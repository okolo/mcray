import math

import numpy as np


def mc_rotate(
    mc_file, jet_half_size_deg, jet_direction_deg, psf_deg=180, rotated_file=None
):
    if rotated_file is None:
        rotated_file = f"{mc_file}_{jet_half_size_deg}_{jet_direction_deg}_{psf_deg}"

    jet_half_size = (
        np.longdouble(jet_half_size_deg) / 180.0 * math.pi
    )  # jet openning angle
    jet_direction = (
        np.longdouble(jet_direction_deg) / 180.0 * math.pi
    )  # angle between direction to source and jet axes

    if jet_half_size > math.pi or jet_half_size <= 0:
        jet_half_size = math.pi
    # here we assume that jet axis is located in y-z plane if angle=0, then jet points towards observer
    jet_axes = np.array(
        [0, np.sin(jet_direction), np.cos(jet_direction)], dtype=np.longdouble
    )

    if jet_direction > math.pi or jet_direction < 0:
        raise ValueError(
            "invalid jet_direction_deg param [0, 180] range expected")

    with open(mc_file, "r") as lines, open(rotated_file, "w") as rotated:
        column_names = []
        for line in lines:
            columns = line.split()
            if len(columns) == 0:
                continue
            if columns[0] == "#":
                if "weight" in columns:
                    column_names = (
                        columns[1:3]
                        + ["Theta", "Phi", "dT_raw/T", "dT_calc/T"]
                        + columns[9:]
                    )
                    print("# " + "\t".join(column_names), file=rotated)
                else:
                    print(line.strip(), file=rotated)
            else:
                particle = np.array(list(map(np.longdouble, columns)))
                beta = phi = 0 #dTcalc = 0

                r = particle[2:5]
                dz = r[
                    2
                    # dz is roughly equal to time delay in units of full path (assuming very small deflection angle)
                ]

                # negative dz may apper in data as a result of rounding errors
                dz = max(dz, 0)  # replacing negative values by zeros

                r[2] = 1.0 - dz  # now r is position of particle

                # direction of momentum
                n = particle[5:8]
                n[2] = 1.0 - n[2]
                # compensate rounding errors for n
                mod_n = np.dot(n, n) ** 0.5
                n = 1 / mod_n * n
                t = 0  # calculated delay
                if dz > 0:  # there is some delay
                    r2 = np.dot(r, r)
                    # assert r2 < 1
                    if r2 > 1:
                        r *= 1.0 / np.sqrt(r2)
                        r2 = 1.0
                    # we continue trajectory in the direction of momentun until it reaches the sphere of radius 1.
                    # to find the delay t, we need to solve the following quadratic equation with respect to variable t:
                    # [(x + V_x)*t]^2 + [(y + V_y)*t]^2 + [(z + V_z)*t]^2 = 1
                    # with extra condition t > 0
                    # below is the solution:
                    rV = np.dot(n, r)
                    t = np.sqrt(1.0 - r2 + rV * rV) - rV

                    # calculating the end point of the trajectory
                    r = r + t * n

                # normalizing r to 1 to fix the rounding errors
                modR = np.dot(r, r) ** 0.5
                r = 1.0 / modR * r

                # now we will construct the rotation of coordinates which converts vector r to 0,0,1
                # to find the starting direction of the particle momentum
                dr2 = r[0] ** 2 + r[1] ** 2
                if dr2 > 0:
                    # building rotation matrix to move x,y,z to 0,0,1
                    norm = 1.0 / (dr2**0.5)
                    # rotation axes:
                    ux = norm * r[1]
                    uy = -norm * r[0]
                    # uz=0
                    cos = r[2]
                    sin = (1 - cos**2) ** 0.5

                    toCenter = np.array(
                        [
                            [cos + ux * ux * (1 - cos), ux *
                             uy * (1 - cos), uy * sin],
                            [ux * uy * (1 - cos), cos + uy *
                             uy * (1 - cos), -ux * sin],
                            [-uy * sin, ux * sin, cos],
                        ]
                    )

                    zAxes = np.array([0.0, 0.0, 1.0], dtype=np.longdouble)
                    # test
                    error = np.fabs(np.dot(toCenter, r) - zAxes).max()
                    assert error < 1e-7
                    # finding starting direction of the momentum
                    nInitial = np.dot(toCenter, zAxes)
                    norm_nInitial = np.sqrt(np.dot(nInitial, nInitial))
                    nInitial *= (
                        1.0 / norm_nInitial
                    )  # normalize after rotation to avoid rounding errors
                    # calculating angle between starting direction and jet axis
                    startCos = np.dot(nInitial, jet_axes)
                    startAngle = 0
                    if (
                        startCos < 1
                    ):  # avoid roundout errors if angle is very close to zero
                        startAngle = np.arccos(startCos)
                    # TODO: here instead of using step function we can modify weight depending
                    # on angular distance from jet axes
                    if (
                        startAngle < jet_half_size
                    ):  # if the particle is inside the cone (passes cut)
                        # calculate angle between arrival direction and direction to source
                        cosBeta = np.dot(
                            r, n
                        )  # this angle is invariant (doesn't depend on turns)
                        if (
                            cosBeta < 1
                        ):  # avoid rounding errors if angle is very close to zero
                            beta = np.arccos(cosBeta)
                            nFinal = np.dot(toCenter, n).transpose()
                            cosPhi = nFinal[0] / np.sin(
                                beta
                                # second angle (could be useful in assymetric tasks) if beta=0 it is undefined (could be any)
                            )
                            if cosPhi < 1 and cosPhi > -1:
                                phi = np.arccos(
                                    cosPhi
                                )  # returns value in interval [0,pi]
                                if nFinal[1] < 0:
                                    phi = 2 * math.pi - phi

                        # s = math.sin(startAngle+beta)
                        # if s>0.:
                        #     dTcalc = (math.sin(startAngle)+math.sin(beta))/s - 1.
                        # dTcalc = (1.-min(modR, 1.))/cosBeta  # condition modR > 1 can only occur due to rounding error
                    else:
                        continue
                elif (
                    jet_direction > jet_half_size
                ):  # no deflection: check if jet is pointing towards observer
                    continue

                if beta * 180.0 / math.pi > psf_deg:
                    continue

                output = [
                    particle[0],
                    particle[1],
                    beta * 180.0 / math.pi,
                    phi * 180.0 / math.pi,
                    dz,
                    t,
                ] + list(particle[8:])
                output = ["{0:g}".format(i) for i in output]
                print(*output, file=rotated)

    return rotated_file
