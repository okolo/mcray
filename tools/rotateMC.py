from __future__ import print_function
__author__ = 'ok'

import sys
import numpy as np
import math

# we assume that in MC all particles were emitted from point (0,0,0) with momenta directed along z axes
# if jet direction angle > 0 we will assume that the jet axis is located in y-z plane

if __name__ == '__main__':
    pass

if len(sys.argv) != 4:
    print('utility to rotate MC output data\nUsage: ' + sys.argv[
        0] + ' <data-file> <jet halh-angle> <jet direction angle>', file=sys.stderr)
    sys.exit(1)
fileName = sys.argv[1]

jet_half_size = np.double(sys.argv[2]) / 180. * math.pi  # jet openning angle
jet_direction = np.double(sys.argv[3]) / 180. * math.pi  # angle between direction to source and jet axes

if len(sys.argv) > 6:
    seed = int(sys.argv[6])
    np.random.seed(seed)

if jet_half_size > math.pi or jet_half_size <= 0:
    jet_half_size = math.pi
# here we assume that jet axis is located in y-z plane if angle=0, then jet points towards observer
jet_axes = [0, np.sin(jet_direction), np.cos(jet_direction)]

if jet_direction > math.pi or jet_direction < 0:
    print('invalid param <jet direction angle> [0, 180] range expected')
    sys.exit(1)

#header_printed = False

# print 'reading from ' + fileName
# data = np.genfromtxt(fileName, dtype='double')
# print str(data.shape[0]) + ' records found'
with open(fileName, "r") as lines:
    column_names = []
    for line in lines:
        columns = line.split()
        if len(columns) == 0:
            continue
        if columns[0] == '#':
            if 'weight' in columns:
                column_names = columns[1:3] + ["Theta","Phi","dT_raw/T","dT_calc/T"] + columns[9:]
                print('# ' + '\t'.join(column_names))
            else:
                print(line.strip())
        else:

            particle = np.array(list(map(np.double, columns)))
            beta = phi = dTcalc = 0

            # if data.shape[0]<1:
            #    sys.exit(0)
            # for particle in data:

            r = (particle[2:5])
            dt = r[2]  # time delay in units of full path (assuming small deflection angle)
            r[2] = 1. - dt  # now r is position of particle
            # normalizing r to 1
            modR = np.dot(r, r) ** 0.5
            r = 1. / modR * r
            # direction of momentum
            n = (particle[5:8])
            n[2] = 1. - n[2]
            dr2 = r[0] ** 2 + r[1] ** 2
            if dr2 > 0:
                # building rotation matrix to move x,y,z to 0,0,1
                norm = 1. / (dr2 ** 0.5)
                # rotation axes:
                ux = norm * r[1]
                uy = -norm * r[0]
                # uz=0
                cos = r[2]
                sin = (1 - cos ** 2) ** 0.5
                toCenter = np.matrix([[cos + ux * ux * (1 - cos), ux * uy * (1 - cos), uy * sin],
                                      [ux * uy * (1 - cos), cos + uy * uy * (1 - cos), -ux * sin],
                                      [-uy * sin, ux * sin, cos]])

                # test
                error = np.fabs( np.dot(toCenter, r) - np.array([0.,0.,1.]) ).max()
                assert error < 1e-7
                zAxes = [0, 0, 1]
                # finding starting direction of the momentum
                nInitial = np.dot(toCenter, zAxes)
                # calculating angle between starting direction and jet axis
                startCos = np.dot(nInitial, jet_axes)
                startAngle = 0
                if startCos < 1:  # avoid roundout errors if angle is very close to zero
                    startAngle = np.arccos(startCos)
                # TODO: here instead of using step function we can modify weight depending
                # on angular distance from jet axes
                if startAngle < jet_half_size:  # if the particle is inside the cone (passes cut)
                    # calculate angle between arrival direction and direction to source
                    cosBeta = np.dot(r, n)  # this angle is invariant (doesn't depend on turns)
                    if cosBeta < 1:  # avoid rounding errors if angle is very close to zero
                        beta = np.arccos(cosBeta)
                        nFinal = np.dot(toCenter, n).transpose()
                        cosPhi = nFinal[0] / np.sin(beta) # second angle (could be useful in assymetric tasks) if beta=0 it is undefined (could be any)
                        if cosPhi < 1 and cosPhi > -1:
                            phi = np.arccos(cosPhi)[0,0]  # returns value in interval [0,pi]
                            if nFinal[1] < 0:
                                phi = 2 * math.pi - phi

                    s = math.sin(startAngle+beta)
                    if s>0.:
                        dTcalc = (math.sin(startAngle)+math.sin(beta))/s - 1.
                else:
                    continue
            elif jet_direction > jet_half_size:  # no deflection: check if jet is pointing towards observer
                continue

            output = [particle[0], particle[1], beta * 180.0 / math.pi, phi * 180.0 / math.pi, dt, dTcalc] + list(particle[8:])
            output = ["{0:g}".format(i) for i in output]
            print(*output)

