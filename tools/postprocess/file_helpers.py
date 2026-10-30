"""Read normalized CRbeam tables. Coordinates exposed by CRbeamFile are in kpc.

Requires NumPy; reading config.yaml also requires PyYAML. No plotting imports.
"""

import io
import math
from pathlib import Path
import re

import numpy as np

NUMBER = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?"


def header_number(header, name):
    match = re.search(r"(?:^|\t)" + re.escape(name) + r"\s*=\s*(" + NUMBER + r")(?=\s|$)", header)
    if match is None:
        raise ValueError("Missing or invalid header field: " + name)
    value = np.longdouble(match.group(1))
    if not np.isfinite(value):
        raise ValueError("Non-finite header field: " + name)
    return value


class CRbeamFile:
    def __init__(self, fname, verbose=False, *, config_path=None):
        self.filename = str(fname)
        path = Path(fname)
        headers, rows = [], []
        with path.open() as stream:
            for line in stream:
                if line.lstrip().startswith("#"):
                    headers.append(line.rstrip("\n"))
                elif line.strip():
                    rows.append(line)
        if not headers:
            raise ValueError("Missing CRbeam header")
        self.header = headers[0]
        if "spatial_format = " in self.header and header_number(self.header, "spatial_format") != 1:
            raise ValueError("Unsupported spatial output version")
        self.scale_Mpc = header_number(self.header, "T/Mpc")
        if self.scale_Mpc <= 0:
            raise ValueError("Coordinate scale must be positive")
        column_headers = [line for line in headers if line.startswith("# E/eV\t")]
        if len(column_headers) != 1:
            raise ValueError("Expected one tab-separated column header")
        self.columns = column_headers[0][2:].split("\t")
        if len(set(self.columns)) != len(self.columns):
            raise ValueError("Duplicate column names")
        self.ftable = (np.loadtxt(io.StringIO("".join(rows)), ndmin=2, dtype=np.longdouble) if rows
                       else np.empty((0, len(self.columns)), dtype=np.longdouble))
        if self.ftable.shape[1] != len(self.columns) or not np.isfinite(self.ftable).all():
            raise ValueError("Invalid table shape or non-finite values")
        required = ["E/eV", "weight", "x/T", "y/T", "1-(z/T)", "Px/P", "Py/P", "1-(Pz/P)"]
        if self.columns[:8] != required:
            raise ValueError("Unsupported coordinate layout")
        self.Norm = self.scale_Mpc * 1000
        self.weights = self["weight"]
        self.x = self["x/T"] * self.Norm
        self.y = self["y/T"] * self.Norm
        self.z = (1 - self["1-(z/T)"]) * self.Norm
        if not np.isfinite(np.column_stack((self.x, self.y, self.z))).all():
            raise ValueError("Reconstructed coordinates are not finite")
        if "in_target" in self.columns and not np.isin(self["in_target"], [0, 1]).all():
            raise ValueError("in_target must contain only zero or one")
        self.Px, self.Py = self["Px/P"], self["Py/P"]
        self.Pz = 1 - self["1-(Pz/P)"]
        self.in_target = (self["in_target"].astype(bool) if "in_target" in self.columns
                          else np.zeros(len(self), dtype=bool))
        self.in_target_num = int(self.in_target.sum())
        self.config, self.params_dict = {}, {}
        directory = path.resolve().parent
        run_directory = directory.parent if directory.name.startswith("z") else directory
        config = Path(config_path) if config_path else run_directory / "config.yaml"
        if config.exists():
            import yaml
            with config.open() as stream:
                self.config = yaml.safe_load(stream)
            source = self.config.get("source", {})
            observer = self.config.get("observer", {})
            self.params_dict["z"] = source.get("redshift")
            for prefix, vector in [("source", source.get("position_Mpc", [0, 0, 0])),
                                   ("observer", observer.get("center_Mpc", [0, 0, 0]))]:
                for name, value in zip("xyz", vector):
                    self.params_dict[prefix + "-" + name] = float(value) * 1000
            if "radius_Mpc" in observer:
                self.params_dict["observer-radius"] = float(observer["radius_Mpc"]) * 1000
            angular = source.get("angular_distribution", {})
            if "half_angle_deg" in angular:
                self.params_dict["jet-angle"] = math.radians(angular["half_angle_deg"])
        else:
            params = run_directory / "params"
            if params.exists():
                for line in params.read_text().splitlines():
                    if line.startswith("--"):
                        key, _, value = line.partition("\t")
                        self.params_dict[key[2:]] = value
        try:
            self.source_position_Mpc = np.array([
                header_number(self.header, "source_" + name + "_Mpc") for name in "xyz"], dtype=np.longdouble)
        except ValueError:
            if "spatial_format = " in self.header:
                raise
            self.source_position_Mpc = np.array([
                float(self.params_dict.get("source-" + name, 0)) / 1000 for name in "xyz"], dtype=np.longdouble)
        if not np.isfinite(self.source_position_Mpc).all():
            raise ValueError("Invalid source position")
        for name, value in zip("xyz", self.source_position_Mpc * 1000):
            setattr(self, "source_" + name, value)
        if verbose:
            print("Read", len(self), "particles from", fname)

    def __getitem__(self, name):
        if name in self.columns:
            return self.ftable[:, self.columns.index(name)]
        return self.params_dict[name]

    def __len__(self):
        return len(self.ftable)

    def keys(self):
        return set(self.columns) | set(self.params_dict)

    def delay_Mpc(self, reference="radial", axis=None):
        """Elapsed simulation flight minus a geometric light-travel reference.

        radial: distance from primary birth position to registration position.
        axis: displacement projected onto the supplied unit-normalized axis.
        This is not an observer-plane intersection or a trajectory extrapolation.
        """
        if "FlightMpc" not in self.columns:
            raise ValueError("Per-particle FlightMpc is required; header T is only a coordinate scale")
        flight = self["FlightMpc"]
        if (flight < 0).any():
            raise ValueError("Negative flight time")
        displacement = np.column_stack((self.x, self.y, self.z)) / 1000 - self.source_position_Mpc
        if reference == "radial":
            distance = np.hypot(np.hypot(displacement[:, 0], displacement[:, 1]), displacement[:, 2])
        elif reference == "axis":
            axis = np.asarray(axis, dtype=np.longdouble)
            if axis.shape != (3,) or not np.isfinite(axis).all() or not np.any(axis):
                raise ValueError("A finite nonzero reference axis is required")
            axis = axis / np.max(np.abs(axis))
            axis /= np.linalg.norm(axis)
            distance = displacement @ axis
        else:
            raise ValueError("Unknown geometric reference: " + reference)
        return flight - distance
