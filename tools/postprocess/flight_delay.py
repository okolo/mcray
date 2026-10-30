#!/usr/bin/env python3
"""Write E/eV, weight and excess simulation flight time in Mpc (c=1)."""
import argparse
import sys

from file_helpers import CRbeamFile
import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("file")
    parser.add_argument("--config", help="Effective config.yaml, normally found beside z*/")
    parser.add_argument("--reference", choices=["radial", "axis"], default="radial")
    parser.add_argument("--axis", type=float, nargs=3, help="Required for the axis reference")
    args = parser.parse_args()
    try:
        data = CRbeamFile(args.file, config_path=args.config)
        delay = data.delay_Mpc(args.reference, args.axis)
    except (ValueError, OSError) as error:
        parser.error(str(error))
    np.savetxt(sys.stdout, np.column_stack((data["E/eV"], data.weights, delay)),
               header="E/eV\tweight\tdelay_Mpc", delimiter="\t", fmt="%.17g")


if __name__ == "__main__":
    main()
