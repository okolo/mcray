# Result postprocessing

This directory contains Python utilities for reading, transforming and
analyzing CRbeam results. Run examples from the repository root. Build/run
scripts and simulation tables remain in bin; input-data preparation tools
remain in tools. Regression checks live in tests.

## Utilities

| Utility | Purpose | Detailed documentation |
| --- | --- | --- |
| [file_helpers.py](file_helpers.py) | CRbeamFile loader: reconstruct coordinates, read run settings, and calculate radial or axis-referenced geometric excess times. | [Spatial tables and Python loading](../../docs/usage/spatial_source_observer.md#python-loading-and-delay) |
| [flight_delay.py](flight_delay.py) | Command that writes energy, weight and geometric excess time in Mpc to stdout. | [Delay definitions and invocation](../../docs/usage/spatial_source_observer.md#python-loading-and-delay) |
| [mc_rotate.py](mc_rotate.py) | Rotate original +z-pencil histories, apply jet/PSF cuts and estimate delays by straight continuation to a sphere. | [Algorithm, input assumptions and block processing](mc_rotate.md) |

The utilities support different input contracts. In particular, mc_rotate
requires the original fixed-flight-time geometry; its requirements differ
from the general spatial-table loader. Check the linked documentation before
processing a new result format.

## Dependencies and use

Install the libraries listed in [requirements.txt](requirements.txt):

```sh
python3 -m pip install -r tools/postprocess/requirements.txt
```

NumPy is used for numerical processing; PyYAML is used by file_helpers.py
when reading config.yaml or run.yaml. These are optional analysis dependencies; CRbeam
itself does not require them.

The command-line utility can be run directly:

```sh
python3 tools/postprocess/flight_delay.py --help
```

For scripts or notebooks, add tools/postprocess to the module search path and
import the required module. Function arguments, units, examples and tests are
listed in the detailed documentation linked above.

The previous modern paths bin/FileHelpers.py and tools/flight_delay.py remain
compatibility entry points. Their maintained implementation is in this folder.
External historical forks are read-only references and are not dependencies.
