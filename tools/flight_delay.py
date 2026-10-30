#!/usr/bin/env python3
"""Compatibility entry point for tools/postprocess/flight_delay.py."""
from pathlib import Path
import runpy

if __name__ == "__main__":
    import sys
    directory = Path(__file__).resolve().parent / "postprocess"
    sys.path.insert(0, str(directory))
    runpy.run_path(str(directory / "flight_delay.py"), run_name="__main__")
