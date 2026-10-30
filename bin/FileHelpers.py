"""Compatibility import; the maintained loader lives in tools/postprocess."""
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools/postprocess"))
from file_helpers import CRbeamFile, header_number
