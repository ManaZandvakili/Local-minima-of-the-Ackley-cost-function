import importlib.util
import os
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
os.environ.setdefault("MPLBACKEND", "Agg")

R_REF = 3.15 / 2.0          # Final.py normalises all positions by this radius
# Final.py hull parameters (normalised) -> physical semi-axes = params * R_REF
FINAL_SHAPES = {
    "sphere": ("sphere", None, (1.0, 1.0, 1.0)),
    "cigar": ("ellipsoid", (1.0, 1.0, 2.5), (1.0, 1.0, 2.5)),
    "oblate": ("ellipsoid", (2.0, 2.0, 0.5), (2.0, 2.0, 0.5)),
}


@pytest.fixture(scope="session")
def final():
    """The unmodified conference-paper module (reference/Final.py)."""
    path = ROOT / "reference" / "Final.py"
    spec = importlib.util.spec_from_file_location("final_reference", path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["final_reference"] = mod
    spec.loader.exec_module(mod)
    return mod
