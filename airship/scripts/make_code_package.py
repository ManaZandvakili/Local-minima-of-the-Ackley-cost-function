"""Build the code package for the journal submission (a zip file).

Contains the Python package, settings, the programs needed to reproduce the
results, the tests, the conference-paper reference code, the MATLAB/Simulink
version, and (unless --code-only) the layouts, tuned gains and trained networks
used in the paper. Videos, figures, one-way case-study scripts and development
logs are left out.

Usage (from the airship/ folder):
    python scripts/make_code_package.py --out ../airship_code.zip
    python scripts/make_code_package.py --out ../airship_code.zip --code-only
"""
import argparse
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SCRIPTS = ["export_reference_layouts.py", "make_hand_layouts.py", "run_optimizer.py", "tune_pid.py",
           "train_rl.py", "evaluate_rl.py", "path_following.py", "main_result_figure.py", "export_matlab.py",
           "make_code_package.py"]
CODE = ["airship/**/*.py", "configs/*.yaml", "tests/*.py", "reference/Final.py", "reference/README.md",
        "matlab/*.m", "matlab/README.md", "matlab/.gitignore", "matlab/data/**/*.mat", "requirements.txt"]
RUNS = ["results/layouts/*.json", "results/layouts/phase1_bidir_20261117_145131/*.json",
        "results/layouts/phase1_bidir_maxN11_20261117_151234/*.json",
        "results/runs/phase1_bidir_20261117_145131/report.md",
        "results/runs/phase1_bidir_maxN11_20261117_151234/report.md",
        "results/pid/20261117_153423/gains/*.json", "results/pid/20261117_153423/config.yaml",
        "results/pid/20261117_153423/report.md", "results/pid/20261117_153423/summary.csv"]
for r in ("20261117_155720_bidir_residual_seed0", "20261117_155720_bidir_plain_seed0"):
    RUNS += [f"results/rl/{r}/{f}" for f in ("model.zip", "vecnormalize.pkl", "config.yaml", "curriculum.json",
                                             "pool.pkl")]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(ROOT.parent / "airship_code.zip"))
    ap.add_argument("--code-only", action="store_true", help="leave out the layouts, gains and networks")
    a = ap.parse_args()
    files = set()
    for pat in CODE + [f"scripts/{s}" for s in SCRIPTS] + ([] if a.code_only else RUNS):
        hits = [p for p in ROOT.glob(pat) if p.is_file() and "__pycache__" not in p.parts]
        if not hits:
            raise SystemExit(f"nothing matches {pat}")
        files.update(hits)
    with zipfile.ZipFile(a.out, "w", zipfile.ZIP_DEFLATED) as z:
        z.write(ROOT / "package" / "README.md", "airship_code/README.md")
        for p in sorted(files):
            z.write(p, f"airship_code/{p.relative_to(ROOT)}")
    size = Path(a.out).stat().st_size / 1e6
    print(f"{a.out}: {len(files) + 1} files, {size:.1f} MB")


if __name__ == "__main__":
    main()
