"""Regenerate CellBender's HTML report for finished samples (e.g. runs whose report step failed).

Usage: python cellbender_report.py <data_dir> <work_dir> [sample ...]
Writes <work_dir>/cellbender/<S>/<S>_cellbender_report.html
"""
import os
import sys
from pathlib import Path

os.environ["PATH"] = os.path.dirname(sys.executable) + os.pathsep + os.environ.get("PATH", "")
from cellbender.remove_background import report  # noqa: E402


def main(data_dir, work_dir, samples):
    data_dir, work_dir = Path(data_dir).resolve(), Path(work_dir).resolve()
    nb = str(Path(report.__file__).with_name("report.ipynb"))
    for s in samples:
        out = work_dir / "cellbender" / s
        os.chdir(out)  # the notebook writes tmp.report.* in the cwd
        os.environ["INPUT_FILE"] = str(data_dir / s / "raw_feature_bc_matrix.h5")
        os.environ["OUTPUT_FILE"] = str(out / f"{s}_cellbender.h5")
        html = report.run_notebook_make_html(file=nb, output=str(out / f"{s}_cellbender_report.html"))
        print(s, html, flush=True)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3:])
