"""Cross-check the MATLAB results against the Python mirror.

Usage (from the repository root):
    python python/compare_matlab_python.py

Reads results/matlab/*.csv and results/python_reference/*.csv and writes
results/comparison/matlab_vs_python.md. Rows are matched on the case keys;
differences should be at round-off level.
"""
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
M_DIR, P_DIR = ROOT / "results" / "matlab", ROOT / "results" / "python_reference"
OUT = ROOT / "results" / "comparison" / "matlab_vs_python.md"

KEYS = ["pcm", "k_scenario", "h", "t_c_min", "bypass", "convection", "shape", "A_K", "P_min"]
VALS = ["w_opt", "w_opt_fs", "w_opt_lo", "w_opt_hi", "E_opt_kJ_per_m", "fl_opt", "dw", "penalty_pct", "dE_pct"]


def load(path):
    df = pd.read_csv(path)
    for col in ("bypass", "convection"):
        df[col] = df[col].astype(bool)
    return df


def main():
    m, p = load(M_DIR / "summary_all.csv"), load(P_DIR / "summary_all.csv")
    j = p.merge(m, on=KEYS, suffixes=("_py", "_m"), how="outer", indicator=True)
    both = j[j["_merge"] == "both"]
    only_m = j[j["_merge"] == "right_only"]
    only_p = j[j["_merge"] == "left_only"]

    lines = ["# MATLAB vs Python cross-check", "",
             f"Rows in MATLAB summary: {len(m)}. Rows in Python summary: {len(p)}. Matched: {len(both)}.", "",
             "| Column | Max absolute difference |", "| --- | --- |"]
    for c in VALS:
        lines.append(f"| {c} | {(both[c + '_m'] - both[c + '_py']).abs().max():.2e} |")

    groups = only_m[["pcm", "k_scenario", "h"]].drop_duplicates()
    lines += ["", f"Rows only in MATLAB: {len(only_m)} "
              f"({', '.join(f'{r.pcm}/{r.k_scenario}/h={r.h}' for r in groups.itertuples()) or 'none'}).",
              f"Rows only in Python: {len(only_p)}.", ""]

    mp = pd.read_csv(M_DIR / "phase_convergence_Pgt_tc.csv")
    pp = pd.read_csv(P_DIR / "phase_convergence_Pgt_tc.csv")
    mp30 = mp[mp["t_c_min"] == 30].assign(shape=lambda d: d["shape_id"].map({1: "sine", 2: "square"}))
    jj = pp.merge(mp30, left_on=["pcm", "k", "shape", "A", "nph"],
                  right_on=["pcm", "k", "shape", "A_K", "n_phases"])
    lines += ["## Phase convergence (P = 60 min, t_c = 30 min)", "",
              f"Matched rows: {len(jj)}. Max penalty difference: "
              f"{(jj['penalty'] - jj['penalty_pct']).abs().max():.2e} percentage points.", ""]

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
