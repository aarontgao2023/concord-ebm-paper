"""Dry run of the packaged comparison on a simulated REF cohort (the ADNI-shaped design).

Simulates one REF_H0 (or a named cell) data set, writes it to CSV exactly as an analyst
would hand it over, and runs the command-line entry point on it. Records timing so the
ADNI budget (n=971, B=599, three estimators) can be planned from measured numbers.
"""
import argparse
import json
import os
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from design_v2 import resolve, simulate  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cell", default="REF_H0")
    ap.add_argument("--seed", type=int, default=31900000)
    ap.add_argument("--B", type=int, default=99)
    ap.add_argument("--stability", type=int, default=20)
    ap.add_argument("--workers", type=int, default=16)
    ap.add_argument("--estimators", default="invariant_min,standard")
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--overrides", default="{}", help="JSON overrides for resolve()")
    a = ap.parse_args()
    a.output.mkdir(parents=True, exist_ok=True)
    df, truth = simulate(resolve(a.cell, **json.loads(a.overrides)), a.seed)
    csv = a.output / f"{a.cell}_{a.seed}.csv"
    df.to_csv(csv, index=False)
    (a.output / f"{a.cell}_{a.seed}_truth.json").write_text(json.dumps(
        {k: truth[k] for k in ("base_order", "group_orderings", "biomarker_names", "group_dx_counts",
                               "pair_truth", "common_order_null")}, indent=1))
    print(f"simulated {a.cell} seed {a.seed}: {len(df)} subjects, truth pair_truth={truth['pair_truth']}",
          flush=True)
    for estimator in a.estimators.split(","):
        out = a.output / f"{a.cell}_{a.seed}_{estimator}.json"
        log = a.output / f"{a.cell}_{a.seed}_{estimator}.log"
        cmd = [sys.executable, str(HERE.parent / "ebm_compare.py"), str(csv), "--group", "APOE",
               "--estimator", estimator, "--B", str(a.B), "--stability", str(a.stability),
               "--workers", str(a.workers), "--seed", str(a.seed), "--output", str(out)]
        print("running:", " ".join(cmd), flush=True)
        t0 = time.monotonic()
        with log.open("w") as stream:
            code = subprocess.call(cmd, stdout=stream, stderr=subprocess.STDOUT,
                                   env=os.environ | {"PYTHONUNBUFFERED": "1"})
        print(f"{estimator}: exit {code} in {(time.monotonic() - t0) / 60:.1f} min", flush=True)
        if out.exists():
            r = json.loads(out.read_text())
            print(json.dumps({"estimator": estimator, "status": r["status"],
                              "distances": r.get("distances", {}).get("pairs"),
                              "p_diagnosis": {k: v["p"] for k, v in r["tests"]["diagnosis"]["pairs"].items()}
                              if r.get("tests") else None,
                              "p_unrestricted": {k: v["p"] for k, v in r["tests"]["unrestricted"]["pairs"].items()}
                              if r.get("tests") else None,
                              "ess": r.get("effective_sample_size"), "fits": r.get("fits"),
                              "timing": r.get("timing")}, indent=1), flush=True)
        print(log.read_text()[-3000:], flush=True)


if __name__ == "__main__":
    main()
