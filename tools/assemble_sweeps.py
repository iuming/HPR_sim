# -*- coding: utf-8 -*-
"""
Assemble the recorded/partial sweep points (produced by run_sweep.py, captured
from console output when the runner was throttled) into the standard per-sweep
CSV/PNG files and a combined JSON, so make_report.py can consume them.

Only fills in sweeps that are missing; real CSV files are left untouched.
"""
import os
import csv
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(os.path.dirname(HERE), "results")

# recorded values: {cavity: {sweep: {parameter: coverage_ratio}}}
RECORDED = {
    "main": {
        # captured from run_sweep output (fast mode)
        "recip": {1: 0.9609, 2: 0.9611, 3: 0.9571, 4: 0.9564},
    },
    "harmonic": {
        "recip": {1: 0.9759, 2: 0.9732, 3: 0.9733, 4: 0.9728},
    },
}

LABELS = {
    "threshold": "cleaning intensity threshold (cos)",
    "tilt": "hole tilt angle [deg]",
    "recip": "reciprocation count",
    "speed": "wand translation speed [mm/s]",
    "rot": "head rotation speed [rpm]",
    "jetangle": "jet spray half-angle [deg]",
    "pressure": "pump pressure [bar]",
}


def load_existing(name, sw):
    p = os.path.join(RES, f"{name}_sweep_{sw}.csv")
    if os.path.exists(p):
        with open(p, encoding="utf-8") as fh:
            r = list(csv.DictReader(fh))
        return [(float(x["parameter"]), float(x["coverage"])) for x in r]
    return None


def main():
    combined = {}
    for name, sweeps in RECORDED.items():
        combined[name] = {}
        for sw, points in sweeps.items():
            existing = load_existing(name, sw)
            pts = existing if existing is not None else sorted((k, v) for k, v in points.items())
            csv_path = os.path.join(RES, f"{name}_sweep_{sw}.csv")
            if existing is None:
                with open(csv_path, "w", newline="", encoding="utf-8") as fh:
                    w = csv.writer(fh)
                    w.writerow(["parameter", "coverage", "cleaned_m2",
                                "reachable_m2", "wall_time_s", "mean_exposure_s"])
                    for par, cov in pts:
                        w.writerow([par, cov, "", "", "", ""])
            # plot
            x = [p for p, _ in pts]
            y = [c * 100 for c, _ in pts]
            fig, ax = plt.subplots(figsize=(7, 5))
            ax.plot(x, y, "o-", color="crimson", lw=2)
            ax.set_xlabel(LABELS.get(sw, sw))
            ax.set_ylabel("effective coverage ratio [%]")
            ax.set_ylim(0, 100)
            ax.set_title(f"{name}: HPR coverage vs {LABELS.get(sw, sw)}")
            fig.tight_layout()
            fig.savefig(os.path.join(RES, f"{name}_sweep_{sw}.png"), dpi=140)
            plt.close(fig)
            combined[name][sw] = [
                dict(parameter=par, coverage=cov) for par, cov in pts
            ]
    with open(os.path.join(RES, "sweeps_recorded.json"), "w",
              encoding="utf-8") as fh:
        json.dump(combined, fh, indent=2, default=float)
    print("assembled recorded sweep files OK")


if __name__ == "__main__":
    main()
