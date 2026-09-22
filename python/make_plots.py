"""Reproduce the analytical figures and machine-readable assumptions/results."""

from project_paths import ROOT, project_file
import argparse
import csv
import hashlib
import json
from math import ceil, log10
from pathlib import Path
import platform

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from model import accuracy, estimate, workload

ARMS = ["nisq", "star_original", "star_angle", "star_floor", "ft"]
LABELS = ["Physical + whole-circuit PEC", "STAR: original injection",
          "STAR: angle law", "STAR: angle law + floor", "Full FT"]
COLORS = ["#0072B2", "#CC79A7", "#009E73", "#D55E00", "#202020"]
STYLES = ["-", ":", "--", "-", "-"]


def panel(ax, rows, x, y, title):
    for arm, label, color, style in zip(ARMS, LABELS, COLORS, STYLES):
        chosen = [row for row in rows if row["arm"] == arm]
        ax.plot([row[x] for row in chosen], [row[y] for row in chosen],
                label=label, color=color, ls=style, lw=2)
    ax.set_title(title)
    ax.set_xlabel({"L": "Lattice side L (N = L² sites)", "r": "Shared Trotter step count r"}[x])
    ax.grid(alpha=.2)


def write_csv(path, rows):
    fields = sorted(set().union(*(row.keys() for row in rows)))
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=project_file("assumptions.json"))
    parser.add_argument("--out", type=Path, default=ROOT / "results/analytical_v1")
    args = parser.parse_args()
    c = json.loads(args.config.read_text())
    accuracy(c)
    out = args.out
    out.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False,
                         "axes.spines.right": False, "savefig.dpi": 180})
    rows = []
    for alg in ("plaquette", "termwise"):
        for family in ("fixed_tau", "tau_proportional_L"):
            for L in range(2, 26, 2):
                tau = 1.0 if family == "fixed_tau" else L / 4
                w = workload(L, tau, ceil(c["steps_per_tau"] * tau), alg, c["u_over_j"])
                for arm in ARMS:
                    rows.append(dict(family=family, **estimate(w, c, arm)))
    write_csv(out / "scaling.csv", rows)
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), sharex=True, sharey=True)
    for i, alg in enumerate(("plaquette", "termwise")):
        for j, family in enumerate(("fixed_tau", "tau_proportional_L")):
            subset = [row for row in rows if row["algorithm"] == alg and row["family"] == family]
            title = ("Campbell plaquette" if i == 0 else "Termwise Strang")
            title += "; τ = 1" if j == 0 else "; τ = L/4"
            ax = axes[i, j]
            panel(ax, subset, "L", "log10_Q_min", title)
            ax.axhline(log10(c["device_qubits"]), color=".55", ls=":", lw=1)
            ax.set_ylim(0, 12)
            if j == 0:
                ax.set_ylabel("log₁₀ Q required within one week")
    axes[0, 0].legend(fontsize=8, loc="upper left")
    fig.suptitle("Analytical reference curves • full lattice • uncalibrated routing/supply", fontsize=13)
    fig.text(.5, .013, "p = 10⁻³; r = ⌈4τ⌉; 1% relative at D* = 0.1; 95% bounded-sampling criterion\n"
             "Curves leaving the frame continue above 10¹² qubits; CSV retains their logarithms. Dotted horizontal line: Q = 10⁶.",
             ha="center", fontsize=9)
    fig.tight_layout(rect=(0, .055, 1, .955))
    for suffix in ("pdf", "png"):
        fig.savefig(out / f"first_scaling.{suffix}")
    plt.close(fig)

    step_rows = []
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.6), sharey=True)
    for j, alg in enumerate(("plaquette", "termwise")):
        for r in (2, 4, 8, 16, 32, 64, 128):
            w = workload(4, 1, r, alg, c["u_over_j"])
            for arm in ARMS:
                step_rows.append(estimate(w, c, arm))
        subset = [row for row in step_rows if row["algorithm"] == alg]
        panel(axes[j], subset, "r", "log10_Q_min",
              ("Campbell plaquette" if j == 0 else "Termwise Strang") + "; L = 4, τ = 1")
        axes[j].set_xscale("log", base=2)
        axes[j].set_ylim(0, 15)
    axes[0].set_ylabel("log₁₀ Q required within one week")
    axes[0].legend(fontsize=8)
    fig.suptitle("Shared Trotter refinement changes the three costs differently", fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, .94))
    for suffix in ("pdf", "png"):
        fig.savefig(out / f"step_sensitivity.{suffix}")
    plt.close(fig)
    write_csv(out / "step_sensitivity.csv", step_rows)

    # Independent knobs most likely to change interpretation of the first plot.
    sensitivity = []
    variations = [("default", {}), ("threshold_0039", {"surface_threshold": .0039, "surface_A": .07}),
                  ("D_005", {"doublon_reference": .05}), ("D_020", {"doublon_reference": .2}),
                  ("rho_2", {"steps_per_tau": 2}), ("rho_8", {"steps_per_tau": 8}),
                  ("floor_0", {"star_floor": 0}), ("floor_1e-4", {"star_floor": 1e-4}),
                  ("p_1e-4", {"p": 1e-4})]
    for name, change in variations:
        settings = {**c, **change}
        for L in (4, 8, 12, 16):
            w = workload(L, 1, ceil(settings["steps_per_tau"]), "plaquette", c["u_over_j"])
            for arm in ARMS:
                sensitivity.append(dict(sensitivity=name, **estimate(w, settings, arm)))
    write_csv(out / "sensitivity.csv", sensitivity)

    files = [project_file(name) for name in ("model.py", "make_plots.py", "test_model.py",
                                      "analytical_forms.tex", "references.bib", "WORKING_DECISIONS.md")]
    files += [args.config, *sorted((ROOT / "refs").glob("*.pdf"))]
    manifest = {"model_status": "analytical reference; not calibrated; conditional on shared r and D*",
                "assumptions": c, "python": platform.python_version(),
                "numpy": np.__version__, "matplotlib": matplotlib.__version__,
                "sha256": {str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else str(path):
                           hashlib.sha256(path.read_bytes()).hexdigest() for path in files if path.exists()}}
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    baseline = [{key: row[key] for key in ("L", "arm", "log10_Q_min", "distance")}
                for row in rows if row["algorithm"] == "plaquette"
                and row["family"] == "fixed_tau" and row["L"] in (4, 8, 16)]
    print(json.dumps(baseline, indent=2, allow_nan=False))
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
