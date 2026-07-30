"""Plot empirical ITRDB Stage B manuscript prototypes from frozen summaries."""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
COLORS = {
    "ca533": "#0072B2",
    "germ012": "#D55E00",
    "fran003": "#009E73",
}
LABELS = {
    "ca533": "ca533",
    "germ012": "GERM012",
    "fran003": "FRAN003",
}


def save_figure(fig: plt.Figure, output_dir: Path, stem: str) -> None:
    fig.savefig(output_dir / f"{stem}.svg", bbox_inches="tight")
    fig.savefig(output_dir / f"{stem}.pdf", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo-root", type=Path, default=REPO_ROOT)
    parser.add_argument("--results-dir", type=Path, default=None)
    parser.add_argument("--output-dir", type=Path, default=None)
    args = parser.parse_args()

    repo_root = args.repo_root.resolve()
    results_dir = (
        args.results_dir.resolve()
        if args.results_dir is not None
        else repo_root / "analysis" / "itrdb_empirical_stage_b_v01" / "results"
    )
    output_dir = (
        args.output_dir.resolve()
        if args.output_dir is not None
        else results_dir / "figures"
    )
    output_dir.mkdir(parents=True, exist_ok=True)

    structure = pd.read_csv(results_dir / "structure_summary_v01.csv")
    primary = structure.loc[structure["detrending_method"] == "AgeDepSpline"].copy()
    subset = pd.read_csv(results_dir / "subset_summary_v01.csv")
    split = pd.read_csv(results_dir / "split_half_summary_v01.csv")
    bootstrap = pd.read_csv(results_dir / "bootstrap_summary_v01.csv")
    loto = pd.read_csv(results_dir / "loto_summary_v01.csv")

    fig, ax = plt.subplots(figsize=(6.3, 4.5))
    for row in primary.itertuples(index=False):
        color = COLORS[str(row.dataset)]
        ax.scatter(
            row.eps,
            row.effective_dimension_fraction,
            s=80,
            color=color,
            edgecolor="black",
            linewidth=0.6,
            zorder=3,
        )
        ax.annotate(
            LABELS[str(row.dataset)],
            (row.eps, row.effective_dimension_fraction),
            xytext=(5, 5),
            textcoords="offset points",
        )
    ax.set_xlabel("Conventional plug-in EPS")
    ax.set_ylabel(r"Effective dimension / $N$")
    ax.set_xlim(0, 1.02)
    ax.set_ylim(0, 1.02)
    ax.grid(alpha=0.25)
    save_figure(fig, output_dir, "figE01_structure_eps_deff_v01")

    fig, ax = plt.subplots(figsize=(7.0, 4.8))
    for dataset, block in subset.groupby("dataset", sort=True):
        block = block.sort_values("subset_size")
        color = COLORS[str(dataset)]
        ax.fill_between(
            block["subset_size"],
            block["correlation_q05"],
            block["correlation_q95"],
            color=color,
            alpha=0.18,
            linewidth=0,
        )
        ax.plot(
            block["subset_size"],
            block["correlation_median"],
            color=color,
            linewidth=2,
            label=LABELS[str(dataset)],
        )
        ax.scatter(
            block["subset_size"],
            block["correlation_min"],
            color=color,
            s=10,
            alpha=0.65,
        )
    ax.set_xlabel("Subset tree count")
    ax.set_ylabel("Subset-to-full correlation")
    ax.set_ylim(-0.05, 1.02)
    ax.grid(alpha=0.25)
    ax.legend(frameon=False)
    save_figure(fig, output_dir, "figE02_subset_profiles_v01")

    split = split.sort_values("correlation_median")
    fig, ax = plt.subplots(figsize=(6.4, 4.3))
    positions = range(len(split))
    for position, row in zip(positions, split.itertuples(index=False)):
        color = COLORS[str(row.dataset)]
        ax.plot(
            [row.correlation_q05, row.correlation_q95],
            [position, position],
            color=color,
            linewidth=4,
            solid_capstyle="round",
        )
        ax.scatter(
            row.correlation_median,
            position,
            color=color,
            edgecolor="black",
            linewidth=0.5,
            s=70,
            zorder=3,
        )
        ax.scatter(
            row.correlation_min,
            position,
            color=color,
            marker="x",
            s=50,
            linewidth=1.5,
            zorder=3,
        )
    ax.set_yticks(list(positions), [LABELS[str(value)] for value in split["dataset"]])
    ax.set_xlabel("Disjoint split-half correlation")
    ax.set_xlim(-0.05, 1.02)
    ax.grid(axis="x", alpha=0.25)
    save_figure(fig, output_dir, "figE03_split_half_intervals_v01")

    bootstrap = bootstrap.sort_values("variance_ratio_median")
    fig, ax = plt.subplots(figsize=(6.4, 4.3))
    positions = range(len(bootstrap))
    for position, row in zip(positions, bootstrap.itertuples(index=False)):
        color = COLORS[str(row.dataset)]
        ax.plot(
            [row.variance_ratio_q05, row.variance_ratio_q95],
            [position, position],
            color=color,
            linewidth=4,
            solid_capstyle="round",
        )
        ax.scatter(
            row.variance_ratio_median,
            position,
            color=color,
            edgecolor="black",
            linewidth=0.5,
            s=70,
            zorder=3,
        )
    ax.axvline(1.0, linestyle="--", linewidth=1, color="black")
    ax.set_yticks(
        list(positions), [LABELS[str(value)] for value in bootstrap["dataset"]]
    )
    ax.set_xlabel("Bootstrap / full chronology variance ratio")
    ax.grid(axis="x", alpha=0.25)
    save_figure(fig, output_dir, "figE04_bootstrap_variance_ratio_v01")

    fig, ax = plt.subplots(figsize=(6.3, 4.5))
    for row in loto.itertuples(index=False):
        color = COLORS[str(row.dataset)]
        ax.scatter(
            row.maximum_absolute_delta_eps,
            row.maximum_standardized_rms_difference,
            s=80,
            color=color,
            edgecolor="black",
            linewidth=0.6,
        )
        ax.annotate(
            LABELS[str(row.dataset)],
            (
                row.maximum_absolute_delta_eps,
                row.maximum_standardized_rms_difference,
            ),
            xytext=(5, 5),
            textcoords="offset points",
        )
    ax.set_xlabel(r"Maximum $|\Delta \mathrm{EPS}|$ under LOTO")
    ax.set_ylabel("Maximum standardized RMS chronology change")
    ax.grid(alpha=0.25)
    save_figure(fig, output_dir, "figE05_loto_domains_v01")

    print(f"Empirical Stage B figures written to: {output_dir}")


if __name__ == "__main__":
    main()
