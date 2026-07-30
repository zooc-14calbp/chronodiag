"""Create manuscript-ready figures from frozen same-EPS Stage A outputs.

The script reads only the fixed correlation matrices and machine-readable tables
written by ``run_same_eps_stage_a_v01.py``.  It does not rerun or alter the
scientific analysis.  Each figure is written as a separate plot to make later
manuscript panel assembly explicit and reversible.

Run from the repository root after installing the optional figure dependencies:

    python -m pip install -e "python[figures]"
    python scripts/plot_same_eps_stage_a_v01.py
"""
from __future__ import annotations

import argparse
from hashlib import sha256
from pathlib import Path
import sys
from typing import Iterable

import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
PYTHON_SRC = REPO_ROOT / "python" / "src"
if str(PYTHON_SRC) not in sys.path:
    sys.path.insert(0, str(PYTHON_SRC))

from chronodiag import eps_value  # noqa: E402

SCENARIO_ORDER = (
    "exchangeable",
    "unequal_loadings",
    "influential_trees",
    "two_clusters",
)
SCENARIO_LABELS = {
    "exchangeable": "Exchangeable",
    "unequal_loadings": "Unequal loadings",
    "influential_trees": "High-loading minority",
    "two_clusters": "Two clusters",
}
SCENARIO_COLORS = {
    "exchangeable": "#1f77b4",
    "unequal_loadings": "#2ca02c",
    "influential_trees": "#ff7f0e",
    "two_clusters": "#9467bd",
}
DEFAULT_FORMATS = ("svg", "png", "pdf")


def _read_tables(analysis_dir: Path) -> dict[str, pd.DataFrame]:
    required = {
        "scenario": "scenario_summary_v01.csv",
        "subset": "subset_summary_v01.csv",
        "split": "split_half_summary_v01.csv",
        "bootstrap": "bootstrap_summary_v01.csv",
        "loto": "loto_results_v01.csv",
    }
    tables: dict[str, pd.DataFrame] = {}
    for name, filename in required.items():
        path = analysis_dir / filename
        if not path.exists():
            raise FileNotFoundError(
                f"required frozen Stage A table was not found: {path}"
            )
        tables[name] = pd.read_csv(path)
    return tables


def _save_figure(
    figure: plt.Figure,
    output_stem: Path,
    *,
    formats: Iterable[str],
    dpi: int,
) -> None:
    output_stem.parent.mkdir(parents=True, exist_ok=True)
    for extension in formats:
        normalized = extension.lower()
        options: dict[str, object] = {"bbox_inches": "tight"}
        if normalized == "png":
            options["dpi"] = dpi
            options["metadata"] = {
                "Software": "chronodiag same-EPS Stage A v01"
            }
        elif normalized == "svg":
            options["metadata"] = {
                "Date": None,
                "Creator": "chronodiag same-EPS Stage A v01",
            }
        elif normalized == "pdf":
            options["metadata"] = {
                "CreationDate": None,
                "ModDate": None,
                "Creator": "chronodiag same-EPS Stage A v01",
            }
        figure.savefig(output_stem.with_suffix(f".{extension}"), **options)
    plt.close(figure)


def _set_common_style() -> None:
    plt.rcParams.update(
        {
            "font.size": 10,
            "axes.titlesize": 11,
            "axes.labelsize": 10,
            "legend.fontsize": 9,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "figure.dpi": 120,
            "savefig.transparent": False,
            "svg.hashsalt": "chronodiag-same-eps-stage-a-v01",
        }
    )


def plot_eps_surface(
    output_dir: Path,
    *,
    formats: Iterable[str],
    dpi: int,
) -> None:
    n_values = np.arange(1, 61)
    r_values = np.linspace(0.0, 0.60, 401)
    n_grid, r_grid = np.meshgrid(n_values, r_values)
    eps_grid = n_grid * r_grid / (1.0 + (n_grid - 1.0) * r_grid)

    figure = plt.figure(figsize=(7.2, 5.4))
    axis = figure.add_axes([0.12, 0.12, 0.72, 0.80])
    image = axis.pcolormesh(
        n_grid,
        r_grid,
        eps_grid,
        shading="auto",
        cmap="viridis",
        norm=Normalize(vmin=0.0, vmax=1.0),
        rasterized=True,
    )
    contour = axis.contour(
        n_grid,
        r_grid,
        eps_grid,
        levels=[0.50, 0.70, 0.85, 0.95],
        linewidths=[0.8, 0.8, 2.0, 0.8],
        colors=["white", "white", "black", "white"],
    )
    axis.clabel(contour, inline=True, fontsize=8, fmt="%.2f")
    stage_a_rbar = 0.23943661971830982
    axis.scatter(
        [18],
        [stage_a_rbar],
        marker="*",
        s=170,
        color="#d62728",
        edgecolor="black",
        linewidth=0.7,
        zorder=5,
        label="Stage A: N = 18, r̄ = 0.2394",
    )
    axis.set_xlabel("Number of biological units, N")
    axis.set_ylabel("Mean interseries correlation, r̄")
    axis.set_xlim(1, 60)
    axis.set_ylim(0, 0.60)
    axis.set_title("Expressed population signal as a joint function of N and r̄")
    axis.legend(loc="lower right", frameon=True)
    color_axis = figure.add_axes([0.87, 0.12, 0.03, 0.80])
    colorbar = figure.colorbar(image, cax=color_axis)
    colorbar.set_label("Conventional plug-in EPS")
    _save_figure(
        figure,
        output_dir / "fig01_eps_surface_v01",
        formats=formats,
        dpi=dpi,
    )


def plot_correlation_matrix(
    scenario: str,
    matrix: np.ndarray,
    effective_dimension: float,
    output_dir: Path,
    *,
    formats: Iterable[str],
    dpi: int,
) -> None:
    display = matrix.copy()
    np.fill_diagonal(display, np.nan)

    figure = plt.figure(figsize=(6.0, 5.4))
    axis = figure.add_axes([0.12, 0.12, 0.72, 0.80])
    cmap = plt.get_cmap("RdBu_r").copy()
    cmap.set_bad("#d9d9d9")
    image = axis.imshow(
        display,
        cmap=cmap,
        vmin=-1.0,
        vmax=1.0,
        interpolation="nearest",
        origin="upper",
    )
    ticks = [0, 2, 5, 8, 11, 14, 17]
    labels = [f"T{index + 1:02d}" for index in ticks]
    axis.set_xticks(ticks, labels=labels, rotation=45, ha="right")
    axis.set_yticks(ticks, labels=labels)
    axis.set_xlabel("Biological unit")
    axis.set_ylabel("Biological unit")
    axis.set_title(
        f"{SCENARIO_LABELS[scenario]}\n"
        f"N = 18, r̄ = 0.2394, EPS = 0.85, Dₑff = {effective_dimension:.2f}"
    )
    color_axis = figure.add_axes([0.87, 0.12, 0.03, 0.80])
    colorbar = figure.colorbar(image, cax=color_axis)
    colorbar.set_label("Interseries correlation")
    _save_figure(
        figure,
        output_dir / f"fig02_matrix_{scenario}_v01",
        formats=formats,
        dpi=dpi,
    )


def plot_subset_profiles(
    subset: pd.DataFrame,
    output_dir: Path,
    *,
    formats: Iterable[str],
    dpi: int,
) -> None:
    figure, axis = plt.subplots(figsize=(7.4, 5.2))
    for scenario in SCENARIO_ORDER:
        block = subset.loc[subset["scenario"] == scenario].sort_values(
            "subset_size"
        )
        x = block["subset_size"].to_numpy(dtype=float)
        median = block["correlation_median"].to_numpy(dtype=float)
        lower = block["correlation_q05"].to_numpy(dtype=float)
        upper = block["correlation_q95"].to_numpy(dtype=float)
        color = SCENARIO_COLORS[scenario]
        axis.plot(
            x,
            median,
            label=SCENARIO_LABELS[scenario],
            color=color,
            linewidth=2.0,
        )
        axis.fill_between(x, lower, upper, color=color, alpha=0.16, linewidth=0)

    theoretical = (
        subset.loc[subset["scenario"] == "exchangeable"]
        .sort_values("subset_size")
    )
    axis.plot(
        theoretical["subset_size"],
        theoretical["theoretical_exchangeable_correlation"],
        color="black",
        linestyle="--",
        linewidth=1.0,
        label="Exchangeable √SSS reference",
    )
    axis.set_xlabel("Subset size, n")
    axis.set_ylabel("Subset-to-full chronology correlation")
    axis.set_xlim(1, 17)
    axis.set_ylim(0.35, 1.01)
    axis.set_xticks(range(1, 18, 2))
    axis.set_title("Composition spread in exhaustive subset-to-full comparisons")
    axis.legend(loc="lower right", frameon=True)
    axis.grid(axis="y", alpha=0.25)
    _save_figure(
        figure,
        output_dir / "fig03_subset_profiles_v01",
        formats=formats,
        dpi=dpi,
    )


def plot_split_half_intervals(
    split: pd.DataFrame,
    output_dir: Path,
    *,
    formats: Iterable[str],
    dpi: int,
) -> None:
    ordered = split.set_index("scenario").loc[list(SCENARIO_ORDER)]
    y = np.arange(len(ordered))

    figure, axis = plt.subplots(figsize=(7.4, 4.4))
    for position, (scenario, row) in enumerate(ordered.iterrows()):
        color = SCENARIO_COLORS[scenario]
        axis.hlines(
            position,
            row["correlation_min"],
            row["correlation_max"],
            color=color,
            linewidth=1.1,
            alpha=0.65,
        )
        axis.hlines(
            position,
            row["correlation_q05"],
            row["correlation_q95"],
            color=color,
            linewidth=6.0,
            alpha=0.75,
        )
        axis.scatter(
            row["correlation_median"],
            position,
            color=color,
            edgecolor="black",
            linewidth=0.5,
            s=55,
            zorder=4,
        )
    reference = float(ordered.iloc[0]["theoretical_exchangeable_half_eps"])
    axis.axvline(
        reference,
        color="black",
        linestyle="--",
        linewidth=1.0,
        label=f"Exchangeable EPS₉ = {reference:.3f}",
    )
    axis.set_yticks(y, [SCENARIO_LABELS[s] for s in SCENARIO_ORDER])
    axis.invert_yaxis()
    axis.set_xlabel("Correlation between disjoint 9-tree chronologies")
    axis.set_xlim(0.0, 0.85)
    axis.set_title(
        "All 24,310 unique balanced split halves\n"
        "thin = full range; thick = 5–95%; point = median"
    )
    axis.grid(axis="x", alpha=0.25)
    axis.legend(loc="upper left", frameon=True)
    _save_figure(
        figure,
        output_dir / "fig04_split_half_intervals_v01",
        formats=formats,
        dpi=dpi,
    )


def plot_bootstrap_intervals(
    bootstrap: pd.DataFrame,
    output_dir: Path,
    *,
    formats: Iterable[str],
    dpi: int,
) -> None:
    ordered = bootstrap.set_index("scenario").loc[list(SCENARIO_ORDER)]
    y = np.arange(len(ordered))

    figure, axis = plt.subplots(figsize=(7.4, 4.4))
    for position, (scenario, row) in enumerate(ordered.iterrows()):
        color = SCENARIO_COLORS[scenario]
        axis.hlines(
            position,
            row["variance_ratio_q01"],
            row["variance_ratio_q99"],
            color=color,
            linewidth=1.0,
            alpha=0.55,
        )
        axis.hlines(
            position,
            row["variance_ratio_q05"],
            row["variance_ratio_q95"],
            color=color,
            linewidth=6.0,
            alpha=0.75,
        )
        axis.scatter(
            row["variance_ratio_median"],
            position,
            color=color,
            edgecolor="black",
            linewidth=0.5,
            s=55,
            zorder=4,
        )
    axis.axvline(1.0, color="black", linestyle="--", linewidth=1.0)
    axis.set_yticks(y, [SCENARIO_LABELS[s] for s in SCENARIO_ORDER])
    axis.invert_yaxis()
    axis.set_xlabel("Bootstrap/full chronology variance ratio")
    axis.set_title(
        "Frozen 18-out-of-18 tree bootstrap\n"
        "thin = 1–99%; thick = 5–95%; point = median"
    )
    axis.grid(axis="x", alpha=0.25)
    _save_figure(
        figure,
        output_dir / "fig05_bootstrap_variance_ratio_v01",
        formats=formats,
        dpi=dpi,
    )


def plot_loto_influence(
    loto: pd.DataFrame,
    output_dir: Path,
    *,
    formats: Iterable[str],
    dpi: int,
) -> None:
    figure, axis = plt.subplots(figsize=(7.2, 5.2))
    for scenario in SCENARIO_ORDER:
        block = loto.loc[loto["scenario"] == scenario]
        color = SCENARIO_COLORS[scenario]
        axis.scatter(
            block["standardized_rms_difference"],
            block["delta_eps"],
            label=SCENARIO_LABELS[scenario],
            color=color,
            alpha=0.82,
            edgecolor="white",
            linewidth=0.5,
            s=48,
        )
    influential = loto.loc[
        (loto["scenario"] == "influential_trees")
        & (loto["omitted_unit_id"].isin(["T01", "T18"]))
    ]
    for _, row in influential.iterrows():
        offset = (5, 5) if row["omitted_unit_id"] == "T18" else (6, -14)
        axis.annotate(
            row["omitted_unit_id"],
            (row["standardized_rms_difference"], row["delta_eps"]),
            xytext=offset,
            textcoords="offset points",
            fontsize=8,
        )
    axis.axhline(0.0, color="black", linewidth=0.8)
    axis.set_xlabel("LOTO standardized RMS chronology difference")
    axis.set_ylabel("Change in conventional plug-in EPS")
    axis.set_title(
        "LOTO influence is metric-specific\n"
        "correlation with the full chronology exceeds 0.9949 in all cases"
    )
    axis.grid(alpha=0.25)
    axis.legend(loc="upper left", frameon=True)
    _save_figure(
        figure,
        output_dir / "fig06_loto_influence_v01",
        formats=formats,
        dpi=dpi,
    )



def write_figure_manifest(output_dir: Path) -> Path:
    rows: list[dict[str, str]] = []
    for path in sorted(output_dir.glob("*.svg")):
        rows.append(
            {
                "filename": path.name,
                "format": "svg",
                "sha256": sha256(path.read_bytes()).hexdigest(),
                "generator": "scripts/plot_same_eps_stage_a_v01.py",
                "status": "frozen_stage_A_prototype",
            }
        )
    manifest = output_dir / "figure_manifest_v01.csv"
    pd.DataFrame(rows).to_csv(manifest, index=False, lineterminator="\n")
    return manifest

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--analysis-dir",
        type=Path,
        default=REPO_ROOT / "analysis" / "same_eps_stage_a_v01",
    )
    parser.add_argument(
        "--fixture-dir",
        type=Path,
        default=REPO_ROOT / "shared_fixtures",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=(
            REPO_ROOT
            / "analysis"
            / "same_eps_stage_a_v01"
            / "figures"
        ),
    )
    parser.add_argument(
        "--formats",
        nargs="+",
        default=list(DEFAULT_FORMATS),
        choices=["svg", "png", "pdf"],
    )
    parser.add_argument("--dpi", type=int, default=300)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    _set_common_style()
    tables = _read_tables(args.analysis_dir)

    plot_eps_surface(args.output_dir, formats=args.formats, dpi=args.dpi)

    scenario_table = tables["scenario"].set_index("scenario")
    for scenario in SCENARIO_ORDER:
        matrix = np.loadtxt(
            args.fixture_dir / f"{scenario}_correlation.csv",
            delimiter=",",
        )
        plot_correlation_matrix(
            scenario,
            matrix,
            float(scenario_table.loc[scenario, "effective_dimension"]),
            args.output_dir,
            formats=args.formats,
            dpi=args.dpi,
        )

    plot_subset_profiles(
        tables["subset"], args.output_dir, formats=args.formats, dpi=args.dpi
    )
    plot_split_half_intervals(
        tables["split"], args.output_dir, formats=args.formats, dpi=args.dpi
    )
    plot_bootstrap_intervals(
        tables["bootstrap"], args.output_dir, formats=args.formats, dpi=args.dpi
    )
    plot_loto_influence(
        tables["loto"], args.output_dir, formats=args.formats, dpi=args.dpi
    )
    manifest = write_figure_manifest(args.output_dir)

    print(f"Wrote same-EPS Stage A figures to {args.output_dir}")
    print(f"Figure manifest: {manifest}")


if __name__ == "__main__":
    main()
