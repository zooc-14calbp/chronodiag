"""Create frozen Yoshino empirical Stage C figures."""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]


def save(fig, output_dir: Path, stem: str) -> None:
    fig.savefig(output_dir / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(output_dir / f"{stem}.svg", bbox_inches="tight")
    fig.clf()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo-root", type=Path, default=REPO_ROOT)
    parser.add_argument("--result-dir", type=Path, default=None)
    parser.add_argument("--output-dir", type=Path, default=None)
    args = parser.parse_args()

    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise SystemExit(
            "matplotlib is required for Yoshino figures; install python[figures]"
        ) from exc

    repo_root = args.repo_root.resolve()
    result_dir = (
        args.result_dir.resolve()
        if args.result_dir is not None
        else repo_root / "analysis" / "yoshino_empirical_stage_c_v01" / "results"
    )
    output_dir = (
        args.output_dir.resolve()
        if args.output_dir is not None
        else result_dir / "figures"
    )
    output_dir.mkdir(parents=True, exist_ok=True)

    subset = pd.read_csv(result_dir / "subset_summary_v01.csv")
    split_raw = pd.read_csv(result_dir / "yoshino_split_half_results_v01.csv.gz")
    bootstrap_raw = pd.read_csv(result_dir / "yoshino_bootstrap_results_v01.csv.gz")
    loto_raw = pd.read_csv(
        result_dir / "yoshino_loto_results_v01.csv.gz",
        dtype={"omitted_unit_id": str},
    )
    structure = pd.read_csv(result_dir / "structure_summary_v01.csv")
    keymap = pd.read_csv(
        repo_root / "data" / "yoshino" / "yoshino_keymap_v01.csv",
        dtype={"tree_id": str, "source_group": str},
    )[["tree_id", "source_group"]]

    # Y1: subset sensitivity profile.
    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    ax.fill_between(
        subset["subset_size"],
        subset["correlation_q05"],
        subset["correlation_q95"],
        alpha=0.25,
        label="5--95%",
    )
    ax.plot(
        subset["subset_size"],
        subset["correlation_median"],
        linewidth=2,
        label="Median",
    )
    ax.plot(
        subset["subset_size"],
        subset["correlation_min"],
        linewidth=1,
        linestyle="--",
        label="Sampled minimum",
    )
    ax.set(xlabel="Subset specimen count", ylabel="Subset-to-full correlation", ylim=(0, 1.01))
    ax.legend(frameon=False)
    save(fig, output_dir, "figY01_subset_profile_v01")

    # Y2: disjoint split-half distribution.
    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    valid_split = split_raw.loc[split_raw["status"] == "ok", "correlation"]
    ax.hist(valid_split, bins=35)
    ax.axvline(valid_split.median(), linewidth=2, label="Median")
    ax.axvline(valid_split.quantile(0.05), linestyle="--", label="5th percentile")
    ax.set(xlabel="Disjoint split-half correlation", ylabel="Number of partitions")
    ax.legend(frameon=False)
    save(fig, output_dir, "figY02_split_half_distribution_v01")

    # Y3: bootstrap variance ratio.
    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    valid_boot = bootstrap_raw.loc[
        bootstrap_raw["status"].isin(["ok", "constant_chronology"])
    ]
    ax.hist(valid_boot["variance_ratio"], bins=35)
    ax.axvline(1.0, linestyle="--", linewidth=1.5, label="Full chronology")
    ax.set(
        xlabel="Bootstrap / full chronology variance ratio",
        ylabel="Number of bootstrap draws",
    )
    ax.legend(frameon=False)
    save(fig, output_dir, "figY03_bootstrap_variance_ratio_v01")

    # Y4: LOTO influence domains, with source-group labels.
    loto = loto_raw.merge(
        keymap,
        left_on="omitted_unit_id",
        right_on="tree_id",
        how="left",
        validate="one_to_one",
    )
    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    for group, block in loto.groupby("source_group", sort=True):
        ax.scatter(
            block["delta_eps"].abs(),
            block["standardized_rms_difference"],
            label=f"Source group {group}",
            alpha=0.8,
        )
    ax.set(
        xlabel=r"$|\Delta \mathrm{EPS}|$ under specimen omission",
        ylabel="Standardized RMS chronology change",
    )
    ax.legend(frameon=False)
    save(fig, output_dir, "figY04_loto_influence_domains_v01")

    # Y5: preprocessing structure sensitivity.
    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    order = ["AgeDepSpline", "ModNegExp", "Spline30"]
    block = structure.set_index("detrending_method").reindex(order)
    ax.plot(order, block["eps"], marker="o", label="Conventional plug-in EPS")
    ax.plot(
        order,
        block["effective_dimension_fraction"],
        marker="s",
        label=r"$D_{\mathrm{eff}}/N$",
    )
    ax.set(xlabel="Detrending method", ylabel="Coefficient or normalized dimension", ylim=(0, 1.02))
    ax.legend(frameon=False)
    save(fig, output_dir, "figY05_preprocessing_structure_v01")

    print(f"Yoshino Stage C figures written to: {output_dir}")


if __name__ == "__main__":
    main()
