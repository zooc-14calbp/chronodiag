"""Generate frozen Yoshino Stage C memberships."""
from __future__ import annotations

import argparse
from hashlib import sha256
import importlib.util
import json
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]


def load_stage_b_generator(repo_root: Path):
    path = repo_root / "scripts" / "generate_itrdb_empirical_memberships_v01.py"
    spec = importlib.util.spec_from_file_location("chronodiag_stage_b_memberships", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"could not load {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def file_sha256(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def stable_manifest_path(path: Path, repo_root: Path) -> str:
    """Return a repository-relative path, or an absolute path outside the repository."""
    resolved = path.resolve()
    try:
        return resolved.relative_to(repo_root.resolve()).as_posix()
    except ValueError:
        return resolved.as_posix()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo-root", type=Path, default=REPO_ROOT)
    parser.add_argument("--preprocessed-dir", type=Path, default=None)
    parser.add_argument("--output-dir", type=Path, default=None)
    args = parser.parse_args()

    repo_root = args.repo_root.resolve()
    analysis_root = repo_root / "analysis" / "yoshino_empirical_stage_c_v01"
    preprocessed_dir = (
        args.preprocessed_dir.resolve()
        if args.preprocessed_dir is not None
        else analysis_root / "preprocessed"
    )
    output_dir = (
        args.output_dir.resolve()
        if args.output_dir is not None
        else analysis_root / "memberships"
    )
    output_dir.mkdir(parents=True, exist_ok=True)

    helper = load_stage_b_generator(repo_root)
    tree_path = preprocessed_dir / "tree_level_all_methods_v01.csv"
    tree_data = pd.read_csv(tree_path, dtype={"tree_id": str})
    primary = tree_data.loc[
        (tree_data["dataset"] == "yoshino")
        & (tree_data["detrending_method"] == "AgeDepSpline")
    ].copy()
    if primary.empty:
        raise ValueError("no primary Yoshino tree-level data were found")

    units = tuple(sorted(primary["tree_id"].astype(str).unique()))
    n_trees = len(units)
    if n_trees != 50:
        raise AssertionError(f"expected 50 Yoshino specimens, observed {n_trees}")
    if len(primary) != n_trees * primary["year"].nunique():
        raise AssertionError("Yoshino tree-level data are not a complete grid")

    subset = helper.generate_subset_table("yoshino", units)
    split = helper.generate_split_table("yoshino", units)
    bootstrap = helper.generate_bootstrap_table("yoshino", units)

    paths = {
        "subset": output_dir / "yoshino_subset_memberships_v01.csv.gz",
        "split": output_dir / "yoshino_split_half_memberships_v01.csv.gz",
        "bootstrap": output_dir / "yoshino_bootstrap_multiplicities_v01.csv.gz",
    }
    helper.write_csv_gzip_reproducible(subset, paths["subset"])
    helper.write_csv_gzip_reproducible(split, paths["split"])
    helper.write_csv_gzip_reproducible(bootstrap, paths["bootstrap"])

    subset_counts = (
        subset[["replicate_id", "subset_size"]]
        .drop_duplicates()
        .groupby("subset_size", sort=True)
        .size()
    )
    design = pd.DataFrame(
        [
            {
                "dataset": "yoshino",
                "n_trees": n_trees,
                "subset_min_size": int(subset_counts.index.min()),
                "subset_max_size": int(subset_counts.index.max()),
                "subset_replicates_total": int(subset["replicate_id"].nunique()),
                "subset_replicates_min_per_size": int(subset_counts.min()),
                "subset_replicates_max_per_size": int(subset_counts.max()),
                "split_group_a_size": n_trees // 2,
                "split_group_b_size": n_trees - n_trees // 2,
                "split_replicates": int(split["replicate_id"].nunique()),
                "bootstrap_replicates": int(bootstrap["replicate_id"].nunique()),
            }
        ]
    )
    design.to_csv(output_dir / "membership_design_summary_v01.csv", index=False)

    manifest = {
        "analysis": "yoshino_empirical_stage_c_memberships",
        "version": "v01",
        "helper_script_sha256": file_sha256(
            repo_root / "scripts" / "generate_itrdb_empirical_memberships_v01.py"
        ),
        "master_seed": int(helper.MASTER_SEED),
        "max_subsets_per_size": int(helper.MAX_SUBSETS_PER_SIZE),
        "split_replicates": int(helper.SPLIT_REPLICATES),
        "bootstrap_replicates": int(helper.BOOTSTRAP_REPLICATES),
        "n_trees": n_trees,
        "tree_ids": list(units),
        "files": {
            label: {
                "path": stable_manifest_path(path, repo_root),
                "sha256": file_sha256(path),
                "rows": int(
                    {
                        "subset": len(subset),
                        "split": len(split),
                        "bootstrap": len(bootstrap),
                    }[label]
                ),
            }
            for label, path in paths.items()
        },
    }
    (output_dir / "membership_manifest_v01.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    print("Yoshino empirical Stage C membership generation completed.")
    print(design.to_string(index=False))


if __name__ == "__main__":
    main()
