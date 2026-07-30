"""Finalize the dual-hash Yoshino RWL provenance manifest."""
from __future__ import annotations

import argparse
from hashlib import sha256
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]


def file_sha256(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def normalized_lf_bytes(path: Path) -> bytes:
    """Return text bytes after CRLF and lone-CR normalization to LF."""
    data = path.read_bytes()
    return data.replace(b"\r\n", b"\n").replace(b"\r", b"\n")


def normalized_lf_sha256(path: Path) -> str:
    return sha256(normalized_lf_bytes(path)).hexdigest()


def normalized_lf_size(path: Path) -> int:
    return len(normalized_lf_bytes(path))


def first_present(row: pd.Series, *names: str) -> object:
    for name in names:
        if name in row.index and pd.notna(row[name]):
            return row[name]
    raise KeyError(f"none of the manifest fields are present: {names}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo-root", type=Path, default=REPO_ROOT)
    parser.add_argument("--raw-path", type=Path, default=None)
    parser.add_argument("--manifest-path", type=Path, default=None)
    args = parser.parse_args()

    repo_root = args.repo_root.resolve()
    raw_path = (
        args.raw_path.resolve()
        if args.raw_path is not None
        else repo_root / "data" / "yoshino" / "raw" / "nr_all.rwl"
    )
    manifest_path = (
        args.manifest_path.resolve()
        if args.manifest_path is not None
        else repo_root / "data" / "yoshino" / "yoshino_raw_manifest_v01.csv"
    )

    if not raw_path.is_file():
        raise FileNotFoundError(f"Yoshino RWL is absent: {raw_path}")
    if not manifest_path.is_file():
        raise FileNotFoundError(f"Yoshino raw manifest is absent: {manifest_path}")

    old = pd.read_csv(manifest_path, encoding="utf-8-sig")
    if len(old) != 1:
        raise AssertionError("Yoshino raw manifest must contain exactly one row")
    row = old.iloc[0]

    source_size = int(first_present(row, "source_size_bytes", "size_bytes"))
    source_sha = str(first_present(row, "source_sha256", "sha256")).strip().lower()
    if len(source_sha) != 64:
        raise AssertionError("source SHA-256 must contain 64 hexadecimal characters")

    revised = pd.DataFrame(
        [
            {
                "filename": str(first_present(row, "filename")),
                "source_path": str(first_present(row, "source_path")),
                "source_size_bytes": source_size,
                "source_sha256": source_sha,
                "normalized_lf_size_bytes": normalized_lf_size(raw_path),
                "normalized_lf_sha256": normalized_lf_sha256(raw_path),
                "input_level": str(first_present(row, "input_level")),
            }
        ]
    )
    revised.to_csv(manifest_path, index=False, lineterminator="\n")

    native_sha = file_sha256(raw_path)
    native_size = raw_path.stat().st_size
    print("Yoshino dual-hash raw manifest finalized.")
    print(revised.to_string(index=False))
    print(f"Current working-tree bytes: {native_size} bytes; SHA-256 {native_sha}")
    source_equal = native_size == source_size and native_sha == source_sha
    print(
        "Source-byte equality: "
        + ("yes" if source_equal else "no (expected after Git EOL normalization)")
    )


if __name__ == "__main__":
    main()
