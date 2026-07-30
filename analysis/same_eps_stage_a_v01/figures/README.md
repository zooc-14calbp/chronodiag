# Same-EPS Stage A figure prototypes

These figures are generated from the frozen machine-readable Stage A tables by:

```text
scripts/plot_same_eps_stage_a_v01.py
```

The repository retains SVG versions for review and exact visual provenance. PDF
and 300-dpi PNG versions are generated locally for manuscript assembly:

```powershell
python -m pip install -e "python[figures]"
python scripts/plot_same_eps_stage_a_v01.py --formats svg pdf png
```

Each plot is intentionally separate. Final multi-panel assembly, panel lettering,
and journal dimensions belong to the manuscript figure pipeline rather than this
scientific-analysis step.

The color mapping is fixed across composition figures:

- blue: exchangeable;
- green: unequal loadings;
- orange: high-loading minority;
- purple: two clusters.

`figure_manifest_v01.csv` records the SHA-256 of the frozen SVG prototypes.
