# Yoshino raw ring-width example

This directory is reserved for the public raw ring-width dataset associated with:

Yonenobu, H., Ohyama, M., Hoshino, Y., Kurihara, S., and Fukuchi, R. (2026).
*Tree-ring evidence from a centuries-long plantation system: Growth trends and
changing hydroclimatic sensitivity of Japanese cedar in Yoshino, western Japan.*
Dendrochronologia 99, 126572. https://doi.org/10.1016/j.dendro.2026.126572

## Deliberately minimal repository metadata

The raw measurement keycode is the authoritative identifier and already embeds the
relevant `site_id`, `specimen_id`, and `tree_id`. The frozen release will document
the substring/parsing rule here, but will not duplicate detailed specimen/site
metadata that are already described in the paper.

The paper reports 60 retained specimens and 72 measured radial directions. Twelve
specimens have two measured directions; those directions are averaged before the
specimen-level chronology is formed. `chronodiag` will reconstruct that aggregation
from the raw keycode rather than distribute a second, manually maintained table of
specimen-mean series.

## Planned file

```text
yoshino_ring_width_raw.csv   # or the original .rwl file, once frozen
```

No climate data or climate-response outputs are required for the chronology-
diagnostic demonstration.
