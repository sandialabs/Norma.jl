# Dynamic Linear OpInf Overlap L2 Error Example

This example mirrors the existing overlap cuboid linear OpInf setup and enables
`compute overlap L2 error: true` on both overlap Schwarz boundary conditions.

Run it with:

```bash
julia --project=. -e 'include("src/Norma.jl"); Norma.run("examples/ahead/overlap/cuboid/dynamic-linear-opinf-rom-overlap-l2-error/cuboids.yaml")'
```

Outputs:

- Standard Exodus results for both subdomains
- `iterations-####.csv` with the Schwarz iteration count per stop
- `overlap-l2-errors-####.csv` with one row per enabled overlap BC:
  `domain,side_set,overlap_l2_error`
