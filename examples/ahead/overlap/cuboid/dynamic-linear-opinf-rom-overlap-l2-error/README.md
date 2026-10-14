# Dynamic Linear OpInf Overlap L2 Error Example

This example is the overlap cuboid linear OpInf setup with
`compute overlap L2 relative error: disp` enabled on both overlap Schwarz boundary conditions.
The field argument can be `disp` (displacement), `velo` (velocity), or `acce` (acceleration).

Run it with:

```bash
julia --project=. -e 'include("src/Norma.jl"); Norma.run("examples/ahead/overlap/cuboid/dynamic-linear-opinf-rom-overlap-l2-error/cuboids.yaml")'
```

Outputs:

- Standard Exodus results for both subdomains
- `iterations-####.csv` with the Schwarz iteration count per stop
- `overlap-l2-disp-rel-errors-####.csv` with one row per enabled overlap BC:
  `domain,side_set,overlap_l2_error`

