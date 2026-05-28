# Overlap L2 Error Additions

## What changed

- Added an optional `compute overlap L2 error` flag for `Schwarz overlap` and `OpInf Schwarz overlap` boundary conditions.
- Built and cached a full overlap-region interpolation map on the wrapped FOM overlap BC, separate from the interface side-set map used for coupling.
- Added a post-convergence diagnostic that computes the displacement mismatch
  `sqrt(sum_i ||u_dst(x_i) - u_src(x_i)||^2)` over mapped destination overlap nodes.
- Kept Schwarz convergence unchanged: the new metric is diagnostic only and does not affect stopping, relaxation, or BC application.

## Reporting

- Logged the overlap-region L2 error once per enabled overlap BC after a converged Schwarz stop.
- Added a dedicated per-stop CSV file, `overlap-l2-errors-####.csv`, instead of overloading `iterations-####.csv`.
- Reused the wrapped FOM overlap BC for `OpInf Schwarz overlap`, so the OpInf wrapper exposes the same diagnostic value without duplicating overlap-region maps or state.

## Coverage

- Added static overlap tests for parsing/default behavior, overlap-region node selection, post-convergence computation, and dedicated CSV output.
- Added a linear OpInf cuboid regression using a new example deck with the diagnostic enabled.
- Added an `OpInf Schwarz overlap` regression on the neural-network OpInf path to confirm the wrapper forwards the option and exposes the stored diagnostic.

## New example

- Added [`examples/ahead/overlap/cuboid/dynamic-linear-opinf-rom-overlap-l2-error`](/Users/ejparis/Research/ahead/Norma.jl/examples/ahead/overlap/cuboid/dynamic-linear-opinf-rom-overlap-l2-error) as a ready-to-run cuboid example for the new diagnostic.
