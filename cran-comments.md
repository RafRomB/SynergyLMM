## CRAN submission

This is a patch update to the SynergyLMM package (1.1.4) fixing the test failure
reported by the additional checks with OpenBLAS for version 1.1.3.

## Fixes

- The failing tests asserted that `lmmSynergy()` emits its "p-values are
  approximated to 0" warning. That warning is only emitted when a simulated
  p-value is exactly 0, which depends on the Monte Carlo draws returned by
  `MASS::mvrnorm()`. Those draws are built from an eigendecomposition and are
  therefore not reproducible across BLAS/LAPACK implementations, even with a
  fixed seed. The tests no longer assert a simulated outcome, and the wording of
  the warning is now covered by a deterministic unit test.
- Other tests that depended on simulated outcomes, on the random number
  generator state left by previously run test files, or on messages emitted by
  other packages, have been made deterministic.

## R CMD check results

0 errors | 0 warnings | 0 notes
