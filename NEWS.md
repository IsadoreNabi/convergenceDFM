# convergenceDFM 0.3.2

## Cross-platform reproducibility of the regression (golden) tests

* **Portable golden tolerances.** The Leave-Cluster-Out and reweighting golden
  tests previously asserted *bit-exact* equality (`expect_identical`) against
  fixtures frozen on the development machine. The underlying computations run
  through floating-point linear algebra (BLAS), whose last few significant
  digits are not portable across CRAN's heterogeneous platforms / BLAS builds,
  so the bit-exact assertions failed on the CRAN check farm at the ~1e-11
  (Leave-Cluster-Out) and ~1e-7 (reweighting) level. The portable invariant is
  now numerical: `expect_equal()` with tolerances (`1e-6` and `1e-5`
  respectively) that absorb the non-portable noise while staying two-to-five
  orders below the scale of any genuine algorithmic regression. The stronger
  bit-exact assertion is retained as a developer-machine guard via
  `skip_on_cran()`. Logical outputs (`robust`) remain compared bit-exactly, as
  they are platform-invariant. No package code changed; this is a test-rigor
  fix that keeps the three rigor layers (algebraic / statistical / numerical)
  strictly separate.

# convergenceDFM 0.3.1

## Wedge diagnostics for the transformation problem

* **New: `compute_wedge()`, `wedge_stationarity()`, `placebo_values()`** (file
  `R/wedge_diagnostics.R`). The wedge `W_i = Phi_i - V_i = K_i * G' - p_i`
  isolates the marxian redistribution of surplus, cancelling the shared cost
  price `k = c + v`, and sums to zero across sectors by construction. The module
  builds the wedge (with a per-period zero-sum check), tests its per-sector and
  panel mean-reversion / unit-root behaviour (AR(1) reversion speed and half-life
  plus an augmented Dickey-Fuller test via `urca`), and constructs
  aggregate-preserving placebo "values" as negative controls (a uniform-profit-
  rate scheme, where the wedge vanishes and value equals the production price, and
  a seeded derangement that mis-pairs capital and surplus while preserving the
  per-period total). Used to move the evidential weight of the gravitation claim
  from the definitional level co-movement to the construction-orthogonal wedge
  dynamics, and to provide negative controls for the internal adversarial
  validation. Fully tested in `tests/testthat/test-wedge.R`.

# convergenceDFM 0.3.0

Single canonical disaggregation engine and a cross-sectional leave-out test for
coupling robustness under input-output dependence.

## Canonical disaggregation (breaking)

* **Local disaggregation duplicate removed.** The deterministic weight-blend
  `run_disaggregation_custom_prior()` (and its private helpers `compute_L_from_P`,
  `spread_likelihood`, `posterior_weighted`) never conditioned on the observed
  CPI and duplicated the purpose of the dedicated package. They are deleted.
  The single source of truth is now the closed-form conjugate engine
  `BayesianDisaggregation::disaggregate_conjugate()` (a Kalman/RTS smoother,
  pure R), added to `Imports`.
* **`test_reweighting_robustness()` rewired.** Each alternative weighting scheme
  is a constant-in-time prior replicated across periods to form the weight
  matrix; the sectoral levels are the conjugate posterior median, genuinely
  conditioned on the CPI, instead of a deterministic blend.
* **Readers retained and relocated.** The hardened `read_cpi()` (clean error on
  an unidentifiable column) and the internal `read_weights_matrix()` move to
  `R/data_io.R`; they are generic Excel input/output, not part of the
  disaggregation engine, and stay in this package.

## Leave-Cluster-Out

* **`test_leave_cluster_out()` (new, exported).** Generalizes the
  delete-one-sector jackknife to delete-one-cluster: each group of sectors is
  dropped from both `X` and `Y` and the coupling is re-estimated. Under
  input-output ("MIP") dependence a naive leave-one-out is optimistic; dropping
  an entire value chain forces the prediction onto the general gravitation. It is
  the cross-sectional complement of the temporal time-shift / block-bootstrap
  nulls and reuses their coupling pipeline (it does not reimplement it).
* **`build_cluster_map()` (new, exported).** Pluggable sector-to-cluster map. A
  genuine map (the Leontief partition) is supplied via `cluster_map`; when it is
  absent a documented fallback (`"correlation"` or `"com"`) is built and flagged
  as a stopgap proxy, not a demand-linkage partition.
* Delete-a-group (block) jackknife `bias`/`se` are reported with the documented
  balanced-cluster caveat; the primary outputs are per-cluster influence and the
  robustness verdict.

## Reproducibility

* `test_reweighting_robustness()` now propagates its `seed` to the data
  diagnosis and the cross-validated PLS component selection, closing a gap where
  the per-scheme couplings depended on the ambient RNG state (call order). A
  pre-existing error-handler bug that dropped the failure marker for a failed
  scheme (`<-` instead of `<<-`) is fixed.

## Tests and documentation

* New always-on tests (`test-disaggregation-import.R`,
  `test-leave-cluster-out.R`) and regression goldens
  (`test-golden-convergence.R`): the Leave-Cluster-Out reproduces a frozen
  result bit-exactly, and the rewired reweighting reproduces frozen couplings
  within numerical tolerance (the layers are kept separate). Goldens and Excel
  fixtures are regenerated by `validacion/make_golden_convergence.R`.
* New vignette `canonical-disaggregation-and-leave-cluster-out` documents both
  decisions.

# convergenceDFM 0.2.0

Robustness-focused overhaul. This release fixes several correctness bugs and
methodological flaws that could bias or invalidate convergence/coupling
inference, and aligns the documentation with the implemented API.

## Methodology (breaking / behavioural)

* **Coupling null reconstructed.** `rotation_null_test()` previously generated
  its null by random orthogonal rotation of the factor space. The four coupling
  statistics (Procrustes, CCA, principal angles, dynamic-beta norm) are
  *invariant* to orthogonal rotation, so the null was degenerate and the test
  could never reject. The default null is now a **circular time-shift /
  moving-block bootstrap** that preserves each series' marginal dynamics while
  breaking the cross-series temporal alignment. The pure-rotation null is kept
  only as an explicitly documented invariance diagnostic (`null_method =
  "rotation"`).
* **Permutation test reconstructed.** `test_permutation_robustness()` previously
  permuted the *columns* of the Y factors; the coupling strength is invariant to
  that permutation, so the null was degenerate. It now permutes the **time
  index** of Y (circular shift) to break the X->Y relationship.
* **Convergence test de-tautologized.** The OU/AR(1) persistence `phi` is no
  longer hard-constrained to `(0, 1)`; it can express unit-root and mildly
  explosive dynamics, so the convergence verdict is now a genuine test of
  `H0: not stationary` rather than a constraint that is always satisfied. The
  prior on `phi` is weakly informative and centred at a neutral value rather
  than at high persistence.
* **MCMC diagnostics.** `estimate_factor_OU()` now returns and (optionally)
  warns on R-hat, bulk/tail effective sample size and divergent transitions, and
  the previously dead `adapt2`/`mtd2` controls now drive an automatic re-run when
  divergences are detected.
* **Real jackknife.** `test_jackknife_sectors()` now performs a leave-one-sector
  -out jackknife with bias and standard-error estimates and influence ranking,
  instead of a single leave-top-k drop.
* **Out-of-sample significance.** Nested out-of-sample comparisons
  (`deltaR2_ou()`, `rescue_short_run_channel()`) now report a Clark-West /
  Diebold-Mariano style test instead of a bare RMSE delta.
* **Multiplicity control.** Benjamini-Hochberg FDR adjustment is applied to the
  per-series unit-root p-values and to the coupling-metric p-values.

## Reproducibility

* All stochastic routines now honour their `seed` argument: `rotation_null_test()`,
  `test_permutation_robustness()`, `sensitivity_analysis_weights()` and the noise
  injection in `diagnose_data()` previously ignored it.

## Bug fixes

* `zzz.R` defined `.onAttach()` twice; the Stan-availability message is restored.
* `.pkgenv` was overwritten with a vector, destroying the package environment; it
  is now a real environment created in `.onLoad()`.
* `read_cpi()` used `which.max()` on a logical, which silently selected column 1
  when no date/CPI column matched; it now errors cleanly.
* `run_convergence_robustness_tests()` used `exists("indices$sensitivity")`
  (always `FALSE`); the sensitivity block is now propagated.
* `to_num_commas()` now parses both European and Anglo-Saxon number formats.
* `calculate_hc_manual()` implements HC4 (it previously fell through to HC0) and
  uses `sandwich`/`lmtest` when available.
* `estimate_DFM()` now returns robust t-statistics and p-values, reports the
  dominant-root half-life, and no longer floors the out-of-sample R^2 at -2.
* `visualize_factor_dynamics()` no longer fails when `save_plot = TRUE` without a
  file name, and graphics are opt-in within the pipeline.
* `test_reweighting_robustness()` no longer references an unassigned object when
  weight reading fails.

## Documentation

* `prepare_marxist_factors()` now actually incorporates the auxiliary series
  (`TMG`, `COM_matrix`, `SPVR_matrix`, `CA`) it accepts.
* README "Basic Usage" used non-existent functions (`prepare_panel_data`,
  `estimate_dfm`, `estimate_ou_process`, `plot_convergence`); it now uses the
  real API.
* `@return` sections of the test functions, `select_optimal_components_safe()`
  and `estimate_DFM()` now match the objects actually returned.
* New vignette section "Methodological notes and design decisions" documents the
  discrete-time OU/AR(1) interpretation, PLS vs PCA factor extraction, the
  generated-regressor caveat in residual-based unit-root tests, and the
  weighting-blend ("disaggregation") method.
* `DESCRIPTION` no longer over-claims `sandwich`/`lmtest` usage; both are now
  Suggests and used when present, with a documented HC fallback.

# convergenceDFM 0.1.4

* Initial package structure.
