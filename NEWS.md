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
