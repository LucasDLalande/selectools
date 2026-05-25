## R CMD check results

0 errors | 0 warnings | 3 note

* This is a new release.

## Notes

- `CITATION.cff` at root level: required by GitHub and Zenodo for citation metadata.
The R citation file can be found at `inst/CITATION` as per R conventions.

- `attach()` usage in `dredge_trajectories()` and `best_trajectory()`: required to 
make threshold datasets available for `MuMIn::dredge()` and `MuMIn::get.models()` 
model re-evaluation. Always cleaned-up via `on.exit()` as recommended in `?attach()` 
'Good practice'.

## Warnings

- Warnings from `phyr::pglmm()` using deprecated `lme4::findbars()` and `lme4::nobars()`: 
upstream issue in `phyr`, not in `selectools::dredge_pglmm()`. 
The warning is displayed once per session.

## Examples

Examples for `estimate_threshold()`, `fit_models()`, `dredge_trajectories()`, 
`best_trajectory()` and `select_trajectory()` are wrapped in `\donttest{}` because 
they involve fitting linear mixed models with `lme4::lmer()` during threshold estimation 
on simulated data, which can produce numerical errors (e.g. "Downdated VtV is not 
positive definite") on certain platforms (observed on macOS ARM64). The examples 
are correct and run successfully on Windows and Linux.
