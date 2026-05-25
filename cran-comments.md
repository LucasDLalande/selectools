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

- Warnings from `phyr::pglmm()` using deprecated `lme4::findbars()` and `lme4::nobars()`: 
upstream issue in `phyr`, not in `selectools::dredge_pglmm()`. 
The warning is displayed once per session.
