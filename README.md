
# selectools

<!-- badges: start -->
<!-- badges: end -->

The goal of `selectools` is to provide and facilitate model selection procedures, 
mostly based on AIC and AICc, and adapting the `dredge` function from the `MuMIn` 
package. `selectools` contains a unique (for now) function (`pglmm_dredge`) that can 
be used to perform model selection of phylogenetic generalized linear mixed 
models (`pglmm` from the `phyr` package). It generates a model selection table 
of models as would do the `dredge` function from `MuMIn` with combinations 
(subsets) of fixed effect terms in the global model, with optional model 
inclusion rules. 

## Installation

You can install the development version of `selectools` from [GitHub](https://github.com/) with:

``` r
install.packages("remotes")
remotes::install_url(
  "https://github.com/LucasDLalande/selectools/archive/refs/heads/master.zip"
)
```

## Usage
The package contains a unique `dredge_pglmm` function:

Define your random effects (`formulaRE`), fixed effects (`fixed`), data, 
model selection criterion (AIC or AICc), distribution family (gaussian , binomial or poisson), 
the covariance matrix of the random effects (`cov_ranef`, i.e. usually an object of class `phylo`)
and whether to display estimates and/or standard-errors.
``` r
dredge_pglmm(
  formulaRE,
  fixed,
  data,
  rank = c("AICc", "AIC"),
  family = c("gaussian", "binomial", "poisson"),
  cov_ranef,
  estimate = T,
  std.err = F,
  round = 3
)
```
It returns a data frame with ordered models and respective fixed effects structure (and estimates and/or standard-errors), based on AIC or AICc.

## Example

This is a basic example of the `dredge_pglmm` function using simulated package data 

``` r
# Write the null model using pglmm() from 'phyr'
mod1 <- pglmm(log_onset ~ 1 + (1|Species__),
             data = senescence, family = "gaussian", cov_ranef = list(Species = tree_ultra),
             REML = TRUE, verbose = FALSE, s2.init = .1)

# Returns a model selection table based on the null model and evaluating all possible models based on the fixed effect structure provided
dredge_pglmm(
        formulaRE = "log_onset ~ 1 + (1 | Species)",
        fixed = c("log_afr", "log_mass"),
        data = senescence,
        rank = "AICc",
        family = "gaussian",
        cov_ranef = list(Species = tree_ultra)
)
```
## Data documentation
`senescence` is a simulated dataset that provides age at first reproduction
(logarithm), onset of senescence (logarithm) and body mass (logarithm) for
33 species. See more details using:
``` r
help(senescence)
```

`tree_ultra` is an ultrametric phylogenetic tree (object of class `phylo`) 
representing evolutionary distances among 33 species. See more details using:
``` r
help(tree_ultra)
```

## Future development

`selectools` will shortly include a function designed for trajectory model selection. 
It will be particularly suited for the model selection of several trajectories as can be done in
senescence research, by automatically computing, testing, comparing and ranking 
constant, linear, quadratic or threshold trajectories. 

## Citation

To cite the `selectools` package in your publications, please use:

  Lalande LD (2025). _selectools: A facilitating model selection procedure_. R package version 0.1.0, 
  <https://github.com/LucasDLalande/selectools.git>.
