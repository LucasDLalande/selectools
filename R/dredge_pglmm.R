#' Model selection for pglmm
#'
#' This is an adaptation of the \code{dredge} function from the \code{MuMIn} package to
#' phylogenetic generalized linear mixed models (\code{pglmm} from the \code{phyr}). It
#' generates a model selection table of models as would do the \code{dredge}
#' function from \code{MuMIn} with combinations (subsets) of fixed effect terms in
#' the global model, with optional model inclusion rules.
#'
#' @param formulaRE A character string specifying the random effects structure (e.g., "~ 1 + (1 | Species)"). Must include an intercept (`1`) as fixed effect.
#' @param fixed A character vector of fixed-effect names to test. Interaction terms are allowed (e.g., "Sepal.Width:Petal.Length"), but some post-processing of tested models might be required.
#' @param data A data.frame containing the variables named in formula.
#' @param rank Criterion used to rank the models. Choose either "AIC" (default) or "AICc".
#' @param family Distribution family to use in model fitting. Options are "gaussian", "binomial", or "poisson".
#' @param cov_ranef A named list of covariance matrices of random terms. The names should be the group variables that are used as random terms with specified covariance matrices (without the two underscores, e.g. list(sp = tree1, site = tree2)). The actual object can be either a phylogeny with class "phylo" or a prepared covariance matrix. If it is a phylogeny, pglmm will prune it and then convert it to a covariance matrix assuming Brownian motion evolution. pglmm will also standardize all covariance matrices to have determinant of one. Group variables will be converted to factors and all covariance matrices will be rearranged so that rows and columns are in the same order as the levels of their corresponding group variables.
#' @param estimate Logical. If `TRUE` (default), includes fixed effect estimates in the output.
#' @param std.err Logical. If `TRUE`, includes standard errors for the fixed effects. Default to `FALSE`.
#' @param round Integer. Number of decimal places to round numeric results to. Default to `3`.
#'
#' @return A data frame summarizing model selection results, ranked by the chosen criterion.
#'
#' @export
#'
#' @examples
#' # Write the null model using pglmm() from 'phyr'
#' mod1 <- phyr::pglmm(log_onset ~ 1 + (1|Species__),
#'                   data = senescence, family = "gaussian", cov_ranef = list(Species = tree_ultra),
#'                   REML = TRUE, verbose = FALSE, s2.init = .1)
#'
#' # Run model selection based
#' dredge_pglmm(
#'        formulaRE = "log_onset ~ 1 + (1 | Species)",
#'        fixed = c("log_afr", "log_mass"),
#'        data = senescence,
#'        rank = "AICc",
#'        family = "gaussian",
#'        cov_ranef = list(Species = tree_ultra)
#'  )
dredge_pglmm <- function(formulaRE, fixed, data,
                         rank=c("AICc", "AIC"),
                         family=c("gaussian", "binomial", "poisson"),
                         cov_ranef,
                         estimate=T, std.err=F, round=3) {

  rank <- match.arg(rank)
  family <- match.arg(family)

  if(!is.data.frame(data)) {
    stop("data must be a dataframe")
  }

  if(!is.list(cov_ranef)) {
    stop("cov_ranef must be a named list of covariance matrices of random terms. For more information check the help page")
  }

  if (!all(vapply(cov_ranef, function(x) is.matrix(x) || inherits(x, "phylo"), logical(1)))) {
    stop("Each element of cov_ranef must be a matrix or a phylo object.")
  }

  if(!is.numeric(round) || round < 0 || round != floor(round)) {
    stop("round must be a non-negative integer")
  }

  null_mod <- .fit_null_model(formulaRE, data, family, cov_ranef)

  fixed_mod <- .fit_fixed_model(null_mod, fixed, data, family, cov_ranef)

  table <- .dredge_table(null_mod, fixed_mod, rank, estimate, std.err, round)

  data_expr <- substitute(data)
  .print_table(null_mod, fixed_mod, data, table, rank, cov_ranef, data_expr)

  invisible(table)
}

