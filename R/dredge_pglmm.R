#' Model selection for pglmm
#'
#' This is an adaptation of \code{\link[MuMIn]{dredge}} to phylogenetic generalized
#' linear mixed models (\code{\link[phyr]{pglmm}}). It generates a model
#' selection table of models as would do \code{\link[MuMIn]{dredge}} with combinations
#' (subsets) of fixed effect terms in the global model, with optional model inclusion rules.
#'
#' @param formulaRE A \code{formula} object specifying the random effects structure
#' (e.g., \code{log_onset ~ 1 + (1 | Species)}). Must include an intercept (\code{1}) as fixed effect.
#' @param fixed A character vector of fixed effects to test.
#' Interaction terms are allowed (e.g., \code{"Sepal.Width:Petal.Length"}),
#' but some post-processing of tested models might be required.
#' @param data A \code{data.frame} containing the variables named in \code{formula}.
#' @param rank Criterion used to rank the models. Choose either \code{"AIC"} (default) or \code{"AICc"}.
#' @param family Distribution family to use in model fitting. Options are \code{"gaussian"},
#' \code{"binomial"}, or \code{"poisson"}.
#' @param cov_ranef A named list of covariance matrices of random terms. The names
#' should be the group variables that are used as random terms with specified covariance
#' matrices (without the two underscores, e.g. \code{list(sp = tree1, site = tree2)}).
#' The actual object can be either a phylogeny with class \code{phylo} or a prepared
#' covariance matrix. If it is a phylogeny, \code{\link[phyr]{pglmm}} will prune it and then convert
#' it to a covariance matrix assuming Brownian motion evolution. \code{\link[phyr]{pglmm}}
#' will also standardize all covariance matrices to have determinant of one.
#' Group variables will be converted to factors and all covariance matrices will
#' be rearranged so that rows and columns are in the same order as the levels of
#' their corresponding group variables.
#' @param estimate Logical. If \code{TRUE} (default), includes fixed effect estimates
#' in the output.
#' @param std.err Logical. If \code{TRUE}, includes standard errors for the fixed effects.
#' Defaults to \code{FALSE}.
#' @param round Integer. Number of decimal places to round numeric results to. Defaults to \code{3}.
#'
#' @return A \code{dredge.table} object containing:
#' \describe{
#'    \item{dredge_table}{A \code{data.frame} object with ranked models according to AICc and AIC.}
#'    \item{full_model_chr}{A \code{character} string displaying the formula of the full model.}
#'    \item{random_display}{A \code{character} string of random terms.}
#'    \item{rank_name}{A \code{character} string specifying the name of the rank criterion chosen.}
#'    \item{data_name}{A \code{character} string specifying the name of the \code{data}.}
#'    \item{cov_ranef_name}{A \code{character} string specifying the name of the covariance matrix of random terms.}
#'    \item{family_name}{A \code{character} string specifying the name of the
#'    family used in \code{pglmm} (\code{"gaussian"}, \code{"binomial"}, or \code{"poisson"}).}
#'    }
#'
#' @export
#'
#' @examples
#' # Run model selection based
#' res <- dredge_pglmm(formulaRE = log_onset ~ 1 + (1 | Species),
#'                     fixed = c("log_afr", "log_mass"),
#'                     data = senescence,
#'                     rank = "AICc",
#'                     family = "gaussian",
#'                     cov_ranef = list(Species = tree_ultra))
#' res
dredge_pglmm <- function(formulaRE, fixed, data,
                         rank=c("AICc", "AIC"),
                         family=c("gaussian", "binomial", "poisson"),
                         cov_ranef,
                         estimate=T, std.err=F, round=3) {


  # ---- data ----
  if (missing(data) || is.null(data) || !is.data.frame(data)) {
    stop("data must be provided and must be an existing data.frame")
  }

  data_name <- deparse(substitute(data))

  # ---- formulaRE ----
  if (!inherits(formulaRE, "formula")) {
    stop("formulaRE must be an object of class formula")
  }

  if (length(formulaRE) != 3) {
    stop("formulaRE must be a two-sided formula")
  }

  # ---- fixed ----
  if (!is.character(fixed)) {
    stop("fixed must be a character string, or a vector of characters")
  }

  if (length(fixed) == 0) {
    stop("fixed must contain at least one variable")
  }

  if (!all(fixed %in% names(data))) {
    stop("All variables in 'fixed' must be variables of 'data'")
  }

  # ---- rank ----
  rank <- match.arg(rank)

  # ---- family ----
  family <- match.arg(family)
  family_name <- family

  # ---- cov_ranef ----
  if (!all(vapply(cov_ranef, function(x) is.matrix(x) || inherits(x, "phylo"), logical(1)))) {
    stop("Each element of cov_ranef must be a matrix or a phylo object.")
  }

  cov_ranef_name <- deparse(substitute(cov_ranef))

  # ---- round ----
  if(!is.numeric(round) || round < 0 || round != floor(round)) {
    stop("round must be a non-negative integer")
  }

  # ---- function ----
  ## null model - random part
  null_mod <- .fit_null_model(formulaRE, data, rank, family, cov_ranef)

  ## fixed part
  fixed_mod <- .fit_fixed_model(null_mod, fixed)

  ## create dredhe table
  table <- .dredge_table(null_mod, fixed_mod, estimate, std.err, round)

  table$data_name <- data_name
  table$cov_ranef_name <- cov_ranef_name
  table$family_name <- family_name

  class(table) <- "dredge.table"

  return(table)
}


#### ---- S3 Methods ---- ####
#' Print dredge table
#'
#' Displays a readable dredge table of \code{\link[phyr]{pglmm}} model obtained by \code{\link{dredge_pglmm}}.
#'
#' @param x An object of class \code{dredge.table}.
#' @param ... Further arguments (currently unused).
#'
#' @return Invisibly returns \code{x}.
#'
#' @method print dredge.table
#' @seealso \code{\link{dredge_pglmm}}
#'
#' @export
print.dredge.table <- function(x, ...) {

  .title_box("DREDGE TABLE FOR PGLMM MODEL")

  cat(paste("Global model: pglmm(", x$full_model_chr,")\n", sep=""))
  cat("With family: ", x$family_name, "\n")
  cat("Data: ", x$data_name, "\n\n")

  .sep_line()
  cat("Model selection table\n\n")
  print(x$dredge_table)

  cat("\nModels ranked by", x$rank_name, "\n")
  cat("Random terms (all models): ", x$random_display, "\n")

  cat("Covariance matrix of random effects: ", x$cov_ranef_name, "\n\n")
  .sep_line()
}








