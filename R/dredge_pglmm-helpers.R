#' Fitting the null model with pglmm
#'
#' Internal helper function. This function fits the null model with pglmm from \pkg{phyr}
#' and extracts the degrees of freedom, log-likelihood, AICc and AIC.
#'
#' @param formulaRE A \code{formula} object specifying the random effects structure
#' (e.g., "log_onset ~ 1 + (1 | Species)"). Must include an intercept (\code{1}) as fixed effect.
#' @param data A \code{data.frame} containing the variables named in \code{formula}.
#' @param rank Criterion used to rank the models. Choose either \code{"AIC"} (default) or \code{"AICc"}.
#' @param family Distribution family to use in model fitting.
#' Options are \code{"gaussian"}, \code{"binomial"}, or \code{"poisson"}.
#' @param cov_ranef A named list of covariance matrices of random terms.
#' The names should be the group variables that are used as random terms with
#' specified covariance matrices (without the two underscores,
#' e.g. \code{list(sp = tree1, site = tree2)}).
#' The actual object can be either a phylogeny with class \code{phylo} or a prepared
#' covariance matrix. If it is a phylogeny, \code{\link[phyr]{pglmm}} will prune it and then convert
#' it to a covariance matrix assuming Brownian motion evolution.
#' \code{\link[phyr]{pglmm}} will also standardize all covariance matrices to have determinant of one.
#' Group variables will be converted to factors and all covariance matrices will
#' be rearranged so that rows and columns are in the same order as the levels of
#' their corresponding group variables.
#'
#' @return A list containing:
#' \describe{
#'    \item{formula.RE}{Random effects formula (class \code{formula}).}
#'    \item{model.RE}{Fitted \code{\code[phyr]{pglmm}} object.}
#'    \item{k.RE}{Number of parameters (class \code{numeric}).}
#'    \item{logLik.RE}{Log-likelihood(class \code{numeric}).}
#'    \item{rank.RE}{AICc or AIC value depending on rank argument (class \code{numeric}).}
#'    \item{data}{A \code{data.frame} containing the variables named in \code{formulaRE}.}
#'    \item{rank}{The criterion used to rank the models (\code{"AIC"} or \code{"AICc"}).}
#'    \item{family}{Distribution family used in model fitting (\code{"gaussian"}, \code{"binomial"}, or \code{"poisson"}).}
#'    \item{cov_ranef}{The named list of covariance matrices of random terms}
#'    }
#'
#' @keywords internal
#'
#' @noRd
.fit_null_model <- function(formulaRE, data, rank, family, cov_ranef) {

  # runs RE model
  mod.RE <- pglmm(formulaRE, data, family,
                  cov_ranef = cov_ranef,
                  REML = FALSE, verbose = FALSE, s2.init = .1)

  # Extracts df, log-likelihood and computes AICc and AIC
  k.RE <- nrow(phyr::fixef(mod.RE))+nrow(phyr::ranef(mod.RE))
  logLik.RE <- mod.RE$logLik

  if (rank == "AICc") {

    rank.RE <- mod.RE$AIC + (2 * k.RE^2 + 2 * k.RE) / (nrow(data) - k.RE - 1)

  } else if (rank == "AIC") {

    rank.RE <- mod.RE$AIC

  }

  return(list(
    formula.RE = formulaRE,
    model.RE = mod.RE,
    k.RE = k.RE,
    logLik.RE = logLik.RE,
    rank.RE = rank.RE,
    data = data,
    rank = rank,
    family = family,
    cov_ranef = cov_ranef
  ))
}

#' Fitting the fixed-effect models with pglmm
#'
#' Internal helper function. This function fits all possible fixed-effect models
#' with \code{\link[phyr]{pglmm}} and extracts the degrees of freedom, log-likelihood, AICc and AIC.
#'
#' @param null_mod An object containing the results of the internal
#' \code{.fit_null_model()} function.
#' @param fixed A character vector of fixed effects to test.
#' Interaction terms are allowed (e.g., \code{"Sepal.Width:Petal.Length"}),
#' but some post-processing of tested models might be required.
#'
#' @return A list containing:
#' \describe{
#'    \item{model_list.fix}{A list of fitted \code{pglmm} object.}
#'    \item{k.fix}{A vector of number of parameters (class \code{numeric}).}
#'    \item{logLik.fix}{A vector of Log-likelihood(class \code{numeric}).}
#'    \item{rank.fix}{A vector of AICc or AIC values depending on rank argument (class \code{numeric}).}
#'    \item{rank_name}{A \code{character} string specifying the name of the rank criterion chosen.}
#'    }
#'
#' @keywords internal
#'
#' @noRd
.fit_fixed_model <- function(null_mod, fixed) {

  data <- null_mod$data
  rank <- null_mod$rank
  family <- null_mod$family
  cov_ranef <- null_mod$cov_ranef

  # creates all possible combination of fixed effects and store it as a simple list
  comb.fix <- unlist(Map(combn, list(fixed), seq_along(fixed), simplify = FALSE), recursive = FALSE)

  # puts the fixed formulas alongside the random effect formula for each combination of fixed effect formula
  form.fix <- lapply(comb.fix, function(vars){ # "vars" is each element of the comb.fix successively
    fixed_part <- paste(vars, collapse = "+")
    update(null_mod$formula.RE, as.formula(paste("~", fixed_part, "+ .")))
  })

  # fits the models, stored into list
  model_list.fix <- lapply(form.fix, function(f) {
    pglmm(f, data, family,
          cov_ranef,
          REML = FALSE, verbose = FALSE, s2.init = .1)
  })

  # extracts df of each model
  k.fix <- vapply(model_list.fix, function(m) {

    nrow(phyr::fixef(m)) + nrow(phyr::ranef(m))

  }, numeric(1))

  # extracts log-likelihood of each model
  logLik.fix <- vapply(model_list.fix, function(m) {

    m$logLik

  }, numeric(1))

  if (rank == "AICc") {
    # extracts AICc of the null model
    rank.fix <- vapply(model_list.fix, function(m) {

      k <- nrow(phyr::fixef(m)) + nrow(phyr::ranef(m))
      n <- nrow(data)
      aic <- m$AIC

      aic + (2 * k^2 + 2 * k) / (n - k - 1)

    }, numeric(1))
  } else if (rank == "AIC") {
    # extracts AIC values for fixed-effects models
    rank.fix <- vapply(model_list.fix, function(m) {

      m$AIC

    }, numeric(1))
  }

return(list(
  model_list.fix = model_list.fix,
  k.fix = k.fix,
  logLik.fix = logLik.fix,
  rank.fix = rank.fix,
  rank_name = rank
))
}

#' Building the dredge table
#'
#' Internal helper function. This function assembles all models tested with the
#' fixed-effects estimate values and standard-error, their number of parameters,
#' AIC(c), loglikelihood and computed delta AIC(c) and weights. Models are order
#' from the smallest delta AIC(c) to the greatest.
#'
#' @param null_mod An object containing the results of the internal
#' \code{.fit_null_model()} function.
#' @param fixed_mod An object containing the results of the internal
#' \code{.fit_fixed_model()} function.
#' @param estimate Logical. If \code{TRUE} (default), includes fixed effect estimates
#' in the output.
#' @param std.err Logical. If \code{TRUE}, includes standard errors for the fixed effects.
#' Defaults to \code{FALSE}.
#' @param round Integer. Number of decimal places to round numeric results to. Defaults to \code{3}.
#'
#' @return  A list containing:
#' \describe{
#'    \item{dredge_table}{A \code{data.frame} object with ranked models according to AICc and AIC}
#'    \item{full_model_chr}{A \code{character} string displaying the formula of the full model}
#'    \item{random_display}{A \code{character} string of random terms}
#'    \item{rank_name}{A \code{character} string specifying the name of the rank criterion chosen}
#'    }
#'
#' @keywords internal
#'
#' @noRd
.dredge_table <- function(null_mod, fixed_mod, estimate, std.err, round) {

  rank_name <- fixed_mod$rank_name
  data_name <- fixed_mod$data_name
  cov_ranef_name <- fixed_mod$cov_ranef_name

  # extracts the full model as character
  full_model_chr <- as.character(as.expression(fixed_mod$model_list.fix[[length(fixed_mod$model_list.fix)]]$formula_original))

  # extracts the random part of the model
  random_terms <- findbars(null_mod$formula.RE)
  random_display <- paste(paste0("(", vapply(random_terms, deparse, character(1)), ")"), collapse = " + ")

  # merge lists for the null and fixed models, dfs, AICcs, AICs and log-likelihoods
  model_list <- c(fixed_mod$model_list.fix, list(null_mod$model.RE))
  k <- c(fixed_mod$k.fix, null_mod$k.RE)
  rank <- c(fixed_mod$rank.fix, null_mod$rank.RE)
  logLik <- c(fixed_mod$logLik.fix, null_mod$logLik.RE)

  # gets the unique fixed effects from all models
  all_fixed_effects <- unique(unlist(lapply(model_list, function(m) rownames(phyr::fixef(m)))))

  # Creates an empty dataframe with columns named after the unique fixed effects
  dredge_table <- data.frame(matrix(ncol = length(all_fixed_effects), nrow = 0))
  colnames(dredge_table) <- all_fixed_effects

  # Iterates through each model and adds a row to the dataframe with fixed effects values
  for (i in seq_along(model_list)) {
    model_fixef_names <- rownames(phyr::fixef(model_list[[i]])) # extracts fixed effects name
    model_values <- as.data.frame(t(phyr::fixef(model_list[[i]])[,1])) # extracts estimates of fixed effects
    model_se <- as.data.frame(t(phyr::fixef(model_list[[i]])[,2])) # extracts standard-errors of fixed effect

    dredge_table[i, ] <- NA

    if (estimate == TRUE & std.err == TRUE) {
      # fills in the values and se for fixed effects present in the model
      dredge_table[i, model_fixef_names] <- paste(round(model_values, round), paste0("(", round(model_se, round), ")"), sep = " ")

    } else {
      if (estimate == TRUE & std.err == FALSE) {
        # fills in the values fixed effects present in the model
        dredge_table[i, model_fixef_names] <- round(model_values, round)

      } else {
        # fills in with "+" instead of estimates values
        dredge_table[i, model_fixef_names] <- "+"

      }
    }
  }

  # adds corresponding k, logLiks and rank to the dataframe and order dataframe
  dredge_table <- cbind(dredge_table, k, logLik, rank)
  dredge_table <- dredge_table[order(dredge_table$rank),]

  # calculates the delta and weight columns
  dredge_table$delta <- dredge_table$rank - min(dredge_table$rank)
  dredge_table$weight <- exp(-0.5 * dredge_table$delta) / sum(exp(-0.5 * dredge_table$delta))

  colnames(dredge_table)[colnames(dredge_table) == "rank"] <- rank_name

  return(list(
    dredge_table = dredge_table,
    full_model_chr = full_model_chr,
    random_display = random_display,
    rank_name = rank_name
    ))
}




