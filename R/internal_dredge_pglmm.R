#' Fitting the null model with pglmm
#'
#' Internal helper function. This function fits the null model with pglmm from \code{phyr}
#' and extracts the degrees of freedom, log-likelihood, AICc and AIC.
#'
#' @param formulaRE A character string specifying the random effects structure (e.g., "~ 1 + (1 | Species)"). Must include an intercept (`1`) as fixed effect.
#' @param data A data.frame containing the variables named in formula.
#' @param family Distribution family to use in model fitting. Options are "gaussian", "binomial", or "poisson".
#' @param cov_ranef A named list of covariance matrices of random terms. The names should be the group variables that are used as random terms with specified covariance matrices (without the two underscores, e.g. list(sp = tree1, site = tree2)). The actual object can be either a phylogeny with class "phylo" or a prepared covariance matrix. If it is a phylogeny, pglmm will prune it and then convert it to a covariance matrix assuming Brownian motion evolution. pglmm will also standardize all covariance matrices to have determinant of one. Group variables will be converted to factors and all covariance matrices will be rearranged so that rows and columns are in the same order as the levels of their corresponding group variables.
#'
#' @return A list containing:
#' \describe{
#'    \item{formula.RE}{Random effects formula (class formula)}
#'    \item{model.RE}{Fitted pglmm object}
#'    \item{df.RE}{Degrees of freedom (class numeric)}
#'    \item{logLik.RE}{Log-likelihood(class numeric)}
#'    \item{AICc.RE}{AICc (class numeric)}
#'    \item{AIC.RE}{AIC (class numeric)}
#'    }
#'
#' @importFrom stats as.formula
#' @importFrom phyr pglmm
#' @importFrom lme4 fixef ranef
#'
#' @keywords internal
#'
#' @noRd
.fit_null_model <- function(formulaRE, data, family, cov_ranef) {

  family <- match.arg(family, c("gaussian", "binomial", "poisson"))

  if (!is.character(formulaRE)) {
    stop("formulaRE must be a character string")
  }

  # converts the character random effect (RE) formula into a formula
  formula.RE <- as.formula(formulaRE)

  # runs RE model
  mod.RE <- pglmm(formula.RE, data, family,
                        cov_ranef = cov_ranef,
                        REML=FALSE, verbose=FALSE, s2.init=.1)

  # Extracts df, log-likelihood and computes AICc and AIC
  df.RE <- nrow(fixef(mod.RE))+nrow(ranef(mod.RE))
  logLik.RE <- mod.RE$logLik
  AICc.RE <- mod.RE$AIC + ((2*df.RE^2)+2*df.RE)/(nrow(data)-df.RE-1)
  AIC.RE <- mod.RE$AIC

  return(list(
    formula.RE = formula.RE,
    model.RE = mod.RE,
    df.RE = df.RE,
    logLik.RE = logLik.RE,
    AICc.RE = AICc.RE,
    AIC.RE = AIC.RE
  ))
}

#' Fitting the fixed-effect models with pglmm
#'
#' Internal helper function. This function fits all possible fixed-effect models
#' with pglmm from \code{phyr} and extracts the degrees of freedom, log-likelihood, AICc and AIC.
#'
#' @param null_mod An object containing the results of the internal \code{.fit_null_model()} function.
#' @param fixed A character vector of fixed effects to test. Interaction terms are allowed (e.g., "Sepal.Width:Petal.Length"), but some post-processing of tested models might be required.
#' @param data A data.frame containing the variables named in formula.
#' @param family Distribution family to use in model fitting. Options are "gaussian", "binomial", or "poisson".
#' @param cov_ranef A named list of covariance matrices of random terms. The names should be the group variables that are used as random terms with specified covariance matrices (without the two underscores, e.g. list(sp = tree1, site = tree2)). The actual object can be either a phylogeny with class "phylo" or a prepared covariance matrix. If it is a phylogeny, pglmm will prune it and then convert it to a covariance matrix assuming Brownian motion evolution. pglmm will also standardize all covariance matrices to have determinant of one. Group variables will be converted to factors and all covariance matrices will be rearranged so that rows and columns are in the same order as the levels of their corresponding group variables.
#'
#' @return A list containing:
#' \describe{
#'    \item{model_list.fix}{A list of fitted pglmm object}
#'    \item{df.fix}{A vector of degrees of freedom (class numeric)}
#'    \item{logLik.fix}{A vector of Log-likelihood(class numeric)}
#'    \item{AICc.fix}{A vector of AICc (class numeric)}
#'    \item{AIC.fix}{A vector of AIC (class numeric)}
#'    }
#'
#' @importFrom utils combn
#' @importFrom stats update as.formula
#' @importFrom phyr pglmm
#' @importFrom lme4 fixef ranef
#'
#' @keywords internal
#'
#' @noRd
.fit_fixed_model <- function(null_mod, fixed, data, family, cov_ranef) {

  family <- match.arg(family, c("gaussian", "binomial", "poisson"))

  if (!is.character(fixed)) {
    stop("fixed must be a character string, or a vector of characters")
  }

  if (length(fixed) == 0) {
    stop("fixed must contain at least one variable")
  }

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
          REML=FALSE, verbose=FALSE, s2.init=.1)
  })

  # extracts df of each model
  df.fix <- vapply(model_list.fix, function(m) {
    nrow(fixef(m))+nrow(ranef(m))
  }, numeric(1))

  # extracts log-likelihood of each model
  logLik.fix <- vapply(model_list.fix, function(m) {
    m$logLik
  }, numeric(1))

  # extracts AICc of the null model
  AICc.fix <- vapply(model_list.fix, function(m) {
    m$AIC + ((2*(nrow(fixef(m))+nrow(ranef(m)))^2)+2*(nrow(fixef(m))+nrow(ranef(m))))/(nrow(data)-(nrow(fixef(m))+nrow(ranef(m)))-1)
  }, numeric(1))

  # extracts AIC values for fixed-effects models
  AIC.fix <- vapply(model_list.fix, function(m) {
    m$AIC
  }, numeric(1))

return(list(
  model_list.fix = model_list.fix,
  df.fix = df.fix,
  logLik.fix = logLik.fix,
  AICc.fix = AICc.fix,
  AIC.fix = AIC.fix
))
}

#' Building the dredge table
#'
#' Internal helper function. This function assembles all models tested with the
#' fixed-effects estimate values and standard-error, their degrees of freedom,
#' AIC(c), loglikelihood and computed delta AIC(c) and weights. Models are order
#' from the smallest delta AIC(c) to the greatest.
#'
#' @param null_mod An object containing the results of the internal \code{.fit_null_model()} function.
#' @param fixed_mod An object containing the results of the internal \code{.fit_fixed_model()} function.
#' @param rank Criterion used to rank the models. Choose either "AIC" (default) or "AICc".
#' @param estimate Logical. If `TRUE` (default), includes fixed effect estimates in the output.
#' @param std.err Logical. If `TRUE`, includes standard errors for the fixed effects. Default is `FALSE`.
#' @param round Integer. Number of decimal places to round numeric results to. Default is `3`.
#'
#' @return  A data frame summarizing model selection results, ranked by the chosen criterion.
#'
#' @importFrom lme4 fixef
#'
#' @keywords internal
#'
#' @noRd
.dredge_table <- function(null_mod, fixed_mod, rank, estimate, std.err, round) {

  rank <- match.arg(rank, c("AICc", "AIC"))

  # merge lists for the null and fixed models, dfs, AICcs, AICs and log-likelihoods
  model_list <- c(fixed_mod$model_list.fix, list(null_mod$model.RE))
  df <- c(fixed_mod$df.fix, null_mod$df.RE)
  AICc <- c(fixed_mod$AICc.fix, null_mod$AICc.RE)
  AIC <- c(fixed_mod$AIC.fix, null_mod$AIC.RE)
  logLik <- c(fixed_mod$logLik.fix, null_mod$logLik.RE)

  all_fixed_effects <- unique(unlist(lapply(model_list, function(m) rownames(fixef(m))))) # gets the unique fixed effects from all models

  # Creates an empty dataframe with columns named after the unique fixed effects
  dredge.table <- data.frame(matrix(ncol = length(all_fixed_effects), nrow = 0))
  colnames(dredge.table) <- all_fixed_effects

  # Iterates through each model and adds a row to the dataframe with fixed effects values
  for (i in seq_along(model_list)) {
    model_fixef_names <- rownames(fixef(model_list[[i]])) # extracts fixed effects name
    model_values <- as.data.frame(t(fixef(model_list[[i]])[,1])) # extracts values of fixed effects
    model_se <- as.data.frame(t(fixef(model_list[[i]])[,2])) # extracts standard-errors of fixed effect

    # Creates a template row with empty values
    template_row <- data.frame(matrix(NA, ncol = length(all_fixed_effects), nrow = 1))
    colnames(template_row) <- all_fixed_effects

    if (estimate==TRUE & std.err==TRUE){
      template_row[, model_fixef_names] <- paste(round(model_values,round), paste0("(", round(model_se,round), ")"), sep = " ") # fills in the values and se for fixed effects present in the model
    }else{
      if (estimate==TRUE & std.err==FALSE){
        template_row[, model_fixef_names] <- round(model_values,round) # fills in the values fixed effects present in the model
      }else{
        template_row[, model_fixef_names] <- "+"
      }
    }
    dredge.table <- rbind(dredge.table, template_row) # adds the row to the dataframe
  }

  if (rank=="AICc") {
    dredge.table <- cbind(dredge.table, df, logLik, AICc) # adds corresponding dfs, logLiks and AICcs to the dataframe

    dredge.table <- dredge.table[order(dredge.table$AICc),] # orders the dataframe according to AICc value
    dredge.table$delta_AICc <- dredge.table$AICc - min(dredge.table$AICc) # calculates the delta_AICc column
    dredge.table$weight <- exp(-0.5*dredge.table$delta_AICc)/sum(exp(-0.5*dredge.table$delta_AICc)) # calculates the AICc weight columns
  } else {
    dredge.table <- cbind(dredge.table, df, logLik, AIC) # adds corresponding dfs, logLiks and AICcs to the dataframe

    dredge.table <- dredge.table[order(dredge.table$AIC),] # orders the dataframe according to AIC value
    dredge.table$delta_AIC <- dredge.table$AIC - min(dredge.table$AIC) # calculates the delta_AIC column
    dredge.table$weight <- exp(-0.5*dredge.table$delta_AIC)/sum(exp(-0.5*dredge.table$delta_AIC)) # calculates the AIC weight columns
  }

  return(dredge.table)
}


#' Displaying the dredge table in the console
#'
#' Internal helper function. This function display a summary of model selection procedure
#' and display the dredge table in the console.
#'
#' @param null_mod An object containing the results of the internal \code{.fit_null_model()} function.
#' @param fixed_mod An object containing the results of the internal \code{.fit_fixed_model()} function.
#' @param data A data.frame containing the variables named in formula.
#' @param rank Criterion used to rank the models. Choose either "AIC" (default) or "AICc".
#' @param table A dataframe resulting from the internal \code{.dredge_table()} function.
#' @param cov_ranef A named list of covariance matrices of random terms. The names should be the group variables that are used as random terms with specified covariance matrices (without the two underscores, e.g. list(sp = tree1, site = tree2)). The actual object can be either a phylogeny with class "phylo" or a prepared covariance matrix. If it is a phylogeny, pglmm will prune it and then convert it to a covariance matrix assuming Brownian motion evolution. pglmm will also standardize all covariance matrices to have determinant of one. Group variables will be converted to factors and all covariance matrices will be rearranged so that rows and columns are in the same order as the levels of their corresponding group variables.
#'
#' @return  A data frame summarizing model selection results, ranked by the chosen criterion.
#'
#' @importFrom lme4 findbars
#'
#' @keywords internal
#'
#' @noRd
.print_table <- function(null_mod, fixed_mod, data, table, rank, cov_ranef, data_expr) {

  rank <- match.arg(rank, c("AICc", "AIC"))

  # extracts the full model as character
  full_model <- as.character(as.expression(fixed_mod$model_list.fix[[length(fixed_mod$model_list.fix)]]$formula_original))

  # extracts the random part of the model
  random_terms <- findbars(null_mod$formula.RE)
  random_display <- vapply(random_terms, deparse, character(1))
  # deparse tranforms an R expression into a character string
  # character(1) assures that each call to deparse returns exactly 1 character type element

  print(paste("Global model: pglmm(", full_model,")", sep=""))
  cat("Data:\n")
  cat(deparse(data_expr), "\n")

  cat("---\nModel selection table\n")
  print(table)
  if (rank=="AICc"){
    cat("Models ranked by AICc\nRandom terms (all models):\n", paste(random_display, collapse = " + "), sep="     ")
  }else{
    cat("Models ranked by AIC\nRandom terms (all models):\n", paste(random_display, collapse = " + "), sep="     ")
  }
  cat("\nCovariance matrix of random effects:\n")
  print(substitute(cov_ranef))
}
