#### ---- internal functions ---- ####
#' Small helper to convert R objects to strings.
#'
#' Internal helper function. This function allows to convert expressions to strings and collapse them.
#'
#' @param x A vector of \code{R} objects to be converted into strings.
#'
#' @noRd
.collapse_deparse <- function(x) {
  txt <- paste(deparse(x), collapse = " ")
  txt <- gsub("\\s+", " ", txt) # replace multiple spaces by single spaces
  trimws(txt)
}

#' Small helper to update name in call
#'
#' Internal helper function. This function allows to modify the data name in a model call for linear models and linear mixed models.
#'
#' @param model An \code{S4} object of class \code{lmerMod} or an \code{S3} object of class \code{lm}.
#' @param data_name_symbol A \code{name} or \code{symbol} representing the dataset to insert in the model call.
#'
#' @noRd
.update_call_data <- function(model, data_name_symbol) {

  if (!inherits(model, c("lm", "lmerMod"))) {
    stop("model must be of class 'lm' or 'lmerMod'")
  }

  # ---- function ----

  if (inherits(model, "lm")) {
    model$call$data <- data_name_symbol
  } else {
    model@call$data <- data_name_symbol
  }

  return(model)
}

#' Small helper to set the model evaluation environment
#'
#' Internal helper function. This function allows to modify the evaluation environment of a linear (mixed) model.
#'
#' @param model An \code{S4} object of class \code{lmerMod} or an \code{S3} object of class \code{lm}.
#' @param data_name_symbol An \code{environment}.
#'
#' @noRd
.set_call_env <- function(model, env) {

  if (!inherits(model, c("lm", "lmerMod"))) {
    stop("model must be of class 'lm' or 'lmerMod'")
  }

  # ---- function ----

  if (inherits(model, "lm")) {
    environment(model$terms) <- env
    environment(model$call) <- env

  } else {
    environment(model@call) <- env
  }
  return(model)
}

#' Small helper to set controls in 'lmerMod' class object
#'
#' Internal helper object. Sets the controls for \code{lmerMod} class objects.
#' Uses Bound Optimization BY Quadratic Approximation, "bobyqa" and 1e5 as
#' maximal number of evaluation.
#'
#' @noRd
lmerControl_default <- lmerControl(optimizer = "bobyqa",
                                   optCtrl = list(maxfun = 1e5))
