#### ---- internal functions ---- ####
#' Creating a copy of \code{data} to be modified.
#'
#' Internal helper function. This function allows to prepare a copy of the original
#' dataset and the variables to be modified according to threshold value for
#' threshold models.
#'
#' @param trajectory_trait A \code{character} string representing the trajectory trait.
#' Must exist in \code{data}.
#' @param data A \code{data.frame}.
#'
#' @return A list of objects including a copy of \code{data} named \code{data2} and \code{var1}
#' and \code{var2} named after \code{trajectory_trait}.
#'
#' @noRd
.make_modification_data <- function(trajectory_trait, data) {

  return(list(
    data2 = data,
    var1 = paste0(trajectory_trait, ".1"),
    var2 = paste0(trajectory_trait, ".2")
  ))
}

#' Adjusting trajectory trait for threshold models (constant-slope).
#'
#' Internal helper function. This function allows to adjust the trajectory trait
#' according to the threshold value for the constant-slope threshold trajectory.
#'
#' @param data_tmp A \code{data.frame} object obtained from original dataset.
#' @param var1 A numeric vector identical to the variable \code{trajectory_trait},
#' to be modified according to threshold value \code{thr}. Here, \code{var1} is unused but
#' required for \code{.adjust_data_SC} and \code{.adjust_data_SS} and kept for
#' consistency across \code{adjust_functions} in \code{\link{estimate_threshold}}.
#' @param var2 A \code{numeric} vector identical to the variable \code{trajectory_trait},
#' to be modified according to threshold value \code{thr}.
#' @param trajectory_trait A \code{character} string representing the trajectory trait.
#' Must exist in \code{data_tmp}.
#' @param thr A \code{numeric} value representing the threshold value tested.
#'
#' @return A \code{data.frame} object similar to \code{data_tmp} with modified \code{var2} variable.
#'
#' @noRd
.adjust_data_CS <- function(data_tmp, var1, var2, trajectory_trait, thr) {

  # ---- var1 and var2 ----
  if (!is.character(var1) || length(var1) != 1) {
    stop("var1 must be a single character string")
  }

  if (!is.character(var2) || length(var2) != 1) {
    stop("var2 must be a single character string")
  }

  # ---- thr ----
  if (!is.numeric(thr)) {
    stop("thr must be a numeric value")
  }

  # ---- function ----

  data_tmp[[var2]] <- ifelse(data_tmp[[trajectory_trait]] < thr, 0, data_tmp[[trajectory_trait]] - thr)

  return(data_tmp)
}

#' Adjusting trajectory trait for threshold models (slope-constant).
#'
#' Internal helper function. This function allows to adjust the trajectory trait
#' according to the threshold value for the slope-constant threshold trajectory.
#'
#' @param data_tmp A \code{data.frame} object obtained from original dataset.
#' @param var1 A numeric vector identical to the variable \code{trajectory_trait},
#' to be modified according to threshold value \code{thr}.
#' @param var2 A \code{numeric} vector identical to the variable \code{trajectory_trait},
#' to be modified according to threshold value \code{thr}. Here, \code{var1} is unused but
#' required for \code{.adjust_data_SC} and \code{.adjust_data_SS} and kept for
#' consistency across \code{adjust_functions} in \code{\link{estimate_threshold}}.
#' @param trajectory_trait A \code{character} string representing the trajectory trait.
#' Must exist in \code{data_tmp}.
#' @param thr A \code{numeric} value representing the threshold value tested.
#'
#' @return A \code{data.frame} object similar to \code{data_tmp} with modified \code{var1} variable.
#'
#' @noRd
.adjust_data_SC <- function(data_tmp, var1, var2, trajectory_trait, thr) {

  # ---- var1 and var2 ----
  if (!is.character(var1) || length(var1) != 1) {
    stop("var1 must be a single character string")
  }

  if (!is.character(var2) || length(var2) != 1) {
    stop("var2 must be a single character string")
  }

  # ---- thr ----
  if (!is.numeric(thr)) {
    stop("thr must be a numeric value")
  }

  # ---- function ----

  data_tmp[[var1]] <- ifelse(data_tmp[[trajectory_trait]] < thr, data_tmp[[trajectory_trait]] - thr, 0)

  return(data_tmp)
}

#' Adjusting trajectory trait for threshold models (slope-slope).
#'
#' Internal helper function. This function allows to adjust the trajectory trait
#' according to the threshold value for the slope-slope threshold trajectory.
#'
#' @param data_tmp A \code{data.frame} object obtained from original dataset.
#' @param var1 A numeric vector identical to the variable \code{trajectory_trait},
#' to be modified according to threshold value \code{thr}.
#' @param var2 A \code{numeric} vector identical to the variable \code{trajectory_trait},
#' to be modified according to threshold value \code{thr}.
#' @param trajectory_trait A \code{character} string representing the trajectory trait.
#' Must exist in \code{data_tmp}.
#' @param thr A \code{numeric} value representing the threshold value tested.
#'
#' @return A \code{data.frame} object similar to \code{data_tmp} with modified \code{var1}
#' and \code{var2} variable.
#'
#' @noRd
.adjust_data_SS <- function(data_tmp, var1, var2, trajectory_trait, thr) {

  # ---- var1 and var2 ----
  if (!is.character(var1) || length(var1) != 1) {
    stop("var1 must be a single character string")
  }

  if (!is.character(var2) || length(var2) != 1) {
    stop("var2 must be a single character string")
  }

  # ---- thr ----
  if (!is.numeric(thr)) {
    stop("thr must be a numeric value")
  }

  # ---- function ----

  data_tmp[[var1]] <- ifelse(data_tmp[[trajectory_trait]] > thr, thr, data_tmp[[trajectory_trait]])
  data_tmp[[var2]] <- ifelse(data_tmp[[trajectory_trait]] <= thr, 0, data_tmp[[trajectory_trait]] - thr)

  return(data_tmp)
}

#' Fitting threshold models for each threshold value
#'
#' Internal function. This function allows to fit threshold models based on
#' threshold formulas constructed with \code{\link{trajectory_formulas}} based on the
#' modified trajectory trait variables for each threshold trajectory as done by
#' \code{.adjust_data_CS}, \code{.adjust_data_SC} and \code{.adjust_data_SS}.
#' Models are fitted as linear model \code{\link[stats]{lm}} or linear mixed models
#' using \code{\link[lme4]{lmer}} depending on whether random effects are absent or present, respectively.
#'
#' @param data_tmp A \code{data.frame} object obtained from original dataset and
#' containing the variables named in \code{threshold_formulas}.
#' @param threshold_formulas A \code{list} of object of class \code{formula}
#' representing the threshold trajectory formulas or a \code{trajectory.formulas}
#' object generated by \code{\link{trajectory_formulas}}.
#' @param random A list of the random effects present in \code{threshold_formulas} obtained by
#' \code{\link[reformulas]{findbars}}.
#' @param lmer_control Sets the controls for \code{lmerMod} class objects. Defaults
#' to internal object \code{lmerControl_default} which uses Bound Optimization BY
#' Quadratic Approximation, "bobyqa" and 1e5 as maximal number of evaluation.
#'
#' @return A \code{list} of \code{S3} or \code{S4} objects.
#'
#' @noRd
.refit_threshold_models <- function(data_tmp, threshold_formulas, random = NULL, lmer_control) {

  # ---- function ----

  if (!is.null(random)) {
    model_refit <- lmer(threshold_formulas, data = data_tmp, control = lmer_control, REML = FALSE)
    model_refit@call[[1]] <- quote(lme4::lmer)
  } else {
    model_refit <- lm(threshold_formulas, data = data_tmp)
    model_refit$call[[1]] <- quote(stats::lm)
  }

  return(model_refit)
}

#### ---- API function ---- ####
#' Estimating threshold for each threshold trajectories
#'
#' This function estimates threshold value for each threshold trajectory based on
#' AIC or AICc values. For each threshold trajectory, \code{.make_modification_data}
#' creates a temporary dataset as a copy of \code{data} and including \code{var1} and \code{var2},
#' copies of \code{trajectory_trait}. \code{var1} and \code{var2} are iteratively recalculated
#' with the function list \code{adjust_functions} along a defined sequence of
#' threshold values. \code{adjust_functions} is a list of \code{.adjust_data_CS},
#' \code{.adjust_data_SC} and \code{.adjust_data_SS} functions that recalculates
#' \code{var1} and \code{var2} according to the threshold value.
#' Threshold estimations are performed following methods used in Douhard et al. (2017),
#' Briga et al. (2019), or Lalande et al. (2024).
#'
#' @param data A \code{data.frame} object containing the variables named in \code{threshold_formulas}.
#' @param threshold_formulas A \code{list} of object of class \code{formula}
#' representing the threshold trajectory formulas or a \code{trajectory.formulas}
#' object generated by \code{\link{trajectory_formulas}}.
#' @param trajectory_trait A \code{character} string representing the trajectory trait.
#' Must exist in \code{data}.
#' @param start,end,iter \code{Numeric} values controlling the threshold range search:
#'   \itemize{
#'     \item \code{start}: Lower bound of the threshold range.
#'     \item \code{end}: Upper bound of the threshold range.
#'     \item \code{iter}: Increment of the sequence.
#'     }
#' @param rank Criterion used to rank the models. Choose either \code{"AIC"} or \code{"AICc"}.
#' @param lmer_control Advanced users can supply custom optimizer controls for
#' \code{lmerMod} class objects. Defaults to \code{NULL}, meaning that \code{\link[lme4]{lmer}} uses the
#' internal object \code{lmerControl_default} whith Bound Optimization BY Quadratic Approximation,
#' "bobyqa" and 1e5 as maximal number of evaluation.
#'
#' @return An object of class \code{threshold.estimation} containing:
#' \describe{
#'    \item{trajectory_trait}{A \code{character} string representing the trajectory trait.}
#'    \item{rank}{The criterion used to rank the models.}
#'    \item{thresholds}{A \code{data.frame} with criterion values for each threshold
#'    value and each threshold trajectory.}
#'    \item{best_thresholds}{A \code{list} of \code{numeric} values representing the threshold
#'    with the lowest criterion value for each threshold trajectory.}
#'    }
#'
#' @export
#'
#' @references
#' Briga, M., Jimeno, B. & Verhulst, S. (2019). Coupling lifespan and aging? The
#' age at onset of body mass decline associates positively with sex-specific lifespan
#' but negatively with environment-specific lifespan. \emph{Experimental Gerontology},
#' 119, 111-119. \doi{10.1016/j.exger.2019.01.030}
#'
#' Douhard, F., Gaillard, J.-M., Pellerin, M., Jacob, L. & Lemaître, J.-F. (2017).
#' The cost of growing large: Costs of post-weaning growth on body mass senescence
#' in a wild mammal. \emph{Oikos}, 126, 1329–1338. \doi{10.1111/oik.04421}
#'
#' Lalande, L. D., Bourgoin, G., Carbillet, J., Cheynel, L., Debias, F., Ferté, H.,
#' Gaillard, J., Garcia, R., Lemaître, J., Palme, R., Pellerin, M., Peroz, C., Rey, B.,
#' Vuarin, P. and Gilot-Fromont, E. (2024). Early-life glucocorticoids accelerate
#' lymphocyte count senescence in roe deer. \emph{General and Comparative Endocrinology},
#' 357, 114595. \doi{10.1016/j.ygcen.2024.114595}
#'
#' @examples
#' # Provide an object resulting from trajectory_formulas() and containing threshold formulas
#' formulas <- trajectory_formulas(formula = log_onset ~ log_mass + (1|Species),
#'                                 trajectory_trait = "log_afr", interactions = "log_mass",
#'                                 data = senescence,
#'                                 trajectories = c("linear", "quadratic", "thresholdCS"))
#'
#' # Pass an object of class "trajectory_formulas" to estimate_threshold()
#' estimate_threshold(data = senescence, threshold_formulas = formulas,
#'                    trajectory_trait = "log_afr", start = -0.5, end = 5.7, iter = 0.1,
#'                    rank = "AICc")
#'
#' # Or pass a named list to estimate_threshold()
#' # Names must be chosen from "thresholdCS", "thresholdSC" and "thresholdSS"
#' # The trajectory trait must have the suffix ".2" if the slope
#' # is after the threshold, and ".1" if the slope is before the threshold.
#' formulas <- list(
#'   thresholdCS = log_onset ~ log_mass + log_afr.2 + log_mass:log_afr.2 + (1|Species),
#'   thresholdSC = log_onset ~ log_mass + log_afr.1 + log_mass:log_afr.1 + (1|Species),
#'   thresholdSS = log_onset ~ log_mass + log_afr.1 + log_afr.2 +
#'                 log_mass:log_afr.1 + log_mass:log_afr.2 + (1|Species))
#'
#' estimate_threshold(data = senescence, threshold_formulas = formulas,
#'                    trajectory_trait = "log_afr", start = -0.5, end = 5.7, iter = 0.1,
#'                    rank = "AICc")
estimate_threshold <- function(data,
                               threshold_formulas,
                               trajectory_trait,
                               start,
                               end,
                               iter,
                               rank = c("AICc", "AIC"),
                               lmer_control = NULL) {

  # ---- data ----
  if (missing(data) || is.null(data) || !is.data.frame(data)) {
    stop("data must be provided and must be an existing data.frame")
  }

  # ---- trajectory_trait ----
  if (!is.character(trajectory_trait) || length(trajectory_trait) != 1) {
    stop("trajectory_trait must be a single character string")
  }

  if (!trajectory_trait %in% names(data)) {
    stop("trajectory_trait must be a variable of 'data'")
  }

  # ---- threshold formulas ----
  if (!inherits(threshold_formulas, "trajectory.formulas")) {
    if (length(threshold_formulas) == 0) {
      stop("No valid threshold formulas provided")
    }
  }

  allowed <- c("thresholdCS", "thresholdSC", "thresholdSS")

  if (!inherits(threshold_formulas, "trajectory.formulas")) {

    if (is.null(names(threshold_formulas)) || !all(names(threshold_formulas) %in% allowed)) {
      stop("Unknown trajectory types in formulas. Names must be from: ", paste(allowed, collapse = ", "))
    }

    threshold_formulas <- threshold_formulas[names(threshold_formulas) %in% allowed]

  } else {
    non_threshold_formulas <- threshold_formulas$formulas$non_threshold_formulas
    threshold_formulas <- threshold_formulas$formulas$threshold_formulas
  }


  if (length(threshold_formulas) > 0) {
    if (!is.list(threshold_formulas) || !all(vapply(threshold_formulas, function(x) inherits(x, "formula"), logical(1)))) {
      stop("threshold_formulas must be a list of formula objects")
    }

    if (!all(vapply(threshold_formulas, function(f)
      length(f) == 3,
      logical(1)
    ))) {
      stop("threshold_formulas must be a list of two-sided formulas")
    }
  }


  types <- vapply(threshold_formulas, function(f) {
    attr(f, "trajectory_type") %||% NA_character_
  }, FUN.VALUE = character(1))

  if (any(is.na(types) | types == "")) {
    warning("Some threshold_formulas were not generated by trajectory_formulas(). Unexpected behaviors can be expected.")
  }

  # ---- sequence arguments ----
  if (!is.numeric(start) || length(start) != 1) {
    stop("start, end and iter must be single numeric values")
  }

  if (!is.numeric(end) || length(end) != 1) {
    stop("start, end and iter must be single numeric values")
  }

  if (!is.numeric(iter) || length(iter) != 1) {
    stop("start, end and iter must be single numeric values")
  }

  # ---- rank ----
  rank <- match.arg(rank)

  # ---- random ----
  random <- NULL
  if (length(threshold_formulas) > 0) {
    random <- findbars(threshold_formulas[[1]])
  } else {
    random <- findbars(non_threshold_formulas[[1]])
  }

  if (!is.null(random)) {
    random_txt <- paste0("(", vapply(random, function(x) .collapse_deparse(x), ""), ")", collapse = " + ")
  } else random_txt <- NULL

  # ---- lmer_control ----
  if (!is.null(lmer_control)) {
    if (!inherits(lmer_control, c("lmerControl", "merControl"))) {
      stop("lmer_control must be NULL or a lmerControl, merControl object. See help page for details")
    }
  }

  if (!is.null(random)) {
    if (is.null(lmer_control)) {
      lmer_control <- lmerControl_default
    }
  } else {
    lmer_control <- NULL
  }

  # ---- data adjusting functions ----
  adjust_functions <- list(
    thresholdCS = function(data_tmp, var1, var2, trajectory_trait, thr)
      .adjust_data_CS(data_tmp, var1, var2, trajectory_trait, thr),
    thresholdSC = function(data_tmp, var1, var2, trajectory_trait, thr)
      .adjust_data_SC(data_tmp, var1, var2, trajectory_trait, thr),
    thresholdSS = function(data_tmp, var1, var2, trajectory_trait, thr)
      .adjust_data_SS(data_tmp, var1, var2, trajectory_trait, thr)
  )

  # create the sequence on which we run threshold models
  sequence <- seq(start, end, by = iter)

  # create a dataframe in which we'll store rank values
  df_threshold <- data.frame(threshold = sequence)

  # create alternative data
  vars <- .make_modification_data(trajectory_trait, data)
  data2 <- vars$data2
  var1 <- vars$var1
  var2 <- vars$var2

  message(paste(names(threshold_formulas), collapse = ", "), " models will be fitted.")

  # modify alternative variable, refit models and store rank values across threshold values sequence
  for (i in seq_along(sequence)) {

    thr <- sequence[i]

    message("Fitting threshold models with breakpoint: ", thr, "\n")

    for (n in names(threshold_formulas)) {

      # reset data variables
      data_tmp <- data2
      data_tmp[[var1]] <- data_tmp[[trajectory_trait]]
      data_tmp[[var2]] <- data_tmp[[trajectory_trait]]

      # adjust variables
      data_tmp <- adjust_functions[[n]](data_tmp, var1, var2, trajectory_trait, thr)

      # refit models
      refitted_model <- .refit_threshold_models(data_tmp, threshold_formulas[[n]], random, lmer_control)

      # store AIC/AICc
      df_threshold[i, paste0(rank, "_", n)] <- if (rank == "AICc") AICc(refitted_model) else AIC(refitted_model)
    }
  }

  # ---- best thresholds ----

  best_thresholds <- setNames(vector("list", length(threshold_formulas)),
                              names(threshold_formulas))

  for (n in names(threshold_formulas)) {
    col_name <- paste0(rank, "_", n)
    best_thresholds[[n]] <- df_threshold$threshold[which.min(df_threshold[[col_name]])]
    attr(best_thresholds[[n]], "best_threshold") <- n
  }

  threshold_estimation <- list(
    trajectory_trait = trajectory_trait,
    rank = rank,
    random = random_txt,
    thresholds = df_threshold,
    best_thresholds = best_thresholds
    )

  class(threshold_estimation) <- "threshold.estimation"
  return(threshold_estimation)
}


#### ---- S3 Methods ---- ####
#' Print threshold estimation
#'
#' Displays a readable summary of the threshold estimations obtained by
#' \code{\link{estimate_threshold}}.
#'
#' @param x An object of class \code{threshold.estimation}.
#' @param ... Further arguments (currently unused).
#'
#' @return Invisibly returns \code{x}.
#'
#' @method print threshold.estimation
#' @seealso \code{\link{estimate_threshold}}
#'
#' @export
print.threshold.estimation <- function(x, ...) {
  .title_box("THRESHOLD ESTIMATION")
  cat("Treshold estimation for trait:", x$trajectory_trait, "\n")
  cat("Ranking criterion:", x$rank, "\n")
  cat("Random term(s):", ifelse(is.null(x$random), "None", x$random), "\n\n")
  .sep_line()

  cat("Best thresholds:\n")
  for (n in names(x$best_thresholds)) {
    cat("- ", n, ":", x$best_thresholds[[n]], "\n")
  }
  cat("\n")

  cat(.sep_line())
  cat("Use $thresholds to consult the data.frame with", x$rank, "values according to threshold.\n")
  cat("For a visual look use plot(x, model =", ifelse(length(x$best_thresholds) > 1, paste0("c('", paste(names(x$best_thresholds), collapse = "', '"), "')).", sep = ""), paste0("'", names(x$best_thresholds), "').")), "\n\n")

  .sep_line()

  return(invisible(x))
}

#' Plot threshold estimation
#'
#' Displays a visual representation of the threshold estimations obtained by
#' \code{\link{estimate_threshold}}.
#'
#' @param x An object of class \code{threshold.estimation}.
#' @param model A \code{character} string for the threshold trajectory results to be
#' displayed. Choose from \code{"thresholdCS"}, \code{"thresholdSC"} or \code{"thresholdSS"}.
#' @param ... Further arguments (currently unused).
#'
#' @return Invisibly returns \code{x}.
#'
#' @method plot threshold.estimation
#' @seealso \code{\link{estimate_threshold}}
#'
#' @export
plot.threshold.estimation <- function(x, model = c("thresholdCS", "thresholdSC", "thresholdSS"), ...) {

  model <- match.arg(model)

  y <- x$thresholds[[paste0(x$rank, "_", model)]]

  plot(x$thresholds$threshold, y,
       xlab = "Threshold",
       ylab = x$rank,
       main = paste0(x$rank, " value according to threshold\nfor the ", model, " model\n"),
       pch = 16,
       ylim = c(min(y)-3,max(y)+3),
       cex.axis = 0.85,
       cex.main = 1.5,
       font.lab = 2)
  abline(h = min(y), col = "navyblue", lwd = 2, lty = 2)
  abline(v = x$thresholds$threshold[y == min(y)], col = "navyblue", lwd = 2, lty = 2)
  abline(h = c(min(y)-2, min(y)+2), col = "lightcoral", lwd = 1, lty = 2)
  mtext(paste("Threshold with minimum", x$rank, "value (blue lines), +/- 2", x$rank, "(red lines)"), side = 3, cex = 0.85)

  return(invisible(x))
}
