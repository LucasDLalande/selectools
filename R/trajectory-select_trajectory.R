#### ---- API function ---- ####
#' Wrapper for trajectory model selection
#'
#' This function allows to build desired trajectory formulas (constant, linear,
#' quadratic, factor and threshold) based on a base formula and a defined trait
#' of interest (the trajectory trait), to estimate threshold points for threshold
#' trajectories, to fit models and to apply \code{\link[MuMIn]{dredge}} on all trajectories,
#' retaining the best model trajectory.
#'
#' This function is a wrapper of \code{\link{trajectory_formulas}}, \code{\link{estimate_threshold}},
#' \code{\link{fit_models}}, \code{\link{dredge_trajectories}} and \code{\link{best_trajectory}}.
#' Please refer to the help page of the functions above for details.
#'
#' @param formula An object of class \code{formula} representing the constant trajectory formula.
#' @param trajectory_trait A character string representing the trajectory trait.
#' Must exist in \code{data}.
#' @param interactions A vector of \code{character} strings representing variables with
#' which \code{trajectory_trait} interacts. Must exist in \code{data}. Defaults to \code{NULL}.
#' @param data A \code{data.frame} object containing the variables named in \code{formula}.
#' @param trajectories A \code{character} vector representing the trajectories to be constructed.
#' Choose among \code{"constant"}, \code{"linear"}, \code{"quadratic"}, \code{"factor"}, \code{"thresholdCS"}, \code{"thresholdSC"}
#' and \code{"thresholdSS"}. Defaults to all trajectories (\code{"all"}). \code{"thresholdCS"} describes
#' a "constant-slope" trajectory where the slope is null before the threshold. Similarly,
#' \code{"thresholdSC"} describes a "slope-constant" trajectory where the slope is null
#' after the threshold, and \code{"thresholdSS"} authorize a non-null slope both before
#' and after the threshold.
#' @param start,end,iter \code{Numeric} values controlling the threshold range search:
#'   \itemize{
#'     \item \code{start}: Lower bound of the threshold range.
#'     \item \code{end}: Upper bound of the threshold range.
#'     \item \code{iter}: Increment of the sequence.
#'     }
#' @param rank Criterion used to rank the models. Choose either \code{"AIC"} or \code{"AICc"}.
#' @param delta_criterion A \code{numeric} value to define the criterion value within which
#' the best model is selected.
#' @param adjust_threshold Logical. If \code{TRUE} (default), adds 1 parameter (df) and
#' 2 AIC/AICc points to threshold models to account for the implicit cost of the
#' threshold parameter as suggested by some methods.
#' @param lmer_control Advanced users can supply custom optimizer controls for
#' \code{lmerMod} class objects. Defaults to \code{NULL}, meaning that \code{\link[lme4]{lmer}} uses the
#' internal object \code{lmerControl_default} whith Bound Optimization BY Quadratic Approximation,
#' "bobyqa" and 1e5 as maximal number of evaluation.
#'
#' @return An object of class \code{best.trajectory} containing:
#' \describe{
#'    \item{best_table}{A \code{data.frame} displaying the rank, number of parameters
#'    and delta rank value for the best model of each trajectory.}
#'    \item{final_traj}{A \code{data.frame} displaying the rank, number of parameters
#'    and delta rank value for the best model across all trajectories and according
#'    to the \code{delta_criterion} argument.}
#'    \item{final_traj_name}{A \code{character} string for the name of the final trajectory.}
#'    \item{summaries}{A \code{list} of \code{summary.merMod} objects (one per trajectory).}
#'    \item{summary}{A \code{list} of \code{summary.merMod} object for the final trajectory.}
#'    }
#'
#' @export
#'
#' @examples
#' # These examples are wrapped in \donttest{} because fitting models on the simulated
#' # data can produce platform-specific numerical errors (macOS ARM64).
#' \donttest{
#' res <- select_trajectory(formula = log_onset ~ log_mass + (1|Species),
#'                          trajectory_trait = "log_afr",
#'                          interactions = "log_mass",
#'                          data = senescence,
#'                          trajectories = "all",
#'                          start = -0.5, end = 5.7, iter = 0.1,
#'                          rank = "AICc",
#'                          delta_criterion = 2,
#'                          adjust_threshold = TRUE,
#'                          lmer_control = NULL)
#'
#' # display an overview of the results
#' res # type = "global" by default
#' print(res, type = "by_trajectory")
#'
#' # display summary (summaries) of retained trajectory (trajectories)
#' summary(res) # type = "global" by default
#' summary(res, type = "by_trajectory")
#' }
select_trajectory <- function(formula,
                              trajectory_trait,
                              interactions = NULL,
                              data,
                              trajectories = "all",
                              start, end, iter,
                              rank = c("AICc", "AIC"),
                              delta_criterion = 2,
                              adjust_threshold = TRUE,
                              lmer_control = NULL) {

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

  # ---- delta_criterion ----
  if (!is.numeric(delta_criterion) || length(delta_criterion) != 1 || is.na(delta_criterion) || delta_criterion < 0) {
    stop("delta_criterion must be a single non-negative numeric value")
  }

  # ---- pipeline ----
  message("\U0001F680 Running select_trajectory()...     \n\n")
  message("\u23F3 1/5 Building formulas...     ")
  res.formulas <- suppressMessages(
    trajectory_formulas(formula = formula,
                        trajectory_trait = trajectory_trait,
                        interactions = interactions,
                        data = data,
                        trajectories = trajectories))
  message("    \u2705 Done!\n")

  message("\u23F3 2/5 Estimating thresholds...     ")
  res.estimate <- suppressMessages(
    estimate_threshold(data = data,
                       threshold_formulas = res.formulas,
                       trajectory_trait = trajectory_trait,
                       start = start,
                       end = end,
                       iter = iter,
                       rank = rank,
                       lmer_control = lmer_control))
  if (length(res.formulas$formulas$threshold_formulas) > 0) message("    \u2728 Check!\n") else message("    \u23ED Threshold trajectories, step skipped.\n")

  message("\u23F3 3/5 Fitting models...     ")
  res.models <- suppressMessages(
    fit_models(formulas = res.formulas,
               trajectory_trait = trajectory_trait,
               data = data,
               best_thresholds = res.estimate,
               lmer_control = lmer_control))
  message("    \U0001F60A Done, keep goin'!\n")

  message("\u23F3 4/5 Dredging models...     ")
  res.dredge <- dredge_trajectories(models = res.models,
                                    trajectory_trait = trajectory_trait,
                                    rank = rank,
                                    delta_criterion = delta_criterion,
                                    data = data)
  message("    \U0001F340 Done, one more to go!\n")

  message("\u23F3 5/5 Selecting best trajectory...     ")
  res.best_trajectory <- suppressMessages(
    best_trajectory(dredge_results = res.dredge,
                    delta_criterion = delta_criterion,
                    adjust_threshold = adjust_threshold))
  message("    \U0001F973 ... and done!\n\n")

  message("Trajectories recap and models summaries can be consulted using:
print(x, type = c('global', 'by_trajectory'))
summary(x, type = c('global', 'by_trajectory'))\n
********************************************************************************
Note: For threshold models, df is incremented by 1 and ", rank, " by 2 to account for
the implicit threshold parameter.
********************************************************************************\n\n")

  message("\U0001F389 Pipeline completed! \U0001F485")

  return(res.best_trajectory)
}


