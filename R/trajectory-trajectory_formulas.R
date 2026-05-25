#### ---- API function ---- ####
#' Building all trajectory formulas
#'
#' This function allows to build desired trajectory formulas (constant, linear,
#' quadratic, factor and threshold) based on a base formula and a defined trait
#' of interest (the trajectory trait).
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
#'
#' @return An object of class \code{trajectory.formulas} containing:
#' \describe{
#'    \item{trajectory_trait}{A \code{character} string representing the trajectory trait.}
#'    \item{interactions}{The interaction terms.}
#'    \item{random}{The random terms.}
#'    \item{data}{The name of the dataset.}
#'    \item{formulas}{A \code{list} of the trajectory formulas splitted into non-threshold
#'    and threshold formulas with a \code{trajectory_type} attribute.}
#'    \item{display}{A readable list of the formulas.}
#'    }
#'
#' @export
#'
#' @examples
#' # Provide a base formula, the trait of interest, its interactions, data and the trajectories desired
#' formulas <- trajectory_formulas(formula = log_onset ~ log_mass + (1|Species),
#'                                 trajectory_trait = "log_afr", interactions = "log_mass",
#'                                 data = senescence,
#'                                 trajectories = c("linear", "quadratic", "thresholdCS")
#'                                 )
trajectory_formulas <- function(formula, trajectory_trait, interactions = NULL, data, trajectories = "all") {

  # ---- data ----
  if (missing(data) || is.null(data) || !is.data.frame(data)) {
    stop("data must be provided and must be an existing data.frame")
  }

  data_name <- substitute(data)

  # ---- formula ----
  if (!inherits(formula, "formula")) {
    stop("formula must be an object of class formula")
  }

  if (length(formula) != 3) {
    stop("formula must be a two-sided formula")
  }

  # ---- trajectory_trait ----
  if (!is.character(trajectory_trait) || length(trajectory_trait) != 1) {
    stop("trajectory_trait must be a single character string")
  }

  if (!trajectory_trait %in% names(data)) {
    stop("trajectory_trait must be a variable of 'data'")
  }

  # ---- interactions ----
  if (!is.null(interactions)) {

    if (!is.character(interactions)) {
      stop("interactions must be NULL or a character vector")
    }

    if (any(is.na(interactions)) || any(interactions == "")) {
      stop("interactions must contain valid variable names")
    }

    if (!all(interactions %in% names(data))) {
      stop("some interactions are not variables in 'data'")
    }
  }

  # ---- trajectories ----
  allowed <- c("constant", "linear", "quadratic", "factor", "thresholdCS", "thresholdSC", "thresholdSS")
  allowed_non_threshold <- c("constant", "linear", "quadratic", "factor")
  allowed_threshold <- c("thresholdCS", "thresholdSC", "thresholdSS")

  if (!is.character(trajectories)) {
    stop("trajectories must be a character object")
  }

  if (length(trajectories) == 1) {
    if (trajectories == "all") {
      trajectories <- allowed
    } else if (trajectories == "non-threshold") {
      trajectories <- allowed_non_threshold
    } else if (trajectories == "threshold") {
      trajectories <- allowed_threshold
    }
  }

  if (!all(trajectories %in% allowed)) {
    stop ("Unknown trajectories: ", paste(setdiff(trajectories, allowed), collapse = ", "), "\nPlease choose among: ", paste(allowed, collapse = ", "), "\n
To select all non-threshold trajectories, choose 'non-threshold'.\n
To select all threshold trajectories, choose 'threshold'.\n
Defaults to 'all', selecting all trajectories.")
  }

  # ---- separate fixed part ----
  fixed <- nobars(formula)
  random <- findbars(formula)
  random_txt <- NULL

  # ---- construct models ----
  ## trajectory term
  terms <- list(
    linear = trajectory_trait,
    quadratic = paste(c(trajectory_trait, paste0("I(", trajectory_trait, "^2)")), collapse = " + "),
    factor = paste0("factor(", trajectory_trait, ")"),
    thresholdCS = paste0(trajectory_trait, ".2"),
    thresholdSC = paste0(trajectory_trait, ".1"),
    thresholdSS = paste(c(paste0(trajectory_trait, ".1"), paste0(trajectory_trait, ".2")), collapse = " + ")
  )

  ## construct fixed formulas without interactions and random effects
  formulas <- lapply(terms, function(term) {
    update(fixed, paste(". ~ . +", term))
  })

  all_formulas <- c(
    list(constant = fixed),
    formulas
  )

  ## models to be displayed
  formulas_display <- list(
    constant = all_formulas$constant,
    linear = all_formulas$linear,
    quadratic = all_formulas$quadratic,
    factor = all_formulas$factor,
    thresholdCS = all_formulas$thresholdCS,
    thresholdSC = all_formulas$thresholdSC,
    thresholdSS = all_formulas$thresholdSS
  )

  # ---- interactions ----
  if (!is.null(interactions)) {

    ## interaction terms
    interaction_terms <- list(
      linear = paste0(terms$linear, ":", interactions),
      quadratic = c(paste0(trajectory_trait, ":", interactions), paste0(paste0("I(", trajectory_trait, "^2)"), ":", interactions)),
      factor = paste0(terms$factor, ":", interactions),
      thresholdCS = paste0(terms$thresholdCS, ":", interactions),
      thresholdSC = paste0(terms$thresholdSC, ":", interactions),
      thresholdSS = c(paste0(paste0(trajectory_trait, ".1"), ":", interactions), paste0(paste0(trajectory_trait, ".2"), ":", interactions))
    )

    ## construct fixed formulas with interactions
    formulas <- mapply(FUN = function(f, t) {
      update(f, paste(". ~ . +", paste(t, collapse = " + ")))
    },
    f = formulas,
    t = interaction_terms,
    SIMPLIFY = FALSE)

    all_formulas <- c(
      list(constant = fixed),
      formulas
    )

    ## models to be displayed
    formulas_display <- list(
      constant = all_formulas$constant,
      linear = all_formulas$linear,
      quadratic = all_formulas$quadratic,
      factor = all_formulas$factor,
      thresholdCS = all_formulas$thresholdCS,
      thresholdSC = all_formulas$thresholdSC,
      thresholdSS = all_formulas$thresholdSS
    )
  }

  # ---- random effects ----
  if (!is.null(random)) {

    ## extract random part
    random_txt <- paste0("(", vapply(random, function(x) .collapse_deparse(x), ""), ")", collapse = " + ")

    ## display
    formulas_display <- lapply(all_formulas, function(f) {
      as.formula(paste(.collapse_deparse(f), "+", random_txt))
    })

    ## construct formulas with interactions and random effects
    all_formulas <- lapply(all_formulas, function(f) {
      update(f, paste(". ~ . +", random_txt))
    })
  }

  # ---- attributes ----
  for (n in names(all_formulas)) {
    attr(all_formulas[[n]], "trajectory_type") <- n
  }

  # ---- filter trajectories ----
  selected_display <- formulas_display[trajectories]

  selected_formulas <- all_formulas[trajectories]
  non_threshold_formulas <- selected_formulas[!grepl("^threshold", names(selected_formulas))]
  threshold_formulas <- selected_formulas[grepl("^threshold", names(selected_formulas))]

  trajectory_formulas <- list(
    trajectory_trait = trajectory_trait,
    interactions = interactions,
    random = random_txt,
    data = data_name,
    formulas = list(
      non_threshold_formulas = non_threshold_formulas,
      threshold_formulas = threshold_formulas),
    display = selected_display
  )

  class(trajectory_formulas) <- "trajectory.formulas"

  return(trajectory_formulas)
}

#### ---- S3 Method ---- ####
#' Print trajectory formulas
#'
#' Displays a readable summary of the trajectory formulas constructed by \code{\link{trajectory_formulas}},
#' including fixed effects, interactions and random effects when present.
#'
#' @param x An object of class \code{trajectory.formulas}.
#' @param ... Further arguments (currently unused).
#'
#' @return Invisibly returns \code{x}.
#'
#' @method print trajectory.formulas
#' @seealso \code{\link{trajectory_formulas}}
#'
#' @export
print.trajectory.formulas <- function(x, ...) {
  .title_box("TRAJECTORY FORMULAS")
  cat("Constructed trajectories based on trait:", x$trajectory_trait, "\n")
  cat("Interaction term(s):", paste(x$interactions, collapse = ", "), "\n")
  cat("Random term(s):", ifelse(is.null(x$random), "None", x$random), "\n")
  cat("Dataset:", x$data, "\n\n")
  .sep_line()

  cat("Trajectories fitted:\n\n")
  cat(toupper("Non-threshold trajectorie(s)\n"))
  if (length(x$formulas$non_threshold_formulas) > 0) {
    for (n in names(x$formulas$non_threshold_formulas)) {
      cat("- ", .first_upp(paste0(n, ":\n")))
      cat(.collapse_deparse(x$display[[n]]), "\n\n")
    }
  } else {
    cat("None\n\n")
  }

  cat(toupper("Threshold trajectorie(s)\n"))
  if (length(x$formulas$threshold_formulas) > 0) {
    for (n in names(x$formulas$threshold_formulas)) {
      cat("- ", .first_upp(paste0(n, ":\n")))
      cat(.collapse_deparse(x$display[[n]]), "\n\n")
    }
  } else {
    cat("None\n\n")
  }

  .sep_line()

  return(invisible(x))
}

