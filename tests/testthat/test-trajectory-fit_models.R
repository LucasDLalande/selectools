#### Test .fit_non_threshold_models() ####
## correct results
test_that(".fit_non_threshold_models correctly fit models", {

  set.seed(123)

  # standardized data to avoid scaling warnings
  df <- data.frame(x = rnorm(100, mean = 0, sd = 1),
                   y = rnorm(100, mean = 0, sd = 1),
                   z = rnorm(100, mean = 0, sd = 1),
                   group = factor(c(rep("A", 25), rep("B", 25), rep("C", 25), rep("D", 25))))

  data_name <- substitute(df)

  ## LMM
  formulas_lmer <- list(
    linear = y ~ x + z + x:z + (1|group),
    quadratic = y ~ x + z + x:z + I(z^2) + x:I(z^2) + (1|group))

  random_lmer <- findbars(formulas_lmer[[1]])

  ## LM
  formulas_lm <- list(
    linear = y ~ x + z + x:z,
    quadratic = y ~ x + z + x:z + I(z^2) + x:I(z^2))

  random_lm <- findbars(formulas_lm[[1]])

  ## control
  control <- lmerControl(optimizer = "bobyqa",
                         optCtrl = list(maxfun = 1e5))

  # ---- lmer ----
  res_lmer <- .fit_non_threshold_models(non_threshold_formulas = formulas_lmer,
                                        data = df,
                                        data_name = data_name,
                                        random = random_lmer,
                                        lmer_control = control)

  ## structure
  expect_type(res_lmer, "list")
  expect_length(res_lmer, 2)

  ## class
  expect_s4_class(res_lmer[[1]], "lmerMod")
  expect_s4_class(res_lmer[[2]], "lmerMod")

  ## formulas preserved
  expect_equal(formula(res_lmer[[1]]), formulas_lmer[[1]], ignore_attr = TRUE)
  expect_equal(formula(res_lmer[[2]]), formulas_lmer[[2]], ignore_attr = TRUE)

  ## fitted on correct data
  expect_equal(res_lmer[[1]]@call$data, df)

  # ---- lm ----
  res_lm <- .fit_non_threshold_models(non_threshold_formulas = formulas_lm,
                                      data = df,
                                      data_name = data_name,
                                      random = random_lm,
                                      lmer_control = control)

  ## structure
  expect_type(res_lm, "list")
  expect_length(res_lm, 2)

  ## class
  expect_s3_class(res_lm[[1]], "lm")
  expect_s3_class(res_lm[[2]], "lm")

  ## formulas preserved
  expect_equal(formula(res_lm[[1]]), formulas_lm[[1]], ignore_attr = TRUE)
  expect_equal(formula(res_lm[[2]]), formulas_lm[[2]], ignore_attr = TRUE)

  ## fitted on correct data
  expect_equal(res_lm[[1]]$call$data, df)
})

#### Test .make_threshold_dataset() ####
## correct results
test_that(".make_threshold_dataset correctly generates a list of datasets with var1 and var2 calculated with best thresholds", {

  # small sample size to facilitate calculation of expected values
  df <- data.frame(x = 11:15, y = 1:5, z = 1:5, group = factor(c(rep("A", 3), rep("B",2))))

  best_thresholds <- list(
    thresholdCS = 3,
    thresholdSC = 4,
    thresholdSS = 2
  )

  res <- .make_threshold_dataset(trajectory_trait = "z",
                                 data = df,
                                 best_thresholds = best_thresholds)

  # structure
  expect_type(res, "list")
  expect_length(res, 3)

  # class
  expect_s3_class(res[[1]], "data.frame")
  expect_s3_class(res[[2]], "data.frame")
  expect_s3_class(res[[3]], "data.frame")

  # var1 and var2 calculations
  expect_equal(res[[1]]$z.2, c(0,0,0,1,2))
  expect_equal(res[[2]]$z.1, c(-3,-2,-1,0,0))
  expect_equal(res[[3]]$z.1, c(1,2,2,2,2))
  expect_equal(res[[3]]$z.2, c(0,0,1,2,3))
})


#### Test .make_model_env() ####
## correct results
test_that(".make_model_env correctly generates an environment with attached datasets", {

  set.seed(123)

  # standardized data to avoid scaling warnings
  df <- data.frame(x = rnorm(100, mean = 0, sd = 1),
                   y = rnorm(100, mean = 0, sd = 1),
                   z = rnorm(100, mean = 0, sd = 1),
                   group = factor(c(rep("A", 25), rep("B", 25), rep("C", 25), rep("D", 25))))

  best_thresholds <- list(
    thresholdCS = 0.3,
    thresholdSC = 0.4,
    thresholdSS = 0.2
  )

  datasets_threshold <- .make_threshold_dataset(trajectory_trait = "z",
                                                data = df,
                                                best_thresholds = best_thresholds)

  res <- .make_model_env(datasets_threshold)

  ## model_env
  expect_type(res, "environment")
  expect_length(res, length(best_thresholds))

  # ---- model environment ----
  ## threshold datasets exist in the models environment
  expect_true(exists(".data_CS", envir = res, inherits = FALSE))
  expect_true(exists(".data_SC", envir = res, inherits = FALSE))
  expect_true(exists(".data_SS", envir = res, inherits = FALSE))

  ## threshold datasets are data.frame objects
  expect_true(is.data.frame(get(".data_CS", envir = res)))
  expect_true(is.data.frame(get(".data_SC", envir = res)))
  expect_true(is.data.frame(get(".data_SS", envir = res)))
})


#### Test .fit_threshold_models() ####
## correct results
test_that(".fit_threshold_models correctly fit models", {

  set.seed(123)

  # standardized data to avoid scaling warnings
  df <- data.frame(x = rnorm(100, mean = 0, sd = 1),
                   y = rnorm(100, mean = 0, sd = 1),
                   z = rnorm(100, mean = 0, sd = 1),
                   group = factor(c(rep("A", 25), rep("B", 25), rep("C", 25), rep("D", 25))))

  df_thr <- list(thresholdCS = data.frame(x = rnorm(100, mean = 0, sd = 1),
                                          y = rnorm(100, mean = 0, sd = 1),
                                          z.2 = rnorm(100, mean = 0, sd = 1),
                                          group = factor(c(rep("A", 25), rep("B", 25), rep("C", 25), rep("D", 25)))),
                 thresholdSC = data.frame(x = rnorm(100, mean = 0, sd = 1),
                                          y = rnorm(100, mean = 0, sd = 1),
                                          z.1 = rnorm(100, mean = 0, sd = 1),
                                          group = factor(c(rep("A", 25), rep("B", 25), rep("C", 25), rep("D", 25)))),
                 thresholdSS = data.frame(x = rnorm(100, mean = 0, sd = 1),
                                          y = rnorm(100, mean = 0, sd = 1),
                                          z.1 = rnorm(100, mean = 0, sd = 1),
                                          z.2 = rnorm(100, mean = 0, sd = 1),
                                          group = factor(c(rep("A", 25), rep("B", 25), rep("C", 25), rep("D", 25))))
  )

  ## LMM
  formulas_lmer <- list(
    thresholdCS = y ~ x + z.2 + x:z.2 + (1|group),
    thresholdSC = y ~ x + z.1 + x:z.1 + (1|group),
    thresholdSS = y ~ x + z.1 + x:z.1 + z.2 + x:z.2 + (1|group))

  random_lmer <- findbars(formulas_lmer[[1]])

  ## LM
  formulas_lm <- list(
    thresholdCS = y ~ x + z.2 + x:z.2,
    thresholdSC = y ~ x + z.1 + x:z.1,
    thresholdSS = y ~ x + z.1 + x:z.1 + z.2 + x:z.2)

  random_lm <- findbars(formulas_lm[[1]])

  ## control
  control <- lmerControl(optimizer = "bobyqa",
                         optCtrl = list(maxfun = 1e5))

  # ---- lmer ----
  res_lmer <- .fit_threshold_models(datasets_threshold = df_thr,
                                    threshold_formulas = formulas_lmer,
                                    random = random_lmer,
                                    lmer_control = control)

  ## structure
  expect_type(res_lmer, "list")
  expect_length(res_lmer, 2)
  expect_type(res_lmer[[1]], "list")
  expect_type(res_lmer[[2]], "environment")

  ## class
  expect_s4_class(res_lmer$models[[1]], "lmerMod")
  expect_s4_class(res_lmer$models[[2]], "lmerMod")
  expect_s4_class(res_lmer$models[[3]], "lmerMod")

  ## formulas preserved
  expect_equal(formula(res_lmer$models[[1]]), formulas_lmer[[1]], ignore_attr = TRUE)
  expect_equal(formula(res_lmer$models[[2]]), formulas_lmer[[2]], ignore_attr = TRUE)
  expect_equal(formula(res_lmer$models[[3]]), formulas_lmer[[3]], ignore_attr = TRUE)

  ## fitted on correct data
  expect_equal(res_lmer$models[[1]]@call$data, quote(.data_CS))
  expect_equal(res_lmer$models[[2]]@call$data, quote(.data_SC))
  expect_equal(res_lmer$models[[3]]@call$data, quote(.data_SS))

  # ---- lm ----
  res_lm <- .fit_threshold_models(datasets_threshold = df_thr,
                                  threshold_formulas = formulas_lm,
                                  random = random_lm,
                                  lmer_control = control)

  ## structure
  expect_type(res_lm, "list")
  expect_length(res_lm, 2)
  expect_type(res_lm[[1]], "list")
  expect_type(res_lm[[2]], "environment")

  ## class
  expect_s3_class(res_lm$models[[1]], "lm")
  expect_s3_class(res_lm$models[[2]], "lm")
  expect_s3_class(res_lm$models[[3]], "lm")

  ## formulas preserved
  expect_equal(formula(res_lm$models[[1]]), formulas_lm[[1]], ignore_attr = TRUE)
  expect_equal(formula(res_lm$models[[2]]), formulas_lm[[2]], ignore_attr = TRUE)
  expect_equal(formula(res_lm$models[[3]]), formulas_lm[[3]], ignore_attr = TRUE)

  ## fitted on correct data
  expect_equal(res_lm$models[[1]]$call$data, quote(.data_CS))
  expect_equal(res_lm$models[[2]]$call$data, quote(.data_SC))
  expect_equal(res_lm$models[[3]]$call$data, quote(.data_SS))
})


#### Test fit_models() ####
## correct results
test_that("fit_models correctly fits models", {

  set.seed(123)

  # standardized data to avoid scaling warnings
  df <- data.frame(x = rnorm(100, mean = 0, sd = 1),
                   y = rnorm(100, mean = 0, sd = 1),
                   z = rnorm(100, mean = 0, sd = 1),
                   group = factor(c(rep("A", 25), rep("B", 25), rep("C", 25), rep("D", 25))))

  ## LMM
  ### formulas generated by trajectory_formulas()
  f_lmer <- trajectory_formulas(formula = y ~ x + (1|group),
                                trajectory_trait = "z",
                                interactions = "x",
                                data = df,
                                trajectories = "all")

  ### threshold estimation generated by estimate_threshold()
  estimation_lmer <- estimate_threshold(data = df,
                                        threshold_formulas = f_lmer,
                                        trajectory_trait = "z",
                                        start = -0.8, end = 0.8, iter = 0.2,
                                        rank = "AICc")
  ## LM
  ### formulas provided manually by users
  formulas_lm <- list(
    linear = y ~ x + z + x:z,
    thresholdCS = y ~ x + z.2 + x:z.2
    )

  ### threshold values provided manually by users
  estimation_threshold_lm <- list(
    thresholdCS = 0.4
  )

  # ---- lmer ----
  ## Objects returned from trajectory_formulas() and estimate_threshold()
  ### convergence, identification or Hessian warnings can occur with lme4 on simulated datasets
  ### we intently ignore those as it does not relate to function issue
  m_lmer <- withCallingHandlers(
    fit_models(formulas = f_lmer,
               trajectory_trait = "z",
               data = df,
               best_thresholds = estimation_lmer),
    warning = function(w) {
      # ignore expected convergence warnings from lme4 on simulated datasets
      if (grepl("failed to converge|nearly unidentifiable|Hessian",
                conditionMessage(w)
      )) {
        invokeRestart("muffleWarning")
      }
    }
  )

  ### general structure
  expect_type(m_lmer, "list")
  expect_length(m_lmer, 7)
  expect_type(c(m_lmer$trajectory_trait, m_lmer$random), "character")
  expect_type(m_lmer$data, "symbol")
  expect_type(m_lmer$model_env, "environment")

  ### datasets_threshold
  expect_type(m_lmer$datasets_threshold, "list")
  expect_length(m_lmer$datasets_threshold, 3)
  expect_s3_class(m_lmer$datasets_threshold[[1]], "data.frame")
  expect_s3_class(m_lmer$datasets_threshold[[2]], "data.frame")
  expect_s3_class(m_lmer$datasets_threshold[[3]], "data.frame")

  ### models
  expect_type(m_lmer$models, "list")
  expect_length(m_lmer$models, 7)
  expect_s4_class(m_lmer$models[[1]], "lmerMod")
  expect_s4_class(m_lmer$models[[2]], "lmerMod")
  expect_s4_class(m_lmer$models[[3]], "lmerMod")
  expect_s4_class(m_lmer$models[[4]], "lmerMod")
  expect_s4_class(m_lmer$models[[5]], "lmerMod")
  expect_s4_class(m_lmer$models[[6]], "lmerMod")
  expect_s4_class(m_lmer$models[[7]], "lmerMod")

  ### best_thresholds
  expect_type(m_lmer$best_thresholds, "list")
  expect_length(m_lmer$best_thresholds, 3)
  expect_equal(class(m_lmer$best_thresholds[[1]]), "numeric")
  expect_equal(class(m_lmer$best_thresholds[[2]]), "numeric")
  expect_equal(class(m_lmer$best_thresholds[[3]]), "numeric")

  ### class
  expect_s3_class(m_lmer, "models")

  ### formulas preserved
  formulas_lmer <- list(
    linear = y ~ x + z + x:z + (1|group),
    thresholdCS = y ~ x + z.2 + x:z.2 + (1|group)
  )

  expect_setequal(
    attr(terms(formula(m_lmer$models[["linear"]])), "term.labels"),
    attr(terms(formulas_lmer[["linear"]]), "term.labels"))

  expect_setequal(
    attr(terms(formula(m_lmer$models[["thresholdCS"]])), "term.labels"),
    attr(terms(formulas_lmer[["thresholdCS"]]), "term.labels"))

  ### fitted on correct data
  expect_equal(m_lmer$models[["thresholdCS"]]@call$data, quote(.data_CS))
  expect_equal(m_lmer$models[["thresholdSC"]]@call$data, quote(.data_SC))
  expect_equal(m_lmer$models[["thresholdSS"]]@call$data, quote(.data_SS))

  # ---- lm ----
  ## objects provided manually by users
  ### we ignore the single expected warning concerning threshold_formulas and best_thresholds not generated by trajectory_formulas() and estimate_threshold() respectively
  models_lm <- withCallingHandlers(
    fit_models(formulas = formulas_lm,
               trajectory_trait = "z",
               data = df,
               best_thresholds = estimation_threshold_lm),
    warning = function(w) {
      # ignore expected warning on arguments not being generated by the package functions
      if (grepl("were not generated by",
                conditionMessage(w)
      )) {
        invokeRestart("muffleWarning")
      }
    }
  )

  ### general structure
  expect_type(models_lm, "list")
  expect_length(models_lm, 7)
  expect_type(c(models_lm$trajectory_trait, models_lm$random), "character")
  expect_type(models_lm$data, "symbol")
  expect_type(models_lm$model_env, "environment")

  ### datasets_threshold
  expect_type(models_lm$datasets_threshold, "list")
  expect_length(models_lm$datasets_threshold, 1)
  expect_s3_class(models_lm$datasets_threshold[[1]], "data.frame")

  ### models
  expect_type(models_lm$models, "list")
  expect_length(models_lm$models, 2)
  expect_s3_class(models_lm$models[[1]], "lm")
  expect_s3_class(models_lm$models[[2]], "lm")

  ### best_thresholds
  expect_type(models_lm$best_thresholds, "list")
  expect_length(models_lm$best_thresholds, 1)
  expect_equal(class(models_lm$best_thresholds[[1]]), "numeric")

  ### class
  expect_s3_class(models_lm, "models")

  ### formulas preserved
  expect_setequal(
    attr(terms(formula(models_lm$models[["linear"]])), "term.labels"),
    attr(terms(formulas_lm[["linear"]]), "term.labels"))

  expect_setequal(
    attr(terms(formula(models_lm$models[["thresholdCS"]])), "term.labels"),
    attr(terms(formulas_lm[["thresholdCS"]]), "term.labels"))

  ### fitted on correct data
  expect_equal(models_lm$models[["thresholdCS"]]$call$data, quote(.data_CS))
})

## errors when invalid arguments
test_that("fit_models errors when invalid arguments", {

  set.seed(123)

  # standardized data to avoid scaling warnings
  df <- data.frame(x = rnorm(100, mean = 0, sd = 1),
                   y = rnorm(100, mean = 0, sd = 1),
                   z = rnorm(100, mean = 0, sd = 1),
                   group = factor(c(rep("A", 25), rep("B", 25), rep("C", 25), rep("D", 25))))

  ## LMM
  ### formulas generated by trajectory_formulas()
  f_lmer <- trajectory_formulas(formula = y ~ x + (1|group),
                                trajectory_trait = "z",
                                interactions = "x",
                                data = df,
                                trajectories = "all")

  ### formulas provided manually by users
  formulas_lmer <- list(
    constant = y ~ x + (1|group),
    linear = y ~ x + z + x:z + (1|group),
    quadratic = y ~ x + z + I(z^2) + x:z + x:I(z^2) + (1|group),
    quadratic = y ~ x + factor(z) + x:factor(z) + (1|group),
    thresholdCS = y ~ x + z.2 + x:z.2 + (1|group),
    thresholdSC = y ~ x + z.1 + x:z.1 + (1|group),
    thresholdSS = y ~ x + z.1 + z.2 + x:z.1 + x:z.2 + (1|group)
  )

  ### threshold estimation generated by estimate_threshold()
  estimation_lmer <- estimate_threshold(data = df, threshold_formulas = f_lmer,
                                        trajectory_trait = "z",
                                        start = -0.8, end = 0.8, iter = 0.2,
                                        rank = "AICc")

  ### threshold values provided manually by users
  estimation_threshold_lmer <- list(
    thresholdCS = 0.4,
    thresholdSC = 0.5,
    thresholdSS = 0.3
  )

  ## LM
  ### formulas provided manually by users
  formulas_lm <- list(
    linear = y ~ x + z + x:z,
    thresholdCS = y ~ x + z.2 + x:z.2
  )

  ### threshold values provided manually by users
  estimation_threshold_lm <- list(
    thresholdCS = 14
  )

  # ---- data ----
  ## not a data.frame
  wrongdata <- "not a data.frame"
  expect_error(
    fit_models(formulas = f_lmer,
               trajectory_trait = "z",
               data = wrongdata,
               best_thresholds = estimation_lmer),
  "data must be provided and must be an existing data.frame")

  # ---- formulas ----
  ## no valid formulas
  expect_error(
    fit_models(formulas = list(),
               trajectory_trait = "z",
               data = df,
               best_thresholds = estimation_lmer),
  "No valid formulas provided")

  ## unknown formula name
  expect_error(
    fit_models(formulas = list(constant = y ~ x + (1|group),
                               formula1 = y ~ x + z + x:z + (1|group),
                               thresholdCS = y ~ x + z.2 + x:z.2 + (1|group)),
               trajectory_trait = "z",
               data = df,
               best_thresholds = estimation_lmer),
  "Unknown trajectory types in formulas. Names must be from: constant, linear, quadratic, factor, thresholdCS, thresholdSC, thresholdSS")

  ## not a formula object
  expect_error(
    fit_models(formulas = list(linear = 123),
               trajectory_trait = "z",
               data = df,
               best_thresholds = estimation_lmer),
  "non_threshold_formulas must be a list of formula objects")

  ## not a two-sided formula
  expect_error(
    fit_models(formulas = list(linear = ~ x + z + x:z + (1|group)),
               trajectory_trait = "z",
               data = df,
               best_thresholds = estimation_lmer),
  "non_threshold_formulas must be a list of two-sided formulas")

  ## not all elements are formula objects
  expect_error(
    fit_models(formulas = list(linear = y ~ x + z + x:z + (1|group),
                               thresholdCS = 123),
               trajectory_trait = "z",
               data = df,
               best_thresholds = estimation_lmer),
  "threshold_formulas must be a list of formula objects")

  ## not all elements are two-sided formulas
  expect_error(
    fit_models(formulas = list(linear = y ~ x + z + x:z + (1|group),
                               thresholdCS = ~ x + z.2 + x:z.2 + (1|group)),
               trajectory_trait = "z",
               data = df,
               best_thresholds = estimation_lmer),
  "threshold_formulas must be a list of two-sided formulas")

  ## warning when formulas provided by users
  expect_snapshot_warning(
    withCallingHandlers(
      fit_models(formulas = formulas_lmer,
                 trajectory_trait = "z",
                 data = df,
                 best_thresholds = estimation_lmer),
       warning = function(w) {
         # ignore expected convergence warnings from lme4 on simulated datasets
         if (grepl("failed to converge|nearly unidentifiable|Hessian",
                   conditionMessage(w)
            )) {
           invokeRestart("muffleWarning")
         }
       }
    )
  )

  # ---- trajectory_trait ----
  ## not a single character string
  expect_error(
    fit_models(formulas = f_lmer,
               trajectory_trait = c(123, "wrong"),
               data = df,
               best_thresholds = estimation_lmer),
  "trajectory_trait must be a single character string")

  ## not a variable from data
  expect_error(
    fit_models(formulas = f_lmer,
               trajectory_trait = "wrong",
               data = df,
               best_thresholds = estimation_lmer),
  "trajectory_trait must be a variable of 'data'")

  # ---- best_thresholds ----
  ## no threshold estimations while threshold formulas are provided
  expect_error(
    fit_models(formulas = f_lmer,
               trajectory_trait = "z",
               data = df,
               best_thresholds = NULL),
  "best_thresholds must be provided when threshold models are included")

  ## wrong threshold estimation name
  expect_error(
    fit_models(formulas = f_lmer,
               trajectory_trait = "z",
               data = df,
               best_thresholds = list(
                 linear = 13,
                 thresholdCS = 14)
               ),
  "Unknown threshold types in best_thresholds. Names must be from: thresholdCS, thresholdSC, thresholdSS")

  ## not a numeric value
  expect_error(
    fit_models(formulas = f_lmer,
               trajectory_trait = "z",
               data = df,
               best_thresholds = list(
                 thresholdCS = "not a numeric")
    ),
  "best_thresholds must be a non-empty list of numeric objects")

  ## not a list
  expect_error(
    fit_models(formulas = f_lmer,
               trajectory_trait = "z",
               data = df,
               best_thresholds = c(
                 thresholdCS = 14)
    ),
  "best_thresholds must be a non-empty list of numeric objects")

  ## warning when threshold estimation is provided by users
  expect_snapshot_warning(
    withCallingHandlers(
      fit_models(formulas = f_lmer,
                 trajectory_trait = "z",
                 data = df,
                 best_thresholds = estimation_threshold_lmer),
      warning = function(w) {
        # ignore expected convergence warnings from lme4 on simulated datasets
        if (grepl("failed to converge|nearly unidentifiable|Hessian",
                  conditionMessage(w)
        )) {
          invokeRestart("muffleWarning")
        }
      }
    )
  )

  ## warning when threshold estimations are provided while no threshold formulas
  ### formulas generated by trajectory_formulas()
  f_lmer_nonthr <- trajectory_formulas(formula = y ~ x + (1|group),
                                       trajectory_trait = "z",
                                       interactions = "x",
                                       data = df,
                                       trajectories = "non-threshold")
  expect_snapshot_warning(
    withCallingHandlers(
      fit_models(formulas = f_lmer_nonthr,
                 trajectory_trait = "z",
                 data = df,
                 best_thresholds = estimation_lmer),
      warning = function(w) {
        # ignore expected convergence warnings from lme4 on simulated datasets
        if (grepl("failed to converge|nearly unidentifiable|Hessian",
                  conditionMessage(w)
        )) {
          invokeRestart("muffleWarning")
        }
      }
    )
  )

  ## error when some threshold estimations are missing
  ### missing threshold values
  estimation_threshold_lmer_missing <- list(
    thresholdCS = 0.4,
    thresholdSC = 0.5
  )
  expect_error(
    withCallingHandlers(
      fit_models(formulas = f_lmer,
                 trajectory_trait = "z",
                 data = df,
                 best_thresholds = estimation_threshold_lmer_missing),
      warning = function(w) {
        # ignore expected convergence warnings from lme4 on simulated datasets
        # and ignore warning concerning best_thresholds not generated by estimate_threshold()
        if (grepl("failed to converge|nearly unidentifiable|Hessian|were not generated by",
                  conditionMessage(w)
        )) {
          invokeRestart("muffleWarning")
        }
      }
    ),
    "Missing threshold values")

  # ---- lmer_control ----
  lmer_control_wrong <- "lmer_control"
  expect_error(
    fit_models(formulas = f_lmer,
               trajectory_trait = "z",
               data = df,
               best_thresholds = estimation_lmer,
               lmer_control = lmer_control_wrong),
  "lmer_control must be NULL or a lmerControl, merControl object. See help page for details")

  # ---- random == NULL, lmer_control == NULL ----
  expect_no_error(
    withCallingHandlers(
      fit_models(formulas = f_lmer,
                 trajectory_trait = "z",
                 data = df,
                 best_thresholds = estimation_lmer,
                 lmer_control = NULL),
      warning = function(w) {
        # ignore expected convergence warnings from lme4 on simulated datasets
        if (grepl("failed to converge|nearly unidentifiable|Hessian",
                  conditionMessage(w)
        )) {
          invokeRestart("muffleWarning")
        }
      }
    )
  )
})

#### Test print.models() ####
test_that("print.models works", {

  set.seed(123)

  # standardized data to avoid scaling warnings
  df <- data.frame(x = rnorm(100, mean = 0, sd = 1),
                   y = rnorm(100, mean = 0, sd = 1),
                   z = rnorm(100, mean = 0, sd = 1),
                   group = factor(c(rep("A", 25), rep("B", 25), rep("C", 25), rep("D", 25))))

  # formulas generated by trajectory_formulas()
  f_lmer <- trajectory_formulas(formula = y ~ x + (1|group),
                                trajectory_trait = "z",
                                interactions = "x",
                                data = df,
                                trajectories = "all")

  # threshold estimation with estimate_threshold()
  estimation_lmer <- estimate_threshold(data = df,
                                        threshold_formulas = f_lmer,
                                        trajectory_trait = "z",
                                        start = -0.8, end = 0.8, iter = 0.2,
                                        rank = "AICc")

  # fit_models()
    res <- withCallingHandlers(
      fit_models(formulas = f_lmer,
                 trajectory_trait = "z",
                 data = df,
                 best_thresholds = estimation_lmer),
      warning = function(w) {
        # ignore expected convergence warnings from lme4 on simulated datasets
        if (grepl("failed to converge|nearly unidentifiable|Hessian",
                  conditionMessage(w)
        )) {
          invokeRestart("muffleWarning")
        }
      }
    )

  expect_invisible(print(res))
})
