#### Test .fit_null_model() ####
test_that(".fit_null_model correctly fits null pglmm, extracts rank, k and logLik and returns correct items", {

  res <- .fit_null_model(
    formulaRE = log_onset ~ 1 + (1|Species),
    data = senescence,
    rank = "AICc",
    family = "gaussian",
    cov_ranef = list(Species = tree_ultra)
  )

  ## structure
  expect_type(res, "list")
  expect_length(res, 9)

  ### formula.RE
  expect_s3_class(res$formula.RE, "formula")

  ### model.RE
  expect_true(inherits(res$model.RE, "pglmm"))
  expect_equal(res$formula.RE, res$model.RE$formula_original)

  ### k.RE
  expect_equal(res$k.RE, nrow(phyr::fixef(res$model.RE)) + nrow(phyr::ranef(res$model.RE)))

  ### logLik.RE
  expect_type(res$logLik.RE, "double")

  ### rank
  expect_type(res$rank.RE, "double")
  expect_equal(res$rank, "AICc")

  ### data
  expect_s3_class(res$data, "data.frame")

  ### family
  expect_type(res$family, "character")
  expect_equal(res$family, "gaussian")

  ### cov_ranef
  expect_type(res$cov_ranef, "list")
  expect_type(res$cov_ranef$Species, "list")
  expect_true(inherits(res$cov_ranef$Species, "phylo"))
})

#### Test .fit_fixed_model() ####
test_that(".fit_fixed_model correctly fits null pglmm, extracts rank, k and logLik and returns correct items", {

  null_mod <- .fit_null_model(
    formulaRE = log_onset ~ 1 + (1|Species),
    data = senescence,
    rank = "AIC",
    family = "gaussian",
    cov_ranef = list(Species = tree_ultra)
  )

  fixed <- c("log_mass", "log_afr")

  res <- .fit_fixed_model(null_mod = null_mod,
                          fixed = c("log_mass", "log_afr"))

  ## number of possible models (excluding Intercept only)
  nm <- 2^length(fixed) - 1

  ## structure
  expect_type(res, "list")
  expect_length(res, 5)

  ### model_list.fix
  expect_length(res$model_list.fix, nm)
  expect_true(all(vapply(res$model.fix, function(m) inherits(m, "pglmm"), logical(1))))
  expect_equal(res$model_list.fix[[1]]$formula_original, log_onset ~ log_mass + (1|Species))
  expect_equal(res$model_list.fix[[2]]$formula_original, log_onset ~ log_afr + (1|Species))
  expect_equal(res$model_list.fix[[3]]$formula_original, log_onset ~ log_mass + log_afr + (1|Species))

  ### k.fix
  expect_length(res$k.fix, nm)
  expect_type(res$k.fix, "double")
  expect_equal(res$k.fix[1], nrow(phyr::fixef(res$model_list.fix[[1]])) + nrow(phyr::ranef(res$model_list.fix[[1]])))
  expect_equal(res$k.fix[2], nrow(phyr::fixef(res$model_list.fix[[2]])) + nrow(phyr::ranef(res$model_list.fix[[2]])))
  expect_equal(res$k.fix[3], nrow(phyr::fixef(res$model_list.fix[[3]])) + nrow(phyr::ranef(res$model_list.fix[[3]])))

  ### logLik.fix
  expect_length(res$logLik.fix, nm)
  expect_type(res$logLik.fix, "double")

  ### rank.fix
  expect_length(res$rank.fix, nm)
  expect_type(res$rank.fix, "double")

  ### rank_name
  expect_equal(res$rank_name, "AIC")
})

#### Test .dredge_table() ####
test_that(".dredge_table correctly generate a dataframe and characters for print", {

  null_mod <- .fit_null_model(
    formulaRE = log_onset ~ 1 + (1|Species),
    data = senescence,
    rank = "AIC",
    family = "gaussian",
    cov_ranef = list(Species = tree_ultra)
  )

  fixed <- c("log_mass", "log_afr")

  fixed_mod <- .fit_fixed_model(null_mod = null_mod,
                          fixed = c("log_mass", "log_afr"))

  ## number of possible models (excluding Intercept only)
  nm <- 2^length(fixed) - 1

  res <- .dredge_table(null_mod,
                       fixed_mod,
                       estimate = TRUE,
                       std.err = TRUE,
                       round = 3)

  ## structure
  expect_type(res, "list")
  expect_length(res, 4)

  ### dredge_table
  expect_s3_class(res$dredge_table, "data.frame")
  expect_equal(nrow(res$dredge_table), nm + 1)
  expect_equal(ncol(res$dredge_table), length(fixed) + 1 + 5) # 2 fixed effects + Intercept(1) + k, logLik, rank, delta, weight(5)
  expect_equal(colnames(res$dredge_table[6]), "AIC")
  expect_equal(sum(res$dredge_table$weight), 1, tolerance = 1e-8)
  expect_true(all(diff(res$dredge_table$AIC) >= 0))

  ### full_model_chr
  expect_type(res$full_model_chr, "character")
  expect_equal(res$full_model_chr, "log_onset ~ log_mass + log_afr + (1 | Species)")
  expect_length(res$full_model_chr, 1)

  ### random_display
  expect_type(res$random_display, "character")
  expect_equal(res$random_display, "(1 | Species)")
  expect_length(res$full_model_chr, 1)

  ### rank_name
  expect_equal(res$rank_name, "AIC")

  ## Estimate and standard-errors
  fixed_cols_est_se <- res$dredge_table[, c("log_mass", "log_afr", "(Intercept)")]
  non_na_est_se <- unlist(fixed_cols_est_se)[!is.na(unlist(fixed_cols_est_se))]
  expect_true(all(grepl("^-?[0-9.]+ \\([0-9.]+\\)$", non_na_est_se)))

  ## Estimates only
  res_est <- .dredge_table(null_mod,
                           fixed_mod,
                           estimate = TRUE,
                           std.err = FALSE,
                           round = 3)
  fixed_cols_est <- res_est$dredge_table[, c("log_mass", "log_afr", "(Intercept)")]
  non_na_est <- unlist(fixed_cols_est)[!is.na(unlist(fixed_cols_est))]
  expect_true(all(vapply(non_na_est, is.numeric, logical(1))))

  ## No estimates or se
  res_no <- .dredge_table(null_mod,
                          fixed_mod,
                          estimate = FALSE,
                          std.err = FALSE,
                          round = 3)
  fixed_cols <- res_no$dredge_table[, c("log_mass", "log_afr", "(Intercept)")]
  non_na <- unlist(fixed_cols)[!is.na(unlist(fixed_cols))]
  expect_true(all(non_na == "+"))
})

