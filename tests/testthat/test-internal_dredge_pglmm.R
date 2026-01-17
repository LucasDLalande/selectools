# Test .fit_null_model() ----
test_that(".fit_null_model errors when invalid random effect formula and family", {
  expect_error(.fit_null_model(13, data = senescence, family = "gaussian", cov_ranef = list(Species = tree_ultra)))
  expect_error(.fit_null_model("log_onset ~ 1 + (1|Species)", data = senescence, family = "neg.binom", cov_ranef = list(Species = tree_ultra)))
})

test_that(".fit_null_model correctly generates and extracts df, log-likelihood, AICc and AIC", {
  res <- .fit_null_model(
    formulaRE = "log_onset ~ 1 + (1|Species)",
    data = senescence,
    family = "gaussian",
    cov_ranef = list(Species = tree_ultra)
  )

  expect_type(res, "list")
  expect_named(res, c("formula.RE", "model.RE", "df.RE", "logLik.RE", "AICc.RE", "AIC.RE"))
  expect_s3_class(res$formula.RE, "formula")
  expect_true(inherits(res$model.RE, "pglmm"))

  expect_equal(res$df.RE, nrow(lme4::fixef(res$model.RE)) + nrow(lme4::ranef(res$model.RE)))
  expect_gt(res$AICc.RE, res$AIC.RE)
})

# Test .fit_fixed_model() ----
test_that(".fit_fixed_model errors when invalid fixed effects and family", {

  null_mod <- .fit_null_model(
    formulaRE = "log_onset ~ 1 + (1|Species)",
    data = senescence,
    family = "gaussian",
    cov_ranef = list(Species = tree_ultra)
  )

  expect_error(.fit_fixed_model(null_mod = null_mod, fixed, data = senescence, family = "gaussian", cov_ranef = list(Species = tree_ultra)))
  expect_error(.fit_null_model(null_mod = null_mod, fixed = c("log_mass", "log_afr"), data = senescence, family = "neg.binom", cov_ranef = list(Species = tree_ultra)))
})

test_that(".fit_fixed_model correctly generates and extracts df, log-likelihood, AICc and AIC", {

  null_mod <- .fit_null_model(
    formulaRE = "log_onset ~ 1 + (1|Species)",
    data = senescence,
    family = "gaussian",
    cov_ranef = list(Species = tree_ultra)
  )

  res <- .fit_fixed_model(null_mod = null_mod,
                  fixed = c("log_mass", "log_afr"),
                  data = senescence,
                  family = "gaussian",
                  cov_ranef = list(Species = tree_ultra))

  expect_type(res, "list")
  expect_named(res, c("model_list.fix", "df.fix", "logLik.fix", "AICc.fix", "AIC.fix"))
  expect_true(inherits(res$model_list.fix[[2]], "pglmm"))

  expect_equal(res$df.fix[[2]], nrow(lme4::fixef(res$model_list.fix[[2]])) + nrow(lme4::ranef(res$model_list.fix[[2]])))
  expect_gt(res$AICc.fix[[2]], res$AIC.fix[[2]])
})

test_that(".fit_fixed_model generates the good number of model to be tested", {

  null_mod <- .fit_null_model(
    formulaRE = "log_onset ~ 1 + (1|Species)",
    data = senescence,
    family = "gaussian",
    cov_ranef = list(Species = tree_ultra)
  )

  res <- .fit_fixed_model(null_mod = null_mod,
                          fixed = c("log_mass", "log_afr"),
                          data = senescence,
                          family = "gaussian",
                          cov_ranef = list(Species = tree_ultra))

  nb_models <- 2^length(c("log_mass", "log_afr"))-1

  expect_equal(length(res$model_list.fix), nb_models)
})


# Test .dredge_table() ----
test_that(".dredge_table errors when invalid rank", {

  null_mod <- .fit_null_model(
    formulaRE = "log_onset ~ 1 + (1|Species)",
    data = senescence,
    family = "gaussian",
    cov_ranef = list(Species = tree_ultra)
  )

  fixed_mod <- .fit_fixed_model(
    null_mod = null_mod,
    fixed = c("log_mass", "log_afr"),
    data = senescence,
    family = "gaussian",
    cov_ranef = list(Species = tree_ultra))

  expect_error(.dredge_table(null_mod, fixed_mod, rank="BIC", estimate=TRUE, std.err=TRUE, round=3))
})

test_that(".dredge_table correctly generate a dataframe", {

  null_mod <- .fit_null_model(
    formulaRE = "log_onset ~ 1 + (1|Species)",
    data = senescence,
    family = "gaussian",
    cov_ranef = list(Species = tree_ultra)
  )

  fixed_mod <- .fit_fixed_model(
    null_mod = null_mod,
    fixed = c("log_mass", "log_afr"),
    data = senescence,
    family = "gaussian",
    cov_ranef = list(Species = tree_ultra))

  table <- .dredge_table(null_mod, fixed_mod, rank="AICc", estimate=TRUE, std.err=TRUE, round=3)

  expect_true(is.data.frame(table))
})

test_that(".dredge_table generate a dataframe ordered by AICc or AIC and that weight sum =1", {

  null_mod <- .fit_null_model(
    formulaRE = "log_onset ~ 1 + (1|Species)",
    data = senescence,
    family = "gaussian",
    cov_ranef = list(Species = tree_ultra)
  )

  fixed_mod <- .fit_fixed_model(
    null_mod = null_mod,
    fixed = c("log_mass", "log_afr"),
    data = senescence,
    family = "gaussian",
    cov_ranef = list(Species = tree_ultra))

  table <- .dredge_table(null_mod, fixed_mod, rank="AICc", estimate=TRUE, std.err=TRUE, round=3)

  expect_true(all(diff(table$AICc) >= 0))

  expect_equal(sum(table$weight), 1, tolerance = 1e-8)
})

test_that(".dredge_table manage correctly estimates and standard errors", {

  null_mod <- .fit_null_model(
    formulaRE = "log_onset ~ 1 + (1|Species)",
    data = senescence,
    family = "gaussian",
    cov_ranef = list(Species = tree_ultra)
  )

  fixed_mod <- .fit_fixed_model(
    null_mod = null_mod,
    fixed = c("log_mass", "log_afr"),
    data = senescence,
    family = "gaussian",
    cov_ranef = list(Species = tree_ultra))

  table <- .dredge_table(null_mod, fixed_mod, rank="AICc", estimate=FALSE, std.err=FALSE)

  vals <- unlist(table[,c("log_mass", "log_afr")])
  vals <- as.character(vals)
  expect_true(all(na.omit(vals) == "+"))
})

# Test .print_table() ----
test_that(".print_table errors when invalid rank", {

  null_mod <- .fit_null_model(
    formulaRE = "log_onset ~ 1 + (1|Species)",
    data = senescence,
    family = "gaussian",
    cov_ranef = list(Species = tree_ultra)
  )

  fixed_mod <- .fit_fixed_model(
    null_mod = null_mod,
    fixed = c("log_mass", "log_afr"),
    data = senescence,
    family = "gaussian",
    cov_ranef = list(Species = tree_ultra))

  table <- .dredge_table(
    null_mod,
    fixed_mod,
    rank="AICc",
    estimate=TRUE,
    std.err=TRUE,
    round=3)

  expect_error(.print_table(null_mod, fixed_mod, data = senescence, table, rank = "BIC", cov_ranef = list(Species = tree_ultra), data_expr = substitute(senescence)))
})

test_that(".print_table displays expected output in the console", {

  null_mod <- .fit_null_model(
    formulaRE = "log_onset ~ 1 + (1|Species)",
    data = senescence,
    family = "gaussian",
    cov_ranef = list(Species = tree_ultra)
  )

  fixed_mod <- .fit_fixed_model(
    null_mod = null_mod,
    fixed = c("log_mass", "log_afr"),
    data = senescence,
    family = "gaussian",
    cov_ranef = list(Species = tree_ultra))

  table <- .dredge_table(
    null_mod,
    fixed_mod,
    rank="AIC",
    estimate=TRUE,
    std.err=TRUE,
    round=3)

  expect_output(.print_table(
    null_mod,
    fixed_mod,
    data = senescence,
    table,
    rank = "AIC",
    cov_ranef = list(Species = tree_ultra),
    data_expr = substitute(senescence)),
    regexp = "Global model: pglmm")

  expect_output(.print_table(
    null_mod,
    fixed_mod,
    data = senescence,
    table,
    rank = "AIC",
    cov_ranef = list(Species = tree_ultra),
    data_expr = substitute(senescence)),
    regexp = "Data:")

  expect_output(.print_table(
    null_mod,
    fixed_mod,
    data = senescence,
    table,
    rank = "AIC",
    cov_ranef = list(Species = tree_ultra),
    data_expr = substitute(senescence)),
    regexp = "Model selection table")

  expect_output(.print_table(
    null_mod,
    fixed_mod,
    data = senescence,
    table,
    rank = "AIC",
    cov_ranef = list(Species = tree_ultra),
    data_expr = substitute(senescence)),
    regexp = "Models ranked by AIC")

  expect_output(.print_table(
    null_mod,
    fixed_mod,
    data = senescence,
    table,
    rank = "AIC",
    cov_ranef = list(Species = tree_ultra),
    data_expr = substitute(senescence)),
    regexp = "Covariance matrix of random effects")
})
