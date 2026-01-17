# Test dredge_pglmm()
test_that("dredge_pglmm runs correctly", {
  expect_invisible(
    dredge_pglmm(
        formulaRE = "log_onset ~ 1 + (1 | Species)",
        fixed = c("log_afr", "log_mass"),
        data = senescence,
        rank = "AICc",
        family = "gaussian",
        cov_ranef = list(Species = tree_ultra),
        estimate = T,
        std.err = T,
        round = 5))
})

test_that("dredge_pglmm errors to invalid arguments", {

  # error to invalid random effect formula
  expect_error(
    dredge_pglmm(
      formulaRE = 123,
      fixed = c("log_afr", "log_mass"),
      data = senescence,
      rank = "AICc",
      family = "gaussian",
      cov_ranef = list(Species = tree_ultra),
      estimate = T,
      std.err = T,
      round = 5))

  # error to null fixed
  expect_error(
    dredge_pglmm(
      formulaRE = "log_onset ~ 1 + (1|Species)",
      fixed = NULL,
      data = senescence,
      rank = "AICc",
      family = "gaussian",
      cov_ranef = list(Species = tree_ultra),
      estimate = T,
      std.err = T,
      round = 5))

  # error to non existing dataset
  expect_error(
    dredge_pglmm(
      formulaRE = "log_onset ~ 1 + (1|Species)",
      fixed = c("log_mass", "log_afr"),
      data = "wrongdata",
      rank = "AICc",
      family = "gaussian",
      cov_ranef = list(Species = tree_ultra),
      estimate = T,
      std.err = T,
      round = 5))

  # error to object that is not a dataframe
  datatest <- c(1,2,3,4,5)

  expect_error(
    dredge_pglmm(
      formulaRE = "log_onset ~ 1 + (1|Species)",
      fixed = c("log_mass", "log_afr"),
      data = datatest,
      rank = "AICc",
      family = "gaussian",
      cov_ranef = list(Species = tree_ultra),
      estimate = T,
      std.err = T,
      round = 5))

  # error to invalid rank
  expect_error(
    dredge_pglmm(
      formulaRE = "log_onset ~ 1 + (1|Species)",
      fixed = c("log_mass", "log_afr"),
      data = senescence,
      rank = "BIC",
      family = "gaussian",
      cov_ranef = list(Species = tree_ultra),
      estimate = T,
      std.err = T,
      round = 5))

  # error to invalid family
  expect_error(
    dredge_pglmm(
      formulaRE = "log_onset ~ 1 + (1|Species)",
      fixed = c("log_mass", "log_afr"),
      data = senescence,
      rank = "AICc",
      family = "neg.binom",
      cov_ranef = list(Species = tree_ultra),
      estimate = T,
      std.err = T,
      round = 5))

  # error to invalid covariance matrix - not a list
  expect_error(
    dredge_pglmm(
      formulaRE = "log_onset ~ 1 + (1|Species)",
      fixed = c("log_mass", "log_afr"),
      data = senescence,
      rank = "AICc",
      family = "gaussian",
      cov_ranef = datatest,
      estimate = T,
      std.err = T,
      round = 5))

  # error to invalid covariance matrix - a list in which elements are not matrix or phylo object
  expect_error(
    dredge_pglmm(
      formulaRE = "log_onset ~ 1 + (1|Species)",
      fixed = c("log_mass", "log_afr"),
      data = senescence,
      rank = "AICc",
      family = "gaussian",
      cov_ranef = list(Species = datatest),
      estimate = T,
      std.err = T,
      round = 5))

  # error to decimal for round
  expect_error(
    dredge_pglmm(
      formulaRE = "log_onset ~ 1 + (1|Species)",
      fixed = c("log_mass", "log_afr"),
      data = senescence,
      rank = "AICc",
      family = "gaussian",
      cov_ranef = list(Species = tree_ultra),
      estimate = T,
      std.err = T,
      round = 5.4))

  # error to negative integer for round
  expect_error(
    dredge_pglmm(
      formulaRE = "log_onset ~ 1 + (1|Species)",
      fixed = c("log_mass", "log_afr"),
      data = senescence,
      rank = "AICc",
      family = "gaussian",
      cov_ranef = list(Species = tree_ultra),
      estimate = T,
      std.err = T,
      round = -5))
})

