#### Test dredge_pglmm() ####
## correct results
test_that("dredge_pglmm runs correctly", {

  res <- dredge_pglmm(formulaRE = log_onset ~ 1 + (1 | Species),
                      fixed = c("log_afr", "log_mass"),
                      data = senescence,
                      rank = "AICc",
                      family = "gaussian",
                      cov_ranef = list(Species = tree_ultra),
                      estimate = TRUE, std.err = FALSE, round = 3)

  ## structure
  expect_type(res, "list")
  expect_length(res, 7)

  ### dredge_table
  #### the content of dredge_table has already been check in test-internal-dredge_pglmm.R
  expect_s3_class(res$dredge_table, "data.frame")

  ### full_model_chr
  expect_type(res$full_model_chr, "character")
  expect_equal(res$full_model_chr, "log_onset ~ log_afr + log_mass + (1 | Species)")
  expect_length(res$full_model_chr, 1)

  ### random_display
  expect_type(res$random_display, "character")
  expect_equal(res$random_display, "(1 | Species)")
  expect_length(res$full_model_chr, 1)

  ### rank_name
  expect_equal(res$rank_name, "AICc")

  ### data_name
  expect_type(res$data_name, "character")
  expect_equal(res$data_name, "senescence")

  ### cov_ranef_name
  expect_type(res$cov_ranef_name, "character")
  expect_equal(res$cov_ranef_name, "list(Species = tree_ultra)")

  ### family_name
  expect_type(res$family_name, "character")
  expect_equal(res$family_name, "gaussian")
})

## errors when invalid arguments
test_that("dredge_pglmm errors to invalid arguments", {

  # ---- formulaRE ----
  ## not a formula
  expect_error(
    dredge_pglmm(formulaRE = "log_onset ~ 1 + (1 | Species)",
                 fixed = c("log_afr", "log_mass"),
                 data = senescence,
                 rank = "AICc",
                 family = "gaussian",
                 cov_ranef = list(Species = tree_ultra),
                 estimate = TRUE, std.err = FALSE, round = 3),
    "formulaRE must be an object of class formula")

  ## not a two-sided formula
  expect_error(
    dredge_pglmm(formulaRE = ~ 1 + (1 | Species),
                 fixed = c("log_afr", "log_mass"),
                 data = senescence,
                 rank = "AICc",
                 family = "gaussian",
                 cov_ranef = list(Species = tree_ultra),
                 estimate = TRUE, std.err = FALSE, round = 3),
    "formulaRE must be a two-sided formula")

  # ---- fixed ----
  ## not a character
  expect_error(
    dredge_pglmm(formulaRE = log_onset ~ 1 + (1 | Species),
                 fixed = 123,
                 data = senescence,
                 rank = "AICc",
                 family = "gaussian",
                 cov_ranef = list(Species = tree_ultra),
                 estimate = TRUE, std.err = FALSE, round = 3),
    "fixed must be a character string, or a vector of characters")

  ## not only character
  expect_error(
    dredge_pglmm(formulaRE = log_onset ~ 1 + (1 | Species),
                 fixed = c(123, "log_mass"),
                 data = senescence,
                 rank = "AICc",
                 family = "gaussian",
                 cov_ranef = list(Species = tree_ultra),
                 estimate = TRUE, std.err = FALSE, round = 3),
    "All variables in 'fixed' must be variables of 'data'")

  ## not a variable of data
  expect_error(
    dredge_pglmm(formulaRE = log_onset ~ 1 + (1 | Species),
                 fixed = c("log_mass", "not a variable of data"),
                 data = senescence,
                 rank = "AICc",
                 family = "gaussian",
                 cov_ranef = list(Species = tree_ultra),
                 estimate = TRUE, std.err = FALSE, round = 3),
    "All variables in 'fixed' must be variables of 'data'")

  # ---- data ----
  ## wrong data
  wrongdata <- "not a data.frame"
  expect_error(
    dredge_pglmm(formulaRE = log_onset ~ 1 + (1 | Species),
                 fixed = c("log_afr", "log_mass"),
                 data = wrongdata,
                 rank = "AICc",
                 family = "gaussian",
                 cov_ranef = list(Species = tree_ultra),
                 estimate = TRUE, std.err = FALSE, round = 3),
    "data must be provided and must be an existing data.frame")

  # ---- cov_ranef ----
  ## not a matrix or phylo object
  expect_error(
    dredge_pglmm(formulaRE = log_onset ~ 1 + (1 | Species),
                 fixed = c("log_afr", "log_mass"),
                 data = senescence,
                 rank = "AICc",
                 family = "gaussian",
                 cov_ranef = 123,
                 estimate = TRUE, std.err = FALSE, round = 3),
    "Each element of cov_ranef must be a matrix or a phylo object.")

  # ---- round ----
  ## not a numeric
  expect_error(
    dredge_pglmm(formulaRE = log_onset ~ 1 + (1 | Species),
                 fixed = c("log_afr", "log_mass"),
                 data = senescence,
                 rank = "AICc",
                 family = "gaussian",
                 cov_ranef = list(Species = tree_ultra),
                 estimate = TRUE, std.err = FALSE, round = "3"),
    "round must be a non-negative integer")

  ## negative
  expect_error(
    dredge_pglmm(formulaRE = log_onset ~ 1 + (1 | Species),
                 fixed = c("log_afr", "log_mass"),
                 data = senescence,
                 rank = "AICc",
                 family = "gaussian",
                 cov_ranef = list(Species = tree_ultra),
                 estimate = TRUE, std.err = FALSE, round = -3),
    "round must be a non-negative integer")

  ## not an integer
  expect_error(
    dredge_pglmm(formulaRE = log_onset ~ 1 + (1 | Species),
                 fixed = c("log_afr", "log_mass"),
                 data = senescence,
                 rank = "AICc",
                 family = "gaussian",
                 cov_ranef = list(Species = tree_ultra),
                 estimate = TRUE, std.err = FALSE, round = 3.2),
    "round must be a non-negative integer")
})

#### Test print.dredge.table() ####
test_that("print.dredge.table works", {

  res_TF <- dredge_pglmm(formulaRE = log_onset ~ 1 + (1 | Species),
                      fixed = c("log_afr", "log_mass"),
                      data = senescence,
                      rank = "AIC",
                      family = "gaussian",
                      cov_ranef = list(Species = tree_ultra),
                      estimate = TRUE, std.err = FALSE, round = 3)

  expect_invisible(print(res_TF))

  res_TT <- dredge_pglmm(formulaRE = log_onset ~ 1 + (1 | Species),
                         fixed = c("log_afr", "log_mass"),
                         data = senescence,
                         rank = "AICc",
                         family = "gaussian",
                         cov_ranef = list(Species = tree_ultra),
                         estimate = TRUE, std.err = TRUE, round = 3)

  expect_invisible(print(res_TT))

  res_TT <- dredge_pglmm(formulaRE = log_onset ~ 1 + (1 | Species),
                         fixed = c("log_afr", "log_mass"),
                         data = senescence,
                         rank = "AIC",
                         family = "gaussian",
                         cov_ranef = list(Species = tree_ultra),
                         estimate = FALSE, std.err = FALSE, round = 3)

  expect_invisible(print(res_TT))

  res_TT5 <- dredge_pglmm(formulaRE = log_onset ~ 1 + (1 | Species),
                         fixed = c("log_afr", "log_mass"),
                         data = senescence,
                         rank = "AICc",
                         family = "gaussian",
                         cov_ranef = list(Species = tree_ultra),
                         estimate = TRUE, std.err = TRUE, round = 5)

  expect_invisible(print(res_TT5))
})
