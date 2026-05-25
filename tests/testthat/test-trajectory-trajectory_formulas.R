#### Test trajectory_formulas() ####
## correct results
test_that("trajectory_formulas correctly creates trajectory formulas with random effects and interactions", {

  set.seed(123)

  # standardized data to avoid scaling warnings
  df <- data.frame(x = rnorm(100, mean = 0, sd = 1),
                   y = rnorm(100, mean = 0, sd = 1),
                   z = rnorm(100, mean = 0, sd = 1),
                   group = factor(c(rep("A", 25), rep("B", 25), rep("C", 25), rep("D", 25))))

  formula <- y ~ x + (1|group)

  res <- trajectory_formulas(formula = formula, trajectory_trait = "z", interactions = "x", data = df)

  res <- list(
    constant = res$formulas$non_threshold_formulas$constant,
    linear = res$formulas$non_threshold_formulas$linear,
    quadratic = res$formulas$non_threshold_formulas$quadratic,
    factor = res$formulas$non_threshold_formulas$factor,
    thresholdCS = res$formulas$threshold_formulas$thresholdCS,
    thresholdSC = res$formulas$threshold_formulas$thresholdSC,
    thresholdSS = res$formulas$threshold_formulas$thresholdSS)

  expected <- list(
    constant = y ~ x + (1|group),
    linear = y ~ x + z + (1|group) + x:z,
    quadratic = y ~ x + z + I(z^2) + (1|group) + x:z + x:I(z^2),
    factor = y ~ x + factor(z) + (1|group) + x:factor(z),
    thresholdCS = y ~ x + z.2 + (1|group) + x:z.2,
    thresholdSC = y ~ x + z.1 + (1|group) + x:z.1,
    thresholdSS = y ~ x + z.1 + z.2 + (1|group) + x:z.1 + x:z.2)

  expect_equal(res, expected, ignore_attr = TRUE)
})

test_that("trajectory_formulas correctly creates trajectory formulas without random effects and interactions", {

  set.seed(123)

  # standardized data to avoid scaling warnings
  df <- data.frame(x = rnorm(100, mean = 0, sd = 1),
                   y = rnorm(100, mean = 0, sd = 1),
                   z = rnorm(100, mean = 0, sd = 1),
                   group = factor(c(rep("A", 25), rep("B", 25), rep("C", 25), rep("D", 25))))

  formula <- y ~ x

  res <- trajectory_formulas(formula = formula, trajectory_trait = "z", data = df)

  res <- list(
    constant = res$formulas$non_threshold_formulas$constant,
    linear = res$formulas$non_threshold_formulas$linear,
    quadratic = res$formulas$non_threshold_formulas$quadratic,
    factor = res$formulas$non_threshold_formulas$factor,
    thresholdCS = res$formulas$threshold_formulas$thresholdCS,
    thresholdSC = res$formulas$threshold_formulas$thresholdSC,
    thresholdSS = res$formulas$threshold_formulas$thresholdSS)

  expected <- list(
    constant = y ~ x,
    linear = y ~ x + z,
    quadratic = y ~ x + z + I(z^2),
    factor = y ~ x + factor(z),
    thresholdCS = y ~ x + z.2,
    thresholdSC = y ~ x + z.1,
    thresholdSS = y ~ x + z.1 + z.2)

  expect_equal(res, expected, ignore_attr = TRUE)
})

test_that("trajectory_formulas correctly creates trajectory formulas with several interaction terms", {

  set.seed(123)

  # standardized data to avoid scaling warnings
  df <- data.frame(x = rnorm(100, mean = 0, sd = 1),
                   y = rnorm(100, mean = 0, sd = 1),
                   z = rnorm(100, mean = 0, sd = 1),
                   a = rnorm(100, mean = 0, sd = 1),
                   group = factor(c(rep("A", 25), rep("B", 25), rep("C", 25), rep("D", 25))))

  formula <- y ~ x + a + (1|group)

  res <- trajectory_formulas(formula = formula, trajectory_trait = "z", interactions = c("x", "a"), data = df)

  res <- list(
    constant = res$formulas$non_threshold_formulas$constant,
    linear = res$formulas$non_threshold_formulas$linear,
    quadratic = res$formulas$non_threshold_formulas$quadratic,
    factor = res$formulas$non_threshold_formulas$factor,
    thresholdCS = res$formulas$threshold_formulas$thresholdCS,
    thresholdSC = res$formulas$threshold_formulas$thresholdSC,
    thresholdSS = res$formulas$threshold_formulas$thresholdSS)

  expected <- list(
    constant = y ~ x + a + (1|group),
    linear = y ~ x + a + z + (1|group) + x:z + a:z,
    quadratic = y ~ x + a + z + I(z^2) + (1|group) + x:z + a:z + x:I(z^2) + a:I(z^2),
    factor = y ~ x + a + factor(z) + (1|group) + x:factor(z) + a:factor(z),
    thresholdCS = y ~ x + a + z.2 + (1|group) + x:z.2 + a:z.2,
    thresholdSC = y ~ x + a + z.1 + (1|group) + x:z.1 + a:z.1,
    thresholdSS = y ~ x + a + z.1 + z.2 + (1|group) + x:z.1 + a:z.1 + x:z.2 + a:z.2)

  expect_equal(res, expected, ignore_attr = TRUE)
})

test_that("trajectory_formulas correctly creates specific trajectories", {

  set.seed(123)

  # standardized data to avoid scaling warnings
  df <- data.frame(x = rnorm(100, mean = 0, sd = 1),
                   y = rnorm(100, mean = 0, sd = 1),
                   z = rnorm(100, mean = 0, sd = 1),
                   a = rnorm(100, mean = 0, sd = 1),
                   group = factor(c(rep("A", 25), rep("B", 25), rep("C", 25), rep("D", 25))))

  formula <- y ~ x + a + (1|group)

  ## specified trajectories
  res1 <- trajectory_formulas(formula = formula, trajectory_trait = "z", interactions = c("x", "a"), data = df, trajectories = c("linear", "thresholdSC"))

  res1 <- list(
    linear = res1$formulas$non_threshold_formulas$linear,
    thresholdSC = res1$formulas$threshold_formulas$thresholdSC)

  expected1 <- list(
    linear = y ~ x + a + z + (1|group) + x:z + a:z,
    thresholdSC = y ~ x + a + z.1 + (1|group) + x:z.1 + a:z.1)

  expect_equal(res1, expected1, ignore_attr = TRUE)

  ## all threshold trajectories
  res2 <- trajectory_formulas(formula = formula, trajectory_trait = "z", interactions = c("x", "a"), data = df, trajectories = "threshold")

  res2 <- list(
    thresholdCS = res2$formulas$threshold_formulas$thresholdCS,
    thresholdSC = res2$formulas$threshold_formulas$thresholdSC,
    thresholdSS = res2$formulas$threshold_formulas$thresholdSS)

  expected2 <- list(
    thresholdCS = y ~ x + a + z.2 + (1|group) + x:z.2 + a:z.2,
    thresholdSC = y ~ x + a + z.1 + (1|group) + x:z.1 + a:z.1,
    thresholdSS = y ~ x + a + z.1 + z.2 + (1|group) + x:z.1 + a:z.1 + x:z.2 + a:z.2)

  expect_equal(res2, expected2, ignore_attr = TRUE)

  ## all non-threshold trajectories
  res3 <- trajectory_formulas(formula = formula, trajectory_trait = "z", interactions = c("x", "a"), data = df, trajectories = "non-threshold")

  res3 <- list(
    constant = res3$formulas$non_threshold_formulas$constant,
    linear = res3$formulas$non_threshold_formulas$linear,
    quadratic = res3$formulas$non_threshold_formulas$quadratic,
    factor = res3$formulas$non_threshold_formulas$factor)

  expected3 <- list(
    constant = y ~ x + a + (1|group),
    linear = y ~ x + a + z + (1|group) + x:z + a:z,
    quadratic = y ~ x + a + z + I(z^2) + (1|group) + x:z + a:z + x:I(z^2) + a:I(z^2),
    factor = y ~ x + a + factor(z) + (1|group) + x:factor(z) + a:factor(z))

  expect_equal(res3, expected3, ignore_attr = TRUE)

  ## all trajectories
  res4 <- trajectory_formulas(formula = formula, trajectory_trait = "z", interactions = c("x", "a"), data = df, trajectories = "all")

  res4 <- list(
    constant = res4$formulas$non_threshold_formulas$constant,
    linear = res4$formulas$non_threshold_formulas$linear,
    quadratic = res4$formulas$non_threshold_formulas$quadratic,
    factor = res4$formulas$non_threshold_formulas$factor,
    thresholdCS = res4$formulas$threshold_formulas$thresholdCS,
    thresholdSC = res4$formulas$threshold_formulas$thresholdSC,
    thresholdSS = res4$formulas$threshold_formulas$thresholdSS)

  expected4 <- list(
    constant = y ~ x + a + (1|group),
    linear = y ~ x + a + z + (1|group) + x:z + a:z,
    quadratic = y ~ x + a + z + I(z^2) + (1|group) + x:z + a:z + x:I(z^2) + a:I(z^2),
    factor = y ~ x + a + factor(z) + (1|group) + x:factor(z) + a:factor(z),
    thresholdCS = y ~ x + a + z.2 + (1|group) + x:z.2 + a:z.2,
    thresholdSC = y ~ x + a + z.1 + (1|group) + x:z.1 + a:z.1,
    thresholdSS = y ~ x + a + z.1 + z.2 + (1|group) + x:z.1 + a:z.1 + x:z.2 + a:z.2)

  expect_equal(res4, expected4, ignore_attr = TRUE)
})

test_that("trajectory_formulas returns objects of correct attributes and class", {

  set.seed(123)

  # standardized data to avoid scaling warnings
  df <- data.frame(x = rnorm(100, mean = 0, sd = 1),
                   y = rnorm(100, mean = 0, sd = 1),
                   z = rnorm(100, mean = 0, sd = 1),
                   a = rnorm(100, mean = 0, sd = 1),
                   group = factor(c(rep("A", 25), rep("B", 25), rep("C", 25), rep("D", 25))))

  formula <- y ~ x + a + (1|group)

  res <- trajectory_formulas(formula = formula, trajectory_trait = "z", interactions = c("x", "a"), data = df, trajectories = c("linear", "thresholdSC"))

  expect_equal(class(res), "trajectory.formulas")
  expect_equal(attributes(res$formulas[[1]][[1]])$trajectory_type, "linear")
  expect_equal(attributes(res$formulas[[2]][[1]])$trajectory_type, "thresholdSC")
})

## errors when invalid arguments
test_that("trajectory_formulas errors when invalid arguments", {

  set.seed(123)

  # standardized data to avoid scaling warnings
  df <- data.frame(x = rnorm(100, mean = 0, sd = 1),
                   y = rnorm(100, mean = 0, sd = 1),
                   z = rnorm(100, mean = 0, sd = 1),
                   a = rnorm(100, mean = 0, sd = 1),
                   group = factor(c(rep("A", 25), rep("B", 25), rep("C", 25), rep("D", 25))))

  formula <- y ~ x + a + (1|group)

  wrongdata <- "not a data.frame"

  # ---- data ----
  expect_error(trajectory_formulas(formula = formula, trajectory_trait = "z", interactions = c("x", "a"), data = wrongdata),
               "data must be provided and must be an existing data.frame")

  # ---- formula ----
  ## not a formula
  expect_error(trajectory_formulas(formula = 123, trajectory_trait = "z", interactions = c("x", "a"), data = df),
               "formula must be an object of class formula")

  ## not a two-sided formula
  expect_error(trajectory_formulas(formula = ~ x + a, trajectory_trait = "z", interactions = c("x", "a"), data = df),
               "formula must be a two-sided formula")

  # ---- trajectory_trait ----
  ## not a single character string
  expect_error(trajectory_formulas(formula = formula, trajectory_trait = c("z", "y"), interactions = c("x", "a"), data = df),
               "trajectory_trait must be a single character string")

  ## not a variable from data
  expect_error(trajectory_formulas(formula = formula, trajectory_trait = "b", interactions = c("x", "a"), data = df),
               "trajectory_trait must be a variable of 'data'")

  # ---- interactions ----
  ## not a character vector or NULL
  expect_error(trajectory_formulas(formula = formula, trajectory_trait = "z", interactions = df$x, data = df),
               "interactions must be NULL or a character vector")

  ## not a variable from data
  expect_error(trajectory_formulas(formula = formula, trajectory_trait = "z", interactions = "b", data = df),
               "some interactions are not variables in 'data'")

  # ---- trajectories ----
  ## errors if not from allowed trajectory names
  expect_error(trajectory_formulas(formula = formula, trajectory_trait = "z", interactions = c("x", "a"), data = df, trajectories = "wrong_traj"),
               "Unknown trajectories")
})

#### Test print.trajectory.formulas() ####
test_that("print.trajectory.formulas works", {

  set.seed(123)

  # standardized data to avoid scaling warnings
  df <- data.frame(x = rnorm(100, mean = 0, sd = 1),
                   y = rnorm(100, mean = 0, sd = 1),
                   z = rnorm(100, mean = 0, sd = 1),
                   a = rnorm(100, mean = 0, sd = 1),
                   group = factor(c(rep("A", 25), rep("B", 25), rep("C", 25), rep("D", 25))))

  formula <- y ~ x + a + (1|group)

  res <- trajectory_formulas(formula = formula, trajectory_trait = "z", interactions = c("x", "a"), data = df)

  expect_invisible(print(res))
})
