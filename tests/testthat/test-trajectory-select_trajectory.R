#### Test select_trajectory() ####
## correct results
test_that("select_trajectory correctly returns the recap table, final trajectory and summaries", {

  set.seed(123)

  # standardized data to avoid scaling warnings
  data <- data.frame(x = rnorm(100, mean = 0, sd = 1),
                     y = rnorm(100, mean = 0, sd = 1),
                     z = rnorm(100, mean = 0, sd = 1),
                     group = factor(c(rep("A", 25), rep("B", 25), rep("C", 25), rep("D", 25))))

  traj <- c("constant", "linear", "thresholdCS")

  res <- select_trajectory(formula = y ~ x + (1|group),
                           trajectory_trait = "z",
                           interactions = "x",
                           data = data,
                           trajectories = traj,
                           start = -0.8, end = 0.8, iter = 0.2,
                           rank = "AICc",
                           delta_criterion = 2,
                           adjust_threshold = FALSE,
                           lmer_control = NULL)

  res_adj <- select_trajectory(formula = y ~ x + (1|group),
                               trajectory_trait = "z",
                               interactions = "x",
                               data = data,
                               trajectories = traj,
                               start = -0.8, end = 0.8, iter = 0.2,
                               rank = "AICc",
                               delta_criterion = 4,
                               adjust_threshold = TRUE,
                               lmer_control = NULL)

  # ---- test expected output structure ----
  expect_equal(class(res), "best.trajectory")
  expect_type(res, "list")
  expect_length(res, 5)

  ## best_table
  expect_s3_class(res$best_table, "data.frame")
  expect_equal(res$best_table$trajectory, traj)
  expect_equal(nrow(res$best_table), length(traj))

  ## final_traj
  expect_s3_class(res$final_traj, "data.frame")

  ## final_traj_name
  expect_type(res$final_traj_name, "character")

  ## summaries
  expect_type(res$summaries, "list")
  expect_equal(names(res$summaries), traj)
  expect_equal(length(res$summaries), length(traj))
  expect_true(all(vapply(res$summaries, function(x) inherits(x, "summary.merMod"), logical(1))))

  ## final_summary
  expect_type(res$final_summary, "list")
  expect_equal(names(res$final_summary), res$final_traj$trajectory)
  expect_equal(length(res$final_summary), nrow(res$final_traj))
  expect_true(all(vapply(res$final_summary, function(x) inherits(x, "summary.merMod"), logical(1))))

  # ---- test results coherence ----
  ## ajust_threshold correctly increment df and rank
  expect_equal(res$best_table$df[grepl("threshold", res$best_table$trajectory)] + 1, res_adj$best_table$df[grepl("threshold", res_adj$best_table$trajectory)])
  expect_equal(res$best_table$AICc[grepl("threshold", res$best_table$trajectory)] + 2, res_adj$best_table$AICc[grepl("threshold", res_adj$best_table$trajectory)])

  ## final_traj have only models under 2 delta AICc
  expect_true(all(res$final_traj$delta <= 2))
  expect_true(all(res_adj$final_traj$delta <= 4))

  ## final_traj_name is a paste of trajectories in final_traj
  expect_equal(paste(res$final_traj$trajectory, collapse = ", "), res$final_traj_name)
})

## errors when invalid argument
test_that("select_rajectory errors when invalid arguments", {

  set.seed(123)

  # standardized data to avoid scaling warnings
  data <- data.frame(x = rnorm(100, mean = 0, sd = 1),
                     y = rnorm(100, mean = 0, sd = 1),
                     z = rnorm(100, mean = 0, sd = 1),
                     group = factor(c(rep("A", 25), rep("B", 25), rep("C", 25), rep("D", 25))))

  traj <- c("constant", "linear", "thresholdCS")

  res <- select_trajectory(formula = y ~ x + (1|group),
                           trajectory_trait = "z",
                           interactions = "x",
                           data = data,
                           trajectories = traj,
                           start = -0.8, end = 0.8, iter = 0.2,
                           rank = "AICc",
                           delta_criterion = 2,
                           adjust_threshold = FALSE,
                           lmer_control = NULL)

  # ---- start/end/iter ----
  ## not numeric values
  ### start
  expect_error(
    select_trajectory(formula = y ~ x + (1|group),
                      trajectory_trait = "z",
                      interactions = "x",
                      data = data,
                      trajectories = traj,
                      start = "-0.8", end = 0.8, iter = 0.2,
                      rank = "AICc",
                      delta_criterion = 2,
                      adjust_threshold = FALSE,
                      lmer_control = NULL),
    "start, end and iter must be single numeric values")

  ### end
  expect_error(
    select_trajectory(formula = y ~ x + (1|group),
                      trajectory_trait = "z",
                      interactions = "x",
                      data = data,
                      trajectories = traj,
                      start = -0.8, end = "0.8", iter = 0.2,
                      rank = "AICc",
                      delta_criterion = 2,
                      adjust_threshold = FALSE,
                      lmer_control = NULL),
    "start, end and iter must be single numeric values")

  ### iter
  expect_error(
    select_trajectory(formula = y ~ x + (1|group),
                      trajectory_trait = "z",
                      interactions = "x",
                      data = data,
                      trajectories = traj,
                      start = -0.8, end = 0.8, iter = "0.2",
                      rank = "AICc",
                      delta_criterion = 2,
                      adjust_threshold = FALSE,
                      lmer_control = NULL),
    "start, end and iter must be single numeric values")

  ## not a single value
  ### start
  expect_error(
    select_trajectory(formula = y ~ x + (1|group),
                      trajectory_trait = "z",
                      interactions = "x",
                      data = data,
                      trajectories = traj,
                      start = c(-0.8, -0.7), end = 0.8, iter = 0.2,
                      rank = "AICc",
                      delta_criterion = 2,
                      adjust_threshold = FALSE,
                      lmer_control = NULL),
    "start, end and iter must be single numeric values")

  ### end
  expect_error(
    select_trajectory(formula = y ~ x + (1|group),
                      trajectory_trait = "z",
                      interactions = "x",
                      data = data,
                      trajectories = traj,
                      start = -0.8, end = c(0.8, 0.7), iter = 0.2,
                      rank = "AICc",
                      delta_criterion = 2,
                      adjust_threshold = FALSE,
                      lmer_control = NULL),
    "start, end and iter must be single numeric values")

  ### iter
  expect_error(
    select_trajectory(formula = y ~ x + (1|group),
                      trajectory_trait = "z",
                      interactions = "x",
                      data = data,
                      trajectories = traj,
                      start = -0.8, end = 0.8, iter = c(0.2, 0.3),
                      rank = "AICc",
                      delta_criterion = 2,
                      adjust_threshold = FALSE,
                      lmer_control = NULL),
    "start, end and iter must be single numeric values")

  # ---- delta_criterion ----
  ## not a numeric
  expect_error(
    select_trajectory(formula = y ~ x + (1|group),
                      trajectory_trait = "z",
                      interactions = "x",
                      data = data,
                      trajectories = traj,
                      start = -0.8, end = 0.8, iter = 0.2,
                      rank = "AICc",
                      delta_criterion = "2",
                      adjust_threshold = FALSE,
                      lmer_control = NULL),
    "delta_criterion must be a single non-negative numeric value")

  ## not a single numeric
  expect_error(
    select_trajectory(formula = y ~ x + (1|group),
                      trajectory_trait = "z",
                      interactions = "x",
                      data = data,
                      trajectories = traj,
                      start = -0.8, end = 0.8, iter = 0.2,
                      rank = "AICc",
                      delta_criterion = c(2, 4),
                      adjust_threshold = FALSE,
                      lmer_control = NULL),
    "delta_criterion must be a single non-negative numeric value")

  ## not a non-negative numeric
  expect_error(
    select_trajectory(formula = y ~ x + (1|group),
                      trajectory_trait = "z",
                      interactions = "x",
                      data = data,
                      trajectories = traj,
                      start = -0.8, end = 0.8, iter = 0.2,
                      rank = "AICc",
                      delta_criterion = -2,
                      adjust_threshold = FALSE,
                      lmer_control = NULL),
    "delta_criterion must be a single non-negative numeric value")
})
