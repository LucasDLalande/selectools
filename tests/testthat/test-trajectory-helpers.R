#### Test .collapse_deparse() ####
test_that(".collapse_deparse correctly convert R objects to strings", {
  x <- response ~ var1 + var2 + (1|random1) + (1|random2)
  random <- findbars(x)

  res <- .collapse_deparse(random)

  expected <- paste(deparse(random), collapse="")

  expect_equal(res, expected)



  y <- quote(1|random1)

  res <- .collapse_deparse(y)

  expected <- paste(deparse(y), collapse="")

  expect_equal(res, expected)
})

#### Test .update_call_data() ####
test_that(".update_call_data correctly update the data name in the model call and errors if not 'lm' or 'lmerMod' class objects", {

  set.seed(123)

  # standardized data to avoid scaling warnings
  df <- data.frame(x = rnorm(100, mean = 0, sd = 1),
                   y = rnorm(100, mean = 0, sd = 1),
                   z = rnorm(100, mean = 0, sd = 1),
                   group = factor(c(rep("A", 25), rep("B", 25), rep("C", 25), rep("D", 25))))

  mod_lm <- lm(y ~ x + z, data = df)
  mod_lmer <- lmer(y ~ x + z + (1|group), data = df, REML = FALSE)

  mydata <- substitute(df)

  updated_lm <- .update_call_data(mod_lm, mydata)
  updated_lmer <- .update_call_data(mod_lmer, mydata)

  expect_equal(updated_lm$call$data, mydata)
  expect_equal(updated_lmer@call$data, mydata)

  expect_error(.update_call_data(123, mydata),
               "model must be of class 'lm' or 'lmerMod'")
})

#### Test .set_call_env() ####
test_that(".set_call_env correctly update evaluation environment and errors if not 'lm' or 'lmer' class objects", {

  set.seed(123)

  # standardized data to avoid scaling warnings
  df <- data.frame(x = rnorm(100, mean = 0, sd = 1),
                   y = rnorm(100, mean = 0, sd = 1),
                   z = rnorm(100, mean = 0, sd = 1),
                   group = factor(c(rep("A", 25), rep("B", 25), rep("C", 25), rep("D", 25))))

  mod_lm <- lm(y ~ x + z, data = df)
  mod_lmer <- lmer(y ~ x + z + (1|group), data = df, REML = FALSE)

  env <- new.env(parent = parent.frame())

  updated_lm <- .set_call_env(mod_lm, env)
  updated_lmer <- .set_call_env(mod_lmer, env)

  expect_equal(environment(updated_lm$call), env)
  expect_equal(environment(updated_lmer@call), env)

  expect_error(.set_call_env(123, env),
               "model must be of class 'lm' or 'lmerMod'")
})


