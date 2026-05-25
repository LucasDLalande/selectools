#### Test .title_box ####
test_that(".title_box produces expected title boxes", {
  res <- .title_box("HELLO WORLD", width = 10, border = "#")

  expected <- cat("\n\n", strrep("#", 10), "\n### HELLO WORLD ###\n", strrep("#", 10), "\n\n")

  expect_equal(res, expected)
})

#### Test .sep_line ####
test_that(".sep_line produces expected separation line", {
  res <- .sep_line(line_type = "-", width = 10)

  expected <- cat(strrep("-", 10), "\n\n")

  expect_equal(res, expected)
})

#### Test .first_upp ####
test_that(".first_upp capitalizes first letter of a character string", {
  res <- .first_upp("first letter should be capitalized")

  expected <- "First letter should be capitalized"

  expect_equal(res, expected)
})
