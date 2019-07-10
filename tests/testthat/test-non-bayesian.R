context("Non-Bayesian")

test_that("preorder works", {
  expect_warning(preorder(y), NA)
  expect_warning(preorder(y, factors = 2), NA)
  expect_warning(preorder(y, factors = 3), NA)
  expect_warning(preorder(y, factors = 2, type = "dynamic"), NA)
  expect_warning(preorder(y, factors = 2, type = "dynamic", transload = exp), NA)
})

test_that("findrestrict works", {
  # default relto
  expect_warning(findrestrict(y, factors = 2), NA)
  expect_warning(findrestrict(y, factors = 3), NA)
  expect_warning(findrestrict(y, factors = 2, transload = exp), NA)
  for (relto in c("all", "none", "current")) {
    # non-default relto
    expect_warning(findrestrict(y, factors = 2, relto = relto), NA)
    expect_warning(findrestrict(y, factors = 3, relto = relto), NA)
    expect_warning(findrestrict(y, factors = 2, relto = relto, transload = exp), NA)
  }
})
