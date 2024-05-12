test_that("entropy_check", {
  bn <- loadHuginNet("../models/alarm.net")
  target <- c("VENTALV", "VENTTUBE")
  query <- c("ARTCO2")
  expect_equal(entropy_check(bn, target, 0.9), TRUE)
  expect_equal(entropy_check(bn, target, 0.8), FALSE)
})

