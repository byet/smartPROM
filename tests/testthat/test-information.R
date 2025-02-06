
test_that("entropy",{
  bn <- loadHuginNet("../testmodels/alarm.net")
  target <- c("VENTALV", "VENTTUBE")
  query <- c("ARTCO2")
  post1 <- querygrain(bn, "VENTALV")
  expect_lt(dist_entropy(post1) - 0.889866244, 0.0000001 )
  post2 <- querygrain(bn, "VENTTUBE")
  expect_lt(dist_entropy(post2) - 0.773995889, 0.0000001 )
  post3 <- querygrain(bn, target,type="joint")
  dist_entropy(post3)
  expect_lt(dist_entropy(post3) - 1.404494867, 0.0000001 )
  
})

test_that("calc_ent",{
  expect_equal(entropy(c(0.6,0.4)),0.673011667)
  expect_equal(entropy(c(1,0,0)),0)
})
