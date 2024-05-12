test_that("entropy",{
  bn <- loadHuginNet("../models/alarm.net")
  target <- c("VENTALV", "VENTTUBE")
  query <- c("ARTCO2")
  marg_ent <- ent(bn, target, "marginal")
  expect_equal(names(marg_ent), target)
  expect_lt(marg_ent[1] - 0.889866244, 0.0000001 )
  expect_lt(marg_ent[2] - 0.773995889, 0.0000001 )
  jnt_ent <- ent(bn, target, "joint")
  expect_lt(jnt_ent - 1.404494867, 0.0000001 )
  
})

test_that("calc_ent",{
  expect_equal(entropy(c(0.6,0.4)),0.673011667)
  expect_equal(entropy(c(1,0,0)),0)
})
