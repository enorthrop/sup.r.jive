

test_that("JIVEpred works with normal outcome", {
  data(SimData.norm)
  train.fit <- JIVE.pred(X=SimData.norm$X,
                     Y=SimData.norm$Y, family = "gaussian",
                     rankJ=1,rankI=c(1,1))
  train.mse <- round(sum((SimData.norm$Y-train.fit$mod.fit$fitted.values)^2),3)

  n <- 300
  p <- 40 #Don't change p unless SimData Changes
  withr::with_seed( 1,
     test.x <- list(matrix(rnorm(n*p), ncol=n),
         matrix(rnorm(n*p), ncol=n)))
  withr::with_seed( 2, test.y <- rnorm(20))
  test.fit <- predict(train.fit, newdata = test.x)
  test.mse <- round(sum((test.y-test.fit$Ypred)^2),3)

  expect_equal(train.mse, 8.607)
  expect_equal(test.mse, 366.887)
})


