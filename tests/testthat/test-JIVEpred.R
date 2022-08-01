

test_that("JIVEpred works with normal outcome", {
  data(SimData.norm)
  train.fit <- JIVE.pred(X=SimData.norm$X,
                     Y=SimData.norm$Y, family = "gaussian",
                     rankJ=1,rankA=c(1,1))
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
  train.fit
  summary(train.fit)
})

test_that("plot JIVEpred is error free - normal", {
  data(SimData.norm)
  fit <- JIVE.pred(X=SimData.norm$X,
              Y=SimData.norm$Y, family = "gaussian",
              rankJ=1,rankA=c(1,1))
  plotHeatmap(fit)
  plotVarExplained(fit)
  plotFittedValues(fit)
})

test_that("JIVEpred works with binary outcome", {
  data(SimData.bin)
  train.fit <- JIVE.pred(X=SimData.bin$X,
                         Y=SimData.bin$Y, family = "binomial",
                         rankJ=1,rankA=c(1,1))
  train.mse <- round(sum((SimData.bin$Y-train.fit$mod.fit$fitted.values)^2),3)

  n <- 300
  p <- 30 #Don't change p unless SimData Changes
  withr::with_seed( 1,
                    test.x <- list(matrix(rnorm(n*p), ncol=n),
                                   matrix(rnorm(n*p), ncol=n)))
  withr::with_seed( 2, test.y <- rnorm(20))
  test.fit <- predict(train.fit, newdata = test.x)
  test.mse <- round(sum((test.y-test.fit$Ypred)^2),3)

  expect_equal(train.mse, 4.136)
  expect_equal(test.mse, 332.003)
  train.fit
  summary(train.fit)
})

test_that("plot JIVEpred is error free - binary", {
  data(SimData.bin)
  fit <- JIVE.pred(X=SimData.bin$X,
                   Y=SimData.bin$Y, family = "binomial",
                   rankJ=1,rankA=c(1,1))
  plotHeatmap(fit)
  suppressWarnings(plotVarExplained(fit))
  plotFittedValues(fit)
})
