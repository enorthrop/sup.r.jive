test_that("sJIVE works", {
  data(SimData.norm)
  train.fit <- sJIVE(X=SimData.norm$X,
                     Y=SimData.norm$Y,
                     rankJ=1,rankA=c(1,1),eta=0.5)
  train.mse <- round(sum((SimData.norm$Y-train.fit$fittedY)^2),3)

  expect_equal(train.mse, 11.132)
  train.fit
  summary(train.fit)
})

test_that("sJIVE prediction works", {
  data(SimData.norm)
  train.fit <- sJIVE(X=SimData.norm$X,
                     Y=SimData.norm$Y,
                     rankJ=1,rankA=c(1,1),eta=0.5)

  n <- 300
  p <- 40 #Don't change p unless SimData Changes
  withr::with_seed( 1,
  test.x <- list(matrix(rnorm(n*p), ncol=n),
                 matrix(rnorm(n*p), ncol=n))
  )
  withr::with_seed( 2,
       test.y <- rnorm(20)
  )
  test.fit <- predict(train.fit, newdata = test.x)
  test.mse <- round(sum((test.y-test.fit$Ypred)^2),3)
  expect_equal(test.mse, 365.342)
})


test_that("plot sJIVE is error free", {
  data(SimData.norm)
  fit <- sJIVE(X=SimData.norm$X,
                     Y=SimData.norm$Y,
                     rankJ=1,rankA=c(1,1),eta=0.5)
  plotHeatmap(fit)
  plotVarExplained(fit)
  plotFittedValues(fit)
})

