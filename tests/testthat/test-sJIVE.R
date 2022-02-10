test_that("sJIVE and sJIVE.predict works", {
  set.seed(1)
  train.x <- list(matrix(rnorm(300), ncol=20),
                  matrix(rnorm(200), ncol=20))
  train.y <- rnorm(20)

  train.fit <- sJIVE(X=train.x,Y=train.y,rankJ=1,rankA=c(1,1),eta=0.5)

  test.x <- list(matrix(rnorm(600), ncol=40),
                 matrix(rnorm(400), ncol=40))
  test.y <- rnorm(20)
  test.fit <- sJIVE.predict(train.fit, newdata = test.x)

  expect_equal(round(sum((test.y-test.fit$Ypred)^2),4), 71.3117)
})
