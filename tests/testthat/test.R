test_that("hmix returns a list with expected structure", {
  result <- hmix(dummy_set$AMZN, horizon = 5, centers = 4, n_hidden = 4, seed = 123, n_tests = 3, warmup = 0.3)

  expect_type(result, "list")
  expect_named(result, c("model", "error_sets"))
  expect_type(result$model, "list")

  model <- result$model

  expect_type(model, "list")
  expect_named(model, c('hmm_model', 'pred_funs'))
  expect_type(model$hmm_model, "list")
  expect_named(model$hmm_model, c('observations', 'hmm', 'trained_hmm'))
  expect_named(model$pred_funs, c('t1', 't2', 't3', 't4', 't5'))
  expect_true(is.numeric(model$pred_funs$t1$rfun(5)) & length(model$pred_funs$t1$rfun(5)) == 5)
  expect_type(model$pred_funs$t1$qfun(c(0, 1)), "double")
  expect_type(model$pred_funs$t1$pfun(tail(ts)), "double")
  expect_type(model$pred_funs$t1$dfun(tail(ts)), "double")
  expect_true(is.numeric(model$pred_funs$t5$rfun(5)) & length(model$pred_funs$t5$rfun(5)) == 5)
  expect_type(model$pred_funs$t5$qfun(c(0, 1)), "double")
  expect_type(model$pred_funs$t5$pfun(tail(ts)), "double")
  expect_type(model$pred_funs$t5$dfun(tail(ts)), "double")

  expect_type(result$error_sets, "list")
  crps_set <- result$error_sets$crps_set
  expect_type(crps_set, "double")
  expect_true(all(crps_set >= 0))
  expect_true(prod(dim(crps_set)) == 15)
  })
