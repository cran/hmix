#' hmix
#'
#' @description hmix function segments the time series with k-means clustering, fits an HMM to model state transitions, and generates future predictions over a specified horizon. It evaluates model accuracy by calculating the Continuous Ranked Probability Score (CRPS) across multiple test points, producing error metrics that assess the model's predictive performance and robustness.
#'
#' @param ts A numeric vector representing the time series data.
#' @param horizon Integer. The prediction horizon, specifying how many future points to forecast.
#' @param centers Integer. Number of clusters for k-means clustering. Default: 10.
#' @param n_hidden Integer. Number of hidden states in the Hidden Markov Model. Default: 4.
#' @param seed Integer. Random seed for reproducibility. Default: 42.
#' @param n_tests Integer. Number of testing points for back-testing. Default: 20.
#' @param warmup Numeric. Proportion of the time series used as the warm-up period before testing. Default: 0.5.
#'
#' @return This function returns a list containing:
#' \itemize{
#'   \item model: The HMM model along with its estimated parameters.
#' \itemize{
#'      \item hmm_model: The object includes classified observations, initial HMM and trained HMM.
#'      \item pred_funs: Prediction functions for each point in horizon (rfun, dfun, pfun, qfun)
#'      }
#'   \item error_sets: A list of error metrics calculated for each testing point (at this time, CRPS).
#' }
#'
#' @import normalp
#' @import glogis
#' @import gld
#' @import edfun
#' @importFrom purrr map transpose map_dfr
#' @importFrom HMM initHMM baumWelch forward viterbiTraining
#' @importFrom mc2d rdirichlet
#' @importFrom cubature adaptIntegrate
#' @importFrom dplyr sample_n
#' @importFrom stats kmeans runif
#' @importFrom utils tail head
#'
#' @export
#'
#' @examples
#' \donttest{
#' # Example usage of hmix function:
#' result <- hmix(dummy_set$AMZN, horizon = 10, centers = 5, n_hidden = 3, n_tests = 2)
#' print(result$model)
#' print(result$error_sets)
#'
#' # Random sampling for each point in predicted horizon
#' result$model$pred_funs$t1$rfun(10)
#'
#' # ICDF for each point in horizon
#' result$model$pred_funs$t5$qfun(c(0, 1))
#'
#' # PDF for each point in horizon
#' result$model$pred_funs$t8$dfun(tail(ts))
#'
#' # CDF for each point in horizon
#' result$model$pred_funs$t10$pfun(tail(ts))
#' }
#'
hmix <- function(ts, horizon, centers = 10, n_hidden = 4, seed = 42, n_tests = 20, warmup = 0.5)
{
  error_sets <- NULL
  if(n_tests > 0)
  {
    n_length <- length(ts)
    testing_points <- unique(round(seq(warmup * n_length, n_length - horizon, length.out = n_tests)))
    if(length(testing_points) < n_tests){warning("adapting testing to available space"); n_tests <- length(testing_points)}

    models <- map(testing_points, ~ hmix_machine(ts = ts[1:.x], horizon = horizon, centers = centers, n_hidden = n_hidden, seed = seed))
    truths <- map(testing_points, ~ tail(ts[1:(.x + horizon)], horizon))

    pred_funs <- transpose(lapply(1:horizon, function(t) tryCatch(map(models, ~ .x$pred_funs[[t]]), error = function(e) NA)))

    crps_set <- t(mapply(function(bc) mapply(function(t) crps(truths, pred_funs, bc, t), t = 1:horizon), bc = 1:n_tests))
    colnames(crps_set) <- paste0("t", 1:horizon)
    rownames(crps_set) <- paste0("test_", 1:n_tests)

    error_sets <- list(crps_set = crps_set)
  }

  model <- hmix_machine(ts = ts, horizon = horizon, centers = centers, n_hidden = n_hidden, seed = seed)
  outcome <- list(model = model, error_sets = error_sets)
  return(outcome)
}



hmix_machine <- function(ts, horizon, centers = 10, n_hidden = 4, n_samples = 1000, seed = 42)
{
  set.seed(seed)

  dts <- tail(ts, -1)/head(ts, -1) - 1
  reframed <- smart_reframer(dts, horizon, horizon)

  groups <- kmeans(reframed, centers, nstart = 10)
  buckets <- map(1:centers, ~ as.data.frame(reframed[groups$cluster == .x,, drop = FALSE]))

  observations <- groups$cluster
  hidden_states <- paste0("s", 1:n_hidden)

  initial_transition_probs <- rdirichlet(n_hidden, rep(1, n_hidden))
  initial_emission_probs <- rdirichlet(n_hidden, rep(1, length(unique(observations))))
  initial_probs <- rdirichlet(1, rep(1, n_hidden))

  hmm <- initHMM(hidden_states, unique(observations), transProbs = initial_transition_probs, emissionProbs = initial_emission_probs, startProbs = initial_probs)

  trained_hmm <- tryCatch(baumWelch(hmm, observations, delta = 1E-9, maxIterations = 1000), error = function(e) viterbiTraining(hmm, observations, maxIterations = 1000, delta = 1E-9, pseudoCount = 1E-9))

  estimated_transition_probs <- NULL
  estimated_emission_probs <- NULL
  forward_probs <- NULL
  next_state_probs <- NULL
  next_observation_probs <- NULL

  if(is.list(trained_hmm))
  {
    estimated_transition_probs <- trained_hmm$hmm$transProbs
    estimated_emission_probs <- trained_hmm$hmm$emissionProbs

    log_forward_probs <- forward(trained_hmm$hmm, observations)
    forward_probs <- apply(log_forward_probs, 2, function(x) exp(x - max(x))/sum(exp(x - max(x))))###NUMERICALLY STABLE
    idx_space <- colSums(forward_probs)
    n_cols <- length(idx_space)
    idx <- sum(idx_space, na.rm = TRUE)
    if(n_cols > idx){warning("computation stability issues")}
    current_state_probs <- forward_probs[, idx]
    next_state_probs <- current_state_probs %*% trained_hmm$hmm$transProbs
    next_observation_probs <- next_state_probs %*% trained_hmm$hmm$emissionProbs

    if(all(!is.finite(next_observation_probs))){next_observation_probs <- rdirichlet(1, rep(1, centers)); colnames(next_observation_probs) <- 1:centers; warning("error with forward, untrained results")}
  }
  else {next_observation_probs <- rdirichlet(1, rep(1, centers)); colnames(next_observation_probs) <- 1:centers; warning("error during training, untrained results")}

  preds <- as.numeric(sample(as.numeric(colnames(next_observation_probs)), n_samples, replace = TRUE, prob = next_observation_probs))
  sampled_set <- map_dfr(preds, ~ sample_n(buckets[[.x]], 1))
  colnames(sampled_set) <- paste0("t", 1:horizon)

  sampled_cumprod <- as.data.frame(t(apply(sampled_set, 1, function(x) tail(ts, 1) * cumprod(x + 1))))
  if(all(ts > 0)){sampled_cumprod <- sampled_cumprod[apply(sampled_cumprod, 1, function(x) all(x > 0)), , drop = FALSE]}
  if(all(ts <= 0)){sampled_cumprod <- sampled_cumprod[apply(sampled_cumprod, 1, function(x) all(x <= 0)), , drop = FALSE]}
  if(all(ts%%1 == 0)){sampled_cumprod <- round(sampled_cumprod, 0)}

  pred_funs <- map(sampled_cumprod, ~ suppressWarnings(mixture(.x)))
  names(pred_funs) <- paste0("t", 1:horizon)

  hmm_model <- list(observations = observations, hmm = hmm, trained_hmm = trained_hmm)

  out <- list(hmm_model = hmm_model, pred_funs = pred_funs)
  return(out)
}

crps <- function(truth, pred_fun, backtesting_cycle, time)
{
  if(is.na(pred_fun[[backtesting_cycle]][time])){return(NA)}
  if(is.na(truth[[backtesting_cycle]][time])){return(NA)}
  ground <- truth[[backtesting_cycle]][time]
  eval_range <- pred_fun[[backtesting_cycle]][[time]]$qfun(c(0, 1))
  cdf <- pred_fun[[backtesting_cycle]][[time]]$pfun
  out <- adaptIntegrate(function(x) (cdf(x) - (x >= ground))^2, eval_range[1], eval_range[2])$integral
  return(out)
}

mixture <- function(x, alpha = 0.025)
{
  gnd_params <- tryCatch(unlist(paramp(x)[c(2, 4, 5)]), error = function(e) NA)
  if(is.numeric(gnd_params)){gnd_set <- qnormp(runif(length(x), alpha, 1-alpha), gnd_params[1], gnd_params[2], gnd_params[3])} else {gnd_set <- NA}

  glogis_params <- tryCatch(glogisfit(x)$parameters, error = function(e) NA)
  if(is.numeric(glogis_params)){glogis_set <- qglogis(runif(length(x), alpha, 1-alpha), glogis_params[1], glogis_params[2], glogis_params[3])} else {glogis_set <- NA}

  gl_params <- tryCatch(fit.fkml(x)$lambda, error = function(e) NA)
  if(is.numeric(gl_params)){gl_set <- qgl(runif(length(x), alpha, 1-alpha), gl_params[1],  gl_params[2], gl_params[3], gl_params[4])} else {gl_set <- NA}

  values <- c(x, gnd_set, glogis_set, gl_set)
  values <- values[is.finite(values)]
  if(all(x >= 0)){values <- values[values >= 0]}
  if(all(x < 0)){values <- values[values < 0]}

  pred_fun <- edfun(values, cut = 0)

  return(pred_fun)
}


smart_reframer <- function(ts, seq_len, stride)
{
  n_length <- length(ts)
  if(seq_len > n_length | stride > n_length){stop("vector too short for sequence length or stride")}
  if(n_length%%seq_len > 0){ts <- tail(ts, - (n_length%%seq_len))}
  n_length <- length(ts)
  idx <- seq(from = 1, to = (n_length - seq_len + 1), by = 1)
  reframed <- t(sapply(idx, function(x) ts[x:(x+seq_len-1)]))
  if(seq_len == 1){reframed <- t(reframed)}
  idx <- rev(seq(nrow(reframed), 1, - stride))
  reframed <- reframed[idx,,drop = FALSE]
  colnames(reframed) <- paste0("t", 1:seq_len)
  return(reframed)
}
