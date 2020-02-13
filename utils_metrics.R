
MetricsPositive <- function(actual, predicted){
  return(sum(actual))
}
MetricsPositive(actual = c(1, 1, 1, 1), predicted = c(1, 0, 0, 0))

MetricsNegative <- function(actual, predicted){
  return(sum(1-actual))
}
MetricsNegative(actual = c(1, 0, 1, 1), predicted = c(1, 0, 0, 0))

MetricsTruePositive <- function(actual, predicted){
  idx_positives <- which(actual == 1)
  return(sum(predicted[idx_positives]))
}
MetricsTruePositive(actual = c(1, 1, 1, 1), predicted = c(1, 0, 0, 0))

MetricsFalsePositive <- function(actual, predicted){
  idx_negatives <- which(actual == 0)
  return(sum(predicted[idx_negatives]))
}
MetricsFalsePositive(actual = c(0, 0, 1, 1), predicted = c(1, 1, 0, 0))

MetricsTrueNegative <- function(actual, predicted){
  idx_negatives <- which(actual == 0)
  return(sum(1-predicted[idx_negatives]))
}
MetricsTrueNegative(actual = c(1, 0, 0, 1), predicted = c(1, 0, 0, 0))

MetricsFalseNegative <- function(actual, predicted){
  idx_positives <- which(actual == 1)
  return(sum(1-predicted[idx_positives]))
}
MetricsFalseNegative(actual = c(1, 1, 1, 1), predicted = c(1, 0, 0, 0))

MetricsSensitivity <- function(actual, predicted){
  p <- MetricsPositive(actual, predicted)
  tp <- MetricsTruePositive(actual, predicted)
  return(tp/p)
}
MetricsSensitivity(actual = c(1, 1, 1, 1), predicted = c(1, 0, 0, 1))

MetricsSpecificity <- function(actual, predicted){
  n <- MetricsNegative(actual, predicted)
  tn <- MetricsTrueNegative(actual, predicted)
  return(tn/n)
}
MetricsSpecificity(actual = c(1, 1, 0, 0), predicted = c(0, 0, 0, 1))

MetricsMissRate <- function(actual, predicted){
  p <- MetricsPositive(actual, predicted)
  fn <- MetricsFalseNegative(actual, predicted)
  return(fn/p)
}
MetricsMissRate(actual = c(1, 1, 0, 0), predicted = c(0, 1, 0, 1))

