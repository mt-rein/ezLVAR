#' Synthetic Example Data
#'
#' @description A synthetic data set with intensive longitudinal data of 60
#'   individuals, with 56 observations per individual. Every individual belongs
#'   to one of two groups which differ in their dynamic processes. At each
#'   time-point, two constructs were assessed with four items each. The items
#'   `v1` to `v4` measure construct 1, and `v5` to `v8` measure construct 2.
#'
#'
#' @format ## `exampledata` A data frame with 3,360 rows and 11 columns:
#' \describe{
#'   \item{v1, v2, v3, v4}{Items measuring construct 1}
#'   \item{v5, v6, v7, v8}{Items measuring construct 2}
#'   \item{id}{ID variable}
#'   \item{obs}{Observation counter per individual}
#'   \item{group}{Grouping variable}
#' }
"exampledata"
