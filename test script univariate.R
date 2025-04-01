library(devtools)
load_all()
# MM <- "f1 =~ v1 + v2 + v3 + v4
#            f2 =~ v5 + v6 + v7 + v8"
#
# output_step1 <- step1(exampledata, measurementmodel = MM)
#
# step2output <- step2(output_step1)
#
# id = "id"
#
# A <- create_A(step2output, random_intercept = FALSE,
#               startvalues = matrix(.2, nrow = 2, ncol = 2))
# Q <- create_Q(step2output, random_intercept = FALSE,
#               startvalues = matrix(.2, nrow = 2, ncol = 2))
#
# test <- step3(step2output, id = "id", A = A, Q = Q)
# summary(test$model)
#

cond <- expand.grid(step1def = c("group", "nogroup"),
                    randomintercept = c(TRUE, FALSE),
                    step3def = c("group", "nogroup"))


do_test <- function(pos, cond){
  MMdef <- cond$MMdef[pos] |> as.character()
  step1def <- cond$step1def[pos] |> as.character()
  randomintercept <- cond$randomintercept[pos] |> as.logical()
  step3def <- cond$step3def[pos] |> as.character()

  MM <- "f1 =~ v1 + v2 + v3 + v4"

  if(step1def == "nogroup"){
    step1output <- step1(exampledata, MM)
  } else {
    step1output <- step1(exampledata, MM, group = "group",
                         invariances = c("loadings", "intercepts"))
  }

  step2output <- step2(step1output)

  A <- create_A(step2output, random_intercept = randomintercept,
                startvalues = matrix(.2, nrow = 1, ncol = 1))
  Q <- create_Q(step2output, random_intercept = randomintercept,
                startvalues = matrix(.2, nrow = 1, ncol = 1))

  if(step3def == "group"){
    step3output <- step3(step2output, id = "id", A = A, Q = Q, step3group = "group")
  } else {
    step3output <- step3(step2output, id = "id", A = A, Q = Q)
  }

  print(paste0("Row ", pos, " finished."))
}

test <- lapply(1:nrow(cond), do_test, cond = cond)
