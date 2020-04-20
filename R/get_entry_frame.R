#' Estimation of all significance levels for entry to binary logistic regression model using current model and external input variables
#'
#' Takes in data frame containing both currently used input variables and external input variables also taking in current binary logistic regression model and names of external input variables; returns result data frame with external names, chi-squared values and corresponding p-values for entry.
#' @param complete_data_frame A data frame containing input variables from the current model and all to be tested external numeric input variables (without NA values) which exist in entry name vector
#' @param current_model A current binary logistic regression model
#' @param entry_name_vector Entry name vector containing names of external numeric input variables (without NA values) to be tested
#' @import stats
#' @import magrittr
#' @importFrom Matrix Matrix Diagonal rankMatrix
#' @return The data frame containing test results for entry of external input variables which includes external names, chi-squared values and p-values for entry
#' @examples
#' # Example 1:
#' library(mlbench)
#' data(PimaIndiansDiabetes2)
#' df <- na.omit(PimaIndiansDiabetes2)
#' df$diabetes <- ifelse(df$diabetes == "pos", 1, 0)
#' m <- glm(diabetes ~ glucose + age + mass + pedigree, df, family = binomial())
#' library(StepwiseBinaryLogisticRegression)
#' print(get_entry_frame(df, m, c("pregnant", "pressure", "triceps", "insulin")))
#'
#' # Example 2:
#' library(StepwiseBinaryLogisticRegression)
#' data(Remission)
#' m2 <- glm(remiss ~ li + temp + cell, Remission, family = binomial())
#' print(get_entry_frame(Remission, m2, c("smear", "infil", "blast")))
#' @export
get_entry_frame <- function(complete_data_frame, current_model, entry_name_vector)
{
  if(family(current_model)$family != "binomial" | family(current_model)$link != "logit")
    stop("The GLM should be binomial logit model.")
  if(length(setdiff(entry_name_vector, names(complete_data_frame))) > 0)
    stop("Some entry names do not exist in complete data frame.")
  entry_frame <- data.frame(var_ = character(), scoreChiSq = numeric(), p.value.entry = numeric(), stringsAsFactors = FALSE)
  for(var_ in entry_name_vector)
  {
    # response variable (which is to be removed from X below) and all current input variables:
    X <- model.frame(current_model)
    # only response variable:
    y <- X[, all.vars(formula(current_model))[1]]
    # predictions with current model:
    p <- predict(current_model, X, type = "response")
    # trying to add the new variable to X matrix:
    X[, var_] <- complete_data_frame[, var_]
    # the unit constant formation for the intercept coefficient by means of replacing/removing response variable:
    X[, 1] <- 1; names(X) <- c("intercept", names(X)[-1]); X <- as.matrix(X)
    # gradient vector with the predictions with current model:
    g <- t(X) %*% (y - p)
    # Hessian matrix with the predictions with current model:
    H <- (Matrix(t(X)) %*% Diagonal(x = -p * (1 - p)) %*% Matrix(X)) %>% as.matrix
    # (-1) * Hessian matrix:
    I <- -H
    # quadratic form for inverse of I matrix:
    if(rankMatrix(I) == min(dim(I)))
      scoreChiSq <- t(g) %*% solve(I) %*% g
    else
      scoreChiSq <- 0
    p.value.entry <- 1 - pchisq(as.numeric(scoreChiSq), df = 1)
    entry_frame %<>% rbind(data.frame(var_ = var_, scoreChiSq = scoreChiSq, p.value.entry = p.value.entry, stringsAsFactors = FALSE))
  }
  return(entry_frame)
}
