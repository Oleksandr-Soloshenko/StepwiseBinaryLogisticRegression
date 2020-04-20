#' Backward variable selection using given full model
#'
#' Takes in any full 'glm' model and significance level for stay (sls) then conducts backward variable selection and returns shrunk model.
#' @param model A full 'glm' model preliminarily built with formula/data/family parameters only
#' @param sls A significance level for stay (sls) applied for each backward elimination step
#' @param printing A logical parameter indicating report information to console output
#' @import stats
#' @return The shrunk 'glm' model where all input variables satisfy significance level for stay
#' @examples
#' # Example 1:
#' library(mlbench)
#' data(PimaIndiansDiabetes2)
#' df <- na.omit(PimaIndiansDiabetes2) # removing rows with any NA values
#' df$diabetes <- ifelse(df$diabetes == "pos", 1, 0)
#' full_model <- glm(diabetes ~ ., data = df, family = binomial())
#' print(formula(full_model))
#' library(StepwiseBinaryLogisticRegression)
#' model <- backward.glm(full_model)
#' print(formula(model))
#' print(summary(model))
#'
#' # Example 2:
#' library(StepwiseBinaryLogisticRegression)
#' data(Remission)
#' full_model2 <- glm(remiss ~ ., data = Remission, family = binomial())
#' print(formula(full_model2))
#' model2 <- backward.glm(full_model2, sls = 0.35)
#' print(formula(model2))
#' print(summary(model2))
#' @export
backward.glm <- function(model, sls = 0.05, printing = TRUE)
{
  while(TRUE)
  {
    p.values <- coef(summary(model))[, "Pr(>|z|)"][-1];
    p.values_bad <- p.values[p.values > sls];
    if(length(p.values_bad) == 0) return(model)
    else {
      # removing the worst variable & re-estimation the model:
      the_worst_variable <- names(p.values_bad[which.max(p.values_bad)]);
      model <- glm(formula = update(formula(model), paste("~ . -", the_worst_variable)), data = model.frame(model), family = family(model));
      if(printing) print(paste0(the_worst_variable, " is REMOVED."))
    }
  }
}
