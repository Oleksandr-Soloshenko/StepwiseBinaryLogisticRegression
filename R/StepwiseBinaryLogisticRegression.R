#' Creation of binary logistic regression model using stepwise variable selection with significance levels for entry (sle) and for stay (sls)
#'
#' Takes in formula with both binary response variable and not empty numeric variables separated with '+' sign or just '.' sign to try all numeric not empty inputs also taking in data frame and significance levels for entry and for stay.
#' @param formula A formula which contains both binary response variable and not empty numeric variables (without NA values) separated with '+' sign or just '.' sign to try all not empty numeric inputs (without NA values). An example: \code{y ~ .} or \code{y ~ x1 + x2 + x3}
#' @param data A data frame containing both binary response variable and not empty (without NA values) numeric input variables (which are mentioned in formula)
#' @param printing A logical parameter indicating report information to console output
#' @param sle A significance level for entry (for each forward entry step)
#' @param sls A significance level for stay (for each backward elimination step)
#' @import stats
#' @import formula.tools
#' @return The binomial/binary 'glm' model with logit link function containing only the input variables included with stepwise variable selection
#' @examples
#' # Example 1:
#' library(mlbench)
#' data(PimaIndiansDiabetes2)
#' df <- na.omit(PimaIndiansDiabetes2) # removing rows with any NA values
#' df$diabetes <- ifelse(df$diabetes == "pos", 1, 0)
#' library(StepwiseBinaryLogisticRegression)
#' m <- StepwiseBinaryLogisticRegression(diabetes ~ pregnant + glucose + pressure +
#' triceps + insulin + mass + pedigree + age, df)
#' print(formula(m))
#' print(summary(m))
#'
#' # Example 2 (actually the same one in this situation):
#' library(mlbench)
#' data(PimaIndiansDiabetes2)
#' df <- na.omit(PimaIndiansDiabetes2) # removing rows with any NA values
#' df$diabetes <- ifelse(df$diabetes == "pos", 1, 0)
#' library(StepwiseBinaryLogisticRegression)
#' m2 <- StepwiseBinaryLogisticRegression(diabetes ~ ., df)
#' print(formula(m2))
#' print(summary(m2))
#'
#' # Example 3:
#' library(StepwiseBinaryLogisticRegression)
#' data(Remission)
#' # here sle <= sls just in case to prevent any possible looping:
#' m3 <- StepwiseBinaryLogisticRegression(remiss ~ ., Remission, sle = 0.30, sls = 0.35)
#' print(formula(m3))
#' print(summary(m3))
#' @export
StepwiseBinaryLogisticRegression <- function(formula, data, sle = 0.05, sls = 0.05, printing = TRUE)
{
  y_name <- lhs.vars(formula) #as.character(formula[[2]])

  if(length(y_name) != 1) # | y_name == ".")
    stop("The left side of formula needs to contain the only one dependent variable.")

  if(!is.numeric(data[, y_name]) | !all(unique(data[, y_name]) %in% 0:1) | !(length(unique(data[, y_name])) == 2))
    stop("The model response needs to be the numeric binary variable with values in {0, 1}. NA values and constant response are not allowed.")

  all_vars <- all.vars(formula)
  if(all_vars[2] == ".") {
    all_vars <- all_vars[1]
    all_vars <- c(all_vars, names(data)[names(data) != all_vars[1]])
  }
  else {
    split_ <- sort(strsplit(gsub(".* ~ ", "", as.character(formula)), " \\+ ")[[1]])
    rhs_vars <- sort(rhs.vars(formula)) #sort(all.vars(formula)[-1])
    err.msg <- "The formula should contain input variables separated with ' + ' sign only; or just use ' . ' to take all columns. An example: y ~ x1 + x2 + x3 + x4 or y ~ ."
    if(length(split_) != length(rhs_vars)) stop(err.msg)
    if(!all(split_ == rhs_vars)) stop(err.msg)
  }
  complete_df  <- data[, all_vars]
  if(sum(is.na(complete_df[, -1])) > 0)
    stop("The input variables need to contain no NA values. Please use 'data <- na.omit(data)' or just don't use columns with NA values (or fill NA values).")
  if(!all(sapply(complete_df[, -1], is.numeric)))
    stop("The input variables need to be numeric only.")
  all_x_names <- all_vars[-1]
  curr_formula <- as.formula(paste(y_name, "~ 1"))

  while(TRUE)
  {
    curr_model <- backward.glm(glm(formula = curr_formula, data = complete_df, family = binomial(link = "logit")), sls = sls, printing = printing);
    entry_names <- setdiff(all_x_names, row.names(coef(summary(curr_model)))[-1]);
    entry_calc <- get_entry_frame(complete_data_frame = complete_df, current_model = curr_model, entry_name_vector = entry_names)
    if(printing) print(entry_calc)
    if(all(entry_calc$p.value.entry > sle)) return(curr_model)
    else {
      variable_to_add <- entry_calc[which.min(entry_calc$p.value.entry), "var_"]
      curr_formula <- as.formula(paste0(y_name, " ~ ", paste0(c(row.names(coef(summary(curr_model)))[-1], variable_to_add), collapse = " + ")))
    }
    if(printing) print(curr_formula)
  }
}
