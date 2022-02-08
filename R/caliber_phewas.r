library(tidyverse)
library(targets)
library(tidymodels)
library(workflowsets)
library(glue)

tar_load(summarised_dataset_gp_only_NEW)
# tar_load(caliber_comorbidities)
tar_load(caliber_comorbidities_merged)

# full dataset
df <- summarised_dataset_gp_only_NEW %>%
  left_join(caliber_comorbidities_merged,
            by = c("eid",
                   "dob")) %>%
  labelled::unlabelled()

# caliber indicator cols
indicator_cols <- caliber_comorbidities_merged %>%
  select(ends_with("indicator")) %>%
  names()

# list of labels/indicator cols
caliber_indicator_dict <- caliber_comorbidities_merged %>%
  select(ends_with("indicator")) %>%
  labelled::generate_dictionary()

# create list of models
outcome <- "COMBINED_SCAR_CY_indicator"

caliber_indicator_dict <- caliber_indicator_dict %>%
  as_tibble() %>%
  select(wflow_id = variable,
         label)

caliber_formulas <- caliber_indicator_dict %>%
  pull(wflow_id) %>%
  set_names() %>%
  map(~ as.formula(rlang::parse_expr(glue("{outcome} ~ age + sex_f31_0_0 + {.x}"))))

# engine
glm_model <- logistic_reg() %>%
  set_engine("glm")

# create workflowset
caliber_models <- workflow_set(preproc = caliber_formulas,
                               models = list(model = glm_model)) %>%
  mutate(wflow_id = str_remove(wflow_id,
                               "_model$"))

# append labels
caliber_models <- caliber_indicator_dict %>%
  right_join(caliber_models,
             by = "wflow_id")

# fit models
# sample df to increase speed
df_10000 <-
  df %>% sample_n(10000)

total_n <- 10

keep_raw_model <- TRUE

caliber_models <- caliber_models[1:10,] %>%
  mutate(fit = imap(info, ~ {
    message(glue("{label[[.y]]}, {.y} of {total_n}"))

    result <- list()

    result$model_raw <- tryCatch(
      fit(.x$workflow[[1]], df_10000),
      error = function(e)
        NULL
    )

    # extract
    if (!is.null(result$model_raw)) {
      result$model_tidy <- broom::tidy(result$model_raw)


      # remove raw model if requested
      if (!keep_raw_model) {
        result$model_raw <- NULL
      }
    } else {
      result <- NULL
    }

    return(result)
  }))

#' Run many models at once
#'
#' @param df Data frame
#' @param dict A data frame with columns 'variable' and 'label'
#' @param model_str A model formula as a string, containing one variable as \code{.x}.
#' @param engine A \code{\link{parsnip}} engine
#'
#' @return A dataframe
#' @export
map_models <- function(df,
                       dict,
                       model_str,
                       engine = logistic_reg() %>%
                         set_engine("glm")) {
  formulas <- dict %>%
    pull(wflow_id) %>%
    set_names() %>%
    map(~ as.formula(rlang::parse_expr(glue(model_str))))

  models <- workflow_set(preproc = formulas,
                                 models = list(model = engine)) %>%
    mutate(wflow_id = str_remove(wflow_id,
                                 "_model$"))

  # append labels
  models <- dict %>%
    right_join(models,
               by = "wflow_id")

  # fit
  total_n <- nrow(models)

  models <- models[1:10,] %>%
    mutate(fit = imap(info, ~ {
      message(glue("{label[[.y]]}, {.y} of {total_n}"))

      result <- list()

      result$model_raw <- tryCatch(
        fit(.x$workflow[[1]], df_10000),
        error = function(e)
          NULL
      )

      # extract
      if (!is.null(result$model_raw)) {
        result$model_tidy <- broom::tidy(result$model_raw)


        # remove raw model if requested
        if (!keep_raw_model) {
          result$model_raw <- NULL
        }
      } else {
        result <- NULL
      }

      return(result)
    }))

  # return results
  return(models)
}
