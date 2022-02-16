#'Create multiple models which vary by one parameter
#'
#'@description Use when creating multiple models which vary by one parameter.
#'  See example below.
#'
#'@details Ideally aim to supply a `df` without missing data. Otherwise, check N
#'  observations for models (in `fit` column, `model_glance`) to see if many (or
#'  any) observations were removed.
#'
#'@param df Data frame.
#'@param model_str A model formula as a string, containing one variable as `.x`.
#'@param engine A [`parsnip`](https://www.tidymodels.org/find/parsnip/) engine.
#'@param params A character vector of values to substitute `.x` in `model_str`.
#'@param rm_raw_model If `TRUE` (the default), the raw models are removed.
#'
#'@return A data frame
#'@export
#'@examples
#' # dummy dataset
#' mtcars2 <- mtcars %>%
#'  tibble::rownames_to_column(var = "make") %>%
#'  dplyr::mutate("is_merc" = ifelse(stringr::str_detect(.data[["make"]], pattern = "^Merc"),
#'                                   yes = 1,
#'                                   no = 0 )) %>%
#'  tibble::as_tibble() %>%
#'  dplyr::mutate(dplyr::across(tidyselect::starts_with("is_"), as.factor))
#'
#' # preview dummy dataset
#' mtcars2
#'
#' # run `map_models()` - logistic regression to predict Mercedes make
#' result <-
#' map_models(
#'   df = mtcars2,
#'   model_str = "is_merc ~ mpg + cyl + disp + {.x}",
#'   params = c(GEAR = 'gear', CARB = 'carb'),
#'   engine = parsnip::set_engine(object = parsnip::logistic_reg(), engine = "glm"),
#'   rm_raw_model = FALSE
#' )
#'
#' # view result
#' result
#'
#' # model outputs are stored in the `fit` column
#' names(result$fit[[1]])
#'
#' # 'tidy' model outputs are under `model_tidy`
#' result$fit[[1]]$model_tidy
#'
#' # 'glance' model outputs are under `model_glance`
#' result$fit[[1]]$model_glance
#'
#' # `model_raw` contains either the raw model under `result`, or an error message under `error`
#' names(result$fit[[1]]$model_raw)
map_models <- function(df,
                       params,
                       model_str,
                       engine =  parsnip::set_engine(object = parsnip::logistic_reg(),
                                             engine = "glm"),
                       rm_raw_model = TRUE) {
  # create list of model formulas
  formulas <- params %>%
    purrr::set_names() %>%
    purrr::map(~ as.formula(rlang::parse_expr(glue::glue(model_str))))

  # combine with parsnip engine to create workflowset
  models <- workflowsets::workflow_set(preproc = formulas,
                         models = list(model = engine)) %>%
    dplyr::mutate("wflow_id" = stringr::str_remove(.data[["wflow_id"]],
                                 "_model$"))

  # append labels
  if (!is.null(names(params))) {
    # create dictionary
    dict <- params %>%
      as.list() %>%
      as.data.frame() %>%
      tidyr::pivot_longer(tidyselect::everything())

    dict <- as.list(stats::setNames(object = dict$name,
                            nm = dict$value))

    # add 'label' col
    models$wflow_label <- dplyr::recode(models$wflow_id,
                                        !!!dict)
  }

  # fit
  total_n <- nrow(models)

  models <- models %>%
    dplyr::mutate("fit" = purrr::imap(.data[["info"]], ~ {
      message(glue::glue("{wflow_id[[.y]]}, {.y} of {total_n}"))

      # set up
      result <- list()
      fit_safely <- purrr::safely(purrr::partial(parsnip::fit, .x$workflow[[1]], df), otherwise = NULL)

      # fit
      result$model_raw <- fit_safely()

      # if error, show as warning
      if (!is.null(result$model_raw$error)) {
        warning("Error! ",result$model_raw$error$message)
      }

      # if successful, tidy and glance with broom. Otherwise, return `NULL`
      if (!is.null(result$model_raw)) {
        result$model_tidy <- broom::tidy(result$model_raw$result)

        result$model_glance <- broom::glance(result$model_raw$result)

        # remove raw model if requested
        if (rm_raw_model) {
          result$model_raw$result <- NA
        }
      } else {
        result <- NULL
      }

      return(result)
    }))

  # add metadata
  models$model_formula <- model_str
  models$model_engine <- engine %>%
    purrr::keep(is.character) %>%
    purrr::reduce(paste, sep = "; ")

  # return results
  return(models)
}
