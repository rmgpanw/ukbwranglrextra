#' Generate multiple models
#'
#' WRT missing data, aim to supply a `df` without missing data. Otherwise, check
#' N observations for models (in `fit` column, `model_raw`, `result`,
#' `model_glance`) to see if many/any were removed.
#'
#' @param df Data frame
#' @param dict A data frame with columns 'variable' and 'label'
#' @param model_str A model formula as a string, containing one variable as
#'   `.x`.
#' @param engine A \code{\link{parsnip}} engine
#'
#' @return A data frame
#' @export
map_models <- function(df,
                       params,
                       model_str,
                       engine =  parsnip::set_engine(object = parsnip::logistic_reg(),
                                             engine = "glm"),
                       butcher_model = TRUE) {
  # create list of model formulas
  formulas <- params %>%
    purrr::set_names() %>%
    purrr::map(~ as.formula(rlang::parse_expr(glue::glue(model_str))))

  # combine with parsnip engine to create workflowset
  models <- workflowsets::workflow_set(preproc = formulas,
                         models = list(model = engine)) %>%
    dplyr::mutate(wflow_id = stringr::str_remove(wflow_id,
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
    dplyr::mutate("fit" = purrr::imap(info, ~ {
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
        if (butcher_model) {
          result$model_raw$result <- butcher::butcher(result$model_raw$result)
        }
      } else {
        result <- NULL
      }

      return(result)
    }))

  # return results
  return(models)
}
