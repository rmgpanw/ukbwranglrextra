#' Create a data frame with cases and matched controls
#'
#'
#'
#' @param df A data frame
#' @param case_control_colname Character. The name of an indicator column in
#'   \code{df} that labels individuals as either cases or controls.
#' @param case_value String Used to determine which individuals are cases. This
#'   value should be present in \code{df[[case_control_colname]]}.
#' @param control_value String Used to determine which individuals are controls
#'   This value should be present in \code{df[[case_control_colname]]}.
#' @param matching_vars_categorical Character vector of column names containing
#'   categorical variables to match on. These columns should all be of type
#'   character.
#' @param matching_vars_continuous Named integer vector. Names should be column
#'   names containing continuous variables to match on (these columns should all
#'   be of type numeric). Integer values are the number of digits to match to.
#'   For example, `c("age_yrs" = 0)` specifies that controls will be matched on
#'   'age_yrs' to the nearest year.
#' @param n_control_ratio Integer. The number of controls to match for each
#'   case.
#' @param seed Integer, for reproducibility.
#'
#' @return A data frame
#' @export
create_case_control_matched_df <- function(df,
                                           case_control_colname,
                                           case_value = "Case",
                                           control_value = "Control",
                                           matching_vars_categorical = NULL,
                                           matching_vars_continuous = NULL,
                                           n_control_ratio = 4,
                                           seed = 42)
{
  # Validate arguments -------

  # check case_control_colname exists in df
  assertthat::assert_that(case_control_colname %in% names(df),
                          msg = paste0("Error! ",
                                       case_control_colname,
                                       " not present in df"))

  # check that expected case/control values are present
  assertthat::assert_that(all(c(case_value, control_value) %in% df[[case_control_colname]]),
                          msg = paste0("Error! Values ",
                                       case_value,
                                       " and/or ",
                                       control_value,
                                       " are not present in df[[",
                                       case_control_colname,
                                       "]]"))

  # check that at least one of categorical or continuous matching vars is supplied
  assertthat::assert_that(!is.null(matching_vars_categorical) | !is.null(matching_vars_continuous),
                          msg = "At least one of either `matching_vars_continuous` or `matching_vars_categorical` must not be NULL")

  # check that matching variables exist in `df`
  if (!is.null(matching_vars_continuous)) {
    assertthat::assert_that(all(names(matching_vars_continuous) %in% names(df)),
                          msg = "Error! Some `matching_vars_continuous` values are not present in `df`")
  }

  if (!is.null(matching_vars_categorical)) {
  assertthat::assert_that(all(matching_vars_categorical %in% names(df)),
                          msg = "Error! Some `matching_vars_categorical` values are not present in `df`")
  }

  # continuous matching vars - create df of arguments to be used later
  if (!is.null(matching_vars_continuous)) {
    matching_vars_continuous <- matching_vars_continuous %>%
      as.list() %>%
      as.data.frame() %>%
      tidyr::pivot_longer(tidyselect::everything(),
                          names_to = "matching_vars_continuous",
                          values_to = "round_n_digits") %>%
      dplyr::mutate("matching_vars_continuous_rounded" = paste(.data[["matching_vars_continuous"]],
                                                               "rounded",
                                                               sep = "..."))

    # TODO - add check no columns in data with these names

    # convert df to list of dfs: one df per matching variable, each with a single
    # row containing a value for `round_n_digits` and a temporary column name
    # (`matching_vars_continuous_rounded`)
    matching_vars_continuous_split <- split(matching_vars_continuous,
                                            matching_vars_continuous$matching_vars_continuous)

    # append cols - round continuous variables to specified n digits
    for (continuous_var in matching_vars_continuous$matching_vars_continuous) {
      # get parameters
      rounded_colname <-
        matching_vars_continuous_split[[continuous_var]]$matching_vars_continuous_rounded
      round_n_digits <-
        matching_vars_continuous_split[[continuous_var]]$round_n_digits

      # append rounded version of continuous variable
      df[[rounded_colname]] <- round(df[[continuous_var]],
                                     round_n_digits)
    }

    # Now generate 'matching ids' from the matching variables

    # filter for cases and select matching variables - if
    # `matching_vars_categorical` is `NULL`, then `all_matching_vars` will
    # simply be `matching_vars_continuous$matching_vars_continuous_rounded`
    all_matching_vars <- c(
      matching_vars_categorical,
      matching_vars_continuous$matching_vars_continuous_rounded
    )
  } else {
    all_matching_vars <- matching_vars_categorical
  }

  matching_id_df <- df[df[[case_control_colname]] == case_value,]

  matching_id_df <- matching_id_df %>%
    dplyr::select(tidyselect::all_of(all_matching_vars))

  # TODO - add step to check there are any missing values for matching variables

  matching_id_df <- matching_id_df %>%
    dplyr::arrange(dplyr::across(tidyselect::all_of(all_matching_vars))) %>%
    dplyr::distinct()

  # matching ids - one per unique combination of matching variables
  matching_id_df[["matching_id"]] <- 1:nrow(matching_id_df)

  # create cases df - append matching id
  cases_df <- df[df[[case_control_colname]] == case_value, ]

  cases_df <- cases_df %>%
    dplyr::left_join(matching_id_df,
                     by = all_matching_vars) %>%
    dplyr::select(-tidyselect::all_of(
      matching_vars_continuous$matching_vars_continuous_rounded
    ))

  # create controls df - append matching id and filter for these
  controls_df <- df[df[[case_control_colname]] == control_value, ]

  controls_df <- controls_df %>%
    dplyr::left_join(matching_id_df,
                     by = all_matching_vars) %>%
    dplyr::select(-tidyselect::all_of(
      matching_vars_continuous$matching_vars_continuous_rounded
    )) %>%
    dplyr::filter(.data[["matching_id"]] %in% cases_df[["matching_id"]])

  # check that controls_df contains all the required 'matching ids'
  missing_match_ids <- cases_df %>%
    dplyr::filter(!.data[["matching_id"]] %in% controls_df[["matching_id"]]) %>%
    dplyr::pull(.data[["matching_id"]]) %>%
    unique()

  if (length(missing_match_ids) > 0) {
    stop(paste0("Error! There are ",
                length(unique(cases_df[["matching_id"]])),
                " unique combinations of the selected matching variables in the cases group, but ",
                length(missing_match_ids),
                " of these are not present in the control group. Try reducing the number of matching variables"))
  }

  ########
  # n cases per matching id
  n_cases_per_matching_id <- cases_df %>%
    dplyr::count(.data[["matching_id"]])

  n_controls_per_matching_id <- n_cases_per_matching_id

  n_controls_per_matching_id$n_controls_to_sample <-
    n_cases_per_matching_id$n * n_control_ratio

  n_controls_per_matching_id$n <- NULL

  # get controls
  controls_df_split <- split(controls_df,
                             controls_df$matching_id)

  # TODO check if there are any matching_ids where there are not enough controls
  # (less than `n_controls_to_sample`)

  control_sampling_mapper <- function(matching_id,
                                      n_controls_to_sample) {
    controls_df_split[[matching_id]] %>%
      dplyr::slice_sample(n = n_controls_to_sample,
                          replace = FALSE)
  }

  set.seed(seed)
  controls_df <- n_controls_per_matching_id %>%
    # for each 'subset', randomly select 4 controls
    purrr::pmap(control_sampling_mapper) %>%
    dplyr::bind_rows()

  # final result - recombine with cases
  result <- dplyr::bind_rows(cases_df,
                             controls_df)

  # add back labels
  for (col_name in names(result)) {
    col_label <- attributes(df[[col_name]])$label

    if (!is.null(col_label)) {
      attributes(result[[col_name]])$label <- col_label
    }
  }

  return(result)
}
