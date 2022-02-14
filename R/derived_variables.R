
# PUBLIC FUNCTIONS --------------------------------------------------------


# Low level utilities -----------------------------------------------------

#' Append a rowwise summary indicator column
#'
#' Low level function that appends a (binary) indicator column to a dataframe,
#' summarising rowwise whether any values in a selection of existing columns
#' match one or more specified values.
#'
#' @param df Data frame
#' @param cols Character vector of columns to summarise
#' @param yes_values Vector. \code{new_indicator_colname} will be "Yes" if any
#'   values in \code{cols} (rowwise) match one of these values.
#' @param new_indicator_colname Name of new indicator column.
#' @param new_indicator_var_label Character string to set as variable label for
#'   new indicator column.
#' @param new_indicator_yes_value Default is "Yes".
#' @param new_indicator_no_value Default is "No"
#'
#' @return A data frame with an additional column, named as per
#'   \code{new_indicator_colname}.
#'
#' @importFrom rlang :=
#' @export
#' @examples
#' # example input df
#' df <- tibble::tribble(
#' ~med_1, ~med_2, ~med_3, ~IGNORED_COL,
#' "Beta blocker", NA, NA, "Beta blocker",
#' "ACE inhibitor", "Beta blocker", NA, "Beta blocker",
#' "Aspirin", "Aspirin", NA, "Beta blocker",
#' "Aspirin", "ACE inhibitor", NA, "Beta blocker"
#' )
#'
#' # add indicator col
#' df <- append_rowwise_summary_indicator(df = df,
#' cols = c("med_1", "med_2", "med_3"),
#' yes_values = c("Beta blocker", "ACE inhibitor"),
#' new_indicator_colname = "on_anti_hypertensive",
#' new_indicator_var_label = "Taking anti-hypertensive medication",
#' new_indicator_yes_value = "Yes",
#' new_indicator_no_value = "No")
#'
#' # print result
#' print(df)
#'
#' # print variable label for new column
#' print(attributes(df$on_anti_hypertensive)$label)
append_rowwise_summary_indicator <- function(df,
                                                 cols,
                                                 yes_values,
                                                 new_indicator_colname,
                                                 new_indicator_var_label,
                                                 new_indicator_yes_value = "Yes",
                                                 new_indicator_no_value = "No") {
  # mutate new indicator column = "Yes" if any of the
  df <- df %>%
    dplyr::mutate(
      !!new_indicator_colname := ifelse(
        rowSums(
          dplyr::across(tidyselect::all_of(cols),
                        ~ .x %in% yes_values),
          na.rm = TRUE
        ) > 0,
        yes = new_indicator_yes_value,
        no = new_indicator_no_value
      )
    )

  # add variable label
  attributes(df[[new_indicator_colname]])$label <-
    new_indicator_var_label

  return(df)
}


# High level utilities ----------------------------------------------------

#' Summarise UK Biobank touchscreen responses for self-reported medications
#' (cholesterol, blood pressure, diabetes or take exogenous hormones)
#'
#' Summarises field IDs
#' [6153](https://biobank.ndph.ox.ac.uk/ukb/field.cgi?id=6153) and
#' [6177](https://biobank.ndph.ox.ac.uk/ukb/field.cgi?id=6177). The returned
#' data frame contains separate indicator columns for each medication type.
#'
#' Note that [field ID
#' 6153](https://biobank.ndph.ox.ac.uk/ukb/field.cgi?id=6153) applies only to
#' female participants, where as [field
#' ID6177](https://biobank.ndph.ox.ac.uk/ukb/field.cgi?id=6177) applies to male
#' participants.
#'
#' @inheritParams ukbwranglr::read_ukb
#' @param instance Character. Optionally specify a single instance to be
#'   summarised e.g. '1'. Default is `NULL`, in which case the appended columns
#'   will reflect a summary across all available instances.
#' @param .drop If `TRUE`, remove input columns for Field IDs 6153 and 6177 from
#'   output.
#'
#' @return A data frame.
#' @export
derive_touchscreen_self_reported_med_bp_dm_chol_exog_hormones <-
  function(ukb_main,
           ukb_data_dict = ukbwranglr::get_ukb_data_dict(),
           instance = NULL,
           .drop = TRUE) {
    # validate args --------

    ## instance must be between 1 and 4
    if (!is.null(instance)) {
      match.arg(instance,
                choices = as.character(1:4))
    }

    # make data dict
    data_dict <- ukbwranglr::make_data_dict(ukb_main,
                                            ukb_data_dict = ukb_data_dict)

    # required cols for FIDs 6153 and 6177 --------
    required_cols <- ukb_data_dict %>%
      dplyr::filter(.data[["FieldID"]] %in% c("6153", "6177"))

    if (!is.null(instance)) {
      required_cols <- required_cols %>%
        dplyr::filter(.data[["instance"]] == !!instance)
    }

    required_cols <- required_cols$descriptive_colnames

    # check required cols are present
    missing_cols <-
      subset(required_cols, !required_cols %in% names(ukb_main))

    assertthat::assert_that(
      length(missing_cols) == 0,
      msg = paste0(
        "Error! The following required cols are missing from `ukb_main`: ",
        stringr::str_c(missing_cols,
                       sep = "",
                       collapse = ", ")
      )
    )

    # df of all parameters ---------
    params <- tibble::tribble(
      ~ yes_values,
      ~ new_indicator_colname,
      ~ new_indicator_var_label,
      "Cholesterol lowering medication",
      "ts_self_reported_med_chol",
      "Self-reported cholesterol-lowering medication (touchscreen, baseline)",
      "Blood pressure medication",
      "ts_self_reported_med_bp",
      "Self-reported anti-hypertensive medication (touchscreen, baseline)",
      "Insulin",
      "ts_self_reported_med_insulin",
      "Self-reported insulin (touchscreen, baseline)",
      "Hormone replacement therapy",
      "ts_self_reported_med_hrt",
      "Self-reported HRT (touchscreen, baseline)",
      "Oral contraceptive pill or minipill",
      "ts_self_reported_med_ocp_minipill",
      "Self-reported OCP/minipill (touchscreen, baseline)",
      "None of the above",
      "ts_self_reported_med_none_of_the_above",
      "Self-reported not taking cholesterol/BP/insulin/HRT/OCP/minipill medication (touchscreen, baseline)",
      "Do not know",
      "ts_self_reported_med_dont_know",
      "Self-reported don't know if taking cholesterol/BP/insulin/HRT/OCP/minipill medication (touchscreen, baseline)",
      "Prefer not to answer",
      "ts_self_reported_med_prefer_not_to_answer",
      "Self-reported prefer not to answer whether taking cholesterol/BP/insulin/HRT/OCP/minipill medication (touchscreen, baseline)"
    )

    if (!is.null(instance)) {
      params$new_indicator_colname <- paste0(params$new_indicator_colname,
                                             "_",
                                             instance)
    }

    params$new_indicator_colname <- paste0(params$new_indicator_colname,
                                           "_indicator")

    # append summary indicator cols ---------
    for (i in 1:nrow(params)) {
      # parameters
      arg_params <- as.list(params[i, ])

      # add indicator col
      ukb_main <- append_rowwise_summary_indicator(
        df = ukb_main,
        cols = required_cols,
        yes_values = arg_params$yes_values,
        new_indicator_colname = arg_params$new_indicator_colname,
        new_indicator_var_label = arg_params$new_indicator_var_label,
        new_indicator_yes_value = "Yes",
        new_indicator_no_value = "No"
      )
    }

    # drop input cols ----------
    if (.drop) {
      ukb_main <- ukb_main %>%
        dplyr::select(-tidyselect::all_of(required_cols))
    }

    # return result --------
    return(ukb_main)
  }

# PRIVATE FUNCTIONS -------------------------------------------------------


