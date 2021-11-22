
# PUBLIC FUNCTIONS --------------------------------------------------------

#' Derive indicator columns for touchscreen self-reported cholesterol, blood
#' pressure, insulin or exogenous hormones medications
#'
#' The data for Field IDs
#' \href{6177}{https://biobank.ndph.ox.ac.uk/ukb/field.cgi?id=6177} (men) and
#' \href{https://biobank.ndph.ox.ac.uk/ukb/field.cgi?id=6153}{6153} (women,
#' additionally includes values for HRT and OCP/minipill) are in wide format.
#' This function separates out the data in arrays 0-3 for instance 0 \emph{only}
#' into separate indicator 'Yes'/'No' columns, one for each medication type.
#'
#' @param ukb_main Main UKB dataset (data frame)
#' @param data_dict A data dictionary for \code{ukb_main}
#' @param .drop If \code{TRUE}, remove original columns from result. Default is
#'   \code{TRUE}.
#'
#' @return A data frame with additional indicator columns for cholesterol, blood
#'   pressure, insulin, HRT an OCP/minipill medications.
#' @export
serive_touchscreen_self_reported_med_bp_dm_chol_exog_hormones_instance_0 <-
  function(ukb_main,
           data_dict = ukbwranglr::get_ukb_data_dict(),
           .drop = TRUE) {

    REQUIRED_FIDS <- c("6177", "6153")

    # validate args
    assertthat::assert_that(
      all(REQUIRED_FIDS %in% data_dict$FieldID),
      msg = paste0(
        "Required (some or both) Field IDs not present in `data_dict`. Required FieldIDs: ",
        stringr::str_c(REQUIRED_FIDS,
                       sep = "",
                       collapse = ", ")
      )
    )

    ukbwranglr:::check_required_cols_exist(df = ukb_main,
                                           REQUIRED_FIDS)

    # get required fields for instance 0 only (men *and* women)
    touchscreen_self_reported_medication_bp_dm_chol_exog_hormones_cols <- data_dict %>%
      dplyr::filter(.data[["FieldID"]] %in% REQUIRED_FIDS) %>%
      dplyr::filter(.data[["instance"]] == "0") %>%
      dplyr::pull(.data[["descriptive_colnames"]])

    # reformat
    reformatted_ts_sr_meds <- ukb_main %>%
      # select only required cols
      dplyr::select(eid,
                    tidyselect::all_of(.data[["touchscreen_self_reported_medication_bp_dm_chol_exog_hormones_cols"]])) %>%

      # to long format
      tidyr::pivot_longer(tidyselect::all_of(.data[["touchscreen_self_reported_medication_bp_dm_chol_exog_hormones_cols"]])) %>%

      # remove 'special' vvalue responses
      dplyr::filter(!.data[["value"]] %in% c("None of the above",
                                  "Do not know",
                                  "Prefer not to answer")) %>%

      # nest
      dplyr::group_by(.data[["eid"]]) %>%
      tidyr::nest() %>%

      # mutate indicator cols, one for each medication type
      dplyr::mutate(
        touchscreen_self_reported_medication_for_cholesterol_instance_0 = purrr::map_chr(
          .data[["data"]],
          ~ ifelse(sum(.x$value == "Cholesterol lowering medication",
                       na.rm = TRUE),
                   yes = "Yes",
                   no = "No")
        ),
        touchscreen_self_reported_medication_for_bp_instance_0 = purrr::map_chr(
          .data[["data"]],
          ~ ifelse(sum(.x$value == "Blood pressure medication",
                       na.rm = TRUE),
                   yes = "Yes",
                   no = "No")
        ),
        touchscreen_self_reported_insulin_instance_0 = purrr::map_chr(
          .data[["data"]],
          ~ ifelse(sum(.x$value == "Insulin",
                       na.rm = TRUE),
                   yes = "Yes",
                   no = "No")
        ),
        touchscreen_self_reported_hrt_instance_0 = purrr::map_chr(
          .data[["data"]],
          ~ ifelse(sum(.x$value == "Hormone replacement therapy",
                       na.rm = TRUE),
                   yes = "Yes",
                   no = "No")
        ),
        touchscreen_self_reported_ocp_or_minipill_instance_0 = purrr::map_chr(
          .data[["data"]],
          ~ ifelse(sum(.x$value == "Oral contraceptive pill or minipill",
                       na.rm = TRUE),
                   yes = "Yes",
                   no = "No")
        )
      ) %>%

      # remove nested data col and unnest
      dplyr::select(-.data[["data"]]) %>%
      tidyr::unnest(cols = c()) %>%
      dplyr::ungroup()

    # remove original cols
    if (.drop) {
      ukb_main <- ukb_main %>%
        dplyr::select(
          -tidyselect::all_of(
            .data[["touchscreen_self_reported_medication_bp_dm_chol_exog_hormones_cols"]]
          )
        )
    }

    # join reformatted cols back
    ukb_main <- ukb_main %>%
      dplyr::full_join(reformatted_ts_sr_meds,
                       by = "eid")

    # return result
    return(ukb_main)
  }


# PRIVATE FUNCTIONS -------------------------------------------------------


