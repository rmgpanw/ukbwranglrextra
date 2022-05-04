
# PUBLIC ------------------------------------------------------------------

#' Map codes in a UKB clinical events table
#'
#' Joins the existing `code` column to a mapping data frame created by
#' [codemapper::get_mapping_df], then filters out any missing values for the new
#' code type. A column is appended containing the newly mapped codes.
#'
#' The `clinical_events` data frame can contain multiple sources, but an error
#' is raised by default if these do not all use the same data coding (see
#' [ukbwranglr::clinical_events_sources()] for recognised data sources and their
#' respective data codings).
#'
#' @param clinical_events A clinical events data frame created by
#'   [ukbwranglr::tidy_clinical_events()]
#' @param from Coding type being mapped from
#' @param to Coding type being mapped to
#' @param from_colname Name of column containing codes to be mapped from. If
#'   `NULL` (default), this is assumed to be named 'code'.
#' @param to_colname Name of new column containing mapped codes. If `NULL`
#'   (default), this will equal the value for argument `to`.
#' @param all_lkps_maps Named list of SQLite database with lookup/mapping tables
#' @param strict_ukb If `TRUE`, extra checks are performed to check that
#'   `source` and `data_coding` match [ukbwranglr::clinical_events_sources()],
#'   and that column names match those expected for a clinical events table.
#' @noRd
#'
#' @return A clinical events data frame.
map_codes_ukb_clinical_events <- function(clinical_events,
                                          from,
                                          to,
                                          from_colname = NULL,
                                          to_colname = NULL,
                                          all_lkps_maps = "all_lkps_maps.db",
                                          strict_ukb = TRUE) {
  # validate args
  codemapper:::check_mapping_args(from = from,
                                  to = to)

  if (!is.null(to_colname)) {
    assertthat::is.string(to_colname)
  }

  if (!is.null(from_colname)) {
    assertthat::is.string(from_colname)
  }

  if (strict_ukb) {
    ukbwranglr:::validate_clinical_events_and_check_type(clinical_events)

    ## check `clinical_events` contains only one `data_coding` type
    assert_sources_match_data_coding(clinical_events = clinical_events,
                                     data_coding = from)
  }

  # all_lkps_maps
  ## connect to database file path
  if (is.character(all_lkps_maps)) {
    con <- codemapper:::check_all_lkps_maps_path(all_lkps_maps)
    all_lkps_maps <- ukbwranglr::db_tables_to_list(con)
    on.exit(DBI::dbDisconnect(con))
  }

  # create mapping df
  mapping_df <- codemapper::get_mapping_df(from = from,
                                           to = to,
                                           all_lkps_maps = all_lkps_maps)

  # map
  if (is.null(from_colname)) {
    from_colname <- "code"
  }

  join_by <- from
  names(join_by) <- from_colname

  result <- clinical_events %>%
    dplyr::inner_join(mapping_df,
                      by = join_by) %>%
    # remove duplicated rows
    dplyr::distinct()

  # rename newly appended column (optionally)
  if (!is.null(to_colname)) {
    result <- ukbwranglr:::rename_cols(df = result,
                                       old_colnames = to,
                                       new_colnames = to_colname)
  }

  # return result
  return(result)
}

#' Filter UK Biobank clinical events for selected data sources
#'
#' @param clinical_events A clinical events table (data frame or `tbl_dbi`
#'   object) created by [ukbwranglr::tidy_clinical_events()].
#' @param sources A character vector of data sources. Must be listed under
#'   `source` in [ukbwranglr::clinical_events_sources()].
#' @param allow_missing_sources If `FALSE` (default), an error is raised if any
#'   values for `sources` are not present in `clinical_events`. If `TRUE`, a
#'   warning is raised instead.
#'
#' @return A data frame
#' @export
get_clinical_events_source <- function(clinical_events,
                                       sources,
                                       allow_missing_sources = FALSE) {
  # validate args
  assertthat::assert_that(all(sources %in% ukbwranglr::clinical_events_sources()$source))
  ukbwranglr:::validate_clinical_events_and_check_type(clinical_events)

  # check selected sources are present
  check_sources <- clinical_events %>%
    dplyr::filter(.data[["source"]] %in% !!sources) %>%
    dplyr::distinct(.data[["source"]],
                    .keep_all = FALSE) %>%
    dplyr::collect()

  # error/warning if any sources are not present
  missing_sources <-
    subset(sources, !sources %in% check_sources$source)

  if (length(missing_sources) != 0) {
    missing_sources_msg <-
      "The following sources are not present in `clinical events`: "
    if (allow_missing_sources) {
      warning(paste0(
        "Warning! ",
        missing_sources_msg,
        stringr::str_c(missing_sources,
                       sep = "",
                       collapse = ", ")
      ))
    } else {
      stop(paste0(
        "Error!",
        missing_sources_msg,
        stringr::str_c(missing_sources,
                       sep = "",
                       collapse = ", ")
      ))
    }
  }

  # filter clinical events table for selected sources
  result <- clinical_events %>%
    dplyr::filter(.data[["source"]] %in% !!sources) %>%
    dplyr::collect()

  # return result
  return(result)
}

# PRIVATE -----------------------------------------------------------------

#' Asssert that all data sources in a `clinical_events` data frame match a
#' single `data_coding`
#'
#' @param clinical_events Clinical events data frame
#' @param data_coding String
#' @noRd
assert_sources_match_data_coding <- function(clinical_events,
                                             data_coding) {
  # validate args
  ukbwranglr:::validate_clinical_events_and_check_type(clinical_events)
  assertthat::is.string(data_coding)
  assertthat::assert_that(
    data_coding %in% ukbwranglr::clinical_events_sources()$data_coding,
    msg = paste0("Error! Unrecognised `data_coding`: ",
                 data_coding)
  )

  # see what sources are in `clinical_events`
  included_sources <- unique(clinical_events$source)

  # get which data codings are associated with these
  included_data_codings <- ukbwranglr::clinical_events_sources() %>%
    dplyr::filter(.data[["source"]] %in% !!included_sources) %>%
    dplyr::pull(.data[["source"]]) %>%
    unique()

  # check only one `data_coding`, and check that this matches what's expected
  assertthat::are_equal(included_data_codings, data_coding)

  invisible(NULL)
}
