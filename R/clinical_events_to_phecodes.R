# EXPORTED ----------------------------------------------------------------

#' Map UK Biobank clinical events to phecodes
#'
#'
#' @inheritParams ukbwranglr::extract_phenotypes
#' @inheritParams codemapper::codes_starting_with
#' @param min_date_only If `TRUE`, result will be filtered for only the earliest
#'   date per eid-phecode pair (date will be recorded as `NA` for cases where
#'   there are no dates).
#'
#' @return A data frame with column names 'eid', 'source', 'index', 'code',
#'   'icd10', 'phecode' and 'date'.
#' @export
map_clinical_events_to_phecodes <- function(clinical_events,
                                            all_lkps_maps = "all_lkps_maps.db",
                                            min_date_only = FALSE) {

  # utility function to strip 'X' from undivided 3 character ICD-10 codes
  strip_x_from_3char_icd10 <- function(df) {
    df %>%
      dplyr::mutate("icd10" = stringr::str_remove(.data[["icd10"]],
                                                  "X$"))
  }

  start_time <- proc.time()
  message("***MAPPING clinical_events TO PHECODES***")

  # Self-reported non-cancer -----------------------------------------------------------
  message("Self-reported non-cancer")
  self_reported_icd10 <-  get_clinical_events_source(
    clinical_events = clinical_events,
    sources = c("f20002_icd10")
  )

  self_reported_icd10 <- self_reported_icd10 %>%
    dplyr::rename("icd10" = .data[["code"]])

  self_reported_icd10 <-
    map_icd10_to_phecode(clinical_events = self_reported_icd10,
                         all_lkps_maps = all_lkps_maps)

  # HES/Cancer ---------------------------------------------------------------------


  ## ICD-10 ------------------------------------------------------------------
  message("HES/cancer register ICD-10")
  hes_cancer_icd10 <-
    get_clinical_events_source(
      clinical_events = clinical_events,
      sources = c("f41270", "f40006")
    )

  hes_cancer_icd10 <- hes_cancer_icd10 %>%
    dplyr::rename("icd10" = .data[["code"]])

  hes_cancer_icd10 <-
    map_icd10_to_phecode(clinical_events = hes_cancer_icd10,
                         all_lkps_maps = all_lkps_maps)

  ## ICD-9 ------------------------------------------------------------------
  message("HES/cancer register ICD-9")
  hes_cancer_icd9 <-
    get_clinical_events_source(
      clinical_events = clinical_events,
      sources = c("f41271", "f40013")
    )

  hes_cancer_icd9 <- map_codes_ukb_clinical_events(
    clinical_events = hes_cancer_icd9,
    from = "icd9",
    to = "icd10",
    all_lkps_maps = all_lkps_maps,
    strict_ukb = FALSE
  ) %>%
    strip_x_from_3char_icd10() %>%
    map_icd10_to_phecode(all_lkps_maps = all_lkps_maps)

  # Death -------------------------------------------------------------------
  message("Death ICD-10")
  death_icd10 <-
    get_clinical_events_source(
      clinical_events = clinical_events,
      sources = c("f40001", "f40002")
    )

  death_icd10 <- death_icd10 %>%
    dplyr::rename("icd10" = .data[["code"]])

  death_icd10 <- map_icd10_to_phecode(clinical_events = death_icd10,
                                      all_lkps_maps = all_lkps_maps)

  # GP - read 3 -------------------------------------------------------------
  message("GP Read 3")
  # note: read 3 codes only for GP data provider 3 (England TPP)

  gp_read3 <- get_clinical_events_source(
    clinical_events = clinical_events,
    sources = c("gpc1_r3", "gpc2_r3", "gpc3_r3", "gpc4_r3"),
    allow_missing_sources = TRUE
  )

  gp_read3 <- gp_read3 %>%
    map_codes_ukb_clinical_events(
      from = "read3",
      to = "icd10",
      all_lkps_maps = all_lkps_maps,
      strict_ukb = FALSE
    ) %>%
    strip_x_from_3char_icd10() %>%
    map_icd10_to_phecode(all_lkps_maps = all_lkps_maps)

  # GP - read 2 -------------------------------------------------------------

  message("GP Read 2")
  # note: no read 2 codes for GP data provider 3 (England TPP)
  gp_read2 <- get_clinical_events_source(
    clinical_events = clinical_events,
    sources = c("gpc1_r2", "gpc2_r2", "gpc3_r2", "gpc4_r2"),
    allow_missing_sources = TRUE
  )

  gp_read2 <- gp_read2 %>%
    map_codes_ukb_clinical_events(
      from = "read2",
      to = "icd10",
      all_lkps_maps = all_lkps_maps,
      strict_ukb = FALSE
    ) %>%
    strip_x_from_3char_icd10() %>%
    map_icd10_to_phecode(all_lkps_maps = all_lkps_maps)

  # Combine -----------------------------------------------------------------

  message("Combining results")
  result <- list(
    self_reported_icd10,
    hes_cancer_icd10,
    hes_cancer_icd9,
    death_icd10,
    gp_read2,
    gp_read3
  ) %>%
    dplyr::bind_rows() %>%
    dplyr::distinct()

  if (min_date_only) {
    result <- result %>%
    # take only earliest date (or `NA` if no dates recorded) per eid-phecode
    dplyr::group_by(.data[["eid"]],
                    .data[["phecode"]]) %>%
    dplyr::arrange(.data[["date"]]) %>%
    dplyr::slice_head(n = 1) %>%
    dplyr::ungroup()
  }

  ukbwranglr:::time_taken_message(start_time)
  return(result)
}

#' Map codes in a UKB clinical events table
#'
#' Joins the existing `code` column to a mapping data frame created by
#' [codemapper::get_mapping_df], then filters out any missing values for the new
#' code type. A column is appended containing the newly mapped codes.
#'
#' The `clinical_events` data frame can contain multiple sources, but an error
#' is raised by default if these do not all use the same data coding (see
#' [ukbwranglr::clinical_events_sources()] for recognised data sources and their
#' respective dagta codings).
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

  # update class of result
  class(result) <- c(stringr::str_c(sort(sources),
                                    sep = "",
                                    collapse = "_"),
                     class(result))

  # return result
  return(result)
}

#' Reverse map from phecodes to Read and ICD 9
#'
#' Requires the output from [map_clinical_events_to_phecodes()], which maps UK
#' Biobank clinical events from Read 2, Read 3 and ICD-9 to ICD-10, then uses
#' [Phecode Map 1.2 with ICD-10 Codes
#' (beta)](https://phewascatalog.org/phecodes_icd10) to map these ICD-10
#' equivalents (and any actual ICD-10 records) to Phecodes. This function uses
#' the output from [map_clinical_events_to_phecodes()] to *reverse* map from
#' phecodes to Read 2, Read 3, ICD-9 and ICD-10. This is useful for checking
#' which raw clinical codes have been used in any phecode-defined phenotypes.
#'
#' @param clinical_events_phecodes A data frame created by
#'   [map_clinical_events_to_phecodes()].
#' @inheritParams codemapper::lookup_codes
#'
#' @return A data frame with columns "phecode", "phecode_description",
#'   "data_coding", "code" "description", "icd10_equivalent" and
#'   "icd10_description". Column 'code' contains the original 'raw' UK clinical
#'   codes, with code type indicated by the 'data_coding' column.
#' @export
make_phecode_reverse_map <- function(clinical_events_phecodes,
                                     all_lkps_maps) {

  start_time <- proc.time()
  # append code types
  clinical_events_phecodes <- ukbwranglr::clinical_events_sources() %>%
    dplyr::select(tidyselect::all_of(c("source",
                                       "data_coding"))) %>%
    dplyr::right_join(clinical_events_phecodes,
                      by = "source")

  # make `code` = `icd10_code` (currently `code` is `NA` if source uses icd10
  # coding)
  clinical_events_phecodes <- clinical_events_phecodes %>%
    dplyr::mutate("code" = dplyr::case_when(
      .data[["data_coding"]] == "icd10" ~ .data[["icd10"]],
      TRUE ~ .data[["code"]]
    ))

  # select required cols only
  clinical_events_phecodes <- clinical_events_phecodes %>%
    dplyr::select(tidyselect::all_of(c(
      "code",
      "data_coding",
      "icd10",
      "phecode"
    )))

  # check there are no columns with missing values (there should be none)
  na_summary <- clinical_events_phecodes %>%
    purrr::map_dbl(~ sum(is.na(.x)))

  cols_with_missing_values <- subset(na_summary,
                                     na_summary > 0)

  assertthat::assert_that(
    length(cols_with_missing_values) == 0,
    msg = paste0(
      "The following columns contain `NA` values (there should be none): ",
      stringr::str_c(names(cols_with_missing_values),
                     sep = "",
                     collapse = ", ")
    )
  )

  # get distinct phecode/code/combinations only
  clinical_events_phecodes <- clinical_events_phecodes %>%
    dplyr::distinct()

  # append code descriptions
  clinical_events_phecodes_split <- split(
    clinical_events_phecodes,
    clinical_events_phecodes$data_coding
  )

  code_descriptions <- clinical_events_phecodes_split %>%
    purrr::imap(~ codemapper::lookup_codes(codes = .x$code,
                                    code_type = .y,
                                    all_lkps_maps = all_lkps_maps,
                                    preferred_description_only = TRUE,
                                    standardise_output = TRUE,
                                    unrecognised_codes = "warning",
                                    col_filters = codemapper::default_col_filters()) %>%
           dplyr::select(tidyselect::all_of(c(
             "code",
             "description"
           ))))

  clinical_events_phecodes <- clinical_events_phecodes_split %>%
    # append code descriptions
    purrr::imap( ~ .x %>%
                   dplyr::left_join(code_descriptions[[.y]],
                                    by = "code")) %>%

    # combine
    dplyr::bind_rows()

  # append phecode_description
  phecode_descriptions <- codemapper::lookup_codes(codes = clinical_events_phecodes$phecode,
                                                   code_type = "phecode",
                                                   all_lkps_maps = all_lkps_maps,
                                                   preferred_description_only = TRUE,
                                                   standardise_output = TRUE,
                                                   unrecognised_codes = "error",
                                                   col_filters = codemapper::default_col_filters()) %>%
    dplyr::select(
      "phecode" = .data[["code"]],
      "phecode_description" = .data[["description"]]
    )

  clinical_events_phecodes <- clinical_events_phecodes %>%
    dplyr::left_join(phecode_descriptions,
              by = c("phecode"))

  # append icd10 equivalent descriptions
  icd10_descriptions <- codemapper::lookup_codes(codes = clinical_events_phecodes$icd10,
                                                 code_type = "icd10",
                                                 all_lkps_maps = all_lkps_maps,
                                                 preferred_description_only = TRUE,
                                                 standardise_output = TRUE,
                                                 unrecognised_codes = "warning",
                                                 col_filters = codemapper::default_col_filters()) %>%
    dplyr::select(
      "icd10" = .data[["code"]],
      "icd10_description" = .data[["description"]]
    ) %>%

    # any warnings will have already appeared above
    suppressWarnings()

  clinical_events_phecodes <- clinical_events_phecodes %>%
    dplyr::left_join(icd10_descriptions,
              by = c("icd10"))

  # reorder cols
  clinical_events_phecodes <- clinical_events_phecodes %>%
    dplyr::select(tidyselect::all_of(
      c(
        "phecode",
        "phecode_description",
        "data_coding",
        "code",
        "description"
      )
    ),
    "icd10_equivalent" = .data[["icd10"]],
    .data[["icd10_description"]])

  # TODO - codes_like() function, use this to make get_undivided_3char_icd10()
  # function (use this to re-append 'X' and avoid warnings for these codes
  # above)

  ukbwranglr:::time_taken_message(start_time)
  return(clinical_events_phecodes)
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



#' Maps ICD10 codes in `clinical_events` to phecodes
#'
#' @param clinical_events Clinical events data frame
#' @param all_lkps_maps Named list of lookup and mapping tables (can be
#'   `tbl_dbi` objects) or the path to an SQLite database containing these
#'
#' @return
#' @noRd
map_icd10_to_phecode <- function(clinical_events,
                                 all_lkps_maps = "all_lkps_maps.db") {
  map_codes_ukb_clinical_events(
    clinical_events = clinical_events,
    from = "icd10",
    to = "phecode",
    from_colname = "icd10",
    to_colname = "phecode",
    all_lkps_maps = all_lkps_maps,
    strict_ukb = FALSE
  )
}
