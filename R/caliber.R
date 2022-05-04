
# PUBLIC ------------------------------------------------------------------

#' Get a mapping table of CALIBER phenotypes/categories
#'
#' Includes a category order and colour for plotting. These were manually
#' assigned based on those used in the [PheWAS::PheWAS] package (see
#' [PheWAS::pheinfo]). The original CALIBER phenotype to categories mapping file
#' was obtained from CALIBER's [github
#' repo](https://github.com/spiros/chronological-map-phenotypes/blob/master/phenotype_categories_mapping.csv).
#'
#' @return A data frame
#' @export
#'
#' @family CALIBER
#' @examples
#' get_caliber_categories_mapping()
get_caliber_categories_mapping <- function() {
  readr::read_csv(file = system.file("extdata", "caliber_phenotype_categories_mapping.csv", package = "ukbwranglrextra"),
                  col_types = list(
                    phenotype = readr::col_character(),
                    category = readr::col_character(),
                    groupnum = readr::col_integer(),
                    colour = readr::col_character()
                  ))
}

#' Download the caliber github repository
#'
#' Downloads to `tempdir()` and unzips, invisibly returning the file path to the
#' unzipped folder.
#'
#' @param url Download URL (e.g.
#'   "https://github.com/spiros/chronological-map-phenotypes/archive/refs/heads/master.zip")
#'
#' @return File path to downloaded (and unzipped) repository, invisibly.
#' @family CALIBER
download_caliber_repo <- function(url) {

  # file paths
  caliber_repo_zip <- tempfile()
  caliber_repo_unzipped <- file.path(tempdir(), tempfile())

  # download zip file
  utils::download.file(url,
                destfile = caliber_repo_zip)

  # unzip
  utils::unzip(caliber_repo_zip,
               exdir = caliber_repo_unzipped)

  # return path to downloaded and unzipped directory
  invisible(list.files(caliber_repo_unzipped))
}

#' Read a local copy of the CALIBER repository into R
#'
#' @description All CALIBER codes are read into a named list containing 3 data
#'   frames: primary care Read 2, secondary care ICD10 and secondary care OPCS4
#'   codes.
#'
#' @details The directory supplied to `caliber_dir_path` is expected to contain
#'   subdirectories `primary_care` and `secondary_care`, each of which contains
#'   csv files with clinical code lists.
#'
#'   Note also that:
#'
#'   1. Medcodes are dropped
#'
#'   2. Minimal reformatting is performed so that primary and secondary care
#'   codes may be combined into a single data frame
#'
#' @param caliber_dir_path Path to a locally downloaded copy of the [CALIBER
#'   github repository](https://github.com/spiros/chronological-map-phenotypes).
#' @param overlapping_disease_categories If 'error' (default), raises an error
#'   if any overlapping disease categories are present after mapping. Specify
#'   'warning' to raise a warning instead.
#'
#' @export
#' @return A named list of data frames.
#' @family CALIBER
#' @examples
#' # read local copy of CALIBER repository into a named list. Note that
#' # (i) Medcodes are dropped and (ii) minimal reformatting is performed so that
#' # primary and secondary care codes may be combined into a single data frame
#' caliber_raw <- read_caliber_raw(dummy_caliber_dir_path())
#' caliber_raw
#'
#' # combine into a single data frame with dplyr
#' dplyr::bind_rows(caliber_raw)
read_caliber_raw <- function(caliber_dir_path,
                             overlapping_disease_categories = "error") {
  # validate args
  match.arg(overlapping_disease_categories,
            c("error", "warning"))

  # Set filepath constants --------------------------------------------------

  CALIBER_PRIMARY <- file.path(caliber_dir_path, "primary_care")
  CALIBER_SECONDARY <- file.path(caliber_dir_path, "secondary_care")
  CSV_REGEX <- "+\\.csv$"

  PRIMARY_CARE_FILES <- list.files(CALIBER_PRIMARY,
                                   pattern = CSV_REGEX)

  SECONDARY_CARE_FILES <- list.files(CALIBER_SECONDARY,
                                     pattern = CSV_REGEX)

  SECONDARY_CARE_FILES_ICD <-
    subset(SECONDARY_CARE_FILES, grepl("^ICD_", SECONDARY_CARE_FILES))
  SECONDARY_CARE_FILES_OPCS <-
    subset(SECONDARY_CARE_FILES, grepl("^OPCS_", SECONDARY_CARE_FILES))

  # Read CALIBER files and reformat -----------------------------------------

  # Note - currently removes medcodes and secondary descriptions for read codes

  # Read files into 3 dataframes - primary care and secondary care (ICD and OPCS)
  caliber <- c(
    "read2",
    "icd10",
    "opcs4"
  ) %>%
    purrr::set_names() %>%
    purrr::map(~ NULL)

  message("Reading CALIBER clinical codes lists into R")
  message("Primary care Read 2 (1 of 3)")
  caliber$read2 <-
    read_csv_to_named_list_and_combine(CALIBER_PRIMARY,
                                       filenames = PRIMARY_CARE_FILES,
                                       standardising_function = standardise_primary_care)

  message("Secondary care ICD10 (2 of 3)")
  caliber$icd10 <-
    read_csv_to_named_list_and_combine(CALIBER_SECONDARY,
                                       filenames = SECONDARY_CARE_FILES_ICD,
                                       standardising_function = standardise_secondary_care_icd10)

  message("Secondary care OPCS4 (3 of 3)")
  caliber$opcs4 <-
    read_csv_to_named_list_and_combine(CALIBER_SECONDARY,
                                       filenames = SECONDARY_CARE_FILES_OPCS,
                                       standardising_function = standardise_secondary_care_opcs4)

  # Check for overlapping disease categories -------
  overlapping_disease_categories_list <- ukbwranglr:::identify_overlapping_disease_categories(dplyr::bind_rows(caliber))

  if (!is.null(overlapping_disease_categories_list)) {
    overlapping_disease_categories_msg <- paste0("The following ",
                                                 length(unique(overlapping_disease_categories_list$clinical_codes$disease)),
                                                 " diseases include categories with non-distinct codes: ",
                                                 stringr::str_c(unique(overlapping_disease_categories_list$clinical_codes$disease),
                                                                sep = "",
                                                                collapse = ", "))

    switch(overlapping_disease_categories,
           error = stop(overlapping_disease_categories_msg),
           warning = warning(overlapping_disease_categories_msg))
  }

  return(caliber)
}

#' Reformat and map CALIBER codes for use with UK Biobank data
#'
#' Reformats Read 2, ICD10 and OPCS4 CALIBER codes to match the format in UK
#' Biobank data, and also maps from Read 2 to Read 3, as well as from ICD10 to
#' ICD9. See `vignette("caliber")` for further details.
#'
#' @param caliber A named list of data frames, created by [read_caliber_raw()].
#' @inheritParams codemapper::codes_starting_with
#' @inheritParams read_caliber_raw
#' @param overlapping_disease_categories_csv File path to a csv containing codes
#'   that are listed under more than one disease category within a disease. This
#'   should have the same format as [ukbwranglr::example_clinical_codes()], with
#'   the author column set to 'caliber' for all rows, plus an additional 'keep'
#'   column with 'Y' values indicating which rows to keep. By default, this is
#'   set to [default_overlapping_disease_categories_csv()].
#'
#' @return A named list of data frames.
#' @export
#' @family CALIBER
#' @examples
#' # read local copy of CALIBER repository into a named list
#' caliber_raw <- read_caliber_raw(dummy_caliber_dir_path())
#'
#' # reformat
#' # reformat_caliber_for_ukb(caliber_raw, all_lkps_maps = codemapper::build_dummy_all_lkps_maps())
reformat_caliber_for_ukb <- function(caliber,
                                     all_lkps_maps,
                                     col_filters = codemapper::default_col_filters(),
                                     overlapping_disease_categories = "error",
                                     overlapping_disease_categories_csv = default_overlapping_disease_categories_csv()) {
  # validate args - TODO


  # connect to database file path if `all_lkps_maps` is a string, or `NULL`
  if (is.character(all_lkps_maps)) {
    con <- codemapper:::check_all_lkps_maps_path(all_lkps_maps)
    all_lkps_maps <- ukbwranglr::db_tables_to_list(con)
    on.exit(DBI::dbDisconnect(con))
  } else if (is.null(all_lkps_maps)) {
    if (Sys.getenv("ALL_LKPS_MAPS_DB") != "") {
      message(paste0("Attempting to connect to ", Sys.getenv("ALL_LKPS_MAPS_DB")))
      con <-
        con <- codemapper:::check_all_lkps_maps_path(Sys.getenv("ALL_LKPS_MAPS_DB"))
      all_lkps_maps <- ukbwranglr::db_tables_to_list(con)
      on.exit(DBI::dbDisconnect(con))
    } else if (file.exists("all_lkps_maps.db")) {
      message("Attempting to connect to all_lkps_maps.db in current working directory")
      con <- codemapper:::check_all_lkps_maps_path("all_lkps_maps.db")
      all_lkps_maps <- ukbwranglr::db_tables_to_list(con)
      on.exit(DBI::dbDisconnect(con))
    } else {
      stop(
        "No/invalid path supplied to `all_lkps_maps` and no file called 'all_lkps_maps.db' found in current working directory. See `?all_lkps_maps_to_db()`"
      )
    }
  }

  # reformat read2, icd10 and opcs4 codes --------
  message("Reformatting Read 2 codes")
  caliber$read2 <-
    reformat_caliber_read2(caliber$read2,
                           all_lkps_maps = all_lkps_maps)

  message("Reformatting ICD10 codes")
  caliber$icd10 <-
    reformat_caliber_icd10(caliber$icd10,
                           all_lkps_maps = all_lkps_maps)

  message("Reformatting OPCS4 codes")
  caliber$opcs4 <-
    reformat_caliber_opcs4(caliber$opcs4,
                           all_lkps_maps = all_lkps_maps)


  # Map codes (read2 - read3,  icd10 - icd9) --------------------------------
  message("Mapping read2 codes to read3")

  caliber$read3 <- map_caliber(df = caliber$read2,
                                                  from = "read2",
                                                  to = "read3",
                                                  all_lkps_maps = all_lkps_maps,
                                                  col_filters = col_filters)

  message("Mapping icd10 to icd9 codes")
  caliber$icd9 <- map_caliber(df = caliber$icd10,
                                                   from = "icd10",
                                                   to = "icd9",
                                                   all_lkps_maps = all_lkps_maps,
                                                   col_filters = col_filters)

  # Remove 'X' from ends of 3 character ICD10 codes
  caliber$icd10 <- caliber$icd10 %>%
    dplyr::mutate("code" = stringr::str_remove(.data[["code"]],
                                                "X$"))

  # Resolve overlapping disease categories ------------

  caliber <- dplyr::bind_rows(caliber)

  if (!is.null(overlapping_disease_categories_csv)) {
    # user-supplied df indicating how to handle overlapping disease categories
    resolved_overlapping_disease_categories <-
      readr::read_csv(
        overlapping_disease_categories_csv,
        progress = FALSE,
        col_types = readr::cols(.default = "c")
      )

    # filter for only rows present in `caliber`
    resolved_overlapping_disease_categories <-
      resolved_overlapping_disease_categories %>%
      dplyr::semi_join(caliber,
                       by = names(caliber))

    # validate
    resolved_overlapping_disease_categories <-
      validate_overlapping_disease_categories_df(resolved_overlapping_disease_categories)

    # resolve overlapping disease categories

    ## setup
    resolved_overlapping_disease_categories <-
      resolved_overlapping_disease_categories %>%
      dplyr::mutate(keep = dplyr::case_when(is.na(.data[["keep"]]) ~ "N",
                                            TRUE ~ .data[["keep"]]))

    resolved_overlapping_disease_categories <-
      split(
        resolved_overlapping_disease_categories,
        resolved_overlapping_disease_categories$keep
      )

    ## remove rows to be removed
    caliber <- caliber %>%
      dplyr::anti_join(resolved_overlapping_disease_categories$N,
                       by = names(caliber))
  }

  # Check for remaining overlapping disease categories ----------
  overlapping_disease_categories_list <- ukbwranglr:::identify_overlapping_disease_categories(caliber)

  if (!is.null(overlapping_disease_categories_list)) {
    overlapping_disease_categories_msg <- paste0("The following ",
                                                 length(unique(overlapping_disease_categories_list$clinical_codes$disease)),
                                                 " diseases include categories with non-distinct codes: ",
                                                 stringr::str_c(unique(overlapping_disease_categories_list$clinical_codes$disease),
                                                                sep = "",
                                                                collapse = ", "))

    switch(overlapping_disease_categories,
           error = stop(overlapping_disease_categories_msg),
           warning = warning(overlapping_disease_categories_msg))
  }

  # return result
  return(caliber)
}


#' CSV file for resolving overlapping disease categories in mapped CALIBER codes
#'
#' Returns the file path to csv included with this package for resolving
#' overlapping disease categories in mapped CALIBER codes. To be used with
#' [reformat_caliber_for_ukb()].
#'
#' @return A file path.
#' @export
#'
#' @family CALIBER
#' @examples
#' # return file path
#' default_overlapping_disease_categories_csv()
#'
#' # read file into R
#' read.csv(default_overlapping_disease_categories_csv())
default_overlapping_disease_categories_csv <- function() {
  system.file("extdata",
              "caliber_mapped_overlapping_disease_categories.csv",
              package = "ukbwranglrextra")
}

# PRIVATE ---------------------------------------------------------------



## Read and reformat CALIBER codes --------------------------------------------------

### Read CALIBER (raw) --------------------

#### Standardising mappers -------------

# standardising functions for primary and secondary care (ICD and OPCS4) csv files
standardise_primary_care <- function(df) {
  df <- df %>%
    tidyr::pivot_longer(
      cols = c("Readcode", "Medcode"),
      names_to = "code_type",
      values_to = "code"
    ) %>%
    dplyr::mutate("author" = "caliber") %>%
    dplyr::select(tidyselect::all_of(
      c(
        "Disease",
        "ReadcodeDescr",
        "Category",
        "code_type",
        "code",
        "author"
      )
    ))

  df <- ukbwranglr:::rename_cols(
    df,
    old_colnames = c("Disease", "ReadcodeDescr", "Category"),
    new_colnames = c("disease", "description", "category")
  )

  return(df)
}

standardise_secondary_care_icd10 <- function(df) {
  df <- df %>%
    dplyr::mutate(code_type = "icd10",
                  author = "caliber") %>%
    dplyr::select(tidyselect::all_of(
      c(
        "Disease",
        "ICD10codeDescr",
        "Category",
        "code_type",
        "ICD10code",
        "author"
      )
    ))

  df <- ukbwranglr:::rename_cols(
    df,
    old_colnames = c("Disease", "ICD10codeDescr", "Category", "ICD10code"),
    new_colnames = c("disease", "description", "category", "code")
  )

  return(df)
}

standardise_secondary_care_opcs4 <- function(df) {
  df <- df %>%
    dplyr::mutate(code_type = "opcs4",
                  author = "caliber") %>%
    dplyr::select(tidyselect::all_of(
      c(
        "Disease",
        "OPCS4codeDescr",
        "Category",
        "code_type",
        "OPCS4code",
        "author"
      )
    ))

  df <- ukbwranglr:::rename_cols(
    df,
    old_colnames = c("Disease", "OPCS4codeDescr", "Category", "OPCS4code"),
    new_colnames = c("disease", "description", "category", "code")
  )

  return(df)
}

#### Read all csvs to named list --------

# reads a list of csv files into a named list, standardises, then combines into single df
read_csv_to_named_list_and_combine <-
  function(# directory where files are located
    directory,
    # vector of file names
    filenames,
    # function to process each file with
    standardising_function) {
    result <- as.list(file.path(directory, filenames))
    names(result) <- filenames

    result %>%
      purrr::map(~ readr::read_csv(.x,
                                   progress = FALSE,
                                   col_types = readr::cols(.default = "c"))) %>%
      purrr::map(standardising_function) %>%
      dplyr::bind_rows()
  }

### Reformat CALIBER --------------------------------------------------------


#' Reformat CALIBER Read 2/Medcode df
#'
#' Filter for Read codes only (i.e. remove Medcodes). Remove last 2 characters from Read codes (these indicate whether a description is primary or not).
#'
#' @param read2_df A caliber df of Read2/Medcodes, read using [read_caliber_raw()].
#' @param all_lkps_maps all_lkps_maps list
#' @param unrecognised_codes Passed to `codemapper::check_codes_exist()`
#'
#' @return Dataframe
#' @noRd
reformat_caliber_read2 <- function(read2_df,
                                   all_lkps_maps,
                                   unrecognised_codes = "warning") {

    read2_df <- read2_df %>%
    # filter for only read2 codes
    dplyr::filter(.data[["code_type"]] == "Readcode") %>%
    # label as 'read2' (ukbwranglr format)
    dplyr::mutate("code_type" = "read2") %>%

    # make column with last 2 characters of `code`
    dplyr::mutate("code_last_2_char" = stringr::str_sub(.data[["code"]],
                                                        start = -2L,
                                                        end = -1L)) %>%

    # remove last 2 characters
    dplyr::mutate("code" = stringr::str_sub(.data[["code"]],
                                          start = 1L,
                                          end = -3L)) %>%

    # take only one description per code, per disease
    dplyr::group_by(.data[["disease"]],
                    # .data[["category"]],
                    .data[["code"]]) %>%
    dplyr::arrange(.data[["code_last_2_char"]]) %>%
    dplyr::slice(1L) %>%
    dplyr::ungroup() %>%
    dplyr::select(-.data[["code_last_2_char"]])

  # check for unrecognised codes - warning if any found
  read_v2_lkp <- all_lkps_maps$read_v2_lkp %>%
    dplyr::collect()

  codemapper:::check_codes_exist(codes = read2_df$code,
                                 lkp_codes = read_v2_lkp$read_code,
                                 table_name = "read_v2_lkp",
                                 code_type = "read2",
                                 return_unrecognised_codes = FALSE,
                                 unrecognised_codes = "warning")

  # return result
  return(read2_df)
}

#' Reformat CALIBER ICD10 df
#'
#' Convert ICD10 codes from ICD10_CODE format to ALT_CODE format (note that
#' undivided 3 character ICD10 codes will be appended by 'X' in the final output
#' from this function). Any unrecognised ICD10 codes will be removed (note that
#' only A90 and A91 are unrecognised). Expand 3 character ICD10 codes to include
#' children (e.g. D25 expanded to include D250, D251, D252 and D529)
#'
#' @param icd10_df A caliber df of Read2/Medcodes, read using
#'   [read_caliber_raw()].
#' @param all_lkps_maps all_lkps_maps list
#' @param unrecognised_codes Passed to [codemapper::reformat_icd10_codes()]
#'
#' @return Dataframe
#' @noRd
reformat_caliber_icd10 <- function(icd10_df,
                                   all_lkps_maps,
                                   unrecognised_codes = "warning") {

  # get vector of all present icd10 codes in ALT_CODE format. Also returns a
  # message listing ICD10 codes with modifiers that will map to >1 ICD10 code in
  # ALT_CODE format. Also raises warning if any unrecognised ICD10 codes are
  # present.
  icd10_codes_in_icd10_df <- codemapper::reformat_icd10_codes(
    icd10_codes = icd10_df$code,
    all_lkps_maps = all_lkps_maps,
    input_icd10_format = "ICD10_CODE",
    output_icd10_format = "ALT_CODE",
    unrecognised_codes = unrecognised_codes,
    strip_x = FALSE
  )

  cols_to_keep <- names(icd10_df)

  icd10_lkp_map <- all_lkps_maps$icd10_lkp %>%
    dplyr::select(tidyselect::all_of(c(
      "ICD10_CODE",
      "ALT_CODE"
    ))) %>%
    dplyr::collect()

  icd10_df <- icd10_lkp_map %>%
    dplyr::filter(.data[["ALT_CODE"]] %in% !!icd10_codes_in_icd10_df) %>%
    dplyr::right_join(icd10_df,
                      by = c("ICD10_CODE" = "code")) %>%
    dplyr::mutate("code" = .data[["ALT_CODE"]]) %>%
    dplyr::select(tidyselect::all_of(cols_to_keep))

  # remove any rows with unrecognised codes
  unrecognised_icd10_df <- icd10_df %>%
    dplyr::filter(is.na(.data[["code"]]))

  warning(paste0("Removing ",
                nrow(unrecognised_icd10_df),
                " rows with unrecognised ICD10 codes. Diseases with unrecognised codes: '",
                stringr::str_c(unique(unrecognised_icd10_df$disease),
                               sep = "",
                               collapse = ", "),
                "'"))

  icd10_df <- icd10_df %>%
    dplyr::filter(!is.na(.data[["code"]]))

  # expand 3 character ICD10 codes (e.g. E10 - see warning note under 'Code
  # list' tab: https://www.caliberresearch.org/portal/show/diabcomp_hes)
  icd10_3_char <- icd10_df %>%
    dplyr::filter(stringr::str_length(.data[["code"]]) == 3)

  icd10_lkp_map_3_char <- icd10_lkp_map %>%
    dplyr::mutate("icd10_3_char" = stringr::str_sub(.data[["ICD10_CODE"]],
                                                    start = 1L,
                                                    end = 3L)) %>%
    dplyr::filter(.data[["icd10_3_char"]] %in% !!icd10_3_char$code) %>%
    dplyr::select(tidyselect::all_of(c(
      "icd10_3_char",
      "ALT_CODE"
    )))

  icd10_3_char <- icd10_3_char %>%
    dplyr::left_join(icd10_lkp_map_3_char,
                     by = c("code" = "icd10_3_char")) %>%
    dplyr::mutate("code" = .data[["ALT_CODE"]]) %>%
    dplyr::select(tidyselect::all_of(cols_to_keep))

  # recombine and remove duplicate rows
  icd10_df <- dplyr::bind_rows(icd10_df,
                               icd10_3_char) %>%
    dplyr::distinct()

  # return result
  return(icd10_df)
}

#' Reformat CALIBER OPCS4 df
#'
#' Filter for Read codes only (i.e. remove Medcodes). Remove last 2 characters
#' from Read codes (these indicate whether a description is primary or not).
#'
#' @param opcs4_df A caliber df of Read2/Medcodes, read using
#'   [read_caliber_raw()].
#' @param all_lkps_maps all_lkps_maps list
#' @param unrecognised_codes Passed to `codemapper:::check_codes_exist()`
#'
#' @return Dataframe
#' @noRd
reformat_caliber_opcs4 <- function(opcs4_df,
                                   all_lkps_maps,
                                   unrecognised_codes = "warning") {

  # remove '.'
  opcs4_df <- opcs4_df %>%
    dplyr::mutate("code" = stringr::str_remove(.data[["code"]],
                                               pattern = "\\."))

  # check for unrecognised codes - warning if any found
  opcs4_lkp <- all_lkps_maps$opcs4_lkp %>%
    dplyr::collect()

  codemapper:::check_codes_exist(codes = opcs4_df$code,
                                 lkp_codes = opcs4_lkp$opcs4_code,
                                 code_type = "read2",
                                 return_unrecognised_codes = FALSE,
                                 unrecognised_codes = "warning")

  # return result
  return(opcs4_df)
}

### Map CALIBER -------------------------------------------------------------

# helper functions for `reformat_caliber_for_ukb()` to map codes from read2 to
# read3 and icd10 to icd9


#' Map CALIBER codelist data frame from one coding system to another
#'
#' Helper function for [reformat_caliber_for_ukb()]. Maps a codes in a
#' reformatted CALIEBR codelist data frame, used for Read 2 to Read 3, and ICD10
#' to ICD9.
#'
#' @param df Reformatted CALIBER code list, data frame
#' @param from Code type to map from
#' @param to Code type to map to
#' @param all_lkps_maps Named list
#' @param col_filters See [codemapper::map_codes()]
#'
#' @return A data frame. The `code_type` column will equal the value for `to`.
map_caliber <- function(df,
                        from,
                        to,
                        all_lkps_maps,
                        col_filters) {
  final_col_order <- names(df)

  # get mapping df
  mapping_df <- codemapper::get_mapping_df(from = from,
                               to = to,
                               all_lkps_maps = all_lkps_maps,
                               col_filters = col_filters,
                               rename_from_to = c(from = "old_code",
                                                  to = "new_code"),
                               reverse_mapping = "warning")

  # join, update `code_type` and rename cols
  result <- df %>%
    dplyr::inner_join(mapping_df,
                     by = c("code" = "old_code")) %>%
    dplyr::mutate("code_type" = !!to) %>%
    dplyr::select(-.data[["code"]]) %>%
    ukbwranglr:::rename_cols(old_colnames = "new_code",
                                    new_colnames = "code")

  # distinct codes only (if there are duplicate codes because of alternative
  # code descriptions, just one code description should be selected)
  result <- result %>%
    dplyr::distinct(dplyr::across(-.data[["description"]]))

  # append code descriptions
  code_descriptions <- codemapper::lookup_codes(codes = result$code,
                                                code_type = to,
                                                all_lkps_maps = all_lkps_maps,
                                                preferred_description_only = TRUE,
                                                standardise_output = TRUE,
                                                unrecognised_codes = "error",
                                                col_filters = col_filters) %>%
    dplyr::select(-.data[["code_type"]])

  result <- result %>%
    dplyr::left_join(code_descriptions,
                     by = "code")

  # reorder
  result %>%
    dplyr::select(tidyselect::all_of(final_col_order))
}

# Validation helpers ------------------------------------------------------

#' Validation helper for [reformat_caliber_for_ukb()]
#'
#' Validate a user-supplied df which indicates how to handle overlapping disease
#' categories.
#'
#' @param df Data frame - columns as for a clinical codes df, plus an additional
#'   'keep' column which should only contain 'Y' and `NA` values.
#'
#' @return The input df, unchanged
validate_overlapping_disease_categories_df <- function(df) {
  # expected names
  expected_colnames <- c("disease",
                         "description",
                         "category",
                         "code_type",
                         "code",
                         "author",
                         "keep")

  assertthat::assert_that(all(names(df) == expected_colnames),
                          msg = paste0("Unexpected colnames detected. Expected colnames: ",
                                       stringr::str_c(expected_colnames,
                                                      sep = "",
                                                      collapse = ", ")))

  # all type character
  coltypes <- df %>%
    purrr::map_chr(class) %>%
    unique()

  assertthat::assert_that(length(unique(coltypes)) == 1 && coltypes == "character",
                          msg = "All columns should be of type character")

  # author column only includes 'caliber'
  assertthat::assert_that(length(unique(df$author)) == 1 && unique(df$author) == "caliber",
                          msg = "Column 'author' should only contain value 'caliber'")

  # keep col only contains 'Y' or NA
  keep_col_values <- df %>%
    dplyr::distinct(.data[["keep"]]) %>%
    dplyr::pull(.data[["keep"]])

  assertthat::assert_that(keep_col_values[1] == c("Y") & is.na(keep_col_values[2]),
                          msg = "Column 'keep' should only contain 'Y' and `NA` values (at least one of each)")

  # each code has one row labelled 'keep', and one or more rows not to keep
  df_code_split <- split(df,
                         paste(df$code_type, df$code, sep = "_"))

  df_code_split_n_y <- df_code_split %>%
    purrr::map_dbl(~ sum(.x$keep == "Y",
                         na.rm = TRUE))

  df_code_split_n_y_gt_1 <- subset(df_code_split_n_y,
                                   df_code_split_n_y > 1)

  df_code_split_n_y_eq_0 <- subset(df_code_split_n_y,
                                   df_code_split_n_y == 0)

  assertthat::assert_that(
    length(df_code_split_n_y_gt_1) == 0 &&
      length(df_code_split_n_y_eq_0) == 0,
    msg = paste0(
      "All codes should have exactly one row flagged as 'Y' in column 'keep'. ",
      length(c(
        df_code_split_n_y_gt_1, df_code_split_n_y_eq_0
      )),
      " codes do not meet this requirement. These are (displaying up to the first 25): ",
      stringr::str_c(utils::head(names(
        c(df_code_split_n_y_gt_1, df_code_split_n_y_eq_0)
      ),
      n = 25),
      sep = "",
      collapse = ", ")
    )
  )

  # no missing values (excluding 'keep' column)
  df_minus_keep_col <- df %>%
    dplyr::select(-.data[["keep"]])

  if (!assertthat::are_equal(df_minus_keep_col,
                             stats::na.omit(df_minus_keep_col))) {
    stop("There should be no missing values (except in column 'keep')")
  }

  return(df)
}
