
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
#' Downloads to `tempdir()` and unzips.
#'
#' @param url
#'
#' @return File path to downloaded (and unzipped) repository.
download_caliber_repo <- function(
  url = "https://github.com/spiros/chronological-map-phenotypes/archive/refs/heads/master.zip") {

  # file paths
  caliber_repo_zip <- tempfile()
  caliber_repo_unzipped <- file.path(tempdir(), tempfile())

  # download zip file
  download.file(url,
                destfile = caliber_repo_zip)

  # unzip
  utils::unzip(caliber_repo_zip,
               exdir = caliber_repo_unzipped)

  # return path to downloaded and unzipped directory
  return(list.files(caliber_repo_unzipped))
}

# OVERVIEW ----------------------------------------------------------------

# Functions to collect all CALIBER code lists into a single standardised table, map
# codes from read2 -> read3 and from icd10 ->icd9, and write this to a SQLite
# database.

# Medcodes are dropped and only primary descriptions are kept for read2 and read3


# FUNCTIONS ---------------------------------------------------------------

map_caliber_for_ukb <- function() {
# read_csv() - all columns read as type character
read_csv_as_character <- purrr::partial(readr::read_csv, col_types = readr::cols(.default = "c"))

# standardising functions for primary and secondary care (ICD and OPCS4) csv files
standardise_primary_care <- purrr::as_mapper(
  ~ .x %>%
    tidyr::pivot_longer(
      cols = c("Readcode", "Medcode"),
      names_to = "code_type",
      values_to = "code"
    ) %>%
    dplyr::mutate(
      author = "caliber"
    ) %>%
    dplyr::select(
      disease = Disease,
      description = ReadcodeDescr,
      category = Category,
      code_type,
      code,
      author
    )
)

standardise_secondary_care_icd10 <- purrr::as_mapper(
  ~ .x %>%
    dplyr::mutate(code_type = "icd10",
                  author = "caliber") %>%
    dplyr::select(
      disease = Disease,
      description = ICD10codeDescr,
      category = Category,
      code_type,
      code = ICD10code,
      author
    )
)

standardise_secondary_care_opcs4 <- purrr::as_mapper(
  ~ .x %>%
    dplyr::mutate(code_type = "opcs4",
                  author = "caliber") %>%
    dplyr::select(
      disease = Disease,
      description = OPCS4codeDescr,
      category = Category,
      code_type,
      code = 	OPCS4code,
      author
    )
)

# functions to reformat codes for UKB

get_icd10_codes_with_modifiers <- function(icd10_lkp) {
  icd10_lkp %>%
    dplyr::collect() %>%
    dplyr::filter(!is.na(MODIFIER_4) | !is.na(MODIFIER_5)) %>%
    # mutate column for codes minus modifiers
    dplyr::mutate(base_code = stringr::str_replace_all(
      ICD10_CODE,
      "\\.[:digit:]*$",
      ""
    ))
}

append_icd10_ALT_CODEs_with_modifiers <- function(df,
                                                  all_lkps_maps) {

  # approach - identify icd10 codes with modifiers and separate these from those
  # that do not. For each disease category (codes should only appear once per disease category)

  # all icd10 codes with modifiers
  icd10_with_modifiers_all <- get_icd10_codes_with_modifiers(all_lkps_maps$icd10_lkp)

  # icd10 codes (base e.g. E10) with modifiers (e.g. E10.0, E10.1 etc) in df
  df_icd10_with_modifiers <- df %>%
    dplyr::filter(code %in% icd10_with_modifiers_all$base_code)

  # ...remove these from df
  df <- df %>%
    dplyr::filter(!code %in% icd10_with_modifiers_all$base_code)

  # append ALT_CODES with modifiers, by disease
  disease_categories <- unique(df_icd10_with_modifiers$category)

  df_icd10_with_modifiers <- disease_categories %>%
    purrr::set_names() %>%
    purrr::map(.f = ~ {
      # get base_codes for single disease category
      base_codes_to_expand <- df_icd10_with_modifiers %>%
        dplyr::filter(category == .x) %>%
        .$code

      # get disease for disease category
      DISEASE <- df_icd10_with_modifiers %>%
        dplyr::filter(category == .x) %>%
        .$disease %>%
        unique() # should be unique

      # get full set of codes + modifiers for these 'base_codes'
      expanded_codes <- icd10_with_modifiers_all %>%
        dplyr::filter(base_code %in% base_codes_to_expand) %>%
        .$ALT_CODE

      # add in base codes
      expanded_codes <- c(expanded_codes, base_codes_to_expand)

      # get descriptions
      codemapper::lookup_codes(
        codes = expanded_codes,
        code_type = "icd10",
        all_lkps_maps = all_lkps_maps,
        preferred_description_only = TRUE,
        standardise_output = TRUE,
        quiet = FALSE
      ) %>%
        dplyr::mutate(
          disease = DISEASE,
          code_type = "icd10",
          author = "caliber",
          category = .x
        ) %>%
        dplyr::select(disease,
                      description,
                      category,
                      code_type,
                      code,
                      author)
    }) %>%
    dplyr::bind_rows() %>%
    dplyr::distinct()

  # now re-join to original df
  result <- dplyr::bind_rows(df, df_icd10_with_modifiers)

  return(result)
}

# TODO - either delete, or keep as test for append_icd10_ALT_CODEs_with_modifiers
expand_3_char_icd10_codes <- function(icd10_df,
                                      all_lkps_maps) {
  # (see caliber warning note
  # here: https://www.caliberresearch.org/portal/show/diabcomp_hes)

  # validate args
  # all `code_type` should be 'icd10', with no missing values
  assertthat::assert_that((length(unique(icd10_df$code_type)) == 1) &&
                            (unique(icd10_df$code_type) == "icd10"),
                          msg = "Error! `icd10_df` contains code types that are not 'icd10' (or possibly `NA` values)")

  # expand all codes - icd10 codes should be either length 3 or 4
  result <- icd10_df %>%
    dplyr::pull(.data[["code"]]) %>%
    codemapper::codes_starting_with(code_type = "icd10",
                                    all_lkps_maps = all_lkps_maps,
                                    codes_only = FALSE,
                                    preferred_description_only = TRUE,
                                    standardise_output = TRUE) %>%
    dplyr::bind_rows(str_len_3_icd10) %>%
    tidyr::fill(tidyselect::everything(),
                .direction = "up")

  assertthat::assert_that(nrow(result) == nrow(na.omit(result)),
                          msg = "Error when expanding 3 character ICD10 codes! Some cells are `NA`")

  return(result)
}

# read2 codes: filter for only primary descriptions and remove last 2 characters
# (the last 2 characters indicate whether description is primary or not for a
# code) and remove "." from ICD-10 codes
reformat_caliber_read2 <- function(read2_df) {
  read2_df %>%
    # filter for only read2 codes
    dplyr::filter(code_type == "Readcode") %>%
    # label as 'read2' (ukbwranglr format)
    dplyr::mutate(code_type = "read2") %>%

    # remove last 2 characters
    dplyr::mutate(code = stringr::str_sub(code,
                                          start = 1L,
                                          end = -3L)) %>%

    # take only one description per code, per disease
    dplyr::group_by(disease, code) %>%
    dplyr::slice(1L) %>%
    dplyr::ungroup()
}

reformat_caliber_icd10 <- function(icd10_df,
                                   all_lkps_maps) {

  # convert to ALT_CODE format - note, some codes are in ALT_CODE format (e.g.
  # diabetes, 'O24'), while others are in ICD_10 format (e.g. acute kidney
  # injury, 'N17').
  icd10_df <- icd10_df %>%
    dplyr::mutate(
      code = dplyr::case_when(
        stringr::str_detect(code, pattern = "\\.") ~ codemapper::reformat_icd10_codes(icd10_codes = .data[["code"]],
                                                                                      all_lkps_maps = all_lkps_maps,
                                                                                      input_icd10_format = "ICD10_CODE",
                                                                                      output_icd10_format = "ALT_CODE",
                                                                                      strip_x = TRUE),
        TRUE ~ .data[["code"]]
      )
    )

  # TODO - raise issue on caliber repo. Temp fix for single code 'N23.X' here (should just be 'N23' on caliber)
  icd10_df <- icd10_df %>%
    dplyr::mutate(code = dplyr::case_when(
      code == "N23.X" ~ "N23",
      TRUE ~ code
    ))

  # 3 character codes with no children (e.g. 'A38', Scarlet fever, 'I10', Essential (primary) hypertension)
  # need 'X' appended, so that mapping to ICD9 will work
  icd10_format_3_char_codes <- all_lkps_maps$icd10_lkp %>%
    dplyr::collect() %>%
    dplyr::filter(stringr::str_detect(.data[["ALT_CODE"]],
                                      pattern = "X$")) %>%
    dplyr::pull(ICD10_CODE)

  icd10_df <- icd10_df %>%
    dplyr::mutate(
      code = dplyr::case_when(.data[["code"]] %in% icd10_format_3_char_codes ~ paste0(.data[["code"]], "X"),
                              TRUE ~ .data[["code"]])
    )

  # expand icd 10 codes to include all modifiers, where appropriate e.g. E10, I70.0 (see caliber warning note
  # here: https://www.caliberresearch.org/portal/show/diabcomp_hes )
  icd10_df <- append_icd10_ALT_CODEs_with_modifiers(df = icd10_df,
                                                    all_lkps_maps = all_lkps_maps)

  return(icd10_df)
}

# functions to map codes from read2 to read3 and icd10 to icd9
map_caliber_single_disease_category <- function(df,
                                                disease,
                                                disease_category,
                                                all_lkps_maps,
                                                from,
                                                to) {
  # process read codes only - drop medcodes

  # map to read3
  df <- df %>%
    purrr::pluck("code") %>%
    codemapper::map_codes(from = from,
                          to = to,
                          all_lkps_maps = all_lkps_maps,
                          codes_only = FALSE,
                          preferred_description_only = TRUE,
                          standardise_output = TRUE)

  # reformat
  if (is.null(df)) {
    return(NULL)
  }

  codemapper:::reformat_standardised_codelist(
    standardised_codelist = df,
    code_type = to,
    disease = disease,
    disease_category = disease_category,
    author = "caliber"
  )
}

map_caliber_multiple_disease_categories <- function(df,
                                                    all_lkps_maps,
                                                    from,
                                                    to,
                                                    verbose = TRUE) {
  # read2 to read3
  result <- unique(df$category) %>%
    purrr::set_names() %>%
    purrr::map( ~ NULL)

  total_n_disease_categories = length(names(result))
  counter <- 1

  # loop through disease categories by disease (nested for loop)
  diseases <- unique(df$disease)
  for (disease in diseases) {
    disease_df <- df[df$disease == disease, ]

    categories <- unique(disease_df$category)

    # for each disease category, map codes
    for (disease_category in categories) {
      if (verbose) {
        message(paste0("Mapping codes for ",
                       disease_category,
                       ". ",
                       counter, " of ", total_n_disease_categories))
      }

      disease_category_df <- disease_df[disease_df$category == disease_category, ]

      result[[disease_category]] <-
        map_caliber_single_disease_category(
          df = disease_category_df,
          disease = disease,
          disease_category = disease_category,
          all_lkps_maps = all_lkps_maps,
          from = from,
          to = to
        )

      counter <- counter + 1
    }
  }

  # combine list of results into a single df and return
  return(dplyr::bind_rows(result))
}

# reads a list of csv files into a named list, standardises, then combines into single df
read_csv_to_named_list_and_combine <- function(
  directory, # directory where files are located
  filenames, # vector of file names
  standardising_function, # function to process each file with
  file_ext = ".csv", # file extension to remove
  read_function = read_csv_as_character # function read files
) {
  paste(directory, filenames, sep = "/") %>%
    purrr::set_names(nm = stringr::str_replace(filenames,
                                               file_ext,
                                               "")) %>% # remove '.csv'
    purrr::map(read_function) %>%
    purrr::map(standardising_function) %>%
    dplyr::bind_rows()
}

# MAIN --------------------------------------------------------------------

get_caliber_codes_standardise_and_map <- function(all_lkps_maps) {


  # Download caliber repo ---------------------------------------------------
  message("Downloading code lists from caliber repo")
  caliber_dir_path <- download_caliber_repo()

  # Set filepath constants --------------------------------------------------

  CALIBER_PRIMARY <- file.path(caliber_dir_path, "primary_care")
  CALIBER_SECONDARY <- file.path(caliber_dir_path, "secondary_care")
  CSV_REGEX <- "+\\.csv$"

  PRIMARY_CARE_FILES <- list.files(CALIBER_PRIMARY,
                                   pattern = CSV_REGEX)

  SECONDARY_CARE_FILES <- list.files(CALIBER_SECONDARY,
                                     pattern = CSV_REGEX)

  SECONDARY_CARE_FILES_ICD <- subset(SECONDARY_CARE_FILES, grepl("^ICD_", SECONDARY_CARE_FILES))
  SECONDARY_CARE_FILES_OPCS <- subset(SECONDARY_CARE_FILES, grepl("^OPCS_", SECONDARY_CARE_FILES))

  # Read CALIBER files and reformat -----------------------------------------

  # Note - currently removes medcodes and secondary descriptions for read codes

  # Read files into 3 dataframes - primary care and secondary care (ICD and OPCS)
  result <- c(
    "primary_care_codes_read2",
    "primary_care_codes_read3",
    "secondary_care_codes_icd10",
    "secondary_care_codes_icd9",
    "secondary_care_codes_opcs4"
  ) %>%
    purrr::set_names() %>%
    purrr::map(~ NULL)

  message("Reading caliber clinical codes lists into R and reformtting")
  result$primary_care_codes_read2 <-
    read_csv_to_named_list_and_combine(CALIBER_PRIMARY,
                                       filenames = PRIMARY_CARE_FILES,
                                       standardising_function = standardise_primary_care) %>%
    reformat_caliber_read2()

  result$secondary_care_codes_icd10 <-
    read_csv_to_named_list_and_combine(CALIBER_SECONDARY,
                                       filenames = SECONDARY_CARE_FILES_ICD,
                                       standardising_function = standardise_secondary_care_icd10) %>%
    reformat_caliber_icd10(all_lkps_maps = all_lkps_maps)

  result$secondary_care_codes_opcs4 <-
    read_csv_to_named_list_and_combine(CALIBER_SECONDARY,
                                       filenames = SECONDARY_CARE_FILES_OPCS,
                                       standardising_function = standardise_secondary_care_opcs4)


  # Map codes (read2 - read3,  icd10 - icd9) --------------------------------
  message("Mapping read2 codes to read3")
  result$primary_care_codes_read3 <-
    map_caliber_multiple_disease_categories(
      result$primary_care_codes_read2,
      all_lkps_maps = all_lkps_maps,
      from = "read2",
      to = "read3"
    )

  message("Mapping icd10 to icd9 codes")
  result$secondary_care_codes_icd9 <-
    map_caliber_multiple_disease_categories(
      result$secondary_care_codes_icd10,
      all_lkps_maps = all_lkps_maps,
      from = "icd10",
      to = "icd9"
    )

  message("Removing 'X' from ends of 3 character ICD10 codes")
  result$secondary_care_codes_icd10 <- codemapper:::strip_x_from_alt_icd10(df = result$secondary_care_codes_icd10,
                                                                           alt_icd10_code_col = "code")

  # combine
  message("Concatenating results")
  dplyr::bind_rows(result)
}

}
