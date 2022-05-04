## code to prepare
## `inst/extdata/mapped_caliber_overlapping_disease_categories.csv` dataset goes
## here

library(codemapper)

# Reformat and map CALIBER repo (remove col_filters to get all possible
# overlaps)
caliber_raw <- read_caliber_raw(Sys.getenv("CALIBER_LOCAL"))

caliber_ukb <- reformat_caliber_for_ukb(
  caliber_raw,
  all_lkps_maps = NULL,
  col_filters = NULL,
  overlapping_disease_categories = "warning"
)

# get overlapping disease categories
overlapping_disease_categories <- caliber_ukb %>%
  dplyr::bind_rows() %>%
  ukbwranglr:::identify_overlapping_disease_categories() %>%
  .$clinical_codes %>%
  dplyr::arrange(disease,
                 code_type,
                 code,
                 category)

# write to csv in inst/extdata - then manually indicate which rows to keep
readr::write_csv(
  overlapping_disease_categories,
  file = file.path(
    "inst",
    "extdata",
    "caliber_mapped_overlapping_disease_categories.csv"
  )
)
