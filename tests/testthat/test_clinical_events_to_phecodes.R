
# CONSTANTS ---------------------------------------------------------------

## Dummy all_lkps_maps --------------------------------------------------------------

all_lkps_maps <- list(
  read_ctv3_icd10 = tibble::tribble(
    ~.rowid, ~read_code, ~icd10_code, ~mapping_status, ~refine_flag, ~add_code_flag, ~element_num, ~block_num, ~icd10_dagger_asterisk,
    75774L,    "XE0Uc",      "I10X",             "D",          "C",            "C",          "0",        "0",                     NA
  )


)

## Dummy data --------------------------------------------------------------

dummy_clinical_events <- tibble::tribble(
  ~eid,    ~source,    ~index,    ~code,    ~date,
  # HTN
  1,       "f41270",  "0_0",    "I10",      "2000-00-00",
  2,       "gpc3_r3",  "0_0",    "XE0Uc",      "2000-00-00",
)

expected_clinical_events_phecodes <- tibble::tribble(
  ~eid,    ~source,    ~index,    ~code,    ~date, ~icd10, ~phecode,
  # HTN
  1,       "f41270",  "0_0",    NA,      "2000-00-00", "I10", "401.1",
  2,       "gpc3_r3",  "0_0",    "XE0Uc",      "2000-00-00", "I10", "401.1"
)



# TESTS -------------------------------------------------------------------

# `map_clinical_events_to_phecodes()` -------------------------------------


test_that("`map_clinical_events_to_phecodes()` returns expected output", {
  result <- map_clinical_events_to_phecodes(clinical_events,
                                            all_lkps_maps = "all_lkps_maps.db",
                                            min_date_only = FALSE)
})
