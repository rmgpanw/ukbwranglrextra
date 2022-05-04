
# CONSTANTS ---------------------------------------------------------------



# TESTS -------------------------------------------------------------------

# `map_clinical_events_to_phecodes()` -------------------------------------

# check:

## - Read code for hypertension maps to 3 character 'I10' icd10, which is
## undivided and therefore ends with 'X' - the 'X' should be removed from icd10
## column in final result

## - 5 character ICD10 codes e.g. M00.9 (Pyogenic arthritis, unspecified) in
## ICD10_CODE format should include in ALT_CODE format M009 + 10 other codes
## starting with M009 (the children specify the location of the arthritis)
test_that("`map_clinical_events_to_phecodes()` returns expected output", {
  result <- map_clinical_events_to_phecodes(DUMMY_CLINICAL_EVENTS,
                                            all_lkps_maps = ALL_LKPS_MAPS,
                                            min_date_only = FALSE)

  expect_equal(
    result,
    EXPECTED_CLINICAL_EVENTS_PHECODES
  )
})
