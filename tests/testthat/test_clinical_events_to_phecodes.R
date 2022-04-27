
# CONSTANTS ---------------------------------------------------------------



# TESTS -------------------------------------------------------------------

# `map_clinical_events_to_phecodes()` -------------------------------------

# check:

## - Read code for hypertension maps to 3 character 'I10' icd10, which is
## undivided and therefore ends with 'X' - the 'X' should be removed from icd10
## column in final result

## - sdfsf
test_that("`map_clinical_events_to_phecodes()` returns expected output", {
  result <- map_clinical_events_to_phecodes2(DUMMY_CLINICAL_EVENTS,
                                            all_lkps_maps = ALL_LKPS_MAPS,
                                            min_date_only = FALSE)

  expect_equal(
    result,
    EXPECTED_CLINICAL_EVENTS_PHECODES
  )
})
