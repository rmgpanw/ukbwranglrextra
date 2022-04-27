
# CONSTANTS ---------------------------------------------------------------



# TESTS -------------------------------------------------------------------

# `map_clinical_events_to_phecodes()` -------------------------------------


test_that("`map_clinical_events_to_phecodes()` returns expected output", {
  result <- map_clinical_events_to_phecodes(dummy_clinical_events,
                                            all_lkps_maps = "all_lkps_maps.db",
                                            min_date_only = FALSE)
})
