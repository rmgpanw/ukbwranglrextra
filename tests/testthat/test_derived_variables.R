
# `append_rowwise_summary_indicator()` ------------------------------------

test_that("`append_rowwise_summary_indicator()` returns expected results", {
  df <- tibble::tribble(
      ~med_1, ~med_2, ~med_3, ~IGNORED_COL,
      "Beta blocker", NA, NA, "Beta blocker",
      "ACE inhibitor", "Beta blocker", NA, "Beta blocker",
      "Aspirin", "Aspirin", NA, "Beta blocker",
      "Aspirin", "ACE inhibitor", NA, "Beta blocker"
      )

  result <- append_rowwise_summary_indicator(df = df,
    cols = c("med_1", "med_2", "med_3"),
    yes_values = c("Beta blocker", "ACE inhibitor"),
    new_indicator_colname = "on_anti_hypertensive",
    new_indicator_var_label = "Taking anti-hypertensive medication",
    new_indicator_yes_value = "Yes",
    new_indicator_no_value = "No")

  expect_equal(as.character(result$on_anti_hypertensive),
               c("Yes", "Yes", "No", "Yes"))

  expect_equal(attributes(result$on_anti_hypertensive)$label,
               "Taking anti-hypertensive medication")
})
