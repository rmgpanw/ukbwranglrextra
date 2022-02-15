
mtcars2 <- mtcars %>%
  tibble::rownames_to_column(var = "make") %>%
  dplyr::mutate(
    "is_merc" = ifelse(
      stringr::str_detect(.data[["make"]],
                          pattern = "^Merc"),
      yes = 1,
      no = 0
    ),
    "is_mazda" = ifelse(
      stringr::str_detect(.data[["make"]],
                          pattern = "^Mazda"),
      yes = 1,
      no = 0
    )
  ) %>%
  tibble::as_tibble() %>%
  dplyr::mutate(dplyr::across(tidyselect::starts_with("is_"),
                              as.factor))

# `map_models()` ----------------------------------------------------------
test_that("`map_models()` returns expected results", {
  # run function
  result <- map_models(df = mtcars2,
                       model_str = "is_merc ~ mpg + cyl + disp + {.x}",
                       params = c(
                         # DRAT = 'drat',
                         #          WT = 'wt',
                         #          QSEC = 'qsec',
                         #          VS = 'vs',
                         #          AM = 'am',
                                  GEAR = 'gear',
                                  CARB = 'carb'),
                       engine = parsnip::set_engine(object = parsnip::logistic_reg(),
                                                    engine = "glm"),
                       rm_raw_model = FALSE)

  # expected colnames
  expect_equal(names(result),
               expected = c("wflow_id",
                            "info",
                            "option",
                            "result",
                            "wflow_label",
                            "fit"))

  # expected names in `fit` column
  expect_equal(names(result$fit[[1]]),
               expected = c("model_raw", "model_tidy", "model_glance"))

  # expected names in `fit` column, `model_raw$result`
  expect_equal(names(result$fit[[1]]$model_raw),
               expected = c("result", "error"))

})

test_that("`map_models()` `rm_raw_model` arg works as expected", {
  # run function
  result_raw_mod_rm <- map_models(df = mtcars2,
                       model_str = "is_merc ~ mpg + cyl + disp + {.x}",
                       params = c(
                         # DRAT = 'drat',
                         #          WT = 'wt',
                         #          QSEC = 'qsec',
                         #          VS = 'vs',
                         #          AM = 'am',
                         GEAR = 'gear',
                         CARB = 'carb'),
                       engine = parsnip::set_engine(object = parsnip::logistic_reg(),
                                                    engine = "glm"),
                       rm_raw_model = TRUE)
})
