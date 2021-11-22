

# CONSTANTS ---------------------------------------------------------------

UKB_SYNTHETIC_DATASET_DETAILS <- tibble::tribble(
  ~ md5, ~ equivalent_dataset,
  'tabular.md5', "ukb_main",
  'medrec.md5', "gp_clinical",
  'genotype.md5', "",
  'rand_chr.md5', "",
  'bulk.md5',  ""
)

# PRIVATE -----------------------------------------------------------------

get_ukb_synthetic_data <- function(dataset = "tabular.md5",
                                   nrows = 10) {
  # validate args
  match.arg(dataset,
            choices = UKB_SYNTHETIC_DATASET_DETAILS$md5)

  start_time <- proc.time()

  # make index df
  ukb_dummy_index_df <- make_ukb_dummy_index_df()
  ukb_dummy_index_df <- split(ukb_dummy_index_df, ukb_dummy_index_df$md5_file)

  # read files
  result <- fread_ukb_dummy_files(urls = ukb_dummy_index_df[[dataset]]$url,
                                  fread_fn = purrr::partial(data.table::fread,
                                                            nrows = nrows))
  ukbwranglr:::time_taken_message(start_time)

  return(result)
}

make_ukb_dummy_index_df <- function() {
  # md5 files
  md5_files <- get_ukb_dummy_html_links_matching("md5$") %>%
    fread_ukb_dummy_files(purrr::partial(data.table::fread, header = FALSE)) %>%
    .$result %>%
    dplyr::bind_rows(.id = "md5")

  names(md5_files) <- c("md5_file", "md5", "file_name")

  # all urls
  all_urls <- tibble::tibble(url = get_all_ukb_dummy_html_links()) %>%
    dplyr::mutate("file_name" = installr::file.name.from.url(.data[["url"]]))

  # combine
  result <- md5_files %>%
    dplyr::left_join(all_urls,
                     by = "file_name")

  return(result)
}

fread_ukb_dummy_files <- function(urls,
                                  fread_fn = data.table::fread) {
  # make safely version of fread_fn
  fread_fn <- purrr::safely(fread_fn)

  # urls can be a named vector. If no names, use file from url link
  names(urls) <- installr::file.name.from.url(urls)

  # try reading urls
  result <- urls %>%
    purrr::map(fread_fn) %>%
    purrr::transpose() %>%
    purrr::map(purrr::compact)

  # warning if any errors
  if (length(result$error) > 0) {
    warning(paste0("Warning! Unable to get files at the following URL(s): ",
                   stringr::str_c(names(result$error),
                                  sep = "",
                                  collapse = ", ")))
  }

  # return results
  return(result)
}

get_ukb_dummy_html_links_matching <- function(pattern = "md5$") {
  # get a subset of urls on the ukb synthetic dataset web page that match a
  # pattern
  result <- get_all_ukb_dummy_html_links()

  result <- subset(result, stringr::str_detect(result,
                                               pattern))


  return(result)
}

get_all_ukb_dummy_html_links <- function() {
  # get all urls on the ukb synthetic dataset web page
  rvest::read_html(
    "https://biobank.ndph.ox.ac.uk/ukb/exinfo.cgi?src=UKB_Synthetic_Dataset.html"
  ) %>%
    rvest::html_elements("a") %>%
    rvest::html_attr("href")
}
