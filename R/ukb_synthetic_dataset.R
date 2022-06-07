

# CONSTANTS ---------------------------------------------------------------

UKB_SYNTHETIC_DATASET_DETAILS <- tibble::tribble(
  ~ md5, ~ equivalent_dataset,
  'tabular.md5', "ukb_main",
  'medrec.md5', "gp_clinical",
  'genotype.md5', "",
  'rand_chr.md5', "",
  'bulk.md5',  ""
)


# PUBLIC ------------------------------------------------------------------

#' Get URL links to from UKB synthetic dataset webpage
#'
#' @param pattern Optionally filter the list of URLs for a subset matching a
#'   pattern.
#'
#' @return A character vector of URLs.
#' @export
#' @examples
#' \dontrun{
#' get_ukb_synthetic_dataset_html_links()
#'}
get_ukb_synthetic_dataset_html_links <- function(pattern = "md5$") {
  # get all urls on the ukb synthetic dataset web page
  result <- rvest::read_html(
    "https://biobank.ndph.ox.ac.uk/ukb/exinfo.cgi?src=UKB_Synthetic_Dataset.html"
  ) %>%
    rvest::html_elements("a") %>%
    rvest::html_attr("href")

  # get a subset of urls on the ukb synthetic dataset web page that match a
  # pattern
  if (!is.null(pattern)) {
    result <- subset(result, stringr::str_detect(result,
                                                 pattern))
  }

  return(result)
}

#' Download the UKB synthetic dataset
#'
#' @param dataset Character. The name of one of the `.md5` files on the [UKB
#'   synthetic data
#'   webpage](https://biobank.ndph.ox.ac.uk/ukb/exinfo.cgi?src=UKB_Synthetic_Dataset.html).
#'    All
#' @param download_dir Directory to download files to.
#' @param timeout Integer. Passed to `options()` to determine how long
#'   `download.file()` attempts to download for.
#'
#' @return Invisibly returns a list with downloaded file paths and errors for
#'   failed downloads.
#' @export
#'
#' @examples
#' \dontrun{
#'  get_ukb_synthetic_data(download_dir = tempdir())
#' }
get_ukb_synthetic_data <- function(
  download_dir,
  dataset = "tabular.md5",
  timeout = 600) {
  # validate args
  match.arg(dataset,
            choices = UKB_SYNTHETIC_DATASET_DETAILS$md5)

  if (!fs::dir_exists(download_dir)) {
    stop("`download_dir` either does not exist or is not a directory. Value supplied:",
         download_dir)
  }

  start_time <- proc.time()

  # make index df
  ukb_dummy_index_df <- make_ukb_dummy_index_df()

  # add download file paths
  ukb_dummy_index_df$path <- file.path(download_dir,
                                                ukb_dummy_index_df$file_name)

  # add download timeout option
  ukb_dummy_index_df$timeout <- timeout

  # split by dataset
  ukb_dummy_index_df <- split(ukb_dummy_index_df, ukb_dummy_index_df$md5_file)

  # download files
  download_file_safely <- purrr::safely(ukbwranglr:::download_file)

  result <- ukb_dummy_index_df[[dataset]][c('download_url', 'path')] %>%
    furrr::future_pmap(ukbwranglr:::download_file)

  ukbwranglr:::time_taken_message(start_time)

  invisible(result)
}

# PRIVATE -----------------------------------------------------------------

make_ukb_dummy_index_df <- function() {
  # md5 files
  md5_files <- get_ukb_synthetic_dataset_html_links(pattern = "md5$") %>%
    fread_ukb_dummy_files(purrr::partial(data.table::fread, header = FALSE)) %>%
    .$result %>%
    dplyr::bind_rows(.id = "md5")

  names(md5_files) <- c("md5_file", "md5", "file_name")

  # tempfix - files listed in `rand_chr.md5` should be suffixed with '.gz'
  md5_files <- md5_files %>%
    dplyr::mutate("file_name" = ifelse(stringr::str_detect(.data[["file_name"]],
                                                    "\\.dat$"),
                                yes = paste0(.data[["file_name"]], ".gz"),
                                no = .data[["file_name"]]))

  # all urls
  all_urls <- tibble::tibble(download_url = get_ukb_synthetic_dataset_html_links(pattern = NULL)) %>%
    dplyr::mutate("file_name" = installr::file.name.from.url(.data[["download_url"]]))

  # combine
  result <- md5_files %>%
    dplyr::left_join(all_urls,
                     by = "file_name") %>%
    tibble::as_tibble()

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
