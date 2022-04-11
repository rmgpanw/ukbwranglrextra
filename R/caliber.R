
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
