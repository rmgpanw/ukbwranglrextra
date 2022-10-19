#' Create a UKB data frame with a unique ID column
#'
#' A convenience function that returns a data frame for a main UK Biobank
#' dataset with a unique ID column, created by concatenating values from a
#' selection of variables. **Manual validation of any subsequent linkage is
#' strongly advised.**
#'
#' @inheritParams ukbwranglr::read_ukb
#' @inheritParams data.table::fread
#' @inheritParams create_unique_id
#'
#' @return A data frame
#' @export
#' @seealso \code{\link{create_unique_id}}
create_unique_id_df <- function(path,
                                delim = "\t",
                                ukb_data_dict = ukbwranglr::get_ukb_data_dict(),
                                ukb_codings = ukbwranglr::get_ukb_codings(),
                                descriptive_colnames = TRUE,
                                label = FALSE,
                                field_ids = c("31",
                                              "52",
                                              "34",
                                              "21000",
                                              "53",
                                              "96",
                                              "50"),
                                instances = "0",
                                id_col = "..unique_id",
                                remove = TRUE,
                                .ignore_duplicate_ids = FALSE) {
  # make data dictionary
  data_dict_full <- ukbwranglr::make_data_dict(path,
                                          ukb_data_dict = ukb_data_dict)

  # checks for link id, and create `data_dict` for `read_ukb()`
  data_dict_read_ukb <-
    create_unique_id_validation_checks(
      data_dict = data_dict_full,
      field_ids = field_ids,
      instances = instances,
      keep_eid_col = TRUE
    )

  # read selected field_ids
  ukb_main <- ukbwranglr::read_ukb(
    path = path,
    delim = delim,
    data_dict = data_dict_read_ukb,
    ukb_data_dict = ukb_data_dict,
    ukb_codings = ukb_codings,
    descriptive_colnames = descriptive_colnames,
    label = label
  )

  # append unique id
  ukb_main <- create_unique_id(
    ukb_main = ukb_main,
    ukb_data_dict = ukb_data_dict,
    field_ids = field_ids,
    instances = instances,
    id_col = id_col,
    remove = remove,
    .ignore_duplicate_ids = .ignore_duplicate_ids
  )
  # return result
  return(ukb_main)
}

#' Create a unique participant ID by combining UKB data field IDs
#'
#' Creates a unique participant ID by concatenating values from a selection of
#' UKB data fields. An error is raised if the final ID column contains
#' non-unique values. **Manual validation of any subsequent linkage is strongly
#' advised.**
#'
#' By default, the following field IDs are used: 31 (Sex), 52 (Month of birth),
#' 34 (Year of birth), 21000 (Ethnic background), 53 (Date of attending
#' assessment centre), 96 (Time since interview start at which blood pressure
#' screen(s) shown), and 50 (Standing height). Any columns of type factor will
#' be converted to type integer.
#'
#' @param ukb_main A data frame - a UKB main dataset.
#' @param ukb_data_dict The UK Biobank data dictionary (available
#'   \href{https://biobank.ctsu.ox.ac.uk/crystal/exinfo.cgi?src=accessing_data_guide}{data
#'    online}).
#' @param field_ids A character vector of fields IDs that will be used to create
#'   the new unique ID column. These should match the values under column
#'   'Field' in the UK Biobank data dictionary.
#' @param id_col Name of the the new column to be created.
#' @param remove If \code{TRUE}, remove input columns from output data frame.
#' @param .ignore_duplicate_ids If \code{TRUE}, allow duplicate ID values and
#'   raise a warning if any are found. May be helpful for debugging. By default
#'   this is \code{FALSE}.
#' @param instances A character vector of instances to include when generating
#'   the new unique ID column. Should contain one or more of the following
#'   digits: '0', '1', '2', '3'. Note that more recent datasets may include
#'   instances that are not present in older datasets. By default only the first
#'   instance is used.
#'
#' @return A data frame with an additional column named as specified by
#'   \code{id_col}.
#' @export
create_unique_id <- function(ukb_main,
                             ukb_data_dict = ukbwranglr::get_ukb_data_dict(),
                             field_ids = c("31",
                                           "52",
                                           "34",
                                           "21000",
                                           "53",
                                           "96",
                                           "50"),
                             instances = "0",
                             id_col = "..unique_id",
                             remove = TRUE,
                             .ignore_duplicate_ids = FALSE) {
  # create data dictionary
  data_dict <- ukb_main %>%
    ukbwranglr::make_data_dict(ukb_data_dict = ukb_data_dict)

  # checks - also filter `data_dict` for specified Field IDs and instances
  data_dict <- create_unique_id_validation_checks(data_dict = data_dict,
                                                  field_ids = field_ids,
                                                  instances = instances,
                                                  keep_eid_col = FALSE)

  # arrange `data_dict` and `ukb_main` by FieldID/instance/array
  data_dict$FieldID_instance_array <- paste(data_dict$FieldID,
                                            data_dict$instance,
                                            data_dict$array,
                                            sep = "_")

  data_dict <- data_dict %>%
    dplyr::arrange(.data[["FieldID_instance_array"]])

  ukb_main <- ukb_main %>%
    dplyr::select(tidyselect::everything(),
                  tidyselect::all_of(data_dict[["colheaders_raw"]])) %>%
    # convert factors to integers
    dplyr::mutate(dplyr::across(tidyselect:::where(is.factor),
                                as.integer))

  # create unique id col - unite values in selected `field_ids` columns
  ukb_main <- tidyr::unite(
    ukb_main,
    col = "..id_col",
    sep = "_",
    remove = remove,
    na.rm = TRUE,
    tidyselect::all_of(data_dict[["colheaders_raw"]])
  )

  # rename unique ID column, as specified
  ukb_main <- ukbwranglr:::rename_cols(ukb_main,
                                       old_colnames = "..id_col",
                                       new_colnames = id_col)

  # add variable label - concatenate columns used to make unique_id, separated
  # by ';'
  attributes(ukb_main[[id_col]])$label <-
    stringr::str_c(data_dict$FieldID_instance_array,
                   sep = "",
                   collapse = "; ")

  # check that unique id col contains only unique values
  n_distinct_ids <- dplyr::n_distinct(ukb_main[[id_col]])

  # if non-unique ids, raise either warning or error
  if (n_distinct_ids != nrow(ukb_main)) {
    msg <- paste0(
      "There are ",
      n_distinct_ids,
      " unique ID values but ",
      nrow(ukb_main),
      " participants in `ukb_main`. Perhaps try a different combination of `field_ids`."
    )

    ifelse(.ignore_duplicate_ids,
           yes = warning(msg),
           no = stop(msg))
  }

  # return result
  return(ukb_main)
}

# my_df <- data.frame(x = letters[1:3],
#                     y = letters[5:7],
#                     z = LETTERS[8:10])
#
# my_df %>%
#   tidyr::unite(
#     col = "unique_id",
#     sep = "_",
#     all_of(c("z", "x")),
#     remove = FALSE,
#     na.rm = FALSE
#   )


#' Validation helper for \code{create_unique_id}
#'
#' See code comments for checks.
#'
#' @inheritParams create_unique_id
#' @param keep_eid_col If `TRUE`, retain eid col in returned data dictionary
#'   data frame.
#'
#' @return Returns \code{data_dict}, filtered for only Field IDs in
#'   \code{field_ids} and instances in `instances`.
#' @noRd
create_unique_id_validation_checks <- function(data_dict,
                                               field_ids,
                                               instances,
                                               keep_eid_col = FALSE) {
  # check that `field_ids` does not include eid's provided by UK Biobank
  assertthat::assert_that(!"Participant identifier ('eid')" %in% field_ids,
                          msg = "Error! `field_ids` cannot include 'Participant identifier ('eid')'")

  # check that 'eid' column is present
  assertthat::assert_that("eid" %in% data_dict[["FieldID"]],
                          msg = "Error! An 'eid' column must be present in `ukb_main`")

  # keep eid row separately
  eid_row <- data_dict %>%
    dplyr::filter(FieldID == "eid")

  # check that all `field_ids` are present in dataset
  missing_fields_ids <-
    subset(field_ids,!field_ids %in% data_dict$FieldID)

  assertthat::assert_that(
    length(missing_fields_ids) == 0,
    msg = paste0(
      "Error! Some fields IDs are not present in `ukb_main`. Try using `ukbwranglr::make_data_dict` to check which fields IDs are available. Missing fields IDs: ",
      stringr::str_c(missing_fields_ids,
                     sep = "",
                     collapse = ", ")
    )
  )

  # filter data dictionary for these field_ids
  data_dict <- data_dict %>%
    dplyr::filter(.data[["FieldID"]] %in% !!field_ids)

  # check that `instances` are valid
  valid_instances_msg <- "Argument `instances` should be a character vector containing one or more of the following digits: '0', '1', '2', '3'."

  assertthat::assert_that(is.character(instances),
                          msg = valid_instances_msg)

  assertthat::assert_that(all(instances %in% c("0", "1", "2", "3")),
                          msg = valid_instances_msg)

  # filter data dictionary for these instances
  data_dict <- data_dict %>%
    dplyr::filter(.data[["instance"]] %in% !!instances)

  # add back eid row, if needed
  if (keep_eid_col) {
    data_dict <- eid_row %>%
      dplyr::bind_rows(data_dict)
  }

  return(data_dict)
}
