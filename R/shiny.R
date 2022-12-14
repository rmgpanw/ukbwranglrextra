#' Interactively create a UK Biobank dataset
#'
#' Launches a shiny app to interactively create a UK Biobank dataset
#'
#' @inheritParams ukbwranglr::read_ukb
#' @inheritParams shiny::shinyApp
#' @param ... Additional parameters passed on to \code{\link[shiny]{shinyApp}}
#' @import shiny
#' @import shinyFiles
#'
#' @return NULL
#'
#' @export
#' @examples
#' \dontrun{
#'   runMakeDataset()
#' }
runMakeDataset <- function(ukb_data_dict = ukbwranglr::get_ukb_data_dict(),
                           ukb_codings = ukbwranglr::get_ukb_codings(),
                           options = list(launch.browser = TRUE),
                           ...) {
  # SETUP -------------------------------------------------------------------

  selected_data_dict_cols <- c(
    "Field",
    "FieldID",
    "instance",
    "array",
    "Path",
    "Category",
    "Participants",
    "Items",
    "Stability",
    "ValueType",
    "Units",
    "ItemType",
    "Strata",
    "Sexed",
    "Instances",
    "Array",
    "Coding",
    "Notes",
    "Link",
    "descriptive_colnames",
    "colheaders_raw",
    "colheaders_processed"
  )

  derived_variables_list <- list(
    "Date of birth" = list(
      derive_function = ukbwranglr::derive_dob,
      required_fids = ukbwranglr::derive_dob(.details_only = TRUE)$required_field_ids
    ),
    "Ethnicity, simplified" = list(
      derive_function = ukbwranglr::derive_ethnic_background_simplified,
      required_fids = ukbwranglr::derive_ethnic_background_simplified(.details_only = TRUE)$required_field_ids
    )
  )

  derived_variables_required_fids_df <- derived_variables_list %>%
    purrr::map("required_fids") %>%
    purrr::map(toString) %>%
    tibble::as_tibble() %>%
    t() %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "Derived variable") %>%
    dplyr::rename("Required FieldIDs" = .data[["V1"]])

  # UI ----------------------------------------------------------------------


  # Conditional UI elements -------------------------------------------------

  download_data_dict_tabs <- tabsetPanel(
    id = "download_data_dict_inputs",
    type = "hidden",
    tabPanel("empty"),
    tabPanel(
      "download_data_dict_inputs",
      tags$hr(),
      h4("Download data dictionary"),
      downloadButton("download_data_dict_full",
                     label = "All variables"),
      downloadButton("download_data_dict_selected",
                     label = "Selected variables")
    )
  )

  make_dataset_tabs <- tabsetPanel(
    id = "make_dataset_inputs",
    type = "hidden",
    tabPanel("empty"),
    tabPanel(
      "make_dataset_inputs",
      tags$hr(),
      h4("Create dataset"),
      numericInput(
        "nrows",
        "N rows to include",
        value = 510000,
        min = 0,
        max = 510000
      ),
      checkboxInput("descriptive_colnames", "Descriptive column names", value = TRUE),
      checkboxInput("labelled", "Labelled", value = TRUE),
      uiOutput("derived_variables_input"),
      selectInput(
        inputId = "summarise_numerical_variables",
        "Summarise numerical variables",
        choices = c("mean", "min", "max", "sum", "n_values"),
        multiple = TRUE
      ),
      textInput(
        "ukb_dataset_outname",
        "Dataset file name",
        value = paste0("my_ukb_dataset_", Sys.Date(), ".rds")
      ),
      downloadButton("ukb_dataset",
                     label = "Make dataset"),
    )
  )

  clinical_events_tabs <- tabsetPanel(
    id = "clinical_events_inputs",
    type = "hidden",
    tabPanel("empty"),
    tabPanel(
      "clinical_events_inputs",
      tags$hr(),
      h4("Clinical events"),
      numericInput(
        "nrows_clinical_events",
        "N rows to include",
        value = 510000,
        min = 0,
        max = 510000
      ),
      textInput(
        "clinical_events_outname",
        "Clinical events file name",
        value = paste0("clinical_events_", Sys.Date(), ".csv")
      ),
      downloadButton("download_clinical_events",
                     label = "Extract clinical events"),
    )
  )

  # Main UI -----------------------------------------------------------------

  ui <- fluidPage(# Application title
    titlePanel("Make a UK Biobank dataset"),

    # Sidebar
    sidebarLayout(
      sidebarPanel(
        width = 2,
        h4("Input file paths"),
        shinyFilesButton(
          "ukb_main_dataset_filepath",
          label = "UKB main dataset",
          title = "Please select a UK Biobank main dataset file",
          multiple = FALSE,
          viewtype = "detail"
        ),
        shinyFilesButton(
          "data_dict_preselected",
          label = "Data dictionary (pre-selected)",
          title = "Please select a UK Biobank main dataset data dictionary file",
          multiple = FALSE,
          viewtype = "detail"
        ),
        download_data_dict_tabs,
        make_dataset_tabs,
        clinical_events_tabs
      ),

      # Show the data dictionary
      mainPanel(
        tabsetPanel(
          type = "pills",
          tabPanel(
            "Interactive data dictionary",
            reactable::reactableOutput("data_dict")
          ),
          tabPanel(
            "Selected variables",
            reactable::reactableOutput("data_dict_selected")
          ),
          tabPanel(
            "Selected input filepaths",
            h4("UKB main dataset"),
            verbatimTextOutput("ukb_main_dataset_filepath"),
            h4("Data dictionary (pre-selected)"),
            verbatimTextOutput("data_dict_preselected_filepath"),
          ),
          tabPanel(
            "Required FieldIDs for derived variables",
            reactable::reactableOutput("derived_variables_required_fids")
          )
        ),
        width = 9
      )
    ))


  # SERVER ------------------------------------------------------------------

  server <- function(input, output, session) {
    # DERIVED VARIABLE REQUIRED FIDS REACTABLE --------------------------------

    output$derived_variables_required_fids <-
      reactable::renderReactable(
        reactable::reactable(
          derived_variables_required_fids_df,
          filterable = TRUE,
          showPageSizeOptions = TRUE,
          pageSizeOptions = c(10, 25, 50, 100, 200),
          paginationType = "jump"
        )
      )

    # SHINYFILES SETUP --------------------------------------------------------

    volumes <-
      c(Home = fs::path_home(),
        "R Installation" = R.home(),
        getVolumes()())

    # Main dataset ------------------------------------------------------------

    # main dataset file
    shinyFileChoose(input,
                    "ukb_main_dataset_filepath",
                    roots = volumes,
                    session = session)

    # print UKB main dataset file selection to console
    observe({
      cat("\ninput$ukb_main_dataset_filepath value:\n\n")
      print(input$ukb_main_dataset_filepath)
    })

    # print file selection to browser
    output$ukb_main_dataset_filepath <- renderPrint({
      if (is.integer(input$ukb_main_dataset_filepath)) {
        cat("Please select a UKB main dataset")
      } else {
        parseFilePaths(volumes, input$ukb_main_dataset_filepath)$datapath %>%
          stringr::str_split(" \\u2022 ") %>%
          .[[1]] %>%
          as.list() %>%
          do.call(file.path, .)
      }
    })

    # get filepath to main dataset
    ukb_main_dataset_filepath <- reactive({
      if (is.integer(input$ukb_main_dataset_filepath)) {
        cat("No files have been selected (shinyFileChoose)")
      } else {
        parseFilePaths(volumes, input$ukb_main_dataset_filepath)$datapath %>%
          stringr::str_split(" \\u2022 ") %>%
          .[[1]] %>%
          as.list() %>%
          do.call(file.path, .)
      }
    })


    # Data dictionary (pre-selected) ------------------------------------------

    # pre-selected data dictionary
    shinyFileChoose(input,
                    "data_dict_preselected",
                    roots = volumes,
                    session = session)

    # print selection to console
    observe({
      cat("\ninput$data_dict_preselected value:\n\n")
      print(input$data_dict_preselected)
    })

    # print file selection to browser
    output$data_dict_preselected_filepath <- renderPrint({
      if (is.integer(input$data_dict_preselected)) {
        cat("Please select a data dictionary (optional)")
      } else {
        parseFilePaths(volumes, input$data_dict_preselected)$datapath %>%
          stringr::str_split(" \\u2022 ") %>%
          .[[1]] %>%
          as.list() %>%
          do.call(file.path, .)
      }
    })

    # get filepath to selected UKB data dictionary file (of pre-selected variables)
    data_dict_preselected_filepath <- reactive({
      if (is.integer(input$data_dict_preselected)) {
        cat("No files have been selected (shinyFileChoose)")
      } else {
        parseFilePaths(volumes, input$data_dict_preselected)$datapath %>%
          stringr::str_split(" \\u2022 ") %>%
          .[[1]] %>%
          as.list() %>%
          do.call(file.path, .)
      }
    })

    # DATA DICTIONARY REACTABLE TABLE -----------------------------------------


    # Full data dictionary reactable ------------------------------------------


    # create data dict
    ukb_main_dataset_filepath_ext <- reactive({
      req(ukb_main_dataset_filepath())

      ext <- ukbwranglr:::extract_file_ext(ukb_main_dataset_filepath())

      if (!ext %in% c("csv", "tsv", "txt", "tab", "dta")) {
        validate(
          "Invalid file extension for main dataset; Please upload a comma (.csv), tab delimited (.tsv/.txt/.tab) or STATA (.dta) file"
        )
      }

      ext
    })

    data_dict <- reactive({
      req(ukb_main_dataset_filepath(),
          ukb_main_dataset_filepath_ext())

      switch(
        ukb_main_dataset_filepath_ext(),
        csv = ukbwranglr::make_data_dict(
          ukb_main_dataset_filepath(),
          delim = ",",
          ukb_data_dict = ukb_data_dict
        ) %>%
          dplyr::select(tidyselect::all_of(selected_data_dict_cols)),
        tsv = ukbwranglr::make_data_dict(
          ukb_main_dataset_filepath(),
          delim = "\t",
          ukb_data_dict = ukb_data_dict
        ) %>%
          dplyr::select(tidyselect::all_of(selected_data_dict_cols)),
        txt = ukbwranglr::make_data_dict(
          ukb_main_dataset_filepath(),
          delim = "\t",
          ukb_data_dict = ukb_data_dict
        ) %>%
          dplyr::select(tidyselect::all_of(selected_data_dict_cols)),
        tab = ukbwranglr::make_data_dict(
          ukb_main_dataset_filepath(),
          delim = "\t",
          ukb_data_dict = ukb_data_dict
        ) %>%
          dplyr::select(tidyselect::all_of(selected_data_dict_cols)),
        dta = ukbwranglr::make_data_dict(
          ukb_main_dataset_filepath(),
          delim = "\t",
          ukb_data_dict = ukb_data_dict
        ) %>%
          dplyr::select(tidyselect::all_of(selected_data_dict_cols))
      )
    })

    # read data dictionary of pre-selected variables (optional)
    data_dict_preselected <- reactive({
      req(data_dict_preselected_filepath())

      ext <- ukbwranglr:::extract_file_ext(data_dict_preselected_filepath())

      if (ext != "csv") {
        validate("Invalid file extension for data dicionary, must be a csv file")
      }

      readr::read_csv(data_dict_preselected_filepath())
    })

    # display data dictionary for selected UKB file
    output$data_dict <- reactable::renderReactable({
      req(data_dict())

      reactable::reactable(
        # data,
        data_dict(),
        filterable = TRUE,
        searchable = TRUE,
        showPageSizeOptions = TRUE,
        pageSizeOptions = c(10, 25, 50, 100, 200),
        groupBy = "Field",
        selection = "multiple",
        onClick = "select",
        defaultSelected = 1,
        paginationType = "jump",
        theme = reactable::reactableTheme(
          rowSelectedStyle = list(backgroundColor = "#eee", boxShadow = "inset 2px 0 0 0 #ffa62d")
        )
      )

    })

    # selected rows
    selected <-
      reactive(reactable::getReactableState("data_dict", "selected"))


    # Update selected variables from pre-selected data dict upload ------------

    # update selected rows to include those in data_dict_preselected
    observeEvent(data_dict_preselected(), {
      # currently selected
      currently_selected_rowids <-
        reactable::getReactableState("data_dict", "selected")

      # get row indices to update
      preselected_rowids <- data_dict() %>%
        tibble::rowid_to_column() %>%
        dplyr::filter(.data[["colheaders_raw"]] %in% data_dict_preselected()$colheaders_raw) %>%
        dplyr::pull(.data[["rowid"]])

      # union
      updated_selection_rowids <-
        unique(c(currently_selected_rowids, preselected_rowids))

      # update
      reactable::updateReactable("data_dict", selected = updated_selection_rowids)
    })


    # Selected variables reactable --------------------------------------------

    # print selected variables as reactable
    output$data_dict_selected <- reactable::renderReactable({
      req(data_dict())
      req(selected())

      reactable::reactable(
        # data,
        data_dict()[selected(),],
        filterable = TRUE,
        searchable = TRUE,
        showPageSizeOptions = TRUE,
        pageSizeOptions = c(10, 25, 50, 100, 200),
        groupBy = "Field",
        paginationType = "jump"
      )

    })

    # print selected rows in console
    observe({
      print(data_dict()[selected(),])
    })


    # DYNAMIC UI ELEMENTS -----------------------------------------------------

    # Download and make dataset tabs ------------------------------------------

    observeEvent(data_dict(), {
      if (!is.null(data_dict())) {
        updateTabsetPanel(inputId = "download_data_dict_inputs", selected = "download_data_dict_inputs")
        updateTabsetPanel(inputId = "make_dataset_inputs", selected = "make_dataset_inputs")
        updateTabsetPanel(inputId = "clinical_events_inputs", selected = "clinical_events_inputs")
      } else {
        updateTabsetPanel(inputId = "download_data_dict_inputs", selected = "empty")
        updateTabsetPanel(inputId = "make_dataset_inputs", selected = "empty")
        updateTabsetPanel(inputId = "clinical_events_inputs", selected = "empty")
      }
    })

    # Derived variables -------------------------------------------------------

    derived_variables_list_available <- reactive({
      req(selected())
      selected_fieldids <- data_dict()[selected(),]$FieldID

      derived_variables_list %>%
        purrr::keep( ~ all(.x$required_fids %in% selected_fieldids))
    })

    output$derived_variables_input <- renderUI({
      selectInput(
        inputId = "derived_variables_input",
        label = "Select variables to derive",
        multiple = TRUE,
        choices = names(derived_variables_list_available())
      )
    })

    # DOWNLOADS ---------------------------------------------------------------

    # Data dictionary - full --------------------------------------------------
    output$download_data_dict_full <- downloadHandler(
      filename = function() {
        ("data_dict_full.csv")
      },
      content = function(file) {
        readr::write_csv(data_dict(), file, na = "")
      }
    )


    # Data dictionary - selected rows only ------------------------------------
    output$download_data_dict_selected <- downloadHandler(
      filename = function() {
        ("data_dict_selected.csv")
      },
      content = function(file) {
        readr::write_csv(data_dict()[selected(),], file, na = "")
      }
    )


    # UKB dataset -------------------------------------------------------------
    output$ukb_dataset <- downloadHandler(
      filename = function() {
        input$ukb_dataset_outname
      },
      content = function(file) {
        # setup
        ext <- ukbwranglr:::extract_file_ext(file)

        if (!ext %in% c("rds", "dta", "csv", "tsv")) {
          showNotification(
            "Invalid file type. File extensionn must be one of '.rds', '.dta', '.csv' or '.tsv'. Please click 'cancel' in dialog box",
            type = "error"
          )
          stop(
            "Invalid file type. File extensionn must be one of '.rds', '.dta', '.csv' or '.tsv'"
          )
        }

        req(data_dict())

        delim <- switch(
          ukb_main_dataset_filepath_ext(),
          csv = ",",
          tsv = "\t",
          txt = "\t",
          tab = "\t",
          dta = ""
        )

        # read selected variables
        notify <- function(msg, id = NULL) {
          showNotification(msg,
                           id = id,
                           duration = NULL,
                           closeButton = FALSE)
        }

        message("Reading UKB data")
        id <- notify("Reading data...")
        on.exit(removeNotification(id), add = TRUE)

        ukb_main <-
          ukbwranglr::read_ukb(
            path = ukb_main_dataset_filepath(),
            delim = delim,
            ukb_data_dict = ukb_data_dict,
            ukb_codings = ukb_codings,
            data_dict = data_dict()[selected(),],
            descriptive_colnames = input$descriptive_colnames,
            label = input$labelled,
            nrows = input$nrows
          )

        # TODO - get ethnic background working
        # derive variables
        if (length(input$derived_variables_input > 0)) {
          notify("Deriving variables...", id = id)
          for (variable in input$derived_variables_input) {
            ukb_main <-
              derived_variables_list_available()[[variable]]$derive_function(ukb_main,
                                                                             ukb_data_dict = ukb_data_dict)
          }
        }

        # summarise numerical variables
        if (length(input$summarise_numerical_variables > 0)) {
          notify("Summarising numerical variables...", id = id)
          for (summary_function in input$summarise_numerical_variables) {
            try(ukb_main <- ukbwranglr::summarise_numerical_variables(ukb_main = ukb_main,
                                                      summary_function = summary_function))
          }
        }

        if (ext == "dta") {
          # STATA colheaders cannot contain '.' and must not be too long
          message("STATA output: renaming columns...\n")
          notify("STATA output: renaming columns...", id = id)
          names(ukb_main) <- ukbwranglr::make_data_dict(ukb_main,
                                            ukb_data_dict = ukb_data_dict)$colheaders_processed
          # names(ukb_main) <- data_dict()[selected(), ]$colheaders_processed
        }

        # write dataset
        message("Writing data\n")
        notify("Writing data...", id = id)

        switch(
          ext,
          rds = saveRDS(ukb_main, file = file),
          dta = haven::write_dta(ukb_main, path = file),
          csv = readr::write_csv(ukb_main, file = file, na = ""),
          tsv = readr::write_tsv(ukb_main, file = file, na = "")
        )

        message("Finished writing data")
        showNotification("Dataset succesfully created!", duration = 10)
      }
    )

    # Clinical events ---------------------------------------------------------
    output$download_clinical_events <- downloadHandler(
      filename = function() {
        input$clinical_events_outname
      },
      content = function(file) {
        # TODO - extract this to function (currently c+p from ukb dataset section)

        # setup
        ext <- ukbwranglr:::extract_file_ext(file)

        if (!ext %in% c("rds", "dta", "csv", "tsv")) {
          showNotification(
            "Invalid file type. File extensionn must be one of '.rds', '.dta', '.csv' or '.tsv'. Please click 'cancel' in dialog box",
            type = "error"
          )
          stop(
            "Invalid file type. File extensionn must be one of '.rds', '.dta', '.csv' or '.tsv'"
          )
        }

        req(data_dict())

        delim <- switch(
          ukb_main_dataset_filepath_ext(),
          csv = ",",
          tsv = "\t",
          txt = "\t",
          tab = "\t",
          dta = ""
        )

        # get available clinical events FIDs
        available_clinical_events <-
          ukbwranglr:::CLINICAL_EVENTS_FIELD_IDS %>%
          purrr::keep(~ all(.x %in% data_dict()$FieldID))

        available_clinical_events_fids <- available_clinical_events %>%
          purrr::reduce(c)

        message(paste0(
          "Available clinical events: ",
          stringr::str_c(
            names(available_clinical_events),
            sep = "",
            collapse = ", "
          )
        ))

        clinical_events_data_dict <- data_dict() %>%
          dplyr::filter(.data[["FieldID"]] %in% c("eid",
                                                  available_clinical_events_fids))

        # read clinical events variables
        notify <- function(msg, id = NULL) {
          showNotification(msg,
                           id = id,
                           duration = NULL,
                           closeButton = FALSE)
        }

        message("Reading UKB data")
        id <- notify("Reading clinical events data...")
        on.exit(removeNotification(id), add = TRUE)

        result <-
          ukbwranglr::read_ukb(
            path = ukb_main_dataset_filepath(),
            delim = delim,
            ukb_data_dict = ukb_data_dict,
            ukb_codings = ukb_codings,
            data_dict = clinical_events_data_dict,
            descriptive_colnames = TRUE,
            label = FALSE,
            nrows = input$nrows_clinical_events
          )

        # tidy
        # write dataset
        message("Tidying clinical events data\n")
        notify("Tidying clinical events...", id = id)
        result <- ukbwranglr::tidy_clinical_events(
          ukb_main = result,
          ukb_data_dict = ukb_data_dict,
          ukb_codings = ukb_codings,
          strict = FALSE
        ) %>%
          dplyr::bind_rows()

        # write dataset
        message("Writing data\n")
        notify("Writing data...", id = id)

        readr::write_csv(result, file = file, na = "")

        message("Finished writing data")
        showNotification("Dataset succesfully created!", duration = 10)
      }
    )
    # RUN APPLICATION ---------------------------------------------------------
  }

  shinyApp(ui = ui,
           server = server,
           options = options,
           ...)
}
