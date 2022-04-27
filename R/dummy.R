## Dummy all_lkps_maps --------------------------------------------------------------

ALL_LKPS_MAPS <- list(
  read_ctv3_icd10 = tibble::tribble(
    ~.rowid, ~read_code, ~icd10_code, ~mapping_status, ~refine_flag, ~add_code_flag, ~element_num, ~block_num, ~icd10_dagger_asterisk,
    75774L,    "XE0Uc",      "I10X",             "D",          "C",            "C",          "0",        "0",                     NA
  ),

  icd10_phecode = tibble::tribble(
    ~ICD10, ~PHECODE, ~Exl..Phecodes, ~Excl..Phenotypes, ~ALT_CODE,
    "E10",  "250.1",   "249-250.99",        "DIABETES",     "E10",
    "I10",  "401.1",   "401-405.99", "HYPERTENSIVE DISEASE",     "I10",
  ),

  read_v2_read_ctv3 = tibble::tribble(
    ~.rowid, ~CHAPTER, ~READV2_CODE,               ~READV2_DESC,                          ~TERMV2_DESC, ~TERMV2_ORDER, ~TERMV2_TYPE, ~READV3_CODE, ~TERMV3_CODE, ~TERMV3_TYPE,                          ~TERMV3_DESC, ~IS_ASSURED,
    60486L,      "C",      "C10E.", "Type 1 diabetes mellitus",            "Type 1 diabetes mellitus",          "00",          "P",      "X40J4",      "Y41PV",          "S",            "Type 1 diabetes mellitus",         "1",
    60487L,      "C",      "C10E.", "Type 1 diabetes mellitus", "Insulin dependent diabetes mellitus",          "12",          "S",      "X40J4",      "Yadjz",          "S", "Insulin-dependent diabetes mellitus",         "1",
    60489L,      "C",      "C10E.", "Type 1 diabetes mellitus",            "Type I diabetes mellitus",          "11",          "S",      "X40J4",      "Yagv5",          "P",            "Type I diabetes mellitus",         "1"
  )
)

## Dummy data --------------------------------------------------------------

DUMMY_CLINICAL_EVENTS <- tibble::tribble(
  ~eid,    ~source,    ~index,    ~code,    ~date,
  # HTN
  1,       "f41270",  "0_0",    "I10",      "2000-00-00",
  2,       "gpc3_r3",  "0_0",    "XE0Uc",      "2000-00-00",
)

EXPECTED_CLINICAL_EVENTS_PHECODES <- tibble::tribble(
  # HTN
  ~eid,   ~source, ~index,   ~code,        ~date, ~icd10, ~phecode,
  2, "gpc3_r3",  "0_0", "XE0Uc", "2000-00-00",  "I10",  "401.1",
  1,  "f41270",  "0_0",      NA, "2000-00-00",  "I10",  "401.1"
)
