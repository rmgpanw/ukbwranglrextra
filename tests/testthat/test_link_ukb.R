
# SET UP ------------------------------------------------------------------

# TODO - need better dummy data for testing

## Dummy data and essentials -----------------------------------------------

# dummy_ukb_main <- ukbwranglr::dummy_main_dataset_clinical_events()
#
# ukb_data_dict <- ukbwranglr::get_ukb_data_dict()
# ukb_codings <- ukbwranglr::get_ukb_codings()
#
# ## Data dictionary and fields for creating unique id
# data_dict <- ukbwranglr::make_data_dict(dummy_ukb_main)
# field_ids <- c(
#   # "Participant identifier ('eid')" = "eid",
#   # "Cancer code, self-reported" = "20001",
#   "Non-cancer illness code, self-reported" = "20002",
#   # "Interpolated Year when cancer first diagnosed" = "20006",
#   "Interpolated Year when non-cancer illness first diagnosed" = "20008",
#   "Diagnoses - ICD10" = "41270",
#   # "Diagnoses - ICD9" = "41271",
#   "Date of first in-patient diagnosis - ICD10" = "41280",
#   # "Date of first in-patient diagnosis - ICD9" = "41281",
#   # "Underlying (primary) cause of death: ICD10" = "40001",
#   # "Contributory (secondary) causes of death: ICD10" = "40002",
#   # "Date of death" = "40000",
#   # "Treatment/medication code" = "20003",
#   # "Date of attending assessment centre" = "53",
#   # "Date of cancer diagnosis" = "40005",
#   # "Type of cancer: ICD10" = "40006",
#   # "Type of cancer: ICD9" = "40013",
#   # "Operative procedures - OPCS4" = "41272",
#   # "Date of first operative procedure - OPCS4" = "41282",
#   # "Operative procedures - OPCS3" = "41273",
#   # "Date of first operative procedure - OPCS3" = "41283",
#   "Operation code" = "20004",
#   "Interpolated Year when operation took place" = "20010"
# )

# TESTS -------------------------------------------------------------------

# create_unique_id_validation_checks(data_dict = data_dict,
#                                    field_ids = field_ids)
