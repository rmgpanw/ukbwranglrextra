
# SETUP -------------------------------------------------------------------

# caliber_raw <- read_caliber_raw(dummy_caliber_dir_path())
#
# all_lkps_maps <- codemapper::build
#
# caliber_ukb <- reformat_caliber_for_ukb(caliber_raw,
#                                         all_lkps_maps = ,
#                                         overlapping_disease_categories_csv = default_overlapping_disease_categories_csv())

# TESTS -------------------------------------------------------------------



# Check codes exist -------------------------------------------------------

## e.g. 4+ character ICD10 codes in ALT code format (majority of codes are in ICD10 code format)

# Undivided 3 character ICD10 codes (e.g. A38, Scarlet fever) -------------



# 3 character ICD10 codes that need expanding -----------------------------

# For CALIBER, where only the parent 3 characters are shown then it implies that
# all the children codes fall under the same disease category. These need
# expanding to include all children codes.

## ICD10 codes with modifiers --------------------------

# (e.g. E10 (has modifier 4), M90 (has modifier 5), F48 (no modifier 4/5, but does have 4 character children))



## Check 3 character ICD10 codes without modifiers are expanded (e.g.  D25, leiomyoma) --------



# Codes that are duplicated within (and between?) diseases ----------------


