vol_for_sequencing(
    sequencing_spreadsheet = "seq_spreadsheet_ex.csv",
    vol_for_sequencing  = 0.00001, # LITERS
    conc_for_sequencing = 5e-9, # MOLES PER LITER
    length_fragment = 176,  # BASE PAIRS (before library prep, including primers)
    colname_sample = "sample_id",
    colname_conc   = "conc_ng_ml", # NANOGRAMS PER MICROLITER
    colname_vol    = "vol_remain_ul", # MICROLITERS
    colname_lib    = "library_index",
    colname_primer = "primer_index", 
    conc_units = "ngml" # options: "ngml" "ngul"
)

# vol_for_sequencing(
    # sequencing_spreadsheet = "/Users/threeprime/Desktop/20160627_elwha_pcr2concentration_JO.csv",
    # vol_for_sequencing  = 0.00001, # LITERS
    # conc_for_sequencing = 5e-9, # MOLES PER LITER
    # length_fragment = 176,  # BASE PAIRS (before library prep, including primers)
    # colname_sample = "sample_id",
    # colname_conc   = "conc_ng_ml", # NANOGRAMS PER MICROLITER
    # colname_vol    = "vol_remain_ul", # MICROLITERS
    # colname_lib    = "Library",
    # colname_primer = "Tag",
    # conc_units = "ngml" # options: "ngml" "ngul"
# )
