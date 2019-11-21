#!/usr/bin/env Rscript

vol_for_sequencing <- function(
  sequencing_spreadsheet,
  length_fragment, # BASE PAIRS
  vol_for_sequencing = 0.00005, # LITERS
  conc_for_sequencing = 5e-9, # MOLES PER LITER
  colname_sample = "sample_id",
  colname_conc   = "conc_ng_ul", # NANOGRAMS PER MICROLITER
  colname_vol    = "vol_remaining_ul", # MICROLITERS
  colname_lib    = "library_index",
  colname_primer = "primer_index",
  conc_units = "ngul" # options: "ngml" "ngul"
  ) {

  ##############################################################################
  # CONSTANTS
  ##############################################################################

  # get this from the core facility
  # vol_for_sequencing  <- 0.00001 # LITERS
  # conc_for_sequencing <- 5e-9 # MOLES PER LITER

  # This is the size of the fragments
  # length_fragment <- 176 # BASE PAIRS (before library prep, including primers)

  # this is the mean grams per mole of 1 base pair double stranded DNA
  # MAKE SURE THIS IS RIGHT!!!
  molar_mass_DNA_bp <- 650 # GRAMS PER MOLE

  # get this from the library prep kit
  # UNUSED FOR NOW
  # minimum_mass_for_library_prep <- 1e-10 # GRAMS (1e-10 = 100 picograms)

  # read in the spreadsheet
  spreadsheet  <- as.data.frame(gs_read(gs_title(sequencing_spreadsheet)))

  # detect number of libraries
  spreadsheet[,colname_lib]<-as.factor(spreadsheet[,colname_lib])
  nlib<-nlevels(spreadsheet[,colname_lib])
  samples.per.lib<-aggregate(spreadsheet[,colname_vol],by=list(spreadsheet[,colname_lib]), table)#Creates a table of number of samples per lib
  b<-match (spreadsheet[,colname_lib], samples.per.lib[,1])
  #Assign number of samples per lib
  ifelse(nlib==1,
        N_samples<- rep(nrow(spreadsheet),nrow(spreadsheet)),
        N_samples<-samples.per.lib[b,2])

  #spreadsheet[,colname_lib]<-as.factor(spreadsheet[,colname_lib])
  #N_samples_per_lib<-rowSums(table(spreadsheet[,colname_lib],spreadsheet[,colname_primer]))#Number of samples per library


  # give the names of the columns of interest
    # colname_sample <- "sample_id"
  # colname_conc   <- "conc_ng_ul" # NANOGRAMS PER MICROLITER
      # colname_vol    <- "vol_remaining_ul" # MICROLITERS
  # colname_lib    <- "library_index"
  # colname_primer <- "primer_index"

  if( conc_units == "ngml" ) {
    conc_conversion <- (1e-9 / 1e-3)
  } else if( conc_units == "ngul" ) {
    conc_conversion <- (1e-9 / 1e-6)
  } else {
    stop("conc_units must be either 'ngml' or 'ngul'")
  }

  sample_concs  <- spreadsheet[,colname_conc]
  sample_concs  <- sample_concs * conc_conversion # CONVERT TO GRAMS PER LITER
  vol_remaining <- spreadsheet[,colname_vol]
  vol_remaining <- spreadsheet[,colname_vol] * 1e-6 # CONVERT TO LITERS
  sample_lib    <- spreadsheet[,colname_lib]
  sample_id     <- spreadsheet[,colname_sample]

  ##############################################################################
  # CALCULATIONS
  ##############################################################################



  moles_for_seq <- conc_for_sequencing * vol_for_sequencing

  moles_per_sample <- moles_for_seq / N_samples
  
  molar_mass_library <- molar_mass_DNA_bp * length_fragment

  vol_per_sample <- molar_mass_library * moles_per_sample / sample_concs

  vol_per_sample_uL <- vol_per_sample * 1e6 #

  volume_warning <- vol_per_sample >= vol_remaining

  volume_warning[is.na(volume_warning)] <- FALSE

  notes <- ifelse(volume_warning, "volume needed exceeds volume remaining", "NA")

  # warning(paste("quadruple check the molar mass of DNA! Currently:", molar_mass_DNA_bp, "grams per mole"))

  output_df <- data.frame(
    sample_id,
    sample_lib,
    vol_per_sample_uL,
    notes,
    stringsAsFactors = FALSE)

  file_base <- sapply(strsplit(sequencing_spreadsheet,".", fixed = TRUE), function(x) paste(x[1:(length(x)-1)], collapse="."))
  output_file <- paste(file_base, "_vol_for_pool.csv", sep = "")
  write.csv(
    x = output_df,
    file = output_file,
    row.names = FALSE
  )
return(output_df)
}

# example:
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
