#!/usr/bin/env Rscript

############################################################################################
# CALCULATION FOR THE NUMBER OF MOLECULES IN A SAMPLE OF DNA
############################################################################################

# how long is your fragment in base pairs?
length_fragment <- 200

# how much mass do you have, in nanograms (ng)?
mass_fragment <- 10



############################################################################################
# DON'T CHANGE THIS STUFF
############################################################################################



num_copies <- function(length_fragment, mass_fragment){
	# average mass of a single base pair (in Daltons, which is equivalent to g/mol aka molar mass)
	# put another way - How much does one mole of base pairs weigh in grams? (650 for double stranded DNA [dsDNA])
	mass_bp <- 650 #(== 1.0794*10^-12 ng)...  (a Dalton is 1.67*10^24 g)

	# Avogadro's number
	avo_num <- 6.022e23
	
	# unit correction (use 1e9 for nanograms)
	unit_correction <- 1e9
	
	copies <- (mass_fragment * avo_num)/(length_fragment * unit_correction * mass_bp)
	return(copies)
	
}

num_copies(length_fragment = 670, mass_fragment = 150)
