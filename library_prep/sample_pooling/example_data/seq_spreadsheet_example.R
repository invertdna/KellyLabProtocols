#!/usr/bin/env Rscript

N_samples   <- 26
N_libraries <- 4

sample_lib_mode <- N_samples %/% N_libraries
samples_per_lib <- rep(sample_lib_mode, N_libraries)
samples_per_lib[N_libraries] <- sample_lib_mode + (N_samples %% N_libraries)



random_index <- function(length, times)
{
  nucleotides <- c("A", "C", "T", "G")
  replicate(times, paste(
    sample(nucleotides,
           size = length,
           replace = TRUE),
   collapse = "")
  )
}

sample_id        <- LETTERS[1:N_samples]
vol_remain_ul <- rep(18, N_samples)
conc_ng_ml        <- rnorm(N_samples, mean = 200, sd = 50)
library_index    <- rep(random_index(8, N_libraries), times = samples_per_lib)
primer_index     <- random_index(6, N_samples)

seq_spreadsheet_ex <- data.frame(
  sample_id, vol_remain_ul, conc_ng_ml, library_index, primer_index
)

write.csv(seq_spreadsheet_ex, file = "seq_spreadsheet_ex.csv", row.names = FALSE)
