library(tidyverse)

encoded <- data # create copy of your data to modify

for (col in colnames(data)) {
  if (!col %in% c("SEQID", "Medium", "Genus", "Strain", "Strain_full", "Phylum")) { # the columns in my dataset that were not ESKAPE columns with bioactivity annotations
    l <- numeric(nrow(data))  
    
    for (i in 1:nrow(data)) {
      counts <- str_count(string = raw_data_filtered[i, col], "1") # number of positive hits in a three-character sequence of 0s and 1s (triplicate scorings, e.g. 110 for two positive hits)
      if (counts > 0) {l[i] <- 1}
      else if (counts == 0) {l[i] <- 0}
    }
    encoded[[col]] <- l
  }
}
