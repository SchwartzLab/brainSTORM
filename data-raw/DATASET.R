## code to prepare `DATASET`
storm_baseCols <- c("chr", "gencoor", "strand", "gene", "txcoor", "pos",
                    "refSeq", "set", "nuc")
storm_baseCoorCols <- c("chr", "gencoor", "strand", "gene", "txcoor", "pos", "refSeq")

RNAmods_vec <- c("Y", "Nm", "m5C", "ac4C", "m1A", "m7G", "m3U")
usethis::use_data(storm_baseCols, overwrite = TRUE)
usethis::use_data(storm_baseCoorCols, overwrite = TRUE)
usethis::use_data(RNAmods_vec, overwrite = TRUE)
