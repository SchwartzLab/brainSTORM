## code to prepare `DATASET`
storm_baseCols <- c("chr", "gencoor", "strand", "gene", "txcoor", "pos",
                    "refSeq", "set", "nuc")
storm_baseCoorCols <- c("chr", "gencoor", "strand", "gene", "txcoor", "pos", "refSeq")

RNAmods_vec <- c("Y", "Nm", "m5C", "ac4C", "m1A", "m7G", "m3U")
usethis::use_data(storm_baseCols, overwrite = TRUE)
usethis::use_data(storm_baseCoorCols, overwrite = TRUE)
usethis::use_data(RNAmods_vec, overwrite = TRUE)

# Metrics to be assigned to specific RNAmods

Y_scores <- c("SRD1bpDS_CMC.TGIRT_Mock.TGIRT",
              "SRD1bpDS_CMC.SSIII_Mock.SSIII",
              "SRlog2FCh1bpDS_CMC.TGIRT_Mock.TGIRT",
              "SRlog2FCh1bpDS_CMC.SSIII_Mock.SSIII",
              "DRD_m5C.TGIRT_Mock.TGIRT",
              "DRD_RBSseqHeatMg.SSIII_Mock.SSIII",
              "DRD_RBSseqHeatMg.TGIRT_Mock.TGIRT")

Nm_scores <- c("NmStopScore_MocklowdNTPs.SSIII_Mock.SSIII",
               "ScoreA3p_Mock.TGIRT",
               "ScoreA3p_Mock.SSIII",
               "ScoreA3p_Mock.RTHIV")

ac4C_scores <- c("CtoT.MRD_Ac4C.TGIRT_Mock.TGIRT",
                 "CtoT.MRD_Ac4C.TGIRT_DeacetylatedAc4C.TGIRT",
                 "CtoT.MRD_Ac4C.SSIII_Mock.SSIII",
                 "CtoT.MRD_Ac4C.SSIII_DeacetylatedAc4C.SSIII")

m1A_scores <- c("SRD1bpDS_Mock.SSIII_Dimroth.SSIII",
                "SRlog2FCh1bpDS_Mock.SSIII_Dimroth.SSIII",
                "MRD_Mock.TGIRT_Dimroth.TGIRT",
                "MRD_Mock.SSIII_RBSseqHeatMg.SSIII",
                "MRD_Mock.TGIRT_RBSseqHeatMg.TGIRT",
                "MRD_Mock.TGIRT_AlkBmix.TGIRT")

m7G_scores <- c("MRD_DeacetylatedAc4C.TGIRT_Ac4C.TGIRT",
                "MRD_DeacetylatedAc4C.SSIII_Ac4C.SSIII",
                "SRD1bpDS_DeacetylatedAc4C.TGIRT_Ac4C.TGIRT",
                "SRD1bpDS_DeacetylatedAc4C.SSIII_Ac4C.SSIII",
                "SRlog2FCh1bpDS_DeacetylatedAc4C.SSIII_Ac4C.SSIII",
                "SRlog2FCh1bpDS_DeacetylatedAc4C.TGIRT_Ac4C.TGIRT",
                "SRD_NaBH4HydBiotin.RTHIV_Mock.RTHIV",
                "SRD_NaBH4HydBiotin.TGIRT_Mock.TGIRT")

m5C_scores <- c("CytPer_m5C.TGIRT_Mock.TGIRT",
                "CytPer_RBSseqHeatMg.SSIII_Mock.SSIII",
                "CytPer_RBSseqHeatMg.TGIRT_Mock.TGIRT")

m3U_scores <- c("MRD_Mock.TGIRT_AlkBmix.TGIRT")

usethis::use_data(Y_scores, overwrite = TRUE)
usethis::use_data(Nm_scores, overwrite = TRUE)
usethis::use_data(ac4C_scores, overwrite = TRUE)
usethis::use_data(m1A_scores, overwrite = TRUE)
usethis::use_data(m7G_scores, overwrite = TRUE)
usethis::use_data(m5C_scores, overwrite = TRUE)
usethis::use_data(m3U_scores, overwrite = TRUE)
