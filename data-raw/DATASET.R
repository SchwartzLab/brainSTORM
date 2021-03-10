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
                "MRD_NaBH4HydBiotin.RTHIV_Mock.RTHIV",
                "MRD_NaBH4HydBiotin.TGIRT_Mock.TGIRT")

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

# Make sample STORM object. ####################################################

# Genome and Gene Annotation
fastaGenome <- "/home/labs/schwartzlab/miguelg/BIGDATA/genome_references/ribosomal/Sc_ribosomal_seqs.fa"
bedAnnotation <- "/home/labs/schwartzlab/miguelg/BIGDATA/gene_annotations/ribosomal/Sc_ribosomal_seqs.bed"
r1Files <- grep(list.files("/home/labs/schwartzlab/joeg/data/lib524/Saccharomyces_cerevisiae_SK1",
                           full.names = TRUE), pattern = "R1", value = TRUE)
OUTDIR <- "/home/labs/schwartzlab/miguelg/WORKSPACE_wexac/storm_seq/lib524/s_cerevisiae"
EXP_NAME <- "S.cerevisiae_SK1-lib524"
NCORES <- 10

vList <- list(organism = c("Yeast"),
              RTase = c("SSIII", "SSIV", "TGIRT", "RTHIV"),
              libTreat = c("Ac4C", "CMC", "DeacetylatedAc4C", "MocklowdNTPs",
                           "Mock", "Dimroth", "m5C", "AlkBmix", "RBSseqHeatMg", "NaBH4HydBiotin"))

META <- storm_META(fileNames = r1Files,
                   varsList = vList,
                   setVars = "organism",
                   idVars = c("organism", "RTase", "libTreat"),
                   groupVars = c("libTreat", "RTase"),
                   outDir = OUTDIR)

# Make transcriptome, make new gene annotation for transcriptome
fastaTxOme <- mkTranscriptome(fastaGenome, bedAnnotation, nCores = NCORES)
bedTxOme <- mkBedFromFastaTxOme(fastaTxOme)

# Creating bisulphite transcriptome
bisTxPath <- bisGenome(fastaTxOme)

# Create STAR genomes
STARGenome <- mkSTARgenome(fastaTxOme, bedTxOme)
STARGenome_bis <- mkSTARgenome(bisTxPath, bedTxOme)

#Loading transcriptome and annotation
GENOME <- txtools::tx_load_genome(fastaTxOme)
TXOME <- txtools::tx_load_bed(bedTxOme)

if(!all(file.exists(META$BAM))){
    # STAR Alignment
    alignSTAR(read1Files = META[libTreat != "m5C" & libTreat != "RBSseqHeatMg", FASTQ],
              nCores = NCORES,
              zipped = TRUE,
              STARgenomeDir = STARGenome,
              alignEndsType = "Local",
              outSAMtype = "BAM Unsorted",
              outDir = OUTDIR)
    alignSTAR(read1Files = META[libTreat == "m5C" | libTreat == "RBSseqHeatMg", FASTQ],
              nCores = NCORES,
              zipped = TRUE,
              STARgenomeDir = STARGenome_bis,
              alignEndsType = "Local",
              outSAMtype = "BAM Unsorted",
              outDir = OUTDIR)
}

if(!all(file.exists(META$RDS))){
    rdsFiles <- parallel::mclapply(mc.cores = NCORES, seq_along(META$FASTQ), function(i){
        bam2TxDT(BAMfile = META$BAM[i],
                 geneAnnot = TXOME,
                 genome = GENOME,
                 dtType = "covNuc",
                 outDir = OUTDIR,
                 nCores = 1,
                 remL = 1000,
                 minR = 0)
    })
}

yeast_STORM <- storm_STORM(META, GENOME, TXOME, nCores = 1)

usethis::use_data(yeast_STORM, overwrite = TRUE)
