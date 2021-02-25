#' Create STORM META - Experiment design matrix
#'
#' STORM META constructor from RDS file names
#'
#' @param fileNames
#' @param varsList
#' @param groupVars
#' @param setVars
#' @param idVars
#' @param outDir
#'
#' @return
#' @export
#'
#' @examples
storm_META <- function(fileNames, varsList, groupVars, setVars, idVars, outDir){
    if(magrittr::not("libTreat" %in% names(varsList))){
        stop("libTreat must be one of the element names in varsList argument")
    }
    DT <- data.table::data.table(FASTQ = fileNames)
    procFileNames <- DT$FASTQ %>% stringr::str_split(pattern = "/") %>%
        lapply(function(x) utils::tail(x, 1)) %>% unlist
    DT <- lapply(names(varsList), function(i){
        patt <- paste0("(", paste(varsList[[i]], collapse = "|"), ")")
        tmp <- stringr::str_extract_all(string = procFileNames, pattern = patt)
        lapply(tmp, function(x) paste(x, collapse = ".")) %>% unlist %>%
            data.table::data.table() %>% magrittr::set_names(i)
    }) %>% do.call(what = cbind) %>% cbind(DT, .)
    DT$group <-  data.frame(DT[, groupVars, with = FALSE]) %>%
        apply(MARGIN = 1, function(x){paste(x, collapse= ".")}) %>% factor()
    DT$set <- data.frame(DT[, setVars, with = FALSE]) %>%
        apply(MARGIN = 1, function(x){paste(x, collapse= ".")}) %>% factor()
    DT$id <- data.frame(DT[, idVars, with = FALSE]) %>%
        apply(MARGIN = 1, function(x){paste(x, collapse= ".")}) %>% as.character()
    if(sum(duplicated(DT$id)) != 0){
        warning("Generated ids are not unique per sample")
    }
    if(sum(duplicated(paste(DT$group, DT$set)))){
        warning("Combinations of group and set should be unique between samples")
    }
    if(sum(DT$libTreat == "") != 0){
        warning("Some samples lack libTreat")
    }
    rootNames <- DT$FASTQ %>% lapply(function(x){strsplit(x, split = "/") %>%
            unlist %>% utils::tail(1)}) %>% unlist %>% gsub(pattern = "_R1(.)+", replacement = "")
    BAM <- file.path(outDir, paste0(rootNames, "_Aligned.out.sorted.bam"))
    RDS <- file.path(outDir, paste0(rootNames, "_Aligned.out.sorted.txDT.rds"))
    DT <- tibble::add_column(.data = DT, .after = "FASTQ", BAM = BAM) %>%
        tibble::add_column(.after = "BAM", RDS = RDS)
    return(DT)
}

#' BAM to txDataTable
#'
#' From BAM file to count data in a table
#'
#' @param BAMfile
#' @param geneAnnot
#' @param genome
#' @param dtType
#' @param paired
#' @param outDir
#' @param remL
#' @param minR
#' @param nCores
#' @param ySize
#' @param verb
#'
#' @return
#' @export
#'
#' @examples
bam2TxDT <- function(BAMfile, geneAnnot, genome, dtType, paired = TRUE,
                     outDir, remL = 10000, minR = 50, nCores = 2,
                     ySize = 100000, verb = TRUE){
    t0 <- Sys.time() # Start time
    # Process arguments
    # Load sequence or not
    if(dtType == "cov"){
        lSeq <- FALSE
    }else if(dtType == "covNuc"){
        lSeq <- TRUE
    }else{
        stop("dtType argument must be one of the options: 'cov' or 'covNuc'")
    }
    # Check yieldSize to be integer
    txtools:::check_integer_arg(ySize, "ySize")
    # Set OUTPUT filename
    outName <- strsplit(BAMfile, split = "/") %>% unlist %>% utils::tail(1) %>%
        gsub(x = ., pattern = ".bam$", replacement = ".txDT.rds", perl = T)
    # Output all genes even with no reads overlapping
    if(minR == 0){makeFULL <- TRUE}else{makeFULL <- FALSE}
    if(minR == 0){minR <- 1}
    # MAIN program
    # Load gene annotation
    gA <- geneAnnot
    if(verb){
        cat("Gene annotation loaded with", length(gA), "gene models.\n")
    }
    GENOME <- genome
    # Load BAM file
    bam <- txtools::tx_load_bam(file = BAMfile,
                                yieldSize = ySize,
                                scanFlag = "default",
                                loadSeq = lSeq,
                                recoverDumpedAligns = FALSE,
                                verbose = verb,
                                pairedEnd = paired)
    t1 <- Sys.time() # Loading BAM file time
    # Convert to transcriptomic
    txReads <- txtools::tx_reads(reads = bam,
                                 geneAnnot = gA,
                                 nCores = nCores,
                                 minReads = minR,
                                 withSeq = lSeq,
                                 verbose = verb)
    # Filter by length
    if(!is.na(remL)){
        txReads <- txtools::tx_filter_maxWidth(x = txReads, thr = remL, nCores = nCores)
    }

    # Data.table generation
    if(dtType == "cov"){
        if(verb){
            cat("Generating coverage data.table")
        }
        OUT <- txtools::tx_makeDT_coverage(x = txReads,
                                           geneAnnot = gA,
                                           nCores = nCores,
                                           fullDT = makeFULL,
                                           genome = GENOME)
    }else if(dtType == "covNuc"){
        if(verb){
            cat("Generating coverage and nucleotide frequency data.table. \n")
        }
        OUT <- txtools::tx_makeDT_covNucFreq(x = txReads,
                                             geneAnnot = gA,
                                             nCores = nCores,
                                             fullDT = makeFULL,
                                             genome = GENOME)
    }else{print("This message should not be printed, let the maintainer know.")}
    t2 <- Sys.time() # Creating data.table time
    # NOTE: In practice using too many cores made for longer processing times
    # newNCores <- min(nCores, 10)

    # Saving file as .rds
    saveRDS(object = OUT, file = file.path(outDir, outName))
    # Report
    timeBam <- t1 - t0 # Total time taken
    timePrc <- t2 - t1 # Total time taken
    timeTot <- t2 - t0 # Total time taken
    reportName <- strsplit(BAMfile, split = "/") %>% unlist %>% utils::tail(1) %>%
        gsub(x = ., pattern = ".bam$", replacement = ".txDT.log", perl = T)
    readsInOut <- parallel::mclapply(mc.cores = nCores, txReads, names) %>% unlist
    uniqReadsInOut <- unique(readsInOut)
    report <- c("BAM file name:", BAMfile,
                "Paired-end reads in BAM file:", length(bam),
                "Output contains:", " ",
                "Number of genes:", length(unique(OUT$gene)),
                "Number of reads in output:", length(readsInOut),
                "Number of unique reads in output:", length(uniqReadsInOut),
                "Fraction of total reads in output:", round(length(uniqReadsInOut)/length(bam), 4),
                "Loading BAM time:", paste(round(timeBam, 2), units(timeBam), sep = " "),
                "Processing time:", paste(round(timePrc, 2), units(timePrc), sep = " "),
                "Total time taken:", paste(round(timeTot, 2), units(timeTot), sep = " ")) %>%
        matrix(ncol = 2, byrow =T)
    utils::write.table(x = report,
                       file = file.path(outDir, reportName),
                       sep = "\t",
                       quote = FALSE,
                       row.names = FALSE,
                       col.names = FALSE)
    file.path(outDir, outName)
}

#' Make transcriptome
#'
#' Extract transcriptome sequences from genome and generate FASTA file
#'
#' @param fastaGenome
#' @param bedAnnotation
#' @param outFile
#' @param nCores
#'
#' @return
#' @export
#'
#' @examples
mkTranscriptome <- function(fastaGenome, bedAnnotation, outFile = "auto", nCores){
    genome <- txtools::tx_load_genome(fastaGenome)
    geneAnnot <- txtools::tx_load_bed(bedAnnotation)
    seqs <- getGeneSeqsfromGenome(genome = genome, geneAnnot = geneAnnot, nCores)
    if(outFile == "auto"){
        mkTmpDir()
        fileN <- strsplit(fastaGenome, split = "/") %>% unlist %>% utils::tail(1)
        outFile <- file.path("STORMtmp_dir", paste0(fileN, ".txOme"))
    }
    names(seqs) <- geneAnnot$name
    Biostrings::writeXStringSet(seqs, filepath = outFile)
    outFile
}

#' Make BED from FASTA Transcriptome
#'
#' Make BED Gene Annotation from a FASTA genome. Considering that the sequences
#' in the FASTA file represent intron-less genes.
#'
#' @param fastaTxOme
#' @param outFile
#'
#' @return
#' @export
#'
#' @examples
mkBedFromFastaTxOme <- function(fastaTxOme, outFile = "auto"){
    fa <- txtools::tx_load_genome(fastaTxOme)
    if(max(Biostrings::width(fa)) > 10000){warning("Maximum sequence length exceeded 10000, make sure sequences represent transcripts.")}
    if(outFile == "auto"){
        mkTmpDir()
        fileN <- strsplit(fastaTxOme, split = "/") %>% unlist %>% utils::tail(1)
        outFile <- file.path("STORMtmp_dir", paste0(fileN, ".bed"))
    }
    tmp <- GenomicRanges::GRanges(seqnames = names(fa),
                                  ranges = IRanges::IRanges(start = 1, width = Biostrings::width(fa)),
                                  strand = "+")
    names(tmp) <- names(fa)
    plyranges::write_bed(tmp, file = outFile)
    outFile
}

#' Add known RNA modifications nucleotide identity
#'
#' Add the nucleotide identity based on a table with columns "pos" and "nuc
#' establishing the position and nucleotide identity
#' Useful to annotate already known RNA modifications. The function will add
#' a column called "nuc" in the "RES" and "CALLS" tables.
#'
#' @param STORM STORM object
#' @param RNAmods data.frame with two columns "pos" for the positions in the
#' same format as the tables inside the STORM object as STORM$RES, a combination
#' in the form of "gene:txcoor"; And "nuc" for the identity of the nucleotide
#' on said position.
#'
#' @return
#' @export
#'
#' @examples
addKnownRNAmods <- function(STORM, RNAmods){
    if("RES" %in% names(STORM)){
        if("refSeq" %in% names(STORM$RES)){
            STORM$RES <- tibble::add_column(nuc = factor(RNAmods$nuc[match(STORM$RES$pos, RNAmods$pos)]),
                               .after = "refSeq",
                               .data = STORM$RES)
        }else{
            STORM$RES$nuc <- factor(RNAmods$nuc[match(STORM$RES$pos, RNAmods$pos)])
        }
    }
    if("CALLS" %in% names(STORM)){
        for(set in unique(STORM$META$set)){
            for(i in seq_along(STORM$CALLS[[set]])){
                tmp <- factor(RNAmods$nuc[match(STORM$CALLS[[set]][[i]]$pos, RNAmods$pos)])
                if("refSeq" %in% names(STORM$CALLS[[set]][[i]])){
                    STORM$CALLS[[set]][[i]] <- tibble::add_column(nuc = tmp,
                                                                  .after = "refSeq",
                                                                  .data = STORM$CALLS[[set]][[i]])
                }else{
                    STORM$CALLS[[set]][[i]] <- tibble::add_column(nuc = tmp,
                                                                  .data = STORM$CALLS[[set]][[i]])
                }
            }
        }
    }else{
        stop("RES is not part of the STORM object structure. ",
             "Metrics have not been calculated.")
    }
    STORM
}
