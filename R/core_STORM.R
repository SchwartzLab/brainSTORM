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
        stop("-d --DT_datatype argument must be one of the options: 'cov' or 'covNuc'")
    }
    # Check yieldSize to be integer
    txtools:::check_integer_arg(ySize, "ySize")
    # Set OUTPUT filename
    outName <- strsplit(BAMfile, split = "/") %>% unlist %>% utils::tail(1) %>%
        gsub(x = ., pattern = ".bam$", replacement = ".txDT.rds", perl = T)
    if(!(dtType == "cov" | dtType == "covNuc")){
        stop("-d parameter is neither 'cov', nor 'covNuc'")
    }
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
    GenomicAlignments::width(txReads) %>% lapply(function(x) stats::quantile(x, 0.99))
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
