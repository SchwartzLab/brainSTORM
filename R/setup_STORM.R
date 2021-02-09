#' STORM META constructor from RDS file names
#'
#' @param varsList
#' @param groupVars
#' @param setVars
#' @param idVars
#' @param fileNames
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
        lapply(function(x) tail(x, 1)) %>% unlist
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
        apply(MARGIN = 1, function(x){paste(x, collapse= ".")}) %>% factor()
    if(sum(duplicated(DT$id)) !=0){
        warning("Generated ids are not unique per sample")
    }
    if(sum(duplicated(paste(DT$group, DT$set)))){
        warning("Combinations of group and set should be unique between samples")
    }
    rootNames <- DT$FASTQ %>% lapply(function(x){strsplit(x, split = "/") %>%
            unlist %>% tail(1)}) %>% unlist %>% gsub(pattern = "_R1(.)+", replacement = "")
    BAM <- file.path(outDir, paste0(rootNames, "_Aligned.out.sorted.bam"))
    RDS <- file.path(outDir, paste0(rootNames, "_Aligned.out.sorted.txDT.rds"))
    DT <- tibble::add_column(.data = DT, .after = "FASTQ", BAM = BAM) %>%
    tibble::add_column(.after = "BAM", RDS = RDS)
    return(DT)
}


# Extract gene sequences from genome and geneAnnotation
getGeneSeqsfromGenome <- function(geneAnnot, genome, nCores = 1){
    txtools:::check_mc_windows(nCores)
    txtools:::check_GA_genome_chrCompat(geneAnnot = geneAnnot, genome = genome)
    parallel::mclapply(mc.cores = nCores, seq_along(geneAnnot), function(i){
        selGene <- geneAnnot[i]
        iGene <- as.character(selGene$name)
        iChr <- as.character(GenomicRanges::seqnames(selGene))
        iStr <- as.character(selGene@strand)
        iGA <- selGene
        iBlocks <- S4Vectors::mcols(iGA)$blocks %>% txtools:::if_IRangesList_Unlist() %>%
            IRanges::shift(IRanges::start(iGA) - 1)
        SEQ <- stringr::str_sub(genome[[iChr]], start = IRanges::start(iBlocks),
                                end = IRanges::end(iBlocks)) %>% paste(collapse = "") %>%
            Biostrings::DNAString()
        if(iStr == "-") {
            SEQ <- Biostrings::reverseComplement(SEQ)
        }
        SEQ
    }) %>% Biostrings::DNAStringSet()
}

# Extract transcriptome sequences from genome and generate FASTA file
mkTranscriptome <- function(fastaGenome, bedAnnotation, outFile = "auto", nCores){
    genome <- txtools::tx_load_genome(fastaGenome)
    geneAnnot <- txtools::tx_load_bed(bedAnnotation)
    seqs <- getGeneSeqsfromGenome(genome = genome, geneAnnot = geneAnnot, nCores)
    if(outFile == "auto"){
        mkTmpDir()
        fileN <- strsplit(fastaGenome, split = "/") %>% unlist %>% tail(1)
        outFile <- file.path("STORMtmp_dir", paste0(fileN, ".txOme"))
    }
    names(seqs) <- geneAnnot$name
    Biostrings::writeXStringSet(seqs, filepath = outFile)
    outFile
}

# Make GeneAnnotation for Transcriptome
mkBedFromFastaTxOme <- function(fastaTxOme, outFile = "auto"){
    fa <- txtools::tx_load_genome(fastaTxOme)
    if(max(Biostrings::width(fa)) > 10000){warning("Maximum sequence length exceeded 10000, make sure sequences represent transcripts.")}
    if(outFile == "auto"){
        mkTmpDir()
        fileN <- strsplit(fastaTxOme, split = "/") %>% unlist %>% tail(1)
        outFile <- file.path("STORMtmp_dir", paste0(fileN, ".bed"))
    }
    tmp <- GenomicRanges::GRanges(seqnames = names(fa),
                                  ranges = IRanges::IRanges(start = 1, width = Biostrings::width(fa)),
                                  strand = "+")
    names(tmp) <- names(fa)
    plyranges::write_bed(tmp, file = outFile)
    outFile
}

