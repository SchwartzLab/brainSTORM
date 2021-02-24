#' STAR alignment
#'
#' Wrapper to perform alignment of sequencing reads to a reference genome
#' using STAR (Dobin) and sorting with Samtools.
#'
#' @param read1Files
#' @param STARgenomeDir
#' @param pairedEnd
#' @param zipped
#' @param nCores
#' @param alignEndsType
#' @param outSAMtype
#' @param outFilterMultimapNmax
#' @param outDir
#'
#' @return
#' @export
#'
#' @examples
alignSTAR <- function(read1Files, STARgenomeDir, pairedEnd = TRUE, zipped = TRUE,
                      nCores = 4, alignEndsType = "Local",
                      outSAMtype = "BAM Unsorted", outFilterMultimapNmax = 10,
                      outDir){
    mkTmpDir()
    if(!dir.exists(outDir)){dir.create(outDir)}
    if(zipped){rFCom <- "zcat"}else if(!zipped){rFCom <- "cat"}
    for(read1F in read1Files){
        if(pairedEnd){
            read2F <- gsub(read1F, pattern = "R1", replacement = "R2")
            if(!file.exists(read2F)){stop(read2F, "does not exist.")}
        }else if(!pairedEnd){
            read2F <- ""
        }else{stop("pairedEnd must be logical either TRUE or FALSE")}
        outFPrefix <- strsplit(read1F, split = "/") %>% unlist %>% tail(1) %>%
            gsub(pattern = "R1.fastq.gz", replacement = "") %>%
            gsub(read1F, pattern = "R1.fastq", replacement = "")
        outFPrefix <- file.path(getwd(), "STORMtmp_dir", outFPrefix)
        com <- paste0("module load STAR/2.7.5c &&",
                      "/apps/RH7U2/general/STAR/2.7.5c/bin/Linux_x86_64/STAR",
                      " --runMode alignReads",
                      " --runThreadN ", nCores,
                      " --genomeDir ", STARgenomeDir,
                      " --readFilesCommand ", rFCom,
                      " --readFilesIn ", read1F, " ", read2F,
                      " --outFileNamePrefix ", outFPrefix,
                      " --outSAMtype ", outSAMtype,
                      " --outFilterMultimapNmax ", outFilterMultimapNmax,
                      " --alignEndsType ", alignEndsType)
        system(com)
    }
    # Alignment Summary Report
    rootNames <- lapply(read1Files, function(x){strsplit(x, split = "/") %>% unlist %>%
            tail(1)}) %>% unlist %>% gsub(pattern = "_R1(.)+", replacement = "")
    logFiles <- file.path("STORMtmp_dir", paste0(rootNames, "_Log.final.out"))
    # all(file.exists(logFiles))
    RES <- lapply(logFiles, function(x){
        utils::read.delim(file = x, header = FALSE, stringsAsFactors = FALSE)
    })
    # Merge in one table
    summary <- lapply(seq_along(RES), function(x){
        RES[[x]][,2]
    }) %>% do.call(what = cbind)
    rownames(summary) <- RES[[1]][,1]
    colnames(summary) <- rootNames
    outReport <- file.path(outDir, "mappingSummary.txt")
    if(file.exists(outReport)){ # Add columns to existing summary report
        tmp <- data.table::fread(outReport, header = TRUE) %>%
            tibble::column_to_rownames("V1")
        utils::write.table(x = cbind(tmp, summary), file = outReport,
                           sep = "\t", quote = F, col.names = NA)
    }else{
        utils::write.table(x = summary, file = outReport, sep = "\t", quote = F,
                           col.names = NA)
    }
    # Sort and index with samtools
    bamFiles <- file.path("STORMtmp_dir", paste0(rootNames, "_Aligned.out.bam"))
    for(file in bamFiles){
        system(paste0("module load samtools/1.9 && ",
                      "/apps/RH7U2/gnu/samtools/1.9/bin/samtools sort -o ",
                      gsub(file, pattern = ".bam$", replacement = ".sorted.bam"),
                      " ", file, " -@", nCores))
        system(paste0("module load samtools/1.9 && ",
                      "/apps/RH7U2/gnu/samtools/1.9/bin/samtools index ",
                      gsub(file, pattern = ".bam$", replacement = ".sorted.bam")))
        system(paste0("rm ", file))
    }
    # Remove all garbage
    garbage <- file.path("STORMtmp_dir",
                         lapply(rootNames, function(x){
                             paste0(x, c("_SJ.out.tab", "_Log.progress.out", "_Log.out",
                                         "_Log.final.out"))
                         }) %>% unlist)
    invisible(file.remove(garbage))
    # Move files to output dir
    BAM <- c(gsub(bamFiles, pattern = ".bam$", replacement = ".sorted.bam"),
             gsub(bamFiles, pattern = ".bam$", replacement = ".sorted.bam.bai"))
    for(file in BAM){
        system(paste("mv", file, outDir))
    }
}

#' STORM object constructor from META table
#'
#' @param META
#' @param genome
#' @param geneAnnot
#' @param nCores
#'
#' @return
#' @export
#'
#' @examples
storm_STORM <- function(META, genome = NULL, geneAnnot = NULL, nCores = 1){
    if(hasDups(META$id)){
        stop("Not allowed duplicated id variables in META")
    }
    DTL <- lapply(META$RDS, function(x) readRDS(x)) %>% magrittr::set_names(META$id)# load RDS files
    # Check if data is in data.table or list format
    if(all(lapply(DTL, class) %>% unlist %>% magrittr::equals("data.frame"))){
        DTL <- lapply(DTL, data.table::data.table) %>% magrittr::set_names(META$id)
    }else if(all(lapply(DTL, class) %>% unlist %>% magrittr::equals("list"))){
        DTL <- lapply(DTL, function(x){
            do.call(x, what = rbind) %>% data.table::data.table()
        }) %>% magrittr::set_names(META$id)
    }
    # Have reference sequence, if not add it
    reqRefSeq <- lapply(DTL, function(x){"refSeq" %in% names(x)}) %>% unlist %>% magrittr::not()
    if(sum(reqRefSeq) > 0 & is.null(genome) | is.null(geneAnnot)){
        stop("Data requires reference sequence, genome and geneAnnot arguments must be provided")
    }
    if(sum(reqRefSeq) > 0){
        DTL[reqRefSeq] <- parallel::mclapply(mc.cores = nCores, DTL[reqRefSeq], function(DT){
            txtools::tx_split_DT(DT) %>% lapply(function(x){
                txtools::tx_add_refSeqDT(DT = x, genome = genome, geneAnnot = geneAnnot)
            }) %>% txtools::tx_merge_DT()
        })
    }
    # Check uniformity of DTs, if unequal equalize
    tmpPos <- lapply(DTL, function(DT){paste(DT$gene, DT$txcoor, sep = ":")})
    identCoors <- lapply(tmpPos[-1], function(x){identical(tmpPos[[1]], x)}) %>% unlist() %>% all()
    if(!identCoors){
        if(is.null(geneAnnot) | is.null(genome)){
            stop("Data needs compatibility adjustment, this requires geneAnnot and
                 genome arguments to be provided")
        }
        tmpGenes <- lapply(DTL, function(x) as.character(x$gene)) %>% unlist() %>% unique()
        tmpGA <- geneAnnot[match(tmpGenes, geneAnnot$name)]
        DTL <- lapply(DTL, function(DT){
            txtools::tx_complete_DT(DT, tmpGA, genome, nCores)
        })
    }
    if(!("pos" %in% names(DTL[[1]]))){
        DTL <- lapply(DTL, function(x) txtools::tx_add_pos(x))
    }
    STORM <- list(META = META,
                  DATA = DTL,
                  RES  = NULL)
    return(STORM)
}
