#' STAR alignment
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

#' STORM-summary
#'
#' Adds STORM$SUMMARY table
#'
#' @param STORM
#'
#' @return
#' @export
#'
#' @examples
storm_summary <- function(STORM){
    RNAmods <- names(STORM$CALLS[[1]])
    OUT <- lapply(RNAmods, function(RNAmod_i){
        tmp_out <- STORM$DATA[[1]][,storm_baseCoorCols, with = FALSE]
        SETS <- levels(STORM$META$set)
        if(is.null(STORM$CALLS[[SETS[1]]][[RNAmod_i]]$logist_Score)){
            warning("logist_Score was not found for ", RNAmod_i, ". Results will",
                    " not be summarized in STORM$SUMMARY. Check if metrics where ",
                    "placed in STORM$CALLS[[1]]$", RNAmod_i)
            return(NULL)
        }else{
            tmp_log <- lapply(SETS, function(x){
                STORM$CALLS[[x]][[RNAmod_i]]$logist_Score
            }) %>% do.call(what = cbind) %>% magrittr::set_colnames(paste("logScore",
                                                                          SETS, sep = "_"))
            tmp_lin <- lapply(SETS, function(x){
                STORM$CALLS[[x]][[RNAmod_i]]$linear_Score
            }) %>% do.call(what = cbind) %>% magrittr::set_colnames(paste("linScore",
                                                                          SETS, sep = "_"))
            tmp_prd <- lapply(SETS, function(x){
                STORM$CALLS[[x]][[RNAmod_i]]$pred
            }) %>% do.call(what = cbind) %>% magrittr::set_colnames(paste("pred",
                                                                          SETS, sep = "_"))
            tmp_out <- cbind(tmp_out, tmp_lin, tmp_log, tmp_prd)
            return(tmp_out)
        }
    }) %>% magrittr::set_names(RNAmods)
    STORM$SUMMARY <- OUT
    STORM
}


# Add results to a STORM object. Remove scores if metric is already present.
hlpr_add_REScols <- function(STORM_RES, REScols){
    iMetric <- unique(REScols[,"metric"]) %>% as.character()
    # remove results for identical metric
    if("metric" %in% names(STORM_RES)){
        STORM_RES <- STORM_RES[STORM_RES$metric != iMetric,]
    }
    STORM_RES <- rbind(STORM_RES, REScols)
    STORM_RES
}

# List files in work dir that match pattern pat
listFilePatt <- function(pattern, path = "."){
    files <- list.files(path)[grep(pattern = pattern, x = list.files(path))]
    return(files)
}


# Notebook 1 ###################################################################

# Tables of files from FASTQ to expected BAM and RDS targets
files_table <- function(META){
    fastq <- META$FASTQ
    bam <- META$BAM
    lce <- gsub(pattern = "Aligned.out.sorted.bam",
                replacement = "Aligned.out.sorted.lce.txt", META$BAM)
    rds <- gsub(pattern = "Aligned.out.sorted.bam",
                replacement = "Aligned.out.sorted.txDT.rds", META$BAM)
    # Table
    tmpDT <- data.table::data.table(FASTQ = c(fastq),
                                    BAM = c(bam),
                                    BAM_ok = file.exists(c(bam)),
                                    lce = c(lce),
                                    lce_ok = file.exists(c(lce)),
                                    rds = c(rds),
                                    rds_ok = file.exists(c(rds)))
    return(tmpDT)
}

# FASTQ Duplication rate (library complexity)
fastq_dupRate <- function(FASTQs_pahts, nCores){
    parallel::mclapply(mc.cores = nCores, FASTQs_pahts, function(file){
        tmp <- ShortRead::readFastq(file)
        tmp2 <- ShortRead::readFastq(gsub(file, pattern = "R1", replacement = "R2"))
        dupR <- paste(tmp@sread, tmp2@sread, sep = "") %>% duplicated() %>% mean
        return(dupR)
    }) %>% unlist
}

#' Nucleotide frequency report
#'
#' @param META
#' @param nCores
#' @param firstN
#'
#' @return
#' @export
#'
#' @examples
fastq_nucFreq <- function(META, nCores, firstN = 1e4){
    parallel::mclapply(mc.cores = nCores, META$FASTQ, function(file){
        tmp <- readLines(file, firstN * 4)[seq(2, firstN*4, 4)] %>% Biostrings::DNAStringSet()
        tmp2 <- readLines(gsub(file, pattern = "R1",
                               replacement = "R2"), firstN * 4)[seq(2, firstN*4, 4)] %>%
            Biostrings::DNAStringSet() %>% Biostrings::complement()
        mR1R2 <- paste(tmp, tmp2, sep = "")
        nucFreq <- mR1R2 %>% stringr::str_split(pattern = "") %>% unlist %>% table
        return(nucFreq[c("A", "C", "G", "T")])
    }) %>% do.call(what = cbind) %>% magrittr::set_colnames(META$id)
}

#' nucleotide frequency plot
#'
#' @param nucF_x
#' @param subtitle
#'
#' @return
#' @export
#'
#' @examples
gg_nucFreq <- function(nucF_x, subtitle){
    tmp <- prop.table(nucF_x, margin = 2) %>% data.frame %>% tibble::rownames_to_column(var = "nuc") %>%
        tidyr::pivot_longer(cols = -nuc, names_to = "Sample", values_to = "Ratio")
    ggplot2::ggplot(tmp, ggplot2::aes(x = Sample, y = Ratio, fill = nuc)) +
        ggplot2::geom_bar(stat = "identity") +
        ggplot2::scale_fill_brewer(palette="Set1") +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) +
        ggplot2::ggtitle("Nucleotide Frequency per library", subtitle = subtitle)
}

#' Alignment report table
#'
#' @param META
#' @param nCores
#'
#' @return
#' @export
#'
#' @examples
reads_report <- function(META, nCores = 4){
    DT <- files_table(META)
    if(all(DT$BAM_ok) & all(DT$lce_ok) & all(DT$rds_ok)){
        res1 <- parallel::mclapply(mc.cores = nCores, DT$FASTQ, function(x){
            tmp <- ShortRead::readFastq(x)
            length(tmp)
        }) %>% unlist
        res2 <- parallel::mclapply(mc.cores = nCores, DT$BAM[DT$BAM_ok], function(x){
            tmp2 <- Rsamtools::scanBam(x)
            tmp2[[1]]$qname %>% unique %>% length
        }) %>% unlist
        res3 <- parallel::mclapply(mc.cores = nCores, DT$rds[DT$rds_ok], function(x){
            tmplog <- data.table::fread(gsub(pattern = "rds", "log", x), header = F)
            tmplog[grep(tmplog$V1, pattern =  "unique reads"),]$V2 %>% as.numeric()
        }) %>% unlist
        tmpDT <- data.table::data.table(sample = META$id,
                                        FASTQ_reads = res1,
                                        BAM_aligns = NA,
                                        tx_starts = NA)
        tmpDT$BAM_aligns[DT$BAM_ok] <- res2
        tmpDT$tx_starts[DT$rds_ok] <- res3
        tmpDT$pC_BAM <- round(tmpDT$BAM_aligns / tmpDT$FASTQ_reads * 100, 2)
        tmpDT$pC_tx <- round(tmpDT$tx_starts / tmpDT$BAM_aligns * 100, 2)
        return(tmpDT)
    }else{
        stop("Some files are missing.")
    }
}

#' Alignment report plots
#'
#' @param rReport
#' @param species
#'
#' @return
#' @export
#'
#' @examples
gg_readStats <- function(rReport, species){
    tmpDT <- rReport[, -c(5, 6)]
    tmpDT$FASTQ <- tmpDT$FASTQ_reads - tmpDT$BAM_aligns
    tmpDT$BAM <- tmpDT$BAM_aligns - tmpDT$tx_starts
    tmpDT$txDT <- tmpDT$tx_starts
    tmpDT$BAM[is.na(tmpDT$BAM)] <- 0; tmpDT$txDT[is.na(tmpDT$txDT)] <- 0
    tmpDT <- tidyr::pivot_longer(data = tmpDT[,c("sample", "FASTQ", "BAM", "txDT")],
                                 cols = c("FASTQ", "BAM", "txDT"), names_to = "Reads")
    tmpDT$Reads <- factor(tmpDT$Reads, levels = c("FASTQ", "BAM", "txDT"))
    tmpDT$value <- magrittr::divide_by(tmpDT$value, 1e6)
    t_GG1 <- ggplot2::ggplot(tmpDT, ggplot2::aes(x = sample, y = value, fill = Reads)) +
        ggplot2::geom_bar(stat="identity") +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) +
        ggplot2::ggtitle(paste("Reads in", species, "libraries by processing step"),) +
        ggplot2::ylab("Million reads") +
        ggplot2::xlab("Samples")
    # Proportion
    tmpDT <- rReport[, -c(5, 6)]
    tmpDT$FASTQ <- tmpDT$FASTQ_reads - tmpDT$BAM_aligns
    tmpDT$BAM <- tmpDT$BAM_aligns - tmpDT$tx_starts
    tmpDT$txDT <- tmpDT$tx_starts
    tmpDT$BAM[is.na(tmpDT$BAM)] <- 0; tmpDT$txDT[is.na(tmpDT$txDT)] <- 0
    tmpDT <- data.frame(sample = tmpDT$sample,
                        apply(tmpDT[,c("FASTQ", "BAM", "txDT")], 1, prop.table) %>% t)
    tmpDT <- tidyr::pivot_longer(data = tmpDT[,c("sample", "FASTQ", "BAM", "txDT")],
                                 cols = c("FASTQ", "BAM", "txDT"), names_to = "Reads")
    tmpDT$Reads <- factor(tmpDT$Reads, levels = c("FASTQ", "BAM", "txDT"))

    t_GG2 <- ggplot2::ggplot(tmpDT, ggplot2::aes(x = sample, y = value, fill = Reads)) +
        ggplot2::geom_bar(stat="identity") +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) +
        ggplot2::ggtitle(paste("Proportion of reads in", species, "libraries by processing step"),) +
        ggplot2::ylab("Proportion") +
        ggplot2::xlab("Samples")
    return(list(t_GG1, t_GG2))
}

# Alignment and transcript data processing efficiency
ggAlignEffPlot <- function(META, rReport){
    tmp <- cbind(META, rReport[,-1])
    tmpGG1 <- ggplot2::ggplot(tmp, ggplot2::aes(x = tmp$libTreat, y = tmp$pC_BAM, colour = tmp$libTreat)) +
        ggplot2::geom_boxplot() + ggplot2::geom_point(colour = "black") +
        ggplot2::ggtitle(paste(META$organism[1], "- Alignment efficiency")) +
        ggplot2::ylab("% Aligned reads") +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))
    tmpGG2 <- ggplot2::ggplot(tmp, ggplot2::aes(x = tmp$bioTreat, y = tmp$pC_BAM, colour = tmp$bioTreat)) +
        ggplot2::geom_boxplot() + ggplot2::geom_point(colour = "black") +
        ggplot2::ggtitle(paste(META$organism[1], "- Alignment efficiency")) +
        ggplot2::ylab("% Aligned reads") +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))
    tmpGG3 <- ggplot2::ggplot(tmp, ggplot2::aes(x = tmp$RTase, y = tmp$pC_BAM, colour = tmp$RTase)) +
        ggplot2::geom_boxplot() + ggplot2::geom_point(colour = "black") +
        ggplot2::ggtitle(paste(META$organism[1], "- Alignment efficiency")) +
        ggplot2::ylab("% Aligned reads") +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))
    list(tmpGG1, tmpGG2, tmpGG3)
}

#' Library complexity extrapolation plot
#'
#' @param META
#' @param tab_name
#' @param speciesName
#'
#' @return
#' @export
#'
#' @examples
gg_lce <- function(META, tab_name, speciesName = ""){
    lceFiles <- gsub(META$BAM, pattern = ".bam$", replacement = ".lce.txt") %>%
        magrittr::set_names(META$id)
    if(!all(file.exists(lceFiles))){
        stop("Report files missing:\n", paste(lceFiles[!file.exists(lceFiles)], collapse = " \n"))
    }
    tmp <- lapply(lceFiles, function(x) data.table::fread(x)) %>%
        magrittr::set_names(META$id)
    tmp <- lapply(names(tmp), function(x){
        cbind(tmp[[x]], id = x)
    }) %>% do.call(what = rbind) %>% data.table::data.table()
    data.table::fwrite(x = tmp, file = tab_name, sep = "\t")
    cat("Library complexity table output:", tab_name)
    tmpT <- table(tmp$TOTAL_READS)
    e_reads <- names(tmpT)[tmpT == max(tmpT)] %>% as.numeric %>% max
    tmp <- tmp[tmp$TOTAL_READS == e_reads,]
    ggOUT <- ggplot2::ggplot(tmp, ggplot2::aes(x = tmp$id, y = tmp$EXPECTED_DISTINCT)) +
        ggplot2::geom_bar(stat="identity", color="black", position= ggplot2::position_dodge()) +
        ggplot2::geom_errorbar(ggplot2::aes(ymin= tmp$LOWER_0.95CI, ymax= tmp$UPPER_0.95CI), width=.2,
                               position= ggplot2::position_dodge(1)) +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) +
        ggplot2::ggtitle(label = paste("Library complexity extrapolation -", speciesName),
                         sub = paste("Expected distinct reads at", e_reads, "depth. CI = 95%")) +
        ggplot2::ylab("Expected distinct reads") +
        ggplot2::xlab("Samples")
    return(ggOUT)
}

# STORM add functions ##########################################################

TGIRT_metrics <- function(STORM){
    STORM %>%
        add_CytPer("m5C", "Mock") %>%
        add_CtoT_MRD("Ac4C", "Mock") %>%
        add_CtoT_MRD("Ac4C", "DeacetylatedAc4C") %>%
        add_MRD("Mock", "Dimroth") %>%
        add_SRD1bpDS("CMC", "Mock") %>%
        add_SRlog2FCh1bpDS("CMC", "Mock") %>%
        storm_add_scoreA3p("Mock", minMedCov = 50) %>%
        add_MRD("DeacetylatedAc4C", "Ac4C")
}

SSIII_metrics <- function(STORM){
    STORM %>%
        add_CtoT_MRD("Ac4C", "Mock") %>%
        add_CtoT_MRD("Ac4C", "DeacetylatedAc4C") %>%
        add_SRD1bpDS("Mock", "Dimroth") %>%
        add_SRD1bpDS("CMC", "Mock") %>%
        add_SRlog2FCh1bpDS("CMC", "Mock") %>%
        add_2OmeScore("MocklowdNTPs", "Mock", perNuc_Zscaling = T) %>%
        add_MRD("DeacetylatedAc4C", "Ac4C") %>%
        storm_add_scoreA3p("Mock") %>%
        add_SRD1bpDS("DeacetylatedAc4C", "Ac4C")
}


# Boxplot of $RES calculated metrics grouping by modified nucleotides
ggMetricsNuc <- function(STORM, title){
    tmpDT <- STORM$RES %>% stats::na.omit()
    ggplot2::ggplot(tmpDT, ggplot2::aes(x = nuc, y = score)) + ggplot2::geom_boxplot(outlier.colour = NA) +
        ggplot2::geom_point(ggplot2::aes(colour = metric), alpha = 0.2) +
        ggplot2::theme_bw() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) +
        ggplot2::facet_wrap(~metric, scales = "free") +
        ggplot2::ggtitle(title)
}

# Plot scatterplot of scores colored by gene by tx coordinate
ggMetricsPos <- function(STORM, title){
    tmpDT <- STORM$RES %>% stats::na.omit()
    ggplot2::ggplot(tmpDT) +
        ggplot2::geom_point(ggplot2::aes(x = tmpDT$pos, y = tmpDT$score, colour = tmpDT$gene), alpha = 0.2) +
        ggplot2::facet_grid(tmpDT$metric~tmpDT$set, scales = "free") +
        ggplot2::theme_bw() +
        ggplot2::theme(axis.text.x = ggplot2::element_blank()) +
        ggplot2::ggtitle(title) + ggplot2::xlab("txCoor")
}

# RNAmod logical vectors from character or factor
is.pseudoU <- function(nucs){nucs %in% c("Y", "Ym")}
is.2Ome <- function(nucs){nucs %in% c("Am", "Gm", "Um", "Cm", "Ym")}
is.m5C <- function(nucs){nucs %in% c("m5C")}
is.ac4C <- function(nucs){nucs %in% c("ac4C")}
is.m1A <- function(nucs){nucs %in% c("m1A")}
is.m7G <- function(nucs){nucs %in% c("m7G")}

#  Balanced groups assignment of two level vectors for K-fold cross validation
CV_balancedGroups <- function(x, k){
    x <- factor(x)
    blocks <- seq(1, k)
    groups <- rep(NA, length(x))
    tmpLog1 <- x == levels(x)[1]
    tmpLog2 <- x == levels(x)[2]
    groups[tmpLog1] <- sample(rep(blocks, ceiling(sum(tmpLog1)/length(blocks))))[seq(groups[tmpLog1])]
    groups[tmpLog2] <- sample(rep(blocks, ceiling(sum(tmpLog2)/length(blocks))))[seq(groups[tmpLog2])]
    groups
}
# TODO: Make a generalized K-fold balanced groups sampler

# Extract metrics from k-CFV res
extract_AUC_sen_spe <- function(res, RNAmod, strategy, organism){
    data.frame(organism = organism,
               strategy = strategy,
               RNAmod = RNAmod,
               AUC = unlist(lapply(seq(res), function(i) res[[i]][RNAmod, "AUC"])),
               sensitivity = unlist(lapply(seq(res), function(i) res[[i]][RNAmod, "Sensitivity"])),
               specificity = unlist(lapply(seq(res), function(i) res[[i]][RNAmod, "Specificity"])))
}

# Matthew Correlation Coefficient
matthewCorrCoeff <- function(scores, successes){
    MCC <- NULL
    known <- rep(0, length(scores)); known[successes] <- 1
    NAs <- which(is.na(scores))
    scores[is.na(scores)] <- 0
    for (i in 1:length(scores)){
        thr <- scores[i]
        calls <- as.numeric(scores >= thr)
        TP <-  sum(calls & known)
        FP <-  sum(calls) - TP
        TN <- sum(!calls & !known)
        FN <- sum(!calls) - TN
        MCC[i] <- (TP * TN - FP * FN) / sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
    }
    MCC[is.nan(MCC)] <- 0
    MCC[is.infinite(MCC)] <- 0
    MCC[NAs] <- NA
    return(MCC)
}

# MCC, using threshold removing NA scores
MCC <- function(scores, successes, thr){
    successes <- successes[!is.na(scores)]
    scores <- scores[!is.na(scores)]
    calls <- as.numeric(scores >= thr)
    known <- rep(0, length(scores)); known[successes] <- 1
    TP <-  sum(calls & known)
    FP <-  sum(calls) - TP
    TN <- sum(!calls & !known)
    FN <- sum(!calls) - TN
    MCC <- (TP * TN - FP * FN) / sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
    if(is.nan(MCC)){MCC <- 0}
    MCC
}

# All combinations of strings in character vector
all_comb_allNelements <- function(x){
    lapply(seq_along(x), function(y) utils::combn(x, y, simplify = FALSE)) %>%
        unlist(recursive = FALSE)
}

# Logistic Model training using STORM$CALLS object.
trainLogModel_RNAmods <- function(STORM, isModFun, modNuc, thresholds){
    CALLS <- STORM$CALLS[[1]][[modNuc]] %>% data.frame()
    CALLS$isMod <- isModFun(CALLS$nuc)
    varNames <- names(CALLS)[lapply(CALLS, is.double) %>% unlist]
    allCombVar <- all_comb_allNelements(varNames)
    combNames <- paste0("Comb_", seq_along(allCombVar))
    OUT <- lapply(allCombVar, function(varNames){
        tData <- CALLS[, c(varNames, "isMod")]
        tData$isMod <- as.numeric(tData$isMod)
        logMod <- stats::glm(isMod ~ ., family = stats::binomial(), data = tData)
        CALLS$logPred <- stats::predict(logMod, newdata = CALLS, type = "response")
        sapply(thresholds, function(x) MCC(CALLS$logPred, CALLS$isMod, x))
    }) %>% do.call(what = cbind) %>% magrittr::set_rownames(thresholds) %>% magrittr::set_colnames(combNames)
    selComb <- which.max(apply(OUT, 2, median))
    selThr <- which.max(OUT[,selComb])
    # FinalModel
    tData <- CALLS[, c(allCombVar[[selComb]], "isMod")]
    tData$isMod <- as.numeric(tData$isMod)
    logMod <- stats::glm(isMod ~ ., family = stats::binomial(), data = tData)
    return(list(logiMod = logMod, thr = thresholds[selThr], vars = allCombVar[[selComb]],
                MCCmat = OUT, MCC= OUT[selThr, selComb]))
}

#' Logistic regressions scores
#'
#' Update STORM$CALLS with logistic and linear scores and prediction
#'
#' @param STORM
#' @param logMods
#'
#' @return
#' @export
#'
#' @examples
call_RNAmods_logRegMods <- function(STORM, logMods){
    STORM$CALLS <- lapply(names(STORM$CALLS), function(callList){
        tmp <- lapply(names(STORM$CALLS[[callList]]), function(RNAmod){
            if(STORM$CALLS[[callList]][[RNAmod]] %>% lapply(is.double) %>% unlist %>% sum == 0){
                return(STORM$CALLS[[callList]][[RNAmod]])
            }
            nData <- data.frame(STORM$CALLS[[callList]][[RNAmod]])
            STORM$CALLS[[callList]][[RNAmod]]$logist_Score <-
                stats::predict.glm(logMods[[RNAmod]]$logiMod, newdata = nData, type = "response")
            STORM$CALLS[[callList]][[RNAmod]]$linear_Score <-
                stats::predict.glm(logMods[[RNAmod]]$logiMod, newdata = nData, type = "link")
            STORM$CALLS[[callList]][[RNAmod]]$pred <-
                STORM$CALLS[[callList]][[RNAmod]]$logist_Score >= logMods[[RNAmod]]$thr
            return(STORM$CALLS[[callList]][[RNAmod]])
        }) %>% magrittr::set_names(names(STORM$CALLS[[callList]]))
        return(tmp)
    }) %>% magrittr::set_names(names(STORM$CALLS))
    STORM
}

# Save as RDS if object is not already RDS
saveAsRDSIfNotAlready <- function(object, fileName){
    if(!file.exists(fileName)){saveRDS(object, fileName)}
}

# Notebook 4 ########################

# Has duplicates?
hasDups <- function(x){
    sum(duplicated(x)) > 0
}

# Limit STORM object to genes in geneAnnot argument
storm_reduceToGA <- function(STORM, geneAnnot){
    STORM$DATA <- lapply(STORM$DATA, function(DT) DT[gene %in% geneAnnot$name,])
    return(STORM)
}

# storm_calls: Makes predictions based on RF and CutPointR for each modification
storm_calls <- function(STORM, RNAmodList, RF_list, CP_models){
    RNAmods_vec <- names(RNAmodList)
    coorSys <- STORM$RES[,names(STORM$RES) %in% c("chr", "gencoor", "strand",
                                                  "gene", "txcoor", "pos",
                                                  "refSeq", "nuc"), with = FALSE]
    iSets <- levels(STORM$META$set)
    STORM$CALLS <- lapply(seq_along(RNAmods_vec), function(i){
        pred_RES <- lapply(seq_along(iSets), function(j){
            lapply(CP_models[[i]], function(CPmod){
                tmpPat <- paste0("^", gsub("_TGIRT", "", CPmod$predictor))
                selVar <- grepl(pattern = tmpPat, names(STORM$RES)) &
                    grepl(pattern = iSets[j], names(STORM$RES))
                tmpDat <- data.frame(x = STORM$RES[[which(selVar)]])
                CPmod$predictor <- "x"
                as.logical(stats::predict(CPmod, tmpDat))
            }) %>% do.call(what = data.frame) %>% unname %>% apply(1, all)
        }) %>% do.call(what = cbind) %>% data.table::data.table() %>%
            magrittr::set_names(paste("pred", RNAmods_vec[i], iSets, sep = "_"))

        scor_RES <- lapply(seq_along(iSets), function(j){
            lapply(CP_models[[i]], function(CPmod){
                tmpPat <- paste0("^", gsub("_TGIRT", "", CPmod$predictor))
                selVar <- grepl(pattern = tmpPat, names(STORM$RES)) &
                    grepl(pattern = iSets[j], names(STORM$RES))
                STORM$RES[,which(selVar), with = FALSE]
            }) %>% do.call(what = cbind)
        }) %>% do.call(what = cbind) %>% data.table::data.table()

        tree_RES <- lapply(seq_along(iSets), function(j){
            RF_mod <- RF_list[[i]]
            RF_vars <- names(RF_mod$forest$xlevels)
            newDat <- lapply(RF_vars, function(RF_v){
                tmpP <- gsub(pattern = "_TGIRT", replacement = "", x = RF_v)
                selVar <- grepl(pattern = tmpP, x = names(STORM$RES))
                selSmp <- grepl(pattern = iSets[j], x = names(STORM$RES)) & selVar
                newDat <- data.frame(x = STORM$RES[, selSmp, with = FALSE]) %>% magrittr::set_names(RF_v)
            }) %>% do.call(what = cbind)
            newDat[is.na(newDat)] <- 0
            out <- stats::predict(RF_mod, newDat, norm.votes = TRUE, type = "vote")
            unname(out[,"TRUE"])
        }) %>% do.call(what = cbind) %>% data.table::data.table() %>%
            magrittr::set_names(paste("votes", RNAmods_vec[i], iSets, sep = "_"))

        data.table::data.table(coorSys, pred_RES, tree_RES, scor_RES)
    }) %>% magrittr::set_names(RNAmods_vec)
    return(STORM)
}

