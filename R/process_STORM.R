# STAR alignment
# Wrapper to perform alignment of sequencing reads to a reference genome
# using STAR (Dobin) and sorting with Samtools.
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
        com <- paste0("/apps/RH7U2/general/STAR/2.7.5c/bin/Linux_x86_64/STAR",
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
        read.delim(file = x, header = FALSE, stringsAsFactors = FALSE)
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
        write.table(x = cbind(tmp, summary), file = outReport,
                    sep = "\t", quote = F, col.names = NA)
    }else{
        write.table(x = summary, file = outReport, sep = "\t", quote = F,
                    col.names = NA)
    }
    # Sort and index with samtools
    bamFiles <- file.path("STORMtmp_dir", paste0(rootNames, "_Aligned.out.bam"))
    for(file in bamFiles){
        system(paste0("/apps/RH7U2/gnu/samtools/1.9/bin/samtools sort -o ",
                      gsub(file, pattern = ".bam$", replacement = ".sorted.bam"),
                      " ", file, " -@", nCores))
        system(paste0("/apps/RH7U2/gnu/samtools/1.9/bin/samtools index ",
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

# STORM object constructor from META table
storm_STORM <- function(META, genome = NULL, geneAnnot = NULL, nCores = 1){
    if(hasDups(META$id)){
        stop("Not allowed duplicated id variables in META")
    }
    DTL <- lapply(META$RDS, function(x) readRDS(x)) # load RDS files
    # Check if data is in data.table or list format
    if(all(lapply(DTL, class) %>% unlist %>% magrittr::equals("data.frame"))){
        DTL <- lapply(DTL, data.table) %>% magrittr::set_names(META$id)
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

# STORM$SUMMARY table
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
            }) %>% do.call(what = cbind) %>% set_colnames(paste("logScore",
                                                                SETS, sep = "_"))
            tmp_lin <- lapply(SETS, function(x){
                STORM$CALLS[[x]][[RNAmod_i]]$linear_Score
            }) %>% do.call(what = cbind) %>% set_colnames(paste("linScore",
                                                                SETS, sep = "_"))
            tmp_prd <- lapply(SETS, function(x){
                STORM$CALLS[[x]][[RNAmod_i]]$pred
            }) %>% do.call(what = cbind) %>% set_colnames(paste("pred",
                                                                SETS, sep = "_"))
            tmp_out <- cbind(tmp_out, tmp_lin, tmp_log, tmp_prd)
            return(tmp_out)
        }
    }) %>% set_names(RNAmods)
    STORM$SUMMARY <- OUT
    STORM
}


# Add results to a STORM object. Remove scores if metric is already present.
hlpr_add_REScols <- function(STORM_RES, REScols){
    iMetric <- unique(REScols[,"metric"]) %>% as.character()
    # remove results for identical metric
    if("metric" %in% names(STORM_RES)){
        STORM_RES <- STORM_RES[metric != iMetric,]
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
files_table <- function(META, outDir){
    fastq <- META$FASTQ
    if(!"BAM" %in% colnames(META)){
        META <- addBAMFileNames(META, outDir)
    }
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
    mclapply(mc.cores = nCores, FASTQs_pahts, function(file){
        tmp <- ShortRead::readFastq(file)
        tmp2 <- ShortRead::readFastq(gsub(file, pattern = "R1", replacement = "R2"))
        dupR <- paste(tmp@sread, tmp2@sread, sep = "") %>% duplicated() %>% mean
        return(dupR)
    }) %>% unlist
}

# FASTQ nucleotide frequency
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

## ggplot nucleotide frequency barplots
gg_nucFreq <- function(nucF_x, subtitle){
    tmp <- prop.table(nucF_x, margin = 2) %>% data.frame %>% tibble::rownames_to_column(var = "nuc") %>%
        tidyr::pivot_longer(cols = -nuc, names_to = "Sample", values_to = "Ratio")
    ggplot2::ggplot(tmp, ggplot2::aes(x = Sample, y = Ratio, fill = nuc)) +
        ggplot2::geom_bar(stat = "identity") +
        ggplot2::scale_fill_brewer(palette="Set1") +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) +
        ggplot2::ggtitle(paste("Nucleotide Frequency per library in", subtitle))
}

# Reads report table
reads_report <- function(DT, META, nCores = 4){
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
}

# read stats plots with ggplots
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
    tmpGG1 <- ggplot(tmp, aes(x = libTreat, y = pC_BAM, colour = libTreat)) +
        geom_boxplot() + geom_point(colour = "black") +
        ggtitle(paste(META$organism[1], "- Alignment efficiency")) +
        ylab("% Aligned reads") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1))
    tmpGG2 <- ggplot(tmp, aes(x = bioTreat, y = pC_BAM, colour = bioTreat)) +
        geom_boxplot() + geom_point(colour = "black") +
        ggtitle(paste(META$organism[1], "- Alignment efficiency")) +
        ylab("% Aligned reads") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1))
    tmpGG3 <- ggplot(tmp, aes(x = RTase, y = pC_BAM, colour = RTase)) +
        geom_boxplot() + geom_point(colour = "black") +
        ggtitle(paste(META$organism[1], "- Alignment efficiency")) +
        ylab("% Aligned reads") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1))
    list(tmpGG1, tmpGG2, tmpGG3)
}

# Library complexity extrapolation barplots and tables
gg_lce <- function(META, tab_name, speciesName = ""){
    lceFiles <- gsub(META$BAM, pattern = ".bam", replacement = ".lce.txt") %>%
        setNames(META$id)
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
    ggOUT <- ggplot2::ggplot(tmp, ggplot2::aes(x = id, y = EXPECTED_DISTINCT)) +
        ggplot2::geom_bar(stat="identity", color="black", position= ggplot2::position_dodge()) +
        ggplot2::geom_errorbar(ggplot2::aes(ymin=LOWER_0.95CI, ymax=UPPER_0.95CI), width=.2,
                               position= ggplot2::position_dodge(1)) +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) +
        ggplot2::ggtitle(label = paste("Library complexity extrapolation -", speciesName),
                         sub = paste("Expected distinct reads at", e_reads, "depth. CI = 95%")) +
        ggplot2::ylab("Expected distinct reads") +
        ggplot2::xlab("Samples")
    return(ggOUT)
}

# Library complexity line plots
gg_lceLines <- function(f_tab){
    require(plotly)
    tmp <- lapply(f_tab$lce, fread) %>%
        set_names(f_tab$id)
    tmp <- lapply(names(tmp), function(x){
        cbind(tmp[[x]], id = x)
    }) %>% do.call(what = rbind) %>% data.table()
    tmpGG <- ggplot(tmp) + geom_line(aes(x = TOTAL_READS,
                                         y = EXPECTED_DISTINCT,
                                         colour = id)) +
        geom_ribbon(aes(x = TOTAL_READS,
                        y = EXPECTED_DISTINCT,
                        ymin=LOWER_0.95CI,
                        ymax=UPPER_0.95CI,
                        fill = id),alpha=0.2) + theme_minimal()
    ggplotly(tmpGG)
}

# Notebook 2 ###################################################################

# Function for cytifine persistence to Bisulphite treatment
detect_m5C <- function(DT){
    tmp <- round(DT$C / (DT$T + DT$C), 4)
    tibble::add_column(DT, det_m5C = tmp)
}

# Function for cytidine persistence to Bisulphite treatment
detect_ac4C <- function(DT){
    DT <- data.table(DT)
    tmp <- round(DT$`T` / (DT$`T` + DT$`C`), 6)
    tibble::add_column(DT, det_ac4C = tmp)
}

# Add column of Different Nucleotide to reference ratio, diffToRef and nucTotal columns are required
add_diffNucToRefRatio <- function(DT){
    DT <- tx_add_diffNucToRef(DT) %>% tx_add_nucTotal()
    tmp <- round(DT$diffToRef / DT$nucTotal, 6)
    tibble::add_column(DT, diffToRef_Ratio = tmp)
}

# Start rate
add_StartRate <- function(DT, minCov = 50){
    tmp <- (DT$start_5p + 1) / (DT$cov + 1)
    tmp[DT$cov < minCov] <- NA
    tibble::add_column(DT, start_Ratio = tmp)
}

# Start rate
add_StartRate_1bpDS <- function(DT, minCov = 50){
    tmp <- (DT$start_5p + 1) / (DT$cov + 1)
    tmp[DT$cov < minCov] <- NA
    DTL <- tibble::add_column(DT, startRate_1bpDS = tmp) %>% tx_split_DT()
    DT <- lapply(DTL, function(DT){
        DT$startRate_1bpDS <- c(tail(DT$startRate_1bpDS, -1), NA)
        DT
    }) %>% tx_merge_DT()
    return(DT)
}

# Start rate 1bp DS with pseudoCount
add_StartRate_1bpDS_pc <- function(DT, minCov = 50, pc = 0.01){
    tmp <- (DT$start_5p + pc) / (DT$cov + pc)
    tmp[DT$cov < minCov] <- NA
    DTL <- tibble::add_column(DT, startRate_1bpDS_pc = tmp) %>% tx_split_DT()
    DT <- lapply(DTL, function(DT){
        DT$startRate_1bpDS_pc <- c(tail(DT$startRate_1bpDS_pc, -1), NA)
        DT
    }) %>% tx_merge_DT()
    return(DT)
}

# 2Omescore calculation: Score A
detect_2OmeScore <- function(DT_low, DT_high, neighFlankSize = 5, minCov = 50){
    DT_low <- DT_low %>% add_StartRate(minCov = minCov)
    DT_high <- DT_high %>% add_StartRate()
    if(all(union(DT_low$gene, DT_high$gene) %in% intersect(DT_low$gene, DT_high$gene))){
        iGenes <- intersect(DT_low$gene, DT_high$gene)
    }else{stop("DT_low and DT_high do not share the same genes")}
    DTL1 <- DT_low %>% tx_split_DT()
    DTL2 <- DT_high %>% tx_split_DT()
    OUT <- lapply(iGenes, function(iG){

        RES <- data.table(DTL1[[iG]][,colnames(DTL1[[iG]]) %in% storm_baseCols, with = FALSE])
        tmpSR <- c(tail(DTL1[[iG]]$start_Ratio / DTL2[[iG]]$start_Ratio, -1), NA)
        sc_range <- (1+neighFlankSize):(nrow(DTL1[[iG]])-neighFlankSize)
        score <- rep(NA, nrow(DTL1[[iG]]))
        for(i in sc_range){
            score[i] <- log2(tmpSR[i]) - log2(mean(c(tmpSR[(i - neighFlankSize):(i - 1)],
                                                     tmpSR[(i + 1):(i + neighFlankSize)])))
        }
        RES$twoOme_score <- score
        return(RES)
    })
    do.call(what = rbind, OUT)
}

# Scaling twoOme_score by nucleotide groups
scale_by_nuc <- function(DT){
    iNucs <- DT$refSeq %>% unique()
    for(iN in iNucs){
        selPos <- which(DT$refSeq == iN)
        DT$twoOme_score[selPos] <- scale(DT$twoOme_score[selPos])
    }
    return(DT)
}

# Add RNAmod nucs and pos columns to STORM$RES
add_Nuc <- function(STORM, rib){
    names(rib) <- gsub(pattern = "(Hs_|Sc_)", "", x = names(rib))
    STORM$RES <- tibble::add_column(STORM$RES,
                                    nuc = rib$nuc[match(STORM$RES$pos, rib$pos)],
                                    .after = "refSeq")
    return(STORM)
}

# Add RNAmod nucs and pos columns to STORM$RES
add_posAndNuc <- function(STORM, rib){
    names(rib) <- gsub(pattern = "(Hs_|Sc_)", "", x = names(rib))
    STORM$RES$pos <- paste(STORM$RES$gene, STORM$RES$txcoor, sep = ":")
    STORM$RES$nuc <- rib$nuc[match(STORM$RES$pos, rib$pos)]
    return(STORM)
}

# Score A Birkedal et al., 2015
scoreA <- function(CHR, d = 6, minMedEnv = 15){
    scoreA <- rep(NA, length(CHR)); names(scoreA) <- names(CHR)
    for (P in (d+1):(length(CHR)-d)){
        ePos <- CHR[P]
        leftFlank <- CHR[(P-d):(P-1)]
        rightFlank <- CHR[(P+1):(P+d)]
        if( min(median(c(rightFlank, leftFlank), na.rm = T)) < minMedEnv ){next()}
        m1 <- mean(leftFlank)
        m2 <- mean(rightFlank)
        sd1 <- sd(leftFlank)
        sd2 <- sd(rightFlank)
        scoreA[P] <- 1 - (((2 * ePos) + 1) /
                              ((0.5*abs(m1-sd1)) + ePos + (0.5*abs(m2-sd2)) + 1))
    }
    scoreA[scoreA < 0] <- 0
    return(scoreA)
}

# Adding score A to multiGene DT
add_scoreA_3p <- function(DT, minMedEnv = 15, nCores = 1){
    DTL <- tx_split_DT(DT)
    parallel::mclapply(mc.cores = nCores, DTL, function(x){
        tmp <- scoreA(x$end_3p, minMedEnv = minMedEnv)
        tibble::add_column(x, scoreA_3p = tmp)
    }) %>% tx_merge_DT()
}

# Adding score A 5p_ends to multiGene DT
add_scoreA_5p <- function(DT, minMedEnv = 15, nCores = 1){
    DTL <- tx_split_DT(DT)
    parallel::mclapply(mc.cores = nCores, DTL, function(x){
        if(!("startRate_1bpDS" %in% names(x))){
            x <- add_StartRate_1bpDS(x)
        }
        tmp <- scoreA(x$start_5p, minMedEnv = minMedEnv)
        tmp <- c(tail(tmp, -1), NA)
        tibble::add_column(x, scoreA_5p = tmp)
    }) %>% tx_merge_DT()
}

# STORM functions ##############################################################





# Special Vectors ##############################################################

storm_baseCols <- c("chr", "gencoor", "strand", "gene", "txcoor", "pos",
                    "refSeq", "set", "nuc")
storm_baseCoorCols <- c("chr", "gencoor", "strand", "gene", "txcoor", "pos", "refSeq")

RNAmods_vec <- c("Y", "Nm", "m5C", "ac4C", "m1A", "m7G", "m3U")

# STORM add functions ##########################################################
# Add C->T misincorporation difference (CtoT_MRD) to STORM object
add_CtoT_MRD <- function(STORM, group_A, group_B, newColName = "auto",
                         onNucs = c("C")){
    if(newColName == "auto"){
        newColName <- paste("CtoT.MRD", group_A, group_B, sep = "_")
    }
    sets <- intersect(STORM$META[group == group_A,]$set, STORM$META[group == group_B,]$set)
    if(length(sets) == 0){stop("No sets found with both group_A and group_B label")}
    OUT <- lapply(sets, function(iSet){
        id_A <- STORM$META[set == iSet,][group == group_A]$id
        id_B <- STORM$META[set == iSet,][group == group_B]$id
        if(isEmpty(id_A)){stop(group_A, " not found in STORM object META data")}
        if(isEmpty(id_B)){stop(group_B, " not found in STORM object META data")}
        DT_A <- STORM$DATA[[id_A]] %>% detect_ac4C()
        DT_B <- STORM$DATA[[id_B]] %>% detect_ac4C()
        tmpRES <- data.table(DT_A$det_ac4C - DT_B$det_ac4C)
        tmpRES <- data.table(iSet, newColName, tmpRES)
        names(tmpRES) <- c("set", "metric", "score")
        tmpRES <- data.table(DT_A[,c("chr", "gencoor", "strand", "gene",
                                     "txcoor", "pos", "refSeq")], tmpRES)
        tmpRES[!(tmpRES$refSeq %in% onNucs), "score"] <- NA
        return(tmpRES)
    }) %>% do.call(what = rbind)
    STORM$RES <- hlpr_add_REScols(STORM$RES, OUT)
    return(STORM)
}

# Add SRD 1bp down-stream to STORM object
add_SRD1bpDS <- function(STORM, group_A, group_B, newColName = "auto",
                         onNucs = c("A", "C", "G", "T")){
    if(newColName == "auto"){
        newColName <- paste("SRD1bpDS", group_A, group_B, sep = "_")
    }
    sets <- intersect(STORM$META[group == group_A,]$set, STORM$META[group == group_B,]$set)
    if(length(sets) == 0){stop("No sets found with both group_A and group_B label")}
    OUT <- lapply(sets, function(iSet){
        id_A <- STORM$META[set == iSet,][group == group_A]$id
        id_B <- STORM$META[set == iSet,][group == group_B]$id
        if(isEmpty(id_A)){stop(group_A, " not found in STORM object META data")}
        if(isEmpty(id_B)){stop(group_B, " not found in STORM object META data")}
        DT_A <- add_StartRate_1bpDS(STORM$DATA[[id_A]])
        DT_B <- add_StartRate_1bpDS(STORM$DATA[[id_B]])
        tmpRES <- data.table(DT_A$startRate_1bpDS - DT_B$startRate_1bpDS)
        tmpRES <- data.table(iSet, newColName, tmpRES)
        names(tmpRES) <- c("set", "metric", "score")
        tmpRES <- data.table(DT_A[,c("chr", "gencoor", "strand", "gene",
                                     "txcoor", "pos", "refSeq")], tmpRES)
        tmpRES[!(tmpRES$refSeq %in% onNucs), "score"] <- NA
        return(tmpRES)
    }) %>% do.call(what = rbind)
    STORM$RES <- hlpr_add_REScols(STORM$RES, OUT)
    return(STORM)
}

# Add Start rate log2-Fold-Change 1bp down-stream to STORM object
add_SRlog2FCh1bpDS <- function(STORM, group_A, group_B, newColName = "auto",
                               onNucs = c("A", "C", "G", "T"), minCov = 50){
    if(newColName == "auto"){
        newColName <- paste("SRlog2FCh1bpDS", group_A, group_B, sep = "_")
    }
    sets <- intersect(STORM$META[group == group_A,]$set, STORM$META[group == group_B,]$set)
    if(length(sets) == 0){stop("No sets found with both group_A and group_B label")}
    OUT <- lapply(sets, function(iSet){
        id_A <- STORM$META[set == iSet,][group == group_A]$id
        id_B <- STORM$META[set == iSet,][group == group_B]$id
        if(isEmpty(id_A)){stop(group_A, " not found in STORM object META data")}
        if(isEmpty(id_B)){stop(group_B, " not found in STORM object META data")}
        DT_A <- add_StartRate_1bpDS(STORM$DATA[[id_A]], minCov = minCov)
        DT_B <- add_StartRate_1bpDS(STORM$DATA[[id_B]], minCov = minCov)
        tmpRES <- data.table(log2(DT_A$startRate_1bpDS / DT_B$startRate_1bpDS))
        tmpRES <- data.table(iSet, newColName, tmpRES)
        names(tmpRES) <- c("set", "metric", "score")
        tmpRES <- data.table(DT_A[,c("chr", "gencoor", "strand", "gene",
                                     "txcoor", "pos", "refSeq")], tmpRES)
        tmpRES[!(tmpRES$refSeq %in% onNucs), "score"] <- NA
        return(tmpRES)
    }) %>% do.call(what = rbind)
    STORM$RES <- hlpr_add_REScols(STORM$RES, OUT)
    return(STORM)
}

# Add SRD to STORM object
add_SRD <- function(STORM, group_A, group_B, newColName = "auto",
                    onNucs = c("A", "C", "G", "T")){
    if(newColName == "auto"){
        newColName <- paste("SRD", group_A, group_B, sep = "_")
    }
    sets <- intersect(STORM$META[group == group_A,]$set, STORM$META[group == group_B,]$set)
    if(length(sets) == 0){stop("No sets found with both group_A and group_B label")}
    OUT <- lapply(sets, function(iSet){
        id_A <- STORM$META[set == iSet,][group == group_A]$id
        id_B <- STORM$META[set == iSet,][group == group_B]$id
        if(isEmpty(id_A)){stop(group_A, " not found in STORM object META data")}
        if(isEmpty(id_B)){stop(group_B, " not found in STORM object META data")}
        DT_A <- add_StartRate(STORM$DATA[[id_A]])
        DT_B <- add_StartRate(STORM$DATA[[id_B]])
        tmpRES <- data.table(DT_A$startRatio - DT_B$startRatio)
        tmpRES <- data.table(iSet, newColName, tmpRES)
        names(tmpRES) <- c("set", "metric", "score")
        tmpRES <- data.table(DT_A[,c("chr", "gencoor", "strand", "gene",
                                     "txcoor", "pos", "refSeq")], tmpRES)
        tmpRES[!(tmpRES$refSeq %in% onNucs), "score"] <- NA
        return(tmpRES)
    }) %>% do.call(what = rbind)
    STORM$RES <- hlpr_add_REScols(STORM$RES, OUT)
    return(STORM)
}

# Add 2-O-methylation score to STORM-seq object
add_2OmeScore <- function(STORM,
                          lib_LowdNTPs,
                          lib_HighdNTPs,
                          newColName = "auto",
                          neighFlankSize = 5,
                          perNuc_Zscaling = F,
                          onNucs = c("A", "C", "G", "T")){
    if(newColName == "auto"){
        newColName <- paste("NmStopScore", lib_LowdNTPs, lib_HighdNTPs, sep = "_")
    }
    sets <- intersect(STORM$META[group == lib_LowdNTPs,]$set, STORM$META[group == lib_HighdNTPs,]$set)
    if(length(sets) == 0){stop("No sets found with both lib_LowdNTPs and lib_HighdNTPs label")}
    OUT <- lapply(sets, function(iSet){
        id_A <- STORM$META[set == iSet,][group == lib_LowdNTPs]$id
        id_B <- STORM$META[set == iSet,][group == lib_HighdNTPs]$id
        if(isEmpty(id_A)){stop(group_A, " not found in STORM object META data")}
        if(isEmpty(id_B)){stop(group_B, " not found in STORM object META data")}
        DT <- detect_2OmeScore(DT_low = STORM$DATA[[id_A]],
                               DT_high = STORM$DATA[[id_B]],
                               neighFlankSize = neighFlankSize)
        tmpRES <- data.table(twoOme_score = DT$twoOme_score)
        tmpRES <- data.table(iSet, newColName, tmpRES)
        names(tmpRES) <- c("set", "metric", "score")
        tmpRES <- data.table(DT[,c("chr", "gencoor", "strand", "gene",
                                   "txcoor", "pos", "refSeq")], tmpRES)
        tmpRES[!(tmpRES$refSeq %in% onNucs), "score"] <- NA
        return(tmpRES)
    })
    if(perNuc_Zscaling == T){
        OUT <- lapply(OUT, function(DT){
            iNucs <- DT$refSeq %>% unique()
            for(iN in iNucs){
                selPos <- which(DT$refSeq == iN)
                DT$score[selPos] <- scale(DT$score[selPos])
            }
            DT
        })
    }
    OUT <- do.call(OUT, what = rbind)
    STORM$RES <- hlpr_add_REScols(STORM$RES, OUT)
    return(STORM)
}

# Add Misincorporation Rate Difference (MRD) treatA-treatB to STORM object
add_MRD <- function(STORM, group_A, group_B, newColName = "auto",
                    onNucs = c("A", "C", "G", "T")){
    if(newColName == "auto"){
        newColName <- paste("MRD", group_A, group_B, sep = "_")
    }
    sets <- intersect(STORM$META[group == group_A,]$set, STORM$META[group == group_B,]$set)
    if(length(sets) == 0){stop("No sets found with both group_A and group_B label")}
    OUT <- lapply(sets, function(iSet){
        id_A <- STORM$META[set == iSet,][group == group_A]$id
        id_B <- STORM$META[set == iSet,][group == group_B]$id
        if(isEmpty(id_A)){stop(group_A, " not found in STORM object META data")}
        if(isEmpty(id_B)){stop(group_B, " not found in STORM object META data")}
        DT_A <- STORM$DATA[[id_A]] %>% add_diffNucToRefRatio()
        DT_B <- STORM$DATA[[id_B]] %>% add_diffNucToRefRatio()
        tmpRES <- data.table(DT_A$diffToRef_Ratio - DT_B$diffToRef_Ratio)
        tmpRES <- data.table(iSet, newColName, tmpRES)
        names(tmpRES) <- c("set", "metric", "score")
        tmpRES <- data.table(DT_A[,c("chr", "gencoor", "strand", "gene",
                                     "txcoor", "pos", "refSeq")], tmpRES)
        tmpRES[!(tmpRES$refSeq %in% onNucs), "score"] <- NA
        return(tmpRES)
    }) %>% do.call(what = rbind)
    STORM$RES <- hlpr_add_REScols(STORM$RES, OUT)
    return(STORM)
}

# Add cytidine persistance index to STORM object
add_CytPer <- function(STORM, group_A, group_B, newColName = "auto",
                       onNucs = c("C")){
    if(newColName == "auto"){
        newColName <- paste("CytPer", group_A, group_B, sep = "_")
    }
    sets <- intersect(STORM$META[group == group_A,]$set, STORM$META[group == group_B,]$set)
    if(length(sets) == 0){stop("No sets found with both group_A and group_B label")}
    OUT <- lapply(sets, function(iSet){
        id_A <- STORM$META[set == iSet,][group == group_A]$id
        id_B <- STORM$META[set == iSet,][group == group_B]$id
        if(isEmpty(id_A)){stop(group_A, " not found in STORM object META data")}
        if(isEmpty(id_B)){stop(group_B, " not found in STORM object META data")}
        DT_A <- STORM$DATA[[id_A]] %>% detect_m5C()
        DT_B <- STORM$DATA[[id_B]] %>% detect_m5C()
        tmpRES <- data.table(DT_A$det_m5C - (1 - DT_B$det_m5C))
        tmpRES <- data.table(iSet, newColName, tmpRES)
        names(tmpRES) <- c("set", "metric", "score")
        tmpRES <- data.table(DT_A[,c("chr", "gencoor", "strand", "gene",
                                     "txcoor", "pos", "refSeq")], tmpRES)
        tmpRES[!(tmpRES$refSeq %in% onNucs), "score"] <- NA
        return(tmpRES)
    }) %>% do.call(what = rbind)
    STORM$RES <- hlpr_add_REScols(STORM$RES, OUT)
    return(STORM)
}

# Add start rate column in STORM$RES
add_SR_1bpDS <- function(STORM, group_A, newColName = "auto", onNucs = c("A", "C", "G", "T")){
    if(newColName == "auto"){
        newColName <- paste("StartRate1bpDS", group_A, sep = "_")
    }
    sets <- STORM$META[group == group_A,]$set
    if(length(sets) == 0){stop("No sets found with both group_A and group_B label")}
    OUT <- lapply(sets, function(iSet){
        id_A <- STORM$META[set == iSet & group == group_A]$id
        if(isEmpty(id_A)){stop(group_A, " not found in STORM object META data")}
        DT_A <- add_StartRate_1bpDS(STORM$DATA[[id_A]])
        tmpRES <- data.table(DT_A$startRate_1bpDS)
        tmpRES <- data.table(iSet, newColName, tmpRES)
        names(tmpRES) <- c("set", "metric", "score")
        tmpRES <- data.table(DT_A[,c("chr", "gencoor", "strand", "gene",
                                     "txcoor", "pos", "refSeq")], tmpRES)
        tmpRES[!(tmpRES$refSeq %in% onNucs), "score"] <- NA
        return(tmpRES)
    }) %>% do.call(what = rbind)
    STORM$RES <- hlpr_add_REScols(STORM$RES, OUT)
    return(STORM)
}

# Add start rate column in STORM$RES
add_MR <- function(STORM, group_A, newColName = "auto", onNucs = c("A", "C", "G", "T")){
    if(newColName == "auto"){
        newColName <- paste("MissIncRate", group_A, sep = "_")
    }
    sets <- STORM$META[group == group_A,]$set
    if(length(sets) == 0){stop("No sets found with both group_A and group_B label")}
    OUT <- lapply(sets, function(iSet){
        id_A <- STORM$META[set == iSet & group == group_A]$id
        if(isEmpty(id_A)){stop(group_A, " not found in STORM object META data")}
        DT_A <- add_diffNucToRefRatio(STORM$DATA[[id_A]])
        tmpRES <- data.table(DT_A$diffToRef_Ratio)
        tmpRES <- data.table(iSet, newColName, tmpRES)
        names(tmpRES) <- c("set", "metric", "score")
        tmpRES <- data.table(DT_A[,c("chr", "gencoor", "strand", "gene",
                                     "txcoor", "pos", "refSeq")], tmpRES)
        tmpRES[!(tmpRES$refSeq %in% onNucs), "score"] <- NA
        return(tmpRES)
    }) %>% do.call(what = rbind)
    STORM$RES <- hlpr_add_REScols(STORM$RES, OUT)
    return(STORM)
}

# Add ScoreA 3prime
storm_add_scoreA3p <- function(STORM, group_A, newColName = "auto",
                               onNucs = c("A", "C", "G", "T"), minMedCov = 15){
    if(newColName == "auto"){
        newColName <- paste("ScoreA3p", group_A, sep = "_")
    }
    sets <- STORM$META[group == group_A,]$set
    if(length(sets) == 0){stop("No sets found with group_A label")}
    OUT <- lapply(sets, function(iSet){
        id_A <- STORM$META[set == iSet & group == group_A]$id
        DT_A <- add_scoreA_3p(STORM$DATA[[id_A]], minMedEnv = minMedCov)
        if(isEmpty(id_A)){stop(group_A, " not found in STORM object META data")}
        tmpRES <- data.table(DT_A$scoreA_3p)
        tmpRES <- data.table(iSet, newColName, tmpRES)
        names(tmpRES) <- c("set", "metric", "score")
        tmpRES <- data.table(DT_A[,c("chr", "gencoor", "strand", "gene",
                                     "txcoor", "pos", "refSeq")], tmpRES)
        tmpRES[!(tmpRES$refSeq %in% onNucs), "score"] <- NA
        return(tmpRES)
    }) %>% do.call(what = rbind)
    STORM$RES <- hlpr_add_REScols(STORM$RES, OUT)
    return(STORM)
}

# Add ScoreA 5prime
storm_add_scoreA5p <- function(STORM, group_A, newColName = "auto",
                               onNucs = c("A", "C", "G", "T"), minMedCov = 15){
    if(newColName == "auto"){
        newColName <- paste("ScoreA5p", group_A, sep = "_")
    }
    sets <- STORM$META[group == group_A,]$set
    if(length(sets) == 0){stop("No sets found with both group_A and group_B label")}
    OUT <- lapply(sets, function(iSet){
        id_A <- STORM$META[set == iSet & group == group_A]$id
        if(isEmpty(id_A)){stop(group_A, " not found in STORM object META data")}
        DT_A <- add_scoreA_5p(STORM$DATA[[id_A]], minMedEnv = minMedCov)
        tmpRES <- data.table(DT_A$scoreA_5p)
        tmpRES <- data.table(iSet, newColName, tmpRES)
        names(tmpRES) <- c("set", "metric", "score")
        tmpRES <- data.table(DT_A[, c("chr", "gencoor", "strand", "gene",
                                      "txcoor", "pos", "refSeq")], tmpRES)
        tmpRES[!(tmpRES$refSeq %in% onNucs), "score"] <- NA
        return(tmpRES)
    }) %>% do.call(what = rbind)
    STORM$RES <- hlpr_add_REScols(STORM$RES, OUT)
    return(STORM)
}

# Notebook 3 ###################################################################

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

SSIIIandTGIRT_metrics <- function(STORM){
    STORM %>%
        add_Y_metrics %>%
        add_Nm_metrics %>%
        add_ac4C_metrics %>%
        add_m1A_metrics %>%
        add_m7G_metrics %>%
        add_m5C_metrics
}

add_Y_metrics <- function(STORM){
    STORM %>%
        add_SRD1bpDS("CMC.TGIRT", "Mock.TGIRT") %>%
        add_SRD1bpDS("CMC.SSIII", "Mock.SSIII") %>%
        add_SRlog2FCh1bpDS("CMC.TGIRT", "Mock.TGIRT") %>%
        add_SRlog2FCh1bpDS("CMC.SSIII", "Mock.SSIII")
}

Y_scores <- c("SRD1bpDS_CMC.TGIRT_Mock.TGIRT",
              "SRD1bpDS_CMC.SSIII_Mock.SSIII",
              "SRlog2FCh1bpDS_CMC.TGIRT_Mock.TGIRT",
              "SRlog2FCh1bpDS_CMC.SSIII_Mock.SSIII")

add_Nm_metrics <- function(STORM){
    STORM %>%
        storm_add_scoreA3p("Mock.TGIRT", minMedCov = 50) %>%
        storm_add_scoreA3p("Mock.SSIII", minMedCov = 50) %>%
        add_2OmeScore("MocklowdNTPs.SSIII", "Mock.SSIII", perNuc_Zscaling = T)
}

Nm_scores <- c("NmStopScore_MocklowdNTPs.SSIII_Mock.SSIII",
               "ScoreA3p_Mock.TGIRT",
               "ScoreA3p_Mock.SSIII")

add_ac4C_metrics <- function(STORM){
    STORM %>%
        add_CtoT_MRD("Ac4C.TGIRT", "Mock.TGIRT") %>%
        add_CtoT_MRD("Ac4C.TGIRT", "DeacetylatedAc4C.TGIRT") %>%
        add_CtoT_MRD("Ac4C.SSIII", "Mock.SSIII") %>%
        add_CtoT_MRD("Ac4C.SSIII", "DeacetylatedAc4C.SSIII")
}

ac4C_scores <- c("CtoT.MRD_Ac4C.TGIRT_Mock.TGIRT",
                 "CtoT.MRD_Ac4C.TGIRT_DeacetylatedAc4C.TGIRT",
                 "CtoT.MRD_Ac4C.SSIII_Mock.SSIII",
                 "CtoT.MRD_Ac4C.SSIII_DeacetylatedAc4C.SSIII")

add_m1A_metrics <- function(STORM){
    STORM %>%
        add_SRD1bpDS("Mock.SSIII", "Dimroth.SSIII") %>%
        add_SRlog2FCh1bpDS("Mock.SSIII", "Dimroth.SSIII") %>%
        add_MRD("Mock.TGIRT", "Dimroth.TGIRT") %>%
        add_SRD1bpDS("Mock.TGIRT", "Dimroth.TGIRT") %>%
        add_SRlog2FCh1bpDS("Mock.TGIRT", "Dimroth.TGIRT")
}

m1A_scores <- c("SRD1bpDS_Mock_Dimroth",
                "MRD_Mock_Dimroth",
                "MRD_Mock.TGIRT_Dimroth.TGIRT",
                "SRD1bpDS_Mock.SSIII_Dimroth.SSIII")

add_m7G_metrics <- function(STORM){
    STORM %>%
        add_SRlog2FCh1bpDS("DeacetylatedAc4C.TGIRT", "Ac4C.TGIRT") %>%
        add_SRlog2FCh1bpDS("DeacetylatedAc4C.SSIII", "Ac4C.SSIII") %>%
        add_SRD1bpDS("DeacetylatedAc4C.TGIRT", "Ac4C.TGIRT") %>%
        add_SRD1bpDS("DeacetylatedAc4C.SSIII", "Ac4C.SSIII") %>%
        add_MRD("DeacetylatedAc4C.TGIRT", "Ac4C.TGIRT") %>%
        add_MRD("DeacetylatedAc4C.SSIII", "Ac4C.SSIII")
}

m7G_scores <- c("MRD_DeacetylatedAc4C.TGIRT_Ac4C.TGIRT",
                "MRD_DeacetylatedAc4C.SSIII_Ac4C.SSIII",
                "SRD1bpDS_DeacetylatedAc4C.SSIII_Ac4C.SSIII")

add_m5C_metrics <- function(STORM){
    STORM %>%
        add_CytPer("m5C.TGIRT", "Mock.TGIRT")
}

m5C_scores <- c("CytPer_m5C_Mock",
                "CytPer_m5C.TGIRT_Mock.TGIRT")

# Scores per RNA modification
add_m3U_metrics <- function(STORM){
    STORM %>%
        add_MRD("Mock", "AlkBmix")
}

m3U_scores <- c("MRD_Mock_AlkBmix")

# Boxplot of $RES calculated metrics grouping by modified nucleotides
ggMetricsNuc <- function(STORM, title){
    tmpDT <- STORM$RES %>% na.omit
    ggplot(tmpDT, aes(x = nuc, y = score)) + geom_boxplot(outlier.colour = NA) +
        geom_point(aes(colour = metric), alpha = 0.2) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        facet_wrap(~metric, scales = "free") +
        ggtitle(title)
}

# Plot scatterplot of scores colored by gene by tx coordinate
ggMetricsPos <- function(STORM, title){
    tmpDT <- STORM$RES %>% na.omit
    ggplot(tmpDT) +
        geom_point(aes(x = pos, y = score, colour = gene), alpha = 0.2) +
        facet_grid(metric~set, scales = "free") +
        theme_bw() +
        theme(axis.text.x=element_blank()) +
        ggtitle(title) + xlab("txCoor")
}

# Initializa STORM$CALLS
hlp_start_CALLS <- function(STORM){
    STORM$CALLS <- lapply(seq_along(unique(STORM$META$set)), function(x) {
        lapply(RNAmods_vec, function(x) NULL) %>% set_names(RNAmods_vec)
    })
    names(STORM$CALLS) <- unique(STORM$META$set)
    STORM
}

# Assign scores to respective RNAmods in STORM$CALLS
hlp_assign_scores <- function(STORM, RNAmod, scores){
    tmp <- STORM$RES %>% pivot_wider(names_from = metric, values_from = score) %>% data.table
    selVars <- c(storm_baseCols, scores)
    tmp <- tmp[,names(tmp) %in% selVars, with = FALSE]
    tmp <- split(tmp, tmp$set)
    for(iN in names(tmp)){
        STORM$CALLS[[iN]][[RNAmod]] <- tmp[[iN]]
    }
    STORM
}

# Abbreviated STORM$CALLS
storm_makeCalls <- function(STORM){
    hlp_start_CALLS(STORM) %>%
        hlp_assign_scores("Y", Y_scores) %>%
        hlp_assign_scores("Nm", Nm_scores) %>%
        hlp_assign_scores("m5C", m5C_scores) %>%
        hlp_assign_scores("ac4C", ac4C_scores) %>%
        hlp_assign_scores("m1A", m1A_scores) %>%
        hlp_assign_scores("m7G", m7G_scores) %>%
        hlp_assign_scores("m3U", m3U_scores)
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

# Predict class based on cutpointr model
predict_cutpointr <- function(cp_model, newdata){
    pr <- cp_model$predictor
    dir <- cp_model$direction
    ocp <- cp_model$optimal_cutpoint
    cuts <- lapply(seq_along(pr), function(x){
        if(dir[x] == "<="){
            newdata[,pr[x]] <= ocp[x]
        }else if(dir[x] == ">="){
            newdata[,pr[x]] >= ocp[x]
        }else{stop("Direction of predictor is not supported")}
    }) %>% do.call(what = "rbind") %>% colSums() %>% equals(length(pr))
    ifelse(cuts, as.character(cp_model$pos_class[1]), as.character(cp_model$neg_class[1]))
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

# pheatmap, no clustering
pheat_NoClust <- function(x, ...){
    pheatmap::pheatmap(x, cluster_cols = FALSE, cluster_rows = FALSE, ...)
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
        logMod <- glm(isMod ~ ., family = binomial(), data = tData)
        CALLS$logPred <- stats::predict(logMod, newdata = CALLS, type = "response")
        sapply(thresholds, function(x) MCC(CALLS$logPred, CALLS$isMod, x))
    }) %>% do.call(what = cbind) %>% set_rownames(thresholds) %>% set_colnames(combNames)
    selComb <- which.max(colMedians(OUT))
    selThr <- which.max(OUT[,selComb])
    # FinalModel
    tData <- CALLS[, c(allCombVar[[selComb]], "isMod")]
    tData$isMod <- as.numeric(tData$isMod)
    logMod <- glm(isMod ~ ., family = binomial(), data = tData)
    return(list(logiMod = logMod, thr = thresholds[selThr], vars = allCombVar[[selComb]],
                MCCmat = OUT, MCC= OUT[selThr, selComb]))
}

# Update STORM$CALLS with logistic and linear scores and prediction
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
        }) %>% set_names(names(STORM$CALLS[[callList]]))
        return(tmp)
    }) %>% set_names(names(STORM$CALLS))
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

# Train RF and CP models
train_RF_CP <- function(STORM, RNAmodsList, varList){
    #RandomFores models
    RES <- STORM$RES
    nucList_Sc <- lapply(RNAmodsList, function(i) RES$nuc %in% i)
    CP_models <- lapply(seq_along(RNAmodsList), function(i){
        selVars <- lapply(varList[[i]], function(i){
            grep(pattern = i, x = names(RES), value = TRUE)}) %>% unlist
        tmpData <- data.frame(set_colnames(RES[, selVars, with = F], varList[[i]]),
                              RNAmod = factor(nucList_Sc[[i]]))
        CPmodel <- lapply(varList[[i]], function(iVar){
            tmpC <- cutpointr(tmpData[, iVar], tmpData[, "RNAmod"],
                              method = maximize_metric, metric = sum_sens_spec,
                              direction = ">=", pos_class = TRUE,
                              use_midpoints = TRUE, na.rm = TRUE)
            tmpC$predictor <- iVar
            tmpC
        })
        CPmodel
    })
    RES[is.na(RES)] <- 0
    RF_models <- lapply(seq_along(varList), function(i){
        selVars <- lapply(varList[[i]], function(i) grep(pattern = i, x = names(RES), value = TRUE)) %>% unlist
        tmpDat <- RES[, selVars, with = F] %>% set_colnames(varList[[i]])
        randomForest(x = tmpDat, y = factor(nucList_Sc[[i]]))
    })
    #CutPointer
    list(RNAmodList = RNAmodsList, RF_models = RF_models, CP_models = CP_models)
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
                as.logical(predict(CPmod, tmpDat))
            }) %>% do.call(what = data.frame) %>% unname %>% apply(1, all)
        }) %>% do.call(what = cbind) %>% data.table() %>%
            set_names(paste("pred", RNAmods_vec[i], iSets, sep = "_"))

        scor_RES <- lapply(seq_along(iSets), function(j){
            lapply(CP_models[[i]], function(CPmod){
                tmpPat <- paste0("^", gsub("_TGIRT", "", CPmod$predictor))
                selVar <- grepl(pattern = tmpPat, names(STORM$RES)) &
                    grepl(pattern = iSets[j], names(STORM$RES))
                STORM$RES[,which(selVar), with = FALSE]
            }) %>% do.call(what = cbind)
        }) %>% do.call(what = cbind) %>% data.table()

        tree_RES <- lapply(seq_along(iSets), function(j){
            RF_mod <- RF_list[[i]]
            RF_vars <- names(RF_mod$forest$xlevels)
            newDat <- lapply(RF_vars, function(RF_v){
                tmpP <- gsub(pattern = "_TGIRT", replacement = "", x = RF_v)
                selVar <- grepl(pattern = tmpP, x = names(STORM$RES))
                selSmp <- grepl(pattern = iSets[j], x = names(STORM$RES)) & selVar
                newDat <- data.frame(x = STORM$RES[, selSmp, with = FALSE]) %>% set_names(RF_v)
            }) %>% do.call(what = cbind)
            newDat[is.na(newDat)] <- 0
            out <- predict(RF_mod, newDat, norm.votes = TRUE, type = "vote")
            unname(out[,"TRUE"])
        }) %>% do.call(what = cbind) %>% data.table() %>%
            set_names(paste("votes", RNAmods_vec[i], iSets, sep = "_"))

        data.table(coorSys, pred_RES, tree_RES, scor_RES)
    }) %>% set_names(RNAmods_vec)
    return(STORM)
}

