# PREDEFINED METRICS ###########################################################

#' Add all default metrics for STORMseq v0.0.2
#'
#' @param STORM list. STORM object as output by \code{\link{storm_STORM}}
#'
#' @return
#' @export
#'
add_default_metrics_v2 <- function(STORM){
    STORM %>%
        add_2OmeScore("MocklowdNTPs.SSIII", "Mock.SSIII", perNuc_Zscaling = T) %>%
        add_CtoT_MRD("Ac4C.SSIII", "DeacetylatedAc4C.SSIII") %>%
        add_CtoT_MRD("Ac4C.SSIII", "Mock.SSIII") %>%
        add_CtoT_MRD("Ac4C.TGIRT", "DeacetylatedAc4C.TGIRT") %>%
        add_CtoT_MRD("Ac4C.TGIRT", "Mock.TGIRT") %>%
        add_CytPer("m5C.TGIRT", "Mock.TGIRT") %>%
        add_CytPer("RBSseqHeatMg.SSIII", "Mock.SSIII") %>%
        add_CytPer("RBSseqHeatMg.TGIRT", "Mock.TGIRT") %>%
        add_DRD("m5C.TGIRT", "Mock.TGIRT") %>%
        add_DRD("RBSseqHeatMg.SSIII", "Mock.SSIII") %>%
        add_DRD("RBSseqHeatMg.TGIRT", "Mock.TGIRT") %>%
        add_MRD("DeacetylatedAc4C.SSIII", "Ac4C.SSIII") %>%
        add_MRD("DeacetylatedAc4C.TGIRT", "Ac4C.TGIRT") %>%
        add_MRD("Mock.SSIII", "RBSseqHeatMg.SSIII") %>%
        add_MRD("Mock.TGIRT", "AlkBmix.TGIRT") %>%
        add_MRD("Mock.TGIRT", "Dimroth.TGIRT") %>%
        add_MRD("Mock.TGIRT", "RBSseqHeatMg.TGIRT") %>%
        add_scoreA3p("Mock.RTHIV") %>%
        add_scoreA3p("Mock.SSIII") %>%
        add_scoreA3p("Mock.TGIRT") %>%
        add_MRD("NaBH4HydBiotin.RTHIV", "Mock.RTHIV") %>%
        add_MRD("NaBH4HydBiotin.TGIRT", "Mock.TGIRT") %>%
        add_SRD1bpDS("CMC.SSIII", "Mock.SSIII") %>%
        add_SRD1bpDS("CMC.TGIRT", "Mock.TGIRT") %>%
        add_SRD1bpDS("DeacetylatedAc4C.SSIII", "Ac4C.SSIII") %>%
        add_SRD1bpDS("DeacetylatedAc4C.TGIRT", "Ac4C.TGIRT") %>%
        add_SRD1bpDS("Mock.SSIII", "Dimroth.SSIII") %>%
        add_SRlog2FCh1bpDS("CMC.SSIII", "Mock.SSIII") %>%
        add_SRlog2FCh1bpDS("CMC.TGIRT", "Mock.TGIRT") %>%
        add_SRlog2FCh1bpDS("DeacetylatedAc4C.SSIII", "Ac4C.SSIII") %>%
        add_SRlog2FCh1bpDS("DeacetylatedAc4C.TGIRT", "Ac4C.TGIRT") %>%
        add_SRlog2FCh1bpDS("Mock.SSIII", "Dimroth.SSIII")
}

#' Add all default metrics for STORMseq v0.0.1
#'
#' @param STORM list. STORM object as output by \code{\link{storm_STORM}}
#'
#' @return
#' @export
#'
#' @examples
add_default_metrics_v1 <- function(STORM){
    STORM %>%
        add_SRD1bpDS("CMC.TGIRT", "Mock.TGIRT") %>%
        add_SRD1bpDS("CMC.SSIII", "Mock.SSIII") %>%
        add_SRlog2FCh1bpDS("CMC.TGIRT", "Mock.TGIRT") %>%
        add_SRlog2FCh1bpDS("CMC.SSIII", "Mock.SSIII") %>%
        add_scoreA3p("Mock.TGIRT", minMedCov = 50) %>%
        add_scoreA3p("Mock.SSIII", minMedCov = 50) %>%
        add_2OmeScore("MocklowdNTPs.SSIII", "Mock.SSIII", perNuc_Zscaling = T) %>%
        add_CtoT_MRD("Ac4C.TGIRT", "Mock.TGIRT") %>%
        add_CtoT_MRD("Ac4C.TGIRT", "DeacetylatedAc4C.TGIRT") %>%
        add_CtoT_MRD("Ac4C.SSIII", "Mock.SSIII") %>%
        add_CtoT_MRD("Ac4C.SSIII", "DeacetylatedAc4C.SSIII") %>%
        add_SRD1bpDS("Mock.SSIII", "Dimroth.SSIII") %>%
        add_SRlog2FCh1bpDS("Mock.SSIII", "Dimroth.SSIII") %>%
        add_MRD("Mock.TGIRT", "Dimroth.TGIRT") %>%
        add_SRD1bpDS("Mock.TGIRT", "Dimroth.TGIRT") %>%
        add_SRlog2FCh1bpDS("Mock.TGIRT", "Dimroth.TGIRT") %>%
        add_SRlog2FCh1bpDS("DeacetylatedAc4C.TGIRT", "Ac4C.TGIRT") %>%
        add_SRlog2FCh1bpDS("DeacetylatedAc4C.SSIII", "Ac4C.SSIII") %>%
        add_SRD1bpDS("DeacetylatedAc4C.TGIRT", "Ac4C.TGIRT") %>%
        add_SRD1bpDS("DeacetylatedAc4C.SSIII", "Ac4C.SSIII") %>%
        add_MRD("DeacetylatedAc4C.TGIRT", "Ac4C.TGIRT") %>%
        add_MRD("DeacetylatedAc4C.SSIII", "Ac4C.SSIII") %>%
        add_CytPer("m5C.TGIRT", "Mock.TGIRT") %>%
        add_MRD("Mock.TGIRT", "AlkBmix.TGIRT")
}

add_Y_metrics <- function(STORM){
    STORM %>%
        add_SRD1bpDS("CMC.TGIRT", "Mock.TGIRT") %>%
        add_SRD1bpDS("CMC.SSIII", "Mock.SSIII") %>%
        add_SRlog2FCh1bpDS("CMC.TGIRT", "Mock.TGIRT") %>%
        add_SRlog2FCh1bpDS("CMC.SSIII", "Mock.SSIII") %>%
        add_DRD("RBSseqHeatMg.TGIRT", "Mock.TGIRT") %>%
        add_DRD("RBSseqHeatMg.SSIII", "Mock.SSIII") %>%
        add_DRD("m5C.TGIRT", "Mock.TGIRT")
}

add_Nm_metrics <- function(STORM){
    STORM %>%
        add_scoreA3p("Mock.TGIRT") %>%
        add_scoreA3p("Mock.SSIII") %>%
        add_scoreA3p("Mock.RTHIV") %>%
        add_2OmeScore("MocklowdNTPs.SSIII", "Mock.SSIII", perNuc_Zscaling = T)
}

add_ac4C_metrics <- function(STORM){
    STORM %>%
        add_CtoT_MRD("Ac4C.TGIRT", "Mock.TGIRT") %>%
        add_CtoT_MRD("Ac4C.SSIII", "Mock.SSIII") %>%
        add_CtoT_MRD("Ac4C.TGIRT", "DeacetylatedAc4C.TGIRT") %>%
        add_CtoT_MRD("Ac4C.SSIII", "DeacetylatedAc4C.SSIII")
}

add_m1A_metrics <- function(STORM){
    STORM %>%
        add_SRD1bpDS("Mock.SSIII", "Dimroth.SSIII") %>%
        add_SRlog2FCh1bpDS("Mock.SSIII", "Dimroth.SSIII") %>%
        add_MRD("Mock.TGIRT", "Dimroth.TGIRT") %>%
        add_MRD("Mock.TGIRT", "RBSseqHeatMg.TGIRT") %>%
        add_MRD("Mock.SSIII", "RBSseqHeatMg.SSIII")
}

add_m7G_metrics <- function(STORM){
    STORM %>%
        add_SRlog2FCh1bpDS("DeacetylatedAc4C.TGIRT", "Ac4C.TGIRT") %>%
        add_SRlog2FCh1bpDS("DeacetylatedAc4C.SSIII", "Ac4C.SSIII") %>%
        add_SRD1bpDS("DeacetylatedAc4C.TGIRT", "Ac4C.TGIRT") %>%
        add_SRD1bpDS("DeacetylatedAc4C.SSIII", "Ac4C.SSIII") %>%
        add_MRD("DeacetylatedAc4C.TGIRT", "Ac4C.TGIRT") %>%
        add_MRD("DeacetylatedAc4C.SSIII", "Ac4C.SSIII") %>%
        add_MRD("NaBH4HydBiotin.TGIRT", "Mock.TGIRT") %>%
        add_MRD("NaBH4HydBiotin.RTHIV", "Mock.RTHIV")
}

add_m5C_metrics <- function(STORM){
    STORM %>%
        add_CytPer("m5C.TGIRT", "Mock.TGIRT") %>%
        add_CytPer("RBSseqHeatMg.TGIRT", "Mock.TGIRT") %>%
        add_CytPer("RBSseqHeatMg.SSIII", "Mock.SSIII")
}

add_m3U_metrics <- function(STORM){
    STORM %>%
        add_MRD("Mock.TGIRT", "AlkBmix.TGIRT")
}

# METRIC ASSIGNMENT FUNCTIONS ##################################################

#' Initialize STORM$CALLS
#'
#' @param STORM list. STORM object as output by \code{\link{storm_STORM}}
#' @param RNAmods character Vector of RNA modification symbols to assign scores
#' to. By default = brainSTORM::RNAmods_vec.
#'
#' @return
#' @export
#'
#' @examples
hlp_start_CALLS <- function(STORM, RNAmods = RNAmods_vec){
    STORM$CALLS <- lapply(seq_along(unique(STORM$META$set)), function(x) {
        lapply(RNAmods, function(x) NULL) %>% magrittr::set_names(RNAmods)
    })
    names(STORM$CALLS) <- unique(STORM$META$set)
    STORM
}

#' Assign scores to respective RNAmods in STORM$CALLS
#'
#' @param STORM list. STORM object as output by \code{\link{storm_STORM}}
#' @param RNAmod
#' @param scores
#'
#' @return
#' @export
#'
#' @examples
hlp_assign_scores <- function(STORM, RNAmod, scores){
    tmp <- STORM$RES %>%
        tidyr::pivot_wider(names_from = metric, values_from = score) %>%
        data.table::data.table()
    selVars <- c(storm_baseCols, scores)
    tmp <- tmp[,names(tmp) %in% selVars, with = FALSE]
    scoresNotInTable <- setdiff(scores, colnames(tmp))
    if(length(scoresNotInTable) > 0){
        warning(paste0(scoresNotInTable, collapse = " "), "metric(s) are not ",
                "present in STORM$RES$metric, and were not added to CALLS. ",
                "This is just a warning.")
    }
    tmp <- split(tmp, tmp$set)
    for(iN in names(tmp)){
        STORM$CALLS[[iN]][[RNAmod]] <- tmp[[iN]]
    }
    STORM
}


#' Assigning default Scores to tables by sets of samples
#'
#' @param STORM list. STORM object as output by \code{\link{storm_STORM}}
#'
#' @return
#' @export
#'
#' @examples
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


# BASE FUNCTIONS ###############################################################

# Between two groups ###########################################################

#' Add Start Ratio Difference
#'
#' Calculates the start ratio difference of group_A minus that
#' of group_B per set of samples in a STORM object.
#'
#' @param STORM list. STORM object as output by \code{\link{storm_STORM}}
#' @param group_A character. Name of group to be compared as found in STORM$META$group
#' @param group_B character. Name of group to be compared as found in STORM$META$group
#' @param newColName character. Name of calculated metric to be stored in STORM$RES,
#' assigned by default based on group_A and group_B
#' @param onNucs character. Nucleotide(s) in which the metric will be calculated
#' @param minCov integer. Minimal coverage in position for metric to be calculated
#'
#' @return
#' @export
#'
add_SRD <- function(STORM, group_A, group_B, newColName = "auto",
                    onNucs = c("A", "C", "G", "T"), minCov = 50){
    if(newColName == "auto"){
        newColName <- paste("SRD", group_A, group_B, sep = "_")
    }
    sets <- intersect(STORM$META[STORM$META$group == group_A,]$set, STORM$META[STORM$META$group == group_B,]$set)
    if(length(sets) == 0){
        warning("No sets found with both ", group_A, " and ", group_B, " group labels.\n",
                newColName, " was not calculated.")
        return(STORM)
    }
    OUT <- lapply(sets, function(iSet){
        id_A <- STORM$META[STORM$META$set == iSet & STORM$META$group == group_A,]$id
        id_B <- STORM$META[STORM$META$set == iSet & STORM$META$group == group_B,]$id
        DT_A <- add_StartRate(STORM$DATA[[id_A]], minCov = minCov)
        DT_B <- add_StartRate(STORM$DATA[[id_B]], minCov = minCov)
        tmpRES <- data.table::data.table(DT_A$startRatio - DT_B$startRatio)
        tmpRES <- data.table::data.table(iSet, newColName, tmpRES)
        names(tmpRES) <- c("set", "metric", "score")
        tmpRES <- data.table::data.table(DT_A[,c("chr", "gencoor", "strand", "gene",
                                                 "txcoor", "pos", "refSeq")], tmpRES)
        tmpRES[!(tmpRES$refSeq %in% onNucs), "score"] <- NA
        return(tmpRES)
    }) %>% do.call(what = rbind)
    STORM$RES <- hlpr_add_REScols(STORM$RES, OUT)
    return(STORM)
}

#' Add Start Ratio Difference 1bp down-stream
#'
#' Calculates the start ratio difference 1bp down-stream of group_A minus that
#' of group_B per set of samples in a STORM object.
#'
#' @param STORM list. STORM object as output by \code{\link{storm_STORM}}
#' @param group_A character. Name of group to be compared as found in STORM$META$group
#' @param group_B character. Name of group to be compared as found in STORM$META$group
#' @param newColName character. Name of calculated metric to be stored in STORM$RES,
#' assigned by default based on group_A and group_B
#' @param onNucs character. Nucleotide(s) in which the metric will be calculated
#' @param minCov integer. Minimal coverage in position for metric to be calculated
#'
#' @return
#' @export
#'
add_SRD1bpDS <- function(STORM, group_A, group_B, newColName = "auto",
                         onNucs = c("A", "C", "G", "T"), minCov = 50){
    if(newColName == "auto"){
        newColName <- paste("SRD1bpDS", group_A, group_B, sep = "_")
    }
    sets <- intersect(STORM$META[STORM$META$group == group_A,]$set, STORM$META[STORM$META$group == group_B,]$set)
    if(length(sets) == 0){
        warning("No sets found with both ", group_A, " and ", group_B, " group labels.\n",
                newColName, " was not calculated.")
        return(STORM)
    }
    OUT <- lapply(sets, function(iSet){
        id_A <- STORM$META[STORM$META$set == iSet & STORM$META$group == group_A,]$id
        id_B <- STORM$META[STORM$META$set == iSet & STORM$META$group == group_B,]$id
        DT_A <- add_StartRate_1bpDS(STORM$DATA[[id_A]], minCov = minCov)
        DT_B <- add_StartRate_1bpDS(STORM$DATA[[id_B]], minCov = minCov)
        tmpRES <- data.table::data.table(DT_A$startRate_1bpDS - DT_B$startRate_1bpDS)
        tmpRES <- data.table::data.table(iSet, newColName, tmpRES)
        names(tmpRES) <- c("set", "metric", "score")
        tmpRES <- data.table::data.table(DT_A[,c("chr", "gencoor", "strand", "gene",
                                                 "txcoor", "pos", "refSeq")], tmpRES)
        tmpRES[!(tmpRES$refSeq %in% onNucs), "score"] <- NA
        return(tmpRES)
    }) %>% do.call(what = rbind)
    STORM$RES <- hlpr_add_REScols(STORM$RES, OUT)
    return(STORM)
}

#' Add Start ratio log2-Fold-Change 1bp down-stream
#'
#' Calculates the log2-Fold-Change of the start ratio difference 1bp down-stream
#' of group_A over that of group_B per set of samples in a STORM object.
#'
#' @param STORM list. STORM object as output by \code{\link{storm_STORM}}
#' @param group_A character. Name of group to be compared as found in STORM$META$group
#' @param group_B character. Name of group to be compared as found in STORM$META$group
#' @param newColName character. Name of calculated metric to be stored in STORM$RES,
#' assigned by default based on group_A and group_B
#' @param onNucs character. Nucleotide(s) in which the metric will be calculated
#' @param minCov integer. Minimal coverage in position for metric to be calculated
#'
#' @return
#' @export
#'
#' @examples
add_SRlog2FCh1bpDS <- function(STORM, group_A, group_B, newColName = "auto",
                               onNucs = c("A", "C", "G", "T"), minCov = 50){
    if(newColName == "auto"){
        newColName <- paste("SRlog2FCh1bpDS", group_A, group_B, sep = "_")
    }
    sets <- intersect(STORM$META[STORM$META$group == group_A,]$set, STORM$META[STORM$META$group == group_B,]$set)
    if(length(sets) == 0){
        warning("No sets found with both ", group_A, " and ", group_B, " group labels.\n",
                newColName, " was not calculated.")
        return(STORM)
    }
    OUT <- lapply(sets, function(iSet){
        id_A <- STORM$META[STORM$META$set == iSet & STORM$META$group == group_A,]$id
        id_B <- STORM$META[STORM$META$set == iSet & STORM$META$group == group_B,]$id
        DT_A <- add_StartRate_1bpDS(STORM$DATA[[id_A]], minCov = minCov)
        DT_B <- add_StartRate_1bpDS(STORM$DATA[[id_B]], minCov = minCov)
        tmpRES <- data.table::data.table(log2(DT_A$startRate_1bpDS / DT_B$startRate_1bpDS))
        tmpRES <- data.table::data.table(iSet, newColName, tmpRES)
        names(tmpRES) <- c("set", "metric", "score")
        tmpRES <- data.table::data.table(DT_A[,c("chr", "gencoor", "strand", "gene",
                                                 "txcoor", "pos", "refSeq")], tmpRES)
        tmpRES[!(tmpRES$refSeq %in% onNucs), "score"] <- NA
        return(tmpRES)
    }) %>% do.call(what = rbind)
    STORM$RES <- hlpr_add_REScols(STORM$RES, OUT)
    return(STORM)
}

#' Add 2-O-methylation score to STORM-seq object
#'
#' Methylation score that considers the log2-FCH of start ratio between
#' two RNA-seq libraries (low-dNTPs and high-dNTPs), subtracting the "noise"
#' from a surrounding window.
#'
#' Incarnato et al. 2017: For Nm Incarnato and collabs. used the fold change
#' as the log2 of the ratio between the stoppage ratio at the low dNTPs
#' concentration and the stoppage ratio in the high dNTPs concentration sample.
#'  Because the stoppage ratio seems to be strongly region specific, they
#'  subtracted the local background defined as the mean ratio of the stoppage
#'  in the neighborhood of a given position (+/- 5 nucleotides, excluding
#'  given position).
#'
#' @param STORM list. STORM object as output by \code{\link{storm_STORM}}
#' @param lib_LowdNTPs character. Name of group with low-dNTPs library preparation
#' found in STORM$META$group
#' @param lib_HighdNTPs character. Name of group with high-dNTPs library preparation
#' found in STORM$META$group
#' @param newColName character. Name of calculated metric to be stored in STORM$RES,
#' assigned by default based on group_A and group_B
#' @param flankSize integer. Length of window to each side of position to
#' consider for metric calculation
#' @param perNuc_Zscaling
#' @param onNucs character. Nucleotide(s) in which the metric will be calculated
#' @param minCov integer. Minimal coverage in position for metric to be calculated
#'
#' @return
#' @export
#'
#' @examples
add_2OmeScore <- function(STORM,
                          lib_LowdNTPs,
                          lib_HighdNTPs,
                          newColName = "auto",
                          flankSize = 5,
                          perNuc_Zscaling = F,
                          onNucs = c("A", "C", "G", "T"), minCov = 50){
    if(newColName == "auto"){
        newColName <- paste("NmStopScore", lib_LowdNTPs, lib_HighdNTPs, sep = "_")
    }
    sets <- intersect(STORM$META[STORM$META$group == lib_LowdNTPs,]$set,
                      STORM$META[STORM$META$group == lib_HighdNTPs,]$set)
    if(length(sets) == 0){warning("No sets found with both ", lib_LowdNTPs,
                                  " and ", lib_HighdNTPs, " group labels.\n",
                                  newColName, " was not calculated.")
        return(STORM)}
    OUT <- lapply(sets, function(iSet){
        id_A <- STORM$META[STORM$META$set == iSet & STORM$META$group == lib_LowdNTPs,]$id
        id_B <- STORM$META[STORM$META$set == iSet & STORM$META$group == lib_HighdNTPs,]$id
        DT <- detect_2OmeScore(DT_low = STORM$DATA[[id_A]],
                               DT_high = STORM$DATA[[id_B]],
                               neighFlankSize = flankSize,
                               minCov = minCov)
        tmpRES <- data.table::data.table(twoOme_score = DT$twoOme_score)
        tmpRES <- data.table::data.table(iSet, newColName, tmpRES)
        names(tmpRES) <- c("set", "metric", "score")
        tmpRES <- data.table::data.table(DT[,c("chr", "gencoor", "strand", "gene",
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

#' Add Misincorporation Rate Difference
#'
#' @param STORM list. STORM object as output by \code{\link{storm_STORM}}
#' @param group_A character. Name of group to be compared as found in STORM$META$group
#' @param group_B character. Name of group to be compared as found in STORM$META$group
#' @param newColName character. Name of calculated metric to be stored in STORM$RES,
#' assigned by default based on group_A and group_B
#' @param onNucs character. Nucleotide(s) in which the metric will be calculated
#' @param minCov integer. Minimal coverage in position for metric to be calculated
#'
#' @return
#' @export
#'
#' @examples
add_MRD <- function(STORM, group_A, group_B, newColName = "auto",
                    onNucs = c("A", "C", "G", "T"), minCov = 50){
    if(newColName == "auto"){
        newColName <- paste("MRD", group_A, group_B, sep = "_")
    }
    sets <- intersect(STORM$META[STORM$META$group == group_A,]$set, STORM$META[STORM$META$group == group_B,]$set)
    if(length(sets) == 0){
        warning("No sets found with both ", group_A, " and ", group_B, " group labels.\n",
                newColName, " was not calculated.")
        return(STORM)
    }
    OUT <- lapply(sets, function(iSet){
        id_A <- STORM$META[STORM$META$set == iSet & STORM$META$group == group_A,]$id
        id_B <- STORM$META[STORM$META$set == iSet & STORM$META$group == group_B,]$id
        DT_A <- STORM$DATA[[id_A]] %>% add_diffNucToRefRatio(minCov = minCov)
        DT_B <- STORM$DATA[[id_B]] %>% add_diffNucToRefRatio(minCov = minCov)
        tmpRES <- data.table::data.table(DT_A$diffToRef_Ratio - DT_B$diffToRef_Ratio)
        tmpRES <- data.table::data.table(iSet, newColName, tmpRES)
        names(tmpRES) <- c("set", "metric", "score")
        tmpRES <- data.table::data.table(DT_A[,c("chr", "gencoor", "strand", "gene",
                                                 "txcoor", "pos", "refSeq")], tmpRES)
        tmpRES[!(tmpRES$refSeq %in% onNucs), "score"] <- NA
        return(tmpRES)
    }) %>% do.call(what = rbind)
    STORM$RES <- hlpr_add_REScols(STORM$RES, OUT)
    return(STORM)
}

#' Add cytidine persistance index to STORM object
#'
#' Cytidine score modified from (Squires et al. 2012). It accounts for positions
#' with Thymine reads in Cytidine positions in the Control/Mock sample.
#' CytPer = (C_bisulphite / (C_bisulphite + T_bisulphite)) - (1 - (C_control / C_control + T_control))
#'
#' @param STORM list. STORM object as output by \code{\link{storm_STORM}}
#' @param group_A character. Name of group to be compared as found in STORM$META$group
#' @param group_B character. Name of group to be compared as found in STORM$META$group
#' @param newColName character. Name of calculated metric to be stored in STORM$RES,
#' assigned by default based on group_A and group_B
#' @param onNucs character. Nucleotide(s) in which the metric will be calculated
#'
#' @return
#' @export
#'
#' @examples
add_CytPer <- function(STORM, group_A, group_B, newColName = "auto",
                       onNucs = c("C"), minCov = 50){
    if(newColName == "auto"){
        newColName <- paste("CytPer", group_A, group_B, sep = "_")
    }
    sets <- intersect(STORM$META[STORM$META$group == group_A,]$set, STORM$META[STORM$META$group == group_B,]$set)
    if(length(sets) == 0){
        warning("No sets found with both ", group_A, " and ", group_B, " group labels.\n",
                newColName, " was not calculated.")
        return(STORM)
    }
    OUT <- lapply(sets, function(iSet){
        id_A <- STORM$META[STORM$META$set == iSet & STORM$META$group == group_A,]$id
        id_B <- STORM$META[STORM$META$set == iSet & STORM$META$group == group_B,]$id
        DT_A <- STORM$DATA[[id_A]] %>% detect_m5C(minCov = minCov)
        DT_B <- STORM$DATA[[id_B]] %>% detect_m5C(minCov = minCov)
        tmpRES <- data.table::data.table(DT_A$det_m5C - (1 - DT_B$det_m5C))
        tmpRES <- data.table::data.table(iSet, newColName, tmpRES)
        names(tmpRES) <- c("set", "metric", "score")
        tmpRES <- data.table::data.table(DT_A[,c("chr", "gencoor", "strand", "gene",
                                                 "txcoor", "pos", "refSeq")], tmpRES)
        tmpRES[!(tmpRES$refSeq %in% onNucs), "score"] <- NA
        return(tmpRES)
    }) %>% do.call(what = rbind)
    STORM$RES <- hlpr_add_REScols(STORM$RES, OUT)
    return(STORM)
}

#' Add C-to-T misincorporation difference
#'
#' Calculates the difference of T conversion at C positions by subtracting the
#' result of "T"/("T"+"C") between groups. Commonly NaCNBH3-Treated vs
#' Control or vs Deacetylated negative control.
#'
#' The reaction of ac4C with sodium cyanoborohydride (NaCNBH3) under acidic
#' conditions forms the reduced nucleobase N4-acetyltetrahydrocytidine.
#' The altered structure of this reduced nucleobase compared with ac4C causes
#' the incorporation of non-cognate deoxynucleotide triphosphates (dNTPs)
#' upon reverse transcription (Sas-Chen et al. 2020).
#'
#' @param STORM list. STORM object as output by \code{\link{storm_STORM}}
#' @param group_A character. Name of group to be compared as found in STORM$META$group
#' @param group_B character. Name of group to be compared as found in STORM$META$group
#' @param newColName character. Name of calculated metric to be stored in STORM$RES,
#' assigned by default based on group_A and group_B
#' @param onNucs character. Nucleotide(s) in which the metric will be calculated
#'
#' @return
#' @export
#'
#' @examples
add_CtoT_MRD <- function(STORM, group_A, group_B, newColName = "auto",
                         onNucs = c("C"), minCov = 50){
    if(newColName == "auto"){
        newColName <- paste("CtoT.MRD", group_A, group_B, sep = "_")
    }
    sets <- intersect(STORM$META[STORM$META$group == group_A,]$set, STORM$META[STORM$META$group == group_B,]$set)
    if(length(sets) == 0){
        warning("No sets found with both ", group_A, " and ", group_B, " group labels.\n",
                newColName, " was not calculated.")
        return(STORM)
    }
    OUT <- lapply(sets, function(iSet){
        id_A <- STORM$META[STORM$META$set == iSet & STORM$META$group == group_A,]$id
        id_B <- STORM$META[STORM$META$set == iSet & STORM$META$group == group_B,]$id
        DT_A <- STORM$DATA[[id_A]] %>% detect_ac4C(minCov = minCov)
        DT_B <- STORM$DATA[[id_B]] %>% detect_ac4C(minCov = minCov)
        tmpRES <- data.table::data.table(DT_A$det_ac4C - DT_B$det_ac4C)
        tmpRES <- data.table::data.table(iSet, newColName, tmpRES)
        names(tmpRES) <- c("set", "metric", "score")
        tmpRES <- data.table::data.table(DT_A[,c("chr", "gencoor", "strand", "gene",
                                                 "txcoor", "pos", "refSeq")], tmpRES)
        tmpRES[!(tmpRES$refSeq %in% onNucs), "score"] <- NA
        return(tmpRES)
    }) %>% do.call(what = rbind)
    STORM$RES <- hlpr_add_REScols(STORM$RES, OUT)
    return(STORM)
}

#' Add deletion rate difference
#'
#' Calculated the deletion rate to nucleotide reads for each position and
#' subtracts it for the rate in the control group.
#'
#' @param STORM list. STORM object as output by \code{\link{storm_STORM}}
#' @param group_A character. Name of group to be compared as found in STORM$META$group
#' @param group_B character. Name of group to be compared as found in STORM$META$group
#' @param newColName character. Name of calculated metric to be stored in STORM$RES,
#' assigned by default based on group_A and group_B
#' @param onNucs character. Nucleotide(s) in which the metric will be calculated
#' @param minReadCounts integer. Minimal number of nucleotide reads (NOT coverage)
#' in position for metric to be calculated
#'
#' @return
#' @export
#'
#' @examples
add_DRD <- function(STORM, group_A, group_B, newColName = "auto",
                    onNucs = c("A", "C", "G", "T"), minReadCounts = 20){
    if(newColName == "auto"){
        newColName <- paste("DRD", group_A, group_B, sep = "_")
    }
    sets <- intersect(STORM$META[STORM$META$group == group_A,]$set, STORM$META[STORM$META$group == group_B,]$set)
    if(length(sets) == 0){
        warning("No sets found with both ", group_A, " and ", group_B, " group labels.\n",
                newColName, " was not calculated.")
        return(STORM)
    }
    OUT <- lapply(sets, function(iSet){
        id_A <- STORM$META[STORM$META$set == iSet & STORM$META$group == group_A,]$id
        id_B <- STORM$META[STORM$META$set == iSet & STORM$META$group == group_B,]$id
        DT_A <- STORM$DATA[[id_A]] %>% deletion_rate(minReadCounts = minReadCounts)
        DT_B <- STORM$DATA[[id_B]] %>% deletion_rate(minReadCounts = minReadCounts)
        tmpRES <- data.table::data.table(DT_A$deletion_Ratio - DT_B$deletion_Ratio)
        tmpRES <- data.table::data.table(iSet, newColName, tmpRES)
        names(tmpRES) <- c("set", "metric", "score")
        tmpRES <- data.table::data.table(DT_A[,c("chr", "gencoor", "strand", "gene",
                                                 "txcoor", "pos", "refSeq")], tmpRES)
        tmpRES[!(tmpRES$refSeq %in% onNucs), "score"] <- NA
        return(tmpRES)
    }) %>% do.call(what = rbind)
    STORM$RES <- hlpr_add_REScols(STORM$RES, OUT)
    return(STORM)
}

# One group metrics ############################################################

#' Add Misincorporation Rate
#'
#' @param STORM list. STORM object as output by \code{\link{storm_STORM}}
#' @param group_A character. Name of group to be compared as found in STORM$META$group
#' @param newColName character. Name of calculated metric to be stored in STORM$RES,
#' assigned by default based on group_A and group_B
#' @param onNucs character. Nucleotide(s) in which the metric will be calculated
#'
#' @return
#' @export
#'
#' @examples
add_MR <- function(STORM, group_A, newColName = "auto",
                   onNucs = c("A", "C", "G", "T"), minCov = 50){
    if(newColName == "auto"){
        newColName <- paste("MR", group_A, sep = "_")
    }
    sets <- STORM$META[STORM$META$group == group_A,]$set
    if(length(sets) == 0){
        warning("No sets found with ", group_A, " group label.\n",
                newColName, " was not calculated.")
        return(STORM)
    }
    OUT <- lapply(sets, function(iSet){
        id_A <- STORM$META[STORM$META$set == iSet & STORM$META$group == group_A,]$id
        DT_A <- add_diffNucToRefRatio(STORM$DATA[[id_A]], minCov = minCov)
        tmpRES <- data.table::data.table(DT_A$diffToRef_Ratio)
        tmpRES <- data.table::data.table(iSet, newColName, tmpRES)
        names(tmpRES) <- c("set", "metric", "score")
        tmpRES <- data.table::data.table(DT_A[,c("chr", "gencoor", "strand", "gene",
                                                 "txcoor", "pos", "refSeq")], tmpRES)
        tmpRES[!(tmpRES$refSeq %in% onNucs), "score"] <- NA
        return(tmpRES)
    }) %>% do.call(what = rbind)
    STORM$RES <- hlpr_add_REScols(STORM$RES, OUT)
    return(STORM)
}

#' Add Read-start to coverage ratio
#'
#' @param STORM list. STORM object as output by \code{\link{storm_STORM}}
#' @param group_A character. Name of group to be compared as found in STORM$META$group
#' @param newColName character. Name of calculated metric to be stored in STORM$RES,
#' assigned by default based on group_A and group_B
#' @param onNucs character. Nucleotide(s) in which the metric will be calculated
#' @param minCov integer. Minimal coverage in position for metric to be calculated
#'
#' @return
#' @export
#'
#' @examples
add_SR <- function(STORM, group_A, newColName = "auto",
                   onNucs = c("A", "C", "G", "T"), minCov = 50){
    if(newColName == "auto"){
        newColName <- paste("SR", group_A, sep = "_")
    }
    sets <- STORM$META[STORM$META$group == group_A,]$set
    if(length(sets) == 0){
        warning("No sets found with ", group_A, " group label.\n",
                newColName, " was not calculated.")
        return(STORM)
    }
    OUT <- lapply(sets, function(iSet){
        id_A <- STORM$META[STORM$META$set == iSet & STORM$META$group == group_A,]$id
        DT_A <- add_StartRate(STORM$DATA[[id_A]])
        tmpRES <- data.table::data.table(DT_A$start_Ratio)
        tmpRES <- data.table::data.table(iSet, newColName, tmpRES)
        names(tmpRES) <- c("set", "metric", "score")
        tmpRES <- data.table::data.table(DT_A[,c("chr", "gencoor", "strand", "gene",
                                                 "txcoor", "pos", "refSeq")], tmpRES)
        tmpRES[!(tmpRES$refSeq %in% onNucs), "score"] <- NA
        return(tmpRES)
    }) %>% do.call(what = rbind)
    STORM$RES <- hlpr_add_REScols(STORM$RES, OUT)
    return(STORM)
}

#' Add read-start to coverage rate 1 bp down-stream
#'
#' @param STORM list. STORM object as output by \code{\link{storm_STORM}}
#' @param group_A character. Name of group to be compared as found in STORM$META$group
#' @param newColName character. Name of calculated metric to be stored in STORM$RES,
#' assigned by default based on group_A and group_B
#' @param onNucs character. Nucleotide(s) in which the metric will be calculated
#' @param minCov integer. Minimal coverage in position for metric to be calculated
#'
#' @return
#' @export
#'
#' @examples
add_SR1bpDS <- function(STORM, group_A, newColName = "auto",
                         onNucs = c("A", "C", "G", "T"), minCov = 50){
    if(newColName == "auto"){
        newColName <- paste("SR1bpDS", group_A, sep = "_")
    }
    sets <- STORM$META[STORM$META$group == group_A,]$set
    if(length(sets) == 0){
        warning("No sets found with ", group_A, " group label.\n",
                newColName, " was not calculated.")
        return(STORM)
    }
    OUT <- lapply(sets, function(iSet){
        id_A <- STORM$META[STORM$META$set == iSet & STORM$META$group == group_A,]$id
        DT_A <- add_StartRate_1bpDS(STORM$DATA[[id_A]], minCov = minCov)
        tmpRES <- data.table::data.table(DT_A$startRate_1bpDS)
        tmpRES <- data.table::data.table(iSet, newColName, tmpRES)
        names(tmpRES) <- c("set", "metric", "score")
        tmpRES <- data.table::data.table(DT_A[,c("chr", "gencoor", "strand", "gene",
                                                 "txcoor", "pos", "refSeq")], tmpRES)
        tmpRES[!(tmpRES$refSeq %in% onNucs), "score"] <- NA
        return(tmpRES)
    }) %>% do.call(what = rbind)
    STORM$RES <- hlpr_add_REScols(STORM$RES, OUT)
    return(STORM)
}

#' Add read-end to coverage ratio
#'
#' @param STORM list. STORM object as output by \code{\link{storm_STORM}}
#' @param group_A character. Name of group to be compared as found in STORM$META$group
#' @param newColName character. Name of calculated metric to be stored in STORM$RES,
#' assigned by default based on group_A and group_B
#' @param onNucs character. Nucleotide(s) in which the metric will be calculated
#' @param minCov integer. Minimal coverage in position for metric to be calculated
#'
#' @return
#' @export
#'
#' @examples
add_ER <- function(STORM, group_A, newColName = "auto",
                   onNucs = c("A", "C", "G", "T"), minCov = 50){
    if(newColName == "auto"){
        newColName <- paste("ER", group_A, sep = "_")
    }
    sets <- STORM$META[STORM$META$group == group_A,]$set
    if(length(sets) == 0){
        warning("No sets found with ", group_A, " group label.\n",
                newColName, " was not calculated.")
        return(STORM)
    }
    OUT <- lapply(sets, function(iSet){
        id_A <- STORM$META[STORM$META$set == iSet & STORM$META$group == group_A,]$id
        DT_A <- add_EndRate(STORM$DATA[[id_A]], minCov = minCov)
        tmpRES <- data.table::data.table(DT_A$end_Ratio)
        tmpRES <- data.table::data.table(iSet, newColName, tmpRES)
        names(tmpRES) <- c("set", "metric", "score")
        tmpRES <- data.table::data.table(DT_A[,c("chr", "gencoor", "strand", "gene",
                                                 "txcoor", "pos", "refSeq")], tmpRES)
        tmpRES[!(tmpRES$refSeq %in% onNucs), "score"] <- NA
        return(tmpRES)
    }) %>% do.call(what = rbind)
    STORM$RES <- hlpr_add_REScols(STORM$RES, OUT)
    return(STORM)
}

#' Add ScoreA accounting for read-ends (3prime)
#'
#' @param STORM list. STORM object as output by \code{\link{storm_STORM}}
#' @param group_A character. Name of group to be compared as found in STORM$META$group
#' @param newColName character. Name of calculated metric to be stored in STORM$RES,
#' assigned by default based on group_A and group_B
#' @param onNucs character. Nucleotide(s) in which the metric will be calculated
#' @param minMedCov integer. Minimal median coverage of window around position,
#' delimited by flankSize
#' @param flankSize integer. Length of window to each side of position to
#' consider for metric calculation
#'
#' @return
#' @export
#' @aliases storm_add_scoreA3p
#'
#' @examples
add_scoreA3p <- function(STORM, group_A, newColName = "auto",
                         onNucs = c("A", "C", "G", "T"), minMedCov = 30,
                         flankSize = 6){
    if(newColName == "auto"){
        newColName <- paste("ScoreA3p", group_A, sep = "_")
    }
    sets <- STORM$META[STORM$META$group == group_A,]$set
    if(length(sets) == 0){
        warning("No sets found with ", group_A, " group label.\n",
                newColName, " was not calculated.")
        return(STORM)
    }
    OUT <- lapply(sets, function(iSet){
        id_A <- STORM$META[STORM$META$set == iSet & STORM$META$group == group_A,]$id
        DT_A <- add_scoreA_3p(STORM$DATA[[id_A]], minMedEnv = minMedCov, flankSize = flankSize)
        tmpRES <- data.table::data.table(DT_A$scoreA_3p)
        tmpRES <- data.table::data.table(iSet, newColName, tmpRES)
        names(tmpRES) <- c("set", "metric", "score")
        tmpRES <- data.table::data.table(DT_A[,c("chr", "gencoor", "strand", "gene",
                                                 "txcoor", "pos", "refSeq")], tmpRES)
        tmpRES[!(tmpRES$refSeq %in% onNucs), "score"] <- NA
        return(tmpRES)
    }) %>% do.call(what = rbind)
    STORM$RES <- hlpr_add_REScols(STORM$RES, OUT)
    return(STORM)
}

#' Add ScoreA accounting for read-starts (5prime)
#'
#' @param STORM list. STORM object as output by \code{\link{storm_STORM}}
#' @param group_A character. Name of group to be compared as found in STORM$META$group
#' @param newColName character. Name of calculated metric to be stored in STORM$RES,
#' assigned by default based on group_A and group_B
#' @param onNucs character. Nucleotide(s) in which the metric will be calculated
#' @param minMedCov integer. Minimal median coverage of window around position,
#' delimited by flankSize
#' @param flankSize integer. Length of window to each side of position to
#' consider for metric calculation
#'
#' @return
#' @export
#'
#' @examples
#' @aliases storm_add_scoreA5p
add_scoreA5p <- function(STORM, group_A, newColName = "auto",
                         onNucs = c("A", "C", "G", "T"), minMedCov = 30,
                         flankSize = 6){
    if(newColName == "auto"){
        newColName <- paste("ScoreA5p", group_A, sep = "_")
    }
    sets <- STORM$META[STORM$META$group == group_A,]$set
    if(length(sets) == 0){
        warning("No sets found with ", group_A, " group label.\n",
                newColName, " was not calculated.")
        return(STORM)
    }
    OUT <- lapply(sets, function(iSet){
        id_A <- STORM$META[STORM$META$set == iSet & STORM$META$group == group_A,]$id
        DT_A <- add_scoreA_5p(STORM$DATA[[id_A]], minMedEnv = minMedCov,
                              flankSize = flankSize)
        tmpRES <- data.table::data.table(DT_A$scoreA_5p)
        tmpRES <- data.table::data.table(iSet, newColName, tmpRES)
        names(tmpRES) <- c("set", "metric", "score")
        tmpRES <- data.table::data.table(DT_A[, c("chr", "gencoor", "strand", "gene",
                                                  "txcoor", "pos", "refSeq")], tmpRES)
        tmpRES[!(tmpRES$refSeq %in% onNucs), "score"] <- NA
        return(tmpRES)
    }) %>% do.call(what = rbind)
    STORM$RES <- hlpr_add_REScols(STORM$RES, OUT)
    return(STORM)
}

#' Add deletion rate
#'
#' Adds the deletion to nucleotide reads ratio
#'
#' @param STORM list. STORM object as output by \code{\link{storm_STORM}}
#' @param group_A character. Name of group to be compared as found in STORM$META$group
#' @param newColName character. Name of calculated metric to be stored in STORM$RES,
#' assigned by default based on group_A and group_B
#' @param onNucs character. Nucleotide(s) in which the metric will be calculated
#' @param minReadCounts integer. Minimal number of nucleotide reads (NOT coverage)
#' in position for metric to be calculated
#'
#' @return
#' @export
#'
#' @examples
add_DR <- function(STORM, group_A, newColName = "auto",
                   onNucs = c("A", "C", "G", "T"), minReadCounts = 30){
    if(newColName == "auto"){
        newColName <- paste("DR", group_A, sep = "_")
    }
    sets <- STORM$META[STORM$META$group == group_A,]$set
    if(length(sets) == 0){
        warning("No sets found with ", group_A, " group label.\n",
                newColName, " was not calculated.")
        return(STORM)
    }
    OUT <- lapply(sets, function(iSet){
        id_A <- STORM$META[STORM$META$set == iSet & STORM$META$group == group_A,]$id
        DT_A <- deletion_rate(STORM$DATA[[id_A]], minReadCounts = minReadCounts)
        tmpRES <- data.table::data.table(DT_A$deletion_Ratio)
        tmpRES <- data.table::data.table(iSet, newColName, tmpRES)
        names(tmpRES) <- c("set", "metric", "score")
        tmpRES <- data.table::data.table(DT_A[,c("chr", "gencoor", "strand", "gene",
                                                 "txcoor", "pos", "refSeq")], tmpRES)
        tmpRES[!(tmpRES$refSeq %in% onNucs), "score"] <- NA
        return(tmpRES)
    }) %>% do.call(what = rbind)
    STORM$RES <- hlpr_add_REScols(STORM$RES, OUT)
    return(STORM)
}

# Helper functions #############################################################
# Function for cytidine persistence to Bisulphite treatment
detect_m5C <- function(DT, minCov = 50){
    tmp <- round(DT$C / (DT$T + DT$C), 4)
    tmp[DT$cov < minCov] <- NA
    tibble::add_column(DT, det_m5C = tmp)
}

# Function for cytidine persistence to Bisulphite treatment
detect_ac4C <- function(DT, minCov = 50){
    DT <- data.table::data.table(DT)
    tmp <- round(DT$`T` / (DT$`T` + DT$`C`), 6)
    tmp[DT$cov < minCov] <- NA
    tibble::add_column(DT, det_ac4C = tmp)
}

# Add column of Different Nucleotides to reference ratio, diffToRef and nucTotal columns are required
add_diffNucToRefRatio <- function(DT, minCov = 50){
    DT <- txtools::tx_add_diffNucToRef(DT) %>% txtools::tx_add_nucTotal()
    tmp <- round(DT$diffToRef / DT$nucTotal, 6)
    tmp[DT$cov < minCov] <- NA
    tibble::add_column(DT, diffToRef_Ratio = tmp)
}

# Start rate
add_StartRate <- function(DT, minCov = 50){
    tmp <- (DT$start_5p + 1) / (DT$cov + 1)
    tmp[DT$cov < minCov] <- NA
    tibble::add_column(DT, start_Ratio = tmp)
}

# Start rate
add_EndRate <- function(DT, minCov = 50){
    tmp <- (DT$end_3p + 1) / (DT$cov + 1)
    tmp[DT$cov < minCov] <- NA
    tibble::add_column(DT, end_Ratio = tmp)
}

# Start rate
add_StartRate_1bpDS <- function(DT, minCov = 50){
    tmp <- (DT$start_5p + 1) / (DT$cov + 1)
    tmp[DT$cov < minCov] <- NA
    DTL <- tibble::add_column(DT, startRate_1bpDS = tmp) %>% txtools::tx_split_DT()
    DT <- lapply(DTL, function(DT){
        DT$startRate_1bpDS <- c(utils::tail(DT$startRate_1bpDS, -1), NA)
        DT
    }) %>% txtools::tx_merge_DT()
    return(DT)
}

# 2Omescore calculation: Score A
detect_2OmeScore <- function(DT_low, DT_high, neighFlankSize = 5, minCov = 50){
    DT_low <- DT_low %>% add_StartRate(minCov = minCov)
    DT_high <- DT_high %>% add_StartRate(minCov = minCov)
    if(all(union(DT_low$gene, DT_high$gene) %in% intersect(DT_low$gene, DT_high$gene))){
        iGenes <- intersect(DT_low$gene, DT_high$gene)
    }else{stop("DT_low and DT_high do not share the same genes")}
    DTL1 <- DT_low %>% txtools::tx_split_DT()
    DTL2 <- DT_high %>% txtools::tx_split_DT()
    OUT <- lapply(iGenes, function(iG){
        RES <- data.table::data.table(DTL1[[iG]][,colnames(DTL1[[iG]]) %in% storm_baseCols, with = FALSE])
        tmpSR <- c(utils::tail(DTL1[[iG]]$start_Ratio / DTL2[[iG]]$start_Ratio, - 1), NA)
        sc_range <- (1 + neighFlankSize):(nrow(DTL1[[iG]])-neighFlankSize)
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

# Score A Birkedal et al., 2015
scoreA <- function(CHR, d = 6, minMedEnv = 30){
    scoreA <- rep(NA, length(CHR)); names(scoreA) <- names(CHR)
    for (P in (d+1):(length(CHR)-d)){
        ePos <- CHR[P]
        leftFlank <- CHR[(P-d):(P-1)]
        rightFlank <- CHR[(P+1):(P+d)]
        if(min(stats::median(c(rightFlank, leftFlank), na.rm = T)) < minMedEnv ){next()}
        m1 <- mean(leftFlank)
        m2 <- mean(rightFlank)
        sd1 <- stats::sd(leftFlank)
        sd2 <- stats::sd(rightFlank)
        scoreA[P] <- 1 - (((2 * ePos) + 1) /
                              ((0.5*abs(m1-sd1)) + ePos + (0.5*abs(m2-sd2)) + 1))
    }
    scoreA[scoreA < 0] <- 0
    return(scoreA)
}

# Adding score A to multiGene DT
add_scoreA_3p <- function(DT, minMedEnv = 30, nCores = 1, flankSize = 6){
    DTL <- txtools::tx_split_DT(DT)
    parallel::mclapply(mc.cores = nCores, DTL, function(x){
        tmp <- scoreA(x$end_3p, minMedEnv = minMedEnv, d = flankSize)
        tibble::add_column(x, scoreA_3p = tmp)
    }) %>% txtools::tx_merge_DT()
}

# Adding score A 5p_ends to multiGene DT
add_scoreA_5p <- function(DT, minMedEnv = 15, nCores = 1, flankSize = 6){
    DTL <- txtools::tx_split_DT(DT)
    parallel::mclapply(mc.cores = nCores, DTL, function(x){
        if(!("startRate_1bpDS" %in% names(x))){
            x <- add_StartRate_1bpDS(x)
        }
        tmp <- scoreA(x$start_5p, minMedEnv = minMedEnv, d = flankSize)
        tmp <- c(utils::tail(tmp, -1), NA)
        tibble::add_column(x, scoreA_5p = tmp)
    }) %>% txtools::tx_merge_DT()
}

# Start rate
add_StartRate <- function(DT, minCov = 50){
    tmp <- (DT$start_5p + 1) / (DT$cov + 1)
    tmp[DT$cov < minCov] <- NA
    tibble::add_column(DT, start_Ratio = tmp)
}

# Calculate deletion rate
deletion_rate <- function(DT, minReadCounts = 30){
    DT <- txtools::tx_add_nucTotal(DT)
    tmp <- round(DT$`-` / DT$nucTotal, 6)
    tmp[DT$nucTotal < minReadCounts] <- NA
    tibble::add_column(DT, deletion_Ratio = tmp)
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
