# PREDEFINED METRICS ###########################################################

#' Add all metrics for STORMseq v0.0.2
#'
#' @param STORM
#'
#' @return
#' @export
#'
#' @examples
add_allRNAmod_metrics <- function(STORM){
    STORM %>%
        add_Nm_metrics() %>%
        add_Y_metrics() %>%
        add_ac4C_metrics() %>%
        add_m1A_metrics() %>%
        add_m7G_metrics() %>%
        add_m5C_metrics() %>%
        add_m3U_metrics()
}

add_Y_metrics <- function(STORM){
    STORM %>%
        try(add_SRD1bpDS("CMC.TGIRT", "Mock.TGIRT"), silent = T) %>%
        try(add_SRD1bpDS("CMC.SSIII", "Mock.SSIII"), silent = T) %>%
        try(add_SRlog2FCh1bpDS("CMC.TGIRT", "Mock.TGIRT"), silent = T) %>%
        try(add_SRlog2FCh1bpDS("CMC.SSIII", "Mock.SSIII"), silent = T) %>%
        try(add_DRD("RBSseqHeatMg.TGIRT", "Mock.TGIRT"), silent = T) %>%
        try(add_DRD("RBSseqHeatMg.SSIII", "Mock.SSIII"), silent = T) %>%
        try(add_DRD("m5C.TGIRT", "Mock.TGIRT"), silent = T) %>%
        try(add_DRD("m5C.SSIII", "Mock.SSIII"), silent = T)
}

add_Nm_metrics <- function(STORM){
    STORM %>%
        storm_add_scoreA3p("Mock.TGIRT", minMedCov = 50) %>%
        storm_add_scoreA3p("Mock.SSIII", minMedCov = 50) %>%
        add_2OmeScore("MocklowdNTPs.SSIII", "Mock.SSIII", perNuc_Zscaling = T)
}

add_ac4C_metrics <- function(STORM){
    STORM %>%
        add_CtoT_MRD("Ac4C.TGIRT", "Mock.TGIRT") %>%
        add_CtoT_MRD("Ac4C.TGIRT", "DeacetylatedAc4C.TGIRT") %>%
        add_CtoT_MRD("Ac4C.SSIII", "Mock.SSIII") %>%
        add_CtoT_MRD("Ac4C.SSIII", "DeacetylatedAc4C.SSIII")
}

add_m1A_metrics <- function(STORM){
    STORM %>%
        add_SRD1bpDS("Mock.SSIII", "Dimroth.SSIII") %>%
        add_SRlog2FCh1bpDS("Mock.SSIII", "Dimroth.SSIII") %>%
        add_MRD("Mock.TGIRT", "Dimroth.TGIRT") %>%
        add_SRD1bpDS("Mock.TGIRT", "Dimroth.TGIRT") %>%
        add_SRlog2FCh1bpDS("Mock.TGIRT", "Dimroth.TGIRT")
}

add_m7G_metrics <- function(STORM){
    STORM %>%
        add_SRlog2FCh1bpDS("DeacetylatedAc4C.TGIRT", "Ac4C.TGIRT") %>%
        add_SRlog2FCh1bpDS("DeacetylatedAc4C.SSIII", "Ac4C.SSIII") %>%
        add_SRD1bpDS("DeacetylatedAc4C.TGIRT", "Ac4C.TGIRT") %>%
        add_SRD1bpDS("DeacetylatedAc4C.SSIII", "Ac4C.SSIII") %>%
        add_MRD("DeacetylatedAc4C.TGIRT", "Ac4C.TGIRT") %>%
        add_MRD("DeacetylatedAc4C.SSIII", "Ac4C.SSIII") %>%
        add_SRD("NaBH4HydBiotin.TGIRT", "Mock.TGIRT") %>%
        add_SRD("NaBH4HydBiotin.RTHIV", "Mock.RTHIV")
}

add_m5C_metrics <- function(STORM){
    STORM %>%
        add_CytPer("m5C.TGIRT", "Mock.TGIRT") %>%
        add_CytPer("RBSseqHeatMg.TGIRT", "Mock.TGIRT")
}

add_m3U_metrics <- function(STORM){
    STORM %>%
        add_MRD("Mock.TGIRT", "AlkBmix.TGIRT")
}

Y_scores <- c("SRD1bpDS_CMC.TGIRT_Mock.TGIRT",
              "SRD1bpDS_CMC.SSIII_Mock.SSIII",
              "SRlog2FCh1bpDS_CMC.TGIRT_Mock.TGIRT",
              "SRlog2FCh1bpDS_CMC.SSIII_Mock.SSIII")

Nm_scores <- c("NmStopScore_MocklowdNTPs.SSIII_Mock.SSIII",
               "ScoreA3p_Mock.TGIRT",
               "ScoreA3p_Mock.SSIII")

ac4C_scores <- c("CtoT.MRD_Ac4C.TGIRT_Mock.TGIRT",
                 "CtoT.MRD_Ac4C.TGIRT_DeacetylatedAc4C.TGIRT",
                 "CtoT.MRD_Ac4C.SSIII_Mock.SSIII",
                 "CtoT.MRD_Ac4C.SSIII_DeacetylatedAc4C.SSIII")

m1A_scores <- c("SRD1bpDS_Mock_Dimroth",
                "MRD_Mock_Dimroth",
                "MRD_Mock.TGIRT_Dimroth.TGIRT",
                "SRD1bpDS_Mock.SSIII_Dimroth.SSIII")

m7G_scores <- c("MRD_DeacetylatedAc4C.TGIRT_Ac4C.TGIRT",
                "MRD_DeacetylatedAc4C.SSIII_Ac4C.SSIII",
                "SRD1bpDS_DeacetylatedAc4C.SSIII_Ac4C.SSIII")

m5C_scores <- c("CytPer_m5C.TGIRT_Mock.TGIRT",
                "CytPer_RBSseqHeatMg.TGIRT_Mock.TGIRT")

m3U_scores <- c("MRD_Mock_AlkBmix")

# METRIC ASSIGNMENT FUNCTIONS ##################################################

#' Initialize STORM$CALLS
#'
#' @param STORM
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
#' @param STORM
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
    tmp <- split(tmp, tmp$set)
    for(iN in names(tmp)){
        STORM$CALLS[[iN]][[RNAmod]] <- tmp[[iN]]
    }
    STORM
}


#' Assigning predefined Scores
#'
#' @param STORM
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

# Add SRD to STORM object
add_SRD <- function(STORM, group_A, group_B, newColName = "auto",
                    onNucs = c("A", "C", "G", "T")){
    if(newColName == "auto"){
        newColName <- paste("SRD", group_A, group_B, sep = "_")
    }
    sets <- intersect(STORM$META[STORM$META$group == group_A,]$set, STORM$META[STORM$META$group == group_B,]$set)
    if(length(sets) == 0){stop("No sets found with both group_A and group_B label")}
    OUT <- lapply(sets, function(iSet){
        id_A <- STORM$META[STORM$META$set == iSet & STORM$META$group == group_A,]$id
        id_B <- STORM$META[STORM$META$set == iSet & STORM$META$group == group_B,]$id
        if(S4Vectors::isEmpty(id_A)){stop(group_A, " not found in STORM object META data")}
        if(S4Vectors::isEmpty(id_B)){stop(group_B, " not found in STORM object META data")}
        DT_A <- add_StartRate(STORM$DATA[[id_A]])
        DT_B <- add_StartRate(STORM$DATA[[id_B]])
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

# Add SRD 1bp down-stream to STORM object
add_SRD1bpDS <- function(STORM, group_A, group_B, newColName = "auto",
                         onNucs = c("A", "C", "G", "T")){
    if(newColName == "auto"){
        newColName <- paste("SRD1bpDS", group_A, group_B, sep = "_")
    }
    sets <- intersect(STORM$META[STORM$META$group == group_A,]$set, STORM$META[STORM$META$group == group_B,]$set)
    if(length(sets) == 0){stop("No sets found with both group_A and group_B label")}
    OUT <- lapply(sets, function(iSet){
        id_A <- STORM$META[STORM$META$set == iSet & STORM$META$group == group_A,]$id
        id_B <- STORM$META[STORM$META$set == iSet & STORM$META$group == group_B,]$id
        if(S4Vectors::isEmpty(id_A)){stop(group_A, " not found in STORM object META data")}
        if(S4Vectors::isEmpty(id_B)){stop(group_B, " not found in STORM object META data")}
        DT_A <- add_StartRate_1bpDS(STORM$DATA[[id_A]])
        DT_B <- add_StartRate_1bpDS(STORM$DATA[[id_B]])
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

# Add Start rate log2-Fold-Change 1bp down-stream to STORM object
add_SRlog2FCh1bpDS <- function(STORM, group_A, group_B, newColName = "auto",
                               onNucs = c("A", "C", "G", "T"), minCov = 50){
    if(newColName == "auto"){
        newColName <- paste("SRlog2FCh1bpDS", group_A, group_B, sep = "_")
    }
    sets <- intersect(STORM$META[STORM$META$group == group_A,]$set, STORM$META[STORM$META$group == group_B,]$set)
    if(length(sets) == 0){stop("No sets found with both group_A and group_B label")}
    OUT <- lapply(sets, function(iSet){
        id_A <- STORM$META[STORM$META$set == iSet & STORM$META$group == group_A,]$id
        id_B <- STORM$META[STORM$META$set == iSet & STORM$META$group == group_B,]$id
        if(S4Vectors::isEmpty(id_A)){stop(group_A, " not found in STORM object META data")}
        if(S4Vectors::isEmpty(id_B)){stop(group_B, " not found in STORM object META data")}
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
    sets <- intersect(STORM$META[STORM$META$group == lib_LowdNTPs,]$set,
                      STORM$META[STORM$META$group == lib_HighdNTPs,]$set)
    if(length(sets) == 0){stop("No sets found with both lib_LowdNTPs and lib_HighdNTPs label")}
    OUT <- lapply(sets, function(iSet){
        id_A <- STORM$META[STORM$META$set == iSet & STORM$META$group == lib_LowdNTPs,]$id
        id_B <- STORM$META[STORM$META$set == iSet & STORM$META$group == lib_HighdNTPs,]$id
        if(S4Vectors::isEmpty(id_A)){stop(group_A, " not found in STORM object META data")}
        if(S4Vectors::isEmpty(id_B)){stop(group_B, " not found in STORM object META data")}
        DT <- detect_2OmeScore(DT_low = STORM$DATA[[id_A]],
                               DT_high = STORM$DATA[[id_B]],
                               neighFlankSize = neighFlankSize)
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

# Add Misincorporation Rate Difference (MRD) treatA-treatB to STORM object
add_MRD <- function(STORM, group_A, group_B, newColName = "auto",
                    onNucs = c("A", "C", "G", "T")){
    if(newColName == "auto"){
        newColName <- paste("MRD", group_A, group_B, sep = "_")
    }
    sets <- intersect(STORM$META[STORM$META$group == group_A,]$set, STORM$META[STORM$META$group == group_B,]$set)
    if(length(sets) == 0){stop("No sets found with both group_A and group_B label")}
    OUT <- lapply(sets, function(iSet){
        id_A <- STORM$META[STORM$META$set == iSet & STORM$META$group == group_A,]$id
        id_B <- STORM$META[STORM$META$set == iSet & STORM$META$group == group_B,]$id
        if(S4Vectors::isEmpty(id_A)){stop(group_A, " not found in STORM object META data")}
        if(S4Vectors::isEmpty(id_B)){stop(group_B, " not found in STORM object META data")}
        DT_A <- STORM$DATA[[id_A]] %>% add_diffNucToRefRatio()
        DT_B <- STORM$DATA[[id_B]] %>% add_diffNucToRefRatio()
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
#' @param STORM
#' @param group_A
#' @param group_B
#' @param newColName
#' @param onNucs
#'
#' @return
#' @export
#'
#' @examples
add_CytPer <- function(STORM, group_A, group_B, newColName = "auto",
                       onNucs = c("C")){
    if(newColName == "auto"){
        newColName <- paste("CytPer", group_A, group_B, sep = "_")
    }
    sets <- intersect(STORM$META[STORM$META$group == group_A,]$set, STORM$META[STORM$META$group == group_B,]$set)
    if(length(sets) == 0){stop("No sets found with both group_A and group_B label")}
    OUT <- lapply(sets, function(iSet){
        id_A <- STORM$META[STORM$META$set == iSet & STORM$META$group == group_A,]$id
        id_B <- STORM$META[STORM$META$set == iSet & STORM$META$group == group_B,]$id
        if(S4Vectors::isEmpty(id_A)){stop(group_A, " not found in STORM object META data")}
        if(S4Vectors::isEmpty(id_B)){stop(group_B, " not found in STORM object META data")}
        DT_A <- STORM$DATA[[id_A]] %>% detect_m5C()
        DT_B <- STORM$DATA[[id_B]] %>% detect_m5C()
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

#' Add C->T misincorporation difference (CtoT_MRD) to STORM object
#'
#' @param STORM
#' @param group_A
#' @param group_B
#' @param newColName
#' @param onNucs
#'
#' @return
#' @export
#'
#' @examples
add_CtoT_MRD <- function(STORM, group_A, group_B, newColName = "auto",
                         onNucs = c("C")){
    if(newColName == "auto"){
        newColName <- paste("CtoT.MRD", group_A, group_B, sep = "_")
    }
    sets <- intersect(STORM$META[STORM$META$group == group_A,]$set, STORM$META[STORM$META$group == group_B,]$set)
    if(length(sets) == 0){stop("No sets found with both group_A and group_B label")}
    OUT <- lapply(sets, function(iSet){
        id_A <- STORM$META[STORM$META$set == iSet & STORM$META$group == group_A,]$id
        id_B <- STORM$META[STORM$META$set == iSet & STORM$META$group == group_B,]$id
        if(S4Vectors::isEmpty(id_A)){stop(group_A, " not found in STORM object META data")}
        if(S4Vectors::isEmpty(id_B)){stop(group_B, " not found in STORM object META data")}
        DT_A <- STORM$DATA[[id_A]] %>% detect_ac4C()
        DT_B <- STORM$DATA[[id_B]] %>% detect_ac4C()
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

#' Add deletion rate difference to STORM object
#'
#' @param STORM
#' @param group_A
#' @param group_B
#' @param newColName
#' @param onNucs
#'
#' @return
#' @export
#'
#' @examples
add_DRD <- function(STORM, group_A, group_B, newColName = "auto",
                    onNucs = c("A", "C", "G", "T")){
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
        DT_A <- STORM$DATA[[id_A]] %>% deletion_rate()
        DT_B <- STORM$DATA[[id_B]] %>% deletion_rate()
        tmpRES <- data.table::data.table(DT_A$deletion_Ratio - DT_B$deletion_Ratio)
        tmpRES <- data.table::data.table(iSet, newColName, tmpRES)
        names(tmpRES) <- c("set", "metric", "score")
        tmpRES <- data.table::data.table(DT_A[,c("chr", "gencoor", "strand", "gene",
                                                 "txcoor", "pos", "refSeq")], tmpRES)
        tmpRES[!(tmpRES$refSeq %in% onNucs), "score"] <- NA
        return(tmpRES)
    }) %>% do.call(what = rbind)
    STORM$RES <- brainSTORM:::hlpr_add_REScols(STORM$RES, OUT)
    return(STORM)
}

# One group metrics ############################################################

#' Add Missincorporation rate
#'
#' @param STORM
#' @param group_A
#' @param newColName
#' @param onNucs
#'
#' @return
#' @export
#'
#' @examples
add_MR <- function(STORM, group_A, newColName = "auto", onNucs = c("A", "C", "G", "T")){
    if(newColName == "auto"){
        newColName <- paste("MissIncRate", group_A, sep = "_")
    }
    sets <- STORM$META[STORM$META$group == group_A,]$set
    if(length(sets) == 0){stop("No sets found with both group_A and group_B label")}
    OUT <- lapply(sets, function(iSet){
        id_A <- STORM$META[STORM$META$set == iSet & STORM$META$group == group_A,]$id
        if(S4Vectors::isEmpty(id_A)){stop(group_A, " not found in STORM object META data")}
        DT_A <- add_diffNucToRefRatio(STORM$DATA[[id_A]])
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

#' Add read-start to coverage ratio
#'
#' @param STORM
#' @param group_A
#' @param newColName
#' @param onNucs
#' @param minCov
#'
#' @return
#' @export
#'
#' @examples
add_SR <- function(STORM, group_A, newColName = "auto", onNucs = c("A", "C", "G", "T"), minCov = 50){
    if(newColName == "auto"){
        newColName <- paste("StartRate", group_A, sep = "_")
    }
    sets <- STORM$META[STORM$META$group == group_A,]$set
    if(length(sets) == 0){stop("No sets found with both group_A and group_B label")}
    OUT <- lapply(sets, function(iSet){
        id_A <- STORM$META[STORM$META$set == iSet & STORM$META$group == group_A,]$id
        if(S4Vectors::isEmpty(id_A)){stop(group_A, " not found in STORM object META data")}
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
#' @param STORM
#' @param group_A
#' @param newColName
#' @param onNucs
#'
#' @return
#' @export
#'
#' @examples
add_SR_1bpDS <- function(STORM, group_A, newColName = "auto", onNucs = c("A", "C", "G", "T")){
    if(newColName == "auto"){
        newColName <- paste("StartRate1bpDS", group_A, sep = "_")
    }
    sets <- STORM$META[STORM$META$group == group_A,]$set
    if(length(sets) == 0){stop("No sets found with both group_A and group_B label")}
    OUT <- lapply(sets, function(iSet){
        id_A <- STORM$META[STORM$META$set == iSet & STORM$META$group == group_A,]$id
        if(S4Vectors::isEmpty(id_A)){stop(group_A, " not found in STORM object META data")}
        DT_A <- add_StartRate_1bpDS(STORM$DATA[[id_A]])
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
#' @param STORM
#' @param group_A
#' @param newColName
#' @param onNucs
#' @param minCov
#'
#' @return
#' @export
#'
#' @examples
add_ER <- function(STORM, group_A, newColName = "auto", onNucs = c("A", "C", "G", "T"), minCov = 50){
    if(newColName == "auto"){
        newColName <- paste("EndRate", group_A, sep = "_")
    }
    sets <- STORM$META[STORM$META$group == group_A,]$set
    if(length(sets) == 0){stop("No sets found with both group_A and group_B label")}
    OUT <- lapply(sets, function(iSet){
        id_A <- STORM$META[STORM$META$set == iSet & STORM$META$group == group_A,]$id
        if(S4Vectors::isEmpty(id_A)){stop(group_A, " not found in STORM object META data")}
        DT_A <- add_EndRate(STORM$DATA[[id_A]])
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

#' Add ScoreA-3prime to a STORM object
#'
#' @param STORM
#' @param group_A
#' @param newColName
#' @param onNucs
#' @param minMedCov
#'
#' @return
#' @export
#' @aliases storm_add_scoreA3p
#'
#' @examples
add_scoreA3p <- function(STORM, group_A, newColName = "auto",
                         onNucs = c("A", "C", "G", "T"), minMedCov = 15){
    if(newColName == "auto"){
        newColName <- paste("ScoreA3p", group_A, sep = "_")
    }
    sets <- STORM$META[STORM$META$group == group_A,]$set
    if(length(sets) == 0){stop("No sets found with group_A label")}
    OUT <- lapply(sets, function(iSet){
        id_A <- STORM$META[STORM$META$set == iSet & STORM$META$group == group_A,]$id
        DT_A <- add_scoreA_3p(STORM$DATA[[id_A]], minMedEnv = minMedCov)
        if(S4Vectors::isEmpty(id_A)){stop(group_A, " not found in STORM object META data")}
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

#' Add ScoreA-3prime to a STORM object
#'
#' @param STORM
#' @param group_A
#' @param newColName
#' @param onNucs
#' @param minMedCov
#'
#' @return
#' @export
#'
#' @examples
#' @aliases storm_add_scoreA5p
add_scoreA5p <- function(STORM, group_A, newColName = "auto",
                         onNucs = c("A", "C", "G", "T"), minMedCov = 15){
    if(newColName == "auto"){
        newColName <- paste("ScoreA5p", group_A, sep = "_")
    }
    sets <- STORM$META[STORM$META$group == group_A,]$set
    if(length(sets) == 0){stop("No sets found with both group_A and group_B label")}
    OUT <- lapply(sets, function(iSet){
        id_A <- STORM$META[STORM$META$set == iSet & STORM$META$group == group_A,]$id
        if(S4Vectors::isEmpty(id_A)){stop(group_A, " not found in STORM object META data")}
        DT_A <- add_scoreA_5p(STORM$DATA[[id_A]], minMedEnv = minMedCov)
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

# Helper functions #############################################################
# Function for cytidine persistence to Bisulphite treatment
detect_m5C <- function(DT){
    tmp <- round(DT$C / (DT$T + DT$C), 4)
    tibble::add_column(DT, det_m5C = tmp)
}

# Function for cytidine persistence to Bisulphite treatment
detect_ac4C <- function(DT){
    DT <- data.table::data.table(DT)
    tmp <- round(DT$`T` / (DT$`T` + DT$`C`), 6)
    tibble::add_column(DT, det_ac4C = tmp)
}

# Add column of Different Nucleotides to reference ratio, diffToRef and nucTotal columns are required
add_diffNucToRefRatio <- function(DT){
    DT <- txtools::tx_add_diffNucToRef(DT) %>% txtools::tx_add_nucTotal()
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

# Start rate 1bp DS with pseudoCount
add_StartRate_1bpDS_pc <- function(DT, minCov = 50, pc = 0.01){
    tmp <- (DT$start_5p + pc) / (DT$cov + pc)
    tmp[DT$cov < minCov] <- NA
    DTL <- tibble::add_column(DT, startRate_1bpDS_pc = tmp) %>% txtools::tx_split_DT()
    DT <- lapply(DTL, function(DT){
        DT$startRate_1bpDS_pc <- c(utils::tail(DT$startRate_1bpDS_pc, -1), NA)
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
        tmpSR <- c(utils::tail(DTL1[[iG]]$start_Ratio / DTL2[[iG]]$start_Ratio, -1), NA)
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

# Score A Birkedal et al., 2015
scoreA <- function(CHR, d = 6, minMedEnv = 15){
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
add_scoreA_3p <- function(DT, minMedEnv = 15, nCores = 1){
    DTL <- txtools::tx_split_DT(DT)
    parallel::mclapply(mc.cores = nCores, DTL, function(x){
        tmp <- scoreA(x$end_3p, minMedEnv = minMedEnv)
        tibble::add_column(x, scoreA_3p = tmp)
    }) %>% txtools::tx_merge_DT()
}

# Adding score A 5p_ends to multiGene DT
add_scoreA_5p <- function(DT, minMedEnv = 15, nCores = 1){
    DTL <- txtools::tx_split_DT(DT)
    parallel::mclapply(mc.cores = nCores, DTL, function(x){
        if(!("startRate_1bpDS" %in% names(x))){
            x <- add_StartRate_1bpDS(x)
        }
        tmp <- scoreA(x$start_5p, minMedEnv = minMedEnv)
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
deletion_rate <- function(DT, minReadCounts = 50){
    DT <- txtools::tx_add_nucTotal(DT)
    tmp <- round(DT$`-` / DT$nucTotal, 6)
    tmp[DT$nucTotal < minReadCounts] <- NA
    tibble::add_column(DT, deletion_Ratio = tmp)
}

# Add column of Different Nucleotides to reference ratio, diffToRef and nucTotal columns are required
add_diffNucToRefRatio <- function(DT){
    DT <- txtools::tx_add_diffNucToRef(DT) %>% txtools::tx_add_nucTotal()
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
