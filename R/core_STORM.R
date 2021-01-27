# STORM object constructor from META table
storm_STORM <- function(META, fastaGenome = NULL, geneAnnot = NULL, nCores = 1){
    if(hasDups(META$id)){
        stop("Not allowed duplicated id variables in META")
    }
    DTL <- lapply(META$RDS, function(x) readRDS(x)) # load RDS files
    # Check if data is in data.table or list format
    if(all(lapply(DTL, class) %>% unlist %>% equals("data.frame"))){
        DTL <- lapply(DTL, data.table) %>% set_names(META$id)
    }else if(all(lapply(DTL, class) %>% unlist %>% equals("list"))){
        DTL <- lapply(DTL, function(x){
            do.call(x, what = rbind) %>% data.table()
        }) %>% set_names(META$id)
    }
    # Have reference sequence, if not add it
    reqRefSeq <- lapply(DTL, function(x){"refSeq" %in% names(x)}) %>% unlist %>% not
    if(sum(reqRefSeq) > 0 & is.null(fastaGenome) | is.null(geneAnnot)){
        stop("Data requires reference sequence, fastaGenome and geneAnnot arguments must be provided")
    }
    if(sum(reqRefSeq) > 0){
        DTL[reqRefSeq] <- mclapply(mc.cores = nCores, DTL[reqRefSeq], function(DT){
            tx_split_DT(DT) %>% lapply(function(x){
                tx_add_refSeqDT(DT = x, genome = fastaGenome, geneAnnot = geneAnnot)
            }) %>% tx_merge_DT()
        })
    }
    # Check uniformity of DTs, if unequal equalize
    tmpPos <- lapply(DTL, function(DT){paste(DT$gene, DT$txcoor, sep = ":")})
    identCoors <- lapply(tmpPos[-1], function(x){identical(tmpPos[[1]], x)}) %>% unlist %>% all
    if(!identCoors){
        if(is.null(geneAnnot) | is.null(fastaGenome)){
            stop("Data needs compatibility adjustment, this requires geneAnnot and
                 fastaGenome arguments to be provided")
        }
        tmpGenes <- lapply(DTL, function(x) as.character(x$gene)) %>% unlist %>% unique
        tmpGA <- geneAnnot[match(tmpGenes, geneAnnot$name)]
        DTL <- lapply(DTL, function(DT){
            complete_DT(DT, tmpGA, fastaGenome, nCores)
        })
    }
    if(!("pos" %in% names(DTL[[1]]))){
        DTL <- lapply(DTL, function(x) tx_add_pos(x))
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
