#' STORM META constructor from RDS file names
#'
#' @param DTfiles
#' @param varsList
#' @param groupVars
#' @param setVars
#' @param idVars
#'
#' @return
#' @export
#'
#' @examples
storm_META <- function(DTfiles, varsList, groupVars, setVars, idVars){
    if(not("libTreat" %in% names(varsList))){
        stop("libTreat must be one of the element names in varsList argument")
    }
    DT <- data.table(RDS = DTfiles)
    procFileNames <- DT$RDS %>% str_split(pattern = "/") %>%
        lapply(function(x) tail(x, 1)) %>% unlist
    DT <- lapply(names(varsList), function(i){
        patt <- paste0("(", paste(varsList[[i]], collapse = "|"), ")")
        tmp <- stringr::str_extract_all(string = procFileNames, pattern = patt)
        lapply(tmp, function(x) paste(x, collapse = ".")) %>% unlist %>%
            data.table() %>% set_names(i)
    }) %>% do.call(what = cbind) %>% cbind(DT, .)
    DT$group <- DT[,groupVars, with = F] %>%
        apply(1, function(x) paste(x, collapse= ".")) %>%
        factor()
    DT$set <- DT[,setVars, with = F] %>%
        apply(1, function(x) paste(x, collapse= ".")) %>%
        factor()
    DT$id <- DT[,idVars, with = F] %>%
        apply(1, function(x) paste(x, collapse= "."))
    if(sum(duplicated(DT$id)) !=0){
        warning("Generated ids are not unique per sample")
    }
    if(sum(duplicated(paste(DT$group, DT$set)))){
        warning("Combinations of group and set should be unique between samples")
    }
    return(DT)
}
