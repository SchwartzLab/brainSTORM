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
storm_META <- function(fileNames, varsList, groupVars, setVars, idVars){
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
    return(DT)
}

