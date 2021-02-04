# mapping Summary Report PDF plots
mapSummPlot <- function(x){
    mapSumm <- data.table::fread(file = file.path(outDir, "mappingSummary.txt"),
                                 header = TRUE, stringsAsFactors = FALSE) %>%
        data.frame
    rownames(mapSumm) <- as.character(mapSumm[,1]) %>% gsub(pattern = " \\|", replacement = "")
    mapSumm <- t(mapSumm[,-1])
    colnames(mapSumm)
}
