#' @export
storm_add_scoreA5p <- function(STORM, group_A, newColName = "auto",
                               onNucs = c("A", "C", "G", "T"), minMedCov = 15){
    # .Deprecated("add_scoreA5p")
    add_scoreA5p(STORM, group_A, newColName, onNucs, minMedCov)
}

#' @export
storm_add_scoreA3p <- function(STORM, group_A, newColName = "auto",
                               onNucs = c("A", "C", "G", "T"), minMedCov = 15){
    # .Deprecated("add_scoreA3p")
    add_scoreA3p(STORM, group_A, newColName, onNucs, minMedCov)
}
