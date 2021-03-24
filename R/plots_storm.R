#' Boxplots of $RES calculated metrics grouping by modified nucleotides
#'
#' @param STORM list. STORM object as output by \code{\link{storm_STORM}}
#' @param metrics character. Selected metrics to plot. By default plots all metrics
#' @param title character. Title to be shown on top of plot
#'
#' @return
#' @export
#'
#' @examples
storm_metricsBoxPlot_byNuc <- function(STORM, metrics = NULL, title = ""){
    if(!"nuc" %in% names(STORM$RES)){stop("'nuc' column with RNAmod identity is not present in STORM$RES")}
    if(!is.null(metrics)){
        tmpDT <- dplyr::filter(STORM$RES, metric %in% metrics) %>% stats::na.omit()
    }else{
        tmpDT <- STORM$RES %>% stats::na.omit()
    }
    ggplot2::ggplot(tmpDT, ggplot2::aes(x = nuc, y = score)) + ggplot2::geom_boxplot(outlier.colour = NA) +
        ggplot2::geom_point(ggplot2::aes(colour = metric), alpha = 0.2) +
        ggplot2::theme_bw() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) +
        ggplot2::facet_wrap(~metric, scales = "free") +
        ggplot2::ggtitle(title)
}

#' STORM scatterplot of scores by txcoor
#'
#' Plot scatterplots of scores in STORM$RES by txcoor and colored by gene
#'
#' @param STORM list. STORM object as output by \code{\link{storm_STORM}}
#' @param metrics character. Selected metrics to plot. By default plots all metrics
#' @param title character. Title to be shown on top of plot
#'
#' @return
#' @export
#'
#' @examples
storm_metricsScatterPlot_byPos <- function(STORM, metrics = NULL, title = "", na.rm = FALSE){
    if(!is.null(metrics)){
        tmpDT <- dplyr::filter(STORM$RES, metric %in% metrics)
    }else{
        tmpDT <- STORM$RES
    }
    if(na.rm){tmpDT <- stats::na.omit(tmpDT)}
    ggplot2::ggplot(tmpDT) +
        ggplot2::geom_point(ggplot2::aes(x = .data$pos, y = .data$score, colour = .data$gene), alpha = 0.2) +
        ggplot2::facet_grid(.data$set ~ .data$metric, scales = "free") +
        ggplot2::theme_bw() +
        ggplot2::theme(axis.text.x = ggplot2::element_blank()) +
        ggplot2::ggtitle(title) + ggplot2::xlab("txCoor") +
        ggplot2::theme_minimal()
}
