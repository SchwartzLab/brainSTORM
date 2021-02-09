SSIIIandTGIRT_metrics <- function(STORM){
    STORM %>%
        add_Y_metrics %>%
        add_Nm_metrics %>%
        add_ac4C_metrics %>%
        add_m1A_metrics %>%
        add_m7G_metrics %>%
        add_m5C_metrics
}
