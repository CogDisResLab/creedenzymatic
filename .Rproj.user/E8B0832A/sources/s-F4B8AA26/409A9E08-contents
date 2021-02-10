#' Plot Quartile Figure
#'
#' Takes KRSA and UKA ranked tables and generate a quartile figure
#'
#' This function takes KRSA and UKA ranked tables and generate a quartile figure
#'
#' @param df1 dataframe, UKA ranked table
#' @param df2 dataframe, KRSA ranked table
#'
#' @return ggplot figure
#'

quartile_figure <- function(df, df2) {


  d1 <- rbind(df, df2)

  d1 %>% pivot_longer(3:ncol(.), names_to = "Tool", values_to = "Quartile") %>%
    mutate(
           present = ifelse(is.na(Quartile), "1", "2"), Quartile = ifelse(present == "1", 2, Quartile)
    ) %>%
    mutate(Tool = gsub("_Quartile", "", Tool ), Family = gsub(" Family", "", `Kinase Family`)) %>%
    ggplot(aes(reorder(Kinase, Family), Tool)) + geom_point(aes(size = Quartile)) +
    coord_flip() +
    scale_size(trans = "reverse", range = c(1,5), breaks = c(1,2,3,4)) +
    #geom_text(position = position_dodge(width = 1), aes(x=`Kinase Family`, y=5, label = `Kinase Family`)) +
    theme_bw() + facet_grid(Family ~ ., scales = "free", space = "free") +
    labs(x = "", y = "")



}
