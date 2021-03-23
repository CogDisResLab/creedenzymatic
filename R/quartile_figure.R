#' Plot Quartile Figure
#'
#' Takes KRSA and UKA ranked tables and generate a quartile figure
#'
#' This function takes KRSA and UKA ranked tables and generate a quartile figure
#'
#' @param df dataframe, combined mapped tables
#' @param kinases dataframe, KRSA ranked table
#'
#' @return ggplot figure
#'

quartile_figure <- function(df, kinases) {


  df %>%
    select(Uniprot_Gene, KinaseFamily, Qrt, Method) %>%
    pivot_wider(names_from = Method, values_from = Qrt) %>%
    pivot_longer(3:4, names_to = "Method", values_to = "Qrt") %>%

    mutate(
      present = ifelse(is.na(Qrt), "1", "2"), Qrt = ifelse(present == "1", 2, Qrt)
    ) %>%

    filter(KinaseFamily %in% kinases) %>%
    ggplot(aes(reorder(Uniprot_Gene, KinaseFamily), Method)) + geom_point(aes(size = Qrt, shape = present)) +

    scale_radius(trans = "reverse") +
    theme_bw() + facet_grid(. ~ KinaseFamily, scales = "free", space = "free") +
    scale_shape_manual(values=c(1, 19)) +
    theme(axis.text.x = element_text(angle = 30, size = 7.5, vjust = 0.7)) +
    labs(x = "", y = "")



}
