#' Plot quartile Figure
#'
#' Takes the combined ranked dataframe (KRSA, UKA, .. etc) and generate a quartile figure
#'
#'
#' @param df dataframe, combined mapped tables
#'
#' @return ggplot figure
#'
#' @export
#'

quartile_figure <- function(df) {

  df %>%
    dplyr::select(Uniprot_Gene, KinaseFamily, Qrt, Method) %>%
    tidyr::pivot_wider(names_from = Method, values_from = Qrt) %>%
    tidyr::pivot_longer(3:4, names_to = "Method", values_to = "Qrt") %>%
    dplyr::mutate(
      present = ifelse(is.na(Qrt), "No", "Yes"), Qrt = ifelse(present == "No", 2, Qrt)
    ) %>%
    #filter(KinaseFamily %in% kinases) %>%
    ggplot2::ggplot(ggplot2::aes(Uniprot_Gene, Method)) +
    ggplot2::geom_point(ggplot2::aes(size = Qrt, shape = present)) +
    ggplot2::scale_radius(trans = "reverse") +
    ggplot2::theme_bw() +
    ggplot2::facet_grid(. ~ KinaseFamily, scales = "free", space = "free") +
    ggplot2::scale_shape_manual(values=c(1, 19)) +
    ggplot2::theme(axis.text.x = element_text(angle = 30, size = 7.5, vjust = 0.7)) +
    ggplot2::labs(x = "", y = "")



}
