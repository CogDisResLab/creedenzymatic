#' Plot quartile Figure
#'
#' Takes the combined ranked dataframe (KRSA, UKA, .. etc) and generate a quartile figure
#'
#'
#' @param df dataframe, combined mapped tables
#' @param grouping character to choose grouping (KinaseFamily, subfamily, or group). Default is KinaseFamily
#'
#' @return ggplot figure
#'
#' @export
#'

quartile_figure <- function(df, grouping = "KinaseFamily") {

  df %>%
    dplyr::select(hgnc_symbol, one_of(grouping), Qrt, Method) %>%
    tidyr::pivot_wider(names_from = Method, values_from = Qrt) %>%
    tidyr::pivot_longer(3:ncol(.), names_to = "Method", values_to = "Qrt", values_fn = unique) %>%
    dplyr::mutate(
      present = ifelse(is.na(Qrt), "No", "Yes"), Qrt = ifelse(present == "No", 2, Qrt),
      present = factor(present, levels = c("Yes", "No")),
      Qrt = factor(Qrt, levels = c(1,2,3,4)),
      Method = factor(Method, levels = c("UKA", "PTM-SEA", "KEA3", "KRSA"))
    ) %>%
    ggplot2::ggplot(ggplot2::aes(hgnc_symbol, Method)) +
    ggplot2::geom_point(ggplot2::aes(size = Qrt, shape = present)) +
    #gplot2::scale_radius(trans = "reverse") +
    ggplot2::scale_size_manual(values = c("4" = 4, "3" = 3, "2" = 2, "1" = 1)) +
    ggplot2::theme_bw() +
    {if(grouping == "subfamily") {ggplot2::facet_grid(. ~ subfamily, scales = "free", space = "free")}
      else if(grouping == "group") {ggplot2::facet_grid(. ~ group, scales = "free", space = "free")}
      else {ggplot2::facet_grid(. ~ KinaseFamily, scales = "free", space = "free")}
      }+
    ggplot2::scale_shape_manual(values=c(Yes = 19, No = 1)) +
    ggplot2::theme(axis.text.x = element_text(angle = 30, size = 7.5, vjust = 0.7)) +
    ggplot2::labs(x = "", y = "")



}
