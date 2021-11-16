#' Extract Top Kinases
#'
#' reads combined dataframe (ranked and quartiled) and extracts top kinases based on adjustable criteria
#'
#' This function takes in the combined dataframe (ranked and quartiled) and extracts top kinases based on adjustable criteria
#'
#' @param combined_df dataframe, Ranked and quartiled dataframe
#' @param min_qrt integer, minimum quartile to count
#' @param min_counts integer, number of minimum hits
#'
#' @return vector, top kinases
#'
#' @export
#'

extract_top_kinases <- function(combined_df, min_qrt, min_counts) {

  combined_df %>% select(hgnc_symbol, Qrt, Method) %>%
    distinct(hgnc_symbol, Method, .keep_all = T) %>%
    pivot_wider(names_from = Method, values_from = Qrt) %>%
    mutate(Avg = rowSums(across(where(is.numeric)) >= min_qrt, na.rm = T)) %>%
    filter(Avg >= min_counts) %>% pull(hgnc_symbol) -> top_kinases

  message(paste0("# of Kinases: ", length(top_kinases)))

  top_kinases


}
