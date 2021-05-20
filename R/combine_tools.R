#' Combine data for quartile figure
#'
#' reads ranked tables from the different tools (KRSA, UKA, ... etc)
#'
#' This function takes in ranked tables from the different tools (KRSA, UKA, ... etc) and map them to the kinome mapping file and return df ready for the quartile figure
#'
#' @param mapping_df kinome mapping df (default is kinome_mp_file_v1)
#' @param KRSA_df dataframe, KRSA table output (requires at least Kinase and Score columns)
#' @param UKA_df dataframe, UKA table output (requires at least Kinase and Score columns)
#' @param KEA3_df dataframe, KEA table output (requires at least Kinase and Score columns)
#' @param PTM_SEA_df dataframe, PTM_SEA table output (requires at least Kinase and Score columns)
#'
#' @return dataframe, ready for quartile figure
#'
#' @export
#'

combine_tools <- function(KRSA_df = NULL, UKA_df = NULL, KEA3_df = NULL, PTM_SEA_df = NULL, mapping_df = kinome_mp_file_v2) {


  my_tibble <- tibble::tibble(
    Kinase = character(),
    Score = numeric(),
    Rank = numeric(),
    Method = character(),
    Perc = numeric(),
    Qrt = numeric(),
    hgnc_symbol = character(),
    group = character(),
    KinaseFamily = character(),
    subfamily = character()

  )

  if(!is.null(KRSA_df)) {
    dplyr::left_join(KRSA_df %>% mutate(Kinase = toupper(Kinase)), dplyr::select(mapping_df, krsa_id, hgnc_symbol, group, family, subfamily), by = c("Kinase" = "krsa_id")) %>%
      dplyr::rename(KinaseFamily = family) %>%
      rbind(my_tibble) -> my_tibble
  }

  if(!is.null(UKA_df)) {
    dplyr::left_join(UKA_df %>% mutate(Kinase = toupper(Kinase)), dplyr::select(mapping_df, uka_id, hgnc_symbol, group, family, subfamily), by = c("Kinase" = "uka_id")) %>%
      dplyr::rename(KinaseFamily = family) %>%
      rbind(my_tibble) -> my_tibble
  }

  if(!is.null(KEA3_df)) {
    dplyr::left_join(KEA3_df %>% mutate(Kinase = toupper(Kinase)), dplyr::select(mapping_df, kea3_id, hgnc_symbol, group, family, subfamily), by = c("Kinase" = "kea3_id")) %>%
      dplyr::rename(KinaseFamily = family) %>%
      rbind(my_tibble) -> my_tibble
  }

  if(!is.null(PTM_SEA_df)) {
    dplyr::left_join(PTM_SEA_df %>% mutate(Kinase = toupper(Kinase)), dplyr::select(mapping_df, ptmsea_id, hgnc_symbol, group, family, subfamily), by = c("Kinase" = "ptmsea_id")) %>%
      dplyr::rename(KinaseFamily = family) %>%
      rbind(my_tibble) -> my_tibble
  }

  my_tibble



}
