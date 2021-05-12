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

combine_tools <- function(KRSA_df = NULL, UKA_df = NULL, KEA3_df = NULL, PTM_SEA_df = NULL, mapping_df = kinome_mp_file_v1) {


  my_tibble <- tibble::tibble(
    Kinase = character(),
    Score = numeric(),
    Rank = numeric(),
    Method = character(),
    Perc = numeric(),
    Qrt = numeric(),
    Uniprot_Gene = character(),
    KinaseFamily = character(),

  )

  if(!is.null(KRSA_df)) {
    dplyr::left_join(KRSA_df, dplyr::select(mapping_df, KRSA, Uniprot_Gene), by = c("Kinase" = "KRSA")) %>%
      dplyr::mutate(KinaseFamily = Kinase) %>%
      rbind(my_tibble) -> my_tibble
  }

  if(!is.null(UKA_df)) {
    dplyr::left_join(UKA_df, dplyr::select(mapping_df, UKA, Uniprot_Gene, KRSA), by = c("Kinase" = "UKA")) %>%
      dplyr::rename(KinaseFamily = KRSA) %>%
      rbind(my_tibble) -> my_tibble
  }

  if(!is.null(KEA3_df)) {
    dplyr::left_join(KEA3_df, dplyr::select(mapping_df, KEA3, Uniprot_Gene, KRSA), by = c("Kinase" = "KEA3")) %>%
      dplyr::rename(KinaseFamily = KRSA) %>%
      rbind(my_tibble) -> my_tibble
  }

  if(!is.null(PTM_SEA_df)) {
    dplyr::left_join(PTM_SEA_df, dplyr::select(mapping_df, PTMSEA, Uniprot_Gene, KRSA), by = c("Kinase" = "PTMSEA")) %>%
      dplyr::rename(KinaseFamily = KRSA) %>%
      rbind(my_tibble) -> my_tibble
  }

  my_tibble



}
