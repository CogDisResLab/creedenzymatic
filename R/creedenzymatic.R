#' Runs Creedenzymatic
#'
#' reads KRSA, UKA, LFC tables and run creedenzymatic
#'
#' This function takes in table and rank and quartile kinases based on the absolute Score values
#'
#' @param KRSA_table dataframe, KRSA table output
#' @param UKA_table dataframe, UKA table output
#' @param LFC_table dataframe, KEA table output
#' @param avg_krsa
#' @param avg_lfc
#' @param ..., arguments passed to other functions
#'
#' @return dataframe, Ranked and quartiled table
#'
#' @export
#'

creedenzymatic <- function(KRSA_table, UKA_table, LFC_table, avg_krsa = T, avg_lfc = T, prefix = "Comp1", ...) {

  suppressMessages(read_delim(KRSA_table, delim = "\t")) %>%
    {if(avg_krsa) {dplyr::select(., Kinase, Score = AvgZ) %>% dplyr::distinct(.keep_all = T)}
  else {dplyr::select(., Kinase, Score = Z)}} %>%
    read_krsa(trns = "abs", sort = "desc") -> krsa_table_ranked

  suppressMessages(read_delim(UKA_table, delim = "\t")) %>%
    select(`Kinase Name`, `Median Final score`) %>%
    rename(Kinase = `Kinase Name`, Score = `Median Final score`) %>%
    read_uka(trns = "abs", sort = "desc") -> uka_table_ranked


  suppressMessages(read_delim(LFC_table, delim = "\t")) %>%
    {if(avg_lfc) {dplyr::select(.,Peptide, Score = totalMeanLFC)} else {dplyr::select(.,Peptide, Score = LFC)}} %>%
    distinct(.keep_all = T) -> lfc_df

  read_kea(lfc_df, filter = T, cutoff_abs = T, sort = "asc", trns = "abs",
           rm_duplicates = T, method = "MeanRank", lib = "kinase-substrate", ...) -> kea_table_ranked


  read_ptmsea(lfc_df) -> ptm_table_ranked

  combine_tools(KRSA_df = krsa_table_ranked, UKA_df = uka_table_ranked, KEA3_df = kea_table_ranked,
                PTM_SEA_df = ptm_table_ranked) -> combined_df_full

  readr::write_delim(combined_df_full, paste0("output_tables/",prefix, "_", "Combined_Table.txt"), delim = "\t")

  combined_df_full


}
