#' Reads a dataframe of Peptides IDs and their Scores and run PTM-SEA
#'
#' reads a dataframe of Peptides IDs and their Scores (LFC, p-value, ... etc) and run PTM-SEA
#'
#' This function a dataframe of Peptides IDs and their Scores (LFC, p-value, ... etc) and run PTM-SEA
#'
#' @param df dataframe, must have at least Peptide and Score columns
#' @param lib searched PTM-SEA libraries "kinase-substrate" or "all" (default is "kinase-substrate" which will return only kinase libraries like ChengKSIN, PTMsigDB, PhosDAll)

#' @param ..., arguments passed to run ptm-sea function
#'
#' @return dataframe, Ranked and quartiled table
#'
#' @export
#'

read_ptmsea <- function(df, ...) {


  if(!is.data.frame(df)) {stop("Make sure your input is a dataframe")}
  if(!all(colnames(df) %in% c("Peptide", "Score"))) {
    stop("check columns names, they should be: Peptide, Score")
  }

  if(!is.character(df$Peptide)) {stop("check that Kinase column is as character column")}
  if(!is.numeric(df$Score)) {stop("check that Score column is as numeric column")}

  df %>% dplyr::mutate(Peptide = gsub("\\s+", "", Peptide)) %>%
    dplyr::left_join(rbind(stk_pamchip_87102_array_layout_ptmsea, ptk_pamchip_86402_array_layout_ptmsea),
                     by = c("Peptide" = "ID")) %>%
    dplyr::select(PTM_SEA_ID, Score) %>% distinct(PTM_SEA_ID, .keep_all = T) %>%
    dplyr::filter(!is.na(PTM_SEA_ID)) %>%
    tibble::column_to_rownames("PTM_SEA_ID") %>%
    as.matrix() -> mat_ds

  (my_new_ds <- new("GCT", mat=mat_ds))

  suppressWarnings(run_ptmsea(my_new_ds, ...) -> ptmsea_results)

  if(!all(is.na(ptmsea_results))) {

    ptmsea_results %>% lapply(`[[`, 1) %>%
      lapply(`[[`, 1) %>% as.data.frame() %>%
      pivot_longer(cols = 1:ncol(.), names_to = "Set", values_to = "ES") %>%
      #mutate(Kinase = str_extract(Set, "_.+"),
      #       Kinase = str_remove(Kinase, "_"),
      #       Kinase = str_replace(Kinase, "\\.", "/")
      dplyr::mutate(Kinase = Set %>% toupper()) %>%
      dplyr::select(Kinase, Score = ES) -> ptmsea_results

    rank_kinases(ptmsea_results, tool = "PTM-SEA", trns = "abs", sort = "desc")
  }

  else {
    message("Couldn't run PTM_SEA successfully")
  }



}
