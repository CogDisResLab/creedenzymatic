#' Reads a dataframe of Peptides IDs and their Scores and run KEA3
#'
#' reads a dataframe of Peptides IDs and their Scores (LFC, p-value, ... etc) and run KEA3 on a subset of these peptides or all of them
#'
#' This function a dataframe of Peptides IDs and their Scores (LFC, p-value, ... etc) and run KEA3 on a subset of these peptides or all of them
#'
#' @param df dataframe, must have at least Peptide and Score columns
#' @param filter boolean to subset peptides or not
#' @param cutoff numeric to act as the cutoff to filter out peptides
#' @param cutoff_abs boolean (use absolute value or not) default is TRUE
#' @param direction ("lower", "higher) filter based on lower than or higher than the cutoff values (default to "higher")
#' @param rm_duplicates boolean (TRUE or FALSE) remove genes duplicates
#' @param method "MeanRank" takes the mean rank across all libraries or "MeanFDR" takes the mean of FDR across all libraries (default is "MeanRank")
#' @param lib searched kea libraries "kinase-substrate" or "all" (default is "kinase-substrate" which will return only kinase libraries like ChengKSIN, PTMsigDB, PhosDAll)

#' @param ..., arguments passed to rank_kinases function
#'
#' @return dataframe, Ranked and quartiled table
#'
#' @export
#'

read_kea <- function(df, filter = T, cutoff = 0.2, cutoff_abs = T, direction = "higher", rm_duplicates = T, method = "MeanRank" , lib = c("kinase-substrate"), ...) {


  if(!is.data.frame(df)) {stop("Make sure your input is a dataframe")}
  if(!all(colnames(df) %in% c("Peptide", "Score"))) {
    stop("check columns names, they should be: Peptide, Score")
  }

  if(!is.character(df$Peptide)) {stop("check that Kinase column is as character column")}
  if(!is.numeric(df$Score)) {stop("check that Score column is as numeric column")}

  if(filter == T) {

    df %>% dplyr::filter(
      if(cutoff_abs == T) {

        if(direction == "higher") {
          abs(Score) >= cutoff
        }

        else if (direction == "lower") {
          abs(Score) <= cutoff
        }


      }

      else {
        if(direction == "higher") {
          Score >= cutoff
        }

        else if (direction == "lower") {
          Score <= cutoff
        }
      }

      ) %>%
      dplyr::distinct(Peptide, .keep_all = T) -> df

  }

  else {
    df %>%
    dplyr::distinct(Peptide, .keep_all = T) -> df

  }

  df %>%
    dplyr::left_join(rbind(stk_pamchip_87102_mapping, ptk_pamchip_86402_mapping), by = c("Peptide" = "ID")) %>%
    dplyr::filter(!is.na(HGNC)) %>%
    dplyr::pull(HGNC) -> pep_hits

  if(rm_duplicates == T) {
    message(paste0("Removed ", length(setdiff(pep_hits, unique(pep_hits))), " gene duplicates"))
    pep_hits <- unique(pep_hits)
  }


  message(paste0("Total # of genes ", length(pep_hits)))

  run_kea(pep_hits, lib) -> kea_results

  if(!all(is.na(kea_results))) {

    if(method == "MeanRank") {
      kea_results <- kea_results$`Integrated--meanRank` %>%
        select(TF, Score) %>%
        rename(Kinase = TF)
    }

    else {
      purrr::map_df(kea_results, base::rbind) %>%
        dplyr::select(TF, FDR) %>%
        dplyr::mutate(FDR = as.numeric(FDR)) %>%
        dplyr::group_by(TF) %>%
        dplyr::summarise(Score = mean(FDR)) %>%
        dplyr::ungroup() %>%
        dplyr::rename(Kinase = TF) -> kea_results
    }

    rank_kinases(kea_results, tool = "KEA3", ...)
  }

  else {
    message("Couldn't connect to KEA3 API successfully")
  }



}
