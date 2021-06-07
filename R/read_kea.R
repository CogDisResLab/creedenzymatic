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
#' @param direction ("lower", "higher) filter based on less than or bigger than the cutoff values (default to "higher")
#' @param lib searched kea libraries (default is "kinases" which will return only kinase libraries like ChengKSIN, PTMsigDB, PhosDAll)

#' @param ..., arguments passed to rank_kinases function
#'
#' @return dataframe, Ranked and quartiled table
#'
#' @export
#'

read_kea <- function(df, filter, cutoff, cutoff_abs = T,direction = "higher", lib = c("kinases"), ...) {


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
    pull(HGNC) -> pep_hits

  pep_hits <- na.omit(pep_hits)

  message(length(pep_hits))

  run_kea(pep_hits, lib) -> kea_results

  if(!all(is.na(kea_results))) {
    purrr::map_df(kea_results, base::rbind) %>%
      dplyr::select(TF, FDR) %>%
      dplyr::mutate(FDR = as.numeric(FDR)) %>%
      dplyr::group_by(TF) %>%
      dplyr::summarise(Score = mean(FDR)) %>%
      dplyr::ungroup() %>%
      dplyr::rename(Kinase = TF) -> kea_results

    rank_kinases(kea_results, tool = "KEA3", ...)
  }

  else {
    message("Couldn't connect to KEA3 API successfully")
  }



}
