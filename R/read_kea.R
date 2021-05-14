#' Reads a dataframe of Peptides IDs and their Scores and run KEA3
#'
#' reads a dataframe of Peptides IDs and their Scores (LFC, p-value, ... etc) and run KEA3 on a subset of these peptides or all of them
#'
#' This function a dataframe of Peptides IDs and their Scores (LFC, p-value, ... etc) and run KEA3 on a subset of these peptides or all of them
#'
#' @param df dataframe, must have at least Peptide and Score columns
#' @param filter boolean to subset peptides or not
#' @param cutoff numeric to act as the cutoff to filter out peptides
#' @param cutoff_abs boolean (use absolute value or not) default to use
#' @param direction ("lower", "higher) filter based on less than or bigger than the cutoff values (default to "higher")
#' @param ..., arguments passed to rank_kinases function
#'
#' @return dataframe, Ranked and quartiled table
#'
#' @export
#'

read_kea <- function(df, ...) {


  if(!is.data.frame(df)) {stop("Make sure your input is a dataframe")}
  if(!all(colnames(df) %in% c("Peptide", "Score"))) {
    stop("check columns names, they should be: Peptide, Score")
  }

  if(!is.character(df$Kinase)) {stop("check that Kinase column is as character column")}
  if(!is.numeric(df$Score)) {stop("check that Score column is as numeric column")}

  if(filter == T) {

    df %>% dplyr::filter(
      if(cutoff_abs == T) {

        if(direction == "higher") {
          abs(score) >= cutoff
        }

        else if (direction == "lower") {
          abs(score) <= cutoff
        }


      }

      else {
        if(direction == "higher") {
          score >= cutoff
        }

        else if (direction == "lower") {
          score <= cutoff
        }
      }

      ) %>%
      dplyr::distinct() -> df

  }

  else {
    df %>%
    dplyr::distinct() -> df

  }

  df %>%
    dplyr::left_join(rbind(stk_pamchip_87102_mapping, ptk_pamchip_86402_mapping), by = c("Peptide" = "ID")) %>%
    pull(HGNC) -> pep_hits

  run_kea(pep_hits, )

  rank_kinases(df, tool = "KRSA", ...)



}
