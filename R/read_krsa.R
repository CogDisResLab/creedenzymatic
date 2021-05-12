#' Reads and Rank KRSA table
#'
#' reads KRSA table and checks for correct format
#'
#' This function takes in table and rank and quartile kinases based on the absolute Score values
#'
#' @param df dataframe, table output (requires at least Kinase and Score columns)
#' @param ..., arguments passed to rank_kinases function
#'
#' @return dataframe, Ranked and quartiled table
#'
#' @export
#'

read_krsa <- function(df, ...) {


  if(!is.data.frame(df)) {stop("Make sure your input is a dataframe")}
  if(!all(colnames(df) %in% c("Kinase", "Score"))) {
    stop("check columns names, they should be: Kinase, Score")
  }

  if(!is.character(df$Kinase)) {stop("check that Kinase column is as character column")}
  if(!is.numeric(df$Score)) {stop("check that Score column is as numeric column")}

  rank_kinases(df, tool = "KRSA", ...)



}
