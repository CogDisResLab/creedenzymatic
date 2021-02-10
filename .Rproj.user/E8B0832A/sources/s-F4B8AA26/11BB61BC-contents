#' Reads and Rank UKA table
#'
#' reads UKA table and checks for correct format
#'
#' This function takes in UKA table and rank and quartile kinases based on the absolute Score values
#'
#' @param df dataframe, UKA table output (requires at least Kinase and Z columns)
#'
#' @return dataframe, Ranked and quartiled UKA table
#'

read_uka <- function(df) {


  if(!is.data.frame(df)) {stop("Make sure your input is a dataframe")}
  if(!all(colnames(df) %in% c("Kinase", "Score"))) {
    stop("check columns names, they should be: Kinase, Score")
  }

  if(!is.character(df$Kinase)) {stop("check that Kinase column is as character column")}
  if(!is.numeric(df$Score)) {stop("check that Score column is as numeric column")}

  rank_kinases(df, trns = "raw", sort = "desc", tool = "UKA")



}
