#' Rank Kinases based on a score
#'
#' This function will scale the scores on a percentile and quartile scales
#'
#'
#'
#' @param df, dataframe with 2 columns: Kinase, Score
#' @param trns, for transformation of the score, the values accepted for this argument are abs and raw (abs: use absolute values of scores, raw: no transformation)
#' @param sort, accepts either asc or desc (ascending and descending)
#' @param tool, specifying the name of the tool
#'
#' @return
#' @export
#'
#' @examples
#'
#'
rank_kinases <- function(df, trns = c("raw", "abs"), sort = c("desc", "asc"), tool = c("KRSA", "UKA")) {

  if(!all(trns %in% c("raw", "abs"))) {
    stop("the trns argument must be either raw or abs")
  }

  if(!all(sort %in% c("desc", "asc"))) {
    stop("the sort argument must be either desc or asc")
  }

  perc <- stats::ecdf(1:length(df$Score))

  df %>%
    dplyr::mutate(Rank = dplyr::dense_rank(
      if(sort == "desc"){
      desc( if(trns == "abs"){abs(Score)}
            else {Score})
      }
      else{
        if(trns == "abs"){abs(Score)}
              else {Score}
      }
      )) %>%
    dplyr::mutate(Method = tool) %>%
    dplyr::mutate(#Qnt = ntile(abs(Score), 4),
           Perc = if(sort == "desc"){
             if(trns == "abs"){dplyr::percent_rank(abs(Score))}
                   else {dplyr::percent_rank(Score)}
           }
           else{
             if(trns == "abs"){dplyr::percent_rank(desc(abs(Score)))}
             else {dplyr::percent_rank(Score)}
           }
             ,
           Qrt = dplyr::case_when(
             1 <= Perc | Perc >= 0.75  ~ 4,
             0.75 < Perc | Perc >= 0.5  ~ 3,
             0.5 < Perc | Perc >= 0.25  ~ 2,
             0.25 < Perc | Perc >= 0  ~ 1
           ))


}
