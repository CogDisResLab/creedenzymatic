#' Rnak Kinases based on a score
#'
#' @param df, dataframe
#'
#' @return
#' @export
#'
#' @examples
#'
#'
rank_kinases <- function(df, trns = c("raw", "abs"), sort = c("desc", "asc"), tool = c("KRSA", "UKA")) {


  perc <- ecdf(1:length(df$Score))

  df %>%
    mutate(Rank = dense_rank(
      if(sort == "desc"){
      desc( if(trns == "abs"){abs(Score)}
            else {Score})
      }
      else{
        if(trns == "abs"){abs(Score)}
              else {Score}
      }
      )) %>%
    mutate(Method = tool) %>%
    mutate(Qnt = perc(Rank),
           Qrt = case_when(
             Qnt <= 0.25 ~ 1,
             Qnt <= 0.5 ~ 2,
             Qnt <= 0.75 ~ 3,
             Qnt <= 1 ~ 4,
           ))


}
