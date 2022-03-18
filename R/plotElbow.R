
#' @title plotElbow
#' @description plotElbow: This function generates an Elbow plot of the scores for each # of features
#' @param results the final results table from computing all HSS combinations
#' @param elbow elbow value outputted from getElbow during HSS. Use to color the automatically computed elbow point.
#' @noRd
plotElbow <- function(results, elbow = NULL){
     dfElbow =split(results, results$no.markers)  %>% lapply(., function(x) { x[which.max(x$score),]})  %>% do.call(rbind, .)
     dfElbow$color = ifelse(dfElbow$no.markers == elbow, "eblow", "other")
     if (is.null(elbow)) {
          p.Elbow = ggplot2::ggplot(dfElbow, ggplot2::aes(x=no.markers, y = score)) +
               theme(legend.position = 'none') +
               ggplot2::geom_line() +
               ggplot2::geom_point()

          return(p.Elbow)
     }
     p.Elbow = ggplot2::ggplot(dfElbow, ggplot2::aes(x=no.markers, y = score)) +
          theme(legend.position = 'none') +
          ggplot2::geom_line() +
          ggplot2::geom_point(aes(color = color))

     return(p.Elbow)
}

