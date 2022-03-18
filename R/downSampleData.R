

#' @title downSampleData
#' @description Downsamples data to one of 3 conditions: 5000 cells in each class label or the class label with the minimum # of cells
#' @param x table with predictors of interest
#' @param y vector of class labels
#' @noRd
downSampleData <- function(x, y){

     df = cbind(x,labels = y)
     rownames(df) = paste0("inputID", rownames(df))
     inputIDs = rownames(df)
     df[["inputIDs"]] = inputIDs

     min.value = min(table(y))
     if (min.value >= 5000) {
          message("Downsampling data to 5000 cells in each class label")
          x = df %>% group_by(labels) %>% sample_n(5000)
     } else  if (min.value <= 5000) {
          message("Downsampling data to ", min.value," cells in each class label, which is the class label with the minimum # of cells")
          print(table(train.y))

          x = df %>% group_by(labels) %>% sample_n(min(table(y)))
     }
     downsampledCells = df
     downsampledCells$inputIDs = downsampledCells$inputIDs %in% x$inputIDs
     message("Downsampling data to ", min.value," cells in each class label, which is the class label with the minimum # of cells")
     return(downsampledCells) ## returns a subsetted dataset with labels relative to original cells that were INPUT
}

