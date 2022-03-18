
#################################
#################################
#################################
## SEPARATION METRIC FUNCTIONS ##
#################################
#################################
#################################

########################
## EUCLIDEAN DISTANCE ##
########################
#' @description getScore_euclidean: An aggregate function to compute LDA and euclidean distance score for HSS.
#' @param x dataframe of training data
#' @param y vector of class labels matching training data rows
#' @param cols vector of column names
#' @noRd
getScore_euclidean <- function(x, y, cols) {
     # performs LDA using columns provided and returns lowest euclidean distance between pop means
     lda.out <- lda(y~as.matrix(x[, cols]))
     eucl_score = min(dist(lda.out$means %*% lda.out$scaling[,1:2]))
     # message(eucl_score)
     return(eucl_score)
}

#' @description getScore_euclidean_mean: An aggregate function to compute LDA and an average euclidean distance score across all classes for HSS
#' @param x dataframe of training data
#' @param y vector of class labels matching training data rows
#' @param cols vector of column names
#' @noRd
getScore_euclidean_mean <- function(x, y, cols) {
     # performs LDA using columns provided and returns lowest euclidean distance between pop means
     lda.out <- lda(y~as.matrix(x[, cols]))
     eucl_score = mean(dist(lda.out$means %*% lda.out$scaling[,1:2]))
     # message(eucl_score)
     return(eucl_score)
}

#' @description getScore_euclidean_2class: An aggregate function to compute LDA and euclidean distance score for HSS specifically for 2-class LDA
#' @param x dataframe of training data
#' @param y vector of class labels matching training data rows
#' @param cols vector of column names
#' @noRd
getScore_euclidean_2class <- function(x, y, cols) {
     # performs LDA using columns provided and returns lowest euclidean distance between pop means
     lda.out <- lda(y~as.matrix(x[, cols]))
     eucl_score = min(dist(lda.out$means %*% lda.out$scaling[,1]))
     # message(eucl_score)
     return(eucl_score)
}
#######################
## SILHOUETTE SCORE  ##
#######################
#' @description getScore_silhouette: An aggregate function to compute LDA and silouette score for HSS.
#' @param x dataframe of training data
#' @param y vector of class labels matching training data rows
#' @param cols vector of column names
#' @noRd
#' @keywords internal
getScore_silhouette  <- function(x, y, cols) {

     ## computes precise silhouette score but is slow since it requires a distance matrix
     df = x[, cols] %>% as.matrix()
     lda.out <- lda(y~as.matrix(x[, cols]))

     dist.matrix = df %>% dist() %>% as.matrix()
     silhoutte.output = silhouette(x = as.numeric(factor(y)), dmatrix = dist.matrix, do.clus.stat = TRUE, do.n.k=TRUE )

     silh.result.all = df %>% cbind(., data.frame(cluster = silhoutte.output[,1], neighbor = silhoutte.output[,2], sil_width = silhoutte.output[,3] ))
     sum.sil = summary(silhoutte.output)
     silhouette.score = mean(sum.sil$clus.avg.widths)
     # message(silhouette.score)
     return(silhouette.score)
}


#' @title create_pixel_grid
#' @description create_pixel_grid: This function generates a pixel grid template, defaults to 10,000 pixels
#' @param xbreaks the # of pixels to break the x-axis into. Defaults to 100.
#' @param ybreaks the # of pixels to break the y-axis into. Defaults to 100.
#' @keywords internal
create_pixel_grid <- function(xbreaks =100, ybreaks = 100) {
     xbreaks <-xbreaks
     ybreaks <-ybreaks
     pixel.grid = expand.grid(1:xbreaks,1:ybreaks) %>% data.frame(.) %>%
          rename(x=Var1, y = Var2) %>%
          mutate(pixel = paste(x, y,sep= "."))
     return(pixel.grid)
}

#' @description: generate_density_map: This function computes the proportion (density) of the class labels in each pixel of the pixel grid generated using create_pixel_grid.
#' @param data: a dataframe with 3 columns: x-axis coordinates (labeled x), y-axis coordinates(labeled y), and the class labels (labeled as labels)
#' @param pixel.grid: the output from the create_pixel_grid function.
#' @param xbreaks: the # of pixels to break the x-axis into. Defaults to 100.
#' @param ybreaks: the # of pixels to break the y-axis into. Defaults to 100.
#' @noRd
generate_density_map <- function(data, pixel.grid = pixel.grid, xbreaks = 100, ybreaks = 100) {
     # data_pixel_counts <- lapply(unique(data$labels), function(class.label) {
     #   message(class.label)

     xbin <- cut(data$x, xbreaks, include.lowest = TRUE)
     ybin <- cut(data$y, ybreaks, include.lowest = TRUE)

     data_pixel_counts <- data %>%
          ungroup() %>%
          mutate(xout = as.numeric(xbin),
                 yout = as.numeric(ybin),
                 pixel = paste(xout,yout,sep=".")) %>%
          group_by(pixel,labels) %>%
          summarize(count = n())  %>%
          ungroup() %>%
          right_join(.,pixel.grid,by="pixel") %>%
          mutate(count.0 = ifelse(is.na(count), 0, count)) %>%
          ungroup() %>%
          group_by(labels) %>%
          mutate(percent.0 = count.0 / sum(count.0))
     return(data_pixel_counts)
}


#' @title calculate_Score
#' @description calculate_Score: This function computes the pixel class entropy score.
#' @param density_metric_output the output from the function, density_metric_output
#' @noRd
calculate_pceScore <- function(data = density_metric_output) {
     data <- na.omit(data) %>% ungroup()
     pce.score <- data %>%
          group_by(pixel) %>%
          summarize(entropy = entropy(count),
                    num.of.labels = n()) %>%
          ungroup() %>%
          mutate(pce.score = 1-(entropy/log2(num.of.labels)),
                 pce.score = case_when(is.na(pce.score) ~ 1,
                                       TRUE ~ pce.score)) # %>%  select(pixel, num.of.labels,entropy, pce.score)
     return(pce.score)
}

#' @title computePCEscore
#' @description computePCEscore: computes the pixel class entropy score for any biaxial dataset with class labels
#' @param data a dataframe with 3 columns: x-axis coordinates (labeled `x`), y-axis coordinates(labeled `y`), and the class labels (labeled as `labels`)
#' @noRd
computePCEscore <- function(data) {
     ## data in format of x(axis1), y(axis2), class label of interest

     if (!exists("pixel.grid")){
          pixel.grid <<- create_pixel_grid()
     }

     density_metric_output <- generate_density_map(data = data, pixel.grid = pixel.grid)
     pce.score_output <- calculate_pceScore(data = density_metric_output)
     pce.score <-  mean(pce.score_output$pce.score)
     return(pce.score)
}

##########################################
## PIXEL CLONALITY ENTROPY (PCE) SCORE  ##
##########################################
#' @description getScore_pce: An aggregate function to compute LDA and the PCE score for HSS.
#' @param x dataframe of training data
#' @param y vector of class labels matching training data rows
#' @param cols vector of column names
#' @noRd
getScore_pce <- function(x, y, cols) {
     ## pixel clonality scoring method
     lda.out <- lda(y~as.matrix(x[, cols]))

     if (!exists("pixel.grid")){
          pixel.grid <- create_pixel_grid()
     }
     data.pixels <- as.matrix(x[, cols]) %*% lda.out$scaling
     data.pixels <- data.pixels[ , c("LD1","LD2" )]%>% data.frame()
     data.pixels["labels"] <- y
     colnames(data.pixels) <- c("x","y","labels")

     pce.score = computePCEscore(data = data.pixels)
     # message(pce.score)
     return(pce.score)
}
###############################
## CUSTOM TEMPLATE FUNCTION  ##
###############################
#' @title getScore_custom
#' @description getScore_custom: placeholder for a custom metric
#' @param x dataframe of training data
#' @param y vector of class labels matching training data rows
#' @param cols vector of column names
#' @param custom.score.method a custom-function that takes in x, y, and cols, and outputs a score where the larger the value, the more optimal your separation criteria is.
getScore_custom <- function(x, y, cols, custom.score.method, ...) {
     df = x[, cols, with=F]
     lda.out <- lda(y~., data=df)
     custom_output = custom.score.method(x, y, cols, lda.out, ...)
     return(custom_output)
}

#####################################
## aggregate function for getScore ##
#####################################
#' @description getScore: wrapper function for each getScore metric
#' @param x dataframe of training data
#' @param y vector of class labels matching training data rows
#' @param cols vector of column names
#' @param score.method the scoring method to use.
#' @noRd
getScore <- function(x , y, cols, score.method) {
     if (length(cols) > 1) {
          if (score.method == "euclidean") {
               # message("euclidean")
               scoreFunction <- getScore_euclidean(x, y, cols)
               return(scoreFunction)
          } else if (score.method == "euclidean_2class") {
               # message("euclidean")
               scoreFunction <- getScore_euclidean_2class(x, y, cols)
               return(scoreFunction)
          } else if (score.method == "euclidean_mean") {
               # message("euclidean")
               scoreFunction <- getScore_euclidean_mean(x, y, cols)
               return(scoreFunction)
          } else if (score.method == "silhouette") {
               # message("silhouette")
               ## pixel clonality scoring method
               scoreFunction <- getScore_silhouette(x, y, cols)
               return(scoreFunction)
          } else if (score.method == "pixel.entropy") {
               # message("pixel.entropy")
               ## pixel clonality scoring method
               scoreFunction <- getScore_pce(x, y, cols)
               return(scoreFunction)
          } else if (score.method == "pixel.density") {
               # message("pixel.density")
               ## pixel clonality scoring method
               scoreFunction <- getScore_pixelDensity(x, y, cols)
               return(scoreFunction)
          } else if (score.method == "custom") {

               if (is.null(custom.score.method)) {
                    stop("Must provide a custom score method")
               }
               # message("custom")
               scoreFunction <- getScore_custom(cols, x, y, custom.score.method)
               return(scoreFunction)
          }
     }
}
