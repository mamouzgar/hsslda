##' Authors: Meelad Amouzgar and David Glass
##' @author Meelad Amouzgar
##' @author David Glass
##' date: July 5th, 2021
##' Description: defines several functions required for hybrid subset selection:
##' 1. hsslda: the hybrid subset selection (HSS)
##' 2.downsampleBalance: downsampling function for faster HSS using the separation metric of choice with balanced downsampling of input class labels
##' 3. various separation metric functions:
##' 3a. Euclidean
##' 3b. Silhouette score
##' 3c. Pixel entropy
##' 3d. Pixel density
##'
##'
##'
##'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 geom_point
#'
#' @importFrom dplyr filter
#' @importFrom dplyr select
#' @importFrom dplyr group_by
#' @importFrom dplyr ungroup
#' @importFrom dplyr rename
#' @importFrom dplyr summarize
#' @importFrom dplyr n
#' @importFrom dplyr %>%
#' @importFrom dplyr mutate
#' @importFrom dplyr rowwise
#' @importFrom dplyr bind_rows
#' @importFrom dplyr left_join
#' @importFrom dplyr right_join
#' @importFrom dplyr full_join
#' @importFrom dplyr case_when
#' @importFrom dplyr sample_n
#' @importFrom tidyr spread
#' @importFrom tidyr gather
#'
#' @importFrom MASS lda
#' @importFrom stats dist
#' @importFrom entropy entropy
#' @importFrom cluster silhouette
#'





#' @title plotElbow
#' @description plotElbow: This function generates an Elbow plot of the scores for each # of features
#' @param results the final results table from computing all HSS combinations
#' @param elbow elbow value outputted from getElbow during HSS. Use to color the automatically computed elbow point.
#' @noRd
plotElbow <- function(results, elbow = NULL){
        dfElbow =split(results, results$no.markers)  %>% lapply(., function(x) { x[which.max(x$score),]})  %>% do.call(rbind, .)

        if (is.null(elbow)) {
                p.Elbow = ggplot2::ggplot(dfElbow, ggplot2::aes(x=no.markers, y = score)) +
                        ggplot2::geom_line() +
                        ggplot2::geom_point()

                return(p.Elbow)
        }
        p.Elbow = ggplot2::ggplot(dfElbow, ggplot2::aes(x=no.markers, y = score)) +
                ggplot2::geom_line() +
                ggplot2::geom_point(color = ifelse(dfElbow$no.markers == elbow, "red", "black"))

        return(p.Elbow)
}





#################################
#################################
#################################
## SEPARATION METRIC FUNCTIONS ##
#################################
#################################
#################################

##############################
## PIXEL ANALYSIS FUNCTIONS ##
##############################
# create_pixel_grid <- function(xbreaks =100, ybreaks = 100) {
#   pixel.grid = base::expand.grid(1:xbreaks,1:ybreaks) %>% data.frame(.)
#   base::colnames(pixel.grid) = c("x","y")
#   pixel.grid$pixel = base::paste(pixel.grid$x,pixel.grid$y, sep =".")
#   return(pixel.grid)
# }
#
# generate_density_map <- function(data, pixel.grid = pixel.grid, xbreaks = 100, ybreaks = 100) {
#   ## returns a density map only of pixels that have cells in them
#   xbin <- base::cut(data$x, xbreaks, include.lowest = TRUE)
#   ybin <- base::cut(data$y, ybreaks, include.lowest = TRUE)
#
#   data_pixel_counts = cbind(data, xout = as.numeric(xbin), yout = as.numeric(ybin) )
#   data_pixel_counts["pixel"] = paste(data_pixel_counts$xout, data_pixel_counts$yout, sep=".")
#   data_pixel_counts = table(pixel = data_pixel_counts$pixel, labels = data_pixel_counts$labels) %>% as.data.frame() %>% .[.$Freq !=0, ] # %>% base::merge(.,pixel.grid, by = "pixel", all.y = TRUE)
#   data_pixel_counts["count.0"] = ifelse(is.na(data_pixel_counts$Freq), 0,data_pixel_counts$Freq )
#   data_pixel_counts = data_pixel_counts %>% split(., .$labels) %>% lapply(., function(dp){dp["percent.0"] = dp$count.0/sum(dp$count.0) ;return(dp) }) %>% do.call("rbind", .)
#   return(data_pixel_counts)
# }


#' @title create_pixel_grid
#' @description create_pixel_grid: This function generates a pixel grid template, defaults to 10,000 pixels
#' @param xbreaks the # of pixels to break the x-axis into. Defaults to 100.
#' @param ybreaks the # of pixels to break the y-axis into. Defaults to 100.
#' @keywords internal
#' @export
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


###############################
## PIXEL CLASS ENTROPY SCORE ##
###############################
# calculate_pceScore <- function(data) {
#   pce.score = split(data, data$pixel) %>% lapply(., function(dp) {data.frame(entropy = entropy(dp$count), num.of.labels = length(dp$labels)) } ) %>% do.call("rbind",.)
#   pce.score["pce.score"] = 1-(pce.score$entropy/log2(pce.score$num.of.labels))
#   pce.score["pce.score"] = ifelse(is.na(pce.score$pce.score), 1 , pce.score$pce.score)
#   pce.score["pixel"] = rownames(pce.score)
#   return(pce.score)
# }

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
#' @export
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

##########################
## PIXEL DENSITY SCORE  ##
##########################
# description: calculate_pixelDensityScore: computes the pixel density score for any biaxial dataset with class labels. Must finish recoding. FYI see lapply commented out for reminder on what needs to change
# param data: a dataframe with 3 columns: x-axis coordinates (labeled `x`), y-axis coordinates(labeled `y`), and the class labels (labeled as `labels`)
# calculate_pixelDensityScore <- function(data = density_metric_output) {
#
#   data <- na.omit(data) %>% ungroup()
#   # message("calculate-pixel-clonality")
#
#   pixelDensity.score <- data %>%
#     dplyr::ungroup() %>%
#     dplyr::mutate(binary.labels = ifelse(labels == class.label, "class.of.interest", "other")) %>%
#     dplyr::select(pixel, x, y, percent.0, binary.labels) %>%
#     dplyr::group_by(pixel, x, y, binary.labels) %>%
#     dplyr::summarize(count = sum(count.0),
#               percent = sum(percent.0)) %>%
#     tidyr::gather(key = "approach", value = "quantity", -pixel,-x,-y,-binary.labels) %>%
#     tidyr::spread(key = "binary.labels", value = "quantity") %>%
#     dplyr::rowwise() %>%
#     dplyr::mutate(class.of.interest = ifelse(is.na(class.of.interest), 0, class.of.interest),
#            other = ifelse(is.na(other), 0, other),
#            density.metric = (class.of.interest / (other+class.of.interest)),
#            # density.metric = ifelse()
#            labels = class.label)  %>%
#     # dplyr::select(pixel,x,y, approach, class.of.interest, other,density.metric, labels) %>%
#     dplyr:: group_by(approach, labels) %>%
#     dplyr::summarize(density.summary = sum(density.metric)) %>%
#     dplyr::mutate(density.summary.normalized = density.summary / 10000)
#   return(pixelDensity.score)
# }


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

##########################################
## PIXEL CLONALITY ENTROPY (PCE) SCORE  ##
##########################################
#' @description getScore_pce: An aggregate function to compute LDA and the PCE score for HSS.
#' @param x dataframe of training data
#' @param y vector of class labels matching training data rows
#' @param cols vector of column names
#' @noRd
#' @keywords internal
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

# getScore_pixelDensity <- function(x, y, cols) {
#   ## pixel clonality scoring method
#   # performs LDA using columns provided and returns lowest euclidean distance between pop means
#   lda.out <- lda(y~., data=x[, cols])
#
#   if (!exists("pixel.grid")){
#     pixel.grid <- create_pixel_grid()
#   }
#   data.pixels <- as.matrix(x[, cols]) %*% lda.out$scaling
#   data.pixels <- data.pixels[ , c("LD1","LD2" )]%>% data.frame()
#   data.pixels["labels"] <- y
#   colnames(data.pixels) <- c("x","y","labels")
#
#   density_metric_output <- generate_density_map(data = data.pixels, pixel.grid = pixel.grid)
#   pixelDensity_output <- calculate_pixelDensityScore(data = density_metric_output)
#   pixelDensity.score <-  mean(pixelDensity_output$density.summary.normalized)
#   # message(pixelDensity.score)
#   return(pixelDensity.score)
# }

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

#' @description hybridSubsetSelection: function that performs hybrid stepwise subset selection
#' @param x dataframe of training data
#' @param y vector of class labels matching training data rows
#' @param score.method the scoring method to use.
#' @param custom.score.method optional input, a custome scoring function (see getScore_custom)
#' @noRd
hybridSubsetSelection <- function(x, y, score.method , custom.score.method = NULL) {
        options(dplyr.summarise.inform = FALSE)

        # performs hybrid stepwise subset selection on LDA reduced dimensions
        # Inputs:
        #   x - data.table to evaluate, rows are cells, columns are columns to evaluate
        #   y - vector of observation classes
        #   two.d - logical if true creates two new axes, if false only one
        # Outputs:
        #   matrix of coefficients where rows are markers and columns are axes

        ### global data structures ###
        keep <- NULL
        channels <- colnames(x)
        n.channels <- length(channels)
        current.score <- 0
        continue <- TRUE
        results <- setNames(data.frame(matrix(nrow=1, ncol=n.channels)), channels)
        results[1,] <- as.list(rep(F, n.channels))
        subtract.log <- results[0,] # record of keep values inputted into subtractOne
        results$score <- 0

        hss.result = list() ## final output

        #############################
        #############################
        #############################
        ##### SCORING FUNCTIONS #####
        #############################
        #############################
        #############################



        # #######################
        # ## SILHOUETTE SCORE  ##
        # #######################
        # getScore_silhouette  <- function(x, y, cols) {
        #   df = x[, cols, with=F]
        #   lda.out <- lda(y~., data=df)
        #   dist.matrix <- df %>% dist() %>% as.matrix()
        #   silhoutte.output <- cluster::silhouette(x = as.numeric(factor(y)), dmatrix = dist.matrix, do.clus.stat = TRUE, do.n.k=TRUE )
        #
        #   silh.result.all<- df %>%
        #     bind_cols(., data.frame(cluster = silhoutte.output[,1], neighbor = silhoutte.output[,2], sil_width = silhoutte.output[,3] ))
        #   sum.sil <- summary(silhoutte.output)
        #   silhouette.score = mean(sum.sil$clus.avg.widths)
        #   message(silhouette.score)
        #   return(silhouette.score)
        # }

        ######################################
        ## PIXEL CLASS ENTROPY (PCE) SCORE  ##
        ######################################
        # getScore_pce <- function(x, y, cols) {
        #   ## pixel clonality scoring method
        #   # performs LDA using columns provided and returns lowest euclidean distance between pop means
        #   lda.out <- lda(y~., data=x[, cols, with=F])
        #   data.pixels <- makeAxes(dt = x[, cols, with=F], co=lda.out$scaling)
        #   data.pixels <- data.pixels[ , c("ld1","ld2" )]
        #   data.pixels$labels <- y
        #   colnames(data.pixels) <- c("x","y","labels")
        #
        #   if (!exists("pixel.grid")){
        #     pixel.grid <<- create_pixel_grid()
        #   }
        #
        #   density_metric_output <- generate_density_map(data = data.pixels, pixel.grid = pixel.grid)
        #   pce.score_output <- calculate_pceScore(data = density_metric_output)
        #   pce.score <-  mean(pce.score_output$pce.score)
        #   message(pce.score)
        #   return(pce.score)
        # }







        ####################
        ## begin analysis ##
        ####################

        ## original, euclidean distance based function
        # } else if (score.method == "pixel.density") {
        #   ## pixel density scoring method
        #   getScore <<- function(cols) {
        #     # performs LDA using columns provided and returns lowest euclidean distance between pop means
        #     message(length(cols))
        #     lda.out <<- lda(y~., data=x[, cols, with=F])
        #     # message(lda.out)
        #     df.density <<- makeAxes(dt = x[, cols, with=F], co=lda.out$scaling)
        #     df.density <<- df.density[ , c("ld1","ld2" )]
        #     df.density$labels <- factor(dat$labels)
        #     colnames(df.density) <- c("x","y","labels")
        #     coordinate.density <- list()
        #     pixel.label.density_sc <- list()
        #
        #     pixel.grid <- create_pixel_grid()
        #     density_metric_output <- generate_density_map(data = df.density, pixel.grid = pixel.grid)
        #     # pixel.clonality_output <- calculate_entropy_metric(data = density_metric_output) %>% mutate(filename = file.name)
        #     # pixel.clonality.metric <<- dplyr::bind_rows(pixel.clonality.metric, pixel.clonality_output)
        #     coordinate.density <- dplyr::bind_rows(coordinate.density, density_metric_output )
        #     density_metric_output <- calculate_pixel_density_metric(data=density_metric_output)
        #     pixel.label.density_sc <-  dplyr::bind_rows(pixel.label.density_sc, density_metric_output)
        #     density_metric_output <- aggregate_pixel_density_metric(data=density_metric_output)
        #     # density_metric_output$filename = file.name
        #     ave.density_metric_output <- density_metric_output %>% dplyr::filter(approach == "percent")
        #     mean.score = mean(ave.density_metric_output$density.summary.normalized)
        #
        #     message(mean.score)
        #     if (two.d) return(mean.score)
        #     return(min(dist(lda.out$means %*% lda.out$scaling[,1])))
        #   }

        ##subset functions
        addOne <- function() {
                # Evaluates the addition of each channel not in keep to keep. Adds best and updates current.score
                temp.results <- results[0,]
                # message(keep)
                # message(channels)
                for (channel in channels[!channels %in% keep]) {
                        temp.keep <- c(keep, channel)
                        temp.score <- getScore(x, y, cols = temp.keep, score.method)
                        temp.results <- rbind(temp.results, as.list(channels %in% temp.keep) %>% append(temp.score))
                        # message(temp.results)
                        # temp.results <<-temp.results
                }
                colnames(temp.results) = colnames(results)
                current.score <<- max(temp.results$score)
                new.keep <- temp.results[temp.results$score == current.score, channels]
                if (nrow(new.keep) > 1) new.keep <- new.keep[sample(.N,1)]
                keep <<- channels[as.logical(new.keep)]
                results <<- unique(rbind(results, temp.results))
        }

        subtractOne <- function() {
                # Evaluates the subtraction of each channel from keep. Removes worst if it improves score and updates current.score
                # If a better subset is found, it calls itself.
                # If this keep has been evaluted before, exits
                # message("test")
                subtract.log <<- rbind(subtract.log, as.list(channels %in% keep))
                if (anyDuplicated(subtract.log) > 0) {
                        subtract.log <<- unique(subtract.log)
                        return()
                }
                temp.results <- results[0,]
                # temp.results <<- temp.results
                for (channel in keep) {
                        temp.keep <- keep[!keep %in% channel]
                        temp.score <- getScore(x, y, cols = temp.keep, score.method)
                        if (is.null(temp.score)){
                             temp.score = 0
                        }
                        temp.results <- rbind(temp.results, as.list(channels %in% temp.keep) %>% append(temp.score))
                }
                # message(colnames(temp.results))
                # message( colnames(results))
                # print(colnames(results))
                colnames(temp.results) = colnames(results)
                # temp.results <<- temp.results
                # current.score <<- current.score

                if (max(temp.results$score) > current.score) {
                        # current.score <<- base::max(temp.results$score)
                        current.score <- max(temp.results$score)
                        new.keep <- temp.results[temp.results$score == current.score, channels]
                        if (nrow(new.keep) > 1) new.keep <- new.keep[1, ]
                        keep <<- channels[as.logical(new.keep)]
                        results <<- unique(rbind(results, temp.results))
                        subtractOne()
                } else results <<- unique(rbind(results, temp.results))
        }


        #############################
        #############################
        #############################
        initializeKeep <- function() {
                # chooses the best scoring pair of markers to initialize keep
                temp.results <- results[0,]

                myChannels = expand.grid(channel.1 = channels, channel.2 = channels) %>% .[.$channel.1 != .$channel.2, ]
                # message(myChannels)
                temp.results = apply(myChannels, 1, function(row.temp.keep){
                        temp.keep = row.temp.keep %>% unlist()
                        temp.score = getScore(x, y, cols = temp.keep, score.method)
                        temp.result = as.list(channels %in% temp.keep) %>% append(temp.score)
                        return(temp.result)
                })
                # temp.results <<- temp.results
                temp.results <- do.call("rbind", lapply(temp.results, unlist)) %>% data.frame()
                colnames(temp.results) = colnames(results)

                current.score <<- max(temp.results$score)
                new.keep <- temp.results[temp.results$score==current.score, channels]
                old.keep <- new.keep
                if (nrow(new.keep) > 1)  new.keep <- new.keep[ which.max(apply(new.keep, 1, sum)), ]
                keep <<- channels[as.logical(new.keep)]
                results <<- unique(rbind(results, temp.results))
        }

        getElbow <- function(res) {
                # takes results and returns the elbow point
                res.lite <- res[ , "no.markers"] %>% unique() %>% .[-1 ] %>% data.frame(no.markers = .)
                res.lite[, "V1"] <- lapply(split(res, res$no.markers) , function(df) {max(df$score)}) %>% unlist(.) %>% .[-1]
                res.lite<-res.lite
                slope <- (res.lite$V1[nrow(res.lite)] - res.lite$V1[1]) / (res.lite$no.markers[nrow(res.lite)] - res.lite$no.markers[1])
                intercept <- res.lite$V1[1] - slope * res.lite$no.markers[1]
                perp.slope <- -1 / slope
                perp.int <- res.lite$V1 - (perp.slope * res.lite$no.markers)
                xcross <- (intercept - perp.int) / (perp.slope - slope)
                ycross <- slope * xcross + intercept
                dists <-  sqrt((res.lite$no.markers - xcross)^2 + (res.lite$V1 - ycross)^2)
                elbowi <- max(which(dists==max(dists))) # if dists are tie, take the largest number of channels
                return(elbowi+1)
        }

        ### main ###
        initializeKeep()
        while(continue) {
                message(paste("Number of markers:", length(keep)))
                addOne()
                message(paste("Number of markers:", length(keep)))
                if (length(keep) > 3) subtractOne()
                if (length(keep)==length(channels)) continue <- FALSE
        }
        results["no.markers"] = apply(results[, channels], 1, sum)
        elbow <- getElbow(res=results)

        markers <- results[results$no.markers==elbow, ] %>%
                .[.$score==max(.$score), colnames(.) %in% channels] %>%
                unlist() %>%
                .[.==1] %>%
                names(.)
        # message(markers)
        lda.out <- lda(y~., data=x[, markers])


        ## save lda.out results
        hss.result[["method"]] = score.method
        hss.result[["finalMarkers"]] = markers
        hss.result[["HSSscores"]] = results
        hss.result[["ElbowPlot"]] = plotElbow(results = results, elbow = elbow)
        hss.result[["HSS-LDA-model"]] = lda.out


        ## restore options
        options(dplyr.summarise.inform = TRUE) ## turn it back on
        # pixel.grid <- NULL

        return(hss.result)
}

#' @title makeAxes
#' @description makeAxes: function that generates new HSS-LDA axes given an LDA model coefficients
#' @param df a dataframe of cells (rows) by markers/genes (columns)
#' @param co the dataframe of LDA coefficients from MASS::lda.
#' @keywords internal
#' @export
makeAxes <- function(df, co=coefficients) {
        # makes new axes based on coefficients
        # Inputs:
        #   df - data.frame of data
        #   co - matrix of coefficients
        #   axis.name - character vector of new axis name (e.g. "ld" results in "ld1" and "ld2")
        # Outputs:
        #   df - data.frame
        x <- as.matrix(df[, rownames(co)])
        HSS.LDs = x %*% co
        colnames(HSS.LDs) = paste0("HSS_", colnames(HSS.LDs), sep = "")
        df = cbind(df, HSS.LDs)
        return(df)
}




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



#' @title runHSS
#' @description This function runs hybrid subset selection.
#' @param x table with predictors of interest
#' @param y vector of class labels
#' @param score.method scoring metric for feature selection using HSS. Options include: 'euclidean', 'silhouette', 'pixel.density', 'pixel.entropy', or 'custom'.
#' @param custom.score.method function for your custom scoring metric. Score.method must be 'custom'
#' @param downsample boolean indicating whether to downsample data. Defaults to TRUE. Downsamples to 5000 cells per class label, or balanced sampling based on the class label with minimum # of cells.
#' @keywords internal
#' @export
runHSS <- function(x, y, score.method, custom.score.method = NULL, downsample = TRUE){
        if (!score.method %in% c("euclidean","euclidean_2class", "euclidean_mean", "silhouette", "pixel.density","pixel.entropy")){
                stop("score.method method must be: 'euclidean', 'euclidean_2class', 'euclidean_mean', silhouette', 'pixel.density', 'pixel.entropy', or 'custom'.")
        }

        ## save original input data
        orig.data = x

        if (downsample == TRUE) {
                downsampledCells = downSampleData(x = x, y = y)
                downsampledCells[["inputIDs"]] = NULL
                cells.included.in.downsampling.analysis = downsampledCells$inputIDs ## a boolean vector of whether the input cells were included in the downsampling analysis
                x = downsampledCells[colnames(x)]
        } else {
                message("Using all input cells for HSS-LDA...")
                cells.included.in.downsampling.analysis = rep(TRUE, nrow(orig.data))
        }

        message("Optimizing dimensionality reduction using ", score.method, " metric...")

        ## HSS results
        hss.results <- hybridSubsetSelection(x, y, score.method = score.method, custom.score.method = custom.score.method)
        coefficients = hss.results[["HSS-LDA-model"]]$scaling
        dat <- makeAxes(df=orig.data, co=coefficients)
        dat[["labels"]] = y
        hss.results[["HSS-LDA-result"]] = dat

        ## add boolean vector of downsampled cell labels
        hss.results[["downsampled.cells.included.in.analysis"]] = cells.included.in.downsampling.analysis

        return(hss.results)
}

# setwd("~/phd-projects")
# library(dplyr)
# library(magrittr)
# library(ggplot2)
# path  <- "~/phd-projects/sc-lda/data/analysis-ready/metabolism-CD8-naive-data_cell-allmarkers.csv"
# channels <- c('GLUT1', 'HK2', 'GAPDH', 'LDHA', 'MCT1', 'PFKFB4', 'IDH2', 'CyclinB1',
#               'GLUD12', 'CS', 'OGDH', 'CytC', 'ATP5A', 'S6_p', 'HIF1A', 'PDK1_p', 'NRF1',
#               'NRF2_p', 'XBP1', 'VDAC1', 'OPA1', 'DRP1', 'ASCT2', 'GLS', 'GOT2', 'CPT1A',
#               'ACADM', 'IdU', 'BrU', 'Puromycin', 'H3_p',
#               "CD45RA",
#               # "CD69","CD25"
#               "CD69","CD3","CD98","CD25","CD27","CD137","CD57"
#
#               # 'DNA','barium'
# )
# dat_input <- data.table::fread(path)
# dat <- dat_input %>%
#   filter(H3_p < 0.05,
#          dead < 0.3)  %>%
#   group_by(labels) %>%
#   sample_n(1000) %>%
#   data.frame()
# # fwrite(dat, file = output_path)
# # training data with appropriate channels
# train.x <- dat[, channels[1:12]]
# # training class labels - must be 3+ unique classes
# train.y <- dat$labels
# #
# # # x=train.x
# # # y=train.y
# #
# start.time = Sys.time()
# # hss.results=runHSS(x = train.x, y = train.y, score.method = "euclidean")
# # hss.results=runHSS(x = train.x, y = train.y, score.method = "silhouette")
# hss.results=runHSS(x = train.x, y = train.y, score.method = "pixel.entropy")
# end.time = Sys.time()
# coefficients <- hybridSubsetSelection(x=train.x, y=train.y, score.method = "pixel.entropy")
# #
# end.time-start.time
#
# start.time.1 = Sys.time()
# hss.results=runHSS(x = train.x, y = train.y, score.method = "euclidean")
# end.time.1 = Sys.time()
# end.time.1-start.time.1
#
# start.time.2 = Sys.time()
# hss.results=runHSS(x = train.x, y = train.y, score.method = "silhouette")
# end.time.2 = Sys.time()
# end.time.2-start.time.2
#
# start.time.3 = Sys.time()
# hss.results=runHSS(x = train.x, y = train.y, score.method = "pixel.entropy")
# end.time.3 = Sys.time()
# end.time.3-start.time.3
# end.time.1-start.time.1
# end.time.2-start.time.2
# end.time.3-start.time.3
# start.time.4 = Sys.time()
# hss.results=runHSS(x = train.x, y = train.y, score.method = "pixel.density")
# end.time.4 = Sys.time()
# hss.results=runHSS(x = train.x, y = train.y, score.method = "silhouette")







# remotes::install_github("mamouzgar/hsslda")







# dat <- makeAxes()
# # cluster::silhouette(x = y)
# # dist.matr = dist(train.x) %>% as.matrix()
# # silhoutte.output <- cluster::silhouette(x = as.numeric(factor(train.y)), dmatrix = dist.matr, do.clus.stat = TRUE, do.n.k=TRUE )
# ggplot(hss.results$`HSS-LDA-result`, aes(x = ld1, y= ld2)) +
#   geom_point(aes(color = labels)) +
#   viridis::scale_color_viridis(discrete=TRUE)


# Once you’ve got your documentation completed, you can simply run:
#
#   devtools::document()
# This will generate the load_mat.Rd file in the man folder:
#
# You will get one .Rd file for each function in your R package.
#
# Each time you add new documentation to your R function, you need to run devtools::document() again to re-generate the .Rd files.

















# ## computes estimate of the silhouette score but is faster since it bypasses the distance matrix
# df = x[, cols] %>% as.matrix()
# lda.out <- lda(y~as.matrix(x[, cols]))
#
# si = cluster::silhouette(x = as.numeric(factor(rownames(as.matrix(dist(lda.out$means %*% lda.out$scaling[ ,1:2]))))), dmatrix = as.matrix(dist(lda.out$means %*% lda.out$scaling[ ,1:2])))
#
#
# summary(silhoutte.output)$clus.avg.wi
# as.matrix(dist(lda.out$means %*% lda.out$scaling[ ,1:2]))
#
# intercluster.distance = apply(as.matrix(dist(lda.out$means %*% lda.out$scaling[ ,1:2])), MARGIN = 1, max)
# avg.intercluster.distance = mean(matrix(dist(lda.out$means %*% lda.out$scaling[ ,1:2]))) ## average distance to each class label centroid
#
# lds = as.matrix(hsslda::TcellHartmann2020_sampleData[1:5]) %*% lda.out$scaling[ ,1:2] %>% data.frame()
# lds$y = factor(hsslda::TcellHartmann2020_sampleData$labels)
# euc.dist <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2))
#
# centroids = lda.out$means %*% lda.out$scaling[ ,1:2] %>% data.frame() %>%mutate(y = rownames(.)) %>% as_tibble() %>% rename(centroid.LD1=LD1, centroid.LD2 = LD2)
# lds = lds %>% left_join(., centroids ) %>%
#         rowwise() %>%
#         mutate(distance = sqrt(((centroid.LD2 - LD2)^2) + ((centroid.LD1 - LD1)^2)),
#                intercluster.distances = lapply(centroids$y, function(centroid.label){
#                        myCentroid = centroids[centroids$y == centroid.label, ]
#                        dist.to.centroid = sqrt(((LD2 - myCentroid$centroid.LD2)^2) + ((LD2 - myCentroid$centroid.LD1)^2))
#                        names(dist.to.centroid) = centroid.label
#                        return(dist.to.centroid)
#                } ) %>%  list(),
#                closest.cluster = unlist(intercluster.distances) %>% .[!names(.) %in% y] %>% min())
#
# silhs = lds %>% rowwise() %>% mutate(silh = (closest.cluster - distance)/max(lds$distance,lds$closest.cluster) )
# nearest.cluster = apply(lds, MARGIN = 1, function(values){values['distance']})
#
#
# lds.mean = lds %>%
#         group_by(y) %>%
#         summarize(intra.cluster.dist = mean(distance))
#
# lds.silh = lds %>% left_join(.,lds.mean)
# Silhouette.Score = lds.silh  %>% mutate(silh = (intra.cluster.dist - distance)/ max(intra.cluster.dist, distance))  # (b-a)/max(a,b)
#
#
# avg.intracluster.distance = mean(lds$intra.cluster.dist ) ## average distance to the centroid
#
#
#
#
#
# silhoutte.output = cluster::silhouette(x = as.numeric(factor(hsslda::TcellHartmann2020_sampleData$labels)), dmatrix = as.matrix(dist(as.matrix(hsslda::TcellHartmann2020_sampleData[1:5]) %*% lda.out$scaling[ ,1:2])))
# sum.sil = summary(silhoutte.output)
# sum.sil$clus.avg.widths
# silhouette.score # 0.02538235
# )
#
#
# lda.out = MASS::lda(hsslda::TcellHartmann2020_sampleData[-(ncol(TcellHartmann2020_sampleData))],hsslda::TcellHartmann2020_sampleData$labels)



# channels = c('GLUT1', 'HK2', 'GAPDH', 'LDHA', 'MCT1', 'PFKFB4', 'IDH2', 'CyclinB1', 'GLUD12', 'CS', 'OGDH', 'CytC', 'ATP5A', 'S6_p', 'HIF1A')
# train.x = TcellHartmann2020_sampleData %>% dplyr::filter(labels %in% c("day0","day4")) %>% .[channels]
# train.y = TcellHartmann2020_sampleData%>% dplyr::filter(labels %in% c("day0","day4")) %>% .[['labels']]
# test = data.table::fread("/Users/meelad/phd-projects/sc-lda/data/original-data/invitro-tcells-c393.csv") %>% as_tibble()
#
# channels.metabolism <- c('GLUT1', 'HK2', 'GAPDH', 'LDHA', 'MCT1', 'PFKFB4', 'IDH2', 'CyclinB1',
#                          'GLUD12', 'CS', 'OGDH', 'CytC', 'ATP5A', 'S6_p', 'HIF1A', 'PDK1_p', 'NRF1',
#                          'NRF2_p', 'XBP1', 'VDAC1', 'OPA1', 'DRP1', 'ASCT2', 'CD98', 'GLS', 'GOT2', 'CPT1A',
#                          'ACADM', 'IdU', 'BrU', 'Puromycin', 'H3_p')
# test$cell.id = paste("cell",1:nrow(test),sep="")
# train = test %>% dplyr::filter(grepl("naive",gate.source),grepl("CD8",gate.source), grepl("day0|day4",gate.source)) %>%
#         group_by(gate.source) %>% sample_n(5000) %>% .$cell.id
#
# train.x = test  %>% .[ .$cell.id %in% train, channels.metabolism]
# train.y = test%>% .[ .$cell.id %in% train, "gate.source"] %>% .$gate.source
# hss.result = runHSS(x = train.x, y = factor(train.y), score.method = 'euclidean_2class', downsample = FALSE)



