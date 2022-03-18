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
#' @importFrom ggplot2 theme
#'
#' @importFrom dplyr filter
#' @importFrom dplyr select
#' @importFrom dplyr group_by
#' @importFrom dplyr ungroup
#' @importFrom dplyr rename
#' @importFrom dplyr summarize
#' @importFrom dplyr n
#' @importFrom dplyr %>%
#' @importFrom magrittr %>%
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



        ## specify HSS_LD outputs ##
        hss.result[["method"]] = score.method
        hss.result[["finalMarkers"]] = markers
        hss.result[["HSSscores"]] = results
        hss.result[["ElbowPlot"]] = plotElbow(results = results, elbow = elbow)
        ## explicitly lda.out results as HSS-LDA for the user
        colnames(lda.out$scaling) = paste0("HSS_", colnames(lda.out$scaling), sep ='')
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
        # colnames(HSS.LDs) = paste0("HSS_", colnames(HSS.LDs), sep = "")
        df = cbind(df, HSS.LDs)
        return(df)
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





