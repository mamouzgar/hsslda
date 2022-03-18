LDA with Hybrid Subset Selection (HSS)
======================================

Authors: Meelad Amouzgar and David Glass

Bendall lab @ Stanford University

This repository is still under development.

Introduction:
-------------

Linear Discriminant Analysis (LDA) is a classification algorithm that
we’ve repurposed for supervised dimensionality reduction of single-cell
data. LDA identifies linear combinations of predictors that optimally
separate a priori labels, enabling users to tailor visualizations to
separate specific aspects of cellular heterogeneity. We implement LDA
with feature selection by Hybrid Subset Selection (HSS) in this R
package, called hsslda.

Installation instructons:
-------------------------

You can install hsslda using:

``` r
remotes::install_github("mamouzgar/hsslda", build_vignettes = FALSE)
```

The github page includes an introduction to the package.
``` r
library(hsslda)
```

Here we read-in the example data

``` r
data(TcellHartmann2020_sampleData)
```

``` r
colnames(TcellHartmann2020_sampleData)
```

    ##  [1] "GLUT1"    "HK2"      "GAPDH"    "LDHA"     "MCT1"     "PFKFB4"  
    ##  [7] "IDH2"     "CyclinB1" "GLUD12"   "CS"       "OGDH"     "CytC"    
    ## [13] "ATP5A"    "S6_p"     "HIF1A"    "labels"

``` r
head(TcellHartmann2020_sampleData,3)
```

    ##        GLUT1        HK2     GAPDH      LDHA       MCT1    PFKFB4      IDH2
    ## 1 0.03409762 0.00000000 0.3481982 0.4066992 0.06882467 0.1601350 0.5051168
    ## 2 0.25528617 0.06572665 0.3005934 0.5078840 0.01220142 0.0241186 0.6933922
    ## 3 0.14032230 0.00000000 0.1665477 0.1122911 0.12135540 0.1192041 0.4822793
    ##      CyclinB1    GLUD12        CS      OGDH       CytC     ATP5A       S6_p
    ## 1 0.063481520 0.4324707 0.5296415 0.5315838 0.20964937 0.3771783 0.03704562
    ## 2 0.102636694 0.4726154 0.5879990 0.4161547 0.02309679 0.6053693 0.16939760
    ## 3 0.007619569 0.3998965 0.4228455 0.1749531 0.08402797 0.2401586 0.03365913
    ##        HIF1A labels
    ## 1 0.02346801   day0
    ## 2 0.08144526   day0
    ## 3 0.09690983   day0

How to run HSS-LDA:
-------------------

You can run hybrid subset selection on a dataset using runHSS().

runHSS takes 3 inputs:

1.  x: a dataframe of cells (rows) and markers/genes (columns)

2.  y: a vector of your class labels of interest.

3.  score.method: A scoring metric to use for feature selection, which
    can include: ‘euclidean’, ‘silhouette’, ‘pixel.entropy’,
    ‘pixel.density’, or ‘custom’

4.  (optional) downsample: Boolean, defaults to TRUE. Downsamples data
    to improve runtime. Set to FALSE to use all input data.

Here we will run HSS using the default scoring method: euclidean
distance.Note that scoring metrics like silhouette or pixel class
entropy (PCE) score will require additional package dependencies.

``` r
channels = c('GLUT1', 'HK2', 'GAPDH', 'LDHA', 'MCT1', 'PFKFB4', 'IDH2', 'CyclinB1', 'GLUD12', 'CS', 'OGDH', 'CytC', 'ATP5A', 'S6_p', 'HIF1A')
train.x = TcellHartmann2020_sampleData[channels]
train.y = TcellHartmann2020_sampleData[['labels']]
hss.result = runHSS(x = train.x, y = train.y, score.method = 'euclidean', downsample = FALSE)
```

    ## Using all input cells for HSS-LDA...

    ## Optimizing dimensionality reduction using euclidean metric...

    ## Number of markers: 2

    ## Number of markers: 3
    ## Number of markers: 3

    ## Number of markers: 4
    ## Number of markers: 4

    ## Number of markers: 5
    ## Number of markers: 5

    ## Number of markers: 6
    ## Number of markers: 6

    ## Number of markers: 7
    ## Number of markers: 7

    ## Number of markers: 8
    ## Number of markers: 8

    ## Number of markers: 9
    ## Number of markers: 9

    ## Number of markers: 10

    ## Number of markers: 9

    ## Number of markers: 10
    ## Number of markers: 10

    ## Number of markers: 11

    ## Number of markers: 10

    ## Number of markers: 11
    ## Number of markers: 11

    ## Number of markers: 12

    ## Number of markers: 11

    ## Number of markers: 12
    ## Number of markers: 12

    ## Number of markers: 13

    ## Number of markers: 12

    ## Number of markers: 13
    ## Number of markers: 13

    ## Number of markers: 14

    ## Number of markers: 13

    ## Number of markers: 14
    ## Number of markers: 14

    ## Number of markers: 15

    ## Number of markers: 14

    ## Number of markers: 15

Output summary:
---------------

The final hss.result object outputted from runHSS contains a few
elements, the most important ones are highlighted below:

1.  HSSscores: a table of all scores for each subsetted model.

2.  ElbowPlot: visualizes the elbow plot and calculated elbow point for
    the optimal feature set.

3.  HSS-LDA-model: The final LDA model using the optimal feature set.

4.  HSS-LDA-result: The final result table containing all initially
    inputted markers, class labels, and linear discriminants from the
    opimal model.

5.  downsampled.cells.included.in.analysis: A boolean vector of rows
    included in the downsampling analysis. If downsampling was not
    performed, all values are TRUE. The final output results includes
    all data projected onto the HSS-LD axes.

You can visualize the elbow plot, which is ggplot configurable:

``` r
hss.result$ElbowPlot
```

![](/private/var/folders/_p/dzrkxwzd30l1jx40q6_p26gj_sqz2b/T/RtmpeTiudj/preview-669969b2ba1d.dir/hsslda-intro-md_files/figure-markdown_github/unnamed-chunk-9-1.png)

The final output object contains a merged dataframe of your markers,
class labels, and newly generated LD axes

``` r
lda.df = hss.result$`HSS-LDA-result`
head(lda.df, 3)
```

    ##        GLUT1        HK2     GAPDH      LDHA       MCT1    PFKFB4      IDH2
    ## 1 0.03409762 0.00000000 0.3481982 0.4066992 0.06882467 0.1601350 0.5051168
    ## 2 0.25528617 0.06572665 0.3005934 0.5078840 0.01220142 0.0241186 0.6933922
    ## 3 0.14032230 0.00000000 0.1665477 0.1122911 0.12135540 0.1192041 0.4822793
    ##      CyclinB1    GLUD12        CS      OGDH       CytC     ATP5A       S6_p
    ## 1 0.063481520 0.4324707 0.5296415 0.5315838 0.20964937 0.3771783 0.03704562
    ## 2 0.102636694 0.4726154 0.5879990 0.4161547 0.02309679 0.6053693 0.16939760
    ## 3 0.007619569 0.3998965 0.4228455 0.1749531 0.08402797 0.2401586 0.03365913
    ##        HIF1A        LD1       LD2      LD3      LD4         LD5 labels
    ## 1 0.02346801  0.1824071 -1.505745 4.121023 3.561247  0.17292992   day0
    ## 2 0.08144526 -1.4912720 -1.858684 4.190646 4.693052 -2.52678835   day0
    ## 3 0.09690983 -1.0026704 -2.758428 2.453262 3.370301  0.09091075   day0

``` r
ggplot2::ggplot(lda.df,ggplot2::aes(x = HSS_LD1, y = HSS_LD2, color = labels)) +
  ggplot2::geom_point() 
```

![](/private/var/folders/_p/dzrkxwzd30l1jx40q6_p26gj_sqz2b/T/RtmpeTiudj/preview-669969b2ba1d.dir/hsslda-intro-md_files/figure-markdown_github/unnamed-chunk-12-1.png)

The final HSS-LDA model is also saved, and can be used like any R model.

``` r
hss.result$`HSS-LDA-model`
```

    ## Call:
    ## lda(y ~ ., data = x[, markers])
    ## 
    ## Prior probabilities of groups:
    ##      day0      day1      day2      day3      day4      day5 
    ## 0.1666667 0.1666667 0.1666667 0.1666667 0.1666667 0.1666667 
    ## 
    ## Group means:
    ##          GLUT1        HK2      IDH2   CyclinB1    GLUD12        CS      OGDH
    ## day0 0.1686856 0.04104809 0.6052301 0.06537171 0.4897904 0.5299507 0.3292659
    ## day1 0.2953167 0.21629789 0.5935709 0.06965901 0.4730912 0.5417745 0.3906114
    ## day2 0.5763669 0.55772316 0.6680269 0.23513566 0.5277285 0.7177145 0.5890289
    ## day3 0.7184115 0.55541671 0.7248705 0.35324105 0.5967461 0.7861387 0.6456323
    ## day4 0.7027508 0.41256773 0.7639317 0.31064447 0.6161448 0.7786976 0.5322304
    ## day5 0.7213829 0.31763930 0.7655344 0.25145793 0.6068985 0.7636594 0.4400868
    ##           CytC     ATP5A     HIF1A
    ## day0 0.1369610 0.4328304 0.1148124
    ## day1 0.1670760 0.5025277 0.2048914
    ## day2 0.2500147 0.6946313 0.2941883
    ## day3 0.3662381 0.7670108 0.3772618
    ## day4 0.3900886 0.7277347 0.2854079
    ## day5 0.4314735 0.6504328 0.2555777
    ## 
    ## Coefficients of linear discriminants:
    ##                 LD1         LD2        LD3       LD4         LD5
    ## GLUT1    -9.3286320 -2.42021822 -4.2104736 -2.427140   0.5830266
    ## HK2       0.4317694  3.60630945 -1.9692387  2.802069   1.0771403
    ## IDH2      0.9131571 -6.23232819  0.2676104  6.611837  -2.5128133
    ## CyclinB1  0.4159760  0.89607252  3.3770898 -1.120539   0.9637185
    ## GLUD12    3.2002486  0.08886794  1.2587234 -2.408197   2.8570995
    ## CS       -4.3611705 -2.34919265  3.6521969  6.055976   4.8876598
    ## OGDH      1.2117074  2.79674730  0.5939495 -1.434667   3.1181463
    ## CytC     -0.8575323 -2.13438292 -0.2152454 -1.984046   0.1923797
    ## ATP5A     1.2688525  4.81834552  3.0803171 -1.425575 -11.1020405
    ## HIF1A    -0.1810920  0.73022516  0.1653555 -3.167533   1.1239380
    ## 
    ## Proportion of trace:
    ##    LD1    LD2    LD3    LD4    LD5 
    ## 0.7906 0.1890 0.0115 0.0066 0.0023

You can use this final model to apply the same dimensionality reduction
to new data using makeAxes.

``` r
lda.df.newData = makeAxes(df = newdata, co =hss.result$`HSS-LDA-model`$scaling)
```

You can also add your own custom separation metric to perform feature
selection with by changing score.method = ‘custom’ and adding custom
function. The custom function must take in as input:

1.  x: dataframe of all your data.

2.  y: vector of class labels matching training data rows.

3.  custom.score.method: the function for your custom separation metric.

You can then input your custom function into the hsslda::runHSS()
function as an argument into custom.score.method, and indicate
score.method = ‘custom’

``` r
hss.resultCustom = runHSS(x = train.x, y = train.y, score.method = 'custom', custom.score.method = myCustomFunction)
```

And that’s how you use LDA with Hybrid Subset Selection!

Questions, comments, or general feedback:

<a href="mailto:amouzgar@stanford.edu" class="email">amouzgar@stanford.edu</a>
