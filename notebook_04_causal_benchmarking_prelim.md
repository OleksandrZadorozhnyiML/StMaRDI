***Description of the dataset***

The Tuebingen cause-effect dataset contains the 108 datasets, every with
2 variables. The description of the data is the following:

the data represented here is a growing database with different data for
testing causal detection algorithms. The goal here is to distinguish
between cause and effect. We searched for data sets with known ground
truth. However, we do not guarantee that all provided ground truths are
correct. The datafiles are .txt-files and contain two variables, one is
the cause and the other the effect. For every example there exists a
description file where you can find the ground truth and how the data
was derived.

Note that not always the first column is the cause and the second the
effect. This is indicated in a meta-data file. There is also a weighting
factor suggested for some pairs which are very similar if you want to
calculate the overall performance.

***Setup of the problem*** We have 89 datasets with variables
*V**a**r*1, *V**a**r*2. The goal is to establish cause-effect relationships between 2 variables using the data. 

``` r
### loading packages
library("liver")
library("qgraph")
library("igraph")
library("bnlearn")
library("dplyr")

list_of_datasets = read.csv("whitelist.csv")
```

__Loading the Tuebingen datasets__

``` r
whitelist = as.list(list_of_datasets)
head(list_of_datasets)
```

    ##   whitelist
    ## 1         1
    ## 2         2
    ## 3         3
    ## 4         4
    ## 5         5
    ## 6         6

``` r
# randomly pick 5 datasets for study comparison 
selected = list_of_datasets[sample(nrow(list_of_datasets),5),]
print(selected)
```

    ## [1] 34 50 61 64  9

__Here we taking a glimpse in the tubingen data__

``` r
# here we extract the selected datasets to perform benchmarking

sel_to_str = as.character(selected)
elem_len= nchar(sel_to_str)
```

__Transcribing the lengths of the array to form the file names__

``` r
result = vector(mode="character", length=5)
for (i in 1:5){
  text_int = "causal_tubingen"
  elem_len = nchar(sel_to_str[i])
  result[i] <- switch(elem_len,
       "1" = paste(text_int,"00",sel_to_str[i],sep = ""),
       "2" = paste(text_int,"0",sel_to_str[i],sep = ""),
       "3" = paste(text_int,"",sel_to_str[i],sep = ""),
)
}
```

``` r
library(ggplot2)

set.seed(1)

# The data used can be obtained from https://webdav.tuebingen.mpg.de/cause-effect/. Here we collect automatically the names and descriptions. 
filenames <- sprintf("https://webdav.tuebingen.mpg.de/cause-effect/pair%04d.txt", seq(1, 108, 1))
descriptions <- sprintf("https://webdav.tuebingen.mpg.de/cause-effect/pair%04d_des.txt", seq(1, 108, 1))
```

We are considering the datasets which have 2 variables: in the following
we use the HSIC criterion to test for independence of the residuals and
covariates

``` r
library("utils")

results <- data.frame(matrix(nrow = 0, ncol = 2))
colnames(results) <- c("pair", "HSIC test")

# dropping the datasets due to 
drop_f = c(17,64:68)
p_val = 0.05
library(krr)
library(dHSIC)
n = 100 

for (i in setdiff(whitelist$whitelist,drop_f)) {
  # Ground truth files

  if (i %in% c(86, 88)) {
    expected = c(1, 2)
  } else {
    description_file  = url(sprintf("https://webdav.tuebingen.mpg.de/cause-effect/pair%04d_des.txt", 1))
    # Pair 86: x to y, 88 x to y
    description1 = paste(readLines(description_file), collapse=" ")
    description1 = gsub(" ", "", description1)
    if (grepl(">[Yy]|[Yy]<", description1)) {
      expected = c(1, 2)
    } else if (grepl(">[Xx]|[Xx]<", description1)) {
      expected = c(2, 1)
    } else {
      print("This shouldn't happen")
    }
    # print(expected)
  }
  
  # data_file <- url(sprintf("https://webdav.tuebingen.mpg.de/cause-effect/pair%04d.txt", i))
    cause_effect_pair = read.table(filenames[i], header=F)
    cause_effect_pair = as.data.frame(as.matrix(cause_effect_pair[1:min(nrow(cause_effect_pair),n),]))
    var1  = as.matrix(cause_effect_pair[,1])
    var2  = as.matrix(cause_effect_pair[,2])
    # LiNGAM
    # Residual analysis based on HSIC test
    
    # Fit the ridge regression model
    fit_1 = krr(var1, var2, sigma = 1, lambda = 1)
    #predict_1 =  predict(fit_1, xnew =x1_seq )

    fit_2 = krr(var2, var1, sigma = 1, lambda = 1)
    #predict_2 = predict(fit_2,xnew  = x2_seq)
    
    # Calculate the residuals
    residuals_1 = residuals(fit_1, type = "deviance")
    residuals_2 = residuals(fit_2, type = "deviance")
    
    
    hsic_12 = dhsic.test(residuals_1,var1)
    hsic_21 = dhsic.test(residuals_2,var2)
    
   #  print(sprintf("Pair%04d",i))
   #  print(hsic_12$p.value)
   #  print(hsic_21$p.value)
    
    if (hsic_12$p.value<hsic_21$p.value){
      result = c(1,2)
    }
    else if (hsic_21$p.value<hsic_12$p.value){
      result = c(2,1)
    }
    else {
      result = c(1,2)
    }
    # Results
    results[nrow(results) + 1,] <- c(i, all(expected == result))
}
```

__Refererences__

\[1\] Pfister, N, and Peter, J,. Independence Testing via Hilbert
Schmidt Independence Criterion. <https://doi.org/10.1111/rssb.12235>

\[2\] <https://vulstats.ucsd.edu/chi-squared.html>

\[3\] A. Gretton, K. Fukumizu, C. H. Teo, L. Song, B. Sch'’olkopf, and
A. Smola. A kernel statistical test of independence. In Advances in
Neural Information Processing Systems 20 (NIPS), 2008. SBN:
978-1-60560-949-2

\[4\] Peters, J., Mooij, J.M., Janzing, D. and Sch"olkopf, B., 2014.
Causal discovery with continuous additive noise models. JMLR
15(58):2009−2053, 2014.
