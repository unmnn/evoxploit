evoxploit
================

<!-- html_document -->

## Overview

evoxploit is a package for augmenting a **longitudinal
dataset**<sup>\*</sup> with so called **evolution features** (Niemann et
al., 2015) and **sequence features** (Hielscher et al., 2014). Theses
features reflect a study participant’s individual change over time,
his/her change over time in comparison with the cluster he/she was
closest to at the baseline examinations, and the change of whole
participant groups.

<sup>\*</sup> The functions of this package require a n\*m-dimensional
dataframe for the input variables and a factor vector of length n with
the target variable (e.g. a medical outcome). The input data frame
should have the following variables names semantics:

  - Each variable name (e.g. `waist_circum_s0`) consists of a *stem* and
    a *suffix*. The stem is the substring until the last underscore `_`
    (`wairst_circum`). The suffix is the substring beginning with the
    last underscore (`_s0`).
  - The suffix indicates the study wave where the measurement was
    conducted. For instance `s0` could represent the baseline
    examinations of a longitudinal study, `s1` the first follow- up
    examinations, `s2` the second follow up examinations and so on.
  - All variables must have a suffix. The suffix must start with "\_"
    and end with at least one digit. The suffix default is `_s`.

## Installation

``` r
# Development version from Github
devtools::install_github("unmnn/evoxploit")
```

## Usage Example

<!-- The example data `epi` is a list with two elements. `data` is a dataframe  -->

<!-- containing 19 input features -->

``` r
# Inspect example data
library(evoxploit)
str(epi)
```

    ## List of 2
    ##  $ data :Classes 'tbl_df', 'tbl' and 'data.frame':   354 obs. of  19 variables:
    ##   ..$ school_s0     : Factor w/ 3 levels "1","2","3": 2 2 1 3 2 2 3 2 2 3 ...
    ##   ..$ coffee_s0     : num [1:354] 7 4 7 4 3 3 4 2 3 7 ...
    ##   ..$ bmi_s0        : num [1:354] 32 22 23 25 27 23 22 21 17 29 ...
    ##   ..$ bmi_s1        : num [1:354] 30 22 21 25 28 25 22 22 17 29 ...
    ##   ..$ bmi_s2        : num [1:354] 31 26 23 24 30 25 24 26 17 31 ...
    ##   ..$ uric_acid_s0  : num [1:354] 330 183 295 241 253 302 178 204 199 305 ...
    ##   ..$ uric_acid_s1  : num [1:354] 288 199 264 247 276 292 174 207 198 296 ...
    ##   ..$ uric_acid_s2  : num [1:354] 361 187 335 192 289 261 163 223 191 480 ...
    ##   ..$ smoke_s0      : Factor w/ 3 levels "1","2","3": 3 1 2 3 1 2 3 2 2 2 ...
    ##   ..$ smoke_s1      : Factor w/ 3 levels "1","2","3": 3 1 3 3 1 1 3 2 2 2 ...
    ##   ..$ smoke_s2      : Factor w/ 3 levels "1","2","3": 3 1 2 3 1 1 3 2 2 2 ...
    ##   ..$ c07aa_s0      : Factor w/ 2 levels "1","2": 1 1 1 1 1 1 1 1 1 1 ...
    ##   ..$ c07aa_s1      : Factor w/ 2 levels "1","2": 1 1 1 1 1 1 1 1 1 1 ...
    ##   ..$ c07aa_s2      : Factor w/ 2 levels "1","2": 1 1 1 1 1 1 1 1 1 1 ...
    ##   ..$ sa75_s0       : Factor w/ 4 levels "1","2","3","4": 1 1 1 1 1 1 3 1 1 1 ...
    ##   ..$ sa75_s2       : Factor w/ 4 levels "1","2","3","4": 2 1 3 3 1 1 1 1 1 1 ...
    ##   ..$ gx_11597390_s0: Factor w/ 3 levels "1","2","3": 2 3 1 3 3 2 2 2 2 3 ...
    ##   ..$ gx_11597390_s1: Factor w/ 3 levels "1","2","3": 2 3 1 3 3 2 2 2 2 3 ...
    ##   ..$ gx_11597390_s2: Factor w/ 3 levels "1","2","3": 2 3 1 3 3 2 2 2 2 3 ...
    ##  $ label:Classes 'tbl_df', 'tbl' and 'data.frame':   354 obs. of  1 variable:
    ##   ..$ label: Factor w/ 3 levels "A","B","C": 1 1 1 1 1 1 1 2 1 1 ...

``` r
# Create an Evoxploit object
epi_evo <- Evoxploit$new(epi$data, epi$label[[1]], wave_suffix = "_s")

# Print summary to console
summary(epi_evo)
```

    ## === Summary ===
    ## Data:
    ##   - observations:  354 
    ##   - columns:  19 
    ##     - #numeric variables:  7 
    ##     - #factor variables:  12 
    ##   - waves:  0, 1, 2 
    ##  
    ## Target: 
    ##   - #classes:  3 
    ##   - predominant class:  A :  60 % 
    ##  
    ## Clustering: 
    ## Cluster  1 
    ##   - attributes:  bmi_s0, uric_acid_s0 
    ##   - minPts:  6 
    ##   - eps:  0.07289373 
    ##   - clusters:  3 
    ##   - predominant cluster:  1 :  95 % 
    ##   - outlier:  5 % 
    ## Cluster  2 
    ##   - attributes:  bmi_s1, uric_acid_s1 
    ##   - minPts:  6 
    ##   - eps:  0.1415658 
    ##   - clusters:  3 
    ##   - predominant cluster:  1 :  93 % 
    ##   - outlier:  7 % 
    ## Cluster  3 
    ##   - attributes:  bmi_s2, uric_acid_s2 
    ##   - minPts:  6 
    ##   - eps:  0.1245425 
    ##   - clusters:  3 
    ##   - predominant cluster:  1 :  80 % 
    ##   - outlier:  20 % 
    ##  
    ## Evo features: 
    ##   - sequence features:  1 
    ##   - descriptive features:  15 
    ##   - other evo features:  106 
    ## ======

``` r
library(dplyr)
library(ggplot2)
# Calculate Gain Ratio of both original and extracted features w.r.t. the 
# target variable
df_evo <- bind_cols(tibble(label = epi_evo$label), epi_evo$all_features)

gain_names <- FSelector::gain.ratio(label ~ ., data = df_evo)
gain_names %>%
  tibble::rownames_to_column("variable") %>%
  arrange(desc(attr_importance)) %>%
  mutate(variable = forcats::fct_reorder(variable, attr_importance)) %>%
  slice(1:20) %>%
  ggplot(aes(variable, attr_importance)) +
  geom_col() +
  coord_flip()
```

![](man/figures/README-unnamed-chunk-3-1.png)<!-- -->

## Bibliography

  - (Niemann et al., 2015) Uli Niemann, Tommy Hielscher, Myra
    Spiliopoulou, Henry Völzke, and Jens-Peter Kühn. “Can we classify
    the participants of a longitudinal epidemiological study from their
    previous evolution?” *Proc. of IEEE Computer-Based Medical Systems*,
    121-126, 2015.
  - (Hielscher et al. 2014) Tommy Hielscher, Myra Spiliopoulou, Henry
    Völzke, and Jens-Peter Kühn. “Mining longitudinal epidemiological
    data to understand a reversible disorder”. *Proc. of Intelligent
    Data Analysis*, 2014.
