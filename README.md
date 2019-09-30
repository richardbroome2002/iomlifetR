# iomlifetR
# add documentation here
An implementation of the IOMLIFET excel spreadsheets in R. The IOMLIFET spreadsheets facilitate the use of life-tables for
assessment of public health risks.

Key functions are:

1. `life_table()` which calculates life expectancy using the Chiang method
1. `burden_le()` which calls `life_table()` to produce life tables for the baseline and reduced exposure scenarios
1. `impact()` which uses leslie matrices to form the future population under the two scenarios. 

The latter function provides estimates of life-years gained in the future from a popoulation-wide reduction in PM2.5 related risk of mortality. 

```
Version: 0.2
Author: Richard Broome, Joshua Horsley, Ivan Hanigan
Maintainer: Ivan Hanigan <ivan.hanigan@gmail.com>
License: GPL-3
```

## Installation notes

- Vignettes are being added. These are not built by devtools by default, need to use 

```r
library(devtools)
install_github("richardbroome2002/iomlifetR", build_vignettes = TRUE)
```
