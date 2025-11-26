
<!-- README.md is generated from README.Rmd. Please edit that file -->

# CustomBoundary

<!-- badges: start -->
<!-- badges: end -->

The goal of CustomBoundary is to generate required input files for the
CASTool when using a custom (non-state) boundary.

## Installation

You can install the development version of CustomBoundary from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("laura-naslund/CustomBoundary")
```

## Generate files

To generate files, provide the file paths 1) the location of the custom
boundary geospatial file (e.g., .shp or .geojson), 2) desired output
directory, 3) the region name that corresponds to the custom boundary.

generateFiles(inputFilePath, outputFolder, region_name)

Additional parameters can be provided to the clustering algorithm script
(e.g., to override the default settings).
