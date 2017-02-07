---

# chromVAR

chromVAR is an R package for the analysis of sparse chromatin accessibility data from single cell or bulk ATAC or DNAse-seq data. 

## Installation


Installation is easiest using the devtools package. The function `install_github` will install the package.

``` r
devtools::install_github("GreenleafLab/chromVAR", auth_token = "my_token")
```

The argument auth\_token takes in your github [personal access token](https://github.com/settings/tokens). This token is needed because at the moment this repository is private.

Two additional packages that are recommended and used in the vignettes:

* motifmatchr - available on [GitHub](https://github.com/GreenleafLab/motifmatchr)
* JASPAR2016  - available from Bioconductor

Depending on your repository settings, the Bioconductor dependencies may fail to install. Use `setRepositories(graphics=FALSE)` to see what repositories you have activated and to add the BioC software repository if need be. You may also have to install some of the BioConductor packages manually, and possible some non-R dependencies for some of the Bioconductor packages.  



## Quickstart

``` r
library(chromVAR)
library(motifmatchr)

### Example of how to read in counts -------------------------------------------

# Caution: FAKE FILENAMES -- Replace with real as appropriate! If you want to 
run on example data in package, start at next wection with "example_counts" data

peakfile <- "mypeaks.bed"
peaks <- get_peaks(peakfile)

bamfiles <- c("mybam1.bam","mybam2.bam")
fragment_counts <- get_counts(bamfiles, peaks, 
                              paired =  TRUE, 
                              by_rg = TRUE, 
                              format = "bam", 
                              colData = DataFrame(celltype = c("GM","K562")))

### ----------------------------------------------------------------------------

### Using example counts from package ------------------------------------------

data(example_counts, package = "chromVAR")
example_counts <- add_gc_bias(example_counts)
counts_filtered <- filter_samples(example_counts, min_depth = 1500,
                                  min_in_peaks = 0.15)
counts_filtered <- filter_peaks(example_counts)
motifs <- get_jaspar_motifs()
motif_ix <- match_motifs(motifs, counts_filtered)

# computing deviations
dev <- compute_deviations(object = counts_filtered, 
                                 annotations = motif_ix)

```

See Vignettes for more information!
