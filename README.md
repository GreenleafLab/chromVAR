---

# chromVAR

chromVAR is an R package for the analysis of sparse chromatin accessibility data from single cell or bulk ATAC or DNAse-seq data. 

## Note on recent function name changes

The chromVAR (and related motifmatchr) function names recently changed to switch over to camelCase from snake_case.  All exported functions now use camelCase, e.g. `compute_deviations` is now `computeDeviations`. If following the current documentation but using an earlier version of the package, either update the package or be aware of the discrepancy. This change was made to comply with Bioconductor naming preferences, as chromVAR will soon be submitted to Bioconductor.    

## Installation


Installation is easiest using the devtools package. The function `install_github` will install the package.

``` r
devtools::install_github("GreenleafLab/chromVAR")
```

A number of needed packages are installed in this process. One of the dependencies has a system requirement for the gsl library, so if this is not installed already it may need to be installed separately.  

For Windows users, some have reported that the S4Vector dependency does not currently function on windows R 3.3.3, but that installation was successful on R 3.3.2

Two additional packages that are recommended and used in the vignettes:

* motifmatchr - available on [GitHub](https://github.com/GreenleafLab/motifmatchr)
* JASPAR2016  - available from Bioconductor

Additionally, the package chromVARmotifs can be useful for loading additional motif collections:

* chromVARmotifs -  available on [GitHub](https://github.com/GreenleafLab/chromVARmotifs)

## Parallelization

Before running chromVAR functions, it is advisable to use the `register` function from BiocParallel to specify your preferred method of parallelization.  For unix systems, 

```r
library(BiocParallel)
register(MultiCoreParam(8)) # Use 8 cores
```

For Windows, `MultiCoreParam` will not work, but you can use SnowParam:

```r
register(SnowParam(SnowParam(workers=1, type = "SOCK")))
```

Even if you don't want to use more than one core, it is recommended to explicitly register that choice:
```r
register()
```

Please see the documentation for [BiocParallel](https://bioconductor.org/packages/release/bioc/html/BiocParallel.html) for more information about the `register` function and the various options for multi-processing. 


## Quickstart

``` r
library(chromVAR)
library(motifmatchr)
library(BSgenome.Hsapiens.UCSC.hg19)

### Example of how to read in counts -------------------------------------------

# Caution: FAKE FILENAMES -- Replace with real as appropriate! If you want to 
# run on example data in package, start at next section with example_counts data

peakfile <- "mypeaks.bed"
peaks <- getPeaks(peakfile)

bamfiles <- c("mybam1.bam","mybam2.bam")
fragment_counts <- getCounts(bamfiles, peaks, 
                              paired =  TRUE, 
                              by_rg = TRUE, 
                              format = "bam", 
                              colData = DataFrame(celltype = c("GM","K562")))

### ----------------------------------------------------------------------------

### Using example counts from package ------------------------------------------

data(example_counts, package = "chromVAR")
example_counts <- addGCBias(example_counts, 
                              genome = BSgenome.Hsapiens.UCSC.hg19)
counts_filtered <- filterSamples(example_counts, min_depth = 1500,
                                  min_in_peaks = 0.15)
counts_filtered <- filterPeaks(counts_filtered)
motifs <- getJasparMotifs()
motif_ix <- matchMotifs(motifs, counts_filtered,
                         genome = BSgenome.Hsapiens.UCSC.hg19)

# computing deviations
dev <- computeDeviations(object = counts_filtered, 
                                 annotations = motif_ix)

```

See [documentation website](https://greenleaflab.github.io/chromVAR/) for more information!
