---
[![Build Status](https://travis-ci.org/GreenleafLab/chromVAR.svg?branch=master)](https://travis-ci.org/GreenleafLab/chromVAR)

# chromVAR

chromVAR is an R package for the analysis of sparse chromatin accessibility data from single cell or bulk ATAC or DNAse-seq data. The package aims to identify motifs or other genomic annotations associated with variability in chromatin accessibility between individual cells or samples.  For a more detail overview of the method, please see the [publication](https://www.nature.com/nmeth/journal/vaop/ncurrent/full/nmeth.4401.html) ([pdf](http://greenleaf.stanford.edu/assets/pdf/nmeth.4401.pdf), [supplement](https://drive.google.com/file/d/0B8eUn6ZURmqvUjBCbE5Hc0p4UFU/view?usp=sharing)). 

For a paper evaluating chromVAR and other methods as a method for enabling clustering of single cells, see [the preprint from Huidong Chen et al](https://www.biorxiv.org/content/10.1101/739011v1). Using kmers + PCA appears to be the best variant of chromVAR for clustering, but newer methods such as [SnapATAC](https://github.com/r3fang/SnapATAC) outperform chromVAR for the clustering tasks evaluated in the paper. chromVAR may be complementary to some other methods, as a way of annotating TF motif usage in cells & clusters rather than cluster identification or embedding.

## Installation

The recommended installation method for `chromVAR` is using BiocManager. You will first have to have installed the [BiocManager package
](https://cran.r-project.org/package=BiocManager).

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("chromVAR")
```

A number of needed packages are installed in this process. One of the dependencies has a system requirement for the gsl library, so if this is not installed already it may need to be installed separately. Several people have reported issues with the GO.db package (a dependency of one of the dependencies) not being installed automatically -- if you see an error relating to that package, try installing it separately first (`BiocManager::install("GO.db")`).

For Windows users, some have reported that the S4Vector dependency does not currently function on windows R 3.3.3, but that installation was successful on R 3.3.2.

For Mac users, some have encountered installation difficulties relating to compiling the C++ code. If you encounter problems, please see the threads and advice in Issues [11](https://github.com/GreenleafLab/chromVAR/issues/11) and [20](https://github.com/GreenleafLab/chromVAR/issues/20).  

Two additional packages that are recommended and used in the vignettes:

* motifmatchr - available on [GitHub](https://github.com/GreenleafLab/motifmatchr) or [development version of Bioconductor](https://bioconductor.org/packages/devel/bioc/html/motifmatchr.html)
* JASPAR2016  - available from Bioconductor

Additionally, the package chromVARmotifs can be useful for loading additional motif collections:

* chromVARmotifs -  available on [GitHub](https://github.com/GreenleafLab/chromVARmotifs)

## Parallelization

Before running chromVAR functions, it is advisable to use the `register` function from BiocParallel to specify your preferred method of parallelization.  For unix systems, 

```r
library(BiocParallel)
register(MulticoreParam(8)) # Use 8 cores
```

For Windows, `MulticoreParam` will not work, but you can use SnowParam:

```r
register(SnowParam(SnowParam(workers=1, type = "SOCK")))
```

Even if you don't want to use more than one core, it is recommended to explicitly register that choice using SerialParam:

```r
register(SerialParam())
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
