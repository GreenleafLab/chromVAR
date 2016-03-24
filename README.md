# chromVAR

##Installation
Installation is easiest using the devtools package.  The function `install_github` will install the package.
```{r}
devtools::install_github("GreenleafLab/chromVAR", auth_token = "my_token")
```
The argument auth_token takes in your github [personal acces token](https://github.com/settings/applications).  This token is needed because at the moment this repository is private.  

A number of needed packages are installed in this process. Note that for functions that require a genome sequence, the package BSgenome.Hsapiens.UCSC.hg19 is used as a default argument. However that package will not be automatically installed -- if using the default argument and that genome build, you will need to install that package.  If using another genome build, the appropraiate BSgenome object for your species should be passed to functions requiring a genome build (e.g. `get_motif_indices`, `get_gc`, and `get_kmer_indices`).

Depending on your repository settings, the Bioconductor dependencies may fail to install.  Use `setRepositories(graphics=F)` to see what repositories you have activated and to add the BioC software repository if need be.  

##Loading the package
Use library or require to load package.
```{r}
library(chromVAR)
```

##Setting multiprocessing options
The package uses BiocParallel to do the multiprocessing.  Check the documentation for BiocParallel to see available options.  The settings can be set using the register function.  For example, to use MulticoreParam with 16 cores:
```{r}

BiocParallel::register(BiocParallel::MulticoreParam(16))

```
##Reading in inputs
```{r}
bed <- "my_bedfile.bed"
bamfiles <- c('bam1.bam','bam2.bam')

peaks <- get_peaks(bed)
counts <- get_counts(bamfiles, peaks, paired =  TRUE, by_rg = TRUE, format = "bam")
```

The function `get_inputs` returns a list with two elements.  The first is a GenomicRanges object with the peaks.  The second is a Matrix of fragment counts per sample/cell for each peak.  The Matrix package is used so that if the matrix is sparse, the matrix will be stored as a sparse Matrix.  If you want to manipulate the counts data in any way in addition to use by chromVAR packages, it is a good idea to load the Matrix package (`library(Matrix)`).

##Filtering inputs

If working with single cell data, it is advisable to filter out samples with insufficient reads or a low proportion of reads in peaks as these may represent empty wells or dead cells. 

```{r}
#compute sequencing depth from bamfiles
depths = get_sample_depths(bamfiles, paired =  TRUE, by_rg = TRUE, format = "bam")

#find indices of samples to keep
samples_to_keep = filter_samples(counts, depths)

#filter counts matrix
counts = counts[,samples_to_keep]

```

Two parameters are used for filtering -- min_in_peaks and min_depth.  If not provided (as above), these cutoffs are estimated based on the medians from the data.  min_in_peaks is set to 0.5 times the median proportion of fragments in peaks.  min_depth is set to the maximum of 500 or 10% of the median library size. 

Unless `plot = FALSE` given as argument to function `filter_samples`, the following type of plot will be generated:
![proportion_in_peaks_vs_depth_plot](example_plot1.png)

For both bulk and single cell data, peaks should be filtered based on having at least a certain number of fragments. At minimum, each peak should have at least one fragment across all the samples (it might be possible to have peaks with zero reads due to using a peak set defined by other data). Otherwise, downstream functions won't work. 

```{r}
peaks_to_keep = filter_peaks(counts)

counts = counts[peaks_to_keep,]
peaks = peaks[peaks_to_keep] #to keep peaks consistent with counts
```

##Finding background peaks
First, compute the gc content from the peaks.  Then use that output as well as the counts to determine background peaks.  
```{r}
bias <- get_gc(inputs$peaks)
bg <- get_background_peaks(counts_mat = inputs$counts, bias = bias)
```

Note that the function `get_gc` also takes in an argument for a BSgenome object.  The default is BSgenome.Hsapiens.UCSC.hg19, so if using a different genome build be sure to provide the correct genome. For example, if using sacCer3 you could do:
```{r}
library(BSgenome.Scerevisiae.UCSC.sacCer3)
bias <- get_gc(inputs$peaks, genome = BSgenome.Scerevisiae.UCSC.sacCer3)
```

Check out `available.genomes` from the BSgenome package for what genomes are available. For making your own BSgenome object, check out `BSgenomeForge`.  

##Get motifs and what peaks contain motifs
The function `get_motifs` fetches motifs from the JASPAR database.  The function `get_motif_indices` finds which peaks contain which motifs.
```{r}
motifs <- get_motifs()
motif_ix <- get_motif_indices(motifs = motifs, peaks = inputs$peaks)
```

The function get_motifs() by default gets human motifs from JASPAR core database.  For other species motifs, change the species argument.  
```{r}
motifs <- get_motifs(species = "Saccharomyces cerevisiae")
```
For using a collection other than core, use the `core` argument.  Options include: "CORE", "CNE", "PHYLOFACTS", "SPLICE", "POLII", "FAM", "PBM", "PBM_HOMEO", "PBM_HLH".

The `get_motifs` function is simply a wrapper around `getMatrixSet` from TFBSTools-- you can also use that function to fetch motifs from JASPAR if you prefer, and/or check out the documentation for that function for more information.  

For the function `get_motif_indices` a genome sequence is again required.  So for sacCer3 for example:

```{r}
motif_ix <- get_motif_indices(motifs = motifs, peaks = inputs$peaks, genome = BSgenome.Scerevisiae.UCSC.sacCer3)
```

Another option is the p.cutoff for determing how stringent motif calling should be. The default value is 0.00005, which tends to give reasonable numbers of motif matches.  

##Compute deviations
```{r}
deviations <- compute_deviations(counts_mat = inputs$counts, background_peaks = bg, peak_indices = motif_ix)
```

The function `compute_deviations` returns a list with two matrices. The first matrix (deviations$z if using command above) will give the deviation z-score for each set of peaks (rows) for each cell or sample (columns).  These z-scores represent how accessible the set of peaks is relative to the expectation based on equal chromatin accessibility profiles across cells/samples, normalized by a set of background peak sets matched for GC and average accessability.   The second matrix (deviations$fc) will give the log2 fold change in accessibility for each set of peaks relative to the expected, again normalized by the set of background peaks.  


