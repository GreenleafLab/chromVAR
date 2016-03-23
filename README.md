
```{r}
devtools::install_github("GreenleafLab/chromVAR")
```

```{r}
library(chromVAR)
BiocParallel::register(BiocParallel::MulticoreParam(16))

bed <- "my_bedfile.bed"

bamfiles <- c('bam1.bam','bam2.bam')

#Reading in inputs
inputs <- get_inputs(bed,bamfiles, by_rg=TRUE, paired=TRUE, format = "bam")

#Finding background peaks
bias <- get_gc(inputs$peaks)
bg <- get_background_peaks(counts_mat = inputs$counts, bias = bias)


#Get motif indices
motifs <- get_motifs()
motif_ix <- get_motif_indices(motifs = motifs, peaks = inputs$peaks)

#Get deviations
deviations <- compute_deviations(counts_mat = inputs$counts, background_peaks = bg, peak_indices = motif_ix)

```