library(chromVAR)
'%ni%' = Negate('%in%')

peak_files <- c("/lab/alicia/chromVARmanuscript/data-raw/scATAC/peaks/GM_peaks.narrowPeak",
                "/lab/alicia/chromVARmanuscript/data-raw/scATAC/peaks/H1_peaks.narrowPeak")

peaks <- lapply(peak_files, read_macs2_narrowpeaks)

cleanup_peaks <- function(peaks){
  peaks <- peaks[which(as.character(seqnames(peaks)) %ni% c("chrM","chrY","chrX"))]
  
  peaks <- sort(peaks)
  
  mito = get_peaks(filename = "/lab/alicia/chromVARmanuscript/data-raw/hg19_mito_blacklist.txt")
  
  blacklisted = overlapsAny(peaks, mito)
  non_mito_peaks = which(!blacklisted)
  peaks <- peaks[non_mito_peaks]
  
  seqlengths(peaks) <- GenomeInfoDb::seqlengths(BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19)[GenomeInfoDb::seqlevels(peaks)]
  within_boundaries = which(GenomicRanges::trim(peaks) == peaks)
  peaks <- peaks[within_boundaries]
  return(peaks)
}
peaks_top <- lapply(peaks, function(x){
  filt = cleanup_peaks(x)
  filt[order(filt$score,decreasing = TRUE)[1:20000]]
})

peaks_merged <- sort(do.call(c,peaks_top))

keep_peaks = 1:length(peaks_merged)
while (!(isDisjoint(peaks_merged[keep_peaks]))){
  chr_names = as.character(seqnames(peaks_merged[keep_peaks]))
  starts = start(peaks_merged[keep_peaks])
  ends = end(peaks_merged[keep_peaks])
  overlap_next = intersect(which(chr_names[1:(length(keep_peaks) -1)] == chr_names[2:(length(keep_peaks))]),
                           which(ends[1:(length(keep_peaks) -1)] >= starts[2:(length(keep_peaks))]))
  overlap_previous = overlap_next + 1
  overlap_comparison = peaks_merged[keep_peaks[overlap_previous]]$qval > peaks_merged[keep_peaks[overlap_next]]$qval
  discard = keep_peaks[c(overlap_previous[!overlap_comparison], overlap_next[overlap_comparison])]
  keep_peaks = keep_peaks[keep_peaks %ni% discard]
}

scatac_peaks <- peaks_merged[keep_peaks]


bam_info = read.table("/lab/alicia/chromVARmanuscript/data-raw/scATAC/bam/scatac_bam_info.txt",sep = "\t", header = TRUE,
                      stringsAsFactors = FALSE)

bam_info = bam_info[c(3,5),]

example_counts <- get_counts(paste0("/lab/alicia/chromVARmanuscript/data-raw/scATAC/bam/",bam_info$Bamfile), scatac_peaks, paired = TRUE, by_rg = TRUE, format = "bam",
                             colData = DataFrame(Cell_Type = bam_info$Cell_Type))

example_counts <- example_counts[,c(which(colData(example_counts)$Cell_Type == "GM")[1:25],
                                    which(colData(example_counts)$Cell_Type == "H1")[1:25])]

example_counts <- filter_peaks(example_counts)



devtools::use_data(example_counts, overwrite = TRUE)

mini_counts <- filter_samples(example_counts, min_depth = 2000,
                              min_in_peaks = 0.15, shiny = FALSE)
mini_counts <- filter_peaks(mini_counts)
mini_counts <- mini_counts[sample(nrow(mini_counts), 1000),]
mini_counts <- add_gc_bias(mini_counts)
devtools::use_data(mini_counts, overwrite = TRUE)

example_counts <- add_gc_bias(example_counts)
counts_filtered <- filter_samples(example_counts, min_depth = 1500,
                                  min_in_peaks = 0.15, shiny = FALSE)
counts_filtered <- filter_peaks(example_counts)
motifs <- get_jaspar_motifs()
motif_ix <- match_motifs(motifs, counts_filtered)

# computing deviations
dev <- compute_deviations(object = counts_filtered, 
                          annotations = motif_ix)



