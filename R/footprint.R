
#' @export
get_footprint <- function(motif_pos, ends, flank = 250){
  
  tmp = resize(motif_pos, width = 1, fix = "center", ignore.strand = FALSE)
  o = findOverlaps(motif_pos, ends, maxgap = flank, ignore.strand = TRUE)
  d = ifelse(as.character(strand(tmp[queryHits(o)])) == "-", 
             start(tmp[queryHits(o)]) - start(ends[subjectHits(o)]), 
             start(tmp[queryHits(o)]) - start(ends[subjectHits(o)]))
  out = tabulate2(d, min_val = -flank, max_val = flank)
  return(out)
}