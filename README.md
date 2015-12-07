# chromVAR
chromatin Variability Across Regions (of the genome!)

[] indicates status {} indicates estimated work time to completion {} indicates priority

Reading in and preprocessing data:
- Read processed bed file [100%]
- Read bam file with RG tags [100%]
- Read multiple bam files [100%]
- read unprocessed peak bed file and find non-overlapping, even-widthed peaks [0%] {1 day} (medium priority)
- find peaks using databank of peaks? [0%] {3+ days} (very low priority)
- run preprocessing of bam files using shell scripts? [0%] {2 days} (low priority)
- get motif indices [100%]
- get kmer indices [100%]
- get gapped kmer indices [80%]

Computing variability
- compute variability for individual peaks (applicable for bulk samples) [90%] {3 hours}
- compute variability for sets of peaks [100%]

Basic extensions for variability
- find samples with high or low accessibility for set of peaks [50%] {1 day} (high priority)
- find peaks with high or low accessibility across set of samples [50%] {1 day}  (high priority)

More complex extensions
- find groups of annotations/motifs/kmers with independent variability [25%] {?}
- for groups of kmers, group kmers that likely make up the same motif [60%] {?}
- clustering of samples based on variability [30%]{?}

Differential variablity
- find sets with greater variability in one sample versus another [0%]{?}

Visualization
- plot variability [90%] {1-2 hours} (medium priority)
- visualize kmers [60%] {2 days}
- heatmaps deviatons across samples [50%] {4 hours}

Documentation 
- functions [60%] {1-3 hours}
- classes [60%] {1-3 hours}
- vignette [0%] {1-3 days}








