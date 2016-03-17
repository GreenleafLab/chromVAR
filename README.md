# chromVAR
chromatin Variability Across Regions (of the genome!)

[] indicates status  

Reading in and preprocessing data:
- Read processed bed file [100%]
- Read bam file with RG tags [100%]
- Read multiple bam files [100%]
- find non-overlappingpeaks [100%]
- get motif indices [100%]
- get kmer indices [100%]
- get gapped kmer indices [100%]

Computing variability
- compute variability for individual peaks (applicable for bulk samples) [90%] 
- compute variability for sets of peaks [100%]

Basic extensions for variability
- find samples with high or low accessibility for set of peaks [90%] 
- find peaks with high or low accessibility across set of samples [90%] 

More complex extensions
- find groups of annotations/motifs/kmers with independent variability [85%] 
- for groups of kmers, group kmers that likely make up the same motif [80%] 
- clustering of samples based on variability [70%]
- correlation with gene expression [0%]

Differential variablity
- find sets with greater variability in one sample versus another [0%]

Visualization
- plot variability [100%] 
- visualize kmers [90%] 
- heatmaps deviatons across samples [80%] 

Documentation 
- functions [70%] 
- classes [70%] 
- vignette [0%] 








