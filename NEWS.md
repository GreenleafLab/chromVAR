# News

## 0.4.0

### S4 Methods

Many functions have been converted to S4 methods to allow for different input formats.

### Motif matching

The motif matching code has been moved to a separate package, motifmatchr.  That code relies on the MOODSR C++ library, and this split was done for ease of licensing.  Additionally, the motif matching functionality could be useful outside of the context of chromVAR

### chromVARDeviations class

A new class has been introduced, which inherits from SummarizedExperiment, to hold output of compute_deviations. This class has a couple of accessor functions (see below).  This addition was mostly for increased transparency.

### Accessor methods

Several accessor methods have been added:
- counts: to access the counts assays of SummarizedExperiment object
- deviations: to access the bias corrected deviations of chromVARDeviations
- deviation_scores: to access the deviation Z-scores of chromVARDeviations