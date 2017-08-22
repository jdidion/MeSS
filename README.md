R functions for computing Methylation Specificity Scores (MeSSes) from WGBS data in multiple tissues.

# Usage

```R
source('entropy.R')
specificity <- cpg_specificity(pos, meth, max.cores=8)
```

Where `pos` is a [GRanges](https://bioconductor.org/packages/release/bioc/html/GenomicRanges.html) object with genomic positions of `N` CpGs, and `meth` is an `N x M` matrix where `M` is the number of tissues and values are Beta estimates extracted from WGBS data.

# Details

Simplification and speed-up of SMART-BS-Seq by Liu et al. 
(https://pypi.python.org/pypi/SMART-BS-Seq).

Rather than starting from multiple .wig files, this script starts from a
matrix of CpG x sample. The output is a data.frame with the following columns:

* A specificity value [0-1] for each CpG, which Liu interpret as follows:
< 0.5 == low specificity; 0.5-0.75 == intermediate specificity;
> 0.75 == high specificity.
* Median methylation, computed by one-step Tukey biweight.
* Distance: geometric difference from previous row).
* DiffEntropy: entropy of difference from previous row).
* MeSS.*: Methylation Specificity Score for each tissue.
* Optionally (see below):
  * p.*: Significance of deviation from median, for each tissue.
  * padj.*: Adjusted p-value for each tissue.
  * ranks: MeSSes converted to ranks.

Minimum beta range threshold: When there are many beta values very close to
each other (this happens frequently at the extremes), but one or a few betas
are slightly different, the standard entropy calculation might yield a high
specificity. While technically correct, these scores cause difficulties for
interpretation. Liu et al. got around this problem by enforcing a minimum
difference of 1.0 (1%) between each beta and the median (biweight). Note that
this is not mentioned in their paper; one must read their code to find it out.
It is difficult to known why they chose 1.0, but we suspect it was an
arbitrary decision. This approach leads to a strange distribution of MeSSes,
with a large spike toward the low end (the value of the spike is determined
by the sample size; for 21 samples it was between 0.18-0.19). We instead
recommend determining this value empirically from the distribution of per-CpG
beta ranges, or by the resolution of the assay (1 / average read depth for
WGBS). We chose to go with the first option and selected 3.0 (3%) as our
threshold; YMMV.

P-values: To determine the specificity of methylation in each tissue,
Liu et al. use the per-sample weights from the Tukey Biweight calculation to
identify samples that belong to the "background" distribution and then run
a one-sample t-test for each non-background sample. They do not perform
multiple test correction. We assert that p-values are problematic for
determining tissue specificity for the reasons described above. Nevertheless,
we show that, after removing CpGs with beta value range below the empirically
determined threshold, there is a high degree of correlation between MeSS and
p-value. We also show that computing p-values based on betas is problematic
due to the non-normality of the data, and we show that first converting betas
to M values leads to higher correlation with MeSS. We also recommend computing
p-values using a 'leave-one-out' strategy (i.e. treat one sample as the mean
and use the remaining samples as the distribution). We provide an option
(disabled by default) to compute p-values.
