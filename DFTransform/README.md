#### README BioScript

I will try to compile here python implementations of known transform functions. I just started with TSS and CSS, maybe some will follow. Inputs are pandas dataframe.

##### Warning

These implementations have been done on my spare time and while I do think they are correct, I do not garantee they are providing the exact similar result compared to their original implementation. Please let me know if you spot aby errors which could be fixed.

##### Concerning CSS

This is a direct translation of Cumulative Sum Scalling (CSS) transform originally defined in [this paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4010126/) and proposed in the R package [metagenomeSeq](https://www.rdocumentation.org/packages/metagenomeSeq/versions/1.14.0). This implementation has been done with numpy. This is a direct translation of [cumNormStat](https://www.rdocumentation.org/packages/metagenomeSeq/versions/1.14.0/topics/cumNormStat) process, and not [cumNormStatFast](https://www.rdocumentation.org/packages/metagenomeSeq/versions/1.14.0/topics/cumNormStatFast), with its associated [calcNormFactors](https://github.com/HCBravoLab/metagenomeSeq/blob/df8a28214fa9cb25870dee0e5cc909c160ce8da2/R/cumNorm.R#L42) function.