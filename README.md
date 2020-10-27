# Evaluation of quality-aware clustering and refinement implemented in GeFaST

This repository covers the comparative evaluation of the quality-aware clustering and refinement methods implemented in GeFaST
and presented in *On the use of sequence-quality information in OTU clustering* (submitted).
GeFaST implements two groups of quality-aware methods.
The *quality-weighted* methods involve alignments with quality-weighted cost functions,
while the *model-supported* methods are inspired by DADA2's consistency model.
The clustering quality of GeFaST is compared to DADA2, USEARCH and VSEARCH on two collections of data sets described in
[*DADA2: High-resolution sample inference from Illumina amplicon data*](https://doi.org/10.1038/nmeth.3869) (Callahan et al.) and
[*Improved OTU-picking using long-read 16S rRNA gene amplicon sequencing and generic hierarchical clustering*](https://doi.org/10.1186/s40168-015-0105-6) (Franzén et al.).

## How to use it

The evaluation is two-staged.
First, the different methods and tools were evaluated separately on the two collections of data sets.
Then, the results were aggregated and the clustering quality of the tools was compared over all data sets.
The workflow of the different parts of evaluation is described in the notebooks contained in the `analyses/` folder. 

The required software (see below) is not installed automatically. 
The analyses expect them two be installed in the `tools/` folder but the paths can be adjusted.

There are several starting points from which the evaluation can be repeated:
1) The code in the `scripts/` folder and the files provided in the subfolders of `analyses/` 
allow to repeat the analyses from scratch by following the workflow described in the respective notebook.
2) The data preparation can be skipped by downloading the respective input files from [PUB](https://doi.org/10.4119/unibi/2946192)
and ignoring the corresponding workflow commands in the notebooks. 
3) The actual evaluation of the clustering results can be rerun by using the aggregated information available as archived CSV files
from [PUB](https://doi.org/10.4119/unibi/2946192).

## Required software
 * [GeFaST](https://github.com/romueller/gefast) (version 2.0.0)
 * [USEARCH](http://www.drive5.com/usearch/download.html) (version 11.0.667)
 * [VSEARCH](https://github.com/torognes/vsearch) (version 2.14.2)
 * [ART](https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm) (`art_illumina`, VanillaIceCream release)
 * [Infernal](eddylab.org/infernal/) (`cmbuild` and `cmalign`, version 1.1.2)
 * GCC (version 4.9.2 or higher)
 * python (version 2.7 or higher; with Biopython, NumPy, matplotlib, pandas, seaborn, statistics)
 * R (version 3.5.3 or higher; with packages dada2, ggplot2, ShortRead)
 * Jupyter Notebook (for viewing or rerunning the evaluations)
 * common Unix tools (gzip, sed, tar, wget)


_Notes:_
Older versions of GCC, python, R etc. might also work but have not been tested.
Version 11 of USEARCH was used to evaluate its clustering quality. To stay as close as possible to the workflow of Callahan et al.,
the preparation of the Callahan data sets used the versions 8.0.1623 and 10.0.240.
