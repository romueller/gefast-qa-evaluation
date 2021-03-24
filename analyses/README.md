### Analyses overview

The evaluation comprises the analyses listed below, involving a specific tool or group of methods and one of the two collections of data sets.

*Clustering quality:*
- `dada2_callahan`: DADA2 on Callahan data sets
- `dada2_franzen`: DADA2 on Franzén data sets
- `model_supported_callahan`: model-supported clustering and refinement methods of GeFaST on Callahan data sets
- `model_supported_franzen`: model-supported clustering and refinement methods of GeFaST on Franzén data sets
- `quality_weighted_callahan`: quality-weighted clustering methods of GeFaST on Callahan data sets
- `quality_weighted_franzen`: quality-weighted clustering methods of GeFaST on Franzén data sets
- `swarm_callahan`: Swarm on Callahan data sets
- `swarm_franzen`: Swarm on Franzén data sets
- `uvsearch_callahan`: USEARCH, VSEARCH and UPARSE on Callahan data sets
- `uvsearch_franzen`: USEARCH, VSEARCH and UPARSE on Franzén data sets

*Performance:*
- `performance`: runtime and memory consumption of GeFaST and the other tools on the largest data set of the quality evaluation

Each such analysis subfolder contains the folder structure and initial files to repeat the respective analysis.
The correponding workflow is described in the eponymous notebook in the subfolder, which also shows a first evaluation of the clustering results.

Above analyses are summarised as follows:

- `model_supported_overall`: The model-supported clustering and refinement methods of GeFaST are compared to Swarm, DADA2, USEARCH, VSEARCH and UPARSE 
    on both the Callahan and Franzén data sets in `model_supported_overall.ipynb`. 
- `quality_weighted_overall`: The quality-weighted clustering methods of GeFaST are compared to Swarm, USEARCH, VSEARCH and UPARSE on both the Callahan and Franzén data sets.
    Initially, the different quality-weighted cost functions are evaluated separately in both quality alignment-score and quality Levenshtein mode in order to pick the best variant of each cost function.
    The corresponding notebooks are stored in the `quality_weighted_overall/` subfolder. 
    Subsequently, the picked variants are compared to the other tools in `quality_weighted_overall.ipynb`.
