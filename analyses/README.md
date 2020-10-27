### Analyses overview

The evaluation comprises the following analyses involving a specific tool or group of methods and one of the two collections of data sets:

- `dada2_callahan`: DADA2 on Callahan data sets
- `dada2_franzen`: DADA2 on Franzén data sets
- `model_supported_callahan`: model-supported clustering and refinement methods of GeFaST on Callahan data sets
- `model_supported_franzen`: model-supported clustering and refinement methods of GeFaST on Franzén data sets
- `quality_weighted_callahan`: quality-weighted clustering methods of GeFaST on Callahan data sets
- `quality_weighted_franzen`: quality-weighted clustering methods of GeFaST on Franzén data sets
- `uvsearch_callahan`: USEARCH and VSEARCH on Callahan data sets
- `uvsearch_franzen`: USEARCH and VSEARCH on Franzén data sets 

Each such analysis subfolder contains the folder structure and initial files to repeat the respective analysis.
The correponding workflow is described in the eponymous notebook in the subfolder, which also shows a first evaluation of the clustering results.

Above analyses are summarised as follows:

- `model_supported_overall`: The model-supported clustering and refinement methods of GeFaST are compared to DADA2, USEARCH and VSEARCH 
    on both the Callahan and Franzén data sets in `model_supported_overall.ipynb`. 
- `quality_weighted_overall`: The quality-weighted clustering methods of GeFaST are compared to USEARCH and VSEARCH on both the Callahan and Franzén data sets.
    Initially, the different quality-weighted cost functions are evaluated separately in both quality alignment-score and quality Levenshtein mode in order to pick the best variant of each cost function.
    The corresponding notebooks are stored in the `quality_weighted_overall/` subfolder. 
    Subsequently, the picked variants are compared to the other tools in `quality_weighted_overall.ipynb`.
