### Evaluation code

The scripts implement the evaluation and are grouped according to their purpose:
- `configuration_handling.py` and `util.py` contain some utility code reused in multiple scripts.
- `analyses/`: The scripts execute the workflows on the two collections of data sets.
The individual steps, e.g. preparating the data or running a specific tool on it, are accessed through different commands.
- `cluster_evaluation/`: Code for assessing the clustering quality from the clustering outputs.
- `data_preparation/`: Computation of ground truths and the reformatting of input files.
- `tool_scripts/`: Wrapper code for conveniently running the examined tools with varying parameters. 