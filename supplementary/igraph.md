# IGraph Setup

# **Input Data**

The input data to create the graph features is done by the rinerator program by executing the following script:

> python2.7 ./RINerator_V0.5.1/Source/get_chains.py data/7KVE_relaxed.pdb data/igraph/in/ data/igraph/auxiliar/chains.txt

The manual setup of the rinerator program is available [here](../supplementary/rinerator_manual_installation.md). Besides that, we made a docker configuration with rinerator already configured.

The outputs are saved in the following path: `igraph/in`. The files necessary for the next step are:

- The ```.sif``` file is a graph delimeted by ```:```. Where in each line there is the type and subtype of the interaction, for instance ```cnt:mc_sc``` means *main chain*, *side chain*.

- The ```_intsc.ea``` file contains a score for every interactions.

- The  ```_nrint.ea``` file contains the number of interactions.



# **How to execute the code**

Inside the **src folder** execute:
```shell
python produce_graph_features.py
```



# **Output Data**

The previous command saves the following files in [this folder: data/igraph/out](`data/igraph/out`):

- 7KVE_relaxed_graph_data_unidirectional.csv
- 7KVE_relaxed_graph_features_unidirectional.csv



