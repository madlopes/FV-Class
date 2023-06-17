# RINerator

The RINerator software can be downloaded [here](https://rinalyzer.de/rinerator.php). The version used is 0.5.1, released on December 2014.

**Execute the following steps in this order:**

- Download the **probe** and **reduce** folder (located in the folder called artifacts) and place them in```~/.local/share/```. Both should be visible by your system therefore one solution is to create a symbolic link to the probe and reduce files in the ``` /usr/bin``` folder.
    - To create a symbolic link: ```ln -s source_file myfile```. Replace source_file with the name of the existing file for which you want to create the symbolic link.

- Install python 2.7.18

- Make sure pip, setuptools and wheel are installed in the latest version available: ```pip install -U pip setuptools wheel```

- Install numpy lib (version 1.16.6): ```pip install numpy==1.16.6```

- Install biopython lib (version 1.59): ```pip install biopython==1.59```

- Create a file which contains the chains IDs, for instance ```chains.txt``` 

    - The file chains.txt contains the chains name of the structure separated by a comma.
    - In order to find the chains IDs go to the chimera software and click in Select $\rightarrow$ Chains $\rightarrow$ Chains IDs.
    
    
- Run the following command: ```python ./RINerator_V0.5.1/Source/get_chains.py <structure.pdb> <output_directory> <chains.txt>```
    - Make sure you are running the correct python version with: ```python --version```.
    - Make sure the dependencies are correct installed with: ```pip list```.

**Output:**

Some files will be generated in the folder previously provided. 

- The ```.sif``` file is a graph delimeted by ```:```. Where in each line there is the type and subtype of the interaction, for instance ```cnt:mc_sc``` means *main chain*, *side chain*.

- The ```_intsc.ea``` file contains a score for every interactions.

- The  ```_nrint.ea``` file contains the number of interactions.


# RODRIGO

Steps:

3 - python2.7 ./RINerator_V0.5.1/Source/get_chains.py data/7KVE_relaxed.pdb data/igraph/in/ data/igraph/auxiliar/chains.txt