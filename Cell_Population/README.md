# Cell Population calculation 

## Notes
This standalone module is compatible with Python 3 (PhylogicNDT is only compatible with Python 2.7)

This module allows you to select the most likely tree and calcualtes the cell population based on the selected tree. (Current PhylogicNDT output calculates cell abundance based on the top tree, which may not always be the case.)

## Input files 
The -sif flag should provide a sif file

The -m flag should provide a mutation ccf file from the clustering result

The -c flag should provide a cluster ccf file from the clustering result

The --tree number should provide the most likely tree number from PhylogicNDT output 

## Output files
Three output files should be generated:

Cell_population_abundances.tsv

Cell_population_mcmc_trace.tsv

Constrained_ccf.tsv

## Run
	python ./PhylogicNDT.py CellPopulation -i Indiv_ID -sif Patient.sif  -m mutation_ccf_file -c cluster_ccf_file --tree_number correct_tree_number
