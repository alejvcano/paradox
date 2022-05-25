README for amidase data from Wrenbeck, et al. in wrenbeck2017_ace.csv

# Source 

Wrenbeck EE, Azouz LR, Whitehead TA. 2017. Single-mutation fitness landscapes for an enzyme on multiple substrates reveal specificity is globally encoded. Nat Commun 8:15695.

The supplementary data actually are on FigShare here: https://figshare.com/articles/dataset/Normalized_fitness_values_for_AmiE_selections/3505901

# Processing

wrenbeck2017_ace.csv is processed from the supplementary data for acetamide selection (file "amiESelectionFitnessData_Acetamide.txt").  This is a well-formed table with 21 rows for each position, covering all 20 amino acids plus stop.  The fitness values are coded as 0 for cases where the mutant amino acid is the wild-type.  These can be used to assign the wild-type amino acid.  I did a sanity check to confirm that the locations are all numbers from 1 to n where n is the number of unique locations. 

From the starting file for acetamide 
* assign the wild-type amino acid for each position
* convert from log2 values using fitness = 2^normalized_fitness
* add titv where this can be inferred (in case we need it)
* output with standardized column names in csv format 

