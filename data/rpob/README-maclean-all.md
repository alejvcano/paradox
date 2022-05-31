README for rpoB data from MacLean, et al. in maclean-all.csv

# Source 

MacLean RC, Perron GG, Gardner A. 2010. Diminishing returns from beneficial mutations and pervasive epistasis shape the fitness landscape for rifampicin resistance in Pseudomonas aeruginosa. Genetics 186:1345-1354.

# Processing

This table combines data on frequency and selection coefficient from Table 1 of MacLean, et al. with the data on mutation from Table 2. 

The relative fitness assigned to a variant is taken to be the average of any measurements available in Table 1 for the 3 different backgrounds A455G, A455T, and C1550T. 

The frequency at which a variant is found in evolution, per Table 1, is defined as the combined occurrence in the 3 backgrounds, where the occurrence is not precluded by the nature of the background genotype (the "NP" values in the original Table 1), and taking into account that the reported frequencies in Table 1 have different denominators.  That is, though 96 trials were conducted in each background, this resulted in different numbers of identifiable rpoB mutations. The discretization of frequency values in Table 1 indicates that the denominators for the frequencies are 70, 78 and 78 in the backgrounds of A455G, A455T, and C1550T, respectvely. 

The mutation value assigned for a variant is the number of counts given in Table 2 (out of a survey of 80 mutations).  We assign a count of 0 to all variants not listed in Table 2, given that the authors state "Any beneficial mutations not shown in this table were not found in this collection of mutants."  

 
