README for rpoB data from MacLean, et al. in maclean-all.csv

# Source 

MacLean RC, Perron GG, Gardner A. 2010. Diminishing returns from beneficial mutations and pervasive epistasis shape the fitness landscape for rifampicin resistance in Pseudomonas aeruginosa. Genetics 186:1345-1354.

# Context for use

MacLean, et al (2010) carried out selection for Rifampicin resistance in replicate cultures of Pseudomonas aeruginosa, with 96 replicates each for 3 slightly different genetic backgrounds.  Their purpose was to study epistasis and diminishing returns.  However, this study was quite unique in that the authors assessed mutation rates. They simply isolated 80 Rifampicin-resistant colonies with rpoB mutations, reporting the count for each type of mutation recovered.  Out of 36 mutations identified in the selection experiment, 11 were identified in the mutation assay, with numbers ranging from 1 to 30.  That is, the observed mutation rates ranged 30-fold, and the range must be much larger if one considers the 25 unobserved mutations. 

Therefore, this is a rare study that reports a triplet of frequency evolved, selection coefficient, and mutation rate for 11 or more variants. 

# Processing

The table included here, derived from MacLean, et al., combines data on frequency and selection coefficient from Table 1 of MacLean, et al. with the data on mutation from Table 2. 

The relative fitness assigned to a variant is taken to be the average of any measurements available in Table 1 for the 3 different backgrounds A455G, A455T, and C1550T. 

The frequency at which a variant is found in evolution, per Table 1, is defined as the combined occurrence in the 3 backgrounds, where the occurrence is not precluded by the nature of the background genotype (the "NP" values in the original Table 1), and taking into account that the reported frequencies in Table 1 have different denominators.  That is, though 96 trials were conducted in each background, this resulted in different numbers of identifiable rpoB mutations. The discretization of frequency values in Table 1 indicates that the denominators for the frequencies are 70, 78 and 78 in the backgrounds of A455G, A455T, and C1550T, respectvely. 

The mutation value assigned for a variant is the number of counts given in Table 2 (out of a survey of 80 mutations).  We assign a count of 0 to all variants not listed in Table 2, given that the authors state "Any beneficial mutations not shown in this table were not found in this collection of mutants."  

 
