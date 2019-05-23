# REvariation
A repository for storing materials of Zou and Zhang (2019, in review) study on 
Relative Exchangeability (RE) variation among 90 clades of species on the Tree of Life.

01_sequences:

This directory contains the pairwise codon sequence alignments of all 90 clades. 
The FASTA-format files are named with the prefix refering to the clade names shown 
in Table S1 of the paper.

02_scripts:

This directory contains the relevant scripts mentioned in Materials and Methods section
of the paper.

01_concat_clean.py: An example script for concatenating coding sequences of individual
                    genes and filtering out codon sites with gaps, missing or ambiguous
                    data.

02_get_matrix.py: An example script for deriving the transition matrix used in sequence
                  simulations.
                  
03_sim_seq.py: An example script for simulating sequence evolution in a clade.

CodonZ.py, TsTv.py: Python module file
