
[![DOI](https://zenodo.org/badge/327256687.svg)](https://zenodo.org/badge/latestdoi/327256687)

# MMEJ_detection


This Matlab code was designed to identify deletions with microhomology-based signatures by analyzing the DNA sequences flanking the deletions. 

The code is written in matlab/R2019b

## Source publication


Tzah Feldman, Akhiad Bercovich, Yoni Moskovitz, Noa Chapal-Ilani, Amanda Mitchell, … & Liran I. Shlush. (2021). Recurrent deletions in clonal hematopoiesis are driven by microhomology-mediated end joining. 

See: https://www.nature.com/articles/s41467-021-22803-y

## Running the code

Main(input_file_name,output_file_name, WITH_MISMATCH):

The input_file_name is a file containing details related to genomic sequences flanking the deletions.  More information of how to generate this input file is described in the source publication under the methods section. 

The following columns are required:

- Column 1 – the name of the cohort from which the deletion was extracted
- Column 2  - gene name
- Column 3 – chromosome
- Column 4 – deletion start position
- Column 5 – deletion end position
- Column 6–  genomic position 19 bp before the start position
- Column 7 – genomic position 20 bp after the end position
- Column 8 – Full genomic coordinates of the flanking sequence (starting from the position at column 6 and ending at the position at column 7)
- Column 9 – Genomic sequence of the flanking region

## Input_file example


| cohort | Gene | chr | start | end | startMinus19 | endPlus20 | fasta_coordinates | fasta | 
| ------ | ------ | ------ | ------ | ------ | ------ | ------ | ------ | ------ |
| BeatAML | SDHA | 5 | 218484 | 218484 | 218464 | 218504 | >5:218464-218504 | GCAGACATGTCGGGGGTCCGGGGCCTGTCGCGGCTGCTGA |
| BeatAML | NLRP6 | 11 | 281553 | 281555 | 281533 | 281575 | >11:281533-281575 | CGAGGACACCGAAGAGCCAGAGGAGGAGGAGGAGGGAGAGGA |
| MDS_cohort | PIGG | 4 | 517287 | 517289 | 517267 | 517309 | >4:517267-517309 | AAGGTGGTCAGAGCTAGACCTTCTTATTCTGTTGGGGACGGC |
| MDS_cohort | C11orf35 | 11 | 559313 | 559313 | 559293 | 559333 | >11:559293-559333 | CCCTCTTTCTCTCTCTTGTCCCCCCAGCTCAGGGCCCTAC |


The output_file_name is a csv file containing the input columns together with the following output columns: 
MH - the sequence of the microhomology (if detected)
Mismatch amount - the number of mismatches that were identified inside the microhomology
MH_len - the length of the microhomology (if detected)
Del_len - the deletion’s length


WITH_MISMATCH flag indicates whether to allow single bp mismatch inside the microhomology. Values are 0 (default) for not allowing and 1 for allowing mismatch.  




