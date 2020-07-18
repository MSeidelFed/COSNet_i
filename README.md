# IntCryOmics
Integration of omics relative changes into Cryo-EM based randomly sampled interaction networks of multiprotein complexes.
The project is composed of independent components:

## Introduction
IntCryOmics is a collection of scripts to randomly select spatial neighborhoods of proteins from a multiprotein complex in order to test whether these neighbors characterize a region within the complex that becomes significantly enriched upon any experimental procedure. The procedure has been detailed in XXXXX and XXXXX publication.

*Workflow*

![Workflow](images/intcryomics_workflow.png)

## Usage

*Describes the workflow from data generation to output*

Python version

```
module add devel/Python-3.6.5
```

## Get PDB files

Use linux terminal to build PDBs from mmCIF files

e.g.:

```
## make directory to store PDB objects

mkdir trial

## extract PDB objects from the mmCIF entity
### python3 [path_to_function] [path_to_mmCIF] [output_path]

python3 Python_Modules/split_cif_by_entity.py Data/4v7e.cif trial/
```

now the PDB objects should be in the trial directory

## Generate FASTA files (optional)

This module is optional within the entire workflow and is meant as a quality control step in order to align the monomer sequences of the constituent molecular entities in the analyzed multimeric complex.

e.g.:

```
## grabbing entity number from original mmCIF

grep "L30" Data/4v7e.cif

## grabbing the sequences and generating the FASTA
### python3 [path_to_function] [path_to_PDB] [path_to_mmCIF] [entity_number] [entity_type] [output_path]

python3 Python_Modules/extract_sequence.py trial/4v7e_L30_60S_ribosomal.pdb Data/4v7e.cif 69 protein .
```

The function will output the sequences of the mmCIF object according to the entity number (69 in the example) and type (protein in the example), it will only consider the sequences a match if every single monomer is listed in both strings in which case no further message will be displayed, otherwise a warning message will tell the user that the sequences are not equal and a visual verification becomes necessary.

## Reindex PDBs (optional)

The residues column inside PDBs can be reindexed if there are inconsistencies in numbering to prevent holes in structures

e.g.:

```
## make new directory to store reindexed PDBs

mkdir trial_redxd

## reindexing objects with a bash script
### names file, if you create the names file with excel you might want to getrid of the carriage returns, e.g.:

tr -d '\r' < file_names.txt

### bash Batch_files/batch_reindex_PDBs.sh [pathtopythoncode] [pathto/infile] [pathtopdbs] [outpath] [startingreindexnumber]

bash Batch_files/batch_reindex_PDBs.sh Python_Modules/ Data/file_names.txt trial/ trial_redxd/ 1


## reindexing one by one
### python3 [path_to_function] [starting_reindex_number] [path_to_PDB] [output_path]

python3 Python_Modules/reindex_pdb.py 1 trial/4v7e_L30_60S_ribosomal.pdb trial_redxd/uL30_rpL7.pdb
```

After iterating through every object the output, reindexed PDBs should be in the trial_redxd folder. At this point you may select which object to include to fit the distance network. In our exemplary case we only used the protein PDBs to fit the network, hiding the remaining PDB files in a separate folder.


## Fit the distance network

Generate file with the combinations of names that will be used to calculate distance matrices between entities

e.g.:

```
## make new text file with the names of reindexed objects that will be used for the network

ls trial_redxd > infile.txt

## Generate file with the combinations of names that will be used to calculate distance matrices between entities
### python3 combination.py [filewithlistoffiles] [combi_num] [outfilename] [outpath]

python3 Python_Modules/combination.py infile.txt 2 combi_names.txt .
```

calculate a distance matrix between entities using reindexed PDBs as input, the file with names conminations as future relationships in the matrix and a bash script to iterate the process between each pair of PDB files. Alternatively one can do each pair individually without using the bash script

e.g.:

```
## make new directory to store distance matrices

mkdir trial_DistMat

## Calculating distance matrices between objects with a bash script
### bash Batch_files/batch_calc_dist.sh [pathtopythoncode] [pathto/infile] [pathtopdbs] [outpath]

bash Batch_files/batch_calc_dist.sh Python_Modules/ combi_names.txt trial_redxd/ trial_DistMat/

## Calculating distance matrices one by one
### python3 calculate_distance.py [pathtoPDBfile1] [pathtoPDBfile2] [outputpath]

python3 Python_Modules/calculate_distance.py trial_redxd/eL13_RPL13.pdb trial_redxd/eL14_RPL14.pdb .

```

Please note that the function expects a naming scheme of the form XX_XX.pdb, which should be arranged in the file_names.txt so that the reindexed objects align with such a scheme. The function then grabs the characters before the floor dash as ID of the resulting .csv matrix.

The resulting distance matrices can be imported into R for visualization.


## Treshold selection and transit probability

select a threshold to define a contact and build a contact matrix based on the distances. Plus output the number of atoms that are in contact to later on define the percentage of coverage and transit probability.

e.g.:

```
## make new directory to store contact matrices

mkdir trial_contacts

## make a list with the generated csv matrices

ls trial_DistMat/ > csv_names_file.txt

## calculate contact number among entities according to a defined treshold
### python3 contacts_from_dist.py [threshold] [pathtoinput] [pathtocsvs] [pathtopdbs] [pathtooutput]

python3 Python_Modules/contacts_from_dist.py 8 csv_names_file.txt trial_DistMat/ trial_redxd/ trial_contacts/
```

The function will take first text protion from the proposed nameing scheme (XX_XX) and map the contacts using the reindexed PDBs, whenever it encounters a redundance in terms of naming, the user will need to input manually the number of the evaluated RP given a list of options. After iteration trough all the files in the csv_names_file.txt several outputs are produced:

Output 1: summary contacts =   eL37 eL39 10 514 161 First two columns determine whether a contact exists. third is the number of residues in contact, fourth is the number of conections between atoms, and last the number of atoms conected.

Output 2: contacts_t8_eL37_eL39_full.dat = 17 51 7.129363 49 18 First is residue # from first protein, second is residue # from second protein and third is distance between them, fourth number of conections between atoms and fifth number of atoms conected. The file named "contacts_t8_eL37_eL39.dat" does not contain the last two columns

The first two columns of Output 1 are used to fit the network, and in order to prepare the edgelist for the random walk, one must also create an edges_with_weights.txt file, containing the numbers of amino acids in contact as well as the nodes in contact. This is simply:

```
awk '{print $1" "$2" "$3}' trial_contacts/summary_contacts_t8.txt > edges_with_weights.txt
```

## Random Walk and Fisher exact test

The user must input a significances file for the hypergeometric or Fisher exact test with the same protein identifiers as those in the edge list and binary calls where 1 is significant for condition X and 0 is not significant, both columns must be separated by a space " ". The module allows repeated protein identifiers, or in other words, multiple significances per identifier in separate lines. This feature is meant to cope with paralogs within protein families, which are characteristic of multiprotein complexes.

e.g.:

```
## python3 intcryomics.py [filewithedgelist] [sigfile] [walklength] [iterationnum]
python3 Python_Modules/intcryomics.py edges_with_weights.txt Data/significance_file 20 10
```

### List of nodes

['eL13', 'eL15', 'eL18', 'uL15', 'uL1', 'uL29', 'uL4', 'eL14', 'eL20', 'eL6', 'uL13', 'uL6', 'eL36', 'eL42', 'eL8', 'uL2', 'uL30', 'eL19', 'eS7', 'uS17', 'eL21', 'uL16', 'eL29', 'uL18', 'eL24', 'eS8', 'uL14', 'uL3', 'eL27', 'eL30', 'eL34', 'eL28', 'eL32', 'eL43', 'eS1', 'uS15', 'eL33', 'eL39', 'eL37', 'uL24', 'uL23', 'eL40', 'uL5', 'eS10', 'eS12', 'uS14', 'uS3', 'eS31', 'eS17', 'RACK1', 'uS2', 'eS19', 'uS13', 'uS9', 'eS26', 'uS11', 'eS21', 'eS27', 'uS4', 'uS5', 'uS8', 'eS24', 'eS4', 'eS6', 'eS25', 'uS7', 'eS28', 'eS30', 'uS12', 'P1', 'P2', 'uL10', 'uL11', 'uS10', 'uS19']

### Minumum covering set

i.e., smallest set that spans the entire node space, gives preference to larger subsets

[{64, 65, 1, 5, 40, 42, 74, 13, 46, 14, 48, 49, 51, 52, 53, 54, 23, 55},
{2, 36, 6, 7, 8, 9, 10, 11, 16, 20, 23, 24, 27},
{35, 68, 60, 46, 48, 50, 18, 56, 57, 59, 28, 29, 30},
{0, 1, 34, 33, 4, 5, 37, 40, 12, 13, 14, 15, 28},
{67, 68, 50, 19, 56, 58, 59, 60, 61, 62, 63},
{73, 43, 44, 45, 46, 47, 48},
{72, 69, 70, 71},
{36, 6, 7, 8, 9, 10, 16, 24, 25, 26, 27, 31},
{0, 33, 1, 3, 37, 38, 39, 40, 5, 6, 14, 15, 16, 28, 29, 30},
{2, 6, 7, 8, 9, 16, 20, 21, 22, 23},
{0, 2, 36, 6, 7, 8, 9, 10, 11, 41, 27, 31},
{35, 68, 48, 17, 18, 50, 19, 56, 57, 59, 60},
{32, 2, 36, 6, 7, 9, 10, 16, 27, 31},
{65, 34, 33, 66, 1, 14, 15, 54, 55, 29}]

### Sampled Regions

i.e., most visited nodes during the random walks starting from the first node

Region 0:	['eS25', 'uS7', 'eL15', 'uL29', 'uL23', 'uL5', 'uS19', 'eL42', 'uS3', 'eL8', 'eS17', 'RACK1', 'eS19', 'uS13', 'uS9', 'eS26', 'uL18', 'uS11']

Region 1:	['eL18', 'eL33', 'uL4', 'eL14', 'eL20', 'eL6', 'uL13', 'uL6', 'uL30', 'eL21', 'uL18', 'eL24', 'uL3']

Region 2:	['uS15', 'uS12', 'uS8', 'uS3', 'eS17', 'uS2', 'eS7', 'eS21', 'eS27', 'uS5', 'eL27', 'eL30', 'eL34']

Region 3:	['eL13', 'eL15', 'eS1', 'eL43', 'uL1', 'uL29', 'eL39', 'uL23', 'eL36', 'eL42', 'eL8', 'uL2', 'eL27']

Region 4:	['eS30', 'uS12', 'uS2', 'uS17', 'eS21', 'uS4', 'uS5', 'uS8', 'eS24', 'eS4', 'eS6']

Region 5:	['uS10', 'eS10', 'eS12', 'uS14', 'uS3', 'eS31', 'eS17']

Region 6:	['uL11', 'P1', 'P2', 'uL10']

Region 7:	['eL33', 'uL4', 'eL14', 'eL20', 'eL6', 'uL13', 'uL30', 'eL24', 'eS8', 'uL14', 'uL3', 'eL28']

Region 8:	['eL13', 'eL43', 'eL15', 'uL15', 'eL39', 'eL37', 'uL24', 'uL23', 'uL29', 'uL4', 'eL8', 'uL2', 'uL30', 'eL27', 'eL30', 'eL34']

Region 9:	['eL18', 'uL4', 'eL14', 'eL20', 'eL6', 'uL30', 'eL21', 'uL16', 'eL29', 'uL18']

Region 10:	['eL13', 'eL18', 'eL33', 'uL4', 'eL14', 'eL20', 'eL6', 'uL13', 'uL6', 'eL40', 'uL3', 'eL28']

Region 11:	['uS15', 'uS12', 'eS17', 'eL19', 'eS7', 'uS2', 'uS17', 'eS21', 'eS27', 'uS5', 'uS8']

Region 12:	['eL32', 'eL18', 'eL33', 'uL4', 'eL14', 'eL6', 'uL13', 'uL30', 'uL3', 'eL28']

Region 13:	['uS7', 'eS1', 'eL43', 'eS28', 'eL15', 'eL8', 'uL2', 'eS26', 'uS11', 'eL30']

### Significances

#### Total number significant proteins: 

36

#### Significant regions Fisher exact test: 

[0.6710013944369708, 1.0, 0.0011576953008030787, 0.12538968815729165, 0.18519788445879923, 0.5219422565692198, 0.0002274210985913904, 0.6162482642592604, 1.0, 0.4110255105920462, 0.31178055617954936, 0.003713939573632823, 0.7839471731468854, 0.139941724115829]

#### Significant regions after Bonferroni correction: 

1. Logical array ([False, False,  True, False, False, False,  True, False, False, False, False, False, False, False])

2. Significances array ([1.        , 1.        , 0.01620773, 1.        , 1.        , 1.        , 0.0031839 , 1.        , 1.        , 1.        ,  1.        , 0.05199515, 1.        , 1.        ])


