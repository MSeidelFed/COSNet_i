# IntCryOmics USAGE.md
This is the usage file for the integration of omics relative changes into Cryo-EM based randomly sampled interaction networks of multiprotein complexes. The project is composed of independent components, written as python scripts (found in [Python_Modules](https://github.com/MSeidelFed/IntCryOmics/tree/master/Python_Modules)), which can be run in batch with bash scripts (found in [Batch_files](https://github.com/MSeidelFed/IntCryOmics/tree/master/Batch_files)). 

Here we document the usage of the various modules, with an example using the yeast and rabbit ribosomal protein complexes, with rRNA removed (PDB ID: [6snt](https://www.rcsb.org/structure/6SNT)) (PDB ID: [6gz5](https://www.rcsb.org/structure/6GZ5)).

## Usage

**In the sections below, we describe the workflow from data generation to output. Code should be run from the command line, Linux-style.**
**Individual python scripts can all be run with the -h or --help flag to get the usage.**

For example, running ```python3 reindex_pdb.py -h``` in the command line gives:
```sh
usage: python3 reindex_pdb.py [-h] <startidx> <pdbfile> <outfile>

Returns PDB file with residues reindexed according to start index

positional arguments:
  startidx    Starting index, usually 1
  pdbfile     PDB file to be reindexed
  outfile     Reindexed PDB file

optional arguments:
  -h, --help  show this help message and exit
```

### _Dependencies_
All code for this project is written for Python version 3 and above. Running on an older version of Python will not work! If you are working with a modular server system, don't forget to load the correct module.
```bash
$ module add devel/Python-3.6.5
```

Otherwise make sure you have Python 3 installed. Additionally, we use the Biopython PDB package to handle PDB files. Installation info can be found [here](https://biopython.org/wiki/Download).

## Extract PDB files from mmCIF files

Many multiprotein complexes are given in mmCIF-formatted files. Here is how to use the command line (linux-like) to extract PDBs from mmCIF files, which we save into the trial directory.
_e.g.:_

**1. Make directory to store PDB objects**
```bash
$ mkdir trial
```
**2. Extract PDB objects from the mmCIF entity**

_Usage:_ ```python3 <path_to_function> <path_to_mmCIF> <output_path>```

```bash
$ python3 Python_Modules/split_cif_by_entity.py Data/4v7e.cif trial/
```

The PDB objects will be saved as individual files in the `trial` directory.
## Generate FASTA files (optional)

This module is optional within the entire workflow and is meant as a quality control step in order to align the monomer sequences of the constituent molecular entities in the analyzed multimeric complex.
_e.g.:_

**1. Grabbing entity number from original mmCIF**
```bash
$ grep "L30" Data/4v7e.cif
```
**2. Grabbing the sequences and generating the FASTA**

_Usage:_ ```python3 <path_to_function> <path_to_PDB> <path_to_mmCIF> <entity_number> <entity_type> <output_path>```
```bash
$ python3 Python_Modules/extract_sequence.py trial/4v7e_L30_60S_ribosomal.pdb Data/4v7e.cif 69 protein .
```

The function will output the sequences of the mmCIF object according to the entity number (69 in the example) and type (protein in the example), it will only consider the sequences a match if every single monomer is listed in both strings in which case no further message will be displayed, otherwise a warning message will tell the user that the sequences are not equal and a visual verification becomes necessary.

## Reindex PDBs (optional)

The residues column inside PDBs can be reindexed (sequential incrementation) if there are inconsistencies in numbering to prevent holes in structures, or if for any other reason the user wishes to start from a certain index.
_e.g.:_

**1. make new directory to store reindexed PDBs**
```bash
$ mkdir trial_redxd
```
**2a. BATCH version: reindexing objects with a bash script**

- Make a 2-column names file, if you create the names file with excel you might want to getrid of the carriage returns, e.g.:
```bash
$ tr -d '\r' < file_names.txt
```
- The file needs to have, per line, the input pdb file, followed by a whitespace and the output file name: ```file.pdb file_reindex.pdb```.

_Usage:_ ```bash Batch_files/batch_reindex_PDBs.sh <pathtopythoncode> <pathto/infile> <pathtopdbs> <outpath> <startingreindexnumber>```

- Make sure your input paths to directories include the final slash `/`

```bash
$ bash Batch_files/batch_reindex_PDBs.sh Python_Modules/ Data/file_names.txt trial/ trial_redxd/ 1
```
**2b. SINGLE version: reindexing one by one**

_Usage:_ ```python3 reindex_pdb.py [-h] <startidx> <pdbfile> <outfile>```
```bash
$ python3 Python_Modules/reindex_pdb.py 1 trial/4v7e_L30_60S_ribosomal.pdb trial_redxd/uL30_rpL7.pdb
```

After iterating through every object the output, reindexed PDBs should be in the trial_redxd folder. At this point you may select which object to include to fit the distance network. In our exemplary case we only used the protein PDBs to fit the network, hiding the remaining PDB files in a separate folder.

## Fit the distance network

Generate file with the combinations of names that will be used to calculate distance matrices between entities
_e.g.:_

**1. Make new text file with the names of reindexed objects that will be used for the network**
```bash
$ ls trial_redxd > infile.txt
```
**2. Generate file with the combinations of names that will be used to calculate distance matrices between entities**

_Usage:_ ```python3 combination.py [-h] <inputfile> <N> <outfile> <outpath```
```bash
$ python3 Python_Modules/combination.py infile.txt 2 combi_names.txt .
```
The N here is your combination base value, so use 2. 

Calculate a distance matrix between entities using reindexed PDBs as input, the file with names conminations as future relationships in the matrix and a bash script to iterate the process between each pair of PDB files. Alternatively one can do each pair individually without using the shell script to run in batch.
_e.g.:_
**3. Make new directory to store distance matrices**
```bash
$ mkdir trial_DistMat
```
**4. Calculating distance matrices between objects with a bash script**

_Usage:_ ```bash Batch_files/batch_calc_dist.sh <pathtopythoncode> <pathto/infile> <pathtopdbs> <outpath>```
```bash
$ bash Batch_files/batch_calc_dist.sh Python_Modules/ combi_names.txt trial_redxd/ trial_DistMat/
```
In this case your ```combi_names.txt``` file should be a two column file where you list ```file1 file2``` with a whitespace inbetween.

**5. Calculating distance matrices one by one**

_Usage:_ ```python3 calculate_distance.py [-h] <pathto/PDBfile1> <pathto/PDBfile2> <outpath>```
```bash
$ python3 Python_Modules/calculate_distance.py trial_redxd/eL13_RPL13.pdb trial_redxd/eL14_RPL14.pdb .
```
Please note that the function expects a naming scheme of the form XX_XX.pdb, which should be arranged in the file_names.txt so that the reindexed objects align with such a scheme. The function then grabs the characters before the underscore as the ID of the resulting .csv matrix. 
The resulting distance matrices can be imported into R for visualization.


## Threshold selection and transit probability

Select a threshold to define a contact and build a contact matrix based on the distances. Plus output the number of atoms that are in contact to later on define the percentage of coverage and transit probability.
_e.g.:_

**1. Make new directory to store contact matrices**
```bash
$ mkdir trial_contacts
```
**2. Make a list with the generated csv matrices**
```bash
$ ls trial_DistMat/ > csv_names_file.txt
```
**3. Calculate contact number among entities according to a defined threshold**

_Usage:_ ```python3 contacts_from_dist.py <threshold> <pathtoinput> <pathtocsvs> <pathtopdbs> <pathtooutput>```
```bash
$ python3 Python_Modules/contacts_from_dist.py 8 csv_names_file.txt trial_DistMat/ trial_redxd/ trial_contacts/
```

The function will take the first text portion from the proposed nameing scheme (XX_XX) and map the contacts using the reindexed PDBs, whenever it encounters a redundance in terms of naming, the user will need to input manually the number of the evaluated RP given a list of options. After iteration through all the files in the csv_names_file.txt several outputs are produced:

Output 1: summary contacts =   eL37 eL39 10 514 161 First two columns determine whether a contact exists. third is the number of residues in contact, fourth is the number of conections between atoms, and last the number of atoms conected.

Output 2: contacts_t8_eL37_eL39_full.dat = 17 51 7.129363 49 18 First is residue # from first protein, second is residue # from second protein and third is distance between them, fourth number of conections between atoms and fifth number of atoms conected. The file named "contacts_t8_eL37_eL39.dat" does not contain the last two columns

The first two columns of Output 1 are used to fit the network, and in order to prepare the edgelist for the random walk, one must also create an edges_with_weights.txt file, containing the numbers of amino acids in contact as well as the nodes in contact. This is simply:

```
$ awk '{print $1" "$2" "$3}' trial_contacts/summary_contacts_t8.txt > Results/edges_with_weights_4v7e_t8.txt
```

## Random Walk and Fisher Exact Test

The user must input a significances file for the hypergeometric or Fisher exact test with the same protein identifiers as those in the edge list and binary calls where 1 is significant for condition X and 0 is not significant, both columns must be separated by a space " ". The module allows repeated protein identifiers, or in other words, multiple significances per identifier in separate lines. This feature is meant to cope with paralogs within protein families, which are characteristic of multiprotein complexes. The example output is listed below:
_e.g.:_

_Usage:_ ```python3 intcryomics.py <filewithedgelist> <sigfile> <walklength> <iterationnum>```
```bash
$ python3 Python_Modules/intcryomics.py edges_with_weights.txt Data/significance_file 20 10 > Results/IntCryOmics_4v7e_t8.txt
```

### List of nodes

['eL13', 'eL15', 'eL18', 'uL15', 'uL1', 'uL29', 'uL4', 'eL14', 'eL20', 'eL6', 'uL13', 'uL6', 'eL36', 'eL42', 'eL8', 'uL2', 'uL30', 'eL19', 'eS7', 'uS17', 'eL21', 'uL16', 'eL29', 'uL18', 'eL24', 'eS8', 'uL14', 'uL3', 'eL27', 'eL30', 'eL34', 'eL28', 'eL32', 'eL43', 'eS1', 'uS15', 'eL33', 'eL39', 'eL37', 'uL24', 'uL23', 'eL40', 'uL5', 'eS10', 'eS12', 'uS14', 'uS3', 'eS31', 'eS17', 'RACK1', 'uS2', 'eS19', 'uS13', 'uS9', 'eS26', 'uS11', 'eS21', 'eS27', 'uS4', 'uS5', 'uS8', 'eS24', 'eS4', 'eS6', 'eS25', 'uS7', 'eS28', 'eS30', 'uS12', 'P1', 'P2', 'uL10', 'uL11', 'uS10', 'uS19']
### Minumum covering set
_i.e., smallest set that spans the entire node space, gives preference to larger subsets_

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

_i.e., most visited nodes during the random walks starting from the first node_

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

### Significances:

#### Total number significant proteins: 

36

#### Significant regions Fisher exact test: 
_given per region, in order_

[0.6710013944369708, 1.0, 0.0011576953008030787, 0.12538968815729165, 0.18519788445879923, 0.5219422565692198, 0.0002274210985913904, 0.6162482642592604, 1.0, 0.4110255105920462, 0.31178055617954936, 0.003713939573632823, 0.7839471731468854, 0.139941724115829]

#### Significant regions after Bonferroni correction: 

1. Logical array ([False, False,  True, False, False, False,  True, False, False, False, False, False, False, False])

2. Significances array ([1.        , 1.        , 0.01620773, 1.        , 1.        , 1.        , 0.0031839 , 1.        , 1.        , 1.        ,  1.        , 0.05199515, 1.        , 1.        ])


## Network drawing and highlight of specific regions

In order to color nodes, the following script allows users to map the subunit belonging of each node in the source and target columns, plus an optional argument is to select all the nodes that interact after the random walk with an specific protein of interest.
_e.g.:_

_Usage:_ ```python3 pimp_my_network.py <Names_file> [Intcryomics_file] [protein_ID] ```
```bash

## names file
ls Results/edges* > Network_names_file.txt

## IntCryOmics names file
ls Results/IntCryOmics_* > IntCryOmics_names_file.txt

$ python3 Python_Modules/pimp_my_network.py Network_names_file.txt IntCryOmics_names_file.txt eL39
```

If an IntCryOmics names file and a protein identifier are not given, the nodes will be named based on the subunit they belong to, either SSU or LSU. This procedure generates as outcome text files with a network structure that can be then visualized. The following is an example that was visualized in [Cytoscape](https://cytoscape.org/) using the protein eL39 as an example to follow the polypeptide exit tunnel (PET) region at different distance thresholds:


_t_ = 5Å ![4v7e_t5](images/out2edges_with_weights_4v7e_t5.txt.png)
_t_ = 8Å ![4v7e_t8](images/out2edges_with_weights_4v7e_t8.txt.png)
_t_ = 12Å ![4v7e_t12](images/out2edges_with_weights_4v7e_t12.txt.png)
_t_ = 20Å ![4v7e_t20](images/out2edges_with_weights_4v7e_t20.txt.png)

