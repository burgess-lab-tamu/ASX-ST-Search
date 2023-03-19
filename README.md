# ASX-ST-Search
Scripts to search &amp; analyze ASX/ST motifs

## Usage

search ASX/ST motifs
```bash
1) open 'motif search' folder
2) download either search_ASX or search_ST.py
3) scripts will need pdb_id and chain_id as inputs.  
   two csv files of unique chains and their corresponding PDB ID can be obtained in the same page.
   chains included here are non-redundant chains filtered from all entries in PDB till 2020 Nov.
   PDB files corresponding to PDB IDs in the csv file must be downloaded to local disk: details of downloading PDB files can be checked in https://www.rcsb.org/downloads
4) scripts are ready to run with two csv files and local PDB files. 
5) output of search_ASX and search_ST gives csv files which include locations of ASX/ST motifs (by PDB ID, chain ID and starting residue ID (Asp/Asn/Ser/Thr)), classfications and at interfaces or not.
6) The output csv files from search_ASX.py or search_ST.py can be used as inputs in sequence_secondary_dihedral.py to aquire sequence information, secondary structures and dihedrals of N' - N4 in each ASX motifs or ST motifs. 
7) The output from sequence_secondary_dihedral.py serves as the datasets of unique ASX or ST motifs in PDB, which can be found in 'nonredundant-ASX/ST,zip'.
8) Following filter scripts can be applied on nonredundant ASX/ST datasets to build childern datasets.  An example of childern dataset, helical N-cap ASX motifs and ST-motifs, were searched by requring motifs in class 5 and at helical N-termini.
9) users can generate different children datasets based on their demands.

```


bioinformatics
```bash
general information
1) open 'bioinformaics' folder
2) in jupyter notebook, run ASX or ST_informatics using generated databases from 'motif search' folder. 
   statistical information such as number of different types of motifs, avearge/standard deviation of dihedral angles and residue abundance can be obtained as output of the last command in the notebook.

hydrophobic interactions
1) open 'count_hydrophobic_interaction' folder
2) run python files inside using generated databases from 'motif search' folder.
3) output is a csv file including locations of ASX motifs, [N',N3,N4], 3-element array [XXX] with 1 for 'have interactions' and 0 for 'none' in the orders of 'N'-N3','N'-N4' and 'N3-N4', and atoms involved in the interactions.

database
generated databases including hydrophobic patterns of ASX/ST motifs (with hydrophobic N' and N4) are provided in zip format.
```

## Citation

to be published
