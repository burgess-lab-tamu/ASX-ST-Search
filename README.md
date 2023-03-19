# ASX-ST-Search
Scripts to search &amp; analyze ASX/ST motifs

## Usage

search ASX/ST motifs
```bash
1) open 'motif search' folder
2) download either serach_ASX or search_ST.py
3) scripts will need pdb_id and chain_id as inputs.  
   two csv files of unique chains and their corresponding PDB ID can be obtained in the same page.
   chains included here are non-redundant chains filtered from all entries in PDB till 2020 Nov.
   PDB files corresponding to PDB IDs in the csv file must be downloaded to local disk: details of downloading PDB files can be checked in https://www.rcsb.org/downloads
4) scripts are ready to run with two csv files and local PDB files. 
5) output is a csv file including the location of ASX/ST motifs (by PDB ID, chain ID and starting residue ID (Asp/Asn/Ser/Thr)) 

database
generated databases of unqiue ASX/ST motifs are provided in zip format.
filters can be implemented on the family database to generate children database.
a children database of alpha-helical ASX/ST motifs is also provided.
```


bioinformatics
```bash
general information
1) open 'bioinformaics' folder
2) in jupyter notebook, run ASX or ST_informatics using generated databases from 'motif search' folder. 
   statistical information such as number of different types of motifs, avearge/standard deviation of dihedral angles and residue abundance can be obtained.

hydrophobic interactions
1) open 'count_hydrophobic_interaction' folder
2) run python files inside using generated databases from 'motif search' folder.

database
generated databases including hydrophobic patterns of ASX/ST motifs (with hydrophobic N' and N4) are provided in zip format.
```

## Citation

to be published
