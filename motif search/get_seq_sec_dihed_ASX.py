#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
import pymol
from pymol import cmd, stored
import numpy as np
import __main__
import pandas as pd
import sys
import datetime
from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP
import math
import collections
import shutil


# In[2]:


__main__.pymol_argv = [ 'pymol', '-qc']
information = []
def get_sequence_secondary_stucture(PDB_path, csv_path, output):
    schellman_csv = pd.read_csv(csv_path + '.csv')
    #get sequence
    pdb_id = list(schellman_csv.iloc[:,1])
    chain_id = list(schellman_csv.iloc[:,2])
    start_id = list(schellman_csv.iloc[:,3])
    loop_type = list(schellman_csv.iloc[:,4])
    surface = list(schellman_csv.iloc[:,5])
    for i in range(len(pdb_id)):
        if i > 0 and i%50 == 0:
            df = pd.DataFrame(information)
            df.to_csv(output + '.csv')
        print(pdb_id[i], chain_id[i], start_id[i], loop_type[i])
        cmd.select('all')
        cmd.delete('all')
        pdb_path = PDB_path + '/' + pdb_id[i] + '.pdb'
        if os.path.isfile(pdb_path):
            cmd.load(PDB_path + '/' + pdb_id[i] + '.pdb', 'test')
            atoms = cmd.get_model('chain ' + chain_id[i])
            sequence = []
            secondary_structure = []
            sub_information = []
            for j in range(8):
                for at in atoms.atom:
                    #print(at.resi)
                    if at.resi.isnumeric():
                        if int(at.resi) == start_id[i]-2+j:
                            #print(at.resn)
                            sequence.append(at.resn)
                            break
            for l in range(start_id[i]-1, start_id[i]+5):
                try:
                    p = PDBParser()
                    structure = p.get_structure("try", PDB_path + '/' + pdb_id[i] + '.pdb') #PDB path
                    model = structure[0]
                    dssp = DSSP(model, PDB_path + '/' + pdb_id[i] + '.pdb', dssp='mkdssp') # PDB path
                    id_x = (chain_id[i], (' ', l, ' '))
                    print(id_x)
                    sec = dssp[id_x][2]
                    print(sec)
                    secondary_structure.append(sec)
                except:
                    secondary_structure.append('Missing')
            secondary_structure = list(dict.fromkeys(secondary_structure))
            try:
                phi_psis = cmd.phi_psi('test & c. ' + chain_id[i] + ' & i. ' + str(start_id[i]-1) + '-' + str(start_id[i] + 5))
                phi_psis = list(phi_psis.items())
                phi_psi = []
                phi_psi.append(phi_psis)
            except:
                phi_psi = ['error']
            sub_information.append([pdb_id[i], chain_id[i], start_id[i], loop_type[i], surface[i], sequence, secondary_structure, phi_psi])
            #print(sub_information)
            information.extend(sub_information)
    data = pd.DataFrame(information)
    data.to_csv(output + '.csv')


# In[ ]:


if __name__== "__main__":
    get_sequence_secondary_stucture(sys.argv[1],sys.argv[2], sys.argv[3])

