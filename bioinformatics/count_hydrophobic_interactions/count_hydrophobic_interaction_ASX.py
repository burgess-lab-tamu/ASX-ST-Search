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


# In[9]:


def read_hydrophobic(csv, name):
    data = pd.read_csv(csv)
    pdb = list(data.iloc[:,1])
    chain = list(data.iloc[:,2])
    start = list(data.iloc[:,3])
    seq = list(data.iloc[:,6])
    final = []
    for i in range(len(pdb)):
        cmd.reinitialize()
        inter = []
        sub_seq = [seq[i][2:5], seq[i][-19:-16], seq[i][-12:-9]]
        pocket = []
        #print(sub_seq)
        mol = cmd.load('PDB_10000_final/' + pdb[i] + '.pdb', 'test')
        cmd.remove('hydrogen')
        #1-6
        cmd.select('sele, test & c. ' + chain[i] + ' & i. ' + str(start[i]-1) + ' & e. C &! n. C+CA+N+O' + ' within 4.5 of test & c. ' + chain[i] + ' & i. ' + str(start[i] + 4) + ' & e. C &! n. C+CA+N+O', state=1)
        sele_list = cmd.index('sele')
        pocket.append(sele_list)
        #print(sele_list)
        if sele_list != []:
            inter.append(1)
        else:
            inter.append(0)
        cmd.deselect()
        #1-5
        cmd.reinitialize()
        mol = cmd.load('PDB_10000_final/' + pdb[i] + '.pdb', 'test')
        cmd.remove('hydrogen')
        cmd.select('sele, test & c. ' + chain[i] + ' & i. ' + str(start[i]-1) + ' & e. C &! n. C+CA+N+O' + ' within 4.5 of test & c. ' + chain[i] + ' & i. ' + str(start[i] + 3) + ' & e. C &! n. C+CA+N+O', state=1)
        sele_list = cmd.index('sele')
        pocket.append(sele_list)
        #print(sele_list)
        if sele_list != []:
            inter.append(1)
        else:
            inter.append(0)
        cmd.deselect()
        #5-6
        cmd.reinitialize()
        mol = cmd.load('PDB_10000_final/' + pdb[i] + '.pdb', 'test')
        cmd.remove('hydrogen')
        cmd.select('sele, test & c. ' + chain[i] + ' & i. ' + str(start[i]+3) + ' & e. C &! n. C+CA+N+O' + ' within 4.5 of test & c. ' + chain[i] + ' & i. ' + str(start[i] + 4) + ' & e. C &! n. C+CA+N+O', state=1)
        sele_list = cmd.index('sele')
        pocket.append(sele_list)
        #print(sele_list)
        if sele_list != []:
            inter.append(1)
        else:
            inter.append(0)
        #print(inter)
        final.append([sub_seq,inter, pocket])
        if i > 0 and i%50 == 0:
            df = pd.DataFrame(final)
            df.to_csv(name + '.csv')
    final = pd.DataFrame(final)
    final.to_csv(name + '.csv')


# In[136]:


if __name__== "__main__":
    read_hydrophobic(sys.argv[1], sys.argv[2])

