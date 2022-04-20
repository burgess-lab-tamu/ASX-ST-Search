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


# In[2]:


def interfaceResidues(cmpx, cA='c. A', cB='c. B', cutoff=1.0, selName="interface"):
        """
        interfaceResidues -- finds 'interface' residues between two chains in a complex.

        PARAMS
            cmpx
                The complex containing cA and cB

            cA
                The first chain in which we search for residues at an interface
                with cB

            cB
                The second chain in which we search for residues at an interface
                with cA

            cutoff
                The difference in area OVER which residues are considered
                interface residues.  Residues whose dASA from the complex to
                a single chain is greater than this cutoff are kept.  Zero
                keeps all residues.

            selName
                The name of the selection to return.

        RETURNS
            * A selection of interface residues is created and named
                depending on what you passed into selName
            * An array of values is returned where each value is:
                ( modelName, residueNumber, dASA )

        NOTES
            If you have two chains that are not from the same PDB that you want
            to complex together, use the create command like:
                create myComplex, pdb1WithChainA or pdb2withChainX
            then pass myComplex to this script like:
                interfaceResidues myComlpex, c. A, c. X

            This script calculates the area of the complex as a whole.  Then,
            it separates the two chains that you pass in through the arguments
            cA and cB, alone.  Once it has this, it calculates the difference
            and any residues ABOVE the cutoff are called interface residues.

        AUTHOR:
            Jason Vertrees, 2009.		
        """
        # Save user's settings, before setting dot_solvent
        oldDS = cmd.get("dot_solvent")
        cmd.set("dot_solvent", 1)

        # set some string names for temporary objects/selections
        tempC, selName1 = "tempComplex", selName+"1"
        chA, chB = "chA", "chB"

        # operate on a new object & turn off the original
        cmd.create(tempC, cmpx)
        cmd.disable(cmpx)

        # remove cruft and inrrelevant chains
        cmd.remove(tempC + " and not (polymer and (%s or %s))" % (cA, cB))

        # get the area of the complete complex
        cmd.get_area(tempC, load_b=1)
        # copy the areas from the loaded b to the q, field.
        cmd.alter(tempC, 'q=b')

        # extract the two chains and calc. the new area
        # note: the q fields are copied to the new objects
        # chA and chB
        cmd.extract(chA, tempC + " and (" + cA + ")")
        cmd.extract(chB, tempC + " and (" + cB + ")")
        cmd.get_area(chA, load_b=1)
        cmd.get_area(chB, load_b=1)

        # update the chain-only objects w/the difference
        cmd.alter( "%s or %s" % (chA,chB), "b=b-q" )

        # The calculations are done.  Now, all we need to
        # do is to determine which residues are over the cutoff
        # and save them.
        stored.r, rVal, seen = [], [], []
        cmd.iterate('%s or %s' % (chA, chB), 'stored.r.append((model,resi,b))')

        cmd.enable(cmpx)
        cmd.select(selName1, 'none')
        for (model,resi,diff) in stored.r:
            key=resi+"-"+model
            if abs(diff)>=float(cutoff):
                if key in seen: continue
                else: seen.append(key)
                rVal.append( (model,resi,diff) )
                # expand the selection here; I chose to iterate over stored.r instead of
                # creating one large selection b/c if there are too many residues PyMOL
                # might crash on a very large selection.  This is pretty much guaranteed
                # not to kill PyMOL; but, it might take a little longer to run.
                cmd.select( selName1, selName1 + " or (%s and i. %s)" % (model,resi))

        # this is how you transfer a selection to another object.
        cmd.select(selName, cmpx + " in " + selName1)
        # clean up after ourselves
        cmd.delete(selName1)
        cmd.delete(chA)
        cmd.delete(chB)
        cmd.delete(tempC)
        # show the selection
        cmd.enable(selName)

        # reset users settings
        cmd.set("dot_solvent", oldDS)
        return rVal


# In[3]:


__main__.pymol_argv = [ 'pymol', '-qc']
ST_total = []
def read_PDB(PDB_path, PDB_csv, chain_csv, start, end, csv_name):
    time_stamp = datetime.datetime.now()
    print("time_stamp       " + time_stamp.strftime('%Y.%m.%d-%H:%M:%S'))
    PDB_list = pd.read_csv(PDB_csv)
    chain_list = pd.read_csv(chain_csv)
    pdb = PDB_list.iloc[:,1]
    chain = chain_list.iloc[:,1]
    fail_pdb = os.listdir('failed_PDB_pdb/')
    fail_pdb = [i.replace('.pdb','') for i in fail_pdb]
    print(len(fail_pdb))
    print(fail_pdb[5])
    start = int(start)
    end = int(end)
    if start == 0:
        pdb = pdb[:end]
        chain = chain[:end]
    else:
        pdb = list(pdb[start-1:end])
        chain = list(chain[start-1:end])
    #print(pdb[:20], chain[:20])
    length = end - start
    for i in range(length):
        cmd.select('all')
        cmd.delete('all')
        if i > 0 and i%50 == 0:
            df = pd.DataFrame(ST_total)
            df.to_csv(csv_name + '.csv')
        #print(pdb[i])
        if pdb[i] not in fail_pdb:
            new_path = PDB_path + '/' + str(pdb[i]) + '.pdb'
            print(new_path)
            if os.path.isfile(new_path):
                cmd.load(new_path,'test')
                cmd.remove('solvent')
                cmd.remove('ino.')
                cmd.remove('org.')
                chain_id = str(chain[i]).replace(' ','')
                chains = cmd.get_chains('test')
                atoms = cmd.get_model('chain ' + chain_id)
                resi_id = []
                for at in atoms.atom:
                    #print(at.resi)
                    if at.resi.isnumeric():
                        resi_id.append(int(at.resi))
                frequency = {x:resi_id.count(x) for x in resi_id}
                L_r, N_r = frequency.keys(), frequency.values()
                average = np.average(list(N_r))
                print(average)
                if average < 5.5:
                    cmd.select('all')
                    cmd.delete('all')
                    continue
                resi_id_final = list(dict.fromkeys(resi_id))
                resi_id_final.sort()
                chain_ST = []
                if len(resi_id_final) >= 8:
                    cmd.h_add(selection = 'test')
                    for k in range(len(resi_id_final)):
                        try:
                            for at in atoms.atom:
                                #print(resi_id_final[k], at.resi)
                                if at.resi == str(resi_id_final[k]):
                                    resi_name = at.resn
                                    break
                            if resi_name == 'SER':
                                #print(k)
                                x = resi_id_final[k]
                                #print(x)
                                if x < resi_id_final[-1] - 4 and x + 4 == resi_id_final[k+4]:
                                    #print('angles for' + str(resi_id_final[k]))
                                    distance_s2 = cmd.distance('chain ' + chain_id + ' & i. ' + str(x) + ' & n. OG', 'chain '+ chain_id +' & i. ' + str(x+2) + ' & n. N')
                                    angle_s2 = cmd.get_angle('chain ' + chain_id + ' & i. ' + str(x) + ' & n. OG', 'chain '+ chain_id +' & i. ' + str(x+2) + ' & h. within 2 of n. N', 'chain '+ chain_id +' & i. ' + str(x+2) + ' & n. N')
                                    distance_m3 = cmd.distance('chain ' + chain_id + ' & i. ' + str(x) + ' & n. O', 'chain ' + chain_id + ' & i. ' + str(x+3) + ' & n. N')
                                    angle_m3 = cmd.get_angle('chain ' + chain_id + ' & i. ' + str(x) + ' & n. O', 'chain ' + chain_id + ' & i. ' + str(x+3) + ' & h. within 2 of n. N', 'chain ' + chain_id + ' & i. ' + str(x+3) + ' & n. N')
                                    distance_s3 = cmd.distance('chain ' + chain_id + ' & i. ' + str(x) + ' & n. OG', 'chain ' + chain_id + ' & i. ' + str(x+3) + ' & n. N')
                                    angle_s3 = cmd.get_angle('chain ' + chain_id + ' & i. ' + str(x) + ' & n. OG', 'chain '+ chain_id +' & i. ' + str(x+3) + ' & h. within 2 of n. N', 'chain '+ chain_id +' & i. ' + str(x+3) + ' & n. N')
                                    distance_m4 = cmd.distance('chain ' + chain_id + ' & i. ' + str(x) + ' & n. O', 'chain ' + chain_id + ' & i. ' + str(x+4) + ' & n. N')
                                    angle_m4 = cmd.get_angle('chain ' + chain_id + ' & i. ' + str(x) + ' & n. O', 'chain '+ chain_id +' & i. ' + str(x+4) + ' & h. within 2 of n. N', 'chain '+ chain_id +' & i. ' + str(x+4) + ' & n. N')
                                    #print('ds2 {}, ds3 {}, dm3 {}, dm4 {}'.format(distance_s2, distance_s3, distance_m3, distance_m4))
                                    #print('angle_s2 {}, angle_s3 {}, angle_m3 {}, angle_m4 {}'.format(angle_s2, angle_s3, angle_m3, angle_m4))
                                    if (distance_s2 < 3.5 and angle_s2 >= 140) and (distance_s3 > 3.5 or angle_s3 < 140) and (distance_m3 < 3.5 and angle_m3 >= 140) and (distance_m4 > 3.5 or angle_m4 < 140):
                                        #print(distance_1, distance_2)
                                        #print(x)
                                        #print('common_shellman_loop')
                                        output_format = [pdb[i], chain_id, x, 'C1']
                                        chain_ST.append(output_format)
                                    elif (distance_s2 < 3.5 and angle_s2 >= 140) and (distance_s3 > 3.5 or angle_s3 < 140) and (distance_m3 < 3.5 and angle_m3 >= 140) and (distance_m4 < 3.5 and angle_m4 >= 140):
                                        #print(distance_1, distance_2)
                                        #print(x)
                                        #print('common_shellman_loop')
                                        output_format = [pdb[i], chain_id, x, 'C2']
                                        chain_ST.append(output_format)
                                    elif (distance_s2 < 3.5 and angle_s2 >= 140) and (distance_s3 < 3.5 and angle_s3 >= 140) and (distance_m3 < 3.5 and angle_m3 >= 140) and (distance_m4 > 3.5 or angle_m4 < 140):
                                        #print(distance_1, distance_2)
                                        #print(x)
                                        #print('common_shellman_loop')
                                        output_format = [pdb[i], chain_id, x, 'C2a']
                                        chain_ST.append(output_format)
                                    elif (distance_s2 < 3.5 and angle_s2 >= 140) and (distance_s3 < 3.5 and angle_s3 >= 140) and (distance_m3 < 3.5 and angle_m3 >= 140) and (distance_m4 < 3.5 and angle_m4 >= 140):
                                        #print(distance_1, distance_2)
                                        #print(x)
                                        #print('common_shellman_loop')
                                        output_format = [pdb[i], chain_id, x, 'C3']
                                        chain_ST.append(output_format)
                                    elif (distance_s2 > 3.5 or angle_s2 < 140) and (distance_s3 < 3.5 and angle_s3 >= 140) and (distance_m3 < 3.5 and angle_m3 >= 140) and (distance_m4 > 3.5 or angle_m4 < 140):
                                        #print(distance_1, distance_2)
                                        #print(x)
                                        #print('common_shellman_loop')
                                        output_format = [pdb[i], chain_id, x, 'C3a']
                                        chain_ST.append(output_format)
                                    elif (distance_s2 < 3.5 and angle_s2 >= 140) and (distance_s3 < 3.5 and angle_s3 >= 140) and (distance_m3 > 3.5 or angle_m3 < 140) and (distance_m4 < 3.5 and angle_m4 >= 140):
                                        #print(distance_1, distance_2)
                                        #print(x)
                                        #print('common_shellman_loop')
                                        output_format = [pdb[i], chain_id, x, 'C4']
                                        chain_ST.append(output_format)
                                    elif (distance_s2 > 3.5 or angle_s2 < 140) and (distance_s3 < 3.5 and angle_s3 >= 140) and (distance_m3 < 3.5 and angle_m3 >= 140) and (distance_m4 < 3.5 and angle_m4 >= 140):
                                        #print(distance_1, distance_2)
                                        #print(x)
                                        #print('common_shellman_loop')
                                        output_format = [pdb[i], chain_id, x, 'C4a']
                                        chain_ST.append(output_format)
                                    elif (distance_s2 > 3.5 or angle_s2 < 140) and (distance_s3 < 3.5 and angle_s3 >= 140) and (distance_m3 > 3.5 or angle_m3 < 140) and (distance_m4 < 3.5 and angle_m4 >= 140):
                                        #print(distance_1, distance_2)
                                        #print(x)
                                        #print('common_shellman_loop')
                                        output_format = [pdb[i], chain_id, x, 'C5']
                                        chain_ST.append(output_format)
                            if resi_name == 'THR':
                                #print(k)
                                x = resi_id_final[k]
                                #print(x)
                                if x < resi_id_final[-1] - 4 and x + 4 == resi_id_final[k+4]:
                                    #print('angles for' + str(resi_id_final[k]))
                                    distance_s2 = cmd.distance('chain ' + chain_id + ' & i. ' + str(x) + ' & n. OG1', 'chain '+ chain_id +' & i. ' + str(x+2) + ' & n. N')
                                    angle_s2 = cmd.get_angle('chain ' + chain_id + ' & i. ' + str(x) + ' & n. OG1', 'chain '+ chain_id +' & i. ' + str(x+2) + ' & h. within 2 of n. N', 'chain '+ chain_id +' & i. ' + str(x+2) + ' & n. N')
                                    distance_m3 = cmd.distance('chain ' + chain_id + ' & i. ' + str(x) + ' & n. O', 'chain ' + chain_id + ' & i. ' + str(x+3) + ' & n. N')
                                    angle_m3 = cmd.get_angle('chain ' + chain_id + ' & i. ' + str(x) + ' & n. O', 'chain ' + chain_id + ' & i. ' + str(x+3) + ' & h. within 2 of n. N', 'chain ' + chain_id + ' & i. ' + str(x+3) + ' & n. N')
                                    distance_s3 = cmd.distance('chain ' + chain_id + ' & i. ' + str(x) + ' & n. OG1', 'chain ' + chain_id + ' & i. ' + str(x+3) + ' & n. N')
                                    angle_s3 = cmd.get_angle('chain ' + chain_id + ' & i. ' + str(x) + ' & n. OG1', 'chain '+ chain_id +' & i. ' + str(x+3) + ' & h. within 2 of n. N', 'chain '+ chain_id +' & i. ' + str(x+3) + ' & n. N')
                                    distance_m4 = cmd.distance('chain ' + chain_id + ' & i. ' + str(x) + ' & n. O', 'chain ' + chain_id + ' & i. ' + str(x+4) + ' & n. N')
                                    angle_m4 = cmd.get_angle('chain ' + chain_id + ' & i. ' + str(x) + ' & n. O', 'chain '+ chain_id +' & i. ' + str(x+4) + ' & h. within 2 of n. N', 'chain '+ chain_id +' & i. ' + str(x+4) + ' & n. N')
                                    #print('ds2 {}, ds3 {}, dm3 {}, dm4 {}'.format(distance_s2, distance_s3, distance_m3, distance_m4))
                                    #print('angle_s2 {}, angle_s3 {}, angle_m3 {}, angle_m4 {}'.format(angle_s2, angle_s3, angle_m3, angle_m4))
                                    if (distance_s2 < 3.5 and angle_s2 >= 140) and (distance_s3 > 3.5 or angle_s3 < 140) and (distance_m3 < 3.5 and angle_m3 >= 140) and (distance_m4 > 3.5 or angle_m4 < 140):
                                        #print(distance_1, distance_2)
                                        #print(x)
                                        #print('common_shellman_loop')
                                        output_format = [pdb[i], chain_id, x, 'C1']
                                        chain_ST.append(output_format)
                                    elif (distance_s2 < 3.5 and angle_s2 >= 140) and (distance_s3 > 3.5 or angle_s3 < 140) and (distance_m3 < 3.5 and angle_m3 >= 140) and (distance_m4 < 3.5 and angle_m4 >= 140):
                                        #print(distance_1, distance_2)
                                        #print(x)
                                        #print('common_shellman_loop')
                                        output_format = [pdb[i], chain_id, x, 'C2']
                                        chain_ST.append(output_format)
                                    elif (distance_s2 < 3.5 and angle_s2 >= 140) and (distance_s3 < 3.5 and angle_s3 >= 140) and (distance_m3 < 3.5 and angle_m3 >= 140) and (distance_m4 > 3.5 or angle_m4 < 140):
                                        #print(distance_1, distance_2)
                                        #print(x)
                                        #print('common_shellman_loop')
                                        output_format = [pdb[i], chain_id, x, 'C2a']
                                        chain_ST.append(output_format)
                                    elif (distance_s2 < 3.5 and angle_s2 >= 140) and (distance_s3 < 3.5 and angle_s3 >= 140) and (distance_m3 < 3.5 and angle_m3 >= 140) and (distance_m4 < 3.5 and angle_m4 >= 140):
                                        #print(distance_1, distance_2)
                                        #print(x)
                                        #print('common_shellman_loop')
                                        output_format = [pdb[i], chain_id, x, 'C3']
                                        chain_ST.append(output_format)
                                    elif (distance_s2 > 3.5 or angle_s2 < 140) and (distance_s3 < 3.5 and angle_s3 >= 140) and (distance_m3 < 3.5 and angle_m3 >= 140) and (distance_m4 > 3.5 or angle_m4 < 140):
                                        #print(distance_1, distance_2)
                                        #print(x)
                                        #print('common_shellman_loop')
                                        output_format = [pdb[i], chain_id, x, 'C3a']
                                        chain_ST.append(output_format)
                                    elif (distance_s2 < 3.5 and angle_s2 >= 140) and (distance_s3 < 3.5 and angle_s3 >= 140) and (distance_m3 > 3.5 or angle_m3 < 140) and (distance_m4 < 3.5 and angle_m4 >= 140):
                                        #print(distance_1, distance_2)
                                        #print(x)
                                        #print('common_shellman_loop')
                                        output_format = [pdb[i], chain_id, x, 'C4']
                                        chain_ST.append(output_format)
                                    elif (distance_s2 > 3.5 or angle_s2 < 140) and (distance_s3 < 3.5 and angle_s3 >= 140) and (distance_m3 < 3.5 and angle_m3 >= 140) and (distance_m4 < 3.5 and angle_m4 >= 140):
                                        #print(distance_1, distance_2)
                                        #print(x)
                                        #print('common_shellman_loop')
                                        output_format = [pdb[i], chain_id, x, 'C4a']
                                        chain_ST.append(output_format)
                                    elif (distance_s2 > 3.5 or angle_s2 < 140) and (distance_s3 < 3.5 and angle_s3 >= 140) and (distance_m3 > 3.5 or angle_m3 < 140) and (distance_m4 < 3.5 and angle_m4 >= 140):
                                        #print(distance_1, distance_2)
                                        #print(x)
                                        #print('common_shellman_loop')
                                        output_format = [pdb[i], chain_id, x, 'C5']
                                        chain_ST.append(output_format)
                        except:
                            print('error')
                            pass
                    print(chain_id, pdb[i])
                    if chain_ST != []:
                        chain_interface_residues = []
                        for l in range(len(chains)):
                            if chains[l] != chain_id:
                                #print(chain_id[l])
                                interface = interfaceResidues('test', cA='c. ' + chain_id, cB='c. ' + chains[l])
                                n_interface = np.shape(interface)[0]
                                for m in range(n_interface):
                                    if interface[m][0] == 'chA':
                                        if interface[m][1].isnumeric():
                                            chain_interface_residues.append(int(interface[m][1]))
                        chain_interface_residues.sort()
                        chain_interface_residues = list(dict.fromkeys(chain_interface_residues))
                        #print(chain_interface_residues)
                        for n in range(len(chain_ST)):
                            #print(n)
                            for o in range(len(chain_interface_residues)):
                                #print(o)
                                if chain_ST[n][3] == 'C1' or chain_ST[n][3] == 'C2a' or chain_ST[n][3] == 'C3a':
                                    if chain_interface_residues[o] >= int(chain_ST[n][2]) and chain_interface_residues[o] <= (int(chain_ST[n][2]) + 3):
                                        chain_ST[n].append('surface')
                                        break
                                else:
                                    if chain_interface_residues[o] >= int(chain_ST[n][2]) and chain_interface_residues[o] <= (int(chain_ST[n][2]) + 4):
                                        chain_ST[n].append('surface')
                                        break
                        print('goood!')
                        ST_total.extend(chain_ST)
                cmd.select('all')
                cmd.delete('all')
    time_stamp = datetime.datetime.now()
    df = pd.DataFrame(ST_total)
    df.to_csv(csv_name + '.csv')
    print("time_stamp       " + time_stamp.strftime('%Y.%m.%d-%H:%M:%S'))
    return ST_total


# In[4]:


#test = read_PDB('PDB_10000_final', 'nonredundant_PDBID.csv', 'nonredundant_CHAINID.csv', 0, 10, 'ST_test')


# In[ ]:


if __name__== "__main__":
    read_PDB(sys.argv[1],sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6])

