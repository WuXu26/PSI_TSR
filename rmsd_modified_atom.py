#program to calculate the similarity between structures using rmsd

import csv
import math
import Bio.PDB
from Bio.PDB import PDBParser
import pandas as pd
import os
from scipy.stats import skew
import numpy as np
from itertools import combinations

df=pd.read_csv('sample_details_et_mix5_name_id_atom.csv') #Change here for different sample details file
PDB_list = df['protein'].to_list()
Chain = df['chain'].to_list()
Drug_name = df['drug_name'].to_list()
Drug_id = df['drug_id'].to_list()
Group = df['group'].to_list()
Atom_list=df['ligand_atom'].to_list()

def calDist(x1, y1, z1, x2, y2, z2): #Calculate distance between two atoms
    return math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2 + (z2 - z1) ** 2)

def Cal_Centroid(res): #Calculates the centriod of structures
    SumX = 0
    SumY = 0
    SumZ = 0
    counter = 0
    for atom1 in residue:
        atomCoord = atom1.get_vector()
        SumX += atomCoord[0]
        SumY += atomCoord[1]
        SumZ += atomCoord[2]
        counter += 1
    Centroid_X_ = SumX / counter
    Centroid_Y_ = SumY / counter
    Centroid_Z_ = SumZ / counter
    return Centroid_X_,Centroid_Y_,Centroid_Z_

def translation(*args): #translates the two structures
    CenX=args[0]
    CenY=args[1]
    CenZ=args[2]
    residue_coord_dict={}
    counter=0
    for atom in residue:
        atomCoord = atom.get_vector()
        X1_New=CenX-atomCoord[0]
        Y1_New=CenY-atomCoord[1]
        Z1_New =CenZ-atomCoord[2]
        residue_coord_dict[counter]=[X1_New,Y1_New,Z1_New]
        counter+=1
    return residue_coord_dict
Dataframe_Index=[]

def kabsch_umeyama(A, B): #the algoritms for alignment by roration and scaling
    assert A.shape == B.shape
    n, m = A.shape

    EA = np.mean(A, axis=0)
    EB = np.mean(B, axis=0)
    VarA = np.mean(np.linalg.norm(A - EA, axis=1) ** 2)

    H = ((A - EA).T @ (B - EB)) / n
    U, D, VT = np.linalg.svd(H)
    d = np.sign(np.linalg.det(U) * np.linalg.det(VT))
    S = np.diag([1] * (m - 1) + [d])

    R = U @ S @ VT
    c = VarA / np.trace(np.diag(D) @ S)
    t = EA - c * R @ EB

    return R, c, t

def rmsd(A,B): #The rmsd method
    sum=0
    for i in range(len(A)):
        dist=calDist(A[i][0],A[i][1],A[i][2],B[i][0],B[i][1],B[i][2])
        sum+=dist**2
    return np.sqrt(sum/len(A))




counter=0
drug_coord_dict={}
Output_Folder_Path="/ddnB/work/wxx6941/TSR/code/code/psi_revision/5pdb/output_chl_atom_rmsd"
OutputFile=open(f'{Output_Folder_Path}/similarity_drugs_rmsd.txt','w')
for i in range(len(PDB_list)):
    PDB_ID=PDB_list[i]
    Chain_Name = Chain[i]
    drug_Name = Drug_name[i]
    drug_Id = Drug_id[i]
    atoms=Atom_list[i]
    if drug_Name == 'CL0':
        #drug_Name='CL0'
        drug_Id=f'0{drug_Id}'

    Drug=f'{PDB_ID}_{Chain_Name}_{drug_Name}_{drug_Id}'
    Dataframe_Index.append(Drug)
    #print(PDB_ID,Chain_Name,drug_Name,drug_Id)
    PDB_File_Path = "/ddnB/work/wxx6941/TSR/code/code/psi_revision/5pdb/PDB_datadir_HRemoved/{}.pdb".format(
        PDB_ID)
    p = Bio.PDB.PDBParser()
    Structure = p.get_structure('PrimaryStructureChain', PDB_File_Path)
    model = Structure[0]
    for chain in model:
        if chain.id == Chain_Name:
            for residue in chain:
                drug_Identifier = str(residue)[17:18].strip()
                resName = str(residue)[9:12].strip()  # residue Name

                numeric_filter = filter(str.isdigit, str(residue.id))
                Res_Id = "".join(numeric_filter)  # Residue ID

                if drug_Identifier == 'H' and resName == drug_Name and Res_Id == str(drug_Id):
                    drug_coord = []
                    for atom in residue:
                        if atom.get_name() in atoms:
                            atomCoord = atom.get_vector()
                            drug_coord.append([atomCoord[0], atomCoord[1], atomCoord[2]])


                    drug_coord_dict[Drug] = drug_coord
                    counter+= 1
                    #Centroid_X,Centroid_Y,Centroid_Z=Cal_Centroid(residue)
                    #print(Centroid_X,Centroid_Y,Centroid_Z)
                    #print(translation(Centroid_X,Centroid_Y,Centroid_Z))

comb_drug = combinations([*drug_coord_dict], 2)
for i in comb_drug:
    A=np.array(drug_coord_dict[i[0]])
    B=np.array(drug_coord_dict[i[1]])
    R, c, t = kabsch_umeyama(A, B)

    B = np.array([t + c * R @ b for b in B])
    rmsd_similarity=rmsd(A,B)
    print(rmsd_similarity)
    OutputFile.write(f'{i[0]}\t{i[1]}\t{rmsd_similarity}')
    OutputFile.write('\n')

#df = pd.DataFrame(Total_drug_dict_features, index=Dataframe_Index)
#df = df.astype('int')
#df.to_csv(f"feature_map_with_header.csv", header=True, index=True)
OutputFile.close()
print('completed successfully')





