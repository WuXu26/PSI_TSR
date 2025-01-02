# Program to calculate the amino acid triplets and key Frequency
# Author:Tarikul Islam Milon
#Created on: 06/27/2023


import csv
import math
import Bio.PDB
from Bio.PDB import PDBParser
import pandas as pd
import os

dTheta = 29
dLen = 35
numOfLabels = 20


df = pd.read_csv('sample_details_psi_ab_mix1.csv')
PDB_list = df['protein'].to_list()
Chain = df['chain'].to_list()
Group_Infor = df['group'].to_list()

atomSeq = {}
atomSeqNumber = open("aminoAcidCode_lexicographic_new.txt", 'r')
lines = atomSeqNumber.readlines()
for line in lines:
    aa_Res = line.split()[0]
    aa_Res_No =line.split()[1]
    atomSeq[aa_Res]=aa_Res_No
atomSeqNumber.close()


# Theta Bin for 3D
def thetaClass_(Theta):
    # classT=0
    if Theta >= 0 and Theta < 12.11:
        classT = 1
    elif Theta >= 12.11 and Theta < 17.32:
        classT = 2
    elif Theta >= 17.32 and Theta < 21.53:
        classT = 3
    elif Theta >= 21.53 and Theta < 25.21:
        classT = 4
    elif Theta >= 25.21 and Theta < 28.54:
        classT = 5
    elif Theta >= 28.54 and Theta < 31.64:
        classT = 6
    elif Theta >= 31.64 and Theta < 34.55:
        classT = 7
    elif Theta >= 34.55 and Theta < 37.34:
        classT = 8
    elif Theta >= 37.34 and Theta < 40.03:
        classT = 9
    elif Theta >= 40.03 and Theta < 42.64:
        classT = 10
    elif Theta >= 42.64 and Theta < 45.17:
        classT = 11
    elif Theta >= 45.17 and Theta < 47.64:
        classT = 12
    elif Theta >= 47.64 and Theta < 50.05:
        classT = 13
    elif Theta >= 50.05 and Theta < 52.43:
        classT = 14
    elif Theta >= 52.43 and Theta < 54.77:
        classT = 15
    elif Theta >= 54.77 and Theta < 57.08:
        classT = 16
    elif Theta >= 57.08 and Theta < 59.38:
        classT = 17
    elif Theta >= 59.38 and Theta < 61.64:
        classT = 18
    elif Theta >= 61.64 and Theta < 63.87:
        classT = 19
    elif Theta >= 63.87 and Theta < 66.09:
        classT = 20
    elif Theta >= 66.09 and Theta < 68.30:
        classT = 21
    elif Theta >= 68.30 and Theta < 70.5:
        classT = 22
    elif Theta >= 70.5 and Theta < 72.69:
        classT = 23
    elif Theta >= 72.69 and Theta < 79.2:
        classT = 24
    elif Theta >= 79.2 and Theta < 81.36:
        classT = 25
    elif Theta >= 81.36 and Theta < 83.51:
        classT = 26
    elif Theta >= 83.51 and Theta < 85.67:
        classT = 27
    elif Theta >= 85.67 and Theta < 87.80:
        classT = 28
    elif Theta >= 87.80 and Theta <= 90.00:
        classT = 29
    return classT

# maxDist bin for 3D
def dist12Class_(dist12):
    #classL=0
    if (dist12<3.83):
        classL=1
    elif dist12>=3.83 and dist12<7.00:
        classL=2
    elif dist12>=7.00 and dist12<9.00:
        classL=3
    elif dist12>=9.00 and dist12<11.00:
        classL=4
    elif dist12>=11.00 and dist12<14.00:
        classL=5
    elif dist12>=14.00 and dist12<17.99:
        classL=6
    elif dist12>=17.99 and dist12<21.25:
        classL=7
    elif dist12>=21.25 and dist12<23.19:
        classL=8
    elif dist12>=23.19 and dist12<24.8:
        classL=9
    elif dist12>=24.8 and dist12<26.26:
        classL=10
    elif dist12>=26.26 and dist12<27.72:
        classL=11
    elif dist12>=27.72 and dist12<28.9:
        classL=12
    elif dist12>=28.9 and dist12<30.36:
        classL=13
    elif dist12>=30.36 and dist12<31.62:
        classL=14
    elif dist12>=31.62 and dist12<32.76:
        classL=15
    elif dist12>=32.76 and dist12<33.84:
        classL=16
    elif dist12>=33.84 and dist12<35.13:
        classL=17
    elif dist12>=35.13 and dist12<36.26:
        classL=18
    elif dist12>=36.26 and dist12<37.62:
        classL=19
    elif dist12>=37.62 and dist12<38.73:
        classL=20
    elif dist12>=38.73 and dist12<40.12:
        classL=21
    elif dist12>=40.12 and dist12<41.8:
        classL=22
    elif dist12>=41.8 and dist12<43.41:
        classL=23
    elif dist12>=43.41 and dist12<45.55:
        classL=24
    elif dist12>=45.55 and dist12<47.46:
        classL=25
    elif dist12>=47.46 and dist12<49.69:
        classL=26
    elif dist12>=49.69 and dist12<52.65:
        classL=27
    elif dist12>=52.65 and dist12<55.81:
        classL=28
    elif dist12>=55.81 and dist12<60.2:
        classL=29
    elif dist12>=60.2 and dist12<64.63:
        classL=30
    elif dist12>=64.63 and dist12<70.04:
        classL=31
    elif dist12>=70.04 and dist12<76.15:
        classL=32
    elif dist12>=76.15 and dist12<83.26:
        classL=33
    elif dist12>=83.26 and dist12<132.45:
        classL=34
    elif dist12>=132.45:
        classL=35
    return classL



def calDist(x1, y1, z1, x2, y2, z2):
    return math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2 + (z2 - z1) ** 2)


Key_Dict_Total = ()
DataFrame_Index = []
Group_Information = []

Output_Folder_Path = '/ddnB/work/wxx6941/TSR/code/code/psi_revision/9pdb/ca_tsr_a_b/output_ca_all'
outputFile3=open('Sample_details2_output.txt','w')
for i in range(len(PDB_list)):
    PDB_ID = PDB_list[i]
    Chain_Name = Chain[i]
    Group = Group_Infor[i]

    PDB_File_Path = "/ddnB/work/wxx6941/TSR/code/code/psi_revision/9pdb/ca_tsr_a_b/pdb/{}.pdb".format(PDB_ID)

    keyDict3D = {}
    xCoord = {}
    yCoord = {}
    zCoord = {}
    Atom = {}
    Res={}
    res_Id={}

    counter2 = 0
    File=open(PDB_File_Path,'r')
    lines_ = File.readlines()
    outputFile1 = open( f'{Output_Folder_Path}/{PDB_ID}_{Chain_Name}.triplets_theta29_dist35', 'w')
    outputFile2 = open(f'{Output_Folder_Path}/key_{PDB_ID}_{Chain_Name}.keys_theta29_dist35', 'w')
    # header
    outputFile1.writelines( 'Residue1   Residue2   Residue3   Edge1  Edge2  Edge3\t   Coor_R1\t           Coor_R2\t         CoorR3\tTheta\tmax_dist\td_3\tkey3D\tLabel1\tLabel2\tLabel3\n')
    outputFile2.writelines('key\t\tfreq\n')

    counter = 0
    for line in lines_:
        if line.startswith('ATOM') or line.startswith('HETATM'):
            if line[21].strip() == Chain_Name:
                if line[13:15].strip() == 'CA':
                    xCoord[counter] = float(line[30:38])
                    yCoord[counter] = float(line[38:46])
                    zCoord[counter] = float(line[46:54])

                    Res[counter] = line[17:20].strip()

                    res_Id[counter] = line[22:27].strip()
                    Atom[counter] = line[13:15].strip()
                    counter += 1

                    print(line[22:27].strip())
        if line.startswith('TER') and line[21].strip() == Chain_Name:
            break
    Max_Dist={}
    Min_Dist={}
    Number_of_Triangles=0
    for i in range(len(xCoord)):
        for j in range(i + 1, len(yCoord)):
            for k in range(j + 1, len(zCoord)):
                L1 = calDist(xCoord[i], yCoord[i], zCoord[i], xCoord[j], yCoord[j], zCoord[j])
                L2 = calDist(xCoord[j], yCoord[j], zCoord[j], xCoord[k], yCoord[k], zCoord[k])
                L3 = calDist(xCoord[i], yCoord[i], zCoord[i], xCoord[k], yCoord[k], zCoord[k])

                l1 = int(atomSeq[Res[j]])
                l2 = int(atomSeq[Res[k]])
                l3 = int(atomSeq[Res[i]])

                Med1 = (1 / 2) * math.sqrt(2 * (L1 ** 2) + 2 * (L2 ** 2) - L3 ** 2)
                Med2 = (1 / 2) * math.sqrt(2 * (L2 ** 2) + 2 * (L3 ** 2) - L1 ** 2)
                Med3 = (1 / 2) * math.sqrt(2 * (L3 ** 2) + 2 * (L1 ** 2) - L2 ** 2)
                Median = [Med1, Med2, Med3]
                Label = [l1, l2, l3]
                index1 = [L3, L1, L2]

                # 1st Condition
                if l1 != l2 != l3:
                    X = [l1, l2, l3]
                    b3 = Median[Label.index(min(l1, l2, l3))]
                    d12 = index1[Label.index(min(l1, l2, l3))]
                    if d12 == L3 and max(l1, l2, l3) == l2:
                        d13 = L2
                    elif d12 == L3 and max(l1, l2, l3) == l3:
                        d13 = L1

                    elif d12 == L2 and max(l1, l2, l3) == l1:
                        d13 = L1
                    elif d12 == L2 and max(l1, l2, l3) == l2:
                        d13 = L3
                    elif d12 == L1 and max(l1, l2, l3) == l1:
                        d13 = L2
                    elif d12 == L1 and max(l1, l2, l3) == l3:
                        d13 = L3
                    X.remove(max(X))
                    X.remove(min(X))
                    Label1 = max(l1, l2, l3)
                    Label2 = X[0]
                    Label3 = min(l1, l2, l3)
                    #outputFile1.writelines(f'{Label1}  {Label2}  {Label3} {max(l1, l2, l3)} {l1} {l2} {l3}')

                # 2nd condition
                elif l1 > l2 == l3:
                    Label1 = l1
                    if L2 > L1:
                        b3 = Med3
                        d13 = L1
                        d12 = L2
                        Label2 = l2
                        Label3 = l3
                    else:
                        b3 = Med2
                        d13 = L2
                        d12 = L1
                        Label2 = l3
                        Label3 = l2

                elif l2 > l1 == l3:
                    Label1 = l2
                    if L3 > L2:
                        b3 = Med1
                        d13 = L2
                        d12 = L3
                        Label2 = l3
                        Label3 = l1
                    else:
                        b3 = Med3
                        d13 = L3
                        d12 = L2
                        Label2 = l1
                        Label3 = l3

                elif l3 > l1 == l2:
                    Label1 = l3
                    if L1 > L3:
                        b3 = Med2
                        d13 = L3
                        d12 = L1
                        Label2 = l1
                        Label3 = l2
                    else:
                        b3 = Med1
                        d13 = L1
                        d12 = L3
                        Label2 = l2
                        Label3 = l1
                # 3rd condition
                elif l1 == l2 > l3:
                    b3 = Med3
                    Label3 = l3
                    if L1 > L3:
                        d13 = L1
                        d12 = L2
                        Label1 = l1
                        Label2 = l2
                    else:
                        d13 = L3
                        d12 = L2
                        Label1 = l2
                        Label2 = l1

                elif l1 == l3 > l2:
                    Label3 = l2
                    b3 = Med2
                    if L2 > L3:
                        d13 = L2
                        d12 = L1
                        Label1 = l1
                        Label2 = l3
                    else:
                        d13 = L3
                        d12 = L1
                        Label1 = l3
                        Label2 = l1
                elif l2 == l3 > l1:
                    Label3 = l1
                    b3 = Med1
                    if L2 > L1:
                        d13 = L2
                        d12 = L3
                        Label1 = l2
                        Label2 = l3
                    else:
                        d13 = L1
                        d12 = L3
                        Label1 = l3
                        Label2 = l2

                # 4th condition
                if l1 == l2 == l3:
                    if L2 >= max(L1, L2, L3):
                        b3 = Med3
                        d13 = L1
                        d12 = L2
                        Label1 = l1
                        Label2 = l2
                        Label3 = l3
                    if L1 >= max(L1, L2, L3):
                        b3 = Med2
                        d13 = L2
                        d12 = L1
                        Label1 = l1
                        Label2 = l3
                        Label3 = l2

                    if L3 >= max(L1, L2, L3):
                        # b3=Med3
                        # d13 =L1
                        # d12 =L2
                        # Corrected
                        b3 = Med1
                        d13 = L1
                        d12 = L3
                        Label1 = l3
                        Label2 = l2
                        Label3 = l1

                a = (d13 ** 2 - (d12 / 2) ** 2 - b3 ** 2)
                b = (2 * (d12 / 2) * b3)
                if L1!=0 and L2!=0 and L3!=0:
                    Theta1 = (math.acos(a / b)) * (180 / math.pi)

                    if Theta1 <= 90:
                        Theta = Theta1
                    else:
                        Theta = abs(180 - Theta1)
                    maxDist = max(L1, L2, L3)
                    minDist = min(L1, L2, L3)

                    Max_Dist[Number_of_Triangles] = maxDist
                    Min_Dist[Number_of_Triangles] = minDist
                    Number_of_Triangles += 1

                    ClassT1 = thetaClass_(Theta)
                    ClassL1 = dist12Class_(maxDist)
                    if maxDist >= 0 and maxDist <= 5000:
                        key3D = dLen * dTheta * (numOfLabels ** 2) * (int(Label1) - 1) + \
                                dLen * dTheta * (numOfLabels) * (int(Label2) - 1) + \
                                dLen * dTheta * (int(Label3) - 1) + \
                                dTheta * (ClassL1 - 1) + \
                                (ClassT1 - 1)

                        if key3D in keyDict3D:
                            keyDict3D[key3D] += 1
                        else:
                            keyDict3D[key3D] = 1

                        outputFile1.write("{}_{}_{}_{}  ".format(Res[i], Chain_Name, res_Id[i], Atom[i]))
                        outputFile1.write("{}_{}_{}_{}  ".format(Res[j], Chain_Name, res_Id[j], Atom[j]))
                        outputFile1.write("{}_{}_{}_{}  ".format(Res[k], Chain_Name, res_Id[k], Atom[k]))
                        outputFile1.write(" {:.2f}  {:.2f}  {:.2f} ".format(L1, L2, L3))
                        outputFile1.write(
                            "  {:.2f},{:.2f},{:.2f}   {:.2f},{:.2f},{:.2f}  ".format(xCoord[i], yCoord[i], zCoord[i],
                                                                                     xCoord[j], yCoord[j], zCoord[j]))

                        outputFile1.write(" {:.2f},{:.2f},{:.2f}  ".format(xCoord[k], yCoord[k], zCoord[k]))
                        outputFile1.write(
                            "{:.2f}   {:.2f}   {:.2f}   {:.0f} {} {} {} \n".format(Theta, maxDist, b3, key3D,
                                                                                   int(Label1),
                                                                                   int(Label2), int(Label3)))
    for value_ in keyDict3D:
        outputFile2.writelines([str(value_), '\t', str(keyDict3D[value_]), '\n'])

    Key_Dict_Total += (keyDict3D,)

    DataFrame_Index.append(f'{PDB_ID}_{Chain_Name}')
    Group_Information.append(Group)

    outputFile3.writelines(
                    f'{PDB_ID}\t{Chain_Name}\t{len(Res)}\t{len(keyDict3D)}\t{sum(keyDict3D.values())}\t{max(Max_Dist.values())}\t{min(Max_Dist.values())}\n')

    outputFile1.close()
    outputFile2.close()



outputFile3.close()
df=pd.DataFrame(Key_Dict_Total,index=Group_Information)
df=df.rename_axis('group')
df=df.fillna(0)
df = df.astype('int')
df.insert(0,'protein',DataFrame_Index)
df.to_csv(f"{Output_Folder_Path}/feature_map_with_header.csv",header=True,index=True)

df_Group=pd.DataFrame(columns=['group'],index=DataFrame_Index)
df_Group=df_Group.rename_axis('protein')
df_Group['group']=Group_Information
df_Group.to_csv(f"{Output_Folder_Path}/sample_details.csv",header=True,index=True)

df_Clustering=pd.DataFrame(Key_Dict_Total,index=DataFrame_Index)
df_Clustering=df_Clustering.fillna(0)
df_Clustering = df_Clustering.astype('int')

first_column = list(df_Clustering.iloc[:, 0])
DataFrame_Index_2=[]
for i in range(len(df_Clustering)):
    DataFrame_Index_2.append(DataFrame_Index[i]+";"+str(int(first_column[i])))
df_Clustering.iloc[:, 0]=DataFrame_Index_2
df_Clustering.to_csv(f"{Output_Folder_Path}/localFeatureVect_theta29_dist35_NoFeatureSelection_keyCombine0.csv",header=False,index=False)
print('completed Successfully')



