#!/usr/bin/env python
import pandas as pd
import pdbio
import os
import sys
import numpy as np
from Bio import Align

three_to_one = {
    'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
    'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
    'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
    'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y',
    'SEC': 'U', 'PYL': 'O'  # Add rare amino acids if needed
}

def get_geometric_CA_center(df,segments=None):
    if type(segments) != type(None):
        #print(segments)
        df_CA=df.loc[(df["atname"]=="CA") & (df["segid"].isin(segments))][['x','y','z']]
    else:
        df_CA=df.loc[df["atname"]=="CA"][['x','y','z']]
    return df_CA.sum(axis=0)/df_CA.shape[0]

def center_by_CA(df,segments=None):
    center_current=get_geometric_CA_center(df,segments)
    #print(center_current)
    df[['x','y','z']]-=center_current
    #return df

# Does not check that both proteins are identical before attempting the join.
def relational(mapper):
    """
    Takes a list of tuples as an argument, and then creates a relational dataframe mapping the tuple values to each other.
    Dataframe columns are "R0" and "R1".
    Each tuple should have a length of 2.
    """
    schema = {'R0':str,'R1':str}
    df = pd.DataFrame(columns=schema.keys()).astype(schema)
    for maptuple in mapper:
        df = df._append({'R0':maptuple[0],'R1':maptuple[1]},ignore_index=True)
    return df


def extract_sequence(df, chainid):
    # Filter the dataframe by the specified chainid and sort by residue number
    filtered_df = df[df['segid'] == chainid].sort_values('resnum').groupby('resnum').first().reset_index()['resname']

    # Convert the 'resname' to one-letter codes using the mapping dictionary
    sequence = ''.join(three_to_one.get(res, 'X') for res in filtered_df)
    #print(sequence)
    return sequence

def map_chains(df1,df2):
    chains = df1.segid.unique()
    chains2 = df2.segid.unique()

    if len(chains)==0 or len(chains2)==0:
        raise ValueError('One or more dataframes has no segids labeled')

    if len(chains) == 1 and len(chains2==1):
        mapper = []
        mapper.append((chains[0],chains2[0]))
    else:
        mapper = []
        for chain in chains:
            for chain2 in chains2:
                if chain==chains2:
                    mapper.append((chain,chain2))

    if len(mapper)==0:
        raise VaueError('No matching segids found')
    return mapper


def different_sequence_RMSD(df1,df2,mapper=None):
    """
    Get RMSD when sequences are not identical. Does a gradient alignment.
    """
#   #print('df1:\n',df1.head())
#   #print('df1:\n',df1.segid.unique())
#   #print('df2:\n',df2.head())
    if mapper is None:
        mapper = map_chains(df1,df2)

    map_df = relational(mapper)
    map_df['alignment'] = np.nan
    map_df['RMSD'] = np.nan


    schema = {'segid1':str,'segid2':str,'index1':int, 'index2':int}
    alignment_df = pd.DataFrame(columns=schema.keys()).astype(schema)



    for index,map_row in map_df.iterrows():
        sequence1 = extract_sequence(df1,map_row['R0'])
        sequence2 = extract_sequence(df2,map_row['R1'])
#       #print('sequence1:\n',sequence1)
#       #print('sequence2:\n',sequence2)

        aligner=Align.PairwiseAligner()
        alignments = aligner.align(sequence1,sequence2)
        alignment = alignments[0]
        normalized_score = alignment.score / max(len(sequence1), len(sequence2))
        map_df.at[index,'alignment'] = normalized_score

        aligned_seq1, aligned_seq2 = alignment.aligned
        for region1, region2 in zip(aligned_seq1, aligned_seq2):
            for i, j in zip(range(region1[0], region1[1]), range(region2[0], region2[1])):
                alignment_df=alignment_df._append({'segid1':map_row['R0'],'segid2':map_row['R1'],'index1':i, 'index2':j},ignore_index=True)

    alignment_df.index1+=1 #Residue number starts count at 1
    alignment_df.index2+=1

    #alignment_df.to_csv('alignment_df.csv')
    #print('alignment_df:\n',alignment_df.head())
    #print('df1 B:\n',df1.loc[df1['segid']=='B'][['segid','resnum']])
    # Perform an inner join between df1 and alignment_df based on segid and resnum = index1
    df1_filtered = pd.merge(df1, alignment_df, how='outer', left_on=['segid', 'resnum'], right_on=['segid1', 'index1'],suffixes=['','_r'])
    df2_filtered = pd.merge(df2, alignment_df, how='outer', left_on=['segid', 'resnum'], right_on=['segid2', 'index2'],suffixes=['','_r'])

    # Drop unnecessary columns if needed (like 'segid' and 'index1' after the merge)
    df1_filtered.rename(columns={'index1': 'alignment_index'}, inplace=True)
    df2_filtered.rename(columns={'index1': 'alignment_index'}, inplace=True)

    df1_filtered = df1_filtered.drop(columns=['index2', 'segid1','segid2']).reset_index(drop=True)
    df2_filtered = df2_filtered.drop(columns=['segid1','segid2']).reset_index(drop=True)
    #print('df2_filtered shape: ',df2_filtered.shape)
   #print('df1_filtered:\n',df1_filtered.head())

    center_by_CA(df1_filtered,map_df['R0'].unique())
   #print('df1_filtered (centered):\n',df1_filtered.head())
    center_by_CA(df2_filtered,map_df['R1'].unique())
    print('centered both dataframes, about to do gradient align')

    return kabsch_align(df1_filtered,df2_filtered,map_df),df1_filtered,df2_filtered




def merge_CA_dataframes(df0,df1,map_df):
    df0_CA=df0.loc[df0["atname"]=="CA"]
    df1_CA=df1.loc[df1["atname"]=="CA"]

    df_mid=pd.merge(df0_CA, map_df, left_on="segid",right_on="R0",how="inner")
    if 'alignment_index' in df1_CA.keys():
        df2=pd.merge(df_mid,df1_CA,left_on=["alignment_index","R1"],right_on=["alignment_index","segid"],how="inner",suffixes=["_l","_r"])
    else:
        df2=pd.merge(df_mid,df1_CA,left_on=["resnum","R1"],right_on=["resnum","segid"],how="inner",suffixes=["_l","_r"])
    return df2

# Assumes both dataframes are centered
def get_CA_RMSD(df0,df1,map_df=None,whitelist_region=None):
    if map_df is None:
        map_df = relational(map_chains(template_struct,modify_struct))

    df2 = merge_CA_dataframes(df0,df1,map_df)

    df2['dx']=df2.x_l-df2.x_r
    df2['dy']=df2.y_l-df2.y_r
    df2['dz']=df2.z_l-df2.z_r
    #map_df.to_csv('map_df.csv')
    #df2.to_csv('RMSD_df2.csv')
    diffs = df2[['dx','dy','dz']]

    # Compute individual RMSD for each group (R0, R1)
    grouped_rmsd = df2.groupby(['R0', 'R1']).apply(lambda group: np.sqrt(np.mean(np.sum(np.square(group[['dx', 'dy', 'dz']]), axis=1))))

    # Update the map_df with RMSD values
    map_df = map_df.set_index(['R0', 'R1'])  # Ensure map_df has R0 and R1 as index
    try:
        map_df['RMSD'] = grouped_rmsd  # Update the map_df with individual RMSD values
    except:
       print('couldn\'t set map_df\'s RMSD value.')
       print('Grouped_rmsd:\n',grouped_rmsd)
       print('Map_df: \n',map_df.head())
       print('df0: \n',df0[['alignment_index',"segid"]].head())
       print('df0: \n',df0.head())
       print('df1: \n',df1.head())
       print('df2: \n',df2[['dx','dy','dz','R0','R1']].head())



    rmsd = np.sqrt(np.sum(np.sum(diffs**2,axis=0),axis=0)/diffs.shape[0])
    if np.isnan(rmsd):
       print(diffs)
       print(diffs.shape)

    return rmsd,map_df

# Define a rotation matrix
def rotation_matrix(axis, theta):
    """
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians.
    """
    axis = np.asarray(axis)
    axis = axis / np.sqrt(np.dot(axis, axis))
    a = np.cos(theta / 2.0)
    b, c, d = -axis * np.sin(theta / 2.0)
    aa, bb, cc, dd = a * a, b * b, c * c, d * d
    bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
    return np.array([[aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
                     [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
                     [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc]])

def kabsch_align(template_struct, modify_struct, map_df=None):
    if map_df is None:
        map_df = relational(map_chains(template_struct,modify_struct))

    df_merge = merge_CA_dataframes(template_struct,modify_struct,map_df)
    template_coords = df_merge[['x_l', 'y_l', 'z_l']].to_numpy()
    modify_coords = df_merge[['x_r', 'y_r', 'z_r']].to_numpy()

    # Center the matched points by their centroids
    centroid_template = np.mean(template_coords, axis=0)
    centroid_modify = np.mean(modify_coords, axis=0)
    centered_template_coords = template_coords - centroid_template
    centered_modify_coords = modify_coords - centroid_modify

    # Compute the covariance matrix and SVD
    covariance_matrix = np.dot(centered_modify_coords.T, centered_template_coords)
    V, S, Wt = np.linalg.svd(covariance_matrix)
    d = np.sign(np.linalg.det(np.dot(V, Wt)))
    rotation_matrix = np.dot(V, np.dot(np.diag([1, 1, d]), Wt))

    # Apply the transformation (rotation + translation) to all atoms in modify_struct
    all_modify_coords = modify_struct[['x', 'y', 'z']].to_numpy() - centroid_modify
    aligned_coords = np.dot(all_modify_coords, rotation_matrix) + centroid_template
    modify_struct[['x', 'y', 'z']] = aligned_coords

    # Calculate the final RMSD using only the matched Cα atoms
    final_rmsd = get_CA_RMSD(template_struct, modify_struct, map_df)
    return final_rmsd

def align_all_structures(project,whitelist_region=None):
    if project[-1]!='/':
        project+='/'
    files = os.listdir(project)
    pdblist = []
    for ufile in files:
        if ufile[-4:]=='.pdb':
            pdblist.append(ufile)
    df0 = pdbio.read_pdb(project+pdblist[0])
    center_by_CA(df0)
    pdbio.write_pdb(df0,project+pdblist[0][:-4]+'_centered.pdb')
    for i in pdblist[1:]:
        df1 = pdbio.read_pdb(project+i)
        center_by_CA(df1)
        kabsch_align(df0,df1)
        print(pdblist[0]+' + '+i+' RMSD:\t',get_CA_RMSD(df0,df1))
        pdbio.write_pdb(df1,project+i[:-4]+'_aligned.pdb')

if __name__ == '__main__':
    align_all_structures(sys.argv[1])
