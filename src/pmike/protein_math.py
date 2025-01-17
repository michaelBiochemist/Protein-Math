#!/usr/bin/env python
import pandas as pd
from pmike import pdbio
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

def get_geometric_CA_center(df,segments=None,selection_criteria=None):
    df_CA = None
    df_CA=df.loc[df["atname"]=="CA"]
    if type(segments) != type(None):
        df_CA=df_CA.loc[(df_CA["segid"].isin(segments))]
    if type(selection_criteria) == str:
        if selection_criteria == 'core-sequence':
            core_df = pd.DataFrame()

            for chain_id, chain_df in df_CA.groupby('segid'):
                min_resnum = chain_df['resnum'].min()
                max_resnum = chain_df['resnum'].max()
                offset = int(max_resnum * 0.15) # Shave off 15% of residues from either side
                start_resnum = min_resnum + offset
                end_resnum = max_resnum - offset

                core_chain_df = chain_df.loc[(chain_df['resnum'] <= end_resnum) & (chain_df['resnum'] >= start_resnum)]
                core_df = pd.concat([core_df,core_chain_df])
            df_CA = core_df

    return df_CA[['x','y','z']].sum(axis=0)/df_CA[['x','y','z']].shape[0]

def get_volume(df):
    maxes = df[['x','y','z']].max()
    mins  = df[['x','y','z']].min()
    diffs = maxes-mins
    return diffs['x']*diffs['y']*diffs['z']

def center_by_CA(df,segments=None,selection_criteria=None):
    center_current=get_geometric_CA_center(df,segments,selection_criteria)
    df[['x','y','z']]-=center_current
    # Does not need a return value as this centers the dataframe passed in the argument (pass by reference)

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
    filtered_df = df[(df['segid'] == chainid) & (df['atname']=="CA")].sort_values('resnum').groupby('resnum').first().reset_index()[['resnum','resname']]
    index_offset = filtered_df.resnum.min()

    # Convert the 'resname' to one-letter codes using the mapping dictionary
    sequence = ''.join(three_to_one.get(res, 'X') for res in filtered_df['resname'])
    if len(sequence)==0:
        return '',0

    # Extract the residue numbers
    residue_numbers = filtered_df['resnum'].values

    # Check if the sequence of residue numbers is contiguous
    is_contiguous = (len(residue_numbers) == (residue_numbers[-1] - residue_numbers[0] + 1))

    # Optionally, find the missing residue numbers
    if not is_contiguous:
        expected_resnums = set(range(residue_numbers[0], residue_numbers[-1] + 1))
        missing_resnums = expected_resnums - set(residue_numbers)
        print("The residue numbers are not contiguous.")
        print(f"Chain {chainid} is missing residue numbers: {sorted(missing_resnums)}")

    return sequence,index_offset

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
                if chain==chain2:
                    mapper.append((chain,chain2))

    if len(mapper)==0:
        raise VaueError('No matching segids found')
    return mapper

def get_alignments(df1,df2,map_df,mode="global",match_score=2,mismatch_score=-1,
        gap_open=-10,gap_extend=-0.5):
    alignment_list = []
    for index,map_row in map_df.iterrows():
        sequence1,index_offset1  = extract_sequence(df1,map_row['R0'])
        sequence2,index_offset2 = extract_sequence(df2,map_row['R1'])

        if len(sequence1)==0 or len(sequence2)==0:
            map_df.at[index,'alignment'] = np.nan
            continue

        #print('sequence1 :',sequence1)
        #print('sequence2 :',sequence2)

        aligner=Align.PairwiseAligner()
        aligner.mode=mode
        aligner.match_score=match_score
        aligner.mismatch_score=mismatch_score
        aligner.open_gap_score=gap_open
        aligner.extend_gap_score=gap_extend

        alignments = aligner.align(sequence1,sequence2)
        alignment = alignments[0]
        normalized_score = alignment.score / max(len(sequence1), len(sequence2))
        map_df.at[index,'alignment'] = normalized_score

        aligned_seq1 = alignment.aligned[0]+index_offset1
        aligned_seq2 = alignment.aligned[1]+index_offset2

        alignment_list.append({'seq1':aligned_seq1,'seq2':aligned_seq2,'score':normalized_score})
    return alignment_list

def get_alignment_df(df1,df2,map_df,mode="global",match_score=2,mismatch_score=-1,
        gap_open=-10,gap_extend=-0.5):
    rows = []
    for index,map_row in map_df.iterrows():
        sequence1,index_offset1 = extract_sequence(df1,map_row['R0'])
        sequence2,index_offset2 = extract_sequence(df2,map_row['R1'])

        if len(sequence1)==0 or len(sequence2)==0:
            map_df.at[index,'alignment'] = np.nan
            continue

        #print('sequence1 :',sequence1)
        #print('sequence2 :',sequence2)

        aligner=Align.PairwiseAligner()
        aligner.mode=mode
        aligner.match_score=match_score
        aligner.mismatch_score=mismatch_score
        aligner.open_gap_score=gap_open
        aligner.extend_gap_score=gap_extend

        alignments = aligner.align(sequence1,sequence2)
        alignment = alignments[0]
        normalized_score = alignment.score / max(len(sequence1), len(sequence2))
        map_df.at[index,'alignment'] = normalized_score

        aligned_seq1 = alignment.aligned[0]+index_offset1
        aligned_seq2 = alignment.aligned[1]+index_offset2
        """
        for region1, region2 in zip(aligned_seq1, aligned_seq2):
            for i, j in zip(range(region1[0], region1[1]), range(region2[0], region2[1])):
                alignment_df=alignment_df._append({'segid1':map_row['R0'],'segid2':map_row['R1'],'index1':i, 'index2':j},ignore_index=True)
        """
        for region1, region2 in zip(aligned_seq1, aligned_seq2):
            for i, j in zip(range(region1[0], region1[1]), range(region2[0], region2[1])):
                rows.append({'segid1': map_row['R0'], 'segid2': map_row['R1'], 'index1': i, 'index2': j})
    alignment_df = pd.DataFrame(rows)
    alignment_df.index1+=1 #Residue number starts count at 1
    alignment_df.index2+=1

    alignment_df['index1'] = alignment_df['index1'].astype('Int64')
    alignment_df['index2'] = alignment_df['index2'].astype('Int64')
    return alignment_df


def different_sequence_RMSD(df1,df2,mapper=None,alignment_params=None,center_selection_criteria=None):
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
    #alignment_df = pd.DataFrame(columns=schema.keys()).astype(schema)
    #print('alignment params: ',alignment_params)
    if alignment_params==None:
        alignment_df=get_alignment_df(df1,df2,map_df)
    else:
        alignment_df=get_alignment_df(df1,df2,map_df,mode=alignment_params['mode'],
                match_score=alignment_params['match_score'],mismatch_score=alignment_params['mismatch_score'],
                gap_open=alignment_params['gap_open'],gap_extend=alignment_params['gap_extend'])

    #alignment_df.to_parquet('alignment_df.parquet')

    df1_filtered = pd.merge(df1, alignment_df, how='left', left_on=['segid', 'resnum'], right_on=['segid1', 'index1'],suffixes=['','_r'])
    df2_filtered = pd.merge(df2, alignment_df, how='left', left_on=['segid', 'resnum'], right_on=['segid2', 'index2'],suffixes=['','_r'])

    # Drop unnecessary columns if needed (like 'segid' and 'index1' after the merge)
    df1_filtered.rename(columns={'index1': 'alignment_index'}, inplace=True)
    df2_filtered.rename(columns={'index1': 'alignment_index'}, inplace=True)

    df1_filtered = df1_filtered.drop(columns=['index2', 'segid1','segid2']).reset_index(drop=True)
    df2_filtered = df2_filtered.drop(columns=['segid1','segid2']).reset_index(drop=True)

    center_by_CA(df1_filtered,segments=map_df['R0'].unique(),selection_criteria=center_selection_criteria)
    center_by_CA(df2_filtered,segments=map_df['R1'].unique(),selection_criteria=center_selection_criteria)

    return kabsch_align(df1_filtered,df2_filtered,map_df,center_selection_criteria=center_selection_criteria),df1_filtered,df2_filtered




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
def get_CA_RMSD(df0,df1,map_df=None,whitelist_region=None,debuggy=False):
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
    normalized_rmsd = grouped_rmsd / np.sqrt(df2.groupby(['R0', 'R1']).size())
    if debuggy:
        df2['dd'] =np.sqrt(df2.dx**2 + df2.dy**2 + df2.dz**2)
        df2.to_parquet('df2.parquet')


    # Update the map_df with RMSD values
    map_df = map_df.set_index(['R0', 'R1'])  # Ensure map_df has R0 and R1 as index
    try:
        map_df['RMSD'] = grouped_rmsd  # Update the map_df with individual RMSD values
        map_df['RMSD (N)'] = normalized_rmsd  # Update the map_df with individual RMSD values
    except:
       print('couldn\'t set map_df\'s RMSD value.')
       print('Grouped_rmsd:\n',grouped_rmsd)
       print('Map_df: \n',map_df.head())
       print('df0: \n',df0[['alignment_index',"segid"]].head())
       print('df0: \n',df0.head())
       print('df1: \n',df1.head())
       print('df2: \n',df2[['dx','dy','dz','R0','R1']].head())



    rmsd = np.sqrt(np.sum(np.sum(diffs**2,axis=0),axis=0)/diffs.shape[0])
    rmsd_n = rmsd / np.sqrt(diffs.shape[0])
    if np.isnan(rmsd):
       print(diffs)
       print(diffs.shape)

    return rmsd,map_df,rmsd_n

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

def kabsch_align(template_struct, modify_struct, map_df=None,center_selection_criteria=None):
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
    final_rmsd = get_CA_RMSD(template_struct, modify_struct, map_df,debuggy=False)
    return final_rmsd

def calculate_chirality_proline_old(group):
    # Extract coordinates for key atoms
    try:
        ca = group[group['atname'] == 'CA'][['x', 'y', 'z']].iloc[0].values
        n = group[group['atname'] == 'N'][['x', 'y', 'z']].iloc[0].values
        c = group[group['atname'] == 'C'][['x', 'y', 'z']].iloc[0].values
        # Extract side-chain atoms for proline ring
        ring_atoms = group[group['atname'].isin(['CD', 'CG', 'CB'])][['x', 'y', 'z']].values
        ring_center = np.mean(ring_atoms, axis=0)  # Geometric center of the ring
    except IndexError:
        # Return None if key atoms are missing
        return None

    # Create vectors
    v_nc = n - c
    v_ca = ca - c
    v_ring = ring_center - ca

    # Calculate determinant
    chirality = np.linalg.det([v_nc, v_ca, v_ring])

    # Assign chirality
    if chirality > 0:
        return "L"
    elif chirality < 0:
        return "D"
    else:
        return None

def calculate_chirality_proline(residue_df):
    """
    Calculate the chirality of proline residues using the CA-HA vector along with the N, CA, and C atoms.
    Returns 'L' or 'D' if chirality can be determined, otherwise returns None if HA is missing.
    """
    # Select atoms for chirality check

    # Get atom coordinates (N, CA, C, HA)
    N = residue_df[residue_df['atname'] == 'N'][['x', 'y', 'z']].values
    CA = residue_df[residue_df['atname'] == 'CA'][['x', 'y', 'z']].values
    C = residue_df[residue_df['atname'] == 'C'][['x', 'y', 'z']].values
    HA = residue_df[residue_df['atname'] == 'HA'][['x', 'y', 'z']].values

    if len(HA) == 0:
        # If HA is not present, return None
        print(f'cannot calculate proline {residue_df["resnum"].min()} chirality as alpha H is missing')
        return None

    # Calculate vectors
    N_CA = CA - N  # N-CA vector
    CA_C = C - CA  # CA-C vector
    CA_HA = HA - CA  # CA-HA vector

    N_CA, CA_C, CA_HA = map(lambda v: v.flatten() if len(v.shape) > 1 else v, [N_CA, CA_C, CA_HA])

    # Compute the cross product of N-CA and CA-C to get the normal vector
    normal = np.cross(N_CA, CA_C)
    # Calculate the scalar triple product of the normal vector with the CA-H vector
    scalar_triple_product = np.dot(normal, CA_HA)

    # Determine chirality (right-handed or left-handed)
    if scalar_triple_product > 0:
        return 'L'  # Left-handed (L-proline)
    else:
        return 'D'  # Right-handed (D-proline)

def calculate_residue_chirality(residue_df):
    if residue_df['resname'].iloc[0] == 'PRO':
        return calculate_chirality_proline(residue_df)
        #return None
    # Extract coordinates for key atoms
    try:
        ca = residue_df[residue_df['atname'] == 'CA'][['x', 'y', 'z']].values[0]
        n = residue_df[residue_df['atname'] == 'N'][['x', 'y', 'z']].values[0]
        c = residue_df[residue_df['atname'] == 'C'][['x', 'y', 'z']].values[0]
        cb = residue_df[residue_df['atname'] == 'CB'][['x', 'y', 'z']].values[0]
    except IndexError:
        # Skip if any key atom is missing
        return None

    # Create vectors from CA
    v_n = n - ca
    v_c = c - ca
    v_cb = cb - ca

    # Calculate determinant
    chirality = np.linalg.det([v_n, v_c, v_cb])

    # Assign chirality
    if chirality > 0:
        return "L"
    elif chirality < 0:
        return "D"
    else:
        print(f'error calculating residue {residue_df["resname"].min()} {residue_df["resnum"].min()} chirality')
        return None
def toggle_residue_chirality(residue_df):
    """
    Convert the chirality of a residue between L and D form by reflecting its side-chain atoms using vectorized operations.

    Parameters:
    residue_df (pd.DataFrame): Dataframe containing the residue's atom data.
                               Columns: 'atname', 'x', 'y', 'z', etc.

    Returns:
    pd.DataFrame: A new dataframe with the modified chirality.
    """
    try:
        # Extract backbone atom coordinates
        ca = residue_df[residue_df['atname'] == 'CA'][['x', 'y', 'z']].iloc[0].values
        n = residue_df[residue_df['atname'] == 'N'][['x', 'y', 'z']].iloc[0].values
        c = residue_df[residue_df['atname'] == 'C'][['x', 'y', 'z']].iloc[0].values
    except IndexError:
        # Return None if any key backbone atom is missing
        print(f'error toggling residue {residue_df["resname"]} {residue_df["resnum"]}')
        return None

    # Define the plane for reflection
    v1 = n - ca  # Vector from CA to N
    v2 = c - ca  # Vector from CA to C

    # Normal vector to the plane (cross product of v1 and v2)
    plane_normal = np.cross(v1, v2)
    plane_normal /= np.linalg.norm(plane_normal)  # Normalize the normal vector

    # Identify side-chain atoms (not part of the backbone)
    side_chain_mask = ~residue_df['atname'].isin(['N', 'CA', 'C', 'O'])
    side_chain_atoms = residue_df[side_chain_mask][['x', 'y', 'z']].values

    # Vectorized reflection for all side-chain atoms
    vec_to_atoms = side_chain_atoms - ca  # Vector from CA to each side-chain atom
    projections = np.dot(vec_to_atoms, plane_normal)[:, np.newaxis] * plane_normal  # Projections onto the plane normal
    reflected_coords = side_chain_atoms - 2 * projections  # Reflect through the plane

    # Update the coordinates in the dataframe
    reflected_residue = residue_df.copy()
    reflected_residue.loc[side_chain_mask, ['x', 'y', 'z']] = reflected_coords

    return reflected_residue

def get_chirality(df):
    chirality_series = (
        df.groupby(['resnum', 'chainid'])
        .apply(calculate_residue_chirality)
        .reset_index(name='chirality')  # Convert output to a dataframe
    )

    # Merge chirality results back into the original dataframe
    df = df.merge(chirality_series, on=['resnum', 'chainid'], how='left')
    return df

def assign_histidine_charge(df, cutoff=4.0):
    """
    Assigns protonation states to histidines based on proximity to charged residues.

    Parameters:
        df (pd.DataFrame): DataFrame containing PDB structure information.
        cutoff (float): Distance cutoff in Ångströms to consider an interaction.

    Returns:
        pd.DataFrame: DataFrame with updated 'resname' for histidines.

    (Assumes hydrogens are not included in the PDB file, as it does not add or remove hydrogen atoms directly)
    """

    # Identify histidines
    histidines = df[df['resname'] == 'HIS']
    if histidines.shape[0]==0:
        return df

    columns = ['atnum','resnum','chainid','resname','atname','x','y','z']
    charged = df.loc[(df['resname'].isin(['ASP','GLU','ARG','LYS'])) & (df['atname'].isin(['OD1','OD2','OE1','OE2','NH1','NH2','NE','NZ'])),columns] # Isolate charged residues
    if charged.shape[0] == 0:
        df.loc[df['resname']=='HIS','resname'] = 'HSD'
        return df

    # Nearby acidic resiudes contribute +1 to protonation, whereas basic ones contribute -1
    charged['interaction'] = np.where(charged['resname'].isin(['ARG', 'LYS']),-1,1)
    merged = pd.merge(histidines,charged,how='cross',suffixes=['','_r'])
    # Get distances and filter where distance is less than the cutoff
    merged['dist']=np.linalg.norm(merged[['x','y','z']].values - merged[['x_r','y_r','z_r']].values, axis=1)
    merged_filter=merged.loc[merged['dist']<=cutoff,['atnum','resnum','chainid','atname','dist','interaction']]

    # Step 1: Aggregate interaction values for ND1 and NE2
    interactions = (
        merged_filter.groupby(['resnum', 'chainid', 'atname'])['interaction']
        .sum()
        .unstack(fill_value=0)  # Pivot 'atname' values (ND1, NE2) into columns
    )

    # Step 2: Compare summed interactions for ND1 and NE2
    interactions['protonation'] = np.where(
        interactions.get('ND1', 0) > interactions.get('NE2', 0),  # Compare ND1 vs NE2
        'HSD',  # ND1 protonated
        'HSE'   # NE2 protonated
    )

    # Step 3: Create mapping of (resnum, chainid) to protonation state
    protonation_map = interactions['protonation'].to_dict()
    # Step 4: Update histidines DataFrame with the new protonation states
    histidines.loc[:,'resname'] = histidines.apply(
        lambda row: protonation_map.get((row['resnum'], row['chainid']), 'HSD'),
        axis=1
    )

    unique_histidines=histidines[['resnum','chainid','resname']].drop_duplicates()

    # Final step: Update original dataframe with new Histidine residue names
    df_merge = pd.merge(df,unique_histidines,how='left',on=['resnum','chainid'],suffixes=['','_r'])
    df_merge.loc[df_merge['resname']=='HIS','resname']=df_merge.loc[df_merge['resname']=='HIS','resname_r']
    df_merge.drop(columns=['resname_r'], inplace=True)

    return df_merge


#if __name__ == '__main__':
