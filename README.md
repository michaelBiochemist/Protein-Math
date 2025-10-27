# Protein-Math
A set of tools in python for analyzing proteins as dataframes.

For those trying to install it on older versions of python, from the shared repo, do the following steps:
	1. Find the correct .whl package
	2. pip install --upgrade --ignore-requires-python dist/pmike-0.0.\<number\>-\*whl

# Explanation of functions and algorithms #

# pdbio #
PDBIO is responsible for file format reading and conversion. The most reliable formats to use are PDB files, though it can also read from .crd files as used in CHARMM MD software. 

The primary functions to do this are "read_coord" and "write_pdb" or "write_crd". The former requires only the file name as argument ard returns a pandas dataframe. The latter requires the dataframe, and optionally accepts filename and comments, otherwise defaulting to output.<file-extension> and "", respectively. 

If a dataframe is longer than 100,000 rows, then the resulting PDB will be broken up into separate files in accordance with the PDB spec. In order to read it back in, one can call the "read_split_pdb" function. Beware, this will read in all pdb files in a given directory. 

# protein_math #
Once the protein coordinate file is loaded into a dataframe, this file contains a set of functions to operate on it. 

get_geometric_CA_center(df,segments=None,selection_criteria=None):
get_volume(df):
center_by_CA(df,segments=None,selection_criteria=None):
relational(mapper):
extract_sequence(df, chainid):
map_chains(df1,df2):
get_alignments(df1,df2,map_df,mode="global",match_score=2,mismatch_score=-1,
get_alignment_df(df1,df2,map_df,mode="global",match_score=2,mismatch_score=-1,
different_sequence_RMSD(df1,df2,mapper=None,alignment_params=None,center_selection_criteria=None):
same_sequence_RMSD(df1,df2,mapper=None):
merge_CA_dataframes(df0,df1,map_df):
get_CA_RMSD(df0,df1,map_df=None,whitelist_region=None,debuggy=False):
# Define a rotation matrix
rotation_matrix(axis, theta):
kabsch_align(template_struct, modify_struct, map_df=None):
pairwise_rmsd(pairwise_tuple):
compute_some_matrix_parallel(args_list, labels, function, max_workers=None,print_count=1000):
compute_rmsd_matrix_parallel(structures, labels, num_jobs=-1):
calculate_chirality_proline_old(group):
calculate_chirality_proline(residue_df):
calculate_residue_chirality(residue_df):
toggle_residue_chirality(residue_df):
get_chirality(df):
assign_histidine_charge(df, cutoff=4.0):
