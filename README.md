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
Once the protein coordinate file is loaded into a dataframe, this file contains a set of functions to operate on it. What follows are a list of some of the most useful capabilities:

- kabsch alignment on two structures
- different-sequence RMSD that aligns and calculates RMSD on two non-identical proteins based on segments with hight sequence similarity
- Parallel computation of an RMSD matrix, though it assumes you are dealing with identical versions of the same protein, can be useful for interpreting output from molecular dynamics simulations 
- Get the chirality of residues in a protein structure for use in setting inputs for some MD simulations
