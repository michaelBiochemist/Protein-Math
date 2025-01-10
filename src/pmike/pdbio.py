#!/usr/bin/env python

import pandas as pd
from numpy import nan
from datetime import datetime
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
import sys
import os

from Bio.PDB import Structure, Model, Chain, Residue, Atom, MMCIFIO
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import three_to_index
import pandas as pd
"""
PDB TOOLS
Created by Michael van Dyk on 07/19/2024
Last modified: 07/19/2024

 Reference PDB column spec described here:
 https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/tutorials/pdbintro.html
 For the column position tuple, since the postion starts count at zero, but is not inclusive on the right, the start position is 1 less than the position in the web page, but the end position is the same
 List is [data_type, (start, end), left_justify_bit]

Updating to use spec from:
    https://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM


TODO:
    Missing support for rows that start with anything other than "ATOM" or "HETATM", such as "HELIX", "SHEET", "TER", etc.
"""

pdb_spec = { # Column specification for PDB files
    "rowtype": ["str", (0,6),1],		# "Row Type (e.g. "ATOM", "REMARK", etc)"
	 "atnum": ["num", (6,11),0],		# "Atom serial number"
	 "atname": ["str", (12,16),1],	# "Atom name"
	 "altloc": ["str", (16,17),0],	# "Alternate location indicator"
	 "resname": ["str", (17,21),1],		# "Residue name"
	 "chainid": ["str",	(21,22),0],	# "Chain identifier"
	 "resnum": ["num", (22,26),0], 		# "Residue sequence number"
	 "inscode": ["str",(26,27),0],		# "Code for insertions of residues"
	 "x": ["num",(30,38),0],		# "X coordinate"
	 "y": ["num",(38,46),0],		# "Y coordinate"
	 "z": ["num",(46,54),0],		# "Z coordinate"
	 "occupancy": ["str",(54,60),0],		# "Occupancy"
	 "tfactor": ["str",(60,66),0],		# "Temperature factor"
	 "segid": ["str",(72,76),0],		# "Segment identifier"
	 "element": ["str",(76,78),0],		# "Element symbol"
	 "charge": ["str",(78,80),0]        # Charge
     }
crd_spec = { # Columns specification for CHARMM extended CRD files
        "atnum":['num',(0,10),0], # Atom no. sequential
        "ires":['num',(11,20),0], # Residue Index - does not revert to 1 on multimers
        "resname":['str',(22,30),1],
        "atname":['str',(32,40),1],
	 "x": ["num",(42,60),0],		# "X coordinate"
	 "y": ["num",(62,80),0],		# "Y coordinate"
	 "z": ["num",(82,100),0],		# "Z coordinate"
     "segid": ["str",(102,106),1],
     "resnum": ["num",(112,126),1],
     "tfactor": ["num",(125,140),0]
     }
def clean_coordinate_df(df,filetype=None):
    if filetype=='PDB':
        mask = df['inscode'].astype(str).str.isdigit()
        df.loc[mask,'resnum']=df.loc[mask,'resnum']*10+df.loc[mask,'inscode'].astype(int)
        df.loc[mask,'inscode'] = nan

    segids = df.segid.unique()
    chainids = df.chainid.unique()

    if len(segids)==1 and pd.isnull(segids[0]) and (len(chainids)!=1 or not pd.isnull(chainids[0])):
        df.segid=df.chainid

    df['resnum']=df.groupby('segid')['resnum'].transform(lambda x: pd.factorize(x)[0] + 1)

    return df

# Reads PDB into dataframe
def read_coordinates(file_name):
    return read_coord(file_name)
def read_coord(file_name):
    read_spec = {'.PDB':read_pdb,
            '.CIF':read_cif,
            '.CRD':read_crd}
    return read_spec[file_name[-4:].upper()](file_name)


def read_pdb(pdbfil):
    df = read_coordinate_file(pdbfil, pdb_spec)
    df = clean_coordinate_df(df,'PDB')
    return df

def read_pdb_directory(directory):
    pdbs = {}
    for a in os.listdir(directory):
        if a[-4:].lower() == '.pdb':
            pdbs[a[:-4]] = read_pdb('/'.join([directory,a]))
    return pdbs

def read_crd(crdfil):
    return clean_coordinate_df(read_coordinate_file(crdfil, crd_spec),'CRD')

def read_cif(cif_file):
    # Parse the CIF file into a dictionary
    mmcif_dict = MMCIF2Dict(cif_file)

    # Extract relevant columns from the dictionary
    column_dict = {
            "_atom_site.Cartn_x":"x",
            "_atom_site.Cartn_y":"y",
            "_atom_site.Cartn_z":"z",
            "_atom_site.group_PDB":"rowtype",
            "_atom_site.id":"atnum",
            "_atom_site.type_symbol":"element",
            "_atom_site.label_atom_id":"atname",
            "_atom_site.label_comp_id":"resname",
        "_atom_site.label_seq_id":"resnum",
        "_atom_site.label_asym_id":"chainid",
        "_atom_site.occupancy":"occupancy",
        "_atom_site.pdbx_formal_charge":"charge"
    }
    atom_site_columns = list(column_dict.keys())
    # Check if column is present and if not throw error or warning depending on the column.
    # Also throw error if present columns are not all of the same length.
    col0_rownum = len(mmcif_dict[atom_site_columns[0]])
    for col in atom_site_columns:
        if col not in mmcif_dict.keys():
            if col in ["_atom_site.occupancy", "_atom_site.pdbx_formal_charge"]:
                print(f'Warning: column {col} not in mmcif dictionary')
                atom_site_columns.remove(col)
            else:
                print(f'Error: important column {col} not in mmcif dictionary, this corresponds to PDB {column_dict[col]}')
                return None # Throw error here when ready
        elif col0_rownum != len(mmcif_dict[col]):
            print(f'Error: row counts do not match between columns {atom_site_columns[0]}: {col0_rownum} and {col}: {len(mmcif_dict[col])}')
            return None

    # Create a pandas DataFrame from these columns
    data = {col: mmcif_dict.get(col, []) for col in atom_site_columns}
    df = pd.DataFrame(data)

    # Convert numeric columns to appropriate types
    make_numeric_list = ["_atom_site.Cartn_x","_atom_site.Cartn_y","_atom_site.Cartn_z",
        "_atom_site.label_seq_id","_atom_site.occupancy","_atom_site.pdbx_formal_charge","_atom_site.id"]
    for make_numeric in make_numeric_list:
        if make_numeric in atom_site_columns:
            df[make_numeric] = pd.to_numeric(df[make_numeric],errors='coerce')
        else:
            df[make_numeric] = nan

    # Rename columns for easier readability
    df.columns = [
        "x", "y", "z", "rowtype", "atnum", "element",
        "atname", "resname", "resnum", "chainid","occupancy","charge"
    ]
    df["segid"]=df["chainid"]
    df["tfactor"]=nan
    df["altloc"]=nan
    df["inscode"]=''

    df[["atnum","resnum"]] = df[["atnum","resnum"]].fillna(0).astype(int)
    df = clean_coordinate_df(df,'CIF')

    return df

def read_coordinate_file(pdbfil, column_spec):
    colspecs = []
    for cname in column_spec.keys():
        colspecs.append(column_spec[cname][1])

    df = pd.read_fwf(pdbfil, colspecs=colspecs,header=None,names=list(column_spec.keys()), dtype=str)
    if "rowtype" in column_spec.keys():
        df=df.loc[df["rowtype"].isin(["ATOM","HETATM"])]
    elif "ires" in column_spec.keys():
        indy=df.loc[df['ires']=='EXT'].index[0]
        df=df.iloc[indy+1:]
    else:
        return df

    # Make numeric columns numeric
    numeric_columns = []
    for label in column_spec.keys():
        if column_spec[label][0] == "num":
            numeric_columns.append(label)
    df[numeric_columns]=df[numeric_columns].apply(pd.to_numeric)
    return df

# Used by write_pdb function
def format_column(column, start, end,prev_end,ljust_bit,file_type):
    width = end - start
    fullwidth=end-prev_end
    if column.name == 'occupancy':
       return column.apply(lambda x: f"{float(x):0.2f}".rjust(fullwidth) if pd.notna(x) else '1.00'.rjust(fullwidth))
    if column.name == 'tfactor':
        return column.apply(lambda x: f"{float(x):0.2f}".rjust(fullwidth) if pd.notna(x) else '0.00'.rjust(fullwidth))
    if column.name in ['x', 'y', 'z']:
        if file_type == 'PDB':
            return column.apply(lambda x: f"{float(x):8.3f}".rjust(fullwidth) if pd.notna(x) else ''.rjust(fullwidth))
        elif file_type == 'CRD':
            return column.apply(lambda x: f"{float(x):8.10f}".rjust(fullwidth) if pd.notna(x) else ''.rjust(fullwidth))
        else:
            pass
    elif column.name == "tfactor" and file_type=='CRD':
            return column.apply(lambda x: f"{float(x):8.10f}".rjust(fullwidth) if pd.notna(x) else ''.rjust(fullwidth))
    else:
        if ljust_bit==1:
            return column.apply(lambda x: str(x).ljust(width).rjust(fullwidth) if pd.notna(x) else ''.ljust(width).rjust(fullwidth))
        else:
            return column.apply(lambda x: str(x).rjust(fullwidth) if pd.notna(x) else ''.rjust(fullwidth))

# Write pdb from a dataframe. It must contain the columns and data types matching column_spec

def write_crd(df,crdout='output.crd',comments=''):
    if 'ires' not in df.keys() or 'rowtype' in df.keys():
        df = pdbdf_to_crddf(df)
    write_coordinate_file(df,crdout,comments,crd_spec)

def write_pdb(df,pdbout='output.pdb',comments=''):
    df = df.copy()
    df['atname'] = df['atname'].apply(lambda x: ' ' + x if len(x) == 1 else x)

    if df.resnum.max() > 9999 or df.atnum.max() > 99999:
        print("WARNING: PDB file is > 99,999 atoms or >9,999 residues. Splitting output into multiple files")
        segids = df.segid.unique()
        i = 0
        for segid in segids:
            df_seg = df.loc[df.segid==segid].copy()
            df_seg.atnum-=(df_seg.atnum.min()-1) # Reset atom number
            if df_seg.resnum.max() <= 9999 and df_seg.atnum.max() <= 99999:
                write_coordinate_file(df_seg,pdbout[:-4]+'_'+str(i)+'_'+segid+'.pdb',comments,pdb_spec)
            else:
                j=0
                while not df_seg.empty:
                    df_select = df_seg.loc[df_seg['resnum'] <= 9999]
                    if df_select.atnum.max() > 99999:
                        max_resnum = df_select.loc[df_select['atnum']==100000].resnum.max()
                        df_select = df_select.loc[df_select['resnum'] < max_resnum]
                    write_coordinate_file(df_select,pdbout[:-4]+'_'+str(i)+'_'+segid+'_'+str(j)+'.pdb',comments,pdb_spec)
                    df_seg = df_seg.loc[df_seg['atnum'] > df_select.atnum.max()].copy()
                    df_seg.atnum-=df_select.atnum.max()
                    df_seg.resnum-=df_select.resnum.max()
                    j+=1
            i+=1
    else:
        write_coordinate_file(df,pdbout,comments,pdb_spec)


def prepare_cif_dataframe(df):
    # Select and rename columns to CIF-like format
    df_cif = df.rename(columns={
        'atnum': '_atom_site.id',
        'atname': '_atom_site.label_atom_id',
        'altloc': '_atom_site.label_alt_id',
        'resname': '_atom_site.label_comp_id',
        'chainid': '_atom_site.label_asym_id',
        'resnum': '_atom_site.label_seq_id',
        'x': '_atom_site.Cartn_x',
        'y': '_atom_site.Cartn_y',
        'z': '_atom_site.Cartn_z',
        'occupancy': '_atom_site.occupancy',
        'tfactor': '_atom_site.B_iso_or_equiv',
        'segid': '_atom_site.auth_asym_id',
        'element': '_atom_site.type_symbol',
        'charge': '_atom_site.pdbx_formal_charge'
    })

    # Ensure correct data types and fill in missing data with defaults if needed
    df_cif['_atom_site.occupancy'] = pd.to_numeric(df_cif['_atom_site.occupancy'], errors='coerce').fillna(1.0)
    df_cif['_atom_site.B_iso_or_equiv'] = pd.to_numeric(df_cif['_atom_site.B_iso_or_equiv'], errors='coerce').fillna(0.0)

    # Format necessary columns for CIF conventions, e.g., coordinates as floats with specified precision
    for col in ['_atom_site.Cartn_x', '_atom_site.Cartn_y', '_atom_site.Cartn_z']:
        df_cif[col] = df_cif[col].map('{:.3f}'.format)

    return df_cif


def df_to_structure(df, structure_id="structure"):
    # Initialize a new structure
    structure = Structure.Structure(structure_id)
    model = Model.Model(0)  # Single model (Model ID = 0)
    structure.add(model)

    # Iterate over the DataFrame grouped by chain
    for chain_id, chain_df in df.groupby('chainid'):
        chain = Chain.Chain(chain_id)
        model.add(chain)

        for res_id, res_df in chain_df.groupby(['resnum', 'inscode']):
            # Use residue name, resnum and inscode to identify residues uniquely
            resname = res_df['resname'].iloc[0]
            het_flag = " "  # Regular residues typically use " "
            resnum = res_id[0]
            inscode = res_id[1] if pd.notna(res_id[1]) else " "

            # Initialize a Residue object
            residue = Residue.Residue((het_flag, resnum, inscode), resname, chain_id)
            chain.add(residue)

            # Add atoms to the residue
            for _, atom_row in res_df.iterrows():
                atom = Atom.Atom(
                    atom_row['atname'],
                    [atom_row['x'], atom_row['y'], atom_row['z']],
                    float(atom_row['tfactor']) if pd.notna(atom_row['tfactor']) else 0.0,
                    float(atom_row['occupancy']) if pd.notna(atom_row['occupancy']) else 1.0,
                    float(atom_row['altloc']) if pd.notna(atom_row['altloc']) else 0.0,
                    atom_row['atname'],
                    atom_row['atnum'],
                    element=atom_row['element'].strip()
                )
                residue.add(atom)
    return structure

def write_structure_to_cif(structure, filename):
    # Use Biopython's MMCIFIO class to write the structure as a .cif file
    io = MMCIFIO()
    io.set_structure(structure)
    io.save(filename)
    print(f"Structure written to CIF file: {filename}")


def write_cif(df_cif, filename):
    # CIF files require headers and loop definitions
    with open(filename, 'w') as f:
        f.write("data_cif_generated\n\n")
        f.write("loop_\n")
        f.write(" " + "\n ".join(df_cif.columns) + "\n")

        # Write all columns as CIF format rows (tab-delimited)
        df_cif.to_csv(f, sep="\t", header=False, index=False)
    print(f"CIF file saved as {filename}")


# Reads in PDB files that were split according to the format from the PDB writer.
# This is for instances where there are >99,999 atoms or >9,999 residues
def read_split_pdb(rootname, filepath="."):
    a = os.listdir(filepath)
    files = []
    for b in a:
        if b.startswith(rootname) and b[-4:] == '.pdb':
            files.append(b)
    files.sort()
    df = read_pdb(filepath+"/"+files[0])
    for ufile in files[1:]:
        df1 = read_pdb(filepath+"/"+ufile)
        segid = df1.segid.max()
        if segid in df.segid.unique():
            df1.resnum+=df.loc[df['segid']==segid,'resnum'].max()
        df1.atnum+=df.atnum.max()
        df = pd.concat([df,df1],ignore_index=True)

    return df

def write_coordinate_file(df,filename,comments,column_spec):
    file_type = filename[-3:].upper()
    if "rowtype" in column_spec.keys():
        file_type='PDB'
    elif "ires" in column_spec.keys():
        file_type="CRD"
    #print("file type is ", file_type)
    prev_end = 0
    formatted_columns=[]
    for label in column_spec:
        start,end=column_spec[label][1]
        ljust_bit=column_spec[label][2]
        formatted_column=format_column(df[label],start,end,prev_end,ljust_bit,file_type)
        formatted_columns.append(formatted_column)
        prev_end=end

    #Concatenate all formatted columns
    formatted_df=pd.concat(formatted_columns,axis=1)

    # Combine formatted columns into a single string per row
    formatted_df['combined'] = formatted_df.apply(lambda row: ''.join(row.values), axis=1)

    # Write the formatted DataFrame to a file
    initfile=open(filename,'w')
    write_comments(initfile,file_type,comments.upper())
    if file_type == "PDB":
        initfile.close()
        formatted_df['combined'].to_csv(filename, header=False, index=False, quoting=3,mode='a')
        appendy = open(filename,'a')
        appendy.write('END')
    elif file_type == "CRD":
        maxatom = df.atnum.max()
        initfile.write(str(maxatom).rjust(10)+'  EXT\n')
        initfile.close()
        formatted_df['combined'].to_csv(filename, header=False, index=False, quoting=3,mode='a')
    else:
        initfile.close()

def pdbdf_to_crddf(df):
    df['ires'] = pd.factorize(df[['resnum', 'segid']].apply(tuple, axis=1))[0] + 1
    return df[crd_spec.keys()]

def write_comments(writefile, file_type,comments):
    comm_format = {"PDB": {"width":72, "comment":'REMARK  '},
        "CRD": {"width":130, "comment":'* '}}
    writefile.write(comm_format[file_type]["comment"]+'WRITTEN BY MIKE\'S PDB READER ON '+datetime.now().isoformat()+'\n')
    commlen = len(comments)
    if commlen != 0:
        i=0
        remark_width=comm_format[file_type]["width"]
        while i*remark_width < commlen:
            writefile.write(comm_format[file_type]['comment']+comments[i*remark_width:(i+1)*remark_width]+'\n')
            i+=1

def split_pdb_by_chain(pdbfile):
    df = read_pdb(pdbfile)
    name_root = pdbfile[:-4]
    chains = df['chainid'].unique()
    for chain in chains:
        write_pdb(df.loc[df['chainid']==chain],name_root+'_'+chain+'.pdb')

# Writes an entre folder of PDB files to parquet files.
# Write single parquet file with df.to_parquet('output.parquet',engine='pyarrow')
# Read from parquet using df = pd.read_parquet('output.parquet',engine='pyarrow')
def folder_to_parquet(directory):
    for cur_file in os.listdir(directory):
        if cur_file[-3:] == 'pdb':
            df = read_pdb(directory+'/'+cur_file)
            df.to_parquet('output/'+cur_file[:-3]+'parquet',engine='pyarrow')

"""
test_file = 'DATA/step5_1.crd'
df = read_crd(test_file)
write_crd(df,'testy.crd')
"""

if __name__ == '__main__':
    split_pdb_by_chain(sys.argv[1])
