import pdbio
def read_pdb(pdbfil):
    return pdbio.read_pdb(pdbfil)
def write_pdb(df,pdbfil='output.pdb',comments=''):
    return pdbio.write_pdb(df,pdbfil='output.pdb',comments='')
def folder_to_parquet(directory):
    pdbio.folder_to_parduet(directory)
