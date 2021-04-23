#!/usr/env/bin/ python3
"""
[split_cif_by_entity.py]

USAGE: python3 split_cif_by_entity.py inputfile pathtooutput

This script is a command-line interface to extract and
save structural data for all entities from a mmCIF file. 
"""

from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from collections import Counter
from pathlib import Path
import argparse
import os

def lower_upper_idx(EntityAtomIDs):
    """
    Gets lower and upper indices (row #s) for each 
    entity in the atomic data block of mmCIF file

    Parameters
    ----------
    EntityAtomIDs: list
        list of entity ids per atom, indicates
        which entity atom belongs to

    Returns
    -------
    idx_dict: dict
        dict with entity atom row ranges
        {entid: (start, stop)}
    """
    ent_instances = dict(Counter(EntityAtomIDs))
    idx_dict = {}
    for i in range(len(ent_instances)):
        if i == 0:
            end = ent_instances['1']-1
            idx_dict[i] = (0,end)
        elif i != 0:
            length = ent_instances[str(i+1)]
            idx_dict[i] = (end+1,end+length)
            end = end+length
    return idx_dict

def writeout_struct_data(ObjName, IDxDict, NamesList, CIFDict, address):
    """
    Takes values from cifdict corresponding to the desired
    structural features
    Saves them to file according to PDB standard format 

    Parameters
    ----------
    ObjName: str
        4-letter cif file identifier
    IDxDict: dict
        dict with entity atom row ranges {entid: (start, stop)}
    NamesList: list
        list of entity names
    CIFDict: Bio.PDB.MMCIF2Dict.MMCIF2Dict
        parsed dict of cif info
    address: str

    Returns
    -------
    None
    """
    print(f'Output PDBs found in: {address}\n')

    atomlist = CIFDict["_atom_site.group_PDB"]
    atomid_num = CIFDict["_atom_site.id"]
    atomid_lbl = CIFDict["_atom_site.label_atom_id"]
    atomtype = CIFDict["_atom_site.type_symbol"]
    atom_comp = CIFDict["_atom_site.label_comp_id"]
    atom_asym = CIFDict["_atom_site.label_asym_id"]
    atom_seqid = CIFDict["_atom_site.label_seq_id"]
    
    atom_Xcoord = CIFDict["_atom_site.Cartn_x"]
    atom_Ycoord = CIFDict["_atom_site.Cartn_y"]
    atom_Zcoord = CIFDict["_atom_site.Cartn_z"]

    atom_occupancy = CIFDict["_atom_site.occupancy"]
    atom_iso = CIFDict["_atom_site.B_iso_or_equiv"]

    for i in range(len(IDxDict)):
        entity_name = NamesList[i].replace(" ","_")
        print(f'{i}: Creating pdb file for {entity_name}\n')
        if 'protein' in entity_name: 
            protid = entity_name.split('_protein_')
            entity_name = f'{protid[1]}_{protid[0]}'
        try:
            pdbfilename = f'{ObjName}_{entity_name}.pdb' 
            structfile = open(f'{address.joinpath(pdbfilename)}', 'w+')
        except:
            print('Error in opening file: Check that naming and path scheme is correct!')
            print('Moving onto next entity.\n')
            continue
        bounds = IDxDict[i]
        k = 1
        for j in range(bounds[0],bounds[1]+1):
            structfile.write('{0:<6}{1:>5} {2:<4} {3:<3}{4:>2}{5:>4}    {6:>8}{7:>8}{8:>8}{9:>6}{10:>6}          {11:>2}'.format(
                            atomlist[j],k,atomid_lbl[j],atom_comp[j],atom_asym[j],
                            atom_seqid[j],atom_Xcoord[j],atom_Ycoord[j],atom_Zcoord[j],
                            atom_occupancy[j],atom_iso[j],atomtype[j]))
            k += 1
            if j != bounds[1]:
                structfile.write('\n')
        structfile.write('\nTER')
        structfile.close()
    return

def check_entityname(entnamelist):
    """
    Checks for non alpha-numeric characters 
    in entity names. Replaces any with underscore.
    Directly modifies input list. 
    Allows: -, (, ), empty space.

    Parameters
    ----------
    entnamelist: list

    Returns 
    -------
    entnamelist: list
    """
    allowedchars = ['-', '(', ')', ' ']
    for i in range(0,len(entnamelist)):
        name = entnamelist[i].split()
        name = ''.join(name)
        if not name.isalnum():
            edited_name = ['_' if (not x.isalnum() and x not in allowedchars) else x for x in list(entnamelist[i])]
            edited_name = ''.join(edited_name)
            entnamelist[i] = edited_name
    return entnamelist
          

def parse_arguments():
    parser = argparse.ArgumentParser(usage="python3 %(prog)s [-h] CIFfile outpath",
                                     description="Splits CIF into separate PDB files per entity")
    parser.add_argument("CIFfile", help="CIF file input to be split", type=str)
    parser.add_argument("outpath",help="Output path", type=str)
    args = parser.parse_args()
    return args
    
def main():
    # get arguments
    Args = parse_arguments()

    # read in data file from command line
    pathtociffile = Path(Args.CIFfile) #pathlib.PosixPath
    obj_name = pathtociffile.stem
    address = Path(Args.outpath)

    if not pathtociffile.exists() or not address.exists():
        raise Exception('Incorrect path specifications in input. Verify!')
    else:
        # save contents as a dictionary
        cif_dict = MMCIF2Dict(pathtociffile) #Bio.PDB.MMCIF2Dict.MMCIF2Dict

        # get lists of identifiers and entities
        ent_id = cif_dict["_entity_poly.entity_id"] #list
        ent_name = cif_dict["_entity.pdbx_description"] #list
        new_ent_name = check_entityname(ent_name) #list
        ent_atom_ids = cif_dict["_atom_site.label_entity_id"] #list

        # get dictionary of indices per entity indicating corresponding rows of the atomic data 
        IDxDict = lower_upper_idx(ent_atom_ids)
        writeout_struct_data(obj_name, IDxDict, new_ent_name, cif_dict, address)
    
## MAIN
if __name__ == "__main__":
    main()
