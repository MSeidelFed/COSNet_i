#!/usr/bin/env python3
"""
[extract_sequence.py]

DESCRIPTION: Extracts sequence corresponding to the documented sequence
in a mmCIF file, also extracts sequence based on structural data in pdb file,
this is to compare sequences for missing residues. 

USAGE: python3 extract_sequence.py PDBfile CIFfile EntityNumber TypeTag 

RETURNS: Fasta file with the two sequences, print in stdout whether sequences
match or not
"""
import os 
import argparse
from pathlib import Path

from Bio.PDB import PDBParser
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
import numpy as np

global THREE_TO_ONE
THREE_TO_ONE = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N',
                'ASP': 'D', 'CYS': 'C', 'GLU': 'E',
                'GLN': 'Q', 'GLY': 'G', 'HIS': 'H',
                'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
                'MET': 'M', 'PHE': 'F', 'PRO': 'P',
                'SER': 'S', 'THR': 'T', 'TRP': 'W',
                'TYR': 'Y', 'VAL': 'V'}
## FUNCTIONS
def three_letter_to_one_letter(chain, code=THREE_TO_ONE):
    """
    Translate a chain of amino acids in three letter code to one letter code.

    Arguments
    ---------
    chain: Bio.PDB.Chain.Chain object or list of Bio.PDB.Residue.Residue
    code:  dict

    Returns
    -------
    translated_chain: str
    """
    translated_chain = []
    for res in list(chain):
        try:
            translated_chain.append(THREE_TO_ONE[res.__dict__['resname']])
        except KeyError:
            warnings.warn(f'Unknown amino acid encountered: {res}, skipping', RuntimeWarning)
    return ''.join(translated_chain)

def get_struct(pathtoPDBfile, PDBname, CIFname):
    """
    Takes in PDB file, parses out structural information
    Returns structure object, original PDB ID of structure, entity name (if 
    file is a portion of the original larger file, say an r protein or rRNA) 

    Note: deals with files of the form: 4v7e_18S_ribosomal_RNA.pdb, 
    or of the form eL13_rpL13.pdb or 4v7e_RACK1.pdb

    Parameters
    ----------
    pathtoPDBfile: str
    PDBname: str
    CIFname: str

    Returns
    -------
    structure: Bio.PDB.Structure.Structure
    struct_id: str
    ent_id: str
    residues: generator
    """
    parser = PDBParser(PERMISSIVE=1)
    struct_id_elems = PDBname.split("_")
    if len(struct_id_elems) == 2:
        if struct_id_elems[0] == CIFname:
            struct_id = struct_id_elems[0]
            ent_id = struct_id_elems[1].split(".")[0]
        else:
            struct_id = CIFname
            ent_id = struct_id_elems[0]
    else:
        struct_id = struct_id_elems[0]
        ent_id = struct_id_elems[1].split(".")[0]
        
    structure = parser.get_structure(ent_id, pathtoPDBfile)
    residues = structure.get_residues()
    atoms = structure.get_atoms()

    return structure, struct_id, ent_id, residues

def get_prot_res_seq_pdb(Residues_obj, tag):
    """
    Gives sequence as taken from the atomic coordinates section

    Parameters
    ----------
    Residues_obj: generator
    tag: str

    Returns
    -------
    seq: str
    """
    ResiList = []
    if (tag == 'protein' or tag == 'prot'):
        print('protein tag')
        for resi in Residues_obj:
            ResiList.append(resi)
        seq = three_letter_to_one_letter(ResiList)
    else:
        for resi in Residues_obj:
            resname = resi.get_resname()
            ResiList.append(resname.strip())
        seq = ''.join(ResiList)
        
    return seq

def get_prot_res_seq_cif(CifDict, EntNum):
    """
    Gives sequence as taken from mmCIF file, entity poly pdbx seq section
    Takes the canonical sequence. 

    Parameters
    ----------
    CifDict: dict
    EntNum: int

    Returns
    -------
    canonseq: str
    """
#    canon_seqlist = CifDict["_entity_poly.pdbx_seq_one_letter_code_can"][EntNum-1].split('\n')
    canon_seqlist = CifDict["_entity_poly.pdbx_seq_one_letter_code"][EntNum-1].split('\n')
    canonseq = ''.join(canon_seqlist)
    return canonseq

def write_out_fasta(fasta_dict, outfile, outpath,flag='w'):
    with open(outfile, flag) as f:
        for k, v in fasta_dict.items():
            f.write(''.join(['>', k, '\n']))
            f.write(''.join([v, '\n']))
    os.rename(outfile, os.path.join(outpath, outfile))

def parse_arguments():
    parser = argparse.ArgumentParser(usage="python3 %(prog)s [-h] pdbfile CIFfile entitynum type outpath",
                                     description="Extracts sequence from CIF file, and also from the corresponding PDB file")
    parser.add_argument("pdbfile", help="PDB file", type=str)
    parser.add_argument("ciffile", help="CIF file", type=str)
    parser.add_argument("entitynum", help="Entity number corresponding to protein, as listed in CIF file (look this up!)", type=int)
    parser.add_argument("type", help="Tag for sequence type: protein, RNA, DNA", type=str)
    parser.add_argument("outputpath", help="Output type", type=str)
    args = parser.parse_args()
    return args

## MAIN
if __name__ == "__main__":
    # parse args
    Args = parse_arguments()
    pathto_pdbfile = Path(Args.pdbfile)
    pdbname = pathto_pdbfile.stem
    pathto_ciffile = Path(Args.ciffile)
    cifname = pathto_ciffile.stem
    entitynum = int(Args.entitynum)
    typetag = str(Args.type) # protein, RNA, DNA
    print(f'----------------------- Type: {typetag}')
    outpath = Args.outputpath
    print(f'Entity Number: {entitynum}')
    #parse out structure and identifiers
    Structure, StructID, EntID, Residues = get_struct(pathto_pdbfile, pdbname, cifname)
    print(f'mmCIF: {StructID} PDB: {EntID}')
    #get sequence from mmCIF file
    cifdict = MMCIF2Dict(pathto_ciffile)
    CIFSeq = get_prot_res_seq_cif(cifdict, entitynum)
    print(f'Sequence from mmCIF: {CIFSeq}')
    #get sequence from pdb file
    PDBSeq = get_prot_res_seq_pdb(Residues, typetag)
    print(f'Sequence from PDB: {PDBSeq}')
  
    #verify sequences are equal
    if CIFSeq != PDBSeq:
        print(f'{StructID} {EntID} {entitynum}: Sequences do not match')

    #print out structural sequence to 
    seqsdict = {}
    seqsdict[f'{EntID}_cifseq'] = CIFSeq
    seqsdict[f'{EntID}_PDBseq'] = PDBSeq
    outfile = f'{pdbname}.fasta'

    #outpath = f'{pdbfile[:-4]}.fasta'
    write_out_fasta(seqsdict, outfile, outpath, 'w')
