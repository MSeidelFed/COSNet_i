#!/usr/env/bin python3
"""
[check_cif_completeness.py]
"""
### IMPORTS
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio import SeqIO, pairwise2
from Bio.pairwise2 import format_alignment
from Bio.PDB import Selection
from pathlib import Path
import argparse

import matplotlib.pyplot as plt
import numpy as np

### VARIABLES
global THREE_TO_ONE
THREE_TO_ONE = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N',
                'ASP': 'D', 'CYS': 'C', 'GLU': 'E',
                'GLN': 'Q', 'GLY': 'G', 'HIS': 'H',
                'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
                'MET': 'M', 'PHE': 'F', 'PRO': 'P',
                'SER': 'S', 'THR': 'T', 'TRP': 'W',
                'TYR': 'Y', 'VAL': 'V'}
### FUNCTIONS
def get_entity_seqs(CIFfile):
    cifdict = MMCIF2Dict(CIFfile)
    parser = MMCIFParser()
    cifstruct = parser.get_structure("ribo", CIFfile)
    chainlst = Selection.unfold_entities(cifstruct, 'C')

    iddict = dict(zip(cifdict["_entity_poly.entity_id"], cifdict["_entity.pdbx_description"]))
    totentries = len(iddict)
    protids = []
    if totentries > 1:
        for i in range(totentries):
            if "polypeptide" in cifdict["_entity_poly.type"][i]:
                protids.append(i)
    else:
        if "polypeptide" in cifdict["_entity_poly.type"]:
            protids.append(int(list(iddict.keys())[0]))

    structseqdict = {}
    canseqdict = {}
    for idx in protids:
        reslist = [res.get_resname() for res in chainlst[idx].get_residues()]
        translated = []
        for res in reslist:
            try:
                translated.append(THREE_TO_ONE[res])
            except:
                translated.append('X')
        structseq = ''.join(translated)
        if len(protids) == 1:
            canonseq = cifdict["_entity_poly.pdbx_seq_one_letter_code_can"] 
        else: 
            canonseq = cifdict["_entity_poly.pdbx_seq_one_letter_code_can"][idx] 
        canonseq = canonseq.strip().split('\n')
        canonseq = ''.join(canonseq)
        structseqdict[idx] = structseq
        canseqdict[idx] = canonseq
    return iddict, structseqdict, canseqdict

def alignseqs(StructSeqD, CanSeqD, IDD):
    pidlist = []
    for key in StructSeqD.keys():
        if len(StructSeqD) == 1:
            print(f'================= Entity {key} name: {IDD[str(key)]} =======================\n')
        else:
            print(f'================= Entity {key} name: {IDD[str(key+1)]} =======================\n')
        print(f'StructSeq: {StructSeqD[key]}\n')
        print(f'CanonSeq: {CanSeqD[key]}\n')
        alns = pairwise2.align.globalms(StructSeqD[key], CanSeqD[key],2, -1, -10,-0.5)
        print(f'Alignment:\n{format_alignment(*alns[0])}')
        percIDval = calc_pid(alns[0][0], alns[0][1])
        pidlist.append(percIDval)
        print(f'Percent ID: {percIDval}')
    return pidlist

def calc_pid(seq1, seq2):
    alnlen = len(seq1)
    identical = 0
    for i in range(0, len(seq1)):
        if seq1[i].isalpha() and seq2[i].isalpha() and seq1[i] == seq2[i]:
            identical += 1
    pid = identical/alnlen
    return pid
    
def parse_cl_arguments():
    parser = argparse.ArgumentParser(usage="python3 %(prog)s [-h] ciffile outpath",
                                     description="Checks for structural completeness of protein entities in CIF")
    parser.add_argument("ciffile", help="CIF file to be checked", type=str)
    parser.add_argument("outpath", help="Output path for results file, default in local dir", type=str)
    args = parser.parse_args()
    return args

def main():
    Args = parse_cl_arguments()
    IDDict, SSDict, CSDict = get_entity_seqs(Args.ciffile)
    print(IDDict)
    print(SSDict)
    print(CSDict)
    name = Args.ciffile.split('.')[0]
    PIDlist = alignseqs(SSDict, CSDict, IDDict)

    with open(f'{Path(Args.outpath)}/{name}_percentage_ids.dat', 'w') as f:
        for pid in PIDlist:
            f.write(str(pid))
            f.write('\n')


### MAIN
if __name__=="__main__":
    main()
