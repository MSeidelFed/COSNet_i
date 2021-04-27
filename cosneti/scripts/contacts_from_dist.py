#!/usr/env/bin/ python3
"""
[contacts_from_dist.py]

This script is a command-line interface to
find contacting residues under a distance threshold and 
print out summary info to the console.
Any contacts are saved in separate data files.

It can take either one distance matrix csv file as input, or 
a text file listing distance matrix csv files, one per line. 
Paths to input and output must be specified.

USAGE: python3 contacts_from_dist.py threshold inputfile pathtocsvs pathtopdbs outpath
"""

import os
import glob
import argparse
from pathlib import Path

from Bio.PDB import PDBParser
import numpy as np

def get_list_of_files(infile):
    """
    Makes list of file names 

    Parameters
    ----------
    infile: pathlib.PosixPath

    Returns
    -------
    filelist: list
        list of file names
    """
    filelist = []
    with open(infile) as f:
        for line in f:
            filelist.append(line.strip())

    return filelist

def read_dist_mtx(inputfile):
    """Reads in csv file as np array

    Parameters
    ----------
    inputfile: str
        path to inputfile

    Returns
    -------
    dist_mtx: numpy.ndarray
        matrix of inter-residue distances
    """
    dist_mtx = np.loadtxt(inputfile, dtype=float, delimiter=',')
    return dist_mtx

def find_contact_idx(dist_mtx, ent1name, ent2name, noncontacts, contacts, threshold, resultspath):
    """Find indices for residues with distances below contact threshold.
       Saves existing contacts to a file, otherwise adds proteins to non_contact file.
       Returns two lists of tuples, where contacts have a tuple within a tuple!

    Parameters
    ----------
    dist_mtx: numpy.ndarray
    ent1name: str
        name of first entity (usually prot1)
    ent2name: str
        name of first entity (usually prot2)
    noncontacts: list
        list of tuples of prots without contacts, e.g. [('prot1', 'prot2'),...]
    contacts: list
        list of tuples of prots with contacts and number of contacts e.g. [(('prot2', 'prot3'), 3),...]
    threshold: int
        distance threshold in Å for contacts
    resultspath: pathlib.PosixPath

    Returns
    -------
    noncontacts: list
    contacts: list
    """
    pair = (ent1name, ent2name)
    if np.all(dist_mtx >= int(threshold)):
        noncontacts.append(pair)
    else:
        idx = np.where(dist_mtx <= int(threshold))
        vals = dist_mtx[idx]
        num_contacts = len(vals)
        contacts.append((pair, num_contacts))
        array1 = idx[0].reshape(idx[0].shape[0],1)
        array2 = idx[1].reshape(idx[1].shape[0],1)
        idx_zipped = np.hstack((array1, array2))
        contact_array = np.hstack((idx_zipped, vals.reshape(vals.shape[0],1)))
        outname = f'contacts_t{threshold}_{ent1name}_{ent2name}.dat'
        outfile = os.path.join(resultspath, outname)
        np.savetxt(outfile, contact_array, delimiter = ' ',fmt='%i %i %f')
    return noncontacts, contacts

def get_structure(protname, datapath):
    """
    Parses in PDB structure, asks for user input
    if multiple options exist.

    Parameters
    ----------
    protname: str
        name for your protein, used to search for pdbs
    datapath: pathlib.PosixPath
        path to dir with pdb files

    Returns
    -------
    structure: Bio.PDB.Structure.Structure
    """
    parser = PDBParser(PERMISSIVE=1)
    pdblist = glob.glob(f'{datapath}/*{protname}[_.]*.pdb') ###
    
    if len(pdblist) > 1:
        print(pdblist)
        fileyouwantdict = {}
        count = 1
     
        print('You have these options for pdb files: ')
        for pdb in pdblist:
            print(f'{count}) {os.path.basename(pdb)}')
            fileyouwantdict[count] = pdb
            count += 1
            
        idx = int(input('Enter number corresponding to correct pdb file: \n'))
        #pdbfile = glob.glob(f'{datapath}/{protname}_*_reidx.pdb')
        #print(os.path.basename(pdbfile[0]))
        pdbfile = fileyouwantdict[idx]
    else:
        pdbfile = pdblist [0]
    structure = parser.get_structure(protname, pdbfile)
    return structure 
    
def count_atoms_under_thresh(struct1, struct2, ent1, ent2, threshold, resultspath):
    """
    Grabs any contact files between ent1 and ent2.
    Compares two structures.
    Counts number of atoms per residue that are falling under the threshold value.
    Also counts number of connections ('links' between atoms) under this threshold.

    Parameters
    ----------
    struct1: Bio.PDB.Structure.Structure
    struct2: Bio.PDB.Structure.Structure
    ent1: str
    ent2: str
    threshold: int
    resultspath: pathlib.PosixPath

    Returns
    -------
    per_res_connections_count: list
    per_res_atoms_count: list
    """
    contactfile = glob.glob(f'{resultspath}/contacts_t{threshold}_{ent1}_{ent2}[._]*')   ### 
    if not contactfile:
        contactfile = glob.glob(f'{resultspath}/contacts_t{threshold}_{ent2}_{ent1}[._]*')    ###
    print('Contact File: ')
    print(os.path.basename(contactfile[0]))
    print('Contact Residues: ')
    filepath = contactfile[0]
    resi_contactlist = []
    with open(filepath, 'r') as f:
        for line in f.readlines():
            contents = line.strip().split()
            resi1 = int(contents[0])+1 # numbering in np arrays starts from 0
            resi2 = int(contents[1])+1 # so add one to match reindexed pdbs
            resi_contactlist.append((resi1, resi2))

    struct1_id = get_model_chain_ids(struct1)
    print(type(struct1_id))
    print(f'(model#, chainID) for {ent1}: {struct1_id}')
    struct2_id = get_model_chain_ids(struct2)
    print(f'(model#, chainID) for {ent2}: {struct2_id}')

    per_res_connections_count = [] # lists per residue pair connections
    per_res_atoms_count = []

    for resipair in resi_contactlist:
        per_res_count = 0
        resi1 = struct1[struct1_id[0]][struct1_id[1]][resipair[0]]
        resi2 = struct2[struct2_id[0]][struct2_id[1]][resipair[1]]
        print(f'------ {resipair}: {resi1.get_resname()} {resi2.get_resname()}')
        atomdict1 = {}
        atomdict2 = {}
        for atom1 in resi1:
            per_atom_count = 0
            atom1_name = f'{atom1.get_serial_number()}_{atom1.get_name()}'
            #print(f'+++++++ Atom 1 {atom1_name}')
            coords1 = np.array(atom1.get_coord())
            atomdict1[atom1_name] = 0
            for atom2 in resi2:
                atom2_name = f'{atom2.get_serial_number()}_{atom2.get_name()}'
                #print(f'++++ Atom 2 {atom2_name}')
                coords2 = np.array(atom2.get_coord())
                dist = np.linalg.norm(coords1-coords2)
                if not atom2_name in atomdict2:
                    atomdict2[atom2_name] = 0
                if dist.item() <= int(threshold):
                    #print(f'Distance Under threshold: {dist}')
                    per_res_count += 1
                    per_atom_count += 1
                    opp_count = atomdict2[atom2_name] + 1
                    atomdict2[atom2_name] = atomdict2[atom2_name] + 1
            atomdict1[atom1_name] = per_atom_count
        print(f'Atom 1 Dict: {atomdict1}')
        print(f'Atom 2 Dict: {atomdict2}')
        print(f'Connections Per Residue: {per_res_count}')
        per_res_connections_count.append(per_res_count)
        print(f'List of connections per residue: {per_res_connections_count}')
        ct = 0
        for val in atomdict1.values():
            if val != 0:
                ct += 1
        for val in atomdict2.values():
            if val != 0:
                ct += 1
        per_res_atoms_count.append(ct) 

    return per_res_connections_count, per_res_atoms_count

def get_model_chain_ids(structure):
    """ 
    Extracts model and chain ids from structure object.

    Parameters
    ----------
    structure: Bio.PDB.Structure.Structure

    Returns
    -------
    struct_id: tuple
        tuple of (modelid, chainid)
    """
    modelnum = []
    chainid = []
    for model in structure:
        modelnum.append(model.get_id())
        for chain in model:
            chainid.append(chain.get_id())

    struct_id = (modelnum[0], chainid[0])
    return struct_id
            
def print_summary_output(contact_list, summed_count_list, filename, resultspath):
    """ 
    Writes contact information to output files.

    Parameters
    ----------
    contact_list: list
    summed_count_list: list
    filename: str
    resultspath: pathlib.PosixPath
        path to dir to write out files
    """
    if contact_list:
        outfile = os.path.join(resultspath, filename)
        contactfile = open(outfile, 'w+')
        for i in range(len(contact_list)):
            if filename.strip().split('_')[1]=='noncontacts':
                contactfile.write(f'{contact_list[i][0]} {contact_list[i][1]}\n')
            else: 
                contactfile.write(f'{contact_list[i][0][0]} {contact_list[i][0][1]} {contact_list[i][1]} {summed_count_list[i][0]} {summed_count_list[i][1]}\n')
        contactfile.close()
    else:
        return

def append_count_info(per_res_count_list, total_atom_count_list, ent1name, ent2name, threshold, resultspath):
    """
    Writes atom and connection count info to 'full' version of contact files.

    Parameters
    ----------
    per_res_count_list: list
    total_atom_count_list: list
    ent1name: str
    ent2name: str
    threshold: int
    resultspath: pathlib.PosixPath
    """
    infile = os.path.join(resultspath, f'contacts_t{threshold}_{ent1name}_{ent2name}.dat')
    outfile = os.path.join(resultspath, f'contacts_t{threshold}_{ent1name}_{ent2name}_full.dat')
    data = open(infile, 'r').readlines()
    output = []
    for i in range(len(data)):
        line = data[i].strip().split()
        res1 = int(line[0])+1
        res2 = int(line[1])+1
        dist = line[2]
        output.append('{} {} {} {} {}'.format(res1, res2, dist, per_res_count_list[i], total_atom_count_list[i]))
    f = open(outfile, 'w')
    f.write('\n'.join(output))
    f.close()
    return

def parse_arguments():
    parser = argparse.ArgumentParser(usage="python3 %(prog)s [-h] threshold inputfile pathtocsvs pathtopdbs outpath",
                                     description="Finds residues in contact under a defined distance threshold")
    parser.add_argument("threshold", help="Value (angstroms) representing distance under which residues are in contact (normally 8)", type=int)
    parser.add_argument("inputfile", help="One CSV or a text file with list of CSV (distance matrices)", type=str)
    parser.add_argument("pathtocsvs", help="Path to where dist matrix CSVs are located", type=str)
    parser.add_argument("pathtopdbs", help="Path to where PDB structure files are located", type=str)
    parser.add_argument("outpath", help="Output path for contact files", type=str)
    args = parser.parse_args()
    return args
 
## MAIN
if __name__ == "__main__":
    # get arguments
    Args = parse_arguments()
    threshold = int(Args.threshold)
    pathtoinfile = Path(Args.inputfile)
    infile = pathtoinfile.name

    pathto_csvs = Path(Args.pathtocsvs)
    pathto_pdbs = Path(Args.pathtopdbs)
    output_path = Path(Args.outpath)

    # is input file one csv or a file listing csvs? 
    if pathtoinfile.suffix != '.csv':
        FilesList = get_list_of_files(pathtoinfile)
        print(FilesList)
    else:
        FilesList = [infile]
        print(FilesList)


    NonContactList = []
    ContactList = []
    SummedCountList = []
    count = 0
    for afile in FilesList:
        print('+++++++++++++++++++++++++++++++++++++')
        print(f'Distance File: {afile}')
        # your dist mtx and pdb files need to have X_x.pdb, Y_y.pdb and dist_mtx_X_Y.csv format!
        fileparts = afile.split('_')
        ent1 = fileparts[2]
        ent2 = fileparts[3].split('.')[0]
        print(f'Comparing: {ent1} {ent2}')
        # Load data
        filename = pathto_csvs.joinpath(afile)
        print(filename)
        try:
            DistMtx = read_dist_mtx(filename)
        except: 
            raise FileNotFoundError(f'WARNING: {filename} not found in {pathto_csvs}')
            continue
        # Find contacts and output to respective files
        NonConList, ConList = find_contact_idx(DistMtx, ent1, ent2, NonContactList, ContactList, threshold, output_path)
        count +=1
        print(f'Running on protein pair {count}')
        pair = (ent1, ent2)
        if not pair in NonConList:
            print('CONTACTS EXIST')
            # Parse out structures
            StructEnt1 = get_structure(ent1, pathto_pdbs)
            StructEnt2 = get_structure(ent2, pathto_pdbs)
            # Calculate atoms in contacts
            PerResCountList, PerResAtomCountList = count_atoms_under_thresh(StructEnt1, StructEnt2, ent1, ent2, threshold, output_path)
            print(PerResCountList)
            print(type(PerResCountList))
            print(PerResAtomCountList)
            print(type(PerResAtomCountList))
            # Print out where you are
            SummedCountList.append((sum(PerResCountList), sum(PerResAtomCountList)))
            append_count_info(PerResCountList, PerResAtomCountList, ent1, ent2, threshold, output_path)
        else:
            print('NO CONTACTS UNDER THRESHOLD')
        print('=====================================\n')

    print_summary_output(NonConList, SummedCountList, f'summary_noncontacts_t{threshold}.txt', output_path)
    print_summary_output(ConList, SummedCountList, f'summary_contacts_t{threshold}.txt', output_path)
