#!/usr/env/bin/python3
"""
[contacts_from_dist.py]: Checks if values in a distance matrix are under
a certain threshold value, if such values exist, prints them out into a separate file
with indices and distances

USAGE: python3 contacts_from_dist.py <filewithcsvs> <thresholdvalue>
"""
## IMPORTS
from Bio.PDB import PDBParser
from sys import argv
import numpy as np
import os
import glob
## FUNCTIONS
def get_files_list(inputfile):
    """ """
    files_list = []
    with open(inputfile, 'r') as f:
        for line in f:
            files_list.append(line.strip())

    return files_list

def read_dist_mtx(inputfile):
    """Reads in csv file as np array"""
    dist_mtx = np.loadtxt(inputfile, dtype=float, delimiter=',')
    return dist_mtx

def find_contact_idx(dist_mtx, ent1name, ent2name, noncontacts, contacts, threshold, resultspath):
    """Find indices and distances below threshold
       outputs to a file if contacts exist otherwise
       adds to non_contact file
       Returns two lists of tuples, where contacts have a tuple within a tuple!
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
    """
    parser = PDBParser(PERMISSIVE=1)
    pdblist = glob.glob(f'{datapath}/*{protname}*.pdb')
    
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
    """Takes in list of interacting residues and counts number of atoms per residue
       that are falling under the threshold value
    """
    contactfile = glob.glob(f'{resultspath}/contacts_t{threshold}*{ent1}**{ent2}*')    
    if not contactfile:
        contactfile = glob.glob(f'{resultspath}/contacts_t{threshold}*{ent2}**{ent1}*')
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
        #print(f'Atom 1 Dict: {atomdict1}')
        #print(f'Atom 2 Dict: {atomdict2}')
        #print(f'Connections Per Residue: {per_res_count}')
        per_res_connections_count.append(per_res_count)
        #print(f'List of connections per residue: {per_res_connections_count}')
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
    """ """
    modelnum = []
    chainid = []
    for model in structure:
        modelnum.append(model.get_id())
        for chain in model:
            chainid.append(chain.get_id())

    struct_id = (modelnum[0], chainid[0])
    return struct_id
            
def print_summary_output(contact_list, summed_count_list, filename, resultspath):
    """ """
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
    """ """
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
## MAIN
if __name__ == "__main__":

    # print usage
    if len(argv) != 6:
        print('#############################################################################################################\n')
        print('  Usage: python3 contacts_from_dist.py [threshold] [pathtoinput] [pathtocsvs] [pathtopdbs] [pathtooutput]    \n')
        print('          The threshold (in angstroms) is the distance under which two residues are in contact               \n')
        print('#############################################################################################################\n')
    else:
        # read in args
        threshold = argv[1]
        pathto_infile = argv[2]
        infile = os.path.basename(pathto_infile)

        pathto_csvs = argv[3]
        pathto_pdbs = argv[4]
        output_path = argv[5]

        # is input file one csv or a file listing csvs? 
        if infile[-4:] != '.csv':
            FilesList = get_files_list(pathto_infile)
        else:
            FilesList = [infile]

        NonContactList = []
        ContactList = []
        SummedCountList = []
        count = 0
        for afile in FilesList:
            print('+++++++++++++++++++++++++++++++++++++')
            print(f'Distance File: {afile}')
            fileparts = afile.split('_')
            ent1 = fileparts[2]
            ent2 = fileparts[3].split('.')[0]
            print(f'Comparing: {ent1} {ent2}')
            # Load data
            filename = os.path.join(pathto_csvs, afile)
            print(filename)
            DistMtx = read_dist_mtx(filename)
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
                # Print out where you are
                SummedCountList.append((sum(PerResCountList), sum(PerResAtomCountList)))
                append_count_info(PerResCountList, PerResAtomCountList, ent1, ent2, threshold, output_path)
            else:
                print('NO CONTACTS UNDER THRESHOLD')
            print('=====================================\n')

        print_summary_output(NonConList, SummedCountList, f'summary_noncontacts_t{threshold}.txt', output_path)
        print_summary_output(ConList, SummedCountList, f'summary_contacts_t{threshold}.txt', output_path)
