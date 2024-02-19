from os import path
from System.sys_funcs.calcs import get_radius
from System.sys_objs.chain import Sol, Chain
from System.sys_objs.atom import make_atom
from System.sys_objs.residue import Residue
import numpy as np
from pandas import DataFrame


def read_pdb(sys):
    # Check to see if the file is provided and use the base file if not
    file = sys.base_file

    # Get the file information and make sure to close the file when done
    with open(file, 'r') as f:
        my_file = f.readlines()
    # Add the system name and reset the atoms and data lists
    sys.name = path.basename(sys.base_file)[:-4] + '_coarse'
    # Set up the atom and the data lists
    atoms, data = [], []
    sys.chains, sys.residues = [], []
    chains, resids = {}, {}
    # Go through each line in the file and check if the first word is the word we are looking for
    reset_checker = 0
    for i in range(len(my_file)):
        # Check to make sure the line isn't empty
        if len(my_file[i]) == 0:
            continue
        # Pull the file line and first word
        line = my_file[i]
        word = line[:4].lower()
        # Check to see if the line is an atom line. If the line is not an atom line store the other data
        if line and word != 'atom':
            data.append(my_file[i].split())
            continue
        # Check for the "m" situation
        if line[76:78] == ' M':
            continue
        # Get the residue sequence of the atom
        res_seq = int(line[22:26])
        if line[22:26] == '    ':
            res_seq = 0
        # Create the atom
        atom = make_atom(location=np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])]),
                         system=sys, element=line[76:78].strip(), res_seq=res_seq, name=line[12:16].strip(),
                         seg_id=line[72:76], index=line[6:11], chain=line[21], residue=line[17:20].strip())
        # Collect the chain type for each atom
        if atom['chain'] == ' ':
            if atom['residue'].lower() in {'sol', 'hoh', 'sod', 'w'}:
                atom['chain'] = 'Z'
            elif atom['residue'].lower() in {'cl', 'mg', 'na', 'k', 'ion'} and 'Z' in chains:
                atom['chain'] = 'X'
            else:
                atom['chain'] = 'A'
        # Create the chain and residue dictionaries
        res_name = atom['chain'] + '_' + atom['residue'] + '_' + str(atom['res_seq']) + '_' + str(reset_checker)
        # If the chain has been made before
        if atom['chain'] in chains:
            # Get the chain from the dictionary and add the atom
            my_chn = chains[atom['chain']]
            my_chn.add_atom(atom['num'])
            atom['chn'] = my_chn
        # Create the chain
        else:
            # If the chain is the sol chain
            if atom['chain'] == 'Z':
                my_chn = Sol(atoms=[atom['num']], residues=[], name=atom['chain'], sys=sys)
                sys.sol = my_chn
            # If the chain is not sol create a regular chain object
            else:
                my_chn = Chain(atoms=[atom['num']], residues=[], name=atom['chain'], sys=sys)
                sys.chains.append(my_chn)
            # Set the chain in the dictionary and give the atom it's chain
            chains[atom['chain']] = my_chn
            atom['chn'] = my_chn

        # Assign the atoms and create the residues
        if res_name in resids:
            my_res = resids[res_name]
            my_res.atoms.append(atom)
        else:
            my_res = Residue(sys=sys, atoms=[atom], name=atom['residue'], sequence=atom['res_seq'], chain=atom['chn'])
            atom['chn'].residues.append(my_res)
            resids[res_name] = my_res
            sys.residues.append(my_res)
        # Assign the residue to the atom
        atom['res'] = my_res

        # Assign the radius
        atom['rad'] = get_radius(atom)
        # If the residue numbers roll over reset the name of the residue to distinguish between the residues
        atoms.append(atom)
        if res_seq == 9999:
            reset_checker += 1
    # Set the colors for the residues based off the default colors for set elements
    res_colors = {'ALA': 'H', 'ARG': "He", 'ASN': 'Li', 'ASP': 'Be', 'ASX': 'B', 'CYS': 'C', 'GLN': 'F', 'GLU': 'O',
                  'GLX': 'S', 'GLY': 'Cl', 'HIS': 'Ar', 'ILE': 'Na', 'LEU': 'Mg', 'LYS': 'Mg', 'MET': 'Al',
                  'PHE': 'Si', 'PRO': 'P', 'SER': 'S', 'THR': 'Cl', 'TRP': 'Ar', 'TYR': 'K', 'VAL': 'Ca',
                  'SOL': 'Ti', 'DA': 'N', 'DC': 'O', 'DG': 'F', 'DT': 'S', 'NA': 'NA', 'CL': 'CL', 'MG': 'MG',
                  'K': 'K'}
    # Set the residues' colors and let the default go to Titanium (Grey)
    for res in sys.residues:
        if res.name == 'ION':
            res.elem_col = res.atoms[0]['name']
        elif res.name in res_colors:
            res.elem_col = res_colors[res.name]
        else:
            res.elem_col = 'Ti'
    # Set the atoms and the data
    sys.atoms, sys.data = DataFrame(atoms), data