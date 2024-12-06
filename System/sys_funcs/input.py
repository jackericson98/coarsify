from os import path
from System.sys_funcs.calcs import get_radius
from System.sys_objs.chain import Sol, Chain
from System.sys_objs.atom import make_atom
from System.sys_objs.residue import Residue
import numpy as np
from pandas import DataFrame
from System.sys_funcs.calcs import calc_dist


def fix_sol(residue):
    """
    Reorganizes atoms within a residue to ensure oxygen and hydrogen atoms are correctly grouped into water molecules.
    Handles cases where there are extra or missing hydrogens.

    Parameters:
    residue (Residue): The residue object containing atoms to be organized.

    Returns:
    list: A list of 'Residue' objects with correctly assigned atoms.
    """

    # Initialize containers for oxygen and hydrogen atoms
    oxy_res = []
    hydrogens = []

    # Separate oxygen and hydrogen atoms into different lists
    for atom in residue.atoms:
        if atom['element'].lower() == 'o':
            # Create a new residue for each oxygen atom
            oxy_res.append(Residue(sys=residue.sys, atoms=[atom], name=atom['residue'],
                                   sequence=atom['res_seq'], chain=atom['chn']))
        elif atom['element'].lower() == 'h':
            hydrogens.append(atom)

    # Assign hydrogens to the nearest oxygen atom to form water molecules
    for h in hydrogens:
        closest_res, min_dist = None, np.inf
        for res in oxy_res:
            dist = calc_dist(res.atoms[0]['loc'], h['loc'])
            if dist < min_dist:
                min_dist = dist
                closest_res = res
        if closest_res:
            closest_res.atoms.append(h)

    # Check the integrity of newly formed residues
    good_resids = []
    incomplete_resids = []
    for res in oxy_res:
        if len(res.atoms) == 3:  # A complete water molecule has 3 atoms: O and 2 H
            good_resids.append(res)
        else:
            incomplete_resids.append(res)

    # Attempt to correct incomplete residues
    for res in incomplete_resids:
        if len(res.atoms) < 3:
            # This block tries to find hydrogens that can be moved to this residue
            for h in hydrogens:
                dist = calc_dist(res.atoms[0]['loc'], h['loc'])
                if dist < 1.5:  # Assumed maximum bond length for O-H
                    res.atoms.append(h)
                    hydrogens.remove(h)
                if len(res.atoms) == 3:
                    break
        if len(res.atoms) != 3:  # Still incomplete after trying to add hydrogens
            good_resids.append(res)  # Optionally, handle still incomplete residues separately

    return good_resids


def read_pdb(sys):
    # Check to see if the file is provided and use the base file if not
    file = sys.base_file

    # Get the file information and make sure to close the file when done
    with open(file, 'r') as f:
        my_file = f.readlines()

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
        word = line[:6].lower().strip()
        # Check to see if the line is an atom line. If the line is not an atom line store the other data
        if line and word not in {'atom', 'hetatm'}:
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
        # if res_name in resids and atom['chain'] == 'Z' and len(resids[res_name].atoms) >= 3:
        #     my_res = Residue(sys=sys, atoms=[atom], name=atom['residue'], sequence=atom['res_seq'], chain=atom['chn'])
        #     atom['chn'].residues.append(my_res)
        #     resids[res_name] = my_res
        #     sys.residues.append(my_res)
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

        # Assign the mass
        atom['mass'] = sys.masses[atom['element'].lower()]

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
    adjusted_residues = []
    for res in sys.residues:
        if res.name == 'SOL' and len(res.atoms) > 3:
            adjusted_residues += fix_sol(res)
        else:
            adjusted_residues.append(res)
    sys.residues = adjusted_residues
    for res in sys.residues:
        if res.name == 'ION':
            res.elem_col = res.atoms[0]['name']
        elif res.name in res_colors:
            res.elem_col = 'Al'
        else:
            res.elem_col = 'Ti'

    # Set the atoms and the data
    sys.atoms, sys.settings = DataFrame(atoms), data