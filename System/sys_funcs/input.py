from os import path
from System.sys_funcs.calcs import get_radius
from System.sys_objs.chain import Sol, Chain
from System.sys_objs.atom import make_atom
from System.sys_objs.residue import Residue
import numpy as np
from pandas import DataFrame
from System.sys_funcs.calcs import calc_dist
from chemistry_interpreter import my_masses


def fix_sol(sys, residue):
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
    for a in residue.atoms:
        atom = sys.atoms.iloc[a]

        if atom['element'].lower() == 'o':
            # Create a new residue for each oxygen atom
            oxy_res.append(Residue(sys=residue.sys, atoms=[a], name='SOL',
                                   sequence=atom['res_seq'], chain=atom['chn']))
            # Add the residue to the atom
            sys.atoms.loc[a, 'res'] = oxy_res[-1]
        elif atom['element'].lower() == 'h':
            hydrogens.append(atom['num'])

    # Assign hydrogens to the nearest oxygen atom to form water molecules
    for h in hydrogens:
        closest_res, min_dist = None, np.inf
        for res in oxy_res:
            dist = calc_dist(sys.atoms['loc'][res.atoms[0]], sys.atoms['loc'][h])
            if dist < min_dist:
                min_dist = dist
                closest_res = res
        if closest_res and min_dist < 2.5:
            closest_res.atoms.append(h)
            # Add the residue to the atom
            sys.atoms.loc[h, 'res'] = closest_res
            hydrogens.remove(h)

    # Check the integrity of newly formed residues
    good_resids = []
    incomplete_resids = []
    for res in oxy_res:
        if len(res.atoms) == 3:  # A complete water molecule has 3 atoms: O and 2 H
            good_resids.append(res)
            for a in res.atoms:
                sys.atoms.loc[a, 'res'] = res
        else:
            incomplete_resids.append(res)

    # Attempt to correct incomplete residues
    for res in incomplete_resids:
        if len(res.atoms) < 3:
            # This block tries to find hydrogens that can be moved to this residue
            for h in hydrogens:
                dist = calc_dist(sys.atoms['loc'][res.atoms[0]], sys.atoms['loc'][h])
                if dist < 2.5:  # Assumed maximum bond length for O-H
                    res.atoms.append(h)
                    # Add the residue to the atom
                    sys.atoms.loc[h, 'res'] = res
                    hydrogens.remove(h)
                if len(res.atoms) == 3:
                    break
            # print([(sys.atoms['name'][_], sys.atoms['res_seq'][_], sys.atoms['loc'][_][0]) for _ in res.atoms])
        good_resids.append(res)

    # Last give the hydrogens their own residues
    for h in hydrogens:
        # Get the hydrogen atoms
        hy = sys.atoms.iloc[h]
        # Create the residue
        good_resids.append(Residue(sys=residue.sys, atoms=[h], name='SOL', sequence=hy['res_seq'], chain=hy['chn'],
                                   element='H'))
        # Add the residue to the atom
        sys.atoms.loc[h, 'res'] = good_resids[-1]

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
    reset_checker, atom_count = 0, 0
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

        name = line[12:16]
        res_seq = line[22:26]
        if line[22:26] == '    ':
            res_seq = 0
        # If no chain is specified, set the chain to 'None'
        res_str, chain_str = line[17:20].strip(), line[21]

        # Create the atom
        mass = my_masses[line[76:78].strip().lower()]
        # Create the atom
        atom = make_atom(location=np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])]),
                         system=sys,
                         element=line[76:78].strip(), res_seq=int(res_seq), res_name=res_str, chn_name=chain_str,
                         name=name.strip(), seg_id=line[72:76], index=atom_count, set_index=int(line[6:11]), mass=mass)
        atom_count += 1

        if chain_str == ' ':
            if res_str.lower() in {'sol', 'hoh', 'sod', 'out', 'cl', 'mg', 'na', 'k', 'ion', 'cla'}:
                chain_str = 'SOL'
            else:
                chain_str = 'A'

        # Create the chain and residue dictionaries
        res_name, chn_name = chain_str + '_' + line[17:20] + str(atom['res_seq']) + '_' + str(reset_checker), chain_str
        # If the chain has been made before
        if chn_name in chains:
            # Get the chain from the dictionary and add the atom
            my_chn = chains[chn_name]
            my_chn.add_atom(atom['num'])
            atom['chn'] = my_chn
        # Create the chain
        else:
            # If the chain is the sol chain
            if res_str.lower() in {'sol', 'hoh', 'sod', 'out', 'cl', 'mg', 'na', 'k', 'ion',
                                   'cla'} or chn_name == 'SOL':
                my_chn = Sol(atoms=[atom['num']], residues=[], name=chn_name, sys=sys)
                sys.sol = my_chn
            # If the chain is not sol create a regular chain object
            else:
                my_chn = Chain(atoms=[atom['num']], residues=[], name=chn_name, sys=sys)
                sys.chains.append(my_chn)
            # Set the chain in the dictionary and give the atom it's chain
            chains[chn_name] = my_chn
            atom['chn'] = my_chn

        # Assign the atoms and create the residues
        # if res_name in resids and atom['chain'] == 'Z' and len(resids[res_name].atoms) >= 3:
        #     my_res = Residue(sys=sys, atoms=[atom], name=atom['residue'], sequence=atom['res_seq'], chain=atom['chn'])
        #     atom['chn'].residues.append(my_res)
        #     resids[res_name] = my_res
        #     sys.residues.append(my_res)
        # Assign the atoms and create the residues
        if res_name in resids:
            my_res = resids[res_name]
            my_res.atoms.append(atom['num'])
        else:
            my_res = Residue(sys=sys, atoms=[atom['num']], name=res_str, sequence=atom['res_seq'],
                             chain=atom['chn'])
            resids[res_name] = my_res
            if res_str.lower() in {'sol', 'hoh', 'sod', 'out', 'cl', 'mg', 'na', 'k', 'ion',
                                   'cla'} or chain_str == 'SOL':
                sys.sol.residues.append(my_res)
            else:
                sys.residues.append(my_res)
                atom['chn'].residues.append(my_res)
        # Assign the residue to the atom
        atom['res'] = my_res

        # Add the atom
        atoms.append(atom)
        # If the residue numbers roll over reset the name of the residue to distinguish between the residues
        if res_seq == 9999:
            reset_checker += 1
    # Set the colors for the residues based off the default colors for set elements
    res_colors = {'ALA': 'H', 'ARG': "He", 'ASN': 'Li', 'ASP': 'Be', 'ASX': 'B', 'CYS': 'C', 'GLN': 'F', 'GLU': 'O',
                  'GLX': 'S', 'GLY': 'Cl', 'HIS': 'Ar', 'ILE': 'Na', 'LEU': 'Mg', 'LYS': 'Mg', 'MET': 'Al',
                  'PHE': 'Si', 'PRO': 'P', 'SER': 'S', 'THR': 'Cl', 'TRP': 'Ar', 'TYR': 'K', 'VAL': 'Ca',
                  'SOL': 'Ti', 'DA': 'N', 'DC': 'O', 'DG': 'F', 'DT': 'S', 'NA': 'NA', 'CL': 'CL', 'MG': 'MG',
                  'K': 'K'}

    # Set the atoms and the data
    sys.atoms, sys.settings = DataFrame(atoms), data
    # Adjust the SOL residues
    adjusted_residues = []
    for res in sys.sol.residues:
        if len(res.atoms) > 3:
            adjusted_residues += fix_sol(sys, res)
        else:
            adjusted_residues.append(res)
    # Add the sys.sol residues
    sys.sol.residues = adjusted_residues
    # Instantiate the separated residues
    new_residues = []
    # Do a final check
    for res in sys.sol.residues:
        # Check that they have neighboring indices
        indices = [sys.atoms['index'][_] for _ in res.atoms]
        # Check for differences greater than 2
        if any([_ > 2 for _ in [abs(indices[0] - __) for __ in indices]]):
            # Get the locations
            locs = [sys.atoms['loc'][_] for _ in res.atoms]
            # Scrutinize their relative distances
            if any([calc_dist(locs[0], _) > 5 for _ in locs]):
                # Remove the residue and split its atoms into separate atom residues.
                sys.sol.residues.remove(res)
                # Split the residue into separate residues
                for a in res.atoms:
                    # Get the atom from the dataframe
                    atom = sys.atoms.iloc[a]
                    # Create the new residue
                    new_residues.append(Residue(sys=res.sys, atoms=[a], name=atom['name'], sequence=atom['res_seq'],
                                                chain=atom['chn'], element=atom['element']))
                    # Add the residue to the atom
                    sys.atoms.loc[a, 'res'] = new_residues[-1]
                # No need to split again
                continue
        # Check for the atom types
        a_types = [sys.atoms['element'][_] for _ in res.atoms]
        # Check that there arent more than one o and 2 h
        if len([_ for _ in a_types if 'o' in _.lower()]) > 1 or len([_ for _ in a_types if 'h' in _.lower()]) > 2:
            # Remove the residue and split its atoms into separate atom residues.
            sys.sol.residues.remove(res)
            # Split the residue into separate residues
            for a in res.atoms:
                # Get the atom from the dataframe
                atom = sys.atoms.iloc[a]
                # Create the new residue
                new_residues.append(Residue(sys=res.sys, atoms=[a], name=atom['name'], sequence=atom['res_seq'],
                                            chain=atom['chn'], element=atom['element']))
                # Add the residue to the atom
                sys.atoms.loc[a, 'res'] = new_residues[-1]

    # Add the residues to the sys.sol
    sys.sol.residues += new_residues
