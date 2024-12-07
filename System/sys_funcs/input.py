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
            oxy_res.append(Residue(sys=residue.sys, atoms=[a], name=atom['residue'],
                                   sequence=atom['res_seq'], chain=atom['chn']))
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
                    hydrogens.remove(h)
                if len(res.atoms) == 3:
                    break
            # print([(sys.atoms['name'][_], sys.atoms['res_seq'][_], sys.atoms['loc'][_][0]) for _ in res.atoms])
        good_resids.append(res)

    # Las sort the hydrogens
    if len(hydrogens) == 1:
        h = hydrogens[0]
        good_resids.append(Residue(sys=residue.sys, atoms=hydrogens, name=sys.atoms['name'][h],
                                   sequence=sys.atoms['res_seq'][h], chain=sys.atoms['chn'][h]))
    elif len(hydrogens) == 2:
        h1, h2 = sys.atoms.iloc[hydrogens[0]], sys.atoms.iloc[hydrogens[1]]
        if calc_dist(h1['loc'], h2['loc']) < 3 and h1['name'] != h2['name']:
            good_resids.append(Residue(sys=residue.sys, atoms=hydrogens, name=h1['residue'],
                                       sequence=h1['res_seq'], chain=h1['chn']))
        else:
            for h in hydrogens:
                good_resids.append(Residue(sys=residue.sys, atoms=[h], name=sys.atoms['name'][h],
                                           sequence=sys.atoms['res_seq'][h], chain=sys.atoms['chn'][h]))
    elif len(hydrogens) == 4:
        # Get the first hydrogen
        h1 = sys.atoms.iloc[hydrogens[0]]
        # Find it's pair
        for h in hydrogens[1:]:
            if sys.atoms['name'][h][:-1] == h1['name'][:-1] and calc_dist(sys.atoms['loc'][h], h1['loc']) < 3:
                hydrogens = [_ for _ in hydrogens[1:] if _ != h]
                good_resids.append(Residue(sys=residue.sys, atoms=[h1['num'], h], name=h1['residue'],
                                           sequence=h1['res_seq'], chain=h1['chn']))
        # Get the first hydrogen
        h1 = sys.atoms.iloc[hydrogens[0]]
        good_resids.append(Residue(sys=residue.sys, atoms=hydrogens, name=h1['residue'],
                                   sequence=h1['res_seq'], chain=h1['chn']))
    # print([(sys.atoms['name'][_], sys.atoms['res_seq'][_], sys.atoms['loc'][_][0]) for _ in hydrogens])
    for res in good_resids:

        if len(res.atoms) > 3:

            print([(sys.atoms['name'][_], sys.atoms['loc'][_], sys.atoms['element'][_], res.seq) for _ in res.atoms])
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
                         name=name.strip(), seg_id=line[72:76], index=atom_count, mass=mass)
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

    # Set the atoms and the data
    sys.atoms, sys.settings = DataFrame(atoms), data
    # Adjust the SOL residues
    adjusted_residues = []
    for res in sys.sol.residues:
        if len(res.atoms) > 3:
            try:
                adjusted_residues += fix_sol(sys, res)
            except TypeError:
                print(res.atoms)
        else:
            adjusted_residues.append(res)

    sys.sol.residues = adjusted_residues
