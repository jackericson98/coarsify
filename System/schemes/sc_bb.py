from System.sys_funcs.calcs import calc_com, calc_dist
from System.sys_objs.ball import Ball


def coarsify_sc_bb_simple(sys, avg_dist=False, therm_cush=0.5, nuc_loc=None, am_loc=None):
    """
    Main coarsify function. Calculates radii and location for residues
    """
    # Check to see if the system has balls yet or not
    if sys.balls is None:
        sys.balls = []
    # Loop through the residues in the system
    for res in sys.residues:

        # If the residue is water we don't need to worry about backbone and side chain
        if res.name in sys.nucleics:
            bb_atoms, sc_atoms = [], []
            # Get the back bone atoms and the side chain atoms
            for atom in res.atoms:
                if atom['name'] in sys.nucleic_bbs:
                    bb_atoms.append(atom)
                elif atom['name'] in sys.nucleic_scs:
                    sc_atoms.append(atom)
                elif atom['name'] in sys.nucleic_ignores:
                    continue
                else:
                    sc_bb_add = input("{} atom type not found in nucleic list (Residue {} # {})! Choose one of the "
                                      "following - 1. Back Bone, 2. Side Chain, 3. Ignore \'{}\' 4. Skip\n>>> "
                                      .format(atom['name'], res.name, res.seq, atom['name']))
                    if sc_bb_add == '1':
                        sys.nucleic_bbs.append(atom['name'])
                        print("{} will be added to nucleic backbones from here on out".format(atom['name']))
                    elif sc_bb_add == '2':
                        sys.nucleic_scs.append(atom['name'])
                        print("{} will be added to nucleic side chains from here on out".format(atom['name']))
                    elif sc_bb_add == 3:
                        sys.nucleic_ignores.append(atom['name'])
                        print("{} will be ignored from here on out".format(atom['name']))
            # Calculate the center of mass for the atoms in a residue
            if nuc_loc is None:
                bb_loc = calc_com([_['loc'] for _ in bb_atoms])
                sc_loc = calc_com([_['loc'] for _ in sc_atoms])
            # Choose the scheme for coarse graining the residues
            if avg_dist:
                # Find the average distance from the center of mass for the backbone
                bb_rad = sum([calc_dist(bb_loc, a['loc']) + a['rad'] for a in bb_atoms]) / len(bb_atoms) + therm_cush
                # Find the average distance from the center of mass for the side chain
                sc_rad = sum([calc_dist(sc_loc, a['loc']) + a['rad'] for a in sc_atoms]) / len(sc_atoms) + therm_cush
            else:
                # Find the maximum of the summations of atom radii and the distance to residue com
                bb_rad = max([calc_dist(bb_loc, a['loc']) + a['rad'] for a in bb_atoms]) + therm_cush
                # Find the maximum of the summations of atom radii and the distance to residue com
                sc_rad = max([calc_dist(sc_loc, a['loc']) + a['rad'] for a in sc_atoms]) + therm_cush
            # Create the ball object
            sys.balls += [
                Ball(loc=bb_loc, rad=bb_rad, element=res.elem_col, residues=[res], atoms=res.atoms, name=res.name,
                     chain=res.chain, seq=res.seq, residue_subsection='bb'),
                Ball(loc=sc_loc, rad=sc_rad, element='pb', residues=[res], atoms=res.atoms, name=res.name,
                     chain=res.chain, seq=res.seq, residue_subsection='sc')]
        elif res.name in sys.aminos:
            bb_atoms, sc_atoms = [], []
            # Get the back bone atoms and the side chain atoms
            for atom in res.atoms:
                if atom['name'] in sys.amino_bbs:
                    bb_atoms.append(atom)
                elif atom['name'] in sys.amino_scs:
                    sc_atoms.append(atom)
                elif atom['name'] in sys.amino_ignores:
                    continue
                else:
                    sc_bb_add = input("{} atom type not found in amino list (Residue {} # {})! Please choose an option "
                                      "-  Add to 1. Back Bone, 2. Side Chain, 3. Ignore \'{}\', 4. Pass\n>>> "
                                      .format(atom['name'], res.name, res.seq, atom['name']))
                    if sc_bb_add == '1':
                        sys.amino_bbs.append(atom['name'])
                        print("{} will be added to amino backbones from here on out".format(atom['name']))
                    elif sc_bb_add == '2':
                        sys.amino_scs.append(atom['name'])
                        print("{} will be added to amino side chains from here on out".format(atom['name']))
                    elif sc_bb_add == '3':
                        sys.amino_ignores.append(atom['name'])
                        print("{} will be ignored from here on out".format(atom['name']))
                    else:
                        continue
            # Calculate the center of mass for the atoms in a residue
            try:
                bb_loc = calc_com([_['loc'] for _ in bb_atoms])
                sc_loc = calc_com([_['loc'] for _ in sc_atoms])
            except IndexError:
                continue
            # Choose the scheme for coarse graining the residues
            if avg_dist:
                # Find the average distance from the center of mass for the backbone
                bb_rad = sum([calc_dist(bb_loc, a['loc']) + a['rad'] for a in bb_atoms]) / len(bb_atoms) + therm_cush
                # Find the average distance from the center of mass for the side chain
                sc_rad = sum([calc_dist(sc_loc, a['loc']) + a['rad'] for a in sc_atoms]) / len(sc_atoms) + therm_cush
            else:
                # Find the maximum of the summations of atom radii and the distance to residue com
                bb_rad = max([calc_dist(bb_loc, a['loc']) + a['rad'] for a in bb_atoms]) + therm_cush
                # Find the maximum of the summations of atom radii and the distance to residue com
                sc_rad = max([calc_dist(sc_loc, a['loc']) + a['rad'] for a in sc_atoms]) + therm_cush
            # Create the ball object
            sys.balls += [
                Ball(loc=bb_loc, rad=bb_rad, element=res.elem_col, residues=[res], atoms=res.atoms, name=res.name,
                     chain=res.chain, seq=res.seq, residue_subsection='bb'),
                Ball(loc=sc_loc, rad=sc_rad, element='pb', residues=[res], atoms=res.atoms, name=res.name,
                     chain=res.chain, seq=res.seq, residue_subsection='sc')]
        else:
            # Calculate the center of mass for the atoms in a residue
            loc = calc_com([_['loc'] for _ in res.atoms])
            # Find the average distance from the center of mass
            if avg_dist:
                rad = sum([calc_dist(loc, a['loc']) + a['rad'] for a in res.atoms]) / len(res.atoms) + therm_cush
            else:
                # Find the maximum of the summations of atom radii and the distance to residue com
                rad = max([calc_dist(loc, a['loc']) + a['rad'] for a in res.atoms]) + therm_cush
            # Create the ball object
            sys.balls.append(Ball(loc=loc, rad=rad, element=res.elem_col, residues=[res], atoms=res.atoms,
                                  name=res.name, chain=res.chain, seq=res.seq))


def coarsify_sc_bb(sys, avg_dist=False, therm_cush=0.5, nuc_loc=None, am_loc=None):
    """
    Main coarsify function. Calculates radii and location for residues
    """
    # Check to see if the system has balls yet or not
    if sys.balls is None:
        sys.balls = []
    # Loop through the residues in the system
    for res in sys.residues:
        # If the residue is water we don't need to worry about backbone and side chain
        if res.name in sys.nucleics:
            bb_atoms, sc_atoms = [], []
            # Get the back bone atoms and the side chain atoms
            for atom in res.atoms:
                if atom['name'] in sys.nucleic_bbs:
                    bb_atoms.append(atom)
                elif atom['name'] in sys.nucleic_scs:
                    sc_atoms.append(atom)
                elif atom['name'] in sys.nucleic_ignores:
                    continue
                else:
                    sc_bb_add = input("{} atom type not found in nucleic list (Residue {} # {})! Choose one of the "
                                      "following - 1. Back Bone, 2. Side Chain, 3. Ignore \'{}\' 4. Skip\n>>> "
                                      .format(atom['name'], res.name, res.seq, atom['name']))
                    if sc_bb_add == '1':
                        sys.nucleic_bbs.append(atom['name'])
                        print("{} will be added to nucleic backbones from here on out".format(atom['name']))
                    elif sc_bb_add == '2':
                        sys.nucleic_scs.append(atom['name'])
                        print("{} will be added to nucleic side chains from here on out".format(atom['name']))
                    elif sc_bb_add == 3:
                        sys.nucleic_ignores.append(atom['name'])
                        print("{} will be ignored from here on out".format(atom['name']))
            # Calculate the center of mass for the atoms in a residue
            if nuc_loc is None:
                bb_loc = calc_com([_['loc'] for _ in bb_atoms])
                sc_loc = calc_com([_['loc'] for _ in sc_atoms])
            # Choose the scheme for coarse graining the residues
            if avg_dist:
                # Find the average distance from the center of mass for the backbone
                bb_rad = sum([calc_dist(bb_loc, a['loc']) + a['rad'] for a in bb_atoms]) / len(bb_atoms) + therm_cush
                # Find the average distance from the center of mass for the side chain
                sc_rad = sum([calc_dist(sc_loc, a['loc']) + a['rad'] for a in sc_atoms]) / len(sc_atoms) + therm_cush
            else:
                # Find the maximum of the summations of atom radii and the distance to residue com
                bb_rad = max([calc_dist(bb_loc, a['loc']) + a['rad'] for a in bb_atoms]) + therm_cush
                # Find the maximum of the summations of atom radii and the distance to residue com
                sc_rad = max([calc_dist(sc_loc, a['loc']) + a['rad'] for a in sc_atoms]) + therm_cush
            # Create the ball object
            sys.balls += [
                Ball(loc=bb_loc, rad=bb_rad, element=res.elem_col, residues=[res], atoms=res.atoms, name=res.name,
                     chain=res.chain, seq=res.seq, residue_subsection='bb'),
                Ball(loc=sc_loc, rad=sc_rad, element='pb', residues=[res], atoms=res.atoms, name=res.name,
                     chain=res.chain, seq=res.seq, residue_subsection='sc')]
        elif res.name in sys.aminos:
            bb_atoms, sc_atoms = [], []
            # Get the back bone atoms and the side chain atoms
            for atom in res.atoms:
                if atom['name'] in sys.amino_bbs:
                    bb_atoms.append(atom)
                elif atom['name'] in sys.amino_scs:
                    sc_atoms.append(atom)
                elif atom['name'] in sys.amino_ignores:
                    continue
                else:
                    sc_bb_add = input("{} atom type not found in amino list (Residue {} # {})! Please choose an option "
                                      "-  Add to 1. Back Bone, 2. Side Chain, 3. Ignore \'{}\', 4. Pass\n>>> "
                                      .format(atom['name'], res.name, res.seq, atom['name']))
                    if sc_bb_add == '1':
                        sys.amino_bbs.append(atom['name'])
                        print("{} will be added to amino backbones from here on out".format(atom['name']))
                    elif sc_bb_add == '2':
                        sys.amino_scs.append(atom['name'])
                        print("{} will be added to amino side chains from here on out".format(atom['name']))
                    elif sc_bb_add == '3':
                        sys.amino_ignores.append(atom['name'])
                        print("{} will be ignored from here on out".format(atom['name']))
                    else:
                        continue
            # Calculate the center of mass for the atoms in a residue
            try:
                bb_loc = calc_com([_['loc'] for _ in bb_atoms])
                sc_loc = calc_com([_['loc'] for _ in sc_atoms])
            except IndexError:
                continue
            # Choose the scheme for coarse graining the residues
            if avg_dist:
                # Find the average distance from the center of mass for the backbone
                bb_rad = sum([calc_dist(bb_loc, a['loc']) + a['rad'] for a in bb_atoms]) / len(bb_atoms) + therm_cush
                # Find the average distance from the center of mass for the side chain
                sc_rad = sum([calc_dist(sc_loc, a['loc']) + a['rad'] for a in sc_atoms]) / len(sc_atoms) + therm_cush
            else:
                # Find the maximum of the summations of atom radii and the distance to residue com
                bb_rad = max([calc_dist(bb_loc, a['loc']) + a['rad'] for a in bb_atoms]) + therm_cush
                # Find the maximum of the summations of atom radii and the distance to residue com
                sc_rad = max([calc_dist(sc_loc, a['loc']) + a['rad'] for a in sc_atoms]) + therm_cush
            # Create the ball object
            sys.balls += [
                Ball(loc=bb_loc, rad=bb_rad, element=res.elem_col, residues=[res], atoms=res.atoms, name=res.name,
                     chain=res.chain, seq=res.seq, residue_subsection='bb'),
                Ball(loc=sc_loc, rad=sc_rad, element='pb', residues=[res], atoms=res.atoms, name=res.name,
                     chain=res.chain, seq=res.seq, residue_subsection='sc')]
        else:
            # Calculate the center of mass for the atoms in a residue
            loc = calc_com([_['loc'] for _ in res.atoms])
            # Find the average distance from the center of mass
            if avg_dist:
                rad = sum([calc_dist(loc, a['loc']) + a['rad'] for a in res.atoms]) / len(res.atoms) + therm_cush
            else:
                # Find the maximum of the summations of atom radii and the distance to residue com
                rad = max([calc_dist(loc, a['loc']) + a['rad'] for a in res.atoms]) + therm_cush
            # Create the ball object
            sys.balls.append(Ball(loc=loc, rad=rad, element=res.elem_col, residues=[res], atoms=res.atoms,
                                  name=res.name, chain=res.chain, seq=res.seq))
