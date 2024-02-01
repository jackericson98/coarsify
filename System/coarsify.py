from System.calcs import calc_dist, calc_com
from System.objects import Ball
from System.cg_designations import proteins, nucleobases, ions, solvents


def coarsify_encapsulate(sys, therm_cush=0.5):
    """
    Main coarsify function. Calculates radii and location for residues
    """
    # Check to see if the system has balls yet or not
    if sys.balls is None:
        sys.balls = []
    # Loop through the residues in the system
    for res in sys.residues:
        # Calculate the center of mass for the atoms in a residue
        loc = calc_com([_['loc'] for _ in res.atoms])
        # Find the maximum of the summations of atom radii and the distance to residue com
        rad = max([calc_dist(loc, a['loc']) + a['rad'] for a in res.atoms]) + therm_cush
        # Create the ball object
        sys.balls.append(Ball(loc=loc, rad=rad, element=res.elem_col, residues=[res], atoms=res.atoms, name=res.name,
                              chain=res.chain, seq=res.seq))


def coarsify_avg_dist(sys, therm_cush=0.5):
    """
    Main coarsify function. Calculates radii and location for residues
    """
    # Check to see if the system has balls yet or not
    if sys.balls is None:
        sys.balls = []
    # Loop through the residues in the system
    for res in sys.residues:
        # Calculate the center of mass for the atoms in a residue
        loc = calc_com([_['loc'] for _ in res.atoms])
        # Find the average distance from the center of mass
        rad = sum([calc_dist(loc, a['loc']) + a['rad'] for a in res.atoms]) / len(res.atoms) + therm_cush
        # Create the ball object
        sys.balls.append(Ball(loc=loc, rad=rad, element=res.elem_col, residues=[res], atoms=res.atoms, name=res.name,
                              chain=res.chain, seq=res.seq))


def coarsify_sc_bb(sys, avg_dist=True, therm_cush=0.5):
    """
    Main coarsify function. Calculates radii and location for residues
    """
    # Check to see if the system has balls yet or not
    if sys.balls is None:
        sys.balls = []
    # Loop through the residues in the system
    for res in sys.residues:
        # If the residue is water we don't need to worry about backbone and side chain
        if res.chain.name in {'Z', 'X'}:
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
            continue
        # Get the back bone atoms and the side chain atoms
        bb_atoms = [_ for _ in res.atoms if _['name'] in sys.bb_names]
        sc_atoms = [_ for _ in res.atoms if _['name'] not in sys.bb_names]
        # Calculate the center of mass for the atoms in a residue
        try:
            bb_loc = calc_com([_['loc'] for _ in bb_atoms])
            sc_loc = calc_com([_['loc'] for _ in sc_atoms])
        except IndexError:
            print(res.name, res.chain.name, len(res.atoms), [_['name'] for _ in res.atoms])
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
        sys.balls += [Ball(loc=bb_loc, rad=bb_rad, element=res.elem_col, residues=[res], atoms=res.atoms, name=res.name,
                           chain=res.chain, seq=res.seq, residue_subsection='bb'),
                      Ball(loc=sc_loc, rad=sc_rad, element='pb', residues=[res], atoms=res.atoms, name=res.name,
                           chain=res.chain, seq=res.seq, residue_subsection='sc')]


def coarsify_c_alpha(sys, therm_cush=0.0):
    """
    C Alpha coarse graining scheme
    """
    # Check to see if the system has balls yet or not
    if sys.balls is None:
        sys.balls = []
    return


def coarsify_martini(sys, therm_cush=0.0):
    """
    CG Martini generated pdb file coarsify function
    """
    # Check to see if the system has balls yet or not
    if sys.balls is None:
        sys.balls = []
    # Go through the atoms in the system
    for i, atom in sys.atoms.iterrows():
        # Proteins
        if atom['residue'] in proteins:
            # Get the radius
            rad = proteins[atom['residue']][atom['name']]['size'] + therm_cush
        # Nucleic acid bases
        elif atom['residue'] in nucleobases:
            # Get the radius
            rad = nucleobases[atom['residue']][atom['name']]['size'] + therm_cush
        # Ions
        elif atom['residue'] in ions:
            # Get the radius
            rad = ions[atom['residue']][atom['name']]['size'] + therm_cush
        # Solvents
        elif atom['residue'] in solvents:
            # Get the radius
            rad = solvents[atom['residue']][atom['name']]['size'] + therm_cush
        else:
            rad = atom['rad']
        # Create the ball
        sys.balls.append(Ball(loc=atom['loc'], rad=rad, element=atom['res'].elem_col, residues=atom['res'],
                              name=atom['residue'], chain=atom['chn'], seq=atom['res'].seq, index=atom['num']))


def coarsify_primo(sys, therm_cush=0.0):
    """
    primo coarse graining model
    :param sys:
    :param therm_cush:
    :return:
    """
    sys.balls = []
    return