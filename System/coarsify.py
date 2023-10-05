from System.calcs import calc_dist, calc_com
from System.objects import Ball
from System.cg_designations import proteins, nucleobases, ions, solvents


def coarsify_basic(sys, average_dist=False, encapsulate=False, therm_cush=0.5):
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
        # Choose the scheme for coarse graining the residues
        if average_dist:
            # Find the average distance from the center of mass
            rad = sum([calc_dist(loc, a['loc']) + a['rad'] for a in res.atoms]) / len(res.atoms) + therm_cush
        elif encapsulate:
            # Find the maximum of the summations of atom radii and the distance to residue com
            rad = max([calc_dist(loc, a['loc']) + a['rad'] for a in res.atoms]) + therm_cush
        # Create the ball object
        sys.balls.append(Ball(loc=loc, rad=rad, element=res.elem_col, residues=[res], atoms=res.atoms, name=res.name,
                              chain=res.chain, seq=res.seq))


def coarsify_CG(sys, therm_cush=0.0):
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
                              name=atom['residue'], chain=atom['chn'], seq=atom['res'].seq))
