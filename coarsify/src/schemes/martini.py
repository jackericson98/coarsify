from coarsify.src.tools.cg_designations import *
from System.sys_objs.ball import Ball


def coarsify_martini(sys):
    """
    CG Martini generated pdb file coarsify function
    """
    therm_cush = sys.therm_cush
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
        sys.balls.append(Ball(loc=atom['loc'], rad=rad, element=atom['res'].element, residues=[atom['res']],
                              name=atom['residue'], chain=atom['chn'], seq=atom['res'].seq, index=atom['num'], mass=1))

