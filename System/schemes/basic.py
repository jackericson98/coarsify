from System.sys_funcs.calcs import calc_com, calc_dist
from System.sys_objs.ball import Ball

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

