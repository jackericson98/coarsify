import numpy as np


def calc_dist(l0, l1):
    """
    Calculate distance function used to simplify code
    :param l0: Point 0 list, array, n-dimensional must match point 1
    :param l1: Point 1 list, array, n-dimensional must match point 0
    :return: float distance between the two points
    """
    # Pythagorean theorem
    return np.sqrt(sum(np.square(l0 - l1)))


def calc_com(points, masses=None):
    """
    Takes in a set of points and returns the weighted center of mass. If no mass list provided all weighted at 1
    :param points: lists of locations in n-dimensions
    :param masses: list of masses corresponding to the points
    :return: Center of mass of the inputs
    """
    # Get the total mass if masses are included
    if masses is None:
        masses = [1 for _ in range(len(points))]
    # Set the running sum for the x, y, z values to 0
    tots = [0 for _ in range(len(points[0]))]
    for j in range(len(points)):
        for i in range(len(points[0])):
            tots[i] += points[j][i] * masses[j]

    # Return the center of mass of inputs
    return [tots[i]/sum(masses) for i in range(len(points[0]))]


# Calculate tetrahedron volume function.
def calc_tetra_vol(p0, p1, p2, p3):
    """
    Calculates the volume of a tetrahedron defined by its vertices
    :param p0: Point 0
    :param p1: Point 1
    :param p2: Point 2
    :param p3: Point 3
    :return: Volume of the tetrahedron made by the points
    """
    # Choose a base point (p0) and find the vectors between it and other points
    r01 = p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]
    r02 = p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]
    r03 = np.array([p3[0] - p0[0], p3[1] - p0[1], p3[2] - p0[2]])

    # Formula for tetrahedron volume: 1/6 * r03 dot (r01 cross r02)
    return (1/6)*abs(np.dot(r03, np.cross(r01, r02)))


def get_time(seconds):
    """
    Turns seconds into hours, minutes and seconds
    :param seconds: Number of seconds in the counter
    :return: hours, minutes, seconds
    """
    # Divide up the values
    hours = seconds // 3600
    minutes = (seconds - (hours * 3600)) // 60
    seconds = seconds - hours * 3600 - minutes * 60
    # Return the values
    return hours, minutes, seconds


def get_radius(atom):
    """
    Finds the radius of the atom from the symbol or vice versa
    :return: The radius of the atom from the symbol or vice versa
    """
    radii, special_radii = atom['sys'].radii, atom['sys'].special_radii
    # Get the radius and the element from the name of the atom
    if atom['res'] is not None and atom['res'].name in special_radii:
        # Check if no atom name exists or its empty
        if atom['name'] is not None and atom['name'] != '':
            for i in range(len(atom['name'])):
                name = atom['name'][:-i]
                # Check the residue name
                if name in special_radii[atom['res'].name]:
                    atom['rad'] = special_radii[atom['res'].name][name]
    # If we have the type and just want the radius, keep scanning until we find the radius
    if atom['rad'] is None and atom['element'].upper() in radii:
        atom['rad'] = radii[atom['element'].upper()]
    # If indicated we return the symbol of atom that the radius indicates
    if atom['rad'] is None or atom['rad'] == 0:
        # Check to see if the radius is in the system
        if atom['element'].upper() in radii:
            atom['rad'] = radii[atom['element'].upper()]
        else:
            # Get the closest atom to it
            min_diff = np.inf
            # Go through the radii in the system looking for the smallest difference
            for radius in radii:
                if radii[radius] - atom['rad'] < min_diff:
                    atom['element'] = radii[radius]
    return atom['rad']


def pdb_line(atom="ATOM", ser_num=0, name="", alt_loc=" ", res_name="", chain="A", res_seq=0, cfir="", x=0, y=0, z=0,
             occ=1, tfact=0, seg_id="", elem="", charge=""):
    """
    Takes in values for a line in a pdb file and places them in the correct locations
    :return: String for each line
    """
    # Write the line for the file
    return "{:<6}{:>5} {:<4}{:1}{:>3} {:^1}{:>4}{:1}   {:>8.3f}{:>8.3f}{:>8.3f}{:>6.2f}{:>6.2f}      {:<4}{:>2}{}\n"\
        .format(atom, ser_num, name, alt_loc, res_name, chain, res_seq, cfir, x, y, z, occ, tfact, seg_id, elem, charge)
