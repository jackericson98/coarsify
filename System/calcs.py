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


def calc_com(points):
    """
    Takes in a set of points and returns the coordinates of the center of mass
    :param points: lists of locations in n-dimensions
    :return: Center of mass of the inputs
    """

    # Set the running sum for the x, y, z values to 0
    tots = [0 for _ in range(len(points[0]))]
    for point in points:
        for i in range(len(points[0])):
            tots[i] += point[i]

    # Return the center of mass of inputs
    return [tots[i]/len(points) for i in range(len(points[0]))]


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
    if atom['rad'] is None and atom['element'].lower() in radii:
        atom['rad'] = radii[atom['element'].lower()]
    # If indicated we return the symbol of atom that the radius indicates
    if atom['rad'] is None or atom['rad'] == 0:
        # Check to see if the radius is in the system
        if atom['rad'] in {radii[_] for _ in radii[1]}:
            atom['element'] = radii[atom['rad']]
        else:
            # Get the closest atom to it
            min_diff = np.inf
            # Go through the radii in the system looking for the smallest difference
            for radius in radii:
                if radii[radius] - atom['rad'] < min_diff:
                    atom['element'] = radii[radius]
    return atom['rad']
