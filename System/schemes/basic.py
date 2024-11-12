from System.sys_funcs.calcs import calc_com, calc_dist, calc_tetra_vol, get_time
import time
from System.sys_objs.ball import Ball
import numpy as np
from itertools import combinations
from System.sys_funcs.verts import calc_vert
from mpl_visualize import plot_verts, plot_balls
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull


def minimum_enclosing_sphere_for_two(sphere1, sphere2):
    """
    Calculates the minimum enclosing sphere for two given spheres.
    Each sphere is defined by a tuple of (center, radius), where `center` is a numpy array.
    """
    center1, radius1 = sphere1
    center2, radius2 = sphere2

    # Calculate the vector from the center of sphere1 to the center of sphere2
    center_vector = center2 - center1
    distance_between_centers = np.linalg.norm(center_vector)

    if distance_between_centers == 0:
        # Both spheres are at the same center
        if radius1 == radius2:
            return center1, radius1
        elif radius1 > radius2:
            return center1, radius1
        else:
            return center2, radius2

    # Normalize the center vector
    direction = center_vector / distance_between_centers

    # Calculate the new center point
    # We use proportion to determine the correct placement between the two centers
    if radius1 > radius2:
        center = center1 + direction * (distance_between_centers + radius2 - radius1) / 2
    else:
        center = center1 + direction * (distance_between_centers - radius2 + radius1) / 2

    # The radius of the enclosing sphere
    enclosing_radius = np.linalg.norm(center - center1) + radius1

    return center, enclosing_radius


def find_circumcenter(a, b, c):
    """Find the circumcenter of the triangle formed by points a, b, c."""
    ac = c - a
    ab = b - a
    ab_cross_ac = np.cross(ab, ac)
    # Squared lengths of the edges
    ab_sq = np.dot(ab, ab)
    ac_sq = np.dot(ac, ac)
    cross_sq = np.dot(ab_cross_ac, ab_cross_ac)

    # Compute the circumcenter
    if cross_sq == 0:
        return None  # Degenerate case: points are collinear or coincident

    s = 0.5 / cross_sq * np.array([
        ac_sq * np.dot(ab, ab_cross_ac),
        ab_sq * np.dot(ac, -ab_cross_ac)
    ])
    circumcenter = a + s[0] * ab + s[1] * ac
    return circumcenter


def enclosing_sphere_radius(center, spheres):
    """Calculate the minimum radius of the sphere that encloses the given spheres with given center."""
    max_radius = 0
    for (sphere_center, radius) in spheres:
        distance = np.linalg.norm(center - sphere_center) + radius
        if distance > max_radius:
            max_radius = distance
    return max_radius


def minimum_enclosing_sphere_3(spheres):
    """Find the minimum enclosing sphere for three spheres."""
    if len(spheres) != 3:
        raise ValueError("Exactly three spheres are required.")

    centers = [s[0] for s in spheres]
    circumcenter = find_circumcenter(np.array(centers[0]), np.array(centers[1]), np.array(centers[2]))

    if circumcenter is None:
        # Handle collinear or coincident points, a simpler enclosing sphere is needed
        return None

    radius = enclosing_sphere_radius(circumcenter, spheres)
    return circumcenter, radius


def find_enclosing_sphere1(spheres):
    max_dist, my_spheres = 0, None
    for i in range(len(spheres)):
        loc_i, rad_i = spheres[i]
        for j in range(i + 1, len(spheres)):
            loc_j, rad_j = spheres[j]
            dist = calc_dist(loc_i, loc_j) + rad_i + rad_i
            if dist > max_dist:
                my_spheres = [spheres[i], spheres[j]]
                max_dist = dist
    s1, s2 = my_spheres
    my_dir = s2[0] - s1[0]
    dhat = my_dir / np.linalg.norm(my_dir)
    loc = s1[0] + (0.5 * max_dist - s1[1]) * dhat
    if encapsulates_all_spheres(loc, max_dist / 2, spheres):
        return loc, max_dist / 2


def find_enclosing_sphere(spheres):
    # spheres is a list of tuples (center, radius), where center is a numpy array
    num_spheres = len(spheres)

    # Calculate the farthest points in the set of spheres to define the initial bounds
    max_distance = 0
    for i in range(num_spheres):
        center_i, radius_i = spheres[i]
        for j in range(i + 1, num_spheres):
            center_j, radius_j = spheres[j]
            distance = np.linalg.norm(center_i - center_j) + radius_i + radius_j
            if distance > max_distance:
                max_distance = distance
                # Initial guess for the enclosing sphere center and radius
                enclosing_center = (center_i + center_j) / 2
                enclosing_radius = distance / 2

    # Iteratively adjust the center of the sphere to minimize the radius
    for _ in range(100):  # Perform 100 iterations for convergence
        max_dist = 0
        for center, radius in spheres:
            dist = np.linalg.norm(enclosing_center - center) + radius
            if dist > max_dist:
                max_dist = dist
                farthest_sphere = (center, radius)

        # Update the sphere to include the farthest sphere
        if max_dist > enclosing_radius:
            enclosing_radius = max_dist
            direction = farthest_sphere[0] - enclosing_center
            direction /= np.linalg.norm(direction)
            enclosing_center += direction * (max_dist - enclosing_radius)

    return enclosing_center, enclosing_radius


def calculate_center_with_radii(spheres):
    """Calculate the center considering both the position and radii of the spheres."""
    total_weight = sum(radius for _, radius in spheres)
    weighted_centers = sum(center * radius for center, radius in spheres)
    center_with_radii = weighted_centers / total_weight
    return center_with_radii


def encapsulates_all_spheres(loc, rad, spheres):
    """Check if a candidate sphere encapsulates all the spheres."""
    for center, radius in spheres:
        if round(calc_dist(center, loc) + radius, 4) > round(abs(rad), 4):
            return False
    return True


def minimum_enclosing_sphere1(spheres):
    """
    Calculate the minimum enclosing sphere for a set of spheres using their convex hull.

    Parameters:
    - spheres: list of tuples (center, radius), where center is a numpy array.

    Returns:
    - Tuple (center, radius) of the minimum enclosing sphere.
    """
    # Extract centers
    centers = np.array([center for center, _ in spheres])

    # Compute the convex hull of the centers
    hull = ConvexHull(centers)
    hull_points = centers[hull.vertices]

    # Use the Ritter's bounding sphere algorithm on hull points (a simple and effective heuristic)
    # Start with an initial sphere encompassing two hull points
    center = (hull_points[0] + hull_points[1]) / 2
    radius = np.linalg.norm(hull_points[0] - center)

    # Adjust sphere to include all hull points
    for point in hull_points:
        d = np.linalg.norm(point - center)
        if d > radius:
            # Expand the sphere to include the new point
            radius = (radius + d) / 2
            direction = (point - center) / d
            center += (d - radius) * direction

    # Adjust radius to include the original spheres' radii
    for center_p, radius_p in spheres:
        distance = np.linalg.norm(center_p - center)
        required_radius = distance + radius_p
        if required_radius > radius:
            radius = required_radius

    return center, radius


def minimum_enclosing_sphere(spheres, plotting=False):
    """
    Find the minimum enclosing sphere by considering combinations of spheres sorted by their distance from the center of
    mass.
    """

    # If there is only one sphere return the sphere
    if len(spheres) == 1:
        return spheres[0][0], spheres[0][1]

    # If there are multiple spheres, the simplest case is the enclosing sphere is defined by the furthest spheres
    my_loc_rad = find_enclosing_sphere1(spheres)
    if my_loc_rad is not None:
        return my_loc_rad

    # In the case of three spheres, the next best case is to do the three case with circumcenter and enclosing radius
    if len(spheres) == 3:
        return minimum_enclosing_sphere_3(spheres)

    # Set up the tetra spheres list
    tetra_spheres = []
    # Loop through the sphere combinations finding the tetrahedron volumes to see which provides the largest
    for four_spheres in combinations(spheres, 4):
        # Calculate the tetrahedron volumes
        vol = calc_tetra_vol(*[_[0] for _ in four_spheres])
        # Add the spheres and the volumes to new lists
        tetra_spheres.append((four_spheres, vol))

    # Sort the list by tetrahedron volume
    sorted_spheres = sorted(tetra_spheres, key=lambda x: x[1], reverse=True)
    sorted_spheres = [_[0] for _ in sorted_spheres]

    # Set up the minimum radius and minimum location for the calculation of vertices
    min_rad, min_loc = np.inf, None
    # Check combinations starting with those spheres furthest from the center of mass
    for four_spheres in sorted_spheres[:int(len(sorted_spheres) / 2)]:
        loc_rad = calc_vert(four_spheres)

        if loc_rad is None:
            continue
        loc, rad = loc_rad
        rad = abs(rad)
        if encapsulates_all_spheres(loc, rad, spheres):
            if rad is not None and rad < min_rad:
                min_rad = rad
                min_loc = loc

    my_loc, my_rad = find_enclosing_sphere(spheres)
    if (min_loc is None or my_rad < min_rad) and encapsulates_all_spheres(my_loc, my_rad, spheres):
        min_loc, min_rad = my_loc, my_rad

    my_loc, my_rad = minimum_enclosing_sphere1(spheres)
    if min_loc is None or my_rad < min_rad and encapsulates_all_spheres(my_loc, my_rad, spheres):
        min_loc, min_rad = my_loc, my_rad

    return min_loc, min_rad


def calculate_average_sphere(spheres):
    # spheres is a list of tuples (center, radius), where center is a numpy array
    centers = np.array([s[0] for s in spheres])

    # Calculate the centroid of the centers
    centroid = np.mean(centers, axis=0)

    # Calculate the average distance from the centroid to the surface of each sphere
    distances = np.array([np.linalg.norm(centroid - center) + radius for center, radius in spheres])
    average_radius = np.mean(distances)

    return centroid, average_radius


def calculate_weighted_average_sphere(spheres):
    """
    Calculates the weighted average sphere given a list of spheres with weights.
    Each element in `spheres` should be a tuple ((center, radius), weight),
    where `center` is a numpy array.

    Parameters:
    - spheres: list of tuples ((center, radius), weight)

    Returns:
    - Tuple (centroid, average_radius): Center and radius of the weighted average sphere.
    """
    if not spheres:
        return None, None

    # Extract centers, radii, and weights
    centers = np.array([s[0][0] for s in spheres])
    radii = np.array([s[0][1] for s in spheres])
    weights = np.array([s[1] for s in spheres])

    # Calculate the weighted centroid of the centers
    total_weight = np.sum(weights)
    weighted_centers = centers * weights[:, np.newaxis]  # Broadcasting weights
    centroid = np.sum(weighted_centers, axis=0) / total_weight

    # Calculate the weighted average distance from the centroid to the surface of each sphere
    distances = np.array([np.linalg.norm(centroid - center) + radius for center, radius in zip(centers, radii)])
    weighted_distances = distances * weights
    average_radius = np.sum(weighted_distances) / total_weight

    return centroid, average_radius


def make_ball(atoms, scheme='Average Distance', mass_weighted=True, include_h=True, therm_cush=0.0):
    # If the length of the atoms is 1 just return that location and radius
    if len(atoms) == 1:
        return atoms[0]['loc'], atoms[0]['rad']
    # If we are including h
    if not include_h:
        atoms = [_ for _ in atoms if _['element'].lower() != 'h']
    # Choose the scheme for coarse graining the residues
    if scheme == 'Average Distance' and mass_weighted:
        loc, rad = calculate_weighted_average_sphere([((a['loc'], a['rad']), a['mass']) for a in atoms])
    elif scheme == 'Average Distance':
        loc, rad = calculate_average_sphere([(a['loc'], a['rad']) for a in atoms])
    else:
        loc, rad = minimum_enclosing_sphere([(_['loc'], _['rad']) for _ in atoms])

    return loc, rad


def coarsify(sys):
    """
    Main coarsify function. Calculates radii and location for residues
    """
    # Check to see if the system has balls yet or not
    if sys.balls is None:
        sys.balls = []
    # Set the start time
    start_time = time.perf_counter()
    # Loop through the residues in the system
    for i, res in enumerate(sys.residues):
        # Get the percentage and print it
        percentage = min((i / (len(sys.residues) - len(sys.sol.residues))) * 100, 100)
        my_time = time.perf_counter() - start_time
        h, m, s = get_time(my_time)
        print("\rRun Time = {}:{:02d}:{:2.2f} - Coarsifying Residue {:>5}, {:<3} {:<5} - {:.2f} %"
              .format(int(h), int(m), round(s, 2), i + 1, res.name, res.seq, percentage), end="")
        # If the residue is water we don't need to worry about backbone and side chain
        if res.name in sys.nucleics and sys.sc_bb:
            sugs, phos, nbase, hs, mass = [], [], [], [], 0
            # Get the back bone atoms and the side chain atoms
            for atom in res.atoms:
                if atom['element'].lower() == 'h':
                    hs.append(atom)
                elif atom['name'] in sys.nucleic_pphte:
                    phos.append(atom)
                elif atom['name'] in sys.nucleic_nbase:
                    sugs.append(atom)
                elif atom['name'] in sys.nucleic_sugrs:
                    nbase.append(atom)
                elif atom['name'] in sys.nucleic_ignores:
                    continue
                else:
                    sc_bb_add = input("{} atom type not found in nucleic list (Residue {} # {})! Choose one of the "
                                      "following - 1. Phosphate, 2. Ribose, 3. Base, 4. Ignore \'{}\' 5. Skip\n>>> "
                                      .format(atom['name'], res.name, res.seq, atom['name']))
                    if sc_bb_add == '1':
                        sys.nucleic_pphte.append(atom['name'])
                        phos.append(atom)
                        print("{} will be added to nucleic phosphate grups from here on out".format(atom['name']))
                    elif sc_bb_add == '2':
                        sys.nucleic_nbase.append(atom['name'])
                        sugs.append(atom)
                        print("{} will be added to nucleic ribose groups from here on out".format(atom['name']))
                    elif sc_bb_add == 3:
                        sys.nucleic_sugr.append(atom['name'])
                        nbase.append(atom)
                        print("{} will be added to nucleic bases from here on out".format(atom['name']))
                    elif sc_bb_add == 4:
                        sys.nucleic_ignores.append(atom['name'])
                        print("{} will be ignored from here on out".format(atom['name']))
                    else:
                        continue

            # Sort the unknown hydrogen into their different groups
            for atom in hs:
                atom_dists = [(calc_dist(atom['loc'], _['loc']), _) for _ in res.atoms if _['element'].lower() != 'h']
                near_atoms = sorted(atom_dists, key=lambda x: x[0])
                close_atom = near_atoms[0][1]
                if close_atom['name'] in sys.nucleic_pphte:
                    phos.append(atom)
                elif close_atom['name'] in sys.nucleic_nbase:
                    sugs.append(atom)
                elif close_atom['name'] in sys.nucleic_sugrs:
                    nbase.append(atom)
                elif close_atom['name'] in sys.nucleic_ignores:
                    continue

            # Create the ball objects
            if len(phos) > 0:
                mass = sum([_['mass'] for _ in phos])
                ph_loc, ph_rad = make_ball(phos, sys.scheme, sys.mass_weighted, sys.include_h, sys.therm_cush)
                sys.balls.append(Ball(loc=ph_loc, rad=ph_rad, element='bi', residues=[res], atoms=phos,
                                      name=res.name, chain=res.chain, seq=res.seq, residue_subsection='phosphate',
                                      mass=mass))
            if len(sugs) > 0:
                mass = sum([_['mass'] for _ in sugs])
                sug_loc, sug_rad = make_ball(sugs, sys.scheme, sys.mass_weighted, sys.include_h, sys.therm_cush)
                sys.balls.append(Ball(loc=sug_loc, rad=sug_rad, element='pb', residues=[res], atoms=sugs, name=res.name,
                     chain=res.chain, seq=res.seq, residue_subsection='sugar', mass=mass))

            if len(nbase) > 0:
                mass = sum([_['mass'] for _ in nbase])
                nbas_loc, nbas_rad = make_ball(nbase, sys.scheme, sys.mass_weighted, sys.include_h, sys.therm_cush)
                sys.balls.append(Ball(loc=nbas_loc, rad=nbas_rad, element='po', residues=[res], atoms=nbase, name=res.name,
                                 chain=res.chain, seq=res.seq, residue_subsection='nbase', mass=mass))

        elif res.name in sys.aminos and sys.sc_bb:
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
            # Get the balls for the aminos
            if len(bb_atoms) > 0:
                mass = sum([_['mass'] for _ in bb_atoms])
                bb_loc, bb_rad = make_ball(bb_atoms, sys.scheme, sys.mass_weighted, sys.include_h, sys.therm_cush)
                sys.balls.append(Ball(loc=bb_loc, rad=bb_rad, element=res.elem_col, residues=[res], atoms=bb_atoms,
                                      name=res.name, chain=res.chain, seq=res.seq, residue_subsection='bb', mass=mass))
            if len(sc_atoms) > 0:
                mass = sum([_['mass'] for _ in sc_atoms])
                sc_loc, sc_rad = make_ball(sc_atoms, sys.scheme, sys.mass_weighted, sys.include_h, sys.therm_cush)
                sys.balls.append(Ball(loc=sc_loc, rad=sc_rad, element='pb', residues=[res], atoms=sc_atoms,
                                      name=res.name, chain=res.chain, seq=res.seq, residue_subsection='sc', mass=mass))
        else:
            # Get the loc and rad
            loc, rad = make_ball(res.atoms, sys.scheme, sys.mass_weighted, sys.include_h, sys.therm_cush)
            # Calculate the mass of the ball
            mass = sum([_['mass'] for _ in res.atoms])
            # Create the ball object
            sys.balls.append(Ball(loc=loc, rad=rad, element=res.elem_col, residues=[res], atoms=res.atoms,
                                  name=res.name, chain=res.chain, seq=res.seq, mass=mass))
