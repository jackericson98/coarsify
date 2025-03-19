from System.sys_funcs.calcs import pdb_line
import os


def set_dir(sys, dir_name=None):
    if not os.path.exists("./Data/user_data"):
        os.mkdir("./Data/user_data")
    # If no outer directory was specified use the directory outside the current one
    if dir_name is None:
        if sys.vpy_dir is not None:
            dir_name = sys.vpy_dir + "/Data/user_data/" + sys.name
        else:
            dir_name = os.getcwd() + "/Data/user_data/" + sys.name
    # Catch for existing directories. Keep trying out directories until one doesn't exist
    i = 0
    while True:
        # Try creating the directory with the system name + the current i_string
        try:
            # Create a string variable for the incrementing variable
            i_str = '_' + str(i)
            # If no file with the system name exists change the string to empty
            if i == 0:
                i_str = ""
            # Try to create the directory
            os.mkdir(dir_name + i_str)
            break
        # If the file exists increment the counter and try creating the directory again
        except FileExistsError:
            i += 1
    # Set the output directory for the system
    sys.dir = dir_name + i_str


def write_pdb(sys):
    # Catch empty atoms cases
    if sys.residues is None or len(sys.residues) == 0:
        return
    # Make note of the starting directory
    start_dir = os.getcwd()
    # Change to the specified directory
    os.chdir(sys.dir)

    # Open the file for writing
    with open(sys.name + ".pdb", 'w') as pdb_file:

        # Write the opening line so vorpy knows what to expect
        pdb_file.write("REMARK coarsify file - Scheme = {}, Thermal Cushion = {}\n"
                       .format(sys.scheme, sys.therm_cush))

        # Go through each atom in the system
        for i, ball in enumerate(sys.balls):
            # Skip any sol with way too big of a radius
            if ball.chain.name in ['SOL', 'Z', 'X', 'W'] and ball.rad > 10:
                print('\nThe Ball {}, {}, {} (a part of SOL) has a radius of {} and was skipped'
                      .format(ball.name, ball.index, ball.seq, ball.rad))
                continue
            # Set the ball index
            if ball.index is None:
                ball.index = i
            # Make sure the ball index is set correctly
            ball.index = str(int(ball.index) % 100000)
            atom_name = ball.element
            if ball.name == 'W':
                ball.name = 'SOL'
            res_name = ball.name
            chain = ball.chain.name
            if chain in ['SOL', 'Z', 'X', 'W']:
                chain = "0"
            res_seq = ball.seq
            x, y, z = ball.loc
            occ = ball.mass
            if ball.sub_section == 'bb':
                charge = 1
            elif ball.sub_section == 'sc':
                charge = 2
            elif ball.sub_section == 'sugar':
                charge = 1
            elif ball.sub_section == 'nbase':
                charge = 2
            elif ball.sub_section == 'phosphate':
                charge = 3
            else:
                charge = 0
            tfact = ball.rad
            elem = ball.element
            if len(ball.atoms) == 1:
                elem = ball.atoms[0]['element']
                atom_name = ball.atoms[0]['name']
            # Write the atom information
            pdb_file.write(pdb_line(ser_num=ball.index, name=atom_name, res_name=res_name, chain=chain, res_seq=res_seq,
                                    x=x, y=y, z=z, occ=occ, tfact=tfact, elem=elem, charge=charge))
    # Change back to the starting directory
    os.chdir(start_dir)


def write_balls(sys, file_name=None):
    """
    Writes a ball.txt file with the pure ball location and radii information along with information about each ball that
    has been cg'd
    """
    # Set up the file name
    if file_name is None:
        file_name = sys.name + '.txt'
    # Open the file
    with open(sys.dir + '/' + file_name, 'w') as wf:
        # Write the first line
        wf.write(f"System: {sys.name}, Info order: Index x y z r # chain res_name res_seq ball_subsection mass element"
                 f" name\n")
        # Loop through the balls in the file
        for i, ball in enumerate(sys.balls):
            # Skip any sol with way too big of a radius
            if ball.chain.name in ['SOL', 'Z', 'X', 'W'] and ball.rad > 10:
                print('\nThe Ball {}, {}, {} (a part of SOL) has a radius of {} and was skipped'
                      .format(ball.name, ball.index, ball.seq, ball.rad))
                continue
            # Set the ball index
            if ball.index is None:
                ball.index = i
            # Make sure the ball index is set correctly
            ball.index = str(int(ball.index) % 100000)
            atom_name = ball.element
            if ball.name == 'W':
                ball.name = 'SOL'
            res_name = ball.name
            chain = ball.chain.name
            if chain in ['SOL', 'Z', 'X', 'W']:
                chain = "0"
            res_seq = ball.seq
            x, y, z = ball.loc
            elem = ball.element
            if len(ball.atoms) == 1:
                elem = ball.atoms[0]['element']
                atom_name = ball.atoms[0]['name']
            # Write the ball information
            wf.write(f"{ball.index} {x} {y} {z} {ball.rad} # {chain} {res_name} {res_seq} {ball.sub_section} "
                     f"{ball.mass} {elem} {atom_name}\n")


def write_pymol_atoms(sys, set_sol=True, file_name=None):
    start_dir = os.getcwd()
    os.chdir(sys.dir)
    if file_name is None:
        file_name = 'set_atoms.pml'
    if not set_sol:
        file_name = file_name[:-4] + '_nosol.pml'
    # Create the file
    with open(file_name, 'w') as file:
        for ball in sys.balls:
            if not set_sol and ball.name.lower() in {'sol', 'hoh'}:
                continue
            name = ball.element
            if len(ball.atoms) == 1:
                name = ball.atoms[0]['name']
            file.write("alter ({} and resn {} and resi {} and name {}), vdw={}\n"
                       .format(sys.name, ball.name, ball.seq, name, round(ball.rad, 3)))
        # Rebuild the system
        file.write("\nrebuild")
    os.chdir(start_dir)


def color_pymol_balls(sys, bb_sc=False):
    """
    A function that creates a script for the pymol balls to be colored by residue or by side chain backbone
    """
    start_dir = os.getcwd()
    os.chdir(sys.dir)
    file_name = 'set_atom_colors.pml'
    sc_bb_colors = {'bb': '0xB8B8B8', 'phosphate': '0xFFA500', 'sugar': '0xB8B8B8'}
    # Create the file
    with open(file_name, 'w') as file:
        for ball in sys.balls:
            if len(ball.atoms) == 1:
                continue

            if not bb_sc:
                file.write("color {}, ({} and resn {} and resi {} and name {})\n".format(ball.residues[0].color, sys.name, ball.name, ball.seq,
                                                                                  ball.element))
            else:
                if ball.sub_section == 'nbase':
                    file.write(
                        "color {}, ({} and resn {} and resi {} and name {})\n".format(ball.residues[0].color, sys.name,
                                                                               ball.name, ball.seq, ball.element))
                elif ball.sub_section is None or ball.sub_section == 'sc':
                    file.write(
                        "color {}, ({} and resn {} and resi {} and name {})\n".format(ball.residues[0].color, sys.name,
                                                                               ball.name, ball.seq, ball.element))
                else:
                    file.write(
                        "color {}, ({} and resn {} and resi {} and name {})\n".format(sc_bb_colors[ball.sub_section], sys.name,
                                                                               ball.name, ball.seq, ball.element))
        # Rebuild the system
        file.write("\nrebuild")
    os.chdir(start_dir)


def set_sys_dir(sys, dir_name=None, root_dir=None):
    """
    Sets the directory for the output data. If the directory exists add 1 to the end number
    :param sys: System to assign the output directory to
    :param dir_name: Name for the directory
    :param root_dir:
    :return:
    """

    # Make sure a user_data path exists
    if sys.vpy_dir is not None and not os.path.exists(sys.vpy_dir + "/Data/user_data"):
        os.mkdir(sys.vpy_dir + "/Data/user_data")
    elif sys.vpy_dir is None and not os.path.exists("./Data/user_data"):
        if not os.path.exists('./Data'):
            os.mkdir(os.path.abspath('.') + '/Data/user_data')
        else:
            os.mkdir(os.path.abspath('./Data') + '/user_data')

    # If no outer directory was specified use the directory outside the current one
    if dir_name is None:
        if root_dir is not None:
            dir_name = root_dir + '/' + sys.name
        elif sys.vpy_dir is not None:

            dir_name = sys.vpy_dir + "/Data/user_data/" + sys.name
        else:
            dir_name = os.getcwd() + "/Data/user_data/" + sys.name
    print(dir_name)
    # Catch for existing directories. Keep trying out directories until one doesn't exist
    i = 0
    while True:
        # Try creating the directory with the system name + the current i_string
        try:
            # Create a string variable for the incrementing variable
            i_str = '_' + str(i)
            # If no file with the system name exists change the string to empty
            if i == 0:
                i_str = ""
            # Try to create the directory
            os.mkdir(dir_name + i_str)
            break
        # If the file exists increment the counter and try creating the directory again
        except FileExistsError:
            i += 1
    # Set the output directory for the system
    sys.dir = dir_name + i_str
