from System.calcs import pdb_line
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
        pdb_file.write("REMARK coarsify file\n")

        # Go through each atom in the system
        for i, ball in enumerate(sys.balls):
            if ball.index is None:
                ball.index = i
            ball.index %= 100000
            atom_name = ball.element
            if ball.name == 'W':
                ball.name = 'SOL'
            res_name = ball.name
            chain = ball.chain.name
            if chain in ['SOL', 'Z', 'X', 'W']:
                chain = " "
            res_seq = ball.seq
            x, y, z = ball.loc
            occ = 0
            if ball.sub_section == 'bb':
                occ = 1
            elif ball.sub_section == 'sc':
                occ = 2
            tfact = ball.rad
            elem = ball.element
            # Write the atom information
            pdb_file.write(pdb_line(ser_num=ball.index, name=atom_name, res_name=res_name, chain=chain, res_seq=res_seq,
                                    x=x, y=y, z=z, occ=occ, tfact=tfact, elem=elem))
    # Change back to the starting directory
    os.chdir(start_dir)


def write_pymol_atoms(sys, set_sol=True):
    start_dir = os.getcwd()
    os.chdir(sys.dir)
    if set_sol:
        file_name = 'set_atoms_all.pml'
    else:
        file_name = 'set_atoms_no_sol.pml'
    # Create the file
    with open(file_name, 'w') as file:
        for ball in sys.balls:
            if not set_sol and ball.name.lower() == 'sol':
                continue
            res_str = "residue {} ".format(ball.seq) if ball != "" else ""
            file.write("alter ({}name {}), vdw={}\n".format(res_str, ball.name, round(ball.rad, 3)))
        # Rebuild the system
        file.write("\nrebuild")
    os.chdir(start_dir)
