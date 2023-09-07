import os.path as path
from pandas import DataFrame
import numpy as np
import os
from System.objects import make_atom, Residue, Sol, Chain
from System.calcs import calc_com, calc_dist, get_radius, pdb_line


class System:
    def __init__(self, file, atoms=None, output_directory=None, root_dir=None, print_actions=False, residues=None,
                 chains=None, segments=None, output=True):
        """
        Class used to import files of all types and return a System
        :param file: Base system file address
        :param atoms: List holding the atom objects
        :param output_directory: Directory for export files to be output to
        """

        # Names
        self.name = None                    # Name                :   Name describing the system

        # Loadable objects
        self.atoms = atoms                  # Atoms               :   List holding the atom objects
        self.residues = residues            # Residues            :   List of residues (lists of atoms)
        self.chains = chains
        self.segments = segments
        self.sol = None                     # Solution            :   List of solution molecules (lists of atoms)

        self.radii = my_radii               # Radii               :   List of atomic radii
        self.special_radii = special_radii  # Special Radii       :   List of special radius situations. Helpful for gro
        self.decimals = None                # Decimals            :   Decimals setting for the whole system

        # Set up the file attributes
        self.data = None                    # Data                :   Additional data provided by the base file
        self.base_file = file               # Base file           :   Primary file address
        self.dir = output_directory         # Output Directory    :   Output directory for the export files
        self.vpy_dir = root_dir             # Vorpy Directory     :   Directory that vorpy is running out of
        self.max_atom_rad = 0               # Max atom rad        :   Largest radius of the system for reference

        # Gui
        self.print_actions = print_actions  # Print actions Bool  :   Tells the system to print or not

        # Run the processes
        self.read_pdb()
        self.print_info()
        self.coarsify()
        if output:
            self.set_dir()
            self.write_pdb()
            self.set_pymol_atoms()

    def read_pdb(self):
        """
        Interprets pdb data into a system of atom objects
        :param self: System to add the pdb information to
        :return: list of tuples of locations and radii
        """
        # Check to see if the file is provided and use the base file if not
        file = self.base_file

        # Get the file information and make sure to close the file when done
        with open(file, 'r') as f:
            my_file = f.readlines()
        # Add the system name and reset the atoms and data lists
        self.name = path.basename(self.base_file)[:-4] + '_coarse'
        # Set up the atom and the data lists
        atoms, data, atom_count = [], [], 0
        self.chains, self.residues = [], []
        chains, resids = {}, {}
        # Go through each line in the file and check if the first word is the word we are looking for
        reset_checker = 0
        for i in range(len(my_file)):
            # Check to make sure the line isn't empty
            if len(my_file[i]) == 0:
                continue
            # Pull the file line and first word
            line = my_file[i]
            word = line[:4].lower()
            # Check to see if the line is an atom line. If the line is not an atom line store the other data
            if line and word != 'atom':
                data.append(my_file[i].split())
                continue
            # Check for the "m" situation
            if line[76:78] == ' M':
                continue
            # Get the residue sequence of the atom
            res_seq = int(line[22:26])
            if line[22:26] == '    ':
                res_seq = 0
            # Create the atom
            atom = make_atom(location=np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])]),
                             system=self, element=line[76:78].strip(), res_seq=res_seq, name=line[12:16].strip(),
                             seg_id=line[72:76], index=atom_count, chain=line[21], residue=line[17:20].strip())

            atom_count += 1
            # Collect the chain type for each atom
            if atom['chain'] == ' ':
                if atom['residue'].lower() in {'sol', 'hoh', 'sod'}:
                    atom['chain'] = 'Z'
                elif atom['residue'].lower() in {'cl', 'mg', 'na', 'k'} and 'Z' in chains:
                    atom['chain'] = 'X'
                else:
                    atom['chain'] = 'A'
            # Create the chain and residue dictionaries
            res_name = atom['chain'] + atom['residue'] + str(atom['res_seq']) + str(reset_checker)
            # If the chain has been made before
            if atom['chain'] in chains:
                # Get the chain from the dictionary and add the atom
                my_chn = chains[atom['chain']]
                my_chn.add_atom(atom['num'])
                atom['chn'] = my_chn
            # Create the chain
            else:
                # If the chain is the sol chain
                if atom['chain'] == 'Z':
                    my_chn = Sol(atoms=[atom['num']], residues=[], name=atom['chain'], sys=self)
                    self.sol = my_chn
                # If the chain is not sol create a regular chain object
                else:
                    my_chn = Chain(atoms=[atom['num']], residues=[], name=atom['chain'], sys=self)
                    self.chains.append(my_chn)
                # Set the chain in the dictionary and give the atom it's chain
                chains[atom['chain']] = my_chn
                atom['chn'] = my_chn

            # Assign the atoms and create the residues
            if res_name in resids:
                my_res = resids[res_name]
                my_res.atoms.append(atom)
            else:
                my_res = Residue(sys=self, atoms=[atom], name=atom['residue'], sequence=atom['res_seq'], chain=atom['chn'])
                atom['chn'].residues.append(my_res)
                resids[res_name] = my_res
                self.residues.append(my_res)
            # Assign the residue to the atom
            atom['res'] = my_res

            # Assign the radius
            atom['rad'] = get_radius(atom)
            # If the residue numbers roll over reset the name of the residue to distinguish between the residues
            atoms.append(atom)
            if res_seq == 9999:
                reset_checker += 1
        # Set the colors for the residues based off the default colors for set elements
        res_colors = {'ALA': 'H', 'ARG': "He", 'ASN': 'Li', 'ASP': 'Be', 'ASX': 'B', 'CYS': 'C', 'GLN': 'F', 'GLU': 'O',
                      'GLX': 'S', 'GLY': 'Cl', 'HIS': 'Ar', 'ILE': 'Na', 'LEU': 'Mg', 'LYS': 'Mg', 'MET': 'Al',
                      'PHE': 'Si', 'PRO': 'P', 'SER': 'S', 'THR': 'Cl', 'TRP': 'Ar', 'TYR': 'K', 'VAL': 'Ca',
                      'SOL': 'Ti', 'DA': 'N', 'DC': 'O', 'DG': 'F', 'DT': 'S'}
        # Set the residues' colors and let the default go to Titanium (Grey)
        for res in self.residues:
            if res.name in res_colors:
                res.elem_col = res_colors[res.name]
            else:
                res.elem_col = 'Ti'
        # Set the atoms and the data
        self.atoms, self.data = DataFrame(atoms), data

    def print_info(self):
        """
        Prints the information for the loaded system
        """
        # Count the number of atoms, residues and chains and print their characteristics
        atoms_var = str(len(self.atoms)) + " Atoms"
        resids_var = str(len(self.residues)) + " Residues"
        chains_var = str(len(self.chains)) + " Chains: " + ", ".join(["{} - {} atoms, {} residues"
                            .format(_.name, len(_.atoms), len(_.residues)) for _ in self.chains])
        # Create the variable for the SOL
        sol_var = ""
        if self.sol is not None:
            sol_var = self.sol.name + " - " + str(len(self.sol.residues)) + " residues"
        # Print everything
        print(atoms_var, resids_var, chains_var, sol_var)

    def coarsify(self):
        """
        Main coarsify function. Calculates radii and location for residues
        """
        # Loop through the residues in the system
        for res in self.residues:
            # Calculate the center of mass for the atoms in a residue
            res.loc = calc_com([_['loc'] for _ in res.atoms])
            # Find the maximum of the summations of atom radii and the distance to residue com
            res.rad = max([calc_dist(res.loc, atom['loc']) + atom['rad'] for atom in res.atoms])

    def set_dir(self, dir_name=None):
        """
        Sets the directory for the output data. If the directory exists add 1 to the end number
        :param self: System to assign the output directory to
        :param dir_name: Name for the directory
        """
        if not os.path.exists("./Data/user_data"):
            os.mkdir("./Data/user_data")
        # If no outer directory was specified use the directory outside the current one
        if dir_name is None:
            if self.vpy_dir is not None:
                dir_name = self.vpy_dir + "/Data/user_data/" + self.name
            else:
                dir_name = os.getcwd() + "/Data/user_data/" + self.name
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
        self.dir = dir_name + i_str

    def write_pdb(self):
        """
        Creates a pdb file type in the current working directory
        :param self: System object used for writing the whole pbd file
        :return: Writes a pdb file for the set of atoms
        """
        # Catch empty atoms cases
        if self.residues is None or len(self.residues) == 0:
            return
        # Make note of the starting directory
        start_dir = os.getcwd()
        # Change to the specified directory
        os.chdir(self.dir)

        # Open the file for writing
        with open(self.name + ".pdb", 'w') as pdb_file:

            # Write the opening line so vorpy knows what to expect
            pdb_file.write("REMARK coarsify file\n")

            # Go through each atom in the system
            for i, res in enumerate(self.residues):

                atom_name = res.elem_col
                res_name = res.name
                chain = res.chain.name
                if chain == 'SOL' or chain == 'Z':
                    chain = " "
                res_seq = res.seq
                x, y, z = res.loc
                occ = 1
                tfact = res.rad
                elem = res.elem_col
                # Write the atom information
                pdb_file.write(pdb_line(ser_num=i, name=atom_name, res_name=res_name, chain=chain, res_seq=res_seq,
                                        x=x, y=y, z=z, occ=occ, tfact=tfact, elem=elem))
        # Change back to the starting directory
        os.chdir(start_dir)

    def set_pymol_atoms(self):
        """
        Creates a script to set the radii of the spheres in pymol
        """
        start_dir = os.getcwd()
        os.chdir(self.dir)
        # Create the file
        with open('set_atoms.pml', 'w') as file:
            for res in self.residues:
                res_str = "residue {} ".format(res.seq) if res != "" else ""
                file.write("alter ({}name {}), vdw={}\n".format(res_str, res.name, res.rad))
            # Rebuild the system
            file.write("\nrebuild")
        os.chdir(start_dir)


##################################################### Atomic Radii #####################################################


my_radii = {'h': 1.30, 'he': 1.40, 'li': 0.76, 'be': 0.45, 'b': 1.92, 'c': 1.80, 'n': 1.60, 'o': 1.50, 'f': 1.33,
            'ne': 1.54, 'na': 1.02, 'mg': 0.72, 'al': 0.60, 'si': 2.10, 'p': 1.90, 's': 1.90, 'cl': 1.81, 'ar': 1.88,
            'k': 1.38, 'ca': 1.00, 'ga': 0.62, 'ge': 0.73, 'as': 0.58, 'se': 1.90, 'br': 1.83, 'kr': 2.02, 'rb': 1.52,
            'sr': 1.18, 'in': 1.93, 'sn': 2.17, 'sb': 2.06, 'te': 2.06, 'i': 2.20, 'xe': 2.16, 'cs': 1.67, 'ba': 1.35,
            'tl': 1.96, 'pb': 2.02, 'bi': 2.07, 'po': 1.97, 'at': 2.02, 'rn': 2.20, 'fr': 3.48, 'ra': 2.83, '': 1.80}
special_radii = {''   : {'C': 1.75, 'CA': 1.90, 'N': 1.70, 'O': 1.49, 'F': 1.33, 'CL': 1.81, 'BR': 1.96, 'I': 2.20},
                 'ALA': {'CB': 1.92},
                 'ARB': {'CB': 1.91, 'CD': 1.88, 'CG': 1.92, 'CZ': 1.80, 'NE': 1.62, 'NH1': 1.62, 'NH2': 1.67},
                 'ASN': {'CB': 1.91, 'CG': 1.81, 'ND2': 1.62, 'OD1': 1.52},
                 'ASP': {'CB': 1.91, 'CG': 1.76, 'OD1': 1.49, 'OD2': 1.49},
                 'CYS': {'CB': 1.91, 'S': 1.88},
                 'GLN': {'CB': 1.91, 'CD': 1.81, 'CG': 1.80, 'NE2': 1.62, 'OE1': 1.52},
                 'GLU': {'CB': 1.91, 'CD': 1.76, 'CG': 1.88, 'OE1': 1.49, 'OE2': 1.49},
                 'HIS': {'CB': 1.91, 'CD': 1.74, 'CE': 1.74, 'CG': 1.80, 'ND1': 1.60, 'ND2': 1.60},
                 'ILE': {'CB': 2.01, 'CD1': 1.92, 'CG1': 1.92, 'CG2': 1.92},
                 'LEU': {'CB': 1.91, 'CD1': 1.92, 'CD2': 1.92, 'CG': 2.01},
                 'LYS': {'CB': 1.91, 'CD': 1.92, 'CE': 1.88, 'CG': 1.92, 'NZ': 1.67},
                 'MET': {'CB': 1.91, 'CE': 1.80, 'CG': 1.92, 'S': 1.94},
                 'PHE': {'CB': 1.91, 'CD': 1.82, 'CE': 1.82, 'CG': 1.74, 'CZ': 1.82},
                 'PRO': {'CB': 1.91, 'CD': 1.92, 'CG': 1.92},
                 'SER': {'CB': 1.91, 'OG': 1.54},
                 'THR': {'CB': 2.01, 'CG2': 1.92, 'OG': 1.54},
                 'TRP': {'CB': 1.91, 'CD': 1.82, 'CE': 1.82, 'CE2': 1.74, 'CG': 1.74, 'CH': 1.82, 'CZ': 1.82, 'NE1': 1.66},
                 'TYR': {'CB': 1.91, 'CD': 1.82, 'CE': 1.82, 'CG': 1.74, 'CZ': 1.80, 'OH': 1.54},
                 'VAL': {'CB': 2.01, 'CG1': 1.92, 'CG2': 1.92}}