from os import path
import tkinter as tk
from tkinter import filedialog
from System.sys_funcs.input import read_pdb
from System.sys_funcs.output import set_dir, write_pdb, write_pymol_atoms
from System.schemes.basic import coarsify_avg_dist, coarsify_encapsulate
from System.schemes.martini import coarsify_martini
from System.schemes.sc_bb import coarsify_sc_bb
from System.schemes.primo import coarsify_primo


class System:
    def __init__(self, file, atoms=None, output_directory=None, root_dir=None, print_actions=False, residues=None,
                 chains=None, segments=None, output=True, scheme=None, thermal_cushion=0.0, include_h=True,
                 mass_weighted=True):
        """
        Class used to import files of all types and return a System
        :param file: Base system file address
        :param atoms: List holding the atom objects
        :param output_directory: Directory for export files to be output to
        """

        # Names
        self.name = None                    # Name                :   Name describing the system
        self.scheme = scheme                # CG Scheme           :   The scheme by which the atoms are coarse grained-
        self.therm_cush = thermal_cushion   # Thermal Cushion     :   How much additional radius is given to each ball
        self.include_h = include_h

        # Loadable objects
        self.atoms = atoms                  # Atoms               :   List holding the atom objects
        self.residues = residues            # Residues            :   List of residues (lists of atoms)
        self.chains = chains
        self.segments = segments
        self.sol = None                     # Solution            :   List of solution molecules (lists of atoms)

        self.radii = my_radii               # Radii               :   List of atomic radii
        self.special_radii = special_radii  # Special Radii       :   List of special radius situations. Helpful for gro
        self.aminos = amino_acids
        self.amino_bbs = amino_bbs
        self.amino_scs = amino_scs
        self.amino_ignores = []
        self.nucleics = nucleic_acids
        self.nucleic_sugrs = nucleic_sugrs
        self.nucleic_pphte = nucleic_pphte
        self.nucleic_nbase = nucleic_nbase
        self.nucleic_ignores = []
        self.decimals = None                # Decimals            :   Decimals setting for the whole system

        self.balls = None                   # Balls               :   Output for the program

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
        self.coarsify(scheme=scheme, therm_cush=thermal_cushion, mass_weighted=mass_weighted)
        if output:
            self.output()

    def read_pdb(self):
        """
        Interprets pdb data into a system of atom objects
        :param self: System to add the pdb information to
        :return: list of tuples of locations and radii
        """
        read_pdb(self)

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

    def coarsify(self, scheme, therm_cush=0.0, mass_weighted=True):
        """
        Main coarsify function. Calculates radii and location for residues
        """
        if scheme == '1':
            coarsify_avg_dist(self, therm_cush, include_h=True, mass_weighted=mass_weighted)
        elif scheme == '2':
            coarsify_encapsulate(self, therm_cush, include_h=self.include_h)
        elif scheme == '3':
            coarsify_sc_bb(self, therm_cush=therm_cush, avg_dist=True, include_h=True, mass_weighted=mass_weighted)
        elif scheme == '4':
            coarsify_sc_bb(self, therm_cush=therm_cush, avg_dist=False, include_h=self.include_h)
        elif scheme == '5':
            coarsify_primo(self, therm_cush)
        elif scheme == '6':
            coarsify_martini(self, therm_cush)
        self.scheme = {'1': 'Average Distance', '2': 'Encapsulate', '3': 'Backbone Side-chain', '4': 'Alpha Carbon',
                       '5': 'Primo', '6': 'CG Martini'}[self.scheme]

    def output(self):
        """
        Outputs the information for the coarsified data
        """
        # Choose whether to output to user_data, the original folder, or some other folder
        output_dir_selection = input("Choose output location: \n\n"
                                     "1. Coarsify User Data     ->      ../coarsify/Data/user_data \n"
                                     "2. Original File Location ->  {}\n"
                                     "3. Other Directory        ->  (Opens Folder Dialog Window)\n"
                                     "  >>>  ".format(path.dirname(self.base_file)))
        # Set the output directory
        if output_dir_selection == '1':
            self.set_dir()
        elif output_dir_selection == '2':
            self.set_dir(path.dirname(self.base_file))
        else:
            root = tk.Tk()
            root.withdraw()
            self.dir = filedialog.askdirectory()
        # Write the pdb
        write_pdb(self)
        # Create the setting script for pymol
        write_pymol_atoms(self)
        write_pymol_atoms(self, set_sol=False)

    def set_dir(self, dir_name=None):
        """
        Sets the directory for the output data. If the directory exists add 1 to the end number
        :param self: System to assign the output directory to
        :param dir_name: Name for the directory
        """
        set_dir(self, dir_name=dir_name)


##################################################### Atomic Radii #####################################################


my_radii = {'h': 1.30, 'he': 1.40, 'li': 0.76, 'be': 0.45, 'b': 1.92, 'c': 1.80, 'n': 1.60, 'o': 1.50, 'f': 1.33,
            'ne': 1.54, 'na': 1.02, 'mg': 0.72, 'al': 0.60, 'si': 2.10, 'p': 1.90, 's': 1.90, 'cl': 1.81, 'ar': 1.88,
            'k': 1.38, 'ca': 1.00, 'ga': 0.62, 'ge': 0.73, 'as': 0.58, 'se': 1.90, 'br': 1.83, 'kr': 2.02, 'rb': 1.52,
            'sr': 1.18, 'in': 1.93, 'sn': 2.17, 'sb': 2.06, 'te': 2.06, 'i': 2.20, 'xe': 2.16, 'cs': 1.67, 'ba': 1.35,
            'tl': 1.96, 'pb': 2.02, 'bi': 2.07, 'po': 1.97, 'at': 2.02, 'rn': 2.20, 'fr': 3.48, 'ra': 2.83, '': 1.80,
            'W': 4.1}
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

amino_acids = {'ALA', 'ARB', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER',
               'THR', 'TRP', 'TYR', 'VAL', 'GLY', 'ARG'}
amino_bbs = ['CA', 'HA', 'HA1', 'HA2', 'N', 'HN', 'H', 'C', 'O', 'OC1', 'OC2', 'OT1', 'OT2', 'H1', 'H2', 'H3']
amino_scs = ['CB', 'HB', 'HB1', 'HB2', 'HB3',
             'SD', 'CD', 'CD1', 'CD2', 'ND1', 'ND2', 'OD1', 'OD2', 'HD1', 'HD2', 'HD3', 'HD11', 'HD12', 'HD13', 'HD21', 'HD22', 'HD23'
             , 'CE', 'CE1', 'CE2', 'CE3', 'OE1', 'OE2', 'NE', 'NE1', 'NE2', 'HE', 'HE1', 'HE2', 'HE3', 'HE21', 'HE22',
             'CG', 'CG1', 'CG2', 'OG', 'SG', 'OG1', 'HG', 'HG1', 'HG2', 'HG11', 'HG12', 'HG13', 'HG21', 'HG22', 'HG23',
             'CH2', 'NH1', 'OH', 'HH', 'HH1', 'HH2', 'HH11', 'HH12', 'NH2', 'HH21', 'HH22',
             'NZ', 'CZ', 'CZ1', 'CZ2', 'CZ3', 'NZ', 'HZ', 'HZ1', 'HZ2', 'HZ3']

nucleic_acids = {'DT', 'DA', 'DG', 'DC', 'DU', 'U', 'G', 'A', 'T', 'C'}
nucleic_sugrs = ['N1', 'N2', 'N3', 'N4', 'N5', 'N6', 'N7', 'N8', 'N9', 'C2', 'C4', 'C5', 'C6', 'C7', 'C8', 'O2', 'O4', 'O6']
nucleic_nbase = ['O3\'', 'O5\'', 'C5\'', 'C4\'', 'O4\'', 'C3\'', 'C2\'', 'C1\'']
nucleic_pphte = ['P', 'O1P', 'O2P', 'OP1', 'OP2']
