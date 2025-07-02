import os
from os import path
import shutil
import tkinter as tk
from tkinter import filedialog
from System.sys_funcs.input import read_pdb
from System.sys_funcs.output import set_dir, write_pdb, write_pymol_atoms, color_pymol_balls, set_sys_dir, write_balls
from System.schemes.martini import coarsify_martini
from System.schemes.basic import coarsify
from GUI import settings_gui


class System:
    def __init__(self, file=None, atoms=None, output_directory=None, root_dir=None, print_actions=False, residues=None,
                 chains=None, segments=None, output=True, scheme=None, thermal_cushion=0.0, include_h=True,
                 mass_weighted=True, sc_bb=None, color_scheme='Shapely', all_methods=False):
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
        self.mass_weighted = mass_weighted
        self.sc_bb = sc_bb
        self.all = all_methods

        # Loadable objects
        self.atoms = atoms                  # Atoms               :   List holding the atom objects
        self.residues = residues            # Residues            :   List of residues (lists of atoms)
        self.chains = chains
        self.segments = segments
        self.sol = None                     # Solution            :   List of solution molecules (lists of atoms)

        self.radii = element_radii               # Radii               :   List of atomic radii
        self.masses = my_masses
        self.special_radii = special_radii  # Special Radii       :   List of special radius situations. Helpful for gro
        self.aminos = amino_acids
        self.amino_bbs = amino_bbs
        self.amino_scs = amino_scs
        self.amino_ignores = []
        self.nucleics = nucleic_acids
        self.nucleic_sugrs = nucleic_nbase
        self.nucleic_pphte = nucleic_pphte
        self.nucleic_nbase = nucleic_sugr
        self.nucleic_ignores = []
        self.decimals = None                # Decimals            :   Decimals setting for the whole system
        self.color_scheme = color_scheme

        self.balls = None                   # Balls               :   Output for the program

        # Set up the file attributes
        self.data = None                    # Data                :   Additional data provided by the base file
        self.base_file = file               # Base file           :   Primary file address
        self.dir = output_directory         # Output Directory    :   Output directory for the export files
        self.vpy_dir = os.getcwd()          # Vorpy Directory     :   Directory that vorpy is running out of
        self.max_atom_rad = 0               # Max atom rad        :   Largest radius of the system for reference

        # Print Actions
        self.print_actions = print_actions  # Print actions Bool  :   Tells the system to print or not

        # Run the processes
        my_vals = settings_gui()
        if my_vals['cg method'] == 'All Schemes':
            self.run_all_schemes(my_vals)
            return

        self.get_vals(my_vals)
        self.set_name()
        self.read_pdb()
        self.print_info()
        self.coarsify()
        if self.dir is None or not os.path.exists(self.dir):
            self.set_sys_dir()
        self.output(self.dir)

    def get_vals(self, my_vals):

        self.include_h = my_vals['include h']
        self.base_file = my_vals['input file']
        self.scheme = my_vals['cg method']
        self.mass_weighted = my_vals['mass weighted']
        self.therm_cush = my_vals['thermal cushion']
        self.sc_bb = my_vals['sc bb']
        self.dir = my_vals['output folder']

    def set_name(self):
        # Set up the sc_bb var
        sc_bb = ''
        if self.sc_bb:
            sc_bb = '_Split'
        mw = ''
        if self.mass_weighted and self.scheme == 'Average Distance':
            mw = '_MW'
        # Add the system name and reset the atoms and data lists
        name = path.basename(self.base_file)[:-4] + '_' + self.scheme + sc_bb + mw
        # Split and rejoin the name
        self.name = '_'.join(name.split(' '))

    def read_pdb(self):
        """
        Interprets pdb data into a system of atom objects
        :param self: System to add the pdb information to
        :return: list of tuples of locations and radii
        """
        read_pdb(self)

    def run_all_schemes(self, my_vals):
        if self.dir is not None:
            rooty_tooty = self.dir
        else:
            rooty_tooty = self.vpy_dir + '/Data/user_data/'
        for scheme in ['Encapsulate', 'Average Distance']:
            for sc_bb_val in [True, False]:
                for mass_weight_val in [True, False]:
                    if mass_weight_val and scheme == 'Encapsulate':
                        continue
                    # Reset everything
                    self.balls, self.name, self.dir, self.scheme = None, None, None, None
                    # Set my_vals values
                    my_vals['cg method'], my_vals['sc bb'], my_vals['mass weighted'] = scheme, sc_bb_val, mass_weight_val
                    self.get_vals(my_vals)
                    self.set_name()
                    self.read_pdb()
                    self.print_info()
                    self.coarsify()
                    self.set_sys_dir(root_dir=rooty_tooty)
                    print(self.dir)
                    self.output(self.dir)

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
        if self.scheme == 'Martini':
            coarsify_martini(self)
        else:
            coarsify(self)

    def set_sys_dir(self, my_dir=None, root_dir=None):
        set_sys_dir(self, my_dir, root_dir=root_dir)

    def output(self, my_dir=None):
        """
        Outputs the information for the coarsified data
        """
        if my_dir is None:
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
        # Write the atom colors
        color_pymol_balls(self, self.sc_bb)
        # Copy the original pdb over
        shutil.copy2(self.base_file, self.dir + '/' + path.basename(self.base_file[:-4]) + '_base.pdb')
        # Create the balls file
        write_balls(self)

    def set_dir(self, dir_name=None):
        """
        Sets the directory for the output data. If the directory exists add 1 to the end number
        :param self: System to assign the output directory to
        :param dir_name: Name for the directory
        """
        set_dir(self, dir_name=dir_name)


##################################################### Atomic Radii #####################################################


my_masses = {'h': 1.008, 'he': 4.003, 'li': 6.941, 'be': 9.012, 'b': 10.811, 'c': 12.011, 'n': 14.007, 'o': 15.999,
             'f': 18.998,'ne': 20.180, 'na': 22.990, 'mg': 24.305, 'al': 26.982, 'si': 28.086, 'p': 30.974,
             's': 32.066, 'cl': 35.453, 'ar': 39.948, 'k': 39.098, 'ca': 40.078, 'ga': 69.723, 'ge': 72.631,
             'as': 74.922, 'se': 78.971, 'br': 79.904, 'kr': 83.798, 'rb': 85.468, 'sr': 87.62, 'in': 114.818,
             'sn': 118.711, 'sb': 121.760, 'te': 27.6, 'i': 126.904, 'xe': 131.293, 'cs': 132.905, 'ba': 137.328,
             'tl': 204.383, 'pb': 207.2, 'bi': 208.980, 'po': 208.982, 'at': 209.987, 'rn': 222.018, 'fr': 223.020,
             'ra': 226.025, '': 1.80, 'W': 4.1}

element_radii = {
    'H': 1.30, 'HE': 1.40, 'LI': 0.76, 'BE': 0.45, 'B': 1.92, 'C': 1.80, 'N': 1.60, 'O': 1.50, 'F': 1.33,
    'NE': 1.54, 'NA': 1.02, 'MG': 0.72, 'AL': 0.60, 'SI': 2.10, 'P': 1.90, 'S': 1.90, 'CL': 1.81, 'AR': 1.88,
    'K': 1.38, 'CA': 1.00, 'GA': 0.62, 'GE': 0.73, 'AS': 0.58, 'SE': 1.90, 'BR': 1.83, 'KR': 2.02, 'RB': 1.52,
    'SR': 1.18, 'IN': 1.93, 'SN': 2.17, 'SB': 2.06, 'TE': 2.06, 'I': 2.20, 'XE': 2.16, 'CS': 1.67, 'BA': 1.35,
    'TL': 1.96, 'PB': 2.02, 'BI': 2.07, 'PO': 1.97, 'AT': 2.02, 'RN': 2.20, 'FR': 3.48, 'RA': 2.83, 'ZN': 1.39
}

special_radii = {
    'ALA': {
        'C': 1.8, 'CA': 1.8, 'CB': 1.92, 'H': 1.3, 'HA': 1.3, 'HB1': 1.3, 'HB2': 1.3, 'HB3': 1.3, 'N': 1.6, 'O': 1.5, 'OC1': 1.5, 'OC2': 1.5
    },
    'ARG': {
        'C': 1.8, 'CA': 1.8, 'CB': 1.91, 'CD': 1.88, 'CG': 1.92, 'CZ': 1.8, 'H': 1.3, 'HA': 1.3, 'HB1': 1.3, 'HB2': 1.3, 'HB3': 1.3, 'HD1': 1.3, 'HD2': 1.3, 'HD3': 1.3, 'HE': 1.3, 'HG1': 1.3, 'HG2': 1.3, 'HH11': 1.3, 'HH12': 1.3, 'HH21': 1.3, 'HH22': 1.3, 'N': 1.6, 'NE': 1.62, 'NH1': 1.62, 'NH2': 1.67, 'O': 1.5
    },
    'ASN': {
        'C': 1.8, 'CA': 1.8, 'CB': 1.91, 'CG': 1.81, 'H': 1.3, 'HA': 1.3, 'HB1': 1.3, 'HB2': 1.3, 'HD21': 1.3, 'HD22': 1.3, 'N': 1.6, 'ND2': 1.62, 'O': 1.5, 'OD1': 1.52
    },
    'ASP': {
        'C': 1.8, 'CA': 1.8, 'CB': 1.91, 'CG': 1.76, 'H': 1.3, 'HA': 1.3, 'HB1': 1.3, 'HB2': 1.3, 'N': 1.6, 'O': 1.5, 'OD1': 1.49, 'OD2': 1.49
    },
    'CYS': {
        'C': 1.8, 'CA': 1.8, 'CB': 1.91, 'H': 1.3, 'HA': 1.3, 'HB1': 1.3, 'HB2': 1.3, 'HG': 1.3, 'N': 1.6, 'O': 1.5, 'S': 1.9, 'SG': 1.9
    },
    'GLN': {
        'C': 1.8, 'CA': 1.8, 'CB': 1.91, 'CD': 1.81, 'CG': 1.8, 'H': 1.3, 'HA': 1.3, 'HB1': 1.3, 'HB2': 1.3, 'HE21': 1.62, 'HE22': 1.62, 'HG1': 1.3, 'HG2': 1.3, 'N': 1.6, 'NE2': 1.62, 'O': 1.5, 'OE1': 1.52
    },
    'GLU': {
        'C': 1.8, 'CA': 1.8, 'CB': 1.91, 'CD': 1.76, 'CG': 1.88, 'H': 1.3, 'HA': 1.3, 'HB1': 1.3, 'HB2': 1.3, 'HG1': 1.3, 'HG2': 1.3, 'N': 1.6, 'O': 1.5, 'OE1': 1.49, 'OE2': 1.49
    },
    'GLY': {
        'C': 1.8, 'CA': 1.8, 'H': 1.3, 'HA1': 1.3, 'HA2': 1.3, 'N': 1.6, 'O': 1.5, 'OC1': 1.5, 'OC2': 1.5
    },
    'HIS': {
        'C': 1.8, 'CA': 1.8, 'CB': 1.91, 'CD': 1.74, 'CE': 1.74, 'CD2': 1.74, 'CE1': 1.74, 'CG': 1.8, 'H': 1.3, 'HA': 1.3, 'HB1': 1.3, 'HB2': 1.3, 'HD2': 1.3, 'HE1': 1.3, 'HE2': 1.3, 'N': 1.6, 'ND1': 1.6, 'ND2': 1.6, 'NE2': 1.6, 'O': 1.5
    },
    'ILE': {
        'C': 1.8, 'CA': 1.8, 'CB': 2.01, 'CD': 1.92, 'CD1': 1.92, 'CG1': 1.92, 'CG2': 1.92, 'H': 1.3, 'HA': 1.3, 'HB': 1.3, 'HD1': 1.3, 'HD2': 1.3, 'HD3': 1.3, 'HD11': 1.3, 'HD12': 1.3, 'HD13': 1.3, 'HG12': 1.3, 'HG13': 1.3, 'HG21': 1.3, 'HG22': 1.3, 'HG23': 1.3, 'N': 1.6, 'O': 1.5
    },
    'LEU': {
        'C': 1.8, 'CA': 1.8, 'CB': 1.91, 'CD1': 1.92, 'CD2': 1.92, 'CG': 2.01, 'H': 1.3, 'HA': 1.3, 'HB1': 1.3, 'HB2': 1.3, 'HD11': 1.3, 'HD12': 1.3, 'HD13': 1.3, 'HD21': 1.3, 'HD22': 1.3, 'HD23': 1.3, 'HG': 1.3, 'N': 1.6, 'O': 1.5
    },
    'LYS': {
        'C': 1.8, 'CA': 1.8, 'CB': 1.91, 'CD': 1.92, 'CE': 1.88, 'CG': 1.92, 'H': 1.3, 'HA': 1.3, 'HB1': 1.3, 'HB2': 1.3, 'HD1': 1.3, 'HD2': 1.3, 'HE1': 1.3, 'HE2': 1.3, 'HG1': 1.3, 'HG2': 1.3, 'HZ1': 1.67, 'HZ2': 1.67, 'HZ3': 1.67, 'N': 1.6, 'NZ': 1.67, 'O': 1.5
    },
    'MET': {
        'C': 1.8, 'CA': 1.8, 'CB': 1.91, 'CE': 1.8, 'CG': 1.92, 'H': 1.3, 'HA': 1.3, 'HB1': 1.3, 'HB2': 1.3, 'HE1': 1.3, 'HE2': 1.3, 'HE3': 1.3, 'HG1': 1.3, 'HG2': 1.3, 'N': 1.6, 'O': 1.5, 'SD': 1.94, 'S': 1.94
    },
    'PHE': {
        'C': 1.8, 'CA': 1.8, 'CB': 1.91, 'CD1': 1.82, 'CD2': 1.82, 'CD': 1.82, 'CE1': 1.82, 'CE2': 1.82, 'CG': 1.74, 'CZ': 1.82, 'H': 1.3, 'HA': 1.3, 'HB1': 1.3, 'HB2': 1.3, 'HD1': 1.3, 'HD2': 1.3, 'HE1': 1.3, 'HE2': 1.3, 'HZ': 1.3, 'N': 1.6, 'O': 1.5
    },
    'PRO': {
        'C': 1.8, 'CA': 1.8, 'CB': 1.91, 'CD': 1.92, 'CG': 1.92, 'H': 1.3, 'HA': 1.3, 'HB1': 1.3, 'HB2': 1.3, 'HD1': 1.3, 'HD2': 1.3, 'HG1': 1.3, 'HG2': 1.3, 'N': 1.6, 'O': 1.5
    },
    'SER': {
        'C': 1.8, 'CA': 1.8, 'CB': 1.91, 'H1': 1.3, 'H2': 1.3, 'H3': 1.3, 'H': 1.3, 'HA': 1.3, 'HB1': 1.3, 'HB2': 1.3, 'HG': 1.54, 'N': 1.6, 'O': 1.5, 'OG': 1.54
    },
    'THR': {
        'C': 1.8, 'CA': 1.8, 'CB': 2.01, 'CG2': 1.92, 'H': 1.3, 'HA': 1.3, 'HB': 1.3, 'HG1': 1.54, 'HG21': 1.3, 'HG22': 1.3, 'HG23': 1.3, 'N': 1.6, 'O': 1.5, 'OG1': 1.54
    },
    'TRP': {
        'C': 1.8, 'CA': 1.8, 'CB': 1.91, 'CD1': 1.82, 'CD2': 1.82, 'CD': 1.82, 'CE': 1.82, 'CE2': 1.74, 'CE3': 1.74, 'CG': 1.74, 'CH': 1.82, 'CH2': 1.82, 'CZ': 1.82, 'CZ1': 1.82, 'CZ2': 1.82, 'H': 1.3, 'HA': 1.3, 'HB1': 1.3, 'HB2': 1.3, 'HD1': 1.3, 'HE1': 1.62, 'HE3': 1.3, 'HZ1': 1.3, 'HZ2': 1.3, 'HH2': 1.3, 'N': 1.6, 'NE1': 1.66, 'O': 1.5
    },
    'TYR': {
        'C': 1.8, 'CA': 1.8, 'CB': 1.91, 'CD': 1.82, 'CD1': 1.82, 'CD2': 1.82, 'CE': 1.82, 'CE1': 1.82, 'CE2': 1.82, 'CG': 1.74, 'CZ': 1.8, 'H': 1.3, 'HA': 1.3, 'HB1': 1.3, 'HB2': 1.3, 'HD1': 1.3, 'HD2': 1.3, 'HE1': 1.3, 'HE2': 1.3, 'HH': 1.3, 'N': 1.6, 'O': 1.5, 'OH': 1.54
    },
    'VAL': {
        'C': 1.8, 'CA': 1.8, 'CB': 2.01, 'CG1': 1.92, 'CG2': 1.92, 'H': 1.3, 'HA': 1.3, 'HB': 1.3, 'HG11': 1.3, 'HG12': 1.3, 'HG13': 1.3, 'HG21': 1.3, 'HG22': 1.3, 'HG23': 1.3, 'N': 1.6, 'O': 1.5
    }
}

amino_acids = {'ALA', 'ARB', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER',
               'THR', 'TRP', 'TYR', 'VAL', 'GLY', 'ARG'}
amino_bbs = ['CA', 'HA', 'HA1', 'HA2', 'N', 'HN', 'H', 'C', 'O', 'OC1', 'OC2', 'OT1', 'OT2', 'H1', 'H2', 'H3']
amino_scs = ['CB', 'HB', 'HB1', 'HB2', 'HB3',
             'SD', 'CD', 'CD1', 'CD2', 'ND1', 'ND2', 'OD1', 'OD2', 'HD1', 'HD2', 'HD3', 'HD11', 'HD12', 'HD13', 'HD21', 'HD22', 'HD23'
             , 'CE', 'CE1', 'CE2', 'CE3', 'OE1', 'OE2', 'NE', 'NE1', 'NE2', 'HE', 'HE1', 'HE2', 'HE3', 'HE21', 'HE22',
             'CG', 'CG1', 'CG2', 'OG', 'SG', 'OG1', 'HG', 'HG1', 'HG2', 'HG11', 'HG12', 'HG13', 'HG21', 'HG22', 'HG23',
             'CH2', 'NH1', 'OH', 'HH', 'HH1', 'HH2', 'HH11', 'HH12', 'NH2', 'HH21', 'HH22',
             'NZ', 'CZ', 'CZ1', 'CZ2', 'CZ3', 'NZ', 'HZ', 'HZ1', 'HZ2', 'HZ3']

nucleic_acids = {'DT', 'DA', 'DG', 'DC', 'DU', 'U', 'G', 'A', 'T', 'C', 'GDP', 'OMC'}

nucleic_nbase = ['N1', 'N2', 'N3', 'N4', 'N5', 'N6', 'N7', 'N8', 'N9', 'C2', 'C4', 'C5', 'C6', 'C7', 'C8', 'O2', 'O4',
                 'O6', 'H2', 'H21', 'H22', 'H3', 'H4', 'H41', 'H42', 'H5', 'H6', 'H61', 'H62', 'H8', 'H71', 'H72', 'H73']
nucleic_sugr = ['O3\'', 'O5\'', 'C5\'', 'C4\'', 'O4\'', 'C3\'', 'C2\'', 'C1\'', 'O2\'', 'CM2', 'H1\'', 'H2\'', 'H2\'\'',
                'H3\'', 'H4\'', 'H5\'', 'H5\'\'', 'H3T', 'H5T', 'H2\'1', 'H2\'2', 'H5\'1', 'H5\'2']
nucleic_pphte = ['P', 'O1P', 'O2P', 'OP1', 'OP2', 'PA', 'PB', 'O1A', 'O1B', 'O2A', 'O2B', 'O3A', 'O3B']

bb_sc_colors = {**{_: 'r' for _ in amino_bbs}, **{_: 'y' for _ in amino_scs}, **{_: 'blue' for _ in nucleic_nbase},
                **{_: 'purple' for _ in nucleic_sugr}, **{_: 'maroon' for _ in nucleic_pphte}}

