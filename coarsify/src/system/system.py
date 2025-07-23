import os
from os import path
import shutil
import tkinter as tk
from tkinter import filedialog
from coarsify.src.input.input import read_pdb
from coarsify.src.output.output import set_dir, write_pdb, write_pymol_atoms, color_pymol_balls, set_sys_dir, write_balls
from coarsify.src.schemes.martini import coarsify_martini
from coarsify.src.schemes.basic import coarsify
from coarsify.src.gui.GUI import settings_gui
from coarsify.src.tools import *


class System:
    def __init__(self, file=None, atoms=None, output_directory=None, root_dir=None, print_actions=False, residues=None,
                 chains=None, segments=None, output=True, scheme=None, thermal_cushion=0.0, include_h=True,
                 mass_weighted=True, sc_bb=None, color_scheme='Shapely', all_methods=False):
        """
        Class used to import files of all types and return a system
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

        self.radii = element_radii          # Radii               :   List of atomic radii
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
        if self.base_file is None:
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
        :param self: system to add the pdb information to
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
        # Print out the directory that the files have been exported to
        print("\nFiles exported to {}".format(self.dir))

    def set_dir(self, dir_name=None):
        """
        Sets the directory for the output data. If the directory exists add 1 to the end number
        :param self: system to assign the output directory to
        :param dir_name: Name for the directory
        """
        set_dir(self, dir_name=dir_name)
