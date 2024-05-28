import tkinter as tk
from tkinter import filedialog
from System.system import System


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    scheme = input("Choose coarsify scheme:\n\n"
                   " (1) Average Distance - Balls located at center of mass of residues with radius average distance of its constisuent atoms\n"
                   " (2) Encapsulate Residues - Balls that minimally encapsulate atoms in a residue\n"
                   " (3) Side Chain/Backbone (Average Distance) - Balls that encapsulate the backbone and the side chain for each residue\n"
                   " (4) Side Chain/Backbone (Encapsulate) - Balls that encapsulate the backbone and the side chain for each residue\n"
                   " (5) Primo - Primo Coarse graining scheme\n"
                   " (6) CG Martini - Pre-processed pdb files using cg-martini\n\n"
                   " >>>   ")
    thermal_cushion = input("Add a thermal cushion? (Y/N)")
    if thermal_cushion in ['Y', 'y', 'Yes', 'yes']:
        thermal_cushion = float(input('Enter cushion (in Angstroms) >>>'))
    else:
        thermal_cushion = 0.0
    mass_weighted = input('Mass Weighted')
    root = tk.Tk()
    root.withdraw()
    root.wm_attributes('-topmost', 1)
    file = filedialog.askopenfilename()
    sys = System(file=file, scheme=scheme, thermal_cushion=thermal_cushion, mass_weighted=mass_weighted)
