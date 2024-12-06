import tkinter as tk
from tkinter import ttk
from tkinter import filedialog

settings = False

"""
Thermal cushion, methodology, sc_bb, average distance, include H, input file, output file
"""


def settings_gui():
    # Function to collect values and print dictionary
    global settings

    def apply_values():
        global settings
        settings = {
            'input file': input_file_var.get(),
            "cg method": cg_method_var.get(),
            "mass weighted": mass_weighted_var.get(),
            "thermal cushion": float(thermal_cushion_var.get()),
            "sc bb": sc_bb_var.get(),
            "include h": include_h_var.get(),
            "output folder": output_folder_var.get(),
        }
        root.destroy()
        return settings  # Or replace this with return data if using within another function or script

    # Function to handle cancel
    def cancel():
        root.destroy()

    # Function to choose an input file
    def choose_input_file():
        choose_input_file_window = tk.Tk()
        choose_input_file_window.withdraw()
        choose_input_file_window.wm_attributes('-topmost', 1)
        input_file = filedialog.askopenfilename()
        input_file_var.set(input_file)
        choose_input_file_window.destroy()

    def choose_output_folder():
        choose_output_folder_window = tk.Tk()
        choose_output_folder_window.withdraw()
        choose_output_folder_window.wm_attributes('-topmost', 1)
        input_file = filedialog.askdirectory()
        output_folder_var.set(input_file)
        choose_output_folder_window.destroy()

    # Main window
    root = tk.Tk()
    root.title("Coarse Grain Settings")

    for i in range(3):
        root.grid_rowconfigure(i, weight=1, minsize=50)
        root.grid_columnconfigure(i, weight=1, minsize=50)

    input_file_var = tk.StringVar(value='Choose File')
    cg_method_var = tk.StringVar(value='Encapsulate')
    thermal_cushion_var = tk.StringVar(value='0.0')
    sc_bb_var = tk.BooleanVar(value=False)
    include_h_var = tk.BooleanVar(value=True)
    mass_weighted_var = tk.BooleanVar(value=True)
    output_folder_var = tk.StringVar(value='Choose Output Folder')

    # Choose File
    tk.Label(root, textvariable=input_file_var).grid(row=0, column=0, columnspan=2)
    tk.Button(root, text='Browse', command=choose_input_file).grid(row=0, column=2)

    # Methodology
    tk.Label(root, text='CG Method').grid(row=1, column=0)
    method_menu = ttk.Combobox(root, textvariable=cg_method_var, values=['Encapsulate', 'Average Distance', 'Martini', 'All Schemes'])
    method_menu.grid(row=1, column=1, columnspan=2)
    method_menu.current(0)

    # Thermal Cushion variable
    tk.Label(root, text='Thermal Cushion').grid(row=2, column=0)
    tk.Entry(root, textvariable=thermal_cushion_var).grid(row=2, column=1)
    tk.Label(root, text='\u212B').grid(row=2, column=2)

    # Sidechain Backbone or Full residue
    tk.Checkbutton(root, text='Split Residue?', variable=sc_bb_var).grid(row=3, column=0, columnspan=3)

    # Include Hydrogen checkbutton
    tk.Checkbutton(root, text='Include Hydrogens?', variable=include_h_var).grid(row=4, column=0, columnspan=3)

    # Mass Weighted button
    tk.Checkbutton(root, text='Mass Weighted', variable=mass_weighted_var).grid(row=5, column=0, columnspan=3)

    tk.Label(root, textvariable=output_folder_var).grid(row=6, column=0, columnspan=2)
    tk.Button(root, text='Browse', command=choose_output_folder).grid(row=6, column=2)

    # Buttons
    tk.Button(root, text="Apply", command=apply_values).grid(row=7, column=2, pady=10, padx=10)
    tk.Button(root, text="Cancel", command=cancel).grid(row=7, column=1)

    # Run the GUI
    root.mainloop()

    return settings
