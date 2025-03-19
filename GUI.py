import tkinter as tk
from tkinter import ttk
from tkinter import filedialog
import tkinter.font as tkfont

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
        return settings

    def cancel():
        root.destroy()

    def show_help():
        help_window = tk.Toplevel(root)
        help_window.title("Coarsify Help")
        help_window.configure(bg='#f0f0f0')
        help_window.minsize(400, 500)
        
        # Make help window modal
        help_window.transient(root)
        help_window.grab_set()
        
        # Create main frame for help content
        help_frame = ttk.Frame(help_window, padding="20")
        help_frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        help_window.grid_columnconfigure(0, weight=1)
        help_window.grid_rowconfigure(0, weight=1)
        
        # Title
        help_title = tk.Label(help_frame, text="Help Guide", 
                            font=tkfont.Font(family="Arial", size=24, weight="bold"),
                            fg="#2196F3", bg='#f0f0f0')
        help_title.grid(row=0, column=0, pady=(0, 20))
        
        # Help content
        help_text = {
            "File Selection": "Choose the input PDB file that you want to coarse-grain. This should be a properly formatted protein structure file.",
            
            "CG Method": """Choose the coarse-graining method to use:
• Encapsulate: Groups atoms based on residue structure
• Average Distance: Uses distance-based clustering
• Martini: Applies the MARTINI force field mapping
• All Schemes: Applies all methods for comparison""",
            
            "Thermal Cushion": "Additional distance (in Angstroms) added to the coarse-graining radius to account for thermal motion. Default is 0.0Å.",
            
            "Split Residue": "When enabled, separates backbone and sidechain atoms into different coarse-grained beads.",
            
            "Include Hydrogens": "When enabled, includes hydrogen atoms in the coarse-graining process. Disable to ignore hydrogens.",
            
            "Mass Weighted": "When enabled, uses mass-weighted averaging for bead positions. Disable for geometric center calculation.",
            
            "Output Folder": "Select the directory where the coarse-grained structure and analysis files will be saved."
        }
        
        # Create text widget for scrollable help content
        help_text_widget = tk.Text(help_frame, wrap=tk.WORD, width=50, height=20,
                                 font=tkfont.Font(family="Arial", size=10),
                                 bg='white', relief="flat")
        help_text_widget.grid(row=1, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        
        # Add scrollbar
        scrollbar = ttk.Scrollbar(help_frame, orient="vertical", command=help_text_widget.yview)
        scrollbar.grid(row=1, column=1, sticky=(tk.N, tk.S))
        help_text_widget.configure(yscrollcommand=scrollbar.set)
        
        # Insert help content with formatting
        help_text_widget.tag_configure("heading", font=tkfont.Font(family="Arial", size=12, weight="bold"))
        help_text_widget.tag_configure("content", font=tkfont.Font(family="Arial", size=10))
        
        for topic, description in help_text.items():
            help_text_widget.insert(tk.END, f"{topic}\n", "heading")
            help_text_widget.insert(tk.END, f"{description}\n\n", "content")
        
        help_text_widget.configure(state='disabled')  # Make text read-only
        
        # Close button
        ttk.Button(help_frame, text="Close", 
                   command=help_window.destroy).grid(row=2, column=0, pady=20)
        
        # Center help window relative to main window
        help_window.update_idletasks()
        width = help_window.winfo_width()
        height = help_window.winfo_height()
        parent_x = root.winfo_x()
        parent_y = root.winfo_y()
        parent_width = root.winfo_width()
        parent_height = root.winfo_height()
        x = parent_x + (parent_width // 2) - (width // 2)
        y = parent_y + (parent_height // 2) - (height // 2)
        help_window.geometry(f'+{x}+{y}')

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
    root.title("Coarsify")
    root.configure(bg='#f0f0f0')
    
    # Set minimum window size
    root.minsize(340, 600)
    
    # Configure style
    style = ttk.Style()
    style.configure('TButton', padding=6, relief="flat", background="#2196F3")
    style.configure('TEntry', padding=6)
    style.configure('TCombobox', padding=6)
    
    # Create main frame
    main_frame = ttk.Frame(root, padding="20")
    main_frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
    root.grid_columnconfigure(0, weight=1)
    root.grid_rowconfigure(0, weight=1)

    # Title Frame with Help Button
    title_frame = ttk.Frame(main_frame)
    title_frame.grid(row=0, column=0, columnspan=3, pady=(0, 20))
    title_frame.grid_columnconfigure(1, weight=1)  # This makes the title center between help button
    
    # Create title with custom font
    title_font = tkfont.Font(family="Arial", size=36, weight="bold")
    title_label = tk.Label(title_frame, text="Coarsify", font=title_font, fg="#2196F3", bg='#f0f0f0')
    title_label.grid(row=0, column=1, pady=20)
    
    # Help button
    help_button = ttk.Button(title_frame, text="?", width=3, command=show_help)
    help_button.grid(row=0, column=2, padx=(10, 0), sticky='ne')

    # Content Frame
    content_frame = ttk.Frame(main_frame, padding="10")
    content_frame.grid(row=1, column=0, columnspan=3, sticky=(tk.W, tk.E))
    
    # Initialize variables
    input_file_var = tk.StringVar(value='Choose File')
    cg_method_var = tk.StringVar(value='Encapsulate')
    thermal_cushion_var = tk.StringVar(value='0.0')
    sc_bb_var = tk.BooleanVar(value=False)
    include_h_var = tk.BooleanVar(value=True)
    mass_weighted_var = tk.BooleanVar(value=True)
    output_folder_var = tk.StringVar(value='Choose Output Folder')

    # File Selection Frame
    file_frame = ttk.LabelFrame(content_frame, text="File Selection", padding="10")
    file_frame.grid(row=0, column=0, columnspan=3, sticky=(tk.W, tk.E), pady=(0, 10))
    
    ttk.Label(file_frame, textvariable=input_file_var).grid(row=0, column=0, columnspan=2, sticky=tk.W, pady=5)
    ttk.Button(file_frame, text='Browse', command=choose_input_file).grid(row=0, column=2, padx=5)

    # Settings Frame
    settings_frame = ttk.LabelFrame(content_frame, text="Settings", padding="10")
    settings_frame.grid(row=1, column=0, columnspan=3, sticky=(tk.W, tk.E), pady=(0, 10))

    # Method Selection
    ttk.Label(settings_frame, text='CG Method:').grid(row=0, column=0, sticky=tk.W, pady=5)
    method_menu = ttk.Combobox(settings_frame, textvariable=cg_method_var, 
                              values=['Encapsulate', 'Average Distance', 'Martini', 'All Schemes'])
    method_menu.grid(row=0, column=1, columnspan=2, sticky=(tk.W, tk.E), pady=5, padx=5)
    method_menu.current(0)

    # Thermal Cushion
    ttk.Label(settings_frame, text='Thermal Cushion:').grid(row=1, column=0, sticky=tk.W, pady=5)
    ttk.Entry(settings_frame, textvariable=thermal_cushion_var).grid(row=1, column=1, sticky=(tk.W, tk.E), pady=5)
    ttk.Label(settings_frame, text='\u212B').grid(row=1, column=2)

    # Checkboxes Frame
    checkbox_frame = ttk.Frame(settings_frame)
    checkbox_frame.grid(row=2, column=0, columnspan=3, pady=10)

    ttk.Checkbutton(checkbox_frame, text='Split Residue?', variable=sc_bb_var).pack(pady=2)
    ttk.Checkbutton(checkbox_frame, text='Include Hydrogens?', variable=include_h_var).pack(pady=2)
    ttk.Checkbutton(checkbox_frame, text='Mass Weighted', variable=mass_weighted_var).pack(pady=2)

    # Output Frame
    output_frame = ttk.LabelFrame(content_frame, text="Output", padding="10")
    output_frame.grid(row=2, column=0, columnspan=3, sticky=(tk.W, tk.E), pady=(0, 10))

    ttk.Label(output_frame, textvariable=output_folder_var).grid(row=0, column=0, columnspan=2, sticky=tk.W, pady=5)
    ttk.Button(output_frame, text='Browse', command=choose_output_folder).grid(row=0, column=2, padx=5)

    # Buttons Frame
    button_frame = ttk.Frame(content_frame)
    button_frame.grid(row=3, column=0, columnspan=3, pady=20)

    ttk.Button(button_frame, text="Apply", command=apply_values, style='TButton').pack(side=tk.RIGHT, padx=5)
    ttk.Button(button_frame, text="Cancel", command=cancel, style='TButton').pack(side=tk.RIGHT, padx=5)

    # Center the window on the screen
    root.update_idletasks()
    width = root.winfo_width()
    height = root.winfo_height()
    x = (root.winfo_screenwidth() // 2) - (width // 2)
    y = (root.winfo_screenheight() // 2) - (height // 2)
    root.geometry(f'{width}x{height}+{x}+{y}')

    root.mainloop()
    return settings


if __name__ == '__main__':
    print(settings_gui())
