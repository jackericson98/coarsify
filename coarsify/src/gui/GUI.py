#!/usr/bin/env python3
"""
GUI entry point for the coarsify graphical user interface.
"""

import sys
import tkinter as tk
from tkinter import ttk
from tkinter import filedialog
import tkinter.font as tkfont

settings = False

"""
Thermal cushion, methodology, sc_bb, average distance, include H, input file, output file
"""


class RoundedButton(tk.Canvas):
  def __init__(self, parent, border_radius, padding, color, text='', command=None):
    tk.Canvas.__init__(self, parent, borderwidth=0,
                       relief="raised", highlightthickness=0, bg=parent["bg"])
    self.command = command
    font_size = 10
    self.font = tkfont.Font(size=font_size, family='Helvetica')
    self.id = None
    height = font_size + (1 * padding)
    width = self.font.measure(text)+(1*padding)

    width = width if width >= 80 else 80

    if border_radius > 0.5*width:
      print("Error: border_radius is greater than width.")
      return None

    if border_radius > 0.5*height:
      print("Error: border_radius is greater than height.")
      return None

    rad = 2*border_radius

    def shape():
      self.create_arc((0, rad, rad, 0),
                      start=90, extent=90, fill=color, outline=color)
      self.create_arc((width-rad, 0, width,
                        rad), start=0, extent=90, fill=color, outline=color)
      self.create_arc((width, height-rad, width-rad,
                        height), start=270, extent=90, fill=color, outline=color)
      self.create_arc((0, height-rad, rad, height), start=180, extent=90, fill=color, outline=color)
      return self.create_polygon((0, height-border_radius, 0, border_radius, border_radius, 0, width-border_radius, 0, width,
                           border_radius, width, height-border_radius, width-border_radius, height, border_radius, height),
                                 fill=color, outline=color)

    id = shape()
    (x0, y0, x1, y1) = self.bbox("all")
    width = (x1-x0)
    height = (y1-y0)
    self.configure(width=width, height=height)
    self.create_text(width/2, height/2,text=text, fill='black', font= self.font)
    self.bind("<ButtonPress-1>", self._on_press)
    self.bind("<ButtonRelease-1>", self._on_release)

  def _on_press(self, event):
      self.configure(relief="sunken")

  def _on_release(self, event):
      self.configure(relief="raised")
      if self.command is not None:
          self.command()


def settings_gui():
    # Function to collect values and print dictionary
    global settings

    def truncate_path(path, num_chars=15):
        """Truncates a path to show only first and last num_chars characters with ... in between"""
        if len(path) <= num_chars * 2 + 3:  # If path is short enough, show full path
            return path
        return f"{path[:num_chars]}...{path[-num_chars:]}"

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
            "Coarsify": "Coarsify is a tool for coarse-graining protein structures. It is designed to be used in conjunction with the molecular visualization tool, PyMOL. ",
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
        if input_file:  # Only update if a file was selected
            input_file_display.set(truncate_path(input_file))
            input_file_var.set(input_file)  # Store full path
        choose_input_file_window.destroy()

    def choose_output_folder():
        choose_output_folder_window = tk.Tk()
        choose_output_folder_window.withdraw()
        choose_output_folder_window.wm_attributes('-topmost', 1)
        output_folder = filedialog.askdirectory()
        if output_folder:  # Only update if a folder was selected
            output_folder_display.set(truncate_path(output_folder))
            output_folder_var.set(output_folder)  # Store full path
        choose_output_folder_window.destroy()

    # Main window
    root = tk.Tk()
    root.wm_attributes('-topmost', 1)
    root.title("Coarsify")
    root.configure(bg='#f0f0f0')
    
    # Set minimum window size
    root.minsize(340, 600)
    
    # Configure style
    style = ttk.Style()
    style.configure('TButton', padding=6, relief="flat", background="#2196F3")
    style.configure('TEntry', padding=6)
    style.configure('TCombobox', padding=6)
    # Configure help button style
    style.configure('Help.TButton', 
                   padding=2,
                   relief="flat",
                   background="#f0f0f0",
                   borderwidth=1,
                   font=('Arial', 8))
    
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
    title_label.grid(row=0, column=1, pady=(20, 0))
    
    # Create subtitle
    subtitle_font = tkfont.Font(family="Arial", size=14)
    subtitle_label = tk.Label(title_frame, text="A molecular coarse graining tool", 
                            font=subtitle_font, fg="#666666", bg='#f0f0f0')
    subtitle_label.grid(row=1, column=1, pady=(0, 20))

    # Content Frame
    content_frame = ttk.Frame(main_frame, padding="10")
    content_frame.grid(row=1, column=0, columnspan=3, sticky=(tk.W, tk.E))
    
    # Initialize variables
    input_file_var = tk.StringVar(value='')  # Stores full path
    input_file_display = tk.StringVar(value='Choose File')  # Displays truncated path
    cg_method_var = tk.StringVar(value='Encapsulate')
    thermal_cushion_var = tk.StringVar(value='0.0')
    sc_bb_var = tk.BooleanVar(value=False)
    include_h_var = tk.BooleanVar(value=True)
    mass_weighted_var = tk.BooleanVar(value=True)
    output_folder_var = tk.StringVar(value='')  # Stores full path
    output_folder_display = tk.StringVar(value='Choose Output Folder')  # Displays truncated path

    # File Selection Frame
    file_frame = ttk.LabelFrame(content_frame, text="File Selection", padding="10")
    file_frame.grid(row=0, column=0, columnspan=3, sticky=(tk.W, tk.E), pady=(0, 10))
    
    # Configure column weights to push button to the right
    file_frame.grid_columnconfigure(0, weight=1)  # Make the label column expand
    file_frame.grid_columnconfigure(1, weight=0)  # Keep button column fixed
    
    ttk.Label(file_frame, textvariable=input_file_display).grid(row=0, column=0, sticky=tk.W, pady=5)
    ttk.Button(file_frame, text='Browse', command=choose_input_file).grid(row=0, column=1, padx=5, sticky=tk.E)

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

    # Configure column weights to push button to the right
    output_frame.grid_columnconfigure(0, weight=1)  # Make the label column expand
    output_frame.grid_columnconfigure(1, weight=0)  # Keep button column fixed
    
    ttk.Label(output_frame, textvariable=output_folder_display).grid(row=0, column=0, sticky=tk.W, pady=5)
    ttk.Button(output_frame, text='Browse', command=choose_output_folder).grid(row=0, column=1, padx=5, sticky=tk.E)

    # Buttons Frame
    button_frame = ttk.Frame(content_frame)
    button_frame.grid(row=3, column=0, columnspan=3, pady=20)
    button_frame.grid_columnconfigure(1, weight=1)  # Center column for main buttons
    
    # Left frame for help button
    left_button_frame = ttk.Frame(button_frame)
    left_button_frame.grid(row=0, column=0, padx=5, sticky=tk.W)
    
    # Center frame for main buttons
    center_button_frame = ttk.Frame(button_frame)
    center_button_frame.grid(row=0, column=1)

    # Main buttons in center
    ttk.Button(center_button_frame, text="Apply", command=apply_values, style='TButton').pack(side=tk.LEFT, padx=5)
    ttk.Button(center_button_frame, text="Cancel", command=cancel, style='TButton').pack(side=tk.LEFT, padx=5)

    # Help button (small and circular)
    help_button = tk.Button(left_button_frame, 
                                text="?",
                                width=2,
                                height=1,
                                font=('Arial', 10),
                                fg='#666666',
                                bg='#f0f0f0',
                                relief='solid',
                                borderwidth=1,
                                command=show_help)
    help_button['border'] = "0"
    help_button.pack(side=tk.RIGHT)
    
    # Make the help button circular
    help_button.bind('<Configure>', lambda e: help_button.configure(width=2 if e.width > e.height else 1))
    
    # Hover effects
    def on_enter(e):
        help_button['bg'] = '#e6e6e6'
    
    def on_leave(e):
        help_button['bg'] = '#f0f0f0'
    
    help_button.bind('<Enter>', on_enter)
    help_button.bind('<Leave>', on_leave)

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
