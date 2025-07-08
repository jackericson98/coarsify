# Coarsify

A Python tool for coarse-graining molecular structures from various file formats (.pdb, .gro, .mol, .cif, .xyz) into simplified representations using different coarse-graining schemes. Designed for researchers working with molecular dynamics simulations who need to approximate molecules as fewer spheres or perform structural analysis.

## Features

- **Multiple Input Formats**: Supports PDB, GRO, MOL, CIF, and XYZ file formats
- **Various Coarse-Graining Schemes**:
  - **Average Distance**: Balls located at center of mass of residues with radius based on average distance of constituent atoms
  - **Encapsulate Residues**: Balls that minimally encapsulate atoms in a residue
  - **Side Chain/Backbone Split**: Separate beads for backbone and sidechain atoms
  - **Martini**: Pre-processed PDB files using CG-MARTINI force field mapping
- **Mass-Weighted Options**: Choose between mass-weighted or geometric center calculations
- **Thermal Cushion**: Add additional radius to account for thermal motion
- **Hydrogen Handling**: Option to include or exclude hydrogen atoms
- **PyMOL Integration**: Automatic generation of PyMOL scripts for visualization
- **Multiple Output Formats**: PDB files with coarse-grained structures and text files with coordinates and radii

## Installation

### Prerequisites
- Python 3.7 or higher
- Required Python packages (see `requirements.txt`):
  - numpy ~= 1.24.4
  - pandas ~= 2.0.3

### Setup
1. Clone or download this repository
2. Install dependencies:
   ```bash
   pip install -r requirements.txt
   ```

## Usage

### Graphical User Interface (Recommended)

The GUI provides an intuitive interface for all coarse-graining operations:

```bash
python GUI.py
```

**GUI Features:**
- File selection dialog for input structures
- Dropdown menus for all parameters
- Real-time parameter validation
- Help system with detailed explanations
- Output folder selection

### Command Line Interface

For batch processing or scripting:

```bash
python coarsify.py
```

**CLI Options:**
- Interactive prompts for all parameters
- File dialog for input selection
- Multiple output location options

## Coarse-Graining Schemes

### 1. Average Distance
- **Description**: Creates balls at the center of mass of residues
- **Radius**: Average distance of constituent atoms from the center
- **Options**: Mass-weighted or geometric center calculation
- **Best for**: General structural analysis and visualization

### 2. Encapsulate Residues
- **Description**: Creates minimal spheres that contain all atoms in a residue
- **Radius**: Minimum radius to encapsulate all residue atoms
- **Best for**: Preserving molecular volume and shape

### 3. Side Chain/Backbone Split
- **Description**: Separates backbone and sidechain atoms into different beads
- **Variants**: 
  - Average Distance with split
  - Encapsulate with split
- **Best for**: Detailed analysis of protein structure and dynamics

### 4. Martini
- **Description**: Applies CG-MARTINI force field mapping
- **Best for**: Molecular dynamics simulations with MARTINI force field

## Parameters

- **Thermal Cushion**: Additional radius (in Angstroms) added to account for thermal motion
- **Mass Weighted**: Use mass-weighted averaging for bead positions (vs. geometric center)
- **Include Hydrogens**: Include hydrogen atoms in coarse-graining calculations
- **Split Residue**: Separate backbone and sidechain into different beads

## Output Files

For each coarse-graining operation, the following files are generated:

1. **`[name]_[scheme].pdb`**: Coarse-grained structure in PDB format
2. **`[name]_base.pdb`**: Copy of the original input structure
3. **`set_atom_colors.pml`**: PyMOL script for sphere visualization
4. **`[name]_[scheme].txt`**: Text file with coordinates and radii of coarse-grained beads

## PyMOL Visualization

The tool automatically generates PyMOL scripts for easy visualization:

1. Load the coarse-grained PDB file in PyMOL
2. Run the generated `set_atom_colors.pml` script
3. The script will:
   - Set appropriate sphere radii for each bead
   - Apply color schemes based on residue types
   - Configure visualization settings

## Example Workflow

1. **Load Structure**: Select your PDB file through the GUI
2. **Choose Scheme**: Select "Average Distance" for general analysis
3. **Set Parameters**: 
   - Thermal cushion: 1.0 Å
   - Mass weighted: Yes
   - Include hydrogens: No
4. **Run Coarse-Graining**: Execute the process
5. **Visualize**: Open the output PDB in PyMOL and run the generated script

## Integration with Vorpy

This tool was originally developed to work in conjunction with [Vorpy](https://github.com/your-username/vorpy) for Voronoi diagram analysis of molecular structures. The coarse-grained representations can be used as input for Voronoi tessellation analysis.

## File Structure

```
coarsify/
├── GUI.py                 # Main GUI application
├── coarsify.py           # CLI interface
├── System/               # Core system modules
│   ├── system.py         # Main system class
│   ├── schemes/          # Coarse-graining algorithms
│   ├── sys_funcs/        # Utility functions
│   └── sys_objs/         # Object definitions
├── Data/                 # Output data directory
│   ├── test_data/        # Test structures
│   └── user_data/        # User-generated outputs
└── requirements.txt      # Python dependencies
```

## Contributing

This tool is actively developed for molecular dynamics research. Contributions are welcome, particularly for:
- Additional coarse-graining schemes
- Support for new file formats
- Performance optimizations
- Documentation improvements

## License

[Add your license information here]

## Citation

If you use this tool in your research, please cite:
[Add citation information when available]

## Contact

For questions, bug reports, or feature requests, please [create an issue](https://github.com/your-username/coarsify/issues) or contact [your-email@domain.com].
