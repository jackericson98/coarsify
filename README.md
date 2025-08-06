
<p align="center">
  <img width="300" height="300" alt="CoarsifyLogo2" src="https://github.com/user-attachments/assets/5e8cd4f7-21ed-4d39-83b7-65efc50ef78f" />
</p>

# Coarsify

**Coarsify** is a Python package for coarse-graining molecular structures into simplified, sphere-based representations. It supports a wide range of input file formats and multiple coarse-graining schemes, allowing researchers to reduce molecular complexity for visualization, analysis, or simulation workflows.

<p align="center">
  <img width="1263" height="397" alt="CoarsifyBanner-page001" src="https://github.com/user-attachments/assets/d21cef57-7d90-4ce3-94c9-096f8f119de4" />
</p>

---

## 🚀 Features

- **Multi-format Input Support**  
  Supports common molecular structure formats: `.pdb`, `.gro`, `.mol`, `.cif`, `.xyz`

- **Flexible Coarse-Graining Schemes**
  - **Average Distance**: Center of mass with radius from atomic dispersion
  - **Encapsulate Residues**: Minimal bounding sphere per residue
  - **Side Chain/Backbone Split**: Separate backbone and sidechain beads
  - **Martini Mapping**: MARTINI-compatible CG mappings

- **PyMOL Integration**  
  Automatically generate `.pml` scripts for visualization with correct radii and color coding

- **Parameter Control**
  - Mass-weighted or geometric centers
  - Optional thermal cushion
  - Hydrogen inclusion toggle
  - Residue splitting for detailed structure

- **Multiple Output Types**  
  Export coarse-grained structures in `.pdb`, along with text-based coordinate and radius files

- **Published on PyPI**  
  Install with `pip install coarsify` and integrate with your Python workflows easily

---

## 🛠️ Installation

### From PyPI (Recommended)

```bash
pip install coarsify
```

### From Source

```bash
git clone https://github.com/jackericson98/coarsify.git
cd coarsify
pip install .
```

### Requirements
- Python 3.7+
- `numpy >= 1.24.0`
- `pandas >= 2.0.0`

---

## 🧪 Usage

### GUI (Recommended)

Launch GUI from command line:

```bash
coarsify
```

Or if installed from source:

```bash
python -m coarsify
```

#### GUI Features
- File selector for input structures
- Dropdowns for all coarse-graining schemes and settings
- Output location selector
- Context-sensitive help

---

### CLI

For scripting or batch processing:

```bash
coarsify
```

Or from source:

```bash
python -m coarsify
```

#### CLI Features
- Interactive prompts for parameters
- Dialog-based file selection
- Compatible with headless environments

---

## ⚙️ Coarse-Graining Schemes

| Scheme                  | Description                                                                 |
|-------------------------|-----------------------------------------------------------------------------|
| **Average Distance**    | Center of mass and average atom distance radius                            |
| **Encapsulate Residues**| Minimal bounding sphere for all atoms in a residue                         |
| **Side Chain / Backbone**| Separate backbone and sidechain beads using average or encapsulation       |
| **Martini Mapping**     | Apply MARTINI force field mapping for simulations                          |

---

## 🔧 Parameters

- `thermal_cushion` – Adds buffer to radii (in Å) to simulate thermal motion
- `mass_weighted` – Use mass-weighted center instead of geometric center
- `include_hydrogens` – Include or exclude hydrogen atoms
- `split_residue` – Generate separate beads for backbone and sidechain atoms

---

## 📁 Output Files

Each run generates:

| File Name                 | Description                                      |
|---------------------------|--------------------------------------------------|
| `[name]_[scheme].pdb`     | Coarse-grained structure                         |
| `[name]_base.pdb`         | Original input structure                         |
| `set_atom_colors.pml`     | PyMOL script for radius/color settings           |
| `[name]_[scheme].txt`     | Text file with (x, y, z, radius) for each bead   |

---

## 🔬 PyMOL Visualization

1. Open the output `.pdb` file in PyMOL  
2. Run the `set_atom_colors.pml` script  
3. The script:
   - Sets sphere radii to match bead size
   - Colors beads by residue or type
   - Applies basic visual styling

---

## 📌 Example Workflow

```text
1. Launch GUI and load a structure
2. Select the "Encapsulate Residues" scheme
3. Set thermal cushion to 1.0 Å, exclude hydrogens
4. Click run to generate output files
5. Visualize in PyMOL using the provided script
```

---

## 🔗 Integration with Vorpy

Coarsify outputs are fully compatible with [Vorpy](https://github.com/your-username/vorpy), a tool for Voronoi diagram analysis. This allows downstream geometric or topological analysis of the coarse-grained structures.

---

## 🤝 Contributing

We welcome contributions for:
- New coarse-graining methods
- Additional input/output format support
- Speed improvements
- Documentation and tutorials

---

## 📄 License

MIT License – see the [LICENSE](LICENSE) file for details.

---

## 📚 Citation

```bibtex
@software{coarsify2024,
  author = {John Ericson},
  title = {Coarsify: A Python package for coarse graining molecular structures},
  year = {2024},
  publisher = {GitHub},
  url = {https://github.com/jackericson98/coarsify}
}
```

---

## 📬 Contact

For bug reports or feature requests, [create an issue](https://github.com/jackericson98/coarsify/issues) or contact [jackericson98@gmail.com](mailto:jackericson98@gmail.com)
