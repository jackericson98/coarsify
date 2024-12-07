class Residue:
    def __init__(self, atoms=None, name=None, sys=None, mol=None, sequence=None, seg_id=None, chain=None, element='Al'):
        """
        Residue Object for holding specific residue information
        :param atoms:
        :param name:
        :param sys:
        :param mol:
        :param sequence:
        :param seg_id:
        """
        self.atoms = atoms
        self.name = name
        self.sys = sys
        self.mol = mol
        self.seq = sequence
        self.id = seg_id
        self.chain = chain
        self.print_name = None
        self.element = element
        try:
            self.color = residue_colors[name][sys.color_scheme]
        except KeyError:
            self.color = '0xFF00FF'

    def add_atom(self, atom):
        self.atoms.append(atom)


residue_colors = {
    # Amino Acids
    "ALA": {"Shapely": "0x8CFF8C"},
    "GLY": {"Shapely": "0xFFFFFF"},
    "LEU": {"Shapely": "0x455E45"},
    "SER": {"Shapely": "0xFF7042"},
    "VAL": {"Shapely": "0xFF8CFF"},
    "THR": {"Shapely": "0xB84C00"},
    "LYS": {"Shapely": "0x4747B8"},
    "ASP": {"Shapely": "0xA00042"},
    "ILE": {"Shapely": "0x004C00"},
    "ASN": {"Shapely": "0xFF7C70"},
    "GLU": {"Shapely": "0x660000"},
    "PRO": {"Shapely": "0x525252"},
    "ARG": {"Shapely": "0x00007C"},
    "PHE": {"Shapely": "0x534C42"},
    "GLN": {"Shapely": "0xFF4C4C"},
    "TYR": {"Shapely": "0x8C704C"},
    "HIS": {"Shapely": "0x7070FF"},
    "CYS": {"Shapely": "0xFFFF70"},
    "MET": {"Shapely": "0xB8A042"},
    "TRP": {"Shapely": "0x4F4600"},
    "ASX": {"Shapely": "0xFF00FF"},  # Assuming ASX and GLX represent Asp/Asn and Glu/Gln ambiguous cases
    "GLX": {"Shapely": "0xFF00FF"},
    "PCA": {"Shapely": "0xFF00FF"},  # Rare in standard use, included for completeness
    "HYP": {"Shapely": "0xFF00FF"},  # Rare in standard use, included for completeness

    # Nucleic Acids
    **{_: {"Shapely": "0xA0A0FF"} for _ in {"DA", "A"}},
    **{_: {"Shapely": "0xFF8C4B"} for _ in {"DC", "C"}},
    **{_: {"Shapely": "0xFF7070"} for _ in {"DG", "G"}},
    **{_: {"Shapely": "0xA0FFA0"} for _ in {"DT", "T"}},
    **{_: {"Shapely": "0xB8B8B8"} for _ in {"DU", "U"}},

    # Other special cases
    "Backbone": {"Shapely": "0xB8B8B8"},
    "Special": {"Shapely": "0x5E005E"},
    "Default": {"Shapely": "0xFF00FF"},
    # Sol colors
    "SOL": {"Shapely": "0x00FFFF"},
    "HOH": {"Shapely": "0x00FFFF"},
    # Single atom colors
    "MG": {"Shapely": "0x00FA6D"},
    "C": {'Shapely': "0xC8C8C8"},
    "O": {'Shapely': "0xF00000"},
    "H": {'Shapely': "0xFFFFFF"},
    "N": {'Shapely': "0x8F8FFF"},
    "S": {'Shapely': "0xFFC832"},
    "P": {'Shapely': "0xFFA500"},
    "CL": {'Shapely': "0x00FF00"},
    "BR": {'Shapely': "0xA52A2A"},
    "ZN": {'Shapely': "0xA52A2A"},
    "NA": {'Shapely': "0x0000FF"},
    "FE": {'Shapely': "0xFFA500"},
    "CA": {'Shapely': "0x808090"}
}

