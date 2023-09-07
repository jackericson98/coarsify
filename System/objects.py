

def make_atom(system=None, location=None, radius=None, index='', name='', residue='', chain='', res_seq="", seg_id="",
              element="", chn=None, res=None):
    atom = {
        # System groups
        'sys': system,           # System       :   Main system object
        'res': res,              # Residue      :   Residue object of which the atom is a part
        'chn': chn,              # Chain        :   Chain object of which the atom is a part

        'loc': location,         # Location     :   Set the location of the center of the sphere
        'rad': radius,           # Radius       :   Set the radius for the sphere object. Default is 1

        # Calculated Traits
        'vol': 0,                # Cell Volume  :   Volume of the voronoi cell for the atom
        'sa': 0,                 # Surface Area :   Surface area of the atom's cell
        'curv': 0,
        'box': [],               # Box          :   The grid location of the atom

        # Network objects
        'averts': [],             # Vertices     :   List of Vertex type objects
        'asurfs': [],             # Surfaces     :   List of Surface type objects
        'aedges': [],             # Edges        :   List of Edge type objects

        # Input traits
        'num': index,            # Number       :   The index from the initial atom file
        'name': name,            # Name         :   Name retrieved from pdb file
        'chain': chain,          # Chain        :   Molecule chain the atom is a part of
        'residue': residue,      # Residue      :   String describing the residue type that the atom is a part of
        'res_seq': res_seq,      # Sequence     :   Sequence of the residue that the atom is a part of
        'seg_id': seg_id,        # Segment ID   :   Segment identifier for the atom
        'element': element,      # Symbol       :   Element of the atom
    }
    return atom


class Residue:
    def __init__(self, atoms=None, name=None, sys=None, mol=None, sequence=None, seg_id=None, chain=None):
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
        self.loc = None
        self.rad = None
        self.elem_col = None

    def add_atom(self, atom):
        self.atoms.append(atom)


class Chain:
    def __init__(self, sys=None, atoms=None, residues=None, name=None):
        self.name = name
        self.atoms = atoms
        self.residues = residues
        self.vol = 0
        self.sa = 0

    def add_atom(self, atom):
        self.atoms.append(atom)


class Sol(Chain):
    def __init__(self, sys=None, atoms=None, residues=None, name="H2O"):
        super().__init__()
        self.atoms = atoms
        self.residues = residues
        self.name = name
        self.vol = 0
        self.sa = 0
