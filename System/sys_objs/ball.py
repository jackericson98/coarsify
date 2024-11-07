class Ball:
    def __init__(self, loc=None, rad=None, element=None, residues=None, atoms=None, name=None, chain=None, seq=None,
                 index=None, residue_subsection=None, mass=None):
        self.loc = loc
        self.rad = rad
        self.element = element
        self.residues = residues
        self.atoms = atoms
        self.name = name
        self.chain = chain
        self.seq = seq
        self.index = index
        self.sub_section = residue_subsection
        self.mass = mass