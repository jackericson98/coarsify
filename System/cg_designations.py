############################################# Amino Acid designations ##################################################

proteins = {
    'ALA': {'BB': {'atoms': {'N', 'HN', 'CA', 'C', 'O'}, 'type': 'polar', 'size': 4.1, 'designation': 'SP2'},
            'SC1': {'atoms': {'CB'}, 'type': 'nonpolar', 'size': 3.4, 'designation': 'TC3'}},
    'ARG': {'BB': {'atoms': {'N', 'HN', 'CA', 'C', 'O'}, 'type': 'polar', 'size': 4.7, 'designation': 'P2'},
            'SC1': {'atoms': {'CB', 'CG', 'CD'}, 'type': 'nonpolar', 'size': 4.1, 'designation': 'SC3'},
            'SC2': {'atoms': {'NE', 'HE', 'CZ', 'NH1', 'HH11', 'HH12', 'NH2', 'HH21', 'HH22'}, 'type': 'charged', 'size': 4.1, 'designation': 'SQ3p'}},
    'ASN': {'BB': {'atoms': {'N', 'HN', 'CA', 'C', 'O'}, 'type': 'polar', 'size': 4.7, 'designation': 'P2'},
            'SC1': {'atoms': {'CB', 'CG', 'OD1', 'ND2', 'HD21', 'HD22'}, 'type': 'polar', 'size': 4.1, 'designation': 'SP5'}},
    'ASP': {'BB': {'atoms': {'N', 'HN', 'CA', 'C', 'O'}, 'type': 'polar', 'size': 4.7, 'designation': 'P2'},
            'SC1': {'atoms': {'CB', 'CG', 'OD1', 'OD2'}, 'type': 'charged', 'size': 4.1, 'designation': 'SQ5n'}},
    'CYS': {'BB': {'atoms': {'N', 'HN', 'CA', 'C', 'O'}, 'type': 'polar', 'size': 4.7, 'designation': 'P2'},
            'SC1': {'atoms': {'CB', 'SG'}, 'type': 'nonpolar', 'size': 3.4, 'designation': 'TC6'}},
    'GLN': {'BB': {'atoms': {'N', 'HN', 'CA', 'C', 'O'}, 'type': 'polar', 'size': 4.7, 'designation': 'P2'},
            'SC1': {'atoms': {'CB', 'CG', 'CD', 'OE1', 'NE2', 'HE21', 'HE22', 'HE11', 'HE12'}, 'type': 'polar', 'size': 4.7, 'designation': 'P5'}},
    'GLU': {'BB': {'atoms': {'N', 'HN', 'CA', 'C', 'O'}, 'type': 'polar', 'size': 4.7, 'designation': 'P2'},
            'SC1': {'atoms': {'CB', 'CG', 'CD', 'OE1', 'OE2'}, 'type': 'charged', 'size': 4.7, 'designation': 'Q5n'}},
    'GLY': {'BB': {'atoms': {'N', 'HN', 'CA', 'C', 'O'}, 'type': 'polar', 'size': 4.1, 'designation': 'SP1'}},
    'HIS': {'BB': {'atoms': {'N', 'HN', 'CA', 'C', 'O'}, 'type': 'polar', 'size': 4.7, 'designation': 'P2'},
            'SC1': {'atoms': {'CB', 'CG'}, 'type': 'nonpolar', 'size': 3.4, 'designation': 'TC4'},
            'SC2': {'atoms': {'CD2', 'HD2', 'NE2', 'HE2'}, 'type': 'intermediate', 'size': 3.4, 'designation': 'TN6d'},
            'SC3': {'atoms': {'ND1', 'HD1', 'CE1', 'HE1'}, 'type': 'intermediate', 'size': 3.4, 'designation': 'TN5a'}},
    'ILE': {'BB': {'atoms': {'N', 'HN', 'CA', 'C', 'O'}, 'type': 'polar', 'size': 4.7, 'designation': 'P2'},
            'SC1': {'atoms': {'CB', 'CG2', 'CG1', 'CD'}, 'type': 'nonpolar', 'size': 4.1, 'designation': 'SC2'}},
    'LEU': {'BB': {'atoms': {'N', 'HN', 'CA', 'C', 'O'}, 'type': 'polar', 'size': 4.7, 'designation': 'P2'},
            'SC1': {'atoms': {'CB', 'CG', 'CD1', 'CD2'}, 'type': 'nonpolar', 'size': 4.1, 'designation': 'SC2'}},
    'LYS': {'BB': {'atoms': {'N', 'HN', 'CA', 'C', 'O'}, 'type': 'polar', 'size': 4.7, 'designation': 'P2'},
            'SC1': {'atoms': {'CB', 'CG', 'CD'}, 'type': 'nonpolar', 'size': 4.1, 'designation': 'SC3'},
            'SC2': {'atoms': {'CE', 'NZ', 'HZ1', 'HZ2', 'HZ3'}, 'type': 'charged', 'size': 4.1, 'designation': 'SQ4p'}},
    'MET': {'BB': {'atoms': {'N', 'HN', 'CA', 'C', 'O'}, 'type': 'polar', 'size': 4.7, 'designation': 'P2'},
            'SC1': {'atoms': {'CB', 'CG', 'SD', 'CE'}, 'type': 'nonpolar', 'size': 4.7, 'designation': 'C6'}},
    'PHE': {'BB': {'atoms': {'N', 'HN', 'CA', 'C', 'O'}, 'type': 'polar', 'size': 4.7, 'designation': 'P2'},
            'SC1': {'atoms': {'CB', 'HB1', 'HB2', 'CG', 'CD1', 'HD1'}, 'type': 'nonpolar', 'size': 4.1, 'designation': 'SC4'},
            'SC2': {'atoms': {'CE1', 'HE1'}, 'type': 'nonpolar', 'size': 3.4, 'designation': 'TC5'},
            'SC3': {'atoms': {'CZ', 'HZ', 'CD2', 'HD2', 'CE2', 'HE2'}, 'type': 'nonpolar', 'size': 3.4, 'designation': 'TC5'}},
    'PRO': {'BB': {'atoms': {'N', 'CA', 'C', 'O'}, 'type': 'polar', 'size': 4.1, 'designation': 'SP2a'},
            'SC1': {'atoms': {'CD', 'CB', 'CG'}, 'type': 'nonpolar', 'size': 4.1, 'designation': 'SC3'}},
    'SER': {'BB': {'atoms': {'N', 'HN', 'CA', 'C', 'O'}, 'type': 'polar', 'size': 4.7, 'designation': 'P2'},
            'SC1': {'atoms': {'CB', 'OG', 'HG1'}, 'type': 'polar', 'size': 3.4, 'designation': 'TP1'}},
    'THR': {'BB': {'atoms': {'N', 'HN', 'CA', 'C', 'O'}, 'type': 'polar', 'size': 4.7, 'designation': 'P2'},
            'SC1': {'atoms': {'CB', 'OG1', 'HG1', 'CG2'}, 'type': 'polar', 'size': 4.1, 'designation': 'SP1'}},
    'TRP': {'BB': {'atoms': {'N', 'HN', 'CA', 'C', 'O'}, 'type': 'polar', 'size': 4.7, 'designation': 'P2'},
            'SC1': {'atoms': {'CB', 'CG'}, 'type': 'nonpolar', 'size': 3.4, 'designation': 'TC4'},
            'SC2': {'atoms': {'CD1', 'HD1', 'NE1', 'HE1'}, 'type': 'intermediate', 'size': 3.4, 'designation': 'TN6d'},
            'SC3': {'atoms': {'CE2', 'CD2'}, 'type': 'nonpolar', 'size': 3.4, 'designation': 'TC5'},
            'SC4': {'atoms': {'CE3', 'HE3', 'CZ3', 'HZ3'}, 'type': 'nonpolar', 'size': 3.4, 'designation': 'TC5'},
            'SC5': {'atoms': {'CZ2', 'HZ2', 'CH2', 'HH2'}, 'type': 'nonpolar', 'size': 3.4, 'designation': 'TC5'}},
    'TYR': {'BB': {'atoms': {'N', 'HN', 'CA', 'C', 'O'}, 'type': 'polar', 'size': 4.7, 'designation': 'P2'},
            'SC1': {'atoms': {'CB', 'CG'}, 'type': 'intermediate', 'size': 3.4, 'designation': 'TC4'},
            'SC2': {'atoms': {'CD1', 'HD1', 'CE1', 'HE1'}, 'type': 'nonpolar', 'size': 3.4, 'designation': 'TC5'},
            'SC3': {'atoms': {'CD2', 'HD2', 'CE2', 'HE2'}, 'type': 'nonpolar', 'size': 3.4, 'designation': 'TC5'},
            'SC4': {'atoms': {'CZ', 'OH', 'HH'}, 'type': 'intermediate', 'size': 3.4, 'designation': 'TN6'}},
    'VAL': {'BB': {'atoms': {'N', 'HN', 'CA', 'C', 'O'}, 'type': 'polar', 'size': 4.1, 'designation': 'SP2'},
            'SC1': {'atoms': {'CB', 'CG1', 'CG2'}, 'type': 'nonpolar', 'size': 4.1, 'designation': 'SC3'}}
}

################################################## Nucleic Acid Base designations ######################################

nucleobases = {
    'ADEN': {'N1': {'atoms': {}, 'type': 'intermediate', 'size': 4.1, 'designation': 'SN1'},
             'N2': {'atoms': {}, 'type': 'intermediate', 'size': 3.4, 'designation': 'TN3A'},
             'N3': {'atoms': {}, 'type': 'intermediate', 'size': 3.4, 'designation': 'TN5a'},
             'N4': {'atoms': {}, 'type': 'polar', 'size': 3.4, 'designation': 'TP1d'},
             'N5': {'atoms': {}, 'type': 'intermediate', 'size': 3.4, 'designation': 'TN3a'},
             'N6': {'atoms': {}, 'type': 'nonpolar', 'size': 3.4, 'designation': 'TC6'}},
    'CYTO': {'N1': {'atoms': {}, 'type': 'intermediate', 'size': 4.1, 'designation': 'SN1'},
             'N2': {'atoms': {}, 'type': 'polar', 'size': 3.4, 'designation': 'TP1a'},
             'N3': {'atoms': {}, 'type': 'intermediate', 'size': 3.4, 'designation': 'TN5a'},
             'N4': {'atoms': {}, 'type': 'polar', 'size': 3.4, 'designation': 'TP1d'}},
    'GUAN': {'N1': {'atoms': {}, 'type': 'intermediate', 'size': 4.1, 'designation': 'SN1'},
             'N2': {'atoms': {}, 'type': 'polar', 'size': 3.4, 'designation': 'TP1d'},
             'N3': {'atoms': {}, 'type': 'polar', 'size': 3.4, 'designation': 'TP1d'},
             'N4': {'atoms': {}, 'type': 'polar', 'size': 3.4, 'designation': 'TP1a'},
             'N5': {'atoms': {}, 'type': 'intermediate', 'size': 3.4, 'designation': 'TN5a'},
             'N6': {'atoms': {}, 'type': 'intermediate', 'size': 3.4, 'designation': 'TN3a'}},
    'THYM': {'N1': {'atoms': {}, 'type': 'intermediate', 'size': 4.1, 'designation': 'SN1'},
             'N2': {'atoms': {}, 'type': 'polar', 'size': 3.4, 'designation': 'TP1a'},
             'N3': {'atoms': {}, 'type': 'intermediate', 'size': 3.4, 'designation': 'TN5d'},
             'N4': {'atoms': {}, 'type': 'polar', 'size': 3.4, 'designation': 'TP1a'},
             'N5': {'atoms': {}, 'type': 'nonpolar', 'size': 3.4, 'designation': 'TC2'}},
    'URAC': {'N1': {'atoms': {}, 'type': 'intermediate', 'size': 4.1, 'designation': 'SN1'},
             'N2': {'atoms': {}, 'type': 'polar', 'size': 3.4, 'designation': 'TP1a'},
             'N3': {'atoms': {}, 'type': 'intermediate', 'size': 3.4, 'designation': 'TN5d'},
             'N4': {'atoms': {}, 'type': 'polar', 'size': 3.4, 'designation': 'TP1a'},
             'N5': {'atoms': {}, 'type': 'nonpolar ', 'size': 3.4, 'designation': 'TC5'}}
}

##################################################### Ions and Solvent designations ####################################

ions = {'ION': {
    'NA': {'atoms': {'NA'}, 'type': 'charged', 'size': 3.4, 'charge': 1.0, 'designation': 'TQ5', 'mass': None},
    'CL': {'atoms': {'CL'}, 'type': 'charged', 'size': 3.4, 'charge': -1.0, 'designation': 'TQ5', 'mass': 35.453},
    'BR': {'atoms': {'BR'}, 'type': 'charged', 'size': 4.1, 'charge': -1.0, 'designation': 'SQ4', 'mass': 79.90},
    'IOD': {'atoms': {'ID'}, 'type': 'charged', 'size': 4.1, 'charge': -1.0, 'designation': 'SQ2', 'mass': 79.90},
    'TMA': {'atoms': {'TMA'}, 'type': 'charged', 'size': 4.7, 'charge': 1.0, 'designation': 'Q2', 'mass': 74.14},
    'ACE': {'atoms': {'CL'}, 'type': 'charged', 'size': 4.1, 'charge': -1.0, 'designation': 'SQ5n', 'mass': 59.044},
    'CA': {'atoms': {'TMA'}, 'type': 'divalent', 'size': 4.1, 'charge': 2.0, 'designation': 'SD', 'mass': 40.078}}
}

solvents = {
    'W': {'W': {'atoms': {}, 'type': 'polar', 'size': 4.1}}
}