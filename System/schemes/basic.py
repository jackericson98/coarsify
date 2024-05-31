from System.sys_funcs.calcs import calc_com, calc_dist
from System.sys_objs.ball import Ball


def make_ball(atoms, scheme='Average Distance', mass_weighted=True, include_h=True, therm_cush=0.0):
    # If we are including h
    if not include_h:
        atoms = [_ for _ in atoms if _['element'].lower() != 'h']

    # Calculate the center of mass for the atoms in a residue
    loc = calc_com([_['loc'] for _ in atoms], [_['mass'] for _ in atoms] if mass_weighted else None)

    # Choose the scheme for coarse graining the residues
    if scheme == 'Average Distance' and mass_weighted:
        rad = (sum([(calc_dist(loc, a['loc']) + a['rad']) * a['mass'] for a in atoms]) / sum([a['mass'] for a in atoms])
               + therm_cush)
    elif scheme == 'Average Distance':
        rad = sum([calc_dist(loc, a['loc']) + a['rad'] for a in atoms]) / len(atoms) + therm_cush
    else:
        rad = max([calc_dist(loc, a['loc']) + a['rad'] for a in atoms]) + therm_cush

    return loc, rad


def coarsify(sys):
    """
    Main coarsify function. Calculates radii and location for residues
    """
    # Check to see if the system has balls yet or not
    if sys.balls is None:
        sys.balls = []
    # Loop through the residues in the system
    for res in sys.residues:
        # If the residue is water we don't need to worry about backbone and side chain
        if res.name in sys.nucleics and sys.sc_bb:
            sugs, phos, nbase, hs = [], [], [], []
            # Get the back bone atoms and the side chain atoms
            for atom in res.atoms:
                if atom['element'].lower() == 'h':
                    hs.append(atom)
                elif atom['name'] in sys.nucleic_pphte:
                    phos.append(atom)
                elif atom['name'] in sys.nucleic_nbase:
                    sugs.append(atom)
                elif atom['name'] in sys.nucleic_sugr:
                    nbase.append(atom)
                elif atom['name'] in sys.nucleic_ignores:
                    continue
                else:
                    sc_bb_add = input("{} atom type not found in nucleic list (Residue {} # {})! Choose one of the "
                                      "following - 1. Phosphate, 2. Ribose, 3. Base, 4. Ignore \'{}\' 5. Skip\n>>> "
                                      .format(atom['name'], res.name, res.seq, atom['name']))
                    if sc_bb_add == '1':
                        sys.nucleic_pphte.append(atom['name'])
                        phos.append(atom)
                        print("{} will be added to nucleic phosphate grups from here on out".format(atom['name']))
                    elif sc_bb_add == '2':
                        sys.nucleic_nbase.append(atom['name'])
                        sugs.append(atom)
                        print("{} will be added to nucleic ribose groups from here on out".format(atom['name']))
                    elif sc_bb_add == 3:
                        sys.nucleic_sugr.append(atom['name'])
                        nbase.append(atom)
                        print("{} will be added to nucleic bases from here on out".format(atom['name']))
                    elif sc_bb_add == 4:
                        sys.nucleic_ignores.append(atom['name'])
                        print("{} will be ignored from here on out".format(atom['name']))
                    else:
                        continue

            # Sort the unknown hydrogen into their different groups
            for atom in hs:
                atom_dists = [(calc_dist(atom['loc'], _['loc']), _) for _ in res.atoms if _['element'].lower() != 'h']
                near_atoms = sorted(atom_dists, key=lambda x: x[0])
                close_atom = near_atoms[0][1]
                if close_atom['name'] in sys.nucleic_pphte:
                    phos.append(atom)
                elif close_atom['name'] in sys.nucleic_nbase:
                    sugs.append(atom)
                elif close_atom['name'] in sys.nucleic_sugr:
                    nbase.append(atom)
                elif close_atom['name'] in sys.nucleic_ignores:
                    continue

            # Create the ball objects
            if len(phos) > 0:
                ph_loc, ph_rad = make_ball(phos, sys.scheme, sys.mass_weighted, sys.include_h, sys.therm_cush)
                sys.balls.append(Ball(loc=ph_loc, rad=ph_rad, element=res.elem_col, residues=[res], atoms=phos,
                                      name=res.name, chain=res.chain, seq=res.seq, residue_subsection='phosphate'))
            sug_loc, sug_rad = make_ball(sugs, sys.scheme, sys.mass_weighted, sys.include_h, sys.therm_cush)
            nbas_loc, nbas_rad = make_ball(nbase, sys.scheme, sys.mass_weighted, sys.include_h, sys.therm_cush)
            sys.balls += [
                Ball(loc=sug_loc, rad=sug_rad, element='pb', residues=[res], atoms=sugs, name=res.name,
                     chain=res.chain, seq=res.seq, residue_subsection='sugar'),
                Ball(loc=nbas_loc, rad=nbas_rad, element='pb', residues=[res], atoms=nbase, name=res.name,
                     chain=res.chain, seq=res.seq, residue_subsection='nbase')
            ]
        elif res.name in sys.aminos and sys.sc_bb:
            bb_atoms, sc_atoms = [], []
            # Get the back bone atoms and the side chain atoms
            for atom in res.atoms:
                if atom['name'] in sys.amino_bbs:
                    bb_atoms.append(atom)
                elif atom['name'] in sys.amino_scs:
                    sc_atoms.append(atom)
                elif atom['name'] in sys.amino_ignores:
                    continue
                else:
                    sc_bb_add = input("{} atom type not found in amino list (Residue {} # {})! Please choose an option "
                                      "-  Add to 1. Back Bone, 2. Side Chain, 3. Ignore \'{}\', 4. Pass\n>>> "
                                      .format(atom['name'], res.name, res.seq, atom['name']))
                    if sc_bb_add == '1':
                        sys.amino_bbs.append(atom['name'])
                        print("{} will be added to amino backbones from here on out".format(atom['name']))
                    elif sc_bb_add == '2':
                        sys.amino_scs.append(atom['name'])
                        print("{} will be added to amino side chains from here on out".format(atom['name']))
                    elif sc_bb_add == '3':
                        sys.amino_ignores.append(atom['name'])
                        print("{} will be ignored from here on out".format(atom['name']))
                    else:
                        continue
            # Get the balls for the aminos
            bb_loc, bb_rad = make_ball(bb_atoms, sys.scheme, sys.mass_weighted, sys.include_h, sys.therm_cush)
            sc_loc, sc_rad = make_ball(sc_atoms, sys.scheme, sys.mass_weighted, sys.include_h, sys.therm_cush)
            # Create the ball object
            sys.balls += [
                Ball(loc=bb_loc, rad=bb_rad, element=res.elem_col, residues=[res], atoms=bb_atoms, name=res.name,
                     chain=res.chain, seq=res.seq, residue_subsection='bb'),
                Ball(loc=sc_loc, rad=sc_rad, element='pb', residues=[res], atoms=sc_atoms, name=res.name,
                     chain=res.chain, seq=res.seq, residue_subsection='sc')]
        else:
            # Get the loc and rad
            loc, rad = make_ball(res.atoms, sys.scheme, sys.mass_weighted, sys.include_h, sys.therm_cush)
            # Create the ball object
            sys.balls.append(Ball(loc=loc, rad=rad, element=res.elem_col, residues=[res], atoms=res.atoms,
                                  name=res.name, chain=res.chain, seq=res.seq))
