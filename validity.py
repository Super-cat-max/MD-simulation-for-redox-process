from collections import Counter
from itertools import product

import numpy as np
import smact
from pymatgen.core import Structure
from smact.screening import pauling_test

chemical_symbols = [
    'X', 'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg',
    'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn',
    'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb',
    'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In',
    'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm',
    'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta',
    'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At',
    'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk',
    'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt',
    'Ds', 'Rg', 'Cn', 'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og'
]


def is_valid(structure, use_pauling_test=True, include_alloys=True):
    try:
        if Structure.is_valid(structure):
            elem_counter = Counter(structure.atomic_numbers)
            composition = [(elem, elem_counter[elem])
                           for elem in sorted(elem_counter.keys())]
            elems, counts = list(zip(*composition))
            counts = np.array(counts)
            counts = counts / np.gcd.reduce(counts)
            comps = tuple(counts.astype('int').tolist())

            def smact_validity(comp, count, use_pauling_test, include_alloys):
                elem_symbols = tuple([chemical_symbols[elem] for elem in comp])
                space = smact.element_dictionary(elem_symbols)
                smact_elems = [e[1] for e in space.items()]
                electronegs = [e.pauling_eneg for e in smact_elems]
                ox_combos = [e.oxidation_states for e in smact_elems]
                if len(set(elem_symbols)) == 1:
                    return True
                if include_alloys:
                    return all(
                        [elem_s in smact.metals for elem_s in elem_symbols])

                threshold = np.max(count)
                compositions = []
                for ox_states in product(*ox_combos):
                    stoichs = [(c,) for c in count]
                    # Test for charge balance
                    cn_e, cn_r = smact.neutral_ratios(ox_states,
                                                      stoichs=stoichs,
                                                      threshold=threshold)
                    # Electronegativity test
                    if cn_e:
                        if use_pauling_test:
                            try:
                                pauling_test(ox_states, electronegs)
                            except TypeError:
                                # if no electronegativity data, assume it is okay
                                for ratio in cn_r:
                                    compositions.append(
                                        tuple([elem_symbols, ox_states,
                                               ratio]))
                        else:
                            for ratio in cn_r:
                                compositions.append(
                                    tuple([elem_symbols, ox_states, ratio]))
                compositions = [(i[0], i[2]) for i in compositions]
                compositions = list(set(compositions))
                return len(compositions) > 0

            return smact_validity(elems, comps, use_pauling_test,
                                  include_alloys)
        else:
            return False
    except:
        return False
