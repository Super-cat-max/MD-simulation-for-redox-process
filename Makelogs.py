# -*- coding: utf-8 -*-
"""
@Time ： 2023/3/8 22:41
@Auth ： max
@File ：Makelogs.py
@IDE ：PyCharm
@Motto：ABC(Always Be Coding)

"""
from ase.io import read,write
from pymatgen.core import Structure
a = read('POSCAR')
a.remove()