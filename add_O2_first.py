# -*- coding: utf-8 -*-
"""
@Time ： 2023/3/10 15:54
@Auth ： max
@File ：add_O2_first.py
@IDE ：PyCharm
@Motto：ABC(Always Be Coding)
"""
import linecache
import os
import subprocess

from ase.build import sort
from ase.io import read, write


class add_O2:
    def __init__(self, path, environment, cmd):
        self.path = path
        self.environment = environment
        self.structure_row = read('{}/row.cif'.format(self.path))
        self.structure = (sort(self.structure_row, self.structure_row.positions[:, 2]))
        self.command = cmd

    def relax_row(self):
        try:
            os.mkdir('{}/relax_row'.format(self.path))
        except OSError:
            pass
        try:
            os.mkdir('{}/step0'.format(self.path))
        except OSError:
            pass
        write('{}/relax_row/input.arc'.format(self.path), self.structure, format='dmol-arc')
        try:
            os.system('cp {}/calculate_files/Potential_files/{} {}/relax_row/{}'.format(
                self.path, self.environment, self.path, self.environment))
        except OSError:
            pass
        try:
            os.system('cp {}/calculate_files/MD_settings/lasp.in {}/relax_row/lasp.in'.format(
                self.path, self.path))
        except OSError:
            pass

    def run_md(self):
        subprocess.call(self.command, shell=True, cwd='{}/relax_row'.format(self.path))
        print('MD finished')

    def get_new_structures(self):
        file = open('{}/relax_row/md.arc'.format(self.path), 'r')
        line1 = file.readlines()
        file2 = open('{}/relax_row/input.arc'.format(self.path), 'r')
        line2 = file2.readlines()
        stru = []
        for i in range(len(line1) - len(line2) + 4, len(line1) + 1):
            tmp = linecache.getline('{}/relax_row/md.arc'.format(self.path), i)
            stru.append(tmp)
        f = open('{}/step0/input0.arc'.format(self.path), 'w')
        f.write('!BIOSYM archive 2')
        f.write('\n')
        # print('!BIOSYM archive 2')
        f.write('PBC=ON')
        f.write('\n')
        for i in stru:
            f.write(i)
        f.close()
        write('{}/step0/input0.cif'.format(self.path),
              read('{}/step0/input0.arc'.format(self.path), format='dmol-arc'))

    def run(self):
        self.relax_row()
        self.run_md()
        self.get_new_structures()
# # try
# path = '/data/home/max/Projects/OxideCuZn/NemethodCu3Zn/Tryauto'
# command = 'mpirun  /data/software/lasp/3.2.0/NN_pro3.2.0_intel18/Src/lasp > run.log 2>&1'
# testmiao = add_O2(path, 'CuZnO.pot', command)
# # testmiao.relax_row()
# # testmiao.run_md()
# # testmiao.get_new_structures()
