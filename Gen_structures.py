# -*- coding: utf-8 -*-
"""
@Time ： 2023/3/8 22:42
@Auth ： max
@File ：Gen_structures.py
@IDE ：PyCharm
@Motto：ABC(Always Be Coding)
"""
import copy
import linecache
import os

import numpy as np
from ase import Atom
from ase.build import sort
from ase.io import read, write
from math import sqrt, pow

import clean_finder as cf


def get_distance(atom1, atom2):
    x1, y1, z1 = atom1.position
    x2, y2, z2 = atom2.position
    return sqrt(pow((x1 - x2), 2) + pow((y1 - y2), 2) + pow((z1 - z2), 2))


def get_distance_2d(atom1, atom2):
    x1, y1, _ = atom1.position
    x2, y2, _ = atom2.position
    return sqrt(pow((x1 - x2), 2) + pow((y1 - y2), 2))


class GenerateStructures:
    """The structure will be like:
        working path
                |___programs
                |___calculate_files
                |   |___MD_settings
                |       |___lasp.in(structure opt)
                |   |___energy_settings
                |       |___lasp.in(MD)
                |   |___Potential_files
                |       |___XXX.pot
                |___step0
                |   |___Top.txt
                |   |___Sites.txt
                |   |___cal_energy
                |   |   |___slab
                |   |       |_input.arc
                |   |       |-XXX.pot
                |   |       |-lasp.in(structure opt)
                |   |   |___00
                |   |       |-input.arc
                |   |       |-XXX.pot
                |   |       |-lasp.in(structure opt)
                |   |___Energy_list.txt
                |   |___cal_MD
                |   |___MD1
                |       |_MD_energy01.txt
                |       |_input.arc
                |       |_XXX.pot
                |       |_lasp.in(MD)
                |   |___MD2
                |       |_MD_energy02.txt
                |       |_input.arc
                |       |-.pot
                |       |-lasp.in(MD)
                |   |___...
                |   |
                |   |___MDN
                |       |_MD_energyNN.txt
                |       |-input.arc
                |       |-XXX.pot
                |       |-lasp.in(MD)
                |___step2
                |...



    """

    def __init__(self, filename, path, simulation_systems, step, run_type='Oxidation'):
        """load paths"""
        self.current_path = path
        self.type = run_type
        self.step = step  # This parameter will indicate the Oxidation/Reduction steps
        if self.type == 'Oxidation':
            self.filename = self.current_path + '/' + 'step{}'.format(self.step) + '/' + filename
        else:
            self.filename = self.current_path + '/' + 'reduction_step{}'.format(self.step) + '/' + filename
        self.structure_row = read(self.filename)
        self.structure = (sort(self.structure_row, self.structure_row.positions[:, 2]))
        self.first_layer = self.find_top_layer()  # Finished: find a method to find layers
        self.simulation_systems = simulation_systems
        self.sitesall = None
        self.files = 0
        self.MD_generate = 0
        self.continue_MD = True
        self.continue_Reduction = True
        self.remove_o_number = 0
        self.continue_add_o = True

        if self.type == 'Oxidation':
            try:
                os.system('mkdir {}/step{}'.format(self.current_path, self.step))
            except OSError:
                print('mkdir {}/step{} exist'.format(self.current_path, self.step))
                pass
            try:
                os.mkdir('{}/step{}/cal_energy'.format(self.current_path, self.step))
            except OSError:
                print('{}/step{}/cal_energy exist'.format(self.current_path, self.step))
                pass
            try:
                os.mkdir('{}/step{}/cal_MD'.format(self.current_path, self.step))
            except OSError:
                print('{}/step{}/cal_MD exist'.format(self.current_path, self.step))
                pass
            try:
                os.mkdir('{}/step{}/cal_MD/MD{}'.format(self.current_path, self.step, self.MD_generate))
            except OSError:
                print('{}/step{}/cal_MD/MD{} exist'.format(self.current_path, self.step, self.MD_generate))
                pass
            if self.step == 0:
                if os.access('{}/step{}/Hollow.txt'.format(self.current_path, + self.step), os.R_OK):
                    print("Sites exists")
                    pass
                else:
                    self.findsites()
            else:
                if os.access('{}/step{}/Top.txt'.format(self.current_path, + self.step), os.R_OK):
                    print("Sites exists")
                    pass
                else:
                    self.findsites()
        else:  # reduction step
            try:
                os.mkdir('{}/reduction_step{}'.format(self.current_path, self.step))
            except OSError:
                print('{}/reduction_step{} exist'.format(self.current_path, self.step))
                pass
            try:
                os.mkdir('{}/reduction_step{}/cal_energy'.format(self.current_path, self.step))
            except OSError:
                print('{}/reduction_step{}/cal_energy exist'.format(self.current_path, self.step))
                pass
            try:
                os.mkdir('{}/reduction_step{}/cal_MD'.format(self.current_path, self.step))
            except OSError:
                print('{}/reduction_step{}/cal_MD exist'.format(self.current_path, self.step))
                pass

    def find_top_layer(self):
        # for i in self.structure:
        tmp_str = read('{}/row.cif'.format(self.current_path))
        tmp_data = set(tmp_str.positions[:, 2])
        return max(tmp_data) - 1.5
        # self.first_layer = max(tmp_data)

    def findsites(self):
        """Find Top sites and save them into current_path/stepX/Top.txt"""
        finder = cf.Adsorb_finderPro(self.filename)
        finder.angle = 45
        finder.After_settings()
        top = finder.re_T
        t, b, h, _ = finder.FindSitesXYZ()
        # step(self.step)/Top.txt
        np.savetxt("{}/step{}/Top.txt".format(self.current_path, self.step), t, fmt='%f', delimiter=',')
        # np.savetxt("{}/Bridge.txt".format(self.current_path), B, fmt='%f', delimiter=',')
        np.savetxt("{}/step{}/Hollow.txt".format(self.current_path, self.step), h, fmt='%f', delimiter=',')

    def make_energy_folders(self):
        """This function will make Energy calculate folders"""
        if self.step == 0:
            sites = np.loadtxt('{}/step{}/Hollow.txt'.format(self.current_path, self.step), delimiter=',')
            print('First add O2')
        else:
            sites = np.loadtxt('{}/step{}/Top.txt'.format(self.current_path, self.step), delimiter=',')
        tmp = 0
        site = []
        for i in sites:
            if i[2] >= self.first_layer:
                tmp_stru = copy.deepcopy(self.structure)
                # print(tmp_stru)
                site.append(i)
                try:
                    os.mkdir('{}/step{}/cal_energy/{}'.format(self.current_path, self.step, tmp))
                except OSError:
                    pass
                if self.step == 0:
                    tmp_stru.append(Atom('O', [i[0], i[1], i[2] + 1.2]))
                    tmp_stru.append(Atom('O', [i[0], i[1], i[2] + 1.2 + 1.20]))
                else:
                    tmp_stru.append(Atom('O', [i[0], i[1], i[2] + 1.2]))
                # print(tmp_stru)

                try:
                    write('{}/step{}/cal_energy/{}/str.cif'.format(self.current_path, self.step, tmp), tmp_stru)
                except OSError:
                    pass
                try:
                    write('{}/step{}/cal_energy/{}/input.arc'.format(self.current_path, self.step, tmp), tmp_stru,
                          format='dmol-arc')
                except OSError:
                    pass
                try:
                    os.system('cp {}/calculate_files/Potential_files/{} {}/step{}/cal_energy/{}/{}'.format(
                        self.current_path, self.simulation_systems, self.current_path, self.step, tmp,
                        self.simulation_systems))
                except OSError:
                    pass
                try:
                    os.system('cp {}/calculate_files/energy_settings/lasp.in {}/step{}/cal_energy/{}/lasp.in'.format(
                        self.current_path, self.current_path, self.step, tmp))
                except OSError:
                    pass
                tmp = tmp + 1
            else:
                pass
        np.savetxt("{}/step{}/Sites.txt".format(self.current_path, self.step), site, fmt='%f', delimiter=',')

        #     np.savetxt("a.txt", a, fmt = '%f', delimiter = ',')
        try:
            os.mkdir('{}/step{}/cal_energy/slab'.format(self.current_path, self.step))
        except OSError:
            pass
        try:
            write('{}/step{}/cal_energy/slab/str.cif'.format(self.current_path, self.step), self.structure)
            write('{}/step{}/cal_energy/slab/input.arc'.format(self.current_path, self.step), self.structure,
                  format='dmol-arc')
        except OSError:
            pass
        try:
            os.system('cp {}/calculate_files/energy_settings/lasp.in {}/step{}/cal_energy/slab/lasp.in'.format(
                self.current_path, self.current_path, self.step))
            os.system('cp {}/calculate_files/Potential_files/{} {}/step{}/cal_energy/slab/{}'.format(
                self.current_path, self.simulation_systems, self.current_path, self.step,
                self.simulation_systems))
        except OSError:
            pass
        print('generated {} structures, good luck researcher !'.format(tmp))
        self.files = tmp

    def take_o_coverage(self):
        """read energy and add O in sequence, and take a best coverage"""
    # todo: add another step to confirm the the O2 cover rate.
    # 伪代码： 得到能量之后，然后通过能量对O原子位点进行排序（读取现成对文件）然后在依次序加入O原子，重新计算具有覆盖度时的吸附能
        pass
    def make_md_folders(self, chosen):
        """This function will make MD folders
        chosen should be the index of sites
        """
        if self.type == 'Oxidation':
            try:
                os.mkdir('{}/step{}/cal_MD/MD{}'.format(self.current_path, self.step, self.MD_generate))
            except OSError:
                pass
        else:  # Reduction
            try:
                os.mkdir('{}/reduction_step{}/cal_MD/MD{}'.format(self.current_path, self.step, self.MD_generate))
            except OSError:
                pass
        # 添加氧原子的部分捏
        if self.type == 'Oxidation':
            self.sitesall = np.loadtxt('{}/step{}/Sites.txt'.format(
                self.current_path, self.step), delimiter=',')
            add = []
            for i in chosen:
                add.append(self.sitesall[i])
            tmp_stru = copy.deepcopy(self.structure)
            for j in add:
                if self.step == 0:
                    tmp_stru.append(Atom('O', [j[0], j[1], j[2] + 1.20]))
                    tmp_stru.append(Atom('O', [j[0], j[1], j[2] + 1.20 + 1.20]))
                else:
                    tmp_stru.append(Atom('O', [j[0], j[1], j[2] + 1.20]))
            try:
                # if os.path.exists('{}/step{}/cal_MD/MD{}/New_structures.cif'.format(self.current_path, self.step,
                #                                                             self.MD_generate)):
                #     pass
                # else:
                write('{}/step{}/cal_MD/MD{}/New_structures.cif'.format(self.current_path, self.step,
                                                                        self.MD_generate), tmp_stru)
            except OSError:
                pass
            try:
                # if os.path.exists('{}/step{}/cal_MD/MD{}/input.arc'.format(self.current_path, self.step,
                # self.MD_generate)): pass else:
                write('{}/step{}/cal_MD/MD{}/input.arc'.format(self.current_path, self.step, self.MD_generate),
                      tmp_stru, format='dmol-arc')
            except OSError:
                pass
            try:
                # if os.path.exists('cp {}/calculate_files/Potential_files/{} {}/step{}/cal_MD/MD{}/{}'.format(
                #         self.current_path, self.simulation_systems, self.current_path, self.step, self.MD_generate,
                #         self.simulation_systems)):
                #     pass
                # else:
                os.system('cp {}/calculate_files/Potential_files/{} {}/step{}/cal_MD/MD{}/{}'.format(
                    self.current_path, self.simulation_systems, self.current_path, self.step, self.MD_generate,
                    self.simulation_systems))
            except OSError:
                pass
            try:
                # if os.path.exists('cp {}/calculate_files/MD_settings/lasp.in {}/step{}/cal_MD/MD{}/lasp.in'.format(
                #         self.current_path, self.current_path, self.step, self.MD_generate)):
                #     pass
                # else:
                os.system('cp {}/calculate_files/MD_settings/lasp.in {}/step{}/cal_MD/MD{}/lasp.in'.format(
                    self.current_path, self.current_path, self.step, self.MD_generate))
            except OSError:
                pass
            # self.MD_generate = self.MD_generate + 1
            print('MD : {} generate'.format(self.MD_generate))
        else:
            self.sitesall = np.loadtxt('{}/reduction_step{}/Removed_O_list.txt'.format(
                self.current_path, self.step), delimiter=',')
            tmp_stru = copy.deepcopy(self.structure)
            remove_list = []
            # Bug fixed: when reduction step is in the end, the remove list will be [], which will cause index error
            try:
                for i in chosen:
                    i = int(i) - 1
                    print(i, self.sitesall[i])
                    remove_list.append(self.sitesall[i])
                tmp_stru = copy.deepcopy(self.structure)
                remove_list.sort(reverse=True)
                print(remove_list)
                # remove_list.sort(revers=True)
                self.remove_o_number = len(remove_list)
                for i in remove_list:
                    i = int(i)
                    del tmp_stru[i]
                try:
                    write('{}/reduction_step{}/cal_MD/MD{}/New_structures.cif'.format(self.current_path, self.step,
                                                                                      self.MD_generate), tmp_stru)
                except OSError:
                    pass
                try:
                    write(
                        '{}/reduction_step{}/cal_MD/MD{}/input.arc'.format(self.current_path, self.step,
                                                                           self.MD_generate),
                        tmp_stru,
                        format='dmol-arc')
                except OSError:
                    pass
                try:
                    os.system('cp {}/calculate_files/Potential_files/{} {}/reduction_step{}/cal_MD/MD{}/{}'.format(
                        self.current_path, self.simulation_systems, self.current_path, self.step, self.MD_generate,
                        self.simulation_systems))
                except OSError:
                    pass
                try:
                    os.system(
                        'cp {}/calculate_files/MD_settings/Reduction/lasp.in {}/reduction_step{}/cal_MD/MD{}/lasp.in'.format(
                            self.current_path, self.current_path, self.step, self.MD_generate))
                except OSError:
                    pass
                # self.MD_generate = self.MD_generate + 1
                print('MD-reduction {} generate'.format(self.MD_generate))
            except IndexError:
                self.remove_o_number = 0
                pass

    def get_new_structures(self):
        if self.type == 'Oxidation':
            file = open('{}/step{}/cal_MD/MD{}/md.arc'.format(self.current_path, self.step, self.MD_generate), 'r')
            line1 = file.readlines()
            file2 = open('{}/step{}/cal_MD/MD{}/input.arc'.format(self.current_path, self.step, self.MD_generate), 'r')
            line2 = file2.readlines()
        else:
            file = open('{}/reduction_step{}/cal_MD/MD{}/md.arc'.format(self.current_path, self.step, self.MD_generate),
                        'r')
            line1 = file.readlines()
            file2 = open(
                '{}/reduction_step{}/cal_MD/MD{}/input.arc'.format(self.current_path, self.step, self.MD_generate), 'r')
            line2 = file2.readlines()
        stru = []
        if self.type == 'Oxidation':
            for i in range(len(line1) - len(line2) + 4, len(line1) + 1):
                tmp = linecache.getline('{}/step{}/cal_MD/MD{}/md.arc'.format(self.current_path, self.step,
                                                                              self.MD_generate), i)
                stru.append(tmp)
            self.MD_generate = self.MD_generate + 1
            try:
                os.mkdir('{}/step{}/cal_MD/MD{}'.format(self.current_path, self.step, self.MD_generate))
            except OSError:
                pass
            f = open('{}/step{}/cal_MD/MD{}/input.arc'.format(self.current_path, self.step, self.MD_generate), 'w')
        else:
            for i in range(len(line1) - len(line2) + 4, len(line1) + 1):
                tmp = linecache.getline('{}/reduction_step{}/cal_MD/MD{}/md.arc'.format(self.current_path, self.step,
                                                                                        self.MD_generate), i)
                stru.append(tmp)
            self.MD_generate = self.MD_generate + 1
            try:
                os.mkdir('{}/reduction_step{}/cal_MD/MD{}'.format(self.current_path, self.step, self.MD_generate))
            except OSError:
                pass
            f = open('{}/reduction_step{}/cal_MD/MD{}/input.arc'.format(self.current_path, self.step, self.MD_generate),
                     'w')
        f.write('!BIOSYM archive 2')
        f.write('\n')
        # print('!BIOSYM archive 2')
        f.write('PBC=ON')
        f.write('\n')
        for i in stru:
            f.write(i)
        f.close()
        # Bug fixed: 2023 03 13: This method will reduce one line in input.arc files each times
        # So we have to read and re-write input.arc files
        if self.type == 'Oxidation':
            tmp = read('{}/step{}/cal_MD/MD{}/input.arc'.format(self.current_path, self.step, self.MD_generate),
                       format='dmol-arc')
            write('{}/step{}/cal_MD/MD{}/input.arc'.format(self.current_path, self.step, self.MD_generate),
                  tmp, format='dmol-arc')
        else:
            tmp = read(
                '{}/reduction_step{}/cal_MD/MD{}/input.arc'.format(self.current_path, self.step, self.MD_generate),
                format='dmol-arc')
            write('{}/reduction_step{}/cal_MD/MD{}/input.arc'.format(self.current_path, self.step, self.MD_generate),
                  tmp, format='dmol-arc')
        if self.MD_generate != 0:
            if self.type == 'Oxidation':
                os.system('cp {}/step{}/cal_MD/MD{}/md.restart {}/step{}/cal_MD/MD{}/md.restart'.format(
                    self.current_path, self.step, self.MD_generate - 1, self.current_path, self.step, self.MD_generate))
                print('cp {}/step{}/cal_MD/MD{}/md.restart {}/step{}/cal_MD/MD{}/md.restart'.format(
                    self.current_path, self.step, self.MD_generate - 1, self.current_path, self.step, self.MD_generate))
            else:
                os.system(
                    'cp {}/reduction_step{}/cal_MD/MD{}/md.restart {}/reduction_step{}/cal_MD/MD{}/md.restart'.format(
                        self.current_path, self.step, self.MD_generate - 1, self.current_path, self.step,
                        self.MD_generate))
                print('cp {}/reduction_step{}/cal_MD/MD{}/md.restart {}/reduction_step{}/cal_MD/MD{}/md.restart'.format(
                    self.current_path, self.step, self.MD_generate - 1, self.current_path, self.step, self.MD_generate))
        else:
            pass
        try:
            if self.type == 'Oxidation':
                os.system('cp {}/calculate_files/Potential_files/{} {}/step{}/cal_MD/MD{}/{}'.format(
                    self.current_path, self.simulation_systems, self.current_path, self.step, self.MD_generate,
                    self.simulation_systems))
            else:
                os.system('cp {}/calculate_files/Potential_files/{} {}/reduction_step{}/cal_MD/MD{}/{}'.format(
                    self.current_path, self.simulation_systems, self.current_path, self.step, self.MD_generate,
                    self.simulation_systems))
        except OSError:
            pass
        try:
            if self.type == 'Oxidation':
                os.system('cp {}/calculate_files/MD_settings/lasp.in {}/step{}/cal_MD/MD{}/lasp.in'.format(
                    self.current_path, self.current_path, self.step, self.MD_generate))
            else:
                os.system(
                    'cp {}/calculate_files/MD_settings/Reduction/lasp.in {}/reduction_step{}/cal_MD/MD{}/lasp.in'.format(
                        self.current_path, self.current_path, self.step, self.MD_generate))
        except OSError:
            pass
        try:  # restart
            if self.type == 'Oxidation':
                os.system('cp {}/calculate_files/MD_settings/lasp.in {}/step{}/cal_MD/MD{}/lasp.in'.format(
                    self.current_path, self.current_path, self.step, self.MD_generate))
            else:
                os.system(
                    'cp {}/calculate_files/MD_settings/Reduction/lasp.in {}/reduction_step{}/cal_MD/MD{}/lasp.in'.format(
                        self.current_path, self.current_path, self.step, self.MD_generate))
        except OSError:
            pass
        # It will save the oxygen Not removed structures
        # It will save in input{step+1}.arc
        if self.continue_MD:  # if continue_MD = True, stop and run next step;
            # if continue_MD = False, means that the MD stop, export files to
            # self.current_path + '/' + 'step{}'.format(self.step) + '/' + filename
            print('MD{} generate'.format(self.MD_generate))
            pass
        else:
            if self.type == 'Oxidation':
                try:
                    os.mkdir('{}/step{}'.format(self.current_path, self.step + 1))
                except OSError:
                    pass
            else:
                try:
                    os.mkdir('{}/reduction_step{}'.format(self.current_path, self.step + 1))
                except OSError:
                    pass
            if self.type == 'Oxidation':
                try:
                    os.mkdir('{}/step{}/cal_energy'.format(self.current_path, self.step + 1))
                except OSError:
                    pass
            else:
                try:
                    os.mkdir('{}/reduction_step{}/cal_energy'.format(self.current_path, self.step + 1))
                except OSError:
                    pass
            if self.type == 'Oxidation':
                try:
                    os.mkdir('{}/step{}/cal_MD'.format(self.current_path, self.step + 1))
                except OSError:
                    pass
            else:
                try:
                    os.mkdir('{}/reduction_step{}/cal_MD'.format(self.current_path, self.step + 1))
                except OSError:
                    pass
            if self.type == 'Oxidation':
                write('{}/step{}/Un_removed_finished.cif'.format(self.current_path, self.step + 1),
                      read('{}/step{}/cal_MD/MD{}/input.arc'.format(self.current_path, self.step, self.MD_generate),
                           format='dmol-arc'))  # after MD, save a cif structure
            else:
                write('{}/reduction_step{}/Un_removed_finished.cif'.format(self.current_path, self.step + 1),
                      read('{}/reduction_step{}/cal_MD/MD{}/input.arc'.format(self.current_path, self.step,
                                                                              self.MD_generate),
                           format='dmol-arc'))
            # It will save the final inputX.cif
            # if self.step == 0:  # finished 此步骤的效果存在疑问，需要进一步验证，在beta-edition中做去掉处理
            #     write('{}/step{}/input{}.cif'.format(self.current_path, self.step + 1, self.step + 1),
            #           read('{}/step{}/Un_removed_finished.cif'.format(self.current_path, self.step + 1)))
            # else:
            if self.type == 'Oxidation':
                row_stru = read('{}/step{}/Un_removed_finished.cif'.format(self.current_path, self.step + 1))
                tmp = self.remove_free_oxygen(row_stru)
                write('{}/step{}/input{}.cif'.format(self.current_path, self.step + 1, self.step + 1), tmp)
            else:
                row_stru = read('{}/reduction_step{}/Un_removed_finished.cif'.format(self.current_path, self.step + 1))
                tmp = row_stru
                write('{}/reduction_step{}/input{}.cif'.format(self.current_path, self.step + 1, self.step + 1), tmp)
            self.MD_generate = self.MD_generate + 1

    def remove_free_oxygen(self, stru):
        """use a triangle matrix to save distance
        0   1   2   ... N
        1   d11
        2   d21 d22
        ... ... ...
        N   dN1 dN2     dNN

        the coordinate like:[1 2 3 4 ... ]
        the elements like: [O Cu Zn ... ]
        test all distance between the O-Cu or O-Zn, if they are all larger than 2.5, remove this O atom
        if Oxygen atom are too close to other O atoms (less than 1.5 * 0.6), the O will be removed
        """
        tmp_strut = stru
        # tmp_strut = read('{}/step{}/input{}.arc'.format(self.current_path, self.step, self.step + 1),
        #                 format='dmol-arc')
        # for test
        # tmp_strut = read('{}/input{}.cif'.format(self.current_path, self.step + 1))

        atoms_list = []
        remove_list = []
        for k in tmp_strut:
            atoms_list.append(k.symbol)
        triangle = np.zeros((len(tmp_strut), len(tmp_strut)))
        for i in range(0, len(tmp_strut)):
            for j in range(0, i):
                triangle[i, j] = get_distance(tmp_strut[i], tmp_strut[j])
        atoms_list = np.array(atoms_list)
        O_list = np.where((atoms_list == 'O'))[0]  # 得到了O的坐标位置
        # print(O_list)
        Metal_list = np.where((atoms_list != 'O'))[0]
        # print(Metal_list)
        # remove too far
        # bug fix: 2023 3 12: for some O. they might adsorb on the bottom of the cell due to the cyclic
        for top in range(0, len(tmp_strut)):
            if (tmp_strut[top].position[2] >= tmp_strut.cell[2][2] - 1.5 or tmp_strut[top].position[2] <= 1.5) and \
                    tmp_strut[top].symbol == 'O':
                # print('too high', top)
                remove_list.append(top)

        for r_O in O_list:
            O_M_distance = []
            for r_M in Metal_list:
                if r_M < r_O:
                    O_M_distance.append(triangle[r_O][r_M])
                elif r_M >= r_O:
                    O_M_distance.append(triangle[r_M][r_O])
            #  ###
            # print(r_O)
            # if r_O == 181:
            #     # print(min(O_M_distance))
            if min(O_M_distance) >= 2.5:
                remove_list.append(r_O)
            # print(O_M_distance)

        # remove too close
        for r_O1 in O_list:
            O_O_distance = []
            for r_O2 in O_list:
                if r_O1 > r_O2:
                    O_O_distance.append(triangle[r_O1][r_O2])
                elif r_O1 < r_O2:
                    O_O_distance.append(triangle[r_O2][r_O1])
            if min(O_O_distance) <= 1.5 * 0.6:
                remove_list.append(r_O1)
            # print(O_O_distance)
        remove_list = list(set(remove_list))
        # print(remove_list)
        # print(O_O_distance)
        # for i in a:
        #     print(i)
        remove_list.sort(reverse=True)
        # sort to avoid change the index of atom
        print('row structure', tmp_strut)
        for remove in remove_list:
            del tmp_strut[remove]
            # print(remove)

        print('removed structure', tmp_strut)
        # write('{}/removed_input.cif'.format(self.current_path), tmp_strut)
        strut_removed = tmp_strut
        # print(remove_list)
        # print(triangle[182][124])
        #     print(k)
        # write('{}/fixed.cif'.format(self.current_path),strut_removed)
        return strut_removed

    def remove_O(self):

        print('aaaaa')

    def make_reduction_folders(self):
        O_list = []
        ####
        # finished: we need a method to find the surface Oxygen
        # 加一个距离矩阵，用来计算距离，对于O原子，其距离较近的一圈内部没有比它高的原子，就是表面O原子
        tmp_atoms = copy.deepcopy(self.structure)

        atoms_list = []
        for tmp_atom in tmp_atoms:
            atoms_list.append(tmp_atom.symbol)

        triangle = np.zeros((len(tmp_atoms), len(tmp_atoms)))
        for tmp_atom in range(0, len(tmp_atoms)):
            for j in range(0, tmp_atom):
                triangle[tmp_atom, j] = get_distance_2d(tmp_atoms[tmp_atom], tmp_atoms[j])
        print(triangle)

        surface_list = []
        atoms_list = np.array(atoms_list)
        O_list = np.where((atoms_list == 'O'))[0]  # 得到了O的坐标位置
        Metal_list = np.where((atoms_list != 'O'))[0]
        # judgement: find atoms close to O atoms and these atoms are lower than O atoms return these O atoms
        for O_atom in O_list:
            tmp_close = []
            for atom in range(0, len(atoms_list)):
                if O_atom >= atom:
                    if triangle[O_atom][atom] <= 4:
                        if tmp_atoms[atom].position[2] >= tmp_atoms[O_atom].position[2] + 1:
                            tmp_close.append(tmp_atoms[atom].symbol)
                else:
                    #             print(triangle[atom][O])
                    if triangle[atom][O_atom] <= 4:
                        if tmp_atoms[atom].position[2] >= tmp_atoms[O_atom].position[2] + 1:
                            tmp_close.append(tmp_atoms[atom].symbol)
            if tmp_close == []:
                surface_list.append(O_atom)
        surface_list.sort(reverse=True)
        # str_tmp = copy.deepcopy(self.structure)
        # for i in surface_list:
        #     del str_tmp[i]
        # write('try.cif', str_tmp)
        ####
        O_list = surface_list
        # for i in self.structure:
        #     if i.symbol == 'O':
        #         O_list.append(i.index)
        # O_list.sort(reverse=True)
        #####

        # if len(O_list) >= 20:
        # for j in O_list:
        #     if O_
        try:
            os.mkdir('{}/reduction_step{}/cal_energy/slab'.format(self.current_path, self.step))
        except OSError:
            pass
        write('{}/reduction_step{}/cal_energy/slab/input.arc'.format(self.current_path, self.step), self.structure,
              format='dmol-arc')
        write('{}/reduction_step{}/cal_energy/slab/str.cif'.format(self.current_path, self.step), self.structure)
        try:
            os.system(
                'cp {}/calculate_files/energy_settings/Reduction/lasp.in {}/reduction_step{}/cal_energy/slab/lasp.in'.format(
                    self.current_path, self.current_path, self.step))
        except OSError:
            pass
        try:
            os.system('cp {}/calculate_files/Potential_files/{} {}/reduction_step{}/cal_energy/slab/{}'.format(
                self.current_path, self.simulation_systems, self.current_path, self.step, self.simulation_systems))
        except OSError:
            pass
        for j in range(0, len(O_list)):
            tmp_stru = copy.deepcopy(self.structure)
            try:
                os.mkdir('{}/reduction_step{}/cal_energy/{}'.format(self.current_path, self.step, j))
            except OSError:
                print('{}/reduction_step{}/cal_energy/{} exist'.format(self.current_path, self.step, j))
                pass
            del tmp_stru[O_list[j]]
            write('{}/reduction_step{}/cal_energy/{}/input.arc'.format(self.current_path, self.step, j), tmp_stru,
                  format='dmol-arc')
            write('{}/reduction_step{}/cal_energy/{}/str.cif'.format(self.current_path, self.step, j), tmp_stru)
            try:
                os.system(
                    'cp {}/calculate_files/energy_settings/Reduction/lasp.in {}/reduction_step{}/cal_energy/{}/lasp.in'.format(
                        self.current_path, self.current_path, self.step, j))
            except OSError:
                pass
            try:
                os.system('cp {}/calculate_files/Potential_files/{} {}/reduction_step{}/cal_energy/{}/{}'.format(
                    self.current_path, self.simulation_systems, self.current_path, self.step, j,
                    self.simulation_systems))
            except OSError:
                pass
        np.savetxt("{}/reduction_step{}/Removed_O_list.txt".format(self.current_path, self.step), O_list, fmt='%f',
                   delimiter=',')
        print('{} structures have been generated, Good luck researcher'.format(len(O_list)))
        self.files = len(O_list)

#
# # test
# path = '/data/home/max/Projects/OxideCuZn/NemethodCu3Zn/tryauto_safe'
#
# # os.mkdir('{}/reduction_step{}'.format(path,1))
# # write('{}/reduction_step{}/Oxide{}.cif'.format(path, 1, 1), read('{}/Oxide.cif'.format(path)))
# try_class = GenerateStructures('Oxide0.cif', path, 'CuZnCHO.pot', 0, type = '114514')
# try_class.make_reduction_folders()
# # try_class.make_md_folders([11,6,9])
# # # try_class.make_energy_folders()make_energy_folders
# # struc = read('/data/home/max/Projects/OxideCuZn/NemethodCu3Zn/Tryauto/input.cif')
# # print('find {}/step{} -name Sites.txt'.format(try_class.current_path, + try_class.step))
# # os.system('find {}/step{} -name Sites.txt'.format(try_class.current_path, + try_class.step))
# # try_class.make_md_folders(list([1, 2, 3, 4, 5, 6, 7, 8]))
# try_class.remove_free_oxygen()
