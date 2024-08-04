# -*- coding: utf-8 -*-
"""
@Time ： 2023/3/8 22:43
@Auth ： max
@File ：Read_energy.py
@IDE ：PyCharm
@Motto：ABC(Always Be Coding)
"""
import subprocess

import numpy as np


def get_slope(data):
    x = []
    # x 应该在10的数量级上，与Y的变化幅度相对应
    fix_large = 1
    if len(data) <= 999:
        fix_large = 10
    elif len(data) <= 9999:
        fix_large = 100
    elif len(data) <= 99999:
        fix_large = 1000
    elif len(data) <= 999999:
        fix_large = 10000
    elif len(data) <= 9999999:
        fix_large = 100000
    for i in range(0, len(data)):
        x.append(i / fix_large)
    slope, intercept = np.polyfit(x, data, 1)
    return slope, intercept


class calculator_lasp:
    def __init__(self, lasp_command, path, step, number, run_type='Oxidation'):
        self.command = lasp_command
        self.path = path
        self.step = step
        self.number = number
        self.energy_data = None
        self.slab_energy = None
        self.o_energy = -4.772286  # O2 energy(calculated by lasp) divide 2   CuZnO.pot : -4.940347
        # CuZnCHO.pot : -4.772286     When you use different potential files, remembered to change this energy
        self.site_list = []
        self.delete_list = ['allfor.arc', 'allkeys.log', 'allstr.arc',
                            'Badstr.arc', '*.pot']
        self.MD_energy = []
        # self.early_stop = False
        self.generate = 0
        self.type = run_type

    def energy_calculate(self):
        for i in range(0, self.number - 1):
            # subprocess.call(self.command, cwd=)
            if self.type == 'Oxidation':
                subprocess.call(self.command, shell=True, cwd=self.path + '/step{}/cal_energy/{}'.format(self.step, i))
            else:
                subprocess.call(self.command, shell=True, cwd=self.path + '/reduction_step{}/cal_energy/{}'.format(
                    self.step, i))
            # remove unneeded files
            try:
                if self.type == 'Oxidation':
                    for cmd in self.delete_list:
                        rm_cmd = 'rm {}'.format(cmd)
                        subprocess.call(rm_cmd, shell=True,
                                        cwd=self.path + '/step{}/cal_energy/{}'.format(self.step, i))
                else:
                    for cmd in self.delete_list:
                        rm_cmd = 'rm {}'.format(cmd)
                        subprocess.call(rm_cmd, shell=True, cwd=self.path + '/reduction_step{}/cal_energy/{}'.
                                        format(self.step, i))
            except OSError:
                pass
            print('str{} finished'.format(i))
        if self.type == 'Oxidation':
            subprocess.call(self.command, shell=True, cwd=self.path + '/step{}/cal_energy/slab'.format(self.step))
        else:
            subprocess.call(self.command, shell=True,
                            cwd=self.path + '/reduction_step{}/cal_energy/slab'.format(self.step))
        if self.type == 'Oxidation':
            try:
                for cmd in self.delete_list:
                    rm_cmd = 'rm {}'.format(cmd)
                    subprocess.call(rm_cmd, shell=True, cwd=self.path + '/step{}/cal_energy/slab'.format(self.step))
            except OSError:
                pass
        else:
            try:
                for cmd in self.delete_list:
                    rm_cmd = 'rm {}'.format(cmd)
                    subprocess.call(rm_cmd, shell=True,
                                    cwd=self.path + '/reduction_step{}/cal_energy/slab'.format(self.step))
            except OSError:
                pass

        print('slab finished')

    def MD_calculate(self):
        pass

    def read_energy(self):
        energy_list = []
        for i in range(0, self.number - 1):
            if self.type == 'Oxidation':
                tmp_path = '{}/step{}/cal_energy/{}/all.arc'.format(self.path, self.step, i)
            else:
                tmp_path = '{}/reduction_step{}/cal_energy/{}/all.arc'.format(self.path, self.step, i)
            with open(tmp_path, 'r') as x:
                energy_line = x.readlines()[2]
                energy = energy_line.split()[3]
                try:  # 2023/3/11 Bug_fix: if the structures are not in the pot_file, the LASP will
                    # return a file : ExceedSym.arc while the energy is '********' which will make this
                    # program crash, we should avoid it and raise a warning
                    energy_list.append([int(i), float(energy)])
                except OSError:
                    if self.type == 'Oxidation':
                        print('Some thing wrong in structure optimization in step{} cal_energy/{}, if you find '
                              '\'ExceedSym.arc\' in it, you should take care because the NN-Potential is not '
                              'accurate'.format(self.step, i))
                    else:
                        print(
                            'Some thing wrong in structure optimization in reduction_step{} cal_energy/{}, if you find '
                            '\'ExceedSym.arc\' in it, you should take care because the NN-Potential is not '
                            'accurate'.format(self.step, i))
                    pass
                # print(energy)
                x.close()
        energy_save = np.array(energy_list)
        if self.type == 'Oxidation':
            np.savetxt("{}/step{}/Energy.txt".format(self.path, self.step), energy_save, fmt='%f', delimiter=',')
        else:
            np.savetxt("{}/reduction_step{}/Energy.txt".format(self.path, self.step), energy_save, fmt='%f',
                       delimiter=',')

    def calculate_adsorb_energy(self):
        if self.type == 'Oxidation':
            self.energy_data = np.loadtxt('{}/step{}/Energy.txt'.format(self.path, self.step), delimiter=',')
        else:
            self.energy_data = np.loadtxt('{}/reduction_step{}/Energy.txt'.format(self.path, self.step), delimiter=',')
        if self.type == 'Oxidation':
            tmp_path2 = '{}/step{}/cal_energy/slab/all.arc'.format(self.path, self.step)
        else:
            tmp_path2 = '{}/reduction_step{}/cal_energy/slab/all.arc'.format(self.path, self.step)
        with open(tmp_path2, 'r') as x:
            energy_line = x.readlines()[2]
            energy = energy_line.split()[3]
            self.slab_energy = float(energy)
            x.close()
        if self.energy_data.shape == (2,):
            self.energy_data = np.array([self.energy_data])
        else:
            pass
        for i in self.energy_data:
            if self.type == 'Oxidation':
                if self.step == 0:
                    i[1] = i[1] - self.slab_energy - self.o_energy * 2
                else:
                    i[1] = i[1] - self.slab_energy - self.o_energy
                i[0] = int(i[0])
            else:
                i[1] = i[1] - self.slab_energy + self.o_energy
                i[0] = int(i[0])
        # print(self.energy_data)
        if self.type == 'Oxidation':
            np.savetxt("{}/step{}/Energy_O_adsorb.txt".format(self.path, self.step), self.energy_data, fmt='%f',
                       delimiter=',')
        else:
            np.savetxt("{}/reduction_step{}/O_formed_adsorb.txt".format(self.path, self.step), self.energy_data,
                       fmt='%f',
                       delimiter=',')

    def get_right_sites(self):
        judgement = np.loadtxt('{}/step{}/Energy_O_adsorb.txt'.format(self.path, self.step), delimiter=',')
        if self.step != 0:
            for i in judgement:
                if i[1] <= 0 and abs(i[1]) <= 3:
                    self.site_list.append(int(i[0]))
        else:  # O_2 steps
            judgement = judgement[judgement[:, 1].argsort()]
            for j in range(0, int(len(judgement) / 5)):
                # print(judgement[j][0])
                self.site_list.append(int(judgement[j][0]))

        self.site_list.sort(reverse=True)
        print('before', self.site_list)
        if len(self.site_list) <= 40:  # 一次加入的O有点太少了，改一下捏捏捏， 采用0.9 的参数来去掉一两个无用的数据
            if not self.site_list[0:int(len(self.site_list) * 0.9)]:
                self.site_list = self.site_list[0:int(len(self.site_list) * 0.9)]
            else:  # 不能为0，即使只有一个元素也要绳之以法！
                self.site_list = self.site_list
        else:
            self.site_list = self.site_list[0:int(len(self.site_list) * 0.5)]
        print('after', self.site_list)

    def get_remove_o_atoms(self):
        judgement = np.loadtxt('{}/reduction_step{}/O_formed_adsorb.txt'.format(self.path, self.step), delimiter=',')
        # judgement = judgement
        # print('judgement judgement', judgement)
        # data = data[np.argsort(data[:, 0])]
        # Bug fixed:2023/03/26
        try:
            if judgement.shape == (2,):
                judgement = judgement
            else:
                try:
                    judgement = judgement[np.argsort(judgement[:, 1])]
                except OSError:
                    judgement = []
        except OSError:
            judgement = []
        except IndexError:
            judgement = []
        # print('sort', judgement)
        # print(judgement)
        # print(judgement)
        # for i in judgement:
        #     print(i[0])
        remove_list = []
        print('judgement', judgement)
        if judgement == []:  # do not remake the type, it will cause some error
            pass
        else:
            if judgement.shape == (2,):
                judgement = np.array([judgement])
            else:
                pass
        if int(len(judgement)) <= 20:
            if 5 <= int(len(judgement)) and judgement != []:
                for i in range(0, 5):
                    remove_list.append(judgement[i][0])
            elif judgement != []:
                for i in range(0, int(len(judgement))):
                    remove_list.append(judgement[i][0])
            elif judgement == []:
                remove_list = []
        else:
            for i in range(0, int(len(judgement) * 0.5)):
                remove_list.append(judgement[i][0])
        # print(remove_list)
        # print('miao')
        print('{} o_atoms will be removed'.format(len(remove_list)))
        self.site_list = remove_list
        # return remove_list

    def run_md(self, generate):
        # finished: (V1.2): I made It! I am a good cat!
        # read lasp.out in cal_MD
        # try to use split method to reduce the time of the MD process
        # each step will cost X ns, and if the energy between two steps are close
        # stop the MD process, it will save lots of times
        if self.type == 'Oxidation':
            subprocess.call(self.command, shell=True,
                            cwd='{}/step{}/cal_MD/MD{}'.format(self.path, self.step, generate))
        else:
            subprocess.call(self.command, shell=True,
                            cwd='{}/reduction_step{}/cal_MD/MD{}'.format(self.path, self.step, generate))
        print('MD finished')

    def run_energy(self):
        self.energy_calculate()
        self.read_energy()
        self.calculate_adsorb_energy()
        if self.type == 'Oxidation':
            self.get_right_sites()
        else:
            self.get_remove_o_atoms()

    # def judge_early_stop(self):
    # finished make judgement function
    # self.early_stop = True
    # self.early_stop = False

    def read_MD_energy(self, generate):
        self.MD_energy = []
        if self.type == 'Oxidation':
            f = open('{}/step{}/cal_MD/MD{}/lasp.out'.format(self.path, self.step, generate))
        else:
            f = open('{}/reduction_step{}/cal_MD/MD{}/lasp.out'.format(self.path, self.step, generate))
        lines = f.readlines()
        start = 0
        for i in range(0, len(lines) - 1):
            if lines[i].find('Equilibrated') != -1:
                start = i
        for j in range(start + 1, len(lines) - 1):
            self.MD_energy.append(float((lines[j].split()[1])))
        f.close()
        if self.type == 'Oxidation':
            np.savetxt("{}/step{}/cal_MD/MD{}/MD_energy.txt".format(self.path, self.step, generate),
                       self.MD_energy, fmt='%f', delimiter=',')
        else:
            np.savetxt("{}/reduction_step{}/cal_MD/MD{}/MD_energy.txt".format(self.path, self.step, generate),
                       self.MD_energy, fmt='%f', delimiter=',')
        print('file saved at')
        if self.type == 'Oxidation':
            print("{}/step{}/cal_MD/MD{}/MD_energy.txt".format(self.path, self.step, generate))
        else:
            print("{}/reduction_step{}/cal_MD/MD{}/MD_energy.txt".format(self.path, self.step, generate))
        if generate >= 2:
            if self.type == 'Oxidation':
                energy1 = np.loadtxt("{}/step{}/cal_MD/MD{}/MD_energy.txt".format(self.path, self.step, generate - 1),
                                     delimiter=',')
                energy2 = np.loadtxt("{}/step{}/cal_MD/MD{}/MD_energy.txt".format(self.path, self.step, generate),
                                     delimiter=',')
            else:
                energy1 = np.loadtxt("{}/reduction_step{}/cal_MD/MD{}/MD_energy.txt".format(self.path, self.step,
                                                                                            generate - 1),
                                     delimiter=',')
                energy2 = np.loadtxt("{}/reduction_step{}/cal_MD/MD{}/MD_energy.txt".format(self.path, self.step,
                                                                                            generate), delimiter=',')
            # Bug fixed :需要类似归一化的处理方法才行，当X数量与Y数量级需要进行匹配时，需要根据Y数量进行归一化才行啊啊啊啊
            # 使用了对于X进行缩放的方法来处理，暂时看是可以停下来的
            mean_energy1 = np.mean(energy1)
            mean_energy2 = np.mean(energy2)
            slope1, intercept1 = get_slope(energy1)
            slope2, intercept2 = get_slope(energy2)
            print(slope1, intercept1)
            print(slope2, intercept2)
            print(abs(mean_energy1 - mean_energy2))
            if abs(slope1) <= 0.5 and abs(slope2) <= 0.5 and abs(mean_energy1 - mean_energy2) <= 0.5:
                return False
            else:
                if generate <= 20:
                    return True
                else:
                    return False
        else:
            return True

# # # # try
# path = '/data/home/max/Projects/OxideCuZn/NemethodCu3Zn/tryauto_safe'
# # # # path = '/data/home/max/Projects/OxideCuZn/NemethodCu3Zn/MDsecond/try'
# command = 'mpirun  /data/software/lasp/3.2.0/NN_pro3.2.0_intel18/Src/lasp > run.log 2>&1'
# laspcal = clculator_lasp(command, path, 0, 14, type='r')
# laspcal.run_energy()
# laspcal.get_remove_o_atoms()
# # laspcal.run_md()
# # # print(laspcal.site_list)
