# -*- coding: utf-8 -*-
"""
@Time ： 2023/3/8 22:43
@Auth ： max
@File ：clean_finder.py
@IDE ：PyCharm
@Motto：ABC(Always Be Coding)

"""
import copy
import json
import os
import random
import sys
from collections import defaultdict

import numpy as np
from ase.build import add_adsorbate
from ase.build import molecule
from ase.constraints import FixAtoms
from ase.io import read
from math import cos
from math import sqrt
from pymatgen.io.vasp import Poscar
from scipy.spatial import Delaunay

sys.path.insert(0, os.path.dirname(os.getcwd()))
random.seed(3018)


class Adsorb_finder:
    """This class can be used to find the adsorb sites on any surface
        you need a poscar file as target structure
        we use alpha-sahpe method so that you can find sites in the hole of the surface
        you should give both the radius of adsorbate and adsorbent,the molecule
        you can adjust to adsorb angle from 0 to 90, the basic vector which detemined the direction
        and the number of the layer of atoms which have interaction with the molecule
        we will add more function in this class


        """

    def __init__(self, filename):
        self.sites_all = None
        self.ans = None
        self.nouse = None
        self.filename = filename
        self.file = self.extract_POSCAR()
        # In fact we need a way to change the D_a and D_b, But I am just a noob, I cannot do anything
        # So I use parameters which can almost fix to CuZnO-systems
        # I'm just a cat ₍˄·͈༝·͈˄₎◞ ̑̑ who can write code, you can't ask for too much! ( ´◔︎ ‸◔︎`)
        self.D_a = 1.5  # 0.86  # default adsorb 被吸附物半径
        self.D_b = 1.5  # 1.39  # default adsorbent 吸附质半径
        self.alpha = self.D_a + self.D_b  # default alpha
        self.judement = np.array([0, 0, 1])  # default basic vector
        self.default_judement = np.array([0, 0, 1])
        self.roundtol = 9  # default number of decimal for rounding
        self.cell = self.file['box_coord']
        self.xyz = self.file['coord_Cartesian']
        self.atoms = molecule('CO', vacuum=3.0)
        self.sub = 2  # default relation atoms
        self.angle = 30
        self.zhegaixishu = 0.5  # 遮盖系数，用来判断是否被遮盖 更新重点
        self.distance_coefficient = 1.25  # 偏移距离系数，是指安装放置吸附分子时的初始距离系数，为吸附质+吸附剂半径乘以该系数
        self.pbcs = np.array(
            [[0, 0, 0], [1, 0, 0], [-1, 0, 0], [0, 1, 0], [1, 1, 0], [-1, 1, 0], [0, -1, 0], [1, -1, 0], [-1, -1, 0],
             [0, 0, 1], [1, 0, 1], [-1, 0, 1], [0, 1, 1], [1, 1, 1], [-1, 1, 1], [0, -1, 1], [1, -1, 1], [-1, -1, 1],
             [0, 0, -1], [1, 0, -1], [-1, 0, -1], [0, 1, -1], [1, 1, -1], [-1, 1, -1], [0, -1, -1], [1, -1, -1],
             [-1, -1, -1]])
        xyz = np.array(self.xyz)
        cell = np.array(self.cell)
        xyz_frac = np.linalg.solve(cell.T, xyz.T).T
        xyz_frac = np.remainder(xyz_frac, 1)
        xyz_frac_pbcs = np.concatenate((self.pbcs.reshape(27, 1, 3) + xyz_frac.reshape(1, -1, 3)), axis=0)
        xyz_pbcs = np.dot(xyz_frac_pbcs, cell)
        self.alpha_array = xyz_pbcs
        self.system = self.extract_POSCAR()
        # self.atom_list = self.make_atomlist()#暂时用不到，该方法需要用到原有的POSCAR函数，但是原有的POSCAR函数不稳定，现在采用ASE代替，暂时没有完全代替
        self.TopXyz, self.BridgeXyz, self.HollowXyz, _ = self.FindSitesXYZ()
        self.database = None  # need a database of json format
        self.map_placetoput = None

    def extract_POSCAR(self):

        a = read(self.filename)
        system = {'coord_Cartesian': a.get_positions(),
                  'box_coord': [a.get_cell()[0], a.get_cell()[1], a.get_cell()[2]]}
        return system

    def make_atomlist(self):
        tag = 0
        tmp = []
        for i in self.system['atom_num']:
            for j in range(0, i):
                tmp.append(self.system['atom'][tag])
            tag = tag + 1
        return tmp

    def display(self):
        print(self.file)

    def alpha_shape_3D(self):
        """
        Compute the alpha shape (concave hull) of a set of 3D points.
        Parameters:
            pos - np.array of shape (n,3) points.
            alpha - alpha value.
            roundtol - number of decimal for rounding
        return
            outer surface vertex indices, edge indices, and triangle indices
        """

        tetra_vertices = Delaunay(self.alpha_array).vertices
        # Find radius of the circumsphere.
        # By definition, radius of the sphere fitting inside the tetrahedral needs
        # to be smaller than alpha value
        # http://mathworld.wolfram.com/Circumsphere.html
        tetrapos = np.take(self.alpha_array, tetra_vertices, axis=0)
        normsq = np.sum(tetrapos ** 2, axis=2)[:, :, None]
        ones = np.ones((tetrapos.shape[0], tetrapos.shape[1], 1))
        a = np.linalg.det(np.concatenate((tetrapos, ones), axis=2))
        Dx = np.linalg.det(np.concatenate((normsq, tetrapos[:, :, [1, 2]], ones), axis=2))
        Dy = -np.linalg.det(np.concatenate((normsq, tetrapos[:, :, [0, 2]], ones), axis=2))
        Dz = np.linalg.det(np.concatenate((normsq, tetrapos[:, :, [0, 1]], ones), axis=2))
        c = np.linalg.det(np.concatenate((normsq, tetrapos), axis=2))
        # Remove bad tetrahedrals. These the ones where volume is zero.
        bad = a == 0
        num = Dx ** 2 + Dy ** 2 + Dz ** 2 - 4 * a * c
        bad[num < 0] = True
        bad = np.where(bad)[0]
        tetra_vertices = np.delete(tetra_vertices, bad, axis=0)
        num = np.delete(num, bad, axis=0)
        a = np.delete(a, bad, axis=0)
        # get radius
        r = np.sqrt(num) / (2 * np.abs(a))
        # Find tetrahedrals
        tetras = tetra_vertices[r < self.alpha, :]
        # triangles
        TriComb = np.array([(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)])
        Triangles = tetras[:, TriComb].reshape(-1, 3)
        Triangles = np.sort(Triangles, axis=1)
        # Remove triangles that occurs twice, because they are within shapes
        TrianglesDict = defaultdict(int)
        for tri in Triangles: TrianglesDict[tuple(tri)] += 1
        Triangles = np.array([tri for tri in TrianglesDict if TrianglesDict[tri] == 1])
        # edges
        if Triangles.size == 0:
            return [], [], []
        EdgeComb = np.array([(0, 1), (0, 2), (1, 2)])
        Edges = Triangles[:, EdgeComb].reshape(-1, 2)
        Edges = np.sort(Edges, axis=1)
        Edges = np.unique(Edges, axis=0)
        Vertices = np.unique(Edges)
        return Vertices, Edges, Triangles

    def FindSitesXYZ(self):  # Find site positions in all directions
        xyz = np.array(self.xyz)
        cell = np.array(self.cell)
        xyz_frac = np.linalg.solve(cell.T, xyz.T).T
        xyz_frac = np.remainder(xyz_frac, 1)
        # account for the pbc condition before applying alpha shape
        xyz_frac_pbcs = np.concatenate((self.pbcs.reshape(27, 1, 3) + xyz_frac.reshape(1, -1, 3)), axis=0)
        xyz_pbcs = np.dot(xyz_frac_pbcs, cell)
        self.alpha_array = xyz_pbcs
        V, E, T = self.alpha_shape_3D()
        if isinstance(V, list) or isinstance(V, np.ndarray) and V.size == 0:
            return [], [], [], 0.0
        # Find ones that are within the original cell
        V = V[V < len(xyz)]
        E = E[np.any(E < len(xyz), axis=1)]
        T = T[np.any(T < len(xyz), axis=1)]
        if V.size == 0:
            return [], [], [], 0.0
        # Remove duplicates
        E = E[np.unique(np.sort(np.remainder(E, len(xyz)), axis=1), axis=0, return_index=True)[1], :]
        T = T[np.unique(np.sort(np.remainder(T, len(xyz)), axis=1), axis=0, return_index=True)[1], :]

        TopXyz = xyz_pbcs[V, :]
        BridgeXyz = np.mean(xyz_pbcs[E, :], axis=1)
        HollowXyz = np.mean(xyz_pbcs[T, :], axis=1)

        temp1 = xyz_pbcs[T[:, [[0, 1], [1, 2], [0, 2]]], :]
        triangle_lengths = np.linalg.norm(temp1[:, :, 0, :] - temp1[:, :, 1, :], axis=2)
        ph = np.sum(triangle_lengths, axis=1) / 2
        area = np.sqrt(
            ph * (ph - triangle_lengths[:, 0]) * (ph - triangle_lengths[:, 1]) * (ph - triangle_lengths[:, 2]))
        # self.TopXyz = TopXyz
        # self.BridgeXyz = BridgeXyz
        # self.HollowXyz = HollowXyz
        return TopXyz, BridgeXyz, HollowXyz, np.sum(area)

    def distance(self, p1, p2, m):
        if m == 2:
            d = sqrt((p1[0] - p2[0]) * (p1[0] - p2[0]) + (p1[1] - p2[1]) * (p1[1] - p2[1]))
        elif m == 3:
            d = sqrt((p1[0] - p2[0]) * (p1[0] - p2[0]) + (p1[1] - p2[1]) * (p1[1] - p2[1]) + (p1[2] - p2[2]) * (
                    p1[2] - p2[2]))
        return d

    def V_minus(self, v1, v2):
        '''向量由V1指向V2'''  # 对于此项目，应该是吸附原子在V2处
        if len(v1) != len(v2):
            print("数据错误")
        else:
            v0 = []
            for i in range(0, len(v1)):
                v0.append(v2[i] - v1[i])
            return v0

    def find_sites(self, TBH=0, Tset=None, Bset=None, Hset=None):
        # 在构造函数中调用了FindesitesXYZ 于是这个部分不需要了应该是这样对的一定是这样的(确信
        # self.FindSitesXYZ()
        # 纯体积版本
        #     surface = []
        ''' TBH = 1 means Top
            TBH = 2 means Bridge
            TBH = 3 means Hollow
            defalt or other input means output all sites
            '''
        if TBH == 1:
            self.sites_all = self.TopXyz
        elif TBH == 2:
            self.sites_all = self.BridgeXyz
        elif TBH == 3:
            self.sites_all = self.HollowXyz
        else:
            tmp = []
            if Tset != None:
                random.shuffle(self.TopXyz)
                if len(self.TopXyz) < Tset:
                    print('warning! the Top site is less than {} !'.format(Tset))
                for iR in range(0, Tset):
                    tmp.append(self.TopXyz[iR])
            else:
                for i in self.TopXyz:
                    tmp.append(i)
            if Bset != None:
                random.shuffle(self.BridgeXyz)
                if len(self.BridgeXyz) < Bset:
                    print('warning! the Bridge site is less than {} !'.format(Bset))
                for jR in range(0, Bset):
                    tmp.append(self.BridgeXyz[jR])
            else:
                for j in self.BridgeXyz:
                    tmp.append(j)
            if Hset != None:
                random.shuffle(self.HollowXyz)
                if len(self.HollowXyz) < Hset:
                    print('warning! the Hollow site is less than {} !'.format(Hset))
                for kR in range(0, Hset):
                    tmp.append(self.HollowXyz[kR])
            else:
                for k in self.HollowXyz:
                    tmp.append(k)
            self.sites_all = tmp
        nouse = []
        ans = []
        # 此处出现重复现象，原因我不知道，解决方法采用np自带的去重复方法处理
        # 该步骤为了计算空间向量，以便于添加吸附原子
        # 该部分用来去除无法吸附的位点
        # 为3.5工作的最终落实
        remove_over = []

        for i in self.sites_all:
            over = 0
            for j in self.xyz:
                if self.distance(i, j, 2) <= (self.D_a + self.D_b) * self.zhegaixishu and i[2] + 0.2 < j[2]:
                    # print(self.distance(i, j, 2))
                    # 出现遮挡同时
                    # print('zhedang')
                    over = over + 1
                else:
                    pass
            # print(over)
            if over == 0:
                remove_over.append(i)

            else:
                over = 0
                # print(self.distance(i, j, 2))
                pass
        # 此处筛选掉了所有存在遮挡的球体

        if self.judement.all() != self.default_judement.all():
            pass
        else:
            self.sites_all = remove_over
        # 此部分寻找空间向量
        for i in self.sites_all:
            tmp = []  # print(tmp)#寻找近邻原子，即有空间接触的原子
            for j in self.xyz:
                if self.distance(i, j, 3) <= self.sub * (self.D_a + self.D_b):
                    tmp.append(j)  # 此处的tmp为近邻原子  注：不考虑晶格空间重复性，因为是整个表面
                else:
                    pass
            tmVector = []  # 计算空间向量
            remakelength = []
            for k in tmp:
                tmpv = self.V_minus(k, i)
                tmVector.append(tmpv)
                remakelength.append(self.distance(k, i, 3))
            # 此处获得所有空间向量，由吸附质指向吸附剂
            if len(tmVector) == 0:
                pass  # 一会再加上
            # 认为这不是个好东西，直接Pass掉，不加了
            else:
                # 20220428修改如下：引入了公式rnew=r*(Da+Db)/distance.用来解决距离远的原子贡献的空间向量反而大的问题
                ntmp = np.array(tmVector)  # 转换为numpy
                for remake_len in range(0, len(ntmp)):
                    # print(ntmp[remake_len],remakelength[remake_len])
                    if remakelength[remake_len] * remakelength[remake_len] <= 0.05:  # 识别进0值防止出现除法出现错误
                        pass
                    else:
                        ntmp[remake_len] = ntmp[remake_len] * self.distance_coefficient * (self.D_a + self.D_b) / \
                                           remakelength[remake_len]
                # print('空间向量',len(ntmp),remakelength)
                V_sum = sum(k1 for k1 in ntmp)
                V_sum = V_sum / sqrt(sum(V_sum * V_sum))  # 得到加和后的单位向量，标准化再此处实现
                if np.dot(self.judement, V_sum) >= cos(self.angle):  # 夹角如果是正数，说明原子在上表面
                    for j2 in self.xyz:
                        if self.distance(i, j2, 3) <= 1.12 * (self.D_a + self.D_b):  # 判断一下原子是否在范围内
                            ans.append([i, V_sum])
                        else:
                            nouse.append(i)
                else:
                    nouse.append(i)
        nouse = np.unique(nouse, axis=0)
        ans = np.array(ans)
        ans = np.unique(ans, axis=0)
        self.ans = ans
        self.nouse = nouse
        return ans, nouse

    def find_idx(self):
        '''用来匹配位点与原子类别'''
        positions = self.xyz
        Top_idx = []
        T, _ = self.find_sites(1)
        # print(T)
        for iT in T:
            T_distance = []
            # T_disdict = {}
            for jT in positions:
                d = self.distance(iT[0], jT, 3)
                T_distance.append(d)
            Top_idx.append(T_distance.index(min(T_distance)) + 1)
        re_T = []
        for i in range(0, len(T)):
            re_T.append([T[i], Top_idx[i]])
        # print(re_T)
        Bridge_idx = []
        B, _ = self.find_sites(2)
        for iB in B:
            B_distance = []
            # T_disdict = {}
            for jB in positions:
                d = self.distance(iB[0], jB, 3)
                B_distance.append(d)
            # print(B_distance)
            tmp = sorted(B_distance)
            # print(tmp)
            # print([B_distance.index(tmp[0]),B_distance.index(tmp[1])])
            Bridge_idx.append([B_distance.index(tmp[0]) + 1, B_distance.index(tmp[1]) + 1])
            # 此处+1是为了与POSCAR文件贴合
        re_B = []
        for i in range(0, len(B)):
            re_B.append([B[i], Bridge_idx[i]])
        # print(re_T)
        Hollow_idx = []
        H, _ = self.find_sites(3)
        for iH in H:
            H_distance = []
            # T_disdict = {}
            for jH in positions:
                d = self.distance(iH[0], jH, 3)
                H_distance.append(d)
            tmp1 = sorted(H_distance)
            Hollow_idx.append(
                [H_distance.index(tmp1[0]) + 1, H_distance.index(tmp1[1]) + 1, H_distance.index(tmp1[2]) + 1])
            # print(H_distance.index(min(H_distance)))
        re_H = []
        for i in range(0, len(T)):
            re_H.append([H[i], Hollow_idx[i]])
        # print(re_T)
        return re_T, re_B, re_H

    def find_moleculeplace(self, Top):
        # self.find_sites()
        tmp = []
        if Top != None:
            for i in self.ans:
                tmp.append([i[0][0], i[0][1], i[0][2] + self.distance_coefficient * (self.D_a + self.D_b)])
        else:
            for i in self.ans:
                tmp.append(i[0] + self.distance_coefficient * (self.D_a + self.D_b) * i[1])
        self.put_place = tmp
        self.map_placetoput = []
        for i in range(0, len(self.put_place)):
            self.map_placetoput.append([self.ans[i], self.put_place[i]])

    def change_molecule(self, moleculename):
        self.atoms = molecule(moleculename, vacuum=3.0)
        self.mol = moleculename

    def get_knowledge(self):
        database = json.read(self.database)
        '''从self.'''

    def add_molecule(self, path, filetype=0, fixed=False, Top=None):
        self.find_moleculeplace(Top)
        tmp = []
        for i in self.xyz:
            tmp.append(i[2])
        height_correct = max(tmp)
        atoms = self.atoms

        j = 0

        # 分别寻找位点
        for i in self.map_placetoput:
            structure_row = read(self.filename)  # ase.io.read
            # get Z "correct the function
            # we can move it to the __init__ to speed up this program
            info = structure_row.info.get('adsorbate_info', {})
            if 'top layer atom index' in info:
                a = info['top layer atom index']
            else:
                a = structure_row.positions[:, 2].argmax()
                if 'adsorbate_info' not in structure_row.info:
                    structure_row.info['adsorbate_info'] = {}
                structure_row.info['adsorbate_info']['top layer atom index'] = a
                # print(type(structure_row))
            z = structure_row.positions[a, 2]
            d = i[1][0:2]
            h = i[1][2] - z
            j = j + 1
            add_adsorbate(structure_row, atoms, h, d, mol_index=0)
            # print(structure_row.positions)
            # print(i)
            if fixed == True:
                distance = []
                for k in structure_row.get_positions():
                    # print(k)
                    dtmp = self.distance(k, i[1], 3)
                    distance.append(dtmp)
                # print(type(structure_row))
                # print(distance)
                aaaaa = np.array(distance)
                index = np.argsort(aaaaa)
                # print(index)
                structure_row = structure_row[index]
                #
                # distance1 = []
                # for l in structure_row.get_positions():
                #     # print(k)
                #     dtmp1 = self.distance(l, i[1], 3)
                #     distance1.append(dtmp1)
                # aaaaa2 = np.array(distance1)
                # index2 = np.argsort(aaaaa2)
                # print(index2)
                # structure_row = structure_row[-index]
                # for k in structure_row.get_positions():
                #     dtmp = self.distance(k, i[1], 3)
                #     distance1.append(dtmp)
                # print(index)
                # print(np.argsort(distance1))
                # print(distance1)
                # stu = ase.io.read(os.path.join(root, name))

                c = FixAtoms(indices=[atom.index for atom in structure_row if
                                      self.distance(atom.position, i[0][0], 3) >= 2 * self.distance_coefficient * (
                                              self.D_a + self.D_b)])
                structure_row.set_constraint(c)

                # structure_row.write('./Tryfit/{}_fixed'.format(name))

            if filetype == 0:
                structure_row.write('{}/POSCAR_{}'.format(path, j))
            else:
                pass
            #     structure_row.write('{}/{}on{}_{}.cif'.format(path, self.mol, self.filename, j))
        # title = './{}/all{}on{}angle{}sub{}d{}V{}'.format(path, self.mol, self.filename, self.angle, self.sub,
        #                                                   self.distance_coefficient, self.judement)


# 与LS-CGCNN对接部分，详情可见jupyter文件

class Make_dataset(Adsorb_finder):
    def FindSitesIdx(self):
        # 相当于找出所有点和每一个点相邻的所有边、三角
        # 号码是从0开始的，跟POSCAR命名顺序有别，注意此点
        xyz = np.array(self.xyz)
        cell = np.array(self.cell)
        xyz_frac = np.linalg.solve(cell.T, xyz.T).T
        xyz_frac = np.remainder(xyz_frac, 1)
        # account for the pbc condition before applying alpha shape
        xyz_frac_pbcs = np.concatenate((self.pbcs.reshape(27, 1, 3) + xyz_frac.reshape(1, -1, 3)), axis=0)
        xyz_pbcs = np.dot(xyz_frac_pbcs, cell)
        V, E, T = self.alpha_shape_3D()
        if isinstance(V, list) or isinstance(V, np.ndarray) and V.size == 0:
            return [], [], [], 0.0
        # Find unique
        V = np.unique(np.remainder(V, len(xyz))).reshape(-1, 1)
        E = np.unique(np.sort(np.remainder(E, len(xyz)), axis=1), axis=0)
        T = np.unique(np.sort(np.remainder(T, len(xyz)), axis=1), axis=0)
        if V.size == 0:
            return [], [], [], 0.0
        self.V = V
        self.E = E
        self.T = T
        return V, E, T

    def set_formation(self):  #
        self.FindSitesIdx()
        Newdata = []
        for site_idx in self.V.tolist() + self.E.tolist() + self.T.tolist():
            tmp_dict = {'cell': self.cell, 'positions': self.xyz, 'pbc': True}
            newdatum = copy.deepcopy(tmp_dict)
        newdatum['site_idx'] = site_idx
        Newdata.append(newdatum)
        json.dump(Newdata, open('./labeled_data{}.json'.format(self.filename), 'w'))


class Make_structures(Make_dataset):
    def __init__(self, filename):
        Adsorb_finder.__init__(self, filename)
        # self.filename = filename
        tmp, _ = self.find_sites()
        self.rowpositions = []
        self.positionvector = []
        for i in tmp:
            self.rowpositions.append(i[0])
            self.positionvector.append(i[1])
        atoms = molecule('CO', vacuum=3.0)
        self.atoms = atoms
        self.path = './structure/{}'
        height_list = []
        for i in self.xyz:
            height_list.append(i[-1])
        self.height_correct = max(height_list)

    def change_molecule(self, moleculename):
        self.atoms = molecule(moleculename, vacuum=3.0)

    def write_structure(self):
        repair_top = []
        for i in self.xyz:
            repair_top.append(i[2])
        self.z_repair = max(repair_top)
        structure_row = read(self.filename)
        poscar = Poscar.from_file(self.filename)

        j = 0

        for i in self.rowpositions:
            structure_row = read(self.filename)
            # get Z "correct the function
            # we can move it to the __init__ to speed up this program
            info = structure_row.info.get('adsorbate_info', {})
            if 'top layer atom index' in info:
                a = info['top layer atom index']
            else:
                a = structure_row.positions[:, 2].argmax()
                if 'adsorbate_info' not in structure_row.info:
                    structure_row.info['adsorbate_info'] = {}
                structure_row.info['adsorbate_info']['top layer atom index'] = a
            z = structure_row.positions[a, 2]

            d = i[0:2]
            h = i[2] - z
            j = j + 1
            add_adsorbate(structure_row, self.atoms, h, d, mol_index=0)

            structure_row.write(self.path + '/POSCAR_{}'.format(j))

    def add_structure(self):
        repair_top = []
        for i in self.xyz:
            repair_top.append(i[2])
        self.z_repair = max(repair_top)
        structure_row = read(self.filename)
        # get Z "correct the function
        # we can move it to the __init__ to speed up this program
        info = structure_row.info.get('adsorbate_info', {})
        if 'top layer atom index' in info:
            a = info['top layer atom index']
        else:
            a = structure_row.positions[:, 2].argmax()
            if 'adsorbate_info' not in structure_row.info:
                structure_row.info['adsorbate_info'] = {}
            structure_row.info['adsorbate_info']['top layer atom index'] = a
        z = structure_row.positions[a, 2]
        poscar = Poscar.from_file(self.filename)
        # atoms = molecule('CO', vacuum=3.0)
        j = 0
        self.rowpositions1 = np.array(self.rowpositions) + 1.15 * self.alpha * np.array(self.positionvector)
        for i in self.rowpositions1:
            structure_row = read(self.filename)
            d = i[0:2]
            h = i[2] - z
            j = j + 1
            add_adsorbate(structure_row, self.atoms, h, d, mol_index=0)
            # structure_row.write('./structure/POSCAR_{}'.format(j))
            structure_row.write(self.path + '/POSCAR_{}'.format(j))


class Adsorb_finderPro(Adsorb_finder):
    def __init__(self, filename):
        Adsorb_finder.__init__(self, filename)
        self.movecoefficient = 0.5
        self.focus = self.find_focus()

    def After_settings(self):
        self.screenedTop, _ = self.find_sites(TBH=1)
        self.screenedBridge, _ = self.find_sites(TBH=2)
        self.screenedHollow, _ = self.find_sites(TBH=3)
        # # self.Topplace, self.Bridgeplace, self.Hollowplace = self.Get_position()
        self.re_T, self.re_B, self.re_H = self.find_idxRe()

        self.Topplace, self.Bridgeplace, self.Hollowplace = self.Get_position()

    def find_focus(self):
        find_f = self.xyz
        return sum(find_f) / len(find_f)

    def find_idxRe(self):
        '''用来匹配位点与原子类别'''
        positions = self.TopXyz
        Top_idx = []
        T = self.screenedTop
        for iT in T:
            T_distance = []
            # T_disdict = {}
            for jT in positions:
                d = self.distance(iT[0], jT, 3)
                T_distance.append(d)
            tmp = sorted(T_distance)
            Top_idx.append([T_distance.index(tmp[1]), T_distance.index(tmp[2]), T_distance.index(tmp[3])])
        re_T = []
        for i in range(0, len(T)):
            re_T.append([T[i], Top_idx[i]])
        # print(re_T)
        Bridge_idx = []
        B = self.screenedBridge
        for iB in B:
            B_distance = []
            # T_disdict = {}
            for jB in positions:
                d = self.distance(iB[0], jB, 3)
                B_distance.append(d)
            # print(B_distance)
            tmp = sorted(B_distance)
            # print(tmp)
            # print([B_distance.index(tmp[0]),B_distance.index(tmp[1])])
            Bridge_idx.append(
                [B_distance.index(tmp[0]), B_distance.index(tmp[1]), B_distance.index(tmp[2]), B_distance.index(tmp[3]),
                 B_distance.index(tmp[4])])
            # 此处+1是为了与POSCAR文件贴合
        re_B = []
        for i in range(0, len(B)):
            re_B.append([B[i], Bridge_idx[i]])
        # print(re_T)
        Hollow_idx = []
        H = self.screenedHollow
        for iH in H:
            H_distance = []
            # T_disdict = {}
            for jH in positions:
                d = self.distance(iH[0], jH, 3)
                H_distance.append(d)
            tmp1 = sorted(H_distance)
            Hollow_idx.append(
                [H_distance.index(tmp1[0]), H_distance.index(tmp1[1]), H_distance.index(tmp1[2])])
            # print(H_distance.index(min(H_distance)))
        re_H = []
        for i in range(0, len(H)):
            re_H.append([H[i], Hollow_idx[i]])
        # print(re_T)

        return re_T, re_B, re_H

    def Get_position(self):
        Topplace = []
        Bridgeplace = []
        Hollowplace = []
        for i in self.re_T:
            # print(i[0][0])
            # print(self.xyz_surface[i[1][0]], self.xyz_surface[i[1][1]], self.xyz_surface[i[1][2]])
            # print(i[0][0])
            # print(a.TopXyz[i[1][0]], a.TopXyz[i[1][1]], a.TopXyz[i[1][2]])
            # print(i[0][0])
            vector_a = self.TopXyz[i[1][1]] - self.TopXyz[i[1][0]]
            vector_b = self.TopXyz[i[1][2]] - self.TopXyz[i[1][0]]

            v_vx = vector_a[1] * vector_b[2] - vector_a[2] * vector_b[1]
            v_vy = vector_a[2] * vector_b[0] - vector_a[0] * vector_b[2]
            v_vz = vector_a[0] * vector_b[1] - vector_a[1] * vector_b[0]
            vector_vertical = [v_vx / sqrt(v_vx * v_vx + v_vy * v_vy + v_vz * v_vz),
                               v_vy / sqrt(v_vx * v_vx + v_vy * v_vy + v_vz * v_vz),
                               v_vz / sqrt(v_vx * v_vx + v_vy * v_vy + v_vz * v_vz)]
            if vector_vertical[2] > 0:
                pass
            else:
                vector_vertical[0] = -vector_vertical[0]
                vector_vertical[1] = -vector_vertical[1]
                vector_vertical[2] = -vector_vertical[2]
            a = [i[0][0][0] - (self.D_b + self.D_a) * vector_vertical[0] * self.movecoefficient,
                 i[0][0][1] - (self.D_b + self.D_a) * vector_vertical[1] * self.movecoefficient,
                 i[0][0][2] - (self.D_b + self.D_a) * vector_vertical[2] * self.movecoefficient]
            b = [i[0][0][0] + (self.D_b + self.D_a) * vector_vertical[0] * self.movecoefficient,
                 i[0][0][1] + (self.D_b + self.D_a) * vector_vertical[1] * self.movecoefficient,
                 i[0][0][2] + (self.D_b + self.D_a) * vector_vertical[2] * self.movecoefficient]
            a_d = self.distance(a, self.focus, 3)
            b_d = self.distance(b, self.focus, 3)
            # print(a_d,b_d)
            if a_d >= b_d:
                Topplace.append(a)
            else:
                Topplace.append(b)

        for j in self.re_B:
            vector_a = self.TopXyz[j[1][2]] - self.TopXyz[j[1][1]]
            vector_b = self.TopXyz[j[1][3]] - self.TopXyz[j[1][1]]
            v_vx = vector_a[1] * vector_b[2] + vector_a[2] * vector_b[1]
            v_vy = vector_a[2] * vector_b[0] + vector_a[0] * vector_b[2]
            v_vz = vector_a[0] * vector_b[1] + vector_a[1] * vector_b[0]
            vector_vertical = [v_vx / sqrt(v_vx * v_vx + v_vy * v_vy + v_vz * v_vz),
                               v_vy / sqrt(v_vx * v_vx + v_vy * v_vy + v_vz * v_vz),
                               v_vz / sqrt(v_vx * v_vx + v_vy * v_vy + v_vz * v_vz)]
            a = [j[0][0][0] - (self.D_b + self.D_a) * vector_vertical[0] * self.movecoefficient,
                 j[0][0][1] - (self.D_b + self.D_a) * vector_vertical[1] * self.movecoefficient,
                 j[0][0][2] - (self.D_b + self.D_a) * vector_vertical[2] * self.movecoefficient]
            b = [j[0][0][0] + (self.D_b + self.D_a) * vector_vertical[0] * self.movecoefficient,
                 j[0][0][1] + (self.D_b + self.D_a) * vector_vertical[1] * self.movecoefficient,
                 j[0][0][2] + (self.D_b + self.D_a) * vector_vertical[2] * self.movecoefficient]
            a_d = self.distance(a, self.focus, 3)
            b_d = self.distance(b, self.focus, 3)
            if a_d >= b_d:
                Bridgeplace.append(a)
            else:
                Bridgeplace.append(b)
            # 加到中间
            # Bridgeplace.append(j)
        for k in self.re_H:
            # k[0][0][2] = k[0][0][2] + (self.D_a + self.D_b) * self.distance_coefficient/1.414
            # print(k[0][0])
            vector_a = self.TopXyz[k[1][1]] - self.TopXyz[k[1][0]]
            vector_b = self.TopXyz[k[1][2]] - self.TopXyz[k[1][0]]
            v_vx = vector_a[1] * vector_b[2] + vector_a[2] * vector_b[1]
            v_vy = vector_a[2] * vector_b[0] + vector_a[0] * vector_b[2]
            v_vz = vector_a[0] * vector_b[1] + vector_a[1] * vector_b[0]
            vector_vertical = [v_vx / sqrt(v_vx * v_vx + v_vy * v_vy + v_vz * v_vz),
                               v_vy / sqrt(v_vx * v_vx + v_vy * v_vy + v_vz * v_vz),
                               v_vz / sqrt(v_vx * v_vx + v_vy * v_vy + v_vz * v_vz)]
            a = [k[0][0][0] - (self.D_b + self.D_a) * vector_vertical[0] * self.movecoefficient,
                 k[0][0][1] - (self.D_b + self.D_a) * vector_vertical[1] * self.movecoefficient,
                 k[0][0][2] - (self.D_b + self.D_a) * vector_vertical[2] * self.movecoefficient]
            b = [k[0][0][0] + (self.D_b + self.D_a) * vector_vertical[0] * self.movecoefficient,
                 k[0][0][1] + (self.D_b + self.D_a) * vector_vertical[1] * self.movecoefficient,
                 k[0][0][2] + (self.D_b + self.D_a) * vector_vertical[2] * self.movecoefficient]
            a_d = self.distance(a, self.focus, 3)
            b_d = self.distance(b, self.focus, 3)
            if a_d >= b_d:
                Hollowplace.append(a)
            else:
                Hollowplace.append(b)
            # Hollowplace.append(k)
            # print(k[0][0])
        return Topplace, Bridgeplace, Hollowplace

    def add_moleculeRe(self, path, fixed=False):
        self.Topplace, self.Bridgeplace, self.Hollowplace = self.Get_position()
        # self.Topplace, self.Bridgeplace, self.Hollowplace
        tmp = []
        for i in self.xyz:
            tmp.append(i[2])
        height_correct = max(tmp)
        atoms = self.atoms
        j = 0
        # 分别寻找位点
        type = ['Top', 'Bridge', 'Hollow']
        for strtype in type:
            target = []
            if strtype == 'Top':
                print('Top')
                for t in self.Topplace:
                    target.append(t[0][0])
                folder = 'str1'
                print(target)
            elif strtype == 'Bridge':
                print('Bridge')
                for y in self.Bridgeplace:
                    target.append(y[0][0])

                folder = 'str2'
            else:
                print('Hollow')
                for u in self.Hollowplace:
                    target.append(u[0][0])
                # target = self.Hollowplace
                folder = 'str3'
            for i in target:
                # print(i[0][0])
                structure_row = read(self.filename)  # ase.io.read
                # get Z "correct the function
                # we can move it to the __init__ to speed up this program
                info = structure_row.info.get('adsorbate_info', {})
                if 'top layer atom index' in info:
                    a = info['top layer atom index']
                else:
                    a = structure_row.positions[:, 2].argmax()
                    if 'adsorbate_info' not in structure_row.info:
                        structure_row.info['adsorbate_info'] = {}
                    structure_row.info['adsorbate_info']['top layer atom index'] = a
                z = structure_row.positions[a, 2]

# a = Adsorb_finderPro('unrelaxed.cif')
# a.angle = 90
# print(a.D_a+a.D_b)
# a.After_settings()
# top = a.re_T
# T,B,H,_ = a.FindSitesXYZ()
# # print(a.re_B)
# for i in H:
#     print(i)
