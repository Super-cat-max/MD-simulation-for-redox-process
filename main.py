# -*- coding: utf-8 -*-
"""
@Time ： 2023/3/8 22:43
@Auth ： max
@File ：Read_energy.py
@IDE ：PyCharm
@Motto：ABC(Always Be Coding)
"""

from Gen_structures import *
from Read_energy import *
from add_O2_first import add_O2


def oxide_step(dic, environment, lasp_command):
    # 初始为input0.cif
    # 然后在表面寻找吸附位点
    # 依次计算吸附能，然后在吸附能小于0的位点加入O，进行MD模拟
    # 得到新的结构input1,在step1中进行新一轮计算
    # 终止条件为没有新的适合吸附的位点/step数量大于预计step数量
    test = False
    if test:
        step = 1  # Test
    else:
        step = 0
    # ____
    if step == 0:
        row = add_O2(dic, environment, lasp_command)
        row.run()
    else:
        pass

    if test:
        continue_Oxide = False  # Test
    else:
        continue_Oxide = True

    while continue_Oxide:
        make_folder = GenerateStructures('input{}.cif'.format(step), dic, environment, step)
        make_folder.make_energy_folders()
        laspcal = calculator_lasp(lasp_command, dic, step, make_folder.files + 1)
        laspcal.run_energy()
        make_folder.make_md_folders(laspcal.site_list)
        o_add = len(laspcal.site_list)
        o_remove = make_folder.remove_o_number
        # finished: we need a method to stop oxidation, I think it can base the Top layers
        judge_structure = read('{}/step{}/input{}.cif'.format(dic, step, step))
        ##
        if step >= 2:
            o_atoms = []
            judge_number = 0
            for atom in judge_structure:
                if atom.symbol == 'O':
                    o_atoms.append(atom)
            for o_atom in o_atoms:
                # print(o_atom.position[2])
                if o_atom.position[2] <= make_folder.first_layer + 0.75:
                    judge_number = judge_number + 1
            print(judge_number / len(o_atoms))
            if judge_number / len(o_atoms) >= 0.025:  # Todo: find a ratio
                continue_Oxide = False
                print(' 氧化物达到预定厚度')
                print('Oxide atoms reached the target location')
            else:
                pass
        else:
            pass
        ##
        if o_add == 0:
            continue_Oxide = False
            rate = 1
            print('cannot adsorb more O')
            print('Oxide process is end')
        else:
            rate = o_remove / o_add
        if rate >= 0.6:
            continue_Oxide = False
            print('cannot adsorb more O')
            print('Oxide process is end')
        else:
            pass

        tmp = 0
        while make_folder.continue_MD:
            tmp = tmp + 1
            laspcal.run_md(make_folder.MD_generate)

            make_folder.continue_MD = laspcal.read_MD_energy(make_folder.MD_generate)
            # # test
            # if tmp <= 2:  # test
            #     make_folder.continue_MD = False  # test

            # ___________
            make_folder.get_new_structures()

            #  _________________________________

        step = step + 1
        # if step == 2:  # Test
        #     continue_Oxide = False  # Test
        # _______________
        #    还原部分捏
        # _______________

    # _____________________reduction
    reduction_step = 0
    if test:
        step = 5  # Test
    else:
        pass
    # reduction_step = 4  # Test
    write('{}/Oxide.cif'.format(dic), read('{}/step{}/input{}.cif'.format(dic, step, step)))
    try:
        os.mkdir('{}/reduction_step{}'.format(dic, reduction_step))
    except OSError:
        pass
    write('{}/reduction_step{}/input{}.cif'.format(dic, reduction_step, reduction_step),
          read('{}/Oxide.cif'.format(dic)))

    continue_Reduction = True
    # continue_Reduction = False  # test

    for_final = 0
    while continue_Reduction:

        reduction_generator = GenerateStructures('input{}.cif'.format(reduction_step), dic, environment, reduction_step,
                                                 run_type='reduction')
        reduction_generator.make_reduction_folders()
        # print('miao', reduction_step, reduction_generator.files)
        las_reeducation_cal = calculator_lasp(lasp_command, dic, reduction_step, reduction_generator.files + 1,
                                              run_type='reduction')
        las_reeducation_cal.run_energy()
        reduction_generator.make_md_folders(las_reeducation_cal.site_list)

        o_remove = reduction_generator.remove_o_number
        print('o_remove', o_remove)
        # if o_remove == 0:
        #     continue_Reduction = False
        #     reduction_generator.continue_MD = False
        #     print('已经结束力！！！！')
        # else:
        #     continue_Reduction = True

        tmp_reduction = 0
        while reduction_generator.continue_MD:
            tmp_reduction = tmp_reduction + 1
            las_reeducation_cal.run_md(reduction_generator.MD_generate)
            try:
                reduction_generator.continue_MD = las_reeducation_cal.read_MD_energy(reduction_generator.MD_generate)
            except ValueError:
                reduction_generator.continue_MD = False
                print('Caught the Value error!')
                pass
            # test
            # if tmp_reduction <= 2:  # test
            #     reduction_generator.continue_MD = False  # test
            try:
                reduction_generator.get_new_structures()
            except FileNotFoundError:
                print('Caught the File not found error!')
                pass
        reduction_step = reduction_step + 1
        for_final = reduction_step
        if o_remove == 0:
            continue_Reduction = False
            reduction_generator.continue_MD = False
            for_final = reduction_step
            print('已经结束力！！！！')
        else:
            continue_Reduction = True
    # make sure all Oxygen atoms have been removed
    # reduction_step = 5  # Test
    # for_final = 5  # test
    print('The final structure is in step:{}'.format(for_final))
    try:
        os.mkdir('{}/reduction_structures'.format(dic))
    except OSError:
        pass
    # Bug fixed 2023 04 03:
    # In some times, the for_final will be bigger than reduction step, so we have to minus 1 to fix it!
    # don't ask me why, I'm just a cat, an Ecat in TJU-ECAT !
    try:
        search_final = read('{}/reduction_step{}/input{}.cif'.format(dic, for_final, for_final))
        print(search_final)
    except FileNotFoundError:
        for_final = for_final - 1
        search_final = read('{}/reduction_step{}/input{}.cif'.format(dic, for_final, for_final))
        print(search_final)
    try:
        final_structure = read('{}/reduction_step{}/input{}.cif'.format(dic, for_final, for_final))
        atoms_list = []
        for atoms in final_structure:
            atoms_list.append(atoms.symbol)
        try:
            if atoms_list.index('O') == -1:
                write('{}/final_structure.cif'.format(dic), final_structure)
            else:
                tmp_structure = copy.deepcopy(final_structure)

                atom_lists = np.array(atoms_list)
                del_list = []
                tmp_list = np.argwhere(atom_lists == 'O')
                for index in tmp_list:
                    del_list.append(int(index))
                del_list.sort(reverse=True)
                # print(del_list)
                for index in del_list:
                    del tmp_structure[index]

                write('{}/reduction_structures/input_final.cif'.format(dic), tmp_structure)
                write('{}/reduction_structures/input.arc'.format(dic), tmp_structure, format='dmol-arc')
                try:
                    os.system(
                        'cp {}/calculate_files/MD_settings/Reduction_long/lasp.in {}/reduction_structures/lasp.in'.format(dic, dic))
                except OSError:
                    pass
                try:
                    os.system('cp {}/calculate_files/Potential_files/{} {}/reduction_structures/{}'.format(dic,
                                                                                        environment, dic, environment))
                except OSError:
                    pass
                subprocess.call(lasp_command, shell=True, cwd='{}/reduction_structures'.format(dic))
                file = open('{}/reduction_structures/md.arc'.format(dic), 'r')
                line1 = file.readlines()
                file2 = open('{}/reduction_structures/md.arc'.format(dic), 'r')
                line2 = file2.readlines()
                stru = []
                for i in range(len(line1) - len(line2) + 4, len(line1) + 1):
                    tmp = linecache.getline('{}/reduction_structures/md.arc'.format(dic), i)
                    stru.append(tmp)
                f = open('{}/reduction_structures/new.arc'.format(dic), 'w')
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
                tmp = read('{}/reduction_structures/new.arc'.format(dic), format='dmol-arc')
                write('{}/reduction_structures/final_structure.arc'.format(dic), tmp, format='dmol-arc')
                write('{}/final_structure.cif'.format(dic), tmp)
        except ValueError:
            write('{}/final_structure.cif'.format(dic), final_structure)
    except FileNotFoundError:
        pass

    print('congratulations!')
    # reduction = GenerateStructures('input{}.cif'.format(step), dic, environment, step, )


lasp_command = 'mpirun  /data/software/lasp/3.2.0/NN_pro3.2.0_intel18/Src/lasp > run.log 2>&1'
path = '/data/home/max/Projects/OxideCuZn/NemethodCu3Zn/tryauto_safe'
oxide_step(path, 'CuZnCHO.pot', lasp_command)
