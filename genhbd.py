# 导入 PyMOL 的 Python 模块
from pymol import cmd

class HdbLine():
    res_name = "MOL"
    def __init__(self,h_type=-1, h_num=-1, h_name=" " , i_atom=" ",j_atom=" ",k_atom=" ", l_atom=" "):
        self.h_type = h_type
        self.h_num  = h_num 
        self.h_name = h_name
        self.i_atom = i_atom
        self.j_atom = j_atom
        self.k_atom = k_atom
        self.l_atom = l_atom
    
    def gen_hdb(self):
        return "{:8s}{:8s}{:8s}{:8s}{:8s}{:8s}{:8s}".format(str(self.h_num), str(self.h_type), self.h_name,self.i_atom, self.j_atom, self.k_atom, self.l_atom)

# cmd.load("cgu_h.pdb")
# weight_atoms_sele = cmd.select("weight_atoms", "resname TTT and not elem H")
# 使用 PyMOL 的 Python API 获取与目标原子相连的氢原子信息
# bond_atoms = cmd.get_model(f'bound_to {my_object}').atom
weight_atoms_sele = cmd.select("weight_atoms", "sele and not elem H")
for atom_index in cmd.index(selection="(weight_atoms)"):
    # print(atom_index)
    single_atom_info = cmd.get_model(f"index {atom_index[1]}").atom[0]
    single_atom = cmd.select("single_atom", f"index {atom_index[1]}")
    
    # 获取bound_to的原子信息
    bound_to_info = cmd.get_model(f"bound_to single_atom").atom
    bound_to_info_sele = cmd.select("bound_to_info_sele",f"bound_to single_atom")
    h_atoms = [ atom for atom in bound_to_info if atom.name[0]=="H"]
    wg_atoms = [ atom for atom in bound_to_info if atom.name[0]!="H"]
    h_nums = len(h_atoms)
    wg_nums = len(wg_atoms)
    bond_type = -1
    h_name = " "
    i_atom = " "
    j_atom = " "
    k_atom = " "
    l_atom = " "
    if h_nums == 2:
        if wg_nums>=2:
            #ok
            bond_type = 6
            i_atom = single_atom_info.name
            j_atom = wg_atoms[0].name
            k_atom = wg_atoms[1].name
            single_hbd = HdbLine(bond_type, h_nums, h_atoms[0].name[0:2],i_atom, j_atom, k_atom, l_atom)
            # print(f"{h_nums} {bond_type} {h_atoms[0].name[0:2]} {i_atom} {j_atom} {k_atom}")
            print(single_hbd.gen_hdb())
        if wg_nums == 1:
            bond_type = 3
            i_atom = single_atom_info.name
            j_atom = wg_atoms[0].name
            temp_single_atom = cmd.select("temp_single_atom", f"name {j_atom} and chain {wg_atoms[0].chain} and resname {wg_atoms[0].resn} and bound_to_info_sele")
            temp_atoms = cmd.get_model("bound_to temp_single_atom").atom
            temp_wg_atoms = [ atom for atom in temp_atoms if atom.name[0]!="H"]
            k_atom = temp_wg_atoms[0].name
            single_hbd = HdbLine(bond_type, h_nums, h_atoms[0].name[0:2],i_atom, j_atom, k_atom, l_atom)
            print(single_hbd.gen_hdb())
    if h_nums == 1 :
        if wg_nums == 2:
            #ok
            bond_type = 1
            i_atom = single_atom_info.name
            j_atom = wg_atoms[0].name
            k_atom = wg_atoms[1].name
            single_hbd = HdbLine(bond_type, h_nums, h_atoms[0].name,i_atom, j_atom, k_atom, l_atom)
            # print(f"{h_nums} {bond_type} {h_atoms[0].name[0:2]} {i_atom} {j_atom} {k_atom}")
            print(single_hbd.gen_hdb())
        if wg_nums >= 3:
            #ok
            bond_type = 5
            i_atom = single_atom_info.name
            j_atom = wg_atoms[0].name
            k_atom = wg_atoms[1].name
            l_atom = wg_atoms[2].name
            single_hbd = HdbLine(bond_type, h_nums, h_atoms[0].name,i_atom, j_atom, k_atom, l_atom)
            print(single_hbd.gen_hdb())
        if wg_nums == 1:
            #ok
            bond_type = 2
            i_atom = single_atom_info.name
            j_atom = wg_atoms[0].name
            temp_single_atom = cmd.select("temp_single_atom", f"name {j_atom} and chain {wg_atoms[0].chain} and resname {wg_atoms[0].resn} and bound_to_info_sele")
            temp_atoms = cmd.get_model("bound_to temp_single_atom").atom
            temp_wg_atoms = [ atom for atom in temp_atoms if atom.name[0]!="H"]
            k_atom = temp_wg_atoms[0].name
            single_hbd = HdbLine(bond_type, h_nums, h_atoms[0].name[0:2],i_atom, j_atom, k_atom, l_atom)
            print(single_hbd.gen_hdb())
    if h_nums == 3:
        i_atom = single_atom_info.name
        j_atom = wg_atoms[0].name
        temp_single_atom = cmd.select("temp_single_atom", f"name {j_atom} and chain {wg_atoms[0].chain} and resname {wg_atoms[0].resn} and bound_to_info_sele")
        temp_atoms = cmd.get_model("bound_to temp_single_atom").atom
        temp_wg_atoms = [ atom for atom in temp_atoms if atom.name[0]!="H"]
        k_atom = temp_wg_atoms[0].name
        single_hbd = HdbLine(bond_type, h_nums, h_atoms[0].name[0:2],i_atom, j_atom, k_atom)
        print(single_hbd.gen_hdb())
            

    # 打印氢原子信息
    # for bond_atom in bound_to_info:
    #     print(f"Hydrogen Name: {bond_atom.name}, Residue: {bond_atom.resn}, Chain: {bond_atom.chain}")
    # print(len(h_atoms))


# 一个H：确认只与一个氢相连
# 如果某重原子周围只与2个重原子相连，则类型必为1
# 如果某重原子周围只与1个重原子相连，则类型必为2，且需要获得相连重原子的另一个相连重原子
# 如果某重原子周围只与3以上个重原子相连，则类型必为5，且需要获得相连重原子的另一个相连重原子
# 确认只与2个氢相连
# 如果某重原子周围与2个以上重原子相连，则类型必为6
# 如果某重原子周围只与1个重原子相连，则类型必为3，且需要获得相连重原子的另一个相连重原子