from rdkit import Chem
from utils.conf_utils import get_coord_from_mol, merge_monomer, set_coord_to_mol
from utils.rigid_utils import compute_rigid_frame_from_three_atoms, to_global_coords, to_local_coords
import rdkit
from typing import Tuple
import numpy as np
import os

def remove_star_atoms(monomer_smiles):
    mol = Chem.MolFromSmiles(monomer_smiles)
    assert 'H' not in [atom.GetSymbol() for atom in mol.GetAtoms()]
    key_point_list = []
    for atom in mol.GetAtoms():
        for neighbor in atom.GetNeighbors():
            if neighbor.GetSymbol() == '*':
                key_point_list.append(atom.GetIdx())
                key_point_list.append(neighbor.GetIdx())

    assert len(key_point_list) == 4, 'Invalid PSMILES: should have two [*] atoms.'

    star_0, star_1 = key_point_list[1], key_point_list[3]
    neighbor_0, neighbor_1 = key_point_list[0], key_point_list[2]
    # 修改原子类型以便还原真实连接
    mol.GetAtomWithIdx(star_0).SetAtomicNum(mol.GetAtomWithIdx(neighbor_1).GetAtomicNum())
    mol.GetAtomWithIdx(star_1).SetAtomicNum(mol.GetAtomWithIdx(neighbor_0).GetAtomicNum())

    # 用真实原子替换 SMILES 中的 [*]
    replacements = [mol.GetAtomWithIdx(neighbor_1).GetSymbol(), mol.GetAtomWithIdx(neighbor_0).GetSymbol()]
    processed_smi = ""
    count = 0
    i = 0
    while i < len(monomer_smiles):
        if monomer_smiles[i:i + 3] == "[*]":
            processed_smi += replacements[count]
            count += 1
            i += 3
        else:
            processed_smi += monomer_smiles[i]
            i += 1
    pre_atoms = [atom.GetSymbol() for atom in mol.GetAtoms()]
    processed_mol = Chem.MolFromSmiles(processed_smi)
    post_atoms = [atom.GetSymbol() for atom in processed_mol.GetAtoms()]
    assert pre_atoms == post_atoms, 'Atom order mismatch after replacement.'
    return key_point_list, processed_mol, processed_smi
def extract_monomers(mol=None,sdf_file=None, psmiles=None):
    if mol is None and sdf_file is not None:
        mol_noH = Chem.MolFromMolFile(sdf_file, removeHs=True)
    else:
        mol_noH = mol
    keypoints, processed_mol, processed_smi = remove_star_atoms(psmiles)
    n1, star1, n2, star2 = keypoints
    monomer_noH = Chem.RemoveHs(Chem.MolFromSmiles(processed_smi))
    n_monomer_atoms = monomer_noH.GetNumAtoms()
    n_poly_atoms = mol_noH.GetNumAtoms()
    n_repeat = (n_poly_atoms-2) //(n_monomer_atoms - 2)
    total_index = list(range(n_poly_atoms))
    filled_index = []
    cursor = 1
    for _ in range(n_repeat):
        idx_map = [-1] * n_monomer_atoms
        for j in range(n_monomer_atoms):
            if j != star1 and j != star2:
                idx_map[j] = total_index[cursor]
                cursor += 1
        filled_index.append(idx_map)
    for i in range(n_repeat):
        if i > 0:
            filled_index[i][star1] = filled_index[i - 1][n2]
        if i < n_repeat - 1:
            filled_index[i][star2] = filled_index[i + 1][n1]
    filled_index[0][star1] = 0
    filled_index[-1][star2] = n_poly_atoms - 1
    # 根据filled_index提取monomer分子列表
    monomer_list = []
    for idx_map in filled_index:
        submol,keypoint_list_new=extract_submol_by_atoms(mol_noH,idx_map,keypoints)
        monomer_list.append(submol)
    
    return keypoints,filled_index,mol_noH,monomer_list

def extract_rigid_from_monomers(monomer_list,keypoints)->Tuple[np.ndarray,np.ndarray,list[rdkit.Chem.rdchem.Mol]]:
    n1,s1,n2,s2 = keypoints

    # for momnomer in monomer_list[1:-1]:
    #     x1 = momnomer.GetConformer().GetAtomPosition(s1)
    #     x2 = momnomer.GetConformer().GetAtomPosition(n1)
    #     x3 = momnomer.GetConformer().GetAtomPosition(s2)
    # x1 = monomer_list[-1].GetConformer().GetAtomPosition(s1)
    # x2 = monomer_list[-1].GetConformer().GetAtomPosition(n1)
    # x3 = monomer_list[-1].GetConformer().GetAtomPosition(n2)
    R_list = []
    t_list = []
    local_coords = []
    global_coords = []
    local_monomers = []

    for idx,monomer in enumerate(monomer_list):
        x1 = np.array([monomer.GetConformer().GetAtomPosition(n1).x,monomer.GetConformer().GetAtomPosition(n1).y,monomer.GetConformer().GetAtomPosition(n1).z])
        x2 = np.array([monomer.GetConformer().GetAtomPosition(s1).x,monomer.GetConformer().GetAtomPosition(s1).y,monomer.GetConformer().GetAtomPosition(s1).z])
        x3 = np.array([monomer.GetConformer().GetAtomPosition(s2).x,monomer.GetConformer().GetAtomPosition(s2).y,monomer.GetConformer().GetAtomPosition(s2).z])
        R, t = compute_rigid_frame_from_three_atoms(x1, x2, x3)
        R_list.append(R)
        t_list.append(t)
        global_coord = np.array(get_coord_from_mol(monomer))
        local_coord = to_local_coords(global_coord,R,t)
        local_coords.append(local_coord)
        global_coords.append(global_coord)
    for monomer, local_coord in zip(monomer_list, local_coords):
        # 创建新分子，复制原子和键
        emol = Chem.RWMol()
        for atom in monomer.GetAtoms():
            new_atom = Chem.Atom(atom.GetAtomicNum())
            emol.AddAtom(new_atom)
        for bond in monomer.GetBonds():
            emol.AddBond(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), bond.GetBondType())
        new_mol = emol.GetMol()
        # 添加局部坐标构象
        conf = Chem.Conformer(monomer.GetNumAtoms())
        for i in range(monomer.GetNumAtoms()):
            x, y, z = local_coord[i]
            conf.SetAtomPosition(i, (float(x), float(y), float(z)))
        new_mol.RemoveAllConformers()
        new_mol.AddConformer(conf, assignId=True)
        local_monomers.append(new_mol)
    return R_list, t_list, local_monomers
def compose_monomers_from_local_monomers(local_monomers:list,R_list:list,t_list:list,keypoints:list)->rdkit.Chem.rdchem.Mol:
    global_monomers = []
    for local_monomer, R, t in zip(local_monomers, R_list, t_list):
        # 这里需要 new 一下，复制一份分子对象
        new_monomer = Chem.Mol(local_monomer)
        local_coord = np.array(get_coord_from_mol(new_monomer))
        global_coord = to_global_coords(local_coord, R, t)
        global_monomer = set_coord_to_mol(new_monomer, global_coord)
        global_monomers.append(global_monomer)
    merged_mol, atom_map = merge_monomer(global_monomers, keypoints)
    return merged_mol

def extract_submol_by_atoms(mol:rdkit.Chem.rdchem.Mol, atom_indices:list,keypoint_list:list)->Tuple[rdkit.Chem.rdchem.Mol,list]:
    """
    从mol中提取由atom_indices指定的子分子（包括键），保留3D构象。
    """
    # 创建原子掩码
    amap = {}
    for i, idx in enumerate(atom_indices):
        amap[idx] = i
    n1_H = keypoint_list[0]
    s1_H = keypoint_list[1]
    n2_H = keypoint_list[2]
    s2_H = keypoint_list[3]
    keypoint_list_new = [-1,-1,-1,-1]
    # 创建新分子
    emol = Chem.RWMol()
    for idx in atom_indices:
      
        atom = mol.GetAtomWithIdx(idx)
 
        new_atom = Chem.Atom(atom.GetAtomicNum())

        new_idx  = emol.AddAtom(new_atom)
        if idx == n1_H:
            keypoint_list_new[0]=new_idx
        if idx == s1_H:
            keypoint_list_new[1]=new_idx
        if idx ==n2_H:
            keypoint_list_new[2]=new_idx
        if idx ==s2_H:
            keypoint_list_new[3]=new_idx 

    bond_set = set()
    for bond in mol.GetBonds():
        begin_idx = bond.GetBeginAtomIdx()
        end_idx = bond.GetEndAtomIdx()
        if begin_idx in amap and end_idx in amap:
            bond_type = bond.GetBondType()
            emol.AddBond(amap[begin_idx], amap[end_idx], bond_type)
            bond_set.add((amap[begin_idx], amap[end_idx]))

    submol = emol.GetMol()

    # 添加构象信息
    conf = mol.GetConformer()
    sub_conf = Chem.Conformer(len(atom_indices))
    for i, idx in enumerate(atom_indices):
        pos = conf.GetAtomPosition(idx)
        sub_conf.SetAtomPosition(i, pos)
    submol.RemoveAllConformers()
    submol.AddConformer(sub_conf)    
    return submol,keypoint_list_new
def save_to_sdf(mol,sdf_path:str):
    writer = Chem.SDWriter(sdf_path)
    writer.write(mol)
    writer.close()
def extract_save_monomers(mol=None,sdf_path=None,psmiles=None,save_dir=None)->list[rdkit.Chem.rdchem.Mol]:
    from rdkit.Chem import AllChem, rdMolAlign
    if mol is not None:
        keypoints,filled_index,mol_noH,monomer_list = extract_monomers(mol=mol,psmiles=psmiles)
    if sdf_path is not None:
        keypoints,filled_index,mol_noH,monomer_list = extract_monomers(sdf_file=sdf_path,psmiles=psmiles)
    if mol is None and sdf_path is None:
        raise ValueError("Either mol or sdf_path must be provided")
    R_list, t_list, local_monomers = extract_rigid_from_monomers(monomer_list,keypoints)
    merged_mol = compose_monomers_from_local_monomers(local_monomers,R_list,t_list,keypoints)
    rmsd = rdMolAlign.GetBestRMS(merged_mol,mol_noH)
    assert rmsd<1e-4, "The merged mol is not the same as the original mol"
    for i,monomer in enumerate(local_monomers):
        save_to_sdf(monomer,os.path.join(save_dir,f"local_monomer_{i}.sdf"))
    for i,monomer in enumerate(monomer_list):
        save_to_sdf(monomer,os.path.join(save_dir,f"global_monomer_{i}.sdf"))
    np.save(os.path.join(save_dir,"R_list.npy"),R_list)
    np.save(os.path.join(save_dir,"t_list.npy"),t_list)
    np.save(os.path.join(save_dir,"keypoints.npy"),keypoints)
    np.save(os.path.join(save_dir,"filled_index.npy"),filled_index)
    return local_monomers,R_list,t_list