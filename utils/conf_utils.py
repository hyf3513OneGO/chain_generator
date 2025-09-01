import numpy as np
import rdkit
from rdkit.Chem import AllChem
import rdkit.Chem as Chem
def get_coord_from_mol(mol:rdkit.Chem.rdchem.Mol)->list:
    conf = mol.GetConformer()
    coord = []
    for i in range(mol.GetNumAtoms()):
        coord.append(conf.GetAtomPosition(i))
    return coord

def set_coord_to_mol(mol:rdkit.Chem.rdchem.Mol,coord:list)->rdkit.Chem.rdchem.Mol:
    conf = mol.GetConformer()
    for i in range(mol.GetNumAtoms()):
        conf.SetAtomPosition(i,coord[i])
    return mol
def merge_monomer(mol_list:list,key_point_list:list)->rdkit.Chem.rdchem.Mol:
    n1,s1,n2,s2 = key_point_list
    exclude_idxs_list = []
    for idx,mol in enumerate(mol_list):
        if idx==0:
            exclude_idxs_list.append([s2])
        elif idx==len(mol_list)-1:
            exclude_idxs_list.append([s1])
        else:
            exclude_idxs_list.append([s1,s2])
    merged = Chem.RWMol()
    atom_map = {} 
    for idx,mol in enumerate(mol_list):
        for atom_idx in range(mol.GetNumAtoms()):
            if atom_idx in exclude_idxs_list[idx]:
                continue
            atom = mol.GetAtomWithIdx(atom_idx)
            new_atom = Chem.Atom(atom)  # 完整复制
            new_idx = merged.AddAtom(new_atom)
            atom_map[(idx, atom_idx)] = new_idx
        for bond in mol.GetBonds():
            a1, a2 = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            if (a1 in exclude_idxs_list[idx] or a2 in exclude_idxs_list[idx]):
                continue
            new_a1 = atom_map[(idx, a1)]
            new_a2 = atom_map[(idx, a2)]
            merged.AddBond(new_a1, new_a2, bond.GetBondType())
            if bond.GetIsAromatic():
                merged.GetBondBetweenAtoms(new_a1, new_a2).SetIsAromatic(True)
        if idx>0:
            # if idx==1 and s1<n2:
            #     prev_n2_idx = n2-1
            # else:
            prev_n2_idx = n2
            new_prev_idx = atom_map[(idx-1,prev_n2_idx)]
            new_curr_idx = atom_map[(idx,n1)]
            merged.AddBond(new_prev_idx,new_curr_idx,Chem.rdchem.BondType.SINGLE)
    merged_mol = merged.GetMol()
    conf = Chem.Conformer(merged_mol.GetNumAtoms())
    for (i, old_idx), new_idx in atom_map.items():
        pos = mol_list[i].GetConformer().GetAtomPosition(old_idx)
        conf.SetAtomPosition(new_idx, pos)
    Chem.SanitizeMol(merged_mol, sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE)
    merged_mol.AddConformer(conf, assignId=True)
    return merged_mol,atom_map