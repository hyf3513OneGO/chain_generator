
import os
import time
import rdkit
import shutil
from rdkit import Chem
from rdkit.Chem import AllChem
from radonpy.sim.preset.single_chain import SingleChainRelax
from radonpy.ff.gaff2 import GAFF2
from radonpy.core import poly
import radonpy.core.utils as rad_utils
from radonpy.sim import qm 
import copy

from utils.chain_split_utils import remove_star_atoms

class HomoChainBuilder:
    def __init__(self,lmps_exec:str,work_dir:str,temp_dir:str,
    omp_psi4=1,mem_psi4=10000,conf_mm_omp=16,conf_mm_mpi=4,conf_mm_gpu=1,conf_mm_mp=1
    ) -> rdkit.Chem.rdchem.Mol:
        self.lmps_exec = lmps_exec
        os.environ["OMPI_ALLOW_RUN_AS_ROOT"] = "1"
        os.environ["OMPI_ALLOW_RUN_AS_ROOT_CONFIRM"] = "1"
        os.environ["LAMMPS_EXEC"] = lmps_exec
        self.work_dir = work_dir
        self.temp_dir = temp_dir
        self.omp_psi4 = omp_psi4
        self.mem_psi4 = mem_psi4
        self.conf_mm_omp = conf_mm_omp
        self.conf_mm_mpi = conf_mm_mpi
        self.conf_mm_gpu = conf_mm_gpu
        self.conf_mm_mp = conf_mm_mp
        if not os.path.exists(temp_dir):
            os.makedirs(temp_dir)
        if not os.path.exists(work_dir):
            os.makedirs(work_dir)
    def create_mol(self,psmiles,n_repeat=None,n_atoms=None,random_walk=False):        
        mol = rad_utils.mol_from_smiles(psmiles)
        key_point_list, processed_mol, processed_smi = remove_star_atoms(psmiles)
        n1,s1,n2,s2 = key_point_list
        ter_tail_symbol = processed_mol.GetAtomWithIdx(n1).GetSymbol()
        # ter_tail_symbol = ""
        ter_head_symbol = processed_mol.GetAtomWithIdx(n2).GetSymbol()
        ter_tail_smi = f'*{ter_tail_symbol}'
        ter_head_smi = f'*{ter_head_symbol}'
        ter_tail = rad_utils.mol_from_smiles(ter_tail_smi)
        ter_head = rad_utils.mol_from_smiles(ter_head_smi)
        params = AllChem.ETKDGv3()
        params.randomSeed = 42
        AllChem.EmbedMolecule(mol, params)
        AllChem.MMFFOptimizeMolecule(mol)
        qm.assign_charges(mol, charge='gasteiger', work_dir=self.work_dir, tmp_dir=self.temp_dir, opt=False, omp=self.omp_psi4, memory=self.mem_psi4)
        qm.assign_charges(ter_head, charge='gasteiger', work_dir=self.work_dir, tmp_dir=self.temp_dir, opt=False, omp=self.omp_psi4, memory=self.mem_psi4)
        qm.assign_charges(ter_tail, charge='gasteiger', work_dir=self.work_dir, tmp_dir=self.temp_dir, opt=False, omp=self.omp_psi4, memory=self.mem_psi4)
        if n_repeat is None and n_atoms is not None:
            n = poly.calc_n_from_num_atoms(mol, n_atoms, terminal1=ter_head,terminal2=ter_tail)
        elif n_repeat is not None:
            n = n_repeat
        n = int(n)
        if random_walk:
            print("use random walk")
            homopoly = poly.polymerize_rw(mol, n)
            homopoly = poly.terminate_rw(homopoly,mol1=ter_head,mol2=ter_tail,terminate_tail=True)
        else:
            print("not use random walk")
            homopoly = poly.polymerize_mols(mol, n)
            homopoly = poly.terminate_mols(homopoly,mol1=ter_head,mol2=ter_tail,terminate_tail=True)
        return homopoly
    def relax_conf(self,mol:rdkit.Chem.rdchem.Mol,save_dir:str,args_dict:dict)->rdkit.Chem.rdchem.Mol:
        confID = args_dict.get('confID',0)
        temp = args_dict.get('temp',298.15)
        high_temp = args_dict.get('high_temp',None)
        prev_nvt_steps = args_dict.get('prev_nvt_steps',5000)
        cool_steps = args_dict.get('cool_steps',100000)
        final_nvt_steps = args_dict.get('final_nvt_steps',100000)
        time_step_low = args_dict.get('time_step_low',0.1)
        time_step = args_dict.get('time_step',1.0)
        box_length = args_dict.get('box_length',120.0)
        comm_cutoff = args_dict.get('comm_cutoff',6.0)
        ff = GAFF2()
        _ = ff.ff_assign(mol)   # 为分子对象写入 GAFF2 参数
        relax = SingleChainRelax(mol, work_dir=self.work_dir, save_dir=save_dir, prefix='poly_')
        mol_relaxed = relax.exec(
            confId=confID,
            temp=temp,
            high_temp=high_temp,          # 若初始构象较差，可设 high_temp=500.0，并保留默认 cool_steps
            pre_nvt_steps=prev_nvt_steps,
            cool_steps=cool_steps,
            final_nvt_steps=final_nvt_steps,
            time_step_low=time_step_low,
            time_step=time_step,           # 若启用 shake=True（默认），1.0–2.0 fs 都可
            box_len=box_length,           # 立方大盒边长（Å），按你的链尺寸可再调大/调小
            comm_cutoff=comm_cutoff,
            omp=self.omp_psi4, mpi=self.conf_mm_mpi, gpu=self.conf_mm_gpu
        )
        return mol_relaxed
    def fromSmiles(self,psmiles,n_repeat=None,n_atoms=None,save_dir:str=None,random_walk=False,args_dict:dict=None):
        import time
        
        # 创建分子
        start_time = time.time()
        mol_init = self.create_mol(psmiles,n_repeat,n_atoms,random_walk)
        create_mol_time = time.time() - start_time
        
        # 复制分子
        start_time = time.time()
        mol = Chem.Mol(mol_init)
        copy_mol_time = time.time() - start_time
        
        # 松弛构象
        start_time = time.time()
        mol_relaxed = self.relax_conf(mol,save_dir,args_dict)
        relax_conf_time = time.time() - start_time
        
        print(f"[时间统计] 创建分子: {create_mol_time:.3f}s, 复制分子: {copy_mol_time:.3f}s, 松弛构象: {relax_conf_time:.3f}s")
        return mol_relaxed,mol_init
    def save_mol(self,mol:rdkit.Chem.rdchem.Mol,save_path:str):
        writer = Chem.SDWriter(save_path)
        writer.write(mol)
        writer.close()
if __name__ == '__main__':
    hcb = HomoChainBuilder(lmps_exec='/data/yifei/lammps-22Jul2025/src/lmp_mpi',work_dir='./work_dir',temp_dir='./temp_dir')
    args_dict={'temp':298,'high_temp':600.0,'prev_nvt_steps':60000,'cool_steps':500000,
    'final_nvt_steps':1000000,'time_step_low':0.5,'time_step':2.0,'box_length':150.0,'comm_cutoff':12.0}
    homopoly_relaxed = hcb.fromSmiles('*OC(*)=O',n_repeat=6,save_dir='./save_dir',
    random_walk=True,args_dict=args_dict)
    hcb.save_mol(homopoly_relaxed,'./homopoly_relaxed.sdf')
