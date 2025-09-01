# -*- coding: utf-8 -*-
# 单链稳定构象预设（NVT-only, 大盒真空）
import os
import datetime
from .. import lammps, preset
from ..md import MD
from ...core import utils

__version__ = '0.1.0'

class SingleChainRelax(preset.Preset):
    """
    目标：生成单条聚合物链在真空/大盒中的稳定构象
    - 仅使用 NVT；不做 NPT、不装填、不 21-step
    - 默认构建大盒，避免周期自相互作用
    - 步数短、可快速完成
    """

    def __init__(self, mol, prefix='', work_dir=None, save_dir=None, solver_path=None, **kwargs):
        super().__init__(mol, prefix=prefix, work_dir=work_dir, save_dir=save_dir, solver_path=solver_path, **kwargs)

        # 文件命名（单阶段）
        self.in_file     = kwargs.get('in_file',  f'{self.prefix}sc_relax.in')
        self.dat_file    = kwargs.get('dat_file', f'{self.prefix}sc_relax.data')
        self.pdb_file    = kwargs.get('pdb_file', f'{self.prefix}sc_relax.pdb')
        self.log_file    = kwargs.get('log_file', f'{self.prefix}sc_relax.log')
        self.dump_file   = kwargs.get('dump_file',f'{self.prefix}sc_relax.dump')
        self.xtc_file    = kwargs.get('xtc_file', f'{self.prefix}sc_relax.xtc')
        self.last_str    = kwargs.get('last_str', f'{self.prefix}sc_relax_last.dump')
        self.last_data   = kwargs.get('last_data',f'{self.prefix}sc_relax_last.data')
        self.pickle_file = kwargs.get('pickle_file',f'{self.prefix}sc_relax_last.pickle')
        self.json_file   = kwargs.get('json_file',  f'{self.prefix}sc_relax_last.json')

    # --- 构建 NVT-only 的流程 ---
    def _build_md(self,
                  temp=300.0,
                  high_temp=None,               # 如需短高温跨势垒：例如 500.0；否则 None
                  min_steps=5000,
                  pre_nvt_steps=5000,           # 300K 快速预热（0.1 fs）
                  cool_steps=100000,            # 若 high_temp!=None，则高温→目标温度的渐冷段
                  final_nvt_steps=100000,       # 目标温度的短采样/稳定
                  time_step_low=0.1,            # 最初阶段（未 shake）的小步长 fs
                  time_step=1.0,                # 之后阶段（shake=True）步长 fs
                  box_len=120.0,                # 立方盒边长（Å），足够大以避免自相互作用
                  comm_cutoff=6.0,              # 通信 cutoff，真空链可设小一点
                  set_init_velocity=False,
                  **kwargs):

        md = MD()
        # 复用父类在准备期设定的相互作用样式（若已被 Preset 配置过）
        md.pair_style = getattr(self, 'pair_style', 'lj/cut')
        md.cutoff_in  = getattr(self, 'cutoff_in', 3.0)
        md.cutoff_out = getattr(self, 'cutoff_out', '')
        md.kspace_style = getattr(self, 'kspace_style', 'none')             # 真空链：通常不需要 PPPM
        md.kspace_style_accuracy = getattr(self, 'kspace_style_accuracy', '')
        md.bond_style = self.bond_style
        md.angle_style = self.angle_style
        md.dihedral_style = self.dihedral_style
        md.improper_style = self.improper_style
        md.neighbor = f'{self.neighbor_dis} bin'
        md.log_file = self.log_file
        md.dat_file = self.dat_file
        md.dump_file = self.dump_file
        md.xtc_file = self.xtc_file
        md.rst = True
        md.outstr = self.last_str
        md.write_data = self.last_data
        md.add.append(f'comm_modify cutoff {float(comm_cutoff)}')

        # 1) 初始最小化（小步长前的“干净起点”）
        md.add_min(min_style='cg')

        # 若你传了 min_add，就把它作为全局额外命令插入到 input 文件（头部）
        if 'min_add' in kwargs and kwargs['min_add']:
            # 支持字符串或字符串列表
            extra = kwargs['min_add']
            if isinstance(extra, (list, tuple)):
                md.add.extend(extra)
            else:
                md.add.append(str(extra))
        # 2) 300K 短预热（小步长、无 shake）
        if pre_nvt_steps and pre_nvt_steps > 0:
            md.add_md('nvt', int(pre_nvt_steps), time_step=time_step_low, shake=False,
                      t_start=temp, t_stop=temp, **kwargs)

        # 3) 可选的短高温/冷却（帮助越过局部极小）
        if high_temp is not None:
            hot_steps = kwargs.get('hot_steps', 50000)  # 50 ps @1 fs
            md.add_md('nvt', int(hot_steps), time_step=time_step, shake=True,
                      t_start=float(high_temp), t_stop=float(high_temp), **kwargs)
            if cool_steps and cool_steps > 0:
                md.add_md('nvt', int(cool_steps), time_step=time_step, shake=True,
                          t_start=float(high_temp), t_stop=float(temp), **kwargs)

        # 4) 目标温度短采样/稳态确认（NVT）
        if final_nvt_steps and final_nvt_steps > 0:
            md.add_md('nvt', int(final_nvt_steps), time_step=time_step, shake=True,
                      t_start=float(temp), t_stop=float(temp), **kwargs)

        # 5) 统一扩盒到大立方体（避免 PBC 自相互作用）
        #    用 deform 设置 x/y/z 到 [-L/2, +L/2]
        if len(md.wf) > 0:
            L = float(box_len)
            md.wf[-1].add_deform(dftype='final', deform_fin_lo=-L/2.0, deform_fin_hi=L/2.0, axis='xyz')

        # 6) 可选：初始化速度
        if set_init_velocity:
            md.set_init_velocity = temp

        return md

    def exec(self, confId=0,
            temp=300.0,
            high_temp=None,
            min_steps=5000, pre_nvt_steps=5000, cool_steps=100000, final_nvt_steps=100000,
            time_step_low=0.1, time_step=1.0,
            box_len=120.0, comm_cutoff=6.0,
            omp=1, mpi=1, gpu=0, intel='auto', opt='auto', **kwargs):
        """
        运行单链稳定构象生成（NVT-only, 大盒），并统计每个阶段耗时
        阶段：min → pre_nvt → hot → cool → final_nvt
        """
        utils.MolToPDBFile(self.mol, os.path.join(self.work_dir, self.pdb_file))
        lmp = lammps.LAMMPS(work_dir=self.work_dir, solver_path=self.solver_path)
        lmp.make_dat(self.mol, file_name=self.dat_file, confId=confId)

        # 统一的“公共字段”初始化器
        def _mk_md_base():
            md = MD()
            md.pair_style = getattr(self, 'pair_style', 'lj/cut')
            md.cutoff_in  = getattr(self, 'cutoff_in', 3.0)
            md.cutoff_out = getattr(self, 'cutoff_out', '')
            md.kspace_style = getattr(self, 'kspace_style', 'none')
            md.kspace_style_accuracy = getattr(self, 'kspace_style_accuracy', '')
            md.bond_style = self.bond_style
            md.angle_style = self.angle_style
            md.dihedral_style = self.dihedral_style
            md.improper_style = self.improper_style
            md.neighbor = f'{self.neighbor_dis} bin'
            md.log_file = self.log_file
            md.dat_file = self.dat_file
            md.dump_file = self.dump_file
            md.xtc_file = self.xtc_file
            md.rst = True
            md.outstr = self.last_str
            md.write_data = self.last_data
            md.add.append(f'comm_modify cutoff {float(comm_cutoff)}')
            return md

        stage_times = {}

        # ---- Stage 1: 最小化 ----
        if min_steps and min_steps > 0:
            md_min = _mk_md_base()
            md_min.add_min(min_style='cg')
            # 可选：把 min_add 插入到 input 头部
            if 'min_add' in kwargs and kwargs['min_add']:
                extra = kwargs['min_add']
                if isinstance(extra, (list, tuple)):
                    md_min.add.extend(extra)
                else:
                    md_min.add.append(str(extra))

            t0 = datetime.datetime.now()
            utils.radon_print('Stage[minimize] running...', level=1)
            self.mol = lmp.run(md_min, mol=self.mol, confId=confId,
                            input_file=self.in_file, last_str=self.last_str, last_data=self.last_data,
                            omp=omp, mpi=mpi, gpu=gpu, intel=intel, opt=opt)
            t1 = datetime.datetime.now()
            stage_times['minimize'] = str(t1 - t0)

        # ---- Stage 2: 300K 预热 ----
        if pre_nvt_steps and pre_nvt_steps > 0:
            md_pre = _mk_md_base()
            md_pre.add_md('nvt', int(pre_nvt_steps), time_step=time_step_low, shake=False,
                        t_start=temp, t_stop=temp, **kwargs)

            t0 = datetime.datetime.now()
            utils.radon_print('Stage[pre_nvt @300K] running...', level=1)
            self.mol = lmp.run(md_pre, mol=self.mol, confId=confId,
                            input_file=self.in_file, last_str=self.last_str, last_data=self.last_data,
                            omp=omp, mpi=mpi, gpu=gpu, intel=intel, opt=opt)
            t1 = datetime.datetime.now()
            stage_times['pre_nvt'] = str(t1 - t0)

        # ---- Stage 3: 高温段（可选） ----
        if high_temp is not None:
            hot_steps = kwargs.get('hot_steps', 50000)
            if hot_steps and hot_steps > 0:
                md_hot = _mk_md_base()
                md_hot.add_md('nvt', int(hot_steps), time_step=time_step, shake=True,
                            t_start=float(high_temp), t_stop=float(high_temp), **kwargs)

                t0 = datetime.datetime.now()
                utils.radon_print(f'Stage[hot @ {high_temp}K] running...', level=1)
                self.mol = lmp.run(md_hot, mol=self.mol, confId=confId,
                                input_file=self.in_file, last_str=self.last_str, last_data=self.last_data,
                                omp=omp, mpi=mpi, gpu=gpu, intel=intel, opt=opt)
                t1 = datetime.datetime.now()
                stage_times['hot'] = str(t1 - t0)

            # ---- Stage 4: 冷却段（可选） ----
            if cool_steps and cool_steps > 0:
                md_cool = _mk_md_base()
                md_cool.add_md('nvt', int(cool_steps), time_step=time_step, shake=True,
                            t_start=float(high_temp), t_stop=float(temp), **kwargs)

                t0 = datetime.datetime.now()
                utils.radon_print('Stage[cool high->target] running...', level=1)
                self.mol = lmp.run(md_cool, mol=self.mol, confId=confId,
                                input_file=self.in_file, last_str=self.last_str, last_data=self.last_data,
                                omp=omp, mpi=mpi, gpu=gpu, intel=intel, opt=opt)
                t1 = datetime.datetime.now()
                stage_times['cool'] = str(t1 - t0)

        # ---- Stage 5: 目标温度稳定 ----
        if final_nvt_steps and final_nvt_steps > 0:
            md_final = _mk_md_base()
            md_final.add_md('nvt', int(final_nvt_steps), time_step=time_step, shake=True,
                            t_start=float(temp), t_stop=float(temp), **kwargs)
            # 在最后一个阶段做统一扩盒（避免 PBC 自相互作用）
            if len(md_final.wf) > 0:
                L = float(box_len)
                md_final.wf[-1].add_deform(dftype='final',
                                        deform_fin_lo=-L/2.0, deform_fin_hi=L/2.0, axis='xyz')

            t0 = datetime.datetime.now()
            utils.radon_print('Stage[final_nvt @target] running...', level=1)
            self.mol = lmp.run(md_final, mol=self.mol, confId=confId,
                            input_file=self.in_file, last_str=self.last_str, last_data=self.last_data,
                            omp=omp, mpi=mpi, gpu=gpu, intel=intel, opt=opt)
            t1 = datetime.datetime.now()
            stage_times['final_nvt'] = str(t1 - t0)

        # ---- 汇总与落盘 ----
        utils.MolToJSON(self.mol, os.path.join(self.save_dir, self.json_file))
        utils.pickle_dump(self.mol, os.path.join(self.save_dir, self.pickle_file))

        # 总耗时
        total = datetime.timedelta(0)
        for k, v in stage_times.items():
            # v like '0:00:05.123456'
            h, m, s = v.split(':')
            # 简易解析（足够用）；或者直接不解析，只展示字符串
            # 这里保留字符串展示更稳妥
            pass
        utils.radon_print(f"Stage time summary: {stage_times}", level=1)
        return self.mol

