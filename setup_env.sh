#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# 交互输入（有默认值）
###############################################################################
read -rp "请输入 LAMMPS 安装目录 [默认: /root/lammps-22Jul2025]: " input_root
LAMMPS_ROOT="${input_root:-/root/lammps-22Jul2025}"

read -rp "请输入 LAMMPS 可执行文件路径 [默认: ${LAMMPS_ROOT}/lmp]: " input_exec
LAMMPS_EXEC="${input_exec:-${LAMMPS_ROOT}/lmp}"

read -rp "请输入 Conda 环境名称 [默认: chain_generator]: " input_env
CONDA_ENV="${input_env:-chain_generator}"

read -rp "请输入 Python 版本 [默认: 3.11]: " input_py
PY_VER="${input_py:-3.11}"

read -rp "是否将 conda-forge 设为严格最高优先级 (strict)? [默认: n] (y/n): " input_strict
STRICT_FLAG="${input_strict:-n}"

###############################################################################
# 工具函数
###############################################################################
log() { printf "\n[%s] %s\n" "$(date +'%F %T')" "$*"; }

detect_conda_sh() {
  local candidates=(
    "$HOME/miniconda3/etc/profile.d/conda.sh"
    "$HOME/anaconda3/etc/profile.d/conda.sh"
    "/opt/conda/etc/profile.d/conda.sh"
    "$HOME/.conda/etc/profile.d/conda.sh"
  )
  for f in "${candidates[@]}"; do
    [[ -f "$f" ]] && { echo "$f"; return 0; }
  done
  if command -v conda >/dev/null 2>&1; then
    local cbase
    cbase="$(dirname "$(dirname "$(command -v conda)")")"
    if [[ -f "$cbase/etc/profile.d/conda.sh" ]]; then
      echo "$cbase/etc/profile.d/conda.sh"; return 0
    fi
  fi
  return 1
}

append_to_bashrc_once() {
  local line="$1"
  local file="$HOME/.bashrc"
  grep -Fqx "$line" "$file" 2>/dev/null || echo "$line" >> "$file"
}

###############################################################################
# 1) 初始化 conda，创建并激活环境
###############################################################################
log "检测并加载 conda..."
CONDA_SH="$(detect_conda_sh || true)"
if [[ -z "${CONDA_SH:-}" ]]; then
  echo "未找到 conda.sh，请先安装 Miniconda/Anaconda 或将 conda 初始化到当前 shell。"; exit 1
fi
# shellcheck source=/dev/null
source "$CONDA_SH"

if [[ "${STRICT_FLAG,,}" == "y" ]]; then
  conda config --add channels conda-forge || true
fi

if ! conda env list | awk '{print $1}' | grep -Fxq "$CONDA_ENV"; then
  log "创建 conda 环境 $CONDA_ENV (python=${PY_VER})..."
  conda create -n "$CONDA_ENV" "python=${PY_VER}" -y
else
  log "环境 $CONDA_ENV 已存在。"
fi

log "激活环境..."
conda activate "$CONDA_ENV"

###############################################################################
# 2) 安装 Python 依赖：RadonPy 可编辑安装 + pip 包
###############################################################################
log "进入 RadonPy 并安装（editable 模式）..."
if [[ -d "RadonPy" ]]; then
  pushd RadonPy >/dev/null
  pip install -e .
  popd >/dev/null
else
  log "未发现 RadonPy 目录，跳过 pip install -e ."
fi

log "安装其他 pip 依赖..."
pip install aio_pika minio tqdm

###############################################################################
# 3) 写入 LAMMPS/OMPI 环境变量到 ~/.bashrc
###############################################################################
log "写入 LAMMPS 与 OMPI 环境变量到 ~/.bashrc..."
append_to_bashrc_once 'export OMPI_ALLOW_RUN_AS_ROOT=1'
append_to_bashrc_once 'export OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1'
append_to_bashrc_once "export LAMMPS_EXEC=${LAMMPS_EXEC}"
append_to_bashrc_once "export PATH=\${PATH}:${LAMMPS_ROOT}"

###############################################################################
# 4) 安装量子化学与相关包（conda-forge / psi4）
###############################################################################
log "安装 psi4 / qcengine / qcelemental / simple-dftd3 / dftd3-python..."
set +u
conda install -c conda-forge -c psi4 -y \
  psi4 rdkit qcengine qcelemental resp mdtraj matplotlib simple-dftd3 dftd3-python
set -u
###############################################################################
# 5) 对齐 C/C++ 运行时并修复 libstdc++ 解析
###############################################################################
log "安装/升级 C/C++ 运行时与 ncurses..."
conda install -c conda-forge "libstdcxx-ng>=14" "libgcc-ng>=14" ncurses -y

log "修复 libstdc++.so.6 符号链接..."
CONDA_LIB_DIR="${CONDA_PREFIX}/lib"
if [[ -f "${CONDA_LIB_DIR}/libstdc++.so.6" ]]; then
  log "已存在 ${CONDA_LIB_DIR}/libstdc++.so.6"
else
  LIBSTD_REAL="$(ls -1 ${CONDA_LIB_DIR}/libstdc++.so.6.* 2>/dev/null | sort -V | tail -n1 || true)"
  if [[ -n "${LIBSTD_REAL}" && -f "${LIBSTD_REAL}" ]]; then
    ln -sf "${LIBSTD_REAL}" "${CONDA_LIB_DIR}/libstdc++.so.6"
    log "已创建链接 ${CONDA_LIB_DIR}/libstdc++.so.6 -> ${LIBSTD_REAL}"
  else
    echo "未找到 ${CONDA_LIB_DIR}/libstdc++.so.6.*，请检查 libstdcxx-ng 安装是否成功。"; exit 1
  fi
fi
append_to_bashrc_once 'export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"'

###############################################################################
# 6) 安装 nvm / node / pm2
###############################################################################
log "安装 nvm..."
export NVM_DIR="$HOME/.nvm"
if [[ ! -d "$NVM_DIR" ]]; then
  curl -o- https://raw.githubusercontent.com/nvm-sh/nvm/v0.40.3/install.sh | bash
fi
# shellcheck source=/dev/null
. "$HOME/.nvm/nvm.sh"

log "安装 Node.js v22 与 pm2..."
nvm install 22
npm install pm2@latest -g

append_to_bashrc_once 'export NVM_DIR="$HOME/.nvm"'
append_to_bashrc_once '[ -s "$NVM_DIR/nvm.sh" ] && \. "$NVM_DIR/nvm.sh"'

###############################################################################
# 7) 创建 configs 目录
###############################################################################
log "创建 configs 目录..."
mkdir -p ./configs || true

###############################################################################
# 完成提示
###############################################################################
log "全部完成。请执行："
echo "1) source ~/.bashrc"
echo "2) conda activate ${CONDA_ENV}"
echo "3) 运行 python chain_relax_worker.py"
