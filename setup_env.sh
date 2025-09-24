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

https_ok() {
  if command -v curl >/dev/null 2>&1; then
    curl -sI https://baidu.com >/dev/null 2>&1
  elif command -v wget >/dev/null 2>&1; then
    wget -q --spider https://baidu.com >/dev/null 2>&1
  else
    return 1
  fi
}

###############################################################################
# 0) 不全局污染系统库路径
###############################################################################
log "不会向 ~/.bashrc 写入 LD_LIBRARY_PATH；仅在 Conda 环境内设置。"

###############################################################################
# 1) 初始化 conda，创建并激活环境
###############################################################################
log "检测并加载 conda..."
CONDA_SH="$(detect_conda_sh || true)"
if [[ -z "${CONDA_SH:-}" ]]; then
  echo "未找到 conda.sh。请先安装 Miniconda/Anaconda 并确保 'conda init' 后可用。"
  exit 1
fi
# shellcheck source=/dev/null
source "$CONDA_SH"

# 使用更快的求解器（可显著加速依赖求解）
log "设置求解器为 libmamba（加速依赖解析）..."
conda config --set solver libmamba || true
# 不修改 channel_priority（保持默认 flexible）
# 不全局添加 conda-forge，避免整体求解变慢

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
log "安装 RadonPy（editable）..."
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
# 3) 写入 LAMMPS/OMPI 到 ~/.bashrc（不涉及 LD_LIBRARY_PATH）
###############################################################################
log "写入 LAMMPS 与 OMPI 环境变量到 ~/.bashrc..."
append_to_bashrc_once 'export OMPI_ALLOW_RUN_AS_ROOT=1'
append_to_bashrc_once 'export OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1'
append_to_bashrc_once "export LAMMPS_EXEC=${LAMMPS_EXEC}"
append_to_bashrc_once "export PATH=\${PATH}:${LAMMPS_ROOT}"

###############################################################################
# 4) 安装量子化学与相关包（只在该命令使用 conda-forge/psi4 渠道）
###############################################################################
log "安装 psi4 / rdkit / qcengine / qcelemental / resp / mdtraj / matplotlib / simple-dftd3 / dftd3-python..."
# 渠道仅对本次安装生效，不改变全局配置，减少求解成本
set +u
conda install -y -c conda-forge -c psi4 \
  psi4 rdkit qcengine qcelemental resp mdtraj matplotlib simple-dftd3 dftd3-python
set -u

###############################################################################
# 5) 对齐 C/C++ 运行时并修复 libstdc++（仅在该环境内）
###############################################################################
log "安装/升级 C/C++ 运行时与 ncurses（环境内）..."
conda install -y -c conda-forge "libstdcxx-ng>=14" "libgcc-ng>=14" ncurses

log "确保 libstdc++.so.6 符号链接存在（环境内）..."
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

###############################################################################
# 6) 仅在该 Conda 环境内注入 LD_LIBRARY_PATH（不写 ~/.bashrc）
###############################################################################
log "将 LD_LIBRARY_PATH 限定为仅在激活 ${CONDA_ENV} 时生效..."
conda env config vars set LD_LIBRARY_PATH="${CONDA_PREFIX}/lib"
conda deactivate
conda activate "$CONDA_ENV"

###############################################################################
# 7) 安装 nvm / node / pm2（仅当 HTTPS 正常时）
###############################################################################
# 检测并安装 curl（如果不存在）
if ! command -v curl &> /dev/null; then
  log "curl 未安装，正在安装 curl..."
  if command -v apt-get &> /dev/null; then
    sudo apt-get update -y && sudo apt-get install -y curl
  elif command -v yum &> /dev/null; then
    sudo yum install -y curl
  elif command -v dnf &> /dev/null; then
    sudo dnf install -y curl
  elif command -v pacman &> /dev/null; then
    sudo pacman -S --noconfirm curl
  else
    log "无法自动安装 curl，请手动安装后重新运行脚本"
    exit 1
  fi
  log "curl 安装完成"
else
  log "curl 已安装"
fi

log "检测 HTTPS 可用性，用于安装 nvm/node..."
if https_ok; then
  log "HTTPS 正常，安装 nvm..."
  export NVM_DIR="$HOME/.nvm"
  if [[ ! -d "$NVM_DIR" ]]; then
    curl -o- https://raw.githubusercontent.com/nvm-sh/nvm/v0.40.3/install.sh | bash
  fi
  # shellcheck source=/dev/null
  . "$HOME/.nvm/nvm.sh"

  log "安装 Node.js v22 与 pm2..."
  # 设定环境变量
  export NVM_NODEJS_ORG_MIRROR=https://mirrors.cernet.edu.cn/nodejs-release/
  # 然后正常使用
  nvm install 22
  npm install pm2@latest -g

  append_to_bashrc_once 'export NVM_DIR="$HOME/.nvm"'
  append_to_bashrc_once '[ -s "$NVM_DIR/nvm.sh" ] && \. "$NVM_DIR/nvm.sh"'
else
  log "HTTPS 不可用，跳过 nvm/node/pm2 安装。修复系统 TLS 后可重跑本段。"
fi

###############################################################################
# 8) 创建 configs 目录
###############################################################################
log "创建 configs 目录..."
mkdir -p ./configs || true

###############################################################################
# 完成提示
###############################################################################
log "全部完成。后续步骤："
echo "1) source ~/.bashrc"
echo "2) conda activate ${CONDA_ENV}"
echo "3) 运行 python chain_relax_worker.py"
echo "说明：未启用 conda-forge strict，且不全局添加 conda-forge；仅在相关安装命令使用对应渠道。已启用 libmamba 求解器以提升速度。"
