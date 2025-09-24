#!/usr/bin/env bash
# 安装并编译 CPU 版本的 LAMMPS（22Jul2025），并写入环境变量到 ~/.bashrc

set -euo pipefail

# ---------- 默认值 ----------
DEFAULT_URL="http://47.117.17.97:8000/lammps-22Jul2025.tar.gz"
DEFAULT_PKG="lammps-22Jul2025.tar.gz"
DEFAULT_SRC="lammps-22Jul2025"
BUILD_DIR="build-cpu"
INSTALL_LMP_NAME="lmp"
JOBS="$(nproc || echo 4)"

# ---------- 用户输入覆盖 ----------
read -r -p "请输入源码下载 URL [默认: $DEFAULT_URL]：" INPUT_URL
URL="${INPUT_URL:-$DEFAULT_URL}"

read -r -p "请输入源码包文件名 [默认: $DEFAULT_PKG]：" INPUT_PKG
PKG="${INPUT_PKG:-$DEFAULT_PKG}"

read -r -p "请输入源码目录名 [默认: $DEFAULT_SRC]：" INPUT_SRC
SRC_DIR="${INPUT_SRC:-$DEFAULT_SRC}"

# ---------- 前置依赖 ----------
echo "[1/6] 安装系统依赖..."
sudo apt-get update -y
sudo apt-get install -y --no-install-recommends \
  build-essential cmake git \
  mpi-default-bin mpi-default-dev \
  libfftw3-dev wget ca-certificates

# ---------- 获取源码 ----------
echo "[2/6] 下载源码..."
wget -O "$PKG" "$URL"

echo "[3/6] 解压源码..."
rm -rf "$SRC_DIR"
tar -xzvf "$PKG"

# ---------- CMake 配置 ----------
echo "[4/6] 配置 CMake..."
cd "$SRC_DIR"
rm -rf "$BUILD_DIR"
mkdir -p "$BUILD_DIR"
cd "$BUILD_DIR"

cmake ../cmake \
  -D CMAKE_BUILD_TYPE=Release \
  -D BUILD_MPI=on \
  -D FFT=FFTW3 \
  \
  -D PKG_KSPACE=on \
  -D PKG_MANYBODY=on \
  -D PKG_MOLECULE=on \
  -D PKG_REAXFF=on \
  -D PKG_CLASS2=on \
  -D PKG_RIGID=on \
  -D PKG_MISC=on \
  -D PKG_EXTRA-MOLECULE=on \
  -D PKG_EXTRA-DUMP=on

# ---------- 编译 ----------
echo "[5/6] 编译..."
cmake --build . -j "${JOBS}"

# ---------- 安装 ----------
echo "[6/6] 安装可执行文件..."
cp -f ./lmp ../"${INSTALL_LMP_NAME}"
chmod +x ../"${INSTALL_LMP_NAME}"
cd ..

LMP_PATH="$(pwd)/${INSTALL_LMP_NAME}"
LMP_DIR="$(pwd)"

echo "编译完成，可执行文件：$LMP_PATH"

# ---------- 写入环境变量 ----------
echo "写入环境变量到 ~/.bashrc ..."
cp -a ~/.bashrc ~/.bashrc.bak.$(date +%F-%H%M%S)

cat >> ~/.bashrc <<EOF
# --- LAMMPS & OpenMPI (auto-added) ---
export OMPI_ALLOW_RUN_AS_ROOT=1
export OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1

export LAMMPS_EXEC=${LMP_PATH}

LAMMPS_SRC_DIR=${LMP_DIR}
case ":\\$PATH:" in
  *":\${LAMMPS_SRC_DIR}:"*) : ;;
  *) export PATH="\\$PATH:\${LAMMPS_SRC_DIR}" ;;
esac
# --- end ---
EOF

# 立即生效
set +u
# shellcheck disable=SC1090
source ~/.bashrc || true
set -u

# ---------- 提示 ----------
cat <<EOF

运行示例：
  mpirun -np 1 ./lmp -in in.shear

多核示例（4个核心）：
  mpirun -np 4 ./lmp -in in.shear

环境变量已写入 ~/.bashrc，并已生效。
可执行文件路径: $LMP_PATH
EOF
