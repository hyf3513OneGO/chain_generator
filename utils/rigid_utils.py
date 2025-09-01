import numpy as np
from typing import Tuple
def gram_schmidt(v1: np.ndarray, v2: np.ndarray) -> np.ndarray:
    """
    使用 Gram-Schmidt 构造右手正交基 (R)
    v1 定义 x 轴方向，v2 协助定义 y 轴
    返回旋转矩阵 R: shape (3,3)
    """
    x = v1 / np.linalg.norm(v1)
    v2_proj = v2 - np.dot(v2, x) * x
    y = v2_proj / np.linalg.norm(v2_proj)
    z = np.cross(x, y)
    return np.stack([x, y, z], axis=1)
def to_local_coords(X: np.ndarray, R: np.ndarray, t: np.ndarray) -> np.ndarray:
    """
    将全局坐标变换为局部坐标
    X: shape (N, 3)
    R: (3, 3)
    t: (3,)
    """
    return (R.T @ (X - t).T).T

def to_global_coords(X_local: np.ndarray, R: np.ndarray, t: np.ndarray) -> np.ndarray:
    """
    将局部坐标映射回全局坐标
    """
    return (R @ X_local.T).T + t
def compute_rigid_frame_from_three_atoms(
    x1: np.ndarray, x2: np.ndarray, x3: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """
    构造局部坐标系：
    - x1, x2: 定义第一个向量 v1 = x2 - x1（前一方向）
    - x3: 当前 monomer 的 anchor，定义原点和方向 v2 = x2 - x3（后一方向）
    """
    v1 = x1 - x2
    v2 = x3 - x2  # 或 x3 - x2，只要保持一致即可
    R = gram_schmidt(v1, v2)
    t = x2
    return R, t