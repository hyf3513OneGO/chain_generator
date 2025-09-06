from minio import Minio
from minio.error import S3Error
import os
import json
from tqdm import tqdm
class MinioConfig:
    def __init__(self, config_path):
        self.config_path = config_path
        with open(self.config_path, "r") as f:
            self.config = json.load(f)
            self.endpoint = self.config["minio"]["endpoint"]
            self.access_key = self.config["minio"]["access_key"]
            self.secret_key = self.config["minio"]["secret_key"]
            self.secure = self.config["minio"]["secure"]

class MinioManager:
    def __init__(self, endpoint, access_key, secret_key, secure=False):
        """
        初始化 Minio 客户端
        :param endpoint: minio 服务地址, 例如 "localhost:9000"
        :param access_key: 访问密钥
        :param secret_key: 密钥
        :param secure: 是否使用 https
        """
        self.client = Minio(
            endpoint,
            access_key=access_key,
            secret_key=secret_key,
            secure=secure
        )

    def make_bucket(self, bucket_name):
        """
        新建桶，如果桶已存在则忽略
        :param bucket_name: 桶名称
        """
        found = self.client.bucket_exists(bucket_name)
        if not found:
            self.client.make_bucket(bucket_name)

    def upload_file(self, bucket_name, object_name, file_path):
        """
        上传单个文件到指定桶
        :param bucket_name: 桶名称`
        :param object_name: 存储在桶中的对象名
        :param file_path: 本地文件路径
        """
        if not self.client.bucket_exists(bucket_name):
            self.client.make_bucket(bucket_name)
        self.client.fput_object(
            bucket_name, object_name, file_path
        )

    def upload_folder(self, bucket_name, folder_path, prefix=""):
        """
        上传文件夹内所有文件到指定桶
        :param bucket_name: 桶名称
        :param folder_path: 本地文件夹路径
        :param prefix: 上传到桶中的前缀路径（可选）
        """
        if not self.client.bucket_exists(bucket_name):
            self.client.make_bucket(bucket_name)
        for root, dirs, files in os.walk(folder_path):
            for file in files:
                local_path = os.path.join(root, file)
                # 计算相对路径，作为object_name
                rel_path = os.path.relpath(local_path, folder_path)
                object_name = os.path.join(prefix, rel_path).replace("\\", "/")
                self.client.fput_object(
                    bucket_name, object_name, local_path
                )
    def remove_bucket(self, bucket_name):
        """
        删除指定的桶（即使桶非空也会删除）
        :param bucket_name: 桶名称
        """
        # 首先删除桶内所有对象
        import concurrent.futures
        objects_to_delete = list(self.client.list_objects(bucket_name, recursive=True))
        with concurrent.futures.ThreadPoolExecutor() as executor:
            list(tqdm(executor.map(lambda obj: self.client.remove_object(bucket_name, obj.object_name), objects_to_delete), total=len(objects_to_delete)))
        # 删除桶
        self.client.remove_bucket(bucket_name)
    def count_folders_in_bucket(self, bucket_name, prefix=""):
        """
        统计指定桶（及可选前缀）下有多少个文件夹（以'/'分隔的目录）
        :param bucket_name: 桶名称
        :param prefix: 前缀路径（可选）
        :return: 文件夹数量
        """
        if not self.client.bucket_exists(bucket_name):
            raise ValueError(f"桶 {bucket_name} 不存在")
        folder_set = set()
        objects = self.client.list_objects(bucket_name, prefix=prefix, recursive=True)
        for obj in objects:
            object_name = obj.object_name
            # 去除前缀
            if prefix and object_name.startswith(prefix):
                rel_path = object_name[len(prefix):]
                if rel_path.startswith("/"):
                    rel_path = rel_path[1:]
            else:
                rel_path = object_name
            # 分析路径中的文件夹
            parts = rel_path.split("/")
            # 只统计有文件夹的对象
            for i in range(1, len(parts)):
                folder_path = "/".join(parts[:i])
                folder_set.add(folder_path)
        return len(folder_set)
    def download_bucket(self, bucket_name, dest_folder, prefix="", num_workers=8):
        """
        下载指定桶（可选前缀）下的所有对象到本地文件夹，支持多线程和进度条
        :param bucket_name: 桶名称
        :param dest_folder: 本地目标文件夹
        :param prefix: 桶内前缀（可选）
        :param num_workers: 并发线程数
        """
        import concurrent.futures
        from tqdm import tqdm

        if not self.client.bucket_exists(bucket_name):
            raise ValueError(f"桶 {bucket_name} 不存在")
        # 获取所有对象列表
        objects = list(self.client.list_objects(bucket_name, prefix=prefix, recursive=True))
        total = len(objects)
        if total == 0:
            return

        def download_one(obj):
            object_name = obj.object_name
            # 计算本地保存路径
            rel_path = object_name[len(prefix):] if prefix and object_name.startswith(prefix) else object_name
            if rel_path.startswith("/"):
                rel_path = rel_path[1:]
            local_path = os.path.join(dest_folder, rel_path)
            local_dir = os.path.dirname(local_path)
            if not os.path.exists(local_dir):
                os.makedirs(local_dir, exist_ok=True)
            self.client.fget_object(bucket_name, object_name, local_path)

        with concurrent.futures.ThreadPoolExecutor(max_workers=num_workers) as executor:
            list(tqdm(executor.map(download_one, objects), total=total, desc="下载进度"))

if __name__ == "__main__":
    import os
    os.chdir("../")
    print("current directory: ", os.getcwd())
    minio_manager = MinioManager(
        endpoint="47.117.17.97:9000",
        access_key="hyf3513MINIO",
        secret_key="hyf3513MINIO",
        secure=False
    )
    minio_manager.remove_bucket("mini")
    # minio_manager.download_bucket("mini", "./results/mini","")