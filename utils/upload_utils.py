from minio import Minio
from minio.error import S3Error
import os
import json

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
        objects_to_delete = self.client.list_objects(bucket_name, recursive=True)
        for obj in objects_to_delete:
            self.client.remove_object(bucket_name, obj.object_name)
        # 删除桶
        self.client.remove_bucket(bucket_name)
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
    # Minio 的 bucket 名称不能包含斜杠（/），即不能有“mini/test-1”这种双层目录结构。
    # 正确做法是 bucket 只用一级名称，比如 "mini"，然后通过 object_name 的前缀来实现“目录”效果。
    # 下面演示如何在 bucket "mini" 下上传到 "test-1/..." 这样的“目录”结构：

    minio_manager.make_bucket("mini")
    minio_manager.upload_folder("mini", "./results/mini/task_0_repeat-5-monomers", prefix="test-1")
    minio_manager.remove_bucket("mini")