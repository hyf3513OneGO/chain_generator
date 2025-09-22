import argparse
import asyncio
import json
import uuid
import time
import os
import shutil
from functools import partial

from rdkit import Chem
from aio_pika import connect_robust, IncomingMessage, Message
from aio_pika.exceptions import AMQPConnectionError, ChannelInvalidStateError
import aiormq

from utils.chain_simulate_utils import HomoChainBuilder
from utils.chain_split_utils import extract_save_monomers
from utils.pika_utils import RMQConfig
from utils.log_utils import node_print, worker_print
from utils.upload_utils import MinioConfig, MinioManager


# =========================
# Configs
# =========================
class ChainRelaxConfig:
    def __init__(self, config_path: str) -> None:
        with open(config_path, "r") as f:
            config = json.load(f)
            self.node = config["chainRelaxation"]["node"]
            self.queue_name = config["chainRelaxation"]["queue_name"]
            self.failed_queue_name = config["chainRelaxation"].get("failed_queue_name", "chainRelaxation_failed")
            self.lmps_exec = config["chainRelaxation"]["lmps_exec"]
            self.work_dir = config["chainRelaxation"]["work_dir"]
            self.temp_dir = config["chainRelaxation"]["temp_dir"]
            self.temp_save_dir = config["chainRelaxation"]["temp_save_dir"]
            self.gpu = config["chainRelaxation"]["gpu"]


# =========================
# Connection Manager
# =========================
class ConnectionManager:
    def __init__(self, rmq_config: RMQConfig, max_retries: int = 5, retry_delay: float = 5.0):
        self.rmq_config = rmq_config
        self.max_retries = max_retries
        self.retry_delay = retry_delay
        self.connection = None
        self.channel = None
        self.queue = None

    async def connect(self):
        """建立连接并重试机制"""
        for attempt in range(self.max_retries):
            try:
                self.connection = await connect_robust(
                    host=self.rmq_config.host,
                    port=self.rmq_config.port,
                    login=self.rmq_config.username,
                    password=self.rmq_config.password,
                    virtualhost=self.rmq_config.virtual_host,
                    heartbeat=180,  # 提高心跳，降低误判
                    client_properties={"connection_name": f"chain_relax_worker_{uuid.uuid4()}"}
                )
                self.channel = await self.connection.channel()
                await self.channel.set_qos(prefetch_count=1)
                return True
            except (AMQPConnectionError, ChannelInvalidStateError, Exception) as e:
                if attempt < self.max_retries - 1:
                    worker_print("ConnectionManager", "main", f"连接失败，{self.retry_delay}秒后重试 (尝试 {attempt + 1}/{self.max_retries}): {e}")
                    await asyncio.sleep(self.retry_delay)
                else:
                    worker_print("ConnectionManager", "main", f"连接失败，已达到最大重试次数: {e}")
                    raise e
        return False

    async def setup_queue(self, queue_name: str):
        """设置主队列"""
        if not self.channel:
            raise RuntimeError("Channel not initialized")
        # 如需单活消费者，可加 arguments={"x-single-active-consumer": True}
        self.queue = await self.channel.declare_queue(queue_name, durable=True)
        return self.queue

    async def ensure_failed_queue(self, queue_name: str):
        """确保失败队列存在（默认交换机按 routing_key 投递到同名队列）"""
        if not self.channel:
            raise RuntimeError("Channel not initialized")
        await self.channel.declare_queue(queue_name, durable=True)

    async def close(self):
        """关闭连接"""
        try:
            if self.connection and not self.connection.is_closed:
                await self.connection.close()
        except Exception as e:
            worker_print("ConnectionManager", "main", f"关闭连接时出错: {e}")

    async def is_connected(self):
        """检查连接状态"""
        return (
            self.connection and not self.connection.is_closed and
            self.channel and not self.channel.is_closed
        )

    async def reconnect(self):
        """重新连接"""
        await self.close()
        await asyncio.sleep(1)  # 等待一秒再重连
        return await self.connect()


# =========================
# 计算与文件处理
# =========================
def process_relaxation_sync(msg: dict, config: ChainRelaxConfig, worker_id: str):
    task_id = msg["id"]
    prefix_folder = msg["prefix"]
    save_folder = os.path.join("results", prefix_folder, task_id)

    lmp_work_dir = os.path.join("lmp_dir", worker_id, config.work_dir)
    lmp_temp_dir = os.path.join("lmp_dir", worker_id, config.temp_dir)
    lmp_temp_save_dir = os.path.join("lmp_dir", worker_id, config.temp_save_dir)

    hcb = HomoChainBuilder(lmps_exec=config.lmps_exec, work_dir=lmp_work_dir, temp_dir=lmp_temp_dir, conf_mm_gpu=config.gpu)

    args_dict = {
        'temp': msg["temp"],
        'high_temp': msg["high_temp"],
        'prev_nvt_steps': msg["prev_nvt_steps"],
        'cool_steps': msg["cool_steps"],
        'final_nvt_steps': msg["final_nvt_steps"],
        'time_step_low': msg["time_step_low"],
        'time_step': msg["time_step"],
        'box_length': msg["box_length"],
        'comm_cutoff': msg["comm_cutoff"]
    }

    psmiles = msg["psmiles"]  # input with [*],[*] format
    homopoly_relaxed, homopoly_init = hcb.fromSmiles(
        psmiles=psmiles,
        n_repeat=msg.get("n_repeat", None),
        n_atoms=msg.get("n_atoms", None),
        save_dir=lmp_temp_save_dir,
        random_walk=msg["random_walk"], args_dict=args_dict
    )

    if not os.path.exists(save_folder):
        os.makedirs(save_folder)

    hcb.save_mol(homopoly_relaxed, os.path.join(save_folder, "homopoly_relaxed.sdf"))
    hcb.save_mol(homopoly_init, os.path.join(save_folder, "homopoly_init.sdf"))
    with open(os.path.join(save_folder, "psmiles.txt"), "w") as f:
        f.write(msg["psmiles"])

    homopoly_relaxed = Chem.RemoveHs(homopoly_relaxed)
    extract_save_monomers(mol=homopoly_relaxed, psmiles=psmiles, save_dir=save_folder)
    return task_id, msg


def run_task_in_subprocess(msg, config, worker_id: str):
    try:
        return process_relaxation_sync(msg, config, worker_id)
    except Exception as e:
        task_id = msg.get("id", "unknown")
        save_folder = os.path.join("results", msg.get("prefix", ""), task_id)
        if os.path.exists(save_folder):
            shutil.rmtree(save_folder)
        return ("error", msg, str(e))


# =========================
# 安全确认封装
# =========================
async def safe_settle(message: IncomingMessage, action: str = "ack", requeue: bool = True) -> str:
    """
    对 ack/nack 做最佳努力，不重试、不抛异常。
    返回: 'done' | 'skipped-closed' | 'failed-closing'
    """
    try:
        ch = message.channel
        if not ch or ch.is_closed:
            return "skipped-closed"
        if action == "ack":
            await message.ack()
        else:
            await message.nack(requeue=requeue)
        return "done"
    except (aiormq.exceptions.ChannelInvalidStateError, ChannelInvalidStateError, AMQPConnectionError):
        return "failed-closing"


# =========================
# 消息处理
# =========================
async def handle_message(message: IncomingMessage, config: ChainRelaxConfig, conn_manager: ConnectionManager, worker_id: str, minio_helper: MinioManager):
    # 提前解析，解析失败直接丢失败队列并 ack
    try:
        msg = json.loads(message.body.decode())
    except Exception:
        try:
            if await conn_manager.is_connected():
                await conn_manager.ensure_failed_queue(config.failed_queue_name)
                await conn_manager.channel.default_exchange.publish(
                    Message(body=message.body, headers=message.headers, content_type="application/json"),
                    routing_key=config.failed_queue_name
                )
        finally:
            await safe_settle(message, "ack")
        return

    task_id = msg.get("id", "unknown")
    node_print(config.node, f"Chain Relax Task <{task_id}> started")

    loop = asyncio.get_running_loop()

    try:
        # 1) 计算密集：线程池/进程池。当前为线程池；如需更稳可切换 ProcessPoolExecutor。
        result = await loop.run_in_executor(None, run_task_in_subprocess, msg, config, worker_id)
        if result[0] == "error":
            raise RuntimeError(result[2])

        # 2) 阻塞 I/O（MinIO 上传）放线程池，避免阻塞事件循环与心跳
        save_folder = os.path.join("results", msg.get("prefix", ""), task_id)
        upload_fn = partial(minio_helper.upload_folder, msg.get("prefix"), save_folder, task_id)
        await loop.run_in_executor(None, upload_fn)
        node_print(config.node, f"Task {task_id} uploaded to minio successfully.")
        node_print(config.node, f"Task {task_id} completed successfully.")

        # 3) 成功路径 ack（最佳努力）
        await safe_settle(message, "ack")

        # 如需幂等标记，可在此处最后一步落 "_DONE" 文件或写数据库

    except asyncio.CancelledError:
        # 取消：尽量 nack 重入队；清理解算产物也放线程池
        await safe_settle(message, "nack", requeue=True)
        try:
            save_folder = os.path.join("results", msg.get("prefix", ""), task_id)
            if os.path.exists(save_folder):
                await loop.run_in_executor(None, shutil.rmtree, save_folder)
        except Exception:
            pass
        raise

    except Exception as e:
        # 失败：清理 → 发失败队列（重连一次）→ nack 重入队（最佳努力）
        try:
            save_folder = os.path.join("results", msg.get("prefix", ""), task_id)
            if os.path.exists(save_folder):
                await loop.run_in_executor(None, shutil.rmtree, save_folder)
        except Exception as cleanup_error:
            node_print(config.node, f"清理任务文件夹时出错: {cleanup_error}")

        # 发送失败消息
        try:
            if await conn_manager.is_connected():
                await conn_manager.ensure_failed_queue(config.failed_queue_name)
                await conn_manager.channel.default_exchange.publish(
                    Message(body=message.body, headers=message.headers, content_type="application/json"),
                    routing_key=config.failed_queue_name
                )
                node_print(config.node, f"Task {task_id} failed with error: {e} -> sent to failed queue.")
            else:
                node_print(config.node, f"连接已断开，无法发送失败消息。任务 {task_id} 失败: {e}")
        except (ChannelInvalidStateError, AMQPConnectionError) as publish_error:
            node_print(config.node, f"发送失败消息时出错: {publish_error}，尝试重连...")
            try:
                await conn_manager.reconnect()
                await conn_manager.setup_queue(config.queue_name)
                await conn_manager.ensure_failed_queue(config.failed_queue_name)
                await conn_manager.channel.default_exchange.publish(
                    Message(body=message.body, headers=message.headers, content_type="application/json"),
                    routing_key=config.failed_queue_name
                )
                node_print(config.node, f"重连后成功发送失败消息。任务 {task_id} 失败: {e}")
            except Exception as retry_error:
                node_print(config.node, f"重连失败，无法发送失败消息。任务 {task_id} 失败: {e}，重连错误: {retry_error}")
        except Exception as publish_error:
            node_print(config.node, f"发送失败消息时出现未知错误: {publish_error}，任务 {task_id} 失败: {e}")

        # 重入队（最佳努力）
        await safe_settle(message, "nack", requeue=True)


# =========================
# Worker 主体
# =========================
async def worker_main(worker_id: str, config_path: str):
    config = ChainRelaxConfig(config_path)
    rmq_config = RMQConfig(config_path)
    minio_config = MinioConfig(config_path)
    minio_helper = MinioManager(minio_config.endpoint, minio_config.access_key, minio_config.secret_key, minio_config.secure)

    worker_print(config.node, worker_id, f"[Starting ChainRelax with aio-pika")

    # 连接管理器
    conn_manager = ConnectionManager(rmq_config)

    # 建立连接与队列
    try:
        await conn_manager.connect()
        await conn_manager.setup_queue(config.queue_name)
        await conn_manager.ensure_failed_queue(config.failed_queue_name)
        worker_print(config.node, worker_id, f"成功连接到 RabbitMQ")
    except Exception as e:
        worker_print(config.node, worker_id, f"连接失败: {e}")
        return

    # 工作目录
    lmp_work_dir = os.path.join("lmp_dir", worker_id, config.work_dir)
    lmp_temp_dir = os.path.join("lmp_dir", worker_id, config.temp_dir)
    lmp_temp_save_dir = os.path.join("lmp_dir", worker_id, config.temp_save_dir)
    for dp in [lmp_work_dir, lmp_temp_dir, lmp_temp_save_dir]:
        os.makedirs(dp, exist_ok=True)

    in_flight = set()
    consumer_tag = None

    async def _wrapped_handler(msg: IncomingMessage):
        # 包装以跟踪在途任务，吞掉取消避免 aio_pika.tools._task_done 噪声
        coro = handle_message(msg, config, conn_manager, worker_id, minio_helper)
        task = asyncio.create_task(coro)
        in_flight.add(task)
        try:
            await task
        except asyncio.CancelledError:
            # 已在 handle_message 内完成 nack/清理
            return
        finally:
            in_flight.discard(task)

    # 启动消费
    consumer_tag = await conn_manager.queue.consume(_wrapped_handler, no_ack=False)
    worker_print(config.node, worker_id, f"开始监听队列: {config.queue_name}")

    # 健康检查循环
    last_health_check = time.time()
    health_check_interval = 30

    try:
        while True:
            await asyncio.sleep(1)
            now = time.time()
            if now - last_health_check > health_check_interval:
                if not await conn_manager.is_connected():
                    worker_print(config.node, worker_id, f"检测到连接断开，尝试重连...")
                    try:
                        await conn_manager.reconnect()
                        await conn_manager.setup_queue(config.queue_name)
                        await conn_manager.ensure_failed_queue(config.failed_queue_name)
                        consumer_tag = await conn_manager.queue.consume(_wrapped_handler, no_ack=False)
                        worker_print(config.node, worker_id, f"重连成功并恢复消费")
                    except Exception as reconnect_error:
                        worker_print(config.node, worker_id, f"重连失败: {reconnect_error}")
                last_health_check = now

    except KeyboardInterrupt:
        worker_print(config.node, worker_id, f"收到中断信号，正在关闭...")
    except Exception as e:
        worker_print(config.node, worker_id, f"Worker 出现未处理的错误: {e}")
    finally:
        # 优雅关停：先停消费，再等在途任务完成，最后关连接
        try:
            if conn_manager.channel and not conn_manager.channel.is_closed and consumer_tag:
                await conn_manager.channel.basic_cancel(consumer_tag)
        except Exception:
            pass
        try:
            if in_flight:
                await asyncio.wait(in_flight, timeout=60)
        except Exception:
            pass
        await conn_manager.close()
        worker_print(config.node, worker_id, f"[Worker {worker_id}] Connection closed.")


# =========================
# main
# =========================
async def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--config', help='Path to configuration file', default='configs/config.json')
    parser.add_argument('--workers', help='Num for workers', default=1)
    args = parser.parse_args()

    tasks = [asyncio.create_task(worker_main(f"worker_{uuid.uuid4()}", args.config)) for _ in range(int(args.workers))]

    try:
        await asyncio.gather(*tasks)
    except asyncio.CancelledError:
        pass


if __name__ == "__main__":
    try:
        asyncio.run(main())
    except KeyboardInterrupt:
        print("Shutdown requested")
