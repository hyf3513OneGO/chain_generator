import argparse
import asyncio
import json
import uuid
from utils.chain_simulate_utils import HomoChainBuilder
from utils.chain_split_utils import extract_save_monomers
from utils.pika_utils import RMQConfig
from aio_pika import connect_robust, IncomingMessage, Message
from utils.log_utils import node_print,worker_print
from rdkit import Chem
from utils.upload_utils import MinioConfig, MinioManager
import os
import shutil
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
def process_relaxation_sync(msg: dict,config:ChainRelaxConfig,worker_id:str):
    task_id = msg["id"]
    prefix_folder = msg["prefix"]
    save_folder = os.path.join("results", prefix_folder, task_id)
    lmp_work_dir = os.path.join("lmp_dir",worker_id,config.work_dir)
    lmp_temp_dir = os.path.join("lmp_dir",worker_id,config.temp_dir)
    lmp_temp_save_dir = os.path.join("lmp_dir",worker_id,config.temp_save_dir)
    hcb = HomoChainBuilder(lmps_exec=config.lmps_exec,work_dir=lmp_work_dir,temp_dir=lmp_temp_dir)
    args_dict={
        'temp':msg["temp"],
        'high_temp':msg["high_temp"],
        'prev_nvt_steps':msg["prev_nvt_steps"],
        'cool_steps':msg["cool_steps"],
        'final_nvt_steps':msg["final_nvt_steps"],
        'time_step_low':msg["time_step_low"],
        'time_step':msg["time_step"],
        'box_length':msg["box_length"],
        'comm_cutoff':msg["comm_cutoff"]
    }
    psmiles = msg["psmiles"] # input with [*],[*] format
    homopoly_relaxed,homopoly_init = hcb.fromSmiles(
    psmiles = psmiles,
    n_repeat=msg.get("n_repeat",None),
    n_atoms=msg.get("n_atoms",None),
    save_dir=lmp_temp_save_dir,
    random_walk=msg["random_walk"],args_dict=args_dict
    )
    if not os.path.exists(save_folder):
        os.makedirs(save_folder)
    hcb.save_mol(homopoly_relaxed,os.path.join(save_folder,"homopoly_relaxed.sdf"))
    hcb.save_mol(homopoly_init,os.path.join(save_folder,"homopoly_init.sdf"))
    with open(os.path.join(save_folder,"psmiles.txt"),"w") as f:
        f.write(msg["psmiles"])
    
    homopoly_relaxed = Chem.RemoveHs(homopoly_relaxed)
    extract_save_monomers(mol=homopoly_relaxed,psmiles=psmiles,save_dir=save_folder)
    return task_id, msg
def run_task_in_subprocess(msg,config,worker_id:str):
    try:
        return process_relaxation_sync(msg,config,worker_id)
    except Exception as e:
        task_id = msg.get("id", "unknown")
        save_folder = os.path.join("results", msg.get("prefix", ""), task_id)
        if os.path.exists(save_folder):
            shutil.rmtree(save_folder)
        return ("error", msg, str(e))
async def handle_message(message: IncomingMessage, config: ChainRelaxConfig, channel,worker_id:str,minio_helper:MinioManager):
    async with message.process():
        try:
            msg = json.loads(message.body.decode())
            task_id = msg["id"]
            node_print(config.node, f"Chain Relax Task <{task_id}> started")

            loop = asyncio.get_running_loop()
            result = await loop.run_in_executor(None, run_task_in_subprocess, msg,config,worker_id)

            if result[0] == "error":
                raise RuntimeError(result[2])

            task_id, msg = result
            
            try:
                save_folder = os.path.join("results", msg.get("prefix", ""), task_id)
                minio_helper.upload_folder(msg.get("prefix"), save_folder, prefix=task_id)
                node_print(config.node, f"Task {task_id} uploaded to minio successfully.")
            except Exception as e:
                raise e
            node_print(config.node, f"Task {task_id} completed successfully.")

        except Exception as e:
            try:
                msg = json.loads(message.body.decode())
                task_id = msg.get("id", "unknown")
                save_folder = os.path.join("results", msg.get("prefix", ""), task_id)
                if os.path.exists(save_folder):
                    shutil.rmtree(save_folder)
            except Exception as e:
                raise e
            await channel.default_exchange.publish(
                Message(body=message.body, headers=message.headers, content_type="application/json"),
                routing_key=config.failed_queue_name
            )
            node_print(config.node, f"Task {task_id} failed with error: {e} and sent to failed queue.")

async def worker_main(worker_id: str, config_path: str):
    config = ChainRelaxConfig(config_path)
    rmq_config = RMQConfig(config_path)
    minio_config = MinioConfig(config_path)
    minio_helper = MinioManager(minio_config.endpoint, minio_config.access_key, minio_config.secret_key, minio_config.secure)
    worker_print(config.node, worker_id, f"[Starting ChainRelax with aio-pika")
    connection = await connect_robust(
        host=rmq_config.host,
        port=rmq_config.port,
        login=rmq_config.username,
        password=rmq_config.password,
        virtualhost=rmq_config.virtual_host,
        heartbeat=30
    )
    channel = await connection.channel()
    await channel.set_qos(prefetch_count=1)
    queue = await channel.declare_queue(config.queue_name, durable=True)
    lmp_work_dir = os.path.join("lmp_dir",worker_id,config.work_dir)
    lmp_temp_dir = os.path.join("lmp_dir",worker_id,config.temp_dir)
    lmp_temp_save_dir = os.path.join("lmp_dir",worker_id,config.temp_save_dir)
    if not os.path.exists(lmp_work_dir):
        os.makedirs(lmp_work_dir)
    if not os.path.exists(lmp_temp_dir):
        os.makedirs(lmp_temp_dir)
    if not os.path.exists(lmp_temp_save_dir):
        os.makedirs(lmp_temp_save_dir)
    await queue.consume(lambda msg: handle_message(msg, config, channel,worker_id,minio_helper), no_ack=False)

    try:
        while True:
            await asyncio.sleep(1)
    finally:
        await connection.close()
        worker_print(config.node, worker_id, f"[Worker {worker_id}] Connection closed.")
async def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--config', help='Path to configuration file', default='configs/config.json')
    parser.add_argument('--workers', help='Num for workers', default=1)
    args = parser.parse_args()

    tasks = [asyncio.create_task(worker_main(f"worker_{uuid.uuid4()}", args.config)) for i in range(int(args.workers))]

    try:
        await asyncio.gather(*tasks)
    except asyncio.CancelledError:
        pass


if __name__ == "__main__":
    try:
        asyncio.run(main())
    except KeyboardInterrupt:
        print("Shutdown requested")