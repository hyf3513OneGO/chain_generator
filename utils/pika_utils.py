import json



class RMQConfig:
    def __init__(self,config_path:str):
        self.config_path = config_path
        with open(self.config_path,"r") as f:
            self.config = json.load(f)
            self.host = self.config["rmq"]["host"]
            self.port = self.config["rmq"]["port"]
            self.username = self.config["rmq"]["username"]
            self.password = self.config["rmq"]["password"]
            self.virtual_host = self.config["rmq"]["virtual_host"]
            