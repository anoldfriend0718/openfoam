import logging
import os
from decorators import singleton

@singleton
class LoggerService(object):
    def __init__(self,case_dir,level=logging.INFO) -> None:
        super().__init__()
        self._logger=logging.getLogger()
        self._formatter=logging.Formatter("[%(asctime)s][pid:%(process)d][%(attr)s]:%(message)s")
        self._logger.setLevel(logging.DEBUG)
        
        self._log_file=os.path.join(case_dir,"workflow.log")
        # print(self._log_file)
        if os.path.exists(self._log_file):
            os.remove(self._log_file)
        self._create_file_logger(self._log_file,level)

    def get_logger(self,extra={}):
        if not "attr" in extra:
            extra={"attr":"manager"}
        return logging.LoggerAdapter(self._logger, extra)
    
    def _create_console_logger(self,level):
        consoleHandler=logging.StreamHandler()
        consoleHandler.setLevel(level)
        consoleHandler.setFormatter(self._formatter)
        self._logger.addHandler(consoleHandler)

    def _create_file_logger(self,log_file,level):
        file_handler=logging.FileHandler(log_file)
        file_handler.setLevel(level)
        file_handler.setFormatter(self._formatter)
        self._logger.addHandler(file_handler)