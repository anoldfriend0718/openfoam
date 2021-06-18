import traceback
import subprocess
import re
import os
import logging

class OFRunner(object):
    def __init__(self,logger_service)->None:
        super().__init__()
        self._logger_service=logger_service
        numeric_const_pattern = r"""
            [-+]? # optional sign
            (?:
                (?: \d* \. \d+ ) # .1 .12 .123 etc 9.1 etc 98.1 etc
                |
                (?: \d+ \.? ) # 1. 12. 123. etc 1 12 123 etc
            )
            # followed by optional exponent part if desired
            (?: [Ee] [+-]? \d+ ) ?
            """
        self._numeric_detecter = re.compile(numeric_const_pattern, re.VERBOSE)
    
    def run(self,case_dir,timeout=None):
        case_name=os.path.basename(case_dir)
        extra={"attr":case_name}
        logger = self._logger_service.get_logger(extra)
        logger.debug(f"change the dir to {case_dir}")

        os.chdir(case_dir)
        logger.debug("about to preProcess")
        try:
            subprocess.check_output("foamCleanTutorials 2>&1 > run.log", shell=True)
            subprocess.check_output("blockMesh 2>&1 >> run.log", shell=True)
            subprocess.check_output("decomposePar 2>&1 >> run.log", shell=True)
        except subprocess.CalledProcessError as exp:
            logger.error(f"Exception happened when preProcessing: {exp} with stack trace {traceback.format_exc()}")  #error 
            raise     
        logger.info("Succeed completing the preProcessing")

        logger.debug("about to trigger computation")        
        np=self._get_processor_num(case_dir,logger)
        worker=subprocess.Popen(f"mpirun -np {np} $FOAM_USER_APPBIN/icoDBSPermFoam -parallel 2>&1 >> run.log", shell=True)
        logger.info("wait for computation completion")
        try:
            worker.wait(timeout=timeout)
            if worker.returncode !=0:
                err_msg="Exception happended when computing: return code !=0"
                logger.error(err_msg)
                raise RuntimeError(err_msg)
        except subprocess.TimeoutExpired:
            err_msg="Time out when computing with limitation of {timeout} seconds"
            logger.error(err_msg)
            raise TimeoutError(err_msg)

        logger.info("about to check the computation result")
        with open("./run.log","r") as fp:
            lines=fp.readlines()
        if "Finalising parallel run" in lines[-1]:
            logger.info("Succeed completing the computing")
        else:
            err_msg = "Exception happended when computing. Please check run.log"
            logger.error(err_msg) 
            raise RuntimeError(err_msg)

        is_converged=False
        for i,line in enumerate(lines):
            if "Permeability computation is converged" in line: 
                logger.info("Computation is converged!")
                is_converged=True
                break 
        if not is_converged:
            raise RuntimeError("Computation is not converged")

        result_line=lines[i+1]
        perm=float(self._numeric_detecter.findall(result_line)[0])
        logger.info(f"Permeability: {perm} mD")

        return perm


    def _get_processor_num(self,case_dir,logger):
        logger.debug("about to get the processor num for computation")
        decomposeParDict_file_path=os.path.join(case_dir,"system/decomposeParDict")
        if not os.path.exists(decomposeParDict_file_path):
            err_msg=f"system/decomposeParDict is not exist in case dir: {case_dir}"
            logger.error(err_msg)
            raise ValueError(err_msg)

        with open(decomposeParDict_file_path,"r") as fp:
            lines=fp.readlines()

        is_get_np=False
        for line in lines:
            if "numberOfSubdomains" in line: 
                np=int(self._numeric_detecter.findall(line)[0])
                is_get_np=True
                break
        if not is_get_np:
            err_msg="numberOfSubdomains is not defined in system/decomposeParDict"
            logger.error(err_msg)
            raise ValueError(err_msg)

        logger.debug(f"Processor num for computation: {np}")
        return np


    


