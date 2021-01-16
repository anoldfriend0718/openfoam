
import sys
sys.path.append(r"/home/anoldfriend/OpenFOAM/anoldfriend-7/utilities/")

import signal
import multiprocessing as mp
import time
from residual_monitor import read_residuals,plot_multiple_residuals,quit

log="run.log"
pressure_name="p_rgh"
nCorrectors=2
interval=5


sample_size=200
m_residuals=[["h","O2","CO2"],["Ux","Uy",pressure_name]]
m_thresholds=[[1e-1,1e-4,1e-5,1e-6,1e-7],[1e-1,1e-6,1e-7,1e-8]]
m_save_files=["residuals1.jpg","residuals2.jpg"]


def process_fun():
    line_offset=0
    iterations_offset=0
    while True:
        df,line_offset,iterations,info=read_residuals(log,line_offset,pressure_name,nCorrectors,sample_size)
        
        if "cum_physical_time" in info.keys():
            physical_time=info["cum_physical_time"]
        else:
            physical_time="not found"

        if "cum_execution_time" in info.keys():
            execution_time=info["cum_execution_time"]
        else:
            execution_time="not found"

        title=f"physical time : {physical_time} s, execution time : {execution_time} s"
        titles=[title]*len(m_residuals)
        
        if "latest_delta_time" in info.keys():
            delta_time=info["latest_delta_time"]
        else:
            delta_time= "not found"

        if "maxCo" in info.keys():
            maxCo=info["maxCo"]
        else:
            maxCo="not found"
        
        if "meanCo" in info.keys():
            meanCo=info["meanCo"]
        else:
            meanCo="not found"
        
        text=f"latest_delta_time: {delta_time} s \n" + \
             f"mean CFL num: {meanCo}\n" + \
             f"max CFL num: {maxCo}"
        texts=[text]*len(m_residuals)
        
        plot_multiple_residuals(df,iterations_offset,m_residuals,m_thresholds,titles,texts,m_save_files)
        iterations_offset+=iterations
        time.sleep(interval)


if __name__=="__main__":
    try:
        signal.signal(signal.SIGINT,quit)
        signal.signal(signal.SIGTERM,quit)

        p=mp.Process(target=process_fun)

        
        p.start()
        p.deamon=True
        
        while True:
            pass
    except Exception as err:
        print(f"Error Message: {err}")

