
import signal
import multiprocessing as mp
import time
from utils import read_residuals,plot_multiple_residuals,quit

log="run.log"
pressure_name="p_rgh"
nCorrectors=2
interval=1


sample_size=200
m_residuals=[["e","Y"],["Ux","Uy",pressure_name]]
m_thresholds=[[1e-4,1e-5,1e-6],[1e-6,1e-8]]
m_save_files=["residuals1.jpg","residuals2.jpg"]


def process_fun():
    line_offset=0
    iterations_offset=0
    while True:
        df,line_offset,iterations,physical_time,execution_time=read_residuals(log,line_offset,pressure_name,nCorrectors,sample_size)
        title=f"physical time : {physical_time} s, execution time : {execution_time} s"
        titles=[title]*len(m_residuals)
        plot_multiple_residuals(df,iterations_offset,m_residuals,m_thresholds,titles,m_save_files)
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

