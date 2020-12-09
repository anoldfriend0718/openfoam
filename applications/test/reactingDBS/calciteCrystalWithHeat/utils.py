import matplotlib.pyplot as plt
import pandas as pd
import numpy as np 
import re

import sys
import time

from IPython.display import clear_output

def read_residuals(log_file,line_offset,pressure_name,nCorrectors,sample_size):
    # s1=datetime.now()
    with open(log_file,"r") as fp:
        fp.seek(line_offset)
        content = fp.readlines()
        line_offset=fp.tell()
    content = [x.strip() for x in content] 
    # e1=datetime.now()
    # print(f"reading time: {(e1-s1).total_seconds()} seconds")

    # s2=datetime.now()
    solver_term="Solving for "
    physical_time_term="Time = "
    execution_time_term="ExecutionTime = "
    objs=set()
    residuals=dict()
    physical_times=list()
    execution_times=list()
    for line in content:
        if re.search(solver_term, line):
            elements=line.split(" ")
            obj=elements[4].rstrip(',')
            residual=elements[8].rstrip(',')
            if obj not in objs:
                residuals[obj]=list()
                objs.add(obj)
            residuals[obj].append(residual)
        elif re.search(execution_time_term, line):
            elements=line.split(" ")
            time=elements[2]
            execution_times.append(time)
        elif re.search(physical_time_term, line):
            elements=line.split(" ")
            time=elements[2]
            physical_times.append(time)
        
    if len(residuals.keys())==0:
        print("Error Message: no residual was read...")
        return pd.DataFrame({}),line_offset,0,0,0

    if pressure_name in residuals.keys():
        residuals[pressure_name]=residuals[pressure_name][1::nCorrectors]
    if "rho" in residuals.keys():
        residuals.pop('rho', None)

    lengths=[len(residuals[k]) for k in residuals.keys()]
    min_length=min(lengths)
    iterations=max(lengths)

    x_start=max(min_length-sample_size,0)
    x_end=min_length
    indexes=np.arange(x_start,x_end)
    
    for obj in residuals.keys():
        residuals[obj]=residuals[obj][x_start:x_end]

    df_sample=pd.DataFrame(residuals,index=indexes).astype(float)

    return df_sample,line_offset,iterations,physical_times[-1],execution_times[-1]


def plot_residuals(df_sample,iterations_offset,residuals,thresholds,title,save_file):
    # clear_output(wait=True)
    indexes=[i+iterations_offset for i in df_sample.index.tolist()]
    if len(indexes)>0:
        x_start=indexes[0]
        x_end=indexes[-1]
        fig,ax=plt.subplots()
        ax.plot(indexes,df_sample.loc[:,residuals])
        x_lims=[x_start,x_end]
        for threshold in thresholds:
            ax.plot(x_lims,[threshold,threshold],"--")

        ax.set_yscale("log")
        ax.set_xlabel("Iterations")
        ax.set_ylabel("Residuals")
        ax.set_title(title)

        threshold_names=[str(value) for value in thresholds]
        legends=np.concatenate([np.array(residuals),threshold_names])
        ax.legend(legends,loc="upper left",bbox_to_anchor=((1.1,0.9)))

        fig.savefig(save_file, bbox_inches='tight')
        plt.clf()

    else:
        print("Error Message: no residual can be plot...")

def plot_multiple_residuals(df,iterations_offset,m_residuals,m_thresholds,titles,m_save_files):
    for i, residuals in enumerate(m_residuals):
        residuals=m_residuals[i]
        thresholds=m_thresholds[i]
        save_file=m_save_files[i]
        title=titles[i]
        plot_residuals(df,iterations_offset,residuals,thresholds,title,save_file)




def quit(signum,frame):
    sys.exit()
