import os
import os.path as path
import sys
sys.path.append(path.dirname(path.abspath(__file__)))
import numpy as np
import pandas as pd
import concurrent.futures
import argparse
import json
import traceback
import tracemalloc
from functools import reduce
import pyFigure



def computeGasPhaseO2Conc(df):
    #df is a grouped dataframe
    df_gas=df[df["eps"]>1.0-1e-6]
    if(df_gas.shape[0]==0):
        return 0.0
    else:
        num=df_gas.shape[0]
        O2Conc=np.sum(df_gas["O2Conc"])/num
        return O2Conc

def computeAverageCokeFraction(df):
    num=df.shape[0]
    df_coke=df[df["coke"]>1e-6]
    if(df_coke.shape[0]==0):
        return 0.0
    else:
        average_coke_fraction=np.sum(df_coke["coke"])/num
        return average_coke_fraction

def computeAverageQdot(df):
    num=df.shape[0]
    df_coke=df[df["coke"]>1e-6]
    if(df_coke.shape[0]==0):
        return 0.0
    else:
        average_Qdot=np.sum(df_coke["Qdot"])/num
        return average_Qdot

def computeTransverselyAverages(df):
    df_group=df.groupby("x")
    df_meanT=df_group["T"].mean()
    df_meanT=df_meanT.reset_index()

    df_O2Conc=df_group.apply(computeGasPhaseO2Conc)
    df_O2Conc=df_O2Conc.reset_index(name = "O2Conc")

    df_mean_coke=df_group.apply(computeAverageCokeFraction)
    df_mean_coke=df_mean_coke.reset_index(name = "coke")

    df_mean_Qdot=df_group.apply(computeAverageQdot)
    df_mean_Qdot=df_mean_Qdot.reset_index(name = "Qdot")

    dfs = [df_meanT, df_O2Conc, df_mean_coke, df_mean_Qdot]

    df_combined=reduce(lambda df_left,df_right: pd.merge(df_left, df_right),dfs)
    return df_combined

def computeTransverselyAveragesAndSave(df,save_path):
    try:
        df_result=computeTransverselyAverages(df)
        df_result.to_csv(save_path,index=False)
    except Exception as e:
        errmsg=f"Unhandled exception happened: {e} with stack trace {traceback.format_exc()}"
        print(errmsg)
        return False
    return True


def readAndcomputeTransverselyAveragesAndSave(data_folder,time,save_folder):
    df=pyFigure.read_data_and_process(data_folder,time)
    save_path=os.path.join(save_folder,f"{time}.csv")
    print(f"save transversely averaged data to :{save_path}")
    ret=computeTransverselyAveragesAndSave(df,save_path)
    return ret

def batchComputeTransverselyAveragesForAll(data_folder,save_folder,worker_num=8):
    time_names=pyFigure.get_times_from_data_folder(data_folder)
    batchComputeTransverselyAverages(data_folder, save_folder, time_names, worker_num)

def batchComputeTransverselyAverages(data_folder, save_folder, time_names,worker_num=8):
    print(f"time names: {time_names}")
    print(f"save folders: {save_folder}")
    if not os.path.exists(save_folder):
        os.mkdir(save_folder)
    futures=[]
    results=[]
    with concurrent.futures.ProcessPoolExecutor(max_workers=worker_num) as executor:
        for time in time_names:
            future=executor.submit(readAndcomputeTransverselyAveragesAndSave,\
                                    data_folder,time,save_folder)
            futures.append(future)
        for _, future in enumerate(futures):
            ret=future.result()
            results.append(ret)
    results=np.array(results)
    print(f"processed time number: {results.shape[0]}, succeed number: {np.sum(results)}")

def computeMaxTemperatureAndOutletO2ConcHistory(min_max_file_path,transverse_data_folder):
    times=pyFigure.get_times_from_data_folder(transverse_data_folder)
    results=[]
    for time in times:
        df=pd.read_csv(f"{transverse_data_folder}/{time}.csv")
        xmax=np.max(df["x"])
        Tmax=np.max(df["T"])
        O2ConcAtOutlet=list((df[df["x"]==xmax])["O2Conc"])[0]
        ret={"Time":float(time),"Transverse_Tmax":Tmax,"O2ConcAtOutlet":O2ConcAtOutlet}
        results.append(ret)
    df_transverse_data=pd.DataFrame(results)

    df_min_max=pyFigure.read_min_max_field(min_max_file_path,1,"T")
    df_min_max=df_min_max[["Time","max"]]

    df_combined=pd.merge(df_min_max,df_transverse_data)

    save_folder=os.path.abspath(os.path.join(transverse_data_folder,"../others"))
    if not os.path.exists(save_folder):
        os.mkdir(save_folder)
    df_combined.to_csv(f"{save_folder}/MaxTemperatureAndOutletO2ConcHistory.csv",index=True)


    return df_combined

def compute_reaction_rate_burning_rate(data_folder,worker_num=8):
    times=pyFigure.get_times_from_data_folder(data_folder)
    cokeRectionRates=[]
    totalCokes=[]
    futures=[]
    with concurrent.futures.ProcessPoolExecutor(max_workers=worker_num) as executor:
        for time in times:
            future=executor.submit(compute_reaction_rate_burning_rate_for_one_time,data_folder, time)
            futures.append(future)
        for _, future in enumerate(futures):
            ret=future.result()
            cokeRectionRates.append(ret[0])
            totalCokes.append(ret[1])
        

    refCoke=totalCokes[0]
    burningRates=[1-coke/refCoke for coke in totalCokes]
    df_rate=pd.DataFrame({"time":[float(time) for time in times],"reaction_rate":cokeRectionRates,"burning_fraction":burningRates})
    
    save_folder=os.path.abspath(os.path.join(data_folder,"./others"))
    if not os.path.exists(save_folder):
        os.mkdir(save_folder)

    df_rate.to_csv(f"{save_folder}/ReactionRateAndBurningRate.csv",index=True)

    return df_rate

def compute_reaction_rate_burning_rate_for_one_time(data_folder, time):
    df=pyFigure.read_postProcess_csv_data(data_folder,time)
    num=df.shape[0]
    cokeRectionRate=np.sum(df["cokeRectionRate"])*-1/num
    totalCoke=np.sum(df["coke"])
    
    return (cokeRectionRate,totalCoke)

if __name__ == "__main__":

    # tracemalloc.start() # 开始跟踪内存分配
    try:
        parser=argparse.ArgumentParser(description="pyCompute")
        parser.add_argument("-d",dest='data_folder',required=True,help='specify the data folder')
        parser.add_argument("-s",dest='save_folder',required=True,help="specify the save folder")
        parser.add_argument("-t",dest='time_names',default="all",help='specify the time names')
        parser.add_argument("-n",dest='worker_num',default=8,type=int,help='specify the worker num')
        
        args=parser.parse_args()
        data_folder=args.data_folder
        worker_num=args.worker_num
        save_folder=args.save_folder
        time_names=args.time_names

        if time_names=="all":
            print(f"specified time names: all")
            batchComputeTransverselyAveragesForAll(data_folder,save_folder,worker_num=worker_num)
        else:
            print(f"specified time names: {time_names}")
            time_names=json.loads(time_names)
            batchComputeTransverselyAverages(data_folder,save_folder,time_names,worker_num=worker_num)

        print("succeed to  batch compute transversely averages")

        compute_reaction_rate_burning_rate(data_folder,worker_num)
        print("succeed to  compute reaction rate and burningrate")

    except Exception as e:
        errmsg=f"Unhandled exception happened: {e} with stack trace {traceback.format_exc()}"
        print(errmsg)
