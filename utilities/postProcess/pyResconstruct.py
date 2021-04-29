# To add a new cell, type '# %%'
# To add a new markdown cell, type '# %% [markdown]'
import os 
import sys 
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl

import os.path as path
sys.path.append(path.dirname(path.abspath(__file__)))
import polyMesh2d as mesh2d

# import sciPyFoam.polyMesh2d as mesh2d
import uuid
import traceback
import concurrent.futures
import json

def get_processor_dirs(caseDir):
    dirs=os.listdir(caseDir)
    procIndexs=[]
    procDirs=[]
    for dir in dirs:
        path=os.path.join(caseDir,dir)
        if(str.startswith(dir,"processor") and os.path.isdir(path) ):
            procIndex=dir.replace("processor",'')
            procIndexs.append(int(procIndex))
    procIndexs=np.sort(np.array(procIndexs))
    for procIndex in procIndexs:
        path=os.path.join(caseDir,f"processor{procIndex}")
        procDirs.append(path)
    return procDirs

def getTimes(caseDir):
    times=os.listdir(caseDir)
    timeDirs=[]
    for t in times:
        if(os.path.isdir(caseDir+'/'+t)):
            if(t.replace('.','',1).isdigit()):
                timeDirs.append(t)
    #         else:
    #             print(t,'is not a directory')
    timeDirs=np.array(timeDirs)
    times=np.array(timeDirs,dtype=float)
    ind=np.argsort(times)
    return timeDirs[ind],times[ind]

def reconstruct_impl(procDirs,timeName, fieldNames,patchName="frontAndBack"):
    dfProcGroup=[]
    for procDir in procDirs:
        meshData=mesh2d.getMesh(procDir, patchName)
        xProcs=meshData["x"]
        yProcs=meshData["y"]

        fields=mesh2d.readCellData_to_pointData(procDir, timeName, fieldNames, meshData)
        fields["pointData"]["x"]=xProcs
        fields["pointData"]["y"]=yProcs
        dfProc=pd.DataFrame(fields["pointData"])
        dfProcGroup.append(dfProc)
        
    df=pd.concat(dfProcGroup) 
    df=df.drop_duplicates(subset=["x","y"], keep="first")
    df=df.reset_index(drop=True)
    df.fillna(0,inplace=True)
    return df

def reconstruct(caseDir,timeName,fieldNames,patchName='frontAndBack'):
    procDirs=get_processor_dirs(caseDir)
    df=reconstruct_impl(procDirs,timeName,fieldNames,patchName)
    return df

def reconstruct_save(caseDir,timeName,fieldNames,saveFolderPath,patchName="frontAndBack",overWrite=False):
    try:
        savedFileName=f"{timeName}.csv"
        savedFilePath=os.path.join(saveFolderPath,savedFileName)
        if(os.path.exists(savedFilePath) and not overWrite):
            print(f"Time {timeName}: {timeName} result already exists")
            isSucceed=1
            return timeName,isSucceed
        df=reconstruct(caseDir,timeName,fieldNames,patchName)
        df.to_csv(savedFilePath,index=False)
        isSucceed=1
    except Exception as e:
         errMsgs=f"Time {timeName}:\n Unhandled exception happened: {e} with stack trace {traceback.format_exc()}\n"
         print(errMsgs)
         isSucceed=0

    return timeName,isSucceed

def reconstruct_list(caseDir,timeNames,fieldNames,worker=8,saveFolder="postProcess",patchName="frontAndBack",overWrite=False):
    saveFolderPath=os.path.join(caseDir,saveFolder)
    if not os.path.exists(saveFolderPath):
        os.mkdir(saveFolderPath)
    if not os.path.isdir(saveFolderPath):
        saveFolderPath="{saveFolderPath}_str(uuid.uuid4())[:8]"
        os.mkdir(saveFolderPath)
    print(f"saving folder path: {saveFolderPath}")

    futures=list()
    dict_status={"timeName":[],"isSucceed":[]}
    with concurrent.futures.ProcessPoolExecutor(max_workers=worker) as executor:
        print("start to map...")
        for timeName in timeNames:
            future=executor.submit(reconstruct_save,caseDir,timeName,fieldNames,saveFolderPath,patchName,overWrite)
            futures.append(future)
        print("start to reduce...")
        for _, future in enumerate(futures):
            returnedTimeName,isSucceed=future.result()
            dict_status["timeName"].append(returnedTimeName)
            dict_status["isSucceed"].append(isSucceed)
        df_status=pd.DataFrame(dict_status)

    df_succeed=df_status[df_status["isSucceed"]==1]
    df_failed=df_status[df_status["isSucceed"]==0]

    succeedNum=np.array(df_succeed["isSucceed"]).shape[0]
    print(f"succeed number: {succeedNum}")
    failedNum=np.array(df_failed["isSucceed"]).shape[0]
    print(f"failed number: {failedNum}")
    if failedNum:
        failedTimeName=df_failed["timeName"]
        print(f"failed time name:\n {failedTimeName.values}")

    return saveFolderPath



def reconstruct_all(caseDir,fieldNames,worker=8,saveFolder="postProcess",sampleRate=1,patchName="frontAndBack",overWrite=False):
    processor0Dir=os.path.join(caseDir,"processor0")
    if not os.path.exists(processor0Dir):
        print("Error: no processor dir exists!")
        return ""

    allTimeNames,_=getTimes(processor0Dir)
    allTimeNames=allTimeNames[allTimeNames!='0']
    indexs=np.arange(0,allTimeNames.shape[0])
    sampleTimeNames=allTimeNames[indexs%sampleRate==0]
    print(f"sample time names: {sampleTimeNames}")

    saveFolderPath=reconstruct_list(caseDir,sampleTimeNames,fieldNames,worker,saveFolder,patchName,overWrite)
    return saveFolderPath



if __name__ == "__main__":
    parser=argparse.ArgumentParser(description="pyResconstruct")
    parser.add_argument("-c",dest='case_dir',required=True,help='specify the test case directory')
    parser.add_argument("-f",dest='field_names',required=True,help='specify the field names')
    parser.add_argument("-t",dest='time_names',required=True,help='specify the time names')
    parser.add_argument("-w",action='store_true',help='if overwrite when the file has already exist')
    parser.add_argument("-s",dest='save_folder',default="postProcess",help="specify the save folder")
    parser.add_argument("-p",dest='patch',default="frontAndBack",help="specify the patch name")
    parser.add_argument("-r",dest='sample_rate',default=1,type=int,help="specify the sampling rate")
    parser.add_argument("-n",dest='worker_num',default=8,type=int,help='specify the worker num')

    args=parser.parse_args()
    caseDir=args.case_dir
    fieldNameTexts=args.field_names
    workerNum=args.worker_num
    saveFolder=args.save_folder
    sampleRate=args.sample_rate
    patchName=args.patch
    overWrite=args.w
    timeNames=args.time_names

    fieldNames=json.loads(fieldNameTexts)
    print(f"overwrite: {overWrite}")
    if timeNames=="all":
        reconstruct_all(caseDir,fieldNames,worker=workerNum,
                        saveFolder=saveFolder,sampleRate=sampleRate,
                        patchName=patchName,overWrite=overWrite)
    else:
        print(f"specified time names: {timeNames}")
        timeNames=json.loads(timeNames)
        reconstruct_list(caseDir,timeNames=timeNames,fieldNames=fieldNames,
                        worker=workerNum,saveFolder=saveFolder,
                        patchName=patchName,overWrite=overWrite)

