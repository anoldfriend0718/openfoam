import os.path as path
import sys
sys.path.append(path.dirname(path.abspath(__file__)))
import numpy as np
import pandas as pd




def read_min_max_field(file_path,sampling_rate,field):
    data=read_field_min_max_file(file_path)
    
    df=data[data["field"].str.contains(field)]
    df.reset_index(inplace=True,drop=True)

    df_sampling=df[df.index%sampling_rate==0]
    df_sampling.reset_index(inplace=True,drop=True)
    df_sampling.loc[:,"max"]=df_sampling.loc[:,"max"].astype(np.float64)
    df_sampling.loc[:,"min"]=df_sampling.loc[:,"min"].astype(np.float64)
    return df_sampling


def read_field_min_max_file(file_path):
    with open(file_path,"r") as fp:
        comment=fp.readline()
        header=fp.readline()
    header=header[1:-1].split()
    indexs_processor=[]
    for i,_ in enumerate(header):
        if header[i]=="processor":
            indexs_processor.append(i)
    indexs_processor.reverse()  

    data=pd.read_csv(file_path,comment='#', sep='\t',header=None)
    data=data.drop(indexs_processor,axis=1)
    data.rename(columns=lambda x:header[x],inplace=True)
    return data