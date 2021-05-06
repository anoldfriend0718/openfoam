# To add a new cell, type '# %%'
# To add a new markdown cell, type '# %% [markdown]'
import os
import os.path as path
import sys
sys.path.append(path.dirname(path.abspath(__file__)))
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import ticker
from matplotlib.ticker import FormatStrFormatter
from colored import fg, attr
import proplot as pplot # there are some nice colormaps in the proplot package
import concurrent.futures
import argparse
import json



# mpl.rcParams['font.family'] = 'Helvetica', #default font family
mpl.rcParams['mathtext.fontset'] = 'cm' #font for math
C_GREEN = fg('green')
C_RED = fg('red')
C_BLUE = fg('blue')
C_DEFAULT = attr('reset')

def use_paper_style(mpl, fontsize=9):
    mpl.rcParams['axes.titlesize']  = fontsize
    mpl.rcParams['axes.labelsize']  = fontsize
    mpl.rcParams['xtick.labelsize'] = fontsize
    mpl.rcParams['ytick.labelsize'] = fontsize
    mpl.rcParams['legend.fontsize'] = fontsize
    mpl.rcParams['axes.titlesize']  = fontsize
    mpl.rcParams['axes.titlesize']  = fontsize

def read_postProcess_csv_data(postProcessDir,timeName):
    fileName=f"{timeName}.csv"
    filePath=os.path.join(postProcessDir,fileName)
    df=pd.read_csv(filePath)
    return df

def tickformatter():
    formatter = ticker.ScalarFormatter(useMathText=True)
    formatter.set_scientific(True) 
    formatter.set_powerlimits((-1,1)) 
    return formatter

def time_str(time_second):
    Scale,name_time=1,'seconds'
    if(time_second>(86400*365)):
        Scale = 3.17e-08
        name_time='years'
    elif(time_second>86400):
        Scale = 1.1574074074074073e-05
        name_time='days'
    elif(time_second>3600):
        Scale = 0.0002777777777777778
        name_time='hours'
    elif(time_second>60):
        Scale = 0.016666666666666666
        name_time='minutes'
    str_time=str('%.3f %s'%(time_second*Scale,name_time))
    return str_time

def read_data_and_process(folder,timeName):
    df=read_postProcess_csv_data(folder,timeName)
    df["UNorm"]=np.sqrt(df['U_0']**2 + df['U_1']**2 + df['U_2']**2)
    MO2=0.032 #g/mol
    df["O2Conc"]=df['rho']*df['eps']*df['O2']/MO2
    MCO2=0.044 #g/mol
    df["CO2Conc"]=df['rho']*df['eps']*df['CO2']/MCO2
    return df



def figsize_cm(w_cm,x,y,w_offset_cm=0.0):
    #Get figure size in unit cm. 
    cm2inch=1/2.54
    w=w_cm*cm2inch

    length_x=x.max()-x.min()
    length_y=y.max()-y.min()
    xyRatio=length_y/length_x
    h=w*xyRatio
    w=w+w_offset_cm*cm2inch
    return (w,h)


def plot_contourf_Impl(X,Y,Z,title,label,cmap=pplot.Colormap('CoolWarm'),levels=250,figwidth=20,vmin=0,vmax=0):
    figsize=figsize_cm(figwidth,X,Y)
    fig=plt.figure(figsize=figsize)
    ax=plt.gca()
    # ax.axis('scaled')
    ax.set_xlim(X.min(),X.max())
    ax.set_ylim(Y.min(),Y.max())
    ax.set_xlabel('X (m)')
    ax.set_ylabel('Y (m)')
    ax.set_title(title)
    formatter=tickformatter()
    ax.xaxis.set_major_formatter(formatter) 
    ax.yaxis.set_major_formatter(formatter) 
    Xi,Yi=np.meshgrid(X, Y)
    if vmin==0:
        vmin=np.min(Z)
    if vmax==0:
        vmax=np.max(Z)
    CS=ax.contourf(Xi,Yi,Z, cmap=cmap, levels=levels,vmin=vmin,vmax=vmax)
    ax_cb = ax.inset_axes([1.04, 0, 0.02,1])
    plt.colorbar(CS,cax=ax_cb,label=label)
    plt.tight_layout()

# def plot_contourf(dfpivot,label,cmap=pplot.Colormap('CoolWarm'),levels=250):
#     X=dfpivot.columns.values
#     Y=dfpivot.index.values
#     Z=dfpivot.values
#     plot_contourf_Impl(X,Y,Z,label,cmap=cmap,levels=levels)

def plot_contourf(df,fieldName,title,label,cmap=pplot.Colormap('CoolWarm'),levels=250,figwidth=20,vmin=0,vmax=0):
    dfpivot=df.pivot("y", "x", fieldName)
    X=dfpivot.columns.values
    Y=dfpivot.index.values
    Z=dfpivot.values
    plot_contourf_Impl(X,Y,Z,title,label,cmap=cmap,levels=levels,figwidth=figwidth,vmin=0,vmax=0)

def plot_contourf_save(df,fieldName,title,label,folder_path,cmap=pplot.Colormap('CoolWarm'),levels=250,figwidth=20,vmin=0,vmax=0,dpi=600):
    dfpivot=df.pivot("y", "x", fieldName)
    X=dfpivot.columns.values
    Y=dfpivot.index.values
    Z=dfpivot.values
    plot_contourf_Impl(X,Y,Z,title,label,cmap=cmap,levels=levels,figwidth=figwidth,vmin=0,vmax=0)
    if not os.path.exists(folder_path):
        os.mkdir(folder_path)
    plt.savefig(f"{folder_path}/{title}.jpg".replace(" ",""),dpi=dpi)

def plot_multiple_contourf_save(df,fields,time_instant,save_folder,xranges={}):
    if "eps" in fields:
        eps_title=f"porosity contour at {time_str(time_instant)}"
        if "eps" in xranges:
            vmin=xranges["eps"]["vmin"]
            vmax=xranges["eps"]["vmax"]
        else:
            vmin=0
            vmax=0
        plot_contourf_save(df,"eps",eps_title,label='porosity',folder_path=save_folder,vmin=vmin,vmax=vmax)

    if "UNorm" in fields:
        unorm_title=f"velocity magnitude contour at {time_str(time_instant)}"
        if "UNorm" in xranges:
            vmin=xranges["UNorm"]["vmin"]
            vmax=xranges["UNorm"]["vmax"]
        else:
            vmin=0
            vmax=0
        plot_contourf_save(df,"UNorm",unorm_title,label='velocity magnitude',folder_path=save_folder,vmin=vmin,vmax=vmax)

    if "O2Conc" in fields:
        O2Conc_title=f"O$_{2}$ concentration contour at {time_str(time_instant)}"
        if "O2Conc" in xranges:
            vmin=xranges["O2Conc"]["vmin"]
            vmax=xranges["O2Conc"]["vmax"]
        else:
            vmin=0
            vmax=0
        plot_contourf_save(df,"O2Conc",O2Conc_title,label='O$_2$ mole concentration (mol/m$^3$)',folder_path=save_folder,vmin=vmin,vmax=vmax)

    if "T" in fields:
        T_title=f"Temperature contour at {time_str(time_instant)}"
        if "T" in xranges:
            vmin=xranges["T"]["vmin"]
            vmax=xranges["T"]["vmax"]
        else:
            vmin=0
            vmax=0
        plot_contourf_save(df,"T",T_title,label='Temperature ($^{\circ}$C)',folder_path=save_folder,vmin=vmin,vmax=vmax)

    if "Qdot" in fields:
        Qdot_title=f"Reaction Heat Rate contour at {time_str(time_instant)}"
        if "Qdot" in xranges:
            vmin=xranges["Qdot"]["vmin"]
            vmax=xranges["Qdot"]["vmax"]
        else:
            vmin=0
            vmax=0
        plot_contourf_save(df,"Qdot",Qdot_title,label='Reaction Heat Rate (J/(m$^3\cdot$s))',folder_path=save_folder,vmin=vmin,vmax=vmax)

def read_plot_multiple_contourf_save(data_folder,fields,time_instant,save_folder,xranges={}):
    df=read_data_and_process(data_folder,time_instant)
    plot_multiple_contourf_save(df,fields,float(time_instant),save_folder,xranges)


def plot_multiple_contourf_save_all_time(fields,data_folder,save_folder,xranges={},time_names=[],worker_num=10):
    if(len(time_names)==0):
        time_names=get_times_from_data_folder(data_folder)
    futures=[]
    with concurrent.futures.ProcessPoolExecutor(max_workers=worker_num) as executor:
        for time in time_names:
            future=executor.submit(read_plot_multiple_contourf_save,data_folder,fields,time,save_folder,xranges)
            futures.append(future)
        for _, future in enumerate(futures):
            future.result()
        
def get_times_from_data_folder(dataFolder):
    times=os.listdir(dataFolder)
    timeNames=[]
    for t in times:
        if(str.endswith(t,".csv")):
            timeNames.append(t.replace(".csv",""))
    timeNames=np.sort(timeNames)
    return timeNames

def plot_temperature_and_coke_evolution(dataFolder,timeNames=[],workerNum=10):
    if(len(timeNames)==0):
        timeNames=get_times_from_data_folder(dataFolder)
    futures=[]
    maxTemperatures=[]
    cokeFractions=[]
    with concurrent.futures.ThreadPoolExecutor(max_workers=workerNum) as executor:
        for timeName in timeNames:
            future=executor.submit(read_postProcess_csv_data,dataFolder,timeName)
            futures.append(future)
        for _, future in enumerate(futures):
            df=future.result()
            maxTemperatures.append(np.max(df["T"]))
            cokeFractions.append(np.mean(df["coke"]))

    cokeFractions=np.array(cokeFractions)
    maxTemperatures=np.array(maxTemperatures)

    with plt.style.context(['science','no-latex', 'std-colors']):
        fig, ax = plt.subplots(figsize=(6,4))

        ax.set_xlabel(f"Time (s)")
        
        ax.set_title(f"Temporal Evolution of Maximum Temperature and Coke Fraction",color="k")
        times=[float(timeStr) for timeStr in timeNames]
        ax.plot(times,maxTemperatures,linestyle="-",label="Maximum Temperature",color="b")
        ax.set_ylabel("Maximum combustion Temperature ($^{\circ}$C)",color="b")
        ax.tick_params(axis='y', labelcolor="b")
        ax.set_xlim([0,np.max(times)])
    
        ax2 = ax.twinx()
        ax2.plot(times,cokeFractions*100,linestyle="-",color="r")
        ax2.set_xlabel('Time (s)',color="r")
        ax2.set_ylabel("Residual coke fraction (%)",color="r")
        ax2.tick_params(axis='y', labelcolor="r")
        ax2.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))

        fig.tight_layout()
        ax.autoscale(tight=True)
        ax2.autoscale(tight=True)


def read_field_min_max_file(file_path,sampling_rate):
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
    data_sampling=data[data.index%sampling_rate==0]
    return data_sampling

def plot_max_temperature(file_path,sampling_rate):
    data_sampling=read_field_min_max_file(file_path,sampling_rate)
    fig, ax = plt.subplots()
    ax.plot(data_sampling["Time"],data_sampling["max"],label="max temperature")
    ax.set_xlabel(f"Time (s)")
    ax.set_ylabel("Temperature  ($^{\circ}$C)")
    ax.set_title("Temporal Evolution of Minimum/Maximum Temperature")
    ax.legend()
    ax.autoscale(tight=True)
    plt.tight_layout()

if __name__ == "__main__":
    parser=argparse.ArgumentParser(description="pyResconstruct")
    parser.add_argument("-f",dest='field_names',required=True,help='specify the field names')
    parser.add_argument("-d",dest='data_folder',required=True,help='specify the data folder')
    parser.add_argument("-s",dest='save_folder',required=True,help="specify the save folder")
    parser.add_argument("-x",dest='xranges',default="{}",help='specify the color bar ranges')
    parser.add_argument("-t",dest='time_names',default="all",help='specify the time names')
    parser.add_argument("-n",dest='worker_num',default=8,type=int,help='specify the worker num')
    
    args=parser.parse_args()
    data_folder=args.data_folder
    field_names_texts=args.field_names
    field_names=json.loads(field_names_texts)
    xranges_texts=args.xranges
    xranges=json.loads(xranges_texts)
    worker_num=args.worker_num
    save_folder=args.save_folder
    time_names=args.time_names
    if time_names=="all":
        plot_multiple_contourf_save_all_time(field_names,data_folder,save_folder,xranges,time_names=[],worker_num=worker_num)
    else:
        print(f"specified time names: {time_names}")
        time_names=json.loads(time_names)
        plot_multiple_contourf_save_all_time(field_names,data_folder,save_folder,xranges,time_names,worker_num=worker_num)

    print("succeed to plot images")
