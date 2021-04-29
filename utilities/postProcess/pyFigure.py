# To add a new cell, type '# %%'
# To add a new markdown cell, type '# %% [markdown]'
import os
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import ticker
from matplotlib.ticker import FormatStrFormatter
from colored import fg, attr
import proplot as pplot # there are some nice colormaps in the proplot package
import concurrent.futures


# mpl.rcParams['font.family'] = 'Helvetica', #default font family
mpl.rcParams['mathtext.fontset'] = 'cm' #font for math
C_GREEN = fg('green')
C_RED = fg('red')
C_BLUE = fg('blue')
C_DEFAULT = attr('reset')

def usePaperStyle(mpl, fontsize=9):
    mpl.rcParams['axes.titlesize']  = fontsize
    mpl.rcParams['axes.labelsize']  = fontsize
    mpl.rcParams['xtick.labelsize'] = fontsize
    mpl.rcParams['ytick.labelsize'] = fontsize
    mpl.rcParams['legend.fontsize'] = fontsize
    mpl.rcParams['axes.titlesize']  = fontsize
    mpl.rcParams['axes.titlesize']  = fontsize

def readData(postProcessDir,timeName):
    fileName=f"{timeName}.csv"
    filePath=os.path.join(postProcessDir,fileName)
    df=pd.read_csv(filePath)
    return df

def tickformatter():
    formatter = ticker.ScalarFormatter(useMathText=True)
    formatter.set_scientific(True) 
    formatter.set_powerlimits((-1,1)) 
    return formatter

def timeStr(time_second):
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
    str_time=str('%.2f %s'%(time_second*Scale,name_time))
    return str_time


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


def plot_contourf_Impl(X,Y,Z,label,cmap=pplot.Colormap('CoolWarm'),levels=250,figwidth=20):
    figsize=figsize_cm(figwidth,X,Y)
    fig=plt.figure(figsize=figsize)
    ax=plt.gca()
    ax.axis('scaled')
    ax.set_xlim(X.min(),X.max())
    ax.set_ylim(Y.min(),Y.max())
    ax.set_xlabel('X (m)')
    ax.set_ylabel('Y (m)')
    formatter=tickformatter()
    ax.xaxis.set_major_formatter(formatter) 
    ax.yaxis.set_major_formatter(formatter) 
    Xi,Yi=np.meshgrid(X, Y)
    CS=ax.contourf(Xi,Yi,Z, cmap=cmap, levels=levels)
    ax.axis('scaled')
    ax_cb = ax.inset_axes([1.04, 0, 0.02,1])
    plt.colorbar(CS,cax=ax_cb,label=label)
    plt.tight_layout()

# def plot_contourf(dfpivot,label,cmap=pplot.Colormap('CoolWarm'),levels=250):
#     X=dfpivot.columns.values
#     Y=dfpivot.index.values
#     Z=dfpivot.values
#     plot_contourf_Impl(X,Y,Z,label,cmap=cmap,levels=levels)

def plot_contourf(df,fieldName,label,cmap=pplot.Colormap('CoolWarm'),levels=250,figwidth=20):
    dfpivot=df.pivot("y", "x", fieldName)
    X=dfpivot.columns.values
    Y=dfpivot.index.values
    Z=dfpivot.values
    plot_contourf_Impl(X,Y,Z,label,cmap=cmap,levels=levels,figwidth=figwidth)


def getTimesFromDataFolder(dataFolder):
    times=os.listdir(dataFolder)
    timeNames=[]
    for t in times:
        if(str.endswith(t,".csv")):
            timeNames.append(t.replace(".csv",""))
    timeNames=np.sort(timeNames)
    return timeNames

def plotTemperatureAndCokeEvolution(dataFolder,timeNames=[],workerNum=10):
    if(len(timeNames)==0):
        timeNames=getTimesFromDataFolder(dataFolder)
    futures=[]
    maxTemperatures=[]
    cokeFractions=[]
    with concurrent.futures.ThreadPoolExecutor(max_workers=workerNum) as executor:
        for timeName in timeNames:
            future=executor.submit(readData,dataFolder,timeName)
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


