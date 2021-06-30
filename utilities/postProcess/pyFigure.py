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
import matplotlib.gridspec as gridspec
from matplotlib.axes._axes import _log as matplotlib_axes_logger
matplotlib_axes_logger.setLevel('ERROR')
from colored import fg, attr
import proplot as pplot # there are some nice colormaps in the proplot package
import concurrent.futures
import argparse
import json
import math
import traceback
import tracemalloc




# mpl.rcParams['font.family'] = 'Helvetica', #default font family
mpl.rcParams['mathtext.fontset'] = 'cm' #font for math
C_GREEN = fg('green')
C_RED = fg('red')
C_BLUE = fg('blue')
C_DEFAULT = attr('reset')

c1 = pplot.scale_luminance('cerulean', 0.5)
c2 = pplot.scale_luminance('red', 0.5)

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
    # if(time_second>(86400*365)):
    #     Scale = 3.17e-08
    #     name_time='years'
    # elif(time_second>86400):
    #     Scale = 1.1574074074074073e-05
    #     name_time='days'
    # elif(time_second>3600):
    #     Scale = 0.0002777777777777778
    #     name_time='hours'
    # elif(time_second>60):
    #     Scale = 0.016666666666666666
    #     name_time='minutes'
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

def clippedcolorbar(CS, **kwargs):
    from matplotlib.cm import ScalarMappable
    from numpy import arange, floor, ceil
    fig = CS.ax.get_figure()
    vmin = CS.get_clim()[0]
    vmax = CS.get_clim()[1]
    m = ScalarMappable(cmap=CS.get_cmap())
    m.set_array(CS.get_array())
    m.set_clim(CS.get_clim())
    step = CS.levels[1] - CS.levels[0]
    cliplower = CS.zmin<vmin
    clipupper = CS.zmax>vmax
    noextend = 'extend' in kwargs.keys() and kwargs['extend']=='neither'
    # set the colorbar boundaries
    boundaries = arange((floor(vmin/step)-1+1*(cliplower and noextend))*step, (ceil(vmax/step)+1-1*(clipupper and noextend))*step, step)
    kwargs['boundaries'] = boundaries
    # if the z-values are outside the colorbar range, add extend marker(s)
    # This behavior can be disabled by providing extend='neither' to the function call
    if not('extend' in kwargs.keys()) or kwargs['extend'] in ['min','max']:
        extend_min = cliplower or ( 'extend' in kwargs.keys() and kwargs['extend']=='min' )
        extend_max = clipupper or ( 'extend' in kwargs.keys() and kwargs['extend']=='max' )
        if extend_min and extend_max:
            kwargs['extend'] = 'both'
        elif extend_min:
            kwargs['extend'] = 'min'
        elif extend_max:
            kwargs['extend'] = 'max'
    return fig.colorbar(m, **kwargs)

def plot_contourf_Impl(X,Y,Z,title,label,cmap=pplot.Colormap('CoolWarm'),levels=250,figwidth=20,vmin=0,vmax=0):
    figsize=figsize_cm(figwidth,X,Y)
    fig, ax = plt.subplots(figsize=figsize)
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
    # cbarticks = np.arange(vmin,vmax,(vmax-vmin)/10)
    # fig.colorbar(CS,cax=ax_cb,ticks=cbarticks,label=label)


    clippedcolorbar(CS,cax=ax_cb,label=label, extend='neither')

    fig.tight_layout()
    return fig,ax

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
    fig,ax=plot_contourf_Impl(X,Y,Z,title,label,cmap=cmap,levels=levels,figwidth=figwidth,vmin=vmin,vmax=vmax)
    return fig,ax

def plot_contourf_save(df,fieldName,title,label,folder_path,cmap=pplot.Colormap('CoolWarm'),levels=250,figwidth=20,vmin=0,vmax=0,dpi=600):
    dfpivot=df.pivot("y", "x", fieldName)
    X=dfpivot.columns.values
    Y=dfpivot.index.values
    Z=dfpivot.values
    fig,_=plot_contourf_Impl(X,Y,Z,title,label,cmap=cmap,levels=levels,figwidth=figwidth,vmin=vmin,vmax=vmax)
    if not os.path.exists(folder_path):
        os.mkdir(folder_path)
    plt.savefig(f"{folder_path}/{title}.jpg".replace(" ","-"),dpi=dpi)
    plt.close(fig)

    # snapshot = tracemalloc.take_snapshot()
    # top_stats = snapshot.statistics('lineno')  # lineno,逐行统计；filename，统计整个文件内存
    # print(top_stats[0])
    # for stat in top_stats[:1]:
       

def plot_multiple_contourf_save(df,fields,time_instant,save_folder,xranges={},dpi=600):
    if "eps" in fields:
        eps_title=f"porosity contour at {time_str(time_instant)}"
        if "eps" in xranges:
            vmin=xranges["eps"]["vmin"]
            vmax=xranges["eps"]["vmax"]
        else:
            vmin=0
            vmax=0
        plot_contourf_save(df,"eps",eps_title,label='porosity',folder_path=save_folder,vmin=vmin,vmax=vmax,dpi=dpi)

    if "UNorm" in fields:
        unorm_title=f"velocity magnitude contour at {time_str(time_instant)}"
        if "UNorm" in xranges:
            vmin=xranges["UNorm"]["vmin"]
            vmax=xranges["UNorm"]["vmax"]
        else:
            vmin=0
            vmax=0
        plot_contourf_save(df,"UNorm",unorm_title,label='velocity magnitude (m/s)',folder_path=save_folder,vmin=vmin,vmax=vmax,dpi=dpi)

    if "O2Conc" in fields:
        O2Conc_title=f"O$_{2}$ concentration contour at {time_str(time_instant)}"
        if "O2Conc" in xranges:
            vmin=xranges["O2Conc"]["vmin"]
            vmax=xranges["O2Conc"]["vmax"]
        else:
            vmin=0
            vmax=0
        plot_contourf_save(df,"O2Conc",O2Conc_title,label='O$_2$ mole concentration (mol/m$^3$)',folder_path=save_folder,vmin=vmin,vmax=vmax,dpi=dpi)

    if "T" in fields:
        T_title=f"Temperature contour at {time_str(time_instant)}"
        if "T" in xranges:
            vmin=xranges["T"]["vmin"]
            vmax=xranges["T"]["vmax"]
        else:
            vmin=0
            vmax=0
        # print(f"T vmax:{vmax}") 
        plot_contourf_save(df,"T",T_title,label='Temperature (K)',folder_path=save_folder,vmin=vmin,vmax=vmax,dpi=dpi)

    if "Qdot" in fields:
        Qdot_title=f"Reaction heat rate contour at {time_str(time_instant)}"
        if "Qdot" in xranges:
            vmin=xranges["Qdot"]["vmin"]
            vmax=xranges["Qdot"]["vmax"]
        else:
            vmin=0
            vmax=0
        # print(f"Qdot vmax:{vmax}") 
        plot_contourf_save(df,"Qdot",Qdot_title,label='Reaction heat rate (J/(m$^3\cdot$s))',folder_path=save_folder,vmin=vmin,vmax=vmax,dpi=dpi)

def read_plot_multiple_field_contourf_save(data_folder,fields,time_instant,save_folder,xranges={},dpi=600):
    df=read_data_and_process(data_folder,time_instant)
    plot_multiple_contourf_save(df,fields,float(time_instant),save_folder,xranges,dpi=dpi)



def read_plot_multiple_field_contourf_save_for_multiple_times(data_folder,fields,time_instants,save_folder,xranges={}):
    if len(time_instants)==0:
        return
    for time in time_instants:
        read_plot_multiple_field_contourf_save(data_folder,fields,time,save_folder,xranges)

def plot_multiple_contourf_save_all_time(fields,data_folder,save_folder,xranges={},time_names=[],worker_num=10):
    if(len(time_names)==0):
        time_names=get_times_from_data_folder(data_folder)

    futures=[]
    chunk_size=math.ceil(len(time_names)/worker_num)
    time_name_chunks=[time_names[x:x+chunk_size] for x in range(0, len(time_names), chunk_size)]
    print(f"time_name_chunks size: {len(time_name_chunks)}")
    print(f"time_name_chunks: {time_name_chunks}")
    with concurrent.futures.ProcessPoolExecutor(max_workers=worker_num) as executor:
        for chunk in time_name_chunks:
            future=executor.submit(read_plot_multiple_field_contourf_save_for_multiple_times,data_folder,fields,\
                chunk,save_folder,xranges)
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
        ax.set_ylabel("Maximum combustion Temperature (K)",color="b")
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

def read_min_max_field(file_path,sampling_rate,field):
    data=read_field_min_max_file(file_path)
    
    df=data[data["field"].str.contains(field)]
    df.reset_index(inplace=True,drop=True)

    df_sampling=df[df.index%sampling_rate==0]
    df_sampling.reset_index(inplace=True,drop=True)
    return df_sampling

def plot_min_max_field(file_path,sampling_rate,field,label,yscale="linear"):
    data_sampling=read_min_max_field(file_path,sampling_rate,f'^{field}')
    fig, ax = plt.subplots()
    ax.plot(data_sampling["Time"],data_sampling["max"],label=f"max {field}")
    ax.plot(data_sampling["Time"],data_sampling["min"],label=f"min {field}")
    ax.set_xlabel(f"Time (s)")
    ax.set_ylabel(label)
    ax.set_yscale(yscale)
    ax.set_title(f"Temporal Evolution of Max/Min {field}")
    ax.legend()
    ax.autoscale(tight=True)
    fig.tight_layout()
    return fig,ax,data_sampling
    # plt.close(fig)

def plot_transverse_averages(transverse_data_folder,time):
    df_transverse=pd.read_csv(os.path.join(transverse_data_folder,f"{time}.csv"))
    fig,ax=plt.subplots()
    c1 = pplot.scale_luminance('cerulean', 0.5)
    c2 = pplot.scale_luminance('red', 0.5)
    ax.plot(df_transverse["x"],df_transverse["O2Conc"],color=c1)
    ax.set_xlabel("X (m)")
    ax.set_ylabel("Tranversely Averaged O$_2$ mole concentration (mol/m$^3$)",color=c1)
    ax.tick_params(axis='y', colors=c1)

    ax2 = ax.twinx()
    ax2.plot(df_transverse["x"],df_transverse["T"],color=c2)
    ax2.set_ylabel("Tranversely Averaged Temperature (K)",color=c2)
    ax2.tick_params(axis='y', colors=c2)

    fig.tight_layout()
    return ax,ax2,fig

def plot_transverse_averages_of_multiple_times(transverse_data_folder,times):
    formatter=tickformatter()
    lines=["-",":","--","-.",(0,(0.01,2))]
    colors=["k","b","g","r"]

    transverse_data={}
    for time in times:
        df=pd.read_csv(f"{transverse_data_folder}/{time}.csv")
        transverse_data[time]=df

    fig = plt.figure(figsize=(9,10))
    outer = gridspec.GridSpec(2, 2, wspace=0.2, hspace=0.2)

    inner00=gridspec.GridSpecFromSubplotSpec(1,1,subplot_spec=outer[0])
    ax = plt.Subplot(fig, inner00[0])
    for i,time in enumerate(transverse_data.keys()):
        df=transverse_data[time]
        ax.plot(df["x"],df["O2Conc"],label=fr"$\mathit{{t}}\ $ = {time} s",linestyle=lines[i],color=colors[i])
        
    ax.set_xlabel("X (m)")
    ax.xaxis.set_major_formatter(formatter) 
    ax.set_ylabel("O$_2$ mole concentration (mol/m$^3$)")
    # ax.legend(loc='best', shadow=True, fancybox=True)
    fig.add_subplot(ax)
    ax0=ax


    inner01=gridspec.GridSpecFromSubplotSpec(1,1,subplot_spec=outer[1], hspace=0)
    ax = plt.Subplot(fig, inner01[0])
    for i,time in enumerate(transverse_data.keys()):
        df=transverse_data[time]
        ax.plot(df["x"],df["T"],label=fr"$\mathit{{t}}\ $ = {time} s",linestyle=lines[i],color=colors[i])
    ax.set_xlabel("X (m)")
    ax.xaxis.set_major_formatter(formatter) 
    ax.set_ylabel("Temperature (K)")
    fig.add_subplot(ax)

    axes = np.empty(shape=(4, 1), dtype=object)
    inner10=gridspec.GridSpecFromSubplotSpec(4,1,subplot_spec=outer[2],wspace=0, hspace=0)
    global_coke_max=0
    for i,time in enumerate(transverse_data.keys()):
        df=transverse_data[time]
        coke_max=df["coke"].max()
        if global_coke_max<coke_max:
            global_coke_max=coke_max

    for i,time in enumerate(transverse_data.keys()):
        ax=plt.Subplot(fig, inner10[i])
        df=transverse_data[time]
        ax.plot(df["x"],df["coke"],label=fr"$\mathit{{t}}\ $ = {time} s",linestyle=lines[i],color=colors[i])
        # ax.legend(loc="upper right")
        ax.set_ylim([-0.001,global_coke_max*1.1])
        axes[i, 0] = ax
        if ax.is_last_row():
            ax.xaxis.set_major_formatter(formatter) 
            ax.set_xlabel("X (m)")
        else:
            plt.setp(ax.get_xticklabels(), visible=False)
        fig.add_subplot(ax)
    plt.setp(axes[2,0], ylabel='coke fraction')



    axes = np.empty(shape=(4, 1), dtype=object)
    inner11=gridspec.GridSpecFromSubplotSpec(4,1,subplot_spec=outer[3],wspace=0, hspace=0)

    for i,time in enumerate(transverse_data.keys()):
        ax=plt.Subplot(fig, inner11[i],sharex=axes[0, 0], sharey=axes[0, 0])
        df=transverse_data[time]
        ax.plot(df["x"],df["Qdot"],label=fr"$\mathit{{t}}\ $ = {time} s",linestyle=lines[i],color=colors[i])
        # ax.legend(loc="upper right")
        # ax.set_ylim([-1e4,maxQdotOfAll*1.2])
        fig.add_subplot(ax)
        axes[i, 0] = ax
        if ax.is_last_row():
            ax.xaxis.set_major_formatter(formatter) 
            ax.set_xlabel("X (m)")
        else:
            plt.setp(ax.get_xticklabels(), visible=False)
    plt.setp(axes[2,0], ylabel='Reaction Heat Rate (J/(m$^3\cdot$s))')


    handles, labels = ax0.get_legend_handles_labels()
    fig.legend(handles, labels, loc=(0.2,0.92), shadow=True, fancybox=True,ncol=4,fontsize=12)
    return fig

def Plot_MaxTemperature_OutletO2ConcHistory(df_combined):
    fig, ax = pplot.subplots( aspect=(4, 3), axwidth=4)

    lns1=ax.plot(df_combined["Time"],df_combined["max"],color=c1,label="Max Point Temperature",linestyle="-")
    lns2=ax.plot(df_combined["Time"],df_combined["Transverse_Tmax"],color=c1,label="Max Transversely Averaged Temperature",linestyle="--")
    ax.format(xlabel="Time (s)",ylabel="Temperature (K)",
                ycolor=c1)

    ax2 = ax.twinx()
    lns3= ax2.plot(df_combined["Time"],df_combined["O2ConcAtOutlet"],color=c2,linestyle="-.",label="O$_2$ Molar Concentration at outlet")
    max_O2=df_combined["O2ConcAtOutlet"].max()
    ax2.format(ylabel="O$_2$ mole concentration At Outlet (mol/m$^3$)",
                ycolor=c2)

    lns = lns1+lns2+lns3
    labs = [l.get_label() for l in lns]
    ax.legend(lns, labs, loc="upper right", ncol=1, fancybox=True)

    # fig.tight_layout()
    return ax,ax2,fig

def plot_reaction_rate_burning_rate(df_rate):
    fig, ax = pplot.subplots( aspect=(4, 3), axwidth=4)
    c1 = pplot.scale_luminance('cerulean', 0.5)
    c2 = pplot.scale_luminance('red', 0.5)
    df_rate.sort_values(by="time",inplace=True)
    
    lns1=ax.plot(df_rate["time"],df_rate["vol_averaged_reaction_rate"],color=c1,
            label="Reaction rate")

    max_rate=df_rate["vol_averaged_reaction_rate"].max()
    ax.format(xlabel="Time (s)",ylabel="Volume-Averaged coke reaction rate (kg/m$^3$/s)",
              ycolor=c1,ylim=(-1,max_rate*1.1))

    ax2 = ax.twinx()
    lns2=ax2.plot(df_rate["time"],df_rate["burning_fraction"]*100,color=c2,
            linestyle="--",label="Conversion")
    ax2.format(xlabel="Time (s)",ylabel="Conversion (%)",ycolor=c2,
               ylim=(-1,100))
    
    lns = lns1+lns2
    labs = [l.get_label() for l in lns]
    ax.legend(lns, labs, loc="upper right", fancybox=True)

    return ax,ax2,fig

def plot_O2_flux_reaction_rate(df_O2_flux_at_inlet,df_rate,pixelResolution,DO2,sampling_rate=5,ylim=(1e-9, 1e-4)):
    MCoke=12
    MO2=32
    df_O2_flux_at_inlet["diffusive_flux"]=np.array(df_O2_flux_at_inlet["O2_diffusive_Flux_By_DO2"])*DO2
    df_O2_flux_at_inlet["advective_flux"]=np.array(df_O2_flux_at_inlet["O2_adv_flux_by_deltaX"])*pixelResolution
    df_O2_flux_at_inlet["total_flux"]=df_O2_flux_at_inlet["diffusive_flux"]+df_O2_flux_at_inlet["advective_flux"]

    df_O2_flux_at_inlet["ratio of diffusion flux"]=df_O2_flux_at_inlet["diffusive_flux"]/df_O2_flux_at_inlet["total_flux"]
    df_O2_flux_at_inlet["ratio of advection flux"]=df_O2_flux_at_inlet["advective_flux"]/df_O2_flux_at_inlet["total_flux"]


    c1 = pplot.scale_luminance('cerulean', 0.5)
    c2 = pplot.scale_luminance('red', 0.5)

    fig, axs = pplot.subplots(ncols=1, nrows=2, aspect=(4, 3), \
        axwidth=4,hspace=(0),sharey=0,sharex=3)
    
    ax=axs[0]
    
    ax.plot(df_O2_flux_at_inlet["time"],df_O2_flux_at_inlet["diffusive_flux"],color=c1,label="Diffusive Flux",linestyle="-.")
    ax.plot(df_O2_flux_at_inlet["time"],df_O2_flux_at_inlet["advective_flux"],color=c1,label="Advective Flux",linestyle="--")
    ax.plot(df_O2_flux_at_inlet["time"],df_O2_flux_at_inlet["total_flux"],color=c1,label="Total Flux",linestyle="-")
    ax.format(xlabel="Time (s)",ylim=ylim, yformatter='sci',yscale='log',
              ylabel="O$_2$ flux at inlet (kg/s)",ycolor=c1)
    ax.legend(loc="best", ncols=1, fancybox=True )
    ax2 = ax.twinx()
    df_sampling=df_rate[df_rate.index%sampling_rate==0]

    ax2.scatter(df_sampling["time"],df_sampling["total_reaction_rate"]/MCoke*MO2*pixelResolution*pixelResolution,label="Total O$_2$ Reaction Rate",color=c2)
    ax2.format(xlabel="Time (s)",ylim=ylim, yformatter='sci', yscale='log', ylabel='O$_2$ reaction rate within domain (kg/s)',ycolor=c2)
 
    df_sampling=df_O2_flux_at_inlet[df_O2_flux_at_inlet.index%sampling_rate==0]
    df_sampling.reset_index(drop=True,inplace=True)
    df_ratios=df_sampling[["ratio of diffusion flux","ratio of advection flux"]]
    ax3=axs[1]
    ax3.format(ylim=(0,1.1),ylabel="Ratio of O$_2$ advection/diffusion flux to total")
    num=df_sampling.shape[0]
    maxTime=df_sampling["time"].iloc[-1]
    ax3.bar(df_sampling["time"],df_ratios,stacked=True, cycle='Blues',
            legend='ur', edgecolor='blue9', width=maxTime/(num*1.2))

    return ax,ax2,fig



if __name__ == "__main__":

    # tracemalloc.start() # 开始跟踪内存分配
    try:
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
    except Exception as e:
        errmsg=f"Unhandled exception happened: {e} with stack trace {traceback.format_exc()}\n"
        print(errmsg)


