from logging import exception
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
mpl.rcParams['mathtext.fontset'] = 'cm'  # font for math
C_GREEN = fg('green')
C_RED = fg('red')
C_BLUE = fg('blue')
C_DEFAULT = attr('reset')

c1 = pplot.scale_luminance('cerulean', 0.5)
c2 = pplot.scale_luminance('red', 0.5)


def use_paper_style(mpl, fontsize=9):
    mpl.rcParams['axes.titlesize'] = fontsize
    mpl.rcParams['axes.labelsize'] = fontsize
    mpl.rcParams['xtick.labelsize'] = fontsize
    mpl.rcParams['ytick.labelsize'] = fontsize
    mpl.rcParams['legend.fontsize'] = fontsize
    mpl.rcParams['axes.titlesize'] = fontsize
    mpl.rcParams['axes.titlesize'] = fontsize


def read_postProcess_csv_data(postProcessDir, timeName):
    fileName = f"{timeName}.csv"
    filePath = os.path.join(postProcessDir, fileName)
    df = pd.read_csv(filePath)
    return df


def tickformatter():
    formatter = ticker.ScalarFormatter(useMathText=True)
    formatter.set_scientific(True)
    formatter.set_powerlimits((-1, 1))
    return formatter


def time_str(time_second):
    Scale, name_time = 1, 'seconds'
    str_time = str('%.3f %s' % (time_second*Scale, name_time))
    return str_time


def read_data_and_process(folder, timeName):
    df = read_postProcess_csv_data(folder, timeName)
    df["UNorm"] = np.sqrt(df['U_0']**2 + df['U_1']**2 + df['U_2']**2)
    MO2 = 0.032  # g/mol
    df["O2Conc"] = df['rho']*df['eps']*df['O2']/MO2
    MCO2 = 0.044  # g/mol
    df["CO2Conc"] = df['rho']*df['eps']*df['CO2']/MCO2
    return df


def figsize_cm(w_cm, x, y, w_offset_cm=0.0):
    # Get figure size in unit cm.
    cm2inch = 1/2.54
    w = w_cm*cm2inch

    length_x = x.max()-x.min()
    length_y = y.max()-y.min()
    xyRatio = length_y/length_x
    h = w*xyRatio
    w = w+w_offset_cm*cm2inch
    return (w, h)


def plot_contourf_Impl(X, Y, Z, title, label, figwidth, cmap=pplot.Colormap('CoolWarm'),
                       nlevels=250, vmin=0, vmax=0):
    figsize = figsize_cm(figwidth, X, Y)
    fig, ax = plt.subplots(figsize=figsize)
    # ax.axis('scaled')
    ax.set_xlim(X.min(), X.max())
    ax.set_ylim(Y.min(), Y.max())
    ax.set_xlabel('X (m)')
    ax.set_ylabel('Y (m)')
    ax.set_title(title)
    formatter = tickformatter()
    ax.xaxis.set_major_formatter(formatter)
    ax.yaxis.set_major_formatter(formatter)
   
    if vmin == 0:
        vmin = np.min(Z)
    if vmax == 0:
        vmax = np.max(Z)
    
    if vmin>np.min(Z) and vmax<np.max(Z):
        extend="both"
    elif vmin>np.min(Z):
        extend="min"
    elif vmax<np.max(Z):
        extend="max"
    else:
        extend="neither"

    levels = np.linspace(vmin, vmax, nlevels)
    CS = ax.tricontourf(X, Y, Z, cmap=cmap, levels=levels, extend=extend)
    ax_cb = ax.inset_axes([1.04, 0, 0.02, 1])
    fig.colorbar(CS, cax=ax_cb, label=label)
    fig.tight_layout()
    return fig, ax


def plot_contourf(df, fieldName, title, label, cmap, nlevels, figwidth, vmin, vmax):
    x=df["x"]
    y=df["y"]
    z = df[fieldName]
    fig, ax = plot_contourf_Impl(x, y, z, title, label, cmap=cmap, nlevels=nlevels,
                                 figwidth=figwidth, vmin=vmin, vmax=vmax)
    return fig, ax


def plot_contourf_save(df, fieldName, title, label, folder_path, cmap, nlevels, figwidth, vmin, vmax, dpi):
    try:
        fig, _ = plot_contourf(df, fieldName, title, label, cmap, nlevels, figwidth, vmin, vmax)
        if not os.path.exists(folder_path):
            os.mkdir(folder_path)
        plt.savefig(f"{folder_path}/{title}.jpg".replace(" ", "-"), dpi=dpi)
        plt.close(fig)
    except Exception as e:
        print(f"Error: failed in plotting {title} with error message: {e}")
        raise

    # snapshot = tracemalloc.take_snapshot()
    # top_stats = snapshot.statistics('lineno')  # lineno,逐行统计；filename，统计整个文件内存
    # print(top_stats[0])
    # for stat in top_stats[:1]:


def plot_multiple_contourf_save(df, fields, time_instant, save_folder, cmap_name, xranges={}, figwidth=20, dpi=600):
    cmap = pplot.Colormap(cmap_name)

    if "eps" in fields:
        eps_title = f"porosity contour at {time_str(time_instant)}"
        vmin, vmax, nlevels = read_colorbar_configuration("eps",xranges)
        if vmax==0:
            vmax=1.001
        plot_contourf_save(df, "eps", eps_title, label='porosity', cmap=cmap, folder_path=save_folder, nlevels=nlevels,
                           vmin=vmin, vmax=vmax, figwidth=figwidth, dpi=dpi)

    if "UNorm" in fields:
        unorm_title = f"velocity magnitude contour at {time_str(time_instant)}"
        vmin, vmax, nlevels = read_colorbar_configuration("UNorm",xranges)
        plot_contourf_save(df, "UNorm", unorm_title, label='velocity magnitude (m/s)', cmap=cmap,
                           folder_path=save_folder, nlevels=nlevels, vmin=vmin, vmax=vmax,
                           figwidth=figwidth, dpi=dpi)

    if "O2Conc" in fields:
        O2Conc_title = f"O$_{2}$ concentration contour at {time_str(time_instant)}"
        vmin, vmax, nlevels = read_colorbar_configuration("O2Conc",xranges)
        plot_contourf_save(
            df, "O2Conc", O2Conc_title, label='O$_2$ mole concentration (mol/m$^3$)', cmap=cmap,
            folder_path=save_folder, nlevels=nlevels, vmin=vmin, vmax=vmax,figwidth=figwidth, dpi=dpi)

    if "T" in fields:
        T_title = f"Temperature contour at {time_str(time_instant)}"
        vmin, vmax, nlevels = read_colorbar_configuration("T",xranges)
        plot_contourf_save(df, "T", T_title, label='Temperature (K)', cmap=cmap, folder_path=save_folder,
                           nlevels=nlevels, vmin=vmin, vmax=vmax, figwidth=figwidth, dpi=dpi)

    if "p" in fields:
        p_title = f"Pressure contour at {time_str(time_instant)}"
        vmin, vmax, nlevels = read_colorbar_configuration("p",xranges)
        plot_contourf_save(df, "p", p_title, label='Pressure (Pa)', cmap=cmap, folder_path=save_folder, nlevels=nlevels,
                           vmin=vmin, vmax=vmax, figwidth=figwidth, dpi=dpi)

    if "Qdot" in fields:
        Qdot_title = f"Reaction heat rate contour at {time_str(time_instant)}"
        vmin, vmax, nlevels = read_colorbar_configuration("Qdot",xranges)
        plot_contourf_save(
            df, "Qdot", Qdot_title, label='Reaction heat rate (J/(m$^3\cdot$s))', cmap=cmap, folder_path=save_folder,
            nlevels=nlevels, vmin=vmin, vmax=vmax, figwidth=figwidth, dpi=dpi)


def read_colorbar_configuration(field,xranges):
    default_nlevels = 100
    if field in xranges:
        vmin = xranges[field]["vmin"]
        vmax = xranges[field]["vmax"]
        nlevels = xranges[field]["nlevels"]
    else:
        vmin = 0
        vmax = 0
        nlevels = default_nlevels
    return vmin, vmax, nlevels


def read_plot_multiple_field_contourf_save(
        data_folder, fields, time_instant, save_folder, cmap_name, xranges, figwidth, dpi):
    df = read_data_and_process(data_folder, time_instant)
    plot_multiple_contourf_save(
        df=df, fields=fields, time_instant=float(time_instant),
        save_folder=save_folder, cmap_name=cmap_name, xranges=xranges, figwidth=figwidth, dpi=dpi)


def read_plot_multiple_field_contourf_save_for_multiple_times(
        data_folder, fields, save_folder, xranges, cmap_name, figwidth, dpi, time_instants=[]):
    if len(time_instants) == 0:
        return
    for time_instant in time_instants:
        read_plot_multiple_field_contourf_save(
            data_folder=data_folder, fields=fields, time_instant=time_instant, save_folder=save_folder,
            cmap_name=cmap_name, xranges=xranges, figwidth=figwidth, dpi=dpi)


def plot_multiple_contourf_save_all_time(
        fields, data_folder, save_folder, time_names, xranges, cmap_name, figwidth, dpi, worker_num=10):
    futures = []
    chunk_size = math.ceil(len(time_names)/worker_num)
    time_name_chunks = [time_names[x:x+chunk_size] for x in range(0, len(time_names), chunk_size)]
    print(f"time_name_chunks size: {len(time_name_chunks)}")
    print(f"time_name_chunks: {time_name_chunks}")
    with concurrent.futures.ProcessPoolExecutor(max_workers=worker_num) as executor:
        for chunk in time_name_chunks:
            future = executor.submit(read_plot_multiple_field_contourf_save_for_multiple_times, data_folder, fields,
                                     save_folder, xranges, cmap_name, figwidth, dpi, chunk)
            futures.append(future)

        for _, future in enumerate(futures):
            future.result()


def get_times_from_data_folder(dataFolder):
    times = os.listdir(dataFolder)
    timeNames = []
    for t in times:
        if(str.endswith(t, ".csv")):
            timeNames.append(t.replace(".csv", ""))
    timeNames = np.sort(timeNames)
    return timeNames


if __name__ == "__main__":

    # tracemalloc.start() # 开始跟踪内存分配
    try:
        parser = argparse.ArgumentParser(description="pyResconstruct")
        parser.add_argument("-f", dest='field_names', required=True, help='specify the field names')
        parser.add_argument("-d", dest='data_folder', required=True, help='specify the data folder')
        parser.add_argument("-s", dest='save_folder', required=True, help="specify the save folder")
        parser.add_argument("-x", dest='xranges', default="{}", help='specify the color bar ranges')
        parser.add_argument("-c", dest='cmap', default="CoolWarm", help='specify the color map')
        parser.add_argument("-p", dest='dpi', default=300, help='specify dpi')
        parser.add_argument("-w", dest='fig_width', default=20, help='specify figure width')
        parser.add_argument("-t", dest='time_names', default="all", help='specify the time names')
        parser.add_argument("-n", dest='worker_num', default=8, type=int, help='specify the worker num')

        args = parser.parse_args()
        data_folder = args.data_folder
        field_names_texts = args.field_names
        field_names = json.loads(field_names_texts)
        xranges_texts = args.xranges
        xranges = json.loads(xranges_texts)
        cmap = args.cmap
        dpi = args.dpi
        fig_width = args.fig_width
        worker_num = args.worker_num
        save_folder = args.save_folder
        time_names = args.time_names

        if time_names == "all":
            all_time_names = get_times_from_data_folder(data_folder)
            plot_multiple_contourf_save_all_time(
                fields=field_names, data_folder=data_folder, save_folder=save_folder, time_names=all_time_names,
                xranges=xranges, cmap_name=cmap, figwidth=fig_width, dpi=dpi, worker_num=worker_num)
        else:
            print(f"specified time names: {time_names}")
            time_names = json.loads(time_names)
            plot_multiple_contourf_save_all_time(
                fields=field_names, data_folder=data_folder, save_folder=save_folder, time_names=time_names,
                xranges=xranges, cmap_name=cmap, figwidth=fig_width, dpi=dpi, worker_num=worker_num)

        print("succeed to plot images")
    except Exception as e:
        errmsg = f"Unhandled exception happened: {e} with stack trace {traceback.format_exc()}\n"
        print(errmsg)
