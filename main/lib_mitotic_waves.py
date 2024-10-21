import os
import glob

import numpy as np
from numpy.lib.scimath import sqrt as csqrt
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
from matplotlib import rc,rcParams,colors
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
from matplotlib.offsetbox import OffsetImage, AnnotationBbox


from scipy.integrate import odeint
from scipy.io import FortranFile
from scipy.optimize import fsolve
from scipy.signal import find_peaks
from scipy.ndimage import gaussian_filter
from scipy.signal import detrend, savgol_filter
from scipy.signal import coherence
from scipy.interpolate import interp1d
from scipy.interpolate import interp2d
from scipy.interpolate import UnivariateSpline
from scipy.special import ndtr
from scipy.stats import norm, gaussian_kde
from scipy.optimize import curve_fit
import statsmodels.api as sm


import pandas as pd
import openpyxl
import rasterio
from rasterio.plot import show
import copy
import pickle
import f90nml
import time
import warnings
from collections import OrderedDict
import subprocess
from functools import reduce



    
plt.rcdefaults()

font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 18}

plt.rcParams.update(
{"figure.figsize": (1.6*4, 4),
"text.usetex": False,
"axes.titlesize" : 18,
"axes.labelsize" : 18,
"legend.fontsize": 18,
"legend.fontsize": 18,
"xtick.labelsize" : 18,
"ytick.labelsize" : 18,
'xtick.major.width' : 2,
'ytick.major.width' : 2,
"savefig.dpi": 300})


class chang2013model:  
    def __init__(self, params_dict = {}):
        self.tf_timeseries = 100
        self.nt_timeseries = 10001
        self.tf_period = 1000
        self.nt_period = 10001

        #Import params from param_dict
        self.params_dict = params_dict
    
    def ode_eqs(self, y, t, ks, a_cdc25, b_cdc25, ec50_cdc25, n_cdc25, a_wee1, b_wee1, ec50_wee1, n_wee1, a_deg, b_deg, ec50_deg, n_deg, Da, Di,tau,alpha,beta):
        
        a,c = y
        p_cdc25 = np.array([ a_cdc25, b_cdc25, ec50_cdc25, n_cdc25 ])
        p_wee1 = np.array([ a_wee1, b_wee1, ec50_wee1, n_wee1 ])
        p_deg = np.array([ a_deg, b_deg, ec50_deg, n_deg ])

        dydt = [tau*(beta*(ks -H_a(p_deg,a)*a) +alpha*(H_a(p_cdc25,a)*(c-a) -H_i(p_wee1,a)*a) ),
                tau*beta*(ks -H_a(p_deg,a)*c)]

        return dydt

        
    def int_ode_eqs(self, y, t): 
        dydt = self.ode_eqs(y,t, **self.params_dict)
        return dydt
    
    def time_series(self, y0 = None,tf=None,nt = None):
        if tf is None:
            tf = self.tf_timeseries
        if nt is None:
            nt = self.nt_timeseries
        if y0 is None:
            y0 = [0 for i in range(self.dim)]
        t = np.linspace(0, tf, nt)
        sol = odeint(self.int_ode_eqs, y0, t)
        return tuple([t] + [sol[:,i] for i in range(sol.shape[1])])
    
    def plot_phaseplane(self, ax, trange = [], xrange = [0,80],yrange = [0,100], color='gray'):
        y0=[0.001,0.001]
        t, a, c = self.time_series(y0=y0)
        ax.plot(c,a,c = color)

        ax.set_xlabel('Total cyclin B - Cdk1 (nM)')
        ax.set_ylabel('Active cyclin B - Cdk1 (nM)')
        ax.set_ylim(yrange)
        ax.set_xlim(xrange)
    
def H_a(param,x):
    a,b,ec50,n = param
    hill = a + b*x**n/(ec50**n + x**n)
    return hill

def H_i(param,x):
    a,b,ec50,n = param
    hill = a + b*ec50**n/(ec50**n + x**n)
    return hill
    
def custom_colormap(palette_colors, index = False):
    if not index:
        index = np.linspace(0,1,len(palette_colors))
    indx_clr = list(zip(index,palette_colors))
    cm = colors.LinearSegmentedColormap.from_list('blue_red', indx_clr, 256)
    return cm

def hex_to_rgb(value):
    value = value.lstrip('#')
    lv = len(value)
    return tuple(int(value[i:i + lv // 3], 16) for i in range(0, lv, lv // 3))

def create_with_timecolormap_expt( x,y,z, cmap="plasma",linewidth=2,linestyle='solid'):
    points = np.array([x[0:-1:100], y[0:-1:100]]).T.reshape(-1,1,2)
    segments = np.concatenate([points[:-1],points[1:]], axis=1)
    norm = plt.Normalize(np.amin(z),np.amax(z))
    lc = LineCollection(segments, cmap=cmap, norm=norm,linewidth=linewidth,linestyle=linestyle)
    lc.set_array(z[0:-1:100])
    return lc

def custom_colormap_palette(palette_colors, index = False):
    if not index:
        index = np.linspace(0,1,len(palette_colors))
    indx_clr = list(zip(index,palette_colors))
    cm = colors.LinearSegmentedColormap.from_list('blue_red', indx_clr, 256)
    return cm

def plot_xyt(ax0, x,y,t, cols, linewidth = 3,label=None):
    xx, yy, tt = x.copy(), y.copy(), t.copy()
    cmap = custom_colormap_palette(cols)
    lc0 = create_with_timecolormap( x, y, t, cmap = cmap,linewidth=linewidth)
    ax0.add_collection(lc0)
    
def plot_xyt1(ax0,ax1,ax2, x,y,t, cols, linewidth = 3,label=None):
    xx, yy, tt = x.copy(), y.copy(), t.copy()
    norm = colors.Normalize(vmin=np.min(tt), vmax=np.max(tt))
    cmap = custom_colormap(cols)
    lc0 = create_with_timecolormap( x, y, t, cmap = cmap,linewidth=linewidth,norm=norm)
    lc1 = create_with_timecolormap( t, x, t, cmap = cmap,linewidth=linewidth,norm=norm)
    lc2 = create_with_timecolormap( t, y, t, cmap = cmap,linewidth=linewidth,norm=norm)
    ax0.add_collection(lc0)
    ax1.add_collection(lc1)
    ax2.add_collection(lc2)

    # ax1.plot(t, x, label=label,c=cols[-1],linewidth=lw)
    # ax2.plot(t, y, label=label,c=cols[-1],linewidth=lw)

def color_scale(palette_colors = ['#402D77','#009ED6', '#009642','#FBE136','#E5382C']):
    index = np.linspace(0,1,len(palette_colors))
    indx_clr = list(zip(index,palette_colors))
    colormap = colors.LinearSegmentedColormap.from_list('blue_red', indx_clr, 256)
    return colormap

def getImage(path, zoom=1):
    return OffsetImage(plt.imread(path), zoom=zoom)

# Simulations
def get_file_folder(file_path):
    ocur = findoccurrences(file_path,'/')
    if not ocur:
        file_name = file_path
        folder = ''    
    else:
        index = ocur[-1]
        file_name = file_path[index+1:]
        folder = file_path[:index+1]
    return file_name, folder

def findoccurrences(s, ch):
    return [i for i, letter in enumerate(s) if letter == ch]

def import_dim(folder, file_name = 'dim.h'): 
    f1 = open(folder + file_name, "r")
    ll = f1.readlines()
    dim = {}
    for el in ll:
        keyval = el.replace('integer, parameter :: ','').replace('=','').replace('!system size.\n','').replace(' ','')
        dim[keyval[0]] = int(keyval[1:])
    return  dim

def get_n_writen(folder):
    model_params, integration_params, init_params = import_param_sim(folder)    
    delta_t, t_write, t_final = integration_params['delta_t'],integration_params['t_write'],integration_params['t_final']
    n_writen = round(t_final/t_write) 
    return n_writen

def import_param_sim(folder, file_name='param.in'): 
    nml = f90nml.read(folder + file_name)
    return dict(nml['model_params'].todict()), dict(nml['integration_params'].todict()), dict(nml['init_params'].todict())

def import_evolution(folder,file_0 = 'a.dat', file_1 = 'i.dat'):
    dims = import_dim(folder)
    n, m = dims['n'], dims['m']
    n_writen = get_n_writen(folder)
    array0 = np.zeros((n_writen,m,n))
    array1 = np.zeros((n_writen,m,n))
    f = FortranFile(folder + file_0,'r')
    g = FortranFile(folder + file_1,'r')

    for k in range(n_writen):
        dataf = f.read_reals('f8')
        datag = g.read_reals('f8')
        array0[k] = np.array([[dataf[i+n*j] for i in range(n)] for j in range(m)])
        array1[k] = np.array([[datag[i+n*j] for i in range(n)] for j in range(m)])

    return array0, array1

def create_spacetime(folder):
    dim = import_dim(folder)
    m,n = dim['m'],dim['n']
    model_params, integration_params, init_params = import_param_sim(folder) 
    dx = integration_params['delta_x']
    dy = integration_params['delta_y']
    dt = integration_params['t_write']
    t_final = integration_params['t_final']
    n_t = int(round(t_final/dt))
    x = np.arange(n)*dx
    y = np.arange(m)*dy
    t = np.arange(n_t)*dt
    return x,y,t

def import_field(file_path): 
    file_name, folder = get_file_folder(file_path) 
    dims = import_dim(folder)
    n, m = dims['n'], dims['m']
    f = FortranFile(file_path,'r')
    data = f.read_reals('f8')
    field = np.array([[data[i+n*j] for i in range(n)] for j in range(m)])
    return field  

def writen_time_array_plus_one(folder):
    model_params, integration_params, init_params = import_param_sim(folder) 
    t_write = integration_params['t_write']
    tt  = np.arange(get_n_writen(folder)+1)*t_write
    return tt

def import_time_array(file_path,n=False): 
    file_name, folder = get_file_folder(file_path) 
    if n:
        array = np.zeros(n)
    else:
        array = create_time_array(folder)
    
    f = FortranFile(file_path,'r')
    for k in range(len(array)): 
        data = f.read_reals('f8')
        array[k] = data
    f.close() 
    return array

# Waves detection
def get_traj(t,x,field):
    it_cycles = get_cycles_array(field) 
    n = field.shape[1]
#     n_cycles = len(it_cycles[-1])
    n_cycles = max([len(el) for el in it_cycles])
    waves = []
    for j in range(n_cycles):
        wave = []
        for i, el in enumerate(it_cycles):
            if j<len(el):
                #print(j, len(el))
                wave.append([t[el[j]],x[i]])
        wave = np.array(wave).T
        waves.append(wave)
    return waves #returns t, x

def get_cycles(time_series,pro_per=0.2,dist=5):
    amplitudes = np.amax(time_series) -  np.amin(time_series)
    zeros = find_peaks(time_series, distance=dist,prominence=pro_per*(amplitudes))[0]
    if len(zeros)<2 or np.amax(amplitudes) <1e-4:
        zeros = []
    return zeros

def get_cycles_array(field):
    cycles_array = []
    for i in range(field.shape[1]):
        cycles_array.append(get_cycles(field[:,i]))
    return cycles_array

def spl_traj(traj0, dt, plot=False):
    traj = traj0[:,~np.isnan(traj0[0])]
    p_knots = np.where(np.diff(traj[0],1) > 0)[0]
    n_knots = np.array(np.where(np.diff(traj[0],1) < 0)[0])+1
    knots = np.unique(np.sort(np.concatenate((p_knots,n_knots,np.array([0,len(traj[0])-1])))))

    xx,yy = traj[1,knots],traj[0,knots]
    if len(knots) == 2: 
        xx,yy = traj[1],traj[0]
    if len(xx)<=2:
        return [np.nan], [np.nan]
    else:
        weights = np.zeros(len(xx))+np.sqrt(12)/dt
        spl = UnivariateSpline(xx,yy,k=2,w = weights,s=0)
        y_spl = spl(traj[1])
        if plot:
            plt.plot(traj[1],traj[0],'.',c='purple')
            plt.plot(traj[1],y_spl,c='red')
            plt.plot(traj[1,p_knots],traj[0,p_knots],'.',c='red')
            plt.plot(traj[1,n_knots],traj[0,n_knots],'.',c='green')
            plt.plot(traj[1,knots],traj[0,knots],'+',c='black')

        for i in range(len(knots)-1):
            j0,j1 = knots[i], knots[i+1]
            # print('aqui:',abs(y_spl[j0] - y_spl[j1]) > 0.9*dt)
            # print('aqui1:',y_spl[j0], y_spl[j1],y_spl[j0] - y_spl[j1] )

            if abs(y_spl[j0] - y_spl[j1]) > 0.9*dt: # Check if knots have different t (Linear interpolation)
                l0,l1 = min(y_spl[j0], y_spl[j1]),max(y_spl[j0], y_spl[j1])
                # print('js:',j0,j1)
                # print(l0,l1)
                # print((y_spl[j0:j1] < l0).any())
                # print((y_spl[j0:j1] > l1).any())
                if (y_spl[j0:j1] < l0).any() or (y_spl[j0:j1] > l1).any(): # Is the curve outside [t,t+-dt]

                    # Plot wrong trajectory
                    x_err, y_err = traj[1,j0:j1],y_spl[j0:j1]
                    if plot:
                        plt.plot(x_err,y_err, c='yellow')
                    # Interpolate new trajectory    
                    xxx = traj[1,j0:j1]
                    x0,x1 = traj[1,j0], traj[1,j1]
                    y0,y1 = y_spl[j0], y_spl[j1]
                    y_spl[j0:j1] = (y1-y0)/(x1-x0)*(xxx-x0)+y0 

            elif abs(y_spl[j0] - y_spl[j1]) < 0.9*dt: # Check if knots have the same t (Quadratic interpolation)

                if abs(traj[0,j0+1] - traj[0,j0]) < (0.9*dt): # Positive parabola
                    l0,l1 = y_spl[j0] - dt, y_spl[j0]
                    if (y_spl[j0:j1] < l0).any() or (y_spl[j0:j1] > l1).any(): # Is the curve outside [t-dt,t]
                        # Plot wrong trajectory
                        x_err, y_err = traj[1,j0:j1],y_spl[j0:j1]
                        if plot:
                            plt.plot(x_err,y_err, c='yellow')
                        # Interpolate new trajectory    
                        xxx = traj[1,j0:j1]
                        x0,x1 = traj[1,j0], traj[1,j1]
                        sign = +1 # Positive parabola
                        a = sign*2*(dt)/(x0 - x1)**2
                        y_spl[j0:j1] = a*(xxx-x0)*(xxx-x1)+ traj[0,j0]

                else: # Negative parabola
                    l0,l1 = y_spl[j0], y_spl[j0] + dt 
                    if (y_spl[j0:j1] < l0).any() or (y_spl[j0:j1] > l1).any(): # Is the curve outside [t-dt,t]
                        # Plot wrong trajectory
                        x_err, y_err = traj[1,j0:j1],y_spl[j0:j1]
                        if plot:
                            plt.plot(x_err,y_err, c='yellow')
                        # Interpolate new trajectory    
                        xxx = traj[1,j0:j1]
                        x0,x1 = traj[1,j0], traj[1,j1]
                        sign = -1# Negative parabola
                        a = sign*2*(dt)/(x0 - x1)**2
                        y_spl[j0:j1] = a*(xxx-x0)*(xxx-x1)+ traj[0,j0]

        if plot:
            plt.show()
        return y_spl,traj[1]

def spl_waves_traj(waves_traj,dt):
    spl_waves = []
    for el in waves_traj:
        tt,xx = spl_traj(el,dt)
        spl_waves.append(np.array([tt,xx]))
    return spl_waves


def extract_data_v1(dir_el, plot_crop=False,plot_traj=False,plot_cur=False,sigma=1,order=150):
    x,t,array,limits = dir_el['x'],dir_el['t'],dir_el['Activity'],dir_el['limits']
    x,t,array = crop(x,t,array,limits)
    if plot_crop:
        plt.pcolormesh(t,x,array)
    clean_array = recursive_fill_nans(array,sx=10,sy=10)
    array0 = detrend(clean_array,axis=1)
    array = gaussian_filter(array0, sigma=sigma)   

    waves_traj = get_traj(t,x,array.T)
    new_waves_traj = reconstruct_waves(waves_traj)
    if plot_traj:
        for el in new_waves_traj:
            plt.plot(el[0],el[1],'.',markersize=2)
    
    #waves = poly_waves_traj(waves_traj,order=order)
    waves = spl_waves_traj(new_waves_traj,t[1]-t[0])
    if plot_cur:
        for el in waves:
            plt.plot(el[0],el[1],c='purple')    
    if plot_traj or plot_crop or plot_cur:
        plt.show()
    data_t = get_period_speed(waves,x,t,array0) 
    return data_t 

def extract_data_plt(ax,dir_el, xticks,flip=False,plot_crop=False,plot_traj=False,plot_cur=False,sigma=1,order=150):
    x,t,array,limits = dir_el['x'],dir_el['t'],dir_el['Activity'],dir_el['limits']
    x,t,array = crop(x,t,array,limits)
    if plot_crop:
        if flip:
            ax.pcolormesh(t,np.flipud(x),array,cmap='jet')
        else:
            ax.pcolormesh(t,x,array,cmap='jet')
    # clean_array = recursive_fill_nans(array,sx=10,sy=10)
    clean_array = fill_nans(array,sx=15,sy=15)
    clean_array1 = fill_nans(clean_array,sx=15,sy=15)
    array0 = detrend(clean_array1,axis=1)
    array = gaussian_filter(array0, sigma=sigma)   

    waves_traj = get_traj(t,x,array.T)
    new_waves_traj = reconstruct_waves(waves_traj)
    if plot_traj:
        for el in new_waves_traj:
            ax.plot(el[0],el[1],'.',markersize=2,c='66C3FF')
    
    #waves = poly_waves_traj(waves_traj,order=order)
    waves = spl_waves_traj(new_waves_traj,t[1]-t[0])
    if plot_cur:
        for el in waves:
            if flip:
                ax.plot(el[0],np.flipud(el[1]),c='w')
            else:
                ax.plot(el[0],el[1],c='w') 
    if plot_traj or plot_crop or plot_cur:
        ax.set_xlabel('Time (min)')
        ax.set_ylabel('Space (mm)')
        # ax.set_yticks([0,2000,4000,6000,8000,10000])
        ax.set_yticklabels(['0.0','2.0','4.0','6.0','8.0','10.0'])
        ax.set_xticks(xticks)
        ax.set_xticklabels(['0','500','1000','1500'])
        # plt.show()
    data_t = get_period_speed(waves,x,t,array0) 
    return data_t

def crop(x,t,array, limits):
    return x[limits[0]:limits[1]],t[limits[2]:limits[3]],array[limits[0]:limits[1],limits[2]:limits[3]]

def recursive_fill_nans(indata,sx=5,sy=5):
    while np.isnan(indata).any():
        indata = fill_nans(indata,sx,sy)
    return indata

def reconstruct_waves(waves_traj,plot=False):
    new_waves_traj = waves_traj.copy()
    for i in range(1,len(waves_traj)):
        el = waves_traj[i].copy()
        prev_wave = waves_traj[i-1].copy()
        if len(el[1])>1:
            dx = np.min(np.diff(el[1]))
            periods = []
            for ix in range(len(el[0])):
                iix = np.where(abs(prev_wave[1] - el[1,ix]) < dx/2)[0] # Find same position than ix element
                if iix.size != 0:
                    periods.append(el[0,ix]- prev_wave[0,iix])
            period = np.average(periods)
            if plot:
                plt.plot(el[1],el[0],'+',c='blue')
            for j in range(1,len(el[0])):
                k = j-1
                prev_x = el[0,k]
                while np.isnan(prev_x):
                    if k == 0:
                        print('Error: Element 0 is NaN')
                        sys.exit(1)
                    else:    
                        k = k-1
                        prev_x = el[0,k]
                if abs(el[0,j] - prev_x) > period/3: # Is the data an outlier
                    repl = []
                    for iw,wave in enumerate(waves_traj): #Iterate for all waves and find a good replacement
                        if plot:
                            if abs(iw-i)<3:
                                plt.plot(wave[1],wave[0],'.',c='gray',markersize=1)
                        jj = np.where(abs(wave[1] - el[1,j]) < dx/2)[0] # Find same position than j element
                        if jj.size != 0:
                            repl.append(wave[0,jj])
                    repl = np.array(repl) 
                    pot_repl = repl[abs(repl - prev_x) < period/3]
                    if len(pot_repl) == 0:
                        el[0,j] = np.nan   
                    else:
                        pot_repl = pot_repl[np.argmin(abs(pot_repl - prev_x) < period/3)] #Find best replacement
                        if plot:
                            plt.plot(el[1,j],pot_repl,'.',c='red')
                        el[0,j] = pot_repl

            new_waves_traj[i] = el.copy()
            if i == len(waves_traj)-1:
                el = waves_traj[i].copy()
                if plot:
                    plt.plot(el[1],el[0],'+',c='blue')
                nn = int((el[1,-1]-el[1,0])/dx) + 1
                xx = np.linspace(el[1,0],el[1,-1],nn)
                new_el = np.zeros((2,nn)) + np.nan
                new_el[1] = xx.copy()
                ix = ((el[1]-el[1,0])/dx).astype(int)
                iix = np.arange(ix[-1])
                for kx,kkx in enumerate(ix):
                    new_el[0,kkx] = el[0,kx]+0
                missing_x = xx[np.setdiff1d(iix,ix)]
                # for line in missing_x:
                #     plt.axvline(line)
                for xj in missing_x: # Iterate over empty points
                    repl = []
                    for iw,wave in enumerate(waves_traj): #Iterate for all waves and find a good replacement
                        if plot:
                            if abs(iw-i)<3:
                                plt.plot(wave[1],wave[0],'.',c='gray',markersize=1)
                        jj = np.where(abs(wave[1] - xj) < dx/2)[0]
                        if jj.size != 0:
                            repl.append(wave[0,jj])
                    repl = np.array(repl) 
                    pot_repl = repl[abs(repl - prev_x) < period/3]
                    if len(pot_repl) != 0:
                        pot_repl = pot_repl[np.argmin(abs(pot_repl - prev_x) < period/3)] #Find best replacement
                        j = np.where(abs(new_el[1] - xj) < dx/2)[0]
                        if plot:
                            plt.plot(new_el[1,j],pot_repl,'.',c='red')
                        new_el[0,j] = pot_repl
                new_waves_traj[i] = new_el.copy()
        if plot:
            plt.plot(el[1],el[0],c='green')
            plt.show()
    return new_waves_traj

def get_period_speed(waves,x,t,array0):
    nn = 0
    for el in waves: 
        nn = nn + el.shape[1]
    data = np.empty((nn,12)) + np.nan # icycle, time, period, speed
    
    dadx,dadt = np.gradient(array0,x,t)
    dadxdt,dadt2 = np.gradient(dadt,x,t)
    dadx2,dadxdt = np.gradient(dadx,x,t)

    j=0
    for i, el in enumerate(waves):
        i_x = np.rint((el[1]-x[0])/(el[1,1]-el[1,0])).astype(int)
        i_t = np.rint((el[0]-t[0])/(t[1]-t[0])).astype(int)
        if i!=0:
            pel = waves[i-1]
            pi_x = np.rint((pel[1]-x[0])/(pel[1,1]-pel[1,0])).astype(int)
            pi_t = np.rint((pel[0]-t[0])/(t[1]-t[0])).astype(int)
            # intersection = np.in1d(pi_x,i_x)
            inter, pind, ind = np.intersect1d(pi_x,i_x,return_indices=True)
            # periods = el[0] - pel[0,intersection]
            periods = el[0,ind] - pel[0,pind]
        else:
            pi_x = np.rint((x-x[0])/(x[1]-x[0])).astype(int)
            pi_t = np.zeros(len(el[0])).astype(int)
            # intersection = np.in1d(pi_x,i_x) 
            inter, pind, ind = np.intersect1d(pi_x,i_x,return_indices=True)
            periods = np.zeros(len(el[0]))
            
        speeds = np.gradient(el[1],1)/np.gradient(el[0],1)
        slope = np.gradient(el[0],1)/np.gradient(el[1],1)
        iaccel = np.gradient(slope,1)/np.gradient(el[1],1)
        loc_std = local_std(speeds)
        
        f_t = el[0,ind]
        avg_t = np.average(f_t)
        f_speeds = abs(speeds[ind])
        f_loc_std = loc_std[ind]
        f_i_x = i_x[ind]
        f_pt = pi_t[pind]
        f_i_t = i_t[ind]
        f_slope = slope[ind]
        f_iaccel = iaccel[ind]
        for ti,pi,si,stdi,ji_x,jpi_t,ji_t,sli,ai in zip(f_t,periods,f_speeds,f_loc_std,f_i_x,f_pt,f_i_t,f_slope,f_iaccel):
            if len(dadt[ji_x,jpi_t:ji_t]) > 0:
                dadti = np.amax(dadt[ji_x,jpi_t:ji_t])
                dadt2i = np.amax(dadt2[ji_x,jpi_t:ji_t])
                dadxi = np.amax(abs(dadx[ji_x,jpi_t:ji_t]))
                dadx2i = np.amax(abs(dadx2[ji_x,jpi_t:ji_t]))
            else:
                dadti, dadt2i, dadxi, dadx2i = np.nan,np.nan,np.nan,np.nan
            data[j,:] = np.array([i,ti,pi,si,stdi,dadti,dadt2i,dadxi,dadx2i,ti-avg_t,sli,ai])
            j+=1         
    return data

def local_std(arr,sigma=1):
    local_avg = gaussian_filter(arr,sigma)
    local_dev = arr - local_avg
    local_std = np.sqrt(gaussian_filter(local_dev**2,sigma))
    return local_std


def plot_xyt_expt(ax0, x,y,t, cols, linewidth = 3,label=None):
    xx, yy, tt = x.copy(), y.copy(), t.copy()
    cmap = custom_colormap_palette(cols)
    lc0 = create_with_timecolormap_expt( x, y, t, cmap = cmap,linewidth=linewidth)
    ax0.add_collection(lc0)
    
    
def group_tubes_v1(dd,depth_dd):
    if depth_dd == 1:
        icol = 0
        for k,v in dd.items(): # Iterate over tubes
            data_ti = extract_data_v1(v)
            if icol==0:
                data_t = data_ti.copy()
                icol=1
            else:
                data_t = np.concatenate((data_t,data_ti))
    elif depth_dd == 2:
        icol = 0
        for kk,vv in dd.items(): # Iterate over positions
            for k,v in vv.items(): # Iterate over tubes
                data_ti = extract_data_v1(v)
                if icol==0:
                    data_t = data_ti.copy()
                    icol=1
                else:
                    data_t = np.concatenate((data_t,data_ti))
    elif depth_dd == 3:
        icol = 0
        for kkk,vvv in dd.items(): # Iterate over experiments
            for kk,vv in vvv.items(): # Iterate over positions
                for k,v in vv.items(): # Iterate over tubes
                    data_ti = extract_data_v1(v)
                    if icol==0:
                        data_t = data_ti.copy()
                        icol=1
                    else:
                        data_t = np.concatenate((data_t,data_ti))
    return data_t


def contour_2scale(ax,xx,yy,zz,labels=['','',r'$Densisty$'],zlims = [0,0.1],n_levels=10,ylim=[0,0.1]):
    
    im = ax.contour(xx,yy,zz,levels=np.linspace(zlims[0],zlims[1],n_levels))
    im1 = ax.contour(xx,yy,zz,cmap=cm.Reds)
    ax.set_xlabel(labels[0])
    ax.set_ylabel(labels[1])
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='5%', pad=0.1)
    cbar = plt.colorbar(im, ax=ax,cax=cax)
    plt.colorbar(im1, ax=ax,cax=cax)
    cbar.ax.set_ylabel(labels[2], rotation=270)
    ax.set_ylim(ylim)
    
def get_sel_elements(*dictionaries):
    new_dict = {}
    i = 0
    for el in dictionaries:
        for k,v in el.items():
            new_dict[str(i)] = v
            i = i + 1
    return new_dict

def get_directory_structure(rootdir,dxdt_path,exceptions):
    """
    Creates a nested dictionary that represents the folder structure of rootdir with tif files
    """
    
    xls = pd.ExcelFile(dxdt_path)
    dir = {}
    rootdir = rootdir.rstrip(os.sep)
    start = rootdir.rfind(os.sep) + 1
    it = sorted([(p, d, f) for p, d, f in os.walk(rootdir)], key=lambda x: x[0])
    for path, dirs, files in it:
        folders = path[start:].split(os.sep)
        files.sort()
        subdir = dict.fromkeys(files)
        for key,value in subdir.items():
            if '.tif' in key and not '.tiff' in key:
                full_path = path+'/'+key
                if not full_path in exceptions:
                    t,x,kymo = import_kymo(full_path)
                    k = folders[-2]
                    kk = folders[-1]
                    dxdt_df = pd.read_excel(xls, k)
                    i_pos = int(kk[kk.find('Pos')+len('Pos')])
                    i_tube = int(key[key.find('Tube')+len('Tube')])
                    aux = dxdt_df[(dxdt_df['Pos'] == i_pos) & (dxdt_df['Tube'] == i_tube)]
                    if not np.isnan(np.array(aux['Length (mm)'])):
                        Lx,dt = float(aux['Length (mm)']), float(aux['dt (min)'])
                        dx = Lx/kymo.shape[0]
                        subdir[key] = {'path':full_path, 't':t*dt,'x':x*dx,'Activity':kymo,'limits':[0,-1,0,-1]}
        parent = reduce(dict.get, folders[:-1], dir)
        parent[folders[-1]] = subdir
        
    my_dir = remove_none(dir)
    return my_dir

def remove_none(obj):
    if isinstance(obj, (list, tuple, set)):
        return type(obj)(remove_none(x) for x in obj if x is not None)
    elif isinstance(obj, dict):
        return type(obj)((remove_none(k), remove_none(v))
            for k, v in obj.items() if k is not None and v is not None)
    else:
        return obj
    

def add_pixels(my_dict,n_p,axis=1,pos=-1):
    array = my_dict['Activity']
    t,x = my_dict['t'],my_dict['x']
    n,m = array.shape
    array[abs(array) == np.inf] = np.nan
    avg_array = np.average(array[~np.isnan(array)])
    if axis == 0 and pos==-1:
        array0 = np.zeros((n_p,m)) + avg_array
        new_array = np.append(array0,array,axis=axis)
        x0 = np.arange(-n_p,0)*(x[1]-x[0])
        new_x = np.append(x0,x)
        new_t = t.copy()
    elif axis == 0 and pos==1:
        array0 = np.zeros((n_p,m)) + avg_array
        new_array = np.append(array,array0,axis=axis)
        x0 = np.arange(1,n_p+1)*(x[1]-x[0]) + x[-1]
        new_x = np.append(x,x0)
        new_t = t.copy()
    elif axis == 1 and pos==-1:
        array0 = np.zeros((n,n_p)) + avg_array
        new_array = np.append(array0,array,axis=axis)
        t0 = np.arange(-n_p,0)*(t[1]-t[0])
        new_t = np.append(t0,t)
        new_x = x.copy()
    elif axis == 1 and pos==1:
        array0 = np.zeros((n,n_p)) + avg_array
        new_array = np.append(array,array0,axis=axis)
        t0 = np.arange(1,n_p+1)*(t[1]-t[0]) + t[-1]
        new_t = np.append(t,t0)
        new_x = x.copy()

    new_dict = my_dict.copy()
    new_dict['Activity'] = new_array
    new_dict['t'] = new_t
    new_dict['x'] = new_x
    return new_dict

def bin_data(data_t,i_bin,i_var,n_bins=20,x=None):
    if x is None:
        x = np.linspace(np.min(data_t[:,i_bin]),np.max(data_t[:,i_bin]),n_bins)
        hist,x = np.histogram(data_t[:,i_bin], bins=n_bins)
        if (hist < 1).any():
            xmin,xmax = x[0],x[-1] 
            x = x[:-1]
            x_r,h_r = x[np.argmax(hist):],hist[np.argmax(hist):]
            x_l,h_l = x[:np.argmax(hist)],hist[:np.argmax(hist)]
            
            if (h_r<1).any():
                xmax = np.amin(x_r[h_r<1])
            if (h_l<1).any():
                xmin = np.amax(x_l[h_l<1])
            x = np.linspace(xmin,xmax,n_bins)

    y = [data_t[np.logical_and(data_t[:,i_bin] >=x[i], data_t[:,i_bin] < x[i+1]) ,i_var] for i in range(len(x)-1)]
    return x[:-1], y

def plot_scattercolor_dens(ax,x,y,datylim ,ylim=[0,1e3],s=5,loglog=[False,False], labels=['','',r'$Densisty$']):

    filter_d = np.where(np.logical_and(y < datylim[1], y > datylim[0]))
    x, y = x[filter_d],y[filter_d]

    xx, yy = np.meshgrid(np.linspace(x.min(),x.max(),100), np.linspace(y.min(),y.max(),100))
    # xx, yy = np.meshgrid(np.linspace(x.min(),x.max(),100), 1/np.linspace(1/y.max(),1/y.min(),100))
    positions = np.vstack([xx.ravel(), (yy).ravel()])
    xy = np.vstack([x,y])
    density = gaussian_kde(xy)
    z = density(xy)
    zz = np.reshape(density(positions).T, xx.shape)

    im = ax.scatter(x, y, c=z,s=s)
    ax.contour(xx,yy,zz)
    ax.set_ylim(ylim[0],ylim[1])
    ax.set_xlabel(labels[0])
    ax.set_ylabel(labels[1])
    cbar = plt.colorbar(im, ax=ax)
    cbar.ax.set_ylabel(labels[2], rotation=270)
    if loglog[0]:
        ax.set_xscale('log')
    if loglog[1]:
        ax.set_yscale('log') 
    return xx,yy,zz

def plot_violin(ax,x,y,labels=[r'$Period$ $(min)$',r'$Speed$ $(\mu m/min)$'],colors = ['#A4C2A5','#6C9D6E',1,'white','k','k']):
    parts = ax.violinplot(y,x,showmeans=False, showmedians=False,showextrema=False,widths=(x[1]-x[0])/2)
    for pc in parts['bodies']:
        pc.set_facecolor(colors[0])
        pc.set_edgecolor(colors[1])
        pc.set_alpha(colors[2])
        
    q1mq3 = []
    for el in y:
        q1mq3.append(np.percentile(el, [25, 50, 75]))
    q1mq3 = np.array(q1mq3)
    quartile1 = q1mq3[:,0]
    medians = q1mq3[:,1]
    quartile3 = q1mq3[:,2]

    whiskers = np.array([
        adjacent_values(sorted_array, q1, q3)
        for sorted_array, q1, q3 in zip(y, quartile1, quartile3)])
    whiskersMin, whiskersMax = whiskers[:, 0], whiskers[:, 1]

    ax.scatter(x, medians, marker='o', color=colors[3], s=30, zorder=3)
    ax.vlines(x, quartile1, quartile3, color=colors[4], linestyle='-', lw=3)
    ax.vlines(x, whiskersMin, whiskersMax, color=colors[5], linestyle='-', lw=1)
    ax.set_xlabel(labels[0])
    ax.set_ylabel(labels[1])
    return q1mq3

def fill_nans(indata,sx,sy):
    indata[np.isinf(indata)]=np.nan
    new_indata = indata.copy()
    # Find the non-NaN indices
    
    ii,jj = np.where(np.isnan(indata))
    for i,j in zip(ii,jj):
        li0,li1 = max(0,i-sy//2), min(i+sy//2,indata.shape[0])
        lj0,lj1 = max(0,j-sx//2), min(j+sx//2,indata.shape[1])
        neighbours = indata[li0:li1,lj0:lj1]
        new_indata[i,j] = np.average(neighbours[~np.isnan(neighbours)])
        #indata[li0:li1,lj0:lj1] =100
        
    return new_indata

def adjacent_values(vals, q1, q3):
    upper_adjacent_value = q3 + (q3 - q1) * 1.5
    upper_adjacent_value = np.clip(upper_adjacent_value, q3, vals[-1])

    lower_adjacent_value = q1 - (q3 - q1) * 1.5
    lower_adjacent_value = np.clip(lower_adjacent_value, vals[0], q1)
    return lower_adjacent_value, upper_adjacent_value

def import_kymo(file_path, print_info=False,plot=False):
    dataset = rasterio.open(file_path,transform= rasterio.Affine(1, 0, 0, 0, 1, 0))
    i_x = np.arange(dataset.width)
    i_y = np.arange(dataset.height)
    array = dataset.read(1)
    if print_info:
        print(dataset.count, dataset.height, dataset.width, dataset.crs)
    if plot:
        show(dataset)
    return i_x,i_y,array

def fill_nans(indata,sx,sy):
    indata[np.isinf(indata)]=np.nan
    new_indata = indata.copy()
    # Find the non-NaN indices
    
    ii,jj = np.where(np.isnan(indata))
    for i,j in zip(ii,jj):
        li0,li1 = max(0,i-sy//2), min(i+sy//2,indata.shape[0])
        lj0,lj1 = max(0,j-sx//2), min(j+sx//2,indata.shape[1])
        neighbours = indata[li0:li1,lj0:lj1]
        new_indata[i,j] = np.average(neighbours[~np.isnan(neighbours)])
        #indata[li0:li1,lj0:lj1] =100
        
    return new_indata


def import_kymo(file_path, print_info=False,plot=False):
    dataset = rasterio.open(file_path,transform= rasterio.Affine(1, 0, 0, 0, 1, 0))
    i_x = np.arange(dataset.width)
    i_y = np.arange(dataset.height)
    array = dataset.read(1)
    if print_info:
        print(dataset.count, dataset.height, dataset.width, dataset.crs)
    if plot:
        show(dataset)
    return i_x,i_y,array

def fill_nans(indata,sx,sy):
    indata[np.isinf(indata)]=np.nan
    new_indata = indata.copy()
    # Find the non-NaN indices
    
    ii,jj = np.where(np.isnan(indata))
    for i,j in zip(ii,jj):
        li0,li1 = max(0,i-sy//2), min(i+sy//2,indata.shape[0])
        lj0,lj1 = max(0,j-sx//2), min(j+sx//2,indata.shape[1])
        neighbours = indata[li0:li1,lj0:lj1]
        new_indata[i,j] = np.average(neighbours[~np.isnan(neighbours)])
        #indata[li0:li1,lj0:lj1] =100
        
    return new_indata

# Figure 1

def custom_guassian_kde(dat_x,dat_y, nx=500,ny=500):
    sy = 0.2
    sx = 50
    xx = np.linspace(dat_x.min(),dat_x.max(),nx)
    yy = np.linspace(dat_y.min(),dat_y.max(),ny)

    dens = np.zeros((len(yy),len(xx)))
    for i,xi in enumerate(xx):
        weights_x = np.exp(-(dat_x - xi)**2/(2*sx*sx))/(np.sqrt(2*np.pi)*sx) 
        # kde = gaussian_kde(dat_y,sy,weights = weights_x)
        kde = gaussian_kde(dat_y,weights = weights_x)
        densy = kde(yy)
        dens[:,i] = densy/np.amax(densy)
    return xx,yy,dens


# Figure 2
def create_with_timecolormap(x,y,z,
                          cmap="plasma",
                          linewidth=3,
                          linestyle='solid',
                          norm = colors.Normalize(vmin=0, vmax=1)):
    
    points = np.array([x, y]).T.reshape(-1,1,2)
    segments = np.concatenate([points[:-1],points[1:]], axis=1)

    lc = LineCollection(segments,
                        cmap=cmap,
                        norm=norm,
                        linewidth=linewidth,
                        linestyle=linestyle)
    lc.set_array(z)
    return lc

def plot_period_dadt_speed(ax,folders,chi,
                           labelx = r'$\alpha$',
                           labelsy = [r'Period (min)',r'$\frac{dA}{dt}$ (a.u.)',r' Wave Speed ($\mu$m/min)'],
                           col = '#6A966B',
                           yaxis=True,twinxaxis=True,
                           xrange = [0,1],
                           yranges = [[0,150],[0,400],[0,150]],
                           inverse_param=False):
    with open(folders[0], 'rb') as f:
        param_scan = np.load(f)
        speeds = np.load(f)
        periods = np.load(f)
         
    with open(folders[1], 'rb') as f:
        p0 = np.load(f)
        period = np.load(f)
        amps = np.load(f)
        max_dadt = np.load(f)
    
    
    twin1 = ax.twinx()
    twin2 = ax.twinx()

    # Offset the right spine of twin2.  The ticks and label have already been
    # placed on the right by twinx above.
    twin2.spines.right.set_position(("axes", 1.3))

    if not inverse_param:
        p1, = ax.plot(param_scan,periods,c= col, label='Period')
        p2, = twin1.plot(p0,max_dadt,c= col,linestyle='dashed',label=r'$\frac{dA}{dt}$')
        p3, = twin2.plot(param_scan,abs(speeds)*chi*1e3,c= col,linestyle='dotted',label='Wave Speed')
    else:
        p1, = ax.plot(1/param_scan,periods,c= col, label='Period')
        p2, = twin1.plot(1/p0,max_dadt,c= col,linestyle='dashed',label=r'$\frac{dA}{dt}$')
        p3, = twin2.plot(1/param_scan,abs(speeds)*chi*1e3,c= col,linestyle='dotted',label='Wave Speed')
        
    ax.set(xlim=xrange, ylim=yranges[0], xlabel=labelx)
    twin1.set(xlim=xrange, ylim=yranges[1], xlabel=labelx)
    twin2.set(xlim=xrange, ylim=yranges[2], xlabel=labelx)

    if yaxis:
        ax.set(ylabel=labelsy[0])
    else:
        ax.axes.yaxis.set_ticklabels([])
    if twinxaxis:
        twin1.set( ylabel=labelsy[1])
        twin2.set( ylabel=labelsy[2])
    else:
        twin1.axes.yaxis.set_ticklabels([])
        twin2.axes.yaxis.set_ticklabels([])
        # twin1.axes.get_yaxis().set_visible(False)
        twin2.axes.get_yaxis().set_visible(False)
        twin2.axis('off')
    ax.legend(handles=[p1, p2, p3],frameon=True,loc='upper left',fontsize=10)


def plot_timeseries_fig2(ins, ax1,color='blue',ylim=[0,70]):
    y0=[10,40]
    # print(ins.params_dict)
    t, a, c = ins.time_series(y0=y0)
    # print(ins.tf_timeseries)
    
    ax1.plot(t,a,c=color,linewidth=2)
    # ax1.set_xlim(0,500)#t[-1])
    ax1.set_ylim(ylim[0],ylim[1])
    ax1.axis('off')

    
def plot_phaseplane_fig2(ins,ax,color='blue',xlim = [40,80],ylim=[0,70],axis=True):
    y0=[10,40]
    # print(ins.params_dict)
    t, a, c = ins.time_series(y0=y0)
    
    ax.plot(c,a,c=color)#,linewidth=3)
    ax.set_ylim(ylim[0],ylim[1])
    ax.set_xlim(xlim[0],xlim[1])
    if axis:
        ax.set_xlabel(r'cyclin B-Cdk1')
        ax.set_ylabel(r'Act. cyclin B-Cdk1')
    else:
        ax.axis('off')
    ax.tick_params(left = False, right = False , labelleft = False , 
                labelbottom = False, bottom = False)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    
    
def plot_kymo_figure2(ax,folder,chi,t_o2,colormap=None):
    phi_a, phi_i = import_evolution(folder)
    x, y, t = create_spacetime(folder)
    tt1 = writen_time_array_plus_one(folder)
    betat = import_time_array(folder  + 'betat.dat',n=len(tt1))
    ksx = import_field(folder  + 'ksx.dat')[0]




    divider = make_axes_locatable(ax)
    xax = divider.append_axes("bottom", size=1.2, pad=0.1)
    yax = divider.append_axes("left", size=1.2, pad=0.1)
    cax = divider.append_axes('right', size='3%', pad=0.05)

    if colormap is None:
        colormap = color_scale()

    im = ax.pcolormesh(t,x*chi,phi_a[:,0,:].T,cmap=colormap,vmax=70,rasterized=True)#, vmin=0,vmax=90)
    
    # colormap = custom_colormap(['white','#E4EBF1','#5c87ad','#19354D'])
    # im = ax.pcolormesh(t,x*chi,phi_a[:,0,:].T,cmap=colormap,vmax=90,rasterized=True)#, vmin=0,vmax=90)
    
    ax.axes.xaxis.set_ticklabels([])
    ax.axes.yaxis.set_ticklabels([])
    ax.set_xlim(t[0],t[-1])
    xax.set_xlim(t[0],t[-1])
    yax.set_ylim(x[0]*chi,x[-1]*chi)

    cbar = plt.colorbar(im, cax=cax, orientation='vertical')
    cbar.set_label(r'$Act.$ $cyclin$ $B-Cdk1$', rotation=270,labelpad=20)

    blueish = '#5c87ad'
    xax.plot(tt1,1/betat[::len(betat)//len(tt1)],c=blueish)
    xax.set_xlabel(r'$t$ $(min)$',fontsize=int(14))
    xax.set_ylabel(r'$1/\beta$ $(Arb.)$',fontsize=int(14))
    xax.yaxis.tick_right()
    xax.yaxis.set_label_position("right")
    xax.axvline(t_o2,c='k')
    xax.text(t_o2*1.05,3.5, r'$t_{o2}$',fontsize=int(14))
    ymin, ymax = xax.get_ylim()
    yfill = np.linspace(ymin, ymax,10)
    xfill = np.zeros(len(yfill))+t_o2
    xax.fill_betweenx(yfill,xfill,facecolor='#DADDDD')
    xax.set_ylim(ymin,ymax)
    
    reddish = '#C44536'
    yax.plot(ksx,x*chi,c=reddish)
    yax.set_ylabel(r'$x$ $(mm)$',fontsize=int(14))
    yax.set_xlabel('$k_s$ $(nM/min)$',fontsize=int(14))
    yax.invert_xaxis()
    yax.text(1.525,6, r'$k_s(1$',rotation='vertical')
    yax.text(1.525,6.5, r'$+P(x)$',rotation='vertical',c='gray')
    yax.text(1.525,7.3, r'$+A_{ks}N(x))$',rotation='vertical',c=reddish)
    yax.axvline(np.average(ksx),c='k',linestyle = (0,(1,10)))

    
    t0_x = 2.5e-3*np.exp(-(x*chi-3)**2/(2*1**2)) \
         + 1.e-3*np.exp(-(x*chi-7.)**2/(2*1**2)) \
         + 1e-2*np.exp(-(x*chi-10)**2/(2*1**2))
    t0_x = t0_x - np.average(t0_x)
    ks_plot = 1.5
    yax.plot(ks_plot*(1+2*t0_x),x*chi,c='gray') 
    
def plot_simple_kymo(ax,folder,chi,plot_traj='None',color='orange',colormap=None):
    phi_a, phi_i = import_evolution(folder)
    x, y, t = create_spacetime(folder)
    tt1 = writen_time_array_plus_one(folder)

    if colormap is None: 
        colormap = color_scale()

    im = ax.pcolormesh(t,x*chi,phi_a[:,0,:].T,cmap=colormap,vmax=70,rasterized=True)#, vmin=0,vmax=90)

    # colormap = custom_colormap(['white','#E4EBF1','#5c87ad'])
    # colormap = custom_colormap(['white','#E4EBF1','#5c87ad','#19354D'])
    # im = ax.pcolormesh(t,x*chi,phi_a[:,0,:].T,cmap=colormap,vmax=90,rasterized=True)#, vmin=0,vmax=90)

    ax.set_xlim(t[0],t[-1])

    ax.set_xlabel(r'$t$ $(min)$',fontsize=18)
    ax.set_ylabel(r'$x$ $(mm)$',fontsize=18)
    ax.set_yticks([0,2,4,6,8])
    
    if plot_traj == 'All':
        waves_traj = get_traj(t,x,phi_a[:,0,:])
        if folder == '../sim/mit-waves18/':
            for ele in waves_traj[:-1]:
                # ax.plot(ele[0],ele[1]*chi,c='white',linewidth=5)
                ax.plot(ele[0],ele[1]*chi,c=color,linewidth=2)
        else:
            for ele in waves_traj:
                # ax.plot(ele[0],ele[1]*chi,c='white',linewidth=5)
                ax.plot(ele[0],ele[1]*chi,c=color,linewidth=2)
    elif type(plot_traj) is int:
        waves_traj = get_traj(t,x,phi_a[:,0,:])
        ele = waves_traj[plot_traj].copy()
        # ax.plot(ele[0],ele[1]*chi,c='white',linewidth=5)
        ax.plot(ele[0],ele[1]*chi,c=color,linewidth=2)
        
def get_dispersion_relation(folder,chi,plot=False,smooth_coeff = 1./3,resize = 4096):
    phi_a, phi_i = import_evolution(folder)
    x, y, t = create_spacetime(folder)
    nx = len(x)
    limits = [0,-1,0,-1]
    array = phi_a[:,0,::nx//resize].T
    x = x[::nx//resize]*chi*1000
    t = t[:]
    data_theo = extract_data_v1({'t':t,'x':x,'Activity':array,'limits':limits},order=150)
    
    wnx_period,wny_period = sm.nonparametric.lowess(data_theo[:,2],data_theo[:,1],  frac=smooth_coeff*1500/t[-1], it=1, return_sorted = True).T
    wnx_slope,wny_slope = sm.nonparametric.lowess(1/(data_theo[:,3]),data_theo[:,1],  frac=smooth_coeff*1500/t[-1], it=1, return_sorted = True).T
    if plot:
        plt.plot(data_theo[:,1],data_theo[:,2],'.',c='gray',alpha=0.5)
        plt.plot(wnx_period,wny_period,c='orange')
        plt.show()
        plt.plot(data_theo[:,1],1/data_theo[:,3],'.',c='gray',alpha=0.5)
        plt.plot(wnx_slope,wny_slope,c='orange')
        plt.show()
        plt.plot(data_theo[:,2],1/data_theo[:,3],'.',c='gray',alpha=0.5)
        plt.plot(wny_period,wny_slope,c='orange')
        plt.show()
    return wnx_period, wny_period, wnx_slope, wny_slope
        

# Figure 5
def get_speed_phase_triggerv0(folder,chi,plot=False,minlengthwave = 40,slopelim = 1e-2):
    '''Returns the speed of phase and trigger waves from a kymograph (simulation) as input (folder)'''
    
    phi_a, phi_i = import_evolution(folder)
    x, y, t = create_spacetime(folder)
    nx = len(x)
    limits = [0,-1,0,-1]
    array = phi_a[:,0,:]
    x = x*chi*1000
    
    waves_traj = get_traj(t,x,array)
    # minlengthwave = int(0.05*len(x))
    waves_traj = [el for el in waves_traj if len(el[0]) > 2*minlengthwave ]    # Ensure the wave is long enough
    waves = spl_waves_traj(waves_traj,t[1]-t[0])

    # fits = np.zeros((len(waves_traj),6))
    fits = []
    for i,el in enumerate(waves):
        xx,tt = el[1,len(el[0])//2:] - el[1,len(el[0])//2] ,el[0,len(el[0])//2:] - el[0,len(el[0])//2]

        slope = abs(np.gradient(tt,1)/np.gradient(xx,1))
        xx,tt = xx[slope > slopelim],tt[slope > slopelim]
        xx,tt = xx[~np.isnan(xx)],tt[~np.isnan(xx)] # Clean x nans
        xx,tt = xx[~np.isnan(tt)],tt[~np.isnan(tt)] # Clean t nans
        xx,tt = xx[~np.isinf(xx)],tt[~np.isinf(xx)] # Clean x infs
        xx,tt = xx[~np.isinf(tt)],tt[~np.isinf(tt)] # Clean t infs
        xx,tt = xx[tt<t[-1]],tt[tt<t[-1]] # Clean t infs
        xx,tt = xx[tt>1e-5],tt[tt>1e-5] # Clean t infs
        xx,tt = xx[:-2],tt[:-2]

        # plt.show()
        if len(xx)>minlengthwave and len(tt)>minlengthwave:
            popt, pcov = curve_fit(func_fit_speed, tt,xx,bounds=([0,0,0], [1000,10,1000]))
            perr = np.sqrt(np.diag(pcov))
            # fits[i] = np.concatenate((np.array(popt),np.array(perr)))
            fits.append(np.concatenate((np.array(popt),np.array(perr))))
        if plot:
            plt.plot(tt,xx,'.')
            ttt = np.linspace(np.amin(tt),np.amax(tt),100)
            if len(xx)>minlengthwave and len(tt)>minlengthwave:
                plt.plot(ttt,func_fit_speed(ttt, *popt), 'r--')
            plt.show()
    fits = np.array(fits).T
    # plt.show()
    return fits

def get_speed_phase_triggerv3(folder,chi,plot=False,slopelim = 1e-2):
    '''Returns the speed of phase and trigger waves from a kymograph (simulation) as input (folder)'''
    phi_a, phi_i = import_evolution(folder)
    x, y, t = create_spacetime(folder)
    nx = len(x)
    limits = [0,-1,0,-1]
    array = phi_a[:,0,:]
    x = x*chi*1000
    
    waves_traj = get_traj(t,x,array)
    waves = spl_waves_traj(waves_traj,t[1]-t[0])

    
    # minlengthwave = int(0.05*len(x))

            
    lastlast_wave = waves[-2]
    # for i,el in enumerate(waves):
    #     plt.plot(el[0],el[1])
    #     plt.axvline(max(el[0]))
    #     plt.axvline(min(lastlast_wave[0]),c='r')
    # plt.show()
    waves = [el for el in waves if max(el[0]) <min(lastlast_wave[0])  ]    # Ensure the wave is long enough

    # fits = np.zeros((len(waves),6))
    fits = []
    for i,el in enumerate(waves):
        # xx,tt = el[1,len(el[0])//2:] - el[1,len(el[0])//2] ,el[0,len(el[0])//2:] - el[0,len(el[0])//2]
        xx,tt = el[1] ,el[0]


        slope = abs(np.gradient(tt,1)/np.gradient(xx,1))
        xx,tt = xx[slope > slopelim],tt[slope > slopelim]
        xx,tt = xx[~np.isnan(xx)],tt[~np.isnan(xx)] # Clean x nans
        xx,tt = xx[~np.isnan(tt)],tt[~np.isnan(tt)] # Clean t nans
        xx,tt = xx[~np.isinf(xx)],tt[~np.isinf(xx)] # Clean x infs
        xx,tt = xx[~np.isinf(tt)],tt[~np.isinf(tt)] # Clean t infs
        xx,tt = xx[tt<t[-1]],tt[tt<t[-1]] # Clean t infs
        xx,tt = xx[tt>1e-5],tt[tt>1e-5] # Clean t infs
        xx,tt = xx[:-2],tt[:-2]

        if len(xx)>1:
            # popt, pcov = curve_fit(func_fit_speed, tt,xx,bounds=([0,0,0], [1000,10,1000]))
            # perr = np.sqrt(np.diag(pcov))
            # fits.append(np.concatenate((np.array(popt),np.array(perr))))
            slope1 = abs(np.gradient(tt,1)/np.gradient(xx,1))
            # fits.append([1/np.average(slope1),0,0,0,0,0])
            fits.append([1/np.average(slope1),1/np.average(slope1)**2*np.std(slope1),np.amin(slope1),np.amax(slope1),np.average(slope1),np.std(slope1)])
            if plot:
                plt.plot(tt+el[0,len(el[0])//2],xx,'.')
                ttt = np.linspace(np.amin(tt),np.amax(tt),100)
                # plt.plot(ttt+el[0,len(el[0])//2],func_fit_speed(ttt, *popt), 'r--')
                # plt.show()
    fits = np.array(fits).T
    if plot:
        plt.show()
    return fits

def func_fit_speed(t,c,gamma,d):
    x = d + c*t**gamma
    return x

def plot_alphascan(ax,axy,folder,chi,alphas,colors= ['red','blue']):
    subfolders = [folder + 'sim-{}/'.format(i) for i in range(len(glob.glob(folder + '*')))]

    fits_p = []
    for i,subf in enumerate(subfolders):
        fits_p.append(get_speed_phase_triggerv0(subf,chi))
    fits_p_pacemaker_gamma = fits_p.copy()
    avg_fits_pacemaker_gamma = np.array([(np.average(el[0],weights = 1/el[3]**2),np.average(el[1],weights = 1/el[4]**2)) for el in fits_p_pacemaker_gamma]).T
    avg_fits_pacemaker_gamma = np.array([(el[0,24],el[1,24]) for el in fits_p_pacemaker_gamma]).T

    cmap = custom_colormap(colors)  
    
    gammalim = [0,2]
    alpha_thr = np.min(alphas[abs(avg_fits_pacemaker_gamma[1] - 1) < 1e-2])
    axy.fill_between(np.linspace(alpha_thr,1.2,10),np.zeros(10)+gammalim[0],np.zeros(10)+gammalim[1],color=cmap(0.1),alpha=0.3)
    axy.axvline(alpha_thr,color=cmap(0.8))
    # axy.plot(np.linspace(alpha_thr,1.2,10),np.zeros(10)+1,color=cmap(0.8))
    
    lc0 = create_with_timecolormap( alphas, avg_fits_pacemaker_gamma[0], alphas , cmap = cmap)
    ax.add_collection(lc0)
    lc1 = create_with_timecolormap( alphas, avg_fits_pacemaker_gamma[1], alphas , cmap = cmap,linestyle='dashed')
    axy.add_collection(lc1)
    
    # cmap1 = custom_colormap(['white',colors[1]]) 
    # lc2 = create_with_timecolormap( alphas, alphas*0+80, alphas**2 , cmap = cmap1,linewidth = 20)
    # ax.add_collection(lc2)
    
    # cmap2 = custom_colormap(['white','#cf3a82']) 
    # lc3 = create_with_timecolormap( alphas, alphas*0+90, alphas , cmap = cmap2,linewidth = 20)
    # ax.add_collection(lc3)
    # ax.plot(np.linspace(0.1+0.0275,1-0.0275,len(alphas)),alphas*0+90,linewidth = 20,c='#F3CEE0')
    
    ax.plot(alphas[0],avg_fits_pacemaker_gamma[0,0],'o',c=cmap(0.3))
    axy.plot(alphas[0],avg_fits_pacemaker_gamma[1,0],'o',c=cmap(0.3), fillstyle='none')
    ax.plot(alphas[-1],avg_fits_pacemaker_gamma[0,-1],'o',c=colors[1])
    axy.plot(alphas[-1],avg_fits_pacemaker_gamma[1,-1],'o',c=colors[1], fillstyle='none')
    
    ax.set_ylim(0,80)
    ax.set_xlim(0,1.03)
    axy.set_ylim(gammalim)
    ax.set_xlabel(r'Timescale sep. $\alpha$',fontsize=16)
    ax.set_ylabel(r'$v$',fontsize=16)
    axy.set_ylabel(r'$\gamma$',fontsize=16)
    line0 = Line2D([0], [0], label=r'$v$', color=cmap(0.5))
    line1 = Line2D([0], [0], label=r'$\gamma$', color=cmap(0.5),linestyle='dashed')

    
    handles, labels = ax.get_legend_handles_labels()

    #add legend
    handles.extend([line0, line1])
    ax.legend(frameon=False,fontsize=16,handles=handles,title=r'$x = d+vt^\gamma$',title_fontsize=18,loc='lower right')

def plot_phaseplane(ax,folder,color='blue',xlim = [40,80],ylim=[0,70],axis=True):
    phi_a, phi_i = import_evolution(folder)
    ax.plot(phi_a[:,0,0] + phi_i[:,0,0],phi_a[:,0,0],c=color,linewidth=3)
    ax.set_ylim(ylim[0],ylim[1])
    ax.set_xlim(xlim[0],xlim[1])
    if axis:
        ax.set_xlabel(r'cyclin B-Cdk1',fontsize = 16)
        ax.set_ylabel(r'Act. cyclin B-Cdk1',fontsize = 16)
    else:
        ax.axis('off') 
        
def plot_phaseplane_color_timescale_sep(ax,folder,color='blue',xlim = [40,80],ylim=[0,70],axis=True,plot_cbar=False):
    phi_a, phi_i = import_evolution(folder)
    # ax.plot(phi_a[:,0,0] + phi_i[:,0,0],phi_a[:,0,0],c=color,linewidth=3)

    cmap = custom_colormap(['#D2DECF','#166D2D','#166D2D'])  
    dadt = abs(np.gradient(phi_a[:,0,0]))
    lc0 = create_with_timecolormap( phi_a[:,0,0] + phi_i[:,0,0], phi_a[:,0,0], dadt , cmap = cmap, norm = colors.Normalize(vmin=0, vmax=5))
    ax.add_collection(lc0)
    ax.set_ylim(ylim[0],ylim[1])
    ax.set_xlim(xlim[0],xlim[1])
    if axis:
        ax.set_xlabel(r'cyclin B-Cdk1',fontsize = 16)
        ax.set_ylabel(r'Act. cyclin B-Cdk1',fontsize = 16)
    else:
        ax.axis('off')
    if plot_cbar:
        divider = make_axes_locatable(ax)
        cax = divider.append_axes('bottom', size='5%', pad=0.05)
        cbar = plt.colorbar(lc0,cax=cax,location='bottom',orientation='horizontal')
        cbar.set_label(r'$\frac{da}{dt}$',fontsize=26)

def plot_evol_speed(ax,folder,chi,color='blue',plot_yaxis=True):
    fit = get_speed_phase_triggerv3(folder,chi)
    ax.plot(fit[0,0:25],c=color)
    # ax.plot(fit[0,0:25],'o-',c=color)
    # ax.fill_between(np.arange(0,26),fit[0,0:26] - fit[1,0:26]/2,fit[0,0:26] + fit[1,0:26]/2,color=color,alpha=0.3)
    ax.set_ylim(0,125)
    ax.set_xlabel(r'Cycle number',fontsize=16)
    if plot_yaxis:
        ax.set_ylabel(r'Avg. Speed ($\mu m/min$)',fontsize=16)
    else:
        ax.axes.yaxis.set_ticklabels([])
        
def plot_timeseries_fig5(ax1,folder,color='blue',xlim = [40,80],ylim=[0,70]):
    phi_a, phi_i = import_evolution(folder)
    x, y, t = create_spacetime(folder)
    
    ax1.plot(t,phi_a[:,0,0],c=color,linewidth=2)
    ax1.set_xlim(0,500)#t[-1])
    ax1.set_ylim(ylim[0],ylim[1])
    ax1.axis('off')
    
def plot_kymo_fig5_1(ax,yax,ax1,cax,folder,chi,color='blue',plot_ictau=True,plot_cbar=True,plot_xaxis=True,plot_yaxis=True,plot_ic=False,plot_waves=False,tlim=1400,colormap=None):
    phi_a, phi_i = import_evolution(folder)
    x, y, t = create_spacetime(folder)
    taux = import_field(folder  + 'taux.dat')[0]

    if colormap is None:
        colormap = color_scale()
    # colormap = custom_colormap(['white','#E4EBF1','#5c87ad','#19354D'])

    
    im = ax.pcolormesh(t,x*chi,phi_a[:,0,:].T,cmap=colormap,vmax=70,rasterized=True)#, vmin=0,vmax=90)
    ax.set_xlim(0,500)
    im1 = ax1.pcolormesh(t,x*chi,phi_a[:,0,:].T,cmap=colormap,vmax=70,rasterized=True)#, vmin=0,vmax=90)
    ax1.set_xlim(tlim-500/6,tlim)
    
    yax.set_ylim(x[0]*chi,x[-1]*chi)
    ax1.axes.yaxis.set_ticklabels([]) 
    ax1.set_xticks([tlim])
    ax.spines['right'].set_visible(False)
    ax1.spines['left'].set_visible(False)
    ax1.tick_params(left = False , labelleft = False)
    # yax1 = yax.twiny()
    if plot_xaxis:
        ax.set_xlabel(r'Time (min)', fontsize=16) 
    else:
        ax.axes.xaxis.set_ticklabels([])
        ax1.set_xticklabels([])
    ax1.set_xticklabels([str(tlim)])

    if plot_yaxis:
        ax.set_ylabel(r'Space (mm)', fontsize=16) 
    else:
        ax.axes.yaxis.set_ticklabels([]) 
    if plot_cbar:
        divider = make_axes_locatable(ax)
        # cax = divider.append_axes('right', size='5%', pad=0.05)
        cbar = plt.colorbar(im, cax=cax, orientation='vertical')
        cbar.set_label(r'Act. cyclin B-Cdk1', rotation=270,labelpad=20, fontsize=16)
    if plot_ictau:
        if np.max(taux) - np.min(taux)< 1e-5:
            yax.plot(0.+np.zeros(len(x)),x*chi,c='black')
        else:
            # print((taux - np.min(taux))/(np.max(taux) - np.min(taux)))
            yax.plot((taux - np.min(taux))/(np.max(taux) - np.min(taux)),x*chi,c='black')
        # yax.set_xlabel(r'$\tau(x)$',fontsize=int(size*2))
        yax.invert_xaxis()
        yax.axis('off')
    if plot_ic:

        ai_max,ai_min = max(np.max(phi_a[0,0,:]),np.max(phi_i[0,0,:])),min(np.min(phi_a[0,0,:]),np.min(phi_i[0,0,:]))
        if ai_max - ai_min < 1e-5:
            yax.plot(1+0.5+np.zeros(len(x)),x*chi,c='#5C87AD')
            yax.plot(1+0.5+np.zeros(len(x)),x*chi,c='#5C87AD',linestyle='dashed')  
        else:
            yax.plot(1+(phi_a[0,0,:]-ai_min)/(ai_max-ai_min),x*chi,c='#5C87AD')
            yax.plot(1+(phi_i[0,0,:]-ai_min)/(ai_max-ai_min),x*chi,c='#5C87AD',linestyle='dashed')
        # yax.set_xlabel(r'$\tau(x)$',fontsize=int(size*2))
        # yax.invert_xaxis()
        yax.axis('off')
    yax.set_xlim(2,0)
    if plot_waves:
        nx = len(x)
        limits = [0,-1,0,-1]
        array = phi_a[:,0,:]
        x_u = x*chi*1000

        waves_traj = get_traj(t,x_u,array)

        lastlast_wave = waves_traj[-2]

        waves_traj = [el for el in waves_traj if max(el[0]) <min(lastlast_wave[0])  ]    # Ensure the wave is long enough

        for i,el in enumerate(waves_traj):
            xx,tt = el[1],el[0]

            slope = abs(np.gradient(tt,1)/np.gradient(xx,1))
            xx,tt = xx[slope > 1e-2],tt[slope > 1e-2]
            xx,tt = xx[~np.isnan(xx)],tt[~np.isnan(xx)] # Clean x nans
            xx,tt = xx[~np.isnan(tt)],tt[~np.isnan(tt)] # Clean t nans
            xx,tt = xx[~np.isinf(xx)],tt[~np.isinf(xx)] # Clean x infs
            xx,tt = xx[~np.isinf(tt)],tt[~np.isinf(tt)] # Clean t infs
            xx,tt = xx[tt<t[-1]],tt[tt<t[-1]] # Clean t infs
            xx,tt = xx[tt>1e-5],tt[tt>1e-5] # Clean t infs
            xx,tt = xx[2:-2],tt[2:-2]

            if len(xx)>1:
                ax.plot(tt,xx/1000,c=color,linewidth=3)
                ax1.plot(tt,xx/1000,c=color,linewidth=3)
                        

    
#Figure S1
def plot_kymo_figure_s1(ax,folder,chi,t_o2,colormap=None):
    phi_a, phi_i = import_evolution(folder)
    x, y, t = create_spacetime(folder)
    tt1 = writen_time_array_plus_one(folder)
    betat = import_time_array(folder  + 'betat.dat',n=len(tt1))
    ksx = import_field(folder  + 'ksx.dat')[0]

    divider = make_axes_locatable(ax)
    xax = divider.append_axes("bottom", size=0.9, pad=0.1)
    yax = divider.append_axes("left", size=0.9, pad=0.1)
    cax = divider.append_axes('right', size='3%', pad=0.05)

    if colormap is None:
        colormap = color_scale()
    im = ax.pcolormesh(t,x*chi,phi_a[:,0,:].T,cmap=colormap,vmax=70,rasterized=True)#, vmin=0,vmax=90)
    
    # colormap = custom_colormap(['white','#E4EBF1','#5c87ad','#19354D'])
    # im = ax.pcolormesh(t,x*chi,phi_a[:,0,:].T,cmap=colormap,vmax=90,rasterized=True)#, vmin=0,vmax=90)
    
    ax.axes.xaxis.set_ticklabels([])
    ax.axes.yaxis.set_ticklabels([])
    ax.set_xlim(t[0],1500)
    xax.set_xlim(t[0],1500)
    yax.set_ylim(x[0]*chi,x[-1]*chi)


    cbar = plt.colorbar(im, cax=cax, orientation='vertical')
    cbar.set_label(r'Act. cyclin B-Cdk1', rotation=270,labelpad=20)

    blueish = '#5c87ad'
    xax.plot(tt1,1/betat[::len(betat)//len(tt1)],c=blueish)
    xax.set_xlabel(r'Time (min)')
    xax.set_ylabel(r'$1/\beta$ (a.u.)')
    xax.yaxis.tick_right()
    xax.yaxis.set_label_position("right")
    xax.axvline(t_o2,c='k')
    # xax.text(t_o2*1.05,3.5, r'$t_{o2}$',fontsize=int(14))
    ymin, ymax = xax.get_ylim()
    yfill = np.linspace(ymin, ymax,10)
    xfill = np.zeros(len(yfill))+t_o2
    xax.fill_betweenx(yfill,xfill,facecolor='#DADDDD')
    xax.set_ylim(ymin,ymax)
    
    reddish = '#C44536'
    yax.plot(ksx,x*chi,c='black')
    yax.set_ylabel(r'Space (mm)')
    yax.set_xlabel('$k_s$ (nM/min)')
    yax.invert_xaxis()
    # yax.text(1.535,6, r'$k_s(1$',rotation='vertical')
    # yax.text(1.525,6.5, r'$+P(x)$',rotation='vertical',c='gray')
    # yax.text(1.535,6.5, r'$+A_{ks}N(x))$',rotation='vertical',c=reddish)
    # yax.axvline(np.average(ksx),c='k',linestyle = (0,(1,10)))
    # yax.set_xticklabels(yax.get_xticklabels(), rotation=45, ha='right')
    yax.set_xlim(1.55,1.45)
    yax.set_xticks([1.45,1.5,1.55])
    yax.set_xticklabels(['1.45','1.5','1.55'],rotation=90)
    


