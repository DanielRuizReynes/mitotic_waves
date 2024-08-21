import os
import f90nml
from collections import OrderedDict
import numpy as np
from scipy.io import FortranFile
import glob
import subprocess
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
from matplotlib import rc,rcParams,colors

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


def export_dim(folder, dims, n_nuclei, file_name = 'dim.h'): 
    if not os.path.isdir(folder):
        os.system('mkdir -p ' + folder)
    f = open(folder + file_name, "w")
    for el in dims:
        f.write('integer, parameter :: {}={}  !system size.\n'.format(el,dims[el]))
    f.write('integer, parameter :: n_nuclei={}  !number of nuclei.\n'.format(n_nuclei))
    
    f.close()

def export_param_sim(folder, params, nuc_params, params_int, params_init, file_name='param.in'):
    if not os.path.isdir(folder):
        os.system('mkdir ' + folder)
    file = folder + file_name 
    nml = f90nml.namelist.Namelist(OrderedDict({
    'model_params': OrderedDict(params),
    'nuclear_params': OrderedDict(nuc_params),
    'integration_params': OrderedDict(params_int),
    'init_params': OrderedDict(params_init)}))
    nml.write(file, force=True)

def import_dim(folder, file_name = 'dim.h'): #method of pde2model
    f1 = open(folder + file_name, "r")
    ll = f1.readlines()
    dim = {}
    for el in ll[:-1]:
        keyval = el.replace('integer, parameter :: ','').replace('=','').replace('!system size.\n','').replace(' ','')
        dim[keyval[0]] = int(keyval[1:])
    return  dim

def import_param_sim(folder, file_name='param.in'): #method of pde2model
    nml = f90nml.read(folder + file_name)
    return dict(nml['model_params'].todict()), dict(nml['nuclear_params'].todict()), dict(nml['integration_params'].todict()), dict(nml['init_params'].todict())

def create_spacetime(folder):
    dim = import_dim(folder)
    m,n = dim['m'],dim['n']
    model_params, nuclear_params, integration_params, init_params = import_param_sim(folder) 
    dx = integration_params['delta_x']
    dy = integration_params['delta_y']
    dt = integration_params['t_write']
    t_final = integration_params['t_final']
    n_t = int(round(t_final/dt))
    x = np.arange(n)*dx
    y = np.arange(m)*dy
    t = np.arange(n_t)*dt
    return x,y,t

def export_nuclei_params(folder, x_i, eps_i, s_i, file_name='nuclei_params.dat'): 
    file_path = folder + file_name
    f = FortranFile(file_path,'w')
    f.write_record(x_i)
    f.write_record(eps_i)
    f.write_record(s_i)
    f.close() 

def import_nuclei_params(folder, file_name='nuclei_params.dat'): 
    file_path = folder + file_name
    f = FortranFile(file_path,'r')
    x_i = f.read_reals('f8')
    eps_i = f.read_reals('f8')
    s_i = f.read_reals('f8')
    f.close()  
    return x_i, eps_i, s_i

def compile_f90(folder, code_path, machine='m',print_command = False):
    executables = glob.glob(folder + '*.x')
    if len(executables) != 0:
        for el in executables:
            os.system('rm ' + el)

    dims = import_dim(folder)
    dims_name = str(dims['n']) + 'x' + str(dims['m'])
    code_name, code_folder = get_file_folder(code_path)
    exe_name = folder + machine + code_name.replace('.f90','') + dims_name +'.x'
    modfiles = '~/Documents/library/modfiles/'
    dranxor_path = '~/Documents/library/dranxor.f90'
    
    # Machine specific command
    nuredduna_command = 'intel; ifort -ipo -O3 -no-prec-div -fp-model fast=2 -march=sandybridge -mtune=core-avx2 -qmkl -I '
    maxwell_command = 'ifort -fast -qmkl -I '
    
    use_command = maxwell_command
    if machine == 'n':
        use_command = nuredduna_command
        
    command = use_command +modfiles+ ' -I ' +folder+ ' -o ' +exe_name+ ' ' +code_path+ ' ' +dranxor_path
    if print_command:
        print(command)
    else:
        result = subprocess.run([command], shell=True, stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        if not result.stdout is None:
            print(result.stdout.decode("utf-8"))
        if not result.stderr is None:
            print(result.stderr.decode("utf-8"))

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

def plot_kymo_xt(folder,chi):
    phi_a, phi_i = import_evolution(folder)
    x, y, t = create_spacetime(folder)
    tt1 = writen_time_array_plus_one(folder)

    size = 7
    fig = plt.figure(figsize=(2*size,size),constrained_layout=True)
    ax = plt.subplot2grid((1, 1), (0, 0), colspan=60)

    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='5%', pad=0.05)

    colormap = color_scale()
    im = ax.pcolormesh(t,x*chi,phi_a[:,0,:].T,cmap=colormap)#, vmin=0,vmax=90)
    ax.set_xlim(t[0],t[-1])
    ax.set_xlim(t[0],t[-1])
    ax.set_ylim(x[0]*chi,x[-1]*chi)

    cbar = plt.colorbar(im, cax=cax, orientation='vertical')
    cbar.set_label(r'$Act.$ $cyclin$ $B-Cdk1$', rotation=270,labelpad=20)

    ax.set_ylabel(r'$x$ $(mm)$')
    ax.set_xlabel(r'$t$ $(min)$')

def import_evolution(folder,file_0 = 'a.dat', file_1 = 'i.dat'): #method of the specific model
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

def writen_time_array_plus_one(folder):
    model_params, nuclear_params, integration_params, init_params = import_param_sim(folder) 
    t_write = integration_params['t_write']
    tt  = np.arange(get_n_writen(folder)+1)*t_write
    return tt

def get_n_writen(folder):
    model_params, nuclear_params, integration_params, init_params = import_param_sim(folder) 
    delta_t, t_write, t_final = integration_params['delta_t'],integration_params['t_write'],integration_params['t_final']
    n_writen = round(t_final/t_write) 
    return n_writen

def color_scale(palette_colors = ['#402D77','#009ED6', '#009642','#FBE136','#E5382C']):
    index = np.linspace(0,1,len(palette_colors))
    indx_clr = list(zip(index,palette_colors))
    colormap = colors.LinearSegmentedColormap.from_list('blue_red', indx_clr, 256)
    return colormap

def nuclei_field_v(x, phi_a, x_i, eps_i, s_i, ec50_neb, n_neb):
    K_neb = ec50_neb**n_neb
    delta_x = x[1] - x[0]
    i_i = (x_i/delta_x).astype(int) + 1
    f_i =  K_neb/(K_neb + phi_a[0,i_i]**n_neb)
    
    v_ji = [- f_i[j]*eps_i[j]*np.exp(-((x-x_i[j])/s_i[j])**2/2) for j in range(len(x_i))]
    return np.sum(v_ji,axis=0) 

def f_i(phi_a,ec50_neb, n_neb):
    K_neb = ec50_neb**n_neb
    return   K_neb/(K_neb + phi_a**n_neb)

def get_random_index_for_nuclei(i0,i1,n, min_dist):
    i_vals = np.zeros(n, dtype='int')
    for j in range(n):
        distance = 0
        while distance < min_dist:
            random_i = np.random.randint(i0,i1,1)[0]
            distance = np.amin(abs(random_i - i_vals))
        i_vals[j] = random_i
    return np.sort(i_vals)