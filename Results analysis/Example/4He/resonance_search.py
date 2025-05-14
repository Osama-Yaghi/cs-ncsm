#! /usr/bin/python3
import sys
import os,errno
import re
import glob
import time
import subprocess
import fileinput
import warnings
import itertools
import pickle
#import matplotlib
import numpy                as np
import matplotlib.pyplot    as plt
from matplotlib             import colors as mcolors,rcParams,cm
from matplotlib.collections import LineCollection
from matplotlib.lines       import Line2D
from matplotlib.ticker      import AutoMinorLocator,MaxNLocator,NullLocator,MultipleLocator
from mpl_toolkits.mplot3d   import Axes3D
from scipy.optimize         import Bounds,minimize,minimize_scalar,root,shgo,dual_annealing,brute,approx_fprime
from scipy                  import ndimage
from scipy.stats import norm
import pandas as pd
rcParams['text.usetex'] = True
#rcParams['text.latex.unicode'] = True

pattern_JpiT  ={'search':re.compile(r"Searching.*J2 pi T2= *([0-9]+) +([-0-9]+) +([0-9]+)"),
                'fit':   re.compile(r"Fitting.*J2 pi T2= *([0-9]+) +([-0-9]+) +([0-9]+) +at E= +([-.0-9]+) +MeV;Gamma= +([-.0-9]+)")}
pattern_hint  ={'search':re.compile(r"E= *([-.0-9]+) *G= *([-.0-9]+)"),
                'fit':   re.compile(r"Compound basis dimension *([0-9]+)")}
pattern_bound =re.compile(r"E_min= *([-.0-9]+) *E_max= *([-.0-9]+) *G_min= *([-.0-9]+)")
pattern_iter  =re.compile(r"G_max= *([-.0-9]+) *Iter_max= *([-.0-9]+)")
pattern_resfil=re.compile(r"( *[-+0-9]+ *)")

ncsmc=None
iter_max=100
warnings.filterwarnings("ignore")

def initialize_ncsmc(mode='search'):
    global ncsmc,iter_max

    while True:
        line  = ncsmc.stdout.readline()
        #print(line.decode("utf-8"))
        result= re.search(pattern_JpiT[mode],line.decode("utf-8"))
        if result:
            J,Pi,T=(int(result.group(i+1)) for i in range(3))
            J=int(J/2)
            T=int(T/2)
            if mode == 'fit':
                fitted_E,fitted_G=(float(result.group(i+4)) for i in range(2))
        if 'Initial guess' in line.decode("utf-8"):
            break
        if 'Compound basis dimension' in line.decode("utf-8"):
            compound_egv=[]
            break
        if '*** all finished' in line.decode("utf-8"):
            return
    if mode == 'search':
        result=re.search(pattern_hint[mode] ,ncsmc.stdout.readline().decode("utf-8"))
        guessed_E,guessed_G=(float(result.group(i+1)) for i in range(2))
        result=re.search(pattern_bound,ncsmc.stdout.readline().decode("utf-8"))
        E_min,E_max,G_min  =(float(result.group(i+1)) for i in range(3))
        result=re.search(pattern_iter ,ncsmc.stdout.readline().decode("utf-8"))
        G_max,iter_max     = float(result.group(1)),int(result.group(2))
        run_dico={"guessed_E":guessed_E,"guessed_G":guessed_G,"E_min":E_min,"E_max":E_max,"G_min":G_min,"G_max":G_max}
    elif mode == 'fit':
        result=re.search(pattern_hint[mode] ,line.decode("utf-8"))
        compound_dimension=int(result.group(1))
        ncsmc.stdout.readline()
        for j in range(compound_dimension):
            compound_egv.append(float(ncsmc.stdout.readline().decode("utf-8")))
        ncsmc.stdout.readline()
        run_dico={'fitted_E':fitted_E,'fitted_G':fitted_G,'compound_egv':compound_egv}
    else:
        print('Unexpected interaction mode')
        raise
    return (J,Pi,T),run_dico

def start_ncsmc():
    global ncsmc

    file_pattern=['resonance*.out']
    files=[]
    for pattern in file_pattern:
        files.extend(glob.glob(pattern))
    for file_to_remove in files:
        os.remove(file_to_remove)

    try:
        ncsmc  = subprocess.Popen(['./ncsm_rgm_Am3_3_plus_Am2_1_1_rmatrix_ortg_omp_mpi.exe'], stdin=subprocess.PIPE,stdout=subprocess.PIPE)
    except:
        raise

def get_ncsmc_value(x,deg_idx=0,Er_fixed=None,Gr_fixed=None):
    global ncsmc
    nn,ny,nx=None,None,None

    if Er_fixed:
        x=[Er_fixed,float(x)]
    if Gr_fixed:
        x=[float(x),Gr_fixed]

    if np.shape(x)==(2,):
        Er_array=[x[0]]
        Gr_array=[x[1]]
    else:
        nn,ny,nx=np.shape(x)
        Er_array,Gr_array=x.reshape((2,nx*ny))

    E_delta=np.array([])
    for Er, Gr in zip(Er_array, Gr_array):
        Er_str='{:.16E}'.format(Er)
        Gr_str='{:.16E}'.format(Gr)
        ncsmc_in = Er_str.replace('E','d')+'  '+Gr_str.replace('E','d')+' \n'

        ncsmc.stdin.write(ncsmc_in.encode('utf-8'))
        ncsmc.stdin.flush()
        while True:
            ncsmc_out=ncsmc.stdout.readline()

            #print(ncsmc_out)
            result = re.findall(r"( *[-+.0-9E]{3,} *)",ncsmc_out.decode("utf-8"))
            if result:
                break

        ncsmc_out=result[deg_idx] if result else np.finfo('d').max

        E_delta=np.append(E_delta,[np.float64(ncsmc_out)])

    if nx and ny:
        E_delta=E_delta.reshape((ny,nx))

    return E_delta

def get_ncsmc_fit(x,l_over=None,mask=None):
    global ncsmc

    if not l_over:
        E_compound=x
    else:
        E_compound=np.zeros(len(mask))
        E_compound.flat[np.where(np.asarray(mask) == True)[0]] = x
        E_compound.flat[np.where(np.asarray(mask) == False)[0]]= l_over

    E_delta=np.array([])
    for E_c in E_compound:
        E_c_str='{:.16E}'.format(E_c)
        ncsmc_in = E_c_str.replace('E','d')+' \n'

        ncsmc.stdin.write(ncsmc_in.encode('utf-8'))
        ncsmc.stdin.flush()

    ncsmc_out=''.encode('utf-8')
    while True:
        ncsmc_tmp=ncsmc.stdout.readline()

        if 'Waiting for new compound energies' in ncsmc_tmp.decode("utf-8"):
            result = ncsmc_out.decode("utf-8")
            break
        else:
            ncsmc_out=ncsmc_tmp
        if ncsmc_out.decode("utf-8")=='': sys.exit()

    ncsmc_out=result if result else np.finfo('d').max

    E_delta=np.append(E_delta,[np.float64(ncsmc_out)])

#    print(E_compound,E_delta)
    return E_delta

def iterate_search():
    global ncsmc

    ncsmc.stdin.write('Next search \n'.encode('utf-8'))
    ncsmc.stdin.flush()

    return

def find_resonance():
    global ncsmc
    def generate_fun_constraint(res_in):
        def new_function(x):
            return np.linalg.norm(x-res_in)-1.5
        return new_function
    def generate_res_constraint(res_list_in):
        for elem in res_list_in:
            yield {'type':'ineq','fun':generate_fun_constraint(elem)}

    start_ncsmc()
    res_search={}


    while True:
        res_iter  =initialize_ncsmc()
        if res_iter:
            res_JPiT,res_Params=res_iter
        else:
            return res_search
        if res_JPiT in res_search.keys():
            if not res_Params in res_search[res_JPiT]['Params']:
                res_search[res_JPiT]['Params'].extend(res_Params)
        else:
            print('Doing J=%d Pi=%d and T=%d'%res_JPiT)
            res_search.update({res_JPiT:{'Params':[res_Params],'Resonances':[]}})

        print('Resonance is guessed in E_r in [%3.1f,%3.1f] and \Gamma_r in [%3.1f,%3.1f]'%(res_Params['E_min'],res_Params['E_max'],res_Params['G_min'],res_Params['G_max']))

        res_found= []
        bounds   = [(res_Params['E_min'],res_Params['E_max']),(res_Params['G_min'],res_Params['G_max'])]
        result   = shgo(get_ncsmc_value, bounds, minimizer_kwargs={'method':'Nelder-Mead'}, options={'local_iter':10,'infty_constraints':False,'f_tol':1e-5},iters=2, sampling_method='sobol')
        for i,elem in enumerate(result.xl):
            if all(bounds[j][0]<elem[j]<bounds[j][1] for j in range(len(elem))) and result.funl[i]<1e-4:
                res_search[res_JPiT]['Resonances'].extend([elem])
                res_found=elem
                print('Resonance was indeed found at E_r=%5.3f and \Gamma_r=%5.3f'%tuple(elem))
                print('The zero is numerically %f'%result.funl[np.where(result.xl==elem)[0][0]])
                break

#        if len(res_found)>1:
#            result  =brute(get_ncsmc_value,(slice(res_found[0]-0.3,res_found[0]+0.3,0.01),slice(0.,res_found[1]-0.3,0.1)),full_output=True,finish=None)
#            print(result[0],result[1])
        iterate_search()
    return

def fit_resonance():
    global ncsmc
    def generate_fun_constraint(res_in):
        def new_function(x):
            return np.linalg.norm(x-res_in)-1.5
        return new_function
    def generate_res_constraint(res_list_in):
        for elem in res_list_in:
            yield {'type':'ineq','fun':generate_fun_constraint(elem)}

    for line in fileinput.input('resonance.in',inplace=True, backup='.bak'):
        result=re.findall(pattern_resfil,line)
        if result:
            if len(result) == 5:
                line=line.replace(result[3].strip(),'4')
                print(line,end='')
            else:
                print(line,end='')
        else:
            print(line,end='')
    fileinput.close()

    start_ncsmc()
    res_fit={}


    while True:
        res_iter  =initialize_ncsmc(mode='fit')
        if res_iter:
            res_JPiT,res_Params=res_iter
        else:
            return res_fit
        if res_JPiT in res_fit.keys():
            if not res_Params in res_fit[res_JPiT]['Params']:
                res_fit[res_JPiT]['Params'].extend(res_Params)
            res_Params['compound_egv'] = np.asarray(res_fit[res_JPiT]['Compound energies'][-1])
        else:
            print('Doing J={0:d} Pi={1:d} and T={2:d}'.format(*res_JPiT))
            res_fit.update({res_JPiT:{'Params':[res_Params],'Compound energies':[],'Jacobian':[],'Mask':[True for i in range(len(res_Params['compound_egv']))]}})
            res_Params['compound_egv']=np.asarray(res_Params['compound_egv'])

        print('Resonance will be fitted at E_r {0:5.3f} and \Gamma_r {1:5.3}'.format(res_Params['fitted_E'],res_Params['fitted_G']))

        guess = res_Params['compound_egv'][res_fit[res_JPiT]['Mask']]
        l_over= list(set(res_Params['compound_egv']) - set(guess))
        result=minimize(get_ncsmc_fit, guess, args=(l_over,res_fit[res_JPiT]['Mask']), method='Nelder-Mead', options={'maxiter':100}, tol=1.e-4)
        if  result.fun<1e-4:
            if not l_over:
                new_compound=result.x
            else:
                new_compound=np.zeros(len(res_fit[res_JPiT]['Mask']))
                new_compound=flat[np.where(np.asarray(res_fit[res_JPiT]['Mask']) == True)[0]] = result.x.tolist()
                new_compound=flat[np.where(np.asarray(res_fit[res_JPiT]['Mask']) == False)[0]]= l_over
            res_fit[res_JPiT]['Compound energies'].append(new_compound.tolist())
            print('Resonance was indeed fitted at '+' and '.join(['E_c={:5.3f}'.format(E_c) for E_c in result.x.tolist()]))
            print('The zero is numerically {0:f}'.format(result.fun))
            grad=approx_fprime(result.x,get_ncsmc_fit,epsilon=np.sqrt(np.finfo(float).eps))
            print('Jacobian is computed to '+' and '.join(['df/dE_{0:1d}={1:5.3e}'.format(i+1,E_c) for i,E_c in enumerate(grad.tolist())]))
            res_fit[res_JPiT]['Jacobian'].append(grad.tolist())
            i_max=np.abs(grad).tolist().index(max(np.abs(grad)))
            res_fit[res_JPiT]['Mask'][i_max]=False

        iterate_search()
    return

def make_plots():
    global ncsmc

    start_ncsmc()
    res_search={}

    while True:
        res_iter  =initialize_ncsmc()
        if res_iter:
            res_JPiT,res_Params=res_iter
        else:
            return 'all done'
        if res_JPiT in res_search.keys():
            if not res_Params in res_search[res_JPiT]['Params']: res_search[res_JPiT]['Params'].extend(res_Params)
        else:
            res_search.update({res_JPiT:{'Params':[res_Params],'Resonances':[]}})
            plot_surface(res_JPiT)

        iterate_search()

def plot_surface(x):

    plt.clf()
    J,Pi,T=x

    plt.rc('text', usetex=True)
    plt.rc('font', **{'family':'sans-serif','sans-serif':['Helvetica']})
    plt.rcParams['text.latex.preamble']=[r"\usepackage{amsmath,amssymb}"]

    x  = np.arange( 0.001, 8,step=0.2)
    y  = np.concatenate((np.geomspace(0.001  , 2.0,num=20),np.arange(2.0  , 16,step=0.2)))
    xgrid, ygrid = np.meshgrid(x, y)
    xy = np.stack([xgrid, ygrid])

    fig = plt.figure()
    #ax = fig.add_subplot(111, projection='3d')
    ax = fig.add_subplot(111)
    ax.margins(0.04,.04)
    #ax.axis('equal')
    #ax.view_init(135, 45)

    zval=get_ncsmc_value(xy)
    #zval[zval>5e-1]= np.nan
    X = ndimage.zoom(xgrid, 4,order=3)
    Y = ndimage.zoom(ygrid, 4,order=3)
    Z = ndimage.zoom(zval , 4,order=3)
    Z[Z>5e-1]= np.nan

    #ax.plot_surface(xgrid, ygrid, zval, cmap='terrain')
    normalize = mcolors.LogNorm(vmin=0.0001, vmax=0.5)
    cf   =plt.contourf(X,Y,Z,levels=np.geomspace(0.0001, 0.5,num=10),norm=normalize, cmap=cm.YlOrBr_r)
    cbar =fig.colorbar(cf)

    ax.set_xlabel(r'$ \mathcal{R}e(E) \, [{\rm MeV}]$')
    ax.set_ylabel(r'$ \mathcal{I}m(E) \, [{\rm MeV}]$')
    #ax.set_zlabel(r'$\delta_E \, [{\rm MeV}]$')
    plt.grid(b=True, which='major', color='#666666', linestyle='-')
    plt.minorticks_on()
    plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
    #plt.show()

    outName='SurfPlot_J-'+str(J)+'_Pi-'+str(Pi)+'_T-'+str(T)

    fig.savefig(outName+'.pdf', bbox_inches='tight',dpi=100)
    fig.savefig(outName+'.eps')
    fig.savefig(outName+'.png', bbox_inches='tight', transparent=True,dpi=100)

    return

#visualization
def fig_size(nx,ny,scale=1):
    latex_width=438.17227*scale
    plot_width=1.0/72.27 *latex_width
    aesthetic_ratio=(5.0**0.5-1)/2*1.1
    plot_height=plot_width*aesthetic_ratio*nx/ny
    figsize=(plot_width,plot_height)
    return figsize
def plot_spectra(spectra_dic,style="width",Thresholds=None,plot_thresholds=False):
   #initialize some parameters
    plt.rcParams.update({'font.size': 12})
    if(style=="width"):
        y2scale=10
    elif(style=="error"):
        y2scale=1
    elif(style=="band"):
        y2scale=1
    if type(spectra_dic) is dict:
        spectra_dic=[spectra_dic]
    num_spectra=len(spectra_dic)
    if(plot_thresholds):
        num_spectra=num_spectra#+1

    plt.clf()

    plt.rc('text', usetex=True)
    plt.rc('font', **{'family':'sans-serif','sans-serif':['Helvetica']})
    #plt.rcParams['text.latex.preamble']=[r"\usepackage{amsmath,amssymb}"]
    plt.rcParams['text.latex.preamble']=r"\usepackage{amsmath,amssymb}"

    bar_width = 0.8

    data_E=np.array([])
    data_G=np.array([])
    data_Err=np.array([])
    data_colors=[]
    data_lcolors=[]
    data_hatch=[]
    data_width=[]
    data_zorder=[]
    data_JPiT=[]

    n_data=len(spectra_dic)
    n_rows=0
    cmap_names=['Greys', 'Purples', 'Blues', 'Greens', 'Oranges',
                'Reds','YlOrBr', 'RdPu','BuGn','YlGn']
    c_names=['grey', 'purple', 'blue', 'green', 'orange', 'red','saddlebrown',
             'pink','cyan','olive']
    patterns = ('-', '+', 'x', '\\', '*', 'o', 'O', '.','.','.')
    JPit_color={}
    max_used_color=-1
    ###
    for j,dic in enumerate(spectra_dic):
        data_E_tmp=np.array([])
        data_G_tmp=np.array([])
        data_Err_tmp=np.array([])
        data_JPiT_tmp=[]
        for JPiT in dic.keys():
            if(isinstance(JPiT,str)):
                continue
            data_JPiT_tmp.append(JPiT)
        data_JPiT_tmp=np.array(data_JPiT_tmp,dtype='i,i,i')
        data_JPiT_tmp=np.sort(data_JPiT_tmp).tolist()
        data_JPiT.extend(data_JPiT_tmp)
        data_colors_tmp=[]
        data_lcolors_tmp=[]
        data_hatch_tmp=[]
        for i,JPiT in enumerate(data_JPiT_tmp):
            if(isinstance(JPiT,str)):
                continue
            if(JPiT not in JPit_color):
                JPit_color[JPiT]=max_used_color+1
                max_used_color=max_used_color+1
            E_tmp=np.array([resonance[0] for resonance in dic[JPiT]['Resonances']])
            G_tmp=np.array([resonance[1] for resonance in dic[JPiT]['Resonances']])
            Err_tmp=np.array([resonance[2] for resonance in dic[JPiT]['Resonances']])
            if len(E_tmp) > 1:
                idx_E_sorted =np.argsort(E_tmp,axis=0)
            else:
                idx_E_sorted =0
            data_E_tmp=np.append(data_E_tmp,E_tmp[idx_E_sorted])
            data_G_tmp=np.append(data_G_tmp,G_tmp[idx_E_sorted])
            data_Err_tmp=np.append(data_Err_tmp,Err_tmp[idx_E_sorted])
            cmaptype=plt.cm.get_cmap(name=cmap_names[JPit_color[JPiT]])
            for rows in range(len(E_tmp)):
                data_colors_tmp.append(cmaptype((0.5+rows)/(len(E_tmp)+1.)))
                data_lcolors_tmp.append(c_names[JPit_color[JPiT]])
                data_hatch_tmp.append(patterns[JPit_color[JPiT]])
        n_rows=max(len(data_E_tmp),n_rows)
        n_last_rows=data_E.shape[0]

        data_width_tmp =[]
        data_zorder_tmp=[]
        data_width_tmp.extend( np.linspace(0.1,bar_width, len(data_E_tmp)))
        data_zorder_tmp.extend(np.linspace(0  ,len(data_E_tmp)-1 , len(data_E_tmp),dtype=np.int))
        tmp_w=np.ndarray(len(data_width_tmp) ,np.float64)
        tmp_o=np.ndarray(len(data_zorder_tmp),np.int)
        tmp_e=np.ndarray(len(data_E_tmp)     ,np.float64)
        tmp_w[:]=data_width_tmp[:]
        tmp_o[:]=data_zorder_tmp[:]
        tmp_e[:]=data_E_tmp[:]
        for i in range(len(data_E_tmp)):
            data_width_tmp[ np.argmin(tmp_e)]=tmp_w[i]
            data_zorder_tmp[np.argmin(tmp_e)]=tmp_o[-i]
            tmp_e[np.argmin(tmp_e)]=10000.

        data_E=np.resize(data_E.T,(j+1,n_rows)).T
        data_G=np.resize(data_G.T,(j+1,n_rows)).T
        data_Err=np.resize(data_Err.T,(j+1,n_rows)).T
        if n_rows > n_last_rows and n_last_rows >0:
            data_colors.extend([[(0,0,0,0)]*j]*(n_rows-n_last_rows))
            data_lcolors.extend(  [["None"]*j]*(n_rows-n_last_rows))
            data_hatch.extend(    [["None"]*j]*(n_rows-n_last_rows))
            data_width.extend(    [[0.]*j]    *(n_rows-n_last_rows))
            data_zorder.extend(   np.array([[100]*j]*(n_rows-n_last_rows)))
            for k in range(j):
                data_E[n_last_rows:,k]=np.array([0.])
                data_G[n_last_rows:,k]=np.array([0.])
                data_Err[n_last_rows:,k]=np.array([0.])
        elif len(data_E_tmp) < n_rows:
            current_rows=len(data_E_tmp)
            data_E_tmp.resize(n_rows)
            data_G_tmp.resize(n_rows)
            data_Err_tmp.resize(n_rows)
            data_E_tmp[current_rows+1:]=np.array([0.])
            data_G_tmp[current_rows+1:]=np.array([0.])
            data_Err_tmp[current_rows+1:]=np.array([0.])
            data_colors_tmp.extend([(0,0,0,0)]*(n_rows-current_rows))
            data_lcolors_tmp.extend(  ["None"]*(n_rows-current_rows))
            data_hatch_tmp.extend(    ["None"]*(n_rows-current_rows))
            data_width_tmp.extend(    [0.]    *(n_rows-current_rows))
            data_zorder_tmp.extend(   np.array([100]*(n_rows-current_rows)))

        data_E[:,-1]=data_E_tmp
        data_G[:,-1]=data_G_tmp
        data_Err[:,-1]=data_Err_tmp
        data_colors =[[*elem1,elem2] if not len(elem1)==4           else [elem1,elem2]  for elem1,elem2 in zip(data_colors,data_colors_tmp)   ] if data_colors  else data_colors_tmp
        data_lcolors=[[*elem1,elem2] if not isinstance(elem1, str)  else [elem1,elem2]  for elem1,elem2 in zip(data_lcolors,data_lcolors_tmp) ] if data_lcolors else data_lcolors_tmp
        data_hatch  =[[*elem1,elem2] if not isinstance(elem1, str)  else [elem1,elem2]  for elem1,elem2 in zip(data_hatch,data_hatch_tmp)     ] if data_hatch   else data_hatch_tmp
        data_width  =[[*elem1,elem2] if not isinstance(elem1, np.float) else [elem1,elem2]  for elem1,elem2 in zip(data_width,data_width_tmp)     ] if data_width   else data_width_tmp
        data_zorder =[[*elem1,elem2] if not isinstance(elem1, np.int64) else [elem1,elem2]  for elem1,elem2 in zip(data_zorder,data_zorder_tmp)   ] if data_zorder  else data_zorder_tmp
        #data_E=np.append(data_E,data_E_tmp)
        #data_G=np.append(data_G,data_G_tmp)

    #data_E=data_E.reshape(n_data,n_rows).T
    #data_G=data_G.reshape(n_data,n_rows).T

    figsize=fig_size(2,1,0.7)
    fig,ax=plt.subplots(1,1,figsize=figsize)
    #fig = plt.figure()
    #ax  = fig.add_subplot(131)
    ax.margins(0.04,.04)
    #ax.axis('equal')
    columns=[]
    #if(plot_thresholds):
    #    columns.append('thresholds')
    for i in range(n_data):
        j=i+1
        label=input('Choose a label for {:2d} th data (enter=none): '.format(i))
        if label: columns.append(r'$\mathrm{'+label+r'}$')
        else:
            columns.append('')

    index = np.arange(num_spectra)+bar_width/2.
    if(plot_thresholds):
        starting_index=0#1
    else:
        starting_index=0
    index=index[::-1]
    columns=columns[::-1]
    line_index= list(zip(index[:num_spectra-starting_index],index[:num_spectra-starting_index]+bar_width))
    lwidth= np.linspace(0.8,1.2, n_rows)

    s=input('do you want to apply shifts?y/n')
    if(s=='y'):
        s2=input('enter levels energy in one line')
        s3=input('enter shifts in one line')
        shifts=np.array([[float(loc),float(shft)] for loc,shft in
                         zip(s2.split(),s3.split())])
        print(shifts)
        data_E_shifts=data_E.copy()
        data_E_shifts[:,:]=0
        for loc,shft in shifts:
            mask=np.zeros(data_E.shape,dtype=bool)
            mask|= np.isclose(loc,data_E,atol=0.05)
            data_E_shifts[mask]=shft
            data_E_shifts[:,0]=0
    else:
        data_E_shifts=np.zeros(data_E.shape)
    #endif
    for row in range(n_rows):
        if(style=='band'):
            shift=0
        else:
            shift=-data_G[row]/(2*y2scale)
        plt.bar(index[:num_spectra-starting_index]+data_E_shifts[row], data_G[row][:num_spectra-starting_index]/y2scale , data_width[row], align='edge',alpha=0.7,lw=0.1,edgecolor=data_lcolors[row],bottom=data_E[row]+shift, color=data_colors[row],zorder=data_zorder[row][-1] if not isinstance(data_zorder[row],np.int64) else data_zorder[row])#, hatch=data_hatch[row])
        ax.errorbar(index[:num_spectra-starting_index]+data_width[row],data_E[row][:num_spectra-starting_index],yerr=data_Err[row][:num_spectra-starting_index],fmt='o',capsize=4,color='grey',markersize=0.01)
        #lines=[[(segment[0],value),(segment[1],value)] for segment,value in zip(line_index,data_E[row]) ]
        lines=[[(segment[0]+shft,value),(segment[1]+shft,value)] for
                   segment,value,shft in
                   zip(line_index,data_E[row],data_E_shifts[row]) ]
        for i,line in enumerate(lines):
            if data_lcolors[row][i] != 'None':
                spectra=LineCollection([line],linewidths=lwidth[row],colors=mcolors.to_rgb(data_lcolors[row][i] if not isinstance(data_lcolors[row],str) else data_lcolors[row]),zorder=n_rows-row)
                collection=ax.add_collection(spectra)


    #plot the thresholds
    if(plot_thresholds):
        threshold_labels=spectra_dic[0]['threshold_labels']
        gauss_range=np.linspace(-3,3,100)
        for j,spect in enumerate(spectra_dic):
            if( 'Thresholds' in spect):
                thresholds=np.array(list(spect['Thresholds']))
                first_threshold=thresholds[0,0]
                thresholds=[[thr[0]-first_threshold,thr[1]] for thr in
                            thresholds]
                for i,threshold in enumerate(thresholds):
                    ax.plot([index[j]-0.3, index[j]+bar_width+0.3],[threshold[0],threshold[0]],
                     color='black', linewidth=0.7,linestyle='--')
                    gaussian_distribution=norm.pdf(gauss_range,0,threshold[1])*0.1
                    if(j==0):
                        plt.annotate(threshold_labels[i],(index[j]+bar_width+0.1,threshold[0]),textcoords='offset points',xytext=(5,4),rotation=90,fontsize=6)
                    #ax.fill_betweenx( threshold[0] + gauss_range,index[j]+bar_width/2 -
                    #          gaussian_distribution,index[1]+bar_width/2,color='black',alpha=0.07)


    # Adjust layout to make room for the table:
    plt.subplots_adjust(left=0.2, bottom=0.2)

    ax.set_ylabel(r'$ E \, [{\rm MeV}]$')
    ax.set_xticks(index)
    ax.set_xticklabels(columns)
    ax.xaxis.set_tick_params(which='both',bottom=False,top=False)
    ax.yaxis.set_tick_params(which='both',direction= 'inout',left=True,right=False)
    #ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
    ax.yaxis.set_major_locator(MultipleLocator(1))
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.set_xlim(left=-0.5)
    ax.set_xlim(-0.5,len(spectra_dic)+0.7)
    ax.set_ylim(bottom=-0.5)
    ax.set_ylim(top=6.8)
    secax = ax.secondary_yaxis('right', functions=(lambda x: y2scale*x, lambda y: y/y2scale))
    if(style=='width'):
        secax.set_ylabel(r'$ \Gamma \, [ {\rm MeV}]$')
    secax.yaxis.set_tick_params(which='both',direction= 'inout',right=True)
    secax.yaxis.set_minor_locator(MultipleLocator(5))
    secax.yaxis.set_major_locator(MultipleLocator(y2scale*2.))
    plt.setp( ax.xaxis.get_majorticklabels(), rotation=-30, ha="left" )

    pi_str=lambda x: '+' if x>0 else '-'
    leg_handle=[Line2D([0], [0], color=c_names[JPit_color[stuff]], lw=1) for i,stuff in enumerate(sorted(set(data_JPiT)))]
#    leg_labels=['$J^\pi = {:2d}'.format(elem[0])+'^{'+pi_str(elem[1])+'}$' for elem in sorted(set(data_JPiT))]
    leg_labels=['${:2d}'.format(elem[0])+'^{'+pi_str(elem[1])+'};'+'{:2d}$'.format(elem[2]) for elem in sorted(set(data_JPiT))]
    legend= ax.legend(leg_handle,leg_labels,loc='upper left',prop={'size':10},frameon=False,ncol=2)
    outName=input("input the name for the output files")
    plt.tight_layout()
    plt.show()
    outName='spectra/Spectra_'+outName
    fig.savefig(outName+'.pdf', bbox_inches='tight', transparent=True,dpi=2000)
    fig.savefig(outName+'.png', bbox_inches='tight', transparent=True,dpi=2000)
    fig.savefig(outName+'.eps', bbox_inches='tight')
    fig.savefig(outName+'.svg', bbox_inches='tight', transparent=True,dpi=2000)

def expt_spectra(nuclei):
    dic_spectra={ '4He':{(0,-1,0):{'Resonances':[(1.195,0.84),(8.825,4.89)]},
                         (0, -1,1):{'Resonances':[(5.465,7.97)]},
                         (0, 1,0):{'Resonances':[(0.395,0.5)]},
                         (1,-1,0):{'Resonances':[(4.435,6.1),(8.555,3.92)]},
                         (1, -1,1):{'Resonances':[(3.825,6.2),(6.135,12.66)]},
                         (1, 1,0):{'Resonances':[(8.495,9.89)]},
                         (2,-1,0):{'Resonances':[(2.025,2.01),(8.575,8.75)]},
                         (2, -1,1):{'Resonances':[(3.515,5.01)]},
                         (2,1,0):{'Resonances':[(7.605,8.69),(8.855,3.78),(10.89,9.72)]},
                         'Thresholds':[[-8.48,0],[-7.718,0],[-4.448,0],[-2.224,0],[0,0]]
                         ,"threshold_labels":[r'$\mathrm{^3H\text{-}p}$',r'$\mathrm{^3He\text{-}n}$',r'$\mathrm{D\text{-}D}$',r'$\mathrm{D\text{-}np}$',r'$\mathrm{2n2p}$']
                        },
                 '4H':{(0,-1,1):{'Resonances':[(5.27,8.9)]},
                         (1,-1,1):{'Resonances':[(3.5,6.7),(6.02,13)]},
                         (2,-1,1):{'Resonances':[(3.19,4.6)]},
                         'Thresholds':[[-8.48,0],[-2.224,0],[0,0]]
                         ,"threshold_labels":[r'$\mathrm{^3H\text{-}n}$',r'$\mathrm{D\text{-}nn}$',r'$\mathrm{3np}$']
                         },
                 '4Li':{(0,-1,1):{'Resonances':[(6.15,9.35)]},
                         (1,-1,1):{'Resonances':[(4.39,7.35),(6.92,13.51)]},
                         (2,-1,1):{'Resonances':[(4.07,6.03)]},
                         'Thresholds':[[-7.718,0],[-2.224,0],[0,0]]
                         ,"threshold_labels":[r'$\mathrm{^3He\text{-}p}$',r'$\mathrm{D\text{-}pp}$',r'$\mathrm{n3p}$']
                         },
                 '10Be':{(0, 1,1):{'Resonances':[(-6.8122,0.0),(-0.6329,0.0)]},
                         (1,-1,1):{'Resonances':[(-0.8523,0.0)]},
                         (1, 1,1):{'Resonances':[( 3.7578,0.0)]},
                         (2, 1,1):{'Resonances':[(-3.4441,0.0),(-0.8538,0.0),(0.7298,0.0063),(2.7478,0.141)]},
                         (2,-1,1):{'Resonances':[(-0.5489,0.0)]},
                         (3,-1,1):{'Resonances':[(0.5588,0.0157),(3.3378,0.296)]},
                         (4,-1,1):{'Resonances':[(2.4578,0.150)]}}}
    return dic_spectra[nuclei]

def paper_spectra(key):
    dic_spectra={ 'Lazauskas_4H_MT18':{(0,-1,1):{'Resonances':[(1.08,4.06,0.03)]},
                         (1, -1,1):{'Resonances':[(1.08,4.06,0.03),(0.88,4.4,0.1)]},
                         (2,-1,1):{'Resonances':[(1.08,4.06,0.05)]},
                         'Thresholds':[[-8.48,0],[-2.224,0],[0,0]]
                         ,"threshold_labels":[r'$\mathrm{^3H\text{-}n}$',r'$\mathrm{D\text{-}nn}$',r'$\mathrm{3np}$']
                        },

                'Lazauskas_4H_N3LO':{(0,-1,1):{'Resonances':[(0.78,7.7,0.15)]},
                         (1, -1,1):{'Resonances':[(0.9,3.8,0.03),(0.2,4.6,0.2)]},
                         (2,-1,1):{'Resonances':[(1.15,3.97,0.05)]},
                         'Thresholds':[[-8.48,0],[-2.224,0],[0,0]]
                         ,"threshold_labels":[r'$\mathrm{^3H\text{-}n}$',r'$\mathrm{D\text{-}nn}$',r'$\mathrm{3np}$']
                        },
                }
    if(not (key in dic_spectra)):
        sys.exit("comparison spectra not found")
    return dic_spectra[key]
def save_2_file(dic_in,file_in):
    the_file = open(file_in+".pkl", "wb")
    pickle.dump(dic_in, the_file)
    the_file.close()

def read_from_file(file_in):
    the_file = open(file_in+".pkl", "rb")
    return  pickle.load(the_file)
##############################
def expected_theta(resonances,nuclei):
    thresholds_exp=np.array(resonances['Thresholds'])[:,0]
    threshold_list=[ el -thresholds_exp[0] for el in thresholds_exp]
    new_dict={}
    for key,value in resonances.items():
        if(isinstance(key,str)):
            continue
        res_list=value['Resonances']
        new_res_list=[]
        print("*********\n",key)
        for res in res_list:
            pos=np.array([res[0]-el for el in threshold_list if el<res[0]])
            theta=np.arctan(res[1]/pos/2)/2
            new_res_list.append((res[0],res[1],theta))
            print(theta)
            new_dict[key]=new_res_list
    return new_dict
def read_resonance_file(file_name):
    import pickle
    with open(file_name+'.pkl', 'rb') as f:
        loaded_data = pickle.load(f)
    print(loaded_data)
    return loaded_data
def plot_expected_thetas(expected_thetas,nucleus):
    colors=['grey', 'purple', 'blue', 'green', 'orange', 'red','saddlebrown',
             'pink','cyan','olive']
    plt.rcParams.update({'font.size': 12})
    figsize=fig_size(2,1,0.7)
    fig,ax=plt.subplots(1,1,figsize=figsize)
    k=0
    all_x=[]
    all_y=[]
    amplitudes=np.linspace(0.5,1,len(expected_thetas.items()))
    for i,(key,value) in enumerate(expected_thetas.items()):
        x=[]
        y=[]
        for resonance in value:
            theta=resonance[2][-1]
            x_value=amplitudes[i]*np.cos(theta)
            y_value=amplitudes[i]*np.sin(theta)
            x.append(x_value)
            y.append(y_value)
            all_x.append(x_value)
            all_y.append(y_value)
            if any([abs(x[-1]-all_x[j])+3*abs(y[-1]-all_y[j])<0.15 for j in range(len(all_x)-1)]):
                continue
            plt.annotate('${:.2f}'.format(resonance[2][-1])+r'$',
                         (x[-1], y[-1]), textcoords="offset points",
                         xytext=(-16,-2), ha='center',fontsize=10)
        pi_str=lambda x: '+' if x>0 else '-'
        label='${:2d}'.format(key[0])+'^{'+pi_str(key[1])+'};'+'{:2d}$'.format(key[2])
        ax.scatter(x,y,label=label,color=colors[i])
    ax.plot([0,1],np.tan(0.3)*np.array([0,1]),ls='--',color='k',label='$\mathrm{SRG\ limit}$')
    plt.fill_between([0,1], [0,np.tan(0.3)], color='gray',
                     alpha=0.2,hatch='\\')
    legend= ax.legend(loc='best',prop={'size': 10},frameon=False)
    ax.set_xlim(0,1)
    ax.set_ylim(0,1)
    plt.tight_layout()
    plt.show()
    outName='spectra/Expected_thetas_'+nucleus
    fig.savefig(outName+'.pdf', bbox_inches='tight', transparent=True,dpi=2000)
    fig.savefig(outName+'.png', bbox_inches='tight', transparent=True,dpi=2000)
    fig.savefig(outName+'.eps', bbox_inches='tight')
    fig.savefig(outName+'.svg', bbox_inches='tight', transparent=True,dpi=2000)
    return
def quick(nuclei1,nuclei2):
    expt_4he=expt_spectra(nuclei1)
    expected_thetas=expected_theta(expt_4he,nuclei1)
    plot_expected_thetas(expected_thetas,nuclei1)
    #plot_spectra([pred_4he,expt_4he])
def energy_levels(file_names,nucleus,bar_style='width',plot_thresholds=False,comparison='exp'):
    split=False
    if(comparison=='exp'):
        comparison_dict=expt_spectra(nucleus)
    else:
        comparison_dict=paper_spectra(comparison)
    if (isinstance(file_names,str)):
        file_names=[file_names]
    pred_dicts=[read_resonance_file(file_name) for file_name in file_names]
    pred_dicts=pred_dicts[::-1]
    # add error values if absent
    for dict_i in pred_dicts:
        for key in dict_i:
            if(isinstance(key,str)):
                continue
            for i in range(len(dict_i[key]['Resonances'])):
                if(len(dict_i[key]['Resonances'][i])==2):
                    dict_i[key]['Resonances'][i]=(dict_i[key]['Resonances'][i][0],dict_i[key]['Resonances'][i][1],0)
    for key in comparison_dict:
        if(isinstance(key,str)):
            continue
        for i in range(len(comparison_dict[key]['Resonances'])):
            if(len(comparison_dict[key]['Resonances'][i])==2):
                comparison_dict[key]['Resonances'][i]=(comparison_dict[key]['Resonances'][i][0],comparison_dict[key]['Resonances'][i][1],0)
    temp_dict={}
    #for filter on R-matrix
    #for key in comparison_dict:
    #    if(isinstance(key,str)):
    #        temp_dict[key]=comparison_dict[key]
    #        continue
    #    if(key[2]==1):
    #        temp_dict[key]=comparison_dict[key]
    #comparison_dict=temp_dict
    #
    #if(plot_thresholds):
    #    thresholds_in,thresholds_errors_in=read_thresholds()
    #    Thresholds=np.column_stack((thresholds_in,thresholds_errors_in))
    #else:
    #    Thresholds=None
    if(split):
        temp_dict={}
        temp_dict2={}
        second_set=[3.9,4,3.74]
        for key in pred_dict:
            temp_list=[]
            temp_list2=[]
            for i in range(len(pred_dict[key]['Resonances'])):
                #if( pred_dict[key]['Resonances'][i][1]<0.6):
                if( any(abs(pred_dict[key]['Resonances'][i][0]-position)<0.1 for position in second_set)):
                    temp_list2.append(pred_dict[key]['Resonances'][i])
                else:
                    temp_list.append(pred_dict[key]['Resonances'][i])
            if(temp_list):
                temp_dict[key]={'Resonances':temp_list.copy()}
            if(temp_list2):
                temp_dict2[key]={'Resonances':temp_list2.copy()}
            print(len(temp_list),len(temp_list2),len(pred_dict[key]['Resonances']))
        print(temp_dict,temp_dict2)
        plot_spectra([comparison_dict,temp_dict,temp_dict2],style=bar_style,plot_thresholds=plot_thresholds)
    else:
        plot_spectra([comparison_dict,*pred_dicts],style=bar_style,plot_thresholds=plot_thresholds)

##############################
def read_thresholds(parameter_filter=None,n=5):
    columns=['E_r','N_max','hw','srg','lambda','NNN','3N_induced','theta','Nf','srg_f','Er_total']

    if(os.path.exists('thresholds.dat')):
        threshold_data=pd.read_csv('thresholds.dat',sep=',')
    else:
        sys.exit("NO threshold file")
    subdata=threshold_data
    if(not(parameter_filter is None)):
        mask=check_params(params,threshold_data)
        subdata=threshold_data[mask]
    num_rows = subdata.shape[0]
    if(num_rows>n):
        print("parameters are :",columns[1:-1])
        answer=input("Enter the paramters of the threshold that you want to use")
        answer=answer.split()
        input_error= [name for name in answer if name not in columns]
        if(len(input_error)>0):
            sys.exit("error in entered parameters")
        values=input("Enter values:")
        values=values.split()
        values=[float(value) for value in values]
        params=dict(zip(answer,values))
        mask=check_params(params,threshold_data)
        subdata=threshold_data[mask]
    print("current thresholds:\n",subdata[['E_r','Er_total']])
    return subdata['E_r'].tolist(),subdata['Er_total'].tolist()
def check_params(params,datai=None):
    if(datai is None):
        datai=data
    mask = pd.Series(True, index=datai.index)
    for column,value in params.items():
        if(column=='lambda'):
            mask=mask & ((datai[column]==value)| (datai['lambda']==0))
        else:
            mask=mask & (datai[column]==value)
    return mask
