#! /usr/bin/python3
import sys
import os, time, errno, shutil
import re
import glob
import numpy as np
import subprocess,functools,socket
import ipyparallel as ipp
import asyncio
from tqdm import tqdm
from scipy.optimize import Bounds,minimize,minimize_scalar
from datetime import datetime
IJCLab_server=[r'ld-theo[0-9]*\..*$',r'ls-theo[0-9]*\..*$',r'ls-phynet[0-9]*\..*$',r'theorie[0-9]*\..*$']
i=1


import psutil
def monitor_cpu_usage(threshold=0.3, check_interval=10,maxtime=3600*24, *args, **kwargs):
    """
    Monitors CPU usage and when CPU usage is below a threshold.
    """
    i=0
    start_time=time.time()
    while True:
        cpu_usage = psutil.cpu_percent(interval=1)
        sys.stdout.write(f"\rCurrent CPU usage: {cpu_usage}%")
        sys.stdout.flush()  # Ensures the output is printed immediately
        if cpu_usage < threshold:
            print(f"\nCPU usage is below {threshold}%. Executing the code...")
            break  # break the loop once the usage is below the threshold

        time.sleep(check_interval)  # Check again after the specified interval
        if (time.time()-start_time>maxtime):
            print("\nmaxtime exceeded. Executing the code...")
            break

def time_format(t):
    seconds=t%60
    t=t//60
    minutes=t%60
    t=t//60
    hours=t%24
    t=t//24
    days=t
    s=""
    if(days>0):
        s=s+" {} days ".format(days)
    if( hours>0):
        s=s+" {} hours  ".format(hours)
    s=s+" {} minutes".format(minutes+1)
    return s
def modify_pless_input_file(**pless_params):
    import fileinput

    filename='pless-params.in'
    for line in fileinput.input(files=filename,inplace=True):
        if 'H0' in line and 'H0'in pless_params:
            print(re.sub(r"(^H0 +)[-.0-9]+e[-+0-9]*( +[.0-9]+)( +[.0-9]+)(.*)",r"\g<1>"+'{:.6e}'.format(pless_params['H0'])+'{:5.1f}'.format(pless_params['power'])+'{:8.1f}'.format(pless_params['cutoff'])+"\g<4>",line) if all(keys in pless_params for keys in ('cutoff','power')) else re.sub(r"(^H0 +)[-.0-9]+e[-+0-9]*( +[.0-9]+)( +[.0-9]+.*)",r"\g<1>"+'{:.6e}'.format(pless_params['H0'])+"\g<2>\g<3>",line),end='')
        elif '3S1' in line and '3S1'in pless_params:
            print(re.sub(r"(^3S1 +)[-.0-9]+e[-+0-9]*( +[.0-9]+)( +[.0-9]+)(.*)",r"\g<1>"+'{:.6e}'.format(pless_params['3S1'])+'{:5.1f}'.format(pless_params['power'])+'{:8.1f}'.format(pless_params['cutoff'])+"\g<4>",line) if all(keys in pless_params for keys in ('cutoff','power')) else re.sub(r"(^3S1 +)[-.0-9]+e[-+0-9]*( +[.0-9]+)( +[.0-9]+.*)",r"\g<1>"+'{:.6e}'.format(pless_params['3S1'])+"\g<2>\g<3>",line),end='')
        elif '1S0' in line and '1S0'in pless_params:
            print(re.sub(r"(^1S0 +)[-.0-9]+e[-+0-9]*( +[.0-9]+)( +[.0-9]+)(.*)",r"\g<1>"+'{:.6e}'.format(pless_params['1S0'])+'{:5.1f}'.format(pless_params['power'])+'{:8.1f}'.format(pless_params['cutoff'])+"\g<4>",line) if all(keys in pless_params for keys in ('cutoff','power')) else re.sub(r"(^1S0 +)[-.0-9]+e[-+0-9]*( +[.0-9]+)( +[.0-9]+.*)",r"\g<1>"+'{:.6e}'.format(pless_params['1S0'])+"\g<2>\g<3>",line),end='')
        elif 'Regulator' in line:
            print(re.sub(r"(^Regulator +)([-a-z]*)(.*)",r"\g<1>"+pless_params['regulator']+"\g<3>",line) if 'regulator' in pless_params.keys() else re.sub(r"(^Regulator +)([a-z]*)(.*)",r"\g<1>"+"\g<2>"+"\g<3>",line),end='')
        else:
            print(line,end='')
    fileinput.close()
    return

def delete_val_pless_input_file(**pless_params_rm):
    import fileinput

    filename='pless-params.in'
    for line in fileinput.input(files=filename,inplace=True):
        if 'H0' in line and 'H0'in pless_params_rm:
            continue
        elif '3S1' in line and '3S1'in pless_params_rm:
            continue
        elif '1S0' in line and '1S0'in pless_params_rm:
            continue
        elif 'Regulator' in line and 'regulator' in pless_params_rm:
            continue
        else:
            print(line,end='')
    fileinput.close()
    return

def modify_manychir_input_file(hw_in):
    import fileinput

    filename='manychir.in'
    for line in fileinput.input(files=filename,inplace=True):
        if 'hbar*Omega' in line:
            new_line=re.sub(r"^[.0-9]+d[-+0-9]*(.*$)",r'{:.2E}'.format(hw_in)+"\g<1>",line)
            new_line=new_line.replace('E','d')
            print(new_line,end='')
        else:
            print(line,end='')
    fileinput.close()
    return

def run_manyeff_light():

    #file_pattern=['*.dat','*.bin','*.tmp','v3b*']
    file_pattern=['v3b*']
    files=[]
    for pattern in file_pattern:
        files.extend(glob.glob(pattern))
    for file_to_remove in files:
        os.remove(file_to_remove)
    with open(os.devnull, 'w') as fp:
        try:
            subprocess.call(['./manyeffv3b.exe'],stdout=fp)
        except:
            raise
    return

def get_gs_energy():

    pattern_2_match=re.compile(r" *Ground-state energy\n *([-.0-9]*)\n")
    many_out =open('many.out','r').readlines()
    gs_energy=re.findall(pattern_2_match,"".join(many_out))[0]

    return float(gs_energy)

def make_run(H0_in,NN_params=None,NN_cutoffs=None):

    H0_in=H0_in*1e-4
    print('Calculating for H_0 value of ','{:.6e}'.format(H0_in))
    triton_be=-8.481

    if NN_params:
        NN_params.update({'H0':H0_in})
    else:
        NN_params={'H0':H0_in}
    if NN_cutoffs: NN_params.update(NN_cutoffs)
    modify_pless_input_file(**NN_params)
    run_manyeff_light()
    new_gs=get_gs_energy()

    print('Triton binding energy is ','{:7.3f}'.format(new_gs))

    return np.abs(triton_be-new_gs)

def fit_gs_energy(pless_params_forfit=None,pless_cutoffs_forfit=None,bounds = (-1.10000, 5.500000)):

    res = minimize_scalar(make_run, method='bounded', bounds=bounds, args=(pless_params_forfit,pless_cutoffs_forfit))

    print(res.x)

    return res.x

def fit_hw_dependance(hw_min,hw_max):
    step_hw =2
    hw      =np.linspace(hw_min ,hw_max ,max([1 ,int((hw_max -hw_min )/step_hw   )+1]),endpoint=True,dtype=float)
    results =np.array([])

    for value in hw:
        modify_manychir_input_file(value)
        results=np.append(results,fit_gs_energy())

    print(hw)
    print(results)

    return hw,results

def get_pless_params(cutoff_in,cutoff_exp,regulator='non-local',physics=None):
    import pandas as pd
    def reform_keys(dic_in,addin):
        for keys,vals in dic_in.items():
            yield (keys,addin),vals

    NN_cutoffs={'power':cutoff_exp,'cutoff':cutoff_in}
    result    =NN_cutoffs
    if regulator == 'local':
        df       = pd.read_excel("LEC-Lorenzo.xlsx",sheet_name="LECS",skiprows=33,usecols=["Cut-Off (fm^-1)", "C0"])
        df       = df[df.C0.apply(lambda x: isinstance(x, (float, int)))]
        df       = df.rename({'Cut-Off (fm^-1)': 'cutoff','C0': '1S0'}, axis='columns')
        df.cutoff= df.cutoff.apply(lambda x: round(x, 1))
        df["3S1"]= df["1S0"]
        df.loc[[not val for val in df.cutoff.isnull()], "regulator"] = "local"

        blanks   = df[df.isnull().all(1)]
        df_list  = np.split(df,blanks.index)
        dic_input= {0:'Unitarity',1:'Neutron matter',2:'LQCD'}

        i=0
        for dfs in df_list:
            dfs= dfs.dropna(how="all",axis=0)
            if not dfs.empty:
                dic_input[dic_input[i]]=dfs
                i+=1
        dic_out   = dic_input["Unitarity"].set_index("cutoff").T.to_dict("dict")
        dic_param = dict(reform_keys(dic_out, 1.0))
    else:
        dic_param={
            (200.0,1.0): {'3S1':-1.44476e-05,'1S0':-7.69537e-06,'H0':-8.3418e-05},
            (400.0,1.0): {'3S1':-5.59689e-06,'1S0':-4.03838e-06,'H0': 2.3342e-04},
            (600.0,1.0): {'3S1':-3.41188e-06,'1S0':-2.73747e-06,'H0': 4.4531e-04},
            (400.0,2.0): {'3S1':-4.42702e-06,'1S0':-3.34979e-06,'H0': 2.5154e-04},
            (600.0,2.0): {'3S1':-2.72698e-06,'1S0':-2.26422e-06,'H0': 3.7016e-04}}

    result.update(dic_param[(cutoff_in,cutoff_exp)])
    return result

def get_pless_params_forfit(cutoff_in,cutoff_exp):
    dic_param={
        (200.0,1.0): {'3S1':-1.44476e-05,'1S0':-7.69537e-06,'H0':None},
        (400.0,1.0): {'3S1':-5.59689e-06,'1S0':-4.03838e-06,'H0':None},
        (600.0,1.0): {'3S1':-3.41188e-06,'1S0':-2.73747e-06,'H0':None},
        (400.0,2.0): {'3S1':-4.42702e-06,'1S0':-3.34979e-06,'H0':None},
        (600.0,2.0): {'3S1':-2.72698e-06,'1S0':-2.26422e-06,'H0':None}}
    return dic_param[(cutoff_in,cutoff_exp)],{'cutoff':cutoff_in,'power':cutoff_exp}

def merge_pless_params(new_params,pless_params_in):
    for keys,items in pless_params_in.items():
        if type(new_params[keys]) is list:
            new_params.setdefault(keys,[]).append(items)
        else:
            new_params[keys]=[new_params[keys],items]
    return new_params

def run_manyeff(nuclei,NN_interaction=['n3lo'],NNN_interaction={'regulator':'lnl','cutoff':[650,500],'c_D':[0.7],'c_E':[-0.06],'E8':[0.0]},run_params={'srg_evolution':True,'NNN_srg':False,'NNNN_srg':False,'lambda_range':[2.0],'srg_f':1.0,'N_m_range':[14],'hw_range':[20],'theta_range':[0.0],'Nf':0},pless_params=None,output_stream=0,J=None,T=None,testing=False):
    import fileinput
    import itertools as it

    NUCLEI_prop={'3He'     :{'protons':2,'neutrons':1,'J':1,'T':1,'T12_max_3NF' :1,'T3_max_3NF':1,'J12_max_3NF':1},
                 '3H'      :{'protons':1,'neutrons':2,'J':1,'T':1,'T12_max_3NF' :1,'T3_max_3NF':1,'J12_max_3NF':1},
                 '4He'     :{'protons':2,'neutrons':2,'J':0,'T':0,'T12_max_3NF' :0,'T3_max_3NF':1,'J12_max_3NF':7},
                 '3A'      :{'protons':2,'neutrons':2,'J':1,'T':1,'T12_max_3NF':1,'T3_max_3NF':3,'J12_max_3NF':13},
                 '4H':{'protons':1,'neutrons':3,'J':0,'T':1,'T12_max_3NF':2,'T3_max_3NF':3,'J12_max_3NF':7},
                 '4Li':{'protons':3,'neutrons':1,'J':0,'T':1,'T12_max_3NF':2,'T3_max_3NF':3,'J12_max_3NF':7},
                 '4N':{'protons':0,'neutrons':4,'J':0,'T':2,'T12_max_3NF':4,'T3_max_3NF':3,'J12_max_3NF':7,'T3_min_3NF':3},
                 '3N':{'protons':0,'neutrons':3,'J':1,'T':3,'T12_max_3NF':1,'T3_max_3NF':3,'J12_max_3NF':7,'T3_min_3NF':3},
                }
    if(J!=None):
        NUCLEI_prop[nuclei]['J']=J
    if(T!=None):
        NUCLEI_prop[nuclei]['T']=T
        NUCLEI_prop[nuclei]['T12_max_3NF']=2*T
        NUCLEI_prop[nuclei]['T3_max_3NF']=3 if T>0 else 0
    if('T3_min_3NF' in NUCLEI_prop[nuclei]):
        NUCLEI_prop['3H']['T']=NUCLEI_prop[nuclei]['T3_min_3NF']
    start_time=time.time()
    def modify_manychir_input_file(nuclei_in,N_max,h_w,theta,NNN_calc,threeff=False,NN_srg=False,NNN_srg=False,NNNN_srg=False,J12_max_3NF=NUCLEI_prop[nuclei]['J12_max_3NF'],Nmax_3N=None,readtmp=False):
        import fileinput
        global f,line
        if not Nmax_3N:
            Nmax_3N=[N_m_range[-1]]*int((J12_max_3NF+1)/2)
        else:
            Nmax_3N_tmp=Nmax_3N
            ln=min(len(Nmax_3N),int((J12_max_3NF+1)/2))
            Nmax_3N=[N_m_range[-1]]*int((J12_max_3NF+1)/2)
            Nmax_3N[:ln]=Nmax_3N_tmp[:ln]
        if(threeff and nuclei_in[0]=='3'):
            print('Nmax for each J: ',Nmax_3N)
        def J1max3eff_Nmax(a):
            global f,line
            new_line=re.sub(r"^[0-9]*(.*$)",r'{:<2d}'.format(J12_max_3NF)+"\g<1>",line)
            new_line=new_line.rstrip()
            new_line=new_line+"\n"
            line=next(f)
            stored_line=''
            for i,J in enumerate(range(1,J12_max_3NF+1,2)):
                if 'Nmax' in line:
                    new_line=new_line+re.sub(r"^[0-9]*(.*$)",r'{:<3d}'.format(Nmax_3N[i])+"\g<1>",line)
                    line=next(f)
                else:
                    new_line=new_line+'{0:<3d}            ! Nmax for J1={1:<d}/2'.format(Nmax_3N[i],J)+'{}'.format("\n")
                    stored_line+=line
                    line=''
            while 'Nmax' in line:
                line=next(f)
            print(new_line+stored_line,end='')
            raise ValueError

        INPUT_ctrl={'protons' :                     ( lambda line : re.sub(r"^[0-9]*(.*$)",r'{:<2d}'.format(NUCLEI_prop[nuclei_in]['protons'])+"\g<1>",line)),
                    'neutrons' :                    ( lambda line : re.sub(r"^[0-9]*(.*$)",r'{:<2d}'.format(NUCLEI_prop[nuclei_in]['neutrons'])+"\g<1>",line)),
                    'hbar*Omega' :                  ( lambda line : re.sub(r"^[.0-9]+[Ed]?[-+0-9]*(.*$)",r'{:.2E}'.format(h_w)+"\g<1>",line)),
                    'Nmax' :                        ( lambda line : re.sub(r"^[0-9]*(.*$)",r'{:<3d}'.format(N_max)+"\g<1>",line)),
                    'J, 2*J for odd A' :            ( lambda line : re.sub(r"^[0-9]*(.*$)",r'{:<2d}'.format(NUCLEI_prop[nuclei_in]['J'])+"\g<1>",line)),
                    'T, 2*T for odd A' :            ( lambda line : re.sub(r"^[0-9]*(.*$)",r'{:<2d}'.format(NUCLEI_prop[nuclei_in]['T'])+"\g<1>",line)),
                    'threeff' :                     ( lambda line :
                                                     re.sub(r"^[TF]*(.*$)",r'{}'.format(not
                                                                                        (NN_only)
                                                                                        or
                                                                                        (False
                                                                                         and NNN_srg)
                                                                                       or
                                                                                       threeff)[0]+"\g<1>",line)),
                    'threeffcal' :                  ( lambda line : re.sub(r"^[TF]*(.*$)",r'{}'.format(run_params.get('threeffcal',False))[0]+"\g<1>",line)),
                    'chirp' :                       ( lambda line : re.sub(r"^[TF](.*$)",r'{}'.format(LEC_NNN[NN]['chirp'])[0]+"\g<1>",line)),
                    'cutnum' :                      ( lambda line : re.sub(r"^[0-9]*(.*$)",r'{:<3d}'.format(icut)+"\g<1>",line)),
                    'iprot3eff' :                   ( lambda line :
                                                     re.sub(r"^[0-9]*(.*$)",r'{:<2d}'.format(NUCLEI_prop[nuclei]['protons'])+"\g<1>",line)),
                    'neut3eff' :                    ( lambda line : re.sub(r"^[0-9]*(.*$)",r'{:<2d}'.format(NUCLEI_prop[nuclei]['neutrons'])+"\g<1>",line)),
                    'itot23eff' :                   ( lambda line : re.sub(r"^[0-9]*(.*$)",r'{:<2d}'.format(NUCLEI_prop[nuclei]['T12_max_3NF'])+"\g<1>",line)),
                    'J1max3eff' :                   J1max3eff_Nmax,
                    'mscheme' :                     ( lambda line : re.sub(r"^[TF]*(.*$)",r'{}'.format(run_params.get('mscheme',False))[0]+"\g<1>",line)),
                    'real3b' :                      ( lambda line : re.sub(r"^[TF]*(.*$)",r'{}'.format(NNN_calc)[0]+"\g<1>",line)),
                    'cutLambda' :                   ( lambda line : re.sub(r"^[.0-9]+[Ed]?[-+0-9]*(.*$)",r'{:.2E}'.format(NNN_cutoff_Lambda)+"\g<1>",line)),
                    'aTM' :                         ( lambda line : re.sub(r"^[-.0-9]+[Ed]?[-+0-9]*(.*$)",r'{:.2E}'.format(LEC_NNN[NN]['c_1'])+"\g<1>",line)),
                    'bTM' :                         ( lambda line : re.sub(r"^[-.0-9]+[Ed]?[-+0-9]*(.*$)",r'{:.2E}'.format(LEC_NNN[NN]['c_3'])+"\g<1>",line)),
                    'dTM' :                         ( lambda line : re.sub(r"^[-.0-9]+[Ed]?[-+0-9]*(.*$)",r'{:.2E}'.format(LEC_NNN[NN]['c_4'])+"\g<1>",line)),
                    'cTM' :                         ( lambda line : re.sub(r"^[-.0-9]+[Ed]?[-+0-9]*(.*$)",r'{:.2E}'.format(c_TM)+"\g<1>",line)),
                    'cON' :                         ( lambda line : re.sub(r"^[-.0-9]+[Ed]?[-+0-9]*(.*$)",r'{:.2E}'.format(c_ON)+"\g<1>",line)),
                    'nonloc cutLambdanonloc' :      ( lambda line : re.sub(r"^[TF]*( *)[.0-9]+[Ed]?[-+0-9]*(.*$)",r'{}'.format(nonlocal_regulator)[0]+"\g<1>"+r'{:.2E}'.format(NNN_cutoff_nonlocal_Lambda)+"\g<2>",line)),
                    'c_D' :                         ( lambda line : re.sub(r"^[-.0-9]+[Ed]?[-+0-9]*(.*$)",r'{:.2E}'.format(c_D)+"\g<1>",line)),
                    'c_E' :                         ( lambda line : re.sub(r"^[-.0-9]+[Ed]?[-+0-9]*(.*$)",r'{:.2E}'.format(c_E)+"\g<1>",line)),
                    'c_LS   x_78' :                 ( lambda line : re.sub(r"^[-.0-9]+[Ed]?[-+0-9]*( *)[-.0-9]+[Ed]?[-+0-9]*(.*$)",r'{:.2E}'.format(c_LS)+"\g<1>"+r'{:.2E}'.format(x_78)+"\g<2>",line)),
                    'srg' :                         ( lambda line : re.sub(r"^[TF]*(.*$)",r'{}'.format(srg_evolution)[0]+"\g<1>",line)),
                    'NN_SRG_for3b_extract' :        ( lambda line : re.sub(r"^[TF]*(.*$)",r'{}'.format(  NN_srg)[0]+"\g<1>",line)),
                    'NNN_SRG' :                     ( lambda line : re.sub(r"^[TF]*(.*$)",r'{}'.format( NNN_srg)[0]+"\g<1>",line)),
                    'NNNN_SRG' :                    ( lambda line : re.sub(r"^[TF]*(.*$)",r'{}'.format(NNNN_srg)[0]+"\g<1>",line)),
                    'it1max3eff' :                  ( lambda line : re.sub(r"^[0-9]*(.*$)",r'{:<2d}'.format(NUCLEI_prop[nuclei]['T3_max_3NF'])+"\g<1>",line)),
                    'lambda' :                      ( lambda line : re.sub(r"^[-.0-9]+[Ed]?[-+0-9]*(.*$)",r'{:.2E}'.format(Lambda)+"\g<1>",line)),
                    'srg_f' :                       ( lambda line :  re.sub(r"^[-.0-9]+[Ed]?[-+0-9]*(.*$)",r'{:.2E}'.format(srg_f)+"\g<1>",line)),
                    'save_NNN_terms' :              ( lambda line : re.sub(r"^[TF]*(.*$)",r'{}'.format(output_NNN_interaction)[0]+"\g<1>",line)),
                    'obd' :                         ( lambda line : re.sub(r"^[TF]*(.*$)",r'{}'.format(observable_obd_calculation)[0]+"\g<1>",line)),
                    'freq_conv' :                   ( lambda line : re.sub(r"^[TF]*(.*$)",r'{}'.format(frequency_conversion)[0]+"\g<1>",line)),
                    'n_Omega' :                     ( lambda line : re.sub(r"^[0-9]*(.*$)",r'{:<2d}'.format(n_omega)+"\g<1>",line)),
                    'delta_omega' :                 ( lambda line : re.sub(r"^[-.0-9]+[Ed]?[-+0-9]*(.*$)",r'{:.2E}'.format(delta_omega)+"\g<1>",line)),
                    'theta':                        ( lambda line : re.sub(r"^[-.0-9]+[Ed]?[-+0-9]*(.*$)",r'{:.5f}'.format(theta)+"\g<1>",line)),
                    'readbas' :                   ( lambda line : re.sub(r"^[TF]*(.*$)",r'{}'.format(readtmp)[0]+"\g<1>",line)),
                    'MTV' :                      ( lambda line :re.sub(r"^[TF]*(.*$)",r'{}'.format(NN=='mtv')[0]+"\g<1>",line)),
                    'Nf for fitting' :                        ( lambda line : re.sub(r"^[0-9]*(.*$)",r'{:<3d}'.format(Nf)+"\g<1>",line)),
                    'v2b_method choose 1 for direct, 2 for fitting' :
                    ( lambda line: re.sub(r"^[0-9]*(.*$)",r'{:<3d}'.format(2 if Nf>0 else 4 if Nf<0 else 1)+"\g<1>",line)),
                    'N2LO':                         ( lambda line : re.sub(r"^[TF]*(.*$)",r'{}'.format(LEC_NNN[NN]['chirp'])[0]+"\g<1>",line))}

        ctrl_pattern  =re.compile(r"! +(.*)")
        f=fileinput.input(files=('manychir_v1.in'),inplace=True)
        line=next(f)
        key_errors=[]
        while True:
            try:
                ctrl_key=re.search(ctrl_pattern,line)
                try:
                    new_line=INPUT_ctrl[ctrl_key.group(1).rstrip()](line)
                except KeyError:
                    new_line=line
                    key_errors.append(ctrl_key.group(1).rstrip())
                except Exception as ex:
                    new_line=line
                    key_errors.append(ctrl_key.group(1).rstrip())
                print(new_line,end='')
                line=next(f)
            except StopIteration:
                break
            except ValueError:
                continue
        f.close()

        with open('key_errors.log','a') as f:
            f.write(' '.join(key_errors))

    def rem(file,afix=''):
        if os.path.isfile(file+afix):
            os.remove(file+afix)

    def ren(file,newname):
       if os.path.isfile(file):
          os.rename(file,newname)

    def rem_auto_read_files():
        files=['Ham_NN_A3_from3eff.bin','Trel_A3_from3eff.bin','Ham_NN_NNN_A3.bin','Trel_A3_from3bre.bin','Ham_NN3N_A3_from3eff.bin','Trel_NN3N_A3_from3eff.bin']
        for file in files:
            rem(file)
#        files=glob.glob("v3b*")
#        for file in files:
#            rem(file)

    def symlink_force(target, link_name):
        try:
            temp_link = link_name+'.new'
            rem(temp_link)
            os.symlink(target, temp_link)
            os.rename(temp_link, link_name)
        except OSError as e:
            raise


    def create_temporary_runspace(dependency_files=None):
        from datetime import datetime
        current_time = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
        folder="run_"+current_time
        os.mkdir(folder)
        parent_dir=os.getcwd()
        files_in_directory=os.listdir(parent_dir)
        if(dependency_files):
            Ham_3b_file,Ham_2b_file,Trel_A3_file=dependency_files
            if (not NN_only or NNN_srg) :symlink_force(os.path.join(parent_dir,Ham_3b_file),os.path.join(folder,'Ham_NN_NNN.bin'))
            symlink_force(os.path.join(parent_dir,Ham_2b_file),os.path.join(parent_dir,folder,'Ham_NN.bin'))
            symlink_force(os.path.join(parent_dir,Ham_2b_file),os.path.join(parent_dir,folder,'Ham_NN_subtr.bin'))
            symlink_force(os.path.join(parent_dir,Trel_A3_file),os.path.join(parent_dir,folder,'Trel.bin'))
            print(dependency_files)
        for file in files_in_directory:
            if('.bin' in file and ('Ham'in file or 'Trel' in file) and not "LaguerrePoly" in file):
                continue
            if('tmp_' in file):
                source=os.path.join(parent_dir,file)
                destination=os.path.join(parent_dir,folder,file)
                symlink_force(source,destination)
                continue
            source=os.path.join(parent_dir,file)
            destination=os.path.join(parent_dir,folder,file)
            if(os.path.isdir(file)):
                if(not file.startswith("run_")):
                    shutil.copytree(source,destination)
            else:
                shutil.copy2(source,destination)
        os.chdir(os.path.join(parent_dir,folder))


    def clean_temporary_runspace(output_file={'many.out':'many.out_{:f}'.format(time.time())},tmp_ext=''):
        curr_dir=os.getcwd()
        parent_dir=os.path.dirname(curr_dir)
#       copy output to parent directory
        shutil.copy('manychir_v1.in',os.path.join(parent_dir,'manychir_run_'+str(i)+'.in'))
        for output,new_name in output_file.items():
            if(os.path.exists(new_name)):
                shutil.copy(new_name,os.path.join(parent_dir,new_name))

        files=glob.glob("*.tmp")
        for file in files:
            ren(file,file+tmp_ext)
        for file in files:
            if(os.path.exists(file+tmp_ext) and tmp_ext!= ''):
                shutil.copy(file+tmp_ext,os.path.join(parent_dir,file+tmp_ext))
        files=glob.glob("LaguerrePoly*.bin")
        for file in files:
            shutil.copy(file,os.path.join(parent_dir,file))
        #go back to the original directory
        os.chdir(parent_dir)
        #rename the run-space folder and remove it if necessary
        os.rename(curr_dir,curr_dir+"_finished")
        temp_run_folders=[rfolder for rfolder in os.listdir() if
                          os.path.isdir(rfolder) and rfolder.startswith("run_")
                         and "finished" in rfolder]
        temp_run_folders.sort()
        if(len(temp_run_folders)>3):
            shutil.rmtree(temp_run_folders[0])



    def run_code(output_file={'many.out':'many.out_{:f}'.format(time.time())}):
        global i
#       go to a temporary working directory
#
        rem_auto_read_files()
        symlink_force('manychir_v1.in','manychir.in')
        with open(os.devnull, 'w') as fp, open("manyeff.log","a") as f:
            try:
                if(output_stream==1):
                    subprocess.call(['./'+exefile])
                elif (output_stream==2):
                    subprocess.call(['./'+exefile],stdout=f)
                else:
                    subprocess.call(['./'+exefile],stdout=fp)

            except OSError:
                sys.exit("Error while running the fortran program")
        rem('manychir','.in')
        shutil.copy('manychir_v1.in','manychir_run_'+str(i)+'.in' )
        i+=1
        for output,new_name in output_file.items():
            if(os.path.exists(output)==False):
                #sys.exit("Error: file:{} does not exist at the end of calculation".format(output))
                print("Error: file:{} does not exist at the end of calculation".format(output))
            ren(output,new_name)
#

    def empty(*args):
        pass

    def make_list(unknown_type):
        return [unknown_type] if type(unknown_type) is not list else unknown_type

    def run_4He_dependance():
        print("running dependence He")
        #extension_out='{0}{1}{2}_{3:d}.{4:d}{5}'.format(nuclear_pot,'-srg' if srg_evolution else '', '{:3.1f}'.format(Lambda) if srg_evolution else '',N_m_range[-1],int(hw),'_theta_{:0.2f}'.format(theta) if theta != 0 else '')
        extension_out='{0}{1}{2}{3}_srgf{12}_{4}-T12_{5}_N{6:d}_hw{7:d}_J{8}_T{9}_theta{10:0.3f}{11}'.format(nuclear_pot,'_srg' if srg_evolution else '', '{:3.1f}'.format(Lambda) if srg_evolution else '','-4Nsrg' if srg_4N_evolution else '',nuclei,NUCLEI_prop['3H']['T12_max_3NF'],N_m,int(hw),NUCLEI_prop[nuclei]['J'],NUCLEI_prop[nuclei]['T'],theta,"_fit_Nf{}".format(int(Nf)) if Nf>0 else '',srg_f)
        Ham_3b_root ='Ham_{0}{1}{2}_A3{3}.bin_'.format('NN3N' if NN_only and
                                                       srg_evolution else
                                                       'NN','_3N-ind' if NNN_srg
                                                       else '','_NNN' if not NN_only else '','_from3eff' if NN_only else '')
        Ham_2b_root ='Ham_{0}_A3_from3eff.bin_'.format('NN3N' if NN_only and srg_evolution else 'NN')
        Trel_A3_root='Trel{0}_A3_from3{1}.bin_'.format('_NN3N' if NN_only and srg_evolution else '','bre' if not NN_only else 'eff')

        Ham_3b_file = Ham_3b_root +extension_out
        Ham_2b_file = Ham_2b_root +extension_out
        Trel_A3_file= Trel_A3_root+extension_out

        tmp_ext="_N{0}_J{1}_T{2}".format(N_m,NUCLEI_prop[nuclei]['J'],NUCLEI_prop[nuclei]['T'])
        Ham_list=glob.glob(Ham_3b_file.replace('_J{:d}'.format(NUCLEI_prop[nuclei]['J']),'_*'))
        Ham_list.sort(reverse=True)
        print("existing ham_list",Ham_list)
        try:
            J_file_list=[int(re.search(r'_J([0-9]*)_hw',Ham_list[-1]).group(1))]
        except:
            J_file_list=[-1]

        if J_file_list[0]>=0:
            Ham_3b_file=Ham_list[0]
            J_file=J_file_list[0]
#            new_extension='{0}{1}{2}_{3:d}.{4:d}{5}'.format(nuclear_pot,'-srg' if srg_evolution else '', '{:3.1f}'.format(Lambda) if srg_evolution else '',N_max_file,int(hw),'_theta_{:0.2f}'.format(theta) if theta != 0 else '')
            new_extension='{0}{1}{2}{3}_srgf{12}_{4}-T12_{5}_N{6:d}_hw{7:d}_J{8}_T{9}_theta{10:0.3f}{11}'.format(nuclear_pot,'_srg' if srg_evolution else '', '{:3.1f}'.format(Lambda) if srg_evolution else '','-4Nsrg' if srg_4N_evolution else '',nuclei,NUCLEI_prop['3H']['T12_max_3NF'],N_max_file,int(hw),J_file,NUCLEI_prop[nuclei]['T'],theta,"_fit_Nf{}".format(int(Nf)) if Nf>0 else '',srg_f)
            Trel_A3_file= Trel_A3_root+new_extension
            if not os.path.isfile(Trel_A3_file) : sys.exit('Incomplete Kinetic file for N_max {:d} Hamiltonian'.format(N_max_file))
        else:
            modify_manychir_input_file('3H',N_m_range[-1],hw,theta,NNN_calc=not
            #NN_only,NNN_srg=any([srg_evolution,NNN_srg])) # modify J12 to not be
            #NN_only,NNN_srg=NNN_srg,Nmax_3N=[56,50,44,42,38]) # modify J12 to not be
            #NN_only,NNN_srg=NNN_srg,Nmax_3N=[20,20,20,20],J12_max_3NF=1,threeff=True
            NN_only,NNN_srg=NNN_srg,Nmax_3N=[56,44,38,36],J12_max_3NF=7,threeff=True
                                      if srg_evolution else False) # modify J12 to not be
            # the default, be larger

            try:
                outputs={'many.out':'many_{0}{1}{2}_{3}-T12_{4}_{5:d}.{6:d}{7}_dep_Ham_3b'.format(nuclear_pot,'-srg' if srg_evolution else '', '{:3.1f}'.format(Lambda) if srg_evolution else '',nuclei,NUCLEI_prop['3H']['T12_max_3NF'],N_m_range[-1],hw,'_theta_{:0.2f}'.format(theta) if theta != 0 else ''),
                         'Ham_{0}{1}_A3{2}.bin'.format('NN3N' if NN_only and srg_evolution else 'NN','_NNN' if not NN_only else '','_from3eff' if NN_only else ''):Ham_3b_file,
                         'Trel{0}_A3_from3{1}.bin'.format('_NN3N' if NN_only and srg_evolution else '','bre' if not NN_only else 'eff')  :Trel_A3_file
        ,'manyeff.log':'manyeff.log'+extension_out}
                #
                create_temporary_runspace()
                run_code(outputs)
                clean_temporary_runspace(outputs)
                print("&&&&&&&&&&&&&&&&&&&&&&\n finished calculating 3N interaction matrix elements")
            except OSError:
                print('Error while running tasks many_{0}{1}{2}_{3}-T12_{4}_{5:d}.{6:d}_dep_Ham_3b'.format(nuclear_pot,'-srg' if srg_evolution else '', '{:3.1f}'.format(Lambda) if srg_evolution else '',nuclei,NUCLEI_prop['3H']['T12_max_3NF'],N_m_range[-1],hw))


        #Ham_list=glob.glob(Ham_2b_file.replace('_{:d}'.format(N_m_range[-1]),'_*'))
        Ham_list=glob.glob(Ham_2b_file.replace('_N{:d}'.format(N_m_range[-1]),'_*'))
        Ham_list.sort(reverse=True)
        try:
            #N_max_file=int(re.search(r'_([0-9]*)$',Ham_list[0]).group(1))
            N_max_infiles=[int(re.search(r'_N([0-9]*)_hw',Ham_file).group(1)) for Ham_file in Ham_list]
        except:
            N_max_infiles=[0]
        if N_m_range[-1] in N_max_infiles:
            N_max_file=N_m_range[-1]
            Ham_2b_file=Ham_list[0]
        else:
            #modify_manychir_input_file('3H',N_m_range[-1],hw,theta,threeff=True,NNN_calc=False,NN_srg=srg_evolution,Nmax_3N=[26,20,20,20])
            #modify_manychir_input_file('3H',N_m_range[-1],hw,theta,threeff=True,NNN_calc=False,NN_srg=True,Nmax_3N=[20,20,20,20],J12_max_3NF=1)
            modify_manychir_input_file('3H',N_m_range[-1],hw,theta,threeff=True,NNN_calc=False,NN_srg=True,Nmax_3N=[56,44,38,36],J12_max_3NF=7)

            try:
                outputs={'many.out':'many_{0}{1}{2}_{3}-T12_{4}_{5:d}.{6:d}{7}_dep_Ham_2b'.format(nuclear_pot,'-srg' if srg_evolution else '', '{:3.1f}'.format(Lambda) if srg_evolution else '',nuclei,NUCLEI_prop['3H']['T12_max_3NF'],N_m_range[-1],hw,'_theta_{:0.2f}'.format(theta) if theta != 0 else ''),
                         #'Ham_{0}_A3_from3eff.bin'.format('NN3N' if NN_only and srg_evolution else 'NN'):Ham_2b_file}
                         'Ham_NN_A3_from3eff.bin':Ham_2b_file,'manyeff.log':'manyeff.log'+extension_out}# because NNN_SRG is false
                         #'Ham_NN_A3.bin':Ham_2b_file}
                create_temporary_runspace()
                run_code(outputs)
                clean_temporary_runspace(outputs)
                print("&&&&&&&&&&&&&&&&&&&&&&\n finished calculating 2N interaction matrix elements")
            except OSError:
                print('Error while running tasks many_{0}{1}{2}_{3}-T12_{4}_{5:d}.{6:d}_dep_Ham_2b'.format(nuclear_pot,'-srg' if srg_evolution else '', '{:3.1f}'.format(Lambda) if srg_evolution else '',nuclei,NUCLEI_prop['3H']['T12_max_3NF'],N_m_range[-1],hw))

        if (not NN_only or NNN_srg) : symlink_force(Ham_3b_file ,'Ham_NN_NNN.bin')
        symlink_force(Ham_2b_file ,'Ham_NN.bin')
        symlink_force(Ham_2b_file ,'Ham_NN_subtr.bin')
        symlink_force(Trel_A3_file,'Trel.bin')
        print("files",Ham_3b_file,Ham_2b_file,Trel_A3_file)
        return Ham_3b_file,Ham_2b_file,Trel_A3_file

    def run_3A_dependance():
#        run_params.update({'mscheme':'T','frequency_conversion':frequency_conversion,'n_omega':n_omega,'delta_omega':delta_omega})
        run_params.update({'mscheme':'T'})
        Ham_3b_file,Ham_2b_file,Trel_A3_file=run_4He_dependance()
        shutil.copy(os.path.join(home_depot+'/v3trans_cJ_omp.exe'),'v3trans.exe')

        input_file={hw:[Ham_3b_file,Ham_2b_file,Trel_A3_file]}
        if frequency_conversion:
            frequency=hw-delta_omega
            while (frequency>hw-n_omega*delta_omega-1.):
                files_in  =[ 'Ham_{0}{1}_A3{2}.bin_{3}_from_{4:d}'.format('NN3N' if NN_only and srg_evolution else 'NN','_NNN' if not NN_only else '','_from3eff' if NN_only else '', frequency,hw)
                            ,'Ham_{0}_A3_from3eff.bin_{1}_from_{2:d}'.format('NN3N' if NN_only and srg_evolution else 'NN',frequency,hw)
                            ,'Trel{0}_A3_from3{1}.bin_{2}_from_{3:d}'.format('_NN3N' if NN_only and srg_evolution else '','bre' if not NN_only else 'eff',frequency,hw)]
                files_save=[Ham_3b_file+'_{0}_from_{1:d}'.format(frequency,hw),
                            Ham_2b_file+'_{0}_from_{1:d}'.format(frequency,hw),
                            Trel_A3_file+'_{0}_from_{1:d}'.format(frequency,hw)]

                for output,new_name in zip(files_in,files_save):
                    ren(output,new_name)
                input_file.update({frequency:files_save})
                frequency-=delta_omega

        return input_file

    def run_3A(files_in):
        if not files_in: return
        for frequency in files_in.keys():
            Ham_3b_file,Ham_2b_file,Trel_A3_file=files_in[frequency]

            if nucleons<=5:
               N_1max=N_m
               N_2max=N_m
            elif nucleons==6:
               N_1max=N_m-1
               N_2max=N_m
            else:
               N_1max=N_m-2
               N_2max=N_m-1

            extension_out='{0}{1}{2}{3}_srgf{13}_{4}-T12_{5}_N{6:d}_hw{7:d}_to{8:d}_J{9}_T{10}_theta{11:0.3f}{12}'.format(nuclear_pot,'_A3srg' if NNN_srg else'_srg' if srg_evolution else '', '{:3.1f}'.format(Lambda) if srg_evolution else '','-4Nsrg' if srg_4N_evolution else '',nuclei,NUCLEI_prop[nuclei]['T12_max_3NF'],N_m,int(hw),int(frequency),NUCLEI_prop[nuclei]['J'],NUCLEI_prop[nuclei]['T'],theta,"_fit_Nf{}".format(int(Nf)) if Nf>0 else '',srg_f)

            #extension_out='{0}{1}{2}{3}_{4}-T12_{5}_{6:d}.{7:d}_to.{8:d}{9}'.format(nuclear_pot,'-srg' if srg_evolution else '', '{:3.1f}'.format(Lambda) if srg_evolution else '','-4Nsrg' if srg_4N_evolution else '',nuclei,NUCLEI_prop[nuclei]['T12_max_3NF'],N_m,int(hw),int(frequency),'_theta_{:0.2f}'.format(theta) if theta != 0 else '')

            v3trans_file_int='v3trans_J3T3.int_'+extension_out

            if not os.path.isfile(v3trans_file_int) and not os.path.isfile(v3trans_file_int+'.gz'):


                run_params.update({'threeffcal':True})
                modify_manychir_input_file('3H',N_m,frequency,theta,NNN_calc=not NN_only,NNNN_srg=srg_4N_evolution,threeff=True,Nmax_3N=[24,24,20,20],J12_max_3NF=1)
                try:
                    outputs={'many.out':'many_'+extension_out,
                             'Density_at0_A4.bin':'Density_at0_A4.bin_'+extension_out,
                             'manyobdvec.tmp':'manyobdvec.tmp_'+extension_out}
                    if not NN_only or NNN_srg : symlink_force(Ham_3b_file ,'Ham_NN_NNN.bin')
                    symlink_force(Ham_2b_file ,'Ham_NN.bin')
                    symlink_force(Trel_A3_file,'Trel.bin')

                    run_code(outputs)

                    key_seq=['protons','neutrons','hbarxomega','N3max','N2max','N1max','complex_scaling']
                    formatted_write={'protons':'{:<13d}! protons\n'.format(NUCLEI_prop[nuclei]['protons']),'neutrons':'{:<13d}! neutrons\n'.format(NUCLEI_prop[nuclei]['neutrons']),'hbarxomega':'{:<13E}! hbar*Omega\n'.format(frequency),'N3max':'{:<13d}! N_3max\n'.format(N_m),'N2max':'{:<13d}! N_2max\n'.format(N_2max),'N1max':'{:<13d}! N_1max\n'.format(N_1max),'complex_scaling':'{} ! complex_scaling\n'.format('.true.' if (theta>0) else '.false') }
                    with open('v3trans.in','w') as f:
                        for items in key_seq:
                            f.write(formatted_write[items])

                    try:
                        tcoef_file='v3trans_tcoef_{0:d}{1:d}{2:d}.sav'.format(N_m,N_2max,N_1max)
                        tcoef_path_file=os.path.join(mother_path,tcoef_file)
                        if os.path.isfile(tcoef_path_file):
                            symlink_force(tcoef_path_file,tcoef_file)
                        with open(os.devnull, 'w') as fp:
                            try:
                                #subprocess.call(['./v3trans.exe'],stdout=fp)
                                if(output_stream):
                                    subprocess.call(['./v3trans.exe'])
                                else:
                                    subprocess.call(['./v3trans.exe'],stdout=fp)
                                if not os.path.isfile(tcoef_path_file):
                                    shutil.copy(tcoef_file,tcoef_path_file)
                            except:
                                raise
                        ren('v3trans.out','v3trans.'+extension_out)
                        ren('v3trans_J3T3.int',v3trans_file_int)
                    except OSError:
                        print('3N force file was not converted for frequency={:d}'.format(frequency))
                        pass
                except OSError:
                    print('Calculation was not performed for frequency={:d}'.format(frequency))
                    pass
    def read_basis_from_file():
        files=glob.glob("*.tmp")
        for file in files:
            rem(file)
        tmp_ext="_{0}{3}_J{1}_T{2}".format(nuclei,NUCLEI_prop[nuclei]['J'],NUCLEI_prop[nuclei]['T'],'_eff'if (NNN_srg or not NN_only) else '')
        files=glob.glob("*.tmp"+tmp_ext+"_N*")
        files.sort(reverse=True)
        readtmp=False
        try:
            N_max_file=np.array([int(re.search(r'_N([0-9]*)$',file).group(1)) for file
            in files])
        except:
            N_max_file=[0]
        if(len(N_max_file)==0):
            N_max_file=[0]
        if (N_m in N_max_file):
            readtmp=True
            files=glob.glob("*.tmp"+tmp_ext+"_N{}".format(N_m))
            for file in files:
                initial_name=file.split(tmp_ext)[0]
                symlink_force(file ,initial_name)
        if( readtmp):
            tmp_ext=""
        else:
            tmp_ext=tmp_ext+"_N{}".format(N_m)
        return readtmp,tmp_ext

    def run_4He(dependency_files=None):
        output_NNN_interaction=False
#
        if((NNN_srg or NNNN_srg )and False):
            readtmp=False
            tmp_ext=''
        else:
            readtmp,tmp_ext=read_basis_from_file()
        print("reading basis params from file set to :",readtmp,tmp_ext)
        #modify_manychir_input_file(nuclei,N_m,hw,theta,NNN_calc=not NN_only,NNN_srg=any([srg_evolution,NNN_srg]),NNNN_srg=srg_4N_evolution)
        modify_manychir_input_file(nuclei,N_m,hw,theta,NNN_calc=not
         NN_only,NNN_srg=NNN_srg,NNNN_srg=srg_4N_evolution,threeff=False if (NN_only and not NNN_srg) or nuclei[0]<'4' else True,readtmp=readtmp)

        #extension_out='{0}{1}{2}{3}_{4}-T12_{5}_{6:d}.{7:d}{8}'.format(nuclear_pot,'-srg' if srg_evolution else '', '{:3.1f}'.format(Lambda) if srg_evolution else '','-4Nsrg' if srg_4N_evolution else '',nuclei,NUCLEI_prop['3H']['T12_max_3NF'],N_m,int(hw),'_theta_{:0.2f}'.format(theta) if theta != 0 else '')
        extension_out='{0}{1}{2}{3}_srgf{12}_{4}-T12_{5}_N{6:d}_hw{7:d}_J{8}_T{9}_theta{10:0.3f}{11}'.format(nuclear_pot,'_A3srg' if NNN_srg else '_srg' if srg_evolution else '', '{:3.1f}'.format(Lambda) if srg_evolution else '','-4Nsrg' if srg_4N_evolution else '',nuclei,NUCLEI_prop['3H']['T12_max_3NF'],N_m,int(hw),NUCLEI_prop[nuclei]['J'],NUCLEI_prop[nuclei]['T'],theta,"_fit_Nf{}".format(int(Nf)) if Nf>0 else '',srg_f)

        try:
            outputs={'many.out':'many_'+extension_out,
                     'Ham_NN_3N_A4.bin':'Ham_NN_3N_A4.bin_'+extension_out,
                     'Trel_A4.bin':'Trel_A4.bin_{0:d}.{1:d}'.format(N_m,hw),
                     'Density_at0_A4.bin':'Density_at0_A4.bin_'+extension_out,
                     'manyobdvec.tmp':'manyobdvec.tmp_'+extension_out,'manyeff.log':'manyeff.log'+extension_out}

            create_temporary_runspace(dependency_files)
            run_code(outputs)
            clean_temporary_runspace(outputs,tmp_ext)
            print("&&&&&&&&&&&&&&&&&&&&&&\nfour-body calculation is finished")
        except OSError:
            print('Error:4He-Calculation was not performed for N_max={:d}'.format(N_m))
            pass

    def run_3H(*args):
        output_NNN_interaction=False
        #modify_manychir_input_file(nuclei,N_m,hw,theta,NNN_calc=not NN_only,NNN_srg=any([srg_evolution,NNN_srg]),NNNN_srg=srg_4N_evolution)
        modify_manychir_input_file(nuclei,N_m,hw,theta,NNN_calc=not NN_only,NNN_srg=NNN_srg,NNNN_srg=srg_4N_evolution)
        print("NNN_srg is:",NNN_srg)
        #extension_out='{0}{1}{2}{3}_{4}-T12_{5}_{6:d}.{7:d}{8}'.format(nuclear_pot,'-srg' if srg_evolution else '', '{:3.1f}'.format(Lambda) if srg_evolution else '','-4Nsrg' if srg_4N_evolution else '',nuclei,NUCLEI_prop['3H']['T12_max_3NF'],N_m,int(hw),'_theta_{:0.2f}'.format(theta) if theta != 0 else '')
        extension_out='{0}{1}{2}{3}_srgf{12}_{4}-T12_{5}_N{6:d}_hw{7:d}_J{8}_T{9}_theta{10:0.3f}{11}'.format(nuclear_pot,'_srg' if srg_evolution else '', '{:3.1f}'.format(Lambda) if srg_evolution else '','-4Nsrg' if srg_4N_evolution else '',nuclei,NUCLEI_prop['3H']['T12_max_3NF'],N_m,int(hw),NUCLEI_prop[nuclei]['J'],NUCLEI_prop[nuclei]['T'],theta,"_fit_Nf{}".format(int(Nf)) if Nf>0 else '',srg_f)


        try:
            outputs={'many.out':'many_'+extension_out,
                     'Ham_NN_3A.bin':'Ham_NN_3A.bin_'+extension_out,
                     'Trel_3A.bin':'Trel_3A.bin_{0:d}.{1:d}'.format(N_m,hw),
                     'Density_at0_A4.bin':'Density_at0_A4.bin_'+extension_out,
                     'manyobdvec.tmp':'manyobdvec.tmp_'+extension_out}

            run_code(outputs)
        except OSError:
            print('Error:3H-Calculation was not performed for N_max={:d}'.format(N_m))
            pass

    LEC_NNN={'n3lo'     :{'icutnum':100,'chirp':True,'c_1':-0.81  ,'c_3':-3.20  ,'c_4':5.40,'NN-only':False},
            'n2loopt'   :{'icutnum':102,'chirp':True,'c_1':-0.9186,'c_3':-3.8887,'c_4':4.3103,'NN-only':False},
            'n2losat'   :{'icutnum':105,'chirp':True,'c_1':-1.121521199632590,'c_3':-3.925005856486820,'c_4':3.765687158585920,'NN-only':False},
            'n4lo500'   :{'icutnum':106,'chirp':True,'c_1':-0.73  ,'c_3':-3.38   ,'c_4':1.69,'NN-only':False},
            'n4lo500new':{'icutnum':107,'chirp':True,'c_1':-1.20  ,'c_3':-4.43   ,'c_4':2.69,'NN-only':False},
            'n2lo500'   :{'icutnum':108,'chirp':True,'c_1':-0.74  ,'c_3':-3.61   ,'c_4':2.49,'NN-only':False},
            'nlo500'    :{'icutnum':109,'chirp':True,'c_1':-0.73  ,'c_3':-3.38   ,'c_4':1.69,'NN-only':True},
            'lo500'     :{'icutnum':110,'chirp':True,'c_1':0.0    ,'c_3': 0.0    ,'c_4':0.09,'NN-only':True},
            'mtv'     :{'icutnum':116,'chirp':False,'c_1':0.0    ,'c_3': 0.0    ,'c_4':0.09,'NN-only':True},
            'pless-lo'  :{'icutnum':116,'chirp':True,'c_1':0.0    ,'c_3': 0.0    ,'c_4':0.0 ,'NN-only':False}}

    mother_path     =os.getcwd()
    home            =os.environ['HOME']
    path_executable =mother_path
    if 'cca' in socket.getfqdn():
        home_depot      =os.path.join(home,'Codes/cs-nn-interaction')
        home_depot      =os.path.join(mother_path,'')
    elif os.getenv('MACHINE_NAME') is not None and 'idris' in os.getenv('MACHINE_NAME'):
        home_depot      =os.path.join(home,'Codes/cs-nn-interaction')
    else:
        home_depot      =os.path.join(home,'Documents/Codes/Ncsm-interactions')
        home_depot      =os.path.join(mother_path,'')

    #exefile='manyeffv3b_SRG_PVPgen.exe'
    exefile='manyeffv3b.exe'

    os.environ["OMP_PLACES"] = "cores"
    os.environ["OMP_NUM_THREADS"] = "240"

    # for genuine NNN switched off set to 1
    if NNN_interaction or ( NNN_interaction == ['pless'] and NNN_interaction is not False):
        NN_only=False
    else:
        NN_only=True

    # for bare interaction with no SRG evolution
    nonlocal_regulator=False
    NNN_cutoff_nonlocal_Lambda=0.0
    lambda_range=[]
    try:
        srg_evolution=run_params['srg_evolution']
    except:
        srg_evolution=False
    try:
        NNN_srg=run_params['NNN_srg']
    except:
        NNN_srg=False
    try:
        NNNN_srg=run_params['NNNN_srg']
    except:
        NNNN_srg=False
    try:
        srg_4N_evolution=run_params['srg_4N_evolution']
    except:
        srg_4N_evolution=False
    try:
        lambda_range =run_params['lambda_range']
    except:
        lambda_range =[0.0]
    try:
        N_m_range    =run_params['N_m_range']
    except:
        raise RuntimeError('No Nmax in input')
    try:
        hw_range     =run_params['hw_range']
    except:
        raise RuntimeError('No hw in input')
    try:
        NNN_cutoff_Lambda=NNN_interaction['cutoff'][0]
    except:
        if NN_interaction == ['pless-lo'] or not NNN_interaction:
            NNN_cutoff_Lambda=0.0
        else:
            raise RuntimeError('No 3N cutoff in input')
    try:
        mscheme=run_params['mscheme']
    except:
        mscheme=False
        run_params.update({'mscheme':'F'})
    try:
        frequency_conversion=run_params['frequency_conversion']
    except:
        frequency_conversion=False
    try:
        n_omega=run_params['n_omega']
    except:
        n_omega=1
    try:
        delta_omega=run_params['delta_omega']
    except:
        delta_omega=1.0
    try:
        HO_trap      =run_params['HO_trap']
    except:
        try:
            hw_trap_range=run_params['hw_trap']
            A_trap_range =run_params['A_trap']
            HO_trap      =True
        except:
            HO_trap      =False
    try:
        theta_range= run_params['theta_range']
    except:
        theta_range=[0.0]
    try:
        Nf=run_params['Nf']
    except:
        Nf=-1

    try:
        srg_f =run_params['srg_f']
    except:
        srg_f =1.0
    if srg_4N_evolution and srg_evolution:
        sys.exit('Error in SRG input')
    if srg_4N_evolution:
        mscheme=True
    if NNN_interaction and NNN_interaction['regulator'] != 'local':
        nonlocal_regulator=True
        NNN_cutoff_nonlocal_Lambda=NNN_interaction['cutoff'][-1]

    c_TM = 0.0
    c_ON = 0.0
    x_78 = 0.0
    c_D_range  =NNN_interaction.get('c_D',[0.0]) if not NN_only else [0.0]
    c_E_range  =NNN_interaction.get('c_E',[0.0]) if not NN_only else [0.0]
    c_LS_range =NNN_interaction.get('E8' ,[0.0]) if not NN_only else [0.0]

    # observable and density calculation
    observable_obd_calculation=False
    density_at_r0=False

    nucleons=NUCLEI_prop[nuclei]['neutrons']+NUCLEI_prop[nuclei]['protons']
    if nucleons>4:
        shell='p'
    else:
        shell=''

    if 'pless-lo' in NN_interaction:
        key_length =len(pless_params.keys())
        item_length=len(list(pless_params.items())[0][1]) if type(list(pless_params.items())[0][1]) is list else 1
        pless_int  =list([ list(zip(keys,values)) for keys,values in zip(np.stack([list(pless_params.keys())]*item_length).tolist(),np.reshape([make_list(pless_params[keys]) for keys in pless_params.keys()],(key_length,item_length)).T.tolist())])
        if NN_only:
            pless_params.pop('power')
            pless_params.pop('H0')

    base_param   = {'hw':hw_range,'N_m':N_m_range,'theta':theta_range,'lambda':lambda_range}
    base_param.update({'hw_trap':hw_trap_range,'A_trap':A_trap_range} if HO_trap else {})
    chiral_param = {'regulator' :[NNN_interaction['regulator'] if not NN_only else 0.0],'cutoff':NNN_interaction['cutoff'][:-1] if not NN_only else [0.0],'c_D':c_D_range,'c_E':c_E_range,'E8':c_LS_range}
    calc_param   = {'n3lo'      :{**chiral_param,**base_param},
                    'n2loopt'   :{**chiral_param,**base_param},
                    'n2losat'  :{**chiral_param,**base_param},
                    'n4lo500'   :{**chiral_param,**base_param},
                    'n4lo500new':{**chiral_param,**base_param},
                    'n2lo500'   :{**chiral_param,**base_param},
                    'nlo500'    :{'hw':hw_range,'N_m':N_m_range,'theta':theta_range},
                    'lo500'     :{'hw':hw_range,'N_m':N_m_range,'theta':theta_range},
                    'mtv'     :{**chiral_param,**base_param},
                    'pless-lo'  :{'interactions':pless_int ,'hw':hw_range,'N_m':N_m_range,'theta':theta_range} if pless_params else None}
    first_run=True
    for NN in NN_interaction:
        param, val   = zip(*calc_param[NN].items())

        for run_dict in tqdm([dict(zip(param, v)) for v in it.product(*val)]):
            if NN == 'pless-lo':
                run_dict.update(dict(run_dict['interactions']))
            cutoff=run_dict['cutoff']
            c_D   =run_dict.get('c_D' ,0.0)
            c_E   =run_dict.get('c_E' ,0.0)
            c_LS  =run_dict.get('E8'  ,0.0)
            hw    =run_dict.get('hw'  ,0.0)
            N_m   =run_dict.get('N_m' ,0.0)
            theta=run_dict.get('theta',0.0)
            Lambda=run_dict.get('lambda',2.0)

            NN_pot  ='NN'+NN
            NNN_pot ='3N' if not NN_only else ''
            NNN_keys={'n3lo'      :NNN_interaction.keys() if not NN_only else [],
                      'n2loopt'   :NNN_interaction.keys() if not NN_only else [],
                      'n2losat'  :NNN_interaction.keys() if not NN_only else [],
                      'n4lo500'   :NNN_interaction.keys() if not NN_only else [],
                      'n4lo500new':NNN_interaction.keys() if not NN_only else [],
                      'n2lo500'   :NNN_interaction.keys() if not NN_only else [],
                      'nlo500'    :[],
                      'lo500'     :[],
                      'mtv'     :[],
                      'pless-lo'  :pless_params.keys() if pless_params else None}
            chirform={'c_D':'{:.1f}','c_E':'{:.1f}','E8':'{:.1f}','cutoff':'{:.0f}','regulator':'{}'}
            pless_form={'1S0':'{:.2e}','3S1':'{:.2e}','H0':'{:.2e}','power':'{:.0f}','cutoff':'{:.1f}','regulator':'{}'}
            pless_unit={'1S0':float,'3S1':float,'H0':float,'power':int,'cutoff':float,'regulator':str}
            dir_form={'n3lo'      :chirform,
                      'n2loopt'   :chirform,
                      'n2losat'  :chirform,
                      'n4lo500'   :chirform,
                      'n4lo500new':chirform,
                      'n2lo500'   :chirform,
                      'nlo500'    :None,
                      'lo500'     :None,
                      'mtv'     :None,
                      'pless-lo'  :pless_form}

            for key in NNN_keys[NN]:
                formatted_str='_{}'+dir_form[NN][key]
                if NN == 'pless-lo': run_dict[key]=pless_unit[key](run_dict[key])
                NNN_pot=NNN_pot+formatted_str.format(key.replace('_','') if all([key != 'regulator',key != 'cutoff']) else '',run_dict[key])
            nuclear_pot=NN_pot+('_'+NNN_pot if NNN_pot else '')

            working_dir=nuclear_pot+('_SRG' if srg_evolution else
            '')+"_"+nuclei+("_testing" if testing else'')#+"_test2"#+("_Nf{}".format(Nf) if Nf>0 else "")
            #new_path='/vol0'
            working_path=os.path.join(mother_path,working_dir)
            #working_path=os.path.join(new_path,working_dir)
            print("moving into directory:",working_path)
            if not os.path.isdir(working_path):
                os.mkdir(working_path)
                os.chdir(working_path)
            else:
                os.chdir(working_path)
            shutil.copy(os.path.join(home_depot+'/manychir-template.in'),'manychir_v1.in')
            #fitting parameters
            source_path=os.path.join(home_depot,"2b_fitting/lmfit_parameters/")
            destination_path=os.path.join("","2b_fitting/lmfit_parameters/")
            #
            #shutil.copy(os.path.join(path_executable,exefile),os.getcwd())
            if(first_run):
                shutil.copy(os.path.join(home_depot,exefile),os.getcwd())
                if(os.path.exists(destination_path)):
                    shutil.rmtree(destination_path)
                shutil.copytree(source_path,destination_path)
                print("copied the execution file")
                first_run=False
            if NN == 'pless-lo':
                if not pless_params: sys.exit('No parameters for pionless interaction')
                shutil.copy(os.path.join(home_depot+'/pless-params-template.in'),'pless-params.in')
                modify_pless_input_file(**run_dict)
            if observable_obd_calculation:
                shutil.copy(os.path.join(path_executable,'manyobdtr.exe'),os.getcwd())
#*************
#            for filename in os.listdir(home_depot):
#                if filename.endswith('.a') or filename.endswith('.so') or filename.endswith('.mod'):
#                    source_path = os.path.join(home_depot, filename)
#                    destination_path = os.path.join('', filename)
#
#            # Copy the file to the destination directory
#                    shutil.copy2(source_path, destination_path)
#                    print(f"Copied: {filename}")

#*************
            if LEC_NNN[NN]:
                icut   =LEC_NNN[NN]['icutnum']
                chirp  =LEC_NNN[NN]['chirp']
                NN_only=LEC_NNN[NN]['NN-only'] if not NN_only else NN_only
            else:
               chirp   =False
               sys.exit('This interaction is not yet setup')
               NN_only =False

            if NN_only:
                output_NNN_interaction=False
                LEC_NNN[NN].update({'c_1':0.0  ,'c_3':0.0   ,'c_4':0.0})
            else:
                output_NNN_interaction=True

            #for Lambda in lambda_range:

            dependancies={'3A'      :run_3A_dependance,
                          '4He'     :run_4He_dependance,
                          '4H'      :run_4He_dependance,
                          '4N'      :run_4He_dependance,
                          '4Li'      :run_4He_dependance
                         }
            run         ={'3A'      :run_3A,
                          '3He'     :run_4He,
                          '3H'      :run_4He,
                          '4He'     :run_4He,
                          '4H'     :run_4He,
                          '4Li'     :run_4He,
                          '4N'     :run_4He,
                         }
            file_info=[]
            if(NNN_srg or not NN_only):
                file_info   = dependancies.get(nuclei, empty)()
            run.get(nuclei, empty)(file_info)

            if observable_obd_calculation:
                with open(os.devnull, 'w') as fp:
                    try:
                        if(output_stream):
                            subprocess.call(['./manyobdtr.exe'])
                        else:
                            subprocess.call(['./manyobdtr.exe'],stdout=fp)
                    except:
                        raise

                ren('manyobdtr.out','manyobdtr_'+extension_out)
                ren('manyobd.tmp','manyobd.tmp_'+extension_out)
                ren('radialdens_trinv.dat','radialdens_trinv_'+extension_out)

            if observable_obd_calculation:
                ren('summary.dat','summary_'+extension_out)

            os.chdir(mother_path)

    end_time=time.time()
    print("This calculation took:"+time_format(end_time-start_time))
    print("--------------------------")
# end of run_manyeff function





def run_manyeff_generator(frequency_conversion=False,n_omega=0,delta_Omega=0.0,NN_interaction=['n3lo'],NNN_interaction={'regulator':'lnl','cutoff':[650,500],'c_D':[0.7],'c_E':[-0.06],'E8':[0.0]},run_params={'srg_evolution':False,'srg_4Nevolution':False,'lambda_range':[0.0],'N_m_range':[16],'hw_range':[20]}):

    run_params.update({'mscheme':'T','frequency_conversion':frequency_conversion,'n_omega':n_omega,'delta_omega':delta_omega})
    run_manyeff('3A',NN_interaction=NN_interaction,NNN_interaction=NNN_interaction,run_params=run_params)

def run_ncsmpnv2beff_rgm(nuclei,NN_interaction=['n3lo'],run_params={'srg_evolution':True,'lambda_range':[2.0],'N_m_range':[14],'hw_range':[20],'theta_range':[0.0]},pless_params=None,tasks_achieved=False,direct_output=False):
    import fileinput
    import tempfile
    import itertools as it

    NUCLEI_prop={'2H' :{'protons':1,'neutrons':1,'T':0,'iso':False},
                 '2A0':{'protons':1,'neutrons':1,'T':0,'iso':True},
                 '2A1':{'protons':1,'neutrons':1,'T':1,'iso':True}}

    if tasks_achieved:
        tasks_done=[]
    else:
        tasks_done=None

    def modify_ncsmpnv2beff_input_file(nuclei_in,N_max,h_w,theta,N_r=80):
        import fileinput
        global f,line

        INPUT_ctrl={'protons' :                     ( lambda line : re.sub(r"^[0-9]*(.*$)",r'{:<2d}'.format(NUCLEI_prop[nuclei_in]['protons'])+"\g<1>",line)),
                    'neutrons' :                    ( lambda line : re.sub(r"^[0-9]*(.*$)",r'{:<2d}'.format(NUCLEI_prop[nuclei_in]['neutrons'])+"\g<1>",line)),
                    'theta':                        ( lambda line : re.sub(r"^[-.0-9]+[Ed]?[-+0-9]*(.*$)",r'{:.5f}'.format(theta)+"\g<1>",line)),
                    'hbar*Omega' :                  ( lambda line : re.sub(r"^[.0-9]+d[-+0-9]*(.*$)",r'{:.2E}'.format(h_w)+"\g<1>",line)),
                    'N_max (really n12_max)' :      ( lambda line : re.sub(r"^[0-9]* +[0-9]*(.*$)",r'{:<3d} {:<3d}'.format(N_max,N_r)+"\g<1>",line)),
                    'iso' :                         ( lambda line : re.sub(r"^[TF](.*$)",r'{}'.format(NUCLEI_prop[nuclei_in]['iso'])[0]+"\g<1>",line)),
                    'T for even A, 2*T for odd A' : ( lambda line : re.sub(r"^[0-9]*(.*$)",r'{:<2d}'.format(NUCLEI_prop[nuclei_in]['T'])+"\g<1>",line)),
                    'cdbpot' :                      ( lambda line : re.sub(r"^[TF](.*$)",r'{}'.format(not LEC_NN[NN]['chirp'])[0]+"\g<1>",line)),
                    'av18' :                        ( lambda line : re.sub(r"^[TF](.*$)",r'{}'.format(not LEC_NN[NN]['chirp'])[0]+"\g<1>",line)),
                    'MTV, Minnesota' :              ( lambda line : re.sub(r"^[TF](.*$)",r'{}'.format(not LEC_NN[NN]['chirp'])[0]+"\g<1>",line)),
                    'chiral potential' :            ( lambda line : re.sub(r"^[TF](.*$)",r'{}'.format(LEC_NN[NN]['chirp'])[0]+"\g<1>",line)),
                    'cutnum - 100 default' :        ( lambda line : re.sub(r"^[0-9]*(.*$)",r'{:<3d}'.format(icut)+"\g<1>",line)),
                    'srg' :                         ( lambda line : re.sub(r"^[TF]*(.*$)",r'{}'.format(srg_evolution)[0]+"\g<1>",line)),
                    'lambda of srg' :               ( lambda line : re.sub(r"^[.0-9]+(.*$)",r'{:.2E}'.format(Lambda)+"\g<1>",line)),
                    'HO_trap' :                     ( lambda line : re.sub(r"^[TF]*(.*$)",r'{}'.format(HO_trap)[0]+"\g<1>",line)),
                    'hbo_trap,nucleonsintrap' :     ( lambda line : re.sub(r"^[.0-9]+d[-+0-9]* +[0-9]*(.*$)",r'{:.2E} {:<2d}'.format(hw_trap,A_trap)+"\g<1>",line))}

        ctrl_pattern  =re.compile(r"! +(.*)")

        f=fileinput.input(files=('ncsmpneff_v1.in'),inplace=True)
        line=next(f)
        while True:
            try:
                ctrl_key=re.search(ctrl_pattern,line)
                try:
                    new_line=INPUT_ctrl[ctrl_key.group(1)](line)
                except KeyError:
                    new_line=line
                print(new_line,end='')
                line=next(f)
            except StopIteration:
                break
            except ValueError:
                continue
        f.close()

    def rem(file,afix=''):
        if os.path.isfile(file+afix):
            os.remove(file+afix)

    def ren(file,newname):
       if os.path.isfile(file):
          os.rename(file,newname)

    def rem_auto_output_files(file_2_avoid):
        import glob
        files=glob.glob("*.dat") + glob.glob("TBME*") + glob.glob("*.int") + glob.glob("*.out*")
        files=filter(lambda  elem:  not any(elem==item for item in file_2_avoid) ,files)
        for file in files:
            rem(file)

    def symlink_force(target, link_name):
        try:
            temp_link = link_name+'.new'
            rem(temp_link)
            os.symlink(target, temp_link)
            os.rename(temp_link, link_name)
        except OSError as e:
            raise

    def run_code(output_file={'ncsmpneff.out':'ncsmpneff.out_{:f}'.format(time.time())}):
        import time
        files_here=glob.glob("*.dat") + glob.glob("TBME*") + glob.glob("*.int*") + glob.glob("*.out*")
        symlink_force('ncsmpneff_v1.in','ncsmpneff.in')
        with open(os.devnull, 'w') as fp:
            try:
                if(direct_output):
                    subprocess.call(['./'+exefile])
                else:
                    subprocess.call(['./'+exefile],stdout=fp)
            except:
                raise
        for output,new_name in output_file.items():
            if not os.path.exists(output):
                print("Error: file:{} does not exist at the end of calculation".format(output))
            ren(output,new_name)
        rem('ncsmpneff','.in')
        rem_auto_output_files([items for keys,items in output_file.items()]+files_here)

    def empty(*args):
        pass

    def make_list(unknown_type):
        return [unknown_type] if type(unknown_type) is not list else unknown_type

    LEC_NN={'n3lo'      :{'icutnum':100,'chirp':True},
            'n2loopt'   :{'icutnum':102,'chirp':True},
            'n2losat'   :{'icutnum':105,'chirp':True},
            'n4lo500'   :{'icutnum':106,'chirp':True},
            'n4lo500new':{'icutnum':107,'chirp':True},
            'n2lo500'   :{'icutnum':108,'chirp':True},
            'nlo500'    :{'icutnum':109,'chirp':True},
            'lo500'     :{'icutnum':110,'chirp':True},
            'pless-lo'  :{'icutnum':116,'chirp':True}}

    mother_path     =os.getcwd()
    home            =os.environ['HOME']
    path_executable =mother_path
    if 'cca' in socket.getfqdn():
        home_depot      =os.path.join(home,'Codes/cs-nn-interaction')
    elif os.getenv('MACHINE_NAME') is not None and 'idris' in os.getenv('MACHINE_NAME'):
        home_depot      =os.path.join(home,'Codes/cs-nn-interaction')
    else:
        home_depot      =os.path.join(home,'Documents/Codes/Ncsm-interactions')
        home_depot      =os.path.join(mother_path,'')

    #exefile='ncsmpnv2beff.exe'
    exefile='ncsmpnv2beff_rgm.exe'

    os.environ["OMP_PLACES"] = "cores"
    os.environ["OMP_NUM_THREADS"] = "10"

    lambda_range=[]

    try:
        srg_evolution=run_params['srg_evolution']
    except:
        srg_evolution=False
    try:
        lambda_range =run_params['lambda_range']
    except:
        lambda_range =[0.0]
    try:
        theta_range= run_params['theta_range']
    except:
        theta_range=[0.0]
    try:
        N_m_range    =run_params['N_m_range']
    except:
        raise RuntimeError('No Nmax in input')
    try:
        hw_range     =run_params['hw_range']
    except:
        raise RuntimeError('No hw in input')
    try:
        HO_trap      =run_params['HO_trap']
        if HO_trap:
            hw_trap_range=run_params['hw_trap']
            A_trap_range =run_params['A_trap']
    except:
        try:
            hw_trap_range=run_params['hw_trap']
            A_trap_range =run_params['A_trap']
            HO_trap      =True
        except:
            HO_trap      =False
    try:
        N_r_max      =run_params['N_r_max']
    except:
        N_r_max      =None

    if 'pless-lo' in NN_interaction:
        if 'H0' in pless_params.keys(): pless_params.pop('H0')
        key_length =len(pless_params.keys())
        item_length=len(list(pless_params.items())[0][1]) if type(list(pless_params.items())[0][1]) is list else 1
        data_interaction=np.asarray([make_list(pless_params[keys]) for keys in pless_params.keys()],dtype=object)
        pless_int  =list([ list(zip(keys,values)) for keys,values in zip(np.stack([list(pless_params.keys())]*item_length).tolist(),np.reshape(data_interaction,(key_length,item_length)).T.tolist())])

    base_param   = {'hw':hw_range,'N_m':N_m_range,'theta':theta_range}
    base_param.update({'hw_trap':hw_trap_range,'A_trap':A_trap_range} if HO_trap else {})
    base_param.update({'hw_trap':hw_trap_range,'A_trap':A_trap_range} if HO_trap else {})
    calc_param   = {'n3lo'      :base_param,
                    'n2loopt'   :base_param,
                    'n2losat0'  :base_param,
                    'n4lo500'   :base_param,
                    'n4lo500new':base_param,
                    'n2lo500'   :base_param,
                    'nlo500'    :base_param,
                    'lo500'     :base_param,
                    'pless-lo'  : dict(**base_param,**{'interactions':pless_int}) if pless_params else {}}

    for NN in NN_interaction:
        param, val   = zip(*calc_param[NN].items())

        for run_dict in [dict(zip(param, v)) for v in it.product(*val)]:
            if NN == 'pless-lo':
                run_dict.update(dict(run_dict['interactions']))
            hw      =run_dict.get('hw' ,0.0)
            N_m     =run_dict.get('N_m',0.0)
            hw_trap =run_dict.get('hw_trap',0.0)
            A_trap  =run_dict.get('A_trap' ,0)
            theta=run_dict.get('theta',0.0)

            NN_pot    ='NN'+NN
            NN_keys   ={'n3lo'      :[],
                        'n2loopt'   :[],
                        'n2losat0'  :[],
                        'n4lo500'   :[],
                        'n4lo500new':[],
                        'n2lo500'   :[],
                        'nlo500'    :[],
                        'lo500'     :[],
                        'pless-lo'  :pless_params.keys() if pless_params else []}
            pless_form={'1S0':'{:.2e}','3S1':'{:.2e}','power':'{:.0f}','cutoff':'{:.1f}','regulator':'{}'}
            pless_unit={'1S0':float,'3S1':float,'power':int,'cutoff':float,'regulator':str}
            dir_form  ={'n3lo'      :None,
                        'n2loopt'   :None,
                        'n2losat0'  :None,
                        'n4lo500'   :None,
                        'n4lo500new':None,
                        'n2lo500'   :None,
                        'nlo500'    :None,
                        'lo500'     :None,
                        'pless-lo'  :pless_form}

            NN_params=''
            for key in NN_keys[NN]:
                formatted_str='_{}'+dir_form[NN][key]
                if NN == 'pless-lo': run_dict[key]=pless_unit[key](run_dict[key])
                NN_params=NN_params+formatted_str.format(key.replace('_',''),run_dict[key])
            nuclear_pot=NN_pot+(NN_params if NN_params else '')

            working_dir=nuclear_pot+('_SRG' if srg_evolution else '')+"_"+nuclei
            working_path=os.path.join(mother_path,working_dir)

            if not os.path.isdir(working_dir):
                try:
                    os.mkdir(working_dir)
                except:
                    pass
            print("going to dir:",working_path)
            os.chdir(working_path)

            with tempfile.TemporaryDirectory() as tmp_path:
                os.chdir(tmp_path)

                shutil.copy(os.path.join(home_depot+'/ncsmpneff-template.in'),'ncsmpneff_v1.in')
                #shutil.copy(os.path.join(path_executable,exefile),os.getcwd())
                shutil.copy(os.path.join(home_depot,exefile),os.getcwd())

                source_path=os.path.join(home_depot,"2b_fitting/lmfit_parameters/")
                destination_path=os.path.join("","2b_fitting/lmfit_parameters/")
                if(os.path.exists(destination_path)):
                    shutil.rmtree(destination_path)
                shutil.copytree(source_path,destination_path)

                if NN == 'pless-lo':
                    if not pless_params: sys.exit('No parameters for pionless interaction')
                    shutil.copy(os.path.join(home_depot+'/pless-params-template.in'),'pless-params.in')
                    delete_val_pless_input_file(**{'H0':None})
                    modify_pless_input_file(**run_dict)

                if LEC_NN[NN]:
                    icut   =LEC_NN[NN]['icutnum']
                    chirp  =LEC_NN[NN]['chirp']
                else:
                   chirp   =False
                   sys.exit('This interaction is not yet setup')

                file_names = []

                for Lambda in lambda_range:

                    inputs=(nuclei,N_m,hw,theta)
                    if(not (N_r_max is None)):
                        inputs=(nuclei,N_m,hw,theta,N_r_max)
                    #modify_ncsmpnv2beff_input_file(*inputs  if N_r_max is None else (inputs,{'N_r':N_r_max}) )
                    modify_ncsmpnv2beff_input_file(*inputs)

                    extension_out='{0}{1}{2}_{3}_{4:.0f}.{5:d}{6}{7}_theta{8:0.3f}'.format(nuclear_pot,'-srg' if srg_evolution else '', '{:3.1f}'.format(Lambda) if srg_evolution else '',nuclei,N_m,hw,'_A{:d}'.format(A_trap) if HO_trap else '','_{:.0f}'.format(hw_trap) if HO_trap else '',theta)


                    try:
                        outputs={'ncsmpneff.out':'ncsmpneff.out_'+extension_out,
                                 'TBME.int':'TBME_'+extension_out+'.int',
                                 'ncsmpneff_v1.in':'ncsmpneff_'+extension_out+'.in'
                                }

                        run_code(outputs)
                    except OSError:
                        print('Calculation was not performed for N_max={:d}'.format(N_m))
                        pass

                    file_names.extend([files for key,files in outputs.items()])
                    if tasks_achieved:
                        tasks_done.append({'tasks': extension_out, 'directory': working_path })

                for file_name in file_names:
                    try:
                        shutil.move(os.path.join(tmp_path, file_name), working_path)
                    except FileNotFoundError:
                        print('File named '+file_name+' was not generated')
                    except OSError as err:
                        if os.path.exists(os.path.join(working_path,file_name)):
                            print('Work already performed',file_name)
                            break
                        else:
                            print(err)

            os.chdir(mother_path)
    return tasks_done

def run_wrapper_TBME(dict_in,nuclei=None,NN_interaction=None,pless_params=None,tasks_returned=None):
    import shlex
    import subprocess
    command = shlex.split("env -i bash -c 'source /import/divers/intel/2017/parallel_studio_xe_2018.2.046/bin/psxevars.sh intel64 && env'")
    proc    = subprocess.Popen(command, stdout = subprocess.PIPE, stderr=subprocess.DEVNULL)
    for line in proc.stdout:
        (key, _, value) = line.decode().partition("=")
        os.environ[key] = value
    proc.communicate()

    tasks_returned.append(run_ncsmpnv2beff_rgm(nuclei,NN_interaction=NN_interaction,run_params=dict_in,pless_params=pless_params,tasks_achieved=True))
    return

async def TBME_generator(NN_interaction=['n3lo'],run_params={'srg_evolution':True,'lambda_range':[2.0],'N_m_range':[14],'hw_range':[20]},pless_params=None,multi_host=False):
    import os,socket
    from  tqdm import tqdm
    from paramiko import SSHClient
    from scp import SCPClient
    from itertools import product
    from joblib import Parallel, delayed
    from multiprocessing import cpu_count,Process,Manager,Pool
    from functools import partial

    def from_list(iterables):
        if isinstance(iterables,list):
            for elems in iterables:
                yield from from_list(elems)
        else:
            yield iterables

    @ipp.interactive
    def identify():
        import os
        import socket
        return {'host': socket.gethostname(), 'pid': os.getpid()}

    @ipp.interactive
    def locate():
        import os
        import socket
        return (socket.gethostname(), os.getcwd())

    @ipp.interactive
    def hello_tasks(tasks_in):
        import os
        import socket
        return ({'host': socket.gethostname(),'status': 'done','directory':os.getcwd()}, tasks_in)

    def make_parallel_wrapper(tasks_in):
        import socket
        from functools import partial
        from multiprocessing import cpu_count,Process,Manager,Pool,get_context

        def make_list(unknown_type):
            return [unknown_type] if type(unknown_type) is not list else unknown_type

        if isinstance(tasks_in, tuple):
            results      =[]
            Lambda,N_m,hw,hw_trap,A_trap = tasks_in
            run_wrapper_TBME({'srg_evolution':run_params['srg_evolution'],'lambda_range':[Lambda],'N_m_range':[N_m],'hw_range':[hw],'hw_trap':[hw_trap],'A_trap':[A_trap] },nuclei=nuclei, NN_interaction=NN_interaction, pless_params=pless_params,tasks_returned=results)
            tasks_in=make_list(tasks_in)
        else:
            manager = Manager()
            results = manager.list()
            with get_context("spawn").Pool(processes=int(cpu_count()/20)) as pool:
                parallel_out = pool.imap( partial(run_wrapper_TBME, nuclei=nuclei, NN_interaction=NN_interaction, pless_params=pless_params,tasks_returned=results) , [{'srg_evolution':run_params['srg_evolution'],'lambda_range':[Lambda],'N_m_range':[N_m],'hw_range':[hw],'hw_trap':[hw_trap],'A_trap':[A_trap] } for Lambda,N_m,hw,hw_trap,A_trap in tasks_in ] )

        directories=[]
        compute_task=[]
        for items in results:
            for subtasks in items:
                directories.append(subtasks['directory'].replace('/vol0/hupin/Runs',''))
                compute_task.append(subtasks['tasks'])
        if not results:
            return list(({'host': socket.gethostname(),'status': 'failed','directory':directories, 'compute_task':compute_task},b) for b in tasks_in)
        else:
            return list(({'host': socket.gethostname(),'status': 'done'  ,'directory':directories, 'compute_task':compute_task},b) for b in tasks_in)

    nuclei='2H'

    try:
        srg_evolution=run_params['srg_evolution']
    except:
        srg_evolution=False
    try:
        lambda_range =run_params['lambda_range']
    except:
        lambda_range =[0.0]
    try:
        N_m_range    =run_params['N_m_range']
    except:
        raise RuntimeError('No Nmax in input')
    try:
        hw_range     =run_params['hw_range']
    except:
        raise RuntimeError('No hw in input')
    try:
        HO_trap      =run_params['HO_trap']
        if HO_trap:
            hw_trap_range=run_params['hw_trap']
            A_trap_range =run_params['A_trap']
    except:
        try:
            hw_trap_range=run_params['hw_trap']
            A_trap_range =run_params['A_trap']
            HO_trap      =True
        except:
            HO_trap      =False

    tasks=list(product(lambda_range,N_m_range,hw_range,hw_trap_range,A_trap_range))

    if not multi_host:
# This is specific to hyperthreaded env.
        n_cores=int(cpu_count()/2)
        #Parallel(n_jobs=n_cores,prefer="threads")(delayed(run_ncsmpnv2beff_rgm) (nuclei,NN_interaction=NN_interaction,run_params={'srg_evolution':run_params['srg_evolution'],'lambda_range':[Lambda],'N_m_range':[N_m],'hw_range':[hw],'hw_trap':[hw_trap],'A_trap':[A_trap]},pless_params=pless_params) for Lambda,N_m,hw,hw_trap,A_trap in tqdm.tqdm(tasks,total=len(lambda_range)*len(N_m_range)*len(hw_range)*len(hw_trap_range)*len(A_trap_range)))
        manager = Manager()
        results = manager.list()

        processes = []
        with Pool(processes=max(int(n_cores/16),2)) as pool:
            with tqdm(total=len(lambda_range)*len(N_m_range)*len(hw_range)*len(hw_trap_range)*len(A_trap_range)) as pbar:
                for i, _ in enumerate(pool.imap( partial(run_wrapper_TBME, nuclei=nuclei, NN_interaction=NN_interaction, pless_params=pless_params,tasks_returned=results) , [{'srg_evolution':run_params['srg_evolution'],'lambda_range':[Lambda],'N_m_range':[N_m],'hw_range':[hw],'hw_trap':[hw_trap],'A_trap':[A_trap] } for Lambda,N_m,hw,hw_trap,A_trap in tasks ] )):
                    pbar.update()
        return results
    else:
        execution_status={}
        cluster = ipp.Cluster(profile="ijclab-ssh")

        rc = await cluster.start_and_connect(n=10)
        hosts_info = rc[:].apply_async(identify)
        hosts      = list(set([val['host'] for val in hosts_info]))

        my_hostname= socket.gethostname()
        #if my_hostname in hosts: hosts.remove(my_hostname)
        paths      = dict(rc[:].apply_async(locate))

        ssh = SSHClient()
        ssh.load_system_host_keys()
        for host in hosts:
            ssh.connect(hostname=host)
            scp = SCPClient(ssh.get_transport())
            for file_put in ['fit.py']:
                scp.put(file_put,paths[host])
            scp.close()

        balanced_view = rc.load_balanced_view()

        achieved_tasks=[]
        try:
            #achieved_tasks = balanced_view.map_async(hello_tasks,tqdm(tasks,total=len(lambda_range)*len(N_m_range)*len(hw_range)*len(hw_trap_range)*len(A_trap_range)))
            #achieved_tasks = balanced_view.map_async(make_parallel_wrapper,tqdm(tasks,total=len(lambda_range)*len(N_m_range)*len(hw_range)*len(hw_trap_range)*len(A_trap_range)))
            with tqdm(total=len(lambda_range)*len(N_m_range)*len(hw_range)*len(hw_trap_range)*len(A_trap_range)) as pbar:
                for i,achieved in enumerate(balanced_view.map_async( make_parallel_wrapper,tasks)):
                    achieved_tasks.append(achieved)
                    pbar.update()
        except Exception as e:
            print(e)
            pass
        achieved_tasks_dict  = {}
        achieved_tasks_dictrv= []
        for items in from_list(list(achieved_tasks)):
            key,val=items
            achieved_tasks_dict.setdefault(key['host'], []).append({'task':val,'status':key['status'],'directory':key['directory']})
            achieved_tasks_dictrv.append(list(reversed(items)))
        achieved_tasks_dictrv=dict(achieved_tasks_dictrv)

        if not len(list(tasks)) == len(achieved_tasks): print('Input/output tasks size mismatch')
        for task in tasks:
            if task in achieved_tasks_dictrv.keys():
                result=achieved_tasks_dictrv[task]
                if 'done' in result['status']:
                    print('Tasks with id={0} has been successfully achieved on {1}'.format(task,achieved_tasks_dictrv[task]['host']))
                else:
                    print('Tasks with id={0} has returned an/ error(s) {1}'.format(task,achieved_tasks_dictrv[task]['status']))
            else:
                print('Tasks with id={0} was not performed'.format(task))

        await cluster.stop_cluster()

        for host in hosts:
            ssh.connect(hostname=host)
            scp = SCPClient(ssh.get_transport())
            try:
                for directory in set([run_dir for val in achieved_tasks_dictrv.values() for run_dir in val['directory']]):
                    scp.get('/vol0/hupin/Runs{}'.format(directory),'./', recursive=True)
            except Exception as e:
                print(e)
                pass
            #except:
            #    pass
            scp.close()
    return

def run_ncsd(nuclei_in,interaction=[('TBME.int')],run_params={'N_m_range':[14],'parity_range':[1]},tasks_achieved=False,multi_host=False,mpi_task=1,mpi_node=1,mpi_task_per_node=1,cpu_per_task=1,gpu=0):
    import fileinput
    import tempfile
    import itertools     as it
    from  tqdm           import tqdm

    class CallableDict(dict):
        def shell_filling(self,i):
            m_state=0
            for n in range(0, int(i/2+1), 1):
                for l in range(0,i+1,1):
                    if 2*n+l != i:
                        pass
                    else:
                        for j in range(abs(2*l-1), 2*l+2, 2):
                            m_state+=j+1
            return m_state
        def n_in_shell(self,particle_in):
            n=0
            lowest_conf =[]
            particle_max=particle_in
            while particle_max > 0:
                lowest_conf.append(min(particle_max,self.shell_filling(n)))
                particle_max-= lowest_conf[n]
                n+=1
            if lowest_conf[-1]==self.shell_filling(n-1):
                result=n-1
            else:
                result=n-2
            return result
        def __getitem__(self, key):
            try:
                val = super().__getitem__(key)
                return val
            except:
                pattern_nucleon = re.compile("^([0-9]+)A$")
                result          = re.match(pattern_nucleon,key)
                return {'protons':0,'neutrons':int(result.group(1)),'N_shell_min':self.n_in_shell(int(result.group(1)))}

    NUCLEI_prop=CallableDict(
                {'2H' :{'protons':1,'neutrons':1,'N_shell_min':0},
                 '3H' :{'protons':1,'neutrons':2,'N_shell_min':0},
                 '3He':{'protons':2,'neutrons':1,'N_shell_min':0},
                 '4He':{'protons':2,'neutrons':2,'N_shell_min':0},
                 '5He':{'protons':2,'neutrons':3,'N_shell_min':1},
                 '6He':{'protons':2,'neutrons':4,'N_shell_min':1},
                 '4Li':{'protons':3,'neutrons':1,'N_shell_min':1},
                 '5Li':{'protons':3,'neutrons':2,'N_shell_min':1},
                 '6Li':{'protons':3,'neutrons':3,'N_shell_min':1},
                 '7Li':{'protons':3,'neutrons':4,'N_shell_min':1},
                 '8Li':{'protons':3,'neutrons':5,'N_shell_min':1},
                 '8Be':{'protons':4,'neutrons':4,'N_shell_min':1},
                 '9Be':{'protons':4,'neutrons':5,'N_shell_min':1},
                '10Be':{'protons':4,'neutrons':6,'N_shell_min':1}})

    if tasks_achieved:
        tasks_done=[]
    else:
        tasks_done=None

    def modify_ncsd_input_file(nuclei_in,N_max,parity,extension):
        import fileinput
        global f,line

        def shell_filling(i):
            m_state=0
            for n in range(0, int(i/2+1), 1):
                for l in range(0,i+1,1):
                    if 2*n+l != i:
                        pass
                    else:
                        for j in range(abs(2*l-1), 2*l+2, 2):
                            m_state+=j+1
            return m_state

        def shell_occupation(a):
            from math import floor
            global f,line
            def n_in_shell(i,particle_in,particle_type):
                if particle_type == 'nucleons':
                    filling=2
                else:
                    filling=1
                if particle_in==0: return 0
                n=0
                lowest_conf =[]
                particle_max=particle_in
                particle    =0
                while particle_max > 0:
                    lowest_conf.append(min(particle_max , filling*shell_filling(n)))
                    particle_max-= lowest_conf[n]
                    n+=1
                quanta=N_max
                parity=-(n-1)%2*particle_in%2
                if parity != parity_calc:
                    for j in range(n-1,-1,-1):
                        if lowest_conf[j] ==  filling*shell_filling(j):
                            break
                    if j==n-1:
                        n+=1
                        lowest_conf.append(min(1 , filling*shell_filling(n)))
                        lowest_conf[j]  -=1
                    else:
                        lowest_conf[j]  -=1
                        lowest_conf[j+1]+=1
                if i<n-1: return filling*shell_filling(i)
                for j in range(n-1,-1,-1):
                    for hole in range(1,lowest_conf[j]+1):
                        if quanta-(i-j) >= 0:
                            quanta-=(i-j)
                            particle+=1
                        else:
                            break
                    else:
                        continue
                    break
                return particle

            neutrons =NUCLEI_prop[nuclei_in]['neutrons']
            protons  =NUCLEI_prop[nuclei_in]['protons']
            nucleons =neutrons+protons
            new_line =line
            new_line =new_line.rstrip()
            new_line =new_line+"\n"
            kept_line=''
            for n_shell in range(0,run_params['interaction'][extension]['N_1_max']+1,1):
                line=next(f)
                neutrons_max=min(n_in_shell(n_shell ,neutrons,'neutrons'),                                  shell_filling(n_shell))
                protons_max =min(n_in_shell(n_shell ,protons ,'protons' ),                                  shell_filling(n_shell))
                nucleons_max=min(n_in_shell(n_shell ,nucleons,'nucleons' if protons != 0 else 'neutrons'),2*shell_filling(n_shell) if protons != 0 else shell_filling(n_shell))
                if 'N=' in line:
                    new_line    =new_line+re.sub(r"^( *[0-9]* *)[0-9]*( *[0-9]* *)[0-9]*( *[0-9]* *)[0-9]*(.*$)","\g<1>"+r'{:<2d}'.format(protons_max)+"\g<2>"+r'{:<2d}'.format(neutrons_max)+"\g<3>"+r'{:<2d}'.format(nucleons_max)+"\g<4>",line)
                else:
                    if kept_line == '': kept_line=line
                    new_line    =new_line+' 0 {0:d}    0 {1:d}    0 {2:d}'.format(protons_max,neutrons_max,nucleons_max)+'     ! N={}'.format(n_shell)+'\n'
            else:
                new_line=new_line+kept_line
            while 'N=' in line:
                line=next(f)
            print(new_line,end='')
            raise ValueError

        def lowest_conf_quanta(particle_in):
            if particle_in==0: return 0
            n=0
            lowest_conf =[]
            particle_max=particle_in
            while particle_max > 0:
                lowest_conf.append(min(particle_max , shell_filling(n)))
                particle_max-= lowest_conf[n]
                n+=1
            quanta=0
            for k,particle in enumerate(lowest_conf):
                quanta+=k*particle
            return quanta

        long_string_key="Major or subshell occupation restrictions for protons, neutrons and (p+n):"
        protons    = NUCLEI_prop[nuclei_in]['protons']
        neutrons   = NUCLEI_prop[nuclei_in]['neutrons']
        nucleons   = protons+neutrons
        Jz         = (nucleons)%2
        parity_calc= int(-(parity-1)/2)
        N_max_calc = N_max+parity_calc-nucleons%2+lowest_conf_quanta(neutrons)+lowest_conf_quanta(protons)

        if N_max+NUCLEI_prop[nuclei_in]['N_shell_min'] > run_params['interaction'][extension]['N_1_max']:
            raise RuntimeError("Interaction file not large for requested run")

        INPUT_ctrl={'Z, N, hbar*Omega' :                    ( lambda line : re.sub(r"^ +[0-9]+ +[0-9]+ +[.0-9]+(.*$)",r'  {:<2d} {:<2d} {:.2f}'.format(protons,neutrons,run_params['interaction'][extension]['hw'])+"\g<1>",line)),
                    'N_min,N_1max,N_2max' :                 ( lambda line : re.sub(r"^ +[0-9]+ +[0-9]+ +[0-9]+(.*$)",r'  {:<2d} {:<2d} {:<2d}'.format(run_params['interaction'][extension]['N_shell_min'],run_params['interaction'][extension]['N_1_max'],run_params['interaction'][extension]['N_12_max'])+"\g<1>",line)),
                    'rank of the hamiltonian' :             ( lambda line : re.sub(r"^[0-9]*(.*$)",r'{:<1d}'.format(run_params['interaction']['rank'])+"\g<1>",line)),
                    'TBMEfile' :                            ( lambda line : re.sub(r"^[a-zA-Z ]*$",'TBME.int',line)),
                    'fileout' :                             ( lambda line : re.sub(r"^[a-zA-Z ]*$",'NCSD.out',line)),
                    'Nhw  =  sum(2n+l),  Parity  (0  for  +,  1  for  -),  Total 2*Jz' :
                                                            ( lambda line : re.sub(r"^ +[0-9]+ +[0-9]+ +[0-9]+(.*$)",r'  {:<2d} {:<2d} {:<2d}'.format(N_max_calc,parity_calc,Jz)+"\g<1>",line)),
                    'ki,kf,nf (Initial, 1st final state, # of final states)' :
                                                            ( lambda line : re.sub(r"^ +[0-9]+ +[0-9]+ +[0-9]+(.*$)",r'  {:<2d} {:<2d} {:<2d}'.format(run_params.get('#1st',1),run_params.get('#last',1),run_params.get('#st',10))+"\g<1>",line)),
                    'Number of iterations required' :       ( lambda line : re.sub(r"^ +[0-9]+(.*$)",r'  {}'.format(run_params.get('iter',80))+"\g<1>",line)),
                    'irest = 4 to restart, = 0 not to.' :   ( lambda line : re.sub(r"^ +[0-9]+(.*$)",r'  {}'.format(run_params.get('restart',0))+"\g<1>",line)),
                    'nhw0, nhw_min (start IT), nhw_max=nhw, nhw_restart (nhw_begin)':
                                                            ( lambda line : re.sub(r"^ +[0-9]+ +[0-9]+ +[0-9]+(.*$)",r'  {:<2d} {:<2d} {:<2d}'.format(N_max_calc,N_max_calc+2,N_max_calc)+"\g<1>",line)),
                    'available memory per MPI process in GB in real(8)':
                                                            ( lambda line : re.sub(r"^ +[.0-9]+ +(.*$)",r'  {:.2f}'.format(run_params.get('mem-per-task',4.8))+"\g<1>",line)),
                    '3Nfile' :                              ( lambda line : re.sub(r"^[0-9a-zA-Z ]*$",'v3N.int' if run_params['interaction'][extension]['rank'] == 3 else '',line)),
                    long_string_key :                       shell_occupation}


        ctrl_pattern  =re.compile(r"! +(.*)|^ *([^\!\n]*) *$")
        str_convert   =lambda x : x or ''
        str_conadd    =lambda x,y : str_convert(x)+str_convert(y)

        f=fileinput.input(files=('mfdp_v1.dat'),inplace=True)
        line=next(f)
        while True:
            try:
                ctrl_key=re.search(ctrl_pattern,line)
                try:
                    new_line=INPUT_ctrl[str_conadd(ctrl_key.group(1),ctrl_key.group(2))](line)
                except KeyError:
                    new_line=line
                print(new_line,end='')
                line=next(f)
            except StopIteration:
                break
            except ValueError:
                continue
        f.close()
        return N_max_calc

    def rem(file,afix=''):
        if os.path.isfile(file+afix):
            os.remove(file+afix)

    def ren(file,newname):
       if os.path.isfile(file):
          os.rename(file,newname)

    def rem_auto_output_files(file_2_avoid):
        import glob
        files=glob.glob("lancz*.tmp") + glob.glob("lzegv.tmp")
        files=filter(lambda  elem:  not any(elem==item for item in file_2_avoid) ,files)
        for file in files:
            rem(file)

    def symlink_force(target, link_name):
        try:
            temp_link = link_name+'.new'
            rem(temp_link)
            os.symlink(target, temp_link)
            os.rename(temp_link, link_name)
        except OSError as e:
            raise

    def run_code(outputs={'mfdp.egv':'mfdp.egv_{:f}'.format(time.time())}):
        import time,shlex,os
        import socket
        from subprocess import Popen,PIPE,DEVNULL,STDOUT
        def joinit(iterable, delimiter):
            it = iter(iterable)
            yield next(it)
            for x in it:
                yield delimiter
                yield x

        files_here=glob.glob("*.dat*") + glob.glob("TBME*") + glob.glob("*.egv*") + glob.glob("*.out*")
        symlink_force('mfdp_v1.dat','mfdp.dat')

        if multi_host:
            command = shlex.split("env -i bash -c 'source /import/divers/intel/2017/parallel_studio_xe_2018.2.046/bin/psxevars.sh intel64 && env'")
            proc    = Popen(command, stdout = PIPE, stderr = DEVNULL)
            for line in proc.stdout:
                (key, _, value) = line.decode().partition("=")
                os.environ[key] = value
            proc.communicate()
            if any([ re.match(server,socket.getfqdn()) for server in IJCLab_server]):
                hostfile_name= 'hostfile.txt'
                if not os.path.isfile(hostfile_name):
                    #hostlist=['theorie1:100','theorie2:100']
                    hostlist=['ls-theo2:28','ls-theo3:28','ls-theo4:32']
                    #hostlist=['ls-theo2:28','ls-theo3:28','ls-theo4:32','ls-phynet01:16','theorie1:100','theorie2:100']
                    hostfile=open(hostfile_name,'w')
                    to_print=list(joinit(hostlist,'\n'))
                    hostfile.writelines(to_print)
                    hostfile.close()
                ntask  =sum([int(nproc) for nproc in dict(items.split(':') for items in hostlist).values()])
                options='-np {:<3d}'.format(ntask)+' -machine hostfile.txt'
            else:
                raise ValueError('Parallel env. not defined')
        else:
            if any([ re.match(server,socket.getfqdn()) for server in IJCLab_server]):
                mpi_dic={'ls-theo2'   : 28,
                         'ls-theo3'   : 28,
                         'ls-theo4'   : 32,
                         'ls-phynet01': 16,
                         'ls-phynet03': 48,
                         'theorie1'   : 64,
                         'theorie2'   : 64}
                ntask  =mpi_dic[socket.gethostname()]
                options='-np {:<3d}'.format(ntask)+' -host '+socket.gethostname()
            elif 'cca' in socket.getfqdn():
                options='--mpi=pmix -n {:<3d} --cpu-bind=cores'.format(mpi_task)
            elif os.getenv('MACHINE_NAME') is not None and 'idris' in os.getenv('MACHINE_NAME'):
                options='--nodes={0:<3d} --ntasks={1:<3d} --ntasks-per-node={2:<3d} --cpus-per-task={3:<3d} --gpus-per-task={4:<3d} --hint=nomultithread'.format(mpi_node,int(mpi_task_per_node*mpi_node),mpi_task_per_node,cpu_per_task,int(gpu/mpi_task_per_node))
            else:
                raise ValueError('Parallel env. not defined')

#        with open(os.devnull, 'w') as fp:
#            try:
#                subprocess.call(['mpirun,'-np {:<3d}'.format(ntask),'./'+exefile],shell=True,stdout=fp,stderr=fp)
#            except:
#                raise
        with open('runtime.log', 'w', buffering=1) as log:
            if any([ re.match(server,socket.getfqdn()) for server in IJCLab_server]):
                log.write('MPI job will start with {} tasks'.format(ntask))
                proc = Popen(shlex.split('mpirun '+options+' -genv I_MPI_PIN_DOMAIN core -genv I_MPI_DEBUG 4 ./'+exefile), stdin=DEVNULL, stdout=PIPE, stderr=STDOUT, encoding='utf-8')
            elif 'cca' in socket.getfqdn():
                log.write('MPI job will start with {} tasks'.format(mpi_task))
                proc = Popen(shlex.split('srun '+options+' ./'+exefile), stdin=DEVNULL, stdout=PIPE, stderr=STDOUT, encoding='utf-8')
            elif os.getenv('MACHINE_NAME') is not None and 'idris' in os.getenv('MACHINE_NAME'):
                log.write('MPI job will start with {0} tasks and {1} GPUs per nodes'.format(mpi_node*mpi_task_per_node,gpu)+'\n')
                proc = Popen(shlex.split('srun '+options+' ./'+exefile), stdin=DEVNULL, stdout=PIPE, stderr=STDOUT, encoding='utf-8')
            else:
                raise ValueError('Parallel env. not defined')
            while proc.poll() is None:
                out = proc.stdout.readline()
                log.write(out)
                log.flush()
        for output,new_name in outputs.items():
            ren(output,new_name)
        rem('mfdp','.dat')
        rem_auto_output_files([items for keys,items in outputs.items()]+files_here)

    def empty(*args):
        pass

    def make_list(unknown_type):
        return [unknown_type] if type(unknown_type) is not list else unknown_type

    mother_path     =os.getcwd()
    home            =os.environ['HOME']
    path_executable =mother_path
    if 'cca' in socket.getfqdn():
        home_depot      =os.path.join(home,'Codes/ncsm')
    elif os.getenv('MACHINE_NAME') is not None and 'idris' in os.getenv('MACHINE_NAME'):
        home_depot      =os.path.join(home,'Codes/ncsm/ncsd-gpu')
    else:
        home_depot      =os.path.join(home,'Documents/Codes/NCSM')

    if os.getenv('MACHINE_NAME') is not None and 'idris' in os.getenv('MACHINE_NAME'):
        exefile='ncsd-it-gpu.exe'
    else:
        exefile='ncsd_light.out'

    nuclei_range     =make_list(nuclei_in)
    interaction      =make_list(interaction)
    try:
        N_m_range    =run_params['N_m_range']
    except:
        raise RuntimeError('No Nmax in input')
    try:
        parity_range =run_params['parity_range']
    except:
        raise RuntimeError('No parity in input')

    base_param   = {'nuclei': nuclei_range,'N_m':N_m_range,'parity':parity_range}

    for interactions_files in tqdm(interaction):

        if type(interactions_files) is tuple:
            NN_file,NNN_file = interactions_files
            NNN_calc =True
        else:
            NN_file,NNN_file = interactions_files,None
            NNN_calc = False


        TBME_pattern =re.compile(r"TBME_NN(?:(.*)(?<!-srg[.0-9]{3,3})|(.*)(-srg)([.0-9]*)?)_([0-9]*[A-Z][a-z]?)_([0-9]*)\.([0-9]*)_?A?([0-9]*)?_?([0-9]*)?.int")
        NNNME_pattern=re.compile(r"v3trans_J3T3.int_(?:NN(.*)_3N(.*)(?<!-srg)|NN(.*)_3N(.*)(-srg)([.0-9]*)?)_([0-9]*[A-Z][a-z]?)_([0-9]*)\.([0-9]*)_?A?([0-9]*)?_?([0-9]*)?")
        NN_file_config  =re.match( TBME_pattern, NN_file)

        try:
            NN_pot       ='NN'+NN_file_config.group(1) if NN_file_config.group(1) else 'NN'+NN_file_config.group(2)
        except:
            raise RuntimeError('No NN_pot info in file')
        if NNN_calc:
            NNN_file_config =re.match(NNNME_pattern,NNN_file)
            try:
                NNN_pot      ='3N'+NNN_file_config.group(1) if NNN_file_config.group(2) else '3N'+NNN_file_config.group(4)
            except:
                raise RuntimeError('3N_pot info unreadable')
        try:
            srg_evolution=True if NN_file_config.group(3) else False
        except:
            srg_evolution=False
        try:
            Lambda       =float(NN_file_config.group(4)) if NN_file_config.group(4) else 0.0
        except:
            Lambda       =0.0
        try:
            N_mshell_min =NUCLEI_prop[NN_file_config.group(5)]['N_shell_min']
        except:
            raise RuntimeError('Cannot get N_mshell_min from file')
        try:
            N_m          =  int(NN_file_config.group(6))
        except:
            raise RuntimeError('No Nmax_1_max in file')
        try:
            hw           =float(NN_file_config.group(7))
        except:
            raise RuntimeError('No hw in file')
        try:
            HO_trap      =True if NN_file_config.group(8) else False
            if HO_trap:
                hw_trap=float(NN_file_config.group(9))
                A_trap =int(  NN_file_config.group(8))
            else:
                hw_trap=0.0
                A_trap =0
        except:
            HO_trap      =False
            hw_trap=0.0
            A_trap =0

        extension =NN_pot+'{}'.format('_'+NNN_pot if NNN_calc else '')

        run_params.setdefault('interaction',{}).update({'interaction' : run_params['interaction'].setdefault(extension,{}).update({'N_1_max':N_m})} if run_params['interaction'] else {extension:{'N_1_max':int(N_m)}})
        run_params['interaction'][extension]['N_shell_min']=int(N_mshell_min)
        run_params['interaction'][extension]['N_12_max']=int(N_m)
        run_params['interaction'][extension]['hw']=float(hw)
        run_params['interaction'][extension]['rank']=3 if NNN_calc else 2
        run_params['interaction'][extension]['2Nfile']=NN_file
        run_params['interaction'][extension]['3Nfile']=NNN_file if NNN_file else None

        working_dir =extension+('_SRG' if srg_evolution else '')
        working_path=os.path.join(mother_path,working_dir)

        if not os.path.isdir(working_dir):
            try:
                os.mkdir(working_dir)
            except:
                pass

        os.chdir(working_path)

        param, val   = zip(*base_param.items())

        for run_dict in [dict(zip(param, v)) for v in it.product(*val)]:
            nuclei  =run_dict.get('nuclei','2H')
            N_m     =run_dict.get('N_m',0.0)
            parity  =run_dict.get('parity' ,0.0)

            if HO_trap:
                if re.match(r'^A[0-9]+$',nuclei): continue
                if A_trap != NUCLEI_prop[nuclei]['neutrons']: continue

            if N_m+NUCLEI_prop[nuclei]['N_shell_min'] > run_params['interaction'][extension]['N_1_max']:
                continue

            extension_out=extension+'{0}{1}_{2}_{3:d}.{4:.0f}{5}{6}{7}'.format('-srg' if srg_evolution else '', '{:3.1f}'.format(Lambda) if srg_evolution else '',nuclei,N_m,hw,'_p' if parity == 1 else '_m','_A{:d}'.format(A_trap) if HO_trap else '','_{:.0f}'.format(hw_trap) if HO_trap else '')
            try:
                stat_info=os.stat('NCSD.out_'+extension_out)
                if stat_info.st_size > 2700:
                    continue
                else:
                    pass
            except:
                pass

            #with tempfile.TemporaryDirectory(dir=None if multi_host == False else {'dir':working_path}) as tmp_path:
            with tempfile.TemporaryDirectory(dir=working_path) as tmp_path:
                os.chdir(tmp_path)
                shutil.copy(os.path.join(home_depot+'/mfdp-template.dat'),'mfdp_v1.dat')
                #shutil.copy(os.path.join(path_executable,exefile),os.getcwd())
                shutil.copy(os.path.join(home_depot,exefile),os.getcwd())

                file_names = []

                N_tot = modify_ncsd_input_file(nuclei,N_m,parity,extension)

                symlink_force(os.path.join(working_path+'/'+run_params['interaction'][extension]['2Nfile']),'TBME.int')
                if NNN_calc: symlink_force(os.path.join(working_path+'/'+run_params['interaction'][extension]['3Nfile']),'v3N.int')

                try:
                    outputs={'NCSD.out':'NCSD.out_'+extension_out,
                             'mfdp_{}.egv'.format(N_tot):'mfdp.egv_'+extension_out,
                             'mfdp.egv'.format(N_tot):'mfdp.egv_'+extension_out,
                             'mfdp_v1.dat':'mfdp.dat_'+extension_out,
                             'mfd.log':'mfd.log_'+extension_out,
                             'runtime.log':'runtime.log'+extension_out}
                    run_code(outputs)
                except OSError:
                    print('Calculation was not performed for N_max={:d}'.format(N_m))
                    pass

                file_names.extend([files for key,files in outputs.items()])
                if tasks_achieved:
                    tasks_done.append({'tasks': extension_out, 'directory': working_path })

                for file_name in set(file_names):
                    try:
                        shutil.move(os.path.join(tmp_path, file_name), working_path)
                    except FileNotFoundError:
                        print('File named '+file_name+' was not generated')
                    except OSError as err:
                        if os.path.exists(os.path.join(working_path,file_name)):
                            shutil.move(os.path.join(working_path, file_name), os.path.join(working_path, file_name+'~'))
                            shutil.move(os.path.join(tmp_path, file_name), working_path)
                            print('Work already performed, file save in archive')
                            #break
                        else:
                            print(err)

        os.chdir(mother_path)
    return tasks_done

def run_wrapper_ncsd(interaction_in,nuclei_in=None,run_params_in=None,tasks_returned=None):
    import shlex
    import subprocess
    command = shlex.split("env -i bash -c 'source /import/divers/intel/2017/parallel_studio_xe_2018.2.046/bin/psxevars.sh intel64 && env'")
    proc    = subprocess.Popen(command, stdout = subprocess.PIPE, stderr=subprocess.DEVNULL)
    for line in proc.stdout:
        (key, _, value) = line.decode().partition("=")
        os.environ[key] = value
    proc.communicate()

    tasks_returned.append(run_ncsd(nuclei_in,interaction=interaction_in,run_params=run_params_in,tasks_achieved=True))
    return

async def ncsd_generator(nuclei=['4He'],interaction=[('TBME.int')],run_params={'N_m_range':[14],'parity_range':[1]},repository=None,multi_host=False):
    import os,socket
    from  tqdm           import tqdm
    from paramiko        import SSHClient
    from scp             import SCPClient
    from itertools       import product
    from joblib          import Parallel,delayed
    from multiprocessing import cpu_count,Process,Manager,Pool
    from functools       import partial

    def from_list(iterables):
        if isinstance(iterables,list):
            for elems in iterables:
                yield from from_list(elems)
        else:
            yield iterables

    def make_list(unknown_type):
        return [unknown_type] if type(unknown_type) is not list else unknown_type

    @ipp.interactive
    def identify():
        import os
        import socket
        return {'host': socket.gethostname(), 'pid': os.getpid()}

    @ipp.interactive
    def locate():
        import os
        import socket
        return (socket.gethostname(), os.getcwd())

    @ipp.interactive
    def hello_tasks(tasks_in):
        import os
        import socket
        return ({'host': socket.gethostname(),'status': 'done','directory':os.getcwd()}, tasks_in)

    def make_parallel_wrapper(tasks_in, repository='.'):
        import socket,os,shutil
        from functools       import partial
        from multiprocessing import cpu_count,Process,Manager,Pool,get_context

        if not isinstance(tasks_in, list):
            results      =[]
            if isinstance(tasks_in, tuple):
                filenames    =tasks_in
                working_dir  =tasks_in[0].replace('.int','').replace('TBME_','')
            else:
                filenames    =(tasks_in,)
                working_dir  =tasks_in.replace('.int','').replace('TBME_','')
            working_dir  =re.sub(r"_([0-9]*[A-Z][a-z]?)_([0-9]*)\.([0-9]*)_?A?([0-9]*)?_?([0-9]*)?","", working_dir)
            if not os.path.isdir(working_dir):
                try:
                    os.mkdir(working_dir)
                except:
                    pass
            for filename in filenames:
                if not os.path.exists(os.path.join(working_dir,filename)):
                        #os.symlink(os.path.join(repository,working_dir,filename),os.path.join(working_dir,filename))
                        shutil.copy(os.path.join(repository,working_dir,filename),os.path.join(working_dir,filename))
            run_wrapper_ncsd(tasks_in, nuclei_in=nuclei_range, run_params_in=run_params, tasks_returned=results)
            tasks_in=make_list(tasks_in)
        else:
            for interaction in tasks_in:
                if isinstance(interaction, tuple):
                    filenames    =tasks_in
                    working_dir  =interaction[0].replace('.int','').replace('TBME_','')
                else:
                    filenames    =(tasks_in,)
                    working_dir  =interaction.replace('.int','').replace('TBME_','')
                working_dir  =interaction.replace('.int','').replace('TBME_','')
                working_dir  =re.sub(r"_([0-9]*[A-Z][a-z]?)_([0-9]*)\.([0-9]*)_?A?([0-9]*)?_?([0-9]*)?","", working_dir)
                if not os.path.isdir(working_dir):
                    try:
                        os.mkdir(working_dir)
                    except:
                        pass
                for filename in filenames:
                    if not os.path.exists(os.path.join(working_dir,filename)):
                        #os.symlink(os.path.join(repository,working_dir,filename),os.path.join(working_dir,filename))
                        shutil.copy(os.path.join(repository,working_dir,filename),os.path.join(working_dir,filename))
            manager = Manager()
            results = manager.list()
            with get_context("spawn").Pool(processes=max(int(cpu_count()/20),2)) as pool:
                parallel_out = pool.imap( partial(run_wrapper_ncsd, nuclei_in=nuclei_range, run_params_in=run_params,tasks_returned=results) , tasks_in )

        directories=[]
        compute_task=[]
        for items in results:
            for subtasks in items:
                directories.append(subtasks['directory'].replace('/vol0/hupin/Runs',''))
                compute_task.append(subtasks['tasks'])
        if not results:
            return list(({'host': socket.gethostname(),'status': 'failed','directory':directories, 'compute_task':compute_task},b) for b in tasks_in)
        else:
            return list(({'host': socket.gethostname(),'status': 'done'  ,'directory':directories, 'compute_task':compute_task},b) for b in tasks_in)

    nuclei_range     =make_list(nuclei)
    interaction      =make_list(interaction)
    try:
        N_m_range    =run_params['N_m_range']
    except:
        raise RuntimeError('No Nmax in input')
    try:
        parity_range =run_params['parity_range']
    except:
        raise RuntimeError('No parity in input')

    tasks=interaction

    if not multi_host:
# This is specific to hyperthreaded env.
        n_cores = int(cpu_count()/2)
        manager = Manager()
        results = manager.list()

        processes = []
        with Pool(processes=max(int(n_cores/16),2)) as pool:
            with tqdm(total=len(tasks)) as pbar:
                for i, _ in enumerate(pool.imap( partial(run_wrapper_ncsd, nuclei_in=nuclei_range, run_params_in=run_params,tasks_returned=results) , tasks )):
                    pbar.update()
        return results
    else:
        execution_status={}
        cluster = ipp.Cluster(profile="ijclab-ssh")

        rc = await cluster.start_and_connect(n=4)
        hosts_info = rc[:].apply_async(identify)
        hosts      = list(set([val['host'] for val in hosts_info]))

        my_hostname= socket.gethostname()
        #if my_hostname in hosts: hosts.remove(my_hostname)
        paths      = dict(rc[:].apply_async(locate))

        ssh = SSHClient()
        ssh.load_system_host_keys()
        for host in hosts:
            ssh.connect(hostname=host)
            scp = SCPClient(ssh.get_transport())
            for file_put in ['fit.py']:
                scp.put(file_put,paths[host])
            scp.close()

        balanced_view = rc.load_balanced_view()

        achieved_tasks=[]
        make_parallel_wrapper_partial=partial(make_parallel_wrapper,repository=repository)
        try:
            make_parallel_wrapper_partial=partial(make_parallel_wrapper,repository=repository)
            with tqdm(total=len(tasks)) as pbar:
                for i,achieved in enumerate(balanced_view.map_async( make_parallel_wrapper_partial, tasks)):
                    achieved_tasks.append(achieved)
                    pbar.update()
        except Exception as e:
            print(e)
            pass
        achieved_tasks_dict  = {}
        achieved_tasks_dictrv= []
        for items in from_list(list(achieved_tasks)):
            key,val=items
            achieved_tasks_dict.setdefault(key['host'], []).append({'task':val,'status':key['status'],'directory':key['directory']})
            achieved_tasks_dictrv.append(list(reversed(items)))
        achieved_tasks_dictrv=dict(achieved_tasks_dictrv)

        if not len(list(tasks)) == len(achieved_tasks): print('Input/output tasks size mismatch')
        for task in tasks:
            if task in achieved_tasks_dictrv.keys():
                result=achieved_tasks_dictrv[task]
                if 'done' in result['status']:
                    print('Tasks with id={0} has been successfully achieved on {1}'.format(task,achieved_tasks_dictrv[task]['host']))
                else:
                    print('Tasks with id={0} has returned an/ error(s) {1}'.format(task,achieved_tasks_dictrv[task]['status']))
            else:
                print('Tasks with id={0} was not performed'.format(task))

        await cluster.stop_cluster()

        for host in hosts:
            ssh.connect(hostname=host)
            scp = SCPClient(ssh.get_transport())
            try:
                for directory in set([run_dir for val in achieved_tasks_dictrv.values() for run_dir in val['directory']]):
                    scp.get('/vol0/hupin/Runs{}'.format(directory),'./', recursive=True)
            except Exception as e:
                print(e)
                pass
            #except:
            #    pass
            scp.close()
    return

def help():
    """Prints information about how to use this script to run manyeffv3b."""

    print("""
======================== HELP: manyeffv3b Interface ========================

This script helps organize input files and run the manyeffv3b nuclear
structure code smoothly. Output files are stored in directories named
after the interaction, SRG parameters, and the nucleus. For example:

    NNn3lo_SRG_4He/

Eigenvalues are stored in files beginning with 'many_'. Example filename:

    many_NNn3lo_srg2.1_srgf1_4He-T12_1_N20_hw20_J0_T0_theta0.300_fit_Nf8

--------------------------------------------------------------------------
Main Callable Function:
--------------------------------------------------------------------------
run_manyeff(nuclei,
            NN_interaction=['n3lo'],
            NNN_interaction={'regulator': 'lnl', 'cutoff': [650,500],
                             'c_D': [0.7], 'c_E': [-0.06], 'E8': [0.0]},
            run_params={'srg_evolution': True, 'NNN_srg': False,
                        'NNNN_srg': False, 'lambda_range': [2.0],
                        'srg_f': 1.0, 'N_m_range': [14], 'hw_range': [20],
                        'theta_range': [0.0], 'Nf': -1},
            pless_params=None,
            output_stream=0,
            J=None,
            T=None,
            testing=False)

Description:
    Executes manyeffv3b runs with specified nuclear force settings, model
    space parameters, and output preferences.

Parameters:
    nuclei (str)           : The nucleus to calculate (e.g., '4He').
    NN_interaction (list)  : List of NN interactions to use (e.g., ['n3lo'], ['mtv']).
    NNN_interaction (dict) : Dictionary of 3N interaction settings or None to disable.
    run_params (dict)      : Dictionary of model space and SRG parameters:
        - srg_evolution    : Apply SRG evolution? (True/False)
        - NNN_srg          : Apply SRG to A=3 space?
        - NNNN_srg         : Apply SRG to A=4 space?
        - lambda_range     : List of SRG λ values (e.g., [2.0, 2.2])
        - srg_f            : SRG generator scaling factor (deprecated)
        - N_m_range        : List of N_max values (e.g., [18, 20])
        - hw_range         : Harmonic oscillator frequencies (MeV)
        - theta_range      : List of complex scaling θ values (radians)
        - Nf               : Number of basis functions in fitting approach.
                             Nf = 0: no fitting, <0: use n3lo_q_c modified regulator

    pless_params (dict)    : Optional dictionary for pless scheme (or None).
    output_stream (int)    : 0=silent, 1=stdout, 2=log file ('manyeff.log...')
    J (int)                : Total angular momentum (optional).
    T (int)                : Isospin (optional).
    testing (bool)         : Run in isolated 'Testing' subdirectory? (True/False)

--------------------------------------------------------------------------
Example Run:
--------------------------------------------------------------------------
Run a calculation for ⁴H with mtv interaction and Nf=12:

    python3 -c 'import fit; fit.run_manyeff(
        nuclei="4H",
        NN_interaction=["mtv"],
        NNN_interaction=None,
        run_params={
            "srg_evolution": True,
            "N_m_range": [21],
            "hw_range": [60],
            "theta_range": [0.58],
            "additional_opt": "4He_th",
            "lambda_range": [2.4],
            "NNN_srg": True,
            "srg_f": 1,
            "Nf": 12
        },
        J=0,
        T=1,
        output_stream=2,
        testing=False
    )'

==========================================================================
""")

if __name__ == "__main__":
    from inspect import getmembers, isfunction
    import fit
    import json
    mod_list=dict(getmembers(fit, isfunction))

    if len(sys.argv) < 2:
       raise RuntimeError('Missing information in inputs')

    command,arguments=sys.argv[1],json.loads(sys.argv[2])

    if command == 'run_ncsd':
        files = glob.glob("*/*TBME*")
        files = [os.path.basename(path) for path in files]
        arguments['interaction']=files
    mod_list[command](**arguments)
