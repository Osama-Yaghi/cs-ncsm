import numpy as np
import pandas as pd
import seaborn as sns
from math import ceil, sqrt
import math
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm,Normalize
from scipy.ndimage.filters import gaussian_filter
from scipy.optimize import curve_fit
from scipy.optimize import least_squares
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes, mark_inset
from matplotlib.ticker import AutoMinorLocator
import os
import sys
import pickle
import inspect
import textwrap
data=pd.DataFrame(columns=['E_r','E_i','N_max','hw','J','T','parity','srg','lambda','NNN','3N_induced','theta','Nf','srg_f'])
error_style={'lambda':'band','N_max':'convergence','hw':'std','Nf':'convergence','theta':'convergence'}
thresholds_exp={"4He":[-8.482,-7.718,-4.448,-2.224,0]
            ,"4H":[-7.718,-2.224,0]
            ,"3H":[-2.224,0]
            ,"4Li":[-2.224,0]
            ,"4N":[0]
           }
thresholds_NN={"4He":[-8.54205,-7.87064,-4.42696,-2.21348,0]
            ,"4H":[-8.54205,-2.21348,0]
            ,"3H":[-2.21348,0]
            ,"4Li":[-7.87064,-2.21348,0]
            ,"4N":[0]
            ,"2H":[0]
           }
thresholds_3Nind={"4He":[-8.04410,-7.38583,-4.42696,-2.21348,0]
            ,"4H":[-8.04410,-2.21348,0]
            ,"3H":[-2.21348,0]
            ,"4Li":[-7.38583,-2.21348,0]
            ,"4N":[0]
           }
thresholds_NN_Nf4={"4He":[-8.049,-7.46,-3.874,-1.937,0]
            ,"4H":[-7.46,-1.937,0]
            ,"3H":[-1.937,0]
            ,"4Li":[-1.937,0]
            ,"4N":[0]
           }
#thresholds_test={"4He":[-10.2,-9.67,-3.874,-1.937,0]}#1.05
#thresholds_test={"4He":[-11.7216,-11.0969,-6.674,-3.337,0]}#1.08
#thresholds_test={"4He":[-12.716,-12.084,-7.450,-3.725,0]}#1.1
thresholds_test={"4He":[-14.259,-13.615,-8.6736,-4.3368,0]}#1.13
#thresholds_test={"4He":[-15.32,-14.67,-9.528,-4.764,0]}#1.15


#Nf12
#thresholds_test={"4He":[-10.8826,-10.184,-7.450,-3.725,0]}#1.05
#thresholds_test={"4He":[-12.4186,-11.706,-7.450,-3.725,0]}#1.08
#thresholds_test={"4He":[-12.612,-11.8956,-7.450,-3.725,0]}#1.08 srg0.2
#thresholds_test={"4He":[-13.484,-12.762,-7.450,-3.725,0]}#1.1
#thresholds_test={"4He":[-16.4763,-15.7295,-10.578,-5.289,0]}#1.15
thresholds_test={"4He":[-25.625,-24.832,-10.578,-5.289,0]}#1.15


symbol={"hw":r"$\hbar\omega$","lambda":r"$\lambda$","theta":r"$\theta$","Nf":r"$N_f$","N_max":r"$\mathrm{N_{max}}$","srg_f":r"$f$","srg":r"$\mathrm{srg}$"}
units={"hw":r"$\mathrm{MeV}$","theta":r"$\mathrm{rad}$"}

curr_nucleus="4He"
plt.rc('text', usetex=True)
plt.rc('font', **{'family':'sans-serif','sans-serif':['Helvetica']})
plt.rcParams['text.latex.preamble']=r"\usepackage{amsmath,amssymb,siunitx,gensymb}"
A=0
curr_nucleus=""
#Reading files
def readfiles(nucleus,name1='',replace=False):
    """
    Reads all output files containing eigenvalue data and stores them 
    in a pandas DataFrame named <nucleus>.

    Parameters:
    - nucleus (str): The name of the nucleus for which data is being read.
    - name1 (str): The prefix used to identify which files to read (e.g., 'many_').
    - replace (bool, optional): If True, replaces existing data in the DataFrame. Defaults to False.

    Returns:
    - None
    """
    global data,curr_nucleus
    curr_nucleus=nucleus
    files=os.listdir(os.getcwd())
    matching_files=[file for file in files if (file.startswith(name1) and nucleus in file)]
    for file in matching_files:
        srg=0
        srglambda=0
        if('_srg' in file):
            extracted=file.split('srg')[1].split('_')[0]
            if(extracted[0]=='f'):
                srg=0
                srglambda=0
            else:
                srg=1
                srglambda=float(extracted)
        if('theta' in file):
            extracted=file.split('theta')[1].split('_')[0]
            theta=float(extracted)
        if('fit_Nf' in file or '_Nf' in file):
            extracted=file.split('Nf')[1].split('_')[0]
            Nf=float(extracted)
        else:
            Nf=0
        if('_srgf' in file):
            extracted=file.split('_srgf')[1].split('_')[0]
            srg_f=float(extracted)
        else:
            srg_f=1
        if('_NNN_' in file):
            NNN=1
        else:
            NNN=0
        if('_A3srg' in file):
            A3srg=1
        else:
            A3srg=0
        if('_N' in file):
            #extracted=file.lstrip(name1).split('_N')[1].split('_')[0]
            extracted=file.split('_N')[1].split('_')[0]
            if(extracted[0]=='N'):
                extracted=file.split('_N')[2].split('_')[0]
            N_max=int(extracted)
        if('_T' in file):
            extracted=file.split('_T')[1].split('_')[0]
            T=int(extracted)
        if('_J' in file):
            extracted=file.split('_J')[1].split('_')[0]
            J=int(extracted)
        if('_hw' in file):
            extracted=file.split('_hw')[1].split('_')[0]
            hw=int(extracted)
        if(N_max%2==0):
            parity=1
        else:
            parity=-1
        mask=check_params({'3N_induced':A3srg,'NNN':NNN,'srg':srg,'lambda':srglambda,'J':J,'T':T,'parity':parity,'hw':hw,'N_max':N_max,'theta':theta,'Nf':Nf,'srg_f':srg_f})
        mask = mask.reset_index(drop=True)
        mask=mask.reindex(data.index,fill_value=False)
        if(not data[mask].empty):
            if(replace):
                data.drop(data[mask].index,inplace=True)
                print("Updating data from file: NNN={8} 3N-induced={9} J={0} T={1} parity={2} hw={3} N_max={4} srg={5} lambda={6} theta={7} Nf={10} srg_f={11}".format(J,T,'+'
                      if(parity==1)else'-',hw,N_max,'On'if(srg==1) else 'Off',srglambda
                      if (srg==1) else '',theta,NNN,A3srg,Nf,srg_f))

            else:
                continue
        else:
            print("reading data from file: NNN={8} 3N-induced={9} J={0} T={1} parity={2} hw={3} N_max={4} srg={5} lambda={6} theta={7} Nf={10} srg_f={11}".format(J,T,'+'
                      if(parity==1)else'-',hw,N_max,'On'if(srg==1) else 'Off',srglambda
                      if (srg==1) else '',theta,NNN,A3srg,Nf,srg_f))

        N=N_max
        with open(file) as f1:
            for line in f1:
                if(not( "Spectrum" in line)):
                    continue
                N=int(line.split('Nmax=')[1])

                #if(not( r"%" in line)):
                #    continue
                E_r_temp=[]
                E_i_temp=[]
                line=next(f1)
                while(not("End of spectrum" in line) and not("Eigen" in line)):
                #while(r"%" in line):
                    s=line.split()
                    E_r_temp.append(float(s[0]))
                    E_i_temp.append(float(s[1]))
                    line=next(f1)
                length=len(E_r_temp)
                theta_temp=[theta for i in range(length)]
                Nf_temp=[Nf for i in range(length)]
                srg_f_temp=[srg_f for i in range(length)]
                N_temp=[N for i in range(length)]
                srg_temp=[srg for i in range(length)]
                NNN_temp=[NNN for i in range(length)]
                A3srg_temp=[A3srg for i in range(length)]
                srglam_temp=[srglambda for i in range(length)]
                J_temp=[J for i in range(length)]
                T_temp=[T for i in range(length)]
                hw_temp=[hw for i in range(length)]
                parity_temp=[parity for i in range(length)]
                new_data={"E_r":E_r_temp,"E_i":E_i_temp,"N_max":N_temp,"hw":hw_temp,'srg':srg_temp,'J':J_temp,'T':T_temp,'parity':parity_temp,'lambda':srglam_temp,'NNN':NNN_temp,'3N_induced':A3srg_temp,'theta':theta_temp,'Nf':Nf_temp,'srg_f':srg_f_temp}
                new_pd=pd.DataFrame(new_data)
                data=pd.concat([data,new_pd]).drop_duplicates()
                #if(N>30):
                #if(N==N_max):
                #    N=N-2
                #    #N=N
                #N=N-2
    data.to_csv(nucleus,sep=',',index=False)

#Data manipulation and UI
def help():
    print("""
==================== HELP: Spectrum Reader =====================

This script reads and visualizes data from the manyeffv3b.f code.
The parameters of the visualization need to be modified directly in the code.
Some debugging and modifications might be necessary to use this script fully.

---------------------------------------------------------------
Examples of Usage:
---------------------------------------------------------------
  1. Reading files for a specific nucleus:
     python3 -c "import spectrum_reader; spectrum_reader.readfiles('4He', 'many_')"

  2. Loading data and plotting a general graph:
     python3 -c "import spectrum_reader; spectrum_reader.load('4He'); spectrum_reader.general_plot('N_max', 0, 'hw')"

---------------------------------------------------------------
Data Manipulation Functions:
---------------------------------------------------------------
""")

    # List and display docstrings for data manipulation functions
    data_functions = [readfiles, load, update, remove_data,write_resonance,write_thresholds,export_resonances]
    for func in data_functions:
        doc=inspect.getdoc(func)
        if doc:
            lines=doc.split('\n',1)
            signature=lines[0]
            description=textwrap.indent(textwrap.fill(lines[1],width=75),prefix='\t\t') if len(lines)>1 else ''
            print(f"  {func.__name__:<20}: {signature}")
            if description:
                print(description)
        print()

    print("---------------------------------------------------------------")
    print("Visualization Functions:")
    print("---------------------------------------------------------------")

    # List and display docstrings for visualization functions
    plotting_functions = [plot_complex_plane, general_plot, track_state_rotation, plot_rotation_hw, resonance_hwconv]
    for func in plotting_functions:
        doc=inspect.getdoc(func)
        if doc:
            lines=doc.split('\n',1)
            signature=lines[0]
            description=textwrap.indent(textwrap.fill(lines[1],width=75),prefix='\t\t') if len(lines)>1 else ''
            print(f"  {func.__name__:<20}: {signature}")
            if description:
                print(description)
        print()

    print("---------------------------------------------------------------")
    print("End of Help")
    print("---------------------------------------------------------------")



def load(nucleus):
    """
    Loads the DataFrame stored in a file corresponding to <nucleus> 
    in preparation for further tasks. Should be called before plotting.

    Parameters:
    - nucleus (str): The name of the nucleus whose data is being loaded.

    Returns:
    - None
    """
    global data,curr_nucleus,A
    data=pd.read_csv(nucleus,sep=',')
    curr_nucleus=nucleus
    A= int(curr_nucleus[0])
    if(nucleus=='2H'):
        data['T']=0
        data.to_csv(nucleus,sep=',',index=False)
def update(nucleus,name1=''):
    """
    Updates the data in the <nucleus> DataFrame by reading new data 
    from output files in the current directory that start with the prefix <name1>.

    Parameters:
    - nucleus (str): The name of the nucleus whose data is being updated.
    - name1 (str): The prefix used to identify which files to update from.

    Returns:
    - None
    """
    global data
    data=pd.read_csv(nucleus,sep=',')
    rep=input("do you want to replace existing data? yes/no")
    replace= (rep=="yes")
    readfiles(nucleus,name1,replace)
def remove_data(nucleus):
    """
    Removes data from the <nucleus> DataFrame. Prompts the user to specify 
    which data to delete interactively.

    Parameters:
    - nucleus (str): The name of the nucleus from which data will be removed.

    Returns:
    - None
    """
    global data
    data=pd.read_csv(nucleus,sep=',')
    curr_nucleus=nucleus
    print("what mode removal do you want?\n 1- a specific parameter combination\n 2-all data with certain parameters\n 3-entire column")
    typ=input("")
    typ=int(typ)
    print(typ)
    if(typ==1):
        params=ask_data_params(['J','T','parity','3N_induced','NNN','hw','srg','lambda','theta','Nf','srg_f'])
        data=data[~(check_params(params))]
    elif(typ==2):
        print("possible parameters are: J, T , parity, hw , N_max, srg, lambda,theta, Nf")
        s=input("which of these you want to specify?")
        params=ask_data_params(s.split())
        data=data[~(check_params(params))]
    elif(typ==3):
        print( 'avaliable columns:',data.columns)
        s=input("which columns do you want to delete? (input the number) ")
        for s1 in s.split():
            i1=int(s1)
            print("are you sure you want to delete?y/n ",data.columns[i1])
            s3=input("")
            if(s3=="y"):
                data=data.drop(data.columns[i1],axis=1)
        print( 'remaining columns:',data.columns)
    else:
        data=data
    s=input("do you still want to delete?y/n")
    if(s=='y'):
        data.to_csv(nucleus,sep=',',index=False)

#Assisting functions
def check_existing_params(nucleus,simple=True):
    global data
    def print_comp_recursive(df,level,s):
        if(level==len(df.columns)):
            print(s)
            return
        col=df.columns[level]
        if(col=="N_max" or (col=="theta" and simple)):
            print_comp_recursive(df,level+1,s)
            return
        vals=df[col].unique()
        for val in vals:
            sub_df=df[df[col]==val]
            print_comp_recursive(sub_df,level+1,s+" "*3+ col+"="+str(val))
    data=pd.read_csv(nucleus,sep=',')
    curr_nucleus=nucleus
    print("enter the parameters you want to check along with their values(i.e 'T 1'")
    print("possible parameters:",data.columns[2:])
    s=input()
    params=s.split()[0::2]
    values=[float(r) for r in s.split()[1::2]]
    mask=data['J']>-5
    for i,param in enumerate(params):
        mask=mask& (data[param]==values[i])
    subdata=data[mask]
    print_comp_recursive(subdata,2,"")
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
def create_name_extension(params):
    s=""
    for param,value in params.items():
        s=s+'_'+param+'{}'.format(value)
    return s
def format_tick_labels(f):
    res=["${:.3f}".format(x).rstrip('0')+"$" for x in f]
    return res
def ask_data_params(requested_params):
    subdata=data
    params={}
    for param in requested_params:
        if(subdata[param].nunique()==0):
            sys.exit("¯\_(¬_¬)_/¯    I don't find this parameter in data:"+param+"\n did you really load the data?")
        if(subdata[param].nunique()==1):
            value=subdata[param].unique()[0]
            print(param+" is set to {}".format(value)+ " (only available value)")
        else:
            print("available values for "+ param+" are:",np.sort(subdata[param].unique()))
            value=input(" choose a value:")
        params[param]=float(value)
        if(param=='lambda'):
            subdata=subdata[(subdata[param]==params[param])|(subdata[param]==0)]
        else:
            subdata=subdata[subdata[param]==params[param]]
        if(subdata.empty):
            sys.exit("¯\_(¬_¬)_/¯    C'est pas possible!")
    print(params)
    return params

def request_input(params):
    result={}
    for param in params:
        value=input("choose a value for {}".format(param))
        result[param]=value
    return result

def save_figure(figure,dir,name):
    os.makedirs(dir, exist_ok=True)
    print("saving figure: "+name)
    figure.savefig(dir+name+".pdf",format='pdf')
    figure.savefig(dir+name+".svg")
    figure.savefig(dir+name+".eps", bbox_inches='tight')
    figure.savefig(dir+name+".png")

#visualization
def fig_size(nx,ny,scale=1):
    latex_width=438.17227*scale
    plot_width=1.0/72.27 *latex_width
    aesthetic_ratio=(5.0**0.5-1)/2*1.1
    plot_height=plot_width*aesthetic_ratio*nx/ny
    figsize=(plot_width,plot_height)
    return figsize

def ask_for_extra_param(extra_param,subdata):
    n_ex_param=1
    if(extra_param!=''):
        ex_param=subdata[extra_param].unique()
        print("available "+extra_param+"  values:",ex_param)
        s=input("choose the desired values ( all in one line ), or -1 for all")
        ex_param_in=[float(s1) for s1 in s.split()]
        if(ex_param_in[0]==-1):
            ex_param_in=ex_param
        ex_param_in=np.sort(ex_param_in)
        n_ex_param=len(ex_param_in)
        print(" working with: ",n_ex_param,ex_param_in)
    else:
        ex_param_in=[]
        n_ex_param=0
    return ex_param_in,n_ex_param
def set_threshold(params):
    if(params['Nf']==4):
        thresholds=thresholds_NN_Nf4
    elif(params['3N_induced']==1):
        thresholds=thresholds_3Nind
    else:
        thresholds=thresholds_NN
    #thresholds=thresholds_test
    return thresholds

def write_thresholds():
    """ Add the parameters of a set of thresholds. They are stored in a thresholds.dat file which can be read when using Export_resonance.

    Parameters:
    - None
    Returns:
    - None
    """
    columns=['E_r','N_max','hw','srg','lambda','NNN','3N_induced','theta','Nf','srg_f','Er_total']

    threshold_data=pd.DataFrame(columns=columns)
    if(os.path.exists('thresholds.dat')):
        threshold_data=pd.read_csv('thresholds.dat',sep=',')
    print("parameters are :",columns[1:-1])
    answer=input("Enter the paramters of the threshold that you want to add")
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
    print("existing thresholds: \n",subdata['E_r'])

    answer=input("do you want to add a new set of thresholds or modify existing ones?(1/2)")
    if(answer=='2'):
        threshold_data=threshold_data[~mask]
    answer=input("enter thresholds in a line:")
    thresholds=np.array([float(item) for item in answer.split()])
    thresholds=np.sort(thresholds)
    answer=input("enter thresholds errors in a line (starting from the first threshold):")
    thresholds_errors=np.array([float(item) for item in answer.split()])
    for threshold,threshold_err in zip(thresholds,thresholds_errors):
        new_data={key:[value] for key,value in params.items()}
        new_data["E_r"]=[threshold]
        new_data["Er_total"]=[threshold_err]
        new_pd=pd.DataFrame(new_data)
        threshold_data=pd.concat([threshold_data,new_pd]).drop_duplicates()

    mask=check_params(params,threshold_data)
    subdata=threshold_data[mask]
    print("new threshold list:\n",subdata)
    answer=input("confirm the edits? y/n")
    if(answer=='y'):
        threshold_data.to_csv('thresholds.dat',sep=',',index=False)

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

def write_resonance():
    """ Add the parameters of an identified resonance. It is stored in resonances.dat file which can be read when using Export_resonance.

    Parameters:
    - None
    Returns:
    - None
    """
    columns=['J','T','parity','E_r','E_i','N_max','hw','srg','lambda','NNN','3N_induced','theta','Nf','srg_f','Er_hw','Er_lambda','Er_Nmax','Er_Nf','Er_theta','Er_total','Ei_hw','Ei_lambda','Ei_Nmax','Ei_Nf','Ei_theta','Ei_total']
    resonance_data=pd.DataFrame(columns=columns)
    if(os.path.exists('resonances.dat')):
        resonance_data=pd.read_csv('resonances.dat',sep=',')
    params_i=request_input(['J','T','parity'])
    params_i={key:float(value) for key,value in params_i.items()}
    mask=check_params(params_i,resonance_data)
    subdata=resonance_data[mask]
    print("existing resonances: \n",subdata)
    answer=input("do you want to add a new resonance or modify an existing one?(-1/ line-to-modify)")
    print("possible resonance parameters to enter:",columns)
    if (int(answer)==-1):
        answer=input("Enter the paramters of the resonance that you want to add")
        answer=answer.split()
        input_error= [name for name in answer if name not in columns]
        if(len(input_error)>0):
            sys.exit("error in entered parameters")
        values=input("Enter values:")
        values=values.split()
        values=[[float(value)] for value in values]
        new_data={**params_i,**dict(zip(answer,values))}
        new_pd=pd.DataFrame(new_data)
        resonance_data=pd.concat([resonance_data,new_pd]).drop_duplicates()
        print(resonance_data)
    else:
        row=int(answer)
        answer=input("which of the parameters you want to add/modify?")
        answer=answer.split()
        input_error= [name for name in answer if name not in columns]
        if(len(input_error)>0):
            sys.exit("error in entered parameters")
        values=input("Enter values:")
        values=values.split()
        values=[float(value) for value in values]
        subset_row=resonance_data.index[mask][row]
        for param,value in zip(answer,values):
            resonance_data.loc[subset_row,param]=value
    #resonance_data['Er_total']=np.sqrt(resonance_data['Er_Nmax']**2+resonance_data['Er_lambda']**2)
    #resonance_data['Ei_total']=np.sqrt(resonance_data['Ei_Nmax']**2+resonance_data['Ei_lambda']**2)
    mask=check_params(params_i,resonance_data)
    subdata=resonance_data[mask]
    print("new resonance list:\n",subdata)
    answer=input("confirm the edits? y/n")
    if(answer=='y'):
        resonance_data.to_csv('resonances.dat',sep=',',index=False)


def export_resonances(curr_nucleus='4He',name_ext='',band='width',threshold_in=None):
    """
    Exports the identified resonances into a .pkl dictionary and provides the input for
    the resonance_search.py program.

    Parameters:
    - curr_nucleus (str): The name of the nucleus whose data is being loaded.
    -name_ext(str): the name extension of the output file.
    -band(str): determines the error used ( 'Er_Nmax,Er_theta,Er_Nf, width, etc.')
    Returns:
    - None
    """
    columns=['J','T','parity','E_r','E_i','N_max','hw','srg','lambda','NNN','3N_induced','theta','Nf','srg_f','Er_hw','Er_lambda','Er_Nmax','Er_Nf','Er_total','Ei_hw','Ei_lambda','Ei_Nmax','Ei_Nf','Ei_total']
    if(os.path.exists('resonances.dat')):
        resonance_data=pd.read_csv('resonances.dat',sep=',')
    else:
        print("No resonances file found")
        return
    print("All available resonances:\n",resonance_data)
    answer=input(" do you want a filter on the exported resonances?y/n")
    if(answer=='y'):
        answer=input("which of the parameters you want to filter?")
        answer=answer.split()
        input_error= [name for name in answer if name not in columns]
        if(len(input_error)>0):
            sys.exit("error in entered parameters")
        values=input("Enter values:")
        values=values.split()
        values=[float(value) for value in values]
        params=dict(zip(answer,values))
        mask=check_params(params,resonance_data)
        subdata=resonance_data[mask]
        input_error=input("extra filter?(n/upper/lower")
        if(input_error=='upper'):
            param1=input('enter parameter')
            value=float(input('enter value'))
            subdata=subdata[subdata[param1]<value]
        if(input_error=='lower'):
            param1=input('enter parameter')
            value=float(input('enter value'))
            subdata=subdata[subdata[param1]>value]
    else:
        subdata=resonance_data
    print(subdata)
    # create the resonance dictionary
    thresholds,thresholds_errors=read_thresholds()
    first_threshold=thresholds[0]
    #if(not(threshold_in is None)):
    #    first_threshold=threshold_in + thresholds[-1]
        #first_threshold=thresholds[0]
    answer=input("continue with the previous thresholds?y/n")
    if(answer=='n'):
        sys.exit()
    resonance_dict={}
#    threshold_dict={}
    for _,row in subdata.iterrows():
        key= (int(row['J']),int(row['parity']),int(row['T']))
        print("adding a resonance in ", key)
        if(band=='width'):
            value= (row['E_r']-first_threshold,2*row['E_i'],0)
        elif('Er' in band):
            value= (row['E_r']-first_threshold,2*row['E_i'],row[band])
        elif(band=='error'):
            value= (row['E_r']-first_threshold,2*row['E_i'],row['Er_Nmax'])
        elif(band=='srg'):
            value= (row['E_r']-first_threshold,row['Er_lambda'])
        else:
            value=(row['E_r']-first_threshold,0,0)
        if( np.isnan(value[1])):
            value=(value[0],0,0)
        if(key in resonance_dict):
               resonance_dict[key]['Resonances'].append(value)
        else:
            resonance_dict[key]={'Resonances':[value]}

 #   threshold_dict['Thresholds']=zip(thresholds,thresholds_errors)
 #   print(threshold_dict)
    resonance_dict['Thresholds']=zip(thresholds,thresholds_errors)
    name="resonances_{}{}".format(curr_nucleus,name_ext)
    with open(name+'.pkl', 'wb') as file:
        pickle.dump(resonance_dict, file)
        print("file written successfully")

#    name="thresholds_{}{}".format(curr_nucleus,name_ext)
#    with open(name+'.pkl', 'wb') as file:
#        pickle.dump(threshold_dict, file)
#        print("file written successfully")
def plot_rotation(num_p,name=''):
    params=ask_data_params(['J','T','parity','hw','N_max','srg','lambda','NNN','3N_induced','Nf','srg_f'])
    subdata=data[check_params(params)]

    thresholds=set_threshold(params)
    n=subdata['theta'].nunique()
    thetas=subdata['theta'].unique()
    print("available Complex Scaling angles:",thetas)
    s=input("choose the desired values ( all in one line ), or -1 for all")
    thetas_in=[float(s1) for s1 in s.split()]
    if(thetas_in[0]==-1):
        thetas_in=thetas
    print(" working with: ",thetas_in)
    thetas_in=np.sort(thetas_in)
    n=len(thetas_in)
    ny=min(2,n)
    nx=ceil(n/2)
    figsize=fig_size(nx,ny,0.7 if ny==1 else 1)
    fig,axs=plt.subplots(nx,ny,figsize=figsize)
    if not isinstance(axs,np.ndarray):
        axs=np.array([[axs]])
    for i in range(n,nx*ny):
        fig.delaxes(axs[i//ny,i%ny])
    plt.rcParams.update({'font.size': 12})
    color=['r','g','b','y','c','tab:brown','m']
    ind=0
    print(axs.shape)
    for l, ax in enumerate(axs.flat):
    #for l in range(len(thetas_in)):
        if(l==n):
            break
        if(subdata[subdata['theta']==thetas_in[l]].empty):
            continue
        E_r=subdata[subdata['theta']==thetas_in[l]]['E_r']
        E_i=subdata[subdata['theta']==thetas_in[l]]['E_i']
        num_p=min(num_p,len(E_r))
        axj=ind%nx
        axi=int(ind/nx)
        ax.set_xlim(max(min(subdata['E_r'].min(),thresholds[curr_nucleus][0])-1,-40),12)
        ax.set_ylim(-9,2)
        ax.scatter(E_r[:num_p],E_i[:num_p],color=color[2],label=r"$2\theta={:0.2f}^\circ$".format(2*thetas_in[l]*180.0/math.pi),linewidths=1,clip_on=True,marker='o',s=50)
        ##########
        for j in thresholds[curr_nucleus]:
            xx=np.linspace(0,40,40)
            yy=-xx*math.tan(2*float(thetas_in[l]))
            xx[:]=xx[:]+j
            ax.plot(xx,yy,color='k',linestyle='--')
        ##########
        ax.axhline(y=0,color='k')
        ax.axvline(x=0,color='k')
        #ax[axi,axj].set_title("θ={:0.3f}".format(thetas[ind]),y=0,pad=-25,fontdict={'fontsize':14})
        ax.legend(loc='lower left', fontsize=9,markerscale=0.7)
        ax.locator_params(axis='both', nbins=5)
        ax.tick_params(axis='both', which='major', labelsize=11)
        ind+=1
    fig.subplots_adjust(left=0.13, bottom=0.15)
    E_name=r"\mathrm{{E^{{{0:}^{{{1:}}}{2:}}}}}".format(int(params['J'])
            if (A%2==0)else r"\frac{"+str(int(params['J']))+"}{2}",'+' if(params['parity']==1) else'-',int(params['T'])
            if (A%2==0) else r"\frac{"+str(int(params['T']))+"}{2}")
    #
    fig.text(0.5,0.03, r"$\mathcal{R}e\left("+E_name+r"\right)[\SI[]{}{\mega\eV}]$",ha='center',va='center')
    fig.text(0.05,0.5,r"$\mathcal{I}m\left("+E_name+r"\right)\ [\SI[]{}{\mega\eV}]$",ha='center',va='center',rotation='vertical')
    #rect = plt.Rectangle((-2.9, -1), 2.5, 2,
    #                     edgecolor='black',linestyle='dashed', linewidth=0.2, fill=False)
    #ax.add_patch(rect)
    #ax.set_xlabel("Real(E) [Mev]",fontsize=22)
    #ax.set_ylabel("Imaginary(E) [Mev]",fontsize=22)
    #ax.tick_params(axis='both', which='major', labelsize=20)
    #plt.tight_layout()
    plt.show()
    extension=create_name_extension(params)
    s=input("do you want to save the figure?y/n")
    if(s=='y'):
        save_figure(fig,'graphs/rotation/',name+curr_nucleus+extension)


#---------------------------------------------------------------------------
#plots the spectrum in the complex plan, each theta has its own sub-figure, with
#the possiblity of multiple points for an additional parameter (hw, lambda, N_max
#                                                              ..etc
def plot_rotation_v2(num_p,extra_param='',name=''):
    #params=ask_data_params(['J','T','parity','hw','N_max','srg','lambda','NNN','3N_induced','Nf','srg_f'])
    a=['J','T','parity','hw','N_max','srg','lambda','NNN','3N_induced','Nf','srg_f']
    if(extra_param!=''):
        a.remove(extra_param)
    params=ask_data_params(a)
    subdata=data[check_params(params)]

    ex_param_in,n_ex_param=ask_for_extra_param(extra_param,subdata)

    thresholds=set_threshold(params)

    n=subdata['theta'].nunique()
    thetas=subdata['theta'].unique()
    print("available Complex Scaling angles:",thetas)
    s=input("choose the desired values ( all in one line ), or -1 for all")
    thetas_in=[float(s1) for s1 in s.split()]
    if(thetas_in[0]==-1):
        thetas_in=thetas
    print(" working with: ",thetas_in)
    thetas_in=np.sort(thetas_in)
    n=len(thetas_in)
    ny=min(2,n)
    nx=ceil(n/2)
    figsize=fig_size(nx,ny)
    fig,axs=plt.subplots(nx,ny,figsize=figsize)
    if not isinstance(axs,np.ndarray):
        axs=np.array([[axs]])
    for i in range(n,nx*ny):
        fig.delaxes(axs[i//ny,i%ny])
    plt.rcParams.update({'font.size': 12})
    color=['r','g','b','y','c','tab:brown','m']
    ind=0
    print(axs.shape)
    for l, ax in enumerate(axs.flat):
    #for l in range(len(thetas_in)):
        if(l==n):
            break
        # extra params
        for i in range(n_ex_param):
            if(extra_param!=''):
                mask= subdata[extra_param]==ex_param_in[i]
                label=symbol[extra_param]+"$ ={}$".format(ex_param_in[i])
            else:
                mask= subdata['J']>-5
                label=""
        #

            if(subdata[(mask)& (subdata['theta']==thetas_in[l])].empty):
                continue
            E_r=subdata[(mask)&                        (subdata['theta']==thetas_in[l])]['E_r']
            E_i=subdata[(mask)& (subdata['theta']==thetas_in[l])]['E_i']
            num_p=min(num_p,len(E_r))
            ax.scatter(E_r[:num_p],E_i[:num_p],label=label,linewidths=1,clip_on=True,marker='+',s=30)
        axj=ind%nx
        axi=int(ind/nx)
        ax.set_xlim(max(min(subdata['E_r'].min(),thresholds[curr_nucleus][0])-1,-40),12)
        #ax.set_ylim(-10,2)
        ax.set_xlim(-10,40)
        ax.set_title(symbol['theta']+r"$ ={:.2f}\degree$".format(2*thetas_in[l]*180.0/math.pi))
        ##########
        for j in thresholds[curr_nucleus]:
            xx=np.linspace(0,40,40)
            yy=-xx*math.tan(2*float(thetas_in[l]))
            xx[:]=xx[:]+j
            ax.plot(xx,yy,color='k',linestyle='--',linewidth=0.5)
        ##########
        ax.axhline(y=0,color='k')
        ax.axvline(x=0,color='k')
        #ax[axi,axj].set_title("θ={:0.3f}".format(thetas[ind]),y=0,pad=-25,fontdict={'fontsize':14})
        ax.legend(loc='lower left', fontsize=9,markerscale=0.7)
        ax.locator_params(axis='both', nbins=5)
        ax.tick_params(axis='both', which='major', labelsize=11)
        ind+=1
    fig.subplots_adjust(left=0.13, bottom=0.15)
    E_name=r"\mathrm{{E^{{{0:}^{{{1:}}}{2:}}}}}".format(int(params['J'])
            if (A%2==0)else r"\frac{"+str(int(params['J']))+"}{2}",'+' if(params['parity']==1) else'-',int(params['T'])
            if (A%2==0) else r"\frac{"+str(int(params['T']))+"}{2}")
    #
    fig.text(0.5,0.05, r"$\mathrm{\Re}\left("+E_name+r"\right)[\SI[]{}{\mega\eV}]$",ha='center',va='center')
    fig.text(0.05,0.5,r"$\mathrm{\Im}\left("+E_name+r"\right)\ [\SI[]{}{\mega\eV}]$",ha='center',va='center',rotation='vertical')
    #rect = plt.Rectangle((-2.9, -1), 2.5, 2,
    #                     edgecolor='black',linestyle='dashed', linewidth=0.2, fill=False)
    #ax.add_patch(rect)
    #ax.set_xlabel("Real(E) [Mev]",fontsize=22)
    #ax.set_ylabel("Imaginary(E) [Mev]",fontsize=22)
    #ax.tick_params(axis='both', which='major', labelsize=20)
    #plt.tight_layout()
    plt.show()
    extension=create_name_extension(params)
    s=input("do you want to save the figure?y/n")
    if(s=='y'):
        save_figure(fig,'graphs/rotation/',name+curr_nucleus+extension)
#---------------------------------------------------------------------------

#---------------------------------------------------------------------------
#plots the spectrum in the complex plan, each theta has its own sub-figure.
#multiple hw values are shown on each sub-figure. the thresholds are estimated
#from the variation with hw
def plot_rotation_hw(num_p,name=''):
    """
    Plots the eigenvalues in the complex plane for various hw values, 
    and estimates thresholds for each scattering state.

    Parameters:
    - num_p (int): The number of states to plot.

    Returns:
    - None
    """
    #params=ask_data_params(['J','T','parity','hw','N_max','srg','lambda','NNN','3N_induced','Nf','srg_f'])
    extra_param='hw'
    a=['J','T','parity','hw','N_max','srg','lambda','NNN','3N_induced','Nf','srg_f']
    if(extra_param!=''):
        a.remove(extra_param)
    params=ask_data_params(a)
    subdata=data[check_params(params)]

    ex_param_in,n_ex_param=ask_for_extra_param(extra_param,subdata)

    thresholds=set_threshold(params)

    n=subdata['theta'].nunique()
    thetas=subdata['theta'].unique()
    print("available Complex Scaling angles:",thetas)
    s=input("choose the desired values ( all in one line ), or -1 for all")
    thetas_in=[float(s1) for s1 in s.split()]
    if(thetas_in[0]==-1):
        thetas_in=thetas
    print(" working with: ",thetas_in)
    thetas_in=np.sort(thetas_in)
    n=len(thetas_in)
    ny=min(2,n)
    nx=ceil(n/2)
    figsize=fig_size(nx,ny,0.7 if ny==1 else 1)
    fig,axs=plt.subplots(nx,ny,figsize=figsize)
    if not isinstance(axs,np.ndarray):
        axs=np.array([[axs]])
    for i in range(n,nx*ny):
        fig.delaxes(axs[i//ny,i%ny])
    plt.rcParams.update({'font.size': 12})
    color=['r','g','b','y','c','tab:brown','m']
    ind=0
    print(axs.shape)
    for l, ax in enumerate(axs.flat):
    #for l in range(len(thetas_in)):
        if(l==n):
            break
        # extra params
        for i in range(n_ex_param):
            if(extra_param!=''):
                mask= subdata[extra_param]==ex_param_in[i]
                label=symbol[extra_param]+"$ ={}$".format(ex_param_in[i])
            else:
                mask= subdata['J']>-5
                label=""
        #

            if(subdata[(mask)& (subdata['theta']==thetas_in[l])].empty):
                continue
            E_r=subdata[(mask)&                        (subdata['theta']==thetas_in[l])]['E_r']
            E_i=subdata[(mask)& (subdata['theta']==thetas_in[l])]['E_i']
            num_p=min(num_p,len(E_r))
            print(ex_param_in[i],num_p)
            ax.scatter(E_r[:num_p],E_i[:num_p],label=label,linewidths=1,clip_on=True,marker='+',s=30)
        axj=ind%nx
        axi=int(ind/nx)
        ax.set_xlim(max(min(subdata['E_r'].min(),thresholds[curr_nucleus][0])-1,-40),12)
        ax.set_ylim(-15,2)
        #ax.set_xlim(-2*7.1-1,-2*7.1+20)
        ax.set_xlim(-2,15)
        #ax.set_ylim(-15,2)
        ax.set_title(r'$2$'+symbol['theta']+r"$={:.2f}^\circ$".format(2*thetas_in[l]*180.0/math.pi),fontsize=10)
        ##########
        #thresholds={curr_nucleus:[-8,-7.3,-4.4,-2.2,0]}
        thresholds={curr_nucleus:[0]}
        #thresholds={curr_nucleus:[0,-7.1,2*-7.1]}
        #rotation_slope=[1.6,1.6,1,1,1]
        #rotation_slope=[1.65,1.6,1.35,1.6,1.0886]
        rotation_slope=[1,1,1,1,1]
        for i,j in enumerate(thresholds[curr_nucleus]):
            xx=np.linspace(0,40,40)
            yy=-xx*math.tan(2*float(thetas_in[l]))*rotation_slope[i]
            xx[:]=xx[:]+j
            ax.plot(xx,yy,color='k',linestyle='--',linewidth=0.5)
        ##########
        ax.axhline(y=0,color='k')
        ax.axvline(x=0,color='k')
        #ax[axi,axj].set_title("θ={:0.3f}".format(thetas[ind]),y=0,pad=-25,fontdict={'fontsize':14})
        ax.legend(loc='lower left', fontsize=9,markerscale=0.7)
        ax.locator_params(axis='both', nbins=5)
        ax.tick_params(axis='both', which='major', labelsize=11)
        ind+=1
    #fig.subplots_adjust(left=0.13, bottom=0.15)
    fig.subplots_adjust(left=0.16, bottom=0.18)
    E_name=r"\mathrm{{E^{{{0:}^{{{1:}}}{2:}}}}}".format(int(params['J'])
            if (A%2==0)else r"\frac{"+str(int(params['J']))+"}{2}",'+' if(params['parity']==1) else'-',int(params['T'])
            if (A%2==0) else r"\frac{"+str(int(params['T']))+"}{2}")
    #
    fig.text(0.5,0.03, r"$\mathcal{R}e\left("+E_name+r"\right)[\SI[]{}{\mega\eV}]$",ha='center',va='center')
    fig.text(0.03,0.5,r"$\mathcal{I}m\left("+E_name+r"\right)\ [\SI[]{}{\mega\eV}]$",ha='center',va='center',rotation='vertical')
    #rect = plt.Rectangle((-2.9, -1), 2.5, 2,
    #                     edgecolor='black',linestyle='dashed', linewidth=0.2, fill=False)
    #ax.add_patch(rect)
    #ax.set_xlabel("Real(E) [Mev]",fontsize=22)
    #ax.set_ylabel("Imaginary(E) [Mev]",fontsize=22)
    #ax.tick_params(axis='both', which='major', labelsize=20)
    #plt.tight_layout()
    plt.show()
    extension=create_name_extension(params)
    s=input("do you want to save the figure?y/n")
    if(s=='y'):
        save_figure(fig,'graphs/rotation/',name+curr_nucleus+extension)


    #threshold fitting
    def curv_resid(a,x,y):
        weights=[i for i in range(1,len(x)+1)]
        return (a[0]+a[1]*x-y)**2
    def linear_func(x,a,b):
        return a+b*x
    fig,axs=plt.subplots(nx,ny,figsize=figsize)
    if not isinstance(axs,np.ndarray):
        axs=np.array([[axs]])
    for i in range(n,nx*ny):
        fig.delaxes(axs[i//ny,i%ny])
    plt.rcParams.update({'font.size': 12})
    res_list=[]
    for l, ax in enumerate(axs.flat):
        intercepts=[]
        norm=0
        weights=[]
        for i in range(num_p):
            subdata2=subdata[subdata['theta']==thetas_in[l]]
            extra_param_values=subdata2[extra_param].unique()
            extra_param_values=ex_param_in
            extra_param_num=subdata2[extra_param].nunique()
            extra_param_num=n_ex_param

        #
            x=np.array([ subdata2.loc[ subdata2[extra_param]
                             ==extra_param_value,'E_r'].iloc[i] for
                        extra_param_value in extra_param_values])
            y=np.array([ subdata2.loc[ subdata2[extra_param]
                             ==extra_param_value,'E_i'].iloc[i] for
               extra_param_value in extra_param_values])
            print(i,extra_param_values)

            slope=-np.tan(2*thetas_in[l])
            intercept_guess= np.sum(y-slope*x)/extra_param_num
            intercept_guess=max(intercept_guess,-50)
            intercept_guess=min(intercept_guess,0)
            x0=[intercept_guess,slope]
            bounds=([-50,slope*1.8],[0,slope*0.6])
            least1=least_squares(curv_resid,x0,loss='soft_l1',f_scale=0.1,bounds=bounds,args=(x[:],y[:]))
            y_pred=linear_func(x,*least1.x)
            #slop=-np.tan(2*thetas_in[l])
            #intercept= np.sum(y-slop*x)/extra_param_num
            interceptx=-least1.x[0]/least1.x[1]
            #y_pred=slop*x+intercept
            y_mean=np.sum(y)/extra_param_num
            ss_res=np.sum((y-y_pred)**2)
            ss_tot=np.sum((y-y_mean)**2)
            r_squared=1-ss_res/ss_tot
            if(ss_tot==0):
                r_squared=1
            print(interceptx,least1.x[1]/slope,r_squared)
            if(r_squared>0.8 and x[0]>-15):
                intercepts.append(interceptx)
                weights.append(np.exp((r_squared-0.97)))
                norm+=np.exp((r_squared-0.97))
            ax.scatter(x,y)
            if(r_squared>0.6):
                xxx=np.linspace(0,40,100)
                yyy=xxx*least1.x[1]
                ax.plot(xxx+interceptx,yyy,color='k')
            plt.draw()
            if(r_squared<0):
                print("possible resonance in state {}".format(i), x[0],y_mean)
                plt.pause(2)
                res_list.append(i)

            plt.pause(0.2)
        plt.pause(10)
        #test
        #intercepts.append(1)
        #weights.append(np.exp((1-0.97)))
        #norm+=np.exp((1-0.97))
        #
        weights=weights/norm
        plt.clf()
        figsize=fig_size(1,1,0.8)
        fig,axs=plt.subplots(1,1,figsize=figsize)
        kde=sns.kdeplot(x=intercepts,weights=weights, bw_adjust=0.14,ax=axs)
        from scipy.signal import find_peaks, peak_widths
        x, y = kde.get_lines()[0].get_data()
        peaks, _ = find_peaks(y)
        peak_centers = x[peaks]
        results_half = peak_widths(y, peaks, rel_height=0.5)
        widths = results_half[0] / 2.355
        left_ips_x = x[np.round(results_half[2]).astype(int)]
        right_ips_x = x[np.round(results_half[3]).astype(int)]
        widths = (right_ips_x-left_ips_x) / 2.355
        axs.plot(x, y)
        axs.plot(x[peaks], y[peaks], "x")
        threshold_labels=[r'$\mathrm{^3H\text{-}p}$',r'$\mathrm{^3He\text{-}n}$',r'$\mathrm{D\text{-}D}$',r'$\mathrm{D\text{-}np}$',r'$\mathrm{2n2p}$']
        for i in range(len(threshold_labels)):
            axs.annotate(threshold_labels[i],(x[peaks][i],y[peaks][i]),textcoords='offset points',xytext=(5,4),fontsize=12)
        axs.hlines(results_half[1],left_ips_x,right_ips_x, color="C2")
        axs.set_xlabel(r'$\mathrm{E\ [MeV]}$',fontsize=12)
        axs.set_ylabel(r'$\mathrm{Probability\ density}$',fontsize=12)
        plt.tight_layout()
        plt.show()
        s=input('do you want to save the graph?(y/n)')
        if(s=='y'):
            save_figure(fig,'graphs/rotation/',name+curr_nucleus+'thresholds')
        for i, center in enumerate(peak_centers):
            print(f"Peak {i+1}: Center = {center}, Approx Width (std dev) = {widths[i]}, height={y[peaks[i]]}")

        figsize=fig_size(1,1)
        fig2,axs2=plt.subplots(1,1,figsize=figsize)
        sns.kdeplot(x=intercepts, weights=weights, bw_adjust=0.18,fill=True,color='blue', alpha=0.6, ax=axs2)
        axs2.set_xlim(-12,2)
        plt.show()
    print(res_list)
            #slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
#---------------------------------------------------------------------------


def plot_convergence(st_num,extra_param='',name='',extrapolate=False):
    #begining of user defined params
    plot_params={'marker':'s','markersize':3,'linewidth':0.7,'linestyle':'-','layout':'h','ylims':'n','ylim1':[-10,-2],'ylim2':[-5,1],'logscale':True}
    #plot_params={'marker':'s','markersize':3,'linewidth':0.7,'linestyle':'-','layout':'h','ylims':'n','ylim1':[-33,-25],'ylim2':[-2,2],'logscale':False}
    #end of user defined params
    a=['J','T','parity','hw','srg','lambda','theta','NNN','3N_induced','Nf','srg_f']
    if(extra_param!=''):
        a.remove(extra_param)
    params=ask_data_params(a)
    subdata=data[check_params(params)]
    if(st_num>0):
        subdata=subdata[subdata['N_max']>max(2,st_num)]
    n_ex_param=1
    if(extra_param!=''):
        ex_param=subdata[extra_param].unique()
        print("available "+extra_param+"  values:",ex_param)
        s=input("choose the desired values ( all in one line ), or -1 for all")
        ex_param_in=[float(s1) for s1 in s.split()]
        if(ex_param_in[0]==-1):
            ex_param_in=ex_param
        ex_param_in=np.sort(ex_param_in)
        n_ex_param=len(ex_param_in)
        print(" working with: ",n_ex_param,ex_param_in)
    n=subdata['N_max'].nunique()
    thresholds=set_threshold(params)
    if(plot_params['layout']=='h'):
        nx=1
        ny=2
    else:
        nx=2
        ny=1
    figsize=fig_size(nx,ny)
    fig,ax=plt.subplots(nx,ny,figsize=figsize,constrained_layout=False)
    #plt.subplots_adjust(wspace=0.2, hspace=0.2)
    plt.rcParams.update({'font.size': 12})
    N_values=subdata['N_max'].unique()
    N_values=np.sort(N_values)
    print(N_values)
    for i in range(n_ex_param):
        if(extra_param!=''):
            mask= subdata[extra_param]==ex_param_in[i]
            label=symbol[extra_param]+"$ ={}$".format(ex_param_in[i])
        else:
            mask= subdata['J']>-5
            label=""
        N_values2=np.sort(subdata[mask]['N_max'].unique())
        x_values=[subdata[(subdata['N_max']==N_value) &
                          (mask)]['E_r'].iloc[st_num] for N_value in N_values2]
        ax[0].plot(N_values2,x_values,marker=plot_params['marker'],markersize=plot_params['markersize'],linestyle=plot_params['linestyle'],linewidth=plot_params['linewidth'],label=label)
        if(i==0 or n_ex_param<3):
            ax[0].axhline(y=x_values[-1],color='k',linewidth=1,linestyle=':')
        #ax[0].tick_params(axis='both', which='major', labelsize=22)
        y_values=[subdata[(subdata['N_max']==N_value)&(mask)]['E_i'].iloc[st_num]
                  for N_value in N_values2]
        ax[1].plot(N_values2,y_values,marker=plot_params['marker'],markersize=plot_params['markersize'],linestyle=plot_params['linestyle'],linewidth=plot_params['linewidth'],label=label)

        if(i==0 or n_ex_param<3):
            ax[1].axhline(y=y_values[-1],color='k',linewidth=1,linestyle=':')
    # add grey box for converged value
    y_min, y_max = ax[0].get_ylim()
    new_label=format_tick_labels([x_values[-1]])[0]
    #ax[0].text(N_values[0]+0.5,x_values[-1]-0.05*(y_max-y_min),new_label,fontsize=12,color='grey', ha='center',va='center',bbox=dict(facecolor='white', edgecolor='grey',alpha=0.5, boxstyle='round,pad=0.5'))
    y_min, y_max = ax[1].get_ylim()
    new_label=format_tick_labels([y_values[-1]])[0]
    #ax[1].text(N_values[0]+0.5,y_values[-1]-0.05*(y_max-y_min),new_label,fontsize=12,color='grey', ha='center',va='center',bbox=dict(facecolor='white', edgecolor='grey',alpha=0.5, boxstyle='round,pad=0.5'))
    ##
    ##
    if(extra_param!=''):
        ax[0].legend(fontsize=9,labelspacing=0.2)#,loc='upper right')
        #ax[1].legend(fontsize=9,labelspacing=0.2)#,loc='upper right')
    E_name=r"\mathrm{{E^{{{0:}^{{{1:}}}{2:}}}_{{{3:d}}}}}".format(int(params['J'])
            if (A%2==0)else r"\frac{"+str(int(params['J']))+"}{2}",'+' if(params['parity']==1) else'-',int(params['T'])
            if (A%2==0) else r"\frac{"+str(int(params['T']))+"}{2}",int(st_num))
    ax[0].set_xlabel(r"$\mathrm{ N_{max}}$",fontsize=12)
    ax[0].set_ylabel(r"$\mathcal{R}e\left("+E_name+r"\right)[\SI[]{}{\mega\eV}]$",fontsize=12)
    ax[1].set_xlabel(r"$\mathrm{N_{max}}$",fontsize=12)
    ax[1].set_ylabel(r"$\mathcal{I}m\left("+E_name+r"\right)\ [\SI[]{}{\mega\eV}]$",fontsize=12)
    ax[1].tick_params(axis='both', which='major', labelsize=11)
    ax[0].tick_params(axis='both', which='major', labelsize=11)
    ax[0].locator_params(axis='x', nbins=5)
    ax[1].locator_params(axis='x', nbins=5)
    ax[0].minorticks_on()
    ax[1].minorticks_on()
    if( plot_params['logscale']):
        ax[0].set_yscale('symlog')
        ax[1].set_yscale('symlog')
    #
    if(plot_params['ylims']=='y'):
        y1=plot_params['ylim1']
        y2=plot_params['ylim2']
        ax[0].set_ylim(y1[0],y1[1])
        ax[1].set_ylim(y2[0],y2[1])
    #
    ax[0].xaxis.set_minor_locator(AutoMinorLocator(2))
    ax[1].xaxis.set_minor_locator(AutoMinorLocator(2))
    plt.tight_layout()
    extension=create_name_extension(params)
    filename='convergence_'+name+extra_param+curr_nucleus+extension+"_{:d}".format(st_num)
    fig.subplots_adjust(left=0.13, bottom=0.23)
    plt.show()
    s=input("do you want to save the figure?y/n")
    if(s=='y'):
        if(extra_param==''):
            save_figure(fig,'graphs/convergence/',filename)
        else:
            save_figure(fig,'graphs/convergence_'+extra_param+'/',filename)

### plot convergence in the complex plane
    fig,ax=plt.subplots(1,1,figsize=(16*0.7,9*0.7))
    plt.rcParams.update({'font.size': 22})
    color_spect = np.linspace(0, 1, len(N_values))
    x_values=[subdata[(subdata['N_max']==N_value)]['E_r'].iloc[st_num] for N_value in N_values]
    y_values=[subdata[(subdata['N_max']==N_value)]['E_i'].iloc[st_num] for N_value in N_values]
    for i in range(n_ex_param):
        if(extra_param!=''):
            mask= subdata[extra_param]==ex_param_in[i]
            label=symbol[extra_param]+"$ ={}$".format(ex_param_in[i])
        else:
            mask= subdata['J']>-5
            label=""
        N_values2=np.sort(subdata[mask]['N_max'].unique())
        x_values=[subdata[(subdata['N_max']==N_value) &
                          (mask)]['E_r'].iloc[st_num] for N_value in N_values2]
        y_values=[subdata[(subdata['N_max']==N_value) &
                          (mask)]['E_i'].iloc[st_num] for N_value in N_values2]
        if(extra_param==''):
            scatter=ax.scatter(x_values,y_values,marker='+',c=N_values,cmap=plt.cm.viridis)
        else:
            scatter=ax.scatter(x_values,y_values,marker='+',label=label)
    for j in thresholds[curr_nucleus]:
        xx=np.linspace(0,40,40)
        yy=-xx*math.tan(2*float(params['theta']))
        xx[:]=xx[:]+j
        ax.plot(xx,yy,color='k',linestyle='--')
        ##########
    ax.axhline(y=0,color='k')
    ax.axvline(x=0,color='k')
    ax.set_xlim(min(min(x_values),thresholds[curr_nucleus][0])-1,12)
    ax.set_ylim(-15,2)
    if(extra_param==''):
        cbar = fig.colorbar(scatter,ax=ax)
        cbar.set_label(r"$\mathrm{ N_{max}}$")
    else:
        ax.legend(fontsize=18)

    ax.set_xlabel("$\mathrm{\Re}"+E_name+r" [\SI[]{}{\mega\eV}]$",fontsize=22)
    ax.set_ylabel(r"$\mathrm{\Im}"+E_name+"\ [\SI[]{}{\mega\eV}]$",fontsize=22)
    plt.show()
### END of plot convergence in the complex plane

### extrapolation section
    if( extrapolate):
        figsize=fig_size(nx,ny)
        fig,ax=plt.subplots(nx,ny,figsize=figsize)
        def curv_resid(a,x,y):
            weights=[i for i in range(1,len(x)+1)]
            return (a[0]+a[1]*np.exp(-a[2]*x)-y)**2
        def curv_func(x,a,b,c):
            return a+b*np.exp(-c*x)
        for i in range(n_ex_param):
            if(extra_param!=''):
                mask= subdata[extra_param]==ex_param_in[i]
                label=symbol[extra_param]+"$ ={}$".format(ex_param_in[i])
            else:
                mask= subdata['J']>-5
                label=""
            x0=np.array([1,0,0.2])
            num_points=4
            bounds=([-40,-150,0.1],[30,150,1])
            N_full=np.linspace(4,30,15)
            x_values=[subdata[(subdata['N_max']==N_value) & (mask)]['E_r'].iloc[st_num] for N_value in N_values]
            y_values=[subdata[(subdata['N_max']==N_value)&(mask)]['E_i'].iloc[st_num] for N_value in N_values]
            #del(x_values[-2])
            #del(y_values[-2])
            #N_values=np.delete(N_values,n-2)
            x0[0]=x_values[-1]
            least1=least_squares(curv_resid,x0,loss='soft_l1',f_scale=0.1,bounds=bounds,args=(N_values[n-num_points:],x_values[n-num_points:]))
            x_pred=curv_func(N_full[:],*least1.x)
            label1=label#+"\n $E_{{\inf}}={:.2f}$".format(least1.x[0])
            ax[0].plot(N_full[4:],x_pred[4:],markersize=10,linestyle=':',label=label1)
            ax[0].plot(N_values,x_values,marker=plot_params['marker'],markersize=plot_params['markersize'],linestyle=plot_params['linestyle'],linewidth=plot_params['linewidth'])


            x0[0]=y_values[-1]
            least2=least_squares(curv_resid,x0,loss='soft_l1',f_scale=0.1,bounds=bounds,args=(N_values[n-num_points:],y_values[n-num_points:]))
            y_pred=curv_func(N_full[:],*least2.x)
            label2=label#+"\n $E_{{\inf}}={:.2f}$".format(least2.x[0])
            ax[1].plot(N_full[4:],y_pred[4:],markersize=10,linestyle=':',label=label2)
            ax[1].plot(N_values,y_values,marker='x',markersize=10)
            residual1=x_pred[n-3:n]-x_values[n-3:]
            residual2=y_pred[n-3:n]-y_values[n-3:]
            d1_error=np.abs(residual1).sum()/3
            d2_error=np.abs(residual2).sum()/3
            E_r_conv=x_pred[n-1]-least1.x[0]
            E_i_conv=y_pred[n-1]-least2.x[0]
            E_r_err=np.sqrt(d1_error**2+E_r_conv**2)
            E_i_err=np.sqrt(d2_error**2+E_i_conv**2)
            #printing results:
            print("******\nfit results for {}:".format(label))
            print(residual1,residual2)
            print("fit params:",least1.x[:],least2.x[:])
            print("E_inf:", least1.x[0],least2.x[0])
            print("final value:",x_values[-1],y_values[-1])
            print("convergence in E_r:",d1_error,E_r_conv,E_r_err)
            print(x_values[-1]-x_values[-2],x_values[-2]-x_values[-3])
            print("convergence in E_i:",d2_error,E_i_conv,E_i_err)
            print(y_values[-1]-y_values[-2],y_values[-2]-y_values[-3])
            if(extra_param!=''):
                ext_params={**params,extra_param:ex_param_in[i]}
            else:
                ext_params=params
            # #
        #ax[0].legend(fontsize=9)
        #ax[1].legend(fontsize=9)
        ax[0].set_xlabel(r"$\mathrm{ N_{max}}$",fontsize=11)
        ax[0].set_ylabel("$\mathrm{\Re}\left("+E_name+r"\right)[\SI[]{}{\mega\eV}]$",fontsize=11)
        ax[1].set_xlabel(r"$\mathrm{N_{max}}$",fontsize=11)
        ax[1].set_ylabel(r"$\mathrm{\Im}\left("+E_name+r"\right)\ [\SI[]{}{\mega\eV}]$",fontsize=11)
        ax[1].tick_params(axis='both', which='major', labelsize=12)
        ax[0].tick_params(axis='both', which='major', labelsize=12)
        ax[0].locator_params(axis='x', nbins=10)
        ax[1].locator_params(axis='x', nbins=10)
        ax[0].minorticks_on()
        ax[1].minorticks_on()
        ax[0].xaxis.set_minor_locator(AutoMinorLocator(2))
        ax[1].xaxis.set_minor_locator(AutoMinorLocator(2))
        plt.tight_layout()
        plt.show()
        s=input("do you want to save the figures?y/n")
        if(s=='y'):
            extension=create_name_extension(params)
            filename='convergence_extrapolation_'+name+curr_nucleus+extension+"_{:d}".format(st_num)
            if(extra_param!=''):
                save_figure(fig,'graphs/convergence_'+extra_param+'/',filename)
            save_figure(fig,'graphs/convergence/',filename)


### resonance tracking section
def ang_v(o_x,o_y,row1,row2):#x1,y1,x2,y2,dt):
    x1=row1['E_r']
    x2=row2['E_r']
    y1=row1['E_i']
    y2=row2['E_i']
    dt=row2['theta']-row1['theta']
    x_avg=(x1+x2)/2.0
    y_avg=(y1+y2)/2.0
    Lx=x_avg-o_x
    Ly=y_avg-o_y
    L=np.sqrt(Ly**2+Lx**2)
    dx=x2-x1
    dy=y2-y1
    d=np.sqrt(dx**2+dy**2)
    d_L=(dx*Lx+dy*Ly)/L
    sign=np.sign((Lx*dy-Ly*dx)*-1)
    dw=np.sqrt(d**2-d_L**2)/L
    return dw/dt*sign
def dir_v(row1,row2):#,x1,y1,x2,y2,dt):
    x1=row1['E_r']
    x2=row2['E_r']
    y1=row1['E_i']
    y2=row2['E_i']
    dt=row2['theta']-row1['theta']
    dx=x2-x1
    dy=y2-y1
    d=np.sqrt(dx**2+dy**2)/dt
    return d

def track_state_rotation(states=[1],name=''):
    """
    Tracks the rotation of a set of states to determine resonances 
    based on predictions from the ABC theorem.

    Parameters:
    - states (list, optional): A list of states to track. Defaults to an empty list.

    Returns:
    - None
    """
    n_res=states[-1]+1
    n_plots=len(states)
    plt.rcParams.update({'font.size': 12})
    params=ask_data_params(['J','T','parity','hw','N_max','srg','lambda','NNN','3N_induced','Nf','srg_f'])
    subdata=data[check_params(params)].copy()
    thresholds=set_threshold(params)[curr_nucleus]

    #
    thetas=subdata['theta'].unique()[:]
    print("available theta values:",thetas)
    s=input("choose the desired values ( all in one line ), or -1 for all")
    theta_in=[float(s1) for s1 in s.split()]
    if(theta_in[0]==-1):
        theta_in=thetas
    theta_in=np.sort(theta_in)
    n_thetas=len(theta_in)
    print(" working with: ",n_thetas,theta_in)
    thetas=theta_in
    #

    ny=min(2,n_plots)
    nx=ceil(n_plots/2)
    figsize=fig_size(nx,ny,0.7 if ny==1 else 1)
    fig,axs=plt.subplots(nx,ny,figsize=figsize)
    if not isinstance(axs,np.ndarray):
        axs=np.array([[axs]])
    for i in range(n_plots,nx*ny):
        fig.delaxes(axs[i//ny,i%ny])
    color=['r','g','b','y','c','tab:brown','m']
    ind=0
    subdata['v']=None
    subdata['E_th']=None
    print(axs.shape)
    n_p=16
    n_res=np.min(np.array([subdata[subdata['theta']==theta]['E_r'].nunique() for theta in thetas]))
    color=['r','g','b','y','c','tab:brown','m']
    subdata.loc[:,'E_th']=subdata['E_r']*np.cos(2*subdata['theta'])-subdata['E_i']*np.sin(2*subdata['theta'])
    #w=np.zeros((n_res,n_thetas,len(thresholds)))
    #threshold=np.zeros(100)
    print(len(axs.flat))
    print(thetas,n_thetas)
    for i in range(n_thetas):
        subdata.loc[subdata['theta']==thetas[i],:]=subdata[subdata['theta']==thetas[i]].sort_values(by='E_th')
    for j in range(n_res):
        for i in range(n_thetas-1):
            row1=subdata[ subdata['theta']==thetas[i]].iloc[j]
            row2=subdata[subdata['theta']==thetas[i+1]].iloc[j]
            index=subdata.index[subdata['theta']==thetas[i]][j]
            subdata.loc[index,'v']=dir_v(row1,row2)
#            for th in range(len(thresholds)):
#                w[:n_res,i,th]=ang_v(thresholds[th],0,E_r_sorted[i][:n_res],E_i_sorted[i][:n_res],E_r_sorted[i+1][:n_res],E_i_sorted[i+1][:n_res],thetas[i+1]-thetas[i])

        #for i in range(n_thetas-1):
            #for th in range(len(thresholds)):
                #w[:n_res,i,th]=v[:n_res,i]/np.sqrt((0.5*E_r_sorted[i+1][:n_res]+0.5*E_r_sorted[i][:n_res]-thresholds[th])**2+(0.5*E_i_sorted[i+1][:n_res]+0.5*E_i_sorted[i][:n_res])**2)


    for j,ax in enumerate(axs.flat):
        if(j>=n_plots):
            break
        st_num=states[j]
        v=[subdata[subdata['theta']==theta]['v'].iloc[st_num] for theta in
           thetas[:n_thetas-1]]
        ax.plot(thetas[:n_thetas-1]*180/math.pi,v[:],color=color[2],label="state_num={0:d}".format(st_num),marker='o')
        ax.axhline(y=0,color='k')
        ax.axvline(x=0,color='k')
        #ax.legend( fontsize=9,labelspacing=0.2)
        ax.tick_params(axis='both', which='major', labelsize=11)
    fig.text(0.5,0.06, r"$\mathrm{\theta\ [deg]}$",ha='center',va='center',fontsize=12)
    fig.text(0.05,0.5, r"$\mathrm{\frac{d}{d\theta}E}\ [\SI[]{}{\mega\eV}]$",ha='center',va='center',rotation='vertical')
    fig.subplots_adjust(left=0.13, bottom=0.23)
    plt.show()
    save_figure(fig,'graphs/resonances/','speed_test')

#    fig.savefig(name+".pdf")
#    fig.savefig(name+".svg")
#    fig,ax=plt.subplots(3,3)
    if(thetas[1]<0.06):
        extracting_thresholds=True
        print("extracting thresholds from {} points".format(n_res))
        n_thresholds=1
        extracted_thresholds=np.ones(n_res)
        for i in range(n_res):
            v=subdata[subdata['theta']==thetas[0]]['v'].iloc[i]
            r=v/2
            extracted_thresholds[i]=subdata[subdata['theta']==thetas[0]]['E_r'].iloc[i]-r
            print("threshold :",extracted_thresholds[i])
    else:
        n_thresholds=len(thresholds)
        extracting_thresholds=False
    #fig,axs=plt.subplots(nx,ny,figsize=(16,int(9*nx/ny)))
    #fig,axs=plt.subplots(nx,ny,figsize=(int(16*ny/nx),5.1)
    figsize=fig_size(nx,ny,0.7 if ny==1 else 1)
    fig,axs=plt.subplots(nx,ny,figsize=figsize)
    #subdata = subdata.assign(vw=pd.Series(dtype='object'))
    subdata['wv'] = None
    subdata['wv'] = subdata['wv'].astype(object)
    if not isinstance(axs,np.ndarray):
        axs=np.array([[axs]])
    for i in range(n_plots,nx*ny):
        fig.delaxes(axs[i//ny,i%ny])
    for j in range(n_res):
        for i in range(n_thetas-1):
            row1=subdata[ subdata['theta']==thetas[i]].iloc[j]
            row2=subdata[subdata['theta']==thetas[i+1]].iloc[j]
            index=subdata.index[subdata['theta']==thetas[i]][j]
            subdata.at[index,'wv']=[]
            if(extracting_thresholds):
                thresholds=[extracted_thresholds[j]]
            for threshold in (thresholds):
                subdata.at[index,'wv'].append(ang_v(threshold,0,row1,row2))

    vw=[subdata[subdata['theta']==theta]['wv'].iloc[0][0] for theta in thetas[:n_thetas-1]]
    E_name=r"$\mathrm{{E^{{{0:}^{{{1:}}}{2:}}}}}$".format(int(params['J'])
            if (A%2==0)else r"\frac{"+str(int(params['J']))+"}{2}",'+' if(params['parity']==1) else'-',int(params['T'])
            if (A%2==0) else r"\frac{"+str(int(params['T']))+"}{2}")
    avg_thetas=(thetas[:n_thetas-1]+thetas[1:])/2.0
    for j,ax in enumerate(axs.flat):
        if(j>=n_plots):
            break
        st_num=states[j]
        if(extracting_thresholds):
            thresholds=[extracted_thresholds[st_num]]
        for i in range(n_thresholds):
            vw=[subdata[subdata['theta']==theta]['wv'].iloc[st_num][i] for theta in thetas[:n_thetas-1]]
            E0=subdata[subdata['theta']==thetas[0]]['E_r'].iloc[st_num]
            E_f=subdata[subdata['theta']==thetas[-1]]['E_r'].iloc[st_num]
            if(abs(vw[0]-2)>100 or (E0+100< thresholds[i] and E0>-10)):
                continue
            if(E0> thresholds[i]):
                linestyle='-'
            else:
                linestyle=':'
            ax.plot(avg_thetas[:]*180/math.pi,vw[:],label="$\mathrm{{th={0:.2f}}}$".format(thresholds[i]),marker='+',linestyle=linestyle)
        ax.axhline(y=0,color='k')
        ax.axvline(x=0,color='k')
        ax.axhline(y=2,color='r')
        ax.set_ylim(-0.5,10)
        #if(j==0):
            #ax.legend( fontsize=9,loc='upper right',labelspacing=0.2)#, bbox_to_anchor=(1.05, 0.5))

       # ax.set_xlabel("$\mathrm{Complex-Scaling\ angle}\
       #               [\SI[]{}{\degree}]$",fontsize=12,labelpad=-5)
       # ax.set_ylabel("$\mathrm{Rotation\ rate\ in}\ $" +
       #      E_name,fontsize=12,labelpad=10)
       # ax.set_ylabel("$\mathrm{Rotation\ rate\ in}\ $" +
       #      E_name,fontsize=18)
        ax.tick_params(axis='both', which='major', labelsize=11)
        ax.locator_params(axis='x', nbins=4)
        #ax[0].minorticks_on()
        #ax[0].xaxis.set_minor_locator(AutoMinorLocator(2))
    fig.text(0.5,0.05, r"$\mathrm{\theta}\ [^\circ]$",ha='center',va='center')
    fig.text(0.01,0.5, r"$\mathrm{\omega_{th}}\ $",ha='center',va='center',rotation='vertical')
    fig.subplots_adjust(left=0.13, bottom=0.23)
    plt.tight_layout()
    plt.show()
    extension=create_name_extension(params)
    filename="resonance_track_"+str(states)+"_"+name+curr_nucleus+extension
    s=input("do you want to save the figure?y/n")
    if(s=='y'):
        save_figure(fig,'graphs/resonances/',filename)
    figsize=fig_size(1,1)
    kde=sns.kdeplot(x=extracted_thresholds, bw_adjust=0.13)
    from scipy.signal import find_peaks, peak_widths
    x, y = kde.get_lines()[0].get_data()
    peaks, _ = find_peaks(y)
    peak_centers = x[peaks]
    results_half = peak_widths(y, peaks, rel_height=0.5)
    widths = results_half[0] / 2.355
    left_ips_x = x[np.round(results_half[2]).astype(int)]
    right_ips_x = x[np.round(results_half[3]).astype(int)]
    widths = (right_ips_x-left_ips_x) / 2.355
    plt.plot(x, y)
    plt.plot(x[peaks], y[peaks], "x")
    plt.hlines(results_half[1],left_ips_x,right_ips_x, color="C2")
    plt.show()
    for i, center in enumerate(peak_centers):
        print(f"Peak {i+1}: Center = {center}, Approx Width (std dev) = {widths[i]}, height={y[peaks[i]]}")
    fig2,axs2=plt.subplots(1,1,figsize=figsize)
    sns.kdeplot(x=extracted_thresholds, bw_adjust=0.13,fill=True,color='blue', alpha=0.6, ax=axs2)
    axs2.set_xlim(-12,2)
    plt.show()


def plot_convergence_inter(st_num,extra_param='',name='',extrapolate=False):
    #begining of user defined params
    name_ext=''
    plot_params={'marker':'s','markersize':3,'linewidth':0.7,'linestyle':'-','layout':'h','ylims':'y','ylim1':[-6,0],'ylim2':[-2,2]}
    #end of user defined params
    a=['J','T','parity','hw','lambda','theta','Nf','srg_f']
    params=ask_data_params(a)
    subdata=data[check_params(params)]
    inter={r"NN\textnormal{-}bare":{"NNN":0,"srg":0,"3N_induced":0},"NN+SRG":{"NNN":0,"srg":1,"3N_induced":0},r"NN+3N\textnormal{-}ind":{"NNN":0,"srg":1,"3N_induced":1},"NN+3N":{"NNN":1,"srg":1,"3N_induced":1}}
    if(st_num>0):
        subdata=subdata[subdata['N_max']>max(2,st_num)]
    n=subdata['N_max'].nunique()
    if(plot_params['layout']=='h'):
        nx=1
        ny=2
    else:
        nx=2
        ny=1
    figsize=fig_size(nx,ny)
    fig,ax=plt.subplots(nx,ny,figsize=figsize,constrained_layout =True)
    plt.rcParams.update({'font.size': 12})
    N_values=subdata['N_max'].unique()
    N_values=np.sort(N_values)
    print(N_values)
    for inter_i,inter_params in inter.items():
        data_i=subdata[check_params(inter_params,subdata)]
        if(data_i.empty):
            print("No ",inter_i," data")
            continue
        else:
            include=input("do you want to add {}?y/n".format(inter_i))
            if(include=='n'):
                continue
        name_ext=name_ext+"I_"
        label=r"$\mathrm{"+inter_i+"}$"
        N_values2=np.sort(data_i['N_max'].unique())
        x_values=[data_i[(data_i['N_max']==N_value)]['E_r'].iloc[st_num] for N_value in N_values2]
        ax[0].plot(N_values2,x_values,marker=plot_params['marker'],markersize=plot_params['markersize'],linestyle=plot_params['linestyle'],linewidth=plot_params['linewidth'],label=label)
        if(True):
            ax[0].axhline(y=x_values[-1],color='k',linestyle=":")
        #ax[0].tick_params(axis='both', which='major', labelsize=22)
        y_values=[data_i[(data_i['N_max']==N_value)]['E_i'].iloc[st_num]
                  for N_value in N_values2]
        ax[1].plot(N_values2,y_values,marker=plot_params['marker'],markersize=plot_params['markersize'],linestyle=plot_params['linestyle'],linewidth=plot_params['linewidth'],label=label)
        if(True):
            ax[1].axhline(y=y_values[-1],color='k',linestyle=":")
    y_min, y_max = ax[0].get_ylim()
    new_label=format_tick_labels([x_values[-1]])[0]
   # ax[0].text(N_values[0]+0.5,x_values[-1]-0.05*(y_max-y_min),new_label,fontsize=18,color='grey', ha='center',va='center',bbox=dict(facecolor='white', edgecolor='grey',alpha=0.5, boxstyle='round,pad=0.5'))
    y_min, y_max = ax[1].get_ylim()
    new_label=format_tick_labels([y_values[-1]])[0]
    #ax[1].text(N_values[0]+0.5,y_values[-1]-0.05*(y_max-y_min),new_label,fontsize=18,color='grey', ha='center',va='center',bbox=dict(facecolor='white', edgecolor='grey',alpha=0.5, boxstyle='round,pad=0.5'))
    ax[0].legend(fontsize=8,loc='upper right',markerscale=0.5)
    ax[1].legend(fontsize=8,loc='lower right',markerscale=0.5)
    E_name=r"\mathrm{{E^{{{0:}^{{{1:}}}{2:}}}_{{{3:d}}}}}".format(int(params['J'])
            if (A%2==0)else r"\frac{"+str(int(params['J']))+"}{2}",'+' if(params['parity']==1) else'-',int(params['T'])
            if (A%2==0) else r"\frac{"+str(int(params['T']))+"}{2}",int(st_num))
    #
    if(plot_params['ylims']=='y'):
        y1=plot_params['ylim1']
        y2=plot_params['ylim2']
        ax[0].set_ylim(y1[0],y1[1])
        ax[1].set_ylim(y2[0],y2[1])
    #
    ax[0].set_xlabel(r"$\mathrm{ N_{max}}$",fontsize=12)
    ax[0].set_ylabel(r"$\mathcal{R}e\left("+E_name+r"\right)[\SI[]{}{\mega\eV}]$",fontsize=12)
    ax[1].set_xlabel(r"$\mathrm{N_{max}}$",fontsize=12)
    ax[1].set_ylabel(r"$\mathcal{I}m\left("+E_name+r"\right)\ [\SI[]{}{\mega\eV}]$",fontsize=12)
    ax[1].tick_params(axis='both', which='major', labelsize=11)
    ax[0].tick_params(axis='both', which='major', labelsize=11)
    ax[0].locator_params(axis='x', nbins=10)
    ax[1].locator_params(axis='x', nbins=10)
    ax[0].minorticks_on()
    ax[1].minorticks_on()
    ax[0].xaxis.set_minor_locator(AutoMinorLocator(2))
    ax[1].xaxis.set_minor_locator(AutoMinorLocator(2))
    #plt.tight_layout()
    extension=create_name_extension(params)
    filename="inter_"+name_ext+name+curr_nucleus+extension+"_{:d}".format(st_num)
    plt.show()
    s=input("do you want to save the figure?y/n")
    if(s=='y'):
        save_figure(fig,'graphs/convergence/',filename)
def general_plot(par_x,st_num,extra_param='',y_shift=[0,0],name=''):
    """
    Plots the energy of a state as a function of a general parameter <par_x>. 
    Optionally, it can plot multiple curves for different values of the <extra_param>.
    Extrapolation is performed depending on the <par_x> value.

    Parameters:
    - par_x (str): The parameter to be plotted on the x-axis.
    - st_num (int): The index of the state to be plotted (in order of increasing real energy).
    - extra_param (str, optional): The parameter to create multiple curves. Defaults to None.

    Returns:
    - None
    """

    #begining of user defined params
    #plot_params={'marker':'s','markersize':3,'linewidth':0.7,'linestyle':'-','layout':'h','ylims':'y','ylim1':[-10,-2],'ylim2':[-5,1]}
    plot_params={'marker':'s','markersize':3,'linewidth':0.7,'linestyle':'-','layout':'h','ylims':'n','ylim1':[-8.6,-7.8],'ylim2':[-0.5,0.1]}
    #end of user defined params
    a=['J','T','parity','hw','srg','lambda','theta','NNN','3N_induced','Nf','srg_f','N_max']
    colors=['blue','orange','green','red','purple','brown']
    a.remove(par_x)
    if(extra_param!=''):
        a.remove(extra_param)
    params=ask_data_params(a)
    subdata=data[check_params(params)]
    if(st_num>0):
        subdata=subdata[subdata['N_max']>2]
    n_ex_param=1
    if(extra_param!=''):
        ex_param=subdata[extra_param].unique()
        print("available "+extra_param+"  values:",ex_param)
        s=input("choose the desired values ( all in one line ), or -1 for all")
        ex_param_in=[float(s1) for s1 in s.split()]
        if(ex_param_in[0]==-1):
            ex_param_in=ex_param
        ex_param_in=np.sort(ex_param_in)
        n_ex_param=len(ex_param_in)
        print(" working with: ",n_ex_param,ex_param_in)

    par_x_values=subdata[par_x].unique()
    print("available "+par_x+"  values:",par_x_values)
    s=input("choose the desired values ( all in one line ), or -1 for all")
    par_x_values_in=[float(s1) for s1 in s.split()]
    if(par_x_values_in[0]==-1):
        par_x_values_in=par_x_values
    par_x_values_in=np.sort(par_x_values_in)
    n_par_x=len(par_x_values_in)
    print(" working with: ",n_par_x,par_x_values_in)

    #n=subdata[par_x].nunique()
    n=n_par_x
    if(plot_params['layout']=='h'):
        nx=1
        ny=2
    else:
        nx=2
        ny=1
    figsize=fig_size(nx,ny)
    fig,ax=plt.subplots(nx,ny,figsize=figsize,constrained_layout=False)
    if par_x in error_style:
        requested_error_style=error_style[par_x] # determines the type of error estimation
    else:
        requested_error_style=None
    #plt.subplots_adjust(wspace=0.2, hspace=0.2)
    plt.rcParams.update({'font.size': 12})
    par_x_values=par_x_values_in
    for i in range(n_ex_param):
        if(extra_param!=''):
            mask= subdata[extra_param]==ex_param_in[i]
            label=symbol[extra_param]+"$ ={}$".format(ex_param_in[i])
            #label=r"$e^{"+"{}".format(int(ex_param_in[i]))+r"i\theta}$"
        else:
            mask= subdata['J']>-5
            label=""
        x_values=np.sort(subdata[mask][par_x].unique())
        x_values=par_x_values
        y_values=[subdata[(subdata[par_x]==x_value) &
                          (mask)]['E_r'].iloc[st_num] for x_value in x_values]
        #y_values=y_values-y_shift[0]

        if(requested_error_style=='std'):
            y_avg=np.sum(y_values)/n
            y_std=np.sqrt(np.sum((y_values-y_avg)**2)/n)
            print("value of the real component vs {} : {} +-{}:".format(par_x,y_avg,y_std))
        elif(requested_error_style=='band'):
            y_avg=np.sum(y_values)/n
            y_band=(np.max(y_values)-np.min(y_values))/2
            print("value of the real component vs {} : {}+-{}:".format(par_x,y_avg,y_band))
        ax[0].plot(x_values,y_values,marker=plot_params['marker'],markersize=plot_params['markersize'],linestyle=plot_params['linestyle'],linewidth=plot_params['linewidth'],label=label)
        if(i==0 or n_ex_param<3):
            #ax[0].axhline(y=y_values[-1],color='k',linewidth=1,linestyle=':')
            #ax[0].axhline(-2.224225,color='k',linewidth=1,linestyle=':')
            print(y_values[0])
        #ax[0].tick_params(axis='both', which='major', labelsize=22)
        y2_values=[subdata[(subdata[par_x]==x_value)&(mask)]['E_i'].iloc[st_num]
                  for x_value in x_values]
        #y2_values=y2_values-y_shift[1]
        if(requested_error_style=='std'):
            y_avg=np.sum(y2_values)/n
            y_std=np.sqrt(np.sum((y2_values-y_avg)**2)/n)
            print("value of the imaginary component vs {} : {} +- {}:".format(par_x,y_avg,y_std))
        elif(requested_error_style=='band'):
            y_avg=np.sum(y2_values)/n
            y_band=(np.max(y2_values)-np.min(y2_values))/2
            print("value of the real component vs {} : {}+-{}:".format(par_x,y_avg,y_band))
        ax[1].plot(x_values,y2_values,marker=plot_params['marker'],markersize=plot_params['markersize'],linestyle=plot_params['linestyle'],linewidth=plot_params['linewidth'],label=label)

        #plot scattering state curve
        if(par_x=='theta'):
            #sm=2.1 0^+;0
            #sm=1.9 #0^-;0
            #sm=1.2 #1^-;0
            x_line=np.linspace(0,0.5,100)
            sm=1.3#slope modifier
            radius=-y2_values[0]/np.sin(sm*2*x_values[0])
            threshold=y_values[0]-radius*np.cos(sm*2*x_values[0])
            c_y_values=[threshold+radius*np.cos(sm*2*x_value) for x_value in
                        x_line]
            c_y2_values=[-radius*np.sin(sm*2*x_value) for x_value in
                         x_line]
            ax[0].plot(x_line,c_y_values,linestyle=':',linewidth=plot_params['linewidth']*1.5,color='r')
            ax[1].plot(x_line,c_y2_values,linestyle=':',linewidth=plot_params['linewidth']*1.5,color='r')

        #

        if(i==0 or n_ex_param<3):
            #ax[1].axhline(y=y2_values[-1],color='k',linewidth=1,linestyle=':')
            #ax[1].axhline(y=-0,color='k',linewidth=1,linestyle=':')
            print(y2_values[0])
    #ax[0].set_ylim(-26.5,-24.5)
    #ax[0].set_ylim(-8.25,-7.25)
    #ax[1].set_ylim(-0.9,0.1)
    # add grey box for converged value
    y_min, y_max = ax[0].get_ylim()
    new_label=format_tick_labels([x_values[-1]])[0]
    #ax[0].text(N_values[0]+0.5,x_values[-1]-0.05*(y_max-y_min),new_label,fontsize=12,color='grey', ha='center',va='center',bbox=dict(facecolor='white', edgecolor='grey',alpha=0.5, boxstyle='round,pad=0.5'))
    y_min, y_max = ax[1].get_ylim()
    new_label=format_tick_labels([y_values[-1]])[0]
    #ax[1].text(N_values[0]+0.5,y_values[-1]-0.05*(y_max-y_min),new_label,fontsize=12,color='grey', ha='center',va='center',bbox=dict(facecolor='white', edgecolor='grey',alpha=0.5, boxstyle='round,pad=0.5'))
    ##
    ##
    if(extra_param!=''):
        ax[0].legend(fontsize=9,labelspacing=0.2,loc='lower right')
        #ax[1].legend(fontsize=9,labelspacing=0.2)#,loc='upper right')
    E_name=r"\mathrm{{E^{{{0:}^{{{1:}}}{2:}}}_{{{3:d}}}}}".format(int(params['J'])
            if (A%2==0)else r"\frac{"+str(int(params['J']))+"}{2}",'+' if(params['parity']==1) else'-',int(params['T'])
            if (A%2==0) else r"\frac{"+str(int(params['T']))+"}{2}",int(st_num))
    x_label=symbol[par_x]+( (" $[$"+units[par_x]+"$]$") if par_x in units else "")
    ax[0].set_xlabel(x_label,fontsize=12)
    #ax[0].set_xlabel(symbol[par_x]+r"$\mathrm{\ [rad]}$",fontsize=12)
    ax[0].set_ylabel(r"$\mathcal{R}e\left("+E_name+r"\right)[\SI[]{}{\mega\eV}]$",fontsize=12)
    ax[1].set_xlabel(x_label,fontsize=12)
    #ax[1].set_xlabel(symbol[par_x]+r"$\mathrm{\ [rad]}$",fontsize=12)
    ax[1].set_ylabel(r"$\mathcal{I}m\left("+E_name+r"\right)\ [\SI[]{}{\mega\eV}]$",fontsize=12)
    #logscale
    #linthresh=100
    #ax[0].set_yscale('symlog', linthresh=linthresh/10)
    #ax[1].set_yscale('symlog', linthresh=linthresh/10)
    #import matplotlib.ticker as ticker
    #ax[1].yaxis.set_major_locator(ticker.SymmetricalLogLocator(base=100,linthresh=linthresh, subs=None))
    #ax[0].yaxis.set_major_locator(ticker.SymmetricalLogLocator(base=100,linthresh=linthresh, subs=None))
    #ax[0].yaxis.set_minor_locator(ticker.NullLocator())
    #ax[1].yaxis.set_minor_locator(ticker.NullLocator())

    ax[1].tick_params(axis='both', which='major', labelsize=11)
    ax[0].tick_params(axis='both', which='major', labelsize=11)
    ax[0].locator_params(axis='y', nbins=4)
    ax[1].locator_params(axis='y', nbins=2)
    ax[0].minorticks_on()
    ax[1].minorticks_on()

    ax[0].locator_params(axis='x', nbins=5)
    ax[1].locator_params(axis='x', nbins=5)
    ax[0].xaxis.set_minor_locator(AutoMinorLocator(2))
    ax[1].xaxis.set_minor_locator(AutoMinorLocator(2))

    if(plot_params['ylims']=='y'):
        y1=plot_params['ylim1']
        y2=plot_params['ylim2']
        ax[0].set_ylim(y1[0],y1[1])
        ax[1].set_ylim(y2[0],y2[1])
    #ax[0].set_xlim(0,0.4)
    #ax[1].set_xlim(0,0.4)
    #ax[0].set_ylim(-2000,0)

    #
    #fig.subplots_adjust(left=0.13, bottom=0.23)
    plt.tight_layout()
    plt.show()
    extension=create_name_extension(params)
    filename='gen_plot_'+name+'_'+par_x+'_'+extra_param+curr_nucleus+extension+"_{:d}".format(st_num)
    s=input("do you want to save the figure?y/n")
    if(s=='y'):
        save_figure(fig,'graphs/general_plot/',filename)
    if(requested_error_style=='convergence'):
        #initializing
        figsize=fig_size(nx,ny)
        fig,ax=plt.subplots(nx,ny,figsize=figsize)
        def curv_resid_N(a,x,y):
            weights=[i for i in range(1,len(x)+1)]
            return (a[0]+a[1]*np.exp(-a[2]*x)-y)**2*weights
        def curv_resid_th_r(a,x,y):
            weights=[i*2 for i in range(1,len(x)+1)]
            #return (a[0]+(threshold+radius*np.cos(2*x)-a[0])*np.exp(-(x/a[1])**a[2])-y)**2*weights
            return (-sm*a[2]*2*radius*np.sin(sm*a[2]*2*x)*np.exp(-(x/a[0])**a[1])-y)**2*weights
        def curv_resid_th_i(a,x,y):
            weights=[i*2 for i in range(1,len(x)+1)]
            return (-sm*a[2]*2*radius*np.cos(sm*a[2]*2*x)*np.exp(-(x/a[0])**a[1])-y)**2*weights
        def curv_func_N(x,a,b,c):
            return a+b*np.exp(-c*x)
        def curv_func_th_r(x,a,b,c):
            #return a+(threshold+radius*np.cos(2*x)-a)*np.exp(-(x/b)**c)
            return -sm*c*2*radius*np.sin(sm*c*2*x)*np.exp(-(x/a)**b)
        def curv_func_th_i(x,a,b,c):
            return -sm*c*2*radius*np.cos(sm*c*2*x)*np.exp(-(x/a)**b)
        def func_integrate(y0,dy,dx):
            y_new=dy.copy()
            y_new[0]=y0
            for i in range(1,len(dy)):
                y_new[i]=y_new[i-1]+dy[i-1]*dx[i-1]
            return y_new
        initial_params_N=np.array([1,10,0.1])
        bounds_N=([-40,-100,0.05],[10,100,1])
        initial_params_th=np.array([0.3,7,1])
        bounds_th=([0.01,5,0.9],[0.9,9,1.1])
        #
        for i in range(n_ex_param):
            #initialize
            if(extra_param!=''):
                mask= subdata[extra_param]==ex_param_in[i]
                label=symbol[extra_param]+"$ ={}$".format(ex_param_in[i])
            else:
                mask= subdata['J']>-5
                label=""
            #x_values=np.sort(subdata[mask][par_x].unique())
            x_values=par_x_values
            x_values_avg=(x_values[:-1]+x_values[1:])/2
            y_r_values=np.array([subdata[(subdata[par_x]==x_value)&(mask)]['E_r'].iloc[st_num] for
                      x_value in x_values])
            y_i_values=np.array([subdata[(subdata[par_x]==x_value)&(mask)]['E_i'].iloc[st_num] for
                      x_value in x_values])

           #
            if(par_x=='theta'):
                num_points=6
                curv_resid=curv_resid_th_r
                curv_func=curv_func_th_r
                x0=initial_params_th.copy()
                bounds=tuple(bounds_th)
            #extrapolate the real component
                x_full=np.linspace(x_values[n-4],0.9,100)
                x_full=np.concatenate((x_full,x_values[n-3:n]))
                x_full=np.sort(x_full)
                x_full=np.unique(x_full)
                y_r_derivatives=np.diff(y_r_values)/np.diff(x_values)
                least_r=least_squares(curv_resid,x0,loss='soft_l1',f_scale=0.1,bounds=bounds,args=(x_values_avg[n-num_points:],y_r_derivatives[n-num_points:]))
                #least1=least_squares(curv_resid,x0,loss='soft_l1',f_scale=0.1,bounds=bounds,args=(x_values[n-num_points:],y1_values[n-num_points:]))

                y_r_der_pred=curv_func(x_full[:],*least_r.x)
                #y_r_der_pred_compare=curv_func(x_values,*least_r.x)
                y_r_pred=func_integrate(y_r_values[n-4],y_r_der_pred[:],np.diff(x_full)[:])
                #y_r_pred_compare=func_integrate(y_r_values[n-num_points],y_r_der_pred_compare[n-num_points:],np.diff(x_values)[n-num_points:])
                #y_r_pred_compare=np.concatenate(([0]*(n-num_points),y_r_pred_compare))
                y_r_pred_compare=y_r_values.copy()
                indices = np.isin(x_full, x_values)
                indices_r = np.isin(x_values, x_full)
                print(indices_r)
                y_r_pred_compare[indices_r]=y_r_pred[indices]
            # extrapolate the imaginary component
                x0=initial_params_th.copy()
                bounds=tuple(bounds_th)
                bounds=(bounds[0],[least_r.x[0]*1.1,bounds[1][1],bounds[1][2]])
                curv_resid=curv_resid_th_i
                curv_func=curv_func_th_i

                y_i_derivatives=np.diff(y_i_values)/np.diff(x_values)
                least_i=least_squares(curv_resid,x0,loss='soft_l1',f_scale=0.1,bounds=bounds,args=(x_values_avg[n-num_points:],y_i_derivatives[n-num_points:]))
                y_i_der_pred=curv_func(x_full[:],*least_i.x)
                y_i_pred=func_integrate(y_i_values[n-4],y_i_der_pred[:],np.diff(x_full)[:])

                y_i_pred_compare=y_i_values.copy()
                y_i_pred_compare[indices_r]=y_i_pred[indices]

                error_bar_shift=0.2+i*0.1
                E_r_inf=y_r_pred[-1]
                E_i_inf=y_i_pred[-1]
            else:
                num_points=6
                x0=initial_params_N.copy()
                bounds=bounds_N
                curv_resid=curv_resid_N
                curv_func=curv_func_N

                st_n_max=1 if params['parity']==-1 else 2
                x_full=np.arange(st_n_max ,st_n_max+24,2)
                x0[0]=y_r_values[-1]
                print(x0)
                least_r=least_squares(curv_resid,x0,loss='soft_l1',f_scale=0.1,bounds=bounds,args=(x_values[n-num_points:],y_r_values[n-num_points:]))

                y_r_pred=curv_func(x_full[:],*least_r.x)
                y_r_pred_compare=curv_func(x_values[:],*least_r.x)
                y_r_pred_compare[:n-4]=y_r_values[:n-4]
            # extrapolate the imaginary component
                num_points=6
                x0=initial_params_N.copy()
                bounds=bounds_N
                curv_resid=curv_resid_N
                curv_func=curv_func_N

                x0[0]=y_i_values[-1]
                x0[1]=-x0[1]
                print(x0)
                least_i=least_squares(curv_resid,x0,loss='soft_l1',f_scale=0.1,bounds=bounds,args=(x_values[n-num_points:],y_i_values[n-num_points:]))
                y_i_pred=curv_func(x_full[:],*least_i.x)
                y_i_pred_compare=curv_func(x_values[:],*least_i.x)
                y_i_pred_compare[:n-4]=y_i_values[:n-4]
                error_bar_shift=i+4
                E_r_inf=least_r.x[0]
                E_i_inf=least_i.x[0]
        #plot the results
            label1=label#+"\n $E_{{\inf}}={:.2f}$".format(least1.x[0])
            ax[0].plot(x_full[2:],y_r_pred[2:],markersize=10,linestyle=':',label=label1,color=colors[i])
            ax[0].plot(x_values,y_r_values,marker=plot_params['marker'],markersize=plot_params['markersize'],linestyle=plot_params['linestyle'],linewidth=plot_params['linewidth'],color=colors[i])


            label2=label#+"\n $E_{{\inf}}={:.2f}$".format(least2.x[0])
            ax[1].plot(x_full[:],y_i_pred[:],markersize=10,linestyle=':',label=label2,color=colors[i])
            ax[1].plot(x_values,y_i_values,marker=plot_params['marker'],markersize=plot_params['markersize'],linestyle=plot_params['linestyle'],linewidth=plot_params['linewidth'],color=colors[i])
            #ax[1].plot(x_full[:],y2_der_pred[:],markersize=10,linestyle=':',label=label2,color='r')
            #ax[1].plot(x_values[:n-1],y2_derivatives,marker='x',markersize=10,color=colors[i])
            residual1=y_r_pred_compare[n-num_points:n]-y_r_values[n-num_points:]
            residual2=y_i_pred_compare[n-num_points:n]-y_i_values[n-num_points:]
            d1_error=np.sqrt(np.abs(residual1**2).sum()/num_points)
            d2_error=np.sqrt(np.abs(residual2**2).sum()/num_points)
            E_r_conv=y_r_pred_compare[n-1]-E_r_inf
            E_i_conv=y_i_values[-1]-E_i_inf
            E_r_err=np.sqrt(d1_error**2+E_r_conv**2)
            E_i_err=np.sqrt(d2_error**2+E_i_conv**2)
            #printing results:
            print("******\nfit results for {}:".format(label))
            print(residual1,residual2)
            print("fit params:",least_r.x[:],least_i.x[:])
            print("E_inf:", y_r_pred[-1],y_i_pred[-1])
            print("final value:",y_r_values[-1],y_i_values[-1])
            print("convergence in E_r:",d1_error,E_r_conv,E_r_err)
            #print(y1_values[-1]-y1_values[-2],y1_values[-2]-y1_values[-3])
            print("convergence in E_i:",d2_error,E_i_conv,E_i_err)
            #print(y2_values[-1]-y2_values[-2],y2_values[-2]-y2_values[-3])
            print("convergence in E_r:",E_r_conv)
            print("convergence in E_i:",E_i_conv)
            ax[0].errorbar(x_values[-1]+error_bar_shift,E_r_inf,yerr=E_r_err,fmt='o',capsize=4+i,color=colors[i],markersize=2)
            ax[1].errorbar(x_values[-1]+error_bar_shift,E_i_inf,yerr=E_i_err,fmt='o',capsize=4+i,color=colors[i],markersize=2)
        #ax[0].set_ylim(-7.5,-5.5)
        #ax[1].set_ylim(-4,0.5)
        ax[0].set_xlabel(symbol[par_x],fontsize=12)
        ax[0].set_ylabel(r"$\mathcal{R}e\left("+E_name+r"\right)[\SI[]{}{\mega\eV}]$",fontsize=12)
        ax[1].set_xlabel(symbol[par_x],fontsize=12)
        ax[1].set_ylabel(r"$\mathcal{I}m\left("+E_name+r"\right)\ [\SI[]{}{\mega\eV}]$",fontsize=12)
        ax[1].tick_params(axis='both', which='major', labelsize=11)
        ax[0].tick_params(axis='both', which='major', labelsize=11)
        ax[0].locator_params(axis='x', nbins=5)
        ax[1].locator_params(axis='x', nbins=5)
        ax[0].minorticks_on()
        ax[1].minorticks_on()
        ax[0].xaxis.set_minor_locator(AutoMinorLocator(2))
        ax[1].xaxis.set_minor_locator(AutoMinorLocator(2))
        if(extra_param!=''):
            ax[0].legend(fontsize=9,labelspacing=0.2)#,loc='upper right')
        plt.tight_layout()
        plt.show()
        filename='gen_plot_extrapolated_'+name+'_'+par_x+'_'+extra_param+curr_nucleus+extension+"_{:d}".format(st_num)
        s=input("do you want to save the figure?y/n")
        if(s=='y'):
            save_figure(fig,'graphs/general_plot/',filename)


def resonance_hwconv(extra_param='',y_shift=[0,0],name=''):
    """
    Plots a set of curves varied by hw for the convergence of the energy 
    vs N_max for multiple states. To specify which states to plot, 
    modifications need to be made directly in the code.

    Parameters:
    - None

    Returns:
    - None
    """
    #states=[[0,0,1,1],[0,0,-1,0],[1,0,-1,0],[2,0,-1,0],[1,1,-1,0],[1,1,-1,1],[2,1,-1,0]]
    states=[[0,0,1,8],[0,0,-1,8],[1,0,-1,8]]
    par_x='N_max'
    #begining of user defined params
    #plot_params={'marker':'s','markersize':3,'linewidth':0.7,'linestyle':'-','layout':'h','ylims':'y','ylim1':[-10,-2],'ylim2':[-5,1]}
    plot_params={'marker':'o','markersize':2,'linewidth':0.7,'linestyle':'-','layout':'h','ylims':'y','ylim1':[-18,1],'ylim2':[-18,1]}
    #end of user defined params
    a=['hw','srg','lambda','theta','NNN','3N_induced','Nf','srg_f']
    if par_x in a:
        a.remove(par_x)
    if(extra_param!=''):
        a.remove(extra_param)
    params=ask_data_params(a)
    subdata=data[check_params(params)]
    subdata=subdata[subdata['N_max']>3]
    n_ex_param=1
    if(extra_param!=''):
        ex_param=subdata[extra_param].unique()
        print("available "+extra_param+"  values:",ex_param)
        s=input("choose the desired values ( all in one line ), or -1 for all")
        ex_param_in=[float(s1) for s1 in s.split()]
        if(ex_param_in[0]==-1):
            ex_param_in=ex_param
        ex_param_in=np.sort(ex_param_in)
        n_ex_param=len(ex_param_in)
        print(" working with: ",n_ex_param,ex_param_in)

    #par_x_values=subdata[par_x].unique()
    #print("available "+par_x+"  values:",par_x_values)
    #s=input("choose the desired values ( all in one line ), or -1 for all")
    #par_x_values_in=[float(s1) for s1 in s.split()]
    #if(par_x_values_in[0]==-1):
    #    par_x_values_in=par_x_values
    #par_x_values_in=np.sort(par_x_values_in)
    #n_par_x=len(par_x_values_in)
    #print(" working with: ",n_par_x,par_x_values_in)

    #n=subdata[par_x].nunique()
    #n=n_par_x
    if(plot_params['layout']=='h'):
        nx=1
        ny=2
    else:
        nx=2
        ny=1
    figsize=fig_size(nx,ny)
    figsize=(figsize[0],figsize[1]*2)
    fig,ax=plt.subplots(nx,ny,figsize=figsize,constrained_layout=False)
    requested_error_style=error_style[par_x] # determines the type of error estimation
    #plt.subplots_adjust(wspace=0.2, hspace=0.2)
    plt.rcParams.update({'font.size': 12})
    colors=['blue','orange','green','red','purple','brown']
    shift1=np.zeros(len(states))
    shift2=np.zeros(len(states))
    for i in range(n_ex_param):
        if(extra_param!=''):
            mask1= subdata[extra_param]==ex_param_in[i]
            label=symbol[extra_param]+"$ ={}$".format(ex_param_in[i])
        else:
            mask1= subdata['J']>-5
            label=""
        for j,state in enumerate(states):
            full_params=params.copy()
            full_params['J']=state[0]
            full_params['T']=state[1]
            full_params['parity']=state[2]
            #mask2=check_params({**params,'J':state[0],'T':state[1],'parity':state[2]})
            mask2=check_params(full_params,subdata)
            st_num=state[3]
            mask=mask1& mask2
            x_values=np.sort(subdata[mask][par_x].unique())
            x_values=x_values[x_values<=20]
            y_values=[subdata[(subdata[par_x]==x_value) &
                              (mask)]['E_r'].iloc[st_num] for x_value in x_values]
            label=None
            if(i==0):
                shift1[j]=-12-j*2-y_values[-1]
                shift1[j]=-7-j*5-y_values[-1]
                label=r"$\mathrm{{E^{{{0:}^{{{1:}}}{2:}}}}}$".format(int(state[0])
            if (A%2==0)else r"\frac{"+str(int(state[0]))+"}{2}",'+'
                                                                   if(state[2]==1)
                                                                   else'-',int(state[1])
            if (A%2==0) else r"\frac{"+str(int(state[1]))+"}{2}")
    #
            y_values=[ value+shift1[j] for value in y_values]
            ax[0].plot(x_values,y_values,marker=plot_params['marker'],markersize=plot_params['markersize'],linestyle=plot_params['linestyle'],linewidth=plot_params['linewidth'],label=label,color=colors[j%6])
            y2_values=[subdata[(subdata[par_x]==x_value)&(mask)]['E_i'].iloc[st_num]
                        for x_value in x_values]
            if(i==0):
                shift2[j]=-1-j*2-y2_values[-1]
                shift2[j]=-2-j*5-y2_values[-1]
            y2_values=[ value+shift2[j] for value in y2_values]
            #y2_values=y2_values-y_shift[1]
            ax[1].plot(x_values,y2_values,marker=plot_params['marker'],markersize=plot_params['markersize'],linestyle=plot_params['linestyle'],linewidth=plot_params['linewidth'],color=colors[j%6])

    # add grey box for converged value
    y_min, y_max = ax[0].get_ylim()
    new_label=format_tick_labels([x_values[-1]])[0]
    #ax[0].text(N_values[0]+0.5,x_values[-1]-0.05*(y_max-y_min),new_label,fontsize=12,color='grey', ha='center',va='center',bbox=dict(facecolor='white', edgecolor='grey',alpha=0.5, boxstyle='round,pad=0.5'))
    y_min, y_max = ax[1].get_ylim()
    new_label=format_tick_labels([y_values[-1]])[0]
    #ax[1].text(N_values[0]+0.5,y_values[-1]-0.05*(y_max-y_min),new_label,fontsize=12,color='grey', ha='center',va='center',bbox=dict(facecolor='white', edgecolor='grey',alpha=0.5, boxstyle='round,pad=0.5'))
    ##
    ##
    #if(extra_param!=''):
    #    ax[0].legend(fontsize=9,labelspacing=0.2)#,loc='upper right')
        #ax[1].legend(fontsize=9,labelspacing=0.2)#,loc='upper right')
    E_name=r"\mathrm{E}"
    ax[0].set_xlabel(symbol[par_x],fontsize=12)
    ax[0].set_ylabel(r"$\mathcal{R}e\left("+E_name+r"\right)[\SI[]{}{\mega\eV}]$",fontsize=12)
    ax[1].set_xlabel(symbol[par_x],fontsize=12)
    ax[1].set_ylabel(r"$\mathcal{I}m\left("+E_name+r"\right)\ [\SI[]{}{\mega\eV}]$",fontsize=12)
    ax[1].tick_params(axis='both', which='major', labelsize=11)
    ax[0].tick_params(axis='both', which='major', labelsize=11)
    ax[0].locator_params(axis='x', nbins=5)
    ax[1].locator_params(axis='x', nbins=5)
    ax[0].minorticks_on()
    ax[1].minorticks_on()
    ax[0].legend(loc='lower left', fontsize=8,markerscale=0.7)
    #
    if(plot_params['ylims']=='y'):
        y1=plot_params['ylim1']
        y2=plot_params['ylim2']
        ax[0].set_ylim(y1[0],y1[1])
        ax[1].set_ylim(y2[0],y2[1])
    #
    ax[0].xaxis.set_minor_locator(AutoMinorLocator(2))
    ax[1].xaxis.set_minor_locator(AutoMinorLocator(2))
    plt.tight_layout()
    extension=create_name_extension(params)
    filename='resonance_hwconv_'+name+'_'+par_x+'_'+extra_param+curr_nucleus+extension+"_{:d}".format(st_num)
    fig.subplots_adjust(left=0.13, bottom=0.23)
    plt.show()
    s=input("do you want to save the figure?y/n")
    if(s=='y'):
        save_figure(fig,'graphs/general_plot/',filename)
#---------------------------------------------------------------------------
#plots the spectrum in the complex plan, each theta has its own sub-figure.
#multiple hw values are shown on each sub-figure. the thresholds are estimated
#from the variation with hw
def  plot_complex_plane(num_p,first_param='',second_param='',name=''):
    """
    Plots multiple subplots of eigenvalues in the complex plane, 
    allowing visualization of different parameters.

    Parameters:
    - num_p (int): The number of states to plot, starting from the lowest one (real component).
    - first_param (str): Defines the sets of data to be plotted in the same subplot (e.g., 'hw').
    - second_param (str): Defines the number of subplots, plotting data for various values of a parameter (e.g., 'lambda').

    Returns:
    - None
    """
    #params=ask_data_params(['J','T','parity','hw','N_max','srg','lambda','NNN','3N_induced','Nf','srg_f'])
    a=['J','T','parity','hw','N_max','srg','lambda','NNN','3N_induced','Nf','srg_f','theta']
    if(first_param!=''):
        a.remove(first_param)
    if(second_param!=''):
        a.remove(second_param)
    params=ask_data_params(a)
    subdata=data[check_params(params)]

    first_param_in,n_first_param=ask_for_extra_param(first_param,subdata)

    thresholds=set_threshold(params)
    second_param_in,n=ask_for_extra_param(second_param,subdata)
    ny=min(2,n)
    nx=ceil(n/2)
    figsize=fig_size(nx,ny,0.7 if ny==1 else 1)
    fig,axs=plt.subplots(nx,ny,figsize=figsize)
    if not isinstance(axs,np.ndarray):
        axs=np.array([[axs]])
    for i in range(n,nx*ny):
        fig.delaxes(axs[i//ny,i%ny])
    plt.rcParams.update({'font.size': 12})
    color=['r','g','b','y','c','tab:brown','m']
    ind=0
    EE=[-8.49,-9.98,-11.5,-13.2,-14.9,-16.7,-18.6]
    EE=[0]*10
    print(axs.shape)
    # indices for the range of states to be plotted
    pi=0
    pe=num_p
    #
    state_list=[6,3,2,1,1,1,1,1,1]
    for l, ax in enumerate(axs.flat):
    #for l in range(len(thetas_in)):
        if(l==n):
            break
        # extra params
        for i in range(n_first_param):
            if(first_param!=''):
                mask= subdata[first_param]==first_param_in[i]
                label=symbol[first_param]+"$ ={}$".format(first_param_in[i])
            else:
                mask= subdata['J']>-5
                label=""
        #

            if(subdata[(mask)&
                       (subdata[second_param]==second_param_in[l])].empty):
                continue
            E_r=subdata[(mask)&
                        (subdata[second_param]==second_param_in[l])]['E_r']
            #
            E_r=E_r-2*EE[l]
            print(EE[l])
            #
            E_i=subdata[(mask)&
                        (subdata[second_param]==second_param_in[l])]['E_i']
            # comment the following to print a default range
            #pi=state_list[i]-1
            #pe=state_list[i]
            #
            pe=min(pe,len(E_r))
            ax.scatter(E_r[pi:pe],E_i[pi:pe],label=label,linewidths=1,clip_on=True,marker='+',s=40)

            thresholds={curr_nucleus:[0,-EE[l],-2*EE[l]]}
            if('theta' in params):
                theta=params['theta']
            elif (first_param=='theta'):
                theta=first_param_in[i]
            else:
                theta=second_param_in[l]
            for i,j in enumerate(thresholds[curr_nucleus]):
                xx=np.linspace(0,40,40)
                yy=-xx*math.tan(2*float(theta))
                xx[:]=xx[:]+j
                ax.plot(xx,yy,color='k',linestyle=':',linewidth=0.8)
        axj=ind%nx
        axi=int(ind/nx)
        ax.set_xlim(max(min(subdata['E_r'].min(),thresholds[curr_nucleus][0])-1,-40),12)
        ax.set_ylim(-10,2)
        ax.set_xlim(-3,10)
        #ax.set_xlim(-3,10)
        #ax.set_ylim(-10,2)
        #ax.set_xlim(-3-l*1.5,10-l*1.5)
        #ax.set_ylim(-12,2)
        #ax.set_xlim(-4,24)
        #ax.set_title(r'$2$'+symbol['theta']+r"$={:.2f}^\circ$".format(2*thetas_in[l]*180.0/math.pi),fontsize=10)
        ax.set_title(symbol[second_param]+r"$={:.2f}$".format(second_param_in[l])+" "+units[second_param],fontsize=10)
        ##########
        #thresholds={curr_nucleus:[-8,-7.3,-4.4,-2.2,0]}
        #thresholds={curr_nucleus:[-7.2,-4.4,-2.29,-0.4,-1.6]}
        #thresholds={curr_nucleus:[0,-EE[l],2*EE[l]]}
        #thresholds={curr_nucleus:[0]}
        thresholds={curr_nucleus:[0,-EE[l],-2*EE[l]]}
        #rotation_slope=[1.6,1.6,1,1,1]
        #rotation_slope=[1.65,1.6,1.35,1.6,1.0886]
        rotation_slope=[1,1,1]
        if('theta' in params):
            theta=params['theta']
        elif (second_param=='theta'):
            theta=second_param_in[l]
        else:
            theta=0.0
        for i,j in enumerate(thresholds[curr_nucleus]):
            xx=np.linspace(0,40,40)
            yy=-xx*math.tan(2*float(theta))*rotation_slope[i]
            xx[:]=xx[:]+j
            #ax.plot(xx,yy,color='k',linestyle='--',linewidth=0.5)
        ##########
        ax.axhline(y=0,color='k')
        ax.axvline(x=0,color='k')
        #ax[axi,axj].set_title("θ={:0.3f}".format(thetas[ind]),y=0,pad=-25,fontdict={'fontsize':14})
        ax.legend(loc='lower left', fontsize=9,markerscale=0.7)
        ax.locator_params(axis='both', nbins=5)
        ax.tick_params(axis='both', which='major', labelsize=11)
        ind+=1
    #fig.subplots_adjust(left=0.13, bottom=0.15)
    fig.subplots_adjust(left=0.20, bottom=0.18,hspace=0.35)
    E_name=r"\mathrm{{E^{{{0:}^{{{1:}}}{2:}}}}}".format(int(params['J'])
            if (A%2==0)else r"\frac{"+str(int(params['J']))+"}{2}",'+' if(params['parity']==1) else'-',int(params['T'])
            if (A%2==0) else r"\frac{"+str(int(params['T']))+"}{2}")
    #
    fig.text(0.5,0.03, r"$\mathcal{R}e\left("+E_name+r"\right)[\SI[]{}{\mega\eV}]$",ha='center',va='center')
    fig.text(0.04,0.5,r"$\mathcal{I}m\left("+E_name+r"\right)\ [\SI[]{}{\mega\eV}]$",ha='center',va='center',rotation='vertical')
    #rect = plt.Rectangle((-2.9, -1), 2.5, 2,
    #                     edgecolor='black',linestyle='dashed', linewidth=0.2, fill=False)
    #ax.add_patch(rect)
    #ax.set_xlabel("Real(E) [Mev]",fontsize=22)
    #ax.set_ylabel("Imaginary(E) [Mev]",fontsize=22)
    #ax.tick_params(axis='both', which='major', labelsize=20)
    #plt.tight_layout()
    plt.show()
    extension=create_name_extension(params)
    s=input("do you want to save the figure?y/n")
    if(s=='y'):
        save_figure(fig,'graphs/rotation/',name+curr_nucleus+extension)

def plot_example(theta=np.pi/4, roots=None, filename='momentum_plane.svg'):
    if roots is None:
        roots = []

    # Lines from origin
    x_line = np.linspace(0, 20, 200)
    y_real = np.zeros_like(x_line)
    y_slope = -np.tan(theta) * x_line

    plt.rcParams.update({'font.size': 12})
    # Create plot
    figsize=fig_size(1,1,0.7)
    fig,ax=plt.subplots(1,1,figsize=figsize)
    ax.axhline(y=0,color='k')
    ax.axvline(x=0,color='k')
    ax.plot(x_line, y_slope,
            label=r'$\mathrm{Integration\ path}$', color='red', linestyle='--')

    # Plot roots
    if roots:
        ax.scatter(roots, np.zeros_like(roots),
                   marker='o',color='b',label='$\mathrm{Roots}$',s=10)

    # Origin

    # Axes settings
    ax.set_xlabel(r'$\mathcal{R}e(p)$')
    ax.set_ylabel(r'$\mathcal{I}m(p)$')
    ax.legend()
    ax.set_title(r'$\mathrm{Momentum\ Complex\ Plane}$')
    ax.set_ylim(-8,3)
    ax.set_xlim(-4,11)

    # Save as SVG
    plt.tight_layout()
    plt.show()
    save_figure(fig,'graphs/rotation/',filename)
    plt.close(fig)


### varius plotings for my thesis (not related to spectra)
def plot_resonance(E_R=4.0, sigma_0=2, Gamma=0.2, E_min=0.0, E_max=10.0, num_points=500):
    E = np.linspace(E_min, E_max, num_points)
    sigma =(sigma_0 * Gamma**2) / ((E - E_R)**2 + Gamma**2)
    conti=2*np.exp(-(E)/4)+0.5
    noise=np.random.normal(0,0.03*np.max(sigma),size=sigma.shape)
    figsize=fig_size(1,1,0.7)
    fig,ax=plt.subplots(1,1,figsize=figsize)
    ax.plot(E, sigma+conti+noise, label=r'$\sigma(E) = \frac{\sigma_0 \Gamma^2}{(E - E_R)^2 + \Gamma^2}$', color='k')
    ax.set_xlabel(r'$E$', fontsize=14)
    ax.set_ylabel(r'$\sigma(E)$', fontsize=14)
    ax.set_ylim(0,4)
    ax.set_xlim(E_min,E_max)
    plt.tight_layout()
    fig.savefig('resonance_cross_section.svg', format='svg')
    plt.show()

def plot_square_well_wave_function(L=10.0, n=1, x_min=0.0, x_max=20.0, num_points=500):
    x = np.linspace(x_min, x_max, num_points)
    k = n * np.pi / L
    psi = np.zeros_like(x)
    psi=np.exp(-1.3*x)*np.sin(3*x**1.5)
    psi[x>4]=psi[x>4]*np.exp(2.5*(x-4))[x>4]
    figsize=fig_size(1,1,0.7)
    fig,ax=plt.subplots(1,1,figsize=figsize)
    ax.plot(x, psi, label=r'$\psi_n(x)$', color='green')
    ax.axhline(y=0,color='k')
    ax.set_xlabel(r'$x$', fontsize=18)
    ax.set_ylabel(r'$\psi(x)$', fontsize=18)
    ax.set_xlim(x_min,x_max)
    plt.tight_layout()
    fig.savefig('siegert_wave_function.svg', format='svg')
    plt.show()
def plot_tetraneutron_data():
    """
    Plot the data for Figure 4 with proper formatting and LaTeX-compatible fonts.
    """
    # Set up LaTeX-compatible fonts
    plt.rcParams.update({'font.size': 12})

    # Create figure and axis
    figsize=fig_size(1,1,0.7)
    fig,ax=plt.subplots(1,1,figsize=figsize)

    # Plot red full symbol data
    red_full_x = [1.75]
    red_full_dx = [0.37]
    red_full_y = [2.37]
    red_full_dy = [0.58]
    ax.errorbar(red_full_x, red_full_y, xerr=red_full_dx, yerr=red_full_dy,
                fmt='o', color='red', markersize=8, capsize=3, capthick=1,
                label=r'$\mathrm{Duer}\ \mathit{et\ al.}$')

    # Plot red open symbol data
    red_open_x = [2.6]
    red_open_dx = [0]
    red_open_y = [0.83]
    red_open_dy = [1.4]
    ax.errorbar(red_open_x, red_open_y, xerr=red_open_dx, yerr=red_open_dy,
                fmt='o', mfc='none', color='red', markersize=8, capsize=3,
                capthick=1,
                label=r'$\mathrm{Kisamori}\ \mathit{et\ al.}$')

    # Plot full stars data
    full_stars_x = [1.4, 0.85, 1.3]
    full_stars_y = [0.8, 0.3, 0.8]
    ax.scatter(full_stars_x, full_stars_y, marker='*', color='blue', s=200,
                label=r'$\mathrm{NCSM}$')

    # Plot open star data
    open_star_x = [2.38]
    open_star_y = [2.64]
    ax.scatter(open_star_x, open_star_y, marker='*', facecolor='none',
               edgecolor='blue', s=200,
                label=r'$\mathrm{NCGSM}$')

    # Plot cross data
    cross_x = [3.7]
    cross_y = [7.22]
    ax.scatter(cross_x, cross_y, marker='x', color='blue', s=200, linewidth=2,
                label=r'$\mathrm{NCGSM}$')

    # Plot band data
    band_x_min = 0
    band_x_max = 4
    band_y = 2.1
    band_dy = 0.2
    ax.fill_betweenx(
        [band_y - band_dy, band_y + band_dy],  # Y-range of the band
        band_x_min, band_x_max,                # X-range (full width)
        color='gray', alpha=0.3,label=r'$\mathrm{QMC}$')

    # Add labels and title
    ax.set_xlabel(r'$\mathrm{Width\ [MeV]}$', fontsize=14)
    ax.set_ylabel(r'$\mathrm{Energy\ [MeV]}$', fontsize=14)

    # Add legend
    ax.legend(loc='best', fontsize=10, framealpha=1,markerscale=0.5)

    # Adjust layout
    plt.tight_layout()

    plt.show()
    fig.savefig('tetraneutron_predictions.svg', format='svg')
    return fig, ax
if __name__=='__main__':
    help()
