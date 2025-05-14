import time
import numpy as np
import scipy.optimize as opt
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from lmfit import Model,Parameters,minimize
import sys
import inspect
import os
import textwrap
#parameters to modify

prange=int(6.0)
Nf=3
tries=2
interaction='n3lo'
#

ee=np.exp(1)
# timing control
import signal
class TimeoutException(Exception):
    print("Function call stoped due to time limit")
def handler(signum,frame):
    raise TimeoutException()
signal.signal(signal.SIGALRM,handler)

#

def time_string(t):
    t=t//60
    minutes=t%60
    t=t//60
    hour=t%24
    t=t//24
    day=t
    s=""
    if(day>0):
        s=s+" {} days ".format(day)
    if( hour>0):
        s=s+" {} hours  ".format(hour)
    s=s+" {} minutes".format(minutes+1)
    return s
def generate_input(x1,num):
    xy = np.zeros((2, num*num))
    for i in range(num):
        for j in range(num):
            ind=(i)*num+j
            xy[0,ind]=x1[i]
            xy[1,ind]=x1[j]
    return xy
def area(v_n,v1,range):
    summ=np.abs(v1).sum()*(range**2)/(v_n**2)
    return summ

def fromfile(v_n,file_name1,file_name2):
    x1=np.zeros(v_n)
    v1=np.zeros((v_n,v_n))
    i=0
    with open(file_name1) as f1:
        for line in f1:
            for r in line.split():
                x1[i]=float(r)
                i+=1
    i=0
    j=0
    with open(file_name2) as f2:
        for line in f2:
            for r in line.split():
                v1[i,j]=float(r)
                j+=1
            j=0
            i+=1
    xy=generate_input(x1,v_n)
    return x1,xy,v1


def poly_gaus(xy,a,xi,xj,wi,wj,ratio_i):
    x, y = xy
    ci=2/(wi*wi)
    cj=2/(wj*wj)
    ni=ci*xi*xi
    nj=cj*xj*xj
    a0=a
    if(xi!=0):
        a0=a0*((xi*xi/ee)**(-ni/2))
    if(xj!=0):
        a0=a0*((xj*xj/ee)**(-nj/2))
    log1=ci*(x**2*-0.5)+cj*(y**2*-0.5)
    log2=np.log(x)*ni+np.log(y)*nj+log1
    #weighted tail
    #log3=np.log(1+np.exp(x-2)+np.exp(y-2))
    #log2=log2+log3
    #
    #protect agains overflow
    log2=np.maximum(log2,-25)
    log2=np.minimum(log2,30)
    #
    nexp=np.exp(log2)
    return a0*nexp
def regulator(xy,cutoff):
    x, y=xy
    nexp=np.exp(-(x/cutoff)**7 -(y/cutoff)**7)
    return nexp
def initial_params():
    #initial_x=np.linspace(0.1,1,Nf)
    #initial_y=np.linspace(0.1,1,Nf)
    initial_x=np.random.random(Nf)+0.1
    initial_y=np.random.random(Nf)+0.1
    initial_x=initial_x*3
    initial_y=initial_y*3
    initial_w = np.random.random(Nf)*2+1
    return initial_x,initial_y,initial_w

def write_parameters(fit_module,scaling,file_name):
    params=fit_module.params.valuesdict()
    f1=open(file_name,"w")
    f1.write(str(Nf)+"\n")
    for i in range(1,Nf+1):
        s=""
        s+="{0:15.9f}  ".format(params["f{:d}_a".format(i)])
        s+="{0:15.9f}  ".format(params["f{:d}_xi".format(i)])
        s+="{0:15.9f}  ".format(params["f{:d}_xj".format(i)])
        s+="{0:15.9f}  ".format(params["f{:d}_wi".format(i)])
        s+="{0:15.9f}  ".format(params["f{:d}_wj".format(i)])
        s+="\n"
        f1.write(s)
        #print(s)
    f1.write(str(scaling)+"\n")
    f1.close()
def plot(x,xy,v1, fit_model,file_name):
    yi, xi = np.meshgrid(x,x,indexing='xy')#:1:30j, -.2:1.2:30j]
    zpred=cost_function(fit_model.params,x=xy)
    z1=zpred.reshape(v1.shape)
    #v1=weighted_tail(x,v1)
    #z1=weighted_tail(x,z1,reverse=True)
    fig= plt.figure(figsize=(20,20))
    ax = fig.add_subplot(221,projection='3d')
    ax2 = fig.add_subplot(222,projection='3d')
    ax.plot_surface(xi, yi, v1, edgecolor='royalblue', lw=0.5, rstride=8, cstride=2,
                alpha=0.3)
    ax.plot_surface(xi, yi, z1, edgecolor='red', lw=0.5, rstride=8, cstride=2,
                alpha=0.3)
    ax2.plot_surface(xi, yi, z1-v1, edgecolor='royalblue', lw=0.5, rstride=8, cstride=2,
                     alpha=0.3)
    #ax.scatter(x, y, c=zobs, s=200, vmin=zpred.min(), vmax=zpred.max())
    #ax.invert_yaxis()
    #plt.show()
    fig.savefig(file_name+".pdf")
    fig.savefig(file_name+".svg")
    plt.close()
def fig_size(nx,ny,scale=1):
    latex_width=438.17227*scale
    plot_width=1.0/72.27 *latex_width
    aesthetic_ratio=(5.0**0.5-1)/2*1.1
    plot_height=plot_width*aesthetic_ratio*nx/ny
    figsize=(plot_width,plot_height)
    return figsize
def plot2(x,xy,v1, fit_model,file_name):
    plt.rc('text', usetex=True)
    plt.rc('font', **{'family':'sans-serif','sans-serif':['Helvetica']})
    plt.rcParams['text.latex.preamble']=r"\usepackage{amsmath,amssymb,siunitx,gensymb}"
    plt.rcParams.update({'font.size': 12})
    yi, xi = np.meshgrid(x,x,indexing='xy')#:1:30j, -.2:1.2:30j]
    figsize=fig_size(1,1)
    fig=plt.figure(figsize=figsize)
    ax = fig.add_subplot(111,projection='3d')
    ax.plot_surface(xi, yi, v1, edgecolor='royalblue', lw=0.5, rstride=8, cstride=2,
                alpha=0.3)
    ax.set_xlabel(r'$p_i[\mathrm{fm^{-1]}}]$', fontsize=12)
    ax.set_ylabel(r'$p_j[\mathrm{fm^{-1}}]$', fontsize=12)
    ax.set_zlabel('$V(p_i,p_j)$', fontsize=12)
    fig.savefig("test/"+file_name+".pdf")
    fig.savefig("test/"+file_name+".eps")
    plt.close()
def weighted_tail(x1,v1,reverse=False):
    v2=np.exp(x1[:,np.newaxis]-2)+np.exp(x1-2)+1
    if(not reverse):
        return v2*v1
    else:
        return v1/v2
def run_fit(v_n,x1,xy,v1):
    #v1=weighted_tail(x1,v1)
    v_mesh=v1.flatten()
    initial_area=area(v_n,v_mesh,4)
    scale=1
    if(initial_area<1 and initial_area>0):
        scale=1/initial_area
    v_mesh=v_mesh*scale
    initial_area=area(v_n,v_mesh,4)
    print("initial area: ",initial_area,'scaled by:',scale)
    smallest_area=initial_area
    mse_min=100
    pref=["f{0:d}_".format(i) for i in range(0,Nf+1)]
    params=Parameters()
    for i in range(1,Nf+1):
        if(i==1):
            Mod=Model(poly_gaus,prefix=pref[i])
        else:
            Mod=Mod+Model(poly_gaus,prefix=pref[i])
        params.add(name=pref[i]+'a',min=-20,max=20)
        params.add(name=pref[i]+'xi',min=0.1,max=4)
        params.add(name=pref[i]+'xj',min=0.1,max=4)
        #params.add(name=pref[i]+'wi',min=0.2,max=4.5)
        #params.add(name=pref[i]+'wj',min=0.2,max=4.5)
        #params.add(name=pref[i]+'ratio_i',expr="f{:d}_xi/f{:d}_wi".format(i,i),min=0,max=2.4)
        #params.add(name=pref[i]+'ratio_j',expr="f{:d}_xj/f{:d}_wj".format(i,i),min=0,max=2.4)
        params.add(name=pref[i]+'ratio_j',min=0.45,max=20,value=1)
        params.add(name=pref[i]+'ratio_i',min=0.45,max=20,value=1)
        params.add(name=pref[i]+'wi',expr="f{:d}_xi*f{:d}_ratio_i".format(i,i),min=0.1,max=4.5)
        params.add(name=pref[i]+'wj',expr="f{:d}_xj*f{:d}_ratio_j".format(i,i),min=0.1,max=4.5)
    #Mod=Mod*Model(regulator)
    #params.add(name='cutoff',min=5,max=8)
    print('********')
    time_limit=1000
    for i in range(tries):
        initial_x,initial_y,initial_w=initial_params()
        for j in range(1,Nf+1):
            params[pref[j]+'a'].value=5 if j<3 else 1 if j<6 else 0.1
            if(initial_area==0):
                params[pref[j]+'a'].value=0

            params[pref[j]+'wi'].value=initial_w[j-1]
            params[pref[j]+'wj'].value=initial_w[j-1]
            params[pref[j]+'xi'].value=initial_x[j-1]
            params[pref[j]+'xj'].value=initial_y[j-1]
        #params['cutoff'].value=6.7
        print('starting the fitting')
        time1=time.time()
        #result=Mod.fit(v_mesh,params,xy=xy,method='trust-constr')
        try:
            signal.alarm(time_limit)
            signal.alarm(0)
            result=Mod.fit(v_mesh,params,xy=xy,method='least_squares')
            signal.alarm(0)
        except TimeoutException:
            continue
        time2=time.time()
        print("$Time elapsed:",time2-time1)
        if((time2-time1)*2.0 < time_limit):
            time_limit=int((time2-time1)*2.0)

        zpred=result.eval(xy=xy)
        final_area=area(v_n,zpred-v_mesh,4)
        if(final_area<=smallest_area):
            smallest_area=final_area
            best_model=result
        print("current area",final_area)
        #for i in range(1,Nf+1):
        #    xi=best_model.params["f{:d}_xi".format(i)].value
        #    wi=best_model.params["f{:d}_wi".format(i)].value
        #    ratio=best_model.params["f{:d}_ratio_i".format(i)].value
        #    if(xi/wi>2.5):
        #        print("xi={0} wi={1} ratio:{2},{3}".format(xi,wi,xi/wi,ratio))
        #    xi=best_model.params["f{:d}_xj".format(i)]
        #    wi=best_model.params["f{:d}_wj".format(i)]
        #    ratio=best_model.params["f{:d}_ratio_j".format(i)].value
        #    if(xi/wi>2.5):
        #        print("xj={0} wj={1} ratio:{2},{3}".format(xi,wi,xi/wi,ratio))
    print('********')
    initial_parameters=best_model.init_values
    best_params=best_model.best_values
    #print("cutoff:",best_model.best_values['cutoff'])
    zpred=best_model.eval(xy=xy)
    #z1=zpred.reshape(v1.shape)
    final_area=area(v_n,zpred-v_mesh,4)
    if(initial_area>0):
        scaling=final_area/initial_area
    else:
        scaling=1
    print("final area",final_area)
    print("scaling=",scaling)
    #
    v1_i=weighted_tail(x1,v1,reverse=False)*scale
    v_mesh_i=v1_i.flatten()
    initial_area=area(v_n,v_mesh_i,4)
    z1_i=zpred.reshape(v1.shape)
    z1_i=weighted_tail(x1,z1_i,reverse=False)
    zpred_i=z1_i.flatten()
    final_area=area(v_n,zpred_i-v_mesh_i,4)
    if(initial_area>0):
        scaling_i=final_area/initial_area
    else:
        scaling_i=1
    print("final true area",final_area)
    print("true scaling=",scaling_i)
    for j in range(1,Nf+1):
        best_model.params[pref[j]+'a'].value/=scale
        best_model.best_values[pref[j]+'a']/=scale
    #print(best_params)
    #print(initial_parameters)
    #v1=weighted_tail(x1,v1,reverse=True)
    return best_model,scaling

def read_params(Nf_in,J_in,ch_in,v_n,inn,prange):
    out_file1,out_file2=output_file_name(Nf_in,J_in,ch_in,v_n,inn,prange)
    old_params={}
    with open(out_file1) as f1:
        next(f1)
        for i in range(1,Nf_in+1):
            line=next(f1)
            s=line.split()
            if (len(s)!=5):
                print("error reading previous fit params")
                break
            param_set=['a','xi','xj','wi','wj']
            for j in range(5):
                old_params["f{0:d}_{1}".format(i,param_set[j])]=float(s[j])
    return old_params
def cost_function(params,x,data=None,eps=None,penalty=1):
    parvals=params.valuesdict()
    par_names=['a','xi','xj','wi','wj','ratio_i']
    ridge=0
    for i in range(1,Nf+1):
        key_pars=[parvals["f{}_".format(i)+par] for par in par_names]
        if(i==1):
            model=poly_gaus(x,*key_pars)
        else:
            model=model+poly_gaus(x,*key_pars)
        ridge=ridge+np.abs(parvals["f{}_a".format(i)]**2)
    if data is None:
        return model
    if eps==None:
        return ((model-data)**2).sum()+penalty*ridge
    return ((model-data)**2).sum()/eps

def run_fit_2(v_n,x1,xy,v1,print_output=True):
    """run_fit_2(v_n, x1, xy, v1)
    Core fitting routine per mesh. Modify for different fitting logic.
    """
    #v1=weighted_tail(x1,v1)
    po=print_output
    v_mesh=v1.flatten()
    initial_area=area(v_n,v_mesh,4)
    scale=1
    if(initial_area>0):
        scale=10/initial_area
    v_mesh=v_mesh*scale
    initial_area=area(v_n,v_mesh,4)
    if (po):
        print("initial area: ",initial_area,'scaled by:',scale)
    smallest_area=initial_area
    mse_min=100
    pref=["f{0:d}_".format(i) for i in range(0,Nf+1)]
    params=Parameters()
    for i in range(1,Nf+1):
        if(i==1):
            Mod=Model(poly_gaus,prefix=pref[i])
        else:
            Mod=Mod+Model(poly_gaus,prefix=pref[i])
        params.add(name=pref[i]+'a',min=-10,max=10)
        params.add(name=pref[i]+'xi',min=0.1,max=4)
        params.add(name=pref[i]+'xj',min=0.1,max=4)
        params.add(name=pref[i]+'ratio_j',min=0.45,max=20,value=1)
        params.add(name=pref[i]+'ratio_i',min=0.45,max=20,value=1)
        params.add(name=pref[i]+'wi',expr="f{:d}_xi*f{:d}_ratio_i".format(i,i),min=0.1,max=4.5)
        params.add(name=pref[i]+'wj',expr="f{:d}_xj*f{:d}_ratio_j".format(i,i),min=0.1,max=4.5)
    #Mod=Mod*Model(regulator)
    #params.add(name='cutoff',min=5,max=8)
    if (po):
        print('********')
    time_limit=1000
    penalty=0
    for i in range(tries):
        if (po):
            print("**")
        initial_x,initial_y,initial_w=initial_params()
        for j in range(1,Nf+1):
            params[pref[j]+'a'].value=5 if j<3 else 1 if j<6 else 0.1
            if(initial_area==0):
                params[pref[j]+'a'].value=0

            params[pref[j]+'wi'].value=initial_w[j-1]
            params[pref[j]+'wj'].value=initial_w[j-1]
            params[pref[j]+'xi'].value=initial_x[j-1]
            params[pref[j]+'xj'].value=initial_y[j-1]
        #params['cutoff'].value=6.7
        if (po):
            print('starting the fitting')
        time1=time.time()
        #result=Mod.fit(v_mesh,params,xy=xy,method='trust-constr')
        try:
            signal.alarm(time_limit)
            signal.alarm(0)
            #result=Mod.fit(v_mesh,params,xy=xy,method='least_squares')
            result=minimize(cost_function,params,method='trust-constr',args=(xy,),kws={'data':v_mesh,'penalty':penalty})
            #result=minimize(cost_function,params,method='least_squares',args=(xy,),kws={'data':v_mesh})
            signal.alarm(0)
        except TimeoutException:
            continue
        time2=time.time()
        if (po):
            print("$Time elapsed:",time2-time1)
        if((time2-time1)*2.0 < time_limit):
            time_limit=int((time2-time1)*2.0)
        cost=cost_function(result.params,xy,v_mesh,v_n**2)
        if (po):
            print("cost=",cost)
        curr_params=result.params
        zpred=cost_function(curr_params,x=xy)
        final_area=area(v_n,zpred-v_mesh,4)
        if (po):
            print("current area",final_area)
        if(i==0):
            penalty=cost*v_n*v_n/(25 if scale <40 else 25)
            if (po):
                print("ridge penalty is:",penalty)
            continue

        if(final_area<=smallest_area):
            smallest_area=final_area
            best_model=result

    if (po):
        print('********')
    best_params=best_model.params
    #print("cutoff:",best_model.best_values['cutoff'])
    zpred=cost_function(best_params,x=xy)
    #z1=zpred.reshape(v1.shape)
    final_area=area(v_n,zpred-v_mesh,4)
    if(initial_area>0):
        scaling=final_area/initial_area
    else:
        scaling=1
    if (po):
        print("final area",final_area)
    if (po):
        print("scaling=",scaling)
    #
    v1_i=weighted_tail(x1,v1,reverse=False)*scale
    v_mesh_i=v1_i.flatten()
    initial_area=area(v_n,v_mesh_i,4)
    z1_i=zpred.reshape(v1.shape)
    z1_i=weighted_tail(x1,z1_i,reverse=False)
    zpred_i=z1_i.flatten()
    final_area=area(v_n,zpred_i-v_mesh_i,4)
    if(initial_area>0):
        scaling_i=final_area/initial_area
    else:
        scaling_i=1
    if (po):
        print("final true area",final_area)
        print("true scaling=",scaling_i)
    cost=cost_function(best_model.params,xy,v_mesh,v_n**2)
    if (po):
        print("cost=",cost)

    for j in range(1,Nf+1):
        best_model.params[pref[j]+'a'].value/=scale
        #best_model.best_values[pref[j]+'a']/=scale
    #print(best_params)
    #print(initial_parameters)
    #v1=weighted_tail(x1,v1,reverse=True)
    #best_model.params.pretty_print()
    #print(best_model.init_values)
    return best_model,scaling



def run_fit_successive(v_n,x1,xy,v1,Nf_in,old_params):
    #v1=weighted_tail(x1,v1)
    v_mesh=v1.flatten()
    initial_area=area(v_n,v_mesh,4)
    scale=1
    if(initial_area>0):
        scale=10/initial_area
    v_mesh=v_mesh*scale
    initial_area=area(v_n,v_mesh,4)
    print("initial area: ",initial_area,'scaled by:',scale)
    smallest_area=initial_area
    mse_min=100
    pref=["f{0:d}_".format(i) for i in range(0,Nf+1)]
    params=Parameters()
    for i in range(1,Nf+1):
        if(i==1):
            Mod=Model(poly_gaus,prefix=pref[i])
        else:
            Mod=Mod+Model(poly_gaus,prefix=pref[i])
        params.add(name=pref[i]+'a',min=-20 if i<5 else -20,max=20 if i<5 else 20)
        params.add(name=pref[i]+'xi',min=0.1,max=4)
        params.add(name=pref[i]+'xj',min=0.1,max=4)
        params.add(name=pref[i]+'ratio_j',min=0.45,max=20)
        params.add(name=pref[i]+'ratio_i',min=0.45,max=20)
        params.add(name=pref[i]+'wi',expr="f{:d}_xi*f{:d}_ratio_i".format(i,i),min=0.1,max=4.5)
        params.add(name=pref[i]+'wj',expr="f{:d}_xj*f{:d}_ratio_j".format(i,i),min=0.1,max=4.5)
    print('********')
    time_limit=1000
    for i in range(tries):
        initial_x,initial_y,initial_w=initial_params()
        for j in range(1,Nf_in+1):
            if(initial_area==0):
                params[pref[j]+'a'].value=0
            param_names=['a','wi','wj','xi','xj']
            for param_name in param_names:
                params[pref[j]+param_name].set(value=old_params[pref[j]+param_name],vary=False)
            params[pref[j]+"ratio_i"].set(value=old_params[pref[j]+"wi"]/old_params[pref[j]+"xi"],vary=False)
            params[pref[j]+"ratio_j"].set(value=old_params[pref[j]+"wj"]/old_params[pref[j]+"xj"],vary=False)
            params[pref[j]+'a'].set(value=old_params[pref[j]+'a']*scale,vary=False)
            params[pref[j]+'wi'].set(expr="f{:d}_xi*f{:d}_ratio_i".format(j,j),vary=False)
            params[pref[j]+'wj'].set(expr="f{:d}_xj*f{:d}_ratio_j".format(j,j),vary=False)
        #params['f1_a'].set(value=0.1)
        for j in range(Nf_in+1,Nf+1):
            params[pref[j]+'a'].value=5 if j<3 else 1 if j<6 else 0.1
            if(initial_area==0):
                params[pref[j]+'a'].value=0

            params[pref[j]+'wi'].value=initial_w[j-1]
            params[pref[j]+'wj'].value=initial_w[j-1]
            params[pref[j]+'xi'].value=initial_x[j-1]
            params[pref[j]+'xj'].value=initial_y[j-1]
        #if(i==0):
        #    params.pretty_print()

        #params['cutoff'].value=6.7
        print('starting the fitting')
        time1=time.time()
        #result=Mod.fit(v_mesh,params,xy=xy,method='trust-constr')
        try:
            signal.alarm(time_limit)
            signal.alarm(0)
            #result=Mod.fit(v_mesh,params,xy=xy,method='least_squares')
            result=minimize(cost_function,params,method='trust-constr',args=(xy,),kws={'data':v_mesh})
            signal.alarm(0)
        except TimeoutException:
            continue
        time2=time.time()
        print("$Time elapsed:",time2-time1)
        if((time2-time1)*2.0 < time_limit):
            time_limit=int((time2-time1)*2.0)
        curr_params=result.params
        zpred=cost_function(curr_params,x=xy)
        final_area=area(v_n,zpred-v_mesh,4)
        if(final_area<=smallest_area):
            smallest_area=final_area
            best_model=result
        print("current area",final_area)
        #
        #for i in range(1,Nf+1):
        #    xi=best_model.params["f{:d}_xi".format(i)].value
        #    wi=best_model.params["f{:d}_wi".format(i)].value
        #    ratio=best_model.params["f{:d}_ratio_i".format(i)].value
        #    if(xi/wi>2.5):
        #        print("xi={0} wi={1} ratio:{2},{3}".format(xi,wi,xi/wi,ratio))
        #    xi=best_model.params["f{:d}_xj".format(i)]
        #    wi=best_model.params["f{:d}_wj".format(i)]
        #    ratio=best_model.params["f{:d}_ratio_j".format(i)].value
        #    if(xi/wi>2.5):
        #        print("xj={0} wj={1} ratio:{2},{3}".format(xi,wi,xi/wi,ratio))
    print('********')
    best_params=best_model.params
    #print("cutoff:",best_model.best_values['cutoff'])
    zpred=cost_function(best_params,x=xy)
    #z1=zpred.reshape(v1.shape)
    final_area=area(v_n,zpred-v_mesh,4)
    if(initial_area>0):
        scaling=final_area/initial_area
    else:
        scaling=1
    print("final area",final_area)
    print("scaling=",scaling)
    #
    v1_i=weighted_tail(x1,v1,reverse=False)*scale
    v_mesh_i=v1_i.flatten()
    initial_area=area(v_n,v_mesh_i,4)
    z1_i=zpred.reshape(v1.shape)
    z1_i=weighted_tail(x1,z1_i,reverse=False)
    zpred_i=z1_i.flatten()
    final_area=area(v_n,zpred_i-v_mesh_i,4)
    if(initial_area>0):
        scaling_i=final_area/initial_area
    else:
        scaling_i=1
    print("final true area",final_area)
    print("true scaling=",scaling_i)
    cost=cost_function(best_model.params,xy,v_mesh,v_n**2)
    print("cost=",cost)

    for j in range(1,Nf+1):
        best_model.params[pref[j]+'a'].value/=scale
        #best_model.best_values[pref[j]+'a']/=scale
    #print(best_params)
    #print(initial_parameters)
    #v1=weighted_tail(x1,v1,reverse=True)
    #best_model.params.pretty_print()
    #print(best_model.init_values)
    return best_model,scaling

def input_file_name(J,ij,v_n,inn,prange):
    os.makedirs("fit_input",exist_ok=True)
    file_name1="fit_input/fit-mesh_v_n{0:d}.dat".format(v_n)
    file_name1="fit_input/fit-mesh_v_n{0:d}_prange{1:.2f}.dat".format(v_n,prange)
    #file_name2="fit_input/interaction_mesh_"+inn+"_n3lo_v_n{0:d}_range{1:.2f}_J{2:d}_ch{3:d}".format(v_n,prange,J,ij)+".dat"
    file_name2="fit_input/interaction_mesh_"+inn+'_'+ interaction+"_v_n{0:d}_range{1:.2f}_J{2:d}_ch{3:d}".format(v_n,prange,J,ij)+".dat"
    return file_name1,file_name2
def output_file_name(Nf,J,ij,v_n,inn,prange):
    os.makedirs("lmfit_parameteres",exist_ok=True)
    os.makedirs("lmfit_graphs",exist_ok=True)
    #root="v_lmfit_params_v5_"+inn+"_n3lo_v_n{0:d}_range{4:.2f}_J{1:d}_ch{2:d}_Nf{3:d}".format(v_n,J,ij,Nf,prange)
    root="v_lmfit_params_v5_"+inn+'_'+interaction+"_v_n{0:d}_range{4:.2f}_J{1:d}_ch{2:d}_Nf{3:d}".format(v_n,J,ij,Nf,prange)
    file_name1="lmfit_parameters/"+root+".fit"
    file_name2="lmfit_graphs/"+root
    return file_name1,file_name2

def fit(J,v_n,inn):
    """fit(J, v_n, inn)
    Fit channels with a single total angular momentum J.

    Parameters:
    - J    : Total angular momentum (int)
    - v_n  : Number of mesh points on each axis (int)
    - inn  : Isospin projection ('pp', 'nn', 'np')
    """
    for ij in range(1,6):
        print("******************")
        print("Channel:{0:d}".format(ij))
        in_file1,in_file2=input_file_name(J,ij,v_n,inn,prange)
        out_file1,out_file2=output_file_name(Nf,J,ij,v_n,inn,prange)
        x1,xy,v1=fromfile(v_n,in_file1,in_file2)
        #old_params=read_params(4,J,ij,v_n,inn,prange)
        #print(old_params)
        #best_fit,scaling=run_fit_successive(v_n,x1,xy,v1,4,old_params)
        best_fit,scaling=run_fit_2(v_n,x1,xy,v1)
        #write_parameters(best_fit,scaling,out_file1)
        plot(x1,xy,v1,best_fit,out_file2)
        #sys.exit()
        #best_fit,scaling=run_fit(v_n,x1,xy,v1)
def fitallj(Jmin,Jmax,v_n,inn,Nf_in=Nf):
    """fitallj(Jmin, Jmax, v_n, inn)
    Fit interactions for all J values in a range and a given isospin.

    Parameters:
    - Jmin : Minimum total angular momentum
    - Jmax : Maximum total angular momentum
    - v_n  : Mesh resolution
    - inn  : Isospin projection ('pp', 'nn', 'np')
    """
    global Nf
    timei=time.time()
    Nf=Nf_in
    print(Nf)
    for J in range(Jmin,Jmax+1):
        print("******************************************")
        print("STARTING J={0:d}".format(J))
        for ij in range(1,7):
            print("\t******************")
            print("\tChannel:{0:d}".format(ij))
            in_file1,in_file2=input_file_name(J,ij,v_n,inn,prange)
            out_file1,out_file2=output_file_name(Nf,J,ij,v_n,inn,prange)
            x1,xy,v1=fromfile(v_n,in_file1,in_file2)
            best_fit,scaling=run_fit_2(v_n,x1,xy,v1)
            write_parameters(best_fit,scaling,out_file1)
            plot(x1,xy,v1,best_fit,out_file2)
            timef=time.time()
            iterations_finished=(J-Jmin)*6+ij
            time_remain=(timef-timei)*((Jmax-Jmin+1)*6-iterations_finished)/iterations_finished
            print("Remaining time ~ "+time_string(time_remain))
            print("\t******************")
    timef=time.time()
    print("Calculation took: "+ time_string(timef-timei))
    print("******************")
def fitall(Jmin,Jmax,v_n,Nf_in=Nf):
    """fitall(Jmin, Jmax, v_n)
    Fit all J values and all isospin projections.

    Parameters:
    - Jmin : Minimum J
    - Jmax : Maximum J
    - v_n  : Mesh resolution
    """
    global Nf
    Nf=Nf_in
    inn=['nn','np','pp']
    timei=time.time()
    for i in range(3):
        print("******************************************")
        print("STARTING "+inn[i])
        for J in range(Jmin,Jmax+1):
            print("\t******************************************")
            print("\tSTARTING J={0:d}".format(J))
            for ij in range(1,7):
                print("\t\t******************")
                print("\t\t Channel:{0:d}".format(ij))
                in_file1,in_file2=input_file_name(J,ij,v_n,inn[i],prange)
                out_file1,out_file2=output_file_name(Nf,J,ij,v_n,inn[i],prange)
                x1,xy,v1=fromfile(v_n,in_file1,in_file2)
                #
                best_fit,scaling=run_fit_2(v_n,x1,xy,v1)
                write_parameters(best_fit,scaling,out_file1)
                plot(x1,xy,v1,best_fit,out_file2)


                timef=time.time()
                iterations_finished= i*(Jmax-Jmin+1)*6+(J-Jmin)*6+ij
                time_remain=(timef-timei)*(3*(Jmax-Jmin+1)*6-iterations_finished)/iterations_finished
                print("Remaining time ~ "+time_string(time_remain))
                print("\t\t******************")

    timef=time.time()
    print("Calculation took: "+ time_string(timef-timei))
    print("******************")


import concurrent.futures

def parallel_func(J,ij,v_n,inn):
    in_file1,in_file2=input_file_name(J,ij,v_n,inn,prange)
    out_file1,out_file2=output_file_name(Nf,J,ij,v_n,inn,prange)
    x1,xy,v1=fromfile(v_n,in_file1,in_file2)
                #
    best_fit,scaling=run_fit_2(v_n,x1,xy,v1,print_output=False)
    write_parameters(best_fit,scaling,out_file1)
    plot(x1,xy,v1,best_fit,out_file2)
    return J,ij,inn,scaling

# Main function to setup and run the parallel tasks
def fit_parallel(Jmin,Jmax,v_n,inn,Nf_in):
    """fit_parallel(Jmin, Jmax, v_n, inn,Nf_in)
    Same as fitallj but runs in parallel for better performance.

    Parameters:
    - Jmin : Minimum J
    - Jmax : Maximum J
    - v_n  : Mesh resolution
    - [inn]  : List of isospin projection ('pp', 'nn', 'np')
    - Nf_in: the number of functions to be used in the fit
    """
    global Nf
    Nf=Nf_in
    timei=time.time()
    with concurrent.futures.ProcessPoolExecutor() as executor:
        # Submit all tasks to the executor

        futures = [executor.submit(parallel_func, J, ij,v_n,inn1) for J in range(Jmin,Jmax+1) for ij in range(1,7) for inn1 in inn]

        # Process the results as they complete
        for future in concurrent.futures.as_completed(futures):
            try:

                J_f,ij_f,inn_f,scaling=future.result()
                print("finished fitting J={0},ij={1},inn={2}: scaling={3}".format(J_f,ij_f,inn_f,scaling))
            except Exception as exc:
                print('error')
    timef=time.time()
    print("Calculation took: "+ time_string(timef-timei))
    print("******************")



input_output_notes = {
    'input_file_name': "Specifies input file path and name.",
    'output_file_name': "Specifies output file path and name."
}

parameters = {
    "prange":       "Momentum range of the mesh",
    "Nf":           "Number of functions used in the fit (2–16)",
    "interaction":  "Label used for the interaction and output files",
    "tries":        "Number of randomized fitting attempts"
}
def help():
    print("""
==================== HELP: Fitting Script =====================

This script fits NN interactions on a momentum-space mesh.
Reads from `fit_input/`, writes to `lmfit_parameters/`.

---------------------------------------------------------------
Required Parameters:
---------------------------------------------------------------
""")
    for k, v in parameters.items():
        print(f"  {k:<14}: {v}")

    print("""
---------------------------------------------------------------
User-Callable Functions:
---------------------------------------------------------------
""")
    user_funcs = [fit, fitallj, fitall, fit_parallel]
    for f in user_funcs:
        doc=inspect.getdoc(f)
        if doc:
            lines=doc.split('\n',1)
            signature=lines[0]
            description=textwrap.indent(textwrap.fill(lines[1],width=75),prefix='\t\t') if len(lines)>1 else ''
            print(f"  {f.__name__:<20}: {signature}")
            if description:
                print(description)
        print()

    print("---------------------------------------------------------------")
    print("Internal Functions (not for user use unless customizing logic):")
    print("---------------------------------------------------------------")
    print(inspect.getdoc(run_fit_2))
    for name, desc in input_output_notes.items():
        print(f"  {name}: {desc}")

    print("""
---------------------------------------------------------------
Usage Example:
---------------------------------------------------------------
  python -c 'import fitting_lmfit_cost as flc; flc.fit_parallel(0, 4, 40, "np")'
===============================================================
""")


if __name__=='__main__':
    help()
