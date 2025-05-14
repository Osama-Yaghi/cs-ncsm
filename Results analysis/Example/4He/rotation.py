import numpy as np
import pandas as pd
import seaborn as sns
from math import ceil, sqrt
import math
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm,Normalize
from scipy.ndimage.filters import gaussian_filter
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes, mark_inset
import os
#plt.rcParams["figure.figsize"] = (2,2)
v=np.zeros((15,50000))
E_r=np.zeros((15,50000))
E_i=np.zeros((15,50000))
E_r_sorted=np.zeros((15,50000))
E_i_sorted=np.zeros((15,50000))
x=[[]]
thetas_in=np.zeros(100)
thresholds=[-8.482,-7.718,-4.448,-2.224,0]
def read(ind,file_name,pos,skip):
    l=0
    with open(file_name) as f:
        for j in range(skip):
            next(f)
        for line in f:
            if('***********' in line):
                continue
            if('theta=' in line):
                ind+=1
                l=0
                thetas_in[ind]=line.split()[1]
                print('theta',thetas_in[ind])
                continue
            s=line.split()
            E_r[ind][l]=float(s[pos])
            E_i[ind][l]=float(s[pos+1])
            l+=1
    #print(E_r[ind][:3])
    #print(E_i[ind][:3])
def ang_v(o_x,o_y,x1,y1,x2,y2,dt):
    x_avg=(x1+x2)/2.0
    y_avg=(y1+y2)/2.0
    Lx=x_avg-o_x
    Ly=y_avg-o_y
    L=np.sqrt(Ly**2+Lx**2)
    dx=x2-x1
    dy=y2-y1
    d=np.sqrt(dx**2+dy**2)
    d_L=(dx*Lx+dy*Ly)/L

    dw=np.sqrt(d**2-d_L**2)/L
    return dw/dt
def dir_v(x1,y1,x2,y2,dt):
    dx=x2-x1
    dy=y2-y1
    d=np.sqrt(dx**2+dy**2)/dt
    return d
def plot_res(n_res,thetas,name):
    plt.rcParams.update({'font.size': 12})
    fig,ax=plt.subplots(3,3)
    n_p=16
    color=['r','g','b','y','c','tab:brown','m']
    n_thetas=len(thetas)
    v=np.zeros((n_res,n_thetas))
    w=np.zeros((n_res,n_thetas,len(thresholds)))
    E_x_temp=np.zeros(100)
    threshold=np.zeros(100)
    for i in range(n_thetas):
        E_x_temp[:n_p]=E_r[i][:n_p]*np.cos(2*thetas[i])-E_i[i][:n_p]*np.sin(2*thetas[i])
        indices=np.argsort(E_x_temp[:n_p])[:n_res]
        E_r_sorted[i][:n_res]=E_r[i][indices]
        E_i_sorted[i][:n_res]=E_i[i][indices]
    for i in range(n_thetas-1):
        #v[:n_res,i]=np.sqrt((E_r_sorted[i+1][:n_res]-E_r_sorted[i][:n_res])**2+(E_i_sorted[i+1][:n_res]-E_i_sorted[i][:n_res])**2)/(thetas[i+1]-thetas[i])
        v[:n_res,i]=dir_v(E_r_sorted[i][:n_res],E_i_sorted[i][:n_res],E_r_sorted[i+1][:n_res],E_i_sorted[i+1][:n_res],thetas[i+1]-thetas[i])
        for th in range(len(thresholds)):
            w[:n_res,i,th]=ang_v(thresholds[th],0,E_r_sorted[i][:n_res],E_i_sorted[i][:n_res],E_r_sorted[i+1][:n_res],E_i_sorted[i+1][:n_res],thetas[i+1]-thetas[i])

    #for i in range(n_thetas-1):
        #for th in range(len(thresholds)):
            #w[:n_res,i,th]=v[:n_res,i]/np.sqrt((0.5*E_r_sorted[i+1][:n_res]+0.5*E_r_sorted[i][:n_res]-thresholds[th])**2+(0.5*E_i_sorted[i+1][:n_res]+0.5*E_i_sorted[i][:n_res])**2)
    for ind in range(n_res):
        axj=ind%3
        axi=int(ind/3)
        ax[axi,axj].plot(thetas[:n_thetas-1],v[ind,:n_thetas-1],color=color[2],label="state_num={0:d}".format(ind),marker='o')
        ax[axi,axj].axhline(y=0,color='k')
        ax[axi,axj].axvline(x=0,color='k')
        #ax[axi,axj].set_xlim(-30,33)
        #ax[axi,axj].set_ylim(-30,10)
        #ax[axi,axj].set_title("θ={:0.3f}".format(thetas[ind]),y=0,pad=-25,fontdict={'fontsize':14})
        ax[axi,axj].legend( fontsize=12)
    ax[2,1].set_xlabel("Complex-Scaling angle [Rad]",fontsize=22)
    ax[1,0].set_ylabel("speed [Mev]",fontsize=22)
    #rect = plt.Rectangle((-2.9, -1), 2.5, 2,
    #                     edgecolor='black',linestyle='dashed', linewidth=0.2, fill=False)
    #ax.add_patch(rect)
    #ax.set_xlabel("Real(E) [Mev]",fontsize=22)
    #ax.set_ylabel("Imaginary(E) [Mev]",fontsize=22)
    #ax.tick_params(axis='both', which='major', labelsize=20)
    plt.show()
    fig.savefig(name+".pdf")
    fig.savefig(name+".svg")
    fig,ax=plt.subplots(3,3)
    for ind in range(n_res):
        axj=ind%3
        axi=int(ind/3)
        for th in range(len(thresholds)):
            ax[axi,axj].plot(thetas[:n_thetas-1],w[ind,:n_thetas-1,th],color=color[2],label="state_num={0:d}\nthreshold num {1:d}".format(ind,th+1),marker='o')
        ax[axi,axj].axhline(y=0,color='k')
        ax[axi,axj].axvline(x=0,color='k')
        ax[axi,axj].axhline(y=2,color='r')
        #ax[axi,axj].set_xlim(-30,33)
        ax[axi,axj].set_ylim(-0.1,4)
        ax[axi,axj].set_title("state number:{:d}".format(ind+1),y=0,pad=-25,fontdict={'fontsize':12})
            #ax[axi,axj].legend(loc='lower right', fontsize=5)
    ax[2,1].set_xlabel("Complex-Scaling angle [Rad]",fontsize=22)
    ax[1,0].set_ylabel("rotation speed in the complex plane",fontsize=22)
    plt.show()
    fig.savefig(name+"_angular"+".pdf")
    fig.savefig(name+"_angular"+".svg")
def plot2(ind_range,num_p,thetas,name):
    fig,ax=plt.subplots(3,3)
    plt.rcParams.update({'font.size': 22})
    color=['r','g','b','y','c','tab:brown','m']
    l=0
    print(thetas,thetas_in)
    for ind in ind_range:
        axj=ind%3
        axi=int(ind/3)
        ax[axi,axj].set_xlim(-33,33)
        ax[axi,axj].set_ylim(-10,4)
        while(not (thetas_in[l] in thetas)):
            print(thetas_in[l])
            l+=1
        ax[axi,axj].scatter(E_r[l][:num_p],E_i[l][:num_p],color=color[2],label="θ={:0.3f}".format(thetas[ind]),linewidths=3,clip_on=True,marker='+')
        ##########
        for j in range(5):
            xx=np.linspace(0,40,40)
            yy=-xx*math.tan(2*float(thetas[ind]))
            xx[:]=xx[:]+thresholds[j]
            ax[axi,axj].plot(xx,yy,color='k',linestyle='--')
        ##########
        ax[axi,axj].axhline(y=0,color='k')
        ax[axi,axj].axvline(x=0,color='k')
        #ax[axi,axj].set_title("θ={:0.3f}".format(thetas[ind]),y=0,pad=-25,fontdict={'fontsize':14})
        ax[axi,axj].legend(loc='lower left', fontsize=14)
        l+=1
    ax[2,1].set_xlabel("Real(E) [Mev]",fontsize=22)
    ax[1,0].set_ylabel("Imaginary(E) [Mev]",fontsize=22)
    #rect = plt.Rectangle((-2.9, -1), 2.5, 2,
    #                     edgecolor='black',linestyle='dashed', linewidth=0.2, fill=False)
    #ax.add_patch(rect)
    #ax.set_xlabel("Real(E) [Mev]",fontsize=22)
    #ax.set_ylabel("Imaginary(E) [Mev]",fontsize=22)
    #ax.tick_params(axis='both', which='major', labelsize=20)
    plt.show()
    fig.savefig(name+".pdf")
    fig.savefig(name+".svg")
def plot(ind_range,num_p,thetas,name):
    fig,ax=plt.subplots(1,1)
    plt.rcParams.update({'font.size': 22})
    color=['r','g','b','y','c','tab:brown','m','tab:orange']
    ax.axhline(y=0,color='k')
    ax.axvline(x=0,color='k')
    axins = zoomed_inset_axes(ax, zoom=8, loc='upper center')  # 'zoom' controls the zoom level
    axins.axhline(y=0,color='k')
    axins.set_xlim(-2.4,-1.4)
    axins.set_ylim(-0.2,0.2)
    for ind in ind_range:
        ax.scatter(E_r[ind][:num_p],E_i[ind][:num_p],color=color[ind],label="θ={:0.2f}".format(thetas[ind]),linewidths=6,clip_on=False)
        axins.scatter(E_r[ind][:1], E_i[ind][:1],color=color[ind])
        xx=np.linspace(0,40,40)
        print(2*float(thetas[ind]))
        yy=-xx*math.tan(2*float(thetas[ind]))
        ax.plot(xx,yy,color=color[ind])
    #rect = plt.Rectangle((-2.9, -1), 2.5, 2,
    #                     edgecolor='black',linestyle='dashed', linewidth=0.2, fill=False)
    mark_inset(ax, axins, loc1=2, loc2=4, fc="none", ec="0.5")
    #ax.add_patch(rect)
    ax.set_xlim(-4,8)
    ax.set_ylim(-14,10)
    ax.legend()
    ax.set_xlabel("Real(E) [Mev]",fontsize=22)
    ax.set_ylabel("Imaginary(E) [Mev]",fontsize=22)
    ax.tick_params(axis='both', which='major', labelsize=20)
    plt.show()
    fig.savefig(name+".pdf")
    fig.savefig(name+".svg")
def rot(dirr,filename,thetas):
    pos=0
    skip=0
    n_th=len(thetas)
    read(-1,dirr+filename,pos,skip)
    plot(np.arange(n_th),8,thetas,dirr+filename)
def CS_spect(dirr,filename,thetas):
    pos=0
    skip=0
    n_th=len(thetas)
    read(-1,dirr+filename,pos,skip)
    print(thetas,thetas_in)
    plot2(np.arange(min(n_th,9)),16,thetas[max(n_th-9,0):n_th],dirr+filename)
    plot_res(8,thetas,dirr+filename+"res_speed")
data=pd.DataFrame(columns=['E_r','E_i','N_max','theta','J','T','parity','srg','lambda'])
def readfiles(name1):
    files=os.listdir(os.getcwd())
    matching_files=[file for file in files if file.startswith(name1)]
    print(matching_files)
    for file in matching_files:
        srg=0
        srglambda=0
        if('srg' in file):
            srg=1
            extracted=file.split('srg')[1].split('_')[0]
            #print(file,extracted)
            srglambda=float(extracted)
        if('theta' in file):
            extracted=file.split('theta')[1]
            theta=float(extracted)
            print(file,theta)
        if('_N' in file.lstrip(name1)):
            extracted=file.lstrip(name1).split('_N')[1].split('_')[0]
            N_max=int(extracted)
        if('_T' in file.lstrip(name1)):
            extracted=file.split('_T')[1].split('_')[0]
            T=int(extracted)
        if('_J' in file.lstrip(name1)):
            extracted=file.split('_J')[1].split('_')[0]
            J=int(extracted)
        if(N_max%2==0):
            parity=0
        else:
            parity=1
        with open(file) as f1:
            for line in f1:
                if(not( "Spectrum" in line)):
                    continue
                N=int(line.split('N_max=')[1])
                E_r_temp=[]
                E_i_temp=[]
                line=next(f1)
                while(not("End of spectrum" in line)):
                    s=line.split()
                    E_r_temp.append(float(s[0]))
                    E_i_temp.append(float(s[1]))
                    line=next(f1)
                length=len(E_r_temp)
                theta_temp=[theta for i in range(length)]
                N_temp=[N for i in range(length)]
                srg_temp=[srg for i in range(length)]
                srglam_temp=[srglambda for i in range(length)]
                J_temp=[J for i in range(length)]
                T_temp=[T for i in range(length)]
                parity_temp=[parity for i in range(length)]
                new_data={"E_r":E_r_temp,"E_i":E_i_temp,"N_max":N_temp,"theta":theta_temp,'srg':srg_temp,'J':J_temp,'T':T_temp,'parity':parity_temp,'lambda':srglam_temp}
                new_pd=pd.DataFrame(new_data)
                data=pd.concat([data,new_pd])
    print(data.head(n=30))
