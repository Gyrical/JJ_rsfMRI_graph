#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 22 08:17:33 2022

@author: tianlu
"""

import numpy as np
import numpy.matlib
import matplotlib.pyplot as plt
from matplotlib import cm
import seaborn as sns
import pandas as pd
import scipy.interpolate
from scipy.stats import pearsonr
import scipy.io
import os


# Load data and prepare settings
Dj = scipy.io.loadmat('fc_graph_data.mat')

# --- 220629 Jessica ---
D = {'CCall_sp' : Dj['pat_all_cc'],
     'CCall_ll' : Dj['pat_dom_cc'],
     'CCall_rl' : Dj['pat_nondom_cc'],
     'CCall_hc' : Dj['con_cc'],
     'pl_sp'    : Dj['pat_all_cp'],
     'pl_ll'    : Dj['pat_dom_cp'],
     'pl_rl'    : Dj['pat_nondom_cp'],
     'pl_hc'    : Dj['con_cp'],
     'CCm1_sp'  : Dj['pat_all_preCGcc'],
     'CCm1_ll'  : Dj['pat_dom_preCGcc'],
     'CCm1_rl'  : Dj['pat_nondom_preCGcc'],
     'CCm1_hc'  : Dj['con_preCGcc'],
     'bc_sp'    : Dj['pat_all_preCGcb'],
     'bc_ll'    : Dj['pat_dom_preCGcb'],
     'bc_rl'    : Dj['pat_nondom_preCGcb'],
     'bc_hc'    : Dj['con_preCGcb'],
     'wb_sp'    : Dj['m0_all'],
     'wb_hc'    : Dj['m0_con'],
     'wb_ll'    : Dj['m0_dom'],
     'wb_rl'    : Dj['m0_nondom'],
     'ho_sp'    : Dj['m4_all'],
     'ho_hc'    : Dj['m4_con'],
     'ho_ll'    : Dj['m4_dom'],
     'ho_rl'    : Dj['m4_nondom']
    }

measl  = 'CCall pl CCm1 bc wb ho'.split() # list(set([x.split('_')[0] for x in keyl]))
groupl = 'sp hc ll rl'.split() #list(set([x.split('_')[1] for x in keyl]))
# ---

# Dl = scipy.io.loadmat(dir_main+'lesionsize.mat')
# lesionsize = Dl['lesionsize']/1000
# lesionsize = lesionsize.astype(np.float64).flatten()
lesionsize = Dj['lesion'].flatten()/1000

sparsl = np.arange(0.1,0.95,0.05) 
keyl   = [x for x in D.keys() if '__' not in x]
labsgm = 'CC$_{avg}$ PL CC$_{preCG}$ BC$_{preCG}$ FC$_{wb}$ FC$_{h}$'.split()
labsgr = ['Stroke patients','Healthy controls','LL','RL']

# https://matplotlib.org/stable/gallery/color/color_cycle_default.html
pltcolors = [plt.rcParams['axes.prop_cycle'].by_key()['color'][x] for x in [4,7,0,1]]#[:4]
pltmarkers = ['^','o','<','>'] # SP HC LL RL

# sublabs  = 'A B C D E F'.split()
sublabs  = ['A','','B','','C','']
sublfont = {'fontsize':14,'fontweight':'bold','fontname':'Arial'}
legendst = {'bbox_to_anchor':(1.01, 1),'frameon':False}
legendsts = {'bbox_to_anchor':(1.01, 0.95),'frameon':False}
scatterst = {'s':30,'zorder':3,'label':'_nolegend_'}
signst1 = {'ha':'center','va':'bottom','fontsize':14}
signst2 = {'ha':'center','va':'baseline','fontsize':14}

plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['font.family'] = 'sans-serif'

sstats = True
pchoic = 1 # 1: uncorrected, 3: bonferroni, 5: fdr

#%% Fig 2 (2-3) Combined SP vs. HC
    
statsm = [['*' if x<.05 else '' for x in Dj['stats_all_'+xx][5,:]] for xx in 'cc cp pcc pcb'.split()] # else '^' if x<.1

fig,axes = plt.subplots(3,2,figsize=(10,9))
for idxm,labm in enumerate(measl):
    ax = axes.flatten()[idxm]
    
    if idxm<4:
        
        for idxg,labg in enumerate(groupl[:2]):
            
            # Prepare data
            data = D[labm+'_'+labg]
            data[data==np.inf] = np.nan
            x = sparsl
            y = np.nanmean(data,0)
            yerr = np.nanstd(data,0)/np.sqrt(np.sum(~np.isnan(data),0))
            
            # Plot data
            if 0:
                x+=idxg*0.01
                ax.errorbar(x,y,yerr,color=pltcolors[idxg],marker=pltmarkers[idxg],capsize=5)
            if 1:
                ysl = scipy.interpolate.interp1d(x,y-yerr)
                ysh = scipy.interpolate.interp1d(x,y+yerr)
                xs = np.linspace(x[0],x[-1],len(sparsl)*5)
                
                ax.plot(x,y,color=pltcolors[idxg],marker=pltmarkers[idxg])
                ax.fill_between(xs,ysl(xs),ysh(xs),alpha=0.25,color=pltcolors[idxg],label='_nolegend_')
        
        if idxm==1: 
            xls,yls = ax.get_xlim(), ax.get_ylim()
            yls = [yls[0],yls[1]*1.022];ax.set_ylim(top = yls[1])
        
        # Show stats markers    
        if sstats:
            # show stats
            xls,yls = ax.get_xlim(), ax.get_ylim()
            for idxs,xs in enumerate(statsm[idxm]): 
                ax.text(sparsl[idxs],yls[1]-np.diff(yls)*0.08,statsm[idxm][idxs],**signst1)
        
        # Display settings
        ax.set_xticks(sparsl[::2])
        ax.set_xticklabels(['{:.2f}'.format(x) for x in sparsl[::2]])
        ax.set_xlabel('sparsity')
        
# ------------------------------------------------------------------------------
# Show FC

    else:
        x = np.arange(1,3)  
        # Prepare data
        data = [[x for x in np.squeeze(D[labm+'_'+groupl[0]]) if x not in [np.nan,np.inf]],
                [x for x in np.squeeze(D[labm+'_'+groupl[1]]) if x not in [np.nan,np.inf]]]
        
        # Plot data    
        violins = ax.violinplot(data,showmeans=False,showmedians=False,showextrema=False)
        for idxv,v in enumerate(violins['bodies']): v.set_facecolor(pltcolors[idxv])
        
        qmq=np.asarray([np.percentile(x,[50,25,75]) for x in data])
        for idx in range(2):
            ax.scatter(x[idx],qmq[idx,0],marker='o',color='white',**scatterst) # pltmarkers[idx]
            ax.plot([idx+1,idx+1],qmq[idx,(1,2)],color='k',label='_nolegend_') #pltcolors[idx]
        
        # Show stats markers
        if sstats and idxm==5:
            xls,yls = ax.get_xlim(), ax.get_ylim()
            yls = [yls[0],yls[1]*1.022]
            ax.set_ylim(top = yls[1])
            esc = np.diff(yls)*0.07
            el = (np.ones((2,))*np.diff(yls)*0.02,np.zeros((2,)))
            ax.errorbar([1,2],np.ones((2,))*yls[1],el,color='k')
            ax.text(1.5,yls[1],'*',**signst2)
            
        # Display settings
        ax.set_xticks([1,2])
        ax.set_xticklabels(labsgr[:2])

    # General display settings
    xls,yls = ax.get_xlim(), ax.get_ylim()
    ax.text(xls[0]-np.diff(xls)*0.25,yls[1],sublabs[idxm],**sublfont)
    
    if idxm==1 or idxm==5: ax.legend(labsgr[:2],**legendsts)
    # if idxm==5: ax.legend(labsgr[:2],**legendst)
     
    ax.set_ylabel(labsgm[idxm])           
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    
fig.tight_layout()
fn_out = os.path.dirname(fn)+ '/fig2_largegroups.png' 
plt.savefig(fn_out, dpi=300)
print(fn_out+' saved!')
plt.show()


#%% Figure 3 (4-5) combined LL vs. RL vs. HC

slabs = ['$LL ≠ RL$','$LL ≠ HC$','$RL ≠ HC$']
statsma = [[['*' if x<.05 else '' for x in Dj['stats_{}_{}'.format(xx,yy)][1,:]] 
              for xx in 'dom_nondom dom nondom'.split()] for yy in 'cc cp pcc pcb'.split()] # else '^' if x<.1 

fig,axes = plt.subplots(3,2,figsize=(10,9))
for idxm,labm in enumerate(measl):
    ax = axes.flatten()[idxm]
    
    if idxm < 4:
        for idxg,labg in enumerate(groupl[1:]):
            
            # Prepare data
            if labm+'_'+labg not in D.keys(): continue
            data = D[labm+'_'+labg]
            data[data==np.inf] = np.nan
            x = sparsl
            y = np.nanmean(data,0)
            yerr = np.nanstd(data,0)/np.sqrt(np.sum(~np.isnan(data),0))
            
            # Plot data
            if 0:
                x+=idxg*0.01
                ax.errorbar(x,y,yerr,color=pltcolors[idxg+1],marker=pltmarkers[idxg+1],capsize=5)
            if 1:
                ysl = scipy.interpolate.interp1d(x,y-yerr)
                ysh = scipy.interpolate.interp1d(x,y+yerr)
                xs = np.linspace(x[0],x[-1],len(sparsl)*5)
                
                ax.plot(x,y,color=pltcolors[idxg+1],marker=pltmarkers[idxg+1],label=labsgr[idxg+1])
                ax.fill_between(xs,ysl(xs),ysh(xs),alpha=0.25,color=pltcolors[idxg+1],label='_nolegend_')
        
        # show stats
        if sstats:
            xls,yls = ax.get_xlim(), ax.get_ylim()
            # yls = [yls[0],yls[1]*1.072]; ax.set_ylim(top = yls[1])
            
            idxs = [idx for idx,x in enumerate(statsma[idxm]) if len(''.join(x))>0]
            for idxx,xx in enumerate(idxs):
                for idxss,ss in enumerate(statsma[idxm][xx]):
                    ax.text(sparsl[idxss],yls[1]-np.diff(yls)*0.08*idxx,ss,color=pltcolors[xx+1],**signst1)
                ax.text(sparsl[np.max([idx for idx,x in enumerate(statsma[idxm][xx]) if len(x)>0])]+0.03,
                        yls[1]+np.diff(yls)*(0.03-0.08*idxx),slabs[xx],va='bottom',color=pltcolors[xx+1],fontsize=8)
                
        # Display settings
        ax.set_xticks(sparsl[::2])
        ax.set_xticklabels(['{:.2f}'.format(x) for x in sparsl[::2]])
        ax.set_xlabel('sparsity')

# ------------------------------------------------------------------------------
# Show FC
            
    else:
        # Prepare data
        x = np.arange(1,4)
        data = [[x for x in np.squeeze(D[labm+'_'+x]) if x not in [np.nan,np.inf]] for x in groupl[1:]]
        
        # Plot data    
        violins = ax.violinplot(data,showmeans=False,showmedians=False,showextrema=False)
        for idxv,v in enumerate(violins['bodies']): 
            v.set_facecolor(pltcolors[idxv+1])
            v.set_label(labsgr[1+idxv])
            # scatter
            # ax.scatter(idxv+.75+np.random.rand(len(data[idxv]),)*.5,data[idxv],3,color=pltcolors[idxv+1])
        
        qmq=np.asarray([np.percentile(x,[50,25,75]) for x in data])
        ax.scatter(x,qmq[:,0],marker='o',color='white',s=30,zorder=3,label='_nolegend_')
        for idx in x: ax.plot([idx,idx],qmq[idx-1,(1,2)].flatten(),color='k',label='_nolegend_')
            
        xls,yls = ax.get_xlim(), ax.get_ylim()
        if idxm==5:
            yls = [yls[0],yls[1]*1.022]
            ax.set_ylim(top = yls[1])
            esc = np.diff(yls)*0.07
            el = (np.ones((2,))*np.diff(yls)*0.02,np.zeros((2,)))
            ax.errorbar([1,3],np.ones((2,))*yls[1],el,color='k')
            ax.text(2,yls[1],'*',ha='center',va='baseline',fontsize=14)
            
            ax.errorbar([1,1.9],np.ones((2,))*yls[1]-esc,el,color='k')
            ax.text(1.5,yls[1]-esc,'*',ha='center',va='baseline',fontsize=14)
            ax.errorbar([2.1,3],np.ones((2,))*yls[1]-esc,el,color='k')
            ax.text(2.5,yls[1]-esc,'*',ha='center',va='baseline',fontsize=14)
            
        # Display settings FC
        ax.set_xticks([1,2,3])
        ax.set_xticklabels(labsgr[1:])
        
    # General display settings
    ax.set_ylabel(labsgm[idxm])
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    xls,yls = ax.get_xlim(), ax.get_ylim()
    ax.text(xls[0]-np.diff(xls)*0.25,yls[1],sublabs[idxm],**sublfont)
    if idxm==1 or idxm==5: 
        handles, labels = ax.get_legend_handles_labels()
        order = [1,2,0]
        ax.legend([handles[idx] for idx in order],[labels[idx] for idx in order],**legendsts)
    # ax.legend(labsgr[1:],**legendst)
    
fig.tight_layout() 
plt.savefig(os.path.dirname(fn)+ '/fig3_subgroups.png', dpi=300)
print('fig3_subgroups.png saved!')
plt.show()


#%% Figure 4 with correlations
import re
if 1:
    # option 1: show single sparsity
    fig,axes = plt.subplots(3,2,figsize=(8,9))
    
    for idxm,labm in enumerate(measl):
        ax = axes.flatten()[idxm]
        
        # Prepare data
        if idxm>3: data = D[labm+'_sp'].flatten()
        else: data = D[labm+'_sp'][:,np.where(np.isclose(sparsl,0.75))[0][0]].flatten() # 4: 0.75
        
        # Plot
        if 1:    
            ldata = pd.DataFrame({'data':data,'lesionsize':lesionsize})
            g = sns.regplot(x='lesionsize',y='data',data=ldata,ax=ax,marker='^',color=pltcolors[0])
            ax.set_xlim([-10,370])
        else:
            ax.scatter(lesionsize,data,marker=pltmarkers[0],c=pltcolors[0])
            ax.plot([np.min(lesionsize),np.max(lesionsize)],[np.mean(data)])
    
        xls,yls = ax.get_xlim(), ax.get_ylim()
        if idxm<4: ax.text(xls[0]+np.diff(xls)*0.05,yls[1],'$sparsity = 0.75$',ha='left',va='bottom',fontsize=8)
            
        # plot stats text
        if sstats:
            filt = ~np.isnan(np.vstack((lesionsize,data))).any(axis=0)
            r,p = pearsonr(lesionsize[filt],data[filt])
            tstr = re.sub('0(?=[.])','','$r = {:.2f}, p = {:.3f}$'.format(r,p))
            if p < 0.05: tstr += '*'
            ax.text(xls[1]-np.diff(xls)*0.1,yls[1]-np.diff(yls)*0.15,tstr,ha='right',va='bottom',fontsize=8)
            
        # Display options
        ax.set_ylabel(labsgm[idxm])
        ax.set_xlabel('lesion size in cm$^3$')
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.text(xls[0]-np.diff(xls)*0.25,yls[1],sublabs[idxm],**sublfont)


# Option 2: show all sparsities

else:
    cmap = np.flip(cm.get_cmap('Purples')(np.arange(0.2,1,0.1)),axis=0)
    
    # Prepare data
    dfs=[]    
    for idxm,labm in enumerate(measl):
        data  = D[labm+'_sp']
        if idxm < 4:
            group = numpy.matlib.repmat(sparsl,data.shape[0],1).T.flatten()
            x     = numpy.matlib.repmat(lesionsize,1,data.shape[1])
        else:
            group = np.min(sparsl)*np.ones(data.size,)
            x     = lesionsize
        dfs.append(pd.DataFrame({'data':data.T.flatten(),'gm':[labm]*data.size,
                              'lesionsize':np.squeeze(x),'sparsity':['{:.02f}'.format(x) for x in group]}))
    dfall = pd.concat(dfs,ignore_index=True)
        
    # Plot graph metrics
    g = sns.lmplot(x='lesionsize',y='data',hue='sparsity',col='gm',data=dfall,height=3,aspect=1.2,
                   col_wrap=2,sharex=False,sharey=False,palette=cmap,markers='^',
                   n_boot=10)
    
    # Display settings
    g.legend.set_bbox_to_anchor([1.08,0.85])
    axes = g.axes.flatten()
    for idxm,ax in enumerate(axes):
        ax.set_title('')
        ax.set_ylabel(labsgm[idxm])
        ax.set_xlim([-5,355])
        ax.set_xlabel('lesion size in cm$^3$')
        xls,yls = ax.get_xlim(), ax.get_ylim()
        ax.text(xls[0]-np.diff(xls)*0.25,yls[1],sublabs[idxm],**sublfont)
    

plt.tight_layout()
plt.savefig(os.path.dirname(fn)+ '/fig4_correlations.png', dpi=300,bbox_inches='tight')
print('fig4_correlations.png saved!')
plt.show()

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
#%% Old versions
# Figure 2
# SP vs. HC; CC, path length, CCm1, BC over 8 sparsity levels
if 0:    
    fig,axes = plt.subplots(2,2,figsize=(11,8))
    
    for idxm,labm in enumerate(measl[:4]):
        ax = axes.flatten()[idxm]
        
        for idxg,labg in enumerate(groupl[:2]):
            
            # Prepare data
            data = D[labm+'_'+labg]
            data[data==np.inf] = np.nan
            x = sparsl
            y = np.nanmean(data,0)
            yerr = np.nanstd(data,0)/np.sqrt(np.sum(~np.isnan(data),0))
            
            # Plot data
            if 0:
                x+=idxg*0.01
                ax.errorbar(x,y,yerr,color=pltcolors[idxg],marker=pltmarkers[idxg],capsize=5)
            if 1:
                ysl = scipy.interpolate.interp1d(x,y-yerr)
                ysh = scipy.interpolate.interp1d(x,y+yerr)
                xs = np.linspace(x[0],x[-1],len(sparsl)*5)
                
                ax.plot(x,y,color=pltcolors[idxg],marker=pltmarkers[idxg])
                ax.fill_between(xs,ysl(xs),ysh(xs),alpha=0.25,color=pltcolors[idxg])
            
        # Display settings
        ax.set_xticks(sparsl)
        ax.set_xticklabels(['{:.2f}'.format(x) for x in sparsl])
        ax.set_xlabel('sparsity')
        ax.set_ylabel(labsgm[idxm])
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        xls,yls = ax.get_xlim(), ax.get_ylim()
        ax.text(xls[0]-np.diff(xls)*0.25,yls[1],sublabs[idxm],**sublfont)
        if idxm==1: ax.legend(labsgr[:2],**legendst)
    
    fig.tight_layout() 
    plt.show()
    
    
    #%% Figure 3
    # SP vs. HC; FC whole brain and homotopic
    
    fig,axes = plt.subplots(1,2,figsize=(10,4))
    
    for idxm,labm in enumerate(measl[4:]):
        ax = axes.flatten()[idxm]
        x = np.arange(1,3)  
        # Prepare data
        data = [[x for x in np.squeeze(D[labm+'_'+groupl[0]]) if x not in [np.nan,np.inf]],
                [x for x in np.squeeze(D[labm+'_'+groupl[1]]) if x not in [np.nan,np.inf]]]
        
        # Plot data    
        violins = ax.violinplot(data,showmeans=False,showmedians=False,showextrema=False)
        for idxv,v in enumerate(violins['bodies']): v.set_facecolor(pltcolors[idxv])
        
        qmq=np.asarray([np.percentile(x,[50,25,75]) for x in data])
        ax.scatter(x,qmq[:,0],marker='o',color='white',s=30,zorder=3,label='_nolegend_')
        for idx in x: ax.plot([idx,idx],qmq[idx-1,(1,2)].flatten(),color='k',label='_nolegend_')
    
        xls,yls = ax.get_xlim(), ax.get_ylim()
        if idxm==1:
            yls = [yls[0],yls[1]*1.022]
            ax.set_ylim(top = yls[1])
            esc = np.diff(yls)*0.07
            el = (np.ones((2,))*np.diff(yls)*0.02,np.zeros((2,)))
            ax.errorbar([1,2],np.ones((2,))*yls[1],el,color='k')
            ax.text(1.5,yls[1],'*',ha='center',va='baseline',fontsize=14)
            
        # Display settings
        ax.set_xticks([1,2])
        ax.set_xticklabels(labsgr[:2])
        
        ax.set_ylabel(labsgm[idxm+4])
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.text(xls[0]-np.diff(xls)*0.25,yls[1],sublabs[idxm],**sublfont)
        if idxm==1: ax.legend(labsgr[:2],bbox_to_anchor=(1.02, 1),frameon=False)
        
    fig.tight_layout() 
    plt.show()
    
    #%% Figure 4
    # LL vs. RL vs. HC; CC, path length, CCm1, BC over 8 sparsity levels
    
    fig,axes = plt.subplots(2,2,figsize=(11,8))
    
    for idxm,labm in enumerate(measl[:4]):
        ax = axes.flatten()[idxm]
        
        for idxg,labg in enumerate(groupl[1:]):
            
            # Prepare data
            data = D[labm+'_'+labg]
            data[data==np.inf] = np.nan
            x = sparsl
            y = np.nanmean(data,0)
            yerr = np.nanstd(data,0)/np.sqrt(np.sum(~np.isnan(data),0))
            
            # Plot data
            if 0:
                x+=idxg*0.01
                ax.errorbar(x,y,yerr,color=pltcolors[idxg+1],marker=pltmarkers[idxg+1],capsize=5)
            if 1:
                ysl = scipy.interpolate.interp1d(x,y-yerr)
                ysh = scipy.interpolate.interp1d(x,y+yerr)
                xs = np.linspace(x[0],x[-1],len(sparsl)*5)
                
                ax.plot(x,y,color=pltcolors[idxg+1],marker=pltmarkers[idxg+1])
                ax.fill_between(xs,ysl(xs),ysh(xs),alpha=0.25,color=pltcolors[idxg+1],label='_nolegend_')
            
        # Display settings
        ax.set_xticks(sparsl)
        ax.set_xticklabels(['{:.2f}'.format(x) for x in sparsl])
        ax.set_xlabel('sparsity')
        ax.set_ylabel(labsgm[idxm])
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        xls,yls = ax.get_xlim(), ax.get_ylim()
        ax.text(xls[0]-np.diff(xls)*0.25,yls[1],sublabs[idxm],**sublfont)
        
        if idxm==1: 
            ax.legend(labsgr[1:],**legendst)
            # handles, labels = ax.get_legend_handles_labels()
            # order = [1,2,0]
            # ax.legend([handles[idx] for idx in order],[labels[idx] for idx in order])
            
    fig.tight_layout() 
    plt.show()
    
    
    #%% Figure 5
    # LL vs. RL vs. HC; FC whole brain and homotopic
    
    fig,axes = plt.subplots(1,2,figsize=(10,4))
    
    for idxm,labm in enumerate(measl[4:]):
        ax = axes.flatten()[idxm]
        x = np.arange(1,4)
            
        # Prepare data
        data = [[x for x in np.squeeze(D[labm+'_'+x]) if x not in [np.nan,np.inf]] for x in groupl[1:]]
        
        # Plot data    
        violins = ax.violinplot(data,showmeans=False,showmedians=False,showextrema=False)
        for idxv,v in enumerate(violins['bodies']): v.set_facecolor(pltcolors[idxv+1])
        
        qmq=np.asarray([np.percentile(x,[50,25,75]) for x in data])
        ax.scatter(x,qmq[:,0],marker='o',color='white',s=30,zorder=3,label='_nolegend_')
        for idx in x: ax.plot([idx,idx],qmq[idx-1,(1,2)].flatten(),color='k',label='_nolegend_')
            
        xls,yls = ax.get_xlim(), ax.get_ylim()
        if idxm==1:
            yls = [yls[0],yls[1]*1.022]
            ax.set_ylim(top = yls[1])
            esc = np.diff(yls)*0.07
            el = (np.ones((2,))*np.diff(yls)*0.02,np.zeros((2,)))
            ax.errorbar([1,3],np.ones((2,))*yls[1],el,color='k')
            ax.text(2,yls[1],'*',ha='center',va='baseline',fontsize=14)
            
            ax.errorbar([1,1.9],np.ones((2,))*yls[1]-esc,el,color='k')
            ax.text(1.5,yls[1]-esc,'*',ha='center',va='baseline',fontsize=14)
            ax.errorbar([2.1,3],np.ones((2,))*yls[1]-esc,el,color='k')
            ax.text(2.5,yls[1]-esc,'*',ha='center',va='baseline',fontsize=14)
            
        # Display settings
        ax.set_xticks([1,2,3])
        ax.set_xticklabels(labsgr[1:])
        
        ax.set_ylabel(labsgm[idxm+4])
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.text(xls[0]-np.diff(xls)*0.25,yls[1],sublabs[idxm],**sublfont)
        if idxm==1: ax.legend(labsgr[1:],bbox_to_anchor=(1.02, 1),frameon=False)
        
    fig.tight_layout() 
    plt.show()







