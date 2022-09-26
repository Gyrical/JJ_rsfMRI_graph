#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 22 08:17:33 2022

@author: tianlu
"""
import os
import numpy as np
import numpy.matlib
import pandas as pd
import scipy.io
import scipy.interpolate
from scipy.stats import pearsonr
from matplotlib import cm
import matplotlib.pyplot as plt
import seaborn as sns

# Load data
dir_main = os.path.realpath(os.path.join(os.path.dirname(__file__),'..'))
D = scipy.io.loadmat(dir_main+'/results/fc_graph_data.mat')
Dstats = scipy.io.loadmat(dir_main+'/results/fc_graph_stats.mat')
lesionsize = scipy.io.loadmat(dir_main+'/data/patients/lesionsize.mat')['lesionsize'].flatten()

# Study settings
sparsl = np.arange(0.1,0.95,0.05) 
measl  = 'cc cp preCGcc preCGcb fcwb fch'.split()
groupl = 'pat con ll rl'.split() 
labsgm = 'CC$_{avg}$ PL CC$_{preCG}$ BC$_{preCG}$ FC$_{wb}$ FC$_{h}$'.split()
labsgr = ['Stroke patients','Healthy controls','LL','RL']
alpha  = 0.05

# Figure settings # https://matplotlib.org/stable/gallery/color/color_cycle_default.html
pltcolors = [plt.rcParams['axes.prop_cycle'].by_key()['color'][x] for x in [4,7,0,1]]#[:4]
pltmarkers = ['^','o','<','>'] # SP HC LL RL
sublabs  = ['A','','B','','C',''] # sublabs  = 'A B C D E F'.split()
sublfont = {'fontsize':14,'fontweight':'bold','fontname':'Arial'}
legendst = {'bbox_to_anchor':(1.01, 1),'frameon':False}
legendsts = {'bbox_to_anchor':(1.01, 1.01),'frameon':False} # 0.95
scatterst = {'s':30,'zorder':3,'label':'_nolegend_'}
signst1 = {'ha':'center','va':'bottom','fontsize':14}
signst2 = {'ha':'center','va':'baseline','fontsize':14}
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['font.family'] = 'sans-serif'

# Stats settings
sstats = True
# pchoic = 1 # 1: uncorrected, 3: bonferroni, 5: fdr


#%% Fig 2 (2-3) Combined SP vs. HC
    
statsm = [['*' if p < alpha else '' for p in Dstats['stats_all_'+x].flatten()] for x in measl[:4]]

fig,axes = plt.subplots(3,2,figsize=(10,9))
for idxm,labm in enumerate(measl):
    ax = axes.flatten()[idxm]
    
    # ------------------------------------------------------------------------------
    # Show graph measures
    
    if idxm<4:
        
        for idxg,labg in enumerate(groupl[:2]):
            
            # Prepare data
            data = D[labg+'_'+labm]
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
        
        # Prepare data
        data = [[x for x in np.squeeze(D[groupl[0]+'_'+labm]) if x not in [np.nan,np.inf]],
                [x for x in np.squeeze(D[groupl[1]+'_'+labm]) if x not in [np.nan,np.inf]]]
        
        # Plot data    
        violins = ax.violinplot(data,showmeans=False,showmedians=False,showextrema=False)
        for idxv,v in enumerate(violins['bodies']): v.set_facecolor(pltcolors[idxv])
        
        qmq=np.asarray([np.percentile(x,[50,25,75]) for x in data])
        x = np.arange(1,3)
        for idx in range(2):
            ax.scatter(x[idx],qmq[idx,0],marker='o',color='white',**scatterst) # pltmarkers[idx]
            ax.plot([idx+1,idx+1],qmq[idx,(1,2)],color='k',label='_nolegend_') #pltcolors[idx]
        
        # Show stats markers
        if sstats and Dstats['stats_all_'+labm] < alpha:
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

    # Show subplot labels
    xls,yls = ax.get_xlim(), ax.get_ylim()
    ax.text(xls[0]-np.diff(xls)*0.25,yls[1],sublabs[idxm],**sublfont)
    
    # Set legend location
    if idxm==1 or idxm==5: ax.legend(labsgr[:2],loc='upper left',**legendsts)
     
    # General display settings
    ax.set_ylabel(labsgm[idxm])           
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    
fig.tight_layout()
fn_out = 'fig2_largegroups.png' 
plt.savefig(fn_out, dpi=300)
print(fn_out+' saved!')
plt.show()


#%% Figure 3 (4-5) combined LL vs. RL vs. HC

slabs = ['$LL ≠ RL$','$LL ≠ HC$','$RL ≠ HC$']
statsma = [[['*' if x < alpha else '' for x in Dstats['stats_{}_{}'.format(xx,yy)]] 
               for xx in 'll_rl con_ll con_rl'.split()] for yy in measl[:4]] # else '^' if x<.1 

fig,axes = plt.subplots(3,2,figsize=(10,9))
for idxm,labm in enumerate(measl):
    ax = axes.flatten()[idxm]
    
    if idxm < 4:
        for idxg,labg in enumerate(groupl[1:]):
            
            # Prepare data
            if labg+'_'+labm not in D.keys(): continue
            data = D[labg+'_'+labm]
            data[data==np.inf] = np.nan
            x = sparsl
            y = np.nanmean(data,0)
            yerr = np.nanstd(data,0)/np.sqrt(np.sum(~np.isnan(data),0))
            
            # Plot data
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
        data = [[x for x in np.squeeze(D[x+'_'+labm]) if x not in [np.nan,np.inf]] for x in groupl[1:]]
        
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
        if Dstats['stats_con_ll_rl_'+labm]<alpha:
            yls = [yls[0],yls[1]*1.022]
            ax.set_ylim(top = yls[1])
            esc = np.diff(yls)*0.07
            el = (np.ones((2,))*np.diff(yls)*0.02,np.zeros((2,)))
            
            if Dstats['stats_con_rl_'+labm] < alpha:
                ax.errorbar([1,3],np.ones((2,))*yls[1],el,color='k')
                ax.text(2,yls[1],'*',ha='center',va='baseline',fontsize=14)

            if Dstats['stats_con_ll_'+labm] < alpha:          
                ax.errorbar([1,1.9],np.ones((2,))*yls[1]-esc,el,color='k')
                ax.text(1.5,yls[1]-esc,'*',ha='center',va='baseline',fontsize=14)
            if Dstats['stats_ll_rl_'+labm] < alpha:
                ax.errorbar([2.1,3],np.ones((2,))*yls[1]-esc,el,color='k')
                ax.text(2.5,yls[1]-esc,'*',ha='center',va='baseline',fontsize=14)
            
        # Display settings FC
        ax.set_xticks([1,2,3])
        ax.set_xticklabels(labsgr[1:])
        
    # General display settings
    ax.set_ylabel(labsgm[idxm])
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    # Show sublabels
    xls,yls = ax.get_xlim(), ax.get_ylim()
    ax.text(xls[0]-np.diff(xls)*0.25,yls[1],sublabs[idxm],**sublfont)
    # Show legend
    if idxm==1 or idxm==5: 
        handles, labels = ax.get_legend_handles_labels()
        order = [1,2,0]
        ax.legend([handles[idx] for idx in order],[labels[idx] for idx in order],loc='upper left',**legendsts)
    
fig.tight_layout() 
plt.savefig('fig3_subgroups.png', dpi=300)
print('fig3_subgroups.png saved!')
plt.show()


#%% Figure 4 with correlations
import re

if 1:
    # Option 1: show single sparsity

    spars = 0.75
    
    fig,axes = plt.subplots(3,2,figsize=(8,9))
    
    for idxm,labm in enumerate(measl):
        ax = axes.flatten()[idxm]
        
        # Prepare data
        if idxm>3: data = D['pat_'+labm].flatten()
        else: data = D['pat_'+labm][:,np.where(np.isclose(sparsl,spars))[0][0]].flatten() 
        
        # Plot
        ldata = pd.DataFrame({'data':data,'lesionsize':lesionsize})
        g = sns.regplot(x='lesionsize',y='data',data=ldata,ax=ax,marker='^',color=pltcolors[0])
        ax.set_xlim([-10,370])

    
        xls,yls = ax.get_xlim(), ax.get_ylim()
        if idxm<4: ax.text(xls[0]+np.diff(xls)*0.05,yls[1],'$sparsity = {:.2f}$'.format(spars),ha='left',va='bottom',fontsize=8)
            
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


else:
    # Option 2: show all sparsities

    cmap = np.flip(cm.get_cmap('Purples')(np.arange(0.2,1,0.1)),axis=0)
    
    # Prepare data
    dfs=[]    
    for idxm,labm in enumerate(measl):
        data  = D['pat_'+labm]
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
plt.savefig('fig4_correlations.png', dpi=300,bbox_inches='tight')
print('fig4_correlations.png saved!')
plt.show()


# %% Graphical abstract

def plot_connectome(ax,meas='cc'):
    from nilearn import plotting
    flip = meas=='bc'
  
    node_coords_ipsi = [[16,48,0],
                       [38,35,0],
                       [12,22,0],
                       [28,16,0],
                       [45,15,0],
                       [18,2,0],
                       [52,-5,0],
                       [35,-20,0], # M1
                       [57,-30,0],
                       [11,-30,0],
                       [45,-45,0],
                       [18,-58,0],
                       [38,-75,0]]
    N_ipsi = len(node_coords_ipsi)
    
    node_coords_contra = [[-x,y,z] for idx,(x,y,z) in enumerate(node_coords_ipsi) if idx in [0,1,2,9,11,12]] # selected nodes
    node_coords = np.vstack((node_coords_ipsi,node_coords_contra))
    N = len(node_coords)
    # print(f'N nodes: {N}')

    # Define groups and connections
    nodegroups = {'m1':7, 
                  'aroundm1':[5,6,8,9,10]}
    nodeconn = {'bc':{'wholebrain':[(0,2),(1,4),(2,5),(3,4),(4,5),(4,6),(5,7),(7,8),(7,9),(8,10),(10,11),(9,11),(11,12)],
                      'homotopic':[(0,N_ipsi),(1,N_ipsi+1),(2,N_ipsi+2),(9,N_ipsi+3),(11,N_ipsi+4)],
                      'ipsilesional':[(0,2),(3,5)]},
                'cc': {'wholebrain':[(0,2),(1,4),(2,5),(3,4),(4,5),(4,6),(5,7),(7,8),(7,9),(8,10),(10,11),(9,11),(11,12)],
                       'homotopic':[(2,N_ipsi+2)],
                       'ipsilesional':[(0,2),(1,2),(2,3),(3,4),(3,5),(4,5)]}
                }
    
    node_stroke = [[-35,-10,0]]
    
    # Flip nodes
    if flip: 
        node_coords = [[-x,y,z] for x,y,z in node_coords]
        node_stroke *= np.asarray([-1,1,1]) 
    
    
    # Grey adjacency matrix
    adjacency_matrix = np.zeros((N,N))
    # Check sequence
    # for idxn in range(N-1): adjacency_matrix[idxn,idxn+1] = 1
    for idxn in nodeconn[meas]['wholebrain']: adjacency_matrix[idxn] = 1
    for idxn in nodeconn[meas]['homotopic']: adjacency_matrix[idxn] = 1
    for idxn in nodeconn[meas]['ipsilesional']: adjacency_matrix[idxn[0]+N_ipsi,idxn[1]+N_ipsi] = 1
    adjacency_matrix += adjacency_matrix.T
    
    nodeconn_hi = {'bc':[(0,2),(2,5),(5,7),(7,8),(8,10),(10,11),(11,12),
                                (1,N_ipsi+1),(1,4),(4,5),(7,9),(9,11),
                                (9,N_ipsi+3),(N_ipsi+3,N_ipsi+5)],
                   'cc':[(4,5),(4,6),(4,9),(5,6),(5,8),(5,9),(5,10),(6,8),(6,9),(6,10),(8,9),(8,10),(9,10)]}
    adj_mat_hi = np.zeros((N,N))
    for idxn in nodeconn_hi[meas]: adj_mat_hi[idxn] = 1
    adj_mat_hi += adj_mat_hi.T
    
    # Display settings
    node_size = [850 if idx==nodegroups['m1'] else 150 for idx in range(N)]
    node_mark = 'o'
    node_color = 'grey'
    
    # Plot stuff
    plotting.plot_connectome(adjacency_matrix=adjacency_matrix, node_coords=node_coords, 
                             node_size = node_size, node_color=node_color, node_kwargs={'edgecolor':'k','marker':node_mark},
                             edge_kwargs = {'color':'#090909','alpha':0.65},
                             display_mode='z', annotate=False,axes=ax)
    
    # Show highlghted connections
    plotting.plot_connectome(adjacency_matrix=adj_mat_hi, node_coords=node_coords, 
                             node_size = node_size, node_color=node_color, node_kwargs={'edgecolor':'k','marker':node_mark},
                             edge_kwargs = {'color':'green','alpha':0.75},
                             display_mode='z', annotate=False,axes=ax)
    
    # Show Stroke lesion    
    plotting.plot_connectome(adjacency_matrix=np.asarray([[1]]), node_coords=node_stroke, 
                             node_size=6500,node_color=node_color, node_kwargs={'edgecolor':'k','marker':(7,1)},
                             display_mode='z', annotate=False,axes=ax)
    
    # Show markers as nodes
    if False:
        plotting.plot_connectome(adjacency_matrix=np.asarray([[1]]), node_coords=node_stroke, 
                                 node_size=500,node_color='k', node_kwargs={'edgecolor':'k','marker':'$S$'},
                                 display_mode='z', annotate=False,axes=ax)
        plotting.plot_connectome(adjacency_matrix=np.asarray([[1]]), node_coords=np.asarray([node_coords[nodegroups['m1']]]),
                                 node_size=250,node_color='k', node_kwargs={'edgecolor':'k','marker':'$M1$'},
                                 display_mode='z', annotate=False,axes=ax)

    
    return ax


fig,axes = plt.subplots(1,2,figsize=(10,6))
plot_connectome(axes[0],'cc')
plot_connectome(axes[1],'bc')
plt.show()
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    








