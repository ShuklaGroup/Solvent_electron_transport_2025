import numpy as np
import matplotlib.pyplot as plt
import os
import mdtraj as md
import pickle
from tqdm import tqdm
import gc
from deeptime.util import energy2d
import warnings
import mdtraj as md
from deeptime.decomposition import TICA
import pandas as pd
import sys
import seaborn as sns
warnings.filterwarnings('ignore')
plt.rcParams.update({'font.size': 24})

from matplotlib.colors import ListedColormap
purple=tuple(['#6247aa','#815ac0','#a06cd5','#b185db','#d2b7e5'])
blue=tuple(['#2c7da0','#468faf','#61a5c2','#89c2d9','#a9d6e5'])
green=tuple(['#718355', '#87986a', '#97a97c', '#a3b18a', '#cfe1b9'])
orange=tuple(['#ffb700','#ffc300','#ffd000','#ffdd00','#ffea00'])
red=tuple(['#f25c54', '#f27059', '#f4845f', '#f79d65', '#f7b267'])
larger=tuple(['#f7b267','#f7b267','#f7b267','#f7b267','#f7b267'])
anh_colors = purple+blue+green+orange+red#+larger
anh_cmap = ListedColormap(anh_colors)


def plot_1d_kde(feats,solv,bond_type,color,weights=None,sys='MAAAM'):
    # Set global font size
    plt.rcParams.update({'font.size': 20, 'axes.grid': False})
    
    print(len(feats),feats[0].shape)
    if weights is None:
        weights = np.ones((len(np.concatenate(feats))))
    
    # Create the dataframe for plotting
    df = pd.DataFrame({'data': np.concatenate(feats).squeeze() * 10, 'weights': weights})

    # Create the plot
    plt.figure(figsize=(8, 6))
    sns.kdeplot(x='data', weights='weights',color=color, data=df, fill=True, thresh=0, bw_adjust=0.8, alpha=0.4)
    
    # Add title and labels with font size explicitly set
    #plt.title(f'{sys} in {solv} ({bond_type}) OH-bond', fontsize=30)
    plt.xlabel("Distance (Å)")
    plt.ylabel('Probability Density')
    plt.ylim(0,1.0)
    plt.xlim(0,12)
    
    # Ensure tick labels also have the correct font size
    plt.tick_params(axis='both', which='major')
    plt.axvline(x = 3.5, linestyle = '--',color='r',linewidth=3)
    
    # Display the plot
    
    plt.tight_layout()
    plt.savefig(f'plots/HO_plots/2_{sys}_{solv}_{bond_type}_OH.png',dpi=300)



def plot_1d_hist(feats, solv, bond_type, color, weights=None, sys='MAAAM'):
    # Set global font size
    plt.rcParams.update({'font.size': 20, 'axes.grid': False})
    
    data = np.concatenate(feats).squeeze() * 10  # Convert to Å

    if weights is None:
        weights = np.ones_like(data)
    
    # Create the histogram plot
    plt.figure(figsize=(8, 6))
    plt.hist(
        data,
        bins=60,
        range=(0, 12),
        weights=weights,
        density=True,       # Normalizes the histogram to be a probability density
        color=color,
        alpha=0.4,
        edgecolor='black'
    )

    # Axes and formatting
    plt.xlabel("Distance (Å)")
    plt.ylabel("Probability Density")
    plt.xlim(0, 12)
    plt.ylim(0, 1.0)
    plt.axvline(x=3.5, linestyle='--', color='r', linewidth=3)

    # Save and show
    plt.tight_layout()
    plt.savefig(f'plots/HO_plots/{sys}_{solv}_{bond_type}_OH_hist.png', dpi=300)


### FEAUTRIZATION 

#solv = ['water']#, 'acn', 'glycerol', 'tfe']
solv = ['acn', 'glycerol', 'tfe']
#solv_mol = ['HOH']#,'ACN','MGL','TFE']
solv_mol = ['ACN','MGL','TFE']

for s,mol in zip(solv,solv_mol):
    feats = []
    H_14, H_15, H_25 = [], [], []

    for rep in tqdm(range(20)):
        parm_path = f'../{s}/MAAAM/parm7/MAAAM_{rep}.parm7' 
        traj_path = f'../{s}/MAAAM/sim/prod2/MAAAM_{rep}.dcd'
        traj = md.load(traj_path,top=parm_path)
        
        pep_atoms = traj.top.select(f'resname != {mol}')
              
        traj.atom_slice(atom_indices=pep_atoms,inplace=True)
        traj.center_coordinates()

        _, psi = md.compute_psi(traj)
        _, phi = md.compute_phi(traj)
        O_i_H_14 = np.hstack((traj.topology.select('name O and resid 0'), traj.topology.select('name H and resid 3'))).reshape(-1, 2)
        O_i_H_15 = np.hstack((traj.topology.select('name O and resid 0'), traj.topology.select('name H and resid 4'))).reshape(-1, 2)
        O_i_H_25 = np.hstack((traj.topology.select('name O and resid 1'), traj.topology.select('name H and resid 4'))).reshape(-1, 2)
        H_14_distance = md.compute_distances(traj, O_i_H_14)
        H_15_distance = md.compute_distances(traj, O_i_H_15)
        H_25_distance = md.compute_distances(traj, O_i_H_25)
        tot_feats = np.hstack((psi,phi))#,H_14_distance,H_15_distance,H_25_distance))
        feats.append(tot_feats)
        H_14.append(H_14_distance)
        H_15.append(H_15_distance)
        H_25.append(H_25_distance)

    pickle.dump(feats,open(f'traj_feats/2feats_{s}.pkl','wb')) 
    pickle.dump(H_14,open(f'traj_feats/2feat_H_14_{s}.pkl','wb'))
    pickle.dump(H_15,open(f'traj_feats/2feat_H_15_{s}.pkl','wb'))
    pickle.dump(H_25,open(f'traj_feats/2feat_H_25_{s}.pkl','wb'))


### SAVING H-BOND PLOTS

#solv = ['water']#, 'acn', 'glycerol', 'tfe']
solv = [ 'acn', 'glycerol', 'tfe']
colors = [sns.color_palette("deep")[1]]#, sns.color_palette("deep")[2],sns.color_palette("deep")[3], sns.color_palette("deep")[0]] 
#water_weights = pickle.load(open("water_weights.pkl","rb"))
acn_weights = pickle.load(open("acn_weights.pkl","rb"))
glycerol_weights = pickle.load(open("glycerol_weights.pkl","rb"))
tfe_weights = pickle.load(open("tfe_weights.pkl","rb"))
weights = [water_weights, acn_weights, glycerol_weights, tfe_weights]

for s,c,w in zip(solv,colors,weights):
    feat = pickle.load(open(f'traj_feats/2feat_H_14_{s}.pkl','rb'))
    plot_1d_kde(feat,s,bond_type='1 to 4',color=c)#, weights=w)
    
    feat = pickle.load(open(f'traj_feats/2feat_H_15_{s}.pkl','rb'))
    plot_1d_kde(feat,s,bond_type='1 to 5',color=c)#, weights=w)

    feat = pickle.load(open(f'traj_feats/2feat_H_25_{s}.pkl','rb'))
    plot_1d_kde(feat,s,bond_type='2 to 5',color=c)#, weights=w)

####SAVING DIHEDRAL PLOTS

for s, mol, c in zip(solv, solv_mol, colors):
    for rs in [1, 2, 3]:
        phi_list = []
        psi_list = []
        
        for rep in tqdm(range(20)):
            parm_path = f'../{s}/MAAAM/parm7/MAAAM_{rep}.parm7'
            traj_path = f'../{s}/MAAAM/sim/prod/MAAAM_{rep}.dcd'
            
            traj = md.load(traj_path, top=parm_path)
            pep_atoms = traj.top.select(f'resname != {mol}')
            traj.atom_slice(atom_indices=pep_atoms, inplace=True)
            traj.center_coordinates()
            
            phi_atoms = traj.top.select(f'(resid {rs-1} and name C) or (resid {rs} and (name N or name CA or name C))')
            psi_atoms = traj.top.select(f'(resid {rs} and (name N or name CA or name C)) or (resid {rs+1} and name N)')
            
            phi_atoms = phi_atoms.reshape((1, 4))
            psi_atoms = psi_atoms.reshape((1, 4))
                            
            phi_feat = md.compute_dihedrals(traj, phi_atoms)
            psi_feat = md.compute_dihedrals(traj, psi_atoms)
            
            phi_list.append(phi_feat.flatten())
            psi_list.append(psi_feat.flatten())
        
        phi = np.concatenate(phi_list) * 180 / np.pi  # Convert to degrees
        psi = np.concatenate(psi_list) * 180 / np.pi

        # --- Plot 1: Free Energy Plot ---
        energies = energy2d(phi, psi, bins=(50, 50), kbt=2.479/4.184, weights=None, shift_energy=True)
        z1lim = [0, 5]  # Energy range
        levels = np.linspace(*z1lim, 20)

        fig, ax = plt.subplots(figsize=(8, 6))
        ax, contour, cbar = energies.plot(contourf_kws=dict(cmap=anh_cmap),levels=levels)
        ticks = np.linspace(*z1lim, 6)
        cbar.ax.set_yticks(ticks)
        cbar.ax.set_yticklabels(ticks)
        cbar.set_label('Free Energy (kcal/mol)')
        ax.set_xlim(-180, 180)
        ax.set_ylim(-180, 180)
        ax.set_xticks(range(-180, 181, 90))
        ax.set_yticks(range(-180, 181, 90))
        ax.set_xlabel(fr'$\phi_{{res {rs}}}$ (°)')
        ax.set_ylabel(fr'$\psi_{{res {rs}}}$ (°)')
        plt.tight_layout()
        plt.savefig(f'plots/res_dihedral_plots/fep_{s}_res{rs}.png', dpi=300)
        plt.close()

        # --- Plot 2: Scatter Plot ---
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.scatter(phi, psi, c=c, s=100, edgecolor='k', alpha=0.8)
        ax.set_xlim(-180, 180)
        ax.set_ylim(-180, 180)
        ax.set_xticks(range(-180, 181, 90))
        ax.set_yticks(range(-180, 181, 90))
        ax.set_xlabel(fr'$\phi_{{res {rs}}}$ (°)')
        ax.set_ylabel(fr'$\psi_{{res {rs}}}$ (°)')
        ax.grid(True)
        plt.tight_layout()
        plt.savefig(f'plots/res_dihedral_plots/scatter_dihed_{s}_res{rs}.png', dpi=300)
        plt.close()


## SAVING TICA PLOTS


for s in solv:
    
    tic_lagtime = 10
    feat_list = pickle.load(open(f'traj_feats/2feats_{s}.pkl','rb'))
    tica = TICA(lagtime=tic_lagtime).fit_fetch(feat_list)
    projection = [tica.transform(x) for x in feat_list]
    concat_projection = np.concatenate(projection)
    
    energies = energy2d(concat_projection[:,0], concat_projection[:,1], bins=(50, 50), kbt=2.479/4.184, weights=None, shift_energy=True)
    z1lim = [0, 5]
    levels=np.linspace(*z1lim,20)
    plt.figure(figsize=(8, 6))
    ax, contour, cbar = energies.plot(contourf_kws=dict(cmap=anh_cmap),levels=levels)
    ticks = np.linspace(*z1lim, 6)
    cbar.ax.set_yticks(ticks)
    cbar.ax.set_yticklabels(ticks)
    cbar.ax.tick_params(axis='y')#,labelsize=14)
    cbar.set_label('Free Energy (kcal/mol)')#,fontsize=14)
    ax.tick_params(axis='both')#,labelsize=14)
    ax.set_xlabel(f'tIC-1')
    ax.set_ylabel(f'tIC-2')
    plt.tight_layout()
    #plt.title(f"{sys} in {solv}")
    plt.savefig(f'plots/tic_plots/2{s}_tics.png',dpi=300)
    plt.close()



### SAVING TICA MAPPED TO H-BOND


num_bins = 50  # You can change this value to set the desired number of bins

bonds = ['14','15','25']
bond_desc = ['1 to 4', '1 to 5', '2 to 5']
for bond,desc in zip(bonds,bond_desc):
    for s in solv:
        metric_list = []  # List to store the separate metric values
        feat_list = pickle.load(open(f'traj_feats/2feats_{s}.pkl','rb'))
        ho_feats = pickle.load(open(f'traj_feats/2feat_H_{bond}_{s}.pkl','rb'))
        h_feat = np.concatenate(ho_feats)
        tic_lagtime = 10

        tica = TICA(lagtime=tic_lagtime).fit_fetch(feat_list)
        projection = [tica.transform(x) for x in feat_list]
        concat_projection = np.concatenate(projection)

        # Assuming you want to create a 2D histogram or density estimation of the separate metric
        metric_values = h_feat * 10
        metric_values = metric_values.flatten()

        # Ensure lengths match
        if len(concat_projection) != len(metric_values):
            metric_values = metric_values[:len(concat_projection)]  # Adjust if necessary

        # Create a 2D histogram for binning
        binned_metric_sum, x_edges, y_edges = np.histogram2d(
            concat_projection[:, 0], concat_projection[:, 1],
            bins=num_bins,
            weights=metric_values
        )
        binned_count, _, _ = np.histogram2d(
            concat_projection[:, 0], concat_projection[:, 1],
            bins=num_bins
        )

        # Calculate average metric value in each bin
        average_metric = np.divide(binned_metric_sum, binned_count, out=np.zeros_like(binned_metric_sum), where=binned_count != 0)

        # Mask zero values
        masked_average_metric = np.ma.masked_where(average_metric == 0, average_metric)

        # Create a contour plot using the average metric
        z1lim = [0, 10]
        levels = np.linspace(*z1lim, 20)
        
        fig, ax = plt.subplots()
        fig.set_size_inches(8, 6)
        # Create a meshgrid for contour plot
        X, Y = np.meshgrid(x_edges[:-1], y_edges[:-1])

        # Create a custom colormap
        cmap = plt.cm.get_cmap(anh_cmap)  # Get the existing colormap
        cmap.set_bad('white')  # Set the color for NaN values to white

        # Transpose masked_average_metric to match the meshgrid
        contour = ax.contourf(X, Y, masked_average_metric.T, levels=levels, cmap='plasma')  # Transpose for correct orientation
        cbar = plt.colorbar(contour, ax=ax)

        # Set ticks from 0 to 10
        ticks = np.linspace(*z1lim, 6)  # Generates 6 evenly spaced ticks from 0 to 10
        cbar.ax.set_yticks(ticks)
        cbar.ax.set_yticklabels(ticks)
        cbar.ax.tick_params(axis='y')
        cbar.set_label(f'{desc} OH-bond distance (Å)')  # Update label to reflect the metric used

        ax.tick_params(axis='both')
        ax.set_xlabel('tIC-1')
        ax.set_ylabel('tIC-2')

        # Set original axis limits based on concat_projection
        ax.set_xlim(np.min(concat_projection[:, 0]), np.max(concat_projection[:, 0]))
        ax.set_ylim(np.min(concat_projection[:, 1]), np.max(concat_projection[:, 1]))

        # Adjust layout
        plt.axis('tight')  # Tighten the axis
        plt.tight_layout()  # Reduce padding

        #plt.title(f"{sys} in {solv}")
        plt.savefig(f'plots/tic_plots/2{s}_{bond}_tic_dis.png', dpi=300)