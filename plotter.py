# -*- coding: utf-8 -*-
"""
Created on Thu Feb 19 11:11:27 2026

    Plot functions for simulator results    

@author: Erik Bashore
"""

from main import integrate
import numpy as np 
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib import cm
from true_diagrams import true_cross_section


# LaTex plot setup
plt.rcParams.update({
    "text.usetex": True,        
    "font.family": "serif",
    "font.size": 12,
    "text.latex.preamble": r"\usepackage{amsmath}"
})

# plt.rcParams.update(plt.rcParamsDefault)

def plot_splices(result, variable = 's'):
    
    s_min, s_max = result['Data']['s min'], result['Data']['s max']
    n_PS = result['Data']['PS points']
    n_angle_points = result['Data']['Angle points']
    
    t_range_rough = np.linspace(0, np.pi, n_angle_points + 1)[:-1]
    s_range_rough = np.linspace(s_min, s_max, n_PS)
    
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize= (8, 5), sharex = True)
    
    if variable == 's':
        for t, theta in enumerate(t_range_rough):
            ax1.scatter(s_range_rough, result['Average'][:, t], color = 'purple', alpha = 1/(t + 1), label = fr'$\theta = {theta}$')
            ax1.plot(np.linspace(s_min, s_max, 1000), result['True s plot'][:, t], 
                     color = 'purple', alpha = 1/(t + 1), linestyle = '--')
            
            ax2.scatter(s_range_rough, result['Average interference'][:, t], 
                        color = 'purple', alpha = 1/(t + 1))
            ax2.plot(np.linspace(s_min, s_max, 1000), result['True s plot interference'][:, t], 
                     color = 'purple', alpha = 1/(t + 1), linestyle = '--')
                    
        ax2.set_xlabel(r'$\sqrt{s}$')
        
    if variable == 'theta':
        for k, s_sqrt in enumerate(s_range_rough):
            ax1.scatter(t_range_rough, result['Average'][k, :], color = 'purple', alpha = 1/(k + 1), label = r'$\sqrt{s}$' + f'= {s_sqrt}')
            ax1.plot(np.linspace(0, np.pi, 1000), result['True theta plot'][k, :], 
                     color = 'purple', alpha = 1/(k + 1), linestyle = '--')
            
            ax2.scatter(t_range_rough, result['Average interference'][k, :], 
                        color = 'purple', alpha = 1/(k + 1))
            ax2.plot(np.linspace(0, np.pi, 1000), result['True theta plot interference'][k, :], 
                     color = 'purple', alpha = 1/(k + 1), linestyle = '--')
                    
        ax2.set_xlabel(r'$\theta$')
    
    
    if result['Data']['Particles'] == 2:
        ax1.set_ylabel(r'$|\mathcal{M}_\gamma + \mathcal{M}_Z|^2$')
        ax2.set_ylabel(r'Int$(\mathcal{M}_\gamma, \mathcal{M}_Z)$')
    else:
        ax1.set_ylabel(r'$|\mathcal{M}_\gamma + \mathcal{M}_Z + \mathcal{M}_{Z^\prime}|^2$')
        ax2.set_ylabel(r'Int$(\mathcal{M}_\gamma + \mathcal{M}_Z, \mathcal{M}_{Z^\prime})$')
    
    ax1.legend()
    plt.tight_layout()
    plt.show()



#%%


def plot_surfaces(result):
    
    # ----- Extract information from result ----- #
    result_ave = result['Average']
    error_matrix = result['Error matrix']
    true_matrix = result['True matrix']
    
    result_ave_int = result['Average interference']
    error_matrix_int = result['Error matrix interference']
    true_matrix_int = result['True matrix interference']
    
    s_min, s_max = result['Data']['s min'], result['Data']['s max']
    n_PS = result['Data']['PS points']
    n_angle_points = result['Data']['Angle points']
    
    t_range_rough = np.linspace(0, np.pi, n_angle_points + 1)[:-1]
    s_range_rough = np.linspace(s_min, s_max, n_PS)
    
    
    # ----- Create figure ----- #
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize = (12, 5), subplot_kw = dict(projection='3d'))
    
    
    # ----- Error matrix plot ----- # 
    
    # Color map
    cmap = plt.cm.Reds.copy()
    cmap.set_bad(color = 'gray', alpha = 0.5)  # Gray for masked values
    
    # Define shifted grid for error matrices   
    s_edges = np.linspace(s_min, s_max, n_PS + 1)               # + shift
    t_edges = np.linspace(0, np.pi, n_angle_points + 2) [:-1]   # + shift
    
    T_edge, S_edge = np.meshgrid(t_edges, s_edges)
    
    # Normalize error for coloring
    norm_err = (error_matrix - np.min(error_matrix)) / np.ptp(error_matrix)
    norm_err_int = (error_matrix_int - np.min(error_matrix_int)) / np.ptp(error_matrix_int)
    
    
    z_offset = np.min(true_matrix) - 0.5 * np.ptp(true_matrix)  # Make errors lie below surface plot
    ax1.plot_surface(S_edge, T_edge, np.full_like(S_edge, z_offset),
                     facecolors = cmap(norm_err),
                     shade=False)
    
    
    z_offset = np.min(true_matrix_int) - 0.5 * np.ptp(true_matrix_int)
    ax2.plot_surface(S_edge, T_edge, np.full_like(S_edge, z_offset),
                     facecolors = cmap(norm_err_int),
                     shade=False)
    
    
    # ----- True surface ----- #
    # Define rough grid
    T, S = np.meshgrid(t_range_rough, s_range_rough)
    ax1.plot_surface(S, T, true_matrix,
                            alpha = 0.2, color = 'black', label = 'True value')
    ax2.plot_surface(S, T, true_matrix_int,
                            alpha = 0.2, color = 'black', label = 'True value')
    

    # ----- Scatter plot ----- #
    mask = np.ma.getmask(error_matrix) == False  # Masked array has true if masked, we want opposite
    mask_inv = mask == False
    
    ax1.scatter(S[mask].ravel(), T[mask].ravel(), result_ave[mask].ravel(),
        c = 'purple', label = 'Circuit output')
    ax1.scatter(S[mask_inv].ravel(), T[mask_inv].ravel(), result_ave[mask_inv].ravel(),
        c = 'black', label = 'Discarded')
    
    
    mask = np.ma.getmask(error_matrix_int) == False
    mask_inv = mask == False
    
    ax2.scatter(S[mask].ravel(), T[mask].ravel(), result_ave_int[mask].ravel(),
        c = 'purple',label = 'Circuit output')
    ax2.scatter(S[mask_inv].ravel(), T[mask_inv].ravel(), result_ave_int[mask_inv].ravel(),
        c = 'black', label = 'Discarded')

    
    # ----- Labels ----- #
    ax1.set_xlabel(r"$\sqrt{s}$" + '  [GeV]')
    ax1.set_ylabel(r"$\theta$")
    
    ax2.set_xlabel(r"$\sqrt{s}$" + '  [GeV]')
    ax2.set_ylabel(r"$\theta$")
    
    
    if result['Data']['Particles'] == 2:
        ax1.set_zlabel(r"$|\mathcal{M}_\gamma + \mathcal{M}_Z|^2$")
        ax2.set_zlabel(r"Int$(\mathcal{M}_\gamma, \mathcal{M}_Z)$")
      
    elif result['Data']['Particles'] == 3:
        ax1.set_zlabel(r"$|\mathcal{M}_\gamma + \mathcal{M}_Z + \mathcal{M}_{Z^\prime}|^2$")
        ax2.set_zlabel(r"Int$(\mathcal{M}_\gamma + \mathcal{M}_Z, \mathcal{M}_{Z^\prime})$")
    
    
    # Error bar on the side
    mappable = cm.ScalarMappable(cmap='Reds')  # Define colormap
    mappable.set_array(error_matrix)
    fig.colorbar(mappable, ax = ax1, shrink = 0.5, pad = 0.1, label = r"Relative error \%")
    
    mappable = cm.ScalarMappable(cmap='Reds')
    mappable.set_array(error_matrix_int)
    fig.colorbar(mappable, ax = ax2, shrink = 0.5, pad = 0.1, label = r"Relative error \%")
    
    fig.tight_layout()
    ax2.legend()
    plt.show()



#%%


def plot_errors(result):
    # Error matrices
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize = (10, 5), sharey = True)
    
    s_min, s_max = result['Data']['s min'], result['Data']['s max']
    
    # Color map
    cmap = plt.cm.Reds.copy()
    cmap.set_bad(color='gray', alpha = 0.5)  # Gray for masked values
    
    # Full term
    error_matrix = result['Error matrix']
    
    mat = ax1.imshow(error_matrix, extent = [0, np.pi, s_max, s_min], aspect = 'auto', cmap = cmap)
    plt.colorbar(mat, label = r'Relative error \%')
    
    
    # Interference
    error_int_matrix = result['Error matrix interference']
    
    mat_int = ax2.imshow(error_int_matrix, extent = [0, np.pi, s_max, s_min], aspect = 'auto', cmap = cmap)
    plt.colorbar(mat_int, label = r'Relative error \%')
    

    # Labels    
    ax1.set_xlabel(r'$\theta$')
    ax2.set_xlabel(r'$\theta$')
    ax1.set_ylabel(r'$\sqrt{s}$' + '  [GeV]')
    
    # Titles
    ax1.set_title('Full term')
    ax2. set_title('Interference')
    
    
    fig.tight_layout()
    plt.show()


#%%



def plot_cross_section(result):
    
    n_PS = result['Data']['PS points']
    s_min, s_max = result['Data']['s min'], result['Data']['s max']
    
    n_particles = result['Data']['Particles']
    
    s_range_rough = np.linspace(s_min, s_max, n_PS)
    s_range_fine = np.linspace(s_min, s_max, 1000)
    
    # Integrate over theta and phi
    I = integrate(result)
    cross_section = 1/(64*np.pi**2) * np.array([(1/(s_range_rough[k]**2)) * I[k] for k in range(n_PS)])
    
    plt.scatter(s_range_rough, cross_section, color = 'purple', label = 'Circuit output')
    
    plt.plot(s_range_fine, [true_cross_section(s_sqrt**2) for s_sqrt in s_range_fine], 
             linestyle = '--', color = 'black', alpha = 0.5, label = 'True value')
    
    plt.xlabel(r'$\sqrt{s}$')
    
    if n_particles == 2:
        plt.ylabel(r'$\sigma\big(q\bar{q}\to \gamma, Z \to \mu^-\mu^+\big)$')
        
    elif n_particles == 3:
        plt.ylabel(r'$\sigma\big(q\bar{q}\to \gamma, Z, Z^\prime \to \mu^-\mu^+\big)$')
    
    plt.yscale('log')
    plt.tight_layout()
    plt.show()
    

    
    
    
    
    
    

