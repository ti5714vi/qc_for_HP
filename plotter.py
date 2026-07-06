# -*- coding: utf-8 -*-
"""
Created on Thu Feb 19 11:11:27 2026

    Plot functions for simulator results    

@author: Erik Bashore
"""

from main import integrate, SDandAveComp
import numpy as np 
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib import cm
from true_diagrams import *
from scipy.interpolate import interp1d, PchipInterpolator

from scipy.integrate import quad

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
    
    t_range_fine = np.linspace(0, np.pi, 51)[:-1]
    s_range_fine = np.linspace(s_min, s_max, 50)
    
    # ----- Create figure ----- #
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize = (12, 5), subplot_kw = dict(projection='3d'), 
                                   constrained_layout = True)
    
    
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
    true_matrix = np.array([[true(s_sqrt**2, theta, n_particles = 2) for theta in t_range_fine] \
                            for s_sqrt in s_range_fine])
    true_matrix_int = np.array([[true_int(s_sqrt**2, theta, n_particles = 2) for theta in t_range_fine] \
                            for s_sqrt in s_range_fine])
    
    # Define fine grid
    T_fine, S_fine = np.meshgrid(t_range_fine, s_range_fine)
    ax1.plot_surface(S_fine, T_fine, true_matrix,
                            alpha = 0.2, color = 'black', label = 'True value')
    
    ax2.plot_surface(S_fine, T_fine, true_matrix_int,
                            alpha = 0.2, color = 'black', label = 'True value')
    

    # ----- Scatter plot ----- #
    # Define rough grid
    T, S = np.meshgrid(t_range_rough, s_range_rough)
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
    vmin = min(error_matrix.min(), error_matrix_int.min())
    vmax = min(error_matrix.max(), error_matrix_int.max())
    
    mappable = cm.ScalarMappable(cmap='Reds')  # Define colormap
    mappable.set_array(np.array([vmin, vmax]))
    fig.colorbar(mappable, ax = [ax1, ax2], shrink = 0.5, pad = 0.05, label = r"Relative error \%")
    
    # mappable = cm.ScalarMappable(cmap='Reds')
    # mappable.set_array(error_matrix_int)
    # fig.colorbar(mappable, ax = ax2, shrink = 0.5, pad = 0.1, label = r"Relative error \%")
    
    # fig.tight_layout()
    ax2.legend()
    plt.show()



#%%


def plot_errors(result):
    # Error matrices
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize = (11, 3.5), sharey = True, 
                                        constrained_layout = True)
    
    s_min, s_max = result['Data']['s min'], result['Data']['s max']
    
    # Color map
    cmap = plt.cm.Reds.copy()
    cmap.set_bad(color='gray', alpha = 0.5)  # Gray for masked values
    
    # Full term
    error_matrix = result['Error matrix']
    print('Full term average:', np.mean(error_matrix))
    
    # Interference
    error_int_matrix = result['Error matrix interference']
    print('Interference average:', np.mean(error_int_matrix))
    
    # Universal min and max
    vmin = min(error_matrix.min(), error_int_matrix.min())
    vmax = min(error_matrix.max(), error_int_matrix.max())
    
    mat = ax1.imshow(error_matrix, extent = [0, np.pi, s_max, s_min], 
                     vmin = vmin, vmax = vmax, aspect = 'auto', cmap = cmap)
    # plt.colorbar(mat, label = r'Relative error \%')
    
    
    mat_int = ax2.imshow(error_int_matrix, extent = [0, np.pi, s_max, s_min], 
                         vmin = vmin, vmax = vmax, aspect = 'auto', cmap = cmap)
    
    # Combined colorbar
    cbar = fig.colorbar(mat_int, ax = [ax1, ax2], pad = 0.05, label = r'Relative error \%')
    
    # Relation between
    ratio_matrix = error_int_matrix/error_matrix
    print('Ratio average:', np.mean(ratio_matrix))
    
    # Color map
    cmap = plt.cm.Blues.copy()
    cmap.set_bad(color='gray', alpha = 0.5)  # Gray for masked values

    ratio_mat = ax3.imshow(ratio_matrix, extent = [0, np.pi, s_max, s_min], aspect = 'auto', cmap = cmap)
    plt.colorbar(ratio_mat, label = r'Ratio \%')

    # Labels    
    ax1.set_xlabel(r'$\theta$')
    ax2.set_xlabel(r'$\theta$')
    ax3.set_xlabel(r'$\theta$')
    ax1.set_ylabel(r'$\sqrt{s}$' + '  [GeV]')
    
    # Titles
    ax1.set_title('Full term')
    ax2.set_title('Interference')
    ax3.set_title('Quotient')
    
    # fig.tight_layout()
    plt.show()


#%%



def plot_cross_section(result):
    
    n_PS = result['Data']['PS points']
    n_angle_points = result['Data']['Angle points']
    s_min, s_max = result['Data']['s min'], result['Data']['s max']
    
    n_particles = result['Data']['Particles']
    
    s_range_rough = np.linspace(s_min, s_max, n_PS)
    s_range_fine = np.linspace(s_min, s_max, 1000)
    
    # Integrate over theta and phi
    I = integrate(result)
    cross_section = 1/(64*np.pi**2) * np.array([(1/(s_range_rough[k]**2)) * I[k] for k in range(n_PS)])
    
    
    
    # Compute error propagation
    SDs = np.zeros((n_PS, n_angle_points))

    for i in range(n_PS):
        for k in range(n_angle_points):
            SDs[i, k] = SDandAveComp(result['Full output'][0, :, i, k])[0]


    dSigma = np.zeros(n_PS)

    theta_range = np.linspace(0, np.pi, n_angle_points + 1)[:-1]
    Dtheta = np.pi/n_angle_points

    for i in range(n_PS):
        frontfactor = 2*np.pi/ (64*np.pi**2*s_range_rough[i]**2)
        
        summation = 0
        for k in range(n_PS):
            summation += (Dtheta*np.sin(theta_range[k]))**2 * SDs[i, k]**2
            
        dSigma[i] = frontfactor * np.sqrt(summation)
    
    # Smooth curve
    f = PchipInterpolator(s_range_rough, cross_section)
    f_upper = PchipInterpolator(s_range_rough, cross_section + dSigma)
    f_lower = PchipInterpolator(s_range_rough, cross_section - dSigma)
    
    central = f(s_range_fine)
    err_upper = f_upper(s_range_fine)
    err_lower = f_lower(s_range_fine)
    
    # plt.fill_between(s_range_fine, err_lower, err_upper, color = 'purple', alpha = 0.15)
    
    # Outputs
    barn_conversion = 3.894 * 1e8
    plt.scatter(s_range_rough, cross_section, color = 'purple')
    plt.errorbar(s_range_rough, cross_section, yerr = dSigma, 
                    fmt = 'o', color = 'purple', ecolor = 'lightgray', capsize = 5, label = 'Circuit output')
    
    
    # True value
    plt.plot(s_range_fine, [true_cross_section(s_sqrt**2) for s_sqrt in s_range_fine], 
             linestyle = '--', color = 'black', alpha = 0.5, label = 'True value')
    
    plt.xlabel(r'$\sqrt{s}$')
    
    if n_particles == 2:
        plt.ylabel(r'$\sigma\big(q\bar{q}\to \gamma/Z \to \mu^-\mu^+\big)$')
        
    elif n_particles == 3:
        plt.ylabel(r'$\sigma\big(q\bar{q}\to \gamma/Z/Z^\prime \to \mu^-\mu^+\big)$')
    
    plt.legend()
    plt.yscale('log')
    plt.tight_layout()
    plt.show()
    

    
    
    
#%%


# Interference cross-section
interference = lambda theta, s: get_photon_diagram(s, theta) * get_Z_diagram(s, theta).conj() + get_photon_diagram(s, theta).conj() * get_Z_diagram(s, theta)

s_range = np.linspace(80, 100, 100)
s_range_rough =  np.linspace(80, 100, 8)
cross_section = []

for s_sqrt in s_range:
    I, err = quad(interference, 0, np.pi, args=(s_sqrt**2))
    cross_section.append(I/s_sqrt**2)

cross_section = 1/(32*np.pi) * np.array(cross_section)


integral = integrate(result, term = 'Int')
circ_outs = np.array([1/(64*np.pi**2*s_range_rough[j]**2)  * integral[j]  for j in range(8)])

# pb convert
circ_outs *= 0.3894 * 1e8
cross_section *= 0.3894 * 1e8


# MadGraph 
MG = 1e3 * np.array([-0.1726, -0.2291, -0.3339, -0.5863, 0.2832, 0.5155, 0.283, 0.1921])

plt.scatter(s_range_rough, MG, label = 'MG*1000')
plt.scatter(s_range_rough, circ_outs, color = 'purple', label = 'Circuit output')
plt.plot(s_range, cross_section, color = 'black', alpha = 0.5, linestyle = '--', label = 'True value')


plt.xlabel(r'$\sqrt{s}$ [GeV]')
plt.ylabel(r'$\sigma_{\gamma/Z}$')
plt.legend()
plt.tight_layout()



#%%


s_range = np.linspace(80, 100, 8)
MG_out = np.array([10.02, 16.79, 36.77, 139.1, 778.1, 114.4, 36.55, 18.31])

num_out = lambda s_sqrt: true_cross_section(s = s_sqrt**2) * 0.3894 * 1e8 /(12*np.pi)

quotient = MG_out/num_out(s_range)


plt.plot(s_range, quotient, label = 'Quotient', color = 'black')
plt.scatter(s_range, MG_out, label = 'MadGraph')
plt.scatter(s_range, num_out(s_range), label = 'Python')
# plt.hlines(y = 1, xmin = s_range[0], xmax = s_range[-1])
plt.yscale('log')

plt.legend()






    
    

