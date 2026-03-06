# -*- coding: utf-8 -*-
"""
Created on Sat Dec 13 11:16:36 2025

    Constant inputs such as quantum numbers and particle masses

@author: erikb
"""

import numpy as np


# Scattering angle
theta = np.pi/3

# Weak mixing angle
theta_W = np.arcsin(0.47143025548)


# Up quark inputs
m_q = 2.2 # GeV
Q_q = 2/3 # eV
I_q = 1/2 # eV

# Muon inputs
m_mu = 0.1057 # GeV
Q_mu = -1     # eV
I_mu = -1/2   # eV

# Electron charge
alpha_inv = 1/1.325070e+02
e = np.sqrt(4*np.pi*alpha_inv) # eV
Q_e = -1
I_e = -1/2 

# Z boson inputs
M_Z = 91.188      # GeV
Gamma_Z = 2.441404 # GeV

# Complex mass
M_Z = M_Z  #- 1j * Gamma_Z/2


# Z prime boson inputs
M_Z_prime = 1000      # GeV
Gamma_Z_prime = 2.2 # GeV

M_Z_prime = M_Z_prime - 1j * Gamma_Z_prime/2

# Couplings as functions of Q, I
C_V_Z = lambda Q, I: 1j*e/(2*np.cos(theta_W)*np.sin(theta_W)) *(I - 2*Q*np.sin(theta_W)**2)  # Vector coupling
C_A_Z = lambda Q, I: 1j*e/(2*np.cos(theta_W)*np.sin(theta_W)) * I                           # Axial coupling

C_L_Z = lambda Q, I: 1j*e/(2*np.cos(theta_W)*np.sin(theta_W))*(2*I - 2*Q*np.sin(theta_W)**2)
C_R_Z = lambda Q, I: 1j*e/(2*np.cos(theta_W)*np.sin(theta_W))*(- 2*Q*np.sin(theta_W)**2)

# Photon inputs
M_A = 0 # GeV

# Couplings
C_V_A = lambda Q: -1j*e*Q
C_A_A = 0

# Gauge parameter
xi = 1


# Collect all couplings
# Cs_q = np.array([[C_V_A(Q_q), C_A_A], [C_V_Z(Q_q, I_q), C_A_Z(Q_q, I_q)]])         # First vertex with quarks

Cs_q = np.array([[C_V_A(e), C_A_A], [C_V_Z(e, I_e), C_A_Z(e, I_e)]])         # First vertex with quarks
Cs_mu = np.array([[C_V_A(Q_mu), C_A_A], [C_V_Z(Q_mu, I_mu), C_A_Z(Q_mu, I_mu)]])   # Second vertex with muons

q_max = np.max(abs(Cs_q))
mu_max = np.max(abs(Cs_mu))

C_roof = np.max([q_max, mu_max])

Cs_q /= C_roof
Cs_mu /= C_roof

# Three particle case in SSM
# # First vertex with quarks
# Cs_q = [[C_V_A(Q_q), C_A_A],                # Photon
#         [C_V_Z(Q_q, I_q), C_A_Z(Q_q, I_q)], # Z
#         [C_V_Z(Q_q, I_q), C_A_Z(Q_q, I_q)]] # Z prime

# # Second vertex with muons    
# Cs_mu = [[C_V_A(Q_mu), C_A_A],                    # Photon
#           [C_V_Z(Q_mu, I_mu), C_A_Z(Q_mu, I_mu)],  # Z
#           [C_V_Z(Q_mu, I_mu), C_A_Z(Q_mu, I_mu)]]  # Z prime



