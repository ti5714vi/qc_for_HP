# -*- coding: utf-8 -*-
"""
Created on Mon Dec 15 14:27:31 2025

    True diagram functions

@author: Erik Bashore
"""

from constant_inputs import *
from main import *
from quantum_gates import gamma_matrices, gamma5_matrix, beta_matrix
from Drell_Yan_circuit import *


def sign(mu):
  if mu == 0:
    return -1
  else:
    return 1


def get_photon_diagram(s, theta = np.pi/3):

  prefactor = 1j*e**2*Q_q*Q_mu/s

  inputs = generate_kinematic_inputs(s, theta = theta)
  prod = lambda gamma, mu: sign(mu) \
                          * (inputs['spinor2'].T @ beta_matrix @ gamma @ inputs['spinor1']) \
                          * (inputs['spinor3'].T @ beta_matrix  @ gamma @ inputs['spinor4'])

  # Sum over all spinor products and mutliply with front-factor
  result = prefactor * ft.reduce(lambda x, y: x+y, map(prod, gamma_matrices, range(4)))

  return result


def get_Z_diagram(s, theta = np.pi/3, M_Z = M_Z):
  E = np.sqrt(s)/2
  spinor1 = get_spinor(E, s = '+', type = 'F', direction = 'in') # quark
  spinor2 = get_spinor(E, s = '-', type = 'A', direction = 'in') # anti-quark
  spinor3 = get_spinor(E, s = '+', type = 'F', direction = 'out', theta = theta) # muon-
  spinor4 = get_spinor(E, s = '-', type = 'A', direction = 'out', theta = theta) # muon+

  prefactor = -1j/(s - M_Z**2)

  # First prod with quark factors
  prod_q = lambda gamma: spinor2.T @ beta_matrix @ (C_V_Z(Q_q, I_q) * gamma + C_A_Z(Q_q, I_q) * gamma @ gamma5_matrix) @ spinor1

  # Second prod with muon factors
  prod_mu = lambda gamma: spinor3.T @ beta_matrix @ (C_V_Z(Q_mu, I_mu) * gamma + C_A_Z(Q_mu, I_mu) * gamma @ gamma5_matrix) @ spinor4

  full_prod = lambda gamma, mu: sign(mu) * prod_q(gamma) * prod_mu(gamma)

  result = prefactor * ft.reduce(lambda x, y: x + y, map(full_prod, gamma_matrices, range(4)))

  return result

def get_interference(M1, M2):
    prod = M1*(M2.conj())
    return 2*prod.real


# True checks
def true(s, theta, n_particles = 2):
    if n_particles == 2:
        return abs(get_photon_diagram(s, theta) + get_Z_diagram(s, theta))**2
    elif n_particles == 3:
        return abs(get_photon_diagram(s, theta) + get_Z_diagram(s, theta) + get_Z_diagram(s, theta, M_Z = M_Z_prime))**2
    else:
        print("Not valid number of particles")


def true_int(s, theta, n_particles = 2):
    if n_particles == 2:
        return get_interference(get_photon_diagram(s, theta), get_Z_diagram(s, theta))
    elif n_particles == 3:
        return get_interference(get_photon_diagram(s, theta) + get_Z_diagram(s, theta), get_Z_diagram(s, theta , M_Z = M_Z_prime))
    else:
        print("Not valid number of particles")


def true_cross_section(s):
    v_q, a_q = C_V_Z(Q_q, I_q), C_A_Z(Q_q, I_q)
    v_mu, a_mu = C_V_Z(Q_mu, I_mu), C_A_Z(Q_mu, I_mu)
    
    chi_Z = s/(s - M_Z**2) / (4* np.sin(theta_W)**2 * np.cos(theta_W)**2)
    
    sigma = 1/(3*s) * (Q_q**2 + 2 * Q_q * v_q * v_mu * chi_Z.real + (v_q**2 + a_q**2)**2 * (v_mu**2 + a_mu**2)**2 * abs(chi_Z)**2)
    
    return sigma.real

    








