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
from itertools import product

def sign(mu):
  if mu == 0:
    return -1
  else:
    return 1


def get_photon_diagram(s, theta = np.pi/3, hel=['+','-','-','+']):

  prefactor = e**2*Q_mu*Q_mu/s

  inputs = generate_kinematic_inputs(s, theta = theta, hel=hel)
  prod = lambda gamma, mu: sign(mu) \
                          * (inputs['spinor2'].T @ beta_matrix @ gamma @ inputs['spinor1']) \
                          * (inputs['spinor3'].T @ beta_matrix  @ gamma @ inputs['spinor4'])

  # Sum over all spinor products and mutliply with front-factor
  result = prefactor * ft.reduce(lambda x, y: x+y, map(prod, gamma_matrices, range(4)))

  return result


def get_Z_diagram(s, theta = np.pi/3, M_Z = M_Z, hel=['+','-','-','+']):
  E = np.sqrt(s)/2
  spinor1 = get_spinor(E, s = hel[0], type = 'F', direction = 'in') # quark
  spinor2 = get_spinor(E, s = hel[1], type = 'A', direction = 'in') # anti-quark
  spinor3 = get_spinor(E, s = hel[2], type = 'F', direction = 'out', theta = theta) # muon-
  spinor4 = get_spinor(E, s = hel[3], type = 'A', direction = 'out', theta = theta) # muon+

  #print('spinors')
  #print(spinor1)
  #print(spinor2)
  #print(spinor3)
  #print(spinor4)
  prefactor = -1/(s - M_Z**2 + 1j*M_Z*Gamma_Z)

  # First prod with quark factors
  prod_el = lambda gamma: spinor2.T @ beta_matrix @ gamma @ (C_L_Z(Q_e, I_e)   * P_L + C_R_Z(Q_e, I_e) * P_R) @ spinor1

  # Second prod with muon factors
  prod_mu = lambda gamma: spinor3.T @ beta_matrix @ gamma @ (C_L_Z(Q_mu, I_mu) * P_L + C_R_Z(Q_mu, I_mu) * P_R) @ spinor4

  full_prod = lambda gamma, mu: sign(mu) * prod_el(gamma) * prod_mu(gamma)

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

   

# define your points
y_vals = np.linspace(0, 3.14159, 16)
x_vals = [85,87.5,90,92.5]

with open("xyz_output.dat", "w") as f:
    for x in x_vals:
        for y in y_vals:
            value = true(x**2, y)
            f.write(f"{y} {x} {value}\n")


##################################################
#    checks with MadGraph output
##################################################

energy=100.0
theta = 2.0
s = 4*energy**2

print('Photon diagram with energy E=',energy,' GeV, theta=',theta,': ')
photon=get_photon_diagram(s,theta)
print(photon*np.conj(photon))
print('Z-boson diagram with energy E=',energy,' GeV, theta=',theta,': ')
z_diag=get_Z_diagram(s,theta)
print(z_diag*np.conj(z_diag))

print('Both diagrams with energy E=',energy,' GeV, theta=',theta,': ')
tot = true(s,theta)
print(tot)

# Define the symbols
symbols = ['+', '-']

# Generate all 4-length combinations
all_combinations = list(product(symbols, repeat=4))

# Convert each tuple to a list
all_combinations = [list(t) for t in all_combinations]

# Check
summ_Z=0.0
summ_A=0.0
for combo in all_combinations:
    z_diag=get_Z_diagram(s, theta, hel=combo)
    summ_Z = summ_Z + z_diag*np.conj(z_diag)
    photon = get_photon_diagram(s, theta, hel=combo)
    summ_A = summ_A + photon*np.conj(photon)

print('Summed over all helcities for Z:',summ_Z)
print('Summed over all helcities for A:',summ_A)


