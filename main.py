# -*- coding: utf-8 -*-
"""
Created on Fri Dec 12 13:07:07 2025

    main file with imports and functions

@author: Erik Bashore
"""

import numpy as np
from qiskit import *
from qiskit.quantum_info import Statevector, Operator
from qiskit.circuit.library import QFTGate, UnitaryGate
import functools as ft
import matplotlib. pyplot as plt
from qiskit_aer import Aer, AerSimulator
from qiskit.transpiler.passes import RemoveBarriers

black = {'backgroundcolor': '#171717',
         'linecolor': '#EEEEEE',
         'textcolor': '#EEEEEE'}


#%% 

# Functions

def binary(x, n):
  "Returns binary representation of x with length n"
  if n != 0:
      return '{0:b}'.format(x).zfill(n)
  else: 
      return ''

def register(regs):
  return sum([list(r) for r in regs])

def get_needed_qubits(n):    
    if n == 1:
        return 1
    else:
        return int(np.ceil(np.log2(n)))

def register(regs):
  qubits = []
  for r in regs:
    qubits += list(r)
  return qubits


# Pole functions
def single(k, M):
  "Takes in four-vector k, mass M and returns m_0"
  return -1j/(k@k - M**2)

def alt_single(s, M):
  "Takes in Mandelstam s, mass M and returns m_0"
  return -1j/(s - M**2)


def double(k, M, xi):
  "Takes in four-vector k, mass M, gauge fix label xi and returns m_1."
  "Feynman gauge: xi = 1, Landau gauge: xi = 0"
  return 1j * (1 - xi) / ((k@k - M**2)*(np.linalg.norm(k)**2 - xi*M**2))




def run(circ, shots = 5000, measure = 'all'):
  # Runs the given circuit for given number of times and returns
  # the counts as dictionary and a total runtime

  # Connect the quantum circuit to a measurement output
  if measure == 'all':
    circ.measure_all()
  else:
    circ.add_register(ClassicalRegister(len(measure)))
    circ.measure([int(s) for s in measure],  # Quantum bits
                 [int(s) for s in range(len(measure))])  # Classical bits

  # Setup backend
  aersim = AerSimulator()

  # Transpile the circuit
  trans_circ = transpile(circ, aersim)

  # Run the circuit
  job = aersim.run(trans_circ, shots = shots)
  counts = job.result().get_counts()
  runtime = job.result().time_taken

  circ.remove_final_measurements()
  print('Runtime:', f'{int(runtime/60)}m {round(runtime%60, 2)}s')
  return counts



# Spinor function

def get_spinor(E, s, type, direction, theta = np.pi/3):
    """
    Parameters
    ----------
    E : Energy - float
    s : Mandelstam variable - float
    type : Type of particle: Fermion/anti-fermion - F/A
    direction : Incoming or outgoing spinor - in/out
    theta : Scattering angle - float [0, 2pi), default is pi/3

    Returns
    -------
    Four-component massless spinor as an array
    
    """
    ct = np.cos(theta/2)
    st = np.sin(theta/2)
    
    if direction == 'in':
      if type == 'F':
        if s == '+':
          return np.sqrt(2*E) * np.array([1, 0, 1, 0])
        else:
          return np.sqrt(2*E) * np.array([0, 1, 0, -1])
      else:
        if s == '+':
          return np.sqrt(2*E) * np.array([0, 1, 0, 1])
        else:
          return np.sqrt(2*E) * np.array([1, 0, -1, 0])
    else:
      if type == 'F':
        if s == '+':
          return np.sqrt(2*E) * np.array([ct, st, ct, st])
        else:
          return np.sqrt(2*E) * np.array([-st, ct, st, -ct])
      else:
        if s == '+':
          return np.sqrt(2*E) * np.array([st, ct, -st, -ct])
        else:
          return np.sqrt(2*E) * np.array([ct, -st, -ct, st])



def SDandAveComp(data):
  # Returns the standard deviation of an array of data points
  N = len(data)
  ave = np.mean(data)

  square_diff = lambda x: (x - ave)**2
  sum = ft.reduce(lambda x, y: x + y, map(square_diff, data))
  SD = np.sqrt(sum / N)
  return SD, ave


def outlier_mask(data):
    for r, row in enumerate(data):
        for e, element in enumerate(row):
            reduced_data = data[data != element]
            
            diam = np.max(reduced_data) - np.min(reduced_data)
            
            dist = np.min(abs(reduced_data - element))
            
            if dist > diam and element != 0:
                print('Outlier found: ', element, '%')
                
                data[r, e] = 0
            
            if element == 100:
                data[r, e] = 0
                
    mask = data == 0
    
    return np.ma.masked_array(data, mask)

                
def integrate(result):
    output = np.zeros(result['Data']['PS points']) 
    
    # Getting the error matrix
    error_matrix = result['Error matrix'] 
    mask = np.ma.getmask(error_matrix)  # Masked array has true if masked, we want opposite
    
    # Removing the outlier points before computing
    result_ave = np.ma.masked_array(result['Average'], mask = mask).filled(0)
    
    # Theta range
    n_angle_points = result['Data']['Angle points']    
    t_range = np.linspace(0, np.pi, n_angle_points + 1)[:-1]
    dtheta = np.pi/n_angle_points
    
    for r, row in enumerate(result_ave):
        output[r] = 2 * np.pi * dtheta * sum([np.sin(t_range[k]) * row[k] for k in range(n_angle_points)])
        
    return output












