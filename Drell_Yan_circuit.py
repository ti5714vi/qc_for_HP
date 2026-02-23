# -*- coding: utf-8 -*-
"""
Created on Sat Dec 13 11:15:54 2025

    Drell-Yan scattering input and circuit generators

@author: Erik Bashore
"""

from main import *
from quantum_gates import *
from constant_inputs import *


def get_m_roof(s, n_particles = 2):
    if n_particles == 2:
        masses = [M_A, M_Z]
    if n_particles == 3:
        masses = [M_A, M_Z, M_Z_prime]
    
    ms = [alt_single(s, M) for M in masses]

    m_roof = max([abs(m) for m in ms])    
    return m_roof


def generate_kinematic_inputs(s = 1000,  theta = np.pi/3, n_particles = 2):
  # Put in energy
  E = np.sqrt(s)/2

  # Four-vectors
  p1 = np.array([E, 0, 0, E])
  p2 = np.array([E, 0, 0, -E])
  p3 = np.array([E, E*np.sin(theta), 0, E*np.cos(theta)])
  p4 = np.array([E, -E*np.sin(theta), 0, -E*np.cos(theta)])

  # Internal momenta
  k = p1 + p2
  # Normalize
  k_hat = k/np.linalg.norm(k)


  # Get associated spinors
  spinor1 = get_spinor(E, s = '+', type = 'F', direction = 'in') # quark
  spinor2 = get_spinor(E, s = '-', type = 'A', direction = 'in') # anti-quark
  spinor3 = get_spinor(E, s = '+', type = 'F', direction = 'out', theta = theta) # muon-
  spinor4 = get_spinor(E, s = '-', type = 'A', direction = 'out', theta = theta) # muon+


  # Colllect them
  spinors = [spinor1, spinor2, spinor3, spinor4]

  # Gauge parameter
  xi = 1

  # Pole terms
  if n_particles == 2:
      ms = np.array([[single(k, M_A), single(k, M_Z)],           # single poles
                    [double(k, M_A, xi), double(k, M_Z, xi)]])   # double poles
      
  elif n_particles == 3:
      ms = np.array([[single(k, M_A), single(k, M_Z), single(k, M_Z_prime)],              # single poles
                    [double(k, M_A, xi), double(k, M_Z, xi), double(k, M_Z_prime, xi)]])  # double poles

  # Stupid way of finding largest m
  func = lambda i: list(map(lambda x: abs(x), ms[i]))
  m_roof = max(ft.reduce(lambda x, y: x + y, map(func, [0,1])))
  
  # Normalize
  ms /= m_roof
      
  return {'spinor1': spinor1, 'spinor2': spinor2, 'spinor3': spinor3, 'spinor4': spinor4,
          'E': E, 'k': k, 'k_hat': k_hat, 'ms': ms, 'm_roof': m_roof, 'theta': theta}


def compensate_calc(s_range, t_range):
    
    s_compensate = np.zeros(len(s_range))
    theta_compensate = np.zeros(len(t_range))
    
    m_roofs = np.array([get_m_roof(s_sqrt**2) for s_sqrt in s_range])
    
    for l, s_sqrt in enumerate(s_range):
        s_compensate[l] = m_roofs[l] * s_sqrt**2 /4
    
    # for t, theta in enumerate(t_range):
    #     spinor3 = get_spinor(E = 1/2, s = '+', type = 'F', direction = 'out', theta = theta)
    #     spinor4 = get_spinor(E = 1/2, s = '-', type = 'A', direction = 'out', theta = theta)
        
    #     theta_compensate[t] = np.linalg.norm(spinor3)**2 * np.linalg.norm(spinor4)**2  
    
    # compensate = np.array([[s_compensate[s] * theta_compensate[t] for s in range(len(s_range))] for t in range(len(t_range))])
    
    compensate = 2*s_compensate
    return compensate



def generate_circuit(inputs, n_particles = 2, 
                     event = 'Interference', 
                     hit_scheme = True, 
                     PS_integration = False, n_PS = 1, s_min = 80, s_max = 100,
                     angle_integration = False, n_angle_points = 1):

  # Size of unitarity register
  n_u = 3

  # Define registers
  registers = []
  
  if PS_integration == True:
      PS = QuantumRegister(get_needed_qubits(n_PS), name = 'PS')
      registers.append(PS)
      
      s_range = np.linspace(s_min, s_max, n_PS)
      # print(f'Phase space integration activated in s range [{s_min}, {s_max}] with {n_PS} point discretization')
      
  else: 
      s_range = np.array([4*inputs['E']**2])
      
  if angle_integration == True:
      Theta = QuantumRegister(get_needed_qubits(n_angle_points), name = 'T')
      registers.append(Theta)
      
      theta_range = np.linspace(0, np.pi, n_angle_points + 1)[:-1]
      
  else:
      theta_range = np.array([inputs['theta']])
      
      
  if hit_scheme == True:
      hit = QuantumRegister(1, name = 'hit')
      registers.append(hit)
  
  v1 = QuantumRegister(2, name = 'v1')
  v2 = QuantumRegister(2, name = 'v2')
  i1 = QuantumRegister(2, name = 'i1')
  i2 = QuantumRegister(2, name = 'i2')
  U = QuantumRegister(n_u, name = 'U')
  p = QuantumRegister(get_needed_qubits(n_particles), name = 'p')
  a1 = QuantumRegister(1, name = 'a1')
  a2 = QuantumRegister(1, name = 'a2')
  a3 = QuantumRegister(1, name = 'a3')
  
      
  # Assemble circuit    
  registers += [v1, v2, i1, i2, U, p, a1, a2, a3]
      
  circ = QuantumCircuit(*registers)
  
  # Find spinors
  if PS_integration == False:
      spinor1 = inputs['spinor1']
      spinor2 = inputs['spinor2']
      spinor3 = inputs['spinor3']
      spinor4 = inputs['spinor4']
      
  else:
      # Energy independent spinors
      spinor1 = get_spinor(E = 1/2, s = '+', type = 'F', direction = 'in')
      spinor2 = get_spinor(E = 1/2, s = '-', type = 'A', direction = 'in')
      spinor3 = get_spinor(E = 1/2, s = '+', type = 'F', direction = 'out')
      spinor4 = get_spinor(E = 1/2, s = '-', type = 'A', direction = 'out')

  # Create spinor gates
  U1, U1_bar = SpinorGate(spinor1, label = '1'), BarredSpinorGate(spinor1, label = '1')
  U2, U2_bar = SpinorGate(spinor2, label = '2'), BarredSpinorGate(spinor2, label = '2')
  U3, U3_bar = SpinorGate(spinor3, label = '3'), BarredSpinorGate(spinor3, label = '3')
  U4, U4_bar = SpinorGate(spinor4, label = '4'), BarredSpinorGate(spinor4, label = '4')

  
  # Check how many particles
  if n_particles == 2:
      # Initialize superposition of photon and Z
      if event == 'Interference':
          circ.h(p)
      elif event == 'Z':
          circ.x(p)
      else: 
          print('Not valid event type')
          
  elif n_particles == 3:
      circ.append(T, p) # Triple state gate 

  # # Initialize indices
  circ.h(i1)
  circ.h(i2)

  if PS_integration == True:
      circ.h(PS)
      
  if angle_integration == True:
      circ.h(Theta)  
      
  circ.barrier()

  # Apply spinor for incoming fermion
  circ.append(U1, v1)

  # Apply first vertex with quarks
  circ.append(V_gate(n_particles = n_particles, Cs = Cs_q, n = n_u), 
              register([i1, v1, U, p, a1]))

  # Apply barred spinor for incoming anti-fermion
  circ.append(U2_bar, v1)

  circ.barrier()

  if PS_integration == False:
      # Apply internal propagator
      circ.append(P_gate(n_particles = n_particles, ms = inputs['ms'], k = inputs['k_hat'], n = n_u), 
                  register([i1, i2, U, p, a2]))

  else:
    # Apply single pole gate controlled by the PS register
    circ.append(Inc_gate(n_u), U)
    circ.append(alt_S_gate(n_particles = n_particles, n_u = n_u, n_PS = n_PS,
                           s_min = s_min, s_max = s_max), 
                register([PS, i1, i2, U, p]))
      
  circ.barrier()
              
  # Apply spinor for outgoing anti-fermion 
  if angle_integration == True:
      
      # Use energy independent spinor gates if PS scheme is activated
      if PS_integration == True:
          E = 1/2

      # Otherwise input s    
      elif PS_integration == False:
          E = inputs['E']
          
      for t, theta in enumerate(theta_range):
          spinor = get_spinor(E = inputs['E'], s = '-', type = 'A', direction = 'outgoing', 
                              theta = theta_range[t])
          gate = SpinorGate(spinor = spinor,
                                  label = str(t))
          
          circ.append(gate.control(len(Theta), ctrl_state = binary(t, len(Theta))), 
                      register([Theta, v2]))
          
  else:
      circ.append(U4, v2)
  
  # Apply second vertex with muons
  circ.append(V_gate(n_particles = n_particles, Cs = Cs_mu, n = n_u), 
              register([i2, v2, U, p, a3]))

  # Apply barred spinor for outgoing fermion
  if angle_integration == True:
      
      # Use energy independent spinor gates if PS scheme is activated
      if PS_integration == True:
          E = 1/2

      # Otherwise input s    
      elif PS_integration == False:
          E = inputs['E']
          
          
      for t, theta in enumerate(theta_range):
          spinor = get_spinor(E = inputs['E'], s = '+', type = 'F', direction = 'outgoing', 
                              theta = theta_range[t])
          gate = BarredSpinorGate(spinor = spinor,
                                  label = str(t))
          
          circ.append(gate.control(len(Theta), ctrl_state = binary(t, len(Theta))), 
                      register([Theta, v2]))
  else:
      circ.append(U3_bar, v2)

  circ.barrier()

  # Final Hadamard
  circ.h(register([i1, i2]))

  # Close particle register
  if event == 'Interference':
      circ.h(p)
  else:
      circ.x(p)
      
  if hit_scheme == True:
      blank_qubits = circ.num_qubits - 1 - len(p) # Number of qubits to project to vacuum
      
      if PS_integration == True:
          blank_qubits -= len(PS)
          
      if angle_integration == True:
          blank_qubits -= len(Theta)
            
      circ.append(X.control(blank_qubits, ctrl_state = '0'*blank_qubits), 
                  register([v1, v2, i1, i2, U, a1, a2, a3, hit]))
      
  # Compensation factor
  compensate = 256 * np.linalg.norm(spinor1) * np.linalg.norm(spinor2) * C_roof \
                   * compensate_calc(s_range, theta_range)
  
  # Extra gate factors  
  if PS_integration == True:              
      compensate *= np.sqrt(2**len(PS)) * 2 # Extra PS Hadamards and the skipping secondary ancilla 
      # compensate /= inputs['m_roof']
      
  if angle_integration == True:
      compensate *= np.sqrt(2**len(Theta)) # Hadamard compensation
      
  if n_particles == 3:
      compensate *= np.sqrt(3) # Hadamard compensation
      
  if event == 'Z':
      compensate /= 2  # Adjust for the missing Hadamard gate on p
      
  return circ, compensate



