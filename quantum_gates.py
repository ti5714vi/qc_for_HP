# -*- coding: utf-8 -*-
"""
Created on Sat Dec 13 10:50:56 2025

    File with quantum gates    

@author: Erik Bashore
"""

from main import *
from constant_inputs import *


#%% 


# Useful gates

circ = QuantumCircuit(1)
circ.y(0)
X = circ.to_gate(label = 'X')

circ = QuantumCircuit(1)
circ.y(0)
Y = circ.to_gate(label = 'Y')

circ = QuantumCircuit(1)
circ.y(0)
Z = circ.to_gate(label = 'Z')


# Three binary state gate
T = UnitaryGate([[np.sqrt(1/3), np.sqrt(1/2) , np.sqrt(1/6), 0],
                 [np.sqrt(1/3), -np.sqrt(1/2), np.sqrt(1/6), 0],
                 [np.sqrt(1/3), 0, -np.sqrt(2/3), 0],
                 [0, 0, 0, 1]], label = 'T')

# Increment gate
def Inc_gate(n):
  circ = QuantumCircuit(n)
  if n == 1:
    circ.x(0)
  else:
    circ.x(0)
    for i in range(1, n):
      circ.append(X.control(i), [j for j in range(0, i+1)])
  return circ.inverse().to_gate(label = '+')


# Value-setting gate
def B_gate(alpha, n):
  circ = QuantumCircuit(n)

  # Compute corresponding a
  theta = 2 * np.arcsin(abs(alpha))  # rotation angle
  phi = np.angle(alpha)             # relative phase

  # Create single-qubit gate
  matrix = [
      [-np.cos(theta/2), np.exp(1j*phi)*np.sin(theta/2)],
      [-np.exp(-1j*phi)*np.sin(theta/2), -np.cos(theta/2)]
      ]

  b = UnitaryGate(matrix)

  # Create controlled gate with '000...0' as control state
  B = b.control(n-1, ctrl_state = '0'*(n-1))
  circ.append(B, [j for j in range(0,n)][::-1])
  return circ.to_gate(label = f'B({alpha})')


# Minkowski signature gate -+++
def metric_gate():
  circ = QuantumCircuit(2)
  circ.x(0)
  circ.cz(1, 0, ctrl_state = '0')
  circ.x(0)

  return circ.to_gate(label = 'M')



def eta_gate():
    circ = QuantumCircuit(4)
    
    circ.h([2, 3])
    circ.cx([0,1], [2, 3])
    circ.x(2)
    circ.cz(3, 2, ctrl_state = '0')
    circ.x(2)
    
    return circ.to_gate(label = r'$\eta$')

eta_gate = eta_gate()




#%% 

# Dirac matrices

gamma0_matrix = np.array([[1,0,0,0],
                          [0,1,0,0],
                          [0,0,-1,0],
                          [0,0,0,-1]])

gamma1_matrix = np.array([[0,0,0,1],
                          [0,0,1,0],
                          [0,-1,0,0],
                          [-1,0,0,0]])

gamma2_matrix = np.array([[0,0,0,-1j],
                          [0,0,1j,0],
                          [0,1j,0,0],
                          [-1j,0,0,0]])

gamma3_matrix = np.array([[0,0,1,0],
                          [0,0,0,-1],
                          [-1,0,0,0],
                          [0,1,0,0]])

gamma5_matrix = np.array([[0,0,1,0],
                          [0,0,0,1],
                          [1,0,0,0],
                          [0,1,0,0]])

beta_matrix = gamma5_matrix

gamma_matrices = [gamma0_matrix, gamma1_matrix, gamma2_matrix, gamma3_matrix]


#%% 

# Dirac quantum gates

# gamma0
circ = QuantumCircuit(2)
circ.z(1)

gamma0 = circ.to_gate(label = r'$\gamma^0$')

# gamma1
circ = QuantumCircuit(2)
circ.x(1)
circ.z(1)
circ.x(0)
gamma1 = circ.to_gate(label = r'$\gamma^1$')

# gamma2
circ = QuantumCircuit(2)
circ.x(1)
circ.z(1)
circ.y(0)
gamma2 = circ.to_gate(label = r'$\gamma^2$')

# gamma3
circ = QuantumCircuit(2)
circ.x(1)
circ.z(1)
circ.z(0)
gamma3 = circ.to_gate(label = r'$\gamma^3$')

gammas = [gamma0, gamma1, gamma2, gamma3]


# gamma5
circ = QuantumCircuit(2)
circ.x(1)
gamma5 = circ.to_gate(label = r'$\gamma^5$')



#%%

# Single pole gate

def S_gate(n_particles, ms, n):
  assert len(ms) == n_particles, 'Number of poles not equal to number of particles'

  i1 = QuantumRegister(2, name = 'i_1')
  i2 = QuantumRegister(2, name = 'i_2')
  U = QuantumRegister(n, name = 'U')
  p = QuantumRegister(get_needed_qubits(n_particles), name = 'p')
  circ = QuantumCircuit(i1, i2, U, p)

  # Apply metric tensor
  circ.append(eta_gate, list(i1) + list(i2))

  # Apply corresponding poles for all relevant particles
  for i, m in enumerate(ms):
    circ.append(B_gate(alpha = m/2, n = n).control(len(p), ctrl_state = binary(i, len(p))), register([p, U]))
  return circ.to_gate(label = r'$\mathcal{S}$')


# Double pole gate

def D_gate(n_particles, ms, k, n):
  assert len(ms) == n_particles, 'Number of poles not equal to number of particles'

  i1 = QuantumRegister(2, name = 'i_1')
  i2 = QuantumRegister(2, name = 'i_2')
  U = QuantumRegister(n, name = 'U')
  p = QuantumRegister(get_needed_qubits(n_particles), name = 'p')
  circ = QuantumCircuit(i1, i2, U, p)

  for i, m in enumerate(ms):
    circ.append(B_gate(alpha = m, n = n).control(len(p), ctrl_state = binary(i, len(p))), register([p, U]))

  for i in [i1, i2]:
    circ.append(Inc_gate(n = n), U)
    for mu in range(4):
      circ.append(B_gate(alpha = k[mu], n = n).control(2, ctrl_state = binary(mu, 2)), register([i, U]))

  circ.append(metric_gate(), i1)
  circ.append(metric_gate(), i2)
  return circ.to_gate(label = r'$\mathcal{D}$')


# Propagator gate

def P_gate(n_particles, ms, k, n):
  # ms = [[m0, ...., m0],
  #       [m1, ...., m1]]
  m0s, m1s = ms[0], ms[1]

  # Define registers
  i1 = QuantumRegister(2, name = 'i_1')
  i2 = QuantumRegister(2, name = 'i_2')
  U = QuantumRegister(n, name = 'U')
  p = QuantumRegister(get_needed_qubits(n_particles), name = 'p')
  a = QuantumRegister(1, name = 'a')

  # Assemble circuit
  circ = QuantumCircuit(i1, i2, U, p, a)

  # Add initial increment and ancilla superposition
  circ.append(Inc_gate(n), U)
  circ.h(a)

  # Add single and double pole gates
  circ.append(S_gate(n_particles = n_particles, ms = m0s, n = n).control(1, ctrl_state = '0'), register([a, i1, i2, U, p]))
  circ.append(D_gate(n_particles = n_particles, ms = m1s, k = k,  n = n).control(1, ctrl_state = '1'), list(a) + list(i1) + list(i2) + list(U) + list(p))

  # Close ancilla
  circ.h(a)

  return circ.to_gate(label = r'$\mathcal{P}$')



#%%

# Alternate ctrl-propagator gate

def alt_S_gate(n_particles, n_u, n_PS, s_min = 80, s_max = 100):

  s_range = np.linspace(s_min, s_max, n_PS)

  PS = QuantumRegister(get_needed_qubits(n_PS), name = 'PS')
  i1 = QuantumRegister(2, name = 'i_1')
  i2 = QuantumRegister(2, name = 'i_2')
  U = QuantumRegister(n_u, name = 'U')
  p = QuantumRegister(get_needed_qubits(n_particles), name = 'p')
  circ = QuantumCircuit(PS, i1, i2, U, p)

  # Apply metric tensor
  circ.append(eta_gate, list(i1) + list(i2))

  # Apply corresponding poles for all relevant particles
  s_labels = [s_label for s_label in range(0, n_PS)]
  for s_label in s_labels:
    s = s_range[s_label]**2
    
    # Pole factors
    ms = np.array([alt_single(s, M) for M in [M_A, M_Z]])
    m_roof = max([abs(m) for m in ms])
    
    ms /= m_roof
    
    for i, m in enumerate(ms):
      circ.append(B_gate(alpha = m/2, n = n_u).control(len(p) + len(PS),
                                                       ctrl_state = binary(i, len(p)) + binary(s_label, len(PS))),
                  register([PS, p, U]))

    # poles = temp_circ.to_gate(label = f'poles({s_label})')

    # circ.append(poles.control(len(PS), ctrl_state = binary(s_label, len(PS))), register([PS, U, p]))
  return circ.to_gate(label = r'$\mathcal{S}$')


def ctrl_P_gate(n_particles, n_u, n_PS, s_min = 80, s_max = 100):

  s_range = np.linspace(s_min, s_max, n_PS)

  # Define registers
  PS = QuantumRegister(get_needed_qubits(n_PS), name = 'PS')
  i1 = QuantumRegister(2, name = 'i_1')
  i2 = QuantumRegister(2, name = 'i_2')
  U = QuantumRegister(n_u, name = 'U')
  p = QuantumRegister(get_needed_qubits(n_particles), name = 'p')
  a = QuantumRegister(1, name = 'a')

  # Assemble circuit
  circ = QuantumCircuit(PS, i1, i2, U, p, a)

  # Add initial increment and ancilla superposition
  circ.append(Inc_gate(n_u), U)
  # circ.h(a)

  # Add single and double pole gates
  circ.append(S_gate(n_particles = n_particles,
                   n_u = n_u, n_PS = n_PS).control(1, ctrl_state = '0'), register([a, PS, i1, i2, U, p]))

  # Close ancilla
  # circ.h(a)

  return circ.to_gate(label = r'$\mathcal{P}$')



#%% 

# Vertex gate

def V_gate(n_particles, Cs, n):
  # Cs = [[C_V1, C_A1], [C_V2, C_A2], ..., [C_Vn, C_An]]
  assert len(Cs) == n_particles, 'Number of couplings (' + str(len(Cs)) + ') not equal to number of particles (' + str(n_particles) + ')'
  i = QuantumRegister(2, name = 'i')
  t = QuantumRegister(2, name = 't')
  U = QuantumRegister(n, name = 'U')
  a = QuantumRegister(1, name = 'a')
  p = QuantumRegister(get_needed_qubits(n_particles), name = 'p')
  circ = QuantumCircuit(t, i, U, p, a)

  circ.h(a)
  circ.append(Inc_gate(n), U)

  for index, pair in enumerate(Cs):
    C_V, C_A = pair[0], pair[1]
    circ.append(B_gate(C_V, n).control(len(p) + 1, ctrl_state = '0' + binary(index, len(p))), register([p, a, U]))
    circ.append(B_gate(C_A, n).control(len(p) + 1, ctrl_state = '1' + binary(index, len(p))), register([p, a, U]))

  circ.append(gamma5.control(1, ctrl_state = '1'), register([a, i]))

  for mu in range(4):
    circ.append(gammas[mu].control(2, ctrl_state = binary(mu, 2)), register([t, i]))

  circ.h(a)

  return circ.to_gate(label = r'$\mathcal{V}$')




#%% 

# Beta gate

circ = QuantumCircuit(2)

circ.x(1)
beta_gate = circ.to_gate(label = r'$\beta$')


#%% 

# Spinor gates

def SpinorGate(spinor, label):
  spinor_hat = spinor/np.linalg.norm(spinor)

  alpha = np.sqrt(spinor_hat[0]**2 + spinor_hat[1]**2)
  beta = np.sqrt(spinor_hat[2]**2 + spinor_hat[3]**2)

  theta0 = 2*np.arccos(spinor_hat[0]/alpha)
  theta1 = 2*np.arccos(spinor_hat[2]/beta)

  circ = QuantumCircuit(2)

  circ.ry(2*np.arccos(alpha), 1)

  circ.cry(theta0, 1, 0, ctrl_state = '0')
  if np.sign(spinor[1]) == -1:
    circ.cz(1, 0, ctrl_state = '0')

  circ.cry(theta1, 1, 0, ctrl_state = '1')
  if np.sign(spinor[3]) == -1:
    circ.cz(1, 0, ctrl_state = '1')

  U = circ.to_gate(label = fr'$U_{label}$')
  return U


def BarredSpinorGate(spinor, label):
  circ = QuantumCircuit(2)

  # Create barred version
  circ = QuantumCircuit(2)
  circ.append(SpinorGate(spinor.conj(), label = label), [0,1])
  circ.append(beta_gate, [0,1])

  U_bar = circ.inverse().to_gate(label = fr'$U_{label}^\dagger$')
  return U_bar


#%% 

def multi_spinor_gate(n_angle_points, energy_dependent = False, E = 0, 
                      spin = '+', type = 'F', direction = 'out'):
    
    theta_range = np.linspace(0, np.pi, n_angle_points + 1)[:-1]
    
    Theta = QuantumRegister(get_needed_qubits(n_angle_points), name = 'T')
    v = QuantumRegister(2, name = 'v')
    
    circ = QuantumCircuit(Theta, v)
    
    
    for t, theta in enumerate(theta_range):
        if energy_dependent == False:
            E = 1/2
        # else: 
            # E = np.sqrt(s)/2
        spinor = get_spinor(E = E, s = spin, type = type, direction = direction, 
                            theta = theta_range[t])
        U = SpinorGate(spinor = spinor,
                       label = str(t))
        
        circ.append(U.control(len(Theta), ctrl_state = binary(t, len(Theta))), register([Theta, v]))
        
    print("\nTheta discretization points: ", end = '')
    for theta in theta_range:
        print(theta, end = ' ')
    print('\n')    
    return circ.to_gate(label = 'MU')












