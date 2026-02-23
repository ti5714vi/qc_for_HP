# -*- coding: utf-8 -*-
"""
Created on Sat Dec 13 11:06:19 2025

        Primary file for simulations

@author: Erik Bashore
"""

from main import alt_single
from Drell_Yan_circuit import *
from true_diagrams import *

plt.rcParams.update({
    "text.usetex": True,        # Enable full LaTeX
    "font.family": "serif",     # Use serif font
    "text.latex.preamble": r"\usepackage{amsmath}"  # Optional LaTeX packages
})

circ = QuantumCircuit(1)
circ.x(0)
X = circ.to_gate(label = 'X')



hit = QuantumRegister(1, name = 'hit')
v1 = QuantumRegister(2, name = 't1')
v2 = QuantumRegister(2, name = 't2')
i1 = QuantumRegister(2, name = 'i_1')
i2 = QuantumRegister(2, name = 'i_2')
U = QuantumRegister(3, name = 'U')
p = QuantumRegister(get_needed_qubits(2), name = 'p')
a1 = QuantumRegister(1, name = 'a1')
a2 = QuantumRegister(1, name = 'a2')
a3 = QuantumRegister(1, name = 'a3')


#%%

# Z cross-section simulation

n_points = 10
n_batches = 10

s_min, s_max = 80, 100

s_rough = np.linspace(s_min, s_max, n_points)
s_fine = np.linspace(s_min, s_max, 1000)


circ_outs = np.zeros((n_points, n_batches))
circ_errs = []



for s_label, s_sqrt in enumerate(s_rough):
  s = s_sqrt**2
  print('s = ', s)

  # Create circuit
  inputs = generate_kinematic_inputs(s)
  circ, compensate = generate_circuit(inputs, event = 'Z')

  circ = RemoveBarriers()(circ)

  Unitary = circ.to_gate(label = r'$\mathcal{M}_Z$')

  circ = QuantumCircuit(hit, v1, v2, i1, i2, U, p, a1, a2, a3)

  circ.append(Unitary, register([v1, v2, i1, i2, U, p, a1, a2, a3]))

  circ.append(X.control(15, ctrl_state = '0'*15), register([v1, v2, i1, i2, U, p, a1, a2, a3, hit]))
  
  for k in range(n_batches):
    print('Batch: ', k + 1, end = ' ')
    shots = 1e7
    counts = run(circ, shots = shots, measure = [0])
    
    output = counts['1']/shots * (compensate)**2
    circ_outs[s_label, k] = output 

  SD, ave = SDandAveComp(circ_outs[s_label])

  print('Run complete')
  print('Output:', ave , '±', SD)
  print('True value:', abs(get_Z_diagram(s))**2)
  print()


#%%

# Plot the batches
for batch in circ_outs.T:
    plt.scatter(s_rough, batch, alpha = 0.2, color = 'gray', s = 30)

# Plot the average of the batches
plt.scatter(s_rough, [np.mean(row) for row in circ_outs], 
            color = 'purple', alpha = 0.75, s = 50, label = 'Circuit output average')

# Plot the true value
plt.plot(s_fine, [abs(get_Z_diagram(s_sqrt**2))**2 for s_sqrt in s_fine], 
         linestyle = '--', alpha = 0.3, color = 'black', label = 'True value')


plt.yscale('log')
plt.legend()

plt.xlabel(r'$\sqrt{s}$')
plt.ylabel(r'$|\mathcal{M}_Z|^2$')
plt.tight_layout()


#%% 

# gamma + Z cross section and interference

n_points = 7
n_batches = 10

s_min, s_max = 90, 92

s_rough = np.linspace(s_min, s_max, n_points)
s_fine = np.linspace(s_min, s_max, 1000)



#%%

circ_outs = np.zeros((n_points, n_batches))      # For the full |gamma + Z|^2
circ_outs_int = np.zeros((n_points, n_batches))  # For the isolated interference

for s_label, s_sqrt in enumerate(s_rough):
    s = s_sqrt**2
    print('s = ', s)
    
    inputs = generate_kinematic_inputs(s)
    circ, compensate = generate_circuit(inputs, event = 'Interference')
    
    for k in range(n_batches):
        print('Batch: ', k + 1, end = ' ')
        shots = 3e7
        counts = run(circ, shots = shots, measure = [0, 12])
        
        if '01' in counts:
            Nplus = counts['01']
        else:
            Nplus = 0
        if '11' in counts:
            Nminus = counts['11']
        else:
            Nminus = 0
    
        # Convert counts to outputs
        amp_output = Nplus * compensate**2 / shots
        interference_output = (Nplus - Nminus)/2 * compensate**2 / shots
    
        circ_outs[s_label, k] = amp_output
        circ_outs_int[s_label, k] = interference_output
    
        
    print('Run complete')
    
    # Normal amp    
    SD, ave = SDandAveComp(circ_outs[s_label])
    print('\nFull term')
    print('-----------')
    print('Output:', ave , '±', SD)
    print('True value:', abs(get_photon_diagram(s) + get_Z_diagram(s))**2)
    print()
    
    # Interference amp    
    SD, ave = SDandAveComp(circ_outs_int[s_label])
    print('Interference')
    print('--------------')
    print('Output:', ave , '±', SD)
    print('True value:', get_interference(get_photon_diagram(s), get_Z_diagram(s)))
    print()


#%%
# Full term plot

# Plot the batches
for batch in circ_outs.T:
    plt.scatter(s_rough, batch, alpha = 0.2, color = 'gray', s = 30)

# Plot the average of the batches
plt.scatter(s_rough, [np.mean(row) for row in circ_outs], 
            color = 'purple', alpha = 0.75, s = 50, label = 'Circuit output average')

# Plot the true value
plt.plot(s_fine, [abs(get_photon_diagram(s_sqrt**2) + get_Z_diagram(s_sqrt**2))**2 for s_sqrt in s_fine], 
         linestyle = '--', alpha = 0.3, color = 'black', label = 'True value')

plt.vlines(x = M_Z.real, ymin = 0, ymax = 5000, color = 'purple', alpha = 0.2, label = r'$M_Z$')

plt.ylim(1000, 4000)
plt.yscale('log')
plt.legend()

plt.xlabel(r'$\sqrt{s}$')
plt.ylabel(r'$|\mathcal{M}_\gamma + \mathcal{M}_Z|^2$')
plt.tight_layout()

#%%
# Interference plot

# Plot the batches
for batch in circ_outs_int.T:
    plt.scatter(s_rough, batch, alpha = 0.2, color = 'gray', s = 30)

# Plot the average of the batches
plt.scatter(s_rough, [np.mean(row) for row in circ_outs_int], 
            color = 'purple', alpha = 0.75, s = 50, label = 'Circuit output average')


# Plot the true value
plt.plot(s_fine, 
         [get_interference(get_photon_diagram(s_sqrt**2), get_Z_diagram(s_sqrt**2)) for s_sqrt in s_fine], 
         linestyle = '--', alpha = 0.3, color = 'black', label = 'True value')

plt.vlines(x = M_Z.real, ymin = -1000, ymax = 1000, color = 'purple', alpha = 0.2, label = r'$M_Z$')

plt.ylim(-500, 700)
plt.legend()
plt.xlabel(r'$\sqrt{s}$')
plt.ylabel(r'Int$(\mathcal{M}_\gamma, \mathcal{M}_Z)$')
plt.tight_layout()


#%%

# Inclusion of Z prime boson

n_points = 10
n_batches = 1

s_min, s_max = 60, 1200

s_rough = np.linspace(s_min, s_max, n_points)
s_fine = np.linspace(s_min, s_max, 1000)


circ_outs = np.zeros((n_points, n_batches))      # For the full |gamma + Z + Z´|^2
circ_outs_int = np.zeros((n_points, n_batches))  # For the isolated interference

for s_label, s_sqrt in enumerate(s_rough):
    s = s_sqrt**2
    print('s = ', s)
    
    inputs = generate_kinematic_inputs(s, n_particles = 3)
    circ, compensate = generate_circuit(inputs, n_particles = 3)
    
    for k in range(n_batches):
        print('Batch: ', k + 1, end = ' ')
        shots = 2e7
        counts = run(circ, shots = shots, measure = [0, 12, 13]) # Measure hit and particle register
        
        if '001' in counts:  # Zero and hit
            Nplus = counts['001']
        else:
            Nplus = 0
        if '101' in counts:  # Two and hit
            Nminus = counts['101']
        else:
            Nminus = 0
    
        # Convert counts to outputs
        amp_output = Nplus * compensate**2 / shots
        interference_output = (Nplus - Nminus)/2 * compensate**2 / shots
    
        circ_outs[s_label, k] = amp_output
        circ_outs_int[s_label, k] = interference_output
    
    print('Run complete')
    
    # Normal amp    
    SD, ave = SDandAveComp(circ_outs[s_label])
    print('\nFull term')
    print('-----------')
    print('Output:', ave , '±', SD)
    print('True value:', 
          abs(get_photon_diagram(s) + get_Z_diagram(s) + get_Z_diagram(s, M_Z = M_Z_prime))**2)
    print()
    
    # Interference amp    
    SD, ave = SDandAveComp(circ_outs_int[s_label])
    print('Interference')
    print('--------------')
    print('Output:', ave , '±', SD)
    print('True value:', 
          get_interference(get_photon_diagram(s) + get_Z_diagram(s), get_Z_diagram(s, M_Z = M_Z_prime)))
    print()


#%%
# Full term plot

# Plot the batches
for batch in circ_outs.T:
    plt.scatter(s_rough, batch, alpha = 0.2, color = 'gray', s = 30)

# Plot the average of the batches
plt.scatter(s_rough, [np.mean(row) for row in circ_outs], 
            color = 'purple', alpha = 0.75, s = 50, label = 'Circuit output average')

# Plot the true value
true = lambda s: abs(get_photon_diagram(s) + get_Z_diagram(s) + get_Z_diagram(s, M_Z = M_Z_prime) )**2
plt.plot(s_fine, [true(s_sqrt**2) for s_sqrt in s_fine], 
         linestyle = '--', alpha = 0.3, color = 'black', label = 'True value')


plt.yscale('log')
plt.legend()

plt.xlabel(r'$\sqrt{s}$')
plt.ylabel(r'$|\mathcal{M}_\gamma + \mathcal{M}_Z + \mathcal{M}_{Z^\prime}|^2$')
plt.tight_layout()

#%%
# Interference plot

# Plot the batches
for batch in circ_outs_int.T:
    plt.scatter(s_rough, batch, alpha = 0.2, color = 'gray', s = 30)

# Plot the average of the batches
plt.scatter(s_rough, [np.mean(row) for row in circ_outs_int], 
            color = 'purple', alpha = 0.75, s = 50, label = 'Circuit output average')

# Plot the true value

true = lambda s: get_interference(get_photon_diagram(s) + get_Z_diagram(s), 
                                  get_Z_diagram(s, M_Z = M_Z_prime))
plt.plot(s_fine, 
         [true(s_sqrt**2) for s_sqrt in s_fine], 
         linestyle = '--', alpha = 0.3, color = 'black', label = 'True value')


plt.legend()
plt.xlabel(r'$\sqrt{s}$')
plt.ylabel(r'Int$(\mathcal{M}_\gamma + \mathcal{M}_Z, \mathcal{M}_{Z^\prime})$')
plt.tight_layout()

#%%

# PS scheme 
n = 3
n_PS = 2**n             # Number of discretization points
s_min, s_max = 90, 92   # Range of s
s_range = np.linspace(s_min, s_max, n_PS)


batches = 10
shots = 1e7 # Number of shots per batch

outputs = np.zeros((n_PS, batches))
outputs_int = np.zeros((n_PS, batches))
errors = np.zeros(n_PS)
errors_int = np.zeros(n_PS)

inputs = generate_kinematic_inputs(s = 1,            # Placeholder s
                                   theta = 3*np.pi/4,  # Scattering angle
                                   n_particles = 2) 

circ, compensate = generate_circuit(inputs, 
                                    event = 'Interference', 
                                    n_particles = 2,
                                    PS_scheme = True, 
                                    n_PS = n_PS, 
                                    s_min = s_min, s_max = s_max)

circ.draw('mpl', style = black)



#%%



print('Batch', end = ' ')
for k in range(batches):
    print(k+1, end = ' ')
    
    # Run the circuit 
    measure_qubits = [i for i in range(n)] + [n, n+12]           # The qubits we want to measure
    counts = run(circ, shots = shots, measure = measure_qubits)
    
    # Loop through the s values
    for s_label in range(0, n_PS):
      s = s_range[s_label]**2
      E = np.sqrt(s)/2
      
      # Find m_roof for later
      ms = [alt_single(s, M) for M in [M_A, M_Z]]
      m_roof = max([abs(m) for m in ms])
      
      # The M1 + M2 terms
      string_plus = '01' + binary(s_label, n)
      if string_plus in counts:
        N_plus = counts[string_plus]
      else:
        N_plus = 0
        print(f'string {string_plus} not found')
    
      # The M1 - M2 terms
      string_minus = '11' + binary(s_label, n)
      if string_minus in counts:
        N_minus = counts[string_minus]
      else:
        N_minus = 0
        print(f'string {string_minus} not found')
      
      # Compute full term output   
      outputs[s_label, k] = m_roof**2 * compensate**2 * E**4 * N_plus/shots
      
      # Compute interference output
      outputs_int[s_label, k] = m_roof**2 * compensate**2 * (E**4) * (N_plus - N_minus)/2/shots
      


for s_label in range(0, n_PS):
    
    s = s_range[s_label]**2
    E = np.sqrt(s)/2
    
    print('s =', s, '\tE =', E)
    print('---------------------------------------------------')
    
    print('Full term:')
    true = abs(get_photon_diagram(s) + get_Z_diagram(s))**2
    print('\tTrue value: ', true, end = '')
    out = np.mean(outputs[s_label, :])
    print('\tCircuit output: ', out)
    rel_error = abs(true - out)/abs(true)*100
    errors[s_label] = rel_error
    print('\tRelative error: ', rel_error, '%')
    print('\n')

    print('Interference:')
    true = get_interference(get_photon_diagram(s), get_Z_diagram(s))
    print('\tTrue value: ', true, end = '')
    out = np.mean(outputs_int[s_label, :])
    print('\tCircuit output: ', out)
    rel_error = abs(true - out)/abs(true)*100
    errors_int[s_label] = rel_error
    print('\tRelative error: ', rel_error, '%')
    print('\n')


#%%

fig, (ax1, ax2) = plt.subplots(2, 1, sharex = True)


# Full term plot

s_rough = np.linspace(s_min, s_max, n_PS)
s_fine = np.linspace(s_min, s_max, 1000)

# Plot the batches
for batch in outputs.T:
    ax1.scatter(s_rough, batch, alpha = 0.1, color = 'gray', s = 30)

# Plot the average of the batches
ax1.scatter(s_rough, [np.mean(out) for out in outputs], 
            color = 'purple', alpha = 0.75, s = 50, label = 'Circuit output average')

# Plot the true value
true = lambda s: abs(get_photon_diagram(s) + get_Z_diagram(s))**2
ax1.plot(s_fine, [true(s_sqrt**2) for s_sqrt in s_fine], 
         linestyle = '--', alpha = 0.3, color = 'black', label = 'True value')


ax1.set_yscale('log')
ax1.legend()

# ax1.set_xlabel(r'$\sqrt{s}$')
ax1.set_ylabel(r'$|\mathcal{M}_\gamma + \mathcal{M}_Z|^2$')
# plt.tight_layout()



# Interference plot

# Plot the batches
for batch in outputs_int.T:
    ax2.scatter(s_rough, batch, alpha = 0.1, color = 'gray', s = 30)

# Plot the average of the batches
ax2.scatter(s_rough, [np.mean(out) for out in outputs_int], 
            color = 'purple', alpha = 0.75, s = 50, label = 'Circuit output average')

# Plot the true value

true = lambda s: get_interference(get_photon_diagram(s), get_Z_diagram(s))
ax2.plot(s_fine, 
         [true(s_sqrt**2) for s_sqrt in s_fine], 
         linestyle = '--', alpha = 0.3, color = 'black', label = 'True value')


# plt.legend()
ax2.set_xlabel(r'$\sqrt{s}$')
ax2.set_ylabel(r'Int$(\mathcal{M}_\gamma, \mathcal{M}_Z)$')
# plt.tight_layout()

fig.tight_layout()



#%% 

# Include angle integration
s = 1000**2
n_angle_points = 8
print(f's = {s}\n')

theta_range = np.linspace(0, 2*np.pi, n_angle_points + 1)[:-1]

n_batches = 10

t_qubits = get_needed_qubits(n_angle_points)
inputs = generate_kinematic_inputs(s = s, theta = 0)

circ, compensate = generate_circuit(inputs, 
                                    n_particles = 2, 
                                    angle_integration = True, n_angle_points = n_angle_points)


circ.draw(output = 'mpl', style = black)

reg_dict = {reg.name : reg for reg in circ.qregs}
indices = lambda reg: [circ.find_bit(qubit).index for qubit in reg]

#%%

circ_outs = np.zeros((n_angle_points, n_batches))
circ_outs_int = np.zeros((n_angle_points, n_batches))

print('Batch: ', end = ' ')
for k in range(n_batches):
    print(k + 1)
    
    shots = 1e7
    measure = indices(reg_dict['T']) + indices(reg_dict['hit']) + indices((reg_dict['p']))
    counts = run(circ, shots = shots, measure = measure)

    for n in range(n_angle_points):
        theta = theta_range[n]
        spinor3 = get_spinor(E = np.sqrt(s)/2, s = '+', type = 'F', direction = 'out', theta = theta)
        spinor4 = get_spinor(E = np.sqrt(s)/2, s = '-', type = 'A', direction = 'out', theta = theta)
        
        # print(f'Theta: {theta_range[n]}')
        
        string_plus = '0' + '1' + binary(n, t_qubits)
        string_minus = '1' + '1' + binary(n, t_qubits)
        
        if string_plus in counts:
            N_plus = counts[string_plus]
        else:
            N_plus = 0.0
            
        if string_minus in counts:
            N_minus = counts[string_minus]     
        else:
            N_minus = 0.0
        
        # Full term
        circ_out = N_plus / shots * compensate**2 \
                        * np.linalg.norm(spinor3)**2 * np.linalg.norm(spinor4)**2
                        
        # Interference
        circ_out_int = (N_plus - N_minus)/2 / shots * compensate**2 \
                        * np.linalg.norm(spinor3)**2 * np.linalg.norm(spinor4)**2
            
        circ_outs[n, k] = circ_out
        circ_outs_int[n,k] = circ_out_int
        
        # print('Circ:', circ_out)
    
        # print('True:',
              # abs(get_photon_diagram(s, theta = theta) + get_Z_diagram(s, theta = theta))**2
              # )
        # print('\n')

#%%

theta_range_fine = np.linspace(0, 2*np.pi, 1000)

fig, (ax1, ax2) = plt.subplots(2, 1, sharex = True)

# Plot the batches
for batch in circ_outs.T:
    ax1.scatter(theta_range, batch, alpha = 0.1, color = 'gray', s = 30)

# Plot the average of the batches
ax1.scatter(theta_range, [np.mean(out) for out in circ_outs], 
            color = 'purple', alpha = 0.75, s = 50, label = 'Circuit output average')

ax1.plot(theta_range_fine, 
         [abs(get_photon_diagram(s, theta = theta) + get_Z_diagram(s, theta = theta))**2 for theta in theta_range_fine], 
         linestyle = '--', color = 'black', alpha = 0.3)

ax1.set_title(r'$\sqrt{s} =$' + str(np.sqrt(s)) + ' GeV')
# ax1.set_title(r'$\sqrt{s} = M_Z$')

# ax1.set_xlabel(r'$\theta$')
ax1.set_ylabel(r'$|\mathcal{M}_\gamma + \mathcal{M}_Z|^2$')


# Plot the batches
for batch in circ_outs_int.T:
    ax2.scatter(theta_range, batch, alpha = 0.1, color = 'gray', s = 30)

# Plot the average of the batches
ax2.scatter(theta_range, [np.mean(out) for out in circ_outs_int], 
            color = 'purple', alpha = 0.75, s = 50, label = 'Circuit output average')

ax2.plot(theta_range_fine, 
         [get_interference(get_photon_diagram(s, theta = theta), get_Z_diagram(s, theta = theta)) for theta in theta_range_fine], 
         linestyle = '--', color = 'black', alpha = 0.3)

ax2.set_xlabel(r'$\theta$')
ax2.set_ylabel(r'$Int(\mathcal{M}_\gamma, \mathcal{M}_Z)$')

fig.tight_layout()


#%%

# Multi-spinor gate implementation test

n_angle_points = 8

E = 40 # Arbitrary energy
t_qubits = get_needed_qubits(n_angle_points)
Theta = QuantumRegister(t_qubits, name = 'T')
v = QuantumRegister(2, name = 'v')
circ = QuantumCircuit(Theta, v)


circ.h(Theta)
circ.append(multi_spinor_gate(n_angle_points = n_angle_points, 
                              energy_dependent = True, 
                              E = E),
            register([Theta, v]))


psi = Statevector(circ)
circ.decompose().draw(output = 'mpl', style = black)

theta_range = np.linspace(0, 2*np.pi, n_angle_points + 1)[:-1]
for t in range(n_angle_points):
    spinor = get_spinor(E = E, 
                        s = '+', type = 'F', direction = 'out', theta = theta_range[t])
    
    label = binary(t, t_qubits)
    print(f"Theta = {theta_range[t]}")
    print("Spinor from circ: ", [np.sqrt(2**t_qubits)*np.linalg.norm(spinor)*psi[binary(j, 2) + label] for j in range(4)])
    print("True spinor: ", spinor)
    print("\n")
    



#%% 

inputs = generate_kinematic_inputs()
circ, compensate = generate_circuit(inputs, PS_integration = True, n_PS = 4, 
                                    angle_integration = True, n_angle_points = 2)

print(compensate)













