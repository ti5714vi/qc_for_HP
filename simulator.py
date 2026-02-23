# -*- coding: utf-8 -*-
"""
Created on Wed Feb  4 15:18:29 2026

    Full simulator function 

@author: Erik Bashore
"""

from main import alt_single
from Drell_Yan_circuit import *
from true_diagrams import *
from plotter import *
from constant_inputs import *


def simulate(n_particles = 2, s = 1000, 
             PS_integration = False, n_PS = 1, s_min = 80, s_max = 100,
             angle_integration = False, n_angle_points = 1, 
             shots = 1e7, n_batches = 10, 
             theta = np.pi/3):
    
    print('\nSIMULATING DY PROCESS')
    print('---------------------------------------------------------------')
    print('---------------------------------------------------------------')
    print(f'- Number of particles: {n_particles}')
    
    if PS_integration == True:
        print(f'- Phase space integration activated with sqrt(s) range [{s_min}, {s_max}] GeV')
        print(f'- Number of discrete PS points: {n_PS}')
        
        s_range_rough = np.linspace(s_min, s_max, n_PS)
        s_range_fine = np.linspace(s_min, s_max, 1000)
        s_qubits = get_needed_qubits(n_PS)
        
    else: 
        print('- sqrt{s} = ', np.sqrt(s))
        s_qubits = 0
        s_range_rough = np.array([s])
        
    
    if angle_integration == True:
        print(f'- Scattering angle integration activated with {n_angle_points} dicrete points')
        
        t_range_rough = np.linspace(0, np.pi, n_angle_points + 1)[:-1]
        t_range_fine = np.linspace(0, np.pi, 1000)
        
        t_qubits = get_needed_qubits(n_angle_points)
        
        
    else:
        print(f'- Scattering angle: {theta}')
        t_qubits = 0
        t_range_rough = np.array([theta])
    
    print(f'- Number of batches: {n_batches}')
    print(f'- Shots per batch: {int(shots)}')
    
    
    # Draw circuit
    # circ.draw(output = 'mpl', style = black)  
    
    print('---------------------------------------------------------------')
    print('---------------------------------------------------------------')
    
    # Kinematic inputs
    inputs = generate_kinematic_inputs(s = s, theta = theta)

    # Assemble circuit
    circ, compensate = generate_circuit(inputs, 
                                        n_particles = n_particles, 
                                        PS_integration = PS_integration, n_PS = n_PS, s_min = s_min, s_max = s_max,
                                        angle_integration = angle_integration, n_angle_points = n_angle_points)
    
    # print('Estimated number of gates in circuit: ', circ.decompose(reps = 2))
    
    # Create dictionary for circuit registers
    reg_dict = {reg.name : reg for reg in circ.qregs}
    
    # Output matrix
    # [full term/interference, batches, s, theta]
    circ_outs = np.zeros((2, n_batches, n_PS, n_angle_points)) 
    
    # Run circuit
    print('\nRunning circuit')
    for b in range(n_batches):
        print('Batch:', b+1, end = ' ')
        
        indices = lambda reg: [circ.find_bit(qubit).index for qubit in reg]
        
        measure = []
        if PS_integration == True:
            measure += indices(reg_dict['PS'])
            
        if angle_integration == True:
            measure += indices(reg_dict['T'])
    
        measure += indices(reg_dict['hit'])
        measure += indices(reg_dict['p'])
        
        # Run circuit        
        counts = run(circ, shots = shots, measure = measure)
        
        # Matrix of binary strings
        strings = np.array([[[sign + '1' + binary(t_label, t_qubits) + binary(s_label, s_qubits) for s_label in range(0, n_PS)] for t_label in range(0, n_angle_points)] for sign in ['0', '1']])
        
        # Ns = [N_plus, N_minus]
        N_plus = np.array([  [counts[string] if string in counts else 0 for string in strings[0, :, j]]  for j in range(n_PS)])
        N_minus = np.array([ [counts[string] if string in counts else 0 for string in strings[1, :, j]]  for j in range(n_PS)])
        
        Ns = np.array([N_plus, N_minus])
        
        for s_label in range(0, n_PS):
            for t_label in range(0, n_angle_points):
                
                # Full term
                circ_outs[0, b, s_label, t_label] = Ns[0, s_label, t_label]/shots * compensate[s_label]**2
                
                # Interference
                circ_outs[1, b, s_label, t_label] = (Ns[0, s_label, t_label] - Ns[1, s_label, t_label])/2/shots * compensate[s_label]**2    
    
    
    
    # ----- Order result ----- #
    print("\nCollecting results...")
    result = {}
    result['Full output'] = circ_outs
    
    # Average over batches
    result['Average'] = sum(circ_outs[0, :, :, :])/n_batches
    result['Average interference'] = sum(circ_outs[1, :, :, :])/n_batches
    
    # Include true matrices
    result['True matrix'] = np.array([[true(s_sqrt**2, theta, n_particles = n_particles) for theta in t_range_rough] \
                            for s_sqrt in s_range_rough])
        
    result['True matrix interference'] = np.array([[true_int(s_sqrt**2, theta, n_particles = n_particles) for theta in t_range_rough] \
                                                   for s_sqrt in s_range_rough])    
    
        
    # True smooth plots
    result['True s plot'] = np.array([[true(s_sqrt**2, theta, n_particles = n_particles) for theta in t_range_rough] \
                            for s_sqrt in s_range_fine])
        
    result['True s plot interference'] = np.array([[true_int(s_sqrt**2, theta, n_particles = n_particles) for theta in t_range_rough] \
                                                   for s_sqrt in s_range_fine]) 
    
    result['True theta plot'] = np.array([[true(s_sqrt**2, theta, n_particles = n_particles) for theta in t_range_fine] \
                            for s_sqrt in s_range_rough])
        
    result['True theta plot interference'] = np.array([[true_int(s_sqrt**2, theta, n_particles = n_particles) for theta in t_range_fine] \
                                                   for s_sqrt in s_range_rough]) 
        
    # Compute error matrix
    error_matrix = abs(result['Average'] - result['True matrix']) / result['True matrix'] * 100
    error_matrix_int = abs( (result['Average interference'] - result['True matrix interference']) / result['True matrix interference'] ) * 100
    
    
    # Mask not found points and outliers
    error_matrix = outlier_mask(error_matrix)
    error_matrix_int = outlier_mask(error_matrix_int)
    
    result['Error matrix'] = error_matrix
    result['Error matrix interference'] = error_matrix_int
    
    result['Data'] = {'PS points': n_PS, 'Angle points': n_angle_points, 's min': s_min, 's max': s_max,
                      'Batches': n_batches, 'Particles': n_particles}
    print("Simulation finished")
    return result
   

#%%


# ----- Run simulation ----- #
n_PS = 4
n_angle_points = 16
s_min = 85
s_max = 100
n_batches = 100
shots = 2e7

result = simulate(PS_integration = True, n_PS = n_PS, s_min = s_min, s_max = s_max,
         angle_integration = True, n_angle_points = n_angle_points,
         n_batches = n_batches, shots = shots)












