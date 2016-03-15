'''
Created on 14/01/2016

@author: iwelsh
'''
import numpy as np
import itertools
import random
import yaml

def dihedral_angle(a, b, c, d, period=[0,2*np.pi], scale='rad'):
    a = np.array(a)
    b = np.array(b)
    c = np.array(c)
    d = np.array(d)
    v1 = (a - b)/np.linalg.norm(a - b)
    v2 = (b - c)/np.linalg.norm(b - c)
    v3 = (c - d)/np.linalg.norm(c - d)
    n1 = np.cross(v1, v2)
    if np.linalg.norm(n1) > 0.0:
        n1 /= np.linalg.norm(n1)
    n2 = np.cross(v2, v3)
    if np.linalg.norm(n1) > 0.0:
        n2 /= np.linalg.norm(n2)
    m = np.cross(n1, v2)
    if np.linalg.norm(m) > 0.0:
        m /= np.linalg.norm(m)
    y, x = np.dot(m, n2), np.dot(n1, n2)
    phi = np.arctan2(y, x)
    if phi <= period[0]:
        phi += 2 * np.pi
    elif phi > period[1]:
        phi -= 2 * np.pi
    if scale == 'deg':
        phi *= 180/np.pi
        if int(round(phi,0)) == 360:
            phi -= 360
    return float(phi)

def kabsch_alignment(A, B, method=2):
    if method == 1:      # I think is wrong. Doesn;t line things up as nice as I would hope/expect
        covariance = np.dot(np.transpose(A), B)
        U, S, V = np.linalg.svd(covariance)
        # correct for right-handed coordinate system
        if (np.linalg.det(U) * np.linalg.det(V) < 0.0):
            S[-1] = -S[-1]
            U[:, -1] = -U[:, -1]
        # calc the rotation matrix
        rot_mat = np.dot(U,V)
        return rot_mat
    elif method == 2:
        covariance = np.dot(np.transpose(A), B)
        V, S, W_transpose = np.linalg.svd(covariance)
        d = np.sign(np.linalg.det(np.dot(np.transpose(W_transpose),np.transpose(V))))
        I = np.identity(3)
        I[2,2] = d
        U = np.dot(np.transpose(W_transpose), I)
        rot_mat = np.dot(U, np.transpose(V))
        return rot_mat
    
def dict_to_array(coordinates):
    A = []
    for atm in sorted(coordinates):
        A.append(np.asarray(coordinates[atm]['vec']))
    return np.asarray(A)
    
def aligned_rmsd(coordinates_A, coordinates_B):  
    # convert dict input to a matrix type thingy
    mat_A = dict_to_array(coordinates_A)
    mat_B = dict_to_array(coordinates_B)
    # centre both coordinate sets on 0,0,0
    mat_A -= sum(mat_A)/len(mat_A)
    mat_B -= sum(mat_B)/len(mat_B)
    # apply rotation to A to minimise RMSD
    mat_A = np.dot(mat_A, kabsch_alignment(mat_A, mat_B))
    # calculate the RMSD
    rmsd = 0.
    for a, b in zip(mat_A, mat_B):
        rmsd += sum([(a[i] - b[i])**2 for i in range(len(a))])
    return float(np.sqrt(rmsd/len(mat_A)))

def perform_fit(matrix, target, fit_phase, fit_method='svd', printme=False):
    if fit_method == 'qr':
        q,r = np.linalg.qr(matrix)
        d = np.dot(np.transpose(q),target)
        x = np.linalg.solve(r, d)
    elif fit_method == 'svd':
        x = np.linalg.lstsq(matrix, target)[0]
    elif fit_method == 'cauchy':
        initial_values = np.linalg.lstsq(matrix, target)[0]
        x = newton_raphson(matrix, target, initial_values)
    
    
    if printme:
        print(x)
    fit = []
    if not fit_phase:
        for k in x:
            if k > 0:
                fit.append(k)
                fit.append(0)
            else:
                fit.append(abs(k))
                fit.append(np.pi)
    else:
        dx, dy = fit_phase[0]*np.pi/180, fit_phase[1]*np.pi/180
        for i in range(0,len(x),2):
            k = (x[i]*np.cos(dx) + x[i+1]*np.cos(dy))**2
            k += (x[i]*np.sin(dx) + x[i+1]*np.sin(dy))**2
            k = np.sqrt(k)
            phase = np.arctan2(x[i]*np.sin(dx) + x[i+1]*np.sin(dy),x[i]*np.cos(dx) + x[i+1]*np.cos(dy))
            fit.append(k)
            fit.append(phase)
        
    return fit

def test_case(coefficient_matrix, target_vector, beta_values):
    term_vector = np.zeros(3)
    term_vector[0] = beta_values[0]**2 + beta_values[1]**2 + beta_values[2]**2 -3
    term_vector[1] = beta_values[0]**2 + beta_values[1]**2 - beta_values[2] - 1
    term_vector[2] = beta_values[0] + beta_values[1] + beta_values[2] - 3
    return -term_vector

def test_case_derivative(coefficient_matrix, target_vector, beta_values):
    jacobian = np.zeros((3,3))
    jacobian[0][0] = 2*beta_values[0]
    jacobian[0][1] = 2*beta_values[1]
    jacobian[0][2] = 2*beta_values[2]
    jacobian[1][0] = 2*beta_values[0]
    jacobian[1][1] = 2*beta_values[1]
    jacobian[1][2] = -1
    jacobian[2][0] = 1
    jacobian[2][1] = 1
    jacobian[2][2] = 1
    return jacobian

def arctan(coefficient_matrix, target_vector, beta_values):
    term_vector = np.zeros(len(beta_values))
    for i in range(len(term_vector)):
        sum_term = 0
        for j in range(len(target_vector)):
            common_term = np.sum(coefficient_matrix[j] * beta_values) - target_vector[j]
            sum_term += (2 * common_term * coefficient_matrix[j][i]) / (common_term ** 4 + 1)
        term_vector[i] = -sum_term
    
    return term_vector

def arctan_derivative(coefficient_matrix, target_vector, beta_values):
    # return the Jacobian
    jacobian = np.zeros((len(beta_values),len(beta_values)))
    for row in range(len(beta_values)):
        for col in range(len(beta_values)):
            sum_term = 0
            for conf in range(len(target_vector)):
                const_term = (np.sum(coefficient_matrix[conf] * beta_values) - target_vector[conf])**4
                multi = coefficient_matrix[conf][row] * coefficient_matrix[conf][col]
                sum_term += (2 * multi)/(const_term + 1) - (8 * multi * const_term)/((const_term + 1)**2)
            jacobian[row][col] = sum_term
            
    return jacobian

def cauchy(coefficient_matrix, target_vector, beta_values):
    # return the jacobian solution
    term_vector = np.zeros(len(beta_values))   # this is going to become b in the new equation to solve
    #print(len(target_vector))
    for i in range(len(term_vector)):        # i give the i-th simultaneous equation in the result
        sum_term = 0
        for j in range(len(target_vector)):  # j give the j-th simultaneous equation in the source
            common_term = np.sum(coefficient_matrix[j] * beta_values) - target_vector[j]
            sum_term += (2 * common_term * coefficient_matrix[j][i]) / (common_term ** 2 + 1)
        term_vector[i] = -sum_term
    
    return term_vector

def cauchy_derivative(coefficient_matrix, target_vector, beta_values):
    # return the Jacobian
    jacobian = np.zeros((len(beta_values),len(beta_values)))
    for row in range(len(beta_values)):
        for col in range(len(beta_values)):
            sum_term = 0
            for conf in range(len(target_vector)):
                const_term = (np.sum(coefficient_matrix[conf] * beta_values) - target_vector[conf])**2
                multi = coefficient_matrix[conf][row] * coefficient_matrix[conf][col]
                sum_term += (2 * multi)/(const_term + 1) - (4 * multi * const_term)/((const_term + 1)**2)
            jacobian[row][col] = sum_term
            
    #print(jacobian)
    return jacobian

def lls(coefficient_matrix, target_vector, beta_values):
    term_vector = np.zeros(len(beta_values))
    for i in range(len(term_vector)):
        sum_term = 0
        for j in range(len(target_vector)):
            common_term = np.sum(coefficient_matrix[j] * beta_values) - target_vector[j]
            sum_term += 2 * common_term * coefficient_matrix[j][i]
        term_vector[i] = -sum_term
    #print(term_vector)
    return term_vector

def lls_derivative(coefficient_matrix, target_vector, beta_values):
    jacobian = np.zeros((len(beta_values),len(beta_values)))
    for row in range(len(beta_values)):
        for col in range(len(beta_values)):
            sum_term = 0
            for conf in range(len(target_vector)):
                sum_term += coefficient_matrix[conf][row] * coefficient_matrix[conf][col] * 2
            jacobian[row][col] = sum_term
    #print(jacobian)
    return jacobian

def smooth_approx(coefficient_matrix, target_vector, beta_values):
    term_vector = np.zeros(len(beta_values))   # this is going to become b in the new equation to solve
    for i in range(len(term_vector)):        # i give the i-th simultaneous equation in the result
        sum_term = 0
        for j in range(len(target_vector)):  # j give the j-th simultaneous equation in the source
            common_term = np.sum(coefficient_matrix[j] * beta_values) - target_vector[j]
            sum_term += (2 * common_term * coefficient_matrix[j][i]) / np.sqrt((common_term ** 2 + 1))
        term_vector[i] = -sum_term
    
    return term_vector

def smooth_approx_derivative(coefficient_matrix, target_vector, beta_values):
    # return the Jacobian
    jacobian = np.zeros((len(beta_values),len(beta_values)))
    for row in range(len(beta_values)):
        for col in range(len(beta_values)):
            sum_term = 0
            for conf in range(len(target_vector)):
                const_term = (np.sum(coefficient_matrix[conf] * beta_values) - target_vector[conf])**2
                multi = coefficient_matrix[conf][row] * coefficient_matrix[conf][col]
                sum_term += (2 * multi)/np.sqrt((const_term + 1)) - (4 * multi * const_term)/np.power((const_term + 1),3/2)
            jacobian[row][col] = sum_term
            
    #print(jacobian)
    return jacobian

def calc_norm_cauchy(matrix, target, values):
    norm = 0
    for row in range(len(target)):
        norm += np.log((np.sum(matrix[row] * values) - target[row])**2 +1)
    return norm

def calc_norm_arctan(matrix, target, values):
    norm = 0
    for row in range(len(target)):
        norm += np.arctan((np.sum(matrix[row] * values) - target[row])**2)
    return norm

def calc_absolute_deviation(matrix, target, values):
    absolute = 0
    for row in range(len(target)):
        absolute += np.abs(np.sum(matrix[row] * values) - target[row])
    return absolute

def newton_raphson(matrix, target, initial, fit_method='cauchy', max_iter=50, tolerance=1e-14):
    #initial = np.array([4.068553366827293, 2.4012157, 3.4419819, 0.4587839, 0.2616956, 0.3234017])
    #initial = np.array([4.0321680, 2.4212157, 3.3019819, 0.2587839, 0.216956, -0.])
    #initial = np.array([-1.63771227, -4.06120612,  3.70036194,  4.00151597,  2.64167739, -1.66428148,
    #                     0.02128043, -0.92143359, -0.1749524,   0.62808971,  0.20873093, -0.37133858])
    
    #print(initial)
    #print(matrix)
    #print(target)
    for _ in range(max_iter):
        #term_vector = smooth_approx(matrix, target, initial)
        #jacobian = smooth_approx_derivative(matrix, target, initial)
        #term_vector = lls(matrix, target, initial)
        #jacobian = lls_derivative(matrix, target, initial)
        term_vector = cauchy(matrix, target, initial)
        jacobian = cauchy_derivative(matrix, target, initial)
        #term_vector = arctan(matrix, target, initial)
        #jacobian = arctan_derivative(matrix, target, initial)
        #term_vector = test_case(matrix, target, initial)
        #jacobian = test_case_derivative(matrix, target, initial)
        iter_s = np.linalg.solve(jacobian, term_vector)
        #print(iter_s)
        #break
        current = initial + iter_s
        
        initial_norm = calc_norm_cauchy(matrix, target, initial)
        current_norm = calc_norm_cauchy(matrix, target, current)
        
        if current_norm > initial_norm:
            # scale the step size down until it does decrease the target function
            while current_norm > initial_norm:# and norm_iter < max_iter:
                norm_ratio = (current_norm/initial_norm)**2
                iter_s *= (np.sqrt(1 + 6 * norm_ratio) - 1) / (3 * norm_ratio)
                current = initial + iter_s
                initial_norm = calc_norm_cauchy(matrix, target, initial)
                current_norm = calc_norm_cauchy(matrix, target, current)
        #elif current_norm < initial_norm:
        #    previous_iter = iter_s
        #    print('scaling up')
        #    # scale the step size up until it no longer decreases the target function
        #    while current_norm < initial_norm:
        #        previous_iter = iter_s
        #        norm_ratio = (current_norm/initial_norm)**2
        #        iter_s /= (np.sqrt(1 + 6 * norm_ratio) - 1) / (3 * norm_ratio)
        #        current = initial + iter_s
        #        initial_norm = calc_norm_cauchy(matrix, target, initial)
        #        current_norm = calc_norm_cauchy(matrix, target, current)
        #        
        #    current = initial + previous_iter
                
        
        allWithinTolerance = True
        for i in range(len(current)):
            if np.abs(current[i] - initial[i]) > tolerance:
                allWithinTolerance = False
                break 
        initial = current
        
        if allWithinTolerance:
            #print("Converged after {} iterations".format(iter))
            break
        
    return initial

def fit_torsion_terms_lls(energies, angles, fit_phase=False, allowed_multiplicites=list(range(1,7)), printme=False, method='svd'):
    
    assert len(energies) == len(angles)
    
    if fit_phase: 
        multiplier = 2
    else: 
        multiplier = 1
    
    # set up the coefficient and energy target matrices
    t, m = 1,len(allowed_multiplicites)*multiplier
    coefficient_matrix = np.zeros((1,m*t))
    target_vector = np.array((0))
    
    # find the  minimum energy
    #min_target = min([energies[x] for x in energies])
    
    for x in sorted(energies):
        #count += 1
        energy_vector = np.zeros((1,m*t))
        for j in range(m):
            #j = j# + m*x
            multiplicity = allowed_multiplicites[int(np.floor(j/multiplier))]
            angle = angles[x]*np.pi/180
            if not fit_phase:
                torsion_value = np.cos(multiplicity*angle)
            else:
                if j % 2 == 0:
                    phase = fit_phase[0] * np.pi / 180
                else:
                    phase = fit_phase[1] * np.pi / 180
                torsion_value = np.cos(multiplicity*angle + phase)
            
            energy_vector[0,j] = torsion_value
        coefficient_matrix = np.vstack((coefficient_matrix, energy_vector))
        target_vector = np.hstack((target_vector,np.array((energies[x]))))
        
    coefficient_matrix = coefficient_matrix[1:]
    target_vector = target_vector[1:]
    
    target_vector -= np.mean(target_vector)
        
    # subtract the mean of cos(mi\thetaj) from each value of coefficient_matrix
    for col in range(len(coefficient_matrix[0])):
        col_vals = [];
        for row in range(len(coefficient_matrix)):
            col_vals.append(coefficient_matrix[row][col])
        mean_col = np.mean(col_vals)
        for row in range(len(coefficient_matrix)):
            coefficient_matrix[row][col] -= mean_col
    
    fit = perform_fit(coefficient_matrix, target_vector, fit_phase, fit_method=method, printme=printme)
    
    
    dict_fit = {}
    for i in range(len(allowed_multiplicites)):
        j = i*2
        dict_fit[allowed_multiplicites[i]] = {'k':float(fit[j]),
                                              'delta':float(fit[j+1])*180/np.pi,
                                              'm':allowed_multiplicites[i]}
    return dict_fit

def two_pass_lls(energies, angles, limit, phase=False):
    return multi_pass_lls(energies, angles, limit, phase=phase)
    if phase:
        assert isinstance(phase, (list,tuple)) and len(phase) == 2
    fit1 = fit_torsion_terms_lls(energies, angles, fit_phase=phase)
    allowed_two = sorted(fit1, key=lambda x:fit1[x]['k'], reverse=True)[0:limit]
        
    fit2 = fit_torsion_terms_lls(energies, angles, fit_phase=phase, allowed_multiplicites=allowed_two)
    
    return fit2

def energy_at_x(fit, x):
    energy = 0
    for y in fit:
        energy += fit[y]['k']*(np.cos(x*fit[y]['m']*np.pi/180 + fit[y]['delta']*np.pi/180))
    return energy

def energy_at_x2(fit, x):
    energy = 0
    for y in fit:
        energy += fit[y]['k']*(1 + np.cos(x*fit[y]['m']*np.pi/180 + fit[y]['delta']*np.pi/180))
    return energy

def fit_deviation(energies, angles, fit):
    errors = []
    mean_energy = np.mean([energies[x] for x in energies])
    for i in angles:
        errors.append((energies[i] -mean_energy - energy_at_x(fit, angles[i]))**2)
    
    rmsd = np.sqrt(np.mean(errors))
    
    return rmsd

def multi_pass_lls(energies, angles, limit, allowed_multiplicites=list(range(1,7)), phase=False):
    if phase:
        assert isinstance(phase, (list,tuple)) and len(phase) == 2
        
    best_fit = {}
    best_deviatons = 1e80
        
    multiplicites_to_use = itertools.combinations(allowed_multiplicites, limit)
    for i in multiplicites_to_use:
        fit = fit_torsion_terms_lls(energies, angles, fit_phase=phase, allowed_multiplicites=i, baseline_shift=False)
        deviation = fit_deviation(energies, angles, fit)#[0][1]**2
       # print("rmsd: ", deviation, "multis: ", i)
        if deviation < best_deviatons:
            #print("New best deviation of {} with {}".format(deviation, i))
            best_deviatons = deviation
            best_fit = fit
            
    return best_fit
            
    
    

def one_pass_lls(energies, angles, limit=None, phase=False, method='svd', printme=False):
    if phase:
        assert isinstance(phase, (list,tuple)) and len(phase) == 2
    fit = fit_torsion_terms_lls(energies, angles, fit_phase=phase, method=method, printme=printme)
    
    if limit is not None:
        sorted_list = sorted(fit, key=lambda x:fit[x]['k'], reverse=True)[:limit]
        for k in set(fit.keys()).difference(sorted_list):
            del fit[k]
    return fit
    
def prune_data(energies, angles, method='fit'): # method can be fit or derivative
    # prune the data to remove outliers.
    # perform a phased, one_pass_lls fit using all points
    # calculate the variance between each point and it's fitted value to find the most variant points
    # if the most variant point is above some threshold, remove it and refit
    # continue until most variant point is below threshold
    start_len = len(energies)
    remove_count = 0
    while True and method == 'fit' and remove_count < 0.25*start_len:
        mean_energy = np.mean([energies[x] for x in energies])
        for x in energies:
            energies[x] -= mean_energy
        fit = one_pass_lls(energies, angles, phase=[0,90], method='cauchy')
        current_variance = 0
        current_variant_point = -1
        total_variance = []
        for ang in sorted(angles):
            variance = np.abs(energy_at_x(fit, angles[ang]) - energies[ang])
            total_variance.append(variance)
            if variance > current_variance:
                current_variance = variance
                current_variant_point = ang
        if current_variance < 5:
            break
        del energies[current_variant_point], angles[current_variant_point]
        remove_count += 1
        #print('removed data at point ', current_variant_point)
        
    
    while True and method == 'derivative' and remove_count < 0.25*start_len:
        derivatives, d_angles = derivative(energies, angles)
        d_mean = np.mean([derivatives[x] for x in derivatives])
        d_std = np.std([derivatives[x] for x in derivatives])
        #print(sorted(derivatives, key=lambda x: derivatives[x])[-1],sorted(derivatives, key=lambda x: derivatives[x])[0])
        largest_d = sorted(derivatives, key=lambda x: derivatives[x])[-1]
        if derivatives[largest_d] < d_mean + 3 * d_std:
            break
        del energies[largest_d], angles[largest_d]
        remove_count += 1
        print('removed data at point ', largest_d)
        
    
    return energies, angles

def rmsd_fit_to_raw(energiesQM, energiesMD, anglesQM, anglesMD, fit):
    fit_raw_energies = {}
    meanQM = np.mean([energiesQM[x] for x in energiesQM])
    meanMD = np.mean([energiesMD[x] for x in energiesMD])
    for i in energiesMD:
        if i not in energiesQM:
            continue
        avg_angle = (anglesQM[i] + anglesMD[i])/2
        fit_raw_energies[i] = energy_at_x(fit, avg_angle) + energiesMD[i] - meanMD
    meanFIT = np.mean([fit_raw_energies[x] for x in fit_raw_energies])
    deviations = []
    for i in fit_raw_energies:
        deviations.append((energiesQM[i] - meanQM - fit_raw_energies[i] - meanFIT)**2)
    rmsd = np.sqrt(np.mean(deviations))
    return rmsd

def rawminusraw(energiesA, energiesB, anglesA, anglesB, limit, phase, fh):
    diffEnergies, diffAngles = {}, {}
    for i in energiesA:
        if i in energiesB:
            diffEnergies[i] = energiesA[i] - energiesB[i]
            diffAngles[i] = anglesA[i]
    fit = one_pass_lls(diffEnergies, diffAngles, limit, phase)
    for term in fit:
        print('k{m:<2} = {k:11.7f}, delta{m:<2} = {delta:11.7f}'.format(**fit[term]), file=fh)
    print("RMSD fit = {:11.7f} \n".format(fit_deviation(diffEnergies, diffAngles, fit)), file=fh)
    print("RMSD total = {:8.4f}".format(rmsd_fit_to_raw(energiesA, energiesB, 
                                                                            anglesA, anglesB, fit)))#, file=fh)
    
def fitminusraw(energiesA, energiesB, anglesA, anglesB, fitset, limit, rawphase, fitphase, fh):
    if fitset == 'A':
        fit = one_pass_lls(energiesA, anglesA, None, rawphase)
    elif fitset == 'B':
        fit = one_pass_lls(energiesB, anglesB, None, rawphase)
        
    energiesDIFF, anglesDIFF = {}, {}
    if fitset == 'B':
        angles = [(anglesA[x],x) for x in anglesA]
    elif fitset == 'A':
        angles = [(anglesB[x],x) for x in anglesB]
    
    for i in angles:
        if fitset == 'A':
            energiesDIFF[i[1]] = energy_at_x(fit, i[0]) - energiesB[i[1]]
            anglesDIFF[i[1]] = i[0]
        elif fitset == 'B':
            energiesDIFF[i[1]] = energiesA[i[1]] - energy_at_x(fit, i[0]) 
            anglesDIFF[i[1]] = i[0]
    difffit = one_pass_lls(energiesDIFF, anglesDIFF, limit, fitphase)
    
    if fitset == 'A':
        print("\nQM Fit:", file=fh)
        for term in fit:
            print('k{m:<2} = {k:11.7f}, delta{m:<2} = {delta:11.7f}'.format(**fit[term]), file=fh)
        print("RMSD fit = {:11.7f} \n".format(fit_deviation(energiesA, anglesA, fit)), file=fh)
    elif fitset == 'B':
        print("\nMD Fit:", file=fh)
        for term in fit:
            print('k{m:<2} = {k:11.7f}, delta{m:<2} = {delta:11.7f}'.format(**fit[term]), file=fh)
        print("RMSD fit = {:11.7f} \n".format(fit_deviation(energiesB, anglesB, fit)), file=fh)
    print("\nDIFF Fit:", file=fh)
    for term in difffit:
        print('k{m:<2} = {k:11.7f}, delta{m:<2} = {delta:11.7f}'.format(**difffit[term]), file=fh)
    print("RMSD fit = {:11.7f} \n".format(fit_deviation(energiesDIFF, anglesDIFF, difffit)), file=fh)
    print("RMSD total = {:8.4f}".format(rmsd_fit_to_raw(energiesA, energiesB, 
                                                                        anglesA, anglesB, difffit)))#, file=fh)
    
    
    
def fitminusfit(energiesA, energiesB, anglesA, anglesB, limit, rawphase, fitphase, fh):
    afit = one_pass_lls(energiesA, anglesA, None, rawphase)
    bfit = one_pass_lls(energiesB, anglesB, None, rawphase)
    energiesDIFF, anglesDIFF = {}, {}
    for i in range(360):
        energiesDIFF[i] = energy_at_x(afit, i) - energy_at_x(bfit, i)
        anglesDIFF[i] = i
    difffit = one_pass_lls(energiesDIFF, anglesDIFF, limit, fitphase)
    
    print("\nQM Fit:", file=fh)
    for term in afit:
        print('k{m:<2} = {k:11.7f}, delta{m:<2} = {delta:11.7f}'.format(**afit[term]), file=fh)
    print("RMSD fit = {:11.7f} \n".format(fit_deviation(energiesA, anglesA, afit)), file=fh)
    
    print("\nMD Fit:", file=fh)
    for term in bfit:
        print('k{m:<2} = {k:11.7f}, delta{m:<2} = {delta:11.7f}'.format(**bfit[term]), file=fh)
    print("RMSD fit = {:11.7f} \n".format(fit_deviation(energiesB, anglesB, bfit)), file=fh)    
        
    print("\nDIFF Fit:", file=fh)
    for term in difffit:
        print('k{m:<2} = {k:11.7f}, delta{m:<2} = {delta:11.7f}'.format(**difffit[term]), file=fh)
    print("RMSD fit = {:11.7f} \n".format(fit_deviation(energiesDIFF, anglesDIFF, difffit)), file=fh)
    print("RMSD total = {:8.4f}".format(rmsd_fit_to_raw(energiesA, energiesB, 
                                                                        anglesA, anglesB, difffit)))#, file=fh)
    
    

def main(f, filePath):
    #import yaml
    
    #with open('/Users/iwelsh/GitHub/ExtractedData/Torsions/AdditionalFixedTorsions/Original/AMINO2/energies.qm.txt', 'r') as fh:
    with open(filePath + "energies.qm.txt", 'r') as fh:
        energiesQM, anglesQM = yaml.load(fh)
    with open(filePath + "energies.aa.txt", 'r') as fh:
        energiesAA, anglesAA = yaml.load(fh)
    with open(filePath + "energies.ua.txt", 'r') as fh:
        energiesUA, anglesUA = yaml.load(fh)
        
    print("Non-phased fitting to raw QM - raw AA:", file=f)
    rawminusraw(energiesQM, energiesAA, anglesQM, anglesAA, 3, False, f)
    print("Phased fitting to raw QM - raw AA:", file=f)
    rawminusraw(energiesQM, energiesAA, anglesQM, anglesAA, 3, [0,90], f)
    print("Non-phased fitting to raw QM - non-phased AA:", file=f)
    fitminusraw(energiesQM, energiesAA, anglesQM, anglesAA, 'B', 3, False, False, f)
    print("Phased fitting to raw QM - non-phased AA:", file=f)
    fitminusraw(energiesQM, energiesAA, anglesQM, anglesAA, 'B', 3, False, [0,90], f)
    print("Non-phased fitting to raw QM - phased AA:", file=f)
    fitminusraw(energiesQM, energiesAA, anglesQM, anglesAA, 'B', 3, [0,90], False, f)
    print("Phased fitting to raw QM - phased AA:", file=f)
    fitminusraw(energiesQM, energiesAA, anglesQM, anglesAA, 'B', 3, [0,90], [0,90], f)
    print("Non-phased fitting to non-phased QM - raw AA:", file=f)
    fitminusraw(energiesQM, energiesAA, anglesQM, anglesAA, 'A', 3, False, False, f)
    print("Phased fitting to non-phased QM - raw AA:", file=f)
    fitminusraw(energiesQM, energiesAA, anglesQM, anglesAA, 'A', 3, False, [0,90], f)
    print("Non-phased fitting to non-phased QM - non-phased AA:", file=f)
    fitminusfit(energiesQM, energiesAA, anglesQM, anglesAA, 3, False , False, f) 
    print("Non-phased fitting to phased QM - raw AA:", file=f)
    fitminusraw(energiesQM, energiesAA, anglesQM, anglesAA, 'A', 3, [0,90], False, f)
    print("Phased fitting to phased QM - raw AA:", file=f)
    fitminusraw(energiesQM, energiesAA, anglesQM, anglesAA, 'A', 3, [0,90], [0,90], f)
    print("Non-phased fitting to phased QM - phased AA:", file=f)
    fitminusfit(energiesQM, energiesAA, anglesQM, anglesAA, 3, [0,90] , False, f) 
    print("Phased fitting to phased QM - phased AA:", file=f)    
    fitminusfit(energiesQM, energiesAA, anglesQM, anglesAA, 3, [0,90] , [0,90], f)
    
    print("\n=======================================\n", file=f)
    print("Non-phased fitting to raw QM - raw UA:", file=f)
    rawminusraw(energiesQM, energiesUA, anglesQM, anglesUA, 3, False, f)
    print("Phased fitting to raw QM - raw UA:", file=f)
    rawminusraw(energiesQM, energiesUA, anglesQM, anglesUA, 3, [0,90], f)
    print("Non-phased fitting to raw QM - non-phased UA:", file=f)
    fitminusraw(energiesQM, energiesUA, anglesQM, anglesUA, 'B', 3, False, False, f)
    print("Phased fitting to raw QM - non-phased UA:", file=f)
    fitminusraw(energiesQM, energiesUA, anglesQM, anglesUA, 'B', 3, False, [0,90], f)
    print("Non-phased fitting to raw QM - phased UA:", file=f)
    fitminusraw(energiesQM, energiesUA, anglesQM, anglesUA, 'B', 3, [0,90], False, f)
    print("Phased fitting to raw QM - phased UA:", file=f)
    fitminusraw(energiesQM, energiesUA, anglesQM, anglesUA, 'B', 3, [0,90], [0,90], f)
    print("Non-phased fitting to non-phased QM - raw UA:", file=f)
    fitminusraw(energiesQM, energiesUA, anglesQM, anglesUA, 'A', 3, False, False, f)
    print("Phased fitting to non-phased QM - raw UA:", file=f)
    fitminusraw(energiesQM, energiesUA, anglesQM, anglesUA, 'A', 3, False, [0,90], f)
    print("Non-phased fitting to non-phased QM - non-phased UA:", file=f)
    fitminusfit(energiesQM, energiesUA, anglesQM, anglesUA, 3, False , False, f)
    print("Non-phased fitting to phased QM - raw UA:", file=f)
    fitminusraw(energiesQM, energiesUA, anglesQM, anglesUA, 'A', 3, [0,90], False, f)
    print("Phased fitting to phased QM - raw UA:", file=f)
    fitminusraw(energiesQM, energiesUA, anglesQM, anglesUA, 'A', 3, [0,90], [0,90], f)
    print("Non-phased fitting to phased QM - phased UA:", file=f)
    fitminusfit(energiesQM, energiesUA, anglesQM, anglesUA, 3, [0,90] , False, f)
    print("Phased fitting to phased QM - phased UA:", file=f)    
    fitminusfit(energiesQM, energiesUA, anglesQM, anglesUA, 3, [0,90] , [0,90], f)
    
    
    
    #energies, angles = prune_data(energies, angles)
    #energies, angles = {}, {}
    #phases = [random.triangular(-np.pi/2, np.pi/2) for _ in range(12)]
    #amplitudes = []
    #larges = [1,3,6]
    #mediums = [4,5,2]
    #random.shuffle(larges)
    #random.shuffle(mediums)
    #for i in range(12):
    #    if i + 1 in larges[0:2]:
    #        amplitudes.append(random.triangular(8.5,20.5))
    #    elif i + 1 in mediums[0:2]:
    #        amplitudes.append(random.triangular(4.5,8.5))
    #    else:
    #        amplitudes.append(random.triangular(0.5,1.0))
    #for i in range(12):
    #    print('k{0:<2} = {1:11.7f}, delta{0:<2} = {2:11.7f}'.format(i+1, amplitudes[i], phases[i]*180/np.pi), file=f)
    #print("___________________________________________", file=f)
    #for i in range(360):
    #    energies[i] = np.sum([amplitudes[x]*(1+np.cos((x+1)*i*np.pi/180 - phases[x])) for x in range(12)]) * random.uniform(0.9, 1.1)
    #    angles[i] = i
    ##angles = list(range(360))
    #print(fit)
#    opnophasefit = one_pass_lls(energies, angles, 3, False)
#    opphasefit = one_pass_lls(energies, angles, 3, [0,90])
#    
#    tpnophasefit = two_pass_lls(energies, angles, 3, False)
#    tpphasefit = two_pass_lls(energies, angles, 3, [0,90])
#    
#    mpnophasefit = multi_pass_lls(energies, angles, 3, phase=False)
#    mpphasefit = multi_pass_lls(energies, angles, 3, phase=[0,90])
#    rmsds = []
#    for fit in [opnophasefit, tpnophasefit, mpnophasefit,opphasefit, tpphasefit, mpphasefit]:
#        rmsds.append(fit_deviation(energies, angles, fit))
#    print("One pass non-phased fitting:", file=f)
#    for term in opnophasefit:
#        print('k{m:<2} = {k:11.7f}, delta{m:<2} = {delta:11.7f}'.format(**opnophasefit[term]), file=f)
#    print("RMSD = {:11.8f}\n".format(rmsds[0]), file=f)
#    
#    print("Two pass non-phased fitting:", file=f)
#    for term in tpnophasefit:
#        print('k{m:<2} = {k:11.7f}, delta{m:<2} = {delta:11.7f}'.format(**tpnophasefit[term]), file=f)
#    print("RMSD = {:11.8f}\n".format(rmsds[1]), file=f)
#    
#    print("Multi pass non-phased fitting:", file=f)
#    for term in mpnophasefit:
#        print('k{m:<2} = {k:11.7f}, delta{m:<2} = {delta:11.7f}'.format(**mpnophasefit[term]), file=f)
#    print("RMSD = {:11.8f}\n".format(rmsds[2]), file=f)
#    
#    print("One pass phased fitting:", file=f)
#    for term in opphasefit:
#        print('k{m:<2} = {k:11.7f}, delta{m:<2} = {delta:11.7f}'.format(**opphasefit[term]), file=f)
#    print("RMSD = {:11.8f}\n".format(rmsds[3]), file=f)
#    
#    print("Two pass phased fitting:", file=f)
#    for term in tpphasefit:
#        print('k{m:<2} = {k:11.7f}, delta{m:<2} = {delta:11.7f}'.format(**tpphasefit[term]), file=f)
#    print("RMSD = {:11.8f}\n".format(rmsds[4]), file=f)
#    
#    print("Multi pass phased fitting:", file=f)
#    for term in mpphasefit:
#        print('k{m:<2} = {k:11.7f}, delta{m:<2} = {delta:11.7f}'.format(**mpphasefit[term]), file=f)
#    print("RMSD = {:11.8f}\n".format(rmsds[5]), file=f)
    #fit = one_pass_lls(energies, angles, 3, phase=[0,90])
    #fit_deviation(energies, angles, fit)
    #fit = two_pass_lls(energies, angles, 3)
    #for i in fit:
    #    print(fit[i])
    #return rmsds

def main2():
    with open('/Users/iwelsh/GitHub/ExtractedData/Torsions/AdditionalFixedTorsions/Original/METH0/energies.ua.txt', 'r') as fh:
        energies, angles = yaml.load(fh)
    energies, angles = prune_data(energies, angles)

def derivative(energies, angles):
    derivatives, d_angles = {},{}
    sorted_derivs = sorted(angles, key=lambda x:angles[x])

    for i in range(len(sorted_derivs)):
        ang = sorted_derivs[i]
        lead_ang = sorted_derivs[i-1]
        try:
            trail_ang = sorted_derivs[i+1]
        except IndexError:
            trail_ang = sorted_derivs[0]
        
        lead_diff = min((np.abs(angles[ang] - angles[lead_ang]),np.abs(angles[ang] - angles[lead_ang] + 360)))
        trail_diff = min((np.abs(angles[ang] - angles[trail_ang]),np.abs(angles[ang] - angles[trail_ang] - 360)))
        if lead_diff > 90:
            print('lead', ang, lead_ang, lead_diff)
        if trail_diff > 90:
            print('trail', ang, lead_ang, trail_diff)
        derivatives[5*i] = (energies[ang] - energies[lead_ang])/(lead_diff)
                            
        d_angles[5*i] = ang
    return derivatives, d_angles

def cycle_iterator(start_idx, d_list):
    for i in range(len(d_list)):
        idx = i + start_idx
        if idx >= len(d_list):
            idx -= len(d_list)
        yield d_list[idx]

if __name__ == '__main__':
    from time import process_time
    import os
    import sys
    mean_rmsds = [0,0,0,0,0,0]
    #main2()
    #sys.exit()
    #iterations = 50000
    #roots = ['30PsiAllFixed', '60PsiBothAllFixed', 'AllHeavyFixed', 'Figures_Pruned', '30PsiBothAllFixed', 
    #         '60PsiOriginal', 'BothAllFixed', 'Figures_Raw', 'Original', 'ThetaAllFixed', '30PsiBothCarbonFixed', 
    #         'BothCarbonFixed', 'FunctionalFixed', 'ThetaCarbonFixed', '30PsiCarbonFixed', 'PsiAllFixed', '30PsiOriginal', 
    #         'AllCarbonFixed', 'PsiCarbonFixed', '60PsiAllFixed', 'AllFixed']
    roots = ['Original']
    mols = ['AMINO','CHLORO','HYDRO','METH','THIO']
    types = [-1,0,1,2]
    levls = ['aa','ua','qm']
    
    i = 0
    
    ks = [4.0, 2.5212157, 3.4019819, 0.2887839, 0.2916956, 0.1534017]
    deltas = [56.6217515, -55.8958340, 5.9586481, -156.6942762, 47.6905846, 76.4370030]
    #deltas = [0,0,0,0,0,180]
    energies, angles = {}, {}
    for i in range(0,360,5):
        energies[i] = 0
        for j in range(len(ks)):
            energies[i] += float(ks[j]*(1 + np.cos((j+1)*i*np.pi/180 - deltas[j]*np.pi/180)))
            #break
        angles[i] = i
    mean_energy = np.mean([energies[x] for x in energies])
    for x in energies:
        energies[x] -= mean_energy
    energies[90] *= -8
    energies[20] *= -6
    
    svd_fit = one_pass_lls(energies, angles, phase=[0,90])
    #print(svd_fit)
    cauchy_fit = fit_torsion_terms_lls(energies, angles, fit_phase=[0,90], method='cauchy')
    fit_angles, svd_energies, cauchy_energies = [], [], []
    for i in range(360):
        fit_angles.append(i)
        svd_energies.append(energy_at_x(svd_fit, i))
        cauchy_energies.append(energy_at_x(cauchy_fit, i))
    
    from plotting import plot_scatters
    plot_scatters([{'x':angles, 'y':energies, 'marker':'b.'}, {'x':fit_angles, 'y':svd_energies, 'marker':'b-'}
                   , {'x':fit_angles, 'y':cauchy_energies, 'marker':'r-'}],
                   '/Users/iwelsh/Desktop/testthing.png', xlim=(0,360))
    
    
    
    print(cauchy_fit)
    sys.exit()
    
    start_time = process_time()
    for r in roots:
        for m in mols:
            for t in types:
                #for l in levls:
                file_path = "/Users/iwelsh/GitHub/ExtractedData/Torsions/AdditionalFixedTorsions/{}/{}{}/".format(r,m,t)
                if not os.path.isfile(file_path + 'energies.aa.txt'):
                    continue
                print('{}{}'.format(m,t))
                with open('/Users/iwelsh/GitHub/ExtractedData/Torsions/differenceMethods_{}{}.txt'.format(m,t), 'w') as fh:
                    main(fh, file_path)
                
#                    for j in range(len(returned_rmsds)):
#                        mean_rmsds[j] += returned_rmsds[j]
#                    if not (i) % 10:
#                        group_time = process_time()
#                        print("On iteration: ", i+1, ", time for last 10 = ", int(group_time - start_time))
#                        start_time = group_time
#        fh.write("Final mean rmsds after {} iterations:\n".format(i))
#        for j in mean_rmsds:
#            fh.write("{}\n".format(j/i))
#        
    