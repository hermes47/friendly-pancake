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

def perform_fit(matrix, target, fit_phase, fit_method='svd'):
    if fit_method == 'qr':
        q,r = np.linalg.qr(matrix)
        d = np.dot(np.transpose(q),target)
        x = np.linalg.solve(r, d)
    elif fit_method == 'svd':
        x = np.linalg.lstsq(matrix, target)[0]
        
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

def fit_torsion_terms_lls(energies, angles, fit_phase=False, allowed_multiplicites=list(range(1,7)), baseline_shift=False):
    
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
    
    for x in energies:
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
    
    # V_{i} - c^{V}, c^{V} = mean(sum of V_i)
    if baseline_shift:
        target_vector -= np.mean(target_vector)
    else:
        target_vector -= np.mean(target_vector)
        
    # subtract the mean of cos(mi\thetaj) from each value of coefficient_matrix
    for col in range(len(coefficient_matrix[0])):
        col_vals = [];
        for row in range(len(coefficient_matrix)):
            col_vals.append(coefficient_matrix[row][col])
        mean_col = np.mean(col_vals)
        for row in range(len(coefficient_matrix)):
            coefficient_matrix[row][col] -= mean_col
    
    fit = perform_fit(coefficient_matrix, target_vector, fit_phase)
    
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
            
    
    

def one_pass_lls(energies, angles, limit=None, phase=False):
    if phase:
        assert isinstance(phase, (list,tuple)) and len(phase) == 2
    fit = fit_torsion_terms_lls(energies, angles, fit_phase=phase)
    
    if limit is not None:
        sorted_list = sorted(fit, key=lambda x:fit[x]['k'], reverse=True)[:limit]
        for k in set(fit.keys()).difference(sorted_list):
            del fit[k]
    return fit
    
def prune_data(energies, angles):
    # prune the data to remove outliers.
    # perform a phased, one_pass_lls fit using all points
    # calculate the variance between each point and it's fitted value to find the most variant points
    # if the most variant point is above some threshold, remove it and refit
    # continue until most variant point is below threshold
    start_len = len(energies)
    while True:
        fit = one_pass_lls(energies, angles, phase=[0,90])
        current_variance = 0
        current_variant_point = -1
        total_variance = []
        for ang in sorted(angles):
            variance = np.abs(energy_at_x(fit, angles[ang]) - energies[ang])
            total_variance.append(variance)
            if variance > current_variance:
                current_variance = variance
                current_variant_point = ang
        mean_variance = np.mean(variance)
        if current_variance < 1:#*mean_variance:
            print('removed: ', start_len - len(energies))
            #print('closing variance: ', current_variant_point, current_variance,mean_variance)
            break
        
        del energies[current_variant_point], angles[current_variant_point]
        #print('removed data at point ', current_variant_point, current_variance,mean_variance)
    
    return energies, angles
    

def fit_torsion_terms_fft(energies, *angles):
    pass

def rawminusraw(energiesA, energiesB, anglesA, limit, phase, fh):
    diffEnergies, diffAngles = {}, {}
    for i in energiesA:
        if i in energiesB:
            diffEnergies[i] = energiesA[i] - energiesB[i]
            diffAngles[i] = anglesA[i]
    fit = one_pass_lls(diffEnergies, diffAngles, limit, phase)
    for term in fit:
        print('k{m:<2} = {k:11.7f}, delta{m:<2} = {delta:11.7f}'.format(**fit[term]), file=fh)
    print("RMSD fit = {:11.7f} \n".format(fit_deviation(diffEnergies, diffAngles, fit)), file=fh)
    
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
        print("\nQM Fit:", file=fh)
        for term in fit:
            print('k{m:<2} = {k:11.7f}, delta{m:<2} = {delta:11.7f}'.format(**fit[term]), file=fh)
        print("RMSD fit = {:11.7f} \n".format(fit_deviation(energiesB, anglesB, fit)), file=fh)
    print("\nDIFF Fit:", file=fh)
    for term in difffit:
        print('k{m:<2} = {k:11.7f}, delta{m:<2} = {delta:11.7f}'.format(**difffit[term]), file=fh)
    print("RMSD fit = {:11.7f} \n".format(fit_deviation(energiesDIFF, anglesDIFF, difffit)), file=fh)
    
    
    
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
    rawminusraw(energiesQM, energiesAA, anglesQM, 3, False, f)
    print("Non-phased fitting to raw QM - raw UA:", file=f)
    rawminusraw(energiesQM, energiesUA, anglesQM, 3, False, f)
    print("Phased fitting to raw QM - raw AA:", file=f)
    rawminusraw(energiesQM, energiesAA, anglesQM, 3, [0,90], f)
    print("Phased fitting to raw QM - raw UA:", file=f)
    rawminusraw(energiesQM, energiesUA, anglesQM, 3, [0,90], f)
    print("===========================================", file=f)
    
    print("Non-phased fitting to raw QM - non-phased AA:", file=f)
    fitminusraw(energiesQM, energiesAA, anglesQM, anglesAA, 'B', 3, False, False, f)
    print("Non-phased fitting to raw QM - non-phased UA:", file=f)
    fitminusraw(energiesQM, energiesUA, anglesQM, anglesUA, 'B', 3, False, False, f)
    print("Non-phased fitting to raw QM - phased AA:", file=f)
    fitminusraw(energiesQM, energiesAA, anglesQM, anglesAA, 'B', 3, [0,90], False, f)
    print("Non-phased fitting to raw QM - phased UA:", file=f)
    fitminusraw(energiesQM, energiesUA, anglesQM, anglesUA, 'B', 3, [0,90], False, f)
    print("Phased fitting to raw QM - non-phased AA:", file=f)
    fitminusraw(energiesQM, energiesAA, anglesQM, anglesAA, 'B', 3, False, [0,90], f)
    print("Phased fitting to raw QM - non-phased UA:", file=f)
    fitminusraw(energiesQM, energiesUA, anglesQM, anglesUA, 'B', 3, False, [0,90], f)
    print("Phased fitting to raw QM - phased AA:", file=f)
    fitminusraw(energiesQM, energiesAA, anglesQM, anglesAA, 'B', 3, [0,90], [0,90], f)
    print("Phased fitting to raw QM - phased UA:", file=f)
    fitminusraw(energiesQM, energiesUA, anglesQM, anglesUA, 'B', 3, [0,90], [0,90], f)
    print("===========================================", file=f)
    
    print("Non-phased fitting to non-phased QM - raw AA:", file=f)
    fitminusraw(energiesQM, energiesAA, anglesQM, anglesAA, 'B', 3, False, False, f)
    print("Non-phased fitting to non-phased QM - raw UA:", file=f)
    fitminusraw(energiesQM, energiesUA, anglesQM, anglesUA, 'B', 3, False, False, f)
    print("Non-phased fitting to phased QM - raw AA:", file=f)
    fitminusraw(energiesQM, energiesAA, anglesQM, anglesAA, 'B', 3, [0,90], False, f)
    print("Non-phased fitting to phased QM - raw UA:", file=f)
    fitminusraw(energiesQM, energiesUA, anglesQM, anglesUA, 'B', 3, [0,90], False, f)
    print("Phased fitting to non-phased QM - raw AA:", file=f)
    fitminusraw(energiesQM, energiesAA, anglesQM, anglesAA, 'B', 3, False, [0,90], f)
    print("Phased fitting to non-phased QM - raw UA:", file=f)
    fitminusraw(energiesQM, energiesUA, anglesQM, anglesUA, 'B', 3, False, [0,90], f)
    print("Phased fitting to phased QM - raw AA:", file=f)
    fitminusraw(energiesQM, energiesAA, anglesQM, anglesAA, 'B', 3, [0,90], [0,90], f)
    print("Phased fitting to phased QM - raw UA:", file=f)
    fitminusraw(energiesQM, energiesUA, anglesQM, anglesUA, 'B', 3, [0,90], [0,90], f)
    print("===========================================", file=f)
    
    print("Non-phased fitting to non-phased QM - non-phased AA:", file=f)
    fitminusfit(energiesQM, energiesAA, anglesQM, anglesAA, 3, False , False, f)   
    print("Non-phased fitting to non-phased QM - non-phased UA:", file=f)
    fitminusfit(energiesQM, energiesUA, anglesQM, anglesUA, 3, False , False, f)    
    print("Non-phased fitting to phased QM - phased AA:", file=f)
    fitminusfit(energiesQM, energiesAA, anglesQM, anglesAA, 3, [0,90] , False, f)
    print("Non-phased fitting to phased QM - phased UA:", file=f)
    fitminusfit(energiesQM, energiesUA, anglesQM, anglesUA, 3, [0,90] , False, f)
    print("Phased fitting to phased QM - phased AA:", file=f)    
    fitminusfit(energiesQM, energiesAA, anglesQM, anglesAA, 3, [0,90] , [0,90], f)
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

if __name__ == '__main__':
    from time import process_time
    import os
    mean_rmsds = [0,0,0,0,0,0]
    #iterations = 50000
    #roots = ['30PsiAllFixed', '60PsiBothAllFixed', 'AllHeavyFixed', 'Figures_Pruned', '30PsiBothAllFixed', 
    #         '60PsiOriginal', 'BothAllFixed', 'Figures_Raw', 'Original', 'ThetaAllFixed', '30PsiBothCarbonFixed', 
    #         'BothCarbonFixed', 'FunctionalFixed', 'ThetaCarbonFixed', '30PsiCarbonFixed', 'PsiAllFixed', '30PsiOriginal', 
    #         'AllCarbonFixed', 'PsiCarbonFixed', '60PsiAllFixed', 'AllFixed']
    roots = ['Original']
    mols = ['AMINO','CHLORO','HYDRO','METH','THIO']
    types = [-1] #[-1,0,1,2]
    levls = ['aa','ua','qm']
    
    i = 0
    
    start_time = process_time()
    for r in roots:
        for m in mols:
            for t in types:
                #for l in levls:
                file_path = "/Users/iwelsh/GitHub/ExtractedData/Torsions/AdditionalFixedTorsions/{}/{}{}/".format(r,m,t)
                if not os.path.isfile(file_path + 'energies.aa.txt'):
                    continue
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
    