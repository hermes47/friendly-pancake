'''
Created on 14/01/2016

@author: iwelsh
'''
import os
import yaml

from vector_math_functions import * #dihedral_angle, two_pass_lls, one_pass_lls, energy_at_x, prune_data, fit_deviation, derivative
from gamess_extraction import extract_conformation, extract_energy
from gromos_extraction import extract_conf, extract_ene
from plotting import plot_scatters
from mrsuper.lib import genlib, inlib, outlib, storage
from numpy import array, mean
import numpy as np

def extract_energies_qm(root, descriptors, molecules, angles_to_measure):
    for b in descriptors:
        for m in molecules:
            energies, angles = {}, {}
            shifts = angles_to_measure[molecules.index(m)]
            formatting = dict(root=root,
                              basis=b,
                              mol=m)
            for a in range(0,360,5):
                formatting['angle2'] = a
                if a == 180:
                    a = 179
                formatting['angle'] = a
                if (not os.path.isfile('{root}/{basis}/{mol}/{mol}.{angle:03}.log'.format(**formatting)) 
                    and not os.path.isfile('{root}/{basis}/{mol}/{mol}.{angle2:03}.log'.format(**formatting))):
                    continue
                try:
                    with open('{root}/{basis}/{mol}/{mol}.{angle:03}.log'.format(**formatting),'r') as fh:
                        file = fh.readlines()
                except FileNotFoundError:
                    with open('{root}/{basis}/{mol}/{mol}.{angle2:03}.log'.format(**formatting),'r') as fh:
                        file = fh.readlines()
                conf = extract_conformation(file)
                if conf is None:
                    continue                
                energies[a] = extract_energy(file)
                angles[a] = dihedral_angle(conf[shifts[0]]['vec'],conf[shifts[1]]['vec'],
                                           conf[shifts[2]]['vec'],conf[shifts[3]]['vec'], scale='deg')
            if not energies:
                print('Nothing extracted for {basis}:{mol}'.format(**formatting))
                continue
            try:
                z_energy = energies[0]
            except KeyError:
                continue
            for i in energies:
                energies[i] -= z_energy
            
            with open('{root}/{basis}/{mol}/energies.qm.txt'.format(**formatting),'w') as fh:
                yaml.dump([energies,angles],fh)
            with open('{root}/{basis}/{mol}/fitted_curve.qm.txt'.format(**formatting),'w') as fh:
                yaml.dump(two_pass_lls(energies, angles, 3, phase=False), fh)
            with open('{root}/{basis}/{mol}/phased_fitted_curve.qm.txt'.format(**formatting),'w') as fh:
                yaml.dump(two_pass_lls(energies, angles, 3, phase=[0,90]), fh)
            with open('{root}/{basis}/{mol}/onepass_fitted_curve.qm.txt'.format(**formatting),'w') as fh:
                yaml.dump(one_pass_lls(energies, angles, phase=False), fh)
            with open('{root}/{basis}/{mol}/onepass_phased_fitted_curve.qm.txt'.format(**formatting),'w') as fh:
                yaml.dump(one_pass_lls(energies, angles, phase=[0,90]), fh)
                
def extract_energies_md(root, descriptors, molecules, angles_to_measure):
    for t in ('aa','ua'):
        for d in descriptors:
            for m in molecules:
                energies, angles = {}, {}
                shifts = angles_to_measure[molecules.index(m)][t]
                rmsds = {}
                formatting = dict(root=root,
                                  des=d,
                                  mol=m,
                                  type=t,
                                  c=1)
                for a in range(0,360,5):
                    formatting['angle2'] = a
                    if a == 180:
                        a = 179
                    formatting['angle'] = a
                    if (not os.path.isfile('{root}/{des}/{mol}/MD_data/{mol}.{angle:03}.{type}.{c}.log'.format(**formatting)) 
                        and not os.path.isfile('{root}/{des}/{mol}/MD_data/{mol}.{angle2:03}.{type}.{c}.log'.format(**formatting))):
                        continue
                    try:
                        with open('{root}/{des}/{mol}/MD_data/{mol}.{angle:03}.{type}.{c}.log'.format(**formatting),'r') as fh:
                            energy_file = fh.readlines()
                    except FileNotFoundError:
                        try:
                            with open('{root}/{des}/{mol}/MD_data/{mol}.{angle2:03}.{type}.{c}.log'.format(**formatting),'r') as fh:
                                energy_file = fh.readlines()
                        except FileNotFoundError:
                            continue
                    try:
                        with open('{root}/{des}/{mol}/MD_data/{mol}.{angle:03}.{type}.min.{c}.cnf'.format(**formatting),'r') as fh:
                            coord_file = fh.readlines()
                        with open('{root}/{des}/{mol}/MD_data/{mol}.{angle:03}.{type}.cnf'.format(**formatting),'r') as fh:
                            init_coord_file = fh.readlines()
                    except FileNotFoundError:
                        try:
                            with open('{root}/{des}/{mol}/MD_data/{mol}.{angle2:03}.{type}.min.{c}.cnf'.format(**formatting),'r') as fh:
                                coord_file = fh.readlines()
                            with open('{root}/{des}/{mol}/MD_data/{mol}.{angle2:03}.{type}.cnf'.format(**formatting),'r') as fh:
                                init_coord_file = fh.readlines()
                        except FileNotFoundError:
                            continue
                    conf = extract_conf(coord_file)
                    init_conf = extract_conf(init_coord_file)
                    if conf is None:
                        continue
                    if init_conf is not None:
                        rmsds[a] = aligned_rmsd(conf, init_conf)
                    energies[a] = extract_ene(energy_file)
                    angles[a] = dihedral_angle(conf[shifts[0]]['vec'], conf[shifts[1]]['vec'],
                                               conf[shifts[2]]['vec'], conf[shifts[3]]['vec'], scale='deg')
                if not energies:
                    print('Nothing extracted for {type}:{des}:{mol}'.format(**formatting))
                    continue
                try:
                    z_energy = energies[0]
                except KeyError:
                    try:
                        z_energy = energies[5]
                    except KeyError:
                        continue
                for i in energies:
                    energies[i] -= z_energy
                with open('{root}/{des}/{mol}/energies.{type}.txt'.format(**formatting), 'w') as fh:
                    yaml.dump([energies, angles],fh)
                with open('{root}/{des}/{mol}/rmsd.{type}.txt'.format(**formatting), 'w') as fh:
                    yaml.dump(rmsds,fh)
                #with open('{root}/{des}/{mol}/fitted_curve.{type}.txt'.format(**formatting),'w') as fh:
                #    yaml.dump(two_pass_lls(energies, angles, 3, phase=False), fh)
                #with open('{root}/{des}/{mol}/phased_fitted_curve.{type}.txt'.format(**formatting),'w') as fh:
                #    yaml.dump(two_pass_lls(energies, angles, 3, phase=[0,90]), fh)
                #with open('{root}/{des}/{mol}/onepass_fitted_curve.{type}.txt'.format(**formatting),'w') as fh:
                #    yaml.dump(one_pass_lls(energies, angles, phase=False), fh)
                #with open('{root}/{des}/{mol}/onepass_phased_fitted_curve.{type}.txt'.format(**formatting),'w') as fh:
                #    yaml.dump(one_pass_lls(energies, angles, phase=[0,90]), fh)

def plot_energies_individual(root, molecules, basis_sets, loaded_data, types):
    for t in types:
        for m in molecules:
            if m not in loaded_data:
                continue
            for d in basis_sets:
                if d not in loaded_data[m]:
                    continue
                if t not in loaded_data[m][d]:
                    continue
                #if t in ["aa",'qm']:
                #    continue
                data = [{'x':loaded_data[m][d][t][1],
                         'y':loaded_data[m][d][t][0],
                         'marker':'b.'}]
                forms = dict(root=root, mol=m, basis=d, t=t)
                if t is None:
                    save_path = '{root}/Figures/{mol}/{basis}.png'.format(**forms)
                else:
                    save_path = '{root}/Figures/{mol}/{basis}.{t}.png'.format(**forms)
                #print(t)
                fit = one_pass_lls(loaded_data[m][d][t][0], loaded_data[m][d][t][1], phase=[0,90])
                cauchy_fit = one_pass_lls(loaded_data[m][d][t][0], loaded_data[m][d][t][1], phase=[0,90], method='cauchy', printme=False)
                fit_data = {'x':[x/10 for x in range(3600)],
                             'y':[energy_at_x(fit, x/10) for x in range(3600)],
                             'marker':'b-'}
                #mean_fit = mean(fit_data['y'])
                #print('mean', mean_fit)
                #for x in fit_data['y']:
                #    x -= mean_fit
                data.append(fit_data)
                data.append({'x':[x/10 for x in range(3600)],
                             'y':[energy_at_x(cauchy_fit, x/10) for x in range(3600)],
                             'marker':'r-'})
                plot_scatters(data, save_path, show_legend=False, xlim=(0,360), 
                      x_label='Dihedral angle (degrees)', y_label=r'Potential (kJ mol$^{-1}$)', 
                      title='Dihedral profile for {mol} under {t}'.format(**forms))
                if t == 'aa' and m == 'HYDRO2':
                    error_term = error_lls = 0
                    for x in loaded_data[m][d][t][1]:
                        error_term += np.log(1 + (energy_at_x(cauchy_fit, loaded_data[m][d][t][1][x]) - loaded_data[m][d][t][0][x])**2)
                        error_lls += (energy_at_x(fit, loaded_data[m][d][t][1][x]) - loaded_data[m][d][t][0][x])**2
                    print('Error:', error_term, error_lls)
                #print(t, cauchy_fit)
                
def plot_energies_all(root, basis_sets, molecules, loaded_data, types, data_set='Combined'):
    for m in molecules:
        if m not in loaded_data:
            continue
        for t in types:
            data = []
            forms = dict(root=root, mol=m, name=data_set, t=t)
            if t is None:
                save_path = '{root}/Figures/{mol}/{name}.png'.format(**forms)
            else:
                save_path = '{root}/Figures/{mol}/{name}.{t}.png'.format(**forms)
            
            for d in basis_sets:
                if d not in loaded_data[m]:
                    continue
                if t not in loaded_data[m][d]:
                    continue
                point_data = {'x':loaded_data[m][d][t][1],
                              'y':loaded_data[m][d][t][0],
                              'label':d}
                data.append(point_data)
            if not data:
                continue
                        
            plot_scatters(data, save_path, show_legend=True, xlim=(0,360), 
                      x_label='Dihedral angle (degrees)', y_label=r'Potential (kJ mol$^{-1}$)', 
                      title='Dihedral profile for {mol} under all setups'.format(**forms))

def derivatives(energies, angles, mol, root):
    save_path = '{}/Figures/Derivatives'.format(root)
    if not os.path.isdir(save_path):
        os.makedirs(save_path)
    save_path += "/{}.png".format(mol)
    derivatives, d_angles = derivative(energies, angles)
    mean = np.mean([derivatives[x] for x in derivatives])
    std = np.std([derivatives[x] for x in derivatives])
    
    # determine if there exists a 10% length with all within xsd of mean
    at_least_1_within = False
    for i in range(len(derivatives)):
        all_within_x = True
        step_count = 0
        for j in cycle_iterator(i, list(sorted(derivatives))):
            step_count += 1
            if np.abs(derivatives[j]) > mean + 1*std:
                all_within_x = False
                break
            if step_count >= 0.1*len(derivatives):
                break
        if all_within_x:
            at_least_1_within = True
            break
    if not at_least_1_within:
        print("Not a series of 10 points within 1 sd of mean")
            
        
            
    
    
    data = [{'x':d_angles, 'y':dict(zip([x for x in sorted(derivatives)],[np.abs(derivatives[x]) for x in sorted(derivatives)])), 'marker':'b.'}]
    data.append({'x':[-1,361], 'y':[mean + std, mean + std], 'marker':'r--'}) # plus 1 sd
    data.append({'x':[-1,361], 'y':[mean + 2*std, mean + 2*std], 'marker':'g--'}) # pls 2 sd
    #data.append({'x':[-1,361], 'y':[mean-std, mean-std], 'marker':'r--'}) # min 1 sd
    #data.append({'x':[-1,361], 'y':[mean-2*std, mean-2*std], 'marker':'g--'}) # min 2 sd
    #data.append({'x':[-1,361], 'y':[mean-3*std, mean-3*std], 'marker':'b--'}) # min 3 sd
    data.append({'x':[-1,361], 'y':[mean+3*std, mean+3*std], 'marker':'b--'}) # plus 3 sd
    data.append({'x':[-1,361], 'y':[mean, mean], 'marker':'k--'}) # mean
    #initial_angle = sorted(d_angles, key=lambda x:d_angles[x])[0]
    #for ang in sorted(d_angles, key=lambda x:d_angles[x])[1:]:
    #    second_d[ang] = (derivatives[ang] - derivatives[initial_angle])/(d_angles[ang] - d_angles[initial_angle])
    #    second_d_angles[ang] = (d_angles[ang] + d_angles[initial_angle]) / 2
    #    initial_angle = ang
    
    plot_scatters(data, save_path, xlim=(0,360),
                  x_label='Dihedral angle (degrees)', y_label=r'd/dx Potential (kJ mol$^{-1}$)',
                  title='Derivative of dihedral profile for {}'.format(mol))
        
def safe_prune(energies, angles, tol=1.5):
    to_remove = []
    for ang in angles:
        if not ang - tol < angles[ang] < ang + tol:
            to_remove.append(ang) 
    for ang in to_remove:
        print("Removing point at", ang)
        del angles[ang], energies[ang]
        
    return energies, angles
       
def process_data(root, basis_sets, molecules, types=['qm', 'aa', 'ua']):
    # load all data
    loaded_data = {}
    for m in molecules:
        #loaded_data = {}
        for t in types:
            formatting = dict(root=root, mol=m)
            if not os.path.isdir('{root}/Figures/{mol}/'.format(**formatting)):
                os.makedirs('{root}/Figures/{mol}/'.format(**formatting))
            for d in basis_sets:
                formatting['basis'] = d
                formatting['type'] = t
                if not os.path.isfile('{root}/{basis}/{mol}/energies.{type}.txt'.format(**formatting)):
                    continue
                with open('{root}/{basis}/{mol}/energies.{type}.txt'.format(**formatting), 'r') as fh:
                    energies, angles = yaml.load(fh)
                print(m, t)
                #if t != 'aa':
                #    continue
                #if t == 'ua' and m in ['METH0','METH-1']:
                derivatives(energies, angles, m + '-' + t, root)
                #energies, angles = safe_prune(energies, angles)
                #derivatives(energies, angles, m + '-' + t + '-Pruned', root)
                mean_energy = mean([energies[x] for x in energies])
                for x in energies:
                    energies[x] -= mean_energy
                
                if m in loaded_data:
                    if d in loaded_data[m]:
                        loaded_data[m][d][t] = energies, angles
                    else:
                        loaded_data[m][d] = {t:(energies, angles)}
                else:
                    loaded_data[m] = {d:{t:(energies, angles)}}
                    
                if t == 'qm' and d == 'Original':# and (m.endswith('-1') or m.endswith('0')):
                    sampling_density(root, m, loaded_data[m][d][t])
                #if m in loaded_data:
                #    loaded_data[m][d] = energies, angles
                #else:
                #    loaded_data[m] = {d:(energies,angles)}
    # plot individual figures
    plot_energies_individual(root, molecules, basis_sets, loaded_data, types)
    # plot all setups on one figure
    #plot_energies_all(root, basis_sets, molecules, loaded_data, types)
    # plot any selected groups of figures
    #plot_fitted_comparisons(root, molecules, basis_sets, loaded_data, types)
    
def generate_new_qm_jobs(root, mols, descrips):
    #############################
    def ifzmat_string(data, values):
        string = []
        for source in values:
            for dihed in data[source]:
                dihed_string = [3] + list(dihed[0:4])
                string.append(','.join(map(str,dihed_string)))
        return 'IFZMAT(1)=\n' + ',\n'.join(string)
    
    def fvalue_string(data, values, fixed_angle):
        if fixed_angle > 180:
            fixed_angle -= 360
        string = []
        for source in values:
            for dihed in data[source]:
                if dihed[4] is not None:
                    string.append('{:.5f}'.format(dihed[4]))
                else:
                    string.append('{:.5f}'.format(fixed_angle))
        
        return 'FVALUE(1)=\n' + ',\n'.join(string)
    
    def gen_files(mol, mol_data, descriptor, descriptor_data, root):
        forms = dict(sed='/usr/local/bin/sed',
                     root=root,
                     mol=mol,
                     des=descriptor)
        try:
            forms['shift'] = int(descriptor[0:2])
            with open('{root}/{mol}_{shift}.template'.format(**forms),'r') as fh:
                mol_template = fh.read()
        except ValueError:
            with open('{root}/{mol}.template'.format(**forms),'r') as fh:
                mol_template = fh.read()
        
        with open('{root}/ll.template'.format(**forms),'r') as fh:
            ll_template = fh.read()
        for a in range(0,360,5):
            forms['ang'] = '{:03}'.format(a)
            if a == 180:
                a = 179
            file_forms = {'ANGLE':a if a < 180 else a - 360,
                          'IFZMAT':ifzmat_string(mol_data, set(['rigid'] + descriptor_data)),
                          'FVALUE':fvalue_string(mol_data, set(['rigid'] + descriptor_data), a)}
            with open('{root}/{des}/{mol}/{mol}.{ang}.inp'.format(**forms),'w') as fh:
                fh.write(mol_template.format(**file_forms))
            with open('{root}/{des}/{mol}/{mol}.{ang}.ll'.format(**forms),'w') as fh:
                fh.write(ll_template.format(**forms))
    ###########################
    
    
    molecules = {'AMINO1':  {'rigid':[(2,1,8,9, None)],
                             'func':[(16,15,14,8,-54.8),(17,15,14,8,59.97)],
                             'thetaC':[(15,14,8,9,-65.08)],
                             'thetaH':[(20,14,8,9,173.09),(21,14,8,9,57.38)],
                             'psiC':[(10,9,8,1,-71.06)],
                             'psiH':[(18,9,8,1,48.94),(19,9,8,1,168.94)],
                             'psiC_30':[(10,9,8,1,-41.06)],
                             'psiH_30':[(18,9,8,1,78.94),(19,9,8,1,-161.06)],
                             'psiC_60':[(10,9,8,1,-11.06)],
                             'psiH_60':[(18,9,8,1,108.94),(19,9,8,1,-131.06)],
                             },
                 'CHLORO1': {'rigid':[(2,1,8,9,None)],
                             'func':[],
                             'thetaC':[(19,16,8,9,-68.01)],
                             'thetaH':[(17,16,8,9,174.00),(18,16,8,9,50.71)],
                             'psiC':[(10,9,8,1,-73.63)],
                             'psiH':[(14,9,8,1,46.37),(15,9,8,1,166.37)],
                             'psiC_30':[(10,9,8,1,-43.63)],
                             'psiH_30':[(14,9,8,1,76.37),(15,9,8,1,-163.63)],
                             'psiC_60':[(10,9,8,1,-13.63)],
                             'psiH_60':[(14,9,8,1,106.37),(15,9,8,1,-133.63)],
                             },
                 'HYDRO1':  {'rigid':[(2,1,8,9, None)],
                             'func':[(16,15,14,8,-63.29)],
                             'thetaC':[(15,14,8,9,-60.87)],
                             'thetaH':[(19,14,8,9,174.98),(20,14,8,9,56.75)],
                             'psiC':[(10,9,8,1,-72.70)],
                             'psiH':[(17,9,8,1,47.30),(18,9,8,1,167.30)],
                             'psiC_30':[(10,9,8,1,-42.70)],
                             'psiH_30':[(17,9,8,1,77.30),(18,9,8,1,-162.70)],
                             'psiC_60':[(10,9,8,1,-12.70)],
                             'psiH_60':[(17,9,8,1,107.30),(18,9,8,1,-132.70)],
                             }, 
                 'METH1':   {'rigid':[(2,1,8,9, None)],
                             'func':[(16,15,14,8,-179.36394),(17,15,14,8,-59.77456),(18,15,14,8,60.574908)],
                             'thetaC':[(15,14,8,9,-73.08)],
                             'thetaH':[(21,14,8,9,165.26),(22,14,8,9,50.02)],
                             'psiC':[(10,9,8,1,-72.73)],
                             'psiH':[(19,9,8,1,47.27),(20,9,8,1,167.27)],
                             'psiC_30':[(10,9,8,1,-42.73)],
                             'psiH_30':[(19,9,8,1,77.27),(20,9,8,1,-162.73)],
                             'psiC_60':[(10,9,8,1,-12.73)],
                             'psiH_60':[(19,9,8,1,107.27),(20,9,8,1,-132.73)],
                             }, 
                 'THIO1':   {'rigid':[(2,1,8,9, None)],
                             'func':[(16,15,14,8,-65.93)],
                             'thetaC':[(15,14,8,9,-67.66)],
                             'thetaH':[(19,14,8,9,169.38),(20,14,8,9,50.38)],
                             'psiC':[(10,9,8,1,-72.12)],
                             'psiH':[(17,9,8,1,47.88),(18,9,8,1,167.88)],
                             'psiC_30':[(10,9,8,1,-42.12)],
                             'psiH_30':[(17,9,8,1,77.88),(18,9,8,1,-162.12)],
                             'psiC_60':[(10,9,8,1,-12.12)],
                             'psiH_60':[(17,9,8,1,107.88),(18,9,8,1,-132.12)],
                             },
                 'AMINO-1':  {'rigid':[(2,1,8,9, None)],
                             },
                 'CHLORO-1':  {'rigid':[(2,1,8,9, None)],
                             },
                 'HYDRO-1':  {'rigid':[(2,1,8,9, None)],
                             },
                 'METH-1':  {'rigid':[(2,1,8,9, None)],
                             },
                 'THIO-1':  {'rigid':[(2,1,8,9, None)],
                             },
                 }
    
    descriptors = {'Original':          ['rigid'],
                   'FunctionalFixed':   ['func'], 
                   'ThetaCarbonFixed':  ['thetaC'],
                   'ThetaAllFixed':     ['thetaC','thetaH'], 
                   'PsiCarbonFixed':    ['psiC'],
                   'PsiAllFixed':       ['psiC','psiH'], 
                   'BothCarbonFixed':   ['thetaC','psiC'],
                   'BothAllFixed':      ['thetaC','thetaH','psiC','psiH'],
                   'AllHeavyFixed':     ['thetaC','psiC','func'],
                   '30PsiCarbonFixed':  ['psiC_30'],
                   '30PsiAllFixed':     ['psiC_30','psiH_30'],
                   '30PsiBothCarbonFixed':['psiC_30','thetaC'],
                   '30PsiBothAllFixed':['psiC_30','thetaC','psiH_30','thetaH'],
                   #'60PsiCarbonFixed':  ['psiC_60'],
                   '60PsiAllFixed':     ['psiC_60','psiH_60'],
                   '60PsiBothAllFixed': ['psiC_60','psiH_60','thetaC','thetaH'],
                   '60PsiOriginal':     ['rigid'],
                   '30PsiOriginal':     ['rigid'],
                   #'AllFixed':         [],  
                   }
    
    for m in mols:
        if m not in molecules:
            continue
        for d in descrips:
            if d not in descriptors:
                continue
            formatting = dict(mol=m, des=d, root=root)
            if not os.path.isdir('{root}/{des}/{mol}'.format(**formatting)):
                os.makedirs('{root}/{des}/{mol}'.format(**formatting))
            if descriptors[d] is None:
                continue
            gen_files(m, molecules[m], d, descriptors[d], root)

def sampling_density(root, molecule, loaded_data):
    if not os.path.isdir('{}/SamplingDensity'.format(root)):
        os.mkdir('{}/SamplingDensity'.format(root))
    energies, angles = loaded_data
    sample_sets = {"05 deg": list(range(0,360,5)),  
                   "10 deg": list(range(0,360,10)), 
                   "15 deg": list(range(0,360,15)), 
                   "20 deg": list(range(0,360,20)), 
                   "30 deg": list(range(0,360,30)), 
                   "expected peaks": [355,0,5,55,60,65,115,120,125,175,180,185,235,240,245,295,300,305],
                   "expected peaks with midpoints": [355,0,5,55,60,65,115,120,125,175,180,185,235,240,
                                                     245,295,300,305,30,90,150,210,270,330],
                   "expected peaks with midpoints half": [355,5,55,65,115,125,175,185,235,
                                                     245,295,305,30,90,150,210,270,330],
                   "expected troughs with midpoints":[0,55,60,65,120,175,180,185,240,295,300,305,30,90,150,210,270,330]
                   }
    formatting = dict(mol=molecule, root=root)
    with open('{root}/SamplingDensity/{mol}.qm.txt'.format(**formatting), 'w') as fh:
        for sample in sorted(sample_sets):
            sampleEnergies, sampleAngles = {}, {}
            for x in sample_sets[sample]:
                if x == 180:
                    try:
                        sampleEnergies[x-1] = energies[x-1]
                        sampleAngles[x-1] = angles[x-1]
                    except KeyError:
                        pass 
                try:
                    sampleEnergies[x] = energies[x]
                    sampleAngles[x] = angles[x]
                except KeyError:
                    pass#print('Missing {}: {}'.format(molecule, x))
            sampleFit = one_pass_lls(sampleEnergies, sampleAngles,  phase=[0,90])
            sampleRMSD = fit_deviation(energies, angles, sampleFit)
            fh.write('{}: {:.4f}\n'.format(sample, sampleRMSD))
            

def plot_fitted_comparisons(root, molecules, descriptors, loaded_data, types):
    for m in molecules:
        if m not in loaded_data:
            continue
        for d in descriptors:
            if d not in loaded_data[m]:
                continue
            if 'qm' not in loaded_data[m][d]:
                continue
            qm_fit = one_pass_lls(loaded_data[m][d]['qm'][0], loaded_data[m][d]['qm'][1], phase=[0,90], method='cauchy')
            for t in types:
                if t not in loaded_data[m][d]:
                    continue
                if t == 'qm':
                    continue
                forms = dict(root=root, mol=m, des=d, type=t, Type=t.upper())
                save_path = '{root}/Figures/{mol}/{des}_Fitting.{type}.png'.format(**forms)
                plotting_data = [{'x':loaded_data[m][d]['qm'][1],
                                  'y':loaded_data[m][d]['qm'][0],
                                  'marker':'b.',
                                  'label':'QM raw'}]
                
                md_fit = one_pass_lls(loaded_data[m][d][t][0], loaded_data[m][d][t][1], phase=[0,90], method='cauchy')
                qm_fit_energies = [energy_at_x(qm_fit, x/10) for x in range(3600)]
                md_fit_energies = [energy_at_x(md_fit, x/10) for x in range(3600)]
                
                diff_energies = list(array(qm_fit_energies) - array(md_fit_energies))
                diff_fit = one_pass_lls(dict(zip(range(3600),diff_energies)),
                                        dict(zip(range(3600),[x/10 for x in range(3600)])), limit=3, phase=False, method='cauchy')
                diff_fit_energies = [energy_at_x(diff_fit, x/10) for x in range(3600)]
                
                plotting_data.append({'x':[x/10 for x in range(3600)],
                                      'y':qm_fit_energies,
                                      'marker':'b-',
                                      'label':'QM fit'})
                plotting_data.append({'x':[x/10 for x in range(3600)],
                                      'y':md_fit_energies,
                                      'marker':'r-',
                                      'label':t.upper() + ' fit'})
                plotting_data.append({'x':[x/10 for x in range(3600)],
                                      'y':diff_fit_energies,
                                      'marker':'g-',
                                      'label':'diff fit'})
                plotting_data.append({'x':[x/10 for x in range(3600)],
                                      'y':list(array(diff_fit_energies)+array(md_fit_energies)),
                                      'marker':'k--',
                                      'label':'MD + diff'})
                
                plotting_data.append({'x':loaded_data[m][d][t][1],
                                      'y':loaded_data[m][d][t][0],
                                      'marker':'r.',
                                      'label':t.upper() + ' raw'})
                plot_scatters(plotting_data, save_path, show_legend=False, xlim=(0,360),
                          x_label='Dihedral angle (degrees)', y_label=r'Potential (kJ mol$^{-1}$)', 
                          title='Fitting to dihedral profile for {mol}:{des} with {Type}'.format(**forms))
                with open('{root}/Figures/{mol}/{des}_diff_fit.{type}.txt'.format(**forms), 'w') as fh:
                    yaml.dump(diff_fit,fh)
            

def run_md_jobs(root, mols, descrips):
    molecules = {
                 'AMINO-1': {'rigid':{'aa':[16,13,7,9],
                                     'ua':[8,7,5,6],
                                     'gamess':[2,1,8,9],
                                     },
                            'name':'G564',
                            },
                 'AMINO0': {'rigid':{'aa':[2,5,8,13],
                                     'ua':[1,2,3,8],
                                     'gamess':[2,1,8,9],
                                     },
                            'name':'N0GW',
                            },
                 'AMINO1': {'rigid':{'aa':[2,5,8,16],
                                     'ua':[1,2,3,8],
                                     'gamess':[2,1,8,9],
                                     },
                            'name':'SHJ0',
                            },
                 'AMINO2': {'rigid':{'aa':[7,10,12,15],
                                     'ua':[5,6,7,8],
                                     'gamess':[2,1,8,9],
                                     },
                            'name':'6AT1',
                            },
                 'CHLORO-1': {'rigid':{'aa':[2,5,8,14],
                                     'ua':[1,2,3,6],
                                     'gamess':[2,1,8,9],
                                     },
                            'name':'PCRN',
                            },
                 'CHLORO0': {'rigid':{'aa':[2,5,8,11],
                                     'ua':[1,2,3,5],
                                     'gamess':[2,1,8,9],
                                     },
                            'name':'682F',
                            },
                 'CHLORO1': {'rigid':{'aa':[2,5,8,14],
                                     'ua':[1,2,3,6],
                                     'gamess':[2,1,8,9],
                                     },
                            'name':'ZB0H',
                            },
                 'CHLORO2': {'rigid':{'aa':[10,8,17,20],
                                     'ua':[5,4,7,8],
                                     'gamess':[2,1,8,9],
                                     },
                            'name':'35MG',
                            },
                 'HYDRO-1': {'rigid':{'aa':[15,12,6,8],
                                     'ua':[7,6,4,5],
                                     'gamess':[2,1,8,9],
                                     },
                            'name':'H5JT',
                            },
                 'HYDRO0': {'rigid':{'aa':[2,5,8,12],
                                     'ua':[1,2,3,7],
                                     'gamess':[2,1,8,9],
                                     },
                            'name':'G096',
                            },
                 'HYDRO1': {'rigid':{'aa':[2,5,8,15],
                                     'ua':[1,2,3,7],
                                     'gamess':[2,1,8,9],
                                     },
                            'name':'W6KL',
                            },
                 'HYDRO2': {'rigid':{'aa':[2,5,8,10],
                                     'ua':[1,2,3,4],
                                     'gamess':[2,1,8,9],
                                     },
                            'name':'5KRK',
                            },
                 'METH-1': {'rigid':{'aa':[2,5,8,14],
                                     'ua':[1,2,3,5],
                                     'gamess':[2,1,8,9],
                                     },
                            'name':'G011',
                            },
                 'METH0': {'rigid':{'aa':[2,5,8,14],
                                     'ua':[1,2,3,5],
                                     'gamess':[2,1,8,9],
                                     },
                            'name':'G011',
                            },
                 'METH1': {'rigid':{'aa':[2,5,8,10],
                                     'ua':[1,2,3,4],
                                     'gamess':[2,1,8,9],
                                     },
                            'name':'4MUW',
                            },
                 'METH2': {'rigid':{'aa':[2,5,8,10],
                                     'ua':[1,2,3,4],
                                     'gamess':[2,1,8,9],
                                     },
                            'name':'YE4E',
                            },
                 'THIO-1': {'rigid':{'aa':[2,5,8,15],
                                     'ua':[1,2,3,7],
                                     'gamess':[2,1,8,9],
                                     },
                            'name':'ZK6W',
                            },
                 'THIO0': {'rigid':{'aa':[2,5,8,12],
                                     'ua':[1,2,3,7],
                                     'gamess':[2,1,8,9],
                                     },
                            'name':'1QIU',
                            },
                 'THIO1': {'rigid':{'aa':[2,5,8,15],
                                     'ua':[1,2,3,7],
                                     'gamess':[2,1,8,9],
                                     },
                            'name':'5G2Q',
                            },
                 'THIO2': {'rigid':{'aa':[6,9,18,21],
                                     'ua':[4,5,8,9],
                                     'gamess':[2,1,8,9],
                                     },
                            'name':'PJLY',
                            },
                 }
    descriptors = {'Original': ['rigid']}
    #'AMINO-1'  : {'gamess':[2,1,8,9], 'aa' : [16,13,7,9], 'ua' : [8,7,5,6], }
    #for m in sorted(molecules):
    #    print("'{}' : {},".format(m, molecules[m]['rigid']))
    #return
    with open(root + '/Parameters/Templates/dihedral.template','r') as fh:
        imd_temp_1 = fh.read()
    with open(root + '/Parameters/Templates/dihedral2.template','r') as fh:
        imd_temp_2 = fh.read()
    with open(root + '/Parameters/Templates/dihedral3.template','r') as fh:
        imd_temp_3 = fh.read()
    with open(root + '/Parameters/Templates/dihedral_constraint.template','r') as fh:
        constraint_temp = fh.read()
    
    for m in mols:
        
        forms = dict(mol=m, root=root)
        aa_mtb = inlib.gromos_mtb_parse(genlib.load_file('{root}/Parameters/{mol}/{mol}.aa.mtb'.format(**forms)))
        ua_mtb = inlib.gromos_mtb_parse(genlib.load_file('{root}/Parameters/{mol}/{mol}.ua.mtb'.format(**forms)))
        
        aa_work_mol = storage.Molecule(molecules[m]['name'], mtb=aa_mtb)
        ua_work_mol = storage.Molecule(molecules[m]['name'], mtb=ua_mtb)
        
        for d in descrips:
            forms['des'] = d
            forms['outdir'] = '{root}/{des}/{mol}/MD_data'.format(**forms)
            forms['name'] = molecules[m]['name']
            if not os.path.exists('{outdir}'.format(**forms)):
                os.makedirs('{outdir}'.format(**forms))
            with open('{outdir}/{mol}.aa.1.imd'.format(**forms),'w') as fh:
                fh.write(imd_temp_1.format(**{'numberOfAtoms':aa_work_mol.atom_count}))
            with open('{outdir}/{mol}.ua.1.imd'.format(**forms),'w') as fh:
                fh.write(imd_temp_1.format(**{'numberOfAtoms':ua_work_mol.atom_count}))
            with open('{outdir}/{mol}.aa.2.imd'.format(**forms),'w') as fh:
                fh.write(imd_temp_2.format(**{'numberOfAtoms':aa_work_mol.atom_count}))
            with open('{outdir}/{mol}.ua.2.imd'.format(**forms),'w') as fh:
                fh.write(imd_temp_2.format(**{'numberOfAtoms':ua_work_mol.atom_count}))
            with open('{outdir}/{mol}.aa.3.imd'.format(**forms),'w') as fh:
                fh.write(imd_temp_3.format(**{'numberOfAtoms':aa_work_mol.atom_count}))
            with open('{outdir}/{mol}.ua.3.imd'.format(**forms),'w') as fh:
                fh.write(imd_temp_3.format(**{'numberOfAtoms':ua_work_mol.atom_count}))
                
                
            os.chdir('{root}'.format(**forms))
            
            os.system(('/usr/local/gromos++/bin/make_top @build Parameters/{mol}/{mol}.aa.mtb Parameters/54A7.mtb '
                       '@param Parameters/54A7.ifp @seq {name} @solv H2O > {outdir}/{mol}.aa.top'.format(**forms)))
            os.system(('/usr/local/gromos++/bin/make_top @build Parameters/{mol}/{mol}.ua.mtb Parameters/54A7.mtb '
                       '@param Parameters/54A7.ifp @seq {name} @solv H2O > {outdir}/{mol}.ua.top'.format(**forms)))
            
            for a in range(0,360,5):
                skip_angle = False
                forms['ang2'] = a
                if a == 180:
                    a = 179
                forms['ang'] = a
                try:
                    with open('{root}/{des}/{mol}/{mol}.{ang:03}.log'.format(**forms), 'r') as fh:
                        gamess = extract_conformation(fh.readlines(), scale=1, key='name', optimised=False)
                except FileNotFoundError:
                    with open('{root}/{des}/{mol}/{mol}.{ang2:03}.log'.format(**forms), 'r') as fh:
                        gamess = extract_conformation(fh.readlines(), scale=1, key='name', optimised=False)     
                for atm in aa_work_mol.atoms:
                    try:
                        atm.xyz = gamess[atm.atm_name.upper()]['vec']
                    except TypeError:
                        print('TypeError', atm.atm_name, d, m, a)
                        skip_angle = True
                        break
                    except KeyError:
                        print('KeyError', atm.atm_name, d, m, a)
                        raise
                if skip_angle:
                    continue
                for atm in ua_work_mol.atoms:
                    atm.xyz = gamess[atm.atm_name.upper()]['vec']
                with open('{outdir}/{mol}.{ang:03}.aa.cnf'.format(**forms), 'w') as fh:
                    fh.write(outlib.print_gromos_cnf(aa_work_mol))
                with open('{outdir}/{mol}.{ang:03}.ua.cnf'.format(**forms), 'w') as fh:
                    fh.write(outlib.print_gromos_cnf(ua_work_mol))
                #with open('{outdir}/{mol}.{ang:03}.aa.pdb'.format(**forms), 'w') as fh:
                #    fh.write(outlib.print_pdb(aa_work_mol))
                #with open('{outdir}/{mol}.{ang:03}.ua.pdb'.format(**forms), 'w') as fh:
                #    fh.write(outlib.print_pdb(ua_work_mol))
                angs_aa, angs_ua = [], []
                for des in descriptors[d]:
                    atm_a = aa_work_mol.atom(molecules[m][des]['aa'][0]).xyz
                    atm_b = aa_work_mol.atom(molecules[m][des]['aa'][1]).xyz
                    atm_c = aa_work_mol.atom(molecules[m][des]['aa'][2]).xyz
                    atm_d = aa_work_mol.atom(molecules[m][des]['aa'][3]).xyz
                    ang = dihedral_angle(atm_a,atm_b,atm_c,atm_d,scale='deg')
                    angs_aa.append('   {1}  {2}  {3}  {4}  1.0  {0:.3f}  0.0   0.0001'.format(ang, *molecules[m][des]['aa']))
                    angs_ua.append('   {1}  {2}  {3}  {4}  1.0  {0:.3f}  0.0   0.0001'.format(ang, *molecules[m][des]['ua']))
                with open('{outdir}/{mol}.{ang:03}.aa.constraints.dat'.format(**forms), 'w') as fh:
                    fh.write(constraint_temp.format(**{'constraints':'\n'.join(angs_aa)}))
                with open('{outdir}/{mol}.{ang:03}.ua.constraints.dat'.format(**forms), 'w') as fh:
                    fh.write(constraint_temp.format(**{'constraints':'\n'.join(angs_ua)}))
                # first pass
                os.system(('/opt/local/gromosXX-1773/bin/md @topo {outdir}/{mol}.aa.top @conf {outdir}/{mol}.{ang:03}.aa.cnf '
                           '@input {outdir}/{mol}.aa.1.imd @fin {outdir}/{mol}.{ang:03}.aa.min.1.cnf @develop '
                           '@dihrest {outdir}/{mol}.{ang:03}.aa.constraints.dat > {outdir}/{mol}.{ang:03}.aa.1.log'.format(**forms) ))
                
                os.system(('/opt/local/gromosXX-1773/bin/md @topo {outdir}/{mol}.ua.top @conf {outdir}/{mol}.{ang:03}.ua.cnf '
                           '@input {outdir}/{mol}.ua.1.imd @fin {outdir}/{mol}.{ang:03}.ua.min.1.cnf @develop '
                           '@dihrest {outdir}/{mol}.{ang:03}.ua.constraints.dat > {outdir}/{mol}.{ang:03}.ua.1.log'.format(**forms) ))
                ## second pass
                #os.system(('/usr/local/md++/bin/md @topo {outdir}/{mol}.aa.top @conf {outdir}/{mol}.{ang:03}.aa.min.1.cnf '
                #           '@input {outdir}/{mol}.aa.2.imd @fin {outdir}/{mol}.{ang:03}.aa.min.2.cnf @develop '
                #           '@dihrest {outdir}/{mol}.{ang:03}.aa.constraints.dat > {outdir}/{mol}.{ang:03}.aa.2.log'.format(**forms) ))
                #
                #os.system(('/usr/local/md++/bin/md @topo {outdir}/{mol}.ua.top @conf {outdir}/{mol}.{ang:03}.ua.min.1.cnf '
                #           '@input {outdir}/{mol}.ua.2.imd @fin {outdir}/{mol}.{ang:03}.ua.min.2.cnf @develop '
                #           '@dihrest {outdir}/{mol}.{ang:03}.ua.constraints.dat > {outdir}/{mol}.{ang:03}.ua.2.log'.format(**forms) ))
                #
                ## third pass
                #os.system(('/usr/local/md++/bin/md @topo {outdir}/{mol}.aa.top @conf {outdir}/{mol}.{ang:03}.aa.min.2.cnf '
                #           '@input {outdir}/{mol}.aa.3.imd @fin {outdir}/{mol}.{ang:03}.aa.min.3.cnf @develop '
                #           '@dihrest {outdir}/{mol}.{ang:03}.aa.constraints.dat > {outdir}/{mol}.{ang:03}.aa.3.log'.format(**forms) ))
                #
                #os.system(('/usr/local/md++/bin/md @topo {outdir}/{mol}.ua.top @conf {outdir}/{mol}.{ang:03}.ua.min.2.cnf '
                #           '@input {outdir}/{mol}.ua.3.imd @fin {outdir}/{mol}.{ang:03}.ua.min.3.cnf @develop '
                #           '@dihrest {outdir}/{mol}.{ang:03}.ua.constraints.dat > {outdir}/{mol}.{ang:03}.ua.3.log'.format(**forms) ))
                

def main():
    import align_structures
    root = os.path.expanduser('~') + '/GitHub/ExtractedData/Torsions/AdditionalFixedTorsions'
    
    descriptors = ['Original', 'FunctionalFixed', 
                   'ThetaCarbonFixed', 'ThetaAllFixed', 
                   'PsiCarbonFixed', 'PsiAllFixed', 
                   'BothAllFixed', 'BothCarbonFixed',
                   '30PsiCarbonFixed', '30PsiAllFixed',
                   '30PsiBothCarbonFixed', '30PsiBothAllFixed',
                   '60PsiBothAllFixed', '60PsiAllFixed',
                   '30PsiOriginal', '60PsiOriginal',
                   #'AllFixed', 'AllHeavyFixed',
                    ]                                           
    mols = {'AMINO-1' : {'gamess': [2, 1, 8, 9], 'aa': [16, 13, 7, 9], 'ua': [8, 7, 5, 6]},
            'AMINO0' : {'gamess': [2, 1, 8, 9], 'aa': [2, 5, 8, 13], 'ua': [1, 2, 3, 8]},
            'AMINO1' : {'gamess': [2, 1, 8, 9], 'aa': [2, 5, 8, 16], 'ua': [1, 2, 3, 8]},
            'AMINO2' : {'gamess': [2, 1, 8, 9], 'aa': [7, 10, 12, 15], 'ua': [5, 6, 7, 8]},
            'CHLORO-1' : {'gamess': [2, 1, 8, 9], 'aa': [2, 5, 8, 14], 'ua': [1, 2, 3, 6]},
            'CHLORO0' : {'gamess': [2, 1, 8, 9], 'aa': [2, 5, 8, 11], 'ua': [1, 2, 3, 5]},
            'CHLORO1' : {'gamess': [2, 1, 8, 9], 'aa': [2, 5, 8, 14], 'ua': [1, 2, 3, 6]},
            'CHLORO2' : {'gamess': [2, 1, 8, 9], 'aa': [10, 8, 17, 20], 'ua': [5, 4, 7, 8]},
            'HYDRO-1' : {'gamess': [2, 1, 8, 9], 'aa': [15, 12, 6, 8], 'ua': [7, 6, 4, 5]},
            'HYDRO0' : {'gamess': [2, 1, 8, 9], 'aa': [2, 5, 8, 12], 'ua': [1, 2, 3, 7]},
            'HYDRO1' : {'gamess': [2, 1, 8, 9], 'aa': [2, 5, 8, 15], 'ua': [1, 2, 3, 7]},
            'HYDRO2' : {'gamess': [2, 1, 8, 9], 'aa': [2, 5, 8, 10], 'ua': [1, 2, 3, 4]},
            'METH-1' : {'gamess': [2, 1, 8, 9], 'aa': [2, 5, 8, 14], 'ua': [5, 3, 2, 1]},
            'METH0' : {'gamess': [2, 1, 8, 9], 'aa': [2, 5, 8, 14], 'ua': [5, 3, 2, 1]},
            'METH1' : {'gamess': [2, 1, 8, 9], 'aa': [2, 5, 8, 10], 'ua': [1, 2, 3, 4]},
            'METH2' : {'gamess': [2, 1, 8, 9], 'aa': [2, 5, 8, 10], 'ua': [1, 2, 3, 4]},
            'THIO-1' : {'gamess': [2, 1, 8, 9], 'aa': [2, 5, 8, 15], 'ua': [1, 2, 3, 7]},
            'THIO0' : {'gamess': [2, 1, 8, 9], 'aa': [2, 5, 8, 12], 'ua': [1, 2, 3, 7]},
            'THIO1' : {'gamess': [2, 1, 8, 9], 'aa': [2, 5, 8, 15], 'ua': [1, 2, 3, 7]},
            'THIO2' : {'gamess': [2, 1, 8, 9], 'aa': [6, 9, 18, 21], 'ua': [4, 5, 8, 9]},                                   
            }
    all_mols = list(mols.keys())
    qm_mols = [x for x in mols if 'gamess' in mols[x]]
    qm_angles = [mols[x]['gamess'] for x in qm_mols]
    #extract_energies_qm(root, descriptors, qm_mols, qm_angles)
    #extract_energies_qm(root, ['Original'], ['AMINO-1', 'CHLORO-1', 'HYDRO-1', 'METH-1', 'THIO-1'], [[2,1,8,9],[2,1,8,9],[2,1,8,9],[2,1,8,9],[2,1,8,9]])
    
    
    #generate_new_qm_jobs(root, ['METH-1'], ['Original'])
    run_md_jobs(root, qm_mols, ['Original'])
    #run_md_jobs(root, ['AMINO1'], ['Original'])
    md_mols = [x for x in mols if 'aa' in mols[x] and 'ua' in mols[x]]
    #md_mols = ['AMINO0']
    md_angles = [mols[x] for x in md_mols]
    #extract_energies_md(root, ['Original'], md_mols, md_angles)
    extract_energies_md(root, ['Original'], md_mols, md_angles)
    #extract_energies_md(root, descriptors, md_mols, md_angles)
    #process_data(root, ['Original'], ['AMINO1', 'CHLORO1', 'HYDRO1', 'METH1', 'THIO1'])
    process_data(root, ['Original'], all_mols)
    #align_structures.main()
    #process_data(root, ['Original'], ['AMINO1'])
    #process_data(root, descriptors, all_mols)
    
if __name__ == '__main__':
    main()
