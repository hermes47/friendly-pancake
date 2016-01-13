'''
Created on 8/12/2015

@author: iwelsh
'''
import yaml
import numpy as np
import os

def dihedral_angle(a, b, c, d, period=[0,2*np.pi]):
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
    return phi

def extract_energies(root, descriptors, molecules, angles_to_measure):
    for b in descriptors:
        for m in molecules:
            energies, angles = {}, {}
            shifts = angles_to_measure[molecules.index(m)]
            formatting = dict(root=root,
                              basis=b,
                              mol=m)
            for a in range(0,360,5):
                if a == 180:
                    a = 179
                formatting['angle'] = a
                string = '      ***** EQUILIBRIUM GEOMETRY LOCATED *****\n'
                if not os.path.isfile('{root}/{basis}/{mol}/{mol}.{angle:03}.log'.format(**formatting)):
                    continue
                with open('{root}/{basis}/{mol}/{mol}.{angle:03}.log'.format(**formatting),'r') as fh:
                    file = fh.readlines()
                    if string not in file:
                        continue
                    start_line = file.index(string)
                    ene_line = file[start_line-3]
                    energies[a] = float(ene_line.split()[3])*2625.5
                    a_coord = list(map(float, file[start_line+shifts[0]].split()[2:5]))
                    b_coord = list(map(float, file[start_line+shifts[1]].split()[2:5]))
                    c_coord = list(map(float, file[start_line+shifts[2]].split()[2:5]))
                    d_coord = list(map(float, file[start_line+shifts[3]].split()[2:5]))
                    angles[a] = float(dihedral_angle(a_coord, b_coord, c_coord, d_coord) * 180/np.pi)
            try:
                z_energy = energies[0]
            except KeyError:
                continue
            for i in energies:
                energies[i] -= z_energy
            with open('{root}/{basis}/{mol}/energies.txt'.format(**formatting),'w') as fh:
                yaml.dump([energies,angles],fh)
                    
def plot_energies_individual(root, descriptors, molecules, ymin=None):
    import matplotlib.pyplot as plt
    for m in molecules:
        formatting = dict(root=root,
                          mol=m)
        if not os.path.isdir('{root}/Figures/{mol}/'.format(**formatting)):
            os.makedirs('{root}/Figures/{mol}/'.format(**formatting))
        for b in descriptors:
            formatting['basis'] = b
            if not os.path.isfile('{root}/{basis}/{mol}/energies.txt'.format(**formatting)):
                continue
            with open('{root}/{basis}/{mol}/energies.txt'.format(**formatting),'r') as fh:
                energies, angles = yaml.load(fh)
                plt.plot([angles[x] for x in sorted(angles)],[energies[x] for x in sorted(energies)], '.')
                plt.xlim((0,360))
        if ymin is not None:
                plt.ylim((ymin,0))
                plt.xlabel('Dihedral angle (degrees)')
                plt.ylabel(r'Potential (kJ mol$^{-1}$)')
                plt.title('Dihedral profile for {mol} with {basis}'.format(**formatting))
                plt.savefig('{root}/Figures/{mol}/{basis}.png'.format(**formatting), papertype='a4')
                plt.close()
                
def plot_energies_all(root, descriptors, molecules, ymin=None):
    import matplotlib.pyplot as plt
    for m in molecules:
        formatting = dict(root=root,
                          mol=m)
        if not os.path.isdir('{root}/Figures/{mol}/'.format(**formatting)):
            os.makedirs('{root}/Figures/{mol}/'.format(**formatting))
        for b in descriptors:
            if b not in ['Original', 
                   'ThetaCarbonFixed', 
                   'PsiCarbonFixed',
                   'FunctionalFixed',
                   'BothCarbonFixed', 
                    ]:
                pass
                #continue
            formatting['basis'] = b
            if not os.path.isfile('{root}/{basis}/{mol}/energies.txt'.format(**formatting)):
                continue
            with open('{root}/{basis}/{mol}/energies.txt'.format(**formatting),'r') as fh:
                energies, angles = yaml.load(fh)
                plt.plot([angles[x] for x in sorted(angles)],[energies[x] for x in sorted(energies)], '.', label=b)
        plt.xlim((0,360))
        if ymin is not None:
            plt.ylim((ymin,0))
        plt.xlabel('Dihedral angle (degrees)')
        plt.ylabel(r'Potential (kJ mol$^{-1}$)')
        plt.title('Dihedral profile for {mol}'.format(**formatting))
        plt.legend(loc=8, fontsize='x-small', numpoints=1, mode='expand',
                   bbox_to_anchor=(0, 0.02, 1, 0.02), ncol=5)
        plt.savefig('{root}/Figures/{mol}.png'.format(**formatting), papertype='a4')
        plt.close()
        
def plot_energies_selected(root, descriptors, molecules, selected=[], ymin=None):
    import matplotlib.pyplot as plt
    for m in molecules:
        formatting = dict(root=root,
                          mol=m)
        if not os.path.isdir('{root}/Figures/{mol}/'.format(**formatting)):
            os.makedirs('{root}/Figures/{mol}/'.format(**formatting))
        for b in descriptors:
            if b not in selected:
                continue
            formatting['basis'] = b
            if not os.path.isfile('{root}/{basis}/{mol}/energies.txt'.format(**formatting)):
                continue
            with open('{root}/{basis}/{mol}/energies.txt'.format(**formatting),'r') as fh:
                energies, angles = yaml.load(fh)
                plt.plot([angles[x] for x in sorted(angles)],[energies[x] for x in sorted(energies)], '.', label=b)
        plt.xlim((0,360))
        if ymin is not None:
            plt.ylim((ymin,0))
        plt.xlabel('Dihedral angle (degrees)')
        plt.ylabel(r'Potential (kJ mol$^{-1}$)')
        plt.title('Dihedral profile for {mol}'.format(**formatting))
        plt.legend(loc=8, fontsize='x-small', numpoints=1, mode='expand',
                   bbox_to_anchor=(0, 0.02, 1, 0.02), ncol=5)
        plt.savefig('{root}/Figures/{mol}.png'.format(**formatting), papertype='a4')
        plt.close()
    
    
                    
def main():
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
    molecules = ['AMINO1', 'CHLORO1', 'HYDRO1', 'METH1', 'THIO1']
    measure_angles = [[5,4,11,12],[5,4,11,12],[5,4,11,12],[5,4,11,12],[5,4,11,12]]   # the index shifts to get a,b,c,d of the dihedral to measure, on a per molecule
    root = '/Users/iwelsh/Documents/AdditionRigidTorsions'
    extract_energies(root, descriptors, molecules, measure_angles)
    plot_energies_all(root, descriptors, molecules)
    #plot_energies_individual(root, descriptors, molecules)
    #plot_energies_selected(root, descriptors, molecules, ['Original','FunctionalFixed', 'ThetaCarbonFixed',
    #                                                      'PsiCarbonFixed','BothCarbonFixed'])
    
    #plot_energies_selected(root, descriptors, molecules, ['Original','30PsiOriginal','60PsiOriginal'])
    #plot_energies_selected(root, descriptors, molecules, ['Original','30PsiBothAllFixed','60PsiBothAllFixed'])
        

if __name__ == '__main__':
    main()