'''
Created on 8/12/2015

@author: iwelsh
'''
import yaml
import numpy as np
import os
import matplotlib.pyplot as plt
            
def plot_energies_individual(root, descriptors, molecules, ymin=None):
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
                plt.plot([angles[x] for x in sorted(angles)],[energies[x] for x in sorted(energies)], '.-', label=b)
        plt.xlim((0,360))
        if ymin is not None:
            plt.ylim((ymin,0))
        plt.xlabel('Dihedral angle (degrees)')
        plt.ylabel(r'Potential (kJ mol$^{-1}$)')
        plt.title('Dihedral profile for {mol}'.format(**formatting))
        plt.legend(loc=8, fontsize='xx-small', numpoints=1, mode='expand',
                   bbox_to_anchor=(0, 0.02, 1, 0.02), ncol=5)
        plt.savefig('{root}/Figures/{mol}.png'.format(**formatting), papertype='a4')
        plt.close()
        
def plot_energies_selected(root, descriptors, molecules, suffix, selected=[], ymin=None):
    for m in molecules:
        formatting = dict(root=root,
                          mol=m,
                          suffix=suffix)
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
        plt.savefig('{root}/Figures/{mol}/{mol}_{suffix}.png'.format(**formatting), papertype='a4')
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
    root = '/Users/ares/Documents/AdditionalRigidTorsions'
    #extract_energies(root, descriptors, molecules, measure_angles)
    plot_energies_all(root, descriptors, molecules, ymin=-60)
    plot_energies_individual(root, descriptors, molecules)
    #plot_energies_selected(root, descriptors, molecules, 'onlycarbons', ['Original','FunctionalFixed', 'ThetaCarbonFixed',
    #                                                      'PsiCarbonFixed','BothCarbonFixed'])
    #plot_energies_selected(root, descriptors, molecules, 'psishift',['Original','30PsiOriginal','60PsiOriginal'])
    #plot_energies_selected(root, descriptors, molecules, 'rigidpsishift',['Original','30PsiBothAllFixed','60PsiBothAllFixed'])
        

if __name__ == '__main__':
    main()