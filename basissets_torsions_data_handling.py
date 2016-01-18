'''
Created on 14/01/2016

@author: iwelsh
'''
import os
import yaml
from gamess_extraction import extract_conformation, extract_energy
from vector_math_functions import dihedral_angle
from plotting import plot_scatters


def extract_energies(sets, molecules, root, angles_to_measure):
    for b in sets:
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
                if not os.path.isfile('{root}/{basis}/{mol}/{mol}.{angle:03}.log'.format(**formatting)):
                    continue
                with open('{root}/{basis}/{mol}/{mol}.{angle:03}.log'.format(**formatting),'r') as fh:
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
            
            with open('{root}/{basis}/{mol}/energies.txt'.format(**formatting),'w') as fh:
                yaml.dump([energies,angles],fh)
                
def process_data(root, basis_sets, molecules):
    # load all data
    loaded_data = {}
    for m in molecules:
        formatting = dict(root=root, mol=m)
        if not os.path.isdir('{root}/Figures/{mol}/'.format(**formatting)):
            os.makedirs('{root}/Figures/{mol}/'.format(**formatting))
        for b in basis_sets:
            formatting['basis'] = b
            if not os.path.isfile('{root}/{basis}/{mol}/energies.txt'.format(**formatting)):
                continue
            with open('{root}/{basis}/{mol}/energies.txt'.format(**formatting), 'r') as fh:
                energies, angles = yaml.load(fh)
            if m in loaded_data:
                loaded_data[m][b] = energies, angles
            else:
                loaded_data[m] = {b:(energies,angles)}
    # plot individual figures
    plot_energies_individual(root, molecules, basis_sets, loaded_data)
    # plot all setups on one figure
    plot_energies_all(root, basis_sets, molecules, loaded_data)
    # plot any selected groups of figures  
    
def plot_energies_individual(root, molecules, basis_sets, loaded_data):
    for m in molecules:
        if m not in loaded_data:
            continue
        for d in basis_sets:
            if d not in loaded_data[m]:
                continue
            data = [{'x':loaded_data[m][d][1],
                     'y':loaded_data[m][d][0],
                     'marker':'k-'}]
            forms = dict(root=root, mol=m, basis=d)
            save_path = '{root}/Figures/{mol}/{basis}.png'.format(**forms)
            plot_scatters(data, save_path, show_legend=False, xlim=(0,360), 
                  x_label='Dihedral angle (degrees)', y_label=r'Potential (kJ mol$^{-1}$)', 
                  title='Dihedral profile for {mol} with {basis} basis set'.format(**forms))
                
def plot_energies_all(root, basis_sets, molecules, loaded_data, data_set='Combined'):
    for m in molecules:
        if m not in loaded_data:
            continue
        data = []
        forms = dict(root=root, mol=m, name=data_set)
        save_path = '{root}/Figures/{mol}/{name}.png'.format(**forms)
        
        for b in basis_sets:
            if b not in loaded_data[m]:
                continue
            point_data = {'x':loaded_data[m][b][1],
                          'y':loaded_data[m][b][0],
                          'label':b}
            data.append(point_data)
        
        plot_scatters(data, save_path, show_legend=True, xlim=(0,360), 
                  x_label='Dihedral angle (degrees)', y_label=r'Potential (kJ mol$^{-1}$)', 
                  title='Dihedral profile for {mol} with all basis sets'.format(**forms),
                  default_point='-')

def main():
    root = os.path.expanduser('~') + '/GitHub/ExtractedData/Torsions/BasisSets'
    sets = ['631Gd','631+Gd', '6311+G2dp', '6311++G2d2p', '6311++G3df3p']
    mols = ['AMINO0', 'CHLORO0', 'HYDRO0', 'METH0', 'THIO0']
    angles_to_measure = [[2,1,8,9],[2,1,8,9],[2,1,8,9],[2,1,8,9],[2,1,8,9]]
    #extract_energies(sets, mols, root, angles_to_measure)
    process_data(root, sets, mols)
    

if __name__ == '__main__':
    main()