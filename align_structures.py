'''
Created on 10/03/2016

@author: iwelsh
'''
from vector_math_functions import kabsch_alignment
import numpy as np
from mrsuper.lib.storage import Molecule
from mrsuper.lib.genlib import load_file
from mrsuper.lib.inlib import gromos_cnf_parse, gromos_topology_parse, gamess_log_parse
from mrsuper.lib.outlib import print_pdb
import os

def main():
    root = '/Users/iwelsh/GitHub/ExtractedData/Torsions/AdditionalFixedTorsions'
    for mol, data in mol_data():
        if not os.path.isdir(root + '/AlignedStructures/' + mol + '/EM_Changes'):
            os.makedirs(root + '/AlignedStructures/' + mol + '/EM_Changes')
        formatting = dict(root=root+'/Original', mol=mol, root2=root)
        topo = gromos_topology_parse(load_file('{root}/{mol}/MD_data/{mol}.ua.top'.format(**formatting)))
        cnf = gromos_cnf_parse(load_file('{root}/{mol}/MD_data/{mol}.000.ua'
                                             '.min.1.cnf'.format(**formatting)))
        cnf_initial = gromos_cnf_parse(load_file('{root}/{mol}/MD_data/{mol}.000.ua.cnf'.format(**formatting)))
        mol_ob = Molecule(mol, cnf=cnf[1], topo=topo)
        mol_ob_initial = Molecule(mol, cnf=cnf_initial[1], topo=topo)
        previous, previous_initial = [], []
        mol_ob.centre_molecule(data['ua'][1])
        mol_ob_initial.centre_molecule(data['ua'][1])
        for atm in data['ua'][:3]:
            previous.append(np.asarray(mol_ob.atom(atm).xyz))
            previous_initial.append(np.asarray(mol_ob_initial.atom(atm).xyz))
        previous = np.asarray(previous)
        previous_initial = np.asarray(previous_initial)
        
        for level, atoms in data.items():
            formatting = dict(root=root+'/Original', level=level, mol=mol, root2=root)
            if level in ['gamess','aa']:
                topo = gromos_topology_parse(load_file('{root}/{mol}/MD_data/{mol}.aa.top'.format(**formatting)))
            else:
                topo = gromos_topology_parse(load_file('{root}/{mol}/MD_data/{mol}.ua.top'.format(**formatting)))
            for i in range(360):
                formatting['i'] = i
                try:
                    if level == 'gamess':
                        gamess = gamess_log_parse(load_file('{root}/{mol}/{mol}.{i:03}.log'.format(**formatting),
                                                             blank_lines=True))
                        mol_ob2 = Molecule(mol, gamess=gamess[max(gamess)], topo=topo)
                        mol_ob2_initial = Molecule(mol, gamess=gamess[min(gamess)], topo=topo)
                    else:
                        cnf = gromos_cnf_parse(load_file('{root}/{mol}/MD_data/{mol}.{i:03}.{level}'
                                                         '.min.1.cnf'.format(**formatting)))
                        cnf_initial = gromos_cnf_parse(load_file('{root}/{mol}/MD_data/{mol}.{i:03}.{level}'
                                                         '.cnf'.format(**formatting)))
                        mol_ob2 = Molecule(mol, cnf=cnf[1], topo=topo)
                        mol_ob2_initial = Molecule(mol, cnf=cnf_initial[1], topo=topo)
                except FileNotFoundError:
                    continue
                    
                current, current_initial = [], []
                if level == 'gamess':
                    mol_ob2.centre_molecule(data['aa'][1])
                    mol_ob2_initial.centre_molecule(data['aa'][1])
                    atms = data['aa']
                else:
                    try:
                        mol_ob2.centre_molecule(atoms[1])
                    except TypeError:
                        print(mol, level, i, atoms, len(mol_ob2.atoms))
                        print(mol_ob2.atom(atoms[1]).xyz)
                        for atm in mol_ob2.atoms:
                            print(atm.data)
                        
                        raise
                    mol_ob2_initial.centre_molecule(atoms[1])
                    atms = atoms
                
                for atm in atms[:3]:
                    current.append(np.asarray(mol_ob2.atom(atm).xyz))
                    current_initial.append(np.asarray(mol_ob2_initial.atom(atm).xyz))
                current = np.asarray(current)
                current_initial = np.asarray(current_initial)
                rot_mat = kabsch_alignment(current, previous)
                rot_mat_initial = kabsch_alignment(current_initial, previous)
                mol_ob2.rotate(rot_mat)
                mol_ob2_initial.rotate(rot_mat_initial)
                with open(('{root2}/AlignedStructures/{mol}/EM_Changes/'
                          '{mol}.{i:03}.{level}.initial.pdb'.format(**formatting)), 'w') as fh:
                    fh.write(print_pdb(mol_ob2_initial))
                with open('{root2}/AlignedStructures/{mol}/{mol}.{i:03}.{level}.pdb'.format(**formatting), 'w') as fh:
                    fh.write(print_pdb(mol_ob2))
                
                    
            
        
def mol_data():
    mols = {'AMINO-1'  : {'gamess': [2, 1, 8, 9], 'aa': [16, 13, 7, 9], 'ua': [8, 7, 5, 6]},
            'AMINO0'   : {'gamess': [2, 1, 8, 9], 'aa': [2, 5, 8, 13],  'ua': [1, 2, 3, 8]},
            'AMINO1'   : {'gamess': [2, 1, 8, 9], 'aa': [2, 5, 8, 16],  'ua': [1, 2, 3, 8]},
            'AMINO2'   : {'gamess': [2, 1, 8, 9], 'aa': [7, 10, 12, 15],'ua': [5, 6, 7, 8]},
            'CHLORO-1' : {'gamess': [2, 1, 8, 9], 'aa': [2, 5, 8, 14],  'ua': [1, 2, 3, 6]},
            'CHLORO0'  : {'gamess': [2, 1, 8, 9], 'aa': [2, 5, 8, 11],  'ua': [1, 2, 3, 5]},
            'CHLORO1'  : {'gamess': [2, 1, 8, 9], 'aa': [2, 5, 8, 14],  'ua': [1, 2, 3, 6]},
            'CHLORO2'  : {'gamess': [2, 1, 8, 9], 'aa': [10, 8, 17, 20],'ua': [5, 4, 7, 8]},
            'HYDRO-1'  : {'gamess': [2, 1, 8, 9], 'aa': [15, 12, 6, 8], 'ua': [7, 6, 4, 5]},
            'HYDRO0'   : {'gamess': [2, 1, 8, 9], 'aa': [2, 5, 8, 12],  'ua': [1, 2, 3, 7]},
            'HYDRO1'   : {'gamess': [2, 1, 8, 9], 'aa': [2, 5, 8, 15],  'ua': [1, 2, 3, 7]},
            'HYDRO2'   : {'gamess': [2, 1, 8, 9], 'aa': [2, 5, 8, 10],  'ua': [1, 2, 3, 4]},
            'METH-1'   : {'gamess': [2, 1, 8, 9], 'aa': [2, 5, 8, 14],  'ua': [5, 3, 2, 1]},
            'METH0'    : {'gamess': [2, 1, 8, 9], 'aa': [2, 5, 8, 14],  'ua': [5, 3, 2, 1]},
            'METH1'    : {'gamess': [2, 1, 8, 9], 'aa': [2, 5, 8, 10],  'ua': [1, 2, 3, 4]},
            'METH2'    : {'gamess': [2, 1, 8, 9], 'aa': [2, 5, 8, 10],  'ua': [1, 2, 3, 4]},
            'THIO-1'   : {'gamess': [2, 1, 8, 9], 'aa': [2, 5, 8, 15],  'ua': [1, 2, 3, 7]},
            'THIO0'    : {'gamess': [2, 1, 8, 9], 'aa': [2, 5, 8, 12],  'ua': [1, 2, 3, 7]},
            'THIO1'    : {'gamess': [2, 1, 8, 9], 'aa': [2, 5, 8, 15],  'ua': [1, 2, 3, 7]},
            'THIO2'    : {'gamess': [2, 1, 8, 9], 'aa': [6, 9, 18, 21], 'ua': [4, 5, 8, 9]},                                   
            }
    for mol in mols:
        yield mol, mols[mol]



if __name__ == '__main__':
    main()