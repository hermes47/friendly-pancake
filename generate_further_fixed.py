'''
Created on 8/01/2016

@author: iwelsh
'''
import os

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

def gen_files(mol, mol_data, descriptor, descriptor_data, root='/Users/iwelsh/Documents/AdditionRigidTorsions'):
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

def main():
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
    
    for m in molecules:
        for d in descriptors:
            formatting = dict(mol=m, des=d, root='/Users/iwelsh/Documents/AdditionRigidTorsions')
            if not os.path.isdir('{root}/{des}/{mol}'.format(**formatting)):
                os.makedirs('{root}/{des}/{mol}'.format(**formatting))
            if descriptors[d] is None:
                continue
            gen_files(m, molecules[m], d, descriptors[d])
                

if __name__ == '__main__':
    main()