'''
Created on 16/12/2015

@author: iwelsh
'''
import shutil

import os

def main():
    root_dir = '/data/welsh/Projects/Headgroup_Parameterisation/PE_Headgroup_Parameterisation'
    sub_template = 'template.sh'
    lipids = ['DPPE','POPE']
    charge_sets = ['Q1','Q2','Q3','Q4','Q5']
    lj_sets = ['LJ1','LJ2','LJ3','LJ4','LJ5','LJ6']
    electro_types = ['TR_RF','SR_RF','PME']
    temps = {'DPPE':338,'POPE':308}
    ff_files = ['{lipid}.itp','ff_dum.itp','ffG54a7bon.itp','ffG54a7.itp','ffG54a7nb.itp','ions.itp','spc.itp']
    ff_files = ['{q}_{lj}/'+x for x in ff_files] + ['{lipid}.ndx','{lipid}.top']
    for lipid in lipids:
        for charge in charge_sets:
            for lj in lj_sets:
                for electro in electro_types:
                    formatting = {'lipid':lipid,
                                  'electro':electro,
                                  'temp':temps[lipid],
                                  'q':charge,
                                  'lj':lj,
                                  'root':root_dir,
                                  'template':sub_template
                                  }
                    
                    # load files
                    with open('{root}/{template}'.format(**formatting), 'r') as fh:
                        template_file = fh.readlines()
                    with open('{root}/JobFiles/{electro}/job1_em_solv.mdp'.format(**formatting), 'r') as fh:
                        job1a = fh.readlines()
                    with open('{root}/JobFiles/{electro}/job1b_em_solv.mdp'.format(**formatting), 'r') as fh:
                        job1b = fh.readlines()
                    with open('{root}/JobFiles/{electro}/job2_anneal.mdp'.format(**formatting), 'r') as fh:
                        job2 = fh.readlines()
                    with open('{root}/JobFiles/{electro}/job3_npt_simul.mdp'.format(**formatting), 'r') as fh:
                        job3 = fh.readlines()
                    
                    if not os.path.isdir('{root}/{electro}/{q}/{lj}/{lipid}'.format(**formatting)):
                        os.makedirs('{root}/{electro}/{q}/{lj}/{lipid}'.format(**formatting))
                    
                    # write files
                    with open('{root}/{electro}/{q}/{lj}/{lipid}/job1_em_solv.mdp'.format(**formatting),'w') as fh:
                        fh.write(''.join((x.format(**formatting) for x in job1a)))
                    with open('{root}/{electro}/{q}/{lj}/{lipid}/job1b_em_solv.mdp'.format(**formatting),'w') as fh:
                        fh.write(''.join((x.format(**formatting) for x in job1b)))
                    with open('{root}/{electro}/{q}/{lj}/{lipid}/job2_anneal.mdp'.format(**formatting),'w') as fh:
                        fh.write(''.join((x.format(**formatting) for x in job2)))
                    with open('{root}/{electro}/{q}/{lj}/{lipid}/job3_npt_simul.mdp'.format(**formatting),'w') as fh:
                        fh.write(''.join((x.format(**formatting) for x in job3)))
                    with open('{root}/{electro}/{q}/{lj}/{lipid}/sub.sh'.format(**formatting),'w') as fh:
                        fh.write(''.join((x.format(**formatting) for x in template_file)))
                        
                    # copy files
                    for file in ff_files:
                        shutil.copy('{root}/ForceFields/'.format(**formatting)+file.format(**formatting), 
                                    '{root}/{electro}/{q}/{lj}/{lipid}/'.format(**formatting))
                    shutil.copy('{root}/Coordinates/{lipid}.gro'.format(**formatting),
                                '{root}/{electro}/{q}/{lj}/{lipid}/'.format(**formatting))
                    
                    #os.chdir('{root}/{electro}/{q}/{lj}/{lipid}/'.format(**formatting))
                    #os.system('sbatch sub.sh')
                        
                    
                    
                    

if __name__ == '__main__':
    main()