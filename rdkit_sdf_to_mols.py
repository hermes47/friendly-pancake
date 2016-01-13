'''
Created on 8/12/2015

@author: iwelsh
'''

def main():
    from rdkit import Chem
    from rdkit.Chem import AllChem
    import os
    source = '/Users/iwelsh/Documents/Lipid_MOL_files/all_mols.sdf'
    suppl = Chem.SDMolSupplier(source)
    for mol in suppl:
        if mol is None:
            continue
        try:
            cat = mol.GetProp('CATEGORY')[:-5]
            cat = cat.replace(' ','_')
        except KeyError:
            cat = 'Other'
        try:
            main_class = mol.GetProp('MAIN_CLASS')[:-7]
            main_class = main_class.replace(' ','_')
        except KeyError:
            main_class = 'Other'
        final_dir = '/Users/iwelsh/Documents/Lipid_MOL_files/{}/{}'.format(cat,main_class)
        if not os.path.isdir(final_dir):
            os.makedirs(final_dir)
        name = mol.GetProp('LM_ID')
        
        
        mol2 = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol2)
        try:
            AllChem.UFFOptimizeMolecule(mol2)
        except ValueError as e:
            print(cat,main_class,name,'failed with UFF')
            try: 
                AllChem.MMFFOptimizeMolecule(mol2)
            except ValueError as e:
                print(cat,main_class,name,'failed with MMFF',e)
        
        
        Chem.MolToPDBFile(mol2, '{}/{}.pdb'.format(final_dir,name))
        Chem.MolToMolFile(Chem.RemoveHs(mol2), '{}/{}.mol'.format(final_dir,name))
        #break
        
        

if __name__ == '__main__':
    main()