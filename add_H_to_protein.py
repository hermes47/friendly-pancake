'''
Created on 9/03/2016

@author: iwelsh
'''
from mrsuper.lib.genlib import load_file
from mrsuper.lib.inlib import gromos_cnf_parse, gromos_topology_parse, pdb_parse
from mrsuper.lib.outlib import print_pdb
from mrsuper.lib.storage import Molecule
from mrsuper.lib.structure import convert_to_aa, neighbor_h_count, unite_atoms, united_h_count

def main():
    root = '/Users/iwelsh/GitHub/ShortScripts/proteinStructures'
    coords = gromos_cnf_parse(load_file(root + '/1d66_ss-dna_54a8.gch'))[1]
    #pdb = pdb_parse(load_file(root + '/3blg.pdb'))
    topo = gromos_topology_parse(load_file(root + '/1d66_ss-dna_54a8.top'))
    mol = Molecule('3BLG', cnf=coords, topo=topo)
    print(len(mol.atoms))
    #for atm in mol.atoms:
    #    print(atm.data['residue_name'])
    
    #unite_atoms(mol)
    conect_graph = mol.labeled_graph
    for atm in mol.atoms:
        atm.hcount = united_h_count(atm)
    with open(root + '/1d66_gch.pdb','w') as fh:
        fh.write(print_pdb(mol))
    with open(root +'/1d66_aa.pdb','w') as fh:
        fh.write(print_pdb(convert_to_aa(mol)))
    
if __name__ == '__main__':
    main()