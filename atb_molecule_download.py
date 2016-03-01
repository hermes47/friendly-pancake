'''
Created on 18/02/2016

@author: iwelsh
'''
import atb_api



def main(api, molid, aa_mtb_out=None, ua_mtb_out=None, pdb_out=None):
    mol = api.Molecules.molid(molid=molid)
    #if aa_mtb_out is not None:
    #    mol.download_file(fnme=aa_mtb_out, format='mtb_aa')
    #if ua_mtb_out is not None:
    #    mol.download_file(fnme=ua_mtb_out, format='mtb_ua')
    if pdb_out is not None:
        mol.download_file(fnme=pdb_out, format='pdb')

if __name__ == '__main__':
    api = atb_api.API(debug=True)
    root = '/Users/iwelsh/GitHub/ExtractedData/Torsions/Parameters'
    mols = {'AMINO0':26640,
            'AMINO2':26642,
            'CHLORO0':26643,
            'CHLORO2':26645,
            'HYDRO0':847,
            'HYDRO2':26647,
            'METH-1':756,
            'METH0':756,
            'METH2':26649,
            'THIO0':26650,
            'THIO2':26652,
            }
    for m in mols:
        ua_out = '{0}/{1}.ua.mtb'.format(root,m)
        aa_out = '{0}/{1}.aa.mtb'.format(root,m)
        pdb_out = '{0}/{1}.pdb'.format(root,m)
        #main(api, mols[m], aa_mtb_out=aa_out, ua_mtb_out=ua_out, pdb_out=pdb_out)
        break
    vals = [[4.18, 4.14, 3.99, 4.81, 4.96],
            [5.51, 4.62, 5.12, 4.81, 5.77],
            [4.04, 3.90, 3.57, 4.35, 4.37],
            [5.02, 4.77, 4.94, 4.66, 4.85],
            [5.49, 5.38, 5.35, 6.05, 5.93],
            [6.23, 6.83, 5.90, 6.05, 6.73],
            [5.36, 5.46, 5.13, 6.02, 5.78],
            [6.35, 6.22, 6.36, 6.26, 6.21]]
    import numpy as np
    for val in vals:
        print(np.std(val))