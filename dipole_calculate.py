'''
Created on 16/12/2015

@author: iwelsh
'''
def main():
    import numpy as np
    ac = 1.602e-19
    debye = 3.336e-30
    charges = {1:{'H':0.4,'N':-0.5,'C':0.3},
               2:{'H':0.35,'N':-0.366,'C':0.316},
               3:{'H':0.311,'N':-0.187,'C':0.254},
               4:{'H':0.248,'N':0.129,'C':0.127},
               5:{'H':0.250,'N':0.050,'C':0.200}}
    positions = {1:{'type':'H',
                    'pos':np.array((0.0, 0.943, 1.033))*1e-10},
                 2:{'type':'H',
                    'pos':np.array((0.816, -0.471, 1.033))*1e-10},
                 3:{'type':'H',
                    'pos':np.array((-0.816, -0.471, 1.033))*1e-10},
                 4:{'type':'N',
                    'pos':np.array((0.0, 0.0, 0.700))*1e-10},
                 5:{'type':'C',
                    'pos':np.array((0.0, 0.000, -0.770))*1e-10}}
    for charge_set in charges:
        tot_dipole = np.array((0,0,0))*1e-10
        for pos in positions:
            tot_dipole = tot_dipole + charges[charge_set][positions[pos]['type']]*positions[pos]['pos']*ac/debye
        print('Charge set: {}, magnitude = {:.3f}, dipole=({:.3f},{:.3f},{:.3f})'
              .format(charge_set,np.linalg.norm(tot_dipole), *tot_dipole))

if __name__ == '__main__':
    main()