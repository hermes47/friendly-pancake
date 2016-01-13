'''
Created on 14/12/2015

@author: iwelsh
'''

import matplotlib.pyplot as plt
from math import sqrt, pow
from copy import copy
from operator import add
import numpy as np

class non_bonded_plotter(object):
    
    # PLOT LIMITS
    V_XMIN = 0.15
    V_XMAX = 1.4
    V_YMIN = -225
    V_YMAX = 75
    
    F_XMIN = 0.15
    F_XMAX = 0.75
    F_YMIN = None
    F_YMAX = None
    
    PLT_START=75
    
    def __init__(self, h_parameters, other_parameters, f=138.935485):
        self.f = f
        self.h = h_parameters
        self.other = copy(other_parameters)
        if self.other['q'] is None:
            self.other['q'] = self.h['qnl']
        self.distances, self.potential, self.forces = [], [], []
        self.pot_coulomb, self.pot_c6s, self.pot_c12s = [], [], []
        self.force_coulomb, self.force_c6s, self.force_c12s = [], [], []
        
        
    def combine_parameters(self):
        self.combined_parameters = {'q1q2':self.f*self.h['qh']*self.other['q'],
                                    'c6':sqrt(self.h['c6']*self.other['c6']),
                                    'c12':sqrt(self.h['c12']*self.other['c12']),
                                    'name':'{}_{}'.format(self.other['name'],self.h['name'])}
        self.name = self.combined_parameters['name']
        
    def calc_interactions(self, table):
        def tri_add(a, b, c):
            return a + b + c
        
        for dist in table:
            x, V_charge, F_charge, V_c6, F_c6, V_c12, F_c12 = dist
            
            self.pot_coulomb.append(V_charge*self.combined_parameters['q1q2'])
            self.pot_c6s.append(V_c6*self.combined_parameters['c6'])
            self.pot_c12s.append(V_c12*self.combined_parameters['c12'])
            
            self.force_coulomb.append(F_charge*self.combined_parameters['q1q2'])
            self.force_c6s.append(F_c6*self.combined_parameters['c6'])
            self.force_c12s.append(F_c12*self.combined_parameters['c12'])
                        
            self.distances.append(x)
            
        self.potential = list(map(tri_add, self.pot_coulomb, self.pot_c6s, self.pot_c12s))
        self.forces = list(map(tri_add, self.force_coulomb, self.force_c6s, self.force_c12s))
            
    def plot_potential(self, outdir):
        plt.plot(self.distances[self.PLT_START:], self.potential[self.PLT_START:], 'b-', label='Total')
        plt.plot([0.8,0.8],[self.V_YMIN,self.V_YMAX],'k:')
        if not (-1e-16 < self.combined_parameters['q1q2'] < 1e-16):
            plt.plot(self.distances[self.PLT_START:], self.pot_coulomb[self.PLT_START:], 'c:', label='Coulomb')
        if not (-1e-16 < self.combined_parameters['c6'] < 1e-16):
            plt.plot(self.distances[self.PLT_START:], self.pot_c6s[self.PLT_START:], 'g:', label=r'C$_6$')
        if not (-1e-16 < self.combined_parameters['c12'] < 1e-16):
            plt.plot(self.distances[self.PLT_START:], self.pot_c12s[self.PLT_START:], 'r:', label=r'C$_{12}$')
        if (not (-1e-16 < self.combined_parameters['c6'] < 1e-16) 
            or not (-1e-16 < self.combined_parameters['c12'] < 1e-16)):
            plt.plot(self.distances[self.PLT_START:], list(map(add, self.pot_c6s, self.pot_c12s))[self.PLT_START:],
                     'm--',label=r'LJ')
        
        plt.title(r'''Potential of {} with {}
{}'''.format(self.other['name'], self.h['name'], self.determine_sigma_epsilon()))
        plt.xlim(xmin=self.V_XMIN, xmax=self.V_XMAX)
        plt.xlabel('Atom separation (nm)')
        plt.ylim(ymin=self.V_YMIN, ymax=self.V_YMAX)
        plt.ylabel(r'V$_{TOT}$ (kJ mol$^{-1}$)')
        plt.legend(loc=4, fontsize='small')
        plt.savefig('{}/POTENTIAL_{}.png'.format(outdir, self.combined_parameters['name']))
        plt.close()
        
    def plot_forces(self, outdir):
        plt.plot(self.distances[self.PLT_START:], self.forces[self.PLT_START:])
        
        plt.xlim(xmin=self.F_XMIN, xmax=self.F_XMAX)
        plt.ylim(ymin=self.F_YMIN, ymax=self.F_YMAX)
        plt.savefig('{}/FORCE_{}.png'.format(outdir, self.combined_parameters['name']))
        plt.close()
        
    def determine_sigma_epsilon(self):
        def find_nearest(l, target=0.):
            l2 = np.array(l)
            idx = (np.abs(l2 - target)).argmin()
            return idx, l2[idx]
        
        def calc_sigma():
            c12_over_q = -1*self.combined_parameters['c12']/self.combined_parameters['q1q2']
            sigma = pow(c12_over_q, 1/11)
            return sigma
        def calc_epsilon():
            c12_over_q = -12*self.combined_parameters['c12']/self.combined_parameters['q1q2']
            epsilon = pow(c12_over_q, 1/11)
            epsilon = self.combined_parameters['q1q2']/epsilon + self.combined_parameters['c12']/(epsilon**12)
            return epsilon
        
        cut_off = self.potential[400]    
        sorted_pot = sorted(self.potential)
        if int(np.sign(sorted_pot[0])) == 0 or int(np.sign(sorted_pot[-1])) == 0:
            return 'cut-off={:10.5f}'.format(cut_off)
        elif int(np.sign(sorted_pot[0])) == int(np.sign(sorted_pot[-1])):
            return 'cut-off={:10.5f}'.format(cut_off)
        else:
            idx, val = find_nearest(self.potential[50:])
            if self.combined_parameters['q1q2'] < 0 and -1e-16 < self.combined_parameters['c6'] < 1e-16:
                sigma = calc_sigma()
                epsilon = calc_epsilon()
            else:
                sigma = self.distances[idx+50]
                epsilon = sorted_pot[0]
            return '$\sigma$={:13.8f}, $\epsilon$={:10.5f}, cut-off={:10.5f}'.format(sigma,epsilon, cut_off)
    
    
    
    
def main():
    def parameter_iter(parameter_set, element='H'):
        q_count = 0
        for q in parameter_set['q']:
            q_count += 1
            c6_count = 0
            for c6 in parameter_set['c6']:
                c6_count += 1
                if c6_count >= 4:
                    continue
                c12_count = 0
                for c12 in parameter_set['c12']:
                    c12_count += 1
                    if c12_count >= 4:
                        continue
                    name = '{}_{:02}{:02}{:02}'.format(element,q_count, c6_count, c12_count)
                    yield dict(qh=q[0],qnl=q[1],c6=c6,c12=c12, name=name)
    
    h_parameters = {'q':[(0.400,-0.500), (0.350,-0.366), (0.248,0.129)],         # first in Hq, second is NLq
                    'c6':[0.000, 3.76996e-05, 8.464e-05],
                    #'c12':[0.000, 4.2999495e-09, 1e-08, 1.5129e-08, 2e-08, 2.5e-08, 5e-08, 7.5e-08, 
#                           1e-07, 2e-07, 2.5e-07, 5e-07, 7.5e-07, 1e-06, 2e-06, 2.5e-06, 5e-06, 7.5e-06, 1e-05,
#                           2e-05, 2.5e-05, 5e-05, 7.5e-05, 1e-04]
                    'c12':[1.5726e-07-.000001e-07*x for x in range(100)]
                    }
    
    non_iter_h_parameters = dict(#H1=dict(qh=0.248, qnl=0.129, c6=0.000, c12=0.000, name='H1l'),
                                 #H2=dict(qh=0.248, qnl=0.129, c6=3.76996e-05, c12=4.2999495e-09, name='H2l'),
                                 #H3=dict(qh=0.248, qnl=0.129, c6=8.464e-05, c12=1.5129e-08, name='H3l'),
                                 H7=dict(qh=0.248, qnl=0.129, c6=0.000, c12=1.0e-07, name='H4l'),
                                 #H9=dict(qh=0.248, qnl=0.129, c6=0.000, c12=1.0e-06, name='H5l'),
                                 #H4=dict(qh=0.400, qnl=-0.500, c6=0.000, c12=0.000, name='H1h'),
                                 #H5=dict(qh=0.400, qnl=-0.500, c6=3.76996e-05, c12=4.2999495e-09, name='H2h'),
                                 #H6=dict(qh=0.400, qnl=-0.500, c6=8.464e-05, c12=1.5129e-08, name='H3h'),
                                 #H8=dict(qh=0.400, qnl=-0.500, c6=0.000, c12=1.0e-07, name='H4h'),
                                 #H10=dict(qh=0.400, qnl=-0.500, c6=0.000, c12=1.0e-06, name='H5h'),
                                 #H16=dict(qh=0.400, qnl=-0.500, c6=0.000, c12=1.0e-03, name='H6h'),
                                 #H11=dict(qh=0.350, qnl=-0.366, c6=0.000, c12=0.000, name='H1m'),
                                 #H12=dict(qh=0.350, qnl=-0.366, c6=0.000, c12=4.2999495e-09, name='H2m'),
                                 #H13=dict(qh=0.350, qnl=-0.366, c6=0.000, c12=1.5129e-08, name='H3m'),
                                 #H14=dict(qh=0.350, qnl=-0.366, c6=0.000, c12=1.0e-07, name='H4m'),
                                 #H15=dict(qh=0.350, qnl=-0.366, c6=0.000, c12=1.0e-06, name='H5m'),
                                 )
    count = 0
    for c12 in h_parameters['c12']:
        non_iter_h_parameters[count] = dict(qh=0.311, qnl=-0.187, c6=0.000, c12=c12, name='H{}'.format(count))
        count += 1
    
    other_parameters = {#'O1':{'q':-0.700,
                        #     'c6':0.0022619536,
                        #     'c12':1e-06,
                        #     'name':'O1'},
                        #'O2':{'q':-0.600,
                        #     'c6':0.0022619536,
                        #     'c12':1e-06,
                        #     'name':'O2'},
                        'OM':{'q':-0.800,
                             'c6':0.0022619536,
                             'c12':7.4149321e-07,
                             'name':'OM'},
                        #'OA1':{'q':-0.800,
                        #     'c6':0.0022619536,
                        #     'c12':1.505529e-06,
                        #     'name':'OA1'},
                        #'OA2':{'q':-0.700,
                        #     'c6':0.0022619536,
                        #     'c12':1.505529e-06,
                        #     'name':'OA2'},
                        #'OE':{'q':-0.700,
                        #     'c6':0.0022619536,
                        #     'c12':1.21e-06,
                        #     'name':'OE'},
                        #'NL':{'q':None,
                        #     'c6':0.0024364096,
                        #     'c12':2.319529e-06,
                        #     'name':'NL'},
                        }
    
    table_path = 'GROMACS_Tables/table6-12.xvg'
    out_path = 'H_Interactions'
    with open(table_path, 'r') as fh:
        table = fh.readlines()
    for i in range(len(table)):
        table[i] = list(map(float,table[i].split()))
    plot_count = 0
    #for H in parameter_iter(h_parameters):
    for H in non_iter_h_parameters:
        H = non_iter_h_parameters[H]
        for other in other_parameters:
            #if plot_count > 10:
            #    break
            plotter = non_bonded_plotter(H, other_parameters[other])
            plotter.combine_parameters()
            plotter.calc_interactions(table)
            plotter.plot_potential(out_path)
            #plotter.plot_forces(out_path)
            plotter.determine_sigma_epsilon()
            plot_count += 1

def test_func():
    #from operator import add
    a = [1, 2, 3, 4]
    b = [1, 1, 1, 1]
    c = [2, 2, 2, 2]
    print(a)
    print(b)
    print(c)
    print(list(map(add, c, list(map(add, a, b)))))

if __name__ == '__main__':
    main()
    #test_func()