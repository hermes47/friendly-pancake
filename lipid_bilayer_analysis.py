'''
Created on 9/03/2016

@author: iwelsh
'''
import numpy as np
from plotting import plot_scatters
from scipy.optimize import curve_fit
import subprocess
import os
import traceback

class BilayerAnalysis(object):
    def __init__(self, mol, electrostaics, charge, lj,
                  root='/data/welsh/Projects/Headgroup_Parameterisation/PE_Headgroup_Parameterisation'):
        self.formatting = dict(mol=mol, electro=electrostaics, q=charge, lj=lj, root=root)
        os.chdir('{root}/{electro}/{q}/{lj}/{mol}'.format(**self.formatting))
        #os.chdir('{root}/{mol}/{electro}/{q}/{lj}'.format(**self.formatting))
        self.formatting['genName'] = '{mol}_{electro}_{q}_{lj}_npt'.format(**self.formatting)
        self.d_time, self.d_x, self.d_y, self.d_z = [], [], [], []
        self.clean_dir = True
        for file in os.listdir('.'):
            if '.pdb' in file:
                self.clean_dir = False
                break
        
        if electrostaics == "SR_RF":
            self.dimensions = [15,16,17,0]
            self.energies = [9,12,0]
        elif electrostaics == "TR_RF":
            self.dimensions = [17,18,19,0]
            self.energies = [11,14,0]
        elif electrostaics == "PME":
            self.dimensions = [15,16,17,0]
            self.energies = [9,12,0]
        else:
            self.dimensions = [15,16,17,0]
            self.energies = [9,12,0]
    
    def check_status(self):
        if not os.path.isfile('{mol}_{electro}_{q}_{lj}.out'.format(**self.formatting)):
            print('No Job Run')
            return False
        grep_command = ['grep','Too many LINCS warnings','{mol}_{electro}_{q}_{lj}.out'.format(**self.formatting)]
        grep = subprocess.Popen(grep_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, _ = grep.communicate()
        if out:
            print('Too many LINCS errors. Job failed')
            return False
        return True
            
    
    def run_gromacs_codes(self):
        g_energy_dimensions = ("g_energy_mpi_d -f {genName}.edr -s {genName}.tpr "
                              "-o dimensions.xvg -xvg none".format(**self.formatting))
        g_energy_energies = ("g_energy_mpi_d -f {genName}.edr -s {genName}.tpr "
                            "-o energies.xvg -xvg none".format(**self.formatting))
        g_density_watermass = ("g_density_mpi_d -f {genName}.gathered.trr -n {mol}.ndx -s {genName}.tpr "
                               "-o watermassdensity.xvg -xvg none -b 50000 -sl 250 -symm".format(**self.formatting))
        trjconv = ("trjconv_mpi_d -f {genName}.trr -n {mol}.ndx -s {genName}.tpr -o {genName}.gathered.trr "
                   " -xvg none -trans 0 0 4.5 -pbc mol -b 50000".format(**self.formatting))
        
        # run energy dimensions
        dimensions = subprocess.Popen(g_energy_dimensions.split(), stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                                      stderr=subprocess.PIPE)
        di_words = ' '.join(map(str,self.dimensions))
        _, _ = dimensions.communicate(di_words.encode())
        
        # run energy energies
        energies = subprocess.Popen(g_energy_energies.split(), stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE)
        en_words = ' '.join(map(str,self.energies))
        _, _ = energies.communicate(en_words.encode())
        
        # run trjconv
        gather = subprocess.Popen(trjconv.split(), stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                                  stderr=subprocess.PIPE)
        _, _ = gather.communicate('1'.encode())
        
        # run density
        density = subprocess.Popen(g_density_watermass.split(), stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)
        _, _ = density.communicate('0'.encode())
        
    def run_anaylsis(self):
        time, x, y, z = [], [], [], []
        analysis_length = 50000    # ps
        all_files = True
        if not os.path.isfile('dimensions.xvg'):
            print('No dimensions.xvg')
            all_files = False
        else:
            with open('dimensions.xvg','r') as fh:
                for line in fh:
                    dat = list(map(float,line.split()))
                    time.append(dat[0])
                    x.append(dat[1])
                    y.append(dat[2])
                    z.append(dat[3])
        
        potential, temp = [], []
        if not os.path.isfile('energies.xvg'):
            print('No energies.xvg')
            all_files = False
        else:
            with open('energies.xvg','r') as fh:
                for line in fh:
                    dat = list(map(float, line.split()))
                    potential.append(dat[1])
                    temp.append(dat[2])
        
        coord, density = [], []
        if not os.path.isfile('watermassdensity.xvg'):
            print('No watermassdensity.xvg')
            all_files = False
        else:
            with open('watermassdensity.xvg','r') as fh:
                for line in fh:
                    dat = list(map(float, line.split()))
                    coord.append(dat[0])
                    density.append(dat[1])
            if density[0] < 10:
                add = coord[-1]
                for i in range(len(coord)//2):
                    coord[i] += add
                density = density[len(density)//2:] + density[:len(density)//2]
                coord = coord[len(coord)//2:] + coord[:len(coord)//2] 
            coord -= np.mean(coord)   # centre on 0
        
        if all_files:
            self.plot_potential(time, potential)
            self.plot_apl(x, y, time)
            self.plot_vpl(x, y, z, time)
            analysis_start = time.index(time[-1] - analysis_length)
            print('APL: {:.3f}'.format(np.mean(self.calc_apl(x[analysis_start:], y[analysis_start:]))))
            print('VPL: {:.3f}'.format(np.mean(self.calc_vpl(x[analysis_start:], y[analysis_start:], z[analysis_start:]))))
            print("Thickness: {:.3f}".format(self.plot_thickness(coord, density)))
        
    def plot_potential(self, time, pot):
        plot_scatters([{'x':np.array(time)/1000, 'y':pot, 'marker':'-', 'label':'Pot'}], './Potential.png', 
                      show_legend=False, x_label="Time (ns)", y_label=r'Potential (kJ mol$^{-1}$)',
                      title='Potential for {genName}'.format(**self.formatting))
        
    def calc_apl(self, x, y, N=64):
        return (np.array(x)*np.array(y)) / N
    
    def plot_apl(self, x, y, time, N=64):
        apl = self.calc_apl(x, y, N)
        plot_scatters([{'x':np.array(time)/1000, 'y':apl, 'marker':'-', 'label':'APL'}], './APL.png',
                      show_legend=False, x_label="Time (ns)", y_label=r'Area per Lipid (nm$^2$)',
                      title='Area per lipid for {genName}'.format(**self.formatting))
    
    def calc_vpl(self, x, y, z, NL=128, NW=5760, VW=0.03177):
        return (np.array(x)*np.array(y)*np.array(z) - NW * VW) / NL
    
    def plot_vpl(self, x, y, z, time, NL=128, NW=5760, VW=0.03177):
        vpl = self.calc_vpl(x, y, z, NL, NW, VW)
        plot_scatters([{'x':np.array(time)/1000, 'y':vpl, 'marker':'-', 'label':'VPL'}], './VPL.png',
                      show_legend=False, x_label="Time (ns)", y_label=r'Volume per Lipid (nm$^3$)',
                      title='Volume per lipid for {genName}'.format(**self.formatting))
        
    def plot_thickness(self, x, y):
        def linear(x,m,c):
            return m*x + c
        def linear_point(target,m,c):
            return (target - c)/m
        b_count = 9
        bulk_count = np.mean((np.mean(y[0:2*b_count]),np.mean(y[-2*b_count:])))
        if not bulk_count > 10:
            bulk_count = np.mean(y[int(len(y)/2)-2*b_count:int(len(y)/2)+2*b_count])
        drop_below, raise_above = 0,0
        
        for i in range(len(y)):
            if y[i] < 0.5*bulk_count and not drop_below and x[i] < 0:
                drop_below = i
            elif y[i] > 0.5*bulk_count and not raise_above and x[i] > 0:
                raise_above = i
        x_low = np.asarray(x[drop_below-b_count:drop_below+b_count])
        x_high = np.asarray(x[raise_above-b_count:raise_above+b_count])
        y_low = np.asarray(y[drop_below-b_count:drop_below+b_count])
        y_high = np.asarray(y[raise_above-b_count:raise_above+b_count])
        popt_low, _ = curve_fit(linear,x_low,y_low,p0=[1,1])
        popt_high, _ = curve_fit(linear,x_high,y_high,p0=[1,1])
        thickness = linear_point(0.5*bulk_count,*popt_high) - linear_point(0.5*bulk_count,*popt_low)
        
        plot_scatters([{'x':x, 'y':y, 'marker':'-', 'label':'water density'}], './WaterDensity.png',
                      show_legend=False, x_label="Position (nm)", y_label=r'Water mass density',
                      title='Water mass density for {genName}'.format(**self.formatting))
        return thickness
    
    def cleanup(self):
        os.remove('{genName}.gathered.trr'.format(**self.formatting))

def testing():
    root = '/Volumes/DataStore/PE_Parameterisation'
    ana = BilayerAnalysis('DPPE','SR_RF','Q1','LJ1',root=root)
    ana.check_status()
    print(ana.clean_dir)
    #ana.run_gromacs_codes()
    ana.run_anaylsis()
    
def run_jobs():
    mols = ['DPPE','POPE']
    qs = ['Q{}'.format(x) for x in range(1,6)]
    ljs = ['LJ{}'.format(x) for x in range(1,7)]
    electros = ['SR_RF', 'TR_RF', 'PME',]
    for m in mols:
        for e in electros:
            for lj in ljs:
                for q in qs:
                    print('{}_{}_{}_{}'.format(m,e,q,lj))
                    ana = BilayerAnalysis(m, e, q, lj, root='/Volumes/DataStore/PE_Parameterisation')
                    if ana.check_status():
                        #try:
                        #    ana.run_gromacs_codes()
                        #except:
                        #    print('Run gromacs failed with:')
                        #    print(traceback.format_exc())
                        #    continue
                        try:
                            ana.run_anaylsis()
                        except:
                            print('Run analysis failed with:')
                            print(traceback.format_exc())
                            continue
                    print('\n')

if __name__ == '__main__':
    run_jobs()
    #testing()