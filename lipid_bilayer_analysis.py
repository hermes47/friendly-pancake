'''
Created on 9/03/2016

@author: iwelsh
'''
import numpy as np
from plotting import plot_scatters
from scipy.optimize import curve_fit

def calc_apl(x, y, N=64):
    return (np.array(x)*np.array(y)) / N

def plot_apl(x, y, time, N=64):
    apl = calc_apl(x, y, N)
    plot_scatters([{'x':np.array(time)/1000, 'y':apl, 'marker':'-', 'label':'APL'}], '/Volumes/DataStore/DPPE/SR_RF/Q1/LJ1/APL.png',
                  show_legend=False, x_label="Time (ns)", y_label=r'Area per Lipid (nm$^2$)',
                  title='Area per lipid for SR_RF-Q1-LJ1-DPPE')

def calc_vpl(x, y, z, NL=128, NW=5760, VW=0.03177):
    return (np.array(x)*np.array(y)*np.array(z) - NW * VW) / NL

def plot_vpl(x, y, z, time, NL=128, NW=5760, VW=0.03177):
    vpl = calc_vpl(x, y, z, NL, NW, VW)
    plot_scatters([{'x':np.array(time)/1000, 'y':vpl, 'marker':'-', 'label':'VPL'}], 
                  '/Volumes/DataStore/DPPE/SR_RF/Q1/LJ1/VPL.png',
                  show_legend=False, x_label="Time (ns)", y_label=r'Volume per Lipid (nm$^3$)',
                  title='Volume per lipid for SR_RF-Q1-LJ1-DPPE')
    
def plot_thickness(x, y):
    def linear(x,m,c):
        return m*x + c
    def linear_point(target,m,c):
        return (target - c)/m
    b_count = 9
    bulk_count = np.mean((np.mean(y[0:2*b_count]),np.mean(y[-2*b_count:])))
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
    
    plot_scatters([{'x':x, 'y':y, 'marker':'-', 'label':'water density'}], 
                  '/Volumes/DataStore/DPPE/SR_RF/Q1/LJ1/WaterDensity.png',
                  show_legend=False, x_label="Position (nm)", y_label=r'Water mass density',
                  title='Water mass density for SR_RF-Q1-LJ1-DPPE')
    return thickness
    

def main():
    time, x, y, z = [], [], [], []
    analysis_length = 50000    # ps
    with open('/Volumes/DataStore/DPPE/SR_RF/Q1/LJ1/dimensions.xvg','r') as fh:
        for line in fh:
            dat = list(map(float,line.split()))
            time.append(dat[0])
            x.append(dat[1])
            y.append(dat[2])
            z.append(dat[3])
    coord, density = [], []
    with open('/Volumes/DataStore/DPPE/SR_RF/Q1/LJ1/watermassdensity.xvg','r') as fh:
        for line in fh:
            dat = list(map(float, line.split()))
            coord.append(dat[0])
            density.append(dat[1])
    coord -= np.mean(coord)   # centre on 0
    plot_apl(x, y, time)
    plot_vpl(x, y, z, time)
    analysis_start = time.index(time[-1] - analysis_length)
    print('APL: {:.3f}'.format(np.mean(calc_apl(x[analysis_start:], y[analysis_start:]))))
    print('VPL: {:.3f}'.format(np.mean(calc_vpl(x[analysis_start:], y[analysis_start:], z[analysis_start:]))))
    print("Thickness: {:.3f}".format(plot_thickness(coord, density)))

if __name__ == '__main__':
    main()