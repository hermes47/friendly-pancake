'''
Created on 14/01/2016

@author: iwelsh
'''
import os
from plotting import plot_scatters
from numpy import array, mean, std

def run_generate_traj(inroot, mol, charge, dispersion, electro, outroot):
    forms = dict(root=inroot, mol=mol, electro=electro, q=charge, lj=dispersion, n=0, out=outroot,
                 g_density='/opt/local/gromacs-4.6.7/bin/g_density_d')
    if not os.path.exists('{root}/{mol}/{electro}/{q}/{lj}/SplitTraj'.format(**forms)):
        os.makedirs('{root}/{mol}/{electro}/{q}/{lj}/SplitTraj'.format(**forms))
    if not os.path.exists('{out}/{mol}/{electro}/{q}/{lj}/Densities'.format(**forms)):
        os.makedirs('{out}/{mol}/{electro}/{q}/{lj}/Densities'.format(**forms))
        
    while os.path.isfile('{root}/{mol}/{electro}/{q}/{lj}/SplitTraj/npt_{n}.trr'.format(**forms)):
        forms['indir'] = '{root}/{mol}/{electro}/{q}/{lj}'.format(**forms)
        forms['outdir'] = '{out}/{mol}/{electro}/{q}/{lj}'.format(**forms)
        forms['lname'] = '{mol}_{electro}_{q}_{lj}_npt'.format(**forms)
        l_electron = ('{g_density} -f {indir}/SplitTraj/npt_{n}.trr -n {indir}/{mol}.ndx -s {indir}/{lname}.tpr '
                          '-ei {out}/{mol}/electrons.dat -o {outdir}/Densities/LipidElectron_{n:03}.xvg -xvg none '
                          '-sl 150 -dens electron'.format(**forms))
        w_electron = ('{g_density} -f {indir}/SplitTraj/npt_{n}.trr -n {indir}/{mol}.ndx -s {indir}/{lname}.tpr '
                          '-ei {out}/{mol}/electrons.dat -o {outdir}/Densities/WaterElectron_{n:03}.xvg -xvg none '
                          '-sl 150 -dens electron'.format(**forms))
        os.system('echo 3 | '+l_electron)
        os.system('echo 4 | '+w_electron)
        forms['n'] += 1
        
    

def load_data(electrostics, charges, dispersions, root, mol, loaded_data={}):
    for e in electrostics:
        if e not in loaded_data:
            loaded_data[e] = {}
        for q in charges:
            if q not in loaded_data[e]:
                loaded_data[e][q] = {}
            for lj in dispersions:
                forms = dict(root=root, elec=e, charge=q, lj=lj, mol=mol)
                if not os.path.isfile('{root}/{mol}/{elec}/{charge}/{lj}/energies_data.xvg'.format(**forms)):
                    continue
                time, x, y, z, l_temp, potential = [], [], [], [], [], []
                with open('{root}/{mol}/{elec}/{charge}/{lj}/energies_data.xvg'.format(**forms), 'r') as fh:
                    for line in fh:
                        split_data = list(map(float, line.split()))
                        time.append(split_data[0])
                        potential.append(split_data[1])
                        x.append(split_data[2])
                        y.append(split_data[3])
                        z.append(split_data[4])
                        l_temp.append(split_data[5])
                if lj not in loaded_data[e][q]:
                    loaded_data[e][q][lj] = {'e_time':array(time),
                                         'x':array(x),
                                         'y':array(y),
                                         'z':array(z),
                                         'l_temp':array(l_temp),
                                         'pot':array(potential),}
                else:
                    loaded_data[e][q][lj]['e_time'] = array(time)
                    loaded_data[e][q][lj]['x'] = array(x)
                    loaded_data[e][q][lj]['y'] = array(y)
                    loaded_data[e][q][lj]['z'] = array(z)
                    loaded_data[e][q][lj]['l_temp'] = array(l_temp)
                    loaded_data[e][q][lj]['pot'] = array(potential)
    return loaded_data
                
def plot_apl(electrostics, charges, dispersions, root, mol, loaded_data):
    for e in electrostics:
        if e not in loaded_data:
            continue
        for q in charges:
            if q not in loaded_data[e]:
                continue
            for lj in dispersions:
                if lj not in loaded_data[e][q]:
                    continue
                forms = dict(root=root, elec=e, charge=q, lj=lj, mol=mol)
                if not os.path.exists('{root}/Figures/AreaPerLipid'.format(**forms)):
                    os.makedirs('{root}/Figures/AreaPerLipid'.format(**forms))
                file_name = '{root}/Figures/AreaPerLipid/{mol}_{elec}_{charge}_{lj}.png'.format(**forms)
                apl = loaded_data[e][q][lj]['x'] * loaded_data[e][q][lj]['y'] / 64
                data = {'x':loaded_data[e][q][lj]['e_time']/1000,
                        'y':apl,
                        'marker':'-'}
                plot_scatters([data], file_name, show_legend=False, 
                  x_label='Time (ns)', y_label=r'Area per lipid nm$^2$', 
                  title='Area per lipid'.format(**forms))

def plot_vpl(electrostics, charges, dispersions, root, mol, loaded_data):
    pass

def plot_potential(electrostics, charges, dispersions, root, mol, loaded_data):
    for e in electrostics:
        if e not in loaded_data:
            continue
        for q in charges:
            if q not in loaded_data[e]:
                continue
            for lj in dispersions:
                if lj not in loaded_data[e][q]:
                    continue
                forms = dict(root=root, elec=e, charge=q, lj=lj, mol=mol)
                if not os.path.exists('{root}/Figures/Potential'.format(**forms)):
                    os.makedirs('{root}/Figures/Potential'.format(**forms))
                file_name = '{root}/Figures/Potential/{mol}_{elec}_{charge}_{lj}.png'.format(**forms)
                data = {'x':loaded_data[e][q][lj]['e_time']/1000,
                        'y':loaded_data[e][q][lj]['pot']/1000,
                        'marker':'-'}
                plot_scatters([data], file_name, show_legend=False, 
                  x_label='Time (ns)', y_label=r'Potential', 
                  title='Potential'.format(**forms))

def plot_isothermal_area_compressibility(electrostics, charges, dispersions, root, mol, loaded_data):
    pass

def calc_properties(electrostics, charges, dispersions, root, mol, loaded_data, time=75):
    for e in electrostics:
        if e not in loaded_data:
            continue
        for q in charges:
            if q not in loaded_data[e]:
                continue
            for lj in dispersions:
                if lj not in loaded_data[e][q]:
                    continue
                forms = dict(root=root, elec=e, charge=q, lj=lj, mol=mol)
                if not os.path.exists('{root}/Data'.format(**forms)):
                    os.makedirs('{root}/Data'.format(**forms))
                data = {}
                first_time = loaded_data[e][q][lj]['e_time'][-1] - time
                count = 0
                while loaded_data[e][q][lj]['e_time'][count] < first_time:
                    count += 1
                    
                apl = loaded_data[e][q][lj]['x'][count:] * loaded_data[e][q][lj]['y'][count:] / 64
                data['apl'] = {'name':'Area per Lipid:',
                               'units':'nm^2',
                               'mean':mean(apl),
                               'sd':std(apl)}
                with open('{root}/Data/{mol}_{elec}_{charge}_{lj}.txt'.format(**forms), 'w') as fh:
                    for x in data:
                        fh.write('{name:<25} {mean:>10.4f} +/- {sd:<10.4f} {units}\n'.format(**data[x]))

def main():
    electrostatics = ['PME','TR_RF','SR_RF']
    charges = ['Q1','Q2','Q3','Q4','Q5']
    dispersions = ['LJ1','LJ2','LJ3','LJ4','LJ5','LJ6']
    root = os.path.expanduser('~') + '/GitHub/ExtractedData/Bilayers'
    gendata_root = '/Volumes/DataStore'
    mol = 'DPPE'
    for e in electrostatics:
        for q in charges:
            for lj in dispersions:
                run_generate_traj(gendata_root, mol, q, lj, e, root)
    #loaded_data = load_data(electrostatics, charges, dispersions, root, mol)
    #plot_apl(electrostatics, charges, dispersions, root, mol, loaded_data)
    #plot_potential(electrostatics, charges, dispersions, root, mol, loaded_data)
    #calc_properties(electrostatics, charges, dispersions, root, mol, loaded_data)

if __name__ == '__main__':
    main()