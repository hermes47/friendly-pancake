#!/usr/bin/env python3
'''
Created on 3/02/2016

@author: iwelsh
'''
import sys

def load_file(path):
    # only load the position block. that's all we're worried about
    with open(path, 'r') as fh:
        data = {}
        box_size = []
        in_position = False
        in_genbox = False
        for line in fh:
            if line.startswith('#'):
                continue
            if line.startswith('POSITION'):
                in_position = True
                continue
            if line.startswith('GENBOX'):
                in_genbox = True
                continue
            if line.startswith('END'):
                in_position = False
                in_genbox = False
                continue
            if in_position:
                frag_num = int(line[0:6])
                frag_name = line[6:12].strip()
                atm_name = line[12:18].strip()
                atm_num = int(line[18:24])
                x_coord = float(line[24:39])
                y_coord = float(line[39:54])
                z_coord = float(line[54:69])
                data[atm_num] = dict(frag_num=frag_num,
                                     frag_name=frag_name,
                                     atm_name=atm_name,
                                     atm_num=atm_num,
                                     x=x_coord,
                                     y=y_coord,
                                     z=z_coord)
            if in_genbox:
                dat = line.split()
                if len(dat) < 3:
                    continue
                box_size = list(map(float,dat))
                in_genbox = False
                
    return data, box_size

def print_file(name, *source_files):
    xmin = source_files[0][1]['x']
    xmax = source_files[0][1]['x']
    ymin = source_files[0][1]['y']
    ymax = source_files[0][1]['y']
    zmin = source_files[0][1]['z']
    zmax = source_files[0][1]['z']
    coordinates = []
    coordinate_line = '{frag_num:>5} {frag_name:<5} {atm_name:<6}{atm_num:>6}{x:15.9f}{y:15.9f}{z:15.9f}'
    # add non_solvent atoms first, then add the solvent atoms
    frag_count = 0
    atom_count = 0
    for frag in ['OTHER','SOLV']:
        for file in source_files:
            previous_frag = None
            previous_frag_count = 0
            for atm in sorted(file):
                dat = file[atm]
                if frag == 'OTHER' and dat['frag_name'] == 'SOLV':
                    continue
                if frag == 'SOLV' and dat['frag_name'] != 'SOLV':
                    continue
                atom_count += 1
                if previous_frag != dat['frag_name'] or previous_frag_count != dat['frag_num']:
                    previous_frag = dat['frag_name']
                    previous_frag_count = dat['frag_num']
                    frag_count += 1
                if dat['x'] < xmin:
                    xmin = dat['x']
                if dat['x'] > xmax:
                    xmax = dat['x']
                if dat['y'] < ymin:
                    ymin = dat['y']
                if dat['y'] > ymax:
                    ymax = dat['y']
                if dat['z'] < zmin:
                    zmin = dat['z']
                if dat['z'] > zmax:
                    zmax = dat['z']
                dat['atm_num'] = atom_count
                dat['frag_num'] = frag_count
                coordinates.append(coordinate_line.format(**dat))
    out_format = dict(xbox = xmax - xmin,
                      ybox = ymax - ymin,
                      zbox = zmax - zmin,
                      coordinates = '\n'.join(coordinates))   
    with open(name, 'w') as fh:
        fh.write("""TITLE
Combined input files
END
POSITION
# first 24 chars ignored
{coordinates}
END
GENBOX
       1
    {xbox:.9f}    {ybox:.9f}    {zbox:.9f}
    90.00000000    90.00000000    90.00000000
    0.000000000    0.000000000    0.000000000
    0.000000000    0.000000000    0.000000000
END
""".format(**out_format))
   
def combine_files(file_1, file_2):
    dat_1, box_1 = load_file(file_1)
    dat_2, box_2 = load_file(file_2)
    center_coordinates(dat_1)
    center_coordinates(dat_2)
    z_shift = (box_2[2] + box_1[2])/2
    for atm in dat_2:
        dat_2[atm]['z'] += z_shift
    print_file('combined_coordinates.g96', dat_1, dat_2)
    
def center_coordinates(data):
    # find the centre shift needed
    x = y = z = 0
    for atm in data:
        x += data[atm]['x']
        y += data[atm]['y']
        z += data[atm]['z']
    x /= len(data)
    y /= len(data)
    z /= len(data)
    
    # shift to the centre
    for atm in data:
        data[atm]['x'] -= x
        data[atm]['y'] -= y
        data[atm]['z'] -= z
    #return data
            
def main():
    file_1 = sys.argv[1]
    file_2 = sys.argv[2]
    combine_files(file_1, file_2)

if __name__ == '__main__':
    main()