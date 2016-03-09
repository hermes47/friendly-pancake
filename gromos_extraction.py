'''
Created on 20/01/2016

@author: iwelsh
'''
from math import isnan

def extract_conf(loaded_file, key='int'):
    string = 'POSITION\n'
    conf = {}
    if string not in loaded_file:
        conf = None
    elif loaded_file[loaded_file.index(string) + 1].startswith('#'):
        start_line = loaded_file.index(string) + 2
    else:
        start_line = loaded_file.index(string) + 1
    count = 0
    while conf is not None and loaded_file[start_line+count] != 'END\n':
        line_data = loaded_file[start_line + count].split()
        if not line_data:
            break
        count += 1
        if key == 'int' or key not in ['int','name']:
            conf[count] = {'atom':  line_data[2],
                           'x':     float(line_data[4]),
                           'y':     float(line_data[5]),
                           'z':     float(line_data[6]),
                           'vec':   [float(line_data[4]),float(line_data[5]),float(line_data[6])]
                           }
            if isnan(conf[count]['x']) or isnan(conf[count]['y']) or isnan(conf[count]['z']):
                conf = None

        elif key == 'name':
            conf[line_data[2]] = {'atom':  count,
                                  'x':     float(line_data[4]),
                                  'y':     float(line_data[5]),
                                  'z':     float(line_data[6]),
                                  'vec':   [float(line_data[4]),float(line_data[5]),float(line_data[6])]
                                  }
    return conf

def extract_ene(loaded_file):
    string = 'MINIMIZED ENERGY\n'
    if string not in loaded_file:
        energy = 0.
    else:
        ene_line = loaded_file[loaded_file.index(string)+1]
        energy = float(ene_line.split()[2])
    return energy


def main():
    with open('/Users/iwelsh/GitHub/ExtractedData/Torsions/AdditionalFixedTorsions/Original/AMINO1/MD_data/AMINO1.240.aa.min.cnf','r') as fh:
        print(extract_conf(fh.readlines()))

if __name__ == '__main__':
    main()