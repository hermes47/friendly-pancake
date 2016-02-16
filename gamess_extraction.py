'''
Functions for extracting data from GAMESS log files

Created on 14/01/2016

@author: iwelsh
'''
def extract_energy(loaded_file, scale=2625.5):
    """"
    Extract the minimised energy from a GAMESS optimisation log file
    """
    string = '      ***** EQUILIBRIUM GEOMETRY LOCATED *****\n'
    if string not in loaded_file:
        energy = 0.
    else:
        ene_line = loaded_file[loaded_file.index(string)-3]
        energy = float(ene_line.split()[3])*scale
    return energy

def extract_conformation(loaded_file, scale=1, key='int'): # or key='name'
    """
    Extracts the minimised conformation from a GAMESS optimisation log file
    """
    string = '      ***** EQUILIBRIUM GEOMETRY LOCATED *****\n'
    backup_string = 'BEGINNING GEOMETRY SEARCH POINT NSERCH=   0 ...\n'
    conf = {}
    if string not in loaded_file and backup_string not in loaded_file:
        conf = None
    elif string not in loaded_file:
        start_line = loaded_file.index(backup_string) + 5
    else:
        start_line = loaded_file.index(string) + 4
    count = 0
    while conf is not None:
        line_data = loaded_file[start_line + count].split()
        if not line_data:
            break
        count += 1
        if key == 'int' or key not in ['int','name']:
            conf[count] = {'atom':  line_data[0],
                           'charge':float(line_data[1]),
                           'x':     float(line_data[2])*scale,
                           'y':     float(line_data[3])*scale,
                           'z':     float(line_data[4])*scale,
                           'vec':   [float(line_data[2])*scale,float(line_data[3])*scale,float(line_data[4])*scale]
                           }
        elif key == 'name':
            conf[line_data[0]] = {'atom':  count,
                                  'charge':float(line_data[1]),
                                  'x':     float(line_data[2])*scale,
                                  'y':     float(line_data[3])*scale,
                                  'z':     float(line_data[4])*scale,
                                  'vec':   [float(line_data[2])*scale,float(line_data[3])*scale,float(line_data[4])*scale]
                                  }
    return conf
            
            
        
        

def main():
    with open('/Users/iwelsh/Documents/AdditionRigidTorsions/30PsiAllFixed/AMINO1/AMINO1.000.log','r') as fh:
        loaded_file = fh.readlines()
    print(extract_energy(loaded_file))
    print(extract_conformation(loaded_file))

if __name__ == '__main__':
    main()