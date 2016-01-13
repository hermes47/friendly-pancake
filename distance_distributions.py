#!/usr/bin/env python3
'''
Created on 14/12/2015

@author: iwelsh
'''
from numpy.linalg import norm
from numpy import array, histogram
def dist(a, b, box_size=[5.76157,5.87282,9.70528]):
    norms = []
    for min_x in [b[0],b[0]+box_size[0],b[0]-box_size[0]]:
        for min_y in [b[1],b[1]+box_size[1],b[1]-box_size[1]]:
            for min_z in [b[2],b[2]+box_size[2],b[2]-box_size[2]]:
                norms.append(norm(array(a)-array((min_x,min_y,min_z))))
    return min(norms)

def main():
    pdb_file = '/Users/iwelsh/Documents/Parametrisation/DPPE_lipids/QM_Calcs/lipids_only.gro'
    with open(pdb_file,'r') as fh:
        file = fh.readlines()
    for i in range(len(file)):
        file[i] = file[i].split()
    print('File loaded and split')
        
    all_coords = {}
    coords = {}
    count = 0
    for line in file:
        if len(line) in [1,7]: # skip headers
            continue
        if len(line) == 3: # create new set of coordinates and save this line as box dimensions
            all_coords[count] = dict(box=list(map(float,line)), coords=coords)
            coords = {}
            count += 1
            continue
        if line[1] in ['H1','H2','H3','O1','O2','O3','O4','O5','O6','O7','O8','N1']:
            if line[1] not in coords:
                coords[line[1]] = []
            coords[line[1]].append([float(line[3]),float(line[4]),float(line[5])])
    print('Coordinates determined')
            
    others = dict(N=['N1'],OA=['O1','O2'],OM=['O3','O4'],OE=['O5','O7'],O=['O6','O8'])
    results = {}
    for other in others:
        distances = []
        for c in range(0,len(all_coords),50):
            coords = all_coords[c]['coords']
            box = all_coords[c]['box']
            for H in ['H1','H2','H3']:
                for O in others[other]:
                    for start in coords[H]:
                        for end in coords[O]:
                            distances.append(dist(start,end, box))
        
        counts, bins = histogram(distances, range=(min(distances),1.4), bins=25)
        results[other] = [list(counts), list(bins)]
        print('done ',other)
        
    for x in results:
        print(x, results[x])
        
                
    
    
        
    

if __name__ == '__main__':
    main()