from __future__ import print_function
'''
Created on 10/03/2016

@author: iwelsh
'''
import yaml

def main(cutoff, files):
    raw_data = []
    for file in files:
        with open(file,'r') as fh:
            raw_data += fh.readlines()
    cut_off = float(cutoff)
    uniqueish_chains = set()
    subs_list = {}
    previous = None
    for line in raw_data:
        comp_value = float(line.split()[1])
        a, b = line.split()[0].split('-')
        current = a.split('_')[0]
        if current != previous:
            uniqueish_chains.update(subs_list.values())
            subs_list = {}
        if a not in subs_list:
            subs_list[a] = a
        if b not in subs_list:
            subs_list[b] = b
        a, b = sorted([a,b], reverse=False)
        if comp_value > cut_off:  # halves are identical by comparison, so sublist update with the first one
            subs_list[b] = subs_list[a]
        previous = current
    uniqueish_chains.update(subs_list.values())
    print(len(uniqueish_chains),'unique-ish chains. List saved to uniqueish.txt')
    with open('uniqueish.txt','w') as fh:
        yaml.dump(sorted(uniqueish_chains),fh)
    print(subs_list)
                
def main2(cutoff, files):
    raw_data = []
    for file in files:
        with open(file,'r') as fh:
            raw_data += fh.readlines()
    cut_off = float(cutoff)
    uniqueish_chains = set()
    subs_list = {}
    previous = None
    for line in raw_data:
        comp_value = float(line.split()[1])
        a, b = line.split()[0].split('-')
        current = a.split('_')[0]
        if current != previous:
            uniqueish_chains.update(subs_list.values())
            subs_list = {}      
        if a in subs_list:
            continue
        if comp_value > cut_off:  # halves are identical by comparison, so sublist update with the first one
            subs_list[a] = b
        previous = current
    uniqueish_chains.update(subs_list.values())
    print(len(uniqueish_chains),'unique-ish chains. List saved to uniqueish.txt')
    with open('uniqueish.txt','w') as fh:
        for c in sorted(uniqueish_chains):
            fh.write(c + '\n')
    print(subs_list)           

if __name__ == '__main__':
    import sys
    main2(sys.argv[1], sys.argv[2:])
    main(sys.argv[1], sys.argv[2:])