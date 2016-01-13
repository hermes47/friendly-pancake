#!/usr/bin/env python3
'''
Created on 8/12/2015

@author: iwelsh
'''


def main():
    import os
    from time import sleep
    from random import shuffle
    source_file = '/Users/iwelsh/Downloads/LMSDSearchResultsDownload13H10M37S07Dec15.csv'
    dest_dir = '/Users/iwelsh/Documents/Lipid_MOL_files/'
    get_cmd = '/opt/local/bin/wget -O "{0}.mol" "http://www.lipidmaps.org/data/LMSDRecord.php?Mode=File&LMID={0}"'
    with open(source_file,'r') as fh:
        file_data = fh.readlines()[1:]
        shuffle(file_data)
    
    for line in file_data:
        s = line.split('","')
        dirt = s[5][:-5]
        dirt = dirt.replace(' ','_')
        sub_dirt = s[6][:-7]
        sub_dirt = sub_dirt.replace(' ','_')
        final_dir = dest_dir+dirt+'/'+sub_dirt
        if not os.path.isdir(final_dir):
            os.makedirs(final_dir)
        os.chdir(final_dir)
        os.system(get_cmd.format(s[0][1:]))
        sleep(15)
        
    

if __name__ == '__main__':
    main()