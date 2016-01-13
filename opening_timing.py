'''
Created on 18/12/2015

@author: iwelsh
'''
import time
import matplotlib.pyplot as plt
def main():
    start_time = time.process_time()
    times = []
    number = 10000
    for i in range(number):
        with open('/Users/iwelsh/Documents/workspace/SmallWorkingScripts/GROMACS_Tables/table6-9.xvg','r') as fh:
            data = fh.readlines()
            for j in range(len(data)):
                data[j] = list(map(float,data[j].split()))
        times.append(time.process_time()-start_time)
        start_time = time.process_time()
        
    plt.plot([x for x in range(number)], times)
    plt.show()    
        

if __name__ == '__main__':
    main()