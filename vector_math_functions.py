'''
Created on 14/01/2016

@author: iwelsh
'''
import numpy as np

def dihedral_angle(a, b, c, d, period=[0,2*np.pi], scale='rad'):
    a = np.array(a)
    b = np.array(b)
    c = np.array(c)
    d = np.array(d)
    v1 = (a - b)/np.linalg.norm(a - b)
    v2 = (b - c)/np.linalg.norm(b - c)
    v3 = (c - d)/np.linalg.norm(c - d)
    n1 = np.cross(v1, v2)
    if np.linalg.norm(n1) > 0.0:
        n1 /= np.linalg.norm(n1)
    n2 = np.cross(v2, v3)
    if np.linalg.norm(n1) > 0.0:
        n2 /= np.linalg.norm(n2)
    m = np.cross(n1, v2)
    if np.linalg.norm(m) > 0.0:
        m /= np.linalg.norm(m)
    y, x = np.dot(m, n2), np.dot(n1, n2)
    phi = np.arctan2(y, x)
    if phi <= period[0]:
        phi += 2 * np.pi
    elif phi > period[1]:
        phi -= 2 * np.pi
    if scale == 'deg':
        phi *= 180/np.pi
        if int(round(phi,0)) == 360:
            phi -= 360
    return float(phi)

def main():
    pass

if __name__ == '__main__':
    main()