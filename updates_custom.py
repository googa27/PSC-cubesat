import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt
import time as tm
import calendar as cld
import quaternion
#from pyquaternion import Quaternion
import tlm_transcription2 as tr
from mpl_toolkits.mplot3d import Axes3D

#------------------------------------------------------------------------

earth_radius = 6378e3 #in meters
earth_magnetic_constant = 30.115e-9 #in tesla

#------------------------------------------------------------------------


#these functions transform vector coordinates

def gc_to_geo(vec, a, b):
    m1 = np.array([[np.cos(b),-np.sin(b)*np.sin(a),-np.sin(b)*np.cos(a)],
                   [0,np.cos(a),-np.sin(a)],
                   [np.sin(a),np.sin(a)*np.cos(b),np.cos(b)*np.cos(a)]])
    return np.dot(m1, vec)

def geo_to_gc(vec, a, b):
    m1 = np.array([[np.cos(b),-np.sin(b)*np.sin(a),-np.sin(b)*np.cos(a)],
                   [0,np.cos(a),-np.sin(a)],
                   [np.sin(a),np.sin(a)*np.cos(b),np.cos(b)*np.cos(a)]])
    return np.linalg.solve(m1, vec)

def geo_to_orb(vec, c):
    m2 = np.array([[np.cos(c),-np.sin(c),0],
                   [np.sin(c),np.cos(c),0],
                   [0,0,1]])
    return np.dot(m2, vec)

def orb_to_geo(vec, c):
    m2 = np.array([[np.cos(c),-np.sin(c),0],
                   [np.sin(c),np.cos(c),0],
                   [0,0,1]])
    return np.linalg.solve(m2, vec)

def gc_to_orb(vec, a, b, c):
    return geo_to_orb(gc_to_geo(vec, a, b), c)

def orb_to_gc(vec, a, b, c):
    return geo_to_gc(orb_to_geo(vec, c), a, b)

#------------------------------------------------------------------------



#------------------------------------------------------------------------
#here are the tests

if __name__ == '__main__':
    vgc = np.array([1, 0, 0])
    vgeo = gc_to_geo(vgc, np.pi, np.pi)
    gc = geo_to_gc(vgeo, np.pi, np.pi)

