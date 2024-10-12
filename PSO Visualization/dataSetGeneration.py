#created by harrison chung 
import math 
import random 
import numpy as np
from scipy.ndimage import gaussian_filter

class Dataset(object):
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z

    def findGlobalMin(self):
        x_min = self.x.ravel()[self.z.argmin()]
        y_min = self.y.ravel()[self.z.argmin()]
        z_min = self.z.ravel()[self.z.argmin()]
        return x_min, y_min, z_min

    def findGlobalMax(self):
        x_max = self.x.ravel()[self.z.argmax()]
        y_max = self.y.ravel()[self.z.argmax()]
        z_max = self.z.ravel()[self.z.argmax()]
        return x_max, y_max, z_max
    
def generateDataset(side,da,db,step):
    x = np.arange(0, side, step)
    y = np.arange(0, side, step)
    x,y = np.meshgrid(x,y)
    x = np.sin(x / 10) * 10
    y = np.sin(y / 10) * 10
    noise = np.random.normal(size = x.shape)
    z = gaussian_filter(noise, sigma=3.5)
    z = (z - np.min(z)) / (np.max(z) - np.min(z)) 
    z = z * (da - db) + db  
    s = Dataset(x, y, z)
    return s
