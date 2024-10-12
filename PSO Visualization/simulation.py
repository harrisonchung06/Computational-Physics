#created by harrison chung 
import dataSetGeneration
from scipy.spatial import cKDTree
import numpy as np
import sys
from PyQt5.QtWidgets import QApplication
import pyqtgraph as pg
import pyqtgraph.opengl as gl
from pyqtgraph.Qt import QtCore, QtGui
from pyqtgraph.opengl import GLViewWidget, GLGridItem, GLLinePlotItem, GLMeshItem, GLScatterPlotItem

class Sim(object):
    def __init__(self, X,Y,Z, X_particles, side, db, da, step):
        #init the window 
        self.app = QApplication(sys.argv)
        self.view = GLViewWidget()
        self.view.setGeometry(0, 110, 1280, 720)
        self.view.show()
        self.view.setWindowTitle('Search Space')
        self.view.setCameraPosition(distance = 20)

        #create grid of size side with spacing of step 
        self.grid = GLGridItem()
        self.view.addItem(self.grid)
        self.grid.setSpacing(x=step, y=step)
        self.grid.setSize(x=side, y=side, z=da-db)  

        self.z=da+2
        self.grid2 = GLGridItem()
        self.view.addItem(self.grid2)
        self.grid2.setSpacing(x=step, y=step)
        self.grid2.setSize(x=side, y=side, z=da-db)  
        self.grid2.translate(0,0,self.z)
        

        #plotting X,Y,Z axes 

        self.x_axis = GLLinePlotItem(pos=np.array([(-side/2, -side/2, 0), (side, -side/2, 0)]), color=(1, 0, 0, 1), width=2, antialias=True)
        self.y_axis = GLLinePlotItem(pos=np.array([(-side/2, -side/2, 0), (-side/2, side, 0)]), color=(0, 1, 0, 1), width=2, antialias=True)
        self.z_axis = GLLinePlotItem(pos=np.array([(-side/2, -side/2, db), (-side/2, -side/2, da)]), color=(0, 0, 1, 1), width=2, antialias=True)

        self.view.addItem(self.x_axis)
        self.view.addItem(self.y_axis)
        self.view.addItem(self.z_axis)

        #generating faces 
        rows, cols = X.shape
        faces = []
        for r in range(rows - 1):
            for c in range(cols - 1):
                p1 = r * cols + c
                p2 = r * cols + (c + 1)
                p3 = (r + 1) * cols + c
                p4 = (r + 1) * cols + (c + 1)
                faces.append([p1, p2, p3])
                faces.append([p2, p4, p3])
        #generating surface 
        self.surface = GLMeshItem(
        vertexes=np.array([X.flatten(), Y.flatten(), Z.flatten()]).T,
        faces=np.array(faces), 
        color=QtGui.QColor(193,68,14), 
        shader = 'edgeHilight',
        smooth=True,
        drawEdges=False,
        drawFaces=True
        )
        self.surface.translate(-side/2, -side/2, 0)
        self.view.addItem(self.surface)

        
        self.X_particles = np.array([list([X_particles[0, i], X_particles[1, i], self.z]) for i in range(X_particles.shape[1])])
        self.particles = GLScatterPlotItem(
        pos=self.X_particles,
        color=QtGui.QColor(255, 0, 0), 
        size=0.1,  
        pxMode=False  
        )
        self.particles.translate(-side/2, -side/2, 0)
        self.view.addItem(self.particles)

    def start(self):
        #check before start 
        if (sys.flags.interactive != 1) or not hasattr(QtCore, 'PYQT_VERSION'):
            QApplication.instance().exec_()

    def update(self):
        update_particles()
        self.X_particles = np.array([list([X_particles[0, i], X_particles[1, i], self.z]) for i in range(X_particles.shape[1])])
        self.particles.setData(pos = self.X_particles)
        
    def animation(self):
        timer = QtCore.QTimer()
        timer.timeout.connect(self.update)
        timer.start(10)
        self.start()
        self.update()

def update_particles():
    global X_particles, V_particles, pbest, pbest_Z, gbest, gbest_Z
    r1 = np.random.rand(n_particles)
    r2 = np.random.rand(n_particles)
    V_particles = (w * V_particles +
            c1 * r1 * (pbest - X_particles) +
            c2 * r2 * (gbest.reshape(-1, 1) - X_particles))
    X_particles += V_particles
    X_particles = np.clip(X_particles, 0, side)
    Z_particles = get_z(X_particles[0], X_particles[1])
    better_mask = Z_particles < pbest_Z
    pbest[:, better_mask] = X_particles[:, better_mask]
    pbest_Z[better_mask] = Z_particles[better_mask]
    new_gbest_idx = Z_particles.argmin()
    if Z_particles[new_gbest_idx] < gbest_Z:
        gbest = X_particles[:, new_gbest_idx]
        gbest_Z = Z_particles[new_gbest_idx]

if __name__ == '__main__':
    side = 5
    db = -2
    da = 2
    step = 0.125 
    c1 = 0.1  
    c2 = 0.1  
    w = 0.8   

    dataSet = dataSetGeneration.generateDataset(side, da, db, step)
    X,Y,Z = dataSet.x, dataSet.y, dataSet.z

    x_flat = X.ravel()
    y_flat = Y.ravel()
    z_flat = Z.ravel()
    tree = cKDTree(np.vstack((x_flat, y_flat)).T)
    n_particles = 100
    #np.random.seed(100)
    X_particles = np.random.rand(2, n_particles) * side 
    V_particles = np.random.randn(2, n_particles) * 0.1 

    def get_z(x, y):
        coords = np.vstack((x, y)).T
        distances, indices = tree.query(coords, k=1)
        nearest_x, nearest_y = x_flat[indices.flatten()], y_flat[indices.flatten()]
        return z_flat[indices.flatten()]
   
    Z_particles = get_z(X_particles[0], X_particles[1])
    pbest = X_particles.copy()
    pbest_Z = Z_particles.copy()

    gbest_idx = Z_particles.argmin()
    gbest = X_particles[:, gbest_idx]
    gbest_Z = Z_particles[gbest_idx]
    
    s = Sim(X,Y,Z, X_particles, side, db, da, step)
    s.animation()

