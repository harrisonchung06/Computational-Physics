import numpy as np
import math 
import sys
from PyQt5 import QtWidgets
from PyQt5.QtWidgets import QApplication
import pyqtgraph as pg
import pyqtgraph.opengl as gl
from pyqtgraph.Qt import QtCore, QtGui
from pyqtgraph.opengl import GLViewWidget, GLGridItem, GLLinePlotItem, GLMeshItem, GLScatterPlotItem
#import matplotlib.pyplot as plt

class SpringWorld:
    def __init__(self, side, height, pos):
        self.side = side
        self.height = height
        self.pos = pos

        self.app = QApplication(sys.argv)
        self.view = GLViewWidget()
        self.view.setGeometry(0, 110, 1280, 720)
        self.view.show()
        self.view.setWindowTitle('Spring World')
        self.view.setCameraPosition(distance = 25, elevation= 50)

        self.grid = GLGridItem()
        self.grid.setSpacing(x=1, y=1)
        self.grid.setSize(x=2*self.side, y=2*self.side)  
        self.grid.translate(self.side/2, self.side/2,0)
        self.view.addItem(self.grid)

        
        self.x_axis = GLLinePlotItem(pos=np.array([(-self.side/2, -self.side/2, 0), (self.side*(3/2), -self.side/2, 0)]), color=(1, 0, 0, 1), width=2, antialias=True) 
        self.y_axis = GLLinePlotItem(pos=np.array([(-self.side/2, -self.side/2, 0), (-self.side/2, self.side*(3/2), 0)]), color=(0, 1, 0, 1), width=2, antialias=True)

        self.view.addItem(self.x_axis)
        self.view.addItem(self.y_axis)

        self.particles = GLScatterPlotItem(
        pos=self.pos,
        color=QtGui.QColor(0, 0, 255), 
        size=0.3,  
        pxMode=False  
        )
        self.view.addItem(self.particles)
        self.energy_graph = EnergyGraph()
        self.energy_graph.show()

    def start(self):
        #check before start 
        if (sys.flags.interactive != 1) or not hasattr(QtCore, 'PYQT_VERSION'):
            QApplication.instance().exec_()

    def update(self):
        update_particles()
        self.pos = pos
        self.particles.setData(pos = pos)

    def vis(self):
        timer = QtCore.QTimer()
        timer.timeout.connect(self.update)
        timer.start(1)
        self.start()
    
class EnergyGraph:
    def __init__(self):
        
        self.window = QtWidgets.QMainWindow()
        self.window.setWindowTitle('Energy Graph')
        self.window.setGeometry(1300, 110, 640, 480)
        
        self.plot_widget = pg.PlotWidget()
        self.window.setCentralWidget(self.plot_widget)
        self.pe = []
        self.pe_curve = self.plot_widget.plot([], pen='b', name='Potential Energy')
        self.ke = []
        self.ke_curve = self.plot_widget.plot([], pen='r', name='Kinetic Energy')
        self.te = []
        self.te_curve = self.plot_widget.plot([], pen='g', name='Total Energy')

        self.timer = QtCore.QTimer()
        self.timer.timeout.connect(self.update_plot)
        self.timer.start(1)  

    def show(self):
        self.window.show()

    def update_plot(self):
        global te, pe, ke, particle_count
        self.te.append(te)
        self.pe.append(np.sum(np.linalg.norm(pe)))
        self.ke.append(np.sum(np.linalg.norm(ke)))
        self.te_curve.setData(self.te)
        self.pe_curve.setData(self.pe)
        self.ke_curve.setData(self.ke)
        

def update_particles():
    global pos,rij, mag_rij, f, fnet, v, a, te, pe,ke, m, r0, spring_const, dt, side, height

    #update loop for distance (rij), mag distance (rig/norm_rij), and forces 
    #forces f[] only contains the positive forces between atom i and atom j  
    #fnet[] is twice the size because it contains both positive and negative combinations of forces between each atom i and j, following newtons third law of motion
    #i is the current atom, j is all other atoms in the array, ahead of i and k is the total number of permutations between i and j   
    k=0
    for i in range(particle_count-1):
        for j in range(i+1, particle_count):
            rij[k] = pos[j]-pos[i]
            mag_rij[k] = math.sqrt((rij[k,0])**2+(rij[k,1])**2+(rij[k,2])**2)
            f[k] = -spring_const*(mag_rij[k]-r0)*(rij[k]/mag_rij[k])
            fnet[i] += -f[k]
            fnet[j] += f[k]
            pe[k] = 0.5*spring_const*(mag_rij[k]-r0)**2
            k+=1
    #kinematics update 
    a = fnet/m
    v = a*dt
    pos += v*dt+0.5*a*dt**2
    ke = 0.5*m*(v**2)
    te = np.sum(np.linalg.norm(ke))+np.sum(np.linalg.norm(pe))

    #bounds 
    for i in range(particle_count):
        if abs(pos[i,0]) >= side:
            v[i,0] = -v[i,0]
        if abs(pos[i,1]) >= side:
            v[i,1] = -v[i,1]
        if abs(pos[i,2]) >= height:
            v[i,2] = -v[i,2]

if __name__ == '__main__':
    #hyperparameters of the spring world 
    side = 10
    height = 10
    dt = 0.001

    np.random.seed(10) 

    #particle hyperparameters 
    particle_count = 8#number of particles 
    m = 1 #kg 
    spring_const = 50 #N/m
    v0 = [0,0,0] #initial force  
    r0 = 5

    #particle inits 
    #pos= np.random.rand(particle_count, 3)*side 
    pos = np.array(([0.0,0.0,0.0],[0.0,0.0,4.0], [0.0,4.0,0.0], [4.0,0.0,0.0], [4.0,4.0,4.0], [4.0,4.0,0.0], [4.0,0.0,4.0], [0.0,4.0,4.0]))  

    #global var inits 
    perm = int((particle_count*(particle_count+1))/2 - particle_count)

    rij = np.zeros((perm, 3))
    mag_rij =  np.zeros((perm, 1))
    f = np.zeros((perm, 3))
    fnet = np.zeros(pos.shape)
    v = np.zeros(fnet.shape)
    a = np.zeros(fnet.shape)
    pe = np.empty(rij.shape)
    ke = np.empty(v.shape)
    te = 0

    #init and start spring world 
    s = SpringWorld(side, height, pos)
    s.vis()
