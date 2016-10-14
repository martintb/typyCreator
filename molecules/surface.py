from molecule import molecule
from typy.creator.util import sphereVol
import numpy as np
from math import sqrt

def create(*args,**kwargs):
  return surface(*args,**kwargs)

class surface(molecule):
  def __init__(self, nx=5, ny=5, nz=3, diameter=1.0, topType='A',botType='B',midType='C',shift=True):
    super(surface,self).__init__()
    self.name='surface'
    self.placed=False
    self.natoms=0
    self.positions=[]
    self.types=[]
    self.diameters=[]
    self.dihedrals=[]
    self.bonds=[]
    self.angles=[]
    self.bodies=[]
    self.beadVol=sphereVol(diameter)*self.natoms
    self.nx = int(nx)
    self.ny = int(ny)
    self.nz = int(nz)
    self.diameter = diameter
    self.radius = diameter/2.0
    self.topType = topType
    self.botType = botType
    self.midType = midType
    for k in range(self.nz):
      for j in range(self.ny):
        for i in range(self.nx):
          x = (2*i + ((j+k)%2))*self.radius
          y = sqrt(3)*(j+(k%2)/3.0)*self.radius
          z = 2.0/3.0 * sqrt(6) * k * self.radius
          self.positions.append([x,y,z])
          self.natoms+=1
          self.diameters.append(diameter)
          self.bodies.append(0)
          if k==0:
            self.types.append(botType)
          elif k==(self.nz-1):
            self.types.append(topType)
          else:
            self.types.append(midType)
    self.positions = np.array(self.positions)
    if shift:
      self.positions[:,2] -= np.average(self.positions,axis=0)[2]


