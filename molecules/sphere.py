from molecule import molecule
from typyCreator.util import sphereVol
import numpy as np

def create(*args,**kwargs):
  return sphere(*args,**kwargs)

class sphere(molecule):
  def __init__(self, diameter=5., type=['A'],bondType='bondA'):
    super(sphere,self).__init__()
    self.name='sphere'
    self.placed=False
    self.natoms=1
    self.positions=np.array([[0.,0.,0.]])
    self.types=type*1
    self.diameters=[diameter]
    self.bonds=[[bondType,0,1]]
    self.angles=[]
    self.dihedrals=[]
    self.bodies=[0]
    self.beadVol=sphereVol(diameter)
  def addBeads(self, radius=2.5, bead_diameter=1, num_beads=None, type=['A']):
    if not num_beads:
      num_beads=int(4*(2*radius)**2.)
    r=2*np.random.rand(num_beads,3)-1 #initial random vectors
    rDist=np.sqrt(np.sum(r**2,axis=1))
    r = r*radius/(rDist[:,np.newaxis])
  
    damp=1; #damping coefficent for equilibration
    notEq=True
    while notEq:
      rPrev=r.copy()
      for i in range(num_beads):
        #difference vector btwn current point and all others
        rDiff = r[i,:]-r
  
        #distance between current point and all others
        rDist = np.sqrt(np.sum(rDiff**2,axis=1))
        rDist[i] = 1 #need to account for self vector
  
        #calculate force between particles
        force=1/(rDist**2)
        force[i] = 0 #need to account for self vector
  
        #normalize difference vectors
        rNormDiff = rDiff/rDist[:, np.newaxis]
  
        #total force on each particle
        rForce = np.sum(rNormDiff*force[:, np.newaxis],axis=0)
  
        #move current particle in direction of force
        rNormForce=rForce/np.linalg.norm(rForce)
        r[i,:]=r[i,:] + damp*rNormForce
  
        #put particle back onto sphere
        r[i,:] = r[i,:]*radius*(np.linalg.norm(r[i,:]))**(-1);
      
      dev=r-rPrev
      dev=sum(np.sqrt(np.sum(dev**2,axis=1)))/num_beads;
  
      if (dev<0.01):
        notEq=False

    self.positions=np.append(self.positions,r,axis=0)
    self.natoms+=num_beads
    self.types.extend(type*num_beads)
    self.diameters.extend([float(bead_diameter)]*num_beads)
    self.bodies.extend([0]*num_beads)
