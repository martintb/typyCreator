from molecule import molecule
import numpy as np

def create(*args,**kwargs):
  return bead(*args,**kwargs)

class bead(molecule):
  def __init__(self, diameter=[5.], type=['A']):
    super(bead,self).__init__()
    self.name='bead'
    self.placed=False
    self.natoms=1
    self.positions=np.array([[0.,0.,0.]])
    self.types=type
    self.diameters=diameter
    self.bonds=[]
    self.angles=[]
    self.dihedrals=[]
    self.bodies=[-1]
    self.beadVol=(4.0/3.0) * np.pi *(diameter[0]/2.0)**(3.0)
