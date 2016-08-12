from molecule import molecule
import numpy as np

def create(*args,**kwargs):
  return bead(*args,**kwargs)

class bead(molecule):
  def __init__(self, 
               bigDiameter=1.0, 
               bigType='A',
               stickyDiameter=1.0,
               stickyType='B',
               stickyBondType='bondA',
               stickyBondLength=1.0/3.0,
               ):
    super(bead,self).__init__()
    self.name='bead'
    self.placed=False
    self.natoms=2
    self.positions=np.array([[0.,0.,0.],[0.,0.,stickyBondLength]])
    self.types=[bigType,stickyType]
    self.diameters=[bigDiameter,stickyDiameter]
    self.bonds=[[stickyBondType, 0, 1]]
    self.dihedrals=[]
    self.angles=[]
    self.bodies=[-1,-1]
    self.beadVol=(4.0/3.0) * np.pi *(bigDiameter/2.0)**(3.0)
