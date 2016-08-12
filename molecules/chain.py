from molecule import molecule
from itertools import cycle
import numpy as np

#class factory is used to make invocation a bit more intuitive
# compare creator.chain.create(*) vs creator.chain.chain(*)
def create(*args,**kwargs):
  return chain(*args,**kwargs)

class chain(molecule):
  def __init__(self, length=20, beadDiameters=[1.0], sequence=['A'], bondTypes=['bondA'],charges=[0], angleType=None,residue_types=None):
    super(chain,self).__init__()
    # self.__dict__.update(**kwargs)
    self.name='chain'
    self.placed=False
    self.natoms=length
    self.dihedrals=[]
    self.positions=[[0,0,0]]*length
    self.diameters=[d for i,d in zip(range(length),cycle(beadDiameters))]
    self.bodies=[-1 for i in range(length)]
    self.types=[x for i,x in zip(range(length),cycle(sequence))]

    if residue_types is not None:
      self.residue_types=[x for i,x in zip(range(length),cycle(residue_types))]
    self.charges=[x for i,x in zip(range(length),cycle(charges))]
    self.bonds = [[k,i,j] for k,i,j in zip(cycle(bondTypes),range(length-1),range(1,length))]
    self.angles=[]
    if angleType is not None:
      self.angles = [[angleType,i,j,k] for i,j,k in zip(range(length-2),range(1,length-1),range(2,length))]
    self.beadVol=0
    for d in self.diameters:
      self.beadVol+=(4.0/3.0)*np.pi*(d/2.0)**(3.0)
  def makeSpiral(self,bondLength=1.4):
    from ..util.geometry import spiralPos
    self.positions=spiralPos(natoms=self.natoms,bondLength=bondLength)
  def makeLinear(self,bondLength=1.4):
    from ..util.geometry import linearPos
    self.positions=linearPos(natoms=self.natoms,bondLength=bondLength)
  def makeZigZag(self,angle=110,bondLength=1.4,center=True):
    angle = np.pi*angle/180
    pos = [[0,0,0]]
    b = bondLength
    L = 2*b*np.sin(angle/2.0)
    h = b * np.cos(angle/2.0)
    for i in range(self.natoms-1):
      newX = pos[-1][0]+L/2.0
      if i%2==0:
        newY = h
      else:
        newY = 0
      pos.append([newX,newY,0])
    com = np.average(pos,axis=0)
    pos -= com
    self.positions = pos
  def makeCompact(self,lead=5):
    from ..util.geometry import compactPos
    self.positions=compactPos(natoms=self.natoms,lead=lead)
  def makeCompact2(self):
    from ..util.geometry import compactPos2
    self.positions=compactPos2(natoms=self.natoms)
  def makeCompact3(self):
    from ..util.geometry import compactPos3
    self.positions=compactPos3(natoms=self.natoms)


