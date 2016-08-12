from molecule import molecule
from itertools import cycle
from copy import deepcopy
import numpy as np

#class factory is used to make invocation a bit more intuitive
# compare creator.chain.create(*) vs creator.chain.chain(*)
def create(*args,**kwargs):
  return stickyBeadChain(*args,**kwargs)

class stickyBeadChain(molecule):
  def __init__(self, 
               length=20, 
               bigDiameter=1.0, 
               bigSequence=['A'], 
               bigBondType='bondA', 
               bigAngleType=None,
               bigBondLength=1.0,
               stickyDiameter=1.0,
               stickySequence=['B'], 
               stickyBondType='bondB', 
               stickyAngleType=None,
               stickyDihedralType=None,
               stickyBondLength=1.0/3.0,
               ):
    super(stickyBeadChain,self).__init__()
    self.name='stickyBeadChain'
    self.placed=False
    self.natoms=int(length*2)
    self.diameters=[bigDiameter]*length + [stickyDiameter]*length
    self.bodies=[-1 for i in range(length*2)]

    self.types  = [x for i,x in zip(range(length),cycle(bigSequence))]
    self.types += [x for i,x in zip(range(length),cycle(stickySequence))]

    self.bonds  = [[bigBondType,i,j] for i,j in zip(range(length-1),range(1,length))]
    self.bonds += [[stickyBondType,i,i+length] for i in range(length)]

    self.angles = []
    if bigAngleType is not None:
      self.angles = [[bigAngleType,i,j,k] for i,j,k in zip(range(length-2),range(1,length-1),range(2,length))]

    for i,j,k in zip(range(length,length*2),range(length),range(1,length+1)):
      if k==length:
        k = j-1
      self.angles.append([stickyAngleType,i,j,k])

    self.dihedrals=[]
    if stickyDihedralType is not None:
      r1 = range(length,2*length-1)
      r2 = range(0,length-1)
      r3 = range(1,length)
      r4 = range(length+1,2*length)
      indexZip = zip(r1,r2,r3,r4)
      self.dihedrals = [[stickyDihedralType,i,j,k,m] for i,j,k,m in indexZip]

    self.beadVol = length * (4.0/3.0)*np.pi*(bigDiameter/2.0)**(3.0)

    #big bead positions
    self.positions=list(self.getSpiralPos(natoms=int(length),bondLength=bigBondLength))

    #sticky bead positions
    stickyPos = deepcopy(self.positions)
    for i in range(length):
      stickyPos[i][2] += stickyBondLength

    self.positions.extend(stickyPos)
    self.positions = np.array(self.positions)

