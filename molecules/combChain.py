from molecule import molecule
from itertools import cycle
import numpy as np

#class factory is used to make invocation a bit more intuitive
# compare creator.combChain.create(*) vs creator.combChain.combChain(*)
def create(*args,**kwargs):
  return combChain(*args,**kwargs)

class combChain(molecule):
  def __init__(self, 
               backbone_length, 
               backbone_sequence,
               sidechain_length,
               sidechain_sequence,
               sidechain_loc,
               beadDiameters=[1.0], 
               bondTypes=[1],
               backboneAngleType=None,
               sidechainAngleType=None,
               dihedralType=None):
    super(combChain,self).__init__()
    # self.__dict__.update(**kwargs)
    self.backbone_length = backbone_length
    self.sidechain_length = sidechain_length
    self.num_sidechains = len(sidechain_loc)
    self.sidechain_loc = sorted(sidechain_loc)
    self.natoms = self.num_sidechains*sidechain_length + backbone_length
    self.name='combChain'
    self.placed=False
    self.dihedrals=[]
    self.positions=[]
    self.diameters=[d for i,d in zip(range(self.natoms),cycle(beadDiameters))]
    self.bodies=[-1 for i in range(self.natoms)]

    self.types=[]
    self.types += [x for i,x in zip(range(backbone_length),cycle(backbone_sequence))]
    for scl in self.sidechain_loc:
      self.types += [x for i,x in zip(range(sidechain_length),cycle(sidechain_sequence))]

    self.bonds = [[k,i,j] for k,i,j in zip(cycle(bondTypes),range(backbone_length-1),range(1,backbone_length))]
    count = backbone_length
    for i,scl in enumerate(self.sidechain_loc):
      self.bonds += [[bondTypes[0],scl,count]]
      self.bonds += [[k,i+count,j+count] for k,i,j in zip(cycle(bondTypes),range(sidechain_length-1),range(1,sidechain_length))]
      count+=sidechain_length

    self.angles=[]
    if backboneAngleType is not None:
      self.angles = [[backboneAngleType,i,j,k] for i,j,k in zip(range(backbone_length-2),range(1,backbone_length-1),range(2,backbone_length))]

    if sidechainAngleType is not None:
      count = backbone_length
      for i,scl in enumerate(self.sidechain_loc):
        self.angles += [[sidechainAngleType,scl,count,count+1]]
        self.angles += [[sidechainAngleType,i+count,j+count,k+count] for i,j,k in zip(range(sidechain_length-2),range(1,sidechain_length-1),range(2,sidechain_length))]
        count+=sidechain_length

    if dihedralType is not None:
      count = backbone_length
      for i,j in enumerate(range(len(self.sidechain_loc))):
        sca1 = self.sidechain_loc[j]
        try:
          sca2 = self.sidechain_loc[j+1]
        except IndexError:
          continue
        scb1 = backbone_length+sidechain_length*(i+0)
        scb2 = backbone_length+sidechain_length*(i+1)
        self.dihedrals += [[dihedralType,scb1,sca1,sca2,scb2]]

    self.beadVol=0
    for d in self.diameters:
      self.beadVol+=(4.0/3.0)*np.pi*(d/2.0)**(3.0)
  def makeExpandedComb(self,bondLength=1.0):
    pos = [[0.,0.,float(i)*bondLength] for i in range(self.backbone_length)]
    for i,scl in enumerate(self.sidechain_loc,start=1):
      zpos = pos[scl][-1]
      sidechain_pos = [[float(i+1)*bondLength,0.,zpos] for i in range(self.sidechain_length)]
      pos.extend(sidechain_pos)
    self.positions = np.array(pos)
  def makeExpandedAlternating(self,bondLength=1.0):
    pos = [[0.,0.,float(i)*bondLength] for i in range(self.backbone_length)]
    for i,scl in enumerate(self.sidechain_loc,start=1):
      zpos = pos[scl][-1]
      if (i%2)==0:
        sidechain_pos = [[-float(i+1)*bondLength,0.,zpos] for i in range(self.sidechain_length)]
      else:
        sidechain_pos = [[float(i+1)*bondLength,0.,zpos] for i in range(self.sidechain_length)]
      pos.extend(sidechain_pos)
    self.positions = np.array(pos)
  def makeLinear(self,bondLength=1.0):
    pos = [[0.,0.,float(i)*bondLength] for i in range(self.backbone_length)]
    for i,scl in enumerate(self.sidechain_loc,start=1):
      zpos = pos[scl][-1]
      sidechain_pos = [[0.,0.,float(i)*bondLength+zpos] for i in range(self.sidechain_length)]
      pos.extend(sidechain_pos)
    self.positions = np.array(pos)
  def makeLinear2(self,bondLength=1.0):
    pos = [[0.,0.,float(i)*bondLength] for i in range(self.backbone_length)]
    for i,scl in enumerate(self.sidechain_loc,start=1):
      zpos = pos[scl][-1]
      if (i%2)==0:
        sidechain_pos = [[0.,bondLength,float(i)*bondLength+zpos] for i in range(self.sidechain_length)]
      else:
        sidechain_pos = [[0.,-bondLength,float(i)*bondLength+zpos] for i in range(self.sidechain_length)]
      pos.extend(sidechain_pos)
    self.positions = np.array(pos)



