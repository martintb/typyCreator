import copy 
import numpy as np
from random import shuffle
from itertools import cycle
from molecule import molecule
from chain import chain
from sphere import sphere
from ..util import sphereVol
from ..util.geom import rotateTo


#class factory is used to make invocation a bit more intuitive
def create(*args,**kwargs):
  return graftedSphere(*args,**kwargs)

class graftedSphere(molecule):
  def __init__(self, 
               chainObj,
               chainKwargs,
               chainPosFunc,
               chainLengths=None,
               diameter=5,
               num_grafts=51, 
               angleType=None,
               anchorType='D',
               surfaceType='E',
               particleType='F',
               shuffleSequence=False,
               makeLinear=True,
               bondType=1,
               makeMolecules=False,
               ):
    super(graftedSphere2,self).__init__()
    self.name='graftedSphere'
    self.placed=False 
    self.natoms=0
    self.positions=[]
    self.types=[]
    self.dihedrals=[]
    self.diameters=[]
    self.bonds=[]
    self.bodies=[]
    self.angles=[]
    self.beadVol=0
    idCount = 0

    core = sphere(diameter=diameter, type=[particleType],bondType=bondType)
    core.addBeads(radius=diameter/2.0-0.5, type=[surfaceType])
    core.addBeads(radius=diameter/2.0-0.5, type=[anchorType], num_beads=num_grafts)
    self.addMolecule(core)
    self.beadVol+=sphereVol(diameter)

    if makeMolecules:
      self.moleculeIDs = [idCount]*self.natoms
      idCount += 1

    anchorIndex=copy.deepcopy(self.natoms)-1
    if chainLengths is not None:
      chainLengthIter=cycle(chainLengths)
    for i in range(num_grafts):
      if shuffleSequence:
        shuffle(chainKwargs['sequence'])
      if chainLengths is not None:
        chainKwargs['length'] = chainLengthIter.next()
      graft = chainObj(**chainKwargs)
      getattr(graft,chainPosFunc)()

      if makeMolecules:
        graft.moleculeIDs = [idCount]*graft.natoms
        idCount += 1

      index=int(-1*(i+1))
      anchorVec=core.positions[index]
      graftVec = graft.positions[1]-graft.positions[0]

      anchorNorm=(anchorVec[0]**2.0+anchorVec[1]**2.0+anchorVec[2]**2.0)**(0.5)
      monomer1Pos=[anchorVec[0]*(diameter/2.0+1.4-0.5)/anchorNorm,
                   anchorVec[1]*(diameter/2.0+1.4-0.5)/anchorNorm,
                   anchorVec[2]*(diameter/2.0+1.4-0.5)/anchorNorm]
      
      graft.positions = rotateTo(graft.positions,graftVec,anchorVec,transPos=monomer1Pos)

      self.bonds.extend([[bondType, anchorIndex,self.natoms]])
      if angleType:
        self.angles.extend([[angleType, anchorIndex,self.natoms,self.natoms+1]])
      anchorIndex-=1
      self.addMolecule(graft)
      self.beadVol+=graft.beadVol
