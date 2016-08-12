import copy 
import numpy as np
from random import shuffle
from itertools import cycle
from molecule import molecule
from chain import chain
from bead import bead
from ..util import sphereVol
from ..util.geometry import rotateTo


#class factory is used to make invocation a bit more intuitive
def create(*args,**kwargs):
  return star(*args,**kwargs)

class star(molecule):
  def __init__(self, 
               chainObj,
               chainKwargs,
               chainPosFunc,
               coreType,
               chainLengths=None,
               diameter=5,
               num_grafts=51, 
               angleType=None,
               bondType=1,
               shuffleSequence=False,
               makeMolecules=False,
               ):
    super(star,self).__init__()
    self.name='star'
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

    core = bead(diameter=[diameter], type=[coreType])
    self.addMolecule(core)
    self.beadVol+=sphereVol(diameter)

    if makeMolecules:
      self.moleculeIDs = [idCount]*self.natoms
      idCount += 1

    anchorIndex=copy.deepcopy(self.natoms)-1
    anchorVecs = np.random.random((num_grafts,3))-0.5
    anchorNorm  = (diameter+1.0)/2.0/np.linalg.norm(anchorVecs,axis=1)
    anchorVecs *= anchorNorm[:,np.newaxis]
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
      anchorVec= anchorVecs[index]
      graftVec = graft.positions[1]-graft.positions[0]
      graft.positions = rotateTo(graft.positions,graftVec,anchorVec,transPos=anchorVec)

      self.bonds.extend([[bondType, 0,self.natoms]])
      if angleType:
        self.angles.extend([[angleType, 0,self.natoms,self.natoms+1]])
      anchorIndex-=1
      self.addMolecule(graft)
      self.beadVol+=graft.beadVol
