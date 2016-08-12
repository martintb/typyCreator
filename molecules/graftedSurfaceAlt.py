from molecule import molecule
from typy.creator.util import sphereVol
from itertools import cycle
from random import shuffle,choice,randint
import copy 
import numpy as np
from chain import chain


#class factory is used to make invocation a bit more intuitive
def create(*args,**kwargs):
  return graftedSurfaceAlt(*args,**kwargs)

class graftedSurfaceAlt(molecule):
  def __init__(self, 
               chainObj,
               chainKwargs,
               chainPosFunc,
               surfObj,
               num_grafts,
               chainLengths=None,
               shuffleSequence=False,
               bondTypes = [1]):
    super(graftedSurfaceAlt,self).__init__()
    self.name='graftedSurfaceAlt'
    self.placed=False 
    self.natoms=0
    self.positions=[]
    self.dihedrals=[]
    self.types=[]
    self.diameters=[]
    self.bonds=[]
    self.bodies=[]
    self.angles=[]
    self.beadVol=0

    if surfObj is None:
      print '.:: Must supply a fully specified surface object!'
      exit(1)
    else:
      self.addMolecule(surfObj)
      self.beadVol += surfObj.beadVol

    anchorPos = copy.deepcopy(surfObj.positions)
    anchorTypes = np.array(surfObj.types)
    topType = surfObj.topType
    botType = surfObj.botType
    topAnchorIdex = list(np.where(anchorTypes==topType)[0])
    botAnchorIdex = list(np.where(anchorTypes==botType)[0])
    shuffle(topAnchorIdex)
    shuffle(botAnchorIdex)


    if chainLengths is not None:
      chainLengthIter=cycle(chainLengths)
    bondTypeIter=cycle(bondTypes)
    for i in range(num_grafts):
      if shuffleSequence:
        shuffle(chainKwargs['sequence'])
      if chainLengths is not None:
        chainKwargs['length'] = chainLengthIter.next()
      graft = chainObj(**chainKwargs)
      getattr(graft,chainPosFunc)()
      topGraft = copy.deepcopy(graft)
      botGraft = copy.deepcopy(graft)

      topIdex = topAnchorIdex.pop()
      botIdex = botAnchorIdex.pop()

      topPos = anchorPos[topIdex]
      botPos = anchorPos[botIdex]
      
      topPos[2] += 1.13*(surfObj.diameter+1.0)/2.0
      botPos[2] -= 1.13*(surfObj.diameter+1.0)/2.0

      topDiff = np.subtract(topPos,topGraft.positions[0])
      botDiff = np.subtract(botPos,botGraft.positions[-1])

      topGraft.positions = np.add(topGraft.positions,topDiff)
      botGraft.positions = np.add(botGraft.positions,botDiff)

      bT = bontTypeIter.next()
      self.bonds.append([bT, topIdex,self.natoms])
      self.addMolecule(topGraft)
      self.bonds.append([bT, botIdex,self.natoms])
      self.addMolecule(botGraft)
      self.beadVol+=2*sphereVol()*graft_length
    self.positions = np.array(self.positions)
