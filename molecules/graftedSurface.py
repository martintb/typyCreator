from molecule import molecule
from itertools import cycle
from random import shuffle,choice,randint
import copy 
import numpy as np
from chain import chain
from ..util import sphereVol
from ..util.geometry import rotateTo


#class factory is used to make invocation a bit more intuitive
def create(*args,**kwargs):
  return graftedSurface(*args,**kwargs)

class graftedSurface(molecule):
  def __init__(self, 
               chainObj,
               chainKwargs,
               chainPosFunc,
               surfObj,
               num_grafts,
               chainLengths=None,
               bondTypes=[1],
               nx=1,
               dx=1,
               ny=1,
               dy=1,
               shuffleSequence=False):
    super(graftedSurface,self).__init__()
    self.name='graftedSurface'
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

    self.addMolecule(surfObj)
    self.beadVol += surfObj.beadVol

    if nx*ny<num_grafts:
      print '.:: Number graft anchors must be at least equal to number of grafts!'
      exit(1)

    topZ = np.unique(surfObj.positions[:,2]).max()
    botZ = np.unique(surfObj.positions[:,2]).min()

    topAnchorList = []
    botAnchorList = []
    for yi in range(ny):
      for xi in range(nx):
        x = xi*dx
        y = yi*dy
        topAnchorList.append([x,y,topZ])
        botAnchorList.append([x,y,botZ])
    shuffle(topAnchorList)
    shuffle(botAnchorList)

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
      # print '=============='
      # print topGraft.positions[np.array(topGraft.types)==2]

      topAnchor = topAnchorList.pop()
      botAnchor = botAnchorList.pop()

      R = np.eye(3)
      R[2,2] = -1
      botGraft.positions = np.dot(botGraft.positions,R)

      # topGraftVec = topGraft.positions[1] - topGraft.positions[0]
      # botGraftVec = botGraft.positions[1] - botGraft.positions[0]
      # topAnchorVec = [0,0,1]
      # botAnchorVec = [0,0,-1]
      # import ipdb; ipdb.set_trace()
      # topGraft.positions = rotateTo(topGraft.positions,topGraftVec,topAnchorVec)
      # botGraft.positions = rotateTo(botGraft.positions,botGraftVec,botAnchorVec)

      topPos = np.array(topAnchor)
      botPos = np.array(botAnchor)
      
      topPos[2] += 1.13*(surfObj.diameter+1.0)/2.0
      botPos[2] -= 1.13*(surfObj.diameter+1.0)/2.0

      topDiff = np.subtract(topPos,topGraft.positions[0])
      botDiff = np.subtract(botPos,botGraft.positions[0])

      topGraft.positions = np.add(topGraft.positions,topDiff)
      botGraft.positions = np.add(botGraft.positions,botDiff)

      bT = bondTypeIter.next()
      self.positions.append(topAnchor)
      self.types.append(surfObj.botType) #botType is on purpose for visualization
      self.diameters.append(surfObj.diameter)
      self.bodies.append(0)
      self.natoms+=1
      self.beadVol+=sphereVol(surfObj.diameter)
      self.bonds.append([bT, self.natoms-1,self.natoms])
      self.addMolecule(topGraft)
      self.beadVol+=sphereVol()*topGraft.natoms

      self.positions.append(botAnchor)
      self.types.append(surfObj.topType) #topType is on purpose for visualization
      self.diameters.append(surfObj.diameter)
      self.bodies.append(0)
      self.natoms+=1
      self.beadVol+=sphereVol(surfObj.diameter)
      self.bonds.append([bT, self.natoms-1,self.natoms])
      self.addMolecule(botGraft)
      self.beadVol+=sphereVol()*botGraft.natoms

    self.positions = np.array(self.positions)
