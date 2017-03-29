import numpy as np
from copy import deepcopy
from collections import defaultdict
import subprocess
from ..util import box
from ..util import overPrint
from ..cy.lib import minDist,minDistCheck

class python_placer(object):
  def placeMolecules(self, placement={},postFitBox=True,wrap=True,overlapCheck=True):
    pstr = 'Placing molecule {{:d}} / {:d}'.format(len(self.molecules))
    for i,trialMol in enumerate(self.molecules,start=1):
      if not trialMol.placed:
        overPrint(pstr.format(i))
        molType = ''.join([str(st) for st in sorted(np.unique(trialMol.types))])
        trialPos = self.getTrialPos(molType, placement)
        trialMol.translateCOM(trialPos, L=self.box.lx)
        if overlapCheck:
          while self.overlaps(trialMol):
            molType = ''.join([str(st) for st in sorted(np.unique(trialMol.types))])
            trialPos = self.getTrialPos(molType, placement)
            trialMol.translateCOM(trialPos, L=self.box.lx)
        trialMol.placed=True
    print ''
    if postFitBox:
      self.fitBox()
    self.molecules_placed=True
  def placeSolvent(self,solvent,num,placement={},postFitBox=True,solventOverlaps=True):
    self.gatherMolData()
    molPos = np.array(self.positions)
    molDiam = np.array(self.diameters)
    self.addMolecule(solvent,num)
    pstr = 'Placing solvent {{:d}} / {:d}'.format(num)
    m=1
    for i,trialMol in enumerate(self.molecules,start=1):
      if not trialMol.placed:
        overPrint(pstr.format(m))
        overLaps = True
        while overLaps:
          molType = ''.join([str(st) for st in sorted(np.unique(trialMol.types))])
          trialPos = self.getTrialPos(molType, placement)
          trialMol.translateCOM(trialPos, L=self.box.lx)
          if not minDistCheck(trialMol.positions,molPos,trialMol.diameters,molDiam,factor=1.13):
            overLaps=False
        trialMol.placed=True
        m+=1
        if not solventOverlaps:
          molPos = np.append(trialMol.positions,molPos,axis=0)
          molDiam = np.append(trialMol.diameters,molDiam)
    print ''
    if postFitBox:
      self.fitBox()
  def getTrialPos(self,molType,placement={}):
    try:
      minX = placement[molType]['minX']
      maxX = placement[molType]['maxX']
    except KeyError:
      minX =  self.box.xlo
      maxX =  self.box.xhi

    try:
      minY = placement[molType]['minY']
      maxY = placement[molType]['maxY']

    except KeyError:
      minY =  self.box.ylo
      maxY =  self.box.yhi

    try:
      minZ = placement[molType]['minZ']
      maxZ = placement[molType]['maxZ']
    except KeyError:
      minZ =  self.box.zlo
      maxZ =  self.box.zhi

    pos = []
    pos.append((maxX-minX)*np.random.rand(1)+minX)
    pos.append((maxY-minY)*np.random.rand(1)+minY)
    pos.append((maxZ-minZ)*np.random.rand(1)+minZ)
    return np.array(pos).flatten()
  def overlaps(self,trialMol=[]):
    if trialMol:
      mol2List = [trialMol]
    else:
      mol2List=self.molecules

    for mol2 in mol2List:
      for mol1 in self.molecules:
        if mol1.placed and mol1 is not mol2:
          if self.softCheck(mol1,mol2):
            if self.hardCheck(mol1,mol2):
              return True
    return False
  def softCheck(self,mol1,mol2):
    '''
    Check COM distance for overlaps
    '''
    dist = minDist(mol1.COM,mol2.COM,0,False)[2]
    if dist<=(mol1.maxR+mol2.maxR+1.13):
      return True
    return False
  def hardCheck(self,mol1,mol2):
    '''
    Check atom by atom for overlaps
    '''
    imin,jmin,dist=minDist(mol1.positions,mol2.positions,0,False)
    if dist<=(1.13)*((mol1.diameters[imin]+mol2.diameters[jmin])/2.0):
      return True
    return False

