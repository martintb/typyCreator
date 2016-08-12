import numpy as np
from ..cy.lib import maxDist,minDist
from copy import deepcopy

class molecule(object):
  '''
  Base molecule class. All molecules should inherit from this class
  and call the stron
  '''
  def __init__(self):
    self.placed=False
    self.natoms=0
    self.positions=[]
    self.names=[]
    self.types=[]
    self.masses=[]
    self.residue_numbers=[]
    self.residue_types=[]
    self.diameters=[]
    self.bonds=[]
    self.angles=[]
    self.dihedrals=[]
    self.bodies=[]
    self.charges=[]
    self.moleculeIDs=None
    self.moleculeNames=[]
    self.name='BaseClass'
  def shiftBondIDs(self,shiftVal):
      newBondCount     = 0
      newAngleCount    = 0
      newDihedralCount = 0
      bondTypeSet     = set()
      angleTypeSet    = set()
      dihedralTypeSet = set()
      newBonds=[]
      for bond in self.bonds:
        newBonds.append([bond[0], bond[1]+shiftVal, bond[2]+shiftVal])
        bondTypeSet.add(bond[0])
        newBondCount+=1
      newAngles = []
      for angle in self.angles:
        newAngles.append([angle[0], angle[1]+shiftVal, angle[2]+shiftVal, angle[3]+shiftVal])
        angleTypeSet.add(angle[0])
        newAngleCount+=1
      newDihedrals = []
      for dihedral in self.dihedrals:
        newDihedrals.append([dihedral[0], dihedral[1]+shiftVal, dihedral[2]+shiftVal, dihedral[3]+shiftVal, dihedral[4]+shiftVal])
        dihedralTypeSet.add(dihedral[0])
        newDihedralCount+=1
      self.bonds=newBonds
      self.angles=newAngles
      self.dihedrals=newDihedrals
      counts = {'bond':newBondCount,'angle':newAngleCount,'dihedral':newDihedralCount}
      types = {'bond':bondTypeSet,'angle':angleTypeSet,'dihedral':dihedralTypeSet}
      return {'counts':counts,'types':types}
  def shiftBody(self,bodyCount):
    newBodies = np.max(self.bodies)+1
    for i,bval in enumerate(self.bodies):
      if bval!=-1:
        self.bodies[i]+=bodyCount
    return newBodies
  def translateCOM(self,newCOM, L=1e5):
    diff =  newCOM - self.COM[0]
    self.positions = self.positions + diff
    self.COM = np.array([newCOM])
  def calcCOM(self):
    self.COM = np.array([np.mean(self.positions,axis=0)])
  def calcMaxR(self):
    self.maxR = maxDist(self.COM,self.positions,0,False)[2]
  def addMolecule(self,newMol):
    newMol = deepcopy(newMol)  
    newMol.shiftBondIDs(self.natoms)
    self.positions.extend(newMol.positions)
    self.types.extend(deepcopy(newMol.types))
    self.residue_types.extend(deepcopy(newMol.residue_types))
    self.diameters.extend(deepcopy(newMol.diameters))
    self.bonds.extend(deepcopy(newMol.bonds))
    self.angles.extend(deepcopy(newMol.angles))
    self.dihedrals.extend(deepcopy(newMol.dihedrals))
    self.bodies.extend(deepcopy(newMol.bodies))
    self.charges.extend(deepcopy(newMol.charges))
    self.natoms+=newMol.natoms
    for d in newMol.diameters:
      self.beadVol += (4.0/3.0) * np.pi * (d/2.0)**(3.0)
    if newMol.moleculeIDs is not None:
      if self.moleculeIDs is None:
        self.moleculeIDs = []
      self.moleculeIDs.extend(deepcopy(newMol.moleculeIDs))
    #self.moleculeNames.extend([newMol.name]*newMol.natoms)
  def toXYZ(self,fname):
    with open(fname,'w') as f:
      f.write('{}\n'.format(self.natoms))
      f.write('TYPY generated XYZ file for molecule: {}\n'.format(self.name))
      for i in range(self.natoms):
        t = str(self.types[i])
        x = self.positions[i,0]
        y = self.positions[i,1]
        z = self.positions[i,2]
        f.write('{:>3s} {:10.5f} {:10.5f} {:10.5f}\n'.format(t,x,y,z))
