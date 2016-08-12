import numpy as np
from copy import deepcopy
from typy.misc import overPrint
from ..util.box import box
from ..placers.python_placer import python_placer
from ..placers.external_placer import external_placer

class system(python_placer,external_placer):
  def __init__(self): 
    self.molecules_placed=False
    self.molecules=[]
    self.uniqueMolecules = []
    self.natoms=0
    self.nbonds=0
    self.nangles=0
    self.ndihedrals=0
    self.nbodies=0
    self.atomTypes=set()
    self.bondTypes=set()
    self.angleTypes=set()
    self.dihedralTypes=set()
    self.box = box()
    self.masses={}
  def addMolecule(self, molecule,num=1):
    molecule.positions=np.array(molecule.positions)
    molecule.calcCOM()
    molecule.calcMaxR()
    molecule.positions = np.array(molecule.positions)
    molecule.diameters = np.array(molecule.diameters)
    self.uniqueMolecules.append([num,molecule])
    pstr = 'Adding molecule {{:d}} / {:d} to system...'.format(num)
    for i in range(num):
      overPrint(pstr.format(i+1))
      newMol=deepcopy(molecule)
      updates = newMol.shiftBondIDs(self.natoms)
      self.natoms  += newMol.natoms
      self.nbonds  += updates['counts']['bond']
      self.nangles += updates['counts']['angle']
      self.ndihedrals += updates['counts']['dihedral']
      if i==0:
        self.atomTypes.update(newMol.types)
        self.bondTypes.update(updates['types']['bond'])
        self.angleTypes.update(updates['types']['angle'])
        self.dihedralTypes.update(updates['types']['dihedral'])
      self.nbodies+=newMol.shiftBody(self.nbodies)
      self.molecules.append(newMol)
      self.box.beadVolume += newMol.beadVol
    print ''
  def gatherMolData(self):
    self.positions=[]
    self.diameters=[]
    self.bodies=[]
    self.types =[]
    self.residue_types =[]
    self.bonds = []
    self.angles=[]
    self.dihedrals=[]
    self.moleculeIDs=[]
    self.images=[]
    self.charges=[]
    self.moleculeNames=[]
    molCount = 0
    for mol in self.molecules:
      self.positions.extend(list(mol.positions))
      self.diameters.extend(list(mol.diameters))
      self.bodies.extend(list(mol.bodies))
      self.types.extend(list(mol.types))
      self.residue_types.extend(list(mol.residue_types))
      self.bonds.extend(list(mol.bonds))
      self.angles.extend(list(mol.angles))
      self.dihedrals.extend(list(mol.dihedrals))
      self.moleculeNames.extend([mol.name]*mol.natoms)
      if mol.charges is None:
        mol.charges = [0]*mol.natoms
      self.charges.extend(list(mol.charges))
      if mol.moleculeIDs is None:
        self.moleculeIDs.extend([molCount]*mol.natoms)
        molCount+=1
      else:
        newMolIDs = np.add(mol.moleculeIDs,molCount)
        self.moleculeIDs.extend(newMolIDs)
        molCount=np.max(newMolIDs)+1
      self.images.extend([[0,0,0]]*mol.natoms)
  def fitBox(self,iso=True):
    positions=[]
    for mol in self.molecules:
      positions.extend(list(mol.positions))
    print '.:: Fitting Box to Contained Atoms!'
    print 'Before fit:', self.box
    self.box.fit(positions,iso)
    print ' After fit:', self.box



