from chain import chain
from copy import deepcopy
import numpy as np

def create(*args,**kwargs):
  return rodCoil(*args,**kwargs)

class rodCoil(chain):
  def __init__(self,
               rodLength=7, 
               rodDiameters=[1.0],
               rodBondLength=1.0, 
               rodBondTypes=[1.0],
               rodSequence=[1], 
               rodAngleType=None,
               numCoils=3,
               coilLength=20, 
               coilDiameters=[2], 
               coilSequence=[2], 
               coilBondLength=1.0,
               coilBondTypes=[1.0],
               coilAngleType=None):
    super(rodCoil,self).__init__()
    self.name='rodCoil'
    self.placed=False
    self.natoms=0
    self.positions=[]
    self.diameters=[]
    self.bodies=[]
    self.types=[]
    self.dihedrals=[]
    self.bonds = []
    self.angles=[]
    self.beadVol=0

    #Make rod
    rod =  chain(length=rodLength,
                 beadDiameters=rodDiameters,
                 bondTypes=rodBondTypes,
                 angleType=rodAngleType,
                 sequence=rodSequence)
    rod.makeLinear(bondLength=rodBondLength)
    

    #make coils
    coil =  chain(length=coilLength,
                 beadDiameters=coilDiameters,
                 bondTypes=coilBondTypes,
                 angleType=coilAngleType,
                 sequence=coilSequence)
    coil.makeLinear(bondLength=coilBondLength)
    coilList = []
    newPos = np.array(rod.positions[0])
    newPos[0] += 1.13*(rod.diameters[0]+coil.diameters[0])/2.0
    th = 2*np.pi/numCoils
    for i in range(numCoils):
      cc = deepcopy(coil)
      diff = np.subtract(newPos,cc.positions[0])
      cc.positions = np.add(cc.positions,diff)
      coilList.append(cc)
      xy = newPos[:2] 
      xy = np.dot(xy,[[np.cos(th),-np.sin(th)],[np.sin(th),np.cos(th)]])
      newPos[:2]=xy
    


    #combine molecules
    self.addMolecule(rod)
    newBondList = []
    for c in coilList:
      newBondList.append([coilBondTypes[0],0,self.natoms])
      self.addMolecule(c)
    self.bonds.extend(newBondList)
    



