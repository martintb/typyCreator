from chain import chain
from copy import deepcopy
import numpy as np

def create(*args,**kwargs):
  return CLP(*args,**kwargs)

class CLP(chain):
  def __init__(self,
               rodLength=15, 
               rodDiameters=[1.0],
               rodBondLength=1.0, 
               rodBondTypes=[1.0],
               rodSequence=[1], 
               rodAngleType=None,
               numRods=3,
               rodCentralDist=2.0,
               numCoils=3,
               coilLength=20, 
               coilDiameters=[2.0], 
               coilSequence=[2], 
               coilBondLength=1.0,
               coilBondTypes=[1],
               coilAngleType=None):
    super(CLP,self).__init__()
    self.name='CLP'
    self.placed=False
    self.natoms=0
    self.positions=[]
    self.dihedrals=[]
    self.diameters=[]
    self.bodies=[]
    self.types=[]
    self.bonds = []
    self.angles=[]
    self.beadVol=0
    self.moleculeIDs=[]

    #Make rod
    rod =  chain(length=rodLength,
                 beadDiameters=rodDiameters,
                 bondTypes=rodBondTypes,
                 angleType=rodAngleType,
                 sequence=rodSequence)
    rod.makeLinear(bondLength=rodBondLength)
    rodList = []
    newPos = np.array([0.,0.,0.]) 
    newPos[0] += rodCentralDist
    th = 2*np.pi/numRods
    for i in range(numRods):
      cc = deepcopy(rod)
      cc.moleculeIDs = [0]*rodLength
      diff = np.subtract(newPos,cc.positions[0])
      cc.positions = np.add(cc.positions,diff)
      rodList.append(cc)
      xy = newPos[:2]
      xy = np.dot(xy,[[np.cos(th),-np.sin(th)],[np.sin(th),np.cos(th)]])
      newPos[:2]=xy
    
    
    #make coils
    coil =  chain(length=coilLength,
                 beadDiameters=coilDiameters,
                 bondTypes=coilBondTypes,
                 angleType=coilAngleType,
                 sequence=coilSequence)
    coil.makeLinear(bondLength=coilBondLength)
    coilList = []
    newPos = np.array([0.,0.,0.])
    newPos[0] += 1.13*(rod.diameters[0]+coil.diameters[0])/2.0 + rodCentralDist
    th = 2*np.pi/numCoils
    for i in range(numCoils):
      cc = deepcopy(coil)
      cc.moleculeIDs = [i+1]*coilLength
      diff = np.subtract(newPos,cc.positions[0])
      cc.positions = np.add(cc.positions,diff)
      coilList.append(cc)
      xy = newPos[:2] 
      xy = np.dot(xy,[[np.cos(th),-np.sin(th)],[np.sin(th),np.cos(th)]])
      newPos[:2]=xy
    


    #combine molecules
    for k,r in enumerate(rodList):
      self.addMolecule(r)
    newBondList = []
    graftIndex = 0
    for c in coilList:
      newBondList.append([coilBondTypes[0],graftIndex,self.natoms])
      graftIndex += rodLength
      self.addMolecule(c)
    self.bonds.extend(newBondList)
    



