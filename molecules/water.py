from molecule import molecule
import numpy as np

def create(*args,**kwargs):
  return water(*args,**kwargs)

class water(molecule):
  def __init__(self,
               typeO=2,
               typeH=1,
               diamO=3.166,
               diamH=0,
               chargeO=-0.8476,
               chargeH=+0.4238,
               rOH = 1.0,
               bondType=1,
               angle=109.47,
               angleType=1,
              ):
    super(water,self).__init__()
    self.name='water'
    self.placed=False
    self.natoms=3
    hx1 = 1
    hy1 = 0
    hx2 = rOH * np.cos(np.pi*angle/180)
    hy2 = rOH * np.sin(np.pi*angle/180)
    self.positions=np.array([
                              [hx1,hy1,0.],
                              [0.,0.,0.],
                              [hx2,hy2,0.],
                            ])
    self.types=[typeH,typeO,typeH]
    self.diameters=[diamH,diamO,diamH]
    self.bodies=[-1,-1,-1]
    self.charges = [chargeH,chargeO,chargeH]
    self.beadVol = (4.0/3.0) * np.pi *(diamO/2.0)**(3.0)
    self.angles = [[angleType,0,1,2]]
    self.bonds = [[bondType,0,1],
                  [bondType,1,2]]
