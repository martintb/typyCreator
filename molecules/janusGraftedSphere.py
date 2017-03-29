from molecule import molecule
from typyCreator.util import sphereVol
import numpy as np
from chain import chain

#class factory is used to make invocation a bit more intuitive
def create(*args,**kwargs):
  return janusGaftedSphere(*args,**kwargs)

class janusGraftedSphere(molecule):
  def __init__(self, diameter=4, graft_length=25, chain_type_list=['A']*12+['B']*12):
    super(janusGraftedSphere,self).__init__()
    self.name='janusGraftedSphere'
    self.placed=False 
    self.natoms=0
    self.dihedrals=[]
    self.positions=[]
    self.types=[]
    self.diameters=[]
    self.bonds=[]
    self.bodies=[]
    self.angles=[]
    self.beadVol=0

    core = sphere(diameter=diameter, type=['F'])
    core.addBeads(radius=diameter/2.0-0.5, type=['E'])
    core.addBeads(radius=diameter/2.0-0.5, type=['E'], num_beads=int(len(chain_type_list)*1.25))
    self.addMolecule(core)
    # self.bonds=[]
    self.beadVol+=sphereVol(diameter)

    anchorIndex=copy.deepcopy(self.natoms)-1
    corePos = list(core.positions)
    coreTypes = list(core.types)
    for chain_type in chain_type_list:
      graft = chain(length=graft_length,sequence=[chain_type],angleName=None)
      graft.makeLinear()

      foundGraftPoint=False
      while not foundGraftPoint:
        select = np.random.randint(1,len(corePos))
        if coreTypes[select] != 'E':
          corePos.pop(select)
          coreTypes.pop(select)
          foundGraftPoint=False
        elif corePos[select][2]>0 and chain_type=='A':
          anchorVec = corePos.pop(select)
          coreTypes.pop(select)
          foundGraftPoint=True
        elif corePos[select][2]<0 and chain_type=='B':
          anchorVec = corePos.pop(select)
          coreTypes.pop(select)
          foundGraftPoint=True
      
      graftVec = graft.positions[1]-graft.positions[0]

      anchorNorm=(anchorVec[0]**2.0+anchorVec[1]**2.0+anchorVec[2]**2.0)**(0.5)
      monomer1Pos=[anchorVec[0]*(diameter/2.0+1.4-0.5)/anchorNorm,
                   anchorVec[1]*(diameter/2.0+1.4-0.5)/anchorNorm,
                   anchorVec[2]*(diameter/2.0+1.4-0.5)/anchorNorm]
      
      graft.rotateTo(graftVec,anchorVec,transPos=monomer1Pos)

      self.bonds.extend([['bondA', anchorIndex,self.natoms]])
      anchorIndex-=1
      self.addMolecule(graft)
      self.beadVol+=sphereVol()*graft_length
  def compressChains(self):
    csys = system()
    csys.addMolecule(self)
    csys.setBOXL(fitToPositions=True)
    csys.HOOMDload()

    b = bond.harmonic()
    b.bond_coeff.set("bondA", k=50, r0=1.4)
    xml1 = dump.xml(filename='pre-compress-PGP.xml', all=True)
    mol1 = dump.mol2(filename='pre-compress-PGP.mol2')
    dcd = dump.dcd(filename='compress-PGP.dcd',period=100,overwrite=True)

    #pair interactions
    allTypes=csys.hoomd.particles.types
    lj = pair.lj(r_cut=2.5)
    lj.set_params(mode="xplor")
    lj.pair_coeff.set(allTypes, allTypes, epsilon = 5.00, sigma=1.0, alpha = 1.0)

    #define integration modes
    int1=integrate.mode_standard(dt=0.001)
    nonRigid=group.nonrigid()
    int2=integrate.nvt(group=nonRigid, tau=1, T=5)

    print '.:: Compressing chains of grafted sphere'
    runobj=run(1e4)

    self.positions=[copy.copy(p.position) for p in csys.hoomd.particles]
    xml2 = dump.xml(filename='post-compress-PGP.xml', all=True)
    mol2 = dump.mol2(filename='post-compress-PGP.mol2')

    del runobj, nonRigid, int1, allTypes,xml1, mol1, xml2, mol2
    int2.disable(); del int2
    dcd.disable(); del dcd
    b.disable(); del b
    lj.disable(); del lj
    del csys.hoomd
    del csys
    del p
    init.reset()
