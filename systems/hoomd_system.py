from hoomd_script import *
from system import system
import numpy as np
from itertools import cycle

def create(*args,**kwargs):
  return hoomd_system(*args,**kwargs)

class hoomd_system(system):
  def HOOMDload(self,beadTypes=[],bondTypes=[],angleTypes=[]):
    init.reset()
    if not beadTypes:
      beadTypes=np.unique([t for mol in self.molecules for t in mol.types])
      numTypes=int(np.max(map(ord,beadTypes))-64)
    if not bondTypes:
      bondTypes=['bondA']
    box = data.boxdim(Lx=self.box.lx,Ly=self.box.ly,Lz=self.box.lz)
    self.hoomd = init.create_empty(N=self.natoms, 
                                   box=box,
                                   particle_types = beadTypes, 
                                   bond_types=bondTypes,
                                   angle_types=angleTypes)
    i=0
    for mol in self.molecules:
      for pos, type, d, b in zip(mol.positions, mol.types, mol.diameters, mol.bodies):
        self.hoomd.particles[i].position=pos
        self.hoomd.particles[i].type=type
        self.hoomd.particles[i].diameter=int(d)
        self.hoomd.particles[i].body=int(b)
        i+=1
      for name,id1,id2 in mol.bonds:
        self.hoomd.bonds.add(name,int(id1),int(id2))
      if angleTypes:
        for name,id1,id2,id3 in mol.angles:
          self.hoomd.angles.add(name,int(id1),int(id2),int(id3))
    self.hoomd.sysdef.getRigidData().initializeData()
  def HOOMDcompress(self, BOXL=None, volumeFraction=0.1,
                    ljData={},
                    angleData=[{'k':5,'t0':np.pi}],
                    bondData =[{'k':50,'r0':1.4},{'k':50,'r0':5}],
                    preMix=True,
                    postMix=True,
                    scale_particles=False,
                    shrink_steps=0.5e6):
    print '.:: Starting BOX xyz',self.box.lx, self.box.ly,self.box.lz
    try:
      numBondTypes = self.hoomd.bonds.bdata.getNBondTypes()
    except AttributeError:
      numBondTypes = 1
    b = bond.harmonic()
    for i,bd in zip(range(numBondTypes),cycle(bondData)):
      b.bond_coeff.set('bond'+str(unichr(i+65)),k=bd['k'],r0=bd['r0'])

    numAngleTypes = self.hoomd.angles.adata.getNAngleTypes()
    if numAngleTypes>0:
      a = angle.harmonic()
      for i,ad in zip(range(numAngleTypes),cycle(angleData)):
        a.set_coeff('angle'+str(unichr(i+65)),k=ad['k'],t0=ad['t0'])
    else:
      a=[]

    allTypes= self.hoomd.particles.types
    lj = pair.lj(r_cut=np.power(2.0,1.0/6.0))
    lj.set_params(mode='xplor')
    lj.pair_coeff.set(allTypes,allTypes, epsilon=1.0, sigma=1.0, alpha=1.0)
    for key in ljData.iterkeys():
      pair0 = key[0]
      pair1 = key[1]
      eps = ljData[key]['eps']
      sigma = ljData[key]['sigma']
      rcut = ljData[key]['rcut']
      lj.pair_coeff.set(pair0,pair0, epsilon=eps, sigma=sigma, alpha=1.0, r_cut=rcut)

    xml1 = dump.xml(filename='pre-compress-sys.xml', all=True)
    mol1 = dump.mol2(filename='pre-compress-sys.mol2')
    dcd = dump.dcd(filename='compress-sys.dcd', period=0.5e5, overwrite=True)
    data=['potential_energy', 'kinetic_energy', 'pressure','temperature','pair_lj_energy','bond_harmonic_energy']
    log=analyze.log(filename='boxCompress.log',quantities=data, period=0.5e5, header_prefix='#')

    int1 = integrate.mode_standard(dt=0.0005)
    randomSeed1=np.random.randint(10,1e7)
    randomSeed2=np.random.randint(10,1e7)
    print '.:: Using the following random seeds for mixing & compression steps:',randomSeed1,randomSeed2
    int2 = integrate.bdnvt(group=group.nonrigid(), T=1, seed=randomSeed1, limit=0.01)
    int3 = integrate.bdnvt_rigid(group=group.rigid(), T=1, seed=randomSeed2)

    if preMix:
      print '.:: Running pre-compression mixing...'
      run1 = run(0.5e6)
      del run1

    initialL = self.hoomd.box[0]
    if not BOXL:
      vol=0
      for mol in self.molecules:
        vol+=mol.beadVol
      finalL=(vol/volumeFraction)**(1.0/3.0)
    else:
      finalL=BOXL
    print '.:: Shrinking box from',initialL,'to',finalL
    shrinkL = variant.linear_interp(points = [(0,initialL),(0.5e6, finalL) ])
    resizer = update.box_resize(Lx=shrinkL, period=500)
    resizer.set_params(scale_particles=scale_particles)

    run2 = run(shrink_steps)
    resizer.disable()

    if postMix:
      print '.:: Running post-compression mixing...'
      run3 = run(0.5e6)
      del run3
    xml2 = dump.xml(filename='post-compress-sys.xml', all=True)
    mol2 = dump.mol2(filename='post-compress-sys.mol2')

    del run2, resizer, shrinkL, initialL, finalL, mol, vol 
    del int1, allTypes, xml1, xml2, mol1,mol2
    if a:
      a.disable(); del a
    int2.disable();del int2
    int3.disable();del int3
    b.disable();   del b
    lj.disable();  del lj
    log.disable;   del log
    dcd.disable(); del dcd
    gc.collect();
    gc.collect();
