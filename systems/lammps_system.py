from system import system
# import ipdb; ist = ipdb.set_trace


def create(*args,**kwargs):
  return lammps_system(*args,**kwargs)

class lammps_system(system):
  def create_datafile(self,fname,atom_style='bond'):
    self.gatherMolData()
    sysDict = self.__dict__
    try:
      numAtomTypes = max(sysDict['atomTypes'])
    except ValueError:
      numAtomTypes = 0
    try:
      numBondTypes = max(sysDict['bondTypes'])
    except ValueError:
      numBondTypes = 0
    try:
      numAngleTypes = max(sysDict['angleTypes'])
    except ValueError:
      numAngleTypes = 0
    try:
      numDihedralTypes = max(sysDict['dihedralTypes'])
    except ValueError:
      numDihedralTypes = 0
    print '.:: Opening {} for writing (will overwrite if already exists)'.format(fname)
    with open(fname,'w') as f:
      f.write('TYPY created LAMMPS data file\n')
      f.write('{natoms:d} atoms \n'.format(**sysDict))
      f.write('{nbonds:d} bonds \n'.format(**sysDict))
      f.write('{:d} atom types \n'.format(numAtomTypes))
      f.write('{:d} bond types \n'.format(numBondTypes))
      if atom_style=='angle':
        f.write('{nangles:d} angles \n'.format(**sysDict))
        f.write('{:d} angle types \n'.format(numAngleTypes))
      if atom_style=='molecular' or atom_style=='full':
        f.write('{nangles:d} angles \n'.format(**sysDict))
        f.write('{:d} angle types \n'.format(numAngleTypes))
        f.write('{ndihedrals} dihedrals \n'.format(**sysDict))
        f.write('{:d} dihedral types \n'.format(numDihedralTypes))
        f.write('{:d} impropers \n'.format(0))
        f.write('{:d} improper types \n'.format(0))


      f.write('{box.xlo:f} {box.xhi:f} xlo xhi \n'.format(**sysDict))
      f.write('{box.ylo:f} {box.yhi:f} ylo yhi \n'.format(**sysDict))
      f.write('{box.zlo:f} {box.zhi:f} zlo zhi \n'.format(**sysDict))
      #f.write('0.0 0.0 0.0 xy xz yz \n')

      f.write('\n# Type Table\n')
      for i,t in enumerate(sysDict['atomTypes'],start=1):
        f.write('# {:d} {}\n'.format(i,t))

      f.write('\n# Bond Table\n')
      for i,t in enumerate(sysDict['bondTypes'],start=1):
        f.write('# {:d} {}\n'.format(i,t))

      f.write('\nMasses \n\n')
      for i,(p,m) in enumerate(sysDict['masses'].items(),start=1):
        f.write('{:d} {:f}\n'.format(p,m))

      f.write('\nAtoms # {}\n\n'.format(atom_style))
      if atom_style == 'full':
        atomZip = zip(sysDict['moleculeIDs'],sysDict['types'],sysDict['charges'],sysDict['positions'],sysDict['images'])
        for i,(m,t,q,p,im) in enumerate(atomZip,start=1):
          lineDict = dict(atomID=int(i),
                          molID=int(m)+1,
                          atomType=t,
                          charge=q,
                          x=p[0],
                          y=p[1],
                          z=p[2],
                          ix=im[0],
                          iy=im[1],
                          iz=im[2])
          f.write('{atomID:d} {molID:d} {atomType:d} {charge:f} {x:f} {y:f} {z:f} {ix:d} {iy:d} {iz:d}\n'.format(**lineDict))
      else:
        atomZip = zip(sysDict['moleculeIDs'],sysDict['types'],sysDict['positions'],sysDict['images'])
        for i,(m,t,p,im) in enumerate(atomZip,start=1):
          lineDict = dict(atomID=int(i),
                          molID=int(m)+1,
                          atomType=t,
                          x=p[0],
                          y=p[1],
                          z=p[2],
                          ix=im[0],
                          iy=im[1],
                          iz=im[2])
          f.write('{atomID:d} {molID:d} {atomType:d} {x:f} {y:f} {z:f} {ix:d} {iy:d} {iz:d}\n'.format(**lineDict))

      if sysDict['nbonds']>0:
        f.write('\nBonds \n\n')
        for i,b in enumerate(sysDict['bonds'],start=1):
          lineDict = dict(bondID=int(i),bondType=int(b[0]),x=b[1]+1,y=b[2]+1)
          f.write('{bondID:d} {bondType:d} {x:d} {y:d} \n'.format(**lineDict))

      if (sysDict['nangles']>0) and (atom_style=='angle' or atom_style=='molecular' or atom_style=='full'):
        f.write('\nAngles \n\n')
        for i,a in enumerate(sysDict['angles'],start=1):
          lineDict = dict(angleID=int(i),angleType=int(a[0]),x=a[1]+1,y=a[2]+1, z=a[3]+1)
          f.write('{angleID:d} {angleType:d} {x:d} {y:d} {z:d} \n'.format(**lineDict))

      if (sysDict['ndihedrals']>0) and (atom_style=='molecular' or atom_style=='full'):
        f.write('\nDihedrals \n\n')
        for i,a in enumerate(sysDict['dihedrals'],start=1):
          lineDict = dict(angleID=int(i),angleType=int(a[0]),w=a[1]+1,x=a[2]+1, y=a[3]+1,z=a[4]+1)
          f.write('{angleID:d} {angleType:d} {w:d} {x:d} {y:d} {z:d} \n'.format(**lineDict))

      f.write('\n\n\n\n\n\n')
    print '.:: Done writing {}!'.format(fname)
