from system import system
from itertools import cycle


def create(*args,**kwargs):
  return gromacs_system(*args,**kwargs)

class gromacs_system(system):
  def create_grofile(self,fname,resID=None):
    self.gatherMolData()
    sysDict = self.__dict__
    if len(self.residue_types)!=len(self.types):
      raise ValueError('Gromacs Output requires that residue types be specified!')
    print '.:: Opening {} for writing (will overwrite if already exists)'.format(fname)
    with open(fname,'w') as f:
      f.write('TYPY created GROMACS GRO file\n')
      f.write('{natoms:5d}\n'.format(**sysDict))

      if resID is None:
        resID = sysDict['moleculeIDs']
      atomZip = zip(resID,sysDict['residue_types'],sysDict['types'],sysDict['positions'])
      for i,(m,r,t,p) in enumerate(atomZip,start=1):
        lineDict = dict(atomID=int(i),
                        resID=int(m),
                        resName=r,
                        atomType=t,
                        x=p[0],
                        y=p[1],
                        z=p[2])
        f.write('{resID:5d}{resName:5s}{atomType:5s}{atomID:5d}{x:8.3f}{y:8.3f}{z:8.3f}\n'.format(**lineDict))
      f.write('{box.lx:10.5f} {box.ly:10.5f} {box.lz:10.5f}'.format(**sysDict))
      f.write('\n') #necessary for vmd molfile plugin
  def create_topfile(self,
                     forcefield='tymos',
                     fname='topol.top',
                     solvent='spce',
                     mol_dicts = []
                    ):
    with open(fname,'w') as f:
      f.write('#include "{}.ff/forcefield.itp"\n\n'.format(forcefield))

    for md in mol_dicts:
      f.write('[ moleculetype ]\n')
      f.write('; Name    nrexcl\n')
      f.write('; {mol.name:} {nrexcl:}\n\n'.format(**md))

      f.write('[ atoms ]\n')
      f.write(';nr  type  resnr residue  atom   cgnr     charge       mass\n')
      mol = md['mol']
      ad = zip(mol.types,
               mol.residue_numbers,
               mol.residue_names,
               mol.names,
               mol.residue_numbers,
               mol.charges,
               mol.masses)
      for i,vals in enumerate(ad):
        f.write('{:d} {:s} {d} {:s} {:s} {:d} {:f} {:f}\n'.format(i,*vals)


      f.write('\n[ bonds ]\n')
      f.write(';ai aj funct    c0\n')
      for i,j,k in mol.bonds:
        f.write('{:d} {:d} {:d} {:s}'.format(i,j,1,k))

      f.write('\n[ angles ]\n')
      f.write(';ai aj ak funct    c0\n')
      for i,j,k,l in mol.angles:
      f.write('{:d} {:d} {:d} {:d} {:s}'.format(i,j,k,1,l))

    f.write('\n\n#include "{}.ff/{}.itp"\n\n'.format(forcefield,solvent))

    f.write('[ system ]\n')
    f.write('TYPY Create Topology\n\n')

    f.write('[ molecules ]\n')
    for md in mol_dicts:
      f.write('{:s} {:d}\n'.format(md['mol'].name,md['count'])
