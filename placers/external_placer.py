import numpy as np
import subprocess
from collections import defaultdict
from functools import partial
from ..util import overPrint

class external_placer(object):
  def placePACKMOL(self,tol=2.0,seed=np.random.randint(999999),randomize=True,useDiam=True):
    packStr = ''
    packStr += 'tolerance {}\n'.format(tol)
    packStr += 'filetype xyz\n'
    packStr += 'add_box_sides 1.0\n'
    packStr += 'output packed.xyz\n'
    packStr += 'seed {:d}\n'.format(seed)
    if randomize:
      packStr += 'randominitialpoint\n'
    packStr += '\n\n'
    boxStr = 'inside box {xlo} {ylo} {zlo} {xhi} {yhi} {zhi}\n'

    if useDiam:
      diamDict = defaultdict(partial(defaultdict,list))
      for i,(num,mol) in enumerate(self.uniqueMolecules,start=1):
        for j,diam in enumerate(mol.diameters,start=1):
          diamDict[i][diam].append(j)

    print '.:: Converting molecules to xyz files...'
    for i,(num,mol) in enumerate(self.uniqueMolecules,start=1):
      molName = mol.name + '.xyz'
      mol.toXYZ(molName)

      packStr += 'structure {}\n'.format(molName)
      packStr += '\tnumber {:d}\n'.format(int(num))
      packStr += '\t' + boxStr.format(**self.box.__dict__) 
      if useDiam:
        for diam,idex_list in diamDict[i].iteritems():
          idex_str = ' '.join(str(idex) for idex in idex_list)
          radius = diam/2.0
          packStr += '\tatoms {}\n'.format(idex_str)
          packStr += '\t\tradius {}\n'.format(radius)
          packStr += '\tend atoms\n'
      packStr += 'end structure\n'.format(molName)
      packStr += '\n\n'


    print '.:: Writing PACKMOL input...'
    with open('packmol.in','w') as f:
      f.write(packStr)

    print '.:: Calling PACKMOL...'
    with open('packmol.in','r') as f:
      subprocess.check_call('packmol',stdin=f)

    print '.:: Unpacking PACKMOL results...'
    packedPos = list(np.loadtxt('packed.xyz',skiprows=2))
    numMol = len(self.molecules)
    pstr = 'Unpacking PACKMOL placed molecule {{:d}} / {:d}'.format(numMol)
    for i,mol in enumerate(self.molecules,start=1):
      molPos = []
      molTypes = []
      overPrint(pstr.format(i))
      for n in range(mol.natoms):
        molData = packedPos.pop(0)
        molPos.append(molData[1:])
        molTypes.append(molData[0])
      if np.all(mol.types == molTypes):
        mol.positions = np.array(molPos)
      else:
        raise ValueError('Non-matching types in system and packed.xyz!')
    print ''



