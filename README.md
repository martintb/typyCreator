# typyCreator #

The goal of typyCreator is to be able to create a wide number of random and customized intial configurations for coarse-grained MD and MC simulations of various molecules and building blocks. As much as possible, smaller building blocks are used to "build-up" larger building blocks e.g. the chain and sphere classes are used to build the graftedSphere molecule. 

In lieu of true documentation (for the time being), see the script below as an example of the usage. 

## Partial Listing of Supported Building Blocks ##
* bead
* chain
* surface
* star
* sphere
* graftedSphere
* graftedSurface
* graftedSurfaceAlt
* janusGraftedSphere
* stickyBead
* stickyBeadChain
* rodCoil
* water

## Example Python Script ##
The following script creates a binary system:

* 1000 copolymers of length 20 with repeating blocky sequence [1-1-1-2-2]*5
* 500 homopolymers of length 10 with repeating blocky sequence [1]*10

The placeMolecules command will randomly place all of chains in the box and ensure no overlaps using a relatively fast cythonized overlap checker. If a large number of molecules are to be placed in the box, typyCreator is also integrated with PACKMOL and just requires that the PACKMOL executable be on the users PATH. See [here ](http://www.ime.unicamp.br/~martinez/packmol/home.shtml)for more information on packmol.
```
import typyCreator
from typyCreator.systems import lammps_system

lsys = lammps_system.create()

lsys.masses[1]=1.12
lsys.masses[2]=1.03

chain1 = creator.chain.create(length=20,sequence=[1,1,1,2,2])
chain1.makeSpiral(diameter=1.0) #give the chain a spiral initial configuration
lsys.addMolecule(chain1,num=1000)

chain2 = creator.chain.create(length=20,sequence=[1])
chain2.makeSpiral(diameter=1.0) #give the chain a spiral initial configuration
lsys.addMolecule(chain2,num=500)

lsys.box.L = 100

lsys.placeMolecules(postFixBox=True)
#  lsys.placePACKMOL(tol=2.0,seed=np.random.randint(999999),randomize=True,useDiam=True):

lsys.atomTypes = set([1,2])
lsys.bondTypes = set([1])
lsys.create_datafile('data.lmpmolecular',atom_style='molecular')

```