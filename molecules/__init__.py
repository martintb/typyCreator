import glob
import os

#Unfortunately, python does everything relative to the users cwd
#this allows for the dynamic importing of all files in this directory
molecule_dir = os.path.dirname(__file__)
molecule_paths = glob.glob(os.path.join(molecule_dir,'[!_]*.py'))
molecule_files = map(os.path.basename,molecule_paths)
molecule_names = [f.replace('.py','') for f in molecule_files]
__all__ = molecule_names

print '.:: Imported molecules:'
for i,f in enumerate(molecule_names):
  print '{:d}) {:s}'.format(i,f)
print ''
