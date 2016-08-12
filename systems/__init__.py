import glob
import os

#Unfortunately, python does everything relative to the users cwd
#this allows for the dynamic importing of all files in this directory
system_dir = os.path.dirname(__file__)
system_paths = glob.glob(os.path.join(system_dir,'[!_]*.py'))
system_files = map(os.path.basename,system_paths)
system_names = [f.replace('.py','') for f in system_files ]
system_names = [f for f in system_names if (f!='system_base' and f!='box')]
__all__ = system_names

print '.:: Imported systems:'
for i,f in enumerate(system_names):
  print '{:d}) {:s}'.format(i,f)
print ''
