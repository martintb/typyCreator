from box import box
import numpy as np
def sphereVol(d=1.0):
  return (4.0/3.0) * np.pi * (d/2.0)**(3.0)

import sys as pysys
pywrite = pysys.stdout.write
pyflush = pysys.stdout.flush
def overPrint(outStr):
  pywrite('\r' + outStr)
  pyflush()
