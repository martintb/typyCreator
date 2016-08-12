#cython: boundscheck=False
# #cython: wraparound=False
#cython: cdivision=True
#cython: nonecheck=False
import numpy as np
cimport numpy as np
cimport cython
from cython.parallel cimport prange
from libc.math cimport fabs as c_fabs
from libc.math cimport sqrt as c_sqrt

DTYPE = np.int
DTYPE2 = np.float #this is a double in cython
DTYPE3 = np.float32 #this is a double in cython
ctypedef np.int_t DTYPE_t
ctypedef np.float_t DTYPE2_t
ctypedef np.float32_t DTYPE3_t


@cython.boundscheck(False)
cdef inline double dist_wrap(double x, double BOXL):
  '''
  It's important to note that this wrapping method only
  works when you only care about the magnitude of the 
  resulting vector.
  '''
  cdef double Bo2=BOXL/2.0
  if x>Bo2:
    return x-BOXL
  return x

@cython.boundscheck(False)
def pairDist(np.ndarray[DTYPE2_t,ndim=2] r1, np.ndarray[DTYPE2_t,ndim=2] r2, double BOXL):
  cdef Py_ssize_t imax = r1.shape[0]
  cdef Py_ssize_t jmax = r2.shape[0]
  cdef Py_ssize_t i, j
  cdef double dx,dy,dz
  cdef np.ndarray[DTYPE2_t,ndim=1] r = np.zeros([imax], dtype=DTYPE2)
  
  if imax!=jmax:
    raise ValueError('First axis of r1 and r2 must be equal')

  for i in range(imax):
    dx=abs(r1[i,0] - r2[i,0])
    dy=abs(r1[i,1] - r2[i,1])
    dz=abs(r1[i,2] - r2[i,2])

    dx=dist_wrap(dx,BOXL)
    dy=dist_wrap(dy,BOXL)
    dz=dist_wrap(dz,BOXL)

    r[i] = (dx**2 + dy**2 + dz**2)**(0.5)
  
  return r

@cython.boundscheck(False)
def allDist(double[:,:] r1, double[:,:] r2, double BOXL, bint wrap):
  cdef Py_ssize_t imax = r1.shape[0]
  cdef Py_ssize_t jmax = r2.shape[0]
  cdef Py_ssize_t i, j
  cdef double dx,dy,dz
  cdef np.ndarray[DTYPE2_t,ndim=2] r = np.zeros([imax,jmax], dtype=DTYPE2)
  cdef double Bo2=BOXL/2.0

  with nogil:
    for i in prange(imax):
      for j in range(jmax):
        dx=c_fabs(r1[i,0] - r2[j,0])
        dy=c_fabs(r1[i,1] - r2[j,1])
        dz=c_fabs(r1[i,2] - r2[j,2])

        if wrap==True:
          if dx>Bo2:
            dx=dx-BOXL
          if dy>Bo2:
            dy=dy-BOXL
          if dz>Bo2:
            dz=dz-BOXL

        r[i,j] = (dx**2 + dy**2 + dz**2)**(0.5)
  
  return r

# @cython.boundscheck(False)
# def allDist(float[:,:] r1, float[:,:] r2, float BOXL, bint wrap):
#   cdef Py_ssize_t imax = r1.shape[0]
#   cdef Py_ssize_t jmax = r2.shape[0]
#   cdef Py_ssize_t i, j
#   cdef float dx,dy,dz
#   cdef np.ndarray[DTYPE3_t,ndim=2] r = np.zeros([imax,jmax], dtype=DTYPE3)
#   cdef float Bo2=BOXL/2.0
# 
#   with nogil:
#     for i in prange(imax):
#       for j in range(jmax):
#         dx=c_fabs(r1[i,0] - r2[j,0])
#         dy=c_fabs(r1[i,1] - r2[j,1])
#         dz=c_fabs(r1[i,2] - r2[j,2])
# 
#         if wrap==True:
#           if dx>Bo2:
#             dx=dx-BOXL
#           if dy>Bo2:
#             dy=dy-BOXL
#           if dz>Bo2:
#             dz=dz-BOXL
# 
#         r[i,j] = (dx**2 + dy**2 + dz**2)**(0.5)
#   
#   return r

@cython.boundscheck(False)
def selfDist(double[:,:] r, double BOXL, bint wrap):
  cdef int imax = r.shape[0]-1
  cdef int jmax = r.shape[0]
  cdef Py_ssize_t i, j
  cdef DTYPE2_t dx,dy,dz
  cdef np.ndarray[DTYPE2_t,ndim=2] dist = np.ones([jmax,jmax], dtype=DTYPE2)*-1
  cdef double Bo2=BOXL/2.0

  with nogil:
    for i in prange(imax):
      for j in range(i+1,jmax):
        dx=c_fabs(r[i,0] - r[j,0])
        dy=c_fabs(r[i,1] - r[j,1])
        dz=c_fabs(r[i,2] - r[j,2])

        if wrap:
          if dx>Bo2:
            dx=dx-BOXL
          if dy>Bo2:
            dy=dy-BOXL
          if dz>Bo2:
            dz=dz-BOXL

        dist[i,j] = (dx**2 + dy**2 + dz**2)**(0.5)
  
  return dist

@cython.boundscheck(False)
def maxDist(double[:,:] r1, double[:,:] r2, double BOXL, bint wrap):
  cdef Py_ssize_t imax = r1.shape[0]
  cdef Py_ssize_t jmax = r2.shape[0]
  cdef Py_ssize_t i, j,ival,jval
  cdef double dx,dy,dz,r
  cdef double output_maxDist = 1e-9

  for i in range(imax):
    for j in range(jmax):
      dx=c_fabs(r1[i,0] - r2[j,0])
      dy=c_fabs(r1[i,1] - r2[j,1])
      dz=c_fabs(r1[i,2] - r2[j,2])

      if wrap:
        dx = dist_wrap(dx,BOXL)
        dy = dist_wrap(dy,BOXL)
        dz = dist_wrap(dz,BOXL)

      r = dx*dx + dy*dy + dz*dz
      if r>output_maxDist: 
        ival=i
        jval=j
        output_maxDist=r
  return ival,jval,output_maxDist**(0.5)

@cython.boundscheck(False)
def minDist(double[:,:] r1, double[:,:] r2, double BOXL, bint wrap):
  cdef Py_ssize_t imax = r1.shape[0]
  cdef Py_ssize_t jmax = r2.shape[0]
  cdef Py_ssize_t i, j
  cdef Py_ssize_t ival,jval
  cdef double dx,dy,dz,r
  cdef double output_minDist = 1e9

  for i in range(imax):
    for j in range(jmax):
      dx=c_fabs(r1[i,0] - r2[j,0])
      dy=c_fabs(r1[i,1] - r2[j,1])
      dz=c_fabs(r1[i,2] - r2[j,2])

      if wrap:
        dx = dist_wrap(dx,BOXL)
        dy = dist_wrap(dy,BOXL)
        dz = dist_wrap(dz,BOXL)

      r = dx*dx + dy*dy + dz*dz
      if r < output_minDist: 
        ival = i
        jval = j
        output_minDist=r
  return ival, jval, output_minDist**(0.5)

@cython.boundscheck(False)
def minDistCheck(double[:,:] r1, double[:,:] r2, double[:] d1, double[:] d2, double factor):
  cdef Py_ssize_t imax = r1.shape[0]
  cdef Py_ssize_t jmax = r2.shape[0]
  cdef Py_ssize_t i, j
  cdef Py_ssize_t ival,jval
  cdef double dx,dy,dz,r
  cdef double minDist
  cdef bint tooClose=False

  for i in range(imax):
    for j in range(jmax):
      dx=c_fabs(r1[i,0] - r2[j,0])
      dy=c_fabs(r1[i,1] - r2[j,1])
      dz=c_fabs(r1[i,2] - r2[j,2])

      r = c_sqrt(dx*dx + dy*dy + dz*dz)
      minDist = factor*(d1[i] + d2[j])/2.0
      if r <= minDist: 
        tooClose=True
        break
    if tooClose:
      break
  return tooClose

