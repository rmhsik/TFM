# distutils: language=c++
# cython: language_level=3
# cython: boundscheck=False

import numpy as np
from cython.parallel import prange

def MatrixMult(double complex[:,:] a,
               double complex[:,:] b,
               int N):
    return np.asarray(cMatrixMult(a,b,N))

cdef double complex[:,:] cMatrixMult(double complex [:,:] a,
                                     double complex [:,:] b,
                                     int N):
    cdef int i,j
    cdef double complex [:,:] out
    
    out = np.zeros((N,N)).astype('complex128')
    
    for i in prange(N,nogil=True):
        for j in range(N):
            out[i,j] = a[i,j]*b[i,j]
    return out