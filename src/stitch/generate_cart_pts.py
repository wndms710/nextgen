#!/usr/bin/env python 

import numpy as np

def generate_cart_1d(n):

  x = np.zeros(n)
  for i in range(0,n):
    x[i] = (float(i) + 0.5)/float(n)

  return x

nx = 8
ny = 8 
nz = 8

x = generate_cart_1d(nx)
y = generate_cart_1d(ny)
z = generate_cart_1d(nz)

for k in range(0,nz):
  for j in range(0,ny):
    for i in range(0,nx):

      print "%12.8f  %12.8f   %12.8f" % (x[i], y[j], z[k])




