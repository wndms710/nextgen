import os
import struct
import numpy as np

def readpcp(pcd):
  le = '<'    # little endian
  be = '>'    # big endian
  fmt_l = 'q' # long long
  fmt_d = 'd' # double
  fmt_f = 'f' # float
  size_l = struct.calcsize(fmt_l)
  size_f = struct.calcsize(fmt_f)
  size_d = struct.calcsize(fmt_d)

  print('reading pcd = %s ...' % pcd)
  pbin = ''
  for s in pcd.split('.')[:-2]: pbin += s+'.'
  pbin += 'pbin'
  print(' > pbin = %s' % pbin)

  # read pbin...
  f = open(pbin,'rb')
  magic = struct.unpack(le+fmt_l,f.read(size_l))[0]
  if magic == 1235813:
    fmt_l = le+fmt_l
    fmt_d = le+fmt_d
    fmt_f = le+fmt_f
  else:
    fmt_l = be+fmt_l
    fmt_d = be+fmt_d
    fmt_f = be+fmt_f
  ver  = struct.unpack(fmt_l,f.read(size_l))[0]
  npts = struct.unpack(fmt_l,f.read(size_l))[0]
  assert(struct.unpack(fmt_l,f.read(size_l))[0]==2)
  print(' > npts = %i' % npts)
  xpts = np.fromfile(f,dtype=[('x',fmt_d),('y',fmt_d),('z',fmt_d)],count=npts)
  ipts = np.fromfile(f,dtype=fmt_l,count=npts)
  f.close()

  # read pcd...
  f = open(pcd,'rb')
  f.seek(2*size_l,os.SEEK_SET) # skip endian,version... (assume same as pbin)
  assert(npts==struct.unpack(fmt_l,f.read(size_l))[0])
  nvar = struct.unpack(fmt_l,f.read(size_l))[0]
  prec = struct.unpack(fmt_l,f.read(size_l))[0]
  print(' > nvar = %i' % nvar) 
  print(' > prec = %i' % prec) 
  if prec == 0:
    myfmt = fmt_f
    mysize = size_f
  else:
    myfmt = fmt_d
    mysize = size_d

  # ??? or we could get var names from README (if present), then read with dtype=[('var',myfmt)]
  data = np.zeros((nvar,npts))
  for i in range(nvar):
    f.seek(5*size_l+i*npts*mysize, os.SEEK_SET)
    data[i] = np.fromfile(f,dtype=myfmt,count=npts)
  f.close()

  # rearrange data to same order as original probe file
  iord = np.argsort(ipts)

  return npts, xpts[iord], data[:,iord]
