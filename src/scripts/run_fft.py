#!/usr/bin/env python 

import numpy as np
import matplotlib.pyplot as plt
import fft_coeffs as fc
import argparse
import sys

#------------------------------------------------------------------------------

# get args from user
parser = argparse.ArgumentParser()
parser.add_argument('--infile',  required=True,      type=str,   help='input data file')
parser.add_argument('--col_x',   default=1,          type=int,   help='0-indexed column number for time')
parser.add_argument('--col_y',   default=3,          type=int,   help='0-indexed column number for pressure')
parser.add_argument('--t_start', default=-1.0,       type=float, help='discard data with t<t_start')
parser.add_argument('--t_chunk', default=0.25,       type=float, help='chunk size (set T_CHUNK<=0 to disable chunking)')
parser.add_argument('--overlap', default=0.5,        type=float, help='chunk overlap')
parser.add_argument('--p_scale', default=.000145038, type=float, help='factor to convert p units (default: pa->psi)')
args = parser.parse_args()

print('input params:')
print('  infile  = %s' % args.infile )
print('  col_x   = %i' % args.col_x  )
print('  col_y   = %i' % args.col_y  )
print('  t_start = %g' % args.t_start)
if (args.t_chunk > 0.0):
  print('  t_chunk = %g' % args.t_chunk)
  print('  overlap = %g' % args.overlap)
print('')

# load data
t_raw,p_raw = np.loadtxt(args.infile, usecols=(args.col_x,args.col_y), unpack=True)
p_raw = p_raw * args.p_scale # convert p units

# discard out-of-order data
nn=len(t_raw) 
t=np.zeros(nn)
p=np.zeros(nn)
t[0]=t_raw[-1]
p[0]=p_raw[-1]
n = 1
for i in np.arange(nn-2,-1,-1):
  if t_raw[i]<t[n-1]:
    n += 1
    t[n-1]=t_raw[i]
    p[n-1]=p_raw[i] 
t=np.flipud(t[:n])
p=np.flipud(p[:n])

ii = np.argmin(np.abs(t-args.t_start))
t  = t[ii:]
p  = p[ii:]
n  = len(t)

print('original signal: n = %i, t = [%g:%g], p = [%g:%g]' % (len(t_raw),t_raw[0],t_raw[-1],min(p_raw),max(p_raw)))
print('trimmed  signal: n = %i, t = [%g:%g], p = [%g:%g]\n' % (n,t[0],t[-1],min(p),max(p)))

if (args.t_chunk <= 0.0): # don't chunk data

  freq,psd,phase = fc.analyze_fft(t,p)
  amp = fc.convertPsdToPktoPk(freq,psd)

  fc.writeTwoColFile(freq,psd,'psd.dat')
  #fc.writeTwoColFile(freq,phase,'phase.dat')
  fc.writeTwoColFile(freq,amp,'amp.dat')

else: # chunk data

  dt      = (t[-1]-t[0])/(n-1)
  n_chunk = int(args.t_chunk/dt) # samples per chunk

  istart  = 0
  ichunk  = 0
  freq    = None
  psd     = None
  first   = True

  while (istart < n): 

    # discard the last bit of the signal... 
    if (istart + n_chunk > n):
      print('done... discarding %i samples' % (n-istart))
      break

    t1,y1 = fc.chunk(t,p,istart,istart+n_chunk)
    if (len(t1) != n_chunk): 
      print('incorrect chunk length = %i' % len(t1))
      raise Error
    print('processing chunk %i: idx = [%i:%i], t = [%g:%g]' % (ichunk,istart,istart+n_chunk,t1[0],t1[-1]))

    f1,psd1,phase = fc.analyze_fft(t1,y1)
    if (first):
      first = False
      freq = f1
      psd = psd1
    else: 
      psd = psd + psd1

    istart = istart + int(n_chunk*(1.0-args.overlap))
    ichunk = ichunk + 1

  psd = psd/float(ichunk)
  fc.writeTwoColFile(freq,psd,'psd.dat')

  amp = fc.convertPsdToPktoPk(freq,psd)
  fc.writeTwoColFile(freq,amp,'amp.dat')

#------------------------------------------------------------------------------
# plot results

plt.figure(1)
plt.plot(t_raw,p_raw,t,p)
plt.xlabel('time [sec]')
plt.ylabel('pressure')

plt.figure(2)
plt.plot(freq,amp)
#plt.loglog(freq,amp)
plt.xlabel('frequency [Hz]')
plt.ylabel('amplitude')
plt.xlim([0,2000])

plt.figure(3)
plt.loglog(freq,np.real(psd))
plt.xlabel('frequency [Hz]')
plt.ylabel('PSD [p$^2$/Hz]')

plt.show()
