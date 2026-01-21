
import numpy as np
import pylab as pyl

def analyze_fft(t,p,fname=None) : 

  dt_samp = t[1] - t[0] 
  print("Dt_samp = ", dt_samp)

  n       = len(p)
  print("n = ", n)
  p_prime = p
  w       = np.hanning(n)
  p_prime = p_prime - np.mean(p_prime)
  print("mean = ", np.mean(p_prime))

  p_prime = p_prime*w
  phat    = np.fft.fft(p_prime) /float(n)
  freq    = np.fft.fftfreq(n,dt_samp)

  df = freq[1]-freq[0]
  print("freq bin = " ,df)

# check our statement of parseval...
  sum_t = 0.0
  for i in range(0,len(p_prime)):
    sum_t = sum_t + p_prime[i]*p_prime[i]
  print("1/N sum t p^2 = ", sum_t/float(len(p_prime)))

  sum_w = 0.0
  for i in range(0,len(phat)):
    sum_w = sum_w + np.real(phat[i]*np.conj(phat[i]))
  print("sum w p^2 = ", sum_w)

# one sided fft (only considering positive freq..)  
  one_sided_correction = 2.0
# energy correction for the hanning window...
  energy_correction = 8.0/3.0
  psd  = phat * np.conj(phat) * energy_correction*one_sided_correction/df
  return freq[1:int(n/2)], psd[1:int(n/2)]

def chunk(t,y,i1,i2): 
  return t[i1:i2], y[i1:i2]

def writePsdFile(f1,psd1,fname) : 
  nn = len(f1)
  tmp = np.zeros([nn,2])
  tmp[:,0] = f1
  tmp[:,1] = psd1
  np.savetxt(fname,tmp)


