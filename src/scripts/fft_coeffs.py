import numpy as np

def analyze_fft(t,y): 

  print('analyze_fft()...')
  n = len(y)
  print('  n = %i' % n)
  dt_samp = (t[-1]-t[0])/(n-1)
  print('  dt_samp = %g' % dt_samp)

  y_prime = y
  w       = np.hanning(n)
  y_prime = y_prime - np.mean(y_prime)
  print('  mean = %g' % np.mean(y_prime))

  y_prime = y_prime*w
  yhat    = np.fft.fft(y_prime)/float(n)
  freq    = np.fft.fftfreq(n,dt_samp)

  df = freq[1]-freq[0]
  print('  freq bin = %g' % df)

  # check our statement of parseval...
  sum_t = 0.0
  for i in range(0,len(y_prime)):
    sum_t = sum_t + y_prime[i]*y_prime[i]
  print('  1/N sum t y^2 = %g' % (sum_t/float(len(y_prime))))

  sum_w = 0.0
  for i in range(0,len(yhat)):
    sum_w = sum_w + np.real(yhat[i]*np.conj(yhat[i]))
  print('  sum w y^2 = %g' % sum_w)

  # one sided fft (only considering positive freq..)  
  one_sided_correction = 2.0
  # energy correction for the hanning window...
  energy_correction = 8.0/3.0
  psd  = yhat * np.conj(yhat) * energy_correction*one_sided_correction/df
  phase = np.angle(yhat)

  print('')
  return freq[1:int(n/2)], psd[1:int(n/2)], phase[1:int(n/2)]

def chunk(t,y,i1,i2): 
  return t[i1:i2], y[i1:i2]

def writeTwoColFile(f,var,fname) : 
  tmp = np.zeros([len(f),2])
  tmp[:,0] = f
  tmp[:,1] = np.real(var)
  np.savetxt(fname,tmp)

def convertPsdToPktoPk(freq,psd):
  df = freq[1]-freq[0]
  pk2pk = np.zeros(len(psd))
  for i in range(0,len(psd)):
    E_om   = np.real(psd[i]*df)
    pk2pk[i] = 2.0*np.sqrt(2.0*E_om)
  return pk2pk 
