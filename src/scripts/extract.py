import numpy as np
import pylab as pyl
import cti_probe
import fft_coeffs2 as fc
import math
import scipy.io as io

def avg_ring(pp,np,n_samp):
  # n_samp is measured from the end...
  f_avg = None
  psd_avg = None
  count   = 0
  for i in range(0,np):
    t_tmp,p_tmp = pp.load('p',i+1)
    n = len(t_tmp)
    ff,psd_tmp  = fc.analyze_fft(t_tmp[n-n_samp:],p_tmp[n-n_samp:])
    if ( count == 0) : 
      f_avg   = ff
      psd_avg = psd_tmp
    else: 
      psd_avg = psd_avg + psd_tmp
    count = count+1
  psd_avg = psd_avg/float(count)
  return f_avg,psd_avg

Pref2 = 4.0e-10 # pa
Pexp  = 101325 # pa
#Ma    = 0.9
Ma = 1.0
gamma = 1.4
pscaling = 1.0 
dBscaling = Ma*pscaling*pscaling/Pref2

def showme(case,fname,title=""): 

   rp5 = cti_probe.Probe(case)
   t6,p6 = rp5.load('p',3)
   print("start, end times: ", t6[0], "   ",t6[-1]) 
   n_samp = len(t6)
   f6,psd6 = avg_ring(rp5,36,n_samp)

   fc.writePsdFile(f6,psd6,fname)

   pyl.figure() 
   pyl.title(title)
   pyl.xlabel("Hz")
   pyl.ylabel("PSD [dB/Hz]")
   pyl.semilogx(f6,10.0*np.log10(psd6*dBscaling/math.pow(Ma,4)), '-xr')
   # set the axis limits...
   #pyl.axis([10,10000,0,130])

   oname = title+".png"
   print(oname)
   pyl.savefig(oname)

#showme('./ring/Ring1','ring1_psd.dat',title="ring1")
#showme('./ring/Ring2','ring2_psd.dat',title="ring2")
#showme('./ring/Ring3','ring3_psd.dat',title="ring3")
#showme('./ring/Ring4','ring4_psd.dat',title="ring4")
#showme('./ring/Ring5','ring5_psd.dat',title="ring5")
#showme('./ring/Ring6','ring6_psd.dat',title="ring6")

showme('./ring/NRing00','ring00_psd.dat',title="ring00")
showme('./ring/NRing01','ring01_psd.dat',title="ring01")
showme('./ring/NRing02','ring02_psd.dat',title="ring02")
showme('./ring/NRing03','ring03_psd.dat',title="ring03")
showme('./ring/NRing04','ring04_psd.dat',title="ring04")
showme('./ring/NRing05','ring05_psd.dat',title="ring05")
showme('./ring/NRing06','ring06_psd.dat',title="ring06")
showme('./ring/NRing07','ring07_psd.dat',title="ring07")
showme('./ring/NRing08','ring08_psd.dat',title="ring08")
showme('./ring/NRing09','ring09_psd.dat',title="ring09")
showme('./ring/NRing10','ring10_psd.dat',title="ring10")
showme('./ring/NRing11','ring11_psd.dat',title="ring11")
showme('./ring/NRing12','ring12_psd.dat',title="ring12")
showme('./ring/NRing13','ring13_psd.dat',title="ring13")
showme('./ring/NRing14','ring14_psd.dat',title="ring14")
showme('./ring/NRing15','ring15_psd.dat',title="ring15")
showme('./ring/NRing16','ring16_psd.dat',title="ring16")
showme('./ring/NRing17','ring17_psd.dat',title="ring17")
showme('./ring/NRing18','ring18_psd.dat',title="ring18")
showme('./ring/NRing19','ring19_psd.dat',title="ring19")
showme('./ring/NRing20','ring20_psd.dat',title="ring20")
showme('./ring/NRing21','ring21_psd.dat',title="ring21")
showme('./ring/NRing22','ring22_psd.dat',title="ring22")
showme('./ring/NRing23','ring23_psd.dat',title="ring23")
showme('./ring/NRing24','ring24_psd.dat',title="ring24")
showme('./ring/NRing25','ring25_psd.dat',title="ring25")
showme('./ring/NRing26','ring26_psd.dat',title="ring26")
showme('./ring/NRing27','ring27_psd.dat',title="ring27")
showme('./ring/NRing28','ring28_psd.dat',title="ring28")
showme('./ring/NRing29','ring29_psd.dat',title="ring29")
showme('./ring/NRing30','ring30_psd.dat',title="ring30")
showme('./ring/NRing31','ring31_psd.dat',title="ring31")
showme('./ring/NRing32','ring32_psd.dat',title="ring32")
showme('./ring/NRing33','ring33_psd.dat',title="ring33")
showme('./ring/NRing34','ring34_psd.dat',title="ring34")
showme('./ring/NRing35','ring35_psd.dat',title="ring35")
showme('./ring/NRing36','ring36_psd.dat',title="ring36")
showme('./ring/NRing37','ring37_psd.dat',title="ring37")
showme('./ring/NRing38','ring38_psd.dat',title="ring38")
showme('./ring/NRing39','ring39_psd.dat',title="ring39")
showme('./ring/NRing40','ring40_psd.dat',title="ring40")
showme('./ring/NRing41','ring41_psd.dat',title="ring41")
showme('./ring/NRing42','ring42_psd.dat',title="ring42")
showme('./ring/NRing43','ring43_psd.dat',title="ring43")
showme('./ring/NRing44','ring44_psd.dat',title="ring44")
showme('./ring/NRing45','ring45_psd.dat',title="ring45")
showme('./ring/NRing46','ring46_psd.dat',title="ring46")
showme('./ring/NRing47','ring47_psd.dat',title="ring47")
showme('./ring/NRing48','ring48_psd.dat',title="ring48")
showme('./ring/NRing49','ring49_psd.dat',title="ring49")
showme('./ring/NRing50','ring50_psd.dat',title="ring50")
showme('./ring/NRing51','ring51_psd.dat',title="ring51")
showme('./ring/NRing52','ring52_psd.dat',title="ring52")
showme('./ring/NRing53','ring53_psd.dat',title="ring53")
showme('./ring/NRing54','ring54_psd.dat',title="ring54")
showme('./ring/NRing55','ring55_psd.dat',title="ring55")
showme('./ring/NRing56','ring56_psd.dat',title="ring56")
showme('./ring/NRing57','ring57_psd.dat',title="ring57")
showme('./ring/NRing58','ring58_psd.dat',title="ring58")
showme('./ring/NRing59','ring59_psd.dat',title="ring59")
showme('./ring/NRing60','ring60_psd.dat',title="ring60")
showme('./ring/NRing61','ring61_psd.dat',title="ring61")
showme('./ring/NRing62','ring62_psd.dat',title="ring62")
showme('./ring/NRing63','ring63_psd.dat',title="ring63")
showme('./ring/NRing64','ring64_psd.dat',title="ring64")
showme('./ring/NRing65','ring65_psd.dat',title="ring65")
showme('./ring/NRing66','ring66_psd.dat',title="ring66")
showme('./ring/NRing67','ring67_psd.dat',title="ring67")
showme('./ring/NRing68','ring68_psd.dat',title="ring68")
showme('./ring/NRing69','ring69_psd.dat',title="ring69")
showme('./ring/NRing70','ring70_psd.dat',title="ring70")
showme('./ring/NRing71','ring71_psd.dat',title="ring71")
showme('./ring/NRing72','ring72_psd.dat',title="ring72")
showme('./ring/NRing73','ring73_psd.dat',title="ring73")
showme('./ring/NRing74','ring74_psd.dat',title="ring74")
showme('./ring/NRing75','ring75_psd.dat',title="ring75")
showme('./ring/NRing76','ring76_psd.dat',title="ring76")
showme('./ring/NRing77','ring77_psd.dat',title="ring77")
showme('./ring/NRing78','ring78_psd.dat',title="ring78")
showme('./ring/NRing79','ring79_psd.dat',title="ring79")
showme('./ring/NRing80','ring80_psd.dat',title="ring80")
showme('./ring/NRing81','ring81_psd.dat',title="ring81")
showme('./ring/NRing82','ring82_psd.dat',title="ring82")
showme('./ring/NRing83','ring83_psd.dat',title="ring83")
showme('./ring/NRing84','ring84_psd.dat',title="ring84")
showme('./ring/NRing85','ring85_psd.dat',title="ring85")
showme('./ring/NRing86','ring86_psd.dat',title="ring86")
showme('./ring/NRing87','ring87_psd.dat',title="ring87")
showme('./ring/NRing88','ring88_psd.dat',title="ring88")
showme('./ring/NRing89','ring89_psd.dat',title="ring89")
showme('./ring/NRing90','ring90_psd.dat',title="ring90")
showme('./ring/NRing91','ring91_psd.dat',title="ring91")
showme('./ring/NRing92','ring92_psd.dat',title="ring92")
showme('./ring/NRing93','ring93_psd.dat',title="ring93")
showme('./ring/NRing94','ring94_psd.dat',title="ring94")
showme('./ring/NRing95','ring95_psd.dat',title="ring95")
showme('./ring/NRing96','ring96_psd.dat',title="ring96")
showme('./ring/NRing97','ring97_psd.dat',title="ring97")
showme('./ring/NRing98','ring98_psd.dat',title="ring98")
showme('./ring/NRing99','ring99_psd.dat',title="ring99")

# uncomment to see everything
#pyl.show()

