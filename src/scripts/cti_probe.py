#!/usr/bin/env python 
# 
# in order to use this, you will need 
# to set your PYTHONPATH appropriately
#  export PYTHONPATH=$CTI_SOLVER_HOME/script:$PYTHONPATH
import numpy as np
import sys 
import os

class Probe: 

  def __init__(self,filename): 
    # load in the readme that stores point locations
    self.fname = filename 
    self.xp  = np.loadtxt(filename+".README") 
    self.np  = len(self.xp[:,0])
    print("number of probe points : " , self.np)
    self.pd  = [] 
    self.var = ""

  def load(self,var,pid,useCache=True): 
    # load the specific variable and probe location
    # cache the entire history of the this particular 
    # probe data unless otherwise requested 
    if ( useCache) : 
      # assuming you have memory available to load this 
      # variable into memory
      if ( self.var != var) : 
        self.var = var
        fpname   = self.fname + "." + var

        f        = open(fpname, 'r') 

        nd       = 0
        for line in f : 
          if (line[0] != '#'):
            nd = nd + 1

        print(" > Number of data lines : " , nd) 
        print(" > Approx mem usage [MB] : " , (nd*(self.np+1)*8.0)/(1024.*1024.))

        f.seek(0)

        self.pd  = []
        self.pd  = np.zeros([nd,self.np+1])
        id       = 0
        step_prev = -1
        dstep = 0
        for line in f : 

          if (line[0] == '#'):
            continue 

          l = line.split()
          #print l  

          step = int(l[0])
          if (step_prev == -1):
            step_prev = step
          elif (dstep == 0):
            dstep = step - step_prev
            step_prev = step
          else:
            if (step - step_prev != dstep):
              print("\n\n***********************\nERROR probe step spacing varies or rewinds\n**********************\n\n") 
              raise UndefinedError
            step_prev = step
                        
          if ( self.np != int(l[2])): 
            raise ProbeLengthError

          self.pd[id,0] = float(l[1]) # time
          for i in range(3,len(l)) :
            self.pd[id,i-2] = float(l[i]) # probe data
          
          id = id + 1 


        f.close() 
    else : 
      # not yet
      raise UndefinedError


    return self.pd[:,0], self.pd[:,pid]
  
   
      


