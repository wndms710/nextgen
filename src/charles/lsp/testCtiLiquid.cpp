#include <iostream>

#include "CTI.hpp"
using namespace CTI;

#include "CtiLiquid.hpp"  


int main(int argc, char * argv[]) {
  
  try {
    
    CTI_Init(argc,argv,"testCtiLiquid.in");

    {
      Param * param = getParam("LSP.MATERIAL");
      if (param == NULL) CERR("could not find param LSP.MATERIAL");

      int iarg = 0;
      string token = param->getString(iarg++);
      if (token == "LIQUID") {
        CtiLiquid * liquid = newCtiLiquid(param,iarg);
        liquid->info();

        double T0 = getDoubleParam("T0",200);
        double T1 = getDoubleParam("T1",700);
        int n = getIntParam("N",200);
        liquid->setInfoPrefix("XXXX");
        liquid->info(T0,T1,n);

      } else {
        CERR("could not find param LSP.MATERIAL LIQUID");
      }
    }
    
    CTI_Finalize();
    
  }
  catch (int e) {
    if (e == 0)
      CTI_Finalize();
    else
      CTI_Abort();
  }
  catch(...) {
    CTI_Abort();
  }
  
  return(0);
  
}


