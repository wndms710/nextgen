#include "CTI.hpp"
using namespace CTI;

// make this part of CTI eventually...

#include "CVM.hpp"
using namespace CVM;

#include "SimpleSurface.hpp"
#include "CtiCanvas.hpp"

void writeDiffImage(SimpleSurface& ss1,SimpleSurface& ss2) {

  // position the pixels...

  double rmax1 = ss1.getBoundingBoxRmax();
  double xc1[3]; ss1.getBoundingBoxCenter(xc1);
  
  double rmax2 = ss2.getBoundingBoxRmax();
  double xc2[3]; ss2.getBoundingBoxCenter(xc2);

  double target[3]; FOR_I3 target[i] = 0.5*(xc1[i]+xc2[i]);
  double width = max(rmax1,rmax2);
  
  double camera[3];
  camera[0] = target[0] + width*0.75;
  camera[1] = target[1] + width*sqrt(3.0)*0.25;
  camera[2] = target[2] + width*0.5;
 
  double up[3];
  up[0] = 0.0;
  up[1] = 0.0;
  up[2] = 1.0;
  
  int size[2] = { 2000, 2000 };
  
  CtiCanvas canvas(target,camera,up,width*2.0,size,NULL);

  for (int ist = 0; ist < ss1.nst; ++ist) {
    canvas.addTriVolVis(ss1.xsp[ss1.spost[ist][0]],ss1.xsp[ss1.spost[ist][1]],ss1.xsp[ss1.spost[ist][2]],1);
  }
  
  for (int ist = 0; ist < ss2.nst; ++ist) {
    canvas.addTriVolVis(ss2.xsp[ss2.spost[ist][0]],ss2.xsp[ss2.spost[ist][1]],ss2.xsp[ss2.spost[ist][2]],2);
  }
  
  canvas.writeImageVolVis("diffgeom.png");
  
}

int main(int argc, char * argv[]) {

  try {

    CTI_Init(argc,argv,"geomdiff.in");

    {

      assert(cvm == NULL);
      cvm = new CtiVarMachine(NULL); // the arg is a solver, but this can be skipped

      {
	
	Param * param1 = getParam("SURF1");
	Param * param2 = getParam("SURF2");
	if (!param1 || !param2) {
	  CERR("missing SURF1 and/or SURF2 param. Usage:\n\n" <<
	       " > diffgeom.exe --SURF1 SBIN <surf1.sbin> --SURF2 SBIN <surf2.sbin>\n");
	}
	
	SimpleSurface ss1;
	ss1.init(param1);
	
	SimpleSurface ss2;
	ss2.init(param2);
	
	writeDiffImage(ss1,ss2);
		
	
      }

      delete cvm;

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
