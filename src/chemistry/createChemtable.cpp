// chemtable generation
#include "createChemtable.hpp"


int main(int argc,char * argv[]) {
  
  try {
    
    CTI_Init(argc,argv, "createChemtable.in"); 
    cout << "Program: createChemtable:" << endl;
    
    // pointer to the abstract class
    AbstractChemtable * Table;
    
    string chemtabletype = getStringParam("CHEMTABLE.TYPE");
    
    if ( chemtabletype == "VIDA_PREMIXED_FPV_CART1D" ) 
      Table = new VidaChemtablePFPVCart1D(chemtabletype);
    else if ( chemtabletype == "VIDA_PREMIXED_FPV_CART2D" ) 
      Table = new VidaChemtablePFPVCart2D(chemtabletype);
    else if ( chemtabletype == "CHARLES_PREMIXED_FPV_CART2D" || chemtabletype == "PREMIXED") 
      Table = new CharlesChemtablePFPVCart2D(chemtabletype);
    else if ( chemtabletype == "VIDA_NON-PREMIXED_FPV_CART3D" ) 
      Table = new VidaChemtableNFPVCart3D(chemtabletype);
    else if ( chemtabletype == "CHRIS_NON-PREMIXED_FPV_CART3D" || chemtabletype == "NONPREMIXED" || chemtabletype == "DIFFUSION")
      Table = new ChrisChemtableNFPVCart3D(chemtabletype);
    else if ( chemtabletype == "VIDA_NON-PREMIXED_FPV_KDT3D" ) 
      Table = new VidaChemtableNFPVKdt3D(chemtabletype);
    else if ( chemtabletype == "CHRIS_NON-PREMIXED_FPV_KDT3D" ) 
      Table = new ChrisChemtableNFPVKdt3D(chemtabletype);
    else if ( chemtabletype == "VIDA_NON-PREMIXED_FPVT_CART4D" ) 
      CERR_S("VIDA_NON-PREMIXED_FPVH_CART4D is not implemented yet")
    else if ( chemtabletype == "VIDA_NON-PREMIXED_FPVH_CART4D" ) 
      CERR_S("VIDA_NON-PREMIXED_FPVH_CART4D is not implemented yet")
    else if ( chemtabletype == "VIDA_NON-PREMIXED_FPVH_KDT4D" )
      CERR_S("VIDA_NON-PREMIXED_FPVH_KDT4D is not implemented yet")
    else if ( chemtabletype == "VIDA_NON-PREMIXED_FPV_KDT5D" ) 
      CERR_S("VIDA_NON-PREMIXED_FPV_KDT5D is not implemented yet")
    else{
      cerr << "unknown CHEMTABLE.TYPE: " << chemtabletype << endl;
      throw(0);
    }

    // initialize, build, and write the table
    Table->init();
    Table->build();
    Table->finalize();
     
    // delete the table
    delete Table; 

    dumpParams();
    dumpParamUsage();
    
  }
  
  catch (int e) {
    if (e == 0) {
      dumpParams();
      dumpParamUsage();
      return(-1);
    }
    else
      return(-1);
  }
  catch(...) {
    return(-1);
  }
  
  return(0);
}


// TODO:
// 1) catch duplicate CHEMTABLE.SPECIES before they hit asserts later on....
// 2) changes CHEMTABLE.SPECIES to CHEMTABLE.VARS to support other types of flamelet data as well as matching solver param syntax
