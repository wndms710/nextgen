#include "Chemtable.hpp"

// ===================
// delete chemtable 1D
// ===================
void deleteChemtable(AbstractChemtable1D * &chemtable) {
  if ( chemtable != NULL ) delete chemtable;
}

// ===================
// delete chemtable 4D
// ===================
void deleteChemtable(AbstractChemtable4D * &chemtable) {
  if ( chemtable != NULL ) delete chemtable;
}




    
