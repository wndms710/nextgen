#ifndef TRILINOSINTERFACE_HPP
#define TRILINOSINTERFACE_HPP

#include <iostream>
#include "Macros.hpp"
#include "MpiStuff.hpp"
using namespace MpiStuff;

// ====================================================================
// the compiler option WITH_TRILINOS must be defined to build these...
// ====================================================================

#ifdef WITH_TRILINOS

// Trilinos
#include "AztecOO_config.h"

//#include "mpi.h"
#include "Epetra_MpiComm.h"

#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_LinearProblem.h"
#include "AztecOO.h"


#include "ml_include.h"
#include "ml_MultiLevelPreconditioner.h"
using namespace Teuchos ; 
using namespace ML_Epetra ; 

class TrilinosSolver { 
public : 
  
  Epetra_MpiComm * g_comm;
  Epetra_Map * g_map;
  Epetra_CrsMatrix * g_A;
  Epetra_Vector * g_x;
  Epetra_Vector * g_b;
  int g_ncv;
  int* g_icv_global;
  int* g_ncols;
  int* g_icv_cols;


  // 
  // solver variables ... 
  // 
  Epetra_LinearProblem * axb;
  AztecOO * solver;
  int max_iters;
  double tol;
  MultiLevelPreconditioner * MLPrec ; 
  
public : 
  TrilinosSolver() { 
    g_comm       = NULL ; 
    g_map        = NULL ; 
    g_A          = NULL ; 
    g_x          = NULL ; 
    g_b          = NULL ; 
    g_ncv        = -1; 
    
    g_icv_global = NULL ; 
    g_ncols      = NULL ; 
    g_icv_cols   = NULL ; 
    

    // solver variables .. 
    axb          = NULL ; 
    solver       = NULL ; 
    max_iters    = -1 ; 
    MLPrec       = NULL ; 
  }


  void setup( int ncv, int * icv_global, int * ncols_i , int * icv_cols, MPI_Comm comm) ; 
  int solve_once (double * values_x, double * values_b , double * values_A, 
		   int option , int max_iters, int max_levels, double tol) { 
    // update this later.. XXXX
    //solve_first(values_x,values_b,values_A,option,max_iters,max_levels,tol) ; 
    assert(0); // deprecated ... 
    return 0; 
  }

  int solve_first( double * values_x, double * values_b, double * values_A, 
                   int option, int max_iters, int max_levels, int repart_level, double tol   ) ;

  int solve_again( double * values_x, double * values_b, bool verbose=true) ; 

  void setMatrixGlobalValues( double * values_A) ; 
  void updateMLPreconditioner() ; 
 
  
  ~TrilinosSolver() { 

    if ( axb != NULL ) { 
      delete axb ; axb = NULL ; 
    }

    if ( solver != NULL) { 
      delete solver ; solver = NULL ; 
    } 

    if ( MLPrec != NULL) { 
      delete MLPrec ; MLPrec = NULL; 
    }

    if ( g_comm != NULL ) { 
      delete g_comm; g_comm = NULL ; 
    } 

    if ( g_map != NULL) { 
      delete g_map ; g_map = NULL; 
    } 

    if ( g_A != NULL) { 
      delete g_A ; g_A = NULL ; 
    }

    if ( g_x != NULL ) { 
      delete g_x ; g_x = NULL ; 
    } 

    if ( g_b != NULL ) { 
      delete g_b ; g_b = NULL ; 
    } 
    
    DELETE(g_icv_global) ; 
    DELETE(g_ncols     ) ; 
    DELETE(g_icv_cols  ) ; 


   
#if defined(BGP_MEMCHECK) || defined(BGQ_MEMCHECK) 
    if ( mpi_rank == 0 ) 
      cout << ">> Memory available after Trilinos cleanup..." << endl; 
    checkMemoryUsage(); 
#endif 
 
    
  }//destructor .. 



};

#else 

// 
// empty trilinos class 
// 
class TrilinosSolver { 
public : 

  
  void setup( int ncv, int * icv_global, int * ncols_i , int * icv_cols, MPI_Comm comm) ; 
  int solve_once (double * values_x, double * values_b , double * values_A, 
		   int option , int max_iters, int max_levels, double tol) { 
    throw(0); 
    return 0; 
  }

  int solve_first( double * values_x, double * values_b, double * values_A, 
                    int option, int max_iters, int max_levels, int repart_level, double tol) ;

  int solve_again( double * values_x, double * values_b, bool verbose=true) ; 

  void setMatrixGlobalValues( double * values_A); 
  void updateMLPreconditioner(); 

};

#endif // #ifdef WITH_TRILINOS
#endif // #ifndef TRILINOSINTERFACE_HPP
