#include "TrilinosInterface.hpp"

#ifdef WITH_TRILINOS

/*
  Epetra_MpiComm * Ept_comm; // Epetra communicator
  Epetra_Map * RowMap; // Map defining numbering and distribution of matrix rows
  
  // Define Epetra matrix and vectors to solve Ax=b
  Epetra_CrsMatrix * Ept_A;
  Epetra_Vector * Ept_x;
  Epetra_Vector * Ept_b;
  
  int nRow; // number of rows owned by calling processor
  int* nEntperRow; // number of non-zero Entries per Row 
  
  //global variables
  int* Row_global; // Array of Row global indices (of length nRow on calling processor)
  int* cnono_v_global; // Array of global indices for non-zero Entries
*/

void TrilinosSolver::setup(int ncv, int* icv_global, int* ncols_i, int* icv_cols, MPI_Comm comm)
{
  // create Communicator object and find first (global) cv of each processor
  g_comm = new Epetra_MpiComm(comm);

 // RCP<const Teuchos::Comm<int> > m_comm (new MpiComm<int> comm);

  //const int myRank == m_comm->getRank();
  //const int numProcs = m_comm->getSize();
 
//  if (myRank == 0) {
//    cout << "Total number of processes: " << numProcs << " " << Tpetra::version() <<  endl;
//  }

  g_ncv = ncv;
  // create global cv index and Map object
  g_icv_global = new int[ncv];
  for (int i = 0; i < ncv; ++ i) {
    //g_icv_global[i] = i + first_cv;
    g_icv_global[i] = icv_global[i];
    //  if (mpi_rank == 0 ) cout << "g_icv_global = " << i << " " << g_icv_global[i] << endl;
  }
  g_map = new Epetra_Map(-1, ncv, g_icv_global, 0, *g_comm);
  
  // copy column index information and create Matrix object
  int j = 0;
  g_ncols = new int[ncv];
  for (int i = 0; i < ncv; ++ i) {
    g_ncols[i] = ncols_i[i+1]-ncols_i[i];
    j += g_ncols[i];
  }
  g_A = new Epetra_CrsMatrix(Copy, *g_map, g_ncols);
    
  // insert 0 into matrix
  g_icv_cols = new int[j];
  j = 0;
  for (int i = 0; i < ncv; ++ i) {
    int icv = g_icv_global[i];
    double values[g_ncols[i]];
    int j0 = j;
    while (j < j0 + g_ncols[i]) {
      g_icv_cols[j] = icv_cols[j];
      values[j - j0] = 0.0;
      ++ j;
    }
    g_A->InsertGlobalValues(icv, g_ncols[i], values, g_icv_cols + j0);
  }
  g_A->FillComplete();
    
  // create Vector objects
  g_x = new Epetra_Vector(*g_map);
  g_b = new Epetra_Vector(*g_map);
    

}

void TrilinosSolver::setMatrixGlobalValues(double * values_A) { 

   // set matrix and vector values
  int j = 0;
  for (int i = 0; i < g_ncv; ++ i) {
    int icv = g_icv_global[i];
    g_A->ReplaceGlobalValues(icv, g_ncols[i], values_A + j,
			     g_icv_cols + j);
    j += g_ncols[i];
  }
}


void TrilinosSolver::updateMLPreconditioner() { 
  // 
  // should be called AFTER setMatrixGlobalValues
  // 
  assert( MLPrec != NULL ) ; 
  MLPrec->ReComputePreconditioner() ; 
} 

  
int TrilinosSolver::solve_first(double* values_x, double* values_b, double* values_A,
                                 int option, int max_iters_, int max_levels, int repart_level, double tol_)
{
  int niters = -1 ; 
    
    
  max_iters = max_iters_;
  tol = tol_;

  setMatrixGlobalValues(values_A) ; 
  
  g_x->ReplaceGlobalValues(g_ncv, values_x, g_icv_global);
  g_b->ReplaceGlobalValues(g_ncv, values_b, g_icv_global);
    
  // create LinearProblem and Solver objects, set default options
  axb = new Epetra_LinearProblem(g_A, g_x, g_b);
  solver = new AztecOO(*axb);
    
  int az_options[AZ_OPTIONS_SIZE];
  double az_params[AZ_PARAMS_SIZE];
  AZ_defaults(az_options, az_params);
  //az_options[AZ_output] = AZ_summary;
  az_options[AZ_output] = AZ_none;
  az_options[AZ_conv] = AZ_rhs;
  solver->SetAllAztecOptions(az_options);
  solver->SetAllAztecParams(az_params);
    
  // set options of the linear solver
  if (option == 0) {
    // BiCGStab with Point Jacobi preconditioner
    solver->SetAztecOption(AZ_solver, AZ_bicgstab);
    solver->SetAztecOption(AZ_precond, AZ_Jacobi);
    solver->SetAztecOption(AZ_poly_ord, 1);
  }
  else if (option == 1) {
    // BiCGStab with Step-2 Block Jacobi preconditioner
    solver->SetAztecOption(AZ_solver, AZ_bicgstab);
    solver->SetAztecOption(AZ_precond, AZ_Jacobi);
    solver->SetAztecOption(AZ_poly_ord, 2);
  }
  else if (option == 2) {
    // BiCGStab with ICC(0) domain decomposition preconditioner
    solver->SetAztecOption(AZ_solver, AZ_bicgstab);
    solver->SetAztecOption(AZ_precond, AZ_dom_decomp);
    solver->SetAztecOption(AZ_subdomain_solve, AZ_icc);
  }
  else if (option == 3) {
    // PCG with ICC(0) domain decomposition preconditioner
    solver->SetAztecOption(AZ_solver, AZ_cg_condnum);
    solver->SetAztecOption(AZ_precond, AZ_dom_decomp);
    solver->SetAztecOption(AZ_subdomain_solve, AZ_icc);
  }
  else if ( option == 4 ) { 
    // 
    // use multi-level preconditioner ... 
    // 
    ParameterList MLList ; 
    ML_Epetra::SetDefaults("DD-ML",MLList); 
    MLList.set("aggregation: type", "Uncoupled-MIS");
    MLList.set("ML output", 0); 
    //MLList.set("max levels", 6); 
    MLList.set("max levels", max_levels) ; 
    MLList.set("increasing or decreasing", "increasing"); // finest level=0 
    //MLList.set("smoother: type", "ILUT"); 
    MLList.set("smoother: type", "Gauss-Seidel"); 
    MLList.set("smoother: sweeps", 3); 
    //MLList.set("coarse: type", "Amesos-Superlu") ; 
    MLList.set("coarse: type", "Gauss-Seidel");  //iterative on coarsest level .. 
    //will also support superLU for direct solve on coarse level 
    //XXX user-defined bicgstab coming .. 
    MLList.set("analyze memory", true); 
    MLList.set("coarse: max size", 32) ;  // keep coarsening until max size is reached (or max levels...) 
      
    // 
    // attempt repartitioning of the ML preconditioner 
    // to avoid communication bottlenecks ... 
    // 
    MLList.set("repartition: enable", 0); // allow repartition of the ML precond
    MLList.set("repartition: start level", repart_level); // start on level =2 by default 
    MLList.set("repartition: partitioner", "ParMETIS"); 
      
    MLPrec = new MultiLevelPreconditioner(*g_A,MLList); 
    solver->SetPrecOperator(MLPrec) ; 
    //solver->SetAztecOption(AZ_solver,AZ_bicgstab); 
    solver->SetAztecOption(AZ_solver,AZ_gmres); 
    solver->SetAztecOption(AZ_poly_ord,1); 
  }
  else if ( option == 5 ) { 
    
    solver->SetAztecOption(AZ_solver,AZ_gmres); 
    solver->SetAztecOption(AZ_kspace,max_levels) ; 
    solver->SetAztecOption(AZ_conv, AZ_noscaled);
    solver->SetAztecOption(AZ_output,AZ_all); 
    //solver->SetAztecOption(AZ_precond,AZ_none); 
    
    solver->SetAztecOption(AZ_precond,AZ_Jacobi) ;
    solver->SetAztecOption(AZ_poly_ord,1) ; 
    
    //solver->SetAztecOption(AZ_precond, AZ_dom_decomp);
    //solver->SetAztecOption(AZ_subdomain_solve, AZ_ilut);
    //solver->SetAztecParam(AZ_ilut_fill,2.0) ;
    //solver->SetAztecOption(AZ_diagnostics,AZ_none) ; 

  }
  else if ( option == 6 ) { 
    // 
    // use multi-level preconditioner ... 
    // 
    ParameterList MLList ; 
    ML_Epetra::SetDefaults("NSSA",MLList); 
    //MLList.set("aggregation: type", "Uncoupled-MIS");
    MLList.set("ML output", 0); 
    MLList.set("max levels", max_levels) ; 
    //MLList.set("repartition: enable",0);
    //MLList.set("coarse: max size", 16) ;  // keep coarsening until max size is reached (or max levels...) 
    //MLList.set("coarse: type", "Amesos-MUMPS");  //iterative on coarsest level .. 


    /*
    ParameterList testList(MLList); 
    testList.set("ML validate parameter list", false); 
    testList.set("test: Aztec",false);
    testList.set("test: Jacobi", true); 
    testList.set("test: Gauss-Seidel",true) ; 
    testList.set("test: symmetric Gauss-Seidel",false); 
    testList.set("test: block Gauss-Seidel",false); 
    testList.set("test: IFPACK", false); 
    testList.set("test: Aztec as solver", true); 
    testList.set("test: ParaSails", false); 
    testList.set("test: ML self smoother", false);
    */

    MLPrec = new MultiLevelPreconditioner(*g_A,MLList);
    /*
    MLPrec->TestSmoothers(testList) ; 
    throw(0); 
    */
    // XXX hack.. 
    const int kspace_size = 50; 
    solver->SetPrecOperator(MLPrec) ; 
    solver->SetAztecOption(AZ_kspace,kspace_size) ; 
    solver->SetAztecOption(AZ_solver,AZ_gmres);
  }
  else {
    assert(0);
  }
    
  // solve the linear system and get the solution
  niters = solver->Iterate(max_iters, tol);
  g_x->ExtractCopy(values_x);
    
  // // check Linf residual
  // g_A->Multiply(false, *g_x, *g_b);
  // g_b->Scale(-1.0);
  // g_b->SumIntoGlobalValues(g_ncv, values_b, g_icv_global);
  // double infNorm;
  // g_b->NormInf(&infNorm);
  // if (g_comm->MyPID() == 0) {
  //     cout << "InfNorm from Trilinos: " << infNorm << endl;
  // }

  if ( niters >= 0 ) // Aztec made it okay.. 
    niters = solver->NumIters() ; 
    
  
  return niters ; 
    
}


int TrilinosSolver::solve_again(double* values_x, double* values_b, bool verbose ) {

  int niters = -1 ; 
  if (verbose) {
    solver->SetAztecOption(AZ_output,1);
  }

  g_x->ReplaceGlobalValues(g_ncv, values_x, g_icv_global);
  g_b->ReplaceGlobalValues(g_ncv, values_b, g_icv_global);
    
  niters = solver->Iterate(max_iters, tol);
  g_x->ExtractCopy(values_x);

  if ( niters >= 0 ) // Aztec made it okay .. 
    niters = solver->NumIters() ; 

  if (verbose) { // return back until requested
    solver->SetAztecOption(AZ_output,AZ_none);
  }
  return niters; 
    
}

#else


void TrilinosSolver::setup(int ncv, int* icv_global, int* ncols_i, int* icv_cols, MPI_Comm comm) {
  
  std::cerr << "Error: to use trilinos you need to compile WITH_TRILINOS." << std::endl;
  throw(0);
  
}
  
int TrilinosSolver::solve_first( double * values_x, double * values_b, double * values_A, 
                                 int option, int max_iters, int max_levels, int repart_level, double tol) {
  
  std::cerr << "Error: to use trilinos you need to compile WITH_TRILINOS." << std::endl;
  throw(0);
  
}

int TrilinosSolver::solve_again(double* values_x, double* values_b, bool verbose) {
  
  std::cerr << "Error: to use trilinos you need to compile WITH_TRILINOS." << std::endl;
  throw(0);
  
}

void TrilinosSolver::setMatrixGlobalValues( double * values_A) {

  std::cerr << "Error: to use trilinos you need to compile WITH_TRILINOS." << std::endl;
  throw(0);
  
}

void TrilinosSolver::updateMLPreconditioner() {
  
  std::cerr << "Error: to use trilinos you need to compile WITH_TRILINOS." << std::endl;
  throw(0);
  
}

#endif
  


