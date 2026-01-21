#include "CTI.hpp"
#include "Macros.hpp"
#include "MpiStuff.hpp"
#include <complex>
#include "observer.hpp"
#include "FFT.hpp"
#include "ByteSwap.hpp"
#include <vector>
#include "TransposeOperator.hpp"
using namespace CTI;

// ==================================================================
//
// The Face contains all the information about the FW-H input surface:
//   - location, normal vector and panel area
//   - flow field data
// It also contains the computed source terms
//
// ==================================================================

class Face {
    
public:
  
  double x,y,z;     // position
  double nx,ny,nz;  // normal vector
  
  // these variables get dimension nt (i.e. # of frames or snapshots) 
  // To reduce memory usage, rho and area not stored anymore 
  double * sourceF1;
  double * sourceF2;
  double * sourceF3;
  double * sourceQ;
  // these variables get dimension nun_freq (i.e. # of frequencies of interest)
  complex<double> * sourceF1hat;
  complex<double> * sourceF2hat;
  complex<double> * sourceF3hat;
  complex<double> * sourceQhat;

  // ----------------------------
  // Zonal Noise Analysis (ZNA)
  int ZNAsize;     // number of zone(s) containing this face
  int * ZNAindex;  // indices of the zone(s) containing this face 

  Face() {
    sourceF1 = NULL;
    sourceF2 = NULL;
    sourceF3 = NULL;
    sourceQ = NULL;
    sourceF1hat = NULL;
    sourceF2hat = NULL;
    sourceF3hat = NULL;
    sourceQhat = NULL;
    ZNAindex = NULL;
  }

  ~Face() {
    DELETE(sourceF1);
    DELETE(sourceF2);
    DELETE(sourceF3);
    DELETE(sourceQ);
    DELETE(sourceF1hat);
    DELETE(sourceF2hat);
    DELETE(sourceF3hat);
    DELETE(sourceQhat);
    DELETE(ZNAindex);
  }

};





class FWH {
  
private:
  
  
  vector<Face> faceVec;
  double gamma,p_ref,rho_ref,c_ref,c2_ref;
 
  // ----------------------------
  // Zonal Noise Analysis (ZNA)
  complex<double> ** p_hat;
  complex<double> ** p_hat_global;

  // parameters for wind-tunnel formulation, i.e. with coflow
  double U0[3];           // wind-tunnel velocity
  double M0;              // wind-tunnel Mach number
  double beta2;           // beta2 = 1-M0*M0

  double p_char;
  double f_char;
  double df_char;

  double oaspl_f_min;
  double oaspl_f_max;

  double wtime;
  double wtime_start;

  // to turn-off thickness and/or loading noise contributions
  double QF_switch[4];

  // debug switch
  int checkRetardedTime;
    
public:
  
  int num_freq_total;     // number of frequency accessible by postprocessing
  int num_freq;           // number of frequency kept for analysis
  Observer * ObserverPtr; // List of observer (i.e., microphone) arrays
  double * window;        // Hanning or None
  double window_coeff;    // for correction of amplitude and energy if window applied
  complex<double> * iomega;       // i times angular frequencies 
  int nt;                 // number of frames used for the FFT
  int nt_total;           // total number of frames (i.e., snapshots)
  double dt;              // time step in LES calculation
  int start,end,delta;    //
  bool dB_calc;

  // ----------------------------
  // FFT overlap for Welch method
  double fft_overlap;    // from input file, typically 0.75 (i.e. 75%)
  double fft_bandwidth;  // from input file, typically St = 0.05
  int nblock;            // total number of blocks with requested bandwidth and overlap
  int dimension;
  bool thetaX_hom;

  // solid vs permeable formulation
  bool solidFWH;

  // to turn on/off convective formulation
  bool convectiveFWH;

  // to apply rotation of the coordinate system if coflow not aligned in x-direction
  bool coord_rot;
  double RyRz[3][3];

  // ----------------------------
  // Zonal Noise Analysis
  int nZNA;                    // number of zone(s)
  int ZNAcoord;                // 0 for x, 1 for y, 2 for z, common to all zones
  double (*ZNAlimit)[2];       // min & max limit for each zone
  
  // ---------------------------- 
  // enable z-aligned arc
  int axis;
 
  FWH() {

    // for timing...
    wtime = MPI_Wtime();
    wtime_start = wtime;

    // initialize    
    if (mpi_rank == 0){
      cout << " ***  FWH::Initialize " << endl;
    }
    
    ZNAlimit = NULL;
    p_hat = NULL;
    p_hat_global = NULL;
    ObserverPtr = NULL;
    window = NULL;
    iomega = NULL;

    // read input file
    initParameters();

    // intialize coordinate rotation (if needed)
    initRotation();   

    // ----------------------------
    // read the surface files
    readData();
    
    // ----------------------------
    // Zonal Noise Analysis (ZNA)
    // initialize FW-H surface breakdown
    initZNAindex();

    // DEBUG *********************************************
    // to look at the time history of the variable on the FW-H surface every checkdata panels
    int checkData = getIntParam("DEBUG_TIME_HISTORY",0);
    if (checkData > 0) DebugCheckDataTimeHistory(checkData);

    // timing
    if (mpi_rank == 0) {
      double wtime0 = wtime;
      wtime = MPI_Wtime();
      cout <<"-------------- Approximate runtime for initialization = " << wtime - wtime0  << " s" << endl;
    }
  }
  


  ~FWH() {
    
    // destroy sources
    if (mpi_rank == 0){
      cout << " ***  FWH::destroySources() ... " << endl;
    }

    // deallocate memory
    if (ObserverPtr != NULL) delete ObserverPtr;   
 
    DELETE(window);
    DELETE(iomega);

    // ----------------------------
    // Zonal Noise Analysis (ZNA)
    if (ZNAlimit != NULL) delete[] ZNAlimit;

    deleteComplexArray2d(p_hat);

    deleteComplexArray2d(p_hat_global);

  }


  
  
  
  
  void run() {
    double x_obs[3];
    
    if (mpi_rank == 0){
      cout << " ***  FWH::Main solver ... " << endl;
    }

    //********************************************************
    // coordinate rotation if coflow not in x-direction
    if (coord_rot){
      // rotate coflow vector
      if (mpi_rank == 0){
	cout << " > coflow vector before rotation: " << U0[0] << "  " << U0[1] << "  " << U0[2] << " ( mag = " << MAG(U0) << " )" << endl;
      }
      applyRotation(RyRz, U0[0], U0[1], U0[2]);
      if (mpi_rank == 0){
	cout << " > coflow vector AFTER rotation (should be x-aligned): " << U0[0] << "  " << U0[1] << "  " << U0[2] << endl;
      }

      // rotate face normal vector and position vector
      for (int ifa = 0; ifa < faceVec.size(); ++ifa) {	
	applyRotation(RyRz, faceVec[ifa].nx,  faceVec[ifa].ny, faceVec[ifa].nz);
	applyRotation(RyRz, faceVec[ifa].x,  faceVec[ifa].y, faceVec[ifa].z);
      }
	
      // if permeable formulation, rotate velocity vector
      if (!solidFWH) {
	for (int ifa = 0; ifa < faceVec.size(); ++ifa) {
	  for (int it = 0; it < nt; ++it) {
	    applyRotation(RyRz, faceVec[ifa].sourceF1[it], faceVec[ifa].sourceF2[it], faceVec[ifa].sourceF3[it]);
	  }
	}
      }
    }


    //********************************************************
    // solid vs permeable (default) formulation
    if (solidFWH){
      if (mpi_rank == 0){
	cout << " > solid FW-H formulation " << endl;
      }
      // No calculation of source term required, only P'
    }
    else{
      if (mpi_rank == 0){
	cout << " > permeable FW-H formulation " << endl;
      }
      // Calculate the source terms that are observer independent  
      // Check if Wind-Tunnel mode i.e., with coflow
      if (U0[0]==0.0){
	calc_sources_obs_indep();
      }
      else{
	calc_sources_obs_indep_WT();
      }
    }
        


    // initializing observer locations
    initObservers();

    // initialize source arrays
    initFFT();

    // loop over the FFT blocks
    for (int iblock = 0; iblock < nblock; ++iblock){

      int block_start=floor(iblock*nt*(1-fft_overlap));
      if (mpi_rank == 0){
	cout << " > Processing block of data # " << iblock +1 << " out of " << nblock << " (snapshots " << block_start + 1 << " to " << block_start + nt << ")" << endl;
      }

      // Remove mean, apply Hanning window and compute FFT of source terms
      fftSources(block_start);

      // Option to turn-off thickness and/or loading noise contributions
      updateSources();

      // Loop over observers
      ObserverPtr->set_this_array(0);
      ObserverPtr->set_this_obs(0);
      for (int i_obs = 0; i_obs < ObserverPtr->get_num_obs(); ++i_obs) {
	if ((mpi_rank == 0)&&(i_obs%10 == 0)){
	  cout << " > working on observer " << i_obs+1 << " of " << ObserverPtr->get_num_obs() << endl;
	}
	
	ObserverPtr->get_x_obs(x_obs);
	
	//********************************************************
	// initialize output to 0 for current observer
	complex<double> z0(0.0,0.0);
	for (int iZNA=0; iZNA<nZNA; ++iZNA){
	  for (int i_freq = 0; i_freq < num_freq; ++i_freq){
	    p_hat[iZNA][i_freq] = z0;
	  }
	}
	
	//********************************************************
	// coordinate rotation if coflow not in x-direction
	if (coord_rot) applyRotation(RyRz, x_obs);
	
	//********************************************************
	// solid (1) vs permeable (default - 0) formulation
	if (solidFWH){
	  // Check if Wind-Tunnel mode (i.e., with coflow) & if convective effect neglected
	  if ((U0[0]==0.0)||(!convectiveFWH)){
	    calc_sources_obs_dep_solidFWH(x_obs);
	  }
	  else{
	    calc_sources_obs_dep_WT_solidFWH(x_obs);
	  }
	}
	else{
	  // Check if Wind-Tunnel mode (i.e., with coflow) & if convective effect neglected
	  if ((U0[0]==0.0)||(!convectiveFWH)){
	    calc_sources_obs_dep(x_obs);
	  }
	  else{
	    calc_sources_obs_dep_WT(x_obs);
	  }
	}
	
	SurfaceIntegration();

	ObserverPtr->save_p_hat(p_hat_global,iblock);

	ObserverPtr->calcPSDandOASPL(iblock);

	ObserverPtr->inc_this_obs();
      }
    }
      
    // timing
    if (mpi_rank == 0) {
      double wtime0 = wtime;
      wtime = MPI_Wtime();
      cout <<"-------------- Approximate runtime for computation = " << wtime - wtime0  << " s" << endl;
    }
  }
  



  void finalize() { 
    
    if (mpi_rank == 0){
      cout << " ***  FWH::Finalize output ... " << endl;
    }
    
    ObserverPtr->set_this_array(0);
    ObserverPtr->set_this_obs(0);

    
    //----------------------------------------------------------------------------
    //----------------------------------------------------------------------------
    //
    // FOR COMPUTATION IN TIME DOMAIN WITH OPTIONS P_TIME
    //    Compute global first and last time of reception over array of observers (default)
    // 
    // Option can be turned off with keyword DEBUG_RETARDED_TIME = 0
    //----------------------------------------------------------------------------
     
    int it_first;
    int it_last;
    double my_tfirst;       // local minimum of reception time (for output in time domain) 
    double my_tlast;        // local maximum of reception time (for output in time domain)

    // Retarded time calculation
    for (int i=0; i<ObserverPtr->get_num_array(); i++){   
      
      if(ObserverPtr->get_array_PLOTtime()) {	
	
	int istart = ObserverPtr->get_array_start();
	ObserverPtr->set_this_obs(istart);
	  
	int iend = ObserverPtr->get_array_end();
	
	double x_obs[3];
	

	// initialize local minimum and local maximum of reception time (for output in time domain)
	my_tfirst = -1000000000.0;
	my_tlast = 1000000000.0;
	
	for (int i_obs=istart; i_obs<iend; i_obs++){
	  ObserverPtr->get_x_obs(x_obs);
	  
	  calc_retarded_time(x_obs,my_tfirst,my_tlast);
	  
	  ObserverPtr->inc_this_obs();
	}
	
	double tfirst_global;
	double tlast_global;
	MPI_Allreduce(&my_tfirst,&tfirst_global,1,MPI_DOUBLE,MPI_MAX,mpi_comm);
	MPI_Allreduce(&my_tlast,&tlast_global,1,MPI_DOUBLE,MPI_MIN,mpi_comm);
	it_first = ceil((tfirst_global/dt-start)/delta);
	it_last = floor((tlast_global/dt-start)/delta);

	// ****************************
	// Recover full time history with retarded time calculation
	if (it_first > it_last) {
	  it_first = 0;
	  it_last  = nt - 1;
	  if (mpi_rank == 0) {
	    cout << " > WARNING !!! OPTIONS P_TIME: not enought time signal to compute retarded time for Array " << ObserverPtr->get_array_name() << endl;
	    cout <<"     The pressure time history will be reported from time " << (start + max(0,it_first)* delta)*dt << " to " << (start + min(nt-1,it_last)* delta)*dt << endl;
	  }
	}
	else {
	  if (mpi_rank == 0) {
	    cout << " > OPTIONS P_TIME: retarded time calculation for Array " << ObserverPtr->get_array_name() << endl;
	    cout <<"       first (global) reception time = " << tfirst_global <<" (discrete time " << ceil((tfirst_global/dt))*dt<< ", time step " << it_first<<") "<< endl;
	    cout <<"        last (global) reception time = " << tlast_global  <<" (discrete time " << floor((tlast_global/dt))*dt << ", time step " << it_last<<") " << endl;
	  }
	  

	  if (checkRetardedTime == 1) {
	    if (mpi_rank == 0) {
	      cout << " > WARNING !!! OPTIONS P_TIME: DEBUG_RETARDED_TIME = 1 , trying to extend end of time history with retarded time calculation ... " << endl;
	    }
	  }
	  else if (checkRetardedTime == 2) {
	    it_last = it_first+nt-1;
	    if (mpi_rank == 0) {
	      cout << " > WARNING !!! OPTIONS P_TIME: DEBUG_RETARDED_TIME = 2, writing full time history (i.e., no retarded time calculation) " << endl;
	    }
	  }
	  else {
	    it_last = nt-1;
	    if (mpi_rank == 0) {
	      cout << " > Default behavior for OPTIONS P_TIME: truncating time history at time step " << it_last +1 << endl;
	    }
	  }
	  
	  if (mpi_rank == 0) {
	    cout <<"     The pressure time history will be reported from time " << (start + it_first* delta)*dt << " to " << (start + (it_last+1)* delta)*dt << endl;
	  }
	}

	// save global first and last time of reception of homogeneous array
	ObserverPtr->set_array_it_first(it_first);
	ObserverPtr->set_array_it_last(it_last);

      }

      // got to next array of observers
      ObserverPtr->inc_this_array();
      
    }
    //----------------------------------------------------------------------------
    //----------------------------------------------------------------------------
         

    // Save all the outputs
    ObserverPtr->finalize_output(dB_calc, window_coeff, start, delta, dt); 

    // timing
    if (mpi_rank == 0) {
      wtime = MPI_Wtime();
      cout <<"-------------- Approximate total runtime = " << wtime - wtime_start  << " s" << endl;
    }
    
  }


      
      
      
  // DEBUG *********************************************
  virtual void ExactSolution(const double *mpole_mic, char *name) {

    if (mpi_rank == 0) {
      cout << " > Computing exact solution ... " << endl;
    } 
  }

    
  // DEBUG *********************************************
  void DebugCheckDataTimeHistory(int checkData){
    if (mpi_rank == 0){
      char filename[128];
      if (solidFWH){
	cout << "DEBUG: extracting time history of P' on solid surface for proc 0 every " << checkData << " panels" << endl;
	// save only pressure
	for (int i = 0; i < faceVec.size(); ++i) {
	  if (i%checkData == 0){
	    sprintf(filename,"%s_%06d.dat","DataTimeHistory_solid",i);
	    FILE * fp=fopen(filename,"w");
	    assert(fp != NULL);
	    fprintf(fp,"VARIABLES =\"X\" \"Y\" \"Z\" \"NT\" \"P'\" \n");
	    
	    for (int j = 0; j < nt; ++j){
	      fprintf(fp,"%12.6f %12.6f %12.6f %8.2f %18.12e\n", faceVec[i].x, faceVec[i].y, faceVec[i].z, (double)j, faceVec[i].sourceQ[j]);
	    }
	    
	    fclose(fp);
	  }
	}
      }
      else{
	cout << "DEBUG: extracting time history of P',u,v,w on permeable surface for proc 0 every " << checkData << " panels" << endl;
	// save pressure and velocity
	for (int i = 0; i < faceVec.size(); ++i) {
	  if (i%checkData == 0){
	    sprintf(filename,"%s_%06d.dat","DataTimeHistory",i);
	    FILE * fp=fopen(filename,"w");
	    assert(fp != NULL);
	    fprintf(fp,"VARIABLES =\"X\" \"Y\" \"Z\" \"NT\" \"P'\" \"U\" \"V\" \"W\" \n");
	    
	    for (int j = 0; j < nt; ++j){
	      fprintf(fp,"%12.6f %12.6f %12.6f %8.2f %18.12e %18.12e %18.12e %18.12e\n", faceVec[i].x, faceVec[i].y, faceVec[i].z, (double)j, faceVec[i].sourceQ[j], faceVec[i].sourceF1[j], faceVec[i].sourceF2[j], faceVec[i].sourceF3[j]);
	    }
	    
	    fclose(fp);
	    
	  }
	}
      }
    } 
  }
  
  // DEBUG *********************************************
  void DebugAddRotation(int nfa, int nsteps_global){
     if (Param * param = getParam("FRAME_ROTATION")) {
       double frame_rotation[3];
       FOR_I3 frame_rotation[i] = 0.0;
       // user has specified frame rotation in rotations/unit time...
       //frame_rotation[0] = getCurrentParamDouble(0);
       frame_rotation[0] = param->getDouble(0);
       if (mpi_rank == 0)
	 cout << " > FRAME_ROTATION active with the following rev/time: " << frame_rotation[0] << " " << frame_rotation[1] << " " << frame_rotation[2] << endl;
       // convert to rads/time unit...
       frame_rotation[0] *= 2.0*M_PI;
     

       for (int ifa = 0; ifa < nfa; ++ifa) {
	 double OmegaxR[3];
	 OmegaxR[0] = frame_rotation[1]*faceVec[ifa].z - frame_rotation[2]*faceVec[ifa].y;
	 OmegaxR[1] = frame_rotation[2]*faceVec[ifa].x - frame_rotation[0]*faceVec[ifa].z;
	 OmegaxR[2] = frame_rotation[0]*faceVec[ifa].y - frame_rotation[1]*faceVec[ifa].x;
	 for (int istep = 0; istep < nsteps_global; ++istep) {
	   faceVec[ifa].sourceF1[istep]+= OmegaxR[0]; //u
	   faceVec[ifa].sourceF2[istep]+= OmegaxR[1]; //v
	   faceVec[ifa].sourceF3[istep]+= OmegaxR[2]; //w
	 }
       }
     }
   }



  



  
private:



  void initParameters(){
    if (mpi_rank == 0){
      cout << " ***  FWH::initParameters() " << endl;
    }


    //********************************************************
    //OPTIONAL: initialize 2D or 3D solver (Default: 3D)
    dimension  = getIntParam("DIMENSION",3);
    if (dimension == 3){
      if (mpi_rank == 0){
	cout << " > Setting up 3D Ffowcs Williams - Hawkings solver ... " << endl;
      }
    }
    else{
      if (mpi_rank == 0)
	cerr << "Error: unrecognized or unavailable option DIMENSION: " << dimension << endl;
      throw(0);
    }


    //********************************************************
    // solid vs permeable (default) formulation
    // OPTIONAL: can force solver into solid formulation with permeable data (velocity field will not used) 
    solidFWH = getBoolParam("SOLID_FWH",false);
    if (solidFWH){
      if (mpi_rank == 0){
	cout << " > Using solid formulation (regardless of type of input data on FWH surface) ... " << endl;
      }
    }


    //********************************************************
    // OPTIONAL: initialize azimuthal direction (i.e., thetaX) as homogeneous (Default: true)
    thetaX_hom = getBoolParam("THETAX_HOM",true);
    if (mpi_rank == 0){
      if (thetaX_hom){
	cout << " > Setting up azimuthal direction as homogeneous (i.e., thetaX averaging on noise results)" << endl;
      }
      else{
	cout << " > No homogeneous direction (i.e., no thetaX averaging on noise results)" << endl;
      }
    }


    //********************************************************
    // initialize the physical variables
    gamma   = getDoubleParam("GAMMA",1.4);
    c_ref   = getDoubleParam("C_REF",1.0);
    rho_ref = getDoubleParam("RHO_REF",1.0);
    dt      = getDoubleParam("DT");
    c2_ref  = c_ref*c_ref;
    p_ref   = rho_ref * c2_ref / gamma;
    
    
    //********************************************************
    // initialize wind-tunnel velocity, i.e., mean coflow assumed positive in +x direction
    FOR_I3 U0[i] = 0.0;
    Param * param = getParam("U_COFLOW");    
    if (param != NULL) {
      int i = 0;
      while ((i < param->size()) && (i < 3)) {
	U0[i] = param->getDouble(i);
	++i;
      }
    }
    M0 = sqrt(DOT_PRODUCT(U0,U0))/c_ref;
    beta2 = 1.0-M0*M0;


    // DEBUG *********************************************
    // to export pressure signal with the option P_TIME
    // 0: retarded time & truncation for start & end of time history respectively (default)
    // 1: retarded time for both start & end of time history (use with caution ...)
    // 2: straight inverse Fourier transform of the pressure signal (not recommended)
    checkRetardedTime = getIntParam("DEBUG_RETARDED_TIME",0); 


    //********************************************************
    // Ability to overwrite p_ref for Helmholtz solver and single-precision database 
    // Set p_ref to 0 in these cases
    param = getParam("P_REF");
    if (param != NULL) {
      p_ref = param->getDouble();
	if (mpi_rank == 0){
	  cout << " > Setting reference pressure to p_ref = " << p_ref << endl;
	}
    }
    else{
      if (mpi_rank == 0){
	cout << " > Using default reference pressure p_ref = " << p_ref << endl;
      }
    }


    //********************************************************
    // ability to turn-off thickness and/or loading noise contributions
    // Default: all on (1)
    FOR_I4 QF_switch[i] = 1.0;
    param = getParam("SOURCE_SWITCH");
    if (param != NULL) {
      FOR_I4 QF_switch[i] = param->getDouble(i);
      if (mpi_rank == 0){
	cout << " > Setting thickness noise calculation = " << QF_switch[0] << endl;
	cout << " > Setting loading noise (x-dir) calculation = " << QF_switch[1] << endl;
	cout << " > Setting loading noise (y-dir) calculation = " << QF_switch[2] << endl;
	cout << " > Setting loading noise (z-dir) calculation = " << QF_switch[3] << endl;
      }
    }

    //********************************************************
    // enable z-axis arc
    param = getParam("AXIS");
    if (param != NULL) {
      axis = param->getInt();
    }
    else {
      axis =1;
    }
    if (mpi_rank == 0){
      cout << " > Observer arc with respect to axis " << axis << " (1 for x, 3 for z, y not implemented yet)" <<endl;
    }
    

    //********************************************************
    // convective formulation
    // OPTIONAL: can force solver to neglect convective effects during propagation even with coflow specified
    convectiveFWH = getBoolParam("CONVECTIVE_EFFECT",TRUE);
    if (!convectiveFWH){
      if (mpi_rank == 0){
	cout << " > Neglecting convective effect during propagation ... " << endl;
      }
    }


    //********************************************************
    // OPTIONAL: initialize units for spectra
    param = getParam("SPECTRA_UNIT");
    if (param != NULL) {
      string name = param->getString();
      if ((name == "dB")||(name == "DB")) {
	dB_calc=true;
	if (mpi_rank == 0){
	  cout << " > Compute spectra in dB rel. 20 10-6 Pa " << endl;
	}
      }
      else if (name == "NONE") {
	dB_calc=false;
	if (mpi_rank == 0){
	  cout << " > Compute spectra in default units" << endl;
	}
      }
      else {      
	if (mpi_rank == 0)
	  cerr << "Error: unrecognized SPECTRA_UNIT: " << name << endl;
	throw(0);
      }
    }
    else {
      dB_calc=false;
      if (mpi_rank == 0)
	cout << " > SPECTRA_UNIT: NONE (default)" << endl;
    }


    //********************************************************
    // OPTIONAL: conversion to dimensional units (if setup in input file)
    // p_char = rho_char u^2_char: characteristic pressure used for non-dimensionalization of p
    // f_char = u_char/l_char: characteristic frequency used for non-dimensionalization of f
    // deltaf_char: characteristic frequency used for non-dimensionalization of frequency bandwidth df (i.e., St)
    //double c_char = getDoubleParam("C_CHAR",1.0); // Removed to clarify non-dimensionalization
    double rho_char = getDoubleParam("RHO_CHAR",1.0);
    double l_char = getDoubleParam("L_CHAR",1.0);
    double u_char = getDoubleParam("U_CHAR",1.0);   
    p_char = rho_char*u_char*u_char;
    f_char = u_char/l_char;
    df_char = getDoubleParam("DF_CHAR",f_char);
    if (mpi_rank == 0){
      cout << " > (Optional) Conversion to characteristic units for the output, as defined in input file: " << endl;
      cout << "      pressure will be multiplied by RHO_CHAR (U_CHAR)^2 = " << p_char << endl;
      cout << "      frequency will be multiplied by U_CHAR / L_CHAR = " << f_char << endl;
      cout << "      PSD bandwidth will be multiplied by DF_CHAR = " << df_char << endl;
    }
    

    //********************************************************
    // OPTIONAL: initialize frequency range for OASPL (Default: all)
    param = getParam("OASPL_RANGE");
    if (param != NULL) {
      int i = 0;
      oaspl_f_min = param->getDouble(i++);
      oaspl_f_max = param->getDouble(i++);
      if (mpi_rank == 0){
	cout << " > Compute OASPL for frequencies from " << oaspl_f_min << " to "<< oaspl_f_max << endl;
      }
    }
    else{
      oaspl_f_min = 0.0;
      oaspl_f_max = 100000000.0; // just a very large number
      if (mpi_rank == 0){
	cout << " > Compute OASPL for all available frequencies (default)" <<  endl;
      }
    }


    //********************************************************
    // OPTIONAL: setup Zonal Noise Analysis (Default: none)
    param = getParam("ZONAL_NOISE_ANALYSIS");
    if (param != NULL) {
      int i = 0;
      string token = param->getString(i++);
      if      (token == "NZONE_X"){
	ZNAcoord = 0;
      }
      else if (token == "NZONE_Y"){
	ZNAcoord = 1;
      }
      else if (token == "NZONE_Z"){
	ZNAcoord = 2;
      }
      else{
	if (mpi_rank == 0)
	  cerr << "\n\n\n*********************************************************\n" <<
	    "Error: unrecognized keyword in ZONAL_NOISE_ANALYSIS: " << token <<
	    ". Example of proper syntax:\n" << 
	    "ZONAL_NOISE_ANALYSIS NZONE_X = 2  RANGE = 0.0 2.0 2.0 10.0 \n" << 
	    "*********************************************************\n\n\n" << endl;
	throw(0);	
      }
      nZNA = param->getInt(i++);
      nZNA++; // first zone (index 0) is alway for the total surface
      assert(nZNA > 1);
      
      // allocate memory - limits for zone 0 not used and set to arbitruary large numbers
      // might read limits from surface file in the future ...
      ZNAlimit=new double[nZNA][2];
      ZNAlimit[0][0]=-1e6;
      ZNAlimit[0][1]= 1e6;

      token = param->getString(i++);
      if (token != "RANGE"){
	if (mpi_rank == 0)
	  cerr << "\n\n\n*********************************************************\n" <<
	    "Error: unrecognized keyword in ZONAL_NOISE_ANALYSIS: " << token <<
	    ". Example of proper syntax:\n" << 
	    "ZONAL_NOISE_ANALYSIS NZONE_X = 2  RANGE = 0.0 2.0 2.0 10.0 \n" << 
	    "*********************************************************\n\n\n" << endl;
	throw(0);	
      }
      for (int iZNA=1; iZNA< nZNA; ++iZNA){
	ZNAlimit[iZNA][0] = param->getDouble(i++);
	ZNAlimit[iZNA][1] = param->getDouble(i++);
	assert (ZNAlimit[iZNA][0] < ZNAlimit[iZNA][1]);
	if (mpi_rank == 0){
	  cout << " > Zonal Noise Analysis for zone " << iZNA << ": direction " << ZNAcoord << " (0 for x, 1 for y, 2 for z); limits: min = "<< ZNAlimit[iZNA][0] << ",  max = " << ZNAlimit[iZNA][1] << endl;
	}
      }
    }
    else{
      nZNA = 1;
      ZNAcoord = 0;
    }


    //********************************************************
    // initialize the number and range of time snapshots 
    if (Param* param = getParam("DATA_RANGE")) {
      start = param->getInt(0);
      delta = param->getInt(1);
      end = param->getInt(2);
    }
    else {
      CERR(" > cannot find DATA_RANGE.");
    }
    //setCurrentParam("DATA_RANGE");
    //start = getCurrentParamInt(0);   
    //delta = getCurrentParamInt(1);
    //end = getCurrentParamInt(2);
    
    if ( (end-start)%delta != 0 ) {     
      end = (int)((end-start)/delta)*delta+start; 
      if (mpi_rank == 0)
	cout << " !!! WARNING: (end-start)%delta failed ... Check DATA_RANGE !!! Frames " << start << " to " << end << " every " << delta << " frame(s) will be used instead !!!" << endl;
    }
    
    nt_total = (end - start)/delta + 1;

    if (nt_total == 1) {
      CERR(" > must have at least 2 data files.");
    }

    if (mpi_rank == 0)
      cout << " > total number of time snapshots: " << nt_total << endl;


    //********************************************************
    // OPTIONAL: setup FFT overlap (Default: none)
    // syntax for 75% overlap: FFT_OVERLAP=0.75 BANDWIDTH=0.05
    param = getParam("FFT_OVERLAP");
    if (param != NULL) {
      int i = 0;
      fft_overlap = param->getDouble(i++);
      string token = param->getString(i++);
      if (token == "BANDWIDTH"){
	fft_bandwidth = param->getDouble(i++);
      }
      else{
	if (mpi_rank == 0)
	  cerr << "\n\n\n*********************************************************\n" <<
	    "Error: unrecognized keyword in FFT_OVERLAP: " << token <<
	    ". Example of proper syntax:\n" << 
	    "FFT_OVERLAP 0.75 BANDWIDTH=0.05\n" << 
	    "*********************************************************\n\n\n" << endl;
	throw(0);	
      }
      // sampling frequency
      double fs = 1.0/(double(delta)*dt);
      // number of snapshots in 1 block of requested bandwidth
      nt = floor(fs/fft_bandwidth);
      // total number of blocks with requested bandwidth and overlap 
      nblock = floor((nt_total - nt)/(nt*(1-fft_overlap))) + 1;
      if (mpi_rank == 0)
      cout << " > (Optional) FFT overlap: total number of fft blocks = " << nblock << ", with bandwidth " << fft_bandwidth << ", and overlap " << fft_overlap << " (fs = " << fs << " , nt = " << nt << ")" << endl;
    }
    else {
      fft_overlap = 1.0;
      fft_bandwidth = 1.0;
      nt = nt_total;
      nblock = 1;
    }


    //********************************************************
    // initialize the window for spectral analysis
    window  = new double[nt];
    param = getParam("WINDOW_TYPE");
    if (param != NULL) {
      string name = param->getString();
      if ((name == "HANNING")||(name == "HANN")) {
	window_coeff = sqrt(8.0/3.0);  // energy preserving with hanning window
	if (mpi_rank == 0){
	  cout << " > periodic Hanning window is initialized (with coefficient "<< window_coeff <<" for energy preservation, and coefficient 2.0 for amplitude preservation)  " << endl;
	}
	for (int it=0; it<nt; ++it){
	  window[it] = window_coeff * 0.5 * ( 1.0 - cos(2.0*M_PI*(double(it)/double(nt))));
	}
      }
      else if (name == "LOCKARD") {
	// ******** 
	// see paper JSV (2000) 229(4), pp 897-911 
	int nfilter=8;
	window_coeff = 0.0;
	int ifilter=floor(nt/nfilter);
	for (int it=0; it<nt; it++){
	  if ((it < ifilter) || (it > nt - ifilter)){
	    window[it] = 0.5 * ( 1.0 - cos(nfilter*M_PI*(double(it)/double(nt))));
	  }
	  else{
	    window[it] = 1.0;
	  }
	  window_coeff +=  window[it] *  window[it];
	}
	window_coeff = 1/sqrt(window_coeff/nt);
	for (int it=0; it<nt; it++){
	  window[it] *= window_coeff;
	}
	if (mpi_rank == 0){
	  cout << " > Hanning + top hat window is initialized, with coefficient "<< window_coeff <<" for energy preservation (see paper JSV (2000) 229(4), pp 897-911)  " << endl;
	}
      }
      else if (name == "NONE") {
	if (mpi_rank == 0){
	  cout << " > NO WINDOW FOR TIME SIGNAL" << endl;
	}
	for (int it=0; it<nt; ++it){
	  window[it] = 1.0;
	}
        window_coeff = 2.0;  // without hanning window
      }
      else {      
	if (mpi_rank == 0)
	  cerr << "Error: unrecognized WINDOW_TYPE: " << name << endl;
	throw(0);
      }
    }
    else {
      if (mpi_rank == 0)
	cerr << "Error: unspecified WINDOW_TYPE = ... The recommended options is HANNING " << endl;
      throw(0);
    }


    //********************************************************
    // OPTIONAL: initialize frequency range of interest for spectra (sampling period is delta*dt)
    num_freq_total = nt/2 + 1; // total number of frequencies accessible by postprocessing
    double maxf = getDoubleParam("MAX_FREQUENCY",-1.0);
    if (maxf != -1.0) {
      num_freq = int(ceil(maxf * (double(nt) * double(delta)*dt))) + 1;
      if (mpi_rank == 0)
	cout << " > (Optional) Reduced number of frequencies: " << num_freq << " (max frequency =  " << double(num_freq-1)/(double(nt) * double(delta)*dt) << ")" <<endl;

    }
    else {
      // Default: keep all the frequencies
      num_freq = num_freq_total;
      if (mpi_rank == 0)
	cout << " > (Default) Total number of frequencies: " << num_freq << " (max frequency =  " << double(num_freq-1)/(double(nt) * double(delta)*dt) << ")" <<endl;
    }


    //********************************************************
    // initialize frequencies, where the sampling period is delta*dt
    //omega = new double[num_freq];
    iomega = new complex<double>[num_freq];
    complex<double> zi(0.0,1.0);
    for (int i_freq=0; i_freq  < num_freq; ++i_freq){
      iomega[i_freq] = -zi*(2.0 * M_PI * double(i_freq)/(double(nt) * double(delta)*dt));
    }


    //********************************************************
    // initialize p_hat
    p_hat = newComplexArray2d(nZNA,num_freq);
    if (mpi_rank == 0){
      p_hat_global = newComplexArray2d(nZNA,num_freq);
    } 
    else 
      p_hat_global = newComplexArray2d(1,1); 
  }
  


  // ############################################################
  // support both solid and permeable formulation
  void readData() {

    if (mpi_rank == 0)
      cout << " ***  readData() " << endl;

    //********************************************************
    // parse what we need for the files...
    
    string dataprefix = getStringParam("DATA_PREFIX");

    //********************************************************
    // OPTIONAL: initialize the file extension for the FW-H data (default: .data)
    Param * param = getParam("DATA_SUFFIX");
    string datasuffix;
    if (param != NULL) {
      datasuffix = param->getString();
    }
    else{
      datasuffix = "data";
    }
     
    if (mpi_rank == 0)
      cout << " > total number of surface files: " << nt_total << endl;

    //********************************************************
    // read the geom file...
    
    char filename[128];
    sprintf(filename,"%s.geom",dataprefix.c_str());


    if (mpi_rank == 0)
      cout << " > reading file " << filename << endl;

    MPI_File fh;
    if (MPI_File_open(mpi_comm,filename,
		      MPI_MODE_RDONLY,
		      MPI_INFO_NULL,&fh) != 0) {
      if (mpi_rank == 0)
	cerr << "Error: could not open file: " << filename << endl;
      throw(0);
    }

    // ------------------------------------------
    // HEADER
    // 0: a magic number (for future byte swapping)
    // 1: version
    // 2: the data layout flag (e.g. one I1, D1, etc...)
    // 3: count
    //
    // data layout:
    //CTIDefs.hpp:#define FFW_DATA_I1                    1
    //CTIDefs.hpp:#define FFW_DATA_D1                    2
    //CTIDefs.hpp:#define FFW_DATA_D2                    3
    //CTIDefs.hpp:#define FFW_DATA_D1D2                  4
    //CTIDefs.hpp:#define FFW_DATA_D2D2                  5
    // ------------------------------------------    
    
    // assume we do not have to byte_swap...
    bool byte_swap = false;

    int itmp[4];
    MPI_Status status;
    if (mpi_rank == 0) { 
      MPI_File_read(fh,itmp,4,MPI_INT,&status);
    }
    MPI_Bcast(itmp,4,MPI_INT,0,mpi_comm); 

    // check byte swap
    if (itmp[0] != UGP_IO_MAGIC_NUMBER) {
      ByteSwap::byteSwap(itmp, 4);
      if (itmp[0] == UGP_IO_MAGIC_NUMBER) {
	byte_swap = true;
	if (mpi_rank == 0)
	  cout << "File requires byte swapping." << endl;
      }
      else{
	if (mpi_rank == 0) 
	  cerr << "Error: file "<< filename << " does not start as expected. Aborting ..."<< endl;
	throw(0);
      }
    }


    // check version
    if (itmp[1] != UGP_IO_VERSION) {
      if (mpi_rank == 0)
	cout << "WARNING: data io version (" << itmp[1] << ") differs from current version (" << UGP_IO_VERSION << "). Should still be compatible ..." << endl;
      /*
	cerr << "Error: io version differs from current version: " << itmp[1] << endl;
      throw(0);
      */
    }


    // check data layout
    int dataLayout = itmp[2];
    if (dataLayout != FFW_DATA_D2D2){
      if (mpi_rank == 0)
	cerr << "Error: file "<< filename << " has wrong type "<< dataLayout << " (Data layout "<< FFW_DATA_D2D2 << " expected). Aborting ..."<< endl;
      throw(0);
    }
    
    
    // check count
    int nfa_global = itmp[3];
    if (mpi_rank == 0)
      cout << " > Total number of faces (i.e. panels): nfa_global = " << nfa_global << endl;




    // uniformly distribute these faces over the available processes...  
    int * faora = NULL; buildUniformXora(faora,nfa_global); 

    // we can get our local face size from this now...
    int nfa = faora[mpi_rank+1]-faora[mpi_rank];
    int nfa3 = 3*nfa;


    // resize the face vector...
    faceVec.resize(nfa);

    
    // Need d2_type (position, normals and velocities) and d1_type (pressure)
    int my_disp3 = 3*faora[mpi_rank];
   
    // and read the positions and normals from the geom file...
    double * dbuf = new double[nfa3];
    MPI_Offset offset = 4*4; // i.e. 4 int
    
    MPI_File_read_at_all(fh, offset+my_disp3*8, dbuf, nfa3, MPI_DOUBLE, &status);
    if(byte_swap) ByteSwap::byteSwap(dbuf, nfa3);
    for (int i = 0; i < nfa; ++i){
      faceVec[i].x = dbuf[3*i  ];
      faceVec[i].y = dbuf[3*i+1];
      faceVec[i].z = dbuf[3*i+2];
    }
    
    offset += 3*8*nfa_global; // i.e. 3*nfa_global double
    
    MPI_File_read_at_all(fh, offset+my_disp3*8, dbuf, nfa3, MPI_DOUBLE, &status);
    if(byte_swap) ByteSwap::byteSwap(dbuf, nfa3);
    for (int i = 0; i < nfa; ++i) {
      faceVec[i].nx = dbuf[3*i  ];
      faceVec[i].ny = dbuf[3*i+1];
      faceVec[i].nz = dbuf[3*i+2];
    }
    
    MPI_File_close(&fh);
    DELETE(dbuf);

    // look at the global normal sum and area...
    double my_buf[4] = { 0.0, 0.0, 0.0, 0.0};
    for (int i = 0; i < nfa; ++i) {
      my_buf[0] += faceVec[i].nx;
      my_buf[1] += faceVec[i].ny;
      my_buf[2] += faceVec[i].nz;
      my_buf[3] += sqrt(faceVec[i].nx*faceVec[i].nx+faceVec[i].ny*faceVec[i].ny+faceVec[i].nz*faceVec[i].nz);
    }
    double buf[4];
    MPI_Reduce(my_buf,buf,4,MPI_DOUBLE,MPI_SUM,0,mpi_comm);
    if (mpi_rank == 0){
      cout << " > surface normal sum: " << buf[0] << " " << buf[1] << " " << buf[2] << ", total area: " << buf[3] << endl;
    }
    
    //********************************************************
    // read the data file...
    // ----------------------------
    // memory management
    int my_ierr = 0;

    // open the serial file... 
    FILE * fhs = NULL; 

    try {
           
      // ###########################################################
      // say we have nsteps_global time steps of data with nfa_global faces in each set. We start
      // by reading in the data in a partitioned way so that each processor has 
      // approximately nsteps_global/np complete sets of face data...
      
      int nsteps_global = nt_total;
      int * stora = NULL; buildUniformXora(stora,nsteps_global); // step-of-rank table
      int nsteps = stora[mpi_rank+1]-stora[mpi_rank]; // we can get the steps we own like this
      
      if (mpi_rank == 0){
	cout << "global size: " << nsteps_global << " steps x " << nfa_global << " faces." << endl;
	cout << " local size: " << nsteps << " steps x " << nfa_global << " faces." << endl;
      }

      int my_solidFWH = 0;
      
      // initially ALL faces for a subset of the steps are read on one processor
      // pressure p is always needed
      double ** p = newArray2d(nsteps,nfa_global);
      // velocity u, v, w might not be needed if solid surface formulation      
      double ** u = NULL;
      double ** v = NULL;
      double ** w = NULL;

#ifdef SERIAL_IO
      FOR_RANK {
        const int nsteps = stora[rank+1]-stora[rank]; // we can get the steps we own like this
        for (int step = 0; step < nsteps; ++step) {
          const int step_global = (step + stora[rank])*delta + start;

          if (mpi_rank == 0) {
            buildIndexedFilename(filename,dataprefix.c_str(),step_global,datasuffix.c_str());
            cout << " > reading data from surface files: step " << step+1 << " out of " << nsteps <<" (filename = " <<  filename <<")"<< endl;

            if ( (fhs=fopen(filename,"r"))==NULL ) {
              cerr << "Error: cannot open file " << filename << endl;
              throw(-1);
            }
          }
	
          // read in data...
          // ------------------------------------------
          // HEADER
          // 0: a magic number (for future byte swapping)
          // 1: version
          // 2: the data layout flag (e.g. one I1, D1, etc...)
          // 3: count
          // ------------------------------------------  
          
          if (mpi_rank == 0) {
            fread(itmp, sizeof(int), 4, fhs);
            if (rank != 0) {
              MPI_Send(itmp,4,MPI_INT,rank,rank,mpi_comm);
            }
          }

          if (mpi_rank == rank) {
            if (rank != 0) {
              MPI_Recv(itmp,4,MPI_INT,0,rank,mpi_comm,MPI_STATUS_IGNORE);
            }

            // check byte swap
            if (itmp[0] != UGP_IO_MAGIC_NUMBER) {
              ByteSwap::byteSwap(itmp, 4);
              if (itmp[0] == UGP_IO_MAGIC_NUMBER) {
                byte_swap = true;
                cout << "File requires byte swapping." << endl;
              }
              else{
                cerr << "Error: file "<< filename << " does not start as expected. Aborting ..." << endl;
                throw(-1);
              }
            }

            // check version
            if (itmp[1] != UGP_IO_VERSION) {
              cout << "WARNING: data io version (" << itmp[1] << ") differs from current version (" << UGP_IO_VERSION << "). Should still be compatible ..." << endl;
              // cerr << "Error: io version differs from current version: " << itmp[1] << endl;
              // throw(-1);
            }

            // check data layout
            dataLayout = itmp[2];
            if ((dataLayout != FFW_DATA_D1) && (dataLayout != FFW_DATA_D1D2)){
              cerr << "Error: file "<< filename << " has wrong type "<< dataLayout << " (Data layout "<< FFW_DATA_D1D2 << " or " << FFW_DATA_D1 << " expected). Aborting ..." << endl;
              throw(-1);
            }

            // check count
            int nfa_global_check = itmp[3];
            if (nfa_global_check != nfa_global) {
              cerr << "Error: data count in file " << filename << " is not as expected: size " << nfa_global_check << " should be " << nfa_global << endl;
              throw(-1);
            }
          }

          if (dataLayout == FFW_DATA_D1) {
            // read pressure only: solid formulation 
            my_solidFWH = 1;
            int input_size = nfa_global;
            double * FWH_data = NULL;
            if (mpi_rank == 0) {
              FWH_data = new double[input_size];
              fread(FWH_data,  sizeof(double), input_size, fhs);
              if (rank != 0) {
                MPI_Send(FWH_data,input_size,MPI_DOUBLE,rank,rank,mpi_comm);
              }
            }
              
            if (mpi_rank == rank) {
              if (rank != 0) {
                FWH_data = new double[input_size];
                MPI_Recv(FWH_data,input_size,MPI_DOUBLE,0,rank,mpi_comm,MPI_STATUS_IGNORE);
              }
              if (byte_swap) ByteSwap::byteSwap(FWH_data, input_size);

              for (int ifa = 0; ifa < nfa_global; ++ifa){
                p[step][ifa] = FWH_data[ifa];
              }
            }
            DELETE(FWH_data);
          }
          // read pressure and velocity, if permeable formulation
          else {
            int input_size = 4*nfa_global;
            double *FWH_data = NULL;
            if (mpi_rank == 0) {
              FWH_data = new double[input_size];
              fread(FWH_data,  sizeof(double), input_size, fhs);
              if (rank != 0) {
                MPI_Send(FWH_data,input_size,MPI_DOUBLE,rank,rank,mpi_comm);
              }
            }

            if (mpi_rank == rank) {
              if (rank != 0) {
                FWH_data = new double[input_size];
                MPI_Recv(FWH_data,input_size,MPI_DOUBLE,0,rank,mpi_comm,MPI_STATUS_IGNORE);
              }
              if (byte_swap) ByteSwap::byteSwap(FWH_data, input_size);

              if (u == NULL) u = newArray2d(nsteps,nfa_global);
              if (v == NULL) v = newArray2d(nsteps,nfa_global);
              if (w == NULL) w = newArray2d(nsteps,nfa_global); 
              for (int ifa = 0; ifa < nfa_global; ++ifa){
                p[step][ifa] = FWH_data[ifa                 ];
                u[step][ifa] = FWH_data[nfa_global + 3*ifa  ];
                v[step][ifa] = FWH_data[nfa_global + 3*ifa+1];
                w[step][ifa] = FWH_data[nfa_global + 3*ifa+2];
              }
            }
            DELETE(FWH_data);
          }

          if (mpi_rank == 0) {
            assert(fhs);
            fclose(fhs);
          }
        }
      }
#else
      // each processor read its data independently. 
      // Cycle through the local steps locally owned ...
      int step_offset = stora[mpi_rank];
      for (int step = 0; step < nsteps; ++step) {
	const int step_global = (step + step_offset)*delta + start;
	buildIndexedFilename(filename,dataprefix.c_str(),step_global,datasuffix.c_str());

	if (mpi_rank == 0)
	  cout << " > reading data from surface files: step " << step+1 << " out of " << nsteps <<" (filename = " <<  filename <<")"<< endl;

	if ( (fhs=fopen(filename,"r"))==NULL ) {
	  cerr << "Error: cannot open file " << filename << endl;
	  throw(0);
	}
	
	// read in data...
	// ------------------------------------------
	// HEADER
	// 0: a magic number (for future byte swapping)
	// 1: version
	// 2: the data layout flag (e.g. one I1, D1, etc...)
	// 3: count
	// ------------------------------------------  
	
	fread(itmp, sizeof(int), 4, fhs);


	// check byte swap
	if (itmp[0] != UGP_IO_MAGIC_NUMBER) {
	  ByteSwap::byteSwap(itmp, 4);
	  if (itmp[0] == UGP_IO_MAGIC_NUMBER) {
	    byte_swap = true;
	    if (mpi_rank == 0)
	      cout << "File requires byte swapping." << endl;
	  }
	  else{
	    CERR("Error: file "<< filename << " does not start as expected. Aborting ...");
	    throw(0);
	  }
	}


	// check version
	if (itmp[1] != UGP_IO_VERSION) {
	  if (mpi_rank == 0)
	    cout << "WARNING: data io version (" << itmp[1] << ") differs from current version (" << UGP_IO_VERSION << "). Should still be compatible ..." << endl;
	  /*
	  CERR("Error: io version differs from current version: " << itmp[1]);
	  throw(0);
	  */
	}
	

	// check data layout
	dataLayout = itmp[2];
	if ((dataLayout != FFW_DATA_D1) && (dataLayout != FFW_DATA_D1D2)){
	  CERR("Error: file "<< filename << " has wrong type "<< dataLayout << " (Data layout "<< FFW_DATA_D1D2 << " or " << FFW_DATA_D1 << " expected). Aborting ...");
	  throw(0);
	}

	// check count
	int nfa_global_check = itmp[3];
	if (nfa_global_check != nfa_global) {
	  CERR("Error: data count in file " << filename << " is not as expected: size " << nfa_global_check << " should be " << nfa_global);
	  throw(0);
	}


	if (dataLayout == FFW_DATA_D1){
	  // read pressure only: solid formulation 
	  my_solidFWH = 1;
	  int input_size = nfa_global;
          double * FWH_data = new double[input_size];
	  fread(FWH_data,  sizeof(double), input_size, fhs);
	  if(byte_swap) ByteSwap::byteSwap(FWH_data, input_size);

	  for (int ifa = 0; ifa < nfa_global; ++ifa){
	    p[step][ifa] = FWH_data[ifa];
	  }
          delete[] FWH_data;
	}
	// read pressure and velocity, if permeable formulation
	else{
	  // allocate arrays
	  if (u == NULL) u = newArray2d(nsteps,nfa_global);
	  if (v == NULL) v = newArray2d(nsteps,nfa_global);
	  if (w == NULL) w = newArray2d(nsteps,nfa_global); 
	  int input_size = 4*nfa_global;
	  double *FWH_data = new double[input_size];
	  fread(FWH_data,  sizeof(double), input_size, fhs);
	  if(byte_swap) ByteSwap::byteSwap(FWH_data, input_size);

	  for (int ifa = 0; ifa < nfa_global; ++ifa){
	    p[step][ifa] = FWH_data[ifa                 ];
	    u[step][ifa] = FWH_data[nfa_global + 3*ifa  ];
	    v[step][ifa] = FWH_data[nfa_global + 3*ifa+1];
	    w[step][ifa] = FWH_data[nfa_global + 3*ifa+2];
	  }
          delete[] FWH_data;
	}

	fclose(fhs);

      }
#endif

      // now we can call build a transpose operator to populate pt from p...
      // Note: the order here is VERY important. provide current i first (the currently distributed
      // index), then current j (the currently full index)...      
      MPI_Barrier(mpi_comm);

      // in case there are more processors than steps, need to set the flag solidFWH 
      int global_solidFWH;
      MPI_Allreduce(&my_solidFWH,&global_solidFWH,1,MPI_INT,MPI_MAX,mpi_comm);
      if (global_solidFWH==1) solidFWH=true;

      if (mpi_rank == 0)
	cout << " Building transpose operator ... "<< endl;

      TransposeOperator * t = new TransposeOperator(stora,faora);       
      
      delete[] stora;
      delete[] faora;

      // the transposed data has ALL steps for a subset of the faces...
      double ** pt = newArray2d(nfa,nsteps_global);
      double ** ut = NULL;
      double ** vt = NULL;
      double ** wt = NULL;
      
      if (mpi_rank == 0)
	cout << " Building time history of pressure ... "<< endl;
      t->apply(pt,p);
      deleteArray2d(p);

      if (!solidFWH){
	ut = newArray2d(nfa,nsteps_global); 
	if (mpi_rank == 0)
	  cout << " Building time history of u velocity ... "<< endl;	
	t->apply(ut,u);
      }
      deleteArray2d(u);

      if (!solidFWH){
	vt = newArray2d(nfa,nsteps_global); 
	if (mpi_rank == 0)
	  cout << " Building time history of v velocity ... "<< endl;
	t->apply(vt,v);
      }
      deleteArray2d(v);
      
      if (!solidFWH){
	wt = newArray2d(nfa,nsteps_global); 
	if (mpi_rank == 0)
	  cout << " Building time history of w velocity ... "<< endl;
	t->apply(wt,w);
      }
      deleteArray2d(w);
      

      // and free the transform
      delete t;

      // set face data size and allocate face data one variable at a time to limit memory usage ... 
      // pressure
      for (int ifa = 0; ifa < nfa; ++ifa) {
	faceVec[ifa].sourceQ = new double[nt_total];
	for (int istep = 0; istep < nsteps_global; ++istep) {
	  // to reduce memory usage, p is stored in faceVec[i].sourceQ
	  faceVec[ifa].sourceQ[istep]= pt[ifa][istep] - p_ref;
	}
      }
      deleteArray2d(pt);
      
      // velocities if permeable surface formulation
      if (!solidFWH){
	// allocate array and save u velocity
	for (int ifa = 0; ifa < nfa; ++ifa) {
	  faceVec[ifa].sourceF1 = new double[nt_total];
	  for (int istep = 0; istep < nsteps_global; ++istep) {
	    // to reduce memory usage, u is stored in faceVec[i].sourceF1
	    faceVec[ifa].sourceF1[istep]= ut[ifa][istep];
	  }
	}
	deleteArray2d(ut);
	
	// allocate array and save v velocity
	for (int ifa = 0; ifa < nfa; ++ifa) {
	  faceVec[ifa].sourceF2 = new double[nt_total];
	  for (int istep = 0; istep < nsteps_global; ++istep) {
	    // to reduce memory usage, v is stored in faceVec[i].sourceF2
	    faceVec[ifa].sourceF2[istep]= vt[ifa][istep];
	  }
	}
	deleteArray2d(vt);
	
	// allocate array and save w velocity
	for (int ifa = 0; ifa < nfa; ++ifa) {
	  faceVec[ifa].sourceF3 = new double[nt_total];
	  for (int istep = 0; istep < nsteps_global; ++istep) {
	    // to reduce memory usage, w is stored in faceVec[i].sourceF3
	    faceVec[ifa].sourceF3[istep]= wt[ifa][istep];
	  }
	}
	deleteArray2d(wt);
      }

      // DEBUG *********************************************
      // DebugAddRotation(nfa,nsteps_global );


      /*
      
      if (mpi_rank == 0)
	cout << "DONE reading ... nfa = " << nfa <<endl;
      
      for (int istep = 0; istep < nsteps_global; ++istep) {
	double test_min = 1000000000.0;
	double test_max = -1000000000.0;
	for (int ifa = 0; ifa < nfa; ++ifa){

	// faceVec[ifa].sourceF1[istep]=faceVec[ifa].sourceQ[istep]-(pt[ifa][istep] - p_ref);
	// test_min = min(test_min,faceVec[ifa].sourceF1[istep]);
	// test_max = max(test_max,faceVec[ifa].sourceF1[istep]);
	  
	  test_min = min(test_min,faceVec[ifa].sourceQ[istep]);
	  test_max = max(test_max,faceVec[ifa].sourceQ[istep]);
	}
	
	double test_min_global;
	double test_max_global;
	MPI_Allreduce(&test_min,&test_min_global,1,MPI_DOUBLE,MPI_MIN,mpi_comm);
	MPI_Allreduce(&test_max,&test_max_global,1,MPI_DOUBLE,MPI_MAX,mpi_comm);
	if (mpi_rank == 0) {
	  cout << " > step " << istep << ", test_min_global = " << test_min_global <<", test_max_global = " << test_max_global << endl;
	}
      }

      */
      

      
    }
    catch(...){
      cout <<"Memory requirement too large on proc " << mpi_rank << endl;
      my_ierr = 1;
    }
      
    int ierr;
    MPI_Allreduce(&my_ierr,&ierr,1,MPI_INT,MPI_MAX,mpi_comm);
    if (ierr != 0) {
      if (mpi_rank == 0)
	cerr<< ERRSTART<< "Error in readData(): more memory/processors needed "<< ERREND<< endl;
      throw(0);
    }
      
  }
  


  void initObservers(){
    if (mpi_rank == 0){
      cout << " ***  FWH::initObservers() " << endl;
    }

    // need omega rather than -i*omega
    double * omega = new double[num_freq];
    for (int i_freq=0; i_freq  < num_freq; ++i_freq){
      omega[i_freq] =  abs(iomega[i_freq]); 
    }

    // initialize list of observers
    ObserverPtr = new Observer(nZNA, ZNAcoord ,ZNAlimit, nt, num_freq, omega, oaspl_f_min, oaspl_f_max, p_char, f_char, df_char, dimension, thetaX_hom, nblock, axis);
    
    // read observer info from input file
    //FOR_PARAM("OBSERVER") ObserverPtr->writeObserver(getCurrentParam());
    FOR_PARAM_MATCHING("OBSERVER") {
      ObserverPtr->writeObserver(&(*param));
    }

    // finalize the locations
    ObserverPtr->finalize_locations();

    DELETE(omega);

  }




    

  // ----------------------------
  // Zonal Noise Analysis (ZNA)
  // initialize ZNA face index
  void initZNAindex() {
    // select all the faces in zone 0 (total surface)
    for (int ifa = 0; ifa < faceVec.size(); ++ifa) {
      faceVec[ifa].ZNAsize = 1;
      faceVec[ifa].ZNAindex = new int[nZNA];
      faceVec[ifa].ZNAindex[0] = 0;
    }

    // select faces in additional zone if ZNA requested
    if (nZNA > 1){
      double faceCoord;
      for (int ifa = 0; ifa < faceVec.size(); ++ifa) {

	// check if zone range in x, y or z direction
	if (ZNAcoord == 0) {
	  faceCoord = faceVec[ifa].x;
	}
	else if (ZNAcoord == 1) {
	  faceCoord = faceVec[ifa].y;
	}
	else {
	  faceCoord = faceVec[ifa].z;
	}

	for (int iZNA = 1; iZNA < nZNA; ++iZNA) {
	  if ((ZNAlimit[iZNA][0] <= faceCoord) &&(faceCoord <= ZNAlimit[iZNA][1] )) {
	    faceVec[ifa].ZNAindex[faceVec[ifa].ZNAsize] = iZNA;
	    faceVec[ifa].ZNAsize += 1;
	  }
	}
      }

      // Debug: compute zone area and write zone to file
      char filename[128];
      if ((mpi_rank==0)&&(cti_verbose))
	cout << " > Zonal Noise Analysis" << endl;
      for (int iZNA = 1; iZNA < nZNA; ++iZNA) {
	sprintf(filename,"FWHzone%06d.plt",iZNA);
	if ((mpi_rank==0)&&(cti_verbose))
	  cout << "write Tecplot ASCII: " << filename << endl;
      
	FILE * fp;
	if ( mpi_rank == 0 ) {
	  if ( (fp=fopen(filename,"w"))==NULL ) {
	    cerr << "Error: cannot open file " << filename << endl;
	    throw(-1);
	  }
	  fprintf(fp,"TITLE = \"flagged FWH faces\"\n");
	  fprintf(fp,"VARIABLES = \"X\"\n");
	  fprintf(fp,"\"Y\"\n");
	  fprintf(fp,"\"Z\"\n");
	}
	else {
	  int dummy;
	  MPI_Status status;
	  MPI_Recv(&dummy,1,MPI_INT,mpi_rank-1,1234,mpi_comm,&status);
	  if ( (fp=fopen(filename,"a"))==NULL ) {
	    cerr << "Error: cannot open file " << filename << endl;
	    throw(-1);
	  }
	}
	

	double my_buf[5] = { 0.0, 0.0, 0.0, 0.0, 0.0};
	double buf[5];
	int lines_written = 0;
	for (int ifa = 0; ifa < faceVec.size(); ++ifa) {
	  for (int izone=1; izone < faceVec[ifa].ZNAsize; ++izone){
	    if (faceVec[ifa].ZNAindex[izone] == iZNA) {
	      fprintf(fp,"%18.15le ",faceVec[ifa].x);
	      fprintf(fp,"%18.15le ",faceVec[ifa].y);
	      fprintf(fp,"%18.15le ",faceVec[ifa].z);
	      if (++lines_written > 5) {
		fprintf(fp,"\n");
		lines_written=0;
	      }
	      if (lines_written != 0)
		fprintf(fp,"\n");

	      // look at the global normal sum and area...
	      my_buf[0] += faceVec[ifa].nx;
	      my_buf[1] += faceVec[ifa].ny;
	      my_buf[2] += faceVec[ifa].nz;
	      my_buf[3] += sqrt(faceVec[ifa].nx*faceVec[ifa].nx+faceVec[ifa].ny*faceVec[ifa].ny+faceVec[ifa].nz*faceVec[ifa].nz);
	      my_buf[4] += 1.0;
	    }
	  }	  
	}
	
        fclose(fp);          
	
	if ( mpi_rank < mpi_size-1 ) {
	  int dummy = 1;
	  MPI_Send(&dummy,1,MPI_INT,mpi_rank+1,1234,mpi_comm);
	}

	MPI_Barrier(mpi_comm);

	MPI_Reduce(my_buf,buf,5,MPI_DOUBLE,MPI_SUM,0,mpi_comm);
	if (mpi_rank == 0){
	  cout << " > FW-H zone " << iZNA <<",number of faces (i.e. panels): "<< (int)buf[4] <<"; normal sum: " << buf[0] << " " << buf[1] << " " << buf[2] << "; total area: " << buf[3] << endl;
	}
      }
    }
  }
  

  // compute source terms

  /*
    ======================================================== 
    Original formulation: NO COFLOW  
  
  */
  // observer independent quantities if U0 = (0, 0, 0), i.e., no coflow
  void calc_sources_obs_indep() {  
    if (mpi_rank == 0){
      cout << " ***  calc_sources_obs_indep(), NO coflow "<< endl;
    }

    int nfa = faceVec.size();

    for (int ifa = 0; ifa < nfa; ++ifa) {
      // To reduce memory usage, area not stored anymore
      // normal vector contains information about area, namely n = n_unit * dA
      double nx = faceVec[ifa].nx;
      double ny = faceVec[ifa].ny;
      double nz = faceVec[ifa].nz;
     
      for (int it = 0; it < nt_total; ++it) {

	double pprime = faceVec[ifa].sourceQ[it];
	double u = faceVec[ifa].sourceF1[it];
	double v = faceVec[ifa].sourceF2[it];
	double w = faceVec[ifa].sourceF3[it];

	
	// compute density
	/*
	  ======================================================== 
	  Choose a formulation:	  
	  1) Formulation 1 (Default option): Spalart and Shur formulation, valid in acoustic region
	  rho = rho_ref + p'/c_ref^2 	  
	  
	  2) Formulation 2: Spalart and Shur generalized formulation, valid in isentropic region
	  rho = rho_ref (1 + p'/P_ref)^(1/gamma)
	  
	  3) Formulation 3: Original formulation, density directly from LES 
	*/ 
	double rho = rho_ref +  pprime/c2_ref;	

	//compute u.n
	double rho_un = rho * (u*nx + v*ny + w*nz);

	faceVec[ifa].sourceQ[it] = rho_un;
	faceVec[ifa].sourceF1[it] = pprime*nx + rho_un*u;
	faceVec[ifa].sourceF2[it] = pprime*ny + rho_un*v;
	faceVec[ifa].sourceF3[it] = pprime*nz + rho_un*w;

      }
    }    
  }


  // observer dependent quantities if U0 = (0, 0, 0), i.e., no coflow
  void calc_sources_obs_dep(double *x_obs) {       
    complex<double> p_hat_local;
   
    for (int ifa = 0; ifa < faceVec.size(); ++ifa) {
      // compute r
      double rx = x_obs[0] - faceVec[ifa].x;
      double ry = x_obs[1] - faceVec[ifa].y;
      double rz = x_obs[2] - faceVec[ifa].z;
      
      // compute rmag
      double rmag = sqrt(rx*rx + ry*ry + rz*rz);
      
      if ( rmag < 1.0e-10 ) 
	{
	  if (mpi_rank == 0)
	    cout << " !!! Error: singular observer ..." << endl;
	  throw(0);
	}
      
      // compute r_unit
      double r_unitx = rx/rmag;
      double r_unity = ry/rmag;
      double r_unitz = rz/rmag;
      
      for (int i_freq = 0; i_freq < num_freq; ++i_freq){
	// iomega now stores -i * omega  
	
	// Green function
	complex<double> G = -exp( iomega[i_freq]*rmag/c_ref) / rmag;
	
	complex<double> coeff_F = (-1.0/rmag + iomega[i_freq]/c_ref);
	
	// compute (-iw Qhat_n G - Fhat_i dG/dy_i) which is equivalent to -iw Qhat_n G + Fhat_i dG/dx_i
	p_hat_local = (faceVec[ifa].sourceQhat[i_freq]  * iomega[i_freq] + 
		      (faceVec[ifa].sourceF1hat[i_freq] * r_unitx + 
		       faceVec[ifa].sourceF2hat[i_freq] * r_unity + 
		       faceVec[ifa].sourceF3hat[i_freq] * r_unitz)* coeff_F) * G ;
      
	// ----------------------------
        // Zonal Noise Analysis (ZNA)
        // store the face contribution in all the different zone(s) containing that face
        for (int j=0; j < faceVec[ifa].ZNAsize; ++j){
	  int iZNA = faceVec[ifa].ZNAindex[j];
	  p_hat[iZNA][i_freq] += p_hat_local;
	}
      }
    }
  }


  // observer dependent quantities if U0 = (0, 0, 0), i.e., no coflow, AND solid surface formulation
  void calc_sources_obs_dep_solidFWH(double *x_obs) {       
    
    complex<double> p_hat_local;
   
    for (int ifa = 0; ifa < faceVec.size(); ++ifa) {
      // compute r
      double rx = x_obs[0] - faceVec[ifa].x;
      double ry = x_obs[1] - faceVec[ifa].y;
      double rz = x_obs[2] - faceVec[ifa].z;
      
      // compute rmag
      double rmag = sqrt(rx*rx + ry*ry + rz*rz);
      
      if ( rmag < 1.0e-10 ) 
	{
	  if (mpi_rank == 0)
	    cout << " !!! Error: singular observer ..." << endl;
	  throw(0);
	}
      
      // compute r_unit
      double r_unitx = rx/rmag;
      double r_unity = ry/rmag;
      double r_unitz = rz/rmag;

      // normal vector contains information about area, namely n = n_unit * dA
      double nx = faceVec[ifa].nx;
      double ny = faceVec[ifa].ny;
      double nz = faceVec[ifa].nz;
      
      for (int i_freq = 0; i_freq < num_freq; ++i_freq){
	// iomega now stores -i * omega  
	
	// Green function
	complex<double> G = -exp( iomega[i_freq]*rmag/c_ref) / rmag;
	
	complex<double> coeff_F = (-1.0/rmag + iomega[i_freq]/c_ref);
	
	// compute - Phat n_i dG/dy_i which is equivalent to Phat n_i dG/dx_i
	p_hat_local = (faceVec[ifa].sourceQhat[i_freq]  * 				   
			  (nx * r_unitx + 
			   ny * r_unity + 
			   nz * r_unitz)* coeff_F) * G ; 

	// ----------------------------
	// Zonal Noise Analysis (ZNA)
	// store the face contribution in all the different zone(s) containing that face
	for (int j=0; j < faceVec[ifa].ZNAsize; ++j){
	  int iZNA = faceVec[ifa].ZNAindex[j];
	  p_hat[iZNA][i_freq] += p_hat_local;
	}
      }
    }
  }



  /*
    ======================================================== 
    Wind-Tunnel formulation: WITH COFLOW  
  
  */
  // observer independent quantities if U0 = (U_COFLOW, 0, 0), i.e., with coflow
  void calc_sources_obs_indep_WT() {  
    if (mpi_rank == 0){
      cout << " ***  calc_sources_obs_indep_WT(), WITH coflow U0 = (" << U0[0] << ", "<< U0[1] << ", "<< U0[2] << ") " <<endl;
    }

    int nfa = faceVec.size();

    for (int ifa = 0; ifa < nfa; ++ifa) {
      // To reduce memory usage, area not stored anymore
      // normal vector contains information about area, namely n = n_unit * dA
      double nx = faceVec[ifa].nx;
      double ny = faceVec[ifa].ny;
      double nz = faceVec[ifa].nz;

      //compute U0.n, i.e. U0[0]*nx since U0[1]=U0[2]=0
      double U0n = U0[0] * nx;

      for (int it = 0; it < nt_total; ++it) {

	double pprime = faceVec[ifa].sourceQ[it];
	double u = faceVec[ifa].sourceF1[it] - U0[0];
	double v = faceVec[ifa].sourceF2[it];
	double w = faceVec[ifa].sourceF3[it];

	// compute density
	/*
	  ======================================================== 
	  Choose a formulation:	  
	  1) Formulation 1 (Default option): Spalart and Shur formulation, valid in acoustic region
	  rho = rho_ref + p'/c_ref^2 	  
	  
	  2) Formulation 2: Spalart and Shur generalized formulation, valid in isentropic region
	  rho = rho_ref (1 + p'/P_ref)^(1/gamma)
	  
	  3) Formulation 3: Original formulation, density directly from LES 
	*/ 
	double rho = rho_ref +  pprime/c2_ref;

	//compute u.n
	double un = u*nx + v*ny + w*nz;

	faceVec[ifa].sourceQ[it] = rho*(un + U0n) - rho_ref*U0n;
	faceVec[ifa].sourceF1[it] = pprime*nx + rho*(u - U0[0])*(un + U0n) + rho_ref*U0[0]*U0n;
	faceVec[ifa].sourceF2[it] = pprime*ny + rho*v*(un + U0n);
	faceVec[ifa].sourceF3[it] = pprime*nz + rho*w*(un + U0n);

      }
    }
  }


  
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  //
  // FOR FUTURE USE (i.e., COMPUTATION IN TIME DOMAIN) WITH OPTIONS P_TIME
  // Compute global first and last time of reception over array of observer
  //
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  //
  // calculate first and last time of reception at observer (for output in time domain)
  // valid for both with and without coflow
  //
  //----------------------------------------------------------------------------
  void calc_retarded_time(double *x_obs, double &my_tfirst, double &my_tlast) {       

    // first time of emission
    double taufirst = start*dt;
    
    // last time of emission
    double taulast = (start + (nt-1) * delta)*dt;
    
    //if (mpi_rank == 0){
    //  cout << "Observer at ( " << x_obs[0] << ", " << x_obs[1]<<", "<< x_obs[2] <<" ),  my_tfirst= " << my_tfirst <<",  my_tlast = " << my_tlast<< endl;
    //}
    
    for (int ifa = 0; ifa < faceVec.size(); ++ifa) {
      // compute Rstar
      double rx = x_obs[0] - faceVec[ifa].x;
      double ry = x_obs[1] - faceVec[ifa].y;
      double rz = x_obs[2] - faceVec[ifa].z;
      
      double Rstar = sqrt(rx*rx + beta2 * (ry*ry + rz*rz));
      if ( Rstar < 1.0e-10 ) 
	{
	  if (mpi_rank == 0)
	    cout << " !!! Error: singular observer ..." << endl;
	  throw(0);
	}
      
      // compute R
      double R = (-M0*rx+Rstar)/beta2;
      
      // compute first time of reception
      double tfirst = taufirst + R/c_ref;
      
      // compute last time of reception
      double tlast = taulast + R/c_ref;
      
      // find local minimum of reception time on proc    
      my_tfirst = max(my_tfirst, tfirst);
      
      // find local maximum of reception time on proc
      my_tlast=min(my_tlast,tlast);
    }    
  }
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  


  // observer dependent quantities if U0 = (U_COFLOW, 0, 0), i.e., with coflow
  void calc_sources_obs_dep_WT(double *x_obs) {       
    complex<double> p_hat_local;

    for (int ifa = 0; ifa < faceVec.size(); ++ifa) {
      // compute Rstar
      double rx = x_obs[0] - faceVec[ifa].x;
      double ry = x_obs[1] - faceVec[ifa].y;
      double rz = x_obs[2] - faceVec[ifa].z;
      
      double Rstar = sqrt(rx*rx + beta2 * (ry*ry + rz*rz));
      if ( Rstar < 1.0e-10 ) 
	{
	  if (mpi_rank == 0)
	    cout << " !!! Error: singular observer ..." << endl;
	  throw(0);
	}
      
      // compute R
      double R = (-M0*rx+Rstar)/beta2;
      
      // compute dRstar/dyi
      double dRstar_x = -rx/Rstar;
      double dRstar_y = -beta2*ry/Rstar;
      double dRstar_z = -beta2*rz/Rstar;

      // compute dR/dyi
      double dR_x = (M0+dRstar_x)/beta2;
      double dR_y = dRstar_y/beta2;
      double dR_z = dRstar_z/beta2;

      for (int i_freq = 0; i_freq < num_freq; ++i_freq){
	// iomega now stores -i * omega  

	// Green function
	complex<double> G = -exp( iomega[i_freq]*R/c_ref) / Rstar;

	// compute (-iw Qhat_n G - Fhat_i dG/dy_i)
	p_hat_local = (faceVec[ifa].sourceQhat[i_freq]  * iomega[i_freq]	
			  - faceVec[ifa].sourceF1hat[i_freq] * (-dRstar_x/Rstar + iomega[i_freq]/c_ref*dR_x) 
                          - faceVec[ifa].sourceF2hat[i_freq] * (-dRstar_y/Rstar + iomega[i_freq]/c_ref*dR_y) 
			  - faceVec[ifa].sourceF3hat[i_freq] * (-dRstar_z/Rstar + iomega[i_freq]/c_ref*dR_z))* G ;
	
	// ----------------------------
        // Zonal Noise Analysis (ZNA)
        // store the face contribution in all the different zone(s) containing that face
        for (int j=0; j < faceVec[ifa].ZNAsize; ++j){
	  int iZNA = faceVec[ifa].ZNAindex[j];
	  p_hat[iZNA][i_freq] += p_hat_local;
	}
      }
    }
  }


  // observer dependent quantities if U0 = (U_COFLOW, 0, 0), i.e., with coflow
  void calc_sources_obs_dep_WT_solidFWH(double *x_obs) {       
    
    complex<double> p_hat_local;
        
    for (int ifa = 0; ifa < faceVec.size(); ++ifa) {
      // compute Rstar
      double rx = x_obs[0] - faceVec[ifa].x;
      double ry = x_obs[1] - faceVec[ifa].y;
      double rz = x_obs[2] - faceVec[ifa].z;
      
      double Rstar = sqrt(rx*rx + beta2 * (ry*ry + rz*rz));
      if ( Rstar < 1.0e-10 ) 
	{
	  if (mpi_rank == 0)
	    cout << " !!! Error: singular observer ..." << endl;
	  throw(0);
	}
      
      // compute R
      double R = (-M0*rx+Rstar)/beta2;
      
      // compute dRstar/dyi
      double dRstar_x = -rx/Rstar;
      double dRstar_y = -beta2*ry/Rstar;
      double dRstar_z = -beta2*rz/Rstar;

      // compute dR/dyi
      double dR_x = (M0+dRstar_x)/beta2;
      double dR_y = dRstar_y/beta2;
      double dR_z = dRstar_z/beta2;
      
      // normal vector contains information about area, namely n = n_unit * dA
      double nx = faceVec[ifa].nx;
      double ny = faceVec[ifa].ny;
      double nz = faceVec[ifa].nz;

      for (int i_freq = 0; i_freq < num_freq; ++i_freq){
	// iomega now stores -i * omega  

	// Green function
	complex<double> G = -exp( iomega[i_freq]*R/c_ref) / Rstar;
	
	// compute - Phat n_i dG/dy_i
	p_hat_local = -faceVec[ifa].sourceQhat[i_freq] * 
	    (nx * (-dRstar_x/Rstar + iomega[i_freq]/c_ref*dR_x) 
	   + ny * (-dRstar_y/Rstar + iomega[i_freq]/c_ref*dR_y) 
	   + nz * (-dRstar_z/Rstar + iomega[i_freq]/c_ref*dR_z))* G ;

	// ----------------------------
        // Zonal Noise Analysis (ZNA)
        // store the face contribution in all the different zone(s) containing that face
        for (int j=0; j < faceVec[ifa].ZNAsize; ++j){
	  int iZNA = faceVec[ifa].ZNAindex[j];
	  p_hat[iZNA][i_freq] += p_hat_local;
	}
      }
    }
  }


  void updateSources(){
    // Option to turn-off thickness and/or loading noise contributions
    if (!solidFWH){
      // thickness (monopole) noise
      if (QF_switch[0] == 0.0){
	for (int ifa = 0; ifa < faceVec.size(); ++ifa) {
	  for (int i_freq = 0; i_freq < num_freq; ++i_freq){
	    faceVec[ifa].sourceQhat[i_freq] = 0.0;
	  }
	}
	if (mpi_rank == 0){
	  cout << " > WARNING: thickness noise not computed" << endl;
	}
      }

      // loading (dipole) noise, x-component
      if (QF_switch[1] == 0.0){
	for (int ifa = 0; ifa < faceVec.size(); ++ifa) {
	  for (int i_freq = 0; i_freq < num_freq; ++i_freq){
	    faceVec[ifa].sourceF1hat[i_freq] = 0.0;
	  }
	}
	if (mpi_rank == 0){
	  cout << " > WARNING: loading noise (x-component) not computed" << endl;
	}
      }
 
      // loading (dipole) noise, y-component
      if (QF_switch[2] == 0.0){
	for (int ifa = 0; ifa < faceVec.size(); ++ifa) {
	  for (int i_freq = 0; i_freq < num_freq; ++i_freq){
	    faceVec[ifa].sourceF2hat[i_freq] = 0.0;
	  }
	}
	if (mpi_rank == 0){
	  cout << " > WARNING: loading noise (y-component) not computed" << endl;
	}
      }

      // loading (dipole) noise, z-component
      if (QF_switch[3] == 0.0){
	for (int ifa = 0; ifa < faceVec.size(); ++ifa) {
	  for (int i_freq = 0; i_freq < num_freq; ++i_freq){
	    faceVec[ifa].sourceF3hat[i_freq] = 0.0;
	  }
	}
	if (mpi_rank == 0){
	  cout << " > WARNING: loading noise (z-component) not computed" << endl;
	}
      } 
    }
  }
 

  void initFFT() {    
    if (!solidFWH){
      for (int ifa = 0; ifa < faceVec.size(); ++ifa) { 
	faceVec[ifa].sourceF1hat = new complex<double>[num_freq];
	faceVec[ifa].sourceF2hat = new complex<double>[num_freq];
	faceVec[ifa].sourceF3hat = new complex<double>[num_freq];
      }
    }
    
    for (int ifa = 0; ifa < faceVec.size(); ++ifa) { 
      faceVec[ifa].sourceQhat = new complex<double>[num_freq];
    } 
  }


  void fftSources(int block_start) {

    assert(block_start+nt-1<=nt_total);

    complex<double> (*output) = new complex<double>[num_freq_total];
    double (*input) = new double[nt];
    
    fftw_plan plan;
    plan = fftw_plan_dft_r2c_1d(nt, input,  reinterpret_cast<fftw_complex*>(output),FFTW_MEASURE);
    
    if (!solidFWH){
      
      for (int ifa = 0; ifa < faceVec.size(); ++ifa) { 
	
	// calculate the mean of sourceF
	double meanF1 = 0.0;
	double meanF2 = 0.0;
	double meanF3 = 0.0;
	for (int it=0; it<nt; it++) {
	  meanF1 += faceVec[ifa].sourceF1[block_start + it];
	  meanF2 += faceVec[ifa].sourceF2[block_start + it];
	  meanF3 += faceVec[ifa].sourceF3[block_start + it];
	}
	meanF1 /= ((double)nt);
	meanF2 /= ((double)nt);
	meanF3 /= ((double)nt);
	
	// calculate fft of sourceF1 ------------- 
	// replace input
	for (int it = 0; it < nt; ++it){ 
	  input[it] = window[it] * (faceVec[ifa].sourceF1[block_start + it] - meanF1);
	}	  
	fftw_execute(plan);
	// replace output  
	for (int i_freq = 0; i_freq < num_freq; ++i_freq){ 
	  faceVec[ifa].sourceF1hat[i_freq] = output[i_freq];
	}
	
	// calculate fft of sourceF2 ------------- 
	// replace input
	for (int it = 0; it < nt; ++it){ 
	  input[it] = window[it] * (faceVec[ifa].sourceF2[block_start + it] - meanF2);
	}	  
	fftw_execute(plan);
	// replace output  
	for (int i_freq = 0; i_freq < num_freq; ++i_freq){ 
	  faceVec[ifa].sourceF2hat[i_freq] = output[i_freq];
	}
	
	// calculate fft of sourceF3 ------------- 
	// replace input
	for (int it = 0; it < nt; ++it){ 
	  input[it] = window[it] * (faceVec[ifa].sourceF3[block_start + it] - meanF3);
	}
	fftw_execute(plan);	  
	// replace output  
	for (int i_freq = 0; i_freq < num_freq; ++i_freq){ 
	  faceVec[ifa].sourceF3hat[i_freq] = output[i_freq];
	}
      }
    }
    
    
    for (int ifa = 0; ifa < faceVec.size(); ++ifa) { 
      
      // calculate the mean of sourceQ
      double meanQ = 0.0;
      for (int it=0; it<nt; it++) {
	meanQ  += faceVec[ifa].sourceQ[block_start + it];
      }
      meanQ  /= ((double)nt);
      
      // calculate fft of sourceQ --------------      
      // replace input
      for (int it = 0; it < nt; ++it){ 
	input[it] = window[it] * (faceVec[ifa].sourceQ[block_start + it] - meanQ);
      }
      fftw_execute(plan);
      // replace output  
      for (int i_freq = 0; i_freq < num_freq; ++i_freq){ 
	faceVec[ifa].sourceQhat[i_freq] = output[i_freq]; 
      }
    }
    
    // clean up
    fftw_destroy_plan(plan);
    DELETE(input);
    DELETE(output);
    
  }
  
  
  // normalize and sum accross processors
  void SurfaceIntegration() {
    
    // 1/4Pi factor and FFTW normalization 
    for (int i_freq = 0; i_freq < num_freq; ++i_freq){
      for (int iZNA=0; iZNA<nZNA; ++iZNA){
	p_hat[iZNA][i_freq] /= (4.0*M_PI*((double)nt));
      }  
    }

    // Sum all
     MPI_Reduce(&p_hat[0][0],&p_hat_global[0][0],2*num_freq*nZNA,MPI_DOUBLE,MPI_SUM,0,mpi_comm);
  }

  

  void initRotation() {
    if ((convectiveFWH) && ((U0[1]!= 0.0) || (U0[2]!= 0.0))){
      coord_rot = true;
      double beta = 0.0;
      double  gam = 0.0;
      if ( U0[0]!= 0){
	// right-hand rule rotation about z axis
	gam = atan(U0[1]/U0[0]);
	double cg = cos(gam);
	double sg = sin(gam);
	// after Rz applied, right-hand rule rotation about new -y' axis (hence - sign)
	double RzU0 = cg*U0[0]+sg*U0[1];
	beta = -atan(U0[2]/RzU0);
	double cb = cos(beta);
	double sb = sin(beta);
	if (mpi_rank == 0){
	  if (gam != 0.0)  cout << " > Initializing coordinates rotation about z-axis of angle  " <<  gam*180/M_PI << " deg." << endl;
	  if (beta != 0.0) cout << " > Initializing coordinates rotation about y'-axis of angle " << beta*180/M_PI << " deg." << endl;
	}
	// Final rotation matrix
	RyRz[0][0]=cb*cg;  RyRz[0][1]=cb*sg;  RyRz[0][2]=-sb;
	RyRz[1][0]=  -sg;  RyRz[1][1]=   cg;  RyRz[1][2]=0.0;
	RyRz[2][0]=sb*cg;  RyRz[2][1]=sb*sg;  RyRz[2][2]= cb;
      }
      else{
	if (mpi_rank == 0){
	  cout << " > Cannot apply coordinates rotation ... Main coflow direction should be X " << endl;
	}
	throw(0);
      }
    }
    else{
      coord_rot = false;
    }
  }



  void applyRotation(const double M[3][3], double &x0, double &x1, double &x2) {
    double v_orig[3];
    
    // build vector
    v_orig[0]= x0;
    v_orig[1]= x1;
    v_orig[2]= x2;
    
    // apply rotation
    double v_rot[3]=MATRIX_PRODUCT(M,v_orig);
    
    // save result
    x0 = v_rot[0];
    x1 = v_rot[1];
    x2 = v_rot[2];
  }



  void applyRotation(const double M[3][3], double x[3]) {
    
    // apply rotation
    double v_rot[3]=MATRIX_PRODUCT(M,x);
    
    // save result
    x[0] = v_rot[0];
    x[1] = v_rot[1];
    x[2] = v_rot[2];
  }
};  

