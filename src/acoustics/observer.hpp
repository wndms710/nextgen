#include "CTI.hpp"
using namespace CTI;
#include <complex>
#include <vector>
#include "FFT.hpp" // wrapper for <fftw3.h>
#include "MiscUtils.hpp"

// ==============================================================
//
// The observers (i.e., microphones) are grouped into arrays, which
// can be:
//   - a single point (see function "addPoint")
//   - an arc (see function "addPointArc")
//   - a rotating arc (see function "addArcMesh")
//   - a cartesian mesh (see function "addCartMesh")
//
// ==============================================================
class Observer {

private:

  int nt, num_freq, num_obs, this_obs, num_obs_hom, num_array, this_array;

  // ----------------------------
  // Zonal Noise Analysis (ZNA)
  int nZNA;                    // number of zone(s)
  int ZNAcoord;                // 0 for x, 1 for y, 2 for z 
  double (*ZNAlimit)[2];       // min & max limit for each zone

  int num_freq_total;
  vector<double> x_obs;        // positions of all the observers
  vector<double> x_obs_hom;    // positions of the homogeneous observer
  vector<int> obs_index;       // indices of observers (size num_obs)
  vector<int> obs_hom_index;   // indices of homogeneous observers (size num_obs_hom)
  
  
  vector<string> array_name;   // observer array name
  vector<int> array_m1;        // array size 1 (for PLOT3D)
  vector<int> array_m2;        // array size 2 (for PLOT3D)
  vector<int> array_m3;        // array size 3 (for PLOT3D) ... not needed ?
  vector<bool> array_PLOT3D;   // Flag to export array in PLOT3D format
  vector<bool> array_PLOTtime; // Flag to export array of p(t) (in time domain)
  vector<bool> array_PLOThat;  // Flag to export array of p_hat(f) (in frequency domain)
  vector<int> array_start;     // start of array in obs_index
  vector<int> array_end;       // end of array in obs_index
  vector<int> array_hom_start; // start of array in obs_hom_index
  vector<int> array_hom_end;   // end of array in obs_hom_index
  vector<int> array_hom_it_first; // index of global first time of reception for array
  vector<int> array_hom_it_last;  // index of global last time of reception for array
  // *******************************************
  // Update on angles alpha & beta calculation 
  // Ability to specify center of coordinate system different than default (0, 0, 0) for angle definition
  vector<double> array_center1;   // x-coordinates of center of coordinate system of array (for angles alpha & beta calculation)
  vector<double> array_center2;   // y-coordinates of center of coordinate system of array (for angles alpha & beta calculation)
  vector<double> array_center3;   // z-coordinates of center of coordinate system of array (for angles alpha & beta calculation)

  // ----------------------------
  // Zonal Noise Analysis (ZNA)
  double * p;                   // acoustic pressure (size nt*num_obs)
  complex<double> ** p_hat;     // Fourier transform of p (size (nZNA,num_freq*num_obs))
  double * f;                   // frequency (size num_freq)
  double * oaspl_f;             // frequency for OASPL calculation (size num_freq)
  double ** psd;                // power spectral density (size (nZNA, num_freq*num_obs))
  double ** oaspl;              // overall sound pressure levels (size (nZNA,num_obs))
  double ** psd_hom;            // power spectral density for homogeneous observers (size (nZNA,num_freq*num_obs_hom))
  double ** oaspl_hom;          // overall sound pressure levels for homogeneous observers (size (nZNA,num_obs_hom))
  
  bool PSD_calc, time_calc;

  double oaspl_f_min;          // OASPL computed for frequencies oaspl_f_min <= f <= oaspl_f_max
  double oaspl_f_max;

  // ----------------------------
  // FFT overlap
  int nblock; // number of blocks with requested bandwidth and overlap
  int nblock_t; // total number of blocks, including average over blocks (i.e., nblock if only 1 block, nblock + 1 otherwise)


  // add z-aligned arc
  int axis; // 1 for x-axis (default), 3 for z-axis, y-axis not implemented yet
 
public:
  
  double p_char;               // coefficient to convert p to characteristic units (if setup in input file)
  double f_char;               // coefficient to convert f to characteristic units (if setup in input file)
  double df_char;              // coefficient to convert df to characteristic units (if setup in input file)
  double p2_rel;               // Reference pressure for dB levels, p2_rel = (20 10-6)^2 Pa^2
      
  double df;

  int dimension;               // Current implementation of 3D for FW-H

  bool thetaX_hom;


  // ----------------------------
  // Zonal Noise Analysis (ZNA)
  // Initialize Observer with input _in passed from class FWH
  Observer(int nZNA_in, int ZNAcoord_in, double (*ZNAlimit_in)[2], int nt_in, int num_freq_in, double *omega_in, double oaspl_f_min_in, double oaspl_f_max_in, double p_char_in, double f_char_in, double df_char_in, int dimension_in, bool thetaX_hom_in, int nblock_in, int axis_in){

    //initialize p_hat, p and psd    
    p = NULL;
    p_hat = NULL;
    psd = NULL;
    oaspl = NULL;
    psd_hom = NULL;
    oaspl_hom = NULL;
    
    axis = axis_in;
    nZNA = nZNA_in;
    ZNAcoord = ZNAcoord_in;
    ZNAlimit = NULL;
    if (nZNA > 1){
      ZNAlimit=new double[nZNA][2];
      for (int iZNA=0; iZNA < nZNA; ++iZNA){
	FOR_I2 ZNAlimit[iZNA][i] = ZNAlimit_in[iZNA][i];
      }
    }

    // copy info
    nt = nt_in;
    num_freq = num_freq_in;
    oaspl_f_min = oaspl_f_min_in;
    oaspl_f_max = oaspl_f_max_in;
    p_char = p_char_in;
    f_char = f_char_in;
    df_char = df_char_in;

    num_freq_total = nt/2 + 1; // total number of frequencies accessible by postprocessing

    nblock = nblock_in;
    if (nblock == 1) {
      nblock_t = nblock;
    }
    else{
      nblock_t = nblock + 1;
    }

    // initialize bandwidth for PSD
    df = (omega_in[1] - omega_in[0]) * df_char /(2.0*M_PI); // convert df to characteristic quantity
    if (mpi_rank == 0){
      cout << " > Frequency bin for narrowband PSD df = " << df << endl;
    }
    
    // initialize frequencies
    f = new double[num_freq];
    oaspl_f = new double[num_freq];
    for (int ifreq=0; ifreq  < num_freq; ++ifreq){
      f[ifreq] =  f_char*omega_in[ifreq]/(2.0*M_PI); // convert f to characteristic quantity
    }
    
    // initialize frequencies range for oaspl
    for (int ifreq=0; ifreq  < num_freq; ++ifreq){
      if ((f[ifreq]>=oaspl_f_min) && (f[ifreq]<=oaspl_f_max)){
	oaspl_f[ifreq] = 1.0;
      }
      else {
	oaspl_f[ifreq] = 0.0;
      }
    }
    
    // initialize the rest ...
    num_obs = 0;
    this_obs = 0;
    num_obs_hom = 0;
    num_array=0;
    this_array=0;

    PSD_calc = false;
    time_calc = false;
    thetaX_hom = thetaX_hom_in;

    dimension = dimension_in;

    p2_rel =4.0E-10;

  }
  


  virtual ~Observer(){
    if (p != NULL) 
      delete[] p;
            
    if (f != NULL) 
      delete[] f;
    
    if (oaspl_f != NULL) 
      delete[] oaspl_f;

    // ----------------------------
    // Zonal Noise Analysis (ZNA)
    if (ZNAlimit != NULL) 
      delete[] ZNAlimit;

    if (p_hat != NULL)
      MiscUtils::deleteComplexArray2d(p_hat);

    if (psd != NULL) 
      MiscUtils::deleteArray2d(psd);
    
    if (oaspl != NULL) 
      MiscUtils::deleteArray2d(oaspl);

    if (psd_hom != NULL) 
      MiscUtils::deleteArray2d(psd_hom);
    
    if (oaspl_hom != NULL) 
      MiscUtils::deleteArray2d(oaspl_hom);  
  }


  /* ==========================================================
  
  Add a single microphone
  ------------------------
  addPoint(xp,name,xc)
     
  xp[3]: coordinate of the microphone
  xc[3]: location of the center of the coordinate system (for calculation of the angles alpha & beta reported in output)
     
  This function  will add one (1) single observer
  */    
  void addPoint(double *x, string name, double *xc) {   
    // add observer
    FOR_I3 x_obs.push_back(x[i]);
    obs_index.push_back(num_obs);

    // store array information
    array_m1.push_back(1);
    array_m2.push_back(1);
    array_m3.push_back(1);
    array_start.push_back(num_obs);
    array_name.push_back(name);
    array_hom_start.push_back(num_obs_hom);
    // *******************************************
    // Update on angles alpha & beta calculation
    // Store center of coordinate system
    array_center1.push_back(xc[0]);
    array_center2.push_back(xc[1]);
    array_center3.push_back(xc[2]);
    
    num_obs++;

    num_array++;
    
    array_end.push_back(num_obs);
    array_hom_end.push_back(num_obs_hom);

    // DO NOT add homogenous observer
    obs_hom_index.push_back(-1);
  }
  


  /* ==========================================================
     Add an array of microphones with uniform azimuthal (theta) spacing. 
     ------------------------   
     addPointArc(xa, ra, ntheta, name, xc, theta_offset)
     
     xa : axial location of microphones
     ra : distance of microphones from the x-axis
     ntheta : number of azimuthal microphones
     xc[3]: location of the center of the coordinate system (for calculation of the angles alpha & beta reported in output)
     theta_offset[ntheta]: array of offset for theta in deg (default: all set to 0)

     theta: if axis = 1, azimuthal angle measured with respect to +x axis 
     (positive in clockwise direction, i.e., theta=90 along -z axis)
     theta: if axis =3, azimuthal angle measured with respect to +z axis 
     (positive in clockwise direction, i.e., theta=90 along -x axis)
     Range:  0 <= theta <= 360 deg
     
     This function will add one (1) homogeneous observer and (ntheta) single observers
  */   
  void addPointArc(double xa, double r, int ntheta, string name, double *xc, double *theta_offset) {  
    double x[3];

    if ((dimension == 2)  && (ntheta > 1)){
      if (mpi_rank == 0){
	cout << "WARNING ... 2D problem ! Parameter ntheta = " << ntheta << " discarded, only one observer in array "<< name << endl;
      }
      ntheta = 1;    
    }
  
    if (thetaX_hom){
      // add homogenous observer
      // Arc along z-axis
      if (axis ==3) {
	x_obs_hom.push_back(0.0);
	x_obs_hom.push_back(r);
	x_obs_hom.push_back(xa);
      }
      // Arc along x-axis
      else {
	x_obs_hom.push_back(xa);
	x_obs_hom.push_back(r);
	x_obs_hom.push_back(0.0);
      }
    }
      
    // store array information
    array_m1.push_back(ntheta);
    array_m2.push_back(1);
    array_m3.push_back(1);
    array_start.push_back(num_obs);
    array_name.push_back(name);
    array_hom_start.push_back(num_obs_hom);
    // *******************************************
    // Update on angles alpha & beta calculation
    // Store center of coordinate system
    array_center1.push_back(xc[0]);
    array_center2.push_back(xc[1]);
    array_center3.push_back(xc[2]);
    

    for (int itheta=0; itheta<ntheta ; itheta++){
      double theta = (((double)itheta)/((double)ntheta) + theta_offset[itheta]/360.0) * 2.0 * M_PI;
      if (axis == 3) {
	x[0] = - r*sin(theta);
	x[1] =   r*cos(theta);
	x[2] = xa;
      }
      else {
	x[0] = xa;
	x[1] =   r*cos(theta);
	x[2] = - r*sin(theta);
      }
      // add observer
      FOR_I3 x_obs.push_back(x[i]);
      obs_index.push_back(num_obs);
      
      if (thetaX_hom) obs_hom_index.push_back(num_obs_hom);
      num_obs++;
    }

    num_array++;        

    array_end.push_back(num_obs);
      
    if (thetaX_hom){
      num_obs_hom++;
    }
    else{
      // DO NOT add homogenous observer
      obs_hom_index.push_back(-1);
    }

    array_hom_end.push_back(num_obs_hom);
  }



  /* ==========================================================  
     Add a rotating arc of microphones
     ------------------------
     addRotatingArc(x_orig, r , alpha_min, alpha_max, nalpha, ntheta, name, xc, theta_offset)
     
     x_orig[3]: coordinates of the center of the arc
     r: radius of the arc
     alpha_min, alpha_max: minimum and maximum INLET angles of the arc microphones in degrees
     nalpha: number of microphones on the arc
     ntheta: number of azimuthal stations
     xc[3]: location of the center of the coordinate system (for calculation of the angles alpha & beta reported in output)
     theta_offset[ntheta]: array of offset for theta in deg (default: all set to 0)
     
     alpha: inlet angle measured with respect to +z axis 
     (positive in clockwise direction, i.e., alpha=180 downstream of jet)
     Range: 0 <= alpha <= 180 deg

     theta: azimuthal angle measured with respect to +x axis 
     (positive in clockwise direction, i.e., theta=90 along -z axis)
     Range: 0 <= theta <= 360 deg
     
     This function will add (nalpha) homogeneous observers and (nalpha*ntheta) single observers 
  */ 
  void addRotatingArc(double *x_orig, double r, double alpha_min, double alpha_max, int nalpha, int ntheta, string name, double *xc, double *theta_offset) {

    if ((dimension == 2) && (ntheta > 1)){
      if (mpi_rank == 0){
	cout << "WARNING ... 2D problem ! Parameter ntheta = " << ntheta << " discarded, only " << nalpha << " observers in array "<< name  << endl;
      }
      ntheta = 1;    
    }
    
    // store array information
    array_m1.push_back(ntheta); 
    array_m2.push_back(nalpha);
    array_m3.push_back(1);
    array_start.push_back(num_obs);
    array_name.push_back(name);
    array_hom_start.push_back(num_obs_hom);
    // *******************************************
    // Update on angles alpha & beta calculation
    // Store center of coordinate system
    array_center1.push_back(xc[0]);
    array_center2.push_back(xc[1]);
    array_center3.push_back(xc[2]);

      
    for (int ialpha=0; ialpha<nalpha ; ialpha++){
      
      double alpha = (alpha_min + ((double)ialpha)/((double)(nalpha - 1)) * (alpha_max - alpha_min)) * M_PI / 180.0;
      
      double xa = x_orig[0] - r * cos(alpha); 
      double ra = x_orig[1] + r * sin(alpha); 

      if (thetaX_hom){
	// add homogenous observer
	x_obs_hom.push_back(xa);
	x_obs_hom.push_back(ra);
	x_obs_hom.push_back(0.0);
      }


      for (int itheta=0; itheta<ntheta ; itheta++){
	double theta = (((double)itheta)/((double)ntheta) + theta_offset[itheta]/360.0) * 2.0 * M_PI;
	
	// add observer
	x_obs.push_back(xa);              // x
        x_obs.push_back(ra*cos(theta));   // y
        x_obs.push_back(-ra*sin(theta));  // z
	obs_index.push_back(num_obs);
	if (thetaX_hom) obs_hom_index.push_back(num_obs_hom);
	
	num_obs++;
      }

      if (thetaX_hom) num_obs_hom++;  
    }

    num_array++;

    array_end.push_back(num_obs);
    array_hom_end.push_back(num_obs_hom);

    if (!thetaX_hom){
      // DO NOT add homogenous observer
      obs_hom_index.push_back(-1);
    }
  }
  


  /* ==========================================================
     Add a cartesian grid of microphones
     ------------------------
     addCartMesh(x0,x1,x2,n01,n02,xc)
     
     x0[3],x1[3],x2[3]: coordinates of the rectangle nodes
     n01,n02: number of intervals on each edge
     xc[3]: location of the center of the coordinate system (for calculation of the angles alpha & beta reported in output)

     
     x2 -----------------------------
     |                               |
     |                               |
    n02                              |
     |                               |
     |                               |
     x0 ------------ n01 ----------- x1
         
     This function will add (n01*n02) single observers     
  */ 
  void addCartMesh(double *x0, double *x1, double *x2, int n01, int n02, string name, double *xc) {
    
    double xi[3], xf[3], x[3];

    if (dimension == 2){
      if (mpi_rank == 0){
	cout << "WARNING ... 2D problem ! Spanwise z=0 for all observers in array "<< name  << endl;
      }
      x0[2] = 0.0;
      x1[2] = 0.0; 
      x2[2] = 0.0;
    }
    
    // store array information
    array_m1.push_back(n02); 
    array_m2.push_back(n01);
    array_m3.push_back(1);
    array_start.push_back(num_obs);
    array_name.push_back(name);
    array_hom_start.push_back(num_obs_hom);
    // *******************************************
    // Update on angles alpha & beta calculation
    // Store center of coordinate system
    array_center1.push_back(xc[0]);
    array_center2.push_back(xc[1]);
    array_center3.push_back(xc[2]);

    for (int i0=0; i0< n01; ++i0){
      FOR_I3 xi[i] = x0[i] + ((double)i0)/((double)(n01 - 1)) * (x1[i] - x0[i]);
      FOR_I3 xf[i] = x2[i] + ((double)i0)/((double)(n01 - 1)) * (x1[i] - x0[i]);
         
      for (int i1=0; i1< n02; ++i1){	
	FOR_I3 x[i] = xi[i] + ((double)i1)/((double)(n02 - 1)) * (xf[i] - xi[i]);

	// add observer
	FOR_I3 x_obs.push_back(x[i]);
	obs_index.push_back(num_obs);

	// DO NOT add homogenous observer
	obs_hom_index.push_back(-1);
	
	num_obs++;
      }
    }

    num_array++;

    array_hom_end.push_back(num_obs_hom);
    array_end.push_back(num_obs);  
  }
  

  
  void finalize_locations() {
    if (mpi_rank == 0){
          
      for (int i=0; i<num_array; ++i){
	cout << " > Observer (Microphone) array "<< i+1 <<", name: " << array_name[i] << endl;
	cout << " >          number of observer(s):" << array_end[i] - array_start[i] << endl;
	if (array_hom_end[i] - array_hom_start[i] > 0){
	  cout << " >          number of homogeneous observer(s): " << array_hom_end[i] - array_hom_start[i]<< endl;
	}
      }
      cout << " > Total number of observer(s): " << num_obs << endl;
      cout << " > Total number of homogeneous observer(s): " << num_obs_hom << endl;
      cout << " > Number of frequencies: " << num_freq << endl;
    
    
      // ----------------------------
      // Zonal Noise Analysis (ZNA)
      p_hat = MiscUtils::newComplexArray2d(nZNA*nblock,num_freq*num_obs);
      psd = MiscUtils::newArray2d(nZNA*nblock_t,num_freq*num_obs);
      oaspl = MiscUtils::newArray2d(nZNA*nblock_t,num_obs);
      for (int i=0; i < nZNA*nblock_t; ++i){
	for (int j=0; j < num_obs; ++j){
	  oaspl[i][j] = 0.0;
	}
	for (int j=0; j < num_freq*num_obs; ++j){
	  psd[i][j] = 0.0;
	}
      }
      psd_hom = MiscUtils::newArray2d(nZNA*nblock_t,num_freq*num_obs_hom);
      oaspl_hom = MiscUtils::newArray2d(nZNA*nblock_t,num_obs_hom);
    }
  }
    
    
  


  int get_num_obs() {
    return num_obs;
  }

  int get_num_obs_hom() {
    return num_obs_hom;
  }

  int get_num_array() {
    return num_array;
  }
    
  void inc_this_obs() {
    this_obs++;
  }

  void set_this_obs(int i) {
    this_obs = i;
  }
 
  void set_this_array(int i) {
    this_array = i;
  }

  void inc_this_array() {
    this_array++;
  }
  
  void get_x_obs(double *x) {
    FOR_I3 x[i] = x_obs[3*this_obs+i];
  }

  bool get_array_PLOT3D() {
    return array_PLOT3D[this_array];
  }

  int get_array_start() {
    return array_start[this_array];
  }
  
  int get_array_end() {
    return array_end[this_array];
  }
  
  string get_array_name() {
    return array_name[this_array];
  }

  void set_array_PLOT3D(bool flag_value) {
    array_PLOT3D.push_back(flag_value);
  }

  bool get_array_PLOTtime() {
    return array_PLOTtime[this_array];
  }

  void set_array_PLOTtime(bool flag_value) {
    array_PLOTtime.push_back(flag_value);
  }

  bool get_array_PLOThat() {
    return array_PLOThat[this_array];
  }

  void set_array_PLOThat(bool flag_value) {
    array_PLOThat.push_back(flag_value);
  }
  
  void set_array_it_first(int it_first) {
    array_hom_it_first.push_back(it_first);
  }

  void set_array_it_last(int it_last) {
    array_hom_it_last.push_back(it_last);
  }

  int get_array_it_first() {
    return array_hom_it_first[this_array];
  }

  int get_array_it_last() {
    return array_hom_it_last[this_array];
  }

  
  // ----------------------------
  // Zonal Noise Analysis (ZNA)
  void save_p_hat(complex<double> (*p_hat_global)) {
    if (mpi_rank == 0){ 
      for (int ifreq = 0; ifreq < num_freq; ++ifreq){ 
	p_hat[0][num_freq*this_obs + ifreq] = p_char*p_hat_global[ifreq]; // convert p_hat to characteristic quantity
      }
    }
  }
  

  void save_p_hat(complex<double> (**p_hat_global), int iblock) {
    if (mpi_rank == 0){
       for (int ifreq = 0; ifreq < num_freq; ++ifreq){ 
	 for (int j=0; j<nZNA; ++j){
	   p_hat[nZNA*iblock + j][num_freq*this_obs + ifreq] = p_char*p_hat_global[j][ifreq]; // convert p_hat to characteristic quantity
	 }
      }
    }
  }


  // ----------------------------
  // Zonal Noise Analysis (ZNA)
  void calcPSDandOASPL(int iblock) {   
    if (mpi_rank == 0) {
      PSD_calc = true;

      complex<double> psd_c;

      for (int iZNA=0; iZNA<nZNA; ++iZNA){
	
	// for OASPL calculation, oaspl_f[ifreq] = 1 if frequency used, 0 otherwise
	for (int ifreq=0; ifreq<num_freq; ++ifreq) {
	  psd_c = p_hat[nZNA*iblock + iZNA][num_freq*this_obs + ifreq] * conj(p_hat[nZNA*iblock + iZNA][num_freq*this_obs + ifreq]);
	  // average over block (
	  psd  [iZNA][num_freq*this_obs+ifreq] += 2.0 * psd_c.real() / df / double(nblock);
	  oaspl[iZNA][this_obs] += 2.0 * oaspl_f[ifreq] * psd_c.real() / double(nblock) ;
	  // block only
	  if (nblock > 1){
	    psd  [nZNA*(iblock+1)+iZNA][num_freq*this_obs+ifreq] = 2.0 * psd_c.real() / df ;
	    oaspl[nZNA*(iblock+1)+iZNA][this_obs] += 2.0 * oaspl_f[ifreq] * psd_c.real() ;
	  }   
	}
      } 
    } 
  } 


  // ----------------------------
  // Zonal Noise Analysis (ZNA)
  void calcPSDandOASPL_hom() {
    if (mpi_rank == 0) {
      if (PSD_calc == false)
	cout << " > Warning : cannot calculate homogeneous PSD and OASPL. Skipping calcPSDandOASPL_hom() " << endl;
      else{
	int *bin_size = new int[num_obs_hom];
	for (int iblock=0; iblock<nblock_t; ++iblock){
	  for (int iZNA=0; iZNA<nZNA; ++iZNA){
	    // initialize
	    for (int iobs_hom=0; iobs_hom<num_obs_hom; ++iobs_hom){	  
	      bin_size[iobs_hom] = 0;
	      oaspl_hom[nZNA*iblock + iZNA][iobs_hom] = 0.0;	  
	      for (int ifreq=0; ifreq<num_freq; ++ifreq){ 
		psd_hom[nZNA*iblock + iZNA][ifreq + num_freq*iobs_hom] = 0.0;
	      }
	    }
	    
	    // sum
	    int iobs_hom = obs_hom_index[0];
	    for (int iobs=0; iobs<num_obs; ++iobs){
	      if (obs_hom_index[iobs] != -1){
		if (iobs_hom == obs_hom_index[iobs]){
		  oaspl_hom[nZNA*iblock + iZNA][iobs_hom] += oaspl[nZNA*iblock + iZNA][iobs];
		  bin_size[iobs_hom] += 1;
		  for (int ifreq=0; ifreq<num_freq; ++ifreq){ 
		    psd_hom[nZNA*iblock + iZNA][ifreq + num_freq*iobs_hom] += psd[nZNA*iblock + iZNA][ifreq + num_freq*iobs];
		  }
		}
		else{
		  iobs_hom = obs_hom_index[iobs];
		  oaspl_hom[nZNA*iblock + iZNA][iobs_hom] += oaspl[nZNA*iblock + iZNA][iobs];
		  bin_size[iobs_hom] += 1;
		  for (int ifreq=0; ifreq<num_freq; ++ifreq){ 
		    psd_hom[nZNA*iblock + iZNA][ifreq + num_freq*iobs_hom] += psd[nZNA*iblock + iZNA][ifreq + num_freq*iobs];
		  }
		}
	      }
	    }
	  
	    //normalize	
	    for (int iobs_hom=0; iobs_hom<num_obs_hom; ++iobs_hom){	  
	    
	      oaspl_hom[nZNA*iblock + iZNA][iobs_hom] /= bin_size[iobs_hom];
	      
	      for (int ifreq=0; ifreq<num_freq; ++ifreq){ 
		psd_hom[nZNA*iblock + iZNA][ifreq + num_freq*iobs_hom] /=  bin_size[iobs_hom];
	      }
	    }
	  }
	}
      
	// clean up
	delete[] bin_size;
      }
    }
  }
  

  // ----------------------------
  // Zonal Noise Analysis (ZNA)
  void calculate_spectra_in_dB(){
    if (mpi_rank == 0) {
      for (int iblock=0; iblock<nblock_t; ++iblock){
        for (int iZNA=0; iZNA<nZNA; ++iZNA){
          // all observer(s)
          for (int iobs=0; iobs<num_obs; ++iobs) {
            for (int ifreq=0; ifreq<num_freq; ++ifreq) {
              psd[nZNA*iblock + iZNA][num_freq*iobs+ifreq] = 10.0*log10(psd[nZNA*iblock + iZNA][num_freq*iobs+ifreq]/p2_rel);
            }
            oaspl[nZNA*iblock + iZNA][iobs]=10.0*log10(oaspl[nZNA*iblock + iZNA][iobs]/p2_rel);
          }
          
          // homogeneous observer(s) (if any)
          for (int iobs_hom=0; iobs_hom<num_obs_hom; ++iobs_hom){
            for (int ifreq=0; ifreq<num_freq; ++ifreq) {
              psd_hom[nZNA*iblock + iZNA][num_freq*iobs_hom+ifreq] = 10.0*log10(psd_hom[nZNA*iblock + iZNA][num_freq*iobs_hom+ifreq]/p2_rel);
            }
            oaspl_hom[nZNA*iblock + iZNA][iobs_hom]=10.0*log10(oaspl_hom[nZNA*iblock + iZNA][iobs_hom]/p2_rel);
          }
        }
      }
    }
  }
  
  
  // ----------------------------
  // Zonal Noise Analysis (ZNA) - time calculation done for total surface only
  void calculate_pressure_in_time(double window_coeff){
    if (!time_calc) time_calc = true;
    
    if (mpi_rank == 0) {

      // Issue warning if hanning window used ...
      if (window_coeff != 2.0){  // with hanning window
	cout << " > !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! " << endl;
	cout << " > !!!                   WARNING                    !!! " << endl;
	cout << " > !!!                                              !!! " << endl;
	cout << " > !!! Periodic Hanning window should NOT be used   !!! " << endl;
	cout << " > !!! with OPTION P_TIME for temporal computations !!! " << endl;
	cout << " > !!!                                              !!! " << endl;
	cout << " > !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! " << endl;
      }
         
      p = new double[nt*num_obs];
      
      complex<double> (*input) = new complex<double>[num_freq_total];
      double (*output) = new double[nt];
      
      fftw_plan plan;
      
      plan = fftw_plan_dft_c2r_1d(nt, reinterpret_cast<fftw_complex*>(input), output, FFTW_MEASURE);
      // ********************************************************************
      // correction for Hanning window:
      // 2 for amplitude p_hat
      // 8/3 for energy p_hat.p_hat*
      // Recall: window_coeff = sqrt(8/3) already applied to p_hat
      //
      // correction for no window:
      // Recall: window_coeff = 2.0 set in initParameters()
      // ********************************************************************
           
      // To take into account the window correction for amplitude, uncomment the next line:
       double coeff = 2.0/window_coeff;
      // Otherwise:
      // double coeff = 1.0;
       
      for (int iobs=0 ; iobs<num_obs; iobs++){	
	// make the input
	for (int ifreq=0; ifreq<num_freq; ifreq++){
	  input[ifreq] = coeff*p_hat[0][ifreq + num_freq*iobs];
	}
	for (int ifreq=num_freq; ifreq<num_freq_total; ifreq++){
	  input[ifreq] = 0.0;
	}
	fftw_execute(plan);
	
	// replace from output
	for (int it=0; it<nt; it++){
	  p[it + nt*iobs] = output[it];
	}
      }
      
      // clean up
      fftw_destroy_plan(plan);
      delete[] input;
      delete[] output;
    }
  }
  

  // ----------------------------
  // Zonal Noise Analysis (ZNA)
  void write_p_hat(double window_coeff) {    
    if (mpi_rank == 0) {
      char real_s[128];
      char imag_s[128];

      // ********************************************************************
      // correction for Hanning window:
      // 2 for amplitude p_hat
      // 8/3 for energy p_hat.p_hat*
      // Recall: window_coeff = sqrt(8/3) already applied in phat
      //
      // correction for no window:
      // Recall: window_coeff = 2.0 already set
      // ********************************************************************
      
      // To take into account the window correction for amplitude, uncomment the next line:
      double coeff = 2.0/window_coeff;
      // Otherwise:
      // double coeff = 1.0;
      for (int iblock=0; iblock<nblock; ++iblock){
	for (int iZNA=0; iZNA<nZNA; ++iZNA){
	  for (int iobs=0; iobs<num_obs; ++iobs){ 
	    for (int ifreq=0; ifreq<num_freq; ++ifreq) {
	      sprintf(real_s,"%20.10e",coeff*p_hat[nZNA*iblock + iZNA][num_freq*iobs+ifreq].real()); 
	      sprintf(imag_s,"%20.10e",coeff*p_hat[nZNA*iblock + iZNA][num_freq*iobs+ifreq].imag()); 
	      
	      cout << " > block " << iblock << ", zone "<< iZNA <<", iobs = " << iobs << ", ifreq = " << ifreq << ", p_hat = " << real_s << " " << imag_s << endl;
	    }
	  }
	}
      }
    }
  }



  // ----------------------------
  // Zonal Noise Analysis (ZNA)
  void write_p_hat_file(double window_coeff) {   
    // rank 0 writes the file...
    if (mpi_rank == 0) {
      char filename[128];
      char extensionname[128];
      int j=0;
      //int k=0;
      int istart = array_start[this_array];
      int iend = array_end[this_array];
      string name = array_name[this_array];
             
      double (*alpha) = new double[iend-istart];
      double (*theta) = new double[iend-istart];

      // *******************************************
      // Update on angles alpha & beta calculation
      double xc[3];
      xc[0] = array_center1[this_array];
      xc[1] = array_center2[this_array];
      xc[2] = array_center3[this_array];

      for (int iZNA=0; iZNA< nZNA; ++iZNA){
	// Check if the array has homogeneous observer and if ZNA requested
	if (obs_hom_index[istart] == -1){
	  if (iZNA > 0){
	    sprintf(extensionname,"%s%02d%s","_zone",iZNA,".dat");
	  }
	  else{
	    sprintf(extensionname,".dat");
	  }
	}
	else{
	  if (iZNA > 0){
	    sprintf(extensionname,"%s%02d%s","_zone",iZNA,"_theta.dat");
	  }
	  else{
	    sprintf(extensionname,"_theta.dat");
	  }
	}
      
	// save p_hat ----------------------------------------------
	sprintf(filename,"%s_%s%s",name.c_str(),"p_hat",extensionname);
	cout << " > writing p_hat in frequency domain: " << filename << endl;    
	FILE * fp=fopen(filename,"w");
	assert(fp != NULL);
	
	// writing definition of zone as a reminder
	if (iZNA > 0){
	  fprintf(fp,"# Zonal Noise Analysis for zone %i: direction %i (0 for x, 1 for y, 2 for z); limits: min = %6.2f,  max = %6.2f\n", iZNA, ZNAcoord, ZNAlimit[iZNA][0], ZNAlimit[iZNA][1]);
	}

	// writing FFT block as a reminder
	if (nblock_t > 1){
	  fprintf(fp,"# FFT overlap: total number of fft blocks %i, each with %i surface data files.\n", nblock, nt);
	}
	
	// writing observer position and angles as a reminder
	compute_alpha_theta(alpha, theta,istart, iend, xc);
	for (int iobs=istart; iobs<iend; ++iobs){
	  fprintf(fp,"# Observer %i, x=%20.10e  y=%20.10e  z=%20.10e  alpha=%6.2f  theta=%6.2f\n", iobs+1, x_obs[3*iobs + 0], x_obs[3*iobs + 1], x_obs[3*iobs + 2], alpha[j], theta[j]);
	  j++;
	}

	// write tecplot header
	fprintf(fp,"VARIABLES = \"Frequency\"\n");
	for (int iblock=0; iblock<nblock; ++iblock){
	  j=0;
	  for (int iobs=istart; iobs<iend; ++iobs){
	    fprintf(fp,"\"p_hat%i_block%i_real\"\n",j+1,iblock+1);
	    fprintf(fp,"\"p_hat%i_block%i_imag\"\n",j+1,iblock+1);
	    j++;
	  }
	}
	fprintf(fp,"ZONE T=\"%s\"\n",name.c_str());
	
	// write frequency and p_hat
	// ********************************************************************
	// correction for Hanning window:
	// 2 for amplitude p_hat
	// 8/3 for energy p_hat.p_hat*
	// Recall: window_coeff = sqrt(8/3) already applied in void initSources()
	//
	// correction for no window:
	// Recall: window_coeff = 2.0 set in void initSources()     
	// ********************************************************************
	
	// To take into account the window correction for amplitude, uncomment the next line:
	double coeff = 2.0/window_coeff;
	// Otherwise:
	// double coeff = 1.0;
	
	for (int ifreq=0; ifreq<num_freq; ++ifreq){
	  fprintf(fp,"%20.10e ",f[ifreq]);
	  for (int iblock=0; iblock<nblock; ++iblock){
	    for (int iobs=istart; iobs<iend; ++iobs){
	      fprintf(fp,"%20.10e %20.10e ", coeff*p_hat[nZNA*iblock + iZNA][num_freq*iobs+ifreq].real(), coeff*p_hat[nZNA*iblock + iZNA][num_freq*iobs+ifreq].imag());
	    }
	  }
	  fprintf(fp,"\n");
	} 
          
	fclose(fp);

      }
      
      // clean up
      delete[] alpha;
      delete[] theta;
    }
  }



  // ----------------------------
  // Zonal Noise Analysis (ZNA)
  void write_psd_and_oaspl_file() {
    if (mpi_rank == 0) {
      char filename[128];
      char extensionname[128];
      int j = 0;
      int istart = array_start[this_array];
      int iend = array_end[this_array];
      string name = array_name[this_array];

      double (*alpha) = new double[iend-istart];
      double (*theta) = new double[iend-istart];

      // *******************************************
      // Update on angles alpha & beta calculation
      double xc[3];
      xc[0] = array_center1[this_array];
      xc[1] = array_center2[this_array];
      xc[2] = array_center3[this_array];
      
      cout << "WRITE array " << name << ", observer(s) " << istart +1 << " to " << iend <<endl;

      for (int iZNA=0; iZNA< nZNA; ++iZNA){
	// ======================================================================
	//
	// First, save all the observer(s) in the array
	//
	// ======================================================================
	// Check if the array has homogeneous observer and if ZNA requested
	if (obs_hom_index[istart] == -1){
	  if (iZNA > 0){
	    sprintf(extensionname,"%s%02d%s","_zone",iZNA,".dat");
	  }
	  else{
	    sprintf(extensionname,".dat");
	  }
	}
	else{
	  if (iZNA > 0){
	    sprintf(extensionname,"%s%02d%s","_zone",iZNA,"_theta.dat");
	  }
	  else{
	    sprintf(extensionname,"_theta.dat");
	  }
	}
	
	// save PSD ----------------------------------------------
	sprintf(filename,"%s_%s%s",name.c_str(),"PSD",extensionname);
	cout << " > writing PSD file: " << filename << endl;       
	FILE * fp=fopen(filename,"w");
	assert(fp != NULL);

	// writing definition of zone as a reminder
	if (iZNA > 0){
	  fprintf(fp,"# Zonal Noise Analysis for zone %i: direction %i (0 for x, 1 for y, 2 for z); limits: min = %6.2f,  max = %6.2f\n", iZNA, ZNAcoord, ZNAlimit[iZNA][0], ZNAlimit[iZNA][1]);
	}

	// writing FFT block as a reminder
	if (nblock_t > 1){
	  fprintf(fp,"# FFT overlap: total number of fft blocks %i, each with %i surface data files.\n", nblock, nt);
	}
	
	// writing observer position and angles as reminder ... to be improved in future releases
	compute_alpha_theta(alpha, theta, istart, iend, xc);
	for (int iobs=istart; iobs<iend; ++iobs){
	  fprintf(fp,"# Observer %i, x=%12.6f  y=%12.6f  z=%12.6f  alpha=%6.2f  theta=%6.2f\n", iobs+1, x_obs[3*iobs + 0], x_obs[3*iobs + 1], x_obs[3*iobs + 2], alpha[j], theta[j]);
	  j++;
	}
		
	// write tecplot header
	j=0;
	fprintf(fp,"VARIABLES = \"Frequency\"\n");
	for (int iobs=istart; iobs<iend; ++iobs){
	  fprintf(fp,"\"PSD%i\"\n",j+1);
	  j++;
	}
	if (nblock_t > 1 ) {
	  for (int iblock=1; iblock<nblock_t; ++iblock){
	    j=0;
	    for (int iobs=istart; iobs<iend; ++iobs){
	      fprintf(fp,"\"PSD%i_block%i\"\n",j+1,iblock);
	      j++;
	    }
	  }
	}
	fprintf(fp,"ZONE T=\"%s\"\n",name.c_str());
	
	// write frequency and PSD
	for (int ifreq=0; ifreq<num_freq; ++ifreq){
	  fprintf(fp,"%20.10e ",f[ifreq]);
	  for (int iblock=0; iblock<nblock_t; ++iblock){
	    for (int iobs=istart; iobs<iend; ++iobs){
	      fprintf(fp,"%20.10e ", psd[nZNA*iblock + iZNA][num_freq*iobs+ifreq]);
	    }
	  }
	  fprintf(fp,"\n");
	}
	
	fclose(fp);
	
	
	// save OASPL ----------------------------------------------
	sprintf(filename,"%s_%s%s",name.c_str(),"OASPL",extensionname);
	cout << " > writing OASPL file: " << filename <<endl;
	
	fp=fopen(filename,"w");
	assert(fp != NULL);

	// writing definition of zone as a reminder
	if (iZNA > 0){
	  fprintf(fp,"# Zonal Noise Analysis for zone %i: direction %i (0 for x, 1 for y, 2 for z); limits: min = %6.2f,  max = %6.2f\n", iZNA, ZNAcoord, ZNAlimit[iZNA][0], ZNAlimit[iZNA][1]);
	}
	
	// writing FFT block as a reminder
	if (nblock_t > 1){
	  fprintf(fp,"# FFT overlap: total number of fft blocks %i, each with %i surface data files.\n", nblock, nt);
	}

	// write tecplot header and frequency range as a reminder ...
	fprintf(fp,"# OASPL computed from frequency %12.6f to %12.6f \n", oaspl_f_min, oaspl_f_max);
	fprintf(fp,"VARIABLES = \"X\" \"Y\" \"Z\" \"ALPHA\" \"THETA\" \"OASPL\" \n");
	if (nblock_t > 1 ) {
	  for (int iblock=1; iblock<nblock_t; ++iblock){
	    fprintf(fp,"\"OASPL_block%i\"\n",iblock);
	  }
	}
	fprintf(fp,"ZONE T=\"%s\"\n",name.c_str());

	j=0;
	// writing observer position and OASPL
	for (int iobs=istart; iobs<iend; ++iobs){
	  fprintf(fp,"%20.10e %20.10e %20.10e %6.2f %6.2f", x_obs[3*iobs + 0], x_obs[3*iobs + 1], x_obs[3*iobs + 2], alpha[j], theta[j]);
	  for (int iblock=0; iblock<nblock_t; ++iblock){
	    fprintf(fp,"%20.10e ", oaspl[nZNA*iblock + iZNA][iobs]);
	  }
	  fprintf(fp,"\n");
	  j++;
	}  
	
	fclose(fp);	
      }
      

      // ======================================================================
      //
      // Now, save the homogeneous observer(s) (if any) associated with the array
      //
      // ======================================================================
      if (obs_hom_index[istart] != -1){
	  
	int istart_hom = array_hom_start[this_array];
	int iend_hom = array_hom_end[this_array];
	cout << "WRITE array " << name << ", homogeneous observer(s) " << istart_hom +1 << " to " << iend_hom << endl;
	
	for (int iZNA=0; iZNA<nZNA; ++iZNA){
	  // Check if ZNA requested
	  if (iZNA > 0){
	    sprintf(extensionname,"%s%02d%s","_zone",iZNA,".dat");
	  }
	  else{
	    sprintf(extensionname,".dat");
	  }
	  
	  // save PSD ----------------------------------------------
	  sprintf(filename,"%s_%s%s",name.c_str(),"PSD",extensionname);
	  cout << " > writing PSD file: " << filename << " for homogeneous observer(s) of the array." <<  endl;   
	  FILE * fp=fopen(filename,"w");
	  assert(fp != NULL);
	  
	  // writing definition of zone as a reminder
	  if (iZNA > 0){
	    fprintf(fp,"# Zonal Noise Analysis for zone %i: direction %i (0 for x, 1 for y, 2 for z); limits: min = %6.2f,  max = %6.2f\n", iZNA, ZNAcoord, ZNAlimit[iZNA][0], ZNAlimit[iZNA][1]);
	  }
	  
	  j=0;
	  // writing observer position as reminder
	  for (int iobs_hom=istart_hom; iobs_hom<iend_hom; ++iobs_hom){ 
	    fprintf(fp,"# Homogeneous Observer %i, x=%12.6f  y=%12.6f  z=%12.6f  alpha=%6.2f  theta=%6.2f \n", iobs_hom+1, x_obs_hom[3*iobs_hom + 0], x_obs_hom[3*iobs_hom + 1], x_obs_hom[3*iobs_hom + 2], alpha[j], theta[j]);
	    j=j+(iend-istart)/(iend_hom-istart_hom);
	  }
	  
	  //write tecplot header
	  j=0;
	  fprintf(fp,"VARIABLES = \"Frequency\"\n");
	  for (int iobs_hom=istart_hom; iobs_hom<iend_hom; ++iobs_hom){
	    fprintf(fp,"\"PSD%i\"\n",j+1);
	    j++;
	  }
	  if (nblock_t > 1 ) {
	    for (int iblock=1; iblock<nblock_t; ++iblock){
	      j=0;
	      for (int iobs_hom=istart_hom; iobs_hom<iend_hom; ++iobs_hom){
		fprintf(fp,"\"PSD%i_block%i\"\n",j+1,iblock);
		j++;
	      }
	    }
	  }
	  fprintf(fp,"ZONE T=\"%s\"\n",name.c_str());
	  
	  // write frequency and PSD
	  for (int ifreq=0; ifreq<num_freq; ++ifreq){
	    fprintf(fp,"%20.10e ",f[ifreq]);
	    for (int iblock=0; iblock<nblock_t; ++iblock){
	      for (int iobs_hom=istart_hom; iobs_hom<iend_hom; ++iobs_hom){
		fprintf(fp,"%20.10e ", psd_hom[nZNA*iblock + iZNA][num_freq*iobs_hom+ifreq]);
	      }
	    }
	    fprintf(fp,"\n");
	  }
	  
	  fclose(fp);
	  
	  // save OASPL ----------------------------------------------
	  sprintf(filename,"%s_%s%s",name.c_str(),"OASPL",extensionname);
	  cout << " > writing OASPL file: " << filename << " for homogeneous observer(s) of the array." <<  endl;

	  fp=fopen(filename,"w");
	  assert(fp != NULL);
	  
	  // writing definition of zone as a reminder
	  if (iZNA > 0){
	    fprintf(fp,"# Zonal Noise Analysis for zone %i: direction %i (0 for x, 1 for y, 2 for z); limits: min = %6.2f,  max = %6.2f\n", iZNA, ZNAcoord, ZNAlimit[iZNA][0], ZNAlimit[iZNA][1]);
	  }

	  // write tecplot header and frequency range as a reminder ...
	  fprintf(fp,"# OASPL computed from frequency %12.6f to %12.6f \n", oaspl_f_min, oaspl_f_max);
	  fprintf(fp,"VARIABLES = \"X\" \"Y\" \"Z\" \"ALPHA\" \"THETA\" \"OASPL\"\n");
	  if (nblock_t > 1 ) {
	    for (int iblock=1; iblock<nblock_t; ++iblock){
	      fprintf(fp,"\"OASPL_block%i\"\n",iblock);
	    }
	  }
	  fprintf(fp,"ZONE T=\"%s\"\n",name.c_str());
	  
	  j=0;
	  // writing observer position and OASPL
	  for (int iobs_hom=istart_hom; iobs_hom<iend_hom; ++iobs_hom){
	    fprintf(fp,"%20.10e %20.10e %20.10e %6.2f %6.2f", x_obs_hom[3*iobs_hom + 0], x_obs_hom[3*iobs_hom + 1], x_obs_hom[3*iobs_hom + 2], alpha[j], theta[j]);
	    for (int iblock=0; iblock<nblock_t; ++iblock){
	      fprintf(fp,"%20.10e ", oaspl_hom[nZNA*iblock + iZNA][iobs_hom]);
	    }
	    fprintf(fp,"\n");
	    j=j+(iend-istart)/(iend_hom-istart_hom);
	  }
	  
	  fclose(fp);
	  
	}
      }

      // clean up
      delete[] alpha;
      delete[] theta;
    }
  }

  // ----------------------------
  // Zonal Noise Analysis (ZNA) - remove time domain option

  /* ==========================================================
     
  Write a file containing the acoustic pressure in the 
  time domain for each observer

  Not efficient, would need time-domain implementation of FW-H solver
 
  */

  
  
  void write_p_files_in_time(int start, int it_first, int it_last, int delta, double dt, double window_coeff) { 
    if (!time_calc) calculate_pressure_in_time(window_coeff);
    
    if (mpi_rank == 0) {
      char filename[128];
      char extensionname[128];
      int j=0;
      string name = array_name[this_array];
      int istart = array_start[this_array];
      int iend = array_end[this_array];

      //cout << name << " > it_first = " <<  it_first << ", it_last = " << it_last << ", nrt = " << it_last - it_first +1 << endl;

      double (*alpha) = new double[iend-istart];
      double (*theta) = new double[iend-istart];

      // *******************************************
      // Update on angles alpha & beta calculation
      double xc[3];
      xc[0] = array_center1[this_array];
      xc[1] = array_center2[this_array];
      xc[2] = array_center3[this_array];

      // Check if the array has homogeneous observer
      if (obs_hom_index[istart] == -1){
	sprintf(extensionname,".dat");
      }
      else{
	sprintf(extensionname,"_theta.dat");
      }
      
      // p_time ----------------------------------------------
      sprintf(filename,"%s_%s%s",name.c_str(),"p_time",extensionname);
      cout << " > writing p in time domain: " <<filename<< endl;
      FILE * fp=fopen(filename,"w");
      assert(fp != NULL);
      
      // writing observer position as reminder
      compute_alpha_theta(alpha, theta, istart, iend, xc);
      for (int iobs=istart; iobs<iend; iobs++){
	fprintf(fp,"# Observer %i, x=%20.10e  y=%20.10e  z=%20.10e  alpha=%6.2f  theta=%6.2f\n", iobs+1, x_obs[3*iobs + 0], x_obs[3*iobs + 1], x_obs[3*iobs + 2], alpha[j], theta[j]);
	j++;
      }

      // write tecplot header
      j=0;
      fprintf(fp,"VARIABLES = \"t\"\n");
      for (int iobs=istart; iobs<iend; iobs++){
	fprintf(fp,"\"p_time%i\"\n",j+1);
	j++;
      }
      fprintf(fp,"ZONE T=\"%s\"\n",name.c_str());
      
      // ****************************
      // recover full time history with retarded time calculation
      // write time and p_time
      int irt = it_first;
      int nrt = it_last - it_first + 1;
      for (int it=0 ; it<nrt; it++) { 
	double time = (start + (it_first+it) * delta)*dt;
	fprintf(fp,"%20.10e ",time);	
	for (int iobs=istart; iobs<iend; iobs++){
          fprintf(fp,"%20.10e ",p[irt + nt*iobs]);
	}
	fprintf(fp,"\n");
	irt +=1;
	if (irt > nt) {
	  irt = 0;
	}
      }

      fclose(fp);

      // clean up
      delete[] alpha;
      delete[] theta;
    }
  }
  


  // Compute angles assuming center of arc at xc (default: xc=(0,0,0)) ... to be improved in future releases
  void compute_alpha_theta(double *alpha, double *theta, const int istart, const int iend){
    const double xc[3]={0.0,0.0,0.0};
    compute_alpha_theta(alpha, theta, istart, iend, xc);
  }


  // Compute angles assuming center of arc at xc  ... to be improved in future releases
  void compute_alpha_theta(double *alpha, double *theta, const int istart, const int iend, const double xc[3]){
    double alpha_offset;
    double theta_offset;
    double zovery;
    int j =0;
    double x,y,z;

    for (int iobs=istart; iobs<iend; ++iobs){
      if (axis == 3) {
	x=x_obs[3*iobs + 2]-xc[2];
	y=x_obs[3*iobs + 1]-xc[1];
        z=x_obs[3*iobs + 0]-xc[0];
      }
      else {
	x=x_obs[3*iobs + 0]-xc[0];
	y=x_obs[3*iobs + 1]-xc[1];
	z=x_obs[3*iobs + 2]-xc[2];
      }

      if (z <= 0.0){
	if (y >= 0.0){
	  theta_offset = 0.0;
	}
	else {
	  theta_offset = 180.0;
	}
      }
      else { 
	if (y <= 0.0){
	  theta_offset = 180.0;
	}
	else {
	  theta_offset = 360.0; 
	}
      }

      if ((z==0) && (y==0)){
	zovery = 0.0;
      }
      else{
	zovery = z/y;
      }
	

      if (x < 0.0){
	alpha_offset = 0.0;
      }
      else{
	alpha_offset = 180.0;
      }
      
      
	   
      theta[j] = theta_offset - atan(zovery)/M_PI*180.0;

      alpha[j] = alpha_offset - atan(sqrt(y*y + z*z)/x)/M_PI*180.0;
      
      j++;
    }
  }




  /* ==========================================================
     
  Write a binary plot3D geometry file
  
  */
  void write_PLOT3D_geom_file() {   
    // rank 0 writes the file...
    if (mpi_rank == 0) {
      char filename[128];
      string name = array_name[this_array];
      int m1 = array_m1[this_array];
      int m2 = array_m2[this_array];
      int m3 = array_m3[this_array];
      int i0 = array_start[this_array];
      sprintf(filename,"%s_%s.grid",name.c_str(),"PLOT3D");
      cout << " > writing PLOT_3D geometry file: " << filename << endl;
      
      FILE * fp=fopen(filename,"w");
      assert(fp != NULL);
	  
      fwrite(&m1, sizeof(int),1,fp);
      fwrite(&m2, sizeof(int),1,fp);
      fwrite(&m3, sizeof(int),1,fp);
	  
      FOR_I3 
	for (int i3=0 ; i3<m3; ++i3){
	  for (int i2=0 ; i2<m2; ++i2){
	    for (int i1=0 ; i1<m1; ++i1)
	      {
		int iobs = obs_index[i0 + i1 + i2 * m1 + i3 * m2 * m1];
		double coord = x_obs[3*iobs + i];
		fwrite(&coord, sizeof(double),1,fp);
	      }
	  }
	}

      fclose(fp);

    }
  }

  

  // ----------------------------
  // Zonal Noise Analysis (ZNA) - remove PLOT3D option

  /* ==========================================================
     
  Write a binary plot3D solution file containing OASPL
  
  */
 
  
  void write_PLOT3D_OASPL_file(){
    if (mpi_rank == 0) {
      char filename[128];
      string name = array_name[this_array];
      sprintf(filename,"%s_%s.dat",name.c_str(),"OASPL_PLOT3D");
      cout << " > writing PLOT_3D OASPL file: " << filename << endl;
      write_PLOT3D_scalar_file(oaspl[0],filename);
    }
  }
  

  /* ==========================================================
     
  Write a time sequence of binary plot3D solution file containing 
  the acoustic pressure in the time domain 
  
  */

  
  void write_PLOT3D_p_files_in_time(const char * const prefix,int start, int it_first, int it_last, int delta, double window_coeff ) { 

    if (!time_calc) calculate_pressure_in_time(window_coeff);
    
    if (mpi_rank == 0) {
      char filename[128];
      string name = array_name[this_array];
      //int pt_size = array_end[this_array] - array_start[this_array];
      double *pt = new double[num_obs];
      
      cout << " > writing PLOT3D p_time files for time steps = " <<start + it_first*delta<< " to " << start + it_last* delta<<", every " << delta<< endl;
      
      for (int it=it_first ; it<=it_last; it++) {
	int step = start + it * delta;
	sprintf(filename,"%s%s_%s_%06d.dat",name.c_str(),prefix,"P_PLOT3D",step);
	
	for (int iobs=array_start[this_array]; iobs<array_end[this_array]; iobs++){
          pt[iobs] = p[it + nt*iobs];
	}
	if ((it+1)%10 == 0){
	  cout << " > writing PLOT_3D pressure file in time domain: " << filename << endl;
	}
	write_PLOT3D_scalar_file(pt,filename);
      }
         
      // clean up
      delete[] pt;
    }
  }
  

  
  void write_PLOT3D_scalar_file(double *input, char *filename){   
    FILE * fp=fopen(filename,"w");
    assert(fp != NULL);

    int m1 = array_m1[this_array];
    int m2 = array_m2[this_array];
    int m3 = array_m3[this_array];
    int i0 = array_start[this_array];
    
    fwrite(&m1, sizeof(int),1,fp);
    fwrite(&m2, sizeof(int),1,fp);
    fwrite(&m3, sizeof(int),1,fp);
    
    int dummy = 1;
    fwrite(&dummy, sizeof(int),1,fp);
    
    for (int i3=0 ; i3<m3; i3++){
      for (int i2=0 ; i2<m2; i2++){
	for (int i1=0 ; i1<m1; i1++){
	  int iobs = obs_index[i0 + i1 + i2 * m1 + i3 * m2 * m1];
	  fwrite(&(input[iobs]), sizeof(double),1,fp);
	}
      }
    }

    fclose(fp);
  }
  

  void writeObserver(Param * param) {
    bool name_flag = false;
    bool geom_flag = false;
    bool option_flag = false;
    bool plot3D_flag = false;
    bool p_time_flag = false;
    bool p_hat_flag = false;
    bool plot3D_locked = false;
    string name;
  
    int i = 0;
    while (i < param->size()) {
      string token = param->getString(i++);
      if (token == "NAME") {
	assert( !name_flag );
	name_flag = true;
	name = param->getString(i++);
      }
      else if (token == "GEOM") {
	assert( !geom_flag);
	geom_flag = true;
	token = param->getString(i++);
	if (token == "POINT") { 
	  double x[3];
	  x[0] = param->getDouble(i++);
	  x[1] = param->getDouble(i++);
	  x[2] = param->getDouble(i++);
	  // *******************************************
	  // Update on angles alpha & beta calculation
	  double xc[3] = {0.0, 0.0,0.0};
	  if (i < param->size()){
	    token = param->getString(i++);
	    if (token == "CENTER") {
	      xc[0] = param->getDouble(i++);
	      xc[1] = param->getDouble(i++);
	      xc[2] = param->getDouble(i++);
	    }
	    else{
	      i--;
	    }
	  }
	  addPoint(x,name,xc);
	  plot3D_locked = true; // PLOT3D not a valid option
	}
	else if (token == "ARC") {
	  double xa = param->getDouble(i++);
	  double ra = param->getDouble(i++);
	  int ntheta = param->getInt(i++);
	  // *******************************************
	  // Update on angles alpha & beta calculation
	  double xc[3] = {0.0, 0.0,0.0};
	  if (i < param->size()){
	    token = param->getString(i++);
	    if (token == "CENTER") {
	      xc[0] = param->getDouble(i++);
	      xc[1] = param->getDouble(i++);
	      xc[2] = param->getDouble(i++);
	    }
	    else{
	      i--;
	    }
	  }
	  // *******************************************
	  // store theta offset (namely start of theta array) in deg 
	  double * theta_offset = new double[ntheta];
	  for (int itheta=0; itheta<ntheta ; itheta++) theta_offset[itheta] = 0.0;
	  if (i < param->size()){
	    token = param->getString(i++);
	    if (token == "THETA_OFFSET") {
	      for (int itheta=0; itheta<ntheta ; itheta++) theta_offset[itheta] = param->getDouble(i++);
	    }
	    else{
	      i--;
	    }
	  }
	  addPointArc(xa, ra, ntheta, name, xc, theta_offset);
	  plot3D_locked = true; // PLOT3D not a valid option
	  delete[] theta_offset; 
	}
	else if ((token == "ROTATINGARC")||(token == "ROTATING_ARC")) {
	  double x_orig[3]; 
	  x_orig[0] = param->getDouble(i++);
	  x_orig[1] = param->getDouble(i++);
	  x_orig[2] = param->getDouble(i++);
	  double r = param->getDouble(i++);
	  double alpha_min = param->getDouble(i++);
	  double alpha_max = param->getDouble(i++);
	  int nalpha = param->getInt(i++);
	  int ntheta = param->getInt(i++);
	  // *******************************************
	  // Update on angles alpha & beta calculation
	  double xc[3] = {0.0, 0.0,0.0};
	  if (i < param->size()){
	    token = param->getString(i++);
	    if (token == "CENTER") {
	      xc[0] = param->getDouble(i++);
	      xc[1] = param->getDouble(i++);
	      xc[2] = param->getDouble(i++);
	    }
	    else{
	      i--;
	    }
	  }
	  // *******************************************
	  // store theta offset (namely start of theta array) in deg 
	  double * theta_offset = new double[ntheta];
	  for (int itheta=0; itheta<ntheta ; itheta++) theta_offset[itheta] = 0.0;
	  if (i < param->size()){
	    token = param->getString(i++);
	    if (token == "THETA_OFFSET") {
	      for (int itheta=0; itheta<ntheta ; itheta++) theta_offset[itheta] = param->getDouble(i++);
	    }
	    else{
	      i--;
	    }
	  }
	  addRotatingArc(x_orig, r, alpha_min, alpha_max, nalpha, ntheta, name, xc, theta_offset);
	  delete[] theta_offset; 
	}
	else if ((token == "CARTMESH")||(token == "CART_MESH")){
	  double x0[3];
	  x0[0] = param->getDouble(i++);
	  x0[1] = param->getDouble(i++);
	  x0[2] = param->getDouble(i++);
	  double x1[3];
	  x1[0] = param->getDouble(i++);
	  x1[1] = param->getDouble(i++);
	  x1[2] = param->getDouble(i++);
	  double x2[3];
	  x2[0] = param->getDouble(i++);
	  x2[1] = param->getDouble(i++);
	  x2[2] = param->getDouble(i++);
	  int n01 = param->getInt(i++);
	  int n02 = param->getInt(i++);
	  // *******************************************
	  // Update on angles alpha & beta calculation
	  double xc[3] = {0.0, 0.0,0.0};
	  if (i < param->size()){
	    token = param->getString(i++);
	    if (token == "CENTER") {
	      xc[0] = param->getDouble(i++);
	      xc[1] = param->getDouble(i++);
	      xc[2] = param->getDouble(i++);
	    }
	    else{
	      i--;
	    }
	  }
	  addCartMesh(x0,x1,x2,n01,n02,name,xc);
	}	
	else {
	  if (mpi_rank == 0)
	    cerr << "\n\n\n*********************************************************\n" <<
	      "Error: unrecognized entry in OBSERVER GEOM=: " << token << ". \n" <<
	      "Valid GEOM include:\n" << 
	      "POINT, ARC, ROTATING_ARC, CART MESH ...\n"
	      "*********************************************************\n\n\n" << endl;
	  throw(0);
	}
      }
      else if (token == "OPTIONS") {
	assert( !option_flag);
	option_flag = true;
	while (i < param->size()) {
	  token = param->getString(i++);
	  if (token == "PLOT3D"){
	    if (!plot3D_locked){
	      plot3D_flag = true;
	    }
	    else{
	      if (mpi_rank == 0){
		cout <<" > PLOT3D not a useful option for current microphone array "<< name <<" - OPTION DISCARDED ..." << endl;
	      }
	    }
	  }
	  else if (token == "P_HAT"){
	    p_hat_flag = true;
	  }
	  else if (token == "P_TIME"){
	    p_time_flag = true;    
	  }
	  
	
	  
	  else {
	    if (mpi_rank == 0)
	      cerr << "\n\n\n*********************************************************\n" <<
		"Error: unrecognized entry in OBSERVER OPTIONS=: " << token <<". \n" <<
		"Valid OPTIONS include:\n" << 
		"PLOT3D, P_HAT, P_TIME ...\n"
		"*********************************************************\n\n\n" << endl;
	    throw(0);
	  }
	}
      }
      else {
	if (mpi_rank == 0)
	  cerr << "\n\n\n*********************************************************\n" <<
	    "Error: unrecognized keyword in OBSERVER: " << token <<
	    ". Example of proper syntax:\n" << 
	    "OBSERVER NAME=Mic GEOM=POINT 100.0 0.0 0.0 OPTION=PLOT3D \n" << 
	    "*********************************************************\n\n\n" << endl;
	throw(0);
      }
    }
  
    // finalize options to default (false) if not defined in input
    set_array_PLOT3D(plot3D_flag);
    set_array_PLOThat(p_hat_flag);
    set_array_PLOTtime(p_time_flag);
  }



  void finalize_output(bool dB_calc, double window_coeff, int start, int delta, double dt) { 
    set_this_array(0);
    set_this_obs(0);
    
    // compute homogeneous observer(s) (if any)
    if (get_num_obs_hom() > 0){
      calcPSDandOASPL_hom();
    }
    
    // compute spectra in dB (if requested in input file)
    if (dB_calc){
      calculate_spectra_in_dB();
    }
    
    for (int i=0; i<get_num_array(); ++i){   
      // Default: always write PSD and OASPL to file
      write_psd_and_oaspl_file();
      
      if(get_array_PLOThat()){
	write_p_hat_file(window_coeff);
      }
      
      if(get_array_PLOT3D()){
  
	write_PLOT3D_geom_file();
	write_PLOT3D_OASPL_file();
	if(get_array_PLOTtime()){
	  int it_first = max(0,get_array_it_first());
	  int it_last = get_array_it_last();
	  write_PLOT3D_p_files_in_time("_time", start, it_first, it_last, delta, window_coeff);
	}
	
      }
      else {

	if(get_array_PLOTtime()){
	  int it_first = max(0,get_array_it_first()); 
	  int it_last = get_array_it_last();
	  write_p_files_in_time(start, it_first, it_last, delta, dt, window_coeff);
	}
	
      }
      
      // got to next array of observers
      inc_this_array();
    }    
  }


protected:
  
  virtual void observerHook() { 
 
  }

  

};
