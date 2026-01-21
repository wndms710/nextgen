#include "PngAnalysis.hpp"
#include "Dmd.hpp"
#include <sstream>
#include <set>

class DmdModeDescriptor { 
public: 
  double freq; 
  double damp;
  double amp_n;  // amplitude at the final instant.. 
  double amp_0;  // amplitude at the first instant
  int idx;
  
  bool operator<(const DmdModeDescriptor& _d) const { 
    return (amp_n < _d.amp_n);
  }
};//class DmdModeDescriptor

class PngDmd : public PngAnalysis {
 
public: 
  complex<double>* Vr ;
  double * lambda_i ;
  double * amp;

  PngDmd(const string& _infile,  
         const bool _bMean) : PngAnalysis(_infile, _bMean), 
                              Vr(NULL), lambda_i(NULL), amp(NULL) { }
  
  PngDmd(const string& _infile, const string &_matfile,  
         const bool _bMean) : PngAnalysis(_infile, _matfile, _bMean), 
                              Vr(NULL), lambda_i(NULL), amp(NULL) { }
  
  void commonInit() { PngAnalysis::commonInit(); }
  
  void dump() const { 

    FILE * fp          = fopen("freq_amp.dat", "w"); 
    fprintf(fp, "#Frequency   \tAmplitude_0 \tGrowth Rate \tAmplitude_n  \tMode #\n");
    //FILE * f_mode_list = fopen("mode_list.dat", "w");

    set<DmdModeDescriptor> descriptors;
    for (int i =0; i < nnz ; ++i) { 
      if ( lambda_i[i] >= 0.0) {
        // also compute and report the a_i*lambda_k^n where n is the 
        // number of time sequences...
        complex<double> lambda_c(lambda_r[i]*dt_samp,lambda_i[i]*2.0*M_PI*dt_samp);
        complex<double> mu_c     = exp(lambda_c);
        complex<double> mu_n     = pow(mu_c,n_snaps);
        const double a_n         = amp[i]*std::abs(mu_n);

        DmdModeDescriptor d; 
        d.freq     = lambda_i[i];
        d.damp     = lambda_r[i];
        d.amp_n    = a_n; 
        d.amp_0    = amp[i];
        d.idx      = i;
        descriptors.insert(d);
        //fprintf(fp, "%12.8g    %12.8g    %12.8g    %12.8g    %d\n", lambda_i[i], amp[i], lambda_r[i],a_n,i); 
      }
 
      /*
      // HACK.. find me a mode... 
      if ( (lambda_i[i] >= 0.0) && (lambda_i[i] <= 2000.0)) { 
        cout << i << "   " << lambda_i[i] << endl; 
	fprintf(f_mode_list, "%d \n", i);
      }
      */
    }//i

    // write out the frequency, amplitude information in sorted order
    for (set<DmdModeDescriptor>::reverse_iterator it = descriptors.rbegin(); it != descriptors.rend(); ++it)  
      fprintf(fp, "%12.8g\t%12.8g\t%12.8g\t%12.8g\t%d\n", it->freq, it->amp_0, it->damp, it->amp_n, it->idx); 

    fclose(fp); 
    //fclose(f_mode_list);
  }//dump()

  void writeMode(const int imo) { 
 
    FILE * fp = NULL;
    char filename[256];

    if ( mpi_rank == 0 ) {  
      cout << " > write dmd modes .... " << endl; 
      cout.flush();
      fp = fopen("mode_key.dat", "a");
    }
  
    // allocate the memory .. 
    const int nx = base->getNx(); 
    const int ny = base->getNy(); 
    const int n  = nx*ny;

    // reconstruct the mode... 
    complex<double> c_coeff; 
    double coeff;

    double my_dotp = 0.0, dotp = 0.0;
    double * my_reconstruct = new double[n]; 
    double * reconstruct    = new double[n]; 

    int * imora = NULL; calcUniformDist(imora, n_snaps-1, mpi_size); 
    for (int iter = 0; iter < 2 ; ++iter ) { 

      for (int ii =0; ii < nx*ny; ++ii) 
        my_reconstruct[ii] = 0.0;

      if (mpi_rank==0){
        if ( iter == 0 ) {  
          //sprintf(filename, "mode.%03d.real.png", imo); 
          sprintf(filename, "mode.%.2f.real.png", lambda_i[imo]);
        } else {  
          //sprintf(filename, "mode.%03d.imag.png", imo); 
          sprintf(filename,"mode.%.2f.imag.png", lambda_i[imo]);
        }

        cout << "Working on reconstructing : " << filename << "; freq : " << lambda_i[imo] <<  endl; 
    
        if ( iter == 0 ) 
          fprintf(fp, "%03d  %12.8g   ", imo, lambda_i[imo]); 
      }

      for (int i = imora[mpi_rank]; i < imora[mpi_rank+1]; ++i) { 
      //for (int i =0; i < n_snaps-1 ; ++i)  
        
        PngData * png2; 
        if ( imageCache.size() > 0 )
          png2 = getImage(i);
        else { 
          png2 = new PngData();
          png2->read(imageList[i].c_str()); 
          png2->finalizeRead(); 
        }
        
        reconstructCoeff(c_coeff, Wk, Vr, nnz, n_snaps, imo, i); 
        coeff = (iter == 0) ? c_coeff.real() : c_coeff.imag(); 

        //cout << i << "   " << coeff << endl; 
        for (int ipx = 0; ipx < n ; ++ipx) { 
          if ( base->pixel_flag[ipx] >= 0 ) { 
            my_reconstruct[ipx] += coeff * png2->pixel_data[ipx]; //coeff * png2->getPhi(ipx); 
          }
        }//ipx

        if ( imageCache.size() > 0 )
          png2->lock = false; 
        else
          delete png2;
      }//i
      
      MPI_Allreduce(my_reconstruct, reconstruct, n, MPI_DOUBLE, MPI_SUM, mpi_comm);

      for (int ipx = 0; ipx < nx*ny; ++ipx) { 
        if ( base->pixel_flag[ipx] >= 0) { 
          my_dotp += reconstruct[ipx] * reconstruct[ipx]/double(base->npx+base->npxS+base->npxP+base->npxI); 
        }
      }

      if (mpi_rank==0){

        PngData * png = new PngData();
        png->read(imageList[0].c_str());
        png->finalizeRead();

        // now set gray
        double range[2]; 
        dumpRange(range, reconstruct, nx*ny, "RECONSTRUCT");
  
        const double eps = 1.0e-08; 
        if ( (fabs(range[0]) < eps) && ( fabs(range[1]) < eps)) { 
          range[0] = -eps; 
          range[1] =  eps;
        }
  
        //cout << " Npx = " << png->npx << endl; 
        fprintf(fp, "  %12.8g   %12.8g ", range[0], range[1]); 
  
        for (int ipx = 0; ipx < nx*ny ; ++ipx ) { 
          if ( png->pixel_flag[ipx] >= 0 ) { 
            //const double val = (reconstruct[ipx] - range[0])/(range[1]-range[0]);
            //assert ( (val>=0.0)&&(val<= 1.0)) ;
            //png->setPhi(ipx,val);
            png->pixel_data[ipx] = reconstruct[ipx];
          } 
        }

        png->rescaleRangesToData();
        png->initializeWrite();
        png->write(filename); 

        delete png;

      }
    }//iter  
  
    MPI_Allreduce(&my_dotp, &dotp, 1, MPI_DOUBLE, MPI_SUM, mpi_comm);
    if (mpi_rank==0){
      fprintf(fp, "\n"); 
      fclose(fp);
      cout << " Mode normalization across " << mpi_size << " ranks: " << dotp/(float)mpi_size << endl; 
    }
   
    delete[] imora;
    delete[] reconstruct;
    delete[] my_reconstruct; 
  }//writeMode() 

  void writeModeAnimation(const int imo, const double dt) { 
    if (mpi_rank == 0)  cout << " > write dmd modes .... " << endl; 
  
    // allocate the memory .. 
    const int nx = base->getNx(); 
    const int ny = base->getNy(); 
    const int n  = nx*ny;

    // reconstruct the mode... 
    complex<double> c_coeff; 
    double coeff;

    double my_dotp = 0.0, dotp = 0.0;
    double (*my_reconstruct)[2] = new double[n][2];
    double (*reconstruct)[2]    = new double[n][2];
    
    int * imora = NULL; calcUniformDist(imora, n_snaps-1, mpi_size); 
    for (int iter = 0; iter < 2 ; ++iter ) { 

      for (int ii =0; ii < nx*ny; ++ii) 
        my_reconstruct[ii][iter] = 0.0;

      if ( mpi_rank == 0 ) {
        if ( iter == 0 ) 
          cout << " reconstructing the real part of the mode; freq : " << lambda_i[imo] << endl;
        else 
          cout << " reconstructing the imag part of the mode; freq : " << lambda_i[imo] << endl;
      }

      for (int i = imora[mpi_rank]; i < imora[mpi_rank+1]; ++i) { 
        
        PngData * png2; 
        if ( imageCache.size() > 0 )
          png2 = getImage(i);
        else { 
          png2 = new PngData();
          png2->read(imageList[i].c_str()); 
          png2->finalizeRead(); 
        }
        
        reconstructCoeff(c_coeff, Wk, Vr, nnz, n_snaps, imo, i); 
        coeff = (iter == 0) ? c_coeff.real() : c_coeff.imag(); 

        //cout << i << "   " << coeff << endl; 
        for (int ipx = 0; ipx < n ; ++ipx) { 
          if ( base->pixel_flag[ipx] >= PngData::VOLUME_DATA_PIXEL ) { 
            my_reconstruct[ipx][iter] += coeff * png2->pixel_data[ipx]; 
          }
        }//ipx

        if ( imageCache.size() > 0 )
          png2->lock = false; 
        else
          delete png2;
      }//i
      
      MPI_Allreduce((double*)(my_reconstruct), (double*)(reconstruct), 2*n, MPI_DOUBLE, MPI_SUM, mpi_comm);

      for (int ipx = 0; ipx < nx*ny; ++ipx) { 
        if ( base->pixel_flag[ipx] >= PngData::VOLUME_DATA_PIXEL) { 
          my_dotp += reconstruct[ipx][iter] * reconstruct[ipx][iter]/double(base->npx+base->npxS+base->npxP+base->npxI); 
        }
      }

    }//iter  
  
    MPI_Allreduce(&my_dotp, &dotp, 1, MPI_DOUBLE, MPI_SUM, mpi_comm);
    if ( mpi_rank==0 ) {
      cout << " Mode normalization across " << mpi_size << " ranks: " << dotp/double(mpi_size) << endl; 
    }
  
    // now reconstruct the time history... 

    if ( mpi_rank == 0 ) { 
      double grange[2] = {HUGE_VAL, -HUGE_VAL};
      
      cout << " check dt: " << dt << endl;
      
      complex<double> lambda_eig(lambda_r[imo]*dt,2.0*M_PI*lambda_i[imo]*dt);
      complex<double> mu = exp(lambda_eig);
    
      cout << " lambda; " << lambda_r[imo] << "   " << lambda_i[imo] << endl;
      cout << " mu: " << mu.real() << "   " << mu.imag() << endl;

      double range[2] = {1e16,-1e16};
      {
        cout << " re-reading an image: " << imageList[0].c_str() << endl;

        PngData * png = new PngData();
        png->read(imageList[0].c_str());
        cout << "file = " << imageList[0].c_str() << "\n";
        png->finalizeRead();

        cout << " setting range ... " << endl;

        for (int ipx = 0; ipx < nx*ny; ++ipx) { 
          if ( png->pixel_flag[ipx] >= PngData::VOLUME_DATA_PIXEL ) { 
            for (int i = 0; i < 2; ++i) { 
              range[0] = min(range[0],reconstruct[ipx][i]);
              range[1] = max(range[1],reconstruct[ipx][i]);
            }
          }
        }

        //if ( range[0] > 0.0) 
        //  range[0] /= 5.0;
        //else 
        //  range[0] *= 5.0;

        //if ( range[1] > 0.0) 
        //  range[1] *= 5.0;
        //else 
        //  range[1] /= 5.0;

        delete png;

        cout << "range: " << "   " << range[0] << "   " << range[1] << endl;
        grange[0] = range[0];
        grange[1] = range[1];

      }

      cout << " here here " << endl;

      // if iter < 2, then we will also write out the imaginary components.. 
      // this seems to be largely just the complex conjugate of the first one
      // so you could skip ahead if necessary
      
      for (int iter = 0; iter < 2; ++iter) { 
        
        for (int i =0; i < n_snaps-1; ++i) { 
  
          char filename[128];
          if ( iter == 0 ) { 
            sprintf(filename, "mode.%.2f.%06d.real.png", lambda_i[imo], i);
          } else { 
            sprintf(filename, "mode.%.2f.%06d.imag.png", lambda_i[imo], i);
          }

          cout << " writing out file: " << filename << endl;

          complex<double> mu_n = pow(mu,i);
          cout << i << "   " << mu_n.real() << "   " << mu_n.imag() << endl;

          PngData * png = new PngData();
          png->read(imageList[0].c_str(),true);
          png->finalizeRead();
          
          for (int ipx = 0; ipx < nx*ny ; ++ipx ) { 
            if ( png->pixel_flag[ipx] >= PngData::VOLUME_DATA_PIXEL ) { 
              double val;
              if (iter == 0 ) 
                val = reconstruct[ipx][0]*mu_n.real() - reconstruct[ipx][1]*mu_n.imag();
              else 
                val = reconstruct[ipx][0]*mu_n.imag() + reconstruct[ipx][1]*mu_n.real();

              png->pixel_data[ipx] = val;
            }
          }
       
          //all data types are operated on as one, so set the same range
          //for all types, initializeWrite will write the applicable ranges
          //for the data types that exist to image metadata.
          png->setRange(grange);
          png->setRangeSurface(grange);
          png->setRangeParticles(grange);
          png->setRangeIso(grange);

          png->initializeWrite();
          png->write(filename); 

          delete png;
        }// i
      }//iter
    }//mpi_rank == 0

    delete[] imora;
    delete[] reconstruct;
    delete[] my_reconstruct; 
  }//writeModeAnimation() 

  void dumpRange(double range[2], const double * buf, const int nn, const char* msg) { 

    range[0] = 1.0e+16; 
    range[1] = -1.0e+16; 

    for (int i =0; i < nn ; ++i) { 
      range[0] = fmin(range[0], buf[i]); 
      range[1] = fmax(range[1], buf[i]); 
    }

    cout << "dumpRange, " << msg << " = " << range[0] << ": " << range[1] << endl;
  }//dumpRange()

  void runDmd() {

    assert( dt_samp != 0.0) ;
    const double dt = dt_samp;
 
    if ( mpi_rank == 0 ) {
      
      assert(amp==NULL);
      assert(Wk==NULL);
      assert(lambda_r==NULL);
      assert(lambda_i==NULL);
      assert(Vr==NULL);
      assert(A!=NULL);

      nnz = dmdAnalysis(n_snaps, A, lambda_r, lambda_i, Wk, Vr, amp, dt);     

      dump();
    }
    else{
      Wk = new double[(n_snaps-1)*(n_snaps-1)];
    }


    // everyone needs the Wk, Vr, nnz
    MPI_Bcast(&nnz, 1, MPI_INT, 0, mpi_comm);     

    if (mpi_rank !=0){  
      Vr = new complex<double>[nnz*nnz];
      lambda_i = new double[nnz];
      lambda_r = new double[nnz];
    }

    // everyone needs the Wk, Vr, nnz, lambda_i for print statement..
    MPI_Bcast(Wk, (n_snaps-1)*(n_snaps-1), MPI_DOUBLE, 0, mpi_comm);
    
    // in general, the c bindings of mpi do not require that the 
    // double precision complex is supported .. therefore, we 
    // change the call into two double broadcats .. 
    
    //MPI_Bcast(Vr, nnz*nnz, MPI_DOUBLE_COMPLEX, 0, mpi_comm);
    
    { 
      double * tmp_r = new double[nnz*nnz];
      double * tmp_i = new double[nnz*nnz];

      if ( mpi_rank == 0) { 

        for (int ii = 0; ii < nnz*nnz; ++ii) { 
          tmp_r[ii] = Vr[ii].real();
          tmp_i[ii] = Vr[ii].imag();
        }

      }

      MPI_Bcast(tmp_r,nnz*nnz,MPI_DOUBLE,0,mpi_comm);
      MPI_Bcast(tmp_i,nnz*nnz,MPI_DOUBLE,0,mpi_comm);

      if ( mpi_rank != 0) { 

        for (int ii = 0; ii < nnz*nnz; ++ii) 
          Vr[ii] = complex<double>(tmp_r[ii],tmp_i[ii]);

      }

      delete[] tmp_r;
      delete[] tmp_i;
    }


    MPI_Bcast(lambda_i, nnz, MPI_DOUBLE, 0, mpi_comm);
    MPI_Bcast(lambda_r, nnz, MPI_DOUBLE, 0, mpi_comm);

    for (int i=0; i < modeList.size(); ++i) {

      // need to find the imo closest to the requested mode frequency...
      const double f_mode = modeList[i].first;

      int imo = -1;
      if ( mpi_rank == 0 ) { 

        double dist_min = HUGE_VAL;
        for (int ii =0; ii < nnz ; ++ii) { 
          const double d = abs(lambda_i[ii] - f_mode);
          if ( d < dist_min) { 
            dist_min = d;
            imo      = ii;
          }
        }
      }

      MPI_Bcast(&imo,1,MPI_INT,0,mpi_comm);
      assert( (imo >= 0) && (imo < nnz));

      writeMode(imo);
      if ( modeList[i].second) 
        writeModeAnimation(imo,dt);
    }
  }//runDmd()

  void run() { runDmd(); }

  ~PngDmd() {

    if(mpi_rank==0) cout << "~PngDmd()" << endl;

    if ( lambda_i != NULL )
      delete[] lambda_i ; 

    if ( amp != NULL)
      delete[] amp;

    if( Vr != NULL)
      delete[] Vr;
  }//~PngDmd()

  template <class T>
  bool from_string(T& t, const std::string &s, std::ios_base& (*f)(std::ios_base&))
  {
     std::istringstream iss(s);
     return !(iss >> f >> t).fail();
  }
};//class PngDmd 


int main(int argc, char* argv[]) { 
  
  try {
    initMpiEnv(argc,argv);
    {//scope

      /* Pixel based cropping (no longer supported)
      Param * param = getParam("CROP");
      if (param){
        crop_start[0] = param->getDouble(0);
        crop_start[1] = param->getDouble(1);
        crop_end[0]   = param->getDouble(2);
        crop_end[1]   = param->getDouble(3);

        if (mpi_rank==0) std::cout << "cropping from ("<< crop_start[0] <<"," << crop_start[1] <<") to (" << crop_end[0] << ","<< crop_end[1] << ")" << std::endl;
      }
      else {
        if (mpi_rank==0) std::cout << "Proceeding without a crop file" << std::endl;
      }

      for(int comp=0; comp<2; ++comp)
        if(crop_start[comp]>crop_end[comp] && crop_end[comp]!=-1) {
          if (mpi_rank==0) std::cerr << "crop start at dimension " << comp << " is larger than crop end...check the crop.in file...aborting." << std::endl;
          abort();
        }
      */

      PngDmd * dmd = NULL; 
  
      if ((argc < 3) || (argc > 5)) { 
        if(mpi_rank == 0) { 
          cerr << "INVALID NUMBER OF ARGUMENTS PROVIDED!" << endl << endl;
          cerr << " ----------------------------------------------------------------------------------------  " << endl;
          cerr << " Usage: ./imageDmd.exe <image-list-file> <dt-samp> [<cached-matrix-file> <mode-list-file>] " << endl; 
          cerr << " ----------------------------------------------------------------------------------------- " << endl;

          cerr << endl;
          cerr << endl;
          cerr << " Performing an image-based dmd requires the specification of list of images that are spaced " << endl;
          cerr << " at a constant time interval (dt-samp). Images may contain one or more data types " << endl;
          cerr << " (planar, surface, iso), but the quantity of interest should be the same for all types " << endl;
          cerr << " (e.g. pressure).  The list of images are presently put in order in a text " << endl;
          cerr << " file with the filenames.  For instance, in a file images.txt" << endl;
          cerr << endl;
          cerr << endl;
          cerr << "       image.00000000.png " << endl;
          cerr << "       image.00000010.png " << endl;
          cerr << "       image.00000020.png " << endl;
          cerr << " where the images are space by a dt = 1.0e-02 (corresponding to a step size of 1.0e-03 in this " << endl;
          cerr << " example. " << endl;
          cerr << endl;
          cerr << endl;
          
          cerr << " Execution of " << endl;
          cerr << endl;
          cerr << " mpirun ./imageDmd.exe images.txt 1.0e-02" << endl;
          cerr << endl;
          cerr << " will compute the dmd and output a file freq_amp.dat.  this contains a list of computed modes " << endl;
          cerr << " with frequency and amplitude information.  It will also output a file A.dat, which contains the " << endl;
          cerr << " time correlation matrix (which is the most time consuming part of the algorithm).  a subsequent " << endl;
          cerr << " execution can extract the mode shapes.  in the example above, if there was a frequency of 100 [1/time unit]" << endl;
          cerr << " present, then in a file mode_list.in, the following can be placed: " << endl;
          cerr << endl;
          cerr << "  100 [true,false] " << endl;
          cerr << endl;
          cerr << endl;
          cerr << " the true false designation alerts the code whether to write an animation of the mode in time. " << endl;
          cerr << " if omitted, it will default to false.  the subsequent execution of the code can be performed with " << endl;
          cerr << endl;
          cerr << " mpirun ./imageDmd.exe images.txt 1.0e-02 A.dat mode_list.in " << endl;
          cerr << endl;
          cerr << " which will now output the mode to a file.  multiple entries (on separate lines) can be provided in " << endl;
          cerr << " the mode_list.in for extracting multiple frequencies. " << endl;
          cerr << endl;
        }//(root)
        throw(0);
      }//(argc < 3) 

      else if (argc == 3) { // ./imageDmd.exe <list.txt> <dt>
        dmd = new PngDmd(argv[1], false);
        dmd->dt_samp = atof(argv[2]);

        // need a way to get the mode list here...
      }//(argc == 3)

      else if ((argc == 4) || (argc == 5)) {// ./imageDmd.exe <list.txt> <dt> <A.dat> <mode_list.txt>
        dmd = new PngDmd(argv[1], argv[3], false); // image file, cached matrix
        dmd->dt_samp = atof(argv[2]);
    
        if (argc == 5) {//mode list is provided
          ifstream fin(argv[4]); // this is a mode file... 
          string line;
          
          while(getline(fin, line)) {
            stringstream ss; ss << line;
            double f_mode;
            string token;
            bool animate = false; // default 

            ss >> f_mode; // should be the mode's frequency

            if (ss >> token) { 
              // then there are additional characters on this line.. 
              if      ((token == "true")  || (token == "TRUE")  || (token == "True"))   animate = true;
              else if ((token == "false") || (token == "FALSE") || (token == "False"))  animate = false;
              else {
                if (mpi_rank == 0) {
                  cout << " > unrecognized token (expecting true/false): " << token << endl;
                  cout << " > defaulting to false" << endl;
                }
              } 
            }//(ss >> token)

            if (mpi_rank == 0) { 
              if (animate)  
                cout << " > adding mode (frequency) : " << f_mode << " with animation " << endl;
              else  
                cout << " > adding mode (frequency) : " << f_mode << " without animation " << endl;
            }//(root)

            dmd->modeList.push_back(pair<double,bool>(f_mode,animate));
          }//(getline(fin, line))
          fin.close(); 
        }//(argc == 5) 
      }//((argc == 4) || (argc == 5))

      else { 
        if(mpi_rank == 0) cerr << "Unhandled " << endl;
        MPI_Abort(mpi_comm,-1);
      } 
  
      if(mpi_rank == 0) cout << " Using dt_samp : " << dmd->dt_samp << endl; 
      dmd->run();  
  
      if(dmd != NULL) {
        delete dmd;
        dmd = NULL;
      }//(dmd != NULL)
    }//scope
  }//try
  catch(...) {
    MPI_Abort(mpi_comm,1);
  }//catch
  MPI_Finalize();
}//main() 
