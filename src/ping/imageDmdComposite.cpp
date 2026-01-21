#include "PngAnalysis.hpp"
#include "Dmd.hpp"
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
};

class PngDmdComposite : public PngAnalysis {
 
public: 
  complex<double>* Vr ;
  double * lambda_i ;
  double * amp;
  string infile0, infile1; 
  Composite* composite_base;

  PngDmdComposite(const string& _infile0, const string& _infile1, const int _crop_start[2], const int _crop_end[2], 
         const bool _bMean) : PngAnalysis(), 
                              Vr(NULL), lambda_i(NULL), amp(NULL), infile0(_infile0), infile1(_infile1),
                              composite_base(NULL) {
    init();
    buildCorrelationMatrix();
  }
  
  PngDmdComposite(const string& _infile0, const string& _infile1, const string &_matfile, const int _crop_start[2], const int _crop_end[2], 
         const bool _bMean) : PngAnalysis(), 
                              Vr(NULL), lambda_i(NULL), amp(NULL), infile0(_infile0), infile1(_infile1), composite_base(NULL)  {
    init();
    readMatrix(_matfile.c_str(),A);
  }
  
  void init() {
    ncache = 15; // >= 2

    int ns0 =0, ns1 = 0;
   
    if (mpi_rank == 0 ) { 
      cout << " init... " << endl;
    }

    { 
      ifstream fin(infile0.c_str()); 
      while ( true ) { 
        string tmp; fin >> tmp; 
        if ( fin.eof() )  break; 
        imageList.push_back(tmp);
        ++ns0;
      }
      fin.close();
    }
    
    {
      ifstream fin(infile1.c_str());
      while ( true ) { 
        string tmp; fin >> tmp;
        if ( fin.eof()) break;
        imageList.push_back(tmp);
        ++ns1;
      }
      fin.close();
    }

    // the stacked vectors must have the same size...
    assert( ns0 == ns1);
    assert( imageList.size() == 2*ns0);
    ns = ns0; 

    if(composite_base==NULL) 
      composite_base = new Composite(imageList[0],imageList[0+ns]);  
    nd = composite_base->phi0->npx + composite_base->phi1->npx; // number of data points.. 

    if ( mpi_rank == 0 ) {  
      cout << " > Ns, nd = " << ns  << " , " << nd << endl;
      cout << " > Mem alloc for matrix (GB) = " << double(ns*ns*8)/double(1024*1024*1024) << endl; 
    }
    
    imageCache.clear(); 
    A = new double[ns*ns] ; // assuming that ns is sufficiently small < 10^4 ..
  }
  
  void buildCorrelationMatrix() { 
  
    if ( mpi_rank == 0) { 
      cout << "building composite correlation matrix... " << endl;
    }

    assert( mean == NULL);
    double * myA = new double[ns*ns]; 
    for (int i =0; i < ns*ns ; ++i) { 
      A[i]   = 0.0; 
      myA[i] = 0.0;
    }//i

    // we have ns^2/2 entries to populate.  the entries will be 
    // distributed to the different ranks.  should be populated 
    // by processing ncache X ncache submatrices in order to make 
    // use of the imageCache... 
    const int bs = max(ncache/2,2); 
    vector<Block*> blocks; 

    for (int j=0; j < ns ; j += bs ) { 
      for (int i =j; i < ns ; i += bs ) { 

        Block* b  = new Block(); 
        b->istart = i; 
        b->iend   = min(i+bs, ns); // NOT inclusive of the end.
        b->jstart = j;
        b->jend   = min(j+bs, ns);
        b->rank   = -1;
        blocks.push_back(b);
      }
    }

    // round-robin the distribution of the blocks now...
    for (int j=0; j < blocks.size(); ++j ) { 
      int r = j%mpi_size; 
      blocks[j]->rank = r;
    }
    
    if ( mpi_rank == 0 ) { 
      cout << " > Building correlation matrix..."; 
      cout.flush(); 
    }


    // operate on your blocks...
    double wtime = MPI_Wtime(); 
    for (int ibl = 0; ibl < blocks.size(); ++ibl) if ( blocks[ibl]->rank == mpi_rank) { 

      Block * b = blocks[ibl]; 
      if ( mpi_rank == 0 ) {  
        cout << "." ; 
        cout.flush(); 
      }

      for (int j = b->jstart; j < b->jend; ++j) { 
        for (int i = b->istart; i < b->iend ; ++i) {

          // for this i,j entry.. we actually need to load in 4 images..
          // this particular way of access ensures that the diagonal
          // entry caching is handled properly... be careful. 
          PngData* pngj0 = getImage(j);
          PngData* pngj1 = getImage(j+ns);
          PngData* pngi0 = getImage(i);
          PngData* pngi1 = getImage(i+ns);

          myA[(j)*ns+i] = dotproduct(pngj0,pngi0) + dotproduct(pngj1,pngi1);
          pngi0->lock = false;
          pngi1->lock = false;
          pngj0->lock = false;
          pngj1->lock = false;
        }//i
      }//j
    }//ibl


    MPI_Allreduce(myA, A, ns*ns, MPI_DOUBLE, MPI_SUM, mpi_comm); 

    // fill out the upper triangle..
    for (int j=0; j < ns ; ++j) { 
      for (int i=0; i < j ; ++i) { 
        A[(j)*ns+i] = A[(i)*ns+j];
      }
    }

    for (int ibl=0; ibl < blocks.size(); ++ibl) 
      delete blocks[ibl];

    delete[] myA; 

    wtime = MPI_Wtime() - wtime; 
    if ( mpi_rank == 0 ) { 
      cout <<  " OK " << endl; 
      cout <<  " > Elapsed time [secs] = " << wtime << endl; 
    }

    if ( mpi_rank == 0) 
      writeMatrix("A.dat", A); 
  }//buildCorrelationMatrix 

  void dump() const { 

    FILE * fp          = fopen("freq_amp.dat", "w"); 
    FILE * f_mode_list = fopen("mode_list.dat", "w");

    set<DmdModeDescriptor> descriptors;
    for (int i =0; i < nnz ; ++i) { 
      if ( lambda_i[i] >= 0.0) {
        // also compute and report the a_i*lambda_k^n where n is the 
        // number of time sequences...
        complex<double> lambda_c(lambda_r[i]*dt_samp,lambda_i[i]*dt_samp);
        complex<double> mu_c     = exp(lambda_c);
        complex<double> mu_n     = pow(mu_c,ns);
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
  
      // HACK.. find me a mode... 
      if ( (lambda_i[i] >= 0.0) && (lambda_i[i] <= 2000.0)) { 
        cout << i << "   " << lambda_i[i] << endl; 
	fprintf(f_mode_list, "%d \n", i);
      }
    }//i

    // write out the frequency, amplitude information in sorted order
    for (set<DmdModeDescriptor>::reverse_iterator it = descriptors.rbegin(); it != descriptors.rend(); ++it)  
      fprintf(fp, "%12.8g    %12.8g    %12.8g    %12.8g    %d\n", it->freq, it->amp_0, it->damp, it->amp_n, it->idx); 

    fclose(fp); 
    fclose(f_mode_list);
  }

  void writeMode(const int imo) { 

    FILE * fp;
    char filename0[128], filename1[128];

    if ( mpi_rank == 0 ) {  
      cout << " > write dmd modes .... " << endl; 
      cout.flush();
      fp = fopen("mode_key.dat", "a");
    }
  
    // allocate the memory .. 
    const int nx0 = composite_base->phi0->getNx(); 
    const int ny0 = composite_base->phi0->getNy();
    const int nx1 = composite_base->phi1->getNx();
    const int ny1 = composite_base->phi1->getNy();
    const int n0  = nx0*ny0;
    const int n1  = nx1*ny1;

    // reconstruct the mode... 
    complex<double> c_coeff; 
    double coeff, max_val; 

    double my_dotp = 0.0, dotp = 0.0;
    double * my_reconstruct0 = new double[n0]; 
    double * my_reconstruct1 = new double[n1];
    double * reconstruct0    = new double[n0]; 
    double * reconstruct1    = new double[n1];

    int * imora = NULL; calcUniformDist(imora, ns-1, mpi_size); 

    for (int iter = 0; iter < 2 ; ++iter ) { 

      for (int ii =0; ii < n0; ++ii) 
        my_reconstruct0[ii] = 0.0;

      for (int ii =0; ii < n1; ++ii) 
        my_reconstruct1[ii] = 0.0;

      if (mpi_rank==0) {
        if ( iter == 0 ) {  
          sprintf(filename0, "mode0.%03d.real.png", imo);
          sprintf(filename1, "mode1.%03d.real.png", imo);
        } else { 
          sprintf(filename0, "mode0.%03d.imag.png", imo);
          sprintf(filename1, "mode1.%03d.imag.png", imo);
        }

        cout << "Working on reconstructing : " << filename0 << "; freq : " << lambda_i[imo] <<  endl; 
    
        if ( iter == 0 ) 
          fprintf(fp, "%03d  %12.8g   ", imo, lambda_i[imo]); 
      }

      //for (int i =0; i < ns-1 ; ++i) { 
      for (int i = imora[mpi_rank]; i < imora[mpi_rank+1]; ++i) { 

        PngData *png0 = getImage(i); 
        PngData *png1 = getImage(i+ns);
        reconstructCoeff(c_coeff, Wk, Vr, nnz, ns, imo, i); 
        coeff = (iter == 0) ? c_coeff.real() : c_coeff.imag(); 

        //cout << i << "   " << coeff << endl; 
        for (int ipx = 0; ipx < n0; ++ipx) { 
          if ( composite_base->phi0->pixel_flag[ipx] >= 0) 
            my_reconstruct0[ipx] += coeff* png0->getPhi(ipx);
        }

        for (int ipx =0; ipx < n1; ++ipx) { 
          if ( composite_base->phi1->pixel_flag[ipx] >= 0) 
            my_reconstruct1[ipx] += coeff* png1->getPhi(ipx);
        }
       
        png0->lock = false;
        png1->lock = false;
      }//i

      MPI_Allreduce(my_reconstruct0, reconstruct0, n0, MPI_DOUBLE, MPI_SUM, mpi_comm);
      MPI_Allreduce(my_reconstruct1, reconstruct1, n1, MPI_DOUBLE, MPI_SUM, mpi_comm);

      for (int ipx = 0; ipx < n0; ++ipx) { 
        if ( composite_base->phi0->pixel_flag[ipx] >= 0)  
          my_dotp += reconstruct0[ipx] * reconstruct0[ipx]/double(composite_base->phi0->npx); 
      }

      for (int ipx = 0; ipx < n1 ; ++ipx) { 
        if ( composite_base->phi1->pixel_flag[ipx] >= 0) 
          my_dotp += reconstruct1[ipx] * reconstruct1[ipx]/double(composite_base->phi1->npx);
      }
      
      if (mpi_rank==0){

        PngData * png0 = new PngData(crop_start, crop_end);
        png0->read(imageList[0].c_str(),true);
        png0->finalizeRead();
        png0->buffer = new unsigned char[nx0*ny0*png0->getStride()];

        PngData * png1 = new PngData(crop_start, crop_end); 
        png1->read(imageList[ns].c_str(),true);
        png1->finalizeRead();
        png1->buffer = new unsigned char[nx1*ny1*png1->getStride()];

        // now set the buffer...
        double range0[2], range1[2]; 
        dumpRange(range0, reconstruct0, nx0*ny0, "RECONSTRUCT0");
        dumpRange(range1, reconstruct1, nx1*ny1, "RECONSTRUCT1");

        const double eps = 1.0e-08;
        if ( (fabs(range0[0]) < eps) && ( fabs(range0[1]) < eps)) { 
          range0[0] = -eps; 
          range0[1] =  eps;
        }

        if ( (fabs(range1[0]) < eps) && ( fabs(range1[1]) < eps)) { 
          range1[0] = -eps; 
          range1[1] =  eps;
        }

        fprintf(fp, "  %12.8g   %12.8g ", range0[0], range0[1]); 
        fprintf(fp, "  %12.8g   %12.8g ", range1[0], range1[1]); 
  
        for (int ipx = 0; ipx < nx0*ny0 ; ++ipx ) { 
          if ( png0->pixel_flag[ipx] >= 0 ) { 
            const int r = int(255.0* (reconstruct0[ipx] - range0[0])/(range0[1]-range0[0]));
            assert ( (r>=0)&&(r<= 255)) ;
            //cout << ipx << "   " << r << endl;
            png0->buffer[png0->getStride()*ipx+0] = (unsigned char)(r) ; 
            png0->buffer[png0->getStride()*ipx+1] = (unsigned char)(r) ; 
            png0->buffer[png0->getStride()*ipx+2] = (unsigned char)(r) ; 
          } else {
            png0->buffer[png0->getStride()*ipx+0] = 73;
            png0->buffer[png0->getStride()*ipx+1] = 175;
            png0->buffer[png0->getStride()*ipx+2] = 205;
          }
        }
        
        png0->setDepthFromPixelFlag();
        png0->range[0] = range0[0];
        png0->range[1] = range0[1];
        png0->write(filename0); 

        for (int ipx = 0; ipx < nx1*ny1 ; ++ipx ) { 
          if ( png1->pixel_flag[ipx] >= 0 ) { 
            const int r = int(255.0* (reconstruct1[ipx] - range1[0])/(range1[1]-range1[0]));
            assert ( (r>=0)&&(r<= 255)) ;
            //cout << ipx << "   " << r << endl;
            png1->buffer[png1->getStride()*ipx+0] = (unsigned char)(r) ; 
            png1->buffer[png1->getStride()*ipx+1] = (unsigned char)(r) ; 
            png1->buffer[png1->getStride()*ipx+2] = (unsigned char)(r) ; 
          } else {
            png1->buffer[png1->getStride()*ipx+0] = 73;
            png1->buffer[png1->getStride()*ipx+1] = 175;
            png1->buffer[png1->getStride()*ipx+2] = 205;
          }
        }
        
        png1->setDepthFromPixelFlag();
        png1->range[0] = range1[0];
        png1->range[1] = range1[1];
        png1->write(filename1); 

        delete png0;
        delete png1;
      }
    }//iter  
  
    MPI_Allreduce(&my_dotp, &dotp, 1, MPI_DOUBLE, MPI_SUM, mpi_comm);
    if ( mpi_rank==0) {
      fprintf(fp, "\n"); 
      fclose(fp);
      cout << " Mode normalization across " << mpi_size << " ranks: " << dotp/(double)mpi_size << endl; 
    }
   
    delete[] imora;
    delete[] reconstruct0;
    delete[] reconstruct1;
    delete[] my_reconstruct0;
    delete[] my_reconstruct1;
  } 

  void dumpRange(double range[2], const double * buf, const int nn, const char* msg) { 
    range[0] = 1.0e+16; 
    range[1] = -1.0e+16; 

    for (int i =0; i < nn ; ++i) { 
      range[0] = fmin(range[0], buf[i]); 
      range[1] = fmax(range[1], buf[i]); 
    }

    cout << "dumpRange, " << msg << " = " << range[0] << ": " << range[1] << endl;
  }

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

      nnz = dmdAnalysis(ns, A, lambda_r, lambda_i, Wk, Vr, amp, dt);     

      dump();
    }
    else{
      Wk = new double[(ns-1)*(ns-1)];
    }

    // everyone needs the Wk, Vr, nnz
    MPI_Bcast(&nnz, 1, MPI_INT, 0, mpi_comm);     

    if (mpi_rank !=0){  
      Vr = new complex<double>[nnz*nnz];
      lambda_i = new double[nnz];
    }

    // everyone needs the Wk, Vr, nnz, lambda_i for print statement..
    MPI_Bcast(Wk, (ns-1)*(ns-1), MPI_DOUBLE, 0, mpi_comm);
    MPI_Bcast(Vr, nnz*nnz, MPI_DOUBLE_COMPLEX, 0, mpi_comm);
    MPI_Bcast(lambda_i, nnz, MPI_DOUBLE, 0, mpi_comm);

    for (int i=0; i < modeList.size(); ++i)
      writeMode(modeList[i]);
  }//runDmd

  void run() { runDmd(); }

  ~PngDmdComposite() {

    if(mpi_rank==0) cout << "~PngDmd()...";

    if ( lambda_i != NULL )
      delete[] lambda_i ; 

    if ( amp != NULL)
      delete[] amp;

    if( Vr != NULL)
      delete[] Vr;

    if(mpi_rank==0) cout << "OK" << endl;
  }

  template <class T>
  bool from_string(T& t, const std::string &s, std::ios_base& (*f)(std::ios_base&))
  {
     std::istringstream iss(s);
     return !(iss >> f >> t).fail();
  }
}; 


int main(int argc, char* argv[]) { 

  initMpiEnv(argc,argv);

  int crop_start[2] = {0, 0};
  int crop_end[2]   = {-1, -1};
  
  PngDmdComposite * dmd = NULL; 
  
  if ( argc < 4 ) { 
    if ( mpi_rank==0) cerr << " Usage: ./imageDmdComposite.exe <image-list-file0> <image-list-file1> <dt-samp> [<cached-matrix-file> <mode-list-file>] " << endl; 
    return -1; 
  } else if (argc == 4) { 
    
    try {
      dmd = new PngDmdComposite(argv[1], argv[2], crop_start, crop_end, false);
      dmd->dt_samp = atof(argv[3]);
    }
    catch(int e){
      return -1;
    }

    // need a way to get the mode list here...
  } else if ( (argc == 5)||(argc == 6) ) { 

    try{
      dmd = new PngDmdComposite(argv[1], argv[2], argv[4], crop_start, crop_end, false); // image file, cached matrix
      dmd->dt_samp = atof(argv[3]);
    }
    catch(int e) {
      return -1;
    }
    
    if ( argc == 6 ) { 
      ifstream fin(argv[5]); // this is a mode file... 
      while ( true) { 
        string s_mode; fin >> s_mode; 

        if (fin.eof())
          break;
        
        const int i_mode = atoi(s_mode.c_str()); 
        if(mpi_rank==0) cout << " Adding mode : " << i_mode << endl; 
        dmd->modeList.push_back(i_mode);
      }
      fin.close(); 
    }
  } else { 
    if(mpi_rank==0) cerr << "Unhandled " << endl;
    return -1;
  } 
  
  if(mpi_rank==0) cout << " Using dt_samp : " << dmd->dt_samp << endl; 
  dmd->run();  
  
  if(dmd!=NULL) {
    delete dmd;
    dmd=NULL;
  }

  MPI_Finalize();
  return 0; 
} 
