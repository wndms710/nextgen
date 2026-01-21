#include "PngAnalysis.hpp"
#include "Pod.hpp"

class ModeLight { 
public: 
  double * phi; 
  double norm ; 
  int n; 
  
  ModeLight() : phi(NULL), norm(0.0), n(0) {} 
  ModeLight(const int _n) : norm(0.0), n(_n) { 
    phi = new double[_n]; 
    for (int i =0; i < _n ; ++i) 
      phi[i] = 0.0;
  }
  ~ModeLight() { 
    if ( phi != NULL ) 
      delete[] phi;
  }
};//class ModeLight

class PngPod : public PngAnalysis {

public:

  PngPod(const string& _infile, const bool _bMean) : PngAnalysis(_infile, _bMean) { }

  PngPod(const string& _infile, const string &_matfile, const bool _bMean) : PngAnalysis(_infile, _matfile, _bMean) { }

  void runPod() { 
    
    
    // lambda_r stores the singular values...
    if (mpi_rank == 0) 
      nnz = podAnalysis(n_snaps, A, lambda_r, Wk);
    else 
      Wk = new double[(n_snaps-1)*(n_snaps-1)];

    // everyone needs the Wk... 
    MPI_Bcast(Wk, (n_snaps-1)*(n_snaps-1), MPI_DOUBLE, 0, mpi_comm);  
    MPI_Bcast(&nnz, 1, MPI_INT, 0, mpi_comm); 

    const int r = 40;
    if (mpi_rank == 0) 
      cout << "podFilter" << endl;
    podFilter(min(nnz,r)); 

  }//runPod()

  void podFilter(const int _r) { 
    
    // reconstruct the time history from the first r pod modes... 
    assert ( Wk != NULL ) ; 
    assert ( A  != NULL ) ; 

    reconstructPodModes(_r); 

    if (mpi_rank == 0) 
      cout << "recontruction complete!" << endl;

    // show me the pod modes and coeffs from rank 0... 
    if ( mpi_rank == 0 ) { 
      cout << "Down here" << endl;   
      // and the coefficients.... 
      double * Tri = new double[(n_snaps-1)*_r]; 
      for (int i = 0; i < n_snaps-1; ++i) { 
        for (int imo=0; imo < _r; ++imo) { 
          
          const int jj = nnz-1-imo;
          Tri[(i)*_r + imo] = 0.0; 
          for (int m =0; m < n_snaps-1; ++m) 
            Tri[(i)*_r + imo] += Wk[(jj)*(n_snaps-1)+m]*A[(m)*n_snaps+i]; 
        }//imo
      }//i
      cout << "Did that" << endl;
      
      FILE * fp = fopen("pod_coeffs.dat", "w"); 
      for (int i=0; i < n_snaps-1; ++i) { 
        for (int imo=0; imo < _r ; ++imo) { 
          fprintf(fp, "%12.8g  ", Tri[(i)*_r + imo]); 
        }
        fprintf(fp, "\n");
      }//i
      
      fclose(fp);
      cout << "DONEZO" << endl;

      delete[] Tri;
    }
 
  }//podFilter() 
  
  void reconstructPodModes(const int _r) { 

    if (mpi_rank == 0) {  
      cout << " > Reconstructing the pod modes .... " << endl; 
      cout.flush();
    }//(root)

    // allocate the memory .. 
    const int nx = base->getNx(); 
    const int ny = base->getNy(); 
    const int n  = nx*ny;

    // its more efficient to aggregate the modes together, but 
    // we may not have the memory to keep all the modes at once.
    // XXX can batch these as well, for now; compute them indv..
    double * my_phi = new double[n];
    double * phi    = new double[n];

    int * imora = NULL; calcUniformDist(imora, n_snaps-1, mpi_size); 

    for (int imo=0; imo < _r; ++imo ) { 

      const int jj = nnz-1-imo; // column of Wk that has the eigenvector coeffs.
      
      for (int i=0; i < n; ++i) 
        my_phi[i] = 0.0; 
      
      for (int i=imora[mpi_rank]; i < imora[mpi_rank+1]; ++i) { 
        PngData * png = getImage(i); 
        for (int ipx =0; ipx < n ; ++ipx) { 
          if (base->pixel_flag[ipx] >= PngData::VOLUME_DATA_PIXEL ) { 
            my_phi[ipx] +=  Wk[(jj)*(n_snaps-1)+i] * png->pixel_data[ipx]; //Wk[(jj)*(n_snaps-1)+i] * png->getPhi(ipx);
          }
        }
        png->lock = false;
      }//i

      MPI_Allreduce(my_phi, phi, n, MPI_DOUBLE, MPI_SUM, mpi_comm); 

      double norm = 0.0; 
      for (int ipx = 0; ipx < n ; ++ipx) {
        //note that bahama blue entries should have 0.0 in their entires...
        norm += phi[ipx] * phi[ipx]/double(base->npx+base->npxS+base->npxP+base->npxI); 
      }
      //assert ( norm > 0.0);
      
      norm = sqrt(norm); 
      for (int ipx =0; ipx < n ; ++ipx) 
        phi[ipx] /= norm;


      if (mpi_rank == 0 ) { //dump...

        cout << " norm  = " << norm << endl; 
        
        char filename[128]; sprintf(filename, "pod_mode.%04d.png", imo); 
        PngData * png = new PngData(); 
        png->read(imageList[0].c_str()); 
        png->finalizeRead(); 
        
        cout << " ... Working on pod mode " << imo  << " of " << _r << endl;
        
        const int nx = png->getNx(); 
        const int ny = png->getNy(); 
        
        double range[2]; dumpRange(range, phi, nx*ny, "PHI"); 
        
        for (int ipx =0; ipx < nx*ny ; ++ipx)  
          if (png->pixel_flag[ipx] >= PngData::VOLUME_DATA_PIXEL)  
            png->pixel_data[ipx] = phi[ipx];


        //all data types are operated on as one, so set the same range
        //for all types, initializeWrite will write the applicable ranges
        //for the data types that exist to image metadata.
        png->setRange(range);
        png->setRangeSurface(range);
        png->setRangeParticles(range);
        png->setRangeIso(range);

        png->initializeWrite();
        png->write(filename); 
        delete png;
      }//(root)
    }//(imo < _r)

    delete[] my_phi; 
    delete[] phi; 


    if (mpi_rank == 0) { 
      cout << " OK. " << endl; 
      cout.flush(); 
    }

    /*
    // check the orthogonality....
    cout << " Checking orthogonality ... " << endl; 
    for (int imo = 0; imo < _r ; ++imo) { 
      for (int imo2 = 0; imo2 < _r ; ++imo2 ) { 

        double dotp = 0.0; 
        for (int ipx = 0; ipx < n ; ++ipx) { 
          dotp += podModes[imo]->phi[ipx] * podModes[imo2]->phi[ipx]/double(npx);
        }//ipx
      
        cout << " < " << imo << " , " << imo2 << " > = " << dotp << endl;
      }
    }
    */


  }//reconstructPodModes 

  void dumpRange(double range[2], const double * buf, const int nn, const char* msg) { 

    range[0] = 1.0e+16; 
    range[1] = -1.0e+16; 

    for (int i =0; i < nn ; ++i) { 
      range[0] = fmin(range[0], buf[i]); 
      range[1] = fmax(range[1], buf[i]); 
    }

    cout << "dumpRange, " << msg << " = " << range[0] << ": " << range[1] << endl;
  }

  void run() { 
    runPod(); 
  }
};//class PngPod



int main(int argc, char* argv[]) { 

  try {

    //string infile = "crop.in";
    // this call to CTI_Init will update the infile based on any command line input
    initMpiEnv(argc,argv);

    {//scope

      /*
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

      PngPod * pod = NULL ; 
      if (argc <= 1) { // No arguments provided
        if ( mpi_rank == 0 ) {
          cerr << "Too few arguments provided!" << endl;
          cerr << " Usage: ./imagePod.exe <image-list-file> [<matrix-file>] " << endl; 
        }//(root)
        throw(0); 
      } 
      else if (argc == 2) { // list provided
        pod = new PngPod(argv[1], true);
      }  
      else if (argc == 3) { // list + matrix provided
        pod = new PngPod(argv[1], argv[2], true); //image file, cached matrix..
      } 
      else { 
        if (mpi_rank == 0) { 
          cerr << "Too many arguments provided!" << endl;
          cerr << " Usage: ./imagePod.exe <image-list-file> [<matrix-file>] " << endl; 
        }//(root)
        throw(0); 
      }

      pod->run(); 
      delete pod;
      if (mpi_rank == 0) cout << "DONZEO! Exiting cleanly" << endl;

    }//scope
  }//try
  catch(...) {
    MPI_Abort(mpi_comm,1);
  }
  MPI_Finalize();
}//main... 

