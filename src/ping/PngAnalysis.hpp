#ifndef __PNG_ANALYSIS__
#define __PNG_ANALYSIS__

#ifdef NO_MPI
#include "NoMpi.hpp"
#else
#include <mpi.h>
#endif

#include <fstream>
#include <ctime>
#include <vector>
#include <deque>
#include <cmath>
#include <sstream>
#include <assert.h>
#include <map>
#include <algorithm>


typedef unsigned short uint2;

namespace MpiEnv {
  int mpi_rank;
  int mpi_size;
  MPI_Comm mpi_comm;
};//MpiEnv

using namespace MpiEnv;

void initMpiEnv(int argc, char* argv[]) {
  MPI_Init(&argc,&argv);
  MPI_Comm_dup(MPI_COMM_WORLD, &mpi_comm);
  MPI_Comm_rank(mpi_comm, &mpi_rank);
  MPI_Comm_size(mpi_comm, &mpi_size);
}//initMpiEnv()

// XXX JOB This needs to move down here (from line 19) b/c it depends on mpi_rank
#include "PngData.hpp"

void calcUniformDist(int * &xod, const int nx, const int ndist) {
  if (xod == NULL) xod = new int[ndist+1];
  xod[0] = 0;
  for (int id = 1; id <= ndist; id++) {
    xod[id] = (int)((double)id/(double)ndist*(double)nx + 0.5);
  }
}//calcUniformDist()

struct PngDataComp {
  int n;
  explicit PngDataComp(int i) : n(i) { }
  inline bool operator()(const PngData* png) const { return png->id == n; }
};//PngDataComp

struct Block {
  int istart,iend;
  int jstart,jend;
  int rank;
};//Block

class PngAnalysis {

public:
  int n_snaps;
  int n_pts;
  double* A;
  string infile ;
  
  // images..
  int n_cache;
  deque<PngData*> imageCache;
  vector<string> imageList;
  PngData* base;
  PngData* mean;
  double * mean_d;
  
  // pod stuff
  double * Wk;
  double * lambda_r ;
  int nnz ;
  double dt_samp ;
  //vector<int> modeList;
  vector<pair<double,bool> > modeList;

  // image metadata
  //double rmin;
  //double rmax;
  //string varName;

  bool bMean;
 
  PngAnalysis() : base(NULL),mean(NULL),mean_d(NULL),Wk(NULL),lambda_r(NULL),nnz(0),dt_samp(0.0),bMean(false) {}

  PngAnalysis(const string& _infile, const bool _bMean) : infile(_infile),base(NULL),mean(NULL),mean_d(NULL),Wk(NULL),lambda_r(NULL),nnz(0),dt_samp(0.0),bMean(_bMean) {
    commonInit();
    buildCorrelationMatrix();
  }//PngAnalysis()

  PngAnalysis(const string& _infile, const string &_matfile, const bool _bMean) : infile(_infile),base(NULL),mean(NULL),mean_d(NULL),Wk(NULL),lambda_r(NULL),nnz(0),dt_samp(0.0),bMean(_bMean)  {
    
    commonInit();
    readMatrix(_matfile.c_str(),A);
  }//PngAnalysis()

  PngData* getImage(const int i) { 

    deque<PngData*>::iterator it = find_if(imageCache.begin(), imageCache.end(), PngDataComp(i)); 
    
    // The image is already in the cache...
    if (it != imageCache.end()) { 
      (*it)->lock = true;
      return *it; 
    }//(it != iamgeCache.end()) 
    // The image is not currently in the cache...
    else {
      if (imageCache.size() == n_cache) { 
        // cache is full, remove the first unlocked element (first in-first out)
        for(it = imageCache.begin(); it != imageCache.end(); ++it) { 
          if (!((*it)->lock)) 
            break;
        }//it

        assert( it != imageCache.end()); // uh oh, the whole cache is locked...
        delete *it; 
        imageCache.erase(it); 
      }//(imageCache.size() == n_cache)
     
      // when we get the image, we will lock it, but 
      // the user is responsible for releasing the lock
      PngData * png = new PngData();

      png->read(imageList[i].c_str());
      png->finalizeRead();
      png->id = i;
      png->lock = true; 
      imageCache.push_back(png);

      assert( imageCache.size() <= n_cache); 

      return png;
    }//(it == imageCache.end())
        
    // shouldnt get here..
    return NULL; 
  }//getImage() 

  void buildMean() { 

    mean = new PngData();
    mean->read(imageList[0].c_str()); 
    mean->finalizeRead(); 

    const int nx = mean->getNx(); 
    const int ny = mean->getNy(); 

    mean_d = new double[nx*ny];
    double * my_mean_d = new double[nx*ny]; 
    for (int ipx=0; ipx < nx*ny; ++ipx)  
      my_mean_d[ipx] = 0.0; 
   
    // stripe the mean (n-1 images) across the ranks and reduce.
    const int nimg = n_snaps-1; 
    int * imora = NULL; calcUniformDist(imora, nimg, mpi_size); 

    if ( mpi_rank == 0 ) { 
      cout << " > Computing mean: " << endl;
      cout.flush(); 
    }

    for (int j = imora[mpi_rank]; j < imora[mpi_rank+1]; ++j) { 

      // Print progress
      if (mpi_rank == 0) {
        int progress = int(round(1000*double(j - imora[mpi_rank])/double(imora[mpi_rank+1] - imora[mpi_rank]))/10.0);
        cout << "                   " << progress << " % complete" << endl;
      }

      PngData* png = getImage(j); 
      for (int ipx = 0; ipx < nx*ny; ++ipx) { 
        if (mean->pixel_flag[ipx] >= PngData::VOLUME_DATA_PIXEL) { 
          assert(png->pixel_flag[ipx] >= PngData::VOLUME_DATA_PIXEL); 
          my_mean_d[ipx]     += png->pixel_data[ipx]; //png->getPhi(ipx);
        }//VOLUME_DATA_PIXEL
      }//ipx
      png->lock = false;
    }//j

    // now reduce the mean_d buffer..
    MPI_Allreduce(my_mean_d, mean_d, nx*ny, MPI_DOUBLE, MPI_SUM, mpi_comm); 

    // and normalize our entries... 
    // and update gray value in mean PngImage
    for (int ipx = 0; ipx < nx*ny ; ++ipx) { 
      if (mean->pixel_flag[ipx] >= PngData::VOLUME_DATA_PIXEL) { 
        mean_d[ipx] /= double(nimg); 
        mean->pixel_data[ipx] = mean_d[ipx];
      }//VOLUME_DATA_PIXEL 
    }//ipx

    //write mean image
    if (mpi_rank == 0) { 
      mean->rescaleRangesToData();
      mean->initializeWrite();
      mean->write("mean.png"); 
      cout << "                   100% complete" << endl;
    }//(root)

    delete[] my_mean_d;
    delete[] imora;
  }//buildMean()

  void commonInit() { 
    // number of images that each rank can hold in memory
    n_cache = 10; // must be >= 2 

    // everybody needs to keep a list of images. 
    ifstream fin(infile.c_str()); 
    while (true) { 
      string tmp; fin >> tmp; 
      if ( fin.eof() )  
        break; 
      imageList.push_back(tmp);
    }//(true) 
    fin.close();

    n_snaps = imageList.size(); // number of snapshots..

    // everybody also keeps a base image so we know
    // how many pixels, which pixels are active, etc...
    if(base==NULL) base = new PngData();
    base->read(imageList[0].c_str()); 
    base->finalizeRead(); 
    n_pts = base->npx+base->npxS+base->npxP+base->npxI;        // number of data points.. 

    if ( mpi_rank == 0 ) {  
      cout << " > n_snaps, n_points/image = " << n_snaps  << " , " << n_pts << endl;
      cout << " > Mem alloc for matrix (GB) = " << double(n_snaps*n_snaps*8)/double(1024*1024*1024) << endl; 
    }
    
    imageCache.clear(); 

    A = new double[n_snaps*n_snaps] ; // assuming that n_snaps is sufficiently small < 10^4 ..

    mean                   = NULL; 
    mean_d                 = NULL; 

    if (bMean) buildMean();  
  }//commonInit()
  
  //Check the following:
  // 1. Image variable matches
  // 2. Image dt matches with some tolerance
  // 3. Image RANGE_MIN and RANGE_MAX match within some tolerance (for now)
  /*
  void metadataProcessing(){
     if (mpi_rank == 0){
       cout << "Checking image metadata..." << endl;
     }

     PngData * pngj;
     string varId;
     double reltol = 0.01;
     double timej0;

     for ( int j = 0; j<n_snaps; j++){
       int myError = 0;
       if (mpi_rank == 0){
         ImageMetadata mdj;
         pngj = new PngData(crop_start, crop_end);
         pngj->read(imageList[j].c_str(),false);
         string varIdj;
         double dtj, timej, rminj, rmaxj;

         if (j==0){
           varId  = pngj->pngGrayMetadata.varId;
           from_string(rmin,pngj->pngGrayMetadata.rangeMin,std::dec);
           from_string(rmax,pngj->pngGrayMetadata.rangeMax,std::dec);
           from_string(timej0,pngj->pngGrayMetadata.time,std::dec);
         }
         else{
           varIdj = pngj->pngGrayMetadata.varId;
           from_string(rminj,pngj->pngGrayMetadata.rangeMin,std::dec);
           from_string(rmaxj,pngj->pngGrayMetadata.rangeMax,std::dec);
           from_string(timej,pngj->pngGrayMetadata.time,std::dec);
           if (j==1){
             dt_samp = timej - timej0;
           }
           dtj = timej - timej0;
           timej0 = timej;

           if (varIdj!=varId){
             cout << " > Error: Non matching variables: " << varId << ", " << varIdj << endl;
             myError++;
           }
           if (rminj!=rmin){
             cout << " > Error: Non matching variable range mininum: " << rmin << ", " << rminj << endl;
             myError++;
           }
           if (rmaxj!=rmax){
             cout << " > Error: Non matching variable range maximum: " << rmax << ", " << rmaxj << endl;
             myError++;
           }
           if (fabs(dtj-dt_samp) > reltol*dt_samp){ //for rounding error, all 1-percent difference in time step
             cout << " > Error: Variable time step found: " << dtj << ", " << dt_samp << endl;
             myError++;
           }

           if (myError>0){
             cout << " > Error found at image " << imageList[j] << endl;
           }
         }
         delete pngj; 
       }

       MPI_Bcast(&myError, 1, MPI_INT, 0, mpi_comm);
       if (myError > 0)
         throw(0);

     }
     
     MPI_Bcast(&dt_samp, 1, MPI_DOUBLE, 0, mpi_comm);
     MPI_Bcast(&rmin, 1, MPI_DOUBLE, 0, mpi_comm);
     MPI_Bcast(&rmax, 1, MPI_DOUBLE, 0, mpi_comm);

     char varIdBuff[128];
     sprintf(varIdBuff,"%s",varId.c_str());
     MPI_Bcast(varIdBuff, 128, MPI_CHAR, 0, mpi_comm);
     varId = varIdBuff;

     if (mpi_rank == 0){
       cout << "Image metadata check successful," << endl; 
       cout << "  found VAR=" << varId << ", DT=" << dt_samp << ", RANGE=" << rmin << "," << rmax << endl;  
     }
        
  }
*/

  virtual void buildCorrelationMatrix() { 
    
    double * myA = new double[n_snaps*n_snaps]; 
    for (int i =0; i < n_snaps*n_snaps ; ++i) { 
      A[i]   = 0.0; 
      myA[i] = 0.0;
    }//(i < n_snaps^2)

    // we have n_snaps^2/2 entries to populate.  the entries will be 
    // distributed to the different ranks.  should be populated 
    // by processing n_cache X n_cache submatrices in order to make 
    // use of the imageCache... 
    const int nx     = base->getNx(); 
    const int ny     = base->getNy(); 
    const int npx    = base->npx+base->npxS+base->npxP+base->npxI; 
//    const int * flag = base->pixel_flag; //DAP

    double mdotm   = 0.0;
    double * mdotp = NULL; 
    if (mean != NULL) {
      // Compute <mean, mean>/npx and store in mdotm
      // I think non-volume pixels break this, b/c npx is unaware whether all pixels are valid or not
      for (int ipx = 0; ipx < nx*ny; ++ipx) 
        mdotm += mean_d[ipx]*mean_d[ipx]/double(npx);

      // stripe the images across the ranks to build mdotp
      mdotp           = new double[n_snaps]; 
      double * my_buf = new double[n_snaps];
      for (int i=0; i < n_snaps ; ++i) 
        my_buf[i] = 0.0;

      int * imora     = NULL; 
      calcUniformDist(imora, n_snaps, mpi_size); 

      // Compute <png, mean>/npx for each png and store in mdotp
      for (int j = imora[mpi_rank]; j < imora[mpi_rank+1]; ++j) { 
        PngData* png = getImage(j); 
        my_buf[j]    = dotproduct(mean_d,png);
        png->lock    = false;
      }

      MPI_Allreduce(my_buf,mdotp,n_snaps, MPI_DOUBLE, MPI_SUM, mpi_comm); 
      delete[] my_buf; 
      delete[] imora;
    }//(mean != NULL)
    
    // build the blocks...
    const int blk_sz = max(n_cache/2,2); 
    vector<Block*> blocks; 

    int my_n_blks = 0; // counter for rank's workload
    for (int j=0; j < n_snaps ; j += blk_sz) { 
      for (int i =j; i < n_snaps ; i += blk_sz) { 
        Block* b  = new Block(); 
        b->istart = i; 
        b->iend   = min(i+blk_sz, n_snaps); // NOT inclusive of the end.
        b->jstart = j;
        b->jend   = min(j+blk_sz, n_snaps);
        //b->rank   = -1;
        b->rank   = blocks.size() % mpi_size;
        // round-robin the distribution of the blocks now...
        if (b->rank == mpi_rank) my_n_blks++;
        blocks.push_back(b);
      }//(i < n_snaps)
    }//(j < n_snaps)

    if (mpi_rank == 0)  
      cout << " > Building correlation matrix..." << endl; 

    // operate on your blocks...
    double wtime = MPI_Wtime(); 
   
    int my_c_blk = 0;// count root's number of blocks to gauge progress
    for (int ibl = 0; ibl < blocks.size(); ++ibl) {
      if (blocks[ibl]->rank == mpi_rank) { 
        Block * b = blocks[ibl]; 

        // Print progress
        if (mpi_rank == 0) {  
          int progress = int(round(1000*double(my_c_blk)/double(my_n_blks))/10.0);
          cout << "                   " << progress << " % complete" << endl;
        }//(root)
        my_c_blk++;

        for (int j = b->jstart; j < b->jend; ++j) { 
          for (int i = b->istart; i < b->iend ; ++i) {
            // this particular way of access ensures that the diagonal
            // entry caching is handled properly... be careful. 
            PngData* pngj = getImage(j); 
            PngData* pngi = getImage(i); 

            double dotp = dotproduct(pngi, pngj);
            // This is equivalent to computing <pngi - mean, pngj - mean> 
            if (mean != NULL) { 
              dotp += mdotm ; 
              dotp -= mdotp[i]; 
              dotp -= mdotp[j];
            }//(mean != NULL)
            myA[(j)*n_snaps+i] = dotp ;
            pngi->lock = false;
            pngj->lock = false;
          }//(i < b->iend)
        }//(j < b->jend)
      }//(blocks[ibl]->rank == mpi_rank)
    }//(ibl < blocks.size())

    MPI_Allreduce(myA, A, n_snaps*n_snaps, MPI_DOUBLE, MPI_SUM, mpi_comm); 
    if (mpi_rank == 0) cout << "                   " << 100 << " % complete" << endl;

    // fill out the upper triangle..
    for (int j=0; j < n_snaps ; ++j) { 
      for (int i=0; i < j ; ++i) { 
        A[(j)*n_snaps+i] = A[(i)*n_snaps+j];
      }
    }

    for (int ibl=0; ibl < blocks.size(); ++ibl) 
      delete blocks[ibl];

    //delete[] blora; 
    delete[] myA; 

    if (mean != NULL) 
      delete[] mdotp; 

    wtime = MPI_Wtime() - wtime; 
    if (mpi_rank == 0) { 
      cout <<  " OK " << endl; 
      cout <<  " > Elapsed time [secs] = " << wtime << endl; 
    }//(root)

    // if the mean has been removed, check the row, col sums... 
    if ((mean != NULL)&&(mpi_rank==0)) { 
      // column sum of the first n-1 entries should be zero...
      for (int j=0; j < n_snaps-1; ++j) { 
        double sum = 0.0;
        for (int i=0; i < n_snaps-1 ; ++i) { 
          sum += A[(j)*n_snaps+i]; 
        }
        if (fabs(sum) > 1.0e-10) { 
          cout << "col sum failed: " << j << "    " << sum/A[(j)*n_snaps+j] << endl;
        }
      }//j
      
      // row sum of the first n-1 entries should be zero...
      for (int i=0; i < n_snaps-1; ++i) { 
        double sum = 0.0;
        for (int j=0; j < n_snaps-1 ; ++j) { 
          sum += A[(j)*n_snaps+i]; 
        }
        if ( fabs(sum) > 1.0e-10) { 
          cout << "row sum failed: " << i << "    " << sum/A[(i)*n_snaps+i] << endl;
        }
      }//i
        
      cout << " > Done with row/col sums tests." << endl;
    }

    if ( mpi_rank == 0) 
      writeMatrix("A.dat", A); 

  }//buildCorrelationMatrix()

  void writeMatrix(const char* filename, const double* A) { 

    FILE * fp = fopen(filename, "wb");
    fwrite(&n_snaps, sizeof(int), 1, fp); 
    fwrite(&n_pts, sizeof(int), 1, fp);
    fwrite(A, sizeof(double), n_snaps*n_snaps, fp); 
    fclose(fp); 
  }//writeMatrix()

  void readMatrix(const char* filename, double* A) { 

    FILE * fp = fopen(filename, "rb"); 
    //int c_fpc_size, c_n_snaps, c_n_pts, c_nnz;  These don't get used
    int c_n_snaps, c_n_pts;

    fread(&c_n_snaps,sizeof(int),1,fp); assert( c_n_snaps == n_snaps); 
    fread(&c_n_pts,sizeof(int),1,fp); assert( c_n_pts == n_pts); 
    fread(A, sizeof(double), n_snaps*n_snaps, fp);
  }//readMatrix()

  template <class T>
  bool from_string(T& t, const std::string &s, std::ios_base& (*f)(std::ios_base&)) {
     std::istringstream iss(s);
     return !(iss >> f >> t).fail();
  }

  virtual ~PngAnalysis() {
      
    if(base!=NULL) delete base;
    if (mean != NULL) delete mean;
    if(mean_d!=NULL) delete[] mean_d;
    
    imageList.clear();
      
    for (int i=0; i < imageCache.size(); ++i) {
      if(imageCache[i]!=NULL)	{
        delete imageCache[i];
        imageCache[i] = NULL;
      }
    }
    
    imageCache.clear();
    
    if ( A        != NULL) delete[] A;
    if ( lambda_r != NULL) delete[] lambda_r;
    if ( Wk       != NULL) delete[] Wk; 
  }//~PngAnalysis()
};//PngAnalysis

#endif // defined(__PNG_ANALYSIS__)
