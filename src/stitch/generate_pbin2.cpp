

#include <iostream>
#include <cstdio>
#include <assert.h>
#include <climits>
#include <cmath>
#include <fstream>
#include "Adt.hpp"

#include "MiscUtils.hpp"
#include "MpiStuff.hpp"
#include "DistributedDataExchanger.hpp"
#include "CommonIo.hpp"
#include "Macros.hpp"

using namespace std;

//
// the format of the pbin file that is passed to stitch 
// is as follows: 
// 
//     int8 ibuf[4] = {POINTS_IO_MAGIC_NUMBER,   // to check endianness 
//                     POINTS_IO_VERSION,        // version number 
//                     np_global,                // number of points, 
//                     1};                       // indicates that the delta record is present
// 
//     double (*xp)[3];                          // point locations
//
//     double *delta;                            // measure of distance between neighboring points
//
// 
// the input file format is a space delimited list of x,y,z coordinates
// 
//       x_0    y_0    z_0 
//       x_1    y_1    z_1 
//             ... 
//       x_{n-1} y_{n-1} z_{n-1}

typedef long long int int8;

#define POINTS_IO_VERSION 2
#define POINTS_IO_MAGIC_NUMBER 1235813
#define HUGE_INT (1<<24)

int getNextToken(char * token,int& pos,const char * buf, const int max_pos) {
  int token_pos = 0;
  while (pos < max_pos) {
    char c = buf[pos++];
    //cout << "getNextToken pos: " << pos << " c: \"" << c << "\"" << endl;
    if ( (c == ' ') || (c == '\t') || (c == '\n') || (c == ',') || (c == '\0') || c == 13 || c == 31 || c == EOF) {
      if (token_pos > 0) {
        // we have a valid token, so quit...
        break;
      } 
      else {
        // nothing yet - keep going...
        continue;
      }
    } 
    else { 
      token[token_pos++] = c;
    }
  }
  // null terminate token...
  if (token_pos != 0) token[token_pos] = '\0';
  return token_pos; // returns size
}

double atod(char * token) {
  double d;
  sscanf(token, "%lf", &d);
  return (d);
}

void readPointsAscii(vector<double>& x_vec,const char* filename) {

  if ( mpi_rank == 0) 
    cout << " > reading points from " << filename << endl;

  // open...

  MPI_File fh;
  int ierr = MPI_File_open(mpi_comm,filename,MPI_MODE_RDONLY,MPI_INFO_NULL,&fh);
  if (ierr != 0) {
    CERR("cannot open restart file: " << filename);
  }

  // read in relevant chunk of file into "chunk",
  // which starts at location in the file globalstart
  // and has size mysize 

  // figure out who reads what
  MPI_Offset filesize;
  MPI_File_get_size(fh, &filesize);
  //filesize--;  /* get rid of text file eof */
  int8 mysize = filesize/mpi_size;
  MPI_Offset globalstart = mpi_rank * mysize;
  MPI_Offset globalend   = globalstart + mysize - 1;
  if (mpi_rank == mpi_size-1) globalend = filesize-1;
  mysize = globalend - globalstart + 1;

  // allocate memory 
  MPI_Barrier(mpi_comm);
  char *chunk = new char[mysize];

  // everyone reads in their part 
  readChunkedData<char>(fh, globalstart, chunk, mysize, false,mpi_comm);
  //MPI_File_read_at_all(fh, globalstart, chunk, mysize, MPI_CHAR, MPI_STATUS_IGNORE);
  if (mpi_rank == 0)
    cout << " > finished initial read via char striping" << endl;

  // count and index entries that you own

  int8 *ofonl = NULL;
  int8 nnl;
  for (int iter = 0; iter < 2; ++iter) {
    nnl = 0;
    for (int i = 0; i < mysize; ++i) {
      if (chunk[i] == '\n' || chunk[i] == EOF) {
        if (iter == 0) 
          nnl++;
        else
          ofonl[nnl++] = i+globalstart;
      }
    }
    if (iter == 0)
      ofonl = new int8[nnl];
  }
  delete[] chunk; chunk = NULL;

  // scan to determine which global entry/newlines you own...

  int8 nnl_global;
  MPI_Scan(&nnl,&nnl_global,1,MPI_INT8,MPI_SUM,mpi_comm);
  int8 my_inl_disp = nnl_global - nnl;
  MPI_Bcast(&nnl_global,1,MPI_INT8,mpi_size-1,mpi_comm); // everybody needs nl_cnt_global to set offset...

  // build striping
  int8 *nlora_striped = NULL;
  MiscUtils::calcThresholdDist(nlora_striped,nnl_global,mpi_size,DIST_THRESHOLD);


  // build dde 
  int8 *my_inl_global = new int8[nnl];
  const int8 nnl_striped = nlora_striped[mpi_rank+1]-nlora_striped[mpi_rank];
  int8 *ofonl_striped = new int8[nnl_striped+1]; ofonl_striped[0] = 0;

  int used_rank = 0;
  if (nnl_striped > 0) used_rank = 1;
  int nrank;
  MPI_Scan(&used_rank,&nrank,1,MPI_INT,MPI_SUM,mpi_comm);
  MPI_Bcast(&nrank,1,MPI_INT,mpi_size-1,mpi_comm); // everybody needs to know used striped ranks ...

  // send offsete of new-line to striped...
  {
    for (int i = 0; i < nnl; ++i) my_inl_global[i] = my_inl_disp+i; 
    DistributedDataExchanger dde(my_inl_global,nnl,nlora_striped);
    dde.push(ofonl_striped+1,ofonl);
  }
  delete[] my_inl_global;
  DELETE(ofonl);

  // send up your last new line...
  MPI_Request send_req,recv_req;
  if (mpi_rank < nrank-1) 
    MPI_Isend(&ofonl_striped[nnl_striped],1,MPI_INT8,mpi_rank+1,1234,MPI_COMM_WORLD,&send_req);
  if ((mpi_rank > 0) && (mpi_rank < nrank))
    MPI_Irecv(&ofonl_striped[0],1,MPI_INT8,mpi_rank-1,1234,MPI_COMM_WORLD,&recv_req);
  if (mpi_rank < nrank-1) 
    MPI_Wait(&send_req,MPI_STATUSES_IGNORE);
  if ((mpi_rank > 0) && (mpi_rank < nrank))
    MPI_Wait(&recv_req,MPI_STATUSES_IGNORE);
  if (mpi_rank == 0) ofonl_striped[0] = ofonl_striped[1];

  // now try reading again

  globalstart = ofonl_striped[0]+1;
  globalend = ofonl_striped[nlora_striped[mpi_rank+1]-nlora_striped[mpi_rank]];
  delete[] ofonl_striped;
  delete[] nlora_striped;
  //if (mpi_rank == mpi_size-1) globalend = filesize-1;
  mysize = globalend - globalstart + 1;
  if (mysize <= 0) {
    mysize = globalend = globalstart = 0;
  }
  assert(chunk == NULL); chunk = new char[mysize];

  // everyone reads in their part 
  readChunkedData<char>(fh, globalstart, chunk, mysize, false,mpi_comm);
  if (mpi_rank == 0)
    cout << " > finished final read using offsets" << endl;

  // get doubles from char tokens...
  char token[256];
  int pos = 0;
  assert(x_vec.size() == 0);
  while (pos < mysize) { 
    if (getNextToken(token,pos,chunk,mysize)) {
      x_vec.push_back(atod(token));
      //if (mpi_rank == 0) cout << mpi_rank << " " << x_vec.back() << endl;
    }
  }
  assert(x_vec.size()%3 == 0); 
  delete[] chunk;

  // and close...
  MPI_File_close(&fh);

}

void initBoundingBoxAndCounts(double bbox[6], double bbox_global[6], int8& np, 
                              int8& np_global, const vector<double>& xp) { 

  assert( xp.size()%3 == 0);
  np = xp.size()/3;

  //
  // this can be done much more intelligently, but for now, we are 
  // going to get a global length scale associated with bbox volumes.
  //
 
  for (int i = 0; i < 6; ++i) 
    bbox[i] = -HUGE_VAL;

  for (int ip = 0; ip < np; ++ip) {
    for (int i = 0; i < 3; ++i) { 
      bbox[i]     = fmax(bbox[i]  , -xp[3*ip+i]); 
      bbox[i+3]   = fmax(bbox[i+3], xp[3*ip+i]);
    }
  }

  MPI_Allreduce(bbox,bbox_global,6,MPI_DOUBLE,MPI_MAX,mpi_comm);

  for (int i = 0; i <3; ++i) { 
    bbox[i] *= -1.0;
    bbox_global[i] *= -1.0;
  }

  MPI_Allreduce(&np,&np_global,1,MPI_INT8,MPI_SUM,mpi_comm);


}


double getLengthScaleEstimateFromBbox(const double bbox_global[6], const int8 np_global) { 

  const double bbox_vol = (bbox_global[3] - bbox_global[0])*
                          (bbox_global[4] - bbox_global[1])*
                          (bbox_global[5] - bbox_global[2]);

  assert( bbox_vol > 0.0);

  const double L_estimate = pow(bbox_vol/double(np_global),1.0/3.0);

  if ( mpi_rank == 0 )
    cout << " > using estimate for local length scale: " << L_estimate << endl;

  return L_estimate;

}

double setDelta_(Adt<double>* xp_adt, const double point[3], const vector<double>& xp, const double L_est,
                 const int ip) {


  double L_guess = L_est;
  vector<int> nbr_vec;
  int iter           = 0;
  bool done          = false;
  const int max_iter = 100;
 
  double this_delta  = HUGE_VAL;

  while ( !done) { 
    
    // copy down the previous solutions 
    
    nbr_vec.clear();
    xp_adt->buildListForSphere(nbr_vec,point,L_guess);

    if ( nbr_vec.empty() || ( (nbr_vec.size() == 1)&&(nbr_vec[0] == ip) ) ) {

      L_guess *= 2.0;

    }
    else {

      for (int ii = 0, nlim = nbr_vec.size(); ii < nlim; ++ii) { 
        
        const int inbr = nbr_vec[ii];
        
        if ( inbr == ip) continue;

        double dist2   = 0.0;
        for (int i = 0 ; i <3 ; ++i) 
          dist2 += (point[i] - xp[3*inbr+i])*(point[i] - xp[3*inbr+i]);
        this_delta     = min(this_delta,sqrt(dist2));

      }

      done = true;

    } 

    ++iter;

    if ( iter == max_iter ) { 

      cout << " unable to find a single neighbor with L_guess = " << L_guess << " at xp: " 
           << COUT_VEC(point) << endl;
      cout.flush();

      assert(0);

    }

  }

  assert( this_delta != HUGE_VAL);
  assert( this_delta != 0.0);

  return this_delta;

}

void setDelta(double* &delta, const vector<double>& xp, const int np, const double L_est,  
              const double bbox[6], const double bbox_global[6],const int nchunks) { 

  if (mpi_rank == 0) cout << " > starting to set delta" << endl;

  assert( delta == NULL);  delta = new double[np];

  // in case the bounding boxes are all overlapping (worst case scenario), then you 
  // must be able to hold np + (np_global-np)/nchunks points in memory at a 
  // given time.

  int nn;
  if ((nchunks == 1)||(np < nchunks))
    nn = np;
  else 
    nn = np/(nchunks-1);

  double (*bbora)[6] = new double[mpi_size][6];
  MPI_Allgather(bbox,6,MPI_DOUBLE,bbora,6,MPI_DOUBLE,mpi_comm);

  // build the adts for the bboxs and the local set of points.. 
  
  Adt<double> *bb_adt = NULL, *xp_adt = NULL;
  
  { 
    double (*bbmin)[3] = new double[mpi_size][3];
    double (*bbmax)[3] = new double[mpi_size][3];

    for (int rank = 0; rank < mpi_size; ++rank) { 

      for (int i = 0; i < 3; ++i) { 
        bbmin[rank][i] = bbora[rank][i];
        bbmax[rank][i] = bbora[rank][i+3];
      }
    }

    delete[] bbora;

    bb_adt = new Adt<double>(mpi_size,bbmin,bbmax);

    delete[] bbmin;
    delete[] bbmax;

  }

  { 

    const double (*xp_tmp)[3] = (const double(*)[3])(&xp[0]);
    xp_adt              = new Adt<double>(np,xp_tmp,xp_tmp);

  }

  for (int ichunk = 0; ichunk < nchunks; ++ichunk ) { 

    const int ip_start  = ichunk * nn;
    const int ip_end    = min((ichunk+1)*nn,np);
    //FOR_RANK {
    //  if (rank == mpi_rank) 
    //    cout << ip_end-ip_start << endl;
    //  MPI_Barrier(mpi_comm);
    //}

    int* send_count     = new int[mpi_size];
    int* send_disp      = new int[mpi_size];

    double * send_buf   = NULL;
    int * ipois         = NULL;

    for (int rank = 0; rank < mpi_size; ++rank) 
      send_count[rank] = 0;

    for (int iter = 0; iter < 2; ++iter) { 

      for (int ip = ip_start ; ip < ip_end; ++ip) {
        
        // set our delta based on our own information .. 
        
        if ( iter == 0) 
          delta[ip] = setDelta_(xp_adt,&xp[3*ip],xp,L_est,ip); 

        // decide if this point needs to go anywhere else .. 

        vector<int> nbr_rank;
        bb_adt->buildListForPoint(nbr_rank,&xp[3*ip]);

        for (int ii  = 0, nlim = nbr_rank.size(); ii < nlim; ++ii) { 
        
          const int rank = nbr_rank[ii];
        
          if ( rank != mpi_rank ) { 

            if ( iter == 0) { 

              send_count[rank] += 3;

            } else { 

              ipois[send_disp[rank]/3]    = ip;
              send_buf[send_disp[rank]+0] = xp[3*ip+0];
              send_buf[send_disp[rank]+1] = xp[3*ip+1];
              send_buf[send_disp[rank]+2] = xp[3*ip+2];
              send_disp[rank] += 3;

            }

          }
        }  
      }
      
      send_disp[0] = 0;
      for (int rank = 1; rank < mpi_size; ++rank) 
        send_disp[rank] = send_disp[rank-1] + send_count[rank-1];
     
      if ( iter == 0 ) { 

        const int send_buf_size = send_disp[mpi_size-1] + send_count[mpi_size-1];
        assert( send_buf_size%3 == 0);
        send_buf = new double[send_buf_size];
        ipois    = new int[send_buf_size/3];

      }
    }

    // send over the count..

    int * recv_count = new int[mpi_size];
    MPI_Alltoall(send_count,1,MPI_INT,recv_count,1,MPI_INT,mpi_comm);

    int * recv_disp  = new int[mpi_size];
    recv_disp[0]     = 0;
    for (int rank = 1; rank < mpi_size; ++rank) 
      recv_disp[rank] = recv_disp[rank-1] + recv_count[rank-1];

    const int recv_buf_size = recv_disp[mpi_size-1] + recv_count[mpi_size-1];
    double * recv_buf       = new double[recv_buf_size];
    

    MPI_Alltoallv(send_buf,send_count,send_disp,MPI_DOUBLE,
                  recv_buf,recv_count,recv_disp,MPI_DOUBLE,mpi_comm);

    // leave the send_buf allocate to receive the data.. 

    assert( recv_buf_size%3 ==0);
    const int np_cand  = recv_buf_size/3;
    
    // we'll overwrite the recv_buf points inline with their respective deltas..

    for (int ip = 0; ip < np_cand; ++ip) { 

      recv_buf[ip] = setDelta_(xp_adt,&recv_buf[3*ip],xp,L_est,-1);

    }

    // compact the counts and disp 

    for (int rank = 0; rank < mpi_size; ++rank) { 

      send_count[rank] /= 3;
      send_disp[rank]  /= 3;
      recv_count[rank] /= 3;
      recv_disp[rank]  /= 3;

    }

    MPI_Alltoallv(recv_buf,recv_count,recv_disp,MPI_DOUBLE,
                  send_buf,send_count,send_disp,MPI_DOUBLE,mpi_comm);



    const int np_remote = send_disp[mpi_size-1] + send_count[mpi_size-1];

    for (int ii = 0; ii < np_remote; ++ii) { 

      const int ip = ipois[ii];
      delta[ip]    = min(delta[ip],send_buf[ii]);

    }

    if (mpi_rank == 0) 
      cout << " > finished chunk " << ichunk+1 << " of " << nchunks << endl;

    delete[] recv_buf;
    delete[] send_buf;
    delete[] recv_count;
    delete[] send_count;
    delete[] recv_disp;
    delete[] send_disp;
    delete[] ipois;

  }

  delete xp_adt;
  delete bb_adt;

}

void writePointsBinary(char* filename, double* delta, vector<double> xp, const int np, const int8 np_global) {

  if ( mpi_rank == 0) 
    cout << " > writing points to " << filename << endl;

  MPI_File_delete(filename,MPI_INFO_NULL);

  MPI_File fh;
  MPI_File_open(mpi_comm,filename,MPI_MODE_WRONLY|MPI_MODE_CREATE,MPI_INFO_NULL,&fh);

  if (mpi_rank == 0) {
    int8 ibuf[4] = { POINTS_IO_MAGIC_NUMBER, POINTS_IO_VERSION, np_global, 1 }; // last indicates delta 
    MPI_File_write_at(fh,0,ibuf,4,MPI_INT8,MPI_STATUS_IGNORE);
  }

  // points...

  int8 my_np = np;
  int8 my_disp;
  MPI_Scan(&my_np,&my_disp,1,MPI_INT8,MPI_SUM,mpi_comm);
  assert( (mpi_rank != mpi_size-1) || (my_disp == np_global) );
  MPI_Offset offset = (my_disp-my_np)*3*8 + 4*8;
  writeChunkedData<double>(fh,offset,&xp[0],np*3,mpi_comm);

  // delta...

  offset = np_global*3*8 + 4*8 + (my_disp-my_np)*8;
  writeChunkedData<double>(fh,offset,delta,np,mpi_comm);

  // and close...

  MPI_File_close(&fh);

}


int main(int argc, char* argv[]) { 

  double L_est = -1.0;
  int nchunks = 25;

  MPI_Init(&argc, &argv);
  initMpiStuff();

  if ( argc < 2 ) { 
    
    CERR(" Usage: ./generate_pbin.exe <file-containing-points> [<approx local length scale> <nchunks>]");
    
  }

  try { 

    vector<double> xp;
    double* delta   = NULL;
    
    // at this point, we read the ascii file in parallel, but the points 
    // are striped across the ranks.. 

    readPointsAscii(xp,argv[1]);

    double bbox[6], bbox_global[6];
    int8 np, np_global;

    initBoundingBoxAndCounts(bbox,bbox_global,np,np_global,xp);

    if ( argc == 2) { 

      L_est = getLengthScaleEstimateFromBbox(bbox_global,np_global);

    } else if ( argc == 3 ) { 

      L_est = atof(argv[2]);
      
      if ( L_est <= 0.0) { 

        CERR( " > must specify a postive length scale if providing one; got : " << L_est);

      }
    } else if ( argc == 4 ) {
      
      nchunks = atoi(argv[3]);

      if ( nchunks <= 0) { 

        CERR( " > must specify a postive number of chunks if providing one; got : " << nchunks);

      }

    }

    setDelta(delta,xp,np,L_est,bbox,bbox_global,nchunks);

    string name(argv[1]);
    size_t period_pos = name.find(".");
    if (period_pos != string::npos);
      name = name.substr(0,period_pos);
    char filename[128];
    sprintf(filename,"%s.pbin",name.c_str());
    writePointsBinary(filename,delta,xp,np,np_global);

    if ( delta != NULL) delete[] delta;

  } catch(...) { 

    CERR(" Unhandled exception -- see error message above.  Exiting... ");

  }

  return 0;

}
