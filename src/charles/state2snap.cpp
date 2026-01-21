#include "../core/common/Common.hpp"
#include "../core/common/ByteSwap.hpp"

int main(const int argc, const char* argv[]) { 
  using namespace std;

  if (argc != 6) { 
    
    cerr << "Usage: ./state2snap.exe <prefix> <nrank> <first step> <delta step> <final step> "<< endl;
    return 1;

  }

  try { 
    
    const string prefix = argv[1];
    const int nrank = atoi(argv[2]);
    const int istep = atoi(argv[3]);
    const int dstep = atoi(argv[4]);
    const int fstep = atoi(argv[5]);
    if ((fstep-istep)%dstep != 0) {
      cerr << "Final step: " << fstep << " cannot be reached from initial step: " << istep << " with interval: " << dstep << endl;
      return 1;
    }

    char filename[128];
    FILE* fp = NULL;
    vector<int8> icvGlobalVec;
    for (int rank = 0; rank < nrank; ++rank) {
      sprintf(filename,"%s.%s.%06d.state",prefix.c_str(),"id",rank);
      fp = fopen(filename,"rb");
      if (fp == NULL) {
        cerr << "Cannot open file: " << string(filename) << endl;
        return 1;
      }
      fseek(fp, 0L, SEEK_END);
      int8 file_size = ftell(fp);
      assert(file_size%(sizeof(int8)) == 0); 
      rewind(fp);
      const int ncv = file_size/sizeof(int8);
      int8* i8buf = new int8[ncv];
      fread(i8buf,sizeof(int8),ncv,fp);
      fclose(fp);
      for (int icv = 0; icv < ncv; ++icv) 
        icvGlobalVec.push_back(i8buf[icv]);
      delete[] i8buf;
    }
    int8 ncv_global = icvGlobalVec.size();
    cout << "number of cvs: " << ncv_global << endl;

    double *p0 = new double[ncv_global];
    double (*u0)[3] = new double[ncv_global][3];
    double *p = new double[ncv_global];
    double (*u)[3] = new double[ncv_global][3];
    for (int step = istep; step <= fstep; step += dstep) {
      cout << "Building sles from state, prefix: " << prefix << ", mpi_size: " << nrank << ", step: " << step << endl;
      ncv_global = 0;
      for (int rank = 0; rank < nrank; ++rank) {
        sprintf(filename,"%s.%s.%06d.%08d.state",prefix.c_str(),"pu",rank,step);
        fp = fopen(filename,"rb");
        if (fp == NULL) {
          cerr << "Cannot open file: " << string(filename) << endl;
          return 1;
        }
        fseek(fp, 0L, SEEK_END);
        int8 file_size = ftell(fp);
        assert(file_size%(4*sizeof(double)) == 0); 
        rewind(fp);
        const int ncv = file_size/(4*sizeof(double));
        fread(p0+ncv_global,sizeof(double),ncv,fp);
        fread((double*)u0+3*ncv_global,sizeof(double),3*ncv,fp);
        fclose(fp);
        ncv_global += ncv;
      }
      assert(ncv_global == icvGlobalVec.size());

      // sort...
      for (int ii = 0; ii < ncv_global; ++ii) {
        assert(icvGlobalVec[ii] < ncv_global);
        p[icvGlobalVec[ii]] = p0[ii];
        for (int i = 0; i < 3; ++i) 
          u[icvGlobalVec[ii]][i] = u0[ii][i];
      }

      sprintf(filename,"%s.%08d.sles",prefix.c_str(),step);
      FILE* fp = fopen(filename,"wb");
      int itmp[2] = {UGP_IO_MAGIC_NUMBER+1, 5}; // convention for data and not mesh data..
      fwrite(itmp,sizeof(int),2,fp); 
      
      // time...
      {
        Header header;
        header.id = UGP_IO_D0;
        sprintf(header.name,"%s","time");
        header.skip = sizeof(Header);
        header.rdata[0] = double(step);
        fwrite(&header,sizeof(Header),1,fp); 
      }

      // step...
      {
	Header header;
	header.id = UGP_IO_I0;
	sprintf(header.name,"%s","step");
	header.skip = sizeof(Header);
	header.idata[0] = step;
        fwrite(&header,sizeof(Header),1,fp); 
      }


      // p...
      {
	Header header;
	header.id = UGP_IO_CV_D1;
	sprintf(header.name,"%s","p");
	header.skip = sizeof(Header) + ncv_global*sizeof(double);
	ByteSwap::setLswMswPairForInt8(header.idata+0,ncv_global);
        fwrite(&header,sizeof(Header),1,fp); 
        fwrite(p,sizeof(double),ncv_global,fp); 
      }

      // u...
      {
	Header header;
	header.id = UGP_IO_CV_D2;
	sprintf(header.name,"%s","u");
	header.skip = sizeof(Header) + 3*ncv_global*sizeof(double);
	ByteSwap::setLswMswPairForInt8(header.idata+0,ncv_global);
	header.idata[2] = 3;
        fwrite(&header,sizeof(Header),1,fp); 
        fwrite((double*)u,sizeof(double),3*ncv_global,fp);
      }

      // eof...
      {
        Header header;
        header.id = UGP_IO_EOF;
        sprintf(header.name,"EOF");
        header.skip = sizeof(header);
        fwrite(&header,sizeof(Header),1,fp); 
      }
      fclose(fp);
    }
    delete[] p0;
    delete[] p;
    delete[] u0;
    delete[] u;

  } catch(...) { 

    cerr << "Unhandled exception -- see error message above.  Exiting... " << endl;
    return 1;

  }

  return 0;
}
