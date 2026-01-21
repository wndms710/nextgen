#ifndef _FWH_SURFACE_HPP_
#define _FWH_SURFACE_HPP_

#include "SimpleGeom.hpp"

class FwhSurface {
public: 

  int status;
  
  string name;
  int interval;
  vector<pair<string,CtiRegister::CtiData*> > var_vec;

  int ncv_pack; // the cvs we are responsible for packing
  int * cv_pack_array;
  
  // ----------------------------------------------------------
  // the following only defined on ranks that buffer and write:
  // Since this buffered data could be large, use one rank per 
  // node using the MpiStuff::rank_of_rank_internode[] array
  // ----------------------------------------------------------

  int nfa; // total number of fwh faces for this surface
  
  int ncv_unpack; // the number of cvs on the write_rank after all above are gathered
  int * count;
  int * disp;
  
  int ndpw; // number of data-probe(fa)-weight triples used to populate the final buffer
  int * cv_unpack_of_dpw; // data
  int * fa_fwh_of_dpw;    // probe
  double * wgt_of_dpw;    // weight

  double * fwh_buf; // the buffered fwh data with total size nfa*fwh_buf_array_count
  int fwh_buf_step;
  int fwh_buf_array_count; // 1 if p, 4 if p,u
  int fwh_buf_type; // e.g. FFW_DATA_D1D2, etc...
  
  FwhSurface() {
    
    //COUT1("[FWH_SURFACE] FwhSurface()");

    status = CTI_STATUS_ACTIVE;

    interval = -1;

    cv_pack_array = NULL;

    count = disp = NULL;

    cv_unpack_of_dpw = NULL;
    fa_fwh_of_dpw = NULL;
    wgt_of_dpw = NULL;
    
    fwh_buf = NULL;
    
  }

  ~FwhSurface() {

    //COUT1("[FWH_SURFACE] ~FwhSurface() name: " << name);
    
    DELETE(cv_pack_array);

    DELETE(count);
    DELETE(disp);

    DELETE(cv_unpack_of_dpw);
    DELETE(fa_fwh_of_dpw);
    DELETE(wgt_of_dpw);

    DELETE(fwh_buf);

  }
  
  void doProbe(const int step,const double time,const int write_rank) {
    
    COUT1("[FWH_SURFACE] " << name << " buffering fwh data");

    // look at what data is requested in var_vec. This should all
    // be CV DN or CV DN3 data (e.g. p and u)...
    
    vector<CtiRegister::CtiData*> data_vec(var_vec.size());
    int array_count = 0;
    for (int ii = 0; ii < var_vec.size(); ++ii) {
      CtiRegister::CtiData * data = var_vec[ii].second;
      if (data == NULL) data = CtiRegister::getUnregisteredCtiData(var_vec[ii].first);
      assert(data);
      assert(data->getUnindexedTopology() == CV_DATA);
      if (data->getType() == DN_DATA) {
	array_count += 1;
      }
      else if (data->getType() == DN3_DATA) {
	array_count += 3;
      }
      else {
	CERR("unsupported data type: " << data->getType());
      }
      // store the CtiData pointer...
      data_vec[ii] = data;
    }

    // we now have the array_count (1 per cv scalar, 3 per cv vector)...
    
    double * pack_buf = new double[ncv_pack*array_count];
    int array_offset = 0;
    for (int ii = 0; ii < data_vec.size(); ++ii) {
      if (data_vec[ii]->getType() == DN_DATA) {
	for (int icv_pack = 0; icv_pack < ncv_pack; ++icv_pack) {
	  const int icv = cv_pack_array[icv_pack];
	  pack_buf[icv_pack*array_count+array_offset] = data_vec[ii]->dn(icv);
	}
	array_offset += 1;
      }
      else if  (data_vec[ii]->getType() == DN3_DATA) {
	for (int icv_pack = 0; icv_pack < ncv_pack; ++icv_pack) {
	  const int icv = cv_pack_array[icv_pack];
	  pack_buf[icv_pack*array_count+array_offset]   = data_vec[ii]->dn3(icv,0);
	  pack_buf[icv_pack*array_count+array_offset+1] = data_vec[ii]->dn3(icv,1);
	  pack_buf[icv_pack*array_count+array_offset+2] = data_vec[ii]->dn3(icv,2);
	}
	array_offset += 3;
      }
      else {
	assert(0); 
      }
    }
    assert(array_offset == array_count);
    
    // unpack on the write_rank...

    double * unpack_buf = NULL;
    if (mpi_rank == write_rank) {
      unpack_buf = new double[ncv_unpack*array_count];
      // count and disp also get the factor of array_count...
      assert(count);
      assert(disp);
      FOR_RANK {
	count[rank] *= array_count;
	disp[rank] *= array_count;
      }
    }
    
    MPI_Gatherv(pack_buf,ncv_pack*array_count,MPI_DOUBLE,unpack_buf,count,disp,MPI_DOUBLE,write_rank,mpi_comm);
    delete[] pack_buf;

    if (mpi_rank == write_rank) {
      // reset count and disp...
      FOR_RANK {
	count[rank] /= array_count;
	disp[rank] /= array_count;
      }
      // the write_rank is managed externally. All we can check is that 
      // our fwh_buf is NULL. fwh_buf is the final buffer that will remain 
      // in memory until flush is called...
      assert(mpi_rank_shared == 0); // must be a shared memory rank 0 
      assert(fwh_buf == NULL);
      fwh_buf = new double[nfa*array_count];
      // zero, because we add into this buf with the weighted data...
      for (int ifa_fwh = 0; ifa_fwh < nfa; ++ifa_fwh) {
	for (int i = 0; i < array_count; ++i) {
	  fwh_buf[ifa_fwh*array_count+i] = 0.0;
	}
      }
      // now loop on unpack stuff and put into this including the weight...
      array_offset = 0;
      for (int ii = 0; ii < data_vec.size(); ++ii) {
	if (data_vec[ii]->getType() == DN_DATA) {
	  for (int idpw = 0; idpw < ndpw; ++idpw) {
	    const int icv_unpack = cv_unpack_of_dpw[idpw]; assert((icv_unpack >= 0)&&(icv_unpack < ncv_unpack));
	    const double wgt     = wgt_of_dpw[idpw];
	    const int ifa_fwh    = fa_fwh_of_dpw[idpw]; assert((ifa_fwh >= 0)&&(ifa_fwh < nfa));
	    fwh_buf[nfa*array_offset+ifa_fwh] += wgt*unpack_buf[icv_unpack*array_count+array_offset];
	  }
	  array_offset += 1;
	}
	else if  (data_vec[ii]->getType() == DN3_DATA) {
	  for (int idpw = 0; idpw < ndpw; ++idpw) {
	    const int icv_unpack = cv_unpack_of_dpw[idpw]; assert((icv_unpack >= 0)&&(icv_unpack < ncv_unpack));
	    const double wgt     = wgt_of_dpw[idpw];
	    const int ifa_fwh    = fa_fwh_of_dpw[idpw]; assert((ifa_fwh >= 0)&&(ifa_fwh < nfa));
	    fwh_buf[nfa*array_offset+ifa_fwh*3  ] += wgt*unpack_buf[icv_unpack*array_count+array_offset  ];
	    fwh_buf[nfa*array_offset+ifa_fwh*3+1] += wgt*unpack_buf[icv_unpack*array_count+array_offset+1];
	    fwh_buf[nfa*array_offset+ifa_fwh*3+2] += wgt*unpack_buf[icv_unpack*array_count+array_offset+2];
	  }
	  array_offset += 3;
	}
	else {
	  assert(0); 
	}
      }
      assert(array_offset == array_count);
      delete[] unpack_buf;
      // finally, on the write_rank, set the other data that remembers what type
      // of buffer to write...
      fwh_buf_array_count = array_count; // 1 if p, 4 if p,u
      fwh_buf_step = step;
      // set the type...
      if ((data_vec.size() == 1)&&(data_vec[0]->getType() == DN_DATA)) {
	fwh_buf_type = FFW_DATA_D1;
      }
      else if ((data_vec.size() == 2)&&(data_vec[0]->getType() == DN_DATA)&&(data_vec[1]->getType() == DN3_DATA)) {
	fwh_buf_type = FFW_DATA_D1D2;
      }
      else {
	// we can support other data patterns in the future
	fwh_buf_type = FFW_DATA_UNKNOWN;
      }
    }

  }
  
  void flush() {
    
    //COUT1("[FWH_SURFACE] " << name << " flush()");
    
    // only flush if we have a defined buffer...
    
    if (fwh_buf != NULL) {

      char filename[256];
      MiscUtils::buildIndexedFilename(filename,name.c_str(),fwh_buf_step,"data");

      // this is coming from different ranks, and can thus mangle the
      // output, so skip...
      //cout << "[FWH_SURFACE] writing " << filename << endl;  
      
      MiscUtils::mkdir_for_file(filename);
      const string tmp_filename = MiscUtils::makeTmpPrefix(filename);
      FILE * fp = fopen(tmp_filename.c_str(),"w");
	
      int ibuf[4];
      ibuf[0] = UGP_IO_MAGIC_NUMBER;
      ibuf[1] = UGP_IO_VERSION;
      ibuf[2] = fwh_buf_type;
      ibuf[3] = nfa;
      fwrite(ibuf,sizeof(int),4,fp);
      
      fwrite(fwh_buf,sizeof(double),nfa*fwh_buf_array_count,fp);
      
      fclose(fp);
      remove(filename);
      rename(tmp_filename.c_str(),filename);

      delete[] fwh_buf;
      fwh_buf = NULL;

    }

  }

};

#endif
