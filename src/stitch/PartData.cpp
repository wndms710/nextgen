#include "PartData.hpp"
#include "SplineStuff.hpp"
#include "FluentReader.hpp"
#include "CuttableVoronoiData.hpp"

namespace PartData {

  // TODO: put this in MpiStuff...

  void MPI_Gather_set(set<int>& intSet,const int rank0,const MPI_Comm& mpi_comm) {
    int mpi_rank; MPI_Comm_rank(mpi_comm,&mpi_rank);
    int mpi_size; MPI_Comm_size(mpi_comm,&mpi_size);
    int n = intSet.size();
    int * nora = NULL;
    if (mpi_rank == rank0) {
      nora = new int[mpi_size];
      n = 0; // do not include anything from rank0. We already have that info
    }
    MPI_Gather(&n,1,MPI_INT,nora,1,MPI_INT,rank0,mpi_comm);
    int * ibuf = new int[n];
    if (mpi_rank != rank0) {
      int ii = 0;
      for (set<int>::iterator it = intSet.begin(); it != intSet.end(); ++it)
        ibuf[ii++] = *it;
      assert(ii == n);
    }
    int * disp = NULL;
    int * ibuf_all = NULL;
    if (mpi_rank == rank0) {
      disp = new int[mpi_size];
      disp[0] = 0;
      for (int rank = 1; rank < mpi_size; ++rank)
        disp[rank] = disp[rank-1] + nora[rank-1];
      const int n_all = disp[mpi_size-1] + nora[mpi_size-1];
      ibuf_all = new int[n_all];
    }
    MPI_Gatherv(ibuf,n,MPI_INT,ibuf_all,nora,disp,MPI_INT,rank0,mpi_comm);
    delete[] ibuf;
    if (mpi_rank == rank0) {
      const int n_all = disp[mpi_size-1] + nora[mpi_size-1];
      delete[] nora;
      delete[] disp;
      for (int ii = 0; ii < n_all; ++ii)
        intSet.insert(ibuf_all[ii]);
      delete[] ibuf_all;
    }
  }

  void MPI_Gather_set(set<pair<int,int> >& intSet,const int rank0,const MPI_Comm& mpi_comm) {
    int mpi_rank; MPI_Comm_rank(mpi_comm,&mpi_rank);
    int mpi_size; MPI_Comm_size(mpi_comm,&mpi_size);
    int n = intSet.size()*2;
    int * nora = NULL;
    if (mpi_rank == rank0) {
      nora = new int[mpi_size];
      n = 0; // do not include anything from rank0. We already have that info
    }
    MPI_Gather(&n,1,MPI_INT,nora,1,MPI_INT,rank0,mpi_comm);
    int * ibuf = new int[n];
    if (mpi_rank != rank0) {
      int ii = 0;
      for (set<pair<int,int> >::iterator it = intSet.begin(); it != intSet.end(); ++it) {
        ibuf[ii] = it->first;
        ibuf[ii+1] = it->second;
        ii += 2;
      }
      assert(ii == n);
    }
    int * disp = NULL;
    int * ibuf_all = NULL;
    if (mpi_rank == rank0) {
      disp = new int[mpi_size];
      disp[0] = 0;
      for (int rank = 1; rank < mpi_size; ++rank)
        disp[rank] = disp[rank-1] + nora[rank-1];
      const int n_all = disp[mpi_size-1] + nora[mpi_size-1];
      ibuf_all = new int[n_all];
    }
    MPI_Gatherv(ibuf,n,MPI_INT,ibuf_all,nora,disp,MPI_INT,rank0,mpi_comm);
    delete[] ibuf;
    if (mpi_rank == rank0) {
      const int n_all = disp[mpi_size-1] + nora[mpi_size-1];
      delete[] nora;
      delete[] disp;
      for (int ii = 0; ii < n_all; ii += 2)
        intSet.insert(pair<int,int>(ibuf_all[ii],ibuf_all[ii+1]));
      delete[] ibuf_all;
    }
  }

  // part containment is done using an integer space that requires
  // these double-to-integer conversion params...
  bool b_dtoi = false;
  double dtoi_x0[3];
  double dtoi_delta;

  void setDtoiStuff() {
    assert(!b_dtoi);
    // to convert doubles to ints (dtoi), use: dtoi_x0, dtoi_delta..
    double bbmin[3],bbmax[3];
    getBbox(bbmin,bbmax);
    FOR_I3 dtoi_x0[i] = 0.5*(bbmin[i]+bbmax[i]);
    dtoi_delta = max(bbmax[0]-bbmin[0],max(bbmax[1]-bbmin[1],bbmax[2]-bbmin[2]));
    b_dtoi = true;
  }

  // if you hit the assert(fabs(dx/delta) < 2.0), it means the delta is not
  // big enough for your domain. Since points are even, surface tris are odd,
  // if you hit the odd, it means you have large tris relative to the current
  // bbox, and need to enlarge. if you hit even problems, it must mean you are checking
  // points that you could eliminate with a surface0 bbox check...

  // currently the moving solver manages its own bbox, frozen at the initial
  // state, then sets the dtoi stuff on its own...

  // note: we had to increase the tolerance here to support periodic copies...

  int getOddInt(const double dx,const double delta) {
    assert(fabs(dx/delta) < 2.0); // see note above
    int value = (int)floor(dx/delta*1000000000); // use 1 billion here?
    if (value%2 == 0)
      value += 1;
    return value;
  }

  int getEvenInt(const double dx,const double delta) {
    assert(fabs(dx/delta) < 2.0); // see note above
    int value = (int)floor(dx/delta*1000000000); // use 1 billion here?
    if (value%2)
      value += 1;
    return value;
  }

  vector<Part*> partVec;
  vector<StZone> zoneVec;
  map<const string,int> zoneMap;
  vector<pair<int,int> > partSurfaceZoneVec;
  Points * hcpPts = NULL;
  bool b_vp = false;

  int getFlattenedSp(const int ipart,const int isp) {
    assert((ipart >= 0)&&(ipart < partVec.size()));
    assert(isp >= 0);
    int sp_disp = 0;
    for (int ii = 0; ii < ipart; ++ii) {
      if (partVec[ii]->surface)
        sp_disp += partVec[ii]->surface->nsp;
    }
    return sp_disp+isp;
  }
  int getFlattenedSt(const int ipart,const int ist) {
    assert((ipart >= 0)&&(ipart < partVec.size()));
    assert(ist >= 0);
    int st_disp = 0;
    for (int ii = 0; ii < ipart; ++ii) {
      if (partVec[ii]->surface)
        st_disp += partVec[ii]->surface->nst;
    }
    return st_disp+ist;
  }
  void getPartStFromFlattenedSt(int& ipart,int& ist,const int flattened_ist) {
    ist = flattened_ist;
    for (ipart = 0; ipart < partVec.size(); ++ipart) {
      if (partVec[ipart]->surface) {
        if (ist >= partVec[ipart]->surface->nst)
          ist -= partVec[ipart]->surface->nst;
        else
          break;
      }
    }
    assert(flattened_ist == getFlattenedSt(ipart,ist));
  }
  void getPartSpFromFlattenedSp(int& ipart,int& isp,const int flattened_isp) {
    isp = flattened_isp;
    for (ipart = 0; ipart < partVec.size(); ++ipart) {
      if (partVec[ipart]->surface) {
        if (isp >= partVec[ipart]->surface->nsp)
          isp -= partVec[ipart]->surface->nsp;
        else
          break;
      }
    }
    assert(flattened_isp == getFlattenedSp(ipart,isp));
  }
  int getFlattenedNsp() {
    int nsp = 0;
    for (int ipart = 0; ipart < partVec.size(); ++ipart) {
      if (partVec[ipart]->surface)
        nsp += partVec[ipart]->surface->nsp;
    }
    return nsp;
  }
  int getFlattenedNst() {
    int nst = 0;
    for (int ipart = 0; ipart < partVec.size(); ++ipart) {
      if (partVec[ipart]->surface)
        nst += partVec[ipart]->surface->nst;
    }
    return nst;
  }
  double getFlattenedVolume() {
    double vol = 0.0;
    for (int ipart = 0; ipart < partVec.size(); ++ipart) {
      if (partVec[ipart]->surface)
        vol -= partVec[ipart]->surface->getVolume(); // NOTE that getVolume() is negative w/ outward normals
    }
    return vol;
  }
  void getFlattenedCentroid(double xc[3]) {
    double vol = 0.0;
    FOR_I3 xc[i] = 0.0;
    for (int ipart = 0; ipart < partVec.size(); ++ipart) {
      if (partVec[ipart]->surface) {
        vol -= partVec[ipart]->surface->getVolume(); // NOTE that getVolume() is negative w/ outward normals
        double tmp[3]; partVec[ipart]->surface->getCentroid(tmp);
        FOR_I3 xc[i] += partVec[ipart]->surface->getVolume()*tmp[i];
      }
    }
    assert(vol > 0.0);
    FOR_I3 xc[i] /= vol;
  }

  bool b_bbminmax = false;
  double bbminmax[6]; // stores bbox as xmin,ymin,zmin,-xmax,-ymax,-zmax

  bool testPointValidForTransform(const double xp[3],const int bits) {
    
    // this is a non-collective call, so make sure we have the b_bbminmax already defined
    // from earlier processing...
    assert(b_bbminmax);
    double xc[3];
    getBoundingBoxCenter(xc); // TODO: this may not be precise enough for certain orientations or thin periodicity
      
    // this routine was developed to eliminate the periodic transform of certain surface tris 
    // or points from consideration. The point we are considering transforming to include amongst our surfaces or
    // nbrs is in xp (this is its original, untransformed location). xc contains the centroid of the 
    // current region. We only want to consider points where the application of the transformation 
    // will switch the orientation of the point wrt the centroid. 
    //
    // This routine was constructed to handle the pathalogical case that emerges with 180 deg
    // periodicity, and both the transform and the reverse transform move the point to the 
    // same location. With the current approach, we should never double transform points, and
    // this was happening before...

    for (int ipt = 0, npt = PeriodicData::periodicTransformVec.size(); ipt < npt; ++ipt) {

      //cout << "ipt: " << ipt << " PeriodicData::periodicTransformVec[ipt].dump(): ";
      //PeriodicData::periodicTransformVec[ipt].dump();

      if (bits&(3<<(ipt*2))) {
	double data[3];
	PeriodicData::periodicTransformVec[ipt].getData(data);
	//cout << " >> DATA: " << COUT_VEC(data) << endl;
	if (bits&(1<<(ipt*2))) {
	  assert(!(bits&(2<<(ipt*2))));
	  switch (PeriodicData::periodicTransformVec[ipt].getKind()) {
	  case PERIODIC_TRANSFORM_CART:
	    // note: tri centroids could not be eliminated in the vki case using this
	    // approach with either sign...
	    //// check if the translation dotted on the dx is positive...
	    //if ((xp[0]-xc[0])*data[0]+(xp[1]-xc[1])*data[1]+(xp[2]-xc[2])*data[2] > 0.0)
	    //  return true;
	    break;
	  case PERIODIC_TRANSFORM_CYL_X:
	    // check if the x-component of xp x xc is signed the same way as the angle stored in data[2]...
	    if ((fabs(data[2]) == 180.0)&&((xp[1]*xc[2]-xp[2]*xc[1])*data[2] < 0.0))
	      return false;
	    break;
	  case PERIODIC_TRANSFORM_CYL_Y:
	    // check if the y-component of xp x xc is signed the same way as the angle stored in data[2]...
	    if ((fabs(data[2]) == 180.0)&&((xp[2]*xc[0]-xp[0]*xc[2])*data[2] < 0.0))
	      return false;
	    break;
	  case PERIODIC_TRANSFORM_CYL_Z:
	    // check if the z-component of xp x xc is signed the same way as the angle stored in data[2]...
	    if ((fabs(data[2]) == 180.0)&&((xp[0]*xc[1]-xp[1]*xc[0])*data[2] < 0.0))
	      return false;
	    break;
	  default:
	    assert(0);
	  }
	}
	else {
	  assert(bits&(2<<(ipt*2)));
	  switch (PeriodicData::periodicTransformVec[ipt].getKind()) {
	  case PERIODIC_TRANSFORM_CART:
	    // see note above...
	    //// check if the translation dotted on the dx is negative...
	    //if ((xp[0]-xc[0])*data[0]+(xp[1]-xc[1])*data[1]+(xp[2]-xc[2])*data[2] < 0.0)
	    //  return true;
	    break;
	  case PERIODIC_TRANSFORM_CYL_X:
	    // check if the x-component of xp x xc is signed the same way as the angle stored in data[2]...
	    if ((fabs(data[2]) == 180.0)&&((xp[1]*xc[2]-xp[2]*xc[1])*data[2] > 0.0))
	      return false;
	    break;
	  case PERIODIC_TRANSFORM_CYL_Y:
	    // check if the y-component of xp x xc is signed the same way as the angle stored in data[2]...
	    if ((fabs(data[2]) == 180.0)&&((xp[2]*xc[0]-xp[0]*xc[2])*data[2] > 0.0))
	      return false;
	    break;
	  case PERIODIC_TRANSFORM_CYL_Z:
	    // check if the z-component of xp x xc is signed the same way as the angle stored in data[2]...
	    if ((fabs(data[2]) == 180.0)&&((xp[0]*xc[1]-xp[1]*xc[0])*data[2] > 0.0))
	      return false;
	    break;
	  default:
	    assert(0);
	  }
	}
      }
    }
      
    // if we were not able to determine invalid, then must be valid...
    return true;
      
  }

  int Part::readBinary(const string& filename) {

    if (mpi_rank == 0) cout << "Part::readBinary: " << filename << endl;

    // parallel read into shm...
    const int zone_name_len_max = 128;

    int ibuf[5];
    int8 i8buf[4];
    int size_buf[12] = { -1, -1, -1, -1, -1,
                         -1, -1, -1, -1, -1,
                         -1, -1 };
    int my_ierr = 0;

    int nzones;
    char * char_buf = NULL;
    double (*xsp_split)[3] = NULL;
    int (*spost_split)[3] = NULL;
    int *znost_split = NULL;

    int ff_nzones;
    char * ff_char_buf = NULL;
    double (*ff_xsp_split)[3] = NULL;
    int (*ff_spost_split)[3] = NULL;
    int *ff_znost_split = NULL;
    double *ff_dxost_split = NULL;

    int my_np;
    double (*xp_split)[3] = NULL;
    double *deltap_split = NULL;

    assert(surface == NULL);
    assert(znosz == NULL);
    assert(ff_surface == NULL);
    assert(ff_surface_dxost == NULL);

    if (mpi_rank_shared == 0) {

      // ---------------------------------------------------------------
      // only one rank per node participates in the parallel read...
      // ---------------------------------------------------------------

      MPI_File fh;
      char dummy[128];
      sprintf(dummy,"%s",filename.c_str());
      my_ierr = MPI_File_open(mpi_comm_internode,dummy,MPI_MODE_RDONLY,MPI_INFO_NULL,&fh);
      if (my_ierr != 0) {
        if (mpi_rank == 0) cout << "Error: cannot open binary file \"" << filename << "\"" << endl;
        my_ierr = -1; // force -1 to ensure the MPI_MIN below catches this._
      }
      else {

        if (mpi_rank == 0) {
          MPI_File_read(fh,ibuf,5,MPI_INT,MPI_STATUS_IGNORE);
        }
        MPI_Bcast(ibuf,5,MPI_INT,0,mpi_comm_internode);

        // version check...
        // ibuf[0] == PART_MAGIC...
        if (mpi_rank == 0) cout << " > version: " << ibuf[1] << endl;
        if (ibuf[1] != 1) {
          if (mpi_rank == 0) cout << "Error: file does not start as expected \"" << filename << "\"" << endl;
          MPI_File_close(&fh);
          my_ierr = -1; // force -1 to ensure the MPI_MIN below catches this.
        }
        else {

          bool b_surf    = (ibuf[2] == 1);
          bool b_ff_surf = (ibuf[3] == 1);
          bool b_pts     = (ibuf[4] == 1);
          if (mpi_rank == 0) cout << " > b_surf: " << b_surf << ", b_ff_surf: " << b_ff_surf << ", b_pts: " << b_pts << endl;

          MPI_Offset offset = int_size*5;

          if (b_surf) {

            // surface is present...
            surface = new SurfaceShm();

            if (mpi_rank == 0) {

              MPI_File_read_at(fh,offset,ibuf,3,MPI_INT,MPI_STATUS_IGNORE);
              offset += int_size*3;

              cout << " > surface nzones: " << ibuf[0] << " nsp: " << ibuf[1] << " nst: " << ibuf[2] << endl;

              nzones = ibuf[0];
              surface->zoneVec.resize(nzones); // mpi_rank 0 only for now
              char_buf = new char[zone_name_len_max*nzones];
              int length;
              for (int izone = 0; izone < nzones; ++izone) {
                MPI_File_read_at(fh,offset,&length,1,MPI_INT,MPI_STATUS_IGNORE);
                assert((length > 0)&&(length+1 <= zone_name_len_max));
                MPI_File_read_at(fh,offset+int_size,char_buf+zone_name_len_max*izone,length,MPI_CHAR,MPI_STATUS_IGNORE);
                offset += int_size + length;
                char_buf[zone_name_len_max*izone+length] = '\0';
                surface->zoneVec[izone].setName(char_buf+zone_name_len_max*izone);
                cout << " > zone: \"" << surface->zoneVec[izone].getName() << "\"" << endl;
              }

              // use a single i8buf to share the counts and offsets...

              i8buf[0] = ibuf[0]; // nzone
              i8buf[1] = ibuf[1]; // nsp
              i8buf[2] = ibuf[2]; // nst
              i8buf[3] = offset;

            }

            // the remaining records are read in striped, then allgather'd into shm later...

            MPI_Bcast(i8buf,4,MPI_INT8,0,mpi_comm_internode);
            nzones       = i8buf[0];
            surface->nsp = i8buf[1];
            surface->nst = i8buf[2];
            offset       = i8buf[3];

            size_buf[0] = nzones;
            size_buf[1] = surface->nsp;
            size_buf[2] = surface->nst;

            int nsp_split = surface->nsp/mpi_size_internode;
            if (surface->nsp%mpi_size_internode != 0) nsp_split += 1;
            size_buf[3] = nsp_split;

            int isp0 = min(surface->nsp,nsp_split*mpi_rank_internode);
            int isp1 = min(surface->nsp,nsp_split*(mpi_rank_internode+1));
            int my_nsp = isp1 - isp0;
            xsp_split = new double[nsp_split][3];
            MPI_File_read_at_all(fh,offset+double_size*nsp_split*mpi_rank_internode*3,
                                 (double*)xsp_split,my_nsp*3,MPI_DOUBLE,MPI_STATUS_IGNORE);
            offset += double_size*surface->nsp*3;

            {
              double my_buf[6] = { HUGE_VAL, HUGE_VAL, HUGE_VAL, HUGE_VAL, HUGE_VAL, HUGE_VAL };
              for (int isp = 0; isp < my_nsp; ++isp) {
                FOR_I3 my_buf[i]   = min(my_buf[i],xsp_split[isp][i]);
                FOR_I3 my_buf[3+i] = min(my_buf[3+i],-xsp_split[isp][i]);
              }
              double buf[6];
              MPI_Reduce(my_buf,buf,6,MPI_DOUBLE,MPI_MIN,0,mpi_comm_internode);
              if (mpi_rank == 0) {
                cout << " > xsp range: 0: " << buf[0] << " " << -buf[3] <<
                  ", 1: " << buf[1] << " " << -buf[4] <<
                  ", 2: " << buf[2] << " " << -buf[5] << endl;
              }
            }

            int nst_split = surface->nst/mpi_size_internode;
            if (surface->nst%mpi_size_internode != 0) nst_split += 1;
            size_buf[4] = nst_split;

            int ist0 = min(surface->nst,nst_split*mpi_rank_internode);
            int ist1 = min(surface->nst,nst_split*(mpi_rank_internode+1));
            int my_nst = ist1 - ist0;
            spost_split = new int[nst_split][3];
            MPI_File_read_at_all(fh,offset+int_size*nst_split*mpi_rank_internode*3,
                                 (int*)spost_split,my_nst*3,MPI_INT,MPI_STATUS_IGNORE);
            offset += int_size*surface->nst*3;

            znost_split = new int[nst_split];
            MPI_File_read_at_all(fh,offset+int_size*nst_split*mpi_rank_internode,
                                 znost_split,my_nst,MPI_INT,MPI_STATUS_IGNORE);
            offset += int_size*surface->nst;

          } // if (b_surf)

          if (b_ff_surf) {

            // ff_surface is present...
            ff_surface = new SurfaceShm();

            if (mpi_rank == 0) {

              MPI_File_read_at(fh,offset,ibuf,3,MPI_INT,MPI_STATUS_IGNORE);
              offset += int_size*3;

              cout << " > ff_surface nzones: " << ibuf[0] << " nsp: " << ibuf[1] << " nst: " << ibuf[2] << endl;

              ff_nzones = ibuf[0];
              ff_surface->zoneVec.resize(ff_nzones); // mpi_rank 0 only for now
              ff_char_buf = new char[zone_name_len_max*ff_nzones];
              int length;
              for (int izone = 0; izone < ff_nzones; ++izone) {
                MPI_File_read_at(fh,offset,&length,1,MPI_INT,MPI_STATUS_IGNORE);
                assert((length > 0)&&(length+1 <= zone_name_len_max));
                MPI_File_read_at(fh,offset+int_size,ff_char_buf+zone_name_len_max*izone,length,MPI_CHAR,MPI_STATUS_IGNORE);
                offset += int_size + length;
                ff_char_buf[zone_name_len_max*izone+length] = '\0';
                ff_surface->zoneVec[izone].setName(ff_char_buf+zone_name_len_max*izone);
                cout << " > zone: \"" << ff_surface->zoneVec[izone].getName() << "\"" << endl;
              }

              // use a single i8buf to share the counts and offsets...

              i8buf[0] = ibuf[0]; // nzone
              i8buf[1] = ibuf[1]; // nsp
              i8buf[2] = ibuf[2]; // nst
              i8buf[3] = offset;

            }

            // the remaining records are read in striped, then allgather'd into shm later...

            MPI_Bcast(i8buf,4,MPI_INT8,0,mpi_comm_internode);
            ff_nzones       = i8buf[0];
            ff_surface->nsp = i8buf[1];
            ff_surface->nst = i8buf[2];
            offset          = i8buf[3];

            size_buf[5] = ff_nzones;
            size_buf[6] = ff_surface->nsp;
            size_buf[7] = ff_surface->nst;

            int nsp_split = ff_surface->nsp/mpi_size_internode;
            if (ff_surface->nsp%mpi_size_internode != 0) nsp_split += 1;
            size_buf[8] = nsp_split;

            int isp0 = min(ff_surface->nsp,nsp_split*mpi_rank_internode);
            int isp1 = min(ff_surface->nsp,nsp_split*(mpi_rank_internode+1));
            int my_nsp = isp1 - isp0;
            ff_xsp_split = new double[nsp_split][3];
            MPI_File_read_at_all(fh,offset+double_size*nsp_split*mpi_rank_internode*3,
                                 (double*)ff_xsp_split,my_nsp*3,MPI_DOUBLE,MPI_STATUS_IGNORE);
            offset += double_size*ff_surface->nsp*3;

            int nst_split = ff_surface->nst/mpi_size_internode;
            if (ff_surface->nst%mpi_size_internode != 0) nst_split += 1;
            size_buf[9] = nst_split;

            int ist0 = min(ff_surface->nst,nst_split*mpi_rank_internode);
            int ist1 = min(ff_surface->nst,nst_split*(mpi_rank_internode+1));
            int my_nst = ist1 - ist0;
            ff_spost_split = new int[nst_split][3];
            MPI_File_read_at_all(fh,offset+int_size*nst_split*mpi_rank_internode*3,
                                 (int*)ff_spost_split,my_nst*3,MPI_INT,MPI_STATUS_IGNORE);
            offset += int_size*ff_surface->nst*3;

            ff_znost_split = new int[nst_split];
            MPI_File_read_at_all(fh,offset+int_size*nst_split*mpi_rank_internode,
                                 ff_znost_split,my_nst,MPI_INT,MPI_STATUS_IGNORE);
            offset += int_size*ff_surface->nst;

            ff_dxost_split = new double[nst_split];
            MPI_File_read_at_all(fh,offset+double_size*nst_split*mpi_rank_internode,
                                 ff_dxost_split,my_nst,MPI_DOUBLE,MPI_STATUS_IGNORE);
            offset += double_size*ff_surface->nst;

          } // if (b_ff_surf)

          if (b_pts) {

            // pts are present...
            pts = new Points();

            int np_global;
            if (mpi_rank == 0) {
              MPI_File_read_at(fh,offset,&np_global,1,MPI_INT,MPI_STATUS_IGNORE);
              cout << " > pts np_global: " << np_global << endl;
            }
            offset += int_size;

            MPI_Bcast(&np_global,1,MPI_INT,0,mpi_comm_internode);
            pts->np_global = np_global;

            // read in the particle locations striped...
            int np_split = np_global/mpi_size_internode;
            if (np_global%mpi_size_internode != 0) np_split += 1;

            int ip0 = min(np_global,np_split*mpi_rank_internode);
            int ip1 = min(np_global,np_split*(mpi_rank_internode+1));
            my_np = ip1-ip0;

            // while we only need to read my_np, we need to allocate the array to be
            // large enough to scatter (rather than scatterv) to all shared nodes...
            int my_np_alloc = my_np;
            if (my_np%mpi_size_shared)
              my_np_alloc += mpi_size_shared - my_np%mpi_size_shared;

            xp_split = new double[my_np_alloc][3];
            MPI_File_read_at_all(fh,offset+double_size*np_split*mpi_rank_internode*3,
                                 (double*)xp_split,my_np*3,MPI_DOUBLE,MPI_STATUS_IGNORE);
            offset += double_size*np_global*3;

            deltap_split = new double[my_np_alloc];
            MPI_File_read_at_all(fh,offset+double_size*np_split*mpi_rank_internode,
                                 deltap_split,my_np,MPI_DOUBLE,MPI_STATUS_IGNORE);
            offset += double_size*np_global;

            {
              double my_buf[8] = { HUGE_VAL, HUGE_VAL, HUGE_VAL,
                                   HUGE_VAL, HUGE_VAL, HUGE_VAL,
                                   HUGE_VAL, HUGE_VAL };
              for (int ip = 0; ip < my_np; ++ip) {
                FOR_I3 my_buf[i]   = min(my_buf[i],xp_split[ip][i]);
                FOR_I3 my_buf[3+i] = min(my_buf[3+i],-xp_split[ip][i]);
                my_buf[6]          = min(my_buf[6],deltap_split[ip]);
                my_buf[7]          = min(my_buf[7],-deltap_split[ip]);
              }
              double buf[8];
              MPI_Reduce(my_buf,buf,8,MPI_DOUBLE,MPI_MIN,0,mpi_comm_internode);
              if (mpi_rank == 0) {
                cout << " > xp range X: " << buf[0] << " " << -buf[3] <<
                  ", Y: " << buf[1] << " " << -buf[4] <<
                  ", Z: " << buf[2] << " " << -buf[5] << endl;
                cout << " > delta range: " << buf[6] << " " << -buf[7] << endl;
              }
            }

            size_buf[10] = np_global;
            size_buf[11] = my_np;

          }

          MPI_File_close(&fh);

        }

      }

    }

    // error checking before going any further...
    int ierr;
    MPI_Allreduce(&my_ierr,&ierr,1,MPI_INT,MPI_MIN,mpi_comm);
    if (ierr != 0)
      return -1;

    // recall...
    //size_buf[0] = nzones;
    //size_buf[1] = surface->nsp;
    //size_buf[2] = surface->nst;
    //size_buf[3] = nsp_split;
    //size_buf[4] = nst_split;
    //size_buf[5] = ff_nzones;
    //size_buf[6] = ff_surface->nsp;
    //size_buf[7] = ff_surface->nst;
    //size_buf[8] = nsp_split;
    //size_buf[9] = nst_split;
    //size_buf[10] = np_global;
    //size_buf[11] = my_np;

    // everyone needs the sizes in size_buf (currently set on only the shared rank 0)...
    MPI_Bcast(size_buf,12,MPI_INT,0,mpi_comm_shared);
    if (mpi_rank_shared != 0) {
      // surface is the first 5...
      if (size_buf[0] != -1) {
        assert(surface == NULL);
        surface = new SurfaceShm();
        nzones = size_buf[0];
        surface->nsp = size_buf[1];
        surface->nst = size_buf[2];
      }
      // then ff_surface...
      if (size_buf[5] != -1) {
        assert(ff_surface == NULL);
        ff_surface = new SurfaceShm();
        ff_nzones = size_buf[5];
        ff_surface->nsp = size_buf[6];
        ff_surface->nst = size_buf[7];
      }
      // finally pts...
      if (size_buf[10] != -1) {
        assert(pts == NULL);
        pts = new Points();
        pts->np_global = size_buf[10];
        my_np = size_buf[11];
      }
    }

    if (surface) {

      // zone names get bcast to EVERYONE...
      if (mpi_rank != 0) {
        assert(char_buf == NULL);
        char_buf = new char[zone_name_len_max*nzones];
      }
      MPI_Bcast(char_buf,zone_name_len_max*nzones,MPI_CHAR,0,mpi_comm);
      if (mpi_rank != 0) {
        surface->zoneVec.resize(nzones);
        for (int izone = 0; izone < nzones; ++izone) {
          surface->zoneVec[izone].setName(char_buf+zone_name_len_max*izone);
        }
      }
      delete[] char_buf;

      // allocate the shared memory for xsp AND spost AND znost

      if (mpi_rank == 0) cout << " > allgather on surface xsp..." << endl;

      // allocate the shared memory for xsp AND spost AND znost

      surface->nsp_max = size_buf[3]*mpi_size_internode; assert(surface->nsp_max >= surface->nsp);
      assert(surface->xsp == NULL);
      CTI_Mmap_rw(surface->xsp,surface->nsp_max); // _rw

      if (mpi_rank_shared == 0) {
        MPI_Allgather((double*)xsp_split,size_buf[3]*3,MPI_DOUBLE,
                      (double*)surface->xsp,size_buf[3]*3,MPI_DOUBLE,mpi_comm_internode);
        delete[] xsp_split;
      }

      if (mpi_rank == 0) cout << " > allgather on surface spost and znost..." << endl;

      surface->nst_max = size_buf[4]*mpi_size_internode; assert(surface->nst_max >= surface->nst);
      assert(surface->spost == NULL);
      CTI_Mmap(surface->spost,surface->nst_max);
      assert(surface->znost == NULL);
      CTI_Mmap(surface->znost,surface->nst_max);

      if (mpi_rank_shared == 0) {
        MPI_Allgather((int*)spost_split,size_buf[4]*3,MPI_INT,
                      (int*)surface->spost,size_buf[4]*3,MPI_INT,mpi_comm_internode);
        delete[] spost_split;
        MPI_Allgather(znost_split,size_buf[4],MPI_INT,
                      surface->znost,size_buf[4],MPI_INT,mpi_comm_internode);
        delete[] znost_split;
      }
      MPI_Barrier(mpi_comm_shared);

    }

    if (ff_surface) {

      // zone names get bcast to EVERYONE...
      if (mpi_rank != 0) {
        assert(ff_char_buf == NULL);
        ff_char_buf = new char[zone_name_len_max*ff_nzones];
      }
      MPI_Bcast(ff_char_buf,zone_name_len_max*ff_nzones,MPI_CHAR,0,mpi_comm);
      if (mpi_rank != 0) {
        ff_surface->zoneVec.resize(ff_nzones);
        for (int izone = 0; izone < ff_nzones; ++izone) {
          ff_surface->zoneVec[izone].setName(ff_char_buf+zone_name_len_max*izone);
        }
      }
      delete[] ff_char_buf;

      // allocate the shared memory for xsp AND spost AND znost

      if (mpi_rank == 0) cout << " > allgather on ff_surface xsp..." << endl;

      ff_surface->nsp_max = size_buf[8]*mpi_size_internode; assert(ff_surface->nsp_max >= ff_surface->nsp);
      assert(ff_surface->xsp == NULL);
      CTI_Mmap_rw(ff_surface->xsp,ff_surface->nsp_max);
      if (mpi_rank_shared == 0) {
        MPI_Allgather((double*)ff_xsp_split,size_buf[8]*3,MPI_DOUBLE,
                      (double*)ff_surface->xsp,size_buf[8]*3,MPI_DOUBLE,mpi_comm_internode);
        delete[] ff_xsp_split;
      }

      if (mpi_rank == 0) cout << " > allgather on ff_surface spost and znost..." << endl;

      ff_surface->nst_max = size_buf[9]*mpi_size_internode; assert(ff_surface->nst_max >= ff_surface->nst);
      assert(ff_surface->spost == NULL);
      CTI_Mmap(ff_surface->spost,ff_surface->nst_max);
      assert(ff_surface->znost == NULL);
      CTI_Mmap(ff_surface->znost,ff_surface->nst_max);
      assert(ff_surface_dxost == NULL);
      CTI_Mmap(ff_surface_dxost,ff_surface->nst_max);
      if (mpi_rank_shared == 0) {
        MPI_Allgather((int*)ff_spost_split,size_buf[9]*3,MPI_INT,
                      (int*)ff_surface->spost,size_buf[9]*3,MPI_INT,mpi_comm_internode);
        delete[] ff_spost_split;
        MPI_Allgather(ff_znost_split,size_buf[9],MPI_INT,
                      ff_surface->znost,size_buf[9],MPI_INT,mpi_comm_internode);
        delete[] ff_znost_split;
        MPI_Allgather(ff_dxost_split,size_buf[9],MPI_DOUBLE,
                      ff_surface_dxost,size_buf[9],MPI_DOUBLE,mpi_comm_internode);
        delete[] ff_dxost_split;
      }
      MPI_Barrier(mpi_comm_shared);

    }

    if (pts) {

      // the points are currently striped across rank0 of the shared mem nodes. Take
      // the data on each rank0 and scatter to all ranks on the node. pts do NOT
      // use shared memory...

      int np_scatter = my_np/mpi_size_shared;
      if (my_np%mpi_size_shared != 0) ++np_scatter;

      int ip0 = min(my_np,np_scatter*mpi_rank_shared);
      int ip1 = min(my_np,np_scatter*(mpi_rank_shared+1));
      pts->np = ip1-ip0; assert(pts->np <= np_scatter);

      // TODO delete this check
      {
        int np_global_check;
        MPI_Allreduce(&pts->np,&np_global_check,1,MPI_INT,MPI_SUM,mpi_comm);
        assert(np_global_check == pts->np_global);
      }

      pts->xp = new double[np_scatter][3];
      MPI_Scatter((double*)xp_split,np_scatter*3,MPI_DOUBLE, // note xp_split is padded to allow this
                  (double*)pts->xp,np_scatter*3,MPI_DOUBLE,
                  0,mpi_comm_shared);
      if (mpi_rank_shared == 0) delete[] xp_split;

      pts->delta = new double[np_scatter];
      MPI_Scatter((double*)deltap_split,np_scatter,MPI_DOUBLE, // note deltap_split is padded to allow this
                  (double*)pts->delta,np_scatter,MPI_DOUBLE,
                  0,mpi_comm_shared);
      if (mpi_rank_shared == 0) delete[] deltap_split;


      //pts->writeTecplot("junk.dat");


    }

    if (mpi_rank == 0) cout << " > done" << endl;

    // for now, hack these settings. Eventually they should be part of the part file format...
    // we limit these parts for now to fans and blade rows in disks where ALL points participate
    // in the calculation...

    // actually starting to change this...

    // TODO: include simple ff's in this file. These should run for the same
    // problem and not blank any points, and will also support the pingpong/bunny cases...

    assert(surface);
    assert(pts);
    assert(ff_surface == NULL);
    assert(ff_type == UNKNOWN_FF);
    ff_type = NO_FF;
    assert(solid_type == UNKNOWN_SOLID);
    //solid_type = NO_SOLID; // we used this for rotor cases
    solid_type = SURFACE_SOLID;
    return 0;

  }

  int Part::writeBinary(const string& filename) {

    if (mpi_rank == 0) cout << "Part::writeBinary: " << filename << endl;

    // the binary part file contains the following:
    // 1. surface
    // 2. ff_surface
    // 3. pts

    // create tmp file to write too...
    MiscUtils::mkdir_for_file_collective(filename,0);
    string tmp_filename = MiscUtils::makeTmpPrefix(filename);

    char dummy[128];
    sprintf(dummy,"%s",tmp_filename.c_str());
    MPI_File fh;
    MPI_Offset offset = 0;
    MPI_File_open(mpi_comm,dummy,MPI_MODE_WRONLY|MPI_MODE_CREATE,MPI_INFO_NULL,&fh);

    // ===============================================================
    // header
    // ===============================================================

    int itmp[5] = { PART_IO_MAGIC_NUMBER, PART_IO_VERSION, 0, 0, 0 };
    if (surface) itmp[2] = 1;
    if (ff_surface) itmp[3] = 1;
    if (pts) itmp[4] = 1;

    if (mpi_rank == 0) {

      cout << " > surface: " << itmp[2] << ", ff_surface: " << itmp[3] << ", pts: " << itmp[4] << endl;

      MPI_File_write(fh,itmp,5,MPI_INT,MPI_STATUS_IGNORE);
      offset += int_size*5;

      if (surface) {

        // ===============================================================
        // surface
        // ===============================================================

        int surface_itmp[3] = { int(surface->zoneVec.size()), surface->nsp, surface->nst };

        cout << " > surface nzones: " << surface_itmp[0] << ", nsp: " << surface_itmp[1] << ", nst: " << surface_itmp[2] << endl;

        MPI_File_write(fh,surface_itmp,3,MPI_INT,MPI_STATUS_IGNORE);
        offset += int_size*3;

        // surface zones...

        for (int izone = 0; izone < surface->zoneVec.size(); ++izone) {
          int length = surface->zoneVec[izone].getName().length();
          MPI_File_write(fh,&length,1,MPI_INT,MPI_STATUS_IGNORE);
          offset += int_size;
          MPI_File_write(fh,(char*)surface->zoneVec[izone].getName().c_str(),length,MPI_CHAR,MPI_STATUS_IGNORE);
          offset += length;
          cout << " > zone " << izone << " name \"" << surface->zoneVec[izone].getName() << "\"" << endl;
        }

        // surface xsp...

        MPI_File_write(fh,(double*)surface->xsp,surface->nsp*3,MPI_DOUBLE,MPI_STATUS_IGNORE);
        offset += double_size*surface->nsp*3;

        // surface spost...

        MPI_File_write(fh,(int*)surface->spost,surface->nst*3,MPI_INT,MPI_STATUS_IGNORE);
        offset += int_size*surface->nst*3;

        // surface znost...

        MPI_File_write(fh,(int*)surface->znost,surface->nst,MPI_INT,MPI_STATUS_IGNORE);
        offset += int_size*surface->nst;

      }

      if (ff_surface) {

        // ===============================================================
        // ff_surface
        // ===============================================================

        int ff_surface_itmp[3] = { int(ff_surface->zoneVec.size()), ff_surface->nsp, ff_surface->nst };

        cout << " > ff_surface nzones: " << ff_surface_itmp[0] << ", nsp: " << ff_surface_itmp[1] << ", nst: " << ff_surface_itmp[2] << endl;

        MPI_File_write(fh,ff_surface_itmp,3,MPI_INT,MPI_STATUS_IGNORE);
        offset += int_size*3;

        // ff_surface zones...

        for (int izone = 0; izone < ff_surface->zoneVec.size(); ++izone) {
          int length = ff_surface->zoneVec[izone].getName().length();
          MPI_File_write(fh,&length,1,MPI_INT,MPI_STATUS_IGNORE);
          offset += int_size;
          MPI_File_write(fh,(char*)ff_surface->zoneVec[izone].getName().c_str(),length,MPI_CHAR,MPI_STATUS_IGNORE);
          offset += length;
          cout << " > zone " << izone << " name \"" << ff_surface->zoneVec[izone].getName() << "\"" << endl;
        }

        // ff_surface xsp...

        MPI_File_write(fh,(double*)ff_surface->xsp,ff_surface->nsp*3,MPI_DOUBLE,MPI_STATUS_IGNORE);
        offset += double_size*ff_surface->nsp*3;

        // ff_surface spost...

        MPI_File_write(fh,(int*)ff_surface->spost,ff_surface->nst*3,MPI_INT,MPI_STATUS_IGNORE);
        offset += int_size*ff_surface->nst*3;

        // ff_surface znost...

        MPI_File_write(fh,(int*)ff_surface->znost,ff_surface->nst,MPI_INT,MPI_STATUS_IGNORE);
        offset += int_size*ff_surface->nst;

        // ff_surface_dxost...

        assert(ff_surface_dxost);
        MPI_File_write(fh,ff_surface_dxost,ff_surface->nst,MPI_DOUBLE,MPI_STATUS_IGNORE);
        offset += double_size*ff_surface->nst;

      }

    }

    MPI_Bcast(&offset,1,MPI_INT8,0,mpi_comm);

    // ===============================
    // finally points...
    // ===============================

    if (pts) {

      // now all ranks count (and then write) all their local points...

      int my_disp;
      MPI_Scan(&pts->np,&my_disp,1,MPI_INT,MPI_SUM,mpi_comm);

      int np_global = my_disp;
      MPI_Bcast(&np_global,1,MPI_INT,mpi_size-1,mpi_comm);
      assert(np_global > 0);

      if (mpi_rank == 0) {
        cout << " > total point count: " << np_global << endl;
        MPI_File_write(fh,&np_global,1,MPI_INT,MPI_STATUS_IGNORE);
      }
      offset += int_size;

      // now xp...

      writeChunkedData<double>(fh,offset+double_size*(my_disp-pts->np)*3,(double*)pts->xp,pts->np*3,mpi_comm);
      offset += double_size*np_global*3;

      // delta...

      writeChunkedData<double>(fh,offset+double_size*(my_disp-pts->np),pts->delta,pts->np,mpi_comm);
      offset += double_size*np_global;

    }

    MPI_File_set_size(fh,offset);
    MPI_File_close(&fh);

    // if we successfully wrote the tmp file, rename/replace the file with it...
    
    if (mpi_rank == 0) {
      remove(filename.c_str());
      rename(tmp_filename.c_str(),filename.c_str());
    }

    return 0;
  }

  void Part::prepareIsInsideFF(const int xy_min[2],const int xy_max[2]) {

    assert(!hasSimpleFF());

    /*
      if (surface) {
      assert(surface_zn_flag); // why?
      assert(ff_surface == NULL);
      }
      else {
      assert(surface_zn_flag == NULL);
      assert(ff_surface);
      }
    */

    {
      if (adt2d_ff != NULL) {
	if ((xy_min[0] >= bbmin_adt2d_ff[0])&&(xy_min[1] >= bbmin_adt2d_ff[1])&&
	    (xy_max[0] <= bbmax_adt2d_ff[0])&&(xy_max[1] <= bbmax_adt2d_ff[1]))
	  return;
	if (mpi_rank == 0) cout << "rebuilding adt2d_ff XXXXXXXXXXXXXXXXXX" << endl;
	delete adt2d_ff;
	adt2d_ff = NULL;
	stobb_ff.clear();
      }

      FOR_I2 bbmin_adt2d_ff[i] = xy_min[i];
      FOR_I2 bbmax_adt2d_ff[i] = xy_max[i];

      // TODO: this should be a part routine, but figure out how to access namspace
      // stuff from the object...

      // now put any FF surface tris that possibly intersect this bbox into stobb_ff...
      assert(stobb_ff.empty());
      vector<int> bbVec;
      int xy0[2],xy1[2],xy2[2];

      if (ff_surface) {
	for (int ist = 0; ist < ff_surface->nst; ++ist) {
	  const double * const xsp0 = ff_surface->xsp[ff_surface->spost[ist][0]];
	  xy0[0] = getOddInt((xsp0[0]-dtoi_x0[0])*e1_ff[0] +
			     (xsp0[1]-dtoi_x0[1])*e1_ff[1] +
			     (xsp0[2]-dtoi_x0[2])*e1_ff[2],dtoi_delta);
	  const double * const xsp1 = ff_surface->xsp[ff_surface->spost[ist][1]];
	  xy1[0] = getOddInt((xsp1[0]-dtoi_x0[0])*e1_ff[0] +
			     (xsp1[1]-dtoi_x0[1])*e1_ff[1] +
			     (xsp1[2]-dtoi_x0[2])*e1_ff[2],dtoi_delta);
	  const double * const xsp2 = ff_surface->xsp[ff_surface->spost[ist][2]];
	  xy2[0] = getOddInt((xsp2[0]-dtoi_x0[0])*e1_ff[0] +
			     (xsp2[1]-dtoi_x0[1])*e1_ff[1] +
			     (xsp2[2]-dtoi_x0[2])*e1_ff[2],dtoi_delta);
	  const int bbmax0 = max(xy0[0],max(xy1[0],xy2[0]));
	  if (bbmax0 < bbmin_adt2d_ff[0])
	    continue;
	  const int bbmin0 = min(xy0[0],min(xy1[0],xy2[0]));
	  if (bbmin0 > bbmax_adt2d_ff[0])
	    continue;
	  xy0[1] = getOddInt((xsp0[0]-dtoi_x0[0])*e2_ff[0] +
			     (xsp0[1]-dtoi_x0[1])*e2_ff[1] +
			     (xsp0[2]-dtoi_x0[2])*e2_ff[2],dtoi_delta);
	  xy1[1] = getOddInt((xsp1[0]-dtoi_x0[0])*e2_ff[0] +
			     (xsp1[1]-dtoi_x0[1])*e2_ff[1] +
			     (xsp1[2]-dtoi_x0[2])*e2_ff[2],dtoi_delta);
	  xy2[1] = getOddInt((xsp2[0]-dtoi_x0[0])*e2_ff[0] +
			     (xsp2[1]-dtoi_x0[1])*e2_ff[1] +
			     (xsp2[2]-dtoi_x0[2])*e2_ff[2],dtoi_delta);
	  const int bbmax1 = max(xy0[1],max(xy1[1],xy2[1]));
	  if (bbmax1 < bbmin_adt2d_ff[1])
	    continue;
	  const int bbmin1 = min(xy0[1],min(xy1[1],xy2[1]));
	  if (bbmin1 > bbmax_adt2d_ff[1])
	    continue;
	  // if we made it here, we need to include this ist...
	  stobb_ff.push_back(ist);
	  bbVec.push_back(bbmin0);
	  bbVec.push_back(bbmin1);
	  bbVec.push_back(bbmax0);
	  bbVec.push_back(bbmax1);
	}
      }
      else if (surface) {
	assert(surface_zn_flag);
	for (int ist = 0; ist < surface->nst; ++ist) {
	  if (surface_zn_flag[surface->znost[ist]] == -1) {
	    const double * const xsp0 = surface->xsp[surface->spost[ist][0]];
	    xy0[0] = getOddInt((xsp0[0]-dtoi_x0[0])*e1_ff[0] +
			       (xsp0[1]-dtoi_x0[1])*e1_ff[1] +
			       (xsp0[2]-dtoi_x0[2])*e1_ff[2],dtoi_delta);
	    const double * const xsp1 = surface->xsp[surface->spost[ist][1]];
	    xy1[0] = getOddInt((xsp1[0]-dtoi_x0[0])*e1_ff[0] +
			       (xsp1[1]-dtoi_x0[1])*e1_ff[1] +
			       (xsp1[2]-dtoi_x0[2])*e1_ff[2],dtoi_delta);
	    const double * const xsp2 = surface->xsp[surface->spost[ist][2]];
	    xy2[0] = getOddInt((xsp2[0]-dtoi_x0[0])*e1_ff[0] +
			       (xsp2[1]-dtoi_x0[1])*e1_ff[1] +
			       (xsp2[2]-dtoi_x0[2])*e1_ff[2],dtoi_delta);
	    const int bbmax0 = max(xy0[0],max(xy1[0],xy2[0]));
	    if (bbmax0 < bbmin_adt2d_ff[0])
	      continue;
	    const int bbmin0 = min(xy0[0],min(xy1[0],xy2[0]));
	    if (bbmin0 > bbmax_adt2d_ff[0])
	      continue;
	    xy0[1] = getOddInt((xsp0[0]-dtoi_x0[0])*e2_ff[0] +
			       (xsp0[1]-dtoi_x0[1])*e2_ff[1] +
			       (xsp0[2]-dtoi_x0[2])*e2_ff[2],dtoi_delta);
	    xy1[1] = getOddInt((xsp1[0]-dtoi_x0[0])*e2_ff[0] +
			       (xsp1[1]-dtoi_x0[1])*e2_ff[1] +
			       (xsp1[2]-dtoi_x0[2])*e2_ff[2],dtoi_delta);
	    xy2[1] = getOddInt((xsp2[0]-dtoi_x0[0])*e2_ff[0] +
			       (xsp2[1]-dtoi_x0[1])*e2_ff[1] +
			       (xsp2[2]-dtoi_x0[2])*e2_ff[2],dtoi_delta);
	    const int bbmax1 = max(xy0[1],max(xy1[1],xy2[1]));
	    if (bbmax1 < bbmin_adt2d_ff[1])
	      continue;
	    const int bbmin1 = min(xy0[1],min(xy1[1],xy2[1]));
	    if (bbmin1 > bbmax_adt2d_ff[1])
	      continue;
	    // if we made it here, we need to include this ist...
	    stobb_ff.push_back(ist);
	    bbVec.push_back(bbmin0);
	    bbVec.push_back(bbmin1);
	    bbVec.push_back(bbmax0);
	    bbVec.push_back(bbmax1);
	  }
	}
      }
      else {
	assert(0);
      }
      const int nstobb_ff = stobb_ff.size();
      int (*bbmin)[2] = new int[nstobb_ff][2];
      int (*bbmax)[2] = new int[nstobb_ff][2];
      for (int ibb = 0; ibb < nstobb_ff; ++ibb) {
	FOR_I2 bbmin[ibb][i] = bbVec[ibb*4+i];
	FOR_I2 bbmax[ibb][i] = bbVec[ibb*4+2+i];
      }
      adt2d_ff = new Adt2d<int>(nstobb_ff,bbmin,bbmax);
      delete[] bbmin;
      delete[] bbmax;
    }

    if (b_link_surface) {
      if (adt2d_ff_surface != NULL) {
	if ((xy_min[0] >= bbmin_adt2d_ff_surface[0])&&(xy_min[1] >= bbmin_adt2d_ff_surface[1])&&
	    (xy_max[0] <= bbmax_adt2d_ff_surface[0])&&(xy_max[1] <= bbmax_adt2d_ff_surface[1]))
	  return;
	if (mpi_rank == 0) cout << "rebuilding adt2d_ff_surface XXXXXXXXXXXXXXXXXX" << endl;
	delete adt2d_ff_surface;
	adt2d_ff_surface = NULL;
	stobb_ff_surface.clear();
      }

      FOR_I2 bbmin_adt2d_ff_surface[i] = xy_min[i];
      FOR_I2 bbmax_adt2d_ff_surface[i] = xy_max[i];

      // now put any surface tris that possibly intersect this bbox into stobb_ff_surface...
      assert(stobb_ff_surface.empty());
      vector<int> bbVec;
      int xy0[2],xy1[2],xy2[2];
      assert(ipart_link.size()>0);
      for (int ii = 0; ii < ipart_link.size(); ++ii) {
	const int ipart_idx = ipart_link[ii];
	const int ist_offset = ist_offset_link[ii];
	SurfaceShm * surface = partVec[ipart_idx]->surface;
	for (int ist = 0; ist < surface->nst; ++ist) {
	  int ipart;
	  partVec[ipart_idx]->getPartForSt(ipart,ist);
	  if (ipart_of_part == ipart) {
	    const double * const xsp0 = surface->xsp[surface->spost[ist][0]];
	    xy0[0] = getOddInt((xsp0[0]-dtoi_x0[0])*e1_ff[0] +
			       (xsp0[1]-dtoi_x0[1])*e1_ff[1] +
			       (xsp0[2]-dtoi_x0[2])*e1_ff[2],dtoi_delta);
	    const double * const xsp1 = surface->xsp[surface->spost[ist][1]];
	    xy1[0] = getOddInt((xsp1[0]-dtoi_x0[0])*e1_ff[0] +
			       (xsp1[1]-dtoi_x0[1])*e1_ff[1] +
			       (xsp1[2]-dtoi_x0[2])*e1_ff[2],dtoi_delta);
	    const double * const xsp2 = surface->xsp[surface->spost[ist][2]];
	    xy2[0] = getOddInt((xsp2[0]-dtoi_x0[0])*e1_ff[0] +
			       (xsp2[1]-dtoi_x0[1])*e1_ff[1] +
			       (xsp2[2]-dtoi_x0[2])*e1_ff[2],dtoi_delta);
	    const int bbmax0 = max(xy0[0],max(xy1[0],xy2[0]));
	    if (bbmax0 < bbmin_adt2d_ff_surface[0])
	      continue;
	    const int bbmin0 = min(xy0[0],min(xy1[0],xy2[0]));
	    if (bbmin0 > bbmax_adt2d_ff_surface[0])
	      continue;
	    xy0[1] = getOddInt((xsp0[0]-dtoi_x0[0])*e2_ff[0] +
			       (xsp0[1]-dtoi_x0[1])*e2_ff[1] +
			       (xsp0[2]-dtoi_x0[2])*e2_ff[2],dtoi_delta);
	    xy1[1] = getOddInt((xsp1[0]-dtoi_x0[0])*e2_ff[0] +
			       (xsp1[1]-dtoi_x0[1])*e2_ff[1] +
			       (xsp1[2]-dtoi_x0[2])*e2_ff[2],dtoi_delta);
	    xy2[1] = getOddInt((xsp2[0]-dtoi_x0[0])*e2_ff[0] +
			       (xsp2[1]-dtoi_x0[1])*e2_ff[1] +
			       (xsp2[2]-dtoi_x0[2])*e2_ff[2],dtoi_delta);
	    const int bbmax1 = max(xy0[1],max(xy1[1],xy2[1]));
	    if (bbmax1 < bbmin_adt2d_ff_surface[1])
	      continue;
	    const int bbmin1 = min(xy0[1],min(xy1[1],xy2[1]));
	    if (bbmin1 > bbmax_adt2d_ff_surface[1])
	      continue;
	    // if we made it here, we need to include this ist...
	    stobb_ff_surface.push_back(ist+ist_offset);
	    bbVec.push_back(bbmin0);
	    bbVec.push_back(bbmin1);
	    bbVec.push_back(bbmax0);
	    bbVec.push_back(bbmax1);
	  }
	}
      }
      const int nstobb_ff_surface = stobb_ff_surface.size();
      int (*bbmin)[2] = new int[nstobb_ff_surface][2];
      int (*bbmax)[2] = new int[nstobb_ff_surface][2];
      for (int ibb = 0; ibb < nstobb_ff_surface; ++ibb) {
	FOR_I2 bbmin[ibb][i] = bbVec[ibb*4+i];
	FOR_I2 bbmax[ibb][i] = bbVec[ibb*4+2+i];
      }
      adt2d_ff_surface = new Adt2d<int>(nstobb_ff_surface,bbmin,bbmax);
      delete[] bbmin;
      delete[] bbmax;
    }

  }

  bool Part::isInsideFF(const double xp[3]) {

    switch(ff_type) {

    case ANNULAR_CYLINDER_FF:
      {
        // =======================================================
        // the annular cylinder has the following in data_ff:
        // x0[3] : ddata_ff[0,1,2]
        // x1[3] : ddata_ff[3,4,5]
        // r0    : ddata_ff[6]
        // r1    : ddata_ff[7];
        // =======================================================
        const double * const x0 = ddata_ff_+0;
        const double * const x1 = ddata_ff_+3;
        const double dx[3] = DIFF(x1,x0);
        const double dxp[3] = DIFF(xp,x0);
        const double dp = DOT_PRODUCT(dxp,dx);
        // bias towards true, i.e. isInsideFF...
        if (dp < 0.0) {
          return false;
        }
        else {
          const double dx_mag2 = DOT_PRODUCT(dx,dx);
          if (dp > dx_mag2) {
            return false;
          }
          else {
            const double r0 = ddata_ff_[6];
            const double r1 = ddata_ff_[7];
            const double r2 = DOT_PRODUCT(dxp,dxp) - dp*dp/dx_mag2;
            if ((r2 < r0*r0)||(r2 > r1*r1))
              return false;
          }
        }
        return true;
      }
      break;

    case CYLINDER_FF:
      {
        // =======================================================
        // the cylinder has the following in data_ff:
        // x0[3] : ddata_ff[0,1,2]
        // x1[3] : ddata_ff[3,4,5]
        // r    : ddata_ff[6]
        // =======================================================
        const double * const x0 = ddata_ff_+0;
        const double * const x1 = ddata_ff_+3;
        const double dx[3] = DIFF(x1,x0);
        const double dxp[3] = DIFF(xp,x0);
        const double dp = DOT_PRODUCT(dxp,dx);
        // bias towards true, i.e. isInsideFF...
        if (dp < 0.0) {
          return false;
        }
        else {
          const double dx_mag2 = DOT_PRODUCT(dx,dx);
          if (dp > dx_mag2) {
            return false;
          }
          else {
            const double r = ddata_ff_[6];
            const double r2 = DOT_PRODUCT(dxp,dxp) - dp*dp/dx_mag2;
            if (r2 > r*r)return false;
          }
        }
        return true;
      }
      break;
    
    case BOX_FF:
      {
        // =======================================================
        // the box has the following in data_ff:
        // x0[3] : ddata_ff[0,1,2]
        // x1[3] : ddata_ff[3,4,5]
        // =======================================================
        FOR_I3 {
          if ((xp[i] < ddata_ff_[i]) || (xp[i] > ddata_ff_[3+i])) return false;
        }
        return true;
      }
      break;
    
    case SPHERE_FF:
      {
        // =======================================================
        // the sphere has the following in data_ff:
        // x[3] : ddata_ff[0,1,2]
        // r    : ddata_ff[3];
        // =======================================================
        const double * const x = ddata_ff_;
        const double r = ddata_ff_[3];
        if (DIST2(xp,x) > r*r)
          return false;
      }
      return true;

    case SURFACE_FAZONE_FF:
      {
        // get our point as an even int...
        assert(b_dtoi);

        int xy[2];
        xy[0] = getEvenInt((xp[0]-dtoi_x0[0])*e1_ff[0] +
                           (xp[1]-dtoi_x0[1])*e1_ff[1] +
                           (xp[2]-dtoi_x0[2])*e1_ff[2],dtoi_delta);
        xy[1] = getEvenInt((xp[0]-dtoi_x0[0])*e2_ff[0] +
                           (xp[1]-dtoi_x0[1])*e2_ff[1] +
                           (xp[2]-dtoi_x0[2])*e2_ff[2],dtoi_delta);
        // any point we check must be inside the current bbmin/max_adt2d_ff:
        // if not, then the prepare was not called or did not do its job...
        assert((xy[0] >= bbmin_adt2d_ff[0])&&(xy[0] <= bbmax_adt2d_ff[0]));
        assert((xy[1] >= bbmin_adt2d_ff[1])&&(xy[1] <= bbmax_adt2d_ff[1]));
        vector<int> bboxVec;
        adt2d_ff->buildListForPoint(bboxVec,xy);
        int closest = 0;
        double dist_closest;
        for (int ibb = 0; ibb < bboxVec.size(); ++ibb) {
          const int ist = stobb_ff[bboxVec[ibb]];
          // and get the tri corners as odd int's...
          const double * const xsp0 = surface->xsp[surface->spost[ist][0]];
          int xy0[2];
          xy0[0] = getOddInt((xsp0[0]-dtoi_x0[0])*e1_ff[0] + (xsp0[1]-dtoi_x0[1])*e1_ff[1] + (xsp0[2]-dtoi_x0[2])*e1_ff[2],dtoi_delta);
          xy0[1] = getOddInt((xsp0[0]-dtoi_x0[0])*e2_ff[0] + (xsp0[1]-dtoi_x0[1])*e2_ff[1] + (xsp0[2]-dtoi_x0[2])*e2_ff[2],dtoi_delta);
          const double * const xsp1 = surface->xsp[surface->spost[ist][1]];
          int xy1[2];
          xy1[0] = getOddInt((xsp1[0]-dtoi_x0[0])*e1_ff[0] + (xsp1[1]-dtoi_x0[1])*e1_ff[1] + (xsp1[2]-dtoi_x0[2])*e1_ff[2],dtoi_delta);
          xy1[1] = getOddInt((xsp1[0]-dtoi_x0[0])*e2_ff[0] + (xsp1[1]-dtoi_x0[1])*e2_ff[1] + (xsp1[2]-dtoi_x0[2])*e2_ff[2],dtoi_delta);
          const double * const xsp2 = surface->xsp[surface->spost[ist][2]];
          int xy2[2];
          xy2[0] = getOddInt((xsp2[0]-dtoi_x0[0])*e1_ff[0] + (xsp2[1]-dtoi_x0[1])*e1_ff[1] + (xsp2[2]-dtoi_x0[2])*e1_ff[2],dtoi_delta);
          xy2[1] = getOddInt((xsp2[0]-dtoi_x0[0])*e2_ff[0] + (xsp2[1]-dtoi_x0[1])*e2_ff[1] + (xsp2[2]-dtoi_x0[2])*e2_ff[2],dtoi_delta);
          const int8 dxy0[2] = { xy0[0]-xy[0], xy0[1]-xy[1] };
          const int8 dxy1[2] = { xy1[0]-xy[0], xy1[1]-xy[1] };
          const int8 dxy2[2] = { xy2[0]-xy[0], xy2[1]-xy[1] };
          const int8 A0 = dxy1[0]*dxy2[1] - dxy1[1]*dxy2[0];
          const int8 A1 = dxy2[0]*dxy0[1] - dxy2[1]*dxy0[0];
          const int8 A2 = dxy0[0]*dxy1[1] - dxy0[1]*dxy1[0];
          // skip slits...
          if (A0+A1+A2 == 0)
            continue;
          if ((A0 >= 0)&&(A1 >= 0)&&(A2 >= 0)) {
            //assert(A0+A1+A2 > 0); // check for slit
            double dist = 0.0;
            FOR_I3 {
              const double dxi = (double(A0)*(xsp0[i]-xp[i]) + double(A1)*(xsp1[i]-xp[i]) + double(A2)*(xsp2[i]-xp[i]))/double(A0+A1+A2);
              dist += dxi*e0_ff[i];
            }
            if ((closest == 0)||(fabs(dist) < fabs(dist_closest))) {
              closest = 1; // i.e. Areas positive
              dist_closest = dist;
            }
          }
          else if ((A0 <= 0)&&(A1 <= 0)&&(A2 <= 0)) {
            //assert(A0+A1+A2 < 0);
            double dist = 0.0;
            FOR_I3 {
              const double dxi = (double(A0)*(xsp0[i]-xp[i]) + double(A1)*(xsp1[i]-xp[i]) + double(A2)*(xsp2[i]-xp[i]))/double(A0+A1+A2);
              dist += dxi*e0_ff[i];
            }
            if ((closest == 0)||(fabs(dist) < fabs(dist_closest))) {
              closest = -1; // i.e. Areas positive
              dist_closest = dist;
            }
          }
        }
        if ( (closest == 0) || ((closest == 1)&&(dist_closest <= 0.0)) || ((closest == -1)&&(dist_closest >= 0.0)) ) return false;
        assert(dist_closest != 0.0);
        return true;
      }

    case FF_SURFACE_FF:
    case FF_SURFACE_FAZONE_FF: // solid check will be handled elsewhere
      {
        // get our point as an even int...
        assert(b_dtoi);

        int xy[2];
        xy[0] = getEvenInt((xp[0]-dtoi_x0[0])*e1_ff[0] +
                           (xp[1]-dtoi_x0[1])*e1_ff[1] +
                           (xp[2]-dtoi_x0[2])*e1_ff[2],dtoi_delta);
        xy[1] = getEvenInt((xp[0]-dtoi_x0[0])*e2_ff[0] +
                           (xp[1]-dtoi_x0[1])*e2_ff[1] +
                           (xp[2]-dtoi_x0[2])*e2_ff[2],dtoi_delta);
        // any point we check must be inside the current bbmin/max_adt2d_ff:
        // if not, then the prepare was not called or did not do its job...
        assert((xy[0] >= bbmin_adt2d_ff[0])&&(xy[0] <= bbmax_adt2d_ff[0]));
        assert((xy[1] >= bbmin_adt2d_ff[1])&&(xy[1] <= bbmax_adt2d_ff[1]));
        vector<int> bboxVec;
        adt2d_ff->buildListForPoint(bboxVec,xy);
	
	// when there is a linked surface ...
	vector<int> bboxVec_surface;
	if (b_link_surface) {
	  assert((xy[0] >= bbmin_adt2d_ff_surface[0])&&(xy[0] <= bbmax_adt2d_ff_surface[0]));
	  assert((xy[1] >= bbmin_adt2d_ff_surface[1])&&(xy[1] <= bbmax_adt2d_ff_surface[1]));	  
	  adt2d_ff_surface->buildListForPoint(bboxVec_surface,xy);
	}

        int closest = 0;
        double dist_closest;
        for (int ibb = 0; ibb < bboxVec.size(); ++ibb) {
          const int ist = stobb_ff[bboxVec[ibb]];
          // and get the tri corners as odd int's...
          const double * const xsp0 = ff_surface->xsp[ff_surface->spost[ist][0]];
          int xy0[2];
          xy0[0] = getOddInt((xsp0[0]-dtoi_x0[0])*e1_ff[0] + (xsp0[1]-dtoi_x0[1])*e1_ff[1] + (xsp0[2]-dtoi_x0[2])*e1_ff[2],dtoi_delta);
          xy0[1] = getOddInt((xsp0[0]-dtoi_x0[0])*e2_ff[0] + (xsp0[1]-dtoi_x0[1])*e2_ff[1] + (xsp0[2]-dtoi_x0[2])*e2_ff[2],dtoi_delta);
          const double * const xsp1 = ff_surface->xsp[ff_surface->spost[ist][1]];
          int xy1[2];
          xy1[0] = getOddInt((xsp1[0]-dtoi_x0[0])*e1_ff[0] + (xsp1[1]-dtoi_x0[1])*e1_ff[1] + (xsp1[2]-dtoi_x0[2])*e1_ff[2],dtoi_delta);
          xy1[1] = getOddInt((xsp1[0]-dtoi_x0[0])*e2_ff[0] + (xsp1[1]-dtoi_x0[1])*e2_ff[1] + (xsp1[2]-dtoi_x0[2])*e2_ff[2],dtoi_delta);
          const double * const xsp2 = ff_surface->xsp[ff_surface->spost[ist][2]];
          int xy2[2];
          xy2[0] = getOddInt((xsp2[0]-dtoi_x0[0])*e1_ff[0] + (xsp2[1]-dtoi_x0[1])*e1_ff[1] + (xsp2[2]-dtoi_x0[2])*e1_ff[2],dtoi_delta);
          xy2[1] = getOddInt((xsp2[0]-dtoi_x0[0])*e2_ff[0] + (xsp2[1]-dtoi_x0[1])*e2_ff[1] + (xsp2[2]-dtoi_x0[2])*e2_ff[2],dtoi_delta);
          const int8 dxy0[2] = { xy0[0]-xy[0], xy0[1]-xy[1] };
          const int8 dxy1[2] = { xy1[0]-xy[0], xy1[1]-xy[1] };
          const int8 dxy2[2] = { xy2[0]-xy[0], xy2[1]-xy[1] };
          const int8 A0 = dxy1[0]*dxy2[1] - dxy1[1]*dxy2[0];
          const int8 A1 = dxy2[0]*dxy0[1] - dxy2[1]*dxy0[0];
          const int8 A2 = dxy0[0]*dxy1[1] - dxy0[1]*dxy1[0];
          // skip slits...
          if (A0+A1+A2 == 0)
            continue;
          if ((A0 >= 0)&&(A1 >= 0)&&(A2 >= 0)) {
            //assert(A0+A1+A2 > 0); // check for slit
            double dist = 0.0;
            FOR_I3 {
              const double dxi = (double(A0)*(xsp0[i]-xp[i]) + double(A1)*(xsp1[i]-xp[i]) + double(A2)*(xsp2[i]-xp[i]))/double(A0+A1+A2);
              dist += dxi*e0_ff[i];
            }
            if ((closest == 0)||(fabs(dist) < fabs(dist_closest))) {
              closest = 1; // i.e. Areas positive
              dist_closest = dist;
            }
          }
          else if ((A0 <= 0)&&(A1 <= 0)&&(A2 <= 0)) {
            //assert(A0+A1+A2 < 0);
            double dist = 0.0;
            FOR_I3 {
              const double dxi = (double(A0)*(xsp0[i]-xp[i]) + double(A1)*(xsp1[i]-xp[i]) + double(A2)*(xsp2[i]-xp[i]))/double(A0+A1+A2);
              dist += dxi*e0_ff[i];
            }
            if ((closest == 0)||(fabs(dist) < fabs(dist_closest))) {
              closest = -1; // i.e. Areas positive
              dist_closest = dist;
            }
          }
        }
	if (b_link_surface) {
	  for (int ii = 0; ii < ipart_link.size(); ++ii) {
	    const int ipart_idx = ipart_link[ii];
	    const int ist_offset = ist_offset_link[ii];
	    SurfaceShm * surface = partVec[ipart_idx]->surface;
	    for (int ibb = 0; ibb < bboxVec_surface.size(); ++ibb) {
	      const int ist = stobb_ff_surface[bboxVec_surface[ibb]] - ist_offset;
	      if (ist >= 0 && ist < surface->nst) {
		// and get the tri corners as odd int's...
		const double * const xsp0 = surface->xsp[surface->spost[ist][0]];
		int xy0[2];
		xy0[0] = getOddInt((xsp0[0]-dtoi_x0[0])*e1_ff[0] + (xsp0[1]-dtoi_x0[1])*e1_ff[1] + (xsp0[2]-dtoi_x0[2])*e1_ff[2],dtoi_delta);
		xy0[1] = getOddInt((xsp0[0]-dtoi_x0[0])*e2_ff[0] + (xsp0[1]-dtoi_x0[1])*e2_ff[1] + (xsp0[2]-dtoi_x0[2])*e2_ff[2],dtoi_delta);
		const double * const xsp1 = surface->xsp[surface->spost[ist][1]];
		int xy1[2];
		xy1[0] = getOddInt((xsp1[0]-dtoi_x0[0])*e1_ff[0] + (xsp1[1]-dtoi_x0[1])*e1_ff[1] + (xsp1[2]-dtoi_x0[2])*e1_ff[2],dtoi_delta);
		xy1[1] = getOddInt((xsp1[0]-dtoi_x0[0])*e2_ff[0] + (xsp1[1]-dtoi_x0[1])*e2_ff[1] + (xsp1[2]-dtoi_x0[2])*e2_ff[2],dtoi_delta);
		const double * const xsp2 = surface->xsp[surface->spost[ist][2]];
		int xy2[2];
		xy2[0] = getOddInt((xsp2[0]-dtoi_x0[0])*e1_ff[0] + (xsp2[1]-dtoi_x0[1])*e1_ff[1] + (xsp2[2]-dtoi_x0[2])*e1_ff[2],dtoi_delta);
		xy2[1] = getOddInt((xsp2[0]-dtoi_x0[0])*e2_ff[0] + (xsp2[1]-dtoi_x0[1])*e2_ff[1] + (xsp2[2]-dtoi_x0[2])*e2_ff[2],dtoi_delta);
		const int8 dxy0[2] = { xy0[0]-xy[0], xy0[1]-xy[1] };
		const int8 dxy1[2] = { xy1[0]-xy[0], xy1[1]-xy[1] };
		const int8 dxy2[2] = { xy2[0]-xy[0], xy2[1]-xy[1] };
		const int8 A0 = dxy1[0]*dxy2[1] - dxy1[1]*dxy2[0];
		const int8 A1 = dxy2[0]*dxy0[1] - dxy2[1]*dxy0[0];
		const int8 A2 = dxy0[0]*dxy1[1] - dxy0[1]*dxy1[0];
                // skip slits...
                if (A0+A1+A2 == 0)
                  continue;
		if ((A0 >= 0)&&(A1 >= 0)&&(A2 >= 0)) {
		  //assert(A0+A1+A2 > 0); // check for slit
		  double dist = 0.0;
		  FOR_I3 {
		    const double dxi = (double(A0)*(xsp0[i]-xp[i]) + double(A1)*(xsp1[i]-xp[i]) + double(A2)*(xsp2[i]-xp[i]))/double(A0+A1+A2);
		    dist += dxi*e0_ff[i];
		  }
		  if ((closest == 0)||(fabs(dist) < fabs(dist_closest))) {
		    closest = 1; // i.e. Areas positive
		    dist_closest = dist;
		  }
		}
		else if ((A0 <= 0)&&(A1 <= 0)&&(A2 <= 0)) {
		  //assert(A0+A1+A2 < 0);
		  double dist = 0.0;
		  FOR_I3 {
		    const double dxi = (double(A0)*(xsp0[i]-xp[i]) + double(A1)*(xsp1[i]-xp[i]) + double(A2)*(xsp2[i]-xp[i]))/double(A0+A1+A2);
		    dist += dxi*e0_ff[i];
		  }
		  if ((closest == 0)||(fabs(dist) < fabs(dist_closest))) {
		    closest = -1; // i.e. Areas positive
		    dist_closest = dist;
		  }
		}
	      }
	    }
	  }
	}
        if ( (closest == 0) || ((closest == 1)&&(dist_closest < 0.0)) || ((closest == -1)&&(dist_closest > 0.0)) ) return false;
        assert(dist_closest != 0.0);
        return true;
      }

    case NO_FF:

      // should never be called...
      assert(0);

    default:

      if (mpi_rank == 0) cout << "Warning: Part::isInsideFF for name \"" << name << "\": ff_type not handled: " << ff_type << endl;
      assert(0);

    }

    return false;

  }
  void Part::prepareIsInsideSolid(const int xy_min[2],const int xy_max[2]) {

    assert(!hasSimpleSolid());
    assert(surface);

    if (adt2d_solid != NULL) {
      if ((xy_min[0] >= bbmin_adt2d_solid[0])&&(xy_min[1] >= bbmin_adt2d_solid[1])&&
          (xy_max[0] <= bbmax_adt2d_solid[0])&&(xy_max[1] <= bbmax_adt2d_solid[1]))
        return;
      if (mpi_rank == 0) cout << "rebuilding adt2d_solid XXXXXXXXXXXXXXXXXX" << endl;
      delete adt2d_solid;
      adt2d_solid = NULL;
      stobb_solid.clear();
    }

    FOR_I2 bbmin_adt2d_solid[i] = xy_min[i];
    FOR_I2 bbmax_adt2d_solid[i] = xy_max[i];


    // now put any FF surface tris that possibly intersect this bbox into stobb_solid...
    assert(stobb_solid.empty());
    vector<int> bbVec;
    int xy0[2],xy1[2],xy2[2];

    for (int ist = 0; ist < surface->nst; ++ist) {
      if ((surface_zn_flag == NULL)||(surface_zn_flag[surface->znost[ist]] != -1)) {
        const double * const xsp0 = surface->xsp[surface->spost[ist][0]];
        xy0[0] = getOddInt((xsp0[0]-dtoi_x0[0])*e1_solid[0] +
                           (xsp0[1]-dtoi_x0[1])*e1_solid[1] +
                           (xsp0[2]-dtoi_x0[2])*e1_solid[2],dtoi_delta);
        const double * const xsp1 = surface->xsp[surface->spost[ist][1]];
        xy1[0] = getOddInt((xsp1[0]-dtoi_x0[0])*e1_solid[0] +
                           (xsp1[1]-dtoi_x0[1])*e1_solid[1] +
                           (xsp1[2]-dtoi_x0[2])*e1_solid[2],dtoi_delta);
        const double * const xsp2 = surface->xsp[surface->spost[ist][2]];
        xy2[0] = getOddInt((xsp2[0]-dtoi_x0[0])*e1_solid[0] +
                           (xsp2[1]-dtoi_x0[1])*e1_solid[1] +
                           (xsp2[2]-dtoi_x0[2])*e1_solid[2],dtoi_delta);
        const int bbmax0 = max(xy0[0],max(xy1[0],xy2[0]));
        if (bbmax0 < bbmin_adt2d_solid[0])
          continue;
        const int bbmin0 = min(xy0[0],min(xy1[0],xy2[0]));
        if (bbmin0 > bbmax_adt2d_solid[0])
          continue;
        xy0[1] = getOddInt((xsp0[0]-dtoi_x0[0])*e2_solid[0] +
                           (xsp0[1]-dtoi_x0[1])*e2_solid[1] +
                           (xsp0[2]-dtoi_x0[2])*e2_solid[2],dtoi_delta);
        xy1[1] = getOddInt((xsp1[0]-dtoi_x0[0])*e2_solid[0] +
                           (xsp1[1]-dtoi_x0[1])*e2_solid[1] +
                           (xsp1[2]-dtoi_x0[2])*e2_solid[2],dtoi_delta);
        xy2[1] = getOddInt((xsp2[0]-dtoi_x0[0])*e2_solid[0] +
                           (xsp2[1]-dtoi_x0[1])*e2_solid[1] +
                           (xsp2[2]-dtoi_x0[2])*e2_solid[2],dtoi_delta);
        const int bbmax1 = max(xy0[1],max(xy1[1],xy2[1]));
        if (bbmax1 < bbmin_adt2d_solid[1])
          continue;
        const int bbmin1 = min(xy0[1],min(xy1[1],xy2[1]));
        if (bbmin1 > bbmax_adt2d_solid[1])
          continue;
        // if we made it here, we need to include this ist...
        stobb_solid.push_back(ist);
        bbVec.push_back(bbmin0);
        bbVec.push_back(bbmin1);
        bbVec.push_back(bbmax0);
        bbVec.push_back(bbmax1);
      }
    }

    const int nstobb_solid = stobb_solid.size();
    int (*bbmin)[2] = new int[nstobb_solid][2];
    int (*bbmax)[2] = new int[nstobb_solid][2];
    for (int ibb = 0; ibb < nstobb_solid; ++ibb) {
      FOR_I2 bbmin[ibb][i] = bbVec[ibb*4+i];
      FOR_I2 bbmax[ibb][i] = bbVec[ibb*4+2+i];
    }
    adt2d_solid = new Adt2d<int>(nstobb_solid,bbmin,bbmax);
    delete[] bbmin;
    delete[] bbmax;
  }

  bool Part::isInsideSolid(const double xp[3],const bool first) {

    // the first part treats the case of zero intersections differently
    // then later parts. If the first part does NOT find an intersection,
    // the point is considered inside the solid, while if a later part does
    // NOT return an intersection, the point is NOT considered inside. This
    // should keep all points within the domain of the first part specified.

    switch(solid_type) {

    case SPHERE_SOLID:
      {
        // =======================================================
        // the sphere has the following in data_solid:
        // x[3] : ddata_solid[0,1,2]
        // r    : ddata_solid[3];
        // =======================================================
        const double * const x = ddata_solid_;
        const double r = ddata_solid_[3];
        if (DIST2(xp,x) > r*r) return false;
        else return true;
      }
      break;
    case CYLINDER_SOLID:
      {
        // =======================================================
        // the cylinder has the following in data_solid:
        // x0[3] : ddata_solid[0,1,2]
        // x1[3] : ddata_solid[3,4,5]
        // r    : ddata_solid[6]
        // =======================================================
        const double * const x0 = ddata_solid_+0;
        const double * const x1 = ddata_solid_+3;
        const double dx[3] = DIFF(x1,x0);
        const double dxp[3] = DIFF(xp,x0);
        const double dp = DOT_PRODUCT(dxp,dx);
        if (dp < 0.0) {
          return true;
        }
        else {
          const double dx_mag2 = DOT_PRODUCT(dx,dx);
          if (dp > dx_mag2) {
            return true;
          }
          else {
            const double r = ddata_solid_[6];
            const double r2 = DOT_PRODUCT(dxp,dxp) - dp*dp/dx_mag2;
            if (r2 > r*r) return true;
          }
        }
        return false;
      }
      break;
    case ANNULAR_SOLID:
      {
        // =======================================================
        // the annular cylinder has the following in data_solid:
        // x0[3] : ddata_solid[0,1,2]
        // x1[3] : ddata_solid[3,4,5]
        // r0    : ddata_solid[6]
        // r1    : ddata_solid[7];
        // =======================================================
        const double * const x0 = ddata_solid_+0;
        const double * const x1 = ddata_solid_+3;
        const double dx[3] = DIFF(x1,x0);
        const double dxp[3] = DIFF(xp,x0);
        const double dp = DOT_PRODUCT(dxp,dx);
        if (dp < 0.0) {
          return true;
        }
        else {
          const double dx_mag2 = DOT_PRODUCT(dx,dx);
          if (dp > dx_mag2) {
            return true;
          }
          else {
            const double r0 = ddata_solid_[6];
            const double r1 = ddata_solid_[7];
            const double r2 = DOT_PRODUCT(dxp,dxp) - dp*dp/dx_mag2;
            if ((r2 < r0*r0)||(r2 > r1*r1)) return true;
          }
        }
        return false;
      }
      break;
    case BOX_SOLID:
      {
        // =======================================================
        // the box has the following in data_solid:
        // x0[3] : ddata_solid[0,1,2]
        // x1[3] : ddata_solid[3,4,5]
        // =======================================================
        FOR_I3 {
          if ((xp[i] < ddata_solid_[i]) || (xp[i] > ddata_solid_[3+i])) return true;
        }
        return false;
      }
      break;
    case SURFACE_SOLID:
      {
        // get our point as an even int...
        assert(b_dtoi);
        assert(adt2d_solid);

        int xy[2];
        xy[0] = getEvenInt((xp[0]-dtoi_x0[0])*e1_solid[0] +
                           (xp[1]-dtoi_x0[1])*e1_solid[1] +
                           (xp[2]-dtoi_x0[2])*e1_solid[2],dtoi_delta);
        xy[1] = getEvenInt((xp[0]-dtoi_x0[0])*e2_solid[0] +
                           (xp[1]-dtoi_x0[1])*e2_solid[1] +
                           (xp[2]-dtoi_x0[2])*e2_solid[2],dtoi_delta);
        // any point we check must be inside the current bbmin/max_adt2d_solid:
        // if not, then the prepare was not called or did not do its job...
        assert((xy[0] >= bbmin_adt2d_solid[0])&&(xy[0] <= bbmax_adt2d_solid[0]));
        assert((xy[1] >= bbmin_adt2d_solid[1])&&(xy[1] <= bbmax_adt2d_solid[1]));
        vector<int> bboxVec;
        adt2d_solid->buildListForPoint(bboxVec,xy);
        int closest = 0;
        double dist_closest;
        for (int ibb = 0; ibb < bboxVec.size(); ++ibb) {
          const int ist = stobb_solid[bboxVec[ibb]];
          // and get the tri corners as odd int's...
          const double * const xsp0 = surface->xsp[surface->spost[ist][0]];
          int xy0[2];
          xy0[0] = getOddInt((xsp0[0]-dtoi_x0[0])*e1_solid[0] + (xsp0[1]-dtoi_x0[1])*e1_solid[1] + (xsp0[2]-dtoi_x0[2])*e1_solid[2],dtoi_delta);
          xy0[1] = getOddInt((xsp0[0]-dtoi_x0[0])*e2_solid[0] + (xsp0[1]-dtoi_x0[1])*e2_solid[1] + (xsp0[2]-dtoi_x0[2])*e2_solid[2],dtoi_delta);
          const double * const xsp1 = surface->xsp[surface->spost[ist][1]];
          int xy1[2];
          xy1[0] = getOddInt((xsp1[0]-dtoi_x0[0])*e1_solid[0] + (xsp1[1]-dtoi_x0[1])*e1_solid[1] + (xsp1[2]-dtoi_x0[2])*e1_solid[2],dtoi_delta);
          xy1[1] = getOddInt((xsp1[0]-dtoi_x0[0])*e2_solid[0] + (xsp1[1]-dtoi_x0[1])*e2_solid[1] + (xsp1[2]-dtoi_x0[2])*e2_solid[2],dtoi_delta);
          const double * const xsp2 = surface->xsp[surface->spost[ist][2]];
          int xy2[2];
          xy2[0] = getOddInt((xsp2[0]-dtoi_x0[0])*e1_solid[0] + (xsp2[1]-dtoi_x0[1])*e1_solid[1] + (xsp2[2]-dtoi_x0[2])*e1_solid[2],dtoi_delta);
          xy2[1] = getOddInt((xsp2[0]-dtoi_x0[0])*e2_solid[0] + (xsp2[1]-dtoi_x0[1])*e2_solid[1] + (xsp2[2]-dtoi_x0[2])*e2_solid[2],dtoi_delta);
          const int8 dxy0[2] = { xy0[0]-xy[0], xy0[1]-xy[1] };
          const int8 dxy1[2] = { xy1[0]-xy[0], xy1[1]-xy[1] };
          const int8 dxy2[2] = { xy2[0]-xy[0], xy2[1]-xy[1] };
          const int8 A0 = dxy1[0]*dxy2[1] - dxy1[1]*dxy2[0];
          const int8 A1 = dxy2[0]*dxy0[1] - dxy2[1]*dxy0[0];
          const int8 A2 = dxy0[0]*dxy1[1] - dxy0[1]*dxy1[0];
          // skip slits...
          if (A0+A1+A2 == 0)
            continue;
          if ((A0 >= 0)&&(A1 >= 0)&&(A2 >= 0)) {
            //assert(A0+A1+A2 > 0); // check for slit
            double dist = 0.0;
            FOR_I3 {
              const double dxi = (double(A0)*(xsp0[i]-xp[i]) + double(A1)*(xsp1[i]-xp[i]) + double(A2)*(xsp2[i]-xp[i]))/double(A0+A1+A2);
              dist += dxi*e0_solid[i];
            }
            if ((closest == 0)||(fabs(dist) < fabs(dist_closest))) {
              closest = 1; // i.e. Areas positive
              dist_closest = dist;
            }
          }
          else if ((A0 <= 0)&&(A1 <= 0)&&(A2 <= 0)) {
            //assert(A0+A1+A2 < 0);
            double dist = 0.0;
            FOR_I3 {
              const double dxi = (double(A0)*(xsp0[i]-xp[i]) + double(A1)*(xsp1[i]-xp[i]) + double(A2)*(xsp2[i]-xp[i]))/double(A0+A1+A2);
              dist += dxi*e0_solid[i];
            }
            if ((closest == 0)||(fabs(dist) < fabs(dist_closest))) {
              closest = -1; // i.e. Areas positive
              dist_closest = dist;
            }
          }
        }
        if ( (first&&(closest == 0)) || ((closest == 1)&&(dist_closest < 0.0)) || ((closest == -1)&&(dist_closest > 0.0)) ) return true;
        else return false;
      }

    case NO_SOLID:
      return false;

    default:

      if (mpi_rank == 0) cout << "Warning: Part::isInsideSolid for name \"" << name << "\": solid_type not handled: " << solid_type << " " << getSolidTypeName() << endl;
      assert(0);

    }

    return false;

  }

  void Part::constrainSmoothing(double dxp[3],const double xp[3]) const {

    // the point at location xp[3] is about to be smoothed to new
    // position xp[3] + dxp[3]. If you want to constrain this motion based on
    // something about the part, do so here. The default constraint is
    // to set all components of dxp[3] to zero (see below).

    switch(ff_type) {
    case ANNULAR_CYLINDER_FF:
      {
        // get the unit radial vector at "x"...
        const double * const x0 = ddata_ff_+0;
        const double * const x1 = ddata_ff_+3;
        const double dx01[3] = DIFF(x1,x0);
        double dx[3] = DIFF(xp,x0);
        double dp = DOT_PRODUCT(dx,dx01);
        double mag2 = DOT_PRODUCT(dx01,dx01);
        FOR_I3 dx[i] -= dp*dx01[i]/mag2;
        // dx should now be a radial vector...
        mag2 = DOT_PRODUCT(dx,dx);
        assert(mag2 > 0.0);
        // now project the proposed lloyd motion onto this radial vector...
        dp = DOT_PRODUCT(dxp,dx);
        // and constrain the passed dxp to this radial motion...
        FOR_I3 dxp[i] = dp*dx[i]/mag2;
      }
      break;
    default:
      // the default constraint is no motion of part points...
      FOR_I3 dxp[i] = 0.0;
    }

  }

  void buildOrRebuildZoneVec() {

    if (mpi_rank == 0) cout << "buildOrRebuildZoneVec()" << endl;

    zoneVec.clear();
    zoneMap.clear();
    partSurfaceZoneVec.clear();
    paozn.clear();
    pzozn.clear();

    for (int ipart = 0; ipart < partVec.size(); ++ipart) {
      if (partVec[ipart]->surface) {
        if (partVec[ipart]->znosz == NULL) partVec[ipart]->znosz = new int[partVec[ipart]->surface->zoneVec.size()];
        for (int isz = 0; isz < partVec[ipart]->surface->zoneVec.size(); ++isz) {
          const int izn = zoneVec.size();
          partVec[ipart]->znosz[isz] = izn;
          string name;
          if (partVec[ipart]->b_name) {
            // if part has been named, then use "." to separate part and zone names...
            name = partVec[ipart]->name+"."+partVec[ipart]->surface->zoneVec[isz].getName();
          }
          else {
            name = partVec[ipart]->surface->zoneVec[isz].getName();
          }
          zoneMap[name] = izn;
          partSurfaceZoneVec.push_back(pair<int,int>(ipart,isz));
          zoneVec.push_back(StZone(partVec[ipart]->surface->zoneVec[isz]));
          zoneVec.back().setName(name);
          paozn.push_back(ipart);
	  pzozn.push_back(pair<int,int>(ipart,isz));
        }
      }
    }

    assert(zoneVec.size() == zoneMap.size());

  }

  int getZoneIndex(const string& name) {

    map<const string,int>::iterator iter = zoneMap.find(name);
    if (iter == zoneMap.end())
      return -1;
    return iter->second;

  }

  pair<int,int> getPartSurfaceZonePair(const int izn) {

    assert((izn >= 0)&&(izn < partSurfaceZoneVec.size()));
    return partSurfaceZoneVec[izn];

  }

  pair<int,int> getPartSurfaceZonePair(const string& name) {

    const int izn = getZoneIndex(name);
    assert((izn >= 0)&&(izn < partSurfaceZoneVec.size()));
    return partSurfaceZoneVec[izn];

  }

  void clear() {

    for (int ipart = 0, npart = partVec.size(); ipart < npart; ++ipart)
      delete partVec[ipart];
    partVec.clear();

    if (hcpPts) {
      delete hcpPts;
      hcpPts = NULL;
    }


  }

  void calcBbminmax() {
    // the bbminmax of this collection of parts is simply the
    // minmax of all surfaces...
    if (mpi_rank == 0) cout << "calcBbminmax(): partVec.size(): " << partVec.size() << endl;
    assert(!b_bbminmax);
    FOR_I6 bbminmax[i] = HUGE_VAL;
    double surface_bbminmax[6];
    for (int ipart = 0; ipart < partVec.size(); ++ipart) {
      if (partVec[ipart]->surface) {
        partVec[ipart]->surface->getBbminmax(surface_bbminmax);
        FOR_I6 bbminmax[i] = min(bbminmax[i],surface_bbminmax[i]);
      }
    }
    if (mpi_rank == 0) cout << " " << COUT_VEC(bbminmax) << " " << COUT_VEC(bbminmax+3) << endl;
    b_bbminmax = true;
  }

  void ensureBbminmax() {
    if (!b_bbminmax) calcBbminmax();
  }

  double getBoundingBoxRmax() {
    if (!b_bbminmax) calcBbminmax();
    // note that dx = -bbminmax[3]-bbminmax[0], but since it is squared, we
    // can flip the sign...
    return 0.5*sqrt((bbminmax[3]+bbminmax[0])*(bbminmax[3]+bbminmax[0]) +
                    (bbminmax[4]+bbminmax[1])*(bbminmax[4]+bbminmax[1]) +
                    (bbminmax[5]+bbminmax[2])*(bbminmax[5]+bbminmax[2]));
  }

  void getBoundingBoxCenter(double bBoxCenter[3]) {
    if (!b_bbminmax) calcBbminmax();
    // recall bbminmax[3,4,5] contains the max as negative...
    bBoxCenter[0] = 0.5*(bbminmax[0]-bbminmax[3]);
    bBoxCenter[1] = 0.5*(bbminmax[1]-bbminmax[4]);
    bBoxCenter[2] = 0.5*(bbminmax[2]-bbminmax[5]);
  }

  void getBbox(double bbmin[3],double bbmax[3]) {
    if (!b_bbminmax) calcBbminmax();
    FOR_I3 bbmin[i] = bbminmax[i];
    FOR_I3 bbmax[i] = -bbminmax[i+3];
  }

  void getBboxForZoneids(double bbmin[3],double bbmax[3],const vector<int>& zoneidVec) {

    // recall that zoneids and surface zoneids are different. zoneids refer to the
    // single list of integers that uniquely define all the parts and surface zones.

    // step 1: turn the set of zoneids into a vector of ipart,isz pairs...
    vector <pair<int,int> > partSurfaceZonePairVec(zoneidVec.size());
    for (int ii = 0; ii < zoneidVec.size(); ++ii) {
      partSurfaceZonePairVec[ii] = getPartSurfaceZonePair(zoneidVec[ii]);
    }

    FOR_I6 bbminmax[i] = HUGE_VAL;
    double surface_bbminmax[6];
    for (int ipart = 0; ipart < partVec.size(); ++ipart) {
      if (partVec[ipart]->surface) {
        // are there any isz's for this part in partSurfaceZonePairVec?...
        vector<int> surfaceZoneidVec;
        for (int ii = 0; ii < partSurfaceZonePairVec.size(); ++ii) {
          if (partSurfaceZonePairVec[ii].first == ipart) {
            surfaceZoneidVec.push_back(partSurfaceZonePairVec[ii].second);
          }
        }
        if (!surfaceZoneidVec.empty()) {
          partVec[ipart]->surface->getBboxForZoneids(surface_bbminmax,surfaceZoneidVec);
          FOR_I6 bbminmax[i] = min(bbminmax[i],surface_bbminmax[i]);
        }
      }
    }

    // and return...
    FOR_I3 bbmin[i] =  bbminmax[i  ];
    FOR_I3 bbmax[i] = -bbminmax[i+3];
  }

  /*
    void getSubzoneBoundingBoxCenterAndDiagonal(double bbmin[3],double bbmax[3],double xbb[3],double &diag,const int * const zone_flag) {

    MPI_Pause("getSubzoneBoundingBoxCenterAndDiagonal");

    }
  */

  // ----------------------------------------------------------------------
  // part output - compresses surfaces and points and writes out a part...
  //               for use in moving solver, for example...
  // ----------------------------------------------------------------------

  void writePart(const string& filename) {

    if (mpi_rank == 0) cout << "PartData::writePart " << filename << "..." << endl;

    // this routine FLATTENS whatever you have done into a single part...

    // the binary part file contains the following:
    // 1. surface
    // 2. ff_surface
    // 3. pts


    // create tmp file to write too...
    MiscUtils::mkdir_for_file_collective(filename,0);
    string tmp_filename = MiscUtils::makeTmpPrefix(filename);

    char dummy[128];
    sprintf(dummy,"%s",tmp_filename.c_str());
    MPI_File fh;
    MPI_Offset offset = 0;
    MPI_File_open(mpi_comm,dummy,MPI_MODE_WRONLY|MPI_MODE_CREATE,MPI_INFO_NULL,&fh);

    // ===============================================================
    // header
    // ===============================================================

    int itmp[5] = { PART_IO_MAGIC_NUMBER, PART_IO_VERSION, 0, 0, 0 };
    for (int ipart = 0; ipart < partVec.size(); ++ipart) {
      if (partVec[ipart]->surface) {
        // look for the keyword FF...
        // TODO: for now, just assume all surface...
        itmp[2] = 1;
      }
      if (partVec[ipart]->ff_surface) {
        // look for the keyword FF...
        // TODO: for now, just assume all surface...
        itmp[3] = 1;
      }
      if (partVec[ipart]->pts) {
        // look for the keyword FF...
        // TODO: for now, just assume all surface...
        itmp[4] = 1;
      }
      if (hcpPts) itmp[4] = 1;
    }
    assert(itmp[2] == 1); // surface
    assert(itmp[3] == 0); // ff_surface
    assert(itmp[4] == 1); // pts (any pts, not just hcpPts)

    if (mpi_rank == 0) {

      cout << " > surface: " << itmp[2] << ", ff_surface: " << itmp[3] << ", pts: " << itmp[4] << endl;

      MPI_File_write(fh,itmp,5,MPI_INT,MPI_STATUS_IGNORE);
      offset += int_size*5;

      if (itmp[2] == 1) {

        // ===============================================================
        // surface
        // ===============================================================

        int surface_itmp[3] = { 0, 0, 0 };
        for (int ipart = 0; ipart < partVec.size(); ++ipart) {
          if (partVec[ipart]->surface) {
            // HACK: to get going, assume only ipart 0 is written...
            if (ipart != 0) {
              cout << "Warning: assuming part: " << partVec[ipart]->name << " is FF" << endl;
              continue;
            }
            SurfaceShm * surface = partVec[ipart]->surface;
            // TODO: here we should separate out anything FF?...
            // for now, just take it all...
            surface_itmp[0] += surface->zoneVec.size();
            surface_itmp[1] += surface->nsp;
            surface_itmp[2] += surface->nst;
          }
        }

        cout << " > surface nzones: " << surface_itmp[0] << ", nsp: " << surface_itmp[1] << ", nst: " << surface_itmp[2] << endl;

        MPI_File_write(fh,surface_itmp,3,MPI_INT,MPI_STATUS_IGNORE);
        offset += int_size*3;

        // compressed surface zones...

        int nzone = 0;
        for (int ipart = 0; ipart < partVec.size(); ++ipart) {
          if (partVec[ipart]->surface) {
            // HACK: to get going, assume only ipart 0 is written...
            if (ipart != 0) {
              if (mpi_rank == 0) cout << "Warning: assuming part: " << partVec[ipart]->name << " is FF" << endl;
              continue;
            }
            SurfaceShm * surface = partVec[ipart]->surface;
            for (int izone = 0; izone < surface->zoneVec.size(); ++izone) {
              int length = surface->zoneVec[izone].getName().length();
              MPI_File_write(fh,&length,1,MPI_INT,MPI_STATUS_IGNORE);
              offset += int_size;
              MPI_File_write(fh,(char*)surface->zoneVec[izone].getName().c_str(),length,MPI_CHAR,MPI_STATUS_IGNORE);
              offset += length;
              cout << " > zone " << nzone << " name \"" << surface->zoneVec[izone].getName() << "\"" << endl;
              ++nzone;
            }
          }
        }
        assert(nzone == surface_itmp[0]);

        // compressed surface xsp...
        // we could do this in parallel, but serial for now...

        int nsp = 0;
        for (int ipart = 0; ipart < partVec.size(); ++ipart) {
          if (partVec[ipart]->surface) {
            // HACK: to get going, assume only ipart 0 is written...
            if (ipart != 0) {
              if (mpi_rank == 0) cout << "Warning: assuming part: " << partVec[ipart]->name << " is FF" << endl;
              continue;
            }
            SurfaceShm * surface = partVec[ipart]->surface;
            MPI_File_write(fh,(double*)surface->xsp,surface->nsp*3,MPI_DOUBLE,MPI_STATUS_IGNORE);
            nsp += surface->nsp;
          }
        }
        assert(nsp == surface_itmp[1]);
        offset += double_size*nsp*3;

        // compressed surface spost...
        int nst = 0;
        for (int ipart = 0; ipart < partVec.size(); ++ipart) {
          if (partVec[ipart]->surface) {
            // HACK: to get going, assume only ipart 0 is written...
            if (ipart != 0) {
              if (mpi_rank == 0) cout << "Warning: assuming part: " << partVec[ipart]->name << " is FF" << endl;
              continue;
            }
            SurfaceShm * surface = partVec[ipart]->surface;
            MPI_File_write(fh,(int*)surface->spost,surface->nst*3,MPI_INT,MPI_STATUS_IGNORE);
            nst += surface->nst;
          }
        }
        assert(nst == surface_itmp[2]);
        offset += int_size*nst*3;

        nst = 0;
        for (int ipart = 0; ipart < partVec.size(); ++ipart) {
          if (partVec[ipart]->surface) {
            // HACK: to get going, assume only ipart 0 is written...
            if (ipart != 0) {
              if (mpi_rank == 0) cout << "Warning: assuming part: " << partVec[ipart]->name << " is FF" << endl;
              continue;
            }
            SurfaceShm * surface = partVec[ipart]->surface;
            MPI_File_write(fh,(int*)surface->znost,surface->nst,MPI_INT,MPI_STATUS_IGNORE);
            nst += surface->nst;
          }
        }
        assert(nst == surface_itmp[2]);
        offset += int_size*nst;

      }

      if (itmp[3] == 1) {

        // ===============================================================
        // ff_surface
        // ===============================================================

        assert(0);

      }

    }

    MPI_Bcast(&offset,1,MPI_INT8,0,mpi_comm);

    // ===============================
    // finally points...
    // ===============================

    if (itmp[4] == 1) {

      // now all ranks count (and then write) all their local points...

      int my_np = 0;
      for (int ipart = 0; ipart < partVec.size(); ++ipart) {
        if (partVec[ipart]->pts) {
          my_np += partVec[ipart]->pts->np;
        }
      }
      if (hcpPts) my_np += hcpPts->np;

      int my_disp;
      MPI_Scan(&my_np,&my_disp,1,MPI_INT,MPI_SUM,mpi_comm);

      int np_global = my_disp;
      MPI_Bcast(&np_global,1,MPI_INT,mpi_size-1,mpi_comm);
      assert(np_global > 0);

      if (mpi_rank == 0) {
        cout << " > total point count: " << np_global << endl;
        MPI_File_write(fh,&np_global,1,MPI_INT,MPI_STATUS_IGNORE);
      }
      offset += int_size;

      // now xp...

      double * dbuf = new double[my_np*3];
      int my_np_check = 0;
      for (int ipart = 0; ipart < partVec.size(); ++ipart) {
        if (partVec[ipart]->pts) {
          for (int ip = 0; ip < partVec[ipart]->pts->np; ++ip) {
            FOR_I3 dbuf[my_np_check*3+i] = partVec[ipart]->pts->xp[ip][i];
            ++my_np_check;
          }
        }
      }
      if (hcpPts) {
        for (int ip = 0; ip < hcpPts->np; ++ip) {
          FOR_I3 dbuf[my_np_check*3+i] = hcpPts->xp[ip][i];
          ++my_np_check;
        }
      }
      assert(my_np_check == my_np);

      writeChunkedData<double>(fh,offset+double_size*(my_disp-my_np)*3,dbuf,my_np*3,mpi_comm);
      offset += double_size*np_global*3;

      // delta...

      my_np_check = 0;
      for (int ipart = 0; ipart < partVec.size(); ++ipart) {
        if (partVec[ipart]->pts) {
          // look for the keyword FF...
          // TODO: for now, just assume all surface...
          for (int ip = 0; ip < partVec[ipart]->pts->np; ++ip) {
            dbuf[my_np_check++] = partVec[ipart]->pts->delta[ip];
          }
        }
      }
      if (hcpPts) {
        for (int ip = 0; ip < hcpPts->np; ++ip) {
          dbuf[my_np_check++] = hcpPts->delta[ip];
        }
      }
      assert(my_np_check == my_np);

      writeChunkedData<double>(fh,offset+double_size*(my_disp-my_np),dbuf,my_np,mpi_comm);
      offset += double_size*np_global;

      delete[] dbuf;

    }

    MPI_File_set_size(fh,offset);
    MPI_File_close(&fh);

    // if we successfully wrote the tmp file, rename/replace the file with it...
    
    if (mpi_rank == 0) {
      remove(filename.c_str());
      rename(tmp_filename.c_str(),filename.c_str());
    }

  }

  bool checkCoplanar(const int ipart0,const int ist0,const int bits0,
                     const int ipart1,const int ist1,const int bits1,const double degrees) {
    // exacly 180 or greater leads to no combining of facets...
    if (degrees >= 180.0)
      return false;
    // degrees is the crease angle. using an angle close to 180 (e.g. 179) makes
    // the check very stringent. I recommend using something like 150...
    const double dp_tol = cos((180.0-degrees)*M_PI/180.0);
    assert(ipart0 == ipart1); // assume we only call this within a single part...
    assert(partVec[ipart0]->surface);
    assert((ist0 >= 0)&&(ist0 < partVec[ipart0]->surface->nst));
    assert((bits0 >= 0)&&(bits0 < (1<<6)));
    assert(partVec[ipart1]->surface);
    assert((ist1 >= 0)&&(ist1 < partVec[ipart1]->surface->nst));
    assert((bits1 >= 0)&&(bits1 < (1<<6)));
    double n0[3] = TRI_NORMAL_2(partVec[ipart0]->surface->xsp[partVec[ipart0]->surface->spost[ist0][0]],
                                partVec[ipart0]->surface->xsp[partVec[ipart0]->surface->spost[ist0][1]],
                                partVec[ipart0]->surface->xsp[partVec[ipart0]->surface->spost[ist0][2]]);
    if (bits0) PeriodicData::periodicRotate(n0,1,bits0);
    double n1[3] = TRI_NORMAL_2(partVec[ipart1]->surface->xsp[partVec[ipart1]->surface->spost[ist1][0]],
                                partVec[ipart1]->surface->xsp[partVec[ipart1]->surface->spost[ist1][1]],
                                partVec[ipart1]->surface->xsp[partVec[ipart1]->surface->spost[ist1][2]]);
    if (bits1) PeriodicData::periodicRotate(n1,1,bits1);
    const double n0_mag = sqrt(DOT_PRODUCT(n0,n0));
    const double n1_mag = sqrt(DOT_PRODUCT(n1,n1));
    if ((n0_mag <= 0.0)||(n1_mag <= 0.0))
      return true; // linear tris return true
    return( DOT_PRODUCT(n0,n1)/(n0_mag*n1_mag) > dp_tol );
  }

  void writeTriDebugFile(const string& filename,const set<pair<int,int> >& partIstSet) {

    // one rank for now...
    assert(mpi_rank == 0);
    cout << " > writeTriDebugFile: " << filename << endl;

    // build the node map...
    map<const pair<int,int>,int> nodeMap;
    int nsp = 0;
    for (set<pair<int,int> >::iterator it = partIstSet.begin(); it != partIstSet.end(); ++it) {
      const int ipart = it->first;
      const int ist = it->second;
      FOR_I3 {
        const int isp = partVec[ipart]->surface->spost[ist][i];
        map<const pair<int,int>,int>::iterator it2 = nodeMap.find(pair<int,int>(ipart,isp));
        if (it2 == nodeMap.end())
          nodeMap[pair<int,int>(ipart,isp)] = nsp++;
      }
    }
    int nst = partIstSet.size();

    FILE * fp = fopen(filename.c_str(),"w");
    assert(fp != NULL);

    fprintf(fp,"TITLE = \"%s\"\n",filename.c_str());
    fprintf(fp,"VARIABLES = \"X\"\n");
    fprintf(fp,"\"Y\"\n");
    fprintf(fp,"\"Z\"\n");

    // zone header
    fprintf(fp,"ZONE T=\"%s\"\n","blah");
    fprintf(fp,"N=%d, E=%d, F=FEPOINT, ET=TRIANGLE\n",nsp,nst);

    // points first...
    nsp = 0;
    for (map<const pair<int,int>,int>::iterator it = nodeMap.begin(); it != nodeMap.end(); ++it) {
      const int ipart = it->first.first;
      const int isp = it->first.second;
      it->second = nsp++; // re-index nodes...
      fprintf(fp,"%18.15le %18.15le %18.15le\n",
              partVec[ipart]->surface->xsp[isp][0],
              partVec[ipart]->surface->xsp[isp][1],
              partVec[ipart]->surface->xsp[isp][2]);
    }

    // then tris...
    for (set<pair<int,int> >::iterator it = partIstSet.begin(); it != partIstSet.end(); ++it) {
      const int ipart = it->first;
      const int ist = it->second;
      int noosp[3];
      FOR_I3 {
        const int isp = partVec[ipart]->surface->spost[ist][i];
        map<const pair<int,int>,int>::iterator it2 = nodeMap.find(pair<int,int>(ipart,isp));
        assert(it2 != nodeMap.end());
        noosp[i] = it2->second;
      }
      fprintf(fp,"%d %d %d\n",noosp[0]+1,noosp[1]+1,noosp[2]+1); // tecplot file is 1-indexed
    }

    fclose(fp);

  }

  void runDiagnostics(Param * param, const bool b_help) {

    if (b_help) {
      WUI(INFO," runDiagnostics options...\n"); // TODO
      return;
    }

    double b_repair = false;
    double dn = 0.0;
    int niter = 100;
    double neg_frac = 0.25;
    bool b_write_displacement = false;
    try { // force user to call repair
      int iarg = 0;
      while (iarg < param->size()) {
        string token = MiscUtils::toUpperCase(param->getString(iarg++));
        if (token == "REPAIR_INTERSECTIONS") {
          b_repair = true;
        }
        else if (token == "DN") {
          dn = param->getDouble(iarg++);
          if (dn <= 0.0)
            throw(1);
        }
        else if (token == "NITER") {
          niter = param->getInt(iarg++);
          if (niter <= 0)
            throw(1);
        }
        else if (token == "NEG_FRAC") {
          neg_frac = param->getDouble(iarg++);
          if ((neg_frac < 0.0)&&(neg_frac > 1.0))
            throw(1);
        }
        else if (token == "WRITE_DISPLACEMENT") {
          if (!b_repair) {
            WUI(WARN,"WRITE_DISPLACEMENT only meaningful when using REPAIR_INTERSECTIONS. Skipping...");
          }
          else {
            b_write_displacement = true;
          }
        }
        else {
          WUI(WARN,"unrecognized RUN_DIAGNOSTICS param " << token << ". Skipping...");
        }
      }
    }
    catch(int e) {
      WUI(WARN,"expecting RUN_DIAGNOSTICS [REPAIR DN <+double> NITER <+int>]");
      return;
    }

    // force user to call repair and dn
    if ((!b_repair)||(dn == 0.0))
      niter = 1;

    if (mpi_rank == 0) {
      cout << "runDiagnostics()" << endl;
      if (b_repair) cout << " > repair: dn " << dn << " niter " << niter << endl;
    }
    double wtime0 = MPI_Wtime();

    // to compute the displacement, store the initial positions of all sp's in
    // surface->dx_sp...

    if (b_repair) {

      for (int ipart = 0; ipart < partVec.size(); ++ipart) {
        if (partVec[ipart]->surface) {
          if (partVec[ipart]->surface->dx_sp == NULL)
            CTI_Mmap_rw(partVec[ipart]->surface->dx_sp,partVec[ipart]->surface->nsp);
          int nsp_split = partVec[ipart]->surface->nsp/mpi_size_shared;
          if (partVec[ipart]->surface->nsp%mpi_size_shared != 0) nsp_split += 1;
          const int isp0 = min(partVec[ipart]->surface->nsp,nsp_split*mpi_rank_shared);
          const int isp1 = min(partVec[ipart]->surface->nsp,nsp_split*(mpi_rank_shared+1));
          for (int isp = isp0; isp < isp1; ++isp) {
            FOR_I3 partVec[ipart]->surface->dx_sp[isp][i] = partVec[ipart]->surface->xsp[isp][i];
          }
        }
      }
      MPI_Barrier(mpi_comm_shared);

    }


    // build a local list of tris and centroids for ALL parts...

    const double angle_tol_degrees = 0.5; // degrees
    const double angle_tol = cos(M_PI * (1.0-(angle_tol_degrees/180.0)));  // radians ~= 0.57 degrees

    const int nst = getFlattenedNst();
    if (mpi_rank == 0) cout << " > flattened nst: " << nst << endl;

    const int ist_begin = int((int8(nst)*int8(mpi_rank))/int8(mpi_size)); assert((ist_begin >= 0)&&(ist_begin <= nst));
    const int ist_end = int((int8(nst)*int8(mpi_rank+1))/int8(mpi_size)); assert((ist_end >= 0)&&(ist_end <= nst));
    int my_nst = ist_end-ist_begin; assert(my_nst >= 0);

    double (*xcc)[3] = new double[my_nst][3];
    int8 * ist_flat = new int8[my_nst];
    for (int ist = ist_begin; ist < ist_end; ++ist) {
      ist_flat[ist-ist_begin] = ist;
      int ipart,ist_in_part;
      getPartStFromFlattenedSt(ipart,ist_in_part,ist);
      assert(partVec[ipart]->surface);
      FOR_I3 xcc[ist-ist_begin][i] = 0.0;
      FOR_J3 {
        const int isp = partVec[ipart]->surface->spost[ist_in_part][j];
        FOR_I3 {
          xcc[ist-ist_begin][i] += partVec[ipart]->surface->xsp[isp][i];
        }
      }
      FOR_I3 xcc[ist-ist_begin][i] /= 3.0;
    }

    if (mpi_rank == 0) cout << " > partitioning surface tris..." << endl;

    repartXcvPadt(xcc,ist_flat,my_nst,mpi_comm);
    delete[] xcc;

    /*
      int nsp_local = 0;
      double *dist_sp_local = NULL;
      uint8 *psposp_local = NULL; // part-isp of isp_local
      int (*spost_local)[3] = NULL;
      if (b_repair) {
      // also need a local array to store signed distance...
      spost_local = new int[my_nst][3];
      map<const pair<int,int>,int> partIspMap;
      for (int ii = 0; ii < my_nst; ++ii) {
      int ipart1,ist1;
      getPartStFromFlattenedSt(ipart1,ist1,ist_flat[ii]);
      FOR_J3 {
      const int isp1 = partVec[ipart1]->surface->spost[ist1][j];
      map<const pair<int,int>,int>::iterator it = partIspMap.find(pair<int,int>(ipart1,isp1));
      if (it == partIspMap.end()) {
      spost_local[ii][j] = nsp_local;
      partIspMap[pair<int,int>(ipart1,isp1)] = nsp_local++;
      }
      else {
      spost_local[ii][j] = it->second;
      }
      }
      }
      psposp_local = new uint8[nsp_local];
      for (map<const pair<int,int>,int>::iterator it = partIspMap.begin(); it != partIspMap.end(); ++it)
      psposp_local[it->second] = ((uint8(it->first.first)<<32)|uint8(it->first.second));
      dist_sp_local = new double[nsp_local];
      }
    */

    bool b_requires_repair = false;
    for (int iter = 0; iter < niter; ++iter) {

      if ((mpi_rank == 0)&&(b_repair))
        cout <<
          "===========================================================" <<
          "\n > repair iter: " << iter <<
          "\n===========================================================" << endl;

      // now each rank has a geometrically local set of tris in ist_flat.
      double bbmin_local[3] = { HUGE_VAL, HUGE_VAL, HUGE_VAL };
      double bbmax_local[3] = { -HUGE_VAL, -HUGE_VAL, -HUGE_VAL };
      for (int ii = 0; ii < my_nst; ++ii) {
        int ipart,ist_in_part;
        getPartStFromFlattenedSt(ipart,ist_in_part,ist_flat[ii]);
        FOR_J3 {
          const int isp = partVec[ipart]->surface->spost[ist_in_part][j];
          FOR_I3 {
            bbmin_local[i] = min(bbmin_local[i],partVec[ipart]->surface->xsp[isp][i]);
            bbmax_local[i] = max(bbmax_local[i],partVec[ipart]->surface->xsp[isp][i]);
          }
        }
      }

      if (b_repair) {
        //FOR_I3 bbmin_local[i] -= overlap_separation;
        //FOR_I3 bbmax_local[i] += overlap_separation;
        // point normals...
        for (int ipart = 0; ipart < partVec.size(); ++ipart) {
          if (partVec[ipart]->surface) {
            if (iter == 0)
              partVec[ipart]->surface->ensureNsp();
            else
              partVec[ipart]->surface->rebuildNsp();
          }
        }
        //for (int isp_local = 0; isp_local < nsp_local; ++isp_local) {
        //  dist_sp_local[isp_local] = overlap_separation;
        //  //cout << "setting dist_sp_local[" << isp_local << "] = " << overlap_separation << endl;
        //}
      }

      // build an adt of ALL surface tris from all parts that touch this bbox...
      vector<pair<int,int> > partIstVec;
      for (int ipart = 0; ipart < partVec.size(); ++ipart) {
        if (partVec[ipart]->surface) {
          for (int ist = 0; ist < partVec[ipart]->surface->nst; ++ist) {
            const int isp0 = partVec[ipart]->surface->spost[ist][0];
            const int isp1 = partVec[ipart]->surface->spost[ist][1];
            const int isp2 = partVec[ipart]->surface->spost[ist][2];
            if ((max(partVec[ipart]->surface->xsp[isp0][0],max(partVec[ipart]->surface->xsp[isp1][0],partVec[ipart]->surface->xsp[isp2][0])) < bbmin_local[0]) ||
                (min(partVec[ipart]->surface->xsp[isp0][0],min(partVec[ipart]->surface->xsp[isp1][0],partVec[ipart]->surface->xsp[isp2][0])) > bbmax_local[0]) ||
                (max(partVec[ipart]->surface->xsp[isp0][1],max(partVec[ipart]->surface->xsp[isp1][1],partVec[ipart]->surface->xsp[isp2][1])) < bbmin_local[1]) ||
                (min(partVec[ipart]->surface->xsp[isp0][1],min(partVec[ipart]->surface->xsp[isp1][1],partVec[ipart]->surface->xsp[isp2][1])) > bbmax_local[1]) ||
                (max(partVec[ipart]->surface->xsp[isp0][2],max(partVec[ipart]->surface->xsp[isp1][2],partVec[ipart]->surface->xsp[isp2][2])) < bbmin_local[2]) ||
                (min(partVec[ipart]->surface->xsp[isp0][2],min(partVec[ipart]->surface->xsp[isp1][2],partVec[ipart]->surface->xsp[isp2][2])) > bbmax_local[2]))
              continue;
            partIstVec.push_back(pair<int,int>(ipart,ist));
          }
        }
      }

      if (mpi_rank == 0) cout << " > building tri adt..." << endl;

      const int nbb = partIstVec.size();
      double (*bbmin)[3] = new double[nbb][3];
      double (*bbmax)[3] = new double[nbb][3];
      for (int ii = 0; ii < nbb; ++ii) {
        const int ipart = partIstVec[ii].first;
        const int ist = partIstVec[ii].second;
        FOR_I3 {
          bbmin[ii][i] = HUGE_VAL;
          bbmax[ii][i] = -HUGE_VAL;
        }
        FOR_J3 {
          const int isp = partVec[ipart]->surface->spost[ist][j];
          FOR_I3 {
            bbmin[ii][i] = min(bbmin[ii][i],partVec[ipart]->surface->xsp[isp][i]);
            bbmax[ii][i] = max(bbmax[ii][i],partVec[ipart]->surface->xsp[isp][i]);
          }
        }
      }

      Adt<double> adt(nbb,bbmin,bbmax);
      delete[] bbmin;
      delete[] bbmax;

      /*
      // this ratio is 1.2 to 1.5 or so...
      int my_ibuf[2] = { my_nst, nbb };
      int ibuf[2];
      MPI_Reduce(my_ibuf,ibuf,2,MPI_INT,MPI_SUM,0,mpi_comm);
      if (mpi_rank == 0) {
      assert(ibuf[0] == nst);
      cout << " > avg nbb/nst: " << double(ibuf[1])/double(ibuf[0]) << endl;
      }
      */

      // and check...

      set<int> edgeMatchSet;
      set<int> unhandledMatchSet;
      set<pair<int,int> > coplanarZoneSet;
      set<pair<int,int> > coplanarTriSet;
      set<pair<int,int> > intersectingZoneSet; // zone-pairs
      set<pair<int,int> > intersectingTriSet;

      double * x0[3];
      double * x1[3];
      int idata[6];
      double ddata[4];
      vector<double> doubleVec;
      int debug_tritri_index = 0;
      vector<int> intVec;
      for (int ii = 0; ii < my_nst; ++ii) {
        int ipart1,ist1;
        getPartStFromFlattenedSt(ipart1,ist1,ist_flat[ii]);
        int (*spost)[3] = partVec[ipart1]->surface->spost;
        double (*xsp)[3] = partVec[ipart1]->surface->xsp;
        const double normal1[3] = TRI_NORMAL_2(xsp[spost[ist1][0]],xsp[spost[ist1][1]],xsp[spost[ist1][2]]);
        const double area1 = MAG(normal1);
        assert(area1 > 0.0);
        FOR_I3 {
          bbmin_local[i] = HUGE_VAL;
          bbmax_local[i] = -HUGE_VAL;
        }
        FOR_J3 {
          const int isp = spost[ist1][j];
          FOR_I3 {
            bbmin_local[i] = min(bbmin_local[i],xsp[isp][i]);
            bbmax_local[i] = max(bbmax_local[i],xsp[isp][i]);
          }
        }
        //if (b_repair) {
        //  FOR_I3 bbmin_local[i] -= overlap_separation;
        //  FOR_I3 bbmax_local[i] += overlap_separation;
        //}
        int edge_count = 0;
        assert(intVec.empty());
        adt.buildListForBBox(intVec,bbmin_local,bbmax_local);
        for (int jj = 0; jj < intVec.size(); ++jj) {
          const int ipart2 = partIstVec[intVec[jj]].first;
          const int ist2 = partIstVec[intVec[jj]].second;
          // at this point, consider anyone but us. Below we limit some work further...
          if ((ipart2 != ipart1)||(ist2 != ist1)) {
            // these tris are different...
            bool b_handled = false;
            bool b_node = false;
            if (ipart2 == ipart1) {
              const int tri_match =
                (spost[ist1][0] == spost[ist2][0]) +
                (spost[ist1][0] == spost[ist2][1])*(1<<1) +
                (spost[ist1][0] == spost[ist2][2])*(1<<2) +
                (spost[ist1][1] == spost[ist2][0])*(1<<3) +
                (spost[ist1][1] == spost[ist2][1])*(1<<4) +
                (spost[ist1][1] == spost[ist2][2])*(1<<5) +
                (spost[ist1][2] == spost[ist2][0])*(1<<6) +
                (spost[ist1][2] == spost[ist2][1])*(1<<7) +
                (spost[ist1][2] == spost[ist2][2])*(1<<8);
              if (tri_match != 0) {
                switch (tri_match) {
		case 1:
		case 2:
		case 4:
		case 8:
		case 16:
		case 32:
		case 64:
		case 128:
		case 256:
		  // single node match. For now leave b_handled false and continue to
		  // intersection checking below...
		  // TODO: is there a more efficient way to check this?...
		  //b_handled = true;
		  b_node = true;
		  break;
		case 10:
		case 20:
		case 33:
		case 68:
		case 80:
		case 129:
		case 160:
		case 258:
		case 264:
		  // matching along a single edge: compute angle when ist2 > ist1...
		  b_handled = true;
		  ++edge_count; // record aligned edge count...
		  if (ist2 > ist1) {
		    const double normal2[3] = TRI_NORMAL_2(xsp[spost[ist2][0]],xsp[spost[ist2][1]],xsp[spost[ist2][2]]);
		    const double area2 = MAG(normal2);
		    const double dp = DOT_PRODUCT(normal1,normal2);
		    // when this normalized dot product is nearly -1, we are completely folded along the edge...
		    if (dp < area1*area2*angle_tol) {
		      assert(partVec[ipart1]->znosz);
		      //assert(partVec[ipart2]->znosz); // ipart1 == ipart2
		      coplanarZoneSet.insert(pair<int,int>(partVec[ipart1]->znosz[partVec[ipart1]->surface->znost[ist1]],
							   partVec[ipart2]->znosz[partVec[ipart2]->surface->znost[ist2]]));
		      coplanarTriSet.insert(pair<int,int>(ipart1,ist1));
		      coplanarTriSet.insert(pair<int,int>(ipart2,ist2));
		    }
		  }
		  break;
		default:
		  // case was not handled...
		  // write some code here and below to handle...
		  FOR_I3 {
		    FOR_J3 {
		      if ((spost[ist1][i] == spost[ist2][(j+1)%3])&&
			  (spost[ist1][(i+1)%3] == spost[ist2][j])&&
			  (spost[ist1][(i+2)%3] != spost[ist2][(j+2)%3])) {
			edgeMatchSet.insert(tri_match);
			assert(!b_handled);
			b_handled = true;
		      }
		    }
		  }
		  if (!b_handled) {
		    unhandledMatchSet.insert(tri_match);
		  }
                }
              }
            }
            /*
	      if (b_repair) {
              if (!(b_handled||b_node)) {
	      double (*xsp2)[3] = partVec[ipart2]->surface->xsp;
	      int (*spost2)[3] = partVec[ipart2]->surface->spost;
	      // loop in tri1's nodes...
	      FOR_I3 {
	      const int isp1_local = spost_local[ii][i];
	      const int isp1 = spost[ist1][i];
	      // normal is in partVec[ipart1]->surface->n_sp[isp1]...
	      double e1[3],e2[3];
	      MiscUtils::getBestE1E2FromE0(e1,e2,partVec[ipart1]->surface->n_sp[isp1]);
	      // build tri2 in 2d...
	      double xtri2d[3][2];
	      FOR_J3 {
	      const double dx[3] = DIFF(xsp2[spost2[ist2][j]],xsp[isp1]);
	      xtri2d[j][0] = DOT_PRODUCT(dx,e1);
	      xtri2d[j][1] = DOT_PRODUCT(dx,e2);
	      }
	      double cp2d[3];
	      double sum = 0.0;
	      FOR_J3 {
	      cp2d[j] = xtri2d[(j+1)%3][0]*xtri2d[(j+2)%3][1] - xtri2d[(j+1)%3][1]*xtri2d[(j+2)%3][0];
	      sum += cp2d[j];
	      }
	      double tol = 0.01*fabs(sum); // generous tol here is fine -- don't want to miss the intersection
	      if ((cp2d[0] < tol)&&(cp2d[1] < tol)&&(cp2d[2] < tol)) {
	      double xint[3] = { 0.0, 0.0, 0.0 };
	      FOR_J3 FOR_K3 xint[j] += cp2d[k]*xsp2[spost2[ist2][k]][j];
	      FOR_J3 xint[j] /= sum;
	      const double dx[3] = DIFF(xint,xsp[isp1]);
	      const double dist = MAG(dx);
	      if (dist < fabs(dist_sp_local[isp1_local])) {
	      const double sign = DOT_PRODUCT(partVec[ipart1]->surface->n_sp[isp1],dx);
	      if (sign > 0.0) {
	      dist_sp_local[isp1_local] = dist;
	      }
	      else {
	      assert(sign < 0.0);
	      dist_sp_local[isp1_local] = -dist;
	      }
	      }
	      }
	      }
              }
	      }
            */
            // if we did not find an edge nbr (aligned or flipped), then run an intersection check...
            if ((!b_handled)&&((ipart2 > ipart1)||((ipart2 == ipart1)&&(ist2 > ist1)))) {
              FOR_I3 x0[i] = xsp[spost[ist1][i]];
              FOR_I3 x1[i] = partVec[ipart2]->surface->xsp[partVec[ipart2]->surface->spost[ist2][i]];
              if (const int n = MiscUtils::calcTriTriIntersection(idata,ddata,x0,x1)) {
                if (n == -2) {
                  // this case is not recognized by the TriTri routine -- save it to gen some new code...
                  char filename[128]; sprintf(filename,"debug/tritri.%06d.%04d.bin",mpi_rank,debug_tritri_index++);
                  MiscUtils::mkdir_for_file(filename);
                  MiscUtils::writeTriTriBin(filename,x0,x1);
                }
                else {
                  if (n == 1) {
                    // a single intersection can be a node-node, which is fine...
                    if ((idata[0] != MiscUtils::NODE_NODE_INT)||(ipart1 != ipart2)||(spost[ist1][idata[1]] != spost[ist2][idata[2]])) {
                      // not node-node...
                      double xi[3]; MiscUtils::getTriTriX(xi,x0,x1,idata[0],idata+1,ddata);
                      doubleVec.push_back(xi[0]);
                      doubleVec.push_back(xi[1]);
                      doubleVec.push_back(xi[2]);
                      assert(partVec[ipart1]->znosz);
                      assert(partVec[ipart2]->znosz);
                      intersectingZoneSet.insert(pair<int,int>(partVec[ipart1]->znosz[partVec[ipart1]->surface->znost[ist1]],
							       partVec[ipart2]->znosz[partVec[ipart2]->surface->znost[ist2]]));
                      intersectingTriSet.insert(pair<int,int>(ipart1,ist1));
                      intersectingTriSet.insert(pair<int,int>(ipart2,ist2));
                    }
                  }
                  else if (n == 2) {
                    // a double intersection is the most common between intersecting tris. We know this is
                    // not along a shared edge with common nodes, because this case is checked for above...
                    double xi0[3]; MiscUtils::getTriTriX(xi0,x0,x1,idata[0],idata+1,ddata);
                    double xi1[3]; MiscUtils::getTriTriX(xi1,x0,x1,idata[3],idata+4,ddata+2);
                    doubleVec.push_back(xi0[0]);
                    doubleVec.push_back(xi0[1]);
                    doubleVec.push_back(xi0[2]);
                    doubleVec.push_back(xi1[0]);
                    doubleVec.push_back(xi1[1]);
                    doubleVec.push_back(xi1[2]);
                    // store the zones...
                    assert(partVec[ipart1]->znosz);
                    assert(partVec[ipart2]->znosz);
                    intersectingZoneSet.insert(pair<int,int>(partVec[ipart1]->znosz[partVec[ipart1]->surface->znost[ist1]],
							     partVec[ipart2]->znosz[partVec[ipart2]->surface->znost[ist2]]));
                    intersectingTriSet.insert(pair<int,int>(ipart1,ist1));
                    intersectingTriSet.insert(pair<int,int>(ipart2,ist2));
                  }
                  else {
                    assert(n == -1); // coplanar tris
                    //MiscUtils::writeTriTri(0,x0,x1);
                  }
                }
              }
            }
          }
        }
        assert(edge_count == 3);
        intVec.clear();
      }

      /*
	if (b_repair) {
        uint8* send_buf_uint8 = NULL;
        double* send_buf_double = NULL;
        int send_count;
        for (int iter = 0; iter < 2; ++iter) {
	send_count = 0;
	for (int isp_local = 0; isp_local < nsp_local; ++isp_local) {
	if (dist_sp_local[isp_local] < overlap_separation) {
	if (iter == 1) {
	send_buf_uint8[send_count] = psposp_local[isp_local];
	send_buf_double[send_count] = dist_sp_local[isp_local];
	}
	++send_count;
	}
	}
	if (iter == 0) {
	send_buf_uint8 = new uint8[send_count];
	send_buf_double = new double[send_count];
	}
        }

        int * nora = NULL;
        if (mpi_rank == 0) nora = new int[mpi_size];
        MPI_Gather(&send_count,1,MPI_INT,nora,1,MPI_INT,0,mpi_comm);

        int * disp = NULL;
        uint8 * recv_buf_uint8 = NULL;
        double * recv_buf_double = NULL;
        int nrecv = 0;
        if (mpi_rank == 0) {
	disp = new int[mpi_size];
	disp[0] = 0;
	for (int rank = 1; rank < mpi_size; ++rank)
	disp[rank] = disp[rank-1]+nora[rank-1];
	nrecv = disp[mpi_size-1]+nora[mpi_size-1];
	recv_buf_uint8 = new uint8[nrecv];
	recv_buf_double = new double[nrecv];
        }
        MPI_Gatherv(send_buf_uint8,send_count,MPI_UINT8,recv_buf_uint8,nora,disp,MPI_UINT8,0,mpi_comm);
        delete[] send_buf_uint8;
        MPI_Gatherv(send_buf_double,send_count,MPI_DOUBLE,recv_buf_double,nora,disp,MPI_DOUBLE,0,mpi_comm);
        delete[] send_buf_double;
        delete[] nora;
        delete[] disp;

        for (int ipart = 0; ipart < partVec.size(); ++ipart) {
	if (partVec[ipart]->surface) {
	if (partVec[ipart]->surface->dist_sp == NULL)
	CTI_Mmap(partVec[ipart]->surface->dist_sp,partVec[ipart]->surface->nsp);
	// we also need flag_sp to record intersections...
	if (partVec[ipart]->surface->flag_sp == NULL)
	CTI_Mmap(partVec[ipart]->surface->flag_sp,partVec[ipart]->surface->nsp);
	if (mpi_rank == 0) {
	for (int isp = 0; isp < partVec[ipart]->surface->nsp; ++isp) {
	partVec[ipart]->surface->dist_sp[isp] = overlap_separation;
	partVec[ipart]->surface->flag_sp[isp] = 0;
	}
	}
	}
        }

        // perform reduction on rank 0 and broadcast across nodes...
        if (mpi_rank == 0) {
	double dist_min = HUGE_VAL;
	double dist_max = -HUGE_VAL;
	for (int irecv = 0; irecv < nrecv; ++irecv) {
	const int ipart = int(recv_buf_uint8[irecv]>>32);
	const int isp   = int(recv_buf_uint8[irecv]&MASK_32BITS);
	assert((ipart >= 0)&&(ipart < int(partVec.size())));
	assert((isp >= 0)&&(isp < partVec[ipart]->surface->nsp));
	partVec[ipart]->surface->dist_sp[isp] = min(partVec[ipart]->surface->dist_sp[isp],recv_buf_double[irecv]);
	dist_min = min(dist_min,recv_buf_double[irecv]);
	dist_max = max(dist_max,recv_buf_double[irecv]);
	}
	cout << " > dist_min, dist_max, overlap_separation: " << dist_min << " " << dist_max << " " << overlap_separation << endl;
	delete[] recv_buf_uint8;
	delete[] recv_buf_double;
        }
        if (mpi_rank_shared == 0) {
	for (int ipart = 0; ipart < partVec.size(); ++ipart) {
	if (partVec[ipart]->surface) {
	MPI_Bcast(partVec[ipart]->surface->dist_sp,partVec[ipart]->surface->nsp,MPI_DOUBLE,0,mpi_comm_internode);
	}
	}
        }
        MPI_Barrier(mpi_comm_shared);

	}
      */

      // code gen -- get rid of this eventually...

      MPI_Gather_set(edgeMatchSet,0,mpi_comm);
      if (mpi_rank == 0) {
        if (!edgeMatchSet.empty()) {
          cout << "WARNING: the following cases need to be added to aligned edge matching:" << endl;
          for (set<int>::iterator it = edgeMatchSet.begin(); it != edgeMatchSet.end(); ++it) {
            cout << "CODE: case " << *it << ":" << endl;
          }
        }
      }

      MPI_Gather_set(unhandledMatchSet,0,mpi_comm);
      if (mpi_rank == 0) {
        if (!unhandledMatchSet.empty()) {
          cout << "WARNING: the following cases were NOT handled:" << endl;
          for (set<int>::iterator it = unhandledMatchSet.begin(); it != unhandledMatchSet.end(); ++it) {
            cout << "ERROR: case " << *it << ":" << endl;
          }
        }
      }

      // coplanar tris...

      MPI_Gather_set(coplanarTriSet,0,mpi_comm);
      MPI_Gather_set(coplanarZoneSet,0,mpi_comm);
      if (mpi_rank == 0) {
        if (!coplanarTriSet.empty()) {
          cout << " > WARNING: got " << coplanarTriSet.size() << " coplanar tris. See file \"debug_coplanar.dat\"" << endl;
          cout << " > zone pairs for coplanar tris:" << endl;
          for (set<pair<int,int> >::iterator it = coplanarZoneSet.begin(); it != coplanarZoneSet.end(); ++it) {
            cout << "   > \"" << zoneVec[it->first].getName() << "\" and \"" << zoneVec[it->second].getName() << "\"" << endl;
          }
          if (iter == (niter-1))
            writeTriDebugFile("debug_coplanar.dat",coplanarTriSet);
        }
        else {
          cout << " > no coplanar tris" << endl;
        }
      }

      if (b_repair) {
        for (int ipart = 0; ipart < partVec.size(); ++ipart) {
          if (partVec[ipart]->surface) {
            //if (partVec[ipart]->surface->dist_sp == NULL)
            //  CTI_Mmap_rw(partVec[ipart]->surface->dist_sp,partVec[ipart]->surface->nsp);
            // we also need flag_sp to record intersections...
            if (partVec[ipart]->surface->flag_sp == NULL)
              CTI_Mmap_rw(partVec[ipart]->surface->flag_sp,partVec[ipart]->surface->nsp);
            int nsp_split = partVec[ipart]->surface->nsp/mpi_size_shared;
            if (partVec[ipart]->surface->nsp%mpi_size_shared != 0) nsp_split += 1;
            const int isp0 = min(partVec[ipart]->surface->nsp,nsp_split*mpi_rank_shared);
            const int isp1 = min(partVec[ipart]->surface->nsp,nsp_split*(mpi_rank_shared+1));
            for (int isp = isp0; isp < isp1; ++isp) {
              //partVec[ipart]->surface->dist_sp[isp] = 0.0;
              partVec[ipart]->surface->flag_sp[isp] = 0;
            }
          }
        }
        MPI_Barrier(mpi_comm_shared);
      }

      int my_count = doubleVec.size();
      int count;
      MPI_Allreduce(&my_count,&count,1,MPI_INT,MPI_SUM,mpi_comm);
      if (count > 0) {
        // ------------------------------------------------------------
        // everyone reduce their pairs of intersections in the zonePairSet
        // to rank 0 for condensing and reporting...
        // ------------------------------------------------------------
        MPI_Gather_set(intersectingTriSet,0,mpi_comm);
        MPI_Gather_set(intersectingZoneSet,0,mpi_comm);
        if (mpi_rank == 0) {
          assert(!intersectingTriSet.empty());
          assert(!intersectingZoneSet.empty());
          cout << " > WARNING: got " << intersectingTriSet.size() << " intersecting tris. See file \"debug_intersecting.dat\"" << endl;
          cout << " > zone pairs for intersecting tris:" << endl;
          for (set<pair<int,int> >::iterator it = intersectingZoneSet.begin(); it != intersectingZoneSet.end(); ++it) {
            cout << "   > \"" << zoneVec[it->first].getName() << "\" and \"" << zoneVec[it->second].getName() << "\"" << endl;
          }
          if (iter == (niter-1)) {
            char filename[128];
            sprintf(filename,"debug_intersecting.%04d.dat",iter);
            writeTriDebugFile(filename,intersectingTriSet);
          }
        }

        if (b_repair) {
          // use the reduced set on rank 0 to flag spost for intersecting tris...
          if (mpi_rank == 0) {
            for (set<pair<int,int> >::iterator it = intersectingTriSet.begin(); it != intersectingTriSet.end(); ++it) {
              const int ipart = it->first;
              const int ist = it->second;
              assert(partVec[ipart]->surface);
              assert(partVec[ipart]->surface->flag_sp);
              // flag the sp's...
              FOR_I3 {
                const int isp = partVec[ipart]->surface->spost[ist][i];
                partVec[ipart]->surface->flag_sp[isp] |= 1;
              }
            }
          }
          // and Bcast...
          if (mpi_rank_shared == 0) {
            for (int ipart = 0; ipart < partVec.size(); ++ipart) {
              if (partVec[ipart]->surface) {
                assert(partVec[ipart]->surface->flag_sp);
                MPI_Bcast(partVec[ipart]->surface->flag_sp,partVec[ipart]->surface->nsp,MPI_INT,0,mpi_comm_internode);
              }
            }
          }
          MPI_Barrier(mpi_comm_shared);
        }

        // ------------------------------------------------------------
        // write the distributed points to a single file using handshaking...
        // ------------------------------------------------------------
        if (iter == (niter-1)) {
          FILE * fp = NULL;
          if (mpi_rank == 0) {
            cout << " > writing intersection points to \"debug_xp.dat\"" << endl;
            fp = fopen("debug_xp.dat","w");
            assert(fp != NULL);
          }
          else {
            int dummy;
            MPI_Status status;
            MPI_Recv(&dummy,1,MPI_INT,mpi_rank-1,4321,mpi_comm,&status);
            if (!doubleVec.empty()) {
              fp = fopen("debug_xp.dat","a");
              assert(fp != NULL);
            }
          }
          if (fp != NULL) {
            for (int ii = 0; ii < doubleVec.size(); ii += 3) {
              fprintf(fp,"%18.15le %18.15le %18.15le\n",doubleVec[ii],doubleVec[ii+1],doubleVec[ii+2]);
            }
            fclose(fp);
          }
          if (mpi_rank < mpi_size-1) {
            int dummy = 1;
            MPI_Send(&dummy,1,MPI_INT,mpi_rank+1,4321,mpi_comm);
          }
          MPI_Barrier(mpi_comm);
        }
      }
      else {
        if (mpi_rank == 0) cout << " > no intersecting tris" << endl;
      }

      if (mpi_rank == 0) cout << " > run time for diagnostics: " << MPI_Wtime()-wtime0 << " [s]" << endl;

      // ===========================================================
      // actually do it -- put this somewhere else eventually?
      // ===========================================================

      if (b_repair) {

        /*
	  if (mpi_rank == 0) {
          int my_count[4] = { 0, 0, 0, 0 };
          double my_buf[4] = { 0.0, 0.0, 0.0, 0.0 };
          for (int ipart = 0; ipart < partVec.size(); ++ipart) {
	  if (partVec[ipart]->surface) {
	  assert(partVec[ipart]->surface->dist_sp);
	  assert(partVec[ipart]->surface->flag_sp);
	  for (int isp = 0; isp < partVec[ipart]->surface->nsp; ++isp) {
	  if (partVec[ipart]->surface->flag_sp[isp]&1) {
	  if (partVec[ipart]->surface->dist_sp[isp] <= 0.0) {
	  ++my_count[0];
	  // record the min negative...
	  my_buf[0] = min(my_buf[0],partVec[ipart]->surface->dist_sp[isp]);
	  }
	  else {
	  ++my_count[1];
	  // also record how positive...
	  my_buf[1] = max(my_buf[1],partVec[ipart]->surface->dist_sp[isp]);
	  }
	  }
	  else {
	  if (partVec[ipart]->surface->dist_sp[isp] <= 0.0) {
	  // this is a negative value on a non-intersecting tri. It is either a
	  // thin fluid region or an overlap away from intersecting tris.
	  ++my_count[2];
	  my_buf[2] = min(my_buf[2],partVec[ipart]->surface->dist_sp[isp]);
	  }
	  else {
	  // positive away from an intersection: should be the norm...
	  ++my_count[3];
	  my_buf[3] = max(my_buf[3],partVec[ipart]->surface->dist_sp[isp]);
	  }
	  }
	  }
	  }
          }
          cout << "OVERLAP_SEPARATION stats" << endl;
          cout << " > user provided length scale: " << overlap_separation << endl;
          cout << " > negative distances on intersecting tris    : " << my_count[0] << " min dist: " << my_buf[0] << endl;
          cout << " > positive distances on intersecting tris    : " << my_count[1] << " max dist: " << my_buf[1] << endl;
          cout << " > negative distances on non-intersecting tris: " << my_count[2] << " min dist: " << my_buf[2] << endl;
          cout << " > positive distances on non-intersecting tris: " << my_count[3] << " max dist: " << my_buf[3] << endl;
	  }
	  // write certain negative values...
	  if (mpi_rank == 0) {
	  FILE * fp = fopen("dist_sp_very_neg.dat","w");
	  assert(fp);
	  for (int isp = 0; isp < partVec[0]->surface->nsp; ++isp) {
	  if ((partVec[0]->surface->dist_sp[isp] < -0.5*overlap_separation)&&(partVec[0]->surface->flag_sp[isp]&1)) {
	  fprintf(fp,"%18.15le %18.15le %18.15le %18.15le\n",
	  partVec[0]->surface->xsp[isp][0],
	  partVec[0]->surface->xsp[isp][1],
	  partVec[0]->surface->xsp[isp][2],
	  partVec[0]->surface->dist_sp[isp]);
	  }
	  }
	  fclose(fp);
	  }
        */

        int nsp_int = 0;
        for (int ipart = 0; ipart < partVec.size(); ++ipart) {
          if (partVec[ipart]->surface) {
            assert(partVec[ipart]->surface->flag_sp);
            for (int isp = 0; isp < partVec[ipart]->surface->nsp; ++isp) {
              if (partVec[ipart]->surface->flag_sp[isp]&1)
                ++nsp_int;
            }
          }
        }
        if ((iter == 0)&&(nsp_int > 0))
          b_requires_repair = true;

        if (mpi_rank == 0)
          cout << " > nsp_int: " << nsp_int << endl;
        if (nsp_int == 0) {
          if (mpi_rank == 0)
            cout << " > no more negative intersection points. leaving..." << endl;
          break;
        }

        // count the number of local nodes associated with intersections...
        int my_nsp_int = 0;
        for (int ipart = 0; ipart < partVec.size(); ++ipart) {
          if (partVec[ipart]->surface) {
            int nsp_split = partVec[ipart]->surface->nsp/mpi_size_shared;
            if (partVec[ipart]->surface->nsp%mpi_size_shared != 0) nsp_split += 1;
            const int isp0 = min(partVec[ipart]->surface->nsp,nsp_split*mpi_rank_shared);
            const int isp1 = min(partVec[ipart]->surface->nsp,nsp_split*(mpi_rank_shared+1));
            assert(partVec[ipart]->surface->flag_sp);
            for (int isp = isp0; isp < isp1; ++isp) {
              if (partVec[ipart]->surface->flag_sp[isp]&1)
                ++my_nsp_int;
            }
          }
        }

        // compute a local node perturbation for our intersection nodes...
        double (*my_dx_sp)[3] = new double[my_nsp_int][3];
        my_nsp_int = 0;
        int my_count[4] = { 0, 0, 0, 0 };
        for (int ipart = 0; ipart < partVec.size(); ++ipart) {
          if (partVec[ipart]->surface) {
            partVec[ipart]->surface->ensureStosp();
            int nsp_split = partVec[ipart]->surface->nsp/mpi_size_shared;
            if (partVec[ipart]->surface->nsp%mpi_size_shared != 0) nsp_split += 1;
            const int isp0 = min(partVec[ipart]->surface->nsp,nsp_split*mpi_rank_shared);
            const int isp1 = min(partVec[ipart]->surface->nsp,nsp_split*(mpi_rank_shared+1));
            for (int isp = isp0; isp < isp1; ++isp) {
              if (partVec[ipart]->surface->flag_sp[isp]&1) {
                const int my_isp_int = my_nsp_int++;
                FOR_I3 my_dx_sp[my_isp_int][i] = 0.0;
                double wgt = 0.0;
                for (int sos = partVec[ipart]->surface->stosp_i[isp]; sos != partVec[ipart]->surface->stosp_i[isp+1]; ++sos) {
                  const int ist = partVec[ipart]->surface->stosp_v[sos];
                  FOR_I3 {
                    const int isp_nbr = partVec[ipart]->surface->spost[ist][i];
                    const double dx[3] = DIFF(partVec[ipart]->surface->xsp[isp_nbr],partVec[ipart]->surface->xsp[isp]);
                    FOR_J3 my_dx_sp[my_isp_int][j] += dx[j];
                    wgt += 1.0;
                  }
                }
                FOR_I3 my_dx_sp[my_isp_int][i] /= wgt;
                // this is a good direction and distance to move the point, but we only want to apply this motion
                // when it is in the OPPOSITE direction of the surface normal. i.e. it tends to move the
                // intersecting surfaces apart...
                const double mag = MAG(my_dx_sp[my_isp_int]);
                if (DOT_PRODUCT(partVec[ipart]->surface->n_sp[isp],my_dx_sp[my_isp_int]) < 0.0) {
                  if (mag > dn) {
                    FOR_I3 my_dx_sp[my_isp_int][i] *= dn/mag;
                    ++my_count[0];
                  }
                  else {
                    ++my_count[1];
                  }
                }
                else {
                  // if "smoothing" moves the normal towards the other surface (i.e. potentially further into
                  // trouble), then leave the node where it is...
                  if (mag > neg_frac*dn) {
                    FOR_I3 my_dx_sp[my_isp_int][i] *= -neg_frac*dn/mag;
                    ++my_count[2];
                  }
                  else {
                    FOR_I3 my_dx_sp[my_isp_int][i] *= -neg_frac;
                    ++my_count[3];
                  }
                }
              }
            }
          }
        }
        int count[4];
        MPI_Reduce(my_count,count,4,MPI_INT,MPI_SUM,0,mpi_comm_shared);
        if (mpi_rank == 0)
          cout << " > perturbations in smoothing direction: limited to dn: " << count[0] <<
            " not limited: " << count[1] <<
            ", non-smoothing direction: limited to -frac*dn: " << count[2] <<
            " limited to -frac: " << count[3] <<
            " total: " << count[0]+count[1]+count[2]+count[3] << endl;

        // and apply...
        my_nsp_int = 0;
        for (int ipart = 0; ipart < partVec.size(); ++ipart) {
          if (partVec[ipart]->surface) {
            int nsp_split = partVec[ipart]->surface->nsp/mpi_size_shared;
            if (partVec[ipart]->surface->nsp%mpi_size_shared != 0) nsp_split += 1;
            const int isp0 = min(partVec[ipart]->surface->nsp,nsp_split*mpi_rank_shared);
            const int isp1 = min(partVec[ipart]->surface->nsp,nsp_split*(mpi_rank_shared+1));
            for (int isp = isp0; isp < isp1; ++isp) {
              if (partVec[ipart]->surface->flag_sp[isp]&1) {
                const int my_isp_int = my_nsp_int++;
                FOR_I3 partVec[ipart]->surface->xsp[isp][i] += my_dx_sp[my_isp_int][i];
              }
            }
          }
        }
        delete[] my_dx_sp;
        MPI_Barrier(mpi_comm_shared);

      }

    }
    DELETE(ist_flat);
    //DELETE(spost_local);
    //DELETE(psposp_local);
    //DELETE(dist_sp_local);

    if (b_repair&&b_requires_repair) {

      char filename[128];
      for (int ipart = 0; ipart < partVec.size(); ++ipart) {
        if (partVec[ipart]->surface) {
          sprintf(filename,"repaired.%02d.sbin",ipart);
          partVec[ipart]->surface->writeBinary(filename);
        }
      }

      // recall the initial xsp on the surface was stored in surface->dx_sp. use the final
      // xsp to compute the local dx, and dot on the local surface unit normal to get the
      // displacement in the normal direction...

      double my_buf[3] = { 0.0, -HUGE_VAL, -HUGE_VAL };
      for (int ipart = 0; ipart < partVec.size(); ++ipart) {
        if (partVec[ipart]->surface) {
          assert(partVec[ipart]->surface->dx_sp);
          assert(partVec[ipart]->surface->n_sp); // normal should be current if we converged (exited without changing node positions on last iter)
          if (partVec[ipart]->surface->dist_sp == NULL)
            CTI_Mmap_rw(partVec[ipart]->surface->dist_sp,partVec[ipart]->surface->nsp);
          int nsp_split = partVec[ipart]->surface->nsp/mpi_size_shared;
          if (partVec[ipart]->surface->nsp%mpi_size_shared != 0) nsp_split += 1;
          const int isp0 = min(partVec[ipart]->surface->nsp,nsp_split*mpi_rank_shared);
          const int isp1 = min(partVec[ipart]->surface->nsp,nsp_split*(mpi_rank_shared+1));
          for (int isp = isp0; isp < isp1; ++isp) {
            FOR_I3 partVec[ipart]->surface->dx_sp[isp][i] = partVec[ipart]->surface->xsp[isp][i] - partVec[ipart]->surface->dx_sp[isp][i];
            partVec[ipart]->surface->dist_sp[isp] = DOT_PRODUCT(partVec[ipart]->surface->dx_sp[isp],partVec[ipart]->surface->n_sp[isp]);
            const double mag = MAG(partVec[ipart]->surface->dx_sp[isp]);
            my_buf[0] = max(my_buf[0],mag);
            my_buf[1] = max(my_buf[1],-partVec[ipart]->surface->dist_sp[isp]);
            my_buf[2] = max(my_buf[2],partVec[ipart]->surface->dist_sp[isp]);
          }
        }
      }
      double buf[3];
      MPI_Reduce(my_buf,buf,3,MPI_DOUBLE,MPI_MAX,0,mpi_comm_shared);
      if (mpi_rank == 0) {
        cout << " > max surface point displacement: " << buf[0] << ", surface normal displacement range: " << -buf[1] << " to " << buf[2] << endl;
      }

      if (b_write_displacement) {

        // now write the tecplot file(s) of displacement...
        for (int ipart = 0; ipart < partVec.size(); ++ipart) {
          if (partVec[ipart]->surface) {
            sprintf(filename,"displacement.%02d.dat",ipart);
            partVec[ipart]->surface->writeTecplot(filename,partVec[ipart]->surface->dist_sp);
          }
        }

      }

    }
  }

  void prepareIsInside() {

    if (mpi_rank == 0) cout << "prepareIsInside()" << endl;

    for (int ipart = 0; ipart < partVec.size(); ++ipart) {
      // --------------------------
      // farfield...
      // --------------------------
      if (partVec[ipart]->hasFF() && (!partVec[ipart]->hasSimpleFF())) {
        const double * const e1_ff = partVec[ipart]->e1_ff;
        const double * const e2_ff = partVec[ipart]->e2_ff;
        // this part may need its isInsideFF adt's rebuilt if they are needed...
        bool b_ff = false;
        int xy_ff_min[2] = { TWO_BILLION, TWO_BILLION };
        int xy_ff_max[2] = { -TWO_BILLION, -TWO_BILLION };
        // every part with points below us will check against our ff...
        for (int ipart1 = 0; ipart1 < ipart; ++ipart1) {
          if (partVec[ipart1]->pts) {
            b_ff = true;
            if (!b_dtoi) setDtoiStuff();
            for (int ip = 0; ip < partVec[ipart1]->pts->np; ++ip) {
              int xy[2];
              xy[0] = getEvenInt((partVec[ipart1]->pts->xp[ip][0]-dtoi_x0[0])*e1_ff[0] +
                                 (partVec[ipart1]->pts->xp[ip][1]-dtoi_x0[1])*e1_ff[1] +
                                 (partVec[ipart1]->pts->xp[ip][2]-dtoi_x0[2])*e1_ff[2],dtoi_delta);
              xy[1] = getEvenInt((partVec[ipart1]->pts->xp[ip][0]-dtoi_x0[0])*e2_ff[0] +
                                 (partVec[ipart1]->pts->xp[ip][1]-dtoi_x0[1])*e2_ff[1] +
                                 (partVec[ipart1]->pts->xp[ip][2]-dtoi_x0[2])*e2_ff[2],dtoi_delta);
              FOR_I2 {
                xy_ff_min[i] = min(xy_ff_min[i],xy[i]);
                xy_ff_max[i] = max(xy_ff_max[i],xy[i]);
              }
            }
          }
        }
        if (b_ff) partVec[ipart]->prepareIsInsideFF(xy_ff_min,xy_ff_max);
      }
      // --------------------------
      // solid...
      // --------------------------
      //
      if (!partVec[ipart]->hasSimpleSolid()) {
        const double * const e1_solid = partVec[ipart]->e1_solid;
        const double * const e2_solid = partVec[ipart]->e2_solid;
        // this part may need its isInsideSolid adt's rebuilt if they are needed...
        bool b_solid = false;
        int xy_solid_min[2] = { TWO_BILLION, TWO_BILLION };
        int xy_solid_max[2] = { -TWO_BILLION, -TWO_BILLION };
        // every part with points above us will check against our solid...
        for (int ipart1 = ipart+1; ipart1 < partVec.size(); ++ipart1) {
          if (partVec[ipart1]->pts) {
            b_solid = true;
            if (!b_dtoi) setDtoiStuff();

            // cout << "part: " << ipart1 << endl;
            // cout << "dtoi_x0: " << COUT_VEC(dtoi_x0) << ", delta: " << dtoi_delta << endl;
            // cout << "e1_solid: " << COUT_VEC(e1_solid) << ", e2_solid: " << COUT_VEC(e2_solid) << endl;

            for (int ip = 0; ip < partVec[ipart1]->pts->np; ++ip) {
              int xy[2];
              xy[0] = getEvenInt((partVec[ipart1]->pts->xp[ip][0]-dtoi_x0[0])*e1_solid[0] +
                                 (partVec[ipart1]->pts->xp[ip][1]-dtoi_x0[1])*e1_solid[1] +
                                 (partVec[ipart1]->pts->xp[ip][2]-dtoi_x0[2])*e1_solid[2],dtoi_delta);
              xy[1] = getEvenInt((partVec[ipart1]->pts->xp[ip][0]-dtoi_x0[0])*e2_solid[0] +
                                 (partVec[ipart1]->pts->xp[ip][1]-dtoi_x0[1])*e2_solid[1] +
                                 (partVec[ipart1]->pts->xp[ip][2]-dtoi_x0[2])*e2_solid[2],dtoi_delta);
              FOR_I2 {
                xy_solid_min[i] = min(xy_solid_min[i],xy[i]);
                xy_solid_max[i] = max(xy_solid_max[i],xy[i]);
              }
            }
          }
        }
        if (!partVec[ipart]->hasFF()) {
          // for parts without a ff, earlier parts will check against our isInsideSolid...
          for (int ipart1 = 0; ipart1 < ipart; ++ipart1) {
            if (partVec[ipart1]->pts) {
              b_solid = true;
              if (!b_dtoi) setDtoiStuff();
              for (int ip = 0; ip < partVec[ipart1]->pts->np; ++ip) {
                int xy[2];
                xy[0] = getEvenInt((partVec[ipart1]->pts->xp[ip][0]-dtoi_x0[0])*e1_solid[0] +
                                   (partVec[ipart1]->pts->xp[ip][1]-dtoi_x0[1])*e1_solid[1] +
                                   (partVec[ipart1]->pts->xp[ip][2]-dtoi_x0[2])*e1_solid[2],dtoi_delta);
                xy[1] = getEvenInt((partVec[ipart1]->pts->xp[ip][0]-dtoi_x0[0])*e2_solid[0] +
                                   (partVec[ipart1]->pts->xp[ip][1]-dtoi_x0[1])*e2_solid[1] +
                                   (partVec[ipart1]->pts->xp[ip][2]-dtoi_x0[2])*e2_solid[2],dtoi_delta);
                FOR_I2 {
                  xy_solid_min[i] = min(xy_solid_min[i],xy[i]);
                  xy_solid_max[i] = max(xy_solid_max[i],xy[i]);
                }
              }
            }
          }
        }
        if (b_solid) partVec[ipart]->prepareIsInsideSolid(xy_solid_min,xy_solid_max);
      }
    }

    if (mpi_rank == 0) cout << " > done prepareIsInside()" << endl;

  }

  void buildDisjointAndCreaseAngleGroups(const double crease_angle_degrees) {

    // just work on part 0's surface...

    assert(!partVec.empty());
    SurfaceShm * surface = partVec[0]->surface;
    assert(surface != NULL);

    if (mpi_rank == 0) cout << "buildDisjointAndCreaseAngleGroups: got surface nsp,nst: " << surface->nsp << " " << surface->nst << endl;

    // building a teost-like structure in parallel is a pain: try rendezvous...

    int nst_split = surface->nst/mpi_size;
    if (surface->nst%mpi_size != 0) nst_split += 1;
    const int ist_f = min(surface->nst,nst_split*mpi_rank);
    const int ist_l = min(surface->nst,nst_split*(mpi_rank+1))-1;
    //const int my_nst = ist_l-ist_f+1;

    int nsp_split = surface->nsp/mpi_size;
    if (surface->nsp%mpi_size != 0) nsp_split += 1;
    //const int isp_f = min(surface->nsp,nsp_split*mpi_rank);
    //const int isp_l = min(surface->nsp,nsp_split*(mpi_rank+1))-1;
    //const int my_nsp = isp_l-isp_f+1;

    // go through our nst's and pack and send ist,isp0,isp1 to
    // ranks for both isp0 and isp1...

    int * send_count = new int[mpi_size];
    FOR_RANK send_count[rank] = 0;

    // count...
    for (int ist = ist_f; ist <= ist_l; ++ist) {
      FOR_I3 {
        const int isp0 = surface->spost[ist][i];
        const int rank0 = isp0/nsp_split;
        assert((rank0 >= 0)&&(rank0 < mpi_size));
        send_count[rank0] += 3;
        const int isp1 = surface->spost[ist][(i+1)%3];
        const int rank1 = isp1/nsp_split;
        assert((rank1 >= 0)&&(rank1 < mpi_size));
        if (rank1 != rank0)
          send_count[rank1] += 3;
      }
    }

    int * send_disp = new int[mpi_size];
    send_disp[0] = 0;
    for (int rank = 1; rank < mpi_size; ++rank)
      send_disp[rank] = send_count[rank-1] + send_disp[rank-1];
    const int send_count_sum = send_disp[mpi_size-1] + send_count[mpi_size-1];

    int * send_buf_int = new int[send_count_sum];

    // pack...
    for (int ist = ist_f; ist <= ist_l; ++ist) {
      FOR_I3 {
        const int isp0 = surface->spost[ist][i];
        const int isp1 = surface->spost[ist][(i+1)%3];
        const int rank0 = isp0/nsp_split;
        send_buf_int[send_disp[rank0]] = ist;
        send_buf_int[send_disp[rank0]+1] = isp0;
        send_buf_int[send_disp[rank0]+2] = isp1;
        send_disp[rank0] += 3;
        const int rank1 = isp1/nsp_split;
        if (rank1 != rank0) {
          send_buf_int[send_disp[rank1]] = ist;
          send_buf_int[send_disp[rank1]+1] = isp0;
          send_buf_int[send_disp[rank1]+2] = isp1;
          send_disp[rank1] += 3;
        }
      }
    }

    // rewind disp...
    send_disp[0] = 0;
    FOR_RANK send_disp[rank+1] = send_disp[rank] + send_count[rank];

    // send->recv...
    int * recv_count = new int[mpi_size];
    MPI_Alltoall(send_count,1,MPI_INT,recv_count,1,MPI_INT,mpi_comm);

    int * recv_disp = new int[mpi_size];
    recv_disp[0] = 0;
    for (int rank = 1; rank < mpi_size; ++rank)
      recv_disp[rank] = recv_count[rank-1] + recv_disp[rank-1];
    const int recv_count_sum = recv_disp[mpi_size-1] + recv_count[mpi_size-1];

    int * recv_buf_int = new int[recv_count_sum];
    MPI_Alltoallv(send_buf_int,send_count,send_disp,MPI_INT,
		  recv_buf_int,recv_count,recv_disp,MPI_INT,mpi_comm);

    // now on the recv side, we have all the tris and their edges...

    map<const pair<int,int>,pair<int,int> > edgeMap;
    for (int irecv = 0; irecv < recv_count_sum; irecv += 3) {
      const int ist = recv_buf_int[irecv];
      const int isp0 = recv_buf_int[irecv+1];
      const int isp1 = recv_buf_int[irecv+2];
      map<const pair<int,int>,pair<int,int> >::iterator iter = edgeMap.find(pair<int,int>(isp1,isp0));
      if (iter == edgeMap.end()) {
        edgeMap.insert(pair<pair<int,int>,pair<int,int> >(pair<int,int>(isp0,isp1),pair<int,int>(ist,irecv/3)));
      }
      else {
        const int ist_nbr = iter->second.first;
        const int index_nbr = iter->second.second;
        recv_buf_int[index_nbr] = ist;
        recv_buf_int[irecv/3] = ist_nbr;
        edgeMap.erase(iter);
      }
    }
    if (!edgeMap.empty()) {
      // these must be open edges...
      cout << "edgeMap.size(): " << edgeMap.size() << endl;
      assert(0);
    }

    // send back...
    FOR_RANK {
      send_count[rank] /= 3;
      send_disp[rank] /= 3;
      recv_count[rank] /= 3;
      recv_disp[rank] /= 3;
    }

    MPI_Alltoallv(recv_buf_int,recv_count,recv_disp,MPI_INT,
                  send_buf_int,send_count,send_disp,MPI_INT,
		  mpi_comm);
    delete[] recv_buf_int;
    delete[] recv_count;
    delete[] recv_disp;

    int (*stost)[3] = new int[nst_split][3];
    for (int ist = ist_f; ist <= ist_l; ++ist) {
      FOR_I3 {
        const int isp0 = surface->spost[ist][i];
        const int rank0 = isp0/nsp_split;
        stost[ist-ist_f][i] = send_buf_int[send_disp[rank0]];
        ++send_disp[rank0];
        const int isp1 = surface->spost[ist][(i+1)%3];
        const int rank1 = isp1/nsp_split;
        if (rank1 != rank0) {
          assert(stost[ist-ist_f][i] == send_buf_int[send_disp[rank1]]);
          ++send_disp[rank1];
        }
      }
    }
    delete[] send_buf_int;
    delete[] send_count;
    delete[] send_disp;

    // ===================================================
    // we now have parallel stost...
    // ===================================================

    // start with each tri in its own group...

    int * group_st = new int[nst_split];
    for (int ist = ist_f; ist <= ist_l; ++ist)
      group_st[ist-ist_f] = ist;

    for (int ist = ist_f; ist <= ist_l; ++ist) {
      FOR_I3 {
        const int ist_nbr = stost[ist-ist_f][i];
        assert((ist_nbr >= 0)&&(ist_nbr < surface->nst));
        if ((surface->znost[ist] == surface->znost[ist_nbr])&&
            checkCoplanar(0,min(ist,ist_nbr),0,0,max(ist,ist_nbr),0,150)) {
          // these tris belong together...

        }
        else {
          // these tris do not belong in the same patch (at least not
          // across this edge, so set stost to -ve 1 indexing...
          stost[ist-ist_f][i] = -stost[ist-ist_f][i]-1;
        }
      }
    }







    MPI_Pause("XXXXXXXXXX OKOK XXXXXXXXXXXXXX");

  }

  void queryLengthScales() {

    buildDisjointAndCreaseAngleGroups(150);

  }

  void queryParts() {
    WUI(INFO,"There are currently " << partVec.size() << " PARTS:");
    for (int ipart = 0; ipart < partVec.size(); ++ipart) {
      WUI(INFO,"  " << ipart << " " << partVec[ipart]->getName());
      if (partVec[ipart]->surface == NULL) {
        WUI(INFO,"    > part has no surface");
      }
      else {
        SurfaceShm * surface = partVec[ipart]->surface;
        WUI(INFO,"    > part has surface with " << surface->zoneVec.size() << " ZONES:");
        for (int isz = 0; isz < surface->zoneVec.size(); ++isz) {
          WUI(INFO,"     " << isz << " " << surface->zoneVec[isz].getName() << " periodic: " << surface->zoneVec[isz].isPeriodic());
        }
      }
    }
  }

  void queryZones() {
    WUI(INFO,"There are currently " << zoneVec.size() << " ZONES:");
    for (int izone = 0; izone < zoneVec.size(); ++izone) {
      WUI(INFO,"  " << izone << " " << zoneVec[izone].getName() << " periodic: " << zoneVec[izone].isPeriodic());
    }
  }

  void queryZones(const string& names_cds) {
    // assume comma-delimited list of zone names...
    int ierr = 0;
    set<int> zoneIdSet;
    vector<string> zone_names;
    MiscUtils::splitCsv(zone_names,names_cds);
    for (int ii = 0; ii < zone_names.size(); ++ii) {
      const int izone = getZoneIndex(zone_names[ii]); // get the unique zoneIndex in the partVec
      if (izone >= 0) {
        zoneIdSet.insert(izone);
      }
      else {
        WUI(WARN,"Unrecognized ZONE: " << zone_names[ii]);
        ierr |= 1;
      }
    }
    if (ierr == 0) {
      assert(!zoneIdSet.empty());
      queryZones(zoneIdSet);
    }
    else {
      WUI(WARN,"Fix Error and try again");
    }
  }

  // formerly PartMaker...

  void processPart(list<Param>::iterator param,const bool b_help) {

    processPart(&(*param),b_help);

  }
  
  void processPart(Param * param,const bool b_help) {

    if (b_help) {
      helpProcessPart();
    }
    else {

      // allocate a new part and make it from the param...
      Part * part = new Part();
      if (makePart(part,param)) {
        // part was sucessfully built...
        // the part's ipart_of_part should contain the index in the partVec...
        if (part->ipart_of_part != partVec.size()) {
          // this is a part replacement. clear the old part...
          delete partVec[part->ipart_of_part];
          partVec[part->ipart_of_part] = part;
        }
        else {
          partVec.push_back(part);
        }
        // and the part's surface_st_bits...
        if (part->surface) {
          part->init_surface_st_bits();
        }
        // after any part is added to the part list, it is necessary to
        // rebuild a few things...
        buildOrRebuildZoneVec();
        b_bbminmax = false; // this will recompute the overall bbox the next time it is req'd
      }
      else {
        delete part;
        WUI(WARN,"Error: PART failed. check syntax");
      }

    }

  }

  void helpProcessPart() {
    WUI(INFO,
        "PART is the driver for including surface and, optionally, points. Some examples:\n" <<
        "PART NAME=main SURF SBIN ../../../../runs/gb_jet_nextgen/gp_jet.sbin PTS HCP\n" <<
        "PART SURF SBIN ../../../../runs/gb_jet_nextgen/gp_jet.sbin PTS HCP\n" <<
        "PART NAME=bl STRAND ZONE_NAMES=main.nozzle-internal-pipe DN 0.002 DT 0.016 N 20 NLAYERS 3,6,12\n" <<
        "PART bl.part");
  }

  bool makePart(Part * part,Param * param) {

    //bool b_hcp = false;
    part->ipart_of_part = partVec.size(); // for routines that require the eventual part index

    int ierr = 0;
    int iarg = 0;
    while (iarg < param->size()) {
      // look for a recognized filename extension first...
      string token = param->getString(iarg++);
      if ((token.size() > 5)&&(token.compare(token.size()-5,5,".sbin") == 0)) {
        // *.sbin file...
        const int this_ierr = makePartFromSbin(part,token);
        if (this_ierr != 0) {
          WUI(WARN,"problem creating PART from sbin file: " << token);
          ierr = -1;
        }
      }
      else if ((token.size() > 5)&&(token.compare(token.size()-5,5,".part") == 0)) {
        // *.part file...
	const int this_ierr = part->readBinary(token);
        if (this_ierr != 0) {
          WUI(WARN,"problem creating PART from part file: " << token);
	  return false;
          ierr = -1;
        }
      }
      else if ((token.size() > 5)&&(token.compare(token.size()-5,5,".mles") == 0)) {
        // *.mles file...
        const int this_ierr = makePartFromMles(part,token);
        if (this_ierr != 0) {
          WUI(WARN,"problem creating PART from mles file: " << token);
          ierr = -1;
        }
      }
      else {
        token = MiscUtils::toUpperCase(token);
        if (token == "NAME") {
          part->b_name = true;
          part->name = param->getString(iarg++);
          // check if this is an existing part...
          for (int ipart = 0; ipart < partVec.size(); ++ipart) {
            if ((partVec[ipart]->b_name)&&(part->name == partVec[ipart]->name)) {
              // replace this part in the partVec...
              part->ipart_of_part = ipart;
              if (mpi_rank == 0) cout << "XXXXX this part is replacing an existing part: " << ipart << endl;
              break;
            }
          }
        }
        else if (token == "SURF") {
          // look for a recognized filename extension first...
          token = param->getString(iarg++);
          if ((token.size() > 5)&&(token.compare(token.size()-5,5,".sbin") == 0)) {
            const int this_ierr = makePartFromSbin(part,token);
            if (this_ierr != 0) {
              WUI(WARN,"problem creating PART SURF for sbin file: " << token);
              ierr = -1;
            }
          }
          else {
            token = MiscUtils::toUpperCase(token);
            if ((token == "SBIN")||(token == "BIN")) {
              token = param->getString(iarg++);
              const int this_ierr = makePartFromSbin(part,token);
              if (this_ierr != 0) {
                WUI(WARN,"problem creating PART SURF for sbin file: " << token);
                ierr = -1;
              }
            }
            else if ((token == "BOX")||(token == "SIMPLE_BOX")) {
              // BOX X0 <x0> <y0> <z0> X1 <x1> <y1> <z1>
              // TODO: parsing robustness...
              double x0[3],x1[3];
              token = param->getString(iarg);
              if (token == "X0") {
                ++iarg;
                FOR_I3 x0[i] = param->getDouble(iarg++);
                token = param->getString(iarg++);
                assert(token == "X1");
                FOR_I3 x1[i] = param->getDouble(iarg++);
              }
              else {
                x0[0] = param->getDouble(iarg++);
                x1[0] = param->getDouble(iarg++);
                x0[1] = param->getDouble(iarg++);
                x1[1] = param->getDouble(iarg++);
                x0[2] = param->getDouble(iarg++);
                x1[2] = param->getDouble(iarg++);
              }
              // also allow the addition of any of PERIODIC_X PERIODIC_Y PERIODIC_Z...
              bool periodic[3] = {false,false,false};
              while (iarg < param->size()) {
                string token = param->getString(iarg);
                if (token == "PERIODIC_X") {
                  ++iarg;
                  periodic[0] = true;
                }
                else if (token == "PERIODIC_Y") {
                  ++iarg;
                  periodic[1] = true;
                }
                else if (token == "PERIODIC_Z") {
                  ++iarg;
                  periodic[2] = true;
                }
                else {
                  break;
                }
              }
              assert(part->surface == NULL);
              part->surface = new SurfaceShm();
              part->surface->makeBox(x0,x1,periodic);
              assert(part->ff_type == UNKNOWN_FF);
              part->ff_type = NO_FF;
              assert(part->solid_type == UNKNOWN_SOLID);

              // specification for box;
              part->solid_type = BOX_SOLID;
              FOR_I3 part->ddata_solid_[i]   = x0[i];
              FOR_I3 part->ddata_solid_[i+3] = x1[i];
            }
            else if (token == "SPHERE") {
              double x[3],r;
              token = param->getString(iarg);
              if ((token == "X0")||(token == "XC")||(token == "X")) {
                // SPHERE X=x y z R=r
                ++iarg;
                FOR_I3 x[i] = param->getDouble(iarg++);
                token = param->getString(iarg++);
                assert(token == "R");
                r = param->getDouble(iarg++);
              }
              else {
                // SPHERE x y z r
                x[0] = param->getDouble(iarg++);
                x[1] = param->getDouble(iarg++);
                x[2] = param->getDouble(iarg++);
                r = param->getDouble(iarg++);
              }
              // if the user has specified an "N" here, this is the number of segments...
              int nseg = 4;
              while (iarg < param->size()) {
                string token = param->getString(iarg);
		if ((token == "N")||(token == "NSEG")) {
                  ++iarg;
                  nseg = param->getInt(iarg++);
                }
                else {
                  break;
                }
              }
              assert(part->surface == NULL);
              part->surface = new SurfaceShm();
	      part->surface->makeSphere("SPHERE",x,r,nseg,true);
              assert(part->ff_type == UNKNOWN_FF);
              part->ff_type = NO_FF;
              assert(part->solid_type == UNKNOWN_SOLID);

              // part->solid_type = CYLINDER_SOLID;
              // FOR_I3 part->ddata_solid[i]   = x0[i];
              // FOR_I3 part->ddata_solid[i+3] = x1[i];
              // part->ddata_solid[6]          = r;

              part->solid_type = SURFACE_SOLID;
            }
            else if ((token == "CYL")||(token == "CYLINDER")||(token == "DISK")) {
              double x0[3],x1[3],r;
              token = param->getString(iarg);
              if (token == "X0") {
                // CYLINDER X0=x y z X1=x y z R=r
                ++iarg;
                FOR_I3 x0[i] = param->getDouble(iarg++);
                token = param->getString(iarg++);
                assert(token == "X1");
                FOR_I3 x1[i] = param->getDouble(iarg++);
                token = param->getString(iarg++);
                assert(token == "R");
                r = param->getDouble(iarg++);
              }
              else {
		// surfer format...
                // CYLINDER x0 y0 z0 x1 y1 z1 r
                x0[0] = param->getDouble(iarg++);
                x0[1] = param->getDouble(iarg++);
                x0[2] = param->getDouble(iarg++);
                x1[0] = param->getDouble(iarg++);
                x1[1] = param->getDouble(iarg++);
                x1[2] = param->getDouble(iarg++);
                r = param->getDouble(iarg++);
              }
              // if the user has specified an "N" here, this is the number of segments...
              int n = 128;
              bool b_periodic = false;
              while (iarg < param->size()) {
                string token = param->getString(iarg);
                if (token == "PERIODIC") {
                  ++iarg;
                  b_periodic = true;
                }
                else if (token == "N") {
                  ++iarg;
                  n = param->getInt(iarg++);
                }
                else {
                  break;
                }
              }
              assert(part->surface == NULL);
              part->surface = new SurfaceShm();
              part->surface->makeCylinder(x0,x1,r,n,true,b_periodic); // true puts on end caps
              assert(part->ff_type == UNKNOWN_FF);
              part->ff_type = NO_FF;
              assert(part->solid_type == UNKNOWN_SOLID);

              // part->solid_type = CYLINDER_SOLID;
              // FOR_I3 part->ddata_solid[i]   = x0[i];
              // FOR_I3 part->ddata_solid[i+3] = x1[i];
              // part->ddata_solid[6]          = r;

              part->solid_type = SURFACE_SOLID;
            }
            else {
              CERR("unrecognized SURF token: " << token);
            }
          }
        }
        else if (token == "PTS") {
          token = MiscUtils::toUpperCase(param->getString(iarg++));
          if (token == "HCP") {
            CWARN("PTS HCP deprecated - this is not necessary");
          }
          else if (token == "CART") {
            CERR("PTS CART not supported. Use PART BOX_WITH_CART_PTS");
          }
          else if (token == "BIN" || token == "PBIN") {
            const string filename = param->getString(iarg++);
            if (part->pts == NULL) {
              part->pts = new Points();
              part->pts->readBinary(filename);
            }
            else {
              part->pts->addFromBinary(filename);
            }
          }
	  else if (token == "MLES") {
            const string filename = param->getString(iarg++);
	    assert(part->pts == NULL);
	    part->pts = new Points();
	    part->pts->readLes(filename); // look for x_vv and r_vv in mles file...
	  }
	  /*
            else if (token == "RESTART") {
            part->pts = new Points();
            const string filename = param->getString(iarg++);
            part->pts->initFromRestart(filename);
            }
          */
          else {
            CERR("PTS " << token << " not supported");
          }
        }
	else if (token == "NCOPY_PTS") {
	  // NCOPY_PTS can be listed after a PTS command and 
	  // will replicate the points based on the specified tokens...
	  const int ncopy = param->getInt(iarg++);
	  assert(ncopy >= 1);
	  token = MiscUtils::toUpperCase(param->getString(iarg++));
	  if (token == "CYL_X") {
	    const double angle_degrees = param->getDouble(iarg++);
	    if (mpi_rank == 0) cout << "copying points CYL_X " << angle_degrees << endl;
	    assert(part->pts);
	    double (*xp_old)[3] = part->pts->xp; assert(xp_old);
	    double *delta_old = part->pts->delta; assert(delta_old);
	    const int np_old = part->pts->np;
	    part->pts->np_global *= ncopy+1;
	    part->pts->np *= ncopy+1;
	    part->pts->xp = new double[part->pts->np][3];
	    part->pts->delta = new double[part->pts->np];
	    for (int ii = 0; ii <= ncopy; ++ii) {
	      const double cos_theta = cos(double(ii)/double(ncopy+1)*M_PI*2.0); 
	      const double sin_theta = sin(double(ii)/double(ncopy+1)*M_PI*2.0); 
	      for (int ip = 0; ip < np_old; ++ip) {
		const int ip_new = ii*np_old + ip;
		part->pts->xp[ip_new][0] = xp_old[ip][0];
		part->pts->xp[ip_new][1] = xp_old[ip][1]*cos_theta + xp_old[ip][2]*sin_theta;
		part->pts->xp[ip_new][2] = xp_old[ip][2]*cos_theta - xp_old[ip][1]*sin_theta;
		part->pts->delta[ip_new] = delta_old[ip];
	      }
	    }
	    delete[] xp_old;
	    delete[] delta_old;
	  }
	  else if (token == "CYL_Z") {
	    const double angle_degrees = param->getDouble(iarg++);
	    if (mpi_rank == 0) cout << "copying points CYL_Z " << angle_degrees << endl;
	    assert(part->pts);
	    double (*xp_old)[3] = part->pts->xp; assert(xp_old);
	    double *delta_old = part->pts->delta; assert(delta_old);
	    const int np_old = part->pts->np;
	    part->pts->np_global *= ncopy+1;
	    part->pts->np *= ncopy+1;
	    part->pts->xp = new double[part->pts->np][3];
	    part->pts->delta = new double[part->pts->np];
	    for (int ii = 0; ii <= ncopy; ++ii) {
	      const double cos_theta = cos(double(ii)/double(ncopy+1)*M_PI*2.0); 
	      const double sin_theta = sin(double(ii)/double(ncopy+1)*M_PI*2.0); 
	      for (int ip = 0; ip < np_old; ++ip) {
		const int ip_new = ii*np_old + ip;
		part->pts->xp[ip_new][0] = xp_old[ip][0]*cos_theta + xp_old[ip][1]*sin_theta;
		part->pts->xp[ip_new][1] = xp_old[ip][1]*cos_theta - xp_old[ip][0]*sin_theta;
		part->pts->xp[ip_new][2] = xp_old[ip][2];
		part->pts->delta[ip_new] = delta_old[ip];
	      }
	    }
	    delete[] xp_old;
	    delete[] delta_old;
	  }
	  else {
	    CERR("NCOPY_PTS " << token << " not supported");
	  }
	}
        else if (token == "MASK_PTS") {
          assert(part->surface != NULL);
          assert(part->pts != NULL);
          double (*xp_old)[3] = part->pts->xp; assert(xp_old);
          double *delta_old = part->pts->delta; assert(delta_old);
          const int np_old = part->pts->np;
          const int8 np_global_old = part->pts->np_global;
          part->pts->np = 0;
          part->surface->initPointIsInside(xp_old,np_old);
          for (int ip = 0; ip < np_old; ++ip) {
            if (!part->surface->pointIsInside(xp_old[ip])) {
              delta_old[ip] = 0.0; // use delta sign as a flag...
            }
            else {
              ++part->pts->np;
            }
          }
          part->pts->xp = new double[part->pts->np][3];
          part->pts->delta = new double[part->pts->np];
          int ip = 0;
          for (int ip_old = 0; ip_old < np_old; ++ip_old) {
            if (delta_old[ip_old] > 0.0) {
              FOR_I3 part->pts->xp[ip][i] = xp_old[ip_old][i];
              part->pts->delta[ip] = delta_old[ip_old];
              ++ip;
            }
          }
          assert(ip == part->pts->np);
          delete[] xp_old;
          delete[] delta_old;
          int8 my_np = part->pts->np;
          MPI_Allreduce(&my_np,&part->pts->np_global,1,MPI_INT8,MPI_SUM,mpi_comm);
          if (mpi_rank == 0) {
            cout << " > MASK_PTS kept " << part->pts->np_global << " of " << np_global_old << " points" << endl;
          }

        }
	else if (token == "STRAND") {
          if (makeStrand(part,param,iarg)) {
            WUI(INFO,"STRAND built successfully");
          }
          else {
            WUI(WARN,"problem creating STRAND");
            ierr = -1;
          }
        }
        else if (token == "TRI_STRAND") {
          if (makeTriStrand(part,param,iarg)) {
            WUI(INFO,"TRI_STRAND mesh built");
          }
          else {
            WUI(WARN,"problem creating TRI_STRAND");
          }
        }
        else if (token == "AIRFOIL_STRAND") {
          if (makeAirfoilStrand(part,param,iarg)) {
            WUI(INFO,"AIRFOIL_STRAND built successfully");
          }
          else {
            WUI(WARN,"problem creating AIRFOIL_STRAND");
            ierr = -1;
          }
        }
	else if (token == "SVP_STRAND") {
          if (makeSvpStrand(part,param,iarg)) {
            WUI(INFO,"SVP_STRAND built successfully");
          }
          else {
            WUI(WARN,"problem creating SVP_STRAND");
            ierr = -1;
          }
        }
        else if (token == "EXTRUDE_SPLINE_WITH_STRAND") {
          if (makeExtrudeSplineWithStrand(part,param,iarg)) {
            WUI(INFO,"EXTRUDE_SPLINE_WITH_STRAND built successfully");
          }
          else {
            WUI(WARN,"problem creating EXTRUDE_SPLINE_WITH_STRAND");
            ierr = -1;
          }
        }
        else if (token == "CYLINDER_WITH_STRAND") {
          if (makeCylinderWithStrand(part,param,iarg)) {
            WUI(INFO,"CYLINDER_WITH_STRAND built successfully");
          }
          else {
            WUI(WARN,"problem creating CYLINDER_WITH_STRAND");
            ierr = -1;
          }
        }
        else if (token == "BOX_WITH_CART") {
          const bool pts_only = false;
          const int this_ierr = makeBoxWithCartPts(part,param,iarg,pts_only);
          if (this_ierr != 0) {
            WUI(WARN,"problem creating BOX_WITH_CART");
            ierr = -1;
          }
        }
        else if (token == "BOX_WITH_CART_PTS") {
          const bool pts_only = true;
          const int this_ierr = makeBoxWithCartPts(part,param,iarg,pts_only);
          if (this_ierr != 0) {
            WUI(WARN,"problem creating BOX_WITH_CART_PTS");
            ierr = -1;
          }
        }
        else if (token == "CYL_HEXCORE") {
          const bool pts_only = false;
          const int this_ierr = makeCylWithHexcorePts(part,param,iarg,pts_only);
          if (this_ierr != 0) {
            WUI(WARN,"problem creating CYL_HEXCORE");
            ierr = -1;
          }
        }
        else if (token == "CYL_HEXCORE_PTS") {
          const bool pts_only = true;
          const int this_ierr = makeCylWithHexcorePts(part,param,iarg,pts_only);
          if (this_ierr != 0) {
            WUI(WARN,"problem creating CYL_HEXCORE_PTS");
            ierr = -1;
          }
        }
        else if (token == "CYL_DLRT") {
          const bool pts_only = false;
          const int this_ierr = makeCylWithDlrtPts(part,param,iarg,pts_only);
          if (this_ierr != 0) {
            WUI(WARN,"problem creating CYL_DLRT");
            ierr = -1;
          }
        }
        else if (token == "CYL_DLRT_PTS") {
          const bool pts_only = true;
          const int this_ierr = makeCylWithDlrtPts(part,param,iarg,pts_only);
          if (this_ierr != 0) {
            WUI(WARN,"problem creating CYL_DLRT_PTS");
            ierr = -1;
          }
        }
        else if (token == "ANNULAR_CYL_PTS") {
          const bool pts_only = true;  // determines if boundaries are farfield or regular surface
          const int this_ierr = makeAnnularCylPts(part,param,iarg,pts_only);
          if (this_ierr != 0) {
            WUI(WARN,"problem creating ANNULAR_CYL_PTS");
            ierr = -1;
          }
        }
        else if (token == "ANNULAR_CYL") {
          const bool pts_only = false;  // determines if boundaries are farfield or regular surface
          const int this_ierr = makeAnnularCylPts(part,param,iarg,pts_only);
          if (this_ierr != 0) {
            WUI(WARN,"problem creating ANNULAR_CYL");
            ierr = -1;
          }
        }
        else if (token == "ANNULAR_DLRT_PTS") {
          const bool pts_only = true;  // determines if boundaries are farfield or regular surface
          const int this_ierr = makeAnnularDlrtPts(part,param,iarg,pts_only);
          if (this_ierr != 0) {
            WUI(WARN,"problem creating ANNULAR_DLRT_PTS");
            ierr = -1;
          }
        }
        else if (token == "ANNULAR_DLRT") {
          const bool pts_only = false;  // determines if boundaries are farfield or regular surface
          const int this_ierr = makeAnnularDlrtPts(part,param,iarg,pts_only);
          if (this_ierr != 0) {
            WUI(WARN,"problem creating ANNULAR_DLRT");
            ierr = -1;
          }
        }
        else if (token == "SPHERE_WITH_BL") {
          const int this_ierr = makeSphereWithBl(part,param,iarg);
          if (this_ierr != 0) {
            WUI(WARN,"problem creating SPHERE_WITH_BL");
            ierr = -1;
          }
        }
        else if (token == "POINTS_1D") {
          const int this_ierr = makePoints1D(part,param,iarg);
          if (this_ierr != 0) {
            WUI(WARN,"problem creating POINTS_1D");
            ierr = -1;
          }
        }
        else if (token == "MEDIAL_AXIS") {
          const int this_ierr = makeMedialAxis(part,param,iarg);
          if (this_ierr != 0) {
            WUI(WARN,"problem creating MEDIAL_AXIS");
            ierr = -1;
          }
        }
        else if (token == "MSH") {
          const string filename = param->getString(iarg++);
          if (mpi_rank == 0) {

            FluentMsh msh(filename);
            msh.writeTecplot("blah.dat");

            int (*noote)[4] = NULL;
            int nte = 0;
            msh.buildTets(nte,noote);

            double vol = 0.0;
            double vol_min = HUGE_VAL;
            double vol_max = -HUGE_VAL;
            const double (*const x_no)[3] = msh.x_no;
            for (int ite = 0; ite < nte; ++ite) {
              const double tet_vol = SIGNED_TET_VOLUME_6(x_no[noote[ite][0]],
                                                         x_no[noote[ite][1]],
                                                         x_no[noote[ite][2]],
                                                         x_no[noote[ite][3]]);
              vol_min = min(vol_min,tet_vol);
              vol_max = max(vol_max,tet_vol);
              vol += tet_vol;
            }

            cout << "GOT nte: " << nte << " total vol: " << vol << " tet vol_min,max: " << vol_min << " " << vol_max << endl;

            //GeomUtils::writeTecplot("blah.dat",spost,nst,xsp);

          }
          MPI_Pause("BLAH");
	  assert(0);
        }
        else if (token == "NLAYERS") {
          part->ff_nlayersString = param->getString(iarg++);
        }
        else if (token == "SCALE") {
          double dx[3];
          FOR_I3 dx[i] = param->getDouble(iarg++);
          if (part->surface) part->surface->scale(dx);
          if (part->pts) part->pts->scale(dx);
          PeriodicData::scalePeriodicTransformVec(dx);
        }
        else if (token == "TRANSLATE") {
          //assert(0); // TODO: need to translate any ddata_ff and ddata_solid, depending on ff_type/solid_type...
          double dx[3];
          FOR_I3 dx[i] = param->getDouble(iarg++);
          if (part->surface) part->surface->translate(dx);
          if (part->ff_surface) part->ff_surface->translate(dx);
          if (part->pts) part->pts->translate(dx);
        }
	else if (token == "ROTATE") {
          // TODO: need to rotate any ddata_ff and ddata_solid, depending on ff_type/solid_type...
          double point[3]; FOR_I3 point[i] = param->getDouble(iarg++);
          double axis[3]; FOR_I3 axis[i] = param->getDouble(iarg++);
          const double angle_deg = param->getDouble(iarg++);
          if (part->surface) part->surface->rotate(point,axis,angle_deg);
          if (part->ff_surface) part->ff_surface->rotate(point,axis,angle_deg);
          if (part->pts) part->pts->rotate(point,axis,angle_deg);
        }
        else if (token == "FLIP") {
          // change orientation of the tris by flipping last 2 nodes...
          if (part->surface) {
            if (mpi_rank_shared == 0) {
              for (int ist = 0; ist < part->surface->nst; ++ist) {
                const int tmp = part->surface->spost[ist][1];
                part->surface->spost[ist][1] = part->surface->spost[ist][2];
                part->surface->spost[ist][2] = tmp;
              }
            }
            MPI_Barrier(mpi_comm_shared);
          }
          else {
            WUI(WARN,"FLIP requires a surface");
            ierr = -1;
          }
        }
        else if (token == "FF") {
          token = MiscUtils::toUpperCase(param->getString(iarg++));
          if (token == "ALL") {
            part->ff_type = SURFACE_FAZONE_FF;
            assert(part->surface_zn_flag == NULL);
            part->surface_zn_flag = new int[part->surface->zoneVec.size()];
            for (int isz = 0; isz < part->surface->zoneVec.size(); ++isz) {
              if (mpi_rank == 0) cout << " > FAZONE " << part->surface->zoneVec[isz].getName() << " set to FF" << endl;
              part->surface_zn_flag[isz] = -1;
            }
          }
          else if ((token == "ZONE")||(token == "FAZONE")||(token == "ZONES")||(token == "FAZONES")) {
            // user should have supplied a comma-delimited list of zone names...
            part->ff_type = SURFACE_FAZONE_FF;
            assert(part->surface);
            token = param->getString(iarg++);
            set<string> zoneNameSet;
            MiscUtils::splitCsv(zoneNameSet,token);
            // default for the surface_zone_flag is 0...
            assert(part->surface_zn_flag == NULL);
            part->surface_zn_flag = new int[part->surface->zoneVec.size()];
            for (int isz = 0; isz < part->surface->zoneVec.size(); ++isz)
              part->surface_zn_flag[isz] = 0;
            // set -1 in every named zone...
            for (int isz = 0; isz < part->surface->zoneVec.size(); ++isz) {
              set<string>::iterator iter = zoneNameSet.find(part->surface->zoneVec[isz].getName());
              if (iter != zoneNameSet.end()) {
                // found in set...
                if (mpi_rank == 0) cout << " > FAZONE " << part->surface->zoneVec[isz].getName() << " set to FF" << endl;
                assert(part->surface_zn_flag[isz] == 0);
                part->surface_zn_flag[isz] = -1;
                zoneNameSet.erase(iter);
              }
            }
            // error checking...
            if (!zoneNameSet.empty()) {
              if (mpi_rank == 0) {
                WUI(WARN,"Warning: FF ZONE does not recognize the following zone names:");
                for (set<string>::iterator iter = zoneNameSet.begin(); iter != zoneNameSet.end(); ++iter)
                  WUI(WARN,"  " << *iter);
              }
              ierr = -1;
            }
          }
          else {
            WUI(WARN,"Unrecognized FF token: " << token << ". Expecting FF ZONE <zonename>, or FF GEOM <geom-data>");
            ierr = -1;
          }
        }
        else if (token == "DX") {
          // TODO: this is for the piston and valves: to build a static version of the moving
          // solver grid in stitch. Add motion parsing to this addPart routine???
          // the DX vector is used in certain motions (e.g. piston). It also sets the e0_* vector used
          // to determine the isInsideSolid funtion...
          part->b_dx = true;
          FOR_I3 part->dx[i] = param->getDouble(iarg++);
          if (mpi_rank == 0) cout << " > DX " << COUT_VEC(part->dx) << endl;
          // also modify the search direction for "isInsideSolid to the "dx" direction...
          const double mag_dx = MAG(part->dx);
          assert(mag_dx > 0.0);
          FOR_I3 part->e0_solid[i] = part->dx[i]/mag_dx;
          MiscUtils::getBestE1E2FromE0(part->e1_solid,part->e2_solid,part->e0_solid);
          // and set the ff to the same...
          FOR_I3 part->e0_ff[i] = part->e0_solid[i];
          FOR_I3 part->e1_ff[i] = part->e1_solid[i];
          FOR_I3 part->e2_ff[i] = part->e2_solid[i];
        }
	else if (token == "NO_WRITE_SURFACE") {
	  part->b_write_surface = false;
	  if (mpi_rank == 0) cout << " > surface will not be written into part file " << endl;
	}
	else if (token == "NO_WRITE_FF") {
	  part->b_write_ff = false;
	  if (mpi_rank == 0) cout << " > ff surface will not be written into part file " << endl;
	}
        else {
          if (mpi_rank == 0) cout << " > WARNING: skipping unrecognized PART token \"" << token << "\"" << endl;
        }
      }
    }

    if (ierr != 0)
      return false;

    // also, the first part should ALWAYS (we think?) have a surface...
    if (part->ipart_of_part == 0) {
      if (part->surface == NULL) {
        WUI(WARN,"The first PART should always have a surface");
        return false;
      }
    }

    return true;

  }

  void processInterface(Param * param,const bool b_help) {
    
    if (b_help) {
      assert(0);
      //helpProcessPart();
    }
    else {
      
      bool b_motion = false;
      bool b_rotation = false;
      bool b_axis = false;
      double axis[3];
      bool b_xc = false;
      double xc[3];

      int iarg = 0;
      const string zoneNamesCsv = param->getString(iarg++);
      
      while (iarg < param->size()) {
	string token = MiscUtils::toUpperCase(param->getString(iarg++));
	if (token == "MOTION") {
	  b_motion = true;
	}
	else if (token == "ROTATION") {
	  b_rotation = true;
	}
	else if (token == "XC") {
	  FOR_I3 xc[i] = param->getDouble(iarg++);
	  b_xc = true;
	}
	else if (token == "AXIS") {
	  FOR_I3 axis[i] = param->getDouble(iarg++);
	  b_axis = true;
	}
	else {
          if (mpi_rank == 0) cout << " > WARNING: skipping unrecognized INTERFACE token \"" << token << "\"" << endl;
        }
      }
      
      assert(b_motion);
      assert(b_rotation);
      assert(b_xc);
      assert(b_axis);
      
      // find the relevant surface...
      assert(partVec.size() == 1);
      assert(partVec[0]->surface != NULL);
      partVec[0]->surface->applyInterface(zoneNamesCsv,axis,xc);
      
    }
  }
  

  int makePartFromSbin(Part * part,const string& filename) {
    // *.sbin file...
    assert(part->surface == NULL);
    part->surface = new SurfaceShm(filename);
    assert(part->ff_type == UNKNOWN_FF);
    part->ff_type = NO_FF;
    assert(part->solid_type == UNKNOWN_SOLID);
    part->solid_type = SURFACE_SOLID;
    return 0;
  }

  int makePartFromMles(Part * part,const string& filename) {

    assert(part->surface == NULL);
    part->surface = new SurfaceShm();
    assert(part->ff_type == UNKNOWN_FF);
    part->ff_type = NO_FF;
    assert(part->ff_surface == NULL);
    assert(part->solid_type == UNKNOWN_SOLID);
    part->solid_type = SURFACE_SOLID;

    assert(part->pts == NULL);
    part->pts = new Points();

    MPI_Barrier(mpi_comm);
    double wtime0 = 0.0; 
    if (mpi_rank == 0)
      wtime0 = MPI_Wtime();

    int8 ncv_global;

    MPI_File fh; 
    char dummy[128];
    sprintf(dummy,"%s",filename.c_str());
    int ierr = MPI_File_open(mpi_comm,dummy,MPI_MODE_RDONLY,MPI_INFO_NULL,&fh);
    
    if (ierr != 0) { 
      //if (mpi_rank == 0) cout << "Error: cannot open restart file: " << filename << endl;
      return -1;
    } 
    
    int itmp[2] = { 0, 0 };
    if (mpi_rank == 0) { 
      // first 2 ints are: 0. magic number, 1. io version
      MPI_File_read(fh,itmp,2,MPI_INT,MPI_STATUS_IGNORE);
    }
    MPI_Bcast(itmp,2,MPI_INT,0,mpi_comm); 

    bool byte_swap = false;
    if (itmp[0] != UGP_IO_MAGIC_NUMBER) {
      ByteSwap::byteSwap(itmp,2);
      if (itmp[0] != UGP_IO_MAGIC_NUMBER) {
	if (mpi_rank == 0) cout << "Error: mles file does not start as expected." << endl;
	return -1;
      }
      if (mpi_rank == 0) cout << " > file requires byte swapping." << endl;
      byte_swap = true;
    }
    
    int io_version = itmp[1];
    if ((io_version != 3)&&(io_version != 5)) {
      if (mpi_rank == 0) cout << "Error: mles file version not 3 or 5: " << io_version << endl;
      return -1;
    }

    MPI_Offset offset = 8; // 2 ints
    Header header; 
    int done = 0;
    set<string> nameSet;

    int8 * cvora = NULL;
    int ncv;

    while (done != 1) {
      
      if (mpi_rank == 0) {
	MPI_File_read_at(fh,offset,&header,1,MPI_Header,MPI_STATUS_IGNORE); 
	if (byte_swap) ByteSwap::byteSwapHeader(&header,1);
      }
      MPI_Bcast(&header,1,MPI_Header,0,mpi_comm); 
      
      switch (header.id) {
	
      case UGP_IO_NO_FA_CV_COUNTS:

        assert(io_version == 3);
	assert(0);
        break;

      case UGP_IO_FA_CHECK:

	assert(io_version == 3);
	assert(0);
        break;

      case UGP_IO_FA_ZONE_HEADER: 

        assert(io_version == 3);
	assert(0);
        break;

      case UGP_IO_FA_ZONE:

        if (mpi_rank == 0) 
          cout << " > boundary and internal face zone: skipping..." << endl;
        
	break;

      case UGP_IO_NOOFA_I_AND_V:

        assert(io_version == 3);
	assert(0);
	break;

      case UGP_IO_CVOFA:

	if (mpi_rank == 0) 
	  cout << " > cvofa_global: skipping..." << endl;
	
	break;

      case UGP_IO_NO_FA_BF_CV_COUNTS:
	
	ncv_global    = ByteSwap::getInt8FromLswMswPair(header.idata+6);
	if (mpi_rank == 0) {
	  cout << " > ncv_global: " << ncv_global << endl;
 	}
	
	assert(cvora == NULL);
	MiscUtils::calcUniformDist(cvora,ncv_global,mpi_size);
	assert(cvora[mpi_rank+1]-cvora[mpi_rank] < TWO_BILLION);
	ncv = cvora[mpi_rank+1]-cvora[mpi_rank];

	part->pts->np_global = ncv_global;
	part->pts->np =  ncv;
	part->pts->xp = new double[part->pts->np][3];
	part->pts->delta = new double[part->pts->np];

	break;

      case UGP_IO_PERIODIC_TRANSFORM:
	
	if (mpi_rank == 0) cout << " > periodic transform(s): skipping... " << endl;

	break;

      case UGP_IO_BF_ZONE_HEADER:

	assert(header.idata[0] == FA_ZONE_BOUNDARY);
	assert(header.idata[1] == part->surface->zoneVec.size());

	// ensure the name is unique!... TODO do we need this (below it is not used for data lookup)
	{
	  string name = header.name;

          //if (nameSet.find(name) != nameSet.end()) {
          //  CERR("non-unique zone name:" << name);
          //}
	  while (nameSet.find(name) != nameSet.end())
	    name += "1";
	  nameSet.insert(name);
	  part->surface->zoneVec.push_back(name);
	}
	
	break;
	
      case UGP_IO_CVOBF_INT8:

	if (mpi_rank == 0) 
	  cout << " > cvobf_global: skipping..." << endl;
	
	break;
	
      case UGP_IO_NOOBF_I_AND_V_INT8:

	if (mpi_rank == 0) 
	  cout << " > noobf_i and noobf_v_global: skipping..." << endl;	
	
	break;
	
      case UGP_IO_SPOBF_I_V_WGT:
	
	if (mpi_rank == 0) 
	  cout << " > spobf_i/v/wgt: skipping..." << endl;
	
	break;

      case UGP_IO_SBOBF_I_AND_V:
	
	if (mpi_rank == 0) 
	  cout << " > sbobf_i and sbobf_v_global: skipping..." << endl;
	
	break;

      case UGP_IO_CVOFA_INT8:

	if (mpi_rank == 0) 
	  cout << " > cvofa_global: skipping..." << endl;
	
	break;

      case UGP_IO_NOOFA_I_AND_V_INT8:
	
	if (mpi_rank == 0) 
	  cout << " > noofa_i and noofa_v_global: skipping..." << endl;
		
	break;

      case UGP_IO_PBI_PNO:

	if (mpi_rank == 0) cout << " > pbi_pno: skipping..." << endl;
          
	break;

      case UGP_IO_X_NO:

	if (mpi_rank == 0) cout << " > x_no: skipping..." << endl;
	
	break;
	
      case UGP_IO_CV_ZONE_HEADER:

	break;
	
      case UGP_IO_CV_D1:

        if (io_version == 3) 
          assert(ncv_global == header.idata[0]);
        else 
          assert(ncv_global == ByteSwap::getInt8FromLswMswPair(header.idata+0));
	
	if ((strcmp(header.name,"vol_cv") == 0)||(strcmp(header.name,"VOL_VV") == 0)) {
	  if (mpi_rank == 0) cout << " > vol_cv: skipping..." << endl;
	}
        else if (strcmp(header.name,"r_vv") == 0) {

	  if (mpi_rank == 0) cout << " > reading r_vv..." << endl;
          readChunkedData<double>(fh,offset+header_size+cvora[mpi_rank]*double_size,part->pts->delta,ncv,byte_swap,mpi_comm);
	  FOR_ICV part->pts->delta[icv] *= 2.0;

	}
	else {
	  if (mpi_rank == 0) cout << " > unrecognized CV_D1: " << header.name << endl;
	}

	break;
	
      case UGP_IO_CV_D2:

        if (io_version == 3) {
          assert(ncv_global == header.idata[0]);
          assert(header.idata[1] == 3);
        }
        else {
          assert(ncv_global == ByteSwap::getInt8FromLswMswPair(header.idata+0));
          assert(header.idata[2] == 3);
        }
	assert(cvora != NULL);
	
	if ((strcmp(header.name,"x_vv") == 0)||(strcmp(header.name,"X_VV") == 0)) {
	  
	  if (mpi_rank == 0) cout << " > reading x_vv..." << endl;
          readChunkedData<double>(fh,offset+header_size+cvora[mpi_rank]*double_size*3,(double*)(part->pts->xp),ncv*3,byte_swap,mpi_comm);

	}
        else if (strcmp(header.name,"x_cv") == 0) {
	  if (mpi_rank == 0) cout << " > x_cv: skipping..." << endl;
	}
	else {
	  if (mpi_rank == 0) cout << " > unrecognized CV_D2: " << header.name << endl;
	}

	break;
	
      case UGP_IO_BF_D1:
	
	if (strcmp(header.name,"area_bf") == 0) {
	  if (mpi_rank == 0) cout << " > area_bf: skipping..." << endl;
	}
	else if (strcmp(header.name,"area_over_delta_bf") == 0) {
	  if (mpi_rank == 0) cout << " > area_over_delta_bf: skipping..." << endl;
	}
	else {
	  if (mpi_rank == 0) cout << " > unrecognized BF_D1: " << header.name << endl;
	}
	
	break;
	
      case UGP_IO_BF_D2:
	
	if (strcmp(header.name,"n_bf") == 0) {
	  if (mpi_rank == 0) cout << " > n_bf: skipping..." << endl;
	}
	else if (strcmp(header.name,"x_bf") == 0) {
	  if (mpi_rank == 0) cout << " > x_bf: skipping..." << endl;
	}
	else {
	  if (mpi_rank == 0) cout << " > unrecognized BF_D2: " << header.name << endl;
	}
	
	break;

      case UGP_IO_BF_DN33:

	if (strcmp(header.name,"Gij_bf") == 0) {
	  if (mpi_rank == 0) cout << " > Gij_bf: skipping..." << endl;
	}
	else {
	  if (mpi_rank == 0) cout << " > unrecognized BF_DN33: " << header.name << endl;
	}
	
	break;

      case UGP_IO_FA_D2:

	if (strcmp(header.name,"n_fa") == 0) {
	  if (mpi_rank == 0) cout << " > n_fa: skipping..." << endl;
	}
	else if (strcmp(header.name,"x_fa") == 0) {
	  if (mpi_rank == 0) cout << " > x_fa: skipping..." << endl;
	}
	else {
	  if (mpi_rank == 0) cout << " > unrecognized FA_D2: " << header.name << endl;
	}
	
	break;

      case UGP_IO_SURFACE:

	if (mpi_rank == 0) cout << " > surface (xp,spost,znost)..." << endl;

	assert(header.idata[0] > 0); // should have at least one tri? what about triple-periodic box?
	assert(header.idata[1] > 0);
        assert(part->surface->xsp == NULL);
        assert(part->surface->spost == NULL);
        assert(part->surface->znost == NULL);

	{
	  const int surface_nsp_global = header.idata[0]; 
	  const int surface_nst_global = header.idata[1];

	  part->surface->init(surface_nsp_global, surface_nst_global);

	  // recall surface is shm...
	  if (mpi_rank_shared == 0) {
#ifdef SERIAL_IO
            if (mpi_rank == 0) {
              MPI_File_read_at(fh,offset+header_size,(double*)(part->surface->xsp),surface_nsp_global*3,MPI_DOUBLE,MPI_STATUS_IGNORE);
              if ( byte_swap ) ByteSwap::byteSwap(part->surface->xsp,surface_nsp_global);

              MPI_File_read_at(fh,offset+header_size+surface_nsp_global*double_size*3,(int*)(part->surface->spost),surface_nst_global*3,MPI_INT,MPI_STATUS_IGNORE);
              if ( byte_swap) ByteSwap::byteSwap(part->surface->spost, surface_nst_global);

              MPI_File_read_at(fh,offset+header_size+surface_nsp_global*double_size*3+surface_nst_global*int_size*3,part->surface->znost,surface_nst_global,MPI_INT,MPI_STATUS_IGNORE);
              if ( byte_swap ) ByteSwap::byteSwap(part->surface->znost, surface_nst_global);
            }
            MPI_Bcast((double*)part->surface->xsp,surface_nsp_global*3,MPI_DOUBLE,0,mpi_comm_internode);
            MPI_Bcast((int*)part->surface->spost,surface_nst_global*3,MPI_INT,0,mpi_comm_internode);
            MPI_Bcast(part->surface->znost,surface_nst_global,MPI_INT,0,mpi_comm_internode);
#else
	    MPI_File_read_at(fh,offset+header_size,(double*)(part->surface->xsp),surface_nsp_global*3,MPI_DOUBLE,MPI_STATUS_IGNORE);
	    if ( byte_swap ) ByteSwap::byteSwap(part->surface->xsp,surface_nsp_global);

	    MPI_File_read_at(fh,offset+header_size+surface_nsp_global*double_size*3,(int*)(part->surface->spost),surface_nst_global*3,MPI_INT,MPI_STATUS_IGNORE);
	    if ( byte_swap) ByteSwap::byteSwap(part->surface->spost, surface_nst_global);

	    MPI_File_read_at(fh,offset+header_size+surface_nsp_global*double_size*3+surface_nst_global*int_size*3,part->surface->znost,surface_nst_global,MPI_INT,MPI_STATUS_IGNORE);
	    if ( byte_swap ) ByteSwap::byteSwap(part->surface->znost, surface_nst_global);
#endif

	  }
	  MPI_Barrier(mpi_comm_shared);
	}

	break;

      case UGP_IO_HASHIDS:

        //b_foundHash = true;
        //read two hash id's from mles file. First id identifies
        //the current file, store in mlesHash.  Ignore the second
        //field for now.  TODO include surface hash as parent
	/*
	  RestartHashUtilities::mlesReadHashes(fh, offset, header); //will bCast to all ranks
	  if (mpi_rank==0){
          cout << " > mles hash: " << RestartHashUtilities::mlesHash << endl;
	  }
	*/

        break;

      case UGP_IO_EOF:
        
	done = 1;
        break;
	
      default:
	
	if (mpi_rank == 0) cout << " > unknown header: " << header.id << " \"" << header.name << "\", skipping." << endl;
	
      }
      
      offset += header.skip;
    
    }

    MPI_File_close(&fh);
    
    MPI_Barrier(mpi_comm);
    if (mpi_rank == 0) {
      const double seconds = MPI_Wtime() - wtime0;
      cout << " > read_mles summary: ncv_global: " << ncv_global << 
	" file size: " << double(offset)/1.0E+9 << " [GB] read rate: " << double(offset)/1.0E+9/seconds << " [GB/s]" << endl;
    }

    delete[] cvora;

    return 0;
  }

  void writeUsedPoints(const string& filename) {
    COUT1("writeUsedPoints(): " << filename);

    // mirror countPoints method to determine who is inside
    // create temporary contiguous memory space for write

    // get hcp points...
    // PartData::hcpPts->ensurePoints();

    // start with the hcp points...
    int8 my_count[2] = {0,0};
    my_count[0] += PartData::hcpPts->np;

    vector<pair<int,int> > part_pts;  // ipart,ip

    // need to add points added from parts...
    PartData::prepareIsInside();
    for (int ipart = 0; ipart < PartData::partVec.size(); ++ipart) {
      if (PartData::partVec[ipart]->pts) {
        for (int ip = 0; ip < PartData::partVec[ipart]->pts->np; ++ip) {
          // we are going to include this point in the build, if and only if it
          // is NOT inside the FF of every part later than us in the list (these parts take priority)
          // AND it is not inside the solid of any part. Since we have checked the FF of the
          // later parts, we can just check the solid of earlier parts...
          bool in = true;
          for (int ipart2 = ipart+1; ipart2 < PartData::partVec.size(); ++ipart2) {
            if (PartData::partVec[ipart2]->hasFF()) {
              if (PartData::partVec[ipart2]->isInsideFF(PartData::partVec[ipart]->pts->xp[ip])) {
                in = false;
                break;
              }
            }
            else {
              if (PartData::partVec[ipart2]->isInsideSolid(PartData::partVec[ipart]->pts->xp[ip],false)) {
                in = false;
                break;
              }
            }
          }
          if (in) {
            for (int ipart2 = 0; ipart2 < ipart; ++ipart2) {
              if (PartData::partVec[ipart2]->isInsideSolid(PartData::partVec[ipart]->pts->xp[ip],ipart2==0)) { // note flag true when ipart2 is the first part
                in = false;
                break;
              }
            }
          }
          if (in) {
            part_pts.push_back(pair<int,int> (ipart,ip));
            ++my_count[1];
          }
        }
      }
    }
    assert(my_count[1] == int8(part_pts.size()));
    int8 buf[2];
    MPI_Allreduce(my_count,buf,2,MPI_INT8,MPI_SUM,mpi_comm);
    if (mpi_rank == 0) {
      const int8 count = buf[0]+buf[1];
      cout << " > global points written: " << count << " (hcp: " << buf[0] << ", part: " << buf[1] << ")" << endl;
    }

    const int8 np_global_write = buf[0]+buf[1];

    // now actually write points
    MiscUtils::mkdir_for_file_collective(filename,0);
    string tmp_filename = MiscUtils::makeTmpPrefix(filename);
    char dummy[128]; assert(tmp_filename.length() < 128);
    sprintf(dummy,"%s",tmp_filename.c_str());
    MPI_File_delete(dummy,MPI_INFO_NULL);

    MPI_File fh;
    MPI_File_open(mpi_comm,dummy,MPI_MODE_WRONLY|MPI_MODE_CREATE,MPI_INFO_NULL,&fh);

    if (mpi_rank == 0) {
      int8 ibuf[4] = { POINTS_IO_MAGIC_NUMBER, POINTS_IO_VERSION, np_global_write, 1 }; // last is include delta or not maybe?
      MPI_File_write_at(fh,0,ibuf,4,MPI_INT8,MPI_STATUS_IGNORE);
    }

    int8 my_np = my_count[0] + my_count[1];
    int8 my_disp;
    MPI_Scan(&my_np,&my_disp,1,MPI_INT8,MPI_SUM,mpi_comm);
    assert( (mpi_rank != mpi_size-1) || (my_disp == np_global_write) );
    MPI_Offset offset = (my_disp-my_np)*3*8 + 4*8;
    if (buf[0]) {
      writeChunkedData<double>(fh,offset,(double *)PartData::hcpPts->xp,3*my_count[0],mpi_comm);
    }
    if (buf[1]) {
      // build temp array to store for writing
      double (*xp_tmp)[3] = new double[my_count[1]][3];
      for (vector<pair<int,int> >::iterator it=part_pts.begin(); it!=part_pts.end(); ++it) {
        FOR_I3 xp_tmp[it-part_pts.begin()][i] = PartData::partVec[it->first]->pts->xp[it->second][i];
      }
      writeChunkedData<double>(fh,offset+3*8*my_count[0],(double *)xp_tmp,3*my_count[1],mpi_comm);
      DELETE(xp_tmp);
    }


    // delta...
    offset = np_global_write*3*8 + 4*8 + (my_disp-my_np)*8;
    if (buf[0]) {
      writeChunkedData<double>(fh,offset,PartData::hcpPts->delta,my_count[0],mpi_comm);
    }
    if (buf[1]) {
      // build temp array to store for writing
      double * delta_tmp = new double[my_count[1]];
      for (vector<pair<int,int> >::iterator it=part_pts.begin(); it!=part_pts.end(); ++it) {
        FOR_I3 delta_tmp[it-part_pts.begin()] = PartData::partVec[it->first]->pts->delta[it->second];
      }
      writeChunkedData<double>(fh,offset+8*my_count[0],delta_tmp,my_count[1],mpi_comm);
      DELETE(delta_tmp);
    }

    // and close...

    MPI_File_close(&fh);
    if (mpi_rank == 0) {
      remove(filename.c_str());
      rename(tmp_filename.c_str(),filename.c_str());
    }

  }

  void processWritePbin(Param * param,const bool b_help) {
    if (b_help) {
      helpProcessWritePbin();
    }
    else {
      try {
        string filename = "points.pbin";
        int iarg = 0;
        while (iarg < param->size() && iarg < 1) {
          string token = param->getString(iarg++);
          // assume this is a filename...
          filename = token;
        }
        if (hcpPts == NULL) throw(1);
        PartData::writeUsedPoints(filename);
      }
      catch(int e) {
        if (e == 1) {
          WUI(WARN,"expecting NSMOOTH <int> prior to WRITE_PBIN <string>");
        }
        else {
          WUI(WARN,"expecting WRITE_PBIN <string>");
        }
        helpProcessWritePbin();
      }
    }
  }
  void helpProcessWritePbin() {
    WUI(INFO,
        "WRITE_PBIN writes the final static mesh points in pbin format. Some examples:\n" <<
        "WRITE_PBIN # write to points.pbin\n" <<
        "WRITE_PBIN my_points.pbin\n" <<
        "for more detail see [$CWIKB:stitch_export]");
  }

  void processWriteHcpPts(Param * param,const bool b_help) {
    if (b_help) {
      helpProcessWriteHcpPts();
    }
    else {
      try {
        string filename = "hcp_points.pbin";
        int iarg = 0;
        while (iarg < param->size() && iarg < 1) {
          string token = param->getString(iarg++);
          // assume this is a filename...
          filename = token;
        }
        if (hcpPts == NULL) throw(1);
        PartData::hcpPts->writeBinary(filename);
      }
      catch(int e) {
        if (e == 1) {
          WUI(WARN,"expecting NSMOOTH <int> prior to WRITE_HCP_PTS <string>");
        }
        else {
          WUI(WARN,"expecting WRITE_HCP_PTS <string>");
        }
        helpProcessWriteHcpPts();
      }
    }
  }
  void helpProcessWriteHcpPts() {
    WUI(INFO,
        "WRITE_HCP_PTS writes HCP generated points in pbin format. Some examples:\n" <<
        "WRITE_HCP_PTS # write to hcp_points.pbin\n" <<
        "WRITE_HCP_PTS my_points.pbin\n" <<
        "for more detail see [$CWIKB:stitch_export]");
  }

  /*
    void rmPart(Param * param) {

    // removes a part based on name...
    //
    // RMPART bl
    // RMPART NAME=bl
    //

    int iarg = 0;
    string name = param->getString(iarg++);
    if (name == "NAME")
    name = param->getString(iarg++);

    // find the index of the part...

    int ipart;
    for (ipart = 0; ipart < partVec.size(); ++ipart)
    if (partVec[ipart]->name == name)
    break;

    if (ipart == partVec.size()) {
    if (mpi_rank == 0) cout << "Warning: RMPART: part not found with name \"" << name << "\". Nothing removed." << endl;
    return;
    }

    if (mpi_rank == 0) cout << "RMPART: found part \"" << name << "\". Removing." << endl;

    partVec[ipart]->clear();

    // copy down any parts above ipart: these are just Part pointer copies...
    for (int ipart2 = ipart+1; ipart2 < partVec.size(); ++ipart2)
    partVec[ipart2-1] = partVec[ipart2];
    partVec.resize(partVec.size()-1);

    // go through surface_st_flag for all other parts and clear any references to this part,
    // and re-index references for parts that have changed...
    if (mpi_rank_shared == 0) {
    for (int ipart2 = 0; ipart2 < partVec.size(); ++ipart2) {
    if (partVec[ipart2]->surface) {
    for (int ist = 0; ist < partVec[ipart2]->surface->nst; ++ist) {
    int this_ipart,this_igr;
    if (partVec[ipart2]->getPartAndGroupForSt(this_ipart,this_igr,ist)) {
    if (this_ipart == ipart) {
    partVec[ipart2]->clearPartAndGroupForSt(ist);
    }
    else if (this_ipart > ipart) {
    assert(0); // what are we doing here?...
    // this reset avoids the warning...
    partVec[ipart2]->clearPartAndGroupForSt(ist);
    partVec[ipart2]->setPartAndGroupForSt(this_ipart-1,this_igr,ist);
    }
    }
    }
    }
    }
    }
    MPI_Barrier(mpi_comm_shared);

    }
  */

  enum GroupKind {
    UNTESTED_GROUP_KIND,
    UNKNOWN_GROUP_KIND,
    CYLINDER_GROUP_KIND,
    LOGICAL_QUAD_GROUP_KIND,
  };

  class GroupData {
  private:
    GroupKind kind;
    int idata[16];
    double ddata[16];
  public:
    double n[3];
    double x[3];
    double area;
    int igr;
    int ipart;
    int nlink;
    int spoli_s;
    int * spoli_i; // surface-point-of-link references isp on ipart
    int * spoli_v; // "
    double * l_li;     // linear length of link
    double (*n_li)[3]; // normal (link is closed)
    double (*x_li)[3]; // edge-weighted center
    double * d_li;     // mean distance // TODO: why do this extra work at this point?
    double * drms_li;  // rms distance
    double spline_crease_tol;
    // other data that gets filled in by certain calls...
    vector<SplineStuff::CubicSpline> csplineVec;
    GroupData() {
      igr = -1;
      ipart = -1; // TODO: why is this needed?
      nlink = 0;
      spoli_s = 0;
      spoli_i = NULL;
      spoli_v = NULL;
      l_li = NULL;
      n_li = NULL;
      x_li = NULL;
      d_li = NULL;
      drms_li = NULL;
      kind = UNTESTED_GROUP_KIND;
      setSplineCreaseTol(45.0);
    }
    void clear() {
      DELETE(spoli_i);
      DELETE(spoli_v);
      DELETE(l_li);
      DELETE(n_li);
      DELETE(x_li);
      DELETE(d_li);
      DELETE(drms_li);
      csplineVec.clear();
      kind = UNTESTED_GROUP_KIND;
    }
    void dump() {
      cout << "GroupData::dump() nlink: " << nlink << endl;
    }
    GroupKind getKind() const { return kind; }
    void setKind(const GroupKind kind) { this->kind = kind; }
    void setSplineCreaseTol(const double angle_d) {
      // specified angle in degrees
      spline_crease_tol = cos(angle_d/180.0*M_PI);
    }
    //================
    // Right Cylinder
    //================
    bool isCylinder(const double tol) {
      // on subsequent calls to this routine, no need to re-compute...
      if (kind != UNTESTED_GROUP_KIND)
        return (kind == CYLINDER_GROUP_KIND);
      // note that tol is treated as an absolute tolerance when +ve, and
      // a relative tolerance when negative (as a fraciton of r_avg)...
      if (nlink != 2) return false;
      // here we use the link centroids as the cylinder axis...
      double e0[3] = DIFF(x_li[1],x_li[0]);
      const double L = MAG(e0);
      FOR_I3 e0[i] /= L;
      SurfaceShm * surface = partVec[ipart]->surface;
      assert(surface);
      // for ALL tris, get the nodal r, rmin, rmax...
      double my_r_avg[2] = { 0.0, 0.0 };
      double my_r_minmax[2] = { HUGE_VAL, HUGE_VAL };
      for (int ist = mpi_rank; ist < surface->nst; ist += mpi_size) {
        int igr_;
        if (partVec[ipart]->getGroupForSt(igr_,ist)) {
          if (igr_ == igr) {
            // compute the area of this tri for weighting...
            const double n2[3] = TRI_NORMAL_2(surface->xsp[surface->spost[ist][0]],
                                              surface->xsp[surface->spost[ist][1]],
                                              surface->xsp[surface->spost[ist][2]]);
            const double wgt = MAG(n2);
            // nodes first...
            FOR_I3 {
              const int isp = surface->spost[ist][i];
              double dxp[3] = DIFF(surface->xsp[isp],x_li[0]);
              const double dp = DOT_PRODUCT(dxp,e0);
              FOR_I3 dxp[i] -= dp*e0[i];
              const double this_r = MAG(dxp);
              my_r_avg[0] += wgt; // only use the nodes for averaging
              my_r_avg[1] += wgt*this_r;
              my_r_minmax[0] = min(my_r_minmax[0],this_r);
              my_r_minmax[1] = min(my_r_minmax[1],-this_r);
            }
            // then edges to determine the effect of segmenting...
            FOR_I3 {
              const int isp0 = surface->spost[ist][i];
              const int isp1 = surface->spost[ist][(i+1)%3];
              double dxp[3]; FOR_I3 dxp[i] = 0.5*(surface->xsp[isp0][i]+surface->xsp[isp1][i])-x_li[0][i];
              const double dp = DOT_PRODUCT(dxp,e0);
              FOR_I3 dxp[i] -= dp*e0[i];
              const double this_r = MAG(dxp);
              my_r_minmax[0] = min(my_r_minmax[0],this_r);
              my_r_minmax[1] = min(my_r_minmax[1],-this_r);
            }
          }
        }
      }
      double r_avg[2];
      MPI_Allreduce(my_r_avg,r_avg,2,MPI_DOUBLE,MPI_SUM,mpi_comm);
      r_avg[1] /= r_avg[0];
      double r_minmax[2];
      MPI_Allreduce(my_r_minmax,r_minmax,2,MPI_DOUBLE,MPI_MIN,mpi_comm);
      // also, check how orthogonal the ends of the cylinder are...
      // do this everywhere, unless we find a case where it makes sense
      // to do this in parallel...
      double dL_buf[4] = { HUGE_VAL, HUGE_VAL, HUGE_VAL, HUGE_VAL };
      for (int il = 0; il < nlink; ++il) {
        for (int sol = spoli_i[il]; sol != spoli_i[il+1]; ++sol) {
          const int isp = spoli_v[sol];
          const double dxp[3] = DIFF(surface->xsp[isp],x_li[il]);
          const double dp = DOT_PRODUCT(dxp,e0);
          dL_buf[il*2] = min(dL_buf[il*2],dp);
          dL_buf[il*2+1] = min(dL_buf[il*2+1],-dp);
        }
      }
      bool r_fails,l_fails;
      if (tol < 0.0) {
        // use -ve tol to indicate a relative tolerance based on a fraction of r...
        if (mpi_rank == 0) {
          cout << " > got r_avg: " << r_avg[1] << ", r_min,max: " << r_minmax[0] << " " << -r_minmax[1] << " dr_min/r_avg: " << (r_minmax[0]-r_avg[1])/r_avg[1] << ", dr_max/r_avg " << (-r_minmax[1]-r_avg[1])/r_avg[1] << endl;
          cout << " > got dL_max/r_avg at x0: " << max(dL_buf[0],-dL_buf[1])/r_avg[1] << ", dL_max/r_avg at x1: " << max(dL_buf[2],-dL_buf[3])/r_avg[1] << endl;
        }
        r_fails = ((fabs(r_minmax[0]-r_avg[1]) > -tol*r_avg[1])||(fabs(-r_minmax[1]-r_avg[1]) > -tol*r_avg[1]));
        l_fails = ((fabs(dL_buf[0]) > -tol*r_avg[1])||(fabs(dL_buf[1]) > -tol*r_avg[1])||(fabs(dL_buf[2]) > -tol*r_avg[1])||(fabs(dL_buf[3]) > -tol*r_avg[1]));
      }
      else {
        assert(tol > 0.0);
        if (mpi_rank == 0) {
          cout << " > got r_avg: " << r_avg[1] << ", r_min,max: " << r_minmax[0] << " " << -r_minmax[1] << " dr_min/tol: " << (r_minmax[0]-r_avg[1])/tol << ", dr_max/tol " << (-r_minmax[1]-r_avg[1])/tol << endl;
          cout << " > got dL_max/r_avg at x0: " << max(dL_buf[0],-dL_buf[1])/r_avg[1] << ", dL_max/r_avg at x1: " << max(dL_buf[2],-dL_buf[3])/r_avg[1] << endl;
          cout << " > got dL at x0: " << dL_buf[0] << " " << -dL_buf[1] << ", at x1: " << dL_buf[2] << " " << -dL_buf[3] << endl;
        }
        r_fails = ((fabs(r_minmax[0]-r_avg[1]) > tol)||(fabs(-r_minmax[1]-r_avg[1]) > tol));
        l_fails = ((fabs(dL_buf[0]) > tol)||(fabs(dL_buf[1]) > tol)||(fabs(dL_buf[2]) > tol)||(fabs(dL_buf[3]) > tol));
      }
      // store the L and the avg radius in real data...
      ddata[0] = L;
      ddata[1] = r_avg[1];
      if (!r_fails) {
        kind = CYLINDER_GROUP_KIND;
        // if this is a right-cylinder, put a 1 into idata[0]....
        if (!l_fails) {
          idata[0] = 1;
        }
        else {
          idata[0] = 0;
        }
        return true;
      }
      return false;
    }

    bool isRightCylinder() {
      // you should only call this once you know you have a cylinder...
      assert(kind == CYLINDER_GROUP_KIND);
      return idata[0] == 1;
    }

    double getCylinderL() const {
      assert(kind == CYLINDER_GROUP_KIND);
      return ddata[0];
    }

    double getCylinderR() const {
      assert(kind == CYLINDER_GROUP_KIND);
      return ddata[1];
    }

    void countRightCylinderBl(int counts[5],const double dn,const double dt,const int n) {
      assert(kind == CYLINDER_GROUP_KIND);
      assert(isRightCylinder());
      assert(dt>dn);
      // pull data
      const double L = ddata[0];
      const double r = ddata[1];
      const double sf = pow(dt/dn,1.0/double(n-1));
      const double delta = (dt*sf-dn)/(sf-1.0);
      if (mpi_rank == 0) {
        cout << " > L: " << L << endl;
        cout << " > BL sf: " << sf << " delta for " << n << " layers: " << delta << endl;
      }
      assert(L > 2.0*delta);
      const int nL = int((L-2.0*delta)/dt); // + 2*n
      const int nTheta = int(2.0*M_PI*r/dt);
      if (mpi_rank == 0) cout << " > total point estimate: nL+2*n: " << nL+2*n << " x nTheta: " << nTheta << " x n: " << n << " = " << (nL+2*n)*nTheta*n << endl;
      // no surface counts...
      counts[0] = counts[1] = 0;
      // ff_surface point count...
      assert(nlink == 2);
      counts[2] = spoli_i[2] + nTheta*(4*n-2);
      // ff_surface tri count...
      counts[3] = spoli_i[2] + nTheta*4*(2*n-1);
      // pts->np
      counts[4] = 0;
      int rank_flag = 0;
      for (int i = 0; i < n; ++i) {
        // the number of points in the azimuthal direction depends on the spacing...
        const double dx = dn*pow(sf,i);
        const double this_nTheta = int(2.0*M_PI*r/dx);
        if (i == n-1) assert(this_nTheta == nTheta);
        for (int k = 0; k < this_nTheta; ++k) {
          for (int j = 0; j <= i; ++j) {
            if (rank_flag == mpi_rank)
              ++counts[4];
            //++rank_flag;
            //if (rank_flag == mpi_size)
            //  rank_flag = 0;
          }
        }
        ++rank_flag;
        if (rank_flag == mpi_size)
          rank_flag = 0;
      }
      for (int i = 0; i < nL; ++i) {
        for (int k = 0; k < nTheta; ++k) {
          for (int j = 0; j < n; ++j) {
            if (rank_flag == mpi_rank)
              ++counts[4];
            //++rank_flag;
            //if (rank_flag == mpi_size)
            //  rank_flag = 0;
          }
        }
        ++rank_flag;
        if (rank_flag == mpi_size)
          rank_flag = 0;
      }
      for (int i = n-1; i >= 0; --i) {
        const double dx = dn*pow(sf,i);
        const double this_nTheta = int(2.0*M_PI*r/dx);
        if (i == n-1) assert(this_nTheta == nTheta);
        for (int k = 0; k < this_nTheta; ++k) {
          for (int j = 0; j <= i; ++j) {
            if (rank_flag == mpi_rank)
              ++counts[4];
            //++rank_flag;
            //if (rank_flag == mpi_size)
            //  rank_flag = 0;
          }
        }
        ++rank_flag;
        if (rank_flag == mpi_size)
          rank_flag = 0;
      }
    }
    void makeRightCylinderBl(Part * part,int disp[5],const double dn,const double dt,const int n) {
      assert(kind == CYLINDER_GROUP_KIND);
      assert(isRightCylinder());
      double e0[3] = DIFF(x_li[1],x_li[0]);
      const double mag_e0 = MAG(e0);
      FOR_I3 e0[i] /= mag_e0;
      double e1[3],e2[3];
      MiscUtils::getBestE1E2FromE0(e1,e2,e0);
      const double L = ddata[0];
      const double r = ddata[1];
      const double sf = pow(dt/dn,1.0/double(n-1));
      const double delta = (dt*sf-dn)/(sf-1.0);
      const int nL = int((L-2.0*delta)/dt); // + 2*n
      const int nTheta = int(2.0*M_PI*r/dt);
      assert(part->ff_surface);
      if (mpi_rank_shared == 0) {
        // ===================================================
        // surface points: xsp_ff...
        // ===================================================
        const int nsp_ff = spoli_i[2] + nTheta*(4*n-2);
        // copy the links into xsp_ff first...
        for (int sol = 0; sol < spoli_i[2]; ++sol) {
          const int isp = spoli_v[sol];
          FOR_I3 part->ff_surface->xsp[disp[2]+sol][i] = partVec[ipart]->surface->xsp[isp][i];
        }
        // then the layers start at spoli_i[1]...
        for (int k = 0; k <= 2*(n-1); ++k) {
          for (int j = 0; j < nTheta; ++j) {
            double xp;
            double rp;
            if (k < n) {
              xp = dn*(pow(sf,k)-1.0)/(sf-1.0);
              rp = r - dn*(pow(sf,k+1)-1.0)/(sf-1.0);
            }
            else {
              const int kp = 2*(n-1)-k;
              xp = L - dn*(pow(sf,kp+1)-1.0)/(sf-1.0);
              rp = r - dn*(pow(sf,kp+1)-1.0)/(sf-1.0);
            }
            //if (j == 0) cout << "p0 k: " << k << " xp: " << xp << " L-xp: " << L-xp << endl;
            const double yp = rp*cos(2.0*M_PI*double(j)/double(nTheta));
            const double zp = rp*sin(2.0*M_PI*double(j)/double(nTheta));
            int isp_ff = spoli_i[2] + nTheta*2*k+j;
            assert(isp_ff < nsp_ff);
            FOR_I3 part->ff_surface->xsp[disp[2]+isp_ff][i] = xp*e0[i] + yp*e1[i] + zp*e2[i];
            if (k < n-1) {
              xp = dn*(pow(sf,k+1)-1.0)/(sf-1.0);
            }
            else {
              const int kp = 2*(n-1)-k;
              xp = L-dn*(pow(sf,kp)-1.0)/(sf-1.0);
            }
            //if (j == 0) cout << "p1 k: " << k << " xp: " << xp << " L-xp: " << L-xp << endl;
            isp_ff = spoli_i[2] + nTheta*(2*k+1)+j;
            assert(isp_ff < nsp_ff);
            FOR_I3 part->ff_surface->xsp[disp[2]+isp_ff][i] = xp*e0[i] + yp*e1[i] + zp*e2[i];
          }
        }
        disp[2] += nsp_ff;
        // ===================================================
        // surface tris: xst_ff...
        // ===================================================
        const int nst_ff = spoli_i[2] + nTheta*4*(2*n-1);
        // 1. the lip tris are special because they have to respect both discretizations...
        int ist_ff = 0;
        {
          int isp_ff_prev = spoli_i[1]-2;
          int isp_ff = spoli_i[1]-1;
          // compute the theta of the middle...
          const double dx[3] = {
            0.5*(part->ff_surface->xsp[isp_ff_prev][0]+part->ff_surface->xsp[isp_ff][0]) - x_li[0][0],
            0.5*(part->ff_surface->xsp[isp_ff_prev][1]+part->ff_surface->xsp[isp_ff][1]) - x_li[0][1],
            0.5*(part->ff_surface->xsp[isp_ff_prev][2]+part->ff_surface->xsp[isp_ff][2]) - x_li[0][2] };
          double yp = DOT_PRODUCT(dx,e1);
          double zp = DOT_PRODUCT(dx,e2);
          double theta = atan2(zp,yp);
          int j_prev = (int)floor(theta/(2.0*M_PI)*double(nTheta)+0.5);
          if (j_prev < 0) j_prev += nTheta;
          assert((j_prev >= 0)&&(j_prev < nTheta));
          isp_ff_prev = isp_ff;
          for (isp_ff = spoli_i[0]; isp_ff < spoli_i[1]; ++isp_ff) {
            // compute the theta...
            const double dx[3] = {
              0.5*(part->ff_surface->xsp[isp_ff_prev][0]+part->ff_surface->xsp[isp_ff][0]) - x_li[0][0],
              0.5*(part->ff_surface->xsp[isp_ff_prev][1]+part->ff_surface->xsp[isp_ff][1]) - x_li[0][1],
              0.5*(part->ff_surface->xsp[isp_ff_prev][2]+part->ff_surface->xsp[isp_ff][2]) - x_li[0][2] };
            double yp = DOT_PRODUCT(dx,e1);
            double zp = DOT_PRODUCT(dx,e2);
            double theta = atan2(zp,yp);
            int j = (int)floor(theta/(2.0*M_PI)*double(nTheta)+0.5);
            if (j < 0) j += nTheta;
            assert((j >= 0)&&(j < nTheta));
            // tris from the nTheta side...
            while (j_prev != j) {
              int j_next = j_prev + 1;
              if (j_next == nTheta) j_next = 0;
              part->ff_surface->spost[ist_ff+disp[3]][0] = spoli_i[2]+j_next;
              part->ff_surface->spost[ist_ff+disp[3]][1] = isp_ff_prev;
              part->ff_surface->spost[ist_ff+disp[3]][2] = spoli_i[2]+j_prev;
              part->ff_surface_dxost[ist_ff+disp[3]] = dn;
              ++ist_ff;
              j_prev = j_next;
            }
            // tris from the lipline side...
            part->ff_surface->spost[ist_ff+disp[3]][0] = isp_ff_prev;
            part->ff_surface->spost[ist_ff+disp[3]][1] = spoli_i[2]+j;
            part->ff_surface->spost[ist_ff+disp[3]][2] = isp_ff;
            part->ff_surface_dxost[ist_ff+disp[3]] = dn;
            ++ist_ff;
            j_prev = j;
            isp_ff_prev = isp_ff;
          }
          assert(ist_ff == spoli_i[1]+nTheta);
        }
        // now the step tris...
        for (int k = 0; k <= 2*(n-1); ++k) {
          double dxr;
          double dxx;
          if (k < (n-1)) {
            dxr = dn*pow(sf,k);
            dxx = dn*pow(sf,k+1);
          }
          else {
            const int kp = 2*(n-1)-k;
            dxr = dn*pow(sf,kp);
            dxx = dn*pow(sf,kp);
          }
          int j_prev = nTheta-1;
          for (int j = 0; j < nTheta; ++j) {
            // 2 tris at the same r...
            assert(ist_ff < nst_ff);
            part->ff_surface->spost[ist_ff+disp[3]][0] = spoli_i[2]+nTheta*2*k+j_prev;
            part->ff_surface->spost[ist_ff+disp[3]][1] = spoli_i[2]+nTheta*(2*k+1)+j;
            part->ff_surface->spost[ist_ff+disp[3]][2] = spoli_i[2]+nTheta*2*k+j;
            part->ff_surface_dxost[ist_ff+disp[3]] = dxr;
            ++ist_ff;
            assert(ist_ff < nst_ff);
            part->ff_surface->spost[ist_ff+disp[3]][0] = spoli_i[2]+nTheta*2*k+j_prev;
            part->ff_surface->spost[ist_ff+disp[3]][1] = spoli_i[2]+nTheta*(2*k+1)+j_prev;
            part->ff_surface->spost[ist_ff+disp[3]][2] = spoli_i[2]+nTheta*(2*k+1)+j;
            part->ff_surface_dxost[ist_ff+disp[3]] = dxr;
            ++ist_ff;
            // 2 tris at the same x...
            if (k < 2*(n-1)) {
              assert(ist_ff < nst_ff);
              part->ff_surface->spost[ist_ff+disp[3]][0] = spoli_i[2]+nTheta*(2*k+1)+j_prev;
              part->ff_surface->spost[ist_ff+disp[3]][1] = spoli_i[2]+nTheta*(2*k+2)+j;
              part->ff_surface->spost[ist_ff+disp[3]][2] = spoli_i[2]+nTheta*(2*k+1)+j;
              part->ff_surface_dxost[ist_ff+disp[3]] = dxx;
              ++ist_ff;
              assert(ist_ff < nst_ff);
              part->ff_surface->spost[ist_ff+disp[3]][0] = spoli_i[2]+nTheta*(2*k+1)+j_prev;
              part->ff_surface->spost[ist_ff+disp[3]][1] = spoli_i[2]+nTheta*(2*k+2)+j_prev;
              part->ff_surface->spost[ist_ff+disp[3]][2] = spoli_i[2]+nTheta*(2*k+2)+j;
              part->ff_surface_dxost[ist_ff+disp[3]] = dxx;
              ++ist_ff;
            }
            // copy j into j_prev for next one...
            j_prev = j;
          }
        }
        // and the second ring of lip tris...
        {
          int isp_ff_prev = spoli_i[2]-2;
          int isp_ff = spoli_i[2]-1;
          // compute the theta of the middle...
          const double dx[3] = {
            0.5*(part->ff_surface->xsp[isp_ff_prev][0]+part->ff_surface->xsp[isp_ff][0]) - x_li[0][0],
            0.5*(part->ff_surface->xsp[isp_ff_prev][1]+part->ff_surface->xsp[isp_ff][1]) - x_li[0][1],
            0.5*(part->ff_surface->xsp[isp_ff_prev][2]+part->ff_surface->xsp[isp_ff][2]) - x_li[0][2] };
          double yp = DOT_PRODUCT(dx,e1);
          double zp = DOT_PRODUCT(dx,e2);
          double theta = atan2(zp,yp);
          int j_prev = (int)floor(theta/(2.0*M_PI)*double(nTheta)+0.5);
          if (j_prev < 0) j_prev += nTheta;
          assert((j_prev >= 0)&&(j_prev < nTheta));
          isp_ff_prev = isp_ff;
          for (isp_ff = spoli_i[1]; isp_ff < spoli_i[2]; ++isp_ff) {
            // compute the theta...
            const double dx[3] = {
              0.5*(part->ff_surface->xsp[isp_ff_prev][0]+part->ff_surface->xsp[isp_ff][0]) - x_li[0][0],
              0.5*(part->ff_surface->xsp[isp_ff_prev][1]+part->ff_surface->xsp[isp_ff][1]) - x_li[0][1],
              0.5*(part->ff_surface->xsp[isp_ff_prev][2]+part->ff_surface->xsp[isp_ff][2]) - x_li[0][2] };
            double yp = DOT_PRODUCT(dx,e1);
            double zp = DOT_PRODUCT(dx,e2);
            double theta = atan2(zp,yp);
            int j = (int)floor(theta/(2.0*M_PI)*double(nTheta)+0.5);
            if (j < 0) j += nTheta;
            assert((j >= 0)&&(j < nTheta));
            // tris from the nTheta side...
            while (j_prev != j) {
              int j_next = j_prev - 1;
              if (j_next == -1) j_next = nTheta-1;
              assert(ist_ff < nst_ff);
              part->ff_surface->spost[ist_ff+disp[3]][0] = spoli_i[2]+(4*n-3)*nTheta+j_next;
              part->ff_surface->spost[ist_ff+disp[3]][1] = isp_ff_prev;
              part->ff_surface->spost[ist_ff+disp[3]][2] = spoli_i[2]+(4*n-3)*nTheta+j_prev;
              part->ff_surface_dxost[ist_ff+disp[3]] = dn;
              ++ist_ff;
              j_prev = j_next;
            }
            // tris from the lipline side...
            assert(ist_ff < nst_ff);
            part->ff_surface->spost[ist_ff+disp[3]][0] = isp_ff_prev;
            part->ff_surface->spost[ist_ff+disp[3]][1] = spoli_i[2]+(4*n-3)*nTheta+j;
            part->ff_surface->spost[ist_ff+disp[3]][2] = isp_ff;
            part->ff_surface_dxost[ist_ff+disp[3]] = dn;
            ++ist_ff;
            j_prev = j;
            isp_ff_prev = isp_ff;
          }
        }
        assert(ist_ff = nst_ff);
        disp[3] += nst_ff;
      }
      MPI_Barrier(mpi_comm_shared);

      // ==================================================================
      // now the points...
      // ==================================================================
      assert(part->pts);
      int rank_flag = 0;
      for (int i = 0; i < n; ++i) {
        const double xp = dn*(pow(sf,0.5+double(i))-1.0)/(sf-1.0);
        // the number of points in the azimuthal direction depends on the spacing...
        const double dx = dn*pow(sf,i);
        const double this_nTheta = int(2.0*M_PI*r/dx);
        if (i == n-1) assert(this_nTheta == nTheta);
        for (int k = 0; k < this_nTheta; ++k) {
          double theta = (double(k)+0.5)/double(this_nTheta)*2.0*M_PI;
          for (int j = 0; j <= i; ++j) {
            double rp = r - dn*(pow(sf,0.5+double(j))-1.0)/(sf-1.0);
            double yp = rp*cos(theta);
            double zp = rp*sin(theta);
            if (rank_flag == mpi_rank) {
              part->pts->xp[disp[4]][0] = xp*e0[0] + yp*e1[0] + zp*e2[0];
              part->pts->xp[disp[4]][1] = xp*e0[1] + yp*e1[1] + zp*e2[1];
              part->pts->xp[disp[4]][2] = xp*e0[2] + yp*e1[2] + zp*e2[2];
              part->pts->delta[disp[4]] = 1.75*dx; // i.e. a little more than sqrt(3)
              ++disp[4];
            }
            // NOTE: we moved it to round robin just on i because it helped speed up imaging by reducing
            // the global pixel buffer from all ranks.
            //++rank_flag;
            //if (rank_flag == mpi_size)
            //  rank_flag = 0;
          }
        }
        ++rank_flag;
        if (rank_flag == mpi_size)
          rank_flag = 0;
      }
      for (int i = 0; i < nL; ++i) {
        double xp = delta + (double(i)+0.5)/double(nL)*(L-2.0*delta);
        for (int k = 0; k < nTheta; ++k) {
          double theta = (double(k)+0.5)/double(nTheta)*2.0*M_PI;
          for (int j = 0; j < n; ++j) {
            double rp = r - dn*(pow(sf,0.5+double(j))-1.0)/(sf-1.0);
            double yp = rp*cos(theta);
            double zp = rp*sin(theta);
            if (rank_flag == mpi_rank) {
              part->pts->xp[disp[4]][0] = xp*e0[0] + yp*e1[0] + zp*e2[0];
              part->pts->xp[disp[4]][1] = xp*e0[1] + yp*e1[1] + zp*e2[1];
              part->pts->xp[disp[4]][2] = xp*e0[2] + yp*e1[2] + zp*e2[2];
              part->pts->delta[disp[4]] = (L-2.0*delta)/nL*1.75;
              ++disp[4];
            }
            //++rank_flag;
            //if (rank_flag == mpi_size)
            //  rank_flag = 0;
          }
        }
        ++rank_flag;
        if (rank_flag == mpi_size)
          rank_flag = 0;
      }
      for (int i = n-1; i >= 0; --i) {
        double xp = L - dn*(pow(sf,0.5+double(i))-1.0)/(sf-1.0);
        const double dx = dn*pow(sf,i);
        const double this_nTheta = int(2.0*M_PI*r/dx);
        if (i == n-1) assert(this_nTheta == nTheta);
        for (int k = 0; k < this_nTheta; ++k) {
          double theta = (double(k)+0.5)/double(this_nTheta)*2.0*M_PI;
          for (int j = 0; j <= i; ++j) {
            double rp = r - dn*(pow(sf,0.5+double(j))-1.0)/(sf-1.0);
            double yp = rp*cos(theta);
            double zp = rp*sin(theta);
            if (rank_flag == mpi_rank) {
              part->pts->xp[disp[4]][0] = xp*e0[0] + yp*e1[0] + zp*e2[0];
              part->pts->xp[disp[4]][1] = xp*e0[1] + yp*e1[1] + zp*e2[1];
              part->pts->xp[disp[4]][2] = xp*e0[2] + yp*e1[2] + zp*e2[2];
              part->pts->delta[disp[4]] = 1.75*dx; // i.e. a little more than sqrt(3)
              ++disp[4];
            }
            //++rank_flag;
            //if (rank_flag == mpi_size)
            //  rank_flag = 0;
          }
        }
        ++rank_flag;
        if (rank_flag == mpi_size)
          rank_flag = 0;
      }
    }

    //==================
    // Logical Quad
    //==================
    bool isLogicalQuad() {
      // if the work has already been done, then skip...
      if (kind != UNTESTED_GROUP_KIND)
        return kind == LOGICAL_QUAD_GROUP_KIND;
      // do the work...
      if (nlink != 1) return false;
      if (spoli_s < 4) return false; // need atleast 4 points
      // got through edges looking for corners (set by spline_crease_tol)
      SurfaceShm * surface = partVec[ipart]->surface;
      int cnt = 0;
      // we also store the 4 approximate edge lengths in ddata[0:3]...
      FOR_I4 ddata[i] = 0.0;
      for (int sol = spoli_i[0]; sol != spoli_i[1]; ++sol) {
        const int isp = spoli_v[sol];
        int isp_prev,isp_next;
        if (sol == spoli_i[0]) {
          isp_prev = spoli_v[spoli_i[1]-1];
          isp_next = spoli_v[sol+1];
        }
        else if (sol == spoli_i[1]-1) {
          isp_prev = spoli_v[sol-1];
          isp_next = spoli_v[spoli_i[0]];
        }
        else {
          isp_prev = spoli_v[sol-1];
          isp_next = spoli_v[sol+1];
        }
        const double dx0[3] = DIFF(surface->xsp[isp],surface->xsp[isp_prev]);
        const double dx1[3] = DIFF(surface->xsp[isp_next],surface->xsp[isp]);
        const double dx0_mag = MAG(dx0);
        const double dx1_mag = MAG(dx1);
        const double dp = DOT_PRODUCT(dx0,dx1);
        // add lengths to ddata[0:3]. We shift this addition so the last length is
        // set first, that way the first link is from node idata[0]->idata[1], etc.
        // (see below sketch)...
        ddata[(cnt+3)%4] += dx0_mag; // i.e. 3,0,1,2,3
        if (dp < spline_crease_tol*dx0_mag*dx1_mag) {
          if (cnt == 4)
            return false;
          idata[cnt++] = sol;
        }
      }
      if (cnt < 4)
        return false;
      assert(cnt == 4);
      // adjust orientation so long dimension is aligned as shown in the below figure...
      if (ddata[0] + ddata[2] < ddata[1] + ddata[3]) {
        const int itmp    = idata[3];
        const double dtmp = ddata[3];
        for (int i = 3; i > 0; --i) {
          idata[i] = idata[i-1];
          ddata[i] = ddata[i-1];
        }
        idata[0] = itmp;
        ddata[0] = dtmp;
      }
      kind = LOGICAL_QUAD_GROUP_KIND;
      return true;
    }

    void countLogicalQuadBl(int counts[5],const double dn,const double dt1,const double dt2,const int n,const int mode) {

      assert(dn > 0.0);
      assert(dt1 > dn);
      assert(dt2 > dn);
      assert(kind == LOGICAL_QUAD_GROUP_KIND);

      // no surface counts...
      counts[0] = counts[1] = 0;

      // the ff_surface gets built over the logical quad using transfinite interpolation from edges.
      //
      // build 4 curves arount the logical quad oriented with the
      // longest dimension along the x'-direction as follows...
      // (arrows indicate the direction of increasing "s")...
      //
      //         <-- spline[2]
      //    3-------------------2
      //    |                   | ^
      //    | spline[3]         | spline[1]
      //    | \/                |
      //    0-------------------1
      //          spline[0]-->
      //
      // estimate the distance along the surface edge between each pair of sols using the edge distance
      // first 3 edges are not split...

      // recall that isLogicalQuad set idata[0,1,2,3] contains the indices of the corners in spoli_v (formerly sols[4]),
      // and distances in ddata[0,1,2,3]...

      SurfaceShm * surface = partVec[ipart]->surface;
      assert(surface);

      // set the 4 edge-based splines with direction as defined in the figure above...
      assert(csplineVec.empty());
      csplineVec.resize(4);

      {
        double (*x)[3] = new double[spoli_s][3]; // definately big enough...
        FOR_I4 {
          int count = 0;
          int sol = idata[i];
          while (sol != idata[(i+1)%4]) {
            const int isp = spoli_v[sol];
            FOR_J3 x[count][j] = surface->xsp[isp][j];
            ++count;
            ++sol;
            if (sol == spoli_i[1])
              sol = 0;
          }
          // add the end...
          const int isp = spoli_v[sol];
          FOR_J3 x[count][j] = surface->xsp[isp][j];
          ++count;
          csplineVec[i].init(x,count);
        }
        delete[] x;
      }

      // try transfinite interpolation...

      // now the grid...

      const double L1 = 0.5*(csplineVec[0].getLength()+csplineVec[2].getLength());
      const double L2 = 0.5*(csplineVec[1].getLength()+csplineVec[3].getLength());
      const int n1 = int(L1/dt1); assert(n1 > 0);
      const int n2 = int(L2/dt2); assert(n2 > 0);

      if (mpi_rank == 0) {
        cout << " > logical quad dimensions: spline L1 " << L1 << " (from line segments: " << 0.5*(ddata[0]+ddata[2]) <<
          "), spline L2: " <<  L2 << " (from line segments: " << 0.5*(ddata[1]+ddata[3]) << ")" << endl;
      }

      if (mode == 0) {

        // ==============================================================================
        // mode 0 maintains the same dn across each layer and has a single dt...
        // ==============================================================================

        if (mpi_rank == 0)
          cout << " > n1: " << n1 << " n2: " << n2 << " n: " << n << endl;

        counts[2] = spoli_s+4 + 2*(n-1)*(n1+n2) + (n1+1)*(n2+1); // surface point count: +4 is because corners get duplicated
        counts[3] = spoli_s + 2*n1+2*n2 + 4*(n-1)*(n1+n2) + 2*n1*n2; // surface tri count
        counts[4] = 0; // local point count
        int rank_flag = 0; // for point round-robin
        for (int i = 0; i < n1; ++i) {
          for (int j = 0; j < n2; ++j) {
            if (mpi_rank == rank_flag)
              counts[4] += n;
            ++rank_flag;
            if (rank_flag == mpi_size)
              rank_flag = 0;
          }
        }

        if (mpi_rank == 0)
          cout << " > ff_nsp: " << counts[2] << " ff_nst: " << counts[3] << " points->np_global: " << n1*n2*n << endl;

      }
      else if (mode == 1) {

        // ==============================================================================
        // mode 1 is a regular square surrounding the cartesian points with point count
        // increasing away from the perimiter until target dn is achieved in the middle...
        // ==============================================================================

        if (mpi_rank == 0) cout << " > n1: " << n1 << " n2: " << n2 << endl;

        // the number of ff_surface nodes is equal to all the nodes along the perimiter
        // plus nodes at the top of the box...

        counts[2] = spoli_s + (n1+1)*(n2+1) + 4; // point count: +4 is because corners get duplicated
        counts[3] = spoli_s + 2*n1 + 2*n2 + 2*n1*n2; // tri count

        // for points....

        const double h = double(n)*0.5*(dt1+dt2);
        const double delta = h;
        const double dt = delta/double(n);

        double dn_current = dt;
        int np = n;
        int i = 0, j = 0;
        counts[4] = 0;
        int rank_flag = 0; // for point round-robin: note rank 0 per node only...
        while (1) {
          // add the perimiter times the current np to the points count...
          const int i_f = i;
          const int i_l = n1-i-1;
          const int j_f = j;
          const int j_l = n2-j-1;
          if ((i_f == i_l)||(i_f+1 == i_l)||(j_f == j_l)||(j_f+1 == j_l)) {
            // last row(s)/column(s)...
            if ((mpi_rank_shared == 0)&&(mpi_rank_internode == rank_flag))
              counts[4] += (i_l-i_f+1)*(j_l-j_f+1)*np;
            break;
          }
          if ((mpi_rank_shared == 0)&&(mpi_rank_internode == rank_flag))
            counts[4] += (2*(i_l-i_f+1) + 2*(j_l-j_f+1) - 4)*np;
          ++i;
          ++j;
          // modify rank_flag for next region...
          ++rank_flag;
          if (rank_flag == mpi_size_internode)
            rank_flag = 0;
          // add another point if we have not reached the target dn...
          if (dn_current > dn*1.01) {
            // add another point to the strand...
            ++np;
            // and recompute dn_current...
            // this starts us on the correct side of the function, half-way towards the
            // zero gradient...
            dn_current = 0.5*exp(log(dt*double(np-1)/(delta-dt))/double(np))*(delta-dt)/double(np-1);
            int iter = 0;
            while (1) {
              ++iter;
              if (iter > 30) {
                if (mpi_rank == 0) cout << "ERROR: cannot converge dn/delta system: delta,dt,np: " << delta << " " << dt << " " << np << endl;
                assert(0);
              }
              const double t1 = pow(dt/dn_current,double(np)/double(np-1));
              const double t2 = pow(dt/dn_current,1.0/double(np-1));
              const double f = dn_current*t1-delta*t2+delta-dn_current;
              const double fp = (delta*t2-dn_current*t1-dn_current*double(np-1))/(dn_current*double(np-1));
              //if (mpi_rank == 0) cout << "XXX np = " << np << " got dn_current: " << dn_current << " f: " << f << " fp: " << fp << " f/fp/dt: " << f/fp/dt << endl;
              if (f > 0.5*dn_current*fp)
                dn_current *= 0.5; // revert to bisection for stability
              else {
                dn_current -= f/fp;
                if (fabs(f/fp) < 1.0E-12*dt)
                  break;
              }
            }
            //const double sf = pow(dt/dn_current,1.0/double(np-1));
            //const double delta_check = (dt*sf-dn_current)/(sf-1.0);
            //if (mpi_rank == 0) cout << " XXXX got sf: " << sf << " dn_current: " << dn_current << " delta_check: " << delta_check << " delta: " << delta << endl;
          }
        }

        if (mpi_rank == 0) cout << " > np_max (in center): " << np << endl;

      }
      else {

        assert(0);

      }

    }

#define IJ(I,J) (J)*(n1+1)+(I)

    inline int getIspFromIJK(const int i, const int j, const int k,const int n1,const int n2,const int n) const {
      assert((k >= 0)&&(k < n));
      assert((j >= 0)&&(j <= n2));

      if (k == 0) {
        assert((i >= 0)&&(i <= n1));
        return IJ(i,j);
      }
      else {
        int base = (n1+1)*(n2+1)+2*(n1+n2)*(k-1);
        if (j == 0) {
          assert((i >= 0)&&(i <= n1));
          return base+i;
        }
        else if (j == n2) {
          assert((i >= 0)&&(i <= n1));
          return base+(n1+1)+2*(n2-1)+i;
        }
        else if (i == 0) {
          return base+(n1+1)+2*(j-1);
        }
        else if (i == n1) {
          return base+(n1+1)+2*(j-1)+1;
        }
        else {
          return -1; // not on perimeter
        }
      }

    }

    void buildTriAdtAndVecForGroup(Adt<double>& adt,vector<int>& istVec,const int igr) {

      SurfaceShm * surface = partVec[ipart]->surface;

      for (int ist = 0; ist < surface->nst; ++ist) {
        int igr_;
        if (partVec[ipart]->getGroupForSt(igr_,ist)) {
          if (igr_ == igr) {
            istVec.push_back(ist);
          }
        }
      }
      const int ntri = istVec.size();
      double (*bbmin)[3] = new double[ntri][3];
      double (*bbmax)[3] = new double[ntri][3];
      for (int itri = 0; itri < ntri; ++itri) {
        const int ist = istVec[itri];
        FOR_I3 bbmin[itri][i] = HUGE_VAL;
        FOR_I3 bbmax[itri][i] = -HUGE_VAL;
        FOR_I3 {
          const double * xsp_ = surface->xsp[surface->spost[ist][i]];
          FOR_J3 bbmin[itri][j] = min(bbmin[itri][j],xsp_[j]);
          FOR_J3 bbmax[itri][j] = max(bbmax[itri][j],xsp_[j]);
        }
      }
      adt.initAdt(ntri,bbmin,bbmax);
      delete[] bbmin;
      delete[] bbmax;
    }

    void buildTransfiniteInterpPoints(double (*xsp_ff)[3],const int n1, const int n2,const double xws[3],const double xes[3],const double xen[3],const double xwn[3],const vector<int>& istVec,Adt<double>& adt) {
      // place points in the xsp_ff array
      SurfaceShm * surface = partVec[ipart]->surface;
      double d2_max = 0.0;
      vector<int> intVec;
      for (int i = 0; i <= n1; ++i) {
        const double xi = double(i)/double(n1);
        // coordinate values along the south and north edge...
        double xs[3]; csplineVec[0].getX(xs,xi);
        double xn[3]; csplineVec[2].getX(xn,1.0-xi);
        for (int j = 0; j <= n2; ++j) {
          const double eta = double(j)/double(n2);
          // coordinate values along the west and east edge...
          double xw[3]; csplineVec[3].getX(xw,1.0-eta);
          double xe[3]; csplineVec[1].getX(xe,eta);
          // transfinite interp...
          FOR_K3 xsp_ff[IJ(i,j)][k] = (1.0-eta)*xs[k] + eta*xn[k] + (1.0-xi)*xw[k] + xi*xe[k] -
            ( xi*eta*xen[k] + xi*(1.0-eta)*xes[k] + (1.0-xi)*eta*xwn[k] + (1.0-xi)*(1.0-eta)*xws[k] );
          // if this is not an edge (or corner) point, then project to the surface...
          if ((i > 0)&&(i < n1)&&(j > 0)&&(j < n2)) {
            // this is an internal point, so project...
            assert(intVec.empty());
            adt.buildListForClosestPoint(intVec,xsp_ff[IJ(i,j)]);
            assert(!intVec.empty());
            double d2_min = HUGE_VAL;
            double xp_min[3];
            for (int ii = 0; ii < intVec.size(); ++ii) {
              const int ist = istVec[intVec[ii]];
              double xp[3]; MiscUtils::getClosestPointOnTriRobust(xp,xsp_ff[IJ(i,j)],
                                                                  surface->xsp[surface->spost[ist][0]],
                                                                  surface->xsp[surface->spost[ist][1]],
                                                                  surface->xsp[surface->spost[ist][2]]);
              const double d2 = DIST2(xp,xsp_ff[IJ(i,j)]);
              if (d2 < d2_min) {
                d2_min = d2;
                FOR_K3 xp_min[k] = xp[k];
              }
            }
            intVec.clear();
            assert(d2_min < HUGE_VAL);
            FOR_K3 xsp_ff[IJ(i,j)][k] = xp_min[k];
            d2_max = max(d2_max,d2_min);
          }
        }
      }
      if (mpi_rank == 0) cout << " > max distance of original transfinite points from surface: " << sqrt(d2_max) << endl;
    }

    double smoothAndReprojectInterpPoints(double (*xsp_ff)[3],const int smooth_iters,const int n1, const int n2,const vector<int>& istVec,Adt<double>& adt) {
      SurfaceShm * surface = partVec[ipart]->surface;
      double (*xsp0_ff)[3] = new double[(n1+1)*(n2+1)][3];

      double d2_max = 0.0;
      for (int iter = 1; iter <= smooth_iters; ++iter) {
        // copy xsp into xsp0...
        for (int i = 0; i <= n1; ++i) {
          for (int j = 0; j <= n2; ++j) {
            FOR_K3 xsp0_ff[IJ(i,j)][k] = xsp_ff[IJ(i,j)][k];
          }
        }
        // now smooth and re-project...
        vector<int> intVec;
        for (int i = 1; i < n1; ++i) {
          for (int j = 1; j < n2; ++j) {
            FOR_K3 xsp_ff[IJ(i,j)][k] =
              0.25*xsp0_ff[IJ(i,j)][k] +
              0.125*(xsp0_ff[IJ(i-1,j)][k]+xsp0_ff[IJ(i+1,j)][k]+xsp0_ff[IJ(i,j-1)][k]+xsp0_ff[IJ(i,j+1)][k]) +
              0.0625*(xsp0_ff[IJ(i-1,j-1)][k]+xsp0_ff[IJ(i+1,j-1)][k]+xsp0_ff[IJ(i-1,j+1)][k]+xsp0_ff[IJ(i+1,j+1)][k]);
            assert(intVec.empty());
            adt.buildListForClosestPoint(intVec,xsp_ff[IJ(i,j)]);
            assert(!intVec.empty());
            double d2_min = HUGE_VAL;
            double xp_min[3];
            for (int ii = 0; ii < intVec.size(); ++ii) {
              const int ist = istVec[intVec[ii]];
              double xp[3]; MiscUtils::getClosestPointOnTriRobust(xp,xsp_ff[IJ(i,j)],
                                                                  surface->xsp[surface->spost[ist][0]],
                                                                  surface->xsp[surface->spost[ist][1]],
                                                                  surface->xsp[surface->spost[ist][2]]);
              const double d2 = DIST2(xp,xsp_ff[IJ(i,j)]);
              if (d2 < d2_min) {
                d2_min = d2;
                FOR_K3 xp_min[k] = xp[k];
              }
            }
            intVec.clear();
            assert(d2_min < HUGE_VAL);
            d2_max = max(d2_max,DIST2(xsp_ff[IJ(i,j)],xp_min));
            FOR_K3 xsp_ff[IJ(i,j)][k] = xp_min[k];
          }
        }
        /*
          if (iter%10 == 0) {
          char filename[128];
          sprintf(filename,"pts.%08d.dat",iter);
          FILE * fp = fopen(filename,"w");
          for (int i = 0; i <= n1; ++i) {
          for (int j = 0; j <= n2; ++j) {
          fprintf(fp,"%18.15le %18.15le %18.15le\n",xsp[IJ(i,j)][0],xsp[IJ(i,j)][1],xsp[IJ(i,j)][2]);
          }
          }
          fclose(fp);
          }
        */
        if (mpi_rank == 0) cout << " > smooth-and-reproject iter: " << iter << " d: " << sqrt(d2_max) << endl;
      }
      delete[] xsp0_ff;
      return d2_max;
    }

    void computeAndSmoothNormals(double (*xsp0_ff)[3],double const (*xsp_ff)[3],const int n1,const int n2,const int niters) {
      for (int i = 0; i <= n1; ++i) {
        for (int j = 0; j <= n2; ++j) {
          FOR_K3 xsp0_ff[IJ(i,j)][k] = 0.0;
        }
      }
      for (int i = 0; i < n1; ++i) {
        for (int j = 0; j < n2; ++j) {
          // consider the rectangle i->i+1, j->j+1...
          const double na[3] = TRI_NORMAL_2(xsp_ff[IJ(i,j)],xsp_ff[IJ(i+1,j)],xsp_ff[IJ(i+1,j+1)]);
          const double nb[3] = TRI_NORMAL_2(xsp_ff[IJ(i,j)],xsp_ff[IJ(i+1,j+1)],xsp_ff[IJ(i,j+1)]);
          FOR_K3 xsp0_ff[IJ(i,j)][k] += na[k]+nb[k];
          FOR_K3 xsp0_ff[IJ(i+1,j)][k] += na[k]+nb[k];
          FOR_K3 xsp0_ff[IJ(i,j+1)][k] += na[k]+nb[k];
          FOR_K3 xsp0_ff[IJ(i+1,j+1)][k] += na[k]+nb[k];
        }
      }
      if (mpi_rank == 0) cout << " > performing " << niters << " smoothing iterations on normal..." << endl;

      // smooth the normals...
      double (*nsp_smooth)[3] = new double[(n1+1)*(n2+1)][3];
      for (int iter = 0; iter < niters; ++iter) {
        // normalize ALL...
        for (int i = 0; i <= n1; ++i) {
          for (int j = 0; j <= n2; ++j) {
            const double mag = MAG(xsp0_ff[IJ(i,j)]);
            FOR_K3 nsp_smooth[IJ(i,j)][k] = xsp0_ff[IJ(i,j)][k]/mag;
          }
        }
        // smooth interior...
        for (int i = 1; i < n1; ++i) {
          for (int j = 1; j < n2; ++j) {
            // just leave this un-normalized...
            FOR_K3 xsp0_ff[IJ(i,j)][k] =
              2.0*(nsp_smooth[IJ(i-1,j)][k]+nsp_smooth[IJ(i+1,j)][k]+nsp_smooth[IJ(i,j-1)][k]+nsp_smooth[IJ(i,j+1)][k]) +
              (nsp_smooth[IJ(i-1,j-1)][k]+nsp_smooth[IJ(i+1,j-1)][k]+nsp_smooth[IJ(i-1,j+1)][k]+nsp_smooth[IJ(i+1,j+1)][k]);
          }
        }
        // and smooth edges...
        for (int i = 1; i < n1; ++i) {
          FOR_K3 xsp0_ff[IJ(i,0)][k] = (nsp_smooth[IJ(i-1,0)][k]+nsp_smooth[IJ(i+1,0)][k]);
          FOR_K3 xsp0_ff[IJ(i,n2)][k] = (nsp_smooth[IJ(i-1,n2)][k]+nsp_smooth[IJ(i+1,n2)][k]);
        }
        for (int j = 1; j < n2; ++j) {
          FOR_K3 xsp0_ff[IJ(0,j)][k] = (nsp_smooth[IJ(0,j-1)][k]+nsp_smooth[IJ(0,j+1)][k]);
          FOR_K3 xsp0_ff[IJ(n1,j)][k] = (nsp_smooth[IJ(n1,j-1)][k]+nsp_smooth[IJ(n1,j+1)][k]);
        }
      }
      delete[] nsp_smooth;
    }

    void computeAndSmoothInPlaneNormals(double (*xsp0_ff)[3],double const (*xsp_ff)[3],const int n1,const int n2,const int niters,double * _n_plane) {

      // compute wgts of normal components to keep
      double n_plane[3];
      FOR_I3 n_plane[i] = _n_plane[i];
      NORMALIZE(n_plane);
      const double in_plane_wgts[3] = {1.0-fabs(n_plane[0]),1.0-fabs(n_plane[1]),1.0-fabs(n_plane[2])};

      for (int i = 0; i <= n1; ++i) {
        for (int j = 0; j <= n2; ++j) {
          FOR_K3 xsp0_ff[IJ(i,j)][k] = 0.0;
        }
      }
      for (int i = 0; i < n1; ++i) {
        for (int j = 0; j < n2; ++j) {
          // consider the rectangle i->i+1, j->j+1...
          double na[3] = TRI_NORMAL_2(xsp_ff[IJ(i,j)],xsp_ff[IJ(i+1,j)],xsp_ff[IJ(i+1,j+1)]);
          double nb[3] = TRI_NORMAL_2(xsp_ff[IJ(i,j)],xsp_ff[IJ(i+1,j+1)],xsp_ff[IJ(i,j+1)]);
          FOR_K3 {
            na[k] *= in_plane_wgts[k];
            nb[k] *= in_plane_wgts[k];
          }
          FOR_K3 xsp0_ff[IJ(i,j)][k] += na[k]+nb[k];
          FOR_K3 xsp0_ff[IJ(i+1,j)][k] += na[k]+nb[k];
          FOR_K3 xsp0_ff[IJ(i,j+1)][k] += na[k]+nb[k];
          FOR_K3 xsp0_ff[IJ(i+1,j+1)][k] += na[k]+nb[k];
        }
      }
      if (mpi_rank == 0) cout << " > performing " << niters << " smoothing iterations on normal..." << endl;

      // smooth the normals...
      double (*nsp_smooth)[3] = new double[(n1+1)*(n2+1)][3];
      for (int iter = 0; iter < niters; ++iter) {
        // normalize ALL...
        for (int i = 0; i <= n1; ++i) {
          for (int j = 0; j <= n2; ++j) {
            const double mag = MAG(xsp0_ff[IJ(i,j)]);
            FOR_K3 nsp_smooth[IJ(i,j)][k] = xsp0_ff[IJ(i,j)][k]/mag;
          }
        }
        // smooth interior...
        for (int i = 1; i < n1; ++i) {
          for (int j = 1; j < n2; ++j) {
            // just leave this un-normalized...
            FOR_K3 xsp0_ff[IJ(i,j)][k] =
              2.0*(nsp_smooth[IJ(i-1,j)][k]+nsp_smooth[IJ(i+1,j)][k]+nsp_smooth[IJ(i,j-1)][k]+nsp_smooth[IJ(i,j+1)][k]) +
              (nsp_smooth[IJ(i-1,j-1)][k]+nsp_smooth[IJ(i+1,j-1)][k]+nsp_smooth[IJ(i-1,j+1)][k]+nsp_smooth[IJ(i+1,j+1)][k]);
          }
        }
        // and smooth edges...
        for (int i = 1; i < n1; ++i) {
          FOR_K3 xsp0_ff[IJ(i,0)][k] = (nsp_smooth[IJ(i-1,0)][k]+nsp_smooth[IJ(i+1,0)][k]);
          FOR_K3 xsp0_ff[IJ(i,n2)][k] = (nsp_smooth[IJ(i-1,n2)][k]+nsp_smooth[IJ(i+1,n2)][k]);
        }
        for (int j = 1; j < n2; ++j) {
          FOR_K3 xsp0_ff[IJ(0,j)][k] = (nsp_smooth[IJ(0,j-1)][k]+nsp_smooth[IJ(0,j+1)][k]);
          FOR_K3 xsp0_ff[IJ(n1,j)][k] = (nsp_smooth[IJ(n1,j-1)][k]+nsp_smooth[IJ(n1,j+1)][k]);
        }
      }
      delete[] nsp_smooth;
    }

    void makeLogicalQuadBl(Part * part,int disp[5],const double dn,const double dt1,const double dt2,const int n,const int mode,const int smooth_iters,double * n_plane) {
      assert(dn > 0.0);
      assert(dt1 > dn);
      assert(dt2 > dn);
      assert(kind == LOGICAL_QUAD_GROUP_KIND);
      assert(csplineVec.size() == 4);

      // recall the grid...

      const double L1 = 0.5*(csplineVec[0].getLength()+csplineVec[2].getLength());
      const double L2 = 0.5*(csplineVec[1].getLength()+csplineVec[3].getLength());
      const double h = double(n)*0.5*(dt1+dt2);
      const double delta = h;
      const double dt = delta/double(n);
      const int n1 = int(L1/dt1); assert(n1 > 0);
      const int n2 = int(L2/dt2); assert(n2 > 0);

      if (mode == 0) {

        const double sf = pow(dt/dn,1.0/double(n-1)); assert(sf > 1.0);
        SurfaceShm * surface = partVec[ipart]->surface;
        assert(surface);

        // build an adt of the group surface tris for surface projection
        vector<int> istVec;
        Adt<double> adt;
        buildTriAdtAndVecForGroup(adt,istVec,igr);
        double xws[3]; FOR_I3 xws[i] = surface->xsp[spoli_v[idata[0]]][i];
        double xes[3]; FOR_I3 xes[i] = surface->xsp[spoli_v[idata[1]]][i];
        double xen[3]; FOR_I3 xen[i] = surface->xsp[spoli_v[idata[2]]][i];
        double xwn[3]; FOR_I3 xwn[i] = surface->xsp[spoli_v[idata[3]]][i];

        // we will also need another array at each point (normals, original surface positions, etc)...
        double (*xsp_ff_tmp)[3] = new double[(n1+1)*(n2+1)][3];
        buildTransfiniteInterpPoints(xsp_ff_tmp,n1,n2,xws,xes,xen,xwn,istVec,adt);

        // now perform a series of smooth-then-reproject steps...
        if (mpi_rank == 0) cout << " > smoothing points: " << smooth_iters << " iterations..." << endl;

        // now compute the normals...
        double (*xsp0_ff)[3] = new double[(n1+1)*(n2+1)][3];
        const int smooth_iters = 10;
        if (n_plane == NULL) computeAndSmoothNormals(xsp0_ff,xsp_ff_tmp,n1,n2,smooth_iters);
        else computeAndSmoothInPlaneNormals(xsp0_ff,xsp_ff_tmp,n1,n2,smooth_iters,n_plane);

        double (*xsp_ff)[3] = part->ff_surface->xsp + disp[2]; assert(xsp_ff);

        if (mpi_rank_shared == 0) {
          // use the normals to lift xsp off the surface, and store the
          // original surface xsp in xsp0_ff...
          for (int i = 0; i <= n1; ++i) {
            for (int j = 0; j <= n2; ++j) {
              const double mag = MAG(xsp0_ff[IJ(i,j)]); assert(mag > 0.0);
              const double nij[3] = {xsp0_ff[IJ(i,j)][0]/mag,xsp0_ff[IJ(i,j)][1]/mag,xsp0_ff[IJ(i,j)][2]/mag};
              // and put the original xsp on the surface into xsp0_ff for use later...
              FOR_K3 xsp0_ff[IJ(i,j)][k] = xsp_ff_tmp[IJ(i,j)][k];
              // points ordered such that the farthest points come first...
              for (int ip = 1; ip <= n; ++ip) {
                const double dx = dn*(pow(sf,ip)-1.0)/(sf-1.0);
                const int isp = getIspFromIJK(i,j,n-ip,n1,n2,n);
                if (isp >= 0)
                  FOR_K3 xsp_ff[isp][k] = xsp0_ff[IJ(i,j)][k] - dx*nij[k];
              }
            }
          }

          /*
          // take a look...
          if (mpi_rank == 0) {
          FILE * fp = fopen("_pts0.dat","w");
          for (int i = 0; i <= n1; ++i) {
          for (int j = 0; j <= n2; ++j) {
          fprintf(fp,"%f %f %f\n",xsp0_ff[IJ(i,j)][0],xsp0_ff[IJ(i,j)][1],xsp0_ff[IJ(i,j)][2]);
          }
          }
          fclose(fp);
          fp = fopen("_pts1.dat","w");
          for (int i = 0; i <= n1; ++i) {
          for (int j = 0; j <= n2; ++j) {
          for (int k = 0; k < n; ++k) {
          //fprintf(fp,"%18.15le %18.15le %18.15le\n",xsp_ff[IJK(i,j,k)][0],xsp_ff[IJK(i,j,k)][1],xsp_ff[IJK(i,j,k)][2]);
          const int isp = getIspFromIJK(i,j,k,n1,n2,n);
          if (isp >= 0)
          fprintf(fp,"%f %f %f\n",xsp_ff[isp][0],xsp_ff[isp][1],xsp_ff[isp][2]);
          }
          }
          }
          fclose(fp);
          }
          */

        }
        else {
          for (int i = 0; i <= n1; ++i) {
            for (int j = 0; j <= n2; ++j) {
              // and put the original xsp on the surface into xsp0_ff for use later...
              FOR_K3 xsp0_ff[IJ(i,j)][k] = xsp_ff_tmp[IJ(i,j)][k];
            }
          }
        }
        delete[] xsp_ff_tmp;

        // build part->ff_surface on mpi_rank_shared == 0...
        if (mpi_rank_shared == 0) {

          int (*spost_ff)[3] = part->ff_surface->spost + disp[3]; assert(spost_ff);
          int *znost_ff = part->ff_surface->znost + disp[3]; assert(znost_ff);
          double *dxost_ff = part->ff_surface_dxost + disp[3]; assert(dxost_ff);

          // add the tris...

          int ist = 0;
          for (int i = 0; i < n1; ++i) {
            for (int j = 0; j < n2; ++j) {
              // consider the rectangle i->i+1, j->j+1...
              spost_ff[ist][0] = IJ(i,j); // we add disp[2] to all points below
              spost_ff[ist][1] = IJ(i+1,j+1);
              spost_ff[ist][2] = IJ(i+1,j);
              dxost_ff[ist]    = dt;
              znost_ff[ist]    = 0;
              ++ist;
              spost_ff[ist][0] = IJ(i,j);
              spost_ff[ist][1] = IJ(i,j+1);
              spost_ff[ist][2] = IJ(i+1,j+1);
              dxost_ff[ist]    = dt;
              znost_ff[ist]    = 0;
              ++ist;
            }
          }
          assert(ist == n1*n2*2);
          for (int k = 0; k < n-1; ++k) {
            const double dx = dn*(pow(sf,n-k)-1.0)/(sf-1.0);
            for (int j = 0; j < n2; ++j) {
              // i == 0
              spost_ff[ist][0] = getIspFromIJK(0,j,k,n1,n2,n); // we add disp[2] to all points below
              spost_ff[ist][1] = getIspFromIJK(0,j+1,k+1,n1,n2,n);
              spost_ff[ist][2] = getIspFromIJK(0,j+1,k,n1,n2,n);
              dxost_ff[ist]    = dx;
              znost_ff[ist]    = 0;
              ++ist;
              spost_ff[ist][0] = getIspFromIJK(0,j,k,n1,n2,n);
              spost_ff[ist][1] = getIspFromIJK(0,j,k+1,n1,n2,n);
              spost_ff[ist][2] = getIspFromIJK(0,j+1,k+1,n1,n2,n);
              dxost_ff[ist]    = dx;
              znost_ff[ist]    = 0;
              ++ist;
              // i == n1
              spost_ff[ist][0] = getIspFromIJK(n1,j,k,n1,n2,n); // we add disp[2] to all points below
              spost_ff[ist][1] = getIspFromIJK(n1,j+1,k,n1,n2,n);
              spost_ff[ist][2] = getIspFromIJK(n1,j+1,k+1,n1,n2,n);
              dxost_ff[ist]    = dx;
              znost_ff[ist]    = 0;
              ++ist;
              spost_ff[ist][0] = getIspFromIJK(n1,j,k,n1,n2,n);
              spost_ff[ist][1] = getIspFromIJK(n1,j+1,k+1,n1,n2,n);
              spost_ff[ist][2] = getIspFromIJK(n1,j,k+1,n1,n2,n);
              dxost_ff[ist]    = dx;
              znost_ff[ist]    = 0;
              ++ist;
            }
            for (int i = 0; i < n1; ++i) {
              // j == 0
              spost_ff[ist][0] = getIspFromIJK(i,0,k,n1,n2,n); // we add disp[2] to all points below
              spost_ff[ist][1] = getIspFromIJK(i+1,0,k+1,n1,n2,n);
              spost_ff[ist][2] = getIspFromIJK(i,0,k+1,n1,n2,n);
              dxost_ff[ist]    = dx;
              znost_ff[ist]    = 0;
              ++ist;
              spost_ff[ist][0] = getIspFromIJK(i,0,k,n1,n2,n);
              spost_ff[ist][1] = getIspFromIJK(i+1,0,k,n1,n2,n);
              spost_ff[ist][2] = getIspFromIJK(i+1,0,k+1,n1,n2,n);
              dxost_ff[ist]    = dx;
              znost_ff[ist]    = 0;
              ++ist;
              // j == n2
              spost_ff[ist][0] = getIspFromIJK(i,n2,k,n1,n2,n); // we add disp[2] to all points below
              spost_ff[ist][1] = getIspFromIJK(i,n2,k+1,n1,n2,n);
              spost_ff[ist][2] = getIspFromIJK(i+1,n2,k+1,n1,n2,n);
              dxost_ff[ist]    = dx;
              znost_ff[ist]    = 0;
              ++ist;
              spost_ff[ist][0] = getIspFromIJK(i,n2,k,n1,n2,n);
              spost_ff[ist][1] = getIspFromIJK(i+1,n2,k+1,n1,n2,n);
              spost_ff[ist][2] = getIspFromIJK(i+1,n2,k,n1,n2,n);
              dxost_ff[ist]    = dx;
              znost_ff[ist]    = 0;
              ++ist;
            }
          }
          assert(ist == (4*(n-1)*(n1+n2) + 2*n1*n2));

          // now add the perimiter points from the group to the end of the ff_surface->xsp...
          // spline0...
          int isp = (n1+1)*(n2+1)+2*(n-1)*(n1+n2);
          vector<int> line0;
          for (int ii = 0; ii < csplineVec[0].getNp(); ++ii) {
            csplineVec[0].getXp(xsp_ff[isp],ii);
            line0.push_back(isp);
            ++isp;
          }
          // this should correspond to the i-line along j=0...
          vector<int> line1;
          for (int i = 0; i <= n1; ++i)
            line1.push_back(getIspFromIJK(i,0,n-1,n1,n2,n));
          ist += GeomUtils::facetGap(spost_ff+ist,line0,line1,xsp_ff,false);
          // spline1...
          line0.clear();
          for (int ii = 0; ii < csplineVec[1].getNp(); ++ii) {
            csplineVec[1].getXp(xsp_ff[isp],ii);
            line0.push_back(isp);
            ++isp;
          }
          line1.clear();
          for (int j = 0; j <= n2; ++j)
            line1.push_back(getIspFromIJK(n1,j,n-1,n1,n2,n));
          ist += GeomUtils::facetGap(spost_ff+ist,line0,line1,xsp_ff,false);
          // spline2...
          line0.clear();
          for (int ii = 0; ii < csplineVec[2].getNp(); ++ii) {
            csplineVec[2].getXp(xsp_ff[isp],ii);
            line0.push_back(isp);
            ++isp;
          }
          line1.clear();
          for (int i = n1; i >= 0; --i)
            line1.push_back(getIspFromIJK(i,n2,n-1,n1,n2,n));
          ist += GeomUtils::facetGap(spost_ff+ist,line0,line1,xsp_ff,false);
          // spline3...
          line0.clear();
          for (int ii = 0; ii < csplineVec[3].getNp(); ++ii) {
            csplineVec[3].getXp(xsp_ff[isp],ii);
            line0.push_back(isp);
            ++isp;
          }
          line1.clear();
          for (int j = n2; j >= 0; --j)
            line1.push_back(getIspFromIJK(0,j,n-1,n1,n2,n));
          ist += GeomUtils::facetGap(spost_ff+ist,line0,line1,xsp_ff,false);

          // offset the spost_ff for all pts...
          for (int ist_ = 0; ist_ < ist; ++ist_)
            FOR_I3 spost_ff[ist_][i] += disp[2];

          // set the znost and dxost for remaining tris...
          for (int ist_ = (4*(n-1)*(n1+n2) + 2*n1*n2); ist_ < ist; ++ist_) {
            znost_ff[ist_] = 0;
            dxost_ff[ist_] = dn;
          }

          // this is only valid if there is one group...
          //assert(isp == part->ff_surface->nsp);
          //assert(ist == part->ff_surface->nst);
          disp[2] += isp;
          disp[3] += ist;

          // take a look...
          //if (mpi_rank == 0) GeomUtils::writeSbin("ff.sbin",spost_ff,ist,xsp_ff);
        }

        MPI_Barrier(mpi_comm_shared);

        // the new points go in part->pts...
        assert(part->pts);
        double (*xp)[3] = part->pts->xp + disp[4];
        double *deltap = part->pts->delta + disp[4];

        int ip = 0;
        int rank_flag = 0; // for point round-robin
        vector<int> intVec;
        for (int i_ = 0; i_ < n1; ++i_) {
          for (int j_ = 0; j_ < n2; ++j_) {
            if (rank_flag == mpi_rank) {
              // grab the points at the bottom and the top of the strand. These
              // are the simple average of the 4 surrounding points...
              double xp0[3]; FOR_K3 xp0[k] = 0.25*(xsp0_ff[IJ(i_,j_)][k]+xsp0_ff[IJ(i_+1,j_)][k]+xsp0_ff[IJ(i_,j_+1)][k]+xsp0_ff[IJ(i_+1,j_+1)][k]);
              // and project this point to the surface...
              assert(intVec.empty());
              adt.buildListForClosestPoint(intVec,xp0);
              assert(!intVec.empty());
              double d2_min = HUGE_VAL;
              double xp_min[3];
              for (int ii = 0; ii < intVec.size(); ++ii) {
                const int ist = istVec[intVec[ii]];
                double xp_[3]; MiscUtils::getClosestPointOnTriRobust(xp_,xp0,
                                                                     surface->xsp[surface->spost[ist][0]],
                                                                     surface->xsp[surface->spost[ist][1]],
                                                                     surface->xsp[surface->spost[ist][2]]);
                const double d2 = DIST2(xp_,xp0);
                if (d2 < d2_min) {
                  d2_min = d2;
                  FOR_K3 xp_min[k] = xp_[k];
                }
              }
              intVec.clear();
              assert(d2_min < HUGE_VAL);
              FOR_K3 xp0[k] = xp_min[k];
              // now use the point at the top to build a unit-normal direction...
              double np0[3]; FOR_K3 np0[k] = 0.25*(xsp_ff[IJ(i_,j_)][k]+xsp_ff[IJ(i_+1,j_)][k]+xsp_ff[IJ(i_,j_+1)][k]+xsp_ff[IJ(i_+1,j_+1)][k]) - xp0[k];
              const double mag = MAG(np0);
              FOR_K3 np0[k] /= mag;
              for (int ip_ = 0; ip_ < n; ++ip_) {
                const double dx = dn*0.5*( (pow(sf,ip_)-1.0)/(sf-1.0) + (pow(sf,ip_+1)-1.0)/(sf-1.0) );
                FOR_K3 xp[ip][k] = xp0[k] + dx*np0[k];
                deltap[ip] = 1.2*sqrt(3*dt*dt); // use largest scale
                ++ip;
              }
            }
            ++rank_flag;
            if (rank_flag == mpi_size)
              rank_flag = 0;
          }
        }

        // store the ip count in disp[4]...
        disp[4] += ip;

        delete[] xsp0_ff;
      }
      else if (mode == 1) {

        // build part->ff_surface on mpi_rank_shared == 0...
        if (mpi_rank_shared == 0) {

          SurfaceShm * surface = partVec[ipart]->surface;
          assert(surface);

          // put the double[(n1+1)*(n2+1)][3] points at the start of the ff_surface->xsp,
          // which has already been allocated...
          // note we offset this by the disp[2] to get a zero-base ptr (because
          // the calling process may be putting multiple objects into this
          // ff_surface)...
          double (*xsp_ff)[3] = part->ff_surface->xsp + disp[2]; assert(xsp_ff);


          // build an adt of the group surface tris for surface projection
          vector<int> istVec;
          Adt<double> adt;
          buildTriAdtAndVecForGroup(adt,istVec,igr);
          double xws[3]; FOR_I3 xws[i] = surface->xsp[spoli_v[idata[0]]][i];
          double xes[3]; FOR_I3 xes[i] = surface->xsp[spoli_v[idata[1]]][i];
          double xen[3]; FOR_I3 xen[i] = surface->xsp[spoli_v[idata[2]]][i];
          double xwn[3]; FOR_I3 xwn[i] = surface->xsp[spoli_v[idata[3]]][i];

          buildTransfiniteInterpPoints(xsp_ff,n1,n2,xws,xes,xen,xwn,istVec,adt);

          // now perform a series of smooth-then-reproject steps...
          if (mpi_rank == 0) cout << " > smoothing points 50 iterations..." << endl;
          double d2_max = smoothAndReprojectInterpPoints(xsp_ff,smooth_iters,n1,n2,istVec,adt);

          // now compute the normals...
          double (*xsp0_ff)[3] = new double[(n1+1)*(n2+1)][3];
          const int smooth_iters = 10;
          computeAndSmoothNormals(xsp0_ff,xsp_ff,n1,n2,smooth_iters);

          // use the normals to lift xsp off the surface, and store the
          // original surface xsp in xsp0_ff...
          for (int i = 0; i <= n1; ++i) {
            for (int j = 0; j <= n2; ++j) {
              const double mag = MAG(xsp0_ff[IJ(i,j)]); assert(mag > 0.0);
              FOR_K3 {
                xsp0_ff[IJ(i,j)][k] /= mag;
                const double tmp = xsp_ff[IJ(i,j)][k];
                xsp_ff[IJ(i,j)][k] -= h*xsp0_ff[IJ(i,j)][k];
                // and put the original xsp on the surface into xsp0_ff for use later...
                xsp0_ff[IJ(i,j)][k] = tmp;
              }
            }
          }
          /*
          // take a look...
          if (mpi_rank == 0) {
	  FILE * fp = fopen("pts0.dat","w");
	  for (int i = 0; i <= n1; ++i) {
	  for (int j = 0; j <= n2; ++j) {
	  fprintf(fp,"%f %f %f\n",xsp0_ff[IJ(i,j)][0],xsp0_ff[IJ(i,j)][1],xsp0_ff[IJ(i,j)][2]);
	  }
	  }
	  fclose(fp);
	  fp = fopen("pts1.dat","w");
	  for (int i = 0; i <= n1; ++i) {
	  for (int j = 0; j <= n2; ++j) {
	  fprintf(fp,"%f %f %f\n",xsp_ff[IJ(i,j)][0],xsp_ff[IJ(i,j)][1],xsp_ff[IJ(i,j)][2]);
	  }
	  }
	  fclose(fp);
          }
          */
          // add the tris...

          int (*spost_ff)[3] = part->ff_surface->spost + disp[3]; assert(spost_ff);
          int *znost_ff = part->ff_surface->znost + disp[3]; assert(znost_ff);

          int ist = 0;
          for (int i = 0; i < n1; ++i) {
            for (int j = 0; j < n2; ++j) {
              // consider the rectangle i->i+1, j->j+1...
              spost_ff[ist][0] = IJ(i,j); // we add disp[2] to all points below
              spost_ff[ist][1] = IJ(i+1,j+1);
              spost_ff[ist][2] = IJ(i+1,j);
              znost_ff[ist]    = 0;
              ++ist;
              spost_ff[ist][0] = IJ(i,j);
              spost_ff[ist][1] = IJ(i,j+1);
              spost_ff[ist][2] = IJ(i+1,j+1);
              znost_ff[ist]    = 0;
              ++ist;
            }
          }
          assert(ist == n1*n2*2);
          // now add the perimiter points from the group to the end of the ff_surface->xsp...
          // spline0...
          int isp = (n1+1)*(n2+1);
          vector<int> line0;
          for (int ii = 0; ii < csplineVec[0].getNp(); ++ii) {
            csplineVec[0].getXp(xsp_ff[isp],ii);
            line0.push_back(isp);
            ++isp;
          }
          // this should correspond to the i-line along j=0...
          vector<int> line1;
          for (int i = 0; i <= n1; ++i)
            line1.push_back(IJ(i,0));
          ist += GeomUtils::facetGap(spost_ff+ist,line0,line1,xsp_ff,false);
          // spline1...
          line0.clear();
          for (int ii = 0; ii < csplineVec[1].getNp(); ++ii) {
            csplineVec[1].getXp(xsp_ff[isp],ii);
            line0.push_back(isp);
            ++isp;
          }
          line1.clear();
          for (int j = 0; j <= n2; ++j)
            line1.push_back(IJ(n1,j));
          ist += GeomUtils::facetGap(spost_ff+ist,line0,line1,xsp_ff,false);
          // spline2...
          line0.clear();
          for (int ii = 0; ii < csplineVec[2].getNp(); ++ii) {
            csplineVec[2].getXp(xsp_ff[isp],ii);
            line0.push_back(isp);
            ++isp;
          }
          line1.clear();
          for (int i = n1; i >= 0; --i)
            line1.push_back(IJ(i,n2));
          ist += GeomUtils::facetGap(spost_ff+ist,line0,line1,xsp_ff,false);
          // spline3...
          line0.clear();
          for (int ii = 0; ii < csplineVec[3].getNp(); ++ii) {
            csplineVec[3].getXp(xsp_ff[isp],ii);
            line0.push_back(isp);
            ++isp;
          }
          line1.clear();
          for (int j = n2; j >= 0; --j)
            line1.push_back(IJ(0,j));
          ist += GeomUtils::facetGap(spost_ff+ist,line0,line1,xsp_ff,false);

          // offset the spost_ff and set the dxost for all pts...
          double *dxost_ff = part->ff_surface_dxost + disp[3]; assert(dxost_ff);
          for (int ist_ = 0; ist_ < ist; ++ist_) {
            FOR_I3 spost_ff[ist_][i] += disp[2];
            dxost_ff[ist_] = dt;
          }

          // this is only valid if there is one group...
          //assert(isp == part->ff_surface->nsp);
          //assert(ist == part->ff_surface->nst);
          disp[2] += isp;
          disp[3] += ist;

          // take a look...

          /*
            if (mpi_rank == 0) GeomUtils::writeSbin("ff.sbin",spost,ist,xsp);
            if (mpi_rank == 0) GeomUtils::writeTecplot("ff.dat",spost,ist,xsp);
            {
            FILE * fp = fopen("line0.dat","w");
            for (int ii = 0; ii < line0.size(); ++ii) {
            const int isp_ = line0[ii];
            fprintf(fp,"%18.15le %18.15le %18.15le\n",xsp[isp_][0],xsp[isp_][1],xsp[isp_][2]);
            }
            fclose(fp);
            fp = fopen("line1.dat","w");
            for (int ii = 0; ii < line1.size(); ++ii) {
            const int isp_ = line1[ii];
            fprintf(fp,"%18.15le %18.15le %18.15le\n",xsp[isp_][0],xsp[isp_][1],xsp[isp_][2]);
            }
            fclose(fp);
            }
          */

          // the new points go in part->pts...
          assert(part->pts);
          double (*xp)[3] = part->pts->xp + disp[4];
          double *deltap = part->pts->delta + disp[4];

          int ip = 0;
          double dn_current = dt;
          int np = n;
          double sf = 1.0;
          int i = 0, j = 0;
          int rank_flag = 0; // for point round-robin: note rank 0 per node only...
          vector<int> intVec;
          while (1) {
            // add the perimiter times the current np to the points count...
            const int i_f = i;
            const int i_l = n1-i-1;
            assert(i_l >= i_f);
            const int j_f = j;
            const int j_l = n2-j-1;
            assert(j_l >= j_f);
            if ((mpi_rank_shared == 0)&&(mpi_rank_internode == rank_flag)) {
              for (int i_ = i_f; i_ <= i_l; ++i_) {
                for (int j_ = j_f; j_ <= j_l; ++j_) {
                  // perimiter only...
                  if ((i_ == i_f)||(i_ == i_l)||(j_ == j_f)||(j_ == j_l)) {
                    // grab the points at the bottom and the top of the strand. These
                    // are the simple average of the 4 surrounding points...
                    double xp0[3]; FOR_K3 xp0[k] = 0.25*(xsp0_ff[IJ(i_,j_)][k]+xsp0_ff[IJ(i_+1,j_)][k]+xsp0_ff[IJ(i_,j_+1)][k]+xsp0_ff[IJ(i_+1,j_+1)][k]);
                    // and project this point to the surface...
                    assert(intVec.empty());
                    adt.buildListForClosestPoint(intVec,xp0);
                    assert(!intVec.empty());
                    double d2_min = HUGE_VAL;
                    double xp_min[3];
                    for (int ii = 0; ii < intVec.size(); ++ii) {
                      const int ist = istVec[intVec[ii]];
                      double xp_[3]; MiscUtils::getClosestPointOnTriRobust(xp_,xp0,
                                                                           surface->xsp[surface->spost[ist][0]],
                                                                           surface->xsp[surface->spost[ist][1]],
                                                                           surface->xsp[surface->spost[ist][2]]);
                      const double d2 = DIST2(xp_,xp0);
                      if (d2 < d2_min) {
                        d2_min = d2;
                        FOR_K3 xp_min[k] = xp_[k];
                      }
                    }
                    intVec.clear();
                    assert(d2_min < HUGE_VAL);
                    d2_max = max(d2_max,DIST2(xsp_ff[IJ(i,j)],xp_min));
                    FOR_K3 xp0[k] = xp_min[k];
                    // now use the point at the top to build a unit-normal direction...
                    double np0[3]; FOR_K3 np0[k] = 0.25*(xsp_ff[IJ(i_,j_)][k]+xsp_ff[IJ(i_+1,j_)][k]+xsp_ff[IJ(i_,j_+1)][k]+xsp_ff[IJ(i_+1,j_+1)][k]) - xp0[k];
                    const double mag = MAG(np0);
                    FOR_K3 np0[k] /= mag;
                    if (np == n) {
                      // this is the uniform section...
                      assert(sf == 1.0);
                      assert(dn_current == dt);
                      for (int ip_ = 0; ip_ < np; ++ip_) {
                        const double dx = delta*(double(ip_)+0.5)/double(np);
                        FOR_K3 xp[ip][k] = xp0[k] + dx*np0[k];
                        deltap[ip] = 1.2*sqrt(3*dt*dt);
                        ++ip;
                      }
                    }
                    else {
                      assert(sf > 1.0);
                      for (int ip_ = 0; ip_ < np; ++ip_) {
                        const double dx = dn_current*0.5*( (pow(sf,ip_)-1.0)/(sf-1.0) + (pow(sf,ip_+1)-1.0)/(sf-1.0) );
                        FOR_K3 xp[ip][k] = xp0[k] + dx*np0[k];
                        deltap[ip] = 1.2*sqrt(3*dt*dt);
                        ++ip;
                      }
                    }
                  }
                }
              }
            }
            if ((i_f == i_l)||(i_f+1 == i_l)||(j_f == j_l)||(j_f+1 == j_l)) {
              // that was the last row/column...
              break;
            }
            // move in one cell in all directions...
            ++i;
            ++j;
            // modify rank_flag for next region...
            ++rank_flag;
            if (rank_flag == mpi_size_internode)
              rank_flag = 0;
            // add another point if we have not reached the target dn...
            if (dn_current > dn*1.01) {
              // add another point to the strand...
              ++np;
              // and recompute dn_current...
              // this starts us on the correct side of the function, half-way towards the
              // zero gradient...
              dn_current = 0.5*exp(log(dt*double(np-1)/(delta-dt))/double(np))*(delta-dt)/double(np-1);
              int iter = 0;
              while (1) {
                ++iter;
                if (iter > 30) {
                  if (mpi_rank == 0) cout << "ERROR: cannot converge dn/delta system: delta,dt,np: " << delta << " " << dt << " " << np << endl;
                  assert(0);
                }
                const double t1 = pow(dt/dn_current,double(np)/double(np-1));
                const double t2 = pow(dt/dn_current,1.0/double(np-1));
                const double f = dn_current*t1-delta*t2+delta-dn_current;
                const double fp = (delta*t2-dn_current*t1-dn_current*double(np-1))/(dn_current*double(np-1));
                //if (mpi_rank == 0) cout << "XXX np = " << np << " got dn_current: " << dn_current << " f: " << f << " fp: " << fp << " f/fp/dt: " << f/fp/dt << endl;
                if (f > 0.5*dn_current*fp)
                  dn_current *= 0.5; // revert to bisection for stability
                else {
                  dn_current -= f/fp;
                  if (fabs(f/fp) < 1.0E-12*dt)
                    break;
                }
              }
              // recompute sf...
              sf = pow(dt/dn_current,1.0/double(np-1));
              //const double delta_check = (dt*sf-dn_current)/(sf-1.0);
              //if (mpi_rank == 0) cout << " XXXX got sf: " << sf << " dn_current: " << dn_current << " delta_check: " << delta_check << " delta: " << delta << endl;
            }
          }
          // store the ip count in disp[4]...
          disp[4] += ip;

          delete[] xsp0_ff;
        }

        // only valid if there is just one group...
        //assert(disp[4] == part->pts->np);
        //GeomUtils::writePtsTecplot("pts.dat",part->pts->xp,part->pts->np);
        //MPI_Pause("SO FAR SO GOOD");

      }
      else {

        assert(0);

      }

    }

#undef IJ

    bool isRectangle(int sols[4],const double tol) {
      cout << "blah" << endl;
      return false;
    }

    void countRectangleBl(int counts[5],const int sols_[4],const double dn,const double dt,const int n,const int mode) {
      assert(dt > dn);
      assert(nlink == 1);
      int sols[4];
      if (mode == 1) {
        FOR_I3 sols[i] = sols_[i+1];
        sols[3] = sols_[0];
      }
      else {
        assert((mode == 0)||(mode == 2));
        FOR_I4 sols[i] = sols_[i];
      }
      SurfaceShm * surface = partVec[ipart]->surface;
      const double * const x0 = surface->xsp[spoli_v[sols[0]]];
      const double u[3] = DIFF(x0,surface->xsp[spoli_v[sols[1]]]);
      const double v[3] = DIFF(x0,surface->xsp[spoli_v[sols[3]]]);
      const double L = MAG(v);
      const double W = MAG(u);
      const double sf = pow(dt/dn,1.0/double(n-1));
      const double delta = (dt*sf-dn)/(sf-1.0);
      if (mpi_rank == 0) {
        cout << " > L: " << L << " W: " << W << endl;
        cout << " > BL sf: " << sf << " delta for " << n << " layers: " << delta << endl;
      }
      // no surface counts...
      counts[0] = counts[1] = 0;
      int rank_flag = 0;
      if (mode == 2) {
        // ff_surface point count...
        counts[2] = spoli_s + 4*n;
        // ff_surface tri count...
        counts[3] = spoli_s + 4 + 2 + 8*(n-1);
        const int nW = int(W/dt);
        const int nL = int(L/dt);
        // pts->np
        counts[4] = 0;
        if (mpi_rank == 0) cout << " > total point estimate: nL: " << nL << " x nW: " << nW << " x n: " << n << " = " << nL*nW*n << endl;
        for (int i = 0; i < nL; ++i) {
          for (int k = 0; k < nW; ++k) {
            for (int j = 0; j < n; ++j) {
              if (rank_flag == mpi_rank)
                ++counts[4];
            }
          }
          ++rank_flag;
          if (rank_flag == mpi_size)
            rank_flag = 0;
        }
      }
      else {
        assert((mode == 0)||(mode == 1));
        assert(L > 2.0*delta);
        // ff_surface point count...
        counts[2] = spoli_s + 4 + 8*(n-1);
        // ff_surface tri count...
        counts[3] = spoli_s + 4 + 2 + 16*(n-1);
        const int nW = int(W/dt);
        const int nL = int((L-2.0*delta)/dt); // + 2*n
        if (mpi_rank == 0) cout << " > total point estimate: nL+2*n: " << nL+2*n << " x nW: " << nW << " x n: " << n << " = " << (nL+2*n)*nW*n << endl;
        // pts->np
        counts[4] = 0;
        for (int i = 0; i < n; ++i) {
          const double dx = dn*pow(sf,i);
          const double this_nW = int(W/dx);
          //cout << i << " " << this_nW << " " << nW << endl;
          //if (i == n-1) assert(this_nW == nW); // i think this is only true for certain dt/dn's (integer?)
          for (int k = 0; k < this_nW; ++k) {
            for (int j = 0; j <= i; ++j) {
              if (rank_flag == mpi_rank)
                ++counts[4];
            }
          }
          ++rank_flag;
          if (rank_flag == mpi_size)
            rank_flag = 0;
        }
        for (int i = 0; i < nL; ++i) {
          for (int k = 0; k < nW; ++k) {
            for (int j = 0; j < n; ++j) {
              if (rank_flag == mpi_rank)
                ++counts[4];
            }
          }
          ++rank_flag;
          if (rank_flag == mpi_size)
            rank_flag = 0;
        }
        for (int i = n-1; i >= 0; --i) {
          const double dx = dn*pow(sf,i);
          const double this_nW = int(W/dx);
          //if (i == n-1) assert(this_nW == nW);
          for (int k = 0; k < this_nW; ++k) {
            for (int j = 0; j <= i; ++j) {
              if (rank_flag == mpi_rank)
                ++counts[4];
            }
          }
          ++rank_flag;
          if (rank_flag == mpi_size)
            rank_flag = 0;
        }
      }
    }
    void makeRectangleBl(Part * part,int disp[5],const int sols_[4],const double dn,const double dt,const int n,const int mode) {
      assert(dt > dn);
      assert(nlink == 1);
      int sols[4];
      if (mode == 1) {
        FOR_I3 sols[i] = sols_[i+1];
        sols[3] = sols_[0];
      }
      else {
        assert((mode == 0)||(mode == 2));
        FOR_I4 sols[i] = sols_[i];
      }
      SurfaceShm * surface = partVec[ipart]->surface;
      const double * const x0 = surface->xsp[spoli_v[sols[0]]];
      const double * const x1 = surface->xsp[spoli_v[sols[1]]];
      const double u[3] = DIFF(x1,x0);
      const double v[3] = DIFF(surface->xsp[spoli_v[sols[3]]],x0);
      const double L = MAG(v);
      const double W = MAG(u);
      const double w[3] = CROSS_PRODUCT(u,v);
      double e0[3],e1[3],e2[3];
      FOR_I3 e0[i] = u[i]/W;
      FOR_I3 e1[i] = v[i]/L;
      FOR_I3 e2[i] = -w[i]/area; // unitize using area
      if (mpi_rank == 0) cout << " > e0: " << COUT_VEC(e0) << " e1: " << COUT_VEC(e1) << " e2: " << COUT_VEC(e2) << endl;
      const double sf = pow(dt/dn,1.0/double(n-1));
      const double delta = (dt*sf-dn)/(sf-1.0);
      assert(part->ff_surface);
      if (mode == 2) {
        const int nL = int(L/dt);
        const int nW = int(W/dt);
        if (mpi_rank_shared == 0) {
          // ===================================================
          // surface points: xsp_ff...
          // ===================================================
          const int nsp_ff = spoli_s + 4*n;
          // copy the links into xsp_ff first...
          for (int sol = 0; sol < spoli_s; ++sol) {
            const int isp = spoli_v[sol];
            FOR_I3 part->ff_surface->xsp[disp[2]+sol][i] = partVec[ipart]->surface->xsp[isp][i];
          }
          // then the layers start at spoli_i[1]...
          int isp_ff = spoli_s;
          for (int k = 0; k < n; ++k) {
            const double dl = dn*(pow(sf,k+1)-1.0)/(sf-1.0);
            FOR_I3 part->ff_surface->xsp[disp[2]+isp_ff][i] = x0[i] +           dl*e2[i]; ++isp_ff;
            FOR_I3 part->ff_surface->xsp[disp[2]+isp_ff][i] = x1[i] +           dl*e2[i]; ++isp_ff;
            FOR_I3 part->ff_surface->xsp[disp[2]+isp_ff][i] = x1[i] + L*e1[i] + dl*e2[i]; ++isp_ff;
            FOR_I3 part->ff_surface->xsp[disp[2]+isp_ff][i] = x0[i] + L*e1[i] + dl*e2[i]; ++isp_ff;
          }
          assert(isp_ff == nsp_ff);
          /*
            if (mpi_rank == 0) {
            for (int isp_ff = 0; isp_ff < nsp_ff; ++isp_ff) {
            cout << part->ff_surface->xsp[disp[2]+isp_ff][0] << ","
            << part->ff_surface->xsp[disp[2]+isp_ff][1] << ","
            << part->ff_surface->xsp[disp[2]+isp_ff][2] << endl;
            }
            }
          */

          // ===================================================
          // surface tris: xst_ff...
          // ===================================================
          const int nst_ff = spoli_s + 4 + 2 + 8*(n-1);
          // 1. the lip tris are special because they have to respect both discretizations...
          int ist_ff = 0;

          // tris from first layer...

          part->ff_surface->spost[ist_ff+disp[3]][0] = disp[2]+sols[0];
          part->ff_surface->spost[ist_ff+disp[3]][1] = disp[2]+spoli_s+4;
          part->ff_surface->spost[ist_ff+disp[3]][2] = disp[2]+spoli_s+0;
          part->ff_surface_dxost[ist_ff+disp[3]] = dn;
          ++ist_ff;
          part->ff_surface->spost[ist_ff+disp[3]][0] = disp[2]+sols[1];
          part->ff_surface->spost[ist_ff+disp[3]][1] = disp[2]+spoli_s+0;
          part->ff_surface->spost[ist_ff+disp[3]][2] = disp[2]+spoli_s+1;
          part->ff_surface_dxost[ist_ff+disp[3]] = dn;
          ++ist_ff;
          part->ff_surface->spost[ist_ff+disp[3]][0] = disp[2]+sols[2];
          part->ff_surface->spost[ist_ff+disp[3]][1] = disp[2]+spoli_s+1;
          part->ff_surface->spost[ist_ff+disp[3]][2] = disp[2]+spoli_s+2;
          part->ff_surface_dxost[ist_ff+disp[3]] = dn;
          ++ist_ff;
          part->ff_surface->spost[ist_ff+disp[3]][0] = disp[2]+sols[3];
          part->ff_surface->spost[ist_ff+disp[3]][1] = disp[2]+spoli_s+2;
          part->ff_surface->spost[ist_ff+disp[3]][2] = disp[2]+spoli_s+3;
          part->ff_surface_dxost[ist_ff+disp[3]] = dn;
          ++ist_ff;

          // tris from lip...

          int isp_ff_prev = sols[0];
          for (int isp_ff = isp_ff_prev+1; isp_ff <= sols[1]; ++isp_ff) {
            part->ff_surface->spost[ist_ff+disp[3]][0] = disp[2]+isp_ff_prev;
            part->ff_surface->spost[ist_ff+disp[3]][1] = disp[2]+spoli_s;
            part->ff_surface->spost[ist_ff+disp[3]][2] = disp[2]+isp_ff;
            part->ff_surface_dxost[ist_ff+disp[3]] = dn;
            ++ist_ff;
            isp_ff_prev = isp_ff;
          }
          isp_ff_prev = sols[1];
          for (int isp_ff = isp_ff_prev+1; isp_ff <= sols[2]; ++isp_ff) {
            part->ff_surface->spost[ist_ff+disp[3]][0] = disp[2]+isp_ff_prev;
            part->ff_surface->spost[ist_ff+disp[3]][1] = disp[2]+spoli_s+1;
            part->ff_surface->spost[ist_ff+disp[3]][2] = disp[2]+isp_ff;
            part->ff_surface_dxost[ist_ff+disp[3]] = dn;
            ++ist_ff;
            isp_ff_prev = isp_ff;
          }
          isp_ff_prev = sols[2];
          for (int isp_ff = isp_ff_prev+1; isp_ff <= sols[3]; ++isp_ff) {
            part->ff_surface->spost[ist_ff+disp[3]][0] = disp[2]+isp_ff_prev;
            part->ff_surface->spost[ist_ff+disp[3]][1] = disp[2]+spoli_s+2;
            part->ff_surface->spost[ist_ff+disp[3]][2] = disp[2]+isp_ff;
            part->ff_surface_dxost[ist_ff+disp[3]] = dn;
            ++ist_ff;
            isp_ff_prev = isp_ff;
          }
          isp_ff_prev = sols[3];
          for (int isp_ff = isp_ff_prev+1; isp_ff < spoli_s; ++isp_ff) {
            part->ff_surface->spost[ist_ff+disp[3]][0] = disp[2]+isp_ff_prev;
            part->ff_surface->spost[ist_ff+disp[3]][1] = disp[2]+spoli_s+3;
            part->ff_surface->spost[ist_ff+disp[3]][2] = disp[2]+isp_ff;
            part->ff_surface_dxost[ist_ff+disp[3]] = dn;
            ++ist_ff;
            isp_ff_prev = isp_ff;
          }
          isp_ff_prev = 0;
          for (int isp_ff = isp_ff_prev+1; isp_ff <= sols[0]; ++isp_ff) {
            part->ff_surface->spost[ist_ff+disp[3]][0] = disp[2]+isp_ff_prev;
            part->ff_surface->spost[ist_ff+disp[3]][1] = disp[2]+spoli_s+3;
            part->ff_surface->spost[ist_ff+disp[3]][2] = disp[2]+isp_ff;
            part->ff_surface_dxost[ist_ff+disp[3]] = dn;
            ++ist_ff;
            isp_ff_prev = isp_ff;
          }
          part->ff_surface->spost[ist_ff+disp[3]][0] = disp[2]+spoli_s-1;
          part->ff_surface->spost[ist_ff+disp[3]][1] = disp[2]+spoli_s+3;
          part->ff_surface->spost[ist_ff+disp[3]][2] = disp[2];
          part->ff_surface_dxost[ist_ff+disp[3]] = dn;
          ++ist_ff;

          for (int k = 0; k < n-1; ++k) {
            int isp_ff = spoli_s+4*k;
            const double dl = dn*pow(sf,k+1);
            FOR_I3 {
              part->ff_surface->spost[ist_ff+disp[3]][0] = disp[2]+isp_ff+i;
              part->ff_surface->spost[ist_ff+disp[3]][1] = disp[2]+isp_ff+i+4;
              part->ff_surface->spost[ist_ff+disp[3]][2] = disp[2]+isp_ff+i+5;
              part->ff_surface_dxost[ist_ff+disp[3]] = dl;
              ++ist_ff;
              part->ff_surface->spost[ist_ff+disp[3]][0] = disp[2]+isp_ff+i;
              part->ff_surface->spost[ist_ff+disp[3]][1] = disp[2]+isp_ff+i+5;
              part->ff_surface->spost[ist_ff+disp[3]][2] = disp[2]+isp_ff+i+1;
              part->ff_surface_dxost[ist_ff+disp[3]] = dl;
              ++ist_ff;
            }
            part->ff_surface->spost[ist_ff+disp[3]][0] = disp[2]+isp_ff+3;
            part->ff_surface->spost[ist_ff+disp[3]][1] = disp[2]+isp_ff+7;
            part->ff_surface->spost[ist_ff+disp[3]][2] = disp[2]+isp_ff+4;
            part->ff_surface_dxost[ist_ff+disp[3]] = dl;
            ++ist_ff;
            part->ff_surface->spost[ist_ff+disp[3]][0] = disp[2]+isp_ff+3;
            part->ff_surface->spost[ist_ff+disp[3]][1] = disp[2]+isp_ff+4;
            part->ff_surface->spost[ist_ff+disp[3]][2] = disp[2]+isp_ff+0;
            part->ff_surface_dxost[ist_ff+disp[3]] = dl;
            ++ist_ff;
            if (k == n-2) {
              // top
              part->ff_surface->spost[ist_ff+disp[3]][0] = disp[2]+isp_ff+4;
              part->ff_surface->spost[ist_ff+disp[3]][1] = disp[2]+isp_ff+6;
              part->ff_surface->spost[ist_ff+disp[3]][2] = disp[2]+isp_ff+5;
              part->ff_surface_dxost[ist_ff+disp[3]] = dl;
              ++ist_ff;
              part->ff_surface->spost[ist_ff+disp[3]][0] = disp[2]+isp_ff+4;
              part->ff_surface->spost[ist_ff+disp[3]][1] = disp[2]+isp_ff+7;
              part->ff_surface->spost[ist_ff+disp[3]][2] = disp[2]+isp_ff+6;
              part->ff_surface_dxost[ist_ff+disp[3]] = dl;
              ++ist_ff;
            }
          }
          //if (mpi_rank == 0) part->ff_surface->writeTecplotNoZones("ff_surface.plt");
          disp[2] += nsp_ff;
          assert(ist_ff == nst_ff);
          disp[3] += nst_ff;
        }
        // ==================================================================
        // now the points...
        // ==================================================================
        assert(part->pts);
        int rank_flag = 0;
        for (int i = 0; i < nL; ++i) {
          const double yp = (double(i)+0.5)/double(nL)*L;
          for (int k = 0; k < nW; ++k) {
            const double xp = (double(k)+0.5)/double(nW)*W;
            for (int j = 0; j < n; ++j) {
              const double zp = dn*(pow(sf,0.5+double(j))-1.0)/(sf-1.0);
              if (rank_flag == mpi_rank) {
                part->pts->xp[disp[4]][0] = x0[0] + xp*e0[0] + yp*e1[0] + zp*e2[0];
                part->pts->xp[disp[4]][1] = x0[1] + xp*e0[1] + yp*e1[1] + zp*e2[1];
                part->pts->xp[disp[4]][2] = x0[2] + xp*e0[2] + yp*e1[2] + zp*e2[2];
                part->pts->delta[disp[4]] = L/(double)nL*1.75;
                ++disp[4];
              }
            }
          }
          ++rank_flag;
          if (rank_flag == mpi_size)
            rank_flag = 0;
        }
        /*
          if (mpi_rank == 0) {
          for (int ii = 0; ii < disp[4]; ++ii) {
          cout << part->pts->xp[ii][0] << ","
          << part->pts->xp[ii][1] << ","
          << part->pts->xp[ii][2] << endl;
          }
          }
        */
      }
      else {
        assert((mode == 0)||(mode == 1));
        assert(L > 2.0*delta);
        const int nL = int((L-2.0*delta)/dt); // + 2*n
        const int nW = int(W/dt);
        if (mpi_rank_shared == 0) {
          // ===================================================
          // surface points: xsp_ff...
          // ===================================================
          const int nsp_ff = spoli_s + 4 + 8*(n-1);
          // copy the links into xsp_ff first...
          for (int sol = 0; sol < spoli_s; ++sol) {
            const int isp = spoli_v[sol];
            FOR_I3 part->ff_surface->xsp[disp[2]+sol][i] = partVec[ipart]->surface->xsp[isp][i];
          }
          // then the layers start at spoli_i[1]...
          int isp_ff = spoli_s;
          for (int k = 0; k < n-1; ++k) {
            const double dl0 = dn*(pow(sf,k)-1.0)/(sf-1.0);
            const double dl1 = dn*(pow(sf,k+1)-1.0)/(sf-1.0);
            FOR_I3 part->ff_surface->xsp[disp[2]+isp_ff][i] = x0[i] +     dl0*e1[i] + dl1*e2[i]; ++isp_ff;
            FOR_I3 part->ff_surface->xsp[disp[2]+isp_ff][i] = x1[i] +     dl0*e1[i] + dl1*e2[i]; ++isp_ff;
            FOR_I3 part->ff_surface->xsp[disp[2]+isp_ff][i] = x1[i] + (L-dl0)*e1[i] + dl1*e2[i]; ++isp_ff;
            FOR_I3 part->ff_surface->xsp[disp[2]+isp_ff][i] = x0[i] + (L-dl0)*e1[i] + dl1*e2[i]; ++isp_ff;
            FOR_I3 part->ff_surface->xsp[disp[2]+isp_ff][i] = x0[i] +     dl1*e1[i] + dl1*e2[i]; ++isp_ff;
            FOR_I3 part->ff_surface->xsp[disp[2]+isp_ff][i] = x1[i] +     dl1*e1[i] + dl1*e2[i]; ++isp_ff;
            FOR_I3 part->ff_surface->xsp[disp[2]+isp_ff][i] = x1[i] + (L-dl1)*e1[i] + dl1*e2[i]; ++isp_ff;
            FOR_I3 part->ff_surface->xsp[disp[2]+isp_ff][i] = x0[i] + (L-dl1)*e1[i] + dl1*e2[i]; ++isp_ff;
          }
          const double dl0 = dn*(pow(sf,n-1)-1.0)/(sf-1.0);
          const double dl1 = dn*(pow(sf,n)-1.0)/(sf-1.0);
          FOR_I3 part->ff_surface->xsp[disp[2]+isp_ff][i] = x0[i] +     dl0*e1[i] + dl1*e2[i]; ++isp_ff;
          FOR_I3 part->ff_surface->xsp[disp[2]+isp_ff][i] = x1[i] +     dl0*e1[i] + dl1*e2[i]; ++isp_ff;
          FOR_I3 part->ff_surface->xsp[disp[2]+isp_ff][i] = x1[i] + (L-dl0)*e1[i] + dl1*e2[i]; ++isp_ff;
          FOR_I3 part->ff_surface->xsp[disp[2]+isp_ff][i] = x0[i] + (L-dl0)*e1[i] + dl1*e2[i]; ++isp_ff;
          assert(isp_ff == nsp_ff);
          /*
            if (mpi_rank == 0) {
            for (int isp_ff = 0; isp_ff < nsp_ff; ++isp_ff) {
            cout << part->ff_surface->xsp[disp[2]+isp_ff][0] << ","
            << part->ff_surface->xsp[disp[2]+isp_ff][1] << ","
            << part->ff_surface->xsp[disp[2]+isp_ff][2] << endl;
            }
            }
          */

          // ===================================================
          // surface tris: xst_ff...
          // ===================================================
          const int nst_ff = spoli_s + 4 + 2 + 16*(n-1);
          // 1. the lip tris are special because they have to respect both discretizations...
          int ist_ff = 0;

          // tris from first layer...

          part->ff_surface->spost[ist_ff+disp[3]][0] = disp[2]+sols[0];
          part->ff_surface->spost[ist_ff+disp[3]][1] = disp[2]+spoli_s+4;
          part->ff_surface->spost[ist_ff+disp[3]][2] = disp[2]+spoli_s+0;
          part->ff_surface_dxost[ist_ff+disp[3]] = dn;
          ++ist_ff;
          part->ff_surface->spost[ist_ff+disp[3]][0] = disp[2]+sols[1];
          part->ff_surface->spost[ist_ff+disp[3]][1] = disp[2]+spoli_s+0;
          part->ff_surface->spost[ist_ff+disp[3]][2] = disp[2]+spoli_s+1;
          part->ff_surface_dxost[ist_ff+disp[3]] = dn;
          ++ist_ff;
          part->ff_surface->spost[ist_ff+disp[3]][0] = disp[2]+sols[1];
          part->ff_surface->spost[ist_ff+disp[3]][1] = disp[2]+spoli_s+1;
          part->ff_surface->spost[ist_ff+disp[3]][2] = disp[2]+spoli_s+5;
          part->ff_surface_dxost[ist_ff+disp[3]] = dn;
          ++ist_ff;
          part->ff_surface->spost[ist_ff+disp[3]][0] = disp[2]+sols[1];
          part->ff_surface->spost[ist_ff+disp[3]][1] = disp[2]+spoli_s+5;
          part->ff_surface->spost[ist_ff+disp[3]][2] = disp[2]+spoli_s+6;
          part->ff_surface_dxost[ist_ff+disp[3]] = dn;
          ++ist_ff;
          part->ff_surface->spost[ist_ff+disp[3]][0] = disp[2]+sols[2];
          part->ff_surface->spost[ist_ff+disp[3]][1] = disp[2]+spoli_s+6;
          part->ff_surface->spost[ist_ff+disp[3]][2] = disp[2]+spoli_s+2;
          part->ff_surface_dxost[ist_ff+disp[3]] = dn;
          ++ist_ff;
          part->ff_surface->spost[ist_ff+disp[3]][0] = disp[2]+sols[3];
          part->ff_surface->spost[ist_ff+disp[3]][1] = disp[2]+spoli_s+2;
          part->ff_surface->spost[ist_ff+disp[3]][2] = disp[2]+spoli_s+3;
          part->ff_surface_dxost[ist_ff+disp[3]] = dn;
          ++ist_ff;
          part->ff_surface->spost[ist_ff+disp[3]][0] = disp[2]+sols[3];
          part->ff_surface->spost[ist_ff+disp[3]][1] = disp[2]+spoli_s+3;
          part->ff_surface->spost[ist_ff+disp[3]][2] = disp[2]+spoli_s+7;
          part->ff_surface_dxost[ist_ff+disp[3]] = dn;
          ++ist_ff;
          part->ff_surface->spost[ist_ff+disp[3]][0] = disp[2]+sols[3];
          part->ff_surface->spost[ist_ff+disp[3]][1] = disp[2]+spoli_s+7;
          part->ff_surface->spost[ist_ff+disp[3]][2] = disp[2]+spoli_s+4;
          part->ff_surface_dxost[ist_ff+disp[3]] = dn;
          ++ist_ff;

          // tris from lip...

          int isp_ff_prev = sols[0];
          for (int isp_ff = isp_ff_prev+1; isp_ff <= sols[1]; ++isp_ff) {
            part->ff_surface->spost[ist_ff+disp[3]][0] = disp[2]+isp_ff_prev;
            part->ff_surface->spost[ist_ff+disp[3]][1] = disp[2]+spoli_s;
            part->ff_surface->spost[ist_ff+disp[3]][2] = disp[2]+isp_ff;
            part->ff_surface_dxost[ist_ff+disp[3]] = dn;
            ++ist_ff;
            isp_ff_prev = isp_ff;
          }
          isp_ff_prev = sols[1];
          for (int isp_ff = isp_ff_prev+1; isp_ff <= sols[2]; ++isp_ff) {
            part->ff_surface->spost[ist_ff+disp[3]][0] = disp[2]+isp_ff_prev;
            part->ff_surface->spost[ist_ff+disp[3]][1] = disp[2]+spoli_s+6;
            part->ff_surface->spost[ist_ff+disp[3]][2] = disp[2]+isp_ff;
            part->ff_surface_dxost[ist_ff+disp[3]] = dn;
            ++ist_ff;
            isp_ff_prev = isp_ff;
          }
          if (mode == 0) {
            isp_ff_prev = sols[2];
            for (int isp_ff = isp_ff_prev+1; isp_ff <= sols[3]; ++isp_ff) {
              part->ff_surface->spost[ist_ff+disp[3]][0] = disp[2]+isp_ff_prev;
              part->ff_surface->spost[ist_ff+disp[3]][1] = disp[2]+spoli_s+2;
              part->ff_surface->spost[ist_ff+disp[3]][2] = disp[2]+isp_ff;
              part->ff_surface_dxost[ist_ff+disp[3]] = dn;
              ++ist_ff;
              isp_ff_prev = isp_ff;
            }
            isp_ff_prev = sols[3];
            for (int isp_ff = isp_ff_prev+1; isp_ff < spoli_s; ++isp_ff) {
              part->ff_surface->spost[ist_ff+disp[3]][0] = disp[2]+isp_ff_prev;
              part->ff_surface->spost[ist_ff+disp[3]][1] = disp[2]+spoli_s+4;
              part->ff_surface->spost[ist_ff+disp[3]][2] = disp[2]+isp_ff;
              part->ff_surface_dxost[ist_ff+disp[3]] = dn;
              ++ist_ff;
              isp_ff_prev = isp_ff;
            }
            isp_ff_prev = 0;
            for (int isp_ff = isp_ff_prev+1; isp_ff <= sols[0]; ++isp_ff) {
              part->ff_surface->spost[ist_ff+disp[3]][0] = disp[2]+isp_ff_prev;
              part->ff_surface->spost[ist_ff+disp[3]][1] = disp[2]+spoli_s+4;
              part->ff_surface->spost[ist_ff+disp[3]][2] = disp[2]+isp_ff;
              part->ff_surface_dxost[ist_ff+disp[3]] = dn;
              ++ist_ff;
              isp_ff_prev = isp_ff;
            }
            part->ff_surface->spost[ist_ff+disp[3]][0] = disp[2]+spoli_s-1;
            part->ff_surface->spost[ist_ff+disp[3]][1] = disp[2]+spoli_s+4;
            part->ff_surface->spost[ist_ff+disp[3]][2] = disp[2];
            part->ff_surface_dxost[ist_ff+disp[3]] = dn;
            ++ist_ff;
          }
          else {
            assert(mode == 1);
            isp_ff_prev = sols[2];
            for (int isp_ff = isp_ff_prev+1; isp_ff < spoli_s; ++isp_ff) {
              part->ff_surface->spost[ist_ff+disp[3]][0] = disp[2]+isp_ff_prev;
              part->ff_surface->spost[ist_ff+disp[3]][1] = disp[2]+spoli_s+2;
              part->ff_surface->spost[ist_ff+disp[3]][2] = disp[2]+isp_ff;
              part->ff_surface_dxost[ist_ff+disp[3]] = dn;
              ++ist_ff;
              isp_ff_prev = isp_ff;
            }
            isp_ff_prev = 0;
            for (int isp_ff = isp_ff_prev+1; isp_ff <= sols[3]; ++isp_ff) {
              part->ff_surface->spost[ist_ff+disp[3]][0] = disp[2]+isp_ff_prev;
              part->ff_surface->spost[ist_ff+disp[3]][1] = disp[2]+spoli_s+4;
              part->ff_surface->spost[ist_ff+disp[3]][2] = disp[2]+isp_ff;
              part->ff_surface_dxost[ist_ff+disp[3]] = dn;
              ++ist_ff;
              isp_ff_prev = isp_ff;
            }
            isp_ff_prev = sols[3];
            for (int isp_ff = isp_ff_prev+1; isp_ff <= sols[0]; ++isp_ff) {
              part->ff_surface->spost[ist_ff+disp[3]][0] = disp[2]+isp_ff_prev;
              part->ff_surface->spost[ist_ff+disp[3]][1] = disp[2]+spoli_s+4;
              part->ff_surface->spost[ist_ff+disp[3]][2] = disp[2]+isp_ff;
              part->ff_surface_dxost[ist_ff+disp[3]] = dn;
              ++ist_ff;
              isp_ff_prev = isp_ff;
            }
            part->ff_surface->spost[ist_ff+disp[3]][0] = disp[2]+spoli_s-1;
            part->ff_surface->spost[ist_ff+disp[3]][1] = disp[2]+spoli_s+2;
            part->ff_surface->spost[ist_ff+disp[3]][2] = disp[2];
            part->ff_surface_dxost[ist_ff+disp[3]] = dn;
            ++ist_ff;
          }
          for (int k = 0; k < n-1; ++k) {
            int isp_ff = spoli_s+8*k;
            const double dl0 = dn*pow(sf,k);
            const double dl1 = dn*pow(sf,k+1);
            // horizontal...
            part->ff_surface->spost[ist_ff+disp[3]][0] = disp[2]+isp_ff;
            part->ff_surface->spost[ist_ff+disp[3]][1] = disp[2]+isp_ff+4;
            part->ff_surface->spost[ist_ff+disp[3]][2] = disp[2]+isp_ff+5;
            part->ff_surface_dxost[ist_ff+disp[3]] = dl0;
            ++ist_ff;
            part->ff_surface->spost[ist_ff+disp[3]][0] = disp[2]+isp_ff;
            part->ff_surface->spost[ist_ff+disp[3]][1] = disp[2]+isp_ff+5;
            part->ff_surface->spost[ist_ff+disp[3]][2] = disp[2]+isp_ff+1;
            part->ff_surface_dxost[ist_ff+disp[3]] = dl0;
            ++ist_ff;
            part->ff_surface->spost[ist_ff+disp[3]][0] = disp[2]+isp_ff+2;
            part->ff_surface->spost[ist_ff+disp[3]][1] = disp[2]+isp_ff+6;
            part->ff_surface->spost[ist_ff+disp[3]][2] = disp[2]+isp_ff+7;
            part->ff_surface_dxost[ist_ff+disp[3]] = dl0;
            ++ist_ff;
            part->ff_surface->spost[ist_ff+disp[3]][0] = disp[2]+isp_ff+2;
            part->ff_surface->spost[ist_ff+disp[3]][1] = disp[2]+isp_ff+7;
            part->ff_surface->spost[ist_ff+disp[3]][2] = disp[2]+isp_ff+3;
            part->ff_surface_dxost[ist_ff+disp[3]] = dl0;
            ++ist_ff;
            // vertical...
            part->ff_surface->spost[ist_ff+disp[3]][0] = disp[2]+isp_ff+4;
            part->ff_surface->spost[ist_ff+disp[3]][1] = disp[2]+isp_ff+8;
            part->ff_surface->spost[ist_ff+disp[3]][2] = disp[2]+isp_ff+9;
            part->ff_surface_dxost[ist_ff+disp[3]] = dl1;
            ++ist_ff;
            part->ff_surface->spost[ist_ff+disp[3]][0] = disp[2]+isp_ff+4;
            part->ff_surface->spost[ist_ff+disp[3]][1] = disp[2]+isp_ff+9;
            part->ff_surface->spost[ist_ff+disp[3]][2] = disp[2]+isp_ff+5;
            part->ff_surface_dxost[ist_ff+disp[3]] = dl1;
            ++ist_ff;
            part->ff_surface->spost[ist_ff+disp[3]][0] = disp[2]+isp_ff+6;
            part->ff_surface->spost[ist_ff+disp[3]][1] = disp[2]+isp_ff+10;
            part->ff_surface->spost[ist_ff+disp[3]][2] = disp[2]+isp_ff+11;
            part->ff_surface_dxost[ist_ff+disp[3]] = dl1;
            ++ist_ff;
            part->ff_surface->spost[ist_ff+disp[3]][0] = disp[2]+isp_ff+6;
            part->ff_surface->spost[ist_ff+disp[3]][1] = disp[2]+isp_ff+11;
            part->ff_surface->spost[ist_ff+disp[3]][2] = disp[2]+isp_ff+7;
            part->ff_surface_dxost[ist_ff+disp[3]] = dl1;
            ++ist_ff;
            if (k < n-2) {
              part->ff_surface->spost[ist_ff+disp[3]][0] = disp[2]+isp_ff+4;
              part->ff_surface->spost[ist_ff+disp[3]][1] = disp[2]+isp_ff+12;
              part->ff_surface->spost[ist_ff+disp[3]][2] = disp[2]+isp_ff+8;
              part->ff_surface_dxost[ist_ff+disp[3]] = dl1;
              ++ist_ff;
              part->ff_surface->spost[ist_ff+disp[3]][0] = disp[2]+isp_ff+4;
              part->ff_surface->spost[ist_ff+disp[3]][1] = disp[2]+isp_ff+7;
              part->ff_surface->spost[ist_ff+disp[3]][2] = disp[2]+isp_ff+12;
              part->ff_surface_dxost[ist_ff+disp[3]] = dl1;
              ++ist_ff;
              part->ff_surface->spost[ist_ff+disp[3]][0] = disp[2]+isp_ff+7;
              part->ff_surface->spost[ist_ff+disp[3]][1] = disp[2]+isp_ff+15;
              part->ff_surface->spost[ist_ff+disp[3]][2] = disp[2]+isp_ff+12;
              part->ff_surface_dxost[ist_ff+disp[3]] = dl1;
              ++ist_ff;
              part->ff_surface->spost[ist_ff+disp[3]][0] = disp[2]+isp_ff+7;
              part->ff_surface->spost[ist_ff+disp[3]][1] = disp[2]+isp_ff+11;
              part->ff_surface->spost[ist_ff+disp[3]][2] = disp[2]+isp_ff+15;
              part->ff_surface_dxost[ist_ff+disp[3]] = dl1;
              ++ist_ff;
              part->ff_surface->spost[ist_ff+disp[3]][0] = disp[2]+isp_ff+5;
              part->ff_surface->spost[ist_ff+disp[3]][1] = disp[2]+isp_ff+9;
              part->ff_surface->spost[ist_ff+disp[3]][2] = disp[2]+isp_ff+13;
              part->ff_surface_dxost[ist_ff+disp[3]] = dl1;
              ++ist_ff;
              part->ff_surface->spost[ist_ff+disp[3]][0] = disp[2]+isp_ff+5;
              part->ff_surface->spost[ist_ff+disp[3]][1] = disp[2]+isp_ff+13;
              part->ff_surface->spost[ist_ff+disp[3]][2] = disp[2]+isp_ff+14;
              part->ff_surface_dxost[ist_ff+disp[3]] = dl1;
              ++ist_ff;
              part->ff_surface->spost[ist_ff+disp[3]][0] = disp[2]+isp_ff+5;
              part->ff_surface->spost[ist_ff+disp[3]][1] = disp[2]+isp_ff+14;
              part->ff_surface->spost[ist_ff+disp[3]][2] = disp[2]+isp_ff+6;
              part->ff_surface_dxost[ist_ff+disp[3]] = dl1;
              ++ist_ff;
              part->ff_surface->spost[ist_ff+disp[3]][0] = disp[2]+isp_ff+14;
              part->ff_surface->spost[ist_ff+disp[3]][1] = disp[2]+isp_ff+10;
              part->ff_surface->spost[ist_ff+disp[3]][2] = disp[2]+isp_ff+6;
              part->ff_surface_dxost[ist_ff+disp[3]] = dl1;
              ++ist_ff;
            }
            else {
              part->ff_surface->spost[ist_ff+disp[3]][0] = disp[2]+isp_ff+4;
              part->ff_surface->spost[ist_ff+disp[3]][1] = disp[2]+isp_ff+11;
              part->ff_surface->spost[ist_ff+disp[3]][2] = disp[2]+isp_ff+8;
              part->ff_surface_dxost[ist_ff+disp[3]] = dl1;
              ++ist_ff;
              part->ff_surface->spost[ist_ff+disp[3]][0] = disp[2]+isp_ff+4;
              part->ff_surface->spost[ist_ff+disp[3]][1] = disp[2]+isp_ff+7;
              part->ff_surface->spost[ist_ff+disp[3]][2] = disp[2]+isp_ff+11;
              part->ff_surface_dxost[ist_ff+disp[3]] = dl1;
              ++ist_ff;
              part->ff_surface->spost[ist_ff+disp[3]][0] = disp[2]+isp_ff+5;
              part->ff_surface->spost[ist_ff+disp[3]][1] = disp[2]+isp_ff+9;
              part->ff_surface->spost[ist_ff+disp[3]][2] = disp[2]+isp_ff+10;
              part->ff_surface_dxost[ist_ff+disp[3]] = dl1;
              ++ist_ff;
              part->ff_surface->spost[ist_ff+disp[3]][0] = disp[2]+isp_ff+5;
              part->ff_surface->spost[ist_ff+disp[3]][1] = disp[2]+isp_ff+10;
              part->ff_surface->spost[ist_ff+disp[3]][2] = disp[2]+isp_ff+6;
              part->ff_surface_dxost[ist_ff+disp[3]] = dl1;
              ++ist_ff;
              // top
              part->ff_surface->spost[ist_ff+disp[3]][0] = disp[2]+isp_ff+8;
              part->ff_surface->spost[ist_ff+disp[3]][1] = disp[2]+isp_ff+11;
              part->ff_surface->spost[ist_ff+disp[3]][2] = disp[2]+isp_ff+10;
              part->ff_surface_dxost[ist_ff+disp[3]] = dl1;
              ++ist_ff;
              part->ff_surface->spost[ist_ff+disp[3]][0] = disp[2]+isp_ff+8;
              part->ff_surface->spost[ist_ff+disp[3]][1] = disp[2]+isp_ff+10;
              part->ff_surface->spost[ist_ff+disp[3]][2] = disp[2]+isp_ff+9;
              part->ff_surface_dxost[ist_ff+disp[3]] = dl1;
              ++ist_ff;
            }
          }
          //if (mpi_rank == 0) part->ff_surface->writeTecplotNoZones("ff_surface.plt");
          disp[2] += nsp_ff;
          assert(ist_ff == nst_ff);
          disp[3] += nst_ff;
        }
        // ==================================================================
        // now the points...
        // ==================================================================
        assert(part->pts);
        int rank_flag = 0;
        for (int i = 0; i < n; ++i) {
          const double yp = dn*(pow(sf,0.5+double(i))-1.0)/(sf-1.0);
          const double dx = dn*pow(sf,i);
          const double this_nW = int(W/dx);
          //if (i == n-1) assert(this_nW == nW);
          for (int k = 0; k < this_nW; ++k) {
            const double xp = (double(k)+0.5)/double(this_nW)*W;
            for (int j = 0; j <= i; ++j) {
              const double zp = dn*(pow(sf,0.5+double(j))-1.0)/(sf-1.0);
              if (rank_flag == mpi_rank) {
                part->pts->xp[disp[4]][0] = x0[0] + xp*e0[0] + yp*e1[0] + zp*e2[0];
                part->pts->xp[disp[4]][1] = x0[1] + xp*e0[1] + yp*e1[1] + zp*e2[1];
                part->pts->xp[disp[4]][2] = x0[2] + xp*e0[2] + yp*e1[2] + zp*e2[2];
                part->pts->delta[disp[4]] = 1.75*dx; // i.e. a little more than sqrt(3)
                ++disp[4];
              }
            }
          }
          ++rank_flag;
          if (rank_flag == mpi_size)
            rank_flag = 0;
        }
        for (int i = 0; i < nL; ++i) {
          const double yp = delta + (double(i)+0.5)/double(nL)*(L-2.0*delta);
          for (int k = 0; k < nW; ++k) {
            const double xp = (double(k)+0.5)/double(nW)*W;
            for (int j = 0; j < n; ++j) {
              const double zp = dn*(pow(sf,0.5+double(j))-1.0)/(sf-1.0);
              if (rank_flag == mpi_rank) {
                part->pts->xp[disp[4]][0] = x0[0] + xp*e0[0] + yp*e1[0] + zp*e2[0];
                part->pts->xp[disp[4]][1] = x0[1] + xp*e0[1] + yp*e1[1] + zp*e2[1];
                part->pts->xp[disp[4]][2] = x0[2] + xp*e0[2] + yp*e1[2] + zp*e2[2];
                part->pts->delta[disp[4]] = (L-2.0*delta)/nL*1.75;
                ++disp[4];
              }
            }
          }
          ++rank_flag;
          if (rank_flag == mpi_size)
            rank_flag = 0;
        }
        for (int i = n-1; i >= 0; --i) {
          const double yp = L - dn*(pow(sf,0.5+double(i))-1.0)/(sf-1.0);
          const double dx = dn*pow(sf,i);
          const double this_nW = int(W/dx);
          //if (i == n-1) assert(this_nW == nW);
          for (int k = 0; k < this_nW; ++k) {
            const double xp = (double(k)+0.5)/double(this_nW)*W;
            for (int j = 0; j <= i; ++j) {
              const double zp = dn*(pow(sf,0.5+double(j))-1.0)/(sf-1.0);
              if (rank_flag == mpi_rank) {
                part->pts->xp[disp[4]][0] = x0[0] + xp*e0[0] + yp*e1[0] + zp*e2[0];
                part->pts->xp[disp[4]][1] = x0[1] + xp*e0[1] + yp*e1[1] + zp*e2[1];
                part->pts->xp[disp[4]][2] = x0[2] + xp*e0[2] + yp*e1[2] + zp*e2[2];
                part->pts->delta[disp[4]] = 1.75*dx; // i.e. a little more than sqrt(3)
                ++disp[4];
              }
            }
          }
          ++rank_flag;
          if (rank_flag == mpi_size)
            rank_flag = 0;
        }
        /*
          if (mpi_rank == 0) {
          for (int ii = 0; ii < disp[4]; ++ii) {
          cout << part->pts->xp[ii][0] << ","
          << part->pts->xp[ii][1] << ","
          << part->pts->xp[ii][2] << endl;
          }
          }
        */
      }
    }

  };

  void buildGroupsForFlaggedZones(vector<GroupData>& groupDataVec,const set<int>& zones) {

    if (mpi_rank == 0) {
      cout << " > buildGroupsForFlaggedZones: ";
      for (set<int>::iterator it = zones.begin(); it != zones.end(); ++it)
        cout << *it << " ";
      cout << endl;
    }

    /*
      {
      
      // HACK to speed up grouping...
      
      // global zone flag indices in zones has to be converted to local surface zone
      // flag for each surface...

      int * zn_flag = new int[zoneVec.size()];
      for (int izn = 0; izn < zoneVec.size(); ++izn)
      zn_flag[izn] = 0;

      for (set<int>::iterator it = zones.begin(); it != zones.end(); ++it) {
      const int izn = *it;
      assert((izn >= 0)&&(izn < zoneVec.size()));
      zn_flag[izn] = 1;
      }

      for (int ipart = 0; ipart < partVec.size(); ++ipart) {
      if (partVec[ipart]->surface) {
      SurfaceShm * surface = partVec[ipart]->surface;
      set<int> surface_zones;
      for (int isz = 0; isz < surface->zoneVec.size(); ++isz) {
      const int izn = partVec[ipart]->znosz[isz]; // convert from surface zone (isz) to global zone (izn)
      assert((izn >= 0)&&(izn < zoneVec.size()));
      if (zn_flag[izn] != 0) {
      surface_zones.insert(isz);
      }
      }
      if (!surface_zones.empty()) {
      surface->flagDisjointGroups(surface_zones);
      }
      }
      }

      delete[] zn_flag;

      //MPI_Pause("HOW WAS THAT?");

      }

      return;
    */



    
    // this routine should build all disjoint groups associated with the
    // passed partData zone indices, and calculate group data for each. It
    // also sets the shm surface_st_flag to the group index...

    int ngr = 0;
    int sizes[2] = {0,0}; // ibuf,dbuf sizes
    int *ibuf;
    double *dbuf;
    if (mpi_rank_shared == 0) {

      vector<pair<int,int> > flagVec;
      map<const pair<int,int>,pair<int,int> > edgeMap;
      vector<int> paogr;
      vector<int> spoli_i;
      spoli_i.push_back(0);
      vector<int> spoli_v;
      vector<int> groli;

      int * zone_flag = new int[zoneVec.size()];
      for (int izn = 0; izn < zoneVec.size(); ++izn)
        zone_flag[izn] = 0;

      for (set<int>::iterator it = zones.begin(); it != zones.end(); ++it) {
        const int izn = *it;
        assert((izn >= 0)&&(izn < zoneVec.size()));
        zone_flag[izn] = 1;
      }

      for (int ipart = 0; ipart < partVec.size(); ++ipart) {
        if (partVec[ipart]->surface) {

          SurfaceShm * surface = partVec[ipart]->surface;
          assert(edgeMap.empty());
          assert(flagVec.empty());
          for (int ist = 0; ist < surface->nst; ++ist) {
            {
              // check that group is clear...
              int igr;
              assert(!partVec[ipart]->getGroupForSt(igr,ist));
            }
            const int isz = surface->znost[ist];
            const int izn = partVec[ipart]->znosz[isz]; // convert from surface zone (isz) to global zone (izn)
            if (zone_flag[izn] != 0) {
              int count = flagVec.size();
              flagVec.push_back(pair<int,int>(ist,count));
              FOR_I3 {
                const int isp0 = surface->spost[ist][i];
                const int isp1 = surface->spost[ist][(i+1)%3];
                map<const pair<int,int>,pair<int,int> >::iterator iter = edgeMap.find(pair<int,int>(isp1,isp0));
                if (iter != edgeMap.end()) {
                  int nbr_count = iter->second.first;
                  edgeMap.erase(iter);
                  while (nbr_count != flagVec[nbr_count].second)
                    nbr_count = flagVec[nbr_count].second;
                  if (nbr_count < count) {
                    assert(flagVec[count].second == count);
                    flagVec[count].second = nbr_count;
                    count = nbr_count;
                  }
                  else if (count < nbr_count) {
                    assert(flagVec[nbr_count].second == nbr_count);
                    flagVec[nbr_count].second = count;
                    nbr_count = count;
                  }
                }
                else {
                  edgeMap[pair<int,int>(isp0,isp1)] = pair<int,int>(count,ist);
                }
              }
            }
          }

          if (mpi_rank == 0) cout << " > ipart, nst flagged: " << ipart << " " << flagVec.size() << endl;

          for (int ii = 0, limit = flagVec.size(); ii < limit; ++ii) {
            int count = flagVec[ii].second;
            while ((count >= 0)&&(count != flagVec[count].second))
              count = flagVec[count].second;
            if (count >= 0) {
              paogr.push_back(ipart);
              ++ngr;
              flagVec[ii].second = -ngr; // -1,-2,...
            }
            else {
              flagVec[ii].second = count;
            }
          }

          for (int ii = 0, limit = flagVec.size(); ii < limit; ++ii) {
            const int ist = flagVec[ii].first; assert((ist >= 0)&&(ist < surface->nst));
            const int igr = -flagVec[ii].second-1; assert((igr >= 0)&&(igr < ngr));
            partVec[ipart]->setGroupForSt(igr,ist);
          }

          // now find the edge loops associated with open edges...
          map<const pair<int,int>,int> linkMap;
          for (map<const pair<int,int>,pair<int,int> >::iterator iter = edgeMap.begin(); iter != edgeMap.end(); ++iter) {
            const int ist = iter->second.second;
            int igr;
            if (partVec[ipart]->getGroupForSt(igr,ist)) {
              assert((igr >= 0)&&(igr < ngr));
              const int isp0 = iter->first.first;
              const int isp1 = iter->first.second;
              linkMap.insert(pair<pair<int,int>,int>(pair<int,int>(igr,isp0),isp1));
            }
          }

          // make space in spoli_v for the nodes...
          int sol = spoli_v.size();
          spoli_v.resize(sol+edgeMap.size());
          edgeMap.clear();

          while (!linkMap.empty()) {

            //double xc[3] = { 0.0, 0.0, 0.0 };
            //double wgt = 0.0;
            map<const pair<int,int>,int>::iterator iter = linkMap.begin();
            assert(iter != linkMap.end());
            const int igr = iter->first.first;

            if (mpi_rank == 0) cout << " > building link for group: " << igr << endl;

            const int start = iter->first.second;
            int next = start;
            do {

              iter = linkMap.find(pair<int,int>(igr,next));
              assert(iter != linkMap.end());
              const int prev = next;
              next = iter->second;
              linkMap.erase(iter);

              spoli_v[sol++] = prev;

            } while (next != start);

            spoli_i.push_back(sol);
            groli.push_back(igr);

          }

          assert(sol == spoli_v.size());
          flagVec.clear();

        }
      }
      delete[] zone_flag;

      if (mpi_rank == 0) cout << "spoli_i.size(): " << spoli_i.size() << endl;

      // allocate group vec on shm rank 0 only...

      assert(groupDataVec.empty());
      groupDataVec.resize(ngr);

      for (int igr = 0; igr < ngr; ++igr) {
        groupDataVec[igr].igr = igr;
        groupDataVec[igr].ipart = paogr[igr];
        assert(groupDataVec[igr].nlink == 0);
        assert(groupDataVec[igr].spoli_s == 0);
      }

      for (int ilink = 0; ilink < groli.size(); ++ilink) {
        const int igr = groli[ilink];
        const int nsp = spoli_i[ilink+1] - spoli_i[ilink];
        ++groupDataVec[igr].nlink;
        groupDataVec[igr].spoli_s += nsp;
      }

      for (int igr = 0; igr < ngr; ++igr) {
        assert(groupDataVec[igr].spoli_i == NULL);
        groupDataVec[igr].spoli_i = new int[groupDataVec[igr].nlink+1];
        groupDataVec[igr].spoli_i[0] = 0;
        groupDataVec[igr].l_li = new double[groupDataVec[igr].nlink];
        groupDataVec[igr].n_li = new double[groupDataVec[igr].nlink][3];
        groupDataVec[igr].x_li = new double[groupDataVec[igr].nlink][3];
        groupDataVec[igr].d_li = new double[groupDataVec[igr].nlink];
        groupDataVec[igr].drms_li = new double[groupDataVec[igr].nlink];
        for (int il = 0; il < groupDataVec[igr].nlink; ++il) {
          groupDataVec[igr].l_li[il] = 0.0;
          FOR_I3 groupDataVec[igr].n_li[il][i] = 0.0;
          FOR_I3 groupDataVec[igr].x_li[il][i] = 0.0;
          groupDataVec[igr].d_li[il] = 0.0;
          groupDataVec[igr].drms_li[il] = 0.0;
        }
        groupDataVec[igr].nlink = 0;
        groupDataVec[igr].spoli_v = new int[groupDataVec[igr].spoli_s];
        groupDataVec[igr].spoli_s = 0;
      }

      for (int ilink = 0; ilink < groli.size(); ++ilink) {
        const int igr = groli[ilink];
        const int nsp = spoli_i[ilink+1] - spoli_i[ilink];
        groupDataVec[igr].spoli_i[groupDataVec[igr].nlink+1] = groupDataVec[igr].spoli_i[groupDataVec[igr].nlink] + nsp;
        ++groupDataVec[igr].nlink;
        for (int sol = spoli_i[ilink]; sol != spoli_i[ilink+1]; ++sol) {
          groupDataVec[igr].spoli_v[groupDataVec[igr].spoli_s] = spoli_v[sol];
          ++groupDataVec[igr].spoli_s;
        }
      }

      // characterize link geometry...
      for (int igr = 0; igr < ngr; ++igr) {
        for (int il = 0; il < groupDataVec[igr].nlink; ++il) {
          SurfaceShm * surface = partVec[groupDataVec[igr].ipart]->surface;

          // first pass to get length and centroid...
          const int last = groupDataVec[igr].spoli_v[groupDataVec[igr].spoli_i[il+1]-1];
          int next = last;
          for (int sol = groupDataVec[igr].spoli_i[il]; sol != groupDataVec[igr].spoli_i[il+1]; ++sol) {
            const int prev = next;
            next = groupDataVec[igr].spoli_v[sol];

            // add segment prev->next to xc,wgt...
            const double dl = DIST(surface->xsp[next],surface->xsp[prev]);
            FOR_I3 groupDataVec[igr].x_li[il][i] += dl*(surface->xsp[prev][i]+surface->xsp[next][i]); // don't forget factor of 2 when normalizing!
            groupDataVec[igr].l_li[il] += dl;

            if ((prev != last)&&(next != last)) {
              const double n2[3] = TRI_NORMAL_2(surface->xsp[last],surface->xsp[prev],surface->xsp[next]);
              FOR_I3 groupDataVec[igr].n_li[il][i] += n2[i];
            }
          }
          FOR_I3 groupDataVec[igr].x_li[il][i] /= 2.0*groupDataVec[igr].l_li[il];
          FOR_I3 groupDataVec[igr].n_li[il][i] *= 0.5;

          // second pass to get mean and variance of distance...
          next = last;
          for (int sol = groupDataVec[igr].spoli_i[il]; sol != groupDataVec[igr].spoli_i[il+1]; ++sol) {
            const int prev = next;
            next = groupDataVec[igr].spoli_v[sol];

            // add segment prev->next to xc,wgt...
            const double dl = DIST(surface->xsp[next],surface->xsp[prev]);
            double xc[3]; FOR_I3 xc[i] = 0.5*(surface->xsp[prev][i]+surface->xsp[next][i]);
            const double d2 = DIST2(xc,groupDataVec[igr].x_li[il]);
            groupDataVec[igr].d_li[il] += dl*sqrt(d2);
            groupDataVec[igr].drms_li[il] += dl*d2;
          }
          groupDataVec[igr].d_li[il] /= groupDataVec[igr].l_li[il];
          groupDataVec[igr].drms_li[il] = sqrt(max(groupDataVec[igr].drms_li[il]/groupDataVec[igr].l_li[il] -
                                                   groupDataVec[igr].d_li[il]*groupDataVec[igr].d_li[il],0.0));
        }
      }

      // populate buffers to communicate across node...

      // ints...
      sizes[0] += 1; // ngr
      sizes[0] += 3*ngr; // ipart,nlink,spoli_s
      for (int igr = 0; igr < ngr; ++igr) {
        sizes[0] += groupDataVec[igr].nlink;   // spoli_i skipping 0
        sizes[0] += groupDataVec[igr].spoli_s; // spoli_v
      }
      ibuf = new int[sizes[0]];
      int count = 0;
      ibuf[count++] = ngr;
      for (int igr = 0; igr < ngr; ++igr) {
        ibuf[count++] = groupDataVec[igr].ipart;
        ibuf[count++] = groupDataVec[igr].nlink;
        ibuf[count++] = groupDataVec[igr].spoli_s;
        for (int il = 0; il < groupDataVec[igr].nlink; ++il)
          ibuf[count++] = groupDataVec[igr].spoli_i[il+1];
        for (int sol = 0; sol < groupDataVec[igr].spoli_s; ++sol)
          ibuf[count++] = groupDataVec[igr].spoli_v[sol];
      }
      assert(count == sizes[0]);
      // doubles...
      for (int igr = 0; igr < ngr; ++igr)
        sizes[1] += 9*groupDataVec[igr].nlink; // l_li,n_li[3],x_li[3],d_li,drms_li
      dbuf = new double[sizes[1]];
      count = 0;
      for (int igr = 0; igr < ngr; ++igr) {
        for (int il = 0; il < groupDataVec[igr].nlink; ++il) {
          dbuf[count++] = groupDataVec[igr].l_li[il];
          FOR_I3 dbuf[count++] = groupDataVec[igr].n_li[il][i];
          FOR_I3 dbuf[count++] = groupDataVec[igr].x_li[il][i];
          dbuf[count++] = groupDataVec[igr].d_li[il];
          dbuf[count++] = groupDataVec[igr].drms_li[il];
        }
      }
      assert(count == sizes[1]);

      if (mpi_rank == 0) cout << " > ngr: " << ngr << endl;

    }
    MPI_Bcast(sizes,2,MPI_INT,0,mpi_comm_shared);
    if (mpi_rank_shared != 0) {
      ibuf = new int[sizes[0]];
      dbuf = new double[sizes[1]];
    }
    MPI_Bcast(ibuf,sizes[0],MPI_INT,0,mpi_comm_shared);
    MPI_Bcast(dbuf,sizes[1],MPI_DOUBLE,0,mpi_comm_shared);
    if (mpi_rank_shared != 0) {
      int count = 0;
      ngr = ibuf[count++];
      assert(groupDataVec.empty());
      groupDataVec.resize(ngr);
      for (int igr = 0; igr < ngr; ++igr) {
        groupDataVec[igr].igr = igr;
        groupDataVec[igr].ipart = ibuf[count++];
        groupDataVec[igr].nlink = ibuf[count++];
        groupDataVec[igr].spoli_i = new int[groupDataVec[igr].nlink+1];
        groupDataVec[igr].spoli_i[0] = 0;
        groupDataVec[igr].spoli_s = ibuf[count++];
        groupDataVec[igr].spoli_v = new int[groupDataVec[igr].spoli_s];
        for (int il = 0; il < groupDataVec[igr].nlink; ++il)
          groupDataVec[igr].spoli_i[1+il] = ibuf[count++];
        for (int sol = 0; sol < groupDataVec[igr].spoli_s; ++sol)
          groupDataVec[igr].spoli_v[sol] = ibuf[count++];
      }
      assert(count == sizes[0]);

      count = 0;
      for (int igr = 0; igr < ngr; ++igr) {
        groupDataVec[igr].l_li = new double[groupDataVec[igr].nlink];
        groupDataVec[igr].n_li = new double[groupDataVec[igr].nlink][3];
        groupDataVec[igr].x_li = new double[groupDataVec[igr].nlink][3];
        groupDataVec[igr].d_li = new double[groupDataVec[igr].nlink];
        groupDataVec[igr].drms_li = new double[groupDataVec[igr].nlink];
        for (int il = 0; il < groupDataVec[igr].nlink; ++il) {
          groupDataVec[igr].l_li[il] = dbuf[count++];
          FOR_I3 groupDataVec[igr].n_li[il][i] = dbuf[count++];
          FOR_I3 groupDataVec[igr].x_li[il][i] = dbuf[count++];
          groupDataVec[igr].d_li[il] = dbuf[count++];
          groupDataVec[igr].drms_li[il] = dbuf[count++];
        }
      }
      assert(count == sizes[1]);
    }
    delete[] ibuf;
    delete[] dbuf;

    // characterize group geometry...
    double (*my_gdogr)[7] = new double[ngr][7]; // n[3],x[3],area
    for (int igr = 0; igr < ngr; ++igr)
      FOR_I7 my_gdogr[igr][i] = 0.0;
    for (int ipart = 0; ipart < partVec.size(); ++ipart) {
      if (partVec[ipart]->surface) {
        SurfaceShm * surface = partVec[ipart]->surface;
        // stride through tris at mpi_size...
        for (int ist = mpi_rank; ist < surface->nst; ist += mpi_size) {
          int igr;
          if (partVec[ipart]->getGroupForSt(igr,ist)) {
            assert((igr >= 0)&&(igr < ngr));
            const double n2[3] = TRI_NORMAL_2(surface->xsp[surface->spost[ist][0]],
                                              surface->xsp[surface->spost[ist][1]],
                                              surface->xsp[surface->spost[ist][2]]);
            FOR_I3 my_gdogr[igr][i] += n2[i];
            const double area = MAG(n2);
            FOR_I3 my_gdogr[igr][3+i] += area*(surface->xsp[surface->spost[ist][0]][i] +
                                               surface->xsp[surface->spost[ist][1]][i] +
                                               surface->xsp[surface->spost[ist][2]][i]); // don't forget about the factor of 3
            my_gdogr[igr][6] += area;
          }
        }
      }
    }
    double (*gdogr)[7] = new double[ngr][7];
    MPI_Allreduce((double*)my_gdogr,(double*)gdogr,7*ngr,MPI_DOUBLE,MPI_SUM,mpi_comm);
    delete[] my_gdogr;
    for (int igr = 0; igr < ngr; ++igr) {
      FOR_I3 groupDataVec[igr].n[i] = 0.5*gdogr[igr][i];
      FOR_I3 groupDataVec[igr].x[i] = gdogr[igr][3+i]/(gdogr[igr][6]*3.0);
      groupDataVec[igr].area = 0.5*gdogr[igr][6];
    }
    delete[] gdogr;
    // note only shared roots have link data...

    if (mpi_rank == 0) {
      for (int igr = 0; igr < ngr; ++igr) {
        const int ipart = groupDataVec[igr].ipart;
        cout << "group: " << igr << " is part of ipart: " << ipart << " and has " << groupDataVec[igr].nlink << " links"
             << ", normal = " << COUT_VEC(groupDataVec[igr].n) << ", centroid = " << COUT_VEC(groupDataVec[igr].x) << ", area = " << groupDataVec[igr].area << endl;
        // take a look...
        for (int il = 0; il < groupDataVec[igr].nlink; ++il) {
          cout << "link: " << il << " has length: " << groupDataVec[igr].l_li[il] << ", centroid: " << COUT_VEC(groupDataVec[igr].x_li[il])
               << ", mean distance: " << groupDataVec[igr].d_li[il] << ", rms: " << groupDataVec[igr].drms_li[il]
               << ", normal: " << COUT_VEC(groupDataVec[igr].n_li[il]) << endl;
          /*
            char filename[128];
            sprintf(filename,"link.%04d.%04d.dat",igr,il);
            FILE * fp = fopen(filename,"w");
            for (int sol = groupDataVec[igr].spoli_i[il]; sol != groupDataVec[igr].spoli_i[il+1]; ++sol) {
            const int isp = groupDataVec[igr].spoli_v[sol];
            fprintf(fp,"%18.15le %18.15le %18.15le\n",partVec[ipart]->surface->xsp[isp][0],partVec[ipart]->surface->xsp[isp][1],partVec[ipart]->surface->xsp[isp][2]);
            }
            fclose(fp);
          */
        }
      }
    }

  }

  void clearGroups(vector<GroupData>& groupDataVec) {

    // the groups have allocated data, so call clear on each
    // member before clearing the vector...

    const int ngr = groupDataVec.size();
    for (int igr = 0; igr < ngr; ++igr)
      groupDataVec[igr].clear();
    groupDataVec.clear();

    // also go through ALL surface ist's and clear their settings if
    // set...

    if (mpi_rank_shared == 0) {
      for (int ipart = 0; ipart < partVec.size(); ++ipart) {
        if (partVec[ipart]->surface) {
          SurfaceShm * surface = partVec[ipart]->surface;
          for (int ist = 0; ist < surface->nst; ++ist) {
            int igr;
            if (partVec[ipart]->getGroupForSt(igr,ist)) {
              assert((igr >= 0)&&(igr < ngr));
              partVec[ipart]->clearGroupForSt(ist);
            }
          }
        }
      }
    }

    MPI_Barrier(mpi_comm_shared);

  }

  void queryZones(const set<int>& zoneIdSet) {

    assert(!zoneIdSet.empty());

    vector<GroupData> groupDataVec;
    buildGroupsForFlaggedZones(groupDataVec,zoneIdSet);

    WUI(INFO,"Specified zone(s) produced " << groupDataVec.size() << " disjoint group(s)...");

    for (int igr = 0; igr < groupDataVec.size(); ++igr) {
      if (groupDataVec[igr].isCylinder(-0.05)) { // passing a negative number for tol means a relative tolerance
        const double L = groupDataVec[igr].getCylinderL();
        const double r = groupDataVec[igr].getCylinderR();
        if (groupDataVec[igr].isRightCylinder()) {
          WUI(INFO,"group " << igr << ": looks like a right cylinder with r: " << r << ", L: " << L << ", L/D: " << L/(2.0*r));
        }
        else {
          WUI(INFO,"group " << igr << ": looks like a cylinder with r: " << r << ", L: " << L << ", L/D: " << L/(2.0*r));
        }
      }
      else if (groupDataVec[igr].isLogicalQuad()) {
        WUI(INFO,"group " << igr << ": looks like logical quad");
      }
      else {
        WUI(INFO,"group " << igr << ": don't know what this looks like");
      }
    }

    clearGroups(groupDataVec);

  }

  class GroupedPoint {
  public:
    double n[3];
    double x[3];
    double dist;
    int kmin,kmax;
    GroupedPoint() {
      // init to nonsense for checking
      FOR_I3 n[i] = HUGE_VAL;
      FOR_I3 x[i] = HUGE_VAL;
      dist = -HUGE_VAL;
      kmin = -1;
      kmax = -1;
    }
  };

  bool makeTriStrand(Part * part,Param * param,int& iarg) {

    set<int> zone_indices;
    bool b_zone_indices = false;
    bool b_dn = false;
    double dn;
    bool b_dt = false;
    double dt; // just used to hold dn_max/dxost
    bool b_n = false;
    int n;
    double factor = 1.0;
    bool b_branch = false;
    bool b_cull = false;

    int ierr = 0;
    while (iarg < param->size()) {
      string token = MiscUtils::toUpperCase(param->getString(iarg++));
      if ((token == "ZONE")||(token == "ZONES")||(token == "ZONE_NAME")||(token == "ZONE_NAMES")) {
        // expect a comma-delimited list of names...
        b_zone_indices = true;
        const string zoneNamesCsv = param->getString(iarg++);
        vector<string> tmp;
        MiscUtils::splitCsv(tmp,zoneNamesCsv);
        for (int ii = 0; ii < tmp.size(); ++ii) {
          const int izone = getZoneIndex(tmp[ii]); // get the unique zoneIndex in the partVec
          if (izone >= 0) {
            zone_indices.insert(izone);
          }
          else {
            WUI(WARN,"TRI_STRAND ZONE not recognized: " << tmp[ii]);
            ierr = -1;
          }
        }
      }
      else if (token == "FACTOR") {
        factor = param->getDouble(iarg++);
      }
      else if (token == "BRANCH") {
        b_branch = true;
      }
      else if (token == "CULL") {
        b_cull = true;
      }
      else if (token == "DN") {
        b_dn = true;
        dn = param->getDouble(iarg++);
      }
      else if ((token == "DT")||(token == "DXOST")) {
        b_dt = true;
        dt = param->getDouble(iarg++);
      }
      else if (token == "N") {
        b_n = true;
        n = param->getInt(iarg++);
      }
      else if (token == "NLAYERS") {
        part->ff_nlayersString = param->getString(iarg++);
      }
      else {
        WUI(WARN,"Unknown TRI_STRAND parameter " << token << ". skipping...");
      }
    }

    if (!b_dn) {
      WUI(WARN,"Error: TRI_STRAND missing min normal spacing DN <double>");
      ierr = -1;
    }
    if (!b_dt) {
      WUI(WARN,"Error: TRI_STRAND missing max normal spacing DT <double>");
      ierr = -1;
    }
    if (!b_n) {
      WUI(WARN,"Error: TRI_STRAND missing normal count N <int>");
      ierr = -1;
    }
    if (!b_zone_indices) {
      WUI(WARN,"Error: TRI_STRAND missing ZONES <comma-delimited-part:zone_names>");
      ierr = -1;
    }

    if (ierr != 0)
      return false;

    // TODO periodicity (can handle by not allowing periodic zone to set |= 2 flag)

    vector<GroupData> groupDataVec;
    buildGroupsForFlaggedZones(groupDataVec,zone_indices);
    if (groupDataVec.size() == 1) {
      const int ipart = groupDataVec[0].ipart;

      SurfaceShm * surface = partVec[ipart]->surface;
      assert(surface);

      int *sp_flag = new int[surface->nsp];
      for (int isp = 0; isp < surface->nsp; ++isp)
        sp_flag[isp] = 0;
      for (int ist = 0; ist < surface->nst; ++ist) {
        int igr;
        if (partVec[ipart]->getGroupForSt(igr,ist)) {
          assert(igr == 0);
          FOR_I3 sp_flag[surface->spost[ist][i]] |= 1;
        }
        else {
          FOR_I3 sp_flag[surface->spost[ist][i]] |= 2;
        }
      }
      vector<int> ispVec;
      for (int isp = 0; isp < surface->nsp; ++isp) {
        // don't seed points from group boundary
        if (sp_flag[isp] == 1) {
          sp_flag[isp] = ispVec.size();
          ispVec.push_back(isp);
        }
        else if (sp_flag[isp] == 3) {
          sp_flag[isp] = ispVec.size();
          ispVec.push_back(-isp-1);
        }
        else {
          sp_flag[isp] = -1;
        }
      }

      if (mpi_rank == 0)
        cout << " > number of original strands: " << ispVec.size() << endl;

      int ngp = ispVec.size();

      double* dist_gp = new double[ngp];
      {

        int * stogp_i = new int[ngp+1];
        for (int igp = 0; igp < ngp; ++igp) stogp_i[igp+1] = 0;
        for (int ist = 0; ist < surface->nst; ++ist) {
          int igr;
          if (partVec[ipart]->getGroupForSt(igr,ist)) {
            assert(igr == 0);
            FOR_I3 stogp_i[sp_flag[surface->spost[ist][i]]+1]++;
          }
        }
        stogp_i[0] = 0;
        for (int igp = 0; igp < ngp; ++igp) stogp_i[igp+1] += stogp_i[igp];
        int *stogp_v = new int[stogp_i[ngp]];
        for (int ist = 0; ist < surface->nst; ++ist) {
          int igr;
          if (partVec[ipart]->getGroupForSt(igr,ist)) {
            assert(igr == 0);
            FOR_I3 stogp_v[stogp_i[sp_flag[surface->spost[ist][i]]]++] = ist;
          }
        }
        for (int igp = ngp-1; igp > 0; --igp) stogp_i[igp] = stogp_i[igp-1];
        stogp_i[0] = 0;
        int* flag_gp = new int[ngp];

        MinHeap trialHeap(ngp);
        for (int igp = 0; igp < ngp; ++igp) {
          if (ispVec[igp] < 0) {
            dist_gp[igp] = 0.0;
            flag_gp[igp] = 1;
            trialHeap.insert(pair<double,int>(0.0,igp));
          }
          else {
            flag_gp[igp] = 0;
            dist_gp[igp] = HUGE_VAL;
          }
        }
        while (trialHeap.size() > 0)
          fmmLoop(trialHeap,flag_gp,dist_gp,stogp_i,stogp_v,sp_flag,surface);

        delete[] stogp_i;
        delete[] stogp_v;
        delete[] flag_gp;
      }

      vector<GroupedPoint> gpVec(ngp);
      for (int igp = 0; igp < ngp; ++igp) {
        const int isp = max(ispVec[igp],-ispVec[igp]-1);
        FOR_I3 gpVec[igp].x[i] = surface->xsp[isp][i];
        FOR_I3 gpVec[igp].n[i] = 0.0;
        gpVec[igp].dist = dist_gp[igp];
        gpVec[igp].kmin = 0;
      }
      dumpRange(dist_gp,ngp,"dist_gp");
      delete[] dist_gp;
      for (int ist = 0; ist < surface->nst; ++ist) {
        int igr;
        if (partVec[ipart]->getGroupForSt(igr,ist)) {
          assert(igr == 0);
          const double * const x0 = surface->xsp[surface->spost[ist][0]];
          const double * const x1 = surface->xsp[surface->spost[ist][1]];
          const double * const x2 = surface->xsp[surface->spost[ist][2]];
          const double n_st[3] = TRI_NORMAL_2(x0,x1,x2);
          FOR_I3 {
            const int igp = sp_flag[surface->spost[ist][i]];
            assert((igp >= 0)&&(igp < ngp));
            FOR_I3 gpVec[igp].n[i] -= n_st[i]; // flip to point in strand direction
          }
        }
      }
      for (int igp = 0; igp < ngp; ++igp) {
        const double n_mag = MAG(gpVec[igp].n); assert(n_mag > 0.0);
        FOR_I3 gpVec[igp].n[i] /= n_mag;
      }

      // TODO manage this array yourself
      vector<set<int> > nbogp(ngp);
      for (int ist = 0; ist < surface->nst; ++ist) {
        int igr;
        if (partVec[ipart]->getGroupForSt(igr,ist)) {
          assert(igr == 0);
          FOR_I3 {
            const int igp0 = sp_flag[surface->spost[ist][i]];
            const int igp1 = sp_flag[surface->spost[ist][(i+1)%3]];
            assert((igp0 >= 0)&&(igp0 < ngp));
            assert((igp1 >= 0)&&(igp1 < ngp));
            // only do 0 to 1. nbr tri will handle 1 to 0 (0 to 1 in their frame)...
            nbogp[igp0].insert(igp1);
          }
        }
      }

      const double sf = pow(dt/dn,1.0/double(n-1));
      const double delta_nom = (dt*sf-dn)/(sf-1.0);
      if (mpi_rank == 0)
        cout << " > delta: " << delta_nom << endl;
      for (int igp = 0; igp < ngp; ++igp) {
        if (ispVec[igp] >= 0) {
          gpVec[igp].kmax = 0;
          double delta = delta_nom;
          const double dn_current = min(dt,max(dn,dt+(dn-dt)/delta_nom*min(gpVec[igp].dist,delta_nom)));
          const double sf_current = (delta_nom-dn_current)/(delta_nom-dt);
          while (delta > 0.0) {
            double dx;
            if (dn_current == dt)
              dx = dn_current*gpVec[igp].kmax;
            else
              dx = dn_current*(pow(sf_current,gpVec[igp].kmax)-1.0)/(sf_current-1.0);
            delta -= dx;
            ++gpVec[igp].kmax;
          }
        }
      }
      ispVec.clear();

      if (b_branch||b_cull) {
        // create branch strands when the tangential length scale exceeds factor*dt...
        int igp = 0;
        while (igp < ngp) {
          if (gpVec[igp].kmax == -1) {
            ++igp;
            continue;
          }
          set<int>::iterator it = nbogp[igp].begin();
          while (it != nbogp[igp].end()) {
            const int inb = *it;
            assert(gpVec[inb].kmin >= 0);
            if (gpVec[inb].kmax == -1) {
              ++it;
              continue;
            }

            const double dx_[3] = DIFF(gpVec[inb].x,gpVec[igp].x);
            const double dx_mag = MAG(dx_); assert(dx_mag > 0.0);

            double n_br[3]; FOR_I3 n_br[i] = gpVec[igp].n[i]+gpVec[inb].n[i];
            double mag = MAG(n_br); assert(mag > 0.0);
            FOR_I3 n_br[i] /= mag;

            double n_perp[3] = CROSS_PRODUCT(n_br,dx_);
            mag = MAG(n_perp); assert(mag > 0.0);
            FOR_I3 n_perp[i] /= mag;

            double tmp = DOT_PRODUCT(gpVec[igp].n,n_perp);
            double np_gp[3]; FOR_I3 np_gp[i] = gpVec[igp].n[i]-tmp*n_perp[i];
            mag = MAG(np_gp); assert(mag > 0.0);
            FOR_I3 np_gp[i] /= mag;

            tmp = DOT_PRODUCT(gpVec[inb].n,n_perp);
            double np_nb[3]; FOR_I3 np_nb[i] = gpVec[inb].n[i]-tmp*n_perp[i];
            mag = MAG(np_nb); assert(mag > 0.0);
            FOR_I3 np_nb[i] /= mag;

            const double alpha = DOT_PRODUCT(np_nb,dx_);
            const double beta = DOT_PRODUCT(np_gp,dx_);
            const double gamma = DOT_PRODUCT(np_nb,np_gp);
            const double v = (gamma*alpha-beta)/(1.0-gamma*gamma);
            const double u = alpha + gamma*v;
            double x0[3]; FOR_I3 x0[i] = gpVec[inb].x[i]-np_nb[i]*u;
            double x_br[3]; FOR_I3 x_br[i] = 0.5*(gpVec[igp].x[i]+gpVec[inb].x[i]);
            double delta = delta_nom;

            if (b_branch&&(u > 0.0)&&(v > 0.0)) {

              double delta_split;
              if (dx_mag > factor*dt) {
                delta_split = 0.0;
              }
              else {
                double w = 0; FOR_I3 w += (x_br[i]-x0[i])*n_br[i];
                delta_split = w*(2.0*factor*dt/dx_mag-1.0);
              }
              if ((delta_split >= 0.0)&&(delta_split < delta)) {
                gpVec.push_back(GroupedPoint());
                FOR_I3 gpVec[ngp].n[i] = n_br[i];
                FOR_I3 gpVec[ngp].x[i] = x_br[i];
                gpVec[ngp].dist = 0.5*(gpVec[igp].dist+gpVec[inb].dist);
                double tmp = delta_split;
                gpVec[ngp].kmin = 0;
                const double dn_current = min(dt,max(dn,dt+(dn-dt)/delta_nom*min(gpVec[ngp].dist,delta_nom)));
                const double sf_current = (delta_nom-dn_current)/(delta_nom-dt);
                while (tmp > 0.0) {
                  double dx;
                  if (dn_current == dt)
                    dx = dn_current*gpVec[ngp].kmin;
                  else
                    dx = dn_current*(pow(sf_current,gpVec[ngp].kmin)-1.0)/(sf_current-1.0);
                  tmp -= dx;
                  ++gpVec[ngp].kmin;
                }
                gpVec[ngp].kmax = gpVec[ngp].kmin;
                while (delta > delta_split) {
                  double dx;
                  if (dn_current == dt)
                    dx = dn_current*gpVec[ngp].kmax;
                  else
                    dx = dn_current*(pow(sf_current,gpVec[ngp].kmax)-1.0)/(sf_current-1.0);
                  delta -= dx;
                  ++gpVec[ngp].kmax;
                }
                assert(gpVec[ngp].kmax > gpVec[ngp].kmin);

                // split edge...

                {
                  set<int>::iterator it2 = nbogp[inb].find(igp);
                  assert(it2 != nbogp[inb].end());
                  nbogp[inb].erase(it2);
                }
                nbogp[igp].erase(it++);

                nbogp.push_back(set<int>());
                nbogp[ngp].insert(inb);
                nbogp[ngp].insert(igp);
                int count = 0;
                for (set<int>::iterator it2 = nbogp[igp].begin(); it2 != nbogp[igp].end(); ++it2) {
                  set<int>::iterator it3 = nbogp[inb].find(*it2);
                  if (it3 != nbogp[inb].end()) {
                    const int inb2 = *it3;
                    nbogp[ngp].insert(inb2);
                    nbogp[inb2].insert(ngp);
                    ++count;
                    if (count == 2)
                      break;
                  }
                }
                assert(count == 2);
                nbogp[igp].insert(ngp);
                nbogp[inb].insert(ngp);
                // TODO should manage memory myself so that I only call this during a rellocation
                it = nbogp[igp].begin();

                ++ngp;
              }
              else {
                ++it;
              }

            }
            else {

              if (b_cull&&(u < 0.0)&&(v < 0.0)) {
                double w = 0; FOR_I3 w -= (x_br[i]-x0[i])*n_br[i]; // flipped sign
                if ((w > 0.0)&&(w < delta)) {
                  // change kmax...
                  int kmax = 0;
                  const double dn_current = min(dt,max(dn,dt+(dn-dt)/delta_nom*min(gpVec[igp].dist,delta_nom)));
                  const double sf_current = (delta_nom-dn_current)/(delta_nom-dt);
                  while (w > 0.0) {
                    double dx;
                    if (dn_current == dt)
                      dx = dn_current*kmax;
                    else
                      dx = dn_current*(pow(sf_current,kmax)-1.0)/(sf_current-1.0);
                    w -= dx;
                    ++kmax;
                  }
                  gpVec[igp].kmax = min(gpVec[igp].kmax,kmax);
                  gpVec[igp].kmin = min(gpVec[igp].kmin,kmax);
                }
              }

              ++it;
            }

          }
          ++igp; // go to next grouped point
        }
      }
      assert(gpVec.size() == ngp);

      if (mpi_rank == 0)
        cout << " > number of final strands: " << ngp << endl;

      int my_ngp_avg = ngp/mpi_size;
      if (ngp%mpi_size) ++my_ngp_avg;
      const int igp0 = min(ngp,mpi_rank*my_ngp_avg);
      const int igp1 = min(ngp,(mpi_rank+1)*my_ngp_avg); // non-inclusive
      assert(igp1-igp0 <= my_ngp_avg);
      //const int my_ngp = igp1-igp0;

      int np = -1;
      for (int iter = 0; iter < 2; ++iter) {
        const int np_check = np;
        np = 0;
        for (int igp = igp0; igp < igp1; ++igp) {
          if (gpVec[igp].kmax == -1)
            continue;
          if (iter == 0) {
            assert((gpVec[igp].kmin >= 0)&&(gpVec[igp].kmax >= 0));
            np += (gpVec[igp].kmax-gpVec[igp].kmin+1);
          }
          else {
            const double dn_current = min(dt,max(dn,dt+(dn-dt)/delta_nom*min(gpVec[igp].dist,delta_nom)));
            const double sf_current = (delta_nom-dn_current)/(delta_nom-dt);
            for (int k = gpVec[igp].kmin; k <= gpVec[igp].kmax; ++k) {
              double dx;
              if (dn_current == dt)
                dx = dn_current*0.5*(2*k+1);
              else
                dx = dn_current*0.5*( (pow(sf_current,k)-1.0)/(sf_current-1.0) + (pow(sf_current,k+1)-1.0)/(sf_current-1.0));
              FOR_I3 part->pts->xp[np][i] = gpVec[igp].x[i] + dx*gpVec[igp].n[i];
              part->pts->delta[np] = 1.2*sqrt(3.0*dt*dt); // use largest scale
              if (k < gpVec[igp].kmax)
                part->pts->flag[np] = 0;
              else
                part->pts->flag[np] = 1;
              ++np;
            }
          }
        }
        if (iter == 0) {
          // pts...
          assert(part->pts == NULL);
          part->pts = new Points();
          part->pts->init(np); // this call takes np (local), does mem allocation, and reduces np_global...
          assert(part->pts->np_global > 0);
          assert(part->pts->flag == NULL);
          part->pts->flag = new int[np]; // added flag points closest to ff_surface
        }
        else {
          assert(np == np_check);
        }
      }
      assert(np >= 0);

      // build farfield surface from nbogp...

      set<Triple<int,int,int> > triSet;
      if (mpi_rank_shared == 0) {
        for (int igp = 0; igp < ngp; ++igp) {
          for (set<int>::iterator it = nbogp[igp].begin(); it != nbogp[igp].end(); ++it) {
            const int inb = *it;
            for (set<int>::iterator it2 = nbogp[igp].begin(); it2 != nbogp[igp].end(); ++it2) {
              set<int>::iterator it3 = nbogp[inb].find(*it2);
              if (it3 != nbogp[inb].end()) {
                const int inb2 = *it3;
                vector<int> tmp(3);
                tmp[0] = inb2;
                tmp[1] = inb;
                tmp[2] = igp;
                sort(tmp.begin(),tmp.end());
                triSet.insert(Triple<int,int,int>(tmp[0],tmp[1],tmp[2]));
              }
            }
          }
        }
      }
      nbogp.clear();
      const int ff_nsp = ngp;
      int ff_nst = triSet.size();
      MPI_Bcast(&ff_nst,1,MPI_INT,0,mpi_comm_shared);

      if (mpi_rank == 0)
        cout << " > farfield surface nsp,nst: " << ff_nsp << " " << ff_nst << endl;

      assert(part->ff_surface == NULL);
      part->ff_surface = new SurfaceShm();
      part->ff_surface->init(ff_nsp,ff_nst);
      assert(part->ff_surface_dxost == NULL);
      CTI_Mmap(part->ff_surface_dxost,part->ff_surface->nst);
      assert(part->ff_type == UNKNOWN_FF);
      part->ff_type = FF_SURFACE_FAZONE_FF;
      assert(part->solid_type == UNKNOWN_SOLID);
      part->solid_type = NO_SOLID;

      if (mpi_rank_shared == 0) {
        for (int isp = 0; isp < ff_nsp; ++isp) {
          const double dn_current = min(dt,max(dn,dt+(dn-dt)/delta_nom*min(gpVec[isp].dist,delta_nom)));
          const double sf_current = (delta_nom-dn_current)/(delta_nom-dt);
          double dx;
          if (dn_current == dt)
            dx = dn_current*(gpVec[isp].kmax+1);
          else
            dx = dn_current*(pow(sf_current,gpVec[isp].kmax+1)-1.0)/(sf_current-1.0);
          FOR_I3 part->ff_surface->xsp[isp][i] = gpVec[isp].x[i] + dx*gpVec[isp].n[i];
        }
        int ist = 0;
        for (set<Triple<int,int,int> >::iterator it = triSet.begin(); it != triSet.end(); ++it) {
          part->ff_surface_dxost[ist] = dt;
          part->ff_surface->znost[ist] = 0; // there is no surface
          const double n[3] = TRI_NORMAL_2(part->ff_surface->xsp[it->first],
					   part->ff_surface->xsp[it->second],part->ff_surface->xsp[it->third]);
          // use average strand normal to check sign
          double n_avg[3];
          FOR_I3 n_avg[i] = gpVec[it->first].n[i];
          FOR_I3 n_avg[i] += gpVec[it->second].n[i];
          FOR_I3 n_avg[i] += gpVec[it->third].n[i];
          if (DOT_PRODUCT(n,n_avg) > 0.0) {
            part->ff_surface->spost[ist][0] = it->first;
            part->ff_surface->spost[ist][1] = it->second;
            part->ff_surface->spost[ist][2] = it->third;
          }
          else {
            part->ff_surface->spost[ist][0] = it->first;
            part->ff_surface->spost[ist][1] = it->third;
            part->ff_surface->spost[ist][2] = it->second;
          }
          ++ist;
        }
        assert(ist == ff_nst);
      }
      gpVec.clear();
      triSet.clear();
      MPI_Barrier(mpi_comm_shared);

      // all grouped st's get their part index set to partVec.size() (note that
      // the part passed into this routine will be pushed onto partVec after we
      // return succesfully)...

      if (mpi_rank_shared == 0) {
        for (int ipart = 0, limit = partVec.size(); ipart < limit; ++ipart) {
          if (partVec[ipart]->surface) {
            SurfaceShm * surface = partVec[ipart]->surface;
            for (int ist = 0; ist < surface->nst; ++ist) {
              int igr;
              if (partVec[ipart]->getGroupForSt(igr,ist)) {
                partVec[ipart]->clearGroupForSt(ist);
                partVec[ipart]->setPartForSt(part->ipart_of_part,ist);
              }
            }
          }
        }
      }
      MPI_Barrier(mpi_comm_shared);

      // and clear the groupDataVec: we don't use the same clear function as in
      // the queryZones routine becasue we have already cleared the ist group above

      for (int igr = 0; igr < groupDataVec.size(); ++igr) {
        groupDataVec[igr].clear();
      }
      groupDataVec.clear();

      //GeomUtils::writePtsTecplot("tri_strand_pts.dat",part->pts->xp,part->pts->np);
      //if (mpi_rank == 0) GeomUtils::writeSbin("ff.sbin",part->ff_surface->spost,part->ff_surface->nst,part->ff_surface->xsp);
      //stringstream temp; temp << "ff" << partVec.size() << ".sbin";
      //string filename = temp.str();
      //if (mpi_rank == 0) GeomUtils::writeSbin(filename,part->ff_surface->spost,part->ff_surface->nst,part->ff_surface->xsp);

      // and return...

      assert(ierr == 0);
      return true;

    }
    else {
      WUI(WARN," > TRI_STRAND currently supports single disjoint group. Returning...");

      if (mpi_rank_shared == 0) {
        for (int ipart = 0, limit = partVec.size(); ipart < limit; ++ipart) {
          if (partVec[ipart]->surface) {
            SurfaceShm * surface = partVec[ipart]->surface;
            for (int ist = 0; ist < surface->nst; ++ist) {
              int igr;
              if (partVec[ipart]->getGroupForSt(igr,ist)) {
                partVec[ipart]->clearGroupForSt(ist);
              }
            }
          }
        }
      }
      MPI_Barrier(mpi_comm_shared);
      for (int igr = 0; igr < groupDataVec.size(); ++igr) {
        groupDataVec[igr].clear();
      }
      groupDataVec.clear();

      return false;

    }

  }

  bool makeStrand(Part * part,Param * param,int& iarg) {

    bool b_zone_indices = false;
    set<int> zone_indices;
    bool b_dn = false;
    double dn;
    bool b_dt = false;
    double dt;
    bool b_dt1 = false;
    double dt1;
    bool b_dt2 = false;
    double dt2;
    bool b_n = false;
    int n;
    bool b_crease = false;
    double s_crease;

    double * n_plane = NULL;
    //bool b_nt1 = false;
    //int nt1;
    //bool b_nt2 = false;
    //int nt2;
    int mode = 0; // changes the stretching behavior in the tangential plane
    int smooth_iters = 50;

    part->ff_type = FF_SURFACE_FAZONE_FF;
    part->solid_type = NO_SOLID;

    int ierr = 0;
    while (iarg < param->size()) {
      string token = MiscUtils::toUpperCase(param->getString(iarg++));
      if (token == "ZONEID") {
        // expect a comma-delimited set of integers...
        b_zone_indices = true;
        const string zonesCsv = param->getString(iarg++);
        vector<int> tmp;
        if (MiscUtils::splitCsv(tmp,zonesCsv)) {
          for (int ii = 0; ii < tmp.size(); ++ii) {
            const int izone = tmp[ii];
            if ((izone >= 0)&&(izone < zoneVec.size())) {
              zone_indices.insert(izone);
            }
            else {
              WUI(WARN,"STRAND ZONEID out of range: " << izone);
              ierr = -1;
            }
          }
        }
        else {
          WUI(WARN,"STRAND ZONEID expecting comma-delimited ints: " << zonesCsv);
          ierr = -1;
        }
      }
      else if ((token == "ZONE")||(token == "ZONES")||(token == "ZONE_NAME")||(token == "ZONE_NAMES")) {
        // expect a comma-delimited list of names...
        b_zone_indices = true;
        const string zoneNamesCsv = param->getString(iarg++);
        vector<string> tmp;
        MiscUtils::splitCsv(tmp,zoneNamesCsv);
        for (int ii = 0; ii < tmp.size(); ++ii) {
          const int izone = getZoneIndex(tmp[ii]); // get the unique zoneIndex in the partVec
          if (izone >= 0) {
            zone_indices.insert(izone);
          }
          else {
            WUI(WARN,"STRAND ZONE not recognized: " << tmp[ii]);
            ierr = -1;
          }
        }
      }
      else if (token == "DN") {
        b_dn = true;
        dn = param->getDouble(iarg++);
      }
      else if (token == "DT") {
        b_dt = true;
        dt = param->getDouble(iarg++);
      }
      else if (token == "DT1") {
        b_dt1 = true;
        dt1 = param->getDouble(iarg++);
      }
      /*
        else if (token == "NT1") {
        b_nt1 = true;
        nt1 = param->getInt(iarg++);
        }
      */
      else if (token == "DT2") {
        b_dt2 = true;
        dt2 = param->getDouble(iarg++);
      }
      /*
        else if (token == "NT2") {
        b_nt2 = true;
        nt2 = param->getInt(iarg++);
        }
      */
      else if (token == "N") {
        b_n = true;
        n = param->getInt(iarg++);
      }
      else if (token == "N_PLANE") {
        n_plane = new double[3];
        FOR_I3 n_plane[i] = param->getDouble(iarg++);
      }
      else if (token == "MODE") {
        mode = param->getInt(iarg++);
      }
      else if (token == "NLAYERS") {
        part->ff_nlayersString = param->getString(iarg++);
      }
      else if (token == "SMOOTH_ITERS") {
        smooth_iters = param->getInt(iarg++);
      }
      else if (token == "SPLINE_CREASE_TOL") {
        b_crease = true;
        s_crease = param->getDouble(iarg++);
      }
      else {
        WUI(WARN,"Unknown STRAND parameter " << token);
        ierr = -1;
      }
    }

    // TODO: nt2

    if (!b_dn) {
      WUI(WARN,"Error: STRAND missing normal spacing DN <double>");
      ierr = -1;
    }
    if (!b_dt) {
      if (b_dt1) {
        b_dt = true;
        dt = dt1;
      }
      else {
        WUI(WARN,"Error: STRAND missing tangential spacing DT <double>");
        ierr = -1;
      }
    }
    else {
      if (!b_dt1) {
        b_dt1 = true;
        dt1 = dt;
      }
      if (!b_dt2) {
        b_dt2 = true;
        dt2 = dt;
      }
    }

    if (!b_n) {
      WUI(WARN,"Error: STRAND missing normal count N <int>");;
      ierr = -1;
    }
    if (!b_zone_indices) {
      WUI(WARN,"Error: STRAND missing ZONES <comma-delimited-indices> or ZONE_NAMES <comma-delimited-part:zone_names>");;
      ierr = -1;
    }

    if (ierr != 0)
      return false;

    vector<GroupData> groupDataVec;
    buildGroupsForFlaggedZones(groupDataVec,zone_indices);

    WUI(INFO,"Specified zone(s) produced " << groupDataVec.size() << " disjoint group(s)...");

    // loop 1: count ff surface tris and points for each valid group...
    int counts[5] = { 0, 0, 0, 0, 0 }; // surface->nsp,nst(both global),ff_surface->nsp,nst(both global),pts->np(local)
    for (int igr = 0; igr < groupDataVec.size(); ++igr) {
      if (b_crease) groupDataVec[igr].setSplineCreaseTol(s_crease);

      if (groupDataVec[igr].isCylinder(dn*0.25)) {
        if (mpi_rank == 0)
          cout << "group " << igr << ": looks like a cylinder with length: " <<
            groupDataVec[igr].getCylinderL() << " and radius: " << groupDataVec[igr].getCylinderR() << endl;
        // right cylinders only...
        if (groupDataVec[igr].isRightCylinder()) {
          if (mpi_rank == 0)
            cout << "... and its a right cylinder!" << endl;
          int my_counts[5];
          groupDataVec[igr].countRightCylinderBl(my_counts,dn,dt,n);
          FOR_I5 counts[i] += my_counts[i];
        }
        else {
          // this is a non-right cylinder, so we cannot do a strand mesh on ANY group...
          ierr = -1;
        }
      }
      else if (groupDataVec[igr].isLogicalQuad()) {
        if (mpi_rank == 0) cout << "group " << igr << ": looks like logical quad" << endl;
        assert(b_dt1 && b_dt2); // TODO: handle the DT case by setting DT1 and DT2 to DT...
        int my_counts[5];
        groupDataVec[igr].countLogicalQuadBl(my_counts,dn,dt1,dt2,n,mode);
        if (mpi_rank == 0) cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXX my_counts[2]: " << my_counts[2] << endl;
        FOR_I5 counts[i] += my_counts[i];
      }
      else {
        // once we have tested, switch to unknown so we don't test again...
        groupDataVec[igr].setKind(UNKNOWN_GROUP_KIND);
        if (mpi_rank == 0) cout << "group: " << igr << " looks like nothing" << endl;
        ierr = -1;
      }
    }

    // only make strand meshes if ALL groups can be made into strand meshes...
    if (ierr != 0) {
      clearGroups(groupDataVec);
      WUI(WARN,"one or more groups are not strandable");
      return false;
    }

    // record the associated parts
    part->b_link_surface = true;

    assert(part->ipart_link.empty());
    for (int igr = 0; igr < groupDataVec.size(); ++igr) {
      const int ipart = groupDataVec[igr].ipart;
      part->ipart_link.push_back(ipart);
    }

    assert(part->ist_offset_link.empty());
    part->ist_offset_link.push_back(0);
    for (int igr = 0; igr < groupDataVec.size()-1; ++igr) {
      const int ipart = groupDataVec[igr].ipart;
      assert(partVec[ipart]->surface);
      const int back = part->ist_offset_link.back();
      part->ist_offset_link.push_back(back+partVec[ipart]->surface->nst);
    }

    if (mpi_rank == 0) cout << "XXXXXXXXXXXXXXXX counts[2]: " << counts[2] << endl;

    // surface: should not be present...
    assert(part->surface == NULL);
    assert(counts[0] == 0);
    assert(counts[1] == 0);
    // ff_surface...
    assert(part->ff_surface == NULL);
    assert(counts[2] > 0);
    assert(counts[3] > 0);
    part->ff_surface = new SurfaceShm();
    part->ff_surface->init(counts[2],counts[3]);
    assert(part->ff_surface_dxost == NULL);
    CTI_Mmap(part->ff_surface_dxost,counts[3]);
    // pts...
    assert(part->pts == NULL);
    part->pts = new Points();
    part->pts->init(counts[4]); // this call takes np (local), does mem allocation, and reduces np_global...
    assert(part->pts->np_global > 0);
    if (mpi_rank == 0) {
      cout << " > overall surface->nsp,nst: " << counts[0] << " " << counts[1] <<
        ", ff_surface->nsp,nst: " << counts[2] << " " << counts[3] <<
        ", pts->np_global: " << part->pts->np_global << endl;
    }

    // now build...
    int counts_check[5]; FOR_I5 counts_check[i] = counts[i];
    FOR_I5 counts[i] = 0;
    for (int igr = 0; igr < groupDataVec.size(); ++igr) {
      if (groupDataVec[igr].isCylinder(dn*0.25)) {
        assert(groupDataVec[igr].isRightCylinder());
        groupDataVec[igr].makeRightCylinderBl(part,counts,dn,dt,n);
      }
      else if (groupDataVec[igr].isLogicalQuad()) {
        groupDataVec[igr].makeLogicalQuadBl(part,counts,dn,dt1,dt2,n,mode,smooth_iters,n_plane);
      }
      else {
        assert(0);
      }
    }
    if (mpi_rank_shared == 0) {
      // these counts only need to match on shared 0...
      assert(counts[0] == counts_check[0]);
      assert(counts[1] == counts_check[1]);
      if (!(counts[2] == counts_check[2]))
        cout << "XXXXXXXX counts[2]: " << counts[2] << " counts_check[2]: " << counts_check[2] << endl;
      assert(counts[2] == counts_check[2]);
      if (!(counts[3] == counts_check[3]))
        cout << "XXXXXXXX counts[3]: " << counts[3] << " counts_check[3]: " << counts_check[3] << endl;
      assert(counts[3] == counts_check[3]);
    }
    // points can be distributed across all processors...
    assert(counts[4] == counts_check[4]);

    // all grouped st's get their part index set to partVec.size() (note that
    // the part passed into this routine will be pushed onto partVec after we
    // return succesfully)...

    if (mpi_rank_shared == 0) {
      for (int ipart = 0, limit = partVec.size(); ipart < limit; ++ipart) {
        if (partVec[ipart]->surface) {
          SurfaceShm * surface = partVec[ipart]->surface;
          for (int ist = 0; ist < surface->nst; ++ist) {
            int igr;
            if (partVec[ipart]->getGroupForSt(igr,ist)) {
              partVec[ipart]->clearGroupForSt(ist);
              partVec[ipart]->setPartForSt(part->ipart_of_part,ist);
            }
          }
        }
      }
    }
    MPI_Barrier(mpi_comm_shared);

    DELETE(n_plane);

    // and clear the groupDataVec: we don't use the same clear function as in
    // the queryZones routine becasue we have already cleared the ist group above

    for (int igr = 0; igr < groupDataVec.size(); ++igr) {
      groupDataVec[igr].clear();
    }
    groupDataVec.clear();

    // and return...

    assert(ierr == 0);
    return true;

  }

  bool makeAirfoilStrand(Part * part,Param * param,int& iarg) {

    set<int> zone_indices;
    bool b_zone_indices = false;
    bool b_name = false;
    string name;
    int ierr = 0;
    bool b_axis = false;
    double axis[3];
    bool b_nr = false;
    int nr; // number of sampling points in spanwise direction, input for spline
    bool b_dn = false;
    double dn; // normal
    bool b_dt = false;
    double dt; // tangential
    bool b_ds = false;
    double ds; // span (== tangential if not specified separately)
    bool b_nn;
    int nn; // normal count: how many spaces to go from dn to dt

    //Part * mypart = new Part();

    while (iarg < param->size()) {
      const string token = MiscUtils::toUpperCase(param->getString(iarg++));
      if (token == "NAME") {
        b_name = true;
        name = param->getString(iarg++);
      }
      /*
	else if ((token == "SBIN")||(token == "BIN")) {
        const string filename = param->getString(iarg++);
        const int this_ierr = makePartFromSbin(mypart,filename);
        if (this_ierr != 0) {
	WUI(WARN,"problem creating AIRFOIL_STRAND SBIN for sbin file: " << token);
	ierr = -1;
        }
	}
      */
      else if ((token == "ZONE")||(token == "ZONES")||(token == "ZONE_NAME")||(token == "ZONE_NAMES")) {
        // expect a comma-delimited list of names...
        b_zone_indices = true;
        const string zoneNamesCsv = param->getString(iarg++);
        vector<string> tmp;
        MiscUtils::splitCsv(tmp,zoneNamesCsv);
        for (int ii = 0; ii < tmp.size(); ++ii) {
          const int izone = getZoneIndex(tmp[ii]); // get the unique zoneIndex in the partVec
          if (izone >= 0) {
            zone_indices.insert(izone);
          }
          else {
            WUI(WARN,"AIRFOIL_STRAND ZONE not recognized: " << tmp[ii]);
            ierr = -1;
          }
        }
      }
      else if (token == "SUBZONES") {
        CERR("token SUBZONES not supported yet.");
      }
      else if (token == "AXIS") {
	b_axis = true;
        FOR_I3 axis[i] = param->getDouble(iarg++);
        const double axis_mag = MAG(axis);
        FOR_I3 axis[i] = axis[i]/axis_mag;
      }
      else if (token == "NR") {
	b_nr = true;
        nr = param->getInt(iarg++); 
      }
      else if (token == "DN") {
        b_dn = true;
        dn = param->getDouble(iarg++);
      }
      else if (token == "NN") {
        b_nn = true;
        nn = param->getInt(iarg++);
        assert(nn > 2);
      }
      else if (token == "DT") {
        b_dt = true;
        dt = param->getDouble(iarg++);
      }
      else if (token == "DS") {
        b_ds = true;
        ds = param->getDouble(iarg++);
      }
      else if (token == "NLAYERS") {
	part->ff_nlayersString = param->getString(iarg++);
      }
    }

    if (!b_axis) {
      if (mpi_rank == 0) cout << "AIRFOIL_STRAND: missing AXIS <double> <double> <double>" << endl;
      ierr = -1;
    }
    if (!b_dn) {
      if (mpi_rank == 0) cout << "AIRFOIL_STRAND: missing DN <double>" << endl;
      ierr = -1;
    }
    if (!b_nn) {
      if (mpi_rank == 0) cout << "AIRFOIL_STRAND: missing NN <int>" << endl;
      ierr = -1;
    }
    if (!b_nr) {
      if (mpi_rank == 0) cout << "AIRFOIL_STRAND: missing NR <int>" << endl;
      ierr = -1;
    }
    if (!b_dt) {
      if (mpi_rank == 0) cout << "AIRFOIL_STRAND: missing DT <double>" << endl;
      ierr = -1;
    }
    if (ierr != 0) {
      return false;
    }

    // default for ds is dt if not specified separately...
    if (!b_ds) {
      ds = dt;
    }

    vector<GroupData> groupDataVec;
    buildGroupsForFlaggedZones(groupDataVec,zone_indices);
    if (groupDataVec.size() == 1) {

      const int ipart = groupDataVec[0].ipart;
      const int izone = *(zone_indices.begin());

      part->b_link_surface = true;
      assert(part->ipart_link.empty());
      part->ipart_link.push_back(ipart);
      assert(part->ist_offset_link.empty());
      part->ist_offset_link.push_back(0);

      SurfaceShm * surface = partVec[ipart]->surface;
      assert(surface);
    
      // build flag for active surface tris
      int * st_flag = new int[surface->nst];
      for (int ist=0; ist<surface->nst; ++ist) {
	const int isz = surface->znost[ist];
	const int izn = partVec[ipart]->znosz[isz]; // convert from surface zone (isz) to global zone (izn)
	if (izn == izone)
	  st_flag[ist] = 1;
	else
	  st_flag[ist] = 0;
      }

      COUT1(" > about to start groupOpenEdges()");
      vector<vector<int> > groupSpVecs;   
      //ierr = surface->groupOpenEdges(groupSpVecs);
      if (mpi_rank == 0) cout << "found " << groupDataVec[0].nlink << " edge loops" << endl;
      for (int il = 0; il < groupDataVec[0].nlink; ++il) {
	vector<int> this_group;
	for (int sol = groupDataVec[0].spoli_i[il]; sol != groupDataVec[0].spoli_i[il+1]; ++sol) {
	  this_group.push_back(groupDataVec[0].spoli_v[sol]);
	}
	if (mpi_rank == 0) cout << "this loop has " << this_group.size() << " points" << endl;
	groupSpVecs.push_back(this_group);
      }

      const int group_index = groupSpVecs.size();
      if (group_index != 2) 
	CERR("There should be 2 and only 2 open edge loops! Check sbin!");
      //for (int igr = 0; igr<group_index; ++igr) {
      //  COUT1("\n > group " << igr);
      //  for (int isp = 0,nsp = groupSpVecs[igr].size(); isp<nsp; ++isp) {
      //    COUT1(COUT_VEC(mypart->surface->xsp[groupSpVecs[igr][isp]]));
      //   }
      //}

      // find leading edge and trailing edge pts based on max distance
      double LE_pts[group_index][3];
      double TE_pts[group_index][3];
      int LE_idx[group_index];
      double ref_vec[3];
      for (int igr = 0; igr<group_index; ++igr) {
	double max_dist = 0.0;
	int idx0, idx1;
	const int np = groupSpVecs[igr].size();
	for (int i=0; i<np; ++i) {
	  for (int j=i+1; j<=np; ++j) {
	    const double dist_ij = DIST(surface->xsp[groupSpVecs[igr][i]],surface->xsp[groupSpVecs[igr][j%np]]);
	    if (dist_ij > max_dist) {
	      max_dist = dist_ij;
	      idx0 = i;
	      idx1 = j;
	    }
	  }
	}
	if (mpi_rank == 0) cout << "igr, idx0, idx1, max_dist = " << igr << " " << idx0 << " " << idx1 << " " << max_dist << endl;
	double pt0[3], pt1[3];
	FOR_I3 pt0[i] = surface->xsp[groupSpVecs[igr][idx0]][i];
	FOR_I3 pt1[i] = surface->xsp[groupSpVecs[igr][idx1]][i];
	double vec01[3];
	FOR_I3 vec01[i] = pt1[i] - pt0[i];
	if (igr == 0) {
	  const double dp = DOT_PRODUCT(vec01, axis);
	  if (dp >= 0) {
	    FOR_I3 LE_pts[igr][i] = pt0[i];
	    FOR_I3 TE_pts[igr][i] = pt1[i];
	    LE_idx[igr] = idx0;
	  }
	  else {
	    FOR_I3 LE_pts[igr][i] = pt1[i];
	    FOR_I3 TE_pts[igr][i] = pt0[i];
	    LE_idx[igr] = idx1;
	  }
	  FOR_I3 ref_vec[i] = TE_pts[igr][i] - LE_pts[igr][i];
	}
	else {
	  const double dp = DOT_PRODUCT(vec01, ref_vec);
	  if (dp >= 0) {
	    FOR_I3 LE_pts[igr][i] = pt0[i];
	    FOR_I3 TE_pts[igr][i] = pt1[i];
	    LE_idx[igr] = idx0;
	  }
	  else {
	    FOR_I3 LE_pts[igr][i] = pt1[i];
	    FOR_I3 TE_pts[igr][i] = pt0[i];
	    LE_idx[igr] = idx1;
	  }
	}
      }
      // check
      if (mpi_rank == 0) {
	for (int igr=0; igr<group_index; ++igr){
	  cout << "igr, LE_pt: " << igr << " " << COUT_VEC(LE_pts[igr]) << endl;
	  cout << "igr, TE_pt: " << igr << " " << COUT_VEC(TE_pts[igr]) << endl;
	}
      }


      // compare the lengths
      // use the longer one as root, the other as tip
      double lengths[group_index];
      for (int igr = 0; igr<group_index; ++igr) {
	lengths[igr] = 0.0;
	const int np = groupSpVecs[igr].size();
	for (int i=0; i<np; ++i) {
	  lengths[igr] += DIST(surface->xsp[groupSpVecs[igr][i]],surface->xsp[groupSpVecs[igr][(i+1)%np]]);
	}
      }
      // check
      if (mpi_rank == 0) cout << "lengths of the two edge loops are " << lengths[0] << " " << lengths[1] << endl;

      int idx_root, idx_tip;
      if (lengths[0] > lengths[1]) {
	idx_root = 0;
	idx_tip  = 1;
      }
      else {
	idx_root = 1;
	idx_tip  = 0;
	// swap LE_pts and TE_pts
	double tmp[3];
	FOR_I3 tmp[i] = LE_pts[0][i];
	FOR_I3 LE_pts[0][i] = LE_pts[1][i];
	FOR_I3 LE_pts[1][i] = tmp[i];
	FOR_I3 tmp[i] = TE_pts[0][i];
	FOR_I3 TE_pts[0][i] = TE_pts[1][i];
	FOR_I3 TE_pts[1][i] = tmp[i];
      }
      // check
      if (mpi_rank == 0) cout << "idx_root, idx_tip = " << idx_root << " " << idx_tip << endl;

      /*
      // now calc ns based on LE length and ds
      const double LE_len = DIST(LE_pts[0], LE_pts[1]);
      const int ns = int(LE_len/ds) + 1;
      cout << "LE_len, ns = " << LE_len << " " << ns << endl;
      */

      // cluster points near the root and the tip. calc s_local
      const double LE_len = DIST(LE_pts[0], LE_pts[1]);
      const double delta = GeomUtils::stretchingFunction(nn,nn,dn,dt);
      const double max_incre = delta - GeomUtils::stretchingFunction(nn-1,nn,dn,dt);
      const int ns_mid = int((LE_len-2.0*delta)/max_incre);
      const double ds_mid = (LE_len-2.0*delta)/double(ns_mid);
      // check
      if (mpi_rank == 0) cout << "ds_mid/max_incre = " << ds_mid/max_incre << endl;

      const int ns = 2*(nn+1) + ns_mid - 1;

      double *s_local = new double[ns];
      for (int is=0; is<ns; ++is) {
	if (is<=nn) {
	  s_local[is] = GeomUtils::stretchingFunction(is,nn,dn,dt)/LE_len;
	}
	else if (is>=nn+ns_mid) {
	  s_local[is] = (LE_len - GeomUtils::stretchingFunction(ns-1-is,nn,dn,dt))/LE_len;
	}
	else {
	  s_local[is] = (delta + double(is-nn)*ds_mid)/LE_len;
	}
      }

      s_local[0] = 0.0;
      s_local[ns-1] = 1.0;

      // smooth s_local
      const int nsmooth1 = 5;
      double * s_copy = new double[ns];
      for (int ismooth = 0; ismooth < nsmooth1; ++ismooth) {
	for (int is = 0; is < ns; ++is)
	  s_copy[is] = s_local[is];
	for (int is = 1; is < ns-1; ++is) {
	  s_local[is] = 0.5*(s_copy[is-1]+s_copy[is+1]);
	}
      }
      delete[] s_copy;
      // check
      //for (int is = 0; is < ns; ++is) {
      //  cout << "s_local[" << is << "] = " << s_local[is] << endl;
      //}

      // find a vector pointing from the root to the tip
      double direction[3];
      // note: LE_pts and TE_pts have already been rearranged, so now 0 is root and 1 is tip...
      FOR_I3 direction[i] = 0.5*(LE_pts[1][i] + TE_pts[1][i] - LE_pts[0][i] - TE_pts[0][i]);
      const double mag_direc = MAG(direction); assert(mag_direc>0.0);
      FOR_I3 direction[i] /= mag_direc;
      // check
      if (mpi_rank == 0) cout << "direction = " << COUT_VEC(direction) << endl;

      // determine orientation of the edge loop
      int np = groupSpVecs[idx_root].size();
      double (*xp)[3] = NULL;
      xp = new double[np][3];
      for (int ip = 0; ip < np; ++ip) {
	// let it start from the leading edge
	FOR_I3 xp[ip][i] = surface->xsp[groupSpVecs[idx_root][(ip+LE_idx[idx_root])%np]][i];
      }

      // check
      double length = 0.0;
      for (int i=0; i<np; ++i) {
	length += DIST(xp[i],xp[(i+1)%np]);
      }
      if (mpi_rank == 0) cout << "length = " << length << endl;
      /*
	if (mpi_rank == 0) {
	char filename[128];
	sprintf(filename,"pts_root.dat");
	FILE * fp = fopen(filename,"w");
	fprintf(fp, "VARIABLES = \"X\" \"Y\" \"Z\"\n");
	for (int i = 0; i < np; ++i) {
	fprintf(fp,"%18.15le %18.15le %18.15le\n", xp[i][0], xp[i][1], xp[i][2]);
	}
	fclose(fp);
	}
      */

      int orientation = getOrientation(np, direction, xp);
      // check
      if (mpi_rank == 0) cout << "orientation = " << orientation << endl;

      double (*xp_root)[3] = new double[np][3];  // for ff_surface. should be ordered in the opposite direction
      for (int ip = 0; ip < np; ++ip) {
	if (orientation == 1) {
	  FOR_I3 xp_root[ip][i] = xp[np-1-ip][i];
	} 
	else {
	  FOR_I3 xp_root[ip][i] = xp[ip][i];
	}
      }
      const int np_root = np;  // record the number of root loop pts as well

      // resample the root loop using spline
      int nt;
      double *t_local_root = NULL;
      double (*xt_root)[3] = NULL;
      double len_root;
      {    
	SplineStuff::CubicSpline cspline;
	cspline.init(xp,np,true); // true == periodic spline
	len_root = cspline.getLength();
	const double fraction = 0.1;
	resampleSpline(cspline, nn, dt, dn, fraction, nt, t_local_root, xt_root);
      }
      // check
      if (mpi_rank == 0) {
	cout << "len_root = " << len_root << endl;
	cout << "nt = " << nt << endl;
      }

      // allocate memory
      double (*surf_pts)[3]  = new double[nt*nr][3];
      double (*surf_mesh)[3] = new double[nt*ns][3];

      for (int it=0; it<nt; ++it) {
	FOR_J3 surf_pts[it][j] = xt_root[it][j];
      }
      delete[] xt_root;

      // resample tip now
      np = groupSpVecs[idx_tip].size();
      double (*xp_tmp)[3] = new double[np][3];
      for (int ip = 0; ip < np; ++ip) {
	FOR_I3 xp_tmp[ip][i] = surface->xsp[groupSpVecs[idx_tip][ip]][i];
      }
      const int myorient = getOrientation(np, direction, xp_tmp);

      // check
      /*
	if (mpi_rank == 0) {
	char filename[128];
	sprintf(filename,"pts_tip.dat");
	FILE * fp = fopen(filename,"w");
	fprintf(fp, "VARIABLES = \"X\" \"Y\" \"Z\"\n");
	for (int i = 0; i < np; ++i) {
	fprintf(fp,"%18.15le %18.15le %18.15le\n", xp_tmp[i][0], xp_tmp[i][1], xp_tmp[i][2]);
	}
	fclose(fp);
	}
      */

      delete[] xp;
      xp = new double[np][3];
      if (myorient == orientation) {
	for (int ip = 0; ip < np; ++ip) {
	  // let it start from the leading edge
	  FOR_I3 xp[ip][i] = xp_tmp[(ip+LE_idx[idx_tip])%np][i];
	}
      }
      else {
	for (int ip = 0; ip < np; ++ip) {
	  // let it start from the leading edge
	  FOR_I3 xp[ip][i] = xp_tmp[(LE_idx[idx_tip]-ip+np)%np][i];
	}
      }
      double (*xp_tip)[3] = new double[np][3];  // for ff_surface. should be ordered in the opposite direction
      for (int ip = 0; ip < np; ++ip) {
	if (orientation == 1) {
	  FOR_I3 xp_tip[ip][i] = xp[ip][i];
	}
	else {
	  FOR_I3 xp_tip[ip][i] = xp[np-1-ip][i];
	}
      }
      const int np_tip = np;  // record the number of tip loop pts as well

      int nt_tip;
      double *t_local_tip = NULL;
      double (*xt_tip)[3] = NULL;
      double len_tip; 
      {
	SplineStuff::CubicSpline cspline;
	cspline.init(xp,np,true);
	len_tip = cspline.getLength();
	double fraction0 = 0.09;
	resampleSpline(cspline, nn, dt, dn, fraction0, nt_tip, t_local_tip, xt_tip);
	if (nt_tip!=nt) {
	  int nt_tip0 = nt_tip;
	  DELETE(t_local_tip);
	  DELETE(xt_tip); 
	  double fraction = 0.11;
	  resampleSpline(cspline, nn, dt, dn, fraction, nt_tip, t_local_tip, xt_tip);

	  // iterate to make nt_tip = nt
	  int nt_tip1;
	  double fraction1;
	  int iter=0;
	  while (nt_tip!=nt) {
	    ++iter;
	    nt_tip1 = nt_tip;
	    fraction1 = fraction;
	    if (nt_tip1 != nt_tip0) {
	      fraction = fraction0 + double(nt-nt_tip0)/double(nt_tip1-nt_tip0)*(fraction1-fraction0);
	    }
	    else {
	      fraction = fraction1 + (fraction1 - fraction0)*1.05;
	    }
	    DELETE(t_local_tip);
	    DELETE(xt_tip);
	    if (mpi_rank == 0) cout << "fraction = " << fraction << endl;
	    resampleSpline(cspline, nn, dt, dn, fraction, nt_tip, t_local_tip, xt_tip);
	    nt_tip0 = nt_tip1;
	    fraction0 = fraction1;
	  }
	  if (mpi_rank == 0) cout << "iter = " << iter << endl;
	}
      }

      for (int it = 0; it < nt; ++it) {
	//cspline.getX(surf_pts[(nr-1)*nt+it],t_local_tip[it]);
	FOR_J3 surf_pts[(nr-1)*nt+it][j] = xt_tip[it][j];
      }
   
      delete[] xt_tip; 
      delete[] xp_tmp;
      delete[] xp;

      /*
	char filename[128];
	sprintf(filename,"pts_root.dat");
	FILE * fp = fopen(filename,"w");
	fprintf(fp, "VARIABLES = \"X\" \"Y\" \"Z\"\n");
	for (int i = 0; i < nt; ++i) {
	fprintf(fp,"%18.15le %18.15le %18.15le\n", xt[i][0], xt[i][1], xt[i][2]);
	}
	fclose(fp);
      */

    
      double *t_local = new double[nt];

      // intersect a tcone, get the profiles and resample
      double LE_line[nr][3];
      {
	// note: LE_pts and TE_pts have already been rearranged, so now 0 is root and 1 is tip...
	FOR_I3 LE_line[0][i] = LE_pts[0][i];
	FOR_I3 LE_line[nr-1][i] = LE_pts[1][i];
	const double dx[3] = {(LE_line[nr-1][0]-LE_line[0][0])/double(nr-1), 
			      (LE_line[nr-1][1]-LE_line[0][1])/double(nr-1), 
			      (LE_line[nr-1][2]-LE_line[0][2])/double(nr-1)};
	for (int i=1; i<=nr-2; ++i) {
	  FOR_J3 LE_line[i][j] = LE_line[0][j] + dx[j]*double(i);
	}
      }

      double TE_line[nr][3];
      {
	// note: LE_pts and TE_pts have already been rearranged, so now 0 is root and 1 is tip...
	FOR_I3 TE_line[0][i] = TE_pts[0][i];
	FOR_I3 TE_line[nr-1][i] = TE_pts[1][i];
	const double dx[3] = {(TE_line[nr-1][0]-TE_line[0][0])/double(nr-1), 
			      (TE_line[nr-1][1]-TE_line[0][1])/double(nr-1), 
			      (TE_line[nr-1][2]-TE_line[0][2])/double(nr-1)};
	for (int i=1; i<=nr-2; ++i) {
	  FOR_J3 TE_line[i][j] = TE_line[0][j] + dx[j]*double(i);
	}
      }

      // check
      /*
	cout << "> leading edge line: " << endl;
	for (int i=0; i<nr; ++i) {
	cout << LE_line[i][0] << " " << LE_line[i][1] << " " << LE_line[i][2] << endl;
	}
	cout << "> trailing edge line: " << endl;
	for (int i=0; i<nr; ++i) {
	cout << TE_line[i][0] << " " << TE_line[i][1] << " " << TE_line[i][2] << endl;
	}
      */

      vector<vector<double> > xp_r;
      for (int ir=1; ir<nr-1; ++ir) {
	//debug
	//for (int ir=1; ir<2; ++ir) {
	const double xp_le[3] = {LE_line[ir][0], LE_line[ir][1], LE_line[ir][2]};
	const double xp_te[3] = {TE_line[ir][0], TE_line[ir][1], TE_line[ir][2]};
	double dp;
	double vec_tmp[3];

	dp = DOT_PRODUCT(xp_le,axis);
	FOR_I3 vec_tmp[i] = xp_le[i] - dp*axis[i];
	const double ra = MAG(vec_tmp);

	dp = DOT_PRODUCT(xp_te,axis);
	FOR_I3 vec_tmp[i] = xp_te[i] - dp*axis[i];
	const double rb = MAG(vec_tmp);

	// check
	if (mpi_rank == 0) cout << "ra, rb = " << ra << " " << rb << endl;

	//vector<array<double, 3> > xp_r;
	surface->intersectWithTcone(axis, ra, rb, xp_le, xp_te, xp_r, st_flag);

	//cout << "points after intersection " << xp_r.size() << endl;
	//for (vector<array<double,3> >::const_iterator i = xp_r.begin(); i != xp_r.end(); ++i) {
	//  cout << (*i)[0] << " " << (*i)[1] << " " << (*i)[2] << endl;
	//}
     
	// check
	/*
	  if (mpi_rank == 0) { 
	  char filename[128];
	  sprintf(filename,"pts_intersect.%03d.dat",ir);
	  FILE * fp = fopen(filename,"w");
	  fprintf(fp, "VARIABLES = \"X\" \"Y\" \"Z\"\n");
	  for (vector<vector<double> >::const_iterator i = xp_r.begin(); i != xp_r.end(); ++i) {
          fprintf(fp,"%18.15le %18.15le %18.15le\n", (*i)[0], (*i)[1], (*i)[2]);
	  }
	  fclose(fp);
	  }
	*/

	const int n_pts = xp_r.size()-1;     //note: the last and first pts in xp_r are the same. remove the last pt
	double (*pts_tmp)[3] = new double[n_pts][3];
	for (int i=0; i<n_pts; ++i) {
	  FOR_J3 pts_tmp[i][j] = xp_r[i][j];
	}

	double max_dist = 0.0;
	int idx0, idx1;
	for (int i=0; i<n_pts; ++i) {
	  for (int j=i+1; j<=n_pts; ++j) {
	    const double dist_ij = DIST(pts_tmp[i],pts_tmp[j%n_pts]);
	    if (dist_ij > max_dist) {
	      max_dist = dist_ij;
	      idx0 = i;
	      idx1 = j;
	    }
	  }
	}
      
	double pt0[3], pt1[3];
	FOR_I3 pt0[i] = pts_tmp[idx0][i];
	FOR_I3 pt1[i] = pts_tmp[idx1][i];
	double vec01[3];
	FOR_I3 vec01[i] = pt1[i] - pt0[i];
	int start_idx;
	dp = DOT_PRODUCT(vec01, ref_vec);
	if (dp >= 0) {
	  start_idx = idx0;
	}
	else {
	  start_idx = idx1;
	}

	// determine the orientation
	double (*pts)[3] = new double[n_pts][3];
	const int myorient = getOrientation(n_pts, direction, pts_tmp);
	if (myorient == orientation) {
	  for (int i=0; i<n_pts; ++i) {
	    FOR_J3 pts[i][j] = pts_tmp[(i+start_idx)%n_pts][j];
	  }
	}
	else {
	  for (int i=0; i<n_pts; ++i) {
	    FOR_J3 pts[i][j] = pts_tmp[(start_idx-i+n_pts)%n_pts][j];
	  }
	}

	// check
	/*
	  if (mpi_rank == 0) { 
	  char filename[128];
	  sprintf(filename,"reordered_pts_intersect.%03d.dat",ir);
	  FILE * fp = fopen(filename,"w");
	  fprintf(fp, "VARIABLES = \"X\" \"Y\" \"Z\"\n");
	  for (int i = 0; i < n_pts; ++i) {
          fprintf(fp,"%18.15le %18.15le %18.15le\n", pts[i][0], pts[i][1], pts[i][2]);
	  }
	  fclose(fp);
	  }
	*/

	SplineStuff::CubicSpline cspline;
	cspline.init(pts,n_pts,true);
	const double fac = double(ir)/double(nr-1);
	for (int it = 0; it < nt; ++it) {
	  t_local[it] = t_local_root[it] + fac*(t_local_tip[it] - t_local_root[it]);
	  cspline.getX(surf_pts[ir*nt+it],t_local[it]);
	}

	delete[] pts_tmp;
	delete[] pts;
	xp_r.clear();
      }

      delete[] st_flag;

      // check
      /*
	if (mpi_rank == 0) {
	char filename[128];
	sprintf(filename,"surf_pts.dat");
	FILE * fp = fopen(filename,"w");
	fprintf(fp, "VARIABLES = \"X\" \"Y\" \"Z\"\n");
	fprintf(fp, "ZONE I=%d, J=%d\n", nt, nr);
	for (int i=0; i<nt*nr; ++i) {
        fprintf(fp,"%18.15le %18.15le %18.15le\n", surf_pts[i][0], surf_pts[i][1], surf_pts[i][2]);
	}
	fclose(fp);
	}
      */

      // final step, resample in the r direction to get the surface mesh
      /*
	double *s_local = new double[ns];
	const double ds_local = 1.0/double(ns-1);
	s_local[0] = 0.0;
	s_local[ns-1] = 1.0;
	for (int i=1; i<ns-1; ++i) {
	s_local[i] = double(i)*ds_local;
	}
      */
    
      xp = new double[nr][3];
      for (int it=0; it<nt; ++it) {
	for (int ir=0; ir<nr; ++ir) {
	  FOR_I3 xp[ir][i] = surf_pts[ir*nt+it][i];
	}
	SplineStuff::CubicSpline cspline;
	cspline.init(xp,nr,false);
	for (int is = 0; is < ns; ++is) {
	  cspline.getX(surf_mesh[is*nt+it],s_local[is]);
	}
      }	

      //char filename[128];
      //FILE * fp = NULL;   
 
      // check
      /*
	if (mpi_rank == 0) {
	sprintf(filename,"surf_mesh.dat");
	fp = fopen(filename,"w");
	fprintf(fp, "VARIABLES = \"X\" \"Y\" \"Z\"\n");
	fprintf(fp, "ZONE I=%d, J=%d\n", nt, ns);
	for (int i=0; i<nt*ns; ++i) {
        fprintf(fp,"%18.15le %18.15le %18.15le\n", surf_mesh[i][0], surf_mesh[i][1], surf_mesh[i][2]);
	}
	fclose(fp);
	}
      */

      // outer pts coordinates
      double (*outer_pts)[3] = new double[nt*ns][3];
      for (int is = 0; is < ns; ++is) {
	for (int it = 0; it < nt; ++it) {
	  const int it_prev = (it-1+nt)%nt;
	  const int it_next = (it+1)%nt;
	  const double dxt[3] = DIFF(surf_mesh[is*nt+it_next],surf_mesh[is*nt+it_prev]);
	  const int is_prev = (is-1<0)?0:is-1;
	  const int is_next = (is+1>=ns)?ns-1:is+1;
	  const double dxs[3] = DIFF(surf_mesh[is_next*nt+it],surf_mesh[is_prev*nt+it]);
	  double cp[3] = CROSS_PRODUCT(dxt,dxs);
	  FOR_I3 cp[i] *= double(orientation);
	  const double cp_mag = MAG(cp); assert(cp_mag > 0.0);
	  const double delta_bl = GeomUtils::stretchingFunction(nn,nn,dn,dt);
	  FOR_I3 outer_pts[is*nt+it][i] = surf_mesh[is*nt+it][i] + delta_bl*cp[i]/cp_mag;
	}
      }

      // smooth the outer points...
      const int nsmooth2 = 200;
      double (*outer_tmp)[3] = new double[nt][3];
      for (int ismooth = 0; ismooth < nsmooth2; ++ismooth) {
	if (mpi_rank==0 && ismooth%10==0) cout << "ismooth = " << ismooth << endl;
	for (int is = 0; is < ns; ++is) {
	  for (int it = 0; it < nt; ++it)
	    FOR_I3 outer_tmp[it][i] = outer_pts[is*nt+it][i];
	  for (int it = 0; it < nt; ++it) {
	    const int it_prev = (it-1+nt)%nt;
	    const int it_next = (it+1)%nt;
	    const double dx1[3] = DIFF(outer_tmp[it_next],outer_tmp[it_prev]);
	    const double dx1_mag2 = DOT_PRODUCT(dx1,dx1);
	    double dx2[3]; FOR_I3 dx2[i] = 0.5*(outer_tmp[it_next][i]+outer_tmp[it_prev][i])-outer_tmp[it][i];
	    const double dp = DOT_PRODUCT(dx2,dx1);
	    FOR_I3 outer_pts[is*nt+it][i] += dp*dx1[i]/dx1_mag2;
	  }
	}
      }
      delete[] outer_tmp;

      // check
      /*
	if (mpi_rank == 0) {
	sprintf(filename,"outer_pts.dat");
	fp = fopen(filename,"w");
	fprintf(fp, "VARIABLES = \"X\" \"Y\" \"Z\"\n");
	fprintf(fp, "ZONE I=%d, J=%d\n", nt, ns);
	for (int i=0; i<nt*ns; ++i) {
        fprintf(fp,"%18.15le %18.15le %18.15le\n", outer_pts[i][0], outer_pts[i][1], outer_pts[i][2]);
	}
	fclose(fp);    
	}
      */

      // and finally build the strand points between the inner and outer layers...

      double (*bl_mesh)[3] = new double[nt*ns*(nn+1)][3];

      for (int ipt = 0; ipt < nt*ns; ++ipt)
	FOR_I3 bl_mesh[ipt][i] = surf_mesh[ipt][i];

      for (int ipt = 0; ipt < nt*ns; ++ipt)
	FOR_I3 bl_mesh[nt*ns*nn+ipt][i] = outer_pts[ipt][i];
    
      for (int is = 0; is < ns; ++is) {
	for (int it = 0; it < nt; ++it) {
	  const int it_prev = (it-1+nt)%nt;
	  const int it_next = (it+1)%nt;
	  const double dxt[3] = DIFF(surf_mesh[is*nt+it_next],surf_mesh[is*nt+it_prev]);
	  const int is_prev = (is-1<0)?0:is-1;
	  const int is_next = (is+1>=ns)?ns-1:is+1;
	  const double dxs[3] = DIFF(surf_mesh[is_next*nt+it],surf_mesh[is_prev*nt+it]);
	  double cp[3] = CROSS_PRODUCT(dxt,dxs);
	  FOR_I3 cp[i] *= double(orientation);
	  const double cp_mag = MAG(cp); assert(cp_mag > 0.0);
	  const double delta_bl = GeomUtils::stretchingFunction(nn,nn,dn,dt);
	  double dx[3]; FOR_I3 dx[i] = outer_pts[is*nt+it][i] - (surf_mesh[is*nt+it][i] + delta_bl*cp[i]/cp_mag);
	  for (int in = 1; in < nn; ++in) {
	    const double delta = GeomUtils::stretchingFunction(in,nn,dn,dt);
	    const double ratio = delta/delta_bl; assert((ratio > 0.0)&&(ratio < 1.0));
	    FOR_I3 bl_mesh[nt*ns*in+is*nt+it][i] = surf_mesh[is*nt+it][i] + ratio*delta_bl*cp[i]/cp_mag + ratio*ratio*dx[i];
	  }
	}
      }

      /*
      // check
      sprintf(filename,"bl_mesh.dat");
      fp = fopen(filename,"w");
      fprintf(fp, "VARIABLES = \"X\" \"Y\" \"Z\"\n");
      fprintf(fp, "ZONE I=%d, J=%d, K=%d\n", nt, ns, nn+1);
      for (int i=0; i<nt*ns*(nn+1); ++i) {
      fprintf(fp,"%18.15le %18.15le %18.15le\n", bl_mesh[i][0], bl_mesh[i][1], bl_mesh[i][2]);
      }
      fclose(fp);
      */

      // no longer need mypart -- delete it
      //delete mypart;


      /*
      // =======================================
      // the surface...
      // =======================================

      assert(part->surface == NULL);
      part->surface = new SurfaceShm();
      part->surface->init(nt*ns,nt*(ns-1)*2); // nsp,nst
      // debug
      if (mpi_rank == 0) cout << "nt*ns = " << nt*ns << endl;
      if (b_name) {
      // if the user gave the NAME param, then use this to set the zone name...
      part->surface->zoneVec.push_back(StZone(name));
      }
      else {
      // otherwise, use the profile name, or a default...
      part->surface->zoneVec.push_back(StZone("EXTRUDE_AIRFOIL"));
      }
      assert(part->solid_type == UNKNOWN_SOLID);
      part->solid_type = SURFACE_SOLID;

      // recall surface is shm...
      if (mpi_rank_shared == 0) {
      for (int is = 0; is < ns; ++is) {
      for (int it = 0; it < nt; ++it) {
      FOR_I3 part->surface->xsp[is*nt+it][i] = surf_mesh[is*nt+it][i];
      }
      }
      for (int is = 0; is < ns-1; ++is) {
      for (int it = 0; it < nt; ++it) {
      const int it_next = (it+1)%nt;
      // pair of tris...
      part->surface->spost[(is*nt+it)*2  ][0] = is*nt+it;
      part->surface->spost[(is*nt+it)*2  ][2] = is*nt+it+nt;
      part->surface->spost[(is*nt+it)*2  ][1] = is*nt+it_next;
      part->surface->znost[(is*nt+it)*2  ] = 0;
      part->surface->spost[(is*nt+it)*2+1][0] = is*nt+it_next;
      part->surface->spost[(is*nt+it)*2+1][2] = is*nt+it+nt;
      part->surface->spost[(is*nt+it)*2+1][1] = is*nt+it_next+nt;
      part->surface->znost[(is*nt+it)*2+1] = 0;
      }
      }
      //GeomUtils::writeTecplot("surface.dat",part->surface->spost,part->surface->nst,part->surface->xsp);
      }
      MPI_Barrier(mpi_comm_shared);
      */


      // =======================================
      // the ff_surface...
      // =======================================

      assert(part->ff_surface == NULL);
      part->ff_surface = new SurfaceShm();
      const int ff_nsp = nt*ns+np_root+np_tip;
      const int ff_nst = nt*(ns-1)*2 + (nt+np_root) + (nt+np_tip);
      part->ff_surface->init(ff_nsp,ff_nst); // nsp,nst
      if (b_name) {
	// if the user gave the NAME param, then use this to set the zone name...
	part->ff_surface->zoneVec.push_back(StZone(name+"_FF"));
      }
      else {
	// otherwise, use the profile name, or a default...
	part->ff_surface->zoneVec.push_back(StZone("EXTRUDE_AIRFOIL_FF"));
      }
      assert(part->ff_surface_dxost == NULL);
      CTI_Mmap(part->ff_surface_dxost,ff_nst);
      assert(part->ff_type == UNKNOWN_FF);
      //part->ff_type = FF_SURFACE_FF;
      part->ff_type = FF_SURFACE_FAZONE_FF;
      assert(part->solid_type == UNKNOWN_SOLID);
      part->solid_type = NO_SOLID;

      // recall surface is shm...
      if (mpi_rank_shared == 0) {
	for (int is = 0; is < ns; ++is) {
	  for (int it = 0; it < nt; ++it) {
	    if (is<=nn) {
	      const int it_prev = (it-1+nt)%nt;
	      const int it_next = (it+1)%nt;
	      const double dxt[3] = DIFF(surf_mesh[is*nt+it_next],surf_mesh[is*nt+it_prev]);
	      const int is_prev = (is-1<0)?0:is-1;
	      const int is_next = (is+1>=ns)?ns-1:is+1;
	      const double dxs[3] = DIFF(surf_mesh[is_next*nt+it],surf_mesh[is_prev*nt+it]);
	      double cp[3] = CROSS_PRODUCT(dxt,dxs);
	      FOR_I3 cp[i] *= double(orientation);
	      const double cp_mag = MAG(cp); assert(cp_mag > 0.0);
	      // lift the points by a small distance in the wall-normal direction
	      FOR_I3 part->ff_surface->xsp[is*nt+it][i] = bl_mesh[nt*ns*is+is*nt+it][i] + dn*double(nn-is)/double(nn)*cp[i]/cp_mag;
	    }
	    else if (is>=nn+ns_mid) {
	      const int it_prev = (it-1+nt)%nt;
	      const int it_next = (it+1)%nt;
	      const double dxt[3] = DIFF(surf_mesh[is*nt+it_next],surf_mesh[is*nt+it_prev]);
	      const int is_prev = (is-1<0)?0:is-1;
	      const int is_next = (is+1>=ns)?ns-1:is+1;
	      const double dxs[3] = DIFF(surf_mesh[is_next*nt+it],surf_mesh[is_prev*nt+it]);
	      double cp[3] = CROSS_PRODUCT(dxt,dxs);
	      FOR_I3 cp[i] *= double(orientation);
	      const double cp_mag = MAG(cp); assert(cp_mag > 0.0);
	      // lift the points by a small distance in the wall-normal direction
	      FOR_I3 part->ff_surface->xsp[is*nt+it][i] = bl_mesh[nt*ns*(ns-1-is)+is*nt+it][i] + dn*double(is-nn-ns_mid)/double(nn)*cp[i]/cp_mag;
	    }
	    else {
	      FOR_I3 part->ff_surface->xsp[is*nt+it][i] = outer_pts[is*nt+it][i];
	    }
	  }
	}
	for (int ip = 0; ip < np_root; ++ip) {
	  FOR_I3 part->ff_surface->xsp[ns*nt+ip][i] = xp_root[ip][i];
	}
	for (int ip = 0; ip < np_tip; ++ip) {
	  FOR_I3 part->ff_surface->xsp[ns*nt+np_root+ip][i] = xp_tip[ip][i];
	}


	for (int is = 0; is < ns-1; ++is) {
	  for (int it = 0; it < nt; ++it) { 
	    const int it_next = (it+1)%nt;
	    double dist0, dist1;
	    // pair of tris...
	    part->ff_surface->spost[(is*nt+it)*2  ][0] = is*nt+it;
	    part->ff_surface->spost[(is*nt+it)*2  ][1] = is*nt+it+nt;
	    part->ff_surface->spost[(is*nt+it)*2  ][2] = is*nt+it_next;
	    part->ff_surface->znost[(is*nt+it)*2  ] = 0;
	    if (is>nn && is<nn+ns_mid) {
	      part->ff_surface_dxost[(is*nt+it)*2  ] = 0.5*dt;
	    }
	    else {
	      dist0 = DIST(part->ff_surface->xsp[is*nt+it], part->ff_surface->xsp[is*nt+it+nt]);
	      dist1 = DIST(part->ff_surface->xsp[is*nt+it], part->ff_surface->xsp[is*nt+it_next]);
	      part->ff_surface_dxost[(is*nt+it)*2  ] = 0.5*(dist0+dist1);
	    }
	    part->ff_surface->spost[(is*nt+it)*2+1][0] = is*nt+it_next;
	    part->ff_surface->spost[(is*nt+it)*2+1][1] = is*nt+it+nt;
	    part->ff_surface->spost[(is*nt+it)*2+1][2] = is*nt+it_next+nt;
	    part->ff_surface->znost[(is*nt+it)*2+1] = 0;
	    if (is>nn && is<nn+ns_mid) {
	      part->ff_surface_dxost[(is*nt+it)*2+1] = 0.5*dt;
	    }
	    else {
	      dist0 = DIST(part->ff_surface->xsp[is*nt+it_next+nt], part->ff_surface->xsp[is*nt+it+nt]);
	      dist1 = DIST(part->ff_surface->xsp[is*nt+it_next+nt], part->ff_surface->xsp[is*nt+it_next]);
	      part->ff_surface_dxost[(is*nt+it)*2+1] = 0.5*(dist0+dist1);
	    }
	  }
	}

	int num_st;
	vector<int> isp0Vec, isp1Vec;
	//double (*xsp_tmp)[3] = new double[nt+np_root][3];
	for (int ip = 0; ip < nt; ++ip) {
	  if (orientation == 1) {
	    isp0Vec.push_back(ip);
	  }
	  else {
	    isp0Vec.push_back(nt-1-ip);
	  }
	  //FOR_I3 xsp_tmp[ip][i] = part->ff_surface->xsp[ip][i];
	}
	for (int ip = 0; ip < np_root; ++ip) {
	  isp1Vec.push_back(ns*nt+ip);
	  //FOR_I3 xsp_tmp[nt+ip][i] = xp_root[ip][i];
	}

	/*
	  sprintf(filename,"ff_root.dat");
	  fp = fopen(filename,"w");
	  fprintf(fp, "VARIABLES = \"X\" \"Y\" \"Z\"\n");
	  for (int i = 0; i < nt; ++i) {
	  fprintf(fp,"%18.15le %18.15le %18.15le\n", part->ff_surface->xsp[isp0Vec[i]][0], part->ff_surface->xsp[isp0Vec[i]][1], part->ff_surface->xsp[isp0Vec[i]][2]);
	  }
	  fclose(fp);
      
	  sprintf(filename,"surf_root.dat");
	  fp = fopen(filename,"w");
	  fprintf(fp, "VARIABLES = \"X\" \"Y\" \"Z\"\n");
	  for (int i = 0; i < np_root; ++i) {
	  fprintf(fp,"%18.15le %18.15le %18.15le\n", part->ff_surface->xsp[isp1Vec[i]][0], part->ff_surface->xsp[isp1Vec[i]][1], part->ff_surface->xsp[isp1Vec[i]][2]);
	  }
	  fclose(fp);
	*/

	num_st = GeomUtils::facetGap(part->ff_surface->spost+nt*(ns-1)*2, isp1Vec, isp0Vec, part->ff_surface->xsp, true);
	assert(num_st == nt+np_root);
	for (int ist = nt*(ns-1)*2; ist < nt*(ns-1)*2 + (nt+np_root); ++ist) {
	  part->ff_surface->znost[ist] = 0;
	  part->ff_surface_dxost[ist] = dn*2.0;
	}
	isp0Vec.clear();
	isp1Vec.clear();
	//delete[] xsp_tmp;

	//xsp_tmp = new double[nt+np_tip][3];
	for (int ip = 0; ip < nt; ++ip) {
	  if (orientation == 1) {
	    isp0Vec.push_back((ns-1)*nt+nt-1-ip);
	  }
	  else {
	    isp0Vec.push_back((ns-1)*nt+ip);
	  }
	  //FOR_I3 xsp_tmp[ip][i] = part->ff_surface->xsp[(ns-1)*nt+ip][i];
	}
	for (int ip = 0; ip < np_tip; ++ip) {
	  isp1Vec.push_back(ns*nt+np_root+ip);
	  //FOR_I3 xsp_tmp[nt+ip][i] = xp_tip[ip][i];
	}

	/*
	  sprintf(filename,"ff_tip.dat");
	  fp = fopen(filename,"w");
	  fprintf(fp, "VARIABLES = \"X\" \"Y\" \"Z\"\n");
	  for (int i = 0; i < nt; ++i) {
	  fprintf(fp,"%18.15le %18.15le %18.15le\n", part->ff_surface->xsp[isp0Vec[i]][0], part->ff_surface->xsp[isp0Vec[i]][1], part->ff_surface->xsp[isp0Vec[i]][2]);
	  }
	  fclose(fp);

	  sprintf(filename,"surf_tip.dat");
	  fp = fopen(filename,"w");
	  fprintf(fp, "VARIABLES = \"X\" \"Y\" \"Z\"\n");
	  for (int i = 0; i < np_tip; ++i) { 
	  fprintf(fp,"%18.15le %18.15le %18.15le\n", part->ff_surface->xsp[isp1Vec[i]][0], part->ff_surface->xsp[isp1Vec[i]][1], part->ff_surface->xsp[isp1Vec[i]][2]);
	  }
	  fclose(fp);
	*/

	num_st = GeomUtils::facetGap(part->ff_surface->spost+nt*(ns-1)*2+nt+np_root, isp1Vec, isp0Vec, part->ff_surface->xsp, true);
	assert(num_st == nt+np_tip);
	for (int ist = nt*(ns-1)*2 + (nt+np_root); ist < ff_nst; ++ist) {
	  part->ff_surface->znost[ist] = 0;
	  part->ff_surface_dxost[ist] = dn*2.0;
	}
	isp0Vec.clear();
	isp1Vec.clear();
	//delete[] xsp_tmp;      

	if (mpi_rank == 0) {
	  //GeomUtils::writeTecplot("ff_surface.dat",part->ff_surface->spost,part->ff_surface->nst,part->ff_surface->xsp);
	  GeomUtils::writeSbin("ff_surface.sbin",part->ff_surface->spost,part->ff_surface->nst,part->ff_surface->xsp);
	}
      }
      MPI_Barrier(mpi_comm_shared);

      // =======================================
      // and finally the points: this part is distributed in parallel...
      // =======================================
      assert(part->pts == NULL);
      part->pts = new Points();
      int nt_local = nt/mpi_size;
      if (nt%mpi_size != 0) nt_local += 1;
      int it_start = min(nt,nt_local*mpi_rank);
      int it_end = min(nt,nt_local*(mpi_rank+1));
      assert(((ns_mid+ns-1)*nn)%2==0);
      const int nsn = (ns_mid+ns+1) * nn / 2;  // number of pts in s and n directions
      part->pts->np = (it_end - it_start) * nsn;
      part->pts->np_global = nt*nsn;

      if (mpi_rank == 0) cout << " > np_global: " << part->pts->np_global << endl;

      part->pts->xp = new double[part->pts->np][3];
      part->pts->delta = new double[part->pts->np];

      int ip = 0;
      for (int it = it_start; it < it_end; ++it) {
	for (int in = 0; in < nn; ++in) {
	  const int s_start = in;
	  const int s_end   = ns-1-in;
	  for (int is = s_start; is < s_end; ++is) {
	    const int it_next = (it+1)%nt;
	    const int is_next = is+1;
	    const int in_next = in+1;
	    double pt0[3]; FOR_I3 pt0[i] = bl_mesh[nt*ns*in      + is     *nt + it     ][i];
	    double pt1[3]; FOR_I3 pt1[i] = bl_mesh[nt*ns*in      + is     *nt + it_next][i];
	    double pt2[3]; FOR_I3 pt2[i] = bl_mesh[nt*ns*in      + is_next*nt + it_next][i];
	    double pt3[3]; FOR_I3 pt3[i] = bl_mesh[nt*ns*in      + is_next*nt + it     ][i];
	    double pt4[3]; FOR_I3 pt4[i] = bl_mesh[nt*ns*in_next + is     *nt + it     ][i];
	    double pt5[3]; FOR_I3 pt5[i] = bl_mesh[nt*ns*in_next + is     *nt + it_next][i];
	    double pt6[3]; FOR_I3 pt6[i] = bl_mesh[nt*ns*in_next + is_next*nt + it_next][i];
	    double pt7[3]; FOR_I3 pt7[i] = bl_mesh[nt*ns*in_next + is_next*nt + it     ][i];
	    FOR_I3 part->pts->xp[ip][i] = (pt0[i] + pt1[i] + pt2[i] + pt3[i] + pt4[i] + pt5[i] + pt6[i] + pt7[i]) / 8.0;
	    const double diag0 = DIST(pt0, pt6);
	    const double diag1 = DIST(pt1, pt7);
	    const double diag2 = DIST(pt2, pt4);
	    const double diag3 = DIST(pt3, pt5);
	    part->pts->delta[ip] = 1.5*MAX4(diag0, diag1, diag2, diag3);  // take stretching into account
	    ++ip;
	  }
	}
      }
      //cout << "ip, part->pts->np = " << ip << " " << part->pts->np << endl;
      assert(ip == part->pts->np);

      //GeomUtils::writePtsTecplot("points.dat",part->pts->xp,part->pts->np);

      /*
      // check
      // only for debugging...
      double (*check_pts)[3] = new double[nt*(ns-1)*nn][3];

      for (int in = 0; in < nn; ++in) {
      for (int is = 0; is < ns-1; ++is) {
      for (int it = 0; it < nt; ++it) {
      const int it_next = (it+1)%nt;
      const int is_next = is+1;
      const int in_next = in+1;
      double pt0[3]; FOR_I3 pt0[i] = bl_mesh[nt*ns*in      + is     *nt + it     ][i];
      double pt1[3]; FOR_I3 pt1[i] = bl_mesh[nt*ns*in      + is     *nt + it_next][i];
      double pt2[3]; FOR_I3 pt2[i] = bl_mesh[nt*ns*in      + is_next*nt + it_next][i];
      double pt3[3]; FOR_I3 pt3[i] = bl_mesh[nt*ns*in      + is_next*nt + it     ][i];
      double pt4[3]; FOR_I3 pt4[i] = bl_mesh[nt*ns*in_next + is     *nt + it     ][i];
      double pt5[3]; FOR_I3 pt5[i] = bl_mesh[nt*ns*in_next + is     *nt + it_next][i];
      double pt6[3]; FOR_I3 pt6[i] = bl_mesh[nt*ns*in_next + is_next*nt + it_next][i];
      double pt7[3]; FOR_I3 pt7[i] = bl_mesh[nt*ns*in_next + is_next*nt + it     ][i];
      FOR_I3 check_pts[nt*(ns-1)*in+is*nt+it][i] = (pt0[i] + pt1[i] + pt2[i] + pt3[i] + pt4[i] + pt5[i] + pt6[i] + pt7[i]) / 8.0;
      }
      }
      }
      sprintf(filename,"check_pts.dat");
      fp = fopen(filename,"w");
      fprintf(fp, "VARIABLES = \"X\" \"Y\" \"Z\"\n");
      fprintf(fp, "ZONE I=%d, J=%d, K=%d\n", nt, ns-1, nn);
      for (int i=0; i<nt*ns*(nn+1); ++i) {
      fprintf(fp,"%18.15le %18.15le %18.15le\n", check_pts[i][0], check_pts[i][1], check_pts[i][2]);
      }
      fclose(fp);

      delete[] check_pts;
      */

      // all grouped st's get their part index set to partVec.size() (note that
      // the part passed into this routine will be pushed onto partVec after we
      // return succesfully)...

      // HEREHERE

      if (mpi_rank_shared == 0) {
	for (int ipart = 0, limit = partVec.size(); ipart < limit; ++ipart) {
	  if (partVec[ipart]->surface) {
	    SurfaceShm * surface = partVec[ipart]->surface;
	    for (int ist = 0; ist < surface->nst; ++ist) {
	      int igr;
	      if (partVec[ipart]->getGroupForSt(igr,ist)) {
		partVec[ipart]->clearGroupForSt(ist);
		partVec[ipart]->setPartForSt(part->ipart_of_part,ist);
	      }
	    }
	  }
	}
      }
      MPI_Barrier(mpi_comm_shared);

      // and clear the groupDataVec: we don't use the same clear function as in
      // the queryZones routine becasue we have already cleared the ist group above

      for (int igr = 0; igr < groupDataVec.size(); ++igr) {
	groupDataVec[igr].clear();
      }
      groupDataVec.clear();

      delete[] xp;
      delete[] xp_root;
      delete[] xp_tip;
      delete[] t_local;
      delete[] s_local;
      delete[] surf_pts;
      delete[] surf_mesh;
      delete[] outer_pts;
      delete[] bl_mesh;

      //MPI_Pause("OK -- lets do this!");

      assert(ierr == 0);
      return true;

    } //if (groupDataVec.size() == 1)
    else {
      WUI(WARN," > AIRFOIL_STRAND currently supports single disjoint group. Returning...");

      if (mpi_rank_shared == 0) {
        for (int ipart = 0, limit = partVec.size(); ipart < limit; ++ipart) {
          if (partVec[ipart]->surface) {
            SurfaceShm * surface = partVec[ipart]->surface;
            for (int ist = 0; ist < surface->nst; ++ist) {
              int igr;
              if (partVec[ipart]->getGroupForSt(igr,ist)) {
                partVec[ipart]->clearGroupForSt(ist);
              }
            }
          }
        }
      }
      MPI_Barrier(mpi_comm_shared);
      for (int igr = 0; igr < groupDataVec.size(); ++igr) {
        groupDataVec[igr].clear();
      }
      groupDataVec.clear();

      return false;
    }

  }

  int getOrientation(const int np, const double direction[3], const double (*xp)[3]) {
    double theta = 0.0;
    for (int ip=1; ip<np+1; ++ip) {
      double dp;
      double vec0[3], vec1[3];
      FOR_I3 vec0[i] = xp[ip%np][i] - xp[ip-1][i];
      FOR_I3 vec1[i] = xp[(ip+1)%np][i] - xp[ip%np][i];
      dp = DOT_PRODUCT(vec0,direction);
      FOR_I3 vec0[i] -= dp*direction[i];
      dp = DOT_PRODUCT(vec1,direction);
      FOR_I3 vec1[i] -= dp*direction[i];
      const double cp[3] = CROSS_PRODUCT(vec1,vec0);
      dp = DOT_PRODUCT(cp,direction);
      const double mag_vec0 = MAG(vec0); assert(mag_vec0>0.0);
      const double mag_vec1 = MAG(vec1); assert(mag_vec1>0.0);
      const double cos_t = DOT_PRODUCT(vec0,vec1);
      double dtheta;
      if (cos_t>=0.0) {
        dtheta = asin(dp/mag_vec0/mag_vec1);
      }
      else {
        if (dp>=0.0) {
          dtheta = M_PI - asin(dp/mag_vec0/mag_vec1);
        }
        else {
          dtheta = asin(dp/mag_vec0/mag_vec1) - M_PI;
        }
      }
      theta += dtheta;
      // check
      //cout << "ip, dtheta = " << ip << " " << dtheta*180.0/M_PI << endl;
    }
    // check
    if (mpi_rank == 0) cout << "theta = " << theta << endl;
    if (theta>0.0) {
      return -1; // right-hand rule, clockwise
    }
    else {
      return 1;  // right-hand rule, counter-clockwise
    }
  }

  void resampleSpline(SplineStuff::CubicSpline &cspline, const int nn, const double dt, const double dn, const double fraction, int &nt, double * &t_local, double (* &xt)[3]) {
    // borrowed from makeExtrudeSplineWithStrand

    // spline length can now be computed...
    const double L = cspline.getLength();
    const double sf = pow(dt/dn,1.0/(double(nn)-1.0)); assert(sf > 1.0);

    // discretize the spline at the smallest resolution: dn and
    // use the standoff distance to estimate the local spacing allowed...
    assert(dn <= dt);
    int np = L/dn;
    double (*xp)[3] = new double[np][3];
    for (int ip = 0; ip < np; ++ip) {
      const double t = double(ip)/double(np);
      cspline.getX(xp[ip],t);
    }
    //const double fraction = 0.1;
    double *dxp = new double[np];
    for (int ip = 0; ip < np; ++ip) {
      const int ip_prev = (ip-1+np)%np;
      const int ip_next = (ip+1)%np;
      // get the stand-off distance...
      double dx[3];
      FOR_I3 dx[i] = xp[ip][i] - 0.5*(xp[ip_prev][i]+xp[ip_next][i]);
      const double standoff = MAG(dx);
      // standoff should vary quadratically with spacing, so...
      if (standoff > fraction*dn) {
        // too much curvature, so just use dn...
        dxp[ip] = dn;
      }
      else {
        const double root = fraction*dn*standoff*(L*L/(double(np)*double(np)) + standoff*standoff - fraction*dn*standoff);
        assert(root > 0.0);
        dxp[ip] = min(dt,sqrt(root)/standoff);
      }
    }
    delete[] xp; xp = NULL;

    // now reduce dxp everywhere based on a stretching ratio...
    bool done = false;
    while (!done) {
      done = true;
      for (int ip = 0; ip < np; ++ip) {
        const int ip_prev = (ip-1+np)%np;
        const int ip_next = (ip+1)%np;
        const double dxp_nbr_min = min(dxp[ip_prev],dxp[ip_next]);
        if (dxp[ip] > dxp_nbr_min*pow(sf,0.5*dn/dxp_nbr_min)) {
          dxp[ip] = dxp_nbr_min*pow(sf,0.5*dn/dxp_nbr_min);
          done = false;
        }
      }
    }

    double factor = 1.0;
    done = false;
    while (!done) {
      double dist = 0.0;
      nt = 1; // first point at start of spline...
      for (int ip = 0; ip < np; ++ip) {
        const int ip_next = (ip+1)%np;
        //in the interval ip to ip_next, we are looking for dist == dx...
        // dist = dist0 + alpha*L/np
        // and assume dx varies linearly in the interval...
        // dx = alpha*dx1 + (1-alpha)*dx0...
        const double alpha = (dist - dxp[ip])/(dxp[ip_next]-dxp[ip]-L/double(np)*factor);
        if (alpha <= 1.0) {
          if (alpha <= 0.0) {
            if (mpi_rank == 0) 
              cout << "WARNING: ip = " << ip << ", alpha = " << alpha << endl;
          }
          //assert(alpha >= 0.0);
          ++nt;
          // reset the distance to the remaining distance in the interval...
          dist = (1.0-alpha)*L/double(np)*factor;
        }
        else {
          // add the whole interval distance to dist...
          dist += L/double(np)*factor;
        }
      }
      if (dist < 0.5*dxp[0]) {
        factor *= (L-dist)/L;
      }
      else {
        factor *= (L-dist)/(L-dxp[0]);
      }
      if (dist < 1.0E-10*dxp[0]) {
        // this means we have added one too many points (i.e. first point was repeated)..
        --nt;
        done = true;
      }
      else if (dxp[0]-dist < 1.0E-10*dxp[0]) {
        // this means we have added with the last dist almost exactly equal to the required dxp[0]...
        done = true;
      }
    }

    if (mpi_rank == 0) cout << " > tangential spacing nt: " << nt << endl;

    // at this point, we build t_local...

    if (t_local != NULL) delete[] t_local;
    t_local = new double[nt];
    double dist = 0.0;
    t_local[0] = 0.0;
    int it = 1; // first point at start of spline...
    for (int ip = 0; ip < np; ++ip) {
      const int ip_next = (ip+1)%np;
      //in the interval ip to ip_next, we are looking for dist == dx...
      // dist = dist0 + alpha*L/np
      // and assume dx varies linearly in the interval...
      // dx = alpha*dx1 + (1-alpha)*dx0...
      const double alpha = (dist - dxp[ip])/(dxp[ip_next]-dxp[ip]-L/double(np)*factor);
      if (alpha <= 1.0) {
        if (alpha <= 0.0) {
          if (mpi_rank == 0)
            cout << "WARNING: ip = " << ip << ", alpha = " << alpha << endl;
        }
        //assert(alpha >= 0.0);
        t_local[it++] = (double(ip)+alpha)/double(np);
        if (it == nt)
          break;
        dist = (1.0-alpha)*L/double(np)*factor;
      }
      else {
        // add the whole interval distance to dist...
        dist += L/double(np)*factor;
      }
    }
    delete[] dxp; dxp = NULL;

    // smooth t_local if requested by the user...

    const int nsmooth1 = 5;
    double * t_copy = new double[nt];
    for (int ismooth = 0; ismooth < nsmooth1; ++ismooth) {
      for (int it = 0; it < nt; ++it)
        t_copy[it] = t_local[it];
      // we want to keep the first point at zero, so its
      // adjustment is applied to all other points...
      const double dt0 = 0.5*(t_copy[nt-1]-1.0+t_copy[1]);
      for (int it = 1; it < nt; ++it) {
        if (it == nt-1) {
          t_local[nt-1] = 0.5*(t_copy[0]+1.0+t_copy[nt-2]) - dt0;
        }
        else {
          t_local[it] = 0.5*(t_copy[it-1]+t_copy[it+1]) - dt0;
        }
      }
    }
    delete[] t_copy;

    // actually set the points...
    if (xt != NULL) delete[] xt;
    xt = new double[nt][3];
    //double (*dxt)[3] = new double[nt][3];
    for (int it = 0; it < nt; ++it) {
      cspline.getX(xt[it],t_local[it]);
    }
  }

  // make surface voronoi point strand...
  bool makeSvpStrand(Part *part,Param *param,int &iarg) {
    srand(2+mpi_rank);

    int ierr = 0;
    set<int> zone_indices;
    int nsmooth = 0;
    double dn = 0.0; // min normal spacing
    double dt = 0.0; // tangential spacing (max normal spacing)
    int nn = 0; // max count in normal direction
    double crease_angle_degrees = 120.0; // increase this if you want to seed more points near large creases
    double area_ratio = 3.0; // increase this if you want to suppress "culling" of strands from divergence
    double seed_factor = 1.0; // increase this if you want less points
    double kappa_limit = HUGE_VAL; 
    bool b_friendly = true;
    bool b_write_dist = false; 
    double taper_factor = 1.0; // decrease this if you want to shrink the strand growth rate from feature
    bool b_const_dt = false; // use wgt'd lloyd when false
    double axis[3];
    bool b_axis = false;

    while (iarg < param->size()) {
      const string token = MiscUtils::toUpperCase(param->getString(iarg++));
      if ((token == "ZONE")||(token == "ZONES")||(token == "ZONE_NAME")||(token == "ZONE_NAMES")) {
        // expect a comma-delimited list of names...
        const string zoneNamesCsv = param->getString(iarg++);
        vector<string> tmp;
        MiscUtils::splitCsv(tmp,zoneNamesCsv);
        for (int ii = 0; ii < tmp.size(); ++ii) {
          const int izone = getZoneIndex(tmp[ii]); // get the unique zoneIndex in the partVec
          if (izone >= 0) {
            zone_indices.insert(izone);
          }
          else {
            WUI(WARN,"SVP_STRAND ZONE not recognized: " << tmp[ii]);
            ierr = -1;
          }
        }
      }
      else if (token == "DN") {
        dn = param->getDouble(iarg++);
      }
      else if (token == "DT") {
        dt = param->getDouble(iarg++); 
      }
      else if (token == "N") {
        nn = param->getInt(iarg++);
      }
      else if (token == "NSMOOTH") {
	nsmooth = param->getInt(iarg++);
      }
      else if (token == "CREASE_ANGLE") {
        crease_angle_degrees = param->getDouble(iarg++);
      }
      else if (token == "AREA_RATIO") {
        area_ratio = param->getDouble(iarg++);
      }
      else if (token == "SEED_FACTOR") {
        seed_factor = param->getDouble(iarg++);
      }
      else if (token == "TAPER_FACTOR") {
        taper_factor = param->getDouble(iarg++);
      }
      else if (token == "CRITICAL_RADIUS") {
        kappa_limit = 1.0/param->getDouble(iarg++);
      }
      else if (token == "CURVATURE_LIMIT") {
        kappa_limit = param->getDouble(iarg++);
      }
      else if (token == "FORCE") {
        b_friendly = false;
      }
      else if (token == "WRITE_DISTANCE") {
        b_write_dist = true;
      }
      else if (token == "CONST_DT") {
        b_const_dt = true;
      }
      else if (token == "AXIS") {
        b_axis = true;
        FOR_I3 axis[i] = param->getDouble(iarg++);
        const double axis_mag = MAG(axis);
        assert(axis_mag > 0.0);
        FOR_I3 axis[i] = axis[i]/axis_mag;
      }
      else if (token == "NLAYERS") {
	part->ff_nlayersString = param->getString(iarg++);
      }
    }

    if (zone_indices.empty()) {
      if (mpi_rank == 0) cout << "SVP_STRAND: missing ZONES <comma-delimited-part:zone_names>" << endl;
      ierr = -1;
    }
    if ((dt <= 0.0)||(dn >= dt)) {
      if (mpi_rank == 0) cout << "SVP_STRAND: missing/bad DT <double greater than DN>" << endl;
      ierr = -1;
    }
    if ((dn <= 0.0)||(dn >= dt)) {
      if (mpi_rank == 0) cout << "SVP_STRAND: missing/bad DN <double less than DT>" << endl;
      ierr = -1;
    }
    if (!(nn > 1)) {
      if (mpi_rank == 0) cout << "SVP_STRAND: missing/bad N <int greater than 1>" << endl;
      ierr = -1;
    }
    if (ierr != 0) {
      return false;
    }

    // default curvature limit:
    if (kappa_limit == HUGE_VAL)
      kappa_limit = 1.0/dt;

    // for isotropic top cells...

    const double ar = dt/dn; // aspect ratio
    const double sf = pow(ar,1.0/double(nn-1)); // size function
    const double dl = (dt*sf-dn)/(sf-1); // strand height

    if (mpi_rank == 0) {
      cout << " > strand min normal spacing: " << dn << ", max normal/tangential spacing: " << dt << ", aspect ratio: " << ar 
           << ", points: " << nn << ", height: " << dl << ", growth rate: " << 100.0*(sf-1.0) << "%" << endl;
    }

    vector<GroupData> groupDataVec;
    buildGroupsForFlaggedZones(groupDataVec,zone_indices);
    if (groupDataVec.size() != 1) {
      if (mpi_rank == 0)
        cout << " > supports a single connected region; break up part into a part per disjoint region." << endl;
      return false;
    }

    // record the associated parts
    part->b_link_surface = true;
    assert(part->ipart_link.empty());
    const int ipart = groupDataVec[0].ipart;
    part->ipart_link.push_back(ipart);
    assert(part->ist_offset_link.empty());
    part->ist_offset_link.push_back(0);
    assert(part->ff_type == UNKNOWN_FF);
    part->ff_type = FF_SURFACE_FAZONE_FF;
    assert(part->solid_type == UNKNOWN_SOLID);
    part->solid_type = NO_SOLID;

    SurfaceShm * surface = partVec[ipart]->surface;
    assert(surface);
    surface->ensureTeost(); // TODO make local

    int nvp = 0; 
    double (*xvp)[3] = NULL;
    double (*normal_vp)[3] = NULL;
    double *nnvp = NULL;
    double *r2vp = NULL;
    if (mpi_rank_internode_actual == 0) {

      if (mpi_rank == 0) 
        cout << "constructing flagged region..." << endl;

      int nst_split = surface->nst/mpi_size_shared;
      if (surface->nst%mpi_size_shared != 0) nst_split += 1;
      const int ist0 = min(surface->nst,nst_split*mpi_rank_shared);
      const int ist1 = min(surface->nst,nst_split*(mpi_rank_shared+1));

      int nft = 0; // number of flagged tris
      int* ftost = NULL;
      CTI_Mmap_rw(ftost,surface->nst);
      for (int ist = ist0; ist < ist1; ++ist)
        ftost[ist] = -1;
      MPI_Barrier(mpi_comm_shared);
      if (mpi_rank_shared == 0) {
        for (int ist = 0; ist < surface->nst; ++ist) {
          int igr;
          if (partVec[ipart]->getGroupForSt(igr,ist)) {
            if (igr == groupDataVec[0].igr) 
              ftost[ist] = nft++;
          }
        }
      }
      MPI_Bcast(&nft,1,MPI_INT,0,mpi_comm_shared);

      int *stoft = NULL;
      CTI_Mmap_rw(stoft,nft);
      for (int ist = ist0; ist < ist1; ++ist) {
        if (ftost[ist] >= 0) {
          assert(ftost[ist] < nft);
          stoft[ftost[ist]] = ist;
        }
        else {
          assert(ftost[ist] == -1);
        }
      }
      MPI_Barrier(mpi_comm_shared);

      int nft_split = nft/mpi_size_shared;
      if (nft%mpi_size_shared != 0) nft_split += 1;
      const int ift0 = min(nft,nft_split*mpi_rank_shared);
      const int ift1 = min(nft,nft_split*(mpi_rank_shared+1));

      int nsp_split = surface->nsp/mpi_size_shared;
      if (surface->nsp%mpi_size_shared != 0) nsp_split += 1;
      const int isp0 = min(surface->nsp,nsp_split*mpi_rank_shared);
      const int isp1 = min(surface->nsp,nsp_split*(mpi_rank_shared+1));

      int nfp = 0; // number of flagged points
      int* fposp = NULL;
      CTI_Mmap_rw(fposp,surface->nsp);
      for (int isp = isp0; isp < isp1; ++isp) 
        fposp[isp] = -1;
      MPI_Barrier(mpi_comm_shared);
      if (mpi_rank_shared == 0) {
        for (int ift = 0; ift < nft; ++ift) {
          const int ist = stoft[ift]; 
          FOR_I3 {
            const int isp = surface->spost[ist][i];
            if (fposp[isp] == -1)
              fposp[isp] = nfp++;
          }
        }
      }
      MPI_Bcast(&nfp,1,MPI_INT,0,mpi_comm_shared);

      int *spofp = NULL;
      CTI_Mmap_rw(spofp,nfp);
      for (int isp = isp0; isp < isp1; ++isp) {
        if (fposp[isp] >= 0) {
          assert(fposp[isp] < nfp);
          spofp[fposp[isp]] = isp;
        }
        else {
          assert(fposp[isp] == -1);
        }
      }
      MPI_Barrier(mpi_comm_shared);

      int nfp_split = nfp/mpi_size_shared;
      if (nfp%mpi_size_shared != 0) nfp_split += 1;
      const int ifp0 = min(nfp,nfp_split*mpi_rank_shared);
      const int ifp1 = min(nfp,nfp_split*(mpi_rank_shared+1));

      if (mpi_rank == 0) 
        cout << " > number of flagged tris: " << nft << ", number of flagged points: " << nfp << ", number of border points: " << groupDataVec[0].spoli_s << endl;

      if (mpi_rank == 0)
        cout << " > building stofp_i/v... " << endl;

      int *stofp_i = NULL;
      CTI_Mmap_rw(stofp_i,nfp+1);
      for (int ifp = ifp0; ifp < ifp1; ++ifp)
        stofp_i[ifp+1] = 0;
      MPI_Barrier(mpi_comm_shared);
      if (mpi_rank_shared == 0) {
        // go through all tris to get non-flagged tri's touching flagged nodes...
        for (int ist = 0; ist < surface->nst; ++ist) {
          FOR_I3 {
            const int isp = surface->spost[ist][i];
            const int ifp = fposp[isp]; 
            if (ifp != -1) {
              assert((ifp >= 0)&&(ifp < nfp));
              ++stofp_i[ifp+1];
            }
          }
        }
        stofp_i[0] = 0;
        for (int ifp = 0; ifp < nfp; ++ifp)
          stofp_i[ifp+1] += stofp_i[ifp];
      }
      MPI_Barrier(mpi_comm_shared);
      int *stofp_v = NULL;
      //cout << nfp << " " << stofp_i[nfp] << endl;
      CTI_Mmap(stofp_v,stofp_i[nfp]);
      if (mpi_rank_shared == 0) {
        for (int ist = 0; ist < surface->nst; ++ist) {
          FOR_I3 {
            const int isp = surface->spost[ist][i];
            const int ifp = fposp[isp]; 
            if (ifp != -1) {
              assert((ifp >= 0)&&(ifp < nfp));
              stofp_v[stofp_i[ifp]++] = ist;
            }
          }
        }
        for (int ifp = nfp-1; ifp > 0; --ifp)
          stofp_i[ifp] = stofp_i[ifp-1];
        stofp_i[0] = 0;
      }
      MPI_Barrier(mpi_comm_shared);

      // build normal on flagged points... 

      if (mpi_rank == 0) 
        cout << " > constructing normal_fp..." << endl;

      double (*normal_fp)[3] = NULL;
      CTI_Mmap_rw(normal_fp,nfp);
      for (int ifp = ifp0; ifp < ifp1; ++ifp) {
        const int isp = spofp[ifp];
        FOR_I3 normal_fp[ifp][i] = 0.0;
        for (int top = stofp_i[ifp]; top != stofp_i[ifp+1]; ++top) {
          const int ist = stofp_v[top];
          double dx[3][3];
          double dx_mag[3];
          FOR_I3 {
            FOR_J3 dx[i][j] = surface->xsp[surface->spost[ist][(i+1)%3]][j] - surface->xsp[surface->spost[ist][i]][j];
            dx_mag[i] = MAG(dx[i]);
          }
          FOR_I3 {
            if (surface->spost[ist][i] == isp) {
              const double cp[3] = CROSS_PRODUCT(dx[(i+2)%3],dx[i]);
              const double cp_mag = MAG(cp);
              if (cp_mag > 1.0E-6*dx_mag[(i+2)%3]*dx_mag[i]) {
                const double dp = -DOT_PRODUCT(dx[(i+2)%3],dx[i]);
                const double theta = atan2(cp_mag,dp);
                assert(theta > 0.0);
                FOR_J3 normal_fp[ifp][j] += theta*cp[j]/cp_mag;
              }
            }
          }
        }
        if (b_axis) {
          const double dp = DOT_PRODUCT(normal_fp[ifp],axis);
          FOR_I3 normal_fp[ifp][i] -= dp*axis[i];
        }
        const double mag = MAG(normal_fp[ifp]);
        assert(mag > 0.0);
        FOR_I3 normal_fp[ifp][i] /= mag; // unit normal
      }
      MPI_Barrier(mpi_comm_shared);

      if (mpi_rank == 0)
        cout << " > constructing dist_fp..." << endl;

      double *dist_fp = NULL;
      CTI_Mmap_rw(dist_fp,nfp);
      int *flag_fp = NULL;
      CTI_Mmap_rw(flag_fp,nfp);
      for (int ifp = ifp0; ifp < ifp1; ++ifp) {
        dist_fp[ifp] = HUGE_VAL;
        flag_fp[ifp] = 0;
        const int isp = spofp[ifp];
        for (int top = stofp_i[ifp]; top != stofp_i[ifp+1]; ++top) {
          const int ist = stofp_v[top];
          if (ftost[ist] >= 0) {
            FOR_I3 {
              if (surface->spost[ist][i] == isp) {
                int ist_nbr,i_nbr,bit;
                surface->unpackTeost(ist_nbr,i_nbr,bit,surface->teost[ist][i]);
                if (ftost[ist_nbr] >= 0) {
                  if (!checkCoplanar(ipart,ist,0,ipart,ist_nbr,0,crease_angle_degrees)) {
                    //cout << "crease: " << endl;
                    dist_fp[ifp] = 0.0;
                    flag_fp[ifp] = 1;
                  }
                }
              }
            }
          }
        }
      }
      MPI_Barrier(mpi_comm_shared);
      double my_kappa_max = 0.0;
      for (int ift = ift0; ift < ift1; ++ift) {
        const int ist = stoft[ift];
        const double normal2[3] = TRI_NORMAL_2(surface->xsp[surface->spost[ist][0]],surface->xsp[surface->spost[ist][1]],surface->xsp[surface->spost[ist][2]]);
        const double area2 = MAG(normal2);
        double xsp_dn[3][3];
        double xsp_mdn[3][3];
        FOR_I3 {
          const int isp = surface->spost[ist][i];
          const int ifp = fposp[isp]; assert((ifp >= 0)&&(ifp < nfp));
          FOR_J3 xsp_dn[i][j] = surface->xsp[isp][j]-normal_fp[ifp][j]*dn;
          FOR_J3 xsp_mdn[i][j] = surface->xsp[isp][j]+normal_fp[ifp][j]*dn;
        }
        const double normal2_dn[3] = TRI_NORMAL_2(xsp_dn[0],xsp_dn[1],xsp_dn[2]);
        const double area2_dn = MAG(normal2_dn);
        const double normal2_mdn[3] = TRI_NORMAL_2(xsp_mdn[0],xsp_mdn[1],xsp_mdn[2]);
        const double area2_mdn = MAG(normal2_mdn);
        const double kappa = max(fabs(area2-area2_dn),fabs(area2-area2_mdn))/(area2*dn);
        my_kappa_max = max(kappa,my_kappa_max);
        if (kappa > kappa_limit) {
          FOR_I3 {
            const int isp = surface->spost[ist][i];
            const int ifp = fposp[isp]; assert((ifp >= 0)&&(ifp < nfp));
            dist_fp[ifp] = 0.0;
            flag_fp[ifp] = 1;
          }
        }
      }
      double kappa_max;
      MPI_Reduce(&my_kappa_max,&kappa_max,1,MPI_DOUBLE,MPI_MAX,0,mpi_comm_shared);
      if (mpi_rank == 0) 
        cout << " > min radius: " << 1.0/kappa_max << ", critical radius: " << 1.0/kappa_limit << endl;
      if (mpi_rank_shared == 0) {
        MinHeap trialHeap(nfp);
        if (!b_axis) {
          for (int sol = 0; sol < groupDataVec[0].spoli_s; ++sol) {
            const int isp = groupDataVec[0].spoli_v[sol];
            const int ifp = fposp[isp]; assert((ifp >= 0)&&(ifp < nfp));
            dist_fp[ifp] = 0.0;
            flag_fp[ifp] = 1;
          }
        }
        for (int ifp = 0; ifp < nfp; ++ifp)
          if (flag_fp[ifp] == 1)
            trialHeap.insert(pair<double,int>(0.0,ifp));
        while (trialHeap.size() > 0)
          fmmLoop(trialHeap,flag_fp,dist_fp,stofp_i,stofp_v,fposp,surface);
      }
      MPI_Barrier(mpi_comm_shared);
      CTI_Munmap(flag_fp,nfp);

      for (int ifp = ifp0; ifp < ifp1; ++ifp) {
        dist_fp[ifp] = fmax(dn,fmin(dl+dn,fabs(dist_fp[ifp])));
        //dist_fp[ifp] = fmax(0.0,fmin(dl,fabs(dist_fp[ifp])));
        //cout << dist_fp[ifp] << endl;
      }
      MPI_Barrier(mpi_comm_shared);

      if (b_write_dist&&(mpi_rank == 0)) {
        FILE * fp = fopen("dist.dat","w");
        for (int ifp = 0; ifp < nfp; ++ifp) {
          const int isp = spofp[ifp];
          fprintf(fp,"%18.15le %18.15le %18.15le %18.15le\n",surface->xsp[isp][0],surface->xsp[isp][1],surface->xsp[isp][2],dist_fp[ifp]);
        }
        fclose(fp);
      }

      if (mpi_rank == 0) 
        cout << " > seeding initial points..." << endl;

      // need to get number of samples using tri area and divergence of point normals...

      //const double area_vor = seed_factor*dt*dt;
      const double area_vor = seed_factor*0.25*M_PI*dt*dt;
      int my_nvp = 0;
      int *vpoft_i = NULL;
      CTI_Mmap_rw(vpoft_i,nft+1);
      double area_left = 0.0;
      for (int ift = ift0; ift < ift1; ++ift) {
        const int ist = stoft[ift];
        const double normal2[3] = TRI_NORMAL_2(surface->xsp[surface->spost[ist][0]],surface->xsp[surface->spost[ist][1]],surface->xsp[surface->spost[ist][2]]);
        const double area = 0.5*MAG(normal2);
        //double xsp_half_dl[3][3];
        //double xsp_half_mdl[3][3];
        double dist_ft = 0.0;
        FOR_I3 {
          const int isp = surface->spost[ist][i];
          const int ifp = fposp[isp];
          //FOR_J3 xsp_half_dl[i][j] = surface->xsp[isp][j]-normal_fp[ifp][j]*0.5*dl;
          //FOR_J3 xsp_half_mdl[i][j] = surface->xsp[isp][j]+normal_fp[ifp][j]*0.5*dl;
          dist_ft += dist_fp[ifp];
        }
        dist_ft /= 3.0;
        //const double normal2_half_dl[3] = TRI_NORMAL_2(xsp_half_dl[0],xsp_half_dl[1],xsp_half_dl[2]);
        //const double area_half_dl = 0.5*MAG(normal2_half_dl);
        //const double normal2_half_mdl[3] = TRI_NORMAL_2(xsp_half_mdl[0],xsp_half_mdl[1],xsp_half_mdl[2]);
        //const double area_half_mdl = 0.5*MAG(normal2_half_mdl);
        // linear
        //const double delta = dn+(dt-dn)*dist_ft/dl;
        // geometric
        //const double n = max(1.0,log(dist_ft/dn*(sf-1.0)+1.0)/log(sf));
        //const double delta = dn*pow(sf,n-1.0);
        //const double wgt = dt*dt/(delta*delta);
        //const double wgt = dt/delta;
        //const double eff_area = MAX3(area,area_half_dl,area_half_mdl)*wgt+area_left; 
        double wgt;
        if (b_const_dt)
          wgt = 1.0;
        else
          wgt = (dl+dn)/dist_ft;
        const double eff_area = area*wgt+area_left; 
        const int nsamples = int(eff_area/area_vor);
        area_left = eff_area-nsamples*area_vor; 
        //cout << nsamples << " " << eff_area << " " << area_vor << " " << area_left << " " << wgt << endl;
        vpoft_i[ift+1] = nsamples;
        my_nvp += nsamples;
      }
      MPI_Allreduce(&my_nvp,&nvp,1,MPI_INT,MPI_SUM,mpi_comm_shared);
      if (mpi_rank == 0)
        cout << " > nvp: " << nvp << endl;

      if (mpi_rank_shared == 0) {
        vpoft_i[0] = 0;
        for (int ift = 0; ift < nft; ++ift) 
          vpoft_i[ift+1] += vpoft_i[ift];
        assert(vpoft_i[nft] == nvp);
      }
      MPI_Barrier(mpi_comm_shared);

      // allocate some helper arrays in shared memory as well...

      int *vpoft_v = NULL;
      CTI_Mmap_rw(vpoft_v,nvp);
      int *ftovp = NULL;  
      CTI_Mmap_rw(ftovp,nvp); 
      assert(r2vp == NULL); CTI_Mmap_rw(r2vp,nvp);
      assert(xvp == NULL); CTI_Mmap_rw(xvp,nvp);
      assert(normal_vp == NULL); CTI_Mmap_rw(normal_vp,nvp);
      double (*dxvp)[3] = normal_vp; // reuse

      // initialize vp data...

      for (int ift = ift0; ift < ift1; ++ift) {
        const int ist = stoft[ift];
        for (int ivp = vpoft_i[ift]; ivp < vpoft_i[ift+1]; ++ivp) {
          const double r0 = 0.999*double(rand())/double(RAND_MAX)+0.0005;
          const double r1 = 0.999*double(rand())/double(RAND_MAX)+0.0005;
          FOR_I3 xvp[ivp][i] = (1.0-sqrt(r0))*surface->xsp[surface->spost[ist][0]][i] +
            sqrt(r0)*(1.0-r1)*surface->xsp[surface->spost[ist][1]][i] +
            sqrt(r0)*r1*surface->xsp[surface->spost[ist][2]][i];
          r2vp[ivp] = dt*dt;
          ftovp[ivp] = ift;
        }
      }
      MPI_Barrier(mpi_comm_shared);
      int nvp_split = nvp/mpi_size_shared;
      if (nvp%mpi_size_shared != 0) nvp_split += 1;
      const int ivp0 = min(nvp,nvp_split*mpi_rank_shared);
      const int ivp1 = min(nvp,nvp_split*(mpi_rank_shared+1));

      CuttableVoronoiData cvd;
      stack<int> ftStack;
      set<pair<double,int> > nbrSet;
      map<const pair<int,int>,int> edgeMap; // can we get around this?

      // are these allocations too large?
      int* flag_st = new int[surface->nst];
      int* flag_sp = new int[surface->nsp]; 
      // no need to set them all...
      for (int ift = 0; ift < nft; ++ift) {
        const int ist = stoft[ift];
        flag_st[ist] = -1;
        FOR_I3 {
          int ist_nbr,i_nbr,bit;
          surface->unpackTeost(ist_nbr,i_nbr,bit,surface->teost[ist][i]);
          flag_st[ist_nbr] = -1;
        }
        FOR_I3 flag_sp[surface->spost[ist][i]] = -1;
      }
      vector<int> ftVec; // for cleaning up flag_st and flag_sp

      if (mpi_rank == 0)
        cout << " > smoothing points..." << endl; 

      double *dist_vp = NULL;
      CTI_Mmap_rw(dist_vp,nvp);
      vector<int>* nbrVec_vp = new vector<int>[ivp1-ivp0];
      set<int> faSet; // to build nbrVec_vp
      map<const int,double*> x0Map;
      assert(nnvp == NULL); CTI_Mmap_rw(nnvp,nvp);

      int ismooth = 0;
      int done = 0;
      while (done != 2) {
        ++ismooth;

        // rebuild vpoft_i/v...
        if (mpi_rank_shared == 0) {
          if (ismooth > 1) {
            for (int ift = 0; ift < nft; ++ift)
              vpoft_i[ift+1] = 0;
            for (int ivp = 0; ivp < nvp; ++ivp) 
              ++vpoft_i[ftovp[ivp]+1];
            assert(vpoft_i[0] == 0);
            for (int ift = 0; ift < nft; ++ift) 
              vpoft_i[ift+1] += vpoft_i[ift];
            assert(vpoft_i[nft] == nvp);
          }
          for (int ivp = 0; ivp < nvp; ++ivp) 
            vpoft_v[vpoft_i[ftovp[ivp]]++] = ivp;
          for (int ift = nft-1; ift > 0; --ift) 
            vpoft_i[ift] = vpoft_i[ift-1];
          vpoft_i[0] = 0;
        }
        MPI_Barrier(mpi_comm_shared);

        // for stats on Lloyd iterations

        double my_dx2_max = 0.0;
        double my_dx2_sum = 0.0;
        int my_count[3] = {0,0,0};

        // calculate weights...
        for (int ivp = ivp0; ivp < ivp1; ++ivp) {
          const int ift = ftovp[ivp]; 
          const int ist = stoft[ift];
          double n_sum = 0.0;
          dist_vp[ivp] = 0.0;
          FOR_I3 {
            const int isp = surface->spost[ist][i];
            const int ifp = fposp[isp]; assert((ifp >= 0)&&(ifp < nfp));
            const double n[3] = TRI_NORMAL_2(xvp[ivp],surface->xsp[surface->spost[ist][(i+1)%3]],surface->xsp[surface->spost[ist][(i+2)%3]]);
            const double n_mag = MAG(n);
            dist_vp[ivp] += n_mag*dist_fp[ifp]; 
            n_sum += n_mag;
          }
          if (n_sum > 0.0)
            dist_vp[ivp] /= n_sum;

        }
        MPI_Barrier(mpi_comm_shared);

        //if (mpi_rank_shared == 0) 
        //  for (int ivp = 0; ivp < nvp; ++ivp)
        //    cout << "dist_vp " << dist_vp[ivp] << " " << ivp << endl;
        //MPI_Pause("ok");

        // now move each ivp...
        for (int ivp = ivp0; ivp < ivp1; ++ivp) {
          cvd.clear();
          nbrSet.clear();
          ftVec.clear(); 
          edgeMap.clear();

          flag_st[stoft[ftovp[ivp]]] = 1;
          ftStack.push(ftovp[ivp]);
          while (!ftStack.empty()) {

            // pop the ift...
            const int ift = ftStack.top(); ftStack.pop();
            const int ist = stoft[ift];
            assert(flag_st[ist] == 1);
            ftVec.push_back(ift);

            // add all nbrs in the same ift -- even the first time...
            for (int vof = vpoft_i[ift]; vof != vpoft_i[ift+1]; ++vof) {
              const int ivp_nbr = vpoft_v[vof];
              if (ivp != ivp_nbr) {
                const double d2 = DIST2(xvp[ivp_nbr],xvp[ivp]);
                if (d2 < r2vp[ivp]) 
                  nbrSet.insert(pair<double,int>(d2,ivp_nbr)); // use a set to sort
              }
            }

            // loop the nodes ccw...
            FOR_I3 {
              const int isp = surface->spost[ist][i];
              if (flag_sp[isp] == -1) {
                const int ino = cvd.new_node();
                flag_sp[isp] = ino;
                FOR_I3 cvd.x_no[ino][i] = surface->xsp[isp][i] - xvp[ivp][i];
              }
            }

            // loop on the edges ccw. By definition, these edges own the i and (i+1)%3 nodes...
            FOR_I3 {
              int ist_nbr,i_nbr,bit;
              surface->unpackTeost(ist_nbr,i_nbr,bit,surface->teost[ist][i]);
              const int ino0 = flag_sp[surface->spost[ist][i]]; assert((ino0 >= 0)&&(ino0 < cvd.nno));
              const int ino1 = flag_sp[surface->spost[ist][(i+1)%3]]; assert((ino1 >= 0)&&(ino1 < cvd.nno));
              map<const pair<int,int>,int>::iterator iter = edgeMap.find(pair<int,int>(ino0,ino1));
              if (iter == edgeMap.end()) {
                const int ied = cvd.new_edge();
                edgeMap[pair<int,int>(ino1,ino0)] = ied;
                cvd.nooed[ied][0] = ino0;
                cvd.nooed[ied][1] = ino1;
                cvd.faoed[ied][0] = ist;
                if (ftost[ist_nbr] == -1) 
                  cvd.faoed[ied][1] = -2; // boundary
                else 
                  cvd.faoed[ied][1] = -1; // look for valid nbr
              }
              else {
                const int ied = iter->second;
                assert(cvd.nooed[ied][0] == ino1);
                assert(cvd.nooed[ied][1] == ino0);
                assert(cvd.faoed[ied][0] != ist); // no edge can have the same tri on both sides
                assert(cvd.faoed[ied][1] == -1);
                cvd.faoed[ied][1] = ist;
                edgeMap.erase(iter);
              }

              // add neighboring tri to stack if partially within radius;
              // only add neighbors if they are part of selected facets
              if ((ftost[ist_nbr] >= 0)&&(flag_st[ist_nbr] == -1)) {
                if (MiscUtils::getPointToEdgeDist2(xvp[ivp],surface->xsp[surface->spost[ist][i]],surface->xsp[surface->spost[ist][(i+1)%3]]) < r2vp[ivp]) {
                  flag_st[ist_nbr] = 1;
                  ftStack.push(ftost[ist_nbr]);
                }
              }
            }
          }

          // reset flags...

          for (int jj = 0,lim2 = ftVec.size(); jj < lim2; ++jj) {
            const int ist = stoft[ftVec[jj]];
            assert(flag_st[ist] == 1);
            flag_st[ist] = -1;
            FOR_I3 {
              const int isp = surface->spost[ist][i];
              flag_sp[isp] = -1;
            }
          }

          cvd.setD2Max();

          // now cut against nbrs...
          for (set<pair<double,int> >::iterator iter = nbrSet.begin(); iter != nbrSet.end(); ++iter) {
            if (iter->first >= 4.0*cvd.d2_max) // factor of 2^2 = 4 here
              break;
            const int ivp_nbr = iter->second;
            const double dn_nbr[3] = DIFF(xvp[ivp_nbr],xvp[ivp]);
            if ((2.0*dn_nbr[0] < cvd.Lmax[0])&&
                (2.0*dn_nbr[1] < cvd.Lmax[1])&&
                (2.0*dn_nbr[2] < cvd.Lmax[2])&&
                (2.0*dn_nbr[0] > cvd.Lmin[0])&&
                (2.0*dn_nbr[1] > cvd.Lmin[1])&&
                (2.0*dn_nbr[2] > cvd.Lmin[2]) )
              cvd.cut_surf(dn_nbr,-ivp_nbr-8); // updates d2_max...
          }

          // calculate the centroid of everything...
          // note that we also populate stVec with JUST the tris remaining in the
          // cut surface here...

          double dxc[3] = { 0.0, 0.0, 0.0 };
          bool has_seed_boundary = false;
          double area2_sum = 0.0;
          x0Map.clear();
          for (int ied = 0; ied < cvd.ned; ++ied) {
            if (cvd.faoed[ied][0] >= 0) {
              map<const int,double*>::const_iterator it = x0Map.find(cvd.faoed[ied][0]);
              if (it == x0Map.end()) {
                // on the first time we visit any face, store a node to use in all
                // the tris with other edges (the face is planar and convex, so this is fine)...
                x0Map[cvd.faoed[ied][0]] = cvd.x_no[cvd.nooed[ied][0]];
              }
              else if (cvd.x_no[cvd.nooed[ied][1]] != it->second) {
                const double normal2[3] = TRI_NORMAL_2(it->second,cvd.x_no[cvd.nooed[ied][0]],cvd.x_no[cvd.nooed[ied][1]]);
                double area2 = MAG(normal2);
                if (!b_const_dt) {
                  const int ivp_nbr = -cvd.faoed[ied][1]-8;
                  if (ivp_nbr >= 0) { 
                    area2 *= pow((dist_vp[ivp]*(sf-1.0)+dn)/(dist_vp[ivp_nbr]*(sf-1.0)+dn),2);
                    // linear
                    //const double delta = dn+(dt-dn)*dist_vp[ivp_nbr]/dl;
                    // geometric
                    //const double n = max(1.0,log(dist_vp[ivp_nbr]/dn*(sf-1.0)+1.0)/log(sf));
                    //const double delta = dn*pow(sf,n-1.0);
                    //area2 *= dn*dn/(delta*delta);
                  }
                }
                FOR_I3 dxc[i] += area2*(it->second[i]+cvd.x_no[cvd.nooed[ied][0]][i]+cvd.x_no[cvd.nooed[ied][1]][i]);
                area2_sum += area2;
              }
            }
            if (cvd.faoed[ied][1] >= 0) {
              map<const int,double*>::const_iterator it = x0Map.find(cvd.faoed[ied][1]);
              if (it == x0Map.end()) {
                x0Map[cvd.faoed[ied][1]] = cvd.x_no[cvd.nooed[ied][1]];
              }
              else if (cvd.x_no[cvd.nooed[ied][0]] != it->second) {
                const double normal2[3] = TRI_NORMAL_2(it->second,cvd.x_no[cvd.nooed[ied][1]],cvd.x_no[cvd.nooed[ied][0]]);
                double area2 = MAG(normal2);
                if (!b_const_dt) {
                  const int ivp_nbr = -cvd.faoed[ied][0]-8;
                  if (ivp_nbr >= 0) {
                    area2 *= pow((dist_vp[ivp]*(sf-1.0)+dn)/(dist_vp[ivp_nbr]*(sf-1.0)+dn),2);
                    // linear
                    //const double delta = dn+(dt-dn)*dist_vp[ivp_nbr]/dl;
                    // geometric
                    //const double n = max(1.0,log(dist_vp[ivp_nbr]/dn*(sf-1.0)+1.0)/log(sf));
                    //const double delta = dn*pow(sf,n-1.0);
                    //area2 *= dn*dn/(delta*delta);
                  }
                }
                FOR_I3 dxc[i] += area2*(it->second[i]+cvd.x_no[cvd.nooed[ied][0]][i]+cvd.x_no[cvd.nooed[ied][1]][i]);
                area2_sum += area2;
              }
            }
            else if (cvd.faoed[ied][1] == -1) {
              has_seed_boundary = true;
              break;
            }
            else {
              // must be a -2: this is ok (flagged region boundary)...
              assert(cvd.faoed[ied][1] == -2);
              //has_boundary = true;
            }
          }

          if (done == 1)
            nnvp[ivp] = double(nn)+0.0001; // add eps so int truncation gives you right number

          if (has_seed_boundary) {
            // we are too small. Increase deltap[ip], and do NOT apply any motion to the point. The remaining
            // patch is arbitrarily large...
            ++my_count[0];
            //assert(done == 0); // b_friendly mode
            r2vp[ivp] *= 2.25;
            FOR_I3 dxvp[ivp][i] = 0.0;
          }
          else {
            // Even though we have cut away all the original seed boundary of this patch, there still
            // may be nbrs out there that can cut some of our corners off if deltap was not large
            // enough...
            if (r2vp[ivp] <= 4.0*cvd.d2_max) {
              // other neighbors might exist that could still cut this volume. Find them in
              // the next iteration...
              ++my_count[1];
              //assert(done == 0); // can get here after done = 1 because we reproject
              r2vp[ivp] = 6.25*cvd.d2_max; // could be just 4.0 here? - use a little larger so face checks pass for sure
              FOR_I3 dxvp[ivp][i] = 0.0;
            }
            else {
              ++my_count[2];
              if (done == 0) {
                // regular lloyd iteration pass
                assert(area2_sum > 0.0);
                r2vp[ivp] = 6.25*cvd.d2_max;
                FOR_I3 dxc[i] /= 3.0*area2_sum;

                // project new point onto closest tri...
                int ist_closest = -1;
                double d2_closest;
                for (int ied = 0; ied < cvd.ned; ++ied) {
                  if (cvd.faoed[ied][0] >= 0) {
                    map<const int,double*>::const_iterator it = x0Map.find(cvd.faoed[ied][0]);
                    assert(it != x0Map.end());
                    if ((cvd.x_no[cvd.nooed[ied][0]] != it->second)&&(cvd.x_no[cvd.nooed[ied][1]] != it->second)) {
                      double this_dxp[3];
                      MiscUtils::getClosestPointOnTriInterior(this_dxp,dxc,it->second,cvd.x_no[cvd.nooed[ied][0]],cvd.x_no[cvd.nooed[ied][1]]);
                      const double this_d2 = DIST2(this_dxp,dxc);
                      if ((ist_closest == -1)||(this_d2 < d2_closest)) {
                        FOR_I3 dxvp[ivp][i] = this_dxp[i];
                        d2_closest = this_d2;
                        ist_closest = cvd.faoed[ied][0];
                      }
                    }
                  }
                  if (cvd.faoed[ied][1] >= 0) {
                    map<const int,double*>::const_iterator it = x0Map.find(cvd.faoed[ied][1]);
                    assert(it != x0Map.end());
                    if ((cvd.x_no[cvd.nooed[ied][0]] != it->second)&&(cvd.x_no[cvd.nooed[ied][1]] != it->second)) {
                      double this_dxp[3];
                      MiscUtils::getClosestPointOnTriInterior(this_dxp,dxc,it->second,cvd.x_no[cvd.nooed[ied][0]],cvd.x_no[cvd.nooed[ied][1]]);
                      const double this_d2 = DIST2(this_dxp,dxc);
                      if ((ist_closest == -1)||(this_d2 < d2_closest)) {
                        FOR_I3 dxvp[ivp][i] = this_dxp[i];
                        d2_closest = this_d2;
                        ist_closest = cvd.faoed[ied][1];
                      }
                    }
                  }
                }
                // make sure we found one!...
                assert(ist_closest != -1);
                ftovp[ivp] = ftost[ist_closest];
              }
              else {
                // repeat to get neighbors...
                FOR_I3 dxvp[ivp][i] = 0.0;

                // store nbrs...
                // nbrVec_vp holds nbr each internal edge...

                faSet.clear();
                for (int ied = 0; ied < cvd.ned; ++ied) {
                  FOR_I2 {
                    if (cvd.faoed[ied][i] <= -8) 
                      faSet.insert(-cvd.faoed[ied][i]-8);
                  }
                }
                nbrVec_vp[ivp-ivp0].clear();
                for (set<int>::iterator iter = faSet.begin(); iter != faSet.end(); ++iter) 
                  nbrVec_vp[ivp-ivp0].push_back(*iter);

                for (int ied = 0; ied < cvd.ned; ++ied) {
                  if (!b_axis) {
                    if ((cvd.faoed[ied][0] >= 0)&&(cvd.faoed[ied][1] == -2)) {
                      nnvp[ivp] = 1.0001;
                      break;
                    }
                  }
                }
              }
              const double this_dx2 = DOT_PRODUCT(dxvp[ivp],dxvp[ivp])/r2vp[ivp];
              my_dx2_max = max(my_dx2_max,this_dx2);
              my_dx2_sum += this_dx2;
            }
          }
        }

        for (int ivp = ivp0; ivp < ivp1; ++ivp) 
          FOR_I3 xvp[ivp][i] += dxvp[ivp][i];

        double dx2_max,dx2_sum;
        MPI_Reduce(&my_dx2_max,&dx2_max,1,MPI_DOUBLE,MPI_MAX,0,mpi_comm_shared);
        MPI_Reduce(&my_dx2_sum,&dx2_sum,1,MPI_DOUBLE,MPI_SUM,0,mpi_comm_shared);

        int count[3];
        MPI_Reduce(my_count,count,3,MPI_INT,MPI_SUM,0,mpi_comm_shared);

        if (mpi_rank == 0)
          cout << "iter: " << ismooth << " counts: " << count[0] << " (has seed), " << count[1] << " (grow delta), " << count[2] << " (passed); l2,linf: " << sqrt(dx2_sum/double(nvp)) << " " << sqrt(dx2_max) << endl;

        if (mpi_rank_shared == 0) {
          if ((done == 0)&&(b_friendly||((count[0]+count[1])==0))&&(ismooth >= nsmooth)) {
            // we need to go through again and compute the nbrs...
            done = 1;
          }
          else if (done == 1) {
            // we are done...
            done = 2;
          }
        }

        MPI_Bcast(&done,1,MPI_INT,0,mpi_comm_shared);

      }
      delete[] flag_st;
      delete[] flag_sp;
      CTI_Munmap(dist_fp,nfp);

      // get neighbors to build triangulation...

      int* vpovp_i = NULL; 
      CTI_Mmap_rw(vpovp_i,nvp+1);
      for (int ivp = ivp0; ivp < ivp1; ++ivp)
        vpovp_i[ivp+1] = nbrVec_vp[ivp-ivp0].size();
      MPI_Barrier(mpi_comm_shared);
      if (mpi_rank_shared == 0) {
        vpovp_i[0] = 0;
        for (int ivp = 0; ivp < nvp; ++ivp)
          vpovp_i[ivp+1] += vpovp_i[ivp];
      }
      MPI_Barrier(mpi_comm_shared);
      int* vpovp_v = NULL;
      CTI_Mmap_rw(vpovp_v,vpovp_i[nvp]);
      for (int ivp = ivp0; ivp < ivp1; ++ivp) {
        const int nnbr = nbrVec_vp[ivp-ivp0].size();
        assert(nnbr == vpovp_i[ivp+1]-vpovp_i[ivp]);
        for (int inbr = 0; inbr < nnbr; ++inbr) {
          const int ivp_nbr = nbrVec_vp[ivp-ivp0][inbr];
          vpovp_v[vpovp_i[ivp]+inbr] = ivp_nbr;
        }
      }
      delete[] nbrVec_vp;
      MPI_Barrier(mpi_comm_shared);

      if (mpi_rank == 0) 
        cout << " > determining number of points in each strand..." << endl;

      for (int ivp = ivp0; ivp < ivp1; ++ivp) {
        double xx = 0.0,xy = 0.0,xz = 0.0;
        double yy = 0.0,yz = 0.0,zz = 0.0;
        double xc[3] = {xvp[ivp][0],xvp[ivp][1],xvp[ivp][2]};
        for (int vov = vpovp_i[ivp]; vov != vpovp_i[ivp+1]; ++vov) {
          const int ivp_nbr = vpovp_v[vov];
          FOR_I3 xc[i] += xvp[ivp_nbr][i];
        }
        FOR_I3 xc[i] /= double(vpovp_i[ivp+1]-vpovp_i[ivp]+1);
        const double dx[3] = DIFF(xvp[ivp],xc);
        const double wgt = exp(-DOT_PRODUCT(dx,dx)/r2vp[ivp]);
        xx += wgt*dx[0]*dx[0];
        xy += wgt*dx[0]*dx[1];
        xz += wgt*dx[0]*dx[2];
        yy += wgt*dx[1]*dx[1];
        yz += wgt*dx[1]*dx[2];
        zz += wgt*dx[2]*dx[2];
        for (int vov = vpovp_i[ivp]; vov != vpovp_i[ivp+1]; ++vov) {
          const int ivp_nbr = vpovp_v[vov];
          const double dx_nbr[3] = DIFF(xvp[ivp_nbr],xc);
          const double wgt = exp(-DOT_PRODUCT(dx_nbr,dx_nbr)/r2vp[ivp]);
          xx += wgt*dx_nbr[0]*dx_nbr[0];
          xy += wgt*dx_nbr[0]*dx_nbr[1];
          xz += wgt*dx_nbr[0]*dx_nbr[2];
          yy += wgt*dx_nbr[1]*dx_nbr[1];
          yz += wgt*dx_nbr[1]*dx_nbr[2];
          zz += wgt*dx_nbr[2]*dx_nbr[2];
        }
        const double det_x = yy*zz - yz*yz;
        const double det_y = xx*zz - xz*xz;
        const double det_z = xx*yy - xy*xy;
        const double det_max = MAX3(det_x,det_y,det_z);
        if (det_max == det_x) {
          normal_vp[ivp][0] = det_x;
          normal_vp[ivp][1] = xz*yz - xy*zz;
          normal_vp[ivp][2] = xy*yz - xz*yy;
        }
        else if (det_max == det_y) {
          normal_vp[ivp][0] = xz*yz - xy*zz;
          normal_vp[ivp][1] = det_y;
          normal_vp[ivp][2] = xy*xz - yz*xx;
        }
        else {
          assert(det_max == det_z);
          normal_vp[ivp][0] = xy*yz - xz*yy;
          normal_vp[ivp][1] = xy*xz - yz*xx;
          normal_vp[ivp][2] = det_z;
        }

        // sign each normal to be into the domain...

        const int ift = ftovp[ivp];
        const int ist = stoft[ift];
        const double normal2[3] = TRI_NORMAL_2(surface->xsp[surface->spost[ist][0]],surface->xsp[surface->spost[ist][1]],surface->xsp[surface->spost[ist][2]]);
        const double dp = DOT_PRODUCT(normal_vp[ivp],normal2);
        if (dp > 0.0)
          FOR_I3 normal_vp[ivp][i] *= -1.0;

        // make unit..

        const double mag = MAG(normal_vp[ivp]);
        if (mag > 0.0)
          FOR_I3 normal_vp[ivp][i] /= mag;
      }
      /*
	for (int ivp = ivp0; ivp < ivp1; ++ivp) {
        const int ift = ftovp[ivp]; 
        const int ist = stoft[ift];
        FOR_I3 normal_vp[ivp][i] = 0.0;
        FOR_I3 {
	const int isp = surface->spost[ist][i];
	const int ifp = fposp[isp]; assert((ifp >= 0)&&(ifp < nfp));
	const double n[3] = TRI_NORMAL_2(xvp[ivp],surface->xsp[surface->spost[ist][(i+1)%3]],surface->xsp[surface->spost[ist][(i+2)%3]]);
	const double n_mag = MAG(n);
	FOR_J3 normal_vp[ivp][j] -= n_mag*normal_fp[ifp][j]; 
        }
        const double mag = MAG(normal_vp[ivp]);
        assert(mag > 0.0);
        FOR_I3 normal_vp[ivp][i] /= mag;
	}
      */
      MPI_Barrier(mpi_comm_shared);

      // now each point gets an "nnvp", the local number of points...

      if (mpi_rank == 0) 
        cout << " > divergence and taper..." << endl;

      // now reduce nn in any tri based on divergence/convergence of the normals...

      //for (int ivp = ivp0; ivp < ivp1; ++ivp) 
      //  nnvp[ivp] = min(nnvp[ivp],log(dist_vp[ivp]/dn*(sf-1.0)+1.0)/log(sf)); 
      //MPI_Barrier(mpi_comm_shared);
      //CTI_Munmap(dist_vp,nvp);

      for (int ivp = ivp0; ivp < ivp1; ++ivp) {
        for (int vov = vpovp_i[ivp]; vov != vpovp_i[ivp+1]; ++vov) {
          const int ivp_nbr = vpovp_v[vov];
          const double dx0[3] = DIFF(xvp[ivp_nbr],xvp[ivp]);
          const double dx0_mag2 = DOT_PRODUCT(dx0,dx0);
          int j;
          for (j = 1; j < nn; ++j) {
            const double delta = dn*(1.0-pow(sf,j))/(1.0-sf);
            double x0[3]; FOR_K3 x0[k] = xvp[ivp][k]+delta*normal_vp[ivp][k]; 
            double x1[3]; FOR_K3 x1[k] = xvp[ivp_nbr][k]+delta*normal_vp[ivp_nbr][k]; 
            // check that the dp is aligned...
            const double dx[3] = DIFF(x1,x0);
            if (DOT_PRODUCT(dx,dx0) <= 0.0)
              break;
            const double rel_size = DOT_PRODUCT(dx,dx)/dx0_mag2;
            if ((rel_size <= 1.0/area_ratio)||(rel_size >= area_ratio)) // how much divergence to allow? (note ratio of squares)...
              break;
          }
          assert((j >= 1)&&(j <= nn));
          nnvp[ivp] = min(nnvp[ivp],double(j)+0.0001);
        }
      }
      MPI_Barrier(mpi_comm_shared);

      // finally, enforce taper...

      done = 0;
      while (done != 1) {

        int my_done = 1;
        for (int ivp = ivp0; ivp < ivp1; ++ivp) {
          for (int vov = vpovp_i[ivp]; vov != vpovp_i[ivp+1]; ++vov) {
            const int ivp_nbr = vpovp_v[vov];
            const double dx0[3] = DIFF(xvp[ivp_nbr],xvp[ivp]);
            const double dnn_limit = taper_factor*MAG(dx0)/dt;
            if (nnvp[ivp] > nnvp[ivp_nbr]+dnn_limit) {
              nnvp[ivp] = nnvp[ivp_nbr]+0.99999*dnn_limit;
              my_done = 0;
            }
          }
        }
        MPI_Allreduce(&my_done,&done,1,MPI_INT,MPI_MIN,mpi_comm_shared);
      }
      CTI_Munmap(vpovp_v,vpovp_i[nvp]);
      CTI_Munmap(vpovp_i,nvp+1);


      //--------------
      // ff_surface...
      //--------------

      if (mpi_rank == 0)
        cout << " > building ff_surface..." << endl;
      
      int* flag_ft = NULL;
      CTI_Mmap_rw(flag_ft,nft);
      for (int ift = ift0; ift < ift1; ++ift)
        flag_ft[ift] = 0; // hold the edges that are split bitwise
      MPI_Barrier(mpi_comm_shared);
      map<pair<int,int>,int> noMap; // isp0/isp1 -> iap
      vector<pair<int,int> > apVec; // iap -> ift/i
      ftVec.clear(); // ifts that are split
      if (mpi_rank_shared == 0) {
        for (int sol = 0; sol < groupDataVec[0].spoli_s; ++sol) {
          const int isp = groupDataVec[0].spoli_v[sol];
          const int ifp = fposp[isp]; assert((ifp >= 0)&&(ifp < nfp));
          for (int top = stofp_i[ifp]; top != stofp_i[ifp+1]; ++top) {
            const int ist = stofp_v[top];
            const int ift = ftost[ist]; 
            if ((ift != -1)&&(flag_ft[ift] == 0)) {
              int cnt = 0;
              bool b_split[3] = {false,false,false};
              FOR_I3 {
                int ist_nbr,i_nbr,bit;
                surface->unpackTeost(ist_nbr,i_nbr,bit,surface->teost[ist][i]);
                const int ift_nbr = ftost[ist_nbr];
                if (ift_nbr >= 0) 
                  b_split[i] = true;
                else 
                  ++cnt;
              }
              if (cnt >= 2) {
                FOR_I3 {
                  if (b_split[i]) {
                    int ist_nbr,i_nbr,bit;
                    surface->unpackTeost(ist_nbr,i_nbr,bit,surface->teost[ist][i]);
                    const int ift_nbr = ftost[ist_nbr];
                    // store split ft's for fast traversal below
                    if (flag_ft[ift] == 0)
                      ftVec.push_back(ift);
                    if (flag_ft[ift_nbr] == 0)
                      ftVec.push_back(ift_nbr);
                    // flag ft's with split edge
                    flag_ft[ift] |= 1<<i;
                    flag_ft[ift_nbr] |= 1<<i_nbr;
                    // look for reverse edge
                    map<pair<int,int>,int>::iterator it = noMap.find(pair<int,int>(surface->spost[ist][(i+1)%3],surface->spost[ist][i]));
                    if (it == noMap.end()) {
                      // add edge to map
                      noMap[pair<int,int>(surface->spost[ist][i],surface->spost[ist][(i+1)%3])] = apVec.size();
                      apVec.push_back(pair<int,int>(ift,i));
                    }
                  }
                }
              }
            }
          }
        }
        assert(apVec.size() == noMap.size());
      }
      MPI_Barrier(mpi_comm_shared);
      assert(part->ff_surface == NULL);
      part->ff_surface = new SurfaceShm();
      int nap = apVec.size();
      part->ff_surface->nsp = nap;
      MPI_Bcast(&part->ff_surface->nsp,1,MPI_INT,0,mpi_comm_shared);
      part->ff_surface->nst_max = part->ff_surface->nst = 2*part->ff_surface->nsp + nft; // split edge: add two tris for each point
      part->ff_surface->nsp += nfp;
      part->ff_surface->nsp_max = part->ff_surface->nsp;
      if (mpi_rank == 0) 
        cout << " > ff_surface nsp: " << part->ff_surface->nsp << ", nst: " << part->ff_surface->nst << endl; 
      assert(part->ff_surface->xsp == NULL);
      CTI_Mmap_rw(part->ff_surface->xsp,part->ff_surface->nsp_max); 
      assert(part->ff_surface->spost == NULL);
      CTI_Mmap_rw(part->ff_surface->spost,part->ff_surface->nst_max);
      assert(part->ff_surface->znost == NULL);
      CTI_Mmap_rw(part->ff_surface->znost,part->ff_surface->nst_max);
      assert(part->ff_surface_dxost == NULL);
      CTI_Mmap_rw(part->ff_surface_dxost,part->ff_surface->nst_max);
      if (part->b_name) {
        // if the user gave the NAME param, then use this to set the zone name...
        part->ff_surface->zoneVec.push_back(StZone(part->name+"_FF"));
      }
      else {
        // otherwise, use the profile name, or a default...
        part->ff_surface->zoneVec.push_back(StZone("SVP_STRAND_FF"));
      }
      if (mpi_rank == 0) 
        cout << " > divergence and taper..." << endl;

      double* nn_fp = NULL;
      CTI_Mmap_rw(nn_fp,nfp);
      /*
	for (int ifp = ifp0; ifp < ifp1; ++ifp) 
	nn_fp[ifp] = double(nn)+0.0001;
	MPI_Barrier(mpi_comm_shared);
	if (mpi_rank_shared == 0) {
	for (int sol = 0; sol < groupDataVec[0].spoli_s; ++sol) {
	const int isp = groupDataVec[0].spoli_v[sol];
	const int ifp = fposp[isp]; assert((ifp >= 0)&&(ifp < nfp));
	for (int top = stofp_i[ifp]; top != stofp_i[ifp+1]; ++top) {
	const int ist = stofp_v[top];
	const int ift = ftost[ist]; 
	if (ift == -1) {
	nn_fp[ifp] = 0.0001; // boundary
	break;
	}
	}
	}
	}
	MPI_Barrier(mpi_comm_shared);
      */
      for (int ifp = ifp0; ifp < ifp1; ++ifp) {
        nn_fp[ifp] = -1.0;
        const int isp = spofp[ifp];
        double dist2_min = HUGE_VAL;
        for (int top = stofp_i[ifp]; top != stofp_i[ifp+1]; ++top) {
          const int ist = stofp_v[top];
          const int ift = ftost[ist]; 
          if (ift == -1) {
            if (!b_axis) {
              nn_fp[ifp] = 0.0001; // boundary
              break;
            }
          }
          else {
            for (int vof = vpoft_i[ift]; vof != vpoft_i[ift+1]; ++vof) {
              const int ivp = vpoft_v[vof];
              const double dist2 = DIST2(xvp[ivp],surface->xsp[isp]);
              if (dist2 < dist2_min) {
                dist2_min = dist2;
                nn_fp[ifp] = nnvp[ivp];
              }
            }
          }
        }
      }
      MPI_Barrier(mpi_comm_shared);
      // fill...
      {
        int done = 0;
        while (done != 1) {

          int my_done = 1;
          for (int ifp = ifp0; ifp < ifp1; ++ifp) {
            if (nn_fp[ifp] < 0.0) {
              my_done = 0;
              const int isp = spofp[ifp];
              double sum_nn_fp = 0.0;
              int cnt = 0;
              for (int top = stofp_i[ifp]; top != stofp_i[ifp+1]; ++top) {
                const int ist = stofp_v[top];
                const int ift = ftost[ist]; 
                if (ift >= 0) {
                  FOR_I3 {
                    const int isp_nbr = surface->spost[ist][i];
                    if (isp != isp_nbr) {
                      const int ifp_nbr = fposp[isp_nbr];
                      if (nn_fp[ifp_nbr] >= 0.0001) {
                        sum_nn_fp += nn_fp[ifp_nbr];
                        ++cnt;
                      }
                    }
                  }
                }
              }
              if (cnt > 0) {
                nn_fp[ifp] = sum_nn_fp/double(cnt);
              }
            }
          }
          MPI_Allreduce(&my_done,&done,1,MPI_INT,MPI_MIN,mpi_comm_shared);
        }
      }
      /*
      // walk it out too to fill in flagged points without any vp nbrs...
      if (mpi_rank_shared == 0) {
      bool b_done = false;
      while (!b_done) {
      b_done = true;
      for (int ift = 0; ift < nft; ++ift) {
      const int ist = stoft[ift];
      int max_nnvp = -1;
      FOR_I3 {
      const int isp = surface->spost[ist][i];
      const int ifp = fposp[isp]; assert((ifp >= 0)&&(ifp < nfp));
      if (nn_fp[ifp] >= 0) {
      max_nnvp = max(max_nnvp,nn_fp[ifp]);
      }
      else {
      assert(nn_fp[ifp] == -1);
      b_done = false;
      }
      }
      // dont propogate 0's...
      if (max_nnvp > 0) {
      FOR_I3 {
      const int isp = surface->spost[ist][i];
      const int ifp = fposp[isp];
      if (nn_fp[ifp] == -1)
      nn_fp[ifp] = max_nnvp;
      }
      }
      }
      }
      }
      */
      //for (int ifp = ifp0; ifp < ifp1; ++ifp) 
      //  cout << nn_fp[ifp] << endl;
      // divergence limiting...
      for (int ifp = ifp0; ifp < ifp1; ++ifp) {
        const int isp = spofp[ifp];
        for (int top = stofp_i[ifp]; top != stofp_i[ifp+1]; ++top) {
          const int ist = stofp_v[top];
          const int ift = ftost[ist]; 
          if (ift >= 0) {
            FOR_I3 {
              const int isp_nbr = surface->spost[ist][i];
              if (isp != isp_nbr) {
                const int ifp_nbr = fposp[isp_nbr];
                const double dx0[3] = DIFF(surface->xsp[isp_nbr],surface->xsp[isp]);
                const double dx0_mag2 = DOT_PRODUCT(dx0,dx0);
                int j;
                for (j = 1; j < nn; ++j) {
                  const double delta = dn*(1.0-pow(sf,j))/(1.0-sf);
                  double x0[3]; FOR_K3 x0[k] = surface->xsp[isp][k]-delta*normal_fp[ifp][k]; 
                  double x1[3]; FOR_K3 x1[k] = surface->xsp[isp_nbr][k]-delta*normal_fp[ifp_nbr][k]; 
                  // check that the dp is aligned...
                  const double dx[3] = DIFF(x1,x0);
                  if (DOT_PRODUCT(dx,dx0) <= 0.0)
                    break;
                  const double rel_size = DOT_PRODUCT(dx,dx)/dx0_mag2;
                  if ((rel_size <= 1.0/area_ratio)||(rel_size >= area_ratio)) 
                    break;
                }
                assert((j >= 1)&&(j <= nn));
                nn_fp[ifp] = min(nn_fp[ifp],double(j)+0.0001);
              }
            }
          }
        }
      }
      MPI_Barrier(mpi_comm_shared);
      // taper...
      {
        int done = 0;
        while (done != 1) {

          int my_done = 1;
          for (int ifp = ifp0; ifp < ifp1; ++ifp) {
            const int isp = spofp[ifp];
            for (int top = stofp_i[ifp]; top != stofp_i[ifp+1]; ++top) {
              const int ist = stofp_v[top];
              const int ift = ftost[ist]; 
              if (ift >= 0) {
                FOR_I3 {
                  const int isp_nbr = surface->spost[ist][i];
                  if (isp != isp_nbr) {
                    const int ifp_nbr = fposp[isp_nbr];
                    const double dx0[3] = DIFF(surface->xsp[isp_nbr],surface->xsp[isp]);
                    const double dnn_limit = taper_factor*MAG(dx0)/dt;
                    if (nn_fp[ifp] > nn_fp[ifp_nbr]+dnn_limit) {
                      nn_fp[ifp] = nn_fp[ifp_nbr]+0.99999*dnn_limit;
                      my_done = 0;
                    }
                  }

                }
              }
            }
          }
          MPI_Allreduce(&my_done,&done,1,MPI_INT,MPI_MIN,mpi_comm_shared);
        }
      }
      // limit nnvp by nn_fp...
      for (int ivp = ivp0; ivp < ivp1; ++ivp) {
        const int ift = ftovp[ivp];
        const int ist = stoft[ift];
        double this_nn = 0.0;
        double n_sum = 0.0;
        FOR_I3 {
          const int isp = surface->spost[ist][i];
          const int ifp = fposp[isp]; assert((ifp >= 0)&&(ifp < nfp));
          const double n[3] = TRI_NORMAL_2(xvp[ivp],surface->xsp[surface->spost[ist][(i+1)%3]],surface->xsp[surface->spost[ist][(i+2)%3]]);
          const double n_mag = MAG(n);
          this_nn += n_mag*nn_fp[ifp]; 
          n_sum += n_mag;
        }
        if (n_sum > 0.0) {
          this_nn /= n_sum;
          nnvp[ivp] = max(1.0001,min(nnvp[ivp],this_nn));
        }
      }
      MPI_Barrier(mpi_comm_shared);
      CTI_Munmap(ftovp,nvp);
      for (int ifp = ifp0; ifp < ifp1; ++ifp) {
        assert(nn_fp[ifp] >= 0.0);
        const double delta = dn*(1.0-pow(sf,(int)nn_fp[ifp]))/(1.0-sf);
        const int isp = spofp[ifp];
        FOR_I3 part->ff_surface->xsp[ifp][i] = surface->xsp[isp][i]-delta*normal_fp[ifp][i];
      }
      MPI_Barrier(mpi_comm_shared);
      CTI_Munmap(spofp,nfp);
      if (mpi_rank_shared == 0) {
        for (int iap = 0; iap < nap; ++iap) {
          const int ift = apVec[iap].first;
          const int ist = stoft[ift]; 
          const int i = apVec[iap].second;
          int ist_nbr,i_nbr,bit;
          surface->unpackTeost(ist_nbr,i_nbr,bit,surface->teost[ist][i]);
          const int ift_nbr = ftost[ist_nbr]; assert(ift_nbr >= 0);
          const int isp0 = surface->spost[ist][i];
          const int ifp0 = fposp[isp0]; assert((ifp0 >= 0)&&(ifp0 < nfp));
          const int isp1 = surface->spost[ist][(i+1)%3];
          const int ifp1 = fposp[isp1]; assert((ifp1 >= 0)&&(ifp1 < nfp));
          FOR_J3 part->ff_surface->xsp[nfp+iap][j] = 0.5*(surface->xsp[isp0][j]+surface->xsp[isp1][j]);
          int nn_ap = 1;
          double dist2_min = HUGE_VAL;
          for (int vof = vpoft_i[ift]; vof != vpoft_i[ift+1]; ++vof) {
            const int ivp = vpoft_v[vof];
            const double dist2 = DIST2(xvp[ivp],part->ff_surface->xsp[nfp+iap]);
            if (dist2 < dist2_min) {
              dist2_min = dist2;
              nn_ap = (int)nnvp[ivp];
            }
          }
          for (int vof = vpoft_i[ift_nbr]; vof != vpoft_i[ift_nbr+1]; ++vof) {
            const int ivp = vpoft_v[vof];
            const double dist2 = DIST2(xvp[ivp],part->ff_surface->xsp[nfp+iap]);
            if (dist2 < dist2_min) {
              dist2_min = dist2;
              nn_ap = (int)nnvp[ivp];
            }
          }
          if (dist2_min == HUGE_VAL) {
            const int ifp_ = fposp[surface->spost[ist_nbr][(i_nbr+2)%3]];
            assert((ifp_ >= 0)&&(ifp_ < nfp));
            if (nn_fp[ifp_] >= 1.0)
              nn_ap = (int)nn_fp[ifp_];
          }
          assert((nn_ap >= 1)&&(nn_ap <= nn));
          const double delta = dn*(1.0-pow(sf,nn_ap))/(1.0-sf);
          double normal_ap[3]; FOR_J3 normal_ap[j] = 0.5*(normal_fp[ifp0][j]+normal_fp[ifp1][j]);
          const double mag = MAG(normal_ap); assert(mag > 0.0);
          FOR_J3 part->ff_surface->xsp[nfp+iap][j] -= delta*normal_ap[j]/mag;
        }
      }
      MPI_Barrier(mpi_comm_shared);
      CTI_Munmap(nn_fp,nfp);
      CTI_Munmap(normal_fp,nfp);
      CTI_Munmap(stofp_v,stofp_i[nfp]);
      CTI_Munmap(stofp_i,nfp+1);

      for (int ift = ift0; ift < ift1; ++ift) {
        const int ist = stoft[ift];
        if (flag_ft[ift] == 0) {
          part->ff_surface->znost[ift] = 0;
          FOR_I3 part->ff_surface->spost[ift][i] = fposp[surface->spost[ist][2-i]];
          part->ff_surface_dxost[ift] = dt;
          for (int vof = vpoft_i[ift]; vof != vpoft_i[ift+1]; ++vof) {
            const int ivp = vpoft_v[vof];
            part->ff_surface_dxost[ift] = min(part->ff_surface_dxost[ift],sqrt(r2vp[ivp]));
          }
        }
      }
      MPI_Barrier(mpi_comm_shared);
      if (mpi_rank_shared == 0) {
        for (int ii = 0, lim = ftVec.size(); ii < lim; ++ii) {
          const int ift = ftVec[ii];
          assert(flag_ft[ift] > 0);
          const int ist = stoft[ift];
          int cnt = 0;
          int i_split[3] = {-1,-1,-1};
          FOR_I3 if (flag_ft[ift]&(1<<i)) {
            if (i_split[0] == -1)
              i_split[0] = i;
            else if (i_split[1] == -1)
              i_split[1] = i;
            else 
              i_split[2] = i;
            ++cnt;
          }
          if (cnt == 1) {
            map<pair<int,int>,int>::iterator it = noMap.find(pair<int,int>(surface->spost[ist][i_split[0]],surface->spost[ist][(i_split[0]+1)%3]));
            int ift_new,iap;
            if (it == noMap.end()) {
              it = noMap.find(pair<int,int>(surface->spost[ist][(i_split[0]+1)%3],surface->spost[ist][i_split[0]]));
              assert(it != noMap.end());
              iap = it->second;
              const int ift_nbr = apVec[iap].first;
              const int ist_nbr = stoft[ift_nbr];
              const int i_nbr = apVec[iap].second;
              int ist_nbr_nbr,i_nbr_nbr,bit;
              surface->unpackTeost(ist_nbr_nbr,i_nbr_nbr,bit,surface->teost[ist_nbr][i_nbr]);
              assert(ist_nbr_nbr == ist);
              assert(i_nbr_nbr == i_split[0]);
              ift_new = nft+2*iap+1;
            }
            else {
              iap = it->second;
              assert(ift == apVec[iap].first);
              assert(i_split[0] == apVec[iap].second);
              ift_new = nft+2*iap;
            }
            assert(ift_new < part->ff_surface->nst);
            double dxost = dt;
            for (int vof = vpoft_i[ift]; vof != vpoft_i[ift+1]; ++vof) {
              const int ivp = vpoft_v[vof];
              dxost = min(dxost,sqrt(r2vp[ivp]));
            }
            part->ff_surface->znost[ift] = 0;
            part->ff_surface_dxost[ift] = dxost;
            part->ff_surface->spost[ift][0] = fposp[surface->spost[ist][(i_split[0]+2)%3]];
            part->ff_surface->spost[ift][1] = fposp[surface->spost[ist][(i_split[0]+1)%3]];
            part->ff_surface->spost[ift][2] = iap+nfp;
            part->ff_surface->znost[ift_new] = 0;
            part->ff_surface_dxost[ift_new] = dxost;
            part->ff_surface->spost[ift_new][0] = fposp[surface->spost[ist][i_split[0]]];
            part->ff_surface->spost[ift_new][1] = fposp[surface->spost[ist][(i_split[0]+2)%3]];
            part->ff_surface->spost[ift_new][2] = iap+nfp;
          }
          else if (cnt == 2) {
            int ift_new[2],iap[2];
            FOR_I2 {
              map<pair<int,int>,int>::iterator it = noMap.find(pair<int,int>(surface->spost[ist][i_split[i]],surface->spost[ist][(i_split[i]+1)%3]));
              if (it == noMap.end()) {
                it = noMap.find(pair<int,int>(surface->spost[ist][(i_split[i]+1)%3],surface->spost[ist][i_split[i]]));
                assert(it != noMap.end());
                iap[i] = it->second;
                const int ift_nbr = apVec[iap[i]].first;
                const int ist_nbr = stoft[ift_nbr];
                const int i_nbr = apVec[iap[i]].second;
                int ist_nbr_nbr,i_nbr_nbr,bit;
                surface->unpackTeost(ist_nbr_nbr,i_nbr_nbr,bit,surface->teost[ist_nbr][i_nbr]);
                assert(ist_nbr_nbr == ist);
                assert(i_nbr_nbr == i_split[i]);
                ift_new[i] = nft+2*iap[i]+1;
              }
              else {
                iap[i] = it->second;
                assert(ift == apVec[iap[i]].first);
                assert(i_split[i] == apVec[iap[i]].second);
                ift_new[i] = nft+2*iap[i];
              }
              assert(ift_new[i] < part->ff_surface->nst);
            }
            double dxost = dt;
            for (int vof = vpoft_i[ift]; vof != vpoft_i[ift+1]; ++vof) {
              const int ivp = vpoft_v[vof];
              dxost = min(dxost,sqrt(r2vp[ivp]));
            }
            part->ff_surface->znost[ift] = 0;
            part->ff_surface_dxost[ift] = dxost;
            FOR_I2 {
              part->ff_surface->znost[ift_new[i]] = 0;
              part->ff_surface_dxost[ift_new[i]] = dxost;
            }
            if (i_split[1] == (i_split[0]+1)) {
              part->ff_surface->spost[ift][0] = iap[0]+nfp;
              part->ff_surface->spost[ift][1] = iap[1]+nfp;
              part->ff_surface->spost[ift][2] = fposp[surface->spost[ist][(i_split[0]+1)%3]];
              part->ff_surface->spost[ift_new[0]][0] = iap[0]+nfp;
              part->ff_surface->spost[ift_new[0]][1] = fposp[surface->spost[ist][(i_split[0]+2)%3]];
              part->ff_surface->spost[ift_new[0]][2] = iap[1]+nfp;
              part->ff_surface->spost[ift_new[1]][0] = fposp[surface->spost[ist][i_split[0]]];
              part->ff_surface->spost[ift_new[1]][1] = fposp[surface->spost[ist][(i_split[0]+2)%3]];
              part->ff_surface->spost[ift_new[1]][2] = iap[0]+nfp;
            }
            else {
              assert(i_split[1] == (i_split[0]+2));
              part->ff_surface->spost[ift][0] = iap[0]+nfp;
              part->ff_surface->spost[ift][1] = fposp[surface->spost[ist][(i_split[0]+2)%3]];
              part->ff_surface->spost[ift][2] = fposp[surface->spost[ist][(i_split[0]+1)%3]];
              part->ff_surface->spost[ift_new[0]][0] = iap[1]+nfp;
              part->ff_surface->spost[ift_new[0]][1] = fposp[surface->spost[ist][(i_split[0]+2)%3]];
              part->ff_surface->spost[ift_new[0]][2] = iap[0]+nfp;
              part->ff_surface->spost[ift_new[1]][0] = fposp[surface->spost[ist][i_split[0]]];
              part->ff_surface->spost[ift_new[1]][1] = iap[1]+nfp;
              part->ff_surface->spost[ift_new[1]][2] = iap[0]+nfp;
            }
          }
          else {
            assert(cnt == 3);
            int ift_new[3],iap[3];
            FOR_I3 {
              assert(i_split[i] == i);
              map<pair<int,int>,int>::iterator it = noMap.find(pair<int,int>(surface->spost[ist][i],surface->spost[ist][(i+1)%3]));
              if (it == noMap.end()) {
                it = noMap.find(pair<int,int>(surface->spost[ist][(i+1)%3],surface->spost[ist][i]));
                assert(it != noMap.end());
                iap[i] = it->second;
                const int ift_nbr = apVec[iap[i]].first;
                const int ist_nbr = stoft[ift_nbr];
                const int i_nbr = apVec[iap[i]].second;
                int ist_nbr_nbr,i_nbr_nbr,bit;
                surface->unpackTeost(ist_nbr_nbr,i_nbr_nbr,bit,surface->teost[ist_nbr][i_nbr]);
                assert(ist_nbr_nbr == ist);
                assert(i_nbr_nbr == i);
                ift_new[i] = nft+2*iap[i]+1;
              }
              else {
                iap[i] = it->second;
                assert(ift == apVec[iap[i]].first);
                assert(i == apVec[iap[i]].second);
                ift_new[i] = nft+2*iap[i];
              }
              assert(ift_new[i] < part->ff_surface->nst);
            }
            double dxost = dt;
            for (int vof = vpoft_i[ift]; vof != vpoft_i[ift+1]; ++vof) {
              const int ivp = vpoft_v[vof];
              dxost = min(dxost,sqrt(r2vp[ivp]));
            }
            part->ff_surface->znost[ift] = 0;
            part->ff_surface_dxost[ift] = dxost;
            part->ff_surface->spost[ift][0] = iap[2]+nfp;
            part->ff_surface->spost[ift][1] = iap[1]+nfp;
            part->ff_surface->spost[ift][2] = iap[0]+nfp;
            FOR_I3 {
              part->ff_surface->znost[ift_new[i]] = 0;
              part->ff_surface_dxost[ift_new[i]] = dxost;
              part->ff_surface->spost[ift_new[i]][0] = iap[i]+nfp;
              part->ff_surface->spost[ift_new[i]][1] = iap[(i+1)%3]+nfp;
              part->ff_surface->spost[ift_new[i]][2] = fposp[surface->spost[ist][(i+1)%3]];
            }
          }
        }
      }
      apVec.clear();
      ftVec.clear();
      noMap.clear();
      CTI_Munmap(ftost,surface->nst);
      CTI_Munmap(flag_ft,nft);
      CTI_Munmap(fposp,surface->nsp);
      CTI_Munmap(vpoft_v,nvp);
      CTI_Munmap(vpoft_i,nft+1);
      CTI_Munmap(stoft,nft);

      if (mpi_rank == 0) {
        char filename[128];
        sprintf(filename,"ff_surface.%02d.sbin",(int)partVec.size());
        GeomUtils::writeSbin(filename,part->ff_surface->spost,part->ff_surface->nst,part->ff_surface->xsp);
      }

    }

    // copy strand base and ff_surface to other nodes...
    if (part->ff_surface == NULL)
      part->ff_surface = new SurfaceShm();
    int buf[3];
    if (mpi_rank == 0) {
      buf[0] = part->ff_surface->nst;
      buf[1] = part->ff_surface->nsp;
      buf[2] = nvp;
    }
    MPI_Bcast(buf,3,MPI_INT,0,mpi_comm);
    if (mpi_rank_internode_actual == 0) {
      assert(buf[0] == part->ff_surface->nst);
      assert(buf[1] == part->ff_surface->nsp);
      assert(buf[2] == nvp);
    }
    else {
      part->ff_surface->nst_max = part->ff_surface->nst = buf[0];
      part->ff_surface->nsp_max = part->ff_surface->nsp = buf[1];
      nvp = buf[2];
      assert(part->ff_surface->xsp == NULL);
      CTI_Mmap_rw(part->ff_surface->xsp,part->ff_surface->nsp_max); 
      assert(part->ff_surface->spost == NULL);
      CTI_Mmap_rw(part->ff_surface->spost,part->ff_surface->nst_max);
      assert(part->ff_surface->znost == NULL);
      CTI_Mmap_rw(part->ff_surface->znost,part->ff_surface->nst_max);
      assert(part->ff_surface_dxost == NULL);
      CTI_Mmap_rw(part->ff_surface_dxost,part->ff_surface->nst_max);
      if (part->b_name) {
        // if the user gave the NAME param, then use this to set the zone name...
        part->ff_surface->zoneVec.push_back(StZone(part->name+"_FF"));
      }
      else {
        // otherwise, use the profile name, or a default...
        part->ff_surface->zoneVec.push_back(StZone("SVP_STRAND_FF"));
      }
      assert(xvp == NULL); CTI_Mmap(xvp,nvp);
      assert(normal_vp == NULL); CTI_Mmap(normal_vp,nvp);
      assert(r2vp == NULL); CTI_Mmap(r2vp,nvp);
      assert(nnvp == NULL); CTI_Mmap(nnvp,nvp);
    }
    if (mpi_rank_shared == 0) {
      MPI_Bcast((double*)part->ff_surface->xsp,3*part->ff_surface->nsp,MPI_DOUBLE,0,mpi_comm_internode);
      MPI_Bcast((int*)part->ff_surface->spost,3*part->ff_surface->nst,MPI_INT,0,mpi_comm_internode);
      MPI_Bcast(part->ff_surface->znost,part->ff_surface->nst,MPI_INT,0,mpi_comm_internode);
      MPI_Bcast(part->ff_surface_dxost,part->ff_surface->nst,MPI_DOUBLE,0,mpi_comm_internode);
      MPI_Bcast((double*)xvp,3*nvp,MPI_DOUBLE,0,mpi_comm_internode);
      MPI_Bcast((double*)normal_vp,3*nvp,MPI_DOUBLE,0,mpi_comm_internode);
      MPI_Bcast(r2vp,nvp,MPI_DOUBLE,0,mpi_comm_internode);
      MPI_Bcast(nnvp,nvp,MPI_DOUBLE,0,mpi_comm_internode);
    }
    MPI_Barrier(mpi_comm_shared);

    /*
      if (mpi_rank == 0) {
      vector<pair<SimplePoint,int> > xpVec(nvp);
      for (int ip = 0; ip < nvp; ++ip)
      xpVec[ip] = pair<SimplePoint,int>(SimplePoint(xvp[ip]),ip);
      sort(xpVec.begin(),xpVec.end());
      int nblank = 0;
      for (int ip = 1; ip < nvp; ++ip) {
      if (DIST2(xpVec[ip].first.x,xpVec[ip-1].first.x) < dn*dn) {
      ++nblank;
      const int ivp = xpVec[ip].second;
      nnvp[ivp] = 0.0;
      //cout << xpVec[ip].second << " " << xpVec[ip-1].second << " " << COUT_VEC(xpVec[ip].first.x) << " " << COUT_VEC(xpVec[ip-1].first.x) << " " << DIST(xpVec[ip].first.x,xpVec[ip-1].first.x) << endl;
      }
      }
      xpVec.clear();
      cout << " > number of blanked coincident surface voronoi points: " << nblank << "/" << nvp << endl;
      }
    */

    //------------
    // points...
    //------------

    if (mpi_rank == 0)
      cout << " > building points..." << endl;

    // now make split based on full comm...

    int nvp_split_2 = nvp/mpi_size;
    if (nvp%mpi_size != 0) nvp_split_2 += 1;
    const int ivp0_2 = min(nvp,nvp_split_2*mpi_rank);
    const int ivp1_2 = min(nvp,nvp_split_2*(mpi_rank+1));

    assert(part->pts == NULL);
    part->pts = new Points();
    part->pts->np_global = 0;
    part->pts->np = 0;
    for (int ivp = ivp0_2; ivp < ivp1_2; ++ivp) 
      part->pts->np += int(nnvp[ivp]);
    assert(part->pts->xp == NULL); part->pts->xp = new double[part->pts->np][3];
    assert(part->pts->delta == NULL); part->pts->delta = new double[part->pts->np];
    int8 my_np = part->pts->np;
    MPI_Allreduce(&my_np,&part->pts->np_global,1,MPI_INT8,MPI_SUM,mpi_comm);
    if (mpi_rank == 0) 
      cout << " > np_global: " << part->pts->np_global << endl;

    // populate points... 

    int ip = 0;
    for (int ivp = ivp0_2; ivp < ivp1_2; ++ivp) {
      for (int jj = 0; jj < int(nnvp[ivp]); ++jj) {
        const double delta = dn*(1.0-pow(sf,double(jj+0.5)))/(1.0-sf);
        FOR_I3 part->pts->xp[ip][i] = xvp[ivp][i]+delta*normal_vp[ivp][i]; 
        part->pts->delta[ip] = sqrt(r2vp[ivp]); 
        ++ip;
      }
    }
    assert(ip == part->pts->np);
    CTI_Munmap(r2vp,nvp);
    CTI_Munmap(nnvp,nvp);
    CTI_Munmap(normal_vp,nvp); 
    CTI_Munmap(xvp,nvp); 

    /*
      double (*xp)[3] = NULL;
      if (mpi_rank == 0) 
      xp = new double[part->pts->np_global][3];
      int * recv_count = NULL;
      if (mpi_rank == 0) recv_count = new int[mpi_size];
      int count = 3*part->pts->np;
      MPI_Gather(&count,1,MPI_INT,recv_count,1,MPI_INT,0,mpi_comm);
      int * recv_disp = NULL;
      if (mpi_rank == 0) {
      recv_disp = new int[mpi_size];
      recv_disp[0] = 0;
      for (int rank = 1; rank < mpi_size; ++rank)
      recv_disp[rank] = recv_disp[rank-1] + recv_count[rank-1];
      }
      MPI_Gatherv((double*)part->pts->xp,3*part->pts->np,MPI_DOUBLE,xp,recv_count,recv_disp,MPI_DOUBLE,0,mpi_comm);
      if (mpi_rank == 0) {
      delete[] recv_count;
      delete[] recv_disp;
      vector<SimplePoint> xpVec(part->pts->np_global);
      for (int ip = 0; ip < part->pts->np_global; ++ip)
      xpVec[ip] = SimplePoint(xp[ip]);
      delete[] xp;
      sort(xpVec.begin(),xpVec.end());
      for (int ip = 1; ip < part->pts->np_global; ++ip) {
      //assert(xpVec[ip] != xpVec[ip-1]);
      if (DIST(xpVec[ip].x,xpVec[ip-1].x) < dn) {
      cout << COUT_VEC(xpVec[ip].x) << " " << COUT_VEC(xpVec[ip-1].x) << " " << DIST(xpVec[ip].x,xpVec[ip-1].x) << endl;
      }
      }
      xpVec.clear();
      }
    */

    //char filename[128];
    //sprintf(filename,"points.%06d.dat",mpi_rank);
    //GeomUtils::writePtsTecplot(string(filename),part->pts->xp,part->pts->np);

    if (mpi_rank_shared == 0) {
      for (int ipart = 0, limit = partVec.size(); ipart < limit; ++ipart) {
        if (partVec[ipart]->surface) {
          SurfaceShm * surface = partVec[ipart]->surface;
          for (int ist = 0; ist < surface->nst; ++ist) {
            int igr;
            if (partVec[ipart]->getGroupForSt(igr,ist)) {
              partVec[ipart]->clearGroupForSt(ist);
              partVec[ipart]->setPartForSt(part->ipart_of_part,ist);
            }
          }
        }
      }
    }
    MPI_Barrier(mpi_comm_shared);
    for (int igr = 0; igr < groupDataVec.size(); ++igr) 
      groupDataVec[igr].clear();
    groupDataVec.clear();

    MPI_Barrier(mpi_comm);
    return true;
  }

  bool makeExtrudeSplineWithStrand(Part * part,Param * param,int& iarg) {

    int profile_type = -1;
    string profile_name;
    vector<double> profileVec;
    bool b_name = false;
    string name;
    bool b_x0 = false;
    double x0[3];
    bool b_x1 = false;
    double x1[3];
    bool b_dn = false;
    double dn; // normal
    bool b_dt = false;
    double dt; // tangential
    bool b_ds = false;
    double ds; // span (== tangential if not specified separately)
    bool b_nn;
    int nn; // normal count: how many spaces to go from dn to dt
    int periodic_type = -1; // 0:x, 1:y, 2:z

    while (iarg < param->size()) {
      const string token = param->getString(iarg++);
      if (token == "NAME") {
        b_name = true;
        name = param->getString(iarg++);
      }
      else if (token == "PROFILE_XY") {
        profile_type = 0; // 0: XY, 1: XZ, 2: YZ
        profile_name = param->getString(iarg++);
        GeomUtils::readXY(profileVec,profile_name);
      }
      else if (token == "SCALE") {
        const double factor = param->getDouble(iarg++);
        // scale everything that has been read...
        if (profile_type != -1) {
          for (int ii = 0; ii < profileVec.size(); ++ii)
            profileVec[ii] *= factor;
        }
        assert(!b_x0);
        assert(!b_x1);
        assert(!b_dn);
        assert(!b_dt);
        assert(!b_ds);
      }
      else if (token == "X0") {
        b_x0 = true;
        FOR_I3 x0[i] = param->getDouble(iarg++);
      }
      else if (token == "X1") {
        b_x1 = true;
        FOR_I3 x1[i] = param->getDouble(iarg++);
      }
      else if (token == "DN") {
        b_dn = true;
        dn = param->getDouble(iarg++);
      }
      else if (token == "NN") {
        b_nn = true;
        nn = param->getInt(iarg++);
        assert(nn > 2);
      }
      else if (token == "DT") {
        b_dt = true;
        dt = param->getDouble(iarg++);
      }
      else if (token == "DS") {
        b_ds = true;
        ds = param->getDouble(iarg++);
      }
      else if (token == "PERIODIC_Z") {
        periodic_type = 2;
      }
      else if (token == "NLAYERS") {
	part->ff_nlayersString = param->getString(iarg++);
      }
      else {
	if (mpi_rank == 0) cout << "unrecognized EXTRUDE_SPLINE_WITH_STRAND token: " << token << endl;
      }
    }

    int ierr = 0;
    if (profile_type == -1) {
      if (mpi_rank == 0) cout << "EXTRUDE_SPLINE_WITH_STRAND: missing PROFILE_XY <filename>" << endl;
      ierr = -1;
    }
    if (!b_dn) {
      if (mpi_rank == 0) cout << "EXTRUDE_SPLINE_WITH_STRAND: missing DN <double>" << endl;
      ierr = -1;
    }
    if (!b_nn) {
      if (mpi_rank == 0) cout << "EXTRUDE_SPLINE_WITH_STRAND: missing NN <int>" << endl;
      ierr = -1;
    }
    if (!b_x0) {
      if (mpi_rank == 0) cout << "EXTRUDE_SPLINE_WITH_STRAND: missing X0 <x> <y> <z>" << endl;
      ierr = -1;
    }
    if (!b_x1) {
      if (mpi_rank == 0) cout << "EXTRUDE_SPLINE_WITH_STRAND: missing X1 <x> <y> <z>" << endl;
      ierr = -1;
    }
    if (!b_dt) {
      if (mpi_rank == 0) cout << "EXTRUDE_SPLINE_WITH_STRAND: missing DT <double>" << endl;
      ierr = -1;
    }
    if (ierr != 0) {
      return false;
    }

    // default for ds is dt if not specified separately...
    if (!b_ds) {
      ds = dt;
    }

    // set up the spline...
    int np;
    double (*xp)[3] = NULL;
    if (profile_type == 0) {
      assert(profileVec.size()%2 == 0);
      np = profileVec.size()/2;
      xp = new double[np][3];
      for (int ip = 0; ip < np; ++ip) {
        xp[ip][0] = profileVec[ip*2  ] + x0[0];
        xp[ip][1] = profileVec[ip*2+1] + x0[1];
        xp[ip][2] = x0[2];
      }
      profileVec.clear();
    }
    else {
      assert(0);
    }
    SplineStuff::CubicSpline cspline;
    cspline.init(xp,np,true); // true == periodic spline
    delete[] xp;

    // spline length can now be computed...
    const double L = cspline.getLength();
    const double sf = pow(dt/dn,1.0/(double(nn)-1.0)); assert(sf > 1.0);

    // discretize the spline at the smallest resolution: dn and
    // use the standoff distance to estimate the local spacing allowed...
    assert(dn <= dt);
    np = max(int(L/dn),1); // 1? or 2?
    xp = new double[np][3];
    for (int ip = 0; ip < np; ++ip) {
      const double t = double(ip)/double(np);
      cspline.getX(xp[ip],t);
    }

    /*
      {
      FILE * fp = fopen("spline0.dat","w");
      for (int ip = 0; ip < np; ++ip) {
      //for (int in = 0; in <= nn_local[it]; ++in) {
      fprintf(fp,"%18.15le %18.15le %18.15le\n",xp[ip][0],xp[ip][1],xp[ip][2]);
      }
      fclose(fp);
      cout << "take a look at spline0.dat" << endl;
      getchar();
      }
    */

    const double fraction = 0.1; // the fraction of dn that we can tolerate in standoff from the curve
    double *dxp = new double[np];
    for (int ip = 0; ip < np; ++ip) {
      const int ip_prev = (ip-1+np)%np;
      const int ip_next = (ip+1)%np;
      // get the stand-off distance...
      double dx[3];
      FOR_I3 dx[i] = xp[ip][i] - 0.5*(xp[ip_prev][i]+xp[ip_next][i]);
      const double standoff = MAG(dx);
      // standoff should vary quadratically with spacing, so...
      if (standoff > fraction*dn) {
        // too much curvature, so just use dn...
        dxp[ip] = dn;
      }
      else {
        const double root = fraction*dn*standoff*(L*L/(double(np)*double(np)) + standoff*standoff - fraction*dn*standoff);
        assert(root > 0.0);
        dxp[ip] = min(dt,sqrt(root)/standoff);
      }
    }
    delete[] xp; xp = NULL;

    // now reduce dxp everywhere based on a stretching ratio...
    bool done = false;
    while (!done) {
      done = true;
      for (int ip = 0; ip < np; ++ip) {
        const int ip_prev = (ip-1+np)%np;
        const int ip_next = (ip+1)%np;
        const double dxp_nbr_min = min(dxp[ip_prev],dxp[ip_next]);
        if (dxp[ip] > dxp_nbr_min*pow(sf,0.5*dn/dxp_nbr_min)) {
          dxp[ip] = dxp_nbr_min*pow(sf,0.5*dn/dxp_nbr_min);
          done = false;
        }
      }
    }

    int nt;
    double factor = 1.0;
    done = false;
    while (!done) {
      double dist = 0.0;
      nt = 1; // first point at start of spline...
      for (int ip = 0; ip < np; ++ip) {
        const int ip_next = (ip+1)%np;
        //in the interval ip to ip_next, we are looking for dist == dx...
        // dist = dist0 + alpha*L/np
        // and assume dx varies linearly in the interval...
        // dx = alpha*dx1 + (1-alpha)*dx0...
        const double alpha = (dist - dxp[ip])/(dxp[ip_next]-dxp[ip]-L/double(np)*factor);
        if (alpha <= 1.0) {
          assert(alpha >= 0.0);
          ++nt;
          // reset the distance to the remaining distance in the interval...
          dist = (1.0-alpha)*L/double(np)*factor;
        }
        else {
          // add the whole interval distance to dist...
          dist += L/double(np)*factor;
        }
      }
      if (dist < 0.5*dxp[0]) {
        factor *= (L-dist)/L;
      }
      else {
        factor *= (L-dist)/(L-dxp[0]);
      }
      if (dist < 1.0E-10*dxp[0]) {
        // this means we have added one too many points (i.e. first point was repeated)..
        --nt;
        done = true;
      }
      else if (dxp[0]-dist < 1.0E-10*dxp[0]) {
        // this means we have added with the last dist almost exactly equal to the required dxp[0]...
        done = true;
      }
    }

    if (mpi_rank == 0) cout << " > tangential spacing nt: " << nt << endl;

    // at this point, we build t_local...

    double * t_local = new double[nt];
    double dist = 0.0;
    t_local[0] = 0.0;
    int it = 1; // first point at start of spline...
    for (int ip = 0; ip < np; ++ip) {
      const int ip_next = (ip+1)%np;
      //in the interval ip to ip_next, we are looking for dist == dx...
      // dist = dist0 + alpha*L/np
      // and assume dx varies linearly in the interval...
      // dx = alpha*dx1 + (1-alpha)*dx0...
      const double alpha = (dist - dxp[ip])/(dxp[ip_next]-dxp[ip]-L/double(np)*factor);
      if (alpha <= 1.0) {
        assert(alpha >= 0.0);
        t_local[it++] = (double(ip)+alpha)/double(np);
        if (it == nt)
          break;
        dist = (1.0-alpha)*L/double(np)*factor;
      }
      else {
        // add the whole interval distance to dist...
        dist += L/double(np)*factor;
      }
    }
    delete[] dxp; dxp = NULL;

    // smooth t_local if requested by the user...

    const int nsmooth1 = 5;
    double * t_copy = new double[nt];
    for (int ismooth = 0; ismooth < nsmooth1; ++ismooth) {
      for (int it = 0; it < nt; ++it)
        t_copy[it] = t_local[it];
      // we want to keep the first point at zero, so its
      // adjustment is applied to all other points...
      const double dt0 = 0.5*(t_copy[nt-1]-1.0+t_copy[1]);
      for (int it = 1; it < nt; ++it) {
        if (it == nt-1) {
          t_local[nt-1] = 0.5*(t_copy[0]+1.0+t_copy[nt-2]) - dt0;
        }
        else {
          t_local[it] = 0.5*(t_copy[it-1]+t_copy[it+1]) - dt0;
        }
      }
    }
    delete[] t_copy;

    const double dx01[3] = DIFF(x1,x0); // spanwise vector

    // actually set the points...
    double (*xt)[3] = new double[nt*(nn+1)][3];
    //double (*dxt)[3] = new double[nt][3];
    for (int it = 0; it < nt; ++it) {
      cspline.getX(xt[it],t_local[it]);
    }

    for (int it = 0; it < nt; ++it) {
      const int it_prev = (it-1+nt)%nt;
      const int it_next = (it+1)%nt;
      const double dxt[3] = DIFF(xt[it_next],xt[it_prev]);
      const double cp[3] = CROSS_PRODUCT(dxt,dx01);
      const double cp_mag = MAG(cp); assert(cp_mag > 0.0);
      // and the outer points...
      const double delta_bl = GeomUtils::stretchingFunction(nn,nn,dn,dt);
      FOR_I3 xt[it+nt*nn][i] = xt[it][i] - delta_bl*cp[i]/cp_mag;
    }

    // smooth the outer points...
    const int nsmooth2 = 100;
    for (int ismooth = 0; ismooth < nsmooth2; ++ismooth) {
      for (int it = 0; it < nt; ++it)
        FOR_I3 xt[it+nt*(nn-1)][i] = xt[it+nt*(nn)][i];
      for (int it = 0; it < nt; ++it) {
        const int it_prev = (it-1+nt)%nt;
        const int it_next = (it+1)%nt;
        const double dx1[3] = DIFF(xt[it_next+nt*(nn-1)],xt[it_prev+nt*(nn-1)]);
        const double dx1_mag2 = DOT_PRODUCT(dx1,dx1);
        double dx2[3]; FOR_I3 dx2[i] = 0.5*(xt[it_next+nt*(nn-1)][i]+xt[it_prev+nt*(nn-1)][i])-xt[it+nt*(nn-1)][i];
        const double dp = DOT_PRODUCT(dx2,dx1);
        FOR_I3 xt[it+nt*(nn)][i] += dp*dx1[i]/dx1_mag2;
      }
    }

    // and finally build the strand points between the first and last points...

    for (int it = 0; it < nt; ++it) {
      const int it_prev = (it-1+nt)%nt;
      const int it_next = (it+1)%nt;
      const double dxt[3] = DIFF(xt[it_next],xt[it_prev]);
      const double cp[3] = CROSS_PRODUCT(dxt,dx01);
      const double cp_mag = MAG(cp); assert(cp_mag > 0.0);
      const double delta_bl = GeomUtils::stretchingFunction(nn,nn,dn,dt);
      double dx[3]; FOR_I3 dx[i] = xt[it+nt*(nn)][i] - (xt[it][i] - delta_bl*cp[i]/cp_mag);
      for (int in = 1; in < nn; ++in) {
        const double delta = GeomUtils::stretchingFunction(in,nn,dn,dt);
        const double ratio = delta/delta_bl; assert((ratio > 0.0)&&(ratio < 1.0));
        FOR_I3 xt[it+nt*in][i] = xt[it][i] - ratio*delta_bl*cp[i]/cp_mag + ratio*ratio*dx[i];
      }
    }

    // for now, just use the full height everywhere...
    int (*nn_local) = new int[nt];
    for (int it = 0; it < nt; ++it)
      nn_local[it] = nn;

    /*
      const double ratio = 1.25;
      for (int it = 0; it < nt; ++it) {
      const int it_prev = (it-1+nt)%nt;
      const int it_next = (it+1)%nt;
      int in;
      for (in = nn; in > 0; --in) {
      const double dt2 = DIST(xt[it_prev+nt*in],xt[it_next+nt*in]);
      if ((dt2 < 2.0*dt*ratio)&&(dt2 > 2.0*dt/ratio))
      break;
      }
      nn_local[it] = in;
      }

      // then ensure there is no single change greater than 1 level...
      done = false;
      while (!done) {
      done = true;
      for (int it = 0; it < nt; ++it) {
      const int it_prev = (it-1+nt)%nt;
      const int it_next = (it+1)%nt;
      const int nn_local_nbr_min = min(nn_local[it_prev],nn_local[it_next]);
      if (nn_local[it] > nn_local_nbr_min+1) {
      nn_local[it] = nn_local_nbr_min+1;
      done = false;
      }
      }
      }
    */

    /*
    // take a look...
    FILE * fp = fopen("spline.dat","w");
    for (int it = 0; it < nt; ++it) {
    for (int in = 0; in <= nn; ++in) {
    //for (int in = 0; in <= nn_local[it]; ++in) {
    fprintf(fp,"%18.15le %18.15le %18.15le\n",xt[it+nt*in][0],xt[it+nt*in][1],xt[it+nt*in][2]);
    }
    }
    fclose(fp);
    cout << "TAKE A LOOK!" << endl;
    getchar();
    */

    // =======================================
    // the surface...
    // =======================================

    assert(part->surface == NULL);
    part->surface = new SurfaceShm();
    part->surface->init(nt*3,nt*4); // nsp,nst
    if (b_name) {
      // if the user gave the NAME param, then use this to set the zone name...
      part->surface->zoneVec.push_back(StZone(name));
    }
    else {
      // otherwise, use the profile name, or a default...
      part->surface->zoneVec.push_back(StZone("EXTRUDE_SPLINE"));
    }
    assert(part->solid_type == UNKNOWN_SOLID);
    part->solid_type = SURFACE_SOLID;

    // recall surface is shm...
    if (mpi_rank_shared == 0) {
      for (int it = 0; it < nt; ++it) {
        FOR_I3 part->surface->xsp[it][i] = xt[it][i];
        FOR_I3 part->surface->xsp[it+nt][i] = xt[it][i] + 0.5*dx01[i]; // split the span into 2
        FOR_I3 part->surface->xsp[it+nt*2][i] = xt[it][i] + dx01[i];
        const int it_next = (it+1)%nt;
        // first pair of tris...
        part->surface->spost[it*4  ][0] = it;
        part->surface->spost[it*4  ][2] = it+nt;
        part->surface->spost[it*4  ][1] = it_next;
        part->surface->znost[it*4  ] = 0;
        part->surface->spost[it*4+1][0] = it_next;
        part->surface->spost[it*4+1][2] = it+nt;
        part->surface->spost[it*4+1][1] = it_next+nt;
        part->surface->znost[it*4+1] = 0;
        // second pair of tris...
        part->surface->spost[it*4+2][0] = it+nt;
        part->surface->spost[it*4+2][2] = it+nt+nt;
        part->surface->spost[it*4+2][1] = it_next+nt;
        part->surface->znost[it*4+2] = 0;
        part->surface->spost[it*4+3][0] = it_next+nt;
        part->surface->spost[it*4+3][2] = it+nt+nt;
        part->surface->spost[it*4+3][1] = it_next+nt+nt;
        part->surface->znost[it*4+3] = 0;
      }
      //GeomUtils::writeTecplot("surface.dat",part->surface->spost,part->surface->nst,part->surface->xsp);
    }
    MPI_Barrier(mpi_comm_shared);

    // =======================================
    // the ff_surface...
    // =======================================

    assert(part->ff_surface == NULL);
    part->ff_surface = new SurfaceShm();
    part->ff_surface->init(nt*4,nt*6); // nsp,nst
    if (b_name) {
      // if the user gave the NAME param, then use this to set the zone name...
      part->ff_surface->zoneVec.push_back(StZone(name+"_FF"));
    }
    else {
      // otherwise, use the profile name, or a default...
      part->ff_surface->zoneVec.push_back(StZone("EXTRUDE_SPLINE_FF"));
    }
    assert(part->ff_surface_dxost == NULL);
    CTI_Mmap(part->ff_surface_dxost,nt*6);
    assert(part->ff_type == UNKNOWN_FF);
    part->ff_type = FF_SURFACE_FF;

    // recall surface is shm...
    if (mpi_rank_shared == 0) {
      for (int it = 0; it < nt; ++it) {
        FOR_I3 part->ff_surface->xsp[it][i]      = xt[it][i];
        FOR_I3 part->ff_surface->xsp[it+nt][i]   = xt[it+nt*nn_local[it]][i];
        FOR_I3 part->ff_surface->xsp[it+2*nt][i] = xt[it+nt*nn_local[it]][i] + dx01[i];
        FOR_I3 part->ff_surface->xsp[it+3*nt][i] = xt[it][i] + dx01[i];
        const int it_next = (it+1)%nt;
        // first pair of tris...
        part->ff_surface->spost[it*6][0] = it;
        part->ff_surface->spost[it*6][1] = it+nt;
        part->ff_surface->spost[it*6][2] = it_next;
        part->ff_surface->znost[it*6] = 0;
        part->ff_surface_dxost[it*6] = dt;
        part->ff_surface->spost[it*6+1][0] = it_next;
        part->ff_surface->spost[it*6+1][1] = it+nt;
        part->ff_surface->spost[it*6+1][2] = it_next+nt;
        part->ff_surface->znost[it*6+1] = 0;
        part->ff_surface_dxost[it*6+1] = dt;
        // second pair of tris...
        part->ff_surface->spost[it*6+2][0] = it+nt;
        part->ff_surface->spost[it*6+2][1] = it+nt+nt;
        part->ff_surface->spost[it*6+2][2] = it_next+nt;
        part->ff_surface->znost[it*6+2] = 0;
        part->ff_surface_dxost[it*6+2] = dt;
        part->ff_surface->spost[it*6+3][0] = it_next+nt;
        part->ff_surface->spost[it*6+3][1] = it+nt+nt;
        part->ff_surface->spost[it*6+3][2] = it_next+nt+nt;
        part->ff_surface->znost[it*6+3] = 0;
        part->ff_surface_dxost[it*6+3] = dt;
        // third pair of tris...
        part->ff_surface->spost[it*6+4][0] = it+nt+nt;
        part->ff_surface->spost[it*6+4][1] = it+nt+nt+nt;
        part->ff_surface->spost[it*6+4][2] = it_next+nt+nt;
        part->ff_surface->znost[it*6+4] = 0;
        part->ff_surface_dxost[it*6+4] = dt;
        part->ff_surface->spost[it*6+5][0] = it_next+nt+nt;
        part->ff_surface->spost[it*6+5][1] = it+nt+nt+nt;
        part->ff_surface->spost[it*6+5][2] = it_next+nt+nt+nt;
        part->ff_surface->znost[it*6+5] = 0;
        part->ff_surface_dxost[it*6+5] = dt;
      }
      //GeomUtils::writeTecplot("ff_surface.dat",part->ff_surface->spost,part->ff_surface->nst,part->ff_surface->xsp);
    }
    MPI_Barrier(mpi_comm_shared);

    // =======================================
    // and finally the points: this part is distributed in parallel...
    // =======================================

    const double dx01_mag = MAG(dx01);
    const int ns = max(1,(int)floor(dx01_mag/ds+0.5));
    assert(part->pts == NULL);
    part->pts = new Points();
    int nt_local = nt/mpi_size;
    if ((nt_local%mpi_size != 0)||(nt_local == 0)) nt_local += 1;
    int it_start = min(nt,nt_local*mpi_rank);
    int it_end = min(nt,nt_local*(mpi_rank+1));
    part->pts->np = 0;
    part->pts->np_global = 0;
    for (int it = 0; it < nt; ++it) {
      // only count our local np...
      if ((it >= it_start)&&(it < it_end))
        part->pts->np += nn_local[it]*ns;
      // all count np_global...
      part->pts->np_global += nn_local[it]*ns;
    }
    if (mpi_rank == 0) cout << " > np_global: " << part->pts->np_global << endl;

    part->pts->xp = new double[part->pts->np][3];
    part->pts->delta = new double[part->pts->np];

    int ip = 0;
    for (int it = it_start; it < it_end; ++it) {
      for (int in = 0; in < nn_local[it]; ++in) {
        double x[3]; FOR_I3 x[i] = 0.5*(xt[it+nt*in][i] + xt[it+nt*(in+1)][i]);
        for (int is = 0; is < ns; ++is) {
          FOR_I3 part->pts->xp[ip][i] = x[i] + (double(is)+0.5)/double(ns)*dx01[i];
          part->pts->delta[ip] = 2.0*dt; // TODO: could be more precise here
          ++ip;
        }
      }
    }
    assert(ip == part->pts->np);

    /*
      GeomUtils::writePtsTecplot("points.dat",part->pts->xp,part->pts->np);

      // take a look...
      FILE * fp = fopen("spline.dat","w");
      for (int it = 0; it < nt; ++it) {
      for (int in = 0; in <= nn_local[it]; ++in) {
      fprintf(fp,"%18.15le %18.15le %18.15le\n",xt[it+nt*in][0],xt[it+nt*in][1],xt[it+nt*in][2]);
      }
      }
      fclose(fp);
    */

    delete[] nn_local;
    delete[] xt;

    // set any periodicity...
    if (periodic_type == 2) {
      // PERIODIC_Z...
      // confirm that the X0->X1 vector is a z-translation...
      //cout << "PERIODIC_Z: dx01: " << COUT_VEC(dx01) << endl;
      if ((dx01[0] != 0.0)||(dx01[1] != 0.0)||(dx01[2] == 0.0)) {
        if (mpi_rank == 0) cout << "WARNING: dx form X0 to X1 not consistent with PERIODIC_Z: " << COUT_VEC(dx01) << endl;
        ierr = -1;
      }
      else {
        // check if there is matching periodic data...
        int ii_match = -1;
        int ii_sign = 0;
        for (int ii = 0; ii < PeriodicData::periodicTransformVec.size(); ++ii) {
          const int kind = PeriodicData::periodicTransformVec[ii].getKind();
          if (kind == PERIODIC_TRANSFORM_CART) {
            double data[3];
            PeriodicData::periodicTransformVec[ii].getData(data);
            if ((data[0] == 0.0)&&(data[1] == 0.0)&&(data[2] != 0.0)) {
              if (data[2] == dx01[2]) {
                ii_match = ii;
                ii_sign = 1;
              }
              else if (data[2] == -dx01[2]) {
                ii_match = ii;
                ii_sign = -1;
              }
              else {
                if (mpi_rank == 0) cout << "WARNING: existing z-periodicity does not match: " << COUT_VEC(data) << endl;
                ierr = -1;
              }
              break;
            }
          }
        }
        if (ii_match == -1) {
          // no matching or conflicting periodicity, so add it ourselves...
          assert(0);
        }
        else {
          if (ii_sign == 1) {
            assert(part->surface->pbi == NULL);
            CTI_Mmap(part->surface->pbi,part->surface->nsp_max);
            if (mpi_rank_shared == 0) {
              // initialize all pbi to the
              for (int it = 0; it < nt; ++it) {
                part->surface->pbi[it]      = BitUtils::packPbiHash(0,it);
                part->surface->pbi[it+nt]   = BitUtils::packPbiHash(0,it+nt); // recall span split into 2
                part->surface->pbi[it+nt*2] = BitUtils::packPbiHash((1<<(2*ii_match)),it);
              }
            }
            MPI_Barrier(mpi_comm_shared);
          }
          else {
            assert(ii_sign == -1);
            assert(0);
          }
        }
      }
    }
    else {
      assert(periodic_type == -1);
    }

    return true;

  }

  bool makeCylinderWithStrand(Part * part,Param * param,int& iarg) {

    bool b_name = false;
    string name;
    bool b_x0 = false;
    double x0[3];
    bool b_x1 = false;
    double x1[3];
    bool b_r = false;
    double r;
    bool b_dn = false;
    double dn; // normal
    bool b_dt = false;
    double dt; // tangential
    bool b_ds = false;
    double ds; // span (== tangential if not specified separately)
    bool b_nn;
    int nn; // normal count: how many spaces to go from dn to dt
    int periodic_type = -1; // 0:x, 1:y, 2:z

    while (iarg < param->size()) {
      const string token = param->getString(iarg++);
      if (token == "NAME") {
        b_name = true;
        name = param->getString(iarg++);
      }
      else if (token == "X0") {
        b_x0 = true;
        FOR_I3 x0[i] = param->getDouble(iarg++);
      }
      else if (token == "X1") {
        b_x1 = true;
        FOR_I3 x1[i] = param->getDouble(iarg++);
      }
      else if (token == "R") {
        b_r = true;
        r = param->getDouble(iarg++);
      }
      else if (token == "DN") {
        b_dn = true;
        dn = param->getDouble(iarg++);
      }
      else if (token == "NN") {
        b_nn = true;
        nn = param->getInt(iarg++);
      }
      else if (token == "DT") {
        b_dt = true;
        dt = param->getDouble(iarg++);
      }
      else if (token == "DS") {
        b_ds = true;
        ds = param->getDouble(iarg++);
      }
      else if (token == "PERIODIC_Z") {
        periodic_type = 2;
      }
      else if (token == "NLAYERS") {
	part->ff_nlayersString = param->getString(iarg++);
      }
      else {
        if (mpi_rank == 0) cout << "unrecognized CYLINDER_WITH_STRAND token: " << token << endl;
      }
    }

    int ierr = 0;
    if (!b_dn) {
      if (mpi_rank == 0) cout << "CYLINDER_WITH_STRAND: missing DN <double>" << endl;
      ierr = -1;
    }
    if (!b_nn) {
      if (mpi_rank == 0) cout << "CYLINDER_WITH_STRAND: missing NN <int>" << endl;
      ierr = -1;
    }
    if (!b_x0) {
      if (mpi_rank == 0) cout << "CYLINDER_WITH_STRAND: missing X0 <x> <y> <z>" << endl;
      ierr = -1;
    }
    if (!b_x1) {
      if (mpi_rank == 0) cout << "CYLINDER_WITH_STRAND: missing X1 <x> <y> <z>" << endl;
      ierr = -1;
    }
    if (!b_r) {
      if (mpi_rank == 0) cout << "CYLINDER_WITH_STRAND: missing R <double>" << endl;
      ierr = -1;
    }
    if (!b_dt) {
      if (mpi_rank == 0) cout << "CYLINDER_WITH_STRAND: missing DT <double>" << endl;
      ierr = -1;
    }
    if (ierr != 0) {
      return false;
    }

    // default for ds is dt if not specified separately...
    if (!b_ds) {
      ds = dt;
    }

    // choose the point count in the tangential direction by using the
    // circumference at r+delta where delta is the full bl thickness...

    const double delta = GeomUtils::stretchingFunction(nn,nn,dn,dt);
    const double L = 2.0*M_PI*(r+delta);
    const int nt = max(int(L/dt),1); // tangential count
    if ((nt < 3)||(nt > 10000000)) {
      if (mpi_rank == 0) cout << "CYLINDER_WITH_STRAND: DT " << dt << " seems unreasonable: circumference L = " << L << endl;
      return false;
    }
    if (mpi_rank == 0) cout << " > nt (~2*pi*(r+delta)/dt) = " << nt << endl;

    const double dx01[3] = DIFF(x1,x0);
    double e1[3],e2[3];
    MiscUtils::getBestE1E2FromE0(e1,e2,dx01);
    double (*xt)[3] = new double[nt*(nn+1)][3];
    for (int it = 0; it < nt; ++it) {
      const double theta = 2.0*M_PI*double(it)/double(nt);
      FOR_I3 xt[it][i] = x0[i] + e1[i]*r*cos(theta) + e2[i]*r*sin(theta);
    }

    // now build the layers...
    for (int it = 0; it < nt; ++it) {
      const int it_prev = (it-1+nt)%nt;
      const int it_next = (it+1)%nt;
      const double dx[3] = DIFF(xt[it_next],xt[it_prev]);
      const double normal[3] = CROSS_PRODUCT(dx,dx01);
      const double normal_mag = MAG(normal);
      assert(normal_mag > 0.0);
      for (int in = 1; in <= nn; ++in) {
        const double delta = GeomUtils::stretchingFunction(in,nn,dn,dt);
        FOR_I3 xt[it+nt*in][i] = xt[it][i] + delta*normal[i]/normal_mag;
      }
    }

    int (*nn_local) = new int[nt];
    for (int it = 0; it < nt; ++it)
      nn_local[it] = nn;

    /*
      const double ratio = 10.5;
      int (*nn_local) = new int[nt];
      for (int it = 0; it < nt; ++it) {
      const int it_prev = (it-1+nt)%nt;
      const int it_next = (it+1)%nt;
      int in;
      for (in = nn; in > 0; --in) {
      const double dt2 = DIST(xt[it_prev+nt*in],xt[it_next+nt*in]);
      if ((dt2 < 2.0*dt*ratio)&&(dt2 > 2.0*dt/ratio))
      break;
      }
      nn_local[it] = in;
      }

      // then ensure there is no single change greater than 1 level...
      bool done = false;
      while (!done) {
      done = true;
      for (int it = 0; it < nt; ++it) {
      const int it_prev = (it-1+nt)%nt;
      const int it_next = (it+1)%nt;
      const int nn_local_nbr_min = min(nn_local[it_prev],nn_local[it_next]);
      if (nn_local[it] > nn_local_nbr_min+1) {
      nn_local[it] = nn_local_nbr_min+1;
      done = false;
      }
      }
      }
    */

    // =======================================
    // the surface...
    // =======================================

    assert(part->surface == NULL);
    part->surface = new SurfaceShm();
    part->surface->init(nt*3,nt*4); // nsp,nst
    if (b_name) {
      // if the user gave the NAME param, then use this to set the zone name...
      part->surface->zoneVec.push_back(StZone(name));
    }
    else {
      // otherwise, use the profile name, or a default...
      part->surface->zoneVec.push_back(StZone("CYLINDER"));
    }
    assert(part->solid_type == UNKNOWN_SOLID);
    part->solid_type = SURFACE_SOLID;

    // recall surface is shm...
    if (mpi_rank_shared == 0) {
      for (int it = 0; it < nt; ++it) {
        FOR_I3 part->surface->xsp[it][i] = xt[it][i];
        FOR_I3 part->surface->xsp[it+nt][i] = xt[it][i] + 0.5*dx01[i]; // split the span into 2
        FOR_I3 part->surface->xsp[it+nt*2][i] = xt[it][i] + dx01[i];
        const int it_next = (it+1)%nt;
        // first pair of tris...
        part->surface->spost[it*4  ][0] = it;
        part->surface->spost[it*4  ][1] = it+nt;
        part->surface->spost[it*4  ][2] = it_next;
        part->surface->znost[it*4  ] = 0;
        part->surface->spost[it*4+1][0] = it_next;
        part->surface->spost[it*4+1][1] = it+nt;
        part->surface->spost[it*4+1][2] = it_next+nt;
        part->surface->znost[it*4+1] = 0;
        // second pair of tris...
        part->surface->spost[it*4+2][0] = it+nt;
        part->surface->spost[it*4+2][1] = it+nt+nt;
        part->surface->spost[it*4+2][2] = it_next+nt;
        part->surface->znost[it*4+2] = 0;
        part->surface->spost[it*4+3][0] = it_next+nt;
        part->surface->spost[it*4+3][1] = it+nt+nt;
        part->surface->spost[it*4+3][2] = it_next+nt+nt;
        part->surface->znost[it*4+3] = 0;
      }
      //GeomUtils::writeTecplot("surface.dat",part->surface->spost,part->surface->nst,part->surface->xsp);
    }
    MPI_Barrier(mpi_comm_shared);

    // =======================================
    // the ff_surface...
    // =======================================

    assert(part->ff_surface == NULL);
    part->ff_surface = new SurfaceShm();
    part->ff_surface->init(nt*4,nt*6); // nsp,nst
    if (b_name) {
      // if the user gave the NAME param, then use this to set the zone name...
      part->ff_surface->zoneVec.push_back(StZone(name+"_FF"));
    }
    else {
      // otherwise, use the profile name, or a default...
      part->ff_surface->zoneVec.push_back(StZone("CYLINDER_FF"));
    }
    assert(part->ff_surface_dxost == NULL);
    CTI_Mmap(part->ff_surface_dxost,nt*6);
    assert(part->ff_type == UNKNOWN_FF);
    part->ff_type = FF_SURFACE_FF;

    if (mpi_rank_shared == 0) {
      for (int it = 0; it < nt; ++it) {
        FOR_I3 part->ff_surface->xsp[it][i]      = xt[it][i];
        FOR_I3 part->ff_surface->xsp[it+nt][i]   = xt[it+nt*nn_local[it]][i];
        FOR_I3 part->ff_surface->xsp[it+2*nt][i] = xt[it+nt*nn_local[it]][i] + dx01[i];
        FOR_I3 part->ff_surface->xsp[it+3*nt][i] = xt[it][i] + dx01[i];
        const int it_next = (it+1)%nt;
        // first pair of tris...
        part->ff_surface->spost[it*6][0] = it;
        part->ff_surface->spost[it*6][2] = it+nt;
        part->ff_surface->spost[it*6][1] = it_next;
        part->ff_surface->znost[it*6] = 0;
        part->ff_surface_dxost[it*6] = dt;
        part->ff_surface->spost[it*6+1][0] = it_next;
        part->ff_surface->spost[it*6+1][2] = it+nt;
        part->ff_surface->spost[it*6+1][1] = it_next+nt;
        part->ff_surface->znost[it*6+1] = 0;
        part->ff_surface_dxost[it*6+1] = dt;
        // second pair of tris...
        part->ff_surface->spost[it*6+2][0] = it+nt;
        part->ff_surface->spost[it*6+2][2] = it+nt+nt;
        part->ff_surface->spost[it*6+2][1] = it_next+nt;
        part->ff_surface->znost[it*6+2] = 0;
        part->ff_surface_dxost[it*6+2] = dt;
        part->ff_surface->spost[it*6+3][0] = it_next+nt;
        part->ff_surface->spost[it*6+3][2] = it+nt+nt;
        part->ff_surface->spost[it*6+3][1] = it_next+nt+nt;
        part->ff_surface->znost[it*6+3] = 0;
        part->ff_surface_dxost[it*6+3] = dt;
        // third pair of tris...
        part->ff_surface->spost[it*6+4][0] = it+nt+nt;
        part->ff_surface->spost[it*6+4][2] = it+nt+nt+nt;
        part->ff_surface->spost[it*6+4][1] = it_next+nt+nt;
        part->ff_surface->znost[it*6+4] = 0;
        part->ff_surface_dxost[it*6+4] = dt;
        part->ff_surface->spost[it*6+5][0] = it_next+nt+nt;
        part->ff_surface->spost[it*6+5][2] = it+nt+nt+nt;
        part->ff_surface->spost[it*6+5][1] = it_next+nt+nt+nt;
        part->ff_surface->znost[it*6+5] = 0;
        part->ff_surface_dxost[it*6+5] = dt;
      }
      //GeomUtils::writeTecplot("ff_surface.dat",part->ff_surface->spost,part->ff_surface->nst,part->ff_surface->xsp);
    }
    MPI_Barrier(mpi_comm_shared);

    // =======================================
    // and finally the points: this part is distributed in parallel...
    // =======================================

    const double dx01_mag = MAG(dx01);
    const int ns = max(1,(int)floor(dx01_mag/ds+0.5));
    assert(part->pts == NULL);
    part->pts = new Points();
    int nt_local = nt/mpi_size;
    if ((nt_local%mpi_size != 0)||(nt_local == 0)) nt_local += 1;
    int it_start = min(nt,nt_local*mpi_rank);
    int it_end = min(nt,nt_local*(mpi_rank+1));
    part->pts->np = 0;
    part->pts->np_global = 0;
    for (int it = 0; it < nt; ++it) {
      // only count our local np...
      if ((it >= it_start)&&(it < it_end))
        part->pts->np += nn_local[it]*ns;
      // all count np_global...
      part->pts->np_global += nn_local[it]*ns;
    }
    if (mpi_rank == 0) cout << " > np_global: " << part->pts->np_global << endl;

    part->pts->xp = new double[part->pts->np][3];
    part->pts->delta = new double[part->pts->np];

    int ip = 0;
    for (int it = it_start; it < it_end; ++it) {
      for (int in = 0; in < nn_local[it]; ++in) {
        double x[3]; FOR_I3 x[i] = 0.5*(xt[it+nt*in][i] + xt[it+nt*(in+1)][i]);
        for (int is = 0; is < ns; ++is) {
          FOR_I3 part->pts->xp[ip][i] = x[i] + (double(is)+0.5)/double(ns)*dx01[i];
          part->pts->delta[ip] = 2.0*dt; // TODO: could be more precise here
          ++ip;
        }
      }
    }
    assert(ip == part->pts->np);

    //GeomUtils::writePtsTecplot("cyl_points.dat",part->pts->xp,part->pts->np);

    /*
    // take a look...
    FILE * fp = fopen("spline.dat","w");
    for (int it = 0; it < nt; ++it) {
    for (int in = 0; in <= nn_local[it]; ++in) {
    fprintf(fp,"%18.15le %18.15le %18.15le\n",xt[it+nt*in][0],xt[it+nt*in][1],xt[it+nt*in][2]);
    }
    }
    fclose(fp);
    */

    delete[] nn_local;
    delete[] xt;

    // set any periodicity...
    if (periodic_type == 2) {
      // PERIODIC_Z...
      // confirm that the X0->X1 vector is a z-translation...
      if ((dx01[0] != 0.0)||(dx01[1] != 0.0)||(dx01[2] == 0.0)) {
        if (mpi_rank == 0) cout << "WARNING: dx form X0 to X1 not consistent with PERIODIC_Z: " << COUT_VEC(dx01) << endl;
        ierr = -1;
      }
      else {
        // check if there is matching periodic data...
        int ii_match = -1;
        int ii_sign = 0;
        for (int ii = 0; ii < PeriodicData::periodicTransformVec.size(); ++ii) {
          const int kind = PeriodicData::periodicTransformVec[ii].getKind();
          //cout << " > got kind: " << kind << " PERIODIC_TRANSFORM_CART: " << PERIODIC_TRANSFORM_CART << endl;
          if (kind == PERIODIC_TRANSFORM_CART) {
            double data[3];
            PeriodicData::periodicTransformVec[ii].getData(data);
            if ((data[0] == 0.0)&&(data[1] == 0.0)&&(data[2] != 0.0)) {
              if (data[2] == dx01[2]) {
                ii_match = ii;
                ii_sign = 1;
              }
              else if (data[2] == -dx01[2]) {
                ii_match = ii;
                ii_sign = -1;
              }
              else {
                if (mpi_rank == 0) cout << "WARNING: existing z-periodicity does not match: " << COUT_VEC(data) << endl;
                ierr = -1;
              }
              break;
            }
          }
        }
        if (ii_match == -1) {
          // no matching or conflicting periodicity, so add it ourselves...
          assert(0);
        }
        else {
          if (ii_sign == 1) {
            assert(part->surface->pbi == NULL);
            CTI_Mmap(part->surface->pbi,part->surface->nsp_max);
            if (mpi_rank_shared == 0) {
              // initialize all pbi to the
              for (int it = 0; it < nt; ++it) {
                part->surface->pbi[it]      = BitUtils::packPbiHash(0,it);
                part->surface->pbi[it+nt]   = BitUtils::packPbiHash(0,it+nt); // recall span split into 2
                part->surface->pbi[it+nt*2] = BitUtils::packPbiHash((1<<(2*ii_match)),it);
              }
            }
            MPI_Barrier(mpi_comm_shared);
          }
          else {
            assert(ii_sign == -1);
            assert(0);
          }
        }
      }
    }
    else {
      assert(periodic_type == -1);
    }

    return true;

  }

  int makeBoxWithCartPts(Part * part,Param * param,int& iarg,const bool pts_only) {

    string part_name;
    if (pts_only) part_name = "BOX_WITH_CART_PTS";
    else part_name = "BOX_WITH_CART";

    bool b_x0 = false;
    double x0[3];
    bool b_x1 = false;
    double x1[3];
    bool b_n = false;
    int n[3];
    bool periodic[3] = {false,false,false};

    if (mpi_rank == 0) cout << " > making " << part_name << "..." << endl;

    while (iarg < param->size()) {
      string token = param->getString(iarg++);
      if (token == "X0") {
        b_x0 = true;
        FOR_I3 x0[i] = param->getDouble(iarg++);
        if (mpi_rank == 0) cout << " > X0: " << COUT_VEC(x0) << endl;
      }
      else if (token == "X1") {
        b_x1 = true;
        FOR_I3 x1[i] = param->getDouble(iarg++);
        if (mpi_rank == 0) cout << " > X1: " << COUT_VEC(x1) << endl;
      }
      else if (token == "N") {
        b_n = true;
        FOR_I3 n[i] = param->getInt(iarg++);
        if (mpi_rank == 0) cout << " > N: " << COUT_VEC(n) << ", " << n[0]*n[1]*n[2] << " pts" << endl;
      }
      else if (token == "PERIODIC_X") {
        periodic[0] = true;
      }
      else if (token == "PERIODIC_Y") {
        periodic[1] = true;
      }
      else if (token == "PERIODIC_Z") {
        periodic[2] = true;
      }
      else {
        CWARN("Unknown " << part_name << " parameter " << token << "; skipping");
      }
    }

    int ierr = 0;
    if (!b_x0) {
      COUT1("Error: " << part_name << " missing X0 <x> <y> <z>");
      ierr = -1;
    }
    if (!b_x1) {
      COUT1("Error: " << part_name << " missing X1 <x> <y> <z>");
      ierr = -1;
    }
    if (!b_n) {
      COUT1("Error: " << part_name << " missing point counts N <nx> <ny> <nz>");
      ierr = -1;
    }

    if (ierr == 0) {
      if (pts_only) {
        // surface...
        assert(part->surface == NULL);
        assert(part->solid_type == UNKNOWN_SOLID);
        part->solid_type = NO_SOLID;

        // ff_surface...
        assert(part->ff_type == UNKNOWN_FF);
        part->ff_type = BOX_FF;
        assert(part->ff_surface == NULL);
        part->ff_surface = new SurfaceShm();
        FOR_I3 periodic[i] = false;  // no periodicity for pts only part
        part->ff_surface->makeBox(x0,x1,periodic);
        FOR_I3 part->ddata_ff_[i] = x0[i];
        FOR_I3 part->ddata_ff_[i+3] = x1[i];
        assert(part->ff_surface_dxost == NULL);
        CTI_Mmap(part->ff_surface_dxost,part->ff_surface->nst);
        if (mpi_rank_shared == 0) {
          const double dx = (x1[0]-x0[0])/double(n[0]);
          const double dy = (x1[1]-x0[1])/double(n[1]);
          const double dz = (x1[2]-x0[2])/double(n[2]);
          const double dxost = max(dx,max(dy,dz));
          for (int ist = 0; ist < part->ff_surface->nst; ++ist) {
            part->ff_surface_dxost[ist] = dxost;
          }
        }
      }
      else {
        // surface...
        assert(part->surface == NULL);
        part->surface = new SurfaceShm();
        part->surface->makeBox(x0,x1,periodic);
        assert(part->solid_type == UNKNOWN_SOLID);
        part->solid_type = BOX_SOLID;
        FOR_I3 part->ddata_solid_[i] = x0[i];
        FOR_I3 part->ddata_solid_[i+3] = x1[i];


        // ff_surface...
        assert(part->ff_type == UNKNOWN_FF);
        part->ff_type = NO_FF;
        assert(part->ff_surface == NULL);
      }
      MPI_Barrier(mpi_comm_shared);

      // pts same for both...
      assert(part->pts == NULL);
      part->pts = new Points();
      part->pts->makeCart(x0,x1,n);
      assert(part->pts->np_global == n[0]*n[1]*n[2]);
    }

    return ierr;

  }

  int makeAnnularCylPts(Part * part,Param * param,int& iarg,const bool pts_only) {

    string part_name;
    if (pts_only) part_name = "ANNULAR_CYL_PTS";
    else part_name = "ANNULAR_CYL";

    bool b_x0 = false;
    double x0[3];
    bool b_x1 = false;
    double x1[3];
    bool b_r0 = false;
    double r0;
    bool b_r1 = false;
    double r1;
    bool b_dl = false;
    double dl;
    bool b_dr = false;
    double dr;
    bool b_dt = false;
    double dt;
    bool b_n = false;
    int n[3];
    bool b_nl = false;
    int nl;
    bool b_nr = false;
    int nr;
    bool b_nt = false;
    int nt;

    if (mpi_rank == 0) cout << " > making " << part_name << "..." << endl;

    int ierr = 0;
    while (iarg < param->size()) {
      string token = param->getString(iarg++);
      if (token == "X0") {
        b_x0 = true;
        FOR_I3 x0[i] = param->getDouble(iarg++);
        if (mpi_rank == 0) cout << " > X0: " << COUT_VEC(x0) << endl;
      }
      else if (token == "X1") {
        b_x1 = true;
        FOR_I3 x1[i] = param->getDouble(iarg++);
        if (mpi_rank == 0) cout << " > X1: " << COUT_VEC(x1) << endl;
      }
      else if (token == "R0") {
        b_r0 = true;
        r0 = param->getDouble(iarg++);
        if (mpi_rank == 0) cout << " > R0: " << r0 << endl;
      }
      else if (token == "R1") {
        b_r1 = true;
        r1 = param->getDouble(iarg++);
        if (mpi_rank == 0) cout << " > R1: " << r1 << endl;
      }
      else if ((token == "DL")||(token == "DELTA")||(token == "DX")) {
        b_dl = true;
        dl = param->getDouble(iarg++);
        if (mpi_rank == 0) cout << " > DL (or DELTA or DX): " << dl << endl;
      }
      else if (token == "DR") {
        b_dr = true;
        dr = param->getDouble(iarg++);
        if (mpi_rank == 0) cout << " > DR: " << dr << endl;
      }
      else if (token == "DT") {
        b_dt = true;
        dt = param->getDouble(iarg++);
        if (mpi_rank == 0) cout << " > DT: " << dt << endl;
      }
      else if (token == "N") {
        b_n = true;
        FOR_I3 n[i] = param->getInt(iarg++);
        if (mpi_rank == 0) cout << " > N: " << COUT_VEC(n) << endl;
      }
      else if ((token == "NL")||(token == "NX")) {
        b_nl = true;
        nl = param->getInt(iarg++);
        if (mpi_rank == 0) cout << " > NL (or NX): " << nl << endl;
      }
      else if (token == "NR") {
        b_nr = true;
        nr = param->getInt(iarg++);
        if (mpi_rank == 0) cout << " > NR: " << nr << endl;
      }
      else if (token == "NT") {
        b_nt = true;
        nt = param->getInt(iarg++);
        if (mpi_rank == 0) cout << " > NT: " << nt << endl;
      }
      else if (token == "NLAYERS") {
        part->ff_nlayersString = param->getString(iarg++);
        if (mpi_rank == 0) cout << " > NLAYERS: " << part->ff_nlayersString << endl;
      }
      else {
        CWARN("Unknown " << part_name << " parameter " << token);
      }
    }

    if (!b_x0) {
      COUT1("Error: " << part_name << " missing X0 <x> <y> <z>");
      ierr = -1;
    }
    if (!b_x1) {
      COUT1("Error: " << part_name << " missing X1 <x> <y> <z>");
      ierr = -1;
    }
    if (!b_r0) {
      COUT1("Error: " << part_name << " missing R0 <double>");
      ierr = -1;
    }
    if (!b_r1) {
      COUT1("Error: " << part_name << " missing R1 <double>");
      ierr = -1;
    }
    // the routine below needs n[3] set...
    if (!b_n) {
      // L direction...
      if (b_nl) {
        n[0] = nl;
      }
      else if (b_dl) {
        n[0] = (int)ceil(DIST(x0,x1)/dl);
      }
      else {
        if (mpi_rank == 0) cout << "Error: " << part_name << " missing N <double> <double> <double> or NL <double>" << endl;
        ierr = -1;
      }
      // R direction...
      if (b_nr) {
        n[1] = nr;
      }
      else if (b_dr) {
        n[1] = (int)ceil((r1-r0)/dr);
      }
      else if (b_dl) {
        n[1] = (int)ceil((r1-r0)/dl);
      }
      else {
        if (mpi_rank == 0) cout << "Error: " << part_name << " missing N <double> <double> <double> or NR <double>" << endl;
        ierr = -1;
      }
      // theta direction...
      if (b_nt) {
        n[2] = nt;
      }
      else if (b_dt) {
        n[2] = (int)ceil((r0+r1)*M_PI/dt);
      }
      else if (b_dl) {
        n[2] = (int)ceil((r0+r1)*M_PI/dl);
      }
      else {
        if (mpi_rank == 0) cout << "Error: " << part_name << " missing N <double> <double> <double> or NT <double>" << endl;
        ierr = -1;
      }
    }

    if (ierr == 0) {

      if (mpi_rank == 0) cout << " > final N (NL,NR,NT): " << COUT_VEC(n) << ", " << n[0]*n[1]*n[2] << " pts" << endl;

      if (pts_only) {
        part->solid_type = NO_SOLID;

        part->ff_type = ANNULAR_CYLINDER_FF;
        part->ddata_ff_[0] = x0[0];
        part->ddata_ff_[1] = x0[1];
        part->ddata_ff_[2] = x0[2];
        part->ddata_ff_[3] = x1[0];
        part->ddata_ff_[4] = x1[1];
        part->ddata_ff_[5] = x1[2];
        part->ddata_ff_[6] = r0;
        part->ddata_ff_[7] = r1;
        // surface...
        assert(part->surface == NULL);
        // ff_surface...
        assert(part->ff_surface == NULL);
        part->ff_surface = new SurfaceShm();
        part->ff_surface->makeAnnularCyl("FF",x0,x1,r0,r1,n[2],false); // for the number of segments, use the ntheta
        assert(part->ff_surface_dxost == NULL);
        CTI_Mmap(part->ff_surface_dxost,part->ff_surface->nst);
        if (mpi_rank_shared == 0) {
          const double dl = DIST(x0,x1)/double(n[0]);
          const double dr = (r1-r0)/double(n[1]);
          const double dt = (r0+r1)*M_PI/double(n[2]);
          const double dxost = max(dl,max(dr,dt));
          for (int ist = 0; ist < part->ff_surface->nst; ++ist) {
            part->ff_surface_dxost[ist] = dxost;
          }
        }
      }
      else {
        part->solid_type = ANNULAR_SOLID;
        part->ddata_solid_[0] = x0[0];
        part->ddata_solid_[1] = x0[1];
        part->ddata_solid_[2] = x0[2];
        part->ddata_solid_[3] = x1[0];
        part->ddata_solid_[4] = x1[1];
        part->ddata_solid_[5] = x1[2];
        part->ddata_solid_[6] = r0;
        part->ddata_solid_[7] = r1;
        // surface...
        assert(part->surface == NULL);
        part->surface = new SurfaceShm();
        part->surface->makeAnnularCyl("annulus",x0,x1,r0,r1,n[2],false); // for the number of segments, use the ntheta

        // ff_surface...
        assert(part->ff_surface == NULL);
        part->ff_type = NO_FF;
      }
      MPI_Barrier(mpi_comm_shared);

      // pts...
      assert(part->pts == NULL);
      part->pts = new Points();
      part->pts->makeAnnularCyl(x0,x1,r0,r1,n);
      assert(part->pts->np_global == n[0]*n[1]*n[2]);
    }

    return ierr;

  }

  int makeAnnularDlrtPts(Part * part,Param * param,int& iarg,const bool pts_only) {

    string part_name;
    if (pts_only) part_name = "ANNULAR_DLRT_PTS";
    else part_name = "ANNULAR_DLRT";

    bool b_x0 = false;
    double x0[3];
    bool b_x1 = false;
    double x1[3];
    bool b_r0 = false;
    double r0;
    bool b_r1 = false;
    double r1;
    bool b_d = false;
    double d[3];
    bool b_dl = false;
    bool b_dr = false;
    bool b_dt = false;
    bool b_nl = false;
    int nl;
    bool b_nr = false;
    int nr;

    if (mpi_rank == 0) cout << " > making " << part_name << "..." << endl;

    int ierr = 0;
    while (iarg < param->size()) {
      string token = param->getString(iarg++);
      if (token == "X0") {
        b_x0 = true;
        FOR_I3 x0[i] = param->getDouble(iarg++);
        if (mpi_rank == 0) cout << " > X0: " << COUT_VEC(x0) << endl;
      }
      else if (token == "X1") {
        b_x1 = true;
        FOR_I3 x1[i] = param->getDouble(iarg++);
        if (mpi_rank == 0) cout << " > X1: " << COUT_VEC(x1) << endl;
      }
      else if (token == "R0") {
        b_r0 = true;
        r0 = param->getDouble(iarg++);
        if (mpi_rank == 0) cout << " > R0: " << r0 << endl;
      }
      else if (token == "R1") {
        b_r1 = true;
        r1 = param->getDouble(iarg++);
        if (mpi_rank == 0) cout << " > R1: " << r1 << endl;
      }
      else if ((token == "DL")||(token == "DELTA")||(token == "DX")) {
        b_dl = true;
        d[0] = param->getDouble(iarg++);
        if (mpi_rank == 0) cout << " > DL (or DELTA or DX): " << d[0] << endl;
      }
      else if (token == "DR") {
        b_dr = true;
        d[1] = param->getDouble(iarg++);
        if (mpi_rank == 0) cout << " > DR: " << d[1] << endl;
      }
      else if (token == "DT") {
        b_dt = true;
        d[2] = param->getDouble(iarg++);
        if (mpi_rank == 0) cout << " > DT: " << d[2] << endl;
      }
      else if (token == "D") {
        b_d = true;
        FOR_I3 d[i] = param->getDouble(iarg++);
        if (mpi_rank == 0) cout << " > D: " << COUT_VEC(d) << endl;
      }
      else if ((token == "NL")||(token == "NX")) {
        b_nl = true;
        nl = param->getInt(iarg++);
        if (mpi_rank == 0) cout << " > NL (or NX): " << nl << endl;
      }
      else if (token == "NR") {
        b_nr = true;
        nr = param->getInt(iarg++);
        if (mpi_rank == 0) cout << " > NR: " << nr << endl;
      }
      else if (token == "NLAYERS") {
        part->ff_nlayersString = param->getString(iarg++);
        if (mpi_rank == 0) cout << " > NLAYERS: " << part->ff_nlayersString << endl;
      }
      else {
        CWARN("Unknown " << part_name << " parameter " << token);
      }
    }

    if (!b_x0) {
      COUT1("Error: " << part_name << " missing X0 <x> <y> <z>");
      ierr = -1;
    }
    if (!b_x1) {
      COUT1("Error: " << part_name << " missing X1 <x> <y> <z>");
      ierr = -1;
    }
    if (!b_r0) {
      COUT1("Error: " << part_name << " missing R0 <double>");
      ierr = -1;
    }
    if (!b_r1) {
      COUT1("Error: " << part_name << " missing R1 <double>");
      ierr = -1;
    }
    // the routine below needs d[3] set...
    if (!b_d) {
      // L direction...
      if (!b_dl) {
        if (b_nl) {
          const double L=DIST(x0,x1);
          d[0] = L/double(nl);
        }
        else {
          if (mpi_rank == 0) cout << "Error: " << part_name << " missing D <double> <double> <double> or DL <double>" << endl;
          ierr = -1;
        }
      }
      // R direction...
      if (!b_dr) {
        if (b_nr) {
          d[1] = (r1-r0)/double(nr);
        }
        else {
          if (mpi_rank == 0) cout << "Error: " << part_name << " missing D <double> <double> <double> or DR <double>" << endl;
          ierr = -1;
        }
      }
      // theta direction...
      if (!b_dt) {
        if (mpi_rank == 0) cout << "Error: " << part_name << " missing D <double> <double> <double> or DT <double>" << endl;
        ierr = -1;
      }
    }

    if (ierr == 0) {

      if (pts_only) {
        part->solid_type = NO_SOLID;

        part->ff_type = ANNULAR_CYLINDER_FF;
        part->ddata_ff_[0] = x0[0];
        part->ddata_ff_[1] = x0[1];
        part->ddata_ff_[2] = x0[2];
        part->ddata_ff_[3] = x1[0];
        part->ddata_ff_[4] = x1[1];
        part->ddata_ff_[5] = x1[2];
        part->ddata_ff_[6] = r0;
        part->ddata_ff_[7] = r1;
        // surface...
        assert(part->surface == NULL);
        // ff_surface...
        assert(part->ff_surface == NULL);
        part->ff_surface = new SurfaceShm();
        const int ntheta = int(ceil(2.0*M_PI*r1/d[2]));
        part->ff_surface->makeAnnularCyl("FF",x0,x1,r0,r1,ntheta,false); // for the number of segments, use the ntheta
        assert(part->ff_surface_dxost == NULL);
        CTI_Mmap(part->ff_surface_dxost,part->ff_surface->nst);
        if (mpi_rank_shared == 0) {
          const double dxost = max(d[0],max(d[1],d[2]));
          for (int ist = 0; ist < part->ff_surface->nst; ++ist) {
            part->ff_surface_dxost[ist] = dxost;
          }
        }
      }
      else {
        part->solid_type = ANNULAR_SOLID;
        part->ddata_solid_[0] = x0[0];
        part->ddata_solid_[1] = x0[1];
        part->ddata_solid_[2] = x0[2];
        part->ddata_solid_[3] = x1[0];
        part->ddata_solid_[4] = x1[1];
        part->ddata_solid_[5] = x1[2];
        part->ddata_solid_[6] = r0;
        part->ddata_solid_[7] = r1;
        // surface...
        assert(part->surface == NULL);
        part->surface = new SurfaceShm();
        const int ntheta = int(ceil(2.0*M_PI*r1/d[2]));
        part->surface->makeAnnularCyl("annulus",x0,x1,r0,r1,ntheta,false); // for the number of segments, use the ntheta

        // ff_surface...
        assert(part->ff_surface == NULL);
        part->ff_type = NO_FF;
      }
      MPI_Barrier(mpi_comm_shared);

      // pts...
      assert(part->pts == NULL);
      part->pts = new Points();
      double radii[2] = {r0,r1};
      part->pts->makeDlrt(x0,x1,radii,radii,d);
    }

    return ierr;

  }

  int makeCylWithHexcorePts(Part * part,Param * param,int& iarg,const bool pts_only) {

    string part_name;
    if (pts_only) part_name = "CYL_HEXCORE_PTS";
    else part_name = "CYL_HEXCORE";

    bool b_x0 = false;
    double x0[3];
    bool b_x1 = false;
    double x1[3];
    bool b_r = false;
    double r;
    bool b_n = false;
    int n[3];
    bool b_nl = false;
    int nl;
    bool b_nr = false;
    int nr;
    bool b_nt = false;
    int nt;

    if (mpi_rank == 0) cout << " > making " << part_name << "..." << endl;

    int ierr = 0;
    while (iarg < param->size()) {
      string token = param->getString(iarg++);
      if (token == "X0") {
        b_x0 = true;
        FOR_I3 x0[i] = param->getDouble(iarg++);
        if (mpi_rank == 0) cout << " > X0: " << COUT_VEC(x0) << endl;
      }
      else if (token == "X1") {
        b_x1 = true;
        FOR_I3 x1[i] = param->getDouble(iarg++);
        if (mpi_rank == 0) cout << " > X1: " << COUT_VEC(x1) << endl;
      }
      else if (token == "R") {
        b_r = true;
        r = param->getDouble(iarg++);
        if (mpi_rank == 0) cout << " > R: " << r << endl;
      }
      else if (token == "N") {
        b_n = true;
        FOR_I3 n[i] = param->getInt(iarg++);
        if (mpi_rank == 0) cout << " > N: " << COUT_VEC(n) << endl;
      }
      else if ((token == "NL")||(token == "NX")) {
        b_nl = true;
        nl = param->getInt(iarg++);
        if (mpi_rank == 0) cout << " > NL (or NX): " << nl << endl;
      }
      else if (token == "NR") {
        b_nr = true;
        nr = param->getInt(iarg++);
        if (mpi_rank == 0) cout << " > NR: " << nr << endl;
      }
      else if (token == "NT") {
        b_nt = true;
        nt = param->getInt(iarg++);
        if (mpi_rank == 0) cout << " > NT: " << nt << endl;
      }
      else if (token == "NLAYERS") {
        part->ff_nlayersString = param->getString(iarg++);
        if (mpi_rank == 0) cout << " > NLAYERS: " << part->ff_nlayersString << endl;
      }
      else {
        CWARN("Unknown " << part_name << " parameter " << token);
      }
    }

    if (!b_x0) {
      COUT1("Error: " << part_name << " missing X0 <x> <y> <z>");
      ierr = -1;
    }
    if (!b_x1) {
      COUT1("Error: " << part_name << " missing X1 <x> <y> <z>");
      ierr = -1;
    }
    if (!b_r) {
      COUT1("Error: " << part_name << " missing R <double>");
      ierr = -1;
    }
    // the routine below needs n[3] set...
    if (!b_n) {
      // L direction...
      if (b_nl) {
        n[0] = nl;
      }
      else {
        if (mpi_rank == 0) cout << "Error: " << part_name << " missing N <double> <double> <double> or NL <double>" << endl;
        ierr = -1;
      }
      // R direction...
      if (b_nr) {
        n[1] = nr;
      }
      else {
        if (mpi_rank == 0) cout << "Error: " << part_name << " missing N <double> <double> <double> or NR <double>" << endl;
        ierr = -1;
      }
      // theta direction...
      if (b_nt) {
        n[2] = nt;
      }
      else {
        if (mpi_rank == 0) cout << "Error: " << part_name << " missing N <double> <double> <double> or NT <double>" << endl;
        ierr = -1;
      }
    }

    if (ierr == 0) {

      if (pts_only) {
        part->solid_type = NO_SOLID;

        part->ff_type = CYLINDER_FF;
        part->ddata_ff_[0] = x0[0];
        part->ddata_ff_[1] = x0[1];
        part->ddata_ff_[2] = x0[2];
        part->ddata_ff_[3] = x1[0];
        part->ddata_ff_[4] = x1[1];
        part->ddata_ff_[5] = x1[2];
        part->ddata_ff_[6] = r;
        // surface...
        assert(part->surface == NULL);
        // ff_surface...
        assert(part->ff_surface == NULL);
        part->ff_surface = new SurfaceShm();
        part->ff_surface->makeCylinder("FF",x0,x1,r,n[2],true);
        assert(part->ff_surface_dxost == NULL);
        CTI_Mmap(part->ff_surface_dxost,part->ff_surface->nst);
        if (mpi_rank_shared == 0) {
          const double dl = DIST(x0,x1)/double(n[0]);
          const double dr = r/double(n[1]);
          const double dt = r*M_PI/double(n[2]);
          const double dxost = max(dl,max(dr,dt));
          for (int ist = 0; ist < part->ff_surface->nst; ++ist) {
            part->ff_surface_dxost[ist] = dxost;
          }
        }
      }
      else {
        part->solid_type = CYLINDER_SOLID;
        part->ddata_solid_[0] = x0[0];
        part->ddata_solid_[1] = x0[1];
        part->ddata_solid_[2] = x0[2];
        part->ddata_solid_[3] = x1[0];
        part->ddata_solid_[4] = x1[1];
        part->ddata_solid_[5] = x1[2];
        part->ddata_solid_[6] = r;
        // surface...
        assert(part->surface == NULL);
        part->surface = new SurfaceShm();
        part->surface->makeCylinder("cylinder",x0,x1,r,n[2],true);

        // ff_surface...
        assert(part->ff_surface == NULL);
        part->ff_type = NO_FF;
      }
      MPI_Barrier(mpi_comm_shared);

      // pts...
      assert(part->pts == NULL);
      part->pts = new Points();
      part->pts->makeHexcore(x0,x1,r,n);
    }

    return ierr;

  }

  int makeCylWithDlrtPts(Part * part,Param * param,int& iarg,const bool pts_only) {

    string part_name;
    if (pts_only) part_name = "CYL_DLRT_PTS";
    else part_name = "CYL_DLRT";

    bool b_x0 = false;
    double x0[3];
    bool b_x1 = false;
    double x1[3];
    bool b_r = false;
    double r;
    bool b_d = false;
    double d[3];
    bool b_nl = false;
    int nl;
    bool b_nr = false;
    int nr;
    bool b_dl = false;
    bool b_dr = false;
    bool b_dt = false;

    if (mpi_rank == 0) cout << " > making " << part_name << "..." << endl;

    int ierr = 0;
    while (iarg < param->size()) {
      string token = param->getString(iarg++);
      if (token == "X0") {
        b_x0 = true;
        FOR_I3 x0[i] = param->getDouble(iarg++);
        if (mpi_rank == 0) cout << " > X0: " << COUT_VEC(x0) << endl;
      }
      else if (token == "X1") {
        b_x1 = true;
        FOR_I3 x1[i] = param->getDouble(iarg++);
        if (mpi_rank == 0) cout << " > X1: " << COUT_VEC(x1) << endl;
      }
      else if (token == "R") {
        b_r = true;
        r = param->getDouble(iarg++);
        if (mpi_rank == 0) cout << " > R: " << r << endl;
      }
      else if (token == "D") {
        b_d = true;
        FOR_I3 d[i] = param->getDouble(iarg++);
        if (mpi_rank == 0) cout << " > D: " << COUT_VEC(d) << endl;
      }
      else if ((token == "NL")||(token == "NX")) {
        b_nl = true;
        nl = param->getInt(iarg++);
        if (mpi_rank == 0) cout << " > NL (or NX): " << nl << endl;
      }
      else if (token == "NR") {
        b_nr = true;
        nr = param->getInt(iarg++);
        if (mpi_rank == 0) cout << " > NR: " << nr << endl;
      }
      else if (token == "DL" || token == "DX") {
        b_dl = true;
        d[0] = param->getDouble(iarg++);
        if (mpi_rank == 0) cout << " > DL: " << d[0] << endl;
      }
      else if (token == "DR") {
        b_dr = true;
        d[1] = param->getDouble(iarg++);
        if (mpi_rank == 0) cout << " > DR: " << d[1] << endl;
      }
      else if (token == "DT") {
        b_dt = true;
        d[2] = param->getDouble(iarg++);
        if (mpi_rank == 0) cout << " > DT: " << d[2] << endl;
      }
      else if (token == "NLAYERS") {
        part->ff_nlayersString = param->getString(iarg++);
        if (mpi_rank == 0) cout << " > NLAYERS: " << part->ff_nlayersString << endl;
      }
      else {
        CWARN("Unknown " << part_name << " parameter " << token);
      }
    }

    if (!b_x0) {
      COUT1("Error: " << part_name << " missing X0 <x> <y> <z>");
      ierr = -1;
    }
    if (!b_x1) {
      COUT1("Error: " << part_name << " missing X1 <x> <y> <z>");
      ierr = -1;
    }
    if (!b_r) {
      COUT1("Error: " << part_name << " missing R <double>");
      ierr = -1;
    }
    // the routine below needs d[3] set...
    if (!b_d) {
      // L direction...
      if (!b_dl) {
        if (b_nl) {
          const double L = DIST(x0,x1);
          d[0] = L/double(nl);
        }
        else {
          if (mpi_rank == 0) cout << "Error: " << part_name << " missing D <double> <double> <double> or DL <double>" << endl;
          ierr = -1;
        }
      }
      // R direction...
      if (!b_dr) {
        if (b_nr) {
          d[1] = r/double(nr);
        }
        else {
          if (mpi_rank == 0) cout << "Error: " << part_name << " missing D <double> <double> <double> or DR <double>" << endl;
          ierr = -1;
        }
      }
      // theta direction...
      if (!b_dt) {
        if (mpi_rank == 0) cout << "Error: " << part_name << " missing D <double> <double> <double> or DT <double>" << endl;
        ierr = -1;
      }
    }

    if (ierr == 0) {

      if (pts_only) {
        part->solid_type = NO_SOLID;

        part->ff_type = CYLINDER_FF;
        part->ddata_ff_[0] = x0[0];
        part->ddata_ff_[1] = x0[1];
        part->ddata_ff_[2] = x0[2];
        part->ddata_ff_[3] = x1[0];
        part->ddata_ff_[4] = x1[1];
        part->ddata_ff_[5] = x1[2];
        part->ddata_ff_[6] = r;
        // surface...
        assert(part->surface == NULL);
        // ff_surface...
        assert(part->ff_surface == NULL);
        part->ff_surface = new SurfaceShm();
        const int ntheta = int(ceil(r*2.0*M_PI/d[2]));
        part->ff_surface->makeCylinder("FF",x0,x1,r,ntheta,true);
        assert(part->ff_surface_dxost == NULL);
        CTI_Mmap(part->ff_surface_dxost,part->ff_surface->nst);
        if (mpi_rank_shared == 0) {
          const double dxost = max(d[0],max(d[1],d[2]));
          for (int ist = 0; ist < part->ff_surface->nst; ++ist) {
            part->ff_surface_dxost[ist] = dxost;
          }
        }
      }
      else {
        part->solid_type = CYLINDER_SOLID;
        part->ddata_solid_[0] = x0[0];
        part->ddata_solid_[1] = x0[1];
        part->ddata_solid_[2] = x0[2];
        part->ddata_solid_[3] = x1[0];
        part->ddata_solid_[4] = x1[1];
        part->ddata_solid_[5] = x1[2];
        part->ddata_solid_[6] = r;
        // surface...
        assert(part->surface == NULL);
        part->surface = new SurfaceShm();
        const int ntheta = int(ceil(r*2.0*M_PI/d[2]));
        part->surface->makeCylinder("cylinder",x0,x1,r,ntheta,true);

        // ff_surface...
        assert(part->ff_surface == NULL);
        part->ff_type = NO_FF;
      }
      MPI_Barrier(mpi_comm_shared);

      // pts...
      assert(part->pts == NULL);
      part->pts = new Points();
      double rad0[2] = {0.0,r};
      part->pts->makeDlrt(x0,x1,rad0,rad0,d);
    }

    return ierr;

  }

  int makeSphereWithBl(Part * part,Param * param,int& iarg) {

    // SPHERE_WITH_BL XC=<x> <y> <z> R=<r> DN=<dn> DT=<dt> N=<n>

    double xc[3] = { 0.0, 0.0, 0.0 };
    bool b_r = false;
    double r;
    bool b_dn = false;
    double dn;
    bool b_dt = false;
    double dt_target;
    bool b_n = false;
    int n;

    if (mpi_rank == 0) cout << " > making SPHERE_WITH_BL..." << endl;

    while (iarg < param->size()) {
      string token = param->getString(iarg++);
      if ((token == "X")||(token == "XC")) {
        FOR_I3 xc[i] = param->getDouble(iarg++);
        if (mpi_rank == 0) cout << " > XC: " << COUT_VEC(xc) << endl;
      }
      else if (token == "R") {
        b_r = true;
        r = param->getDouble(iarg++);
        if (mpi_rank == 0) cout << " > R: " << r << endl;
        // once we have everything, break...
        if (b_r && b_dn && b_dt && b_n)
          break;
      }
      else if (token == "DN") {
        b_dn = true;
        dn = param->getDouble(iarg++);
        if (mpi_rank == 0) cout << " > DN: " << dn << endl;
        // once we have everything, break...
        if (b_r && b_dn && b_dt && b_n)
          break;
      }
      else if (token == "DT") {
        b_dt = true;
        dt_target = param->getDouble(iarg++);
        if (mpi_rank == 0) cout << " > DT: " << dt_target << endl;
        // once we have everything, break...
        if (b_r && b_dn && b_dt && b_n)
          break;
      }
      else if (token == "N") {
        b_n = true;
        n = param->getInt(iarg++);
        if (mpi_rank == 0) cout << " > N: " << n << endl;
        // once we have everything, break...
        if (b_r && b_dn && b_dt && b_n)
          break;
      }
      else {
        CWARN("Unknown SPHERE_WITH_BL parameter " << token << "; skipping");
      }
    }

    int ierr = 0;
    if (!b_r) {
      COUT1("Error: SPHERE_WITH_BL missing R <double>");
      ierr = -1;
    }
    if (!b_dn) {
      COUT1("Error: SPHERE_WITH_BL missing DN <double>");
      ierr = -1;
    }
    if (!b_dt) {
      COUT1("Error: SPHERE_WITH_BL missing DT <double>");
      ierr = -1;
    }
    if (!b_n) {
      COUT1("Error: SPHERE_WITH_BL missing point counts N <nx> <ny> <nz>");
      ierr = -1;
    }

    if (ierr == 0) {

      // the edge segmentation...
      const int nseg = max(1,int(1.208*r/dt_target));
      const double dt = 1.208*r/double(nseg);
      if (mpi_rank == 0) cout << " > sphere edge segmentation: " << nseg << ", actual DT: " << dt << endl;

      // now do some BL calcs...
      // TODO: used should probably set the FF dt surface and the

      // for the ff calculation, use the constant stretching function
      // formulation. We need to solve the coupled non-linear
      // problem:
      //
      // (eq1) delta = (dtff*sf-dn)/(sf-1.0)
      // (eq2) dtff = 1.208*(r+delta)/double(nseg)
      //
      // where dtff is the tangential AND normal spacing (last BL cell is isotropic)
      // at the maximum radius, i.e. the FF...

      double delta = 0.5*double(n)*(dn+dt);
      double dtff = 1.208*(r+delta)/double(nseg);
      double sf,ddtff;
      int iter = 0;
      do {
        ++iter;
        assert(iter < 100);
        sf = pow(dtff/dn,1.0/double(n-1));
        delta = (dtff*sf-dn)/(sf-1.0);
        ddtff = 1.208*(r+delta)/double(nseg) - dtff;
        dtff += ddtff;
      }
      while (fabs(ddtff) > 1.0E-10*dtff);

      if (mpi_rank == 0) cout << " > N: " << n << " BL thickness delta: " << delta << " DT at r+delta: " << dtff << endl;


      part->solid_type = SPHERE_SOLID;
      part->ddata_solid_[0] = xc[0];
      part->ddata_solid_[1] = xc[1];
      part->ddata_solid_[2] = xc[2];
      part->ddata_solid_[3] = r;

      part->ff_type = SPHERE_FF;
      part->ddata_ff_[0] = xc[0];
      part->ddata_ff_[1] = xc[1];
      part->ddata_ff_[2] = xc[2];
      part->ddata_ff_[3] = r+delta;

      // -------------------------
      // surface...
      // -------------------------
      assert(part->surface == NULL);
      part->surface = new SurfaceShm();
      part->surface->makeSphere(part->name,xc,r,nseg,true); // true: fluid outside

      // -------------------------
      // ff_surface...
      // -------------------------
      assert(part->ff_surface == NULL);
      part->ff_surface = new SurfaceShm();
      part->ff_surface->makeSphere("FF",xc,r+delta,nseg,false); // false: fluid inside
      assert(part->ff_surface_dxost == NULL);
      CTI_Mmap(part->ff_surface_dxost,part->ff_surface->nst);
      if (mpi_rank_shared == 0) {
        for (int ist = 0; ist < part->ff_surface->nst; ++ist) {
          part->ff_surface_dxost[ist] = dtff;
        }
      }
      MPI_Barrier(mpi_comm_shared);

      // -------------------------
      // pts...
      // -------------------------
      assert(part->pts == NULL);
      part->pts = new Points();
      part->pts->np_global = part->surface->nsp*n;

      // count...
      part->pts->np = 0;
      int rr_rank = 0;
      for (int isp = 0; isp < part->surface->nsp; ++isp) {
        if (rr_rank == mpi_rank) {
          part->pts->np += n;
        }
        ++rr_rank;
        if (rr_rank == mpi_size)
          rr_rank = 0;
      }

      part->pts->xp = new double[part->pts->np][3];
      part->pts->delta = new double[part->pts->np];

      int ip = 0;
      rr_rank = 0;
      for (int isp = 0; isp < part->surface->nsp; ++isp) {
        if (rr_rank == mpi_rank) {
          // extrude points from this xsp in the xsp-xc direction...
          const double dx[3] = DIFF(part->surface->xsp[isp],xc);
          const double dx_mag = MAG(dx);
          for (int i = 0; i < n; ++i) {
            const double rp = r + dn*(pow(sf,0.5+double(i))-1.0)/(sf-1.0);
            FOR_J3 part->pts->xp[ip][j] = xc[j] + rp*dx[j]/dx_mag;
            part->pts->delta[ip] = 1.208*rp/double(nseg)*1.5; // TODO: add a correction factor here
            ++ip;
          }
        }
        ++rr_rank;
        if (rr_rank == mpi_size)
          rr_rank = 0;
      }

      // this is only valid if there is one group...
      //assert(ip == part->pts->np);

      // add suggested views for this part...
      // TODO

      // this a a sphere, with BL r+delta (width: 1.2*(r+delta)), centered on X (target 0,0,0), y-plane, y=0



    }

    return ierr;

  }

  int makePoints1D(Part * part,Param * param,int& iarg) {

    set<int> zone_indices;

    bool b_dx = false;
    double dx;

    while (iarg < param->size()) {
      string token = param->getString(iarg++);
      if ((token == "ZONE")||(token == "ZONES")) {
        // expect a comma-delimited set of integers...
        const string zonesCsv = param->getString(iarg++);
        vector<int> tmp;
        MiscUtils::splitCsv(tmp,zonesCsv);
        for (int ii = 0; ii < tmp.size(); ++ii) {
          const int izone = tmp[ii];
          if ((izone >= 0)&&(izone < zoneVec.size())) {
            zone_indices.insert(izone);
          }
          else {
            CWARN("Unrecognized zone index: " << izone << "; skipping");
          }
        }
      }
      else if ((token == "ZONE_NAME")||(token == "ZONE_NAMES")) {
        // expect a comma-delimited list of names...
        const string zoneNamesCsv = param->getString(iarg++);
        vector<string> zone_names;
        MiscUtils::splitCsv(zone_names,zoneNamesCsv);
        for (int ii = 0; ii < zone_names.size(); ++ii) {
          const int izone = getZoneIndex(zone_names[ii]); // get the unique zoneIndex in the partVec
          if (izone >= 0) {
            zone_indices.insert(izone);
          }
          else {
            CWARN("Unrecognized zone name: " << zone_names[ii] << "; skipping");
          }
        }
      }
      else if (token == "DX") {
        b_dx = true;
        dx = param->getDouble(iarg++);
      }
      else if (token == "NLAYERS") {
        //part->ff_nlayers = param->getInt(iarg++);
        part->ff_nlayersString = param->getString(iarg++);
      }
      else {
        CWARN("Unknown POINTS_1D parameter " << token << "; skipping");
      }
    }

    if (zone_indices.empty()) {
      COUT1("Error: makePoints1D missing ZONES <comma-delimited-indices> or ZONE_NAMES <comma-delimited-part:zone_names>");
      return -1;
    }

    assert(partVec[0]->surface);  // this is a HACK. we assume the zones come from one part
    partVec[0]->surface->flagDisjointGroups(zone_indices);
    
    const int ipart = 0; // HACK the points 1D are added to the first part for now

    part->b_link_surface = true;
    assert(part->ipart_link.empty());
    part->ipart_link.push_back(ipart);
    assert(part->ist_offset_link.empty());
    part->ist_offset_link.push_back(0); // HACK offset of the ist's in part ipart, which is the first part now, so zero
    
    SurfaceShm * surface = partVec[0]->surface; // HACK insist on partVec[0]'s surface...

    // we are going to send the loops to rank0 as a Gatherv, but just the
    // regular loops. The periodic loops will only be used locally to
    // build the points...
    int * send_buf_int = new int[surface->spogr_i[surface->ngr]];
    int send_count = 0;
    
    vector<double> xpVec;
    double my_average = 0.0;
    double my_minmax[2] = { HUGE_VAL, HUGE_VAL };

    for (int igr = 0; igr < surface->ngr; ++igr) {
      const int nloops = surface->nlogr[igr].first;
      const int nloops_periodic = surface->nlogr[igr].second;
      assert(nloops == 2);
      assert((nloops_periodic == 0)||(nloops_periodic == 2));
      // loops use a repeated point to mark the end...

      // loop1...
      int sog = surface->spogr_i[igr];
      int isp_start = surface->spogr_v[sog++];
      int isp_next = isp_start;
      send_buf_int[send_count++] = isp_next;
      double xc0[3] = { 0.0, 0.0, 0.0 };
      double length = 0.0;
      do {
        const int isp_prev = isp_next;
        isp_next = surface->spogr_v[sog++];
        send_buf_int[send_count++] = isp_next;
        // work with link isp_prev->isp_next...
        //if (mpi_rank == 0) cout << "igr: " << igr << " loop1: " << isp_prev << " " << isp_next << endl;
	const double dlength = DIST(surface->xsp[isp_prev],surface->xsp[isp_next]);
	length += dlength;
	FOR_I3 xc0[i] += dlength*(surface->xsp[isp_prev][i]+surface->xsp[isp_next][i]);
      } while (isp_next != isp_start);
      FOR_I3 xc0[i] /= length*2.0;
      
      // loop2...
      isp_start = surface->spogr_v[sog++];
      isp_next = isp_start;
      send_buf_int[send_count++] = isp_next;
      double xc1[3] = { 0.0, 0.0, 0.0 };
      length = 0.0;
      do {
        const int isp_prev = isp_next;
        isp_next = surface->spogr_v[sog++];
        send_buf_int[send_count++] = isp_next;
        // work with link isp_prev->isp_next...
        //if (mpi_rank == 0) cout << "igr: " << igr << " loop2: " << isp_prev << " " << isp_next << endl;
	const double dlength = DIST(surface->xsp[isp_prev],surface->xsp[isp_next]);
	length += dlength;
	FOR_I3 xc1[i] += dlength*(surface->xsp[isp_prev][i]+surface->xsp[isp_next][i]);
      } while (isp_next != isp_start);
      FOR_I3 xc1[i] /= length*2.0;

      if (nloops_periodic == 0) {
        
        // non-periodic treatment...
        const double dist = DIST(xc0,xc1);
        my_average += dist;
        my_minmax[0] = min(my_minmax[0],dist);
        my_minmax[1] = min(my_minmax[1],-dist);
	const int n = max(1,int(dist/dx)); // always atleast 1 point
	double ds[3] = DIFF(xc1,xc0);
	for (int i = 0; i < n; ++i) {
	  FOR_J3 xpVec.push_back(xc0[j] + ds[j]/double(n)/2.0*double(2*i+1));
	}
        
      }
      else {

        // periodic part...
        assert(nloops_periodic == 2);
        
        // periodic loop1...
        isp_start = surface->spogr_v[sog++];
        isp_next = isp_start;
	double xcp0[3] = { 0.0, 0.0, 0.0 };
	length = 0.0;
        do {
          const int isp_prev = isp_next;
          isp_next = surface->spogr_v[sog++];
	  const double dlength = DIST(surface->xsp[isp_prev],surface->xsp[isp_next]);
	  length += dlength;
	  FOR_I3 xcp0[i] += dlength*(surface->xsp[isp_prev][i]+surface->xsp[isp_next][i]);
        } while (isp_next != isp_start);
	FOR_I3 xcp0[i] /= length*2.0;
        
        // periodic loop2...
        isp_start = surface->spogr_v[sog++];
        isp_next = isp_start;
	double xcp1[3] = { 0.0, 0.0, 0.0 };
	length = 0.0;
        do {
          const int isp_prev = isp_next;
          isp_next = surface->spogr_v[sog++];
	  const double dlength = DIST(surface->xsp[isp_prev],surface->xsp[isp_next]);
	  length += dlength;
	  FOR_I3 xcp1[i] += dlength*(surface->xsp[isp_prev][i]+surface->xsp[isp_next][i]);
        } while (isp_next != isp_start);
        FOR_I3 xcp1[i] /= length*2.0;
	
        // use distance to ensure xcp0 is associated with xc0, and 
        // xcp1 is associated with xc1. There may be a better way to do this
        // in the future: for example, the loops are associated in the 
        // current order during creation...

        const double d1 = DIST(xc0,xcp0);
	const double d2 = DIST(xc0,xcp1);
	if ( d1 > d2 ) {
          assert( DIST(xc1,xcp1) > DIST(xc1,xcp0) );
          // swap xcp0 and xcp1...
          FOR_I3 {
            const double tmp = xcp0[i];
            xcp0[i] = xcp1[i];
            xcp1[i] = tmp;
          }
        }
        else {
          assert( DIST(xc1,xcp1) < DIST(xc1,xcp0) );
        }

        // now xcp0,xc0 and xcp1,xc1 are properly asssociated...
        
        const double dist0 = DIST(xc0,xcp0);
        const double dist1 = DIST(xc1,xcp1);
        const double dist = dist0+dist1; 
        my_average += dist;
        my_minmax[0] = min(my_minmax[0],dist);
        my_minmax[1] = min(my_minmax[1],-dist);
	const int n = max(1,int(dist/dx)); // always atleast 1 point
        const double dr = dist/double(n);
        // xc0 to xcp0 part
        const int n1 = int(floor(dist0*2.0/dr-1.0)/2) + 1;
        assert((n1 >= 0)&&(n1 <= n));
        double ds1[3] = DIFF(xcp0,xc0);
        double mag = MAG(ds1);
        assert(mag > 0.0);
        FOR_I3 ds1[i] /= mag;
        for (int i = 0; i < n1; ++i) {
          FOR_J3 xpVec.push_back(xc0[j] + 0.5*dr*ds1[j]*double(2*i+1));
        }
        // the other part
        double ds2[3] = DIFF(xcp1,xc1);
        mag = MAG(ds2);
        FOR_I3 ds2[i] /= mag;
        for (int i = 0; i < n-n1; ++i) {
          FOR_J3 xpVec.push_back(xc1[j] + 0.5*dr*ds2[j]*double(2*i+1));
        }
	
      }

      // we should be through the spogr_i...
      assert(sog == surface->spogr_i[igr+1]);

    }

    int np = xpVec.size()/3;
    int np_global;
    MPI_Allreduce(&np,&np_global,1,MPI_INT,MPI_SUM,mpi_comm);
    if (mpi_rank == 0) cout << " > np_global = " << np_global << endl;
    
    // pipe length stats...
    double average;
    MPI_Reduce(&my_average,&average,1,MPI_DOUBLE,MPI_SUM,0,mpi_comm);
    double minmax[2];
    MPI_Reduce(my_minmax,minmax,2,MPI_DOUBLE,MPI_MIN,0,mpi_comm);
    if (mpi_rank == 0) 
      cout << " > length stats [min,avg,max]: " << minmax[0] << " " << average/double(surface->ngr_global) << " " << -minmax[1] << endl;
    
    // -------------------------
    // pts...
    // -------------------------
    assert(part->pts == NULL);
    part->pts = new Points();
    part->pts->np_global = np_global;

    part->pts->np = np;
    part->pts->xp = new double[part->pts->np][3];
    part->pts->delta = new double[part->pts->np];

    for (int ip = 0; ip < np; ++ip) {
      FOR_I3 part->pts->xp[ip][i] = xpVec[ip*3+i];
      part->pts->delta[ip] = 1.2*dx; // TODO: what factor to use here?
    }
    xpVec.clear();

    // check that we didn't run out of memory here. Should be 
    // less than or equal to the surface->spogr_v size because
    // we are not sending periodic loops...
    
    assert(send_count <= surface->spogr_i[surface->ngr]);
    
    // gather the loop data at rank0...
    
    int * recv_count = NULL;
    if (mpi_rank == 0) recv_count = new int[mpi_size];
    MPI_Gather(&send_count,1,MPI_INT,recv_count,1,MPI_INT,0,mpi_comm);

    int * recv_disp = NULL;
    int recv_count_sum;
    int * recv_buf_int = NULL;
    if (mpi_rank == 0) {
      recv_disp = new int[mpi_size];
      recv_disp[0] = 0;
      for (int rank = 1; rank < mpi_size; ++rank)
        recv_disp[rank] = recv_disp[rank-1] + recv_count[rank-1];
      recv_count_sum = recv_disp[mpi_size-1] + recv_count[mpi_size-1];
      recv_buf_int = new int[recv_count_sum];
    }

    MPI_Gatherv(send_buf_int,send_count,MPI_INT,recv_buf_int,recv_count,recv_disp,MPI_INT,0,mpi_comm);
    delete[] send_buf_int;

    // recall that the SurfaceShm has a full/complete copy of the surface tris in shared memory on
    // every node. To avoid messaging the other nodes later, just communicate this recv_buf_int to 
    // everyone...
    
    int buf[2];
    if (mpi_rank_shared == 0) {
      
      // get recv_count_sum and recv_buf_int on all nodes. It is currently
      // only on the first...

      MPI_Bcast(&recv_count_sum,1,MPI_INT,0,mpi_comm_internode);
      if (mpi_rank != 0) {
        assert(recv_buf_int == NULL);
        recv_buf_int = new int[recv_count_sum];
      }
      MPI_Bcast(recv_buf_int,recv_count_sum,MPI_INT,0,mpi_comm_internode);
      
      // now the first rank on every node has the recv_buf_int...

      // the number of surface points should be exactly the recv_count_sum, 
      // because each loop requires a center point, however each loop is
      // identified by a duplicate node.
      
      buf[0] = recv_count_sum; // ff_nsp
      
      // the number of tris should be the number of points around the
      // loops. Since we know the global number of groups, this is just...
      
      buf[1] = recv_count_sum - 2*surface->ngr_global; // ff_nst (2 loops per group)
      
    }

    MPI_Bcast(buf,2,MPI_INT,0,mpi_comm_shared);
    
    // now everyone has the ff_surface tri counts and we can allocate the ff_surface... 

    const int ff_nsp = buf[0];
    const int ff_nst = buf[1];
    
    // -------------------------
    // ff_surface...
    // -------------------------
    assert(part->ff_surface == NULL);
    part->ff_surface = new SurfaceShm();
    part->ff_surface->init(ff_nsp,ff_nst);
    part->ff_surface->zoneVec.push_back(StZone("FF"));
    assert(part->ff_surface_dxost == NULL);
    CTI_Mmap(part->ff_surface_dxost,ff_nst);

    assert(part->ff_type == UNKNOWN_FF);
    part->ff_type = FF_SURFACE_FAZONE_FF;
    assert(part->solid_type == UNKNOWN_SOLID);
    part->solid_type = NO_SOLID;

    SurfaceShm * ff_surface = part->ff_surface;

    // build on shared rank 0 of every node...
    
    if (mpi_rank_shared == 0) {
      
      int nloop = 0;
      int irecv = 0;
      int ist = 0;
      while (irecv < recv_count_sum) {
        // pull the next loop...
        const int isp0 = irecv;
        const int isp_start = recv_buf_int[irecv++];
        FOR_I3 ff_surface->xsp[irecv][i] = surface->xsp[isp_start][i];
        int isp_next = isp_start;
        double xc[3] = { 0.0, 0.0, 0.0 };
        double length = 0.0;
        do {
          const int isp_prev = isp_next;
          isp_next = recv_buf_int[irecv++];
          const double dlength = DIST(surface->xsp[isp_prev],surface->xsp[isp_next]);
          length += dlength;
          FOR_I3 xc[i] += dlength*(surface->xsp[isp_prev][i]+surface->xsp[isp_next][i]);
          if (isp_next != isp_start) {
            FOR_I3 ff_surface->xsp[irecv][i] = surface->xsp[isp_next][i];
            ff_surface->spost[ist][0] = isp0;
            ff_surface->spost[ist][1] = irecv;
            ff_surface->spost[ist][2] = irecv-1;
	    part->ff_surface_dxost[ist] = dx;
          }
          else {
            ff_surface->spost[ist][0] = isp0;
            ff_surface->spost[ist][1] = isp0+1;
            ff_surface->spost[ist][2] = irecv-1;
	    part->ff_surface_dxost[ist] = dx;
          }
          ++ist;
        } while (isp_next != isp_start);
        FOR_I3 xc[i] /= length*2.0;
        FOR_I3 ff_surface->xsp[isp0][i] = xc[i];
        ++nloop;
      }
      assert(irecv == recv_count_sum);
      assert(nloop == 2*surface->ngr_global);
      assert(ist == ff_nst);
    }
    MPI_Barrier(mpi_comm_shared);

    /*
      if (mpi_rank == 0) {
      GeomUtils::writeTecplot("ff_surface.dat",part->ff_surface->spost,part->ff_surface->nst,part->ff_surface->xsp);
      GeomUtils::writeSbin("ff_surface.sbin",part->ff_surface->spost,part->ff_surface->nst,part->ff_surface->xsp);
      }
      GeomUtils::writePtsTecplot("pts.dat",part->pts->xp,part->pts->np);
    */

    // finally, set the surface flags. HACK for now use the set zone_indices
    // used to flagDisjointGroups above...
    
    if (mpi_rank_shared == 0) {
      
      const int nsz = surface->zoneVec.size();
      int * sz_flag = new int[nsz];
      for (int isz = 0; isz < nsz; ++isz)
        sz_flag[isz] = 0;
      
      for (set<int>::iterator it = zone_indices.begin(); it != zone_indices.end(); ++it) {
        const int isz = *it;
        assert((isz >= 0)&&(isz < nsz));
        sz_flag[isz] = 1;
      }
      
      for (int ist = 0; ist < surface->nst; ++ist) {
        const int isz = surface->znost[ist];
        assert((isz >= 0)&&(isz < nsz));
        if (sz_flag[isz]) {
          partVec[ipart]->setPartForSt(part->ipart_of_part,ist);
        }
      }
      
    }
    MPI_Barrier(mpi_comm_shared);
    
    return 0;

  }

  void writeTecplotForFlagedSt(string& fileName, int& ipart, int * st_flag, const int& izone, double* sp_buf) {
    SurfaceShm * surface = partVec[ipart]->surface;

    int num_st = 0;
    for (int ist = 0; ist < surface->nst; ist++) {
      if (st_flag[ist] == 1) {
        num_st++;
      }
    }

    FILE * fp = fopen(fileName.c_str(),"w");
    fprintf(fp,"TITLE = \"title\"\n");
    fprintf(fp,"VARIABLES = \"X\"\n");
    fprintf(fp,"\"Y\"\n");
    fprintf(fp,"\"Z\"\n");
    fprintf(fp,"\"VAL\"\n");
    fprintf(fp,"ZONE T=\"flag_zone_%d\"\n", izone);
    fprintf(fp,"N=%d, E=%d, F=FEPOINT, ET=TRIANGLE\n",num_st*3,num_st);


    for (int ist = 0; ist < surface->nst; ist++) {
      if (st_flag[ist] == 1) {
        fprintf(fp,"%lf %lf %lf %lf\n",surface->xsp[surface->spost[ist][0]][0],
                surface->xsp[surface->spost[ist][0]][1],
                surface->xsp[surface->spost[ist][0]][2],
                (sp_buf[surface->spost[ist][0]]+sp_buf[surface->spost[ist][1]]+sp_buf[surface->spost[ist][2]])/3.);
        fprintf(fp,"%lf %lf %lf %lf\n",surface->xsp[surface->spost[ist][1]][0],
                surface->xsp[surface->spost[ist][1]][1],
                surface->xsp[surface->spost[ist][1]][2],
                (sp_buf[surface->spost[ist][0]]+sp_buf[surface->spost[ist][1]]+sp_buf[surface->spost[ist][2]])/3.);
        fprintf(fp,"%lf %lf %lf %lf\n",surface->xsp[surface->spost[ist][2]][0],
                surface->xsp[surface->spost[ist][2]][1],
                surface->xsp[surface->spost[ist][2]][2],
                (sp_buf[surface->spost[ist][0]]+sp_buf[surface->spost[ist][1]]+sp_buf[surface->spost[ist][2]])/3.);
      }
    }

    int ii = 1;
    for (int ist = 0; ist < surface->nst; ist++) {
      if (st_flag[ist] == 1) {
        fprintf(fp,"%d %d %d\n", ii, ii+1, ii+2);
        ii+=3;
      }
    }

    fclose(fp);

  }

  int makeMedialAxis(Part * part,Param * param,int& iarg) {

    assert(0); // need code review

    for (int iz = 0, n = zoneVec.size(); iz < n; iz++) {
      COUT1(" > zone " << iz << " name: " << zoneVec[iz].getName());
    }
    set<int> zone_indices; // TODO screwed this up need to make set<pair<int,int> > ?
    vector<int> zone_npoints;

    while (iarg < param->size()) {
      string token = param->getString(iarg++);
      if ((token == "ZONE")||(token == "ZONES")) {
        // expect a comma-delimited set of integers...
        const string zonesCsv = param->getString(iarg++);
        vector<int> tmp;
        MiscUtils::splitCsv(tmp,zonesCsv);
        for (int ii = 0; ii < tmp.size(); ++ii) {
          const int izone = tmp[ii];
          if ((izone >= 0)&&(izone < zoneVec.size())) {
            zone_indices.insert(izone);
          }
          else {
            CWARN("Unrecognized zone index: " << izone << "; skipping");
          }
        }
      }
      else if ((token == "ZONE_NAME")||(token == "ZONE_NAMES")) {
        // expect a comma-delimited list of names...
        const string zoneNamesCsv = param->getString(iarg++);
        vector<string> zone_names;
        MiscUtils::splitCsv(zone_names,zoneNamesCsv);
        for (int ii = 0; ii < zone_names.size(); ++ii) {
          const int izone = getZoneIndex(zone_names[ii]); // get the unique zoneIndex in the partVec
          if (izone >= 0) {
            zone_indices.insert(izone);
          }
          else {
            CWARN("Unrecognized zone name: " << zone_names[ii] << "; skipping");
          }
        }
      }
      else if ((token == "NPOINTS")) {
        // expect a comma-delimited set of integers...
        const string ptsCsv = param->getString(iarg++);
        vector<int> tmp;
        MiscUtils::splitCsv(tmp,ptsCsv);
        for (int ii = 0; ii < tmp.size(); ++ii) {
          int pts = tmp[ii];
          if (pts > 0) {
            zone_npoints.push_back(pts);
          }
          else {
            CERR("Wrong number of points for position: " << ii);
          }
        }
      }
      else {
        CWARN("Unknown MEDIAL_AXIS parameter " << token << "; skipping");
      }
    }

    if (zone_indices.size() != zone_npoints.size())
      CERR("Number of zones doesn't match the number of points");

    //for (int i = 0; i<zone_indices.size(); i++) {
    //  cout << "zone index: " << zone_indices[i] << " ";
    //}
    //cout << endl;
    //for (int i = 0; i<zone_npoints.size(); i++) {
    //  cout << "zone npoints: " << zone_npoints[i] << " ";
    //}
    //cout << endl;
    //getchar();

    vector<GroupData> groupDataVec;
    buildGroupsForFlaggedZones(groupDataVec,zone_indices);

    // keep track of Geo distance calc for each part
    vector<bool> b_geoDistVec(partVec.size(), false);
    vector<double*> sp_bufVec(partVec.size(), NULL);

    // cubic splines for each zone
    SplineStuff::CubicSpline * cspline = new SplineStuff::CubicSpline[groupDataVec.size()];
    SplineStuff::CubicSplineXy * cspline_sd = new SplineStuff::CubicSplineXy[groupDataVec.size()];

    int ff_nsp = 0;
    int ff_nst = 0;
    int np_global = 0;
    int np = 0;
    int rr_rank = 0;
    for (int izone = 0; izone < zone_npoints.size(); izone++) {
      np_global += zone_npoints[izone];
      for (int ip = 0; ip < zone_npoints[izone]; ip++) {
        if (rr_rank == mpi_rank) {
          np++;
        }
        rr_rank = (rr_rank+1) % mpi_size;
      }
    }

    for (int igr = 0; igr < groupDataVec.size(); igr++) {
      COUT1(" > processing group " << igr);
      SurfaceShm * surface = partVec[groupDataVec[igr].ipart]->surface;

      if (mpi_rank == 0) {
        cout << " >> group " << igr << " link 0 normal: " << COUT_VEC(groupDataVec[igr].n_li[0]) << " , link 0 x: " << COUT_VEC(groupDataVec[igr].x_li[0]) << endl;
        cout << " >> group " << igr << " link 1 normal: " << COUT_VEC(groupDataVec[igr].n_li[1]) << " , link 1 x: " << COUT_VEC(groupDataVec[igr].x_li[1]) << endl;
      }

      // get the link information for ff_surface
      ff_nsp += groupDataVec[igr].spoli_s + 2;
      ff_nst += groupDataVec[igr].spoli_s;

      // set distance function in sp_buf
      // NOTE: This is the most sensitive part of the algorithm,
      // right now the goemetric distance is computed from a point.
      // you may get different solutions with different geo distance functions.
      // The start point can be important, link 0/1.
      double * sp_buf = NULL;
      if (not b_geoDistVec[groupDataVec[igr].ipart]) {
        assert(sp_bufVec[groupDataVec[igr].ipart] == NULL);
        vector<int> ispVec;
        int ilink = 1;
        for (int sol = groupDataVec[igr].spoli_i[ilink]; sol != groupDataVec[igr].spoli_i[ilink+1]; sol++) {
          int isp = groupDataVec[igr].spoli_v[sol];
          ispVec.push_back(isp);
        }
        //double x_start[3];
        //FOR_I3 x_start[i] = groupDataVec[igr].x_li[1][i];
        sp_buf = new double [surface->nsp];
        //calcGeoDistFromPoint(x_start, groupDataVec[igr].ipart, sp_buf);
        calcGeoDistFromPoint(ispVec, groupDataVec[igr].ipart, sp_buf);
        sp_bufVec[groupDataVec[igr].ipart] = sp_buf;
        b_geoDistVec[groupDataVec[igr].ipart] = true;
      }
      else {
        sp_buf = sp_bufVec[groupDataVec[igr].ipart];
        COUT1(" > using existing Geo distance function for part " << groupDataVec[igr].ipart);
      }

      int * st_flag = new int [surface->nst];
      for (int ist = 0; ist < surface->nst; ist++) {
        st_flag[ist] = -1;
      }

      FILE * fpts;
      if (mpi_rank == 0) {
        stringstream tmp;
        tmp << "points_gr_" << igr << ".dat";
        string pts_name = tmp.str();
        fpts = fopen(pts_name.c_str(), "w");
        fprintf(fpts,"VARIABLES = \"X\" \"Y\" \"Z\" \"NX\" \"NY\" \"NZ\"\n");
      }

      // start from link 0
      double x_cur[3]; FOR_I3 x_cur[i] = groupDataVec[igr].x_li[0][i];
      double n_cur[3]; FOR_I3 n_cur[i] = groupDataVec[igr].n_li[0][i];
      double mag_n = MAG(n_cur);
      FOR_I3 n_cur[i] /= mag_n;
      double dx_cur;
      bool done = false;
      int iteration = 0;
      vector<double*> xVec;
      vector<double> surfDistVec;
      while (not done) {
        // find the closest surface point to x_cur
        COUT1(" > iter: " << iteration);

        int isp_closest = -1;
        double dist_closest = HUGE_VAL;
        for (int ist = 0; ist < surface->nst; ist++) {
          int this_igr;
          if (partVec[groupDataVec[igr].ipart]->getGroupForSt(this_igr,ist)) {
            if (this_igr == igr) {
              FOR_I3 {
                const int isp = surface->spost[ist][i];
                const double this_dist = DIST(surface->xsp[isp],x_cur);
                if (this_dist < dist_closest) {
                  isp_closest = isp;
                  dist_closest = this_dist;
                }
              }
            }
          }
        }
        assert(isp_closest >= 0);

        dx_cur = dist_closest;

        double g_min = sp_buf[isp_closest] - dx_cur;
        double g_max = sp_buf[isp_closest] + dx_cur;

        int n_flag = 0;
        for (int ist = 0; ist < surface->nst; ist++) {
          int this_igr;
          if (partVec[groupDataVec[igr].ipart]->getGroupForSt(this_igr,ist)) {
            if (this_igr == igr) {
              double this_g = (sp_buf[surface->spost[ist][0]]+sp_buf[surface->spost[ist][1]]+sp_buf[surface->spost[ist][2]])/3.;
              if ( (this_g >= g_min) && (this_g < g_max)) {
                st_flag[ist] = 1;
                n_flag++;
              }
            }
          }
        }
        if (n_flag == 0) {
          CERR(" > Did not find any triangles for this iteration, exit...");
        }
        //stringstream temp; temp << "tris_for_gr_" << igr << "_iter_" << iteration << ".dat";
        //string fileName = temp.str();
        //writeTecplotForFlagedSt(fileName, groupDataVec[igr].ipart, st_flag, iteration, sp_buf);
        double xc[3] = {0,0,0};
        double nc[3] = {0,0,0};
        calcClosestPassageForFlaggedTris(xc,nc,groupDataVec[igr].ipart,st_flag,n_cur);
        if (mpi_rank == 0) fprintf(fpts, "%lf %lf %lf %lf %lf %lf\n",xc[0],xc[1],xc[2],nc[0],nc[1],nc[2]);
        xVec.push_back(new double[3]);
        FOR_I3 xVec[xVec.size()-1][i] = xc[i];
        double dist_farthest = -HUGE_VAL;
        for (int ist = 0; ist < surface->nst; ist++) {
          if (st_flag[ist] == 1) {
            const double * const x0 = surface->xsp[surface->spost[ist][0]];
            const double * const x1 = surface->xsp[surface->spost[ist][1]];
            const double * const x2 = surface->xsp[surface->spost[ist][2]];
            double x_st[3]; FOR_I3 x_st[i] = (x0[i]+x1[i]+x2[i])/3.0;
            double this_dist = DIST(xc, x_st);
            if (this_dist > dist_farthest) {
              dist_farthest = this_dist;
            }
          }
        }
        surfDistVec.push_back(dist_farthest);
        // update x and n
        FOR_I3 x_cur[i] = xc[i] + n_cur[i]*2.*dx_cur;
        FOR_I3 n_cur[i] = nc[i];
        // check crossing with link 1
        double dx0[3] = DIFF(xc, groupDataVec[igr].x_li[1]);
        double dx1[3] = DIFF(x_cur, groupDataVec[igr].x_li[1]);
        double dot_prod = DOT_PRODUCT(dx0, dx1);
        if (dot_prod < 0.0) {
          COUT1(" > crossed the other side, exit...");
          done = true;
        }
        // set the flag back
        for (int ist = 0; ist < surface->nst; ist++) {
          st_flag[ist] = -1;
        }
        iteration++;
      }
      if (mpi_rank == 0) fclose(fpts);

      for (int i = 0; i < xVec.size(); i++) {
        COUT1(" > xVec["<<i<<"]: "<<xVec[i][0]<<" "<<xVec[i][1]<<" "<<xVec[i][2]);
      }

      assert(surfDistVec.size() == xVec.size());

      // build cubic spline for this group
      double (*x_cspl)[3] = new double[xVec.size()+2][3];
      FOR_I3 x_cspl[0][i] = groupDataVec[igr].x_li[0][i];
      for (int ip = 1; ip < xVec.size()+1; ip++) {
        FOR_I3 x_cspl[ip][i] = xVec[ip-1][i];
      }
      FOR_I3 x_cspl[xVec.size()+1][i] = groupDataVec[igr].x_li[1][i];

      cspline[igr].init(x_cspl,xVec.size()+2);
      delete [] x_cspl;

      // build cubic spline for surface distance
      double * surfDist = new double [surfDistVec.size()];
      double * sp       = new double [surfDistVec.size()];
      for (int ip = 0; ip < surfDistVec.size(); ip++) {
        surfDist[ip] = surfDistVec[ip];
        sp[ip] = cspline[igr].getSp(ip+1);
      }
      cspline_sd[igr].init(sp, surfDist, surfDistVec.size());

      //if (mpi_rank == 0) {
      //  stringstream tem3;
      //  tem3 << "Nodes_surfDist_gr_" << igr << ".dat";
      //  string fn3 = tem3.str();
      //  FILE * spfp3 = fopen(fn3.c_str(),"w");
      //  for (int ii = 0; ii < surfDistVec.size(); ii++) {
      //    fprintf(spfp3,"%18.15le %18.15le\n",sp[ii] ,surfDist[ii]);
      //  }
      //  fclose(spfp3);
      //}

      delete [] surfDist;
      delete [] sp;

      //if (mpi_rank == 0) {
      //  stringstream tem;
      //  tem << "cspline_surfDist_gr_" << igr << ".dat";
      //  string fn2 = tem.str();
      //  FILE * spfp2 = fopen(fn2.c_str(),"w");
      //  for (int ii = 0; ii <= zone_npoints[igr]; ii++) {
      //    const double s = double(ii)/double(zone_npoints[igr]);
      //    double dis = cspline_sd[igr].getY(s);
      //    fprintf(spfp2,"%18.15le %18.15le\n",s ,dis);
      //  }
      //  fclose(spfp2);
      //}

      if (mpi_rank == 0) {
        stringstream temp;
        temp << "cspline_gr_" << igr << ".dat";
        string fn = temp.str();
        FILE * spfp = fopen(fn.c_str(),"w");
        fprintf(spfp,"TITLE = \"title\"\n");
        fprintf(spfp,"VARIABLES = \"X\"\n");
        fprintf(spfp,"\"Y\"\n");
        fprintf(spfp,"\"Z\"\n");
        fprintf(spfp,"ZONE T=\"spline_zone\"\n");
        for (int ii = 0; ii <= zone_npoints[igr]; ii++) {
          const double s = double(ii)/double(zone_npoints[igr]);
          double x[3];
          cspline[igr].getX(x,s);
          fprintf(spfp,"%18.15le %18.15le %18.15le\n",x[0],x[1],x[2]);
        }
        fclose(spfp);
      }

      int num_st = 0;
      for (int ist = 0; ist < surface->nst; ist++) {
        int this_igr;
        if (partVec[groupDataVec[igr].ipart]->getGroupForSt(this_igr,ist)) {
          if (this_igr == igr) {
            num_st++;
          }
        }
      }

      if (mpi_rank == 0) {
        stringstream temp2;
        temp2 << "subSurf_gr_" << igr << ".dat";
        string fileName = temp2.str();
        FILE * fp2 = fopen(fileName.c_str(),"w");
        fprintf(fp2,"TITLE = \"title\"\n");
        fprintf(fp2,"VARIABLES = \"X\"\n");
        fprintf(fp2,"\"Y\"\n");
        fprintf(fp2,"\"Z\"\n");
        fprintf(fp2,"\"VAL\"\n");
        fprintf(fp2,"ZONE T=\"sub_zone\"\n");
        fprintf(fp2,"N=%d, E=%d, F=FEPOINT, ET=TRIANGLE\n",num_st*3,num_st);


        for (int ist = 0; ist < surface->nst; ist++) {
          int this_igr;
          if (partVec[groupDataVec[igr].ipart]->getGroupForSt(this_igr,ist)) {
            if (this_igr == igr) {
              fprintf(fp2,"%lf %lf %lf %lf\n",surface->xsp[surface->spost[ist][0]][0],
                      surface->xsp[surface->spost[ist][0]][1],
                      surface->xsp[surface->spost[ist][0]][2],
                      (sp_buf[surface->spost[ist][0]]+sp_buf[surface->spost[ist][1]]+sp_buf[surface->spost[ist][2]])/3.);
              fprintf(fp2,"%lf %lf %lf %lf\n",surface->xsp[surface->spost[ist][1]][0],
                      surface->xsp[surface->spost[ist][1]][1],
                      surface->xsp[surface->spost[ist][1]][2],
                      (sp_buf[surface->spost[ist][0]]+sp_buf[surface->spost[ist][1]]+sp_buf[surface->spost[ist][2]])/3.);
              fprintf(fp2,"%lf %lf %lf %lf\n",surface->xsp[surface->spost[ist][2]][0],
                      surface->xsp[surface->spost[ist][2]][1],
                      surface->xsp[surface->spost[ist][2]][2],
                      (sp_buf[surface->spost[ist][0]]+sp_buf[surface->spost[ist][1]]+sp_buf[surface->spost[ist][2]])/3.);
            }
          }
        }

        int ii = 1;
        for (int ist = 0; ist < surface->nst; ist++) {
          int this_igr;
          if (partVec[groupDataVec[igr].ipart]->getGroupForSt(this_igr,ist)) {
            if (this_igr == igr) {
              fprintf(fp2,"%d %d %d\n", ii, ii+1, ii+2);
              ii+=3;
            }
          }
        }

        fclose(fp2);
      }

      delete [] st_flag;
      for (int i = 0; i < xVec.size(); i++) {
        delete [] xVec[i];
      }

    } // group loop


      /*
        SurfaceShm * surface = partVec[0]->surface;
        double * sp_buf = sp_bufVec[0];

        FILE * fp3 = fopen("TotSurf_total.dat","w");
        fprintf(fp3,"TITLE = \"title\"\n");
        fprintf(fp3,"VARIABLES = \"X\"\n");
        fprintf(fp3,"\"Y\"\n");
        fprintf(fp3,"\"Z\"\n");
        fprintf(fp3,"\"VAL\"\n");
        fprintf(fp3,"ZONE T=\"total_zone\"\n");
        fprintf(fp3,"N=%d, E=%d, F=FEPOINT, ET=TRIANGLE\n",surface->nst*3,surface->nst);


        for (int ist = 0; ist < surface->nst; ist++) {
        fprintf(fp3,"%lf %lf %lf %lf\n",surface->xsp[surface->spost[ist][0]][0],
        surface->xsp[surface->spost[ist][0]][1],
        surface->xsp[surface->spost[ist][0]][2],
        (sp_buf[surface->spost[ist][0]]+sp_buf[surface->spost[ist][1]]+sp_buf[surface->spost[ist][2]])/3.);
        fprintf(fp3,"%lf %lf %lf %lf\n",surface->xsp[surface->spost[ist][1]][0],
        surface->xsp[surface->spost[ist][1]][1],
        surface->xsp[surface->spost[ist][1]][2],
        (sp_buf[surface->spost[ist][0]]+sp_buf[surface->spost[ist][1]]+sp_buf[surface->spost[ist][2]])/3.);
        fprintf(fp3,"%lf %lf %lf %lf\n",surface->xsp[surface->spost[ist][2]][0],
        surface->xsp[surface->spost[ist][2]][1],
        surface->xsp[surface->spost[ist][2]][2],
        (sp_buf[surface->spost[ist][0]]+sp_buf[surface->spost[ist][1]]+sp_buf[surface->spost[ist][2]])/3.);
        }

        int ii = 1;
        for (int ist = 0; ist < surface->nst; ist++) {
        fprintf(fp3,"%d %d %d\n", ii, ii+1, ii+2);
        ii+=3;
        }

        fclose(fp3);
      */
    for (int i = 0; i < sp_bufVec.size(); i++) {
      delete [] sp_bufVec[i];
    }


    // -------------------------
    // surface should not exist
    // -------------------------
    assert(part->surface == NULL);
    // -------------------------
    // ff_surface...
    // -------------------------
    assert(part->ff_surface == NULL);
    part->ff_surface = new SurfaceShm();
    part->ff_surface->init(ff_nsp,ff_nst);
    part->ff_surface->zoneVec.push_back(StZone("FF"));
    assert(part->ff_surface_dxost == NULL);
    CTI_Mmap(part->ff_surface_dxost,part->ff_surface->nst);

    COUT1(" > ff_nsp: " << ff_nsp << " ff_nst: " << ff_nst << " np_global: " << np_global);

    // -------------------------
    // pts...
    // -------------------------
    assert(part->pts == NULL);
    part->pts = new Points();
    part->pts->np_global = np_global;

    part->pts->np = np;
    part->pts->xp = new double[part->pts->np][3];
    part->pts->delta = new double[part->pts->np];

    ff_nsp = 0;
    ff_nst = 0;
    np_global = 0;
    np = 0;
    rr_rank = 0;
    for (int igr = 0; igr < groupDataVec.size(); ++igr) {
      SurfaceShm * surface = partVec[groupDataVec[igr].ipart]->surface; assert(surface);
      assert(groupDataVec[igr].nlink == 2);
      if (mpi_rank_shared == 0) {
        // link 0 first...
        // ff_surface points...
        int ff_isp0 = ff_nsp++;
        FOR_I3 part->ff_surface->xsp[ff_isp0][i] = groupDataVec[igr].x_li[0][i];
        for (int sol = groupDataVec[igr].spoli_i[0]; sol != groupDataVec[igr].spoli_i[1]; ++sol) {
          const int isp = groupDataVec[igr].spoli_v[sol];
          FOR_I3 part->ff_surface->xsp[ff_nsp][i] = surface->xsp[isp][i];
          ++ff_nsp;
        }
        // ff_surface tris...
        for (int ff_isp1 = ff_isp0+1; ff_isp1 < ff_nsp; ++ff_isp1) {
          int ff_isp2 = ff_isp1+1;
          if (ff_isp2 == ff_nsp)
            ff_isp2 = ff_isp0+1;
          part->ff_surface->spost[ff_nst][0] = ff_isp0;
          part->ff_surface->spost[ff_nst][1] = ff_isp2;
          part->ff_surface->spost[ff_nst][2] = ff_isp1;
          part->ff_surface->znost[ff_nst] = 0;
          part->ff_surface_dxost[ff_nst] = groupDataVec[igr].d_li[0]; // TODO find the correct value
          ++ff_nst;
        }
        // then link 1...
        // ff_surface points...
        ff_isp0 = ff_nsp++;
        FOR_I3 part->ff_surface->xsp[ff_isp0][i] = groupDataVec[igr].x_li[1][i];
        for (int sol = groupDataVec[igr].spoli_i[1]; sol != groupDataVec[igr].spoli_i[2]; ++sol) {
          const int isp = groupDataVec[igr].spoli_v[sol];
          FOR_I3 part->ff_surface->xsp[ff_nsp][i] = surface->xsp[isp][i];
          ++ff_nsp;
        }
        // ff_surface tris...
        for (int ff_isp1 = ff_isp0+1; ff_isp1 < ff_nsp; ++ff_isp1) {
          int ff_isp2 = ff_isp1+1;
          if (ff_isp2 == ff_nsp)
            ff_isp2 = ff_isp0+1;
          part->ff_surface->spost[ff_nst][0] = ff_isp0;
          part->ff_surface->spost[ff_nst][1] = ff_isp2;
          part->ff_surface->spost[ff_nst][2] = ff_isp1;
          part->ff_surface->znost[ff_nst] = 0;
          part->ff_surface_dxost[ff_nst] = groupDataVec[igr].d_li[1]; // TODO find the correct value
          ++ff_nst;
        }
      }

      // points...
      const double dx_pts = cspline[igr].getLength()/double(zone_npoints[igr]);
      np_global += zone_npoints[igr]; // add atleast 1 point...
      // distribute the points...
      for (int ip = 0; ip < zone_npoints[igr]; ++ip) {
        if (rr_rank == mpi_rank) {
          const double s = (double(ip)+0.5)/double(zone_npoints[igr]);
          double x[3];
          cspline[igr].getX(x,s);
          FOR_I3 part->pts->xp[np][i] = x[i];
          double surfDist = cspline_sd[igr].getY(s);
          part->pts->delta[np] = 1.1*max(dx_pts,surfDist);
          ++np;
        }
        rr_rank = (rr_rank+1) % mpi_size;
      }
    }

    delete[] cspline;
    delete[] cspline_sd;

    if (mpi_rank_shared == 0) {
      assert(ff_nsp == part->ff_surface->nsp);
      assert(ff_nst == part->ff_surface->nst);
    }

    assert(np_global == part->pts->np_global);
    assert(np == part->pts->np);

    for (int igr = 0; igr < groupDataVec.size(); ++igr)
      groupDataVec[igr].clear();

    return 0;

  }

  void calcClosestPassageForFlaggedTris(double xc[3],double nc[3],const int& ipart,const int* st_flag, const double dir[3]) {

    const SurfaceShm * const surface = partVec[ipart]->surface;
    double LS_mat[3][3] = { { 0.0, 0.0, 0.0 },
                            { 0.0, 0.0, 0.0 },
                            { 0.0, 0.0, 0.0 } };

    for (int ist = 0; ist < partVec[ipart]->surface->nst; ist++) {
      if (st_flag[ist] == 1) {
        double unit_n_st[3];
        const double * const x0 = surface->xsp[surface->spost[ist][0]];
        const double * const x1 = surface->xsp[surface->spost[ist][1]];
        const double * const x2 = surface->xsp[surface->spost[ist][2]];
        const double this_n[3] = TRI_NORMAL_2(x0,x1,x2);
        const double mag = MAG(this_n);
        assert(mag > 0.0);
        FOR_I3 unit_n_st[i] = this_n[i]/mag;

        for (int i=0; i<3; i++) {
          for (int j=0; j<3; j++) {
            LS_mat[i][j] += unit_n_st[i] * unit_n_st[j];
          }
        }

      }
    }

    double eig_vec[3][3];
    double eig_val[3];
    MiscUtils::eigenDecomposition(LS_mat, eig_vec, eig_val);
    //FOR_I3 cout << "eig_val i: " << i << " " << eig_val[i] << endl;
    //FOR_I3
    //  FOR_J3
    //  cout << "eig_vec[" << i << "][" << j << "]: "<< eig_vec[i][j] << endl;

    // cylinder normal
    FOR_I3 nc[i] = eig_vec[i][0];
    if (DOT_PRODUCT(nc,dir) < 0.0) {
      FOR_I3 nc[i] *= -1;
    }

    // cylinder other axes
    double uc[3]; FOR_I3 uc[i] = eig_vec[i][1];
    double vc[3] = CROSS_PRODUCT(nc,uc);
    assert(MAG(nc) - 1.0 < 1e-6);
    assert(MAG(uc) - 1.0 < 1e-6);
    assert(MAG(vc) - 1.0 < 1e-6);
    assert(DOT_PRODUCT(nc,uc) < 1e-6);
    assert(DOT_PRODUCT(uc,vc) < 1e-6);
    assert(DOT_PRODUCT(vc,nc) < 1e-6);
    COUT1("   > nc: " << COUT_VEC(nc) << " ,uc: " << COUT_VEC(uc) << " ,vc: " << COUT_VEC(vc));

    FOR_I3 xc[3] = 0.0;
    double sum_wgt = 0.0;

    for (int ist = 0; ist < surface->nst; ist++) {
      if (st_flag[ist] == 1){
        const double * const x0 = surface->xsp[surface->spost[ist][0]];
        const double * const x1 = surface->xsp[surface->spost[ist][1]];
        const double * const x2 = surface->xsp[surface->spost[ist][2]];
        const double this_n[3] = TRI_NORMAL_2(x0,x1,x2);
        const double area = MAG(this_n);
        assert(area > 0.0);
        FOR_I3 xc[i] += area*(x0[i] + x1[i] + x2[i]);
        sum_wgt += area;
      }
    }
    FOR_I3 xc[i] /= sum_wgt*3.0;

    // for the radius, just take the mean radius fromthe guess...
    double r_cur = 0.0;
    sum_wgt = 0.0;
    for (int ist = 0; ist < surface->nst; ist++) {
      if (st_flag[ist] == 1){
        const double * const x0 = surface->xsp[surface->spost[ist][0]];
        const double * const x1 = surface->xsp[surface->spost[ist][1]];
        const double * const x2 = surface->xsp[surface->spost[ist][2]];
        const double this_n[3] = TRI_NORMAL_2(x0,x1,x2);
        const double area = MAG(this_n);
        assert(area > 0.0);
        double x0_2d[2] = {DOT_PRODUCT(x0,uc), DOT_PRODUCT(x0,vc)};
        r_cur += area*DIST_2D(xc,x0_2d);
        double x1_2d[2] = {DOT_PRODUCT(x1,uc), DOT_PRODUCT(x1,vc)};
        r_cur += area*DIST_2D(xc,x1_2d);
        double x2_2d[2] = {DOT_PRODUCT(x2,uc), DOT_PRODUCT(x2,vc)};
        r_cur += area*DIST_2D(xc,x2_2d);
        sum_wgt += area;
      }
    }
    r_cur /= sum_wgt*3.0;
    assert(r_cur > 0.0);

    COUT1("   > mean location: xc: " << COUT_VEC(xc) << " , mean radius r: " << r_cur);

  }

  // solve Eikonal equation starting from a provided point
  void calcGeoDistFromPoint(const double x[3], const int& ipart, double* sp_buf) {
    COUT1(" > calcGeoDistFromPoint for ipart " << ipart);

    SurfaceShm * surface = partVec[ipart]->surface;

    // build surface-tri-of-surface-point to allow use to go from isp to ied

    int * stosp_i = new int[surface->nsp+1]; // node,edge-of-surface-point
    for (int isp = 0; isp < surface->nsp; ++isp) stosp_i[isp+1] = 0;

    // count number of tris for each node
    for (int ist = 0; ist < surface->nst; ++ist) {
      FOR_I3 {
        stosp_i[surface->spost[ist][i]+1]++;
      }
    }

    // scan to get disp from counts...
    stosp_i[0] = 0;
    for (int isp = 0; isp < surface->nsp; ++isp) stosp_i[isp+1] += stosp_i[isp];
    assert(stosp_i[surface->nsp] == surface->nst*3);

    // populate stosp_v....
    int *stosp_v = new int[stosp_i[surface->nsp]];
    for (int ist = 0; ist < surface->nst; ++ist) {
      FOR_I3 {
        stosp_v[stosp_i[surface->spost[ist][i]]++] = ist;
      }
    }

    // rewind...
    for (int isp = surface->nsp-1; isp > 0; --isp) stosp_i[isp] = stosp_i[isp-1];
    stosp_i[0] = 0;

    // find closest point and put its g/isp into CloseVec
    int isp_closest = -1;
    double dist_closest = HUGE_VAL;
    for (int isp = 0; isp < surface->nsp; ++isp) {
      const double this_dist = DIST(surface->xsp[isp],x);
      if (this_dist < dist_closest) {
        isp_closest = isp;
        dist_closest = this_dist;
      }
    }
    assert(isp_closest >= 0);

    IntFlag sp_flag;

    // store g values...
    if (sp_buf == NULL) sp_buf = new double[surface->nsp];
    sp_flag.setLength(surface->nsp);
    for (int isp = 0; isp < surface->nsp; ++isp) {
      sp_buf[isp] = HUGE_VAL; // g
      sp_flag[isp] = 0; // 1: alive, 0: far, < 0: close
    }

    // following notation from Sethian Fast Marching Method, 1998...
    MinHeap trialHeap(surface->nsp); // nbrs to alive pts
    sp_buf[isp_closest] = dist_closest;
    sp_flag[isp_closest] = 1; // first heap element (-1 indexing)
    trialHeap.insert(pair<double,int>(dist_closest,isp_closest)); // negated because heap is descending

    // perform fmm loops until finished...
    while (trialHeap.size() > 0) {
      //cout << trialHeap.size() << endl;
      fmmLoop(stosp_i,stosp_v,trialHeap,sp_flag,sp_buf,surface);
    }

    // negate g values (or add invert to color map)?
    //for (int isp = 0; isp < surface->nsp; ++isp) sp_buf[isp] = -sp_buf[isp];

    delete[] stosp_i;
    delete[] stosp_v;
    sp_flag.clear();

  }

  // solve Eikonal equation starting from a provided point
  void calcGeoDistFromPoint(vector<int> ispVec, const int& ipart, double* sp_buf) {
    COUT1(" > calcGeoDistFromPoint for ipart " << ipart);

    SurfaceShm * surface = partVec[ipart]->surface;

    // build surface-tri-of-surface-point to allow use to go from isp to ied

    int * stosp_i = new int[surface->nsp+1]; // node,edge-of-surface-point
    for (int isp = 0; isp < surface->nsp; ++isp) stosp_i[isp+1] = 0;

    // count number of tris for each node
    for (int ist = 0; ist < surface->nst; ++ist) {
      FOR_I3 {
        stosp_i[surface->spost[ist][i]+1]++;
      }
    }

    // scan to get disp from counts...
    stosp_i[0] = 0;
    for (int isp = 0; isp < surface->nsp; ++isp) stosp_i[isp+1] += stosp_i[isp];
    assert(stosp_i[surface->nsp] == surface->nst*3);

    // populate stosp_v....
    int *stosp_v = new int[stosp_i[surface->nsp]];
    for (int ist = 0; ist < surface->nst; ++ist) {
      FOR_I3 {
        stosp_v[stosp_i[surface->spost[ist][i]]++] = ist;
      }
    }

    // rewind...
    for (int isp = surface->nsp-1; isp > 0; --isp) stosp_i[isp] = stosp_i[isp-1];
    stosp_i[0] = 0;

    IntFlag sp_flag;

    // store g values...
    if (sp_buf == NULL) sp_buf = new double[surface->nsp];
    sp_flag.setLength(surface->nsp);
    for (int isp = 0; isp < surface->nsp; ++isp) {
      sp_buf[isp] = HUGE_VAL; // g
      sp_flag[isp] = 0; // 1: alive, 0: far, < 0: close
    }

    // following notation from Sethian Fast Marching Method, 1998...
    MinHeap trialHeap(surface->nsp); // nbrs to alive pts

    for (int ip = 0; ip < ispVec.size(); ip++) {
      int isp = ispVec[ip];
      sp_buf[isp] = 0;
      sp_flag[isp] = 1; // first heap element (-1 indexing)
      trialHeap.insert(pair<double,int>(0,isp)); // negated because heap is descending
    }

    // perform fmm loops until finished...
    while (trialHeap.size() > 0) {
      //cout << trialHeap.size() << endl;
      fmmLoop(stosp_i,stosp_v,trialHeap,sp_flag,sp_buf,surface);
    }

    // negate g values (or add invert to color map)?
    //for (int isp = 0; isp < surface->nsp; ++isp) sp_buf[isp] = -sp_buf[isp];

    delete[] stosp_i;
    delete[] stosp_v;
    sp_flag.clear();

  }

  void fmmLoop(MinHeap& trialHeap,int* flag_gp,double* dist_gp,const int* stogp_i,const int* stogp_v,const int* gposp,SurfaceShm * surface) {

    // get trial point from close vector...
    pair<double,int> trial = trialHeap.extractMin();
    const int trial_igp = trial.second;
    trialHeap.deleteMin(); // remove from trial
    assert(flag_gp[trial_igp] == 1);
    flag_gp[trial_igp] = 2; // tag as alive

    // mark trial_igp neighbors that aren't alive as close...
    for (int tog = stogp_i[trial_igp]; tog != stogp_i[trial_igp+1]; ++tog) {
      const int ist = stogp_v[tog];
      FOR_I3 {
        const int igp = gposp[surface->spost[ist][i]];
        // some tri's touch non grouped points...
        if (igp != -1) {
          if (flag_gp[igp] == 0) {
            flag_gp[igp] = 1; // tag as close
            trialHeap.insert(pair<double,int>(HUGE_VAL,igp)); // push onto heap
          }
        }
      }
    }

    // recompute g's for all close trial_igp neighbors...
    for (int tog = stogp_i[trial_igp]; tog != stogp_i[trial_igp+1]; ++tog) {
      const int ist = stogp_v[tog];
      // count number of alive points
      int count = 0;
      FOR_I3 {
        const int igp = gposp[surface->spost[ist][i]];
        if (igp != -1) {
          if (flag_gp[igp] == 2) ++count;
        }
      }
      if (count == 3) {
        continue; // skip tris where all the pts are alive
      }
      else if (count == 2) {
        // find two alive points
        int isp_a,isp_b;
        double g_a = HUGE_VAL, g_b = HUGE_VAL; // a smallest alive value
        FOR_I3 {
          const int isp = surface->spost[ist][i];
          const int igp = gposp[isp];
          if (igp != -1) {
            const double this_g = dist_gp[igp];
            if (flag_gp[igp] == 2) {
              if (this_g < g_a) {
                isp_b = isp_a;
                g_b = g_a;
                isp_a = isp;
                g_a = this_g;
              }
              else if (this_g < g_b) {
                isp_b = isp;
                g_b = this_g;
              }
            }
          }
        }
        // update close points g
        FOR_I3 {
          const int isp = surface->spost[ist][i];
          const int igp = gposp[isp];
          if (igp != -1) {
            if (flag_gp[igp] == 1) {
              // calc g using  g_b and g_a
              assert(g_a <= g_b); // make sure our assumption is correct
              const double dx_ca[3] = DIFF(surface->xsp[isp_a],surface->xsp[isp]);
              const double dx_cb[3] = DIFF(surface->xsp[isp_b],surface->xsp[isp]);
              const double a2 = DOT_PRODUCT(dx_cb,dx_cb);
              const double b2 = DOT_PRODUCT(dx_ca,dx_ca);
              // const double a = sqrt(a2);
              const double b = sqrt(b2);
              const double a_cos_theta = DOT_PRODUCT(dx_ca,dx_cb)/b;
              const double a2_sin2_theta = a2-a_cos_theta*a_cos_theta;
              const double u = g_b-g_a;
              // solve quadratic
              const double q2 = (a2+b2-2.0*b*a_cos_theta);
              const double q1 = 2.0*b*u*(a_cos_theta-b);
              const double q0 = b2*(u*u-a2_sin2_theta);
              const double disc = q1*q1-2.0*q2*q0;
              double t;
              if (disc >= 0.0) {
                t = max((-q1+sqrt(disc))/(2.0*q2),(-q1-sqrt(disc))/(2.0*q2));
              }
              // update g with best available estimate
              if ( (disc >= 0.0) && (u < t) && ( (b*(t-u)/t) > a_cos_theta) && ( (b*(t-u)/t) < (a2/a_cos_theta) ) ) {
                dist_gp[igp] = min(dist_gp[igp],t+g_a);
              }
              else {
                // calc g using simple dist
                dist_gp[igp] = min(dist_gp[igp],g_a+MAG(dx_ca));
                dist_gp[igp] = min(dist_gp[igp],g_b+MAG(dx_cb));
              }
              // update trialHeap
              if (trialHeap.index(igp) >= 0)
                trialHeap.updateDouble(trialHeap.index(igp),dist_gp[igp]);
            }
          }
        }
      }
      else {
        assert(count == 1); // we need at least one alive point
        // find alive point
        int isp_a;
        double g_a;
        FOR_I3 {
          const int isp = surface->spost[ist][i];
          const int igp = gposp[isp];
          if (igp != -1) {
            const double this_g = dist_gp[igp];
            if (flag_gp[igp] == 2) {
              isp_a = isp;
              g_a = this_g;
            }
          }
        }
        // update close points g
        FOR_I3 {
          const int isp = surface->spost[ist][i];
          const int igp = gposp[isp];
          if (igp != -1) {
            if (flag_gp[igp] == 1) {
              // calc g using simple dist
              dist_gp[igp] = min(dist_gp[igp],g_a+DIST(surface->xsp[isp],surface->xsp[isp_a]));
              // update trialHeap
              if (trialHeap.index(igp) >= 0)
                trialHeap.updateDouble(trialHeap.index(igp),dist_gp[igp]);
            }
          }
        }
      }
    }
  }

  void fmmLoop(const int * stosp_i, const int * stosp_v, MinHeap& trialHeap, IntFlag& sp_flag, double* sp_buf, SurfaceShm * surface) {

    // get trial point from close vector...
    pair<double,int> trial = trialHeap.extractMin();
    const int trial_isp = trial.second;
    trialHeap.deleteMin(); // remove from trial
    assert(sp_flag[trial_isp] == 1);
    sp_flag[trial_isp] = 2; // tag as alive
    //cout << "trial " << trial_isp << " g " << sp_buf[trial_isp] << " g again " << trial.first << endl;

    // mark trial_isp neighbors that aren't alive as close...
    for (int tos = stosp_i[trial_isp]; tos != stosp_i[trial_isp+1]; ++tos) {
      const int ist = stosp_v[tos];
      FOR_I3 {
        const int isp = surface->spost[ist][i];
        if (sp_flag[isp] == 0) {
          sp_flag[isp] = 1; // tag as close
          trialHeap.insert(pair<double,int>(HUGE_VAL,isp)); // push onto heap
        }
      }
    }

    /*
    // XXX Can we do this better using up/downHeap???
    for (int index = 0, limit = trialHeap.size(); index < limit; ++index) {
    // marked as close and maintain back-pointers
    sp_flag[trialHeap.heap[index].second] = -index-1; // maintain back pointers;
    }
    */

    // recompute g's for all close trial_isp neighbors...
    for (int tos = stosp_i[trial_isp]; tos != stosp_i[trial_isp+1]; ++tos) {
      const int ist = stosp_v[tos];
      // count number of alive points
      int count = 0;
      FOR_I3 {
        const int isp = surface->spost[ist][i];
        //cout << "isp " << isp << " flag " << sp_flag[isp] << " g " << sp_buf[isp] << endl;
        if (sp_flag[isp] == 2) ++count;
      }
      if (count == 3) {
        continue; // skip tris where all the pts are alive
      }
      else if (count == 2) {
        // find two alive points
        int isp_a, isp_b;
        double g_a = HUGE_VAL, g_b = HUGE_VAL; // a smallest alive value
        FOR_I3 {
          const int isp = surface->spost[ist][i];
          const double this_g = sp_buf[isp];
          if (sp_flag[isp] == 2) {
            if (this_g < g_a) {
              isp_b = isp_a;
              g_b = g_a;
              isp_a = isp;
              g_a = this_g;
            }
            else if (this_g < g_b) {
              isp_b = isp;
              g_b = this_g;
            }
          }
        }
        // update close points g
        FOR_I3 {
          const int isp = surface->spost[ist][i];
          if (sp_flag[isp] == 1) {
            // calc g using  g_b and g_a
            assert(g_a <= g_b); // make sure our assumption is correct
            const double dx_ca[3] = DIFF(surface->xsp[isp_a],surface->xsp[isp]);
            const double dx_cb[3] = DIFF(surface->xsp[isp_b],surface->xsp[isp]);
            const double a2 = DOT_PRODUCT(dx_cb,dx_cb);
            const double b2 = DOT_PRODUCT(dx_ca,dx_ca);
            // const double a = sqrt(a2);
            const double b = sqrt(b2);
            const double a_cos_theta = DOT_PRODUCT(dx_ca,dx_cb)/b;
            const double a2_sin2_theta = a2-a_cos_theta*a_cos_theta;
            const double u = g_b-g_a;
            // solve quadratic
            const double q2 = (a2+b2-2.0*b*a_cos_theta);
            const double q1 = 2.0*b*u*(a_cos_theta-b);
            const double q0 = b2*(u*u-a2_sin2_theta);
            const double disc = q1*q1-2.0*q2*q0;
            double t;
            if (disc >= 0.0) {
              t = max((-q1+sqrt(disc))/(2.0*q2),(-q1-sqrt(disc))/(2.0*q2));
            }
            // update g with best available estimate
            if ( (disc >= 0.0) && (u < t) && ( (b*(t-u)/t) > a_cos_theta) && ( (b*(t-u)/t) < (a2/a_cos_theta) ) ) {
              sp_buf[isp] = min(sp_buf[isp],t+g_a);
            }
            else {
              // calc g using simple dist
              sp_buf[isp] = min(sp_buf[isp], g_a+MAG(dx_ca));
              sp_buf[isp] = min(sp_buf[isp], g_b+MAG(dx_cb));
            }
            // update trialHeap
            if (trialHeap.index(isp) >= 0)
              trialHeap.updateDouble(trialHeap.index(isp),sp_buf[isp]);
          }
        }
      }
      else {
        assert(count == 1); // we need at least one alive point
        // find alive point
        int isp_a;
        double g_a;
        FOR_I3 {
          const int isp = surface->spost[ist][i];
          const double this_g = sp_buf[isp];
          if (sp_flag[isp] == 2) {
            isp_a = isp;
            g_a = this_g;
          }
        }
        // update close points g
        FOR_I3 {
          const int isp = surface->spost[ist][i];
          if (sp_flag[isp] == 1) {
            // calc g using simple dist
            sp_buf[isp] = min(sp_buf[isp], g_a+DIST(surface->xsp[isp],surface->xsp[isp_a]));
            // update trialHeap
            trialHeap.updateDouble(trialHeap.index(isp),sp_buf[isp]);
          }
        }
      }
    }
  }

  void processPistonFF(Param * param) {

    assert(0); // April 2019

    //bool b_suggest = false;
    set<int> zone_indices;

    double dr0=0.0,dn0=0.0;
    double dr1=0.0,dn1=0.0;
    double dr2=0.0,dn2=0.0;
    double dn3=0.0;

    int iarg = 0;
    while (iarg < param->size()) {
      string token = param->getString(iarg++);
      if ((token == "ZONE")||(token == "ZONES")) {
        // expect a comma-delimited set of integers...
        const string zonesCsv = param->getString(iarg++);
        vector<int> tmp;
        MiscUtils::splitCsv(tmp,zonesCsv);
        for (int ii = 0; ii < tmp.size(); ++ii) {
          const int izone = tmp[ii];
          if ((izone >= 0)&&(izone < zoneVec.size())) {
            zone_indices.insert(izone);
          }
          else {
            CWARN("ZONE index out of range: " << izone << "; skipping");
          }
        }
      }
      else if ((token == "ZONE_NAME")||(token == "ZONE_NAMES")||(token == "ZONENAME")||(token == "ZONENAMES")) {
        // expect a comma-delimited list of names...
        const string zoneNamesCsv = param->getString(iarg++);
        vector<string> zone_names;
        MiscUtils::splitCsv(zone_names,zoneNamesCsv);
        for (int ii = 0; ii < zone_names.size(); ++ii) {
          const int izone = getZoneIndex(zone_names[ii]); // get the unique zoneIndex in the partVec
          if (izone >= 0) {
            zone_indices.insert(izone);
          }
          else {
            CWARN("Unrecognized ZONENAME: " << zone_names[ii] << "; skipping");
          }
        }
      }
      else if (token == "ALL") {
        // flag all zones...
        for (int izone = 0; izone < zoneVec.size(); ++izone)
          zone_indices.insert(izone);
      }
      //else if (token == "SUGGEST") {
      //  b_suggest = true;
      //}
      else if (token == "DR0") {
        dr0 = param->getDouble(iarg++);
      }
      else if (token == "DN0") {
        dn0 = param->getDouble(iarg++);
      }
      else if (token == "DR1") {
        dr1 = param->getDouble(iarg++);
      }
      else if (token == "DN1") {
        dn1 = param->getDouble(iarg++);
      }
      else if (token == "DR2") {
        dr2 = param->getDouble(iarg++);
      }
      else if (token == "DN2") {
        dn2 = param->getDouble(iarg++);
      }
      else if (token == "DN3") {
        dn3 = param->getDouble(iarg++);
      }
      else {
        CWARN("Unknown PISTON_FF parameter " << token << "; skipping");
      }
    }

    if (zone_indices.empty()) {
      COUT1("Error: PISTON_FF missing ZONES <comma-delimited-indices> or ZONE_NAMES <comma-delimited-part:zone_names>");;
    }
    else {

      vector<GroupData> groupDataVec;
      buildGroupsForFlaggedZones(groupDataVec,zone_indices);

      if (groupDataVec.size() != 1) {
        COUT1("Error: PISTON_FF requires just one zone at a time for now.");
      }
      else {
        const int igr = 0;
        // for a group to be a piston, it needs just one link...
        if (groupDataVec[igr].nlink != 1) {
          if (mpi_rank == 0) cout << "group " << igr << " has more than one link. Cannot build PISTON_FF" << endl;
        }
        else {
          // we have one link. Now try figuring out the orientation...
          if (mpi_rank == 0) cout << "group " << igr << " has one link, orientation: " <<
                               COUT_VEC(groupDataVec[igr].n_li[0]) << endl;
          double e0[3];
          if (fabs(groupDataVec[igr].n_li[0][0]) > 10.0*max(fabs(groupDataVec[igr].n_li[0][1]),fabs(groupDataVec[igr].n_li[0][2]))) {
            e0[0] = groupDataVec[igr].n_li[0][0];
            e0[1] = 0.0;
            e0[2] = 0.0;
          }
          else if (fabs(groupDataVec[igr].n_li[0][1]) > 10.0*max(fabs(groupDataVec[igr].n_li[0][2]),fabs(groupDataVec[igr].n_li[0][0]))) {
            e0[0] = 0.0;
            e0[1] = groupDataVec[igr].n_li[0][1];
            e0[2] = 0.0;
          }
          else if (fabs(groupDataVec[igr].n_li[0][2]) > 10.0*max(fabs(groupDataVec[igr].n_li[0][0]),fabs(groupDataVec[igr].n_li[0][1]))) {
            e0[0] = 0.0;
            e0[1] = 0.0;
            e0[2] = groupDataVec[igr].n_li[0][2];
          }
          else {
            e0[0] = groupDataVec[igr].n_li[0][0];
            e0[1] = groupDataVec[igr].n_li[0][1];
            e0[2] = groupDataVec[igr].n_li[0][2];
          }
          // and normalize...
          double e0_mag = MAG(e0);
          FOR_I3 e0[i] /= e0_mag;
          if (mpi_rank == 0) cout << " > looks like piston points in direction: " << COUT_VEC(e0) << endl;
          // now build the centroid in that direction...
          SurfaceShm * surface = partVec[groupDataVec[igr].ipart]->surface;
          double x_piston[3] = { 0.0, 0.0, 0.0 };
          double wgt_piston = 0.0;
          int sol_prev = groupDataVec[igr].spoli_i[1]-1;
          for (int sol = groupDataVec[igr].spoli_i[0]; sol != groupDataVec[igr].spoli_i[1]; ++sol) {
            const int isp_prev = groupDataVec[igr].spoli_v[sol_prev];
            sol_prev = sol;
            const int isp = groupDataVec[igr].spoli_v[sol];
            const double dxp_prev[3] = DIFF(surface->xsp[isp_prev],groupDataVec[igr].x_li[0]);
            const double dxp[3] = DIFF(surface->xsp[isp],groupDataVec[igr].x_li[0]);
            // figure out the area projected in the e0 direction...
            const double wgt = CROSS_DOT(dxp,dxp_prev,e0);
            FOR_I3 x_piston[i] += wgt*(surface->xsp[isp_prev][i]+surface->xsp[isp][i]+groupDataVec[igr].x_li[0][i]);
            wgt_piston += wgt;
          }
          FOR_I3 x_piston[i] /= 3.0*wgt_piston;
          // now change some of these to zero's...
          if (mpi_rank == 0) cout << " > exact piston center: " << COUT_VEC(x_piston) << endl;
          const double tol = 1.0E-5*sqrt(wgt_piston);
          FOR_I3 {
            if (fabs(x_piston[i]) < tol)
              x_piston[i] = 0.0;
          }
          if (mpi_rank == 0) cout << " > piston center after tol: " << COUT_VEC(x_piston) << endl;

          if ((dr0 != 0.0)||(dn0 != 0.0)) {
            if (mpi_rank == 0) cout << " > expanding outside loop by DR0: " << dr0 << " DN0: " << dn0 << endl;
            if (mpi_rank_shared == 0) {
              for (int sol = groupDataVec[igr].spoli_i[0]; sol != groupDataVec[igr].spoli_i[1]; ++sol) {
                const int isp = groupDataVec[igr].spoli_v[sol];
                double dxp[3];
                FOR_I3 dxp[i] = surface->xsp[isp][i] - x_piston[i];
                const double dp = DOT_PRODUCT(dxp,e0);
                FOR_I3 dxp[i] -= dp*e0[i];
                const double r = MAG(dxp);
                FOR_I3 surface->xsp[isp][i] += dr0*dxp[i]/r + dn0*e0[i];
              }
            }
            MPI_Barrier(mpi_comm_shared);
          }

          // report the max and min radius of the piston edge...

          double r2_min = HUGE_VAL;
          double r2_max = -HUGE_VAL;
          sol_prev = groupDataVec[igr].spoli_i[1]-1;
          for (int sol = groupDataVec[igr].spoli_i[0]; sol != groupDataVec[igr].spoli_i[1]; ++sol) {
            const int isp = groupDataVec[igr].spoli_v[sol];
            double dxp[3];
            FOR_I3 dxp[i] = surface->xsp[isp][i] - x_piston[i];
            const double dp = DOT_PRODUCT(dxp,e0);
            FOR_I3 dxp[i] -= dp*e0[i];
            double r2 = DOT_PRODUCT(dxp,dxp);
            r2_min = min(r2_min,r2);
            r2_max = max(r2_max,r2);
            // and the edge midpoint...
            const int isp_prev = groupDataVec[igr].spoli_v[sol_prev];
            sol_prev = sol;
            FOR_I3 dxp[i] = 0.5*(surface->xsp[isp][i]+surface->xsp[isp_prev][i]) - x_piston[i];
            r2 = DOT_PRODUCT(dxp,dxp);
            r2_min = min(r2_min,r2);
            r2_max = max(r2_max,r2);
          }
          if (mpi_rank == 0) cout << " > r_min/max: " << sqrt(r2_min) << " " << sqrt(r2_max) << endl;

          // and add the new points and tris...

          const int izone = surface->zoneVec.size();
          surface->zoneVec.push_back(StZone("PISTON_FF"));

          const int nsp0 = surface->nsp;
          const int nst0 = surface->nst;
          const int ned = groupDataVec[igr].spoli_i[1];
          partVec[groupDataVec[igr].ipart]->resize_surface(nsp0+2*ned+1,nst0+5*ned);

          if (mpi_rank_shared == 0) {

            int isp_new = nsp0;
            for (int sol = groupDataVec[igr].spoli_i[0]; sol != groupDataVec[igr].spoli_i[1]; ++sol) {
              const int isp = groupDataVec[igr].spoli_v[sol];
              double dxp[3];
              FOR_I3 dxp[i] = surface->xsp[isp][i] - x_piston[i];
              const double dp = DOT_PRODUCT(dxp,e0);
              FOR_I3 dxp[i] -= dp*e0[i];
              const double r = MAG(dxp);
              // first and second row of points...
              FOR_I3 surface->xsp[isp_new][i] = surface->xsp[isp][i] + dr1*dxp[i]/r + dn1*e0[i];
              FOR_I3 surface->xsp[isp_new+ned][i] = surface->xsp[isp_new][i] + (dr1+dr2)*dxp[i]/r + (dn1+dn2)*e0[i];
              ++isp_new;
            }
            // and the center point...
            FOR_I3 surface->xsp[nsp0+2*ned][i] = x_piston[i] + (dn0+dn1+dn2+dn3)*e0[i];

            // and the triangulation...
            int ist_new = nst0;
            isp_new = nsp0;
            sol_prev = groupDataVec[igr].spoli_i[1]-1;
            for (int sol = groupDataVec[igr].spoli_i[0]; sol != groupDataVec[igr].spoli_i[1]; ++sol) {
              const int isp = groupDataVec[igr].spoli_v[sol];
              const int isp_prev = groupDataVec[igr].spoli_v[sol_prev];
              sol_prev = sol;
              int isp_new_prev = isp_new-1;
              if (isp_new_prev < nsp0) isp_new_prev += ned;
              surface->spost[ist_new][1] = isp_prev;
              surface->spost[ist_new][0] = isp;
              surface->spost[ist_new][2] = isp_new;
              surface->znost[ist_new] = izone;
              ++ist_new;
              surface->spost[ist_new][1] = isp_prev;
              surface->spost[ist_new][0] = isp_new;
              surface->spost[ist_new][2] = isp_new_prev;
              surface->znost[ist_new] = izone;
              ++ist_new;
              // and second layer...
              surface->spost[ist_new][1] = isp_new_prev;
              surface->spost[ist_new][0] = isp_new;
              surface->spost[ist_new][2] = isp_new+ned;
              surface->znost[ist_new] = izone;
              ++ist_new;
              surface->spost[ist_new][1] = isp_new_prev;
              surface->spost[ist_new][0] = isp_new+ned;
              surface->spost[ist_new][2] = isp_new_prev+ned;
              surface->znost[ist_new] = izone;
              ++ist_new;
              // and last...
              surface->spost[ist_new][1] = isp_new_prev+ned;
              surface->spost[ist_new][0] = isp_new+ned;
              surface->spost[ist_new][2] = nsp0+2*ned;
              surface->znost[ist_new] = izone;
              ++ist_new;
              ++isp_new;
            }
            assert(ist_new == surface->nst);

          }
          MPI_Barrier(mpi_comm_shared);

          buildOrRebuildZoneVec();
          b_bbminmax = false; // this will recompute the overall bbox the next time it is req'd

        }

      }

      clearGroups(groupDataVec);

    }


  }

  void processWriteSbin(Param * param,const bool b_help) {
    if (b_help) {
      WUI(INFO,
          "WRITE_SBIN <name> writes the current surface as an sbin file\n" <<
          "for more detail see [$CWIKB:stitch_export]");
      return;
    }

    if (param->size() == 0) {
      WUI(WARN,"expecting WRITE_SBIN <name>");
    }
    else {
      if (PartData::partVec.empty()) {
        WUI(WARN,"no surface to write");
      }
      else if (PartData::partVec.size() > 1) {
        WUI(WARN,"multiple surface writing not implemented");
      }
      else if (PartData::partVec[0]->surface == NULL) {
        WUI(WARN,"first part has no surface");
      }
      else {
        string name = param->getString();
        PartData::partVec[0]->surface->writeBinary(name);
        WUI(INFO,"sbin file " << name << " written");
      }
    }
  }

  // ----------------------------------------------
  // start of namespace-based SolverData...
  // ----------------------------------------------
  
  int step = 0;
  int check_interval = 1;

  int8 ncv_global = 0;
  int ncv = 0,ncv_g = 0;
  int * flag_cv = NULL;
  double (*x_vd)[3] = NULL; // the voronoi forming points, if available
  double *delta_vd = NULL; // the nbr sphere radius that guarantees all voronoi nbrs (and surface tri patches) are present 

  double * vol_cv = NULL;
  double (*x_cv)[3] = NULL;
  
  vector<int> cvopa_i;
  vector<int> paozn; // TODO: replace with the below.first when required
  vector<pair<int,int> > pzozn;
  
  map<const uint8,int> rbiMap;
  vector<int8> rbi_g;

  void clear_new() {

    DELETE(flag_cv);
    DELETE(x_vd);
    DELETE(delta_vd);
    DELETE(vol_cv);
    DELETE(x_cv);
    
  }


  // moved over from VoronoiBuilder...
  
  enum LoadBalanceMode {
    POINTS_LOAD_BALANCE_MODE,
    SURFACE_LOAD_BALANCE_MODE,
    ALL_LOAD_BALANCE_MODE,
  };
  LoadBalanceMode load_balance_mode = POINTS_LOAD_BALANCE_MODE;
  int load_balance_bits = 7; // x y and z (111) directions

  void processLoadBalanceMode(Param * param,const bool b_help) {
    if (b_help) {
      helpLoadBalanceMode();
    }
    else {
      try {
        int iarg = 0;
        int bits = 0;
        while (iarg < param->size()) {
          const string token = MiscUtils::toUpperCase(param->getString(iarg++));
          if (token == "POINTS") {
            load_balance_mode = POINTS_LOAD_BALANCE_MODE;
            WUI(INFO,"LOAD_BALANCE_MODE set to POINTS (default)");
          }
          else if (token == "SURFACE") {
            load_balance_mode = SURFACE_LOAD_BALANCE_MODE;
            WUI(INFO,"LOAD_BALANCE_MODE set to SURFACE");
          }
          else if (token == "ALL") {
            load_balance_mode = ALL_LOAD_BALANCE_MODE;
            WUI(INFO,"LOAD_BALANCE_MODE set to ALL");
          }
          else if (token == "X") {
            bits |= 1;
          }
          else if (token == "Y") {
            bits |= 2;
          }
          else if (token == "Z") {
            bits |= 4;
          }
          else {
            WUI(WARN,"unrecognized LOAD_BALANCE_MODE token: " << token);
            helpLoadBalanceMode();
          }
        }
        if (bits > 0)
          load_balance_bits = bits;
      }
      catch(int e) {
        WUI(WARN,"Error parsing LOAD_BALANCE_MODE tokens");
        helpLoadBalanceMode();
      }
    }
  }

  void helpLoadBalanceMode() {
    WUI(INFO,
        "LOAD_BALANCE_MODE sets the base PADT partition of space. Examples:\n" <<
        "LOAD_BALANCE_MODE POINTS # (default) partitions the domain using hcp and part points.\n" <<
        "LOAD_BALANCE_MODE SURFACE X # partitions the domain in x using surface points.\n" <<
        "LOAD_BALANCE_MODE ALL X Y # partitions the domain in x and y using surface, hcp and part points."
	);
  }

  void clearPoints() {

    DELETE(x_vd);
    DELETE(delta_vd);
    DELETE(flag_cv);

  }

  void ensurePoints() {

    // leave if we alread have forming points...
    if (x_vd) { 
      assert(delta_vd);
      assert(flag_cv);
      return;
    }

    COUT1("setPoints()");

    // count the points...

    int np = 0;
    if (checkParam("SKIP_IS_INSIDE_CHECK")) {

      if (mpi_rank == 0) cout << "WARNING: SKIP_IS_INSIDE_CHECK found. Assuming all part points are inside" << endl;

      for (int ipart = 0; ipart < partVec.size(); ++ipart) {
        if (partVec[ipart]->pts) {
          // leave delta positive, and include all points...
          np += partVec[ipart]->pts->np;
        }
      }

    }
    else {

      prepareIsInside();

      for (int ipart = 0; ipart < partVec.size(); ++ipart) {
        if (partVec[ipart]->pts) {
          for (int ip = 0; ip < partVec[ipart]->pts->np; ++ip) {
            // ensure delta is positive and non-zero...
            const double r = partVec[ipart]->pts->delta[ip];
            assert(r > 0.0);
            // we are going to include this point in the build, if and only if it
            // is NOT inside the FF of every part later than us in the list (these parts take priority)
            // AND it is not inside the solid of any part. Since we have checked the FF of the
            // later parts, we can just check the solid of earlier parts...
            for (int ipart2 = ipart+1; ipart2 < partVec.size(); ++ipart2) {
              if (partVec[ipart2]->hasFF()) {
                if (partVec[ipart2]->isInsideFF(partVec[ipart]->pts->xp[ip])) {
                  partVec[ipart]->pts->delta[ip] = -r;
                  break;
                }
              }
              else {
                if (partVec[ipart2]->isInsideSolid(partVec[ipart]->pts->xp[ip],false)) { // ipart2 != 0
                  partVec[ipart]->pts->delta[ip] = -r;
                  break;
                }
              }
            }
            if (partVec[ipart]->pts->delta[ip] > 0.0) {
              for (int ipart2 = 0; ipart2 < ipart; ++ipart2) {
                // the passed bool (ipart2==0) controls how this function behaves when no
                // intersection with the solid is found. See details in Part::isInsideSolid...
                if (partVec[ipart2]->isInsideSolid(partVec[ipart]->pts->xp[ip],ipart2==0)) {
                  partVec[ipart]->pts->delta[ip] = -r;
                  break;
                }
              }
            }
            if (partVec[ipart]->pts->delta[ip] > 0.0) ++np;
          }
        }
      }

    }

    // all the hcp points get included...
    if (hcpPts) 
      np += hcpPts->np;

    // build the striped distribution as input to the load balancer...
    int8 * ipora = NULL;
    buildXora(ipora,np);

    ncv_global = ipora[mpi_size];

    if (mpi_rank == 0) cout << " > ncv_global: " << ncv_global << endl;

    // take a copy of the points...

    double (*xp)[3] = new double[np][3];
    double *deltap = new double[np];
    int *flagp = new int[np];
    int ip = 0;
    for (int ipart = 0; ipart < partVec.size(); ++ipart) {
      if (partVec[ipart]->pts) {
        for (int ip_ = 0; ip_ < partVec[ipart]->pts->np; ++ip_) {
          if (partVec[ipart]->pts->delta[ip_] > 0.0) {
            FOR_I3 xp[ip][i] = partVec[ipart]->pts->xp[ip_][i];
            deltap[ip]       = partVec[ipart]->pts->delta[ip_];
            flagp[ip]        = ipart;
            ++ip;
          }
          else {
            // flip back to positive...
            partVec[ipart]->pts->delta[ip_] = -partVec[ipart]->pts->delta[ip_];
          }
        }
      }
    }
    if (hcpPts) {
      for (int ip_ = 0; ip_ < hcpPts->np; ++ip_) {
        FOR_I3 xp[ip][i] = hcpPts->xp[ip_][i];
        deltap[ip]       = hcpPts->delta[ip_];
        flagp[ip] = partVec.size(); // 1 bigger than the last part
        ++ip;
      }
    }
    assert(ip == np);

    // load balance...

    int8 * ip_global = NULL;
    if (load_balance_mode == POINTS_LOAD_BALANCE_MODE) {
      ip_global = new int8[np];
      for (int ip = 0; ip < np; ++ip)
        ip_global[ip] = ipora[mpi_rank] + ip;

      int np_new = np;
      repartXcvPadt(xp,ip_global,np_new,mpi_comm,load_balance_bits);
      assert(x_vd == NULL); x_vd = xp; xp = NULL;
      ncv = ncv_g = np_new;
    }
    else {
      assert((load_balance_mode == SURFACE_LOAD_BALANCE_MODE)||(load_balance_mode == ALL_LOAD_BALANCE_MODE));
      int nsp;
      double (*xsp)[3] = NULL;
      for (int iter = 0; iter < 2; ++iter) {
        nsp = 0;
        for (int ipart = 0; ipart < partVec.size(); ++ipart) {
          if (partVec[ipart]->surface) {
            SurfaceShm* surface = partVec[ipart]->surface;
            for (int isp = mpi_rank; isp < surface->nsp; isp += mpi_size) {
              if (iter == 1)
                FOR_I3 xsp[nsp][i] = surface->xsp[isp][i];
              ++nsp;
            }
          }
        }
        if (load_balance_mode == ALL_LOAD_BALANCE_MODE) {
          if (iter == 0) {
            nsp += np;
          }
          else {
            for (int ip = 0; ip < np; ++ip) {
              FOR_I3 xsp[nsp][i] = xp[ip][i];
              ++nsp;
            }
          }
        }
        if (iter == 0)
          xsp = new double[nsp][3];
      }
      int * part_sp = new int[nsp];
      calcCvPartPadt(part_sp,xsp,nsp,mpi_comm,load_balance_bits);

      double (*my_bbminmax)[6] = new double[mpi_size][6];
      FOR_RANK FOR_I6 my_bbminmax[rank][i] = HUGE_VAL;
      for (int isp = 0; isp < nsp; ++isp) {
        const int rank = part_sp[isp]; assert((part_sp[isp] >= 0)&&(part_sp[isp] < mpi_size));
        FOR_I3 my_bbminmax[rank][i  ] = min(my_bbminmax[rank][i  ], xsp[isp][i]);
        FOR_I3 my_bbminmax[rank][3+i] = min(my_bbminmax[rank][3+i],-xsp[isp][i]);
      }
      delete[] part_sp;
      delete[] xsp;

      double (*bbminmax)[6] = new double[mpi_size][6];
      MPI_Allreduce((double*)my_bbminmax,(double*)bbminmax,6*mpi_size,MPI_DOUBLE,MPI_MIN,mpi_comm);
      delete[] my_bbminmax;

      // now bbminmax stores the bounding box of the load-balanced partition...

      int* part_ip = new int[np];
      for (int ip = 0; ip < np; ++ip) {
        part_ip[ip] = -1;
        double min_dist = HUGE_VAL;
        FOR_RANK {
          double x_bb[3]; FOR_I3 x_bb[i] = 0.5*(bbminmax[rank][i]-bbminmax[rank][3+i]);
          double this_dist = DIST(x_bb,xp[ip]);
          if (this_dist < min_dist) {
            min_dist = this_dist;
            part_ip[ip] = rank;
          }
        }
        assert((part_ip[ip] >= 0)&&(part_ip[ip] < mpi_size));
      }
      delete[] bbminmax;

      int* send_count = new int[mpi_size];
      FOR_RANK send_count[rank] = 0;
      int* send_disp = new int[mpi_size];
      int* send_buf_int = NULL;
      double* send_buf_double = NULL;
      for (int iter = 0; iter < 2; ++iter) {
        for (int ip = 0; ip < np; ++ip) {
          const int rank = part_ip[ip];
          if (iter == 0) {
            ++send_count[rank];
          }
          else {
            send_buf_int[send_disp[rank]] = ip;
            FOR_I3 send_buf_double[send_disp[rank]*3+i] = xp[ip][i];
            ++send_disp[rank];
          }
        }

        send_disp[0] = 0;
        for (int rank = 1; rank < mpi_size; ++rank)
          send_disp[rank] = send_count[rank-1] + send_disp[rank-1];

        if (iter == 0) {
          const int send_count_sum = send_disp[mpi_size-1] + send_count[mpi_size-1];
          send_buf_int    = new int[send_count_sum];
          send_buf_double = new double[3*send_count_sum];
        }
      }
      delete[] xp; xp = NULL;
      delete[] part_ip;

      // set up recv-side stuff...

      int *recv_count = new int[mpi_size];
      MPI_Alltoall(send_count,1,MPI_INT,recv_count,1,MPI_INT,mpi_comm);

      int *recv_disp = new int[mpi_size];
      recv_disp[0] = 0;
      for (int rank = 1; rank < mpi_size; ++rank)
        recv_disp[rank] = recv_count[rank-1] + recv_disp[rank-1];
      const int recv_count_sum = recv_disp[mpi_size-1] + recv_count[mpi_size-1];

      int * recv_buf_int = new int[recv_count_sum];
      MPI_Alltoallv(send_buf_int,send_count,send_disp,MPI_INT,
		    recv_buf_int,recv_count,recv_disp,MPI_INT,mpi_comm);
      delete[] send_buf_int;

      FOR_RANK {
        send_count[rank] *= 3;
        send_disp[rank]  *= 3;
        recv_count[rank] *= 3;
        recv_disp[rank]  *= 3;
      }

      double * recv_buf_double = new double[recv_count_sum*3];
      MPI_Alltoallv(send_buf_double,send_count,send_disp,MPI_DOUBLE,
		    recv_buf_double,recv_count,recv_disp,MPI_DOUBLE,mpi_comm);
      delete[] send_count;
      delete[] send_disp;
      delete[] send_buf_double;;

      FOR_RANK {
        recv_count[rank] /= 3;
        recv_disp[rank]  /= 3;
      }

      ncv = 0;
      assert(ip_global == NULL); ip_global = new int8[recv_count_sum];
      assert(x_vd == NULL); x_vd = new double[recv_count_sum][3];
      FOR_RANK {
        for (int irecv = recv_disp[rank], lim = recv_disp[rank]+recv_count[rank]; irecv < lim; ++irecv) {
          ip_global[ncv] = recv_buf_int[irecv]+ipora[rank];
          FOR_I3 x_vd[ncv][i] = recv_buf_double[irecv*3+i];
          ++ncv;
        }
      }
      assert(ncv == recv_count_sum);
      ncv_g = ncv;
      delete[] recv_count;
      delete[] recv_disp;
      delete[] recv_buf_int;
      delete[] recv_buf_double;
    }

    reorderXcv(x_vd,ip_global,ncv);

    // use dde to pull the other data from the striped into the load balanced partition...
    DistributedDataExchanger * dde = new DistributedDataExchanger(ip_global,ncv,ipora);
    delete[] ipora;
    delete[] ip_global;

    assert(delta_vd == NULL); delta_vd = new double[ncv];
    dde->pull(delta_vd,deltap);
    delete[] deltap; deltap = NULL;

    assert(flag_cv == NULL); flag_cv = new int[ncv];
    dde->pull(flag_cv,flagp);
    delete[] flagp; flagp = NULL;

    delete dde;

    // populate the global arrays...

    const int npart = partVec.size();
    cvopa_i.clear();
    cvopa_i.resize(npart+2); // last for hcp points
    cvopa_i[0] = 0;
    for (int ipart = 0; ipart <= npart; ++ipart)
      cvopa_i[ipart+1] = 0;
    FOR_ICV {
      const int ipart = flag_cv[icv];
      assert((ipart >= 0)&&(ipart <= npart)); 
      ++cvopa_i[ipart+1];
    }
    for (int ipart = 0; ipart <= npart; ++ipart)
      cvopa_i[ipart+1] += cvopa_i[ipart];
    assert(cvopa_i[npart+1] == ncv);

    // need to reorder points to cvopa ordered...

    assert(xp == NULL); xp = x_vd; x_vd = new double[ncv][3];
    assert(deltap == NULL); deltap = delta_vd; delta_vd = new double[ncv];
    assert(flagp == NULL); flagp = flag_cv; flag_cv = new int[ncv];
    FOR_ICV {
      const int ipart = flagp[icv];
      assert((ipart >= 0)&&(ipart <= npart)); 
      flag_cv[cvopa_i[ipart]] = ipart;
      delta_vd[cvopa_i[ipart]] = deltap[icv];
      FOR_I3 x_vd[cvopa_i[ipart]][i] = xp[icv][i];
      ++cvopa_i[ipart];
    }
    delete[] xp;
    delete[] deltap;
    delete[] flagp;

    for (int ipart = npart; ipart > 0; --ipart) 
      cvopa_i[ipart] = cvopa_i[ipart-1];
    cvopa_i[0] = 0;

  }

} // namespace PartData
