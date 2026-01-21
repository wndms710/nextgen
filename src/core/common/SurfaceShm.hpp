#ifndef _SURFACE_SHM_HPP_
#define _SURFACE_SHM_HPP_

#include "MiscUtils.hpp"
#include "GaussQuadrature.hpp"
#include "StZone.hpp"
#include "GeodesicSphere.hpp"
#include "GeomUtils.hpp"
#include "Utils.hpp"
#include "PeriodicData.hpp"
//#include "Adt2dShmNew.hpp"
#include "Adt2d.hpp"
#include <stack>

class SurfaceShm {
private:

  bool b_bbminmax;
  double bbminmax[6]; // stored as xmin,ymin,zmin,-xmax,-ymax,-zmax

  bool b_volume_and_xc;
  double volume; // the volume is computed at the same time as the centroid.
  double xc[3];

  bool b_mom;
  double mom_diag[3],mom_offd[3];

  bool b_area;
  double area;

public:

  // data stored when groups are built/rebuilt... 
  
  set<int> group_zones;
  int ngr, ngr_global;
  // these store loops by repeating the value of the first point...
  vector<int> spogr_i; // index of first isp in group.
  vector<int> spogr_v; // value, i.e. isp
  vector<pair<int,int> > nlogr; // loop and periodic loop count for each group...

  // see comment below...
  // variables for point-is-inside surface methods
  /*
    Adt2dShmNew * adt2d;
    double surface_x0[3],surface_dx;
    double surface_e0[3],surface_e1[3],surface_e2[3];
    class XpIntersection {
    public:
    double xp;
    int level;
    int sign;
    bool operator<(const XpIntersection& rhs) const { return (xp < rhs.xp); }
    };
  */

  // FH comments: Adt2dShmNew.hpp not initializing properly in all cases, so 
  // moving to local sub-surface approach that requires initialization by
  // the full bounding box of points being considered...

  Adt2d<double> * adt2d;
  vector<int> adt2d_st_vec; // lookup table for the surface tris in the local adt2d
  double adt2d_bbmin[2]; // 2d bounding box for the local adt2d: can be much larger
  double adt2d_bbmax[2]; 
  
  // we carry these sizes so we can use standard
  // allgather from the nodes...
  int nsp_max;
  int nst_max;

  int nsp;
  double (*xsp)[3];
  int nst;
  int (*spost)[3];
  int *znost;
  vector<StZone> zoneVec;
  uint8 *pbi;

  // inverse of spost:
  // build this with ensureStosp...
  int *stosp_i;
  int *stosp_v;

  double (*n_sp)[3]; // nodal normal (not built by default)
  double (*dx_sp)[3]; // set during RUN_DIAGNOSTICS if requested
  double *dist_sp; // set during RUN_DIAGNOSTICS if requested
  double *dist_st; // set during RUN_DIAGNOSTICS if requested
  int *flag_st;
  int *flag_sp;

  // not built by default. stost holds the
  // 3 neighboring tris...
  int (*stost)[3];

  uint8 (*teost)[3]; // replace ALL structures above eventually

  // zone-of-zone graph for walking the surface...
  int *znozn_i;
  int *znozn_v;
  
  SurfaceShm() {
    init();
  }

  SurfaceShm(const string& filename) {
    init();
    const int ierr = readBinary(filename);
    if (ierr != 0)
      throw(0);
  }

  void init() {
    xsp = NULL;
    spost = NULL;
    znost = NULL;
    b_bbminmax = false;
    b_volume_and_xc = false;
    b_mom = false;
    pbi = NULL;
    stosp_i = stosp_v = NULL;
    n_sp = NULL;
    dx_sp = NULL;
    dist_sp = NULL;
    dist_st = NULL;
    flag_st = NULL;
    flag_sp = NULL;
    stost = NULL;
    adt2d = NULL;
    teost = NULL;
    // default latice unit vectors: x,y,z
    //surface_e0[0] = 1.0; surface_e0[1] = 0.0; surface_e0[2] = 0.0;
    //surface_e1[0] = 0.0; surface_e1[1] = 1.0; surface_e1[2] = 0.0;
    //surface_e2[0] = 0.0; surface_e2[1] = 0.0; surface_e2[2] = 1.0;
    znozn_i = NULL;
    znozn_v = NULL;
  }    

  void init(const int nsp,const int nst) {
    this->nsp = nsp_max = nsp;
    assert(xsp == NULL);
    CTI_Mmap_rw(xsp,nsp_max); // _rw
    this->nst = nst_max = nst;
    assert(spost == NULL);
    CTI_Mmap(spost,nst_max);
    assert(znost == NULL);
    CTI_Mmap(znost,nst_max);
    MPI_Barrier(mpi_comm_shared);
  }
  
  void resize(const int nsp_new,const int nst_new) {

    // clear any stosp_i/v...
    // rebuild this later if required...

    if (stosp_i) {
      assert(stosp_v);
      CTI_Munmap(stosp_v,stosp_i[nsp]);
      CTI_Munmap(stosp_i,nsp+1);
      stosp_i = stosp_v = NULL;
    }

    // clear any teost as well.
    // rebuild this later if required...
    
    if (teost) {
      CTI_Munmap(teost,nst);
      teost = NULL;
    }

    // clear any stost as well.
    // rebuild this later if required...

    if (stost) {
      CTI_Munmap(stost,nst);
      stost = NULL;
    }
    if (n_sp) {
      CTI_Munmap(n_sp,nsp);
      n_sp = NULL;
    }
    if (dx_sp) {
      CTI_Munmap(dx_sp,nsp);
      dx_sp = NULL;
    }
    if (dist_sp) {
      CTI_Munmap(dist_sp,nsp);
      dist_sp = NULL;
    }
    if (dist_st) {
      CTI_Munmap(dist_st,nst);
      dist_st = NULL;
    }
    if (flag_st) {
      CTI_Munmap(flag_st,nst);
      flag_st = NULL;
    }
    if (flag_sp) {
      CTI_Munmap(flag_sp,nsp);
      flag_sp = NULL;
    }
    
    // sp data...

    if (nsp_new > nsp_max) {

      double (*xsp_new)[3] = NULL;
      CTI_Mmap_rw(xsp_new,nsp_new); // _rw
      if (xsp != NULL) {
        if (mpi_rank_shared == 0)
          memcpy(xsp_new,xsp,sizeof(double)*nsp*3);
        MPI_Barrier(mpi_comm_shared);
        CTI_Munmap(xsp,nsp_max);
      }
      xsp = xsp_new;

      assert(pbi == NULL); // no resize on periodic (yet)

      nsp_max = nsp_new;

    }
    nsp = nsp_new;

    // st data...

    if (nst_new > nst_max) {

      int (*spost_new)[3] = NULL;
      CTI_Mmap(spost_new,nst_new);
      if (spost != NULL) {
        if (mpi_rank_shared == 0)
          memcpy(spost_new,spost,sizeof(int)*nst*3);
        MPI_Barrier(mpi_comm_shared);
        CTI_Munmap(spost,nst_max);
      }
      spost = spost_new;

      int *znost_new = NULL;
      CTI_Mmap(znost_new,nst_new);
      if (znost != NULL) {
        if (mpi_rank_shared == 0)
          memcpy(znost_new,znost,sizeof(int)*nst);
        MPI_Barrier(mpi_comm_shared);
        CTI_Munmap(znost,nst_max);
      }
      znost = znost_new;

      nst_max = nst_new;

    }
    nst = nst_new;


  }

  uint8 packTeost(const int ist,const int i,const int bit) {
    assert((bit == 0)||(bit == (1<<0))||(bit == (1<<1))); // will fail for multiple periodicity: just add individual bits up to 6
    assert((i >= 0)&&(i < 3));
    return ((uint8(bit)<<52)|(uint8(i)<<50)|uint8(ist)); 
  }

  void unpackTeost(int& ist,int& i,int& bit,const uint8 teost) {
    ist = (teost&MASK_32BITS);
    i   = ((teost>>50)&MASK_2BITS);
    bit = ((teost>>52)&MASK_6BITS);
    assert(packTeost(ist,i,bit) == teost);
  }

  ~SurfaceShm() {
    if (xsp != NULL) CTI_Munmap(xsp,nsp_max);
    if (spost != NULL) CTI_Munmap(spost,nst_max);
    if (znost != NULL) CTI_Munmap(znost,nst_max);
    if (pbi != NULL) CTI_Munmap(pbi,nsp_max);
    if (stosp_i) {
      assert(stosp_v);
      CTI_Munmap(stosp_v,stosp_i[nsp]);
      CTI_Munmap(stosp_i,nsp+1);
      stosp_i = stosp_v = NULL;
    }
    if (stost) {
      CTI_Munmap(stost,nst);
      stost = NULL;
    }
    if (n_sp) {
      CTI_Munmap(n_sp,nsp);
      n_sp = NULL;
    }
    if (dx_sp) {
      CTI_Munmap(dx_sp,nsp);
      dx_sp = NULL;
    }
    if (dist_sp) {
      CTI_Munmap(dist_sp,nsp);
      dist_sp = NULL;
    }
    if (dist_st) {
      CTI_Munmap(dist_st,nst);
      dist_st = NULL;
    }
    if (flag_st) {
      CTI_Munmap(flag_st,nst);
      flag_st = NULL;
    }
    if (flag_sp) {
      CTI_Munmap(flag_sp,nsp);
      flag_sp = NULL;
    }
    if (adt2d) {
      delete adt2d;
      adt2d = NULL;
    }
    // could have these shared as well...
    if (znozn_i) delete[] znozn_i;
    if (znozn_v) delete[] znozn_v;
  }

  bool checkPbi() const { return( pbi != NULL ); }

  uint8 getPbi(const int isp) const {
    assert((isp >= 0)&&(isp < nsp));
    if (pbi) {
      return( pbi[isp] );
    }
    else {
      return uint8(isp);
    }
  }

  void flagZonesMatching(const string& names) {
    for (int izone = 0; izone < zoneVec.size(); ++izone)
      zoneVec[izone].flag = 0;
    // expect a comma-delimited string...
    vector<string> zoneNameVec;
    MiscUtils::splitCsv(zoneNameVec,names);
    // still allow for the possibility that each token uses a wildcard...
    for (int ii = 0; ii < zoneNameVec.size(); ++ii) {
      for (int izone = 0; izone < zoneVec.size(); ++izone) {
        if (zoneVec[izone].flag == 0) {
          if (MiscUtils::strcmp_wildcard(zoneVec[izone].getName(),zoneNameVec[ii]))
            zoneVec[izone].flag = 1;
        }
      }
    }
  }

  int readBinary(const string& filename) {

    if (mpi_rank == 0) cout << "SurfaceShm::readBinary \"" << filename << "\"..." << endl;

    // parallel read into shm...
    const int zone_name_len_max = 128;

    int size_buf[8] = { 0, 0, 0, 0, 0, 0, 0, 0 };
    int my_ierr = 0;
    double (*xsp_split)[3] = NULL;
    int (*spost_split)[3] = NULL;
    int *znost_split = NULL;
    char *char_buf = NULL;
    int nzones = -1;

    int npt = 0;
    int periodic_kind[3];
    double periodic_data[3][3];
    int npz = 0;
    int *periodic_zone_array = NULL;
    uint2 *periodic_zone_bits = NULL;
    uint8 *pbi_split = NULL;

    if (mpi_rank_shared == 0) {

      // ---------------------------------------------------------------
      // only one rank per node participates in the parallel read...
      // ---------------------------------------------------------------

      bool byte_swap = false;
      MPI_File fh;
      char dummy[128];
      sprintf(dummy,"%s",filename.c_str());
      my_ierr = MPI_File_open(mpi_comm_internode,dummy,MPI_MODE_RDONLY,MPI_INFO_NULL,&fh);
      if (my_ierr != 0) {
	if (mpi_rank == 0) cout << "Error: cannot open binary file \"" << filename << "\"" << endl;
        my_ierr = -1; // force -1 to ensure the MPI_MIN below catches this._
      }
      else {

        int itmp[2];
        if (mpi_rank == 0) {
          // first 2 ints are: 0. version, 1. nzones
          MPI_File_read(fh,itmp,2,MPI_INT,MPI_STATUS_IGNORE);
        }
        MPI_Bcast(itmp,2,MPI_INT,0,mpi_comm_internode);

        if ( (itmp[0] < 1)||(itmp[0] > 10)) {

          ByteSwap::byteSwap(itmp,2);
          if ( (itmp[0] < 1) ||(itmp[0] > 10)) {
            if (mpi_rank == 0) cout << "Error: file does not start as expected \"" << filename << "\"" << endl;
            MPI_File_close(&fh);
            my_ierr = -1; // force -1 to ensure the MPI_MIN below catches this.
          } else {

            byte_swap = true;

          }
        }

        if ( my_ierr == 0) {

          MPI_Offset offset = int_size*2;

          // set nzones...
          nzones = size_buf[0] = itmp[1];
          if (mpi_rank == 0) cout << " > nzones: " << nzones << endl;

          // zone names are stored using length (int), then chars.
          // Read these into a char_buf on rank 0 only...
          if (mpi_rank == 0) {
            zoneVec.resize(nzones); // mpi_rank 0 only for now
            char_buf = new char[zone_name_len_max*nzones];
            memset(char_buf,0,zone_name_len_max*nzones); 
            int length;
            for (int izone = 0; izone < nzones; ++izone) {
              MPI_File_read(fh,&length,1,MPI_INT,MPI_STATUS_IGNORE);
              if ( byte_swap) length = ByteSwap::byteSwap(length);
              assert((length > 0)&&(length+1 <= zone_name_len_max));
              MPI_File_read(fh,char_buf+zone_name_len_max*izone,length,MPI_CHAR,MPI_STATUS_IGNORE);
              offset += int_size + length;
              char_buf[zone_name_len_max*izone+length] = '\0';
              zoneVec[izone].setName(char_buf+zone_name_len_max*izone);
              cout << " > zone: " << izone << " \"" << zoneVec[izone].getName() << "\"" << endl;
            }
          }

          // next is the number of surface points...
          int8 i8tmp[2];
          if (mpi_rank == 0) {
            MPI_File_read(fh,&nsp,1,MPI_INT,MPI_STATUS_IGNORE);
            if ( byte_swap) nsp = ByteSwap::byteSwap(nsp);
            offset += int_size;
            cout << " > nsp: " << nsp << endl;
            i8tmp[0] = nsp;
            i8tmp[1] = offset;
          }
          MPI_Bcast(i8tmp,2,MPI_INT8,0,mpi_comm_internode);
          nsp = size_buf[1] = i8tmp[0];
          offset = i8tmp[1];

          // read in the points in parallel using an equal split with the last
          // reader padded if necessary...

          int nsp_split = nsp/mpi_size_internode;
          if (nsp%mpi_size_internode != 0) nsp_split += 1;
          size_buf[2] = nsp_split;
          size_buf[3] = mpi_size_internode;
          const int isp0 = min(nsp,nsp_split*mpi_rank_internode);
          const int isp1 = min(nsp,nsp_split*(mpi_rank_internode+1));
          const int my_nsp = isp1 - isp0;
          assert(my_nsp <= nsp_split);
          xsp_split = new double[nsp_split][3];
          MPI_File_read_at_all(fh,offset+double_size*nsp_split*mpi_rank_internode*3,
                               (double*)xsp_split,my_nsp*3,MPI_DOUBLE,MPI_STATUS_IGNORE);
          if ( byte_swap) ByteSwap::byteSwap(xsp_split,my_nsp);
          offset += double_size*nsp*3;

          if (mpi_rank == 0) {
            MPI_File_read_at(fh,offset,&nst,1,MPI_INT,MPI_STATUS_IGNORE);
            if ( byte_swap) nst = ByteSwap::byteSwap(nst);
            cout << " > nst: " << nst << endl;
          }
          offset += int_size;
          MPI_Bcast(&nst,1,MPI_INT,0,mpi_comm_internode);
          size_buf[4] = nst;

          int nst_split = nst/mpi_size_internode;
          if (nst%mpi_size_internode != 0) nst_split += 1;
          size_buf[5] = nst_split;
          int ist0 = min(nst,nst_split*mpi_rank_internode);
          int ist1 = min(nst,nst_split*(mpi_rank_internode+1));
          int my_nst = ist1 - ist0;
          spost_split = new int[nst_split][3];
          MPI_File_read_at_all(fh,offset+int_size*nst_split*mpi_rank_internode*3,
                               (int*)spost_split,my_nst*3,MPI_INT,MPI_STATUS_IGNORE);
          if ( byte_swap) ByteSwap::byteSwap(spost_split,my_nst);
          offset += int_size*nst*3;

          znost_split = new int[nst_split];
          MPI_File_read_at_all(fh,offset+int_size*nst_split*mpi_rank_internode,
                               znost_split,my_nst,MPI_INT,MPI_STATUS_IGNORE);
          if ( byte_swap) ByteSwap::byteSwap(znost_split,my_nst);
          offset += int_size*nst;

          // any periodic?...
          // NOTE: a filesize check seems to be the only robust way to do this. I tried to read the
          // next int (npt) from the file as done in the readBinary routine of Surface,
          // and this sets status.MPI_ERROR in some way that is not obvious nor general. Eventually
          // the sbin file needs to indicate this more robustly.

          MPI_Offset offset_check;
          MPI_File_get_size(fh,&offset_check); // TODO: should everyone do this? or read and bcast?

          if (offset_check > offset) {

            if (mpi_rank == 0) {
              MPI_File_read_at(fh,offset,&npt,1,MPI_INT,MPI_STATUS_IGNORE);
              if ( byte_swap) npt = ByteSwap::byteSwap(npt);
              cout << " > sbin has periodic data: " << npt << endl;
              assert((npt > 0)&&(npt <= 3));
            }
            offset += int_size;
            MPI_Bcast(&npt,1,MPI_INT,0,mpi_comm_internode);
            size_buf[6] = npt;

            // we have "npt" pairs of periodic transforms. Since the periodic
            // transform is managed by the PartData now, we need to add there or check
            // for consistency...
            // note: periodic data is one int and 3 doubles...
            if (mpi_rank == 0) {
              for (int ipt = 0; ipt < npt; ++ipt) {
                MPI_File_read_at(fh,offset,periodic_kind+ipt,1,MPI_INT,MPI_STATUS_IGNORE);
                if (byte_swap) periodic_kind[ipt] = ByteSwap::byteSwap(periodic_kind[ipt]);
                offset += int_size;
                MPI_File_read_at(fh,offset,periodic_data[ipt],3,MPI_DOUBLE,MPI_STATUS_IGNORE);
                if ( byte_swap) ByteSwap::byteSwap(periodic_data[ipt],3);
                offset += double_size*3;
                cout << " > periodic kind: " << periodic_kind[ipt] << " periodic_data: " << COUT_VEC(periodic_data[ipt]) << endl;
              }
            }
            else {
              offset += (int_size+double_size*3)*npt;
            }
            MPI_Bcast(periodic_kind,npt,MPI_INT,0,mpi_comm_internode);
            MPI_Bcast(periodic_data,npt*3,MPI_DOUBLE,0,mpi_comm_internode);

            if (mpi_rank == 0) {
              MPI_File_read_at(fh,offset,&npz,1,MPI_INT,MPI_STATUS_IGNORE);
              if ( byte_swap) npz = ByteSwap::byteSwap(npz);
              cout << " > npz: " << npz << endl;
              // assert((npz >= 2*npt)&&(npz <= nzones));  // edge-based periodicity will not have any flagged zones, so need to adjust criterion
              assert(npz <= nzones);
            }
            offset += int_size;
            MPI_Bcast(&npz,1,MPI_INT,0,mpi_comm_internode);
            size_buf[7] = npz;

            periodic_zone_array = new int[npz];
            periodic_zone_bits = new uint2[npz];
            if (mpi_rank == 0) {
              MPI_File_read_at(fh,offset,periodic_zone_array,npz,MPI_INT,MPI_STATUS_IGNORE);
              if (byte_swap) ByteSwap::byteSwap(periodic_zone_array,npz);
              MPI_File_read_at(fh,offset+int_size*npz,periodic_zone_bits,npz,MPI_UINT2,MPI_STATUS_IGNORE);
              if (byte_swap) ByteSwap::byteSwap(periodic_zone_bits,npz);
              for (int ipz = 0; ipz < npz; ++ipz) {
                cout << " > zone: " << periodic_zone_array[ipz] << " bits: ";
                for (int ibit = 5; ibit >= 0; --ibit) {
                  if (periodic_zone_bits[ipz]&(uint2(1)<<ibit)) {
                    cout << "1";
                  }
                  else {
                    cout << "0";
                  }
                }
                cout << endl;
              }
            }
            offset += (int_size+sizeof(uint2))*npz;
            MPI_Bcast(periodic_zone_array,npz,MPI_INT,0,mpi_comm_internode);
            MPI_Bcast(periodic_zone_bits,npz,MPI_UINT2,0,mpi_comm_internode);

            int npbi;
            if (mpi_rank == 0) {
              MPI_File_read_at(fh,offset,&npbi,1,MPI_INT,MPI_STATUS_IGNORE);
              if ( byte_swap) npbi = ByteSwap::byteSwap(npbi);
              cout << " > npbi: " << npbi << endl;
            }
            offset += int_size;
            MPI_Bcast(&npbi,1,MPI_INT,0,mpi_comm_internode);

            int * isp_array = new int[npbi];
            int * isp_p_array = new int[npbi];
            uint2 * isp_bits_array = new uint2[npbi];
            if (mpi_rank == 0) {
              MPI_File_read_at(fh,offset,isp_array,npbi,MPI_INT,MPI_STATUS_IGNORE);
              if ( byte_swap) ByteSwap::byteSwap(isp_array,npbi);
              MPI_File_read_at(fh,offset+int_size*npbi,isp_p_array,npbi,MPI_INT,MPI_STATUS_IGNORE);
              if ( byte_swap) ByteSwap::byteSwap(isp_p_array,npbi);
              MPI_File_read_at(fh,offset+int_size*npbi*2,isp_bits_array,npbi,MPI_UINT2,MPI_STATUS_IGNORE);
              if ( byte_swap) ByteSwap::byteSwap(isp_bits_array,npbi);
            }
            offset += (int_size*2 + sizeof(uint2))*npbi;
            MPI_Bcast(isp_array,npbi,MPI_INT,0,mpi_comm_internode);
            MPI_Bcast(isp_p_array,npbi,MPI_INT,0,mpi_comm_internode);
            MPI_Bcast(isp_bits_array,npbi,MPI_UINT2,0,mpi_comm_internode);

            assert(pbi_split == NULL);
            pbi_split = new uint8[nsp_split];
            for (int isp = isp0; isp < isp1; ++isp)
              pbi_split[isp-isp0] = isp;

            for (int ipbi = 0; ipbi < npbi; ++ipbi) {
              const int isp = isp_array[ipbi];
              if ((isp >= isp0)&&(isp < isp1)) {
                const int isp_parent = isp_p_array[ipbi];
                const uint2 bits     = isp_bits_array[ipbi];
                assert(pbi_split[isp-isp0] == uint8(isp));
                pbi_split[isp-isp0] = BitUtils::packPbiHash(bits,isp_parent);
              }
            }

            delete[] isp_array;
            delete[] isp_p_array;
            delete[] isp_bits_array;
          }

          // either with or without periodic data...
	  // note: this does not work as part of stitch_s with NO_MPI on onyx...
          //if (!(offset == offset_check))
	  //  cout << "offset: " << offset << " offset_check: " << offset_check << endl;
          //assert(offset == offset_check);
          MPI_File_close(&fh);

        }

      }

    } // if (mpi_rank_shared == 0)...

    // error checking before going any further...
    int ierr;
    MPI_Allreduce(&my_ierr,&ierr,1,MPI_INT,MPI_MIN,mpi_comm);
    if (ierr != 0)
      return -1;

    // everyone needs the sizes in size_buf (currently set on only the shared rank 0...
    MPI_Bcast(size_buf,8,MPI_INT,0,mpi_comm_shared);
    if (mpi_rank_shared != 0) {
      nzones = size_buf[0];
      nsp = size_buf[1];
      nst = size_buf[4];
      npt = size_buf[6];
      npz = size_buf[7];
    }

    // if periodicity was found, set the periodicTransformVec...
    // this is a member of the partData namespace now, so if
    // it already exists, we will have to do some work reconciling
    // bits. For now, just assume/insist it does not exist, or that the
    // values and bits match...
    if (npt > 0) {
      MPI_Bcast(periodic_kind,npt,MPI_INT,0,mpi_comm_shared);
      MPI_Bcast(periodic_data,npt*3,MPI_DOUBLE,0,mpi_comm_shared);
      if (PeriodicData::periodicTransformVec.empty()) {
        PeriodicData::periodicTransformVec.resize(npt);
        for (int ipt = 0; ipt < npt; ++ipt) {
          PeriodicData::periodicTransformVec[ipt].setKindAndData(periodic_kind[ipt],periodic_data[ipt]);
        }
      }
      else {
        cout << "WARNING: Assuming existing periodicity matches" << endl;
      }
    }

    // zone names were read into the char_buf on rank 0 only.
    // bcast to everyone (ALL ranks)...
    if (mpi_rank != 0) {
      assert(char_buf == NULL);
      char_buf = new char[zone_name_len_max*nzones];
    }
    MPI_Bcast(char_buf,zone_name_len_max*nzones,MPI_CHAR,0,mpi_comm); // bcast to everyone
    if (mpi_rank != 0) {
      zoneVec.resize(nzones); // nzones
      for (int izone = 0; izone < nzones; ++izone) {
        zoneVec[izone].setName(char_buf+zone_name_len_max*izone);
      }
    }
    delete[] char_buf;

    // if some zones were periodic, then set this now...
    if (npz > 0) {
      if (mpi_rank_shared != 0) {
        assert(periodic_zone_array == NULL);
        periodic_zone_array = new int[npz];
        assert(periodic_zone_bits == NULL);
        periodic_zone_bits = new uint2[npz];
      }
      else {
        assert(periodic_zone_array != NULL);
        assert(periodic_zone_bits != NULL);
      }
      MPI_Bcast(periodic_zone_array,npz,MPI_INT,0,mpi_comm_shared);
      MPI_Bcast(periodic_zone_bits,npz,MPI_UINT2,0,mpi_comm_shared);
      for (int ipz = 0; ipz < npz; ++ipz) {
        const int izone = periodic_zone_array[ipz]; assert((izone >= 0)&&(izone < zoneVec.size()));
        zoneVec[izone].setPeriodicBits(periodic_zone_bits[ipz]);
      }
      delete[] periodic_zone_array;
      delete[] periodic_zone_bits;
    }

    // allocate the shared memory for xsp AND spost AND znost

    //if (mpi_rank == 0) cout << " > allgather on xsp..." << endl;
    nsp_max = size_buf[2]*size_buf[3]; assert(nsp_max >= nsp);
    assert(xsp == NULL);
    CTI_Mmap_rw(xsp,nsp_max); // _rw
    if (mpi_rank_shared == 0) {
      MPI_Allgather((double*)xsp_split,size_buf[2]*3,MPI_DOUBLE,
                    (double*)xsp,size_buf[2]*3,MPI_DOUBLE,mpi_comm_internode);
      delete[] xsp_split;
    }

    if (npt > 0) {
      //if (mpi_rank == 0) cout << " > allgather on pbi..." << endl;
      assert(pbi == NULL);
      CTI_Mmap(pbi,nsp_max);
      if (mpi_rank_shared == 0) {
        assert(pbi_split != NULL);
        MPI_Allgather(pbi_split,size_buf[2],MPI_UINT8,
                      pbi,size_buf[2],MPI_UINT8,mpi_comm_internode);
        delete[] pbi_split;
      }
    }

    //if (mpi_rank == 0) cout << " > allgather on spost and znost..." << endl;
    nst_max = size_buf[5]*size_buf[3]; assert(nst_max >= nst);
    assert(spost == NULL);
    CTI_Mmap(spost,nst_max);
    assert(znost == NULL);
    CTI_Mmap(znost,nst_max);
    if (mpi_rank_shared == 0) {
      MPI_Allgather((int*)spost_split,size_buf[5]*3,MPI_INT,
                    (int*)spost,size_buf[5]*3,MPI_INT,mpi_comm_internode);
      delete[] spost_split;
      MPI_Allgather(znost_split,size_buf[5],MPI_INT,
                    znost,size_buf[5],MPI_INT,mpi_comm_internode);
      delete[] znost_split;
    }

    MPI_Barrier(mpi_comm_shared);

    // check pbi...
    // should do this in pieces...
    if (pbi) {
      double my_d2_max = 0.0;
      int count = 0;
      for (int isp = 0; isp < nsp; ++isp) {
        if (pbi[isp] != uint8(isp)) {
          ++count;
          int bits,isp_parent;
          BitUtils::unpackPbiHash(bits,isp_parent,pbi[isp]);
          assert((isp_parent >= 0)&&(isp_parent < nsp));
          double xsp_t[3];
          FOR_I3 xsp_t[i] = xsp[isp_parent][i];
          PeriodicData::periodicTranslate(xsp_t,1,bits);
          const double d2 = DIST2(xsp_t,xsp[isp]);
          my_d2_max = max(my_d2_max,d2);
        }
      }
      double d2_max;
      MPI_Reduce(&my_d2_max,&d2_max,1,MPI_DOUBLE,MPI_MAX,0,mpi_comm);
      if (mpi_rank == 0)
        cout << " > check pbi (should be zero): " << sqrt(d2_max) << " for " << count << " child nodes." << endl;
    }

    if (mpi_rank == 0) cout << " > done" << endl;
    return 0;

  }

  int writeSbin(const string& filename) { return writeBinary(filename); }
  
  int writeBinary(const string& filename) {

    if (mpi_rank == 0) {
      cout << "SurfaceShm::writeBinary \"" << filename << "\"..." << endl;

      MiscUtils::mkdir_for_file(filename);
      string tmp_filename = MiscUtils::makeTmpPrefix(filename);
      FILE * fp = fopen(tmp_filename.c_str(),"wb");
      assert(fp != NULL);

      const int version = 1;
      fwrite(&version,sizeof(int),1,fp);

      const int count = zoneVec.size();
      fwrite(&count,sizeof(int),1,fp);
      cout << " > writing " << count << " zones" << endl;
      for (int izone = 0, limit = zoneVec.size(); izone < limit; ++izone) {
        const int length = zoneVec[izone].getName().length();
        fwrite(&length,sizeof(int),1,fp);
        fwrite(zoneVec[izone].getName().c_str(),sizeof(char),length,fp);
        cout << " > zone " << zoneVec[izone].getName() << endl;
      }

      cout << " > nsp = " << nsp << endl;
      fwrite(&nsp,sizeof(int),1,fp);
      fwrite(xsp,sizeof(double),nsp*3,fp);

      cout << " > nst = " << nst << endl;
      fwrite(&nst,sizeof(int),1,fp);
      fwrite(spost,sizeof(int),nst*3,fp);
      fwrite(znost,sizeof(int),nst,fp);

      // include periodic data if present...
      if (pbi) {
        const int npt = PeriodicData::periodicTransformVec.size();
        cout << " > periodic transforms: " << npt << endl;
        fwrite(&npt,sizeof(int),1,fp);
        for (int ipt = 0; ipt < npt; ++ipt)
          PeriodicData::periodicTransformVec[ipt].writeBinary(fp);
        // write periodic bits for periodic zones only...
        int npz = 0;
        for (int izone = 0, limit = zoneVec.size(); izone < limit; ++izone)
          if (zoneVec[izone].getPeriodicBits())
            ++npz;
        cout << " > npz: " << npz << endl;
        fwrite(&npz,sizeof(int),1,fp);
        int * periodic_zone_array = new int[npz];
        uint2 * periodic_zone_bits = new uint2[npz];
        int ipz = 0;
        for (int izone = 0, limit = zoneVec.size(); izone < limit; ++izone) {
          if (zoneVec[izone].getPeriodicBits()) {
            periodic_zone_array[ipz] = izone;
            periodic_zone_bits[ipz]  = zoneVec[izone].getPeriodicBits();
            ++ipz;
          }
        }
        assert(ipz == npz);
        fwrite(periodic_zone_array,sizeof(int),npz,fp);
        fwrite(periodic_zone_bits,sizeof(uint2),npz,fp);
        delete[] periodic_zone_array;
        delete[] periodic_zone_bits;
        // only write the pbi's that contain transform data...
        int npbi = 0;
        for (int isp = 0; isp < nsp; ++isp)
          if (pbi[isp] != uint8(isp))
            ++npbi;
        cout << " > npbi: " << npbi << endl;
        fwrite(&npbi,sizeof(int),1,fp);
        int * isp_array = new int[npbi];
        int * isp_p_array = new int[npbi];
        uint2 * isp_bits_array = new uint2[npbi];
        int ipbi = 0;
        for (int isp = 0; isp < nsp; ++isp) {
          if (pbi[isp] != uint8(isp)) {
            isp_array[ipbi]      = isp;
            isp_p_array[ipbi]    = int(pbi[isp]&MASK_52BITS);
            isp_bits_array[ipbi] = uint2(pbi[isp]>>52);
            ++ipbi;
          }
        }
        assert(ipbi == npbi);
        fwrite(isp_array,sizeof(int),npbi,fp);
        fwrite(isp_p_array,sizeof(int),npbi,fp);
        fwrite(isp_bits_array,sizeof(uint2),npbi,fp);
        delete[] isp_array;
        delete[] isp_p_array;
        delete[] isp_bits_array;
      }

      fclose(fp);
      remove(filename.c_str());
      rename(tmp_filename.c_str(),filename.c_str());

    }

    return 0;
  }


  void calcZoneBbminmax(double zone_bbminmax[6], const string &zoneNames) {

    if (mpi_rank == 0) cout << "SurfaceShm::calcZoneBbminmax: " << zoneNames << endl;
    double my_bbminmax[6] = { HUGE_VAL, HUGE_VAL, HUGE_VAL, HUGE_VAL, HUGE_VAL, HUGE_VAL };

    //use striped distribution of surface tris, include only points associated
    //with tris in zoneNames.
    int my_nst_avg = nst/mpi_size;
    if (nst%mpi_size) ++my_nst_avg;
    const int ist0 = min(nst,mpi_rank*my_nst_avg);
    const int ist1 = min(nst,(mpi_rank+1)*my_nst_avg);
    assert(ist1-ist0 <= my_nst_avg);

    set<string> zoneNameSet;
    MiscUtils::splitCsv(zoneNameSet,zoneNames);
    for (int izn = 0; izn < zoneVec.size(); ++izn) {
      if (zoneNameSet.find(zoneVec[izn].getName()) != zoneNameSet.end()) {
        for (int ist = ist0; ist < ist1; ++ist) {
          if (znost[ist] == izn) {
            FOR_I3 {
              const int isp = spost[ist][i];
              FOR_I3 my_bbminmax[i] = min(my_bbminmax[i],xsp[isp][i]);
              FOR_I3 my_bbminmax[i+3] = min(my_bbminmax[i+3],-xsp[isp][i]);
            }
          }
        }
      }
    }

    MPI_Allreduce(my_bbminmax,zone_bbminmax,6,MPI_DOUBLE,MPI_MIN,mpi_comm);
    if (mpi_rank == 0) cout << "  Bounding Box: " <<
                         zone_bbminmax[0] << ":" << -zone_bbminmax[3] << " " <<
                         zone_bbminmax[1] << ":" << -zone_bbminmax[4] << " " <<
                         zone_bbminmax[2] << ":" << -zone_bbminmax[5] << endl;

  }

  void clearGeom() {
    b_bbminmax = b_volume_and_xc = b_mom = b_area = false;
  }

  void calcBbminmax() {
    if (mpi_rank == 0) cout << "SurfaceShm::calcBbminmax: ";
    assert(!b_bbminmax);
    // to compute the bbox fast, use a striped ditribution of surface points...
    double my_bbminmax[6] = { HUGE_VAL, HUGE_VAL, HUGE_VAL, HUGE_VAL, HUGE_VAL, HUGE_VAL };
    int my_nsp_avg = nsp/mpi_size;
    if (nsp%mpi_size) ++my_nsp_avg;
    const int isp0 = min(nsp,mpi_rank*my_nsp_avg);
    const int isp1 = min(nsp,(mpi_rank+1)*my_nsp_avg);
    assert(isp1-isp0 <= my_nsp_avg);
    for (int isp = isp0; isp < isp1; ++isp) {
      FOR_I3 my_bbminmax[i] = min(my_bbminmax[i],xsp[isp][i]);
      FOR_I3 my_bbminmax[i+3] = min(my_bbminmax[i+3],-xsp[isp][i]);
    }
    MPI_Allreduce(my_bbminmax,bbminmax,6,MPI_DOUBLE,MPI_MIN,mpi_comm);
    b_bbminmax = true;
    if (mpi_rank == 0) cout << bbminmax[0] << ":" << -bbminmax[3] << " " <<
                         bbminmax[1] << ":" << -bbminmax[4] << " " <<
                         bbminmax[2] << ":" << -bbminmax[5] << endl;
  }

  void getBbminmax(double bbminmax[6]) {
    if (!b_bbminmax) calcBbminmax();
    FOR_I6 bbminmax[i] = this->bbminmax[i];
  }

  void getBbox(double bbmin[3],double bbmax[3]) {
    if (!b_bbminmax) calcBbminmax();
    FOR_I3 bbmin[i] = bbminmax[i];
    FOR_I3 bbmax[i] = -bbminmax[i+3];
  }

  void getBboxForZoneids(double bbminmax[6],const vector<int>& zoneidVec) {

    std::vector<int> sz_flag(zoneVec.size());
    for (int isz = 0; isz < zoneVec.size(); ++isz)
      sz_flag[isz] = 0;
    for (int ii = 0; ii < zoneidVec.size(); ++ii) {
      const int isz = zoneidVec[ii]; assert((isz >= 0)&&(isz < zoneVec.size()));
      sz_flag[isz] = 1;
    }

    double my_bbminmax[6] = { HUGE_VAL, HUGE_VAL, HUGE_VAL, HUGE_VAL, HUGE_VAL, HUGE_VAL };
    int my_nst_avg = nst/mpi_size;
    if (nst%mpi_size) ++my_nst_avg;

    const int ist0 = min(nst,mpi_rank*my_nst_avg);
    const int ist1 = min(nst,(mpi_rank+1)*my_nst_avg);
    assert(ist1-ist0 <= my_nst_avg);

    for (int ist = ist0; ist < ist1; ++ist) {
      const int isz = znost[ist]; assert((isz >= 0)&&(isz < zoneVec.size()));
      if (sz_flag[isz] == 1) {
        // this tri contributes to the bbox of the surface zone...
        FOR_I3 {
          const int isp = spost[ist][i];
          FOR_J3 my_bbminmax[j] = min(my_bbminmax[j],xsp[isp][j]);
          FOR_J3 my_bbminmax[j+3] = min(my_bbminmax[j+3],-xsp[isp][j]);
        }
      }
    }
    MPI_Allreduce(my_bbminmax,bbminmax,6,MPI_DOUBLE,MPI_MIN,mpi_comm);

  }

  void getCentroid(double xc_[3]) {
    if (!b_volume_and_xc) calcVolumeAndCentroid();
    assert(b_volume_and_xc);
    FOR_I3 xc_[i] = xc[i];
  }

  double getVolume() {
    if(!b_volume_and_xc) calcVolumeAndCentroid();
    assert(b_volume_and_xc);
    return volume;
  }

  void calcVolumeAndCentroid() {
    if (mpi_rank == 0) cout << "SurfaceShm::calcVolumeAndCentroid: ";
    assert(!b_volume_and_xc);
    int nst_split = nst/mpi_size;
    if (nst%mpi_size != 0) nst_split += 1;
    const int ist0 = min(nst,nst_split*mpi_rank);
    const int ist1 = min(nst,nst_split*(mpi_rank+1));
    double my_buf[4] = { 0.0, 0.0, 0.0, 0.0 };
    // all tet's use the first point to compute volume from...
    const double * const xsp0 = xsp[0];
    for (int ist = ist0; ist < ist1; ++ist) {
      const double * const xsp1 = xsp[spost[ist][0]];
      const double * const xsp2 = xsp[spost[ist][1]];
      const double * const xsp3 = xsp[spost[ist][2]];
      const double this_vol6 = SIGNED_TET_VOLUME_6(xsp0,xsp1,xsp2,xsp3);
      FOR_I3 my_buf[i] += this_vol6*(xsp0[i]+xsp1[i]+xsp2[i]+xsp3[i]);
      my_buf[3] += this_vol6;
    }
    double buf[4];
    MPI_Allreduce(my_buf,buf,4,MPI_DOUBLE,MPI_SUM,mpi_comm);
    assert(buf[3] != 0.0);
    volume=-buf[3]/6.0;
    FOR_I3 xc[i] = 0.25*buf[i]/buf[3];
    // NOTE: the volume is computed as the volume of a solid (i.e. silver
    // out, gold in) mainly because of the interpretation of solids in 
    // the moving solver. This is not consistent with the volume inside and 
    // the isInside checks, so report -volume here, and revisit when 
    // next using the moving solver...
    // TODO: change volume to be positive when silver points in, gold out...
    if (mpi_rank == 0) cout << -volume << " " << COUT_VEC(xc) << endl;
    b_volume_and_xc = true;
  }

  double getArea() {
    if(!b_area) calcArea();
    assert(b_area);
    return area;
  }

  void calcArea() {
    if (mpi_rank == 0) cout << "SurfaceShm::calcArea: ";
    assert(!b_area);
    int nst_split = nst/mpi_size;
    if (nst%mpi_size != 0) nst_split += 1;
    const int ist0 = min(nst,nst_split*mpi_rank);
    const int ist1 = min(nst,nst_split*(mpi_rank+1));
    double my_area = 0.0;
    // all tet's use the first point to compute volume from...
    for (int ist = ist0; ist < ist1; ++ist) {
      const double * const x0 = xsp[spost[ist][0]];
      const double * const x1 = xsp[spost[ist][1]];
      const double * const x2 = xsp[spost[ist][2]];
      const double n_st[3] = TRI_NORMAL_2(x0,x1,x2);
      my_area += MAG(n_st);
    }
    MPI_Allreduce(&my_area,&area,1,MPI_DOUBLE,MPI_SUM,mpi_comm);
    area *= 0.5;
    if (mpi_rank == 0) cout << area << endl;
    b_area = true;
  }

  void translate(const double dx[3]) {
    if (mpi_rank == 0) cout << "SurfaceShm::translate: " << COUT_VEC(dx) << endl;
    if (mpi_rank_shared == 0) {
      for (int isp = 0; isp < nsp; ++isp) {
        FOR_I3 xsp[isp][i] += dx[i];
      }
    }
    clearGeom();
    MPI_Barrier(mpi_comm_shared);
  }

  void rotate(const double point[3], const double _axis[3], const double angle_deg) {
    if (mpi_rank_shared == 0) {
      double axis[3] = {_axis[0],_axis[1],_axis[2]};
      NORMALIZE(axis);
      const double angle = angle_deg/180.0*M_PI;  // radians

      COUT1("SurfaceShm::rotate()");
      COUT1(" > axis point: " << point[0] << " " << point[1] << " " << point[2]);
      COUT1(" > axis unit vector: " << axis[0] << " " << axis[1] << " " << axis[2]);
      COUT1(" > theta: " << angle_deg << " (degrees)");

      const double cp_mat[3][3] = {
        {0.0,-axis[2],axis[1]},
        {axis[2],0.0,-axis[0]},
        {-axis[1],axis[0],0.0}
      };

      // create rotation matrix
      // R = cos0*I + sin0*[axis]_x + (1-cos0)*(u tensor-product u)
      double R[3][3];
      FOR_I3 {
        R[i][i] = cos(angle);
        FOR_J3 {
          if (i != j) R[i][j] = sin(angle)*cp_mat[i][j];
          R[i][j] += (1-cos(angle))*axis[i]*axis[j];
        }
      }

      FOR_ISP {
        double _xsp[3];
        FOR_J3 _xsp[j] = xsp[isp][j] - point[j];
        const double temp[3] = MATRIX_PRODUCT(R,_xsp);
        FOR_J3 xsp[isp][j] = temp[j] + point[j];
      }
    }
    clearGeom();
    MPI_Barrier(mpi_comm_shared);
  }

  void scale(const double dx[3]) {
    if (mpi_rank == 0) cout << "SurfaceShm::scale: " << COUT_VEC(dx) << endl;
    if (mpi_rank_shared == 0) {
      for (int isp = 0; isp < nsp; ++isp) {
        FOR_I3 xsp[isp][i] *= dx[i];
      }
    }
    clearGeom();
    //HERE
    MPI_Barrier(mpi_comm_shared);
  }

  void getMoments(double diag_[3],double offd_[3]) {
    // notes: corrdinates are wrt xc...
    // diag[0]: int(x*x,dV);
    // diag[1]: int(y*y,dV);
    // diag[2]: int(z*z,dV);
    // offd[0]: int(y*z,dV);
    // offd[1]: int(x*z,dV);
    // offd[2]: int(x*y,dV);
    if (!b_mom) calcMoments();
    FOR_I3 diag_[i] = mom_diag[i];
    FOR_I3 offd_[i] = mom_offd[i];
  }

  void calcMoments() {
    using GaussQuadrature::tet4;
    if (!b_volume_and_xc) calcVolumeAndCentroid();
    if (mpi_rank == 0) cout << "SurfaceShm::calcMoments: ";
    int nst_split = nst/mpi_size;
    if (nst%mpi_size != 0) nst_split += 1;
    const int ist0 = min(nst,nst_split*mpi_rank);
    const int ist1 = min(nst,nst_split*(mpi_rank+1));
    double my_buf[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    // all tet's use the centroid, since it is available...
    double xtet[4][3];
    FOR_I3 xtet[0][i] = 0.0;
    for (int ist = ist0; ist < ist1; ++ist) {
      FOR_I3 xtet[1][i] = xsp[spost[ist][0]][i] - xc[i];
      FOR_I3 xtet[2][i] = xsp[spost[ist][1]][i] - xc[i];
      FOR_I3 xtet[3][i] = xsp[spost[ist][2]][i] - xc[i];
      const double this_vol6 = SIGNED_TET_VOLUME_6(xtet[0],xtet[1],xtet[2],xtet[3]);
      // quadrature loop...
      for (int ip = 0; ip < 4; ++ip) {
        double x[3] = { 0.0, 0.0, 0.0 };
        for (int i = 0; i < 4; ++i) FOR_J3 x[j] += tet4[ip][i]*xtet[i][j];
        // recall wgts that sum to 1 are stored in tet4[ip][4]...
        my_buf[0] += this_vol6*tet4[ip][4]*x[0]*x[0];
        my_buf[1] += this_vol6*tet4[ip][4]*x[1]*x[1];
        my_buf[2] += this_vol6*tet4[ip][4]*x[2]*x[2];
        my_buf[3] += this_vol6*tet4[ip][4]*x[1]*x[2];
        my_buf[4] += this_vol6*tet4[ip][4]*x[0]*x[2];
        my_buf[5] += this_vol6*tet4[ip][4]*x[0]*x[1];
        //my_buf[6] += this_vol6;
      }
    }
    double buf[6];
    MPI_Allreduce(my_buf,buf,6,MPI_DOUBLE,MPI_SUM,mpi_comm);
    // here we account for the handedness to produce +ve moments...
    FOR_I3 mom_diag[i] = -buf[i]/6.0; // recall vol6 used above
    FOR_I3 mom_offd[i] = -buf[i+3]/6.0; // recall vol6 used above
    if (mpi_rank == 0) cout << "diag: " << COUT_VEC(mom_diag) << ", offd: " << COUT_VEC(mom_offd) << endl;
    b_mom = true;
  }

  void writeTecplot(const string& filename,double * data = NULL) const {

    if (mpi_rank == 0) {

      COUT1("SurfaceShm::writeTecplot(" << filename << ")");

      FILE * fp = fopen(filename.c_str(),"w");
      fprintf(fp,"TITLE = \"%s\"\n",filename.c_str());
      fprintf(fp,"VARIABLES = \"X\"\n");
      fprintf(fp,"\"Y\"\n");
      fprintf(fp,"\"Z\"\n");
      if (data) fprintf(fp,"\"DATA\"\n");

      int * sp_flag = new int[nsp];
      for (int izone = 0; izone < zoneVec.size(); ++izone) {

        for (int isp = 0; isp < nsp; ++isp)
          sp_flag[isp] = -1;

        int nt = 0;
        for (int ist = 0; ist < nst; ++ist) {
          if (znost[ist] == izone) {
            ++nt;
            FOR_I3 {
              const int isp = spost[ist][i];
              sp_flag[isp] = 0;
            }
          }
        }

        int np = 0;
        for (int isp = 0; isp < nsp; ++isp)
          if (sp_flag[isp] == 0)
            sp_flag[isp] = np++;

        if (nt > 0) {

          fprintf(fp,"ZONE T=\"%s\"\n",zoneVec[izone].getName().c_str());
          fprintf(fp,"N=%d, E=%d, F=FEPOINT, ET=TRIANGLE\n",np,nt);

          // data...
          if (data) {
            for (int isp = 0; isp < nsp; ++isp) {
              if (sp_flag[isp] >= 0) {
                fprintf(fp,"%18.15e %18.15e %18.15e %18.15e\n",xsp[isp][0],xsp[isp][1],xsp[isp][2],data[isp]);
              }
            }
          }
          else {
            for (int isp = 0; isp < nsp; ++isp) {
              if (sp_flag[isp] >= 0) {
                fprintf(fp,"%18.15e %18.15e %18.15e\n",xsp[isp][0],xsp[isp][1],xsp[isp][2]);
              }
            }
          }

          for (int ist = 0; ist < nst; ++ist) {
            if (znost[ist] == izone) {
              fprintf(fp,"%d %d %d\n",
                      sp_flag[spost[ist][0]]+1,
                      sp_flag[spost[ist][1]]+1,
                      sp_flag[spost[ist][2]]+1);
            }
          }

        }

      }

      fclose(fp);

      delete[] sp_flag;

    }

  }

  void writeTecplotNoZones(const string& filename) const {

    if (mpi_rank == 0) {

      COUT1("SurfaceShm::writeTecplot(" << filename << ")");

      FILE * fp = fopen(filename.c_str(),"w");
      fprintf(fp,"TITLE = \"%s\"\n",filename.c_str());
      fprintf(fp,"VARIABLES = \"X\"\n");
      fprintf(fp,"\"Y\"\n");
      fprintf(fp,"\"Z\"\n");
      fprintf(fp,"ZONE T=\"%s\"\n",filename.c_str());
      fprintf(fp,"N=%d, E=%d, F=FEPOINT, ET=TRIANGLE\n",nsp,nst);

      // data...
      for (int isp = 0; isp < nsp; ++isp) {
        fprintf(fp,"%18.15e %18.15e %18.15e\n",xsp[isp][0],xsp[isp][1],xsp[isp][2]);
      }

      for (int ist = 0; ist < nst; ++ist) {
        fprintf(fp,"%d %d %d\n",
                spost[ist][0]+1,
                spost[ist][1]+1,
                spost[ist][2]+1);
      }

      fclose(fp);

    }

  }

  void makeBox(const double x0[3],const double x1[3],const bool periodic[3]) {

    // x0,x1 are WSB and ENT points respectively.

    //
    //        7------6
    //       /      /|
    //      4------5 |
    //      | 3    | 2
    //      |      |/
    //      0------1
    //
    // x is 0->1, y is 0->3, z is 0->4

    assert(zoneVec.empty());
    assert(xsp == NULL);
    assert(spost == NULL);
    assert(znost == NULL);
    assert(pbi == NULL);

    // everyone gets the zones...

    zoneVec.push_back(StZone("x0"));
    zoneVec.push_back(StZone("x1"));
    zoneVec.push_back(StZone("y0"));
    zoneVec.push_back(StZone("y1"));
    zoneVec.push_back(StZone("z0"));
    zoneVec.push_back(StZone("z1"));

    // shm allocation...

    init(8,12);

    if (mpi_rank_shared == 0) {

      xsp[0][0] = x0[0]; xsp[0][1] = x0[1]; xsp[0][2] = x0[2];
      xsp[1][0] = x1[0]; xsp[1][1] = x0[1]; xsp[1][2] = x0[2];
      xsp[2][0] = x1[0]; xsp[2][1] = x1[1]; xsp[2][2] = x0[2];
      xsp[3][0] = x0[0]; xsp[3][1] = x1[1]; xsp[3][2] = x0[2];
      xsp[4][0] = x0[0]; xsp[4][1] = x0[1]; xsp[4][2] = x1[2];
      xsp[5][0] = x1[0]; xsp[5][1] = x0[1]; xsp[5][2] = x1[2];
      xsp[6][0] = x1[0]; xsp[6][1] = x1[1]; xsp[6][2] = x1[2];
      xsp[7][0] = x0[0]; xsp[7][1] = x1[1]; xsp[7][2] = x1[2];

      const int zone_x0 = 0;
      const int zone_x1 = 1;
      const int zone_y0 = 2;
      const int zone_y1 = 3;
      const int zone_z0 = 4;
      const int zone_z1 = 5;

      // x0...
      spost[0][0] = 0; spost[0][1] = 4; spost[0][2] = 7; znost[0] = zone_x0;
      spost[1][0] = 0; spost[1][1] = 7; spost[1][2] = 3; znost[1] = zone_x0;

      // x1...
      spost[2][0] = 1; spost[2][1] = 6; spost[2][2] = 5; znost[2] = zone_x1;
      spost[3][0] = 1; spost[3][1] = 2; spost[3][2] = 6; znost[3] = zone_x1;

      // y0...
      spost[4][0] = 0; spost[4][1] = 1; spost[4][2] = 5; znost[4] = zone_y0;
      spost[5][0] = 0; spost[5][1] = 5; spost[5][2] = 4; znost[5] = zone_y0;

      // y1...
      spost[6][0] = 2; spost[6][1] = 3; spost[6][2] = 6; znost[6] = zone_y1;
      spost[7][0] = 3; spost[7][1] = 7; spost[7][2] = 6; znost[7] = zone_y1;

      // z0...
      spost[8][0] = 0; spost[8][1] = 3; spost[8][2] = 2; znost[8] = zone_z0;
      spost[9][0] = 0; spost[9][1] = 2; spost[9][2] = 1; znost[9] = zone_z0;

      // z1...
      spost[10][0] = 4; spost[10][1] = 5; spost[10][2] = 6; znost[10] = zone_z1;
      spost[11][0] = 4; spost[11][1] = 6; spost[11][2] = 7; znost[11] = zone_z1;

    }
    MPI_Barrier(mpi_comm_shared);

    if ( periodic[0] && !periodic[1] && !periodic[2] ) {
      if (mpi_rank == 0) cout << " > setting x periodicity..." << endl;
      assert(PeriodicData::periodicTransformVec.empty());
      PeriodicData::periodicTransformVec.push_back(PeriodicTransform(PERIODIC_TRANSFORM_CART,x1[0]-x0[0],0.0,0.0));
      CTI_Mmap(pbi,8);
      if (mpi_rank_shared == 0) {
        pbi[0] = BitUtils::packPbiHash(0,0);
        pbi[3] = BitUtils::packPbiHash(0,3);
        pbi[4] = BitUtils::packPbiHash(0,4);
        pbi[7] = BitUtils::packPbiHash(0,7);
        pbi[1] = BitUtils::packPbiHash((1<<0),0);
        pbi[2] = BitUtils::packPbiHash((1<<0),3);
        pbi[5] = BitUtils::packPbiHash((1<<0),4);
        pbi[6] = BitUtils::packPbiHash((1<<0),7);
      }
      zoneVec[0].setPeriodicBits(1<<1); // x0
      zoneVec[1].setPeriodicBits(1<<0); // x1
    }
    else if ( !periodic[0] && periodic[1] && !periodic[2] ) {
      if (mpi_rank == 0) cout << " > setting y periodicity..." << endl;
      assert(PeriodicData::periodicTransformVec.empty());
      PeriodicData::periodicTransformVec.push_back(PeriodicTransform(PERIODIC_TRANSFORM_CART,0.0,x1[1]-x0[1],0.0));
      CTI_Mmap(pbi,8);
      if (mpi_rank_shared == 0) {
        pbi[0] = BitUtils::packPbiHash(0,0);
        pbi[4] = BitUtils::packPbiHash(0,4);
        pbi[5] = BitUtils::packPbiHash(0,5);
        pbi[1] = BitUtils::packPbiHash(0,1);
        pbi[3] = BitUtils::packPbiHash((1<<0),0);
        pbi[7] = BitUtils::packPbiHash((1<<0),4);
        pbi[6] = BitUtils::packPbiHash((1<<0),5);
        pbi[2] = BitUtils::packPbiHash((1<<0),1);
      }
      zoneVec[2].setPeriodicBits(1<<1); // y0
      zoneVec[3].setPeriodicBits(1<<0); // y1
    }
    else if ( !periodic[0] && !periodic[1] && periodic[2] ) {
      if (mpi_rank == 0) cout << " > setting z periodicity..." << endl;
      assert(PeriodicData::periodicTransformVec.empty());
      PeriodicData::periodicTransformVec.push_back(PeriodicTransform(PERIODIC_TRANSFORM_CART,0.0,0.0,x1[2]-x0[2]));
      CTI_Mmap(pbi,8);
      if (mpi_rank_shared == 0) {
        pbi[0] = BitUtils::packPbiHash(0,0);
        pbi[1] = BitUtils::packPbiHash(0,1);
        pbi[2] = BitUtils::packPbiHash(0,2);
        pbi[3] = BitUtils::packPbiHash(0,3);
        pbi[4] = BitUtils::packPbiHash((1<<0),0);
        pbi[5] = BitUtils::packPbiHash((1<<0),1);
        pbi[6] = BitUtils::packPbiHash((1<<0),2);
        pbi[7] = BitUtils::packPbiHash((1<<0),3);
      }
      zoneVec[4].setPeriodicBits(1<<1); // z0
      zoneVec[5].setPeriodicBits(1<<0); // z1
    }
    else if ( periodic[0] && periodic[1] && !periodic[2] ) {
      if (mpi_rank == 0) cout << " > setting xy periodicity..." << endl;
      assert(PeriodicData::periodicTransformVec.empty());
      PeriodicData::periodicTransformVec.push_back(PeriodicTransform(PERIODIC_TRANSFORM_CART,x1[0]-x0[0],0.0,0.0));
      PeriodicData::periodicTransformVec.push_back(PeriodicTransform(PERIODIC_TRANSFORM_CART,0.0,x1[1]-x0[1],0.0));
      CTI_Mmap(pbi,8);
      if (mpi_rank_shared == 0) {
        pbi[0] = BitUtils::packPbiHash(0,0);
        pbi[4] = BitUtils::packPbiHash(0,4);
        pbi[1] = BitUtils::packPbiHash((1<<0),0);
        pbi[5] = BitUtils::packPbiHash((1<<0),4);
        pbi[3] = BitUtils::packPbiHash((1<<2),0); // 2,3 is the second bit pair
        pbi[7] = BitUtils::packPbiHash((1<<2),4);
        pbi[2] = BitUtils::packPbiHash((1<<0)|(1<<2),0);
        pbi[6] = BitUtils::packPbiHash((1<<0)|(1<<2),4);
      }
      zoneVec[0].setPeriodicBits(1<<1); // x0
      zoneVec[1].setPeriodicBits(1<<0); // x1
      zoneVec[2].setPeriodicBits(1<<3); // y0
      zoneVec[3].setPeriodicBits(1<<2); // y1
    }
    else if ( periodic[0] && !periodic[1] && periodic[2] ) {
      if (mpi_rank == 0) cout << " > setting xz periodicity..." << endl;
      assert(PeriodicData::periodicTransformVec.empty());
      PeriodicData::periodicTransformVec.push_back(PeriodicTransform(PERIODIC_TRANSFORM_CART,x1[0]-x0[0],0.0,0.0));
      PeriodicData::periodicTransformVec.push_back(PeriodicTransform(PERIODIC_TRANSFORM_CART,0.0,0.0,x1[2]-x0[2]));
      CTI_Mmap(pbi,8);
      if (mpi_rank_shared == 0) {
        pbi[0] = BitUtils::packPbiHash(0,0);
        pbi[3] = BitUtils::packPbiHash(0,3);
        pbi[1] = BitUtils::packPbiHash((1<<0),0);
        pbi[2] = BitUtils::packPbiHash((1<<0),3);
        pbi[4] = BitUtils::packPbiHash((1<<2),0); // 2,3 is the second bit pair
        pbi[7] = BitUtils::packPbiHash((1<<2),3);
        pbi[5] = BitUtils::packPbiHash((1<<0)|(1<<2),0);
        pbi[6] = BitUtils::packPbiHash((1<<0)|(1<<2),3);
      }
      zoneVec[0].setPeriodicBits(1<<1); // x0
      zoneVec[1].setPeriodicBits(1<<0); // x1
      zoneVec[4].setPeriodicBits(1<<3); // z0
      zoneVec[5].setPeriodicBits(1<<2); // z1
    }
    else if ( !periodic[0] && periodic[1] && periodic[2] ) {
      if (mpi_rank == 0) cout << " > setting yz periodicity..." << endl;
      assert(PeriodicData::periodicTransformVec.empty());
      PeriodicData::periodicTransformVec.push_back(PeriodicTransform(PERIODIC_TRANSFORM_CART,0.0,x1[1]-x0[1],0.0));
      PeriodicData::periodicTransformVec.push_back(PeriodicTransform(PERIODIC_TRANSFORM_CART,0.0,0.0,x1[2]-x0[2]));
      CTI_Mmap(pbi,8);
      if (mpi_rank_shared == 0) {
        pbi[0] = BitUtils::packPbiHash(0,0);
        pbi[1] = BitUtils::packPbiHash(0,1);
        pbi[3] = BitUtils::packPbiHash((1<<0),0);
        pbi[2] = BitUtils::packPbiHash((1<<0),1);
        pbi[4] = BitUtils::packPbiHash((1<<2),0); // 2,3 is the second bit pair
        pbi[5] = BitUtils::packPbiHash((1<<2),1);
        pbi[7] = BitUtils::packPbiHash((1<<0)|(1<<2),0);
        pbi[6] = BitUtils::packPbiHash((1<<0)|(1<<2),1);
      }
      zoneVec[2].setPeriodicBits(1<<1); // y0
      zoneVec[3].setPeriodicBits(1<<0); // y1
      zoneVec[4].setPeriodicBits(1<<3); // z0
      zoneVec[5].setPeriodicBits(1<<2); // z1
    }
    else if ( periodic[0] && periodic[1] && periodic[2] ) {
      if (mpi_rank == 0) cout << " > setting xyz periodicity..." << endl;
      assert(PeriodicData::periodicTransformVec.empty());
      PeriodicData::periodicTransformVec.push_back(PeriodicTransform(PERIODIC_TRANSFORM_CART,x1[0]-x0[0],0.0,0.0));
      PeriodicData::periodicTransformVec.push_back(PeriodicTransform(PERIODIC_TRANSFORM_CART,0.0,x1[1]-x0[1],0.0));
      PeriodicData::periodicTransformVec.push_back(PeriodicTransform(PERIODIC_TRANSFORM_CART,0.0,0.0,x1[2]-x0[2]));
      CTI_Mmap(pbi,8);
      if (mpi_rank_shared == 0) {
        pbi[0] = BitUtils::packPbiHash(0,0);
        pbi[1] = BitUtils::packPbiHash((1<<0),0);
        pbi[3] = BitUtils::packPbiHash((1<<2),0);
        pbi[4] = BitUtils::packPbiHash((1<<4),0);
        pbi[2] = BitUtils::packPbiHash((1<<0)|(1<<2),0); // 2,3 is the second bit pair
        pbi[5] = BitUtils::packPbiHash((1<<0)|(1<<4),0); // 4,5 is the third bit pair
        pbi[7] = BitUtils::packPbiHash((1<<2)|(1<<4),0);
        pbi[6] = BitUtils::packPbiHash((1<<0)|(1<<2)|(1<<4),0);
      }
      zoneVec[0].setPeriodicBits(1<<1); // x0
      zoneVec[1].setPeriodicBits(1<<0); // x1
      zoneVec[2].setPeriodicBits(1<<3); // y0
      zoneVec[3].setPeriodicBits(1<<2); // y1
      zoneVec[4].setPeriodicBits(1<<5); // z0
      zoneVec[5].setPeriodicBits(1<<4); // z1
    }
    else {
      assert(!periodic[0]);
      assert(!periodic[1]);
      assert(!periodic[2]);
    }

    MPI_Barrier(mpi_comm_shared);

  }

  void makeSphere(const string& name,const double xc[3],const double r,const int nseg,const bool flip) {

    assert(zoneVec.empty());
    assert(xsp == NULL);
    assert(spost == NULL);
    assert(znost == NULL);
    assert(pbi == NULL);

    // everyone gets one zone...

    zoneVec.push_back(StZone(name));

    // shm allocation...

    const int nsp = GeodesicSphere::getSphereNodeCount(nseg);
    const int nst = GeodesicSphere::getSphereTriCount(nseg);
    init(nsp,nst);

    if (mpi_rank_shared == 0) {

      GeodesicSphere::addSphere(xsp,spost,xc,r,nseg,flip); // convention: false: fluid is inside, true: fluid is outside

      // the above doesn't add any znost...
      for (int ist = 0; ist < nst; ++ist) znost[ist] = 0;

    }

    MPI_Barrier(mpi_comm_shared);

  }

  void makeAnnularCyl(const string& name,const double x0[3],const double x1[3],const double r0,const double r1,const int nseg,const bool flip) {

    assert(zoneVec.empty());
    assert(xsp == NULL);
    assert(spost == NULL);
    assert(znost == NULL);
    assert(pbi == NULL);

    // everyone gets one zone...

    zoneVec.push_back(StZone(name));

    // shm allocation...

    const int nsp = GeomUtils::getAnnularCylNodeCount(nseg);
    const int nst = GeomUtils::getAnnularCylTriCount(nseg);
    init(nsp,nst);

    if (mpi_rank_shared == 0) {

      GeomUtils::addAnnularCyl(xsp,spost,x0,x1,r0,r1,nseg,flip); // convention: false: fluid is inside, true: fluid is outside

      // the above doesn't add any znost...
      for (int ist = 0; ist < nst; ++ist) znost[ist] = 0;

    }

    MPI_Barrier(mpi_comm_shared);

  }

  void makeTcone(const double x0[3],const double r0,const double x1[3],const double r1,const int n,const int zone_tcone) {

    // shm allocation...

    init(n*2,n*2);

    if (mpi_rank_shared == 0) {

      // n points around base...
      const double dx[3] = DIFF(x1,x0);
      const double dx_mag = MAG(dx);
      assert(dx_mag > 0.0);

      // absolute value ensures no singularity in the middle of the tcone as well as positive volume
      assert(r0 > 0.0);
      assert(r1 > 0.0);

      // first and second nodes are cap centers
      //FOR_I3 xsp[0][i] = x0[i];
      //FOR_I3 xsp[1][i] = x1[i];

      // create nodes around each cap, cap0 first [0:n-1],
      // then cap1 [n:2n-1]
      GeomUtils::createCirclePts(xsp  ,n,x0,dx,r0); // note that axis (dx) need not be normalized
      GeomUtils::createCirclePts(xsp+n,n,x1,dx,r1);

      // tcone
      GeomUtils::facetCircleToCircle(spost,znost,0,n,n,zone_tcone,true,true); // bools are: close, flip

    }

    MPI_Barrier(mpi_comm_shared);

  }

  void makeTcone(const double x0[3],const double r0,const double x1[3],const double r1,const int n,const int zone_x0,const int zone_x1,const int zone_tcone) {

    // shm allocation...

    init(2+n*2,n*4);

    if (mpi_rank_shared == 0) {

      // n points around base...
      const double dx[3] = DIFF(x1,x0);
      const double dx_mag = MAG(dx);
      assert(dx_mag > 0.0);

      // absolute value ensures no singularity in the middle of the tcone as well as positive volume
      assert(r0 > 0.0);
      assert(r1 > 0.0);

      // first and second nodes are cap centers
      FOR_I3 xsp[0][i] = x0[i];
      FOR_I3 xsp[1][i] = x1[i];

      // create nodes around each cap, cap0 first [1:n],
      // then cap1 [n+1:2n]
      GeomUtils::createCirclePts(xsp+2  ,n,x0,dx,r0); // note that axis (dx) need not be normalized
      GeomUtils::createCirclePts(xsp+2+n,n,x1,dx,r1);

      // cap0, cap1
      GeomUtils::facetCircleToPoint(spost,    znost,    0,2,  n,zone_x0,false);
      GeomUtils::facetCircleToPoint(spost+3*n,znost+3*n,1,2+n,n,zone_x1,true);

      // tcone
      GeomUtils::facetCircleToCircle(spost+n,znost+n,2,2+n,n,zone_tcone,true,true); // bools are: close, flip

    }

    MPI_Barrier(mpi_comm_shared);

  }

  void makeCylinder(const double x0[3],const double x1[3],const double r,const int n,const bool b_caps = true,const bool b_periodic = false) {

    // x0,x1 are at the centers of the left and right cylinder cap. r is the radius

    assert(zoneVec.empty());
    assert(xsp == NULL);
    assert(spost == NULL);
    assert(znost == NULL);
    assert(pbi == NULL);

    // everyone gets the zones...

    if (b_caps) {
      const int zone_x0 = 0;  zoneVec.push_back(StZone("x0"));
      const int zone_x1 = 1;  zoneVec.push_back(StZone("x1"));
      const int zone_cyl = 2; zoneVec.push_back(StZone("cyl"));
      makeTcone(x0,r,x1,r,n,zone_x0,zone_x1,zone_cyl);
      if (b_periodic) {
        const double axis[3] = DIFF(x1,x0);;
        PeriodicData::periodicTransformVec.push_back(PeriodicTransform(PERIODIC_TRANSFORM_CART,axis[0],axis[1],axis[2]));
        //pbi = new uint8[nsp];
        CTI_Mmap(pbi,nsp);
        if (mpi_rank_shared == 0) {
          pbi[0] = BitUtils::packPbiHash(0,1);
          pbi[1] = BitUtils::packPbiHash((1<<0),0);
          for(int i = 2; i < n+2; ++i) {
            pbi[i]=BitUtils::packPbiHash(0,i);
            pbi[n+i]=BitUtils::packPbiHash((1<<0),i);
          }
        }
        MPI_Barrier(mpi_comm_shared);
	zoneVec[0].setPeriodicBits(1<<1); 
	zoneVec[1].setPeriodicBits(1<<0); 
      }
    }
    else {
      const int zone_cyl = 0; zoneVec.push_back(StZone("cyl"));
      makeTcone(x0,r,x1,r,n,zone_cyl);
    }
  }

  void makeCylinder(const string& zonename,const double x0[3],const double x1[3],const double r,const int n,const bool b_caps = true) {

    // x0,x1 are at the centers of the left and right cylinder cap. r is the radius

    assert(zoneVec.empty());
    assert(xsp == NULL);
    assert(spost == NULL);
    assert(znost == NULL);
    assert(pbi == NULL);

    // everyone gets the zone...
    zoneVec.push_back(StZone(zonename));

    if (b_caps) {
      makeTcone(x0,r,x1,r,n,0,0,0);
    }
    else {
      makeTcone(x0,r,x1,r,n,0);
    }

  }

  void makeTcone(const double x0[3],const double r0,const double x1[3],const double r1,const int n,const bool b_caps = true) {

    // x0,x1 are at the centers of the left and right tcone cap. r0,r1 are the radii of the left and right cap.

    assert(zoneVec.empty());
    assert(xsp == NULL);
    assert(spost == NULL);
    assert(znost == NULL);
    assert(pbi == NULL);

    // everyone gets the zones...

    if (b_caps) {
      const int zone_x0 = 0;    zoneVec.push_back(StZone("x0"));
      const int zone_x1 = 1;    zoneVec.push_back(StZone("x1"));
      const int zone_tcone = 2; zoneVec.push_back(StZone("tcone"));
      makeTcone(x0,r0,x1,r1,n,zone_x0,zone_x1,zone_tcone);
    }
    else {
      const int zone_tcone = 0; zoneVec.push_back(StZone("tcone"));
      makeTcone(x0,r0,x1,r1,n,zone_tcone);
    }

  }

  void ensureNsp() {
    if (n_sp == NULL) {
      if (mpi_rank == 0) cout << " > allocating n_sp..." << endl;
      CTI_Mmap_rw(n_sp,nsp);
      rebuildNsp();
    }
  }

  void rebuildNsp() {
    // called after node movement, for example...
    assert(n_sp);
    if (mpi_rank == 0) cout << " > computing n_sp..." << endl;
    int nsp_split = nsp/mpi_size_shared;
    if (nsp%mpi_size_shared != 0) nsp_split += 1;
    const int isp0 = min(nsp,nsp_split*mpi_rank_shared);
    const int isp1 = min(nsp,nsp_split*(mpi_rank_shared+1));
    for (int isp = isp0; isp < isp1; ++isp) 
      FOR_I3 n_sp[isp][i] = 0.0;

    // node-based...
    ensureStosp();
    double dx[3][3];
    double dx_mag[3];
    for (int isp = isp0; isp < isp1; ++isp) {
      for (int top = stosp_i[isp]; top != stosp_i[isp+1]; ++top) {
        const int ist = stosp_v[top];
        FOR_I3 {
          FOR_J3 dx[i][j] = xsp[spost[ist][(i+1)%3]][j] - xsp[spost[ist][i]][j];
          dx_mag[i] = MAG(dx[i]);
        }
        FOR_I3 {
          if (spost[ist][i] == isp) {
            const double cp[3] = CROSS_PRODUCT(dx[(i+2)%3],dx[i]);
            const double cp_mag = MAG(cp);
            if (cp_mag > 1.0E-6*dx_mag[(i+2)%3]*dx_mag[i]) {
              const double dp = -DOT_PRODUCT(dx[(i+2)%3],dx[i]);
              const double theta = atan2(cp_mag,dp);
              assert(theta > 0.0);
              const int isp = spost[ist][i];
              FOR_J3 n_sp[isp][j] += theta*cp[j]/cp_mag;
            }
          }
        }
      }
    }

    /*
    // tri-based...
    for (int ist = 0; ist < nst; ++ist) {
    FOR_I3 {
    FOR_J3 dx[i][j] = xsp[spost[ist][(i+1)%3]][j] - xsp[spost[ist][i]][j];
    dx_mag[i] = MAG(dx[i]);
    }
    FOR_I3 {
    const double cp[3] = CROSS_PRODUCT(dx[(i+2)%3],dx[i]);
    const double cp_mag = MAG(cp);
    if (cp_mag > 1.0E-6*dx_mag[(i+2)%3]*dx_mag[i]) {
    const double dp = -DOT_PRODUCT(dx[(i+2)%3],dx[i]);
    const double theta = atan2(cp_mag,dp);
    assert(theta > 0.0);
    const int isp = spost[ist][i];
    FOR_J3 n_sp[isp][j] += theta*cp[j]/cp_mag;
    }
    }
    }
    */

    for (int isp = isp0; isp < isp1; ++isp) {
      const double mag = MAG(n_sp[isp]);
      assert(mag > 0.0);
      FOR_I3 n_sp[isp][i] /= mag; // unit normal
    }
    MPI_Barrier(mpi_comm_shared);
  }

  void ensureStost() {
    if (stost == NULL) {
      ensureStosp();
      CTI_Mmap(stost,nst);
      // XXXX TODO: must be a fast way to build this in parallel. For now, just do it
      // in serial on shared 0...
      if (mpi_rank_shared == 0) {
        for (int ist = 0; ist < nst; ++ist) {
          FOR_I3 {
            // stost takes -1 unless we find a match...
            stost[ist][i] = -1;
            const int isp0 = spost[ist][i];
            const int isp1 = spost[ist][(i+1)%3];
            for (int sos = stosp_i[isp0]; sos != stosp_i[isp0+1]; ++sos) {
              const int ist_nbr = stosp_v[sos];
              if (ist_nbr != ist) {
                if ( ((spost[ist_nbr][0] == isp1)&&(spost[ist_nbr][1] == isp0)) ||
                     ((spost[ist_nbr][1] == isp1)&&(spost[ist_nbr][2] == isp0)) ||
                     ((spost[ist_nbr][2] == isp1)&&(spost[ist_nbr][0] == isp0)) ) {
                  stost[ist][i] = ist_nbr;
                  break;
                }
              }
            }
          }
        }
      }
      MPI_Barrier(mpi_comm_shared);
    }
  }

  void ensureStosp() {
    // only take action if stosp_i/v is NULL...
    if (stosp_i == NULL) {
      assert(stosp_v == NULL);
      CTI_Mmap(stosp_i,nsp+1);
      if (mpi_rank_shared == 0) {
        for (int isp = 0; isp < nsp; ++isp)
          stosp_i[isp+1] = 0;
        for (int ist = 0; ist < nst; ++ist) {
          FOR_I3 {
            const int isp = spost[ist][i];
            ++stosp_i[isp+1]; // count
          }
        }
        // set csr...
        stosp_i[0] = 0;
        for (int isp = 0; isp < nsp; ++isp)
          stosp_i[isp+1] += stosp_i[isp];
      }
      MPI_Barrier(mpi_comm_shared);
      CTI_Mmap(stosp_v,stosp_i[nsp]);
      if (mpi_rank_shared == 0) {
        for (int ist = 0; ist < nst; ++ist) {
          FOR_I3 {
            const int isp = spost[ist][i];
            stosp_v[stosp_i[isp]++] = ist;
          }
        }
        // rewind csr...
        for (int isp = nsp-1; isp > 0; --isp)
          stosp_i[isp] = stosp_i[isp-1];
        stosp_i[0] = 0;
      }
      MPI_Barrier(mpi_comm_shared);
    }
  }

  void gatherAll(vector<int>& allVec,vector<int>& myVec) {
    
    // TODO: put this somewhere else... 
    
    int my_count = myVec.size();
    int * count = NULL;
    if (mpi_rank == 0) count = new int[mpi_size];
    MPI_Gather(&my_count,1,MPI_INT,count,1,MPI_INT,0,mpi_comm);
    
    int * disp = NULL;
    int * int_buf = NULL;
    int count_sum;
    if (mpi_rank == 0) {
      // check for overflow...
      disp = new int[mpi_size];
      disp[0] = 0;
      for (int rank = 1; rank < mpi_size; ++rank)
	disp[rank] = disp[rank-1] + count[rank-1];
      count_sum = disp[mpi_size-1] + count[mpi_size-1];
      int_buf = new int[count_sum];
    }
    int * my_ptr = (myVec.empty() ? NULL : &myVec.front());
    MPI_Gatherv(my_ptr,my_count,MPI_INT,int_buf,count,disp,MPI_INT,0,mpi_comm);
    
    assert(allVec.empty());
    if (mpi_rank == 0) {
      delete[] count;
      delete[] disp;
      set<int> intSet;
      for (int ii = 0; ii < count_sum; ++ii) {
	intSet.insert(int_buf[ii]);
      }
      delete[] int_buf;
      for (set<int>::iterator it = intSet.begin(); it != intSet.end(); ++it) {
	allVec.push_back(*it);
      }
    }
    
    int n = allVec.size();
    MPI_Bcast(&n,1,MPI_INT,0,mpi_comm);
    if (n > 0) {
      if (mpi_rank != 0) allVec.resize(n);
      MPI_Bcast(&allVec.front(),n,MPI_INT,0,mpi_comm);
    }
  
  }
  
  void gather(set<pair<int,int> >& globalPairSet,const set<pair<int,int> >& pairSet) {
    
    // TODO: put this somewhere else... 
    
    int my_count = pairSet.size()*2;
    int * count = NULL;
    if (mpi_rank == 0) count = new int[mpi_size];
    MPI_Gather(&my_count,1,MPI_INT,count,1,MPI_INT,0,mpi_comm);

    int * my_int_buf = new int[my_count];
    my_count = 0;
    for (set<pair<int,int> >::const_iterator it = pairSet.begin(); it != pairSet.end(); ++it) {
      my_int_buf[my_count] = it->first;
      my_int_buf[my_count+1] = it->second;
      my_count += 2;
    }
    
    int * disp = NULL;
    int * int_buf = NULL;
    int count_sum;
    if (mpi_rank == 0) {
      disp = new int[mpi_size];
      disp[0] = 0;
      for (int rank = 1; rank < mpi_size; ++rank)
	disp[rank] = disp[rank-1] + count[rank-1];
      count_sum = disp[mpi_size-1] + count[mpi_size-1];
      int_buf = new int[count_sum];
    }
    MPI_Gatherv(my_int_buf,my_count,MPI_INT,int_buf,count,disp,MPI_INT,0,mpi_comm);
    delete[] my_int_buf;
    
    assert(globalPairSet.empty());
    if (mpi_rank == 0) {
      delete[] count;
      delete[] disp;
      for (int ii = 0; ii < count_sum; ii += 2) {
	globalPairSet.insert(pair<int,int>(int_buf[ii],int_buf[ii+1]));
      }
      delete[] int_buf;
    }

  }
  
  void buildZnozn() {
    if (mpi_rank == 0) cout << "buildZnozn()" << endl;
    assert(znozn_i == NULL);
    assert(znozn_v == NULL);
    ensureTeost();
    // split the st's up so we all visit a few...
    int nst_split = nst/mpi_size;
    if (nst%mpi_size != 0) nst_split += 1;
    const int ist0 = min(nst,nst_split*mpi_rank);
    const int ist1 = min(nst,nst_split*(mpi_rank+1));
    set<pair<int,int> > pairSet;
    for (int ist = ist0; ist < ist1; ++ist) {
      const int izn = znost[ist];
      assert((izn >= 0)&&(izn < zoneVec.size()));
      FOR_I3 {
	int ist_nbr,i_nbr,bit_nbr;
	unpackTeost(ist_nbr,i_nbr,bit_nbr,teost[ist][i]);
	if ((ist_nbr != ist)&&(znost[ist_nbr] != izn)) {
	  pairSet.insert(pair<int,int>(izn,znost[ist_nbr]));
	}
      }
    }
    set<pair<int,int> > globalPairSet;
    gather(globalPairSet,pairSet);
    znozn_i = new int[zoneVec.size()+1];
    if (mpi_rank == 0) {
      for (int izn = 0; izn < zoneVec.size(); ++izn)
	znozn_i[izn+1] = 0;
      znozn_v = new int[globalPairSet.size()];
      int zoz = 0;
      for (set<pair<int,int> >::iterator it = globalPairSet.begin(); it != globalPairSet.end(); ++it) {
	const int izn = it->first;
	assert((izn >= 0)&&(izn < zoneVec.size()));
	++znozn_i[izn+1];
	const int izn_nbr = it->second;
	assert((izn_nbr >= 0)&&(izn_nbr < zoneVec.size()));
	znozn_v[zoz++] = izn_nbr;
      }
      assert(zoz == globalPairSet.size());
      znozn_i[0] = 0;
      for (int izn = 0; izn < zoneVec.size(); ++izn)
	znozn_i[izn+1] += znozn_i[izn];
      assert(znozn_i[zoneVec.size()] == zoz);
      for (int izn = 0; izn < zoneVec.size(); ++izn) {
	cout << " > izn: " << izn << " \"" << zoneVec[izn].getName() << "\" is connected to zones: ";
	for (int zoz = znozn_i[izn]; zoz != znozn_i[izn+1]; ++zoz) {
	  cout << znozn_v[zoz] << ",";
	}
	cout << endl;
      }
    }
    // finally Bcast the structures...
    // TODO: this could/should be shared mem?...
    MPI_Bcast(znozn_i,zoneVec.size()+1,MPI_INT,0,mpi_comm);
    if (mpi_rank != 0) {
      assert(znozn_v == NULL);
      znozn_v = new int[znozn_i[zoneVec.size()]];
    }
    MPI_Bcast(znozn_v,znozn_i[zoneVec.size()],MPI_INT,0,mpi_comm);
  }

  void ensureZnozn() {
    if (znozn_i == NULL) {
      assert(znozn_v == NULL);
      buildZnozn();
    }
  }
  
  /*
    inline void setJK(int jk[2],const double x[3]) const {
    
    // j is along surface_e1...
    jk[0] = (int)floor(surface_dx*((x[0]-surface_x0[0])*surface_e1[0] + 
    (x[1]-surface_x0[1])*surface_e1[1] + 
    (x[2]-surface_x0[2])*surface_e1[2]));
    if (jk[0]%2 == 0) jk[0] += 1; // ensure odd
    
    // and k is along surface_e2...
    jk[1] = (int)floor(surface_dx*((x[0]-surface_x0[0])*surface_e2[0] + 
    (x[1]-surface_x0[1])*surface_e2[1] + 
    (x[2]-surface_x0[2])*surface_e2[2]));
    if (jk[1]%2 == 0) jk[1] += 1; // ensure odd
    
    }

    void addXpIntersectionsFromSurface(vector<XpIntersection>& xpIntersectionVec,const int jk[2],const Adt2dShmNew * const adt2d,const int level) {
    
    vector<int> bboxVec;
    adt2d->addBboxesForPoint(bboxVec,jk);
    
    if (!bboxVec.empty()) {
      
    for (int ibb = 0; ibb < bboxVec.size(); ++ibb) {

    const int ist = bboxVec[ibb];

    // calculate the integer areas...
    int8 dij[3][2];
    FOR_I3 {
    int jk_sp[2]; setJK(jk_sp,xsp[spost[ist][i]]);
    FOR_J2 dij[i][j] = jk_sp[j] - jk[j];
    }
        
    int8 A[3];
    int8 Asum = 0;
    FOR_I3 {
    A[i] = dij[(i+1)%3][0]*dij[(i+2)%3][1] - dij[(i+1)%3][1]*dij[(i+2)%3][0];
    Asum += A[i];
    }
        
    if (Asum > 0) {
    // this is oriented towards the ray...
    // there is only an intersection if all A[i]'s are also positive (or zero)...
    if ((A[0] >= 0)&&(A[1] >= 0)&&(A[2] >= 0)) {
    //there should only be zero or one equal to zero...
    assert( ((A[0] != 0)||(A[1] != 0)||(A[2] != 0)) || ((A[0] == 0)||(A[1] != 0)||(A[2] != 0)) || ((A[0] != 0)||(A[1] == 0)||(A[2] != 0)) || ((A[0] != 0)||(A[1] != 0)||(A[2] == 0)) );
    // the double x in the e0 direction...
    double * xsp0 = xsp[spost[ist][0]];
    const double xp0 = 
    (xsp0[0]-surface_x0[0])*surface_e0[0] +  
    (xsp0[1]-surface_x0[1])*surface_e0[1] +  
    (xsp0[2]-surface_x0[2])*surface_e0[2];
    double * xsp1 = xsp[spost[ist][1]];
    const double xp1 = 
    (xsp1[0]-surface_x0[0])*surface_e0[0] +  
    (xsp1[1]-surface_x0[1])*surface_e0[1] +  
    (xsp1[2]-surface_x0[2])*surface_e0[2];
    double * xsp2 = xsp[spost[ist][2]];
    const double xp2 = 
    (xsp2[0]-surface_x0[0])*surface_e0[0] +  
    (xsp2[1]-surface_x0[1])*surface_e0[1] +  
    (xsp2[2]-surface_x0[2])*surface_e0[2];
    xpIntersectionVec.push_back(XpIntersection());
    xpIntersectionVec.back().xp = ( double(A[0])*xp0 + double(A[1])*xp1 + double(A[2])*xp2 )/double(Asum);
    xpIntersectionVec.back().level = level;
    if ((A[0] == 0)||(A[1] == 0)||(A[2] == 0)) {
    // this is an edge intersection. It gets a sign of just +1...
    xpIntersectionVec.back().sign = 1;
    }
    else {
    // all areas are positive, so sign = +2...
    xpIntersectionVec.back().sign = 2;
    }
    }
    }
    else if (Asum < 0) {
    // this is oriented AWAY FROM  the ray...
    // there is only an intersection if all A[i]'s are also negative (or zero)...
    if ((A[0] <= 0)&&(A[1] <= 0)&&(A[2] <= 0)) {
    //there should only be zero or one equal to zero...
    assert( ((A[0] != 0)||(A[1] != 0)||(A[2] != 0)) || ((A[0] == 0)||(A[1] != 0)||(A[2] != 0)) || ((A[0] != 0)||(A[1] == 0)||(A[2] != 0)) || ((A[0] != 0)||(A[1] != 0)||(A[2] == 0)) );
    // the double x in the e0 direction...
    double * xsp0 = xsp[spost[ist][0]];
    const double xp0 = 
    (xsp0[0]-surface_x0[0])*surface_e0[0] +  
    (xsp0[1]-surface_x0[1])*surface_e0[1] +  
    (xsp0[2]-surface_x0[2])*surface_e0[2];
    double * xsp1 = xsp[spost[ist][1]];
    const double xp1 = 
    (xsp1[0]-surface_x0[0])*surface_e0[0] +  
    (xsp1[1]-surface_x0[1])*surface_e0[1] +  
    (xsp1[2]-surface_x0[2])*surface_e0[2];
    double * xsp2 = xsp[spost[ist][2]];
    const double xp2 = 
    (xsp2[0]-surface_x0[0])*surface_e0[0] +  
    (xsp2[1]-surface_x0[1])*surface_e0[1] +  
    (xsp2[2]-surface_x0[2])*surface_e0[2];
    xpIntersectionVec.push_back(XpIntersection());
    xpIntersectionVec.back().xp = ( double(A[0])*xp0 + double(A[1])*xp1 + double(A[2])*xp2 )/double(Asum);
    xpIntersectionVec.back().level = level;
    if ((A[0] == 0)||(A[1] == 0)||(A[2] == 0)) {
    // this is an edge intersection. It gets a sign of just +1...
    xpIntersectionVec.back().sign = -1;
    }
    else {
    // all areas are positive, so sign = +2...
    xpIntersectionVec.back().sign = -2;
    }
    }
    }
        
    }

    }

    }

    void initPointIsInside() {

    if (adt2d)
    return;

    // set integer precision using surface bbox...

    const int BIT_MAX = 28; // do not go above 30

    if (!b_bbminmax) calcBbminmax();
    FOR_I3 surface_x0[i] = 0.5*(bbminmax[i]-bbminmax[3+i]);
    const double delta_max = max(-bbminmax[3]-bbminmax[0],max(-bbminmax[4]-bbminmax[1],-bbminmax[5]-bbminmax[2]));
    surface_dx = double(1<<BIT_MAX)/delta_max;

    // build 2d shm int adt...

    int my_nst_avg = nst/mpi_size;
    if (nst%mpi_size) ++my_nst_avg;
    const int ist0 = min(nst,mpi_rank*my_nst_avg);
    const int ist1 = min(nst,(mpi_rank+1)*my_nst_avg);
    assert(ist1-ist0 <= my_nst_avg);
    
    int my_nbb = ist1-ist0;
    vector<int> mybbVec(my_nbb*5); // 5 ints per bbox: original index (in this case ist), jmin, kmin, -jmax, -kmax
    
    for (int ibb = 0; ibb < my_nbb; ++ibb) {
    const int ist = ist0 + ibb;
    mybbVec[5*ibb] = ist;
    int jk0[2]; setJK(jk0,xsp[spost[ist][0]]);
    int jk1[2]; setJK(jk1,xsp[spost[ist][1]]);
    int jk2[2]; setJK(jk2,xsp[spost[ist][2]]);
    FOR_I2 {
    mybbVec[5*ibb+1+i]   = min( jk0[i],min( jk1[i], jk2[i]));
    mybbVec[5*ibb+1+i+2] = min(-jk0[i],min(-jk1[i],-jk2[i]));
    }
    }
    
    adt2d = new Adt2dShmNew(mybbVec,my_nbb);

    }

    bool pointIsInside(const double xp[3]) {

    if (adt2d == NULL) {
    CERR("call initIsInisde before using pointIsInside()");
    }

    // convert the passed point to j,k coords...

    int jk[2]; setJK(jk,xp);
    vector<XpIntersection> xpIntersectionVec;
    addXpIntersectionsFromSurface(xpIntersectionVec,jk,adt2d,0);

    // sort along xp...
    sort(xpIntersectionVec.begin(),xpIntersectionVec.end());

    int count = 0;
    for (int ii = 0, lim = xpIntersectionVec.size(); ii < lim; ++ii) {
    if (xpIntersectionVec[ii].xp > xp[0])
    break;
    assert(xpIntersectionVec[ii].level == 0);
    count += xpIntersectionVec[ii].sign;
    }

    if (count == 0) {
    return false;
    }
    else {
    return true;
    }
    }
  */
  
  // =================================================================
  // new point-is-inside logic for SurfaceShm based on sign of nearest 
  // tri in the x-direction...
  // =================================================================
  
  void initPointIsInside(const double bbmin3d[3],const double bbmax3d[3]) {

    // here we project everything in x by default, so the important 
    // dimensions of the bbox are [1] and [2]...
    
    // if we already have an adt2d, then check its stored bbox range. If sufficient
    // then just return...
    
    if (adt2d) {
      if ((bbmin3d[1] >= adt2d_bbmin[0])&&(bbmin3d[2] >= adt2d_bbmin[1])&&
          (bbmax3d[1] <= adt2d_bbmax[0])&&(bbmax3d[2] <= adt2d_bbmax[1])) {
        return;
      }
      else {
        delete adt2d;
        adt2d_st_vec.clear();
      }
    }
    
    // take the new bbox as the y and z components of the passed bbox3d...  
    
    FOR_I2 {
      adt2d_bbmin[i] = bbmin3d[i+1];
      adt2d_bbmax[i] = bbmax3d[i+1];
    }
    
    // adt2d_st_vec stores the lookup back to the surface tri...
    assert(adt2d_st_vec.empty());
    for (int ist = 0; ist < nst; ++ist) {
      if ((max(xsp[spost[ist][0]][1],max(xsp[spost[ist][1]][1],xsp[spost[ist][2]][1])) < adt2d_bbmin[0])||
          (min(xsp[spost[ist][0]][1],min(xsp[spost[ist][1]][1],xsp[spost[ist][2]][1])) > adt2d_bbmax[0])) 
        continue;
      if ((max(xsp[spost[ist][0]][2],max(xsp[spost[ist][1]][2],xsp[spost[ist][2]][2])) < adt2d_bbmin[1])||
          (min(xsp[spost[ist][0]][2],min(xsp[spost[ist][1]][2],xsp[spost[ist][2]][2])) > adt2d_bbmax[1])) 
        continue;
      adt2d_st_vec.push_back(ist);
    }
    
    int nbb = adt2d_st_vec.size();
    double (*bbmin)[2] = new double[nbb][2];
    double (*bbmax)[2] = new double[nbb][2];
    for (int ii = 0; ii < nbb; ++ii) {
      const int ist = adt2d_st_vec[ii];
      FOR_I2 bbmin[ii][i] = min(xsp[spost[ist][0]][i+1],min(xsp[spost[ist][1]][i+1],xsp[spost[ist][2]][i+1]));
      FOR_I2 bbmax[ii][i] = max(xsp[spost[ist][0]][i+1],max(xsp[spost[ist][1]][i+1],xsp[spost[ist][2]][i+1]));
    }
    
    adt2d = new Adt2d<double>(nbb,bbmin,bbmax);
    
    delete[] bbmin;
    delete[] bbmax;
    
  }

  void initPointIsInside(const double (* const xp)[3],const int np) {
  
    // another way to initialize the point-is-inside function is to 
    // pass ALL the points to be queried. This simply builds the 
    // bbox and then calls the regular bbox-based initialization...  
  
    double bbmin3d[3] = { HUGE_VAL, HUGE_VAL, HUGE_VAL };
    double bbmax3d[3] = { -HUGE_VAL, -HUGE_VAL, -HUGE_VAL };
    
    for (int ip = 0; ip < np; ++ip) {
      FOR_I3 {
        bbmin3d[i] = min(bbmin3d[i],xp[ip][i]);
        bbmax3d[i] = max(bbmax3d[i],xp[ip][i]);
      }
    }
    
    initPointIsInside(bbmin3d,bbmax3d);
    
  }
  
  bool pointIsInside(const double xp[3]) const {

    if (!adt2d) {
      CERR("SurfaceShm::pointIsInside: must call initPointIsInside with the full range of points to be queried");
    }
    
    // check bbox -- recall it is a 2d box in y-z:
    if ((xp[1] < adt2d_bbmin[0])||(xp[1] > adt2d_bbmax[0])||(xp[2] < adt2d_bbmin[1])||(xp[2] > adt2d_bbmax[1])) {
      CERR("Surface::pointIsInside: point is outside of adt2d. check initPointIsInside");
    }

    // note: to handle robustly the problem of a concave cusp, we insist that 
    // there is finite positive wgt associated with every sp of st being considered.
    // This should break ties associated with the cusp edge intersection...
    
    const double eps = 1.0E-10; // minimum normalized weight of a point on a tri  
    
    // grab candidates...
    double dx_min = HUGE_VAL;
    bool b_in = false; // default for any point is false
    vector<int> intVec;
    adt2d->buildListForPoint(intVec,xp+1);
    for (int ii = 0; ii < intVec.size(); ++ii) {
      const int ist = adt2d_st_vec[intVec[ii]]; assert((ist >= 0)&&(ist < nst));
      double A[3];
      double Asum = 0.0;
      FOR_I3 {
        // consider the normal associated with the point and each edge opposite the node...
        // here we compute 3 different areas to handle double precision round-off robustly
        const double A0 = (xsp[spost[ist][(i+2)%3]][1]-xp[1])*(xsp[spost[ist][(i+1)%3]][2]-xp[2]) - 
          (xsp[spost[ist][(i+2)%3]][2]-xp[2])*(xsp[spost[ist][(i+1)%3]][1]-xp[1]);
        const double A1 = (xp[1]-xsp[spost[ist][(i+1)%3]][1])*(xsp[spost[ist][(i+2)%3]][2]-xsp[spost[ist][(i+1)%3]][2]) - 
          (xp[2]-xsp[spost[ist][(i+1)%3]][2])*(xsp[spost[ist][(i+2)%3]][1]-xsp[spost[ist][(i+1)%3]][1]);
        const double A2 = (xsp[spost[ist][(i+1)%3]][1]-xsp[spost[ist][(i+2)%3]][1])*(xp[2]-xsp[spost[ist][(i+2)%3]][2]) - 
          (xsp[spost[ist][(i+1)%3]][2]-xsp[spost[ist][(i+2)%3]][2])*(xp[1]-xsp[spost[ist][(i+2)%3]][1]);
        if ( ((A0 > 0.0)&&(A1 > 0.0)&&(A2 > 0.0)) || ((A0 < 0.0)&&(A1 < 0.0)&&(A2 < 0.0)) ) {
          A[i] = A0 + A1 + A2;
          Asum += A[i];
        }
        else {
          A[i] = 0.0;
        }
      }
      // only consider tris that project in x...
      if (fabs(Asum) > 0.0) {
        // now find the closest tri that has an intersection...
        if ((A[0] >= 0.0)&&(A[1] >= 0.0)&&(A[2] >= 0.0)) {
          // all areas positive or zero...
          double Asum2 = 0.0;
          double dx = 0.0;
          double yz_check[2] = { 0.0, 0.0 }; // TODO: get rid of this
          FOR_I3 {
            const double wgt = max(A[i],eps*Asum);
            dx += wgt*(xsp[spost[ist][i]][0]-xp[0]);
            yz_check[0] += wgt*(xsp[spost[ist][i]][1]-xp[1]);
            yz_check[1] += wgt*(xsp[spost[ist][i]][2]-xp[2]);
            Asum2 += wgt;
          }
          dx /= Asum2;
          yz_check[0] /= Asum2;
          yz_check[1] /= Asum2;
          //assert(fabs(yz_check[0]) < 1.0E-10);
          //assert(fabs(yz_check[1]) < 1.0E-10);
          if (fabs(dx) < dx_min) {
            dx_min = fabs(dx);
            b_in = (dx <= 0.0);
          }
        }
        else if ((A[0] <= 0.0)&&(A[1] <= 0.0)&&(A[2] <= 0.0)) {
          // all areas negative or zero...
          double Asum2 = 0.0;
          double dx = 0.0;
          double yz_check[2] = { 0.0, 0.0 }; // get rid of this
          FOR_I3 {
            const double wgt = min(A[i],eps*Asum);
            dx += wgt*(xsp[spost[ist][i]][0]-xp[0]);
            yz_check[0] += wgt*(xsp[spost[ist][i]][1]-xp[1]);
            yz_check[1] += wgt*(xsp[spost[ist][i]][2]-xp[2]);
            Asum2 += wgt;
          }
          dx /= Asum2;
          yz_check[0] /= Asum2;
          yz_check[1] /= Asum2;
          //assert(fabs(yz_check[0]) < 1.0E-10);
          //assert(fabs(yz_check[1]) < 1.0E-10);
          if (fabs(dx) < dx_min) {
            dx_min = fabs(dx);
            b_in = (dx >= 0.0);
          }
        }
      }
    }
    
    return b_in;

  }
  

  void intersectWithTcone(const double axis[3], const double ra, const double rb, const double xa[3], const double xb[3], vector<vector<double> >& xp_r, int * st_flag) {
    //not quite general now, assume the axis is X or Z

    ensureStost();
    
    int ist_start;
    int isp0_start, isp1_start;
    int count = 0;
    double xp[3];

    FOR_IST {
      if (st_flag[ist] == 1) {
	FOR_I3 {
          const int isp0 = spost[ist][i];
          const int isp1 = spost[ist][(i+1)%3];

          /*
	    const double rsp0 = sqrt(xsp[isp0][1]*xsp[isp0][1] + xsp[isp0][2]*xsp[isp0][2]);
	    const double r0 = ra + (xsp[isp0][0] - xa[0])/(xb[0] - xa[0])*(rb - ra);
	    const double rsp1 = sqrt(xsp[isp1][1]*xsp[isp1][1] + xsp[isp1][2]*xsp[isp1][2]);
	    const double r1 = ra + (xsp[isp1][0] - xa[0])/(xb[0] - xa[0])*(rb - ra);
          */

	  if ( isAnIntersection(axis, ra, rb, xa, xb, isp0, isp1) ) {
	    //cout << "found an intersecting tri: " << endl;
	    //cout << "xsp[isp0] = " << COUT_VEC(xsp[isp0]) << ", xsp[isp1] = " << COUT_VEC(xsp[isp1]) << endl;
	    ++count;

            /*
	      const double c0 = (rb-ra)/(xb[0]-xa[0]);
	      const double x0 = xsp[isp0][0];
	      const double y0 = xsp[isp0][1];
	      const double z0 = xsp[isp0][2];
	      const double x1 = xsp[isp1][0];
	      const double y1 = xsp[isp1][1];
	      const double z1 = xsp[isp1][2];
	      const double c1 = (y1 - y0)/(x1 - x0);
	      const double c2 = (z1 - z0)/(x1 - x0);

	      const double a = c1*c1 + c2*c2 - c0*c0;
	      const double b = 2.0*(y0*c1 - x0*c1*c1 + z0*c2 - x0*c2*c2 - ra*c0 + xa[0]*c0*c0);
	      const double c = (c1*x0-y0)*(c1*x0-y0) + (c2*x0-z0)*(c2*x0-z0) - (c0*xa[0]-ra)*(c0*xa[0]-ra);

	      const double root = b*b - 4.0*a*c;
	      assert(root>0);

	      const double root0 = (-b + sqrt(root)) /2.0 /a;
	      const double root1 = (-b - sqrt(root)) /2.0 /a;

	      cout << "roots = " << root0 << " " << root1 << endl;

	      double xp[3];

	      if ( (root0>=x0 && root0<=x1) || (root0>=x1 && root0<=x0) ) {
	      xp[0] = root0;
	      }
	      else {
	      xp[0] = root1;
	      assert((root1>=x0 && root1<=x1) || (root1>=x1 && root1<=x0));
	      }

	      xp[1] = y0 + (xp[0]-x0)/(x1-x0)*(y1-y0);
	      xp[2] = z0 + (xp[0]-x0)/(x1-x0)*(z1-z0);
            */
            getIntersectionPt(axis, ra, rb, xa, xb, isp0, isp1, xp);
            //cout << "xp = " << COUT_VEC(xp) << endl;

            ist_start = ist;
            isp0_start = isp0;
            isp1_start = isp1;
	    vector<double> x_r;
	    FOR_I3 x_r.push_back(xp[i]);
            xp_r.push_back(x_r);
	  }
	}
	if (count == 2) break;
      }
    }
    if (count == 0) CERR("Cannot find any intersecting tri!");

    int ist_current = ist_start;
    int isp0_current = isp0_start;
    int isp1_current = isp1_start;
    int ist_next = -1;

    while ( ist_current != ist_start || count <= 2 ) {
      if (ist_next != ist_current) {
        FOR_I3 {
          ist_next = stost[ist_current][i];
	  bool found = false;
	  FOR_J3 {
	    const int isp0 = spost[ist_next][j];
	    const int isp1 = spost[ist_next][(j+1)%3];
	    if (isp0 == isp1_current && isp1 == isp0_current) {
	      //cout << "ist_next = " << ist_next <<endl;
	      found = true;
	      break;
	    }
	  }
	  if (found) {
	    //cout << "ist_currnet, spost, isp currents = " << ist_current << " " << COUT_VEC(spost[ist_current]) << " " << isp0_current << " " << isp1_current << endl;
	    ist_current = ist_next;
	    //cout << "now ist_current, spost = " << ist_current << " " << COUT_VEC(spost[ist_current]) << endl;
	    break;
	  }
          //}
        }
        if (ist_next != ist_current) CERR("Cannot find next tri. Searching failed!");
      }
      else {
        ist_next = -1;
        FOR_I3 {
          const int isp0 = spost[ist_current][i];
          const int isp1 = spost[ist_current][(i+1)%3];
          
          if ( isp0 != isp1_current || isp1 !=isp0_current ) {
            if ( isAnIntersection(axis, ra, rb, xa, xb, isp0, isp1) ) {
              ++count;
              getIntersectionPt(axis, ra, rb, xa, xb, isp0, isp1, xp);
              //cout << "xp = " << COUT_VEC(xp) << endl;
	      vector<double> x_r;
	      FOR_I3 x_r.push_back(xp[i]);
              xp_r.push_back(x_r);
              //cout << "old isp0_current, isp1_current = " << isp0_current << " " << isp1_current << endl;
              isp0_current = isp0;
              isp1_current = isp1;
              //cout << "now ist_current, isp0_current, isp1_current = " << ist_current << " " << isp0_current << " " << isp1_current << endl;
              break;
            }
          }
        }
      }
    }
 
    if (mpi_rank == 0) cout << "count = " << count << endl;
  }

  bool isAnIntersection(const double axis[3], const double ra, const double rb, const double xa[3], const double xb[3], const int isp0, const int isp1) {

    const double x_axis[3] = {1.0, 0.0, 0.0};
    const double z_axis[3] = {0.0, 0.0, 1.0};
    if (DIST(axis,x_axis) < 1.0e-14) {
      // axis is in x direction
      const double rsp0 = sqrt(xsp[isp0][1]*xsp[isp0][1] + xsp[isp0][2]*xsp[isp0][2]);
      const double r0 = ra + (xsp[isp0][0] - xa[0])/(xb[0] - xa[0])*(rb - ra);
      const double rsp1 = sqrt(xsp[isp1][1]*xsp[isp1][1] + xsp[isp1][2]*xsp[isp1][2]);
      const double r1 = ra + (xsp[isp1][0] - xa[0])/(xb[0] - xa[0])*(rb - ra);
      return ((rsp0 >= r0) && (rsp1 < r1)) || ((rsp0 < r0) && (rsp1 >= r1));
    }
    else if (DIST(axis,z_axis) < 1.0e-14) {
      // axis is in z direction
      const double rsp0 = sqrt(xsp[isp0][0]*xsp[isp0][0] + xsp[isp0][1]*xsp[isp0][1]);
      const double r0 = ra + (xsp[isp0][2] - xa[2])/(xb[2] - xa[2])*(rb - ra);
      const double rsp1 = sqrt(xsp[isp1][0]*xsp[isp1][0] + xsp[isp1][1]*xsp[isp1][1]);
      const double r1 = ra + (xsp[isp1][2] - xa[2])/(xb[2] - xa[2])*(rb - ra);
      return ((rsp0 >= r0) && (rsp1 < r1)) || ((rsp0 < r0) && (rsp1 >= r1));
    }
    else {
      CERR("isAnIntersection: unsupported axis!")
	return false;
    }
  }

  void getIntersectionPt(const double axis[3], const double ra, const double rb, const double xa[3], const double xb[3], const int isp0, const int isp1, double xp[3]) {

    const double x0 = xsp[isp0][0];
    const double y0 = xsp[isp0][1];
    const double z0 = xsp[isp0][2];
    const double x1 = xsp[isp1][0];
    const double y1 = xsp[isp1][1];
    const double z1 = xsp[isp1][2];

    const double x_axis[3] = {1.0, 0.0, 0.0};
    const double z_axis[3] = {0.0, 0.0, 1.0};

    if (DIST(axis,x_axis) < 1.0e-14) {
      // axis is in x direction
      if (fabs(x1-x0)<1.0e-13) {
        const double r = ra + (x0 - xa[0])/(xb[0] - xa[0])*(rb - ra);
        const double r0 = sqrt(y0*y0 + z0*z0);
        const double r1 = sqrt(y1*y1 + z1*z1);
        const double factor = (r-r0)/(r1-r0);
        xp[0] = x0;
        xp[1] = y0 + factor*(y1-y0);
        xp[2] = z0 + factor*(z1-z0);
        // check
        //cout << "edge case, xsp[isp0], xsp[isp1] = " << COUT_VEC(xsp[isp0]) << " " << COUT_VEC(xsp[isp1]) << endl;
        //cout << "xp = " << COUT_VEC(xp) << endl;
      }
      else {
        const double c0 = (rb-ra)/(xb[0]-xa[0]);
        const double c1 = (y1 - y0)/(x1 - x0);
        const double c2 = (z1 - z0)/(x1 - x0);

        const double a = c1*c1 + c2*c2 - c0*c0;
        const double b = 2.0*(y0*c1 - x0*c1*c1 + z0*c2 - x0*c2*c2 - ra*c0 + xa[0]*c0*c0);
        const double c = (c1*x0-y0)*(c1*x0-y0) + (c2*x0-z0)*(c2*x0-z0) - (c0*xa[0]-ra)*(c0*xa[0]-ra);

        const double root = b*b - 4.0*a*c;
        assert(root>=0.0);

        const double root0 = (-b + sqrt(root)) /2.0 /a;
        const double root1 = (-b - sqrt(root)) /2.0 /a;

        //cout << "roots = " << root0 << " " << root1 << endl;

        if ( (root0>=x0 && root0<=x1) || (root0>=x1 && root0<=x0) ) {
          xp[0] = root0;
        }
        else {
          xp[0] = root1;
          assert((root1>=x0 && root1<=x1) || (root1>=x1 && root1<=x0));
        }

        if (fabs(xp[0]-x0) < 1.0e-13 || fabs(xp[0]-x1) < 1.0e-13 ) {
          if (mpi_rank == 0) cout << "getIntersectionPt: algorithm may fail!" << endl;
        }
        xp[1] = y0 + (xp[0]-x0)/(x1-x0)*(y1-y0);
        xp[2] = z0 + (xp[0]-x0)/(x1-x0)*(z1-z0);
      }
    }
    else if (DIST(axis,z_axis) < 1.0e-14) {
      // axis is in z direction
      if (fabs(z1-z0)<1.0e-13) {
        const double r = ra + (z0 - xa[2])/(xb[2] - xa[2])*(rb - ra);
        const double r0 = sqrt(x0*x0 + y0*y0);
        const double r1 = sqrt(x1*x1 + y1*y1);
        const double factor = (r-r0)/(r1-r0);
        xp[0] = x0 + factor*(x1-x0);
        xp[1] = y0 + factor*(y1-y0);
        xp[2] = z0;
        // check
        //cout << "edge case, xsp[isp0], xsp[isp1] = " << COUT_VEC(xsp[isp0]) << " " << COUT_VEC(xsp[isp1]) << endl;
        //cout << "xp = " << COUT_VEC(xp) << endl;
      }
      else {
        const double c0 = (rb-ra)/(xb[2]-xa[2]);
        const double c1 = (x1 - x0)/(z1 - z0);
        const double c2 = (y1 - y0)/(z1 - z0);

        const double a = c1*c1 + c2*c2 - c0*c0;
        const double b = 2.0*(x0*c1 - z0*c1*c1 + y0*c2 - z0*c2*c2 - ra*c0 + xa[2]*c0*c0);
        const double c = (c1*z0-x0)*(c1*z0-x0) + (c2*z0-y0)*(c2*z0-y0) - (c0*xa[2]-ra)*(c0*xa[2]-ra);

        const double root = b*b - 4.0*a*c;
        assert(root>=0.0);

        const double root0 = (-b + sqrt(root)) /2.0 /a;
        const double root1 = (-b - sqrt(root)) /2.0 /a;

        //cout << "roots = " << root0 << " " << root1 << endl;

        if ( (root0>=z0 && root0<=z1) || (root0>=z1 && root0<=z0) ) {
          xp[2] = root0;
        }
        else {
          xp[2] = root1;
          assert((root1>=z0 && root1<=z1) || (root1>=z1 && root1<=z0));
        }

        if (fabs(xp[2]-z0) < 1.0e-13 || fabs(xp[2]-z1) < 1.0e-13 ) {
          if (mpi_rank == 0) cout << "getIntersectionPt: algorithm may fail!" << endl;
        }
        xp[0] = x0 + (xp[2]-z0)/(z1-z0)*(x1-x0);
        xp[1] = y0 + (xp[2]-z0)/(z1-z0)*(y1-y0);
      }
    }
    else {
      CERR("getIntersectionPt: unsupported axis!")
	}
  }
  
  void ensureTeost() {
    
    if (teost == NULL) buildTeost();

  }

  void clearTeost() {
    
    if (teost) {
      CTI_Munmap(teost,nst);
      teost = NULL;
    }

  }
  
  void buildTeost() {

    // START = new delete buildTeost XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

    if (mpi_rank == 0) cout << "buildTeost(): nst: " << nst << " nsp: " << nsp << endl;

    assert(teost == NULL);
    
    // partition nodes...
    
    const int nsp_local_target = (nsp + mpi_size - 1)/mpi_size;
    const int isp_begin = min(nsp,mpi_rank*nsp_local_target);
    const int isp_end = min(nsp,(mpi_rank+1)*nsp_local_target);
    int nsp_local = isp_end-isp_begin;
    
    int * stosp_local_i = new int[nsp_local+1];
    for (int isp_local = 0; isp_local < nsp_local; ++isp_local)
      stosp_local_i[isp_local+1] = 0;

    // this is painful, but do a FULL ist loop here, storing
    // the participating ist's in a vec to avoid doing this twice...

    int nst_boundary = 0;
    vector<int> st_vec;
    for (int ist = 0; ist < nst; ++ist) {
      // conside boundary faces only: i.e. skip periodic...
      const int isz = znost[ist];
      assert((isz >= 0)&&(isz < zoneVec.size()));
      if (zoneVec[isz].isBoundary()) {
        ++nst_boundary;
        bool b_in = false;
        FOR_I3 {
          const int isp = spost[ist][i];
          if ((isp >= isp_begin)&&(isp < isp_end)) {
            const int isp_local = isp-isp_begin;
            ++stosp_local_i[isp_local+1];
            b_in = true;
          }
        }
        if (b_in) st_vec.push_back(ist);
      }
    }
    
    stosp_local_i[0] = 0;
    for (int isp_local = 0; isp_local < nsp_local; ++isp_local)
      stosp_local_i[isp_local+1] += stosp_local_i[isp_local];
    const int stosp_local_s = stosp_local_i[nsp_local];
    int * stosp_local_v = new int[stosp_local_s];
    
    for (int ii = 0; ii < st_vec.size(); ++ii) {
      const int ist = st_vec[ii];
      FOR_I3 {
        const int isp = spost[ist][i];
        if ((isp >= isp_begin)&&(isp < isp_end)) {
          const int isp_local = isp-isp_begin;
          stosp_local_v[stosp_local_i[isp_local]++] = ist;
        }
      }
    }
    
    for (int isp_local = nsp_local-1; isp_local > 0; --isp_local)
      stosp_local_i[isp_local] = stosp_local_i[isp_local-1];
    stosp_local_i[0] = 0;
    
    int * send_count = new int[mpi_size];
    FOR_RANK send_count[rank] = 0;
    vector<int> intVec;
    vector<int> periodicIntVec;
    
    // now figure out which tris are paired with which...
    
    int my_count[4] = { 0, 0, 0, 0 };
    for (int isp0_local = 0; isp0_local < nsp_local; ++isp0_local) {
      const int isp0 = isp0_local + isp_begin;
      for (int sos = stosp_local_i[isp0_local]; sos != stosp_local_i[isp0_local+1]; ++sos) {
        const int ist = stosp_local_v[sos];
        // find the connecting node AND the next node...
        int i;
        for (i = 0; i < 3; ++i)
          if (spost[ist][i] == isp0)
            break;
        assert(i < 3); // must be found
        const int isp1 = spost[ist][(i+1)%3];
        assert(isp0 != isp1);
        // now search for a tri that shares these 2 nodes...
        int sos2,ist2,j;
        for (sos2 = stosp_local_i[isp0_local]; sos2 != stosp_local_i[isp0_local+1]; ++sos2) {
          ist2 = stosp_local_v[sos2];
          if (ist2 != ist) {
            for (j = 0; j < 3; ++j) {
              if ((spost[ist2][j] == isp1)&&(spost[ist2][(j+1)%3] == isp0))
                break;
            }
            if (j < 3) {
              // this is a match...
              break;
            }
          }
        }
        if (sos2 != stosp_local_i[isp0_local+1]) {
          // this is a match: 
          // ist,edge i is paired with ist2,edge j
          if (ist < ist2) {
            ++my_count[0];
            intVec.push_back(ist);
            intVec.push_back(i);
            intVec.push_back(ist2);
            intVec.push_back(j);
          }
        }
        else {
          // there was no match. This could be a tri along a periodic boundary.
          // For half these tris, this is indicated by pbi...
          assert(pbi);
          int bits0,isp0_parent;
          BitUtils::unpackPbiHash(bits0,isp0_parent,pbi[isp0]);
          int bits1,isp1_parent;
          BitUtils::unpackPbiHash(bits1,isp1_parent,pbi[isp1]);
          if ((bits0 != 0)&&(bits1 != 0)) {
            assert(bits0 == 1);
            assert(bits1 == 1);
            ++my_count[1];
            periodicIntVec.push_back(ist);
            periodicIntVec.push_back(i);
            // we are going to send this to the isp0_parent and look for isp1_parent->isp0_parent
            const int rank = isp0_parent/nsp_local_target;
            assert((rank >= 0)&&(rank < mpi_size));
            send_count[rank] += 5; // message to rank: ist,i,bit,isp1_parent,isp0_parent
          }
          else {
            // this is (probably) a periodic face, but on the boundary
            // with no pbi defined...
            assert(bits0 == 0);
            assert(bits1 == 0);
            ++my_count[2];
          }
        }
      }
    }
    
    // for periodic surfaces, process the periodicIntVec...

    int * send_disp = new int[mpi_size];
    int * send_buf_int = NULL;
    int * recv_count = new int[mpi_size];
    int * recv_disp = new int[mpi_size];
    int * recv_buf_int = NULL;

    if (pbi) {
      
      // for parallel, use alltoall communication to build additional 
      // edge links on the parent ranks...

      send_disp[0] = 0;
      for (int rank = 1; rank < mpi_size; ++rank)
        send_disp[rank] = send_count[rank-1] + send_disp[rank-1];
      const int send_count_sum = send_disp[mpi_size-1] + send_count[mpi_size-1];
      
      send_buf_int = new int[send_count_sum];
      
      for (int ii = 0; ii < periodicIntVec.size(); ii += 2) {
        const int ist = periodicIntVec[ii]; assert((ist >= 0)&&(ist < nst));
        const int i = periodicIntVec[ii+1]; assert((i >= 0)&&(i < 3));
        const int isp0 = spost[ist][i];
        const int isp1 = spost[ist][(i+1)%3];
        int bits0,isp0_parent;
        BitUtils::unpackPbiHash(bits0,isp0_parent,pbi[isp0]);
        int bits1,isp1_parent;
        BitUtils::unpackPbiHash(bits1,isp1_parent,pbi[isp1]);
        assert((bits0 != 0)&&(bits1 != 0));
        assert(bits0 == 1);
        assert(bits1 == 1);
        const int bit = 1; 
        // NOTE: for multiple periodicity (i.e. double or triple) you will need
        // to determine the bit that is common between bits0 and bits1. 
        const int rank = isp0_parent/nsp_local_target;
        send_buf_int[send_disp[rank]  ] = ist;
        send_buf_int[send_disp[rank]+1] = i;
        send_buf_int[send_disp[rank]+2] = bit;
        send_buf_int[send_disp[rank]+3] = isp0_parent;
        send_buf_int[send_disp[rank]+4] = isp1_parent;
        send_disp[rank] += 5;
      }
      periodicIntVec.clear();
      
      // rewind...
      
      send_disp[0] = 0;
      for (int rank = 1; rank < mpi_size; ++rank)
        send_disp[rank] = send_count[rank-1] + send_disp[rank-1];
      
      // set up recv-side stuff...
      
      MPI_Alltoall(send_count,1,MPI_INT,recv_count,1,MPI_INT,mpi_comm);

      recv_disp[0] = 0;
      for (int rank = 1; rank < mpi_size; ++rank)
        recv_disp[rank] = recv_count[rank-1] + recv_disp[rank-1];
      const int recv_count_sum = recv_disp[mpi_size-1] + recv_count[mpi_size-1];
      
      recv_buf_int = new int[recv_count_sum];
      MPI_Alltoallv(send_buf_int,send_count,send_disp,MPI_INT,
                    recv_buf_int,recv_count,recv_disp,MPI_INT,mpi_comm);
      delete[] send_buf_int; send_buf_int = NULL;
      
      for (int irecv = 0; irecv < recv_count_sum; irecv += 5) {
        const int ist = recv_buf_int[irecv];
        const int i   = recv_buf_int[irecv+1];
        const int bit = recv_buf_int[irecv+2];
        const int isp0 = recv_buf_int[irecv+3];
        const int isp1 = recv_buf_int[irecv+4];
        // isp1 should be local...
        assert((isp0 >= isp_begin)&&(isp0 < isp_end));
        const int isp0_local = isp0-isp_begin;
        // loop through the tris that have isp1_local and find the one with isp0...
        int sos2,ist2,j;
        for (sos2 = stosp_local_i[isp0_local]; sos2 != stosp_local_i[isp0_local+1]; ++sos2) {
          ist2 = stosp_local_v[sos2];
          if (ist2 != ist) {
            for (j = 0; j < 3; ++j) {
              if ((spost[ist2][j] == isp1)&&(spost[ist2][(j+1)%3] == isp0))
                break;
            }
            if (j < 3) {
              // this is a match...
              break;
            }
          }
        }
        if (sos2 != stosp_local_i[isp0_local+1]) {
          // we matched!
          ++my_count[3];
          intVec.push_back(ist);
          intVec.push_back(i|(bit<<2));
          intVec.push_back(ist2);
          intVec.push_back(j);
        }
        else {
          cout << "expected a match!" << endl;
          assert(0);
        }
      }
      
      // cleanup...
      
      delete[] recv_buf_int; recv_buf_int = NULL;
      
    }
    else {
      assert(periodicIntVec.empty());
    }

    delete[] stosp_local_i; stosp_local_i = NULL;
    delete[] stosp_local_v; stosp_local_v = NULL;
    
    // report and check...
    
    int count[4];
    MPI_Reduce(my_count,count,4,MPI_INT,MPI_SUM,0,mpi_comm);
    if (mpi_rank == 0) {
      
      cout << " > got internal edge links: " << count[0] << " and periodic edge links: " << count[1] << " " << count[2] << " " << count[3] << endl;
      cout << " > total nst: " << nst << " boundary (i.e. non-periodic) nst: " << nst_boundary << endl;
      
      // "Euler" check...
      assert(nst_boundary*3 == 2*(count[0]+count[1]));

      cout << " > passed Euler check. Surface is properly walkable" << endl;
      
    }

    // partition the tris...
    
    const int nst_local_target = (nst + mpi_size - 1)/mpi_size;
    const int ist_begin = min(nst,mpi_rank*nst_local_target);
    const int ist_end = min(nst,(mpi_rank+1)*nst_local_target);
    int nst_local = ist_end-ist_begin;

    // and send 
    
    assert(send_count);
    FOR_RANK send_count[rank] = 0;
    assert(intVec.size()%4 == 0);
    for (int ii = 0; ii < intVec.size(); ii += 4) {
      const int ist =       intVec[ii]; assert((ist >= 0)&&(ist < nst));
      const int i_and_bit = intVec[ii+1];
      const int i = (i_and_bit&3); assert((i >= 0)&&(i < 3));
      const int bit = (i_and_bit>>2); assert((bit >= 0)&&(bit <= 1)); // for now 0 or 1
      const int ist2 =      intVec[ii+2]; assert((ist2 >= 0)&&(ist2 < nst));
      const int j =         intVec[ii+3]; assert((j >= 0)&&(j < 3));
      const int rank = ist/nst_local_target;
      send_count[rank] += 4;
      const int rank2 = ist2/nst_local_target;
      if (rank2 != rank)
        send_count[rank2] += 4;
    }
    
    assert(send_disp);
    send_disp[0] = 0;
    for (int rank = 1; rank < mpi_size; ++rank)
      send_disp[rank] = send_count[rank-1] + send_disp[rank-1];
    const int send_count_sum = send_disp[mpi_size-1] + send_count[mpi_size-1];

    assert(send_buf_int == NULL); send_buf_int = new int[send_count_sum];
    
    // pack...
    
    for (int ii = 0; ii < intVec.size(); ii += 4) {
      const int ist =       intVec[ii]; assert((ist >= 0)&&(ist < nst));
      const int i_and_bit = intVec[ii+1];
      const int i = (i_and_bit&3); assert((i >= 0)&&(i < 3));
      const int bit = (i_and_bit>>2); assert((bit >= 0)&&(bit <= 1)); // for now 0 or 1
      const int ist2 =      intVec[ii+2]; assert((ist2 >= 0)&&(ist2 < nst));
      const int j =         intVec[ii+3]; assert((j >= 0)&&(j < 3));
      const int rank = ist/nst_local_target;
      send_buf_int[send_disp[rank]  ] = ist;
      send_buf_int[send_disp[rank]+1] = i_and_bit;
      send_buf_int[send_disp[rank]+2] = ist2;
      send_buf_int[send_disp[rank]+3] = j;
      send_disp[rank] += 4;
      const int rank2 = ist2/nst_local_target;
      if (rank2 != rank) {
        send_buf_int[send_disp[rank2]  ] = ist2;
        send_buf_int[send_disp[rank2]+1] = (j|(bit<<3)); // NOTE: shift one more for this periodicity
        send_buf_int[send_disp[rank2]+2] = ist;
        send_buf_int[send_disp[rank2]+3] = i;
        send_disp[rank2] += 4;
      }
    }
    intVec.clear();
    
    // rewind...
    
    send_disp[0] = 0;
    for (int rank = 1; rank < mpi_size; ++rank)
      send_disp[rank] = send_count[rank-1] + send_disp[rank-1];
    
    MPI_Alltoall(send_count,1,MPI_INT,recv_count,1,MPI_INT,mpi_comm);

    recv_disp[0] = 0;
    for (int rank = 1; rank < mpi_size; ++rank)
      recv_disp[rank] = recv_count[rank-1] + recv_disp[rank-1];
    const int recv_count_sum = recv_disp[mpi_size-1] + recv_count[mpi_size-1];
    
    assert(recv_buf_int == NULL); recv_buf_int = new int[recv_count_sum];
    
    MPI_Alltoallv(send_buf_int,send_count,send_disp,MPI_INT,
                  recv_buf_int,recv_count,recv_disp,MPI_INT,mpi_comm);
    delete[] send_buf_int; send_buf_int = NULL;
    delete[] recv_disp; 
    delete[] recv_count; 
    delete[] send_disp; 
    delete[] send_count; 
    
    // allocate and initialize teost...
    
    CTI_Mmap_rw(teost,nst);
    
    // initialize the shared memory range that our rank will set...

    for (int ist = ist_begin; ist < ist_end; ++ist) 
      FOR_I3 teost[ist][i] = packTeost(ist,0,0);
    
    // no need to barrier here. We are only writing our part of teost...
    // unpack into our local nst's...
   
    for (int irecv = 0; irecv < recv_count_sum; irecv += 4) {
      const int ist =       recv_buf_int[irecv]; assert((ist >= ist_begin)&&(ist < ist_end)); // should be local
      const int i_and_bit = recv_buf_int[irecv+1];
      const int i = (i_and_bit&3); assert((i >= 0)&&(i < 3));
      const int bit = (i_and_bit>>2); assert((bit >= 0)&&(bit <= 2)); // for now 0 or 1/2
      const int ist2 =      recv_buf_int[irecv+2]; assert((ist2 >= 0)&&(ist2 < nst));
      const int j =         recv_buf_int[irecv+3]; assert((j >= 0)&&(j < 3));
      assert(ist != ist2);
      int ist_check,i_check,bit_check;
      unpackTeost(ist_check,i_check,bit_check,teost[ist][i]);
      assert(ist_check == ist);
      assert(i_check == 0);
      assert(bit_check == 0);
      teost[ist][i] = packTeost(ist2,j,bit);
      // check if ist2 is also local. If it is, set the teost in ist2 to point to ist...
      if ((ist2 >= ist_begin)&&(ist2 < ist_end)) {
        unpackTeost(ist_check,i_check,bit_check,teost[ist2][j]);
        assert(ist_check == ist2);
        assert(i_check == 0);
        assert(bit_check == 0);
        // this will have to be smarter with multiple periodicity...
        if (bit == 0) {
          teost[ist2][j] = packTeost(ist,i,0);
        }
        else if (bit == 1) {
          teost[ist2][j] = packTeost(ist,i,2);
        }
        else {
          assert(bit == 2);
          teost[ist2][j] = packTeost(ist,i,1);
        }
      }
    }
    
    delete[] recv_buf_int; recv_buf_int = NULL;   

    MPI_Barrier(mpi_comm_shared);
    
    if ((mpi_rank_shared == 0)&&(mpi_size_internode > 1)) {
      
      for (int rank_internode = 0; rank_internode < mpi_size_internode; ++rank_internode) {
        
        // Bcast the data that the ranks on this node wrote into teost to all other nodes...
        
        const int rank = rank_of_rank_internode[rank_internode];
        const int ist_begin = min(nst,rank*nst_local_target);
        const int rank_p1 = ( rank_internode < mpi_size_internode-1 ? rank_of_rank_internode[rank_internode+1] : mpi_size );
        const int ist_end = min(nst,rank_p1*nst_local_target);

        MPI_Bcast((uint8*)(teost+ist_begin),(ist_end-ist_begin)*3,MPI_UINT8,rank_internode,mpi_comm_internode); 
        
      }

    }

    MPI_Barrier(mpi_comm_shared);

    /*
    // check teost...
    {
    double my_d2_max = 0.0;
    for (int ist = 0; ist < nst; ++ist) {
    const int isz = znost[ist];
    assert((isz >= 0)&&(isz < zoneVec.size()));
    if (zoneVec[isz].isBoundary()) {
    FOR_I3 {
    int ist_nbr,i_nbr,bit;
    unpackTeost(ist_nbr,i_nbr,bit,teost[ist][i]);
    assert(ist_nbr != ist);
    assert((i_nbr >= 0)&&(i_nbr < 3)); 
    // our 2 xsp's on this edge are...
    const int isp0 = spost[ist][i];
    const int isp1 = spost[ist][(i+1)%3];
    const int isp0_nbr = spost[ist_nbr][i_nbr];
    const int isp1_nbr = spost[ist_nbr][(i_nbr+1)%3];
    if (bit == 0) {
    // this should be a direct connection to ist_nbr...
    assert(isp0 == isp1_nbr);
    assert(isp1 == isp0_nbr);
    }
    else {
    // we should have pbi set in some of these nodes...
    if (pbi[isp0] == uint8(isp0)) {
    assert(pbi[isp1] == uint8(isp1));
    assert(pbi[isp0_nbr] != uint8(isp0_nbr)); 
    assert(pbi[isp1_nbr] != uint8(isp1_nbr)); 
    int isp0_nbr_bits,isp0_nbr_parent;
    BitUtils::unpackPbiHash(isp0_nbr_bits,isp0_nbr_parent,pbi[isp0_nbr]);
    int isp1_nbr_bits,isp1_nbr_parent;
    BitUtils::unpackPbiHash(isp1_nbr_bits,isp1_nbr_parent,pbi[isp1_nbr]);
    assert(isp0 == isp1_nbr_parent);
    assert(isp1 == isp0_nbr_parent);
    }
    else {
    assert(pbi[isp1] != uint8(isp1));
    assert(pbi[isp0_nbr] == uint8(isp0_nbr)); 
    assert(pbi[isp1_nbr] == uint8(isp1_nbr)); 
    int isp0_bits,isp0_parent;
    BitUtils::unpackPbiHash(isp0_bits,isp0_parent,pbi[isp0]);
    int isp1_bits,isp1_parent;
    BitUtils::unpackPbiHash(isp1_bits,isp1_parent,pbi[isp1]);
    assert(isp0_parent == isp1_nbr);
    assert(isp1_parent == isp0_nbr);
    }
    // also, this illustrates how to translate the nbr to match our 
    // coordinates...
    double xsp0_nbr_t[3]; FOR_J3 xsp0_nbr_t[j] = xsp[isp0_nbr][j];
    PeriodicData::periodicTranslate(xsp0_nbr_t,1,bit);
    double xsp1_nbr_t[3]; FOR_J3 xsp1_nbr_t[j] = xsp[isp1_nbr][j];
    PeriodicData::periodicTranslate(xsp1_nbr_t,1,bit);
    const double d2a = DIST2(xsp[isp0],xsp1_nbr_t);
    const double d2b = DIST2(xsp[isp1],xsp0_nbr_t);
    my_d2_max = max(my_d2_max,max(d2a,d2b));
    }
    }
    }
    else {
    // for periodic, should just be self...
    FOR_I3 {
    int ist_nbr,i_nbr,bit;
    unpackTeost(ist_nbr,i_nbr,bit,teost[ist][i]);
    assert(ist_nbr == ist);
    assert(i_nbr == 0);
    assert(bit == 0);
    }
    }
    }
    double d2_max;
    MPI_Reduce(&my_d2_max,&d2_max,1,MPI_DOUBLE,MPI_MAX,0,mpi_comm);
    if (mpi_rank == 0)
    cout << " > periodic node match distance (should be zero): " << sqrt(my_d2_max) << endl;
    }
    */

    if (mpi_rank == 0)
      cout << " > teost built successfully" << endl;

    // END = new delete buildTeost XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    
  }

  void ensureFlagSt() {
    if (flag_st == NULL) CTI_Mmap_rw(flag_st,nst);
  }
  
  void flagDisjointGroups(const set<int>& zones) {

    // START = new delete flagDisjointGroups XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    
    if (mpi_rank == 0) {
      cout << "SurfaceShm::flagDisjointGroups()" << endl;
      for (set<int>::iterator it = zones.begin(); it != zones.end(); ++it) {
        const int isz = *it;
        assert((isz >= 0)&&(isz < zoneVec.size()));
        cout << " > surface zone " << isz << " \"" << zoneVec[isz].getName() << "\"" << endl;
      }
    }

    if (zones == group_zones) {
      if (mpi_rank == 0) 
        cout << " > ngr_global: " << ngr_global << " groups already built for requested zones" << endl;
      return;
    }
    
    // clear any current group data, and start to build new group data
    // as part of this routine...

    group_zones = zones;
    ngr = 0;
    ngr_global = 0;
    spogr_i.clear(); // index of first isp in group.
    spogr_i.push_back(0);
    spogr_v.clear(); // value, i.e. isp
    nlogr.clear(); // loop and periodic loop count for each group...
    
    // we need teost...

    ensureTeost();
    
    // and a quick way to check if an ist is to be considered... 

    const int nsz = zoneVec.size();
    int * sz_flag = new int[nsz];
    for (int isz = 0; isz < nsz; ++isz)
      sz_flag[isz] = 0;
    
    for (set<int>::iterator it = zones.begin(); it != zones.end(); ++it) {
      const int isz = *it;
      assert((isz >= 0)&&(isz < nsz));
      sz_flag[isz] = 1;
    }

    const int nst_local_target = (nst + mpi_size - 1)/mpi_size;
    int ist_begin = min(nst,mpi_rank*nst_local_target);
    int ist_end = min(nst,(mpi_rank+1)*nst_local_target);
    int nst_local = ist_end-ist_begin;

    // load balance the relevant tris geometrically...
    
    int my_nst_flagged = 0; 
    for (int ist = ist_begin; ist < ist_end; ++ist) {
      const int isz = znost[ist];
      assert((isz >= 0)&&(isz < nsz));
      if (sz_flag[isz])
	++my_nst_flagged;
    }

    // report...
    // for small counts, it may be better to NOT load balance on points.
    // just not worth it...

    {
      int nst_flagged;
      MPI_Reduce(&my_nst_flagged,&nst_flagged,1,MPI_INT,MPI_SUM,0,mpi_comm);
      if (mpi_rank == 0)
        cout << " > considering " << nst_flagged << " tris..." << endl;
    }
    
    double (*xst_flagged)[3] = new double[my_nst_flagged][3];
    int8 *ist_flagged = new int8[my_nst_flagged];

    my_nst_flagged = 0; 
    for (int ist = ist_begin; ist < ist_end; ++ist) {
      const int isz = znost[ist];
      assert((isz >= 0)&&(isz < nsz));
      if (sz_flag[isz]) {
        FOR_I3 {
          xst_flagged[my_nst_flagged][i] = 0.0;
          FOR_J3 xst_flagged[my_nst_flagged][i] += xsp[spost[ist][j]][i];
        }
        ist_flagged[my_nst_flagged] = ist;
	++my_nst_flagged;
      }
    }

    MpiStuff::repartXcvPadt(xst_flagged,ist_flagged,my_nst_flagged,mpi_comm);

    delete[] xst_flagged;

    // we need a shared-memory tri flag...
    
    ensureFlagSt();

    if (mpi_rank_shared == 0)
      for (int ist = 0; ist < nst; ++ist)
        flag_st[ist] = -1;
    MPI_Barrier(mpi_comm_shared);

    // now put a positive rank in our local tris...
    
    for (int ii = 0; ii < my_nst_flagged; ++ii) {
      const int ist = ist_flagged[ii];
      assert((ist >= 0)&&(ist < nst));
      assert(flag_st[ist] == -1);
      flag_st[ist] = mpi_rank;
    }

    // include a second int here to account for the 
    // bits of the tris, although this is ignored for now...

    set<pair<int,int> > istBitsSet;
    stack<pair<int,int> > istBitsStack;
    vector<pair<int,int> > linkVec;
    vector<pair<int,int> > periodicLinkVec;
    map<const int,int> spLocalMap;
      
    for (int ii = 0; ii < my_nst_flagged; ++ii) {
      const int ist_seed = ist_flagged[ii];
      const int isz_seed = znost[ist_seed];
      assert((isz_seed >= 0)&&(isz_seed < nsz));
      assert(sz_flag[isz_seed]); // should be the case
      if (flag_st[ist_seed] == mpi_rank) {

	const int igroup = ngr++;
        
        // we decide who owns a particular group by who owns the 
        // maximum tri index in the group...
	int ist_max = ist_seed;
        bool ist_max_is_local = true;
        flag_st[ist_seed] = -mpi_rank-2; // flip to -2-indexed rank

	// this is a valid seed...
	assert(istBitsSet.empty());
	assert(istBitsStack.empty());
	assert(linkVec.empty());
	assert(periodicLinkVec.empty());
	istBitsSet.insert(pair<int,int>(ist_seed,0));
	istBitsStack.push(pair<int,int>(ist_seed,0));
	while (!istBitsStack.empty()) {
	  pair<int,int> ist_bits_pair = istBitsStack.top(); istBitsStack.pop();
	  const int ist = ist_bits_pair.first;
	  const int bits = ist_bits_pair.second;
	  // loop our nbrs and add them AND stack them if they are not already in the set...
	  FOR_I3 {
	    int ist_nbr,i_nbr,bit_nbr;
	    unpackTeost(ist_nbr,i_nbr,bit_nbr,teost[ist][i]);
	    assert(ist_nbr != ist); // will fail if user has asked to group a periodic surface -- figure this out later
	    // if this tri nbr is in the zone(s) being considered...
	    const int isz_nbr = znost[ist_nbr];
	    assert((isz_nbr >= 0)&&(isz_nbr < nsz));
	    if (sz_flag[isz_nbr]) {
              if (bit_nbr != 0) {
                // this nbr is across a periodic boundary. push both links into 
                // periodicLinkVec...
                periodicLinkVec.push_back(pair<int,int>(spost[ist][i],spost[ist][(i+1)%3]));
              }
	      if (flag_st[ist_nbr] == mpi_rank) {
                istBitsSet.insert(pair<int,int>(ist_nbr,0));
                istBitsStack.push(pair<int,int>(ist_nbr,0));
                flag_st[ist_nbr] = -mpi_rank-2;
                if (ist_nbr > ist_max) {
                  ist_max = ist_nbr;
                  ist_max_is_local = true;
                }
	      }
	      else if (flag_st[ist_nbr] != -mpi_rank-2) {
                if (istBitsSet.find(pair<int,int>(ist_nbr,0)) == istBitsSet.end()) {
                  istBitsSet.insert(pair<int,int>(ist_nbr,0));
                  istBitsStack.push(pair<int,int>(ist_nbr,0));
                  if (ist_nbr > ist_max) {
                    ist_max = ist_nbr;
                    ist_max_is_local = false;
                  }
                }
	      }
	    }
	    else {
	      // this edge is an open link associated with this disjoint group...
	      linkVec.push_back(pair<int,int>(spost[ist][i],spost[ist][(i+1)%3]));
	    }
	  }
	}

	if (ist_max_is_local) {

          // we are responsible for this group...
          
	  // decide how many loops there are in the linkVec...
          assert(spLocalMap.empty());
	  int nsp_local = 0;
	  for (int ii = 0; ii < linkVec.size(); ++ii) {
	    const int isp0 = linkVec[ii].first;
	    map<const int,int>::iterator it = spLocalMap.find(isp0);
	    if (it == spLocalMap.end()) {
	      spLocalMap[isp0] = linkVec[ii].first = nsp_local++;
	    }
	    else {
	      linkVec[ii].first = it->second;
	    }
	    const int isp1 = linkVec[ii].second;
	    it = spLocalMap.find(isp1);
	    if (it == spLocalMap.end()) {
	      spLocalMap[isp1] = linkVec[ii].second = nsp_local++;
	    }
	    else {
	      linkVec[ii].second = it->second;
	    }
	  }
          
	  int * sposp_local = new int[nsp_local];
	  for (map<const int,int>::iterator it = spLocalMap.begin(); it != spLocalMap.end(); ++it) 
	    sposp_local[it->second] = it->first;
	  spLocalMap.clear();
          
	  int * sp_next = new int[nsp_local];
	  for (int isp_local = 0; isp_local < nsp_local; ++isp_local)
	    sp_next[isp_local] = 0;
          
	  for (int ii = 0; ii < linkVec.size(); ++ii) {
	    const int isp0_local = linkVec[ii].first;
	    ++sp_next[isp0_local];
	    const int isp1_local = linkVec[ii].second;
	    ++sp_next[isp1_local];
	  }
          
	  for (int isp_local = 0; isp_local < nsp_local; ++isp_local) {
	    if (sp_next[isp_local] != 2) {
	      cout << "Complex link connectivity!" << endl;
	      assert(0);
	    }
	    sp_next[isp_local] = -1;
	  }
          
	  for (int ii = 0; ii < linkVec.size(); ++ii) {
	    const int isp0_local = linkVec[ii].first;
	    const int isp1_local = linkVec[ii].second;
	    assert(sp_next[isp0_local] == -1); // will fail for places where more than 2 links share the same node
	    sp_next[isp0_local] = isp1_local;
	  }

	  for (int isp_local = 0; isp_local < nsp_local; ++isp_local) 
	    assert(sp_next[isp_local] >= 0);
          
	  int nloop = 0;
	  for (int isp_local = 0; isp_local < nsp_local; ++isp_local) {
	    if (sp_next[isp_local] >= 0) {
	      const int iloop = nloop++;
	      const int isp_local_start = isp_local;
	      int isp_local_next = sp_next[isp_local_start];
	      sp_next[isp_local_start] = -2;
	      int isp_start = sposp_local[isp_local_start];
              spogr_v.push_back(isp_start);
	      int isp_next = sposp_local[isp_local_next];
              spogr_v.push_back(isp_next);
	      double length = DIST(xsp[isp_start],xsp[isp_next]);
	      double xc[3]; FOR_I3 xc[i] = length*(xsp[isp_start][i]+xsp[isp_next][i]);
	      double normal[3] = { 0.0, 0.0, 0.0 };
	      while (isp_local_next != isp_local_start) {
		const int isp_local_prev = isp_local_next;
		const int isp_prev = isp_next;
		isp_local_next = sp_next[isp_local_next];
		assert((isp_local_next >= 0)&&(isp_local_next < nsp_local));
		sp_next[isp_local_prev] = -2;
		isp_next = sposp_local[isp_local_next];
                spogr_v.push_back(isp_next);
		const double dlength = DIST(xsp[isp_prev],xsp[isp_next]);
		FOR_I3 xc[i] += dlength*(xsp[isp_prev][i]+xsp[isp_next][i]);
		length += dlength;
		const double n2[3] = TRI_NORMAL_2(xsp[isp_start],xsp[isp_prev],xsp[isp_next]);
		FOR_I3 normal[i] += n2[i];
	      }
	      FOR_I3 xc[i] /= length*2.0;
	      FOR_I3 normal[i] *= 0.5;
	      const double area = MAG(normal);
              //cout << " > iloop: " << iloop << " length: " << length << " xc: " << COUT_VEC(xc) << " normal: " << COUT_VEC(normal) << " area: " << area << endl;
	    }
	  }
	  for (int isp_local = 0; isp_local < nsp_local; ++isp_local) 
	    assert(sp_next[isp_local] == -2);

	  delete[] sposp_local; sposp_local = NULL;
	  delete[] sp_next; sp_next = NULL;

          // ----------------------------------------------
          // handle any periodic loops if present...
          // ----------------------------------------------
          
	  int nloop_periodic = 0;
          if (!periodicLinkVec.empty()) {

            if ((mpi_rank == 120) && (igroup == 1)) GeomUtils::writeTecplotEdges("edges_periodic.dat",xsp,periodicLinkVec);

            // decide how many loops there are in the linkVec...
            assert(spLocalMap.empty());
            int nsp_local = 0;
            for (int ii = 0; ii < periodicLinkVec.size(); ++ii) {
              const int isp0 = periodicLinkVec[ii].first;
              map<const int,int>::iterator it = spLocalMap.find(isp0);
              if (it == spLocalMap.end()) {
                spLocalMap[isp0] = periodicLinkVec[ii].first = nsp_local++;
              }
              else {
                periodicLinkVec[ii].first = it->second;
              }
              const int isp1 = periodicLinkVec[ii].second;
              it = spLocalMap.find(isp1);
              if (it == spLocalMap.end()) {
                spLocalMap[isp1] = periodicLinkVec[ii].second = nsp_local++;
              }
              else {
                periodicLinkVec[ii].second = it->second;
              }
            }
          
            assert(sposp_local == NULL); sposp_local = new int[nsp_local];
            for (map<const int,int>::iterator it = spLocalMap.begin(); it != spLocalMap.end(); ++it) 
              sposp_local[it->second] = it->first;
            spLocalMap.clear();
          
            assert(sp_next == NULL); sp_next = new int[nsp_local];
            for (int isp_local = 0; isp_local < nsp_local; ++isp_local)
              sp_next[isp_local] = 0;
          
            for (int ii = 0; ii < periodicLinkVec.size(); ++ii) {
              const int isp0_local = periodicLinkVec[ii].first;
              ++sp_next[isp0_local];
              const int isp1_local = periodicLinkVec[ii].second;
              ++sp_next[isp1_local];
            }
          
            for (int isp_local = 0; isp_local < nsp_local; ++isp_local) {
              if (sp_next[isp_local] != 2) {
                cout << "mpi_rank, igroup: " << mpi_rank << " " << igroup << " --Complex periodic link connectivity!" << endl;
                assert(0);
              }
              sp_next[isp_local] = -1;
            }
            
            for (int ii = 0; ii < periodicLinkVec.size(); ++ii) {
              const int isp0_local = periodicLinkVec[ii].first;
              const int isp1_local = periodicLinkVec[ii].second;
              assert(sp_next[isp0_local] == -1); // will fail for places where more than 2 links share the same node
              sp_next[isp0_local] = isp1_local;
            }

            for (int isp_local = 0; isp_local < nsp_local; ++isp_local) 
              assert(sp_next[isp_local] >= 0);
          
            assert(nloop_periodic == 0); // set outside loop
            for (int isp_local = 0; isp_local < nsp_local; ++isp_local) {
              if (sp_next[isp_local] >= 0) {
                const int iloop = nloop_periodic++;
                const int isp_local_start = isp_local;
                int isp_local_next = sp_next[isp_local_start];
                sp_next[isp_local_start] = -2;
                int isp_start = sposp_local[isp_local_start];
                spogr_v.push_back(isp_start);
                int isp_next = sposp_local[isp_local_next];
                spogr_v.push_back(isp_next);
                double length = DIST(xsp[isp_start],xsp[isp_next]);
                double xc[3]; FOR_I3 xc[i] = length*(xsp[isp_start][i]+xsp[isp_next][i]);
                double normal[3] = { 0.0, 0.0, 0.0 };
                while (isp_local_next != isp_local_start) {
                  const int isp_local_prev = isp_local_next;
                  const int isp_prev = isp_next;
                  isp_local_next = sp_next[isp_local_next];
                  assert((isp_local_next >= 0)&&(isp_local_next < nsp_local));
                  sp_next[isp_local_prev] = -2;
                  isp_next = sposp_local[isp_local_next];
                  spogr_v.push_back(isp_next);
                  const double dlength = DIST(xsp[isp_prev],xsp[isp_next]);
                  FOR_I3 xc[i] += dlength*(xsp[isp_prev][i]+xsp[isp_next][i]);
                  length += dlength;
                  const double n2[3] = TRI_NORMAL_2(xsp[isp_start],xsp[isp_prev],xsp[isp_next]);
                  FOR_I3 normal[i] += n2[i];
                }
                FOR_I3 xc[i] /= length*2.0;
                FOR_I3 normal[i] *= 0.5;
                const double area = MAG(normal);
                //cout << " > iloop: " << iloop << " length: " << length << " xc: " << COUT_VEC(xc) << " normal: " << COUT_VEC(normal) << " area: " << area << endl;
              }
            }
            for (int isp_local = 0; isp_local < nsp_local; ++isp_local) 
              assert(sp_next[isp_local] == -2);
           
            delete[] sposp_local; sposp_local = NULL;
            delete[] sp_next; sp_next = NULL;
            
          }
          
          nlogr.push_back(pair<int,int>(nloop,nloop_periodic));
          spogr_i.push_back(spogr_v.size());
          
	}
	else {
          // rewind: this group belonds to another rank...
	  --ngr;
	}
          
	//GeomUtils::writeTecplot("subzone.dat",spost,nst,xsp,istBitsSet);
	//GeomUtils::writeTecplotEdges("edges.dat",xsp,linkVec);

	istBitsSet.clear();
	linkVec.clear();
        periodicLinkVec.clear();
          
      }
    }
    
    assert(nlogr.size() == ngr);
    
    delete[] ist_flagged;
    delete[] sz_flag;
    
    int my_buf[4] = { ngr, 0, 0, 0 };
    for (int igr = 0; igr < ngr; ++igr) {
      if ((nlogr[igr].first == 2)&&(nlogr[igr].second == 0)) {
        ++my_buf[1];
      }
      else if ((nlogr[igr].first == 2)&&(nlogr[igr].second == 2)) {
        ++my_buf[2];
      }
      else {
        ++my_buf[3];
      }
    }

    int buf[4];
    MPI_Reduce(my_buf,buf,4,MPI_INT,MPI_SUM,0,mpi_comm);
    
    // set ngr_global...
    if (mpi_rank == 0) ngr_global = buf[0];
    MPI_Bcast(&ngr_global,1,MPI_INT,0,mpi_comm);
    
    if (mpi_rank == 0) {
      cout << " > got " << ngr_global << " disjoint groups." << endl;
      cout << "   > " << buf[1] << " groups with 2 links, no periodicity" << endl;
      cout << "   > " << buf[2] << " groups with 2 links, with periodicity" << endl;
      cout << "   > " << buf[3] << " groups with other link counts" << endl;
    }
    
    // END = new delete flagDisjointGroups XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  }
  
  void applyInterface(const string& zoneNamesCsv,const double axis[3],const double xc[3]) {
    
    ensureTeost();
    flagZonesMatching(zoneNamesCsv);
    
    double axis_mag = MAG(axis);
    double unit_axis[3];
    FOR_I3 unit_axis[i] = axis[i]/axis_mag;

    if (mpi_rank_shared == 0) {
      // step 1: look for open edge loops...
      int * sp_flag = new int[nsp];
      double (*dx_sp)[3] = new double[nsp][3];
      for (int isp = 0; isp < nsp; ++isp) {
	sp_flag[isp] = -1;
	FOR_I3 dx_sp[isp][i] = 0.0;
      }
      for (int ist = 0; ist < nst; ++ist) {
	if (zoneVec[znost[ist]].flag == 1) {
	  FOR_I3 {
	    int ist_nbr,i_nbr,bit_nbr;
	    unpackTeost(ist_nbr,i_nbr,bit_nbr,teost[ist][i]);
	    assert(ist_nbr != ist);
	    if (zoneVec[znost[ist_nbr]].flag == 0) {
	      // this is an edge of the flagged region: make sure it is as we expect...
	      assert(spost[ist][i] == spost[ist_nbr][(i_nbr+1)%3]);
	      assert(spost[ist][(i+1)%3] == spost[ist_nbr][i_nbr]);
	      assert(sp_flag[spost[ist][i]] == -1);
	      sp_flag[spost[ist][i]] = spost[ist][(i+1)%3];
	      //cout << "LOOKS GOOD!" << endl;
	    }
	  }
	}
      }
      for (int isp = 0; isp < nsp; ++isp) {
	if (sp_flag[isp] >= 0) {
	  // this is the start of a loop...
	  double xsum = 0.0;
	  double rsum = 0.0;
	  double wsum = 0.0;
	  double xcsum[3] = { 0.0, 0.0, 0.0 };
	  double nsum[3] = { 0.0, 0.0, 0.0 };
	  const int isp_start = isp;
	  double dx[3] = DIFF(xsp[isp_start],xc);
	  double x_start = DOT_PRODUCT(dx,unit_axis);
	  double r_start = sqrt(max(0.0,DOT_PRODUCT(dx,dx)-x_start*x_start));
	  int isp_next = isp_start;
	  double x_next = x_start;
	  double r_next = r_start;
	  int count = 0;
	  while (1) {
	    ++count;
	    const int isp_prev = isp_next;
	    const double x_prev = x_next;
	    const double r_prev = r_next;
	    isp_next = sp_flag[isp_prev];
	    assert(isp_next >= 0);
	    double dx[3] = DIFF(xsp[isp_next],xc);
	    x_next = DOT_PRODUCT(dx,unit_axis);
	    r_next = sqrt(max(0.0,DOT_PRODUCT(dx,dx)-x_next*x_next));
	    const double wgt = DIST(xsp[isp_prev],xsp[isp_next]);
	    xsum += wgt*(x_prev+x_next);
	    rsum += wgt*(r_prev+r_next);
	    wsum += wgt;
	    // and the center and normal...
	    FOR_I3 xcsum[i] += wgt*(xsp[isp_prev][i]+xsp[isp_next][i]);
	    double n2[3] = TRI_NORMAL_2(xsp[isp_start],xsp[isp_prev],xsp[isp_next]);
	    FOR_I3 nsum[i] += n2[i];
	    if (isp_next == isp_start)
	      break;
	  }
	  xsum /= 2.0*wsum;
	  rsum /= 2.0*wsum;
	  FOR_I3 xcsum[i] /= 2.0*wsum;
	  double nmag = MAG(nsum);
	  FOR_I3 nsum[i] /= nmag;
	  // loop again, setting dx_sp to snap to exact x,r...
	  double d2_max = 0.0;
	  double x_max[3];
	  isp_next = isp_start;
	  while (1) {
	    const int isp_prev = isp_next;
	    isp_next = sp_flag[isp_prev];
	    assert(isp_next >= 0);
	    sp_flag[isp_prev] = -2;
	    double dx[3] = DIFF(xsp[isp_next],xc);
	    x_next = DOT_PRODUCT(dx,unit_axis);
	    r_next = sqrt(max(0.0,DOT_PRODUCT(dx,dx)-x_next*x_next));
	    double dx_next = x_next-xsum;
	    double dr_next = r_next-rsum;
	    // the unit_axis gives the dx correction...
	    FOR_I3 dx_sp[isp_next][i] = -dx_next*unit_axis[i];
	    // the r is...
	    //double dr[3]; FOR_I3 dr[i] = dx[i] - x_next*unit_axis[i];
	    FOR_I3 dx_sp[isp_next][i] -= dr_next*(dx[i] - x_next*unit_axis[i])/r_next;
	    const double d2 = DOT_PRODUCT(dx_sp[isp],dx_sp[isp]);
	    if (d2 > d2_max) {
	      d2_max = d2;
	      FOR_I3 x_max[i] = xsp[isp][i];
	    }
	    if (isp_next == isp_start)
	      break;
	  }
	  if (mpi_rank == 0) cout << " > open edge loop with " << count << 
			       " edges: avg r: " << rsum << " n: " << COUT_VEC(nsum) << " xc: " << COUT_VEC(xcsum) << 
			       " max dx: " << sqrt(d2_max) << " at x: " << COUT_VEC(x_max) << endl;
	}
      }
      // now look for axis points: i.e. points with high st complexity...
      for (int isp = 0; isp < nsp; ++isp) sp_flag[isp] = 0;
      for (int ist = 0; ist < nst; ++ist) {
	if (zoneVec[znost[ist]].flag == 1) {
	  FOR_I3 {
	    const int isp = spost[ist][i];
	    ++sp_flag[isp];
	  }
	}
      }
      for (int isp = 0; isp < nsp; ++isp) {
	if (sp_flag[isp] >= 32) {
	  double dx[3] = DIFF(xsp[isp],xc);
	  double x = DOT_PRODUCT(dx,unit_axis);
	  double dr = sqrt(max(0.0,DOT_PRODUCT(dx,dx)-x*x));
	  FOR_I3 dx_sp[isp][i] = -(dx[i] - x*unit_axis[i]);
	  if (mpi_rank == 0) cout << " > got axis sp with " << sp_flag[isp] << " st: max dx: " << MAG(dx_sp[isp]) << " at x: " << COUT_VEC(xsp[isp]) << endl;
	}
      }
      // now we need to smooth dx_sp...
      // ...
      // ...
      // ...
      // then add it to xsp...
      for (int isp = 0; isp < nsp; ++isp) {
	FOR_I3 xsp[isp][i] += dx_sp[isp][i];
      }
      delete[] dx_sp;
      delete[] sp_flag; 
    }
    MPI_Barrier(mpi_comm_shared);
  }
  
};

#endif
