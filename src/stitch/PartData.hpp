#ifndef _PART_DATA_HPP_
#define _PART_DATA_HPP_

#include "CTI.hpp"
using namespace CTI;

#include "SurfaceShm.hpp"
#include "Points.hpp"
#include "Adt2d.hpp"
#include "StZone.hpp"
#include "PeriodicData.hpp"
#include "VoxelGrid.hpp"
#include "WebUI.hpp"
#include "MinHeap.hpp"
#include "IntFlag.hpp"
#include "SplineStuff.hpp"
#include "SubSurface.hpp"

namespace PartData {

  // part containment is done using an integer space that requires
  // these double-to-integer conversion params...
  extern bool b_dtoi;
  extern double dtoi_x0[3];
  extern double dtoi_delta;

  void setDtoiStuff();
  int getOddInt(const double dx,const double delta);
  int getEvenInt(const double dx,const double delta);

#define LEVEL_MASK_BITS 31  // max level < 2^5-1 (reserve 31 for part points)

  // TODO: eliminate partkind and use FarFieldType instead: this
  // is what kind is used for anyways...

  enum FarFieldType {
    UNKNOWN_FF,
    NO_FF,
    ANNULAR_CYLINDER_FF,
    CYLINDER_FF,
    SURFACE_FAZONE_FF,     // ff checks use one or more zones from the part's surface (NOT ff_surface)...
    FF_SURFACE_FF,         // ff uses ff_surface
    FF_SURFACE_FAZONE_FF,  // a part like a strand mesh that has both its own ff_surface AND flagged tris on some other surface
    SPHERE_FF,
    BOX_FF,
    //SIMPLE_DISK_Z_FF,
    //SIMPLE_BOX_FF,
    //PART_FF,
  };

  enum SolidType {
    UNKNOWN_SOLID,
    NO_SOLID,
    SURFACE_SOLID, // inside solid checks use the part's surface (unless they are flagged ff)...
    BOX_SOLID,
    ANNULAR_SOLID,
    CYLINDER_SOLID,
    SPHERE_SOLID,
  };

  bool testPointValidForTransform(const double xp[3],const int bits);
    
  class Part {
  private:

    uint * surface_st_bits;

  public:

    // this stores the (eventual) index of this part in the partVec. It is used by some routines before
    // the part is pushed into partVec...
    int ipart_of_part;

    bool b_name;
    string name;

    SurfaceShm * surface;
    int * surface_local_sp_flag;
    int * surface_local_st_flag;
    int * znosz;
    int * surface_zn_flag;

    // the far-field surface of a part can come from either the following
    // "ff_surface" which can store a variable ff_surface_dxost at all of its
    // tris to help stitch match the resolution of nearby hcp points, AND/OR
    // from zones in the "surface" flagged with surface_zn_flag allocated and
    // == -1

    SurfaceShm * ff_surface;
    double * ff_surface_dxost;
    int ff_level_max;
    int ff_nlayers_default;
    string ff_nlayersString;
    int ff_nlayers[LEVEL_MASK_BITS];

    Points * pts;

    // other part data...

    FarFieldType ff_type;
    double ddata_ff_[16];
    double e0_ff[3];
    double e1_ff[3];
    double e2_ff[3];
    Adt2d<int> * adt2d_ff;
    vector<int> stobb_ff; // surface-tri-of-bbox
    int bbmin_adt2d_ff[2];
    int bbmax_adt2d_ff[2];

    SolidType solid_type;
    double ddata_solid_[16];
    double e0_solid[3];
    double e1_solid[3];
    double e2_solid[3];
    Adt2d<int> * adt2d_solid;
    vector<int> stobb_solid; // surface-tri-of-bbox
    int bbmin_adt2d_solid[2];
    int bbmax_adt2d_solid[2];

    bool b_dx;
    double dx[3];

    // there may be some surfaces linked to the part
    // allow it to associate with multiple parts and surfaces
    bool b_link_surface;
    vector<int> ipart_link;
    vector<int> ist_offset_link;

    Adt2d<int> * adt2d_ff_surface;
    vector<int> stobb_ff_surface; // surface-tri-of-bbox
    int bbmin_adt2d_ff_surface[2];
    int bbmax_adt2d_ff_surface[2];

    // For WRITE_PART (moving solver)
    // flag surfaces or ff surfaces in some parts if we don't want
    // to write them into the part file
    bool b_write_surface;
    bool b_write_ff;

    // stuff from moving renamed as vd for now....
    Adt<double> * x_vd_bbox_adt;
    Adt<double> * x_vd_adt;
    
    //SurfaceShm * surface;
    SubSurface * subSurface; // also from moving

    Part() {

      b_name = false;

      surface = NULL;
      surface_st_bits = NULL;
      surface_local_sp_flag = NULL;
      surface_local_st_flag = NULL;
      znosz = NULL;
      surface_zn_flag = NULL;

      ff_surface = NULL;
      ff_surface_dxost = NULL;
      ff_nlayers_default = 10;
      ff_level_max = -1;

      pts = NULL;

      ff_type = UNKNOWN_FF;
      e0_ff[0] = 1.0; e0_ff[1] = 0.0; e0_ff[2] = 0.0;
      e1_ff[0] = 0.0; e1_ff[1] = 1.0; e1_ff[2] = 0.0;
      e2_ff[0] = 0.0; e2_ff[1] = 0.0; e2_ff[2] = 1.0;
      adt2d_ff = NULL;

      b_link_surface = false;
      adt2d_ff_surface = NULL;

      solid_type = UNKNOWN_SOLID;
      e0_solid[0] = 1.0; e0_solid[1] = 0.0; e0_solid[2] = 0.0;
      e1_solid[0] = 0.0; e1_solid[1] = 1.0; e1_solid[2] = 0.0;
      e2_solid[0] = 0.0; e2_solid[1] = 0.0; e2_solid[2] = 1.0;
      adt2d_solid = NULL;

      b_dx = false;

      b_write_surface = true;
      b_write_ff = true;
      
      x_vd_bbox_adt = NULL;
      x_vd_adt = NULL;
    
      subSurface = NULL;
      
    }

    void clear() {

      if (surface_st_bits) {
        CTI_Munmap(surface_st_bits,surface->nst);
        surface_st_bits = NULL;
      }
      if (surface) {
        delete surface;
        surface = NULL;
      }
      DELETE(znosz);
      DELETE(surface_zn_flag);
      DELETE(surface_local_sp_flag);
      DELETE(surface_local_st_flag);
      if (ff_surface_dxost) {
        CTI_Munmap(ff_surface_dxost,ff_surface->nst);
        ff_surface_dxost = NULL;
      }
      if (ff_surface) {
        delete ff_surface;
        ff_surface = NULL;
      }
      if (pts) {
        delete pts;
        pts = NULL;
      }
      if (adt2d_ff) {
        delete adt2d_ff;
        adt2d_ff = NULL;
      }
      if (adt2d_ff_surface) {
        delete adt2d_ff_surface;
        adt2d_ff_surface = NULL;
      }
      if (adt2d_solid) {
        delete adt2d_solid;
        adt2d_solid = NULL;
      }
      
      if (x_vd_bbox_adt) {
        delete x_vd_bbox_adt;
        x_vd_bbox_adt = NULL;
      }
      if (x_vd_adt) {
        delete x_vd_adt;
        x_vd_adt = NULL;
      }
      if (subSurface) {
        delete subSurface;
        subSurface = NULL;
      }

      // new solver data clearing
      

    }

    ~Part() {
      clear();
    }


    string getName() const {
      if (b_name) return name;
      else return "UNNAMED";
    }

    string getFfTypeName() const {
      switch (ff_type) {
      case UNKNOWN_FF:
        return "UNKNOWN_FF";
      case NO_FF:
        return "NO_FF";
      case ANNULAR_CYLINDER_FF:
        return "ANNULAR_CYLINDER_FF";
      case CYLINDER_FF:
        return "CYLINDER_FF";
      case SURFACE_FAZONE_FF:
        return "SURFACE_FAZONE_FF";
      case FF_SURFACE_FF:
        return "FF_SURFACE_FF";
      case FF_SURFACE_FAZONE_FF:
        return "FF_SURFACE_FAZONE_FF";
      case SPHERE_FF:
        return "SPHERE_FF";
      case BOX_FF:
        return "BOX_FF";
      default:
        return "UNRECOGNIZED FF TYPE: update PartData::Part::getFfTypeName()";
      }
    }

    string getSolidTypeName() const {
      switch (solid_type) {
      case UNKNOWN_SOLID:
        return "UNKNOWN_SOLID";
      case NO_SOLID:
        return "NO_SOLID";
      case SURFACE_SOLID:
        return "SURFACE_SOLID";
      case ANNULAR_SOLID:
        return "ANNULAR_SOLID";
      case BOX_SOLID:
        return "BOX_SOLID";
      case CYLINDER_SOLID:
        return "CYLINDER_SOLID";
      case SPHERE_SOLID:
        return "SPHERE_SOLID";
      default:
        return "UNRECOGNIZED SOLID TYPE: update PartData::Part::getSolidTypeName()";
      }
    }

    bool hasSimpleFF() const {
      assert(ff_type != UNKNOWN_FF);
      return !((ff_type == SURFACE_FAZONE_FF)||(ff_type == FF_SURFACE_FF)||(ff_type == FF_SURFACE_FAZONE_FF));
    }

    bool hasFF() const {
      assert(ff_type != UNKNOWN_FF);
      return ff_type != NO_FF;
    }

    bool hasSimpleSolid() const {
      assert(solid_type != UNKNOWN_SOLID);
      return !(solid_type == SURFACE_SOLID);
    }

    void resize_surface(const int nsp_new,const int nst_new) {
      assert(surface);
      if (surface_st_bits) {
        assert(nst_new > surface->nst);
        uint * surface_st_bits_new = NULL;
        CTI_Mmap(surface_st_bits_new,nst_new);
        if (mpi_rank_shared == 0) {
          for (int ist = 0; ist < surface->nst; ++ist)
            surface_st_bits_new[ist] = surface_st_bits[ist];
          for (int ist = surface->nst; ist < nst_new; ++ist)
            surface_st_bits_new[ist] = 0;
        }
        MPI_Barrier(mpi_comm_shared);
        CTI_Munmap(surface_st_bits,surface->nst);
        surface_st_bits = surface_st_bits_new;
      }
      DELETE(znosz); // gets reallocated
      assert(surface_zn_flag == NULL); //DELETE(surface_zn_flag);
      assert(surface_local_sp_flag == NULL); //DELETE(surface_local_sp_flag);
      assert(surface_local_st_flag == NULL); //DELETE(surface_local_st_flag);
      // and call the surface resize too...
      surface->resize(nsp_new,nst_new);
    }

    // functions associated with surface_st_bits...

    void init_surface_st_bits() {
      assert(surface_st_bits == NULL);
      assert(surface);
      CTI_Mmap(surface_st_bits,surface->nst);
      if (mpi_rank_shared == 0) {
        for (int ist = 0; ist < surface->nst; ++ist)
          surface_st_bits[ist] = 0;
      }
      MPI_Barrier(mpi_comm_shared);
    }

    void setGroupForSt(const int igr,const int ist) {
      assert(mpi_rank_shared == 0);
      assert(surface_st_bits);
      int igr_check;
      assert(!getGroupForSt(igr_check,ist));
      assert((igr >= 0)&&(igr < 65535));
      surface_st_bits[ist] |= (uint(igr+1)<<16); // 1-indexing of group
      // check...
      const bool works = getGroupForSt(igr_check,ist);
      assert(works);
      assert(igr_check == igr);
    }

    void clearGroupForSt(const int ist) {
      assert(mpi_rank_shared == 0);
      assert(surface_st_bits);
      surface_st_bits[ist] &= 65535u;
    }

    bool getGroupForSt(int& igr,const int ist) const {
      assert(surface_st_bits);
      igr = int(surface_st_bits[ist]>>16)-1;
      if (igr >= 0)
        return true;
      return false;
    }

    void setPartForSt(const int ipart,const int ist) {
      assert(mpi_rank_shared == 0);
      assert(surface_st_bits);
      int ipart_check;
      assert((!getPartForSt(ipart_check,ist))||(ipart_check == ipart));
      //assert((ipart >= 0)&&(ipart < 31));
      //surface_st_bits[ist] |= (uint(ipart+1)<<11); // 1-indexing of part
      assert((ipart >= 0)&&(ipart < 63));
      surface_st_bits[ist] |= (uint(ipart+1)<<10); // 1-indexing of part
      // check...
      const bool works = getPartForSt(ipart_check,ist);
      assert(works);
      assert(ipart_check == ipart);
    }

    bool getPartForSt(int& ipart,const int ist) const {
      assert(surface_st_bits);
      //ipart = int((surface_st_bits[ist]&65535u)>>11)-1;
      ipart = int((surface_st_bits[ist]&65535u)>>10)-1;
      if (ipart >= 0)
        return true;
      return false;
    }

    void clearPartForSt(const int ist) {
      assert(mpi_rank_shared == 0);
      assert(surface_st_bits);
      //surface_st_bits[ist] &= ~(31u<<11);
      surface_st_bits[ist] &= ~(63u<<10);
    }

    void clearSurfaceStBitsBitsForPart(const int ipart_clear) {
      if (mpi_rank_shared == 0) {
        for (int ist = 0; ist < surface->nst; ++ist) {
          int ipart;
          if (getPartForSt(ipart,ist)) {
            if (ipart == ipart_clear)
              clearPartForSt(ist);
          }
        }
      }
      MPI_Barrier(mpi_comm_shared);
    }

    void setWindowForSt(const int iwd,const int ist) {
      assert(mpi_rank_shared == 0);
      assert(surface_st_bits);
      //assert((iwd >= 0)&&(iwd < 2047)); // 2047 = 2^11-1
      //surface_st_bits[ist] &= ~2047u;
      assert((iwd >= 0)&&(iwd < 1023)); // 1023 = 2^10-1
      surface_st_bits[ist] &= ~1023u;
      surface_st_bits[ist] |= uint(iwd+1); // 1-indexing of the window
      // check...
      int iwd_check;
      const bool works = getWindowForSt(iwd_check,ist);
      assert(works);
      assert(iwd_check == iwd);
    }

    bool getWindowForSt(int& iwd,const int ist) const {
      assert(surface_st_bits);
      //iwd = int(surface_st_bits[ist]&2047u)-1;
      iwd = int(surface_st_bits[ist]&1023u)-1;
      if (iwd >= 0)
        return true;
      return false;
    }

    void buildOrRebuildSubSurface(const double bbmin[3],const double bbmax[3]) {

      assert(surface);
      if (subSurface == NULL) {
        subSurface = new SubSurface();
      }
      else {
        subSurface->clear();
	// rebuilding the subsurface is expensive and means it was not big enough the first time. If we
	// hit this regularly, we should consider increasing the initial tolerances on the subsurface 
	// bbox...
        cout << "[" << mpi_rank << "] WARNING: buildOrRebuildSubSurface: " << COUT_VEC(bbmin) << " " << COUT_VEC(bbmax) << endl;
      }

      // set the bbmin/bbmax...
      FOR_I3 subSurface->bbmin[i] = bbmin[i];
      FOR_I3 subSurface->bbmax[i] = bbmax[i];

      vector<int> bitVec;
      PeriodicData::buildBitVec(bitVec);

      vector<pair<int,int> > stVec; // ist, bits
      double xsp_t[3][3]; 
      for (int ist = 0; ist < surface->nst; ++ist) {
	// skip farfield surfaces: these are flagged as -1 in surface_zn_flag...
	if (surface_zn_flag) {
	  const int izn = surface->znost[ist];
	  assert((izn >= 0)&&(izn < surface->zoneVec.size()));
	  if (surface_zn_flag[izn] == -1)
	    continue;
	}
	// skip periodic...
	if ((surface->znost[ist] >= 0)&&(surface->zoneVec[surface->znost[ist]].isBoundary())) { 
	  for (int ib = 0, nib = bitVec.size(); ib < nib; ++ib) {
	    const int bits = bitVec[ib];
	    double xsp_avg[3] = { 0.0, 0.0, 0.0 }; 
            FOR_I3 {
              const int isp = surface->spost[ist][i];
              FOR_J3 {
		xsp_t[i][j] = surface->xsp[surface->spost[ist][i]][j];
		xsp_avg[j] += surface->xsp[surface->spost[ist][i]][j]/3.0;
	      }
	    }
            if (bits) {
	      // to handle the pathological cases associated with 180 degree periodicity, and also limit
	      // the tris we check for, check that the xsp_avg is atleast on the correct side of the
	      // bbox centroid for consideration...
	      if (!testPointValidForTransform(xsp_avg,bits))
		continue;
	      // otherwise, transform the corner point...
	      PeriodicData::periodicTranslate(xsp_t,3,bits);
      	    }

            if (((xsp_t[0][0] < bbmin[0])&&(xsp_t[1][0] < bbmin[0])&&(xsp_t[2][0] < bbmin[0]))||
                ((xsp_t[0][0] > bbmax[0])&&(xsp_t[1][0] > bbmax[0])&&(xsp_t[2][0] > bbmax[0]))||
                ((xsp_t[0][1] < bbmin[1])&&(xsp_t[1][1] < bbmin[1])&&(xsp_t[2][1] < bbmin[1]))||
                ((xsp_t[0][1] > bbmax[1])&&(xsp_t[1][1] > bbmax[1])&&(xsp_t[2][1] > bbmax[1]))||
                ((xsp_t[0][2] < bbmin[2])&&(xsp_t[1][2] < bbmin[2])&&(xsp_t[2][2] < bbmin[2]))||
                ((xsp_t[0][2] > bbmax[2])&&(xsp_t[1][2] > bbmax[2])&&(xsp_t[2][2] > bbmax[2])))
              continue;
	    stVec.push_back(pair<int,int>(ist,bits));
	  }
	}
      }

      subSurface->nst = stVec.size();
      assert(subSurface->ist_global_and_bits == NULL);
      subSurface->ist_global_and_bits = new int8[subSurface->nst];
      assert(subSurface->spost == NULL);
      subSurface->spost = new int[subSurface->nst][3];
      
      map<const pair<int,int>,int> pbiMap;
      subSurface->nsp = 0;
      for (int ii = 0, nii = stVec.size(); ii < nii; ++ii) {
	const int ist = stVec[ii].first;
	const int bits = stVec[ii].second;
	subSurface->setIstGlobalAndBits(ist,bits,ii);
	FOR_I3 {
	  int isp = surface->spost[ist][i];
	  int bits_active,isp_active;
	  BitUtils::unpackPbiHash(bits_active,isp_active,surface->getPbi(isp));
	  int bits_plus_bits_active = 0;
	  if (bits || bits_active) bits_plus_bits_active = BitUtils::addPeriodicBits(bits,bits_active);
	  // TODO: retain the check to test 180 periodicity on complex case, then remove because
	  // other cases (e.g. vki case with one tri across span I think) trigger this...
	  //assert(bits_plus_bits_active < (1<<6)); // the elimination above should ensure we never get double periodicity...
	  // we have the unique representation of this isp...
	  map<const pair<int,int>,int>::iterator iter = pbiMap.find(pair<int,int>(bits_plus_bits_active,isp_active));
	  if (iter == pbiMap.end()) {
	    // not found, so add...
	    pbiMap[pair<int,int>(bits_plus_bits_active,isp_active)] = subSurface->spost[ii][i] = subSurface->nsp++;
	  }
	  else {
	    subSurface->spost[ii][i] = iter->second;
	  }
	}
      }
      stVec.clear();

      // subSurface nodes...
      assert(subSurface->xp == NULL);
      subSurface->xp = new double[subSurface->nsp][3];
      for (map<const pair<int,int>,int>::iterator iter = pbiMap.begin(); iter != pbiMap.end(); ++iter) {
	const int bits = iter->first.first;
	const int isp = iter->first.second;
	const int ii = iter->second;
	FOR_I3 subSurface->xp[ii][i] = surface->xsp[isp][i];
	if (bits) PeriodicData::periodicTranslate(subSurface->xp[ii],1,bits);
      }
      pbiMap.clear();
      
      // finally, build the subsurface adt...
      double (*bbmin_adt)[3] = new double[subSurface->nst][3];
      double (*bbmax_adt)[3] = new double[subSurface->nst][3];
      for (int ist_ss = 0; ist_ss < subSurface->nst; ++ist_ss) {
        FOR_I3 bbmin_adt[ist_ss][i] = min(subSurface->xp[subSurface->spost[ist_ss][0]][i],
            min(subSurface->xp[subSurface->spost[ist_ss][1]][i],
              subSurface->xp[subSurface->spost[ist_ss][2]][i]));
        FOR_I3 bbmax_adt[ist_ss][i] = max(subSurface->xp[subSurface->spost[ist_ss][0]][i],
            max(subSurface->xp[subSurface->spost[ist_ss][1]][i],
              subSurface->xp[subSurface->spost[ist_ss][2]][i]));
      }
      assert(subSurface->adt == NULL);
      subSurface->adt = new Adt<double>(subSurface->nst,bbmin_adt,bbmax_adt);
      delete[] bbmin_adt;
      delete[] bbmax_adt;

    }

    void buildOrRebuildSubSurfaceOld(const double bbmin[3],const double bbmax[3]) {

      assert(0);

      // note: this one allocated a flag for ALL surface->nsp and surface->nst, which
      // can run out of memory for large surfaces. Use the new version above. This also
      // does not contain the point elimination check...

      assert(surface);
      if (subSurface == NULL) {
        subSurface = new SubSurface();
      }
      else {
        subSurface->clear();
        //cout << "[" << mpi_rank << "] buildOrRebuildSubSurface: " << COUT_VEC(bbmin) << " " << COUT_VEC(bbmax) << endl;
      }

      // set the bbmin/bbmax...
      FOR_I3 subSurface->bbmin[i] = bbmin[i];
      FOR_I3 subSurface->bbmax[i] = bbmax[i];

      int *st_flag = new int[surface->nst];
      int *sp_flag = new int[surface->nsp];
      for (int isp = 0; isp < surface->nsp; ++isp) 
        sp_flag[isp] = 0;
      for (int ist = 0; ist < surface->nst; ++ist) {
        st_flag[ist] = 0;
        // skip farfield surfaces: these are flagged as -1 in surface_zn_flag...
        if (surface_zn_flag) {
          const int izn = surface->znost[ist];
          assert((izn >= 0)&&(izn < surface->zoneVec.size()));
          if (surface_zn_flag[izn] == -1)
            continue;
        }
        // skip periodic...
        if ((surface->znost[ist] >= 0)&&(surface->zoneVec[surface->znost[ist]].isBoundary())) { 
          st_flag[ist] = 1;
          FOR_I3 {
            const int isp = surface->spost[ist][i];
            sp_flag[isp] = 1;
          }
        }
      }

      // also set bit 2 of points that participate in periodic boundaries. These
      // points need to be mapped when the surface is built...

      if (surface->checkPbi()) {
        for (int isp = 0; isp < surface->nsp; ++isp) {
          if (sp_flag[isp] & 1) {
            if (surface->getPbi(isp) != uint8(isp)) {
              sp_flag[isp] |= 2;
              // also flag the referenced point...
              int bits_active,isp_active;
              BitUtils::unpackPbiHash(bits_active,isp_active,surface->getPbi(isp));
              assert((isp_active >= 0)&&(isp_active < surface->nsp));
              assert(surface->getPbi(isp_active) == uint(isp_active));
              sp_flag[isp_active] |= 2;
            }
          }
        }
      }

      vector<int> bitVec;
      PeriodicData::buildBitVec(bitVec);

      const int nb = bitVec.size();
      assert(nb <= 29); // should be 27 at most
      subSurface->nsp = 0;
      subSurface->nst = 0;
      for (int ib = 0; ib < nb; ++ib) {
        const int bits = bitVec[ib];
        const int this_bit = (1<<(ib+2)); assert(this_bit > 2); // cannot use 1 or 2...

        // look for tris that could cross this bbox...
        for (int ist = 0; ist < surface->nst; ++ist) {
          if (st_flag[ist] & 1) {
            // if all points of this tri are outside the bounding box in any dimension,
            // then we can skip that tri...
            double xsp_t[3][3]; 
            FOR_I3 {
              const int isp = surface->spost[ist][i]; assert(sp_flag[isp] & 1);
              FOR_J3 xsp_t[i][j] = surface->xsp[surface->spost[ist][i]][j];
            }
            if (bits) PeriodicData::periodicTranslate(xsp_t,3,bits);
            if (((xsp_t[0][0] < bbmin[0])&&(xsp_t[1][0] < bbmin[0])&&(xsp_t[2][0] < bbmin[0]))||
                ((xsp_t[0][0] > bbmax[0])&&(xsp_t[1][0] > bbmax[0])&&(xsp_t[2][0] > bbmax[0]))||
                ((xsp_t[0][1] < bbmin[1])&&(xsp_t[1][1] < bbmin[1])&&(xsp_t[2][1] < bbmin[1]))||
                ((xsp_t[0][1] > bbmax[1])&&(xsp_t[1][1] > bbmax[1])&&(xsp_t[2][1] > bbmax[1]))||
                ((xsp_t[0][2] < bbmin[2])&&(xsp_t[1][2] < bbmin[2])&&(xsp_t[2][2] < bbmin[2]))||
                ((xsp_t[0][2] > bbmax[2])&&(xsp_t[1][2] > bbmax[2])&&(xsp_t[2][2] > bbmax[2])))
              continue;
            // flag and count this tri...
            assert(!(st_flag[ist] & this_bit));
            st_flag[ist] |= this_bit;
            ++subSurface->nst;
            // flag and count the points...
            FOR_I3 {
              const int isp = surface->spost[ist][i]; 
              if (!(sp_flag[isp] & this_bit)) {
                sp_flag[isp] |= this_bit;
                ++subSurface->nsp;
              }
            }
          }
        }
      }

      assert(subSurface->xp == NULL);
      subSurface->xp = new double[subSurface->nsp][3];
      assert(subSurface->ist_global_and_bits == NULL);
      subSurface->ist_global_and_bits = new int8[subSurface->nst];
      assert(subSurface->spost == NULL);
      subSurface->spost = new int[subSurface->nst][3];

      // now consider each possible set of bits...
      // we will need a map to handle common points...

      int isp_ss = 0;
      int ist_ss = 0;
      int* sp_index = new int[surface->nsp];
      map<const pair<int,int>,int> pbiMap;
      for (int ib = 0; ib < nb; ++ib) {
        const int bits = bitVec[ib];
        const int this_bit = (1<<(ib+2)); assert(this_bit > 2); // cannot use 1 or 2...

        // nodes first, so we have the numbering...

        for (int isp = 0; isp < surface->nsp; ++isp) {
          if ((sp_flag[isp] & 1)&&(sp_flag[isp] & this_bit)) {
            // we want this node. It may, however, altready be indexed...
            if (sp_flag[isp] & 2) {
              assert(surface->checkPbi()); // must be the case -- the presence of pbi determines the 2 bit setting above.
              int bits_active,isp_active;
              BitUtils::unpackPbiHash(bits_active,isp_active,surface->getPbi(isp));
              int bits_plus_bits_active = 0;
              if (bits || bits_active) bits_plus_bits_active = BitUtils::addPeriodicBits(bits,bits_active);
              map<const pair<int,int>,int>::iterator iter = pbiMap.find(pair<int,int>(bits_plus_bits_active,isp_active));
              if (iter == pbiMap.end()) {
                // not found, so add...
                pbiMap[pair<int,int>(bits_plus_bits_active,isp_active)] = isp_ss;
                sp_index[isp] = isp_ss;
                FOR_I3 subSurface->xp[isp_ss][i] = surface->xsp[isp][i];
                if (bits) PeriodicData::periodicTranslate(subSurface->xp[isp_ss],1,bits); // note that we only use bits here to transform
                ++isp_ss;
              }
              else {
                // we have already added this node...
                sp_index[isp] = iter->second;
              }
            }
            else {
              // this node is always new, but there is no need to map it...
              sp_index[isp] = isp_ss;
              FOR_I3 subSurface->xp[isp_ss][i] = surface->xsp[isp][i];
              if (bits) PeriodicData::periodicTranslate(subSurface->xp[isp_ss],1,bits); // note that we only use bits here to transform
              ++isp_ss;
            }
          }
        }

        // and tris second -- the nodes should all be numbered...

        for (int ist = 0; ist < surface->nst; ++ist) {
          if ((st_flag[ist] & 1)&&(st_flag[ist] & this_bit)) {
            // we want this tri...
            FOR_I3 subSurface->spost[ist_ss][i] = sp_index[surface->spost[ist][i]];
            subSurface->setIstGlobalAndBits(ist,bits,ist_ss);
            ++ist_ss;
          }
        }

      }
      assert(isp_ss <= subSurface->nsp);
      subSurface->nsp = isp_ss;
      assert(ist_ss == subSurface->nst);
      delete[] sp_index;
      delete[] sp_flag;
      delete[] st_flag;
      pbiMap.clear();

      // finally, build the subsurface adt...
      double (*bbmin_adt)[3] = new double[subSurface->nst][3];
      double (*bbmax_adt)[3] = new double[subSurface->nst][3];
      for (int ist_ss = 0; ist_ss < subSurface->nst; ++ist_ss) {
        FOR_I3 bbmin_adt[ist_ss][i] = min(subSurface->xp[subSurface->spost[ist_ss][0]][i],
            min(subSurface->xp[subSurface->spost[ist_ss][1]][i],
              subSurface->xp[subSurface->spost[ist_ss][2]][i]));
        FOR_I3 bbmax_adt[ist_ss][i] = max(subSurface->xp[subSurface->spost[ist_ss][0]][i],
            max(subSurface->xp[subSurface->spost[ist_ss][1]][i],
              subSurface->xp[subSurface->spost[ist_ss][2]][i]));
      }
      assert(subSurface->adt == NULL);
      subSurface->adt = new Adt<double>(subSurface->nst,bbmin_adt,bbmax_adt);
      delete[] bbmin_adt;
      delete[] bbmax_adt;
    }

    int readBinary(const string& filename);
    int writeBinary(const string& filename);

    // isInside routines...
    void prepareIsInsideFF(const int xy_min[2],const int xy_max[2]);
    bool isInsideFF(const double xp[3]);
    void prepareIsInsideSolid(const int xy_min[2],const int xy_max[2]);
    bool isInsideSolid(const double xp[3],const bool first);

    // for the case of some parts, it is useful to constrain the
    // general smoothing to certain directions, or typically to zero...
    void constrainSmoothing(double dxp[3],const double xp[3]) const;

  };

  // the vector of parts. Note that each part can have
  // surface - the boundary tris
  // ff_surface - tris that blank the points
  // pts - part points: normally structured in some way
  extern vector<Part*> partVec;

  // zoneVec contains a global zone name, representing the
  // contiguous list of surface zones...
  extern vector<StZone> zoneVec;

  // hcp points: the ones generated by HcpPointBuilderNew...
  extern Points * hcpPts;
  extern bool b_vp; // hcpPts are set from voronoi points

  int getFlattenedSp(const int ipart,const int isp);
  int getFlattenedSt(const int ipart,const int ist);
  void getPartStFromFlattenedSt(int& ipart,int& ist,const int flattened_ist);
  void getPartSpFromFlattenedSp(int& ipart,int& isp,const int flattened_isp);
  int getFlattenedNsp();
  int getFlattenedNst();
  double getFlattenedVolume();
  void getFlattenedCentroid(double xc[3]);

  extern bool b_bbminmax;
  extern double bbminmax[6]; // stores bbox as xmin,ymin,zmin,-xmax,-ymax,-zmax

  void buildOrRebuildZoneVec();
  int getZoneIndex(const string& name);
  pair<int,int> getPartSurfaceZonePair(const int izn);
  pair<int,int> getPartSurfaceZonePair(const string& name);

  void clear();

  // global part bbox...

  void ensureBbminmax();
  double getBoundingBoxRmax();
  void getBoundingBoxCenter(double bBoxCenter[3]);
  void getBbox(double bbmin[3],double bbmax[3]);
  void getBboxForZoneids(double bbmin[3],double bbmax[3],const vector<int>& zoneidVec);
  //void getSubzoneBoundingBoxCenterAndDiagonal(double bbmin[3],double bbmax[3],double xbb[3],double &diag,const int * const zone_flag);

  // part output - compresses surfaces and points and writes out a part...

  void writePart(const string& filename);

  bool checkCoplanar(const int ipart0,const int ist0,const int bits0,
                     const int ipart1,const int ist1,const int bits1,const double degrees);

  void runDiagnostics(Param* param, const bool b_help);

  void prepareIsInside();

  void queryLengthScales();
  void queryParts();
  void queryZones();
  void queryZones(const string& names_cds); // cds == comma-delimited-string

  // formerly PartMaker, but now condensed under the PartData namespace...

  void processPart(list<Param>::iterator param,const bool b_help = false);
  void processPart(Param * param,const bool b_help = false);
  void helpProcessPart();

  void processInterface(Param *param,const bool b_help = false);

  //void rmPart(Param * param);

  // TODO: the return "bool" in these makePart* routines below is more current and preferred...

  bool makePart(Part * part,Param * param);
  int makePartFromSbin(Part * part,const string& filename);
  int makePartFromMles(Part * part,const string& filename);
  void writeUsedPoints(const string& filename);

  void queryZones(const set<int>& zoneIdSet);
  bool makeStrand(Part * part,Param * param,int& iarg);
  bool makeTriStrand(Part * part,Param * param,int& iarg);
  bool makeAirfoilStrand(Part * part,Param * param,int& iarg);
  bool makeSvpStrand(Part * part,Param * param,int& iarg);
  bool makeExtrudeSplineWithStrand(Part * part,Param * param,int& iarg);
  bool makeCylinderWithStrand(Part * part,Param * param,int& iarg);
  int makeBoxWithCartPts(Part * part,Param * param,int& iarg,const bool pts_only);
  int makeAnnularCylPts(Part * part,Param * param,int& iarg,const bool pts_only);
  int makeAnnularDlrtPts(Part * part,Param * param,int& iarg,const bool pts_only);
  int makeCylWithHexcorePts(Part * part,Param * param,int& iarg,const bool pts_only);
  int makeCylWithDlrtPts(Part * part,Param * param,int& iarg,const bool pts_only);
  int makeSphereWithBl(Part * part,Param * param,int& iarg);
  int makePoints1D(Part * part,Param * param,int& iarg);
  int makeMedialAxis(Part * part,Param * param,int& iarg);
  void calcClosestPassageForFlaggedTris(double xc[3],double nc[3],const int& ipart,const int* st_flag,const double dir[3]);
  void calcGeoDistFromPoint(const double x[3], const int& ipart, double* sp_buf);
  void calcGeoDistFromPoint(vector<int> ispVec, const int& ipart, double* sp_buf);
  void fmmLoop(MinHeap& trialHeap,int* flag_gp,double* dist_gp,const int* stogp_i,const int* stogp_v,const int* gposp,SurfaceShm * surface);
  void fmmLoop(const int * stosp_i, const int * stosp_v, MinHeap& trialHeap, IntFlag& sp_flag, double* sp_buf, SurfaceShm * surface);

  // move this to PartData, but it requires the geometric data processing
  // stuff in part maker for now...

  void processPistonFF(Param * param);

  void processWritePbin(Param * param,const bool b_help);
  void helpProcessWritePbin();
  void processWriteHcpPts(Param * param,const bool b_help);
  void helpProcessWriteHcpPts();

  void processWriteSbin(Param * param,const bool b_help);

  // helper functions for airfoil strand...
  int getOrientation(const int np, const double direction[3], const double (*xp)[3]);
  void resampleSpline(SplineStuff::CubicSpline &cspline, const int nn, const double dt, const double dn, const double fraction, int &nt, double * &t_local, double (* &xt)[3]);

  // ----------------------------------------------
  // start of namespace-based SolverData...
  // ----------------------------------------------

  extern int step;
  extern int check_interval;

  extern int8 ncv_global;
  extern int ncv,ncv_g;
  extern int * flag_cv;
  extern double (*x_vd)[3]; // the voronoi forming points, if available
  extern double *delta_vd; // the nbr sphere radius that guarantees all voronoi nbrs (and surface tri patches) are present 

  extern double * vol_cv;
  extern double (*x_cv)[3];
  
  extern vector<int> cvopa_i; // sometimes the cv arrays are arranged contiguous by part - any hcp points are in the final range

  extern vector<int> paozn; // part-of-zone - short vector that relates parts to flattened zones 
  extern vector<pair<int,int> > pzozn;  // part-isz-of-zone - short vector that relates parts to flattened zones - replace above 

  extern map<const uint8,int> rbiMap;
  extern vector<int8> rbi_g;

  void rebuildZoneVecStuff(); // copied over from Move -- reconcile with buildOrRebuildZoneVec above
  void clear_new();
  
  // moved over from VoronoiBuilder...
  void processLoadBalanceMode(Param *param,const bool b_help);
  void helpLoadBalanceMode();
  void clearPoints();
  void ensurePoints();

} // namespace PartData

#endif
