#ifndef _VOXEL_GRID_HPP_
#define _VOXEL_GRID_HPP_

#include "CTI.hpp"
using namespace CTI;

#include "SurfaceShm.hpp"
#include "SimplePoint.hpp"
#include "Adt.hpp"
#include "DistributedDataExchanger.hpp"

// size of (internal) block in each dimension, must be a power of 2 to use bit shifting
#define BLOCK_BITS  4
#define BLOCK_SIZE  (1<<BLOCK_BITS) 
#define BLOCK_MASK  (BLOCK_SIZE-1)
#define BLOCK_SIZE3 (BLOCK_SIZE*BLOCK_SIZE*BLOCK_SIZE)
// ghost halo...
// local  = (global&BLOCK_MASK)+N_GHOSTS
// block  = (global>>BLOCK_BITS)
// global = (block<<BLOCK_BITS)+local-N_GHOSTS;
#define N_GHOSTS 1
#define BLOCK_SIZE_G (BLOCK_SIZE+2*N_GHOSTS)
#define BLOCK_SIZE3_G (BLOCK_SIZE_G*BLOCK_SIZE_G*BLOCK_SIZE_G)
#define ZONE_MASK 65535

class VoxelGrid { 
public:

  class VoxelBlockData {
  public:
    double data[BLOCK_SIZE3_G];
    double data0[BLOCK_SIZE3_G];
    double data00[BLOCK_SIZE3_G];
    int zone[BLOCK_SIZE3_G];
    int flag[BLOCK_SIZE3_G];
    VoxelBlockData() {
      for (int ijk = 0; ijk < BLOCK_SIZE3_G; ++ijk) data[ijk] = HUGE_VAL; 
      for (int ijk = 0; ijk < BLOCK_SIZE3_G; ++ijk) data0[ijk] = HUGE_VAL; 
      for (int ijk = 0; ijk < BLOCK_SIZE3_G; ++ijk) data00[ijk] = HUGE_VAL;
      for (int ijk = 0; ijk < BLOCK_SIZE3_G; ++ijk) zone[ijk] = -1;
      for (int ijk = 0; ijk < BLOCK_SIZE3_G; ++ijk) flag[ijk] = 0; // no bits set 
    }
  };

  VoxelBlockData ** vbd;
  int *vbd_rank;
  vector<int> vbdVec;
    
  /* --------
   * |\      \
   * | \______\
   * | |      |
   * \ |      |
   *  \|______|
   * x0
   */

  double x0[3];     // front bottom left corner of image in world coords
  double tol;   // coincident point relative tolerance
  double dx,inv_dx; // width/ni,ni/width
  double width;     // image width in world coords
  int ni,nj,nk;     // image voxel size: width,height,depth
  int nI,nJ,nK;     // block size (image divided into BLOCK_SIZExBLOCK_SIZExBLOCK_SIZE blocks)
  // marching cube tables which are set outside class
  static const int edgeTable[256];
  static const int triTable[256][16];
  vector<StZone> zoneVec; // need to store zones from initialized surface

  enum UpdateBits {
    UPDATE_DATA   = 1,
    UPDATE_DATA0  = 2,
    UPDATE_DATA00 = 4,
    UPDATE_ZONE   = 8,
    UPDATE_FLAG   = 16,
    UPDATE_ALL    = 65535, // all bits
  };

  enum FlagBits {
    SURFACE_BIT  = 1,
  };

  VoxelGrid() {
    vbd = NULL;
    vbd_rank = NULL;
  }

  ~VoxelGrid() {
    if (vbd) {
      for (int IJK = 0; IJK < nI*nJ*nK; ++IJK) {
        if (vbd[IJK] != NULL) delete vbd[IJK];
      }
      delete[] vbd;
    }
    DELETE(vbd_rank);
  }


  void init(const double bbmin[3],const double bbmax[3],const int nx) {

    tol = 1.0E-12;

    // record x0 and width... 
    
    x0[0] = bbmin[0];
    x0[1] = bbmin[1];
    x0[2] = bbmin[2];
    width = bbmax[0]-bbmin[0];

    // and image size...

    ni = nx; assert((ni >= 1)&&(ni < 16000));
    
    dx     = width/double(ni);
    inv_dx = double(ni)/width;

    nj = int((bbmax[1]-bbmin[1])*inv_dx); assert((nj >= 1)&&(nj < 16000));
    nk = int((bbmax[2]-bbmin[2])*inv_dx); assert((nk >= 1)&&(nk < 16000));

    // and allocate the image voxel storage...

    nI = ((ni-1)>>BLOCK_BITS)+1; assert((nI >= 1)&&(nI < 2000));
    nJ = ((nj-1)>>BLOCK_BITS)+1; assert((nJ >= 1)&&(nJ < 2000));
    nK = ((nk-1)>>BLOCK_BITS)+1; assert((nK >= 1)&&(nK < 2000));
    // assume there is an integer number of blocks...
    assert((nI*nJ*nK >= 0)&&(nI*nJ*nK < TWO_BILLION));

    vbd = new VoxelBlockData*[nI*nJ*nK];
    vbd_rank = new int[nI*nJ*nK];
    for (int IJK = 0; IJK < nI*nJ*nK; ++IJK) {
      vbd[IJK] = NULL;
      vbd_rank[IJK] = -1;
    }

    if (mpi_rank == 0) {
      cout << " > initialized VoxelGrid with voxel size: " << ni << "x" << nj << "x" << nk << "=" << 
        int8(ni)*int8(nj)*int8(nk) << ", block size: " << nI << "x" << nJ << "x" << nK << "=" << 
        nI*nJ*nK << ", width: " << width << ", and halo bloat: " << 
        double(BLOCK_SIZE3_G)/double(BLOCK_SIZE3)*100.0 << "%" << endl;
    }
  }
  
  double getPointToEdgeSignedDist(const double xp[3],const double v0[3],const double v1[3],const double n0[3],const double n1[3]) {
    
    // returns the square of the euclidean distance from point xp to edge v0->v1.

    const double dx[3] = DIFF(v1,v0);
    const double dxp[3] = DIFF(xp,v0);
    const double dp = DOT_PRODUCT(dx,dxp);
    double xc[3],nc[3];
    if (dp <= 0.0) {
      // we are closest to the first point v0. Note this includes the case when dx is zero, and/or dxp is zero...
      FOR_I3 xc[i] = v0[i];
      FOR_I3 nc[i] = n0[i];
    }
    else {
      const double dx2 = DOT_PRODUCT(dx,dx);
      if (dp >= dx2) {
	// we are closest to the second point v1...
        FOR_I3 xc[i] = v1[i];
        FOR_I3 nc[i] = n1[i];
      }
      else {
        FOR_I3 xc[i] = dp/dx2*v1[i] + (1.0-dp/dx2)*v0[i];
        FOR_I3 nc[i] = dp/dx2*n1[i] + (1.0-dp/dx2)*n0[i];
      }
    } 
    //cout << "c " << COUT_VEC(nc) << endl;

    const double dxc[3] = DIFF(xp,xc);
    //return DOT_PRODUCT(dxc,nc);
    return MAG(dx)*SGN(DOT_PRODUCT(dxc,nc));

  }

  double getPointToTriSignedDist(const double xp[3],
      const double v0[3],const double v1[3],const double v2[3],
      const double n0[3],const double n1[3],const double n2[3]) {
    
    // this one robust to bad tris...
    // should eventually modify the below routine in the same way...

    // compute Facet Triangle Edges: T(s,t) = v0 + s*edge0 + t*edge1
    const double edge0[3] = DIFF(v1,v0);
    const double edge1[3] = DIFF(v2,v0);

    // compute coefficients of distance squared function, Q(s,t)
    const double qa = DOT_PRODUCT(edge0,edge0);
    const double qc = DOT_PRODUCT(edge1,edge1);

    double xc[3],nc[3];
    if ((qa <= 0.0)&&(qc <= 0.0)) {
      // this is a zero tri: all points must be the same, so return distance to the
      // first point...
      FOR_I3 xc[i] = v0[i];
      FOR_I3 nc[i] = n0[i];
      //cout << "b " << COUT_VEC(nc) << endl;
      const double dx[3] = DIFF(xp,xc);
      return MAG(dx)*SGN(DOT_PRODUCT(dx,nc));
    }

    const double qb = DOT_PRODUCT(edge0,edge1);

    // compute coordinates Qmin = Q(sbar/det,tbar/det)
    double det = qa*qc - qb*qb;
    if (det <= 0.0) {
      const double edge01[3] = DIFF(v2,v1);
      if (DOT_PRODUCT(edge01,edge01) > max(qa,qc)) {
        return getPointToEdgeSignedDist(xp,v1,v2,n1,n2);
      }
      else if (qa >= qc) {
        return getPointToEdgeSignedDist(xp,v0,v1,n0,n1);
      }
      else {
        return getPointToEdgeSignedDist(xp,v0,v2,n0,n2);
      }
    }

    double diff_v0xp[3];
    FOR_I3 diff_v0xp[i] = (v0[i]) - (xp[i]);
    const double qd = DOT_PRODUCT(edge0,diff_v0xp);
    const double qe = DOT_PRODUCT(edge1,diff_v0xp);
    //const double qf = DOT_PRODUCT(diff_v0xp,diff_v0xp);
    double sbar = qb*qe - qc*qd;
    double tbar = qb*qd - qa*qe;

    /***********************************
     * Identify region containing xp
     *        t
     *      \2|
     *       \|
     *        \
     *        |\
     *        | \
     *      3 | 0\  1
     *     ---|---\-----s
     *        |    \
     *      4 | 5   \  6
     ************************************/
    int region = -1;
    if (sbar + tbar <= det){
      if (sbar<0){
        if (tbar<0) region = 4;
        else region = 3;
      }
      else if (tbar<0) region = 5;
      else region = 0;
    }
    else{
      if (sbar<0) region = 2;
      else if (tbar<0) region = 6;
      else region = 1;
    }
    assert(region>-1);

    // facet coordinates that give min distance to xp
    // will be assigned to sbar, tbar

    double invDet, numer, denom, tmp0, tmp1;
    switch (region) {
    case 0:
      invDet = 1/det;
      sbar *= invDet;
      tbar *= invDet;
      break;
    case 1:
      numer = (qc+qe-qb-qd);
      if (numer<=0) sbar = 0;
      else {
	denom = qa - 2*qb + qc; // positive
	sbar = ( numer >= denom ? 1 : numer/denom );
      }
      tbar = 1 - sbar;
      break;
    case 2:
      tmp0 = qb + qd;
      tmp1 = qc + qe;
      if (tmp1 > tmp0){ // minimum on edge s+t = 1
	numer = tmp1 - tmp0;
	denom = qa - 2*qb + qc;
	sbar = ( numer >= denom ? 1 : numer/denom);
	tbar = 1-sbar;
      }
      else { // minimum on edge s = 0
	sbar = 0;
	tbar = ( tmp1 <= 0 ? 1 : ( qe >= 0 ? 0 : -qe/qc ) );
      }
      break;
    case 3:
      sbar = 0;
      tbar = ( qe >= 0 ? 0 : ( -qe >= qc ? 1 : -qe/qc ) );
      break;
    case 4:
      if (qe < 0){ // minimum on edge s = 0
	sbar = 0;
	tbar = ( qc <= -qe ? 1 : -qe/qc );
      }
      else { // minimum on edge t = 0
	tbar = 0;
	sbar = ( qa <= -qd ? 1 : ( qd >= 0 ? 0 : -qd/qa ) );
      }
      break;
    case 5:
      tbar = 0;
      sbar = ( qd >= 0 ? 0 : ( -qd >= qa ? 1 : -qd/qa ) );
      break;
    case 6:
      tmp0 = qa + qd;
      tmp1 = qb + qe;
      if (tmp0 > tmp1){ // minimum on edge s+t = 1
	numer = tmp0 - tmp1;
	denom = qa - 2*qb + qc;
	tbar = ( numer >= denom ? 1 : numer/denom );
	sbar = 1-tbar;
      }
      else{ // minimum on edge t = 0
	tbar = 0;
	sbar = ( tmp0 <= 0 ? 1 : ( qd >= 0 ? 0 : -qd/qa ) );
      }
      break;
    }

    // return point coordinates in physical space...
    FOR_I3 xc[i] = sbar*edge0[i] + tbar*edge1[i] + v0[i];

    // perform same interpolation for vertex normals...
    FOR_I3 nc[i] = sbar*n1[i] + tbar*n2[i] + (1.0-sbar-tbar)*n0[i];
    //cout << "a " <<  COUT_VEC(nc) << endl;

    const double dx[3] = DIFF(xp,xc);
    return MAG(dx)*SGN(DOT_PRODUCT(dx,nc));

  }

  void initVoxelDataFromSurface(SurfaceShm* surface,const double band_width_ = 0.0,const bool b_write = false) {

    double band_width;
    if (band_width_ == 0.0)
      band_width = 10.0*dx;
    else 
      band_width = band_width_;

    if (mpi_rank == 0) 
      cout << "initVoxelDataFromSurface()" << endl;

    // copy over zoneVec...
    assert(zoneVec.empty());
    zoneVec.resize(surface->zoneVec.size());
    for (int iz = 0; iz < surface->zoneVec.size(); ++iz)
      zoneVec[iz] = surface->zoneVec[iz];

    int my_nst_avg = surface->nst/mpi_size;
    if (surface->nst%mpi_size) ++my_nst_avg;
    const int ist0 = min(surface->nst,mpi_rank*my_nst_avg);
    const int ist1 = min(surface->nst,(mpi_rank+1)*my_nst_avg);
    assert(ist1-ist0 <= my_nst_avg);

    // compute point normals using angle weight surface normals...

    map<const int,int> pointMap;
    int nsp_tmp = 0;
    for (int ist = ist0; ist < ist1; ++ist) {
      FOR_I3 {
        const int isp = surface->spost[ist][i];
        map<const int,int>::iterator iter = pointMap.find(isp);
        if (iter == pointMap.end()) 
          pointMap[isp] = nsp_tmp++;
      }
    }

    int *isp_of_isp_tmp = new int[nsp_tmp];
    double (*n_sp_tmp)[3] = new double[nsp_tmp][3];
    for (int isp_tmp = 0; isp_tmp < nsp_tmp; ++isp_tmp) 
      FOR_I3 n_sp_tmp[isp_tmp][i] = 0.0;
    for (int ist = ist0; ist < ist1; ++ist) {
      double xsps[3][3];
      FOR_I3 FOR_J3 xsps[i][j] = surface->xsp[surface->spost[ist][i]][j];
      const double dx01[3] = DIFF(xsps[1],xsps[0]);
      const double dx02[3] = DIFF(xsps[2],xsps[0]);
      double n_st[3] = CROSS_PRODUCT(dx01,dx02);
      const double mag = MAG(n_st); assert(mag > 0.0);
      const double angle = atan2(mag,DOT_PRODUCT(dx01,dx02));
      FOR_I3 n_st[i] *= angle/mag;
      FOR_I3 {
        const int isp = surface->spost[ist][i];
        map<const int,int>::iterator iter = pointMap.find(isp);
        assert(iter != pointMap.end());
        const int isp_tmp = iter->second;
        isp_of_isp_tmp[isp_tmp] = isp;
        FOR_J3 n_sp_tmp[isp_tmp][j] += n_st[j];
      }
    }

    {
      int *spora = NULL;
      MiscUtils::calcUniformDist(spora,surface->nsp,mpi_size);
      DistributedDataExchanger dde(isp_of_isp_tmp,nsp_tmp,spora);
      double (*n_sp)[3] = new double[spora[mpi_rank+1]-spora[mpi_rank]][3];
      for (int isp = 0; isp < spora[mpi_rank+1]-spora[mpi_rank]; ++isp) 
        FOR_I3 n_sp[isp][i] = 0.0;
      dde.push(n_sp,n_sp_tmp,ADD_DATA);
      for (int isp = 0; isp < spora[mpi_rank+1]-spora[mpi_rank]; ++isp) {
        const double mag = MAG(n_sp[isp]); assert(mag > 0.0);
        FOR_I3 n_sp[isp][i] /= mag;
      }
      dde.pull(n_sp_tmp,n_sp);
      delete[] n_sp;
      delete[] spora;
    }
    delete[] isp_of_isp_tmp; 

    vector<pair<pair<int,int>,pair<double,int> > > blockLocalDistSignedZoneVec;
    //const double halfsize[3] = {0.5*dx,0.5*dx,0.5*dx};
    const double halfsize[3] = {dx*(1.0-tol),dx*(1.0-tol),dx*(1.0-tol)}; // lets double the box size to fix nbrs as well. 

    //const double halfsize[3] = {(0.5+tol)*dx,(0.5+tol)*dx,(0.5+tol)*dx};
    //const double halfsize[3] = {(1.0+tol)*dx,(1.0+tol)*dx,(1.0+tol)*dx};
    for (int ist = ist0; ist < ist1; ++ist) {

      // get bbox of tri and bloat it with prescribed buffer...
      int isps[3]; 
      FOR_I3 isps[i] = surface->spost[ist][i];
      double xsps[3][3];
      FOR_I3 FOR_J3 xsps[i][j] = surface->xsp[isps[i]][j];
      double n_sps[3][3];
      FOR_I3 {
        map<const int,int>::iterator iter = pointMap.find(isps[i]);
        assert(iter != pointMap.end());
        FOR_J3 n_sps[i][j] = n_sp_tmp[iter->second][j];
      }

      double bbmin[3] = {HUGE_VAL,HUGE_VAL,HUGE_VAL};
      double bbmax[3] = {-HUGE_VAL,-HUGE_VAL,-HUGE_VAL};
      FOR_I3 {
        bbmin[i] = MIN4(bbmin[i],xsps[0][i],xsps[1][i],xsps[2][i])-band_width;
        bbmax[i] = MAX4(bbmax[i],xsps[0][i],xsps[1][i],xsps[2][i])+band_width;
      }
      // convert to global integer coordinates
      int ijk_min[3],ijk_max[3];
      const int nijk[3] = {ni,nj,nk};
      FOR_I3 {
        ijk_min[i] = max(0,int(ceil((bbmin[i]-x0[i])*inv_dx)));
        ijk_max[i] = min(nijk[i]-1,int(floor((bbmax[i]-x0[i])*inv_dx)));
      }

      double n_st[3] = TRI_NORMAL_2(xsps[0],xsps[1],xsps[2]);
      const double mag = MAG(n_st); 
      FOR_L3 n_st[l] /= mag;
      double x[3];
      for (int i = ijk_min[0]; i < ijk_max[0]; ++i) {
        x[0] = i*dx+x0[0];
        for (int j = ijk_min[1]; j < ijk_max[1]; ++j) {
          x[1] = j*dx+x0[1];
          for (int k = ijk_min[2]; k < ijk_max[2]; ++k) {
            x[2] = k*dx+x0[2];
            const int IJK = ((k>>BLOCK_BITS)*nJ+(j>>BLOCK_BITS))*nI + (i>>BLOCK_BITS);
            assert((IJK >= 0)&&(IJK < nI*nJ*nK));

            // include ghost shift...
            const int ijk = (((k&BLOCK_MASK)+N_GHOSTS)*BLOCK_SIZE_G+((j&BLOCK_MASK)+N_GHOSTS))*BLOCK_SIZE_G+((i&BLOCK_MASK)+N_GHOSTS);
            assert((ijk >= N_GHOSTS)&&(ijk < BLOCK_SIZE3_G));

            // set vbd_rank[IJK] to -2 so that we know it will be activated...
            vbd_rank[IJK] = -2; 

            // get the signed distance and push into vec... 
            //cout << COUT_VEC(n_st) << " ";
            const double signed_dist = getPointToTriSignedDist(x,xsps[0],xsps[1],xsps[2],n_sps[0],n_sps[1],n_sps[2]);
            const double dist = fabs(signed_dist);

            //const double dist = sqrt(MiscUtils::getPointToTriDist2(x,xsps[0],xsps[1],xsps[2]));
            //const double dx0[3] = DIFF(xsps[0],x);
            
            int zone_sign_surface = surface->znost[ist];
            //if (DOT_PRODUCT(n_st,dx0) < 0.0)
            if (signed_dist > 0.0) 
              zone_sign_surface |= (1<<16); // zone limited to 65535
            if (triBoxOverlap(x,halfsize,xsps)) {
              zone_sign_surface |= (1<<17); 
              blockLocalDistSignedZoneVec.push_back(pair<pair<int,int>,pair<double,int> >
                  (pair<int,int>(IJK,ijk),pair<double,int>(dist,zone_sign_surface))); 
            }
          }
        }
      }

    }
    pointMap.clear();
    delete[] n_sp_tmp;
    
    MPI_Allreduce(MPI_IN_PLACE,vbd_rank,nI*nJ*nK,MPI_INT,MPI_MIN,mpi_comm);

    // at this point vbd_rank = -2 if it is going to be used. now perform a simple striping...
   
    int rank = 0;
    int8 nvoxels_global = 0;
    for (int IJK = 0; IJK < nI*nJ*nK; ++IJK) {
      if (vbd_rank[IJK] == -2) {
        nvoxels_global += BLOCK_SIZE3_G;
        if (rank == mpi_rank) {
          vbdVec.push_back(IJK);
          assert(vbd[IJK] == NULL); 
          vbd[IJK] = new VoxelBlockData();
        }
        vbd_rank[IJK] = rank++;
        if (rank == mpi_size) 
          rank = 0;
      }
    }

    int my_nblocks = vbdVec.size();
    int max_nblocks;
    MPI_Reduce(&my_nblocks,&max_nblocks,1,MPI_INT,MPI_MAX,0,mpi_comm);
    if (mpi_rank == 0) {
      int nblocks_global = nvoxels_global/BLOCK_SIZE3_G;
      cout << " > global active voxels (with ghosts): " << nvoxels_global << 
        ", global active blocks: " << nblocks_global << 
        ", load imbalance: " << max_nblocks << "/" << double(nblocks_global)/double(mpi_size) << 
        "=" << double(max_nblocks)*double(mpi_size)/double(nblocks_global) << endl;
    }

    // now sort the vector, and locally prune it...
    
    sort(blockLocalDistSignedZoneVec.begin(),blockLocalDistSignedZoneVec.end());

    int *send_count = new int[mpi_size];
    FOR_RANK send_count[rank] = 0;
    int *send_disp = new int[mpi_size];
    int *send_buf_int = NULL;
    double *send_buf_double = NULL;
    for (int iter = 0; iter < 2; ++iter) {
      int current_IJK = -1;
      int current_ijk = -1;
      for (int ii = 0, limit = blockLocalDistSignedZoneVec.size(); ii < limit; ++ii) {
        const int IJK = blockLocalDistSignedZoneVec[ii].first.first;
        const int ijk = blockLocalDistSignedZoneVec[ii].first.second;
        if ((IJK == current_IJK)&&(ijk == current_ijk)) 
          continue;
        
        // new voxel, which we know is the min distance based on sort so we will send it...
        const int rank = vbd_rank[IJK]; assert((rank >= 0)&&(rank < mpi_size));
        if (iter == 0) {
          ++send_count[rank];
        }
        else {
          send_buf_int[send_disp[rank]*3  ] = IJK;
          send_buf_int[send_disp[rank]*3+1] = ijk;
          send_buf_int[send_disp[rank]*3+2] = blockLocalDistSignedZoneVec[ii].second.second; // zone_sign_surface
          send_buf_double[send_disp[rank]] = blockLocalDistSignedZoneVec[ii].second.first;
          ++send_disp[rank];
        }

        // update indices...
        current_IJK = IJK;
        current_ijk = ijk;
      }

      send_disp[0] = 0;
      for (int rank = 1; rank < mpi_size; ++rank)
        send_disp[rank] = send_count[rank-1] + send_disp[rank-1];

      if (iter == 0) {
        const int send_count_sum = send_disp[mpi_size-1] + send_count[mpi_size-1];
        send_buf_int = new int[3*send_count_sum];
        send_buf_double = new double[send_count_sum];
      }

    }
    blockLocalDistSignedZoneVec.clear(); 

    // exchange...

    int *recv_count = new int[mpi_size];
    MPI_Alltoall(send_count,1,MPI_INT,recv_count,1,MPI_INT,mpi_comm);

    int *recv_disp = new int[mpi_size];
    recv_disp[0] = 0;
    for (int rank = 1; rank < mpi_size; ++rank)
      recv_disp[rank] = recv_count[rank-1] + recv_disp[rank-1];
    const int recv_count_sum = recv_disp[mpi_size-1] + recv_count[mpi_size-1];

    double * recv_buf_double = new double[recv_count_sum];
    MPI_Alltoallv(send_buf_double,send_count,send_disp,MPI_DOUBLE,
                  recv_buf_double,recv_count,recv_disp,MPI_DOUBLE,mpi_comm);
    delete[] send_buf_double; 

    FOR_RANK {
      send_count[rank] *= 3;
      send_disp[rank] *= 3;
      recv_count[rank] *= 3;
      recv_disp[rank] *= 3;
    }

    int * recv_buf_int = new int[3*recv_count_sum];
    MPI_Alltoallv(send_buf_int,send_count,send_disp,MPI_INT,
                  recv_buf_int,recv_count,recv_disp,MPI_INT,mpi_comm);
    delete[] send_buf_int; 
    delete[] send_count;
    delete[] send_disp;
    delete[] recv_count;
    delete[] recv_disp;
    
    // now set our data...
    
    for (int irecv = 0; irecv < recv_count_sum; ++irecv) {
      const int IJK = recv_buf_int[irecv*3  ]; assert(vbd_rank[IJK] == mpi_rank); // I own it!
      const int ijk = recv_buf_int[irecv*3+1];
      const int zone = (recv_buf_int[irecv*3+2]&65535);
      double signed_dist = recv_buf_double[irecv];
      if (recv_buf_int[irecv*3+2]&(1<<16))
        signed_dist *= -1.0;
      assert(vbd[IJK] != NULL);
      if (fabs(signed_dist) < fabs(vbd[IJK]->data[ijk])) {
        vbd[IJK]->data[ijk] = signed_dist;
        vbd[IJK]->zone[ijk] = zone;
        if (recv_buf_int[irecv*3+2]&(1<<17))
          vbd[IJK]->flag[ijk] |= SURFACE_BIT;
      }
    }
    delete[] recv_buf_int;
    delete[] recv_buf_double;

    updateVoxelData(UPDATE_DATA|UPDATE_ZONE|UPDATE_FLAG,true);

    // check...
    //updateVoxelDataReverse(UPDATE_DATA|UPDATE_ZONE|UPDATE_FLAG,true);
    if (b_write)
      writeVoxelDataToTecplot("init.dat");

  }

  void updateVoxelData(const int update_bits = UPDATE_DATA,const bool b_full_halo = true) {

    int int_width = 2; // IJK,ijk
    if (update_bits&UPDATE_ZONE)
      ++int_width;
    if (update_bits&UPDATE_FLAG)
      ++int_width;

    int double_width = 0;
    if (update_bits&UPDATE_DATA)
      ++double_width;
    if (update_bits&UPDATE_DATA0)
      ++double_width;
    if (update_bits&UPDATE_DATA00)
      ++double_width;

    // send your internal data to your nbrs ghosts if they exist...

    int* send_count = new int[mpi_size];
    int* send_disp  = new int[mpi_size];
    FOR_RANK send_count[rank] = 0;
    int* send_buf_int = NULL;
    double* send_buf_double = NULL;
    for (int iter = 0; iter < 2; ++iter) {
      for (int ivbd = 0, limit = vbdVec.size(); ivbd < limit; ++ivbd) {
        const int IJK = vbdVec[ivbd];
        assert(vbd[IJK] != NULL); 
        const int I = IJK%nI;
        const int J = (IJK/nI)%nJ;
        const int K = IJK/(nI*nJ);
        assert(IJK == (K*nJ+J)*nI+I);

        for (int kk = 0; kk < 3; ++kk) {
          for (int jj = 0; jj < 3; ++jj) {
            for (int ii = 0; ii < 3; ++ii) {
              if (!b_full_halo) {
                if ( (abs(ii-1)+abs(jj-1)+abs(kk-1)) != 1)
                  continue;
              }
              else {
                if ((ii == 1)&&(jj == 1)&&(kk == 1))
                  continue;
              }
              const int IJK_nbr = ((K+kk-1)*nJ+J+jj-1)*nI+I+ii-1;
              if ((IJK_nbr >= 0)&&(IJK_nbr < nI*nJ*nK)&&(vbd_rank[IJK_nbr] >= 0)) {
                const int rank = vbd_rank[IJK_nbr];
                if (iter == 0) {
                  int count = 1;
                  if (ii == 1) 
                    count *= BLOCK_SIZE;
                  else 
                    count *= N_GHOSTS;
                  if (jj == 1) 
                    count *= BLOCK_SIZE;
                  else 
                    count *= N_GHOSTS;
                  if (kk == 1) 
                    count *= BLOCK_SIZE;
                  else 
                    count *= N_GHOSTS;
                  send_count[rank] += count;
                }
                else {
                  int kmin,kmax,kshift,jmin,jmax,jshift,imin,imax,ishift;
                  if (ii == 0) {
                    imin = N_GHOSTS;
                    imax = 2*N_GHOSTS;
                    ishift = BLOCK_SIZE;
                  }
                  else if (ii == 1) {
                    imin = N_GHOSTS;
                    imax = BLOCK_SIZE_G-N_GHOSTS;
                    ishift = 0;
                  }
                  else {
                    assert(ii == 2);
                    imin = BLOCK_SIZE;
                    imax = N_GHOSTS+BLOCK_SIZE;
                    ishift = -BLOCK_SIZE;
                  }
                  if (jj == 0) {
                    jmin = N_GHOSTS;
                    jmax = 2*N_GHOSTS;
                    jshift = BLOCK_SIZE;
                  }
                  else if (jj == 1) {
                    jmin = N_GHOSTS;
                    jmax = BLOCK_SIZE_G-N_GHOSTS;
                    jshift = 0;
                  }
                  else {
                    assert(jj == 2);
                    jmin = BLOCK_SIZE;
                    jmax = N_GHOSTS+BLOCK_SIZE;
                    jshift = -BLOCK_SIZE;
                  }
                  if (kk == 0) {
                    kmin = N_GHOSTS;
                    kmax = 2*N_GHOSTS;
                    kshift = BLOCK_SIZE;
                  }
                  else if (kk == 1) {
                    kmin = N_GHOSTS;
                    kmax = BLOCK_SIZE_G-N_GHOSTS;
                    kshift = 0;
                  }
                  else {
                    assert(kk == 2);
                    kmin = BLOCK_SIZE;
                    kmax = N_GHOSTS+BLOCK_SIZE;
                    kshift = -BLOCK_SIZE;
                  }
                  for (int k = kmin; k < kmax; ++k) {
                    for (int j = jmin; j < jmax; ++j) {
                      for (int i = imin; i < imax; ++i) {
                        const int ijk = (k*BLOCK_SIZE_G+j)*BLOCK_SIZE_G+i;
                        const int ijk_nbr = ((k+kshift)*BLOCK_SIZE_G+j+jshift)*BLOCK_SIZE_G+i+ishift;
                        int disp = 0;
                        send_buf_int[send_disp[rank]*int_width+disp++] = IJK_nbr;
                        send_buf_int[send_disp[rank]*int_width+disp++] = ijk_nbr;
                        if (update_bits&UPDATE_ZONE) 
                          send_buf_int[send_disp[rank]*int_width+disp++] = vbd[IJK]->zone[ijk];
                        if (update_bits&UPDATE_FLAG) 
                          send_buf_int[send_disp[rank]*int_width+disp++] = vbd[IJK]->flag[ijk];
                        disp = 0;
                        if (update_bits&UPDATE_DATA) 
                          send_buf_double[send_disp[rank]*double_width+disp++] = vbd[IJK]->data[ijk];
                        if (update_bits&UPDATE_DATA0) 
                          send_buf_double[send_disp[rank]*double_width+disp++] = vbd[IJK]->data0[ijk];
                        if (update_bits&UPDATE_DATA00) 
                          send_buf_double[send_disp[rank]*double_width+disp++] = vbd[IJK]->data00[ijk];
                        ++send_disp[rank];
                      }
                    }
                  }
                }
              }
            }
          }
        }

      }

      send_disp[0] = 0;
      for (int rank = 1; rank < mpi_size; ++rank)
        send_disp[rank] = send_count[rank-1] + send_disp[rank-1];

      if (iter == 0) {
        const int send_count_sum = send_disp[mpi_size-1] + send_count[mpi_size-1];
        send_buf_int = new int[int_width*send_count_sum];
        send_buf_double = new double[double_width*send_count_sum];
      }

    }

    // exchange...

    int *recv_count = new int[mpi_size];
    MPI_Alltoall(send_count,1,MPI_INT,recv_count,1,MPI_INT,mpi_comm);

    int *recv_disp = new int[mpi_size];
    recv_disp[0] = 0;
    for (int rank = 1; rank < mpi_size; ++rank)
      recv_disp[rank] = recv_count[rank-1] + recv_disp[rank-1];
    const int recv_count_sum = recv_disp[mpi_size-1] + recv_count[mpi_size-1];

    FOR_RANK {
      send_count[rank] *= double_width;
      send_disp[rank] *= double_width;
      recv_count[rank] *= double_width;
      recv_disp[rank] *= double_width;
    }

    double * recv_buf_double = new double[double_width*recv_count_sum];
    MPI_Alltoallv(send_buf_double,send_count,send_disp,MPI_DOUBLE,
                  recv_buf_double,recv_count,recv_disp,MPI_DOUBLE,mpi_comm);
    delete[] send_buf_double; 

    FOR_RANK {
      send_count[rank] = (send_count[rank]/double_width)*int_width;;
      send_disp[rank] = (send_disp[rank]/double_width)*int_width;;
      recv_count[rank] = (recv_count[rank]/double_width)*int_width;;
      recv_disp[rank] = (recv_disp[rank]/double_width)*int_width;;
    }

    int * recv_buf_int = new int[int_width*recv_count_sum];
    MPI_Alltoallv(send_buf_int,send_count,send_disp,MPI_INT,
                  recv_buf_int,recv_count,recv_disp,MPI_INT,mpi_comm);
    delete[] send_buf_int; 
    delete[] send_count;
    delete[] send_disp;
    delete[] recv_count;
    delete[] recv_disp;
    
    // now set our data...

    for (int irecv = 0; irecv < recv_count_sum; ++irecv) {
      int disp = 0;
      const int IJK = recv_buf_int[irecv*int_width+disp++]; assert(vbd_rank[IJK] == mpi_rank);
      const int ijk = recv_buf_int[irecv*int_width+disp++];
      assert(vbd[IJK] != NULL);
      if (update_bits&UPDATE_ZONE) 
        vbd[IJK]->zone[ijk] = recv_buf_int[irecv*int_width+disp++];
      if (update_bits&UPDATE_FLAG) 
        vbd[IJK]->flag[ijk] = recv_buf_int[irecv*int_width+disp++];
      disp = 0;
      if (update_bits&UPDATE_DATA) 
        vbd[IJK]->data[ijk] = recv_buf_double[irecv*double_width+disp++];
      if (update_bits&UPDATE_DATA0) 
        vbd[IJK]->data0[ijk] = recv_buf_double[irecv*double_width+disp++];
      if (update_bits&UPDATE_DATA00) 
        vbd[IJK]->data00[ijk] = recv_buf_double[irecv*double_width+disp++];
    }
    delete[] recv_buf_int;
    delete[] recv_buf_double;

  }

  void updateVoxelDataReverse(const int update_bits = UPDATE_DATA,const bool b_full_halo = true) {

    int int_width = 2; // IJK,ijk
    if (update_bits&UPDATE_ZONE)
      ++int_width;
    if (update_bits&UPDATE_FLAG)
      ++int_width;

    int double_width = 0;
    if (update_bits&UPDATE_DATA)
      ++double_width;
    if (update_bits&UPDATE_DATA0)
      ++double_width;
    if (update_bits&UPDATE_DATA00)
      ++double_width;

    int* send_count = new int[mpi_size];
    int* send_disp  = new int[mpi_size];
    FOR_RANK send_count[rank] = 0;
    int* send_buf_int = NULL;
    double* send_buf_double = NULL;
    for (int iter = 0; iter < 2; ++iter) {
      for (int ivbd = 0, limit = vbdVec.size(); ivbd < limit; ++ivbd) {
        const int IJK = vbdVec[ivbd];
        assert(vbd[IJK] != NULL); 
        const int I = IJK%nI;
        const int J = (IJK/nI)%nJ;
        const int K = IJK/(nI*nJ);
        assert(IJK == (K*nJ+J)*nI+I);

        for (int kk = 0; kk < 3; ++kk) {
          for (int jj = 0; jj < 3; ++jj) {
            for (int ii = 0; ii < 3; ++ii) {
              if (!b_full_halo) {
                if ( (abs(ii-1)+abs(jj-1)+abs(kk-1)) != 1)
                  continue;
              }
              else {
                if ((ii == 1)&&(jj == 1)&&(kk == 1))
                  continue;
              }
              const int IJK_nbr = ((K+kk-1)*nJ+J+jj-1)*nI+I+ii-1;
              if ((IJK_nbr >= 0)&&(IJK_nbr < nI*nJ*nK)&&(vbd_rank[IJK_nbr] >= 0)) {
                const int rank = vbd_rank[IJK_nbr];
                if (iter == 0) {
                  int count = 1;
                  if (ii == 1) 
                    count *= BLOCK_SIZE;
                  else 
                    count *= N_GHOSTS;
                  if (jj == 1) 
                    count *= BLOCK_SIZE;
                  else 
                    count *= N_GHOSTS;
                  if (kk == 1) 
                    count *= BLOCK_SIZE;
                  else 
                    count *= N_GHOSTS;
                  send_count[rank] += count;
                }
                else {
                  int kmin,kmax,kshift,jmin,jmax,jshift,imin,imax,ishift;
                  if (ii == 0) {
                    imin = 0;
                    imax = N_GHOSTS;
                    ishift = BLOCK_SIZE;
                  }
                  else if (ii == 1) {
                    imin = N_GHOSTS;
                    imax = BLOCK_SIZE_G-N_GHOSTS;
                    ishift = 0;
                  }
                  else {
                    assert(ii == 2);
                    imin = N_GHOSTS+BLOCK_SIZE;
                    imax = BLOCK_SIZE_G; 
                    ishift = -BLOCK_SIZE;
                  }
                  if (jj == 0) {
                    jmin = 0;
                    jmax = N_GHOSTS;
                    jshift = BLOCK_SIZE;
                  }
                  else if (jj == 1) {
                    jmin = N_GHOSTS;
                    jmax = BLOCK_SIZE_G-N_GHOSTS;
                    jshift = 0;
                  }
                  else {
                    assert(jj == 2);
                    jmin = N_GHOSTS+BLOCK_SIZE;
                    jmax = BLOCK_SIZE_G;
                    jshift = -BLOCK_SIZE;
                  }
                  if (kk == 0) {
                    kmin = 0;
                    kmax = N_GHOSTS;
                    kshift = BLOCK_SIZE;
                  }
                  else if (kk == 1) {
                    kmin = N_GHOSTS;
                    kmax = BLOCK_SIZE_G-N_GHOSTS;
                    kshift = 0;
                  }
                  else {
                    assert(kk == 2);
                    kmin = N_GHOSTS+BLOCK_SIZE;
                    kmax = BLOCK_SIZE_G;
                    kshift = -BLOCK_SIZE;
                  }
                  for (int k = kmin; k < kmax; ++k) {
                    for (int j = jmin; j < jmax; ++j) {
                      for (int i = imin; i < imax; ++i) {
                        const int ijk = (k*BLOCK_SIZE_G+j)*BLOCK_SIZE_G+i;
                        const int ijk_nbr = ((k+kshift)*BLOCK_SIZE_G+j+jshift)*BLOCK_SIZE_G+i+ishift;
                        int disp = 0;
                        send_buf_int[send_disp[rank]*int_width+disp++] = IJK_nbr;
                        send_buf_int[send_disp[rank]*int_width+disp++] = ijk_nbr;
                        if (update_bits&UPDATE_ZONE) 
                          send_buf_int[send_disp[rank]*int_width+disp++] = vbd[IJK]->zone[ijk];
                        if (update_bits&UPDATE_FLAG) 
                          send_buf_int[send_disp[rank]*int_width+disp++] = vbd[IJK]->flag[ijk];
                        disp = 0;
                        if (update_bits&UPDATE_DATA) 
                          send_buf_double[send_disp[rank]*double_width+disp++] = vbd[IJK]->data[ijk];
                        if (update_bits&UPDATE_DATA0) 
                          send_buf_double[send_disp[rank]*double_width+disp++] = vbd[IJK]->data0[ijk];
                        if (update_bits&UPDATE_DATA00) 
                          send_buf_double[send_disp[rank]*double_width+disp++] = vbd[IJK]->data00[ijk];
                        ++send_disp[rank];
                      }
                    }
                  }
                }
              }
            }
          }
        }

      }

      send_disp[0] = 0;
      for (int rank = 1; rank < mpi_size; ++rank)
        send_disp[rank] = send_count[rank-1] + send_disp[rank-1];

      if (iter == 0) {
        const int send_count_sum = send_disp[mpi_size-1] + send_count[mpi_size-1];
        send_buf_int = new int[int_width*send_count_sum];
        send_buf_double = new double[double_width*send_count_sum];
      }

    }

    // exchange...

    int *recv_count = new int[mpi_size];
    MPI_Alltoall(send_count,1,MPI_INT,recv_count,1,MPI_INT,mpi_comm);

    int *recv_disp = new int[mpi_size];
    recv_disp[0] = 0;
    for (int rank = 1; rank < mpi_size; ++rank)
      recv_disp[rank] = recv_count[rank-1] + recv_disp[rank-1];
    const int recv_count_sum = recv_disp[mpi_size-1] + recv_count[mpi_size-1];

    FOR_RANK {
      send_count[rank] *= double_width;
      send_disp[rank] *= double_width;
      recv_count[rank] *= double_width;
      recv_disp[rank] *= double_width;
    }

    double * recv_buf_double = new double[double_width*recv_count_sum];
    MPI_Alltoallv(send_buf_double,send_count,send_disp,MPI_DOUBLE,
                  recv_buf_double,recv_count,recv_disp,MPI_DOUBLE,mpi_comm);
    delete[] send_buf_double; 

    FOR_RANK {
      send_count[rank] = (send_count[rank]/double_width)*int_width;;
      send_disp[rank] = (send_disp[rank]/double_width)*int_width;;
      recv_count[rank] = (recv_count[rank]/double_width)*int_width;;
      recv_disp[rank] = (recv_disp[rank]/double_width)*int_width;;
    }

    int * recv_buf_int = new int[int_width*recv_count_sum];
    MPI_Alltoallv(send_buf_int,send_count,send_disp,MPI_INT,
                  recv_buf_int,recv_count,recv_disp,MPI_INT,mpi_comm);
    delete[] send_buf_int; 
    delete[] send_count;
    delete[] send_disp;
    delete[] recv_count;
    delete[] recv_disp;
    
    // now set our data...

    for (int irecv = 0; irecv < recv_count_sum; ++irecv) {
      int disp = 0;
      const int IJK = recv_buf_int[irecv*int_width+disp++]; assert(vbd_rank[IJK] == mpi_rank);
      const int ijk = recv_buf_int[irecv*int_width+disp++];
      assert(vbd[IJK] != NULL);
      if (update_bits&UPDATE_ZONE) 
        vbd[IJK]->zone[ijk] = recv_buf_int[irecv*int_width+disp++];
      if (update_bits&UPDATE_FLAG) 
        vbd[IJK]->flag[ijk] = recv_buf_int[irecv*int_width+disp++];
      disp = 0;
      if (update_bits&UPDATE_DATA) 
        vbd[IJK]->data[ijk] = recv_buf_double[irecv*double_width+disp++];
      if (update_bits&UPDATE_DATA0) 
        vbd[IJK]->data0[ijk] = recv_buf_double[irecv*double_width+disp++];
      if (update_bits&UPDATE_DATA00) 
        vbd[IJK]->data00[ijk] = recv_buf_double[irecv*double_width+disp++];
    }
    delete[] recv_buf_int;
    delete[] recv_buf_double;

  }

  void setVoxelData(const int i,const int j,const int k,const double data) {

    assert((i >= 0)&&(i < ni));
    assert((j >= 0)&&(j < nj));
    assert((k >= 0)&&(k < nk));
    const int IJK = ((k>>BLOCK_BITS)*nJ+(j>>BLOCK_BITS))*nI + (i>>BLOCK_BITS);
    assert((IJK >= 0)&&(IJK < nI*nJ*nK));
    if (vbd[IJK] == NULL) vbd[IJK] = new VoxelBlockData();
    const int ijk = (((k&BLOCK_MASK)+N_GHOSTS)*BLOCK_SIZE_G+((j&BLOCK_MASK)+N_GHOSTS))*BLOCK_SIZE_G+((i&BLOCK_MASK)+N_GHOSTS);
    assert((ijk >= N_GHOSTS)&&(ijk < BLOCK_SIZE3_G));

    vbd[IJK]->data[ijk] = data;

  }

  bool getVoxelData(double& data,const int i,const int j,const int k) const {

    assert((i >= 0)&&(i < ni));
    assert((j >= 0)&&(j < nj));
    assert((k >= 0)&&(k < nk));
    const int IJK = ((k>>BLOCK_BITS)*nJ+(j>>BLOCK_BITS))*nI + (i>>BLOCK_BITS);
    if (vbd[IJK] == NULL) return false;
    assert((IJK >= 0)&&(IJK < nI*nJ*nK));
    const int ijk = (((k&BLOCK_MASK)+N_GHOSTS)*BLOCK_SIZE_G+((j&BLOCK_MASK)+N_GHOSTS))*BLOCK_SIZE_G+((i&BLOCK_MASK)+N_GHOSTS);
    assert((ijk >= N_GHOSTS)&&(ijk < BLOCK_SIZE3_G));

    data = vbd[IJK]->data[ijk];
    return true;
  }

  void writeVoxelDataToTecplot(const string& filename) const {

    if (mpi_rank == 0) cout << "writeVoxelDataToTecplot: " << filename << endl;

    // do round-robin writing of points to tecplot...

    FILE * fp;
    if ( mpi_rank == 0 ) {
      fp = fopen(filename.c_str(),"w");
      assert(fp != NULL);
      fprintf(fp,"TITLE = \"Voxel Data\"\n");
      fprintf(fp,"VARIABLES = \"X\"\n");
      fprintf(fp,"\"Y\"\n");
      fprintf(fp,"\"Z\"\n");
      fprintf(fp,"\"DATA\"\n");
      fprintf(fp,"\"ZONE\"\n");
      fprintf(fp,"\"FLAG\"\n");
    }
    else {
      int dummy;
      MPI_Status status;
      MPI_Recv(&dummy,1,MPI_INT,mpi_rank-1,1234,mpi_comm,&status);
      fp = fopen(filename.c_str(),"a");
      assert(fp != NULL);
    }

    for (int ivbd = 0, limit = vbdVec.size(); ivbd < limit; ++ivbd) {
      const int IJK = vbdVec[ivbd];
      assert(vbd[IJK] != NULL); 
      const int I = IJK%nI;
      const int J = (IJK/nI)%nJ;
      const int K = IJK/(nI*nJ);
      assert(IJK == (K*nJ+J)*nI+I);
      double x[3];
      for (int k = N_GHOSTS; k < BLOCK_SIZE_G-N_GHOSTS; ++k) {
        x[2] = ((K<<BLOCK_BITS)+k-N_GHOSTS)*dx+x0[2];
        for (int j = N_GHOSTS; j < BLOCK_SIZE_G-N_GHOSTS; ++j) {
          x[1] = ((J<<BLOCK_BITS)+j-N_GHOSTS)*dx+x0[1];
          for (int i = N_GHOSTS; i < BLOCK_SIZE_G-N_GHOSTS; ++i) {
            x[0] = ((I<<BLOCK_BITS)+i-N_GHOSTS)*dx+x0[0];
            const int ijk = (k*BLOCK_SIZE_G+j)*BLOCK_SIZE_G+i;
            assert((ijk >= N_GHOSTS)&&(ijk < BLOCK_SIZE3_G));
            //if (vbd[IJK]->flag[ijk]&SURFACE_BIT) 
            if (vbd[IJK]->data[ijk] != HUGE_VAL) 
              fprintf(fp,"%18.15le %18.15le %18.15le %18.15le %18.15le %18.15le\n",
                  x[0],x[1],x[2],vbd[IJK]->data[ijk],double(vbd[IJK]->zone[ijk]),double(vbd[IJK]->flag[ijk]));
          }
        }
      }
    }

    fclose(fp);

    if ( mpi_rank < mpi_size-1 ) {
      int dummy = 1;
      MPI_Send(&dummy,1,MPI_INT,mpi_rank+1,1234,mpi_comm);
    }

    MPI_Barrier(mpi_comm);

  }

  void vertexInterp(double p[3],const double iso_level,const double p1[3],const double p2[3],
      const double valp1,const double valp2) const {

    if (fabs(iso_level-valp1) < tol*dx) {
      FOR_I3 p[i] = p1[i];
    }
    else if (fabs(iso_level-valp2) < tol*dx) {
      FOR_I3 p[i] = p2[i];
    }
    else if (fabs(valp1-valp2) < tol*dx) {
      FOR_I3 p[i] = p1[i];
    }
    else {
      const double mu = (iso_level - valp1) / (valp2 - valp1);
      FOR_I3 p[i] = p1[i] + mu*(p2[i]-p1[i]);
    }

  }
  
  // here we are staggering the voxel, where il,jl,kl are the min coordinates...
  // we do this to avoid interpolation...

  void triangulateStaggeredVoxel(vector<SimplePoint>& vertexVec,vector<int>& triZoneVec,
      const int IJK,const int il,const int jl,const int kl,const double iso_level) const {
    assert(N_GHOSTS >= 1);
    // make sure this is an internal location...
    assert((IJK >= 0)&&(IJK < nI*nJ*nK)&&(vbd_rank[IJK] == mpi_rank));
    assert((il >= N_GHOSTS)&&(il < BLOCK_SIZE_G-N_GHOSTS));
    assert((jl >= N_GHOSTS)&&(jl < BLOCK_SIZE_G-N_GHOSTS));
    assert((kl >= N_GHOSTS)&&(kl < BLOCK_SIZE_G-N_GHOSTS));

    // only build if the your next pt in a direction is internal or, if it is a ghost, your peer rank exists

    const int I = IJK%nI;
    const int J = (IJK/nI)%nJ;
    const int K = IJK/(nI*nJ);
    assert(IJK == (K*nJ+J)*nI+I);
    if (!((((il+1) >= N_GHOSTS)&&((il+1) < BLOCK_SIZE_G-N_GHOSTS))||  
          (((il+1) < N_GHOSTS)&&
           (((((K+0)*nJ+J+0)*nI+I-1) >= 0)&&((((K+0)*nJ+J+0)*nI+I-1) < nI*nJ*nK)&&(vbd_rank[((K+0)*nJ+J+0)*nI+I-1] >= 0)))||
          (((il+1) >= BLOCK_SIZE_G-N_GHOSTS)&&
           (((((K+0)*nJ+J+0)*nI+I+1) >= 0)&&((((K+0)*nJ+J+0)*nI+I+1) < nI*nJ*nK)&&(vbd_rank[((K+0)*nJ+J+0)*nI+I+1] >= 0))))) {
      return;
    }
    if (!((((jl+1) >= N_GHOSTS)&&((jl+1) < BLOCK_SIZE_G-N_GHOSTS))||  
          (((jl+1) < N_GHOSTS)&&
           (((((K+0)*nJ+J-1)*nI+I+0) >= 0)&&((((K+0)*nJ+J-1)*nI+I+0) < nI*nJ*nK)&&(vbd_rank[((K+0)*nJ+J-1)*nI+I+0] >= 0)))||
          (((jl+1) >= BLOCK_SIZE_G-N_GHOSTS)&&
           (((((K+0)*nJ+J+1)*nI+I+0) >= 0)&&((((K+0)*nJ+J+1)*nI+I+0) < nI*nJ*nK)&&(vbd_rank[((K+0)*nJ+J+1)*nI+I+0] >= 0))))) {
      return;
    }
    if (!((((kl+1) >= N_GHOSTS)&&((kl+1) < BLOCK_SIZE_G-N_GHOSTS))||  
          (((kl+1) < N_GHOSTS)&&
           (((((K-1)*nJ+J+0)*nI+I+0) >= 0)&&((((K-1)*nJ+J+0)*nI+I+0) < nI*nJ*nK)&&(vbd_rank[((K-1)*nJ+J+0)*nI+I+0] >= 0)))||
          (((kl+1) >= BLOCK_SIZE_G-N_GHOSTS)&&
           (((((K+1)*nJ+J+0)*nI+I+0) >= 0)&&((((K+1)*nJ+J+0)*nI+I+0) < nI*nJ*nK)&&(vbd_rank[((K+1)*nJ+J+0)*nI+I+0] >= 0))))) {
      return;
    }

    // get grid values...

    double grid_data[8];
    grid_data[0] = vbd[IJK]->data[((kl+0)*BLOCK_SIZE_G+jl+0)*BLOCK_SIZE_G+il+0]; if (grid_data[0] == HUGE_VAL) return;
    grid_data[1] = vbd[IJK]->data[((kl+0)*BLOCK_SIZE_G+jl+0)*BLOCK_SIZE_G+il+1]; if (grid_data[1] == HUGE_VAL) return;
    grid_data[2] = vbd[IJK]->data[((kl+0)*BLOCK_SIZE_G+jl+1)*BLOCK_SIZE_G+il+1]; if (grid_data[2] == HUGE_VAL) return;
    grid_data[3] = vbd[IJK]->data[((kl+0)*BLOCK_SIZE_G+jl+1)*BLOCK_SIZE_G+il+0]; if (grid_data[3] == HUGE_VAL) return;
    grid_data[4] = vbd[IJK]->data[((kl+1)*BLOCK_SIZE_G+jl+0)*BLOCK_SIZE_G+il+0]; if (grid_data[4] == HUGE_VAL) return;
    grid_data[5] = vbd[IJK]->data[((kl+1)*BLOCK_SIZE_G+jl+0)*BLOCK_SIZE_G+il+1]; if (grid_data[5] == HUGE_VAL) return;
    grid_data[6] = vbd[IJK]->data[((kl+1)*BLOCK_SIZE_G+jl+1)*BLOCK_SIZE_G+il+1]; if (grid_data[6] == HUGE_VAL) return;
    grid_data[7] = vbd[IJK]->data[((kl+1)*BLOCK_SIZE_G+jl+1)*BLOCK_SIZE_G+il+0]; if (grid_data[7] == HUGE_VAL) return;

    // get global voxel location...

    const int i = (I<<BLOCK_BITS)+il-N_GHOSTS;
    const int j = (J<<BLOCK_BITS)+jl-N_GHOSTS;
    const int k = (K<<BLOCK_BITS)+kl-N_GHOSTS;

    // grid corner positions...
    
    double grid_x[8][3];
    grid_x[0][0] = (double(i)+0.0)*dx+x0[0];
    grid_x[0][1] = (double(j)+0.0)*dx+x0[1];
    grid_x[0][2] = (double(k)+0.0)*dx+x0[2];
    grid_x[1][0] = (double(i)+1.0)*dx+x0[0];
    grid_x[1][1] = (double(j)+0.0)*dx+x0[1];
    grid_x[1][2] = (double(k)+0.0)*dx+x0[2];
    grid_x[2][0] = (double(i)+1.0)*dx+x0[0];
    grid_x[2][1] = (double(j)+1.0)*dx+x0[1];
    grid_x[2][2] = (double(k)+0.0)*dx+x0[2];
    grid_x[3][0] = (double(i)+0.0)*dx+x0[0];
    grid_x[3][1] = (double(j)+1.0)*dx+x0[1];
    grid_x[3][2] = (double(k)+0.0)*dx+x0[2];
    grid_x[4][0] = (double(i)+0.0)*dx+x0[0];
    grid_x[4][1] = (double(j)+0.0)*dx+x0[1];
    grid_x[4][2] = (double(k)+1.0)*dx+x0[2];
    grid_x[5][0] = (double(i)+1.0)*dx+x0[0];
    grid_x[5][1] = (double(j)+0.0)*dx+x0[1];
    grid_x[5][2] = (double(k)+1.0)*dx+x0[2];
    grid_x[6][0] = (double(i)+1.0)*dx+x0[0];
    grid_x[6][1] = (double(j)+1.0)*dx+x0[1];
    grid_x[6][2] = (double(k)+1.0)*dx+x0[2];
    grid_x[7][0] = (double(i)+0.0)*dx+x0[0];
    grid_x[7][1] = (double(j)+1.0)*dx+x0[1];
    grid_x[7][2] = (double(k)+1.0)*dx+x0[2];

    // determine the index into the edge table which tells us 
    // which vertices are inside of the surface

    int cube_index = 0;
    if (grid_data[0] < iso_level) cube_index |= 1;
    if (grid_data[1] < iso_level) cube_index |= 2;
    if (grid_data[2] < iso_level) cube_index |= 4;
    if (grid_data[3] < iso_level) cube_index |= 8;
    if (grid_data[4] < iso_level) cube_index |= 16;
    if (grid_data[5] < iso_level) cube_index |= 32;
    if (grid_data[6] < iso_level) cube_index |= 64;
    if (grid_data[7] < iso_level) cube_index |= 128;
 
    // cube is entirely in/out of the surface 
   
    if (edgeTable[cube_index] == 0) 
      return;

    // find the vertices where the surface intersects the cube 

    double vertlist[12][3];
    if (edgeTable[cube_index] & 1)
      vertexInterp(vertlist[ 0],iso_level,grid_x[0],grid_x[1],grid_data[0],grid_data[1]);
    if (edgeTable[cube_index] & 2)
      vertexInterp(vertlist[ 1],iso_level,grid_x[1],grid_x[2],grid_data[1],grid_data[2]);
    if (edgeTable[cube_index] & 4)
      vertexInterp(vertlist[ 2],iso_level,grid_x[2],grid_x[3],grid_data[2],grid_data[3]);
    if (edgeTable[cube_index] & 8)
      vertexInterp(vertlist[ 3],iso_level,grid_x[3],grid_x[0],grid_data[3],grid_data[0]);
    if (edgeTable[cube_index] & 16)
      vertexInterp(vertlist[ 4],iso_level,grid_x[4],grid_x[5],grid_data[4],grid_data[5]);
    if (edgeTable[cube_index] & 32)
      vertexInterp(vertlist[ 5],iso_level,grid_x[5],grid_x[6],grid_data[5],grid_data[6]);
    if (edgeTable[cube_index] & 64)
      vertexInterp(vertlist[ 6],iso_level,grid_x[6],grid_x[7],grid_data[6],grid_data[7]);
    if (edgeTable[cube_index] & 128)
      vertexInterp(vertlist[ 7],iso_level,grid_x[7],grid_x[4],grid_data[7],grid_data[4]);
    if (edgeTable[cube_index] & 256)
      vertexInterp(vertlist[ 8],iso_level,grid_x[0],grid_x[4],grid_data[0],grid_data[4]);
    if (edgeTable[cube_index] & 512)
      vertexInterp(vertlist[ 9],iso_level,grid_x[1],grid_x[5],grid_data[1],grid_data[5]);
    if (edgeTable[cube_index] & 1024)
      vertexInterp(vertlist[10],iso_level,grid_x[2],grid_x[6],grid_data[2],grid_data[6]);
    if (edgeTable[cube_index] & 2048)
      vertexInterp(vertlist[11],iso_level,grid_x[3],grid_x[7],grid_data[3],grid_data[7]);
   
    // get zones at corners...

    int grid_zone[8];
    grid_zone[0] = vbd[IJK]->zone[((kl+0)*BLOCK_SIZE_G+jl+0)*BLOCK_SIZE_G+il+0]; 
    grid_zone[1] = vbd[IJK]->zone[((kl+0)*BLOCK_SIZE_G+jl+0)*BLOCK_SIZE_G+il+1]; 
    grid_zone[2] = vbd[IJK]->zone[((kl+0)*BLOCK_SIZE_G+jl+1)*BLOCK_SIZE_G+il+1]; 
    grid_zone[3] = vbd[IJK]->zone[((kl+0)*BLOCK_SIZE_G+jl+1)*BLOCK_SIZE_G+il+0]; 
    grid_zone[4] = vbd[IJK]->zone[((kl+1)*BLOCK_SIZE_G+jl+0)*BLOCK_SIZE_G+il+0]; 
    grid_zone[5] = vbd[IJK]->zone[((kl+1)*BLOCK_SIZE_G+jl+0)*BLOCK_SIZE_G+il+1]; 
    grid_zone[6] = vbd[IJK]->zone[((kl+1)*BLOCK_SIZE_G+jl+1)*BLOCK_SIZE_G+il+1]; 
    grid_zone[7] = vbd[IJK]->zone[((kl+1)*BLOCK_SIZE_G+jl+1)*BLOCK_SIZE_G+il+0]; 

    // add triangles...

    for (int l = 0; triTable[cube_index][l] != -1; l += 3) {
      vertexVec.push_back(SimplePoint(vertlist[triTable[cube_index][l  ]]));
      vertexVec.push_back(SimplePoint(vertlist[triTable[cube_index][l+1]]));
      vertexVec.push_back(SimplePoint(vertlist[triTable[cube_index][l+2]]));

      // take zone from closest corner to centroid...
      double xc[3];
      FOR_M3 xc[m]  = vertlist[triTable[cube_index][l  ]][m]; 
      FOR_M3 xc[m] += vertlist[triTable[cube_index][l+1]][m]; 
      FOR_M3 xc[m] += vertlist[triTable[cube_index][l+2]][m]; 
      FOR_M3 xc[m] /= 3.0;
      double dist2 = HUGE_VAL;
      int zone = -1;
      for (int ii = 0; ii < 8; ++ii) {
        const double this_dist2 = DIST2(xc,grid_x[ii]);
        if (this_dist2 <= dist2) {
          zone = grid_zone[ii];
          dist2 = this_dist2;
        }
      }
      assert(zone != -1);
      triZoneVec.push_back(zone);
    }

  }

  void triangulateVoxelData(SurfaceShm *&surface,const double iso_level = 0.0,const bool b_write = false) const {

    if (mpi_rank == 0) cout << "triangulateVoxelData()" << endl;

    // build triangles on your rank...
    
    vector<SimplePoint> vertexVec;
    vector<int> triZoneVec;
    for (int ivbd = 0, limit = vbdVec.size(); ivbd < limit; ++ivbd) {
      const int IJK = vbdVec[ivbd];
      assert(vbd[IJK] != NULL); 
      for (int k = N_GHOSTS; k < BLOCK_SIZE_G-N_GHOSTS; ++k) {
        for (int j = N_GHOSTS; j < BLOCK_SIZE_G-N_GHOSTS; ++j) {
          for (int i = N_GHOSTS; i < BLOCK_SIZE_G-N_GHOSTS; ++i) {
            triangulateStaggeredVoxel(vertexVec,triZoneVec,IJK,i,j,k,iso_level);
          }
        }
      }
    }
    assert(triZoneVec.size()*3 == vertexVec.size());

    // in parallel, just gather everyone at root...
  
    int *recv_count = NULL;
    if (mpi_rank == 0) recv_count = new int[mpi_size];
    int send_count = triZoneVec.size();
    MPI_Gather(&send_count,1,MPI_INT,recv_count,1,MPI_INT,0,mpi_comm);

    int *recv_disp = NULL;
    int recv_count_sum = -1;
    int * recv_buf_int = NULL;
    if (mpi_rank == 0) {
      recv_disp = new int[mpi_size];
      recv_disp[0] = 0;
      for (int rank = 1; rank < mpi_size; ++rank)
        recv_disp[rank] = recv_disp[rank-1] + recv_count[rank-1];
      recv_count_sum = recv_disp[mpi_size-1] + recv_count[mpi_size-1];
      recv_buf_int = new int[recv_count_sum];
    }
    
    MPI_Gatherv(&triZoneVec[0],send_count,MPI_INT,recv_buf_int,recv_count,recv_disp,MPI_INT,0,mpi_comm);
    triZoneVec.clear();

    double *recv_buf_double = NULL;
    send_count *= 9;
    if (mpi_rank == 0) {
      recv_buf_double = new double[9*recv_count_sum]; 
      FOR_RANK {
        recv_count[rank] *= 9;
        recv_disp[rank] *= 9;
      }
    }

    MPI_Gatherv((double*)(&vertexVec[0].x[0]),send_count,MPI_DOUBLE,recv_buf_double,recv_count,recv_disp,MPI_DOUBLE,0,mpi_comm);
    vertexVec.clear();

    if (mpi_rank == 0) {
      delete[] recv_count;
      delete[] recv_disp;
      for (int irecv = 0; irecv < recv_count_sum; ++irecv) {
        triZoneVec.push_back(recv_buf_int[irecv]);
        vertexVec.push_back(SimplePoint(&recv_buf_double[irecv*9  ]));
        vertexVec.push_back(SimplePoint(&recv_buf_double[irecv*9+3]));
        vertexVec.push_back(SimplePoint(&recv_buf_double[irecv*9+6]));
      }
      delete[] recv_buf_int;
      delete[] recv_buf_double;
    }

    if (mpi_rank == 0) 
      cout << " > nsp and nst after triangulation: " << vertexVec.size() << " " << vertexVec.size()/3 << endl;

    // check....
    //if (mpi_rank == 0) writeVertexVecTriplesToTecplot(vertexVec,"triangulate.dat");
    
    // build simple shared memory surface...
    
    assert(surface == NULL);
    surface = new SurfaceShm();

    int *flag = NULL;
    int nt = -1,nv = -1;
    if (mpi_rank == 0) {

      // now index and build bbox
      const int nsp = vertexVec.size();
      double (*bbmin)[3] = new double[nsp][3];
      double (*bbmax)[3] = new double[nsp][3];
      int *sp_flag = new int[nsp];
      const double merge_tol = tol*dx;
      FOR_ISP {
        FOR_I3 {
          bbmin[isp][i] = vertexVec[isp].x[i] - merge_tol;
          bbmax[isp][i] = vertexVec[isp].x[i] + merge_tol;
        }
        sp_flag[isp] = isp;  // reset sp_flag to own index - used next for compress
      }
      Adt<double> *nodesAdt = new Adt<double>(nsp,bbmin,bbmax);
      delete[] bbmin;
      delete[] bbmax;

      // match candidates
      vector<int> candidateVec;
      FOR_ISP {
        candidateVec.clear();
        nodesAdt->buildListForPoint(candidateVec,vertexVec[isp].x);
        for (vector<int>::iterator it=candidateVec.begin(); it!=candidateVec.end(); ++it) {
          const int c_isp = *it;
          if (c_isp < isp) {
            sp_flag[isp] = c_isp;
          }
          else if (c_isp > isp) {
            sp_flag[c_isp] = isp;
          }
          else {
            assert(c_isp == isp);
            continue; // don't process myself
          }
        }
      }
      delete nodesAdt;

      // sp_flag should now be indexed such that nodes who are being compressed index their
      // lower partner

      // push sp_flag to unique value
      FOR_ISP {
        if (sp_flag[isp] != isp) {
          int isp_t = sp_flag[isp];
          while (sp_flag[isp_t] != isp_t) isp_t = sp_flag[isp_t];
          sp_flag[isp] = isp_t;
        }
      }

      // overwrite data with "equivalent" lowest index...
      FOR_ISP {
        if (sp_flag[isp] != isp) 
          FOR_I3 vertexVec[isp].x[i] = vertexVec[sp_flag[isp]].x[i];
      }
      delete[] sp_flag;

      nv = 0;
      map<const pair<SimplePoint,SimplePoint>,int> edgeMap;
      const int ntri = vertexVec.size()/3;
      vector<int> edgeVec(vertexVec.size());
      for (int itri = 0; itri < ntri; ++itri) {
        const int iv0 = nv++;
        //const int iv1 = nv++;
        //const int iv2 = nv++;
        nv += 2;

        // look for edge pairs...
        for (int i = 0; i < 3; ++i) {
          map<const pair<SimplePoint,SimplePoint>,int>::iterator iter = 
            edgeMap.find(pair<SimplePoint,SimplePoint>(vertexVec[iv0+i],vertexVec[iv0+(i+1)%3]));

          if (iter != edgeMap.end()) {
            assert(edgeVec[iter->second] == -1);
            edgeVec[iter->second] = iv0+i;
            edgeVec[iv0+i] = iter->second;
            // found a match...
            edgeMap.erase(iter);
          }
          else {
            // not found, so add the reverse edge to edgeMap...
            edgeMap[pair<SimplePoint,SimplePoint>(vertexVec[iv0+(i+1)%3],vertexVec[iv0+i])] = iv0+i;
            edgeVec[iv0+i] = -1;
          }
        }
      }
      edgeMap.clear();

      // check...

      assert(vertexVec.size()%3 == 0);
      nt = vertexVec.size()/3;
      assert(edgeVec.size() == vertexVec.size());

      // compress nodes...

      flag = new int[vertexVec.size()];
      for (int iv = 0, nv_all = vertexVec.size(); iv < nv_all; ++iv) flag[iv] = iv;

      for (int ie = 0, ne = edgeVec.size(); ie < ne; ++ie) {
        if (edgeVec[ie] != -1) {
          const int ie_match = edgeVec[ie];
          assert(ie_match != ie);
          assert((ie_match >= 0)&&(ie_match < ne));
          assert(edgeVec[ie_match] == ie);
          if (ie_match > ie) {
            {
              int iv0 = ie;
              int iv1_match = ie_match+1; if (iv1_match%3 == 0) iv1_match -= 3;
              assert(vertexVec[iv0] == vertexVec[iv1_match]);
              int iv0_ = flag[iv0];
              while (iv0_ != flag[iv0_]) iv0_ = flag[iv0_];
              int iv1_match_ = flag[iv1_match];
              while (iv1_match_ != flag[iv1_match_]) iv1_match_ = flag[iv1_match_];
              flag[iv0_] = flag[iv1_match_] = min(iv0_,iv1_match_);
            }
            {
              int iv1 = ie+1; if (iv1%3 == 0) iv1 -= 3;
              int iv0_match = ie_match;
              assert(vertexVec[iv1] == vertexVec[iv0_match]);
              int iv1_ = flag[iv1];
              while (iv1_ != flag[iv1_]) iv1_ = flag[iv1_];
              int iv0_match_ = flag[iv0_match];
              while (iv0_match_ != flag[iv0_match_]) iv0_match_ = flag[iv0_match_];
              flag[iv1_] = flag[iv0_match_] = min(iv1_,iv0_match_);
            }
          }
        }
      }
      edgeVec.clear();

      nv = 0;
      for (int iv = 0, nv_all = vertexVec.size(); iv < nv_all; ++iv) {
        if (flag[iv] == iv) {
          ++nv;
          flag[iv] = -nv; // -1 indexing
        }
        else {
          int iv_ = flag[iv];
          while (iv_ >= 0) iv_ = flag[iv_];
          flag[iv] = iv_;
        }
      }

      cout << " > nsp and nst after reduction: " << nv << " " << nt << endl;
    }

    // ========================================
    // set the simple surface stuff...
    // ========================================

    // broadcast sizes across nodes...
    int size_buf[2];
    if (mpi_rank == 0) {
      size_buf[0] = nv;
      size_buf[1] = nt;
    }
    MPI_Bcast(size_buf,2,MPI_INT,0,mpi_comm); 

    surface->init(size_buf[0],size_buf[1]);

    if (mpi_rank_shared == 0) {
      assert(surface->zoneVec.empty());
      surface->zoneVec.resize(zoneVec.size());
      for (int iz = 0; iz < zoneVec.size(); ++iz)
        surface->zoneVec[iz] = zoneVec[iz];

      if (mpi_rank == 0) {
        nv = 0;
        for (int iv = 0, nv_all = vertexVec.size(); iv < nv_all; ++iv) {
          if (flag[iv] == -nv-1) {
            FOR_I3 surface->xsp[nv][i] = vertexVec[iv].x[i];
            ++nv;
          }
        }
        assert(nv == surface->nsp);
        vertexVec.clear();

        for (int it = 0; it < nt; ++it) { 
          FOR_I3 surface->spost[it][i] = -flag[it*3+i]-1;
          surface->znost[it] = triZoneVec[it];
        }
        triZoneVec.clear();
        delete[] flag;

      }
      MPI_Bcast((double*)surface->xsp,3*surface->nsp,MPI_DOUBLE,0,mpi_comm_internode); 
      MPI_Bcast((int*)surface->spost,3*surface->nst,MPI_INT,0,mpi_comm_internode); 
      MPI_Bcast(surface->znost,surface->nst,MPI_INT,0,mpi_comm_internode); 

    }
    MPI_Barrier(mpi_comm_shared);

    // check....
    if (b_write)
      surface->writeTecplot("triangulate.dat");

  }

  void writeVertexVecTriplesToTecplot(const vector<SimplePoint>& vertexVec,const string& filename) const {
    assert(mpi_rank == 0);

    FILE * fp = fopen(filename.c_str(),"w");
    assert(fp != NULL);

    fprintf(fp,"TITLE = \"%s\"\n",filename.c_str());
    fprintf(fp,"VARIABLES = \"X\"\n");
    fprintf(fp,"\"Y\"\n");
    fprintf(fp,"\"Z\"\n");

    assert(vertexVec.size()%3 == 0);
    const int ntri = vertexVec.size()/3;
    fprintf(fp,"ZONE T=\"%s\"\n",filename.c_str());
    fprintf(fp,"N=%d, E=%d, F=FEPOINT, ET=TRIANGLE\n",3*ntri,ntri);

    for (int itri = 0; itri < ntri; ++itri) {
      fprintf(fp,"%lf %lf %lf\n",vertexVec[3*itri  ].x[0],vertexVec[3*itri  ].x[1],vertexVec[3*itri  ].x[2]);
      fprintf(fp,"%lf %lf %lf\n",vertexVec[3*itri+1].x[0],vertexVec[3*itri+1].x[1],vertexVec[3*itri+1].x[2]);
      fprintf(fp,"%lf %lf %lf\n",vertexVec[3*itri+2].x[0],vertexVec[3*itri+2].x[1],vertexVec[3*itri+2].x[2]);
    }
    for (int itri = 0; itri < ntri; ++itri) {
      fprintf(fp,"%d %d %d\n",3*itri+1,3*itri+2,3*itri+3); // 1-indexed
    }

    fclose(fp);
  }

#define FINDMINMAX(x0,x1,x2,min,max) \
  min = max = x0;                    \
  if (x1 < min) min = x1;            \
  if (x1 > max) max = x1;            \
  if (x2 < min) min = x2;            \
  if (x2 > max) max = x2;

#define AXISTEST_X01(a,b,fa,fb)                                 \
  p0 = a*v0[1] - b*v0[2];                                       \
  p2 = a*v2[1] - b*v2[2];                                       \
  if (p0 < p2) {min = p0; max = p2;} else {min = p2; max = p0;} \
  rad = fa * boxhalfsize[1] + fb * boxhalfsize[2];              \
  if (min > rad || max < -rad) return 0;

#define AXISTEST_X2(a,b,fa,fb)                                  \
  p0 = a*v0[1] - b*v0[2];                                       \
  p1 = a*v1[1] - b*v1[2];                                       \
  if (p0 < p1) {min = p0; max = p1;} else {min = p1; max = p0;} \
  rad = fa * boxhalfsize[1] + fb * boxhalfsize[2];              \
  if (min > rad || max < -rad) return 0;

#define AXISTEST_Y02(a,b,fa,fb)                                 \
  p0 = -a*v0[0] + b*v0[2];                                      \
  p2 = -a*v2[0] + b*v2[2];                                      \
  if (p0 < p2) {min = p0; max = p2;} else {min = p2; max = p0;} \
  rad = fa * boxhalfsize[0] + fb * boxhalfsize[2];              \
  if (min > rad || max < -rad) return 0;

#define AXISTEST_Y1(a,b,fa,fb)                                  \
  p0 = -a*v0[0] + b*v0[2];                                      \
  p1 = -a*v1[0] + b*v1[2];                                      \
  if (p0 < p1) {min = p0; max = p1;} else {min = p1; max = p0;} \
  rad = fa * boxhalfsize[0] + fb * boxhalfsize[2];              \
  if (min > rad || max < -rad) return 0;

#define AXISTEST_Z12(a,b,fa,fb)                                 \
  p1 = a*v1[0] - b*v1[1];                                       \
  p2 = a*v2[0] - b*v2[1];                                       \
  if (p2 < p1) {min = p2; max = p1;} else {min = p1; max = p2;} \
  rad = fa * boxhalfsize[0] + fb * boxhalfsize[1];              \
  if (min > rad || max < -rad) return 0;

#define AXISTEST_Z0(a,b,fa,fb)                                  \
  p0 = a*v0[0] - b*v0[1];                                       \
  p1 = a*v1[0] - b*v1[1];                                       \
  if (p0 < p1) {min = p0; max = p1;} else {min = p1; max = p0;} \
  rad = fa * boxhalfsize[0] + fb * boxhalfsize[1];              \
  if (min > rad || max < -rad) return 0;

  bool planeBoxOverlap(const double normal[3],const double vert[3],const double maxbox[3]) {

    double vmin[3],vmax[3];
    for (int q = 0; q < 3; q++) {
      const double v = vert[q];      
      if (normal[q] > 0.0) {
        vmin[q] = -maxbox[q] - v; 
        vmax[q] =  maxbox[q] - v;  
      }
      else {
        vmin[q] =  maxbox[q] - v; 
        vmax[q] = -maxbox[q] - v; 
      }
    }

    if (DOT_PRODUCT(normal,vmin) >  0.0) return false; 
    if (DOT_PRODUCT(normal,vmax) >= 0.0) return true; 

    return false;

  }

  bool triBoxOverlap(const double boxcenter[3],const double boxhalfsize[3],const double triverts[3][3]) {
   
    // algo by Tomas Akenine-Moller 
    // use separating axis theorem to test overlap between triangle and box 
    //   need to test for overlap in these directions:
    //   1) the {x,y,z}-directions (actually,since we use the AABB of the triangle 
    //      we do not even need to test these) 
    //   2) normal of the triangle 
    //   3) crossproduct(edge from tri,{x,y,z}-direction) 
    //      this gives 3x3=9 more tests 

    double min,max,p0,p1,p2,rad,fex,fey,fez;

    // move everything so that the boxcenter is in (0,0,0) 

    const double v0[3] = DIFF(triverts[0],boxcenter);
    const double v1[3] = DIFF(triverts[1],boxcenter);
    const double v2[3] = DIFF(triverts[2],boxcenter);

    // compute triangle edges 
    
    const double e0[3] = DIFF(v1,v0);
    const double e1[3] = DIFF(v2,v1);
    const double e2[3] = DIFF(v0,v2); 

    // Bullet 3:
    // test the 9 tests first (this was faster) 
    
    fex = fabs(e0[0]);
    fey = fabs(e0[1]);
    fez = fabs(e0[2]);

    AXISTEST_X01(e0[2],e0[1],fez,fey);
    AXISTEST_Y02(e0[2],e0[0],fez,fex);
    AXISTEST_Z12(e0[1],e0[0],fey,fex);

    fex = fabs(e1[0]);
    fey = fabs(e1[1]);
    fez = fabs(e1[2]);

    AXISTEST_X01(e1[2],e1[1],fez,fey);
    AXISTEST_Y02(e1[2],e1[0],fez,fex);
    AXISTEST_Z0(e1[1],e1[0],fey,fex);

    fex = fabs(e2[0]);
    fey = fabs(e2[1]);
    fez = fabs(e2[2]);

    AXISTEST_X2(e2[2],e2[1],fez,fey);
    AXISTEST_Y1(e2[2],e2[0],fez,fex);
    AXISTEST_Z12(e2[1],e2[0],fey,fex);

    // Bullet 1: 
    //   first test overlap in the {x,y,z}-directions 
    //   find min,max of the triangle each direction, and test for overlap in 
    //   that direction -- this is equivalent to testing a minimal AABB around 
    //   the triangle against the AABB 

    // test in X-direction 

    FINDMINMAX(v0[0],v1[0],v2[0],min,max);
    if (min >  boxhalfsize[0] || max < -boxhalfsize[0]) return false;

    // test in Y-direction
    
    FINDMINMAX(v0[1],v1[1],v2[1],min,max);
    if (min > boxhalfsize[1] || max < -boxhalfsize[1]) return false;

    // test in Z-direction 
    
    FINDMINMAX(v0[2],v1[2],v2[2],min,max);
    if (min > boxhalfsize[2] || max < -boxhalfsize[2]) return false;

    // Bullet 2: 
    //   test if the box intersects the plane of the triangle 
    //   compute plane equation of triangle: normal*x+d=0

    const double normal[3] = CROSS_PRODUCT(e0,e1);
    if (!planeBoxOverlap(normal,v0,boxhalfsize)) return false;       

    return true; // box and triangle overlaps

  }

#undef FINDMINMAX
#undef AXISTEST_X01
#undef AXISTEST_X2
#undef AXISTEST_Y02
#undef AXISTEST_Y1
#undef AXISTEST_Z12
#undef AXISTEST_Z0

  int diffuseVoxelData(const double zero = 1.0E-8,const double relax = 1.0,const int maxiter = 100,
      const bool b_verbose = false,const bool b_write = false) {

    if (mpi_rank == 0) 
      cout << "diffuseVoxelData()" << endl;

    map<const int,int> zoneCountMap;
    map<const int,int>::iterator zoneCountMap_iter;

    int iter = 0;
    int done = 0;
    while (done == 0) {
      ++iter;

      double my_res_max = 0.0;
      for (int ivbd = 0, limit = vbdVec.size(); ivbd < limit; ++ivbd) {
        const int IJK = vbdVec[ivbd];
        assert(vbd[IJK] != NULL); 

        // apply filter (to voxels not containing tri, flag = 0)...
        for (int k = N_GHOSTS; k < BLOCK_SIZE_G-N_GHOSTS; ++k) {
          for (int j = N_GHOSTS; j < BLOCK_SIZE_G-N_GHOSTS; ++j) {
            for (int i = N_GHOSTS; i < BLOCK_SIZE_G-N_GHOSTS; ++i) {
              const int ijk = (k*BLOCK_SIZE_G+j)*BLOCK_SIZE_G+i;
              if ((vbd[IJK]->flag[ijk]&SURFACE_BIT) == 0) {
                assert(zoneCountMap.empty());

                // its data becomes the simple average of its nbrs non-HUGE_DATA data.
                // update in place to act like a Gauss-Seidel. its zone becomes the mode of its nbrs zones...
                double data_new = 0.0;
                int count = 0;
                for (int ii = 0; ii < 3; ++ii) {
                  for (int jj = 0; jj < 3; ++jj) {
                    for (int kk = 0; kk < 3; ++kk) {
                      // don't include yourself...
                      if ((ii == 1)&&(jj == 1)&&(kk == 1))
                        continue;
                      const int ijk_nbr = ((k+(kk-1))*BLOCK_SIZE_G+j+(jj-1))*BLOCK_SIZE_G+i+(ii-1);
                      if (vbd[IJK]->data[ijk_nbr] != HUGE_VAL) {
                        data_new += vbd[IJK]->data[ijk_nbr];
                        ++count;
                      }
                      if (vbd[IJK]->zone[ijk_nbr] >= 0) {
                        zoneCountMap_iter = zoneCountMap.find(vbd[IJK]->zone[ijk_nbr]);
                        if (zoneCountMap_iter == zoneCountMap.end()) {
                          zoneCountMap[vbd[IJK]->zone[ijk_nbr]] = 1;
                        }
                        else {
                          ++zoneCountMap_iter->second;
                        }
                      }
                    }
                  }
                }
                if (count > 0) {
                  if (vbd[IJK]->data[ijk] != HUGE_VAL) {
                    data_new = relax*data_new/(double)count + (1.0-relax)*vbd[IJK]->data[ijk];
                  }
                  else { 
                    data_new = data_new/(double)count;
                  }
                  my_res_max = max(my_res_max,fabs(vbd[IJK]->data[ijk]-data_new));
                  vbd[IJK]->data[ijk] = data_new;
                }
                else {
                  my_res_max = HUGE_VAL;
                }
                int max_zone_count = 0;
                for (zoneCountMap_iter = zoneCountMap.begin(); zoneCountMap_iter != zoneCountMap.end(); ++zoneCountMap_iter) {
                  if (zoneCountMap_iter->second > max_zone_count) {
                    max_zone_count = zoneCountMap_iter->second;
                    vbd[IJK]->zone[ijk] = zoneCountMap_iter->first;
                  }
                }
                zoneCountMap.clear();
              }
            }
          }
        }
      }

      double res_max;
      MPI_Reduce(&my_res_max, &res_max, 1, MPI_DOUBLE, MPI_MAX, 0, mpi_comm);
      if (mpi_rank == 0) {
        res_max *= inv_dx; // assume data is distance...
        if ((b_verbose)||(iter > maxiter/2))
          cout << " > diffuseVoxelData iter, res_max: " << iter << " " << res_max << endl;
        if (res_max < zero) {
          done = 1;
        }
        else if (iter > maxiter) {
          cout << " > Warning: diffuseVoxelData did not converge after " << maxiter << " iters, res_max: " << res_max << endl;
          done = 2;
        }
      }
      MPI_Bcast(&done, 1, MPI_INT, 0, mpi_comm);

      updateVoxelData(UPDATE_DATA|UPDATE_ZONE,true);

    }

    if (b_write)
      writeVoxelDataToTecplot("diffuse.dat");

    return( done == 1 );
  }

  int redistanceVoxelData(const double zero = 1.0E-8,const double cfl = 0.5,const int maxiter = 100,
      const bool b_verbose = false,const bool b_write = false) {
    assert(N_GHOSTS >= 1);
    const double dt  = cfl*dx; // speed is at most 1

    // commented sections allow for you to switch b/w euler and rk3 abd weno5 with gudonov...

    if (mpi_rank == 0) 
      cout << "redistanceVoxelData()" << endl;

    int iter = 0;
    int done = 0;
    while (done == 0) {
      ++iter;

      // need to fill all buffers, communicate with closest periodic nbr in each direction...
      updateVoxelDataClosestPeriodicFaceNbr();

      for (int ivbd = 0, limit = vbdVec.size(); ivbd < limit; ++ivbd) {
        const int IJK = vbdVec[ivbd];
        assert(vbd[IJK] != NULL); 
        for (int ijk = 0; ijk < BLOCK_SIZE3_G; ++ijk) {
          vbd[IJK]->data00[ijk] = vbd[IJK]->data[ijk]; 
        }
      }

      double my_res_max = 0.0;
      for (int irk = 0; irk < 3; ++irk) {

        for (int ivbd = 0, limit = vbdVec.size(); ivbd < limit; ++ivbd) {
          const int IJK = vbdVec[ivbd];
          assert(vbd[IJK] != NULL); 
          for (int ijk = 0; ijk < BLOCK_SIZE3_G; ++ijk) {
            switch (irk) { 
              case 0:
              case 1:
                vbd[IJK]->data0[ijk] = vbd[IJK]->data[ijk];
                break;
              case 2:
                vbd[IJK]->data0[ijk] = 0.75*vbd[IJK]->data00[ijk] + 0.25*vbd[IJK]->data[ijk];
                break;
              default:
                assert(0);
            }
          }
        }

        for (int ivbd = 0, limit = vbdVec.size(); ivbd < limit; ++ivbd) {
          const int IJK = vbdVec[ivbd];
          assert(vbd[IJK] != NULL); 

          for (int k = N_GHOSTS; k < BLOCK_SIZE_G-N_GHOSTS; ++k) {
            for (int j = N_GHOSTS; j < BLOCK_SIZE_G-N_GHOSTS; ++j) {
              for (int i = N_GHOSTS; i < BLOCK_SIZE_G-N_GHOSTS; ++i) {
                const int ijk = (k*BLOCK_SIZE_G+j)*BLOCK_SIZE_G+i;
                if (vbd[IJK]->flag[ijk]&SURFACE_BIT)
                  continue;

                double ddx_mp[2][3]; // minus/plus,direction
                ddx_mp[0][0] = (vbd[IJK]->data0[(k*BLOCK_SIZE_G+j)*BLOCK_SIZE_G+i]-
                                   vbd[IJK]->data0[(k*BLOCK_SIZE_G+j)*BLOCK_SIZE_G+i-1])*inv_dx;
                ddx_mp[0][1] = (vbd[IJK]->data0[(k*BLOCK_SIZE_G+j)*BLOCK_SIZE_G+i]-
                                   vbd[IJK]->data0[(k*BLOCK_SIZE_G+j-1)*BLOCK_SIZE_G+i])*inv_dx;
                ddx_mp[0][2] = (vbd[IJK]->data0[(k*BLOCK_SIZE_G+j)*BLOCK_SIZE_G+i]-
                                   vbd[IJK]->data0[((k-1)*BLOCK_SIZE_G+j)*BLOCK_SIZE_G+i])*inv_dx;
                ddx_mp[1][0] = (vbd[IJK]->data0[(k*BLOCK_SIZE_G+j)*BLOCK_SIZE_G+i+1]-
                                   vbd[IJK]->data0[(k*BLOCK_SIZE_G+j)*BLOCK_SIZE_G+i])*inv_dx;
                ddx_mp[1][1] = (vbd[IJK]->data0[(k*BLOCK_SIZE_G+j+1)*BLOCK_SIZE_G+i]-
                                   vbd[IJK]->data0[(k*BLOCK_SIZE_G+j)*BLOCK_SIZE_G+i])*inv_dx;
                ddx_mp[1][2] = (vbd[IJK]->data0[((k+1)*BLOCK_SIZE_G+j)*BLOCK_SIZE_G+i]-
                                   vbd[IJK]->data0[(k*BLOCK_SIZE_G+j)*BLOCK_SIZE_G+i])*inv_dx;

                double mag_grad2;
                if (vbd[IJK]->data0[ijk] <= 0.0) {
                  mag_grad2 = max( min(ddx_mp[0][0],0.0)*min(ddx_mp[0][0],0.0) , max(ddx_mp[1][0],0.0)*max(ddx_mp[1][0],0.0) )+
                              max( min(ddx_mp[0][1],0.0)*min(ddx_mp[0][1],0.0) , max(ddx_mp[1][1],0.0)*max(ddx_mp[1][1],0.0) )+
                              max( min(ddx_mp[0][2],0.0)*min(ddx_mp[0][2],0.0) , max(ddx_mp[1][2],0.0)*max(ddx_mp[1][2],0.0) );
                }
                else {
                  mag_grad2 = max( max(ddx_mp[0][0],0.0)*max(ddx_mp[0][0],0.0) , min(ddx_mp[1][0],0.0)*min(ddx_mp[1][0],0.0) )+
                              max( max(ddx_mp[0][1],0.0)*max(ddx_mp[0][1],0.0) , min(ddx_mp[1][1],0.0)*min(ddx_mp[1][1],0.0) )+
                              max( max(ddx_mp[0][2],0.0)*max(ddx_mp[0][2],0.0) , min(ddx_mp[1][2],0.0)*min(ddx_mp[1][2],0.0) );
                }

                vbd[IJK]->data[ijk] = vbd[IJK]->data0[ijk]*
                  (1.0 - dt*(sqrt(mag_grad2)-1.0)/sqrt(vbd[IJK]->data0[ijk]*vbd[IJK]->data0[ijk]+dx*dx*mag_grad2));

              }
            }
          }
        }

        updateVoxelDataClosestPeriodicFaceNbr();
      }

      for (int ivbd = 0, limit = vbdVec.size(); ivbd < limit; ++ivbd) {
        const int IJK = vbdVec[ivbd];
        assert(vbd[IJK] != NULL); 
        for (int ijk = 0; ijk < BLOCK_SIZE3_G; ++ijk) {
          vbd[IJK]->data[ijk] = 1.0/3.0*vbd[IJK]->data00[ijk]+2.0/3.0*vbd[IJK]->data[ijk]; 
          my_res_max = max(my_res_max,fabs(vbd[IJK]->data[ijk]-vbd[IJK]->data00[ijk]));
        }
      }

      double res_max;
      MPI_Reduce(&my_res_max, &res_max, 1, MPI_DOUBLE, MPI_MAX, 0, mpi_comm);
      if (mpi_rank == 0) {
        res_max *= inv_dx; // assume data is distance...
        if ((b_verbose)||(iter > maxiter/2))
          cout << " > redistanceVoxelData iter, res_max: " << iter << " " << res_max << endl;
        if (res_max < zero) {
          done = 1;
        }
        else if (iter > maxiter) {
          cout << " > Warning: redistanceVoxelData did not converge after " << maxiter << " iters, res_max: " << res_max << endl;
          done = 2;
        }
      }
      MPI_Bcast(&done, 1, MPI_INT, 0, mpi_comm);

    }

    // fill all the ghosts when you exit...
    updateVoxelData(UPDATE_DATA,true); // get in all ghosts

    if (b_write)
      writeVoxelDataToTecplot("redistance.dat");

    return( done == 1 );
  }

  int redistanceVoxelDataWeno(const double zero = 1.0E-8,const double cfl = 0.5,const int maxiter = 100,
      const bool b_verbose = false,const bool b_write = false) {
    assert(N_GHOSTS >= 3);
    const double dt  = cfl*dx; // speed is at most 1

    // commented sections allow for you to switch b/w euler and rk3 abd weno5 with gudonov...

    if (mpi_rank == 0) 
      cout << "redistanceVoxelData()" << endl;

    int iter = 0;
    int done = 0;
    while (done == 0) {
      ++iter;

      // need to fill all buffers, communicate with closest periodic nbr in each direction...
      updateVoxelDataClosestPeriodicFaceNbr();

      for (int ivbd = 0, limit = vbdVec.size(); ivbd < limit; ++ivbd) {
        const int IJK = vbdVec[ivbd];
        assert(vbd[IJK] != NULL); 
        for (int ijk = 0; ijk < BLOCK_SIZE3_G; ++ijk) {
          vbd[IJK]->data00[ijk] = vbd[IJK]->data[ijk]; 
        }
      }

      double my_res_max = 0.0;
      for (int irk = 0; irk < 3; ++irk) {

        for (int ivbd = 0, limit = vbdVec.size(); ivbd < limit; ++ivbd) {
          const int IJK = vbdVec[ivbd];
          assert(vbd[IJK] != NULL); 
          for (int ijk = 0; ijk < BLOCK_SIZE3_G; ++ijk) {
            switch (irk) { 
              case 0:
              case 1:
                vbd[IJK]->data0[ijk] = vbd[IJK]->data[ijk];
                break;
              case 2:
                vbd[IJK]->data0[ijk] = 0.75*vbd[IJK]->data00[ijk] + 0.25*vbd[IJK]->data[ijk];
                break;
              default:
                assert(0);
            }
          }
        }

        for (int ivbd = 0, limit = vbdVec.size(); ivbd < limit; ++ivbd) {
          const int IJK = vbdVec[ivbd];
          assert(vbd[IJK] != NULL); 

          for (int k = N_GHOSTS; k < BLOCK_SIZE_G-N_GHOSTS; ++k) {
            for (int j = N_GHOSTS; j < BLOCK_SIZE_G-N_GHOSTS; ++j) {
              for (int i = N_GHOSTS; i < BLOCK_SIZE_G-N_GHOSTS; ++i) {
                const int ijk = (k*BLOCK_SIZE_G+j)*BLOCK_SIZE_G+i;
                if (vbd[IJK]->flag[ijk]&SURFACE_BIT)
                  continue;

                double d_mp[2][5][3]; // minus/plus,finite difference,direction
                for (int ind = 0; ind < 5; ++ind) {
                  d_mp[0][ind][0] = (vbd[IJK]->data0[(k*BLOCK_SIZE_G+j)*BLOCK_SIZE_G+i-2+ind]-
                                     vbd[IJK]->data0[(k*BLOCK_SIZE_G+j)*BLOCK_SIZE_G+i-3+ind])*inv_dx;
                  d_mp[0][ind][1] = (vbd[IJK]->data0[(k*BLOCK_SIZE_G+j-2+ind)*BLOCK_SIZE_G+i]-
                                     vbd[IJK]->data0[(k*BLOCK_SIZE_G+j-3+ind)*BLOCK_SIZE_G+i])*inv_dx;
                  d_mp[0][ind][2] = (vbd[IJK]->data0[((k-2+ind)*BLOCK_SIZE_G+j)*BLOCK_SIZE_G+i]-
                                     vbd[IJK]->data0[((k-3+ind)*BLOCK_SIZE_G+j)*BLOCK_SIZE_G+i])*inv_dx;
                }
                for (int ind = 0; ind < 5; ++ind) {
                  if (ind == 0) {
                    d_mp[1][ind][0] = (vbd[IJK]->data0[(k*BLOCK_SIZE_G+j)*BLOCK_SIZE_G+i-2+(5-ind)]-
                                       vbd[IJK]->data0[(k*BLOCK_SIZE_G+j)*BLOCK_SIZE_G+i-3+(5-ind)])*inv_dx;
                    d_mp[1][ind][1] = (vbd[IJK]->data0[(k*BLOCK_SIZE_G+j-2+(5-ind))*BLOCK_SIZE_G+i]-
                                       vbd[IJK]->data0[(k*BLOCK_SIZE_G+j-3+(5-ind))*BLOCK_SIZE_G+i])*inv_dx;
                    d_mp[1][ind][2] = (vbd[IJK]->data0[((k-2+(5-ind))*BLOCK_SIZE_G+j)*BLOCK_SIZE_G+i]-
                                       vbd[IJK]->data0[((k-3+(5-ind))*BLOCK_SIZE_G+j)*BLOCK_SIZE_G+i])*inv_dx;
                  }
                  else {
                    d_mp[1][ind][0] = d_mp[0][5-ind][0];
                    d_mp[1][ind][1] = d_mp[0][5-ind][1];
                    d_mp[1][ind][2] = d_mp[0][5-ind][2];
                  }
                }

                double ddx_mp[2][3];
                for (int ii = 0; ii < 2; ++ii) {
                  for (int jj = 0; jj < 3; ++jj) {
                    const double eps = 1.0E-6*MAX5(d_mp[ii][0][jj]*d_mp[ii][0][jj],d_mp[ii][1][jj]*d_mp[ii][1][jj],
                          d_mp[ii][2][jj]*d_mp[ii][2][jj],d_mp[ii][3][jj]*d_mp[ii][3][jj],d_mp[ii][4][jj]*d_mp[ii][4][jj])+1.0E-99;
                    double omega[3];
                    omega[0]  = 13.0/12.0*(    d_mp[ii][0][jj]-2.0*d_mp[ii][1][jj]+    d_mp[ii][2][jj])
                                         *(    d_mp[ii][0][jj]-2.0*d_mp[ii][1][jj]+    d_mp[ii][2][jj])+
                                     0.25*(    d_mp[ii][0][jj]-4.0*d_mp[ii][1][jj]+3.0*d_mp[ii][2][jj])
                                         *(    d_mp[ii][0][jj]-4.0*d_mp[ii][1][jj]+3.0*d_mp[ii][2][jj]);
                    omega[0] = 0.1/((omega[0]+eps)*(omega[0]+eps));
                    omega[1] = 13.0/12.0*(    d_mp[ii][1][jj]-2.0*d_mp[ii][2][jj]+    d_mp[ii][3][jj])
                                        *(    d_mp[ii][1][jj]-2.0*d_mp[ii][2][jj]+    d_mp[ii][3][jj])+
                                    0.25*(    d_mp[ii][1][jj]                    -    d_mp[ii][3][jj])
                                        *(    d_mp[ii][1][jj]                    -    d_mp[ii][3][jj]);
                    omega[1] = 0.6/((omega[1]+eps)*(omega[1]+eps));
                    omega[2] = 13.0/12.0*(    d_mp[ii][2][jj]-2.0*d_mp[ii][3][jj]+    d_mp[ii][4][jj])
                                        *(    d_mp[ii][2][jj]-2.0*d_mp[ii][3][jj]+    d_mp[ii][4][jj])+
                                    0.25*(3.0*d_mp[ii][2][jj]-4.0*d_mp[ii][3][jj]+    d_mp[ii][4][jj])
                                        *(3.0*d_mp[ii][2][jj]-4.0*d_mp[ii][3][jj]+    d_mp[ii][4][jj]);
                    omega[2] = 0.3/((omega[2]+eps)*(omega[2]+eps));
                    const double tot = omega[0]+omega[1]+omega[2]; 
                    assert(tot > 0.0);
                    omega[0] /= tot;
                    omega[1] /= tot;
                    omega[2] /= tot;
                    ddx_mp[ii][jj] = (omega[0]*( 2.0*d_mp[ii][0][jj]-7.0*d_mp[ii][1][jj]+11.0*d_mp[ii][2][jj])+
                                      omega[1]*(-1.0*d_mp[ii][1][jj]+5.0*d_mp[ii][2][jj]+ 2.0*d_mp[ii][3][jj])+
                                      omega[2]*( 2.0*d_mp[ii][2][jj]+5.0*d_mp[ii][3][jj]- 1.0*d_mp[ii][4][jj]))/6.0;
                  }
                }

                double mag_grad2;
                if (vbd[IJK]->data0[ijk] <= 0.0) {
                  mag_grad2 = max( min(ddx_mp[0][0],0.0)*min(ddx_mp[0][0],0.0) , max(ddx_mp[1][0],0.0)*max(ddx_mp[1][0],0.0) )+
                              max( min(ddx_mp[0][1],0.0)*min(ddx_mp[0][1],0.0) , max(ddx_mp[1][1],0.0)*max(ddx_mp[1][1],0.0) )+
                              max( min(ddx_mp[0][2],0.0)*min(ddx_mp[0][2],0.0) , max(ddx_mp[1][2],0.0)*max(ddx_mp[1][2],0.0) );
                }
                else {
                  mag_grad2 = max( max(ddx_mp[0][0],0.0)*max(ddx_mp[0][0],0.0) , min(ddx_mp[1][0],0.0)*min(ddx_mp[1][0],0.0) )+
                              max( max(ddx_mp[0][1],0.0)*max(ddx_mp[0][1],0.0) , min(ddx_mp[1][1],0.0)*min(ddx_mp[1][1],0.0) )+
                              max( max(ddx_mp[0][2],0.0)*max(ddx_mp[0][2],0.0) , min(ddx_mp[1][2],0.0)*min(ddx_mp[1][2],0.0) );
                }

                vbd[IJK]->data[ijk] = vbd[IJK]->data0[ijk]*
                  (1.0 - dt*(sqrt(mag_grad2)-1.0)/sqrt(vbd[IJK]->data0[ijk]*vbd[IJK]->data0[ijk]+dx*dx*mag_grad2));

              }
            }
          }
        }

        updateVoxelDataClosestPeriodicFaceNbr();
      }

      for (int ivbd = 0, limit = vbdVec.size(); ivbd < limit; ++ivbd) {
        const int IJK = vbdVec[ivbd];
        assert(vbd[IJK] != NULL); 
        for (int ijk = 0; ijk < BLOCK_SIZE3_G; ++ijk) {
          vbd[IJK]->data[ijk] = 1.0/3.0*vbd[IJK]->data00[ijk]+2.0/3.0*vbd[IJK]->data[ijk]; 
          my_res_max = max(my_res_max,fabs(vbd[IJK]->data[ijk]-vbd[IJK]->data00[ijk]));
        }
      }

      double res_max;
      MPI_Reduce(&my_res_max, &res_max, 1, MPI_DOUBLE, MPI_MAX, 0, mpi_comm);
      if (mpi_rank == 0) {
        res_max *= inv_dx; // assume data is distance...
        if ((b_verbose)||(iter > maxiter/2))
          cout << " > redistanceVoxelData iter, res_max: " << iter << " " << res_max << endl;
        if (res_max < zero) {
          done = 1;
        }
        else if (iter > maxiter) {
          cout << " > Warning: redistanceVoxelData did not converge after " << maxiter << " iters, res_max: " << res_max << endl;
          done = 2;
        }
      }
      MPI_Bcast(&done, 1, MPI_INT, 0, mpi_comm);

    }

    // fill all the ghosts when you exit...
    updateVoxelData(UPDATE_DATA,true); // get in all ghosts

    if (b_write)
      writeVoxelDataToTecplot("redistance.dat");

    return( done == 1 );
  }

  int orientVoxelData(const int maxiter = 1000,const double alpha = 1.0,const double beta = 0.5,
      const bool b_verbose = false,const bool b_write = false) {
    assert(N_GHOSTS >= 1);

    if (mpi_rank == 0) 
      cout << "orientVoxelData()" << endl;

    // consensus of signed distance field...

    int iter = 0;
    int done = 0;
    while (done == 0) {
      ++iter;

      for (int ivbd = 0, limit = vbdVec.size(); ivbd < limit; ++ivbd) {
        const int IJK = vbdVec[ivbd];
        assert(vbd[IJK] != NULL); 
        for (int ijk = 0; ijk < BLOCK_SIZE3_G; ++ijk) {
          vbd[IJK]->data0[ijk] = vbd[IJK]->data[ijk];
        }
      }

      int8 my_nflip = 0;
      for (int ivbd = 0, limit = vbdVec.size(); ivbd < limit; ++ivbd) {
        const int IJK = vbdVec[ivbd];
        assert(vbd[IJK] != NULL); 

        // apply filter (to voxels not containing tri, flag = 0)...
        for (int k = N_GHOSTS; k < BLOCK_SIZE_G-N_GHOSTS; ++k) {
          for (int j = N_GHOSTS; j < BLOCK_SIZE_G-N_GHOSTS; ++j) {
            for (int i = N_GHOSTS; i < BLOCK_SIZE_G-N_GHOSTS; ++i) {
              const int ijk = (k*BLOCK_SIZE_G+j)*BLOCK_SIZE_G+i;
              //if (vbd[IJK]->flag[ijk]&X_RAY_IN_BIT) {
              if (vbd[IJK]->data0[ijk] != HUGE_VAL) {
                int n1 = 0,n2 = 0,n3 = 0,n4 = 0;
                for (int ii = 0; ii < 3; ++ii) {
                  for (int jj = 0; jj < 3; ++jj) {
                    for (int kk = 0; kk < 3; ++kk) {
                      // don't include yourself...
                      if ((ii == 1)&&(jj == 1)&&(kk == 1))
                        continue;
                      const int ijk_nbr = ((k+(kk-1))*BLOCK_SIZE_G+j+(jj-1))*BLOCK_SIZE_G+i+(ii-1);
                      if (vbd[IJK]->data0[ijk_nbr] != HUGE_VAL) {
                        const double delta = dx*sqrt(double((ii-1)*(ii-1)+(jj-1)*(jj-1)+(kk-1)*(kk-1)));
                        if (fabs(vbd[IJK]->data0[ijk_nbr]-vbd[IJK]->data0[ijk]) <= (alpha-tol)*delta) 
                          ++n1;
                        else 
                          ++n2;
                        if (fabs(-vbd[IJK]->data0[ijk_nbr]-vbd[IJK]->data0[ijk]) <= (alpha-tol)*delta) 
                          ++n3;
                        else 
                          ++n4;
                      }
                    }
                  }
                }
                if ((n2+n3) > beta*(n1+n2+n3+n4)) {
                  vbd[IJK]->data[ijk] *= -1;
                  ++my_nflip;
                }
              }
            }
          }
        }
      }

      int8 nflip;
      MPI_Reduce(&my_nflip,&nflip,1,MPI_INT8,MPI_SUM,0,mpi_comm);
      if (mpi_rank == 0) {
        if ((b_verbose)||(iter > maxiter/2))
          cout << " > orientVoxelData iter, nflip: " << iter << " " << nflip << endl;
        if (nflip == 0) {
          done = 1;
        }
        else if (iter > maxiter) {
          cout << " > Warning: orientVoxelData did not converge after " << maxiter << " iters, nflip: " << nflip << endl;
          done = 2;
        }
      }
      MPI_Bcast(&done,1,MPI_INT,0,mpi_comm);

      updateVoxelData(UPDATE_DATA,true);
    }

    if (b_write)
      writeVoxelDataToTecplot("flip.dat");

    return( done == 1 );
  }

  // finds your closet periodic face nbr in each direction. assumes we are communicating data only... 
  void updateVoxelDataClosestPeriodicFaceNbr() {

    // send your internal data to your nbrs ghosts if they exist...

    int* send_count = new int[mpi_size];
    int* send_disp  = new int[mpi_size];
    FOR_RANK send_count[rank] = 0;
    int* send_buf_int = NULL;
    double* send_buf_double = NULL;
    for (int iter = 0; iter < 2; ++iter) {
      for (int ivbd = 0, limit = vbdVec.size(); ivbd < limit; ++ivbd) {
        const int IJK = vbdVec[ivbd];
        assert(vbd[IJK] != NULL); 
        const int I = IJK%nI;
        const int J = (IJK/nI)%nJ;
        const int K = IJK/(nI*nJ);
        assert(IJK == (K*nJ+J)*nI+I);

        // top...
        
        {
          int Kp = -1,IJKp = -1;
          for (Kp = K+1; Kp < nK; ++Kp) {
            IJKp = (Kp*nJ+J)*nI+I;
            assert((IJKp >= 0)&&(IJKp < nI*nJ*nK));
            if (vbd_rank[IJKp] >= 0) 
              break;
          }
          if (Kp == nK) {
            for (Kp = 0; Kp <= K; ++Kp) {
              IJKp = (Kp*nJ+J)*nI+I;
              assert((IJKp >= 0)&&(IJKp < nI*nJ*nK));
              if (vbd_rank[IJKp] >= 0) 
                break;
            }
          }
          assert(vbd_rank[IJKp] >= 0);
          const int rank = vbd_rank[IJKp];
          if (iter == 0) {
            send_count[rank] += N_GHOSTS*BLOCK_SIZE*BLOCK_SIZE;
          }
          else {
            for (int k = BLOCK_SIZE; k < N_GHOSTS+BLOCK_SIZE; ++k) {
              for (int j = N_GHOSTS; j < BLOCK_SIZE_G-N_GHOSTS; ++j) {
                for (int i = N_GHOSTS; i < BLOCK_SIZE_G-N_GHOSTS; ++i) {
                  const int ijk = (k*BLOCK_SIZE_G+j)*BLOCK_SIZE_G+i;
                  const int ijkp = ((k-BLOCK_SIZE)*BLOCK_SIZE_G+j)*BLOCK_SIZE_G+i;
                  send_buf_int[send_disp[rank]*2  ] = IJKp;
                  send_buf_int[send_disp[rank]*2+1] = ijkp;
                  send_buf_double[send_disp[rank]] = vbd[IJK]->data[ijk];
                  ++send_disp[rank];
                }
              }
            }
          }
        }

        // bottom...
        
        {
          int Km = -1,IJKm = -1;
          for (Km = K-1; Km >= 0; --Km) {
            IJKm = (Km*nJ+J)*nI+I;
            assert((IJKm >= 0)&&(IJKm < nI*nJ*nK));
            if (vbd_rank[IJKm] >= 0)
              break;
          }
          if (Km == -1) {
            for (Km = nK-1; Km >= K; --Km) {
              IJKm = (Km*nJ+J)*nI+I;
              assert((IJKm >= 0)&&(IJKm < nI*nJ*nK));
              if (vbd_rank[IJKm] >= 0)
                break;
            }
          }
          assert(vbd_rank[IJKm] >= 0);
          const int rank = vbd_rank[IJKm];
          if (iter == 0) {
            send_count[rank] += N_GHOSTS*BLOCK_SIZE*BLOCK_SIZE;
          }
          else {
            for (int k = N_GHOSTS; k < 2*N_GHOSTS; ++k) {
              for (int j = N_GHOSTS; j < BLOCK_SIZE_G-N_GHOSTS; ++j) {
                for (int i = N_GHOSTS; i < BLOCK_SIZE_G-N_GHOSTS; ++i) {
                  const int ijk = (k*BLOCK_SIZE_G+j)*BLOCK_SIZE_G+i;
                  const int ijkm = ((k+BLOCK_SIZE)*BLOCK_SIZE_G+j)*BLOCK_SIZE_G+i;
                  send_buf_int[send_disp[rank]*2  ] = IJKm;
                  send_buf_int[send_disp[rank]*2+1] = ijkm;
                  send_buf_double[send_disp[rank]] = vbd[IJK]->data[ijk];
                  ++send_disp[rank];
                }
              }
            }
          }
        }

        // north...
        
        {
          int Jp = -1,IJpK = -1;
          for (Jp = J+1; Jp < nJ; ++Jp) {
            IJpK = (K*nJ+Jp)*nI+I;
            assert((IJpK >= 0)&&(IJpK < nI*nJ*nK));
            if (vbd_rank[IJpK] >= 0) 
              break;
          }
          if (Jp == nJ) {
            for (Jp = 0; Jp <= J; ++Jp) {
              IJpK = (K*nJ+Jp)*nI+I;
              assert((IJpK >= 0)&&(IJpK < nI*nJ*nK));
              if (vbd_rank[IJpK] >= 0) 
                break;
            }
          }
          assert(vbd_rank[IJpK] >= 0);
          const int rank = vbd_rank[IJpK];
          if (iter == 0) {
            send_count[rank] += N_GHOSTS*BLOCK_SIZE*BLOCK_SIZE;
          }
          else {
            for (int k = N_GHOSTS; k < BLOCK_SIZE_G-N_GHOSTS; ++k) {
              for (int j = BLOCK_SIZE; j < N_GHOSTS+BLOCK_SIZE; ++j) {
                for (int i = N_GHOSTS; i < BLOCK_SIZE_G-N_GHOSTS; ++i) {
                  const int ijk = (k*BLOCK_SIZE_G+j)*BLOCK_SIZE_G+i;
                  const int ijpk = (k*BLOCK_SIZE_G+j-BLOCK_SIZE)*BLOCK_SIZE_G+i; 
                  send_buf_int[send_disp[rank]*2  ] = IJpK;
                  send_buf_int[send_disp[rank]*2+1] = ijpk;
                  send_buf_double[send_disp[rank]] = vbd[IJK]->data[ijk];
                  ++send_disp[rank];
                }
              }
            }
          }
        }

        // south...
        
        {
          int Jm = -1,IJmK = -1;
          for (Jm = J-1; Jm >= 0; --Jm) {
            IJmK = (K*nJ+Jm)*nI+I;
            assert((IJmK >= 0)&&(IJmK < nI*nJ*nK));
            if (vbd_rank[IJmK] >= 0)
              break;
          }
          if (Jm == -1) {
            for (Jm = nJ-1; Jm >= J; --Jm) {
              IJmK = (K*nJ+Jm)*nI+I;
              assert((IJmK >= 0)&&(IJmK < nI*nJ*nK));
              if (vbd_rank[IJmK] >= 0)
                break;
            }
          }
          assert(vbd_rank[IJmK] >= 0);
          const int rank = vbd_rank[IJmK];
          if (iter == 0) {
            send_count[rank] += N_GHOSTS*BLOCK_SIZE*BLOCK_SIZE;
          }
          else {
            for (int k = N_GHOSTS; k < BLOCK_SIZE_G-N_GHOSTS; ++k) {
              for (int j = N_GHOSTS; j < 2*N_GHOSTS; ++j) {
                for (int i = N_GHOSTS; i < BLOCK_SIZE_G-N_GHOSTS; ++i) {
                  const int ijk = (k*BLOCK_SIZE_G+j)*BLOCK_SIZE_G+i;
                  const int ijmk = (k*BLOCK_SIZE_G+j+BLOCK_SIZE)*BLOCK_SIZE_G+i;
                  send_buf_int[send_disp[rank]*2  ] = IJmK;
                  send_buf_int[send_disp[rank]*2+1] = ijmk;
                  send_buf_double[send_disp[rank]] = vbd[IJK]->data[ijk];
                  ++send_disp[rank];
                }
              }
            }
          }
        }

        // east...
        
        {
          int Ip = -1,IpJK = -1;
          for (Ip = I+1; Ip < nI; ++Ip) {
            IpJK = (K*nJ+J)*nI+Ip;
            assert((IpJK >= 0)&&(IpJK < nI*nJ*nK));
            if (vbd_rank[IpJK] >= 0) 
              break;
          }
          if (Ip == nI) {
            for (Ip = 0; Ip <= I; ++Ip) {
              IpJK = (K*nJ+J)*nI+Ip;
              assert((IpJK >= 0)&&(IpJK < nI*nJ*nK));
              if (vbd_rank[IpJK] >= 0) 
                break;
            }
          }
          assert(vbd_rank[IpJK] >= 0);
          const int rank = vbd_rank[IpJK];
          if (iter == 0) {
            send_count[rank] += N_GHOSTS*BLOCK_SIZE*BLOCK_SIZE;
          }
          else {
            for (int k = N_GHOSTS; k < BLOCK_SIZE_G-N_GHOSTS; ++k) {
              for (int j = N_GHOSTS; j < BLOCK_SIZE_G-N_GHOSTS; ++j) {
                for (int i = BLOCK_SIZE; i < N_GHOSTS+BLOCK_SIZE; ++i) {
                  const int ijk = (k*BLOCK_SIZE_G+j)*BLOCK_SIZE_G+i;
                  const int ipjk = (k*BLOCK_SIZE_G+j)*BLOCK_SIZE_G+i-BLOCK_SIZE;
                  send_buf_int[send_disp[rank]*2  ] = IpJK;
                  send_buf_int[send_disp[rank]*2+1] = ipjk;
                  send_buf_double[send_disp[rank]] = vbd[IJK]->data[ijk];
                  ++send_disp[rank];
                }
              }
            }
          }
        }

        // west...
   
        {
          int Im = -1,ImJK = -1;
          for (Im = I-1; Im >= 0; --Im) {
            ImJK = (K*nJ+J)*nI+Im;
            assert((ImJK >= 0)&&(ImJK < nI*nJ*nK));
            if (vbd_rank[ImJK] >= 0)
              break;
          }
          if (Im == -1) {
            for (Im = nI-1; Im >= I; --Im) {
              ImJK = (K*nJ+J)*nI+Im;
              assert((ImJK >= 0)&&(ImJK < nI*nJ*nK));
              if (vbd_rank[ImJK] >= 0)
                break;
            }
          }
          assert(vbd_rank[ImJK] >= 0);
          const int rank = vbd_rank[ImJK];
          if (iter == 0) {
            send_count[rank] += N_GHOSTS*BLOCK_SIZE*BLOCK_SIZE;
          }
          else {
            for (int k = N_GHOSTS; k < BLOCK_SIZE_G-N_GHOSTS; ++k) {
              for (int j = N_GHOSTS; j < BLOCK_SIZE_G-N_GHOSTS; ++j) {
                for (int i = N_GHOSTS; i < 2*N_GHOSTS; ++i) {
                  const int ijk = (k*BLOCK_SIZE_G+j)*BLOCK_SIZE_G+i;
                  const int imjk = (k*BLOCK_SIZE_G+j)*BLOCK_SIZE_G+i+BLOCK_SIZE;
                  send_buf_int[send_disp[rank]*2  ] = ImJK;
                  send_buf_int[send_disp[rank]*2+1] = imjk;
                  send_buf_double[send_disp[rank]] = vbd[IJK]->data[ijk];
                  ++send_disp[rank];
                }
              }
            }
          }
        }
      }

      send_disp[0] = 0;
      for (int rank = 1; rank < mpi_size; ++rank)
        send_disp[rank] = send_count[rank-1] + send_disp[rank-1];

      if (iter == 0) {
        const int send_count_sum = send_disp[mpi_size-1] + send_count[mpi_size-1];
        send_buf_int = new int[2*send_count_sum];
        send_buf_double = new double[send_count_sum];
      }

    }

    // exchange...

    int *recv_count = new int[mpi_size];
    MPI_Alltoall(send_count,1,MPI_INT,recv_count,1,MPI_INT,mpi_comm);

    int *recv_disp = new int[mpi_size];
    recv_disp[0] = 0;
    for (int rank = 1; rank < mpi_size; ++rank)
      recv_disp[rank] = recv_count[rank-1] + recv_disp[rank-1];
    const int recv_count_sum = recv_disp[mpi_size-1] + recv_count[mpi_size-1];

    double * recv_buf_double = new double[recv_count_sum];
    MPI_Alltoallv(send_buf_double,send_count,send_disp,MPI_DOUBLE,
                  recv_buf_double,recv_count,recv_disp,MPI_DOUBLE,mpi_comm);
    delete[] send_buf_double; 

    FOR_RANK {
      send_count[rank] *= 2;
      send_disp[rank] *= 2;
      recv_count[rank] *= 2;
      recv_disp[rank] *= 2;
    }

    int * recv_buf_int = new int[2*recv_count_sum];
    MPI_Alltoallv(send_buf_int,send_count,send_disp,MPI_INT,
                  recv_buf_int,recv_count,recv_disp,MPI_INT,mpi_comm);
    delete[] send_buf_int; 
    delete[] send_count;
    delete[] send_disp;
    delete[] recv_count;
    delete[] recv_disp;
    
    // now set our data...

    for (int irecv = 0; irecv < recv_count_sum; ++irecv) {
      const int IJK = recv_buf_int[irecv*2  ]; assert(vbd_rank[IJK] == mpi_rank);
      const int ijk = recv_buf_int[irecv*2+1];
      assert(vbd[IJK] != NULL);
      vbd[IJK]->data[ijk] = recv_buf_double[irecv];
    }
    delete[] recv_buf_int;
    delete[] recv_buf_double;
  }
};

#endif 
