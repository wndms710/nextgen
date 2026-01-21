#include "StitchNS.hpp"
#include "CuttableVoronoiData.hpp"
#include "Prcomm.hpp"
#include "WriteImageData.hpp"

// exact copy from Move.hpp, so be careful...

#define PART_MASK_BITS  65535
#define REBUILD_BIT    (1<<16)
#define NEW_BIT        (1<<17)
#define ACTIVE_BIT     (1<<18)
#define ACTIVE0_BIT    (1<<19)
#define TMP1_BIT       (1<<20)
#define TMP2_BIT       (1<<21)
#define TMP3_BIT       (1<<22)
#define REBUILD2_BIT   (1<<23)

#define ORPHAN_FACE_OFFSET 10000 // use to be 1000

// 4 types of layers for now, may want to adjust their widths; currently uniform
// 8 bits for surface layers, resolution transition layers, periodic zone adacent layers, and part-boundary layers
#define SMOOTH_FLAG_MASK 255
#define SMOOTH_FLAG_S_OFFSET 0
#define SMOOTH_FLAG_T_OFFSET 8
#define SMOOTH_FLAG_P_OFFSET 16
#define SMOOTH_FLAG_PB_OFFSET 24

#define ZERO_LOCAL_FACE_BIT      (1<<0)
#define NEAR_ZERO_LOCAL_FACE_BIT (1<<1)
#define VALID_LOCAL_FACE_BIT     (1<<2)
#define NO_NBR_FACE_BIT          (1<<3)
#define ZERO_NBR_FACE_BIT        (1<<4)
#define NEAR_ZERO_NBR_FACE_BIT   (1<<5)
#define VALID_NBR_FACE_BIT       (1<<6)

namespace StitchNS {

  HcpPointBuilder hcp;
  
  // process "SHOW" commands from killfile...
  bool b_show_all = false;
  bool b_show_touching = false;
  bool b_show_zones = false;
  vector<int> show_hiddenSubzonesVec; // different meanings depending on boolean above
  
  // lloyd iteration stuff...
  double (*xp0_and_delta0)[4] = NULL;
  double smooth_limit = 0.5;
  int nsmooth = 0;
  int smoothing_step = 0;
  bool b_force_smooth_all = false;
  uint smooth_mode_idata[4] = {3,  // BLAYER smoothing: 3 layer default
                               0,  // TRANSITION smoothing: default is off
                               1,  // PERIODIC smoothing: default is 1 cell within delta/2
                               3}; // PART_LAYER smoothing: default is 3 layers
  double crease_angle_degrees = 175;
  int8 rc_local_threshold = 0; // default is to rebalance all Voronoi points
  int8 rc_limit = 100000000;

  // maybe put in common namespace?...

  Adt<double> * hcp_x_vd_bbox_adt = NULL;
  Adt<double> * hcp_x_vd_adt = NULL;

  enum CvdReturnDataStatus {
    VD_UNKNOWN,
    VD_FAILED_SEED_AND_NO_NBRS,
    VD_FAILED_SEED,
    VD_FAILED_DELTA_TOO_SMALL,
    VD_PASSED,
  };

  class VdReturnData {
    public:
      // ints...
      int icv;
      int ned,nno,nfa,nbf,ngr;
      int (*nooed)[2];
      // when faoed is positive (or 0), look into this to get the associated part and surface tri
      // when faoed is negative, this is -1 indexing on the cvofa array to get the nbr...
      int (*faoed)[2];
      int * cvofa;
      int * bits_fa; // face bits
      int (*spbobf)[2]; // ist and part-bits of boundary face: uncondensed bf's, so always 1:1
      int * edogr_i;
      // doubles...
      double delta;
      double (*x_no)[3];
      double (*r2_fa); // max radius squared of any face node
      double *area2_fa; // squared area
      int orphan_chunk_data;
      VdReturnData() {
        icv = -1;
        nooed = NULL;
        faoed = NULL;
        cvofa = NULL;
        bits_fa = NULL;
        spbobf = NULL;
        edogr_i = NULL;
        x_no = NULL;
        r2_fa = NULL;
        area2_fa = NULL;
        orphan_chunk_data = -1;
      }
      // NOTE: DO NOT write a fancy destructor here PLEASE.
      // Leave destructor blank because of allocated arrays. Use a clear()
      // memory management model instead for the vector of VdReturnData...
      void clear() {
        icv = -1;
        nno = ned = nfa = ngr = nbf = 0;
        DELETE(nooed);
        DELETE(faoed);
        DELETE(cvofa);
        DELETE(bits_fa);
        DELETE(spbobf);
        DELETE(edogr_i);
        DELETE(x_no);
        DELETE(r2_fa);
        DELETE(area2_fa);
        orphan_chunk_data = -1;
      }
      bool faceIsActive(const int ifa) const {
        assert(bits_fa);
        return ((bits_fa[ifa]&VALID_LOCAL_FACE_BIT)&&(bits_fa[ifa]&VALID_NBR_FACE_BIT));
      }
      int getByteCount() const {
        return sizeof(VdReturnData) + (ngr+1+ned*4+nfa*2+nbf*2)*sizeof(int) + (nno*3+nfa*2)*sizeof(double);
      }
      bool hasPeriodicNbrs() const {
        for (int ifa = 0; ifa < nfa; ++ifa) {
          const int icv_nbr = cvofa[ifa];
          if (icv_nbr >= ncv) {
            int rank,bits,index;
            BitUtils::unpackRankBitsIndex(rank,bits,index,rbi_g[icv_nbr-ncv]);
            if (bits) return true;
          }
        }
        return false;
      }
      void writeTecplot(const string& filename) const {
        double x0[3] = { 0.0, 0.0, 0.0 };
        writeTecplot(filename,x0);
      }
      void writeTecplot(const string& filename,const double x0[3]) const {
        FILE * fp = fopen(filename.c_str(),"w");
        fprintf(fp,"TITLE = \"%s\"\n","debug.dat");
        fprintf(fp,"VARIABLES = \"X\"\n");
        fprintf(fp,"\"Y\"\n");
        fprintf(fp,"\"Z\"\n");
        fprintf(fp,"ZONE T=\"%s\"\n",filename.c_str());
        fprintf(fp,"N=%d, E=%d, F=FEPOINT, ET=TRIANGLE\n",nno,ned);
        for (int ino = 0; ino < nno; ++ino) {
          fprintf(fp,"%18.15e %18.15e %18.15e\n",x_no[ino][0]+x0[0],x_no[ino][1]+x0[1],x_no[ino][2]+x0[2]);
        }
        for (int ied = 0; ied < ned; ++ied) {
          fprintf(fp,"%d %d %d\n",nooed[ied][0]+1,nooed[ied][1]+1,nooed[ied][1]+1);
        }
        fclose(fp);
      }
      void extractGroupAsCvd(CuttableVoronoiData& cvd,const int igr) {
        assert(cvd.empty());
        assert((igr >= 0)&&(igr < ngr));
        assert(edogr_i);
        cvd.resize_ned(edogr_i[igr+1]-edogr_i[igr]);
        int ino_min = nno;
        int ino_max = 0;
        for (int ied = edogr_i[igr]; ied != edogr_i[igr+1]; ++ied) {
          ino_min = min(ino_min,min(nooed[ied][0],nooed[ied][1]));
          ino_max = max(ino_max,max(nooed[ied][0],nooed[ied][1]));
        }
        for (int ied = edogr_i[igr]; ied != edogr_i[igr+1]; ++ied) {
          const int ied_cvd = ied-edogr_i[igr];
          cvd.nooed[ied_cvd][0] = nooed[ied][0]-ino_min;
          cvd.nooed[ied_cvd][1] = nooed[ied][1]-ino_min;
          // note: positive values here reference ist's (actually an index in spbobf), and
          // negatives refer to internal nbrs (in cvofa, -1 indexed for now): may need to revisit
          cvd.faoed[ied_cvd][0] = faoed[ied][0];
          cvd.faoed[ied_cvd][1] = faoed[ied][1];
        }
        cvd.resize_nno(ino_max-ino_min+1);
        for (int ino = ino_min; ino <= ino_max; ++ino) {
          FOR_I3 cvd.x_no[ino-ino_min][i] = x_no[ino][i];
        }
      }
      void calcVolumeAndCentroidForGroup(double &volume,double centroid[3],const int igr) const {
        assert((igr >= 0)&&(igr < ngr));
        volume = 0.0;
        FOR_I3 centroid[i] = 0.0;
        map<const int,const double *> firstNodeMap;
        for (int ied = edogr_i[igr]; ied != edogr_i[igr+1]; ++ied) {
          FOR_I2 {
            const int ifa = faoed[ied][i];
            map<const int,const double *>::iterator iter = firstNodeMap.find(ifa);
            if (iter == firstNodeMap.end()) {
              firstNodeMap[ifa] = x_no[nooed[ied][i]];
            }
            else {
              const double * const x1 = x_no[nooed[ied][i]];
              const double * const x2 = x_no[nooed[ied][1-i]];
              if ((iter->second != x1)&&(iter->second != x2)) {
                const double this_vol = CROSS_DOT(iter->second,x1,x2);
                volume += this_vol;
                FOR_I3 centroid[i] += this_vol*(iter->second[i]+x1[i]+x2[i]); // average 4 tet corners: note one corner is (0,0,0)
              }
            }
          }
        }
        FOR_I3 centroid[i] /= volume*4.0;
        volume /= 6.0;
      }
      void removeTwoEdgeNodes() {

        cout << "removeTwoEdgeNodes()" << endl;

        int (*no_flag)[2] = new int[nno][2];

        for (int ino = 0; ino < nno; ++ino)
          FOR_I2 no_flag[ino][i] = -1;

        for (int ied = 0; ied < ned; ++ied) {
          FOR_I2 {
            const int ino = nooed[ied][i];
            if (no_flag[ino][0] == -1) {
              no_flag[ino][0] = ied;
            }
            else if (no_flag[ino][1] == -1) {
              no_flag[ino][1] = ied;
            }
            else {
              no_flag[ino][0] = no_flag[ino][1] = -2;
            }
          }
        }

        {

          writeTecplot("before.dat");


          FILE * fp = fopen("nodes.dat","w");
          for (int ino = 0; ino < nno; ++ino) {
            if ((no_flag[ino][0] >= 0)&&(no_flag[ino][1] >= 0)) {
              fprintf(fp,"%18.15e %18.15e %18.15e\n",x_no[ino][0],x_no[ino][1],x_no[ino][2]);
              cout << " faoed for edges touching ino: " << ino << endl;
              FOR_I2 {
                const int ied = no_flag[ino][i];
                cout << " > ied: " << ied << " faoed: " << faoed[ied][0] << " " << faoed[ied][1] << endl;
                MiscUtils::writeEdge(ied,x_no[nooed[ied][0]],x_no[nooed[ied][1]]);
              }
              cout << "take a look" << endl;
            }
          }
          fclose(fp);
          cout << "take a look at 2-edge nodes in nodes.dat" << endl;
          getchar();
        }

        int nno_new = 0;
        for (int ino = 0; ino < nno; ++ino) {
          if ((no_flag[ino][0] == -2)&&(no_flag[ino][1] == -2)) {
            // this node has 3 or more edges: keep, but re-index...
            const int ino_new = nno_new++;
            // use no_flag[ino][0] to store the new node index...
            no_flag[ino][0] = -ino_new-3;
            if (ino_new < ino) {
              FOR_I3 x_no[ino_new][i] = x_no[ino][i];
            }
          }
          else if ((no_flag[ino][0] >= 0)&&(no_flag[ino][1] >= 0)) {
            // remove the second edge, reconnecting with the first...
            const int ied0 = no_flag[ino][0];
            const int ied1 = no_flag[ino][1];
            if (nooed[ied0][0] == ino) {
              if (nooed[ied1][0] == ino) {
                assert(faoed[ied0][0] == faoed[ied1][1]);
                assert(faoed[ied0][1] == faoed[ied1][0]);
                nooed[ied0][0] = nooed[ied1][1];
                if (no_flag[nooed[ied1][1]][0] == ied1) {
                  no_flag[nooed[ied1][1]][0] = ied0;
                }
                else if (no_flag[nooed[ied1][1]][1] == ied1) {
                  no_flag[nooed[ied1][1]][1] = ied0;
                }
                nooed[ied1][0] = nooed[ied1][1] = -1;
              }
              else {
                assert(nooed[ied1][1] == ino);
                assert(faoed[ied0][0] == faoed[ied1][0]);
                assert(faoed[ied0][1] == faoed[ied1][1]);
                nooed[ied0][0] = nooed[ied1][0];
                if (no_flag[nooed[ied1][0]][0] == ied1) {
                  no_flag[nooed[ied1][0]][0] = ied0;
                }
                else if (no_flag[nooed[ied1][0]][1] == ied1) {
                  no_flag[nooed[ied1][0]][1] = ied0;
                }
                nooed[ied1][0] = nooed[ied1][1] = -1;
              }
            }
            else {
              assert(nooed[ied0][1] == ino);
              if (nooed[ied1][0] == ino) {
                assert(faoed[ied0][0] == faoed[ied1][0]);
                assert(faoed[ied0][1] == faoed[ied1][1]);
                nooed[ied0][1] = nooed[ied1][1];
                if (no_flag[nooed[ied1][1]][0] == ied1) {
                  no_flag[nooed[ied1][1]][0] = ied0;
                }
                else if (no_flag[nooed[ied1][1]][1] == ied1) {
                  no_flag[nooed[ied1][1]][1] = ied0;
                }
                nooed[ied1][0] = nooed[ied1][1] = -1;
              }
              else {
                assert(nooed[ied1][1] == ino);
                assert(faoed[ied0][0] == faoed[ied1][1]);
                assert(faoed[ied0][1] == faoed[ied1][0]);
                nooed[ied0][1] = nooed[ied1][0];
                if (no_flag[nooed[ied1][0]][0] == ied1) {
                  no_flag[nooed[ied1][0]][0] = ied0;
                }
                else if (no_flag[nooed[ied1][0]][1] == ied1) {
                  no_flag[nooed[ied1][0]][1] = ied0;
                }
                nooed[ied1][0] = nooed[ied1][1] = -1;
              }
            }
          }
          else {
            assert((no_flag[ino][0] == -1)&&(no_flag[ino][1] == -1));
          }
        }

        int ned_new = 0;
        for (int ied = 0; ied < ned; ++ied) {
          if ((nooed[ied][0] >= 0)&&(nooed[ied][1] >= 0)) {
            const int ied_new = ned_new++;
            FOR_I2 {
              const int ino = nooed[ied][i];
              assert(no_flag[ino][0] <= -3);
              nooed[ied_new][i] = -no_flag[ino][0]-3;
            }
            faoed[ied_new][0] = faoed[ied][0];
            faoed[ied_new][1] = faoed[ied][1];
          }
          else {
            assert((nooed[ied][0] == -1)&&(nooed[ied][1] == -1));
          }
        }

        delete[] no_flag;

        nno = nno_new;
        assert(ngr == 1);
        ned = edogr_i[1] = ned_new;

      }
      void packVdVecs(vector<int>& intVec,vector<double>& doubleVec) {
        assert(ngr == 1);

        intVec.push_back(nno);
        intVec.push_back(ned);
        intVec.push_back(nfa);
        intVec.push_back(ngr);
        intVec.push_back(nbf);
        for (int igr = 0; igr < ngr; ++igr)
          intVec.push_back(edogr_i[igr+1]-edogr_i[igr]);
        for (int ied = 0; ied < ned; ++ied) {
          FOR_I2 intVec.push_back(nooed[ied][i]);
          FOR_I2 intVec.push_back(faoed[ied][i]);
        }
        for (int ifa = 0; ifa < nfa; ++ifa) {
          const int icv_nbr = cvofa[ifa];
          assert(icv_nbr != icv);
          int rank,index,bits;
          if (icv_nbr < ncv) {
            rank = mpi_rank;
            bits = 0;
            index = icv_nbr;
          }
          else {
            // this is an off-rank nbr that needs rebuilding...
            assert(icv_nbr-ncv < rbi_g.size());
            BitUtils::unpackRankBitsIndex(rank,bits,index,rbi_g[icv_nbr-ncv]);
          }
          intVec.push_back(BitUtils::packGroupRankBits(0,BitUtils::packRankBits(rank,bits)));
          intVec.push_back(index);
        }
        for (int ibf = 0; ibf < nbf; ++ibf)
          FOR_I2 intVec.push_back(spbobf[ibf][i]);
        for (int ino = 0; ino < nno; ++ino)
          FOR_I3 doubleVec.push_back(x_no[ino][i]);
      }

      /*
         void writeEdge(const int ied) {
         char filename[128];
         sprintf(filename,"edge.%06d.dat",ied);
         FILE * fp = fopen(filename,"w");
         const int np = 20;
         for (int ip = 0; ip <= np; ++ip) {
         const double w = double(ip)/double(np); // 0..1
         double xp[3];
         FOR_I3 xp[i] = w*w*x_no[nooed[ied][0]][i] + (1.0-w*w)*x_no[nooed[ied][1]][i];
         fprintf(fp,"%18.15e %18.15e %18.15e\n",xp[0],xp[1],xp[2]);
         }
         fclose(fp);
         }
         */

  };

  // the goal of stitch is to build the vector of VdReturnData, i.e.
  // the consistent set of voronoi diagrams returned from the
  // load balanced rebuildVd* routines...

  vector<VdReturnData> vdReturnDataVec;
  int * ivrdocv = NULL; // relationship between above and the cv's

  class OrphanChunkData {
    // this ocd is DIFFERENT than the original one in VoronoiBuilder,
    // so be careful...
    public:
      // ints...
      int nno,ned,next;
      int8 rbi;
      int (*nooed)[2];
      // when faoed is positive (or 0), look into this to get the associated part and surface tri
      // when faoed is negative, this is -1 indexing on the cvofa array to get the nbr...
      int8 (*faoed)[2];
      double (*x_no)[3];
      OrphanChunkData() {
        nno = 0;
        ned = 0;
        rbi = -1;
        next = -1;
        nooed = NULL;
        faoed = NULL;
        x_no = NULL;
      }
      // NOTE: DO NOT write a fancy destructor here PLEASE.
      // Leave destructor blank because of allocated arrays. Use a clear()
      // memory management model instead for the vector of OrphanChunkData...
      void clear() {
        nno = 0;
        ned = 0;
        next = -1;
        rbi = -1;
        DELETE(nooed);
        DELETE(faoed);
        DELETE(x_no);
      }
      void writeTecplot(const string& filename) const {
        double x0[3] = { 0.0, 0.0, 0.0 };
        writeTecplot(filename,x0);
      }
      void writeTecplot(const string& filename,const double x0[3]) const {
        FILE * fp = fopen(filename.c_str(),"w");
        fprintf(fp,"TITLE = \"%s\"\n","debug.dat");
        fprintf(fp,"VARIABLES = \"X\"\n");
        fprintf(fp,"\"Y\"\n");
        fprintf(fp,"\"Z\"\n");
        fprintf(fp,"ZONE T=\"%s\"\n",filename.c_str());
        fprintf(fp,"N=%d, E=%d, F=FEPOINT, ET=TRIANGLE\n",nno,ned);
        for (int ino = 0; ino < nno; ++ino) {
          fprintf(fp,"%18.15e %18.15e %18.15e\n",x_no[ino][0]+x0[0],x_no[ino][1]+x0[1],x_no[ino][2]+x0[2]);
        }
        for (int ied = 0; ied < ned; ++ied) {
          fprintf(fp,"%d %d %d\n",nooed[ied][0]+1,nooed[ied][1]+1,nooed[ied][1]+1);
        }
        fclose(fp);
      }
  };

  vector<OrphanChunkData> ocdVec;

  void finalize() {
    if (hcp_x_vd_bbox_adt) {
      delete hcp_x_vd_bbox_adt;
      hcp_x_vd_bbox_adt = NULL;
    }
    if (hcp_x_vd_adt) {
      delete hcp_x_vd_adt;
      hcp_x_vd_adt = NULL;
    }
    for (int ipart = 0; ipart < partVec.size(); ++ipart) {
      if (partVec[ipart]->x_vd_bbox_adt) {
        delete partVec[ipart]->x_vd_bbox_adt;
        partVec[ipart]->x_vd_bbox_adt = NULL;
        assert(partVec[ipart]->x_vd_adt);
        delete partVec[ipart]->x_vd_adt;
        partVec[ipart]->x_vd_adt = NULL;
      }
    }
    // clean up
    for (int ii = 0, limit = vdReturnDataVec.size(); ii < limit; ++ii) {
      vdReturnDataVec[ii].clear();
    }
    vdReturnDataVec.clear();
    if (ivrdocv) {
      delete[] ivrdocv;
      ivrdocv = NULL;
    }
    for (int ii = 0; ii < ocdVec.size(); ++ii) {
      ocdVec[ii].clear();
    }
    ocdVec.clear();
  }

  // used to need this but then skipped this part in setVdReturnData3 because
  // it can be expensive and there is no point in doing it if you are just Lloyd
  // iterating...

  void ensureXvdAdt() {

    // offsets into the part should be available. This should be managed differently eventually
    // because this vector could legitimately be emty on some processor...

    assert(!cvopa_i.empty());

    double my_bbmin[3],my_bbmax[3];
    double (*bbmin)[3] = NULL;
    double (*bbmax)[3] = NULL;
    const int npart = partVec.size();
    for (int ipart = 0; ipart < npart; ++ipart) {
      if (partVec[ipart]->pts) {
        if (partVec[ipart]->x_vd_adt == NULL) {
          assert(partVec[ipart]->x_vd_bbox_adt == NULL); // both should be NULL
          // build the local 3D double adt...
          partVec[ipart]->x_vd_adt = new Adt<double>(cvopa_i[ipart+1]-cvopa_i[ipart],x_vd+cvopa_i[ipart],x_vd+cvopa_i[ipart]);
          // the top leaf of this local adt stores the local bbox...
          partVec[ipart]->x_vd_adt->getBbox(my_bbmin,my_bbmax);
          if (bbmin == NULL) {
            bbmin = new double[mpi_size][3];
            bbmax = new double[mpi_size][3];
          }
          MPI_Allgather(my_bbmin,3,MPI_DOUBLE,bbmin,3,MPI_DOUBLE,mpi_comm);
          MPI_Allgather(my_bbmax,3,MPI_DOUBLE,bbmax,3,MPI_DOUBLE,mpi_comm);
          partVec[ipart]->x_vd_bbox_adt = new Adt<double>(mpi_size,bbmin,bbmax);
        }
      }
    }

    // finally, if there are hcp points AFTER any points that may be associated with the parts,
    // they need to be in a separate adt...
    if (hcp_x_vd_adt == NULL) {
      assert(hcp_x_vd_bbox_adt == NULL);
      assert(cvopa_i.size() == npart+2);
      // build the local 3D double adt...
      hcp_x_vd_adt = new Adt<double>(cvopa_i[npart+1]-cvopa_i[npart],x_vd+cvopa_i[npart],x_vd+cvopa_i[npart]);
      // the top leaf of this local adt stores the local bbox...
      hcp_x_vd_adt->getBbox(my_bbmin,my_bbmax);
      if (bbmin == NULL) {
        bbmin = new double[mpi_size][3];
        bbmax = new double[mpi_size][3];
      }

      // TODO: one all gather? here or for all parts with pts?...
      MPI_Allgather(my_bbmin,3,MPI_DOUBLE,bbmin,3,MPI_DOUBLE,mpi_comm);
      MPI_Allgather(my_bbmax,3,MPI_DOUBLE,bbmax,3,MPI_DOUBLE,mpi_comm);
      hcp_x_vd_bbox_adt = new Adt<double>(mpi_size,bbmin,bbmax);
    }

    // cleanup...
    if (bbmin) {
      delete[] bbmin;
      delete[] bbmax;
    }

  }

  int addNearbySurfaceTrisAndTransforms(vector<pair<int,int> >& storp_v,const double xp[3],const double delta) {

    // this routine adds the surface tris that may cut the cubic voronoi diagram.

    // recall that the seed will be constructed as a centered cube with sides 1.01*delta -- i.e. 1%
    // larger than the maximum size necessary to hold the Voronoi diagram. So surface tris must
    // be included that are within 2% of this box, i.e. 1.02*delta, or +/-0.51*delta...

    const double bbmin[3] = {
      xp[0]-0.51*delta,
      xp[1]-0.51*delta,
      xp[2]-0.51*delta
    };
    const double bbmax[3] = {
      xp[0]+0.51*delta,
      xp[1]+0.51*delta,
      xp[2]+0.51*delta
    };

    vector<int> bboxVec;
    int part_count = 0;
    for (int ipart = 0; ipart < partVec.size(); ++ipart) {
      if (partVec[ipart]->surface) {
        assert(partVec[ipart]->subSurface);
        SubSurface * ss = partVec[ipart]->subSurface;
        // grab the nearby tris from the SubSurface. Recall SubSurface is the part of the
        // surface that may intersect the local points being rebuilt. It should be guaranteed to
        // service this request...
        assert(ss->bbmin[0] <= bbmin[0]);
        assert(ss->bbmin[1] <= bbmin[1]);
        assert(ss->bbmin[2] <= bbmin[2]);
        assert(ss->bbmax[0] >= bbmax[0]);
        assert(ss->bbmax[1] >= bbmax[1]);
        assert(ss->bbmax[2] >= bbmax[2]);
        assert(bboxVec.empty());
        ss->adt->buildListForBBox(bboxVec,bbmin,bbmax);
        // if we got some, add the actual ist's and any associated periodic
        // transforms to the storp_v pair vec. This information can be
        // then sent anywhere to build the starting seed with size +/- 0.505*delta.
        if (!bboxVec.empty()) {
          const size_t size0 = storp_v.size();
          storp_v.resize(size0 + bboxVec.size());
          for (int ii = 0; ii < bboxVec.size(); ++ii) {
            // recall the ss is only a partial copy of the surface on this rank,
            // but its ist's are related to the full shared-memory copy that everyone has...
            const int ist_ss = bboxVec[ii];
            int ist,bits;
            ss->getIstGlobalAndBits(ist,bits,ist_ss);
            storp_v[size0+ii].first = ist;
            assert((bits >= 0)&&(bits < (1<<6))); // Note: allow for 6 bits: have lots of room for more
            storp_v[size0+ii].second = ((ipart<<6)|bits);
          }
          bboxVec.clear();
          ++part_count;
        }
      }
      else {
        assert(partVec[ipart]->subSurface == NULL);
      }
    }

    return part_count;

  }

  void addNbrRanksAndTransforms(vector< pair<int,int> >& rborp_v,const double xp[3],const double delta) {

    // here we are going to build a list of all ranks that we need to send xp (or a possibly transformed xp)
    // to to get the complete set of nbrs within "delta" of xp...

    vector<int> bitVec;
    PeriodicData::buildBitVec(bitVec);

    assert(x_vd);
    vector<int> rankVec;
    set<pair<int,int> > rankSet;
    for (int ipart = 0; ipart < partVec.size(); ++ipart) {
      // not all parts have points...
      if (partVec[ipart]->pts) {
        assert(partVec[ipart]->x_vd_bbox_adt);
        for (int ibv = 0; ibv < bitVec.size(); ++ibv) {
          assert(rankVec.empty());
          const int bits = bitVec[ibv];
          double xp_t[3]; FOR_I3 xp_t[i] = xp[i];
          if (bits != 0) {
	    // try to handle the 180 deg case generally...
	    assert(bits < (1<<6));
	    //if (!testPointValidForTransform(xp,BitUtils::flipPeriodicBits(bits)))
	    if (!testPointValidForTransform(xp,bits))
	      continue;
	    PeriodicData::periodicTranslate(xp_t,1,bits);
	  }
          partVec[ipart]->x_vd_bbox_adt->buildListForSphere(rankVec,xp_t,delta);
          for (int ir = 0; ir < rankVec.size(); ++ir) {
            if ((rankVec[ir] != mpi_rank)||(bits != 0))
              rankSet.insert(pair<int,int>(rankVec[ir],bits));
          }
          rankVec.clear();
        }
      }
    }

    // plus the hcp points...
    assert(hcp_x_vd_bbox_adt);
    for (int ibv = 0; ibv < bitVec.size(); ++ibv) {
      assert(rankVec.empty());
      const int bits = bitVec[ibv];
      double xp_t[3]; FOR_I3 xp_t[i] = xp[i];
      if (bits != 0) {
	// try to handle the 180 deg case generally...
	assert(bits < (1<<6));
	//if (!testPointValidForTransform(xp,BitUtils::flipPeriodicBits(bits)))
	if (!testPointValidForTransform(xp,bits))
	  continue;
	PeriodicData::periodicTranslate(xp_t,1,bits);
      }
      hcp_x_vd_bbox_adt->buildListForSphere(rankVec,xp_t,delta);
      for (int ir = 0; ir < rankVec.size(); ++ir) {
        if ((rankVec[ir] != mpi_rank)||(bits != 0))
          rankSet.insert(pair<int,int>(rankVec[ir],bits));
      }
      rankVec.clear();
    }

    for (set<pair<int,int> >::iterator iter = rankSet.begin(); iter != rankSet.end(); ++iter)
      rborp_v.push_back(*iter);

  }

  void addNearbyNbrs(vector<pair<int,int> >& nborp_v,const int icv,const double xp[3],const double delta,const int bits) {

    assert((icv == -1)||(bits == 0));

    bool found_icv = false;
    vector<int> nbrVec;
    const int npart = partVec.size();
    for (int ipart = 0; ipart < npart; ++ipart) {
      // do all parts have points?...
      if (partVec[ipart]->pts) {
        assert(partVec[ipart]->x_vd_adt);
        assert(nbrVec.empty());
        partVec[ipart]->x_vd_adt->buildListForSphere(nbrVec,xp,delta);
        if (!nbrVec.empty()) {
          for (int ii = 0; ii < nbrVec.size(); ++ii) {
            const int icv_nbr = nbrVec[ii] + cvopa_i[ipart];
            if (icv_nbr == icv) {
              // this is us. We must be active!...
              assert(flag_cv[icv]&ACTIVE_BIT);
              // and we must be being rebuilt...
              assert(flag_cv[icv]&REBUILD_BIT);
              assert(!found_icv);
              found_icv = true;
            }
            else if (flag_cv[icv_nbr]&ACTIVE_BIT) {
              // only add active nbrs...
              nborp_v.push_back( pair<int,int>(icv_nbr,bits) );
            }
          }
          nbrVec.clear();
        }
      }
    }

    assert(hcp_x_vd_adt);
    assert(nbrVec.empty());
    hcp_x_vd_adt->buildListForSphere(nbrVec,xp,delta);
    if (!nbrVec.empty()) {
      for (int ii = 0; ii < nbrVec.size(); ++ii) {
        const int icv_nbr = nbrVec[ii] + cvopa_i[npart];
        if (icv_nbr == icv) {
          // this is us. We must be active!...
          assert(flag_cv[icv]&ACTIVE_BIT);
          // and we must be being rebuilt...
          assert(flag_cv[icv]&REBUILD_BIT);
          assert(!found_icv);
          found_icv = true;
        }
        else if (flag_cv[icv_nbr]&ACTIVE_BIT) {
          // only add active nbrs...
          nborp_v.push_back( pair<int,int>(icv_nbr,bits) );
        }
      }
      nbrVec.clear();
    }

    // make sure we found ourselves...
    if (!(found_icv || (icv == -1))) {
      cout << "rank: " << mpi_rank << " FAILED to find icv: " << icv << endl;
    }
    assert(found_icv || (icv == -1));

  }

  void addNearbyNbrs(vector<pair<int,int> >& nborc_v,const double xp[3],const double delta,const int bits) {

    // this one for nbrs with respect to a point, but no ip known...

    addNearbyNbrs(nborc_v,-1,xp,delta,bits);

  }

  class VdTri {
    public:
      double* x[3];
      VdTri(double * x0,double * x1,double * x2) {
        x[0] = x0;
        x[1] = x1;
        x[2] = x2;
      }
  };

  void setVdReturnData3(VdReturnData& vd,const CuttableVoronoiData& cvd,const double x_vd[3],const vector<pair<double,pair<int,int> > >& nbrVec,const int * const ist_and_ipart_bits) {

    // This is a lighter version than in setVdReturnData2 (in moving), and the whole VdReturnData class was changed

    int * no_flag = new int[cvd.nno];
    for (int ino = 0; ino < cvd.nno; ++ino)
      no_flag[ino] = ino;

    for (int ied = 0; ied < cvd.ned; ++ied) {
      int ino0 = no_flag[cvd.nooed[ied][0]];
      while (ino0 != no_flag[ino0])
        ino0 = no_flag[ino0];
      int ino1 = no_flag[cvd.nooed[ied][1]];
      while (ino1 != no_flag[ino1])
        ino1 = no_flag[ino1];
      no_flag[ino0] = no_flag[ino1] = min(ino0,ino1);
    }

    // now count the groups...
    int ngr = 0;
    for (int ino = 0; ino < cvd.nno; ++ino) {
      if (ino == no_flag[ino]) {
        ++ngr;
        no_flag[ino] = -ngr; // use -ve indexing
      }
      else {
        int ino_ = no_flag[ino];
        while (ino_ >= 0)
          ino_ = no_flag[ino_];
        no_flag[ino] = ino_;
      }
    }

    int * gr_flag = NULL;

    // for the case of multiple groups, we may need to combine groups
    // together beause of "donut faces"...
    if (ngr > 1) {

      // first, figure out if any of the faces are associated with "donut" faces
      // by looking at the total internal face volume assocated with each
      // group. A negative volume indicates a loop that should be associated with
      // another group...

      double * vol_gr = new double[ngr];
      for (int igr = 0; igr < ngr; ++igr)
        vol_gr[igr] = 0.0;

      map<const pair<int,int>,int> faceMap;
      for (int ied = 0; ied < cvd.ned; ++ied) {
        const int igr = -no_flag[cvd.nooed[ied][0]]-1;
        assert(igr == -no_flag[cvd.nooed[ied][1]]-1);
        assert((igr >= 0)&&(igr < ngr));
        // only necessary to consider internal faces in this volume check...
        if (cvd.faoed[ied][0] < 0) {
          // negative is an internal face. Index should be offset
          // by 8 for the original seed surfaces, which should be
          // cut away at this point...
          assert(cvd.faoed[ied][0] <= -8);
          map<const pair<int,int>,int>::iterator iter = faceMap.find(pair<int,int>(igr,cvd.faoed[ied][0]));
          if (iter == faceMap.end()) {
            faceMap[pair<int,int>(igr,cvd.faoed[ied][0])] = cvd.nooed[ied][0];
          }
          else if ((cvd.nooed[ied][0] != iter->second)&&(cvd.nooed[ied][1] != iter->second)) {
            vol_gr[igr] += SIGNED_TET_VOLUME_6_VEC(cvd.x_no[iter->second],cvd.x_no[cvd.nooed[ied][0]],cvd.x_no[cvd.nooed[ied][1]]);
          }
        }
        if (cvd.faoed[ied][1] < 0) {
          // negative is an internal face. Index should be offset
          // by 8 for the original seed surfaces, which should be
          // cut away at this point...
          assert(cvd.faoed[ied][1] <= -8);
          map<const pair<int,int>,int>::iterator iter = faceMap.find(pair<int,int>(igr,cvd.faoed[ied][1]));
          if (iter == faceMap.end()) {
            faceMap[pair<int,int>(igr,cvd.faoed[ied][1])] = cvd.nooed[ied][1];
          }
          else if ((cvd.nooed[ied][0] != iter->second)&&(cvd.nooed[ied][1] != iter->second)) {
            vol_gr[igr] += SIGNED_TET_VOLUME_6_VEC(cvd.x_no[iter->second],cvd.x_no[cvd.nooed[ied][1]],cvd.x_no[cvd.nooed[ied][0]]);
          }
        }
      }

      // negative group volumes (likely) indicate a face penetration by some
      // non-convex geometry (donut), and should be grouped with the main volume. When
      // there are both negative volumes AND orphans, we have more work to do...

      const int ngr_old = ngr;
      assert(gr_flag == NULL);
      gr_flag = new int[ngr_old];
      ngr = 0;
      for (int igr = 0; igr < ngr_old; ++igr) {
        if (vol_gr[igr] < 0.0) {
          gr_flag[igr] = -1; // use -1 indexing for the donut groups
        }
        else {
          gr_flag[igr] = ngr++;
        }
      }

      //cout << "ngr_old: " << ngr_old << ", ngr: " << ngr << endl;
      //for (int igr_old = 0; igr_old < ngr_old; ++igr_old)
      //cout << "igr_old: " << igr_old << ", gr_flag[igr_old]: " << gr_flag[igr_old] << " vol_gr[igr_old]: " << vol_gr[igr_old] << endl;
      //getchar();

      // now set the no_flag[ino] to the positive group index (it is currently -1 indexed
      // to the ngr_old groups)...

      if (ngr < ngr_old) {

        // for the case of one main group, this is trivial, and we just combine all
        // groups into the one main group, group 0...

        if (ngr == 1) {

          // one main group, so everyone should be part of this...
          for (int ino = 0; ino < cvd.nno; ++ino)
            no_flag[ino] = 0;

        }
        else {

          assert(ngr > 1);

          // if you get here, your mesh has both "donut" faces and legitimate orphans, and
          // we will need to figure out which donut goes with which group. This is not that hard:
          // we need to just choose a point on the inside face (the one associated with the
          // negative volume) and see if when we sum the angles subtended by the outside face
          // edges, we get 360, or zero. We should only get 360 for one of the faces, and that
          // is the group we should associate with...
          //
          // It is also possible that we have small negative regions that are completely
          // separated from the main volume -- i.e. floating particles. These will be
          // a little more difficult to handle

          // now loop on edges. When we hit one in one of these donut groups, record the
          // internal face as a negative index in the gr_flag...

          for (int ied = 0; ied < cvd.ned; ++ied) {
            const int igr = -no_flag[cvd.nooed[ied][0]]-1;
            assert(igr == -no_flag[cvd.nooed[ied][1]]-1);
            if (gr_flag[igr] < 0) {
              int ifa_check = -1;
              if (cvd.faoed[ied][0] < 0) {
                assert(cvd.faoed[ied][0] < -1);
                ifa_check = cvd.faoed[ied][0];
              }
              if (cvd.faoed[ied][1] < 0) {
                assert(ifa_check == -1);
                assert(cvd.faoed[ied][1] < -1);
                ifa_check = cvd.faoed[ied][1];
              }
              if (ifa_check != -1) {
                assert(ifa_check <= -8);
                if (gr_flag[igr] == -1) {
                  gr_flag[igr] = ifa_check; // set the group flag to the -8 indexed nbr face
                }
                // does this matter? I think you can just use the first face to decide on the main group it is a part of
                //else if (gr_flag[igr] != ifa_check) {
                //cout << gr_flag[igr] << " " << ifa_check << endl;
                //cerr << "Negative group seems to have more than one face -- could be a tube through another cv?" << endl;
                //for (int ino = 0; ino < cvd.nno; ++ino) if (-no_flag[ino]-1 == igr)
                //  cout << cvd.x_no[ino][0] << "," << cvd.x_no[ino][1] << "," << cvd.x_no[ino][2] << endl;
                //cvd.writeBinary("cvd_debug_tube.bin");
                //assert(0);
                //}
              }
            }
          }

          // make sure all gr_flags were set with the internal face they are part of...

          for (int igr_donut = 0; igr_donut < ngr_old; ++igr_donut) {
            assert((gr_flag[igr_donut] >= 0)||(gr_flag[igr_donut] <= -8)); // if you hit this, maybe a floating/subgrid piece of geometry?
            if (gr_flag[igr_donut] <= -8) {
              //cout << "Working on donut group " << igr_donut << endl;
              // -------------------------------------------------------------------
              // first check is if there is a single main that has the face of the donut
              // now find edges that share this face in a main group...
              // -------------------------------------------------------------------
              int igr_main_match = -1;
              for (int igr_main = 0; igr_main < ngr_old; ++igr_main) {
                assert((gr_flag[igr_main] >= 0)||(gr_flag[igr_main] < -1));
                if (gr_flag[igr_main] >= 0) {
                  //cout << "Got main group: " << gr_flag[igr_main] << endl;
                  for (int ied = 0; ied < cvd.ned; ++ied) {
                    const int this_igr = -no_flag[cvd.nooed[ied][0]]-1;
                    assert(this_igr == -no_flag[cvd.nooed[ied][1]]-1);
                    if (this_igr == igr_main) {
                      FOR_I2 {
                        if (cvd.faoed[ied][i] == gr_flag[igr_donut]) {
                          //double x_ed[3];
                          //FOR_I3 x_ed[i] = 0.5*(cvd.x_no[cvd.nooed[ied][0]][i]+cvd.x_no[cvd.nooed[ied][1]][i]);
                          //cout << " > found edge touching same face " << COUT_VEC(x_ed) << endl;
                          if (igr_main_match == -1) {
                            igr_main_match = igr_main;
                          }
                          else if (igr_main_match != igr_main) {
                            igr_main_match = -2;
                            break;
                          }
                        }
                      }
                    }
                  }
                  if (igr_main_match == -2)
                    break;
                }
              }
              if (igr_main_match == -1) {
                cerr << "ERROR: did not find any main with this donuts face." << endl;
                cvd.writeBinary("cvd_debug_no_match.bin");
                assert(0);
              }
              else if (igr_main_match >= 0) {
                // found one and only one match for this donut. Set its gr_flag to the
                // new group associated with this match...
                assert((gr_flag[igr_main_match] >= 0)&&(gr_flag[igr_main_match] < ngr));
                assert(gr_flag[igr_donut] <= -8);
                gr_flag[igr_donut] = -gr_flag[igr_main_match]-1; // -1 based new group index
              }
              else {
                assert(igr_main_match == -2);
                // find a node on this group that touches the face...
                int ino_center = -1;
                for (int ied = 0; ied < cvd.ned; ++ied) {
                  const int this_igr = -no_flag[cvd.nooed[ied][0]]-1;
                  assert(this_igr == -no_flag[cvd.nooed[ied][1]]-1);
                  if (this_igr == igr_donut) {
                    if ((cvd.faoed[ied][0] == gr_flag[igr_donut])||(cvd.faoed[ied][1] == gr_flag[igr_donut])) {
                      ino_center = cvd.nooed[ied][0];
                      break;
                    }
                  }
                }
                assert(ino_center >= 0);
                // angle sum stuff...
                vector<pair<double,int> > sweepVec(ngr_old);
                for (int igr_main = 0; igr_main < ngr_old; ++igr_main) {
                  sweepVec[igr_main] = pair<double,int>(0.0,igr_main);
                  assert((gr_flag[igr_main] >= 0)||(gr_flag[igr_main] < -1));
                  if (gr_flag[igr_main] >= 0) {
                    //cout << "Got main group: " << gr_flag[igr_main] << endl;
                    for (int ied = 0; ied < cvd.ned; ++ied) {
                      const int this_igr = -no_flag[cvd.nooed[ied][0]]-1;
                      assert(this_igr == -no_flag[cvd.nooed[ied][1]]-1);
                      if (this_igr == igr_main) {
                        FOR_I2 {
                          if (cvd.faoed[ied][i] == gr_flag[igr_donut]) {
                            const double dx0[3] = DIFF(cvd.x_no[cvd.nooed[ied][0]],cvd.x_no[ino_center]);
                            const double dx1[3] = DIFF(cvd.x_no[cvd.nooed[ied][1]],cvd.x_no[ino_center]);
                            sweepVec[igr_main].first += acos(DOT_PRODUCT(dx0,dx1)/(MAG(dx0)*MAG(dx1)));
                            //double x_ed[3];
                            //FOR_I3 x_ed[i] = 0.5*(cvd.x_no[cvd.nooed[ied][0]][i]+cvd.x_no[cvd.nooed[ied][1]][i]);
                            //cout << " > found edge touching same face " << COUT_VEC(x_ed) << endl;
                          }
                        }
                      }
                    }
                  }
                  sweepVec[igr_main].first = fabs(sweepVec[igr_main].first);
                }
                sort(sweepVec.begin(),sweepVec.end());
                // take the largest abs sweep angle as the group (should be very close to 360 deg)
                gr_flag[igr_donut] = -gr_flag[sweepVec[ngr_old-1].second]-1; // -1 based new group index
                //for (int ii = 0; ii < ngr_old; ++ii)
                //  cout << sweepVec[ii].second << " " << sweepVec[ii].first*180.0/M_PI << endl;

                //cerr << "Donut face occurs on more than one main cv. Need to implement swept angle test." << endl;
                //cvd.writeBinary("cvd_debug_swept_angle.bin");
                //assert(0);
              }
            }
          } // for (int igr_donut...

          // check that donut groups were all assigned...

          for (int igr_old = 0; igr_old < ngr_old; ++igr_old) {
            const int igr_new = gr_flag[igr_old];
            if (igr_new >= 0) {
              if (!((igr_new >= 0)&&(igr_new < ngr)))
                cvd.writeBinary("cvd_debug_wtf.bin");
              assert((igr_new >= 0)&&(igr_new < ngr));
            }
            else {
              assert((-igr_new-1 >= 0)&&(-igr_new-1 < ngr));
            }
          }

          // and modify node to reflect new main group number or, in the
          // case of donut groups, to joint the main grouop they were determined to
          // be part of...

          for (int ino = 0; ino < cvd.nno; ++ino) {
            const int igr_old = -no_flag[ino]-1;
            assert((igr_old >= 0)&&(igr_old < ngr_old));
            const int igr_new = gr_flag[igr_old];
            if (igr_new >= 0) {
              assert((igr_new >= 0)&&(igr_new < ngr));
              no_flag[ino] = igr_new;
            }
            else {
              assert((-igr_new-1 >= 0)&&(-igr_new-1 < ngr));
              no_flag[ino] = -igr_new-1;
            }
          }

        }

      }
      else {

        assert(ngr == ngr_old);

        // there was no group combination, so just flip no_flag to reflect
        // a positive group index 0,1,2...

        for (int ino = 0; ino < cvd.nno; ++ino)
          no_flag[ino] = -no_flag[ino]-1;

      }

      delete[] vol_gr;

      // if the number of groups is still larger than 1, we need to reorder groups
      // so the first group is the primary group. We determine the primary group
      // by computing proximity to the boundary tris. The boundary we are closest to
      // will be in the closest group...

      if (ngr > 1) {

        // use gr_flag to ensure that all groups have boundary tris. If they don't
        // we have made a mistake with the donut merge above?...
        for (int igr = 0; igr < ngr; ++igr)
          gr_flag[igr] = 0;

        int igr_closest = -1;
        double d2_closest;
        map<const int,int> firstNodeMap;
        for (int ied = 0; ied < cvd.ned; ++ied) {
          if (cvd.faoed[ied][0] >= 0) {
            // a positive index in faoed at this point refers to an index of a cut
            // boundary face. It is a local indexing based on recv_bit_buf at this point.
            map<const int,int>::iterator iter = firstNodeMap.find(cvd.faoed[ied][0]);
            if (iter == firstNodeMap.end()) {
              firstNodeMap[cvd.faoed[ied][0]] = cvd.nooed[ied][0];
              const int igr = no_flag[cvd.nooed[ied][0]];
              assert((igr >= 0)&&(igr < ngr));
              ++gr_flag[igr];
            }
            else if ((iter->second != cvd.nooed[ied][0])&&(iter->second != cvd.nooed[ied][1])) {
              // here we could have a simplified function that assumes x0 is (0,0,0)...
              const double x0[3] = { 0.0, 0.0, 0.0 };
              const double this_d2 = getPointToTriDist2(x0,
                  cvd.x_no[iter->second],
                  cvd.x_no[cvd.nooed[ied][0]],
                  cvd.x_no[cvd.nooed[ied][1]]);
              if ((igr_closest == -1)||(this_d2 < d2_closest)) {
                igr_closest = no_flag[iter->second];
                d2_closest = this_d2;
                assert((igr_closest >= 0)&&(igr_closest < ngr));
                assert(no_flag[cvd.nooed[ied][0]] == igr_closest);
                assert(no_flag[cvd.nooed[ied][1]] == igr_closest);
              }
            }
          }
          if (cvd.faoed[ied][1] >= 0) {
            map<const int,int>::iterator iter = firstNodeMap.find(cvd.faoed[ied][1]);
            if (iter == firstNodeMap.end()) {
              firstNodeMap[cvd.faoed[ied][1]] = cvd.nooed[ied][1];
              const int igr = no_flag[cvd.nooed[ied][1]];
              assert((igr >= 0)&&(igr < ngr));
              ++gr_flag[igr];
            }
            else if ((iter->second != cvd.nooed[ied][0])&&(iter->second != cvd.nooed[ied][1])) {
              // here we could have a simplified function that assumes x0 is (0,0,0)...
              const double x0[3] = { 0.0, 0.0, 0.0 };
              const double this_d2 = getPointToTriDist2(x0,
                  cvd.x_no[iter->second],
                  cvd.x_no[cvd.nooed[ied][1]],
                  cvd.x_no[cvd.nooed[ied][0]]);
              if ((igr_closest == -1)||(this_d2 < d2_closest)) {
                igr_closest = no_flag[iter->second];
                d2_closest = this_d2;
                assert((igr_closest >= 0)&&(igr_closest < ngr));
                assert(no_flag[cvd.nooed[ied][0]] == igr_closest);
                assert(no_flag[cvd.nooed[ied][1]] == igr_closest);
              }
            }
          }
        }

        assert(igr_closest >= 0);
        int ngr_other = 1;
        for (int igr = 0; igr < ngr; ++igr) {
          assert(gr_flag[igr] > 0);
          if (igr == igr_closest)
            gr_flag[igr] = 0;
          else
            gr_flag[igr] = ngr_other++;
        }
        assert(ngr_other == ngr);

        for (int ino = 0; ino < cvd.nno; ++ino) {
          const int igr_old = no_flag[ino];
          assert((igr_old >= 0)&&(igr_old < ngr));
          no_flag[ino] = gr_flag[igr_old];
        }

        // check that no_flag is set and all groups are represented
        for (int igr = 0; igr < ngr; ++igr)
          gr_flag[igr] = 0;

        for (int ino = 0; ino < cvd.nno; ++ino) {
          const int igr = no_flag[ino];
          assert((igr >= 0)&&(igr < ngr));
          ++gr_flag[igr];
        }

        int start = 0;
        for (int igr = 0; igr < ngr; ++igr) {
          assert(gr_flag[igr] > 0);
          const int inc = gr_flag[igr];
          gr_flag[igr] = start;
          start += inc;
        }

      }
      else {

        assert(gr_flag);
        gr_flag[0] = 0;

      }

    }
    else {

      assert(ngr == 1);
      assert(gr_flag == NULL);
      gr_flag = new int[1];
      gr_flag[0] = 0;

      for (int ino = 0; ino < cvd.nno; ++ino)
        no_flag[ino] = 0;
    }

    // build edogr_i and reorder edges, nodes...

    // we also need stogr_i/v and the face structures...

    // and nodes/edges, but no x_fa, x_cv, vol, ...

    vd.ngr = ngr;
    assert(vd.edogr_i == NULL);
    vd.edogr_i = new int[ngr+1];
    for (int igr = 0; igr < ngr; ++igr)
      vd.edogr_i[igr+1] = 0;
    for (int ied = 0; ied < cvd.ned; ++ied) {
      const int igr = no_flag[cvd.nooed[ied][0]];
      assert(igr == no_flag[cvd.nooed[ied][1]]);
      assert((igr >= 0)&&(igr < ngr));
      ++vd.edogr_i [igr+1];
    }
    // make vd.edogr_i csr...
    vd.edogr_i[0] = 0;
    for (int igr = 0; igr < ngr; ++igr)
      vd.edogr_i[igr+1] += vd.edogr_i[igr];
    assert(vd.edogr_i[ngr] == cvd.ned);

    // now reorder edges into vd.nooed,vd.faoed and count
    // nfa,nbf...
    map<const int,pair<int,double> > fa_r2_map;
    map<const pair<int,int>,int> bf_spb_map;
    assert(vd.nooed == NULL);
    assert(vd.faoed == NULL);
    vd.nfa = 0;
    vd.nbf = 0;
    vd.ned = cvd.ned;
    vd.nooed = new int[cvd.ned][2];
    vd.faoed = new int[cvd.ned][2];
    for (int ied = 0; ied < cvd.ned; ++ied) {
      const int igr = no_flag[cvd.nooed[ied][0]];
      assert(igr == no_flag[cvd.nooed[ied][1]]);
      assert((igr >= 0)&&(igr < ngr));
      const int ied_vd = vd.edogr_i[igr]++;
      FOR_I2 vd.nooed[ied_vd][i] = cvd.nooed[ied][i];
      FOR_I2 {
        const int ifa = cvd.faoed[ied][i];
        if (ifa < 0) {
          const double r2_ed = max(DOT_PRODUCT(cvd.x_no[cvd.nooed[ied][0]],cvd.x_no[cvd.nooed[ied][0]]),
              DOT_PRODUCT(cvd.x_no[cvd.nooed[ied][1]],cvd.x_no[cvd.nooed[ied][1]]));
          const int inbr = -ifa-8; assert((inbr >= 0)&&(inbr < nbrVec.size()));
          map<const int,pair<int,double> >::iterator it = fa_r2_map.find(inbr);
          if (it == fa_r2_map.end()) {
            fa_r2_map[inbr] = pair<int,double>(vd.nfa,r2_ed);
            ++vd.nfa;
            vd.faoed[ied_vd][i] = -vd.nfa; // simple -1 indexing into fa_r2_map
          }
          else {
            vd.faoed[ied_vd][i] = -it->second.first-1;
            it->second.second = max(it->second.second,r2_ed);
          }
        }
        else {
          // this face refers to a boundary surface tri...
          const int ist = ist_and_ipart_bits[ifa];
          const int ipart_and_bits = ist_and_ipart_bits[ifa+1];
          //const int bits = (ist_and_ipart_bits[ifa+1]&MASK_6BITS);
          map<const pair<int,int>,int>::iterator it = bf_spb_map.find(pair<int,int>(ist,ipart_and_bits));
          if (it == bf_spb_map.end()) {
            bf_spb_map[pair<int,int>(ist,ipart_and_bits)] = vd.faoed[ied_vd][i] = vd.nbf++;
          }
          else {
            vd.faoed[ied_vd][i] = it->second;
          }
        }
      }
    }

    // reset csr...
    for (int igr = ngr-1; igr > 0; --igr)
      vd.edogr_i[igr] = vd.edogr_i[igr-1];
    vd.edogr_i[0] = 0;

    // use the maps to build the spbobf and cvofa...
    assert(vd.cvofa == NULL);
    assert(fa_r2_map.size() == vd.nfa);
    vd.cvofa = new int[vd.nfa];
    vd.r2_fa = new double[vd.nfa];
    for (map<const int,pair<int,double> >::iterator it = fa_r2_map.begin(); it != fa_r2_map.end(); ++it) {
      const int ifa = it->second.first;
      assert((ifa >= 0)&&(ifa < vd.nfa));
      const int inbr = it->first; assert((inbr >= 0)&&(inbr < nbrVec.size()));
      vd.cvofa[ifa] = nbrVec[inbr].second.first;
      vd.r2_fa [ifa] = it->second.second;
    }

    assert(vd.spbobf == NULL);
    vd.spbobf = new int[vd.nbf][2];
    for (map<const pair<int,int>,int>::iterator it = bf_spb_map.begin(); it != bf_spb_map.end(); ++it) {
      const int ibf = it->second;
      vd.spbobf[ibf][0] = it->first.first;
      vd.spbobf[ibf][1] = it->first.second;
    }

    // finally node-based quantites...
    assert(gr_flag);
    assert(gr_flag[0] == 0);
    assert(vd.x_no == NULL);
    vd.nno = cvd.nno;
    vd.x_no = new double[cvd.nno][3];
    for (int ino = 0; ino < cvd.nno; ++ino) {
      const int igr = no_flag[ino];
      const int ino_vd = gr_flag[igr]++;
      no_flag[ino] = ino_vd;
      FOR_I3 vd.x_no[ino_vd][i] = cvd.x_no[ino][i];
    }

    // finally, update the node indices...
    for (int ied = 0; ied < cvd.ned; ++ied)
      FOR_I2 vd.nooed[ied][i] = no_flag[vd.nooed[ied][i]];

    delete[] gr_flag;
    delete[] no_flag;

  }

  void _rebuildFlaggedVds(vector<VdReturnData>& vdReturnDataVec,const bool b_surface,
      int * send_count,int * send_disp,int * recv_count,int * recv_disp) {

    double my_time_buf[2] = { 0.0, 0.0 };

    //auto t0 = std::chrono::high_resolution_clock::now();
    //double my_time_buf_sum[2] = { 0.0, 0.0 };
    //double my_time_buf_min[2] = { HUGE_VAL, HUGE_VAL };
#ifdef WITH_CHRONO2
    vector<VdTimingData> timingVec;
    static int timing_index = 0;
    ++timing_index;
#endif

    static int debug_count = 0;
    //const int debug_count_target = -1; //64488; //-1

    // we will need a current x_vd_adt always. We may or may not need a
    // current surface adt. If all rebuilding is occuring internally, then
    // we don't. See how bool b_surface is used below...

    // Note: if this gets to be dominant in terms of cost, we could
    // do it on a part basis, and even do it once after load balance if
    // only SBM is considered. But for now, this is the most general...

#ifdef WITH_CHRONO
    auto ta = std::chrono::high_resolution_clock::now();
#endif

    ensureXvdAdt();

#ifdef WITH_CHRONO
    auto tb = std::chrono::high_resolution_clock::now();
#endif

    // define a cost where everyone below will not move...

    int8 * rcorp = NULL; // rebuild-cost of rebuild-point

    // these kept out of loop to avoid reallocation...

    vector<int> cvorp;
    int8 * rcora = NULL;
    int * wrorp = NULL;

    int * storp_i = NULL;
    vector<pair<int,int> > storp_v;

    int * rborp_i = NULL;
    vector<pair<int,int> > rborp_v; // rank and bits (0 for now) as a pair...

    int send_count_int_sum_max = 0;
    int * send_buf_int = NULL;
    int send_count_double_sum_max = 0;
    double * send_buf_double = NULL;

    int recv_count_int_sum_max = 0;
    int * recv_buf_int = NULL;
    int recv_count_double_sum_max = 0;
    double * recv_buf_double = NULL;

    int * send_count_double = new int[mpi_size];
    int * send_disp_double  = new int[mpi_size];

    int * nborc_i = NULL;
    int nrecv_p1_max = 0;
    vector<pair<int,int> > nborc_v; // nbr and bits of recv-block...

    int * nborp_i = NULL; // nbr of rebuild points
    vector<pair<int,int> > nborp_v; // nbr (and bits) of rebuild point...

    vector<pair<int,int> > vdVec;
    multimap<pair<int,int>,pair<int,int> > nbrMultimap;

    CuttableVoronoiData cvd;
    vector<int> bbVec;
    map<const uint8,int> nodeMap;
    map<const pair<int,int>,int> edgeMap;
    map<const pair<uint8,uint8>,int> edgeMap2;
    vector<pair<double,pair<int,int> > > nbrVec; // double is the dist2, ints is nbr index
    uint8 pbsost[3];
    int edost[3];
    int noost[3];

    vector<int> intVec;
    vector<double> doubleVec;
    vector<VdReturnData> vrdVec;

    // ======================================================================
    // now iterate because the delta_vd may not be large enough to
    // build the voronoi diagrams correctly...
    // ======================================================================

    double bbmin[3],bbmax[3];
    double bbmin2[3],bbmax2[3];

#ifdef WITH_CHRONO
    auto tc = std::chrono::high_resolution_clock::now();
#endif

    int iter = 0;
    bool done = false;
    while (!done) {

      ++iter;
      if (mpi_rank == 0) cout << " > iter: " << iter << endl;

#ifdef WITH_CHRONO2
      auto t1 = std::chrono::high_resolution_clock::now();
#endif

      // store the initial size of vdReturnDataVec for this iteration...

      const int vdReturnDataVec_size0 = vdReturnDataVec.size();

      // figure out how many vd's each rank has to build. At the same time,
      // ensure our fast local copy of the surfaces (the subsurfaces) are
      // large enough...

      // TODO - here we could use part-based bounding boxes as well, and build the subsurface
      // in an additive way?...

      FOR_I3 bbmin[i] = HUGE_VAL;
      FOR_I3 bbmax[i] = -HUGE_VAL;
      FOR_I3 bbmin2[i] = HUGE_VAL;
      FOR_I3 bbmax2[i] = -HUGE_VAL;
      assert(cvorp.empty());
      FOR_ICV {
        if (flag_cv[icv]&REBUILD_BIT) {
          assert(flag_cv[icv]&ACTIVE_BIT);
          cvorp.push_back(icv);
          // figure out the point bbox at 2 sizes. The first size is the bare minimum
          // required for this iteration. The second is one that is larger and will
          // potentially accomodate some growth in delta (or motion for that matter,
          // if the surface is not moving)...
          FOR_I3 bbmin[i] = min(bbmin[i],x_vd[icv][i]-0.51*delta_vd[icv]);
          FOR_I3 bbmax[i] = max(bbmax[i],x_vd[icv][i]+0.51*delta_vd[icv]);
          FOR_I3 bbmin2[i] = min(bbmin2[i],x_vd[icv][i]-0.51*delta_vd[icv]*5.0); // was: 1.5*1.5 = 2.25+eps -- i.e. 2 growths in delta
          FOR_I3 bbmax2[i] = max(bbmax2[i],x_vd[icv][i]+0.51*delta_vd[icv]*5.0); // increase to 5.0 to avoid rebuild
        }
      }
      int nrp = cvorp.size(); // number of rebuild points

      // make sure the local subSurfaces are large enough for the local points being rebuilt...

      if ((nrp > 0)&&(b_surface)) {
        for (int ipart = 0; ipart < partVec.size(); ++ipart) {
          if (partVec[ipart]->surface) {
            if (partVec[ipart]->subSurface == NULL) {
              partVec[ipart]->buildOrRebuildSubSurface(bbmin2,bbmax2);
            }
            else if ((bbmin[0] < partVec[ipart]->subSurface->bbmin[0])||
                (bbmin[1] < partVec[ipart]->subSurface->bbmin[1])||
                (bbmin[2] < partVec[ipart]->subSurface->bbmin[2])||
                (bbmax[0] > partVec[ipart]->subSurface->bbmax[0])||
                (bbmax[1] > partVec[ipart]->subSurface->bbmax[1])||
                (bbmax[2] > partVec[ipart]->subSurface->bbmax[2])) {
              partVec[ipart]->buildOrRebuildSubSurface(bbmin2,bbmax2);
            }
          }
        }
      }

      // count, pack, send first set of messages...
      // this first set is a message to other processors that may have to send
      // the nbrs they own to work-rank...

      if (storp_i == NULL) storp_i = new int[nrp+1]; // surface-tri of rebuild points
      assert(storp_v.empty());
      storp_i[0] = 0;
      if (rcorp == NULL) rcorp = new int8[nrp]; // rebuild-cost of rebuild points
      int8 rc = 0;
      if (rborp_i == NULL) rborp_i = new int[nrp+1]; // rank-bits of rebuild points
      assert(rborp_v.empty());
      rborp_i[0] = 0;
      const int nrp0 = nrp;
      nrp = 0;
      for (int irp = 0; irp < nrp0; ++irp) {
        const int icv = cvorp[irp];
        int part_count = 0;
        if (b_surface) part_count = addNearbySurfaceTrisAndTransforms(storp_v,x_vd[icv],delta_vd[icv]);
        storp_i[irp+1] = storp_v.size();
        // COST MODEL:
        // the "rebuild cost" depends on both the number of tris grabbed AND the part_count, because
        // this will cause the intersection part of the cutting to be performed...
        rcorp[irp] = 10 + part_count*(storp_i[irp+1]-storp_i[irp]); // TODO: could use the previous "cost" here if the diagram has been built before
        rc += rcorp[irp];
        // and the non-local nbr requests within the radius=delta-sphere...
        addNbrRanksAndTransforms(rborp_v,x_vd[icv],delta_vd[icv]);
        rborp_i[irp+1] = rborp_v.size();
        ++nrp;
        if (rc >= rc_limit) 
          break;
      }
      MiscUtils::dumpRange(&rc,1,"RC");

#ifdef WITH_CHRONO2
      auto t2 = std::chrono::high_resolution_clock::now();
#endif

      // now based on cost decide where everyone goes...

      MiscUtils::buildXora(rcora,rc);
      assert(rcora[mpi_rank+1]-rcora[mpi_rank] == rc);

      // allocate on the first time only. nrp will only get smaller on subsequent iterations...
      if (wrorp == NULL) wrorp = new int[nrp0]; // work-rank of rebuild-point. use un-limited nrp to ensure it gets smaller on each iter

      // HACK: use a simple cost model that considers integer message size...
      /*
      int8 my_cost = (4+2*30)*nrp + 2*storp_i[nrp]; // assume 30 nbrs
      int8 cost_sum;
      MPI_Scan(&my_cost,&cost_sum,1,MPI_INT8,MPI_SUM,mpi_comm);
      int8 cost_total = cost_sum;
      MPI_Bcast(&cost_total,1,MPI_INT8,mpi_size-1,mpi_comm);
      cost_sum -= my_cost;
      for (int irp = 0; irp < nrp; ++irp) {
        wrorp[irp] = cost_sum*mpi_size/cost_total; assert((wrorp[irp] >= 0)&&(wrorp[irp] < mpi_size));
        cost_sum += 4+2*30 + 2*(storp_i[irp+1]-storp_i[irp]);
      }
      */

      // sum all of the cost before us that gets redistributed, i.e. not the
      // cost kept local...
      int8 rc_sum = 0;
      for (int rank = 0; rank < mpi_rank; ++rank)
        rc_sum += max(int8(0),rcora[rank+1]-rcora[rank]-rc_local_threshold);
      const int8 rc_target = rcora[mpi_size]/mpi_size+1;
      int8 rc_local_sum = 0;
      for (int irp = 0; irp < nrp; ++irp) {
        if (rc_local_sum < rc_local_threshold) {
          // this guy gets rebuilt locally
          wrorp[irp] = mpi_rank;
          rc_local_sum += rcorp[irp];
        }
        else {
          // this is the part of the cost that gets redistributed...
          // part of this logic can get pushed out of the loop, but move on for now...
          int8 rc_sum_check = rc_sum; // full redistributed sum
          int rank;
          for (rank = 0; rank < mpi_size-1; ++rank) {
            // the amount of extra room on this rank was...
            rc_sum_check -= rc_target - min(rcora[rank+1]-rcora[rank],rc_local_threshold);
            if (rc_sum_check < 0)
              break;
          }
          assert(rank < mpi_size);
          wrorp[irp] = rank;
          rc_sum += rcorp[irp];
        }
      }

      // =============================================================================
      // we now know where each of our rp's is going to be built. We need to send it
      // there along with any nbrs we own, and also send a message to other ranks who
      // may have other nbrs to send their nbrs to rp's rank...
      // =============================================================================

      // at this point, storp_i/v contains the ist and ipart/bits for each rp, and
      // rborp_i/v contains the ranks (besides our own) for each rp that may have nbrs...

      FOR_RANK send_count[rank] = 0;

      for (int irp = 0; irp < nrp; ++irp) {
        // the wrank, or work-rank is easily calculated, but we store in an array just in case
        // it is handled more complexly in the future -- e.g. better balancing communication and work, etc...
        //const int wrank = wrorp[irp];
        //assert((wrank >= 0)&&(wrank < mpi_size));
        // we also want to check the rank and bits where we are heading. If this
        // is our rank, then no message is required.
        for (int ror = rborp_i[irp]; ror != rborp_i[irp+1]; ++ror) {
          // we need to create a message for every request that involves looking
          // up nbr points that are not on the current rank. We could also delay lookups
          // associated with the wrank, but these additional self-messages should be fast,
          // and the algorithm is much simpler without. Also, we could be redistributing to
          // gpu eventually, so do not consider any special actions for wrank...
          const int rank = rborp_v[ror].first; // rank is where the nbr-point-query has to happen
          if (rank != mpi_rank) {
            send_count[rank] += 4; // rank, bits, icv, wrank
          }
        }
      }

      send_disp[0] = 0;
      for (int rank = 1; rank < mpi_size; ++rank)
        send_disp[rank] = send_count[rank-1] + send_disp[rank-1];
      const int send_count_sum = send_disp[mpi_size-1] + send_count[mpi_size-1];

      if (send_buf_int == NULL) {
        send_count_int_sum_max = send_count_sum;
        send_buf_int = new int[send_count_int_sum_max];
      }
      else if (send_count_sum > send_count_int_sum_max) {
        send_count_int_sum_max = send_count_sum;
        delete[] send_buf_int;
        send_buf_int = new int[send_count_int_sum_max];
      }

      // in this case, the double count and the int count are the same...

      if (send_buf_double == NULL) {
        send_count_double_sum_max = send_count_sum;
        send_buf_double = new double[send_count_double_sum_max];
      }
      else if (send_count_sum > send_count_double_sum_max) {
        send_count_double_sum_max = send_count_sum;
        delete[] send_buf_double;
        send_buf_double = new double[send_count_double_sum_max]; // same size (4): x_vd, delta_vd
      }

      // and pack the send buffers...

      for (int irp = 0; irp < nrp; ++irp) {
        const int wrank = wrorp[irp];
        const int icv = cvorp[irp];
        for (int ror = rborp_i[irp]; ror != rborp_i[irp+1]; ++ror) {
          const int rank = rborp_v[ror].first;
          if (rank != mpi_rank) {
            const int bits = rborp_v[ror].second;
            send_buf_int[send_disp[rank]  ] = mpi_rank;
            send_buf_int[send_disp[rank]+1] = bits;
            send_buf_int[send_disp[rank]+2] = icv;
            send_buf_int[send_disp[rank]+3] = wrank; // where the resulting query has to go eventually
            send_buf_double[send_disp[rank]  ] = x_vd[icv][0];
            send_buf_double[send_disp[rank]+1] = x_vd[icv][1];
            send_buf_double[send_disp[rank]+2] = x_vd[icv][2];
            if (bits) PeriodicData::periodicTranslate(send_buf_double+send_disp[rank],1,bits);
            send_buf_double[send_disp[rank]+3] = delta_vd[icv];
            send_disp[rank] += 4;
          }
        }
      }

      // rewind send_disp...

      send_disp[0] = 0;
      for (int rank = 1; rank < mpi_size; ++rank)
        send_disp[rank] = send_count[rank-1] + send_disp[rank-1];

      // set up recv-side stuff...

      MPI_Alltoall(send_count,1,MPI_INT,recv_count,1,MPI_INT,mpi_comm);

      recv_disp[0] = 0;
      for (int rank = 1; rank < mpi_size; ++rank)
        recv_disp[rank] = recv_count[rank-1] + recv_disp[rank-1];
      int recv_count_sum = recv_disp[mpi_size-1] + recv_count[mpi_size-1];

      if (recv_buf_int == NULL) {
        recv_count_int_sum_max = recv_count_sum;
        recv_buf_int = new int[recv_count_int_sum_max];
      }
      else if (recv_count_sum > recv_count_int_sum_max) {
        recv_count_int_sum_max = recv_count_sum;
        delete[] recv_buf_int;
        recv_buf_int = new int[recv_count_int_sum_max];
      }

      MPI_Alltoallv(send_buf_int,send_count,send_disp,MPI_INT,
          recv_buf_int,recv_count,recv_disp,MPI_INT,mpi_comm);

      if (recv_buf_double == NULL) {
        recv_count_double_sum_max = recv_count_sum;
        recv_buf_double = new double[recv_count_double_sum_max];
      }
      else if (recv_count_sum > recv_count_double_sum_max) {
        recv_count_double_sum_max = recv_count_sum;
        delete[] recv_buf_double;
        recv_buf_double = new double[recv_count_double_sum_max];
      }

      MPI_Alltoallv(send_buf_double,send_count,send_disp,MPI_DOUBLE,
          recv_buf_double,recv_count,recv_disp,MPI_DOUBLE,mpi_comm);

      // -----------------------------------------------------------------
      // the second set of messages sends everything required to compute
      // the voronoi diagrams to the work ranks.
      //
      // unpack the recv side and start the nbr message counts. These
      // first messages are simply as follows:
      // int: rank,icv,nnbr,(rb_nbr0,index_nbr0),(rb_nbr1,index_nbr1)...
      // double: xp_nbr0,xp_nbr1...
      // -----------------------------------------------------------------

      FOR_RANK send_count[rank]        = 0; // used to count ints
      FOR_RANK send_count_double[rank] = 0;

      // these messages are the ones that result from the previous send/recv...

      assert(recv_count_sum%4 == 0); // the recv buffers are divied into blocks of 4
      const int nrecv = recv_count_sum/4;
      if (nborc_i == NULL) {
        nrecv_p1_max = nrecv+1;
        nborc_i = new int[nrecv_p1_max]; // nbr of recv-block
      }
      else if (nrecv+1 > nrecv_p1_max) {
        nrecv_p1_max = nrecv+1;
        delete[] nborc_i;
        nborc_i = new int[nrecv_p1_max]; // nbr of recv-block
      }
      nborc_i[0] = 0;
      assert(nborc_v.empty()); // nbr and bits of recv-block...
      for (int irecv = 0; irecv < nrecv; ++irecv) {
        const int rank  = recv_buf_int[irecv*4];
        assert(rank != mpi_rank);
        const int bits  = recv_buf_int[irecv*4+1];
        //const int index = recv_buf_int[irecv*4+2];
        const int wrank = recv_buf_int[irecv*4+3];
        if (bits) {
          // the point is already translated in send_* above...
          addNearbyNbrs(nborc_v,recv_buf_double+irecv*4,recv_buf_double[irecv*4+3],BitUtils::flipPeriodicBits(bits));
        }
        else {
          addNearbyNbrs(nborc_v,recv_buf_double+irecv*4,recv_buf_double[irecv*4+3],0);
        }
        nborc_i[irecv+1] = nborc_v.size();
        // now we can count the message sizes. Here we only add a message if some
        // nbrs were found...
        if (nborc_i[irecv+1]-nborc_i[irecv] > 0) {
          // ints...
          send_count[wrank] += 3 + 2*(nborc_i[irecv+1]-nborc_i[irecv]); // move to 3?
          // doubles: just coordinates...
          send_count_double[wrank] += 3*(nborc_i[irecv+1]-nborc_i[irecv]);
        }
      }

      // ----------------------------------------------------------
      // the second type of message will look as follows...
      // ints: rank,icv,nst,nnbr,(ist0,ist0p),(ist1,ist1p),...,inbr0,inbr1,inbr2 (note these are just indices, rank is known, bits are zero)...
      // doubles: xp,delta,xp_nbr0,xp_nbr1,xp_nbr2...
      // ----------------------------------------------------------

      if (nborp_i == NULL) nborp_i = new int[nrp0+1]; // nbr of rebuild points. use un-limited to ensure size decreases on subsequent iters
      nborp_i[0] = 0;
      assert(nborp_v.empty());
      for (int irp = 0; irp < nrp; ++irp) {
        const int wrank = wrorp[irp]; // work-rank
        const int icv = cvorp[irp];
        // this active point has a work-rank different than the current rank (i.e.
        // mpi_rank where the active point lives), so add its local nbrs to the
        // nborp_v array. These don't require rank/bits...
        addNearbyNbrs(nborp_v,icv,x_vd[icv],delta_vd[icv],0);
        // there may also be some local
        for (int ror = rborp_i[irp]; ror != rborp_i[irp+1]; ++ror) {
          const int rank = rborp_v[ror].first;
          if (rank == mpi_rank) {
            // must be periodicity...
            const int bits = rborp_v[ror].second; assert(bits != 0);
            double xp_t[3]; FOR_I3 xp_t[i] = x_vd[icv][i];
            PeriodicData::periodicTranslate(xp_t,1,bits);
            addNearbyNbrs(nborp_v,xp_t,delta_vd[icv],BitUtils::flipPeriodicBits(bits));
          }
        }
        nborp_i[irp+1] = nborp_v.size();
        // and the message sizes...
        send_count[wrank] += 4 + 2*(storp_i[irp+1]-storp_i[irp]) + 2*(nborp_i[irp+1]-nborp_i[irp]);
        send_count_double[wrank] += 4 + 3*(nborp_i[irp+1]-nborp_i[irp]); // x_vd, delta_vd, then nbrs like above message
      }
      rborp_v.clear();

      //MiscUtils::dumpRange(nborp_i[nrp],"nborp_i[nrp]");

      // counting is done. Now pack...

      send_disp[0] = 0;
      for (int rank = 1; rank < mpi_size; ++rank)
        send_disp[rank] = send_count[rank-1] + send_disp[rank-1];
      int send_count_int_sum = send_disp[mpi_size-1] + send_count[mpi_size-1];

      assert(send_buf_int);
      if (send_count_int_sum > send_count_int_sum_max) {
        send_count_int_sum_max = send_count_int_sum;
        delete[] send_buf_int;
        send_buf_int = new int[send_count_int_sum_max];
      }

      //MiscUtils::dumpRange(send_count_int_sum_max,"send_count_int_sum_max");

      //cout << "rank: " << mpi_rank << " send_count_int_sum_max: " << send_count_int_sum_max << " send_count_int_sum: " << send_count_int_sum << endl;

      send_disp_double[0] = 0;
      for (int rank = 1; rank < mpi_size; ++rank)
        send_disp_double[rank] = send_count_double[rank-1] + send_disp_double[rank-1];
      int send_count_double_sum = send_disp_double[mpi_size-1] + send_count_double[mpi_size-1];

      assert(send_buf_double);
      if (send_count_double_sum > send_count_double_sum_max) {
        send_count_double_sum_max = send_count_double_sum;
        delete[] send_buf_double;
        send_buf_double = new double[send_count_double_sum_max]; // same size (4): x_vd, delta_vd
      }

      // and pack the send buffers...

      for (int irecv = 0; irecv < nrecv; ++irecv) {
        if (nborc_i[irecv+1]-nborc_i[irecv] > 0) {
          const int rank  = recv_buf_int[irecv*4];
          assert((rank >= 0)&&(rank < mpi_size));
          assert(rank != mpi_rank);
          const int bits  = recv_buf_int[irecv*4+1];
          const int icv   = recv_buf_int[irecv*4+2];
          const int wrank = recv_buf_int[irecv*4+3];
          send_buf_int[send_disp[wrank]  ] = -rank-1; // -1 indexing to indicate it is a nborc
          send_buf_int[send_disp[wrank]+1] = icv;
          send_buf_int[send_disp[wrank]+2] = nborc_i[irecv+1]-nborc_i[irecv];
          send_disp[wrank] += 3;
          for (int nor = nborc_i[irecv]; nor != nborc_i[irecv+1]; ++nor) {
            const int icv_nbr = nborc_v[nor].first;
            const int icv_bits = nborc_v[nor].second; assert(icv_bits == BitUtils::flipPeriodicBits(bits));
            send_buf_int[send_disp[wrank]  ] = BitUtils::packRankBits(mpi_rank,icv_bits);
            send_buf_int[send_disp[wrank]+1] = icv_nbr;
            send_disp[wrank] += 2;
            FOR_I3 send_buf_double[send_disp_double[wrank]+i] = x_vd[icv_nbr][i];
            if (icv_bits) PeriodicData::periodicTranslate(send_buf_double+send_disp_double[wrank],1,icv_bits); // i.e. translate into ip's reference frame
            send_disp_double[wrank] += 3;
          }
        }
      }
      nborc_v.clear();

      for (int irp = 0; irp < nrp; ++irp) {
        const int wrank = wrorp[irp]; // work-rank
        assert((wrank >= 0)&&(wrank < mpi_size));
        const int icv = cvorp[irp];
        send_buf_int[send_disp[wrank]  ] = mpi_rank; // + indexing to indicate it is a nborp
        send_buf_int[send_disp[wrank]+1] = icv;
        send_buf_int[send_disp[wrank]+2] = storp_i[irp+1]-storp_i[irp];
        send_buf_int[send_disp[wrank]+3] = nborp_i[irp+1]-nborp_i[irp];
        send_disp[wrank] += 4;
        // point xp and delta...
        FOR_I3 send_buf_double[send_disp_double[wrank]+i] = x_vd[icv][i];
        send_buf_double[send_disp_double[wrank]+3]        = delta_vd[icv];
        send_disp_double[wrank] += 4;
        // local tris...
        for (int sor = storp_i[irp]; sor != storp_i[irp+1]; ++sor) {
          send_buf_int[send_disp[wrank]  ] = storp_v[sor].first; // ist
          send_buf_int[send_disp[wrank]+1] = storp_v[sor].second; // ipart/bits
          send_disp[wrank] += 2;
        }
        // local nbr coordinates ()...
        for (int nor = nborp_i[irp]; nor != nborp_i[irp+1]; ++nor) {
          const int icv_nbr = nborp_v[nor].first;
          const int bits   = nborp_v[nor].second; // already contains the inverse (i.e. how we got here)
          assert((icv_nbr != icv)||(bits != 0));
          send_buf_int[send_disp[wrank]  ] = BitUtils::packRankBits(mpi_rank,bits);
          send_buf_int[send_disp[wrank]+1] = icv_nbr;
          send_disp[wrank] += 2;
          FOR_I3 send_buf_double[send_disp_double[wrank]+i] = x_vd[icv_nbr][i];
          if (bits) PeriodicData::periodicTranslate(send_buf_double+send_disp_double[wrank],1,bits);
          assert(DIST(x_vd[icv],send_buf_double+send_disp_double[wrank]) < delta_vd[icv]*1.000000001);
          send_disp_double[wrank] += 3;
        }
      }
      storp_v.clear();
      nborp_v.clear();
      cvorp.clear();

      // reset the disps...

      send_disp[0] = 0;
      for (int rank = 1; rank < mpi_size; ++rank)
        send_disp[rank] = send_count[rank-1] + send_disp[rank-1];

      send_disp_double[0] = 0;
      for (int rank = 1; rank < mpi_size; ++rank)
        send_disp_double[rank] = send_count_double[rank-1] + send_disp_double[rank-1];

      // now send to workers...
      // ints first...

      MPI_Alltoall(send_count,1,MPI_INT,recv_count,1,MPI_INT,mpi_comm);

      recv_disp[0] = 0;
      for (int rank = 1; rank < mpi_size; ++rank)
        recv_disp[rank] = recv_count[rank-1] + recv_disp[rank-1];
      int recv_count_int_sum = recv_disp[mpi_size-1] + recv_count[mpi_size-1];

      assert(recv_buf_int != NULL);
      if (recv_count_int_sum > recv_count_int_sum_max) {
        recv_count_int_sum_max = recv_count_int_sum;
        delete[] recv_buf_int;
        recv_buf_int = new int[recv_count_int_sum_max];
        if (recv_buf_int == NULL)
          cout << "Warning: new failed: " << recv_count_int_sum_max << endl;
      }

      //MiscUtils::dumpRange(recv_count_int_sum_max,"recv_count_int_sum_max");

      MPI_Alltoallv(send_buf_int,send_count,send_disp,MPI_INT,
          recv_buf_int,recv_count,recv_disp,MPI_INT,mpi_comm);

      // doubles...

      MPI_Alltoall(send_count_double,1,MPI_INT,recv_count,1,MPI_INT,mpi_comm);

      recv_disp[0] = 0;
      for (int rank = 1; rank < mpi_size; ++rank)
        recv_disp[rank] = recv_count[rank-1] + recv_disp[rank-1];
      int recv_count_double_sum = recv_disp[mpi_size-1] + recv_count[mpi_size-1];

      assert(recv_buf_double != NULL);
      if (recv_count_double_sum > recv_count_double_sum_max) {
        recv_count_double_sum_max = recv_count_double_sum;
        delete[] recv_buf_double;
        recv_buf_double = new double[recv_count_double_sum_max];
      }

      MPI_Alltoallv(send_buf_double,send_count_double,send_disp_double,MPI_DOUBLE,
          recv_buf_double,recv_count,recv_disp,MPI_DOUBLE,mpi_comm);

      // now the data for building the vds is load balanced across the wrank's...

      // everybody ensure they can look up their faces quickly based
      // on the global surface bbox expanded to the size required by the
      // surface lookups. We do the surface part first because the surface count
      // is required for the Voronoi cost model...

      // now the data for building the vds is load balanced across the wrank's...

      assert(vdVec.empty());
      assert(nbrMultimap.empty());
      int irecv_int = 0;
      int irecv_double = 0;
      while (irecv_int < recv_count_int_sum) {
        const int rank = recv_buf_int[irecv_int];
        const int icv  = recv_buf_int[irecv_int+1];
        if (rank < 0) {
          const int nnbr = recv_buf_int[irecv_int+2];
          // this is a set of neighbors that need to be cut against the
          // point (rank,icv) which should be elsewhere in this message...
          nbrMultimap.insert(pair<pair<int,int>,pair<int,int> >(pair<int,int>(-rank-1,icv),pair<int,int>(irecv_int,irecv_double)));
          irecv_int += 3 + 2*nnbr;
          irecv_double += 3*nnbr;
        }
        else {
          const int nst  = recv_buf_int[irecv_int+2];
          const int nnbr = recv_buf_int[irecv_int+3];
          // this is a vd that needs to be cut, including some tris and nbrs...
          vdVec.push_back(pair<int,int>(irecv_int,irecv_double));
          irecv_int += 4 + 2*nst + 2*nnbr;
          irecv_double += 4 + 3*nnbr;
        }
      }
      assert(irecv_int == recv_count_int_sum);
      assert(irecv_double == recv_count_double_sum);

      // now vdVec stores the pair of irecv offsets for all voronoi diagrams that
      // need to be built. The multimap stores other nbrs from other ranks that
      // need to be cut against the vd. We should use all of these in building the
      // vd...

      // check...
      {
        int8 my_nnbr_sum = 0;
        int my_nnbr_max = 0;
        int multimap_count = 0;
        for (int ii = 0; ii < vdVec.size(); ++ii) {
          irecv_int = vdVec[ii].first;
          irecv_double = vdVec[ii].second;
          double this_x_vd[3]; FOR_I3 this_x_vd[i] = recv_buf_double[irecv_double+i];
          const double this_delta_vd               = recv_buf_double[irecv_double+3];
          const int rank = recv_buf_int[irecv_int];
          const int icv  = recv_buf_int[irecv_int+1];
          int nnbr = recv_buf_int[irecv_int+3];
          for (int inbr = 0; inbr < nnbr; ++inbr) {
            double this_x_vd_nbr[3]; FOR_I3 this_x_vd_nbr[i] = recv_buf_double[irecv_double+4+inbr*3+i];
            assert(DIST(this_x_vd,this_x_vd_nbr) < this_delta_vd*1.000000001);
          }
          multimap<pair<int,int>,pair<int,int> >::iterator iter = nbrMultimap.find(pair<int,int>(rank,icv));
          while ((iter != nbrMultimap.end())&&(iter->first.first == rank)&&(iter->first.second == icv)) {
            ++multimap_count;
            irecv_int = iter->second.first;
            irecv_double = iter->second.second;
            assert(rank == -recv_buf_int[irecv_int]-1);
            assert(icv  == recv_buf_int[irecv_int+1]);
            const int this_nnbr = recv_buf_int[irecv_int+2];
            for (int inbr = 0; inbr < this_nnbr; ++inbr) {
              double this_x_vd_nbr[3]; FOR_I3 this_x_vd_nbr[i] = recv_buf_double[irecv_double+inbr*3+i]; // no offset for icv
              assert(DIST(this_x_vd,this_x_vd_nbr) < this_delta_vd*1.000000001);
            }
            nnbr += this_nnbr;
            ++iter;
          }
          my_nnbr_sum += nnbr;
          my_nnbr_max = max(my_nnbr_max,nnbr);
        }
        assert(multimap_count == nbrMultimap.size());

        if (step%check_interval == 0) {
          int8 my_i8buf[2] = { (int8)vdVec.size(), my_nnbr_sum };
          int8 i8buf[2];
          MPI_Reduce(my_i8buf,i8buf,2,MPI_INT8,MPI_SUM,0,mpi_comm);
          int my_ibuf[1]= { my_nnbr_max };
          int ibuf[1];
          MPI_Reduce(my_ibuf,ibuf,1,MPI_INT,MPI_MAX,0,mpi_comm);
          if (mpi_rank == 0) {
            cout << "   > build vd stats: count=" << i8buf[0] <<
              " nnbr/vd (avg:max): " << double(i8buf[1])/double(i8buf[0]) << " " << ibuf[0] << endl;
          }
        }
      }

#ifdef WITH_CHRONO2
      auto t3 = std::chrono::high_resolution_clock::now();
#endif

      // ==============================================================
      // now build all vd's on the work ranks...
      // ==============================================================

      assert(intVec.empty());
      assert(doubleVec.empty());
      assert(vrdVec.empty());

      double my_max_time = 0.0;
      double my_time_sum = 0.0;

      // for debugging...
      //int one = 1;
      //int zero = 0;
      //int global_max;

      int8 my_count[4] = { 0, 0, 0, 0 };
      double wtime0 = MPI_Wtime();
      for (int iv = 0; iv < vdVec.size(); ++iv) {

        // this is not formally a synchronization point, but we can make it one by
        // allreducing 1 once we are out...
        /*
           MPI_Allreduce(&one,&global_max,1,MPI_INT,MPI_MAX,mpi_comm);
           if (mpi_rank == 0)
           cout << " > iv: " << iv << endl;
           */

#ifdef WITH_CHRONO2
        auto tstart = std::chrono::high_resolution_clock::now();
#endif

        // offsets into the recv arrays...
        int irecv_int = vdVec[iv].first;
        int irecv_double = vdVec[iv].second;

        // the geometry for this VD...
        const double * const this_x_vd = recv_buf_double + irecv_double;
        const double this_delta_vd     = recv_buf_double[irecv_double+3];
        const int rank = recv_buf_int[irecv_int];
        const int icv  = recv_buf_int[irecv_int+1];
        const int nst  = recv_buf_int[irecv_int+2];
        const int nnbr_local = recv_buf_int[irecv_int+3];

        ++debug_count;
        bool debug = false;
        bool debug2 = false;
        if (Param * param = getParam("DEBUG_RANK_ICV")) {
          const int debug_rank = param->getInt(0);
          const int debug_icv = param->getInt(1);
          if ((mpi_rank == debug_rank)&&(icv == debug_icv)) {
            cout << "DEBUG_RANK_ICV: got match at this_x_vd: " << COUT_VEC(this_x_vd) << " this_delta_vd: " << this_delta_vd << " debug_count: " << debug_count << endl;
          }
        }
        if (Param * param = getParam("DEBUG_RANK_COUNT")) {
          const int debug_rank = param->getInt(0);
          const int count = param->getInt(1);
          if ((mpi_rank == debug_rank)&&(count == debug_count)) {
            cout << "DEBUG_RANK_COUNT: got match at this_x_vd: " << COUT_VEC(this_x_vd) << " this_delta_vd: " << this_delta_vd << " icv: " << icv << endl;
            debug2 = true;
          }
        }

        if (debug) cout << "DEBUG: got cv at x_vd_debug: " << COUT_VEC(this_x_vd) << " mpi_rank, debug_count: " << mpi_rank << " " << debug_count << endl;

        // =============================================
        // start with the local surface patch...
        // XXXXX TODO: use surface adt here eventually...
        // =============================================

        cvd.clear();
        assert(cvd.nno == 0);
        assert(cvd.ned == 0);

        double t1 = 0.0;
        double t2 = 0.0;

        if (b_surface) {

          // build the patch...

          assert(nodeMap.empty());
          assert(edgeMap2.empty());

          for (int ist = 0; ist < nst; ++ist) {
            const int ist_global = recv_buf_int[irecv_int+4+ist*2]; // the actual ist on the partVec[ipart]->surface
            const int ipart_bits = recv_buf_int[irecv_int+4+ist*2+1];
            const int ipart = (ipart_bits>>6);
            const int bits = (ipart_bits&MASK_6BITS);
            assert(partVec[ipart]->surface);
            // step 1: set pbsost (part-bits-isp of st)...
            FOR_I3 {
              const int isp_global = partVec[ipart]->surface->spost[ist_global][i];
              int bits_active,isp_active;
              BitUtils::unpackPbiHash(bits_active,isp_active,partVec[ipart]->surface->getPbi(isp_global));
              // this is the total transform applied to isp_active
              const int bits_tot = BitUtils::addPeriodicBits(bits,bits_active);
              pbsost[i] = ( (uint8((ipart<<12)|bits_tot)<<32) | uint8(isp_active) );
            }
            // step 2: use edge-matching to set as much of the noost and edost that we can...
            FOR_I3 noost[i] = -1;
            FOR_I3 {
              edost[i] = -1;
              // here the edgeMap2 is used to find edges in terms of surface node pairs...
              map<const pair<uint8,uint8>,int>::iterator iter = edgeMap2.find(pair<uint8,uint8>(pbsost[(i+1)%3],pbsost[i]));
              if (iter != edgeMap2.end()) {
                edost[i] = iter->second;
                edgeMap2.erase(iter);
                if (noost[i] == -1) noost[i] = cvd.nooed[edost[i]][1];
                else assert(noost[i] == cvd.nooed[edost[i]][1]);
                if (noost[(i+1)%3] == -1) noost[(i+1)%3] = cvd.nooed[edost[i]][0];
                else assert(noost[(i+1)%3] == cvd.nooed[edost[i]][0]);
              }
            }
            // step 3: now noost and edost are set from edges where possible. At this point, complete any nodes not set...
            FOR_I3 {
              if (noost[i] == -1) {
                map<const uint8,int>::iterator iter = nodeMap.find(pbsost[i]);
                if (iter == nodeMap.end()) {
                  const int ino = cvd.new_node();
                  noost[i] = ino;
                  nodeMap[pbsost[i]] = ino;
                  const int ipart_bits = (pbsost[i]>>32);
                  const int bits = (ipart_bits&MASK_12BITS);
                  const int isp_active = (pbsost[i]&MASK_32BITS);

                  if (bits) {
                    double xp_t[3]; FOR_J3 xp_t[j] = partVec[ipart]->surface->xsp[isp_active][j];
                    PeriodicData::periodicTranslate(xp_t,1,bits);
                    FOR_J3 cvd.x_no[ino][j] = xp_t[j] - this_x_vd[j];
                  }
                  else {
                    FOR_J3 cvd.x_no[ino][j] = partVec[ipart]->surface->xsp[isp_active][j] - this_x_vd[j];
                  }
                }
                else {
                  noost[i] = iter->second;
                }
              }
            }
            // step 4: finally complete the edges...
            FOR_I3 {
              if (edost[i] == -1) {
                const int ied = cvd.new_edge();
                edgeMap2[pair<uint8,uint8>(pbsost[i],pbsost[(i+1)%3])] = ied;
                cvd.nooed[ied][0] = noost[i];
                cvd.nooed[ied][1] = noost[(i+1)%3];
                cvd.faoed[ied][0] = irecv_int+4+ist*2; // reference to recv_buf_int ist_global/part/bits...
                cvd.faoed[ied][1] = -1;
              }
              else {
                const int ied = edost[i];
                assert(cvd.nooed[ied][1] == noost[i]);
                assert(cvd.nooed[ied][0] == noost[(i+1)%3]);
                assert(cvd.faoed[ied][0] != irecv_int+4+ist*2); // no edge can have the same tri on both sides
                assert(cvd.faoed[ied][1] == -1);
                cvd.faoed[ied][1] = irecv_int+4+ist*2;
              }
            }

          }

          nodeMap.clear();
          edgeMap2.clear();

        }

        // =====================================================================================
        // now we may or may not have a surface in cvd. Time to build the seed. If there is no
        // surface because surface flag is false, or no tris were grabbed, this just returns a
        // cube. Otherwise it does potentially a lot of work to return a valid positive volume
        // =====================================================================================

        try {

          cvd.buildSeed(this_delta_vd,recv_buf_int,debug2);

          assert(!cvd.empty());
          for (int ied = 0; ied < cvd.ned; ++ied) {
            FOR_I2 {
              if (cvd.faoed[ied][i] == -1)
                throw(-1);
            }
          }

        }
        catch(...) {

          // if we fail, rebuild the case and write it out in binary to directly diagnose using the
          // CVD_DEBUG override in the moving solver...

          cout << "FAILED cvd.buildSeed at x_vd: " << COUT_VEC(this_x_vd) << " rank: " << mpi_rank << endl;
          cvd.clear();
          assert(cvd.nno == 0);
          assert(cvd.ned == 0);
          if (b_surface) {
            // build the patch...
            assert(nodeMap.empty());
            assert(edgeMap2.empty());
            for (int ist = 0; ist < nst; ++ist) {
              const int ist_global = recv_buf_int[irecv_int+4+ist*2]; // the actual ist on the partVec[ipart]->surface
              const int ipart_bits = recv_buf_int[irecv_int+4+ist*2+1];
              const int ipart = (ipart_bits>>6);
              const int bits = (ipart_bits&MASK_6BITS);
              assert(partVec[ipart]->surface);
              // step 1: set pbsost (part-bits-isp of st)...
              FOR_I3 {
                const int isp_global = partVec[ipart]->surface->spost[ist_global][i];
                int bits_active,isp_active;
                BitUtils::unpackPbiHash(bits_active,isp_active,partVec[ipart]->surface->getPbi(isp_global));
                // this is the total transform applied to isp_active
                const int bits_tot = BitUtils::addPeriodicBits(bits,bits_active);
                pbsost[i] = ( (uint8((ipart<<12)|bits_tot)<<32) | uint8(isp_active) );
              }
              // step 2: use edge-matching to set as much of the noost and edost that we can...
              FOR_I3 noost[i] = -1;
              FOR_I3 {
                edost[i] = -1;
                // here the edgeMap2 is used to find edges in terms of surface node pairs...
                map<const pair<uint8,uint8>,int>::iterator iter = edgeMap2.find(pair<uint8,uint8>(pbsost[(i+1)%3],pbsost[i]));
                if (iter != edgeMap2.end()) {
                  edost[i] = iter->second;
                  edgeMap2.erase(iter);
                  if (noost[i] == -1) noost[i] = cvd.nooed[edost[i]][1];
                  else assert(noost[i] == cvd.nooed[edost[i]][1]);
                  if (noost[(i+1)%3] == -1) noost[(i+1)%3] = cvd.nooed[edost[i]][0];
                  else assert(noost[(i+1)%3] == cvd.nooed[edost[i]][0]);
                }
              }
              // step 3: now noost and edost are set from edges where possible. At this point, complete any nodes not set...
              FOR_I3 {
                if (noost[i] == -1) {
                  map<const uint8,int>::iterator iter = nodeMap.find(pbsost[i]);
                  if (iter == nodeMap.end()) {
                    const int ino = cvd.new_node();
                    noost[i] = ino;
                    nodeMap[pbsost[i]] = ino;
                    const int ipart_bits = (pbsost[i]>>32);
                    const int bits = (ipart_bits&MASK_12BITS);
                    const int isp_active = (pbsost[i]&MASK_32BITS);

                    if (bits) {
                      double xp_t[3]; FOR_J3 xp_t[j] = partVec[ipart]->surface->xsp[isp_active][j];
                      PeriodicData::periodicTranslate(xp_t,1,bits);
                      FOR_J3 cvd.x_no[ino][j] = xp_t[j] - this_x_vd[j];
                    }
                    else {
                      FOR_J3 cvd.x_no[ino][j] = partVec[ipart]->surface->xsp[isp_active][j] - this_x_vd[j];
                    }
                  }
                  else {
                    noost[i] = iter->second;
                  }
                }
              }
              // step 4: finally complete the edges...
              FOR_I3 {
                if (edost[i] == -1) {
                  const int ied = cvd.new_edge();
                  edgeMap2[pair<uint8,uint8>(pbsost[i],pbsost[(i+1)%3])] = ied;
                  cvd.nooed[ied][0] = noost[i];
                  cvd.nooed[ied][1] = noost[(i+1)%3];
                  cvd.faoed[ied][0] = irecv_int+4+ist*2; // reference to recv_buf_int ist_global/part/bits...
                  cvd.faoed[ied][1] = -1;
                }
                else {
                  const int ied = edost[i];
                  assert(cvd.nooed[ied][1] == noost[i]);
                  assert(cvd.nooed[ied][0] == noost[(i+1)%3]);
                  assert(cvd.faoed[ied][0] != irecv_int+4+ist*2); // no edge can have the same tri on both sides
                  assert(cvd.faoed[ied][1] == -1);
                  cvd.faoed[ied][1] = irecv_int+4+ist*2;
                }
              }
            }
            nodeMap.clear();
            edgeMap2.clear();
          }
          // write the failing problem...
          char filename[128];
          sprintf(filename,"cvd_debug-%06d.bin",mpi_rank);
          cvd.writeBinary(filename);
          double zero[3] = { 0.0, 0.0, 0.0 };
          cvd.writeTecplot(mpi_rank,zero);
          sprintf(filename,"cvd_debug_support-%06d.bin",mpi_rank);
          FILE * fp = fopen(filename,"wb");
          fwrite(&this_delta_vd,sizeof(double),1,fp);
          fwrite(&recv_count_int_sum,sizeof(int),1,fp);
          fwrite(recv_buf_int,sizeof(int),recv_count_int_sum,fp);
          fclose(fp);
          // die always...
          throw(-1);
        }

        double t3 = MPI_Wtime()-wtime0;

        // d2 max is the maximum distance of any node relative the x_vd...
        cvd.setD2Max();

#ifdef WITH_CHRONO2
        auto tseed = std::chrono::high_resolution_clock::now();
#endif

        // start with the nbrs that are local to this VD back on the owning rank (recall
        // we are on the work-rank now)...
        assert(nbrVec.empty());
        assert(nnbr_local == recv_buf_int[irecv_int+3]);
        for (int inbr = 0; inbr < nnbr_local; ++inbr) {
          const double * const this_x_vd_nbr = recv_buf_double + irecv_double+4+inbr*3;
          const double d2 = DIST2(this_x_vd,this_x_vd_nbr);
          assert(sqrt(d2) < this_delta_vd*1.000000001);
          if (!(sqrt(d2) > this_delta_vd*1.0E-5)) {
            cout << "sqrt(d2): " << sqrt(d2) << " this_x_vd: " << COUT_VEC(this_x_vd) << " this_x_vd_nbr: " << COUT_VEC(this_x_vd_nbr) << endl;
          }
          assert(sqrt(d2) > this_delta_vd*1.0E-5); // if we hit this, a nbr is right on top of cv -- check it out
          nbrVec.push_back( pair<double,pair<int,int> >(d2,pair<int,int>(irecv_int+4+nst*2+inbr*2,irecv_double+4+inbr*3)) );
        }
        multimap<pair<int,int>,pair<int,int> >::iterator iter = nbrMultimap.find(pair<int,int>(rank,icv));
        while ((iter != nbrMultimap.end())&&(iter->first.first == rank)&&(iter->first.second == icv)) {
          irecv_int = iter->second.first;
          irecv_double = iter->second.second;
          assert(rank == -recv_buf_int[irecv_int]-1);
          assert(icv  == recv_buf_int[irecv_int+1]);
          const int this_nnbr = recv_buf_int[irecv_int+2];
          for (int inbr = 0; inbr < this_nnbr; ++inbr) {
            const double * const this_x_vd_nbr = recv_buf_double + irecv_double+inbr*3; // no offset for icv
            const double d2 = DIST2(this_x_vd,this_x_vd_nbr);
            assert(sqrt(d2) < this_delta_vd*1.000000001);
            assert(sqrt(d2) > this_delta_vd*1.0E-5); // if we hit this, a nbr is right on top of cv -- check it out
            nbrVec.push_back( pair<double,pair<int,int> >(d2,pair<int,int>(irecv_int+3+inbr*2,irecv_double+inbr*3)) );
          }
          ++iter;
        }

        // now sort based on proximity...
        sort(nbrVec.begin(),nbrVec.end());

        // we now have a sorted list of nbrs. Cut the cvd against the nbrs until it cannot possibly
        // be cut anymore...
        for (int in = 0; in < nbrVec.size(); ++in) {
          // recall that the .first contains the d2 nbr distance. and cvd.d2_max contains
          // the maximum node radius. as soon as d2 >= 4*cvd.d2_max, there cannot be any
          // more nbrs that cut this cv...
          if (nbrVec[in].first > 4.0*cvd.d2_max)
            break; // done: all remaining nbrs are too far away to possibly cut the current vd...
          double * this_x_vd_nbr = recv_buf_double + nbrVec[in].second.second;
          double dn[3]; FOR_I3 dn[i] = 0.5*(this_x_vd_nbr[i] - this_x_vd[i]);
          // only cut if we are inside the 6 paraboloids...
          if ((2.0*dn[0] > cvd.Lmax[0])||
              (2.0*dn[1] > cvd.Lmax[1])||
              (2.0*dn[2] > cvd.Lmax[2])||
              (2.0*dn[0] < cvd.Lmin[0])||
              (2.0*dn[1] < cvd.Lmin[1])||
              (2.0*dn[2] < cvd.Lmin[2]) ) {
            continue;
          }
          // if we got here, we are cutting...
          //cvd.cut(dn,-8-nbrVec[in].second.first); // note that local nbrs use a -8 convention. cvd.d2_max updated automatically
          cvd.cut(dn,-8-in); // note that local nbrs use a -8 convention. cvd.d2_max updated automatically
        }

        try {
          cvd.check();
        }
        catch(...) {
          char filename[128];
          sprintf(filename,"cvd_debug.%06d.bin",mpi_rank);
          cvd.writeBinary(filename);
          double zero[3] = { 0.0, 0.0, 0.0 };
          cvd.writeTecplot(mpi_rank,zero);
          cout << "Error: failed second check: --DEBUG_RANK_COUNT " << mpi_rank << " " << debug_count << " this_x_vd: " << COUT_VEC(this_x_vd) << endl;
          throw(-1);
        }

        double t4 = MPI_Wtime()-wtime0;

        // ----------------------------------------------------------------------------------
        // we may be done if both of the following are true:
        // 1. there is no part of the original seed surface in the current vd, and
        // 2. we have considered all possible nbrs -- i.e. delta is big enough...
        // It may still not be a complete voronoi diagram for other reasons related to
        // the surface, but this we will figure out when we try to pair faces...
        // ----------------------------------------------------------------------------------

        int kind = 0;
        if (cvd.hasSeedBoundary()) {

          // if the volume still has some of the original seed faces on it, it needs to be
          // cut by more nbrs, or made larger and cut by more neighbors... Here we could save a bit of
          // work in some cases by rebuilding it larger and cutting against the current nbrs,
          // but for now, just throw all this work out and start again...
          // if we didn't even find a single nbr, then we increase delta more rapidly than
          // the regular rate...
          if (nbrVec.size() == 0) {
            //if (debug) cout << "DEBUG: failed with no nbrs" << endl;
            ++my_count[0];
            if (mpi_rank == rank) {
              // is wrank and rank are the same, immediately process without
              // any messaging...
              assert(flag_cv[icv]&REBUILD_BIT);
              assert(flag_cv[icv]&ACTIVE_BIT);
              delta_vd[icv] *= 1.5; // increase delta by 50%, LEAVE rebuild bit set

              kind = 1;

            }
            else {
              // otherwise we need to send a message back to rank and do
              // the same thing...
              intVec.push_back(VD_FAILED_SEED_AND_NO_NBRS);
              intVec.push_back(rank);
              intVec.push_back(icv);

              kind = 2;

            }
          }
          else {
            //if (debug) cout << "DEBUG: failed with nbrs" << endl;
            ++my_count[1];
            if (mpi_rank == rank) {
              assert(flag_cv[icv]&REBUILD_BIT);
              assert(flag_cv[icv]&ACTIVE_BIT);
              delta_vd[icv] *= 1.5; // increase delta by 50%, LEAVE rebuild bit set

              kind = 3;

            }
            else {
              intVec.push_back(VD_FAILED_SEED);
              intVec.push_back(rank);
              intVec.push_back(icv);

              kind = 4;

            }
          }
        }
        else {
          // the volume has no seed boundaries left -- this is good. They were all cut away.
          // now we can guarantee that there are no other nbrs that can cut any part of
          // the volume if they are atleast 2x as far away as the furthest point on the
          // current hull is (from xv)...
          const double delta_half = sqrt(cvd.d2_max); // cvd.d2_max stores the furthest distance of any node
          if (this_delta_vd <= 2.0*delta_half) {
            // other neighbors might exist that could still cut this volume. Find them in
            // the next iteration...
            // note that we can still have problems with non-convex elements not
            // cutting enough surface. It may be that delta has to be quite a bit larger
            //if (debug) cout << "DEBUG: failed: delta not big enough" << endl;
            ++my_count[2];
            if (mpi_rank == rank) {
              // leave the status in flag_cv and adjust delta...
              assert(flag_cv[icv]&REBUILD_BIT);
              assert(flag_cv[icv]&ACTIVE_BIT);
              delta_vd[icv] = 2.001*delta_half; // could be just 2.0 here? - use a little larger so face checks pass for sure

              kind = 5;

            }
            else {
              intVec.push_back(VD_FAILED_DELTA_TOO_SMALL);
              intVec.push_back(rank);
              intVec.push_back(icv);
              doubleVec.push_back(2.001*delta_half); // see comment above -- could be just 2.0

              kind = 6;

            }
          }
          else {
            // this one is done! (unless the face matching fails, then we are back here)...
            ++my_count[3];
            // calc the geometry required...
            // if this is on our rank, then put the cvd right into the vdReturnDataVec...

            // for timing...
            //VdReturnData * vrd;

            if (debug2) {
              cout << "AFTER EVERYTHING: icv: " << icv << endl;
              cvd.writeTecplot(icv,this_x_vd);
            }

            if (mpi_rank == rank) {
              assert(flag_cv[icv]&REBUILD_BIT);
              assert(flag_cv[icv]&ACTIVE_BIT);
              flag_cv[icv] &= ~REBUILD_BIT;
              flag_cv[icv] |= NEW_BIT;
              // this is our rank, so push directly into the vdReturnDataVec...
              vdReturnDataVec.push_back(VdReturnData());
              // note: this routine translates the nbrVec indices currently in the cvd (using -8 indexing)
              // into irecv_int indices that can be converted into nbr info from the int buffer below. This is
              // not done immediately above because this routine also eliminates small faces that do
              // not need to be included using the common length scale of nbr distance (squared) available
              // in nbrVec...
              setVdReturnData3(vdReturnDataVec.back(),cvd,this_x_vd,nbrVec,recv_buf_int);
              vdReturnDataVec.back().icv = icv;
              delta_vd[icv] = 2.001*sqrt(cvd.d2_max);
              // still need to adjust the referenced nbrs - they currently refer to buffer offsets
              // do this below while we still have the buffers...

              //vrd = &vdReturnDataVec.back();

              kind = 7;

            }
            else {
              // for off-rank icv's, we use the same routine to pack separate buffer...
              intVec.push_back(VD_PASSED);
              intVec.push_back(rank);
              vrdVec.push_back(VdReturnData());
              // see comment above...
              setVdReturnData3(vrdVec.back(),cvd,this_x_vd,nbrVec,recv_buf_int);
              vrdVec.back().icv = icv;
              vrdVec.back().delta = 2.001*sqrt(cvd.d2_max);
              // no need to push_back icv into intVec -- it is in the vdReturnData...

              //vrd = &vrdVec.back();

              kind = 8;

            }

          }
        }

        {
          const double wtime = MPI_Wtime();
          const double elapsed_time = wtime-wtime0;
          my_time_buf[0] += elapsed_time;
          my_time_sum += elapsed_time;
          my_max_time = max(my_max_time,elapsed_time);
          wtime0 = wtime;
          if (elapsed_time > 30.0) {
            cout << "WRITING CVD: elapsed_time: " << elapsed_time << " t1,t2,t3,t4: " << t1 << " " << t2 << " " << t3 << " " << t4 << endl;
            //double zero[3] = { 0.0, 0.0, 0.0 };
            //cvd.writeTecplot(iv,zero);
          }
        }

#ifdef WITH_CHRONO2
        auto tend = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed_seconds = tend-tstart;
        std::chrono::duration<double> seed_seconds = tseed-tstart;
        timingVec.push_back(VdTimingData(elapsed_seconds.count(),seed_seconds.count(),nst,kind));
#endif

        // clear nbrVec for next vd...
        nbrVec.clear();

      } // vdVec loop


      {
        MPI_Barrier(mpi_comm);
        my_time_buf[1] += MPI_Wtime()-wtime0;
      }

      {
        double my_buf[2] = { (double)vdVec.size(), my_time_sum };
        double buf[2];
        MPI_Reduce(my_buf,buf,2,MPI_DOUBLE,MPI_SUM,0,mpi_comm);
        double max_time;
        MPI_Reduce(&my_max_time,&max_time,1,MPI_DOUBLE,MPI_MAX,0,mpi_comm);
        if (mpi_rank == 0) cout << "   > vd build timing: count: " << buf[0] << " avg: " << buf[1]/buf[0] << " max: " << max_time << endl;
      }

#ifdef WITH_CHRONO2
      auto t4 = std::chrono::high_resolution_clock::now();
#endif

      // at this point, all vd's have been built and assessed.

      vdVec.clear();
      nbrMultimap.clear();

      // report the status for this iteration...

      int8 count[4];
      MPI_Allreduce(my_count,count,4,MPI_INT8,MPI_SUM,mpi_comm);
      // if no one failed, we are done this iteration...
      if (count[0] + count[1] + count[2] == 0)
        done = true;
      if ((mpi_rank == 0)&&(step%check_interval == 0)) {
        cout << "   > vd failed seed and no nbrs: " << count[0] << " vd failed seed: " << count[1] <<
          " vd failed delta too small: " << count[2] << " vd passed: " << count[3] << " sum: " <<
          count[0]+count[1]+count[2]+count[3] << endl;
      }

#ifdef WITH_CHRONO2
      auto t5start = std::chrono::high_resolution_clock::now();
#endif

      // while the recv_* stuff is still available, we need to change the
      // nbr references in local vdReturnDataVec...

      for (int ivd = vdReturnDataVec_size0; ivd < vdReturnDataVec.size(); ++ivd) {
        // cycle through the nfa cvofa references and modify to an actual icv...
        for (int ifa = 0; ifa < vdReturnDataVec[ivd].nfa; ++ifa) {
          const int irecv_int = vdReturnDataVec[ivd].cvofa[ifa];
          int rank,bits;
          BitUtils::unpackRankBits(rank,bits,recv_buf_int[irecv_int]);
          const int icv_nbr   = recv_buf_int[irecv_int+1];
          if ((rank == mpi_rank)&&(bits == 0)) {
            // this is a local icv_nbr, so set it in vdReturnDataVec...
            vdReturnDataVec[ivd].cvofa[ifa] = icv_nbr;
          }
          else {
            // this is a ghost (either periodic, and/or off-rank)...
            const uint8 rbi = BitUtils::packRankBitsIndex(rank,bits,icv_nbr);
            map<const uint8,int>::iterator iter = rbiMap.find(rbi);
            if (iter == rbiMap.end()) {
              // this is a ghost not currently in the rbi...
              rbiMap[rbi] = vdReturnDataVec[ivd].cvofa[ifa] = rbi_g.size() + ncv;
              rbi_g.push_back(rbi);
            }
            else {
              // we found this one -- TODO should flag somehow...
              vdReturnDataVec[ivd].cvofa[ifa] = iter->second;
            }
          }
        }
      }

      // now send back everything...

      FOR_RANK send_count[rank] = 0;
      FOR_RANK send_count_double[rank] = 0;
      for (int iter = 0; iter < 2; ++iter) {
        int ii = 0;
        int id = 0;
        int ivrd = 0;
        // use the intVec to drive this packing...
        while (ii < intVec.size()) {
          const int status = intVec[ii++];
          const int rank = intVec[ii++];
          assert(rank != mpi_rank);
          switch (status) {
            case VD_FAILED_SEED_AND_NO_NBRS:
            case VD_FAILED_SEED:
              if (iter == 0) {
                ++ii; // skip icv when counting
                send_count[rank] += 2; // status,icv
              }
              else {
                send_buf_int[send_disp[rank]  ] = status; // VD_FAILED_SEED_AND_NO_NBRS or VD_FAILED_SEED
                send_buf_int[send_disp[rank]+1] = intVec[ii++]; // icv
                send_disp[rank] += 2;
              }
              break;
            case VD_FAILED_DELTA_TOO_SMALL:
              if (iter == 0) {
                ++ii;
                send_count[rank] += 2; // status,icv
                ++id;
                send_count_double[rank] += 1; // new delta_vd
              }
              else {
                send_buf_int[send_disp[rank]  ] = status; // VD_FAILED_DELTA_TOO_SMALL
                send_buf_int[send_disp[rank]+1] = intVec[ii++]; // icv
                send_disp[rank] += 2;
                send_buf_double[send_disp_double[rank]] = doubleVec[id++]; // new delta_vd
                send_disp_double[rank] += 1;
              }
              break;
            case VD_PASSED:
              // for now just pack the status & icv...
              if (iter == 0) {
                // recall data in vrdVec...
                /*
                   int icv,ned,nno,nfa,nbf,ngr;
                   int (*nooed)[2];
                   int (*faoed)[2];
                   int * cvofa;
                   int (*spbobf)[2];
                   int * edogr_i;
                   double delta;
                   double (*x_no)[3];
                   double (*r2_fa);
                   */
                // ints: status,icv,ned,nno,nfa,nbf,ngr, for each ned: 4, for each nfa: 2 (rank,bits),icv_nbr, for each nbf: 2, for each ngr: 1
                send_count[rank] += 7 + vrdVec[ivrd].ned*4 + vrdVec[ivrd].nfa*2 + vrdVec[ivrd].nbf*2 + vrdVec[ivrd].ngr;
                // doubles: delta, for each nno: 3 (x_no), for each nfa: 1 (r2_fa)...
                send_count_double[rank] += 1 + vrdVec[ivrd].nno*3 + vrdVec[ivrd].nfa;
                ++ivrd;
              }
              else {
                // 7 ints...
                send_buf_int[send_disp[rank]  ] = status; // VD_PASSED
                send_buf_int[send_disp[rank]+1] = vrdVec[ivrd].icv; // icv
                send_buf_int[send_disp[rank]+2] = vrdVec[ivrd].ned;
                send_buf_int[send_disp[rank]+3] = vrdVec[ivrd].nno;
                send_buf_int[send_disp[rank]+4] = vrdVec[ivrd].nfa;
                send_buf_int[send_disp[rank]+5] = vrdVec[ivrd].nbf;
                send_buf_int[send_disp[rank]+6] = vrdVec[ivrd].ngr;
                send_disp[rank] += 7;
                // 1 double...
                send_buf_double[send_disp_double[rank]] = vrdVec[ivrd].delta;
                send_disp_double[rank] += 1;
                // edge data...
                for (int ied = 0; ied < vrdVec[ivrd].ned; ++ied) {
                  send_buf_int[send_disp[rank]+ied*4  ] = vrdVec[ivrd].faoed[ied][0];
                  send_buf_int[send_disp[rank]+ied*4+1] = vrdVec[ivrd].faoed[ied][1];
                  send_buf_int[send_disp[rank]+ied*4+2] = vrdVec[ivrd].nooed[ied][0];
                  send_buf_int[send_disp[rank]+ied*4+3] = vrdVec[ivrd].nooed[ied][1];
                }
                send_disp[rank]        += vrdVec[ivrd].ned*4;
                // no data...
                for (int ino = 0; ino < vrdVec[ivrd].nno; ++ino) {
                  FOR_I3 send_buf_double[send_disp_double[rank]+ino*3+i] = vrdVec[ivrd].x_no[ino][i];
                }
                send_disp_double[rank] += vrdVec[ivrd].nno*3;
                // fa data...
                for (int ifa = 0; ifa < vrdVec[ivrd].nfa; ++ifa) {
                  // ints...
                  const int irecv_int = vrdVec[ivrd].cvofa[ifa];
                  send_buf_int[send_disp[rank]+ifa*2  ] = recv_buf_int[irecv_int]; // rank_bits
                  send_buf_int[send_disp[rank]+ifa*2+1] = recv_buf_int[irecv_int+1]; // icv_nbr
                  // doubles...
                  send_buf_double[send_disp_double[rank]+ifa] = vrdVec[ivrd].r2_fa[ifa];
                }
                send_disp[rank]        += vrdVec[ivrd].nfa*2;
                send_disp_double[rank] += vrdVec[ivrd].nfa;
                // bf data...
                for (int ibf = 0; ibf < vrdVec[ivrd].nbf; ++ibf) {
                  // ints...
                  send_buf_int[send_disp[rank]+2*ibf]   = vrdVec[ivrd].spbobf[ibf][0];
                  send_buf_int[send_disp[rank]+2*ibf+1] = vrdVec[ivrd].spbobf[ibf][1];
                }
                send_disp[rank]        += vrdVec[ivrd].nbf*2;
                // gr data...
                for (int igr = 0; igr < vrdVec[ivrd].ngr; ++igr) {
                  send_buf_int[send_disp[rank]+igr] = vrdVec[ivrd].edogr_i[igr+1] - vrdVec[ivrd].edogr_i[igr];
                }
                send_disp[rank]        += vrdVec[ivrd].ngr;
                ++ivrd; // index in the packed vrdVec
              }
              break;
            default:
              assert(0);
          }
        }
        assert(ii == intVec.size());
        assert(id == doubleVec.size());
        assert(ivrd == vrdVec.size());
        // rebuild disp on both iters...
        send_disp[0] = 0;
        for (int rank = 1; rank < mpi_size; ++rank)
          send_disp[rank] = send_count[rank-1] + send_disp[rank-1];
        send_disp_double[0] = 0;
        for (int rank = 1; rank < mpi_size; ++rank)
          send_disp_double[rank] = send_count_double[rank-1] + send_disp_double[rank-1];
        // allocate memory on first iter...
        if (iter == 0) {
          // ints...
          send_count_int_sum = send_count[mpi_size-1] + send_disp[mpi_size-1];
          assert(send_buf_int != NULL);
          if (send_count_int_sum > send_count_int_sum_max) {
            send_count_int_sum_max = send_count_int_sum;
            delete[] send_buf_int;
            send_buf_int = new int[send_count_int_sum_max];
          }
          // doubles...
          send_count_double_sum = send_count_double[mpi_size-1] + send_disp_double[mpi_size-1];
          assert(send_buf_double != NULL);
          if (send_count_double_sum > send_count_double_sum_max) {
            send_count_double_sum_max = send_count_double_sum;
            delete[] send_buf_double;
            send_buf_double = new double[send_count_double_sum_max];
          }
        }
      }

      intVec.clear();
      doubleVec.clear();
      for (int ii = 0, limit = vrdVec.size(); ii < limit; ++ii)
        vrdVec[ii].clear();
      vrdVec.clear();

      // ints...

      MPI_Alltoall(send_count,1,MPI_INT,recv_count,1,MPI_INT,mpi_comm);

      recv_disp[0] = 0;
      for (int rank = 1; rank < mpi_size; ++rank)
        recv_disp[rank] = recv_count[rank-1] + recv_disp[rank-1];
      recv_count_int_sum = recv_disp[mpi_size-1] + recv_count[mpi_size-1];

      assert(recv_buf_int);
      if (recv_count_int_sum > recv_count_int_sum_max) {
        recv_count_int_sum_max = recv_count_int_sum;
        delete[] recv_buf_int;
        recv_buf_int = new int[recv_count_int_sum_max];
      }

      MPI_Alltoallv(send_buf_int,send_count,send_disp,MPI_INT,
          recv_buf_int,recv_count,recv_disp,MPI_INT,mpi_comm);

      // doubles...

      MPI_Alltoall(send_count_double,1,MPI_INT,recv_count,1,MPI_INT,mpi_comm);

      recv_disp[0] = 0;
      for (int rank = 1; rank < mpi_size; ++rank)
        recv_disp[rank] = recv_count[rank-1] + recv_disp[rank-1];
      recv_count_double_sum = recv_disp[mpi_size-1] + recv_count[mpi_size-1];

      assert(recv_buf_double);
      if (recv_count_double_sum > recv_count_double_sum_max) {
        recv_count_double_sum_max = recv_count_double_sum;
        delete[] recv_buf_double;
        recv_buf_double = new double[recv_count_double_sum];
      }

      MPI_Alltoallv(send_buf_double,send_count_double,send_disp_double,MPI_DOUBLE,
          recv_buf_double,recv_count,recv_disp,MPI_DOUBLE,mpi_comm);

      // and unpack on the owning ranks...

      irecv_int = 0;
      irecv_double = 0;
      while (irecv_int < recv_count_int_sum) {
        const int status = recv_buf_int[irecv_int++];
        const int icv    = recv_buf_int[irecv_int++];
        assert((icv >= 0)&&(icv < ncv));
        assert(flag_cv[icv]&REBUILD_BIT);
        assert(flag_cv[icv]&ACTIVE_BIT);
        switch (status) {
          case VD_FAILED_SEED_AND_NO_NBRS:
          case VD_FAILED_SEED:
            delta_vd[icv] *= 1.5; // increase delta by 50%, LEAVE rebuild bit set
            break;
          case VD_FAILED_DELTA_TOO_SMALL:
            delta_vd[icv] = recv_buf_double[irecv_double++]; // reset delta to passed value
            break;
          case VD_PASSED:
            flag_cv[icv] &= ~REBUILD_BIT;
            flag_cv[icv] |= NEW_BIT;
            {
              vdReturnDataVec.push_back(VdReturnData());
              VdReturnData * vrd = &vdReturnDataVec.back();
              vrd->icv = icv;
              vrd->ned = recv_buf_int[irecv_int++];
              vrd->nno = recv_buf_int[irecv_int++];
              vrd->nfa = recv_buf_int[irecv_int++];
              vrd->nbf = recv_buf_int[irecv_int++];
              vrd->ngr = recv_buf_int[irecv_int++];
              delta_vd[icv] = recv_buf_double[irecv_double++];
              // edge data...
              vrd->faoed = new int[vrd->ned][2];
              vrd->nooed = new int[vrd->ned][2];
              for (int ied = 0; ied < vrd->ned; ++ied) {
                FOR_I2 vrd->faoed[ied][i] = recv_buf_int[irecv_int++];
                FOR_I2 vrd->nooed[ied][i] = recv_buf_int[irecv_int++];
              }
              // no data...
              vrd->x_no = new double[vrd->nno][3];
              for (int ino = 0; ino < vrd->nno; ++ino) {
                FOR_I3 vrd->x_no[ino][i] = recv_buf_double[irecv_double++];
              }
              // fa data...
              vrd->cvofa = new int[vrd->nfa];
              vrd->r2_fa = new double[vrd->nfa];
              for (int ifa = 0; ifa < vrd->nfa; ++ifa) {
                int rank,bits;
                BitUtils::unpackRankBits(rank,bits,recv_buf_int[irecv_int++]);
                const int icv_nbr = recv_buf_int[irecv_int++];
                if ((rank == mpi_rank)&&(bits == 0)) {
                  // this is a local icv_nbr, so set it directly...
                  vrd->cvofa[ifa] = icv_nbr;
                }
                else {
                  // this is a ghost (either periodic, and/or off-rank)...
                  const uint8 rbi = BitUtils::packRankBitsIndex(rank,bits,icv_nbr);
                  map<const uint8,int>::iterator iter = rbiMap.find(rbi);
                  if (iter == rbiMap.end()) {
                    // this is a ghost not currently in the rbi...
                    rbiMap[rbi] = vrd->cvofa[ifa] = rbi_g.size() + ncv;
                    rbi_g.push_back(rbi);
                  }
                  else {
                    // we found this one -- TODO should flag somehow...
                    vrd->cvofa[ifa] = iter->second;
                  }
                }
                // doubles...
                vrd->r2_fa[ifa] = recv_buf_double[irecv_double++];
              }
              // bf data...
              vrd->spbobf = new int[vrd->nbf][2];
              for (int ibf = 0; ibf < vrd->nbf; ++ibf) {
                FOR_I2 vrd->spbobf[ibf][i] = recv_buf_int[irecv_int++];
              }
              // gr data...
              vrd->edogr_i = new int[vrd->ngr+1];
              vrd->edogr_i[0] = 0;
              for (int igr = 0; igr < vrd->ngr; ++igr) {
                vrd->edogr_i[igr+1] = vrd->edogr_i[igr] + recv_buf_int[irecv_int++];
              }
            }
            break;
          default:
            assert(0);
        }
      }
      assert(irecv_int == recv_count_int_sum);
      assert(irecv_double == recv_count_double_sum);

#ifdef WITH_CHRONO2
      auto t5 = std::chrono::high_resolution_clock::now();

      double my_buf[4];
      std::chrono::duration<double> t12 = t2-t1; my_buf[0] = t12.count();
      std::chrono::duration<double> t23 = t3-t2; my_buf[1] = t23.count();
      std::chrono::duration<double> t34 = t4-t3; my_buf[2] = t34.count();
      std::chrono::duration<double> t55 = t5-t5start; my_buf[3] = t55.count();
      double buf_min[4];
      MPI_Reduce(my_buf,buf_min,4,MPI_DOUBLE,MPI_MIN,0,mpi_comm);
      double buf_max[4];
      MPI_Reduce(my_buf,buf_max,4,MPI_DOUBLE,MPI_MAX,0,mpi_comm);
      if (mpi_rank == 0) {
        cout << "[TIMING]: t12: " << buf_min[0] << " " << buf_max[0] <<
          " t23: " << buf_min[1] << " " << buf_max[1] <<
          " t34: " << buf_min[2] << " " << buf_max[2] <<
          " t55: " << buf_min[3] << " " << buf_max[3] << endl;
      }
#endif

      } // while (!done)

      // TODO cleanup memory, set cost, do coloring, blah blah blah...

      if (send_buf_int)    delete[] send_buf_int;
      if (send_buf_double) delete[] send_buf_double;
      if (recv_buf_int)    delete[] recv_buf_int;
      if (recv_buf_double) delete[] recv_buf_double;

      delete[] rcorp;
      delete[] wrorp;
      delete[] rcora;
      delete[] storp_i;
      delete[] rborp_i;
      delete[] nborc_i;
      delete[] nborp_i;

      delete[] send_count_double;
      delete[] send_disp_double;

      {
        int8 my_buf[2] = { (int8)vdReturnDataVec.size(), 0 };
        for (int ii = 0; ii < vdReturnDataVec.size(); ++ii) {
          my_buf[1] += vdReturnDataVec[ii].nbf;
        }
        int8 buf[2];
        MPI_Reduce(my_buf,buf,2,MPI_INT8,MPI_SUM,0,mpi_comm);
        if (mpi_rank == 0)
          cout << " > total vds built: " << buf[0] << " total nbf: " << buf[1] << endl;

        double time_buf[2];
        MPI_Reduce(my_time_buf,time_buf,2,MPI_DOUBLE,MPI_SUM,0,mpi_comm);
        if (mpi_rank == 0) cout << " > work time: " << time_buf[0]/double(mpi_size) << " wasted time: " << time_buf[1]/double(mpi_size) << endl;

        double my_minmax[6] = { my_time_buf[0], -my_time_buf[0], my_time_buf[1], -my_time_buf[1], my_time_buf[0]+my_time_buf[1], -my_time_buf[0]-my_time_buf[1] };
        double minmax[6];
        MPI_Reduce(my_minmax,minmax,6,MPI_DOUBLE,MPI_MIN,0,mpi_comm);
        if (mpi_rank == 0) {
          cout << " > work time  min/max: " << minmax[0] << " " << -minmax[1] <<
            " wasted time min/max: " << minmax[2] << " " << -minmax[3] <<
            " total time min/max: " << minmax[4] << " " << -minmax[5] << endl;
        }
      }

    }

    void rebuildFlaggedVds(const bool b_surface,const bool b_final,
        int * send_count,int * send_disp,int * recv_count,int * recv_disp) {

      // this started from an exact copy from moving...

#ifdef WITH_CHRONO
      auto t0 = std::chrono::high_resolution_clock::now();
#endif

      if ((mpi_rank == 0)&&(step%check_interval == 0)) cout << "rebuildFlaggedVds(), b_surface: " << b_surface << " final: " << b_final << endl;

      vector<int> intVec;
      vector<double> doubleVec;

      int send_count_sum_max = 0;
      int * send_buf_int = NULL;
      double * send_buf_double = NULL;

      int recv_count_sum_max = 0;
      int * recv_buf_int = NULL;
      double * recv_buf_double = NULL;

      assert(rbi_g.size() == rbiMap.size());
      const int ncv_g0 = rbi_g.size() + ncv;

      // -----------------------------------------------------------------
      // iterate to rebuild all required cvs...
      // at the end of this loop, we should have a consistent vd made up
      // of potentially old (transformed) and new (just constructed) vds.
      // -----------------------------------------------------------------

      int done = 0;
      while (done == 0) {

#ifdef WITH_CHRONO
        auto t00 = std::chrono::high_resolution_clock::now();
#endif

        // assume we are done...
        int my_done = 1;

        // store the current size. We only need to go through the NEW
        // vd's produced this iteration...
        const int vdReturnDataVec_size0 = vdReturnDataVec.size();

        // see rebuildFlaggedVds.hpp...
        _rebuildFlaggedVds(vdReturnDataVec,b_surface,
            send_count,send_disp,recv_count,recv_disp);

#ifdef WITH_CHRONO
        auto t01 = std::chrono::high_resolution_clock::now();
#endif

        // loop to decide if any of the nbrs of new vd's need to be rebuilt. A nbr of a new vd
        // should be rebuilt if:
        // 1. the nbr's delta is smaller than 2 x the maximum radius of our (shared) face nodes
        // 2. the new vd is appearing: all nbrs should be new
        // 3. the new vd is not appearing, but the nbr is new (i.e. not in our old nbr list)
        FOR_RANK send_count[rank] = 0;
        assert(intVec.empty());
        assert(doubleVec.empty());
        for (int ii = vdReturnDataVec_size0, limit = vdReturnDataVec.size(); ii < limit; ++ii) { // cycle through new vd's only...
          const int icv = vdReturnDataVec[ii].icv; assert((icv >= 0)&&(icv < ncv));
          for (int ifa = 0; ifa < vdReturnDataVec[ii].nfa; ++ifa) {
            const int icv_nbr = vdReturnDataVec[ii].cvofa[ifa];
            assert(icv_nbr != icv);
            if (icv_nbr < ncv) {
              assert(flag_cv[icv_nbr]&ACTIVE_BIT); // your nbr exists at the current time level
              // NEW or OLD, look at the nbr's delta and make sure it is large
              // enough to capture our face...
              if (delta_vd[icv_nbr]*delta_vd[icv_nbr] < 4.0*vdReturnDataVec[ii].r2_fa[ifa]) {
                delta_vd[icv_nbr] = 2.001*sqrt(vdReturnDataVec[ii].r2_fa[ifa]);
                flag_cv[icv_nbr] |= REBUILD_BIT;
                my_done = 0;
              }
              else if (!(flag_cv[icv_nbr]&NEW_BIT)) {
                // if the nbr has not been rebuilt, we need to flag it if either:
                // 1. we are appearing (then all our new nbs need to be rebuilt)...
                // 2. there is no existing fa0 between icv,icv_nbr
                if (!(flag_cv[icv]&ACTIVE0_BIT)) {
                  // we are an appearing cell, all our nbrs need to be rebuilt...
                  //assert(0); // DO WE EVER HIT THIS? - yes, in ICE for example
                  flag_cv[icv_nbr] |= REBUILD_BIT;
                  my_done = 0;
                }
                else {
                  // the moving solver uses the cvocv_i0 structure here. For stitch, this
                  // is not built because it is not needed. Instead, we allow the previous
                  // entry in vdReturnDataVec to persist and use its connections...
                  assert(ivrdocv);
                  const int ii_prev = ivrdocv[icv];
                  assert((ii_prev >= 0)&&(ii_prev < ii));
                  assert(vdReturnDataVec[ii_prev].icv == icv);
                  int ifa_prev;
                  for (ifa_prev = 0; ifa_prev < vdReturnDataVec[ii_prev].nfa; ++ifa_prev) {
                    const int icv_nbr_prev = vdReturnDataVec[ii_prev].cvofa[ifa_prev];
                    assert((icv_nbr_prev >= 0)&&(icv_nbr_prev < ncv_g0));
                    if (icv_nbr_prev == icv_nbr)
                      break;
                  }
                  if (ifa_prev == vdReturnDataVec[ii_prev].nfa) {
                    // no match found, so rebuild the nbr
                    flag_cv[icv_nbr] |= REBUILD_BIT;
                    my_done = 0;
                  }
                }
              }
            }
            else {
              // we want the off-rank behavior to exacly match the on-rank behavior
              // above. We always send one double (the fa r2) to the cv on the other rank...
              assert(icv_nbr-ncv < rbi_g.size());
              int rank,bits,index;
              BitUtils::unpackRankBitsIndex(rank,bits,index,rbi_g[icv_nbr-ncv]);
              ++send_count[rank]; // we need to send just send the index - this guy has to be rebuilt
              doubleVec.push_back(vdReturnDataVec[ii].r2_fa[ifa]);
              intVec.push_back(rank);
              // for the index, use -ve 1 indexing to indicate when the nbr absolutely
              // must be rebuilt: because we are either appearing, or this nbr is new for
              // us (i.e. NOT in our cvocv_i0/v0)...
              if (icv_nbr < ncv_g0) {
                // this nbr was in the old ghost range. Set the flag to keep it.
                //ghost_flag[icv_nbr-ncv] = true;
                if (!(flag_cv[icv]&ACTIVE0_BIT)) {
                  // we are an appearing cell, all our nbrs need to be rebuilt...
                  intVec.push_back(-index-1);
                }
                else {
                  // see notes above on moving solver use of cvocv_i0.
                  assert(ivrdocv);
                  const int ii_prev = ivrdocv[icv];
                  assert((ii_prev >= 0)&&(ii_prev < ii));
                  assert(vdReturnDataVec[ii_prev].icv == icv);
                  int ifa_prev;
                  for (ifa_prev = 0; ifa_prev < vdReturnDataVec[ii_prev].nfa; ++ifa_prev) {
                    const int icv_nbr_prev = vdReturnDataVec[ii_prev].cvofa[ifa_prev];
                    assert((icv_nbr_prev >= 0)&&(icv_nbr_prev < ncv_g0));
                    if (icv_nbr_prev == icv_nbr)
                      break;
                  }
                  if (ifa_prev == vdReturnDataVec[ii_prev].nfa) {
                    // no match found, so rebuild the nbr
                    // no match found, so rebuild the nbr
                    intVec.push_back(-index-1);
                  }
                  else {
                    // we found a match. Just send the index...
                    intVec.push_back(index);
                  }
                }
              }
              else {
                // this is a new ghost nbr. It could not possibly be in our old nbr
                // list. No matter if we are appearing or not, it needs to
                // be new...
                intVec.push_back(-index-1);
              }
            }
          }
        }
        assert(intVec.size() == doubleVec.size()*2);

        send_disp[0] = 0;
        for (int rank = 1; rank < mpi_size; ++rank)
          send_disp[rank] = send_disp[rank-1] + send_count[rank-1];
        const int send_count_sum = send_disp[mpi_size-1] + send_count[mpi_size-1];
        assert(send_count_sum == doubleVec.size());

        if (send_count_sum > send_count_sum_max) {
          send_count_sum_max = send_count_sum;
          if (send_buf_int != NULL) delete[] send_buf_int;
          send_buf_int = new int[send_count_sum_max];
          if (send_buf_double != NULL) delete[] send_buf_double;
          send_buf_double = new double[send_count_sum_max];
        }

        for (int isend = 0; isend < send_count_sum; ++isend) {
          const int rank = intVec[isend*2];
          assert((rank >= 0)&&(rank < mpi_size));
          send_buf_int[send_disp[rank]] = intVec[isend*2+1]; // may be positive or -1 indexed
          send_buf_double[send_disp[rank]] = doubleVec[isend];
          ++send_disp[rank];
        }
        intVec.clear();
        doubleVec.clear();

        // reset the disp...

        send_disp[0] = 0;
        for (int rank = 1; rank < mpi_size; ++rank)
          send_disp[rank] = send_disp[rank-1] + send_count[rank-1];

        // set up recv side and send...

        MPI_Alltoall(send_count,1,MPI_INT,recv_count,1,MPI_INT,mpi_comm);

        recv_disp[0] = 0;
        for (int rank = 1; rank < mpi_size; ++rank)
          recv_disp[rank] = recv_count[rank-1] + recv_disp[rank-1];

        const int recv_count_sum = recv_disp[mpi_size-1] + recv_count[mpi_size-1];
        if (recv_count_sum > recv_count_sum_max) {
          recv_count_sum_max = recv_count_sum;
          if (recv_buf_int != NULL) delete[] recv_buf_int;
          recv_buf_int = new int[recv_count_sum_max];
          if (recv_buf_double != NULL) delete[] recv_buf_double;
          recv_buf_double = new double[recv_count_sum_max];
        }

        MPI_Alltoallv(send_buf_int,send_count,send_disp,MPI_INT,
            recv_buf_int,recv_count,recv_disp,MPI_INT,
            mpi_comm);

        MPI_Alltoallv(send_buf_double,send_count,send_disp,MPI_DOUBLE,
            recv_buf_double,recv_count,recv_disp,MPI_DOUBLE,
            mpi_comm);

        // process...

        for (int irecv = 0; irecv < recv_count_sum; ++irecv) {
          if (recv_buf_int[irecv] >= 0) {
            // a positive index means we only do the r2 comparison. The sending nbr was not
            // appearing and we were already in its nbr list...
            const int icv = recv_buf_int[irecv];
            assert((icv >= 0)&&(icv < ncv));
            assert(flag_cv[icv]&ACTIVE_BIT);
            if (delta_vd[icv]*delta_vd[icv] < 4.0*recv_buf_double[irecv]) {
              //cout << "XXXXXXXXXXXXXXXXXXXXX got small delta_vd. increasing by: " << 2.001*sqrt(recv_buf_double[irecv])/delta_vd[icv] << endl;
              delta_vd[icv] = 2.001*sqrt(recv_buf_double[irecv]);
              flag_cv[icv] |= REBUILD_BIT;
              my_done = 0;
            }
          }
          else {
            // a negative index means automatic rebuild if we are not new: still check r2 first...
            const int icv = -recv_buf_int[irecv]-1;
            assert((icv >= 0)&&(icv < ncv));
            assert(flag_cv[icv]&ACTIVE_BIT);
            if (delta_vd[icv]*delta_vd[icv] < 4.0*recv_buf_double[irecv]) {
              // got small delta, so increase it...
              delta_vd[icv] = 2.001*sqrt(recv_buf_double[irecv]);
              flag_cv[icv] |= REBUILD_BIT;
              my_done = 0;
            }
            else if (!(flag_cv[icv]&NEW_BIT)) {
              // automatically rebuild if not new. The send side did the
              // checking and set the -ve 1 indexing of icv...
              flag_cv[icv] |= REBUILD_BIT;
              my_done = 0;
            }
          }
        }

        MPI_Allreduce(&my_done,&done,1,MPI_INT,MPI_MIN,mpi_comm);

        // make sure rebuild wasn't set in a cv touching the boundary when
        // skip_surface_rebuild has been set...
        assert(b_surface);

#ifdef WITH_CHRONO
        auto t02 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed01 = t01-t00;
        std::chrono::duration<double> elapsed02 = t02-t01;
        my_buf[1] += elapsed01.count();
        my_buf[2] += elapsed02.count();
#endif

      } // while (done == 0)...

#ifdef WITH_CHRONO
      auto t1 = std::chrono::high_resolution_clock::now();
#endif

      delete[] send_buf_double; send_buf_double = NULL;
      delete[] recv_buf_double; recv_buf_double = NULL;

      delete[] send_buf_int; send_buf_int = NULL;
      delete[] recv_buf_int; recv_buf_int = NULL;

      // ===========================================================
      // the vdReturnDataVec is complete, so complete grid...
      // ===========================================================

      // at this point, there may be duplicate icv entries in vdReturnDataVec. Cycle through this
      // in reverse so we give the last ones (i.e. built last) priority and do not consider earlier ones...

      if (ivrdocv == NULL) ivrdocv = new int[ncv]; // first orphan index of cv
      FOR_ICV ivrdocv[icv] = -1;

      int my_buf[2] = { 0, 0 }; // count 0: cvs with orphans, and 1: orphans
      bool duplicates = false;
      for (int ii = vdReturnDataVec.size()-1; ii >= 0; --ii) {
        const int icv = vdReturnDataVec[ii].icv; assert((icv >= 0)&&(icv < ncv));
        assert(flag_cv[icv] & ACTIVE_BIT); // must be active
        assert(!(flag_cv[icv] & REBUILD_BIT));
        if (ivrdocv[icv] == -1) {
          ivrdocv[icv] = ii;
          if (vdReturnDataVec[ii].ngr > 1) {
            ++my_buf[0];
            my_buf[1] += vdReturnDataVec[ii].ngr-1;
          }
          assert(vdReturnDataVec[ii].orphan_chunk_data == -1); // should be no orphans referenced
        }
        else {
          // this is an earlier duplicate, so clear and set icv to -1...
          duplicates = true;
          vdReturnDataVec[ii].clear();
          assert(vdReturnDataVec[ii].icv == -1); // clear sets -1
        }
      }

      // duplicates is an entirely local thing: if present,
      // roll forward, compressing vdReturnDataVec. Note that
      // we have been careful NOT to write a proper destructor for
      // VdReturnDataVec which allows simple POD data copies here including
      // the pointers...
      if (duplicates) {
        int ii_new = 0;
        for (int ii = 0; ii < vdReturnDataVec.size(); ++ii) {
          const int icv = vdReturnDataVec[ii].icv;
          if (icv >= 0) {
            if (ii_new != ii) {
              vdReturnDataVec[ii_new] = vdReturnDataVec[ii];
              ivrdocv[icv] = ii_new;
            }
            ++ii_new;
          }
        }
        assert(ii_new < vdReturnDataVec.size());
        vdReturnDataVec.resize(ii_new);
      }

      assert(ocdVec.empty());

      int * cv_check = new int[ncv];
      FOR_ICV cv_check[icv] = 0;

      for (int ii = 0; ii < vdReturnDataVec.size(); ++ii) {
        const int icv = vdReturnDataVec[ii].icv;
        assert(icv >= 0);
        cv_check[icv] = 1;
        double x_cv_vd[3];
        vdReturnDataVec[ii].calcVolumeAndCentroidForGroup(vol_cv[icv],x_cv_vd,0);
        FOR_I3 x_cv[icv][i] = vol_cv[icv]*x_cv_vd[i]; // leave x_cv in volume-weighted local coordinates relative to x_vd
      }

      // check: TODO: get rid of this eventually...
      FOR_ICV FOR_I3 assert(x_cv[icv][i] == x_cv[icv][i]);

      // if ANY orphans are present, they will modify the centroids so we ALL have to handle them...

      int buf[2];
      MPI_Allreduce(my_buf,buf,2,MPI_INT,MPI_SUM,mpi_comm);

      if (buf[0] == 0) {

        assert(buf[1] == 0);
        if (mpi_rank == 0) cout << " > got 0 orphans" << endl;

      }
      else {

        const int nor = my_buf[1];
        int my_nor_minmax[2] = { nor, -nor };
        int nor_minmax[2];
        MPI_Reduce(my_nor_minmax,nor_minmax,2,MPI_INT,MPI_MIN,0,mpi_comm);

        if (mpi_rank == 0) {
          cout << " > got " << buf[0] << " cvs with " <<
            buf[1] << " orphans: orphans/cv avg: " << double(buf[1])/double(buf[0]) << endl;
          cout << " > orphans per rank [min,avg,max]: " <<
            nor_minmax[0] << " " <<
            double(buf[1])/double(mpi_size) << " " <<
            -nor_minmax[1] << endl;
        }

        // we need to link the orphans so they reach out through nbring orphans into their
        // set of primary cvs...

        int (*cgoor)[2] = new int[nor][2]; // cv-and-group-of-orphan

        // connect the orphans...
        // step 1: build the varios csr structures required by the orphan graph...

        int * oroii_i = new int[vdReturnDataVec.size()+1];
        oroii_i[0] = 0;
        set<int> intSet;
        set<pair<int,int> > intPairSet;
        int * nboor_i = new int[nor+1]; // nbr-of-orphan: index
        vector<int> nboor_v;            // this value should eventually hold the ior_nbr, but holds icv at first
        vector<int> stonboor_i; // from a particular nbr (orphan face) this indexes into the ist's (if any)
        vector<int> stonboor_v; // "flattened" ist values
        nboor_i[0] = 0;
        stonboor_i.push_back(0);
        int ior_check = 0;
        for (int ii = 0; ii < vdReturnDataVec.size(); ++ii) {
          const int icv = vdReturnDataVec[ii].icv;
          // only work in this vd if he is active (icv >= 0), AND he has orphans...
          if ((icv >= 0)&&(vdReturnDataVec[ii].ngr > 1)) {
            for (int igr = 1; igr < vdReturnDataVec[ii].ngr; ++igr) {
              // store the icv and group info for this orphan
              cgoor[ior_check][0] = icv;
              cgoor[ior_check][1] = igr;
              ++ior_check;
              assert(intSet.empty());
              assert(intPairSet.empty());
              for (int ied = vdReturnDataVec[ii].edogr_i[igr]; ied != vdReturnDataVec[ii].edogr_i[igr+1]; ++ied) {
                FOR_I2 {
                  const int ifa = vdReturnDataVec[ii].faoed[ied][i];
                  if (ifa < 0) {
                    // a face less than zero references an internal face cut by
                    // a cv in cvofa. At this point we build 2 sets:
                    // intSet: the set of unique faces (associated with unique icv_nbr's)
                    // intPairSet: face/ist pairs that we come across during the edge checking...
                    intSet.insert(-ifa-1); // 0..nfa-1
                    if (vdReturnDataVec[ii].faoed[ied][1-i] >= 0) {
                      // this is an edge of a face touching a surface tri...
                      const int ifa_opp = vdReturnDataVec[ii].faoed[ied][1-i];
                      assert(ifa_opp < vdReturnDataVec[ii].nbf);
                      // this face references a particular ist and part_bits...
                      const int ist = vdReturnDataVec[ii].spbobf[ifa_opp][0];
                      const int ipart_bits = vdReturnDataVec[ii].spbobf[ifa_opp][1];
                      //const int bits = (ipart_bits&MASK_6BITS);
                      const int ipart = (ipart_bits>>6);
                      assert((ipart >= 0)&&(ipart < partVec.size()));
                      intPairSet.insert(pair<int,int>(-ifa-1,getFlattenedSt(ipart,ist)));
                    }
                  }
                }
              }
              // at this point we have all the faces in intSet, and any face-ist pairs in intPairVec...
              set<pair<int,int> >::iterator it2 = intPairSet.begin();
              for (set<int>::iterator it = intSet.begin(); it != intSet.end(); ++it) {
                const int ifa = *it;
                assert((ifa >= 0)&&(ifa < vdReturnDataVec[ii].nfa));
                nboor_v.push_back(vdReturnDataVec[ii].cvofa[ifa]); // for now, push in the local cv value that cut this face
                while ((it2 != intPairSet.end())&&(it2->first == ifa)) {
                  stonboor_v.push_back(it2->second);
                  ++it2;
                }
                stonboor_i.push_back(stonboor_v.size());
              }
              nboor_i[ior_check] = nboor_v.size();
              assert(it2 == intPairSet.end());
              intSet.clear();
              intPairSet.clear();
            }
          }
          oroii_i[ii+1] = ior_check;
        }
        assert(ior_check == nor);

        // Now implement the matching rules to set the nboor_v's to actual orphan nbr
        // index when the connection is with an orphan, or a -1-indexed icv when it is
        // a main cv connection.
        //
        // Use the idea that a face that has any (or all) of the same ist's must be
        // a valid connection. If no valid connection can be found with matching ist's,
        // then there must be one and only one face that connects these vd groups.

        //int debug_rank = getIntParam("DEBUG_RANK",-1);
        //int debug_ior = getIntParam("DEBUG_IOR",-1);

        double my_unmatched_n2_max = 0.0;
        set<pair<int,int> > orphanLinkSet;
        set<int>* cvoor_set = new set<int>[nor];

        int* send_count2 = new int[mpi_size]; // for return message
        FOR_RANK {
          send_count[rank] = 0;
          send_count2[rank] = 0;
        }
        for (int ior = 0; ior < nor; ++ior) {
          const int icv = cgoor[ior][0];
          // loop on our face nbrs...
          for (int noo = nboor_i[ior]; noo != nboor_i[ior+1]; ++noo) {
            // use a -ve indexing to handle a particular connection just once...
            const int icv_nbr = nboor_v[noo];
            assert(icv_nbr >= 0); // for now?
            if (icv_nbr < ncv) {
              // this is a local nbr. Loop through the orphans of the nbring vd, if any,
              // looking for a match...
              const int ii_nbr = ivrdocv[icv_nbr];
              assert((ii_nbr >= 0)&&(ii_nbr < vdReturnDataVec.size()));
              assert(vdReturnDataVec[ii_nbr].icv == icv_nbr);
              // if/when the above fails: note that for partial rebuilds the vd for this icv may not
              // be available, in which case it will not have orphans, and can be
              // used simply as a primary cv contributing to the orphan cutting...
              // check the oroii_i index structure...
              assert(vdReturnDataVec[ii_nbr].ngr-1 == oroii_i[ii_nbr+1]-oroii_i[ii_nbr]);
              // to match an orphan, or a main vd for that matter, to this orphan nbr,
              // we use the following rules (in priority order):
              // 1. if the connection shares ANY (and normally ALL) of the surface tris
              // associated with our face, then it is the correct connection, and is unique
              // 2. if the connection is the ONLY connection to another vd, then
              // that is the unique connection.
              // 3. if there is no connection, the orphan must be tiny, and can be discarded
              //
              // for fast ist matching, put the tris (if any) associated with our face into intSet...
              assert(intSet.empty());
              for (int sonoo = stonboor_i[noo]; sonoo != stonboor_i[noo+1]; ++sonoo)
                intSet.insert(stonboor_v[sonoo]);
              // we now have any tris this orphan face touches in intSet...
              bool got_any_match = false;
              bool got_internal_match = false;
              if (vdReturnDataVec[ii_nbr].ngr > 1) {
                // this nbr has orphans, so check the orphan nbrs first...
                for (int ior_nbr = oroii_i[ii_nbr]; ior_nbr != oroii_i[ii_nbr+1]; ++ior_nbr) {
                  // we already know that ior has a connection with the ior_nbr's ii. loop
                  // on ior_nbr's connections and see if there is a connection back to icv (us)...
                  for (int noo_nbr = nboor_i[ior_nbr]; noo_nbr != nboor_i[ior_nbr+1]; ++noo_nbr) {
                    // use a -ve indexing to handle a particular connection just once...
                    const int icv_nbr_nbr = nboor_v[noo_nbr];
                    if (icv_nbr_nbr == icv) {
                      // this orphan has a connection back. check if any of the st's associated
                      // with this nbr's face match the set of st's in intSet associated with our
                      // face...
                      if (intSet.empty() && (stonboor_i[noo_nbr+1]-stonboor_i[noo_nbr] == 0)) {
                        // this is a face match all internal edges. There should be just 1...
                        assert(!got_internal_match);
                        assert(!got_any_match);
                        orphanLinkSet.insert(pair<int,int>(ior,ior_nbr));
                        got_internal_match = got_any_match = true;
                      }
                      else {
                        for (int sonoo_nbr = stonboor_i[noo_nbr]; sonoo_nbr != stonboor_i[noo_nbr+1]; ++sonoo_nbr) {
                          if (intSet.find(stonboor_v[sonoo_nbr]) != intSet.end()) {
                            orphanLinkSet.insert(pair<int,int>(ior,ior_nbr));
                            got_any_match = true;
                            break;
                          }
                        }
                      }
                    }
                  }
                }
              }
              // and finish by checking the main vd at ii_nbr to either determine the match if
              // not yet found, or to ensure there are no conflicts with the matching rules...
              assert(vdReturnDataVec[ii_nbr].icv == icv_nbr);
              bool got_icv = false;
              bool got_ist = false;
              bool got_ist_match = false;
              for (int ied = 0; ied < vdReturnDataVec[ii_nbr].edogr_i[1]; ++ied) {
                FOR_I2 {
                  const int ifa = vdReturnDataVec[ii_nbr].faoed[ied][i];
                  if (ifa < 0) {
                    assert(-ifa-1 < vdReturnDataVec[ii_nbr].nfa);
                    const int icv_nbr_nbr = vdReturnDataVec[ii_nbr].cvofa[-ifa-1];
                    if (icv_nbr_nbr == icv) {
                      // this is a potential connection back to the orphan...
                      got_icv = true;
                      if (vdReturnDataVec[ii_nbr].faoed[ied][1-i] >= 0) {
                        got_ist = true;
                        // this is an edge of the face touching a surface tri...
                        const int ifa_opp = vdReturnDataVec[ii_nbr].faoed[ied][1-i];
                        assert(ifa_opp < vdReturnDataVec[ii_nbr].nbf);
                        // this face references a particular ist and part_bits...
                        const int ist = vdReturnDataVec[ii_nbr].spbobf[ifa_opp][0];
                        const int ipart_bits = vdReturnDataVec[ii_nbr].spbobf[ifa_opp][1];
                        //const int bits = (ipart_bits&MASK_6BITS);
                        const int ipart = (ipart_bits>>6);
                        assert((ipart >= 0)&&(ipart < partVec.size()));
                        if (intSet.find(getFlattenedSt(ipart,ist)) != intSet.end())
                          got_ist_match = true;
                      }
                    }
                  }
                }
              }
              if (got_icv) {
                if (intSet.empty() && (!got_ist)) {
                  // internal face connection with parent icv_nbr. It should be the ONLY connection...
                  assert(!got_internal_match);
                  assert(!got_any_match);
                  cvoor_set[ior].insert(icv_nbr);
                  got_internal_match = got_any_match = true;
                }
                else if (got_ist_match) {
                  assert(!got_internal_match);
                  cvoor_set[ior].insert(icv_nbr);
                  got_any_match = true;
                }
              }
              if (!got_any_match) {
                // no matches anywhere: do small face check...
                double n_fa[3] = { 0.0, 0.0, 0.0 };
                // check that this face is small...
                const int ii = ivrdocv[icv];
                assert(vdReturnDataVec[ii].icv == icv);
                const int igr = cgoor[ior][1];
                assert((igr >= 1)&&(igr < vdReturnDataVec[ii].ngr));
                double * x0 = NULL;
                for (int ied = vdReturnDataVec[ii].edogr_i[igr]; ied < vdReturnDataVec[ii].edogr_i[igr+1]; ++ied) {
                  FOR_I2 {
                    const int ifa = vdReturnDataVec[ii].faoed[ied][i];
                    if (ifa < 0) {
                      assert(-ifa-1 < vdReturnDataVec[ii].nfa);
                      if (vdReturnDataVec[ii].cvofa[-ifa-1] == icv_nbr) {
                        if (x0 == NULL) {
                          x0 = vdReturnDataVec[ii].x_no[vdReturnDataVec[ii].nooed[ied][i]];
                        }
                        else {
                          double * x1 = vdReturnDataVec[ii].x_no[vdReturnDataVec[ii].nooed[ied][i]];
                          double * x2 = vdReturnDataVec[ii].x_no[vdReturnDataVec[ii].nooed[ied][1-i]];
                          double this_n[3] = TRI_NORMAL_2(x0,x1,x2);
                          FOR_I3 n_fa[i] += this_n[i];
                        }
                      }
                    }
                  }
                }
                my_unmatched_n2_max = max(my_unmatched_n2_max,DOT_PRODUCT(n_fa,n_fa));
              }
              intSet.clear();
            }
            else {

              // parallel...
              assert(icv_nbr >= ncv);
              int rank,bits,index;
              BitUtils::unpackRankBitsIndex(rank,bits,index,rbi_g[icv_nbr-ncv]);
              assert((bits != 0)||(rank != mpi_rank));
              // need to send over everything the nbr needs to check for a connection...
              // index on rank, rank_bits here, icv here, ior here, nst, ist0, ist1...
              send_count[rank] += 5 + stonboor_i[noo+1]-stonboor_i[noo];
              // also, prepare the send_count2 to get back the result as 3 ints...
              // first in is the status:
              // if status == 0, this is an orphan connection, rank_bits and or index returned
              // if status == 1, this is a parent cv connection, rank_bits and cv index returned
              // if status == -1: no connection found
              send_count2[rank] += 3;
            }
          }
        }

        // now the parallel...

        send_disp[0] = 0;
        for (int rank = 1; rank < mpi_size; ++rank)
          send_disp[rank] = send_disp[rank-1] + send_count[rank-1];
        int send_count_sum = send_disp[mpi_size-1] + send_count[mpi_size-1];

        assert(send_buf_int == NULL);
        send_buf_int = new int[send_count_sum];
        for (int ior = 0; ior < nor; ++ior) {
          const int icv = cgoor[ior][0];
          // loop on our face nbrs...
          for (int noo = nboor_i[ior]; noo != nboor_i[ior+1]; ++noo) {
            // use a -ve indexing to handle a particular connection just once...
            const int icv_nbr = nboor_v[noo];
            assert(icv_nbr >= 0); // for now?
            if (icv_nbr >= ncv) {
              // parallel...
              assert(icv_nbr >= ncv);
              int rank,bits,index;
              BitUtils::unpackRankBitsIndex(rank,bits,index,rbi_g[icv_nbr-ncv]);
              assert((bits != 0)||(rank != mpi_rank));
              // need to send over everything the nbr needs to check for a connection...
              // index on rank, rank_bits here, icv here, ior here, nst, ist0, ist1...
              send_buf_int[send_disp[rank]++] = index;
              send_buf_int[send_disp[rank]++] = BitUtils::packRankBits(mpi_rank,BitUtils::flipPeriodicBits(bits));
              send_buf_int[send_disp[rank]++] = icv;
              send_buf_int[send_disp[rank]++] = ior; // send instead of group
              send_buf_int[send_disp[rank]++] = stonboor_i[noo+1]-stonboor_i[noo];
              for (int sonoo = stonboor_i[noo]; sonoo != stonboor_i[noo+1]; ++sonoo)
                send_buf_int[send_disp[rank]++] = stonboor_v[sonoo];
            }
          }
        }

        // rewind...

        send_disp[0] = 0;
        for (int rank = 1; rank < mpi_size; ++rank)
          send_disp[rank] = send_disp[rank-1] + send_count[rank-1];

        MPI_Alltoall(send_count,1,MPI_INT,recv_count,1,MPI_INT,mpi_comm);

        recv_disp[0] = 0;
        for (int rank = 1; rank < mpi_size; ++rank)
          recv_disp[rank] = recv_count[rank-1] + recv_disp[rank-1];
        int recv_count_sum = recv_disp[mpi_size-1] + recv_count[mpi_size-1];

        assert(recv_buf_int == NULL);
        recv_buf_int = new int[recv_count_sum];
        MPI_Alltoallv(send_buf_int,send_count,send_disp,MPI_INT,
            recv_buf_int,recv_count,recv_disp,MPI_INT,mpi_comm);
        //delete[] send_buf_int; send_buf_int = NULL;

        assert(recv_buf_double == NULL); recv_buf_double = new double[recv_count_sum]; // should be large enough

        // need to expand your orphans to include off rank orphans that have nbr connections on this rank...
        int nor_g = nor;
        map<const uint8,int> rboMap; // rank-bits-orphan
        int irecv = 0;
        int irecv2 = 0;
        FOR_RANK recv_count[rank] = 0; // recount for the response
        while (irecv < recv_count_sum) {
          const int icv_nbr = recv_buf_int[irecv++]; // on rank is the nbr of orphan
          assert((icv_nbr >= 0)&&(icv_nbr < ncv));
          const int rank_bits = recv_buf_int[irecv++];
          const int index = recv_buf_int[irecv++];
          const int ior_there = recv_buf_int[irecv++]; // for send back
          const int nsonoo = recv_buf_int[irecv++];
          int rank,bits;
          BitUtils::unpackRankBits(rank,bits,rank_bits);
          const int flip_bits = BitUtils::flipPeriodicBits(bits);
          recv_count[rank] += 3; // response size
          // connection associated with a small face on the send side
          map<const uint8,int>::iterator it = rbiMap.find(BitUtils::packRankBitsIndex(rank,bits,index));
          if (it == rbiMap.end()) {
            // rbi not found in our current ghost data. This must mean the face send must be small...
            FOR_I3 recv_buf_double[irecv2+i] = HUGE_VAL;
            recv_buf_int[irecv2++] = -1;
            recv_buf_int[irecv2++] = 0;
            recv_buf_int[irecv2++] = 0;
            irecv += nsonoo;
          }
          else {
            // found connection, so group ghost orphans...
            const int icv = it->second; assert(icv >= ncv);
            const uint8 rbo = BitUtils::packRankBitsIndex(rank,bits,ior_there);
            map<const uint8,int>::iterator it = rboMap.find(rbo);
            int ior;
            if (it == rboMap.end()) {
              rboMap[rbo] = ior = nor_g++;
            }
            else {
              ior = it->second;
            }
            assert((ior >= nor)&&(ior < nor_g));
            const int ii_nbr = ivrdocv[icv_nbr];
            assert((ii_nbr >= 0)&&(ii_nbr < vdReturnDataVec.size()));
            assert(vdReturnDataVec[ii_nbr].icv == icv_nbr);
            // if/when the above fails: note that for partial rebuilds the vd for this icv may not
            // be available, in which case it will not have orphans, and can be
            // used simply as a primary cv contributing to the orphan cutting...
            // check the oroii_i index structure...
            assert(vdReturnDataVec[ii_nbr].ngr-1 == oroii_i[ii_nbr+1]-oroii_i[ii_nbr]);

            // for fast ist matching, put the tris (if any) associated with our face into intSet...
            assert(intSet.empty());
            for (int sonoo = 0; sonoo < nsonoo; ++sonoo)
              intSet.insert(recv_buf_int[irecv++]);

            int ior_nbr_match = -1;
            bool got_any_match = false;
            bool got_internal_match = false;
            if (vdReturnDataVec[ii_nbr].ngr > 1) {
              // this nbr has orphans, so check the orphan nbrs first...
              for (int ior_nbr = oroii_i[ii_nbr]; ior_nbr != oroii_i[ii_nbr+1]; ++ior_nbr) {
                // we already know that ior has a connection with the ior_nbr's ii. loop
                // on ior_nbr's connections and see if there is a connection back to icv (us)...
                for (int noo_nbr = nboor_i[ior_nbr]; noo_nbr != nboor_i[ior_nbr+1]; ++noo_nbr) {
                  // use a -ve indexing to handle a particular connection just once...
                  const int icv_nbr_nbr = nboor_v[noo_nbr];
                  if (icv_nbr_nbr == icv) {
                    // this orphan has a connection back. check if any of the st's associated
                    // with this nbr's face match the set of st's in intSet associated with our
                    // face...
                    if (intSet.empty() && (stonboor_i[noo_nbr+1]-stonboor_i[noo_nbr] == 0)) {
                      // this is a face match all internal edges. There should be just 1...
                      assert(!got_internal_match);
                      assert(!got_any_match);
                      orphanLinkSet.insert(pair<int,int>(ior_nbr,ior));
                      ior_nbr_match = ior_nbr;
                      got_internal_match = got_any_match = true;
                    }
                    else {
                      for (int sonoo_nbr = stonboor_i[noo_nbr]; sonoo_nbr != stonboor_i[noo_nbr+1]; ++sonoo_nbr) {
                        if (intSet.find(stonboor_v[sonoo_nbr]) != intSet.end()) {
                          orphanLinkSet.insert(pair<int,int>(ior_nbr,ior));
                          ior_nbr_match = ior_nbr;

                          got_any_match = true;
                          break;
                        }
                      }
                    }
                  }
                }
              }
            }
            // and finish by checking the main vd at ii_nbr to either determine the match if
            // not yet found, or to ensure there are no conflicts with the matching rules...
            assert(vdReturnDataVec[ii_nbr].icv == icv_nbr);
            bool got_icv = false;
            bool got_ist = false;
            bool got_ist_match = false;
            for (int ied = 0; ied < vdReturnDataVec[ii_nbr].edogr_i[1]; ++ied) {
              FOR_I2 {
                const int ifa = vdReturnDataVec[ii_nbr].faoed[ied][i];
                if (ifa < 0) {
                  assert(-ifa-1 < vdReturnDataVec[ii_nbr].nfa);
                  const int icv_nbr_nbr = vdReturnDataVec[ii_nbr].cvofa[-ifa-1];
                  if (icv_nbr_nbr == icv) {
                    // this is a potential connection back to the orphan...
                    got_icv = true;
                    if (vdReturnDataVec[ii_nbr].faoed[ied][1-i] >= 0) {
                      got_ist = true;
                      // this is an edge of the face touching a surface tri...
                      const int ifa_opp = vdReturnDataVec[ii_nbr].faoed[ied][1-i];
                      assert(ifa_opp < vdReturnDataVec[ii_nbr].nbf);
                      // this face references a particular ist and part_bits...
                      const int ist = vdReturnDataVec[ii_nbr].spbobf[ifa_opp][0];
                      const int ipart_bits = vdReturnDataVec[ii_nbr].spbobf[ifa_opp][1];
                      //const int bits = (ipart_bits&MASK_6BITS);
                      const int ipart = (ipart_bits>>6);
                      assert((ipart >= 0)&&(ipart < partVec.size()));
                      if (intSet.find(getFlattenedSt(ipart,ist)) != intSet.end())
                        got_ist_match = true;
                    }
                  }
                }
              }
            }
            bool b_irecv2_set = false;
            if (got_icv) {
              if (intSet.empty() && (!got_ist)) {
                // internal face connection with parent icv_nbr. It should be the ONLY connection...
                assert(!got_internal_match);
                assert(!got_any_match);
                FOR_I3 recv_buf_double[irecv2+i] = x_vd[icv_nbr][i];
                if (bits) PeriodicData::periodicTranslate(recv_buf_double+irecv2,1,flip_bits);
                recv_buf_int[irecv2++] = 1;
                recv_buf_int[irecv2++] = BitUtils::packRankBits(mpi_rank,flip_bits);
                recv_buf_int[irecv2++] = icv_nbr;
                got_internal_match = got_any_match = true;
                b_irecv2_set = true;
              }
              else if (got_ist_match) {
                assert(!got_internal_match);
                FOR_I3 recv_buf_double[irecv2+i] = x_vd[icv_nbr][i];
                if (bits) PeriodicData::periodicTranslate(recv_buf_double+irecv2,1,flip_bits);
                recv_buf_int[irecv2++] = 1;
                recv_buf_int[irecv2++] = BitUtils::packRankBits(mpi_rank,flip_bits);
                recv_buf_int[irecv2++] = icv_nbr;
                got_any_match = true;
                b_irecv2_set = true;
              }
            }
            if (!b_irecv2_set) {
              if (got_any_match) {
                // any match we already found should trigger the response that
                // we connected one or more ghosts. Here we just send back the
                // last connection if there were multiple, so the checking is
                // not exaustive...
                FOR_I3 recv_buf_double[irecv2+i] = HUGE_VAL;
                recv_buf_int[irecv2++] = 0; // status 0
                recv_buf_int[irecv2++] = BitUtils::packRankBits(mpi_rank,flip_bits);
                assert(ior_nbr_match >= 0); // should have an ior_nbr_match in here...
                recv_buf_int[irecv2++] = ior_nbr_match;
              }
              else {
                // no matches anywhere: this will trigger a small face check on the send side...
                FOR_I3 recv_buf_double[irecv2+i] = HUGE_VAL;
                recv_buf_int[irecv2++] = -1; // status
                recv_buf_int[irecv2++] = 0; // not used
                recv_buf_int[irecv2++] = 0; // not used
              }
            }
            intSet.clear();
          }
        }
        assert(irecv == recv_count_sum);
        uint8 *rbo_g = new uint8[nor_g-nor];
        for (map<const uint8,int>::iterator it = rboMap.begin(); it != rboMap.end(); ++it) {
          const int ior = it->second; assert((ior >= nor)&&(ior < nor_g));
          rbo_g[ior-nor] = it->first;
        }

        // at this point, we have been through the recv data, repacking the response
        // in a structured way. It is different, so we need to rebuild the
        // disp's...

        recv_disp[0] = 0;
        for (int rank = 1; rank < mpi_size; ++rank)
          recv_disp[rank] = recv_count[rank-1] + recv_disp[rank-1];
        if (!(irecv2 == recv_disp[mpi_size-1] + recv_count[mpi_size-1]))
          cout << "ERROR: " << mpi_rank << " " << irecv2 << " " << recv_disp[mpi_size-1] + recv_count[mpi_size-1] << endl;
        assert(irecv2 == recv_disp[mpi_size-1] + recv_count[mpi_size-1]);

        // note the use of send_count2...

        // check: get rid of this eventually...
        {
          MPI_Alltoall(recv_count,1,MPI_INT,send_count,1,MPI_INT,mpi_comm);
          FOR_RANK assert(send_count[rank] == send_count2[rank]);
        }

        send_disp[0] = 0;
        for (int rank = 1; rank < mpi_size; ++rank)
          send_disp[rank] = send_disp[rank-1] + send_count2[rank-1];

        // reverse send
        MPI_Alltoallv(recv_buf_int,recv_count,recv_disp,MPI_INT,
            send_buf_int,send_count2,send_disp,MPI_INT,mpi_comm);
        delete[] recv_buf_int; recv_buf_int = NULL;
        assert(send_buf_double == NULL); send_buf_double = new double[send_disp[mpi_size-1]+send_count2[mpi_size-1]];
        MPI_Alltoallv(recv_buf_double,recv_count,recv_disp,MPI_DOUBLE,
            send_buf_double,send_count2,send_disp,MPI_DOUBLE,mpi_comm);
        delete[] recv_buf_double; recv_buf_double = NULL;

        // unpack and respond with checking. But before we do, allocate the set's associated
        // with each group.

        // finish by processing the returned send_buf_int. We could have sent just 1 int
        // back per request, but we return 3 to allow some checking for now...

        map<const uint8,int> rbiOrMap;
        vector<uint8> rbi_or_g;
        vector<double> xcv_or_g;
        for (int ior = 0; ior < nor; ++ior) {
          const int icv = cgoor[ior][0];
          for (int noo = nboor_i[ior]; noo != nboor_i[ior+1]; ++noo) {
            const int icv_nbr = nboor_v[noo];
            assert(icv_nbr >= 0); // for now?
            if (icv_nbr >= ncv) {
              int rank,bits,index_nbr;
              BitUtils::unpackRankBitsIndex(rank,bits,index_nbr,rbi_g[icv_nbr-ncv]);
              assert((bits != 0)||(rank != mpi_rank));
              // index on rank, rank_bits, icv here, nst, ist0, ist1...
              const int status    = send_buf_int[send_disp[rank]++];
              const int rank_bits = send_buf_int[send_disp[rank]++];
              const int index     = send_buf_int[send_disp[rank]++];
              assert((status >= -1)&&(status <= 1));
              if (status == 1) {
                // this is a returned icv link. It should correspond to the rbi associated with  the ghost icv above...
                assert(rbi_g[icv_nbr-ncv] == BitUtils::packRankBitsIndex(rank_bits,index));
                // add this nbr to ior's cv set w/ -1 idx...
                map<const uint8,int>::iterator it = rbiOrMap.find(rbi_g[icv_nbr-ncv]);
                if (it == rbiOrMap.end()) {
                  // this is a ghost not currently in the rbi...
                  cvoor_set[ior].insert(-rbi_or_g.size()-1);
                  rbiOrMap[rbi_g[icv_nbr-ncv]] = rbi_or_g.size();
                  rbi_or_g.push_back(rbi_g[icv_nbr-ncv]);
                  FOR_I3 xcv_or_g.push_back(send_buf_double[(send_disp[rank]-3)+i]);
                }
                else {
                  cvoor_set[ior].insert(-it->second-1);
                }
              }
              else if (status == 0) {
                // produced one or more ghost connections. By symmetry, these should already
                // be in the orphan ghost data map. We can check the one that was sent back, but
                // others may have been made that can't be checked, but this is better than
                // nothing...
                map<const uint8,int>::iterator it = rboMap.find(BitUtils::packRankBitsIndex(rank_bits,index));
                assert(it != rboMap.end()); // should exist
                // AND a connection between this ghost value and the local ior's group should also exist...
                assert(orphanLinkSet.find(pair<int,int>(ior,it->second)) != orphanLinkSet.end());
              }
              else {
                assert(status == -1);
                // no matches anywhere: do small face check...
                double n_fa[3] = { 0.0, 0.0, 0.0 };
                // check that this face is small...
                const int ii = ivrdocv[icv];
                assert(vdReturnDataVec[ii].icv == icv);
                const int igr = cgoor[ior][1];
                assert((igr >= 1)&&(igr < vdReturnDataVec[ii].ngr));
                double * x0 = NULL;
                for (int ied = vdReturnDataVec[ii].edogr_i[igr]; ied < vdReturnDataVec[ii].edogr_i[igr+1]; ++ied) {
                  FOR_I2 {
                    const int ifa = vdReturnDataVec[ii].faoed[ied][i];
                    if (ifa < 0) {
                      assert(-ifa-1 < vdReturnDataVec[ii].nfa);
                      if (vdReturnDataVec[ii].cvofa[-ifa-1] == icv_nbr) {
                        if (x0 == NULL) {
                          x0 = vdReturnDataVec[ii].x_no[vdReturnDataVec[ii].nooed[ied][i]];
                        }
                        else {
                          double * x1 = vdReturnDataVec[ii].x_no[vdReturnDataVec[ii].nooed[ied][i]];
                          double * x2 = vdReturnDataVec[ii].x_no[vdReturnDataVec[ii].nooed[ied][1-i]];
                          double this_n[3] = TRI_NORMAL_2(x0,x1,x2);
                          FOR_I3 n_fa[i] += this_n[i];
                        }
                      }
                    }
                  }
                }
                my_unmatched_n2_max = max(my_unmatched_n2_max,DOT_PRODUCT(n_fa,n_fa));
              }
            }
          }
        }
        rboMap.clear();
        delete[] send_buf_int; send_buf_int = NULL;
        delete[] send_buf_double; send_buf_double = NULL;

        double unmatched_n2_max;
        MPI_Reduce(&my_unmatched_n2_max,&unmatched_n2_max,1,MPI_DOUBLE,MPI_MAX,0,mpi_comm);
        if (mpi_rank == 0)
          cout << " > unmatched_n2_max: " << unmatched_n2_max << endl;

        int* send_disp2 = new int[mpi_size];
        int* recv_count2 = new int[mpi_size];
        int* recv_disp2 = new int[mpi_size];

        // new orphan set-building...
        int done = 0;
        set<int>* potential_cvoor_set = new set<int>[nor];
        set<int>* final_cvoor_set = new set<int>[nor];
        while (done == 0) {

          // assume we are done...
          int my_done = 1;

          FOR_RANK {
            send_count[rank] = 0;
            send_count2[rank] = 0;
          }

          // push all icv's from current cvoor set into potential nbrs of our orphan nbrs...
          set<pair<int,int> >::iterator it = orphanLinkSet.begin();
          for (int ior = 0; ior < nor; ++ior) {
            // cycle through the links for this orphan, pushing everything to nbrs. It is the
            // nbr that does the checking before pushing potential into final...
            while ((it != orphanLinkSet.end())&&(it->first == ior)) {
              const int ior_nbr = it->second; assert((ior_nbr >= 0)&&(ior_nbr < nor_g)); // nbrs could be non-local
              ++it;
              if (ior_nbr < nor) {
                for (set<int>::iterator it2 = cvoor_set[ior].begin(); it2 != cvoor_set[ior].end(); ++it2) {
                  // only add this to nbr's potential if it is not in final or in the set about to be added to final...
                  if ((final_cvoor_set[ior_nbr].find(*it2) == final_cvoor_set[ior_nbr].end())&&
                      (cvoor_set[ior_nbr].find(*it2) == cvoor_set[ior_nbr].end())) {
                    potential_cvoor_set[ior_nbr].insert(*it2);
                  }
                }
              }
              else if (!cvoor_set[ior].empty()) {
                // send over cvoor_set[ior]...
                int rank,bits,index;
                BitUtils::unpackRankBitsIndex(rank,bits,index,rbo_g[ior_nbr-nor]);
                send_count[rank] += 1+cvoor_set[ior].size();
                send_count2[rank] += 3*cvoor_set[ior].size();
              }
            }
          }
          assert(it == orphanLinkSet.end());

          // now the parallel...

          send_disp[0] = 0;
          for (int rank = 1; rank < mpi_size; ++rank)
            send_disp[rank] = send_disp[rank-1] + send_count[rank-1];
          int send_count_sum = send_disp[mpi_size-1] + send_count[mpi_size-1];
          uint8* send_buf_uint8 = new uint8[send_count_sum];

          send_disp2[0] = 0;
          for (int rank = 1; rank < mpi_size; ++rank)
            send_disp2[rank] = send_disp2[rank-1] + send_count2[rank-1];
          int send_count_sum2 = send_disp2[mpi_size-1] + send_count2[mpi_size-1];
          double* send_buf_double = new double[send_count_sum2];

          it = orphanLinkSet.begin();
          for (int ior = 0; ior < nor; ++ior) {
            // cycle through the links for this orphan, pushing everything to nbrs. It is the
            // nbr that does the checking before pushing potential into final...
            while ((it != orphanLinkSet.end())&&(it->first == ior)) {
              const int ior_nbr = it->second; assert((ior_nbr >= 0)&&(ior_nbr < nor_g)); // nbrs could be non-local
              ++it;
              if ((!cvoor_set[ior].empty())&&(ior_nbr >= nor)) {
                // send over cvoor_set[ior]...
                int rank,bits,index;
                BitUtils::unpackRankBitsIndex(rank,bits,index,rbo_g[ior_nbr-nor]);
                const int flip_bits = BitUtils::flipPeriodicBits(bits);
                send_buf_uint8[send_disp[rank]++] = BitUtils::packRankIndex((int)cvoor_set[ior].size(),index); // note put count not rank
                for (set<int>::iterator it2 = cvoor_set[ior].begin(); it2 != cvoor_set[ior].end(); ++it2) {
                  if (*it2 >= 0) {
                    send_buf_uint8[send_disp[rank]++] = BitUtils::packRankBitsIndex(mpi_rank,flip_bits,*it2);
                    FOR_I3 send_buf_double[send_disp2[rank]++] = x_vd[*it2][i];
                    if (flip_bits) PeriodicData::periodicTranslate(send_buf_double+send_disp2[rank]-3,1,flip_bits);
                  }
                  else {
                    if (flip_bits) {
                      int rank0,bits0,index0;
                      BitUtils::unpackRankBitsIndex(rank0,bits0,index0,rbi_or_g[-(*it2)-1]);
                      send_buf_uint8[send_disp[rank]++] = BitUtils::packRankBitsIndex(rank0,BitUtils::addPeriodicBits(bits0,flip_bits),index0);
                    }
                    else {
                      send_buf_uint8[send_disp[rank]++] = rbi_or_g[-(*it2)-1];
                    }
                    FOR_I3 send_buf_double[send_disp2[rank]++] = xcv_or_g[(-(*it2)-1)*3+i];
                    if (flip_bits) PeriodicData::periodicTranslate(send_buf_double+send_disp2[rank]-3,1,flip_bits);
                  }
                }
              }
            }
            // now move cvoor_set into final...
            for (set<int>::iterator it2 = cvoor_set[ior].begin(); it2 != cvoor_set[ior].end(); ++it2)
              final_cvoor_set[ior].insert(*it2);
            cvoor_set[ior].clear();
          }
          assert(it == orphanLinkSet.end());

          // rewind...

          send_disp[0] = 0;
          for (int rank = 1; rank < mpi_size; ++rank)
            send_disp[rank] = send_disp[rank-1] + send_count[rank-1];

          MPI_Alltoall(send_count,1,MPI_INT,recv_count,1,MPI_INT,mpi_comm);

          recv_disp[0] = 0;
          for (int rank = 1; rank < mpi_size; ++rank)
            recv_disp[rank] = recv_count[rank-1] + recv_disp[rank-1];
          int recv_count_sum = recv_disp[mpi_size-1] + recv_count[mpi_size-1];

          uint8 * recv_buf_uint8 = new uint8[recv_count_sum];
          MPI_Alltoallv(send_buf_uint8,send_count,send_disp,MPI_UINT8,
              recv_buf_uint8,recv_count,recv_disp,MPI_UINT8,mpi_comm);
          delete[] send_buf_uint8;

          // rewind...

          send_disp2[0] = 0;
          for (int rank = 1; rank < mpi_size; ++rank)
            send_disp2[rank] = send_disp2[rank-1] + send_count2[rank-1];

          MPI_Alltoall(send_count2,1,MPI_INT,recv_count2,1,MPI_INT,mpi_comm);

          recv_disp2[0] = 0;
          for (int rank = 1; rank < mpi_size; ++rank)
            recv_disp2[rank] = recv_count2[rank-1] + recv_disp2[rank-1];
          int recv_count_sum2 = recv_disp2[mpi_size-1] + recv_count2[mpi_size-1];

          double * recv_buf_double = new double[recv_count_sum2];
          MPI_Alltoallv(send_buf_double,send_count2,send_disp2,MPI_DOUBLE,
              recv_buf_double,recv_count2,recv_disp2,MPI_DOUBLE,mpi_comm);
          delete[] send_buf_double;

          // expand cvoor_set, rbi_or_g, xcv_or_g...
          int irecv = 0;
          int irecv2 = 0;
          while (irecv < recv_count_sum) {
            assert(irecv2 < recv_count_sum2);
            int count,ior_nbr;
            BitUtils::unpackRankIndex(count,ior_nbr,recv_buf_uint8[irecv++]);
            assert((ior_nbr >= 0)&&(ior_nbr < nor));
            for (int jj = 0; jj < count; ++jj) {
              const uint8 rbi = recv_buf_uint8[irecv++];
              int rank,bits,index;
              BitUtils::unpackRankBitsIndex(rank,bits,index,rbi);
              if ((bits != 0)||(rank != mpi_rank)) {
                // add this nbr to ior's cv set w/ -1 idx...
                map<const uint8,int>::iterator it = rbiOrMap.find(rbi);
                if (it == rbiOrMap.end()) {
                  // this is a ghost not currently in the rbi...
                  potential_cvoor_set[ior_nbr].insert(-rbi_or_g.size()-1);
                  rbiOrMap[rbi] = rbi_or_g.size();
                  rbi_or_g.push_back(rbi);
                  //const int icv = cgoor[ior_nbr][0];
                  //cout << COUT_VEC(x_vd[icv]) << " " << COUT_VEC(recv_buf_double+irecv) << endl;
                  FOR_I3 xcv_or_g.push_back(recv_buf_double[irecv2++]);
                }
                else {
                  // only add this to nbr's potential if it is not in final...
                  const int idx = -it->second-1;
                  assert(cvoor_set[ior_nbr].empty());
                  if (final_cvoor_set[ior_nbr].find(idx) == final_cvoor_set[ior_nbr].end())
                    potential_cvoor_set[ior_nbr].insert(idx);
                  irecv2 += 3;
                }
              }
              else {
                assert(index < ncv);
                if (final_cvoor_set[ior_nbr].find(index) == final_cvoor_set[ior_nbr].end())
                  potential_cvoor_set[ior_nbr].insert(index);
                irecv2 += 3;
              }
            }
          }
          assert(irecv == recv_count_sum);
          assert(irecv2 == recv_count_sum2);
          delete[] recv_buf_uint8;
          delete[] recv_buf_double;

          // now loop again moving potential to final when potential potentially cuts the orphan...
          double xmid[3];
          double nmid[3];
          for (int ior = 0; ior < nor; ++ior) {
            // and move potential into cvoor_set in preparation for the next iteration...
            // anyone that (potentially) cuts the orphan should be included...
            // the parent cv of this orphan ior...
            const int icv = cgoor[ior][0];
            const int igr = cgoor[ior][1];
            const int ii = ivrdocv[icv];
            assert(vdReturnDataVec[ii].icv == icv);
            assert((igr >= 1)&&(igr < vdReturnDataVec[ii].ngr));
            // the node range for this group is not explicitly stored, but is contiguous. Compute it
            // once outside the loop...
            int ino_min = vdReturnDataVec[ii].nno;
            int ino_max = 0;
            for (int ied = vdReturnDataVec[ii].edogr_i[igr]; ied < vdReturnDataVec[ii].edogr_i[igr+1]; ++ied) {
              ino_min = min(ino_min,min(vdReturnDataVec[ii].nooed[ied][0],vdReturnDataVec[ii].nooed[ied][1]));
              ino_max = max(ino_max,max(vdReturnDataVec[ii].nooed[ied][0],vdReturnDataVec[ii].nooed[ied][1]));
            }
            for (set<int>::iterator it_pot = potential_cvoor_set[ior].begin(); it_pot != potential_cvoor_set[ior].end(); ++it_pot) {
              const int icv_nbr_pot = *it_pot;
              double * this_x_vd_pot;
              if (icv_nbr_pot >= 0) {
                assert(icv_nbr_pot < ncv);
                this_x_vd_pot = x_vd[icv_nbr_pot];
              }
              else {
                assert((-icv_nbr_pot-1) < int(rbi_or_g.size()));
                this_x_vd_pot = &xcv_or_g[(-icv_nbr_pot-1)*3];
              }
              // we can eliminate a nbr if all nodes of this group are positive wrt the bisector
              // for ANY of the final nbrs...
              bool b_final_closer = false;
              // only do this checking on cvs that are not the orphan's cv. For the case of the orphan cv, if
              // that is part of the final_cvoor_set we need to add it for sure...
              if (icv_nbr_pot != icv) {
                for (set<int>::iterator it_fin = final_cvoor_set[ior].begin(); it_fin != final_cvoor_set[ior].end(); ++it_fin) {
                  const int icv_nbr_fin = *it_fin;
                  assert(icv_nbr_pot != icv_nbr_fin);
                  double * this_x_vd_fin;
                  if (icv_nbr_fin >= 0) {
                    assert(icv_nbr_fin < ncv);
                    this_x_vd_fin = x_vd[icv_nbr_fin];
                  }
                  else {
                    assert((-icv_nbr_fin-1) < int(rbi_or_g.size()));
                    this_x_vd_fin = &xcv_or_g[(-icv_nbr_fin-1)*3];
                  }
                  // the midpoint position and normal in the reference frame of the original orphan is...
                  FOR_I3 xmid[i] = 0.5*(this_x_vd_pot[i]+this_x_vd_fin[i]) - x_vd[icv][i];
                  FOR_I3 nmid[i] = this_x_vd_pot[i] - this_x_vd_fin[i];
                  b_final_closer = true;
                  // using edges to get to nodes, this will be 2x as expensive, but hopefully not the controlling cost...
                  for (int ino = ino_min; ino <= ino_max; ++ino) {
                    if ((vdReturnDataVec[ii].x_no[ino][0]-xmid[0])*nmid[0] +
                        (vdReturnDataVec[ii].x_no[ino][1]-xmid[1])*nmid[1] +
                        (vdReturnDataVec[ii].x_no[ino][2]-xmid[2])*nmid[2] > 0.0) {
                      b_final_closer = false;
                      break;
                    }
                  }
                  // if we find even a single final that is closer to all nodes...
                  if (b_final_closer)
                    break;
                }
              }
              // if there was a case where the for now just add it to cvoor_set for the next iteration...
              if (!b_final_closer) {
                cvoor_set[ior].insert(*it_pot);
                my_done = 0;
              }
            }
            potential_cvoor_set[ior].clear();
          }
          MPI_Allreduce(&my_done,&done,1,MPI_INT,MPI_MIN,mpi_comm);
        }

        // clean up intermediate sets...
        delete[] rbo_g;
        rbiOrMap.clear();
        orphanLinkSet.clear();
        delete[] cvoor_set;
        delete[] potential_cvoor_set;

        if (false) {
          FOR_RANK {
            if (mpi_rank == rank) {
              for (int ior = 0; ior < nor; ++ior) {
                cout << "XXXXX: " << rank << " " << ior << " " << cgoor[ior][0] << " " << COUT_VEC(x_vd[cgoor[ior][0]]) << endl;
                for (set<int>::iterator it2 = final_cvoor_set[ior].begin(); it2 != final_cvoor_set[ior].end(); ++it2) {
                  cout << *it2 << " ";
                  if (*it2 <= -1) {
                    int rank_,bits,index;
                    BitUtils::unpackRankBitsIndex(rank_,bits,index,rbi_or_g[-(*it2)-1]);
                    cout << index << " " << rank_ << " " << COUT_VEC(&xcv_or_g[(-(*it2)-1)*3]) << endl;
                  }
                  else {
                    cout << COUT_VEC(x_vd[*it2]) << endl;
                  }
                }
                cout.flush();
              }
            }
            MPI_Barrier(mpi_comm);
          }
        }

        // check that every orphan has atleast one valid nbr to join. If this is not
        // the case, then it means that a collection of orphans are isolated completely
        // from the rest of the domain, and can probably be eliminated...
        int my_ierr = 0;
        for (int ior = 0; ior < nor; ++ior) {
          if (final_cvoor_set[ior].empty()) {
            my_ierr = -1;
            break;
          }
        }
        int ierr;
        MPI_Reduce(&my_ierr,&ierr,1,MPI_INT,MPI_MIN,0,mpi_comm);
        if ((mpi_rank == 0)&&(ierr != 0))
          cout << " > WARNING: some orphans had no main nbrs to make orphan chunks." << endl;

        delete[] oroii_i; oroii_i = NULL;
        stonboor_i.clear();
        stonboor_v.clear();
        delete[] nboor_i;
        nboor_v.clear();

        int nxvd = rbi_or_g.size();

        CuttableVoronoiData cvd,cvd_copy;

        // update the current centroids using parent vd's in x_cv...

        assert(intVec.empty());
        assert(doubleVec.empty());

        double my_vol_err = 0.0;
        vector<pair<double,int> > nbrVec;
        for (int ior = 0; ior < nor; ++ior) {
          if (final_cvoor_set[ior].empty())
            continue;

          const int icv_or = cgoor[ior][0];
          const int igr_or = cgoor[ior][1];
          const int ii_or = ivrdocv[icv_or];
          assert(igr_or < vdReturnDataVec[ii_or].ngr);
          assert(vdReturnDataVec[ii_or].icv == icv_or);
          vdReturnDataVec[ii_or].extractGroupAsCvd(cvd,igr_or);

          const double vol = cvd.calcVolume();
          //cout << "Extracting ior " << ior << " volume: " << vol << endl;

          // to reduce the set we need to cut against, we compute the
          // approximate location of the center of the orphan and find
          // the closest site...
          double x_no_avg[3] = { 0.0, 0.0, 0.0 };
          for (int ino = 0; ino < cvd.nno; ++ino)
            FOR_I3 x_no_avg[i] += cvd.x_no[ino][i];
          FOR_I3 x_no_avg[i] = x_vd[icv_or][i] + x_no_avg[i]/double(cvd.nno);
          // now build the nbrs to cut against...
          assert(nbrVec.empty());
          int ii_min = -1;
          double d2_min;
          bool b_self_orphan = false;
          for (set<int>::iterator it = final_cvoor_set[ior].begin(); it != final_cvoor_set[ior].end(); ++it) {
            const int ixvd = *it;
            if (ixvd >= 0) {
              // local...
              const double d2 = DIST2(x_no_avg,x_vd[ixvd]);
              if ((ii_min == -1)||(d2 < d2_min)) {
                ii_min = nbrVec.size();
                d2_min = d2;
              }
              nbrVec.push_back(pair<double,int>(d2,ixvd));
              if (ixvd == icv_or) {
                b_self_orphan = true;
                cout << " > local self orphan " << endl;
                break;
              }
            }
            else {
              const double d2 = DIST2(x_no_avg,&xcv_or_g[(-ixvd-1)*3]);
              if ((ii_min == -1)||(d2 < d2_min)) {
                ii_min = nbrVec.size();
                d2_min = d2;
              }
              nbrVec.push_back(pair<double,int>(d2,ixvd));
              // if it is a self orphan it needs to be on your rank. maybe through the
              // orphan linking process it is stored as an rbi_or_g with rank = mpi_rank
              int rank,bits,index;
              BitUtils::unpackRankBitsIndex(rank,bits,index,rbi_or_g[-ixvd-1]);
              //assert(bits == 0);
              if ((index == icv_or)&&(rank == mpi_rank)) {
                b_self_orphan = true;
                cout << " > non-local self orphan " << endl;
                break;
              }
            }
          }
          final_cvoor_set[ior].clear();

          // if this is a "self-orphan", then no cutting is necessary, and we should join
          // the main cv: icv_or. We could do this by refactoring the main voronoi diagram, but
          // recall that we don't want to change these fundamentally until the user has
          // requested writing of the mles (i.e. the fuseOrphanChunks routine), so we treat this
          // as a chunk...

          if (b_self_orphan) {

            double vol_chunk,x_chunk[3];
            cvd.calcVolumeAndCentroid(vol_chunk,x_chunk);

            // this is a local nbr, so update the volume and centroid with our
            // contribution directly...
            assert(cv_check[icv_or] == 1);
            vol_cv[icv_or] += vol_chunk;
            FOR_I3 x_cv[icv_or][i] += vol_chunk*x_chunk[i]; // recall x_cv is currently vol*dx_cv relative to x_vd

            // we will fuse this at the node level later once all orphans have been collected
            // at a particular cv. For now, just add it as an orphan chunk...
            if (b_final) {

              assert(vdReturnDataVec[ii_or].icv == icv_or);

              const int iocd = ocdVec.size();
              ocdVec.push_back(OrphanChunkData());
              OrphanChunkData * ocd = &ocdVec.back();

              // insert this into vdReturnDataVec[ii_or]...
              ocd->next = vdReturnDataVec[ii_or].orphan_chunk_data; // might be -1
              vdReturnDataVec[ii_or].orphan_chunk_data = iocd;

              // we also need the originating rbi for the orphan chunk for fuse
              ocd->rbi = BitUtils::packRankBitsIndex(mpi_rank,0,icv_or);

              // and copy the cvd info into ocd...
              ocd->nno = cvd.nno;
              ocd->x_no = new double[cvd.nno][3];
              for (int ino = 0; ino < cvd.nno; ++ino)
                FOR_I3 ocd->x_no[ino][i] = cvd.x_no[ino][i];
              ocd->ned = cvd.ned;
              ocd->nooed = new int[cvd.ned][2];
              for (int ied = 0; ied < cvd.ned; ++ied)
                FOR_I2 ocd->nooed[ied][i] = cvd.nooed[ied][i];
              ocd->faoed = new int8[ocd->ned][2];
              for (int ied = 0; ied < cvd.ned; ++ied) {
                FOR_I2 {
                  const int ifa = cvd.faoed[ied][i];
                  if (ifa >= 0) {
                    const int ist = vdReturnDataVec[ii_or].spbobf[ifa][0];
                    const int ipart_bits = vdReturnDataVec[ii_or].spbobf[ifa][1];
                    // TODO: use a pack routine here?...
                    ocd->faoed[ied][i] = ((int8(ipart_bits)<<32)|(int8(ist)));
                    assert(ocd->faoed[ied][i] >= 0);
                  }
                  else {
                    assert(ifa > -ORPHAN_FACE_OFFSET);
                    // this ifa is -1 index'd into cvofa...
                    const int icv_nbr = vdReturnDataVec[ii_or].cvofa[-ifa-1];
                    //assert((icv_nbr >= 0)&&(icv_nbr < ncv));
                    if (icv_nbr >= 0) {
                      if (icv_nbr < ncv) {
                        ocd->faoed[ied][i] = -int8(BitUtils::packRankBitsIndex(mpi_rank,0,icv_nbr))-1;
                      }
                      else {
                        assert((icv_nbr-ncv) < rbi_g.size());
                        ocd->faoed[ied][i] = -int8(rbi_g[icv_nbr-ncv])-1;
                      }
                    }
                    else {
                      assert((-icv_nbr-1) < nxvd);
                      ocd->faoed[ied][i] = -int8(rbi_or_g[-icv_nbr-1])-1;
                    }
                    assert(ocd->faoed[ied][i] < 0);
                  }
                }
              }

            }

          }
          else { // !b_self_orphan

            // !b_self_orphan, so cut into chunks...

            assert(!nbrVec.empty());
            assert(ii_min != -1);
            double * x_closest;
            if (nbrVec[ii_min].second >= 0) {
              x_closest = x_vd[nbrVec[ii_min].second];
            }
            else {
              x_closest = &xcv_or_g[(-nbrVec[ii_min].second-1)*3];
            }

            // now formally compute the max distance to the furthest orphan node
            // from this closest site. Because of the Euclidean definition of
            // Voronoi diagrams, it will not be possible for anyone outside
            // twice this distance to cut this orphan, no matter which reference
            // frame it is moved into...
            double d2_max = 0.0;
            for (int ino = 0; ino < cvd.nno; ++ino)  {
              double d2 = 0;
              FOR_I3 {
                const double dx = cvd.x_no[ino][i] + (x_vd[icv_or][i]-x_closest[i]);
                d2 += dx*dx;
              }
              d2_max = max(d2_max,d2);
            }
            for (int ii = 0; ii < nbrVec.size(); ++ii) {
              // now change the first distance = the distance to the closest guy...
              if (nbrVec[ii].second >= 0) {
                nbrVec[ii].first = DIST2(x_closest,x_vd[nbrVec[ii].second]);
              }
              else {
                nbrVec[ii].first = DIST2(x_closest,&xcv_or_g[(-nbrVec[ii].second-1)*3]);
              }
            }

            // and sort. We will retain this order throughout the cutting. It is actually
            // worse to re-sort based on distance from the current zero.
            sort(nbrVec.begin(),nbrVec.end());
            int nii;
            for (nii = 0; nii < nbrVec.size(); ++nii) {
              if (nbrVec[nii].first > 4.0*d2_max)
                break;
            }

            double vol_sum = 0.0;

            for (int ii = 0; ii < nii; ++ii) {
              //cout << "working on ii: " << ii << endl;
              double * x_vd_ii;
              if (nbrVec[ii].second >= 0) {
                x_vd_ii = x_vd[nbrVec[ii].second];
              }
              else {
                x_vd_ii = &xcv_or_g[(-nbrVec[ii].second-1)*3];
              }

              // take a copy
              cvd_copy.init(cvd);
              for (int ino = 0; ino < cvd_copy.nno; ++ino) {
                FOR_I3 cvd_copy.x_no[ino][i] += (x_vd[icv_or][i]-x_vd_ii[i]);
              }

              // hack: change -ve ifa's back to -8 indexing - this prevents cutting
              // problems with the low indexing like -1, etc. TODO: should get rid
              // of this behavior in the cvd.cut routine and make it faster and general.
              for (int ied = 0; ied < cvd_copy.ned; ++ied) {
                FOR_I2 {
                  const int ifa = cvd_copy.faoed[ied][i];
                  if (ifa < 0) {
                    cvd_copy.faoed[ied][i] = ifa-7;
                  }
                }
              }

              cvd_copy.setD2Max();
              for (int jj = 0; jj < nii; ++jj) {
                if (jj != ii) {

                  double * x_vd_jj;
                  if (nbrVec[jj].second >= 0) {
                    x_vd_jj = x_vd[nbrVec[jj].second];
                  }
                  else {
                    x_vd_jj = &xcv_or_g[(-nbrVec[jj].second-1)*3];
                  }

                  double dn[3]; FOR_I3 dn[i] = 0.5*(x_vd_jj[i] - x_vd_ii[i]);
                  // only cut if we are inside the 6 paraboloids...
                  if ((2.0*dn[0] > cvd_copy.Lmax[0])||
                      (2.0*dn[1] > cvd_copy.Lmax[1])||
                      (2.0*dn[2] > cvd_copy.Lmax[2])||
                      (2.0*dn[0] < cvd_copy.Lmin[0])||
                      (2.0*dn[1] < cvd_copy.Lmin[1])||
                      (2.0*dn[2] < cvd_copy.Lmin[2]) ) {
                    continue;
                  }
                  // if we got here, we are cutting...
                  cvd_copy.cut(dn,-jj-ORPHAN_FACE_OFFSET); // cut orphans so they don't interfere with existing cuts...
                  if (cvd_copy.empty())
                    break;
                }
              }

              if (!cvd_copy.empty()) {

                // this guy's volume and centroid has to be sent to "ii":

                double vol_chunk,x_chunk[3];
                cvd_copy.calcVolumeAndCentroid(vol_chunk,x_chunk);

                if (nbrVec[ii].second >= 0) {

                  // this is a local nbr, so update the volume and centroid with our
                  // contribution directly...
                  const int icv_fp = nbrVec[ii].second; // fp = foster-parent: the new parent that gets the cut chunk
                  assert(icv_fp != icv_or);
                  assert((icv_fp >= 0)&&(icv_fp < ncv));
                  assert(cv_check[icv_fp] == 1);
                  vol_cv[icv_fp] += vol_chunk;
                  FOR_I3 x_cv[icv_fp][i] += vol_chunk*x_chunk[i]; // recall x_cv is currently vol*dx_cv relative to x_vd
                  // we will fuse this at the node level later once all orphans have been collected
                  // at a particular cv. For now, just add it as an orphan chunk...
                  if (b_final) {

                    const int ii_fp = ivrdocv[icv_fp];
                    assert(vdReturnDataVec[ii_fp].icv == icv_fp);

                    const int iocd = ocdVec.size();
                    ocdVec.push_back(OrphanChunkData());
                    OrphanChunkData * ocd = &ocdVec.back();

                    // insert this into vdReturnDataVec[ii_fp]...
                    ocd->next = vdReturnDataVec[ii_fp].orphan_chunk_data; // might be -1
                    vdReturnDataVec[ii_fp].orphan_chunk_data = iocd;

                    // we also need the originating rbi for the orphan chunk for fuse
                    ocd->rbi = BitUtils::packRankBitsIndex(mpi_rank,0,icv_or);

                    // and copy the cvd_copy info into ocd...
                    ocd->nno = cvd_copy.nno;
                    ocd->x_no = new double[cvd_copy.nno][3];
                    for (int ino = 0; ino < cvd_copy.nno; ++ino)
                      FOR_I3 ocd->x_no[ino][i] = cvd_copy.x_no[ino][i];
                    ocd->ned = cvd_copy.ned;
                    ocd->nooed = new int[cvd_copy.ned][2];
                    for (int ied = 0; ied < cvd_copy.ned; ++ied)
                      FOR_I2 ocd->nooed[ied][i] = cvd_copy.nooed[ied][i];
                    ocd->faoed = new int8[ocd->ned][2];
                    for (int ied = 0; ied < cvd_copy.ned; ++ied) {
                      FOR_I2 {
                        const int ifa = cvd_copy.faoed[ied][i];
                        if (ifa >= 0) {
                          const int ist = vdReturnDataVec[ii_or].spbobf[ifa][0];
                          const int ipart_bits = vdReturnDataVec[ii_or].spbobf[ifa][1];
                          // TODO: use a pack routine here?...
                          ocd->faoed[ied][i] = ((int8(ipart_bits)<<32)|(int8(ist)));
                          assert(ocd->faoed[ied][i] >= 0);
                        }
                        else if (ifa <= -ORPHAN_FACE_OFFSET) {
                          // this face was cut during the current orphan-chunk cutting. That means its icv associated with
                          // cutting is actually the foster parent...
                          const int icv_nbr = nbrVec[-ifa-ORPHAN_FACE_OFFSET].second;
                          if (icv_nbr >= 0) {
                            assert(icv_nbr < ncv);
                            ocd->faoed[ied][i] = -int8(BitUtils::packRankBitsIndex(mpi_rank+mpi_size,0,icv_nbr))-1;
                          }
                          else {
                            assert((-icv_nbr-1) < nxvd);
                            int rank,bits,index;
                            BitUtils::unpackRankBitsIndex(rank,bits,index,rbi_or_g[-icv_nbr-1]);
                            //assert(bits == 0);
                            //if (bits) {
                            //  cout << " rbi_or_g " << endl;
                            //  double fwd[3],bwd[3]; FOR_I3 fwd[i] = bwd[i] = xcv_or_g[3*(-icv_nbr-1)+i];
                            //  cout << COUT_VEC(fwd) << "->";
                            //  PeriodicData::periodicTranslate(fwd,1,bits);
                            //  cout << COUT_VEC(fwd) << endl;
                            //  cout << COUT_VEC(bwd) << "<-";
                            //  PeriodicData::periodicTranslate(bwd,1,BitUtils::flipPeriodicBits(bits));
                            //  cout << COUT_VEC(bwd) << endl;
                            //}
                            ocd->faoed[ied][i] = -int8(BitUtils::packRankBitsIndex(rank+mpi_size,bits,index))-1;
                          }
                          assert(ocd->faoed[ied][i] < 0);
                        }
                        else {
                          // this ifa is -8 index'd into cvofa: recall it was adjusted above to
                          // allow recutting with the cut routine...
                          const int icv_nbr = vdReturnDataVec[ii_or].cvofa[-ifa-8];
                          assert(icv_nbr >= 0); // original ghosts (not from orphan linking)
                          if (icv_nbr < ncv) {
                            ocd->faoed[ied][i] = -int8(BitUtils::packRankBitsIndex(mpi_rank,0,icv_nbr))-1;
                          }
                          else {
                            assert(icv_nbr-ncv < rbi_g.size());
                            ocd->faoed[ied][i] = -int8(rbi_g[icv_nbr-ncv])-1;
                          }
                          assert(ocd->faoed[ied][i] < 0);
                        }
                      }
                    }

                  }

                }
                else {

                  // parallel...

                  const int ixvd = -nbrVec[ii].second-1;
                  assert((ixvd >= 0)&&(ixvd < nxvd));
                  int rank,bits,index;
                  BitUtils::unpackRankBitsIndex(rank,bits,index,rbi_or_g[ixvd]);
                  const int flip_bits = BitUtils::flipPeriodicBits(bits);
                  if (flip_bits) PeriodicData::periodicRotate(x_chunk,1,flip_bits);

                  intVec.push_back(ixvd);
                  doubleVec.push_back(vol_chunk);
                  FOR_I3 doubleVec.push_back(vol_chunk*x_chunk[i]);

                  if (b_final) {


                    intVec.push_back(BitUtils::packRankBits(mpi_rank,flip_bits));
                    intVec.push_back(icv_or);

                    intVec.push_back(cvd_copy.nno);
                    if (flip_bits) PeriodicData::periodicRotate(cvd_copy.x_no,cvd_copy.nno,flip_bits);
                    for (int ino = 0; ino < cvd_copy.nno; ++ino)
                      FOR_I3 doubleVec.push_back(cvd_copy.x_no[ino][i]);
                    intVec.push_back(cvd_copy.ned);
                    for (int ied = 0; ied < cvd_copy.ned; ++ied)
                      FOR_I2 intVec.push_back(cvd_copy.nooed[ied][i]);
                    for (int ied = 0; ied < cvd_copy.ned; ++ied) {
                      FOR_I2 {
                        const int ifa = cvd_copy.faoed[ied][i];
                        if (ifa >= 0) {
                          intVec.push_back(vdReturnDataVec[ii_or].spbobf[ifa][0]); // ist
                          const int part_st = (vdReturnDataVec[ii_or].spbobf[ifa][1]>>6);
                          const int bits_st = (vdReturnDataVec[ii_or].spbobf[ifa][1]&MASK_6BITS);
                          const int bits_tot = BitUtils::addPeriodicBits(flip_bits,bits_st);
                          assert((bits_tot>>6) == 0); // should be <= 1 translation
                          intVec.push_back((part_st<<6)|bits_tot);
                        }
                        else if (ifa <= -ORPHAN_FACE_OFFSET) {
                          // this face was cut during the current orphan-chunk cutting. That means its icv associated with
                          // cutting is actually the foster parent...
                          const int icv_nbr = nbrVec[-ifa-ORPHAN_FACE_OFFSET].second;
                          if (icv_nbr >= 0) {
                            assert(icv_nbr < ncv);
                            intVec.push_back(-BitUtils::packRankBits(mpi_rank+mpi_size,flip_bits)-1); // -1 index for ifa
                            intVec.push_back(icv_nbr);
                          }
                          else {
                            assert((-icv_nbr-1) < nxvd);
                            int rank0,bits0,index0;
                            BitUtils::unpackRankBitsIndex(rank0,bits0,index0,rbi_or_g[-icv_nbr-1]);
                            const int bits_tot = BitUtils::addPeriodicBits(flip_bits,bits0);
                            intVec.push_back(-BitUtils::packRankBits(rank0+mpi_size,bits_tot)-1); // -1 index for ifa
                            intVec.push_back(index0);
                          }
                        }
                        else {
                          // this ifa is -8 index'd into cvofa: recall it was adjusted above to
                          // allow recutting with the cut routine...
                          const int icv_nbr = vdReturnDataVec[ii_or].cvofa[-ifa-8];
                          assert(icv_nbr >= 0); // original ghost (not from orphan linking)
                          if (icv_nbr < ncv) {
                            intVec.push_back(-BitUtils::packRankBits(mpi_rank,flip_bits)-1); // -1 index for ifa
                            intVec.push_back(icv_nbr);
                          }
                          else {
                            assert((icv_nbr-ncv) < rbi_g.size());
                            int rank0,bits0,index0;
                            BitUtils::unpackRankBitsIndex(rank0,bits0,index0,rbi_g[icv_nbr-ncv]);
                            const int bits_tot = BitUtils::addPeriodicBits(flip_bits,bits0);
                            intVec.push_back(-BitUtils::packRankBits(rank0,bits_tot)-1); // -1 index for ifa
                            intVec.push_back(index0);
                          }
                        }
                      }
                    }

                  }

                }

                // this is ii's part...
                vol_sum += vol_chunk;

                //cout << " > cvd_copy volume: " << ii << " " << vol_chunk << endl;

                try {
                  cvd_copy.check();
                }
                catch(...) {
                  cout << "cvd_copy failed check!" << endl;
                  assert(0);
                  /*
                     char filename[128];
                     sprintf(filename,"cvd_debug.%06d.bin",mpi_rank);
                     cvd_copy.writeBinary(filename);
                     double zero[3] = { 0.0, 0.0, 0.0 };
                     cvd_copy.writeTecplot(mpi_rank,zero);
                     throw(-1);
                     */
                }

              }

            }

            //cout << "vol diff: " << vol-vol_sum << endl;
            my_vol_err += vol-vol_sum;

          }

          nbrVec.clear();
          cvd.clear();

        } // ior loop
        delete[] final_cvoor_set;
        delete[] cgoor; cgoor = NULL;

        double vol_err;
        MPI_Reduce(&my_vol_err,&vol_err,1,MPI_DOUBLE,MPI_SUM,0,mpi_comm);
        if (mpi_rank == 0) cout << " > vol_err = " << vol_err << endl;

        // now send out the centroid corrections...

        FOR_RANK {
          send_count[rank] = 0;
          send_count2[rank] = 0;
        }

        {
          int ii = 0;
          while (ii < intVec.size()) {
            const int ixvd = intVec[ii++];
            assert((ixvd >= 0)&&(ixvd < nxvd));
            int rank,bits,index;
            BitUtils::unpackRankBitsIndex(rank,bits,index,rbi_or_g[ixvd]);
            send_count[rank]  += 1; // index
            send_count2[rank] += 4; // vol,vol*dx_cv
            if (b_final) {
              ii += 2; // rank_bits,icv_or
              send_count[rank] += 2; // rank_bits,icv_or
              const int nno = intVec[ii++];
              const int ned = intVec[ii++];
              ii += 6*ned;
              send_count[rank] += 2+6*ned; // nno,ned,nooed's(x2),faoed's(x4)
              send_count2[rank] += 3*nno; // x_no's
            }
          }
        }

        send_disp[0] = 0;
        for (int rank = 1; rank < mpi_size; ++rank)
          send_disp[rank] = send_disp[rank-1] + send_count[rank-1];
        send_count_sum = send_disp[mpi_size-1] + send_count[mpi_size-1];
        assert(send_count_sum == intVec.size());

        send_disp2[0] = 0;
        for (int rank = 1; rank < mpi_size; ++rank)
          send_disp2[rank] = send_disp2[rank-1] + send_count2[rank-1];
        const int send_count_sum2 = send_disp2[mpi_size-1] + send_count2[mpi_size-1];
        assert(send_count_sum2 == doubleVec.size());

        assert(send_buf_int == NULL); send_buf_int = new int[send_count_sum];
        assert(send_buf_double == NULL); send_buf_double = new double[send_count_sum2];

        {
          int ii = 0;
          int ii2 = 0;
          while (ii < intVec.size()) {
            assert(ii2 < doubleVec.size());
            const int ixvd = intVec[ii++];
            assert((ixvd >= 0)&&(ixvd < nxvd));
            int rank,bits,index;
            BitUtils::unpackRankBitsIndex(rank,bits,index,rbi_or_g[ixvd]);
            //assert(bits == 0);
            send_buf_int[send_disp[rank]++] = index;
            send_buf_double[send_disp2[rank]++] = doubleVec[ii2++];
            send_buf_double[send_disp2[rank]++] = doubleVec[ii2++];
            send_buf_double[send_disp2[rank]++] = doubleVec[ii2++];
            send_buf_double[send_disp2[rank]++] = doubleVec[ii2++];
            if (b_final) {
              send_buf_int[send_disp[rank]++] = intVec[ii++]; // rank_bits
              send_buf_int[send_disp[rank]++] = intVec[ii++]; // icv_or
              const int nno = intVec[ii++];
              const int ned = intVec[ii++];
              send_buf_int[send_disp[rank]++] = nno;
              for (int ino = 0; ino < nno; ++ino) {
                FOR_I3 send_buf_double[send_disp2[rank]++] = doubleVec[ii2++]; // x_no
              }
              send_buf_int[send_disp[rank]++] = ned;
              for (int ied = 0; ied < ned; ++ied)
                FOR_I2 send_buf_int[send_disp[rank]++] = intVec[ii++]; // nooed
              for (int ied = 0; ied < ned; ++ied)
                FOR_I4 send_buf_int[send_disp[rank]++] = intVec[ii++]; // faoed
            }
          }
          assert(ii == intVec.size());
          assert(ii2 == doubleVec.size());
        }

        // rewind...

        send_disp[0] = 0;
        for (int rank = 1; rank < mpi_size; ++rank)
          send_disp[rank] = send_disp[rank-1] + send_count[rank-1];

        send_disp2[0] = 0;
        for (int rank = 1; rank < mpi_size; ++rank)
          send_disp2[rank] = send_disp2[rank-1] + send_count2[rank-1];

        MPI_Alltoall(send_count,1,MPI_INT,recv_count,1,MPI_INT,mpi_comm);
        MPI_Alltoall(send_count2,1,MPI_INT,recv_count2,1,MPI_INT,mpi_comm);

        recv_disp[0] = 0;
        for (int rank = 1; rank < mpi_size; ++rank)
          recv_disp[rank] = recv_count[rank-1] + recv_disp[rank-1];
        recv_count_sum = recv_disp[mpi_size-1] + recv_count[mpi_size-1];

        recv_disp2[0] = 0;
        for (int rank = 1; rank < mpi_size; ++rank)
          recv_disp2[rank] = recv_count2[rank-1] + recv_disp2[rank-1];
        const int recv_count_sum2 = recv_disp2[mpi_size-1] + recv_count2[mpi_size-1];

        assert(recv_buf_int == NULL); recv_buf_int = new int[recv_count_sum];
        assert(recv_buf_double == NULL); recv_buf_double = new double[recv_count_sum2];

        MPI_Alltoallv(send_buf_int,send_count,send_disp,MPI_INT,
            recv_buf_int,recv_count,recv_disp,MPI_INT,mpi_comm);
        delete[] send_buf_int; send_buf_int = NULL;

        MPI_Alltoallv(send_buf_double,send_count2,send_disp2,MPI_DOUBLE,
            recv_buf_double,recv_count2,recv_disp2,MPI_DOUBLE,mpi_comm);
        delete[] send_buf_double; send_buf_double = NULL;
        delete[] send_count2;
        delete[] send_disp2;
        delete[] recv_count2;
        delete[] recv_disp2;

        {
          int irecv = 0;
          int irecv2 = 0;
          while (irecv < recv_count_sum) {
            assert(irecv2 < recv_count_sum2);
            const int icv_fp = recv_buf_int[irecv++];
            assert((icv_fp >= 0)&&(icv_fp < ncv));
            assert(cv_check[icv_fp] == 1);
            // add the orphan contribution to the volume and centroid correction...
            vol_cv[icv_fp] += recv_buf_double[irecv2++];
            FOR_I3 x_cv[icv_fp][i] += recv_buf_double[irecv2++];
            FOR_I3 assert(x_cv[icv_fp][i] == x_cv[icv_fp][i]);
            if (b_final) {
              const int ii_fp = ivrdocv[icv_fp];
              assert(vdReturnDataVec[ii_fp].icv == icv_fp);

              const int iocd = ocdVec.size();
              ocdVec.push_back(OrphanChunkData());
              OrphanChunkData * ocd = &ocdVec.back();

              // insert this into vdReturnDataVec[ii_fp]...
              ocd->next = vdReturnDataVec[ii_fp].orphan_chunk_data; // might be -1
              vdReturnDataVec[ii_fp].orphan_chunk_data = iocd;

              // need rbi for fuse...
              const int rank_bits = recv_buf_int[irecv++];
              const int icv_or = recv_buf_int[irecv++];
              ocd->rbi = BitUtils::packRankBitsIndex(rank_bits,icv_or);

              // and copy the cvd_copy info into ocd...
              ocd->nno = recv_buf_int[irecv++];
              ocd->x_no = new double[ocd->nno][3];
              for (int ino = 0; ino < ocd->nno; ++ino)
                FOR_I3 ocd->x_no[ino][i] = recv_buf_double[irecv2++];
              ocd->ned = recv_buf_int[irecv++];
              ocd->nooed = new int[ocd->ned][2];
              for (int ied = 0; ied < ocd->ned; ++ied)
                FOR_I2 ocd->nooed[ied][i] = recv_buf_int[irecv++];
              ocd->faoed = new int8[ocd->ned][2];
              for (int ied = 0; ied < ocd->ned; ++ied) {
                FOR_I2 {
                  const int first  = recv_buf_int[irecv++];
                  const int second = recv_buf_int[irecv++];
                  if (first >= 0) {
                    const int ist = first;
                    const int ipart_bits = second;
                    ocd->faoed[ied][i] = ((int8(ipart_bits)<<32)|(int8(ist)));
                    assert(ocd->faoed[ied][i] >= 0);
                  }
                  else {
                    int rank_bits = -first-1; // note rank = rank+mpi_size indicates cut face (in fuse below)
                    const int icv_nbr = second;
                    ocd->faoed[ied][i] = -int8(BitUtils::packRankBitsIndex(rank_bits,icv_nbr))-1;
                    assert(ocd->faoed[ied][i] < 0);
                  }
                }
              }
            }

          }
          assert(irecv == recv_count_sum);
          assert(irecv2 == recv_count_sum2);
          delete[] recv_buf_double;
          delete[] recv_buf_int;
        }

      }

      // and return x_cv back to the true centroid...
      for (int ii = 0; ii < vdReturnDataVec.size(); ++ii) {
        const int icv = vdReturnDataVec[ii].icv;
        assert(cv_check[icv] == 1);
        // put x_cv back to a true centroid in physical coords...
        assert(vol_cv[icv] > 0.0);
        FOR_I3 x_cv[icv][i] = x_cv[icv][i]/vol_cv[icv] + x_vd[icv][i];
        FOR_I3 assert(x_cv[icv][i] == x_cv[icv][i]);
      }

      delete[] cv_check; cv_check = NULL;

      {
        // normalized MAG(dx_cv)/delta_vd...
        double my_sum[3] = { 0.0, 0.0, 0.0 };
        double my_minmax[2] = { HUGE_VAL, HUGE_VAL };
        FOR_ICV {
          assert(delta_vd[icv] > 0.0);
          // normalize by 1/2*delta_vd...
          const double val = 2.0*DIST(x_cv[icv],x_vd[icv])/delta_vd[icv];
          my_sum[0] += 1.0;
          my_sum[1] += val;
          if (val < 1.0E-12) my_sum[2] += 1.0;
          my_minmax[0] = min(my_minmax[0],val);
          my_minmax[1] = min(my_minmax[1],-val);

        }
        double sum[3];
        MPI_Reduce(my_sum,sum,3,MPI_DOUBLE,MPI_SUM,0,mpi_comm);
        double minmax[2];
        MPI_Reduce(my_minmax,minmax,2,MPI_DOUBLE,MPI_MIN,0,mpi_comm);
        if (mpi_rank == 0) {
          cout << " > normalized vd centroid displacement [min,avg,max]: " << minmax[0] <<
            " " << sum[1]/sum[0] << " " << -minmax[1] << endl;
          cout << " > fraction of vds that are centroidal: " << sum[2]/sum[0] << endl;
        }

      }

    }

    void flagSmoothedVds(uint *smooth_flag,int * send_count,int * send_disp,int * recv_count,int * recv_disp) {
      assert(smooth_flag);
      if (mpi_rank == 0)
        cout << "flagSmoothedVds()" << endl;

      if (smooth_mode_idata[0]+smooth_mode_idata[1]+smooth_mode_idata[2]+smooth_mode_idata[3] == 0) {
        // no layers have been set for any subgroups, so this means all cells should be iterated
        // smooth_mode determines whether FORCE inside structured parts or not
        // flag it all
        FOR_ICV smooth_flag[icv] = 1;
      }
      else {
        // TODO: this is perhaps made a little more complex than necessary, because
        // now lloyd iteration can be called in any part points and these are automatically
        // constrained by zero, unless another part-based constraint has been implemented
        // in PartData...

        assert(send_count);
        assert(send_disp);
        assert(recv_count);
        assert(recv_disp);

        FOR_ICV smooth_flag[icv] = 0; // originally do not smooth

        const uint n_blayers = min((uint)256,smooth_mode_idata[0]);
        const uint n_tlayers = min((uint)256,smooth_mode_idata[1]);
        const uint n_players = min((uint)256,smooth_mode_idata[2]);
        const uint n_pblayers = min((uint)256,smooth_mode_idata[3]);

        if (mpi_rank == 0) {
          cout << " > smoothing layer counts" << endl;
          cout << "    > surface layers: " << n_blayers << ", resolution transition layers: " << n_tlayers << ", periodicity-adjacent layers: " << n_players << ", part-boundary layers: " << n_pblayers << endl;
        }

        // flag different types of cells appropriately
        for (int ii = 0; ii < vdReturnDataVec.size(); ++ii) {
          const int icv = vdReturnDataVec[ii].icv;
          assert((icv >= 0)&&(icv < ncv));
          if (n_blayers && (vdReturnDataVec[ii].nbf > 0)) {
            if ((flag_cv[icv]&PART_MASK_BITS) == partVec.size()) {
              smooth_flag[icv] |= (1 << SMOOTH_FLAG_S_OFFSET); // this is a boundary hcp cv
            }
            else {
              // also add any points that are in parts that touch boundaries that they are not stranded from...
              bool smooth_cv = true;
              for (int ibf = 0; ibf < vdReturnDataVec[ii].nbf; ++ibf) {
                const int ist = vdReturnDataVec[ii].spbobf[ibf][0];
                const int ipart = (vdReturnDataVec[ii].spbobf[ibf][1] >> 6);
                //const int bits = (vdReturnDataVec[ii].spbobf[ibf][1]&MASK_6BITS);
                assert(partVec[ipart]);
                int ipart_of_ist;
                if (partVec[ipart]->getPartForSt(ipart_of_ist,ist)) {
                  if ((flag_cv[icv]&PART_MASK_BITS) == ipart_of_ist) {
                    // this part is touching at least one face from which it is extruded, so it should NOT be smoothed...
                    smooth_cv = false;
                    break;
                  }
                }
              }
              // if smooth_icv is still true, then we did not hit any boundary flagged with our part index...
              if (smooth_cv) smooth_flag[icv] |= (1 << SMOOTH_FLAG_S_OFFSET);
            }
          }

          if (n_players) {
            // only worry about HCP periodic points for now; assume structured ones are perfect...
            if ((flag_cv[icv]&PART_MASK_BITS) == partVec.size()) {
              if (vdReturnDataVec[ii].hasPeriodicNbrs())
                smooth_flag[icv] |= (1 << SMOOTH_FLAG_P_OFFSET); // this is a periodic-adjacent hcp cv
            }
          }
        }

        // also, any hcp cvs touching nbrs in parts get flagged. Here we do this in reverse so we
        // only need one exchange. We loop on the cvs from parts and flag any of our hcp nbrs...

        set<uint8> rbiSet;
        set<uint8> rbiNbrSet;
        if (n_pblayers) {
          for (int ii = 0; ii < vdReturnDataVec.size(); ++ii) {
            const int icv = vdReturnDataVec[ii].icv;
            assert((icv >= 0)&&(icv < ncv));
            if ((flag_cv[icv]&PART_MASK_BITS) < partVec.size()) {
              for (int ifa = 0; ifa < vdReturnDataVec[ii].nfa; ++ifa) {
                const int icv_nbr = vdReturnDataVec[ii].cvofa[ifa];
                if (icv_nbr < ncv) {
                  if ((flag_cv[icv_nbr]&PART_MASK_BITS) == partVec.size())
                    smooth_flag[icv_nbr] |= (1 << SMOOTH_FLAG_PB_OFFSET);
                }
                else {
                  int rank,bits,index;
                  BitUtils::unpackRankBitsIndex(rank,bits,index,rbi_g[icv_nbr-ncv]);
                  if (rank == mpi_rank) {
                    assert((index >= 0)&&(index < ncv));
                    if ((flag_cv[index]&PART_MASK_BITS) == partVec.size())
                      smooth_flag[index] |= (1 << SMOOTH_FLAG_PB_OFFSET);
                  }
                  else {
                    rbiSet.insert(rbi_g[icv_nbr-ncv]);
                  }
                }
              }
            }
          }

          FOR_RANK send_count[rank] = 0;

          int * send_buf_int = new int[rbiSet.size()];
          int isend = 0;
          for (set<uint8>::const_iterator iter = rbiSet.begin(); iter != rbiSet.end(); ++iter) {
            int rank,bits,index;
            BitUtils::unpackRankBitsIndex(rank,bits,index,*iter);
            assert(rank != mpi_rank);
            send_buf_int[isend++] = index;
            ++send_count[rank];
          }

          send_disp[0] = 0;
          for (int rank = 1; rank < mpi_size; ++rank)
            send_disp[rank] = send_count[rank-1] + send_disp[rank-1];
          assert(rbiSet.size() == send_disp[mpi_size-1] + send_count[mpi_size-1]);
          rbiSet.clear();

          // set up recv-side stuff...

          MPI_Alltoall(send_count,1,MPI_INT,recv_count,1,MPI_INT,mpi_comm);

          recv_disp[0] = 0;
          for (int rank = 1; rank < mpi_size; ++rank)
            recv_disp[rank] = recv_count[rank-1] + recv_disp[rank-1];
          const int recv_count_sum = recv_disp[mpi_size-1] + recv_count[mpi_size-1];

          int * recv_buf_int = new int[recv_count_sum];
          MPI_Alltoallv(send_buf_int,send_count,send_disp,MPI_INT,
              recv_buf_int,recv_count,recv_disp,MPI_INT,mpi_comm);
          delete[] send_buf_int; send_buf_int = NULL;

          for (int irecv = 0; irecv < recv_count_sum; ++irecv) {
            const int icv = recv_buf_int[irecv];
            assert((icv >= 0)&&(icv < ncv));
            if ((flag_cv[icv]&PART_MASK_BITS) == partVec.size())
              smooth_flag[icv] |= (1 << SMOOTH_FLAG_PB_OFFSET);
          }
          delete[] recv_buf_int; recv_buf_int = NULL;
        }

        if (n_tlayers) {
          double * send_buf_double = NULL;
          int * send_buf_int = NULL;
          FOR_RANK send_count[rank] = 0;
          for (int iter = 0; iter < 2; ++iter) {
            for (int ii = 0; ii < vdReturnDataVec.size(); ++ii) {
              const int icv = vdReturnDataVec[ii].icv;
              assert((icv >= 0)&&(icv < ncv));
              // only on hcp transitions...
              if ((flag_cv[icv]&PART_MASK_BITS) == partVec.size()) {
                for (int ifa = 0; ifa < vdReturnDataVec[ii].nfa; ++ifa) {
                  const int icv_nbr = vdReturnDataVec[ii].cvofa[ifa];
                  if (iter == 0) {
                    if (icv_nbr < ncv) {
                      if ((flag_cv[icv_nbr]&PART_MASK_BITS) == partVec.size()) {
                        const double factor = delta_vd[icv_nbr]/delta_vd[icv];
                        if ((factor > 1.5)||(factor < 0.75))
                          smooth_flag[icv_nbr] |= (1 << SMOOTH_FLAG_T_OFFSET);
                      }
                    }
                    else {
                      int rank,bits,index;
                      BitUtils::unpackRankBitsIndex(rank,bits,index,rbi_g[icv_nbr-ncv]);
                      if (rank == mpi_rank) {
                        assert((index >= 0)&&(index < ncv));
                        if ((flag_cv[index]&PART_MASK_BITS) == partVec.size()) {
                          const double factor = delta_vd[index]/delta_vd[icv];
                          if ((factor > 1.5)||(factor < 0.75))
                            smooth_flag[index] |= (1 << SMOOTH_FLAG_T_OFFSET);
                        }
                      }
                      else {
                        ++send_count[rank];
                      }
                    }
                  }
                  else if (icv_nbr >= ncv) {
                    int rank,bits,index;
                    BitUtils::unpackRankBitsIndex(rank,bits,index,rbi_g[icv_nbr-ncv]);
                    if (rank != mpi_rank) {
                      send_buf_double[send_disp[rank]] = delta_vd[icv];
                      send_buf_int[send_disp[rank]] = index;
                      ++send_disp[rank];
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
              send_buf_double = new double[send_count_sum];
              send_buf_int = new int[send_count_sum];
            }
          }

          // set up recv-side stuff...

          MPI_Alltoall(send_count,1,MPI_INT,recv_count,1,MPI_INT,mpi_comm);

          recv_disp[0] = 0;
          for (int rank = 1; rank < mpi_size; ++rank)
            recv_disp[rank] = recv_count[rank-1] + recv_disp[rank-1];
          const int recv_count_sum = recv_disp[mpi_size-1] + recv_count[mpi_size-1];

          int * recv_buf_int = new int[recv_count_sum];
          MPI_Alltoallv(send_buf_int,send_count,send_disp,MPI_INT,
              recv_buf_int,recv_count,recv_disp,MPI_INT,mpi_comm);
          delete[] send_buf_int;

          double * recv_buf_double = new double[recv_count_sum];
          MPI_Alltoallv(send_buf_double,send_count,send_disp,MPI_DOUBLE,
              recv_buf_double,recv_count,recv_disp,MPI_DOUBLE,mpi_comm);
          delete[] send_buf_double;

          for (int irecv = 0; irecv < recv_count_sum; ++irecv) {
            const int icv = recv_buf_int[irecv];
            assert((icv >= 0)&&(icv < ncv));
            if ((flag_cv[icv]&PART_MASK_BITS) == partVec.size()) {
              const double delta = recv_buf_double[irecv];
              const double factor = delta/delta_vd[icv];
              if ((factor > 1.5)||(factor < 0.75))
                smooth_flag[icv] |= (1 << SMOOTH_FLAG_T_OFFSET);
            }
          }
          delete[] recv_buf_int;
          delete[] recv_buf_double;
        }

        // now set nbrs of flag == 1 cvs with flag == 2 for each info type

        uint layer_max = 0;
        FOR_I4 layer_max = max(layer_max,smooth_mode_idata[i]);
        if (mpi_rank == 0) cout << " > growing layers" << endl;

        for (uint layer = 1; layer < layer_max; ++layer) {

          // 1. set local nbrs and store off processor nbrs...
          // here we loop through the parts because we only grow layers within the same part...

          for (int ipart = 0; ipart <= partVec.size(); ++ipart) {

            assert(rbiSet.empty());
            if (layer < n_blayers) {
              for (int ii = 0; ii < vdReturnDataVec.size(); ++ii) {
                const int icv = vdReturnDataVec[ii].icv;
                assert((icv >= 0)&&(icv < ncv));
                const uint flag_layer = ((smooth_flag[icv] >> SMOOTH_FLAG_S_OFFSET)&SMOOTH_FLAG_MASK);
                if ((flag_layer == layer) && ((flag_cv[icv]&PART_MASK_BITS) == ipart)) {
                  for (int ifa = 0; ifa < vdReturnDataVec[ii].nfa; ++ifa) {
                    const int icv_nbr = vdReturnDataVec[ii].cvofa[ifa];
                    if (icv_nbr < ncv) {
                      const uint flag_nbr = ((smooth_flag[icv_nbr] >> SMOOTH_FLAG_S_OFFSET)&SMOOTH_FLAG_MASK);
                      if ((flag_nbr == 0)&&((flag_cv[icv_nbr]&PART_MASK_BITS) == ipart))
                        smooth_flag[icv_nbr] |= (layer+1)<<SMOOTH_FLAG_S_OFFSET;
                    }
                    else {
                      int rank,bits,index;
                      BitUtils::unpackRankBitsIndex(rank,bits,index,rbi_g[icv_nbr-ncv]);
                      if (rank == mpi_rank) {
                        // this is local...
                        assert((index >= 0)&&(index < ncv));
                        const uint flag_nbr = ((smooth_flag[index] >> SMOOTH_FLAG_S_OFFSET)&SMOOTH_FLAG_MASK);
                        if ((flag_nbr == 0)&&((flag_cv[index]&PART_MASK_BITS) == ipart))
                          smooth_flag[index] |= (layer+1)<<SMOOTH_FLAG_S_OFFSET;
                      }
                      else {
                        rbiSet.insert(rbi_g[icv_nbr-ncv]);
                      }
                    }
                  }
                }
              }
            }

            if (layer < n_tlayers) {
              for (int ii = 0; ii < vdReturnDataVec.size(); ++ii) {
                const int icv = vdReturnDataVec[ii].icv;
                assert((icv >= 0)&&(icv < ncv));
                const uint flag_layer = ((smooth_flag[icv] >> SMOOTH_FLAG_T_OFFSET)&SMOOTH_FLAG_MASK);
                if ((flag_layer == layer) && ((flag_cv[icv]&PART_MASK_BITS) == ipart)) {
                  for (int ifa = 0; ifa < vdReturnDataVec[ii].nfa; ++ifa) {
                    const int icv_nbr = vdReturnDataVec[ii].cvofa[ifa];
                    if (icv_nbr < ncv) {
                      const uint flag_nbr = ((smooth_flag[icv_nbr] >> SMOOTH_FLAG_T_OFFSET)&SMOOTH_FLAG_MASK);
                      if ((flag_nbr == 0)&&((flag_cv[icv_nbr]&PART_MASK_BITS) == ipart))
                        smooth_flag[icv_nbr] |= (layer+1)<<SMOOTH_FLAG_T_OFFSET;
                    }
                    else {
                      int rank,bits,index;
                      BitUtils::unpackRankBitsIndex(rank,bits,index,rbi_g[icv_nbr-ncv]);
                      if (rank == mpi_rank) {
                        // this is local...
                        assert((index >= 0)&&(index < ncv));
                        const uint flag_nbr = ((smooth_flag[index] >> SMOOTH_FLAG_T_OFFSET)&SMOOTH_FLAG_MASK);
                        if ((flag_nbr == 0)&&((flag_cv[index]&PART_MASK_BITS) == ipart))
                          smooth_flag[index] |= (layer+1)<<SMOOTH_FLAG_T_OFFSET;
                      }
                      else {
                        rbiSet.insert(rbi_g[icv_nbr-ncv]);
                      }
                    }
                  }
                }
              }
            }

            if (layer < n_players) {
              for (int ii = 0; ii < vdReturnDataVec.size(); ++ii) {
                const int icv = vdReturnDataVec[ii].icv;
                assert((icv >= 0)&&(icv < ncv));
                const uint flag_layer = ((smooth_flag[icv] >> SMOOTH_FLAG_P_OFFSET)&SMOOTH_FLAG_MASK);
                if ((flag_layer == layer) && ((flag_cv[icv]&PART_MASK_BITS) == ipart)) {
                  for (int ifa = 0; ifa < vdReturnDataVec[ii].nfa; ++ifa) {
                    const int icv_nbr = vdReturnDataVec[ii].cvofa[ifa];
                    if (icv_nbr < ncv) {
                      const uint flag_nbr = ((smooth_flag[icv_nbr] >> SMOOTH_FLAG_P_OFFSET)&SMOOTH_FLAG_MASK);
                      if ((flag_nbr == 0)&&((flag_cv[icv_nbr]&PART_MASK_BITS) == ipart))
                        smooth_flag[icv_nbr] |= (layer+1)<<SMOOTH_FLAG_P_OFFSET;
                    }
                    else {
                      int rank,bits,index;
                      BitUtils::unpackRankBitsIndex(rank,bits,index,rbi_g[icv_nbr-ncv]);
                      if (rank == mpi_rank) {
                        // this is local...
                        assert((index >= 0)&&(index < ncv));
                        const uint flag_nbr = ((smooth_flag[index] >> SMOOTH_FLAG_P_OFFSET)&SMOOTH_FLAG_MASK);
                        if ((flag_nbr == 0)&&((flag_cv[index]&PART_MASK_BITS) == ipart))
                          smooth_flag[index] |= (layer+1)<<SMOOTH_FLAG_P_OFFSET;
                      }
                      else {
                        rbiSet.insert(rbi_g[icv_nbr-ncv]);
                      }
                    }
                  }
                }
              }
            }

            if (layer < n_pblayers) {
              for (int ii = 0; ii < vdReturnDataVec.size(); ++ii) {
                const int icv = vdReturnDataVec[ii].icv;
                assert((icv >= 0)&&(icv < ncv));
                const uint flag_layer = ((smooth_flag[icv] >> SMOOTH_FLAG_PB_OFFSET)&SMOOTH_FLAG_MASK);
                if ((flag_layer == layer) && ((flag_cv[icv]&PART_MASK_BITS) == ipart)) {
                  for (int ifa = 0; ifa < vdReturnDataVec[ii].nfa; ++ifa) {
                    const int icv_nbr = vdReturnDataVec[ii].cvofa[ifa];
                    if (icv_nbr < ncv) {
                      const uint flag_nbr = ((smooth_flag[icv_nbr] >> SMOOTH_FLAG_PB_OFFSET)&SMOOTH_FLAG_MASK);
                      if ((flag_nbr == 0)&&((flag_cv[icv_nbr]&PART_MASK_BITS) == ipart))
                        smooth_flag[icv_nbr] |= (layer+1)<<SMOOTH_FLAG_PB_OFFSET;
                    }
                    else {
                      int rank,bits,index;
                      BitUtils::unpackRankBitsIndex(rank,bits,index,rbi_g[icv_nbr-ncv]);
                      if (rank == mpi_rank) {
                        // this is local...
                        assert((index >= 0)&&(index < ncv));
                        const uint flag_nbr = ((smooth_flag[index] >> SMOOTH_FLAG_PB_OFFSET)&SMOOTH_FLAG_MASK);
                        if ((flag_nbr == 0)&&((flag_cv[index]&PART_MASK_BITS) == ipart))
                          smooth_flag[index] |= (layer+1)<<SMOOTH_FLAG_PB_OFFSET;
                      }
                      else {
                        rbiSet.insert(rbi_g[icv_nbr-ncv]);
                      }
                    }
                  }
                }
              }
            }

            // it is possible that rbiSet is empty everywhere. If this is the case, this avoids
            // the all-2-all communication...

            uint8 my_count = rbiSet.size();
            uint8 count;
            MPI_Allreduce(&my_count,&count,1,MPI_UINT8,MPI_SUM,mpi_comm);
            if (count == 0)
              continue;

            // because the set is already in rank-order, we can immediately pack and
            // count at the same time...

            FOR_RANK send_count[rank] = 0;

            int * send_buf_int = new int[rbiSet.size()];
            int isend = 0;
            for (set<uint8>::iterator iter = rbiSet.begin(); iter != rbiSet.end(); ++iter) {
              int rank,bits,index;
              BitUtils::unpackRankBitsIndex(rank,bits,index,*iter);
              assert(rank != mpi_rank);
              send_buf_int[isend++] = index;
              ++send_count[rank];
            }

            send_disp[0] = 0;
            for (int rank = 1; rank < mpi_size; ++rank)
              send_disp[rank] = send_count[rank-1] + send_disp[rank-1];
            assert(rbiSet.size() == send_disp[mpi_size-1] + send_count[mpi_size-1]);
            rbiSet.clear();

            // set up recv-side stuff...

            MPI_Alltoall(send_count,1,MPI_INT,recv_count,1,MPI_INT,mpi_comm);

            recv_disp[0] = 0;
            for (int rank = 1; rank < mpi_size; ++rank)
              recv_disp[rank] = recv_count[rank-1] + recv_disp[rank-1];
            const int recv_count_sum = recv_disp[mpi_size-1] + recv_count[mpi_size-1];

            int * recv_buf_int = new int[recv_count_sum];
            MPI_Alltoallv(send_buf_int,send_count,send_disp,MPI_INT,
                recv_buf_int,recv_count,recv_disp,MPI_INT,mpi_comm);
            delete[] send_buf_int;

            for (int irecv = 0; irecv < recv_count_sum; ++irecv) {
              const int index = recv_buf_int[irecv];
              assert((index >= 0)&&(index < ncv));
              uint flag_tmp = ((smooth_flag[index] >> SMOOTH_FLAG_S_OFFSET)&SMOOTH_FLAG_MASK);
              if ((flag_tmp == 0)&&((flag_cv[index]&PART_MASK_BITS) == ipart))
                smooth_flag[index] |= (layer+1)<<SMOOTH_FLAG_S_OFFSET;

              flag_tmp = ((smooth_flag[index] >> SMOOTH_FLAG_T_OFFSET)&SMOOTH_FLAG_MASK);
              if ((flag_tmp == 0)&&((flag_cv[index]&PART_MASK_BITS) == ipart))
                smooth_flag[index] |= (layer+1)<<SMOOTH_FLAG_T_OFFSET;

              flag_tmp = ((smooth_flag[index] >> SMOOTH_FLAG_P_OFFSET)&SMOOTH_FLAG_MASK);
              if ((flag_tmp == 0)&&((flag_cv[index]&PART_MASK_BITS) == ipart))
                smooth_flag[index] |= (layer+1)<<SMOOTH_FLAG_P_OFFSET;

              flag_tmp = ((smooth_flag[index] >> SMOOTH_FLAG_PB_OFFSET)&SMOOTH_FLAG_MASK);
              if ((flag_tmp == 0)&&((flag_cv[index]&PART_MASK_BITS) == ipart))
                smooth_flag[index] |= (layer+1)<<SMOOTH_FLAG_PB_OFFSET;
            }
            delete[] recv_buf_int;

          }

        }
      }

    }

    bool fuseOrphanChunksTol() {

      MPI_Sync("fuseOrphanChunksTol start");

      // tolerance-base fuse...

      assert(vdReturnDataVec.size() == ncv);

      const double d2_tol = 1.0E-20; // 1.0E-16 was too large
      const double angle_tol = 1.0*M_PI/180.0; // 1 degree tolerance
      const double dp_tol = 1.0-cos(angle_tol);

      int * no_flag = NULL;
      int nno_size = 0;
      double (*x_no)[3] = NULL;

      int (*nooed)[2] = NULL;
      int (*faoed)[2] = NULL;
      int *edono_v = NULL;
      int ned_size = 0;

      int *faofa0 = NULL;
      int nfa_size = 0;

      // augment ghosts with orphan neighbors...

      for (int iocd = 0; iocd < ocdVec.size(); ++iocd) {
        int rank,bits,index;
        BitUtils::unpackRankBitsIndex(rank,bits,index,ocdVec[iocd].rbi);
        if ((bits != 0)||(rank != mpi_rank)) {
          map<const uint8,int>::iterator it = rbiMap.find(ocdVec[iocd].rbi);
          if (it == rbiMap.end()) {
            rbiMap[ocdVec[iocd].rbi] = rbi_g.size() + ncv;
            rbi_g.push_back(ocdVec[iocd].rbi);
          }
        }
        for (int ied = 0; ied < ocdVec[iocd].ned; ++ied) {
          FOR_I2 {
            const int8 ifa = ocdVec[iocd].faoed[ied][i];
            if (ifa < 0) {
              // this should be a -1 index'd rbi...
              int rank,bits,index;
              uint8 rbi_nbr = -ifa-1;
              BitUtils::unpackRankBitsIndex(rank,bits,index,rbi_nbr);
              //assert(bits == 0);
              if (rank >= mpi_size) {
                rank -= mpi_size;
                rbi_nbr = BitUtils::packRankBitsIndex(rank,bits,index);
              }
              if ((bits != 0)||(rank != mpi_rank)) {
                map<const uint8,int>::iterator it = rbiMap.find(rbi_nbr);
                if (it == rbiMap.end()) {
                  rbiMap[rbi_nbr] = rbi_g.size() + ncv;
                  rbi_g.push_back(rbi_nbr);
                }
              }
            }
          }
        }
      }

      // get ghost coordinates...

      int8 * cvora = NULL;
      MiscUtils::buildXora(cvora,ncv);
      assert(cvora[mpi_rank+1]-cvora[mpi_rank] == ncv);
      const int ncv_g_only = rbi_g.size();
      int8 * icv_global_g = new int8[ncv_g_only];
      for (int icv_g = 0; icv_g < ncv_g_only; ++icv_g) {
        int rank,bits,index;
        BitUtils::unpackRankBitsIndex(rank,bits,index,rbi_g[icv_g]);
        icv_global_g[icv_g] = cvora[rank] + index;
      }
      DistributedDataExchanger * dde = new DistributedDataExchanger(icv_global_g,ncv_g_only,cvora);
      delete[] icv_global_g;
      delete[] cvora;
      double (*x_vd_g)[3] = new double[ncv_g_only][3];
      dde->pull(x_vd_g,x_vd);
      delete dde;
      for (int icv_g = 0; icv_g < ncv_g_only; ++icv_g) {
        int rank,bits,index;
        BitUtils::unpackRankBitsIndex(rank,bits,index,rbi_g[icv_g]);
        if (bits) PeriodicData::periodicTranslate(x_vd_g[icv_g],1,bits);
      }

      map<const pair<int,int>,int> stMap;
      map<const pair<int,int>,int> cvMap;

      double my_d2_close_min = HUGE_VAL;
      double my_d2_max = 0.0;

      const int debug_rank = getIntParam("DEBUG_FUSE_RANK",-1);
      const int debug_icv = getIntParam("DEBUG_FUSE_ICV",-1);

      int my_ierr = 0;
      for (int ii = 0; ii < vdReturnDataVec.size(); ++ii) {

        const int icv = vdReturnDataVec[ii].icv;
        assert((icv >= 0)&&(icv < ncv));

        bool debug = (mpi_rank == debug_rank && icv == debug_icv);

        if (debug) {
          cout << "WRITING ii: " << ii << " corresponding to icv: " << vdReturnDataVec[ii].icv << " " << COUT_VEC(x_vd[icv]) << endl;
          CuttableVoronoiData cvd_tmp;
          vdReturnDataVec[ii].extractGroupAsCvd(cvd_tmp,0);
          cvd_tmp.writeTecplot(icv,x_vd[icv]);
          // this guy has some orphan chunks. Join them up...
          int iocd = vdReturnDataVec[ii].orphan_chunk_data;
          while (iocd != -1) {
            cout << "WRITING orphan chunk: " << iocd << endl;
            char filename[128];
            sprintf(filename,"ocd.%06d.dat",iocd);
            ocdVec[iocd].writeTecplot(filename,x_vd[icv]);
            iocd = ocdVec[iocd].next;
          }
        }

        //cout << "working on ii: " << ii << " vdReturnDataVec[ii].icv: " << vdReturnDataVec[ii].icv << " vdReturnDataVec[ii].orphan_chunk_data: " << vdReturnDataVec[ii].orphan_chunk_data << endl;

        int nno_max = 0;
        int ned_max = 0;
        if (vdReturnDataVec[ii].ngr == 1) {
          // if there is just one group, then the nno_max is just nno...
          nno_max = vdReturnDataVec[ii].nno;
          ned_max = vdReturnDataVec[ii].ned;
          assert(vdReturnDataVec[ii].edogr_i[1] == vdReturnDataVec[ii].ned);
        }
        else {
          // otherwise we need to cycle through first group of edges to get the
          // last node in the first group...
          for (int ied = 0; ied < vdReturnDataVec[ii].edogr_i[1]; ++ied) {
            FOR_I2 nno_max = max(nno_max,vdReturnDataVec[ii].nooed[ied][i]);
          }
          nno_max += 1;
          ned_max = vdReturnDataVec[ii].edogr_i[1];
        }
        const int nno_main = nno_max;
        int iocd = vdReturnDataVec[ii].orphan_chunk_data;
        while (iocd != -1) {
          nno_max += ocdVec[iocd].nno;
          ned_max += ocdVec[iocd].ned;
          iocd = ocdVec[iocd].next;
        }

        // reallocate memory if req'd...

        if (nno_max > nno_size) {
          nno_size = nno_max;
          if (no_flag) delete[] no_flag;
          no_flag = new int[nno_size*2]; // needed for several things
          if (x_no) delete[] x_no;
          x_no = new double[nno_size][3];
        }

        if (ned_max > ned_size) {
          ned_size = ned_max;
          if (nooed) delete[] nooed;
          nooed = new int[ned_size][2];
          if (faoed) delete[] faoed;
          faoed = new int[ned_size][2];
          if (edono_v) delete[] edono_v;
          edono_v = new int[ned_size*2];
        }

        for (int ino = 0; ino < nno_max; ++ino)
          no_flag[ino] = ino;

        // only copy the group 0 nodes...
        for (int ino = 0; ino < nno_main; ++ino) {
          FOR_I3 x_no[ino][i] =  vdReturnDataVec[ii].x_no[ino][i];
        }

        // this is N^2 for now...
        double d2_close = HUGE_VAL;
        double d2_max = 0.0;

        set<int> nodes;
        for (int ino = 0; ino < nno_main; ++ino) {
          // cycle through ino_new nodes and check for closeness...
          for (int ino_prev = 0; ino_prev < ino; ++ino_prev) {
            const double d2 = DIST2(x_no[ino_prev],x_no[ino]);
            if (d2 < delta_vd[icv]*delta_vd[icv]*d2_tol) {
              d2_max = max(d2_max,d2);
              // this node ino corresponds to ino_prev...
              if (no_flag[ino] != ino) {
                // node ino has already been set to node no_flag[ino]...
                assert(nodes.empty());
                nodes.insert(ino);
                int ino_set = ino;
                while (ino_set != no_flag[ino_set]) {
                  ino_set = no_flag[ino_set];
                  nodes.insert(ino_set);
                }
                nodes.insert(ino_prev);
                ino_set = ino_prev;
                while (ino_set != no_flag[ino_set]) {
                  ino_set = no_flag[ino_set];
                  nodes.insert(ino_set);
                }
                int ino_last = -1;
                for (set<int>::iterator it = nodes.begin(); it != nodes.end(); ++it) {
                  const int ino = *it;
                  if (ino_last == -1) {
                    assert(no_flag[ino] == ino);
                  }
                  else {
                    no_flag[ino] = ino_last;
                  }
                  ino_last = ino;
                }
                nodes.clear();
              }
              else {
                no_flag[ino] = ino_prev;
              }
            }
            else {
              d2_close = min(d2_close,d2);
            }
          }
        }
        iocd = vdReturnDataVec[ii].orphan_chunk_data;
        int nno_offset = nno_main;
        while (iocd != -1) {
          for (int ino = 0; ino < ocdVec[iocd].nno; ++ino) {
            FOR_I3 x_no[nno_offset+ino][i] = ocdVec[iocd].x_no[ino][i];
            for (int ino_prev = 0; ino_prev < nno_offset+ino; ++ino_prev) {
              const double d2 = DIST2(x_no[ino_prev],x_no[nno_offset+ino]);
              if (d2 < delta_vd[icv]*delta_vd[icv]*d2_tol) {
                d2_max = max(d2_max,d2);
                // this node ino corresponds to ino_prev...
                if (no_flag[nno_offset+ino] != nno_offset+ino) {
                  // node ino has already been set to node no_flag[ino]...
                  assert(nodes.empty());
                  nodes.insert(nno_offset+ino);
                  int ino_set = nno_offset+ino;
                  while (ino_set != no_flag[ino_set]) {
                    ino_set = no_flag[ino_set];
                    nodes.insert(ino_set);
                  }
                  nodes.insert(ino_prev);
                  ino_set = ino_prev;
                  while (ino_set != no_flag[ino_set]) {
                    ino_set = no_flag[ino_set];
                    nodes.insert(ino_set);
                  }
                  int ino_last = -1;
                  for (set<int>::iterator it = nodes.begin(); it != nodes.end(); ++it) {
                    const int ino = *it;
                    if (ino_last == -1) {
                      assert(no_flag[ino] == ino);
                    }
                    else {
                      no_flag[ino] = ino_last;
                    }
                    ino_last = ino;
                  }
                  nodes.clear();
                }
                else {
                  no_flag[nno_offset+ino] = ino_prev;
                }
              }
              else {
                d2_close = min(d2_close,d2);
              }
            }
          }
          nno_offset += ocdVec[iocd].nno;
          iocd = ocdVec[iocd].next;
        }
        assert(nno_offset == nno_max);

        // update the d2's...

        my_d2_close_min = min(my_d2_close_min,d2_close/(delta_vd[icv]*delta_vd[icv]));
        my_d2_max = max(my_d2_max,d2_max/(delta_vd[icv]*delta_vd[icv]));

        int nno_new = 0;
        for (int ino = 0; ino < nno_max; ++ino) {
          if (no_flag[ino] == ino) {
            const int ino_new = nno_new++;
            FOR_I3 x_no[ino_new][i] = x_no[ino][i];
            no_flag[ino] = -nno_new;
          }
          else {
            int ino_ = no_flag[ino];
            while (ino_ >= 0)
              ino_ = no_flag[ino_];
            no_flag[ino] = ino_;
          }
        }

        // no_flag is now -1 indexed with the new node index...

        // now copy the edges...

        assert(stMap.empty());
        int nbf_new = 0;
        assert(cvMap.empty());
        int nfa_new = 0;

        for (int ied = 0; ied < vdReturnDataVec[ii].edogr_i[1]; ++ied) {
          FOR_I2 {
            const int ino = vdReturnDataVec[ii].nooed[ied][i];
            nooed[ied][i] = -no_flag[ino]-1;
            assert((nooed[ied][i] >= 0)&&(nooed[ied][i] < nno_new));
          }
          // for the face, start a new indexing for the remaining edges and the surfaces (positive)
          // or local cvs (negative) they reference...
          FOR_I2 {
            const int ifa = vdReturnDataVec[ii].faoed[ied][i];
            if (ifa >= 0) {
              // a positive face references surface tri geometry...
              const int ist = vdReturnDataVec[ii].spbobf[ifa][0];
              const int ipart_bits = vdReturnDataVec[ii].spbobf[ifa][1];
              map<const pair<int,int>,int>::iterator it2 = stMap.find(pair<int,int>(ipart_bits,ist));
              if (it2 == stMap.end()) {
                stMap[pair<int,int>(ipart_bits,ist)] = faoed[ied][i] = nbf_new++;
              }
              else {
                faoed[ied][i] = it2->second;
              }
            }
            else {
              const int icv_nbr = vdReturnDataVec[ii].cvofa[-ifa-1];
              assert((icv_nbr >= 0)&&(icv_nbr < (rbi_g.size()+ncv)));
              // this internal plane is associated with a pair of cvs...
              pair<int,int> cv_pair = pair<int,int>(min(icv,icv_nbr),max(icv,icv_nbr));
              map<const pair<int,int>,int>::iterator it2 = cvMap.find(cv_pair);
              if (it2 == cvMap.end()) {
                cvMap[cv_pair] = nfa_new++;
                if (debug) {
                  cout << " fp pair " << COUT_VEC(x_vd[icv]);
                  if (icv_nbr < ncv_g)
                    cout << " internal: " << COUT_VEC(x_vd[icv_nbr]) << endl;
                  else
                    cout << " ghost: " << COUT_VEC(x_vd_g[icv_nbr-ncv]) << endl;
                }
                faoed[ied][i] = -nfa_new;
              }
              else {
                faoed[ied][i] = -it2->second-1;
              }
            }
          }
        }
        nno_offset = nno_main;
        int ned_offset = vdReturnDataVec[ii].edogr_i[1];
        iocd = vdReturnDataVec[ii].orphan_chunk_data;
        while (iocd != -1) {
          for (int ied = 0; ied < ocdVec[iocd].ned; ++ied) {
            FOR_I2 {
              const int ino = ocdVec[iocd].nooed[ied][i];
              nooed[ned_offset+ied][i] = -no_flag[nno_offset+ino]-1;
              assert((nooed[ned_offset+ied][i] >= 0)&&(nooed[ned_offset+ied][i] < nno_new));
            }
            FOR_I2 {
              const int8 ifa = ocdVec[iocd].faoed[ied][i];
              if (ifa >= 0) {
                // this is a ipart_bits and ist...
                const int ipart_bits = int(ifa>>32);
                const int ist = int(ifa&MASK_32BITS);
                map<const pair<int,int>,int>::iterator it2 = stMap.find(pair<int,int>(ipart_bits,ist));
                if (it2 == stMap.end()) {
                  stMap[pair<int,int>(ipart_bits,ist)] = faoed[ned_offset+ied][i] = nbf_new++;
                }
                else {
                  faoed[ned_offset+ied][i] = it2->second;
                }
              }
              else {
                // this should be a -1 index'd rbi...
                int rank,bits,index;
                uint8 rbi_nbr = -ifa-1;
                BitUtils::unpackRankBitsIndex(rank,bits,index,rbi_nbr);
                //assert(bits == 0);
                //bool b_orphan_cut;
                int icv0;
                if (rank >= mpi_size) {
                  rank -= mpi_size;
                  //b_orphan_cut = true;
                  icv0 = icv;
                  rbi_nbr = BitUtils::packRankBitsIndex(rank,bits,index);
                }
                else {
                  //b_orphan_cut = false;
                  int rank0,bits0,index0;
                  BitUtils::unpackRankBitsIndex(rank0,bits0,index0,ocdVec[iocd].rbi);
                  if ((bits0 == 0)&&(rank0 == mpi_rank)) {
                    icv0 = index0;
                  }
                  else {
                    map<const uint8,int>::iterator it = rbiMap.find(ocdVec[iocd].rbi);
                    assert(it != rbiMap.end());
                    icv0 = it->second;
                  }
                }
                int icv_nbr;
                if ((bits == 0)&&(rank == mpi_rank)) {
                  assert((index >= 0)&&(index < ncv));
                  icv_nbr = index;
                }
                else {
                  map<const uint8,int>::iterator it = rbiMap.find(rbi_nbr);
                  assert(it != rbiMap.end());
                  icv_nbr = it->second;
                }

                pair<int,int> cv_pair = pair<int,int>(min(icv0,icv_nbr),max(icv0,icv_nbr));
                map<const pair<int,int>,int>::iterator it = cvMap.find(cv_pair);
                if (it == cvMap.end()) {
                  cvMap[cv_pair] = nfa_new++;
                  if (debug) {
                    cout << " or pair ";
                    if (icv0 < ncv)
                      cout << " internal: " << COUT_VEC(x_vd[icv0]);
                    else
                      cout << " ghost: " << COUT_VEC(x_vd_g[icv0-ncv]);
                    if (icv_nbr < ncv_g)
                      cout << " internal: " << COUT_VEC(x_vd[icv_nbr]) << endl;
                    else
                      cout << " ghost: " << COUT_VEC(x_vd_g[icv_nbr-ncv]) << endl;
                  }
                  faoed[ned_offset+ied][i] = -nfa_new;
                }
                else {
                  faoed[ned_offset+ied][i] = -it->second-1;
                }

              }
            }
          }
          nno_offset += ocdVec[iocd].nno;
          ned_offset += ocdVec[iocd].ned;
          const int iocd_next = ocdVec[iocd].next;
          ocdVec[iocd].clear();
          iocd = iocd_next;
        }
        assert(ned_offset == ned_max);

        // now look for edge pairs to eliminate...

        int * edono_i = no_flag;

        for (int ino = 0; ino < nno_new; ++ino)
          edono_i[ino+1] = 0;

        for (int ied = 0; ied < ned_max; ++ied) {
          const int ino0 = nooed[ied][0]; assert((ino0 >= 0)&&(ino0 < nno_new));
          const int ino1 = nooed[ied][1]; assert((ino1 >= 0)&&(ino1 < nno_new));
          if (ino0 == ino1) {
            // if the 2 nodes are identical, the edge can be eliminated right away...
            nooed[ied][0] = -1;
            nooed[ied][1] = -1;
          }
          else {
            ++edono_i[ino0+1];
            ++edono_i[ino1+1];
          }
        }

        edono_i[0] = 0;
        for (int ino = 0; ino < nno_new; ++ino)
          edono_i[ino+1] += edono_i[ino];

        const int edono_s = edono_i[nno_new];
        assert(edono_s <= ned_max*2);

        for (int ied = 0; ied < ned_max; ++ied) {
          const int ino0 = nooed[ied][0];
          const int ino1 = nooed[ied][1];
          if ((ino0 >= 0)&&(ino1 >= 0)) {
            edono_v[edono_i[ino0]++] = ied;
            edono_v[edono_i[ino1]++] = ied;
          }
        }

        for (int ino = nno_new-1; ino > 0; --ino)
          edono_i[ino] = edono_i[ino-1];
        edono_i[0] = 0;

        // now look for edge pair matches...

        for (int ied = 0; ied < ned_max; ++ied) {
          const int ino0 = nooed[ied][0];
          if (ino0 >= 0) {
            const int ino1 = nooed[ied][1];
            assert(ino1 >= 0);
            for (int eon = edono_i[ino0]; eon != edono_i[ino0+1]; ++eon) {
              const int ied_nbr = edono_v[eon];
              if ((ied_nbr > ied)&&(nooed[ied_nbr][0] >= 0)) {
                assert(nooed[ied_nbr][1] >= 0);
                if (nooed[ied_nbr][0] == ino0) {
                  if (nooed[ied_nbr][1] == ino1) {
                    if ((faoed[ied][0] == faoed[ied_nbr][1])&&(faoed[ied][1] == faoed[ied_nbr][0])) {
                      nooed[ied][0] = nooed[ied][1] = nooed[ied_nbr][0] = nooed[ied_nbr][1] = -1;
                    }
                    else if (faoed[ied][1] == faoed[ied_nbr][0]) {
                      // transfer info from ied_nbr to ied...
                      faoed[ied][1] = faoed[ied_nbr][1];
                      // and eliminate ied_nbr...
                      nooed[ied_nbr][0] = nooed[ied_nbr][1] = -1;
                    }
                    else if (faoed[ied][0] == faoed[ied_nbr][1]) {
                      // transfer info from ied_nbr to ied...
                      faoed[ied][0] = faoed[ied_nbr][0];
                      // and eliminate ied_nbr...
                      nooed[ied_nbr][0] = nooed[ied_nbr][1] = -1;
                    }
                    //else {
                    //  cout << "F2: faoed[ied][0,1]: " << faoed[ied][0] << " " << faoed[ied][1] << " faoed[ied_nbr][1,0]: " << faoed[ied_nbr][1] << " " << faoed[ied_nbr][0] << endl;
                    //  if ((faoed[ied][0] < 0)&&(faoed[ied][1] < 0)) MiscUtils::writeEdge(ied,x_no[ino0],x_no[ino1]);
                    //  debug = true;
                    //}
                  }
                }
                else {
                  assert(nooed[ied_nbr][1] == ino0);
                  if (nooed[ied_nbr][0] == ino1) {
                    if ((faoed[ied][0] == faoed[ied_nbr][0])&&(faoed[ied][1] == faoed[ied_nbr][1])) {
                      nooed[ied][0] = nooed[ied][1] = nooed[ied_nbr][0] = nooed[ied_nbr][1] = -1;
                    }
                    else if (faoed[ied][0] == faoed[ied_nbr][0]) {
                      faoed[ied][0] = faoed[ied_nbr][1];
                      nooed[ied_nbr][0] = nooed[ied_nbr][1] = -1;
                    }
                    else if (faoed[ied][1] == faoed[ied_nbr][1]) {
                      faoed[ied][1] = faoed[ied_nbr][0];
                      nooed[ied_nbr][0] = nooed[ied_nbr][1] = -1;
                    }
                    //else {
                    //  cout << "R2: faoed[ied][0,1]: " << faoed[ied][0] << " " << faoed[ied][1] << " faoed[ied_nbr][0,1]: " << faoed[ied_nbr][0] << " " << faoed[ied_nbr][1] << endl;
                    //  if ((faoed[ied][0] < 0)&&(faoed[ied][1] < 0)) MiscUtils::writeEdge(ied,x_no[ino0],x_no[ino1]);
                    //  debug = true;
                    //}
                  }
                }
              }
            }
          }
        }


        // manage face mem...

        if (nfa_new > nfa_size) {
          nfa_size = nfa_new;
          if (faofa0) delete[] faofa0;
          faofa0 = new int[nfa_size];
        }

        // look for coplanar degeneracies...

        const int nfa_new0 = nfa_new;
        for (int ifa0 = 0; ifa0 < nfa_new0; ++ifa0)
          faofa0[ifa0] = ifa0;
        for (map<const pair<int,int>,int>::iterator it = cvMap.begin(); it != cvMap.end(); ++it) {
          if ((it->first.first == icv)||(it->first.second == icv))
            continue;
          double x_fa_or[3];
          double n_fa_or[3];
          if (it->first.first < ncv) {
            FOR_I3 n_fa_or[i] = -x_vd[it->first.first][i];
            FOR_I3 x_fa_or[i] = x_vd[it->first.first][i];
          }
          else {
            FOR_I3 n_fa_or[i] = -x_vd_g[it->first.first-ncv][i];
            FOR_I3 x_fa_or[i] = x_vd_g[it->first.first-ncv][i];
          }
          //if (icv == 4) cout << "or: " << COUT_VEC(x_fa_or) << endl;
          if (it->first.second < ncv) {
            FOR_I3 n_fa_or[i] += x_vd[it->first.second][i];
            FOR_I3 x_fa_or[i] += x_vd[it->first.second][i];
          }
          else {
            FOR_I3 n_fa_or[i] += x_vd_g[it->first.second-ncv][i];
            FOR_I3 x_fa_or[i] += x_vd_g[it->first.second-ncv][i];
          }
          FOR_I3 x_fa_or[i] *= 0.5;
          //if (icv == 4) cout << "or2: " << COUT_VEC(x_fa_or) << endl;
          // look for a coplanar foster parent face...
          for (map<const pair<int,int>,int>::iterator it2 = cvMap.begin(); it2 != cvMap.end(); ++it2) {
            if ((it2->first.first != icv)&&(it2->first.second != icv))
              continue;
            double n_fa_fp[3];
            double x_fa_fp[3];
            FOR_I3 {
              n_fa_fp[i] = -x_vd[icv][i];
              x_fa_fp[i] = x_vd[icv][i];
            }
            //if (icv == 4) cout << "fp: " << COUT_VEC(x_fa_fp) << endl;
            if (icv == it2->first.first) {
              if (it2->first.second < ncv) {
                FOR_I3 {
                  n_fa_fp[i] += x_vd[it2->first.second][i];
                  x_fa_fp[i] += x_vd[it2->first.second][i];
                }
              }
              else {
                FOR_I3 {
                  n_fa_fp[i] += x_vd_g[it2->first.second-ncv][i];
                  x_fa_fp[i] += x_vd_g[it2->first.second-ncv][i];
                }
              }
            }
            else {
              assert(icv == it2->first.second);
              if (it2->first.first < ncv) {
                FOR_I3 {
                  n_fa_fp[i] += x_vd[it2->first.first][i];
                  x_fa_fp[i] += x_vd[it2->first.first][i];
                }
              }
              else {
                FOR_I3 {
                  n_fa_fp[i] += x_vd_g[it2->first.first-ncv][i];
                  x_fa_fp[i] += x_vd_g[it2->first.first-ncv][i];
                }
              }
            }
            FOR_I3 x_fa_fp[i] *= 0.5;
            //if (icv == 4) cout << "fp2: " << COUT_VEC(x_fa_fp) << endl;

            const double dx[3] = DIFF(x_fa_or,x_fa_fp);
            const double phi = DOT_PRODUCT(n_fa_fp,dx);
            const double mag_n2 = DOT_PRODUCT(n_fa_fp,n_fa_fp);
            if (debug)
              cout << phi*phi/(mag_n2*delta_vd[icv]*delta_vd[icv]) << " <? " << d2_tol << " " << fabs(DOT_PRODUCT(n_fa_or,n_fa_fp))/(sqrt(mag_n2)*MAG(n_fa_or)) << " >? " << (1.0-dp_tol) << endl;
            if (phi*phi < mag_n2*delta_vd[icv]*delta_vd[icv]*d2_tol) {
              if (fabs(DOT_PRODUCT(n_fa_or,n_fa_fp)) >= (1.0-dp_tol)*sqrt(mag_n2)*MAG(n_fa_or)) {
                faofa0[it->second] = it2->second;
                --nfa_new;
                break;
              }
            }
          }
        }
        // reindex faces...
        int ifa_new = 0;
        for (int ifa0 = 0; ifa0 < nfa_new0; ++ifa0) {
          if (faofa0[ifa0] == ifa0) {
            ifa_new++;
            faofa0[ifa0] = -ifa_new;
          }
        }
        assert(ifa_new == nfa_new);
        for (map<const pair<int,int>,int>::iterator it = cvMap.begin(); it != cvMap.end(); ++it) {
          const int ifa0 = it->second;
          if (faofa0[ifa0] >= 0) {
            assert(faofa0[ifa0] != ifa0);
            int ifa = faofa0[ifa0];
            while (ifa >= 0)
              ifa = faofa0[ifa];
            faofa0[ifa0] = ifa;
            it->second = -1;
          }
          else {
            it->second = -faofa0[ifa0]-1; // update map
          }
          assert(faofa0[ifa0] < 0);
        }

        map<const pair<int,int>,int> edgeMap;
        for (int ied = 0; ied < ned_max; ++ied) {
          if (nooed[ied][0] == -1) {
            assert(nooed[ied][1] == -1);
            continue;
          }
          FOR_I2 {
            if (faoed[ied][i] < 0)
              faoed[ied][i] = faofa0[-faoed[ied][i]-1];
          }
          map<const pair<int,int>,int>::iterator it = edgeMap.find(pair<int,int>(nooed[ied][0],nooed[ied][1]));
          if (it != edgeMap.end()) {
            // I think this is a collapsed face, lets kill these edge pairs...
            const int ied_pair = it->second;
            if (((nooed[ied][0] == nooed[ied_pair][0])&&
                  (nooed[ied][1] == nooed[ied_pair][1])&&
                  (faoed[ied][0] == faoed[ied_pair][1])&&
                  (faoed[ied][1] == faoed[ied_pair][0]))||
                ((nooed[ied][0] == nooed[ied_pair][1])&&
                 (nooed[ied][1] == nooed[ied_pair][0])&&
                 (faoed[ied][0] == faoed[ied_pair][0])&&
                 (faoed[ied][1] == faoed[ied_pair][1])))
              nooed[ied][0] = nooed[ied][1] = nooed[ied_pair][0] = nooed[ied_pair][1] = -1;
          }
          else {
            edgeMap[pair<int,int>(nooed[ied][1],nooed[ied][0])] = ied;
            edgeMap[pair<int,int>(nooed[ied][0],nooed[ied][1])] = ied;
          }
          if (faoed[ied][0] == faoed[ied][1])
            nooed[ied][0] = nooed[ied][1] = -1;
        }
        edgeMap.clear();

        if (debug) {
          cout << " faces: " << endl;
          for (int ifa0 = 0; ifa0 < nfa_new0; ++ifa0)
            cout << ifa0 << " " << faofa0[ifa0] << " " << endl;
          cout << " edges: " << endl;
          for (int ied = 0; ied < ned_max; ++ied) {
            if (faoed[ied][0] == -2 || faoed[ied][1] == -2)
              cout << ied << " " << faoed[ied][0] << " " << faoed[ied][1] << " " << nooed[ied][0] << " " << nooed[ied][1] << endl;
          }
          cout.flush();
        }

        // at this point, eliminated edges have their nooed set to -1. There is still more
        // elimination to do. We need to detect and eliminate double-edge nodes. These nodes
        // cannot add any information to the vd and can be eliminated...

        /*
           for (int ino = 0; ino < nno_new*2; ++ino)
           no_flag[ino] = -1;

           for (int ied = 0; ied < ned_max; ++ied) {
           if (nooed[ied][0] >= 0) {
           assert(nooed[ied][1] >= 0);
           FOR_I2 {
           const int ino = nooed[ied][i];
           assert((ino >= 0)&&(ino < nno_new));
           if (no_flag[ino*2] == -1) {
        // the first position in no_flag is available to store the edge...
        no_flag[ino*2] = ied;
        }
        else if (no_flag[ino*2+1] == -1) {
        // the second position in no_flag is available to store the edge...
        no_flag[ino*2+1] = ied;
        }
        else {
        // there are more than 2 edges touching this node, so set to -2...
        no_flag[ino*2] = no_flag[ino*2+1] = -2;
        }
        }
        }
        }
        */
        for (int ino = 0; ino < nno_new*2; ++ino)
          no_flag[ino] = -2;

        int nno_new_final = 0;
        for (int ino = 0; ino < nno_new; ++ino) {
          if (no_flag[ino*2] == -2) {
            assert(no_flag[ino*2+1] == -2);
            // this node has 3 or more edges: keep, but re-index...
            const int ino_new = nno_new_final++;
            // use no_flag[ino] to store the new node index...
            no_flag[ino] = -ino_new-3;
            if (ino_new < ino) {
              FOR_I3 x_no[ino_new][i] = x_no[ino][i];
            }
          }
          else if (no_flag[ino*2] >= 0) {
            assert(no_flag[ino*2+1] >= 0);
            // remove the second edge, reconnecting with the first...
            const int ied0 = no_flag[ino*2];
            const int ied1 = no_flag[ino*2+1];
            if (nooed[ied0][0] == ino) {
              if (nooed[ied1][0] == ino) {
                assert(faoed[ied0][0] == faoed[ied1][1]);
                assert(faoed[ied0][1] == faoed[ied1][0]);
                nooed[ied0][0] = nooed[ied1][1];
                if (no_flag[nooed[ied1][1]*2] == ied1) {
                  no_flag[nooed[ied1][1]*2] = ied0;
                }
                else if (no_flag[nooed[ied1][1]*2+1] == ied1) {
                  no_flag[nooed[ied1][1]*2+1] = ied0;
                }
                nooed[ied1][0] = nooed[ied1][1] = -1;
              }
              else {
                assert(nooed[ied1][1] == ino);
                if (!(faoed[ied0][0] == faoed[ied1][0])) {
                  cout << " YYYYYY " << ii <<  " " << ied0 << " " << ied1 << endl;
                  cout << faoed[ied0][0] << " " << faoed[ied0][1] << endl;
                  cout << faoed[ied1][0] << " " << faoed[ied1][1] << endl;
                  cout.flush();
                }
                assert(faoed[ied0][0] == faoed[ied1][0]);
                assert(faoed[ied0][1] == faoed[ied1][1]);
                nooed[ied0][0] = nooed[ied1][0];
                if (no_flag[nooed[ied1][0]*2] == ied1) {
                  no_flag[nooed[ied1][0]*2] = ied0;
                }
                else if (no_flag[nooed[ied1][0]*2+1] == ied1) {
                  no_flag[nooed[ied1][0]*2+1] = ied0;
                }
                nooed[ied1][0] = nooed[ied1][1] = -1;
              }
            }
            else {
              assert(nooed[ied0][1] == ino);
              if (nooed[ied1][0] == ino) {
                if (!(faoed[ied0][0] == faoed[ied1][0])) {
                  cout << " YYYYYY " << ii <<  " " << ied0 << " " << ied1 << endl;
                  cout << faoed[ied0][0] << " " << faoed[ied0][1] << endl;
                  cout << faoed[ied1][0] << " " << faoed[ied1][1] << endl;
                  cout.flush();
                }
                assert(faoed[ied0][0] == faoed[ied1][0]);
                assert(faoed[ied0][1] == faoed[ied1][1]);
                nooed[ied0][1] = nooed[ied1][1];
                if (no_flag[nooed[ied1][1]*2] == ied1) {
                  no_flag[nooed[ied1][1]*2] = ied0;
                }
                else if (no_flag[nooed[ied1][1]*2+1] == ied1) {
                  no_flag[nooed[ied1][1]*2+1] = ied0;
                }
                nooed[ied1][0] = nooed[ied1][1] = -1;
              }
              else {
                assert(nooed[ied1][1] == ino);
                assert(faoed[ied0][0] == faoed[ied1][1]);
                assert(faoed[ied0][1] == faoed[ied1][0]);
                nooed[ied0][1] = nooed[ied1][0];
                if (no_flag[nooed[ied1][0]*2] == ied1) {
                  no_flag[nooed[ied1][0]*2] = ied0;
                }
                else if (no_flag[nooed[ied1][0]*2+1] == ied1) {
                  no_flag[nooed[ied1][0]*2+1] = ied0;
                }
                nooed[ied1][0] = nooed[ied1][1] = -1;
              }
            }
          }
          else {
            assert((no_flag[ino*2] == -1)&&(no_flag[ino*2+1] == -1));
          }
        }

        // reindex faces one more time...
        int nfa_new_final = 0;
        map<const int,int> faMap;
        for (int ied = 0; ied < ned_max; ++ied) {
          if ((nooed[ied][0] >= 0)&&(nooed[ied][1] >= 0)) {
            FOR_I2 {
              if (faoed[ied][i] < 0) {
                const int ifa = -faoed[ied][i]-1;
                map<const int,int>::iterator it = faMap.find(ifa);
                if (it == faMap.end()) {
                  faMap[ifa] = faofa0[ifa] = nfa_new_final++; // reuse faofa0 (+ is keeper)
                  faoed[ied][i] = -faofa0[ifa]-1;
                }
                else {
                  faoed[ied][i] = -it->second-1;
                }
              }
            }
          }
        }
        faMap.clear();
        assert(nfa_new_final <= nfa_new);

        int ned_new_final = 0;
        for (int ied = 0; ied < ned_max; ++ied) {
          if ((nooed[ied][0] >= 0)&&(nooed[ied][1] >= 0)) {
            const int ied_new = ned_new_final++;
            FOR_I2 {
              const int ino = nooed[ied][i];
              assert(no_flag[ino] <= -3);
              nooed[ied_new][i] = -no_flag[ino]-3;
            }
            faoed[ied_new][0] = faoed[ied][0];
            faoed[ied_new][1] = faoed[ied][1];
          }
          else {
            assert((nooed[ied][0] == -1)&&(nooed[ied][1] == -1));
          }
        }

        if (nno_new_final > vdReturnDataVec[ii].nno) {
          delete[] vdReturnDataVec[ii].x_no;
          vdReturnDataVec[ii].x_no = new double[nno_new_final][3];
        }
        vdReturnDataVec[ii].nno = nno_new_final;
        for (int ino = 0; ino < nno_new_final; ++ino)
          FOR_I3 vdReturnDataVec[ii].x_no[ino][i] = x_no[ino][i];

        if (ned_new_final > vdReturnDataVec[ii].ned) {
          delete[] vdReturnDataVec[ii].nooed;
          delete[] vdReturnDataVec[ii].faoed;
          vdReturnDataVec[ii].nooed = new int[ned_new_final][2];
          vdReturnDataVec[ii].faoed = new int[ned_new_final][2];
        }
        vdReturnDataVec[ii].ngr = 1;
        vdReturnDataVec[ii].ned = vdReturnDataVec[ii].edogr_i[1] = ned_new_final;
        for (int ied = 0; ied < ned_new_final; ++ied) {
          FOR_I2 vdReturnDataVec[ii].nooed[ied][i] = nooed[ied][i];
          FOR_I2 vdReturnDataVec[ii].faoed[ied][i] = faoed[ied][i];
        }

        // stMap...
        assert(nbf_new == stMap.size());
        if (nbf_new > vdReturnDataVec[ii].nbf) {
          delete[] vdReturnDataVec[ii].spbobf;
          vdReturnDataVec[ii].spbobf = new int[nbf_new][2];
        }
        vdReturnDataVec[ii].nbf = nbf_new;
        for (map<const pair<int,int>,int>::iterator it = stMap.begin(); it != stMap.end(); ++it) {
          const int ibf = it->second;
          assert((ibf >= 0)&&(ibf < nbf_new));
          vdReturnDataVec[ii].spbobf[ibf][0] = it->first.second; // ist first
          vdReturnDataVec[ii].spbobf[ibf][1] = it->first.first; // ipart_bits second
        }
        stMap.clear();

        if (nfa_new_final > vdReturnDataVec[ii].nfa) {
          delete[] vdReturnDataVec[ii].cvofa;
          vdReturnDataVec[ii].cvofa = new int[nfa_new_final];
        }
        vdReturnDataVec[ii].nfa = nfa_new_final;
        for (map<const pair<int,int>,int>::iterator it = cvMap.begin(); it != cvMap.end(); ++it) {
          if (it->second >= 0) {
            const int ifa0 = it->second; assert(ifa0 < nfa_new);
            const int ifa = faofa0[ifa0];
            if (ifa >= 0) {
              assert(ifa < nfa_new_final);
              if (it->first.first == vdReturnDataVec[ii].icv) {
                assert(it->first.second != vdReturnDataVec[ii].icv);
                vdReturnDataVec[ii].cvofa[ifa] = it->first.second;
              }
              else if (it->first.second == vdReturnDataVec[ii].icv) {
                vdReturnDataVec[ii].cvofa[ifa] = it->first.first;
              }
              else {
                //vdReturnDataVec[ii].cvofa[ifa] = -1; // this guy should never be called...
                vdReturnDataVec[ii].cvofa[ifa] = -it->first.first-1; // -1 idx'n to indicate face mismatch
              }
            }
          }
        }
        cvMap.clear();

        if (debug) {
          cout << COUT_VEC(x_vd[icv]) << endl;
          cout << " faces: " << endl;
          for (int ifa0 = 0; ifa0 < nfa_new0; ++ifa0) {
            const int ifa = faofa0[ifa0];
            cout << ifa0 << " " << ifa;
            if (ifa >= 0)
              cout << " " << vdReturnDataVec[ii].cvofa[ifa];
            cout << endl;
          }
          cout << " edges: " << endl;
          for (int ied = 0; ied < ned_new_final; ++ied) {
            cout << ied << " " << faoed[ied][0] << " " << faoed[ied][1] << " " << nooed[ied][0] << " " << nooed[ied][1] << endl;
          }
          cout.flush();
          cout << " debug: rank,ii,icv: " << mpi_rank << " " << ii << " " << icv << endl;
          char filename[128];
          sprintf(filename,"final.%02d.%06d.dat",mpi_rank,icv);
          vdReturnDataVec[ii].writeTecplot(filename,x_vd[icv]);
        }

        // make sure no edge references a negative cvofa. These should all have been fused...
        for (int ied = 0; ied < vdReturnDataVec[ii].ned; ++ied) {
          FOR_I2 {
            const int ifa = faoed[ied][i];
            if (ifa < 0) {
              assert((-ifa-1 >= 0)&&(-ifa-1 < vdReturnDataVec[ii].nfa));
              if (vdReturnDataVec[ii].cvofa[-ifa-1] < 0) {
                my_ierr = -1;
                vdReturnDataVec[ii].cvofa[-ifa-1] = -vdReturnDataVec[ii].cvofa[-ifa-1]-1; // need a valid icv_nbr
              }
              //assert(vdReturnDataVec[ii].cvofa[-ifa-1] >= 0);
            }
          }
        }
        if (my_ierr != 0) {
          cout << " > failed to fuse. try --DEBUG_FUSE_RANK " << mpi_rank << " --DEBUG_FUSE_ICV " << icv << endl;
          // cleanup and leave...
          for (int jj = ii; jj < vdReturnDataVec.size(); ++jj)
            vdReturnDataVec[jj].orphan_chunk_data = -1;
          for (int iocd = 0; iocd < ocdVec.size(); ++iocd)
            ocdVec[iocd].clear();
          break;
        }

        vdReturnDataVec[ii].orphan_chunk_data = -1;
      }

      // there shouldn't be any orphans now...
      ocdVec.clear();

      delete[] x_vd_g;
      if (faofa0 != NULL) delete[] faofa0;
      if (no_flag != NULL) delete[] no_flag;
      if (x_no != NULL) delete[] x_no;
      if (nooed != NULL) delete[] nooed;
      if (faoed != NULL) delete[] faoed;
      if (edono_v != NULL) delete[] edono_v;

      int ierr;
      MPI_Allreduce(&my_ierr,&ierr,1,MPI_INT,MPI_MIN,mpi_comm);
      if (ierr != 0) {
        if (mpi_rank == 0)
          cout << " > increasing nsmooth..." << endl;
        return false;
      }

      double my_buf[2] = {my_d2_max,-my_d2_close_min};
      double buf[2];
      MPI_Reduce(my_buf,buf,2,MPI_DOUBLE,MPI_MAX,0,mpi_comm);
      if (mpi_rank == 0)
        cout << " > d2_max: " << buf[0] << " d2_close_min: " << -buf[1] << " d2_tol: " << d2_tol << endl;

      MPI_Sync("fuseOrphanChunksTol end");

      return true;
    }

    void rebuildVd(const int nsmooth_inc) {

      MPI_Sync("rebuildVd");

      // there could of been a previous call to nsmooth
      const int nsmooth0 = nsmooth;
      nsmooth += nsmooth_inc;
      smoothing_step = nsmooth0;

      if (mpi_rank == 0)
        cout << "StitchNS::rebuildVd(): nsmooth: " << nsmooth_inc << endl;

      assert(&ncv == &ncv);
      assert(ncv >= 0); // could be zero I guess
      assert(x_vd != NULL);
      assert(delta_vd != NULL);
      assert(flag_cv != NULL);

      // either first call to rebuildVd or after clearSmoothing...

      if (xp0_and_delta0 == NULL) {
        // used for SMOOTH_LIMIT
        xp0_and_delta0 = new double[ncv][4];
        FOR_ICV {
          FOR_I3 xp0_and_delta0[icv][i] = x_vd[icv][i];
          xp0_and_delta0[icv][3] = delta_vd[icv];
          // set everyone to active and rebuild...
          assert(flag_cv[icv] == (flag_cv[icv]&PART_MASK_BITS));
          flag_cv[icv] |= ACTIVE_BIT|REBUILD_BIT;
        }
      }

      // we also need the vol_cv and x_cv...

      if (vol_cv == NULL) vol_cv = new double[ncv];
      if (x_cv == NULL) x_cv = new double[ncv][3];

      int * send_count = new int[mpi_size];
      int * send_disp = new int[mpi_size];
      int * recv_count = new int[mpi_size];
      int * recv_disp = new int[mpi_size];

      if (smoothing_step == 0) {
        if (mpi_rank == 0) {
          cout << "=====================================================" << endl;
          cout << "first Voronoi build: " << endl;
          cout << "=====================================================" << endl;
        }
	PartData::ensureBbminmax();
	rebuildFlaggedVds(true,smoothing_step==nsmooth,send_count,send_disp,recv_count,recv_disp); // first bool is surface flag
      }

      while (1) {

	MPI_Sync("1");

        ++smoothing_step;
        if (smoothing_step > nsmooth) {
          // if you successfully fuse the vd's you can leave and write, otherwise continue smoothing...
          if (fuseOrphanChunksTol())
            break;
          else
            ++nsmooth;
        }

        if (mpi_rank == 0) {
          cout << "=====================================================" << endl;
          cout << "smoothing step: " << smoothing_step << " of " << nsmooth << endl;
          cout << "=====================================================" << endl;
        }

        uint* smooth_flag = new uint[ncv];
        flagSmoothedVds(smooth_flag,send_count,send_disp,recv_count,recv_disp);

        // move x_vd to x_cv in certain cvs and rebuild them. Here we
        // only set the rebuild bit in cvs that moved or in the cvs they touch. More cvs may have to
        // be rebuilt, and this gets enforced by clearing the NEW_BIT, meaning that certain
        // vds may have picked up new nbrs that were not in the current  set...

        double my_dx_minmax[2] = { HUGE_VAL, HUGE_VAL };
        double my_buf[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
        assert(x_cv);
        const int npart = partVec.size();
        FOR_ICV {
          ++my_buf[0]; // global count
          assert(flag_cv[icv]&ACTIVE_BIT);
          flag_cv[icv] |= ACTIVE0_BIT; // try and make this seem like the moving solver as much as possible
          flag_cv[icv] &= ~NEW_BIT;
          assert(!(flag_cv[icv]&REBUILD_BIT));
          assert(!(flag_cv[icv]&TMP1_BIT)); // used below in nbr flagging
          // skip unsmoothed cells...
          if (smooth_flag[icv] != 0) {
            // clear the NEW_BIT, because we are now
            const double val = 2.0*DIST(x_cv[icv],x_vd[icv])/delta_vd[icv];
            if (val > 1.0E-12) {
              // this cv is not centroidal.

              // Move its x_vd to its current centroid
              //FOR_I3 x_vd[icv][i] = x_cv[icv][i];

              // limit total motion to some fraction of original delta...
              double dx0[3];
              FOR_I3 dx0[i] = x_cv[icv][i] - xp0_and_delta0[icv][i];
              const double mag = MAG(dx0);
              if (mag > smooth_limit*xp0_and_delta0[icv][3]) {
                ++my_buf[4]; // limiting smooth count
                FOR_I3 dx0[i] *= smooth_limit*xp0_and_delta0[icv][3]/mag;
              }

              // if this is a seed point that came from a part, it will have additional constraints
              // (e.g. typically setting the motion to zero). In b_force_smooth_all, we skip any
              // constraints and apply looyd iteration everywhere, possibly limited by SMOOTH_LIMIT...
              if ((!b_force_smooth_all)&&((flag_cv[icv]&PART_MASK_BITS) < npart)) {
                ++my_buf[5]; // is Part constraints
                partVec[flag_cv[icv]&PART_MASK_BITS]->constrainSmoothing(dx0,x_vd[icv]);
              }
              FOR_I3 x_vd[icv][i] += dx0[i] - (x_vd[icv][i] - xp0_and_delta0[icv][i]);

              // and flag as moved...
              flag_cv[icv] |= REBUILD_BIT;
              ++my_buf[1];
              my_dx_minmax[0] = min(my_dx_minmax[0],val);
              my_dx_minmax[1] = min(my_dx_minmax[1],-val);
              my_buf[2] += val;
            }
          }
        }
        delete[] smooth_flag;

        // we also need to flag anyone touching a rebuild cv...
        // use TMP1_BIT to do this...

        vector<int> intVec;
        FOR_RANK send_count[rank] = 0;

        for (int ii = 0; ii < vdReturnDataVec.size(); ++ii) {
          const int icv = vdReturnDataVec[ii].icv;
          assert(ivrdocv[icv] == ii);
          if ((icv >= 0)&&(flag_cv[icv]&REBUILD_BIT)) {
            for (int ifa = 0; ifa < vdReturnDataVec[ii].nfa; ++ifa) {
              const int icv_nbr = vdReturnDataVec[ii].cvofa[ifa];
              assert(icv_nbr != icv);
              if (icv_nbr < ncv) {
                // this nbr is local, so it needs its rebuild bit set. Note that
                // we do NOT just set the local rebuild bit here because that could
                // lead to creep in this loop because it is looking for REBUILD. Use
                // TMP1_BIT and then flip it to rebuild outside this loop...
                if (!(flag_cv[icv_nbr]&REBUILD_BIT))
                  flag_cv[icv_nbr] |= TMP1_BIT;
              }
              else {
                // this is an off-rank nbr that needs rebuilding...
                assert(icv_nbr-ncv < rbi_g.size());
                int rank,bits,index;
                BitUtils::unpackRankBitsIndex(rank,bits,index,rbi_g[icv_nbr-ncv]);
                ++send_count[rank]; // we need to send just send the index - this guy has to be rebuilt
                intVec.push_back(rank);
                intVec.push_back(index);
              }
            }
          }
        }

        send_disp[0] = 0;
        for (int rank = 1; rank < mpi_size; ++rank)
          send_disp[rank] = send_disp[rank-1] + send_count[rank-1];
        const int send_count_sum = send_disp[mpi_size-1] + send_count[mpi_size-1];
        assert(intVec.size() == send_count_sum*2);

        int * send_buf_int = new int[send_count_sum];

        // use the intVec to pack the rebuild cvs...

        for (int ii = 0; ii < send_count_sum; ++ii) {
          const int rank = intVec[2*ii  ];
          const int icv  = intVec[2*ii+1];
          send_buf_int[send_disp[rank]++] = icv;
        }
        intVec.clear();

        // rewind...

        send_disp[0] = 0;
        for (int rank = 1; rank < mpi_size; ++rank)
          send_disp[rank] = send_disp[rank-1] + send_count[rank-1];

        MPI_Alltoall(send_count,1,MPI_INT,recv_count,1,MPI_INT,mpi_comm);

        recv_disp[0] = 0;
        for (int rank = 1; rank < mpi_size; ++rank)
          recv_disp[rank] = recv_count[rank-1] + recv_disp[rank-1];
        const int recv_count_sum = recv_disp[mpi_size-1] + recv_count[mpi_size-1];

        int * recv_buf_int = new int[recv_count_sum];

        MPI_Alltoallv(send_buf_int,send_count,send_disp,MPI_INT,
            recv_buf_int,recv_count,recv_disp,MPI_INT,
            mpi_comm);
        delete[] send_buf_int;

        for (int irecv = 0; irecv < recv_count_sum; ++irecv) {
          const int icv = recv_buf_int[irecv];
          assert((icv >= 0)&&(icv < ncv));
          if (!(flag_cv[icv]&REBUILD_BIT))
            flag_cv[icv] |= TMP1_BIT;
        }

        delete[] recv_buf_int; recv_buf_int = NULL;

        // finally flip the TMP1_BITS to REBUILD_BITS, and count these
        // as well...
        FOR_ICV {
          if (flag_cv[icv]&TMP1_BIT) {
            assert(!(flag_cv[icv]&REBUILD_BIT));
            flag_cv[icv] &= ~TMP1_BIT;
            flag_cv[icv] |= REBUILD_BIT;
            ++my_buf[3]; // count of cvs where the x_vd did not change but the current nbr did, so need to rebuilt too
          }
        }

        {
          double dx_minmax[2];
          MPI_Reduce(my_dx_minmax,dx_minmax,2,MPI_DOUBLE,MPI_MIN,0,mpi_comm);
          double buf[6];
          MPI_Reduce(my_buf,buf,6,MPI_DOUBLE,MPI_SUM,0,mpi_comm);
          if (mpi_rank == 0) {
            cout << " > smoothing " << buf[1] << " out of " << buf[0] << " cvs. normalized dx [min,avg,max]: " <<
              dx_minmax[0] << " " << buf[2]/buf[1] << " " << -dx_minmax[1] << endl;
            cout << " > rebuild flag set in " << buf[1] << " cvs and " << buf[3] << " nbrs: " << buf[1]+buf[3] << " total" << endl;
            cout << " > limiting smooth in " << buf[4] << " cvs" << endl;
            cout << " > applying part constraints in " << buf[5] << " cvs" << endl;
          }
        }

        // check if any parts moved, and clear their point lookup adts appropriately. For
        // now, just clear the main point adts...

        assert(hcp_x_vd_bbox_adt); delete hcp_x_vd_bbox_adt; hcp_x_vd_bbox_adt = NULL;
        assert(hcp_x_vd_adt); delete hcp_x_vd_adt; hcp_x_vd_adt = NULL;
        if (b_force_smooth_all) {
          for (int ipart = 0; ipart < npart; ++ipart) {
            if (partVec[ipart]->pts) {
              assert(partVec[ipart]->x_vd_bbox_adt);
              delete partVec[ipart]->x_vd_bbox_adt;
              partVec[ipart]->x_vd_bbox_adt = NULL;
              assert(partVec[ipart]->x_vd_adt);
              delete partVec[ipart]->x_vd_adt;
              partVec[ipart]->x_vd_adt = NULL;
            }
          }
        }

        rebuildFlaggedVds(true,smoothing_step==nsmooth,send_count,send_disp,recv_count,recv_disp); // surface flag
      }

      MPI_Sync("2");

      // flag inactive faces...

      double my_area2_max = 0.0;
      const double area2_ratio_zero = 1.0E-24; // copied over for VoronoiBuilder
      int nfa_max = 0;
      double (*n_fa)[3] = NULL;
      for (int ii = 0; ii < ncv; ++ii) {
        assert(vdReturnDataVec[ii].nfa > 0);
        DELETE(vdReturnDataVec[ii].bits_fa);
        vdReturnDataVec[ii].bits_fa = new int[vdReturnDataVec[ii].nfa]; // hold first node of face for now
        if (vdReturnDataVec[ii].nfa > nfa_max) {
          nfa_max = vdReturnDataVec[ii].nfa;
          DELETE(n_fa);
          n_fa = new double[nfa_max][3];
        }
        DELETE(vdReturnDataVec[ii].r2_fa);
        vdReturnDataVec[ii].r2_fa = new double[vdReturnDataVec[ii].nfa];
        for (int ifa = 0; ifa < vdReturnDataVec[ii].nfa; ++ifa) {
          FOR_J3 n_fa[ifa][j] = 0.0;
          vdReturnDataVec[ii].bits_fa[ifa] = -1;
          vdReturnDataVec[ii].r2_fa[ifa] = 0.0; // overwrite with r2_max from VoronoiBuilder
        }
        for (int ied = 0; ied < vdReturnDataVec[ii].ned; ++ied) {
          FOR_I2 {
            const int ifa = -vdReturnDataVec[ii].faoed[ied][i]-1;
            if (ifa >= 0) {
              assert(ifa < vdReturnDataVec[ii].nfa);
              vdReturnDataVec[ii].r2_fa[ifa] = max(vdReturnDataVec[ii].r2_fa[ifa],
                  max(DOT_PRODUCT(vdReturnDataVec[ii].x_no[vdReturnDataVec[ii].nooed[ied][0]],
                                  vdReturnDataVec[ii].x_no[vdReturnDataVec[ii].nooed[ied][0]]),
                      DOT_PRODUCT(vdReturnDataVec[ii].x_no[vdReturnDataVec[ii].nooed[ied][1]],
                                  vdReturnDataVec[ii].x_no[vdReturnDataVec[ii].nooed[ied][1]])));
              if (vdReturnDataVec[ii].bits_fa[ifa] == -1) {
                vdReturnDataVec[ii].bits_fa[ifa] = vdReturnDataVec[ii].nooed[ied][i];
              }
              else {
                assert((vdReturnDataVec[ii].bits_fa[ifa] >= 0)&&(vdReturnDataVec[ii].bits_fa[ifa] < vdReturnDataVec[ii].nno));
                const double this_n[3] = TRI_NORMAL_2(vdReturnDataVec[ii].x_no[vdReturnDataVec[ii].bits_fa[ifa]],
                                                      vdReturnDataVec[ii].x_no[vdReturnDataVec[ii].nooed[ied][i]],
                                                      vdReturnDataVec[ii].x_no[vdReturnDataVec[ii].nooed[ied][1-i]]);
                FOR_J3 n_fa[ifa][j] += this_n[j];
              }
            }
          }
        }
        DELETE(vdReturnDataVec[ii].area2_fa);
        vdReturnDataVec[ii].area2_fa = new double[vdReturnDataVec[ii].nfa];
        for (int ifa = 0; ifa < vdReturnDataVec[ii].nfa; ++ifa) {
          vdReturnDataVec[ii].area2_fa[ifa] = DOT_PRODUCT(n_fa[ifa],n_fa[ifa]);
          if (vdReturnDataVec[ii].area2_fa[ifa] > vdReturnDataVec[ii].r2_fa[ifa]*vdReturnDataVec[ii].r2_fa[ifa]*area2_ratio_zero) {
            vdReturnDataVec[ii].bits_fa[ifa] = VALID_LOCAL_FACE_BIT;
          }
          else {
            vdReturnDataVec[ii].bits_fa[ifa] = ZERO_LOCAL_FACE_BIT;
            my_area2_max = max(my_area2_max,vdReturnDataVec[ii].area2_fa[ifa]/(vdReturnDataVec[ii].r2_fa[ifa]*vdReturnDataVec[ii].r2_fa[ifa]));
          }
        }
      }
      DELETE(n_fa);

      FOR_RANK send_count[rank] = 0;
      for (int ii = 0; ii < ncv; ++ii) {
        const int icv = vdReturnDataVec[ii].icv;
        assert((icv >= 0)&&(icv < ncv));
        for (int ifa = 0; ifa < vdReturnDataVec[ii].nfa; ++ifa) {
          const int icv_nbr = vdReturnDataVec[ii].cvofa[ifa];
          // try to find duplicated face on nbr...
          if (icv_nbr < ncv) {
            const int ii_nbr = ivrdocv[icv_nbr];
            assert((ii_nbr >= 0)&&(ii_nbr < vdReturnDataVec.size()));
            assert(vdReturnDataVec[ii_nbr].icv == icv_nbr);
            int ifa_nbr;
            for (ifa_nbr = 0; ifa_nbr < vdReturnDataVec[ii_nbr].nfa; ++ifa_nbr) {
              const int icv_nbr_nbr = vdReturnDataVec[ii_nbr].cvofa[ifa_nbr];
              if (icv_nbr_nbr == icv)
                break;
            }
            if (ifa_nbr != vdReturnDataVec[ii_nbr].nfa) {
              // we found our nbr so we can flag as valid if the face is large enough...
              if (vdReturnDataVec[ii_nbr].bits_fa[ifa_nbr]&VALID_LOCAL_FACE_BIT) {
                vdReturnDataVec[ii].bits_fa[ifa] |= VALID_NBR_FACE_BIT;
              }
              else {
                assert(vdReturnDataVec[ii_nbr].bits_fa[ifa_nbr]&ZERO_LOCAL_FACE_BIT);
                vdReturnDataVec[ii].bits_fa[ifa] |= ZERO_NBR_FACE_BIT;
              }
            }
            else {
              vdReturnDataVec[ii].bits_fa[ifa] |= NO_NBR_FACE_BIT;
            }
            if ((vdReturnDataVec[ii].bits_fa[ifa]&VALID_NBR_FACE_BIT) == 0)
              my_area2_max = max(my_area2_max,vdReturnDataVec[ii].area2_fa[ifa]/(vdReturnDataVec[ii].r2_fa[ifa]*vdReturnDataVec[ii].r2_fa[ifa]));
          }
          else {
            int rank,bits,index;
            BitUtils::unpackRankBitsIndex(rank,bits,index,rbi_g[icv_nbr-ncv]);
            send_count[rank] += 3; // index there, mpi_rank/bits_fa, icv
          }
        }
      }

      send_disp[0] = 0;
      for (int rank = 1; rank < mpi_size; ++rank)
        send_disp[rank] = send_disp[rank-1] + send_count[rank-1];
      const int send_count_sum = send_disp[mpi_size-1] + send_count[mpi_size-1];

      int * send_buf_int = new int[send_count_sum];
      for (int ii = 0; ii < ncv; ++ii) {
        const int icv = vdReturnDataVec[ii].icv;
        assert((icv >= 0)&&(icv < ncv));
        for (int ifa = 0; ifa < vdReturnDataVec[ii].nfa; ++ifa) {
          const int icv_nbr = vdReturnDataVec[ii].cvofa[ifa];
          // try to find duplicated face on nbr...
          if (icv_nbr >= ncv) {
            int rank,bits,index;
            BitUtils::unpackRankBitsIndex(rank,bits,index,rbi_g[icv_nbr-ncv]);
            send_buf_int[send_disp[rank]++] = index;
            send_buf_int[send_disp[rank]++] = BitUtils::packRankBits(mpi_rank,BitUtils::flipPeriodicBits(bits));
            send_buf_int[send_disp[rank]++] = icv;
          }
        }
      }

      // rewind...

      send_disp[0] = 0;
      for (int rank = 1; rank < mpi_size; ++rank)
        send_disp[rank] = send_disp[rank-1] + send_count[rank-1];

      MPI_Alltoall(send_count,1,MPI_INT,recv_count,1,MPI_INT,mpi_comm);

      recv_disp[0] = 0;
      for (int rank = 1; rank < mpi_size; ++rank)
        recv_disp[rank] = recv_count[rank-1] + recv_disp[rank-1];
      int recv_count_sum = recv_disp[mpi_size-1] + recv_count[mpi_size-1];

      int * recv_buf_int = new int[recv_count_sum];
      MPI_Alltoallv(send_buf_int,send_count,send_disp,MPI_INT,
          recv_buf_int,recv_count,recv_disp,MPI_INT,mpi_comm);

      assert(recv_count_sum%3 == 0);
      recv_count_sum /= 3;
      for (int irecv = 0; irecv < recv_count_sum; ++irecv) {
        const int icv = recv_buf_int[irecv*3]; assert((icv >= 0)&&(icv < ncv));
        int rank_nbr,flip_bits_nbr; BitUtils::unpackRankBits(rank_nbr,flip_bits_nbr,recv_buf_int[irecv*3+1]);
        const int icv_nbr = recv_buf_int[irecv*3+2];
        const int ii = ivrdocv[icv];
        assert((ii >= 0)&&(ii < vdReturnDataVec.size()));
        assert(vdReturnDataVec[ii].icv == icv);
        int ifa;
        for (ifa = 0; ifa < vdReturnDataVec[ii].nfa; ++ifa) {
          if (vdReturnDataVec[ii].cvofa[ifa] >= ncv) {
            int rank,bits,index;
            BitUtils::unpackRankBitsIndex(rank,bits,index,rbi_g[vdReturnDataVec[ii].cvofa[ifa]-ncv]);
            if ((rank_nbr == rank)&&(flip_bits_nbr == bits)&&(icv_nbr == index))
              break;
          }
        }
        if (ifa != vdReturnDataVec[ii].nfa) {
          // we found our nbr so we can flag as valid if large enough...
          if (vdReturnDataVec[ii].bits_fa[ifa]&VALID_LOCAL_FACE_BIT) {
            recv_buf_int[irecv] = VALID_NBR_FACE_BIT;
          }
          else {
            assert(vdReturnDataVec[ii].bits_fa[ifa]&ZERO_LOCAL_FACE_BIT);
            recv_buf_int[irecv] = ZERO_NBR_FACE_BIT;
          }
        }
        else {
          recv_buf_int[irecv] = NO_NBR_FACE_BIT;
        }
      }

      // reverse send of face bits...

      FOR_RANK {
        send_count[rank] /= 3;
        send_disp[rank] /= 3;
        recv_count[rank] /= 3;
        recv_disp[rank] /= 3;
      }

      MPI_Alltoallv(recv_buf_int,recv_count,recv_disp,MPI_INT,
          send_buf_int,send_count,send_disp,MPI_INT,mpi_comm);
      delete[] recv_buf_int;

      for (int ii = 0; ii < ncv; ++ii) {
        const int icv = vdReturnDataVec[ii].icv;
        assert((icv >= 0)&&(icv < ncv));
        for (int ifa = 0; ifa < vdReturnDataVec[ii].nfa; ++ifa) {
          const int icv_nbr = vdReturnDataVec[ii].cvofa[ifa];
          // try to find duplicated face on nbr...
          if (icv_nbr >= ncv) {
            int rank,bits,index;
            BitUtils::unpackRankBitsIndex(rank,bits,index,rbi_g[icv_nbr-ncv]);
            vdReturnDataVec[ii].bits_fa[ifa] |= send_buf_int[send_disp[rank]++];
            if ((vdReturnDataVec[ii].bits_fa[ifa]&VALID_NBR_FACE_BIT) == 0)
              my_area2_max = max(my_area2_max,vdReturnDataVec[ii].area2_fa[ifa]/(vdReturnDataVec[ii].r2_fa[ifa]*vdReturnDataVec[ii].r2_fa[ifa]));
          }
        }
      }
      delete[] send_buf_int;

      double area2_max;
      MPI_Reduce(&my_area2_max,&area2_max,1,MPI_DOUBLE,MPI_MAX,0,mpi_comm);
      if (mpi_rank == 0)
        cout << " > max relative area2 for mismatched/ignored faces (should be small): " << area2_max << endl;

      // cleanup...
      delete[] send_count;
      delete[] send_disp;
      delete[] recv_count;
      delete[] recv_disp;

    }

    double reconcileNodeLoops(vector<pair<double*,double*> >& loopVec,
        vector<double*>& x0Vec,const int i0_f,const int i0_l,
        vector<double*>& x1Vec,const int i1_f,const int i1_l,
        const double dx[3],const double delta,const double d2_tol) {

      static int counter = 0;
      ++counter;

      const int debug_rank = getIntParam("DEBUG_RANK",-1);
      const int debug_count = getIntParam("DEBUG_COUNT",-1);

      // find the closest nodes (not already flagged)...

      assert(x0Vec[i0_l] == NULL);
      assert(x1Vec[i1_l] == NULL);

      // step 0: how many loops:

      int nloop = 0; // should be same for both loops, so just name nloop
      for (int i0 = i0_f; i0 <= i0_l; ++i0) {
        if (x0Vec[i0] == NULL) {
          ++nloop;
        }
      }
      int nloop_check = 0;
      for (int i1 = i1_f; i1 <= i1_l; ++i1) {
        if (x1Vec[i1] == NULL) {
          ++nloop_check;
        }
      }
      if (nloop != nloop_check) {
        cout << "x0Vec" << endl;
        for (int i0 = i0_f; i0 <= i0_l; ++i0) {
          if (x0Vec[i0] != NULL)
            cout << "[ " << x0Vec[i0][0]-dx[0] << "," << x0Vec[i0][1]-dx[1] << "," << x0Vec[i0][2]-dx[2] << " ] " << endl;
          else
            cout << "NULL" << endl;
        }
        cout << "x1Vec" << endl;
        for (int i1 = i1_f; i1 <= i1_l; ++i1) {
          if (x1Vec[i1] != NULL)
            cout << COUT_VEC(x1Vec[i1]) << endl;
          else
            cout << "NULL" << endl;
        }
        cout.flush();
      }
      assert(nloop == nloop_check); // if you hit this, good luck
      int * loop1_flag = new int[nloop];
      int * loop0_f = new int[nloop];
      int * loop1_f = new int[nloop];
      int * loop0_l = new int[nloop];
      int * loop1_l = new int[nloop];
      for (int iloop = 0; iloop < nloop; ++iloop) {
        loop1_flag[iloop] = 0; // set to 1 when this loop1 has been used
      }
      {
        int iloop = 0;
        loop0_f[0] = i0_f;
        for (int i0 = i0_f; i0 <= i0_l; ++i0) {
          if (x0Vec[i0] == NULL) {
            loop0_l[iloop] = i0-1;
            ++iloop;
            if (iloop < nloop) loop0_f[iloop] = i0+1;
          }
        }
        assert(iloop == nloop);
        iloop = 0;
        loop1_f[0] = i1_f;
        for (int i1 = i1_f; i1 <= i1_l; ++i1) {
          if (x1Vec[i1] == NULL) {
            loop1_l[iloop] = i1-1;
            ++iloop;
            if (iloop < nloop) loop1_f[iloop] = i1+1;
          }
        }
        assert(iloop == nloop);
      }

      /*
         cout << "i0_f: " << i0_f << " i0_l: " << i0_l << " loop0_f[0]: " << loop0_f[0] << " loop0_l[0]: " << loop0_l[0] << endl;
         cout << "i1_f: " << i1_f << " i1_l: " << i1_l << " loop1_f[0]: " << loop1_f[0] << " loop1_l[0]: " << loop1_l[0] << endl;
         getchar();
         */

      for (int iloop0 = 0; iloop0 < nloop; ++iloop0) {
        int i0_closest = -1,i1_closest=0,iloop1_closest=0; // 0 to avoid compiler warnings
        double d2_min = HUGE_VAL;
        for (int i0 = loop0_f[iloop0]; i0 <= loop0_l[iloop0]; ++i0) {
          for (int iloop1 = 0; iloop1 < nloop; ++iloop1) {
            // only consider loop1's we have not used yet...
            if (loop1_flag[iloop1] == 0) {
              for (int i1 = loop1_f[iloop1]; i1 <= loop1_l[iloop1]; ++i1) {
                double d2 = 0;
                FOR_I3 {
                  const double dxi = x1Vec[i1][i] - x0Vec[i0][i] + dx[i];
                  d2 += dxi*dxi;
                }
                if ((i0_closest == -1)||(d2 < d2_min)) {
                  i0_closest = i0;
                  i1_closest = i1;
                  iloop1_closest = iloop1;
                  d2_min = d2;
                }
              }
            }
          }
        }
        assert(i0_closest != -1);
        loop1_flag[iloop1_closest] = 1; // don't use this loop again

        // if this is not the first loop, we need to add the very first pair
        // to the vecLoop again...

        //if (iloop0 > 0) loopVec.push_back(loopVec.front()); // remove when using NULL pairs to segment loops
        //const int loopVec_size0 = loopVec.size();

        int i0 = i0_closest;
        int i1 = i1_closest;
        assert(x0Vec[i0]);
        assert(x1Vec[i1]);

        double dist0 = 0.0; // sum arc-lengths
        double dist1 = 0.0;
        bool i0_looped = false;
        bool i1_looped = false;
        while ((!i0_looped)||(!i1_looped)) {

          // if coincident, add to loopVec
          double d2 = 0;
          FOR_I3 {
            const double dxi = x1Vec[i1][i] - x0Vec[i0][i] + dx[i];
            d2 += dxi*dxi;
          }

          if ((mpi_rank == debug_rank)&&(counter == debug_count)) cout << "XXXXXXXXXXXXXX  i0-i0_f: " << i0-i0_f << " i1-i1_f: " << i1-i1_f << " d2/(delta*delta): " << d2/(delta*delta) << endl;

          if (d2 < d2_tol*delta*delta) {
            loopVec.push_back(pair<double*,double*>(x0Vec[i0],x1Vec[i1]));
            //cout << "got coincident point, i0: " << i0 << " i1: " << i1 << " dist: " << sqrt(d2) << " dist2/(delta*delta): " << d2/(delta*delta) << endl;
            // and reset the cumulative distances...
            dist0 = dist1 = 0.0;
          }

          // increment
          int i0_p = i0;
          if (i0 == loop0_l[iloop0]) i0 = loop0_f[iloop0];
          else ++i0;

          int i1_p = i1;
          if (i1 == loop1_f[iloop1_closest]) i1 = loop1_l[iloop1_closest];
          else --i1;

          assert(x0Vec[i0]);
          assert(x1Vec[i1]);

          if (i0_looped) {
            i0 = i0_p;
            i1_looped = (i1 == i1_closest);
          }
          else if (i1_looped) {
            i1 = i1_p;
            i0_looped = (i0 == i0_closest);
          }
          else {
            // compare lengths...
            const double ddist0 = DIST(x0Vec[i0],x0Vec[i0_p]);
            const double ddist1 = DIST(x1Vec[i1],x1Vec[i1_p]);
            if ( (dist0 + ddist0) > (dist1 + ddist1) )  {
              i0 = i0_p; // keep i1, send i0 back
              dist1 += ddist1;
              i1_looped = (i1 == i1_closest);
            }
            else {
              i1 = i1_p; // keep i0, send i1 back
              dist0 += ddist0;
              i0_looped = (i0 == i0_closest);
            }
          }
          //cout << "iter: " << iter << " i0 " << i0 << " " << i1 << endl;

        }

        //assert(loopVec.size()-loopVec_size0 >= 3); // atleast need a tri -- NOPE. segments can occur and are ok at this point
        //if (iloop0 > 0) loopVec.push_back(loopVec[loopVec_size0]); // remove when using NULL pair to segment loop
        loopVec.push_back(pair<double*,double*>(NULL,NULL));

      } // iloop


      for (int iloop1 = 0; iloop1 < nloop; ++iloop1) {
        assert(loop1_flag[iloop1] == 1);
      }

      delete[] loop1_l;
      delete[] loop0_l;
      delete[] loop1_f;
      delete[] loop0_f;
      delete[] loop1_flag;

      // for geometry, check normals...

      double n0[3] = { 0.0, 0.0, 0.0 };
      assert(x0Vec[i0_f]);
      int i00 = i0_f;
      for (int i0 = i0_f+1; i0 < i0_l; ++i0) {
        if (x0Vec[i0] == NULL) {
          i00 = i0+1;
        }
        else if (x0Vec[i0+1]) {
          const double this_n2[3] = TRI_NORMAL_2(x0Vec[i00],x0Vec[i0],x0Vec[i0+1]);
          FOR_I3 n0[i] += this_n2[i];
        }
      }

      double n1[3] = { 0.0, 0.0, 0.0 };
      assert(x1Vec[i1_f]);
      int i10 = i1_f;
      for (int i1 = i1_f+1; i1 < i1_l; ++i1) {
        if (x1Vec[i1] == NULL) {
          i10 = i1+1;
        }
        else if (x1Vec[i1+1]) {
          const double this_n2[3] = TRI_NORMAL_2(x1Vec[i10],x1Vec[i1],x1Vec[i1+1]);
          FOR_I3 n1[i] += this_n2[i];
        }
      }

      double n_loop[3] = { 0.0, 0.0, 0.0 };
      double x0[3]; FOR_I3 x0[i] = 0.5*(loopVec[0].first[i]+loopVec[0].second[i]);
      for (int iloop = 1; iloop < loopVec.size()-1; ++iloop) {
        if (loopVec[iloop].first == NULL) {
          assert(loopVec[iloop].second == NULL);
          FOR_I3 x0[i] = 0.5*(loopVec[iloop+1].first[i]+loopVec[iloop+1].second[i]);
        }
        else if (loopVec[iloop+1].first) {
          assert(loopVec[iloop+1].second);
          double x1[3]; FOR_I3 x1[i] = 0.5*(loopVec[iloop].first[i]+loopVec[iloop].second[i]);
          double x2[3]; FOR_I3 x2[i] = 0.5*(loopVec[iloop+1].first[i]+loopVec[iloop+1].second[i]);
          const double this_n2[3] = TRI_NORMAL_2(x0,x1,x2);
          FOR_I3 n_loop[i] += this_n2[i];
        }
      }

      //double dn[3] = { n0[0]+n1[0], n0[1]+n1[1], n0[2]+n1[2] };
      double dn_loop[3] = {
        n_loop[0]-0.5*(n0[0]-n1[0]),
        n_loop[1]-0.5*(n0[1]-n1[1]),
        n_loop[2]-0.5*(n0[2]-n1[2])
      };

      const double dn_loop_rel_err_2 = DOT_PRODUCT(dn_loop,dn_loop)/(delta*delta*delta*delta);

      /*
         if (dn_loop_rel_err_2 > 1.0E-24) {
         cout << "PROBLEM: mpi_rank: " << mpi_rank << " counter: " << counter << " " << dn_loop_rel_err_2 << endl;


         FILE * fp = fopen("x0Vec.dat","w");
         for (int i0 = i0_f; i0 <= i0_l; ++i0) {
         if (x0Vec[i0]) {
         fprintf(fp,"%18.15le %18.15le %18.15le\n",
         x0Vec[i0][0],
         x0Vec[i0][1],
         x0Vec[i0][2]);
         }
         }
         fclose(fp);

         fp = fopen("x1Vec.dat","w");
         for (int i1 = i1_f; i1 <= i1_l; ++i1) {
         if (x1Vec[i1]) {
         fprintf(fp,"%18.15le %18.15le %18.15le\n",
         x0Vec[i1][0],
         x0Vec[i1][1],
         x0Vec[i1][2]);
         }
         }
         fclose(fp);

         fp = fopen("xpVec.dat","w");
         for (int iloop = 0; iloop < loopVec.size(); ++iloop) {
         double x[3]; FOR_I3 x[i] = 0.5*(loopVec[iloop].first[i]+loopVec[iloop].second[i]);
         fprintf(fp,"%18.15le %18.15le %18.15le\n",x[0],x[1],x[2]);
         }
         fclose(fp);

         assert(0);
         }
         */

      //return pair<int,double>(nloop,dn_loop_rel_err_2);
      return dn_loop_rel_err_2;

      /*
         cout << "compare: difference between loops: " << DOT_PRODUCT(dn,dn)/(delta*delta*delta*delta) << " final loop - mean loop: " << DOT_PRODUCT(dnloop,dnloop)/(delta*delta*delta*delta) << endl;
         getchar();
         */

    }

    class FaData {
      public:
        uint8 rbi_nbr;
        double n[3];
        double x[3];
        FaData(const uint8 rbi_nbr,const FaceGeometryData& fgd) {
          this->rbi_nbr = rbi_nbr;
          FOR_I3 n[i] = fgd.normal[i];
          FOR_I3 x[i] = fgd.xc[i];
        }
    };

    class BfData {
      public:
        int zone;
        double n[3];
        double x[3];
        double area;
        double area_over_delta;
        double Gij[3][3];
        BfData(const int zone,const BfGeometryData& fgd) {
          this->zone = zone;
          FOR_I3 n[i] = fgd.normal[i];
          FOR_I3 x[i] = fgd.xc[i];
          area = fgd.area;
          area_over_delta = fgd.area_over_delta;
          FOR_I3 FOR_J3 Gij[i][j] = fgd.Gij[i][j];
        }
    };

    void averageFaDupData(vector<FaData>& fa_dup_data_vec,int * faocv_dup_i) {

      double dn2_max = 0.0;
      double dx2_max = 0.0;

      int* send_count = new int[mpi_size];
      FOR_RANK send_count[rank] = 0;

      for (int icv = 0; icv < ncv; ++icv) {
        for (int ifa_dup = faocv_dup_i[icv]; ifa_dup != faocv_dup_i[icv+1]; ++ifa_dup) {
          int rank,bits,icv_nbr;
          BitUtils::unpackRankBitsIndex(rank,bits,icv_nbr,fa_dup_data_vec[ifa_dup].rbi_nbr);
          if ((rank != mpi_rank)||(bits)) {
            ++send_count[rank]; // 6 doubles
          }
        }
      }

      int* send_disp = new int[mpi_size];
      send_disp[0] = 0;
      for (int rank = 1; rank < mpi_size; ++rank)
        send_disp[rank] = send_count[rank-1] + send_disp[rank-1];
      int send_count_sum = send_disp[mpi_size-1] + send_count[mpi_size-1];

      int * send_int_buf = new int[send_count_sum*4];
      double * send_double_buf = new double[send_count_sum*6];

      for (int icv = 0; icv < ncv; ++icv) {
        for (int ifa_dup = faocv_dup_i[icv]; ifa_dup != faocv_dup_i[icv+1]; ++ifa_dup) {
          const uint8 rbi_nbr = fa_dup_data_vec[ifa_dup].rbi_nbr;
          int rank,bits,icv_nbr;
          BitUtils::unpackRankBitsIndex(rank,bits,icv_nbr,rbi_nbr);
          if ((rank == mpi_rank)&&(bits == 0)) {
            // consider only when icv < icv_nbr...
            if (icv < icv_nbr) {
              const uint8 rbi = BitUtils::packRankBitsIndex(mpi_rank,0,icv);
              int ifa_dup_nbr;
              for (ifa_dup_nbr = faocv_dup_i[icv_nbr]; ifa_dup_nbr != faocv_dup_i[icv_nbr+1]; ++ifa_dup_nbr) {
                if (fa_dup_data_vec[ifa_dup_nbr].rbi_nbr == rbi)
                  break;
              }
              if (ifa_dup_nbr == faocv_dup_i[icv_nbr+1]) {
                cout << "ifa_dup_nbr != faocv_dup_i[icv_nbr+1] for rank0,cv0: "
                  << rank << " " << icv_nbr << ", rank1,cv1: " << mpi_rank << " " << icv << " "
                  << "area: " << 0.5*MAG(fa_dup_data_vec[ifa_dup].n) << endl;
                cout.flush();
              }
              assert(ifa_dup_nbr != faocv_dup_i[icv_nbr+1]);
              double dn[3]; FOR_I3 dn[i] = fa_dup_data_vec[ifa_dup].n[i] + fa_dup_data_vec[ifa_dup_nbr].n[i];
              dn2_max = max(dn2_max,DOT_PRODUCT(dn,dn));
              double x_fa[3];     FOR_I3 x_fa[i] = fa_dup_data_vec[ifa_dup].x[i] + x_vd[icv][i];
              double x_fa_nbr[3]; FOR_I3 x_fa_nbr[i] = fa_dup_data_vec[ifa_dup_nbr].x[i] + x_vd[icv_nbr][i];
              const double dx2 = DIST2(x_fa,x_fa_nbr);
              dx2_max = max(dx2_max,dx2);
              // complete internal geometry...
              FOR_I3 fa_dup_data_vec[ifa_dup].n[i] = 0.5*fa_dup_data_vec[ifa_dup].n[i] - 0.5*fa_dup_data_vec[ifa_dup_nbr].n[i];
              FOR_I3 fa_dup_data_vec[ifa_dup_nbr].n[i] = -fa_dup_data_vec[ifa_dup].n[i];
              FOR_I3 fa_dup_data_vec[ifa_dup].x[i] = 0.5*(x_fa[i] + x_fa_nbr[i]);
              FOR_I3 fa_dup_data_vec[ifa_dup_nbr].x[i] = fa_dup_data_vec[ifa_dup].x[i];
            }
          }
          else {
            const int inv_bits = BitUtils::flipPeriodicBits(bits);
            send_int_buf[send_disp[rank]*4  ] = icv_nbr; // icv on rank
            send_int_buf[send_disp[rank]*4+1] = mpi_rank;
            send_int_buf[send_disp[rank]*4+2] = inv_bits;
            send_int_buf[send_disp[rank]*4+3] = icv;
            // and the doubles...
            FOR_I3 send_double_buf[send_disp[rank]*6+i]   = fa_dup_data_vec[ifa_dup].n[i];
            FOR_I3 send_double_buf[send_disp[rank]*6+3+i] = fa_dup_data_vec[ifa_dup].x[i] + x_vd[icv][i];
            if (bits) {
              PeriodicData::periodicRotate(send_double_buf+send_disp[rank]*6,1,inv_bits);
              PeriodicData::periodicTranslate(send_double_buf+send_disp[rank]*6+3,1,inv_bits);
            }
            send_disp[rank]++;
          }
        }
      }

      // rewind and exchange....

      send_disp[0] = 0;
      for (int rank = 1; rank < mpi_size; ++rank)
        send_disp[rank] = send_count[rank-1] + send_disp[rank-1];

      int* recv_count = new int[mpi_size];
      MPI_Alltoall(send_count,1,MPI_INT,recv_count,1,MPI_INT,mpi_comm);

      int* recv_disp = new int[mpi_size];
      recv_disp[0] = 0;
      for (int rank = 1; rank < mpi_size; ++rank)
        recv_disp[rank] = recv_count[rank-1] + recv_disp[rank-1];
      int recv_count_sum = recv_disp[mpi_size-1] + recv_count[mpi_size-1];

      FOR_RANK {
        send_count[rank] *= 4;
        send_disp[rank] *= 4;
        recv_count[rank] *= 4;
        recv_disp[rank] *= 4;
      }

      int * recv_int_buf = new int[recv_count_sum*4];
      MPI_Alltoallv(send_int_buf,send_count,send_disp,MPI_INT,
          recv_int_buf,recv_count,recv_disp,MPI_INT,mpi_comm);
      delete[] send_int_buf; send_int_buf = NULL;

      FOR_RANK {
        send_count[rank] = (send_count[rank]/4)*6;
        send_disp[rank] = (send_disp[rank]/4)*6;
        recv_count[rank] = (recv_count[rank]/4)*6;
        recv_disp[rank] = (recv_disp[rank]/4)*6;
      }

      double * recv_double_buf = new double[recv_count_sum*6];
      MPI_Alltoallv(send_double_buf,send_count,send_disp,MPI_DOUBLE,
          recv_double_buf,recv_count,recv_disp,MPI_DOUBLE,mpi_comm);
      delete[] send_double_buf;
      delete[] send_count;
      delete[] send_disp;
      delete[] recv_count;
      delete[] recv_disp;

      for (int irecv = 0; irecv < recv_count_sum; irecv++) {
        const int icv = recv_int_buf[irecv*4];
        assert((icv >= 0)&&(icv < ncv));
        int rank_nbr = recv_int_buf[irecv*4+1];
        int bits = recv_int_buf[irecv*4+2];
        int icv_on_nbr = recv_int_buf[irecv*4+3];
        const uint8 rbi_nbr = BitUtils::packRankBitsIndex(rank_nbr,bits,icv_on_nbr);
        double n_nbr[3]; FOR_I3 n_nbr[i] = recv_double_buf[irecv*6+i];
        double x_fa_nbr[3]; FOR_I3 x_fa_nbr[i] = recv_double_buf[irecv*6+3+i];
        int ifa_dup;
        for (ifa_dup = faocv_dup_i[icv]; ifa_dup != faocv_dup_i[icv+1]; ++ifa_dup) {
          if (fa_dup_data_vec[ifa_dup].rbi_nbr == rbi_nbr)
            break;
        }
        assert(ifa_dup != faocv_dup_i[icv+1]);
        double dn[3]; FOR_I3 dn[i] = fa_dup_data_vec[ifa_dup].n[i] + n_nbr[i];
        dn2_max = max(dn2_max,DOT_PRODUCT(dn,dn));
        double x_fa[3];     FOR_I3 x_fa[i] = fa_dup_data_vec[ifa_dup].x[i] + x_vd[icv][i];
        const double dx2 = DIST2(x_fa,x_fa_nbr);
        dx2_max = max(dx2_max,dx2);
        // complete internal geometry...
        FOR_I3 fa_dup_data_vec[ifa_dup].n[i] = 0.5*fa_dup_data_vec[ifa_dup].n[i] - 0.5*n_nbr[i];
        FOR_I3 fa_dup_data_vec[ifa_dup].x[i] = 0.5*(x_fa[i] + x_fa_nbr[i]);
      }

      delete[] recv_int_buf;
      delete[] recv_double_buf;

      double my_errs[2] = {dn2_max,dx2_max};
      double errs[2];
      MPI_Allreduce(my_errs,errs,2,MPI_DOUBLE,MPI_MAX,mpi_comm);
      if (mpi_rank == 0)
        cout << " > averageFaDupData dn2_max: " << errs[0] << " dx2_max: " << errs[1] << endl;

    }

    namespace PeriodicFunc {

      // periodic nodes are arranged in groups at the start of the node list,
      // followed by boundary nodes, followed by regular internal nodes. Note that
      // some boundary nodes can be periodic.
      //
      // in this routine, the 6 bits 0-5 (i.e. 3 bit pairs) are used to determine
      // groups, bit (1<<6) is used to indicate a boundary node.

      void dumpBits(const int bits) {
        for (int i = 0; i <= 6; ++i) {
          if (bits & (1<<i)) cout << "1";
          else cout << "0";
        }
      }

      inline static int getFaceZoneForBits(const int bits) {

        // should be just bit pairs...
        assert(bits <= ((1<<5)|(1<<3)|(1<<1)));

        switch (bits) {
          case (0):
            return 26;
          case (1<<0):
            return 0;
          case (1<<1):
            return 1;
          case (1<<2):
            return 2;
          case (1<<3):
            return 3;
          case (1<<2)|(1<<0):
            return 4;
          case (1<<2)|(1<<1):
            return 5;
          case (1<<3)|(1<<0):
            return 6;
          case (1<<3)|(1<<1):
            return 7;
          case (1<<4):
            return 8;
          case (1<<5):
            return 9;
          case (1<<4)|(1<<0):
            return 10;
          case (1<<4)|(1<<1):
            return 11;
          case (1<<5)|(1<<0):
            return 12;
          case (1<<5)|(1<<1):
            return 13;
          case (1<<4)|(1<<2):
            return 14;
          case (1<<4)|(1<<3):
            return 15;
          case (1<<5)|(1<<2):
            return 16;
          case (1<<5)|(1<<3):
            return 17;
          case (1<<4)|(1<<2)|(1<<0):
            return 18;
          case (1<<4)|(1<<2)|(1<<1):
            return 19;
          case (1<<4)|(1<<3)|(1<<0):
            return 20;
          case (1<<4)|(1<<3)|(1<<1):
            return 21;
          case (1<<5)|(1<<2)|(1<<0):
            return 22;
          case (1<<5)|(1<<2)|(1<<1):
            return 23;
          case (1<<5)|(1<<3)|(1<<0):
            return 24;
          case (1<<5)|(1<<3)|(1<<1):
            return 25;
          default:
            cout << "getBitsFromFaceZone: cannot find bits: " << bits << " ";
            for (int i = 5; i >= 0; --i)
              if (bits & (1<<i))
                cout << "1";
              else
                cout << "0";
            cout << endl;
            throw(0);

            return -1;
        }
      }

      inline static int getBitsForFaceZone(const int zone) {

        switch (zone) {
          case (26):
            return (0);
          case (0):
            return (1<<0);
          case (1):
            return (1<<1);
          case (2):
            return (1<<2);
          case (3):
            return (1<<3);
          case (4):
            return (1<<2)|(1<<0);
          case (5):
            return (1<<2)|(1<<1);
          case (6):
            return (1<<3)|(1<<0);
          case (7):
            return (1<<3)|(1<<1);
          case (8):
            return (1<<4);
          case (9):
            return (1<<5);
          case (10):
            return (1<<4)|(1<<0);
          case (11):
            return (1<<4)|(1<<1);
          case (12):
            return (1<<5)|(1<<0);
          case (13):
            return (1<<5)|(1<<1);
          case (14):
            return (1<<4)|(1<<2);
          case (15):
            return (1<<4)|(1<<3);
          case (16):
            return (1<<5)|(1<<2);
          case (17):
            return (1<<5)|(1<<3);
          case (18):
            return (1<<4)|(1<<2)|(1<<0);
          case (19):
            return (1<<4)|(1<<2)|(1<<1);
          case (20):
            return (1<<4)|(1<<3)|(1<<0);
          case (21):
            return (1<<4)|(1<<3)|(1<<1);
          case (22):
            return (1<<5)|(1<<2)|(1<<0);
          case (23):
            return (1<<5)|(1<<2)|(1<<1);
          case (24):
            return (1<<5)|(1<<3)|(1<<0);
          case (25):
            return (1<<5)|(1<<3)|(1<<1);
          default:
            // should be one of the 27 possible zones (0:26)...
            cout << "getBitsFromFaceZone: cannot find zone: " << zone << endl;
            throw(0);

            return -1;
        }
      }

      inline static int getNodeGroup(const int bits) {

        // node groups are reordered so the largest (8 members) are first,
        // then the 4's, then the 2's. This simplifies some of the bit
        // calculations used

        switch (bits) {
          case 0: // regular internal node
            return 8;

          case (1<<6): // boundary node with no periodicity
            return 7;

            // xyz...
          case (1<<4)|(1<<2)|(1<<0):
          case (1<<4)|(1<<2)|(1<<1):
          case (1<<4)|(1<<3)|(1<<0):
          case (1<<4)|(1<<3)|(1<<1):
          case (1<<5)|(1<<2)|(1<<0):
          case (1<<5)|(1<<2)|(1<<1):
          case (1<<5)|(1<<3)|(1<<0):
          case (1<<5)|(1<<3)|(1<<1):
          case (1<<6)|(1<<4)|(1<<2)|(1<<0):
          case (1<<6)|(1<<4)|(1<<2)|(1<<1):
          case (1<<6)|(1<<4)|(1<<3)|(1<<0):
          case (1<<6)|(1<<4)|(1<<3)|(1<<1):
          case (1<<6)|(1<<5)|(1<<2)|(1<<0):
          case (1<<6)|(1<<5)|(1<<2)|(1<<1):
          case (1<<6)|(1<<5)|(1<<3)|(1<<0):
          case (1<<6)|(1<<5)|(1<<3)|(1<<1):
            return 0;

            // xy...
          case (1<<2)|(1<<0):
          case (1<<2)|(1<<1):
          case (1<<3)|(1<<0):
          case (1<<3)|(1<<1):
          case (1<<6)|(1<<2)|(1<<0):
          case (1<<6)|(1<<2)|(1<<1):
          case (1<<6)|(1<<3)|(1<<0):
          case (1<<6)|(1<<3)|(1<<1):
            return 1;

            // xz...
          case (1<<4)|(1<<0):
          case (1<<4)|(1<<1):
          case (1<<5)|(1<<0):
          case (1<<5)|(1<<1):
          case (1<<6)|(1<<4)|(1<<0):
          case (1<<6)|(1<<4)|(1<<1):
          case (1<<6)|(1<<5)|(1<<0):
          case (1<<6)|(1<<5)|(1<<1):
            return 2;

            // yz...
          case (1<<4)|(1<<2):
          case (1<<4)|(1<<3):
          case (1<<5)|(1<<2):
          case (1<<5)|(1<<3):
          case (1<<6)|(1<<4)|(1<<2):
          case (1<<6)|(1<<4)|(1<<3):
          case (1<<6)|(1<<5)|(1<<2):
          case (1<<6)|(1<<5)|(1<<3):
            return 3;

            // x...
          case (1<<0):
          case (1<<1):
          case (1<<6)|(1<<0): // boundary
          case (1<<6)|(1<<1):
            return 4;

            // y...
          case (1<<2):
          case (1<<3):
          case (1<<6)|(1<<2):
          case (1<<6)|(1<<3):
            return 5;

            // z...
          case (1<<4):
          case (1<<5):
          case (1<<6)|(1<<4):
          case (1<<6)|(1<<5):
            return 6;

          default:
            cout << "cannot find group for bits: " << bits << endl;
            assert(0);
            return -1;
        }
      }

      inline static int getNodeSubGroup(const int bits) {

        switch (bits) {

          case 0:
          case (1<<6):
          case (1<<0):
          case (1<<2):
          case (1<<2)|(1<<0):
          case (1<<4):
          case (1<<4)|(1<<0):
          case (1<<4)|(1<<2):
          case (1<<4)|(1<<2)|(1<<0):
          case (1<<6)|(1<<0):
          case (1<<6)|(1<<2):
          case (1<<6)|(1<<2)|(1<<0):
          case (1<<6)|(1<<4):
          case (1<<6)|(1<<4)|(1<<0):
          case (1<<6)|(1<<4)|(1<<2):
          case (1<<6)|(1<<4)|(1<<2)|(1<<0):
            return 0;

          case (1<<1):
          case (1<<3):
          case (1<<2)|(1<<1):
          case (1<<5):
          case (1<<4)|(1<<1):
          case (1<<4)|(1<<3):
          case (1<<4)|(1<<2)|(1<<1):
          case (1<<6)|(1<<1):
          case (1<<6)|(1<<3):
          case (1<<6)|(1<<2)|(1<<1):
          case (1<<6)|(1<<5):
          case (1<<6)|(1<<4)|(1<<1):
          case (1<<6)|(1<<4)|(1<<3):
          case (1<<6)|(1<<4)|(1<<2)|(1<<1):
            return 1;

          case (1<<3)|(1<<0):
          case (1<<5)|(1<<0):
          case (1<<5)|(1<<2):
          case (1<<4)|(1<<3)|(1<<0):
          case (1<<6)|(1<<3)|(1<<0):
          case (1<<6)|(1<<5)|(1<<0):
          case (1<<6)|(1<<5)|(1<<2):
          case (1<<6)|(1<<4)|(1<<3)|(1<<0):
            return 2;

          case (1<<3)|(1<<1):
          case (1<<5)|(1<<1):
          case (1<<5)|(1<<3):
          case (1<<4)|(1<<3)|(1<<1):
          case (1<<6)|(1<<3)|(1<<1):
          case (1<<6)|(1<<5)|(1<<1):
          case (1<<6)|(1<<5)|(1<<3):
          case (1<<6)|(1<<4)|(1<<3)|(1<<1):
            return 3;

          case (1<<5)|(1<<2)|(1<<0):
          case (1<<6)|(1<<5)|(1<<2)|(1<<0):
            return 4;

          case (1<<5)|(1<<2)|(1<<1):
          case (1<<6)|(1<<5)|(1<<2)|(1<<1):
            return 5;

          case (1<<5)|(1<<3)|(1<<0):
          case (1<<6)|(1<<5)|(1<<3)|(1<<0):
            return 6;

          case (1<<5)|(1<<3)|(1<<1):
          case (1<<6)|(1<<5)|(1<<3)|(1<<1):
            return 7;

          default:
            cout << "cannot find sub group for bits: " << bits << endl;
            assert(0);
            return -1;
        }
      }

      inline static int getNodeGroupSize(const int ing) {

        switch (ing) {
          case 0:
            return 8;
          case 1:
          case 2:
          case 3:
            return 4;
          case 4:
          case 5:
          case 6:
            return 2;
          case 7: // boundary only
          case 8: // internal
            return 1;

          default:
            cout << "cannot find size for node group: " << ing << endl;
            assert(0);
            return -1;
        }
      }

      static int8 calcGlobalNodeNumber(const int8 ino_global_even,const int bits,const int8 count[]) {

        const int ing = getNodeGroup(bits);
        int8 ing_offset = 0;
        for (int ing_prev = 0; ing_prev < ing; ++ing_prev)
          ing_offset += count[ing_prev];

        if (!((ino_global_even >= ing_offset)&&(ino_global_even < ing_offset+count[ing]))) {
          cout << "bits: " << bits << " ing: " << ing << " ino_global_even: " << ino_global_even << " ing_offset: " << ing_offset << endl;
        }
        assert((ino_global_even >= ing_offset)&&(ino_global_even < ing_offset+count[ing]));

        // also make sure the node number is staggered appropriately in the group...
        // it should be the first in a set of getNodeGroupSize...

        const int ing_size = getNodeGroupSize(ing);
        assert((ino_global_even-ing_offset)%ing_size == 0);

        return ino_global_even + getNodeSubGroup(bits);

      }

      static int getFirstSixBitsForNodeGroupMember(const int ing,const int index) {

        switch (ing) {

          case 0: // xyz
            switch (index) {
              case 0:
                return (1<<0)|(1<<2)|(1<<4);
              case 1:
                return (1<<1)|(1<<2)|(1<<4);
              case 2:
                return (1<<0)|(1<<3)|(1<<4);
              case 3:
                return (1<<1)|(1<<3)|(1<<4);
              case 4:
                return (1<<0)|(1<<2)|(1<<5);
              case 5:
                return (1<<1)|(1<<2)|(1<<5);
              case 6:
                return (1<<0)|(1<<3)|(1<<5);
              case 7:
                return (1<<1)|(1<<3)|(1<<5);
              default:
                assert(0);
            }

          case 1: // xy
            switch (index) {
              case 0:
                return (1<<0)|(1<<2);
              case 1:
                return (1<<1)|(1<<2);
              case 2:
                return (1<<0)|(1<<3);
              case 3:
                return (1<<1)|(1<<3);
              default:
                assert(0);
            }

          case 2: // xz
            switch (index) {
              case 0:
                return (1<<0)|(1<<4);
              case 1:
                return (1<<1)|(1<<4);
              case 2:
                return (1<<0)|(1<<5);
              case 3:
                return (1<<1)|(1<<5);
              default:
                assert(0);
            }

          case 3: // yz
            switch (index) {
              case 0:
                return (1<<2)|(1<<4);
              case 1:
                return (1<<3)|(1<<4);
              case 2:
                return (1<<2)|(1<<5);
              case 3:
                return (1<<3)|(1<<5);
              default:
                assert(0);
            }

          case 4: // x
            switch (index) {
              case 0:
                return 1<<0;
              case 1:
                return 1<<1;
              default:
                assert(0);
            }

          case 5: // y
            switch (index) {
              case 0:
                return 1<<2;
              case 1:
                return 1<<3;
              default:
                assert(0);
            }

          case 6: // z
            switch (index) {
              case 0:
                return 1<<4;
              case 1:
                return 1<<5;
              default:
                assert(0);
            }

          case 7: // boundary only
          case 8: // internal
            assert(index == 0);
            return 0;
          default:
            cout << "cannot find node group: " << ing << endl;
            assert(0);
            return -1;
        }

      }

      // test conistency of periodic transforms...
      static void testGroupNumberingRoutines() {

        // consistency of node group numbering...
        for (int ing = 0; ing < 7; ++ing) { // x,y,z,xy,xz,yz,xyz
          const int nm = getNodeGroupSize(ing);
          for (int ii = 0; ii < nm; ++ii) {
            const int bits = getFirstSixBitsForNodeGroupMember(ing,ii);
            assert(ing == getNodeGroup(bits));
          }
        }

      }

    };

    bool secondPairSecondInt(const pair<double*,pair<int,int> >& nv0, const pair<double*,pair<int,int> >& nv1) {
      return (nv0.second.second < nv1.second.second);
    }

    // recursive function to perform stost depth-first comparison
    void checkTriNbrUnitN(vector<VdTri>& vdBfTriVec, vector<int>& vdBfColorVec, int *seed_st, const int (*stost)[3], const int ist, const int nst,
        const int iseed, const int ist_seed, const double n_st_seed[3], const double n_mag_seed, const double dp_tol, int depth) {
      if (depth < 10000) {
        FOR_I3 {
          const int ist_nbr = stost[ist][i];
          // if internal and of the same color
          if ( (ist_nbr >= 0) && (vdBfColorVec[ist] == vdBfColorVec[ist_nbr]) ) {
            // if not linked to seed
            if (seed_st[ist_nbr] < 0) {
              const double n_st_nbr[3] = TRI_NORMAL_2(vdBfTriVec[ist_nbr].x[0],vdBfTriVec[ist_nbr].x[1],vdBfTriVec[ist_nbr].x[2]);
              const double n_mag_nbr = MAG(n_st_nbr);
              if ( (n_mag_nbr <= 1.0E-24) || (DOT_PRODUCT(n_st_nbr,n_st_seed)/(n_mag_nbr*n_mag_seed) > dp_tol) ) {
                seed_st[ist_nbr] = iseed;
                checkTriNbrUnitN(vdBfTriVec,vdBfColorVec,seed_st,stost,ist_nbr,nst,iseed,ist_seed,n_st_seed,n_mag_seed,dp_tol,++depth);
              }
            }
          }
        }
      }
    }

  int *force_planar_zone_flag = NULL;

    int buildBfColorMap(map<const pair<pair<int,int>,int>,int>& colorMap,const int ii) {

      // this routine would be better as a class member in VdArray[ip], except we also
      // need access to the ophans that are in ocdVec[next]...

      vector<pair<int,int> > colorPairVec;
      assert(colorMap.empty());
      int n_color = 0;
      for (int ied = 0; ied < vdReturnDataVec[ii].ned; ++ied) {
        // each edge has 2 faces...
        const int ibf0 = vdReturnDataVec[ii].faoed[ied][0];
        int color0 = -1;
        int ipart0,ist0,bits0;
        if (ibf0 >= 0) {
          ist0 = vdReturnDataVec[ii].spbobf[ibf0][0];
          bits0 = (vdReturnDataVec[ii].spbobf[ibf0][1]&MASK_6BITS);
          ipart0 = (vdReturnDataVec[ii].spbobf[ibf0][1]>>6);
          assert((ipart0 >= 0)&&(ipart0 < partVec.size()));
          map<const pair<pair<int,int>,int>,int>::iterator iter = colorMap.find(pair<pair<int,int>,int>(pair<int,int>(ipart0,ist0),bits0));
          if (iter == colorMap.end()) {
            colorMap[pair<pair<int,int>,int>(pair<int,int>(ipart0,ist0),bits0)] = color0 = n_color++;
          }
          else {
            color0 = iter->second;
          }
        }
        const int ibf1 = vdReturnDataVec[ii].faoed[ied][1];
        int color1 = -1;
        int ipart1,ist1,bits1;
        if (ibf1 >= 0) {
          ist1 = vdReturnDataVec[ii].spbobf[ibf1][0];
          bits1 = (vdReturnDataVec[ii].spbobf[ibf1][1]&MASK_6BITS);
          ipart1 = (vdReturnDataVec[ii].spbobf[ibf1][1]>>6);
          assert((ipart1 >= 0)&&(ipart1 < partVec.size()));
          map<const pair<pair<int,int>,int>,int>::iterator iter = colorMap.find(pair<pair<int,int>,int>(pair<int,int>(ipart1,ist1),bits1));
          if (iter == colorMap.end()) {
            colorMap[pair<pair<int,int>,int>(pair<int,int>(ipart1,ist1),bits1)] = color1 = n_color++;
          }
          else {
            color1 = iter->second;
          }
        }
        // edges between 2 surface tris can combine the 2 colors they touch...
        if ((ibf0 >= 0)&&(ibf1 >= 0)) {
          assert(color0 >= 0);
          assert(color1 >= 0);
          // this is an boundary edge between 2 positive, different boundary tris...
          assert((ipart0 != ipart1)||(ist0 != ist1)||(bits0 != bits1)); // the same tri can be self-periodic if it is in the r=0 corner of a slice
          // if this edge is below the crease angle, record the link between these colors...
          const int izone0 = partVec[ipart0]->znosz[partVec[ipart0]->surface->znost[ist0]];
          const int izone1 = partVec[ipart1]->znosz[partVec[ipart1]->surface->znost[ist1]];
          if (izone0 == izone1) {
	    // do not pair if this zone has FORCE_PLANAR...
	    if ((force_planar_zone_flag == NULL)||(force_planar_zone_flag[izone0] == 0)) {
	      // same zone, so check angle...
	      if (checkCoplanar(ipart0,ist0,bits0,ipart1,ist1,bits1,crease_angle_degrees)) {
		// same zone and coplanar, so traverse to bottom and joint...
		colorPairVec.push_back(pair<int,int>(color0,color1));
	      }
	    }
          }
        }
      }

      int * bf_color = new int[n_color];
      for (int ic = 0; ic < n_color; ++ic) bf_color[ic] = ic;

      for (int ii = 0, limit = colorPairVec.size(); ii < limit; ++ii) {
        int ic0 = bf_color[colorPairVec[ii].first];
        while (ic0 != bf_color[ic0]) ic0 = bf_color[ic0];
        int ic1 = bf_color[colorPairVec[ii].second];
        while (ic1 != bf_color[ic1]) ic1 = bf_color[ic1];
        bf_color[ic0] = bf_color[ic1] = min(ic0,ic1);
      }

      int n_color_final = 0;
      for (int ic = 0; ic < n_color; ++ic) {
        if (bf_color[ic] == ic) {
          ++n_color_final;
          bf_color[ic] = -n_color_final;
        }
        else {
          int ic_ = bf_color[ic];
          while (ic_ >= 0) ic_ = bf_color[ic_];
          bf_color[ic] = ic_;
        }
      }

      // flip back to positive...

      for (int ic = 0; ic < n_color; ++ic) {
        bf_color[ic] = -bf_color[ic]-1;
        assert((bf_color[ic] >= 0)&&(bf_color[ic] < n_color_final));
      }

      // and modify the map...
      for (map<const pair<pair<int,int>,int>,int>::iterator it = colorMap.begin(); it != colorMap.end(); ++it) {
        const int ic = it->second;
        assert((ic >= 0)&&(ic < n_color));
        it->second = bf_color[ic];
      }

      // and check...
      for (int ii = 0, limit = colorPairVec.size(); ii < limit; ++ii) {
        assert(bf_color[colorPairVec[ii].first] == bf_color[colorPairVec[ii].second]);
      }

      delete[] bf_color;

      return n_color_final;

    }

#include "tetIntegrateAbsSignedDist.hpp"

    void writeMles(const string& filename,const int io_version,const bool write_surface) {
      COUT1("writeMles: " << filename << " io_version=" << io_version);

      const double d2_tol = 1.0E-16; // tolerance used in VoronoiBuilder...

      // clean up part memory that is not required...
      for (int ipart = 0; ipart < partVec.size(); ++ipart) {
        if (partVec[ipart]->surface)
          partVec[ipart]->surface->clearTeost();
	if (partVec[ipart]->subSurface) {
	  partVec[ipart]->subSurface->clear();
	  delete partVec[ipart]->subSurface;
	  partVec[ipart]->subSurface = NULL;
	}
      }

      PeriodicFunc::testGroupNumberingRoutines(); // just in case someone changed something

      static int debug_counter = 0;
      const int debug_rank = getIntParam("DEBUG_RANK",-1);
      const int debug_count = getIntParam("DEBUG_COUNT",-1);

      vector<double*> noobf_double_ptr_vec;
      vector<double*> noofa_dup_double_ptr_vec;
      int * bfocv_i = new int[ncv+1]; bfocv_i[0] = 0;
      int * faocv_dup_i = new int[ncv+1]; faocv_dup_i[0] = 0;
      vector<int> noobf_i_vec; noobf_i_vec.push_back(0);
      vector<int> spobf_i_vec; spobf_i_vec.push_back(0);
      vector<pair<uint8,double> > spobf_v_wgt_vec;
      vector<int> sbobf_i_vec; sbobf_i_vec.push_back(0);
      vector<pair<pair<int,int>,int> > psbobf_v_vec; // ipart-ist-bits pair of v_vec
      vector<int> noofa_dup_i_vec; noofa_dup_i_vec.reserve(ncv*6+1); noofa_dup_i_vec.push_back(0);
      vector<FaData> fa_dup_data_vec;
      fa_dup_data_vec.reserve(ncv*6); // is this worth it?
      vector<BfData> bf_data_vec;

      // ==================================================================================
      // step 1 -- build face loops? that include primary nodes and cut orphan data nodes...
      // these are loops associated with the nodes of each primary and cut orphan
      // voronoi diagram...
      // ==================================================================================

      double my_buf[4] = { 0.0, 0.0, 0.0, 0.0 };

      multimap<const pair<uint8,double*>,double*> faMultimap;
      multimap<const pair<int,double*>,double*> bfMultimap;

      int nfa_i0 = 0; // count of winning internal faces (icv0,icv1 on the same processor, icv0<icv1)...
      int nfa_i1 = 0; // count of winning inter-processor faces (icv0 on this rank,icv1 on another rank, rank(icv0) < rank(icv1), with some extra tie-breaking for self-periodic)
      // after this routine these counts are added to get the total ACTUAL faces (before any further reduction due to node compression)...

      const int nzone = zoneVec.size();
      vector<int> bf_zone_of_color_final;
      map<const pair<pair<int,int>,int>,int> colorMap; // boundary face color with ((ipart,ist),bits) key

      // maps for local internal and face geometry associated with a single ip...
      map<const uint8,FaceGeometryData> internalFaceGeometryMap; // map of face geometry with rbi key
      map<const pair<pair<int,int>,int>,BfGeometryData> boundaryFaceGeometryMap; // map of face geometry with ((ipart,ist),bits) key
      vector<BfGeometryData> boundaryFaceDataVec;

      // the local cv count is simply the number of local Voronoi diagrams, i.e. the number of forming points
      assert(vol_cv);
      assert(x_cv);
      double * r_vv = new double[ncv];

      vector<VdTri> vdFaTriVec; // internal faces

      // when this routine computes delta (for 1/delta required by boundary viscous closure), it
      // uses either a n2 algorithm, or CI's retriangulation of approx planar patches
      // approach. We found some vd's with huge nbf counts can take literally hours if you
      // just use the n2 algo blindly...

      const int delta_n2_algo_limit = 500; // not rigorous at this point -- adjust if necessary
      int nst_cap = 5000;
      int nseed_cap = 100;
      map<pair<const double*,const double*>,int> edgeMap; // used to build stost
      vector<VdTri> vdBfTriVec; vdBfTriVec.reserve(nst_cap); // boundary faces tri
      vector<int> vdBfColorVec; vdBfColorVec.reserve(nst_cap); // boundary faces color
      int (*stost)[3] = new int[nst_cap][3]; // st of st (-1 for boundary)
      int *seed_st = new int[nst_cap]; // seed of an st
      double (*x_seed)[3] = new double[nseed_cap][3]; // center point of seed
      vector<int> color_seed; color_seed.reserve(nseed_cap); // color of seed
      vector<int> x1x0_seed_i_vec; x1x0_seed_i_vec.reserve(nseed_cap+1); // boundary of seed edge coord indices
      vector<pair<double*,double*> > x1x0_seed_v_vec; x1x0_seed_v_vec.reserve(nst_cap); // boundary of seed edge coord

      MPI_Sync("start");

      int8 my_tot_tri[2] = { 0 , 0 };
      double my_max_err[4] = { 0.0, 0.0, 0.0, 0.0 };
      assert(vdReturnDataVec.size() == ncv);
      // needs to be icv ordered...
      for (int icv = 0; icv < ncv; ++icv) {
        const int ii = ivrdocv[icv];
        assert(icv == vdReturnDataVec[ii].icv);
        assert(vdReturnDataVec[ii].ngr == 1);

        ++debug_counter;
        const bool debug = ((mpi_rank == debug_rank)&&(debug_counter == debug_count));

        if (debug) vdReturnDataVec[ii].writeTecplot("vd_debug.dat",x_vd[icv]);

        //cout << "WORKING on icv: " << icv << endl;

        // step 1.1: populate face map with edges. We have to do this so we can work
        // on the edges associated with a particular face. Right now, the edges are in
        // no partular order (no ed-of-fa structure maintained for Voronoi, for example)...

        const int n_color_final = buildBfColorMap(colorMap,ii);

        // now lets setup bf_zone_of_color_final
        bf_zone_of_color_final.resize(n_color_final);
        for (int ic = 0; ic < n_color_final; ++ic) bf_zone_of_color_final[ic] = -1;

        // =========================================================
        // now build the internal and boundary face maps...
        // =========================================================

        // check -- it was emptied last time
        assert(faMultimap.empty());
        assert(bfMultimap.empty());

        // and internal and boundary face geometry...

        assert(internalFaceGeometryMap.empty());
        assert(boundaryFaceGeometryMap.empty());
        vol_cv[icv] = 0.0;
        FOR_I3 x_cv[icv][i] = 0.0;
        r_vv[icv] = 0.0;

        // and build the tri vec -- basically the tris that, when integrated with the
        // center point -- which is zero -- could be used to compute the volume of
        // the Voronoi diagram.

        assert(vdFaTriVec.empty());
        assert(vdBfTriVec.empty());
        assert(vdBfColorVec.empty());

        if (debug) cout << "1" << endl;

        // no groups or orphans...
        for (int ied = 0; ied < vdReturnDataVec[ii].ned; ++ied) {
          // set r_vv first (as r2_cv for now)...
          FOR_I2 {
            double * x = vdReturnDataVec[ii].x_no[vdReturnDataVec[ii].nooed[ied][i]];
            r_vv[icv] = max(r_vv[icv],DOT_PRODUCT(x,x));
          }
          // each edge has 2 faces...
          FOR_I2 {
            const int ibf_or_ifa = vdReturnDataVec[ii].faoed[ied][i];
            // a negative ibf_or_ifa indicates an internal nbr (positive is a boundary)...
            if (ibf_or_ifa < 0) {
              // recall faoed uses -1 indexing for faces...
              assert((-ibf_or_ifa-1 >= 0)&&(-ibf_or_ifa-1 < vdReturnDataVec[ii].nfa));
              // if face is active...
              if (vdReturnDataVec[ii].faceIsActive(-ibf_or_ifa-1)) {
                uint8 rbiHash; // rbi uniquely represents the nbr we have been cut against
                if (vdReturnDataVec[ii].cvofa[-ibf_or_ifa-1] < ncv)
                  rbiHash = BitUtils::packRankBitsIndex(mpi_rank,0,vdReturnDataVec[ii].cvofa[-ibf_or_ifa-1]);
                else
                  rbiHash = rbi_g[vdReturnDataVec[ii].cvofa[-ibf_or_ifa-1]-ncv];
                int rank,bits,index; BitUtils::unpackRankBitsIndex(rank,bits,index,rbiHash);
                // recall i = 0,1, so edge goes form node 0->1 or 1->0, depending on the face
                double * x0 = vdReturnDataVec[ii].x_no[vdReturnDataVec[ii].nooed[ied][i]];
                double * x1 = vdReturnDataVec[ii].x_no[vdReturnDataVec[ii].nooed[ied][1-i]];
                faMultimap.insert(pair<pair<uint8,double*>,double*>(pair<uint8,double*>(rbiHash,x0),x1));
                // also, build the geometric stuff...
                map<const uint8,FaceGeometryData>::iterator iter = internalFaceGeometryMap.find(rbiHash);
                if (iter == internalFaceGeometryMap.end()) {
                  internalFaceGeometryMap[rbiHash] = FaceGeometryData(x0);
                }
                else if ((x0 != iter->second.x0)&&(x1 != iter->second.x0)) {
                  assert(iter->second.x0 != NULL);
                  const double this_n[3] = TRI_NORMAL_2(iter->second.x0,x0,x1);
                  const double this_vol = DOT_PRODUCT(iter->second.x0,this_n);
                  vol_cv[icv] += this_vol;
                  vdFaTriVec.push_back(VdTri(iter->second.x0,x0,x1));
                  iter->second.area += this_vol;
                  FOR_J3 {
                    const double tmp = this_vol*(iter->second.x0[j]+x0[j]+x1[j]);
                    iter->second.xc[j] += tmp;
                    x_cv[icv][j] += tmp;
                    iter->second.normal[j] += this_n[j];
                  }
                }
              }

            }
            else {

              // a positive ibf_or_ifa is the ibf index, meaning it is an edge of a boundary face.
              // If the other face is an internal face, and
              // it has a connection with an orphan (i.e. group > 0), we skip it...

              assert((ibf_or_ifa >= 0)&&(ibf_or_ifa < vdReturnDataVec[ii].nbf));
              const int ist = vdReturnDataVec[ii].spbobf[ibf_or_ifa][0];
              const int bits = (vdReturnDataVec[ii].spbobf[ibf_or_ifa][1]&MASK_6BITS);
              const int ipart = (vdReturnDataVec[ii].spbobf[ibf_or_ifa][1]>>6);
              assert((ipart >= 0)&&(ipart < partVec.size()));

              map<const pair<pair<int,int>,int>,int>::iterator iter0 = colorMap.end();
              if (vdReturnDataVec[ii].faoed[ied][1-i] >= 0) {

                // the opposite face is a boundary face too. if the color of the face is
                // the same as our color, also skip...

                const int ibf_opp = vdReturnDataVec[ii].faoed[ied][1-i];
                assert((ibf_opp >= 0)&&(ibf_opp < vdReturnDataVec[ii].nbf));
                const int ist_opp = vdReturnDataVec[ii].spbobf[ibf_opp][0];
                const int bits_opp = (vdReturnDataVec[ii].spbobf[ibf_opp][1]&MASK_6BITS);
                const int ipart_opp = (vdReturnDataVec[ii].spbobf[ibf_opp][1]>>6);
                assert((ipart_opp >= 0)&&(ipart_opp < partVec.size()));

                iter0 = colorMap.find(pair<pair<int,int>,int>(pair<int,int>(ipart,ist),bits));
                assert(iter0 != colorMap.end());
                map<const pair<pair<int,int>,int>,int>::iterator iter1 =
                  colorMap.find(pair<pair<int,int>,int>(pair<int,int>(ipart_opp,ist_opp),bits_opp)); assert(iter1 != colorMap.end());
                if (iter0->second == iter1->second) {
                  // we are about to throw this edge out because it connects 2 regions of the same final color,
                  // but we still need to include it in the ifa-based boundary face data...
                  double * x0 = vdReturnDataVec[ii].x_no[vdReturnDataVec[ii].nooed[ied][i]]; // the coordinates of this node would be x0[0],x0[1],x0[2].
                  double * x1 = vdReturnDataVec[ii].x_no[vdReturnDataVec[ii].nooed[ied][1-i]];
                  map<const pair<pair<int,int>,int>,BfGeometryData>::iterator iter =
                    boundaryFaceGeometryMap.find(pair<pair<int,int>,int>(pair<int,int>(ipart,ist),bits));
                  if (iter == boundaryFaceGeometryMap.end()) {
                    boundaryFaceGeometryMap[pair<pair<int,int>,int>(pair<int,int>(ipart,ist),bits)] = BfGeometryData(x0);
                  }
                  else if ((x0 != iter->second.x0)&&(x1 != iter->second.x0)) {
                    assert(iter->second.x0 != NULL);
                    const double this_n[3] = TRI_NORMAL_2(iter->second.x0,x0,x1);
                    const double this_vol = DOT_PRODUCT(iter->second.x0,this_n);
                    vol_cv[icv] += this_vol;
                    vdBfTriVec.push_back(VdTri(iter->second.x0,x0,x1));
                    vdBfColorVec.push_back(iter0->second);
                    const double this_area = MAG(this_n); // because these cut unique surface tris can only be cut by voronoi faces, they are convex
                    iter->second.area += this_area;
                    FOR_J3 {
                      iter->second.xc[j]     += this_area*(iter->second.x0[j]+x0[j]+x1[j]);
                      x_cv[icv][j]            += this_vol*(iter->second.x0[j]+x0[j]+x1[j]);
                      iter->second.normal[j] += this_n[j];
                    }
                  }
                  continue;
                }
              }

              // if we have not "continue'd" by now, this edge should be included in bfMultimap...
              if (iter0 == colorMap.end()) iter0  = colorMap.find(pair<pair<int,int>,int>(pair<int,int>(ipart,ist),bits));
              assert(iter0 != colorMap.end());
              // recall i = 0,1, so edge goes form node 0->1 or 1->0, depending on the face
              double * x0 = vdReturnDataVec[ii].x_no[vdReturnDataVec[ii].nooed[ied][i]]; // the coordinates of this node would be x0[0],x0[1],x0[2].
              double * x1 = vdReturnDataVec[ii].x_no[vdReturnDataVec[ii].nooed[ied][1-i]];
              //assert(bfMultimap.find(pair<int,double*>(iter0->second,x0)) == bfMultimap.end()); // if you hit this, you probably have a "figure-8"
              const int color_final = iter0->second;
              bfMultimap.insert(pair<pair<int,double*>,double*>(pair<int,double*>(color_final,x0),x1));
              if (bf_zone_of_color_final[color_final] == -1) {
                bf_zone_of_color_final[color_final] = partVec[ipart]->znosz[partVec[ipart]->surface->znost[ist]];
              }
              else {
                assert(bf_zone_of_color_final[color_final] == partVec[ipart]->znosz[partVec[ipart]->surface->znost[ist]]);
              }
              // geometry...
              map<const pair<pair<int,int>,int>,BfGeometryData>::iterator iter =
                boundaryFaceGeometryMap.find(pair<pair<int,int>,int>(pair<int,int>(ipart,ist),bits));
              if (iter == boundaryFaceGeometryMap.end()) {
                boundaryFaceGeometryMap[pair<pair<int,int>,int>(pair<int,int>(ipart,ist),bits)] = BfGeometryData(x0);
              }
              else if ((x0 != iter->second.x0)&&(x1 != iter->second.x0)) {
                assert(iter->second.x0 != NULL);
                const double this_n[3] = TRI_NORMAL_2(iter->second.x0,x0,x1);
                const double this_vol = DOT_PRODUCT(iter->second.x0,this_n);
                vol_cv[icv] += this_vol;
                vdBfTriVec.push_back(VdTri(iter->second.x0,x0,x1));
                vdBfColorVec.push_back(iter0->second);
                const double this_area = MAG(this_n); // because these cut surface tris can only be cut by voronoi faces, they are convex
                iter->second.area += this_area;
                FOR_J3 {
                  iter->second.xc[j]     += this_area*(iter->second.x0[j]+x0[j]+x1[j]);
                  x_cv[icv][j]            += this_vol*(iter->second.x0[j]+x0[j]+x1[j]);
                  iter->second.normal[j] += this_n[j];
                }
              }
            }

          }
        }

        if (debug) cout << "2" << endl;

        // normalize the x_cv and vol_cv for this ip...

        FOR_I3 x_cv[icv][i] = x_cv[icv][i]/(4.0*vol_cv[icv]);
        vol_cv[icv] /= 6.0;
        r_vv[icv] = sqrt(r_vv[icv]);

        //cout << icv << " " << r_vv[icv] << " " << vol_cv[icv] << " " << COUT_VEC(x_cv[icv]) << endl;
        //getchar();

        // check that the volume we get from vdFaTriVec/vdBfTriVec is the same...

        /*
           double vol_check = 0.0;
           for (int jj = 0; jj < vdFaTriVec.size(); ++jj) {
           vol_check += SIGNED_TET_VOLUME_6_VEC(vdFaTriVec[jj].x[0],vdFaTriVec[jj].x[1],vdFaTriVec[jj].x[2]);
           }
           for (int jj = 0; jj < vdBfTriVec.size(); ++jj) {
           vol_check += SIGNED_TET_VOLUME_6_VEC(vdBfTriVec[jj].x[0],vdBfTriVec[jj].x[1],vdBfTriVec[jj].x[2]);
           }
           vol_check /= 6.0;
           cout << "vol_check: " << vol_check << " vol_cv[icv]-vol_check: " << vol_cv[icv]-vol_check << endl;
           getchar();
           */

        // compare to canned geometry routine -- remove in future...
        {

          double vol_cv_check,x_cv_check[3];
          vdReturnDataVec[ii].calcVolumeAndCentroidForGroup(vol_cv_check,x_cv_check,0);

          my_max_err[0] = max(my_max_err[0],fabs(1.0-vol_cv[icv]/vol_cv_check));
          my_max_err[1] = max(my_max_err[1],fabs(x_cv_check[0]-x_cv[icv][0]));
          my_max_err[2] = max(my_max_err[2],fabs(x_cv_check[1]-x_cv[icv][1]));
          my_max_err[3] = max(my_max_err[3],fabs(x_cv_check[2]-x_cv[icv][2]));

        }

        // should we reduce the number of boundary tries for the O(n^2) algo to calc delta
        if (vdBfTriVec.size() > delta_n2_algo_limit) {

          if (debug) cout << "3a" << endl;

          const int nst = vdBfTriVec.size();
          //cout << "nst orig: " << nst << endl;
          if (nst > nst_cap) {
            delete[] stost; stost = new int[nst][3];
            delete[] seed_st; seed_st = new int[nst];
            nst_cap = nst;
          }
          map<pair<const double*,const double*>,int>::iterator iter;
          for (int ist = 0; ist < nst; ++ist) {
            // look for edge pairs...
            FOR_I3 {
              const double * x0 = vdBfTriVec[ist].x[i];
              const double * x1 = vdBfTriVec[ist].x[(i+1)%3];
              const int ied = ist*3+i;
              iter = edgeMap.find(pair<const double*,const double*>(x0,x1));
              if (iter != edgeMap.end()) {
                const int ist_nbr = iter->second/3;
                const int i_nbr = iter->second-ist_nbr*3;
                assert(stost[ist_nbr][i_nbr] == -1);
                stost[ist][i] = iter->second/3;
                stost[ist_nbr][i_nbr] = ist;
                // found a match...
                edgeMap.erase(iter);
              }
              else {
                // not found, so add the reverse edge to edgeMap...
                edgeMap[pair<const double*,const double*>(x1,x0)] = ied;
                stost[ist][i] = -1;
              }
            }
          }
          edgeMap.clear();
          // edgeMap is not empty due the edge boundaries

          // seed index of surface tri
          for (int ist = 0; ist < nst; ++ist) seed_st[ist] = -1;
          int nseed = 0;
          const double dp_tol = cos((180.0-crease_angle_degrees)*M_PI/180.0);
          for (int ist = 0; ist < nst; ++ist) {
            // find next seed
            if (seed_st[ist] < 0) {
              const double n_st_seed[3] = TRI_NORMAL_2(vdBfTriVec[ist].x[0],vdBfTriVec[ist].x[1],vdBfTriVec[ist].x[2]);
              const double n_mag_seed = MAG(n_st_seed);
              if (n_mag_seed > 1.0E-24) {
                seed_st[ist] = nseed;
                color_seed.push_back(vdBfColorVec[ist]);
                checkTriNbrUnitN(vdBfTriVec,vdBfColorVec,seed_st,stost,ist,nst,nseed,ist,n_st_seed,n_mag_seed,dp_tol,0);
                ++nseed;
              }
            }
          }
          //cout << vdBfTriVec.size() << " " << nseed << endl;
          my_tot_tri[0] += vdBfTriVec.size();

          // get counts for bounding edges for each seed...
          x1x0_seed_i_vec.resize(nseed+1);
          for (int iseed = 0; iseed <= nseed; ++iseed) x1x0_seed_i_vec[iseed] = 0;
          for (int ist = 0; ist < nst; ++ist) {
            if (seed_st[ist] >= 0) {
              FOR_I3 {
                // if neighbors have different seed or at boundary
                if ((stost[ist][i] < 0)||((stost[ist][i] >= 0)&&(seed_st[stost[ist][i]] != seed_st[ist])&&(seed_st[stost[ist][i]] >= 0))) {
                  int iseed = seed_st[ist];
                  ++x1x0_seed_i_vec[iseed+1]; // increment seed counts
                }
              }
            }
          }
          // scan...
          for (int iseed = 1; iseed <= nseed; ++iseed) x1x0_seed_i_vec[iseed] += x1x0_seed_i_vec[iseed-1];
          // populate bounding edge points for each seed...
          x1x0_seed_v_vec.resize(x1x0_seed_i_vec[nseed]);
          //cout << "nseed nst" << nseed << " " << nst << endl;
          for (int ist = 0; ist < nst; ++ist) {
            if (seed_st[ist] >= 0) {
              FOR_I3 {
                // if neighbors have different seed or at boundary
                if ((stost[ist][i] < 0)||((stost[ist][i] >= 0)&&(seed_st[stost[ist][i]] != seed_st[ist])&&(seed_st[stost[ist][i]] >= 0))) {
                  int iseed = seed_st[ist];
                  //cout << iseed << " " << ist << endl;
                  double * x0 = vdBfTriVec[ist].x[i];
                  double * x1 = vdBfTriVec[ist].x[(i+1)%3];
                  x1x0_seed_v_vec[x1x0_seed_i_vec[iseed]++] = pair<double*,double*>(x1,x0);
                }
              }
            }
          }
          // rewind...
          for (int iseed = nseed-1; iseed > 0; --iseed) x1x0_seed_i_vec[iseed] = x1x0_seed_i_vec[iseed-1];
          x1x0_seed_i_vec[0] = 0;

          // calc planar area/normal and center point for retriangulation
          if (nseed > nseed_cap) {
            delete[] x_seed; x_seed = new double[nseed][3]; // center point
            nseed_cap = nseed;
          }
          //double (*n_seed)[3] = new double[nseed][3]; // 2x area weighted normal
          for (int iseed = 0; iseed < nseed; ++iseed) {
            FOR_I3 x_seed[iseed][i] = 0.0;
            // edge (x1/x0) of seed
            if (x1x0_seed_i_vec[iseed+1]-x1x0_seed_i_vec[iseed] > 3) {
              double L = 0.0;
              for (int eos = x1x0_seed_i_vec[iseed]; eos != x1x0_seed_i_vec[iseed+1]; ++eos) {
                double * x1 = x1x0_seed_v_vec[eos].first;
                double * x0 = x1x0_seed_v_vec[eos].second;
                const double dl = DIST(x1,x0);
                L += dl;
                FOR_I3 x_seed[iseed][i] += dl*(x1[i]+x0[i]);
              }
              FOR_I3 x_seed[iseed][i] /= (2.0*L);
            }
          }
          vdBfTriVec.clear();
          vdBfColorVec.clear();
          for (int iseed = 0; iseed < nseed; ++iseed) {
            //FOR_I3 n_seed[iseed][i] = 0.0;
            // edge (x1/x0) of seed
            if (x1x0_seed_i_vec[iseed+1]-x1x0_seed_i_vec[iseed] > 3) {
              for (int eos = x1x0_seed_i_vec[iseed]; eos != x1x0_seed_i_vec[iseed+1]; ++eos) {
                double * x1 = x1x0_seed_v_vec[eos].first;
                double * x0 = x1x0_seed_v_vec[eos].second;
                //double this_n[3] = TRI_NORMAL_2(x0,x1,x_seed[iseed]);
                //FOR_I3 n_seed[iseed][i] += this_n[i];
                vdBfTriVec.push_back(VdTri(x0,x1,x_seed[iseed]));
                vdBfColorVec.push_back(color_seed[iseed]);
              }
            }
            else {
              // make sure we get valid index range ...
              if (x1x0_seed_i_vec[iseed+1] - x1x0_seed_i_vec[iseed] == 3) {
                const int eos = x1x0_seed_i_vec[iseed];
                double * x0 = x1x0_seed_v_vec[eos].second;
                double * x1 = x1x0_seed_v_vec[eos+1].second;
                double * x2 = x1x0_seed_v_vec[eos+2].second;
                vdBfTriVec.push_back(VdTri(x0,x1,x2));
                vdBfColorVec.push_back(color_seed[iseed]);
              }
            }
          }
          //cout << "nst red: " << vdBfTriVec.size() << endl;
          my_tot_tri[1] += vdBfTriVec.size();

          // cleanup
          color_seed.clear();
          x1x0_seed_i_vec.clear();
          x1x0_seed_v_vec.clear();

          // reduce the boundary face geometry to the colors...

          boundaryFaceDataVec.resize(n_color_final);
          for (int icolor = 0; icolor < n_color_final; ++icolor)
            boundaryFaceDataVec[icolor].zero();

          for (map<const pair<pair<int,int>,int>,BfGeometryData>::iterator iter =
              boundaryFaceGeometryMap.begin(); iter != boundaryFaceGeometryMap.end(); ++iter) {
            // get the final color associated with this patch of boundary face...
            const pair<pair<int,int>,int> ipart_ist_bits_pair = iter->first;
            map<const pair<pair<int,int>,int>,int>::iterator iter0 = colorMap.find(ipart_ist_bits_pair); assert(iter0 != colorMap.end());
            const int color_final = iter0->second;
            assert((color_final >= 0)&&(color_final < n_color_final));
            // store the part index in index1 for this color. There should never be multiple
            // parts per color...
            if (boundaryFaceDataVec[color_final].index1 == -1)
              boundaryFaceDataVec[color_final].index1 = ipart_ist_bits_pair.first.first;
            assert(boundaryFaceDataVec[color_final].index1 == ipart_ist_bits_pair.first.first);
            // and sum area, area-weighted centroid, and area-magnitude normal...
            boundaryFaceDataVec[color_final].area += iter->second.area;
            FOR_J3 {
              boundaryFaceDataVec[color_final].xc[j] += iter->second.xc[j];
              boundaryFaceDataVec[color_final].normal[j] += iter->second.normal[j];
            }
          }
          for (int ist = 0, limit = vdBfTriVec.size(); ist < limit; ++ist) {
            // also, for each surface tri associated with each boundary face, compute the contribution
            // to the viscous area/delta closure. Compute the normal using the actual surface tri for
            // accuracy/robustness...
            const int color_final = vdBfColorVec[ist];
            const double *x0 = vdBfTriVec[ist].x[0];
            const double n_st[3] = TRI_NORMAL_2(vdBfTriVec[ist].x[0],vdBfTriVec[ist].x[1],vdBfTriVec[ist].x[2]);
            const double n_st_mag = MAG(n_st);
            // only perform calculation on non-zero area tris
            if (n_st_mag > 0.0) {
              // this should be the same (direction) as iter->second.normal
              // recall the first point found on each boundary face is stored as x0. Do the integration
              // of the abs signed distance function to get delta...
              const double zero[3] = { 0.0, 0.0, 0.0 };
              double delta_avg = 0.0;
              for (int ii = 0, limit2 = vdFaTriVec.size(); ii < limit2; ++ii) {
                // compute the signed distance at the corners...
                const double phi0 = -(x0[0]*n_st[0]+x0[1]*n_st[1]+x0[2]*n_st[2])/n_st_mag;
                const double phi1 = ((vdFaTriVec[ii].x[0][0]-x0[0])*n_st[0]+
                    (vdFaTriVec[ii].x[0][1]-x0[1])*n_st[1]+
                    (vdFaTriVec[ii].x[0][2]-x0[2])*n_st[2])/n_st_mag;
                const double phi2 = ((vdFaTriVec[ii].x[1][0]-x0[0])*n_st[0]+
                    (vdFaTriVec[ii].x[1][1]-x0[1])*n_st[1]+
                    (vdFaTriVec[ii].x[1][2]-x0[2])*n_st[2])/n_st_mag;
                const double phi3 = ((vdFaTriVec[ii].x[2][0]-x0[0])*n_st[0]+
                    (vdFaTriVec[ii].x[2][1]-x0[1])*n_st[1]+
                    (vdFaTriVec[ii].x[2][2]-x0[2])*n_st[2])/n_st_mag;
                // note first point is zero...
                delta_avg += tetIntegrateAbsSignedDist(phi0,phi1,phi2,phi3,zero,vdFaTriVec[ii].x[0],vdFaTriVec[ii].x[1],vdFaTriVec[ii].x[2]);
              }
              for (int ii = 0, limit2 = vdBfTriVec.size(); ii < limit2; ++ii) {
                // compute the signed distance at the corners...
                const double phi0 = -(x0[0]*n_st[0]+x0[1]*n_st[1]+x0[2]*n_st[2])/n_st_mag;
                const double phi1 = ((vdBfTriVec[ii].x[0][0]-x0[0])*n_st[0]+
                    (vdBfTriVec[ii].x[0][1]-x0[1])*n_st[1]+
                    (vdBfTriVec[ii].x[0][2]-x0[2])*n_st[2])/n_st_mag;
                const double phi2 = ((vdBfTriVec[ii].x[1][0]-x0[0])*n_st[0]+
                    (vdBfTriVec[ii].x[1][1]-x0[1])*n_st[1]+
                    (vdBfTriVec[ii].x[1][2]-x0[2])*n_st[2])/n_st_mag;
                const double phi3 = ((vdBfTriVec[ii].x[2][0]-x0[0])*n_st[0]+
                    (vdBfTriVec[ii].x[2][1]-x0[1])*n_st[1]+
                    (vdBfTriVec[ii].x[2][2]-x0[2])*n_st[2])/n_st_mag;
                // note first point is zero...
                delta_avg += tetIntegrateAbsSignedDist(phi0,phi1,phi2,phi3,zero,vdBfTriVec[ii].x[0],vdBfTriVec[ii].x[1],vdBfTriVec[ii].x[2]);
              }
              delta_avg /= (vol_cv[icv]*24.0);
              //cout << "flat: st_area " << 0.5*n_st_mag << " delta_avg " << delta_avg << endl;
              boundaryFaceDataVec[color_final].area_over_delta += 0.5*n_st_mag/delta_avg; // the area here has a factor of 2
            }
          }

          //DELETE(n_seed);
        }
        else {

          if (debug) cout << "3b" << endl;

          // -----------------------------------
          // use the "n2" algo...
          // -----------------------------------

          // reduce the boundary face geometry to the colors...

          boundaryFaceDataVec.resize(n_color_final);
          for (int icolor = 0; icolor < n_color_final; ++icolor)
            boundaryFaceDataVec[icolor].zero();

          for (map<const pair<pair<int,int>,int>,BfGeometryData>::iterator iter =
              boundaryFaceGeometryMap.begin(); iter != boundaryFaceGeometryMap.end(); ++iter) {
            // get the final color associated with this patch of boundary face...
            const pair<pair<int,int>,int> ipart_ist_bits_pair = iter->first;
            map<const pair<pair<int,int>,int>,int>::iterator iter0 = colorMap.find(ipart_ist_bits_pair); assert(iter0 != colorMap.end());
            const int color_final = iter0->second;
            assert((color_final >= 0)&&(color_final < n_color_final));
            // store the part index in index1 for this color. There should never be multiple
            // parts per color...
            if (boundaryFaceDataVec[color_final].index1 == -1)
              boundaryFaceDataVec[color_final].index1 = ipart_ist_bits_pair.first.first;
            assert(boundaryFaceDataVec[color_final].index1 == ipart_ist_bits_pair.first.first);
            // and sum area, area-weighted centroid, and area-magnitude normal...
            boundaryFaceDataVec[color_final].area += iter->second.area;
            FOR_J3 {
              boundaryFaceDataVec[color_final].xc[j] += iter->second.xc[j];
              boundaryFaceDataVec[color_final].normal[j] += iter->second.normal[j];
            }
            // also, for each surface tri associated with each boundary face, compute the contribution
            // to the viscous area/delta closure...
            const double * const n_st = iter->second.normal;
            const double n_st_mag = MAG(n_st);
            if (n_st_mag > 0.0) {
              // this should be the same (direction) as iter->second.normal
              // recall the first point found on each boundary face is stored as x0. Do the integration
              // of the abs signed distance function to get delta...
              const double zero[3] = { 0.0, 0.0, 0.0 };
              double delta_avg = 0.0;
              for (int ii = 0; ii < vdFaTriVec.size(); ++ii) {
                // compute the signed distance at the corners...
                const double phi0 = -(iter->second.x0[0]*n_st[0]+iter->second.x0[1]*n_st[1]+iter->second.x0[2]*n_st[2])/n_st_mag;
                const double phi1 = ((vdFaTriVec[ii].x[0][0]-iter->second.x0[0])*n_st[0]+
                    (vdFaTriVec[ii].x[0][1]-iter->second.x0[1])*n_st[1]+
                    (vdFaTriVec[ii].x[0][2]-iter->second.x0[2])*n_st[2])/n_st_mag;
                const double phi2 = ((vdFaTriVec[ii].x[1][0]-iter->second.x0[0])*n_st[0]+
                    (vdFaTriVec[ii].x[1][1]-iter->second.x0[1])*n_st[1]+
                    (vdFaTriVec[ii].x[1][2]-iter->second.x0[2])*n_st[2])/n_st_mag;
                const double phi3 = ((vdFaTriVec[ii].x[2][0]-iter->second.x0[0])*n_st[0]+
                    (vdFaTriVec[ii].x[2][1]-iter->second.x0[1])*n_st[1]+
                    (vdFaTriVec[ii].x[2][2]-iter->second.x0[2])*n_st[2])/n_st_mag;
                // note first point is zero...
                delta_avg += tetIntegrateAbsSignedDist(phi0,phi1,phi2,phi3,zero,vdFaTriVec[ii].x[0],vdFaTriVec[ii].x[1],vdFaTriVec[ii].x[2]);
              }
              for (int ii = 0; ii < vdBfTriVec.size(); ++ii) {
                // compute the signed distance at the corners...
                const double phi0 = -(iter->second.x0[0]*n_st[0]+iter->second.x0[1]*n_st[1]+iter->second.x0[2]*n_st[2])/n_st_mag;
                const double phi1 = ((vdBfTriVec[ii].x[0][0]-iter->second.x0[0])*n_st[0]+
                    (vdBfTriVec[ii].x[0][1]-iter->second.x0[1])*n_st[1]+
                    (vdBfTriVec[ii].x[0][2]-iter->second.x0[2])*n_st[2])/n_st_mag;
                const double phi2 = ((vdBfTriVec[ii].x[1][0]-iter->second.x0[0])*n_st[0]+
                    (vdBfTriVec[ii].x[1][1]-iter->second.x0[1])*n_st[1]+
                    (vdBfTriVec[ii].x[1][2]-iter->second.x0[2])*n_st[2])/n_st_mag;
                const double phi3 = ((vdBfTriVec[ii].x[2][0]-iter->second.x0[0])*n_st[0]+
                    (vdBfTriVec[ii].x[2][1]-iter->second.x0[1])*n_st[1]+
                    (vdBfTriVec[ii].x[2][2]-iter->second.x0[2])*n_st[2])/n_st_mag;
                // note first point is zero...
                delta_avg += tetIntegrateAbsSignedDist(phi0,phi1,phi2,phi3,zero,vdBfTriVec[ii].x[0],vdBfTriVec[ii].x[1],vdBfTriVec[ii].x[2]);
              }
              delta_avg /= (vol_cv[icv]*24.0);
              //cout << "orig: st_area " << 0.5*n_st_mag << " delta_avg " << delta_avg << endl;
              boundaryFaceDataVec[color_final].area_over_delta += 0.5*iter->second.area/delta_avg; // the area here has a factor of 2
            }
          }

        } // vdBfTriVec.size() > ?

        if (debug) cout << "4" << endl;

        // ================================================================================
        // for the cv-based gradient approximation, we need to compute the part of the
        // cv-grad tensor associated with boundary faces. It has to be done here while
        // we still have access to the full surface.
        //
        // also add spobf_i/v/wgt stuff here...
        // ================================================================================

        map<const pair<int,int>,double> * spobfMap = new map<const pair<int,int>,double>[n_color_final];

        for (map<const pair<pair<int,int>,int>,BfGeometryData>::iterator iter = boundaryFaceGeometryMap.begin();
            iter != boundaryFaceGeometryMap.end(); ++iter) {
          // get the final color associated with this patch of boundary face...
          const pair<pair<int,int>,int> ipart_ist_bits_pair = iter->first;
          map<const pair<pair<int,int>,int>,int>::iterator iter0 = colorMap.find(ipart_ist_bits_pair); assert(iter0 != colorMap.end());
          const int color_final = iter0->second;
          assert((color_final >= 0)&&(color_final < n_color_final));
          // note that the normal has area magnitude here, and normal and area are basically the same,
          // but we need the unit normal in the bf direction to get the projection...
          if (iter->second.area > 0.0) {
            double projected_area = 0.0;
            double fa_area    = 0.0;
            double * two_n_bf = iter->second.normal;
            double x_bf[3]; FOR_I3 x_bf[i] = iter->second.xc[i]/(3.0*iter->second.area);
            for (map<const uint8,FaceGeometryData>::iterator iter_fgd = internalFaceGeometryMap.begin(); iter_fgd != internalFaceGeometryMap.end(); ++iter_fgd) {
              double * two_n_fa = iter_fgd->second.normal;
              projected_area += abs(DOT_PRODUCT(two_n_fa,two_n_bf));
              fa_area        += sqrt(max(0.0,DOT_PRODUCT(two_n_fa,two_n_fa)));
            }
            // in some cases this projected area will be less than the area of the face. It can
            // even be zero -- i.e. there are no faces with a component in the direction of the
            // projection. so we scale the correction of the A matrix with this area ratio...
            const double r[3] = DIFF(x_bf,x_cv[icv]); // note x_cv is still in the Voronoi reference frame here
            //FOR_I3 FOR_J3 boundaryFaceDataVec[color_final].Gij[i][j] -= r[j]*0.5*two_n_bf[i]*min(1.0,projected_area/(iter->second.area*iter->second.area));
            //FOR_I3 FOR_J3 boundaryFaceDataVec[color_final].Gij[i][j] -= r[j]*0.5*two_n_bf[i]*min(1.0,projected_area/(iter->second.area*fa_area));

            // no projection handled here -- we are going to store the full Gij.  the singularities that
            // might arise are handled in the computation of the operator in StripedMesh/StaticSolver downstream.

            FOR_I3 FOR_J3 boundaryFaceDataVec[color_final].Gij[i][j] -= r[j]*0.5*two_n_bf[i];

            // ========================================
            // now the spobf_i/v/wgt stuff...
            //
            // note that we are using area*x to calculate the centroid of the bf patch. This is
            // fine until the patch is cut by another patch and the bf is no longer convex. In
            // this case, we will need to use n-dot-nbf to get the signed area to ensure we
            // account for the possible non-convex bf...
            // ========================================
            const int ipart = ipart_ist_bits_pair.first.first;
            const int ist = ipart_ist_bits_pair.first.second;
            const int bits = ipart_ist_bits_pair.second;
            const int isp0 = partVec[ipart]->surface->spost[ist][0];
            const int isp1 = partVec[ipart]->surface->spost[ist][1];
            const int isp2 = partVec[ipart]->surface->spost[ist][2];
            // build the barycentric weights from cross product magnitudes...
            double dx0[3];
            double dx1[3];
            double dx2[3];
            if (bits == 0) {
              FOR_I3 dx0[i] = partVec[ipart]->surface->xsp[isp0][i]-x_vd[icv][i] - x_bf[i];
              FOR_I3 dx1[i] = partVec[ipart]->surface->xsp[isp1][i]-x_vd[icv][i] - x_bf[i];
              FOR_I3 dx2[i] = partVec[ipart]->surface->xsp[isp2][i]-x_vd[icv][i] - x_bf[i];
            }
            else {
              // need to apply periodic transforms here to get gometry correct...
              FOR_I3 dx0[i] = partVec[ipart]->surface->xsp[isp0][i];
              PeriodicData::periodicTranslate(dx0,1,bits);
              FOR_I3 dx0[i] -= x_vd[icv][i] + x_bf[i];
              FOR_I3 dx1[i] = partVec[ipart]->surface->xsp[isp1][i];
              PeriodicData::periodicTranslate(dx1,1,bits);
              FOR_I3 dx1[i] -= x_vd[icv][i] + x_bf[i];
              FOR_I3 dx2[i] = partVec[ipart]->surface->xsp[isp2][i];
              PeriodicData::periodicTranslate(dx2,1,bits);
              FOR_I3 dx2[i] -= x_vd[icv][i] + x_bf[i];
            }
            const double w0 = CROSS_DOT(dx1,dx2,iter->second.normal);
            const double w1 = CROSS_DOT(dx2,dx0,iter->second.normal);
            const double w2 = CROSS_DOT(dx0,dx1,iter->second.normal);
            const double wsum = w0 + w1 + w2;
            // add weight*area to each isp...
            // isp0...
            //map<const pair<int,int>,double>::iterator witer = spobfMap[color_final].find(pair<int,int>(isp0,bits));
            map<const pair<int,int>,double>::iterator witer = spobfMap[color_final].find(pair<int,int>(getFlattenedSp(ipart,isp0),bits));
            if (witer == spobfMap[color_final].end())
              spobfMap[color_final][pair<int,int>(getFlattenedSp(ipart,isp0),bits)] = iter->second.area*w0/wsum;
            else
              witer->second += iter->second.area*w0/wsum;
            // isp0...
            //witer = spobfMap[color_final].find(pair<int,int>(isp1,bits));
            witer = spobfMap[color_final].find(pair<int,int>(getFlattenedSp(ipart,isp1),bits));
            if (witer == spobfMap[color_final].end())
              spobfMap[color_final][pair<int,int>(getFlattenedSp(ipart,isp1),bits)] = iter->second.area*w1/wsum;
            else
              witer->second += iter->second.area*w1/wsum;
            // isp2...
            //witer = spobfMap[color_final].find(pair<int,int>(isp2,bits));
            witer = spobfMap[color_final].find(pair<int,int>(getFlattenedSp(ipart,isp2),bits));
            if (witer == spobfMap[color_final].end())
              spobfMap[color_final][pair<int,int>(getFlattenedSp(ipart,isp2),bits)] = iter->second.area*w2/wsum;
            else
              witer->second += iter->second.area*w2/wsum;
          }
        }

        if (debug) cout << "5" << endl;

        // now convert the spobfMap into CSR spobf_i/v/wgt...

        const int nbf0 = spobf_i_vec.size()-1;
        spobf_i_vec.resize(nbf0+1+n_color_final);
        const int spobf_s0 = spobf_v_wgt_vec.size(); assert(spobf_s0 == spobf_i_vec[nbf0]);
        int spobf_s = spobf_s0;
        for (int icolor = 0; icolor < n_color_final; ++icolor)
          spobf_s += spobfMap[icolor].size();
        spobf_v_wgt_vec.resize(spobf_s);
        spobf_s = spobf_s0;
        for (int icolor = 0; icolor < n_color_final; ++icolor) {
          assert(boundaryFaceDataVec[icolor].area > 0.0);
          const double inv_area = 1.0/boundaryFaceDataVec[icolor].area;
          for (map<const pair<int,int>,double>::const_iterator iter = spobfMap[icolor].begin(); iter != spobfMap[icolor].end(); ++iter) {
            const int isp = iter->first.first;
            const int bits = iter->first.second;
            const double wgt = iter->second;
            spobf_v_wgt_vec[spobf_s++] = pair<uint8,double>((uint8(bits)<<52)|uint8(isp),wgt*inv_area);
          }
          spobf_i_vec[nbf0+1+icolor] = spobf_s;
        }
        assert(spobf_s == spobf_v_wgt_vec.size());
        delete[] spobfMap;

        if (debug) cout << "6" << endl;

        // we need to store sbobf_i/v in order to index the surface using our bf this is a local bf to global st csr relationship
        assert(nbf0 == sbobf_i_vec.size()-1);
        const int size0 = psbobf_v_vec.size(); assert(size0 == sbobf_i_vec[nbf0]);
        sbobf_i_vec.resize(nbf0+1+n_color_final);
        for (int icolor = 0; icolor < n_color_final; ++icolor) sbobf_i_vec[nbf0+1+icolor] = 0; // use sbobf_i_vec to store the count first
        // count...
        for (map<const pair<pair<int,int>,int>,BfGeometryData>::iterator iter =
            boundaryFaceGeometryMap.begin(); iter != boundaryFaceGeometryMap.end(); ++iter) {
          const pair<pair<int,int>,int> ipart_ist_bits_pair = iter->first;
          map<const pair<pair<int,int>,int>,int>::iterator iter0 = colorMap.find(ipart_ist_bits_pair);
          assert(iter0 != colorMap.end());
          const int color_final = iter0->second;
          assert((color_final >= 0)&&(color_final < n_color_final));
          ++sbobf_i_vec[nbf0+1+color_final];
        }
        // scan...
        for (int icolor = 0; icolor < n_color_final; ++icolor) sbobf_i_vec[nbf0+1+icolor] += sbobf_i_vec[nbf0+icolor];
        psbobf_v_vec.resize(sbobf_i_vec[nbf0+n_color_final]);
        for (int sob = size0; sob < psbobf_v_vec.size(); ++sob)
          psbobf_v_vec[sob].first.first = -1;
        // populate...
        for (map<const pair<pair<int,int>,int>,BfGeometryData>::iterator iter =
            boundaryFaceGeometryMap.begin(); iter != boundaryFaceGeometryMap.end(); ++iter) {
          const pair<pair<int,int>,int> ipart_ist_bits_pair = iter->first;
          map<const pair<pair<int,int>,int>,int>::iterator iter0 = colorMap.find(ipart_ist_bits_pair);
          assert(iter0 != colorMap.end());
          const int color_final = iter0->second;
          assert((color_final >= 0)&&(color_final < n_color_final));
          assert(psbobf_v_vec[sbobf_i_vec[nbf0+color_final]].first.first == -1);
          psbobf_v_vec[sbobf_i_vec[nbf0+color_final]++] = ipart_ist_bits_pair;
        }
        // rewind...
        for (int icolor = n_color_final-1; icolor > 0; --icolor) sbobf_i_vec[nbf0+icolor] = sbobf_i_vec[nbf0+icolor-1];
        sbobf_i_vec[nbf0] = size0;
        for (int sob = size0; sob < psbobf_v_vec.size(); ++sob) {
          // first.first is ipart...
          assert((psbobf_v_vec[sob].first.first >= 0)&&(psbobf_v_vec[sob].first.first < partVec.size()));
          assert((psbobf_v_vec[sob].second >= 0)&&(psbobf_v_vec[sob].second < (1<<6)));
        }

        if (debug) cout << "7" << endl;

        // cleanup
        boundaryFaceGeometryMap.clear();
        colorMap.clear();
        vdFaTriVec.clear();
        vdBfTriVec.clear();
        vdBfColorVec.clear();

        for (int icolor = 0; icolor < n_color_final; ++icolor) boundaryFaceDataVec[icolor].complete();

        // some cv-based geometry. Here we simply add the Vornooi point to x_cv...

        FOR_I3 x_cv[icv][i] += x_vd[icv][i];

        /*
           cout << "icv: " << icv << " x_cv: " << COUT_VEC(x_cv[icv]) << endl;
           for (int icolor = 0; icolor < n_color_final; ++icolor) {
           cout << "boundaryFace normal: " << COUT_VEC(boundaryFaceDataVec[icolor].normal) <<
           " boundaryFace xc: " << COUT_VEC(boundaryFaceDataVec[icolor].xc) <<
           " area: " << boundaryFaceDataVec[icolor].area <<
           " area_over_delta: " << boundaryFaceDataVec[icolor].area_over_delta <<
           " zone: " << surface->zoneVec[bf_zone_of_color_final[icolor]].getName() << endl;
           }
           getchar();
           */

        // ===================================================
        // at this point the bfMultimap and faMultimap has been filled with the node
        // coordinate pointers for the faces associated with voronoi diagram ip.
        // Also the bf_zone_of_color_final should have been set to the correct zone.
        //
        // pull edges from map and organize...
        // ===================================================

        if (debug) cout << "8" << endl;

        // first work with 1 boundary bf at a time...

        {

          int color_check = -1;

          multimap<const pair<int,double*>,double*>::iterator iter_end = bfMultimap.begin();
          while (iter_end != bfMultimap.end()) {

            // the previous end becomes the begin...
            multimap<const pair<int,double*>,double*>::iterator iter_begin = iter_end;
            assert(iter_begin == bfMultimap.begin());

            const int this_color = iter_begin->first.first; // -1 indexed, but leave in this form
            assert((this_color >= 0)&&(this_color < n_color_final));
            assert(bf_zone_of_color_final[this_color] >= 0);

            // ensure monotonic continuous color...
            assert(this_color == color_check+1);
            ++color_check;

            bf_data_vec.push_back(BfData(bf_zone_of_color_final[this_color],boundaryFaceDataVec[this_color]));

            // advance iter_end to the end of this bf's edges...
            while ((iter_end != bfMultimap.end())&&(iter_end->first.first == this_color)) ++iter_end;

            // We now have the bfMultimap members associated with a particular color (one boundary region
            // with a common zone and basically planar) between iter_begin and iter_end...

            // the list of neighbors here should form 1 or more continuous links...
            multiset<double*> x0Multiset;
            for (multimap<const pair<int,double*>,double*>::iterator iter = iter_begin; iter != iter_end; ++iter) {
              // confirm we are associated with the current bf...
              assert(iter->first.first == iter_begin->first.first);
              // put x0 in the set (or take out if already in)...
              double * x0 = iter->first.second;
              //multiset<double*>::iterator iter_x0 = x0Multiset.find(x0);
              //assert(iter_x0 == x0Multiset.end()); // x0's (and x1's) should be unique: MAY fail here -- use multimap!!!
              x0Multiset.insert(x0);
            }
            multiset<double*> x1Multiset;
            for (multimap<const pair<int,double*>,double*>::iterator iter = iter_begin; iter != iter_end; ++iter) {
              double * x1 = iter->second;
              set<double*>::iterator iter_x0 = x0Multiset.find(x1);
              if (iter_x0 == x0Multiset.end()) {
                x1Multiset.insert(x1);
              }
              else {
                x0Multiset.erase(iter_x0);
              }
            }
            assert(x0Multiset.size() == x1Multiset.size());
            if (x0Multiset.empty()) {

              // this is the normal case. All edges close, so building the edge loop (or loops) is easy...

              while (bfMultimap.begin() != iter_end) {
                multimap<const pair<int,double*>,double*>::iterator iter_first = bfMultimap.begin();
                assert(iter_first != bfMultimap.end());
                double * x0 = iter_first->first.second;
                noobf_double_ptr_vec.push_back(x0);
                //cout << "starting at x0: " << x0 << endl;
                double * x1 = iter_first->second;
                bfMultimap.erase(iter_first);
                while (x1 != x0) {
                  //cout << "next node is: " << x1 << endl;
                  noobf_double_ptr_vec.push_back(x1);
                  iter_first = bfMultimap.find(pair<int,double*>(this_color,x1));
                  assert(iter_first != bfMultimap.end());
                  x1 = iter_first->second;
                  bfMultimap.erase(iter_first);
                }
                noobf_double_ptr_vec.push_back(NULL);
              }

            }
            else {

              // we need to pair the x0 and x1 with additional bfMultimap members...

              multimap<double*,double*> x1x0Multimap;
              double max_winner = 0.0;
              double min_loser = HUGE_VAL;
              for (multiset<double*>::iterator iter_x0 = x0Multiset.begin(); iter_x0 != x0Multiset.end(); ++iter_x0) {
                if (debug) cout << "working on iter_x0: " << *iter_x0 << endl;
                multiset<double*>::iterator iter_x1_closest = x1Multiset.end();
                double closest_d2 = HUGE_VAL;
                for (multiset<double*>::iterator iter_x1 = x1Multiset.begin(); iter_x1 != x1Multiset.end(); ++iter_x1) {
                  const double this_d2 = DIST2(*iter_x0,*iter_x1);
                  if (debug) cout << " > comparing to iter_x1: " << *iter_x1 << " dist: " << sqrt(this_d2) << endl;
                  if (this_d2 < closest_d2) {
                    if (closest_d2 < min_loser) min_loser = closest_d2;
                    closest_d2 = this_d2;
                    iter_x1_closest = iter_x1;
                  }
                  else if (this_d2 < min_loser) {
                    min_loser = this_d2;
                  }
                }
                assert(iter_x1_closest != x1Multiset.end());
                if (closest_d2 > max_winner) {
                  max_winner = closest_d2;
                }
                // store to close bfMultimap
                //x1x0Map[*iter_x1_closest] = *iter_x0;
                x1x0Multimap.insert(pair<double*,double*>(*iter_x1_closest,*iter_x0));
                assert(*iter_x0 != *iter_x1_closest);
                x1Multiset.erase(iter_x1_closest);
              }
              assert(x1Multiset.empty());
              if (min_loser < 2.0*max_winner) {
                cout << "About to fail bf loser/winner check: rank=" << mpi_rank << " icv: " << icv << " debug_counter: " << debug_counter << endl;
                cout.flush();
                if (debug) vdReturnDataVec[ii].writeTecplot("vd.dat",x_vd[icv]);
              }
              assert(min_loser >= 2.0*max_winner);
              my_buf[0] = max(my_buf[0],max_winner);
              my_buf[1] = max(my_buf[1],max_winner/min_loser);

              if (debug) cout << "The largest winner, " << max_winner << ", should be much less than the smallest loser, " << min_loser << endl;

              while (bfMultimap.begin() != iter_end) {
                multimap<const pair<int,double*>,double*>::iterator iter_first = bfMultimap.begin();
                assert(iter_first != bfMultimap.end());
                double * x0 = iter_first->first.second;
                noobf_double_ptr_vec.push_back(x0);
                double * x1 = iter_first->second;
                bfMultimap.erase(iter_first);
                //cout << "starting at x0: " << x0 << endl;
                while (x1 != x0) {
                  //cout << "next node is: " << iter_first->second << endl;
                  noobf_double_ptr_vec.push_back(x1);
                  iter_first = bfMultimap.find(pair<int,double*>(this_color,x1));
                  if (iter_first == bfMultimap.end()) {
                    // hop the bridge b/w main and orphan
                    multimap<double*,double*>::iterator iter_x0 = x1x0Multimap.find(x1);
                    assert(iter_x0 != x1x0Multimap.end());
                    x1 = iter_x0->second;
                    x1x0Multimap.erase(iter_x0);
                  }
                  else {
                    x1 = iter_first->second;
                    bfMultimap.erase(iter_first);
                  }
                }
                noobf_double_ptr_vec.push_back(NULL);
              }
              assert(x1x0Multimap.empty());

            }

            noobf_i_vec.push_back(noobf_double_ptr_vec.size());

          } // bfMultimap loop

          // in either case (orphan or not) we should be empty

          assert(bfMultimap.empty());

          bfocv_i[icv+1] = noobf_i_vec.size()-1;

        }

        if (debug) cout << "9" << endl;

        // ==================================================================
        // now work the internal face work...
        // ==================================================================

        {

          map<const uint8,FaceGeometryData>::iterator iter_fgd = internalFaceGeometryMap.begin();

          map<const pair<uint8,double*>,double*>::iterator iter_end = faMultimap.begin();
          while (iter_end != faMultimap.end()) {

            // the previous end becomes the begin...
            map<const pair<uint8,double*>,double*>::iterator iter_begin = iter_end;
            assert(iter_begin == faMultimap.begin());

            const uint8 this_rbiHash = iter_begin->first.first;

            // the internalFaceGeometryMap should be sorted on rbi as well...
            assert(iter_fgd->first ==  this_rbiHash);
            iter_fgd->second.complete(); // normalize
            fa_dup_data_vec.push_back(FaData(this_rbiHash,iter_fgd->second));
            ++iter_fgd; // for the next time
            //rbi_nbr_of_fa_dup_vec.push_back(this_rbiHash);

            // we can also count the ACTUAL faces here...
            int rank,bits,index;
            BitUtils::unpackRankBitsIndex(rank,bits,index,this_rbiHash);
            if ( (rank == mpi_rank)&&(bits == 0)&&(icv < index) ) {
              ++nfa_i0;
            }
            else if ((bits > BitUtils::flipPeriodicBits(bits)) || ((bits == 0)&&(mpi_rank < rank))) {
              ++nfa_i1;
            }

            // advance iter_end to the end of this face's edges...
            while ((iter_end != faMultimap.end())&&(iter_end->first.first == this_rbiHash)) ++iter_end;

            // We now have the faMultimap members associated with a particular rbiHash (one nbr) between
            // iter_begin and iter_end...

            // the list of neighbors here should form 1 or more continuous links...
            multiset<double*> x0Multiset;
            for (multimap<const pair<uint8,double*>,double*>::iterator iter = iter_begin; iter != iter_end; ++iter) {
              // confirm we are associated with the current face...
              assert(iter->first.first == iter_begin->first.first);
              // put x0 in the set (or take out if already in)...
              double * x0 = iter->first.second;
              //set<double*>::iterator iter_x0 = x0Set.find(x0);
              //assert(iter_x0 == x0Set.end()); // x0's (and x1's) should be unique...
              x0Multiset.insert(x0);
            }
            multiset<double*> x1Multiset;
            for (multimap<const pair<uint8,double*>,double*>::iterator iter = iter_begin; iter != iter_end; ++iter) {
              double * x1 = iter->second;
              multiset<double*>::iterator iter_x0 = x0Multiset.find(x1);
              if (iter_x0 == x0Multiset.end()) {
                x1Multiset.insert(x1);
              }
              else {
                x0Multiset.erase(iter_x0);
              }
            }

            assert(x0Multiset.size() == x1Multiset.size());
            if (x0Multiset.size() == 0) {

              // this is the normal case. All edges close, so building the edge loop (or loops) is easy...

              double* x1_prev = NULL;
              while (faMultimap.begin() != iter_end) {
                // start at repeated key if it exists...
                multimap<const pair<uint8,double*>,double*>::iterator iter_first;
                if (x1_prev == NULL) {
                  iter_first = faMultimap.begin();
                  while (iter_first != iter_end) {
                    multimap<const pair<uint8,double*>,double*>::iterator iter_next = iter_first; ++iter_next;
                    if ((iter_next != iter_end)&&(iter_next->first == iter_first->first))
                      break;
                    ++iter_first;
                  }
                  // start at beginning if key DNE...
                  if (iter_first == iter_end)
                    iter_first = faMultimap.begin();
                }
                else {
                  iter_first = faMultimap.find(pair<uint8,double*>(this_rbiHash,x1_prev));
                  // start at beginning if key DNE...
                  if (iter_first == faMultimap.end())
                    iter_first = faMultimap.begin();
                }
                double * x0 = iter_first->first.second;
                if (x0 == x1_prev)
                  noofa_dup_double_ptr_vec.back() = x0; // figure 8
                else
                  noofa_dup_double_ptr_vec.push_back(x0);
                double * x1 = iter_first->second;
                faMultimap.erase(iter_first);
                while (x1 != x0) {
                  noofa_dup_double_ptr_vec.push_back(x1);
                  iter_first = faMultimap.find(pair<uint8,double*>(this_rbiHash,x1));
                  assert(iter_first != faMultimap.end());
                  x1 = iter_first->second;
                  faMultimap.erase(iter_first);
                }
                noofa_dup_double_ptr_vec.push_back(NULL);
                x1_prev = x1;
              }

            }
            else {

              // we need to pair the x0 and x1 with additional faMultimap members...

              multimap<double*,double*> x1x0Multimap;
              double max_winner = 0.0;
              double min_loser = HUGE_VAL;
              for (multiset<double*>::iterator iter_x0 = x0Multiset.begin(); iter_x0 != x0Multiset.end(); ++iter_x0) {
                if (debug) cout << "working on iter_x0: " << *iter_x0 << endl;
                set<double*>::iterator iter_x1_closest = x1Multiset.end();
                double closest_d2 = HUGE_VAL;
                for (multiset<double*>::iterator iter_x1 = x1Multiset.begin(); iter_x1 != x1Multiset.end(); ++iter_x1) {
                  const double this_d2 = DIST2(*iter_x0,*iter_x1);
                  if (debug) cout << " > comparing to iter_x1: " << *iter_x1 << " dist: " << sqrt(this_d2) << endl;
                  if (this_d2 < closest_d2) {
                    if (closest_d2 < min_loser) min_loser = closest_d2;
                    closest_d2 = this_d2;
                    iter_x1_closest = iter_x1;
                  }
                  else if (this_d2 < min_loser) {
                    min_loser = this_d2;
                  }
                }
                assert(iter_x1_closest != x1Multiset.end());
                if (closest_d2 > max_winner) {
                  max_winner = closest_d2;
                }
                // store to close faMultimap
                x1x0Multimap.insert(pair<double*,double*>(*iter_x1_closest,*iter_x0));
                assert(*iter_x0 != *iter_x1_closest);
                x1Multiset.erase(iter_x1_closest);
              }
              assert(x1Multiset.empty());
              if (min_loser < 2.0*max_winner) {
                if (min_loser < 2.0*max_winner) {
                  cout << "About to fail fa loser/winner check: rank=" << mpi_rank << " icv: " << icv << " debug_counter: " << debug_counter << endl;
                  cout.flush();
                  if (debug) vdReturnDataVec[ii].writeTecplot("vd.dat",x_vd[icv]);
                }
              }
              assert(min_loser >= 2.0*max_winner);
              my_buf[2] = max(my_buf[2],max_winner);
              my_buf[3] = max(my_buf[3],max_winner/min_loser);

              if (debug) cout << "The largest winner, " << max_winner << ", should be much less than the smallest loser, " << min_loser << endl;

              //for (map<double*,double*>::iterator iter = x1x0Multimap.begin(); iter != x1x0Multimap.end(); ++iter) {
              //  cout << "Proximity pair (x1,x0): (" << iter->first << "," << iter->second << ") = (" << COUT_VEC(iter->first) << "," << COUT_VEC(iter->second) << ")" << endl;
              //}
              //MPI_Pause("SHOULD BE ORPHAN!");

              while (faMultimap.begin() != iter_end) {
                // start at repeated key if it exists...
                multimap<const pair<uint8,double*>,double*>::iterator iter_first = faMultimap.begin();
                while (iter_first != iter_end) {
                  multimap<const pair<uint8,double*>,double*>::iterator iter_next = iter_first; ++iter_next;
                  if (iter_next->first == iter_first->first)
                    break;
                  ++iter_first;

                }
                // start at beginning if key DNE...
                if (iter_first == iter_end)
                  iter_first = faMultimap.begin();
                assert(iter_first != faMultimap.end());
                double * x0 = iter_first->first.second;
                noofa_dup_double_ptr_vec.push_back(x0);
                double * x1 = iter_first->second;
                faMultimap.erase(iter_first);
                //cout << "starting at x0: " << x0 << endl;
                while (x1 != x0) {
                  //cout << "next node is: " << iter_first->second << endl;
                  noofa_dup_double_ptr_vec.push_back(x1);
                  iter_first = faMultimap.find(pair<uint8,double*>(this_rbiHash,x1));
                  if (iter_first == faMultimap.end()) {
                    // hop the bridge b/w main and orphan
                    multimap<double*,double*>::iterator iter_x0 = x1x0Multimap.find(x1);
                    assert(iter_x0 != x1x0Multimap.end());
                    x1 = iter_x0->second;
                    x1x0Multimap.erase(iter_x0);
                  }
                  else {
                    x1 = iter_first->second;
                    faMultimap.erase(iter_first);
                  }
                }
                noofa_dup_double_ptr_vec.push_back(NULL);
              }
              assert(x1x0Multimap.empty());

              //cout << "done" << endl;
            }

            noofa_dup_i_vec.push_back(noofa_dup_double_ptr_vec.size());

          } // faMultimap face loop

          // in either case (orphan or not) we should be empty

          assert(faMultimap.empty());

          faocv_dup_i[icv+1] = noofa_dup_i_vec.size()-1;

          // also, make sure we are through the geomtry...

          assert(iter_fgd == internalFaceGeometryMap.end());
          internalFaceGeometryMap.clear();

        }

        if (debug) cout << "10" << endl;

      } // icv loop

      MPI_Sync("end of first icv loop");

      delete[] x_seed;
      delete[] seed_st;
      delete[] stost;

      int8 tot_tri[2];
      MPI_Reduce(my_tot_tri,tot_tri,2,MPI_INT8,MPI_SUM,0,mpi_comm);
      if (mpi_rank == 0) {
        cout << " > surface triangle reduction: unreduced " << tot_tri[0] << " reduced " << tot_tri[1] << endl;
      }

      double max_err[4];
      MPI_Reduce(my_max_err,max_err,4,MPI_DOUBLE,MPI_MAX,0,mpi_comm);
      if (mpi_rank == 0) {
        cout << " > max geometry check errors: vol_cv " << max_err[0] << " x_cv: " << COUT_VEC(max_err+1) << endl;
      }

      bf_zone_of_color_final.clear();
      vector<int>().swap(bf_zone_of_color_final);

      double buf[4];
      MPI_Reduce(my_buf,buf,4,MPI_DOUBLE,MPI_MAX,0,mpi_comm);
      if (mpi_rank == 0) {
        cout << " > boundary max dist face loop (should be zero): " << sqrt(buf[0]) << " max ratio (should be small): " << sqrt(buf[1]) << endl;
        cout << " > internal max dist face loop (should be zero): " << sqrt(buf[2]) << " max ratio (should be small): " << sqrt(buf[3]) << endl;
      }

      // now do the hard work of connecting common face pairs...

      const int nbf_local = bf_data_vec.size();
      assert(nbf_local == bfocv_i[ncv]);
      assert(bf_data_vec.size()+1 == noobf_i_vec.size());
      assert(spobf_i_vec.size() == nbf_local+1);
      assert(sbobf_i_vec.size() == nbf_local+1);

      int nno_local_b = 0;

      // nodeVec links x[3] pointer from vdData or ocdVec or recv_int8 bufs TO a
      // (ino_local,ip) pair.
      vector<pair<double*,pair<int,int> > > nodeVec;

      // =======================================================================================
      // =======================================================================================
      // =======================================================================================
      // At this point we have the noofa_dup_double_ptr_vec/faocv_dup_i and the
      // noobf_double_ptr_vec/bfocv_i for the entire grid built from the vdReturnData using
      // the bfMultimap and faMultimap routines above. We should not have to revisit these
      // individual vds in what follows.
      //
      // Node reduction is going to be based on the nodeVec and use the "lowest-index-array" algo
      // to reduce the nodes, with boundary nodes first. At this point, there has been no consideration
      // of parallelism. We have built the data for the vd's we own.
      // =======================================================================================
      // =======================================================================================
      // =======================================================================================

      int8 * cvora = NULL;
      buildXora(cvora,ncv);
      assert(int8(ncv) == cvora[mpi_rank+1]-cvora[mpi_rank]);

      // average face normals and coordinates so they are the same for
      // both members of the fa_dup associated with the same face...

      averageFaDupData(fa_dup_data_vec,faocv_dup_i);

      // To get the boundary nodes listed first, populate the start of nodeVec with boundary nodes.
      // We need these to appear first in the result/restart file to make certain viz algos efficient.

      vector<int> noobf_v_vec(noobf_double_ptr_vec.size());
      int * cvobf_local = new int[nbf_local];
      for (int icv = 0; icv < ncv; ++icv) {
        for (int ibf = bfocv_i[icv]; ibf != bfocv_i[icv+1]; ++ibf) {
          cvobf_local[ibf] = icv;
          // includes NULL for loop separation
          for (int nob = noobf_i_vec[ibf]; nob != noobf_i_vec[ibf+1]; ++nob) {
            if (noobf_double_ptr_vec[nob]) {
              nodeVec.push_back(pair<double*,pair<int,int> >(noobf_double_ptr_vec[nob],pair<int,int>(nno_local_b,icv)));
              noobf_v_vec[nob] = nno_local_b++;
            }
            else {
              noobf_v_vec[nob] = -1;
            }
          }
        }
      }
      noobf_double_ptr_vec.clear();
      vector<double*>().swap(noobf_double_ptr_vec);

      // at this point noobf_i_vec and noobf_v_vec are consistent

      int nno_local = nno_local_b;

      MPI_Sync("end of second icv loop");

      // before we can push in the internal face nodes, we need to pair them
      // together into a single represetation...

      assert(fa_dup_data_vec.size()+1 == noofa_dup_i_vec.size());
      int my_max_nloop = 0;
      double my_max_dn2 = 0.0;

      // to handle inter-processor (and periodic eventually) faces, we need to exchange
      // their noofa_dup_double_ptr_vec info as well. We only build each face once, even through
      // it is represented twice in the current structures. Pack and send the noofa_dup_double_ptr_vec
      // from the higher rank to the lower rank, and build the face (i.e. loopVec)on the lower rank.

      int * send_int_count = new int[mpi_size];
      int * send_double_count = new int[mpi_size];
      int * send_int8_count = new int[mpi_size];
      FOR_RANK {
        send_int_count[rank] = send_double_count[rank] = send_int8_count[rank] = 0;
      }

      // we know the count of local faces...

      int nfa_local = nfa_i0 + nfa_i1;
      vector<int> noofa_i_vec(nfa_local+1); noofa_i_vec[0] = 0;
      vector<int> noofa_v_vec; // hard to estimate -- NOT noofa_dup_double_ptr_vec.size(), because this includes dup faces
      vector<pair<double*,double*> > loopVec;
      int *cv0ofa_local = new int[nfa_local];
      int *cv1ofa_local_i0 = new int[nfa_i0];
      int *fa_dup_of_fa_local = new int[nfa_local];

      int ifa_local = 0;
      for (int icv = 0; icv < ncv; ++icv) {

        // here duplicate faces include all faces in the current mesh -- i.e. not paired yet...

        for (int ifa_dup = faocv_dup_i[icv]; ifa_dup != faocv_dup_i[icv+1]; ++ifa_dup) {

          int rank,bits,icv_nbr;
          BitUtils::unpackRankBitsIndex(rank,bits,icv_nbr,fa_dup_data_vec[ifa_dup].rbi_nbr);
          if ((rank == mpi_rank)&&(bits == 0)) {

            // this is a face with a local nbr...
            // only build it when icv_nbr > icv...

            if (icv_nbr > icv) {

              // find the relevant face in icv_nbr...

              const uint8 rbi_nbr_nbr = BitUtils::packRankBitsIndex(mpi_rank,0,icv);
              int ifa_dup_nbr;
              for (ifa_dup_nbr = faocv_dup_i[icv_nbr]; ifa_dup_nbr != faocv_dup_i[icv_nbr+1]; ++ifa_dup_nbr) {
                if (fa_dup_data_vec[ifa_dup_nbr].rbi_nbr == rbi_nbr_nbr)
                  break;
              }
              assert(ifa_dup_nbr != faocv_dup_i[icv_nbr+1]); // make sure we found ifa_dup_nbr

              const double dx[3] = DIFF(x_vd[icv_nbr],x_vd[icv]);
              const double delta = 0.5*(delta_vd[icv]+delta_vd[icv_nbr]);

              assert(loopVec.empty());
              //cout << icv << " " << icv_nbr << endl;
              const double my_dn2 = reconcileNodeLoops(loopVec,
                  noofa_dup_double_ptr_vec,noofa_dup_i_vec[ifa_dup],noofa_dup_i_vec[ifa_dup+1]-1,
                  noofa_dup_double_ptr_vec,noofa_dup_i_vec[ifa_dup_nbr],noofa_dup_i_vec[ifa_dup_nbr+1]-1,
                  dx,delta,d2_tol);
              //if (mpi_rank == 331) {
              //  if (my_dn2 > 1.0E-16) {
              //    cout << my_dn2 << " " << icv << " " << icv_nbr << endl;
              //    writeFullCvTecplot(icv,x_vd[icv]);
              //    writeFullCvTecplot(icv_nbr,x_vd[icv_nbr]);
              //    assert(0);
              //  }
              //}
              my_max_dn2 = max(my_max_dn2,my_dn2);

              size_t noofa_v_vec_size0 = noofa_v_vec.size();
              noofa_v_vec.resize(noofa_v_vec_size0+loopVec.size());
              int my_nloop = 0;
              for (int ii = 0, limit = loopVec.size(); ii < limit; ++ii) {
                if (loopVec[ii].first) {
                  assert(loopVec[ii].second);
                  nodeVec.push_back(pair<double*,pair<int,int> >(loopVec[ii].first,pair<int,int>(nno_local,icv)));
                  nodeVec.push_back(pair<double*,pair<int,int> >(loopVec[ii].second,pair<int,int>(nno_local,icv_nbr)));
                  noofa_v_vec[noofa_v_vec_size0+ii] = nno_local++;
                }
                else {
                  assert(loopVec[ii].second == NULL);
                  ++my_nloop;
                  noofa_v_vec[noofa_v_vec_size0+ii] = -1;
                }
              }
              loopVec.clear();

              my_max_nloop = max(my_max_nloop,my_nloop);

              // and add the face...

              cv0ofa_local[ifa_local] = icv;
              fa_dup_of_fa_local[ifa_local] = ifa_dup;
              cv1ofa_local_i0[ifa_local] = icv_nbr;
              ++ifa_local;
              noofa_i_vec[ifa_local] = noofa_v_vec.size();

            }

          }
          else if ((bits < BitUtils::flipPeriodicBits(bits)) || ((bits ==0)&&(mpi_rank > rank))) {

            // we are going to build the loopVec on the winner, so the loser needs to send their nodes over...

            send_int_count[rank] += 4; // icv-(rank-bits)-index and nnof on rank
            send_double_count[rank] += 4 + (noofa_dup_i_vec[ifa_dup+1]-noofa_dup_i_vec[ifa_dup])*3; // xp[3], delta, and 3 doubles per node
            // for the ino globals, skip the nulls...
            for (int nof = noofa_dup_i_vec[ifa_dup]; nof != noofa_dup_i_vec[ifa_dup+1]; ++nof) {
              if (noofa_dup_double_ptr_vec[nof])
                ++send_int8_count[rank]; // ino_local
            }

          }

        }

      }
      assert(ifa_local == nfa_i0);

      MPI_Sync("end of third icv loop");

      int * send_int_disp = new int[mpi_size];
      send_int_disp[0] = 0;
      for (int rank = 1; rank < mpi_size; ++rank)
        send_int_disp[rank] = send_int_count[rank-1] + send_int_disp[rank-1];
      const int send_int_count_sum = send_int_disp[mpi_size-1] + send_int_count[mpi_size-1];
      int * send_int_buf = new int[send_int_count_sum];

      int * send_double_disp = new int[mpi_size];
      send_double_disp[0] = 0;
      for (int rank = 1; rank < mpi_size; ++rank) {
        send_double_disp[rank] = send_double_count[rank-1] + send_double_disp[rank-1];
      }
      const int send_double_count_sum = send_double_disp[mpi_size-1] + send_double_count[mpi_size-1];
      double * send_double_buf = new double[send_double_count_sum];

      // pack the noofa_dup_double_ptr_vec and its index info (i.e. which face it connects to on the recv side)
      // on the loser rank (i.e. higher rank) of any inter-processor face pair...

      for (int icv = 0; icv < ncv; ++icv) {

        for (int ifa_dup = faocv_dup_i[icv]; ifa_dup != faocv_dup_i[icv+1]; ++ifa_dup) {

          int rank,bits,icv_nbr;
          BitUtils::unpackRankBitsIndex(rank,bits,icv_nbr,fa_dup_data_vec[ifa_dup].rbi_nbr);
          const int inv_bits = BitUtils::flipPeriodicBits(bits);
          if ((bits < inv_bits) || ((bits ==0)&&(mpi_rank > rank))) {

            // we are going to build the loopVec on the lower rank, so the higher rank needs to
            // send their nodes over...

            send_int_buf[send_int_disp[rank]  ] = icv_nbr; // icv on the recv side
            send_int_buf[send_int_disp[rank]+1] = BitUtils::packRankBits(mpi_rank,inv_bits);
            send_int_buf[send_int_disp[rank]+2] = icv;
            send_int_buf[send_int_disp[rank]+3] = noofa_dup_i_vec[ifa_dup+1]-noofa_dup_i_vec[ifa_dup];
            send_int_disp[rank] += 4;

            // exchanged doubles start with the xp and delta of icv...

            send_double_buf[send_double_disp[rank]  ] = x_vd[icv][0];
            send_double_buf[send_double_disp[rank]+1] = x_vd[icv][1];
            send_double_buf[send_double_disp[rank]+2] = x_vd[icv][2];

            // periodic transform of above - translate
            if (bits) PeriodicData::periodicTranslate(send_double_buf+send_double_disp[rank],1,inv_bits);

            send_double_buf[send_double_disp[rank]+3] = delta_vd[icv];
            send_double_disp[rank] += 4;

            for (int nof = noofa_dup_i_vec[ifa_dup]; nof != noofa_dup_i_vec[ifa_dup+1]; ++nof) {
              if (noofa_dup_double_ptr_vec[nof]) {
                FOR_I3 send_double_buf[send_double_disp[rank]+i] = noofa_dup_double_ptr_vec[nof][i];
                // periodic transform - rotate
                if (bits) PeriodicData::periodicRotate(send_double_buf+send_double_disp[rank],1,inv_bits);
                send_double_disp[rank] += 3;
              }
              else {
                // for a NULL, we need a unique way to indicate. Use HUGE_VAL...
                FOR_I3 send_double_buf[send_double_disp[rank]++] = HUGE_VAL;
              }
            }

          }

        }

      }

      MPI_Sync("end of 4th icv loop");

      // rewind disp...

      send_int_disp[0] = 0;
      for (int rank = 1; rank < mpi_size; ++rank)
        send_int_disp[rank] = send_int_count[rank-1] + send_int_disp[rank-1];

      int * recv_int_count = new int[mpi_size];
      MPI_Alltoall(send_int_count,1,MPI_INT,recv_int_count,1,MPI_INT,mpi_comm);

      int * recv_int_disp = new int[mpi_size];
      recv_int_disp[0] = 0;
      for (int rank = 1; rank < mpi_size; ++rank)
        recv_int_disp[rank] = recv_int_count[rank-1] + recv_int_disp[rank-1];
      const int recv_int_count_sum = recv_int_disp[mpi_size-1] + recv_int_count[mpi_size-1];

      int * recv_int_buf = new int[recv_int_count_sum];
      MPI_Alltoallv(send_int_buf,send_int_count,send_int_disp,MPI_INT,
          recv_int_buf,recv_int_count,recv_int_disp,MPI_INT,
          mpi_comm);

      delete[] send_int_buf; send_int_buf = NULL;

      delete[] send_int_count;
      delete[] send_int_disp;
      delete[] recv_int_count;
      delete[] recv_int_disp;

      // and doubles...

      send_double_disp[0] = 0;
      for (int rank = 1; rank < mpi_size; ++rank)
        send_double_disp[rank] = send_double_count[rank-1] + send_double_disp[rank-1];

      int * recv_double_count = new int[mpi_size];
      MPI_Alltoall(send_double_count,1,MPI_INT,recv_double_count,1,MPI_INT,mpi_comm);

      int * recv_double_disp = new int[mpi_size];
      recv_double_disp[0] = 0;
      for (int rank = 1; rank < mpi_size; ++rank)
        recv_double_disp[rank] = recv_double_count[rank-1] + recv_double_disp[rank-1];
      const int recv_double_count_sum = recv_double_disp[mpi_size-1] + recv_double_count[mpi_size-1];

      double * recv_double_buf = new double[recv_double_count_sum];
      MPI_Alltoallv(send_double_buf,send_double_count,send_double_disp,MPI_DOUBLE,
          recv_double_buf,recv_double_count,recv_double_disp,MPI_DOUBLE,
          mpi_comm);

      delete[] send_double_buf; send_double_buf = NULL; // useda again below

      delete[] send_double_count;
      delete[] send_double_disp;
      delete[] recv_double_count;
      delete[] recv_double_disp;

      // and int8's, except no exchange...

      int * send_int8_disp = new int[mpi_size];
      send_int8_disp[0] = 0;
      for (int rank = 1; rank < mpi_size; ++rank) {
        send_int8_disp[rank] = send_int8_count[rank-1] + send_int8_disp[rank-1];
      }
      const int send_int8_count_sum = send_int8_disp[mpi_size-1] + send_int8_count[mpi_size-1];

      int * recv_int8_count = new int[mpi_size];
      MPI_Alltoall(send_int8_count,1,MPI_INT,recv_int8_count,1,MPI_INT,mpi_comm);

      int * recv_int8_disp = new int[mpi_size];
      recv_int8_disp[0] = 0;
      for (int rank = 1; rank < mpi_size; ++rank)
        recv_int8_disp[rank] = recv_int8_count[rank-1] + recv_int8_disp[rank-1];
      const int recv_int8_count_sum = recv_int8_disp[mpi_size-1] + recv_int8_count[mpi_size-1];

      int * recv_int8_ino_local = new int[recv_int8_count_sum];
      for (int irecv_int8 = 0; irecv_int8 < recv_int8_count_sum; ++irecv_int8)
        recv_int8_ino_local[irecv_int8] = -1; // not used

      // on the recv side, build the last loops associated with these faces...

      vector<pair<int,int> > sameInoLocalVec;

      uint8 *rbi1ofa_local_i1 = new uint8[nfa_i1];
      int irecv_int = 0;
      int irecv_double = 0;
      int irecv_int8 = 0;
      vector<double*> local_noofa_double_ptr_vec;
      vector<int> irecv_int8_vec;
      while (irecv_int < recv_int_count_sum) {
        const int icv       = recv_int_buf[irecv_int++]; assert((icv >= 0)&&(icv < ncv));
        int rank_nbr,bits_nbr;
        BitUtils::unpackRankBits(rank_nbr,bits_nbr,recv_int_buf[irecv_int++]);
        const int icv_nbr   = recv_int_buf[irecv_int++];
        const int nnof     = recv_int_buf[irecv_int++];

        // recall first 4 doubles are xp[icv_nbr][3] and delta[icv_nbr]...

        const double dx[3] = DIFF(recv_double_buf+irecv_double,x_vd[icv]);
        const double delta = 0.5*(delta_vd[icv]+recv_double_buf[irecv_double+3]);
        irecv_double += 4;

        local_noofa_double_ptr_vec.resize(nnof);
        irecv_int8_vec.resize(nnof);
        for (int nof = 0; nof < nnof; ++nof) {
          if (recv_double_buf[irecv_double+nof*3] == HUGE_VAL) {
            assert(recv_double_buf[irecv_double+nof*3+1] == HUGE_VAL);
            assert(recv_double_buf[irecv_double+nof*3+2] == HUGE_VAL);
            local_noofa_double_ptr_vec[nof] = NULL;
            irecv_int8_vec[nof] = -1;
          }
          else {
            local_noofa_double_ptr_vec[nof] = recv_double_buf+irecv_double+nof*3;
            irecv_int8_vec[nof] = irecv_int8++;
          }
        }
        irecv_double += nnof*3;

        // find the ifa_dup in icv's face list...

        const uint8 rbi_nbr = BitUtils::packRankBitsIndex(rank_nbr,bits_nbr,icv_nbr);
        int ifa_dup;
        for (ifa_dup = faocv_dup_i[icv]; ifa_dup != faocv_dup_i[icv+1]; ++ifa_dup) {
          if (fa_dup_data_vec[ifa_dup].rbi_nbr == rbi_nbr)
            break;
        }
        assert(ifa_dup != faocv_dup_i[icv+1]); // make sure we found ifa_dup

        assert(loopVec.empty());
        //cout << "2: " << icv << " " << icv_nbr << endl;
        const double my_dn2 = reconcileNodeLoops(loopVec,
            noofa_dup_double_ptr_vec,noofa_dup_i_vec[ifa_dup],noofa_dup_i_vec[ifa_dup+1]-1,
            local_noofa_double_ptr_vec,0,nnof-1,
            dx,delta,d2_tol);
        my_max_dn2 = max(my_max_dn2,my_dn2);

        size_t noofa_v_vec_size0 = noofa_v_vec.size();
        noofa_v_vec.resize(noofa_v_vec_size0+loopVec.size());

        int my_nloop = 0;
        for (int ii = 0, limit = loopVec.size(); ii < limit; ++ii) {
          if (loopVec[ii].first) {
            assert(loopVec[ii].second);
            nodeVec.push_back(pair<double*,pair<int,int> >(loopVec[ii].first,pair<int,int>(nno_local,icv)));
            //nodeVec.push_back(pair<double*,pair<int,int> >(loopVec[ii].second,pair<int,int>(nno_local,-1))); // this pairing handled with sameInoLocalVec below
            // get the associated location in the recv_int8_ino_local from the
            // lookup table using some simple pointer arithmetic...
            assert((loopVec[ii].second-local_noofa_double_ptr_vec.front())%3 == 0);
            const int nof = ((loopVec[ii].second-local_noofa_double_ptr_vec.front()))/3;
            assert((nof >= 0)&&(nof < nnof));
            const int this_irecv_int8 = irecv_int8_vec[nof];
            // for a valid node, this should NOT be -1, but instead a location in the global recv_int8_ino_local...
            assert(this_irecv_int8 >= 0);
            if (recv_int8_ino_local[this_irecv_int8] == -1) {
              recv_int8_ino_local[this_irecv_int8] = nno_local;
            }
            else {
              // recv_int8_ino_local[this_irecv_int8] already contains a global index, and it
              // should already be linked to nno_local. check this for now...
              sameInoLocalVec.push_back(pair<int,int>(recv_int8_ino_local[this_irecv_int8],nno_local));
            }
            noofa_v_vec[noofa_v_vec_size0+ii] = nno_local++;
          }
          else {
            assert(loopVec[ii].second == NULL);
            ++my_nloop;
            noofa_v_vec[noofa_v_vec_size0+ii] = -1;
          }
        }
        loopVec.clear();

        my_max_nloop = max(my_max_nloop,my_nloop);

        // add a new ifa_local...

        cv0ofa_local[ifa_local] = icv;
        fa_dup_of_fa_local[ifa_local] = ifa_dup;
        rbi1ofa_local_i1[ifa_local-nfa_i0] = rbi_nbr;
        ++ifa_local;
        noofa_i_vec[ifa_local] = noofa_v_vec.size();

      }
      assert(irecv_int == recv_int_count_sum);
      assert(irecv_double == recv_double_count_sum);
      assert(irecv_int8 == recv_int8_count_sum);
      assert(ifa_local == nfa_local);

      MPI_Sync("5");

      irecv_int8_vec.clear();
      vector<int>().swap(irecv_int8_vec);
      local_noofa_double_ptr_vec.clear();
      vector<double*>().swap(local_noofa_double_ptr_vec);
      delete[] recv_int_buf; recv_int_buf = NULL; // used again below
      delete[] recv_double_buf; recv_double_buf = NULL; // used again below

      // at this point, the recv_int8_ino_local will have a ino_local in it, with potentially
      // some -1's where passed nodes were not matched in the loopVec. Send these back to
      // the send side, where the indexing will be meaningless, BUT the -1's indicate which
      // nodes did not match and shouldn't be indexed...

      // sorry about the "int8" stuff here. These are ints...

      int * send_int8_ino_local = new int[send_int8_count_sum];
      MPI_Alltoallv(recv_int8_ino_local,recv_int8_count,recv_int8_disp,MPI_INT,
          send_int8_ino_local,send_int8_count,send_int8_disp,MPI_INT,
          mpi_comm);

      // record the face bits associated with the node pairings...

      int * send_int8_faono_bits = new int[send_int8_count_sum];

      // pass through the loser side and index loopVec nodes...

      for (int icv = 0; icv < ncv; ++icv) {

        for (int ifa_dup = faocv_dup_i[icv]; ifa_dup != faocv_dup_i[icv+1]; ++ifa_dup) {

          int rank,bits,icv_nbr;
          BitUtils::unpackRankBitsIndex(rank,bits,icv_nbr,fa_dup_data_vec[ifa_dup].rbi_nbr);
          if ((bits < BitUtils::flipPeriodicBits(bits))||((bits == 0)&&(mpi_rank > rank))) { // on-proc non-periodic losers handled above

            for (int nof = noofa_dup_i_vec[ifa_dup]; nof != noofa_dup_i_vec[ifa_dup+1]; ++nof) {
              if (noofa_dup_double_ptr_vec[nof]) {
                if (send_int8_ino_local[send_int8_disp[rank]] >= 0) {
                  nodeVec.push_back(pair<double*,pair<int,int> >(noofa_dup_double_ptr_vec[nof],pair<int,int>(nno_local,icv)));
                  send_int8_ino_local[send_int8_disp[rank]] = nno_local++;
                  send_int8_faono_bits[send_int8_disp[rank]] = bits;
                }
                else {
                  assert(send_int8_ino_local[send_int8_disp[rank]] == -1);
                  send_int8_faono_bits[send_int8_disp[rank]] = 0;
                }
                ++send_int8_disp[rank];
              }
            }

          }

        }

      }

      MPI_Sync("6");

      //rbi_nbr_of_fa_dup_vec.clear();
      noofa_dup_double_ptr_vec.clear();
      vector<double*>().swap(noofa_dup_double_ptr_vec);
      noofa_dup_i_vec.clear();
      vector<int>().swap(noofa_dup_i_vec);
      delete[] faocv_dup_i;

      // rewind -- it gets used below...

      send_int8_disp[0] = 0;
      for (int rank = 1; rank < mpi_size; ++rank)
        send_int8_disp[rank] = send_int8_count[rank-1] + send_int8_disp[rank-1];

      int max_nloop;
      MPI_Reduce(&my_max_nloop,&max_nloop,1,MPI_INT,MPI_MAX,0,mpi_comm);
      double max_dn2;
      MPI_Reduce(&my_max_dn2,&max_dn2,1,MPI_DOUBLE,MPI_MAX,0,mpi_comm);
      if (mpi_rank == 0)
        cout << " > max number of loops per face: " << max_nloop << " max area-weighted face normal difference over delta^2 for reconciled face (should be small): " << sqrt(max_dn2) << endl;

      // sort by double* so that global indices associated with the same double* can be set
      // the same...

      sort(nodeVec.begin(),nodeVec.end());

      assert(nno_local < TWO_BILLION);
      int * no_flag = new int[nno_local];
      for (int ino = 0; ino < nno_local; ++ino) no_flag[ino] = ino;

      {
        int iend = 0;
        while (iend < nodeVec.size()) {
          const int istart = iend++;
          double* current = nodeVec[istart].first;
          while ((iend < nodeVec.size())&&(current == nodeVec[iend].first)) ++iend;
          if (iend > istart+1) {
            //cout << "found some common doubles: " << iend-i << endl;
            int ino_smallest = nno_local; // can't be this large
            for (int ii = istart; ii < iend; ++ii) {
              int ino = no_flag[nodeVec[ii].second.first];
              while (ino != no_flag[ino]) ino = no_flag[ino];
              ino_smallest = min(ino_smallest,ino);
              //cout << "ii: " << ii << " ino: " << ino << " ino_smallest: " << ino_smallest << endl;
            }
            assert(ino_smallest == no_flag[ino_smallest]);
            for (int ii = istart; ii < iend; ++ii) {
              int ino = no_flag[nodeVec[ii].second.first];
              while (ino != no_flag[ino]) ino = no_flag[ino];
              no_flag[ino] = ino_smallest;
            }
          }
        }
      }

      // include sameInoLocalVec...
      for (int ii = 0,limit = sameInoLocalVec.size(); ii < limit; ++ii) {
        int ino0 = sameInoLocalVec[ii].first;
        while (ino0 != no_flag[ino0]) ino0 = no_flag[ino0];
        int ino1 = sameInoLocalVec[ii].second;
        while (ino1 != no_flag[ino1]) ino1 = no_flag[ino1];
        no_flag[ino0] = no_flag[ino1] = min(ino0,ino1);
      }
      sameInoLocalVec.clear();
      vector<pair<int,int> >().swap(sameInoLocalVec);

      // stable_sort based on icv to get nodes in the same vd. We use
      // stable sort here because this sorted order is also used below
      // to reduce node locations, which are sensitive to the order of
      // operations (at machine precision).
      //stable_sort(nodeVec.begin(),nodeVec.end(),secondPairSecondInt);
      // actually, the memory can be different in each run
      sort(nodeVec.begin(),nodeVec.end(),secondPairSecondInt);

      {
        int iend = 0;
        int prev_icv = -1;
        while (iend < nodeVec.size()) {
          const int istart = iend++;
          const int current_icv = nodeVec[istart].second.second;
          assert((prev_icv+1) == current_icv);
          prev_icv = current_icv;
          while ((iend < nodeVec.size())&&(current_icv == nodeVec[iend].second.second)) ++iend;
          // we have our range, now compare dists to see if we can reduce more with N^2 algo
          for (int ii = istart; ii < iend-1; ++ii) {
            for (int jj = ii+1; jj < iend; ++jj) {
              if (nodeVec[ii].first != nodeVec[jj].first) {
                if (DIST2(nodeVec[ii].first,nodeVec[jj].first)/(delta_vd[current_icv]*delta_vd[current_icv]) < d2_tol) {
                  int ino0 = no_flag[nodeVec[ii].second.first];
                  while (ino0 != no_flag[ino0]) ino0 = no_flag[ino0];
                  int ino1 = no_flag[nodeVec[jj].second.first];
                  while (ino1 != no_flag[ino1]) ino1 = no_flag[ino1];
                  no_flag[ino0] = no_flag[ino1] = min(ino0,ino1) ;
                }
              }
            }
          }
        }
      }

      int nno_local_reduced = 0;
      int nno_local_b_reduced = 0;
      for (int ino = 0; ino < nno_local; ++ino) {
        if (no_flag[ino] == ino) {
          ++nno_local_reduced;
          if (ino < nno_local_b)
            nno_local_b_reduced = nno_local_reduced;
          no_flag[ino] = -nno_local_reduced;
        }
        else {
          int ino_ = no_flag[ino];
          while (ino_ >= 0) ino_ = no_flag[ino_];
          no_flag[ino] = ino_;
        }
      }

      // At this point, all nodes are -1 indexed in no_flag.
      // before we throw out no_flag, modify any references to nodes...

      for (int nob = 0, limit = noobf_v_vec.size(); nob < limit; ++nob) {
        if (noobf_v_vec[nob] >= 0) {
          noobf_v_vec[nob] = -no_flag[noobf_v_vec[nob]]-1; // note no_flag currently -1 indexed
        }
        else {
          assert(noobf_v_vec[nob] == -1);
        }
      }

      for (int nof = 0, limit = noofa_v_vec.size(); nof < limit; ++nof) {
        if (noofa_v_vec[nof] >= 0) {
          noofa_v_vec[nof] = -no_flag[noofa_v_vec[nof]]-1; // note no_flag currently -1 indexed
        }
        else {
          assert(noofa_v_vec[nof] == -1);
        }
      }

      for (int isend = 0; isend < send_int8_count_sum; ++isend) {
        const int ino = send_int8_ino_local[isend];
        if (ino == -1) continue;
        assert((ino >= 0)&&(ino < nno_local));
        // these guys don't know their loops
        send_int8_ino_local[isend] = -no_flag[ino]-1; // note no_flag currently -1 indexed
      }

      for (int irecv = 0; irecv < recv_int8_count_sum; ++irecv) {
        const int ino = recv_int8_ino_local[irecv];
        if (ino == -1) continue;
        assert((ino >= 0)&&(ino < nno_local));
        recv_int8_ino_local[irecv] = -no_flag[ino]-1;
      }

      for (int ii = 0, limit = nodeVec.size(); ii < limit; ++ii) {
        // the first of the pair is "ino_local"
        nodeVec[ii].second.first = -no_flag[nodeVec[ii].second.first]-1;
        assert((nodeVec[ii].second.first >= 0)&&(nodeVec[ii].second.first < nno_local_reduced));
      }

      delete[] no_flag;

      {

        int8 vd_mem = 0;
        int8 vd_node_count = 0;
        for (int ii = 0; ii < ncv; ++ii) {
          vd_mem += vdReturnDataVec[ii].getByteCount();
          vd_node_count += vdReturnDataVec[ii].nno;
        }

        int8 my_count[5] = { int8(nodeVec.size()*sizeof(pair<double*,pair<int,int> >)),
          int8(noobf_v_vec.size()*sizeof(int)),
          int8(noofa_v_vec.size()*sizeof(int)),
          vd_mem,
          vd_node_count };
        int8 count[5];
        MPI_Reduce(my_count,count,5,MPI_INT8,MPI_MAX,0,mpi_comm);
        if (mpi_rank == 0) cout << " > memory max: nodeVec MB: " << double(count[0])*1.0E-6 <<
          " noobf_v_vec MB: " << double(count[1])*1.0E-6 <<
            " noofa_v_vec MB: " << double(count[2])*1.0E-6 <<
            " vdReturnDataVec MB: " << double(count[3])*1.0E-6 <<
            " vd node count: " << count[4] << endl;

        MPI_Reduce(my_count,count,5,MPI_INT8,MPI_SUM,0,mpi_comm);
        if (mpi_rank == 0) cout << " > memory avg: nodeVec MB: " << double(count[0])/double(mpi_size)*1.0E-6 <<
          " noobf_v_vec MB: " << double(count[1])/double(mpi_size)*1.0E-6 <<
            " noofa_v_vec MB: " << double(count[2])/double(mpi_size)*1.0E-6 <<
            " vdReturnDataVec MB: " << double(count[3])/double(mpi_size)*1.0E-6 <<
            " vd node count: " << double(count[4])/double(mpi_size) << endl;

      }

      // =======================================
      // set the no_bits to contain the bits of any
      // touching face. Boundary faces are (1<<6)...
      // =======================================

      int * no_bits = new int[nno_local_reduced];

      // set boundary. These are first right now...
      // TODO XXXXXXXX : note that this does not consider periodicity at this point. Is this correct?
      for (int ino_local = 0; ino_local < nno_local_b_reduced; ++ino_local)
        no_bits[ino_local] = (1<<6);

      // and the rest start with no bits...
      for (int ino_local = nno_local_b_reduced; ino_local < nno_local_reduced; ++ino_local)
        no_bits[ino_local] = 0;

      // add and exchange periodic bits. Any periodic bits have been stored on loser
      // faces in the array send_int8_faono_bits (built above) and are added now...

      assert(send_int_buf == NULL); send_int_buf = new int[send_int8_count_sum];
      assert(recv_int_buf == NULL); recv_int_buf = new int[recv_int8_count_sum];

#include "compressNoBits.hpp" // not really a "compress", but similar in format to compressNodesNew used below

      delete[] send_int_buf; send_int_buf = NULL;
      delete[] recv_int_buf; recv_int_buf = NULL;
      delete[] send_int8_faono_bits; send_int8_faono_bits = NULL;

      /*
         {

         cout << "nno_local_reduced: " << nno_local_reduced << " nno_local_b_reduced: " << nno_local_b_reduced << endl;

         int iend = 0;
         int prev_icv = -1;
         while (iend < nodeVec.size()) {
         const int istart = iend++;
         const int current_icv = nodeVec[istart].second.second;
         const double * const current_xp = x_vd[current_icv];
         assert((prev_icv+1) == current_icv);
         prev_icv = current_icv;
         while ((iend < nodeVec.size())&&(current_icv == nodeVec[iend].second.second)) ++iend;
         for (int ii = istart; ii < iend; ++ii) {
         const int ino_local = nodeVec[ii].second.first;
         cout << "got a node: ino_local: " << ino_local << " this_x_i: ";
         FOR_I3 {
         const double this_x_i = nodeVec[ii].first[i] + current_xp[i];
         cout << this_x_i << " ";
         }
         cout << " bits: ";
         for (int i = 0; i <= 6; ++i) {
         if (no_bits[ino_local]&(1<<i)) cout << "1";
         else cout << "0";
         }
         cout << endl;
         }
         }

         }
         */

      // ==================================================
      // at this point the nodes have their actual bits, not
      // neccessarily only even bits. So that each node lives in
      // the appropriate group, we accept any sets of bits. In
      // the future, when we consider axisymmetric singularities,
      // I think we can just add a group to getNodeGroup above.
      // ==================================================

      // index all nodes as follows:
      // all even periodic in the first 7 places (0-6)
      // then actual boundary nodes   7
      // then other (internal)        8

      // the maximum group size. We need to make room for up to this
      // many nodes in each group. It turns out that if one node is
      // present in a group, it is not guaranteed that all nodes will
      // be present (larger groups only). There seems to be 4 of
      // 8 present in certain triply periodic grids. This spacing is
      // corrected in the node ordering at the end.

      const int nng = 9;
      //const int ing_size[9] = { 2, 2, 4, 2, 4, 4, 8, 1, 1 }; // old node grouping
      const int ing_size[9] = { 8, 4, 4, 4, 2, 2, 2, 1, 1 }; // regrouped like this to allow 8 alignment to align all

      int8 my_count[nng];
      for (int i = 0; i < nng; ++i) {
        assert(ing_size[i] == PeriodicFunc::getNodeGroupSize(i)); // double-check the above array is consistent
        my_count[i] = 0;
      }

      for (int ino_local = 0; ino_local < nno_local_reduced; ++ino_local) {
        const int ing = PeriodicFunc::getNodeGroup(no_bits[ino_local]);
        assert((ing >= 0)&&(ing < nng));
        my_count[ing] += ing_size[ing];
      }

      int8 count[nng];
      MPI_Allreduce(my_count,count,nng,MPI_INT8,MPI_SUM,mpi_comm);

      int8 nno_global_early = 0;
      for (int igr = 0; igr < nng; ++igr)
        nno_global_early += count[igr];

      //if (mpi_rank == 0) {
      // cout << " > node group counts:" << endl;
      // for (int i = 0; i < nng; ++i)
      //   cout << "    > ing: " << i << " count: " << count[i] << endl;
      //}

      int8 my_disp[nng];
      int8 my_disp_copy[nng]; // for checking
      MPI_Scan(my_count,my_disp,nng,MPI_INT8,MPI_SUM,mpi_comm);
      for (int i = 0; i < nng; ++i) {
        my_disp[i] -= my_count[i];
        my_disp_copy[i] = my_disp[i];
      }

      // index ALL nodes

      int8 * ino_global = new int8[nno_local_reduced];

      for (int ino_local = 0; ino_local < nno_local_reduced; ++ino_local) {
        const int ing = PeriodicFunc::getNodeGroup(no_bits[ino_local]);
        assert((ing >= 0)&&(ing < nng));
        ino_global[ino_local] = my_disp[ing]; my_disp[ing] += ing_size[ing]; // note addition of group size to leave room for periodic
        for (int ing_prev = 0; ing_prev < ing; ++ing_prev)
          ino_global[ino_local] += count[ing_prev];
      }

      // rewind my_disp...

      for (int i = 0; i < nng; ++i) {
        my_disp[i] -= my_count[i];
        assert(my_disp_copy[i] == my_disp[i]); // check
      }

      // ======================================================
      // node reduction...
      // ======================================================

      // compressNodesNew requires these 2 buffers...

      int8 * send_int8_buf = new int8[send_int8_count_sum];
      int8 * recv_int8_buf = new int8[recv_int8_count_sum];

#include "compressNodes.hpp"

      // identify and eliminate face loops (and their nodes potentially) that have less than 3 unique points...
      // also use this as a chance to reorder the boundary faces...

      // ===============================================
      // bfs first...
      // ===============================================

      vector<pair<int,int> > zone_local_pair_vec(nbf_local);
      for (int ibf_local = 0; ibf_local < nbf_local; ++ibf_local) {
        zone_local_pair_vec[ibf_local].first = bf_data_vec[ibf_local].zone;
        zone_local_pair_vec[ibf_local].second = ibf_local;
      }
      //zone_of_bf_vec.clear();

      sort(zone_local_pair_vec.begin(),zone_local_pair_vec.end());

      int8 * my_nbf_zone = new int8[nzone];
      int8 * my_noobf_zone = new int8[nzone];
      int8 * my_spobf_zone = new int8[nzone];
      int8 * my_sbobf_zone = new int8[nzone];

      // also the geomtry...
      double (*my_geom_bf_zone)[8] = new double[nzone][8];

      for (int izone = 0; izone < nzone; ++izone) {
        my_nbf_zone[izone] = 0;
        my_noobf_zone[izone] = 0;
        my_spobf_zone[izone] = 0;
        my_sbobf_zone[izone] = 0;
        // zero geom stuff...
        my_geom_bf_zone[izone][0] = 0.0; // sum(area_bf)
        my_geom_bf_zone[izone][1] = 0.0; // sum(area_over_delta_bf)
        my_geom_bf_zone[izone][2] = 0.0; // sum(n_bf[0])
        my_geom_bf_zone[izone][3] = 0.0; // sum(n_bf[1])
        my_geom_bf_zone[izone][4] = 0.0; // sum(n_bf[2])
        my_geom_bf_zone[izone][5] = 0.0; // x_bf[0] avg (eventually)
        my_geom_bf_zone[izone][6] = 0.0; // x_bf[1] avg (eventually)
        my_geom_bf_zone[izone][7] = 0.0; // x_bf[2] avg (eventually)
      }

      int * noobf_i = new int[nbf_local+1];
      noobf_i[0] = 0;
      vector<int> noobf_v; noobf_v.reserve(noobf_v_vec.size());

      // surface-point-of-bf...
      int * spobf_i = new int[nbf_local+1];
      spobf_i[0] = 0;
      vector<uint8> spobf_v; spobf_v.reserve(spobf_v_wgt_vec.size()); // we go to uint8 here to support global ist > 2BILLION in the future
      vector<double> spobf_wgt; spobf_wgt.reserve(spobf_v_wgt_vec.size()); // "

      // surface-tri-of-bf...
      // TODO: we can extract the part index from this. part is a bf-property, like zone. part should be the same for ALL surface tris of a single bf,
      // even if intersections are incolrporated, coloring should not cross the intersections of parts...
      int * sbobf_i = new int[nbf_local+1];
      sbobf_i[0] = 0;
      vector<uint8> psbobf_v; psbobf_v.reserve(psbobf_v_vec.size()); // we go to uint8 here to support global ist > 2BILLION in the future
      int8 * cvobf_global = new int8[nbf_local]; // this is atleast big enough

      double (*n_bf)[3]          = new double[nbf_local][3];
      double (*x_bf)[3]          = new double[nbf_local][3];
      double *area_bf            = new double[nbf_local];
      double *area_over_delta_bf = new double[nbf_local];
      double (*Gij_bf)[3][3]     = new double[nbf_local][3][3];

      int nbf = 0; // final local nbf
      vector<int> nodeLoopVec;

      for (int ii = 0, limit = zone_local_pair_vec.size(); ii < limit; ++ii) {
        const int zone = zone_local_pair_vec[ii].first;
        const int ibf_local = zone_local_pair_vec[ii].second;
        // includes NULL for loop separation, including when there is just one loop.
        const int nob_final_first = noobf_v.size();
        const int spob_final_first = spobf_v.size();
        assert(spob_final_first == spobf_wgt.size());
        const int stob_final_first = psbobf_v.size(); // ipart/ist/bits now
        int nob1 = noobf_i_vec[ibf_local];
        int loop_count_final = 0;
        while (nob1 < noobf_i_vec[ibf_local+1]) {
          assert(nodeLoopVec.empty());
          const int nob0 = nob1++; // can we increment? a loop with nothing is possible?
          assert(noobf_v_vec[nob0] >= 0); // check above
          while (noobf_v_vec[nob1] >= 0) ++nob1;
          assert(noobf_v_vec[nob1] == -1); // should end with -1
          assert(nodeLoopVec.empty());
          for (int nob = nob0; nob < nob1; ++nob) {
            assert(noobf_v_vec[nob] >= 0);
            const int ino_local = noobf_v_vec[nob]; assert(ino_local < nno_local_reduced);
            if ( nodeLoopVec.empty() ||
                (max(ino_global[ino_local],-ino_global[ino_local]-2) !=
                 max(ino_global[nodeLoopVec.back()],-ino_global[nodeLoopVec.back()]-2)) ||
                (no_bits[ino_local] != no_bits[nodeLoopVec.back()]) ) {
              nodeLoopVec.push_back(ino_local);
            }
          }
          // if the last node and first node are the same, pop the back...
          if ( (max(ino_global[nodeLoopVec.front()],-ino_global[nodeLoopVec.front()]-2) ==
                max(ino_global[nodeLoopVec.back()],-ino_global[nodeLoopVec.back()]-2)) &&
              (no_bits[nodeLoopVec.front()] == no_bits[nodeLoopVec.back()]) ) {
            nodeLoopVec.pop_back(); // remove the last element -- it matched the first
          }
          // advance nob1 -- skip trailing -1...
          ++nob1;
          // only add this loop if it is 3 or larger...
          if (nodeLoopVec.size() >= 3) {
            // we got atleast 3 unique elements, so add...
            ++loop_count_final;
            // on the second and later loops, add the very first
            if (loop_count_final > 1)
              noobf_v.push_back(noobf_v[nob_final_first]);
            for (int ii = 0; ii < nodeLoopVec.size(); ++ii) {
              noobf_v.push_back(nodeLoopVec[ii]);
              // and flip touched nodes to -2 indexing...
              if (ino_global[nodeLoopVec[ii]] >= 0)
                ino_global[nodeLoopVec[ii]] = -ino_global[nodeLoopVec[ii]]-2;
            }
            // on the second and later loops, close the loop...
            if (loop_count_final > 1)
              noobf_v.push_back(nodeLoopVec[0]);
          }
          nodeLoopVec.clear();
        }
        assert(nob1 == noobf_i_vec[ibf_local+1]);
        // if we got any new loops, then this bf is valid...
        if (loop_count_final >= 1) {
          cvobf_global[nbf] = cvora[mpi_rank] + cvobf_local[ibf_local];
          my_nbf_zone[zone] += 1;
          my_noobf_zone[zone] += noobf_v.size() - nob_final_first;
          // copy over surface points and weights for this bf..
          for (int spob = spobf_i_vec[ibf_local]; spob != spobf_i_vec[ibf_local+1]; ++spob) {
            // note here we break them into two separate vectors for
            // efficient writing...
            spobf_v.push_back(spobf_v_wgt_vec[spob].first);
            spobf_wgt.push_back(spobf_v_wgt_vec[spob].second);
          }
          my_spobf_zone[zone] += spobf_v.size() - spob_final_first;
          // copy over surface tris for this bf...
          for (int stob = sbobf_i_vec[ibf_local]; stob != sbobf_i_vec[ibf_local+1]; ++stob) {
            const pair<pair<int,int>,int> ipart_ist_bits_pair = psbobf_v_vec[stob];
            assert(ipart_ist_bits_pair.second < (1<<6)); // (1<<12) is actually fine here some day: 52 bit shift anticipates the possibility of double (or triple) transforms on ist's
            //psbobf_v.push_back( (uint8(ipart_ist_bits_pair.first.second)) |
            //                    (uint8(ipart_ist_bits_pair.first.first)<<32) |  // ipart shift 32
            //                    (uint8(ipart_ist_bits_pair.second)<<52) ); // shift any bits 52 to allow 12 bits of periodicity (some day)
            // store the flattened surface index
            psbobf_v.push_back( uint8(getFlattenedSt(ipart_ist_bits_pair.first.first,ipart_ist_bits_pair.first.second)) |
                (uint8(ipart_ist_bits_pair.second)<<52) ); // shift any bits 52 to allow 12 bits of periodicity (some day)
          }
          my_sbobf_zone[zone] += psbobf_v.size() - stob_final_first;
          // geometry...
          FOR_I3 n_bf[nbf][i] = bf_data_vec[ibf_local].n[i];
          FOR_I3 x_bf[nbf][i] = bf_data_vec[ibf_local].x[i] + x_vd[cvobf_local[ibf_local]][i];
          area_bf[nbf] = bf_data_vec[ibf_local].area;
          area_over_delta_bf[nbf] = bf_data_vec[ibf_local].area_over_delta;
          FOR_I3 FOR_J3 Gij_bf[nbf][i][j] = bf_data_vec[ibf_local].Gij[i][j];
          // build the partial geometric sums in the bfZoneVec...
          // 0: sum(area_bf)
          // 1: sum(area_over_delta_bf)
          // 2: sum(n_bf[0])
          // 3: sum(n_bf[1])
          // 4: sum(n_bf[2])
          // 5: x_bf[0] avg (eventually)
          // 6: x_bf[1] avg (eventually)
          // 7: x_bf[2] avg (eventually)
          my_geom_bf_zone[zone][0] += area_bf[nbf];
          my_geom_bf_zone[zone][1] += area_over_delta_bf[nbf];
          my_geom_bf_zone[zone][2] += n_bf[nbf][0];
          my_geom_bf_zone[zone][3] += n_bf[nbf][1];
          my_geom_bf_zone[zone][4] += n_bf[nbf][2];
          my_geom_bf_zone[zone][5] += area_bf[nbf]*x_bf[nbf][0];
          my_geom_bf_zone[zone][6] += area_bf[nbf]*x_bf[nbf][1];
          my_geom_bf_zone[zone][7] += area_bf[nbf]*x_bf[nbf][2];
          // and increment...
          ++nbf;
          noobf_i[nbf] = noobf_v.size();
          spobf_i[nbf] = spobf_v.size();
          assert(spobf_i[nbf] == spobf_wgt.size());
          sbobf_i[nbf] = psbobf_v.size();
        }
      }
      delete[] bfocv_i; // no longer needed
      delete[] cvobf_local;
      bf_data_vec.clear();
      vector<BfData>().swap(bf_data_vec);
      noobf_v_vec.clear();
      vector<int>().swap(noobf_v_vec);
      noobf_i_vec.clear();
      vector<int>().swap(noobf_i_vec);
      psbobf_v_vec.clear();
      vector<pair<pair<int,int>,int> >().swap(psbobf_v_vec);
      sbobf_i_vec.clear();
      vector<int>().swap(sbobf_i_vec);
      spobf_v_wgt_vec.clear();
      vector<pair<uint8,double> >().swap(spobf_v_wgt_vec);
      spobf_i_vec.clear();
      vector<int>().swap(spobf_i_vec);

      assert(nbf <= nbf_local);

      for (int nob = 0; nob < noobf_v.size(); ++nob)
        assert(ino_global[noobf_v[nob]] <= -2);

      // ===============================================
      // faces next...
      // ===============================================

      zone_local_pair_vec.clear(); // we don't know the size exactly because of periodic faces, so...

      // the first faces are all internal, no bits...
      for (int ifa_local = 0; ifa_local < nfa_i0; ++ifa_local) {
        const int zone = PeriodicFunc::getFaceZoneForBits(0); // bits == 0 for these faces
        zone_local_pair_vec.push_back(pair<int,int>(zone,ifa_local));
      }

      // the next faces are interprocessor and have an associated rbi...
      for (int ifa_local = nfa_i0; ifa_local < nfa_local; ++ifa_local) {
        int rank,bits,index;
        BitUtils::unpackRankBitsIndex(rank,bits,index,rbi1ofa_local_i1[ifa_local-nfa_i0]);
        const int zone = PeriodicFunc::getFaceZoneForBits(bits);
        zone_local_pair_vec.push_back(pair<int,int>(zone,ifa_local));
        // if this is a periodic face, push its periodic loser...
        if (bits) {
          const int inv_bits = BitUtils::flipFirstSixPeriodicBits(bits);
          assert(bits > inv_bits);
          const int zone_shadow = PeriodicFunc::getFaceZoneForBits(inv_bits);
          assert(zone > zone_shadow);
          zone_local_pair_vec.push_back(pair<int,int>(zone_shadow,ifa_local));
        }
      }

      // now we have the final face count including periodic...
      const int nfa_local_with_periodic = zone_local_pair_vec.size();

      sort(zone_local_pair_vec.begin(),zone_local_pair_vec.end());

      int8 * my_nfa_zone = new int8[27];
      int8 * my_noofa_zone = new int8[27];
      for (int i = 0; i < 27; ++i) {
        my_nfa_zone[i] = 0;
        my_noofa_zone[i] = 0;
      }

      int * noofa_i = new int[nfa_local_with_periodic+1];
      noofa_i[0] = 0;
      vector<pair<int,int> > noofa_v; noofa_v.reserve(noofa_v_vec.size());

      int8 (*cvofa_global)[2] = new int8[nfa_local_with_periodic][2]; // this is atleast big enough

      // face geometry...

      double (*n_fa)[3] = new double[nfa_local_with_periodic][3];
      double (*x_fa)[3] = new double[nfa_local_with_periodic][3];

      int nfa = 0; // final local nfa
      for (int ii = 0, limit = zone_local_pair_vec.size(); ii < limit; ++ii) {
        const int zone = zone_local_pair_vec[ii].first;
        const int bits = PeriodicFunc::getBitsForFaceZone(zone);
        const int flip_bits = BitUtils::flipFirstSixPeriodicBits(bits);
        const int ifa_local = zone_local_pair_vec[ii].second;
        // includes NULL for loop separation, including when there is just one loop.
        const int nof_final_first = noofa_v.size();
        int nof1 = noofa_i_vec[ifa_local];
        int loop_count_final = 0;
        while (nof1 < noofa_i_vec[ifa_local+1]) {
          assert(nodeLoopVec.empty());
          const int nof0 = nof1++; // can we increment? a loop with nothing is possible?
          assert(noofa_v_vec[nof0] >= 0); // check above
          while (noofa_v_vec[nof1] >= 0) ++nof1;
          assert(noofa_v_vec[nof1] == -1); // should end with -1
          assert(nodeLoopVec.empty());
          for (int nof = nof0; nof < nof1; ++nof) {
            assert(noofa_v_vec[nof] >= 0);
            const int ino_local = noofa_v_vec[nof]; assert(ino_local < nno_local_reduced);
            if ( nodeLoopVec.empty() ||
                (max(ino_global[ino_local],-ino_global[ino_local]-2) !=
                 max(ino_global[nodeLoopVec.back()],-ino_global[nodeLoopVec.back()]-2)) ||
                (no_bits[ino_local] != no_bits[nodeLoopVec.back()]) ) {
              nodeLoopVec.push_back(ino_local);
            }
          }
          // if the last node and first node are the same, pop the back...
          if ( (max(ino_global[nodeLoopVec.front()],-ino_global[nodeLoopVec.front()]-2) ==
                max(ino_global[nodeLoopVec.back()],-ino_global[nodeLoopVec.back()]-2)) &&
              (no_bits[nodeLoopVec.front()] == no_bits[nodeLoopVec.back()]) ) {
            nodeLoopVec.pop_back(); // remove the last element -- it matched the first
          }
          // advance nof1 -- skip trailing -1...
          ++nof1;
          // only add this loop if it is 3 or larger...
          if (nodeLoopVec.size() >= 3) {
            // we got atleast 3 unique elements, so add...
            ++loop_count_final;
            // on the second and later loops, add the very first
            if (loop_count_final > 1) {
              noofa_v.push_back(noofa_v[nof_final_first]);
            }
            const int nof_final_first_loop = noofa_v.size();
            for (int ii = 0; ii < nodeLoopVec.size(); ++ii) {
              const int this_ino_local = nodeLoopVec[ii];
              if (flip_bits > bits) { // TODO: some other way to determine if this face is a loser...
                // this is a loser face, so we need to use the face bits to modify the node bits. We
                // can also check that this does not modify the nodes in a way that results in a
                // double transform. It should not, because the node bits were set by the face bits
                // originally...
                int this_no_bits = no_bits[this_ino_local];
                FOR_I3 {
                  if (bits & (1<<(2*i))) {
                    assert(!(bits & (1<<(2*i+1))));
                    assert(this_no_bits & (1<<(2*i+1)));
                    this_no_bits &= ~(1<<(2*i+1));
                    this_no_bits |= (1<<(2*i));
                  }
                  else if (bits & (1<<(2*i+1))) {
                    assert(!(bits & (1<<(2*i))));
                    assert(this_no_bits & (1<<(2*i)));
                    this_no_bits &= ~(1<<(2*i));
                    this_no_bits |= (1<<(2*i+1));
                  }
                }
                noofa_v.push_back(pair<int,int>(this_ino_local,this_no_bits));
              }
              else {
                noofa_v.push_back(pair<int,int>(this_ino_local,no_bits[this_ino_local]));
              }
              // and flip touched nodes to -2 indexing...
              if (ino_global[this_ino_local] >= 0)
                ino_global[this_ino_local] = -ino_global[this_ino_local]-2;
            }
            // on the second and later loops, close the loop...
            if (loop_count_final > 1) {
              noofa_v.push_back(noofa_v[nof_final_first_loop]);
            }
          }
          nodeLoopVec.clear();
        }
        assert(nof1 == noofa_i_vec[ifa_local+1]);
        // if we got any new loops, then this fa is valid...
        if (loop_count_final >= 1) {
          // if this is a loser face copy, we need to make sure everything is flipped appropriately...
          if (flip_bits > bits) {
            assert(ifa_local >= nfa_i0);
            int rank,bits_check,index;
            BitUtils::unpackRankBitsIndex(rank,bits_check,index,rbi1ofa_local_i1[ifa_local-nfa_i0]);
            assert(bits_check == flip_bits);
            cvofa_global[nfa][0] = (cvora[rank]+index);
            cvofa_global[nfa][1] = (cvora[mpi_rank]+cv0ofa_local[ifa_local]) | (int8(bits)<<55); // TODO: why 55?
            // zone counts...
            my_nfa_zone[zone] += 1;
            // for noofa_v, completely flip the order...
            const int nnof = noofa_v.size() - nof_final_first;
            my_noofa_zone[zone] += nnof;
            for (int i = 0; i < nnof/2; ++i) {
              const int nof0_final = nof_final_first+i;
              const int nof1_final = nof_final_first+nnof-i-1;
              const pair<int,int> tmp = noofa_v[nof0_final];
              noofa_v[nof0_final] = noofa_v[nof1_final];
              noofa_v[nof1_final] = tmp;
            }
            // geometry (already averaged)...
            FOR_I3 n_fa[nfa][i] = -fa_dup_data_vec[fa_dup_of_fa_local[ifa_local]].n[i]; // note sign change
            FOR_I3 x_fa[nfa][i] = fa_dup_data_vec[fa_dup_of_fa_local[ifa_local]].x[i];
            // since this is the loser, we need to transform these points/normals...
            PeriodicData::periodicRotate(n_fa[nfa],1,bits);
            PeriodicData::periodicTranslate(x_fa[nfa],1,bits);
            ++nfa;
            noofa_i[nfa] = noofa_v.size();
          }
          else {
            cvofa_global[nfa][0] = cvora[mpi_rank] + cv0ofa_local[ifa_local];
            if (ifa_local < nfa_i0) {
              cvofa_global[nfa][1] = cvora[mpi_rank] + cv1ofa_local_i0[ifa_local];
            }
            else {
              int rank,bits,index;
              BitUtils::unpackRankBitsIndex(rank,bits,index,rbi1ofa_local_i1[ifa_local-nfa_i0]);
              cvofa_global[nfa][1] = (cvora[rank]+index) | (int8(bits)<<55); // TODO: why 55?
            }
            my_nfa_zone[zone] += 1;
            my_noofa_zone[zone] += noofa_v.size() - nof_final_first;
            // geometry (already averaged)...
            FOR_I3 n_fa[nfa][i] = fa_dup_data_vec[fa_dup_of_fa_local[ifa_local]].n[i];
            FOR_I3 x_fa[nfa][i] = fa_dup_data_vec[fa_dup_of_fa_local[ifa_local]].x[i];
            ++nfa;
            noofa_i[nfa] = noofa_v.size();
          }
        }
      }
      delete[] fa_dup_of_fa_local;
      delete[] cv0ofa_local;
      delete[] cv1ofa_local_i0;
      delete[] rbi1ofa_local_i1;
      zone_local_pair_vec.clear();
      vector<pair<int,int> >().swap(zone_local_pair_vec);
      noofa_v_vec.clear();
      vector<int>().swap(noofa_v_vec);
      noofa_i_vec.clear();
      vector<int>().swap(noofa_i_vec);
      fa_dup_data_vec.clear();
      vector<FaData>().swap(fa_dup_data_vec);

      for (int nof = 0; nof < noofa_v.size(); ++nof)
        assert(ino_global[noofa_v[nof].first] <= -2);

      // at this point, any nodes that are used by active fa/bf's present in the final grid
      // are flipped to -2 indexing on the rank where we know they are active. So we can do a
      // reduction to reduce this negative information and make the node numbering consistent (all
      // negative if active)...

#include "compressNodes.hpp"

      for (int nof = 0; nof < noofa_v.size(); ++nof)
        assert((-ino_global[noofa_v[nof].first]-2 >= 0)&&(-ino_global[noofa_v[nof].first]-2 < nno_global_early));

      // now in the new node numbering, decide which ones we own, and index them...

      int8 my_count_reduced[nng];
      for (int i = 0; i < nng; ++i)
        my_count_reduced[i] = 0;

      for (int ino_local = 0; ino_local < nno_local_reduced; ++ino_local) {
        const int ing = PeriodicFunc::getNodeGroup(no_bits[ino_local]);
        int8 this_ino_global = my_disp[ing]; my_disp[ing] += ing_size[ing];
        for (int ing_prev = 0; ing_prev < ing; ++ing_prev)
          this_ino_global += count[ing_prev];
        if (this_ino_global == -ino_global[ino_local]-2) {
          // this is an original node of a surviving face that was selected to be active...
          my_count_reduced[ing] += ing_size[ing];
        }
      }

      // rewind my_disp and check...

      for (int i = 0; i < nng; ++i) {
        my_disp[i] -= my_count[i];
        assert(my_disp_copy[i] == my_disp[i]);
      }

      int8 count_reduced[nng];
      MPI_Allreduce(my_count_reduced,count_reduced,nng,MPI_INT8,MPI_SUM,mpi_comm);

      // the following counts may be reduced slightly if there are
      // certain node members of the large periodic groups that do not
      // appear. There is no good way that I can think of to figure this
      // out here, so we will just do a cleanup at the end, and these
      // limits should be adjusted by the removed node count. We can also
      // use this final reduction as an opportunity to remove nodes that
      // add no information (i.e. they split an edge)...

      // the periodic count...
      int8 nno_p_global = 0;
      for (int i = 0; i < 7; ++i)
        nno_p_global += count_reduced[i];

      // the boundary node count includes group 7...
      int8 nno_pb_global = nno_p_global + count_reduced[7];

      // and the total node count includes group 8...
      int8 nno_global = nno_pb_global + count_reduced[8];

      int8 my_disp_reduced[nng];
      MPI_Scan(my_count_reduced,my_disp_reduced,nng,MPI_INT8,MPI_SUM,mpi_comm);

      for (int i = 0; i < nng; ++i) {
        my_disp_reduced[i] -= my_count_reduced[i];
        my_disp_copy[i] = my_disp_reduced[i]; // switch copy check to the reduce counts
      }

      //if (mpi_rank == 0) {
      //  cout << " > nno_global: " << nno_global << " nno_p_global: " << nno_p_global << " nno_pb_global: " << nno_pb_global << endl;
      //  int8 first = 0;
      //  for (int i = 0; i < nng; ++i) {
      //    cout << " >> node group: " << i << " count: " << count_reduced[i] << " index range: " << first << " : " << first+count_reduced[i]-1 << endl;
      //    first += count_reduced[i];
      //  }
      //}

      // finally, re-index the nodes...

      for (int ino_local = 0; ino_local < nno_local_reduced; ++ino_local) {
        const int ing = PeriodicFunc::getNodeGroup(no_bits[ino_local]);
        assert(ing >= 0);
        int8 this_ino_global = my_disp[ing]; my_disp[ing] += ing_size[ing];
        for (int ing_prev = 0; ing_prev < ing; ++ing_prev)
          this_ino_global += count[ing_prev];
        if (ino_global[ino_local] >= 0) {
          // this node was not needed by any bf/fa on ANY processor...
          ino_global[ino_local] = -nno_global-10;
        }
        else if (this_ino_global == -ino_global[ino_local]-2) {
          // this is an original node that was selected to be active...
          ino_global[ino_local] = my_disp_reduced[ing]; my_disp_reduced[ing] += ing_size[ing];
          for (int ing_prev = 0; ing_prev < ing; ++ing_prev)
            ino_global[ino_local] += count_reduced[ing_prev];
          assert((ino_global[ino_local] >= 0)&&(ino_global[ino_local] < nno_global));
        }
        else {
          ino_global[ino_local] = nno_global; // large value - should not survive reduction
        }
      }

      for (int i = 0; i < nng; ++i) {
        my_disp_reduced[i] -= my_count_reduced[i];
        assert(my_disp_copy[i] == my_disp_reduced[i]);
      }

#include "compressNodes.hpp"

      delete[] send_int8_buf; send_int8_buf = NULL;
      delete[] recv_int8_buf; recv_int8_buf = NULL;
      delete[] send_int8_ino_local;
      delete[] recv_int8_ino_local;

      // =========================================================
      // build x_no. Start by building the avg,rms and count
      // at each of the reduced locations.
      // =========================================================

      double (*x_local)[3] = new double[nno_local_reduced][3];
      double (*x2_local)[3] = new double[nno_local_reduced][3];
      int *x_local_count = new int[nno_local_reduced];
      for (int ino_local = 0; ino_local < nno_local_reduced; ++ino_local) {
        FOR_I3 x_local[ino_local][i] = 0.0;
        FOR_I3 x2_local[ino_local][i] = 0.0;
        x_local_count[ino_local] = 0;
      }

      // recall nodeVec is currently sorted by icv. Since we have to add x_vd[icv] to
      // the local vdReturnData doubles in nodeVec, we can do this efficiently as follows...

      {
        int iend = 0;
        int prev_icv = -1;
        while (iend < nodeVec.size()) {
          const int istart = iend++;
          const int current_icv = nodeVec[istart].second.second;
          const double * const current_xp = x_vd[current_icv];
          assert((prev_icv+1) == current_icv);
          prev_icv = current_icv;
          while ((iend < nodeVec.size())&&(current_icv == nodeVec[iend].second.second)) ++iend;
          for (int ii = istart; ii < iend; ++ii) {
            const int ino_local = nodeVec[ii].second.first;
            FOR_I3 {
              const double this_x_i = nodeVec[ii].first[i] + current_xp[i];
              x_local[ino_local][i] += this_x_i;
              x2_local[ino_local][i] += this_x_i*this_x_i;
            }
            ++x_local_count[ino_local];
          }
        }
      }
      nodeVec.clear();
      vector<pair<double*,pair<int,int> > >().swap(nodeVec);

      // normalize local values...
      for (int ino_local = 0; ino_local < nno_local_reduced; ++ino_local) {
        assert(x_local_count[ino_local] > 0);
        const double tmp = 1.0/(double)x_local_count[ino_local];
        FOR_I3 x_local[ino_local][i] *= tmp;
        FOR_I3 x2_local[ino_local][i] = sqrt(max(0.0,x2_local[ino_local][i]*tmp - x_local[ino_local][i]*x_local[ino_local][i]));
      }

      dumpRange(x_local,nno_local_reduced,"x_local");
      dumpRange(x2_local,nno_local_reduced,"x_local RMS (should be small)");
      delete[] x2_local; x2_local = NULL;

      // finally, put the nodes into a balanced striped distribution across processors. This
      // is JUST a balanced, single striping. Here we force an alignment of 8 to ensure that
      // no individual periodic node group is broken across processors. This allows us
      // to determine pbi locally including the effects of eliminating no-touch nodes later...

      int8 * noora = NULL; calcUniformDistAlign8(noora,nno_global,mpi_size);
      int nno = noora[mpi_rank+1]-noora[mpi_rank]; // our final local node count
      // also confirm the 8 alignment for our partition...
      assert((noora[mpi_rank]%8 == 0)||(noora[mpi_rank] == noora[mpi_size]));

      double (*x_no)[3]  = new double[nno][3];
      double (*x2_no)[3]  = new double[nno][3];
      int *x_no_count = new int[nno];
      for (int ino = 0; ino < nno; ++ino) {
        FOR_I3 x_no[ino][i] = 0.0;
        FOR_I3 x2_no[ino][i] = 0.0;
        x_no_count[ino] = 0;
      }

      FOR_RANK send_int8_count[rank] = 0;

      for (int ino_local = 0; ino_local < nno_local_reduced; ++ino_local) {
        assert(x_local_count[ino_local] > 0);
        if (ino_global[ino_local] >= 0) {
          assert(ino_global[ino_local] < nno_global);
          // need to send to all permutations of this nodes bit pairs....
          const int ing = PeriodicFunc::getNodeGroup(no_bits[ino_local]);
          const int nng = PeriodicFunc::getNodeGroupSize(ing);
          for (int ii = 0; ii < nng; ++ii) {
            int bits = PeriodicFunc::getFirstSixBitsForNodeGroupMember(ing,ii);
            bits |= ((~MASK_6BITS)&no_bits[ino_local]); // need to put back any other bits (e.g., boundary)
            const int8 this_ino_global = PeriodicFunc::calcGlobalNodeNumber(ino_global[ino_local],bits,count_reduced);
            const int rank = getRankInXora(this_ino_global,noora);
            if (rank == mpi_rank) {
              // for on-rank points, just add here...
              int dbits = (bits&MASK_6BITS); // recall transforms support 12 bits, but here we use (1<<6) to indicate boundary
              FOR_I3 {
                const int bo_0 = no_bits[ino_local] & (1<<(2*i));
                const int bo_1 = no_bits[ino_local] & (1<<(2*i+1));
                const int bp_0 = bits & (1<<(2*i));
                const int bp_1 = bits & (1<<(2*i+1));
                assert(!(bo_0 && bo_1));
                assert(!(bp_0 && bp_1));
                assert((bo_0 || bo_1) == (bp_0 || bp_1) );
                // if they are diff take the permutations bit pair (bp)
                if ((bp_0 == bo_0)&&(bp_1 == bo_1))
                  dbits &= ~(MASK_2BITS<<(2*i));
              }
              double x[3]; FOR_I3 x[i] = x_local[ino_local][i];
              PeriodicData::periodicTranslate(x,1,dbits);
              const int ino = this_ino_global - noora[mpi_rank]; assert((ino >= 0)&&(ino < nno));
              FOR_I3 x_no[ino][i] += x[i]*(double)x_local_count[ino_local];
              FOR_I3 x2_no[ino][i] += x[i]*x[i]*(double)(x_local_count[ino_local]);
              x_no_count[ino] += x_local_count[ino_local];
            }
            else {
              send_int8_count[rank] += 2; // ino,count
            }
          }
        }
        else {
          assert(ino_global[ino_local] == -nno_global-10);
        }
      }

      send_int8_disp[0] = 0;
      for (int rank = 1; rank < mpi_size; ++rank)
        send_int8_disp[rank] = send_int8_count[rank-1] + send_int8_disp[rank-1];
      int send_count_sum = send_int8_disp[mpi_size-1] + send_int8_count[mpi_size-1];

      assert(send_int_buf == NULL); send_int_buf       = new int[send_count_sum];   // ino,count
      assert(send_double_buf == NULL); send_double_buf = new double[send_count_sum*3]; // x,x2

      for (int ino_local = 0; ino_local < nno_local_reduced; ++ino_local) {
        assert(x_local_count[ino_local] > 0);
        if (ino_global[ino_local] >= 0) {
          assert(ino_global[ino_local] < nno_global);
          // need to send to all permutations of this nodes bit pairs....
          const int ing = PeriodicFunc::getNodeGroup(no_bits[ino_local]);
          const int nng = PeriodicFunc::getNodeGroupSize(ing);
          for (int ii = 0; ii < nng; ++ii) {
            int bits = PeriodicFunc::getFirstSixBitsForNodeGroupMember(ing,ii);
            bits |= ((~MASK_6BITS)&no_bits[ino_local]); // need to put back any other bits (e.g., boundary)
            const int8 this_ino_global = PeriodicFunc::calcGlobalNodeNumber(ino_global[ino_local],bits,count_reduced);
            const int rank = getRankInXora(this_ino_global,noora);
            if (rank != mpi_rank) {
              // for off-rank points, fill buffers...
              int dbits = (bits&MASK_6BITS); // recall transforms support 12 bits, but here we use (1<<6) to indicate boundary
              FOR_I3 {
                const int bo_0 = no_bits[ino_local] & (1<<(2*i));
                const int bo_1 = no_bits[ino_local] & (1<<(2*i+1));
                const int bp_0 = bits & (1<<(2*i));
                const int bp_1 = bits & (1<<(2*i+1));
                assert(!(bo_0 && bo_1));
                assert(!(bp_0 && bp_1));
                assert((bo_0 || bo_1) == (bp_0 || bp_1) );
                // if they are diff take the permutations bit pair (bp)
                if ((bp_0 == bo_0)&&(bp_1 == bo_1))
                  dbits &= ~(MASK_2BITS<<(2*i));
              }
              double x[3]; FOR_I3 x[i] = x_local[ino_local][i];
              PeriodicData::periodicTranslate(x,1,dbits);
              // ino,count
              const int ino = this_ino_global - noora[rank]; assert((ino >= 0)&&(ino < noora[rank+1]-noora[rank]));
              send_int_buf[send_int8_disp[rank]  ] = ino;
              send_int_buf[send_int8_disp[rank]+1] = x_local_count[ino_local];
              //  x,x2
              FOR_I3 send_double_buf[send_int8_disp[rank]*3+i] = x[i]*(double)x_local_count[ino_local];
              FOR_I3 send_double_buf[send_int8_disp[rank]*3+3+i] = x[i]*x[i]*(double)(x_local_count[ino_local]);
              // and inc disp...
              send_int8_disp[rank] += 2;
            }
          }
        }
        else {
          assert(ino_global[ino_local] == -nno_global-10);
        }
      }

      delete[] x_local_count;
      delete[] x_local;

      // rewind first...

      send_int8_disp[0] = 0;
      for (int rank = 1; rank < mpi_size; ++rank)
        send_int8_disp[rank] = send_int8_count[rank-1] + send_int8_disp[rank-1];

      // recv-side counts....

      MPI_Alltoall(send_int8_count,1,MPI_INT,recv_int8_count,1,MPI_INT,mpi_comm);

      recv_int8_disp[0] = 0;
      for (int rank = 1; rank < mpi_size; ++rank)
        recv_int8_disp[rank] = recv_int8_count[rank-1] + recv_int8_disp[rank-1];
      int recv_count_sum = recv_int8_disp[mpi_size-1] + recv_int8_count[mpi_size-1];

      assert(recv_int_buf == NULL); recv_int_buf    = new int[recv_count_sum];   // ino,count
      MPI_Alltoallv(send_int_buf,send_int8_count,send_int8_disp,MPI_INT,
          recv_int_buf,recv_int8_count,recv_int8_disp,MPI_INT,
          mpi_comm);
      delete[] send_int_buf; send_int_buf = NULL;

      FOR_RANK {
        send_int8_count[rank] *= 3;
        send_int8_disp[rank]  *= 3;
        recv_int8_count[rank] *= 3;
        recv_int8_disp[rank]  *= 3;
      }

      assert(recv_double_buf == NULL); recv_double_buf = new double[recv_count_sum*3];   // ino,count
      MPI_Alltoallv(send_double_buf,send_int8_count,send_int8_disp,MPI_DOUBLE,
          recv_double_buf,recv_int8_count,recv_int8_disp,MPI_DOUBLE,
          mpi_comm);

      delete[] send_double_buf;

      for (int irecv = 0; irecv < recv_count_sum; irecv += 2) {
        const int ino = recv_int_buf[irecv]; assert((ino >= 0)&&(ino < nno));
        const int x_count = recv_int_buf[irecv+1]; assert(x_count > 0);
        x_no_count[ino] += x_count;
        FOR_I3 x_no[ino][i] += recv_double_buf[irecv*3+i];
        FOR_I3 x2_no[ino][i] += recv_double_buf[irecv*3+3+i];
      }
      delete[] recv_int_buf; recv_int_buf = NULL;
      delete[] recv_double_buf;

      // normalize x_no, and compute the rms...
      for (int ino = 0; ino < nno; ++ino) {
        assert(x_no_count[ino] > 0);
        const double tmp = 1.0/double(x_no_count[ino]);
        FOR_I3 x_no[ino][i] *= tmp;
        FOR_I3 x2_no[ino][i] = sqrt(max(0.0,x2_no[ino][i]*tmp - x_no[ino][i]*x_no[ino][i]));
      }
      dumpRange(x_no,nno,"x_no");
      dumpRange(x2_no,nno,"x_no RMS (should be small)");

      delete[] x_no_count; x_no_count = NULL;
      delete[] x2_no; x2_no = NULL;

      // =======================================================================
      // it is still possible for some of the possible and allowed-for
      // periodic nodes to be missing. i.e. not all 8 nodes of the 010101 variety
      // need be present, we could have just 4 of the 8, or maybe just 2?
      // so we need one last condense on the nodes. Also, we want to take this
      // opportunity to eliminate any nodes in loops that do not add any information
      // to the mesh. For boundary faces this will reduce the amount of information
      // in the nodes and could affect mesh metrics when reconstructed from these
      // nodes. In V5 io and later, the mesh metrics are written to the mles file,
      // so this is not a problem. Either with or without these additional nodes
      // in the boundary is associated with various visualization artifacts.
      // =======================================================================

      FOR_RANK send_int8_count[rank] = 0;

      assert(noobf_i[nbf] == noobf_v.size());
      int8 * noobf_v_global = new int8[noobf_i[nbf]];

      set<int8> uniqueSet;
      for (int ibf = 0; ibf < nbf; ++ibf) {
        assert(noobf_i[ibf+1] - noobf_i[ibf] >= 3);
        assert(uniqueSet.empty());
        for (int nob = noobf_i[ibf]; nob != noobf_i[ibf+1]; ++nob) {
          const int ino_local = noobf_v[nob];
          assert((ino_local >= 0)&&(ino_local < nno_local_b_reduced));
          const int bits = no_bits[ino_local];
          assert(bits & (1<<6));
          const int8 this_ino_global = PeriodicFunc::calcGlobalNodeNumber(ino_global[ino_local],bits,count_reduced);
          assert((this_ino_global >= 0)&&(this_ino_global < nno_pb_global));
          noobf_v_global[nob] = (this_ino_global | (int8(bits)<<55)); // include the bits in noobf_v_global for now
          assert((noobf_v_global[nob]&MASK_55BITS) == this_ino_global);
          uniqueSet.insert(this_ino_global);
          const int rank = getRankInXora(this_ino_global,noora);
          send_int8_count[rank] += 3; // boundaries send 1 the global node number and bits, prev and next...
        }
        assert(uniqueSet.size() >= 3);
        uniqueSet.clear();
      }

      noobf_v.clear();
      vector<int>().swap(noobf_v);
      delete[] no_bits; no_bits = NULL;

      // and build the noofa_v_global...

      assert(noofa_i[nfa] < TWO_BILLION); // use more processors if you hit this
      assert(noofa_i[nfa] == noofa_v.size());
      int8 * noofa_v_global = new int8[noofa_i[nfa]];

      for (int ifa = 0; ifa < nfa; ++ifa) {
        assert(noofa_i[ifa+1] - noofa_i[ifa] >= 3);
        assert(uniqueSet.empty());
        for (int nof = noofa_i[ifa]; nof != noofa_i[ifa+1]; ++nof) {
          const int ino_local = noofa_v[nof].first;
          assert((ino_local >= 0)&&(ino_local < nno_local_reduced));
          const int bits = noofa_v[nof].second;
          const int8 this_ino_global = PeriodicFunc::calcGlobalNodeNumber(ino_global[ino_local],bits,count_reduced);
          assert((this_ino_global >= 0)&&(this_ino_global < nno_global));
          noofa_v_global[nof] = (this_ino_global | (int8(bits)<<55)); // include the bits in noofa_v_global for now
          assert((noofa_v_global[nof]&MASK_55BITS) == this_ino_global);
          uniqueSet.insert(this_ino_global);
          const int rank = getRankInXora(this_ino_global,noora);
          send_int8_count[rank] += 3; // internal faces send 3 values, the global node number and bits, and the same for prev and next...
        }
        assert(uniqueSet.size() >= 3);
        uniqueSet.clear();
      }

      noofa_v.clear();
      vector<pair<int,int> >().swap(noofa_v);
      delete[] ino_global; ino_global = NULL;

      // size buffers for packing...

      send_int8_disp[0] = 0;
      for (int rank = 1; rank < mpi_size; ++rank)
        send_int8_disp[rank] = send_int8_count[rank-1] + send_int8_disp[rank-1];
      send_count_sum = send_int8_disp[mpi_size-1] + send_int8_count[mpi_size-1];

      assert(send_int8_buf == NULL); send_int8_buf = new int8[send_count_sum]; // ino/bits, ino_global_prev,ino_global_next

      // pack...
      // for now we pack the adjacent nbrs for bf's AND fa's, although we
      // do not allow bf node elimination at this time, but we could in the
      // future...

      for (int ibf = 0; ibf < nbf; ++ibf) {
        for (int nob = noobf_i[ibf]; nob != noobf_i[ibf+1]; ++nob) {
          const int8 this_ino_global = (noobf_v_global[nob]&MASK_55BITS);
          assert((this_ino_global >= 0)&&(this_ino_global < nno_pb_global));
          const int rank = getRankInXora(this_ino_global,noora);
          // the previous entry...
          int nob_p = nob-1; if (nob == noobf_i[ibf]) nob_p = noobf_i[ibf+1]-1;
          // the next entry...
          int nob_n = nob+1; if (nob == (noobf_i[ibf+1]-1)) nob_n = noobf_i[ibf];
          // TODO: revisit node elimination later: check uniqueness: had to remove this...
          //assert(this_ino_global != (noobf_v_global[nob_p]&MASK_55BITS));
          //assert(this_ino_global != (noobf_v_global[nob_n]&MASK_55BITS));
          //assert( ((noobf_v_global[nob_p]&MASK_55BITS) != (noobf_v_global[nob_n]&MASK_55BITS))||(noobf_i[ibf+1]-noobf_i[ibf] > 3));
          // pack...
          send_int8_buf[send_int8_disp[rank]  ] = noobf_v_global[nob]; // includes the bits
          send_int8_buf[send_int8_disp[rank]+1] = noobf_v_global[nob_p];
          send_int8_buf[send_int8_disp[rank]+2] = noobf_v_global[nob_n];
          // inc disp...
          send_int8_disp[rank] += 3;
        }
      }

      for (int ifa = 0; ifa < nfa; ++ifa) {
        for (int nof = noofa_i[ifa]; nof != noofa_i[ifa+1]; ++nof) {
          const int8 this_ino_global = noofa_v_global[nof]&MASK_55BITS;
          assert((this_ino_global >= 0)&&(this_ino_global < nno_global));
          const int rank = getRankInXora(this_ino_global,noora);
          // the previous entry...
          int nof_p = nof-1; if (nof == noofa_i[ifa]) nof_p = noofa_i[ifa+1]-1;
          // the next entry...
          int nof_n = nof+1; if (nof == (noofa_i[ifa+1]-1)) nof_n = noofa_i[ifa];
          // check uniqueness...
          assert(this_ino_global != (noofa_v_global[nof_p]&MASK_55BITS));
          assert(this_ino_global != (noofa_v_global[nof_n]&MASK_55BITS));
          //assert((noofa_v_global[nof_p]&MASK_55BITS) != (noofa_v_global[nof_n]&MASK_55BITS));
          // pack...
          send_int8_buf[send_int8_disp[rank]  ] = noofa_v_global[nof]; // includes the bits
          send_int8_buf[send_int8_disp[rank]+1] = noofa_v_global[nof_p];
          send_int8_buf[send_int8_disp[rank]+2] = noofa_v_global[nof_n];
          // inc disp...
          send_int8_disp[rank] += 3;
        }
      }

      // rewind...

      send_int8_disp[0] = 0;
      for (int rank = 1; rank < mpi_size; ++rank)
        send_int8_disp[rank] = send_int8_count[rank-1] + send_int8_disp[rank-1];

      // recv-side counts....

      MPI_Alltoall(send_int8_count,1,MPI_INT,recv_int8_count,1,MPI_INT,mpi_comm);

      recv_int8_disp[0] = 0;
      for (int rank = 1; rank < mpi_size; ++rank)
        recv_int8_disp[rank] = recv_int8_count[rank-1] + recv_int8_disp[rank-1];
      recv_count_sum = recv_int8_disp[mpi_size-1] + recv_int8_count[mpi_size-1];

      assert(recv_int8_buf == NULL); recv_int8_buf = new int8[recv_count_sum]; // pbi,ino_global_p,ino_global_n
      MPI_Alltoallv(send_int8_buf,send_int8_count,send_int8_disp,MPI_INT8,
          recv_int8_buf,recv_int8_count,recv_int8_disp,MPI_INT8,
          mpi_comm);

      // two_no_edge_nbrs will store the first two node nbrs associated with
      // this particular node. If there are more than 2 unique nodes, then these
      // values get flipped to -2...

      assert(no_bits == NULL); no_bits = new int[nno];
      int8 (*two_no_edge_nbrs)[2] = new int8[nno][2];
      for (int ino = 0; ino < nno; ++ino) {
        no_bits[ino] = -1;
        FOR_I2 two_no_edge_nbrs[ino][i] = -1;
      }

      for (int irecv = 0; irecv < recv_count_sum; irecv += 3) {
        const int8 this_ino_global = (recv_int8_buf[irecv]&MASK_55BITS);
        assert((this_ino_global >= noora[mpi_rank])&&(this_ino_global < noora[mpi_rank+1]));
        const int ino = this_ino_global - noora[mpi_rank];
        const int bits = (recv_int8_buf[irecv]>>55);
        if (no_bits[ino] == -1) no_bits[ino] = bits;
        else assert(no_bits[ino] == bits);
        const int8 this_ino_global_p = (recv_int8_buf[irecv+1]&MASK_55BITS);
        const int8 this_ino_global_n = (recv_int8_buf[irecv+2]&MASK_55BITS);
        // TODOO: see comment above - uniqueness check...
        //assert(this_ino_global != this_ino_global_p);
        //assert(this_ino_global != this_ino_global_n);
        // two-edge check...
        if (two_no_edge_nbrs[ino][0] == -1) {
          assert(two_no_edge_nbrs[ino][1] == -1);
          two_no_edge_nbrs[ino][0] = this_ino_global_p;
          two_no_edge_nbrs[ino][1] = this_ino_global_n;
        }
        else if (two_no_edge_nbrs[ino][0] >= 0) {
          assert(two_no_edge_nbrs[ino][1] >= 0);
          //assert(this_ino_global_p != this_ino_global_n);
          if (!(((this_ino_global_p == two_no_edge_nbrs[ino][0])&&(this_ino_global_n == two_no_edge_nbrs[ino][1]))||
                ((this_ino_global_p == two_no_edge_nbrs[ino][1])&&(this_ino_global_n == two_no_edge_nbrs[ino][0])))) {
            two_no_edge_nbrs[ino][0] = -2;
            two_no_edge_nbrs[ino][1] = -2;
          }
        }
        else {
          assert(two_no_edge_nbrs[ino][0] == -2); // have atleast 3 nbrs!
          assert(two_no_edge_nbrs[ino][1] == -2);
        }
      }

      // now cycle on nodes and see which ones have survived...

      for (int ing = 0; ing < nng; ++ing)
        my_count_reduced[ing] = 0; // redo node group counts

      // this is the node elimination loop. boundary nodes are known here from no_bits[ino]
      // having (1<<6) so we could avoid boundary node elimination if we wanted...

      uint8 * pbi_no = NULL;
      if (PeriodicData::isPeriodic()) {
        pbi_no = new uint8[nno];
        for (int ino = 0; ino < nno; ++ino)
          pbi_no[ino] = 0; // initialize with a zero everywhere
      }

      double my_d2_max = 0.0;
      int nno_final = 0;
      int ino_final_first = -1;
      for (int ino = 0; ino < nno; ++ino) {
        // we can always use the global node index to determine what group we are in...
        const int8 this_ino_global = noora[mpi_rank]+ino;
        int ing;
        int8 first = 0;
        for (ing = 0; ing < nng; ++ing) {
          if ((this_ino_global >= first)&&(this_ino_global < first+count_reduced[ing]))
            break;
          first += count_reduced[ing];
        }
        assert(ing < nng);
        const int igm = (this_ino_global-first)%PeriodicFunc::getNodeGroupSize(ing); // the group member in igr...
        if (igm == 0) // reset the ino_final_first: it is used to  this is the first group member
          ino_final_first = -1;

        if (no_bits[ino] == -1) { // this node was not touched by anyone...
          // the only way this should happen is if it is double or triple periodicity...
          // Actually this happens in single periodicity when you had symmetry between forming
          // points and the periodic surface (cartesian) due to degenerate cutting.
          //assert((ing == 0)||(ing == 1)||(ing == 2)||(ing == 3)); // xyz,xy,xz,yz
          assert(two_no_edge_nbrs[ino][0] == -1); // this -1 will force the skip below...
          assert(two_no_edge_nbrs[ino][1] == -1);
        }
        else {
          // this node has one or more entries in noofa_v_global somewhere, but may
          // be a 2-edge node. The no_bits have been passed from the noobf_v/noofa_v and checked
          // for self-consistency. Here we use the node group stuff as a check...
          assert(two_no_edge_nbrs[ino][0] != -1);
          assert(two_no_edge_nbrs[ino][1] != -1);
          const int bits = PeriodicFunc::getFirstSixBitsForNodeGroupMember(ing,igm);
          assert(bits == (no_bits[ino]&MASK_6BITS)); // this is a good node group consistency check
          // this node appears in the face loops, but will still get eliminated if it is a 2-edge node...
          // HACK HACK HACK to keep these nodes...
          two_no_edge_nbrs[ino][0] = -2;
          two_no_edge_nbrs[ino][1] = -2;
          // end of HACK HACK HACK
          if (two_no_edge_nbrs[ino][0] >= 0) {
            assert(two_no_edge_nbrs[ino][1] >= 0);
            // this node has two and only two valid node nbrs, meaning it is a two-edge node. Use
            // -1 to skip this node in the noobf_v_global/noofa_v_global loops below
            two_no_edge_nbrs[ino][0] = -1;
          }
          else {
            // this is a valid node. It must be in several bf or fa loops...
            assert(two_no_edge_nbrs[ino][0] == -2);
            assert(two_no_edge_nbrs[ino][1] == -2);
            // reuse two_no_edge_nbrs to store the new node index...
            const int ino_final = nno_final++;
            two_no_edge_nbrs[ino][0] = ino_final; // gets updated with offset below
            ++my_count_reduced[ing]; // update node group counts
            // copy down coordinates and bits...
            if (ino_final != ino) {
              assert(ino_final < ino);
              FOR_I3 x_no[ino_final][i] = x_no[ino][i];
              no_bits[ino_final]        = no_bits[ino];
            }
            // finally, if ino_final_first is -1, then this is the first node...
            if (ino_final_first == -1) {
              //cout << "SETTING ino_final_first: bits: "; dumpBits(bits); cout << " position: " << COUT_VEC(x_no[ino_final]) << endl;
              ino_final_first = ino_final;
            }
            else {
              //cout << " > CHILD of ino_final:   bits: "; dumpBits(bits); cout << " position: " << COUT_VEC(x_no[ino_final]) << endl;
              // the only time we should get children is with periodicity...
              assert(pbi_no); assert(pbi_no[ino_final] == 0);
              assert(bits);
              // the pbi for each node child stores a reference to the
              // parent node, and the transform necessary to get FROM the
              // parent location TO the child...
              int pbi_bits = 0;
              FOR_I3 {
                const int parent_bit_pair = ((no_bits[ino_final_first]>>(2*i))&MASK_2BITS);
                const int child_bit_pair = ((bits>>(2*i))&MASK_2BITS);
                if (parent_bit_pair == 1) {
                  if (child_bit_pair == 2) {
                    pbi_bits |= (1<<(2*i+1));
                  }
                  else {
                    assert(child_bit_pair == 1);
                  }
                }
                else if (parent_bit_pair == 2) {
                  if (child_bit_pair == 1) {
                    pbi_bits |= (1<<(2*i));
                  }
                  else {
                    assert(child_bit_pair == 2);
                  }
                }
                else {
                  assert(parent_bit_pair == 0);
                  assert(child_bit_pair == 0);
                }
              }
              // check what you get when you apply transform to parent...
              double xp_t[3]; FOR_I3 xp_t[i] = x_no[ino_final_first][i];
              PeriodicData::periodicTranslate(xp_t,1,pbi_bits);
              const double this_d2 = DIST2(xp_t,x_no[ino_final]);
              my_d2_max = max(my_d2_max,this_d2);
              pbi_no[ino_final] = BitUtils::packPbiHash(pbi_bits,ino_final_first); // use the local node for now and adjust below when we have the offset
            }
          }
        }
      }

      // the final node distribution...
      int8 * noora_final = NULL;
      buildXora(noora_final,nno_final);
      assert(nno_final == noora_final[mpi_rank+1]-noora_final[mpi_rank]);

      // add local offset to the active nodes and prepare to send back to
      // finalize noobf_v_global and noofa_v_global...

      int ino_check = 0;
      for (int ino = 0; ino < nno; ++ino) {
        if (two_no_edge_nbrs[ino][0] >= 0) {
          assert(two_no_edge_nbrs[ino][0] == ino_check);
          two_no_edge_nbrs[ino][0] += noora_final[mpi_rank];
          ++ino_check;
        }
        else {
          assert(two_no_edge_nbrs[ino][0] == -1);
        }
      }
      assert(ino_check == nno_final);

      // adjust pbi and report its match...

      if (pbi_no) {

        double d2_max;
        MPI_Reduce(&my_d2_max,&d2_max,1,MPI_DOUBLE,MPI_MAX,0,mpi_comm);
        if (mpi_rank == 0)
          cout << " > pbi_no match (should be zero): " << sqrt(d2_max) << endl;

        for (int ino = 0; ino < nno_final; ++ino) {
          const int8 this_ino_global = noora_final[mpi_rank] + ino;
          if (pbi_no[ino] == 0) {
            pbi_no[ino] = BitUtils::packPbiHash(0,this_ino_global);
            assert(pbi_no[ino] == uint8(this_ino_global));
            //cout << this_ino_global << " " << this_ino_global << " " << endl;
          }
          else {
            int bits;
            int8 ino_parent;
            BitUtils::unpackPbiHash(bits,ino_parent,pbi_no[ino]);
            assert(bits);
            assert(ino_parent < ino);
            const int8 ino_parent_global = noora_final[mpi_rank] + ino_parent;
            pbi_no[ino] = BitUtils::packPbiHash(bits,ino_parent_global);
            //cout << ino_parent_global << " " << this_ino_global << " " << endl;
          }

        }

      }

      // put the new global index back into the recv buffer...

      for (int irecv = 0; irecv < recv_count_sum; irecv += 3) {
        const int8 this_ino_global = (recv_int8_buf[irecv]&MASK_55BITS);
        assert((this_ino_global >= noora[mpi_rank])&&(this_ino_global < noora[mpi_rank+1]));
        const int ino = this_ino_global - noora[mpi_rank];
        recv_int8_buf[irecv/3] = two_no_edge_nbrs[ino][0]; // could be -1
      }

      delete[] two_no_edge_nbrs;
      delete[] no_bits;

      FOR_RANK {
        send_int8_count[rank] /= 3;
        send_int8_disp[rank]  /= 3;
        recv_int8_count[rank] /= 3;
        recv_int8_disp[rank]  /= 3;
      }

      // and back to the send side...

      MPI_Alltoallv(recv_int8_buf,recv_int8_count,recv_int8_disp,MPI_INT8,
          send_int8_buf,send_int8_count,send_int8_disp,MPI_INT8,
          mpi_comm);
      delete[] recv_int8_buf;
      delete[] send_int8_count; // still need send_int8_disp to unpack below...
      delete[] recv_int8_count;
      delete[] recv_int8_disp;

      // now modify the global references, compressing as we go...

      for (int izone = 0; izone < nzone; ++izone)
        my_noobf_zone[izone] = 0;

      int noobf_s = 0;
      int ibf_end = 0;
      for (int izone = 0; izone < nzone; ++izone) {
        const int ibf_begin = ibf_end;
        ibf_end += my_nbf_zone[izone];
        for (int ibf = ibf_begin; ibf != ibf_end; ++ibf) {
          const int nob_f = noobf_i[ibf];
          noobf_i[ibf] = noobf_s;
          for (int nob = nob_f; nob != noobf_i[ibf+1]; ++nob) {
            const int8 this_ino_global = noobf_v_global[nob]&MASK_55BITS;
            assert((this_ino_global >= 0)&&(this_ino_global < nno_pb_global));
            const int rank = getRankInXora(this_ino_global,noora);
            const int8 this_ino_global_new = send_int8_buf[send_int8_disp[rank]++];
            if (this_ino_global_new >= 0) {
              assert(this_ino_global_new <= this_ino_global);
              noobf_v_global[noobf_s++] = this_ino_global_new;
              ++my_noobf_zone[izone];
            }
            else {
              assert(this_ino_global_new == -1);
            }
          }
        }
      }
      noobf_i[nbf] = noobf_s;

      for (int izone = 0; izone < 27; ++izone)
        my_noofa_zone[izone] = 0;

      int noofa_s = 0;
      int ifa_end = 0;
      for (int izone = 0; izone < 27; ++izone) {
        const int ifa_begin = ifa_end;
        ifa_end += my_nfa_zone[izone];
        for (int ifa = ifa_begin; ifa != ifa_end; ++ifa) {
          const int nof_f = noofa_i[ifa];
          noofa_i[ifa] = noofa_s;
          for (int nof = nof_f; nof != noofa_i[ifa+1]; ++nof) {
            const int8 this_ino_global = noofa_v_global[nof]&MASK_55BITS;
            assert((this_ino_global >= 0)&&(this_ino_global < nno_global));
            const int rank = getRankInXora(this_ino_global,noora);
            const int8 this_ino_global_new = send_int8_buf[send_int8_disp[rank]++];
            if (this_ino_global_new >= 0) {
              assert(this_ino_global_new <= this_ino_global);
              noofa_v_global[noofa_s++] = this_ino_global_new;
              ++my_noofa_zone[izone];
            }
            else {
              assert(this_ino_global_new == -1);
            }
          }
        }
      }
      noofa_i[nfa] = noofa_s;

      delete[] send_int8_buf;
      delete[] send_int8_disp;

      // update noora, nno...
      delete[] noora;
      noora = noora_final;
      nno = nno_final;
      assert(nno == noora[mpi_rank+1]-noora[mpi_rank]);

      // reduce node group counts one last time...
      MPI_Allreduce(my_count_reduced,count_reduced,nng,MPI_INT8,MPI_SUM,mpi_comm);

      // the periodic count...
      nno_p_global = 0;
      for (int i = 0; i < 7; ++i)
        nno_p_global += count_reduced[i];

      // the boundary node count includes group 7...
      nno_pb_global = nno_p_global + count_reduced[7];

      // and the total node count includes group 8...
      nno_global = nno_pb_global + count_reduced[8];
      assert(nno_global == noora[mpi_size]);

      MPI_Scan(my_count_reduced,my_disp_reduced,nng,MPI_INT8,MPI_SUM,mpi_comm);
      for (int i = 0; i < nng; ++i) {
        my_disp_reduced[i] -= my_count_reduced[i];
      }
      if (mpi_rank == 0) {
        cout << " > final nno_global: " << nno_global << " nno_p_global: " << nno_p_global << " nno_pb_global: " << nno_pb_global << endl;
        int8 first = 0;
        for (int i = 0; i < nng; ++i) {
          cout << " > node group: " << i << " count: " << count_reduced[i] << endl;
          first += count_reduced[i];
        }
      }

      // ======================================
      // serial tecplot hack...
      // ======================================

      /*
         if (mpi_size == 1) {

         assert(mpi_size == 1);

         char filename[128];
         sprintf(filename,"debug.dat");

         FILE * fp = fopen(filename,"w");
         fprintf(fp,"TITLE = \"%s\"\n","debug.dat");
         fprintf(fp,"VARIABLES = \"X\"\n");
         fprintf(fp,"\"Y\"\n");
         fprintf(fp,"\"Z\"\n");

         assert(nzone == surface->zoneVec.size());
         int ibf_end = 0;
         for (int izone = 0; izone < nzone; ++izone) {
         if (my_nbf_zone[izone] > 0) {
         fprintf(fp,"ZONE T=\"%s\"\n",surface->zoneVec[izone].getName().c_str());
         fprintf(fp,"N=%lld, E=%lld, F=FEPOINT, ET=TRIANGLE\n",my_noobf_zone[izone]+my_nbf_zone[izone],my_noobf_zone[izone]);
         const int ibf_begin = ibf_end;
         ibf_end += my_nbf_zone[izone];
      // nodes...
      for (int ibf = ibf_begin; ibf != ibf_end; ++ibf) {
      assert(noobf_i[ibf+1]-noobf_i[ibf] >= 3);
      // write the x_bf...
      fprintf(fp,"%18.15e %18.15e %18.15e\n",x_bf[ibf][0],x_bf[ibf][1],x_bf[ibf][2]);
      for (int nob = noobf_i[ibf]; nob != noobf_i[ibf+1]; ++nob) {
      const int8 this_ino_global = noobf_v_global[nob]&MASK_55BITS;
      const int ino = this_ino_global; assert((ino >= 0)&&(ino < nno)); // one rank only...
      fprintf(fp,"%18.15e %18.15e %18.15e\n",x_no[ino][0],x_no[ino][1],x_no[ino][2]);
      }
      }
      // tris...
      int nno_offset = 0;
      for (int ibf = ibf_begin; ibf != ibf_end; ++ibf) {
      int ino_bf = nno_offset;
      int ino_0 = ino_bf + (noobf_i[ibf+1]-noobf_i[ibf]);
      for (int nob = noobf_i[ibf]; nob != noobf_i[ibf+1]; ++nob) {
      const int ino_1 = ino_bf + 1 + (nob-noobf_i[ibf]);
      fprintf(fp,"%d %d %d\n",ino_bf+1,ino_0+1,ino_1+1); // 1-indexed
      ino_0 = ino_1;
      }
      nno_offset += 1 + (noobf_i[ibf+1]-noobf_i[ibf]);
      }
      }
      }

      int ifa_end = 0;
      for (int izone = 0; izone < 26; ++izone) { // i.e. skip internal
      if (my_nfa_zone[izone] > 0) {
      fprintf(fp,"ZONE T=\"%d\"\n",izone);
      fprintf(fp,"N=%lld, E=%lld, F=FEPOINT, ET=TRIANGLE\n",my_noofa_zone[izone]+my_nfa_zone[izone],my_noofa_zone[izone]);
      const int ifa_begin = ifa_end;
      ifa_end += my_nfa_zone[izone];
      // nodes...
      for (int ifa = ifa_begin; ifa != ifa_end; ++ifa) {
      assert(noofa_i[ifa+1]-noofa_i[ifa] >= 3);
      // write the x_fa...
      fprintf(fp,"%18.15e %18.15e %18.15e\n",x_fa[ifa][0],x_fa[ifa][1],x_fa[ifa][2]);
      for (int nof = noofa_i[ifa]; nof != noofa_i[ifa+1]; ++nof) {
      const int8 this_ino_global = noofa_v_global[nof]&MASK_55BITS;
      const int ino = this_ino_global; assert((ino >= 0)&&(ino < nno)); // one rank only...
      fprintf(fp,"%18.15e %18.15e %18.15e\n",x_no[ino][0],x_no[ino][1],x_no[ino][2]);
      }
      }
      // tris...
      int nno_offset = 0;
      for (int ifa = ifa_begin; ifa != ifa_end; ++ifa) {
      int ino_fa = nno_offset;
      int ino_0 = ino_fa + (noofa_i[ifa+1]-noofa_i[ifa]);
      for (int nof = noofa_i[ifa]; nof != noofa_i[ifa+1]; ++nof) {
        const int ino_1 = ino_fa + 1 + (nof-noofa_i[ifa]);
        fprintf(fp,"%d %d %d\n",ino_fa+1,ino_0+1,ino_1+1); // 1-indexed
        ino_0 = ino_1;
      }
      nno_offset += 1 + (noofa_i[ifa+1]-noofa_i[ifa]);
    }
    }
    }

    fclose(fp);

    }
    */

      // ======================================
      // start to write...
      // ======================================

      // copy LegacyIO from dev May 8, 2017 for start of the following...

      // convert noobf_i to a noobf_count...

      assert(noobf_i[nbf] == noobf_s);
    for (int ibf = 0; ibf < nbf; ++ibf)
      noobf_i[ibf] = noobf_i[ibf+1]-noobf_i[ibf];
    int * noobf_count = noobf_i;

    // convert spobf_i to a spobf_count...

    const int spobf_s = spobf_i[nbf];
    assert(spobf_s == spobf_v.size());
    assert(spobf_s == spobf_wgt.size());
    for (int ibf = 0; ibf < nbf; ++ibf)
      spobf_i[ibf] = spobf_i[ibf+1]-spobf_i[ibf];
    int * spobf_count = spobf_i;

    // convert sbobf_i to a sbobf_count...

    const int sbobf_s = sbobf_i[nbf];
    assert(sbobf_s == psbobf_v.size());
    for (int ibf = 0; ibf < nbf; ++ibf)
      sbobf_i[ibf] = sbobf_i[ibf+1]-sbobf_i[ibf];
    int * sbobf_count = sbobf_i;

    // same with noofa_i...

    assert(noofa_i[nfa] == noofa_s);
    for (int ifa = 0; ifa < nfa; ++ifa)
      noofa_i[ifa] = noofa_i[ifa+1]-noofa_i[ifa];
    int * noofa_count = noofa_i;

    // =====================================================================
    // we now have SOME the following ready for writing...
    //
    // int ncv = mainPart->ncv;
    // int8 * cvora
    //
    // int * noobf_count = new int[nbf];
    // int8 * noobf_v_global = new int8[noobf_s];
    // int8 * cvobf_global = new int8[nbf];
    //
    // NOTE that not all zones in the surface->zoneVec may be represented, and
    // also orphans may or may not be present
    //
    // int * noofa_count = new int[nfa_pi];
    // int8 * noofa_v_global = new int8[noofa_pi_s];
    // int8 (*cvofa_global)[2] = new int8[nfa_pi][2];
    //
    // int8 * noora
    // int8 nno_global (== noora[mpi_size])
    // int8 nno_p_global // the periodic nodes -- it would be good to include this in io version 4
    // int8 nno_pb_global // the periodic and boundary nodes -- also good to include in io version 4
    // int nno // == noora[mpi_rank+1]-noora[mpi_rank];
    // double (*x_no)[3] = new double[nno][3];
    // ...
    // other?...
    // =====================================================================

    // switch the meaning of the disp to stripe the bf's and fa's...

    int8 * my_nbf_disp = new int8[nzone];
    MPI_Scan(my_nbf_zone,my_nbf_disp,nzone,MPI_INT8,MPI_SUM,mpi_comm);
    int8 * nbf_zone = new int8[nzone];
    if (mpi_rank == mpi_size-1)
      for (int izone = 0; izone < nzone; ++izone)
        nbf_zone[izone] = my_nbf_disp[izone];
    MPI_Bcast(nbf_zone,nzone,MPI_INT8,mpi_size-1,mpi_comm);
    for (int izone = 0; izone < nzone; ++izone)
      my_nbf_disp[izone] -= my_nbf_zone[izone];

    /*
       cout << "nzone: " << nzone << endl;
       for (int izone = 0; izone < nzone; ++izone)
       cout << " > my_nbf_zone[izone]: " << my_nbf_zone[izone] << endl;
       getchar();
       */

    int8 * my_noobf_disp = new int8[nzone];
    MPI_Scan(my_noobf_zone,my_noobf_disp,nzone,MPI_INT8,MPI_SUM,mpi_comm);
    int8 * noobf_zone = new int8[nzone];
    if (mpi_rank == mpi_size-1)
      for (int izone = 0; izone < nzone; ++izone)
        noobf_zone[izone] = my_noobf_disp[izone];
    MPI_Bcast(noobf_zone,nzone,MPI_INT8,mpi_size-1,mpi_comm);
    for (int izone = 0; izone < nzone; ++izone)
      my_noobf_disp[izone] -= my_noobf_zone[izone];

    int8 * my_spobf_disp = NULL;
    int8 * spobf_zone = NULL;
    if (write_surface) {
      my_spobf_disp = new int8[nzone];
      MPI_Scan(my_spobf_zone,my_spobf_disp,nzone,MPI_INT8,MPI_SUM,mpi_comm);
      spobf_zone = new int8[nzone];
      if (mpi_rank == mpi_size-1)
        for (int izone = 0; izone < nzone; ++izone)
          spobf_zone[izone] = my_spobf_disp[izone];
      MPI_Bcast(spobf_zone,nzone,MPI_INT8,mpi_size-1,mpi_comm);
      for (int izone = 0; izone < nzone; ++izone)
        my_spobf_disp[izone] -= my_spobf_zone[izone];
    }

    int8 * my_sbobf_disp = NULL;
    int8 * sbobf_zone = NULL;
    if (write_surface) {
      my_sbobf_disp = new int8[nzone];
      MPI_Scan(my_sbobf_zone,my_sbobf_disp,nzone,MPI_INT8,MPI_SUM,mpi_comm);
      sbobf_zone = new int8[nzone];
      if (mpi_rank == mpi_size-1)
        for (int izone = 0; izone < nzone; ++izone)
          sbobf_zone[izone] = my_sbobf_disp[izone];
      MPI_Bcast(sbobf_zone,nzone,MPI_INT8,mpi_size-1,mpi_comm);
      for (int izone = 0; izone < nzone; ++izone)
        my_sbobf_disp[izone] -= my_sbobf_zone[izone];
    }

    int8 nbf_global = 0;
    int8 noobf_global = 0;
    int8 spobf_global = 0;
    int8 sbobf_global = 0;
    for (int izone = 0; izone < nzone; ++izone) {
      nbf_global += nbf_zone[izone];
      noobf_global += noobf_zone[izone];
      if (write_surface) {
        spobf_global += spobf_zone[izone];
        sbobf_global += sbobf_zone[izone];
      }
    }

    // nfa zone counting...

    int8 * my_nfa_disp = new int8[27];
    MPI_Scan(my_nfa_zone,my_nfa_disp,27,MPI_INT8,MPI_SUM,mpi_comm);
    int8 * nfa_zone = new int8[27];
    if (mpi_rank == mpi_size-1)
      for (int izone = 0; izone < 27; ++izone)
        nfa_zone[izone] = my_nfa_disp[izone];
    MPI_Bcast(nfa_zone,27,MPI_INT8,mpi_size-1,mpi_comm);
    for (int izone = 0; izone < 27; ++izone)
      my_nfa_disp[izone] -= my_nfa_zone[izone];

    int8 * my_noofa_disp = new int8[27];
    MPI_Scan(my_noofa_zone,my_noofa_disp,27,MPI_INT8,MPI_SUM,mpi_comm);
    int8 * noofa_zone = new int8[27];
    if (mpi_rank == mpi_size-1)
      for (int izone = 0; izone < 27; ++izone)
        noofa_zone[izone] = my_noofa_disp[izone];
    MPI_Bcast(noofa_zone,27,MPI_INT8,mpi_size-1,mpi_comm);
    for (int izone = 0; izone < 27; ++izone)
      my_noofa_disp[izone] -= my_noofa_zone[izone];

    int8 nfa_global = 0;
    int8 noofa_global = 0;
    for (int izone = 0; izone < 27; ++izone) {
      nfa_global += nfa_zone[izone];
      noofa_global += noofa_zone[izone];
    }

    if (mpi_rank == 0)
      cout << " > nbf_global: " << nbf_global << " avg noobf: " << double(noobf_global)/double(nbf_global) <<
        " nfa_global: " << nfa_global << " avg noofa: " << double(noofa_global)/double(nfa_global) << endl;

    // create tmp file to write too...
    string tmp_filename = MiscUtils::makeTmpPrefix(filename);
    char fname[128];
    sprintf(fname,"%s",tmp_filename.c_str());
    MPI_File_delete(fname,MPI_INFO_NULL);

    if (io_version >= 4) {

      dumpRange(area_over_delta_bf,nbf,"area_over_delta_bf (should be non-negative and bounded)");

      // HACK to rotate mesh prior to write...
      if (Param * param = getParam("ROTATE_MLES")) {

        double axis[3] = {0.0,0.0,0.0};
        double point[3] = {0.0,0.0,0.0};
        double angle = 0.0;

        int iarg = 0;
        while (iarg < param->size()) {
          const string token = param->getString(iarg++);
          if (token == "AXIS") {
            FOR_I3 axis[i] = param->getDouble(iarg++);
            NORMALIZE(axis);
          }
          else if (token == "X") {
            axis[0] = 1.0;
            axis[1] = 0.0;
            axis[2] = 0.0;
          }
          else if (token == "Y") {
            axis[0] = 0.0;
            axis[1] = 1.0;
            axis[2] = 0.0;
          }
          else if (token == "Z") {
            axis[0] = 0.0;
            axis[1] = 0.0;
            axis[2] = 1.0;
          }
          else if (token == "POINT") {
            FOR_I3 point[i] = param->getDouble(iarg++);
          }
          else if (token == "ANGLE") {
            angle = param->getDouble(iarg++)/180.0*M_PI; // radians;
          }
        }

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

        // no stuff...

        FOR_INO {
          const double dx[3] = DIFF(x_no[ino],point);
          const double temp[3] = MATRIX_PRODUCT(R,dx);
          FOR_J3 x_no[ino][j] = temp[j] + point[j];
        }

        // cv stuff...

        FOR_ICV {
          const double dx[3] = DIFF(x_cv[icv],point);
          const double temp[3] = MATRIX_PRODUCT(R,dx);
          FOR_J3 x_cv[icv][j] = temp[j] + point[j];
        }
        FOR_ICV {
          const double dx[3] = DIFF(x_vd[icv],point);
          const double temp[3] = MATRIX_PRODUCT(R,dx);
          FOR_J3 x_vd[icv][j] = temp[j] + point[j];
        }

        // fa stuff...

        FOR_IFA {
          const double dx[3] = DIFF(x_fa[ifa],point);
          const double temp[3] = MATRIX_PRODUCT(R,dx);
          FOR_J3 x_fa[ifa][j] = temp[j] + point[j];
        }
        FOR_IFA {
          const double temp[3] = MATRIX_PRODUCT(R,n_fa[ifa]);
          FOR_J3 n_fa[ifa][j] = temp[j];
        }

        // bf stuff...

        FOR_IBF {
          const double dx[3] = DIFF(x_bf[ibf],point);
          const double temp[3] = MATRIX_PRODUCT(R,dx);
          FOR_J3 x_bf[ibf][j] = temp[j] + point[j];
        }
        FOR_IBF {
          const double temp[3] = MATRIX_PRODUCT(R,n_bf[ibf]);
          FOR_J3 n_bf[ibf][j] = temp[j];
        }
        FOR_IBF {
          // G' = RGR^T
          const double temp[3][3] = MAT_TMAT_MULT(Gij_bf[ibf],R);
          const double temp2[3][3] = MAT_MAT_MULT(R,temp);
          FOR_I3 FOR_J3 Gij_bf[ibf][i][j] = temp2[i][j];
        }

        // surface stuff...
        if (mpi_rank_shared == 0) {
          for (int ipart = 0; ipart < partVec.size(); ++ipart) {
            if (partVec[ipart]->surface) {

              SurfaceShm * surface = partVec[ipart]->surface;
              for (int isp = 0; isp < surface->nsp; ++isp) {
                const double _xsp[3] = DIFF(surface->xsp[isp],point);
                const double temp[3] = MATRIX_PRODUCT(R,_xsp);
                FOR_J3 surface->xsp[isp][j] = temp[j] + point[j];
              }
            }
          }
        }
        MPI_Barrier(mpi_comm_shared);

        // periodic stuff...
        for (int i = 0; i < PeriodicData::periodicTransformVec.size(); ++i) {
          const int kind = PeriodicData::periodicTransformVec[i].getKind();
          double data[3]; PeriodicData::periodicTransformVec[i].getData(data);
          if (kind == PERIODIC_TRANSFORM_CART) {
            const double temp[3] = MATRIX_PRODUCT(R,data);
            PeriodicData::periodicTransformVec[i].setData(temp);
          }
          else {
            CERR(" > currently only CART periodicity permits rotation");
          }
        }

      }

      // ================================================================================
      // ================================================================================
      // ================================================================================
      // io version 4
      // full support for periodicity, 2B, etc...
      // ================================================================================
      // ================================================================================
      // ================================================================================

      MPI_File fh;
      MPI_Offset offset = 0;

      MiscUtils::mkdir_for_file_collective(filename,0);

      MPI_File_open(mpi_comm,fname,MPI_MODE_WRONLY|MPI_MODE_CREATE,MPI_INFO_NULL,&fh);

      if ( mpi_rank == 0 )  {
        int itmp[2] = {UGP_IO_MAGIC_NUMBER, io_version};
        cout << " > UGP_IO_VERSION: " << itmp[1] << endl;
        MPI_File_write(fh,itmp,2,MPI_INT,MPI_STATUS_IGNORE);
      }

      offset += int_size*2;

      //=============================================
      // start with a bunch of stuff on rank0 only...
      //=============================================

      if (mpi_rank == 0) {

        //-----------
        // time...
        //-----------

        {

          time_t raw_time = std::time(NULL);
          tm * now = gmtime(&raw_time);

          char timestamp[CTI_TIMESTAMP_SIZE];
          CTI_Make_timestamp(timestamp,now);
          cout << " > timestamp (GMT) " << timestamp << endl;

          Header header;
          sprintf(header.name,"UGP_IO_TIMESTAMP");
          header.id = UGP_IO_TIMESTAMP;
          header.skip = header_size;
          header.idata[0] = now->tm_sec;   // seconds of minutes from 0 to 61
          header.idata[1] = now->tm_min;   // minutes of hour from 0 to 59
          header.idata[2] = now->tm_hour;  // hours of day from 0 to 24
          header.idata[3] = now->tm_mday;  // day of month from 1 to 31
          header.idata[4] = now->tm_mon;   // month of year from 0 to 11
          header.idata[5] = now->tm_year;  // year since 1900
          header.idata[6] = now->tm_wday;  // days since sunday
          header.idata[7] = now->tm_yday;  // days since January 1st
          header.idata[8] = now->tm_isdst; // hours of daylight savings time

          //Store the approximate runtime in the header
          //header.rdata[0] = (MPI_Wtime()-mpi_wtime_0); // MPI_Wtime reports seconds as a double...

          MPI_File_write(fh,&header,1,MPI_Header,MPI_STATUS_IGNORE);

          offset += header.skip;

        }

        //-----------
        // params...
        //-----------

        {

          cout << " > params" << endl;

          std::stringstream ss;
          CTI_Dump_solver_info(ss);
          dumpParams(ss);

          Header header;
          sprintf(header.name,"UGP_IO_PARAMS");
          header.id = UGP_IO_PARAMS;
          header.idata[0] = ss.str().length();
          header.skip = header_size + header.idata[0];

          MPI_File_write(fh,&header,1,MPI_Header,MPI_STATUS_IGNORE);
          MPI_File_write(fh,(void*)ss.str().c_str(),header.idata[0],MPI_CHAR,MPI_STATUS_IGNORE);

          offset += header.skip;

        }

        //-----------
        // counts
        //-----------

        {

          cout << " > counts" << endl;

          Header header;
          sprintf(header.name,"no_fa_bf_cv_counts");
          header.id     = UGP_IO_NO_FA_BF_CV_COUNTS;
          header.skip   = header_size + int8_size*54;
          ByteSwap::setLswMswPairForInt8(header.idata+0,noora[mpi_size]);
          ByteSwap::setLswMswPairForInt8(header.idata+2,nbf_global);
          ByteSwap::setLswMswPairForInt8(header.idata+4,nfa_global);
          ByteSwap::setLswMswPairForInt8(header.idata+6,cvora[mpi_size]);
          ByteSwap::setLswMswPairForInt8(header.idata+8,nno_p_global); // first nodes are periodic
          ByteSwap::setLswMswPairForInt8(header.idata+10,nno_pb_global); // and boundary
          // use the next idata to indicate if the surface is written...
          header.idata[12] = 0;
          if (write_surface)
            header.idata[12] = 1;
          // and write...
          MPI_File_write(fh,&header,1,MPI_Header,MPI_STATUS_IGNORE);
          // add the 27*2 int8's that describe the periodic distribution of faces...
          MPI_File_write(fh,nfa_zone,27,MPI_INT8,MPI_STATUS_IGNORE);
          MPI_File_write(fh,noofa_zone,27,MPI_INT8,MPI_STATUS_IGNORE);

          offset += header_size + int8_size*54;

        }

        //-----------
        // hashing...
        //   permissive for now, in the sense that we should be able to use an sles
        //   file with any mles as long as the counts mesh. If we generate an mles
        //   on different resources/core counts it will get the same hash if its
        //   counts match.
        //   DP -- this may be too permissive. what about date for mles?
        //-----------

        {
          RestartHashUtilities::clearHashes();
          std::stringstream ss;
          //CTI_Dump_solver_info(ss);  //code version, build date, core count dependency...
          //dumpParams(ss);            //includes all parameters...inputs, image writing...
          //only include counts in hash for now
          ss << noora[mpi_size] << nbf_global << nfa_global << cvora[mpi_size] << nno_p_global << nno_pb_global;
          RestartHashUtilities::myHash.init(ss,RestartHashUtilities::sha1hashlength);
          RestartHashUtilities::mlesHash.init(ss,RestartHashUtilities::sha1hashlength); //TODO, hash surface properties for mles parent
          cout << " > hash " << RestartHashUtilities::myHash << endl;
          offset += RestartHashUtilities::writeHashToRestartHeader(fh); //write and update offset
        }

      }

      // now sync everyone else at this offset...

      MPI_Bcast(&offset,1,MPI_OFFSET_DATATYPE,0,mpi_comm);

      //================================
      // boundary face records.
      // recall that some of the zones have not been touched.
      //================================

      {

        // globally reduce the bf zone geomtric data...

        assert(nzone == zoneVec.size());
        double (*geom_bf_zone)[8] = NULL;
        if (mpi_rank == 0) geom_bf_zone = new double[nzone][8];
        MPI_Reduce(my_geom_bf_zone,geom_bf_zone,nzone*8,MPI_DOUBLE,MPI_SUM,0,mpi_comm);

        int nzone_final = 0;
        int8 ibf_end = 0;
        int8 noobf_end = 0;
        int8 spobf_end = 0;
        int8 sbobf_end = 0;
        for (int izone = 0; izone < nzone; ++izone) {
          // could have an entire zone without noobf -- e.g. 1D circular pipe, periodic zone
          if (nbf_zone[izone] > 0) {
            assert(noobf_zone[izone] > 0);
            // this zone needs to be written...
            if ( mpi_rank == 0 ) {
              const int8 ibf_begin = ibf_end; ibf_end += nbf_zone[izone];
              const int8 noobf_begin = noobf_end; noobf_end += noobf_zone[izone];
              Header header;
              header.id       = UGP_IO_BF_ZONE_HEADER;
              header.skip     = header_size;
              header.idata[0] = FA_ZONE_BOUNDARY;
              //sprintf(header.name,"%s",zoneVec[izone].getName().c_str());
              snprintf(header.name,UGP_IO_HEADER_NAME_LEN-1,"%s",zoneVec[izone].getName().c_str());
              cout << " > bf zone \"" << header.name << "\"" << endl;
              header.idata[1] = nzone_final; // can be different than izone if zones gets eliminated - e.g. periodic zones
              ByteSwap::setLswMswPairForInt8(header.idata+2,nbf_zone[izone]);
              ByteSwap::setLswMswPairForInt8(header.idata+4,noobf_zone[izone]);
              ByteSwap::setLswMswPairForInt8(header.idata+6,ibf_begin);
              ByteSwap::setLswMswPairForInt8(header.idata+8,noobf_begin);
              // add count and begin for the surface-tri-of-boundary-face...
              if (write_surface) {
                const int8 sbobf_begin = sbobf_end; sbobf_end += sbobf_zone[izone];
                ByteSwap::setLswMswPairForInt8(header.idata+10,sbobf_zone[izone]);
                ByteSwap::setLswMswPairForInt8(header.idata+12,sbobf_begin);
                // for the spobf stuff, use the uint8 space...
                const int8 spobf_begin = spobf_end; spobf_end += spobf_zone[izone];
                header.ui8data[0] = spobf_zone[izone];
                header.ui8data[1] = spobf_begin;
              }
              // recall...
              // 0: sum(area_bf)
              // 1: sum(area_over_delta_bf)
              // 2: sum(n_bf[0])
              // 3: sum(n_bf[1])
              // 4: sum(n_bf[2])
              // 5: x_bf[0] avg (eventually)
              // 6: x_bf[1] avg (eventually)
              // 7: x_bf[2] avg (eventually)
              header.rdata[0] = geom_bf_zone[izone][0];
              header.rdata[1] = geom_bf_zone[izone][1];
              header.rdata[2] = geom_bf_zone[izone][2];
              header.rdata[3] = geom_bf_zone[izone][3];
              header.rdata[4] = geom_bf_zone[izone][4];
              assert(geom_bf_zone[izone][0] > 0.0); // area should be non-zero
              header.rdata[5] = geom_bf_zone[izone][5]/geom_bf_zone[izone][0];
              header.rdata[6] = geom_bf_zone[izone][6]/geom_bf_zone[izone][0];
              header.rdata[7] = geom_bf_zone[izone][7]/geom_bf_zone[izone][0];
              MPI_File_write_at(fh,offset,&header,1,MPI_Header,MPI_STATUS_IGNORE);
              offset += header_size;
            }
            // everyone needs to know how many zones to update their offsets...
            ++nzone_final;
          }
          else {
            // this should be a periodic zone. Another possibility is it is a tiny zone
            // not present in the final mesh for some reason (e.g. discarded because of tolerances)
            // but present in the surface...
            assert(nbf_zone[izone] == 0);
            assert(noobf_zone[izone] == 0);
            if (write_surface) {
              assert(spobf_zone[izone] == 0);
              assert(sbobf_zone[izone] == 0);
            }
            // this could happen for a zone that is FF
            //assert(zoneVec[izone].isPeriodic());
          }
        }

        if (mpi_rank == 0) delete[] geom_bf_zone;

        //offset += header_size*nzone_final;

        // now sync everyone else at this offset...

        MPI_Bcast(&offset,1,MPI_OFFSET_DATATYPE,0,mpi_comm);

      }

      // =================================
      // cvobf_global
      // =================================

      if ( mpi_rank == 0 ) {
        cout << " > cvobf" << endl;
        Header header;
        sprintf(header.name,"%s","cvobf_global");
        header.id   = UGP_IO_CVOBF_INT8;
        header.skip = header_size + nbf_global*sizeof(int8);
        ByteSwap::setLswMswPairForInt8(header.idata+0,nbf_global);
        MPI_File_write_at(fh,offset,&header,1,MPI_Header,MPI_STATUS_IGNORE);
      }

      offset += header_size;

      int8 my_disp_local = 0;
      for (int izone = 0; izone < nzone; ++izone) {
        if (nbf_zone[izone] > 0) {

          /*
             {
          // check...
          assert(my_nbf_zone[izone] == nbf_zone[izone]);
          for (int ibf = my_disp_local; ibf < my_disp_local+nbf_zone[izone]; ++ibf) {
          cout << "zone: " << izone << " ibf: " << ibf << " cvobf_global[ibf]: " << cvobf_global[ibf] << endl;
          }
          getchar();
          }
          */

          writeChunkedData<int8>( fh,offset+my_nbf_disp[izone]*sizeof(int8),cvobf_global+my_disp_local,my_nbf_zone[izone],mpi_comm);
          offset += nbf_zone[izone]*sizeof(int8);
          my_disp_local += my_nbf_zone[izone];
        }
      }
      assert(my_disp_local == nbf);

      // =================================
      // noobf_i/v global
      // =================================

      if ( mpi_rank == 0 ) {
        cout << " > noobf_i/v " << endl;
        Header header;
        sprintf(header.name,"noobf_i_and_v_int8");
        header.id   = UGP_IO_NOOBF_I_AND_V_INT8;
        header.skip = header_size + nbf_global*int_size + noobf_global*int8_size;
        ByteSwap::setLswMswPairForInt8(header.idata+0,nbf_global);
        ByteSwap::setLswMswPairForInt8(header.idata+2,noobf_global);
        MPI_File_write_at(fh,offset,&header,1,MPI_Header,MPI_STATUS_IGNORE);
      }

      offset += header_size;

      my_disp_local = 0;
      for (int izone = 0; izone < nzone; ++izone) {
        if (nbf_zone[izone] > 0) {
          writeChunkedData<int>(fh,offset+my_nbf_disp[izone]*sizeof(int),noobf_count+my_disp_local,my_nbf_zone[izone],mpi_comm);
          offset += nbf_zone[izone]*sizeof(int);
          my_disp_local += my_nbf_zone[izone];
        }
      }
      assert(my_disp_local == nbf);
      MiscUtils::dumpRange(noobf_count,nbf,"noobf_count");

      my_disp_local = 0;
      for (int izone = 0; izone < nzone; ++izone) {
        if (nbf_zone[izone] > 0) {
          writeChunkedData<int8>(fh,offset+my_noobf_disp[izone]*sizeof(int8),noobf_v_global+my_disp_local,my_noobf_zone[izone],mpi_comm);
          offset += noobf_zone[izone]*sizeof(int8);
          my_disp_local += my_noobf_zone[izone];
        }
      }

      if (write_surface) {

        // =================================
        // spobf_i/v/wgt global
        // =================================

        if ( mpi_rank == 0 ) {
          cout << " > spobf_i/v/wgt " << endl;
          Header header;
          sprintf(header.name,"spobf_i_v_wgt");
          header.id   = UGP_IO_SPOBF_I_V_WGT;
          header.skip = header_size + nbf_global*int_size + spobf_global*(uint8_size+double_size);
          header.ui8data[0] = nbf_global;
          header.ui8data[1] = spobf_global;
          MPI_File_write_at(fh,offset,&header,1,MPI_Header,MPI_STATUS_IGNORE);
        }

        offset += header_size;

        my_disp_local = 0;
        for (int izone = 0; izone < nzone; ++izone) {
          if (nbf_zone[izone] > 0) {
            writeChunkedData<int>(fh,offset+my_nbf_disp[izone]*sizeof(int),spobf_count+my_disp_local,my_nbf_zone[izone],mpi_comm);
            offset += nbf_zone[izone]*sizeof(int);
            my_disp_local += my_nbf_zone[izone];
          }
        }

        my_disp_local = 0;
        for (int izone = 0; izone < nzone; ++izone) {
          if (nbf_zone[izone] > 0) {
            writeChunkedData<uint8>(fh,offset+my_spobf_disp[izone]*sizeof(uint8),&(spobf_v[0])+my_disp_local,my_spobf_zone[izone],mpi_comm);
            offset += spobf_zone[izone]*sizeof(uint8);
            my_disp_local += my_spobf_zone[izone];
          }
        }

        my_disp_local = 0;
        for (int izone = 0; izone < nzone; ++izone) {
          if (nbf_zone[izone] > 0) {
            writeChunkedData<double>(fh,offset+my_spobf_disp[izone]*sizeof(double),&(spobf_wgt[0])+my_disp_local,my_spobf_zone[izone],mpi_comm);
            offset += spobf_zone[izone]*sizeof(double);
            my_disp_local += my_spobf_zone[izone];
          }
        }

        delete[] my_spobf_disp;
        delete[] spobf_zone;

        // =================================
        // sbobf_i/v global :: TODO: change to stobf ("sb" confusing, bits are always implied)
        // =================================

        if ( mpi_rank == 0 ) {
          cout << " > sbobf_i/v " << endl;
          Header header;
          sprintf(header.name,"sbobf_i_and_v_int");
          header.id   = UGP_IO_SBOBF_I_AND_V;
          header.skip = header_size + nbf_global*int_size + sbobf_global*uint8_size;
          ByteSwap::setLswMswPairForInt8(header.idata+0,nbf_global);
          ByteSwap::setLswMswPairForInt8(header.idata+2,sbobf_global);
          MPI_File_write_at(fh,offset,&header,1,MPI_Header,MPI_STATUS_IGNORE);
        }

        offset += header_size;

        my_disp_local = 0;
        for (int izone = 0; izone < nzone; ++izone) {
          if (nbf_zone[izone] > 0) {
            writeChunkedData<int>(fh,offset+my_nbf_disp[izone]*sizeof(int),sbobf_count+my_disp_local,my_nbf_zone[izone],mpi_comm);
            offset += nbf_zone[izone]*sizeof(int);
            my_disp_local += my_nbf_zone[izone];
          }
        }

        my_disp_local = 0;
        for (int izone = 0; izone < nzone; ++izone) {
          if (nbf_zone[izone] > 0) {
            writeChunkedData<uint8>(fh,offset+my_sbobf_disp[izone]*sizeof(uint8),&(psbobf_v[0])+my_disp_local,my_sbobf_zone[izone],mpi_comm);
            offset += sbobf_zone[izone]*sizeof(uint8);
            my_disp_local += my_sbobf_zone[izone];
          }
        }

        delete[] my_sbobf_disp;
        delete[] sbobf_zone;

      }

      // =================================
      // cvofa_global
      // =================================

      if ( mpi_rank == 0 ) {
        cout << " > cvofa (int+periodic)" << endl;
        Header header;
        sprintf(header.name,"cvofa_global");
        header.id    = UGP_IO_CVOFA_INT8;
        header.skip  = header_size + nfa_global*int8_size*2;
        ByteSwap::setLswMswPairForInt8(header.idata+0,nfa_global);
        header.idata[2] = 2;
        MPI_File_write_at(fh,offset,&header,1,MPI_Header,MPI_STATUS_IGNORE);
      }

      offset += header_size;

      my_disp_local = 0;
      for (int izone = 0; izone < 27; ++izone) {
        if (nfa_zone[izone] > 0) {
          writeChunkedData<int8>( fh,offset+my_nfa_disp[izone]*sizeof(int8)*2,(int8*)(cvofa_global+my_disp_local),my_nfa_zone[izone]*2,mpi_comm);
          offset += nfa_zone[izone]*sizeof(int8)*2;
          my_disp_local += my_nfa_zone[izone];
        }
      }

      // =================================
      // noofa_i/v global
      // =================================

      if ( mpi_rank == 0 ) {
        cout << " > noofa_i/v " << endl;
        Header header;
        sprintf(header.name,"noofa_i_and_v_int8");
        header.id   = UGP_IO_NOOFA_I_AND_V_INT8;
        header.skip = header_size + nfa_global*int_size + noofa_global*int8_size;
        ByteSwap::setLswMswPairForInt8(header.idata+0,nfa_global);
        ByteSwap::setLswMswPairForInt8(header.idata+2,noofa_global);
        MPI_File_write_at(fh,offset,&header,1,MPI_Header,MPI_STATUS_IGNORE);
      }

      offset += header_size;

      my_disp_local = 0;
      for (int izone = 0; izone < 27; ++izone) {
        if (nfa_zone[izone] > 0) {
          writeChunkedData<int>(  fh,offset+my_nfa_disp[izone]*sizeof(int),noofa_count+my_disp_local,my_nfa_zone[izone],mpi_comm);
          offset += nfa_zone[izone]*sizeof(int);
          my_disp_local += my_nfa_zone[izone];
        }
      }
      assert(my_disp_local == nfa);
      MiscUtils::dumpRange(noofa_count,nfa,"noofa_count");

      my_disp_local = 0;
      for (int izone = 0; izone < 27; ++izone) {
        if (nfa_zone[izone] > 0) {
          writeChunkedData<int8>( fh,offset+my_noofa_disp[izone]*sizeof(int8),noofa_v_global+my_disp_local,my_noofa_zone[izone],mpi_comm);
          offset += noofa_zone[izone]*sizeof(int8);
          my_disp_local += my_noofa_zone[izone];
        }
      }

      //===================================
      // periodic transforms
      //===================================

      if (PeriodicData::isPeriodic()) {
        assert(pbi_no);
        // periodicity is written differently depending on version.
        if (io_version == 4) {
          // for v4 io, the
          vector<int> bitsVec;
          if (PeriodicData::checkPeriodicBitPair(2)) {
            // triple periodicity...
            bitsVec.push_back((1<<0));
            bitsVec.push_back((1<<1));
            bitsVec.push_back((1<<2));
            bitsVec.push_back((1<<3));
            bitsVec.push_back((1<<4));
            bitsVec.push_back((1<<5));
            bitsVec.push_back((1<<0)|(1<<2));
            bitsVec.push_back((1<<0)|(1<<3));
            bitsVec.push_back((1<<0)|(1<<4));
            bitsVec.push_back((1<<0)|(1<<5));
            bitsVec.push_back((1<<1)|(1<<2));
            bitsVec.push_back((1<<1)|(1<<3));
            bitsVec.push_back((1<<1)|(1<<4));
            bitsVec.push_back((1<<1)|(1<<5));
            bitsVec.push_back((1<<2)|(1<<4));
            bitsVec.push_back((1<<2)|(1<<5));
            bitsVec.push_back((1<<3)|(1<<4));
            bitsVec.push_back((1<<3)|(1<<5));
            bitsVec.push_back((1<<0)|(1<<2)|(1<<4));
            bitsVec.push_back((1<<0)|(1<<2)|(1<<5));
            bitsVec.push_back((1<<0)|(1<<3)|(1<<4));
            bitsVec.push_back((1<<0)|(1<<3)|(1<<5));
            bitsVec.push_back((1<<1)|(1<<2)|(1<<4));
            bitsVec.push_back((1<<1)|(1<<2)|(1<<5));
            bitsVec.push_back((1<<1)|(1<<3)|(1<<4));
            bitsVec.push_back((1<<1)|(1<<3)|(1<<5));
            // that's 26 entries!!
          }
          else if (PeriodicData::checkPeriodicBitPair(1)) {
            // double periodicity...
            bitsVec.push_back((1<<0));
            bitsVec.push_back((1<<1));
            bitsVec.push_back((1<<2));
            bitsVec.push_back((1<<3));
            bitsVec.push_back((1<<2)|(1<<0));
            bitsVec.push_back((1<<3)|(1<<0));
            bitsVec.push_back((1<<2)|(1<<1));
            bitsVec.push_back((1<<3)|(1<<1));
          }
          else {
            // single periodicity...
            assert(PeriodicData::checkPeriodicBitPair(0));
            bitsVec.push_back((1<<0));
            bitsVec.push_back((1<<1));
          }
          // ----------------------------------------------------
          // now write a header for each possible transform with
          // the transform data stored as R[9] and t[3]...
          // ----------------------------------------------------
          for (int ii = 0; ii < bitsVec.size(); ++ii) {
            const int bits = bitsVec[ii];
            // cycle through bits pairs...
            // SB: "lucky the identity matrix is symmetric" -- inside joke re column-major format
            // of R -- hee hee
            double R[9];
            const bool has_R = PeriodicData::getPeriodicR(R,bits);
            double t[3];
            const bool has_t = PeriodicData::getPeriodicT(t,bits);
            // write on rank 0...
            if (mpi_rank == 0) { // root writes the transform record
              Header header;
              header.id     = UGP_IO_BIT_TRANSFORM_DATA;
              header.skip   = header_size;
              sprintf(header.name,"bit-%d",bits);
              if (has_R) {
                header.idata[0] = 1;
                for (int i = 0; i < 9; ++i)
                  header.rdata[i] = R[i];
              }
              else {
                header.idata[0] = -1; // R is NULL
              }
              if (has_t) {
                header.idata[1] =  1;
                for (int i = 0; i < 3 ; ++i)
                  header.rdata[i+9] = t[i];
              }
              else {
                header.idata[1] = -1; // t is NULL
              }

              header.idata[2] = bits;
              MPI_File_write_at(fh,offset,&header,1,MPI_Header,MPI_STATUS_IGNORE);
            }
            // everybody increment...
            offset += header_size;
          }
        }
        else if (io_version == 5) {
          if (mpi_rank == 0) { // root writes the transform record
            cout << " > periodic transform(s)" << endl;
            Header header;
            header.id     = UGP_IO_PERIODIC_TRANSFORM;
            header.skip   = header_size;
            sprintf(header.name,"periodic transforms");
            header.idata[0] = PeriodicData::periodicTransformVec.size();
            assert(PeriodicData::periodicTransformVec.size() <= 3); // allow 3 transforms max to be written
            for (int i = 0; i < PeriodicData::periodicTransformVec.size(); ++i) {
              header.idata[1+i] = PeriodicData::periodicTransformVec[i].getKind();
              FOR_J3 header.rdata[i*3+j] = PeriodicData::periodicTransformVec[i].getData(j);
            }
            MPI_File_write_at(fh,offset,&header,1,MPI_Header,MPI_STATUS_IGNORE);
          }
          // everybody increment...
          offset += header_size;

          // TODO XXXXX need to write pbi eventually

        }
        else {
          CERR("unsupported io version: " << io_version);
        }
      }
      else {
        assert(pbi_no == NULL);
      }

      //===================================
      // nodal records
      //===================================

      // pbi_pno (pno == periodic nodes)...

      if ( mpi_rank == 0 ) {
        cout << " > pbi_pno" << endl;
        Header header;
        sprintf(header.name,"pbi_pno");
        header.id     = UGP_IO_PBI_PNO;
        header.skip   = header_size + nno_p_global*int8_size;
        ByteSwap::setLswMswPairForInt8(header.idata+0,nno_p_global);
        MPI_File_write_at(fh,offset,&header,1,MPI_Header,MPI_STATUS_IGNORE);
      }

      int nno_p = 0;
      for (int i = 0; i < 7; ++i)
        nno_p += my_count_reduced[i];

      int8 my_nno_p = nno_p;
      int8 my_nno_p_disp;
      MPI_Scan(&my_nno_p,&my_nno_p_disp,1,MPI_INT8,MPI_SUM,mpi_comm);
      my_nno_p_disp -= my_nno_p;
      if (mpi_rank == mpi_size-1) assert(my_nno_p+my_nno_p_disp == nno_p_global);

      writeChunkedData<uint8>(fh,offset+header_size+my_nno_p_disp*int8_size,pbi_no,nno_p,mpi_comm);

      offset += header_size + nno_p_global*int8_size;

      // nodal coordinates...

      if ( mpi_rank == 0 ) {
        cout << " > x_no" << endl;
        Header header;
        sprintf(header.name,"x_no");
        header.id     = UGP_IO_X_NO;
        header.skip   = header_size + noora[mpi_size]*double_size*3;
        ByteSwap::setLswMswPairForInt8(header.idata+0,noora[mpi_size]);
        header.idata[2] = 3;
        MPI_File_write_at(fh,offset,&header,1,MPI_Header,MPI_STATUS_IGNORE);
      }

      writeChunkedData<double>(fh,offset+header_size+noora[mpi_rank]*double_size*3,(double*)x_no,nno*3,mpi_comm);

      offset += header_size + noora[mpi_size]*double_size*3;

      //==================================
      // cv based records
      //==================================

      if (mpi_rank == 0 ) { // one cv zone named "fluid"...
        Header header;
        sprintf(header.name,"fluid");
        header.id       = UGP_IO_CV_ZONE_HEADER;
        header.skip     = header_size;
        header.idata[0] = CV_ZONE_FLUID;
        header.idata[1] = 0;
        ByteSwap::setLswMswPairForInt8(header.idata+2,cvora[mpi_size]);
        int8 zero = 0;
        ByteSwap::setLswMswPairForInt8(header.idata+4,zero); // offset of this cv_zone in full range of cvs
        MPI_File_write_at(fh,offset,&header,1,MPI_Header,MPI_STATUS_IGNORE);
      }
      offset += header_size;

      // vol_cv...

      if ( mpi_rank == 0 ) {
        cout << " > vol_cv" << endl;
        Header header;
        sprintf(header.name,"vol_cv");
        header.id     = UGP_IO_CV_D1;
        header.skip   = header_size + cvora[mpi_size]*double_size;
        ByteSwap::setLswMswPairForInt8(header.idata+0,cvora[mpi_size]);
        MPI_File_write_at(fh,offset,&header,1,MPI_Header,MPI_STATUS_IGNORE);
      }

      writeChunkedData<double>(fh,offset+header_size+cvora[mpi_rank]*double_size,vol_cv,ncv,mpi_comm);

      offset += header_size + cvora[mpi_size]*double_size;

      // x_vv

      if ( mpi_rank == 0 ) {
        cout << " > x_vv" << endl;
        Header header;
        sprintf(header.name,"x_vv");
        header.id       = UGP_IO_CV_D2;
        header.skip     = header_size + cvora[mpi_size]*double_size*3;
        ByteSwap::setLswMswPairForInt8(header.idata+0,cvora[mpi_size]);
        header.idata[2] = 3;
        MPI_File_write_at(fh,offset,&header,1,MPI_Header,MPI_STATUS_IGNORE);
      }

      writeChunkedData<double>(fh,offset+header_size+cvora[mpi_rank]*double_size*3,(double*)x_vd,ncv*3,mpi_comm);

      offset += header_size + cvora[mpi_size]*double_size*3;

      if (io_version >= 5) {

        // r_vv...

        if ( mpi_rank == 0 ) {
          cout << " > r_vv" << endl;
          Header header;
          sprintf(header.name,"r_vv");
          header.id     = UGP_IO_CV_D1;
          header.skip   = header_size + cvora[mpi_size]*double_size;
          ByteSwap::setLswMswPairForInt8(header.idata+0,cvora[mpi_size]);
          MPI_File_write_at(fh,offset,&header,1,MPI_Header,MPI_STATUS_IGNORE);
        }

        writeChunkedData<double>(fh,offset+header_size+cvora[mpi_rank]*double_size,r_vv,ncv,mpi_comm);

        offset += header_size + cvora[mpi_size]*double_size;

        // ---------------------------------------------------------
        // also write ncv per rank and x_vv+r_vv bounding boxes...
        // ---------------------------------------------------------

        double * bbox_global = NULL;
        if ( mpi_rank == 0 ) {

          cout << " > vv_bbox" << endl;
          Header header;
          sprintf(header.name,"vv_bbox");
          header.id     = UGP_IO_VV_BBOX;
          header.skip   = header_size + mpi_size*int_size + mpi_size*double_size*6;
          header.idata[0] = mpi_size;
          MPI_File_write_at(fh,offset,&header,1,MPI_Header,MPI_STATUS_IGNORE);

          int * ncvora = new int[mpi_size];
          FOR_RANK {
            assert(cvora[rank+1]-cvora[rank] < TWO_BILLION);
            ncvora[rank] = cvora[rank+1]-cvora[rank];
          }
          MPI_File_write_at(fh,offset+header_size,ncvora,mpi_size,MPI_INT,MPI_STATUS_IGNORE);
          delete[] ncvora;

          bbox_global = new double[mpi_size*6];

        }

        offset += header_size + mpi_size*int_size;

        // bbox is xmin,ymin,zmin,xmax,ymax,zmax...

        double bbox[6] = { 1.0E+20, 1.0E+20, 1.0E+20, -1.0E+20, -1.0E+20, -1.0E+20 };
        FOR_ICV {
          FOR_I3 bbox[i] = min(bbox[i],x_vd[icv][i]-r_vv[icv]);
          FOR_I3 bbox[3+i] = max(bbox[3+i],x_vd[icv][i]+r_vv[icv]);
        }

        MPI_Gather(bbox,6,MPI_DOUBLE,bbox_global,6,MPI_DOUBLE,0,mpi_comm);

        if ( mpi_rank == 0 ) {

          MPI_File_write_at(fh,offset,bbox_global,mpi_size*6,MPI_DOUBLE,MPI_STATUS_IGNORE);
          delete[] bbox_global;

        }

        offset += mpi_size*double_size*6;

        // =========================================
        // cv centroids -- v5 and greater
        // =========================================

        // x_cv

        if ( mpi_rank == 0 ) {
          cout << " > x_cv" << endl;
          Header header;
          sprintf(header.name,"x_cv");
          header.id       = UGP_IO_CV_D2;
          header.skip     = header_size + cvora[mpi_size]*double_size*3;
          ByteSwap::setLswMswPairForInt8(header.idata+0,cvora[mpi_size]);
          header.idata[2] = 3;
          MPI_File_write_at(fh,offset,&header,1,MPI_Header,MPI_STATUS_IGNORE);
        }

        // HACK...
        /*
           {
           assert(ncv == cvora[mpi_size]);
           cout << "ncv: " << ncv << endl;
           for (int icv = 0; icv < ncv; ++icv) {
           cout << "icv: " << icv << " x_cv: " << COUT_VEC(x_cv[icv]) << endl;
           }
           getchar();
           }
           */

        writeChunkedData<double>(fh,offset+header_size+cvora[mpi_rank]*double_size*3,(double*)x_cv,ncv*3,mpi_comm);

        offset += header_size + cvora[mpi_size]*double_size*3;

        // =================================
        // boundary face-based geometry records...
        // =================================

        // n_bf...

        if ( mpi_rank == 0 ) {
          cout << " > n_bf" << endl;
          Header header;
          sprintf(header.name,"n_bf");
          header.id   = UGP_IO_BF_D2;
          header.skip = header_size + nbf_global*sizeof(double)*3;
          ByteSwap::setLswMswPairForInt8(header.idata+0,nbf_global);
          header.idata[2] = 3;
          MPI_File_write_at(fh,offset,&header,1,MPI_Header,MPI_STATUS_IGNORE);
        }

        offset += header_size;

        my_disp_local = 0;
        for (int izone = 0; izone < nzone; ++izone) {
          if (nbf_zone[izone] > 0) {
            writeChunkedData<double>( fh,offset+my_nbf_disp[izone]*sizeof(double)*3,(double*)(n_bf+my_disp_local),my_nbf_zone[izone]*3,mpi_comm);
            offset += nbf_zone[izone]*sizeof(double)*3;
            my_disp_local += my_nbf_zone[izone];
          }
        }
        assert(my_disp_local == nbf);

        // x_bf...

        if ( mpi_rank == 0 ) {
          cout << " > x_bf" << endl;
          Header header;
          sprintf(header.name,"x_bf");
          header.id   = UGP_IO_BF_D2;
          header.skip = header_size + nbf_global*sizeof(double)*3;
          ByteSwap::setLswMswPairForInt8(header.idata+0,nbf_global);
          header.idata[2] = 3;
          MPI_File_write_at(fh,offset,&header,1,MPI_Header,MPI_STATUS_IGNORE);
        }

        offset += header_size;

        my_disp_local = 0;
        for (int izone = 0; izone < nzone; ++izone) {
          if (nbf_zone[izone] > 0) {

            /*
               {
            // check...
            assert(my_nbf_zone[izone] == nbf_zone[izone]);
            for (int ibf = my_disp_local; ibf < my_disp_local+nbf_zone[izone]; ++ibf) {
            cout << "zone: " << izone << " ibf: " << ibf << " x_bf[ibf]: " << COUT_VEC(x_bf[ibf]) <<
            " cvobf_global[ibf]: " << cvobf_global[ibf] << " x_cv: " << COUT_VEC(x_cv[cvobf_global[ibf]]) << endl;
            }
            getchar();
            }
            */

            writeChunkedData<double>( fh,offset+my_nbf_disp[izone]*sizeof(double)*3,(double*)(x_bf+my_disp_local),my_nbf_zone[izone]*3,mpi_comm);
            offset += nbf_zone[izone]*sizeof(double)*3;
            my_disp_local += my_nbf_zone[izone];
          }
        }
        assert(my_disp_local == nbf);

        // area_bf...

        if ( mpi_rank == 0 ) {
          cout << " > area_bf" << endl;
          Header header;
          sprintf(header.name,"area_bf");
          header.id   = UGP_IO_BF_D1;
          header.skip = header_size + nbf_global*sizeof(double);
          ByteSwap::setLswMswPairForInt8(header.idata+0,nbf_global);
          MPI_File_write_at(fh,offset,&header,1,MPI_Header,MPI_STATUS_IGNORE);
        }

        offset += header_size;

        my_disp_local = 0;
        for (int izone = 0; izone < nzone; ++izone) {
          if (nbf_zone[izone] > 0) {
            writeChunkedData<double>( fh,offset+my_nbf_disp[izone]*sizeof(double),(double*)(area_bf+my_disp_local),my_nbf_zone[izone],mpi_comm);
            offset += nbf_zone[izone]*sizeof(double);
            my_disp_local += my_nbf_zone[izone];
          }
        }
        assert(my_disp_local == nbf);

        // area_over_delta_bf...

        if ( mpi_rank == 0 ) {
          cout << " > area_over_delta_bf" << endl;
          Header header;
          sprintf(header.name,"area_over_delta_bf");
          header.id   = UGP_IO_BF_D1;
          header.skip = header_size + nbf_global*sizeof(double);
          ByteSwap::setLswMswPairForInt8(header.idata+0,nbf_global);
          MPI_File_write_at(fh,offset,&header,1,MPI_Header,MPI_STATUS_IGNORE);
        }

        offset += header_size;

        my_disp_local = 0;
        for (int izone = 0; izone < nzone; ++izone) {
          if (nbf_zone[izone] > 0) {
            writeChunkedData<double>( fh,offset+my_nbf_disp[izone]*sizeof(double),(double*)(area_over_delta_bf+my_disp_local),my_nbf_zone[izone],mpi_comm);
            offset += nbf_zone[izone]*sizeof(double);
            my_disp_local += my_nbf_zone[izone];
          }
        }
        assert(my_disp_local == nbf);

        // Gij_bf...

        if ( mpi_rank == 0 ) {
          cout << " > Gij_bf" << endl;
          Header header;
          sprintf(header.name,"Gij_bf");
          header.id   = UGP_IO_BF_DN33;
          header.skip = header_size + nbf_global*sizeof(double)*9;
          ByteSwap::setLswMswPairForInt8(header.idata+0,nbf_global);
          header.idata[2] = 3;
          header.idata[3] = 3;
          MPI_File_write_at(fh,offset,&header,1,MPI_Header,MPI_STATUS_IGNORE);
        }

        offset += header_size;

        my_disp_local = 0;
        for (int izone = 0; izone < nzone; ++izone) {
          if (nbf_zone[izone] > 0) {
            writeChunkedData<double>( fh,offset+my_nbf_disp[izone]*sizeof(double)*9,(double*)(Gij_bf+my_disp_local),my_nbf_zone[izone]*9,mpi_comm);
            offset += nbf_zone[izone]*sizeof(double)*9;
            my_disp_local += my_nbf_zone[izone];
          }
        }
        assert(my_disp_local == nbf);

        // =================================
        // compact face-based geometry records...
        // =================================

        // n_fa...

        if ( mpi_rank == 0 ) {
          cout << " > n_fa (int+periodic) " << endl;
          Header header;
          sprintf(header.name,"n_fa");
          header.id    = UGP_IO_FA_D2;
          header.skip  = header_size + nfa_global*double_size*3;
          ByteSwap::setLswMswPairForInt8(header.idata+0,nfa_global);
          header.idata[2] = 3;
          MPI_File_write_at(fh,offset,&header,1,MPI_Header,MPI_STATUS_IGNORE);
        }

        offset += header_size;

        my_disp_local = 0;
        for (int izone = 0; izone < 27; ++izone) {
          if (nfa_zone[izone] > 0) {
            writeChunkedData<double>( fh,offset+my_nfa_disp[izone]*sizeof(double)*3,(double*)(n_fa+my_disp_local),my_nfa_zone[izone]*3,mpi_comm);
            offset += nfa_zone[izone]*sizeof(double)*3;
            my_disp_local += my_nfa_zone[izone];
          }
        }

        // x_fa...

        if ( mpi_rank == 0 ) {
          cout << " > x_fa (int+periodic) " << endl;
          Header header;
          sprintf(header.name,"x_fa");
          header.id    = UGP_IO_FA_D2;
          header.skip  = header_size + nfa_global*double_size*3;
          ByteSwap::setLswMswPairForInt8(header.idata+0,nfa_global);
          header.idata[2] = 3;
          MPI_File_write_at(fh,offset,&header,1,MPI_Header,MPI_STATUS_IGNORE);
        }

        offset += header_size;

        my_disp_local = 0;
        for (int izone = 0; izone < 27; ++izone) {
          if (nfa_zone[izone] > 0) {
            writeChunkedData<double>( fh,offset+my_nfa_disp[izone]*sizeof(double)*3,(double*)(x_fa+my_disp_local),my_nfa_zone[izone]*3,mpi_comm);
            offset += nfa_zone[izone]*sizeof(double)*3;
            my_disp_local += my_nfa_zone[izone];
          }
        }

      }

      // ====================
      // write the surface...
      // Note: for some reason on excalibur, the shared memory and the
      // MPI_FILE_write_at_all did not play well together, so the striped
      // part of the SHM surface is copied into a temporary buffer before
      // writing...
      // ====================

      if (!write_surface) {

        if ( mpi_rank == 0 ) {
          cout << " > skipping surface write (NO_SURFACE in WRITE_RESTART)" << endl;
        }

      }
      else {

        if ( mpi_rank == 0 ) {
          cout << " > surface (xp,spost,znost)" << endl;
          Header header;
          sprintf(header.name,"surface (xp,spost,znost)");
          header.id    = UGP_IO_SURFACE;
          header.skip  = header_size + getFlattenedNsp()*double_size*3 + getFlattenedNst()*int_size*(3+1);
          header.idata[0] = getFlattenedNsp();
          header.idata[1] = getFlattenedNst();
          MPI_File_write_at(fh,offset,&header,1,MPI_Header,MPI_STATUS_IGNORE);
        }

        offset += header_size;

        // xsp...

        for (int ipart = 0; ipart < partVec.size(); ++ipart) {
          SurfaceShm* surface = partVec[ipart]->surface;
          if (surface != NULL) {

            int * Xora = NULL;
            MiscUtils::calcThresholdDist(Xora,surface->nsp,mpi_size,DIST_THRESHOLD);

            // check...
            FOR_RANK assert(Xora[rank+1] >= Xora[rank]);
            assert(Xora[mpi_size] == surface->nsp);

            int nsp_check = surface->nsp;
            MPI_Bcast(&nsp_check,1,MPI_INT,0,mpi_comm);
            if (nsp_check != surface->nsp)
              cout << "rank: " << mpi_rank << " surface->nsp: " << surface->nsp << endl;

            int my_nsp = Xora[mpi_rank+1]-Xora[mpi_rank];

            assert(surface->xsp);
            double (*my_dbuf)[3] = new double[my_nsp][3];
            for (int isp = 0; isp < my_nsp; ++isp) {
              FOR_I3 my_dbuf[isp][i] = surface->xsp[Xora[mpi_rank]+isp][i];
            }

            writeChunkedData<double>(fh,offset+Xora[mpi_rank]*double_size*3,(double*)my_dbuf,my_nsp*3,mpi_comm);

            delete[] my_dbuf;
            delete[] Xora;

            offset += double_size*surface->nsp*3;
          }
        }

        // spost...

        for (int ipart = 0; ipart < partVec.size(); ++ipart) {
          SurfaceShm* surface = partVec[ipart]->surface;
          if (surface != NULL) {

            int * Xora = NULL;
            MiscUtils::calcThresholdDist(Xora,surface->nst,mpi_size,DIST_THRESHOLD);

            // check...
            FOR_RANK assert(Xora[rank+1] >= Xora[rank]);
            assert(Xora[mpi_size] == surface->nst);

            int my_nst = Xora[mpi_rank+1]-Xora[mpi_rank];

            assert(surface->spost);
            int (*my_ibuf)[3] = new int[my_nst][3];
            for (int ist = 0; ist < my_nst; ++ist) {
              FOR_I3 my_ibuf[ist][i] = getFlattenedSp(ipart,surface->spost[Xora[mpi_rank]+ist][i]);
            }

            writeChunkedData<int>(fh,offset+Xora[mpi_rank]*int_size*3,(int*)my_ibuf,my_nst*3,mpi_comm);

            delete[] my_ibuf;
            delete[] Xora;

            offset += int_size*surface->nst*3;
          }
        }


        // if periodic reorder zones putting periodic zones at the end (this is for efficient serial restart reading)

        int npz = 0;
        const int nz = zoneVec.size();
        int *new_zone = new int[nz];
        for (int izone = 0; izone < nz; ++izone)
          if (zoneVec[izone].getPeriodicBits())
            ++npz;

        int iz = 0;
        int izp = nz-npz;
        for (int izone = 0; izone < nz; ++izone) {
          if (zoneVec[izone].getPeriodicBits()) {
            //new_zone[izp++] = izone;
            new_zone[izone] = izp++;
          }
          else {
            //new_zone[iz++] = izone;
            new_zone[izone] = iz++;
          }
        }
        assert(iz == nz-npz);
        assert(izp == nz);

        // znost...

        for (int ipart = 0; ipart < partVec.size(); ++ipart) {
          SurfaceShm* surface = partVec[ipart]->surface;
          if (surface != NULL) {

            int * Xora = NULL;
            MiscUtils::calcThresholdDist(Xora,surface->nst,mpi_size,DIST_THRESHOLD);

            assert(surface->znost);
            assert(partVec[ipart]->znosz);

            // check...
            FOR_RANK assert(Xora[rank+1] >= Xora[rank]);
            assert(Xora[mpi_size] == surface->nst);

            int my_nst = Xora[mpi_rank+1]-Xora[mpi_rank];

            int *my_ibuf = new int[my_nst];
            assert(surface->znost);
            for (int ist = 0; ist < my_nst; ++ist) {
              my_ibuf[ist] = new_zone[partVec[ipart]->znosz[surface->znost[Xora[mpi_rank]+ist]]];
            }

            writeChunkedData<int>(fh,offset+Xora[mpi_rank]*int_size,my_ibuf,my_nst,mpi_comm);

            delete[] my_ibuf;
            delete[] Xora;

            offset += int_size*surface->nst;
          }
        }

        // include periodic data if present...
        const int npt = PeriodicData::periodicTransformVec.size();
        assert((npt >= 0)&&(npt <= 3));

        if (npt > 0) {

          // only write the pbi's that contain transform data...
          int npbi = 0;
          for (int ipart = 0; ipart < partVec.size(); ++ipart) {
            SurfaceShm* surface = partVec[ipart]->surface;
            if ((surface != NULL)&&(surface->pbi != NULL)) {
              for (int isp = 0; isp < surface->nsp; ++isp) {
                if (surface->pbi[isp] != uint8(isp))
                  ++npbi;
              }
            }
          }

          if ( mpi_rank == 0 ) {
            cout << " > surface periodic info" << endl;
            Header header;
            sprintf(header.name,"surface periodic info");
            header.id    = UGP_IO_SURFACE_PERIODIC_INFO;
            header.skip  = header_size + npz*(int_size+sizeof(uint2)) + npbi*(2*int_size+sizeof(uint2));
            header.idata[0] = npt;
            header.idata[1] = npz;
            header.idata[2] = npbi;

            // throw the transform info into the header...

            for (int ipt = 0; ipt < npt; ++ipt) {
              header.idata[3+ipt] = PeriodicData::periodicTransformVec[ipt].getKind();
              double temp[3];
              PeriodicData::periodicTransformVec[ipt].getData(temp);
              FOR_I3 header.rdata[ipt*3+i] = temp[i];
            }

            MPI_File_write_at(fh,offset,&header,1,MPI_Header,MPI_STATUS_IGNORE);
          }

          offset += header_size;

          // periodic zone indices and bits...
          // npz should be small... so lets just get rank 0 to write it...

          if (mpi_rank == 0) {
            int *periodic_zone_array = new int[npz];
            uint2 *periodic_zone_bits = new uint2[npz];
            int ipz = 0;
            for (int izone = 0; izone < nz; ++izone) {
              if (zoneVec[izone].getPeriodicBits()) {
                periodic_zone_array[ipz] = new_zone[izone]; // new zone index
                periodic_zone_bits[ipz]  = zoneVec[izone].getPeriodicBits();
                ++ipz;
              }
            }
            assert(ipz == npz);

            MPI_File_write_at(fh,offset,periodic_zone_array,npz,MPI_INT,MPI_STATUS_IGNORE);
            MPI_File_write_at(fh,offset+npz*int_size,periodic_zone_bits,npz,MPI_UINT2,MPI_STATUS_IGNORE);
            delete[] periodic_zone_array;
            delete[] periodic_zone_bits;

          }

          offset += npz*(int_size+sizeof(uint2));

          // pbi...

          int * Xora = NULL;
          MiscUtils::calcThresholdDist(Xora,npbi,mpi_size,DIST_THRESHOLD);
          FOR_RANK assert(Xora[rank+1] >= Xora[rank]);
          assert(Xora[mpi_size] == npbi);
          const int my_npbi = Xora[mpi_rank+1]-Xora[mpi_rank];
          int *isp_array = new int[my_npbi];
          int *isp_p_array = new int[my_npbi];
          uint2 *isp_bits_array = new uint2[my_npbi];
          int ipbi = 0;
          for (int ipart = 0; ipart < partVec.size(); ++ipart) {
            SurfaceShm* surface = partVec[ipart]->surface;
            if ((surface != NULL)&&(surface->pbi != NULL)) {
              for (int isp = 0; isp < surface->nsp; ++isp) {
                if (surface->pbi[isp] != uint8(isp)) {
                  // only pack your own...
                  if (ipbi >= Xora[mpi_rank] && ipbi < Xora[mpi_rank+1]) {
                    isp_array[ipbi] = getFlattenedSp(ipart,isp);
                    isp_p_array[ipbi] = getFlattenedSp(ipart,int(surface->pbi[isp]&MASK_52BITS));
                    isp_bits_array[ipbi] = uint2(surface->pbi[isp]>>52);
                  }
                  ++ipbi;
                }
              }
            }
          }
          assert(ipbi == npbi);
          writeChunkedData<int>(fh,offset+Xora[mpi_rank]*int_size,isp_array,my_npbi,mpi_comm);
          delete[] isp_array;
          writeChunkedData<int>(fh,offset+Xora[mpi_rank]*int_size+npbi*int_size,isp_p_array,my_npbi,mpi_comm);
          delete[] isp_p_array;
          writeChunkedData<uint2>(fh,offset+Xora[mpi_rank]*sizeof(uint2)+2*npbi*int_size,isp_bits_array,my_npbi,mpi_comm);
          delete[] isp_bits_array;
          delete[] Xora;

          offset += npbi*(2*int_size+sizeof(uint2));
        }
        delete[] new_zone;

      }

      //======================
      // and we're done
      //======================

      if ( mpi_rank == 0 ) {
        cout << " > EOF" << endl;
        Header header;
        header.id = UGP_IO_EOF;
        sprintf(header.name,"EOF");
        header.skip = header_size;
        MPI_File_write_at(fh,offset,&header,1,MPI_Header,MPI_STATUS_IGNORE);
      }

      offset += header_size;

      MPI_File_set_size(fh,offset);
      MPI_File_close(&fh);

    }
    else if (io_version == 3) {

      // ================================================================================
      // ================================================================================
      // ================================================================================
      // io version 3
      // ================================================================================
      // ================================================================================
      // ================================================================================

      MPI_File fh;
      MPI_Offset offset = 0;
      MPI_File_open(mpi_comm,fname,MPI_MODE_WRONLY|MPI_MODE_CREATE,MPI_INFO_NULL,&fh);

      if ( mpi_rank == 0 )  {
        int itmp[2] = {UGP_IO_MAGIC_NUMBER, io_version};
        cout << " > UGP_IO_VERSION: " << itmp[1] << endl;
        MPI_File_write(fh,itmp,2,MPI_INT,MPI_STATUS_IGNORE);
      }

      offset += int_size*2;

      //=========================
      // counts
      //=========================

      if ( mpi_rank == 0 ) {
        Header header;
        sprintf(header.name,"no_fa_bf_cv_counts");
        header.id     = UGP_IO_NO_FA_CV_COUNTS;
        header.skip   = header_size;
        ByteSwap::setLswMswPairForInt8(header.idata+0,noora[mpi_size]);
        ByteSwap::setLswMswPairForInt8(header.idata+2,nbf_global+nfa_global);
        ByteSwap::setLswMswPairForInt8(header.idata+4,cvora[mpi_size]);
        MPI_File_write(fh,&header,1,MPI_Header,MPI_STATUS_IGNORE);
      }

      offset += header_size;

      // ============================
      // fa_check...
      // ============================

      if ( mpi_rank == 0 ) {
        Header header;
        sprintf(header.name,"%s","fa_check");
        header.id   = UGP_IO_FA_CHECK;
        header.skip = header_size + int_size*(nbf_global+nfa_global);
        header.idata[0] = nbf_global+nfa_global;
        MPI_File_write(fh,&header,1,MPI_Header,MPI_STATUS_IGNORE);
      }

      offset += header_size;

      {

        int my_nfa_zone_max = 0;
        for (int izone = 0; izone < nzone; ++izone)
          my_nfa_zone_max = max(my_nfa_zone_max,(int)my_nbf_zone[izone]);
        for (int izone = 0; izone < 27; ++izone)
          my_nfa_zone_max = max(my_nfa_zone_max,(int)my_nfa_zone[izone]);

        int * fa_flag = new int[my_nfa_zone_max];

        int ifa_offset = 0;
        for (int izone = 0; izone < nzone; ++izone) {
          if (nbf_zone[izone] > 0) {
            for (int ifa = 0; ifa < my_nbf_zone[izone]; ++ifa)
              fa_flag[ifa] = ifa_offset+my_nbf_disp[izone]+ifa;
            MPI_File_write_at_all(fh,offset+my_nbf_disp[izone]*sizeof(int),
                fa_flag,my_nbf_zone[izone],
                MPI_INT,MPI_STATUS_IGNORE);
            offset += nbf_zone[izone]*sizeof(int);
            ifa_offset += nbf_zone[izone];
          }
        }
        for (int izone = 0; izone < 27; ++izone) {
          if (nfa_zone[izone] > 0) {
            for (int ifa = 0; ifa < my_nfa_zone[izone]; ++ifa)
              fa_flag[ifa] = ifa_offset+my_nfa_disp[izone]+ifa;
            MPI_File_write_at_all(fh,offset+my_nfa_disp[izone]*sizeof(int),
                fa_flag,my_nfa_zone[izone],
                MPI_INT,MPI_STATUS_IGNORE);
            offset += nfa_zone[izone]*sizeof(int);
            ifa_offset += nfa_zone[izone];
          }
        }

        delete[] fa_flag;

      }

      //=========================
      // face zone headers:
      // boundary zones first...
      //=========================

      double (*geom_bf_zone)[8] = NULL;
      if (mpi_rank == 0) geom_bf_zone = new double[nzone][8];
      MPI_Reduce(my_geom_bf_zone,geom_bf_zone,nzone*8,MPI_DOUBLE,MPI_SUM,0,mpi_comm);

      int nzone_final = 0;
      //int8 ibf_end = 0;
      //int8 noobf_end = 0;
      for (int izone = 0; izone < nzone; ++izone) {
        if (nbf_zone[izone] > 0) {
          assert(noobf_zone[izone] > 0);
          //const int8 ibf_begin = ibf_end;
          //ibf_end += nbf_zone[izone];
          //const int8 noobf_begin = noobf_end; noobf_end += noobf_zone[izone];
          // this zone needs to be written...
          if ( mpi_rank == 0 ) {
            Header header;
            header.id       = UGP_IO_FA_ZONE_HEADER;
            header.skip     = header_size;
            header.idata[0] = FA_ZONE_BOUNDARY;
            sprintf(header.name,"%s",zoneVec[izone].getName().c_str());
            cout << " > FA_ZONE \"" << header.name << "\" nfa: " << nbf_zone[izone] << endl;
            header.idata[1] = nzone_final;
            ByteSwap::setLswMswPairForInt8(header.idata+2,nbf_zone[izone]);
            // this data is not required/used in v4 restarts, but allows lesinfo to report
            // geom for v4 restart files...
            // recall...
            // 0: sum(area_bf)
            // 1: sum(area_over_delta_bf)
            // 2: sum(n_bf[0])
            // 3: sum(n_bf[1])
            // 4: sum(n_bf[2])
            // 5: x_bf[0] avg (eventually)
            // 6: x_bf[1] avg (eventually)
            // 7: x_bf[2] avg (eventually)
            header.rdata[0] = geom_bf_zone[izone][0];
            header.rdata[1] = geom_bf_zone[izone][1];
            header.rdata[2] = geom_bf_zone[izone][2];
            header.rdata[3] = geom_bf_zone[izone][3];
            header.rdata[4] = geom_bf_zone[izone][4];
            assert(geom_bf_zone[izone][0] > 0.0); // area should be non-zero
            header.rdata[5] = geom_bf_zone[izone][5]/geom_bf_zone[izone][0];
            header.rdata[6] = geom_bf_zone[izone][6]/geom_bf_zone[izone][0];
            header.rdata[7] = geom_bf_zone[izone][7]/geom_bf_zone[izone][0];
            // and write...
            MPI_File_write_at(fh,offset,&header,1,MPI_Header,MPI_STATUS_IGNORE);
          }
          ++nzone_final;
          offset += header_size;
        }
        else {
          assert(nbf_zone[izone] == 0);
          assert(noobf_zone[izone] == 0);
        }
      }

      if (mpi_rank == 0) delete[] geom_bf_zone;

      // ============================
      // periodic zones...
      // NOT supported for io version 3
      // ============================

      for (int izone = 0; izone < 26; ++izone) {
        if (nfa_zone[izone] > 0) {
          assert(0);
        }
      }

      // ============================
      // default internal...
      // ============================

      if (nfa_zone[26] > 0) {

        if ( mpi_rank == 0 ) {
          Header header;
          header.id       = UGP_IO_FA_ZONE_HEADER;
          header.skip     = header_size;
          header.idata[0] = FA_ZONE_INTERNAL;
          sprintf(header.name,"default-internal"); //.. protected name
          cout << " > FA_ZONE \"" << header.name << "\" nfa: " << nfa_zone[26] << endl;
          header.idata[1] = nzone_final;
          ByteSwap::setLswMswPairForInt8(header.idata+2,nfa_zone[26]);
          MPI_File_write_at(fh,offset,&header,1,MPI_Header,MPI_STATUS_IGNORE);
        }
        ++nzone_final;
        offset += header_size;

      }

      // ============================
      // fa_zone...
      // ============================

      if ( mpi_rank == 0 ) {
        Header header;
        sprintf(header.name,"%s","fa_zone");
        header.id   = UGP_IO_FA_ZONE;
        header.skip = header_size + int_size*(nbf_global+nfa_global);
        header.idata[0] = nbf_global+nfa_global;
        MPI_File_write_at(fh,offset,&header,1,MPI_Header,MPI_STATUS_IGNORE);
      }

      offset += header_size;

      {

        int my_nfa_zone_max = 0;
        for (int izone = 0; izone < nzone; ++izone)
          my_nfa_zone_max = max(my_nfa_zone_max,(int)my_nbf_zone[izone]);
        for (int izone = 0; izone < 27; ++izone)
          my_nfa_zone_max = max(my_nfa_zone_max,(int)my_nfa_zone[izone]);

        int * fa_zone = new int[my_nfa_zone_max];

        int izone_final = 0;
        for (int izone = 0; izone < nzone; ++izone) {
          if (nbf_zone[izone] > 0) {
            for (int ifa = 0; ifa < my_nbf_zone[izone]; ++ifa)
              fa_zone[ifa] = izone_final;
            MPI_File_write_at_all(fh,offset+my_nbf_disp[izone]*sizeof(int),
                fa_zone,my_nbf_zone[izone],
                MPI_INT,MPI_STATUS_IGNORE);
            offset += nbf_zone[izone]*sizeof(int);
            ++izone_final;
          }
        }
        for (int izone = 0; izone < 27; ++izone) {
          if (nfa_zone[izone] > 0) {
            for (int ifa = 0; ifa < my_nfa_zone[izone]; ++ifa)
              fa_zone[ifa] = izone_final;
            MPI_File_write_at_all(fh,offset+my_nfa_disp[izone]*sizeof(int),
                fa_zone,my_nfa_zone[izone],
                MPI_INT,MPI_STATUS_IGNORE);
            offset += nfa_zone[izone]*sizeof(int);
            ++izone_final;
          }
        }
        assert(izone_final == nzone_final);

        delete[] fa_zone;

      }

      //=========================
      // cv zone header:
      // - just 1 fluid zone
      //=========================

      if (mpi_rank == 0 ) { // one cv zone named "fluid"...
        Header header;
        sprintf(header.name,"fluid");
        header.id       = UGP_IO_CV_ZONE_HEADER;
        header.skip     = header_size;
        header.idata[0] = CV_ZONE_FLUID;
        header.idata[1] = 0;
        ByteSwap::setLswMswPairForInt8(header.idata+2,cvora[mpi_size]);
        MPI_File_write_at(fh,offset,&header,1,MPI_Header,MPI_STATUS_IGNORE);
      }

      offset += header_size;

      if ( mpi_rank == 0 ) {
        cout << " > cv_zone... " << endl;
        Header header;
        sprintf(header.name,"cv_zone");
        header.id     = UGP_IO_CV_ZONE;
        header.skip   = header_size + cvora[mpi_size]*int_size;
        header.idata[0] = cvora[mpi_size];
        MPI_File_write_at(fh,offset,&header,1,MPI_Header,MPI_STATUS_IGNORE);
      }

      {
        // cv_zone == 0
        int * cv_zone = new int[ncv];

        for (int icv = 0; icv < ncv; ++icv)
          cv_zone[icv] = 0;

        MPI_File_write_at_all(fh,offset+header_size+cvora[mpi_rank]*int_size,
            cv_zone,ncv,MPI_INT,MPI_STATUS_IGNORE);

        delete[] cv_zone;
      }

      offset += header_size + cvora[mpi_size]*int_size;

      // ============================
      // noofa_i/v...
      // ============================

      if ( mpi_rank == 0 ) {
        Header header;
        sprintf(header.name,"NOOFA_I_AND_V");
        header.id = UGP_IO_NOOFA_I_AND_V;
        header.skip = header_size + int_size*(nbf_global+nfa_global) + int_size*(noobf_global+noofa_global);
        header.idata[0] = nbf_global+nfa_global;
        // noofa_count gets written as 2 ints...
        ByteSwap::setLswMswPairForInt8(header.idata+1,noobf_global+noofa_global);
        MPI_File_write_at(fh,offset,&header,1,MPI_Header,MPI_STATUS_IGNORE);
      }

      offset += header_size;

      int my_disp_local = 0;
      for (int izone = 0; izone < nzone; ++izone) {
        if (nbf_zone[izone] > 0) {
          MPI_File_write_at_all(fh,offset+my_nbf_disp[izone]*sizeof(int),
              noobf_count+my_disp_local,my_nbf_zone[izone],
              MPI_INT,MPI_STATUS_IGNORE);
          offset += nbf_zone[izone]*sizeof(int);
          my_disp_local += my_nbf_zone[izone];
        }
      }
      my_disp_local = 0;
      for (int izone = 0; izone < 27; ++izone) {
        if (nfa_zone[izone] > 0) {
          MPI_File_write_at_all(fh,offset+my_nfa_disp[izone]*sizeof(int),
              noofa_count+my_disp_local,my_nfa_zone[izone],
              MPI_INT,MPI_STATUS_IGNORE);
          offset += nfa_zone[izone]*sizeof(int);
          my_disp_local += my_nfa_zone[izone];
        }
      }

      {

        int * noofa_v = new int[max(noobf_s,noofa_s)];

        for (int nob = 0; nob < noobf_s; ++nob) {
          assert(noobf_v_global[nob] == int(noobf_v_global[nob]));
          noofa_v[nob] = noobf_v_global[nob];
        }
        my_disp_local = 0;
        for (int izone = 0; izone < nzone; ++izone) {
          if (nbf_zone[izone] > 0) {
            MPI_File_write_at_all(fh,offset+my_noobf_disp[izone]*sizeof(int),
                noofa_v+my_disp_local,my_noobf_zone[izone],
                MPI_INT,MPI_STATUS_IGNORE);
            offset += noobf_zone[izone]*sizeof(int);
            my_disp_local += my_noobf_zone[izone];
          }
        }

        for (int nof = 0; nof < noofa_s; ++nof) {
          assert(noofa_v_global[nof] == int(noofa_v_global[nof]));
          noofa_v[nof] = noofa_v_global[nof];
        }
        my_disp_local = 0;
        for (int izone = 0; izone < 27; ++izone) {
          if (nfa_zone[izone] > 0) {
            MPI_File_write_at_all(fh,offset+my_noofa_disp[izone]*sizeof(int),
                noofa_v+my_disp_local,my_noofa_zone[izone],
                MPI_INT,MPI_STATUS_IGNORE);
            offset += noofa_zone[izone]*sizeof(int);
            my_disp_local += my_noofa_zone[izone];
          }
        }

        delete[] noofa_v;

      }

      // ============================
      // cvofa...
      // ============================

      if ( mpi_rank == 0 ) {
        Header header;
        sprintf(header.name,"CVOFA");
        header.id = UGP_IO_CVOFA;
        header.skip = header_size + int_size*(nbf_global+nfa_global)*2;
        header.idata[0] = nbf_global+nfa_global;
        header.idata[1] = 2;
        MPI_File_write_at(fh,offset,&header,1,MPI_Header,MPI_STATUS_IGNORE);
      }

      offset += header_size;

      {
        int (*cvofa)[2] = new int[max(nbf,nfa)][2];

        for (int ibf = 0; ibf < nbf; ++ibf) {
          assert(cvobf_global[ibf] == int(cvobf_global[ibf]));
          cvofa[ibf][0] = (int)cvobf_global[ibf];
          cvofa[ibf][1] = -1;
        }

        my_disp_local = 0;
        for (int izone = 0; izone < nzone; ++izone) {
          if (nbf_zone[izone] > 0) {
            MPI_File_write_at_all(fh,offset+my_nbf_disp[izone]*sizeof(int)*2,
                cvofa+my_disp_local,my_nbf_zone[izone]*2,
                MPI_INT,MPI_STATUS_IGNORE);
            offset += nbf_zone[izone]*sizeof(int)*2;
            my_disp_local += my_nbf_zone[izone];
          }
        }
        assert(my_disp_local == nbf);

        // XXXXX fix this for periodic zones eventually...

        for (int ifa = 0; ifa < nfa; ++ifa) {
          assert(cvofa_global[ifa][0] == int(cvofa_global[ifa][0]));
          assert(cvofa_global[ifa][1] == int(cvofa_global[ifa][1]));
          cvofa[ifa][0] = (int)cvofa_global[ifa][0];
          cvofa[ifa][1] = (int)cvofa_global[ifa][1];
        }

        my_disp_local = 0;
        for (int izone = 0; izone < 27; ++izone) {
          if (nfa_zone[izone] > 0) {
            MPI_File_write_at_all(fh,offset+my_nfa_disp[izone]*sizeof(int)*2,
                cvofa+my_disp_local,my_nfa_zone[izone]*2,
                MPI_INT,MPI_STATUS_IGNORE);
            offset += nfa_zone[izone]*sizeof(int)*2;
            my_disp_local += my_nfa_zone[izone];
          }
        }
        assert(my_disp_local == nfa);

        delete[] cvofa;

      }

      //===================================
      // nodal records
      //===================================

      // nodal coordinates...
      // double (*x_no)[3] = new double[nno][3];

      if ( mpi_rank == 0 ) {
        cout << " > x_no... " << endl;
        Header header;
        sprintf(header.name,"X_NO");
        header.id     = UGP_IO_X_NO;
        header.skip   = header_size + noora[mpi_size]*double_size*3;
        header.idata[0] = noora[mpi_size];
        header.idata[1] = 3;
        MPI_File_write_at(fh,offset,&header,1,MPI_Header,MPI_STATUS_IGNORE);
      }

      MPI_File_write_at_all(fh,offset+header_size+noora[mpi_rank]*double_size*3,
          x_no,nno*3,MPI_DOUBLE,MPI_STATUS_IGNORE);

      offset += header_size + noora[mpi_size]*double_size*3;

      //==================================
      // cv based records
      //==================================

      {

        if ( mpi_rank == 0 ) {
          cout << " > vol_vv... " << endl;
          Header header;
          sprintf(header.name,"VOL_VV");
          header.id     = UGP_IO_CV_D1;
          header.skip   = header_size + cvora[mpi_size]*double_size;
          header.idata[0] = cvora[mpi_size];
          MPI_File_write_at(fh,offset,&header,1,MPI_Header,MPI_STATUS_IGNORE);
        }

        MPI_File_write_at_all(fh,offset+header_size+cvora[mpi_rank]*double_size,
            vol_cv,ncv,MPI_DOUBLE,MPI_STATUS_IGNORE);
      }

      offset += header_size + cvora[mpi_size]*double_size;

      // x_vv
      if ( mpi_rank == 0 ) {
        cout << " > x_vv" << endl;
        Header header;
        sprintf(header.name,"X_VV");
        header.id       = UGP_IO_CV_D2;
        header.skip     = header_size + cvora[mpi_size]*double_size*3;
        header.idata[0] = cvora[mpi_size];
        header.idata[1] = 3;
        MPI_File_write_at(fh,offset,&header,1,MPI_Header,MPI_STATUS_IGNORE);
      }

      MPI_File_write_at_all(fh,offset+header_size+cvora[mpi_rank]*double_size*3,
          x_vd,ncv*3,MPI_DOUBLE,MPI_STATUS_IGNORE);

      offset += header_size + cvora[mpi_size]*double_size*3;

      // ============================
      // EOF...
      // ============================

      if ( mpi_rank == 0 ) {
        cout << " > EOF" << endl;
        Header header;
        header.id = UGP_IO_EOF;
        sprintf(header.name,"EOF");
        header.skip = header_size;
        MPI_File_write_at(fh,offset,&header,1,MPI_Header,MPI_STATUS_IGNORE);
      }

      offset += header_size;

      MPI_File_set_size(fh,offset);
      MPI_File_close(&fh);

    }
    else {

      CERR("unsupported io_version: " << io_version << ". Versions supported are 4 (default) and 3.");

    }

    // if we successfully wrote the tmp file, rename/replace the file with it...

    if (mpi_rank == 0) {
      remove(filename.c_str());
      rename(tmp_filename.c_str(),filename.c_str());
    }


    // cleanup...

    if (pbi_no) delete[] pbi_no;

    delete[] n_fa;
    delete[] x_fa;

    delete[] n_bf;
    delete[] x_bf;
    delete[] area_over_delta_bf;
    delete[] area_bf;
    delete[] Gij_bf;

    delete[] r_vv;

    delete[] noora;
    delete[] cvora;

    delete[] noobf_v_global;
    delete[] noofa_v_global;

    delete[] x_no;

    delete[] noofa_zone;
    delete[] my_noofa_disp;
    delete[] nfa_zone;
    delete[] my_nfa_disp;

    delete[] noobf_zone;
    delete[] my_noobf_disp;
    delete[] nbf_zone;
    delete[] my_nbf_disp;

    delete[] my_geom_bf_zone;

    delete[] cvofa_global;
    delete[] noofa_i;

    delete[] my_noofa_zone;
    delete[] my_nfa_zone;

    delete[] cvobf_global;
    delete[] noobf_i;

    delete[] my_noobf_zone;
    delete[] my_nbf_zone;

    delete[] sbobf_i;
    psbobf_v.clear();
    delete[] my_sbobf_zone;

    delete[] spobf_i;
    spobf_v.clear();
    spobf_wgt.clear();
    delete[] my_spobf_zone;

    }

    void clearSmoothing() {

      if (nsmooth != 0) {

        assert(xp0_and_delta0);
        nsmooth = 0;

        // now reset points...

        for (int icv = 0; icv < ncv; ++icv) {
          FOR_I3 x_vd[icv][i] = xp0_and_delta0[icv][i];
          delta_vd[icv]       = xp0_and_delta0[icv][3];
          flag_cv[icv]        = (flag_cv[icv]&PART_MASK_BITS);
        }

        delete[] xp0_and_delta0; xp0_and_delta0 = NULL;

      }

      b_vp = false;

    }

    void writeAutoCompleteJson(FILE * fp) {
      fprintf(fp,",\n\"siData\": {");
      fprintf(fp,"\n  \"kbVersion\": \"%s\"",CTI::cti_docs_version);
      const string siData =
#include "stitch.siData"
        ;
      fprintf(fp,",\n  %s\n",siData.c_str());
      fprintf(fp,"}");
    }

    bool b_skip_image = false;
    string skip_image_name;

    class View {
      public:
        string group;
        string name;

        double view_xp[3];
        double view_np[3];
        double view_up[3];
        double view_width;

        double plane_xp[3];
        double plane_np[3];

        bool operator==(const View &other) const {
          if (group!=other.group)
            return false;
          if (name!=other.name)
            return false;
          if (view_width!=other.view_width)
            return false;
          FOR_I3{
            if (view_xp[i]!=other.view_xp[i])
              return false;
            if (view_np[i]!=other.view_np[i])
              return false;
            if (view_up[i]!=other.view_up[i])
              return false;
            if (plane_xp[i]!=other.plane_xp[i])
              return false;
            if (plane_np[i]!=other.plane_np[i])
              return false;
          }
          return true;
        }

    };

    bool sortViewByGroup(View a,View b) {
      if (a.group==b.group)
        return (a.name<b.name);
      else
        return (a.group<b.group);
    }

    vector<View> viewVec;

  void processShow(Param * param,const bool b_help) {

    // look for help...
    if (b_help) {
      WUI(INFO,
          "SHOW changes the hidden zone settings\n" <<
          "  example:\n" <<
          "    SHOW ALL\n"
          "for more detail see [$CWIKB:stitch_export]");
      return;
    }

    assert(b_show_all == false);
    assert(b_show_touching == false);
    assert(b_show_zones == false);
    
    int iarg = 0;
    while (iarg < param->size()) {
      string token = MiscUtils::toUpperCase(param->getString(iarg++));
      if (token == "ALL") {
	WUI(INFO,"SHOW ALL processed");
	WebUI::webUIOutput.ensureImage();
	b_show_all = true;
      }
      else if (token == "TOUCHING") {
	// in this case, the current zones will be passed in the WRITE_JSON request...
	WUI(INFO,"SHOW TOUCHING processed");
	WebUI::webUIOutput.ensureImage();
	b_show_touching = true;
      }
      else if ((token == "FAZONE")||(token == "ZONE")||(token == "ZONES")) {
	// in this case, the zones are in the next param... 
	for (int izn = 0; izn < PartData::zoneVec.size(); ++izn)
	  PartData::zoneVec[izn].flag = 0;
        vector<string> zoneNameVec;
        MiscUtils::splitCsv(zoneNameVec,param->getString(iarg++));
	int ierr = 0;
        for (int ii = 0; ii < zoneNameVec.size(); ++ii) {
          // allow for a wildcard...
          bool found = false;
          for (int izn = 0; izn < PartData::zoneVec.size(); ++izn) {
            if (strcmp_wildcard(PartData::zoneVec[izn].getName(),zoneNameVec[ii])) {
	      PartData::zoneVec[izn].flag = 1;
	      found = true;
            }
          }
          if (!found) {
	    // maybe it was a zone index...
	    int izn;
	    if (from_string<int>(izn,zoneNameVec[ii],std::dec)&&(izn >= 0)&&(izn < PartData::zoneVec.size())) {
	      PartData::zoneVec[izn].flag = 1;
	    }
	    else {
	      WUI(WARN,"no matching zone found for " << zoneNameVec[ii]);
	      ierr = -1;
	    }
	  }
	}
	if (ierr != 0) {
	  WUI(WARN,"fix syntax and try again");
	}
	else {
	  assert(show_hiddenSubzonesVec.empty());
	  for (int izn = 0; izn < PartData::zoneVec.size(); ++izn) {
	    if (zoneVec[izn].flag == 0) {
	      show_hiddenSubzonesVec.push_back(izn);
	    }
	  }
	  WUI(INFO,"SHOW ZONE processed");
	  WebUI::webUIOutput.ensureImage();
	  b_show_zones = true;
	}
      }
      else if (token == "WITHOUT") {
	token = MiscUtils::toUpperCase(param->getString(iarg++));
	if (token == "HCP_WINDOW") {
	  // show the zones that DO NOT have an associated hcp_window...
	  set<int> intSet;
	  hcp.getHcpWindowZoneIds(intSet);
	  for (int izn = 0; izn < PartData::zoneVec.size(); ++izn)
	    PartData::zoneVec[izn].flag = 0;
	  for (set<int>::iterator it = intSet.begin(); it != intSet.end(); ++it) {
	    const int izn = *it;
	    assert((izn >= 0)&&(izn < PartData::zoneVec.size()));
	    PartData::zoneVec[izn].flag = 1;
	  }
	  assert(show_hiddenSubzonesVec.empty());
	  for (int izn = 0; izn < PartData::zoneVec.size(); ++izn) {
	    if (zoneVec[izn].flag == 1) {
	      show_hiddenSubzonesVec.push_back(izn);
	    }
	  }
	  WUI(INFO,"SHOW WITHOUT HCP_WINDOW processed");
	  WebUI::webUIOutput.ensureImage();
	  b_show_zones = true;
	}
	else {
	  WUI(WARN,"unrecognized SHOW WITHOUT token: " << token);
	}
      }
      else if (token == "WITH") {
	token = MiscUtils::toUpperCase(param->getString(iarg++));
	if (token == "HCP_WINDOW") {
	  // show the zones that DO NOT have an associated hcp_window...
	  set<int> intSet;
	  hcp.getHcpWindowZoneIds(intSet);
	  for (int izn = 0; izn < PartData::zoneVec.size(); ++izn)
	    PartData::zoneVec[izn].flag = 0;
	  for (set<int>::iterator it = intSet.begin(); it != intSet.end(); ++it) {
	    const int izn = *it;
	    assert((izn >= 0)&&(izn < PartData::zoneVec.size()));
	    PartData::zoneVec[izn].flag = 1;
	  }
	  assert(show_hiddenSubzonesVec.empty());
	  for (int izn = 0; izn < PartData::zoneVec.size(); ++izn) {
	    if (zoneVec[izn].flag == 0) {
	      show_hiddenSubzonesVec.push_back(izn);
	    }
	  }
	  WUI(INFO,"SHOW WITH HCP_WINDOW processed");
	  WebUI::webUIOutput.ensureImage();
	  b_show_zones = true;
	}
	else {
	  WUI(WARN,"unrecognized SHOW WITH token: " << token);
	}
      }
      else {
	WUI(WARN,"unrecognized SHOW token: " << token);
      }
    }
    
  
  }
  
  void processWriteJson(Param * param,const bool b_help) {

      if (b_help) {
        helpWriteJson();
      }
      else {

        try {

          //static bool first = true;

          // send user some recommendations if they haven't set a background resolution...
          // TODO: in some cases, all the points are provided and this is not required:
          // e.g. combining parts where all points are set...
          /*
             if (!hcp.b_hcp_delta) {
             if (first) {
             WebUI::webUIOutput.add(UIMessage(WARN,"You need to set a background resolution by setting HCP_DELTA."));
             hcp.reportHcpDeltaHelp();
             first = false;
             }
             }
             */

          bool b_hide_subzones = false;
          string hiddenSubzones;

          int iarg = 0;
          string prefix;
          uint json_bits = 0;
          while (iarg < param->size()) {
            string token = MiscUtils::toUpperCase(param->getString(iarg++));
            if (token == "NAME") {
              prefix = param->getString(iarg++);
            }
            else if (token == "HIDE_SUBZONES") {
              b_hide_subzones = true;
              hiddenSubzones = param->getString(iarg++);
            }
            else if (token == "HIDE_ZONES") {
              if (mpi_rank == 0) cout << "Warning: HIDE_ZONES ignored by stitch" << endl;
              ++iarg;
            }
            else if (token == "BITS") {
              json_bits = uint(param->getInt(iarg++));
            }
            else {
              if (mpi_rank == 0) cout << "unrecognized WRITE_JSON token \"" << token << "\"" << endl;
              throw(1);
            }
          }
          if (prefix == "") {
            if (mpi_rank == 0) cout << "expective WRITE_JSON NAME <string>" << endl;
            throw(1);
          }

          // this is implemented as per StaticSolver writeJson...
          // may not need this additional guard: i.e. could use the "newImage" behavior
          // without any messages some day...
          if (WebUI::webUIOutput.message_flag) {
            // the processParam routine that may have pushed messages into this JSON
            // may have set the newImage to false. If this is the case, we can skip
            // the image writing, which should share the same prefix as this JSON...
            if (!WebUI::webUIOutput.newImage) {
              b_skip_image = true;
              skip_image_name = prefix;
            }
          }

          // compute the bounding box including the effects of hiding...
          double xbb[3],diag,bbmin[3],bbmax[3];
          if (!b_hide_subzones) {
            getBoundingBoxCenter(xbb);
            diag = 2.0*getBoundingBoxRmax();
            getBbox(bbmin,bbmax);
          }
          else {
            // for now use full bbox -- TODO: needs implementation like surfer
            getBoundingBoxCenter(xbb);
            diag = 2.0*getBoundingBoxRmax();
            getBbox(bbmin,bbmax);
            /*
               vector<int> intVec;
               MiscUtils::splitCsv(intVec,hiddenSubzones);
               getBboxForZoneids(bbmin,bbmax,intVec);
               diag = sqrt( (bbmax[0]-bbmin[0])*(bbmax[0]-bbmin[0]) +
               (bbmax[1]-bbmin[1])*(bbmax[1]-bbmin[1]) +
               (bbmax[2]-bbmin[2])*(bbmax[2]-bbmin[2]) );
               xbb[0] = 0.5*(bbmin[0]+bbmax[0]);
               xbb[1] = 0.5*(bbmin[1]+bbmax[1]);
               xbb[2] = 0.5*(bbmin[2]+bbmax[2]);
               */
          }
	  
	  // ------------------------------------------------------
	  // show/hide commands may result in changes to the 
	  // next wid/image and need to be partially processed in parallel...
	  // ------------------------------------------------------
	  if (b_show_touching) {
	    assert(b_hide_subzones);
	    vector<int> intVec;
	    MiscUtils::splitCsv(intVec,hiddenSubzones);
	    for (int izn = 0; izn < zoneVec.size(); ++izn)
	      zoneVec[izn].flag = 1;
	    // turn off anything in the intVec...
	    for (int ii = 0; ii < intVec.size(); ++ii) {
	      const int izn = intVec[ii];
	      assert((izn >= 0)&&(izn < zoneVec.size()));
	      zoneVec[izn].flag = 0;
	    }
	    for (int izn = 0; izn < zoneVec.size(); ++izn) {
	      if (zoneVec[izn].flag) {
		const int ipart = pzozn[izn].first; assert((ipart >= 0)&&(ipart < partVec.size()));
		const int isz = pzozn[izn].second; assert(partVec[ipart]->surface); assert((isz >= 0)&&(isz < partVec[ipart]->surface->zoneVec.size()));
		partVec[ipart]->surface->ensureZnozn();
		for (int zoz = partVec[ipart]->surface->znozn_i[isz]; zoz != partVec[ipart]->surface->znozn_i[isz+1]; ++zoz) {
		  const int isz_nbr = partVec[ipart]->surface->znozn_v[zoz];
		  // the index in zoneVec of this nbr should share the same offset as izn-isz...
		  zoneVec[isz_nbr-isz+izn].flag = 1;
		}
	      }
	    }
	    assert(show_hiddenSubzonesVec.empty());
	    for (int izn = 0; izn < zoneVec.size(); ++izn) {
	      if (zoneVec[izn].flag == 0) {
		show_hiddenSubzonesVec.push_back(izn);
	      }
	    }
	  }
	  
          // only rank 0 write the json...

          if (mpi_rank == 0) {

            const string filename = prefix + ".json";
            const string tmp_filename = MiscUtils::makeTmpPrefix(filename);

            FILE * fp = fopen(tmp_filename.c_str(),"w");
            assert(fp);

            bool b_surface_loaded = !partVec.empty();
            bool b_volumeData = b_surface_loaded;

            fprintf(fp,"{\n\"solver\":\"stitch\"");
            fprintf(fp,",\n\"b_surface_loaded\":%s",(b_surface_loaded?"true":"false"));
            if (b_surface_loaded) {
              fprintf(fp,",\n\"BoundingBox\":[%f,%f,%f,%f,%f,%f]",bbmin[0],bbmax[0],bbmin[1],bbmax[1],bbmin[2],bbmax[2]);
              fprintf(fp,",\n\"BoundingBoxCentroid\":[%f,%f,%f]",xbb[0],xbb[1],xbb[2]);
              fprintf(fp,",\n\"BoundingBoxDiagonal\":%f",diag);
            }
            else {
              fprintf(fp,",\n\"BoundingBox\":[%f,%f,%f,%f,%f,%f]",0.0,0.0,0.0,0.0,0.0,0.0);
              fprintf(fp,",\n\"BoundingBoxCentroid\":[%f,%f,%f]",0.0,0.0,0.0);
              fprintf(fp,",\n\"BoundingBoxDiagonal\":%f",0.0);
            }
            fprintf(fp,",\n\"LightLocation\":[0.039503,0,1.01189,-0.530497,-0.57,0.581893]");

            if (!WebUI::webUIOutput.empty()) {
              assert(WebUI::webUIOutput.message_flag);
              WebUI::webUIOutput.writeJson(fp);
              //WebUI::webUIOutput.clear(); // done below on all ranks...
            }

            // bit-based criterion
            if (json_bits & (1<<0)) {
              // siData for ui autocompletion
              writeAutoCompleteJson(fp);
            }

            // need this to set "nwindows" in app, so that we can unpack index into level and part/window
            fprintf(fp,",\n\"HcpWindowCount\":%d",(int)hcp.hcpWindowVec.size());

            fprintf(fp,",\n\"HcpWindowParams\":[\n");
            for (int i = 0; i < hcp.hcpWindowVec.size(); ++i) {
              if (i == 0) fprintf(fp," \"%s\"",hcp.hcpWindowVec[i].param_str.c_str());
              else fprintf(fp,",\n \"%s\"",hcp.hcpWindowVec[i].param_str.c_str());
            }
            fprintf(fp,"\n]");

            // DP/ME: what is filetype
            //fprintf(fp,",\n\"filetype\": \"%s\"", filetype.c_str());

            double (*zone_x_area)[4] = new double[zoneVec.size()][4];
            for (int i = 0; i < zoneVec.size(); ++i) FOR_J4 zone_x_area[i][j] = 0.0;

            // don't do this here: expensive. Do it in parallel in a function that doesn't redo it
            // every time this is called...
            /*
               for (int ist = 0; ist < surface->nst; ++ist) {
               const double n_st[3] = TRI_NORMAL_2(surface->xp[surface->spost[ist][0]],surface->xp[surface->spost[ist][1]],surface->xp[surface->spost[ist][2]]);
               const double n_st_mag = MAG(n_st);
               assert(n_st_mag > 0.0);
               FOR_I3  zone_x_area[surface->znost[ist]][i] += n_st_mag*(surface->xp[surface->spost[ist][0]][i]+
               surface->xp[surface->spost[ist][1]][i]+
               surface->xp[surface->spost[ist][2]][i]);
               zone_x_area[surface->znost[ist]][3] += n_st_mag;
               }
               for (int i = 0; i < surface->zoneVec.size(); ++i) {
               FOR_J3 zone_x_area[i][j] /= 3.0*zone_x_area[i][3];
               zone_x_area[i][3] *= 0.5;
               }
               */

            fprintf(fp,",\n\"zoneAreas\":[\n");
            for (int i = 0; i < zoneVec.size(); ++i) {
              if (i == 0) fprintf(fp," %f",zone_x_area[i][3]);
              else fprintf(fp,",\n %f",zone_x_area[i][3]);
            }
            fprintf(fp,"\n]");
            fprintf(fp,",\n\"zoneCentroids\":[\n");
            for (int i = 0; i < zoneVec.size(); ++i) {
              if (i == 0) fprintf(fp," [%f,%f,%f]",zone_x_area[i][0],zone_x_area[i][1],zone_x_area[i][2]);
              else fprintf(fp,",\n [%f,%f,%f]",zone_x_area[i][0],zone_x_area[i][1],zone_x_area[i][2]);
            }
            fprintf(fp,"\n]");
            delete[] zone_x_area;
            fprintf(fp,",\n\"zoneIds\":[\n");
            for (int i = 0; i < zoneVec.size(); ++i) {
              if (i == 0) fprintf(fp," %d",i);
              else fprintf(fp,",\n %d",i);
            }
            fprintf(fp,"\n]");
            fprintf(fp,",\n\"zoneNames\":[\n");
            for (int i = 0, limit = zoneVec.size(); i < limit; ++i) {
              if (i == 0) fprintf(fp," \"%s\"",zoneVec[i].getName().c_str());
              else fprintf(fp,",\n \"%s\"",zoneVec[i].getName().c_str());
            }
            fprintf(fp,"\n]");
            // for subzones, just use the zone index -- i.e. no subzoning
            // DP/ME is this correct?
            fprintf(fp,",\n\"zoneSubzones\":[\n");

            /*
               for (int i = 0, limit = ss.szozn_i.getLength(); i < limit; ++i) {
               if (i == 0) fprintf(fp," \"%d\"",ss.szozn_i[i]);
               else fprintf(fp,",\n \"%d\"",ss.szozn_i[i]);
               }
               fprintf(fp,"\n]");
               */

            for (int i = 0; i <= zoneVec.size(); ++i) {
              if (i == 0) fprintf(fp," \"%d\"",i);
              else fprintf(fp,",\n \"%d\"",i);
            }
            fprintf(fp,"\n]");
            /*
               if (ss.selectedSubzoneVec.size()) {
               fprintf(fp,",\n\"selectedSubzones\":[\n");
               for (int i = 0, limit = ss.selectedSubzoneVec.size(); i < limit; ++i) {
               if (i == 0) fprintf(fp," \"%d\"",ss.selectedSubzoneVec[i]);
               else fprintf(fp,",\n \"%d\"",ss.selectedSubzoneVec[i]);
               }
               fprintf(fp,"\n]");
               }
               */
            fprintf(fp,",\n\"zonePeriodicPairs\":[]");

            if (b_volumeData) {
              fprintf(fp,",\n\"volumeVars\":{\n");
              fprintf(fp," \"D1\":[\n");
              fprintf(fp,"  \"mesh\"");
              fprintf(fp,",\n  \"%s\"","level");
              fprintf(fp,"\n ],\n");
              fprintf(fp," \"D2\":[\n");
              fprintf(fp,"\n  ]");
              fprintf(fp,"\n }");
            }

	    if (b_show_all) {
	      // show all has been requested...
	      fprintf(fp,",\n\"hiddenSubzones\":[\n");
	      fprintf(fp,"\n]");
	    }
	    else if ((b_show_touching)||(b_show_zones)) {
	      // for the case of touching, we did the necessary processing in parallel above...
	      fprintf(fp,",\n\"hiddenSubzones\":[\n");
	      for (int i = 0, limit = show_hiddenSubzonesVec.size(); i < limit; ++i) {
		if (i == 0) {
		  fprintf(fp," \"%d\"",show_hiddenSubzonesVec[i]);
		}
		else {
		  fprintf(fp,",\n \"%d\"",show_hiddenSubzonesVec[i]);
		}
	      }
	      fprintf(fp,"\n]");
	    }
	    
            if (!viewVec.empty()) {
              //sort by group;
              sort(viewVec.begin(), viewVec.end(), sortViewByGroup);

              fprintf(fp,",\n\"animateViews\":{\n");

              string groupVal = "";
              for (int ii = 0, size = viewVec.size(); ii < size; ++ii) {
                if (viewVec[ii].group != groupVal){
                  groupVal = viewVec[ii].group;
                  fprintf(fp,"\"%s\":[\n",groupVal.c_str());
                }
                fprintf(fp,"  {\n");
                fprintf(fp,"    \"name\": \"%s\",\n",viewVec[ii].name.c_str());
                fprintf(fp,"    \"view\": {\n");
                fprintf(fp,"      \"target\": {\n");
                fprintf(fp,"        \"%d\": %f,\n",0,viewVec[ii].view_xp[0]);
                fprintf(fp,"        \"%d\": %f,\n",1,viewVec[ii].view_xp[1]);
                fprintf(fp,"        \"%d\": %f\n",2,viewVec[ii].view_xp[2]);
                fprintf(fp,"      },\n");
                fprintf(fp,"      \"np\": {\n");
                fprintf(fp,"        \"%d\": %f,\n",0,viewVec[ii].view_np[0]);
                fprintf(fp,"        \"%d\": %f,\n",1,viewVec[ii].view_np[1]);
                fprintf(fp,"        \"%d\": %f\n",2,viewVec[ii].view_np[2]);
                fprintf(fp,"      },\n");
                fprintf(fp,"      \"up\": {\n");
                fprintf(fp,"        \"%d\": %f,\n",0,viewVec[ii].view_up[0]);
                fprintf(fp,"        \"%d\": %f,\n",1,viewVec[ii].view_up[1]);
                fprintf(fp,"        \"%d\": %f\n",2,viewVec[ii].view_up[2]);
                fprintf(fp,"      },\n");
                fprintf(fp,"      \"width\": %f\n",viewVec[ii].view_width);
                fprintf(fp,"    },\n");
                fprintf(fp,"    \"plane\": [\n");
                fprintf(fp,"      \"%f\",\n",viewVec[ii].plane_xp[0]);
                fprintf(fp,"      \"%f\",\n",viewVec[ii].plane_xp[1]);
                fprintf(fp,"      \"%f\",\n",viewVec[ii].plane_xp[2]);
                fprintf(fp,"      \"%f\",\n",viewVec[ii].plane_np[0]);
                fprintf(fp,"      \"%f\",\n",viewVec[ii].plane_np[1]);
                fprintf(fp,"      \"%f\"\n",viewVec[ii].plane_np[2]);
                fprintf(fp,"    ]\n");
                if (ii<size-1) {
                  if (groupVal==viewVec[ii+1].group)
                    fprintf(fp,"  },\n");
                  else{
                    groupVal=viewVec[ii+1].group;
                    fprintf(fp,"  }\n],\n\"%s\":[\n",groupVal.c_str());
                  }
                }
                else
                  fprintf(fp,"  }\n]\n");
              }
              fprintf(fp,"}");
            }
            fprintf(fp,"\n}\n");

            fclose(fp);
            remove(filename.c_str());
            rename(tmp_filename.c_str(),filename.c_str());

          }

          // all need to clear...
          WebUI::webUIOutput.clear();


        }
        catch(int e) {
          WUI(WARN,"expecting WRITE_JSON NAME <string>");
          helpWriteJson();
        }

      }

    }

    void helpWriteJson() {
      WUI(INFO,
          "WRITE_JSON writes some solver-specific information (in JavaScript Object Notation) to be used by the app. Example:\n" <<
          "WRITE_JSON NAME my_filename # write JSON data to file my_filename.json");
    }

    class WeightedCentroid {
      public:
        double x[3],wgt;
        WeightedCentroid() {
          x[0] = 0.0;
          x[1] = 0.0;
          x[2] = 0.0;
          wgt = 0.0;
        }
    };

    void addSimpleTrisInternal(vector<SimpleTri>& triVec,const int ii,const double dx[3]) {

      vector<WeightedCentroid> wcVec;
      map<const uint8,int> internalWCMap;

      for (int cycle = 0; cycle <= 1; ++cycle) {

        for (int ied = 0; ied < vdReturnDataVec[ii].ned; ++ied) {
          FOR_I2 {
            const int ifa = vdReturnDataVec[ii].faoed[ied][i];
            if (ifa < 0) {
              // recall vdArray[icv] uses -1 indexing for faces...
              assert((-ifa-1 >= 0)&&(-ifa-1 < vdReturnDataVec[ii].nfa));
              // this is an active face - add it to the faceMap...
              if (vdReturnDataVec[ii].faceIsActive(-ifa-1)) {
                const int icv_nbr = vdReturnDataVec[ii].cvofa[-ifa-1];
                uint8 rbiHash;
                if (icv_nbr < ncv)
                  rbiHash = BitUtils::packRankBitsIndex(mpi_rank,0,icv_nbr);
                else
                  rbiHash = rbi_g[icv_nbr-ncv];
                map<const uint8,int>::iterator iter = internalWCMap.find(rbiHash);
                int iwc;
                if (iter == internalWCMap.end()) {
                  assert(cycle == 0);
                  internalWCMap[rbiHash] = iwc =  wcVec.size();
                  wcVec.push_back(WeightedCentroid());
                }
                else {
                  iwc = iter->second;
                }
                const double * const x0 = vdReturnDataVec[ii].x_no[vdReturnDataVec[ii].nooed[ied][i]];
                const double * const x1 = vdReturnDataVec[ii].x_no[vdReturnDataVec[ii].nooed[ied][1-i]];

                if (cycle == 0) {
                  const double length = DIST(x0,x1);
                  FOR_J3 wcVec[iwc].x[j] += length*(x0[j]+x1[j]);
                  wcVec[iwc].wgt += length;
                }
                else {
                  triVec.push_back(SimpleTri(x1,x0,wcVec[iwc].x,dx));
                }
              }
            }
          }
        }

        // on the first time through, normalize the centroid...

        if (cycle == 0) {
          for (int iwc = 0; iwc < wcVec.size(); ++iwc) {
            assert(wcVec[iwc].wgt > 0.0);
            FOR_I3 wcVec[iwc].x[i] /= 2.0*wcVec[iwc].wgt;
          }
        }
      }
    }

    void addSimpleTrisBoundary(vector<pair<SimpleTriWithData,int> >& triVec,const int ii,const double dx[3]) {

      map<const pair<pair<int,int>,int>,int> colorMap; // boundary face color with ((ipart,ist),bits) key
      buildBfColorMap(colorMap,ii);

      vector<WeightedCentroid> wcVec;
      vector<WeightedCentroid> ifa_normal;  // use container to store averaged normal of decimated faces
      map<const int,int> boundaryWCMap;

      for(int cycle = 0; cycle <= 2; ++cycle) {

        for (int ied = 0; ied < vdReturnDataVec[ii].ned; ++ied) {
          FOR_I2 {
            const int ifa = vdReturnDataVec[ii].faoed[ied][i];
            if (ifa >= 0) {
              // a positive face indicates an ibf reference to an ist/bit pair...
              assert(ifa < vdReturnDataVec[ii].nbf);
              const int ist = vdReturnDataVec[ii].spbobf[ifa][0];
              const int bits = (vdReturnDataVec[ii].spbobf[ifa][1]&MASK_6BITS);
              const int ipart = (vdReturnDataVec[ii].spbobf[ifa][1]>>6);
              assert((ipart >= 0)&&(ipart < partVec.size()));
              map<const pair<pair<int,int>,int>,int>::iterator iter =
                colorMap.find(pair<pair<int,int>,int>(pair<int,int>(ipart,ist),bits));
              assert(iter != colorMap.end());
              // only image this edge if its other edge is also valid and has a different unit normal...
              const int ifa_nbr = vdReturnDataVec[ii].faoed[ied][1-i];
              if (ifa_nbr >= 0) {
                assert(ifa_nbr < vdReturnDataVec[ii].nbf);
                const int ist_nbr = vdReturnDataVec[ii].spbobf[ifa_nbr][0];
                const int bits_nbr = (vdReturnDataVec[ii].spbobf[ifa_nbr][1]&MASK_6BITS);
                const int ipart_nbr = (vdReturnDataVec[ii].spbobf[ifa_nbr][1]>>6);
                assert((ipart_nbr >= 0)&&(ipart_nbr < partVec.size()));
                map<const pair<pair<int,int>,int>,int>::iterator iter_nbr =
                  colorMap.find(pair<pair<int,int>,int>(pair<int,int>(ipart_nbr,ist_nbr),bits_nbr));
                assert(iter_nbr != colorMap.end());
                if (iter->second == iter_nbr->second)
                  continue;
              }

              // this is a boundary face pointing to a surface ist...
              map<const int,int>::iterator iter_w = boundaryWCMap.find(iter->second);
              int iwc;
              if (iter_w == boundaryWCMap.end()) {
                assert(cycle == 0);
                boundaryWCMap[iter->second] = iwc = wcVec.size();
                wcVec.push_back(WeightedCentroid());
                ifa_normal.push_back(WeightedCentroid());
              }
              else {
                iwc = iter_w->second;
              }
              const double * const x0 = vdReturnDataVec[ii].x_no[vdReturnDataVec[ii].nooed[ied][i]];
              const double * const x1 = vdReturnDataVec[ii].x_no[vdReturnDataVec[ii].nooed[ied][1-i]];

              if (cycle == 0) {
                // on the first time through, we are building the face centroid...
                const double length = DIST(x0,x1);
                FOR_J3 wcVec[iwc].x[j] += length*(x0[j]+x1[j]);
                wcVec[iwc].wgt += length;
              }
              else if (cycle == 1) {
                double n_fa[3] = TRI_NORMAL_2(x0,x1,wcVec[iwc].x);
                FOR_J3 {
                  n_fa[j] *= 0.5/wcVec[iwc].wgt;
                  ifa_normal[iwc].x[j] += n_fa[j];
                }
              }
              else {
                triVec.push_back(pair<SimpleTriWithData,int>(SimpleTriWithData(x0,x1,wcVec[iwc].x,dx,ifa_normal[iwc].x[0],ifa_normal[iwc].x[1],ifa_normal[iwc].x[2]),
                      partVec[ipart]->znosz[partVec[ipart]->surface->znost[ist]]));
              }
            }
          }
        }

        if (cycle == 0) {
          // on the first time through, normalize the centroid...
          for (int iwc = 0, end=wcVec.size(); iwc<end; ++iwc) {
            assert(wcVec[iwc].wgt > 0.0);
            const double one_o_wgt = 1.0/wcVec[iwc].wgt;
            FOR_I3 wcVec[iwc].x[i] *= 0.5*one_o_wgt;
          }
        }
        else if (cycle == 1) {
          // second time through average the summed normals
          for (int iwc = 0, end=wcVec.size(); iwc<end; ++iwc) {
            assert(wcVec[iwc].wgt > 0.0);
            const double one_o_wgt = 1.0/wcVec[iwc].wgt;
            FOR_I3 ifa_normal[iwc].x[i] *= 0.5*one_o_wgt;
          }
        }
      }
    }

    void processWriteImage(Param * param,const bool b_help) {

      if (b_help) {
        helpWriteImage();
      }
      else {

        try {

          // here we use a 2-step initialization of the scene, where we
          // initialize a "wid" first using the wid.init(param,step) method
          // which will return non-zero if step%interval != 0, and we
          // can return quickly without further instantiation...

          // note that we hard-code step = 0 here but some day when this
          // is part of stitch lloyd iterations, this may be useful...

#ifdef WITH_CHRONO
          auto t0 = std::chrono::high_resolution_clock::now();
#endif

          int step = 0;
          WriteImageData wid;
          if (wid.init(param,step) != 0)
            return;

          // when interacting with the client, sometimes an image is not required,
          // even though one was requested.
          if (b_skip_image && wid.b_name && (skip_image_name == wid.name)) {
            // skip this immage...
            if (mpi_rank == 0) cout << " > skipping this image" << endl;
            b_skip_image = false;
            return;
          }
	  
	  // check for app-driven show/hide...
	  // NOTE: this approach to selection/show is probably less desirable than letting the
	  // app drive, because it puts stitch in a certain "state" that assumes we write a json
	  // and then write an image. There may however be necessary for certain complex 
	  // show requests that the app cannot figure out...
	  if (b_show_all) {
	    wid.b_hidesubzones = false;
            wid.hiddenSubzonesSet.clear();
	    b_show_all = false;
	  }
	  else if ((b_show_touching)||(b_show_zones)) {
	    wid.b_hidesubzones = true;
            wid.hiddenSubzonesSet.clear(); // TODO: why is this a set? slows down parsing
	    for (int ii = 0; ii < show_hiddenSubzonesVec.size(); ++ii) {
	      wid.hiddenSubzonesSet.insert(show_hiddenSubzonesVec[ii]);
	    }
	    show_hiddenSubzonesVec.clear();
	    b_show_touching = false;
	    b_show_zones = false;
	  }
	  
	  CtiScene * scene = new CtiScene(wid);

          // pass details about the model geometry to the scene->..
          // TODO: make sure the following are also parallel...

          scene->setRmax(getBoundingBoxRmax());
          double bbmin[3],bbmax[3];
          getBbox(bbmin,bbmax);
          scene->setBoundingBbox(bbmin[0],bbmax[0],bbmin[1],bbmax[1],bbmin[2],bbmax[2]);
          double center[3];
          getBoundingBoxCenter(center);
          scene->setCenter(center[0],center[1],center[2]);
          //scene->setMaxIndex(35);

          // and draw...

          for (int izone = 0; izone < zoneVec.size(); ++izone) scene->addZoneName(zoneVec[izone].getName());
          scene->convertHiddenZoneNamesToIndices();
          scene->initCanvas();

          // put this following into init...

          if (scene->hasGeomPlane()) {
            scene->setGeomPlaneAsDataPlane();
            if (scene->blankDataPlane())
              scene->addGeomPlaneAsBlankPlane(); // should be changed to back in scene/canvas
          }

          const double wtime0 = MPI_Wtime();
          double wtime_prev = wtime0;

          if (!(scene->hasVar()&&(scene->getVar() == "3d_mesh"))) {

            // need full surface (for now) to get properly masked canvas, similar to "points" viz in stitch...

            for (int ipart = 0; ipart < partVec.size(); ++ipart) {
              if (partVec[ipart]->surface) {
                // divide up surface tris evenly and contiguously amongst all ranks...
                int my_nst_avg = partVec[ipart]->surface->nst/mpi_size;
                if (partVec[ipart]->surface->nst%mpi_size) ++my_nst_avg;
                const int ist0 = min(partVec[ipart]->surface->nst,mpi_rank*my_nst_avg);
                const int ist1 = min(partVec[ipart]->surface->nst,(mpi_rank+1)*my_nst_avg);
                assert(ist1-ist0 <= my_nst_avg);
                // for the zone information, we provide a local array that maps the surface tri zone to the
                // solver zone using the part's znosz (zone-of-surface-zone)...
                int * zone = new int[ist1-ist0];
                int (*spost)[3] = new int[ist1-ist0][3];
                int nst_render = 0;
                for (int ist = ist0; ist < ist1; ++ist) {
                  // skip surface tris flagged as FF...
                  if ((partVec[ipart]->surface_zn_flag == NULL)||(partVec[ipart]->surface_zn_flag[partVec[ipart]->surface->znost[ist]] != -1)) {
                    zone[nst_render] = partVec[ipart]->znosz[partVec[ipart]->surface->znost[ist]];
                    FOR_I3 spost[nst_render][i] = partVec[ipart]->surface->spost[ist][i];
                    ++nst_render;
                  }
                }
                scene->addSurfaceTris(partVec[ipart]->surface->xsp,spost,zone,nst_render);
                delete[] zone;
                delete[] spost;
              }
            }

            MPI_Barrier(mpi_comm);
            if (mpi_rank == 0) {
              double wtime = MPI_Wtime();
              cout << " tri time: " << wtime-wtime_prev << " [s], total time so far: " << wtime-wtime0 << " [s]" << endl;
              wtime_prev = wtime;
            }

            if (scene->hasGeomPlane()) {

              if (!scene->hasVar()) throw(1);
              if (!((scene->getVar() == "mesh")||(scene->getVar() == "level"))) throw(2);

              double xp[3],np[3];
              scene->getGeomPlaneXpAndNp(xp,np);
              const double mag = MAG(np); assert(mag > 0.0);
              const double unit_np[3] = {np[0]/mag,np[1]/mag,np[2]/mag};
              if (mpi_rank == 0) cout << "getGeomPlaneXpAndNp: " << COUT_VEC(xp) << " " << COUT_VEC(np) << endl;

              int count;
              if (b_vp) {
                assert(x_vd != NULL);
                count = ncv;
              }
              else {
                hcp.ensurePoints(xp,np);
                count = hcp.vertexVec.size();
              }
              double (*x_vv)[3] = new double[count][3];
              double *r_vv = new double[count];
              int * i_vv = NULL;
              double * v_vv = NULL;
              if (scene->getVar() == "mesh") {
                i_vv = new int[count];
              }
              else {
                assert(scene->getVar() == "level");
                v_vv = new double[count];
              }
              if (b_vp) {
                // need to prune points to those within thickened plane...
                if (scene->getVar() == "mesh") {
                  assert(i_vv != NULL);
                  count = 0;
                  for (int ip = 0; ip < ncv; ++ip) {
                    bool b_included = false;
                    if (wid.b_transform_cyl_x) {
                      const double xp0[3] = {x_vd[ip][0], x_vd[ip][1], x_vd[ip][2]};
                      const double r = sqrt(xp0[1]*xp0[1]+xp0[2]*xp0[2]);
                      if (fabs(r-xp[1]) <= delta_vd[ip]) b_included = true;
                    }
                    else if (wid.b_transform_cyl_z) {
                      const double xp0[3]= {x_vd[ip][0], x_vd[ip][1], x_vd[ip][2]};
                      const double r = sqrt(xp0[0]*xp0[0]+xp0[1]*xp0[1]);
                      if (fabs(r-xp[0]) <= delta_vd[ip]) b_included = true;
                    }
                    else {
                      const double dx[3] = DIFF(x_vd[ip],xp);
                      const double r = delta_vd[ip];
                      if (fabs(DOT_PRODUCT(dx,unit_np)) <= r) b_included = true;
                    }

                    if (b_included) {
                      FOR_I3 x_vv[count][i] = x_vd[ip][i];
                      r_vv[count] = delta_vd[ip];
                      int level = 0;
                      // its possible to have no hcp delta (e.g. BOX_WITH_CART_PTS)
                      if (hcp.b_hcp_delta) {
                        while (delta_vd[ip] < hcp.getPointsRvvForLevel(level)) ++level;
                      }
                      i_vv[count] = level;
                      ++count;
                    }
                  }
                }
                else {
                  assert(scene->getVar() == "level");
                  assert(v_vv != NULL);
                  count = 0;
                  for (int ip = 0; ip < ncv; ++ip) {
                    const double dx[3] = DIFF(x_vd[ip],xp);
                    const double r = delta_vd[ip];
                    if (fabs(DOT_PRODUCT(dx,unit_np)) <= r) {
                      FOR_I3 x_vv[count][i] = x_vd[ip][i];
                      r_vv[count] = r;
                      int level = 0;
                      if (hcp.b_hcp_delta) {
                        while (delta_vd[ip] < hcp.getPointsRvvForLevel(level)) ++level;
                      }
                      v_vv[count] = (double)level;
                      ++count;
                    }
                  }
                }
              }
              else if (hcp.b_plane_mesh_points) {
                // put all points into write...
                for (int iv = 0; iv < hcp.vertexVec.size(); ++iv) {
                  FOR_I3 x_vv[iv][i] = hcp.vertexVec[iv].x[i];
                  r_vv[iv] = hcp.getPointsRvvForLevel(hcp.vertexVec[iv].level);
                }
                if (scene->getVar() == "mesh") {
                  for (int iv = 0; iv < hcp.vertexVec.size(); ++iv) {
                    i_vv[iv] = (hcp.vertexVec[iv].window<<HCP_WINDOW_SHIFT_BITS)|hcp.vertexVec[iv].level;
                    // HACK for load balance visualization
                    //i_vv[iv] = ((mpi_rank+1)<<HCP_WINDOW_SHIFT_BITS)|LEVEL_MASK_BITS;
                  }
                }
                else {
                  assert(scene->getVar() == "level");
                  for (int iv = 0; iv < hcp.vertexVec.size(); ++iv)
                    v_vv[iv] = (double)hcp.vertexVec[iv].level;
                }
              }
              else if (hcp.b_all_mesh_points) {
                // need to prune points to those within thickened plane...
                if (scene->getVar() == "mesh") {
                  assert(i_vv != NULL);
                  count = 0;
                  for (int iv = 0; iv < hcp.vertexVec.size(); ++iv) {
                    bool b_included = false;
                    if (wid.b_transform_cyl_x) {
                      const double xp0[3]= {hcp.vertexVec[iv].x[0], hcp.vertexVec[iv].x[1], hcp.vertexVec[iv].x[2]};
                      const double r = sqrt(xp0[1]*xp0[1]+xp0[2]*xp0[2]);
                      const double delta = hcp.getPointsRvvForLevel(hcp.vertexVec[iv].level);
                      if (fabs(r-xp[1]) <= delta) b_included = true;
                    }
                    else if (wid.b_transform_cyl_z) {
                      const double xp0[3]= {hcp.vertexVec[iv].x[0], hcp.vertexVec[iv].x[1], hcp.vertexVec[iv].x[2]};
                      const double r = sqrt(xp0[0]*xp0[0]+xp0[1]*xp0[1]);
                      const double delta = hcp.getPointsRvvForLevel(hcp.vertexVec[iv].level);
                      if (fabs(r-xp[0]) <= delta) b_included = true;
                    }
                    else {
                      const double dx[3] = DIFF(hcp.vertexVec[iv].x,xp);
                      const double r = hcp.getPointsRvvForLevel(hcp.vertexVec[iv].level);
                      if (fabs(DOT_PRODUCT(dx,unit_np)) <= r) b_included = true;
                    }
                    if (b_included) {
                      FOR_I3 x_vv[count][i] = hcp.vertexVec[iv].x[i];
                      const double r = hcp.getPointsRvvForLevel(hcp.vertexVec[iv].level);
                      r_vv[count] = r;
                      i_vv[count] = (hcp.vertexVec[iv].window<<HCP_WINDOW_SHIFT_BITS)|hcp.vertexVec[iv].level;
                      ++count;
                    }
                  }
                }
                else {
                  assert(scene->getVar() == "level");
                  assert(v_vv != NULL);
                  count = 0;
                  for (int iv = 0; iv < hcp.vertexVec.size(); ++iv) {
                    const double dx[3] = DIFF(hcp.vertexVec[iv].x,xp);
                    const double r = hcp.getPointsRvvForLevel(hcp.vertexVec[iv].level);
                    if (fabs(DOT_PRODUCT(dx,unit_np)) <= r) {
                      FOR_I3 x_vv[count][i] = hcp.vertexVec[iv].x[i];
                      r_vv[count] = r;
                      v_vv[count] = (double)hcp.vertexVec[iv].level;
                      ++count;
                    }
                  }
                }
              }

              MPI_Barrier(mpi_comm);

              if (mpi_rank == 0) {
                double wtime = MPI_Wtime();
                if (b_vp)
                  cout << " vb point build time: " << wtime-wtime_prev << " [s], total time so far: " << wtime-wtime0 << " [s]" << endl;
                else
                  cout << " hcp point build time: " << wtime-wtime_prev << " [s], total time so far: " << wtime-wtime0 << " [s]" << endl;
                wtime_prev = wtime;
              }

              // use sign of delta as a flag...
              int full_count = count;
              int indexed_count = count;
              if (!b_vp) {
                prepareIsInside();
                for (int ipart = 0; ipart < partVec.size(); ++ipart) {
                  if (partVec[ipart]->pts) {
                    full_count += partVec[ipart]->pts->np;
                    for (int ip = 0; ip < partVec[ipart]->pts->np; ++ip) {
                      const double dx[3] = DIFF(partVec[ipart]->pts->xp[ip],xp);
                      const double r = partVec[ipart]->pts->delta[ip];
                      assert(r > 0.0);
                      if (fabs(DOT_PRODUCT(dx,unit_np)) <= r) {
                        // we are going to include this point in the build, if and only if it
                        // is NOT inside the FF of every part later than us in the list (these parts take priority)
                        // AND it is not inside the solid of any part. Since we have checked the FF of the
                        // later parts, we can just check the solid of earlier parts...
                        for (int ipart2 = ipart+1; ipart2 < partVec.size(); ++ipart2) {
                          if (partVec[ipart2]->hasFF()) {
                            if (partVec[ipart2]->isInsideFF(partVec[ipart]->pts->xp[ip])) {
                              // HACKHACK
                              partVec[ipart]->pts->delta[ip] = -r;
                              break;
                            }
                          }
                          else {
                            if (partVec[ipart2]->isInsideSolid(partVec[ipart]->pts->xp[ip],false)) { // ipart2!=0
                              partVec[ipart]->pts->delta[ip] = -r;
                              break;
                            }
                          }
                        }
                        if (partVec[ipart]->pts->delta[ip] > 0.0) {
                          for (int ipart2 = 0; ipart2 < ipart; ++ipart2) {
                            if (partVec[ipart2]->isInsideSolid(partVec[ipart]->pts->xp[ip],ipart2==0)) {
                              partVec[ipart]->pts->delta[ip] = -r;
                              break;
                            }
                          }
                        }

                        if (partVec[ipart]->pts->delta[ip] > 0.0)
                          ++indexed_count;

                      }
                      else {
                        // this point is on the other side of the plane...
                        partVec[ipart]->pts->delta[ip] = -r;
                      }
                    }
                  }
                }
              }

              int* ivv_global = NULL;
              if (scene->getVar() == "mesh") {
                assert(i_vv != NULL);

                // build ivv_global...
                ivv_global = new int[full_count];
                int8 * vvora = NULL;
                buildXora(vvora,indexed_count);
                if (mpi_rank == 0) cout << " points within thickened plane: " << vvora[mpi_size] << endl;
                // actually needs to be <= 7 digits to not saturate float buffer
                assert(vvora[mpi_size] < TWO_BILLION);
                int index = 0;
                for (int iv = 0; iv < count; ++iv) {
                  ivv_global[iv] = (int)vvora[mpi_rank]+index;
                  ++index;
                }
                if (!b_vp) {
                  int disp = 0;
                  for (int ipart = 0; ipart < partVec.size(); ++ipart) {
                    if (partVec[ipart]->pts) {
                      for (int ip = 0; ip < partVec[ipart]->pts->np; ++ip) {
                        ivv_global[count+disp+ip] = (int)vvora[mpi_rank]+index;
                        if (partVec[ipart]->pts->delta[ip] > 0.0)
                          ++index;
                      }
                      disp += partVec[ipart]->pts->np;
                    }
                  }
                }
                delete[] vvora;

                scene->addPlaneMeshRvvPositive(x_vv,r_vv,i_vv,ivv_global,count,2.0);
                delete[] i_vv;
              }
              else {
                assert(scene->getVar() == "level");
                assert(v_vv != NULL);
                scene->addPlaneDataRvvPositive(v_vv,x_vv,r_vv,count,2.0);
                delete[] v_vv;
              }
              delete[] x_vv;
              delete[] r_vv;
              MPI_Barrier(mpi_comm);
              if (mpi_rank == 0) {
                double wtime = MPI_Wtime();
                if (b_vp)
                  cout << " vb scene build time: " << wtime-wtime_prev << " [s], total time so far: " << wtime-wtime0 << " [s]" << endl;
                else
                  cout << " hcp scene build time: " << wtime-wtime_prev << " [s], total time so far: " << wtime-wtime0 << " [s]" << endl;
                wtime_prev = wtime;
              }

              // also add any part-based pts that we have...
              if (!b_vp) {
                int disp = 0;
                for (int ipart = 0; ipart < partVec.size(); ++ipart) {
                  if (partVec[ipart]->pts) {
                    if (scene->getVar() == "mesh") {
                      assert(ivv_global != NULL);
                      assert(i_vv != NULL);
                      int * flag = new int[partVec[ipart]->pts->np];
                      for (int ip = 0; ip < partVec[ipart]->pts->np; ++ip) {
                        flag[ip] = ((hcp.hcpWindowVec.size()+ipart+1)<<HCP_WINDOW_SHIFT_BITS)|LEVEL_MASK_BITS;
                        // HACK for load balance visualization
                        //flag[ip] = ((mpi_rank+1)<<HCP_WINDOW_SHIFT_BITS)|LEVEL_MASK_BITS;
                      }
                      scene->addPlaneMeshRvvPositive(partVec[ipart]->pts->xp,partVec[ipart]->pts->delta,flag,ivv_global+count+disp,partVec[ipart]->pts->np,1.002);
                      delete[] flag;
                      disp += partVec[ipart]->pts->np;
                    }
                    else {
                      assert(scene->getVar() == "level");
                      assert(v_vv != NULL);
                      double * var = new double[partVec[ipart]->pts->np];
                      for (int ip = 0; ip < partVec[ipart]->pts->np; ++ip) {
                        int level = 0;
                        if (partVec[ipart]->pts->delta[ip] > 0.0) {
                          while (partVec[ipart]->pts->delta[ip] < hcp.getPointsRvvForLevel(level))
                            ++level;
                        }
                        var[ip] = (double)level;
                      }
                      scene->addPlaneDataRvvPositive(var,partVec[ipart]->pts->xp,partVec[ipart]->pts->delta,partVec[ipart]->pts->np,2.0);
                      delete[] var;
                    }
                  }
                }
                for (int ipart = 0; ipart < partVec.size(); ++ipart) {
                  if (partVec[ipart]->pts) {
                    for (int ip = 0; ip < partVec[ipart]->pts->np; ++ip)
                      partVec[ipart]->pts->delta[ip] = fabs(partVec[ipart]->pts->delta[ip]);
                  }
                }
                // timing...
                MPI_Barrier(mpi_comm);
                if (mpi_rank == 0) {
                  double wtime = MPI_Wtime();
                  cout << " part scene build time: " << wtime-wtime_prev << " [s], total time so far: " << wtime-wtime0 << " [s]" << endl;
                  wtime_prev = wtime;
                }
              }
              DELETE(ivv_global);
            }
            else {

              // look at ff surfaces? note this can only be used when we are NOT doing the mesh viz above (aux taken)...
              // TODO expand zoneVec to include these and write out as surface data
              if (checkParam("DISPLAY_FF_SURFACES")) {
                for (int ipart = 0; ipart < partVec.size(); ++ipart) {
                  if (partVec[ipart]->ff_surface&&partVec[ipart]->pts) {

                    Points *pts = partVec[ipart]->pts;
                    SurfaceShm* ff_surface = partVec[ipart]->ff_surface;

                    // TODO: can we reduce this at all to a local subsurface?
                    int * zone = new int[ff_surface->nst]; // must be less than number of local points
                    for (int ist = 0; ist < ff_surface->nst; ++ist) {
                      //zone[ist] = 65534-32768; // which zone should I use?
                      zone[ist] = zoneVec.size();
                    }

                    scene->addSurfaceTris(ff_surface->xsp,ff_surface->spost,zone,ff_surface->nst);

                    double (*xp)[3] = new double[pts->np][3];
                    double *delta = new double[pts->np];
                    int *index = new int[pts->np];
                    int np_render = 0;
                    for (int ip = 0; ip < pts->np; ++ip) {
                      //if (pts->flag[ip] == 1) {
                      //  assert(pts->flag[ip] < ff_surface->nst);
                      FOR_I3 xp[np_render][i] = pts->xp[ip][i];
                      delta[np_render] = pts->delta[ip];
                      index[np_render] = zone[0];
                      ++np_render;
                      //}
                      //else {
                      //  assert(pts->flag[ip] == 0);
                      //}
                    }
                    delete[] zone;

                    int* ip_global = new int[np_render];
                    int8 * ipora = NULL;
                    buildXora(ipora,np_render);
                    for (int ip = 0; ip < np_render; ++ip)
                      ip_global[ip] = (int)ipora[mpi_rank]+ip;
                    delete[] ipora;

                    scene->addSurfaceMeshRvvPositive(xp,delta,index,ip_global,np_render,2.001);
                    delete[] xp;
                    delete[] delta;
                    delete[] ip_global;

                  }
                  }
                }
              }

            }
            else if (scene->hasGeomPlane()&&scene->hasVar()) {

              if (!b_vp)
                throw(3);

              vector<SimpleTri> meshTriVec;
              double geom_plane_xp[3],geom_plane_np[3];
              scene->getGeomPlaneXpAndNp(geom_plane_xp,geom_plane_np);

              vector<pair<SimpleTriWithData,int> > boundaryTriVec; // int is zone

              // 3d Voronoi mesh imaging
              assert(scene->getVar() == "3d_mesh"); // from above
              scene->removeGeomPlaneAsBlankPlane();  // only blank based on rendered tris, not the geom plane itself

              assert(vdReturnDataVec.size() == ncv);
              for (int ii = 0; ii < ncv; ++ii) {
                const int icv = vdReturnDataVec[ii].icv;
                assert((icv >= 0)&&(icv < ncv));
                // if this coordinate is behind the plane, render its boundary cvs, and if it is
                // within its delta of the plane, render its internal tris...
                const double dn =
                  (x_vd[icv][0]-geom_plane_xp[0])*geom_plane_np[0] +
                  (x_vd[icv][1]-geom_plane_xp[1])*geom_plane_np[1] +
                  (x_vd[icv][2]-geom_plane_xp[2])*geom_plane_np[2];
                if ((dn <= 0)||(!scene->blankDataPlane())) {
                  // everyone behind the requested plane gets their boundary faces rendered...
                  addSimpleTrisBoundary(boundaryTriVec,ii,x_vd[icv]);
                }
                if ((dn <= 0)&&(vdReturnDataVec[ii].hasPeriodicNbrs()||(dn >= -2.0*delta_vd[icv]))) {
                  // anyone close enought to the plane to potentially have nbrs over the
                  // plane gets their internal faces rendered as well...
                  addSimpleTrisInternal(meshTriVec,ii,x_vd[icv]);
                }
              }
              // draw on first edge of tris only
              scene->addSimpleInternalTrisWithMesh(meshTriVec,1,true);
              scene->addSimpleSurfaceTrisWithMesh(boundaryTriVec,1,true); // prevent data plane surface blanking

            }

#ifdef WITH_CHRONO
            auto t1 = std::chrono::high_resolution_clock::now();
#endif

            // --------------------------------------
            // finally, write the image...
            // --------------------------------------

            if (scene->isClientRequest()) {
              // don't add step index to filename when interacting with Cascade client...
              scene->writeImage();
            }
            else {
              scene->writeImage(step);
            }
            MPI_Barrier(mpi_comm);
            if (mpi_rank == 0) {
              double wtime = MPI_Wtime();
              cout << " write image time: " << wtime-wtime_prev << " [s], total time so far: " << wtime-wtime0 << " [s]" << endl;
              wtime_prev = wtime;
            }

            delete scene;

#ifdef WITH_CHRONO
            auto t2 = std::chrono::high_resolution_clock::now();
            {
              std::chrono::duration<double> elapsed1 = t1-t0;
              std::chrono::duration<double> elapsed2 = t2-t1;
              double my_buf[2];
              my_buf[0] = elapsed1.count();
              my_buf[1] = elapsed2.count();
              double buf[2];
              MPI_Reduce(my_buf,buf,2,MPI_DOUBLE,MPI_SUM,0,mpi_comm);
              if (mpi_rank == 0) {
                cout << "[TIMING]: " << buf[0]/double(mpi_size) << " " << buf[1]/double(mpi_size) << endl;
              }
            }
#endif

          }
          catch(int e) {
            if (e == 3) {
              WUI(WARN,"expecting NSMOOTH <int> prior to WRITE_IMAGE VAR 3d_mesh");
            }
            else {
              WUI(WARN,"WRITE_IMAGE GEOM_PLANE <x> <y> <z> <nx> <ny> <nz> expects VAR <mesh/level>");
            }
            helpWriteImage();
          }
        }
      }

      void helpWriteImage() {
        WUI(INFO,
            "WRITE_IMAGE is highly-customizable png-writer that can sample your data in many ways. Some examples:\n" <<
            "WRITE_IMAGE NAME plane_x_mesh GEOM PLANE 0.5 0.5 0.5 1 0 0 VAR mesh\n" <<
            "WRITE_IMAGE NAME plane_y_level GEOM PLANE 0.5 0.5 0.5 0 1 0 VAR level\n" <<
            "for more detail see [$CWIKB:write_image]");
      }

  void processForcePlanar(Param * param,const bool b_help) {
    if (b_help) {
      WUI(INFO,"FORCE_PLANAR zone1,zone2...");
    }
    else {
      COUT1("processing FORCE_PLANAR...");
      assert(force_planar_zone_flag == NULL);
      const int nzone = PartData::zoneVec.size();
      force_planar_zone_flag = new int[nzone];
      for (int izn = 0; izn < nzone; ++izn) {
	force_planar_zone_flag[izn] = 0;
      }
      vector<string> nameVec;
      MiscUtils::splitCsv(nameVec,param->getString());
      if (nameVec.empty()) {
	CERR("FORCE_PLANAR requires atleast one valid zone name");
      }
      for (int ii = 0, nii = nameVec.size(); ii < nii; ++ii) {
	int izn;
	for (izn = 0; izn < nzone; ++izn) {
	  if (nameVec[ii] == PartData::zoneVec[izn].getName()) {
	    COUT1(" > zone " << nameVec[ii] << " will be written as planar bfs");
	    force_planar_zone_flag[izn] = 1;
	    break;
	  }
	}
	if (izn == nzone) {
	  CERR("FORCE_PLANAR could not match zone " << nameVec[ii]);
	}
      }
    }
  }
  
  void processWriteMles(Param * param,const bool b_help) {
        if (b_help) {
          helpWriteMles();
        }
        else {
          try {

            // need to build voronoi first...
            if (!b_vp)
              throw(1);

            int version;
            string filename;
            bool write_surface;
            const string token = MiscUtils::toUpperCase(param->getName());
            if (MiscUtils::toUpperCase(param->getName()) == "WRITE_RESTART_V3") {
              filename = "restart_v3.les";
              version  = 3;
              write_surface = false;
            }
            else if (MiscUtils::toUpperCase(param->getName()) == "WRITE_RESTART_V4") {
              filename = "restart_v4.les";
              version  = 4;
              write_surface = false;
            }
            else {
              filename = "restart.mles";
              version  = 5;
              write_surface = true;
              bool got_filename = false;
              int iarg = 0;
              while (iarg < param->size()) {
                string token = param->getString(iarg++);
                if (MiscUtils::toUpperCase(token) == "NO_SURFACE") {
                  write_surface = false;
                }
                else {
                  // assume this is a filename...
                  if (got_filename) {
                    CWARN("Filename already set. Unrecognized parameter: " << token);
                  }
                  else {
                    filename = token;
                    got_filename = true;
                  }
                }
              }
            }

            WUI(INFO,"WRITE_MLES filename=" << filename << ", version=" << version << ", surface=" << write_surface);

            writeMles(filename,version,write_surface);
          }
          catch(int e) {
            if (e == 1) {
              WUI(WARN,"expecting NSMOOTH <int> prior to WRITE_MLES <string>");
            }
            else {
              WUI(WARN,"expecting WRITE_MLES <string>");
            }
            helpWriteMles();
          }
        }

      }

      void helpWriteMles() {
        WUI(INFO,
            "WRITE_MLES writes the current voronoi points to disk. Examples:\n" <<
            "WRITE_MLES output.mles # also accepts WRITE_RESTART and WRITE_RESTART_V5\n" <<
            "WRITE_MLES # filename assumed to be restart.mles\n" <<
            "WRITE_RESTART_V3 # legacy version 3 format, filename assumed to be restart_v3.mles\n" <<
            "WRITE_RESTART_V4 # legacy version 4 format, filename assumed to be restart_v4.mles\n" <<
            "for more detail see [$CWIKB:stitch_export]"); // TODO change link and knowledge base to WRITE_MLES
      }

      void processClearSmoothing(Param *param, const bool b_help) {
        if (b_help) {
          WUI(INFO,
              "CLEAR_SMOOTHING returns points to original seed locations.\n"<<
              "for more detail see [$CWIKB:stitch_smoothing]");
        }
        else {
          b_vp = false; // cleared below
        }
      }

      void processSmoothLimit(Param * param,const bool b_help) {
        if (b_help) {
          helpSmoothLimit();
        }
        else {
          try {
            const double value = param->getDouble();
            if (value < 0.0) throw(1);
            smooth_limit = value;
            b_vp = false; // reset smoothing
            WUI(INFO,"SMOOTH_LIMIT set to " << smooth_limit);
          }
          catch(int e) {
            WUI(WARN,"expecting SMOOTH_LIMIT <+double>");
            helpSmoothLimit();
          }
        }
      }

      void helpSmoothLimit() {
        WUI(INFO,
            "SMOOTH_LIMIT sets the max normalized displacement of the seed points from their original location. Example:\n" <<
            "SMOOTH_LIMIT 0.5 # (default) allows points to move 0.5*delta\n" <<
            "for more detail see [$CWIKB:stitch_smoothing]");
      }

      void processSmoothMode(Param * param,const bool b_help) {
        if (b_help) {
          helpSmoothMode();
        }
        else {
          try {
            int iarg = 0;
            const string token = MiscUtils::toUpperCase(param->getString(iarg++));
            if (token == "ALL") {
              // smooth_mode = ALL_SMOOTH_MODE;
              FOR_I4 smooth_mode_idata[i] = 0;
              WUI(INFO,"SMOOTH-ing will be applied to ALL cells");
            }
            else if (token == "BLAYER") {
              const int n = param->getInt(iarg++);
              if (n < 0) {
                WUI(WARN,"SMOOTH_MODE BLAYER requires <n> >= 0; setting to 0");
              }
              else {
                // smooth_mode = BLAYER_SMOOTH_MODE;
                smooth_mode_idata[0] = n;
                WUI(INFO,"SMOOTH-ing will be applied to " << smooth_mode_idata[0] << " BLAYER cells");
              }
            }
            else if (token == "TRANSITION") {
              const int n = param->getInt(iarg++);
              if (n < 0) {
                WUI(WARN,"SMOOTH_MODE TRANSITION requires <n> >= 0; setting to 0");
              }
              else {
                // smooth_mode = BLAYER_SMOOTH_MODE;
                smooth_mode_idata[1] = n;
                WUI(INFO,"SMOOTH-ing will be applied to " << smooth_mode_idata[1] << " TRANSITION cells");
              }
            }
            else if (token == "PERIODIC") {
              const int n = param->getInt(iarg++);
              if (n < 0) {
                WUI(WARN,"SMOOTH_MODE PERIODIC requires <n> >= 0; setting to 0");
              }
              else {
                // smooth_mode = BLAYER_SMOOTH_MODE;
                smooth_mode_idata[2] = n;
                WUI(INFO,"SMOOTH-ing will be applied to " << smooth_mode_idata[2] << " PERIODIC cells");
              }
            }
            else if (token == "PART_LAYER") {
              const int n = param->getInt(iarg++);
              if (n < 0) {
                WUI(WARN,"SMOOTH_MODE PART_LAYER requires <n> >= 0; setting to 0");
              }
              else {
                smooth_mode_idata[3] = n;
                WUI(INFO,"SMOOTH-ing will be applied to " << smooth_mode_idata[3] << " PART_LAYER cells");
              }
            }
            else if (token == "BLAYER_TRANSITION") {
              // keep backwards compatible
              const int nblayer = param->getInt(iarg++);
              const int ntransition = param->getInt(iarg++);
              if (nblayer < 0) {
                WUI(WARN,"SMOOTH_MODE BLAYER_TRANSITION requires <nblayer> >= 0");
              }
              else if (ntransition < 0) {
                WUI(WARN,"SMOOTH_MODE BLAYER_TRANSITION requires <ntransition> >= 0");
              }
              else {
                // smooth_mode = TRANSITION_SMOOTH_MODE;
                smooth_mode_idata[0] = smooth_mode_idata[2] = smooth_mode_idata[3] = nblayer;
                smooth_mode_idata[1] = ntransition;
                WUI(INFO,"SMOOTH-ing will be applied to " << smooth_mode_idata[0] << " BLAYER and " << smooth_mode_idata[1] << " TRANSITIONS cells");
              }
            }
            else if (token == "HCP_ONLY") {
              FOR_I3 smooth_mode_idata[i] = -1;
              WUI(INFO,"SMOOTH_MODE set to HCP_ONLY (recommended)");
            }
            else if (token == "FORCE_ALL") {
              b_force_smooth_all = true;
              FOR_I4 smooth_mode_idata[i] = 0;
              WUI(INFO,"SMOOTH_MODE set to FORCE_ALL (not recommended)");
            }
            else {
              WUI(WARN,"unrecognized SMOOTH_MODE token: " << token);
              helpSmoothMode();
            }
          }
          catch(int e) {
            WUI(WARN,"Error parsing SMOOTH_MODE tokens");
            helpSmoothMode();
          }
        }
      }

      void helpSmoothMode() {
        WUI(INFO,
            "SMOOTH_MODE sets the way Voronoi seed points are selected for smoothing. Note that smoothing will be\n" <<
            "constrained for seed points coming from parts, unless you choose FORCE_ALL mode (not recommended). Examples:\n" <<
            "SMOOTH_MODE BLAYER 3 # (default) unconstrained seed points within 3 cells of any boundary will be smoothed\n" <<
            "SMOOTH_MODE TANSITION 3 1 # unconstrained seed points within 3 cells of any boundary and 1 cell of the transtion will be smoothed\n" <<
            "SMOOTH_MODE ALL # all unconstrained seed points will be smoothed\n" <<
            "SMOOTH_MODE FORCE_ALL # smooth all points, including constrained points. Good for illustrating the stability seeding patterns to smoothing\n"
            "for more detail see [$CWIKB:stitch_smoothing]");
      }

      void processCreaseAngleDegrees(Param * param,const bool b_help) {
        if (b_help) {
          helpCreaseAngleDegrees();
        }
        else {
          try {
            const double value = param->getDouble();
            if (value < 0.0) throw(1);
            crease_angle_degrees = value;
            WUI(INFO,"CREASE_ANGLE_DEGREES set to " << crease_angle_degrees);
          }
          catch(int e) {
            WUI(WARN,"expecting CREASE_ANGLE_DEGREES <+double>");
            helpCreaseAngleDegrees();
          }
        }
      }

      void helpCreaseAngleDegrees() {
        WUI(INFO,
            "CREASE_ANGLE_DEGREES sets the triangle crease angle at which we create new boundary face. Example:\n" <<
            "CREASE_ANGLE_DEGREES 175 # (default)\n" <<
            "for more detail see [$CWIKB:stitch_smoothing]");
      }

      void processRcLocalThreshold(Param * param,const bool b_help) {
        if (b_help) {
          helpRcLocalThreshold();
        }
        else {
          try {
            const int8 value = param->getInt8();
            if (value < 0) 
              throw(1);
            rc_local_threshold = value;
            WUI(INFO,"RC_LOCAL_THRESHOLD set to " << rc_local_threshold);
          }
          catch(int e) {
            WUI(WARN,"expecting RC_LOCAL_THRESHOLD <+int>");
            helpRcLocalThreshold();
          }
        }
      }
      void helpRcLocalThreshold() {
        WUI(INFO,
            "RC_LOCAL_THRESHOLD sets the cost threshold for locally built Voronoid diagrams. Example:\n" <<
            "RC_LOCAL_THRESHOLD 0 # (default) load balance all Voronoi diagrams\n" <<
            "for more detail see [$CWIKB:stitch_smoothing]");
      }

      void processRcLimit(Param * param,const bool b_help) {
        if (b_help) {
          helpRcLimit();
        }
        else {
          try {
            const int8 value = param->getInt8();
            if (value <= 0) 
              throw(1);
            rc_limit = value;
            WUI(INFO,"RC_LIMIT set to " << rc_limit);
          }
          catch(int e) {
            WUI(WARN,"expecting RC_LIMIT <+int>");
            helpRcLimit();
          }
        }
      }
      void helpRcLimit() {
        WUI(INFO,
            "RC_LIMIT sets the cost limit for Voronoi points to be cut during this iteration. Example:\n" <<
            "RC_LIMIT 100000000 # (default)\n" <<
            "for more detail see [$CWIKB:stitch_smoothing]");
      }

    } // namespace StitchNS
