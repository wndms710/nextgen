#include "CTI.hpp"
using namespace CTI;

#include "FluentReader.hpp"
#include "SimpleFunc.hpp"
#include "Adt.hpp"
#include "DistributedDataExchanger.hpp"
#include "Prcomm.hpp"
#include "Utils.hpp"
#include "GaussQuadrature.hpp" // tet4 and tri3

class Fea {
  private:
    
    int *send_count,*send_disp;
    int *recv_count,*recv_disp;
    Adt<double> * adt;
    Adt<double> * bbox_adt;
    uint8 * rbi_g; // nodal ghost rbi_g
    vector<CvPrcomm> noPrcommVec; // TODO: call CvPrcomm->Prcomm, and Prcomm->PrcommFull
    map<const void*,MpiRequestStuff*> mpiRequestMap;

    // for striping io
    int8 nno_global;
    int8* ino_global;
    int8* noora_striped;
    DistributedDataExchanger* dde_striped;

    int * noono_i;
    int * noono_v;

    int nzn; // number of volume zones...
    int nno0; // active corner nodes
    int nno; // active nodes
    int nno_g; // active+ghost nodes (nno_g >= nno)
    int nno_i; // active internal nodes: edge nbrs are also active (nno_i <= nno). Used for latency hiding
    int nte; // all tets
    int nte_i; // tets touching local nodes only (nte_i <= nte_ib)
    int nte_ib; // tets touching local and ghost nodes that we own (nte_ib <= nte)
    int (*noote)[10];
    int *znote;
    double (*x_no)[3];
    double *T_no;
    double (*s_no_d)[3]; // sxx,syy,szz
    double (*s_no_od)[3]; // sxy,sxz,syz

    // for setting pin constraints
    double dxminmax[2]; // min and max edge lengths
    double bbminmax[6]; // xmin,ymin,zmin,xmax,ymax,zmax

    int tet10otri6[4][6]; // mapping from 6 node tri + face to 10 node tet
    double T0; // reference temperature

    vector<SimpleFunc*> EFuncVec; // Young's modulus
    vector<SimpleFunc*> nuFuncVec; // Poisson's ratio
    vector<SimpleFunc*> alphaFuncVec; // thermal expansion coefficient

    // Kd = f...
    double (*K_no)[3][3]; // global stiffness matrix (noono_i[nno]])
    double (*d_no)[3]; // global displacements (nno_g)
    double (*f_no)[3]; // global external forces (nno_g)

  public:

    Fea() {

      send_count = NULL;
      send_disp = NULL;
      recv_count = NULL;
      recv_disp = NULL;
      adt = NULL;
      bbox_adt = NULL;
      rbi_g = NULL;

      nno_global = 0;
      ino_global = NULL;
      noora_striped = NULL;
      dde_striped = NULL;

      noono_i = NULL;
      noono_v = NULL;

      nno0 = nno = nno_g = nno_i = 0;
      nte = nte_i = nte_ib = 0;
      nzn = 0;

      noote = NULL;
      znote = NULL;
      x_no = NULL;
      T_no = NULL;
      s_no_d = NULL;
      s_no_od = NULL;

      K_no = NULL;
      d_no = NULL;
      f_no = NULL;

      FOR_I2 dxminmax[i] = HUGE_VAL;
      FOR_I6 bbminmax[i] = HUGE_VAL;

      T0 = getDoubleParam("T0",0.0);

      tet10otri6[0][0] = 0; tet10otri6[0][1] = 1; tet10otri6[0][2] = 2; 
      tet10otri6[0][3] = 4; tet10otri6[0][4] = 5; tet10otri6[0][5] = 6;

      tet10otri6[1][0] = 0; tet10otri6[1][1] = 3; tet10otri6[1][2] = 1; 
      tet10otri6[1][3] = 7; tet10otri6[1][4] = 8; tet10otri6[1][5] = 4;

      tet10otri6[2][0] = 3; tet10otri6[2][1] = 0; tet10otri6[2][2] = 2; 
      tet10otri6[2][3] = 7; tet10otri6[2][4] = 6; tet10otri6[2][5] = 9;

      tet10otri6[3][0] = 1; tet10otri6[3][1] = 3; tet10otri6[3][2] = 2; 
      tet10otri6[3][3] = 8; tet10otri6[3][4] = 9; tet10otri6[3][5] = 5;
    }

    virtual ~Fea() {
      DELETE(send_count);
      DELETE(send_disp);
      DELETE(recv_count);
      DELETE(recv_disp);
      if (adt != NULL) {
        delete adt;
        adt = NULL;
      }
      if (bbox_adt != NULL) {
        delete bbox_adt;
        bbox_adt = NULL;
      }
      DELETE(rbi_g);

      DELETE(ino_global);
      DELETE(noora_striped);
      if (dde_striped) {
        delete dde_striped;
        dde_striped = NULL;
      }

      DELETE(noono_i);
      DELETE(noono_v);

      DELETE(noote);
      DELETE(znote);
      DELETE(x_no);
      DELETE(T_no);
      DELETE(s_no_d);
      DELETE(s_no_od);

      DELETE(K_no);
      DELETE(d_no);
      DELETE(f_no);

      for (int ii = 0, size = EFuncVec.size(); ii < size; ++ii) 
        delete EFuncVec[ii];
      for (int ii = 0, size = nuFuncVec.size(); ii < size; ++ii) 
        delete nuFuncVec[ii];
      for (int ii = 0, size = alphaFuncVec.size(); ii < size; ++ii) 
        delete alphaFuncVec[ii];
    }

  // ===========================================================
  // ghost node data update routines...
  // ===========================================================

#define UPDATE_START updateNoDataStart
#define UPDATE_FINISH updateNoDataFinish
#define UPDATE updateNoData
#define PRCOMM_VEC noPrcommVec

  // updateNoData(int * s...

#define T int
#define PACK_BUF pack_buf_int
#define MPI_T MPI_INT
#define UPDATE_TAG 12121
#include "updateDN.hpp"
#undef T
#undef PACK_BUF
#undef MPI_T
#undef UPDATE_TAG

  // updateCvData(double * s and double (*s)[3]...

#define T double
#define PACK_BUF pack_buf_double
#define MPI_T MPI_DOUBLE
#define UPDATE_TAG 12122
#include "updateDN.hpp"
#include "updateDN3.hpp"
#undef T
#undef PACK_BUF
#undef MPI_T
#undef UPDATE_TAG

#undef UPDATE_START
#undef UPDATE_FINISH
#undef UPDATE
#undef PRCOMM_VEC

  // ===========================================================
  // end of ghost node data update routines
  // ===========================================================
  
  // -----------------------------------
  // 6 node tri:
  // 
  //    2
  //  5/ \4
  //  /_ _\
  // 0  3  1
  // 
  // 10 node tet (outward normal):
  // 
  //    2        1        2        2  
  //  6/ \5    4/ \8    9/ \6    5/ \9 
  //  /_ _\    /_ _\    /_ _\    /_ _\
  // 0  4  1  0  7  3  3  7  0  1  8  3
  // -----------------------------------
  
    void getNTri(double NTri[3][3*6],const double ksi[4]) {
      FOR_I3 for (int j = 0; j < 3*10; ++j) NTri[i][j] = 0.0;

      FOR_I3 FOR_J3 NTri[j][i*3+j] = 2.0*ksi[i]*(ksi[i]-0.5); 

      FOR_J3 NTri[j][3*3+j] = 4.0*ksi[0]*ksi[1];
      FOR_J3 NTri[j][4*3+j] = 4.0*ksi[1]*ksi[2];
      FOR_J3 NTri[j][5*3+j] = 4.0*ksi[0]*ksi[2];
    }

    void getdNTri(double dNTri[2][6],const double ksi[4]) {
      FOR_I2 FOR_J6 dNTri[i][j] = 0.0;
      
      FOR_I2 dNTri[i][i] = 4.0*ksi[i]-1.0;
      FOR_I2 dNTri[i][2] = 1.0-4.0*ksi[2];

      dNTri[0][3] = 4.0*ksi[1];
      dNTri[1][3] = 4.0*ksi[0];

      dNTri[0][4] = -4.0*ksi[1];
      dNTri[1][4] = 4.0*(ksi[2]-ksi[1]);

      dNTri[0][5] = 4.0*(ksi[2]-ksi[0]);
      dNTri[1][5] = -4.0*ksi[0];
    }

    double getdSTri(const double dNTri[2][6],const int ite,const int itri_local) {
      double dxdksi[2][3] = {{0.0,0.0,0.0},{0.0,0.0,0.0}};
      FOR_I6 {
        const int ino = noote[ite][tet10otri6[itri_local][i]];
        FOR_J2 FOR_K3 dxdksi[j][k] += dNTri[j][i]*x_no[ino][k];
      }
      return CROSS_PRODUCT_MAG(dxdksi[0],dxdksi[1]);
    }

    void getN(double N[3][3*10],const double ksi[5]) {
      FOR_I3 for (int j = 0; j < 3*10; ++j) N[i][j] = 0.0;

      FOR_I4 FOR_J3 N[j][i*3+j] = 2.0*ksi[i]*(ksi[i]-0.5); 

      FOR_J3 N[j][4*3+j] = 4.0*ksi[0]*ksi[1];
      FOR_J3 N[j][5*3+j] = 4.0*ksi[1]*ksi[2];
      FOR_J3 N[j][6*3+j] = 4.0*ksi[0]*ksi[2];
      FOR_J3 N[j][7*3+j] = 4.0*ksi[0]*ksi[3];
      FOR_J3 N[j][8*3+j] = 4.0*ksi[1]*ksi[3];
      FOR_J3 N[j][9*3+j] = 4.0*ksi[2]*ksi[3];
    }

    void getdN(double dN[3][10],const double ksi[5]) {
      FOR_I3 FOR_J10 dN[i][j] = 0.0;
      
      FOR_I3 dN[i][i] = 4.0*ksi[i]-1.0;
      FOR_I3 dN[i][3] = 1.0-4.0*ksi[3];

      dN[0][4] = 4.0*ksi[1];
      dN[1][4] = 4.0*ksi[0];

      dN[1][5] = 4.0*ksi[2];
      dN[2][5] = 4.0*ksi[1];

      dN[0][6] = 4.0*ksi[2];
      dN[2][6] = 4.0*ksi[0];

      dN[0][7] = 4.0*(ksi[3]-ksi[0]);
      dN[1][7] = -4.0*ksi[0];
      dN[2][7] = -4.0*ksi[0];

      dN[0][8] = -4.0*ksi[1];
      dN[1][8] = 4.0*(ksi[3]-ksi[1]);
      dN[2][8] = -4.0*ksi[1];

      dN[0][9] = -4.0*ksi[2];
      dN[1][9] = -4.0*ksi[2];
      dN[2][9] = 4.0*(ksi[3]-ksi[2]);
    }

    void getJ(double J[3][3],const double dN[3][10],const int ite) {
      FOR_I3 FOR_J3 {
        J[i][j] = 0.0;
        for (int ino_local = 0; ino_local < 10; ++ino_local)
          J[i][j] += dN[i][ino_local]*x_no[noote[ite][ino_local]][j];
      }
    }

    void getB(double B[6][3*10],const double J[3][3],const double dN[3][10],const int ite) {
      double Jinv[3][3]; MiscUtils::invertMat(Jinv,J);
      double der[3][10]; 
      FOR_I3 FOR_J10 for (int ino_local = 0; ino_local < 10; ++ino_local) {
        der[i][ino_local] = 0.0;
        FOR_J3 der[i][ino_local] += Jinv[i][j]*dN[j][ino_local];
      }
      FOR_I6 for (int j = 0; j < 3*10; ++j) B[i][j] = 0.0;
      for (int ino_local = 0; ino_local < 10; ++ino_local) {
        const int k = 3*ino_local;
        FOR_J3 B[j][k+j] = der[j][ino_local];
        
        B[3][k] = der[1][ino_local];
        B[3][k+1] = der[0][ino_local];

        B[4][k+1] = der[2][ino_local];
        B[4][k+2] = der[1][ino_local];

        B[5][k] = der[2][ino_local];
        B[5][k+2] = der[0][ino_local];
      }

    }

    void getKeAndFe(double Ke[3*10][3*10],double fe[3*10],const int ite) {

      // constiutive matrix...
      double D[6][6]; 
      FOR_I6 FOR_J6 D[i][j] = 0.0;

      for (int i = 0; i < 3*10; ++i)
        for (int j = 0; j < 3*10; ++j)
          Ke[i][j] = 0.0;
      for (int i = 0; i < 3*10; ++i)
        fe[i] = 0.0;
      double J[3][3]; // Jacobian matrix
      double B[6][3*10]; // strain-displacement matrix
      double N[3][3*10]; // shape matrix
      double dN[3][10]; // derivative of shape matrix
      //double dV = 0.0;
      for (int iq = 0; iq < 4; ++iq) {
        getdN(dN,GaussQuadrature::tet4[iq]);
        getJ(J,dN,ite);
        double detJ = fabs(DETERMINANT(J));
        getB(B,J,dN,ite);

        // temperature at quadrature point...
        getN(N,GaussQuadrature::tet4[iq]);
        double T_q = 0.0;
        FOR_J10 T_q += N[0][3*j]*T_no[noote[ite][j]];

	const int izn = ( (znote == NULL) ? 0 : znote[ite] );
        double nu_q;
	nuFuncVec[izn]->eval(&nu_q,&T_q,1);
        double alpha_q;
	alphaFuncVec[izn]->eval(&alpha_q,&T_q,1);

        D[0][0] = 1.0-nu_q;
        D[0][1] = nu_q;
        D[0][2] = nu_q;

        D[1][0] = nu_q;
        D[1][1] = 1.0-nu_q;
        D[1][2] = nu_q;

        D[2][0] = nu_q;
        D[2][1] = nu_q;
        D[2][2] = 1.0-nu_q;

        D[3][3] = (0.5-nu_q);
        D[4][4] = (0.5-nu_q);
        D[5][5] = (0.5-nu_q);

        // need to include these when we have loads other than thermal...
        //const double coeff = E/((1.0+nu)*(1.0-2.0*nu));
        //FOR_I6 FOR_J6 D[i][j] *= coeff;

        const double wgt = GaussQuadrature::tet4[iq][4]/6.0; 

        // stiffness matrix...
        for (int i = 0; i < 3*10; ++i)
          for (int j = 0; j < 3*10; ++j)
            for (int k = 0; k < 6; ++k)
              for (int l = 0; l < 6; ++l)
                Ke[i][j] += B[k][i]*D[k][l]*B[l][j]*detJ*wgt;

        // virtual thermal load...
        const double dT = T_q-T0;
        for (int i = 0; i < 3*10; ++i) {
          for (int k = 0; k < 6; ++k)
            for (int l = 0; l < 3; ++l) // [1 1 1 0 0 0]^T
              fe[i] += B[k][i]*D[k][l]*alpha_q*dT*detJ*wgt;
        }
        //dV += detJ*wgt;
      }
      //cout << ": " << dV << endl;

      // apply pin constraints...
      for (int i = 0; i < 3*10; ++i) {
        if (d_no[noote[ite][i/3]][i%3] == HUGE_VAL) {
          fe[i] = 0.0;
          for (int j = 0; j < 3*10; ++j) {
            if (i == j) {
              Ke[i][j] = 1.0;
            }
            else {
              Ke[i][j] = 0.0;
              Ke[j][i] = 0.0;
            }
          }
        }
      }

      /*
      {
        double dNTri[2][6];
        for (int it = 0; it < 4; ++it) {
          double dS = 0.0;
          for (int iq = 0; iq < 3; ++iq) {
            const double wgt = 0.5*GaussQuadrature::tri3[iq][3];
            getdNTri(dNTri,GaussQuadrature::tri3[iq]);
            dS += wgt*getdSTri(dNTri,ite,it);
          }
          cout << ": " << dS << endl;
        }
      }
      */

    }

    void getNoStress() {

      double D[6][6]; // constiutive matrix...
      FOR_I6 FOR_J6 D[i][j] = 0.0;
      double J[3][3]; // Jacobian matrix
      double B[6][3*10]; // strain-displacement matrix
      double N[3][3*10]; // shape matrix
      double dN[3][10]; // derivative of shape matrix
      double Se[6]; // symmetric stress...

      assert(s_no_d == NULL); s_no_d = new double[nno][3]; 
      assert(s_no_od == NULL); s_no_od = new double[nno][3]; 
      double *wgt_no = new double[nno];
      FOR_INO {
        FOR_I3 {
          s_no_d[ino][i] = 0.0;
          s_no_od[ino][i] = 0.0;
        }
        wgt_no[ino] = 0.0;
      }

      FOR_ITE {

        // get element weighted stress and weight...

        FOR_I6 Se[i] = 0.0;
        double Ve = 0.0;
        for (int iq = 0; iq < 4; ++iq) {
          getdN(dN,GaussQuadrature::tet4[iq]);
          getJ(J,dN,ite);
          double detJ = fabs(DETERMINANT(J));
          getB(B,J,dN,ite);

          // temperature at quadrature point...
          getN(N,GaussQuadrature::tet4[iq]);
          double T_q = 0.0;
          FOR_J10 T_q += N[0][3*j]*T_no[noote[ite][j]];

          const int izn = ( (znote == NULL) ? 0 : znote[ite] );
          double nu_q;
          nuFuncVec[izn]->eval(&nu_q,&T_q,1);
          double alpha_q;
          alphaFuncVec[izn]->eval(&alpha_q,&T_q,1);
          double E_q;
          EFuncVec[izn]->eval(&E_q,&T_q,1);

          D[0][0] = 1.0-nu_q;
          D[0][1] = nu_q;
          D[0][2] = nu_q;

          D[1][0] = nu_q;
          D[1][1] = 1.0-nu_q;
          D[1][2] = nu_q;

          D[2][0] = nu_q;
          D[2][1] = nu_q;
          D[2][2] = 1.0-nu_q;

          D[3][3] = (0.5-nu_q);
          D[4][4] = (0.5-nu_q);
          D[5][5] = (0.5-nu_q);

          const double coeff = E_q/((1.0+nu_q)*(1.0-2.0*nu_q));
          FOR_I6 FOR_J6 D[i][j] *= coeff;

          const double wgt = GaussQuadrature::tet4[iq][4]/6.0;

          // integrated stress...
          for (int k = 0; k < 6; ++k) 
            for (int l = 0; l < 6; ++l) 
              FOR_I10 FOR_J3 Se[k] += D[k][l]*B[l][i*3+j]*d_no[noote[ite][i]][j]*detJ*wgt;

          // integrated volume...
          Ve += detJ*wgt;

        }

        // add contribution to all nodes...
        FOR_I10 {
          const int ino = noote[ite][i];
          if (ino < nno) {
            FOR_J3 {
              s_no_d[ino][j] += Se[j];
              s_no_od[ino][j] += Se[j+3];
            }
            wgt_no[ino] += Ve;
          }
        }

      }

      MiscUtils::dumpRange(wgt_no,nno,"wgt_no");

      // normalize...
      FOR_INO {
        assert(wgt_no[ino] > 0.0);
        FOR_I3 {
          s_no_d[ino][i] /= wgt_no[ino];
          s_no_od[ino][i] /= wgt_no[ino];
        }
      }
      delete[] wgt_no;

      MiscUtils::dumpRange(s_no_d,nno,"s_no_d");
      MiscUtils::dumpRange(s_no_od,nno,"s_no_od");
    }

    void writeTecplot(const string& filename) { 

      // reduce active node set to corners...
      int * flag_no = new int[nno_g];
      FOR_INO flag_no[ino] = -1;
      FOR_ITE {
        FOR_I4 {
          const int ino = noote[ite][i];
          if (ino < nno) 
            flag_no[ino] = 0;
        }
      }
      int nno_flag = 0; 
      FOR_INO {
        if (flag_no[ino] == 0)
          flag_no[ino] = nno_flag++;
        else
          assert(flag_no[ino] == -1);
      }

      // handshake version of ascii tecplot write...

      int my_buf[2] = {nno_flag, nte_ib};
      int buf[2];
      MPI_Reduce(my_buf,buf,2,MPI_INT,MPI_SUM,0,mpi_comm);

      int offset;
      MPI_Scan(&nno_flag,&offset,1,MPI_INT,MPI_SUM,mpi_comm);
      offset -= nno_flag;

      FOR_INO {
        if (flag_no[ino] >= 0)
          flag_no[ino] += offset;
        else
          assert(flag_no[ino] == -1);
      }
      updateNoData(flag_no);

      // ====================
      // header and nodes first...
      // ====================

      FILE * fp;
      if ( mpi_rank == 0 ) {
        fp = fopen(filename.c_str(),"w");
        assert(fp != NULL);
        fprintf(fp,"TITLE = \"%s\"\n",filename.c_str());
        fprintf(fp,"VARIABLES = \"X\"\n");
        fprintf(fp,"\"Y\"\n");
        fprintf(fp,"\"Z\"\n");
        if (s_no_d) {
          assert(s_no_od);
          assert(d_no);
          fprintf(fp,"\"X0\"\n");
          fprintf(fp,"\"Y0\"\n");
          fprintf(fp,"\"Z0\"\n");
          fprintf(fp,"\"T\"\n");
          fprintf(fp,"\"SXX\"\n");
          fprintf(fp,"\"SYY\"\n");
          fprintf(fp,"\"SZZ\"\n");
          fprintf(fp,"\"SXY\"\n");
          fprintf(fp,"\"SXZ\"\n");
          fprintf(fp,"\"SYZ\"\n");
        }
        fprintf(fp,"ZONE T=\"%s\"\n",filename.c_str());
        fprintf(fp,"N=%d, E=%d, F=FEPOINT, ET=TETRAHEDRON\n",buf[0],buf[1]);
        cout << " > writing: " << buf[0] << " nodes, " << buf[1] << " tets..." << endl;
      }
      else {
        int dummy;
        MPI_Recv(&dummy,1,MPI_INT,mpi_rank-1,1234,mpi_comm,MPI_STATUS_IGNORE);
        fp = fopen(filename.c_str(),"a");
        assert(fp != NULL);
      }

      for (int ino = 0; ino < nno; ++ino) {
        if (flag_no[ino] >= 0) {
          if (s_no_d) {
            assert(s_no_od);
            fprintf(fp,"%18.15le %18.15le %18.15le %18.15le %18.15le %18.15le %18.15le %18.15le %18.15le %18.15le %18.15le %18.15le %18.15le\n",
                x_no[ino][0],x_no[ino][1],x_no[ino][2],x_no[ino][0]-d_no[ino][0],x_no[ino][1]-d_no[ino][1],x_no[ino][2]-d_no[ino][2],
                T_no[ino],s_no_d[ino][0],s_no_d[ino][1],s_no_d[ino][2],s_no_od[ino][0],s_no_od[ino][1],s_no_od[ino][2]);
          }
          else {
            fprintf(fp,"%18.15le %18.15le %18.15le\n",x_no[ino][0],x_no[ino][1],x_no[ino][2]);
          }
        }
      }

      fclose(fp);

      if ( mpi_rank < mpi_size-1 ) {
        int dummy = 1;
        MPI_Send(&dummy,1,MPI_INT,mpi_rank+1,1234,mpi_comm);
      }
      MPI_Barrier(mpi_comm);

      // ====================
      // then noote...
      // ====================

      if ( mpi_rank == 0 ) {
        fp = fopen(filename.c_str(),"a");
        assert(fp != NULL);
      }
      else {
        int dummy;
        MPI_Recv(&dummy,1,MPI_INT,mpi_rank-1,1235,mpi_comm,MPI_STATUS_IGNORE);
        fp = fopen(filename.c_str(),"a");
        assert(fp != NULL);
      }

      for (int ite = 0; ite < nte_ib; ++ite) {
        assert(flag_no[noote[ite][0]] >= 0);
        assert(flag_no[noote[ite][1]] >= 0);
        assert(flag_no[noote[ite][2]] >= 0);
        assert(flag_no[noote[ite][3]] >= 0);
        fprintf(fp,"%d %d %d %d\n",
            flag_no[noote[ite][0]]+1,
            flag_no[noote[ite][1]]+1,
            flag_no[noote[ite][2]]+1,
            flag_no[noote[ite][3]]+1);
      }
      delete[] flag_no;

      fclose(fp);

      if ( mpi_rank < mpi_size-1 ) {
        int dummy = 1;
        MPI_Send(&dummy,1,MPI_INT,mpi_rank+1,1235,mpi_comm);
      }
      MPI_Barrier(mpi_comm);

    }

    void writeTecplotByZone(const string prefix) {

      // handshake version of ascii tecplot write by zone...

      int * flag_no = new int[nno_g];
      char filename[128];
      for (int izn = 0; izn < nzn; ++izn) {

        if (nzn > 1) 
          sprintf(filename,"%s.%d.dat",prefix.c_str(),izn);
        else
          sprintf(filename,"%s.dat",prefix.c_str());

        if (mpi_rank == 0)
          cout << " > writing zone: " << izn << " in tecplot ascii format..." << endl;

        // reduce node set to corner nodes for this zone...
        int my_buf[2] = {0,0}; // nno, nte
        FOR_INO flag_no[ino] = -1; 
        FOR_ITE {
          if ( izn == ((znote == NULL) ? 0 : znote[ite]) ) {
            if (ite < nte_ib) 
              ++my_buf[1]; 
            FOR_I4 { 
              const int ino = noote[ite][i]; 
              if (ino < nno) {
                if (flag_no[ino] == -1) {
                  ++my_buf[0];
                  flag_no[ino] = 0;
                }
              }
            }
          }
        }
        int buf[2];
        MPI_Reduce(my_buf,buf,2,MPI_INT,MPI_SUM,0,mpi_comm);

        int offset;
        MPI_Scan(&my_buf[0],&offset,1,MPI_INT,MPI_SUM,mpi_comm);
        offset -= my_buf[0];

        int nno_check = 0;
        FOR_INO {
          if (flag_no[ino] == 0)
            flag_no[ino] = offset+nno_check++;
          else
            assert(flag_no[ino] == -1);
        }
        assert(nno_check == my_buf[0]);
        updateNoData(flag_no);

        // ====================
        // header and nodes first...
        // ====================

        FILE * fp;
        if ( mpi_rank == 0 ) {
          fp = fopen(filename,"w");
          assert(fp != NULL);
          fprintf(fp,"TITLE = \"%s\"\n",filename);
          fprintf(fp,"VARIABLES = \"X\"\n");
          fprintf(fp,"\"Y\"\n");
          fprintf(fp,"\"Z\"\n");
          if (s_no_d) {
            assert(s_no_od);
            assert(d_no);
            fprintf(fp,"\"X0\"\n");
            fprintf(fp,"\"Y0\"\n");
            fprintf(fp,"\"Z0\"\n");
            fprintf(fp,"\"T\"\n");
            fprintf(fp,"\"SXX\"\n");
            fprintf(fp,"\"SYY\"\n");
            fprintf(fp,"\"SZZ\"\n");
            fprintf(fp,"\"SXY\"\n");
            fprintf(fp,"\"SXZ\"\n");
            fprintf(fp,"\"SYZ\"\n");
          }
          fprintf(fp,"ZONE T=\"%s\"\n",filename);
          fprintf(fp,"N=%d, E=%d, F=FEPOINT, ET=TETRAHEDRON\n",buf[0],buf[1]);
          cout << " > writing: " << buf[0] << " nodes, " << buf[1] << " tets..." << endl;
        }
        else {
          int dummy;
          MPI_Recv(&dummy,1,MPI_INT,mpi_rank-1,1234,mpi_comm,MPI_STATUS_IGNORE);
          fp = fopen(filename,"a");
          assert(fp != NULL);
        }

        FOR_INO {
          if (flag_no[ino] >= 0) {
            if (s_no_d) {
              assert(s_no_od);
              fprintf(fp,"%18.15le %18.15le %18.15le %18.15le %18.15le %18.15le %18.15le %18.15le %18.15le %18.15le %18.15le %18.15le %18.15le\n",
                  x_no[ino][0],x_no[ino][1],x_no[ino][2],x_no[ino][0]-d_no[ino][0],x_no[ino][1]-d_no[ino][1],x_no[ino][2]-d_no[ino][2],
                  T_no[ino],s_no_d[ino][0],s_no_d[ino][1],s_no_d[ino][2],s_no_od[ino][0],s_no_od[ino][1],s_no_od[ino][2]);
            }
            else {
              fprintf(fp,"%18.15le %18.15le %18.15le\n",x_no[ino][0],x_no[ino][1],x_no[ino][2]);
            }
          }
        }

        fclose(fp);

        if ( mpi_rank < mpi_size-1 ) {
          int dummy = 1;
          MPI_Send(&dummy,1,MPI_INT,mpi_rank+1,1234,mpi_comm);
        }
        MPI_Barrier(mpi_comm);

        // ====================
        // then noote...
        // ====================

        if ( mpi_rank == 0 ) {
          fp = fopen(filename,"a");
          assert(fp != NULL);
        }
        else {
          int dummy;
          MPI_Recv(&dummy,1,MPI_INT,mpi_rank-1,1235,mpi_comm,MPI_STATUS_IGNORE);
          fp = fopen(filename,"a");
          assert(fp != NULL);
        }

        for (int ite = 0; ite < nte_ib; ++ite) {
          if ( izn == ((znote == NULL) ? 0 : znote[ite]) ) {
            FOR_I4 assert(flag_no[noote[ite][i]] >= 0); 
            // tecplot uses 1 indexing
            fprintf(fp,"%d %d %d %d\n",
                flag_no[noote[ite][0]]+1, 
                flag_no[noote[ite][1]]+1,
                flag_no[noote[ite][2]]+1,
                flag_no[noote[ite][3]]+1);
          }
        }

        fclose(fp);

        if ( mpi_rank < mpi_size-1 ) {
          int dummy = 1;
          MPI_Send(&dummy,1,MPI_INT,mpi_rank+1,1235,mpi_comm);
        }
        MPI_Barrier(mpi_comm);

      }
      delete[] flag_no;
    }

    void buildPrcomm(vector<CvPrcomm>& cvPrcommVec,const int ncv,const int ncv_g,uint8 * rbi_g) {

      // this routine uses the names "cv, ncv, icv, cvPrcommVec", but works for nodes, etc. or any other data
      // topology where paired communicators need to be built to update ghost data...

      if (mpi_rank == 0) cout << "buildPrcomm()" << endl;

      // ------------------------------------------------
      // build the paired communicator -- this should be
      // symmetric in terms of rank, but possibly not count...
      // ------------------------------------------------

      assert(cvPrcommVec.empty());

      int rank_current = -1;
      int bits_current = -1;
      CvPrcomm * prcomm = NULL;
      for (int icv = ncv;  icv < ncv_g; ++icv) {
        int rank,bits,index;
        BitUtils::unpackRankBitsIndex(rank,bits,index,rbi_g[icv-ncv]);
        if (rank > rank_current) {
          // we just switched ranks. If there was a prcomm, then
          // complete the size...
          if (prcomm)
            prcomm->unpack_size = icv-prcomm->unpack_offset;
          rank_current = rank;
          bits_current = -1;
          cvPrcommVec.push_back(CvPrcomm());
          prcomm = &cvPrcommVec.back();
          prcomm->rank = rank;
          prcomm->unpack_offset = icv;
        }
        else {
          assert(rank_current == rank);
          assert(prcomm);
        }
        // just tests monotoncity of the bits...
        if (bits > bits_current) {
          bits_current = bits;
        }
        else {
          assert(bits_current == bits);
        }
      }
      // we just finished. If there was a prcomm, then
      // complete the size...
      if (prcomm) prcomm->unpack_size = ncv_g-prcomm->unpack_offset;


      // finally, we need to send/recv the indices and bits to the pack side
      // and build the packVecs. Note that these are not necessarily symmetric by construction
      // because internal faces may have been removed based on tolerance on one vd,
      // and not on the other...

      MPI_Request * sendRequestArray = new MPI_Request[cvPrcommVec.size()];
      MPI_Request * recvRequestArray = new MPI_Request[cvPrcommVec.size()];

      for (int ii = 0, ii_end=cvPrcommVec.size(); ii < ii_end; ++ii) {

        // post irecv...
        MPI_Irecv(&(cvPrcommVec[ii].pack_size),1,MPI_INT,cvPrcommVec[ii].rank,12345,mpi_comm,&(recvRequestArray[ii]));

        // and the send...
        MPI_Issend(&(cvPrcommVec[ii].unpack_size),1,MPI_INT,cvPrcommVec[ii].rank,12345,mpi_comm,&(sendRequestArray[ii]));

      }

      MPI_Waitall(cvPrcommVec.size(),sendRequestArray,MPI_STATUSES_IGNORE);
      MPI_Waitall(cvPrcommVec.size(),recvRequestArray,MPI_STATUSES_IGNORE);

      // now send from the unpack side to the pack side...

      int pack_size = 0;
      for (int ii = 0, ii_end=cvPrcommVec.size(); ii < ii_end; ++ii)
        pack_size += cvPrcommVec[ii].pack_size;

      uint8 * recv_rbi_g = new uint8[pack_size];

      pack_size = 0;
      for (int ii = 0, ii_end=cvPrcommVec.size(); ii < ii_end; ++ii) {

        // post irecv...
        MPI_Irecv(recv_rbi_g+pack_size,cvPrcommVec[ii].pack_size,MPI_UINT8,
            cvPrcommVec[ii].rank,12346,mpi_comm,&(recvRequestArray[ii]));
        pack_size += cvPrcommVec[ii].pack_size;

        // and the send...
        MPI_Issend(rbi_g+(cvPrcommVec[ii].unpack_offset-ncv),cvPrcommVec[ii].unpack_size,MPI_UINT8,
            cvPrcommVec[ii].rank,12346,mpi_comm,&(sendRequestArray[ii]));

      }

      MPI_Waitall(cvPrcommVec.size(),sendRequestArray,MPI_STATUSES_IGNORE);
      MPI_Waitall(cvPrcommVec.size(),recvRequestArray,MPI_STATUSES_IGNORE);

      delete[] sendRequestArray;
      delete[] recvRequestArray;

      // now build the packVec and periodicity...

      double R[9], t[3];
      pack_size = 0;
      for (int ii = 0, ii_end=cvPrcommVec.size(); ii < ii_end; ++ii) {
        assert(cvPrcommVec[ii].packVec.empty());
        cvPrcommVec[ii].packVec.resize(cvPrcommVec[ii].pack_size);
        CvPrcomm::Transform * transform = NULL;
        int bits_current = -1;
        for (int i = 0; i < cvPrcommVec[ii].pack_size; ++i) {
          int rank,bits,index;
          BitUtils::unpackRankBitsIndex(rank,bits,index,recv_rbi_g[pack_size+i]);
          assert(rank == mpi_rank);
          if (bits > bits_current) {
            // bits are about to change, so complete any translate and rotate...
            if (transform) transform->end = i;
            // look for new rotations/translations...
            const bool has_R = false; // getPeriodicR(R,bits); <- TODO
            const bool has_t = false; // getPeriodicT(t,bits);
            cvPrcommVec[ii].transformVec.push_back(CvPrcomm::Transform(has_R,R,has_t,t,bits,i));
            transform = &(cvPrcommVec[ii].transformVec.back());
            bits_current = bits;
          }
          else {
            assert(bits_current == bits);
          }
          assert((index >= 0)&&(index < ncv));
          cvPrcommVec[ii].packVec[i] = index;
        }
        if (transform) transform->end = cvPrcommVec[ii].pack_size;
        pack_size += cvPrcommVec[ii].pack_size;
      }

      delete[] recv_rbi_g;

    }

    void reorderPrcommPackVec(vector<CvPrcomm>& cvPrcommVec,const int * const order,const int ncv) {

      for (int ii = 0; ii < cvPrcommVec.size(); ++ii) {
        for (int i = 0, limit = cvPrcommVec[ii].packVec.size(); i < limit; ++i) {
          const int icv = cvPrcommVec[ii].packVec[i];
          assert((icv >= 0)&&(icv < ncv));
          cvPrcommVec[ii].packVec[i] = order[icv];
        }
      }

    }

    void buildNoono() {

      if (mpi_rank == 0) cout << "Fea::buildNoono()..." << endl;

      assert(noote != NULL);
      assert(noono_i == NULL);
      assert(noono_v == NULL);

      // build teono_i/v...

      int * teono_i = new int[nno+1];
      FOR_INO teono_i[ino+1] = 0;
      FOR_ITE {
        for (int i = 0; i < 10; ++i) {
          const int ino = noote[ite][i];
          if (ino < nno)
            ++teono_i[ino+1];
        }
      }
      teono_i[0] = 0;
      FOR_INO teono_i[ino+1] += teono_i[ino];
      const int teono_s = teono_i[nno];
      int * teono_v = new int[teono_s];
      FOR_ITE {
        for (int i = 0; i < 10; ++i) {
          const int ino = noote[ite][i];
          if (ino < nno)
            teono_v[teono_i[ino]++] = ite;
        }
      }
      for (int ino = nno-1; ino > 0; --ino)
        teono_i[ino] = teono_i[ino-1];
      teono_i[0] = 0;

      // now use that to build noono_i/v...

      noono_i = new int[nno+1];
      FOR_INO noono_i[ino+1] = 0;

      int * no_flag = new int[nno_g];
      for (int ino = 0; ino < nno_g; ++ino) no_flag[ino] = -1;

      FOR_INO {
        ++noono_i[ino+1]; // diagonal...
        no_flag[ino] = ino;
        for (int ton = teono_i[ino]; ton != teono_i[ino+1]; ++ton) {
          const int ite = teono_v[ton];
          for (int i = 0; i < 10; ++i) {
            const int ino_nbr = noote[ite][i];
            if (no_flag[ino_nbr] != ino) {
              ++noono_i[ino+1];
              no_flag[ino_nbr] = ino;
            }
          }
        }
      }
      noono_i[0] = 0;
      FOR_INO noono_i[ino+1] += noono_i[ino];
      const int noono_s = noono_i[nno];
      noono_v = new int[noono_s];

      FOR_INO {
        noono_v[noono_i[ino]++] = ino; // diagonal...
        no_flag[ino] = ino+nno; // nno offset trick
        for (int ton = teono_i[ino]; ton != teono_i[ino+1]; ++ton) {
          const int ite = teono_v[ton];
          for (int i = 0; i < 10; ++i) {
            const int ino_nbr = noote[ite][i];
            if (no_flag[ino_nbr] != ino+nno) {
              noono_v[noono_i[ino]++] = ino_nbr;
              no_flag[ino_nbr] = ino+nno;
            }
          }
        }
      }

      for (int ino = nno-1; ino > 0; --ino)
        noono_i[ino] = noono_i[ino-1];
      noono_i[0] = 0;

      delete[] no_flag;
      delete[] teono_i;
      delete[] teono_v;

      // checks...

      // [0:nno_i) nodes only have [0:nno) nbrs...

      for (int ino = 0; ino < nno_i; ++ino) {
        for (int non = noono_i[ino]; non != noono_i[ino+1]; ++non) {
          const int ino_nbr = noono_v[ino];
          assert((ino_nbr >= 0)&&(ino_nbr < nno));
        }
      }

    }

    void redistReorder() {
      if (mpi_rank == 0) 
        cout << "Fea::redistReorder()..." << endl;

      assert(nno_global == 0);
      nno_global = nno;
      MPI_Bcast(&nno_global,1,MPI_INT8,0,mpi_comm);
      assert(noora_striped == NULL);
      MiscUtils::calcThresholdDist(noora_striped,nno_global,mpi_size,DIST_THRESHOLD);

      int8 * noora = NULL;
      MiscUtils::buildXora(noora,nno);

      assert(noora[0] == 0);
      assert(noora[mpi_rank+1]-noora[mpi_rank] == nno);

      int8 *ino_global2 = new int8[nno];
      FOR_INO ino_global2[ino] = noora[mpi_rank] + ino; // global no index -- need this

      int nno_orig = nno;
      repartXcvPadt(x_no,ino_global2,nno,mpi_comm);

      // x_no should be balanced now. This reorder puts nearby data closer in memory...
      reorderXcv(x_no,ino_global2,nno); // note: new nno

      DistributedDataExchanger* dde = new DistributedDataExchanger(ino_global2,nno,noora);
      delete[] ino_global2; ino_global2 = NULL;

      // rbi == "rank-bits-index"...
      uint8* rbi = new uint8[nno];
      FOR_INO rbi[ino] = BitUtils::packRankBitsIndex(mpi_rank,0,ino);

      uint8 * rbi_orig = new uint8[nno_orig];
      dde->push(rbi_orig,rbi);
      delete[] rbi;
      delete dde; dde = NULL;

      // we now have the rank and index to send all nno_orig stuff in rbi_orig...
      // now pull the rbi_orig out to the tets...

      assert(ino_global2 == NULL); ino_global2 = new int8[nte*10];
      FOR_ITE FOR_I10 ino_global2[ite*10+i] = noote[ite][i];
      delete[] noote; noote = NULL;

      assert(dde == NULL); dde = new DistributedDataExchanger(ino_global2,nte*10,noora);
      delete[] ino_global2; ino_global2 = NULL;
      delete[] noora; // reused by nst below if present

      uint8 (*rbiote)[10] = new uint8[nte][10];
      dde->pull((uint8*)rbiote,rbi_orig);
      delete dde; dde = NULL;

      // and send the tets as 10 rbi's in the final distribution. This
      // also set the ghost data. We also use this exchange to send the
      // tet volume zone, it present. Otherwise, we send zero...

      if (send_count == NULL) send_count = new int[mpi_size];
      FOR_RANK send_count[rank] = 0;

      // here we temporarily use recv_count as a flag...
      if (recv_count == NULL) recv_count = new int[mpi_size];
      FOR_RANK recv_count[rank] = -1;
      FOR_ITE {
        FOR_I10 {
          int rank,bits,index;
          BitUtils::unpackRankBitsIndex(rank,bits,index,rbiote[ite][i]);
          if (recv_count[rank] != ite) {
            send_count[rank] += 11;
            recv_count[rank] = ite;
          }
        }
      }

      if (send_disp == NULL) send_disp = new int[mpi_size];
      send_disp[0] = 0;
      for (int rank = 1; rank < mpi_size; ++rank)
        send_disp[rank] = send_disp[rank-1] + send_count[rank-1];
      int send_count_sum = send_disp[mpi_size-1] + send_count[mpi_size-1];

      uint8 * send_buf_uint8 = new uint8[send_count_sum];
      FOR_ITE {
        FOR_I10 {
          int rank,bits,index;
          BitUtils::unpackRankBitsIndex(rank,bits,index,rbiote[ite][i]);
          if (recv_count[rank] != nte+ite) { // use nte offset here to save resetting flag to -1 ;)
            FOR_J10 send_buf_uint8[send_disp[rank]+j] = rbiote[ite][j];
            send_buf_uint8[send_disp[rank]+10] = ((znote == NULL) ? 0 : znote[ite]);
            send_disp[rank] += 11;
            recv_count[rank] = nte+ite;
          }
        }
      }
      delete[] rbiote;
      delete[] znote; znote = NULL;

      // rewind...
      send_disp[0] = 0;
      for (int rank = 1; rank < mpi_size; ++rank)
        send_disp[rank] = send_disp[rank-1] + send_count[rank-1];

      MPI_Alltoall(send_count,1,MPI_INT,recv_count,1,MPI_INT,mpi_comm);

      if (recv_disp == NULL) recv_disp = new int[mpi_size];
      recv_disp[0] = 0;
      for (int rank = 1; rank < mpi_size; ++rank)
        recv_disp[rank] = recv_count[rank-1] + recv_disp[rank-1];
      int recv_count_sum = recv_disp[mpi_size-1] + recv_count[mpi_size-1];
      assert(recv_count_sum%11 == 0);

      uint8 * recv_buf_uint8 = new uint8[recv_count_sum];
      MPI_Alltoallv(send_buf_uint8,send_count,send_disp,MPI_UINT8,
          recv_buf_uint8,recv_count,recv_disp,MPI_UINT8,mpi_comm);
      delete[] send_buf_uint8; send_buf_uint8 = NULL;

      map<const uint8,int> rbiMap;
      nno_g = nno; // recall nno stores the new node distribution
      nte = recv_count_sum/11;
      int * flag_te = new int[nte];
      FOR_ITE flag_te[ite] = 0;
      assert(noote == NULL); noote = new int[nte][10];
      assert(znote == NULL); znote = new int[nte];
      for (int irecv = 0; irecv < recv_count_sum; irecv += 11) {
        const int ite = irecv/11;
        FOR_I10 {
          int rank,bits,index;
          BitUtils::unpackRankBitsIndex(rank,bits,index,recv_buf_uint8[irecv+i]);
          if ((rank == mpi_rank)&&(bits == 0)) {
            assert((index >= 0)&&(index < nno));
            noote[ite][i] = index;
          }
          else {
            // mark this tet as a boundary tet:
            // 2: a lower rank also has this tet
            // 1: we are the lowest rank with this tet...
            if (rank < mpi_rank) flag_te[ite] = 2;
            else flag_te[ite] = max(flag_te[ite],1);
            map<const uint8,int>::iterator iter = rbiMap.find(recv_buf_uint8[irecv+i]);
            if (iter == rbiMap.end()) {
              // this is the first time we have seen this node...
              rbiMap[recv_buf_uint8[irecv+i]] = noote[ite][i] = nno_g++;
            }
            else {
              // this node is already indexed...
              noote[ite][i] = iter->second;
            }
          }
        }
        znote[ite] = recv_buf_uint8[irecv+10];
      }
      delete[] recv_buf_uint8; recv_buf_uint8 = NULL;

      // reorder the tets in 3 ranges:
      // 0: internal,
      // 1: internal-boundary-owned,
      // 2: internal-boundary-owned-by-another-rank

      int counts[3] = { 0, 0, 0 };
      FOR_ITE ++counts[flag_te[ite]];

      int nte_old = nte;
      int (*noote_old)[10] = noote;
      noote = new int[nte][10];
      int *znote_old = znote;
      znote = new int[nte];
      nte_i = 0;
      nte_ib = counts[0];
      nte = counts[0]+counts[1];
      // also we will be flagging the nodes for local reorder...
      int * order = new int[nno];
      FOR_INO order[ino] = 0;
      for (int ite = 0; ite < nte_old; ++ite) {
        switch (flag_te[ite]) {
          case 0:
            FOR_I10 noote[nte_i][i] = noote_old[ite][i];
            znote[nte_i] = znote_old[ite];
            ++nte_i;
            break;
          case 1:
            FOR_I10 {
              noote[nte_ib][i] = noote_old[ite][i];
              if (noote_old[ite][i] < nno) {
                order[noote_old[ite][i]] = 1;
              }
            }
            znote[nte_ib] = znote_old[ite];
            ++nte_ib;
            break;
          case 2:
            FOR_I10 {
              noote[nte][i] = noote_old[ite][i];
              if (noote_old[ite][i] < nno) {
                order[noote_old[ite][i]] = 1;
              }
            }
            znote[nte] = znote_old[ite];
            ++nte;
            break;
          default:
            assert(0);
        }
      }
      assert(nte == nte_old);
      delete[] noote_old;
      delete[] znote_old;
      delete[] flag_te;

      // at this point, the active nodes (ino < nno) on the boundary are flagged with
      // order == 1, and order == 0 elsewhere...
      counts[0] = 0;
      counts[1] = 0;
      FOR_INO ++counts[order[ino]];

      int nno_old = nno;
      double (*x_no_old)[3] = x_no;
      x_no = new double[nno_g][3];
      nno_i = 0;
      nno = counts[0];
      for (int ino = 0; ino < nno_old; ++ino) {
        switch(order[ino]) {
          case 0:
            FOR_I3 x_no[nno_i][i] = x_no_old[ino][i];
            order[ino] = nno_i++;
            break;
          case 1:
            FOR_I3 x_no[nno][i] = x_no_old[ino][i];
            order[ino] = nno++;
            break;
          default:
            assert(0);
        }
      }
      assert(nno == nno_old);
      delete[] x_no_old;

      // reorder ghost data so it is contiguous in rank/b/index. This is exactly the current
      // order of the map...

      int * order_g = new int[nno_g-nno];
      assert(rbi_g == NULL); rbi_g = new uint8[nno_g-nno];
      int ino_g = nno;
      for (map<const uint8,int>::iterator iter = rbiMap.begin(); iter != rbiMap.end(); ++iter) {
        order_g[iter->second-nno] = ino_g;
        rbi_g[ino_g-nno] = iter->first;
        ++ino_g;
      }
      assert(ino_g == nno_g);
      rbiMap.clear();

      // change the ghost order in noote and noost...
      FOR_ITE {
        FOR_I10 {
          const int ino = noote[ite][i];
          if (ino < nno) 
            noote[ite][i] = order[ino];
          else 
            noote[ite][i] = order_g[ino-nno];
        }
      }
      delete[] order_g;

      // and finally use rbi_g to build the node communicator...

      buildPrcomm(noPrcommVec,nno,nno_g,rbi_g);
      reorderPrcommPackVec(noPrcommVec,order,nno);

      // reorder the index in rbi_g...
      {

        int* send_count = new int[mpi_size];
        int* send_disp = new int[mpi_size];
        FOR_RANK send_count[rank] = 0;
        int* send_buf = new int[nno_g-nno];
        for (int iter = 0; iter < 2; ++iter) {
          for (int ino = nno; ino < nno_g; ++ino) {
            int rank,bits,index;
            BitUtils::unpackRankBitsIndex(rank,bits,index,rbi_g[ino-nno]);
            if (iter == 0)
              ++send_count[rank];
            else
              send_buf[send_disp[rank]++] = index;
          }
          send_disp[0] = 0;
          for (int rank = 1; rank < mpi_size; ++rank)
            send_disp[rank] = send_count[rank-1] + send_disp[rank-1];
        }

        int * recv_count = new int[mpi_size];
        MPI_Alltoall(send_count,1,MPI_INT,recv_count,1,MPI_INT,mpi_comm);

        int * recv_disp = new int[mpi_size];
        recv_disp[0] = 0;
        for (int rank = 1; rank < mpi_size; ++rank)
          recv_disp[rank] = recv_count[rank-1] + recv_disp[rank-1];
        int recv_count_sum = recv_disp[mpi_size-1] + recv_count[mpi_size-1];

        int * recv_buf = new int[recv_count_sum];
        MPI_Alltoallv(send_buf,send_count,send_disp,MPI_INT,
            recv_buf,recv_count,recv_disp,MPI_INT,mpi_comm);

        for (int irecv = 0; irecv < recv_count_sum; ++irecv) {
          const int ino0 = recv_buf[irecv];
          recv_buf[irecv] = order[ino0];
        }

        MPI_Alltoallv(recv_buf,recv_count,recv_disp,MPI_INT,
            send_buf,send_count,send_disp,MPI_INT,mpi_comm);
        delete[] recv_buf;
        delete[] recv_count;
        delete[] recv_disp;

        for (int ino = nno; ino < nno_g; ++ino) {
          int rank,bits,index;
          BitUtils::unpackRankBitsIndex(rank,bits,index,rbi_g[ino-nno]);
          rbi_g[ino-nno] = BitUtils::packRankBitsIndex(rank,bits,send_buf[send_disp[rank]++]);
        }
        delete[] send_buf;
        delete[] send_count;
        delete[] send_disp;

      }

      // update x_no into ghost in prep for operator build...

      updateNoData(x_no);

      // finally, bring over T_no and ino_global from the original global index striping....

      FOR_RANK send_count[rank] = 0;
      for (int ino = 0; ino < nno_orig; ++ino) {
        int rank,bits,index;
        BitUtils::unpackRankBitsIndex(rank,bits,index,rbi_orig[ino]);
        ++send_count[rank];
      }
      send_disp[0] = 0;
      for (int rank = 1; rank < mpi_size; ++rank)
        send_disp[rank] = send_disp[rank-1] + send_count[rank-1];
      const int send_count_sum2 = send_disp[mpi_size-1] + send_count[mpi_size-1];
      int * send_buf_int = new int[send_count_sum2*2];
      double * send_buf_double = new double[send_count_sum2];
      for (int ino = 0; ino < nno_orig; ++ino) {
        int rank,bits,index;
        BitUtils::unpackRankBitsIndex(rank,bits,index,rbi_orig[ino]);
        assert(bits == 0);
        send_buf_int[send_disp[rank]*2  ] = index;
        send_buf_int[send_disp[rank]*2+1] = ino;
        send_buf_double[send_disp[rank]] = T_no[ino];
        ++send_disp[rank];
      }
      delete[] rbi_orig;
      if (T_no != NULL) delete[] T_no; T_no = NULL;
      // rewind...
      send_disp[0] = 0;
      for (int rank = 1; rank < mpi_size; ++rank)
        send_disp[rank] = send_disp[rank-1] + send_count[rank-1];
      // prepare recv side...
      MPI_Alltoall(send_count,1,MPI_INT,recv_count,1,MPI_INT,mpi_comm);
      recv_disp[0] = 0;
      for (int rank = 1; rank < mpi_size; ++rank)
        recv_disp[rank] = recv_count[rank-1] + recv_disp[rank-1];
      const int recv_count_sum2 = recv_disp[mpi_size-1] + recv_count[mpi_size-1];

      int * recv_buf_int = new int[2*recv_count_sum2];
      FOR_RANK {
        send_count[rank] *= 2;
        send_disp[rank] *= 2;
        recv_count[rank] *= 2;
        recv_disp[rank] *= 2;
      }
      MPI_Alltoallv(send_buf_int,send_count,send_disp,MPI_INT,
          recv_buf_int,recv_count,recv_disp,MPI_INT,mpi_comm);
      delete[] send_buf_int;
      FOR_RANK {
        send_count[rank] /= 2;
        send_disp[rank] /= 2;
        recv_count[rank] /= 2;
        recv_disp[rank] /= 2;
      }

      double * recv_buf_double = new double[recv_count_sum2];
      MPI_Alltoallv(send_buf_double,send_count,send_disp,MPI_DOUBLE,
          recv_buf_double,recv_count,recv_disp,MPI_DOUBLE,mpi_comm);
      delete[] send_buf_double;
      // now unpack into T_no...
      T_no = new double[nno_g];
      assert(ino_global == NULL); ino_global = new int8[nno]; // need ghosts?
      assert(noora_striped != NULL);
      for (int ino = 0; ino < nno_g; ++ino) T_no[ino] = HUGE_VAL;
      FOR_RANK {
        int irecv = recv_disp[rank];
        while (irecv < recv_disp[rank]+recv_count[rank]) {
          const int ino = recv_buf_int[irecv*2  ]; assert((ino >= 0)&&(ino < nno));
          const int ino_striped = recv_buf_int[irecv*2+1]; assert((ino >= 0)&&(ino < nno));
          const int ino_new = order[ino]; assert((ino_new >= 0)&&(ino_new < nno));
          assert(T_no[ino_new] == HUGE_VAL);
          T_no[ino_new] = recv_buf_double[irecv];
          ino_global[ino_new] = ino_striped + noora_striped[rank];
          ++irecv;
        }
      }
      delete[] recv_buf_int;
      delete[] recv_buf_double;
      delete[] order;
      assert(dde_striped == NULL); dde_striped = new DistributedDataExchanger(ino_global,nno,noora_striped);

      MiscUtils::dumpRange(&nno,1,"cht nodes");

      // check on tet ordering...

      for (int ite = 0; ite < nte_i; ++ite) {
        FOR_I10 {
          const int ino = noote[ite][i];
          assert((ino >= 0)&&(ino < nno));
        }
      }

      buildNoono();

      // user defined temperature profile for testing...

      if (Param * param = getParam("T_PROFILE")) {
        const string name = param->getString(0);
        if (name == "QUADRATIC") {
          if (mpi_rank == 0)
            cout << " > applying T = x*x..." << endl;
          for (int ino = 0; ino < nno_g; ++ino) 
            T_no[ino] = x_no[ino][0]*x_no[ino][0];
        }
        else if (name == "LINEAR") {
          if (mpi_rank == 0)
            cout << " > applying T = x..." << endl;
          for (int ino = 0; ino < nno_g; ++ino) 
            T_no[ino] = x_no[ino][0];
        }
        else if (name == "CONSTANT") {
          if (mpi_rank == 0)
            cout << " > applying T = 1..." << endl;
          for (int ino = 0; ino < nno_g; ++ino) 
            T_no[ino] = 1.0;
        }
        else {
          CERR(" > must specify T_PROFILE QUADRATIC, LINEAR, or CONSTANT...");
        }
      }
      updateNoData(T_no);
      MiscUtils::dumpRange(T_no,nno_g,"T_no");

      // build no adt and get global bounding box...
      
      {

        if (mpi_rank == 0) 
          cout << " > building x_no adt..." << endl;

        // put the local nodes in an adt...

        assert(adt == NULL);
        assert(bbox_adt == NULL);
        adt = new Adt<double>(nno,x_no,x_no);

        // the top leaf of this local adt stores the local bbox...
        double my_bbmin[3],my_bbmax[3];
        adt->getBbox(my_bbmin,my_bbmax);
        double (*bbmin)[3] = new double[mpi_size][3];
        double (*bbmax)[3] = new double[mpi_size][3];
        MPI_Allgather(my_bbmin,3,MPI_DOUBLE,bbmin,3,MPI_DOUBLE,mpi_comm);
        MPI_Allgather(my_bbmax,3,MPI_DOUBLE,bbmax,3,MPI_DOUBLE,mpi_comm);
        bbox_adt = new Adt<double>(mpi_size,bbmin,bbmax);
        delete[] bbmin;
        delete[] bbmax;
        bbox_adt->getBbox(bbminmax,bbminmax+3);

        if (mpi_rank == 0) {
          cout << " > bounding box: " <<
            bbminmax[0] << ":" << bbminmax[3] << " " <<
            bbminmax[1] << ":" << bbminmax[4] << " " <<
            bbminmax[2] << ":" << bbminmax[5] << endl;
        }

      }


      // user defined zone profile for testing...
      if (checkParam("SPLIT_ZONE")) {
        if (mpi_rank == 0)
          cout << " > setting zone to 1 x > bbmox mid point..." << endl;
        FOR_ITE {
          double x_avg[3] = {0.0,0.0,0.0};
          FOR_I4 FOR_J3 x_avg[j] += x_no[noote[ite][i]][j];
          FOR_I3 x_avg[i] *= 0.25;
          if (x_avg[0] > 0.5*(bbminmax[0]+bbminmax[3]))
            znote[ite] = 1;
        }
      }

      // =====================================================
      // load in the properties and build constitutive matrix...
      // =====================================================
      
      assert(EFuncVec.empty());
      assert(nuFuncVec.empty());
      assert(alphaFuncVec.empty());
      int my_nzn = 0;
      FOR_ITE {
        assert(znote[ite] >= 0);
        my_nzn = max(my_nzn,znote[ite]+1);
      }
      assert(nzn == 0);
      MPI_Allreduce(&my_nzn,&nzn,1,MPI_INT,MPI_MAX,mpi_comm);
      if (nzn == 1) {
        // for one zone only, look for
        // CHT.RHO, CHT.CP, CHT.K
        if (mpi_rank == 0) cout << " > CHT.E..." << endl;
        EFuncVec.push_back(processSimpleFunc(getParam("CHT.E")));
        if (mpi_rank == 0) cout << " > CHT.NU..." << endl;
        nuFuncVec.push_back(processSimpleFunc(getParam("CHT.NU")));
        if (mpi_rank == 0) cout << " > CHT.ALPHA..." << endl;
        alphaFuncVec.push_back(processSimpleFunc(getParam("CHT.ALPHA")));
      }
      else {
        // for multiple zones, use indexing...
        // CHT.0.E, CHT.0.NU, CHT.0.ALPHA
        // CHT.1.E, CHT.1.NU, CHT.1.ALPHA
        char name[128];
        for (int izn = 0; izn < nzn; ++izn) {
          sprintf(name,"CHT.%d.E",izn);
          if (mpi_rank == 0) cout << " > " << name << "..." << endl;
          EFuncVec.push_back(processSimpleFunc(getParam(name)));
          sprintf(name,"CHT.%d.NU",izn);
          if (mpi_rank == 0) cout << " > " << name << "..." << endl;
          nuFuncVec.push_back(processSimpleFunc(getParam(name)));
          sprintf(name,"CHT.%d.ALPHA",izn);
          if (mpi_rank == 0) cout << " > " << name << "..." << endl;
          alphaFuncVec.push_back(processSimpleFunc(getParam(name)));
        }
      }

      double my_dxminmax[2] = {HUGE_VAL,HUGE_VAL};
      for (int ite = 0; ite < nte; ++ite) {
        for (int i = 0; i < 4; ++i) {
          for (int j = i+1; j < 4; ++j) {
            const double dx = 0.5*DIST(x_no[noote[ite][i]],x_no[noote[ite][j]]);
            my_dxminmax[0] = min(dx,my_dxminmax[0]);
            my_dxminmax[1] = min(-dx,my_dxminmax[1]);
          }
        }
      }
      MPI_Allreduce(my_dxminmax,dxminmax,2,MPI_DOUBLE,MPI_MIN,mpi_comm);
      dxminmax[1] *= -1.0;

      if (mpi_rank == 0.0)
        cout << " > min/max dx0: " << dxminmax[0]<< "/" << dxminmax[1] << endl; 

    }

    class FacePin {
      public:
        string face; // xmin,ymin,zmin,xmax,ymax,zmax
        bool pin[3];
    };

    class PointPin {
      public:
        double x[3];
        bool pin[3];
    };

    void setConstraints() {
      if (mpi_rank == 0) 
        cout << "Fea::setConstraints()..." << endl;

      bool b_pin_dir[3] = {false,false,false};  // x,y,z
      vector<PointPin> ppVec;
      vector<FacePin> fpVec;

      FOR_PARAM_MATCHING("PIN") {
        int iarg = 0;
        while ( iarg < param->size()) {
          string token = param->getString(iarg++);
          if (token == "FACE") {
            fpVec.push_back(FacePin());
            string face = param->getString(iarg++);
            if ((face == "X0")||(face == "XMIN")) {
              fpVec.back().face = "XMIN";
            }
            else if ((face == "Y0")||(face == "YMIN")) {
              fpVec.back().face = "YMIN";
            }
            else if ((face == "Z0")||(face == "ZMIN")) {
              fpVec.back().face = "ZMIN";
            }
            else if ((face == "X1")||(face == "XMAX")) {
              fpVec.back().face = "XMAX";
            }
            else if ((face == "Y1")||(face == "YMAX")) {
              fpVec.back().face = "YMAX";
            }
            else if ((face == "Z1")||(face == "ZMAX")) {
              fpVec.back().face = "ZMAX";
            }
            else {
              CERR( " > unrecognized token in PIN FACE: " << face);
            }
            // get direction to pin...
            FOR_I3 fpVec.back().pin[i] = param->getBool(iarg++);
            if (mpi_rank == 0) 
              cout << " > pin nodes on face: " << fpVec.back().face << " in directions: " << COUT_VEC(fpVec.back().pin) << endl;
          } 
          else if ((token == "DIR")||(token == "DIRECTION")) {
            string dir = param->getString(iarg++);
            if ((dir == "0")||(dir == "X")) {
              b_pin_dir[0] = true;
              if (mpi_rank == 0) 
                cout << " > pin all nodes in x direction..." << endl;
            }
            else if ((dir == "1")||(dir == "Y")) {
              b_pin_dir[1] = true;
              if (mpi_rank == 0) 
                cout << " > pin all nodes in y direction..." << endl;
            }
            else if ((dir == "2")||(dir == "Z")) {
              b_pin_dir[2] = true;
              if (mpi_rank == 0) 
                cout << " > pin all nodes in z direction..." << endl;
            }
            else {
              CERR( " > unrecognized token in PIN DIRECTION: " << dir);
            }
          } 
          else if (token == "POINT") {
            ppVec.push_back(PointPin());
            FOR_I3 ppVec.back().x[i] = param->getDouble(iarg++);
            FOR_I3 ppVec.back().pin[i] = param->getBool(iarg++);
            if (mpi_rank == 0) 
              cout << " > pin node closest too point: " << COUT_VEC(ppVec.back().x) << " in directions: " << COUT_VEC(ppVec.back().pin) << endl;
          }
          else {
            CERR( " > unrecognized token in PIN : " << token);
          }
        }
      }

      // use displacement to flag/describe pins...

      assert(d_no == NULL); d_no = new double[nno_g][3];
      for (int ino = 0; ino < nno_g; ++ino) 
        FOR_I3 d_no[ino][i] = 0.0;

      // direction constraints... 

      for (int ino = 0; ino < nno; ++ino) 
        FOR_I3 if (b_pin_dir[i]) d_no[ino][i] = HUGE_VAL;

      // face constraints...

      const double tol = 0.1*dxminmax[0];
      for (int ifp = 0, nfp = fpVec.size(); ifp < nfp; ++ifp) {
        FOR_INO {
          if ( ((fpVec[ifp].face == "XMIN")&&(x_no[ino][0] < (bbminmax[0]+tol)))||
               ((fpVec[ifp].face == "YMIN")&&(x_no[ino][1] < (bbminmax[1]+tol)))||
               ((fpVec[ifp].face == "ZMIN")&&(x_no[ino][2] < (bbminmax[2]+tol)))||
               ((fpVec[ifp].face == "XMAX")&&(x_no[ino][0] > (bbminmax[4]-tol)))||
               ((fpVec[ifp].face == "YMAX")&&(x_no[ino][1] > (bbminmax[5]-tol)))||
               ((fpVec[ifp].face == "ZMAX")&&(x_no[ino][2] > (bbminmax[6]-tol))) ) {
            FOR_I3 if (fpVec[ifp].pin[i]) d_no[ino][i] = HUGE_VAL;
          }
        }
      }

      // point constraints...

      const int npp = ppVec.size();
      DoubleInt * my_di = new DoubleInt[npp]; // rank/distance
      int * ino_closest = new int[npp];  
      for (int ipp = 0; ipp < npp; ++ipp) {
        ino_closest[ipp] = -1;
        my_di[ipp].this_double = HUGE_VAL; 
        // decided against adt here, so that we gauruntee that a pin always gets a node...
        FOR_INO {
          const double d2 = DIST2(x_no[ino],ppVec[ipp].x);
          if ((ino_closest[ipp] == -1)||(d2 < my_di[ipp].this_double)) {
            ino_closest[ipp] = ino;
            my_di[ipp].this_double = d2;
          }
        }
        if (ino_closest[ipp] != -1) 
          my_di[ipp].this_int = mpi_rank;
        else
          my_di[ipp].this_int = -1; // set the int to -1 unless we found one...
      }
      DoubleInt * di = new DoubleInt[npp];
      MPI_Allreduce(my_di,di,npp,MPI_DOUBLE_INT,MPI_MINLOC,mpi_comm);
      delete[] my_di;

      for (int ipp = 0; ipp < npp; ++ipp) {
        if (di[ipp].this_int == mpi_rank) {
          assert(ino_closest[ipp] >= 0);
          FOR_I3 {
            if (ppVec[ipp].pin[i]) 
              d_no[ino_closest[ipp]][i] = HUGE_VAL;
          }
        }
      }
      ppVec.clear();
      delete[] di;
      delete[] ino_closest;

      updateNoDataStart(d_no);

      int8 my_nc = 0;
      FOR_INO FOR_I3 if (d_no[ino][i] == HUGE_VAL) ++my_nc;

      updateNoDataFinish(d_no);

      int8 nc;
      MPI_Allreduce(&my_nc,&nc,1,MPI_INT8,MPI_SUM,mpi_comm);
      if (nc < 6) {
        // should also check that they are not all collinear constraints...
        CERR( " > must add PIN FACE, PIN DIRECTION and PIN POINT constraints to remove rigid body motion...");
      }

    }

    void buildKdf() {
      if (mpi_rank == 0) 
        cout << "Fea::buildKdf()..." << endl;

      assert(noono_i);
      assert(noono_v);
      assert(d_no);

      assert(f_no == NULL); f_no = new double[nno][3];
      FOR_INO FOR_I3 f_no[ino][i] = 0.0;
      assert(K_no == NULL); K_no = new double[noono_i[nno]][3][3];
      for (int non = 0; non < noono_i[nno]; ++non) FOR_I3 FOR_J3 K_no[non][i][j] = 0.0;

      // build stiffness matrix and rhs..
      double Ke[3*10][3*10]; // element stiffness matrix
      double fe[3*10]; // element thermal load
      FOR_ITE {
        getKeAndFe(Ke,fe,ite);
        FOR_I10 {
          const int ino = noote[ite][i];
          if (ino < nno) {
            FOR_J3 f_no[ino][j] += fe[i*3+j];
            const int non = noono_i[ino]; assert(noono_v[non] == ino);
            FOR_J10 {
              const int ino_nbr = noote[ite][j];
              int non_nbr = non;
              while (noono_v[non_nbr] != ino_nbr)
                ++non_nbr;
              assert(non_nbr < noono_i[ino+1]); // make sure we found ino_nbr...
              FOR_K3 FOR_L3 K_no[non_nbr][k][l] += Ke[i*3+k][j*3+l];
            }
          }
        }
      }
      MiscUtils::dumpRange(f_no,nno,"f_no");

      // set initial guess to zero...
      for (int ino = 0; ino < nno_g; ++ino)
        FOR_I3 d_no[ino][i] = 0.0;

    }

    int solveNoCg(double (*phi)[3],const double (*const A)[3][3],const double (* const rhs)[3],
        const double zero,const int maxiter,const bool verbose) {

      // we need the following work arrays...
      double (*res)[3]      = new double[nno][3];
      double (*v)[3]        = new double[nno][3];
      double (*p)[3]        = new double[nno_g][3];
      double (*inv_diag)[3][3] = new double[nno][3][3];

      // initialize...
      FOR_INO {
        assert(noono_v[noono_i[ino]] == ino); // confirm diagonal first in noono_i/v CSR
        FOR_I3 assert(A[noono_i[ino]][i][i] != 0.0); // diag cannot be zero for diag-preconditioned cg.
        MiscUtils::invertMat(inv_diag[ino],A[noono_i[ino]]);
      }

      for (int ino = 0; ino < nno_g; ++ino) {
        FOR_I3 {
          assert(phi[ino][i] == 0.0); // assumes initial guess is 0 
          p[ino][i] = 0.0; 
        }
      }
      double rho = 1.0;

      // calculate the residual in rhs format...
      // because the initial guess is zero, this is just the rhs...
      FOR_INO FOR_I3 res[ino][i] = rhs[ino][i];

      // diagonal precon/compute normalized residual...
      double my_res_max = 0.0;
      FOR_INO {
        FOR_I3 {
          v[ino][i] = 0.0;
          FOR_J3 v[ino][i] += inv_diag[ino][i][j]*res[ino][j];
          my_res_max = max( my_res_max, fabs(v[ino][i]) );
        }
      }
      double res_max = 0.0;
      if (verbose) {
        MPI_Reduce(&my_res_max,&res_max,1,MPI_DOUBLE,MPI_MAX,0,mpi_comm);
        if ((mpi_rank == 0)&&verbose) cout << " > initial res_max: " << res_max << endl;
      }

      int iter = 0;
      int done = 0;
      while (done == 0) {

        ++iter;

        assert(rho != 0.0);
        const double rho_prev = rho;

        double my_rho = 0.0;
        FOR_INO FOR_I3 my_rho += res[ino][i]*v[ino][i];
        MPI_Allreduce(&my_rho,&rho,1,MPI_DOUBLE,MPI_SUM,mpi_comm);

        assert(rho_prev != 0.0);
        const double beta = rho/rho_prev;
        FOR_INO FOR_I3 p[ino][i] = v[ino][i] + beta*p[ino][i];

        // v = [Ap]{p}...
        // with some attempt to hide latency...

        updateNoDataStart(p);

        for (int ino = 0; ino < nno_i; ++ino) {
          FOR_I3 {
            v[ino][i] = 0.0;
            FOR_J3 v[ino][i] += A[noono_i[ino]][i][j]*p[ino][j];
          }
          for (int non = noono_i[ino]+1; non < noono_i[ino+1]; ++non) {
            const int ino_nbr = noono_v[non];
            assert(ino_nbr < nno);
            FOR_I3 FOR_J3 v[ino][i] += A[non][i][j]*p[ino_nbr][j];
          }
        }

        updateNoDataFinish(p);

        for (int ino = nno_i; ino < nno; ++ino) {
          FOR_I3 {
            v[ino][i] = 0.0;
            FOR_J3 v[ino][i] += A[noono_i[ino]][i][j]*p[ino][j];
          }
          for (int non = noono_i[ino]+1; non < noono_i[ino+1]; ++non) {
            const int ino_nbr = noono_v[non];
            FOR_I3 FOR_J3 v[ino][i] += A[non][i][j]*p[ino_nbr][j];
          }
        }

        double my_gamma = 0.0;
        FOR_INO FOR_I3 my_gamma += p[ino][i]*v[ino][i];
        double gamma;
        MPI_Allreduce(&my_gamma,&gamma,1,MPI_DOUBLE,MPI_SUM,mpi_comm);
        assert(gamma != 0.0);
        assert(gamma == gamma); // nan check

        const double alpha = rho/gamma;

        // check if we are done...
        if (iter%10 == 0) {

          for (int ino = 0; ino < nno_g; ++ino) {
            FOR_I3 {
              assert(p[ino][i] == p[ino][i]); // nan check
              phi[ino][i] += alpha*p[ino][i];
            }
          }

          // recompute the residual...
          FOR_INO {
            FOR_I3 {
              res[ino][i] = rhs[ino][i];
              FOR_J3 res[ino][i] -= A[noono_i[ino]][i][j]*phi[ino][j]; 
            }
            for (int non = noono_i[ino]+1; non < noono_i[ino+1]; ++non) {
              const int ino_nbr = noono_v[non];
              FOR_I3 FOR_J3 res[ino][i] -= A[non][i][j]*phi[ino_nbr][j];
            }
          }

          // compute the max (L-infinity) normalized residual...
          my_res_max = 0.0;
          FOR_INO {
            FOR_I3 {
              v[ino][i] = 0.0;
              FOR_J3 v[ino][i] += inv_diag[ino][i][j]*res[ino][j];
              my_res_max = max( my_res_max, fabs(v[ino][i]) );
            }
          }
          MPI_Reduce(&my_res_max,&res_max,1,MPI_DOUBLE,MPI_MAX,0,mpi_comm);
          if (mpi_rank == 0) {
            // only share the last half of the convergence behaviour...
            if (verbose || (iter > maxiter/2))
              cout << " > solveCgLocal iter " << iter << " res_max " << res_max << endl;
            if (res_max <= zero) {
              if (verbose) cout << " > Successfully converged error to " << res_max << endl;
              done = 1;
            }
            else if (iter > maxiter) {
              cout << "Warning: solveCgLocal did not converge after " << maxiter <<
                " iters, res_max: " << res_max << endl;
              done = 2;
            }
          }
          MPI_Bcast(&done,1,MPI_INT,0,mpi_comm);

        }
        else {

          // update full phi...
          for (int ino = 0; ino < nno_g; ++ino)
            FOR_I3 phi[ino][i] += alpha*p[ino][i];

          FOR_INO {
            // on the other iterations, use this approximation to update
            // the unreduced residual...
            FOR_I3 res[ino][i] -= alpha*v[ino][i];
            // still need to compute v, diag precon for next iteration...
            FOR_I3 {
              v[ino][i] = 0.0;
              FOR_J3 v[ino][i] += inv_diag[ino][i][j]*res[ino][j];
            }
          }

        }


      }

      delete[] res;
      delete[] v;
      delete[] p;
      delete[] inv_diag;

      // let the calling routine know if we were successful...
      return done;

    }

    void init() {
      if ( mpi_rank == 0) 
        cout << "Fea::init() " << endl;

      // noote,x_no,T_no...
      if (mpi_rank == 0) {
        if (Param * param = getParam("RESTART")) {
          const string mesh_filename = param->getString(0);
          // records for 4 node tets...
          int (*noote0)[4] = NULL;
          double (*x_no0)[3] = NULL;
          nno0 = -1;
          if ((mesh_filename.size() > 8)&&(mesh_filename.compare(mesh_filename.size()-8,8,".tet_bin") == 0)) {
            // simple binary tet file...
            FILE * fp = fopen(mesh_filename.c_str(),"rb");
            if (fp == NULL) {
              cout << "Error: cannot open cht tet_bin file: " << mesh_filename << endl;
              nno0 = -1;
            }
            else {
              bool b_zones = false;
              size_t nread;
              nread = fread(&nno0,sizeof(int),1,fp);
              assert(nread == 1);
              if (nno0 == -1) {
                b_zones = true;
                nread = fread(&nno0,sizeof(int),1,fp);
                assert(nread == 1);
              }
              nread = fread(&nte,sizeof(int),1,fp);
              assert(nread == 1);
              cout << " > cht tet_bin file: " << mesh_filename << " has nno: " << nno0 << ", nte: " << nte << " b_zones: " << b_zones << endl;
              assert(x_no0 == NULL);
              x_no0 = new double[nno0][3];
              nread = fread(x_no0,sizeof(double),nno0*3,fp);
              assert(nread == (3*nno0));
              assert(noote0 == NULL);
              noote0 = new int[nte][4];
              nread = fread(noote0,sizeof(int),nte*4,fp);
              assert(nread == (4*nte));
              assert(znote == NULL);
              if (b_zones) {
                znote = new int[nte];
                nread = fread(znote,sizeof(int),nte,fp);
                assert(nread == nte);
                set<int> zoneSet;
                for (int ite = 0; ite < nte; ++ite) {
                  const int izn = znote[ite];
                  zoneSet.insert(izn);
                }
                cout << " > cht tet_bin file has these volume zones: ";
                for (set<int>::iterator it = zoneSet.begin(); it != zoneSet.end(); ++it) {
                  cout << " " << *it;
                }
                cout << endl;
              }
              fclose(fp);
            }
          }
          else {
            // assume coming from fluent...
            FluentMsh * msh = new FluentMsh(mesh_filename);
            // for now, everything on rank 0...
            nno0 = msh->nno;
            assert(noote0 == NULL);
            msh->buildTets(nte,noote0);
            assert(x_no0 == NULL); x_no0 = new double[nno0][3];
            for (int ino = 0; ino < nno0; ++ino) 
              FOR_I3 x_no0[ino][i] = msh->x_no[ino][i];
            delete msh;
          }

          // add nodes for each edge...

          map<pair<int,int>,int> edgeMap;
          vector<pair<int,int> > nodeVec(4);
          nno = nno0;
          for (int ite = 0; ite < nte; ++ite) {
            for (int ino_local = 0; ino_local < 4; ++ino_local) 
              nodeVec[ino_local] = pair<int,int>(noote0[ite][ino_local],ino_local);
            sort(nodeVec.begin(),nodeVec.end());
            for (int ino_local0 = 0; ino_local0 < 4; ++ino_local0) {
              for (int ino_local1 = ino_local0+1; ino_local1 < 4; ++ino_local1) {
                pair<int,int> edge = pair<int,int>(nodeVec[ino_local0].first,nodeVec[ino_local1].first);
                map<pair<int,int>,int>::iterator it = edgeMap.find(edge);
                if (it == edgeMap.end()) 
                  edgeMap[edge] = nno++;
              }
            }
          }
          cout << " > expanding node set to include edge midpoints: " << nno0 << " -> " << nno << endl;

          assert(x_no == NULL);
          x_no = new double[nno][3];
          for (int ino = 0; ino < nno0; ++ino)
            FOR_I3 x_no[ino][i] = x_no0[ino][i];
          delete[] x_no0;

          // if the restart file has a thermal solution, then set it...

          assert(T_no == NULL);
          T_no = new double[nno];
          if (param->size() > 1) {
            assert(param->size() == 2);
            const string soln_filename = param->getString(1);
            FILE * fp = fopen(soln_filename.c_str(),"rb");
            if (fp == NULL) {
              cout << "Error: cannot open cht soln file: " << soln_filename << endl;
            }
            else {
              int8 nno_check = -1;
              size_t nread;
              nread = fread(&nno_check,sizeof(int8),1,fp);
              assert(nread == 1);
              if (nno_check != nno0) {
                cout << "Error: mesh nno does no match that of cht soln file: " << soln_filename << endl;
              }
              else {
                cout << " > read T_no from cht soln file: " << soln_filename << endl;
                assert(T_no);
                nread = fread(T_no,sizeof(double),nno0,fp);
                assert(nread == nno0);
              }
            }
          }
          else {
            for (int ino = 0; ino < nno0; ++ino) 
              T_no[ino] = T0;
          }

          assert(noote == NULL);
          noote = new int[nte][10];
          for (int ite = 0; ite < nte; ++ite) {
            FOR_I4 noote[ite][i]   = noote0[ite][i];
            FOR_I6 noote[ite][i+4] = -1;
          }
          delete[] noote0;
          for (int ite = 0; ite < nte; ++ite) {
            for (int ino_local = 0; ino_local < 4; ++ino_local) 
              nodeVec[ino_local] = pair<int,int>(noote[ite][ino_local],ino_local);
            sort(nodeVec.begin(),nodeVec.end());
            int nno_local = 4;
            for (int ino_local0 = 0; ino_local0 < 4; ++ino_local0) {
              for (int ino_local1 = ino_local0+1; ino_local1 < 4; ++ino_local1) {
                const int ino0 = nodeVec[ino_local0].first;
                const int ino1 = nodeVec[ino_local1].first;
                pair<int,int> edge = pair<int,int>(ino0,ino1);
                map<pair<int,int>,int>::iterator it = edgeMap.find(edge);
                assert(it != edgeMap.end());
                // take average node temperature as edge midpoint temperature
                T_no[it->second] = 0.5*(T_no[ino0]+T_no[ino1]);
                FOR_I3 x_no[it->second][i] = 0.5*(x_no[ino0][i]+x_no[ino1][i]);
                // sort edge midpoints to follow shape function definitions
                const int imin = min(nodeVec[ino_local0].second,nodeVec[ino_local1].second);
                const int imax = max(nodeVec[ino_local0].second,nodeVec[ino_local1].second);
                if ((imin == 0)&&(imax == 1))
                  noote[ite][4] = it->second;
                else if ((imin == 1)&&(imax == 2))
                  noote[ite][5] = it->second;
                else if ((imin == 0)&&(imax == 2))
                  noote[ite][6] = it->second;
                else if ((imin == 0)&&(imax == 3))
                  noote[ite][7] = it->second;
                else if ((imin == 1)&&(imax == 3))
                  noote[ite][8] = it->second;
                else if ((imin == 2)&&(imax == 3))
                  noote[ite][9] = it->second;
                else
                  assert(0);
                ++nno_local;
              }
            }
            FOR_I10 assert((noote[ite][i] >= 0)&&(noote[ite][i] < nno));
            assert(nno_local == 10);
          }
          edgeMap.clear();

        }
        else if (checkParam("SINGLE_TET")) {
          nno0 = 4;
          nno = 10;
          nte = 1;
          assert(x_no == NULL);
          x_no = new double[nno][3];
          x_no[0][0] = 0.0; x_no[0][1] = 0.0; x_no[0][2] = 0.0;
          x_no[1][0] = 1.0; x_no[1][1] = 0.0; x_no[1][2] = 0.0;
          x_no[2][0] = 0.0; x_no[2][1] = 1.0; x_no[2][2] = 0.0;
          x_no[3][0] = 0.0; x_no[3][1] = 0.0; x_no[3][2] = 1.0;
          x_no[4][0] = 0.5; x_no[4][1] = 0.0; x_no[4][2] = 0.0;
          x_no[5][0] = 0.5; x_no[5][1] = 0.5; x_no[5][2] = 0.0;
          x_no[6][0] = 0.0; x_no[6][1] = 0.5; x_no[6][2] = 0.0;
          x_no[7][0] = 0.0; x_no[7][1] = 0.0; x_no[7][2] = 0.5;
          x_no[8][0] = 0.5; x_no[8][1] = 0.0; x_no[8][2] = 0.5;
          x_no[9][0] = 0.0; x_no[9][1] = 0.5; x_no[9][2] = 0.5;
          assert(noote == NULL);
          noote = new int[1][10];
          for (int i = 0; i < 10; ++i) 
            noote[0][i] = i;
          assert(T_no == NULL);
          T_no = new double[nno];
          for (int ino = 0; ino < nno; ++ino) 
            T_no[ino] = T0;
        }
        else {
          cout << "Error: must specify RESTART mesh.tet_bin [soln.cht]" << endl;
        }


        // if table is english units and tet_bin/cht are in metric...

        if (checkParam("METRIC_TO_ENGLISH")) {
          for (int ino = 0; ino < nno; ++ino) {
            T_no[ino] = (T_no[ino]-273.15)*9.0/5.0+32.0; // K to F
            FOR_I3 x_no[ino][i] *= 39.3701; // m to in
          }
        }

      }

      // load balances the above striping
      
      redistReorder(); 

      // write out initial state...
     
      if (checkParam("WRITE_INITIAL_MESH"))
        writeTecplotByZone("mesh");

      setConstraints();

      buildKdf();

    }

    void run() {
      if ( mpi_rank == 0) 
        cout << "Fea::run() " << endl;

      // solve K*d=f
      const bool verbose = getBoolParam("VERBOSE",true);
      const int maxiter = getIntParam("MAXITER",10000);
      const double zero = getDoubleParam("ZERO",1E-6);
      solveNoCg(d_no,K_no,f_no,zero,maxiter,verbose);
      getNoStress();

      // update nodal positions (must happen after stress comp)...
      for (int ino = 0; ino < nno; ++ino) 
        FOR_I3 x_no[ino][i] += d_no[ino][i];

      MiscUtils::dumpRange(d_no,nno,"d_no");
      MiscUtils::dumpRange(x_no,nno,"x_no final");

      writeTecplotByZone("soln");
    }

};

int main(int argc, char* argv[]) { 

  try { 

    CTI_Init(argc,argv,"fea.in");

    { 

      Fea solver;
      solver.init();
      solver.run();

    }

    CTI_Finalize();

  }
  catch( int e) { 

    if ( e == 0) 
      CTI_Finalize();
    else 
      CTI_Abort();

  }
  catch(...) { 

    CTI_Abort();

  }

  return 0;
} 
