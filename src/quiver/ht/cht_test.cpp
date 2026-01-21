
#include "Common.hpp"
#include "CTI.hpp"
using namespace CTI;
#include "Utils.hpp"
#include "Adt.hpp"
#include "DistributedDataExchanger.hpp"
#include "Cht.hpp"

class ChtTest : public Cht { 
public:

  void init() { 

    const string tles_file = getStringParam("RESTART_TLES");
    read_tles(tles_file);

    // HACK initial condition

    assert(T_no == NULL); T_no = new double[nno];
    for (int ino = 0; ino < nno; ++ino) 
      T_no[ino] = x_no[ino][0]+x_no[ino][1]+x_no[ino][2]; 

    Cht::init();

  }
  
  // titantium thermal conductivity: 24.5 W/m/K

  void run() { 

    const double bc_tol    = getDoubleParam("BC_XTOL", 1.0e-8);
    const double thickness = getDoubleParam("METAL_THICKNESS");
    const double xlo       = getDoubleParam("XLO");
    const double Tlo       = getDoubleParam("TLO");
    const double xhi       = getDoubleParam("XHI");
    const double Thi       = getDoubleParam("THI");
    //const double Qonlt     = 492.31; // W/m^3  
    //const double Qonlt     = 837.0; // W/m^3  
    const double Qonlt     = 899.0; // W/m^3  

    assert( rhs_no != NULL);
    assert( laplace_noono != NULL);
    assert( A != NULL);
    assert( T_no != NULL);
    assert( vol_no != NULL);

    if ( mpi_rank == 0) 
      cout << " ... looks like everyone is initialized " << endl;

    // set a constant heat source rate inside the metal.. 

    double * rhs = new double[nno];

    FOR_INO rhs[ino] = vol_no[ino]*Qonlt*k;

    // assemble the operator .. 
    
    for (int non = 0; non < noono_i[nno]; ++non) 
      A[non] = -k*laplace_noono[non];

    // inject boundary conditions at the two ends .. 

    int * no_flag = new int[nno];
    FOR_INO no_flag[ino] = 0;

    // dont have boundary condition structures so we need to 
    // grab the tris based on coordinates.  assume that the 
    // nodes for the two different boundary conditions are 
    // distinct and uniquely specified (ie no corners.. )

    int8 my_buf[2] = {0,0};

    for (int ist = 0; ist < nst; ++ist) { 

      double x_c[3] = {0.0,0.0,0.0};
      FOR_I3 FOR_J3 x_c[j] += x_no[noost[ist][i]][j]/3.0;

      if ( abs(x_c[0] - xlo) < bc_tol) { 
        
        ++my_buf[0];

        FOR_I3 { 
          const int ino = noost[ist][i];
          assert( (no_flag[ino] ==0) || (no_flag[ino] == 1));
          no_flag[ino] = 1;
        }
      
      } else if ( abs(x_c[0] - xhi) < bc_tol) { 

        ++my_buf[1];

        FOR_I3 { 
          const int ino = noost[ist][i];
          assert( (no_flag[ino] == 0) || (no_flag[ino] == 2));
          no_flag[ino] = 2;
        }

      }
    }

    // report the number of boundary condition triangles found .. 

    int8 buf[2];
    MPI_Reduce(my_buf,buf,2,MPI_INT8,MPI_SUM,0,mpi_comm);

    if ( mpi_rank == 0) 
      cout << " found nst for x0, x1 bcs : " << buf[0] << "    " << buf[1] << endl;

    FOR_INO { 
      
      if ( no_flag[ino] == 1) { 

        const double diag = pow(vol_no[ino],1.0/3.0);
        A[noono_i[ino]]   = diag;

        for (int non = noono_i[ino]+1; non != noono_i[ino+1]; ++non) {
            A[non] = 0.0;
        
            // also modify the matrix coeff in row ino_nbr if it is an internal 
            const int ino_nbr = noono_v[non];
            if (no_flag[ino_nbr] == 0) {
              // look for ino in ino_nbr's row...
              int non_nbr = noono_i[ino_nbr]+1; // no need to consider the diagonal
              while (noono_v[non_nbr] != ino)
                ++non_nbr;
              assert(non_nbr < noono_i[ino_nbr+1]);
              rhs[ino_nbr] -= A[non_nbr]*Tlo;
              A[non_nbr] = 0.0;
            }
          }
          rhs[ino] = diag*Tlo;

      } else if ( no_flag[ino] == 2) { 

        const double diag = pow(vol_no[ino],1.0/3.0);
        A[noono_i[ino]]   = diag;

        for (int non = noono_i[ino]+1; non != noono_i[ino+1]; ++non) {
            A[non] = 0.0;
        
            // also modify the matrix coeff in row ino_nbr if it is an internal 
            const int ino_nbr = noono_v[non];
            if (no_flag[ino_nbr] == 0) {
              // look for ino in ino_nbr's row...
              int non_nbr = noono_i[ino_nbr]+1; // no need to consider the diagonal
              while (noono_v[non_nbr] != ino)
                ++non_nbr;
              assert(non_nbr < noono_i[ino_nbr+1]);
              rhs[ino_nbr] -= A[non_nbr]*Thi;
              A[non_nbr] = 0.0;
            }
          }
          rhs[ino] = diag*Thi;
      
      }
    }

    // set an initial guess for the temperature .. 

    FOR_INO T_no[ino] = 0.5*(Tlo + Thi);

    solveNoCg(T_no,A,rhs,1.0e-8,2000,true);

    MiscUtils::dumpRange(T_no,nno,"final temperature");


    // try to verify the streamwise temperature profile .. 

    for (int ist = 0; ist < nst; ++ist) { 

      double x_c[3] = {0.0,0.0,0.0};
      FOR_I3 FOR_J3 x_c[j] += x_no[noost[ist][i]][j]/3.0;

      if ( abs(x_c[1]) < bc_tol) { 

        double T_avg = 0.0;
        FOR_I3 T_avg += T_no[noost[ist][i]]/3.0;
        cout << " XXXX x T = " << x_c[0] << "    " << T_avg << endl;

      }
    }

    delete[] rhs;
    delete[] no_flag;

  }
};

int main(int argc, char* argv[]) { 


  CTI_Init(argc,argv,"cht_test.in");
  ChtTest solver;
  solver.init();
  solver.run();
  CTI_Finalize();
  return 0;

}
