#ifndef CTIINTERPOLATETRANSFORM_HPP
#define CTIINTERPOLATETRANSFORM_HPP

#include "CTI.hpp"
#include "CtiRegister.hpp"
#include "MiscUtils.hpp"

using namespace CTI;
using namespace CtiRegister; 
using namespace MiscUtils;

inline bool myInterpComp(const std::pair<double, double> &a, const std::pair<double, double> &b) {
      return a.first < b.first;
}

class CtiInterpolateTransform : public CtiDataProducer { 
public: 

  double * x;
  double * y;
  int n;
  string func_name;

  CtiInterpolateTransform(const string& _func_name, const string& filename) { 

    addCtiDataProducer(this);

    func_name = _func_name;

    x = NULL; 
    y = NULL;
    n = -1;

    if ( mpi_rank == 0) { 

      ifstream infile(filename.c_str());

      if ( infile.is_open()) { 

        vector<pair<double,double> > tmp_x;
        double xx,yy;
        while ( infile >> xx >> yy)  
          tmp_x.push_back(pair<double,double>(xx,yy));

        infile.close();

        sort(tmp_x.begin(),tmp_x.end(),myInterpComp);

        n = tmp_x.size();
        x = new double[n];
        y = new double[n];

        for (int i = 0; i < n; ++i) {
          x[i] = tmp_x[i].first;
          y[i] = tmp_x[i].second;

        }

      }


    }

    MPI_Bcast(&n,1,MPI_INT,0,mpi_comm);

    if ( n == -1) { 
      CERR( " > unable to open file: " << filename );
    }

    if ( mpi_rank != 0) { 

      x = new double[n];
      y = new double[n];

    } 

    MPI_Bcast(x,n,MPI_DOUBLE,0,mpi_comm);
    MPI_Bcast(y,n,MPI_DOUBLE,0,mpi_comm);

  }


  virtual ~CtiInterpolateTransform() { 

    DELETE(x);
    DELETE(y);

  }


  double interp1d(const double& _x) const { 

    if ( _x <= x[0]) { 
      
      return y[0]; 
    
    } else if ( _x >= x[n-1]) { 
      
      return y[n-1];
    
    } else { 
   
      int i; findInterval(i,_x,x,n);
      assert( (i >=0) && (i < n-1));
      assert( (_x >= x[i]) && (_x <= x[i+1]));

      return y[i] + (_x - x[i])/(x[i+1]-x[i])*(y[i+1] - y[i]);

    }
  }

  CtiRegister::CtiDataError funcEvalCtiData(CtiRegister::CtiData& v,
                                            const string& name,
                                            list<CtiRegister::CtiData>& args,
                                            const bool b_eval_func) {

    
    if ( name == func_name) { 
      
      
      if ( args.size() != 1) 
        return CTI_DATA_ARG_COUNT;
      
      list<CtiRegister::CtiData>::iterator arg = args.begin();
      if ( arg->getType() == DN_DATA) { 
        
        v.new_dn(*arg);
        
        if ( b_eval_func) { 
          
          for (int i = 0, lim = v.size(); i < lim; ++i)  
            v.dn(i) = interp1d(arg->dn(i));
        
        }
        
        return CTI_DATA_OK;

      } else { 

        return CTI_DATA_NOT_VALID;

      }
    
    }

    return CTI_DATA_NOT_FOUND;

  }
};

#endif
