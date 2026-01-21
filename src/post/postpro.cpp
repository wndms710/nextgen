#include "CTI.hpp"
using namespace CTI;

#include "PostproSolver.hpp"

class CharlesPostpro : public PostproSolver {
public:

  int step;
  double time,dt;
  double (*u)[3];
  double *rho;
  double *Z;
  double *T;
  
  CharlesPostpro() {
  
    COUT1("CharlesPostpro()");
    
    registerData(step,"step",READWRITE_DATA);
    registerData(time,"time",READWRITE_DATA);
    registerData(dt,"dt",READWRITE_DATA);
    
    u = NULL; registerCvData(u,"u",READWRITE_DATA);
    rho = NULL; registerCvData(rho,"rho",READWRITE_DATA);
    Z = NULL; registerCvData(Z,"Z",READWRITE_DATA);
    T = NULL; registerCvData(T,"T",READWRITE_DATA);
    
  }

  void initData() {
    assert(u == NULL); u = new double[ncv_g][3];
    assert(rho == NULL); rho = new double[ncv_g];
    assert(Z == NULL); Z = new double[ncv_g];
    assert(T == NULL); T = new double[ncv_g];
  }
 
  ~CharlesPostpro() {
    
    COUT1("~CharlesPostpro()");
    
    DELETE(u);
    DELETE(rho);
    DELETE(Z);
    DELETE(T);
    
  }
 
  virtual void temporalHook() {
    
    COUT1("CharlesPostpro::temporalHook()");
    
    if (mpi_rank == 0) {
      cout << " > step: " << step << endl;
      cout << " > time: " << time << endl;
      cout << " > dt:   " << dt << endl;
    }
    
    MiscUtils::dumpRange(u,ncv,"u");
    MiscUtils::dumpRange(rho,ncv,"rho");
    MiscUtils::dumpRange(Z,ncv,"Z");
    MiscUtils::dumpRange(T,ncv,"T");
    
  }

  virtual CtiRegister::CtiDataError funcEvalCtiData(CtiRegister::CtiData& v,const string& name,list<CtiRegister::CtiData>& args, const bool b_eval_func) {
    
    if (mpi_rank == 0) cout << "CharlesPostpro::funcEvalCtiData() being asked to evaluate: \"" << name << "\"" << endl;
    
    if (name == "mdot_fa") {
      
      if (!args.empty()) 
	return(CtiRegister::CTI_DATA_ARG_COUNT);

      v.new_dn(FA_DATA,nfa);

      updateCvData(u);
      updateCvData(rho);
      
      FOR_IFA {
	const int icv0 = cvofa[ifa][0];
	const int icv1 = cvofa[ifa][1];
	const double rho_avg = 0.5*(rho[icv0]+rho[icv1]);
	const double v_avg   = rho_avg/(rho[icv0]*rho[icv1]); // specific volume average ...
	v.dn(ifa) = 0.5*( (u[icv0][0]+u[icv1][0])*n_fa[ifa][0] + 
			  (u[icv0][1]+u[icv1][1])*n_fa[ifa][1] + 
			  (u[icv0][2]+u[icv1][2])*n_fa[ifa][2] )/v_avg;
      }
      
      return CtiRegister::CTI_DATA_OK;
      
    }
    if (name == "mdotZ_fa") {

      if (!args.empty()) 
	return(CtiRegister::CTI_DATA_ARG_COUNT);

      v.new_dn(FA_DATA,nfa);

      updateCvData(u);
      updateCvData(rho);
      updateCvData(Z);
      
      FOR_IFA {
	const int icv0 = cvofa[ifa][0];
	const int icv1 = cvofa[ifa][1];
	const double rho_avg = 0.5*(rho[icv0]+rho[icv1]);
	const double v_avg   = rho_avg/(rho[icv0]*rho[icv1]); // specific volume average ...
	const double mdot = 0.5*( (u[icv0][0]+u[icv1][0])*n_fa[ifa][0] + 
				  (u[icv0][1]+u[icv1][1])*n_fa[ifa][1] + 
				  (u[icv0][2]+u[icv1][2])*n_fa[ifa][2] )/v_avg;
	v.dn(ifa) = mdot*0.5*(Z[icv0]+Z[icv1]);
	
      }
      
      return CtiRegister::CTI_DATA_OK;
      
    }
    else if (name == "diff_fa") {
      
      if (args.size() != 1)
	return(CtiRegister::CTI_DATA_ARG_COUNT);
      
      list<CtiRegister::CtiData>::iterator arg = args.begin();
      
      if ((arg->getTopology() == CV_DATA)&&(arg->getType() == DN_DATA)) {

        //double * dphi_fa = new double[nfa];
        v.new_dn(FA_DATA,nfa);
        FOR_INTERPROC_IFA {
          const int icv0 = cvofa[ifa][0]; assert((icv0 >= 0)&&(icv0 < ncv));
          //dphi_fa[ifa] = -arg->dn(icv0);
          v.dn(ifa) = -arg->dn(icv0);
        }
        //updateFaDataStart(dphi_fa,SUBTRACT_DATA);
        updateFaDataStart(v.getDNptr(),SUBTRACT_DATA);
        FOR_INTERNAL_IFA {
          const int icv0 = cvofa[ifa][0]; assert((icv0 >= 0)&&(icv0 < ncv));
          const int icv1 = cvofa[ifa][1]; assert((icv1 >= 0)&&(icv1 < ncv));
          //dphi_fa[ifa] = arg->dn(icv1) - arg->dn(icv0);
          v.dn(ifa) = arg->dn(icv1) - arg->dn(icv0);
        }
        //updateFaDataFinish(dphi_fa,SUBTRACT_DATA);
        updateFaDataFinish(v.getDNptr(),SUBTRACT_DATA);

        //v.new_dn(FA_DATA,nfa);
        //FOR_IFA {
        //  v.dn(ifa) = dphi_fa[ifa];
        //}
        //delete[] dphi_fa;
	
        return CtiRegister::CTI_DATA_OK;
	
      }
      else {

	return CtiRegister::CTI_DATA_NOT_VALID;
      
      }
      
    }
    else if (name == "avg_fa") {
      
      if (args.size() != 1)
	return(CtiRegister::CTI_DATA_ARG_COUNT);
      
      list<CtiRegister::CtiData>::iterator arg = args.begin();
    
      if (arg->getTopology() == CV_DATA) {
        if (arg->getType() == DN_DATA) {

          v.new_dn(FA_DATA,nfa);
          FOR_INTERPROC_IFA {
            const int icv0 = cvofa[ifa][0]; assert((icv0 >= 0)&&(icv0 < ncv));
            v.dn(ifa) = arg->dn(icv0);
          }
          updateFaDataStart(v.getDNptr(),ADD_DATA);
          FOR_INTERNAL_IFA {
            const int icv0 = cvofa[ifa][0]; assert((icv0 >= 0)&&(icv0 < ncv));
            const int icv1 = cvofa[ifa][1]; assert((icv1 >= 0)&&(icv1 < ncv));
            v.dn(ifa) = arg->dn(icv1) + arg->dn(icv0);
          }
          updateFaDataFinish(v.getDNptr(),ADD_DATA);
          FOR_IFA v.dn(ifa) *= 0.5; 
          
          return CtiRegister::CTI_DATA_OK;
          
        }
        else if (arg->getType() == DN3_DATA) {

          v.new_dn3(FA_DATA,nfa);
          FOR_INTERPROC_IFA {
            const int icv0 = cvofa[ifa][0]; assert((icv0 >= 0)&&(icv0 < ncv));
            FOR_I3 v.dn3(ifa,i) = arg->dn3(icv0,i);
          }
          updateFaDataStart(v.getDN3ptr(),ADD_DATA);
          FOR_INTERNAL_IFA {
            const int icv0 = cvofa[ifa][0]; assert((icv0 >= 0)&&(icv0 < ncv));
            const int icv1 = cvofa[ifa][1]; assert((icv1 >= 0)&&(icv1 < ncv));
            FOR_I3 v.dn3(ifa,i) = arg->dn3(icv1,i) + arg->dn3(icv0,i);
          }
          updateFaDataFinish(v.getDN3ptr(),ADD_DATA);
          FOR_IFA FOR_I3 v.dn3(ifa,i) *= 0.5; 

          return CtiRegister::CTI_DATA_OK;
          
        }
        else {

          return CtiRegister::CTI_DATA_NOT_VALID;

        }
      }
      else {

	return CtiRegister::CTI_DATA_NOT_VALID;

      }
	
    }
    
    return StaticSolver::funcEvalCtiData(v,name,args,b_eval_func);
    
  }

};


int main(int argc, char * argv[]) {

  try {
    
    CTI_Init(argc,argv,"postpro.in");
    
    {
      
      //PostproSolver solver;
      CharlesPostpro solver;
      solver.init();
      solver.run();
      
    }
    
    CTI_Finalize();
    
  }
  catch (int e) {
    if (e == 0)
      CTI_Finalize();
    else
      CTI_Abort();
  }
  catch(...) {
    CTI_Abort();
  }

  return(0);

}

