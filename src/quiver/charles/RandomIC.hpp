#ifndef RANDOMIC_HPP
#define RANDOMIC_HPP

class RandomIc : public IdealGasSolver {
public:

  RandomIc() {}

  void initialHook() { 
    
    if ( checkParam("DEBUG_IC")) { 
      const double rho_ref = getDoubleParam("RHO_REF");
      const double p_ref   = getDoubleParam("P_REF");
      const double gamma   = getDoubleParam("GAMMA", 1.4);
      
      const double eps = 1.0E-12;
      
      for (int icv = 0; icv < ncv; ++icv) { 
        const double sos0      = sqrt(gamma*p_ref/rho_ref);
        const double rho_prime = eps*double(rand())/double(RAND_MAX);
        rho[icv] = rho_ref + rho_prime;
        FOR_I3 { 
          u[icv][i] = eps*2.0*(2.0*double(rand())/double(RAND_MAX) - 1.0);
        }
        rhoE[icv]  = (p_ref + sos0*sos0*rho_prime)/(gamma-1.0);
        rhoE[icv] += 0.5*rho[icv]*DOT_PRODUCT(u[icv],u[icv]);
      }
      
      dumpRange(rho,ncv,"rho_init");
      dumpRange(u,ncv,"u_init");
      dumpRange(rhoE,ncv,"rhoE_init");
    }
  }
};

#endif
