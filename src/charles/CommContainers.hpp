#ifndef COMMCONTAINERS_HPP
#define COMMCONTAINERS_HPP

// the following class is a container to allow 
// for the packing of the state variables for 
// a batched message to be sent to update the ghosts..

//===================================
// non-reacting ideal gas solver
//==================================

class IgComm { 
public: 
  double *rho;
  double (*u)[3];
  double *rhoE;

  IgComm(double* _rho, double (*_u)[3], double* _rhoE) : rho(_rho), u(_u), rhoE(_rhoE) {} 
  static int data_size() { return 5;}
};

inline void pack_class(double* buf, const IgComm* s, const int icv) { 
  buf[0] = s->rho[icv];
  buf[1] = s->u[icv][0];
  buf[2] = s->u[icv][1];
  buf[3] = s->u[icv][2];
  buf[4] = s->rhoE[icv];
}

inline void pack_class(double* buf, const IgComm* s, const int icv, const double * R) { 
  buf[0] = s->rho[icv];
  MiscUtils::matVecMult(&buf[1],R,s->u[icv]);
  buf[4] = s->rhoE[icv];
}

inline void unpack_class(double* buf, const IgComm* s, const int icv) { 
  s->rho[icv]  = buf[0];
  s->u[icv][0] = buf[1];
  s->u[icv][1] = buf[2];
  s->u[icv][2] = buf[3];
  s->rhoE[icv] = buf[4];
}

//===================================
// premixed reacting solver
//===================================

class PremixedComm { 
public: 
  double * rho;
  double (*u)[3];
  double * rhoE;
  double * Z;
  double * C;

  PremixedComm(double* _rho, double (*_u)[3], double* _rhoE, double* _Z, double* _C) : 
    rho(_rho), u(_u), rhoE(_rhoE), Z(_Z), C(_C) {} 
  static int data_size() { return 7;}
};

inline void pack_class(double* buf, const PremixedComm* s, const int icv) { 
  buf[0] = s->rho[icv];
  buf[1] = s->u[icv][0];
  buf[2] = s->u[icv][1];
  buf[3] = s->u[icv][2];
  buf[4] = s->rhoE[icv];
  buf[5] = s->Z[icv];
  buf[6] = s->C[icv];
}

inline void pack_class(double* buf, const PremixedComm* s, const int icv, const double *R ) { 
  buf[0] = s->rho[icv];
  MiscUtils::matVecMult(&buf[1],R,s->u[icv]);
  buf[4] = s->rhoE[icv];
  buf[5] = s->Z[icv];
  buf[6] = s->C[icv];
}


inline void unpack_class(double* buf, const PremixedComm* s, const int icv) { 
  s->rho[icv]  = buf[0];
  s->u[icv][0] = buf[1];
  s->u[icv][1] = buf[2];
  s->u[icv][2] = buf[3];
  s->rhoE[icv] = buf[4];
  s->Z[icv]    = buf[5];
  s->C[icv]    = buf[6];
}

//===============================================
// nonpremixed communicator 
// this is for use with an algebraic closure of 
// the mixture fraction variance s.t. the mixture 
// fraction is updated separately
//===============================================

class NonpremixedComm { 
public: 

  double * rho;
  double (*u)[3];
  double * rhoE;
  double * C;

  NonpremixedComm(double* _rho, double (*_u)[3], double* _rhoE, double* _C) : 
    rho(_rho), u(_u), rhoE(_rhoE), C(_C) {} 
  static int data_size() { return 6;}

};

inline void pack_class(double* buf, const NonpremixedComm* s, const int icv) { 
  buf[0] = s->rho[icv];
  buf[1] = s->u[icv][0];
  buf[2] = s->u[icv][1];
  buf[3] = s->u[icv][2];
  buf[4] = s->rhoE[icv];
  buf[5] = s->C[icv];
}

inline void pack_class(double* buf, const NonpremixedComm* s, const int icv, const double *R ) { 
  buf[0] = s->rho[icv];
  MiscUtils::matVecMult(&buf[1],R,s->u[icv]);
  buf[4] = s->rhoE[icv];
  buf[5] = s->C[icv];
}


inline void unpack_class(double* buf, const NonpremixedComm* s, const int icv) { 
  s->rho[icv]  = buf[0];
  s->u[icv][0] = buf[1];
  s->u[icv][1] = buf[2];
  s->u[icv][2] = buf[3];
  s->rhoE[icv] = buf[4];
  s->C[icv]    = buf[5];
}


//====================================
// sgs communicator (compressible solver)
//====================================

class SgsComm { 
public: 
  double * mu_sgs;
  double * a_sgs;
  double * loc_sgs;

  SgsComm(double *_mu_sgs, double *_a_sgs, double* _loc_sgs) : 
    mu_sgs(_mu_sgs), a_sgs(_a_sgs), loc_sgs(_loc_sgs) {}
  static int data_size() { return 3;}
};

inline void pack_class(double* buf, const SgsComm* s, const int icv) { 
  buf[0] = s->mu_sgs[icv];
  buf[1] = s->a_sgs[icv];
  buf[2] = s->loc_sgs[icv];
}

inline void pack_class(double* buf, const SgsComm* s, const int icv, const double* R) { 
  buf[0] = s->mu_sgs[icv];
  buf[1] = s->a_sgs[icv];
  buf[2] = s->loc_sgs[icv];
}

inline void unpack_class(double* buf, const SgsComm* s, const int icv) { 
  s->mu_sgs[icv] = buf[0];
  s->a_sgs[icv]  = buf[1];
  s->loc_sgs[icv] = buf[2];
}

#endif
