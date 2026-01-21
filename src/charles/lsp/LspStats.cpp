#include "LspStats.hpp"
#include "CtiLiquid.hpp"

LspStats::LspStats(Param * param) {

  name = "";
  param_line = "# param line: ";
  param_line += getParamLine(param);
  interval = -1;

  nbins = nbins_pdf = -1;
  dmin = dmax = -1;
  dt_sum = 0.0;

  np_sum = NULL;
  nm_sum = NULL;
  mp_sum = NULL;
  mm_sum = NULL;
  mup_sum = NULL;
  mum_sum = NULL;
  d2p_sum = NULL;
  d2m_sum = NULL;
  d3p_sum = NULL;
  d3m_sum = NULL;

  got_pdf = false;
  diam = NULL;
  pdf_mass_p = NULL;
  pdf_mass_m = NULL;
  pdf_num_p = NULL;
  pdf_num_m = NULL;
}

LspStats::~LspStats() {
  DELETE(np_sum);
  DELETE(nm_sum);
  DELETE(mp_sum);
  DELETE(mm_sum);
  DELETE(mup_sum);
  DELETE(mum_sum);
  DELETE(d2p_sum);
  DELETE(d2m_sum);
  DELETE(d3p_sum);
  DELETE(d3m_sum);
  if (pdf_mass_p != NULL) DELETE(pdf_mass_p);
  if (pdf_mass_m != NULL) DELETE(pdf_mass_m);
  if (pdf_num_p != NULL) DELETE(pdf_num_p);
  if (pdf_num_m != NULL) DELETE(pdf_num_m);
  if (diam != NULL) DELETE(diam);
}

LspStatsDisk::LspStatsDisk(Param * param) : LspStats(param) {

  FOR_I3 xc[i] = 0.0;
  FOR_I3 normal[i] = 0.0;
  FOR_I3 t1[i] = 0.0;
  FOR_I3 t2[i] = 0.0;

  r0 = r1 = theta0 = theta1 = -1;
  nr = ntheta = -1;
}

LspStatsRect::LspStatsRect(Param * param) : LspStats(param) {

  FOR_I3 xc[i] = 0.0;
  FOR_I3 normal[i] = 0.0;
  FOR_I3 t1[i] = 0.0;
  FOR_I3 t2[i] = 0.0;

  s1 = s2 = -1;
  nbins_s1 = nbins_s2 = -1;
}

string LspStats::getParamLine(Param * param) {
  string param_line = "LSP.STATS ";
  int iarg = 0;
  while (iarg < param->size()) {
    param_line += param->getString(iarg++);
    param_line += " ";
  }
  param_line += "\n";
  return param_line;
}

void LspStatsDisk::init(Param * param) {

  bool got_name = false;
  bool got_interval = false;
  bool got_geom = false;

  int iarg = 0;
  while (iarg < param->size()) {
    string token = param->getString(iarg++);
    if (token == "NAME") {
      if (got_name) CERR("in LSP.STATS wrong use of NAME!");
      got_name = true;
      name = param->getString(iarg++);
    }
    else if (token == "INTERVAL") {
      if (got_interval) CERR("in LSP.STATS wrong use of INTERVAL!");
      got_interval = true;
      interval = param->getInt(iarg++);
      if (interval < 1) CERR("in LSP.STATS wrong value of INTERVAL!");
    }
    else if (token == "PDF") {
      if (got_pdf) CERR("in LSP.STATS wrong use of PDF!");
      this->got_pdf = true;
      token = param->getString(iarg++);
      if (token != "NBINS") CERR("in LSP.STATS wrong position of NBINS after PDF!");
      nbins_pdf = param->getInt(iarg++);

      token = param->getString(iarg++);
      if (token != "DRANGE") CERR("in LSP.STATS wrong position of DRANGE after PDF!");
      dmin = param->getDouble(iarg++);
      dmax = param->getDouble(iarg++);
      if (dmin >= dmax) CERR("in LSP.STATS wrong values for dmin and dmax");
    }  
    else if (token == "GEOM") {
      if (got_geom) CERR("in LSP.STATS wrong use of GEOM!");
      got_geom = true;
      token = param->getString(iarg++);
      assert(token=="DISK");
      if (param->getString(iarg++) != "X") CERR("in LSP.STATS expect X in DISK!");
      xc[0] = param->getDouble(iarg++);
      xc[1] = param->getDouble(iarg++);
      xc[2] = param->getDouble(iarg++);
      if (param->getString(iarg++) != "N") CERR("in LSP.STATS expect N in DISK!");
      normal[0] = param->getDouble(iarg++);
      normal[1] = param->getDouble(iarg++);
      normal[2] = param->getDouble(iarg++);
      // make this a unit normal...
      NORMALIZE(normal);
      double tmp[3];
      if (fabs(normal[0]) <= min(fabs(normal[1]),fabs(normal[2]))) {
        tmp[0] = 1.0;
        tmp[1] = 0.0;
        tmp[2] = 0.0;
      }
      else if (fabs(normal[1]) <= min(fabs(normal[0]),fabs(normal[2]))) {
        tmp[0] = 0.0;
        tmp[1] = 1.0;
        tmp[2] = 0.0;
      }
      else {
        tmp[0] = 0.0;
        tmp[1] = 0.0;
        tmp[2] = 1.0;
      }
      t1[0] = CROSS_PRODUCT_0(tmp,normal);
      t1[1] = CROSS_PRODUCT_1(tmp,normal);
      t1[2] = CROSS_PRODUCT_2(tmp,normal);
      NORMALIZE(t1);
      t2[0] = CROSS_PRODUCT_0(normal,t1);
      t2[1] = CROSS_PRODUCT_1(normal,t1);
      t2[2] = CROSS_PRODUCT_2(normal,t1);
      assert(DOT_PRODUCT(t1,t2)<1e-10);
      assert(DOT_PRODUCT(t1,normal)<1e-10);
      assert(DOT_PRODUCT(t2,normal)<1e-10);
      if (param->getString(iarg++) != "R0") CERR("in LSP.STATS expect R0 in DISK!");
      r0 = param->getDouble(iarg++);
      if (param->getString(iarg++) != "R1") CERR("in LSP.STATS expect R1 in DISK!");
      r1 = param->getDouble(iarg++);
      if (param->getString(iarg++) != "NR") CERR("in LSP.STATS expect NR in DISK!");
      nr = param->getInt(iarg++);
      if (param->getString(iarg++) != "THETA0") CERR("in LSP.STATS expect THETA0 in DISK!");
      theta0 = M_PI/180.0*param->getDouble(iarg++);
      if (param->getString(iarg++) != "THETA1") CERR("in LSP.STATS expect THETA1 in DISK!");
      theta1 = M_PI/180.0*param->getDouble(iarg++);
      if (param->getString(iarg++) != "NTHETA") CERR("in LSP.STATS expect NTHETA in DISK!");
      ntheta = param->getInt(iarg++);

      nbins = nr*ntheta;
    }
    else {
      CERR("in LSP.STATS unrecognized token: "+token);
    }
  }

  //int ierr = 0;
  if (!got_name)
    CERR("Error: LSP.STATS missing NAME");
  if (!got_interval)
    CERR("Error: LSP.STATS missing INTERVAL");
  if (!got_geom) 
    CERR("Error: LSP.STATS missing GEOM");

  assert(nbins>0);
  np_sum  =  new double[nbins];
  nm_sum  =  new double[nbins];
  mp_sum  =  new double[nbins];
  mm_sum  =  new double[nbins];
  mup_sum  =  new double[nbins];
  mum_sum  =  new double[nbins];
  d2p_sum =  new double[nbins];
  d2m_sum =  new double[nbins];
  d3p_sum =  new double[nbins];
  d3m_sum =  new double[nbins];

  if (this->got_pdf) {
    assert(nbins_pdf>0);
    pdf_mass_p = new double[nbins_pdf];
    pdf_mass_m = new double[nbins_pdf];
    pdf_num_p = new double[nbins_pdf];
    pdf_num_m = new double[nbins_pdf];
    diam = new double[nbins_pdf];
    for (int ibin = 0; ibin < nbins_pdf; ++ibin) {
      diam[ibin]  =  getD_pdf(ibin);
      pdf_mass_p[ibin]  =  0.0;
      pdf_mass_m[ibin]  =  0.0;
      pdf_num_p[ibin]  =  0.0;
      pdf_num_m[ibin]  =  0.0;
    }
  }

  // and zero everything. Here we could read from restart file some day...

  for (int bin = 0; bin < nbins; ++bin) {
    np_sum[bin]  =  0.0;
    nm_sum[bin]  =  0.0;
    mp_sum[bin]  =  0.0;
    mm_sum[bin]  =  0.0;
    mup_sum[bin]  =  0.0;
    mum_sum[bin]  =  0.0;
    d2p_sum[bin] =  0.0;
    d2m_sum[bin] =  0.0;
    d3p_sum[bin] =  0.0;
    d3m_sum[bin] =  0.0;
  }

  // also initialize the computational geometry helper stuff...

  // report 
  IF_RANK0 {
    cout << "LSP.STATS: " << name << " x: " << COUT_VEC(xc) << " n: " << COUT_VEC(normal) << " r: " << r0 << "-" << r1 << " theta: " << theta0 << "-" << theta1 << " nr: " << nr << " ntheta: " << ntheta << " nbins: " << nbins << endl;
    if (this->got_pdf) cout << "LSP.STATS got pdf, nbins: " << nbins_pdf << " dmin: " << dmin << " dmax: " << dmax << endl;
  }
}

void LspStatsRect::init(Param * param) {

  bool got_name = false;
  bool got_interval = false;
  bool got_geom = false;

  int iarg = 0;
  while (iarg < param->size()) {
    string token = param->getString(iarg++);
    if (token == "NAME") {
      if (got_name) CERR("in LSP.STATS wrong use of NAME!");
      got_name = true;
      name = param->getString(iarg++);
    }
    else if (token == "INTERVAL") {
      if (got_interval) CERR("in LSP.STATS wrong use of INTERVAL!");
      got_interval = true;
      interval = param->getInt(iarg++);
      if (interval < 1) CERR("in LSP.STATS wrong value of INTERVAL!");
    }
    else if (token == "PDF") {
      if (got_pdf) CERR("in LSP.STATS wrong use of PDF!");
      this->got_pdf = true;
      token = param->getString(iarg++);
      if (token != "NBINS") CERR("in LSP.STATS wrong position of NBINS after PDF!");
      nbins_pdf = param->getInt(iarg++);

      token = param->getString(iarg++);
      if (token != "DRANGE") CERR("in LSP.STATS wrong position of DRANGE after PDF!");
      dmin = param->getDouble(iarg++);
      dmax = param->getDouble(iarg++);
      if (dmin >= dmax) CERR("in LSP.STATS wrong values for dmin and dmax");
    }  
    else if (token == "GEOM") {
      if (got_geom) CERR("in LSP.STATS wrong use of GEOM!");
      got_geom = true;
      token = param->getString(iarg++);
      assert(token=="RECTANGLE");
      if (param->getString(iarg++) != "X") CERR("in LSP.STATS expect X in RECTANGLE!");
      xc[0] = param->getDouble(iarg++);
      xc[1] = param->getDouble(iarg++);
      xc[2] = param->getDouble(iarg++);
      if (param->getString(iarg++) != "N") CERR("in LSP.STATS expect N in RECTANGLE!");
      normal[0] = param->getDouble(iarg++);
      normal[1] = param->getDouble(iarg++);
      normal[2] = param->getDouble(iarg++);
      // make this a unit normal...
      NORMALIZE(normal);
      if (param->getString(iarg++) != "T") CERR("in LSP.STATS expect T in RECTANGLE!");
      t1[0] = param->getDouble(iarg++);
      t1[1] = param->getDouble(iarg++);
      t1[2] = param->getDouble(iarg++);
      // make this a unit normal...
      NORMALIZE(t1);
      // t1 and normal should be orthogonal
      if (DOT_PRODUCT(t1,normal)>1e-10) CERR("in LSP.STATS expect T and N be orthogonal in RECTANGLE!");
      const double tmp[3] = CROSS_PRODUCT(normal,t1);
      FOR_J3 t2[j] = tmp[j];
      assert(DOT_PRODUCT(t1,t2)<1e-10);
      assert(DOT_PRODUCT(t1,normal)<1e-10);
      assert(DOT_PRODUCT(t2,normal)<1e-10);
      if (param->getString(iarg++) != "SIDE1") CERR("in LSP.STATS expect SIDE1 in RECTANGLE!");
      s1 = param->getDouble(iarg++);
      if (param->getString(iarg++) != "SIDE2") CERR("in LSP.STATS expect SIDE2 in RECTANGLE!");
      s2 = param->getDouble(iarg++);
      if (param->getString(iarg++) != "NBINS1") CERR("in LSP.STATS expect NBINS1 in RECTANGLE!");
      nbins_s1 = param->getInt(iarg++);
      if (param->getString(iarg++) != "NBINS2") CERR("in LSP.STATS expect NBINS2 in RECTANGLE!");
      nbins_s2 = param->getInt(iarg++);

      nbins = nbins_s1*nbins_s2;
    }
    else {
      CERR("in LSP.STATS unrecognized token: "+token);
    }
  }

  //int ierr = 0;
  if (!got_name)
    CERR("Error: LSP.STATS missing NAME");
  if (!got_interval)
    CERR("Error: LSP.STATS missing INTERVAL");
  if (!got_geom) 
    CERR("Error: LSP.STATS missing GEOM");

  assert(nbins>0);
  np_sum  =  new double[nbins];
  nm_sum  =  new double[nbins];
  mp_sum  =  new double[nbins];
  mm_sum  =  new double[nbins];
  mup_sum  =  new double[nbins];
  mum_sum  =  new double[nbins];
  d2p_sum =  new double[nbins];
  d2m_sum =  new double[nbins];
  d3p_sum =  new double[nbins];
  d3m_sum =  new double[nbins];

  if (this->got_pdf) {
    assert(nbins_pdf>0);
    pdf_mass_p = new double[nbins_pdf];
    pdf_mass_m = new double[nbins_pdf];
    pdf_num_p = new double[nbins_pdf];
    pdf_num_m = new double[nbins_pdf];
    diam = new double[nbins_pdf];
    for (int ibin = 0; ibin < nbins_pdf; ++ibin) {
      diam[ibin]  =  getD_pdf(ibin);
      pdf_mass_p[ibin]  =  0.0;
      pdf_mass_m[ibin]  =  0.0;
      pdf_num_p[ibin]  =  0.0;
      pdf_num_m[ibin]  =  0.0;
    }
  }

  // and zero everything. Here we could read from restart file some day...

  for (int bin = 0; bin < nbins; ++bin) {
    np_sum[bin]  =  0.0;
    nm_sum[bin]  =  0.0;
    mp_sum[bin]  =  0.0;
    mm_sum[bin]  =  0.0;
    mup_sum[bin]  =  0.0;
    mum_sum[bin]  =  0.0;
    d2p_sum[bin] =  0.0;
    d2m_sum[bin] =  0.0;
    d3p_sum[bin] =  0.0;
    d3m_sum[bin] =  0.0;
  }

  // also initialize the computational geometry helper stuff...

  // report 
  IF_RANK0 {
    cout << "LSP.STATS: " << name << " x: " << COUT_VEC(xc) << " n: " << COUT_VEC(normal) << " t1: " << COUT_VEC(t1) << " t2: " << COUT_VEC(t2) << " side1: " << s1 << " side2: " << s2 << " nbins_s1: " << nbins_s1 << " nbins_s2: " << nbins_s2 << endl;
    if (this->got_pdf) cout << "LSP.STATS got pdf, nbins: " << nbins_pdf << " dmin: " << dmin << " dmax: " << dmax << endl;
  }
}

void LspStats::update(const LspState * lsp,const int np,CtiMaterial* fuel,const double dt,const bool verbose) {

  dt_sum += dt;

  for (int ip = 0; ip < np; ++ip) {

    // check if this particle crossed the plane...

    double dx0[3] = DIFF(lsp[ip].xp0,xc);
    double d0 = DOT_PRODUCT(dx0,normal);

    double dx1[3] = DIFF(lsp[ip].xp,xc);
    double d1 = DOT_PRODUCT(dx1,normal);

    if ((d0 < 0)&&(d1 >= 0)) {
      // the particle crossed from the minus side to the plus side of the plane...
      // the crossing point is...
      double dxc[3];
      FOR_I3 dxc[i] = (-d0*dx1[i] + d1*dx0[i])/(d1-d0);
      int bin = getBin(dxc); // could return -1 if not in a bin
      if (bin >= 0) {
        np_sum[bin] += lsp[ip].npar;
        mp_sum[bin] += lsp[ip].mp;
        const double up_norm = DOT_PRODUCT(lsp[ip].up,normal);
        //assert(up_norm > 0);
        mup_sum[bin] += lsp[ip].mp*up_norm;
        const double rhop  = fuel->calcRho(lsp[ip].Tp);
        const double dp3   = lsp[ip].mp/(M_PI/6.0*rhop*lsp[ip].npar);
        d2p_sum[bin] += lsp[ip].npar * pow(dp3,2.0/3.0);
        d3p_sum[bin] += lsp[ip].npar * dp3;
        // add to pdf only if it passes within the geom
        if (got_pdf) {
          int bin = getBin_pdf(lsp[ip].dp);
          if (bin >= 0) {
            pdf_mass_p[bin] += lsp[ip].mp;
            pdf_num_p[bin] += 1;
          }
        }
      }
    }
    else if ((d0 >= 0)&&(d1 < 0)) {
      // the particle crossed from the plus side to the minus side of the plane...
      // the crossing point is...
      double dxc[3];
      FOR_I3 dxc[i] = (d0*dx1[i] - d1*dx0[i])/(d0-d1);
      int bin = getBin(dxc); // could return -1 if not in a bin
      if (bin >= 0) {
        nm_sum[bin] += lsp[ip].npar;
        mm_sum[bin] += lsp[ip].mp;
        const double up_norm = DOT_PRODUCT(lsp[ip].up,normal);
        //assert(up_norm < 0);
        mum_sum[bin] += lsp[ip].mp*(-up_norm);
        const double rhop  = fuel->calcRho(lsp[ip].Tp);
        const double dp3   = lsp[ip].mp/(M_PI/6.0*rhop*lsp[ip].npar);
        d2m_sum[bin] += lsp[ip].npar * pow(dp3,2.0/3.0);
        d3m_sum[bin] += lsp[ip].npar * dp3;
        // add to pdf only if it passes within the geom
        if (got_pdf) {
          int bin = getBin_pdf(lsp[ip].dp);
          if (bin >= 0) {
            pdf_mass_m[bin] += lsp[ip].mp;
            pdf_num_m[bin] += 1;
          }
        }
      }
    }
  }
}

void LspStats::reset(const bool verbose) {

  if (verbose and (mpi_rank==0)) cout << " > reset all stats to zero!" << endl;

  dt_sum = 0.0;
  for (int bin = 0; bin < nbins; ++bin) {
    np_sum[bin]  =  0.0;
    nm_sum[bin]  =  0.0;
    mp_sum[bin]  =  0.0;
    mm_sum[bin]  =  0.0;
    mup_sum[bin]  =  0.0;
    mum_sum[bin]  =  0.0;
    d2p_sum[bin] =  0.0;
    d2m_sum[bin] =  0.0;
    d3p_sum[bin] =  0.0;
    d3m_sum[bin] =  0.0;
  }
  
  if (got_pdf) {
    for (int bin = 0; bin < nbins_pdf; ++bin) {
      pdf_mass_p[bin] = 0.0;
      pdf_mass_m[bin] = 0.0;
      pdf_num_p[bin] = 0.0;
      pdf_num_m[bin] = 0.0;
    }
  }

}

int LspStatsDisk::getBin(const double dxc[]) {

  double r = sqrt(DOT_PRODUCT(dxc,dxc));
  int r_bin = int((r-r0)/(r1-r0)*double(nr));

  if ((r_bin < 0)||(r_bin >= nr)) {
    return(-1);
  }

  int theta_bin;
  if (r == 0.0)
    theta_bin = 0;
  else {
    double dx = DOT_PRODUCT(dxc,t1);
    double dy = DOT_PRODUCT(dxc,t2);
    double theta = atan2(dy,dx);
    theta_bin = int((theta-theta0)/(theta1-theta0)*double(ntheta));
    if ((theta_bin < 0)||(theta_bin >= ntheta)) {
      return(-1);
    }
  }

  return(r_bin*ntheta + theta_bin);

}

int LspStatsRect::getBin(const double dxc[]) {

  double x1 = s1/2. + DOT_PRODUCT(dxc,t1);
  double x2 = s2/2. + DOT_PRODUCT(dxc,t2);

  int bin_x1 = int(x1/s1*double(nbins_s1));
  int bin_x2 = int(x2/s2*double(nbins_s2));
  
  if ((bin_x1 < 0) or (bin_x1 >= nbins_s1))
    return -1;
  if ((bin_x2 < 0) or (bin_x2 >= nbins_s2))
    return -1;

  return (bin_x2*nbins_s1 + bin_x1);

}

int LspStats::getBin_pdf(const double d) {

  const int d_bin = (int)floor( (d-dmin)/(dmax-dmin)*double(nbins_pdf) );

  if ((d_bin < 0)||(d_bin >= nbins_pdf)) {
    return(-1);
  }

  return(d_bin);

}

double LspStats::getD_pdf(const int ibin) {
  assert((ibin>=0)and(ibin<nbins_pdf));
  const double dD  =  (dmax-dmin)/double(nbins_pdf);
  const double d   =  dmin + (double(ibin) + 0.5)*dD;
  return(d);
}

void LspStats::writePdf(const double time,const int step,const bool verbose) {

  assert(pdf_mass_p);
  assert(pdf_mass_m);
  assert(pdf_num_p);
  assert(pdf_num_m);
  assert(nbins_pdf>0);

  double * pdf_mass_p_global = NULL;
  double * pdf_mass_m_global = NULL;
  double * pdf_num_p_global = NULL;
  double * pdf_num_m_global = NULL;

  if (mpi_rank == 0) {
    pdf_mass_p_global = new double[nbins_pdf];
    pdf_mass_m_global = new double[nbins_pdf];
    pdf_num_p_global = new double[nbins_pdf];
    pdf_num_m_global = new double[nbins_pdf];
  }

  MPI_Reduce(pdf_mass_p,pdf_mass_p_global,nbins_pdf,MPI_DOUBLE,MPI_SUM,0,mpi_comm);
  MPI_Reduce(pdf_mass_m,pdf_mass_m_global,nbins_pdf,MPI_DOUBLE,MPI_SUM,0,mpi_comm);
  MPI_Reduce(pdf_num_p,pdf_num_p_global,nbins_pdf,MPI_DOUBLE,MPI_SUM,0,mpi_comm);
  MPI_Reduce(pdf_num_m,pdf_num_m_global,nbins_pdf,MPI_DOUBLE,MPI_SUM,0,mpi_comm);

  if (mpi_rank == 0) {

    char filename[128];
    buildIndexedFilename(filename,(name+".PDF"+".lspstats").c_str(),step,"dat");
    mkdir_for_file(filename); // allow the user to have specified a subdirectory...
    COUT1(" > writeLSPStats: " << filename);
    FILE * fp = fopen(filename,"w");

    const double dD = (dmax-dmin)/double(nbins_pdf);
    double sum_mass_p = 0.0;
    double sum_mass_m = 0.0;
    double sum_num_p = 0.0;
    double sum_num_m = 0.0;

    for (int ibin = 0; ibin < nbins_pdf; ++ibin) {
      sum_mass_p += pdf_mass_p_global[ibin];
      sum_mass_m += pdf_mass_m_global[ibin];
      sum_num_p += pdf_num_p_global[ibin];
      sum_num_m += pdf_num_m_global[ibin];
    }

    if (sum_mass_p > 0) {
      for (int ibin = 0; ibin < nbins_pdf; ++ibin) 
        pdf_mass_p_global[ibin] /= (dD*sum_mass_p);

      // check the integral of pdf
      double sum_pdf = 0.0;
      for (int ibin = 0; ibin < nbins_pdf; ++ibin) 
        sum_pdf += pdf_mass_p_global[ibin]*dD;
      cout << " > integral pdf_mass_p - 1.0 (should be 0): " << sum_pdf-1.0 << endl;
      assert((sum_pdf-1.0)<1e-12);
    }

    if (sum_mass_m > 0) {
      for (int ibin = 0; ibin < nbins_pdf; ++ibin) 
        pdf_mass_m_global[ibin] /= (dD*sum_mass_m);

      // check the integral of pdf
      double sum_pdf = 0.0;
      for (int ibin = 0; ibin < nbins_pdf; ++ibin) 
        sum_pdf += pdf_mass_m_global[ibin]*dD;
      cout << " > integral pdf_mass_m - 1.0 (should be 0): " << sum_pdf-1.0 << endl;
      assert((sum_pdf-1.0)<1e-12);
  
    }

    if (sum_num_p > 0) {
      for (int ibin = 0; ibin < nbins_pdf; ++ibin) 
        pdf_num_p_global[ibin] /= (dD*sum_num_p);

      // check the integral of pdf
      double sum_pdf = 0.0;
      for (int ibin = 0; ibin < nbins_pdf; ++ibin) 
        sum_pdf += pdf_num_p_global[ibin]*dD;
      cout << " > integral pdf_num_p - 1.0 (should be 0): " << sum_pdf-1.0 << endl;
      assert((sum_pdf-1.0)<1e-12);
    }

    if (sum_num_m > 0) {
      for (int ibin = 0; ibin < nbins_pdf; ++ibin) 
        pdf_num_m_global[ibin] /= (dD*sum_num_m);

      // check the integral of pdf
      double sum_pdf = 0.0;
      for (int ibin = 0; ibin < nbins_pdf; ++ibin) 
        sum_pdf += pdf_num_m_global[ibin]*dD;
      cout << " > integral pdf_num_m - 1.0 (should be 0): " << sum_pdf-1.0 << endl;
      assert((sum_pdf-1.0)<1e-12);
    }

    fprintf(fp,"%s",param_line.c_str());
    fprintf(fp,"# simulation time: %8.16le , stats average time:  %8.16le\n",time,dt_sum);
    fprintf(fp,"# sum number plus: %18.16le , sum mass plus: %18.16le , sum number minus: %18.16le , sum mass minus: %18.16le\n",sum_num_p,sum_mass_p,sum_num_m,sum_mass_m);
    fprintf(fp,"# 0: diameter 1: number pdf plus dir 2: mass pdf plus dir 3: number pdf minus dir 4: mass pdf minus dir\n");
    fprintf(fp,"#             5: number cdf plus dir 6: mass cdf plus dir 7: number cdf minus dir 8: mass cdf minus dir\n");
    double cdf_num_p = 0.0;
    double cdf_num_m = 0.0;
    double cdf_mass_p = 0.0;
    double cdf_mass_m = 0.0;
    for (int ibin = 0; ibin < nbins_pdf; ++ibin) {
      cdf_num_p += pdf_num_p_global[ibin]*dD;
      cdf_num_m += pdf_num_m_global[ibin]*dD;
      cdf_mass_p += pdf_mass_p_global[ibin]*dD;
      cdf_mass_m += pdf_mass_m_global[ibin]*dD;
      fprintf(fp,"%18.15le %18.15le %18.15le %18.15le %18.15le %18.15le %18.15le %18.15le %18.15le\n",diam[ibin], pdf_num_p_global[ibin], pdf_mass_p_global[ibin], pdf_num_m_global[ibin], pdf_mass_m_global[ibin], cdf_num_p, cdf_mass_p, cdf_num_m, cdf_mass_m);
    }
    fclose(fp);

    delete[] pdf_mass_p_global;
    delete[] pdf_mass_m_global;
    delete[] pdf_num_p_global;
    delete[] pdf_num_m_global;

  }

}

void LspStatsDisk::write(const double time,const int step,const bool verbose) {

  // reduce to rank0 and plot...

  double * np_sum_global = NULL;
  double * nm_sum_global = NULL;
  double * mp_sum_global = NULL;
  double * mm_sum_global = NULL;
  double * mup_sum_global = NULL;
  double * mum_sum_global = NULL;
  double * d2p_sum_global = NULL;
  double * d3p_sum_global = NULL;
  double * d2m_sum_global = NULL;
  double * d3m_sum_global = NULL;

  if (mpi_rank == 0) {
    np_sum_global  =  new double[nbins];
    nm_sum_global  =  new double[nbins];
    mp_sum_global  =  new double[nbins];
    mm_sum_global  =  new double[nbins];
    mup_sum_global  =  new double[nbins];
    mum_sum_global  =  new double[nbins];
    d2p_sum_global =  new double[nbins];
    d3p_sum_global =  new double[nbins];
    d2m_sum_global =  new double[nbins];
    d3m_sum_global =  new double[nbins];
  }

  MPI_Reduce(np_sum,np_sum_global,nbins,MPI_DOUBLE,MPI_SUM,0,mpi_comm);
  MPI_Reduce(nm_sum,nm_sum_global,nbins,MPI_DOUBLE,MPI_SUM,0,mpi_comm);
  MPI_Reduce(mp_sum,mp_sum_global,nbins,MPI_DOUBLE,MPI_SUM,0,mpi_comm);
  MPI_Reduce(mm_sum,mm_sum_global,nbins,MPI_DOUBLE,MPI_SUM,0,mpi_comm);
  MPI_Reduce(mup_sum,mup_sum_global,nbins,MPI_DOUBLE,MPI_SUM,0,mpi_comm);
  MPI_Reduce(mum_sum,mum_sum_global,nbins,MPI_DOUBLE,MPI_SUM,0,mpi_comm);
  MPI_Reduce(d2p_sum,d2p_sum_global,nbins,MPI_DOUBLE,MPI_SUM,0,mpi_comm);
  MPI_Reduce(d2m_sum,d2m_sum_global,nbins,MPI_DOUBLE,MPI_SUM,0,mpi_comm);
  MPI_Reduce(d3p_sum,d3p_sum_global,nbins,MPI_DOUBLE,MPI_SUM,0,mpi_comm);
  MPI_Reduce(d3m_sum,d3m_sum_global,nbins,MPI_DOUBLE,MPI_SUM,0,mpi_comm);

  if (mpi_rank == 0) {

    char filename[128];
    buildIndexedFilename(filename,(name+".lspstats").c_str(),step,"dat");
    mkdir_for_file(filename); // allow the user to have specified a subdirectory...
    COUT1(" > writeLSPStats: " << filename);
    FILE * fp = fopen(filename,"w");
    
    fprintf(fp,"%s",param_line.c_str());
    fprintf(fp,"# simulation time: %8.16le , stats average time:  %8.16le\n",time,dt_sum);
    fprintf(fp,"title = \"lsp.stats\"\n");
    fprintf(fp,"variables = \"X\"\n\"Y\"\n\"Z\"\n\"R\"\n\"THETA\"\n\"Np\"\n\"Nm\"\n\"N\"\n\"Mp\"\n\"Mm\"\n\"M\"\n\"Vp\"\n\"Vm\"\n\"V\"\n\"MUp\"\n\"MUn\"\n\"MU\"\n\"SMDp\"\n\"SMDm\"\n\"SMD\"\n");
    fprintf(fp,"zone i = %d, j = %d, DATAPACKING=POINT\n",nr,ntheta);

    double mdotp = 0.0;
    double mdotm = 0.0;
    double mdot  = 0.0;
    double vdotp = 0.0;
    double vdotm = 0.0;
    double vdot = 0.0;
    double totalarea = 0.0;
    double np = 0.0;
    double nm = 0.0;

    for (int itheta = 0; itheta < ntheta; ++itheta) {
      for (int ir = 0; ir < nr; ++ir) {
        const int bin = ir*ntheta + itheta;
        double r = (double(ir)+0.5)/double(nr)*(r1-r0)+r0;
        double theta = (double(itheta)+0.5)/double(ntheta)*(theta1-theta0)+theta0;
        // coords...
        double r_cos_theta = r*cos(theta);
        double r_sin_theta = r*sin(theta);
        double x = xc[0] + r_cos_theta*t1[0] + r_sin_theta*t2[0];
        double y = xc[1] + r_cos_theta*t1[1] + r_sin_theta*t2[1];
        double z = xc[2] + r_cos_theta*t1[2] + r_sin_theta*t2[2];
        double area = r*(r1-r0)/double(nr)*(theta1-theta0)/double(ntheta);
        fprintf(fp,"%18.15le %18.15le %18.15le %18.15le %18.15le ",x,y,z,r,theta*180.0/M_PI);
        fprintf(fp,"%18.15le %18.15le %18.15le ",np_sum_global[bin]/dt_sum/area,-nm_sum_global[bin]/dt_sum/area,(np_sum_global[bin]-nm_sum_global[bin])/dt_sum/area);
        fprintf(fp,"%18.15le %18.15le %18.15le ",mp_sum_global[bin]/dt_sum/area,-mm_sum_global[bin]/dt_sum/area,(mp_sum_global[bin]-mm_sum_global[bin])/dt_sum/area);
        fprintf(fp,"%18.15le %18.15le %18.15le ",M_PI/6.0*d3p_sum_global[bin]/dt_sum/area,-M_PI/6.0*d3m_sum_global[bin]/dt_sum/area,M_PI/6.0*(d3p_sum_global[bin]-d3m_sum_global[bin])/dt_sum/area);
        fprintf(fp,"%18.15le %18.15le %18.15le ",mup_sum_global[bin]/dt_sum/area,-mum_sum_global[bin]/dt_sum/area,(mup_sum_global[bin]-mum_sum_global[bin])/dt_sum/area);

        if (d2p_sum_global[bin]+d2m_sum_global[bin] == 0.0)
          fprintf(fp,"0.0 0.0 0.0\n");
        else if (d2p_sum_global[bin] == 0.0)
          fprintf(fp,"0.0 %18.15le %18.15le\n",d3m_sum_global[bin]/d2m_sum_global[bin],d3m_sum_global[bin]/d2m_sum_global[bin]);
        else if (d2m_sum_global[bin] == 0.0)
          fprintf(fp,"%18.15le 0.0 %18.15le\n",d3p_sum_global[bin]/d2p_sum_global[bin],d3p_sum_global[bin]/d2p_sum_global[bin]);
        else
          fprintf(fp,"%18.15le %18.15le %18.15le\n",d3p_sum_global[bin]/d2p_sum_global[bin],d3m_sum_global[bin]/d2m_sum_global[bin],(d3p_sum_global[bin]+d3m_sum_global[bin])/(d2p_sum_global[bin]+d2m_sum_global[bin]));


        mdotp += mp_sum_global[bin];
        mdotm -= mm_sum_global[bin];

        vdotp += M_PI/6.0*d3p_sum_global[bin];
        vdotm -= M_PI/6.0*d3m_sum_global[bin];
        totalarea += area;

        np += np_sum_global[bin];
        nm += nm_sum_global[bin];
      }
    }

    fclose(fp);

    mdotp /= dt_sum;
    mdotm /= dt_sum;
    mdot   = mdotp + mdotm;
    vdotp /= dt_sum;
    vdotm /= dt_sum;
    vdot   = vdotp + vdotm;

    COUT1(" > particle mass through plane " << name << ": mdotp = " << mdotp << "  mdotm = " << mdotm << "  mdot = " << mdot);
    COUT1(" > particle volume flux through plane " << name << ": vdotp = " << vdotp/totalarea << "  vdotm = " << vdotm/totalarea << "  vdot = " << vdot/totalarea);
    COUT1(" > total area and dt_sum " << name << ": area = " << totalarea << " " <<"dt_sum =" << " " << dt_sum );
    COUT1(" > total number of particles " << name << ": np = " << np << " " << " nm = " << nm );

    double SUM_d3m = 0.0;
    double SUM_d3p = 0.0;
    double SUM_d2m = 0.0;
    double SUM_d2p = 0.0;
    for (int bin = 0; bin < nbins; bin++) {
      SUM_d3m += d3m_sum_global[bin];
      SUM_d2m += d2m_sum_global[bin];
      SUM_d3p += d3p_sum_global[bin];
      SUM_d2p += d2p_sum_global[bin];
    }
    if (SUM_d2p+SUM_d2m == 0.0) {
      COUT1(" > particle SMDp: " << "0.0 " << "SMDm: " << "0.0 " << "SMD: " << "0.0");
    } else if (SUM_d2p == 0.0) {
      COUT1(" > particle SMDp: " << "0.0 " << "SMDm: " << SUM_d3m/SUM_d2m << " SMD: " << SUM_d3m/SUM_d2m);
    } else if (SUM_d2m == 0.0) {
      COUT1(" > particle SMDp: " << SUM_d3p/SUM_d2p << " SMDm: " << "0.0 " << "SMD: " << SUM_d3p/SUM_d2p);
    } else {
      COUT1(" > particle SMDp: " << SUM_d3p/SUM_d2p << " SMDm: " << SUM_d3m/SUM_d2m << " SMD: " << (SUM_d3p+SUM_d3m)/(SUM_d2p+SUM_d2m));
    }

    delete[] np_sum_global;
    delete[] nm_sum_global;
    delete[] mp_sum_global;
    delete[] mm_sum_global;
    delete[] mup_sum_global;
    delete[] mum_sum_global;
    delete[] d2p_sum_global;
    delete[] d3p_sum_global;
    delete[] d2m_sum_global;
    delete[] d3m_sum_global;
  
  }

  if (got_pdf) LspStats::writePdf(time,step,verbose);
}

void LspStatsRect::write(const double time,const int step,const bool verbose) {

  // reduce to rank0 and plot...

  double * np_sum_global = NULL;
  double * nm_sum_global = NULL;
  double * mp_sum_global = NULL;
  double * mm_sum_global = NULL;
  double * mup_sum_global = NULL;
  double * mum_sum_global = NULL;
  double * d2p_sum_global = NULL;
  double * d3p_sum_global = NULL;
  double * d2m_sum_global = NULL;
  double * d3m_sum_global = NULL;

  if (mpi_rank == 0) {
    np_sum_global  =  new double[nbins];
    nm_sum_global  =  new double[nbins];
    mp_sum_global  =  new double[nbins];
    mm_sum_global  =  new double[nbins];
    mup_sum_global  =  new double[nbins];
    mum_sum_global  =  new double[nbins];
    d2p_sum_global =  new double[nbins];
    d3p_sum_global =  new double[nbins];
    d2m_sum_global =  new double[nbins];
    d3m_sum_global =  new double[nbins];
  }

  MPI_Reduce(np_sum,np_sum_global,nbins,MPI_DOUBLE,MPI_SUM,0,mpi_comm);
  MPI_Reduce(nm_sum,nm_sum_global,nbins,MPI_DOUBLE,MPI_SUM,0,mpi_comm);
  MPI_Reduce(mp_sum,mp_sum_global,nbins,MPI_DOUBLE,MPI_SUM,0,mpi_comm);
  MPI_Reduce(mm_sum,mm_sum_global,nbins,MPI_DOUBLE,MPI_SUM,0,mpi_comm);
  MPI_Reduce(mup_sum,mup_sum_global,nbins,MPI_DOUBLE,MPI_SUM,0,mpi_comm);
  MPI_Reduce(mum_sum,mum_sum_global,nbins,MPI_DOUBLE,MPI_SUM,0,mpi_comm);
  MPI_Reduce(d2p_sum,d2p_sum_global,nbins,MPI_DOUBLE,MPI_SUM,0,mpi_comm);
  MPI_Reduce(d2m_sum,d2m_sum_global,nbins,MPI_DOUBLE,MPI_SUM,0,mpi_comm);
  MPI_Reduce(d3p_sum,d3p_sum_global,nbins,MPI_DOUBLE,MPI_SUM,0,mpi_comm);
  MPI_Reduce(d3m_sum,d3m_sum_global,nbins,MPI_DOUBLE,MPI_SUM,0,mpi_comm);

  if (mpi_rank == 0) {

    char filename[128];
    buildIndexedFilename(filename,(name+".lspstats").c_str(),step,"dat");
    mkdir_for_file(filename); // allow the user to have specified a subdirectory...
    COUT1(" > writeLSPStats: " << filename);
    FILE * fp = fopen(filename,"w");
    
    fprintf(fp,"%s",param_line.c_str());
    fprintf(fp,"# simulation time: %8.16le , stats average time:  %8.16le\n",time,dt_sum);
    fprintf(fp,"title = \"lsp.stats\"\n");
    fprintf(fp,"variables = \"X\"\n\"Y\"\n\"Z\"\n\"Np\"\n\"Nm\"\n\"N\"\n\"Mp\"\n\"Mm\"\n\"M\"\n\"Vp\"\n\"Vm\"\n\"V\"\n\"MUp\"\n\"MUn\"\n\"MU\"\n\"SMDp\"\n\"SMDm\"\n\"SMD\"\n");
    fprintf(fp,"zone i = %d, j = %d, DATAPACKING=POINT\n",nbins_s1,nbins_s2);

    double mdotp = 0.0;
    double mdotm = 0.0;
    double mdot  = 0.0;
    double vdotp = 0.0;
    double vdotm = 0.0;
    double vdot = 0.0;
    double totalarea = 0.0;
    double np = 0.0;
    double nm = 0.0;

    for (int ix1 = 0; ix1 < nbins_s1; ++ix1) {
      for (int ix2 = 0; ix2 < nbins_s2; ++ix2) {
        const int bin = ix2*nbins_s1 + ix1;
        const double dx1 = (double(ix1)/double(nbins_s1) - 0.5)*s1; 
        const double dx2 = (double(ix2)/double(nbins_s2) - 0.5)*s2; 
        // coords...
        double x[3]; 
        FOR_I3 x[i] = xc[i] + dx1*t1[i] + dx2*t2[i];
        double area = (s1/double(nbins_s1))*(s2/double(nbins_s2));
        fprintf(fp,"%18.15le %18.15le %18.15le ",x[0],x[1],x[2]);
        fprintf(fp,"%18.15le %18.15le %18.15le ",np_sum_global[bin]/dt_sum/area,-nm_sum_global[bin]/dt_sum/area,(np_sum_global[bin]-nm_sum_global[bin])/dt_sum/area);
        fprintf(fp,"%18.15le %18.15le %18.15le ",mp_sum_global[bin]/dt_sum/area,-mm_sum_global[bin]/dt_sum/area,(mp_sum_global[bin]-mm_sum_global[bin])/dt_sum/area);
        fprintf(fp,"%18.15le %18.15le %18.15le ",M_PI/6.0*d3p_sum_global[bin]/dt_sum/area,-M_PI/6.0*d3m_sum_global[bin]/dt_sum/area,M_PI/6.0*(d3p_sum_global[bin]-d3m_sum_global[bin])/dt_sum/area);
        fprintf(fp,"%18.15le %18.15le %18.15le ",mup_sum_global[bin]/dt_sum/area,-mum_sum_global[bin]/dt_sum/area,(mup_sum_global[bin]-mum_sum_global[bin])/dt_sum/area);

        if (d2p_sum_global[bin]+d2m_sum_global[bin] == 0.0)
          fprintf(fp,"0.0 0.0 0.0\n");
        else if (d2p_sum_global[bin] == 0.0)
          fprintf(fp,"0.0 %18.15le %18.15le\n",d3m_sum_global[bin]/d2m_sum_global[bin],d3m_sum_global[bin]/d2m_sum_global[bin]);
        else if (d2m_sum_global[bin] == 0.0)
          fprintf(fp,"%18.15le 0.0 %18.15le\n",d3p_sum_global[bin]/d2p_sum_global[bin],d3p_sum_global[bin]/d2p_sum_global[bin]);
        else
          fprintf(fp,"%18.15le %18.15le %18.15le\n",d3p_sum_global[bin]/d2p_sum_global[bin],d3m_sum_global[bin]/d2m_sum_global[bin],(d3p_sum_global[bin]+d3m_sum_global[bin])/(d2p_sum_global[bin]+d2m_sum_global[bin]));


        mdotp += mp_sum_global[bin];
        mdotm -= mm_sum_global[bin];

        vdotp += M_PI/6.0*d3p_sum_global[bin];
        vdotm -= M_PI/6.0*d3m_sum_global[bin];
        totalarea += area;

        np += np_sum_global[bin];
        nm += nm_sum_global[bin];
      }
    }

    fclose(fp);

    mdotp /= dt_sum;
    mdotm /= dt_sum;
    mdot   = mdotp + mdotm;
    vdotp /= dt_sum;
    vdotm /= dt_sum;
    vdot   = vdotp + vdotm;

    COUT1(" > particle mass through plane " << name << ": mdotp = " << mdotp << "  mdotm = " << mdotm << "  mdot = " << mdot);
    COUT1(" > particle volume flux through plane " << name << ": vdotp = " << vdotp/totalarea << "  vdotm = " << vdotm/totalarea << "  vdot = " << vdot/totalarea);
    COUT1(" > total area and dt_sum " << name << ": area = " << totalarea << " " <<"dt_sum =" << " " << dt_sum );
    COUT1(" > total number of particles " << name << ": np = " << np << " " << " nm = " << nm );

    double SUM_d3m = 0.0;
    double SUM_d3p = 0.0;
    double SUM_d2m = 0.0;
    double SUM_d2p = 0.0;
    for (int bin = 0; bin < nbins; bin++) {
      SUM_d3m += d3m_sum_global[bin];
      SUM_d2m += d2m_sum_global[bin];
      SUM_d3p += d3p_sum_global[bin];
      SUM_d2p += d2p_sum_global[bin];
    }
    if (SUM_d2p+SUM_d2m == 0.0) {
      COUT1(" > particle SMDp: " << "0.0 " << "SMDm: " << "0.0 " << "SMD: " << "0.0");
    } else if (SUM_d2p == 0.0) {
      COUT1(" > particle SMDp: " << "0.0 " << "SMDm: " << SUM_d3m/SUM_d2m << " SMD: " << SUM_d3m/SUM_d2m);
    } else if (SUM_d2m == 0.0) {
      COUT1(" > particle SMDp: " << SUM_d3p/SUM_d2p << " SMDm: " << "0.0 " << "SMD: " << SUM_d3p/SUM_d2p);
    } else {
      COUT1(" > particle SMDp: " << SUM_d3p/SUM_d2p << " SMDm: " << SUM_d3m/SUM_d2m << " SMD: " << (SUM_d3p+SUM_d3m)/(SUM_d2p+SUM_d2m));
    }

    delete[] np_sum_global;
    delete[] nm_sum_global;
    delete[] mp_sum_global;
    delete[] mm_sum_global;
    delete[] mup_sum_global;
    delete[] mum_sum_global;
    delete[] d2p_sum_global;
    delete[] d3p_sum_global;
    delete[] d2m_sum_global;
    delete[] d3m_sum_global;
  
  }

  if (got_pdf) LspStats::writePdf(time,step,verbose);
}

LspStats * newLspStats(Param * param) {
  
  string geomType;
  bool got_geom = false;

  for (int iarg = 0; iarg < param->size(); iarg++) {
    string token = param->getString(iarg);
    if (token == "GEOM") {
      got_geom = true;
      geomType = param->getString(iarg+1);
      break;
    }
  }

  if (not got_geom) CERR("could not find GEOM in LSP.STATS!");

  if (geomType == "DISK") {
    LspStats * stats = new LspStatsDisk(param);
    stats -> init(param);
    return stats;
  } else if (geomType == "RECTANGLE") {
    LspStats * stats = new LspStatsRect(param);
    stats -> init(param);
    return stats;
  } else {
    CERR("in LSP.STATS unsupported GEOM: "+geomType);
  }
  
  return NULL; // should not end up here
}

LspPDF::LspPDF() {
  pdf = NULL;
  diam = NULL;
}

LspPDF::~LspPDF() {
  DELETE(pdf);
  DELETE(diam);
}


void LspPDF::init(Param * param) {

  // a typical stats line...
  // LSP.PDF NAME=r4 INTERVAL=10 START=10 GEOM=BOX <xmin> <xmax> <ymin> <ymax> <zmin> <zmax> NBINS <nbins> DRANGE <dmin> <dmax>
  // LSP.PDF NAME=r4 INTERVAL=10 START=10 GEOM=CONE X <x> <y> <z> N <nx> <ny> <nz> H <h> R0 <r0> R1 <r1> NBINS <nbins> DRANGE <dmin> <dmax>

  bool got_name = false;
  bool got_interval = false;
  bool got_start = false;
  bool got_geom = false;

  int iarg = 0;
  while (iarg < param->size()) {
    string token = param->getString(iarg++);
    if (token == "NAME") {
      assert(!got_name);
      got_name = true;
      name = param->getString(iarg++);
    } else if (token == "INTERVAL") {
      assert(!got_interval);
      got_interval = true;
      interval = param->getInt(iarg++);
      assert(interval >= 1);
    } else if (token == "START") {
      assert(!got_start);
      got_start = true;
      start = param->getInt(iarg++);
      assert(start >= 1);
    } else if (token == "GEOM") {
      assert(!got_geom);
      got_geom = true;
      token = param->getString(iarg++);
      if (token == "BOX") {
        geom = BOX;
        xmin[0] = param->getDouble(iarg++);
        xmax[0] = param->getDouble(iarg++);
        xmin[1] = param->getDouble(iarg++);
        xmax[1] = param->getDouble(iarg++);
        xmin[2] = param->getDouble(iarg++);
        xmax[2] = param->getDouble(iarg++);
      } else if (token == "CONE") { 
        geom = CONE;
        assert(param->getString(iarg) == "X");
        ++iarg;
        xc[0] = param->getDouble(iarg++);
        xc[1] = param->getDouble(iarg++);
        xc[2] = param->getDouble(iarg++);
        assert(param->getString(iarg) == "N");
        ++iarg;
        normal[0] = param->getDouble(iarg++);
        normal[1] = param->getDouble(iarg++);
        normal[2] = param->getDouble(iarg++);
        // make this a unit normal...
        const double mag = MAG(normal);
        assert(mag > 0.0);
        FOR_J3 normal[j] /= mag;
        assert(param->getString(iarg) == "H");
        ++iarg;
        hight = param->getDouble(iarg++);
        assert(param->getString(iarg) == "R0");
        ++iarg;
        r0 = param->getDouble(iarg++);
        assert(param->getString(iarg) == "R1");
        ++iarg;
        r1 = param->getDouble(iarg++);
      } else {
        CERR("unsupported GEOM: " << token);
      }
    } else if (token == "NBINS") {
        nbins = param->getInt(iarg++);
    } else if (token == "DRANGE") {
        dmin = param->getDouble(iarg++);
        dmax = param->getDouble(iarg++);
    } else {
      CERR("unrecognized token: " << token);
    }
  }

  int ierr = 0;
  if (!got_name) {
    if (mpi_rank == 0)
      cerr << "Error: LSP.PDF missing NAME" << endl;
    ierr = -1;
  }
  if (!got_interval) {
    if (mpi_rank == 0)
      cerr << "Error: LSP.PDF missing INTERVAL" << endl;
    ierr = -1;
  }
  if (!got_start) {
    if (mpi_rank == 0)
      cerr << "Error: LSP.PDF missing START" << endl;
    ierr = -1;
  }
  if (!got_geom) {
    if (mpi_rank == 0)
      cerr << "Error: LSP.PDF missing GEOM" << endl;
    ierr = -1;
  }

  if (ierr != 0)
    throw(0);

  if (mpi_rank == 0)
    cout << "LSP.PDF: " << name << " nbins=" << nbins << endl;

  pdf   =  new double[nbins];
  diam  =  new double[nbins];

  // and zero everything. Here we could read from restart file some day...

  for (int ibin = 0; ibin < nbins; ++ibin) {
    diam[ibin]  =  getD(ibin);
    pdf[ibin]  =  0.0;
  }

}

int LspPDF::getBin(const double d) {

  const int d_bin = (int)floor( (d-dmin)/(dmax-dmin)*double(nbins) );

  if ((d_bin < 0)||(d_bin >= nbins)) {
    return(-1);
  }

  return(d_bin);

}

double LspPDF::getD(const int ibin) {
  const double dD  =  (dmax-dmin)/nbins;
  const double d   =  dmin + (double(ibin) + 0.5)*dD;
  return(d);
}

bool LspPDF::inside(const LspState& lsp) {
  bool is_inside = false;
  const double* xp = lsp.xp;
  switch(geom) {
    case BOX:
    if ((xp[0] >= xmin[0]) && (xp[0] <= xmax[0]))
      if ((xp[1] >= xmin[1]) && (xp[1] <= xmax[1]))
        if ((xp[2] >= xmin[2]) && (xp[2] <= xmax[2]))
          is_inside = true;
    break;
    case CONE:
      double xc1[3]; // the center of the other base
      FOR_I3 xc1[i] = xc[i] + normal[i]*hight;
      const double N[3] = TRI_NORMAL_2(xc,xc1,xp); // reference normal,
      const double dir[3] = CROSS_PRODUCT(N,normal);
      const double magdir = MAG(dir);
      
      double xR[3], xR1[3];
      FOR_I3 xR[i] = xc[i] + r0*dir[i]/magdir;
      FOR_I3 xR1[i] = xc1[i] + r1*dir[i]/magdir;

      const double n0[3] = TRI_NORMAL_2(xc1,xR1,xp);
      if (DOT_PRODUCT(N,n0) >= 0) {
        const double n1[3] = TRI_NORMAL_2(xR1,xR,xp);
        if (DOT_PRODUCT(N,n1) >= 0) {
          const double n2[3] = TRI_NORMAL_2(xR,xc,xp);
          if (DOT_PRODUCT(N,n2) >= 0) {
            is_inside = true;
          }
        }
      }
    break;
  }

  return is_inside;
}

void LspPDF::update(const LspState * lsp,const int np,const double dt,const bool verbose) {
  reset(verbose);

  for (int ip = 0; ip < np; ++ip) {

    // check if the particle inside...
    
    if (inside(lsp[ip])) {
      const int ibin = getBin(lsp[ip].dp);
      if (ibin < 0) continue;
      assert(ibin<nbins);
      pdf[ibin] += lsp[ip].mp;
    } else { continue;}
    
  }

}

void LspPDF::reset(const bool verbose) {
  for (int ibin = 0; ibin < nbins; ++ibin) 
    pdf[ibin] = 0.0;
}

void LspPDF::write(const double time,const int step,const bool verbose) {

  // reduce the pdf to rank 0 then write

  double * pdf_global = NULL;

  IF_RANK0 {
    pdf_global = new double [nbins];
  }

  MPI_Reduce(pdf,pdf_global,nbins,MPI_DOUBLE,MPI_SUM,0,mpi_comm);

  IF_RANK0 {
    const double dD = (dmax-dmin)/nbins;
    double sum_mass = 0.0;
    for (int ibin = 0; ibin < nbins; ++ibin)
      sum_mass += pdf_global[ibin];
    if (sum_mass == 0.0) {
      IF_RANK0 cout << "No particles present in the geometry, skipping writing pdf file!" << endl;
      return;
    }
    for (int ibin = 0; ibin < nbins; ++ibin)
      pdf_global[ibin] /= (dD*sum_mass);
    double sum_pdf = 0.0;
    for (int ibin = 0; ibin < nbins; ++ibin) 
      sum_pdf += pdf_global[ibin]*dD;
    cout << "sum_pdf-1.0 (should be 0): " << sum_pdf-1.0 << endl;
    assert((sum_pdf-1.0)<1e-12);

    char filename[128];
    buildIndexedFilename(filename,(name+".lsppdf").c_str(),step,"dat");
    mkdir_for_file(filename); // allow the user to have specified a subdirectory...

    COUT1("writeLSPPDF: " << filename);

    FILE * fp = fopen(filename,"w");

    fprintf(fp,"title: \"mass PDF\"\n");
    fprintf(fp,"variables: \"D\"\n \"PDF\"\n");
    fprintf(fp,"zone: i=%d, DATAPACKING=POINT\n",nbins);
    fprintf(fp,"total mass: %18.15le\n",sum_mass);
  
    for (int ibin = 0; ibin < nbins; ++ibin) {
      fprintf(fp,"%18.15le %18.15le\n",diam[ibin],pdf_global[ibin]);
    }
    fclose(fp);
    
  }

  delete[] pdf_global;

}

void LspPDF::report() {
  IF_RANK0 {
    cout << "name: " << name << endl
        << "interval: " << interval << endl
        << "geometry: " << geom << " , (BOX: 0, CONE: 1) " << endl
        << "xmin: " << COUT_VEC(xmin) << " , xmax: " << COUT_VEC(xmax) << endl
        << "xc: " << COUT_VEC(xc) << " , normal: " << COUT_VEC(normal) << " , hight: " << hight << endl
        << "r0: " << r0 << " , r1: " << r1 << endl
        << "nbins: " << nbins << " , dmin: " << dmin << " , dmax: " << dmax << endl;
  }
}


