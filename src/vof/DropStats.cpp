#include "DropStats.hpp"
#include "CtiLiquid.hpp"

DropPDF::DropPDF() {
  pdf = NULL;
  num = NULL;
  diam = NULL;
}

DropPDF::~DropPDF() {
  DELETE(pdf);
  DELETE(num);
  DELETE(diam);
}


void DropPDF::init(Param * param) {

  // a typical stats line...
  // LSP.PDF NAME=r4 INTERVAL=10 START=10 GEOM=BOX <xmin> <xmax> <ymin> <ymax> <zmin> <zmax> NBINS <nbins> DRANGE <dmin> <dmax>
  // LSP.PDF NAME=r4 INTERVAL=10 START=10 GEOM=CONE X <x> <y> <z> N <nx> <ny> <nz> H <h> R0 <r0> R1 <r1> NBINS <nbins> DRANGE <dmin> <dmax>

  bool got_name = false;
  bool got_interval = false;
  bool got_start = false;
  bool got_geom = false;

  int i = 0;
  while (i < param->size()) {
    string token = param->getString(i++);
    if (token == "NAME") {
      assert(!got_name);
      got_name = true;
      name = param->getString(i++);
    } else if (token == "INTERVAL") {
      assert(!got_interval);
      got_interval = true;
      interval = param->getInt(i++);
      assert(interval >= 1);
    } else if (token == "START") {
      assert(!got_start);
      got_start = true;
      start = param->getInt(i++);
      assert(start >= 1);
    } else if (token == "GEOM") {
      assert(!got_geom);
      got_geom = true;
      token = param->getString(i++);
      if (token == "BOX") {
        geom = BOX;
        xmin[0] = param->getDouble(i++);
        xmax[0] = param->getDouble(i++);
        xmin[1] = param->getDouble(i++);
        xmax[1] = param->getDouble(i++);
        xmin[2] = param->getDouble(i++);
        xmax[2] = param->getDouble(i++);
      } else if (token == "CONE") { 
        geom = CONE;
        assert(param->getString(i++) == "X");
        xc[0] = param->getDouble(i++);
        xc[1] = param->getDouble(i++);
        xc[2] = param->getDouble(i++);
        assert(param->getString(i++) == "N");
        normal[0] = param->getDouble(i++);
        normal[1] = param->getDouble(i++);
        normal[2] = param->getDouble(i++);
        // make this a unit normal...
        const double mag = MAG(normal);
        assert(mag > 0.0);
        FOR_J3 normal[j] /= mag;
        assert(param->getString(i++) == "H");
        hight = param->getDouble(i++);
        assert(param->getString(i++) == "R0");
        r0 = param->getDouble(i++);
        assert(param->getString(i++) == "R1");
        r1 = param->getDouble(i++);
      } else {
        CERR("unsupported GEOM: " << token);
      }
    } else if (token == "NBINS") {
        nbins = param->getInt(i++);
    } else if (token == "DRANGE") {
        dmin = param->getDouble(i++);
        dmax = param->getDouble(i++);
    } else if (token == "VAR") {
      varname = param->getString(i++);
    } else {
      CERR("unrecognized token: " << token);
    }
  }

  int ierr = 0;
  if (!got_name) {
    if (mpi_rank == 0)
      cerr << "Error: DROP.PDF missing NAME" << endl;
    ierr = -1;
  }
  if (!got_interval) {
    if (mpi_rank == 0)
      cerr << "Error: DROP.PDF missing INTERVAL" << endl;
    ierr = -1;
  }
  if (!got_start) {
    if (mpi_rank == 0)
      cerr << "Error: DROP.PDF missing START" << endl;
    ierr = -1;
  }
  if (!got_geom) {
    if (mpi_rank == 0)
      cerr << "Error: DROP.PDF missing GEOM" << endl;
    ierr = -1;
  }

  if (ierr != 0)
    throw(0);

  if (mpi_rank == 0)
    cout << "DROP.PDF: " << name << " nbins=" << nbins << endl;

  pdf   =  new double[nbins];
  diam  =  new double[nbins];
  num = new double[nbins];

  // and zero everything. Here we could read from restart file some day...

  for (int ibin = 0; ibin < nbins; ++ibin) {
    diam[ibin]  =  getD(ibin);
    pdf[ibin]  =  0.0;
    num[ibin] = 0.0;
  }

}

int DropPDF::getBin(const double d) {

  const int d_bin = (int)floor( (d-dmin)/(dmax-dmin)*double(nbins) );

  if ((d_bin < 0)||(d_bin >= nbins)) {
    return(-1);
  }

  return(d_bin);

}

double DropPDF::getD(const int ibin) {
  const double dD  =  (dmax-dmin)/nbins;
  const double d   =  dmin + (double(ibin) + 0.5)*dD;
  return(d);
}

bool DropPDF::inside(const double xp[3]) {
  bool is_inside = false;
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

void DropPDF::update(const double* dropD, const double (*dropX)[3], const double (*dropU)[3],  const int ndrops, const double time, const double dt) {

  int my_ndrops = 0;
  //reset(1);

  for (int ib = 0; ib < ndrops; ++ib) {
    if (inside(dropX[ib])) {
      if (varname == "D") {
        const double dp = dropD[ib];
        const int ibin = getBin(dp);
        if (ibin < 0) continue;
        assert(ibin<nbins);
        pdf[ibin] += M_PI/6.0*pow(dp,3);
        num[ibin] += 1;
      }
      else if (varname == "U_X") {
        const double up = dropU[ib][0];
        const double dp = dropD[ib];
        const int ibin = getBin(up);
        if (ibin < 0) continue;
        assert(ibin<nbins);
        pdf[ibin] += M_PI/6.0*pow(dp,3);
        num[ibin] += 1;
      }
      else if (varname == "U_YZ") {
        const double vp = dropU[ib][1];
        const double wp = dropU[ib][2];
        const double ut = sqrt(vp*vp + wp*wp);
        const double dp = dropD[ib];
        const int ibin = getBin(ut);
        if (ibin < 0) continue;
        assert(ibin<nbins);
        pdf[ibin] += M_PI/6.0*pow(dp,3);
        num[ibin] += 1;
      }
      else if (varname == "U_R") {
        const double yc = 0.1197;
        const double zc = 0.1037;
        const double yp = dropX[ib][1] - yc;
        const double zp = dropX[ib][2] - zc;
        double theta = atan(fabs(zp)/fabs(yp));

        if ( yp>= 0.0 && zp >= 0.0) theta = theta;
        else if (yp< 0.0 && zp >= 0.0) theta = M_PI-theta;
        else if (yp< 0.0 && zp < 0.0) theta = M_PI+ theta;
        else if (yp>= 0.0 && zp < 0.0) theta = 2.0*M_PI-theta;
        else CERR("Critical Error");

        const double vp = dropU[ib][1];
        const double wp = dropU[ib][2];

        const double ur = vp*cos(theta) + wp*sin(theta);
        const double utheta = vp*sin(theta) + wp*cos(theta); 
        const double dp = dropD[ib];
        const int ibin = getBin(ur);
        if (ibin < 0) continue;
        assert(ibin<nbins);
        pdf[ibin] += M_PI/6.0*pow(dp,3);
        num[ibin] += 1;
      }
      else if (varname == "U_THETA") {
        const double yc = 0.1197;
        const double zc = 0.1037;
        const double yp = dropX[ib][1] - yc;
        const double zp = dropX[ib][2] - zc;
        double theta = atan(fabs(zp)/fabs(yp));

        if ( yp>= 0.0 && zp >= 0.0) theta = theta;
        else if (yp< 0.0 && zp >= 0.0) theta = M_PI-theta;
        else if (yp< 0.0 && zp < 0.0) theta = M_PI+ theta;
        else if (yp>= 0.0 && zp < 0.0) theta = 2.0*M_PI-theta;
        else CERR("Critical angle  Error");

        const double vp = dropU[ib][1];
        const double wp = dropU[ib][2];

        const double ur = vp*cos(theta) + wp*sin(theta);
        const double utheta = vp*sin(theta) + wp*cos(theta); 
        const double dp = dropD[ib];
        const int ibin = getBin(utheta);
        if (ibin < 0) continue;
        assert(ibin<nbins);
        pdf[ibin] += M_PI/6.0*pow(dp,3);
        num[ibin] += 1;
      }
      else CERR("unknown var");
    } 
    else continue;
  }
  
}

void DropPDF::update(const double* blob_surf, const double* blob_vol, const double (*blobX)[3], const int nblobs, const double time, const double dt,const bool verbose) {

  reset(verbose);

  double my_dropSurf = 0.0;
  double my_dropVol = 0.0;
  int my_ndrops = 0;

  // count vof blobs first ...
  for (int ib = 0; ib < nblobs; ++ib) {
    if (inside(blobX[ib])) {
      const double dp = pow(6.0*blob_vol[ib]/M_PI,1.0/3.0);
      const int ibin = getBin(dp);
      if (ibin < 0) continue;
      assert(ibin<nbins);
      pdf[ibin] += blob_vol[ib];
      num[ibin] += 1.0;
      my_dropSurf += M_PI*dp*dp;
      my_dropVol  += blob_vol[ib];
      my_ndrops += 1;
    }
    else {continue;}
  }
  
  double dropSurf = 0.0;
  double dropVol  = 0.0;
  int ndrops = 0;
  MPI_Reduce(&my_dropSurf,&dropSurf,1,MPI_DOUBLE,MPI_SUM,0,mpi_comm);
  MPI_Reduce(&my_dropVol,&dropVol,1,MPI_DOUBLE,MPI_SUM,0,mpi_comm);
  MPI_Reduce(&my_ndrops,&ndrops,1,MPI_INT,MPI_SUM,0,mpi_comm);

  if (mpi_rank == 0) cout << "Drop STATS = " <<  name << " time= " << time << " " <<  "ndrops = " << ndrops << " " <<   "SMD23 = " <<  6.0*dropVol/dropSurf << endl; 
  if (mpi_rank == 0 ) cout << "Time, Surface, Volume = " << time << " " << dropSurf << " " << dropVol << endl;
}



void DropPDF::update(const double* blob_surf, const double* blob_vol, const double (*blobX)[3], const int nblobs, 
    const LspState * lsp,const int np, CtiMaterial* fuel, const double time, const double dt,const bool verbose) {

  reset(verbose);


  double my_dropSurf = 0.0;
  double my_dropVol = 0.0;
  double my_lspSurf = 0.0;
  double my_lspVol = 0.0;
  int my_lsps = 0;
  int my_drops =0;
  // count vof blobs first ...
  for (int ib = 0; ib < nblobs; ++ib) {
    if (inside(blobX[ib])) {
      const double dp = pow(6.0*blob_vol[ib]/M_PI,1.0/3.0);
      const int ibin = getBin(dp);
      if (ibin < 0) continue;
      assert(ibin<nbins);
      pdf[ibin] += blob_vol[ib];
      num[ibin] += 1.0;
      my_dropSurf += M_PI*dp*dp;
      my_dropVol  += blob_vol[ib];
      my_drops += 1;
    }
    else {continue;}
  }
 
  
  for (int ip = 0; ip < np; ++ip) {
    if (inside(lsp[ip].xp)) {
      const int ibin = getBin(lsp[ip].dp);
      if (ibin < 0) continue;
      assert(ibin<nbins);
      const double rhop  = fuel->calcRho(lsp[ip].Tp);
      pdf[ibin] += lsp[ip].mp/rhop; // volume pdf
      num[ibin] += 1.0; // number pdf 
      my_lspSurf += M_PI*lsp[ip].dp*lsp[ip].dp;
      my_lspVol  += M_PI*pow(lsp[ip].dp,3.0)/6.0;
      my_lsps   += 1;
    } else { continue;}
    
  }
  double dropSurf = 0.0;
  double dropVol  = 0.0;
  double lspSurf = 0.0;
  double lspVol  = 0.0;

  int nlsps = 0;
  int ndrops = 0;
  int np_total = 0;
  MPI_Reduce(&my_dropSurf,&dropSurf,1,MPI_DOUBLE,MPI_SUM,0,mpi_comm);
  MPI_Reduce(&my_dropVol,&dropVol,1,MPI_DOUBLE,MPI_SUM,0,mpi_comm);
  MPI_Reduce(&my_lspSurf,&lspSurf,1,MPI_DOUBLE,MPI_SUM,0,mpi_comm);
  MPI_Reduce(&my_lspVol,&lspVol,1,MPI_DOUBLE,MPI_SUM,0,mpi_comm);
  MPI_Reduce(&my_lsps,&nlsps,1,MPI_INT,MPI_SUM,0,mpi_comm);
  MPI_Reduce(&my_drops,&ndrops,1,MPI_INT,MPI_SUM,0,mpi_comm);
  MPI_Reduce(&np,&np_total,1,MPI_INT,MPI_SUM,0,mpi_comm);


  if (mpi_rank == 0) cout << "Drop STATS = " <<  name << " time= " << time << " " << "nblobs= " << ndrops << " " <<  "nlsps= " << nlsps << " " <<   "SMD23 = " <<  6.0*(dropVol+lspVol)/(dropSurf+lspSurf) << endl; 
  if (mpi_rank == 0 ) cout << "Time, blob Surface, lsp Surface, blob Volume, lsp Volume = " << " " << time << " " << dropSurf << " " << lspSurf << " " <<  dropVol << " " << lspVol <<  endl;
}

void DropPDF::reset(const bool verbose) {
  for (int ibin = 0; ibin < nbins; ++ibin)  {
    pdf[ibin] = 0.0;
    num[ibin] = 0.0;
  }
 
}

void DropPDF::write(const double time,const int step,const bool verbose) {

  // reduce the pdf to rank 0 then write

  double * pdf_global = NULL;
  double * num_global = NULL;

  IF_RANK0 {
    pdf_global = new double [nbins];
    num_global = new double [nbins];
  }

  MPI_Reduce(pdf,pdf_global,nbins,MPI_DOUBLE,MPI_SUM,0,mpi_comm);
  MPI_Reduce(num,num_global,nbins,MPI_DOUBLE,MPI_SUM,0,mpi_comm);

  IF_RANK0 {
    const double dD = (dmax-dmin)/nbins;
    double sum_mass = 0.0;
    double sum_num = 0.0;
    for (int ibin = 0; ibin < nbins; ++ibin) {
      sum_mass += pdf_global[ibin];
      sum_num += num_global[ibin];
    }
    if (sum_mass == 0.0 || sum_num == 0.0 ) {
      IF_RANK0 cout << "No particles present in the geometry, skipping writing pdf file!" << endl;
      return;
    }
    for (int ibin = 0; ibin < nbins; ++ibin) {
      pdf_global[ibin] /= (dD*sum_mass);
      num_global[ibin] /= (dD*sum_num);
    }
    double sum_pdf = 0.0;
    double num_pdf = 0.0;
    for (int ibin = 0; ibin < nbins; ++ibin) {
      sum_pdf += pdf_global[ibin]*dD;
      num_pdf += num_global[ibin]*dD;
    }
    //cout << "sum_pdf-1.0 (should be 0): " << sum_pdf-1.0 << endl;
    assert((sum_pdf-1.0)<1e-12);
    assert((num_pdf-1.0)<1e-12);

    char filename[128];
    char filename2[128];
    buildIndexedFilename(filename,(name+".pdf").c_str(),step,"dat");
    buildIndexedFilename(filename2,(name+".num").c_str(),step,"dat");
    mkdir_for_file(filename); // allow the user to have specified a subdirectory...
    mkdir_for_file(filename2); // allow the user to have specified a subdirectory...

    COUT1("writeDropPDF: " << filename);
    COUT1("writeDropNumber: " << filename2);

    FILE * fp = fopen(filename,"w");

    //fprintf(fp,"VARIABLES =  "D," " PDF"\n");
  
    for (int ibin = 0; ibin < nbins; ++ibin) {
      fprintf(fp,"%18.15le %18.15le\n",diam[ibin],pdf_global[ibin]);
    }
    fclose(fp);

    FILE * fp2 = fopen(filename2,"w");

  
    for (int ibin = 0; ibin < nbins; ++ibin) {
      fprintf(fp2,"%18.15le %18.15le\n",diam[ibin],num_global[ibin]);
    }
    fclose(fp2);
    
  }

  delete[] pdf_global;
  delete[] num_global;

}

void DropPDF::report() {
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


