#include "LiquidFlux.hpp"

LiquidFlux::LiquidFlux() {
  //pointer null
  liquid_flux = NULL;
  yp = NULL;
  zp = NULL;
  projected_area = NULL;
}

LiquidFlux::~LiquidFlux() {
  DELETE(liquid_flux);
  DELETE(yp);
  DELETE(zp);
  DELETE(projected_area);
}

void LiquidFlux::init(Param * param) {

  // a typical stats line...
  // LIQUID.FLUX NAME=n1 SAMPLE_INTERVAL=10 WRITE_INTERVAL=100 XP <x> <y> <z> GEOM <ymin> <ymax> <zmin> <zmax> NY <ny> NZ <nz> 

  bool got_name = false;
  bool got_interval = false;
  bool got_geom = false;

  int i = 0;
  while (i < param->size()) {
    string token = param->getString(i++);
    if (token == "NAME") {
      assert(!got_name);
      got_name = true;
      name = param->getString(i++);
    } else if (token == "SAMPLE_INTERVAL") {
      assert(!got_interval);
      got_interval = true;
      sample_interval = param->getInt(i++);
      assert(sample_interval >= 1);
    } else if (token == "WRITE_INTERVAL") {
      write_interval = param->getInt(i++);
      assert(write_interval >= 1);
    } else if (token == "XP") {
      xp[0] = param->getDouble(i++);
      xp[1] = param->getDouble(i++);
      xp[2] = param->getDouble(i++);
    } else if (token == "GEOM") {
      assert(!got_geom);
      got_geom = true;
      ymin = param->getDouble(i++);
      ymax = param->getDouble(i++);
      zmin = param->getDouble(i++);
      zmax = param->getDouble(i++);
    } else if (token == "NY") {
      ny = param->getInt(i++);
    } else if (token == "NZ") {
      nz = param->getInt(i++);
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
  if (!got_geom) {
    if (mpi_rank == 0)
      cerr << "Error: DROP.PDF missing GEOM" << endl;
    ierr = -1;
  }

  if (ierr != 0)
    throw(0);

  nbins = ny*nz;
  if (mpi_rank == 0)
    cout << "LIQUID FLUX: " << name << " nbins=" << nbins << endl;

  liquid_flux   =  new double[nbins];
  projected_area  =  new double[nbins];

  //reset 
  for (int ibin=0; ibin < nbins; ibin++) {
    liquid_flux[ibin] = 0.0;
    projected_area[ibin] = 00.0;
    //cout << "init area = " << ibin << " " << projected_area[ibin] << endl;
  }

  yp = new double[ny];
  zp = new double[nz];

  const double dy = (ymax-ymin)/double(ny);
  const double dz = (zmax-zmin)/double(nz);

  for (int ybin=0; ybin<ny; ybin++) {
    yp[ybin] = ymin+0.5*dy  + ybin*dy;
  }

  for (int zbin=0; zbin<nz; zbin++) {
    zp[zbin] = zmin+0.5*dz  + zbin*dz;
  }


 // FOR_ICV {
 //   cout << x_cv[icv][0] << endl;
 // }

  // and zero everything. Here we could read from restart file some day...

}

void LiquidFlux::calcProjectedArea(double* my_area) {
  if (mpi_rank == 0 ) 
    cout << "calc projected area for each bin" << endl;


  MPI_Reduce(my_area,projected_area, nbins, MPI_DOUBLE, MPI_SUM, 0, mpi_comm);

  if (mpi_rank == 0 ) {
    double sum_area = 0.0;
    for (int ibin=0; ibin<nbins; ibin++){
      sum_area += projected_area[ibin];
      cout << "projected area = " << ibin << " " << projected_area[ibin] << endl;
    }
    cout << "total flux area = " << sum_area << endl;
  }

}

void LiquidFlux::calcLiquidFlux(double* local_flux) {
  if (mpi_rank == 0 ) 
    cout << "calc liquid flux for each bin" << endl;


  MPI_Reduce(local_flux, liquid_flux, nbins, MPI_DOUBLE, MPI_SUM, 0, mpi_comm);

  if (mpi_rank == 0 ) {
    double sum_flux = 0.0;
    for (int ibin=0; ibin<nbins; ibin++){
      sum_flux += liquid_flux[ibin];
  //    cout << "liquid flux = " << ibin << " " << liquid_flux[ibin] << endl;
    }
  //  cout << "total liquid flux = " << sum_flux << endl;
  }

}


void LiquidFlux::report() {

  if (mpi_rank == 0 ) 
    cout << "LIQUID FLUX "  << name << " " << ymin << " " << ymax << " " << zmin << " " << zmax << endl;

//  for (int ibin =0; ibin<nbins; ibin++) {
//   cout << "report " << ibin << " " <<  projected_area[ibin] << endl;
//  }


}

void LiquidFlux::write(const double time, const int step) {

  double * liquid_flux_global = NULL;
  double * projected_area_global = NULL;

  IF_RANK0 {
    liquid_flux_global = new double [nbins];
    projected_area_global = new double [nbins];
  }

  MPI_Reduce(liquid_flux,liquid_flux_global,nbins,MPI_DOUBLE,MPI_SUM,0,mpi_comm);
  MPI_Reduce(projected_area,projected_area_global,nbins,MPI_DOUBLE,MPI_SUM,0,mpi_comm);


  if (mpi_rank == 0 ) {

    char filename[128];
    buildIndexedFilename(filename,(name+".lmf").c_str(),step,"dat");
    mkdir_for_file(filename); // allow the user to have specified a subdirectory...

    COUT1("writeLiquidFlux: " << filename);

    FILE * fp = fopen(filename,"w");

    //fprintf(fp,"VARIABLES =  "D," " PDF"\n");
    fprintf(fp, "VARIABLES = X, Y, liquid_mass, area \n");
    fprintf(fp, "ZONE I=");
    fprintf(fp, "%6d", ny);
    fprintf(fp, ", J=");
    fprintf(fp, "%6d", nz);
    fprintf(fp, ",  DATAPACKING=POINT \n");

    for (int iz = 0; iz < nz; ++iz) {
      for (int iy=0; iy < ny; ++iy) {
        int ibin = iy + iz*ny;
        if (projected_area_global[ibin] > 1.0E-16) 
          fprintf(fp,"%18.7e %18.7e %18.7e %18.7e\n", yp[iy], zp[iz], liquid_flux_global[ibin]/projected_area_global[ibin], projected_area_global[ibin]);
        else
          fprintf(fp,"%18.7e %18.7e %18.7e %18.7e\n", yp[iy], zp[iz], 0.0, projected_area_global[ibin]);
      }
    }
    fclose(fp);
  }

  delete[] liquid_flux_global;
  delete[] projected_area_global;


}


/*
int LiquidFlux::getBin(const double d) {

  const int d_bin = (int)floor( (d-dmin)/(dmax-dmin)*double(nbins) );

  if ((d_bin < 0)||(d_bin >= nbins)) {
    return(-1);
  }

  return(d_bin);

}


bool LiquidFlux::inside(const double xp[3]) {
  bool is_inside = false;
  if ((xp[0] >= xmin[0]) && (xp[0] <= xmax[0]))
    if ((xp[1] >= xmin[1]) && (xp[1] <= xmax[1]))
      if ((xp[2] >= xmin[2]) && (xp[2] <= xmax[2]))
        is_inside = true;
  return is_inside;
}

*/

