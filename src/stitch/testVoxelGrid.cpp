#include "CTI.hpp"
using namespace CTI;

#include "VoxelGrid.hpp"

int main(int argc, char * argv[]) {

  try {
      
    CTI_Init(argc,argv,"testVoxelGrid.in");
    
    {
      VoxelGrid vg;
      SurfaceShm* surface = NULL; 
      int nx = 100;
      int nband = 10;
      bool b_write = false;
      if (Param* param = getParam("VOXEL_GRID")) {
        int iarg = 0;
        while (iarg < param->size()) {
          string token = param->getString(iarg++);
          if (token == "SURF") {
            token = param->getString(iarg++);
            if ((token == "SBIN")||(token == "BIN")) {
              const string filename = param->getString(iarg++);
              surface = new SurfaceShm(filename);
            }
          }
          else if (token == "NX") {
            nx = param->getInt(iarg++);
          }
          else if (token == "NBAND") {
            nband = param->getInt(iarg++);
          }
          else if (token == "WRITE") {
            b_write = true;
          }
        }
        double bbmin[3],bbmax[3];
        surface->getBbox(bbmin,bbmax);
        const double width = bbmax[0]-bbmin[0];
        const double dx = width/double(nx);
        const double band_width = double(nband)*dx;
        FOR_I3 {
          bbmin[i] -= band_width;
          bbmax[i] += band_width;
        }
        vg.init(bbmin,bbmax,nx+2*nband);
        vg.initVoxelDataFromSurface(surface,band_width,b_write);
        delete surface; surface = NULL;
      }

      if (Param* param = getParam("ORIENT")) {
        int maxiter = 1000;
        double alpha = 1.0; 
        double beta = 0.5;
        bool b_verbose = false;
        bool b_write = false;
        int iarg = 0;
        while (iarg < param->size()) {
          string token = param->getString(iarg++);
          if (token == "MAXITER") {
            maxiter = param->getInt(iarg++);
          }
          else if (token == "ALPHA") {
            alpha = param->getDouble(iarg++);
            if (alpha < 1.0)
              COUT1("Warning: alpha should be >= 1");
          }
          else if (token == "BETA") {
            beta = param->getDouble(iarg++);
            if (beta < 0.5)
              COUT1("Warning: beta should be >= 0.5");
          }
          else if (token == "VERBOSE") {
            b_verbose = true;
          }
          else if (token == "WRITE") {
            b_write = true;
          }
        }

        // repair surface using diffusion...
        vg.orientVoxelData(maxiter,alpha,beta,b_verbose,b_write);
      }

      if (Param* param = getParam("DIFFUSE")) {
        int maxiter = 1000;
        double zero = 1.0E-6;
        double relax = 1.0;
        bool b_verbose = false;
        bool b_write = false;
        int iarg = 0;
        while (iarg < param->size()) {
          string token = param->getString(iarg++);
          if (token == "MAXITER") {
            maxiter = param->getInt(iarg++);
          }
          else if (token == "ZERO") {
            zero = param->getDouble(iarg++);
          }
          else if (token == "RELAX") {
            relax = param->getDouble(iarg++);
          }
          else if (token == "VERBOSE") {
            b_verbose = true;
          }
          else if (token == "WRITE") {
            b_write = true;
          }
        }

        // repair surface using diffusion...
        vg.diffuseVoxelData(zero,relax,maxiter,b_verbose,b_write);
      }

      if (Param* param = getParam("REDISTANCE")) {
        int maxiter = 1000;
        double zero = 1.0E-6;
        double cfl = sqrt(3.0);
        bool b_verbose = false;
        bool b_write = false;
        int iarg = 0;
        while (iarg < param->size()) {
          string token = param->getString(iarg++);
          if (token == "MAXITER") {
            maxiter = param->getInt(iarg++);
          }
          else if (token == "ZERO") {
            zero = param->getDouble(iarg++);
          }
          else if (token == "CFL") {
            cfl = param->getDouble(iarg++);
          }
          else if (token == "VERBOSE") {
            b_verbose = true;
          }
          else if (token == "WRITE") {
            b_write = true;
          }
        }

        // repair surface using diffusion...
        vg.redistanceVoxelData(zero,cfl,maxiter,b_verbose,b_write);
      }

      // overwrite surface...
      if (Param* param = getParam("TRIANGULATE")) {
        double iso = 0.0;
        bool b_write = false;
        int iarg = 0;
        while (iarg < param->size()) {
          string token = param->getString(iarg++);
          if (token == "ISO") {
            iso = param->getDouble(iarg++);
          }
          else if (token == "WRITE") {
            b_write = true;
          }
        }
        vg.triangulateVoxelData(surface,iso,b_write);
      }

      if (Param* param = getParam("WRITE_SBIN")) {
        string name = "surface.sbin";
        int iarg = 0;
        while (iarg < param->size()) {
          string token = param->getString(iarg++);
          if (token == "NAME") {
            name = param->getString(iarg++);
          }
        }
        surface->writeBinary(name);

      }

      // cleanup...
      if (surface) delete surface;

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

