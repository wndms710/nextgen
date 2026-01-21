#include "StaticSolver.hpp"
#include "LspInjector.hpp"
#include "MiscUtils.hpp"
#include "SubSurface.hpp"

double error_function(const double eta);
double sampleNormalDist(const double mean, const double std, const double xmin, const double xmax);
void project_to_plane(double * x_pro, const double x[3], const double xp[3], const double np[3]);
bool CampFunction(pair<double, int> p1, pair<double, int> p2) { return (p1.second < p2.second); }

LspInjector::LspInjector() {
  IF_RANK0 cout << "LspInjector()" << endl;

  InjectorType = -1;
  name = "";
  Tp = 0.0;
  rhop = 0.0;
  mdot = 0.0;
  np = 0;

  DistType = -1;
  for (int i = 0; i < 5; ++i) DistData[i] = 0.0;

  cdfType = NA;
  countCdf_x = NULL;
  countCdf_y = NULL;
  countCdf_n = 500;
  c0 = c1 = -1;

  GeomType = -1;
  for (int i = 0; i < 11; ++i) GeomData[i] = 0.0;

  ResidualMass = 0.0;

  injectTime = -1.0;

  nst = 0;
  nsp = 0;
  spost = NULL;
  xp = NULL;
  cdf_area = NULL;

  FOR_I3 {
    e1[i] = 0;
    e2[i] = 0;
  }

  triVec.clear();
  cdf_area_tri = NULL;
  cdf_area_rank = NULL;

  solver = NULL;
  u = NULL;

  overWriteVelocity_b = false;
  FOR_I2 OWVData[i] = -1;

  pShape = -1;
  axis_std = axis_e = axis_f = -1;

}

LspInjector::~LspInjector() {
  IF_RANK0 cout << "~LspInjector()" << endl;
  
  if (spost != NULL) delete[] spost;
  if (xp != NULL) delete[] xp;
  if (cdf_area != NULL) delete[] cdf_area;
  if (countCdf_x != NULL) delete[] countCdf_x;
  if (countCdf_y != NULL) delete[] countCdf_y;
  if (cdf_area_tri != NULL) delete[] cdf_area_tri;
  if (cdf_area_rank != NULL) delete[] cdf_area_rank;
}

void LspInjector::initFromParams(Param * param, CtiMaterial * fuel, StaticSolver * solver, double (*u) [3]) {

  this->solver = solver;
  this->u = u;

  // add the injector using the Current Parameter...

  assert(param->getName() == "LSP.INJECTOR" );
  name = param->getString();

  if (mpi_rank == 0)
    cout << " > initializing LSP.INJECTOR: " << name << endl;

  bool TP_flag = false;

  pShape = SPHERE_SHAPE; // shpere is default, unless otherwise is specified in SHAPE

  int iarg = 1;
  while (iarg < param->size()) {
    string token = param->getString(iarg++);
    if (token == "TYPE") {
      // ==========================================
      // TYPE of injector:
      // ==========================================
      token = param->getString(iarg++);
      if (token == "COUNT") {
        assert( InjectorType == -1 );
        InjectorType = COUNT_TYPE;
        np = param->getInt(iarg++);
        if (mpi_rank == 0)
          cout << "   > injector TYPE = COUNT " << np << endl;
      }
      else if ( token == "MDOT" ) {
        assert( InjectorType == -1 );
        InjectorType = MDOT_TYPE;
        mdot =  param->getDouble(iarg++);
        if (mpi_rank == 0)
          cout << "   > injector TYPE = MDOT " << mdot << endl;
      }
      else {
        if (mpi_rank == 0)
          cerr << "\n\n\n*********************************************************\n" <<
            "Error: unrecognized TYPE in LSP.INJECTOR: " << token << ". Valid types are: \n" <<
            "TYPE COUNT <n>   -- injects n particles once\n" <<
            "TYPE MDOT <mdot> -- injects at rate mdot\n" <<
            "*********************************************************\n\n\n" << endl;
        throw(0);
      }
    }
    else if (token == "DIST") {
      // ==========================================
      // DIST of particles:
      // ==========================================
      token = param->getString(iarg++);
      if (token == "UNIFORM") {
        assert( DistType == -1 );
        DistType = UNIFORM_DIST;
        DistData[0] = param->getDouble(iarg++); // <diam>
        if (mpi_rank == 0)
          cout << "   > injector DIST = UNIFORM " << " Particle diam: " << DistData[0] << endl;
      }
      else if (token == "RR") {
        assert( DistType == -1 );
        DistType = RR_DIST;
        string token = param->getString(iarg++);
        if (token=="COUNT_DIST") {
          cdfType = CDF_COUNT_DIST;
        }
        else if (token=="MASS_DIST") {
          cdfType = CDF_MASS_DIST;
        } else {
          CERR("RR cdf distribution type is missing. Options: COUNT_DIST, MASS_DIST");
        }
        DistData[0] = param->getDouble(iarg++); // <dmin>
        DistData[1] = param->getDouble(iarg++); // <dmax>
        DistData[2] = param->getDouble(iarg++); // <dbar>
        if (DistData[1] <= DistData[0]) CERR("Wrong diameter values in LSP.INJECTOR->DIST");
        DistData[3] = param->getDouble(iarg++); // <nrr>
        if (mpi_rank == 0)
          cout << "   > injector DIST = RR, type: " << cdfType << 
            " (count_dist: " << CDF_COUNT_DIST << ", mass_dist: " << CDF_MASS_DIST << ")" << 
            " dmin: " << DistData[0] << " dmax: " << DistData[1] << " davg: " <<  DistData[2] <<
            " n: " << DistData[3] << endl;
      }
      else if (token == "GAUSSIAN") {
        assert( DistType == -1);
        DistType = GAUSSIAN_DIST;
        string token = param->getString(iarg++);
        if (token=="COUNT_DIST") {
          cdfType = CDF_COUNT_DIST;
        }
        else if (token=="MASS_DIST") {
          cdfType = CDF_MASS_DIST;
        } else {
          CERR("Gaussian cdf distribution type is missing. Options: COUNT_DIST, MASS_DIST");
        }
        DistData[0] = param->getDouble(iarg++); // min diameter
        DistData[1] = param->getDouble(iarg++); // max diameter
        DistData[2] = param->getDouble(iarg++); // mean diameter
        DistData[3] = param->getDouble(iarg++); // standard deviation
        if (DistData[1] <= DistData[0]) CERR("Wrong diameter values in LSP.INJECTOR->DIST");
        if (mpi_rank == 0)
          cout << "   > injector DIST = GAUSSIAN type: " << cdfType << 
            " (count_dist: " << CDF_COUNT_DIST << ", mass_dist: " << CDF_MASS_DIST << ")" << 
            " dmin: " << DistData[0] << " dmax: " << DistData[1] << " davg: " << DistData[2] <<
            " dstd: " << DistData[3] << endl;
      }
      else {
        if (mpi_rank == 0)
          cerr << "\n\n\n*********************************************************\n" <<
            "Error: unrecognized DIST in LSP.INJECTOR: " << token << ". Valid dists are: \n" <<
            "DIST UNIFORM <diam>                                            -- injects constant diameter particles\n" <<
            "DIST RR COUNT_DIST/MASS_DIST <dmin> <dmax> <dbar> <nrr>        -- Rosin-Rammler distribution\n" <<
            "DIST GAUSSIAN COUNT_DIST/MASS_DIST <dmin> <dmax> <davg> <std>  -- Gaussian distribution\n" <<
            "*********************************************************\n\n\n" << endl;
        throw(0);
      }
      build_countCdf_stuff();
    }
    else if (token == "GEOM") {
      // ==========================================
      // GEOM of particles:
      // ==========================================
      token = param->getString(iarg++);
      if (token == "POINT") {
        assert( GeomType == -1 );
        GeomType = POINT_GEOM;
        token = param->getString(iarg++);
        if (token != "X")
          CERR("LSP.INJECTOR expects X");
        GeomData[0] = param->getDouble(iarg++); // <x>
        GeomData[1] = param->getDouble(iarg++); // <y>
        GeomData[2] = param->getDouble(iarg++); // <z>
        token = param->getString(iarg++);
        if (token != "N")
          CERR("LSP.INJECTOR expects N");
        GeomData[3] = param->getDouble(iarg++);
        GeomData[4] = param->getDouble(iarg++);
        GeomData[5] = param->getDouble(iarg++);
        NORMALIZE((GeomData+3));
        token = param->getString(iarg++);
        if (token != "VMAG")
          CERR("LSP.INJECTOR expects VMAG");
        GeomData[6] = param->getDouble(iarg++);
        if (mpi_rank == 0)
          cout << "   > injector GEOM = POINT  X=" << GeomData[0] << " " << GeomData[1] << " " << GeomData[2] <<
            " N=" << GeomData[3] << " " << GeomData[4] << " " << GeomData[5] << " VMAG=" << GeomData[6] << endl;
      }
      else if (token == "POINT3") {
        assert( GeomType == -1 );
        GeomType = POINT3_GEOM;
        token = param->getString(iarg++);
        if (token != "X0")
          CERR("LSP.INJECTOR expects X0");
        GeomData[0] = param->getDouble(iarg++); // <x>
        GeomData[1] = param->getDouble(iarg++); // <y>
        GeomData[2] = param->getDouble(iarg++); // <z>
        token = param->getString(iarg++);
        if (token != "X1")
          CERR("LSP.INJECTOR expects X1");
        GeomData[3] = param->getDouble(iarg++); // <x>
        GeomData[4] = param->getDouble(iarg++); // <y>
        GeomData[5] = param->getDouble(iarg++); // <z>
        token = param->getString(iarg++);
        if (token != "X2")
          CERR("LSP.INJECTOR expects X2");
        GeomData[6] = param->getDouble(iarg++); // <x>
        GeomData[7] = param->getDouble(iarg++); // <y>
        GeomData[8] = param->getDouble(iarg++); // <z>
        token = param->getString(iarg++);
        if (token != "VMAG")
          CERR("LSP.INJECTOR expects VMAG");
        GeomData[9] = param->getDouble(iarg++);
        if (mpi_rank == 0)
          cout << "   > injector GEOM = POINT3  X0=" << GeomData[0] << " " << GeomData[1] << " " << GeomData[2] <<
            " X1=" << GeomData[3] << " " << GeomData[4] << " " << GeomData[5] <<
            " X2=" << GeomData[6] << " " << GeomData[7] << " " << GeomData[8] << 
            " VMAG=" << GeomData[9] << endl;
      }
      else if (token == "RING") {
        assert( GeomType == -1 );
        GeomType = RING_GEOM;
        token = param->getString(iarg++);
        if (token != "X")
          CERR("LSP.INJECTOR expects X");
        GeomData[0] = param->getDouble(iarg++); // <x>
        GeomData[1] = param->getDouble(iarg++); // <y>
        GeomData[2] = param->getDouble(iarg++); // <z>
        token = param->getString(iarg++);
        if (token != "N")
          CERR("LSP.INJECTOR expects N");
        GeomData[3] = param->getDouble(iarg++); // <nx>
        GeomData[4] = param->getDouble(iarg++); // <ny>
        GeomData[5] = param->getDouble(iarg++); // <nz>
        token = param->getString(iarg++);
        if (token != "R")
          CERR("LSP.INJECTOR expects R");
        GeomData[6] = param->getDouble(iarg++); // <r>
        token = param->getString(iarg++);
        if (token != "ANGLE")
          CERR("LSP.INJECTOR expects ANGLE");
        GeomData[7] = param->getDouble(iarg++); // <cone angle>
        token = param->getString(iarg++);
        if (token != "V")
          CERR("LSP.INJECTOR expects V");
        GeomData[8] = param->getDouble(iarg++); // <injection speed>
        token = param->getString(iarg++);
        if (token != "SF")
          CERR("LSP.INJECTOR expects SF");
        GeomData[9] = param->getDouble(iarg++); // <swirl fraction: 1 == 45 degrees>
        token = param->getString(iarg++);
        if (token != "SPREAD")
          CERR("LSP.INJECTOR expects SPREAD");
        GeomData[10] = param->getDouble(iarg++); // <cone angle spread fraction>
        // make sure that the normal has a unity length
        double mag = 0; FOR_I3 mag += GeomData[i+3]*GeomData[i+3];
        FOR_I3 GeomData[i+3] /= sqrt(mag);
        if (mpi_rank == 0)
          cout << "   > injector GEOM = RING X=" << GeomData[0] << " " << GeomData[1] << " " << GeomData[2] <<
            " V=" << GeomData[3] << " " << GeomData[4] << " " << GeomData[5] << " R=" << GeomData[6] <<
            " ANGLE=" << GeomData[7] << " V=" << GeomData[8] << " SF=" << GeomData[9] << " SPREAD=" << GeomData[10] << endl;
      }
      else if (token == "BOX") {
        assert( GeomType == -1 );
        GeomType = LPBOX_GEOM;
        token = param->getString(iarg++);
        if (token != "X")
          CERR("LSP.INJECTOR expects X");
        GeomData[0] = param->getDouble(iarg++); // <xmin>
        GeomData[1] = param->getDouble(iarg++); // <xmax>
        GeomData[2] = param->getDouble(iarg++); // <ymin>
        GeomData[3] = param->getDouble(iarg++); // <ymax>
        GeomData[4] = param->getDouble(iarg++); // <zmin>
        GeomData[5] = param->getDouble(iarg++); // <zmax>
        token = param->getString(iarg++);
        if (token != "V")
          CERR("LSP.INJECTOR expects V");
        GeomData[6] = param->getDouble(iarg++); // <vx>
        GeomData[7] = param->getDouble(iarg++); // <vy>
        GeomData[8] = param->getDouble(iarg++); // <vz>
        if (mpi_rank == 0)
          cout << "   > injector GEOM = BOX X=" << GeomData[0] << ":" << GeomData[1] << 
            " Y=" << GeomData[2] << ":" << GeomData[3] <<
            " Z=" << GeomData[4] << ":" << GeomData[5] << 
            " V=" << GeomData[6] << " " << GeomData[7] << " " << GeomData[8] << endl;
        if ( (GeomData[0] > GeomData[1]) or (GeomData[2] > GeomData[3]) or (GeomData[4] > GeomData[5]) ) 
          CERR("BOX is defined incorrectly in LSP.STATIC_INJECTOR. Check box dimensions...");
      } 
      else if (token == "CIRCLE") {
        assert( GeomType == -1 );
        GeomType = CIRCLE_GEOM;
        token = param->getString(iarg++);
        if (token != "X")
          CERR("LSP.INJECTOR expects X");
        GeomData[0] = param->getDouble(iarg++); // <x>
        GeomData[1] = param->getDouble(iarg++); // <y>
        GeomData[2] = param->getDouble(iarg++); // <z>
        token = param->getString(iarg++);
        if (token != "N")
          CERR("LSP.INJECTOR expects N");
        double n[3], mag_n;
        n[0] = param->getDouble(iarg++);
        n[1] = param->getDouble(iarg++);
        n[2] = param->getDouble(iarg++);
        mag_n = MAG(n);
        FOR_J3 n[j] /= mag_n;
        GeomData[3] = n[0];
        GeomData[4] = n[1];
        GeomData[5] = n[2];
        token = param->getString(iarg++);
        if (token != "R")
          CERR("LSP.INJECTOR expects R");
        GeomData[6] = param->getDouble(iarg++);
        token = param->getString(iarg++);
        if (token != "V")
          CERR("LSP.INJECTOR expects V");
        GeomData[7] = param->getDouble(iarg++);
        if (mpi_rank == 0)
          cout << "   > injector GEOM = CIRCLE X=" << GeomData[0] << " " << GeomData[1] << " " << GeomData[2] << " ,N= " << GeomData[3] << " " << GeomData[4] << " " << GeomData[5] << " ,R= " << GeomData[6] << 
          " ,V= " << GeomData[7] << endl;
      }
      else if (token == "ZONE") {
        assert( GeomType == -1 );
        GeomType = ZONE_GEOM;

        string zone_name = param->getString(iarg++);

        token = param->getString(iarg++);
        if (token != "VMAG")
          CERR("LSP.INJECTOR expects VMAG");
        GeomData[0] = param->getDouble(iarg++);

        int ierr = initZoneInjector(zone_name);
        if (ierr == 1)
          CERR("Problem in ZONE injector initiation, check the zone name")

      }
      else if (token == "PLANE") {
        assert( GeomType == -1 );
        GeomType = LPPLANE_GEOM;

        token = param->getString(iarg++);
        if (token != "X")
          CERR("LSP.INJECTOR expects X");
        GeomData[0] = param->getDouble(iarg++);
        GeomData[1] = param->getDouble(iarg++);
        GeomData[2] = param->getDouble(iarg++);
        
        token = param->getString(iarg++);
        if (token != "N")
          CERR("LSP.INJECTOR expects N");
        GeomData[3] = param->getDouble(iarg++);
        GeomData[4] = param->getDouble(iarg++);
        GeomData[5] = param->getDouble(iarg++);
        NORMALIZE((GeomData+3));

        token = param->getString(iarg++);
        if (token != "VMAG")
          CERR("LSP.INJECTOR expects VMAG");
        GeomData[6] = param->getDouble(iarg++);

        int ierr = initPlaneInjector();
        if (ierr == 1) 
          CERR("Problem in PLANE injector initiation, check the plane location and normal");

      }
      else if (token == "CONE") {
        assert( GeomType == -1 );
        GeomType = CONE_GEOM;
        
        bool got_x = false;
        bool got_n = false;
        bool got_a = false;
        bool got_v = false;
        int nParam = 4;
        try {
          for (int iter = 0; iter < nParam; iter++) {
            string token = param->getString(iarg++);
            if (token == "X") {
              got_x = true;
              GeomData[0] = param->getDouble(iarg++);
              GeomData[1] = param->getDouble(iarg++);
              GeomData[2] = param->getDouble(iarg++);
            }
            else if (token == "N") {
              got_n = true;
              GeomData[3] = param->getDouble(iarg++);
              GeomData[4] = param->getDouble(iarg++);
              GeomData[5] = param->getDouble(iarg++);
              NORMALIZE((GeomData+3));

              // form e1 and e2
              MiscUtils::getBestE1E2FromE0(e1, e2, (GeomData+3));
            }
            else if (token == "CONE_ANGLE") {
              got_a = true;
              GeomData[6] = param->getDouble(iarg++);
            }
            else if (token == "VMAG") {
              got_v = true;
              GeomData[7] = param->getDouble(iarg++);
            }
            else {
              throw token;
            }
          }
          if (not got_x) throw 1;
          if (not got_n) throw 2;
          if (not got_a) throw 3;
          if (not got_v) throw 4;
        }
        catch (int ierr) {
          if (ierr == 1) {
            CERR("couldn't find X for CONE");
          } else if (ierr == 2) {
            CERR("couldn't find N for CONE");
          } else if (ierr == 3) {
            CERR("couldn't find CONE_ANGLE for CONE");
          } else if (ierr == 4) {
            CERR("couldn't find VMAG for CONE");
          } else {
            CERR("unrecognized parameter for LSP.INJECTOR GEOM CONE");
          }
        }
        catch (string token) {
          CERR("expected parameter for CONE: X, N, CONE_ANGLE, VMAG. Make sure you have them all. reached to unrecognized token: "+token);
        }
        catch(...) {
          CERR("unknown error in GEOM=CONE");
        }

        IF_RANK0 {
          cout << "   > injector GEOM = CONE X=" << GeomData[0] << " " << GeomData[1] << " " << GeomData[2] << " ,N= " << GeomData[3] << " " << GeomData[4] << " " << GeomData[5] << " ,CONE_ANGLE= " << GeomData[6] << " ,VMAG= " << GeomData[7] << " e1: " << COUT_VEC(e1) << " e2: " << COUT_VEC(e2) << endl;
        }
      }
      else {
        if (mpi_rank == 0)
          cerr << "\n\n\n*********************************************************\n" <<
            "Error: unrecognized GEOM in LSP.INJECTOR: " << token << ". Valid geoms are: \n" <<
            "GEOM POINT X <x> <y> <z> N <nx> <ny> <nz> VMAG <vmag>\n" <<
            "   -- point source at (x,y,z), injection direction (nx,ny,nz) and velocity magnitude vmag\n" <<
            "GEOM POINT3 X0 <x> <y> <z> X1 <x> <y> <z> X2 <x> <y> <z> VMAG <vmag> \n" <<
            "   -- point source the center of X0-X1-X2 with velocity orthogonal to the tri with magnitude vmag\n" <<
            "GEOM RING <x> <y> <z> <nx> <ny> <nz> <r> <cone angle> <injection speed> <swirl ratio> <spread>\n" <<
            "   -- ring source with radius (r), center (x,y,z) and normal (nx,ny,nz), \n" <<
            "   -- injecting particels with normal distribution through a cone  \n" <<
            "GEOM BOX X <xmin> <xmax> <ymin> <ymax> <zmin> <zmax> V <vx> <vy> <vz> \n" <<
            "GEOM CIRCLE X <x> <y> <z> N <nx> <ny> <nz> R <r> V <vmag> \n" <<
            "GEOM ZONE <zone-name> VMAG <vmag> \n" <<
            "GEOM CONE X <x> <y> <z> N <nx> <ny> <nz> CONE_ANGLE <angle> VMAG <vmag> \n" <<
            "GEOM PLANE X <x> <y> <z> N <nx> <ny> <nz> VMAG <vmag> \n" <<
            "*********************************************************\n\n\n" << endl;
        throw(0);
      }
    }
    else if ( token == "TP" ) {
      // ==========================================
      // TP <temperature>
      // ==========================================
      assert( !TP_flag );
      TP_flag = true;
      Tp = param->getDouble(iarg++);
      rhop = fuel->calcRho(Tp);
      if (mpi_rank == 0)
        cout << "   > injector temperature: " << Tp << ", density: " << rhop << endl;
    }
    else if (token == "TIME_SPAN") {
      injectTime = param->getDouble(iarg++);
      if (injectTime <= 0) CERR("TIME_SPAN should be positive.");
      if (mpi_rank == 0)
        cout << "   > injector time span:  " << injectTime << endl;
    }
    else if ( token == "OVER_WRITE_VELOCITY" ) {
      overWriteVelocity_b = true;
      try {
        token = param->getString(iarg++);
      }
      catch(...) {
        CERR( "Error: OVER_WRITE_VELOCITY needs a MODE.\nOptions are: DEFAULT/MULT/MAG/MAX. Example: OVER_WRITE_VELOCITY MODE MULT 0.7");
      }
      if (token != "MODE") CERR( "Error: we need a MODE after OVER_WRITE_VELOCITY in LSP.INJECTOR");
      token = param->getString(iarg++);
      if (token == "DEFAULT") {
        overWriteVelocity_mode = DEFAULT;
      }
      else if (token == "MULT") {
        overWriteVelocity_mode = MULT;
        OWVData[0] = param->getDouble(iarg++);
      }
      else if (token == "MAG") {
        overWriteVelocity_mode = MAG;
        OWVData[0] = param->getDouble(iarg++);
        if (OWVData[0]<0) CERR( "Error: MAG value must be positive. In LSP.INJECTOR");
      }
      else if (token == "MAX") {
        overWriteVelocity_mode = MAX;
        OWVData[0] = param->getDouble(iarg++);
        if (OWVData[0]<0) CERR( "Error: MAX value must be positive. In LSP.INJECTOR");
      }
      else
        CERR( "Error: MODE: "<<token<<" after OVER_WRITE_VELOCITY was not recognized in LSP.INJECTOR");
      IF_RANK0 {
        stringstream ss0; ss0<<OWVData[0];
        cout << "   > Over writing the injection velocity to the fluid velocity with MODE: " << overWriteVelocity_mode << 
                ", values: " << ((overWriteVelocity_mode != DEFAULT) ? ss0.str():"") << " " << endl <<
                "     ("<<DEFAULT<<":DEFAULT, "<<MULT<<":MULT, "<<MAG<<":MAG, "<<MAX<<":MAX)"<<endl;
      }
    }
    else if ( token == "SHAPE") {
      token = param->getString(iarg++);
      if ( token == "SPHERE")
        pShape = SPHERE_SHAPE;
      else if ( token == "ELLIPSOID") {
        pShape = ELLIPSOID_SHAPE;
        token = param->getString(iarg++);
        if (token == "STD") {
          axis_std = param->getDouble(iarg++);
          if (axis_std < 0) CERR( "Error: axis std should be larger than 0, in LSP.INJECTOR");
        }
        else if (token == "EF") {
          axis_e = param->getDouble(iarg++);
          axis_f = param->getDouble(iarg++);
          if ((axis_e < 0)||(axis_e > 1)) CERR( "Error: elongation value should be in range 0:1, in LSP.INJECTOR");
          if ((axis_f < 0)||(axis_f > 1)) CERR( "Error: flatness value should be in range 0:1, in LSP.INJECTOR");
        }
        else if (token == "DEFAULT") 
          axis_std = 0.1;
        else 
          CERR( "Error: unrecognized token: "<<token<<" after SHAPE ELLIPSOID");
      }
      else 
        CERR( "Error: unrecognized SHAPE: " << token << " in LSP.INJECTOR. Options are: SPHERE/ELLIPSOID");
    }
    else {
      if (mpi_rank == 0)
        cerr << "\n\n\n*********************************************************\n" <<
          "Error: unrecognized keyword in LSP.INJECTOR: " << token <<
          "\n*********************************************************\n\n\n" << endl;
      throw(0);
    }
  }

  // now check...
  if ( InjectorType == -1 ) {
    if (mpi_rank == 0)
      cerr << "\n\n\n*********************************************************\n" <<
        "Error: LSP.INJECTOR " << name << ": TYPE not specified. Valid types are: \n" <<
        "TYPE COUNT <n>   -- injects n particles per time step\n" <<
        "TYPE MDOT <mdot> -- injects at rate mdot\n" <<
          "*********************************************************\n\n\n" << endl;
    throw(0);
  }

  if ( DistType == -1 ) {
    if (mpi_rank == 0)
      cerr << "\n\n\n*********************************************************\n" <<
        "Error: LSP.INJECTOR: " << name << ": DIST not specified. Valid dists are: \n" <<
        "DIST UNIFORM <diam>                                            -- injects constant diameter particles\n" <<
        "DIST RR COUNT_DIST/MASS_DIST <dmin> <dmax> <dbar> <nrr>        -- Rosin-Rammler distribution\n" <<
        "DIST GAUSSIAN COUNT_DIST/MASS_DIST <dmin> <dmax> <davg> <std>  -- Gaussian distribution\n" <<
        "*********************************************************\n\n\n" << endl;
    throw(0);
  }

  if ( GeomType == -1 ) {
    if (mpi_rank == 0)
      cerr << "\n\n\n*********************************************************\n" <<
        "Error: LSP.INJECTOR: " << name << " did not specified GEOM: Valid geoms are: \n" <<
        "GEOM POINT X <x> <y> <z> N <nx> <ny> <nz> VMAG <vmag>\n" <<
        "   -- point source at (x,y,z) with velocity (vx,vy,vz)\n" <<
        "GEOM POINT3 X0 <x> <y> <z> X1 <x> <y> <z> X2 <x> <y> <z> VMAG <vmag> \n" <<
        "   -- point source the center of X0-X1-X2 with velocity orthogonal to the tri with magnitude vmag\n" <<
        "GEOM RING <x> <y> <z> <nx> <ny> <nz> <r> <cone angle> <injection speed> <swirl ratio> \n" <<
        "   -- ring source with radius (r), center (x,y,z) and normal (nx,ny,nz), \n" <<
        "   -- injecting particels with normal distribution through a cone  \n" <<
        "GEOM BOX X <xmin> <xmax> <ymin> <ymax> <zmin> <zmax> V <vx> <vy> <vz> \n" <<
        "GEOM CIRCLE X <x> <y> <z> N <nx> <ny> <nz> R <r> V <vmag> \n" <<
        "GEOM ZONE <zone-name> VMAG <vmag> \n" <<
        "GEOM CONE X <x> <y> <z> N <nx> <ny> <nz> CONE_ANGLE <angle> VMAG <vmag> \n" <<
        "GEOM PLANE X <x> <y> <z> N <nx> <ny> <nz> VMAG <vmag> \n" <<
        "*********************************************************\n\n\n" << endl;
    throw(0);
  }

  if ( !TP_flag ) {
    if (mpi_rank == 0)
      cerr << "\n\n\n*********************************************************\n" <<
        "Error: LSP.INJECTOR: " << name << ": TP (temperature) not specified.\n" <<
        "*********************************************************\n\n\n" << endl;
    throw(0);
  }

}

void LspInjector::addNewLsp(vector<LspData>& lspDataVec,
   const double time,const double dt,const bool verbose) {

  // time: the relative time from starting the simulation
  // Don't add particles if larger than the injection time
  if ((injectTime > 0) && (time > injectTime)) return;

  // here we use the LspDataVec to add this injector's new particles

  int from_ind = lspDataVec.size();

  int npnew = 0;                         // number of particles injected at each time step
  double * mpnew = NULL;                 // particle mass
  double * Tpnew = NULL;                 // particle temperature
  double * my_mpnew = NULL;              // particle mass for each rank
  double * my_Tpnew = NULL;              // particle temperature for each rank
  double (*xpnew)[3] = NULL;             // particle position
  double (*upnew)[3] = NULL;             // particle velocity
  int * np_count_rank = NULL;            // np for each rank for plane injectors

  // ============================================
  // on rank 0 build the injected particle list:
  // This is necessary because random numbers are involved
  // in some of the particle distributions, and we want
  // all ranks to start with the same distribution.
  // ============================================

  IF_RANK0 {

    // this would cout every time step, and this routine doesn't know about step and check_intervall, so for now turned it off
    //if (verbose && (step%check_interval == 0))
    //  cout << " > injector " << getName() << ": ";

    // if the injector is MDOT modify the residual mass

    if ( InjectorType == MDOT_TYPE )
      ResidualMass += mdot*dt;

    // -----------------------------
    // distribution...
    // -----------------------------

    if ( DistType == UNIFORM_DIST ) {

      double dp   = DistData[0];
      double mp   = rhop*M_PI/6.0*pow(dp,3);

      if ( InjectorType == MDOT_TYPE ){
        double npcheck = ResidualMass/mp;
        if (npcheck > 1000000.0) {
          cerr << "Error: LspInjector " << getName() << " is trying to inject " << npcheck << " particles this time step. Check LSP.INJECT settings." << endl;
          throw(-1);
        }
        npnew = (int)npcheck;
        ResidualMass -= npnew*mp;
      }
      else if ( InjectorType == COUNT_TYPE ) {
        // for COUNT injectors, the np stores the count. Use it once, then
        // reset to zero...
        npnew = np;
        np = 0;
      }
      else{
        assert(0);
      }

      // allocate data 
      if ( npnew > 0 ) {
        mpnew    = new double[npnew];
        Tpnew    = new double[npnew];
        for (int j = 0; j < npnew; ++j) {
          mpnew[j]   = mp;
          Tpnew[j]   = Tp;
        }
      }
    }

    else if (DistType == RR_DIST ) {

      /// for RR injection, the dist_data is the
      // user-specified dmin, dmax, dbar, nrr, so add particles 
      // until the stored mass is negative...

      vector<double> dpVec; // manage memory with a vector
      bool done = (InjectorType==MDOT_TYPE) ? (ResidualMass<0) : (np<=0);
      while (!done) {

        // get the next particle diameter...
        double dp = 0;
        if (cdfType == CDF_COUNT_DIST) {
          dp = diam_rr_countCdf(DistData[0],DistData[1],DistData[2],DistData[3]);
        } else if (cdfType == CDF_MASS_DIST) {
          dp = diam_rr_massCdf(DistData[0],DistData[1],DistData[2],DistData[3]);
        } else {
          assert(0);
        }
        const double mp = rhop*M_PI/6.0*pow(dp,3);

        dpVec.push_back(dp);

        switch(InjectorType) {
        case MDOT_TYPE:
          ResidualMass -= mp;
          done = (ResidualMass<0);
          break;
        case COUNT_TYPE:
          // for COUNT injectors, we inject a count of particles...
          np--;
          done = (np<=0);
          break;
        default:
          assert(0);
        }

      }

      npnew = dpVec.size();

      if ( npnew > 0 ) {
        mpnew    = new double[npnew];
        Tpnew    = new double[npnew];
        for (int j = 0; j < npnew; ++j) {
          double dp = dpVec[j];
          double mp = rhop*M_PI/6.0*pow(dp,3);
          mpnew[j]   = mp;
          Tpnew[j]   = Tp;
        }
      }

    }
    else if (DistType == GAUSSIAN_DIST ) {

      // for Gaussian injection, the dist_data is the
      // user-specified dmin, dmax, davg, dstd, so add particles 
      // until the stored mass is smaller than one particle mass...

      vector<double> dpVec; // manage memory with a vector
      bool done = (InjectorType==MDOT_TYPE) ? (ResidualMass<0) : (np<=0);
      while (!done) {

        // get the next particle diameter...
        const double dp = diam_gaussian(DistData[0],DistData[1],DistData[2],DistData[3]);
        const double mp = rhop*M_PI/6.0*pow(dp,3);

        dpVec.push_back(dp);

        switch(InjectorType) {
        case MDOT_TYPE:
          ResidualMass -= mp;
          done = (ResidualMass<0);
          break;
        case COUNT_TYPE:
          // for COUNT injectors, we inject a count of particles...
          np--;
          done = (np<=0);
          break;
        default:
          assert(0);
        }

      }

      npnew = dpVec.size();

      if ( npnew > 0 ) {
        mpnew    = new double[npnew];
        Tpnew    = new double[npnew];
        for (int j = 0; j < npnew; ++j) {
          double dp = dpVec[j];
          double mp = rhop*M_PI/6.0*pow(dp,3);
          mpnew[j]   = mp;
          Tpnew[j]   = Tp;
        }
      }
    }

    else {
      assert(0);
    }

    // -----------------------
    // then geometry...
    // -----------------------

    if ( npnew > 0 ) {

      xpnew = new double[npnew][3];
      upnew = new double[npnew][3];

      switch (GeomType) {

      case POINT_GEOM:
      {
        for (int j = 0; j < npnew; ++j) {
          FOR_I3 xpnew[j][i] = GeomData[i];
          FOR_I3 upnew[j][i] = GeomData[i+3]*GeomData[6];
        }
        break;
      }
      case POINT3_GEOM:
        {
          double X[3][3];
          for (int ix = 0; ix < 3; ++ix) {
            FOR_I3 X[ix][i] = GeomData[ix*3 + i];
          }
          const double X0X1[3] = DIFF(X[1],X[0]);
          const double X0X2[3] = DIFF(X[2],X[0]);
          double n[3] = CROSS_PRODUCT(X0X1,X0X2);
          const double magn = MAG(n);
          FOR_I3 n[i] /= magn;
          for (int j = 0; j < npnew; ++j) {
            FOR_I3 xpnew[j][i] = (X[0][i] + X[1][i] + X[2][i])/3.;
            FOR_I3 upnew[j][i] = GeomData[9]*n[i];
            //cout << "!!! checking velocity: " << COUT_VEC(upnew[j]) << ", checking n: " << COUT_VEC(n) << " , checking MAG(n): " << MAG(n) << endl;
          }
        }
        break;

      case RING_GEOM:
        {
          // copy the data
          double XInjector[3];
          int i = 0;
          XInjector[0] = GeomData[i++];
          XInjector[1] = GeomData[i++];
          XInjector[2] = GeomData[i++];
          double Vinjector[3];
          Vinjector[0] = GeomData[i++];
          Vinjector[1] = GeomData[i++];
          Vinjector[2] = GeomData[i++];
          double RInjector        = GeomData[i++];
          double ConeAngleDegrees = GeomData[i++];
          double InjectSpeed      = GeomData[i++];
          double SwirlRatio       = GeomData[i++];
          double spreadFrac       = GeomData[i++];
          assert( i==11 );

          // generate random numbers for both location and injection angle
          double * RandX = new double[npnew];
          for ( int j =0; j < npnew; ++j) {
            RandX[j] = 2.0*M_PI*(double)(rand())/(double)(RAND_MAX); // 0 to 2pi
          }

          // Box-Muller
          double * RandA = new double[npnew];
          for ( int j =0; j < npnew; ++j) {
            // not sure id rand() can return 0, so make sure it doesn't as follows...
            double u1 = (double(rand())+1.0)/(double(RAND_MAX)+1.0);
            double u2 = (double(rand())+1.0)/(double(RAND_MAX)+1.0);
            double normal_rand_number = sqrt(-2.0*log(u1))*cos(2.0*M_PI*u2);
            // a second normal (not used here) would be 
            // = sqrt(-2.0*log(u1))*sin(2.0*M_PI*u2);
            // scatter the points about the user-defined cone half-angle according to spreadFrac
            RandA[j] = ConeAngleDegrees*M_PI/180.0*(1.0 + spreadFrac*normal_rand_number);
            // should be very rare, but limit between -90 and 90...
            RandA[j] = max(-M_PI/2.0,min(M_PI/2.0,RandA[j]));
          }

          //dumpHist0(RandA,npnew,"RandA");
          //dumpHist0(RandX,npnew,"RandX");
          //getchar();

          // build a set of orthogonal vectors
          // first vector is just the injection direction (already normalized)...
          double e1[3]; FOR_I3 e1[i] = Vinjector[i];

          // for e2, first choose the min component of e1...
          double e3[3] = { 0.0, 0.0, 0.0 };
                if      (fabs(e1[0]) <= min(fabs(e1[1]),fabs(e1[2])))
                  e3[0] = 1.0;
                else if (fabs(e1[1]) <= min(fabs(e1[0]),fabs(e1[2])))
                  e3[1] = 1.0;
                else
                  e3[2] = 1.0;

          // e2 is the cross product 
          double e2[3];
          e2[0] = e1[1]*e3[2] - e1[2]*e3[1];
          e2[1] = e1[2]*e3[0] - e1[0]*e3[2];
          e2[2] = e1[0]*e3[1] - e1[1]*e3[0];
          double inv_mag = 1.0/sqrt(DOT_PRODUCT(e2,e2));
          FOR_I3 e2[i] *= inv_mag;

          // e3 is e1 x e2... 
          e3[0] = e1[1]*e2[2] - e1[2]*e2[1];
          e3[1] = e1[2]*e2[0] - e1[0]*e2[2];
          e3[2] = e1[0]*e2[1] - e1[1]*e2[0];

          // should already be normal...
          assert( fabs(1.0-sqrt(DOT_PRODUCT(e3,e3))) < 1.0E-10);

          // break the velocity into a tangential part and a "normal" part. The normal
          // part is the projection in the "e1-R" plane...
          const double u_n = InjectSpeed/pow((1.0+pow(SwirlRatio,2)),0.5);
          const double u_t = u_n*SwirlRatio;

          // for all the new particles
          for (int j = 0; j < npnew; ++j) {

            // first randomly generate particles on a 360 ring
            double Rdirection[3]; FOR_I3 Rdirection[i] = cos(RandX[j])*e2[i] + sin(RandX[j])*e3[i];
            FOR_I3 xpnew[j][i] = XInjector[i] + RInjector*Rdirection[i];

            // swirl direction T = e1 x R...
            double Tdirection[3];
            Tdirection[0] = e1[1]*Rdirection[2] - e1[2]*Rdirection[1];
            Tdirection[1] = e1[2]*Rdirection[0] - e1[0]*Rdirection[2];
            Tdirection[2] = e1[0]*Rdirection[1] - e1[1]*Rdirection[0];

            FOR_I3 upnew[j][i] = u_t*Tdirection[i] + u_n*(cos(RandA[j])*e1[i] + sin(RandA[j])*Rdirection[i]);

          }

          delete[] RandX;
          delete[] RandA;

        }
        break;
    
        case LPBOX_GEOM:
        {
          for (int j = 0; j < npnew; ++j) {
            double rand_x[3];
            FOR_I3 rand_x[i] = rand()/double(RAND_MAX);
            FOR_I3 xpnew[j][i] = GeomData[2*i] + (GeomData[2*i+1] - GeomData[2*i])*rand_x[i];
            FOR_I3 upnew[j][i] = GeomData[i+6];
          }
        }
        break;

        case CIRCLE_GEOM:
        {
          double X[3], n[3], R;
          FOR_I3 X[i] = GeomData[i]; 
          FOR_I3 n[i] = GeomData[i+3];
          R = GeomData[6];
          int i_min = -1;
          double n_min = 1e20;
          FOR_I3 {
            if (abs(n[i]) < n_min) {
              i_min = i;
              n_min = n[i];
            }
          }
          double dir[3];
          dir[i_min] = 0.0;
          if (i_min == 0) {
            dir[1] = n[2];
            dir[2] = -n[1];
          } else if (i_min == 1) {
            dir[0] = n[2];
            dir[2] = -n[0];
          } else if (i_min == 2) {
            dir[0] = n[1];
            dir[1] = -n[0];
          } else {assert(0);}
          double mag_dir = MAG(dir);
          assert(mag_dir>0);
          FOR_I3 dir[i] /= mag_dir;
          assert(DOT_PRODUCT(dir,n) < 1e-12);
          for (int j = 0; j < npnew; ++j) {
            double rand_r_ = rand()/double(RAND_MAX);
            double rand_t = rand()/double(RAND_MAX);   // random variables for r and theta
            double r = R * sqrt(rand_r_);
            double theta = rand_t * 2. * M_PI;
            double Rot[3][3];
            // build three dimensional rotation matrix
            FOR_I3 Rot[i][i] = cos(theta) + n[i]*n[i]*(1.-cos(theta));
            Rot[0][1] = n[0]*n[1]*(1.-cos(theta)) - n[2]*sin(theta);
            Rot[1][0] = n[0]*n[1]*(1.-cos(theta)) + n[2]*sin(theta);
            Rot[0][2] = n[0]*n[2]*(1.-cos(theta)) + n[1]*sin(theta);
            Rot[2][0] = n[0]*n[2]*(1.-cos(theta)) - n[1]*sin(theta);
            Rot[2][1] = n[1]*n[2]*(1.-cos(theta)) + n[0]*sin(theta);
            Rot[1][2] = n[1]*n[2]*(1.-cos(theta)) - n[0]*sin(theta);
            double rand_dir[3] = MATRIX_PRODUCT(Rot,dir);
            assert(MAG(rand_dir) - 1. < 1e-12);
            assert(DOT_PRODUCT(rand_dir,n) < 1e-12);
            FOR_I3 xpnew[j][i] = X[i] + r * rand_dir[i];
            FOR_I3 upnew[j][i] = n[i] * GeomData[7];
          }
        }
        break;

        case ZONE_GEOM:
        {
          assert(this->nst>0);
          assert(this->nsp>0);
          assert(spost);
          assert(xp);
          assert(cdf_area);
          for (int j = 0; j < npnew; ++j) {
            // first randomly choose a tri then choose the point in the tri randomly
            double rand_val[3];
            FOR_I3 rand_val[i] = rand()/double(RAND_MAX);
            int ist_rand = 0;
            for (int ist = 1; ist < this->nst+1; ist++) {
              if (rand_val[0] < cdf_area[ist]) {
                ist_rand = ist-1;
                break;
              }
            }
            //IF_RANK0 cout << "rand_val: " << COUT_VEC(rand_val) << " ist_rand: " << ist_rand << endl;
            assert((ist_rand>=0)&&(ist_rand<this->nst));

            const double * const x0 = xp[spost[ist_rand][0]];
            const double * const x1 = xp[spost[ist_rand][1]];
            const double * const x2 = xp[spost[ist_rand][2]];

            // set the normal to be inward
            double normal[3] = TRI_NORMAL_2(x0,x2,x1);
            double mag = MAG(normal);
            FOR_I3 normal[i] /= mag;
            
            // locate a random point inside the tri by partitioning 1 into three segments
            double rand_min = min(rand_val[1],rand_val[2]);
            double rand_max = max(rand_val[1],rand_val[2]);
            // move the particle slightly inside
            double small = 1e-7;
            FOR_I3 xpnew[j][i] = rand_min*x0[i] + (rand_max-rand_min)*x1[i] + (1. - rand_max)*x2[i] + small*normal[i];
            FOR_I3 upnew[j][i] = GeomData[0]*normal[i];
          }
        }
        break;

        case LPPLANE_GEOM:
        {
          assert(cdf_area_rank);
          assert(np_count_rank==NULL);
          np_count_rank = new int [mpi_size];
          FOR_RANK np_count_rank[rank] = 0;
  
          int * np_rank = new int [npnew];
  
          for (int ip = 0; ip < npnew; ++ip) {
            double rand_rank = double(rand())/(double(RAND_MAX) + 1.0);
            FOR_RANK {
              if (rand_rank < cdf_area_rank[rank]) {
                assert(rank < mpi_size);
                np_count_rank[rank]++;
                np_rank[ip] = rank;
                break;
              }
            }
          }

          int np_sum = 0;
          FOR_RANK np_sum += np_count_rank[rank];
          assert(np_sum == npnew);

  
          // sort the mpnew and Tpnew based on their rank
          vector<pair<double, int> > mpSortVec; // first: mp, second: rank
          vector<pair<double, int> > TpSortVec; // first: Tp, second: rank
          for (int ip = 0; ip < npnew; ++ip) {
            pair<double, int> tempMp(mpnew[ip], np_rank[ip]);
            pair<double, int> tempTp(Tpnew[ip], np_rank[ip]);
            mpSortVec.push_back(tempMp);
            TpSortVec.push_back(tempTp);
          }
  
          std::sort(mpSortVec.begin(), mpSortVec.end(), CampFunction);
          std::sort(TpSortVec.begin(), TpSortVec.end(), CampFunction);
          
          for (int ip = 0; ip < npnew; ++ip) {
            mpnew[ip] = mpSortVec[ip].first;
            Tpnew[ip] = TpSortVec[ip].first;
          }

          delete[] np_rank;

          //assert(this->nst>0);
          //assert(this->nsp>0);
          //assert(spost);
          //assert(xp);
          //assert(cdf_area);
          //for (int j = 0; j < npnew; ++j) {
          //  // first randomly choose a tri then choose the point in the tri randomly
          //  double rand_val[3];
          //  FOR_I3 rand_val[i] = rand()/double(RAND_MAX);
          //  int ist_rand = 0;
          //  for (int ist = 1; ist < this->nst+1; ist++) {
          //    if (rand_val[0] < cdf_area[ist]) {
          //      ist_rand = ist-1;
          //      break;
          //    }
          //  }
          //  //IF_RANK0 cout << "rand_val: " << COUT_VEC(rand_val) << " ist_rand: " << ist_rand << endl;
          //  assert((ist_rand>=0)&&(ist_rand<this->nst));

          //  const double * const x0 = xp[spost[ist_rand][0]];
          //  const double * const x1 = xp[spost[ist_rand][1]];
          //  const double * const x2 = xp[spost[ist_rand][2]];

          //  const double * const np_plane = GeomData+3;

          //  // locate a random point inside the tri by partitioning 1 into three segments
          //  double rand_min = min(rand_val[1],rand_val[2]);
          //  double rand_max = max(rand_val[1],rand_val[2]);
          //  FOR_I3 xpnew[j][i] = rand_min*x0[i] + (rand_max-rand_min)*x1[i] + (1. - rand_max)*x2[i];
          //  FOR_I3 upnew[j][i] = GeomData[6]*np_plane[i];
          //}
        }
        break;

        case CONE_GEOM:
        {
          for (int j = 0; j < npnew; ++j) {
            double rand_val[2];
            FOR_I2 rand_val[i] = rand()/double(RAND_MAX);

            double theta = (rand_val[0] - 0.5) * 2.0 * (GeomData[6]/180.*M_PI);
            double phi = rand_val[1] * 2.0 * M_PI;
            double sin_theta = sin(theta);
            double cos_theta = cos(theta);
            double sin_phi = sin(phi);
            double cos_phi = cos(phi);
            // if we take e1, e2, and n as the main axes, 
            // we rotate in the e1-n plane with theta, then rotate in e1-e2 plane with phi
            // so n = (0,0,1) goes to...
            double normal[3] = {0,0,0};
            FOR_I3 normal[i] += -cos_phi*sin_theta*e1[i];
            FOR_I3 normal[i] += -sin_phi*sin_theta*e2[i];
            FOR_I3 normal[i] += cos_theta*GeomData[3+i];
            assert((DOT_PRODUCT(normal,normal)-1.)<1e-6);

            FOR_I3 xpnew[j][i] = GeomData[i];
            FOR_I3 upnew[j][i] = GeomData[7]*normal[i];
          }
          
        }
        break;

      default:
        assert(0);
      }

    }

    // report...

    // this would cout every time step, and this routine doesn't know about step and check_intervall, so for now turned it off
    //if (verbose && (step%check_interval == 0)) {
    //  if (npnew > 0) {
    //    double dp3_min = 1.0E+10;
    //    double dp3_max = 0.0;
    //    double mass_sum = 0.0;
    //    for (int j = 0; j < npnew; ++j) {
    //      double dp3 = mpnew[j]/(rhop*M_PI/6.0);
    //      dp3_min = min(dp3_min,dp3);
    //      dp3_max = max(dp3_max,dp3);
    //      mass_sum += mpnew[j];
    //    }
    //    cout << "injecting " << npnew <<
    //      " particles, diam range: " << pow(dp3_min,1.0/3.0) << ":" << pow(dp3_max,1.0/3.0) <<
    //      ", total mass: " << mass_sum << endl;
    //  }
    //  else {
    //    cout << " no parcels" << endl;
    //  }
    //}

  }

  // ======================================
  // rank 0 bcast to everyone...
  // every rank should deside which particles
  // should be injected in their rank.
  // for plane injector each rank should 
  // determine the location of injection as well.
  // ======================================

  MPI_Bcast(&npnew,1,MPI_INT,0,mpi_comm);
  if (npnew == 0) return;

  switch (GeomType) {

    case LPPLANE_GEOM:    // in case of a plane injector the location of the points are not still determined
                          // so each rank prepares particles for itself
    {

    //const double * const xp_plane = GeomData;
    const double * const np_plane = GeomData+3;

    assert(cdf_area_tri);
    IF_RANK0 {
      assert(np_count_rank);
      assert(cdf_area_rank);
    }

    int my_npnew;
    MPI_Scatter(np_count_rank, 1, MPI_INT, &my_npnew, 1, MPI_INT, 0, mpi_comm);
    
    //if (mpi_rank != 0) {
      my_mpnew    = new double[my_npnew];
      my_Tpnew    = new double[my_npnew];
    //}

    int * send_disp;
    IF_RANK0 {
      send_disp = new int [mpi_size];
      send_disp[0] = 0;
      for (int irank = 1; irank < mpi_size; irank++) {
        send_disp[irank] = send_disp[irank-1] + np_count_rank[irank-1];
      }
      assert((send_disp[mpi_size-1]+np_count_rank[mpi_size-1])==npnew);
    }
    MPI_Scatterv(mpnew, np_count_rank, send_disp, MPI_DOUBLE, my_mpnew, my_npnew, MPI_DOUBLE, 0, mpi_comm);
    MPI_Scatterv(Tpnew, np_count_rank, send_disp, MPI_DOUBLE, my_Tpnew, my_npnew, MPI_DOUBLE, 0, mpi_comm);
    
    IF_RANK0 delete[] send_disp;

    int npold = lspDataVec.size();
    lspDataVec.resize(npold+my_npnew);

    int ipnew = 0;
    while (ipnew < my_npnew) {
      
      int itri = -1;
      double rand_tri = double(rand()) / (double(RAND_MAX) + 1.0);
      for (int ii = 0; ii < triVec.size(); ++ii) {
        if (rand_tri < cdf_area_tri[ii]) {
          itri = ii;
          break;
        }
      }
      assert((itri>=0)&&(itri<triVec.size()));

      // locate a random point inside the tri by partitioning 1 into three segments
      double rand_val[2]; 
      FOR_I2 rand_val[i] = double(rand()) / (double(RAND_MAX) + 1.0);
      double rand_min = min(rand_val[0],rand_val[1]);
      double rand_max = max(rand_val[0],rand_val[1]);
  
      const double * const x0 = triVec[itri].first.x0;
      const double * const x1 = triVec[itri].first.x1;
      const double * const x2 = triVec[itri].first.x2;

      double xp_seed[3];
      FOR_I3 xp_seed[i] = rand_min*x0[i] + (rand_max-rand_min)*x1[i] + (1. - rand_max)*x2[i];

      int icv_closest;
      const int isInside = pointIsInside(xp_seed, icv_closest);

      if ( isInside ) {
        assert( (icv_closest >= 0) && (icv_closest < solver->ncv));
        lspDataVec[npold+ipnew].icv = icv_closest;
        lspDataVec[npold+ipnew].mp = my_mpnew[ipnew];
        lspDataVec[npold+ipnew].Tp = my_Tpnew[ipnew];
        FOR_I3 lspDataVec[npold+ipnew].xp[i] = xp_seed[i];
        FOR_I3 lspDataVec[npold+ipnew].up[i] = GeomData[6]*np_plane[i];
        if (pShape == ELLIPSOID_SHAPE) {
          if ( (axis_e>0)&&(axis_f>0) ) {
            lspDataVec[npold+ipnew].e = axis_e;
            lspDataVec[npold+ipnew].f = axis_f;
          } 
          else if (axis_std>0) {
            const double dmin = 1e-10;
            lspDataVec[npold+ipnew].e = sampleNormalDist(1, axis_std, dmin, 1);
            lspDataVec[npold+ipnew].f = sampleNormalDist(1, axis_std, dmin, 1);
          }
        }
        else {
          lspDataVec[npold+ipnew].e = 1.;
          lspDataVec[npold+ipnew].f = 1.;
        }
        ipnew++;
      }
    }

    int np_count;
    MPI_Allreduce(&ipnew, &np_count, 1, MPI_INT, MPI_SUM, mpi_comm);
    assert(np_count == npnew);

    delete[] np_count_rank;
    }
    break;

    default:      // in the default case rank 0 sets all mp, Tp, xp, up and we just decide which rank owns which particles
    {
    if (mpi_rank != 0) {
      mpnew    = new double[npnew];
      Tpnew    = new double[npnew];
      xpnew    = new double[npnew][3];
      upnew    = new double[npnew][3];
    }

    MPI_Bcast(mpnew,            npnew,MPI_DOUBLE,0,mpi_comm);
    MPI_Bcast(Tpnew,            npnew,MPI_DOUBLE,0,mpi_comm);
    MPI_Bcast((double*)xpnew, 3*npnew,MPI_DOUBLE,0,mpi_comm);
    MPI_Bcast((double*)upnew, 3*npnew,MPI_DOUBLE,0,mpi_comm);

    // add to the lspDataVec...

    int my_np_count = 0;
    for (int ip = 0; ip < npnew; ++ip) {
      const double xp_seed[3] = {xpnew[ip][0], xpnew[ip][1], xpnew[ip][2]};
      int icv_closest;
      const int isInside = pointIsInside(xp_seed,icv_closest);
      if ( isInside ) {
        assert( (icv_closest >= 0) && (icv_closest < solver->ncv));
        lspDataVec.resize(lspDataVec.size()+1);
        int iback = lspDataVec.size()-1;
        lspDataVec[iback].icv = icv_closest;
        lspDataVec[iback].mp = mpnew[ip];
        lspDataVec[iback].Tp = Tpnew[ip];
        FOR_I3 lspDataVec[iback].xp[i] = xp_seed[i];
        FOR_I3 lspDataVec[iback].up[i] = upnew[ip][i];
        if (pShape == ELLIPSOID_SHAPE) {
          if ( (axis_e>0)&&(axis_f>0) ) {
            lspDataVec[iback].e = axis_e;
            lspDataVec[iback].f = axis_f;
          } 
          else if (axis_std>0) {
            const double dmin = 1e-10;
            lspDataVec[iback].e = sampleNormalDist(1, axis_std, dmin, 1);
            lspDataVec[iback].f = sampleNormalDist(1, axis_std, dmin, 1);
          }
        }
        else {
          lspDataVec[iback].e = 1.;
          lspDataVec[iback].f = 1.;
        }
        my_np_count++;
      }
    }

    int np_count;
    MPI_Allreduce(&my_np_count, &np_count, 1, MPI_INT, MPI_SUM, mpi_comm);
    assert(np_count == npnew);

    }
    break;
  }

  if (overWriteVelocity_b) overWriteVelocity(lspDataVec,from_ind);

  delete[] mpnew;
  delete[] Tpnew;
  delete[] xpnew;
  delete[] upnew;
  delete[] my_mpnew;
  delete[] my_Tpnew;

}

double LspInjector::diam_rr_countCdf(const double &dmin,const double &dmax,const double &davg,const double &q) {

  // clipped Rosin-Rammler (Weibull) distribution
  //
  // dmin --> minimum diameter
  // dmax --> maximum diameter
  // davg --> scale parameter (not a true mean diameter)
  // q    --> shape parameter (exponent in distribution)

  assert((c0>=0)&&(c0<1));
  assert((c1>0)&&(c1<=1));
  const double r = double(rand())/double(RAND_MAX); // 0 < r < 1
  const double rand_val = c0 + r*(c1-c0);

  double d_sel = davg*pow(-log(1.-rand_val),1.0/q);

  d_sel = min(d_sel,dmax);
  d_sel = max(d_sel,dmin);

  return d_sel;
}

double LspInjector::diam_rr_massCdf(const double &dmin,const double &dmax,const double &davg,const double &q) {

  // clipped Rosin-Rammler (Weibull) distribution
  //
  // dmin --> minimum diameter
  // dmax --> maximum diameter
  // davg --> scale parameter (not a true mean diameter)
  // q    --> shape parameter (exponent in distribution)

  assert((c0>=0)&&(c0<1));
  assert((c1>0)&&(c1<=1));
  const double r = double(rand())/double(RAND_MAX); // 0 < r < 1
  const double rand_val = c0 + r*(c1-c0);

  int bin_sel = -1;

  for (int ibin = 1; ibin < countCdf_n; ibin++) {

    if (countCdf_y[ibin] >= rand_val) {
      bin_sel = ibin - 1;
      break;
    }
      
  }
  assert(bin_sel >= 0);

  // add randomness inside each bin
  const double dbin = countCdf_x[2] - countCdf_x[1];
  const double r_inbin = double(rand())/double(RAND_MAX); // 0 < r_inbin < 1
  double d_sel = countCdf_x[bin_sel] + r_inbin*dbin;
  d_sel = min(d_sel,dmax);
  d_sel = max(d_sel,dmin);

  return d_sel;
  
}

double LspInjector::diam_gaussian(const double &dmin, const double &dmax, const double &davg, const double &d_std) {
  
  assert((c0>=0)&&(c0<1));
  assert((c1>0)&&(c1<=1));
  const double r = double(rand())/double(RAND_MAX); // 0 < r < 1
  const double rand_val = c0 + (c1 - c0)*r;

  int bin_sel = -1;

  for (int ibin = 1;ibin < countCdf_n; ibin++) { 

    if (countCdf_y[ibin] >= rand_val) {
      bin_sel = ibin - 1;
      break;
    }

  }
  assert(bin_sel >= 0);

  // add randomness inside each bin
  const double dbin = countCdf_x[2] - countCdf_x[1];
  const double r_inbin = double(rand())/double(RAND_MAX); // 0 < r_inbin < 1
  double d_sel = countCdf_x[bin_sel] + r_inbin*dbin;
  d_sel = min(d_sel,dmax);
  d_sel = max(d_sel,dmin);

  return d_sel;
}

void LspInjector::build_countCdf_stuff() {
  
  IF_RANK0 cout << "   > build_countCdf_stuff()" << endl;

  assert(countCdf_x==NULL);
  assert(countCdf_y==NULL);
  assert(countCdf_n>0);

  switch (DistType) {
  case UNIFORM_DIST:
    break;
  case RR_DIST:
    {
      const double dmin = DistData[0];
      const double dmax = DistData[1];
      const double davg = DistData[2];
      const double n = DistData[3];

      switch (cdfType) {
      case CDF_COUNT_DIST:
        // count cdf is analytical and only need the beginning and the end of y axis
        c0 = 1. - exp(-pow((dmin/davg),n));
        c1 = 1. - exp(-pow((dmax/davg),n));
        break;
      case CDF_MASS_DIST:
        // count cdf should be build from analytical mass cdf
        {
        countCdf_x = new double [countCdf_n];
        countCdf_y = new double [countCdf_n];
        for (int ibin = 0; ibin < countCdf_n; ibin++) {
          countCdf_x[ibin] = 0.;
          countCdf_y[ibin] = 0.;
        }

        const double dmax_ext = dmax*5.;
        const double dmin_ext = 0.;
        const double dbin = (dmax_ext-dmin_ext)/double(countCdf_n);

        const int dmin_ind = int(floor(dmin / dbin)); 
        const int dmax_ind = int(floor(dmax / dbin)); 
        
        for (int ibin = 0; ibin < countCdf_n; ibin++) {
          countCdf_x[ibin] = dmin_ext + (double(ibin))*dbin;
          assert(countCdf_y[ibin]==0.);
          if (ibin > 0) {
            countCdf_y[ibin] = countCdf_y[ibin-1] + n*pow(countCdf_x[ibin]/davg,n)/pow(countCdf_x[ibin],4)*exp(-pow(countCdf_x[ibin]/davg,n));
          }
        }
        
        // normalize the cdf with the last element of cdf
        for (int ibin = 0; ibin < countCdf_n; ibin++)
          countCdf_y[ibin] /= countCdf_y[countCdf_n-1];

        c0 = countCdf_y[dmin_ind];
        c1 = countCdf_y[dmax_ind];
        }
        break;
      default:
        assert(0);
      }
    }
    break;
  case GAUSSIAN_DIST:
    {
      switch (cdfType) {
      case CDF_COUNT_DIST:
        {
        const double dmin = DistData[0];
        const double dmax = DistData[1];
        const double davg = DistData[2];
        const double std  = DistData[3];

        // only count cdf is supported for now
        countCdf_x = new double [countCdf_n];
        countCdf_y = new double [countCdf_n];
        for (int ibin = 0; ibin < countCdf_n; ibin++) {
          countCdf_x[ibin] = 0.;
          countCdf_y[ibin] = 0.;
        }

        const double dmax_ext = dmax*5.;
        const double dmin_ext = 0.;
        const double dbin = (dmax_ext-dmin_ext)/double(countCdf_n);

        const int dmin_ind = int(floor(dmin / dbin)); 
        const int dmax_ind = int(floor(dmax / dbin)); 
        
        for (int ibin = 0; ibin < countCdf_n; ibin++) {
          countCdf_x[ibin] = dmin_ext + (double(ibin))*dbin;
          countCdf_y[ibin] = 0.5*(1.0 + error_function((countCdf_x[ibin]-davg)/(sqrt(2)*std)));
        }
        
        // normalize the cdf with the last element of cdf
        for (int ibin = 0; ibin < countCdf_n; ibin++)
          countCdf_y[ibin] /= countCdf_y[countCdf_n-1];

        c0 = countCdf_y[dmin_ind];
        c1 = countCdf_y[dmax_ind];
        }
        break;
      case CDF_MASS_DIST:
        {
        const double dmin = DistData[0];
        const double dmax = DistData[1];
        const double davg = DistData[2];
        const double std  = DistData[3];

        // only count cdf is supported for now
        countCdf_x = new double [countCdf_n];
        countCdf_y = new double [countCdf_n];
        for (int ibin = 0; ibin < countCdf_n; ibin++) {
          countCdf_x[ibin] = 0.;
          countCdf_y[ibin] = 0.;
        }

        const double dmax_ext = dmax*5.;
        const double dmin_ext = 0.;
        const double dbin = (dmax_ext-dmin_ext)/double(countCdf_n);

        const int dmin_ind = int(floor(dmin / dbin)); 
        const int dmax_ind = int(floor(dmax / dbin)); 
        
        for (int ibin = 0; ibin < countCdf_n; ibin++) {
          countCdf_x[ibin] = dmin_ext + (double(ibin))*dbin;
          assert(countCdf_y[ibin]==0.);
          if (ibin > 0) {
            countCdf_y[ibin] = countCdf_y[ibin-1] + exp(-pow(countCdf_x[ibin]-davg,2)/(2.*pow(std,2))) / pow(countCdf_x[ibin],3);
          }
        }
        
        // normalize the cdf with the last element of cdf
        for (int ibin = 0; ibin < countCdf_n; ibin++)
          countCdf_y[ibin] /= countCdf_y[countCdf_n-1];

        c0 = countCdf_y[dmin_ind];
        c1 = countCdf_y[dmax_ind];
        }
        break;
      default:
        assert(0);
      }
    }
    break;
  default:
    assert(0);
  }

}

int LspInjector::initZoneInjector(string zone_name) {

  map<const string, int>::iterator iter = solver->bfZoneNameMap.find(zone_name);
  if (iter == solver->bfZoneNameMap.end()) {
    return 1;
  }

  int zone_ind = solver->bfZoneNameMap[zone_name];
  if (mpi_rank == 0) {
    cout << "   > injector GEOM = ZONE name: " << zone_name << " id (bfZone index): " << zone_ind << " " << endl;
  }

  int my_nst = 0; // number of triangles in this zone
  for (int ist = 0; ist < solver->subSurface->nst; ist++) {
    if (solver->subSurface->znost[ist] == zone_ind)
      my_nst++;
  }

  int (* my_spost)[3] = new int [my_nst][3];
  map<int,int> spMap; // a map from subsurface isp to the zone isp
  int my_nsp = 0; // number of distinct points
  int ist = 0;
  if (my_nst > 0) {
    for (int ist_ss = 0; ist_ss < solver->subSurface->nst; ist_ss++) {
      if (solver->subSurface->znost[ist_ss] == zone_ind) {
        FOR_I3 {
          int isp_ss = solver->subSurface->spost[ist][i];
          map<int,int>::iterator it = spMap.find(isp_ss);
          if (it == spMap.end()) { // ip_old is a new point, add it to the map
            spMap[isp_ss] = my_nsp++;
          }
          my_spost[ist][i] = spMap[isp_ss];
        }
        ist++;
      }
    }
  }
  assert(ist==my_nst);

  double (* my_xp)[3] = new double [my_nsp][3];
  for (map<int,int>::iterator it = spMap.begin(); it!=spMap.end(); it++) {
    int isp_ss = it->first;
    int isp = it->second;
    FOR_I3 my_xp[isp][i] = solver->subSurface->xp[isp_ss][i];
  }

  // communicate nst, nsp
  int * nst_ranks = new int [mpi_size]; // nst for each processor
  int * nsp_ranks = new int [mpi_size]; // nsp for each processor
  MPI_Allgather(&my_nst,1,MPI_INT,nst_ranks,1,MPI_INT,mpi_comm);
  MPI_Allgather(&my_nsp,1,MPI_INT,nsp_ranks,1,MPI_INT,mpi_comm);

  this->nst = 0;
  this->nsp = 0;
  for (int iproc = 0; iproc < mpi_size; iproc++) {
    this->nst += nst_ranks[iproc];
    this->nsp += nsp_ranks[iproc];
  }

  // communicate spost
  // correct my_spost, shift the isp to accomodate for other processors
  int shift = 0;
  for (int iproc = 0; iproc < mpi_rank; iproc++) 
    shift += nsp_ranks[iproc];
  for (int ist = 0; ist < my_nst; ist++)
    FOR_I3 my_spost[ist][i] += shift;

  int * disp = new int [mpi_size];
  disp[0] = 0;
  for (int iproc = 1; iproc < mpi_size; iproc++)
    disp[iproc] = disp[iproc-1] + nst_ranks[iproc-1]*3;
  assert((disp[mpi_size-1]+nst_ranks[mpi_size-1]*3)==(this->nst*3));

  int * rcv_count = new int [mpi_size];
  for (int iproc = 0; iproc < mpi_size; iproc++)
    rcv_count[iproc] = nst_ranks[iproc]*3;
  
  assert(spost == NULL);
  spost = new int [this->nst][3];
  MPI_Allgatherv((int*)my_spost,my_nst*3,MPI_INT,(int*)spost,rcv_count,disp,MPI_INT,mpi_comm);

  // communicate xp
  disp[0] = 0;
  for (int iproc = 1; iproc < mpi_size; iproc++)
    disp[iproc] = disp[iproc-1] + nsp_ranks[iproc-1]*3;
  assert((disp[mpi_size-1]+nsp_ranks[mpi_size-1]*3)==(this->nsp*3));

  for (int iproc = 0; iproc < mpi_size; iproc++)
    rcv_count[iproc] = nsp_ranks[iproc]*3;

  assert(xp == NULL);
  xp = new double[this->nsp][3]; 
  MPI_Allgatherv((double*)my_xp,my_nsp*3,MPI_DOUBLE,(double*)xp,rcv_count,disp,MPI_DOUBLE,mpi_comm);

  assert(cdf_area == NULL);
  cdf_area = new double [this->nst+1];
  cdf_area[0] = 0.;
  double avg_normal[3] = {0,0,0};
  for (int ist = 0; ist < this->nst; ist++) {
    const double * const x0 = xp[spost[ist][0]];
    const double * const x1 = xp[spost[ist][1]];
    const double * const x2 = xp[spost[ist][2]];
    
    // twice the tri normal
    double normal[3] = TRI_NORMAL_2(x0,x2,x1);
    double area = MAG(normal) / 2.;
    cdf_area[ist+1] = cdf_area[ist] + area;

    FOR_I3 avg_normal[i] += normal[i];
  }

  FOR_I3 avg_normal[i] /= (2.*cdf_area[this->nst]);

  for (int ist = 0; ist < this->nst+1; ist++)
    cdf_area[ist] /= cdf_area[this->nst];

  IF_RANK0 {
    cout << "   > zone info, nst: " << this->nst << " nsp: " << this->nsp << " average normal: " << COUT_VEC(avg_normal) << " inject velocity: " << GeomData[0] << endl;
  //  // dump as tecplot...
  //  FILE * fp = fopen("tris.dat","w");
  //  assert(fp != NULL);
  //  fprintf(fp,"TITLE = \"tris\"\n");
  //  fprintf(fp,"VARIABLES = \"X\"\n");
  //  fprintf(fp,"\"Y\"\n");
  //  fprintf(fp,"\"Z\"\n");
  //  // zone header
  //  fprintf(fp,"ZONE T=\"blah\"\n");
  //  fprintf(fp,"N=%d, E=%d, F=FEPOINT, ET=TRIANGLE\n",this->nsp,this->nst);
  //
  //  // order should be good...
  //  for (int isp = 0; isp < this->nsp; ++isp) {
  //    fprintf(fp,"%lf %lf %lf\n",xp[isp][0],xp[isp][1],xp[isp][2]);
  //  }
  //
  //  for (int ist = 0; ist < this->nst; ++ist) {
  //    fprintf(fp,"%d %d %d\n",spost[ist][0]+1,spost[ist][1]+1,spost[ist][2]+1);
  //  }
  //
  //  fclose(fp);

  //  for (int ist = 0; ist < this->nst; ist++) 
  //    cout << "ist: " << ist << " isp: " << spost[ist][0] << " " << spost[ist][1] << " " <<spost[ist][2] << " " << COUT_VEC(xp[spost[ist][0]]) << " " << COUT_VEC(xp[spost[ist][1]]) << " " << COUT_VEC(xp[spost[ist][2]]) << endl;
  //  for (int ist = 0; ist < this->nst+1; ist++) 
  //    cout << "ist: " << ist << " cdf_area: " << cdf_area[ist] << endl;
  }
  
  delete[] my_spost;
  delete[] my_xp;
  delete[] nst_ranks;
  delete[] nsp_ranks;
  delete[] rcv_count;
  delete[] disp;
  return 0;

}

int LspInjector::initPlaneInjector() {

  const double * const xp_plane = GeomData;
  const double * const np_plane = GeomData+3;

  MiscUtils::getBestE1E2FromE0(e1, e2, np_plane);

  IF_RANK0 cout << "   > injector GEOM = PLANE xp: " << COUT_VEC(xp_plane) << " np: " << COUT_VEC(np_plane) << endl;
 
  double * iso_var_no = new double[solver->nno];

  for (int ino = 0; ino < solver->nno; ++ino) {
    iso_var_no[ino] =
      (solver->x_no[ino][0]-xp_plane[0])*np_plane[0] +
      (solver->x_no[ino][1]-xp_plane[1])*np_plane[1] +
      (solver->x_no[ino][2]-xp_plane[2])*np_plane[2];
  }
  
  assert(triVec.empty()); 
  solver->buildIsoWithIcv(triVec,iso_var_no,0.0);

  delete[] iso_var_no;

  // check if any triangles are selected
  int my_ntri = triVec.size();
  int ntri_global;
  MPI_Allreduce(&my_ntri, &ntri_global, 1, MPI_INT, MPI_SUM, mpi_comm);
  if (ntri_global == 0) return 1;

  //GeomUtils::writeTecplot("triVec.dat",triVec);
  //writeTecplot("triVec.dat",triVec);

  // put the selected cvs in groups based on their participation in trivec...

  // start by putting every cv in its own group. Recall the cv associated with
  // every tri in the triVec is stored in the second...
  int * cv_flag = new int[solver->ncv_g];
  for (int icv = 0; icv < solver->ncv; ++icv) cv_flag[icv] = -1;
  int ngr = 0;
  for (int ii = 0; ii < triVec.size(); ++ii) {
    const int icv = triVec[ii].second;
    assert((icv >= 0)&&(icv < solver->ncv));
    if (cv_flag[icv] == -1) cv_flag[icv] = ngr++;
  }

  // offset the cv_flag so we are all globally unique...
  int offset;
  MPI_Scan(&ngr,&offset,1,MPI_INT,MPI_SUM,mpi_comm);
  offset -= ngr; assert(offset >= 0);

  int icv_closest = -1;
  double d2_closest = HUGE_VAL;
  for (int icv = 0; icv < solver->ncv; ++icv) {
    if (cv_flag[icv] >= 0) {
      cv_flag[icv] += offset;
      const double d2 = DIST2(solver->x_vv[icv],xp_plane);
      if (d2 < d2_closest) {
        icv_closest = icv;
        d2_closest = d2;
      }
    }
  }
  
  int done = 0;
  while (done == 0) {
    
    int my_done = 1;
    
    solver->updateCvData(cv_flag);
    
    // TODO: this makes no use of ordering and internal-only, so 100's of iterations
    // may be necessary for large problems...
    // e.g. loop through interproc faces and see if any local cvs need to be updated,etc...

    for (int icv = 0; icv < solver->ncv; ++icv) {
      if (cv_flag[icv] >= 0) {
        for (int coc = solver->cvocv_i[icv]+1; coc != solver->cvocv_i[icv+1]; ++coc) {
          const int icv_nbr = solver->cvocv_v[coc];
          if ((cv_flag[icv_nbr] >= 0)&&(cv_flag[icv_nbr] < cv_flag[icv])) {
            cv_flag[icv] = cv_flag[icv_nbr];
            my_done = 0;
          }
        }
      }
    }
    
    MPI_Allreduce(&my_done,&done,1,MPI_INT,MPI_MIN,mpi_comm);
    
  }
  
  // at this point the group indexing will be unique, but not continuous.
  // this does not matter. We just want to save the group of the 
  // closest icv...
  
  DoubleInt my_di, di;
  my_di.this_double = d2_closest;
  my_di.this_int    = mpi_rank;
  MPI_Allreduce(&my_di,&di,1,MPI_DOUBLE_INT,MPI_MINLOC, mpi_comm);

  int igr_closest;
  if (di.this_int == mpi_rank) {
    assert(icv_closest >= 0);
    igr_closest = cv_flag[icv_closest];
    assert(igr_closest >= 0);
  }
  MPI_Bcast(&igr_closest,1,MPI_INT,di.this_int,mpi_comm);
          
  // clean the triVec...
  int ii_new = 0;
  for (int ii = 0; ii < triVec.size(); ++ii) {
    const int icv = triVec[ii].second;
    if (cv_flag[icv] == igr_closest) {
      if (ii_new != ii) {
        FOR_I3 triVec[ii_new].first.x0[i] = triVec[ii].first.x0[i];
        FOR_I3 triVec[ii_new].first.x1[i] = triVec[ii].first.x1[i];
        FOR_I3 triVec[ii_new].first.x2[i] = triVec[ii].first.x2[i];
        triVec[ii_new].second = icv;
      }
      ++ii_new;
    }
  }
  triVec.resize(ii_new);

  delete[] cv_flag;

  //GeomUtils::writeTecplot("triVec2.dat",triVec);
  //writeTecplot("triVec3.dat",triVec);
  
  // =====================================
  // also build a cdf of tris and ranks...
  // =====================================
  
  const double np_mag = MAG(np_plane); assert(np_mag > 0.0); assert(fabs(np_mag-1.0)<1e-6);
  assert(cdf_area_tri==NULL);
  cdf_area_tri = new double [triVec.size()];
  double my_area = 0.0;
  for (int ii = 0; ii < triVec.size(); ++ii) {
    //const int icv = triVec[ii].second;
    const double this_n[3] = TRI_NORMAL_2(triVec[ii].first.x0,triVec[ii].first.x1,triVec[ii].first.x2);
    const double area = 0.5*DOT_PRODUCT(this_n,np_plane);
    my_area += area;
    cdf_area_tri[ii] = fabs(area);
    if (ii > 0) cdf_area_tri[ii] += cdf_area_tri[ii-1];
  }

  for (int ii = 0; ii < triVec.size(); ++ii) {
    cdf_area_tri[ii] /= cdf_area_tri[triVec.size()-1];
  }

  assert(cdf_area_rank == NULL);
  if (mpi_rank == 0) {
    cdf_area_rank = new double[mpi_size];
  }
  MPI_Gather(&my_area,1,MPI_DOUBLE,cdf_area_rank,1,MPI_DOUBLE,0,mpi_comm);
  if (mpi_rank == 0) {
    double area_sum = 0.0;
    FOR_RANK {
      area_sum += cdf_area_rank[rank];
      cdf_area_rank[rank] = area_sum;
    }
    assert(area_sum>0);
    FOR_RANK cdf_area_rank[rank] /= area_sum;
    cout << "   > plane area: " << area_sum << endl;
  }

  return 0;

  //if (mpi_rank==12) {
  //  for (int ii = 0; ii < triVec.size(); ++ii) {
  //    if (ii > 0) assert(cdf_area_tri[ii]>=cdf_area_tri[ii-1]);
  //    cout << ii << " cdf_area_tri: " << cdf_area_tri[ii] << endl;
  //  }
  //}

  //IF_RANK0 {
  //  FOR_RANK {
  //    cout << rank << " cdf_area_rank: " << cdf_area_rank[rank] << endl;
  //  }
  //}
  //MPI_Pause("HERE");
  
}

StaticLspInjector::StaticLspInjector() : LspInjector() {
  IF_RANK0 cout << "StaticLspInjector()" << endl;
  
  npx = npy = npz = -1;

  total_mass = -1;

  np_in_cell = -1;

}

StaticLspInjector::~StaticLspInjector() {
  IF_RANK0 cout << "~StaticLspInjector()" << endl;
}

void StaticLspInjector::initFromParams(Param * param, CtiMaterial * fuel, StaticSolver * solver, double (*u) [3]) {

  this->solver = solver;
  this->u = u;

  // add the injector using the Current Parameter...

  assert(param->getName() == "LSP.STATIC_INJECTOR" );

  pShape = SPHERE_SHAPE; // shpere is default, unless otherwise is specified in SHAPE

  // expect a name. This is not used right now, but could be used to 
  // uniquely define this injector and add a tag to a particles from this injector 
  // for analysis/vis.

  try {
    
    name = param->getString();
    
    if (mpi_rank == 0)
      cout << " > initializing LSP.STATIC_INJECTOR: " << name << endl;

    bool TP_flag = false;

    int iarg = 1;
    while (iarg < param->size()) {
      string token = param->getString(iarg++);
      if (token == "TYPE") {
        // ==========================================
        // TYPE of injector:
        // ==========================================
        token = param->getString(iarg++);
        if (token == "COUNT") {
          assert( InjectorType == -1 );
          InjectorType = COUNT_TYPE;
          np = param->getInt(iarg++);
          if (mpi_rank == 0)
            cout << "   > static injector TYPE = COUNT " << np << endl;
        }
        else if ( token == "MASS" ) {
          assert( InjectorType == -1 );
          InjectorType = MASS_TYPE;
          total_mass = param->getDouble(iarg++);
          if (mpi_rank == 0)
            cout << "   > static injector TYPE = MASS, total mass: " << total_mass << endl;
        }
        else if ( token == "STRUCTURED" ) {
          assert( InjectorType == -1 );
          InjectorType = STRUCTURED_TYPE;
          npx = param->getInt(iarg++);
          npy = param->getInt(iarg++);
          npz = param->getInt(iarg++);
          np = npx * npy * npz;
          if (mpi_rank == 0)
            cout << "   > static injector TYPE = STRUCTURED npx: " << npx << " npy: " << npy << " npz: " << npz << endl;
        }
        else if (token == "EACH_CELL") {
          assert( InjectorType == -1 );
          InjectorType = EACH_CELL_TYPE;
          assert(param->getString(iarg) == "N");
          ++iarg;
          np_in_cell = param->getInt(iarg++); // number of particles in each cell
          int np_local = solver->ncv * np_in_cell;
          MPI_Reduce(&np_local,&np,1,MPI_INT,MPI_SUM,0,mpi_comm);
          if (mpi_rank == 0) {
            cout << "   > static injector TYPE = EACH_CELL N: " << np_in_cell << endl; 
            cout << "   > this injector injects " << np_in_cell << " particles in each cell with the specified size distribution" << endl;
            cout << "   > the actuall number of particles depends on the GEOM parameter. Maximum possible number of particles: " << np << endl; 
          }
        }
        else {
          CERR("unrecognized TYPE in LSP.STATIC_INJECTOR: " << token << ". Valid types are: \n" <<
               "TYPE COUNT <n>                     -- injects n particles once\n" <<
               "TYPE MASS <mass>                   -- injects total particles of specified mass\n" <<
               "TYPE STRUCTURED <npx> <npy> <npz>  -- injects npx x npy x npz equidistance particles\n" <<
               "TYPE EACH_CELL N <np in cell>      -- injects np particles in each cell");
        }
      }
      else if (token == "DIST") {
        // ==========================================
        // DIST of particles:
        // ==========================================
        token = param->getString(iarg++);
        if (token == "UNIFORM") {
          assert( DistType == -1 );
          DistType = UNIFORM_DIST;
          DistData[0] = param->getDouble(iarg++); // <diam>
          if (mpi_rank == 0)
            cout << "   > static injector DIST = UNIFORM " << " Particle diam: " << DistData[0] << endl;
        }
        else if (token == "RR") {
          assert( DistType == -1 );
          DistType = RR_DIST;
          string token = param->getString(iarg++);
          if (token=="COUNT_DIST") {
            cdfType = CDF_COUNT_DIST;
          }
          else if (token=="MASS_DIST") {
            cdfType = CDF_MASS_DIST;
          } else {
            CERR("RR cdf distribution type is missing. Options: COUNT_DIST, MASS_DIST");
          }
          DistData[0] = param->getDouble(iarg++); // <dmin>
          DistData[1] = param->getDouble(iarg++); // <dmax>
          DistData[2] = param->getDouble(iarg++); // <dbar>
          if (DistData[1] <= DistData[0]) CERR("Wrong diameter values in LSP.INJECTOR->DIST");
          DistData[3] = param->getDouble(iarg++); // <nrr>
          if (mpi_rank == 0)
            cout << "   > static injector DIST = RR, type: " << cdfType << 
              " (count_dist: " << CDF_COUNT_DIST << ", mass_dist: " << CDF_MASS_DIST << ")" << 
              " dmin: " << DistData[0] << " dmax: " << DistData[1] << " davg: " <<  DistData[2] <<
              " n: " << DistData[3] << endl;
        }
        else if (token == "GAUSSIAN") {
          assert( DistType == -1);
          DistType = GAUSSIAN_DIST;
          string token = param->getString(iarg++);
          if (token=="COUNT_DIST") {
            cdfType = CDF_COUNT_DIST;
          }
          else if (token=="MASS_DIST") {
            cdfType = CDF_MASS_DIST;
          } else {
            CERR("Gaussian cdf distribution type is missing. Options: COUNT_DIST, MASS_DIST");
          }
          DistData[0] = param->getDouble(iarg++); // min diameter
          DistData[1] = param->getDouble(iarg++); // max diameter
          DistData[2] = param->getDouble(iarg++); // mean diameter
          DistData[3] = param->getDouble(iarg++); // standard deviation
          if (DistData[1] <= DistData[0]) CERR("Wrong diameter values in LSP.INJECTOR->DIST");
          if (mpi_rank == 0)
            cout << "   > static injector DIST = GAUSSIAN type: " << cdfType << 
              " (count_dist: " << CDF_COUNT_DIST << ", mass_dist: " << CDF_MASS_DIST << ")" << 
              " dmin: " << DistData[0] << " dmax: " << DistData[1] << " davg: " << DistData[2] <<
              " dstd: " << DistData[3] << endl;
        }
        else {
          CERR("unrecognized DIST in LSP.STATIC_INJECTOR: " << token << ". Valid dists are: \n" <<
               "DIST UNIFORM <diam>                                           -- injects constant diameter particles\n" <<
               "DIST RR COUNT_DIST/MASS_DIST <dmin> <dmax> <dbar> <nrr>       -- Rosin-Rammler distribution\n" <<
               "DIST GAUSSIAN COUNT_DIST/MASS_DIST<dmin> <dmax> <davg> <std>  -- Gaussian distribution");
        }
        build_countCdf_stuff();
      }
      else if (token == "GEOM") {
        // ==========================================
        // GEOM of particles:
        // ==========================================
        token = param->getString(iarg++);
        if (token == "BOX") {
          assert( GeomType == -1 );
          GeomType = LPBOX_GEOM;
          assert(param->getString(iarg) == "X");
          ++iarg;
          GeomData[0] = param->getDouble(iarg++); // <xmin>
          GeomData[1] = param->getDouble(iarg++); // <xmax>
          GeomData[2] = param->getDouble(iarg++); // <ymin>
          GeomData[3] = param->getDouble(iarg++); // <ymax>
          GeomData[4] = param->getDouble(iarg++); // <zmin>
          GeomData[5] = param->getDouble(iarg++); // <zmax>
          assert(param->getString(iarg) == "V");
          ++iarg;
          GeomData[6] = param->getDouble(iarg++); // <vx>
          GeomData[7] = param->getDouble(iarg++); // <vy>
          GeomData[8] = param->getDouble(iarg++); // <vz>
          if (mpi_rank == 0)
            cout << "   > injector GEOM = BOX X=" << GeomData[0] << ":" << GeomData[1] << 
              " Y=" << GeomData[2] << ":" << GeomData[3] <<
              " Z=" << GeomData[4] << ":" << GeomData[5] << 
              " V=" << GeomData[6] << " " << GeomData[7] << " " << GeomData[8] << endl;
          if ( (GeomData[0] > GeomData[1]) or (GeomData[2] > GeomData[3]) or (GeomData[4] > GeomData[5]) ) 
            CERR("BOX is defined incorrectly in LSP.STATIC_INJECTOR. Check box dimensions...");
        } 
        else {
          CERR("unrecognized GEOM in LSP.STATIC_INJECTOR: " << token << ". Valid geoms are: \n" <<
               "GEOM BOX X <xmin> <xmax> <ymin> <ymax> <zmin> <zmax> V <vx> <vy> <vz>");
        }
      }
      else if ( token == "TP" ) {
        // ==========================================
        // TP <temperature>
        // ==========================================
        assert( !TP_flag );
        TP_flag = true;
        Tp = param->getDouble(iarg++);
        // use LiquidProps to get the density at this temperature...
        rhop = fuel->calcRho(Tp);
        if (mpi_rank == 0)
          cout << "   > injector temperature: " << Tp << ", density: " << rhop << endl;
      }
      else if ( token == "OVER_WRITE_VELOCITY" ) {
        overWriteVelocity_b = true;
        try {
          token = param->getString(iarg++);
        }
        catch(...) {
          CERR( "Error: OVER_WRITE_VELOCITY needs a MODE.\nOptions are: DEFAULT/MULT/MAG/MAX. Example: OVER_WRITE_VELOCITY MODE MULT 0.7");
        }
        if (token != "MODE") CERR( "Error: we need a MODE after OVER_WRITE_VELOCITY in LSP.INJECTOR");
        token = param->getString(iarg++);
        if (token == "DEFAULT") {
          overWriteVelocity_mode = DEFAULT;
        }
        else if (token == "MULT") {
          overWriteVelocity_mode = MULT;
          OWVData[0] = param->getDouble(iarg++);
        }
        else if (token == "MAG") {
          overWriteVelocity_mode = MAG;
          OWVData[0] = param->getDouble(iarg++);
          if (OWVData[0]<0) CERR( "Error: MAG value must be positive. In LSP.INJECTOR");
        }
        else if (token == "MAX") {
          overWriteVelocity_mode = MAX;
          OWVData[0] = param->getDouble(iarg++);
          if (OWVData[0]<0) CERR( "Error: MAX value must be positive. In LSP.INJECTOR");
        }
        else
          CERR( "Error: MODE: "<<token<<" after OVER_WRITE_VELOCITY was not recognized in LSP.INJECTOR");
        IF_RANK0 {
          stringstream ss0; ss0<<OWVData[0];
          cout << "   > Over writing the injection velocity to the fluid velocity with MODE: " << overWriteVelocity_mode << 
                  ", values: " << ((overWriteVelocity_mode != DEFAULT) ? ss0.str():"") << " " << endl <<
                  "     ("<<DEFAULT<<":DEFAULT, "<<MULT<<":MULT, "<<MAG<<":MAG, "<<MAX<<":MAX)"<<endl;
        }
      }
      else if ( token == "SHAPE") {
        token = param->getString(iarg++);
        if ( token == "SPHERE")
          pShape = SPHERE_SHAPE;
        else if ( token == "ELLIPSOID") {
          pShape = ELLIPSOID_SHAPE;
          token = param->getString(iarg++);
          if (token == "STD") {
            axis_std = param->getDouble(iarg++);
            if (axis_std < 0) CERR( "Error: axis std should be larger than 0, in LSP.INJECTOR");
          }
          else if (token == "EF") {
            axis_e = param->getDouble(iarg++);
            axis_f = param->getDouble(iarg++);
            if ((axis_e < 0)||(axis_e > 1)) CERR( "Error: elongation value should be in range 0:1, in LSP.INJECTOR");
            if ((axis_f < 0)||(axis_f > 1)) CERR( "Error: flatness value should be in range 0:1, in LSP.INJECTOR");
          }
          else if (token == "DEFAULT") 
            axis_std = 0.1;
          else 
            CERR( "Error: unrecognized token: "<<token<<" after SHAPE ELLIPSOID");
        }
        else 
          CERR( "Error: unrecognized SHAPE: " << token << " in LSP.INJECTOR. Options are: SPHERE/ELLIPSOID");
      }
      else {
        CERR("unrecognized keyword in LSP.STATIC_INJECTOR: " << token);
      }
    }
    
    // now check...
    if ( InjectorType == -1 ) {
      CERR("LSP.STATIC_INJECTOR " << name << ": TYPE not specified. Valid types are: \n" <<
           "TYPE COUNT <n>                             -- injects n particles per time step\n" <<
           "TYPE MASS <mass>                           -- injects total particles of specified mass\n" <<
           "TYPE STRUCTURED <npx> <npy> <npz>          -- injects npx x npy x npz equidistance particles\n" <<
           "TYPE EACH_CELL N <np in cell>              -- injects np particles in each cell");
    }
    
    if ( DistType == -1 ) {
      CERR("LSP.STATIC_INJECTOR: " << name << ": DIST not specified. Valid dists are: \n" <<
           "DIST UNIFORM <diam>                                            -- injects constant diameter particles\n" <<
           "DIST RR COUNT_DIST/MASS_DIST <dmin> <dmax> <dbar> <nrr>        -- Rosin-Rammler distribution\n" <<
           "DIST GAUSSIAN COUNT_DIST/MASS_DIST <dmin> <dmax> <davg> <std>  -- Gaussian distribution"); 
    }
    
    if ( GeomType == -1 ) {
      CERR("LSP.STATIC_INJECTOR: " << name << " did not specified GEOM: Valid geoms are: \n" <<
           "GEOM BOX X <xmin> <xmax> <ymin> <ymax> <zmin> <zmax> V <vx> <vy> <vz>");
    }
    
    if ( !TP_flag ) {
      CERR("LSP.STATIC_INJECTOR: " << name << ": TP (temperature) not specified.");
    }
    
  }
  catch (int ierr) {
    CERR("problem parsing LSP.STATIC_INJECTOR. Examples:\n" <<
         "LSP.STATIC_INJECTOR my_inj_1 GEOM=BOX X <x0> <x1> <y0> <y1> <z0> <z1> V <u> <v> <w> TYPE=COUNT <n> DIST=UNIFORM <diam> TP <temp>");
  }

}

void StaticLspInjector::addNewLsp(vector<LspData>& lspDataVec,
   const double time,const double dt,const bool verbose) {

  assert(solver->cvAdt);

  int from_ind = lspDataVec.size();

  int npnew = 0;                         // number of particles injected at each time step
  double * mpnew = NULL;                 // particle mass
  double * Tpnew = NULL;                 // particle temperature
  //double (*xpnew)[3] = NULL;             // particle position
  //double (*upnew)[3] = NULL;             // particle velocity

  // ============================================
  // on rank 0 build the injected particle list:
  // This is necessary because random numbers are involved
  // in some of the particle distributions, and we want
  // all ranks to start with the same distribution.
  // ============================================

  IF_RANK0 {

    if ( DistType == UNIFORM_DIST ) {

      double dp   = DistData[0];
      double mp   = rhop*M_PI/6.0*pow(dp,3);

      if ( InjectorType == COUNT_TYPE ) {
        npnew = np;
      } 
      else if ( InjectorType == MASS_TYPE ) {
        npnew = int(total_mass / mp);
        if (npnew > 1e7) {
          CWARN("Injecting more than 1e7 particles for LSP.STATIC_INJECTOR name: " << name << " np: " << npnew);
        }
      }
      else if ( InjectorType == STRUCTURED_TYPE ) {
        npnew = np;
      }
      else if ( InjectorType == EACH_CELL_TYPE ) {
        npnew = np;
      }
      else {
        assert(0);
      }

      // allocate data 
      mpnew    = new double[npnew];
      Tpnew    = new double[npnew];
      for (int j = 0; j < npnew; ++j) {
        mpnew[j]   = mp;
        Tpnew[j]   = Tp;
      }

    } else if (DistType == RR_DIST ) {

      /// for RR injection, the dist_data is the
      // user-specified dmin, dmax, dbar, nrr

      vector<double> dpVec; 
      bool done = false;
      while (!done) {

        // get the next particle diameter...
        double dp = 0;
        if (cdfType == CDF_COUNT_DIST) {
          dp = diam_rr_countCdf(DistData[0],DistData[1],DistData[2],DistData[3]);
        } else if (cdfType == CDF_MASS_DIST) {
          dp = diam_rr_massCdf(DistData[0],DistData[1],DistData[2],DistData[3]);
        } else {
          assert(0);
        }
        const double mp = rhop*M_PI/6.0*pow(dp,3);

        dpVec.push_back(dp);

        switch(InjectorType) {
        case COUNT_TYPE:
          np--;
          done = (np<=0);
          break;
        case MASS_TYPE:
          total_mass -= mp;
          done = (total_mass<=0);
          break;
        case STRUCTURED_TYPE:
          np--;
          done = (np<=0);
          break;
        case EACH_CELL_TYPE:
          np--;
          done = (np<=0);
          break;
        default:
          assert(0);
        }

      }

      npnew = dpVec.size();

      mpnew    = new double[npnew];
      Tpnew    = new double[npnew];
      for (int j = 0; j < npnew; ++j) {
        double dp = dpVec[j];
        double mp = rhop*M_PI/6.0*pow(dp,3);
        mpnew[j]   = mp;
        Tpnew[j]   = Tp;
      }

    } else if (DistType == GAUSSIAN_DIST ) {

      // for Gaussian injection, the dist_data is the
      // user-specified dmin, dmax, davg, dstd

      vector<double> dpVec; // manage memory with a vector
      bool done = false;
      while (!done) {

        // get the next particle diameter...
        const double dp = diam_gaussian(DistData[0],DistData[1],DistData[2],DistData[3]);
        const double mp = rhop*M_PI/6.0*pow(dp,3);

        dpVec.push_back(dp);

        switch(InjectorType) {
        case COUNT_TYPE:
          np--;
          done = (np<=0);
          break;
        case MASS_TYPE:
          total_mass -= mp;
          done = (total_mass<=0);
          break;
        case STRUCTURED_TYPE:
          np--;
          done = (np<=0);
          break;
        case EACH_CELL_TYPE:
          np--;
          done = (np<=0);
          break;
        default:
          assert(0);
        }

      }

      npnew = dpVec.size();

      mpnew    = new double[npnew];
      Tpnew    = new double[npnew];
      for (int j = 0; j < npnew; ++j) {
        double dp = dpVec[j];
        double mp = rhop*M_PI/6.0*pow(dp,3);
        mpnew[j]   = mp;
        Tpnew[j]   = Tp;
      }

    }

  }

  // ======================================
  // rank 0 bcast to everyone...
  // ======================================

  MPI_Bcast(&npnew,1,MPI_INT,0,mpi_comm);
  if (npnew == 0) return;

  if (mpi_rank != 0) {
    mpnew    = new double[npnew];
    Tpnew    = new double[npnew];
  }

  MPI_Bcast(mpnew,npnew,MPI_DOUBLE,0,mpi_comm);
  MPI_Bcast(Tpnew,npnew,MPI_DOUBLE,0,mpi_comm);

  //////////////////////////////////////////////////
  // Now set the location and velocity of particles
  //////////////////////////////////////////////////

  int my_np_count = 0;

  if (GeomType == LPBOX_GEOM) {

    double x_min = GeomData[0]; // <xmin>
    double x_max = GeomData[1]; // <xmax>
    double y_min = GeomData[2]; // <ymin>
    double y_max = GeomData[3]; // <ymax>
    double z_min = GeomData[4]; // <zmin>
    double z_max = GeomData[5]; // <zmax>

    if ( InjectorType == STRUCTURED_TYPE ) {

      // the position for this injector type is not random

      int np_total = npx * npy * npz;
      double dx = (x_max - x_min)/double(npx);
      double dy = (y_max - y_min)/double(npy);
      double dz = (z_max - z_min)/double(npz);

      for (int ip = 0; ip < np_total; ip++) {
        // set xp_seed
        double xp_seed[3];
        int iz = ip / (npx*npy);
        int iy = (ip-npx*npy*iz) / npx;
        int ix = ip % npx;
        xp_seed[0] = x_min + (double(ix)+0.5)*dx;
        xp_seed[1] = y_min + (double(iy)+0.5)*dy;
        xp_seed[2] = z_min + (double(iz)+0.5)*dz;
        assert(xp_seed[0]>=x_min);
        assert(xp_seed[0]<=x_max);
        assert(xp_seed[1]>=y_min);
        assert(xp_seed[1]<=y_max);
        assert(xp_seed[2]>=z_min);
        assert(xp_seed[2]<=z_max);

        int icv_closest;
        const int isInside = pointIsInside(xp_seed, icv_closest);
        if ( isInside ) {
          assert( (icv_closest >= 0) && (icv_closest < solver->ncv));
          lspDataVec.resize(lspDataVec.size()+1);
          int iback = lspDataVec.size()-1;
          lspDataVec[iback].icv = icv_closest;
          lspDataVec[iback].mp = mpnew[ip];
          lspDataVec[iback].Tp = Tpnew[ip];
          FOR_I3 lspDataVec[iback].xp[i] = xp_seed[i];
          FOR_I3 lspDataVec[iback].up[i] = GeomData[6+i];
          if (pShape == ELLIPSOID_SHAPE) {
            if ( (axis_e>0)&&(axis_f>0) ) {
              lspDataVec[iback].e = axis_e;
              lspDataVec[iback].f = axis_f;
            } 
            else if (axis_std>0) {
              const double dmin = 1e-10;
              lspDataVec[iback].e = sampleNormalDist(1, axis_std, dmin, 1);
              lspDataVec[iback].f = sampleNormalDist(1, axis_std, dmin, 1);
            }
          }
          else {
            lspDataVec[iback].e = 1.;
            lspDataVec[iback].f = 1.;
          }
          my_np_count++;
        }

      }

      // we are finished
      int count;
      MPI_Reduce(&my_np_count,&count,1,MPI_INT,MPI_SUM,0,mpi_comm);
      IF_RANK0 cout << "   > injecting " << count << " particles for this static injector" << endl; 
      
    } else if ( InjectorType == EACH_CELL_TYPE ) {

      // this injector seeds particles on the center of each cv

      assert(np_in_cell>0);

      int8 * cvora = NULL;
      buildXora(cvora,solver->ncv);
      assert( cvora[mpi_rank+1] - cvora[mpi_rank] == solver->ncv );
      int ip_rank = cvora[mpi_rank] * np_in_cell; // start from the location for your rank

      for (int icv = 0; icv < solver->ncv; icv++) {
        if ((solver->x_vv[icv][0] >= x_min) and (solver->x_vv[icv][0] <= x_max)) {
          if ((solver->x_vv[icv][1] >= y_min) and (solver->x_vv[icv][1] <= y_max)) {
            if ((solver->x_vv[icv][2] >= z_min) and (solver->x_vv[icv][2] <= z_max)) {
              int icv_closest;
              const int isInside = pointIsInside(solver->x_vv[icv], icv_closest);
              if ( isInside ) {
                assert(icv == icv_closest);
                for (int ip = 0; ip < np_in_cell; ip++) {
                  lspDataVec.resize(lspDataVec.size()+1);
                  int iback = lspDataVec.size()-1;
                  lspDataVec[iback].icv = icv_closest;
                  lspDataVec[iback].mp = mpnew[ip_rank];
                  lspDataVec[iback].Tp = Tpnew[ip_rank];
                  FOR_I3 lspDataVec[iback].xp[i] = solver->x_vv[icv][i];
                  FOR_I3 lspDataVec[iback].up[i] = GeomData[6+i];
                  if (pShape == ELLIPSOID_SHAPE) {
                    if ( (axis_e>0)&&(axis_f>0) ) {
                      lspDataVec[iback].e = axis_e;
                      lspDataVec[iback].f = axis_f;
                    } 
                    else if (axis_std>0) {
                      const double dmin = 1e-10;
                      lspDataVec[iback].e = sampleNormalDist(1, axis_std, dmin, 1);
                      lspDataVec[iback].f = sampleNormalDist(1, axis_std, dmin, 1);
                    }
                  }
                  else {
                    lspDataVec[iback].e = 1.;
                    lspDataVec[iback].f = 1.;
                  }
                  ip_rank++;
                  my_np_count++;
                }
              }
            }
          }
        }
      }
      delete[] cvora;

      // we are finished
      int count;
      MPI_Reduce(&my_np_count,&count,1,MPI_INT,MPI_SUM,0,mpi_comm);
      IF_RANK0 cout << "   > injecting " << count << " particles for this static injector" << endl; 

    } else if ( (InjectorType == COUNT_TYPE) or (InjectorType == MASS_TYPE) ) {

      // for all the other InjectorTypes: COUNT or MASS, we have a random process

      //////////////////////////////////////////////////////
      // find out how many particles goes to each processor
      //////////////////////////////////////////////////////

      // compute what portion of the box volume belongs to this processor
      double my_total_volume = 0;
      const int n_max = 100;
      double dx = (x_max - x_min)/double(n_max);
      double dy = (y_max - y_min)/double(n_max);
      double dz = (z_max - z_min)/double(n_max);
      double d_max = max(max(dx,dy),dz);
      int nx = int(floor((x_max-x_min)/d_max));
      int ny = int(floor((y_max-y_min)/d_max));
      int nz = int(floor((z_max-z_min)/d_max));
      for (int ix = 0; ix <= nx; ix++) {
        double x = x_min + double(ix)*d_max;
        double x_next = min( (x + d_max) , x_max );
        for (int iy = 0; iy <= ny; iy++) {
          double y = y_min + double(iy)*d_max;
          double y_next = min( (y + d_max) , y_max );
          for (int iz = 0; iz <= nz; iz++) {
            double z = z_min + double(iz)*d_max;
            double z_next = min( (z + d_max) , z_max );

            double X_mid[3] = {(x+x_next)/2., (y+y_next)/2., (z+z_next)/2.};
            vector<int> cvList;
            solver->cvAdt->buildListForPoint(cvList, X_mid);
            const int icv = getClosestICVFromList(cvList,X_mid);

            if (icv >= 0) my_total_volume += d_max*d_max*d_max;
          }
        }
      }
              
      double * total_volume_ranks = NULL;
      IF_RANK0 total_volume_ranks = new double [mpi_size];
      MPI_Gather(&my_total_volume, 1, MPI_DOUBLE, total_volume_ranks, 1, MPI_DOUBLE, 0, mpi_comm);

      int * npnew_ranks;
      IF_RANK0 {

        double total_volume_sum = 0;
        for (int irank = 0; irank < mpi_size; irank++) {
          total_volume_sum += total_volume_ranks[irank];
        }
        
        npnew_ranks = new int [mpi_size];
        int npnew_rem = npnew;
        for (int irank = 0; irank < mpi_size; irank++) {
          npnew_ranks[irank] = (int)floor(total_volume_ranks[irank] / total_volume_sum * double(npnew));
          npnew_rem -= npnew_ranks[irank];
        }
        
        for (int irank = 0; irank < mpi_size; irank++) {
          if (npnew_ranks[irank] > 0) {
            npnew_ranks[irank] += npnew_rem;
            break;
          }
        }

        delete[] total_volume_ranks;
      }
  
      int my_npnew;
      MPI_Scatter(npnew_ranks, 1, MPI_INT, &my_npnew, 1, MPI_INT, 0, mpi_comm);
      assert(my_npnew >= 0);
      IF_RANK0 delete[] npnew_ranks;

      int ip_new = 0;
      while (ip_new < my_npnew) {
        // set xp_seed
        double xp_seed[3];
        xp_seed[0] = x_min + (x_max - x_min) * double(rand()) / double(RAND_MAX);
        xp_seed[1] = y_min + (y_max - y_min) * double(rand()) / double(RAND_MAX);
        xp_seed[2] = z_min + (z_max - z_min) * double(rand()) / double(RAND_MAX);
        
        int icv_closest;
        const int isInside = pointIsInside(xp_seed, icv_closest);
        if ( isInside ) {
          lspDataVec.resize(lspDataVec.size()+1);
          int iback = lspDataVec.size()-1;
          lspDataVec[iback].icv = icv_closest;
          lspDataVec[iback].mp = mpnew[ip_new];
          lspDataVec[iback].Tp = Tpnew[ip_new];
          FOR_I3 lspDataVec[iback].xp[i] = xp_seed[i];
          FOR_I3 lspDataVec[iback].up[i] = GeomData[6+i];
          if (pShape == ELLIPSOID_SHAPE) {
            if ( (axis_e>0)&&(axis_f>0) ) {
              lspDataVec[iback].e = axis_e;
              lspDataVec[iback].f = axis_f;
            } 
            else if (axis_std>0) {
              const double dmin = 1e-10;
              lspDataVec[iback].e = sampleNormalDist(1, axis_std, dmin, 1);
              lspDataVec[iback].f = sampleNormalDist(1, axis_std, dmin, 1);
            }
          }
          else {
            lspDataVec[iback].e = 1.;
            lspDataVec[iback].f = 1.;
          }
          ip_new++;
          my_np_count++;
        }
      
      }

      // we are finished
      int count;
      MPI_Reduce(&my_np_count,&count,1,MPI_INT,MPI_SUM,0,mpi_comm);
      IF_RANK0 cout << "   > injecting " << count << " particles for this static injector" << endl; 

    } else {
      CERR("Injector type wasn't found");
    }

  } else {
    CERR("GeomType wasn't found");
  }

  if (overWriteVelocity_b) overWriteVelocity(lspDataVec,from_ind,true);

  delete[] mpnew;
  delete[] Tpnew;
  return;
}

int StaticLspInjector::getClosestICVFromList(vector<int> &cvList, const double* xp) {
  double dist_closest = 1e20;
  int icv_closest = -1;
  for (int i = 0; i < cvList.size(); ++i) {
    const int icv = cvList[i];
    double dist = DIST(xp,solver->x_cv[icv]);
    if (dist < dist_closest) {
      dist_closest = dist;
      icv_closest = icv;
    }
  }
  return(icv_closest);
}

void LspInjector::overWriteVelocity(vector<LspData>& lspDataVec,int from_ind, bool verbose) {
    
  IF_RANK0 
    if (verbose)
      cout << "   > over writing injector velocity..." << endl;

  switch (overWriteVelocity_mode) {
    case DEFAULT:
    {
      for (int ivec = from_ind; ivec < lspDataVec.size(); ivec++) {
        FOR_I3 lspDataVec[ivec].up[i] = u[lspDataVec[ivec].icv][i];
      }
      break;
    }
    case MULT:
    {
      for (int ivec = from_ind; ivec < lspDataVec.size(); ivec++) {
        FOR_I3 lspDataVec[ivec].up[i] = u[lspDataVec[ivec].icv][i] * OWVData[0];
      }
      break;
    }
    case MAG:
    {
      for (int ivec = from_ind; ivec < lspDataVec.size(); ivec++) {
        double dir[3]; 
        FOR_I3 dir[i] = u[lspDataVec[ivec].icv][i];
        NORMALIZE(dir);
        FOR_I3 lspDataVec[ivec].up[i] = dir[i] * OWVData[0];
      }
      break;
    }
    case MAX:
    {
      for (int ivec = from_ind; ivec < lspDataVec.size(); ivec++) {
        double u_mag = MAG(u[lspDataVec[ivec].icv]);
        FOR_I3 lspDataVec[ivec].up[i] = u[lspDataVec[ivec].icv][i] * ( (u_mag > OWVData[0]) ? (OWVData[0]/u_mag):1.0);
      }
      break;
    }
  }


}

int LspInjector::pointIsInside(const double xp[3], int &icv_ret) {

  icv_ret = -1;

  // returns:
  // 1: is inside   icv_ret: icv_closest
  // 0: is outside  icv_ret: -1
  solver->ensureCvAdt();

  vector<int> cvList;
  solver->cvAdt->buildListForPoint(cvList, xp);

  int icv_closest = -1;
  double d2_closest;
  for (int ivec = 0; ivec < cvList.size(); ivec++) {
    const int icv = cvList[ivec];
    const double d2 = DIST2(solver->x_vv[icv], xp);
    if ((icv_closest==-1)||(d2<d2_closest)) {
      d2_closest = d2;
      icv_closest = icv;
    }
  }
  if (icv_closest==-1) {
    // could not find a cell to own the point
    return 0;
  }

  // find the closest cell to the point
  assert((icv_closest>=0)&&(icv_closest<solver->ncv));
  int done = 0;
  int counter = 0;
  int counter_max = 1000;
  while (!done) {
    
    done = 1;
    int icv_nbr_closest = -1;
    double d2_nbr_closest;
    for (int coc = solver->cvocv_i[icv_closest]+1; coc != solver->cvocv_i[icv_closest+1]; coc++) {
      const int icv_nbr = solver->cvocv_v[coc];
      assert(icv_nbr!=icv_closest);
      const double d2_nbr = DIST2(solver->x_vv[icv_nbr], xp);
      if ((icv_nbr_closest==-1)||(d2_nbr < d2_nbr_closest)) {
        d2_nbr_closest = d2_nbr;
        icv_nbr_closest = icv_nbr;
      }
    }

    if (d2_nbr_closest < d2_closest) {
      // we got a closer neighbor
      done = 0;
      if (icv_nbr_closest >= solver->ncv) {
        // the neighbor is a ghost so we are outside
        return 0;
      } else {
        icv_closest = icv_nbr_closest;
        d2_closest = d2_nbr_closest;
      }
    }

    ++counter;
    if (counter>counter_max) cout << "WARNING, reached maximum iteration for rank: " << mpi_rank << " counter: " << counter << " d2_closest: " << d2_closest << " icv_closest: " << icv_closest << " xp: " << COUT_VEC(xp) << endl;

  }

  // now we should have icv_closest with the point inside of it
  assert((icv_closest>=0)&&(icv_closest<solver->ncv));

  // check if the closest cell has a boundary surface and decide if it is inside
  const int my_nbf = solver->bfocv_i[icv_closest+1] - solver->bfocv_i[icv_closest];
  if (my_nbf==0) {
    // no boundary faces so the point is inside
    icv_ret = icv_closest;
    return 1;
  }
  else {
    
    // and we need the closest ibf,ist_ss...
    // We need a length scale...
    const double tol2 = 1.0E-10*pow(solver->vol_cv[icv_closest],2.0/3.0);
    int ibf_closest = -1;
    int ist_ss_closest;
    double d2_closest,dist_closest;
    //double dx_closest[3],normal_closest[3];
    for (int boc = solver->bfocv_i[icv_closest]; boc != solver->bfocv_i[icv_closest+1]; ++boc) {
      const int ibf = solver->bfocv_v[boc];
      for (int sob = solver->sstobf_i[ibf]; sob != solver->sstobf_i[ibf+1]; ++sob) {
        const int ist_ss = solver->sstobf_v[sob];
              
        // the corner coords of this sub-surface tri are (recall that the
        // subsurface tri already has any periodic transfrom taken care of)...
        const double * const x0 = solver->subSurface->xp[solver->subSurface->spost[ist_ss][0]];
        const double * const x1 = solver->subSurface->xp[solver->subSurface->spost[ist_ss][1]];
        const double * const x2 = solver->subSurface->xp[solver->subSurface->spost[ist_ss][2]];
        
        double xp_tri[3]; MiscUtils::getClosestPointOnTriRobust(xp_tri,xp,x0,x1,x2);
        
        const double dx[3] = DIFF(xp_tri, xp);
        const double d2 = DOT_PRODUCT(dx,dx);
        
        if ((ibf_closest == -1)||(d2 < d2_closest+tol2)) {
          const double normal[3] = TRI_NORMAL_2(x0,x1,x2);
          const double normal_mag = MAG(normal);
          assert(normal_mag != 0); // the subsurface should not have zero-area tris!
          const double dist = DOT_PRODUCT(dx,normal)/normal_mag;
          if ((ibf_closest == -1)||(d2 < d2_closest-tol2)||(fabs(dist) > fabs(dist_closest))) {
            ibf_closest = ibf;
            ist_ss_closest = ist_ss;
            d2_closest = d2;
            //FOR_I3 dx_closest[i] = dx[i];
            //FOR_I3 normal_closest[i] = normal[i];
            dist_closest = dist;
          }
        }
      }
    }
    assert(ist_ss_closest >= 0);
    if (dist_closest >= 0.0) {
      // dist with positive sign means we are inside the fluid volume...
      icv_ret = icv_closest;
      return 1;
    } else {
      // dist with negative sign means we are outside the fluid volume...
      return 0;
    }
  }

  // we should never reach here
  assert(0);
  return 0;
}

double error_function(const double eta) {

  double z = fabs(eta);
  double t = 1.0/(1.0+0.5*z);
  double erfcc = t*exp(-z*z-1.26551223+t*(1.00002368+t*(0.37409196+
              t*(0.09678418+t*(-0.18628806+t*(0.27886807+t*(-1.13520398+
                        t*(1.48851587+t*(-0.82215223+t*0.17087277)))))))));
  if (eta < 0.0) erfcc = 2.0-erfcc;

  double err_func = 1.0 - erfcc;
  return err_func;

}

void project_to_plane(double * x_pro, const double x[3], const double xp[3], const double np[3]) {
  
  double dx[3] = DIFF(x,xp);
  double d = DOT_PRODUCT(dx,np);
  double mag2_np = DOT_PRODUCT(np,np);
  FOR_I3 x_pro[i] = x[i] - d * np[i] / mag2_np;
}
//void LspInjector::getRotationMatrix(double (* Rot)[3], const double n[3], const double theta) {

//}

double sampleNormalDist(const double mean, const double std, const double xmin, const double xmax) {

  // using 0.5(1+erf([x-mu]/[sigma*sqrt(2)])) as the cdf

  double err = 1e20;
  double tol = 1e-7;
  // centralize and normalize x to have 0 mean and variance 0.5
  double eta_min = (xmin-mean)/(std*sqrt(2.));
  double eta_max = (xmax-mean)/(std*sqrt(2.));
  double eta_mid;

  const double f_min = 0.5 * ( 1.0 + error_function(eta_min));
  const double f_max = 0.5 * ( 1.0 + error_function(eta_max));
  double rand_val = (f_max-f_min) * double(rand())/(double(RAND_MAX)+1.) + f_min; // -f_min < rand_val < f_max

  int iter = 0;
  while ( err > tol ) {
    
    eta_mid = (eta_min+eta_max)/2.;
    double f_mid = 0.5 * ( 1.0 + error_function(eta_mid));

    if (f_mid > rand_val) {
      eta_max = eta_mid;  
    }
    else {
      eta_min = eta_mid;
    }

    err = fabs(f_mid - rand_val);
    iter++;
    if (iter>1000)
      cout << " > WARNING: sampleNormalDist cannot converge. rand val: " << rand_val << " iteration: " << iter << endl;
  }

  // scale and shift eta
  double x_rand = eta_mid*(std*sqrt(2.)) + mean;
  x_rand = min(x_rand, xmax);
  x_rand = max(x_rand, xmin);
  return x_rand;

}

