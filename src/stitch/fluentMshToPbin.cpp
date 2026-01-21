#include "CTI.hpp"
using namespace CTI;

#include "DataContainer.hpp"

void processMsh() {

  FOR_PARAM_MATCHING("MSH") {

    string msh_filename = param->getString();

    if (mpi_rank == 0) cout << " > working on msh/cas file: " << msh_filename << endl;

    DataContainer * dc = new DataContainer();

    // look jsut for byte-swap
    int byte_swap = 0;  // 0: unset, 1: no swap, 2: swap
    int iarg = 1;
    while (iarg < param->size()) {
      string token = MiscUtils::toUpperCase(param->getString(iarg++));
      if (token == "BYTE_SWAP") {
        string tok = MiscUtils::toUpperCase(param->getString(iarg++));
        if (tok == "TRUE") byte_swap = 2;
        else if (tok == "FALSE") byte_swap = 1;
      }
    }

    dc->initPointsFromFluent(msh_filename,byte_swap);

    // parse remaining args...
    iarg = 1;
    while (iarg < param->size()) {
      string token = MiscUtils::toUpperCase(param->getString(iarg++));
      if (token == "SCALE") {
        const double factor = param->getDouble(iarg++);
        COUT1(" > SCALE currently not supported; skipping");
        // for (int ino = 0; ino < msh->nno; ++ino) {
        //   FOR_I3 msh->x_no[ino][i] *= factor;
        // }
      }
      else if (token == "TRANSLATE") {
        double dx[3];
        FOR_I3 dx[i] = param->getDouble(iarg++);
        COUT1(" > applying TRANSLATE " << dx[0] << " " << dx[1] << " " << dx[2]);
        dc->translate(dx[0],dx[1],dx[2]);
      }
      else if (token == "ROTATE") {
        COUT1(" > applying ROTATE");
        double x0[3];
        FOR_I3 x0[i] = param->getDouble(iarg++);
        COUT1("    > ref. point " << x0[0] << " " << x0[1] << " " << x0[2]);
        double n0[3];
        FOR_I3 n0[i] = param->getDouble(iarg++);
        COUT1("    > axis dir. " << n0[0] << " " << n0[1] << " " << n0[2]);
        const double alpha_deg = param->getDouble(iarg++);
        dc->rotate(x0,n0,alpha_deg);
      }
      else if (token == "MIRROR") {
        COUT1(" > applying MIRROR");
        double x0[3];
        FOR_I3 x0[i] = param->getDouble(iarg++);
        COUT1("    > ref. point " << x0[0] << " " << x0[1] << " " << x0[2]);
        double n0[3];
        FOR_I3 n0[i] = param->getDouble(iarg++);
        COUT1("    > plane normal " << n0[0] << " " << n0[1] << " " << n0[2]);
        dc->mirror(x0,n0);
      }
      else {
        CWARN("unrecognized MSH token: " << token);
      }
    }

    // write to pbin
    size_t last_dot = msh_filename.find_last_of('.');
    string pbin_filename = msh_filename.substr(0,last_dot) + ".pbin";
    dc->writePointsBinary(pbin_filename);

    delete dc;
  }
  cout << "done" << endl;

}

int main(int argc, char * argv[]) {

  try {

    CTI_Init(argc,argv,"fluentMshToPbin.in");

    {

      processMsh();

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
