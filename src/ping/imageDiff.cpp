#include "CTI.hpp"
using namespace CTI;
#include "MiscUtils.hpp"
#include "PngData.hpp"

int main(int argc, char * argv[]) {

  try {
    if (argc < 4) throw(-1);

    int iarg = 1;
    string metric = "AE";
    string files[2] = {"",""};
    double tol = -1.0;
    bool force_image = false;
    while (iarg < argc) {
      const string token=argv[iarg++];
      if (token == "METRIC") {
        metric = argv[iarg++];
      }
      else if (token == "IMAGES" || token == "FILES") {
        files[0] = argv[iarg++];
        files[1] = argv[iarg++];
      }
      else if (token == "FORCE_IMAGE") {
        force_image = true;
      }
      else if (token == "DIFF_TOL") {
        tol = atof(argv[iarg++]);
      }
      else {
        CWARN("unrecognized imageDiff parameter \"" << token << "\"; skipping");
      }
    }

    PngData * img0;
    PngData * img1;
    try {
      img0 = new PngData();
      img0->read(files[0].c_str());
      img0->finalizeRead();

      img1 = new PngData();
      img1->read(files[1].c_str());
      img1->finalizeRead();
    }
    catch (...) {
      throw(-2);
    }

    // use reference filename for prefix
    string outputPrefix = MiscUtils::getPrefix(files[0]);
    COUT2(" > output file prefix will be \"" << outputPrefix << "\"");

    return img1->diffWithStats(img0,metric,outputPrefix,tol,force_image);

  }
  catch (int e) {
    switch (e) {
      case -1:
        cout << "Error: expecting 2 arguments: ./imageDiff.exe FILES <image-file> <image-file> [optional args]" << endl;
        break;
      case -2:
        cout << "Error: could not properly read input files" << endl;
        break;
      default:
        cout << "Error: unknown type" << endl;
    }
    return e;
  }
  catch(...) {
    cout << "Error: unexpected exception was hit; exiting" << endl;
    return -10;
  }
}
