#include "SimpleSurface.hpp"
#include "UgridReader.hpp"

int SimpleSurface::addUgrid(const string& filename, const string& bc_filename) {

  // current surface values
  const int ss_nsp0 = nsp;
  const int ss_nst0 = nst;

  UgridReader * fr;

  // parse file pre-suffix for endianness and binary length
  std::size_t dots[2] = {std::string::npos,std::string::npos};
  std::size_t found = filename.rfind(".");
  if (found!=std::string::npos) {
    dots[1] = found-1;  // position before the period
    found = filename.rfind(".",(dots[1]));
    if (found!=std::string::npos) {
      dots[0] = found+1;  // position after dot
    }
  }

  // initialize either ascii or binary reader based on prefix existence
  if (dots[0] == std::string::npos || dots[1] == std::string::npos) {
    WUI(WARN," ! did not find file type pre-suffix; assuming ASCII formatting");
    try {
      fr = new UgridAsciiReader(filename);
    }
    catch (...) {
      return -1;
    }
  }
  else {
    // binary
    try {
      fr = new UgridBinaryReader(filename);

      const string type_string = filename.substr(dots[0],(dots[1]-dots[0]+1));
      COUT1(" > format type: " << type_string);

      if (type_string == "r4") {
        fr->setFortran();
        fr->setByteSwap();
      }
      else if (type_string == "r8") {
        fr->setFortran();
        fr->setByteSwap();
        fr->setBytes(8);
      }
      else if (type_string == "lr4") {
        fr->setFortran();
      }
      else if (type_string == "lr8") {
        fr->setFortran();
        fr->setBytes(8);
      }
      else if (type_string == "b4") {
        fr->setByteSwap();
      }
      else if (type_string == "b8") {
        fr->setByteSwap();
        fr->setBytes(8);
      }
      else if (type_string == "lb4") {
      }
      else if (type_string == "lb8") {
        fr->setBytes(8);
      }
      else {
        WUI(WARN,"unrecognized file type pre-suffix " << type_string << "; skipping");
        return -1;
      }
    }
    catch (...) {
      return -1;
    }
  }

  fr->readHeader(nsp,nst);

  if (ss_nsp0 == 0) {
    COUT1(" > nsp: " << nsp);  // just the new points
  }
  else {
    nsp += ss_nsp0;
    COUT1(" > nsp: " << nsp << " (+" << (nsp-ss_nsp0) << ")");
  }

  if (ss_nst0 == 0) {
    COUT1(" > nst: " << nst);  // just the new points
  }
  else {
    nst += ss_nst0;
    COUT1(" > nst: " << nst << " (+" << (nst-ss_nst0) << ")");
  }

  // read in new nodes
  growNspData(nsp,ss_nsp0);
  fr->readNodes(&xsp[ss_nsp0]);

  // read in new tris
  growNstData(nst,ss_nst0);
  fr->readTrisAndQuads(&spost[ss_nst0],&znost[ss_nst0]);

  // properly offset and account for 1-indexing in ugrid files
  for (int ist=ss_nst0; ist<nst; ++ist) {
    FOR_I3 spost[ist][i] += ss_nsp0-1;
  }

  // create a map from specified ugrid zone indices to new contiguous ones
  map<int,int> ugridzone2zone;
  map<int,int>::iterator it;

  // no map file provided; use indices in file for zone names
  for (int ist=ss_nst0; ist<nst; ++ist) {
    it = ugridzone2zone.find(znost[ist]);

    if (it != ugridzone2zone.end()) {
      znost[ist] = it->second;
    }
    else {
      // not found, so add to zoneVec and map
      ugridzone2zone.insert(pair<int,int> (znost[ist],zoneVec.size()));
      stringstream oss; oss << znost[ist];
      const string zonename = "ugridzone_" + oss.str();
      znost[ist] = zoneVec.size();

      zoneVec.push_back(SurfaceZone(zonename));
    }
  }

  if (bc_filename.size()) {
    COUT1(" > mapping zone names from bc file: " << bc_filename << " (assumed mcfd.bc file)");

    BcAsciiReader bcr;
    try {
      bcr.init(bc_filename);
    }
    catch (...) {
      WUI(WARN,"problem parsing zone-mapping file; skipping");
      return -1;  // TODO just fail this and go to else{} instead of return...
    }

    map<int,string> ugridzone2name;
    bcr.getZonenameMap(ugridzone2name,"SLUGG");  // extensible for other filetypes eventually

    for (map<int,string>::iterator it2=ugridzone2name.begin(); it2!=ugridzone2name.end(); ++it2) {
      it = ugridzone2zone.find(it2->first);
      if (it != ugridzone2zone.end()) {
        // found this ugrid-zone
        const int izn = it->second;
        zoneVec[izn].setName(it2->second);
      }
      else {
        WUI(WARN,"ugrid zone " << it2->first << " wasn't found; skipping renaming of this zone to " << it2->second);
      }
    }
  }

  // ====================================================
  // make sure all zones were set...
  for (int ist = 0; ist < nst; ++ist) {
    assert((znost[ist] >= 0)&&(znost[ist] < int(zoneVec.size())));
  }

  return 0;
}
