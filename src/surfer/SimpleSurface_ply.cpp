#include "SimpleSurface.hpp"
#include "PlyReader.hpp"

int SimpleSurface::addPly(const string& filename) {

  // FILE * fp = fopen(filename.c_str(),"rb");
  // if (fp == NULL) {
  //   cerr << "Error: cannot open file: " << filename << endl;
  //   return -1;
  // }

  // current surface values
  const int ss_nsp0 = nsp;
  const int ss_nst0 = nst;

  PlyReader fr;
  try {
    fr.init(filename);
  }
  catch (...) {
    return -1;
  }

  string token;

  COUT1(" > reading header information");
  token = fr.getNextTokenAsString();
  if (token != "ply") {
    CWARN("this file doesn't appear to be a ply file");
    return -1;
  }

  // format specification
  token = fr.getNextTokenAsString();
  if (token != "format") {
    CWARN("this file doesn't specify its format properly");
    return -1;
  }
  else {
    token = fr.getNextTokenAsString();
    if (token == "ascii") {
      COUT2(" > ascii file");
    }
    else if (token == "binary_little_endian") {
      COUT2(" > binary file (little endian)");
      fr.setBinary();
    }
    else if (token == "binary_big_endian") {
      fr.setBinary();
      fr.setByteSwap();
      COUT2(" > binary file (big endian)");
    }
    else {
      CWARN("unrecognized format: " << token);
      return -1;
    }
  }
  token = fr.getNextTokenAsString();
  assert(token == "1.0");  // all ply files use hardcoded version

  token = fr.getNextTokenAsString();
  while (token != "end_header") {

    if (token == "comment") {
      fr.goToNextLine();
    }
    else if (token == "element") {
      const string name = fr.getNextTokenAsString();
      const int count = fr.getNextTokenAsInt();
      fr.addElement(name,count);

      // special treatment for surface mesh vars
      if (name == "vertex") {
        nsp += count;
      }
      if (name == "face") {
        nst += count;
      }
    }
    else if (token == "property") {
      string name;
      string var_type = fr.getNextTokenAsString();

      if (var_type == "list") {
        const string list_var_type = fr.getNextTokenAsString();
        var_type = fr.getNextTokenAsString();
        name = fr.getNextTokenAsString();
        fr.addProperty(list_var_type,var_type,name);
      }
      else {
        // not a list
        name = fr.getNextTokenAsString();
        fr.addProperty(var_type,name);
      }
    }
    else {
      COUT2(" > unrecognized ply header token: " << token);
    }
    token = fr.getNextTokenAsString();
  }

  // header summary
  for (int e=0, lim=fr.elementsVec.size(); e<lim; ++e) {
    cout << " > element \""<< fr.elementsVec[e].getName() << "\" has properties: ";
    cout << fr.elementsVec[e].propertiesVec[0].getName();
    for (int p=1, lim=fr.elementsVec[e].propertiesVec.size(); p<lim; ++p) {
      cout << ","<< fr.elementsVec[e].propertiesVec[p].getName();
    }
    cout << endl;
  }

  // allocate new sizes for nsp/nst surface objects
  growNspData(nsp,ss_nsp0);
  growNstData(nst,ss_nst0);

  for (int e=0, lim=fr.elementsVec.size(); e<lim; ++e) {
    if (fr.elementsVec[e].getName() == "vertex") {
      COUT2(" > reading " << fr.elementsVec[e].getCount() << " vertices");
      for (int isp=ss_nsp0; isp<nsp; ++isp) {
        for (int p=0, lim=fr.elementsVec[e].propertiesVec.size(); p<lim; ++p) {
          if (fr.elementsVec[e].propertiesVec[p].getName() == "x") {
            xsp[isp][0] = fr.elementsVec[e].propertiesVec[p].getDouble();
          }
          else if (fr.elementsVec[e].propertiesVec[p].getName() == "y") {
            xsp[isp][1] = fr.elementsVec[e].propertiesVec[p].getDouble();
          }
          else if (fr.elementsVec[e].propertiesVec[p].getName() == "z") {
            xsp[isp][2] = fr.elementsVec[e].propertiesVec[p].getDouble();
          }
          // skip unknown properties
          else if (fr.elementsVec[e].propertiesVec[p].isInt()) {
            fr.elementsVec[e].propertiesVec[p].getInt();
          }
          else {
            fr.elementsVec[e].propertiesVec[p].getDouble();
          }
        }
      }
    }
    if (fr.elementsVec[e].getName() == "face") {
      COUT2(" > reading " << fr.elementsVec[e].getCount() << " faces");
      for (int ist=ss_nst0; ist<nst; ++ist) {
        for (int p=0, lim=fr.elementsVec[e].propertiesVec.size(); p<lim; ++p) {

          if (fr.elementsVec[e].propertiesVec[p].getName() == "vertex_indices") {
            const int n_vals = fr.elementsVec[e].propertiesVec[p].getListInt();
            if (n_vals != 3) {
              CWARN("file specified polygon greater than triangle; currently not supported");
              return -1;
            }
            for (int i=0; i<n_vals; ++i) {
              spost[ist][i] = fr.elementsVec[e].propertiesVec[p].getInt();
            }
          }
          else {
            // skip through other data
            if (fr.elementsVec[e].propertiesVec[p].isList()) {
              const int n_vals = fr.elementsVec[e].propertiesVec[p].getListInt();
              if (fr.elementsVec[e].propertiesVec[p].isInt()) {
                for (int i=0; i<n_vals; ++i) fr.elementsVec[e].propertiesVec[p].getInt();
              }
              else {
                for (int i=0; i<n_vals; ++i) fr.elementsVec[e].propertiesVec[p].getDouble();
              }
            }
            else {
              // not a list
              if (fr.elementsVec[e].propertiesVec[p].isInt()) {
                fr.elementsVec[e].propertiesVec[p].getInt();
              }
              else {
                fr.elementsVec[e].propertiesVec[p].getDouble();
              }
            }
          }
        }
      }
    }
  }

  const int new_zone = zoneVec.size();
  size_t filename_pos = filename.find_last_of("/\\");
  if (filename_pos!=string::npos) zoneVec.push_back(SurfaceZone(filename.substr(filename_pos+1)));
  else zoneVec.push_back(SurfaceZone(filename));

  // no zones specified in ply yet, so all get new zone based on filename
  for (int ist=ss_nst0; ist < nst; ++ist) {
    znost[ist] = new_zone;
  }

  if (ss_nst0 > 0) {
    // if adding surface, properly offset node indices
    for (int ist = ss_nst0; ist < nst; ++ist) {
      FOR_I3 spost[ist][i] += ss_nsp0;
    }
  }

  COUT1(" > initFromPly done");
  // ====================================================
  // make sure all zones were set...
  for (int ist = 0; ist < nst; ++ist) {
    assert((znost[ist] >= 0)&&(znost[ist] < int(zoneVec.size())));
  }

  return 0;  // success
}
