#include "SimpleSurface.hpp"

#define FAST_OBJ_BUFFER_SIZE 65535
#define FAST_OBJ_TOKEN_SIZE 65535

char * getNextObjToken(bool& endOfLine, bool& endOfFile, FILE * fp,int& pos,int& max_pos,char * buf, bool& max_is_eof) {

  char * token = NULL;
  while (1) {

    if (pos >= max_pos) {

      assert(pos == max_pos);
      if (token) {
        // if the token is active, but not completed, then we need to shift the token part of
        // the current buf (i.e. the end) to the start of the buf, and read the next part in...
        //cout << "pos >= max_pos: " << pos << " " << max_pos << " token[0]: " << token[0] << " max_pos-token+buf: " << max_pos-int(token-buf) << endl;
        pos = max_pos - int(token-buf);
        if (token != buf) {
          memmove(buf,token,pos);
          token = buf; // reset to start...
        }
      }
      else {
        pos = 0;
      }

      max_pos = pos + fread(buf+pos, sizeof(char), FAST_OBJ_BUFFER_SIZE, fp);
      if (feof(fp)) max_is_eof = true;

      if (max_pos == pos) {
        buf[pos] = '\0';
        endOfLine = false;
        if (max_is_eof) endOfFile = true;
        return token;
      }

      /*
      cout << "buf is: \"";
      for (int i = 0; i < max_pos; ++i)
      cout << buf[i];
      cout << "\"" << endl;
      getchar();
      */
    }

    //cout << "getNextStlToken pos: " << pos << " c: \"" << buf[pos] << "\"" << endl;
    const char c = buf[pos++];

    if ( (c == ' ') || (c == '\t') || (c == '\n') || (c == 13) ) {
      // c is a whitespace character - this represents either the termination
      // of a token, or some space at the start of the next token...
      if ((c == '\n') || (c == 13)) endOfLine = true;
      else endOfLine = false;

      if (token || endOfLine) break;
    }
    else if (!token) {
      // any other character is the start of the token...
      token = buf+pos-1;
      endOfLine = false;
    }

  }

  // terminate string and return...
  buf[pos-1] = '\0';
  return token;

}

void moveToStartOfNextLine(bool& endOfFile,FILE * fp,int& pos,int& max_pos,char * buf,bool& max_is_eof) {

  char * token = NULL;
  while (1) {

    if (pos >= max_pos) {

      assert(pos == max_pos);
      if (token) {
        // if the token is active, but not completed, then we need to shift the token part of
        // the current buf (i.e. the end) to the start of the buf, and read the next part in...
        //cout << "pos >= max_pos: " << pos << " " << max_pos << " token[0]: " << token[0] << " max_pos-token+buf: " << max_pos-int(token-buf) << endl;
        pos = max_pos - int(token-buf);
        if (token != buf) {
          memmove(buf,token,pos);
          token = buf; // reset to start...
        }
      }
      else {
        pos = 0;
      }

      max_pos = pos + fread(buf+pos, sizeof(char), FAST_OBJ_BUFFER_SIZE, fp);
      if (feof(fp)) max_is_eof = true;

      if (max_pos == pos) {
        buf[pos] = '\0';
        if (max_is_eof) endOfFile = true;
        return;
      }

    }

    const char c = buf[pos++];

    if ((c == '\n') || (c == 13)) {
      // c is a line end or carraige return
      if (token) break;
    }
    else if (!token) {
      // any other character is the start of the token...
      token = buf+pos-1;
    }
  }

  // terminate string and return...
  buf[pos-1] = '\0';
}

struct TriNodes {
public:
  int nodes[3];
  int zone;
  TriNodes() {
    FOR_I3 nodes[i] = -1;
    zone = -1;
  }
};

// obj files are always ascii. Their binary counterparts are .mod files.
// This is not a complete obj file reader. We ignore vertex normals and vertex texture data.
int SimpleSurface::addObj(const string& filename) {

  FILE * fp = NULL;
  const int file_err = MiscUtils::openFile(&fp,filename,"rb");
  if (file_err != 0) return file_err;

  // current surface values
  const int ss_nsp0 = nsp;
  const int ss_nst0 = nst;
  vector<DoubleVertex> vertexVec;
  vector<TriNodes> TriNodesVec;

  COUT1(" > only ascii obj files are supported (does not support binary .mod files)");

  fseek(fp,0,SEEK_SET);
  COUT1(" > scanning file for number of vertices and tris");

  int pos = 0;
  int max_pos = 0;
  char * buf = new char[FAST_OBJ_BUFFER_SIZE+FAST_OBJ_TOKEN_SIZE];
  char * firstTokenOfLine = NULL;
  bool endOfLine = false;
  bool endOfFile = false;
  bool max_is_eof = false;
  int n_zones = -1;
  int n_faces = 0;
  int f_split = 0;

  while (!endOfFile) {
    firstTokenOfLine = getNextObjToken(endOfLine,endOfFile,fp,pos,max_pos,buf,max_is_eof);

    if (firstTokenOfLine != NULL) {
      if (strcmp(firstTokenOfLine,"v") == 0) {
        nsp++;
        vertexVec.push_back(DoubleVertex());
        sscanf(getNextObjToken(endOfLine,endOfFile,fp,pos,max_pos,buf,max_is_eof),"%lf",&vertexVec.back().x);
        sscanf(getNextObjToken(endOfLine,endOfFile,fp,pos,max_pos,buf,max_is_eof),"%lf",&vertexVec.back().y);
        sscanf(getNextObjToken(endOfLine,endOfFile,fp,pos,max_pos,buf,max_is_eof),"%lf",&vertexVec.back().z);

        //DEBUG
        // cout << "vertex: " << vertexVec.back().x << " " << vertexVec.back().y << " " << vertexVec.back().z << endl;
      }
      else if (strcmp(firstTokenOfLine,"f") == 0) {
        ++n_faces;
        vector<int> my_vertices;
        char * token = NULL;
        char * str = NULL;

        while (!endOfLine && !endOfFile) {
          token = getNextObjToken(endOfLine,endOfFile,fp,pos,max_pos,buf,max_is_eof);// This contains vertex number/texture/normal; we only want the vertex
          if (token != NULL) {
            str = strtok(token,"/");  // split the string at the first "/"
            my_vertices.push_back(atoi(str));
          }
        }

        if (my_vertices.size() > 3) ++f_split;
        nst += my_vertices.size()-2;  // assumes face can be built using a fan rather than centroid splitting; avoids node insertion

        //DEBUG
        // cout << "face:";
        // for (int tri=0,ntri=my_vertices.size(); tri<ntri; ++tri) {
        //   cout << " " << my_vertices[tri];
        // }
        // cout << endl;

        for (int tri=0,ntri=my_vertices.size()-2; tri<ntri; ++tri) {
          TriNodes thisTri;
          thisTri.nodes[0] = my_vertices[0];
          thisTri.nodes[1] = my_vertices[tri+1];
          thisTri.nodes[2] = my_vertices[tri+2];
          thisTri.zone = n_zones;
          TriNodesVec.push_back(thisTri);
        }
        my_vertices.clear();
      }
      else if (strcmp(firstTokenOfLine,"o") == 0 || strcmp(firstTokenOfLine,"g") == 0) {
        ++n_zones;
        char * token = getNextObjToken(endOfLine,endOfFile,fp,pos,max_pos,buf,max_is_eof);
        zoneVec.push_back(string(token));
      }

      if (!endOfLine) moveToStartOfNextLine(endOfFile,fp,pos,max_pos,buf,max_is_eof); // ignore unknown lines and possible empty space at end of line

    }
  }

  COUT2(" > imported vertex count: " << nsp);
  if (f_split) {COUT2(" > processed face count " << n_faces << "; required splitting " << f_split << " non-triangular faces");}
  else COUT2(" > imported tri count: " << nst);
  if (n_zones != -1) COUT2(" > imported " << (++n_zones) << " objects");

  // allocate new sizes for nsp/nst surface objects
  // allocate new sizes for nsp/nst surface objects
  growNspData(nsp,ss_nsp0);
  growNstData(nst,ss_nst0);

  int nv = 0;
  for (int iv = 0, nv_all = vertexVec.size(); iv < nv_all; ++iv) {
    xsp[nv+ss_nsp0][0] = vertexVec[iv].x;
    xsp[nv+ss_nsp0][1] = vertexVec[iv].y;
    xsp[nv+ss_nsp0][2] = vertexVec[iv].z;
    ++nv;
  }

  assert(nv+ss_nsp0 == nsp);
  vertexVec.clear();


  int newst = 0;
  for (int ist = 0, st_all = TriNodesVec.size(); ist < st_all; ++ist) {
    spost[ist+ss_nst0][0] = ss_nsp0 + TriNodesVec[ist].nodes[0] - 1;
    spost[ist+ss_nst0][1] = ss_nsp0 + TriNodesVec[ist].nodes[1] - 1;
    spost[ist+ss_nst0][2] = ss_nsp0 + TriNodesVec[ist].nodes[2] - 1;
    ++newst;
  }

  assert(newst+ss_nst0 == nst);

  if (n_zones == -1) {
    // no groups were in the file, so simply use filename for all zones
    const int new_zone = zoneVec.size();
    zoneVec.push_back(SurfaceZone(filename));

    for (int ist=ss_nst0; ist < nst; ++ist) {
      znost[ist] = new_zone;
    }
  }
  else {
    // use groups definedin obj file
    for (int ist = 0, st_all = TriNodesVec.size(); ist < st_all; ++ist) {
      znost[ist+ss_nst0] = TriNodesVec[ist].zone;
    }
  }
  TriNodesVec.clear();

  COUT1(" > initFromObj done");
  // ====================================================
  // make sure all zones were set...
  for (int ist = 0; ist < nst; ++ist) {
    assert((znost[ist] >= 0)&&(znost[ist] < int(zoneVec.size())));
  }

  return 0;  // success
}
