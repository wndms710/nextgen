#ifndef NEW_FLUENT_READER_HPP
#define NEW_FLUENT_READER_HPP

#include <math.h>
#include <assert.h>
#include "ByteSwap.hpp"
#include "MiscUtils.hpp"

using namespace std;

#define FLUENT_BUFFER_SIZE     65536
#define FLUENT_TOKEN_SIZE      65536
//#define FLUENT_BUFFER_SIZE     25
//#define FLUENT_TOKEN_SIZE      25

class NewFluentReader {

private:

  int pos, max_pos;
  char buf[FLUENT_BUFFER_SIZE+FLUENT_TOKEN_SIZE];

protected:

  FILE * fp;
  int level;
  bool byteSwap;

  int atox(char * token) {
    int i;
    sscanf(token, "%x", &i);
    return (i);
  }

  double atod(char * token) {
    double d;
    sscanf(token, "%lf", &d);
    return (d);
  }

public:

  NewFluentReader() {
    // reset tokenizing stuff...
    pos = max_pos = 0;
    level = 0;
    fp = NULL;
    byteSwap = false;
  }

  NewFluentReader(const string& filename) {
    init(filename);
  }

  ~NewFluentReader() {
    if (fp != NULL) fclose(fp);
  }

  void init(const string& filename) {
    // reset tokenizing stuff...
    pos = max_pos = 0;
    level = 0;
    byteSwap = false;
    const int file_err = MiscUtils::openFile(&fp,filename,"rb");
    if (file_err != 0) throw(file_err);
    // fp = fopen(filename.c_str(),"rb");
    // if (fp == NULL) {
    //   cerr << "Error: cannot open file: "<< filename << endl;
    //   throw(20);
    // }
  }

  void setByteSwap(const bool _byteSwap) {
    byteSwap = _byteSwap;
  }

  inline int getLevel() const {return level;}

  void rewindFile() {
    rewind(fp);
  }

  template<typename T>
  void readBinaryNodes(T * d, const int dim) {
    if (pos+int(sizeof(T))*dim > max_pos) {
      assert(pos > 0);
      const int n = max_pos-pos;
      if (n > 0)
        memmove(buf,buf+pos,n);
      max_pos = n + fread(buf+n, sizeof(char), FLUENT_BUFFER_SIZE, fp);
      pos = 0;
    }
    assert(pos+int(sizeof(T))*dim <= max_pos);
    memcpy(d,buf+pos,sizeof(T)*dim);
    pos += sizeof(T)*dim;
  }

  void readBinaryInt(int& d) {
    d = readBinary<int,int>();
  }

  void readBinaryFaces(int * d, const int nNodes) {

    // node list + two cells
    int nvars = nNodes + 2;

    for (int i=0; i<nvars; ++i) {
      d[i] = readBinary<int,int>();
    }
  }

  template<typename T_To, typename T_From>
  T_To readBinary() {
    T_From d;
    // we need to read sizeof(T_From) bytes from the current buf+pos. If the remaining
    // bytes are not sufficient to read, then realign before reading...

    if (pos+int(sizeof(T_From)) > max_pos) {
      assert(pos > 0);
      const int n = max_pos-pos;
      if (n > 0)
        memmove(buf,buf+(pos),n);

      max_pos = n + fread(buf+n, sizeof(char), FLUENT_BUFFER_SIZE, fp);
      pos = 0;
    }
    assert(pos+int(sizeof(T_From)) <= max_pos);
    memcpy(&d,buf+pos,sizeof(T_From));
    if (byteSwap) d = ByteSwap::byteSwap(d);
    pos += sizeof(T_From);

    return T_To(d);  // cast to the appropriate format
  };

  int getNextToken(char * &token) {

    token = NULL;
    int token_level = level;
    int quote_mode = 0;

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

        max_pos = pos + fread(buf+pos, sizeof(char), FLUENT_BUFFER_SIZE, fp);
        if (max_pos == pos) {
          buf[pos] = '\0';
          cerr << "Error at end of file: \"" << buf << "\"" << endl;
          assert(0);
        }

        /*
          cout << "buf is: \"";
          for (int i = 0; i < max_pos; ++i)
          cout << buf[i];
          cout << "\"" << endl;
          getchar();
        */
      }

      //cout << "getNextToken pos: " << pos << " c: \"" << buf[pos] << "\"" << endl;
      const char c = buf[pos++];

      if (c == '"') {
        if (quote_mode == 1)
          break;

        if (token != NULL) cout << "offending token that already exists: " << token << endl;
        assert(token == NULL);
        token = buf+pos; // first char after quote
        quote_mode = 1;
      }
      else if (quote_mode == 0) {
        if ( (c == ' ') || (c == '\t') || (c == '\n') || (c == 13) || (c == '\r') ) {
          // c is a whitespace character - this represents either the termination
          // of a token, or some space at the start of the next token...
          if (token)
            break;
        }
        else if (c == '(') {
          ++level;
          if (token)
            break;

          ++token_level;
        }
        else if (c == ')') {
          --level;
          if (token)
            break;

          --token_level;
        }
        else if (!token) {
          // any other character is the start of the token...
          token = buf+pos-1;
        }
      }
    }

    // terminate string and return...
    buf[pos-1] = '\0';
    return (token_level);
  }

  int getNextTokenAsInt(const int level_check) {
    char * token;
    const int level = getNextToken(token); assert(level == level_check);
    return atoi(token);
  }

  double getNextTokenAsDouble(const int level_check) {
    char * token;
    const int level = getNextToken(token); assert(level == level_check);
    return atod(token);
  }

  char * getNextTokenAsString(const int level_check) {
    char * token;
    const int level = getNextToken(token); assert(level == level_check);
    return token;
  }

  int getNextTokenAsHex(const int level_check) {
    char * token;
    const int level = getNextToken(token); assert(level == level_check);
    return atox(token);
  }

  bool advanceToLevel(const int target) {

    int quote_mode = 0;

    while (1) {

      if (pos >= max_pos) {
        max_pos = fread(buf, sizeof(char), FLUENT_BUFFER_SIZE, fp);
        //cout << "got max_pos: " << max_pos << endl;
        if (max_pos == 0)
          return false;

        pos = 0;
      }

      //cout << "advanceToLevel: " << target << " pos: " << pos << " c: \"" << buf[pos] << "\"" << endl;
      char c = buf[pos++];

      if (c == '"') {
        quote_mode = 1 - quote_mode;
      }
      else if (quote_mode == 0) {
        if (c == '(') {
          ++level;
          if (level == target)
            return true;
        }
        else if (c == ')') {
          --level;
        }
      }

    }

  }

  //void readMsh(


};

class FluentMsh {
public:

  vector<string> zoneVec;
  int nno,nfa,ncv;
  double (*x_no)[3];
  int *noofa_i,*noofa_v;
  int (*cvofa)[2];
  int *znofa;
  int *no_flag;
  int *fa_flag;
  
  FluentMsh() {
    x_no = NULL;
    noofa_i = NULL;
    noofa_v = NULL;
    cvofa = NULL;
    znofa = NULL;
    no_flag = NULL;
    fa_flag = NULL;
  }

  FluentMsh(const string& filename) {
    x_no = NULL;
    noofa_i = NULL;
    noofa_v = NULL;
    cvofa = NULL;
    znofa = NULL;
    no_flag = NULL;
    fa_flag = NULL;
    init(filename);
  }

  ~FluentMsh() {
    DELETE(x_no);
    DELETE(noofa_i);
    DELETE(noofa_v);
    DELETE(cvofa);
    DELETE(znofa);
    DELETE(no_flag);
    DELETE(fa_flag);
  }
  
  void writeSurfaceSbin(const string& filename) const {
    
    cout << "writeSurfaceSbin: \"" << filename << "\"" << endl;

    // boundary zones only...

    const int nzn = zoneVec.size();
    int * zone_count = new int[nzn];
    int * zone_new = new int[nzn];

    for (int izn = 0; izn < nzn; ++izn) {
      zone_count[izn] = 0;
    }

    int * no_flag = new int[nno];
    FOR_INO no_flag[ino] = -1;
    
    for (int ifa = 0; ifa < nfa; ++ifa) {
      if (cvofa[ifa][1] < 0) { 
        const int izn = znofa[ifa];
        assert((izn >= 0)&&(izn < nzn));
        // make sure this is a tri...
        const int nnof = noofa_i[ifa+1]-noofa_i[ifa];
        assert(nnof == 3); // otherwise, we need to split
        for (int nof = noofa_i[ifa]; nof != noofa_i[ifa+1]; ++nof) {
          const int ino = noofa_v[nof];
          assert((ino >= 0)&&(ino < nno));
          no_flag[ino] = 0;
        }
        ++zone_count[izn];
      }
    }
    
    int nzn_surface = 0;
    for (int izn = 0; izn < nzn; ++izn) {
      if (zone_count[izn] > 0) {
        zone_new[izn] = nzn_surface++;
      }
      else {
        assert(zone_count[izn] == 0);
        zone_new[izn] = -1;
      }
    }
    cout << " > surface nzn: " << nzn_surface << endl;

    int nno_surface = 0;
    for (int ino = 0; ino < nno; ++ino) {
      if (no_flag[ino] == 0) {
        no_flag[ino] = nno_surface++;
      }
      else {
        assert(no_flag[ino] == -1);
      }
    }
    cout << " > surface nno: " << nno_surface << endl;
    
    // turn zone_count into a displacement....
    int * zone_disp = new int[nzn];
    zone_disp[0] = 0;
    for (int izn = 1; izn < nzn; ++izn)
      zone_disp[izn] = zone_disp[izn-1] + zone_count[izn-1];
    const int nfa_surface = zone_disp[nzn-1] + zone_count[nzn-1];
    cout << " > surface nfa: " << nfa_surface << endl;
    
    int * fa_flag = new int[nfa_surface];
    for (int ifa = 0; ifa < nfa; ++ifa) {
      if (cvofa[ifa][1] < 0) { 
        const int izn = znofa[ifa];
        fa_flag[zone_disp[izn]++] = ifa;
      }
    }

    // rewind...
    zone_disp[0] = 0;
    for (int izn = 1; izn < nzn; ++izn)
      zone_disp[izn] = zone_disp[izn-1] + zone_count[izn-1];
    
    // now write the sbin file...
    FILE * fp = fopen(filename.c_str(),"wb");
    assert(fp != NULL);
    
    const int version = 1;
    cout << " > sbin version: " << version << endl;
    fwrite(&version,sizeof(int),1,fp);
    
    const int count = nzn_surface;
    fwrite(&count,sizeof(int),1,fp);
    
    cout << " > writing " << count << " zones:" << endl;
    for (int izn = 0; izn < nzn; ++izn) {
      if (zone_count[izn] > 0) {
        const int length = zoneVec[izn].length();
        fwrite(&length,sizeof(int),1,fp);
        fwrite(zoneVec[izn].c_str(),sizeof(char),length,fp);
        cout << "    > \"" << zoneVec[izn] << "\"" << endl;
      }
    }
      
    cout << " > nsp: " << nno_surface << endl;
    fwrite(&nno_surface,sizeof(int),1,fp);
    for (int ino = 0; ino < nno; ++ino) {
      if (no_flag[ino] >= 0) {
        fwrite(x_no[ino],sizeof(double),3,fp);
      }
    }
    
    cout << " > nst: " << nfa_surface << endl;
    fwrite(&nfa_surface,sizeof(int),1,fp);

    // spost first...
    int nof[3];
    for (int izn = 0; izn < nzn; ++izn) {
      if (zone_count[izn] > 0) {
        for (int ii = zone_disp[izn]; ii < zone_disp[izn]+zone_count[izn]; ++ii) {
          const int ifa = fa_flag[ii];
          assert(znofa[ifa] == izn);
          const int nnof = noofa_i[ifa+1]-noofa_i[ifa];
          assert(nnof == 3);
          nof[0] = no_flag[noofa_v[noofa_i[ifa]  ]];
          nof[1] = no_flag[noofa_v[noofa_i[ifa]+1]];
          nof[2] = no_flag[noofa_v[noofa_i[ifa]+2]];
          fwrite(nof,sizeof(int),3,fp);
        }
      }
    }

    // then the NEW zone index...
    for (int izn = 0; izn < nzn; ++izn) {
      if (zone_count[izn] > 0) {
        const int izn_surface = zone_new[izn];
        for (int ii = zone_disp[izn]; ii < zone_disp[izn]+zone_count[izn]; ++ii) {
          const int ifa = fa_flag[ii];
          assert(znofa[ifa] == izn);
          fwrite(&izn_surface,sizeof(int),1,fp);
        }
      }
    }

    fclose(fp);

    delete[] zone_count;
    delete[] zone_new;
    delete[] zone_disp;
    delete[] fa_flag;
    delete[] no_flag;
  
  }

  void writeTecplot(const string& filename) const {
    
    cout << "writeTecplot: \"" << filename << "\"" << endl;

    // boundary zones only...

    const int nzn = zoneVec.size();
    int * zone_count = new int[nzn];
    int * zone_nnof_max = new int[nzn];

    for (int izn = 0; izn < nzn; ++izn) {
      zone_count[izn] = 0;
      zone_nnof_max[izn] = 0;
    }

    for (int ifa = 0; ifa < nfa; ++ifa) {
      if (cvofa[ifa][1] < 0) { 
        const int izn = znofa[ifa];
        assert((izn >= 0)&&(izn < nzn));
        ++zone_count[izn];
        zone_nnof_max[izn] = max(zone_nnof_max[izn],noofa_i[ifa+1]-noofa_i[ifa]);
      }
    }

    for (int izn = 0; izn < nzn; ++izn) {
      cout << "zone " << izn << " zone_count: " << zone_count[izn] << endl;
    }

    // turn zone_count into a displacement....

    int * zone_disp = new int[nzn];
    zone_disp[0] = 0;
    for (int izn = 1; izn < nzn; ++izn)
      zone_disp[izn] = zone_disp[izn-1] + zone_count[izn-1];
    const int nfa_b = zone_disp[nzn-1] + zone_count[nzn-1];

    cout << " > nfa_b: " << nfa_b << endl;

    int * fa_flag = new int[nfa_b];
    for (int ifa = 0; ifa < nfa; ++ifa) {
      if (cvofa[ifa][1] < 0) { 
        const int izn = znofa[ifa];
        fa_flag[zone_disp[izn]++] = ifa;
      }
    }

    // rewind...
    zone_disp[0] = 0;
    for (int izn = 1; izn < nzn; ++izn)
      zone_disp[izn] = zone_disp[izn-1] + zone_count[izn-1];

    int * no_flag = new int[nno];
    FOR_INO no_flag[ino] = -1;

    FILE * fp = fopen(filename.c_str(),"w");
    fprintf(fp,"TITLE = \"%s\"\n",filename.c_str());
    fprintf(fp,"VARIABLES = \"X\"\n");
    fprintf(fp,"\"Y\"\n");
    fprintf(fp,"\"Z\"\n");

    MPI_Pause("about to do it!");

    for (int izn = 0; izn < nzn; ++izn) {
      
      cout << "working on zone: " << izn << " zone_count[izn]: " << zone_count[izn] << endl;
      
      if (zone_count[izn] > 0) {
        for (int ii = zone_disp[izn]; ii < zone_disp[izn]+zone_count[izn]; ++ii) {
          const int ifa = fa_flag[ii];
          assert(znofa[ifa] == izn);
          for (int nof = noofa_i[ifa]; nof != noofa_i[ifa+1]; ++nof) {
            const int ino = noofa_v[nof];
            assert((ino >= 0)&&(ino < nno));
            no_flag[ino] = -2;
          }
        }
        int nno_zone = 0;
        FOR_INO {
          if (no_flag[ino] == -2) {
            no_flag[ino] = nno_zone++;
          }
          else {
            no_flag[ino] = -1;
          }
        }
        if (zone_nnof_max[izn] == 3) {
          // triangle zone...
          fprintf(fp,"ZONE T=\"%s\"\n",zoneVec[izn].c_str());
          fprintf(fp,"N=%d, E=%d, F=FEPOINT, ET=TRIANGLE\n",nno_zone,zone_count[izn]);
          FOR_INO {
            if (no_flag[ino] >= 0) {
              fprintf(fp,"%18.15e %18.15e %18.15e\n",x_no[ino][0],x_no[ino][1],x_no[ino][2]);
            }
          }
          for (int ii = zone_disp[izn]; ii < zone_disp[izn]+zone_count[izn]; ++ii) {
            const int ifa = fa_flag[ii];
            assert(znofa[ifa] == izn);
            const int nof_f = noofa_i[ifa];
            assert(noofa_i[ifa+1]-nof_f == 3);
            fprintf(fp,"%d %d %d\n",no_flag[noofa_v[nof_f]]+1,no_flag[noofa_v[nof_f+1]]+1,no_flag[noofa_v[nof_f+2]]+1);
          }
        }
        else {
          // need to implement quads and polygons...
          assert(0);
        }
      }
    }

    fclose(fp);

    delete[] zone_count;
    delete[] zone_nnof_max;
    delete[] zone_disp;
    delete[] fa_flag;
    delete[] no_flag;

  }

  void buildTets(int& nte,int (*&noote)[4]) const {

    // start with collecting face counts and node counts at the cvs...

    int * faocv = new int[ncv];
    int * noofaocv = new int[ncv];
    FOR_ICV {
      faocv[icv] = 0;
      noofaocv[icv] = 0;
    }

    FOR_IFA {
      FOR_I2 {
        const int icv = cvofa[ifa][i];
        if (icv >= 0) {
          assert(icv < ncv);
          ++faocv[icv];
          noofaocv[icv] += noofa_i[ifa+1]-noofa_i[ifa];
        }
      }
    }

    nte = 0;
    FOR_ICV {
      if ((faocv[icv] == 4)&&(noofaocv[icv] == 12)) {
        // this is a single tet...
        // switch faocv to store the first tet, and
        // noofaocv to store the type...
        faocv[icv] = nte++;
        noofaocv[icv] = 1; // type 1: 1 tet
      }
      else {
        cout << "HOW TO DIVIDE INTO TETS: faocv[icv]: " << faocv[icv] << " noofaocv[icv]: " << noofaocv[icv] << endl;
        assert(0);
      }
    }

    assert(noote == NULL);
    noote = new int[nte][4];
    for (int ite = 0; ite < nte; ++ite)
      FOR_I4 noote[ite][i] = -1;

    // now push faces into the elements...

    FOR_IFA {
      const int nof_f = noofa_i[ifa];
      const int nnof = noofa_i[ifa+1]-nof_f;
      {
        const int icv0 = cvofa[ifa][0];
        assert(icv0 >= 0);
        // the starting (and only) tet index is in faocv...
        const int ite = faocv[icv0];
        // noofaocv stores the type...
        if (noofaocv[icv0] == 1) {
          // this is a tet...
          assert(nnof == 3);
          const int ino2 = noofa_v[nof_f];
          const int ino1 = noofa_v[nof_f+1];
          const int ino0 = noofa_v[nof_f+2];
          if ((noote[ite][0] == -1)&&(noote[ite][1] == -1)&&(noote[ite][2] == -1)) {
            noote[ite][0] = ino0;
            noote[ite][1] = ino1;
            noote[ite][2] = ino2;
          }
          else if ((noote[ite][1] == ino0)&&(noote[ite][0] == ino1)&&((noote[ite][3] == -1)||(noote[ite][3] == ino2))) {
            noote[ite][1] = ino0;
            noote[ite][0] = ino1;
            noote[ite][3] = ino2;
          }
          else if ((noote[ite][1] == ino1)&&(noote[ite][0] == ino2)&&((noote[ite][3] == -1)||(noote[ite][3] == ino0))) {
            noote[ite][1] = ino1;
            noote[ite][0] = ino2;
            noote[ite][3] = ino0;
          }
          else if ((noote[ite][1] == ino2)&&(noote[ite][0] == ino0)&&((noote[ite][3] == -1)||(noote[ite][3] == ino1))) {
            noote[ite][1] = ino2;
            noote[ite][0] = ino0;
            noote[ite][3] = ino1;
          }
          else if ((noote[ite][0] == ino0)&&(noote[ite][2] == ino1)&&((noote[ite][3] == -1)||(noote[ite][3] == ino2))) {
            noote[ite][0] = ino0;
            noote[ite][2] = ino1;
            noote[ite][3] = ino2;
          }
          else if ((noote[ite][0] == ino1)&&(noote[ite][2] == ino2)&&((noote[ite][3] == -1)||(noote[ite][3] == ino0))) {
            noote[ite][0] = ino1;
            noote[ite][2] = ino2;
            noote[ite][3] = ino0;
          }
          else if ((noote[ite][0] == ino2)&&(noote[ite][2] == ino0)&&((noote[ite][3] == -1)||(noote[ite][3] == ino1))) {
            noote[ite][0] = ino2;
            noote[ite][2] = ino0;
            noote[ite][3] = ino1;
          }
          else if ((noote[ite][2] == ino0)&&(noote[ite][1] == ino1)&&((noote[ite][3] == -1)||(noote[ite][3] == ino2))) {
            noote[ite][2] = ino0;
            noote[ite][1] = ino1;
            noote[ite][3] = ino2;
          }
          else if ((noote[ite][2] == ino1)&&(noote[ite][1] == ino2)&&((noote[ite][3] == -1)||(noote[ite][3] == ino0))) {
            noote[ite][2] = ino1;
            noote[ite][1] = ino2;
            noote[ite][3] = ino0;
          }
          else if ((noote[ite][2] == ino2)&&(noote[ite][1] == ino0)&&((noote[ite][3] == -1)||(noote[ite][3] == ino1))) {
            noote[ite][2] = ino2;
            noote[ite][1] = ino0;
            noote[ite][3] = ino1;
          }
          else {
            assert(0);
          }
        }
        else {
          assert(0);
        }
      }
      {
        const int icv1 = cvofa[ifa][1];
        if (icv1 >= 0) {
          // the starting (and only) tet index is in faocv...
          const int ite = faocv[icv1];
          // noofaocv stores the type...
          if (noofaocv[icv1] == 1) {
            // this is a tet...
            assert(nnof == 3);
            // reverse the node order for this cvofa 1...
            const int ino0 = noofa_v[nof_f];
            const int ino1 = noofa_v[nof_f+1];
            const int ino2 = noofa_v[nof_f+2];
            if ((noote[ite][0] == -1)&&(noote[ite][1] == -1)&&(noote[ite][2] == -1)) {
              noote[ite][0] = ino0;
              noote[ite][1] = ino1;
              noote[ite][2] = ino2;
            }
            else if ((noote[ite][1] == ino0)&&(noote[ite][0] == ino1)&&((noote[ite][3] == -1)||(noote[ite][3] == ino2))) {
              noote[ite][1] = ino0;
              noote[ite][0] = ino1;
              noote[ite][3] = ino2;
            }
            else if ((noote[ite][1] == ino1)&&(noote[ite][0] == ino2)&&((noote[ite][3] == -1)||(noote[ite][3] == ino0))) {
              noote[ite][1] = ino1;
              noote[ite][0] = ino2;
              noote[ite][3] = ino0;
            }
            else if ((noote[ite][1] == ino2)&&(noote[ite][0] == ino0)&&((noote[ite][3] == -1)||(noote[ite][3] == ino1))) {
              noote[ite][1] = ino2;
              noote[ite][0] = ino0;
              noote[ite][3] = ino1;
            }
            else if ((noote[ite][0] == ino0)&&(noote[ite][2] == ino1)&&((noote[ite][3] == -1)||(noote[ite][3] == ino2))) {
              noote[ite][0] = ino0;
              noote[ite][2] = ino1;
              noote[ite][3] = ino2;
            }
            else if ((noote[ite][0] == ino1)&&(noote[ite][2] == ino2)&&((noote[ite][3] == -1)||(noote[ite][3] == ino0))) {
              noote[ite][0] = ino1;
              noote[ite][2] = ino2;
              noote[ite][3] = ino0;
            }
            else if ((noote[ite][0] == ino2)&&(noote[ite][2] == ino0)&&((noote[ite][3] == -1)||(noote[ite][3] == ino1))) {
              noote[ite][0] = ino2;
              noote[ite][2] = ino0;
              noote[ite][3] = ino1;
            }
            else if ((noote[ite][2] == ino0)&&(noote[ite][1] == ino1)&&((noote[ite][3] == -1)||(noote[ite][3] == ino2))) {
              noote[ite][2] = ino0;
              noote[ite][1] = ino1;
              noote[ite][3] = ino2;
            }
            else if ((noote[ite][2] == ino1)&&(noote[ite][1] == ino2)&&((noote[ite][3] == -1)||(noote[ite][3] == ino0))) {
              noote[ite][2] = ino1;
              noote[ite][1] = ino2;
              noote[ite][3] = ino0;
            }
            else if ((noote[ite][2] == ino2)&&(noote[ite][1] == ino0)&&((noote[ite][3] == -1)||(noote[ite][3] == ino1))) {
              noote[ite][2] = ino2;
              noote[ite][1] = ino0;
              noote[ite][3] = ino1;
            }
            else {
              assert(0);
            }
          }
          else {
            assert(0);
          }
        }
      }
    }

    delete[] faocv;
    delete[] noofaocv;

    // make sure all tets are fully specified...
    for (int ite = 0; ite < nte; ++ite) {
      FOR_I4 {
        assert(noote[ite][i] >= 0);
      }
    }

  }

private:

  void init(const string& filename) {

    // here we use try/throw/catch error management beacause
    // we are calling from the constructor...

    NewFluentReader fr;
    fr.init(filename);

    nno = 0;
    nfa = 0;
    ncv = 0;
    assert(x_no == NULL);
    assert(noofa_i == NULL);
    assert(noofa_v == NULL);
    assert(cvofa == NULL);
    assert(znofa == NULL);
    assert(no_flag == NULL);
    assert(fa_flag == NULL);

    vector<int> noofa_v_counter;
    int * noofa_i_counter = NULL;

    int ifa_counter = 0;
    int icv_current = 0;

    int zn_header_count = 0;
    int no_zone_count = 0;
    int cv_zone_count = 0;
    int fa_zone_count = 0;

    cout << " > reading msh file \"" << filename << "\"..." << endl;

    while (fr.advanceToLevel(1)) {

      char * token;
      int level = fr.getNextToken(token);
      assert(level == 1);
      //cout << " > level 1 token \"" << token << "\"..." << endl;

      if (strcmp(token,"2")==0) {
        // format:
        // (2 ND)
        cout << "    > grid dimensions declaration:" << endl;
        const int nd = fr.getNextTokenAsInt(1);
        cout << "      - grid is " << nd << "-D" << endl;
        if (nd != 3) {
          cerr << "Error: only 3-D grids are supported" << endl;
          throw(-1);
        }
      }
      else if (strcmp(token,"4")==0) {
        // machine configuration; will tell you if little or big Endian
        const int endianFlag = fr.getNextTokenAsInt(2);
        cout << "    > endianFlag: " << endianFlag << endl;
        if (endianFlag != 60) {
          fr.setByteSwap(true);
          cout << "   > byte-swapping binary reads" << endl;
        }
        // for now don't bother with other flags since we don't know what they are....
      }
      else if (strcmp(token,"10")==0) {
        // format:
        // (10 (id start end type [ND]) (x y [z]))
        // [.] indicate optional parameters, depending on dimension
        const int id    = fr.getNextTokenAsHex(2);
        const int start = fr.getNextTokenAsHex(2);
        const int end   = fr.getNextTokenAsHex(2);
        const int type  = fr.getNextTokenAsHex(2);

        if (id == 0) {
          assert(start == 1);
          // declaration section, ND is omitted
          // set the node count for checking
          cout << "    > global node declaration: " << end << " nodes" << endl;
          // intialize new surface points data
          assert(nno == 0); nno = end;
          assert(x_no == NULL); x_no = new double[nno][3];
          assert(no_flag == NULL); no_flag = new int[nno];
          for (int ino = 0; ino < nno; ++ino) no_flag[ino] = 0;
        }
        else {
          // data present section
          const int noND  = fr.getNextTokenAsHex(2); // future: make ND param reading optional
          assert (noND == 3);

          cout << "    > nodal data block:" << endl;
          cout << "      - zone-id     : " << id << endl;
          cout << "      - start-index : " << start << endl;
          cout << "      - end-index   : " << end << endl;
          cout << "      - node-type   : " << type << endl;
          cout << "      - node-dim    : " << noND << endl;

          ++no_zone_count;
          fr.advanceToLevel(2);

          // all nodes in all zones must be read in for now...
          int ino_current = start-1;
          assert(nno >= end);
          assert(x_no != NULL);
          for (int i = start; i <= end; ++i) {
            assert(no_flag[ino_current] == 0);
            no_flag[ino_current] = 1;
            FOR_J3 x_no[ino_current][j] = fr.getNextTokenAsDouble(2);
            ++ino_current;
          }
        }
      }
      else if (strcmp(token,"12")==0) {
        // format:
        // (12 (id start end type [elem-type]) ())
        // [.] indicate optional parameters, depending on dimension
        const int id     = fr.getNextTokenAsHex(2);
        const int start  = fr.getNextTokenAsHex(2);
        const int end    = fr.getNextTokenAsHex(2);
        const int type   = fr.getNextTokenAsHex(2);

        if (id == 0) {
          // declaration section, ND is omitted
          if (end != 0) {
            cout << "    > global cell declaration: " << end << " cells" << endl;
            assert(start == 1);
            ncv = end;
          }
          else {
            cout << "    > global cell declaration: boundary mesh (no cell data)" << endl;
          }
        }
        else {
          // data present section
          const int elType   = fr.getNextTokenAsHex(2);
          cout << "    > cell data block:" << endl;
          cout << "      - zone-id     : " << id << endl;
          cout << "      - start-index : " << start << endl;
          cout << "      - end-index   : " << end << endl;
          cout << "      - cell-type   : " << type << endl;
          cout << "      - elem-type   : " << elType << endl;

          assert(icv_current == start-1);
          icv_current += end-start+1;
          ++cv_zone_count;
        }

      }
      else if (strcmp(token,"13")==0) {
        // format:
        // (13 (id start end bc-type fa-type) ())
        const int id     = fr.getNextTokenAsHex(2);
        const int start  = fr.getNextTokenAsHex(2);
        const int end    = fr.getNextTokenAsHex(2);
        const int type   = fr.getNextTokenAsHex(2);

        if (id == 0) {
          // declaration section, bc-type is omitted
          // set the global face count for checking
          assert(start == 1);
          cout << "    > file face declaration: " << end << " faces" << endl;

          // if we are reading everything...
          assert(nfa == 0);
          nfa = end;
          assert(ifa_counter == 0);
          assert(noofa_i == NULL); // leave this one for now
          assert(noofa_i_counter == NULL); noofa_i_counter = new int[nfa+1];
          noofa_i_counter[0] = 0;
          assert(noofa_v_counter.empty());
          assert(cvofa == NULL); cvofa = new int[nfa][2];
          assert(znofa == NULL); znofa = new int[nfa];
          assert(fa_flag == NULL); fa_flag = new int[nfa];
          for (int ifa = 0; ifa < nfa; ++ifa)
            fa_flag[ifa] = -1;
        }
        else {
          // data present section
          const int faType   = fr.getNextTokenAsHex(2);

          cout << "    > face data block:" << endl;
          cout << "      - zone-id     : " << id << endl;
          cout << "      - start-index : " << start << endl;
          cout << "      - end-index   : " << end << endl;
          cout << "      - bc-type     : " << type << endl;
          cout << "      - fa-type     : " << faType << endl;

          fr.advanceToLevel(2);

          // create unique Zone identifier
          const int izone = zoneVec.size();
          ++fa_zone_count;
          ostringstream oss;
          oss << "X" << id;
          zoneVec.push_back(oss.str());

          if ((faType >= 2)&&(faType <=4)) {
            // list all nodes, then two cells
            cout << "      > no_per_fa: " << faType << endl;
            int ifa_current = start-1;
            assert(nfa >= end);
            for (int ifa = start; ifa <= end; ++ifa) {
              assert(fa_flag[ifa_current] == -1);
              fa_flag[ifa_current] = ifa_counter++;
              for (int i = 0; i < faType; ++i) {
                const int ino = fr.getNextTokenAsHex(2)-1;
                assert((ino >= 0)&&(ino < nno));
                noofa_v_counter.push_back(ino);
              }
              noofa_i_counter[ifa_counter] = noofa_v_counter.size();
              znofa[ifa_current] = izone;
              cvofa[ifa_current][0] = fr.getNextTokenAsHex(2)-1;
              cvofa[ifa_current][1] = fr.getNextTokenAsHex(2)-1;
              ++ifa_current;
            }
          }
          else if (faType == 0 || faType == 5) {
            // list all nodes, then two cells
            int ifa_current = start-1;
            assert(nfa >= end);
            for (int ifa = start; ifa <= end; ++ifa) {
              assert(fa_flag[ifa_current] == -1);
              fa_flag[ifa_current] = ifa_counter++;
              const int no_per_fa = fr.getNextTokenAsHex(2);
              for (int i = 0; i < no_per_fa; ++i) {
                const int ino = fr.getNextTokenAsHex(2)-1;
                assert((ino >= 0)&&(ino < nno));
                noofa_v_counter.push_back(ino);
              }
              noofa_i_counter[ifa_counter] = noofa_v_counter.size();
              znofa[ifa_current] = izone;
              cvofa[ifa_current][0] = fr.getNextTokenAsHex(2)-1;
              cvofa[ifa_current][1] = fr.getNextTokenAsHex(2)-1;
              ++ifa_current;
            }
          }
          else {
            cout << "Unsupported boundary face type detected (faType = " << faType << ")" << endl;
          }
        }
      }
      else if (strcmp(token,"3010")==0 || strcmp(token,"2010")==0) {
        // binary node block...
        //const int id    = fr.getNextTokenAsHex(2);
        //const int start = fr.getNextTokenAsHex(2);
        //const int end   = fr.getNextTokenAsHex(2);
        //const int type  = fr.getNextTokenAsHex(2);
        //const int noND  = fr.getNextTokenAsHex(2);
        fr.getNextTokenAsHex(2);
        fr.getNextTokenAsHex(2);
        fr.getNextTokenAsHex(2);
        fr.getNextTokenAsHex(2);
        fr.getNextTokenAsHex(2);

        assert(0);
        /*

          nno_count += end - start + 1;
          ++no_zone_count;

          cout << "    > nodal data block:" << endl;
          cout << "      - zone-id     : " << id << endl;
          cout << "      - start-index : " << start << endl;
          cout << "      - end-index   : " << end << endl;
          cout << "      - node-type   : " << type << endl;
          cout << "      - node-dim    : " << noND << endl;

          fr.advanceToLevel(2);

          // ------------------------------------------
          // set surface node stuff...
          assert(noND == 3);
          assert((start == 1) && (end == nno)); // ensures there is one node zone -- support multiple in the future
          assert(no_zone_count == 1);
          assert((type == 1) || (type == 2)); // regular or boundary nodes

          if (strcmp(token,"3010")==0) {
          if (x_no.size() < end) x_no.resize(end);
          double x[3];
          for (int i = start; i <= end; ++i) {
          fr.readBinaryNodes<double>(x,3);
          // nodes have 0 indexing, Fluent uses 1-indexing
          x_no[i-1].first = x[0];
          x_no[i-1].second = x[1];
          x_no[i-1].third = x[2];
          }
          }
          else {
          assert(strcmp(token,"2010")==0);
          if (x_no.size() < end) x_no.resize(end);
          float x[3];
          for (int i = start; i <= end; ++i) {
          fr.readBinaryNodes<float>(x,3);
          // nodes have 0 indexing, Fluent uses 1-indexing
          x_no[i-1].first = (double)x[0];
          x_no[i-1].second = (double)x[1];
          x_no[i-1].third = (double)x[2];
          }
          }
          // ------------------------------------------
          */

      }
      else if (strcmp(token,"2012")==0 || strcmp(token,"3012")==0) {
        // Need to properly stride across binary sections
        // (2012 (id start end bc-type fa-type) ())
        const int id     = fr.getNextTokenAsHex(2);
        const int start  = fr.getNextTokenAsHex(2);
        const int end    = fr.getNextTokenAsHex(2);
        const int type   = fr.getNextTokenAsHex(2);
        const int elemType = fr.getNextTokenAsHex(2);
        UNUSED(id);
        UNUSED(type);

        // only need to read cell description if of mixed-type
        if (elemType == 0) {
          fr.advanceToLevel(2);

          int cellType;
          for (int i = start; i <= end; ++i) {
            fr.readBinaryInt(cellType); // TODO: just advance if not using data...
          }
        }
        // otherwise cell block is empty

        // ------------------------------------------

      }
      else if (strcmp(token,"2013")==0 || strcmp(token,"3013")==0) {
        // binary faces block...
        // (2013 (id start last bc-type face-type) no n1 n2 c0 c1)
        const int id     = fr.getNextTokenAsHex(2);
        const int start  = fr.getNextTokenAsHex(2);
        const int end    = fr.getNextTokenAsHex(2);
        const int bcType = fr.getNextTokenAsHex(2);
        const int faType = fr.getNextTokenAsHex(2);

        cout << "    > face data block:" << endl;
        cout << "      - zone-id     : " << id << endl;
        cout << "      - start-index : " << start << endl;
        cout << "      - end-index   : " << end << endl;
        cout << "      - bc-type     : " << bcType << endl;
        cout << "      - fa-type     : " << faType << endl;

        fr.advanceToLevel(2);

        assert(0);

        /*

          bool skipZone = false;
          set<int>::iterator it;
          it = skipFaceZoneIds.find(id);
          if (it!=skipFaceZoneIds.end()) {
          // this zone should be skipped b/c it is internal
          skipZone = true;
          }

          // ----------------------------------------------
          // update surface
          const int izone = zoneVec.size();
          if (!skipZone) {
          nfa_count += end - start + 1;
          ++fazone_count;
          ostringstream oss;
          oss << "X" << id;
          zoneVec.push_back(oss.str());
          }
          // ----------------------------------------------

          if (faType >= 2 && faType <= 4) {
          int x[faType+2];
          for (int ifa = start; ifa <= end; ++ifa) {
          fr.readBinaryFaces(x,faType);

          for (int i = 0; i < faType; ++i) {
          assert((x[i] >= 1) && (x[i] <= nno));
          }
          // for (int i = faType; i < (faType+2); ++i) {
          //   assert(x[i] == 0); // true if boundary mesh
          // }

          if (!skipZone) {
          // create the tris for this polygon
          for (int ist = 0; ist < faType-2; ++ist) {
          spost[current_tri][0] = x[0]-1;     // faces and nodes are 1-indexed
          spost[current_tri][1] = x[ist+1]-1; // faces and nodes are 1-indexed
          spost[current_tri][2] = x[ist+2]-1; // faces and nodes are 1-indexed

          assert(znost[current_tri] == -1);
          znost[current_tri] = izone;
          ++current_tri;
          }
          }
          }
          }
          else if (faType == 0 || faType == 5) {
          // first variable declares number of points per face
          // still only support tris
          int no_per_fa;
          for (int i = start; i <= end; ++i) {
          fr.readBinaryInt(no_per_fa);

          int x[no_per_fa+2];
          fr.readBinaryFaces(x,no_per_fa);
          for (int i = 0; i < no_per_fa; ++i) {
          assert((x[i] >= 1) && (x[i] <= nno));
          }
          // for (int i = no_per_fa; i < (no_per_fa+2); ++i) {
          //   assert(x[i] == 0); // true if boundary mesh
          // }

          if (!skipZone) {
          // create the tris for this polygon
          for (int ist = 0; ist < no_per_fa-2; ++ist) {
          spost[current_tri][0] = x[0]-1;     // faces and nodes are 1-indexed
          spost[current_tri][1] = x[ist+1]-1; // faces and nodes are 1-indexed
          spost[current_tri][2] = x[ist+2]-1; // faces and nodes are 1-indexed

          assert(znost[current_tri] == -1);
          znost[current_tri] = izone;
          ++current_tri;
          }
          }
          }
          }
        */

      }
      else if (strcmp(token,"2041")==0) {
        // Need to properly stride across binary sections
        // (2041 (start end) ())
        const int start  = fr.getNextTokenAsHex(2);
        const int end    = fr.getNextTokenAsHex(2);

        fr.advanceToLevel(2);

        // assumes an int flag
        int tmp_flag;
        for (int i = start; i <= end; ++i) {
          fr.readBinaryInt(tmp_flag);
        }
      }
      else if (strcmp(token,"2059")==0 || strcmp(token,"2058")==0) {
        // Need to properly stride across binary sections
        // (2059/58 (start end parent-zone child-zone) ())
        // where zones are either face or cell blocks, depending
        const int start = fr.getNextTokenAsHex(2);
        const int end = fr.getNextTokenAsHex(2);
        const int parentZone = fr.getNextTokenAsHex(2);
        const int childZone = fr.getNextTokenAsHex(2);
        UNUSED(parentZone);
        UNUSED(childZone);

        fr.advanceToLevel(2);

        int n_children;
        for (int i = start; i <= end; ++i) {
          fr.readBinaryInt(n_children);
          for (int j=0; j<n_children; ++j) fr.readBinary<int,int>();  // don't need to store this info
        }
      }
      else if ((strcmp(token,"45")==0) || (strcmp(token,"39")==0))  {
        // format:
        // (45 (id zone-type zone-name) ())
        const int id     = fr.getNextTokenAsInt(2);
        const char * znType = fr.getNextTokenAsString(2);
        const char * znName = fr.getNextTokenAsString(2);

        cout << "    > zone header block:" << endl;
        cout << "      - zone-id   : " << id << endl;
        cout << "      - zone-type : " << znType << endl;
        cout << "      - zone-name : " << znName << endl;

        ostringstream oss;
        oss << "X" << id;
        const string oss_str = oss.str();

        // look for this zone in zoneVec...
        int i;
        for (i = 0; i < int(zoneVec.size()); ++i) {
          if (zoneVec[i] == oss_str) {
            zoneVec[i] = znName;
            break;
          }
        }
        // make sure we found it...
        //assert(i < this->zoneVec.size()); # assumes sequential ordering from 1...

        if (i< int(zoneVec.size())) ++zn_header_count;

      }
      else {

        cout << "    > skipping unrecognized block token \"" << token << "\"" << endl;

      }
    }
    
    // check that al nodes were read...
    for (int ino = 0; ino < nno; ++ino)
      assert(no_flag[ino] == 1);
    
    // at this point, the fa_flag contains the index to the
    // faces as read in file order. We need to reorder the nodes
    // and populate noofa_v properly...
    assert(ifa_counter == nfa);
    assert(noofa_i == NULL);
    noofa_i = new int[nfa+1];
    noofa_i[0] = 0;
    assert(noofa_v == NULL);
    noofa_v = new int[noofa_v_counter.size()];
    // loop through in face-order, pulling nodes from the file order...
    int nof = 0; 
    for (int ifa = 0; ifa < nfa; ++ifa) {
      ifa_counter = fa_flag[ifa];
      assert((ifa_counter >= 0)&&(ifa_counter < nfa));
      for (int nof_counter = noofa_i_counter[ifa_counter]; nof_counter != noofa_i_counter[ifa_counter+1]; ++nof_counter) {
        noofa_v[nof++] = noofa_v_counter[nof_counter];
      }
      noofa_i[ifa+1] = nof;
    }
    assert(nof == noofa_v_counter.size());
    delete[] noofa_i_counter;
    noofa_v_counter.clear();
    
    assert(icv_current == ncv);

    // Count Checking
    cout << " > initFromMsh Summary:" << endl;
    cout << "     " << nno << " nodes" << endl;
    cout << "     " << nfa << " faces" << endl;
    cout << "     " << ncv << " cvs" << endl;
    cout << "     " << no_zone_count << " node zones" << endl;
    cout << "     " << fa_zone_count << " face zones" << endl;
    cout << "     " << cv_zone_count << " cell zones" << endl;
    cout << "     " << zn_header_count << " zone headers" << endl;

    for (int izone = 0,nzn=zoneVec.size(); izone < nzn; ++izone) cout << " zone: " << izone << ": \"" << zoneVec[izone] << "\"" << endl;
    
    cout << "average nodes per face: " << double(nof)/double(nfa) << endl;

    // ====================================================
    // make sure all zones were set...
    for (int ifa = 0; ifa < nfa; ++ifa) {
      assert((znofa[ifa] >= 0)&&(znofa[ifa] < int(zoneVec.size())));
    }

    // flip all faces because Fluent convention is opposite of ours...
    bool b_flip_cvofa = false;
    bool b_do_not_flip_cvofa = false;
    for (int ifa = 0; ifa < nfa; ++ifa) {
      const int nnof = noofa_i[ifa+1] - noofa_i[ifa];
      for (int i = 0; i < nnof/2; ++i) {
        const int tmp = noofa_v[noofa_i[ifa]+i];
        noofa_v[noofa_i[ifa]+i] = noofa_v[noofa_i[ifa+1]-1-i];
        noofa_v[noofa_i[ifa+1]-1-i] = tmp;
      }
      // also flip cvofa if ordering is wrong...
      // this is the case in one of the ennova files...
      if (cvofa[ifa][1] < 0) {
        assert(cvofa[ifa][1] == -1);
        assert(cvofa[ifa][0] >= 0);
        b_do_not_flip_cvofa = true;
      }
      if (cvofa[ifa][0] < 0) {
        assert(cvofa[ifa][0] == -1);
        assert(cvofa[ifa][1] >= 0);
        b_flip_cvofa = true;
      }
    }
    
    if (b_do_not_flip_cvofa && b_flip_cvofa) {
      CERR("cvofa has mixed orientation. What to do?");
    }
    
    if (b_flip_cvofa) {
      cout << " > WARNING: flipping cvofa" << endl;
      for (int ifa = 0; ifa < nfa; ++ifa) {
        const int tmp = cvofa[ifa][0];
        cvofa[ifa][0] = cvofa[ifa][1];
        cvofa[ifa][1] = tmp;
      }
    }
    
    // finally test that all nodes are touched by faces...
    FOR_INO no_flag[ino] = -1;
    for (int nof = 0; nof < noofa_i[nfa]; ++nof) {
      const int ino = noofa_v[nof];
      assert((ino >= 0)&&(ino < nno));
      no_flag[ino] = 0;
    }

    int n_unused=0;
    int count = 0;
    FOR_INO {
      if (no_flag[ino] != 0) {
        ++n_unused;
        COUT1(" ! free floating node at: " << COUT_VEC(x_no[ino]));
      }
      else no_flag[ino] = count++;  // index map
     // assert(no_flag[ino] == 0);
     //
    }

    if (n_unused) {
    CWARN(n_unused << " nodes not matched to face; pruning (may indicate missing tet(s))");

      double (*tmp_no)[3] = x_no; x_no = NULL;
      x_no = new double [nno-n_unused][3];

      FOR_INO {
        if (no_flag[ino] >= 0) {
          FOR_I3 x_no[no_flag[ino]][i] = tmp_no[ino][i];
        }
      }

      for (int nof = 0; nof < noofa_i[nfa]; ++nof) {
        const int ino = noofa_v[nof];
        assert((ino >= 0)&&(ino < nno));
        noofa_v[nof] = no_flag[ino];
      }

      nno -= n_unused;  // prune
      delete[] tmp_no; tmp_no = NULL;
    }

    delete[] no_flag; no_flag = NULL;
  }
  
};

#endif
