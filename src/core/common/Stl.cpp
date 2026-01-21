#include "Stl.hpp"
#include "Params.hpp"
#include <stack>

void Stl::initFromPlot3DFile(const char * filename, int istart, int iend, int jstart, int jend, int kstart, int kend) {

  FILE * fp;
  if ( (fp = fopen(filename,"r")) == NULL ) {
    cerr << "Error: plot3d file " << filename << " cannot be found/opened." << endl;
    throw(-1);
  }


  // assuming to extract one surface from a single block grid....
  int ni,nj,nk;
  fscanf(fp, "%i %i %i \n ", &ni, &nj, &nk );

  /*
    int istart = 1;
    int iend = ni;
    int jstart = 1;
    int jend = nj;
    int kstart = 63;
    int kend = 63;
  */

  // number of facets (per each quad on the structured grid -> 2 facets)
  if (istart==iend) {
    nf = 2*(jend-jstart)*(kend-kstart);
  }
  if (jstart==jend) {
    nf = 2*(iend-istart)*(kend-kstart);
  }
  if (kstart==kend) {
    nf = 2*(iend-istart)*(jend-jstart);
  }

#ifdef DEBUG
  cout << " reading Plot3d file " << filename << endl;
  cout << " block dimensions " << ni << "x" << nj << "x" << nk << endl;
  cout << " extracting surface ("<<istart<<":"<<iend<<"),("<<jstart<<":"<<jend<<"),("<<kstart<<":"<<kend<<")" << endl;
  cout << " generating an STL surface with "<< nf << " facets " << endl;
#endif

  // allocate temporary array for the entire grid
  float (*xtemp)[3] = new float[ni*nj*nk][3];
  float xmax[3],xmin[3];

  int id;
  for (id=0;id<3;id++) {
    xmax[id] = -10000.;
    xmin[id] =  10000.;
  }
  int i,j,k,ijk;
  for (id=0;id<3;id++) {
    for (k=0;k<nk;k++) {
      for (j=0;j<nj;j++) {
        for (i=0;i<ni;i++) {
          ijk = i+j*(ni)+k*(ni)*(nj);
          fscanf(fp, "%f", &xtemp[ijk][id]);
          //                cout << "read xtemp["<<ijk<<"]["<<id<<"]= " << xtemp[ijk][id] << endl;
          xmax[id] = max(xmax[id],xtemp[ijk][id]);
          xmin[id] = min(xmin[id],xtemp[ijk][id]);
        }
      }
    }
  }
  //      for (ijk=0;ijk<ni*nj*nk;ijk++) {
  //        cout << " x["<<ijk<<"]=("<<xtemp[ijk][0]<<","<< xtemp[ijk][1]<<","<< xtemp[ijk][2]<<")"<< endl;
  //      }
#ifdef DEBUG
  cout << " bouding box (" << xmin[0] << ":" << xmax[0] << "),("<< xmin[1] << ":" << xmax[1] << "),(" << xmin[2] << ":" << xmax[2] << ")" << endl;

  cout << " finished reading the Plot3D file ... " << endl;
#endif
  fclose(fp);

  // for now assume 3 unique vertices per facet...
  nv = 3*nf;

  // allocate...
  normal = new float[nf][3];
  v_of_f = new int[nf][3];
  xv     = new float[nv][3];

  // now construct the in a j plane facets....

  int ifacet = 0;
  int iv = 0;


  if (istart==iend) {
    i = istart-1;
    for (k=kstart-1;k<kend-1;k++) {
      for (j=jstart-1;j<jend-1;j++) {
        // first facet
        // first vertex
        ijk = (i  )+(j  )*(ni)+(k  )*(ni)*(nj);
        for (id=0;id<3;id++) { xv[iv][id] = xtemp[ijk][id]; }
        v_of_f[ifacet][0] = iv; iv++;
        // second vertex
        ijk = (i  )+(j+1)*(ni)+(k  )*(ni)*(nj);
        for (id=0;id<3;id++) { xv[iv][id] = xtemp[ijk][id]; }
        v_of_f[ifacet][1] = iv; iv++;
        // third vertex
        ijk = (i  )+(j+1)*(ni)+(k+1)*(ni)*(nj);
        for (id=0;id<3;id++) { xv[iv][id] = xtemp[ijk][id]; }
        v_of_f[ifacet][2] = iv; iv++; ifacet++;
        // second facet
        // first vertex
        ijk = (i  )+(j  )*(ni)+(k  )*(ni)*(nj);
        for (id=0;id<3;id++) { xv[iv][id] = xtemp[ijk][id]; }
        v_of_f[ifacet][0] = iv; iv++;
        // second vertex
        ijk = (i  )+(j+1)*(ni)+(k+1)*(ni)*(nj);
        for (id=0;id<3;id++) { xv[iv][id] = xtemp[ijk][id]; }
        v_of_f[ifacet][1] = iv; iv++;
        // third vertex
        ijk = (i  )+(j  )*(ni)+(k+1)*(ni)*(nj);
        for (id=0;id<3;id++) { xv[iv][id] = xtemp[ijk][id]; }
        v_of_f[ifacet][2] = iv; iv++; ifacet++;
      }
    }
  }

  if (jstart==jend) {
    j = jstart-1;
    for (k=kstart-1;k<kend-1;k++) {
      for (i=istart-1;i<iend-1;i++) {
        // first facet
        // first vertex
        ijk = (i  )+j*(ni)+(k  )*(ni)*(nj);
        for (id=0;id<3;id++) { xv[iv][id] = xtemp[ijk][id]; }
        v_of_f[ifacet][0] = iv; iv++;
        // second vertex
        ijk = (i+1)+j*(ni)+(k  )*(ni)*(nj);
        for (id=0;id<3;id++) { xv[iv][id] = xtemp[ijk][id]; }
        v_of_f[ifacet][1] = iv; iv++;
        // third vertex
        ijk = (i+1)+j*(ni)+(k+1)*(ni)*(nj);
        for (id=0;id<3;id++) { xv[iv][id] = xtemp[ijk][id]; }
        v_of_f[ifacet][2] = iv; iv++; ifacet++;
        // second facet
        // first vertex
        ijk = (i  )+j*(ni)+(k  )*(ni)*(nj);
        for (id=0;id<3;id++) { xv[iv][id] = xtemp[ijk][id]; }
        v_of_f[ifacet][0] = iv; iv++;
        // second vertex
        ijk = (i+1)+j*(ni)+(k+1)*(ni)*(nj);
        for (id=0;id<3;id++) { xv[iv][id] = xtemp[ijk][id]; }
        v_of_f[ifacet][1] = iv; iv++;
        // third vertex
        ijk = (i  )+j*(ni)+(k+1)*(ni)*(nj);
        for (id=0;id<3;id++) { xv[iv][id] = xtemp[ijk][id]; }
        v_of_f[ifacet][2] = iv; iv++; ifacet++;
      }
    }
  }

  if (kstart==kend) {
    k = kstart-1;
    for (j=jstart-1;j<jend-1;j++) {
      for (i=istart-1;i<iend-1;i++) {
        // first facet
        // first vertex
        ijk = (i  )+(j  )*(ni)+(k  )*(ni)*(nj);
        for (id=0;id<3;id++) { xv[iv][id] = xtemp[ijk][id]; }
        v_of_f[ifacet][0] = iv; iv++;
        // second vertex
        ijk = (i+1)+(j  )*(ni)+(k  )*(ni)*(nj);
        for (id=0;id<3;id++) { xv[iv][id] = xtemp[ijk][id]; }
        v_of_f[ifacet][1] = iv; iv++;
        // third vertex
        ijk = (i+1)+(j+1)*(ni)+(k  )*(ni)*(nj);
        for (id=0;id<3;id++) { xv[iv][id] = xtemp[ijk][id]; }
        v_of_f[ifacet][2] = iv; iv++; ifacet++;
        // second facet
        // first vertex
        ijk = (i  )+(j  )*(ni)+(k  )*(ni)*(nj);
        for (id=0;id<3;id++) { xv[iv][id] = xtemp[ijk][id]; }
        v_of_f[ifacet][0] = iv; iv++;
        // second vertex
        ijk = (i+1)+(j+1)*(ni)+(k  )*(ni)*(nj);
        for (id=0;id<3;id++) { xv[iv][id] = xtemp[ijk][id]; }
        v_of_f[ifacet][1] = iv; iv++;
        // third vertex
        ijk = (i  )+(j+1)*(ni)+(k  )*(ni)*(nj);
        for (id=0;id<3;id++) { xv[iv][id] = xtemp[ijk][id]; }
        v_of_f[ifacet][2] = iv; iv++; ifacet++;
      }
    }
  }



  //cout << iv << " " << nv << " " << ifacet << " " << nf << endl;

  // checks...
  if (iv != nv) {
    cerr << "Error: problem in stl file 3." << endl;
    throw(-1);
  }

  if (ifacet != nf) {
    cerr << "Error: problem in stl file 4." << endl;
    throw(-1);
  }

  // get a "fake" normal

  for (ifacet=0;ifacet<nf;ifacet++) {
    for (id=0;id<3;id++) {
      normal[ifacet][id] = 1.;
    }
  }

  // cleanup

  delete[] xtemp;

}

void Stl::asciiCountObjects(const char * filename) {
  cout << " > opening stl file " << filename << endl;
  FILE * fp;
  if ( (fp = fopen(filename,"r")) == NULL ) {
    cerr << "    > ERROR: stl file " << filename << " cannot be found/opened." << endl;
    throw(-1);
  }
  asciiCountObjects(fp);

  fclose(fp);
}

void Stl::asciiCountObjects(FILE * fp) {
  cout << "    > counting elements" << endl;
  rewind(fp); // ensure we are at top of file
  // count the number of facets, solids
  nf = 0;
  ns = 0;
  char line[128], s[16];

  while (fgets(line,128,fp) != NULL) {
    if (sscanf(line,"%s",s) == 1) {
      if (strcmp(s,"facet") == 0) {
        nf++;
      }
      else if (strcmp(s,"solid") == 0) {
        ns++;
      }
    }
  }
  cout << "      solids: " << ns << endl;
  cout << "      facets: " << nf << endl;
  nv = nf*3; //assume this for now
}

void Stl::asciiReadObjects(const char * filename) {
        cout << " > opening stl file " << filename << endl;
        FILE * fp;
        if ( (fp = fopen(filename,"r")) == NULL ) {
          cerr << "    > ERROR: stl file " << filename << " cannot be found/opened." << endl;
          throw(-1);
        }

        asciiReadObjects(fp);

        fclose(fp);
}

void Stl::asciiReadObjects(FILE * fp) {
  cout << "    > reading in data" << endl;
  int ifacet = 0;
  int isolid = 0;
  int iv = 3;
  rewind(fp); // ensure read from top of file

  int nlin = 0;
  char line[128], s[16];
        while (fgets(line,128,fp) != NULL) {
                nlin++;
                if (sscanf(line,"%s",s) == 1) {
                        if (strcmp(s,"facet") == 0) {
                                if (iv != 3) {
                                        cerr << "    > ERROR: problem with number of vertices listed on line: " << nlin << endl;
                                        throw(-1);
                                }
                                sscanf(line,"%*s %*s %f %f %f\n",&(normal[ifacet][0]),
                                &(normal[ifacet][1]),&(normal[ifacet][2]));

                                sof[ifacet] = isolid-1; // maintain zero-indexing

                                ifacet++;
                                iv = 0;

                        }
                        else if (strcmp(s,"vertex") == 0) {
                                if (iv >= 3) {
                                        cerr << "    > ERROR: problem with number of vertices listed on line: " << nlin << endl;
                                        throw(-1);
                                }
                                sscanf(line,"%*s %f %f %f\n",&(xv[3*(ifacet-1)+iv][0]),
                                &(xv[3*(ifacet-1)+iv][1]),&(xv[3*(ifacet-1)+iv][2]));
                                v_of_f[ifacet-1][iv] = 3*(ifacet-1) + iv;
                                iv++;
                        }
                        else if (strcmp(s,"solid") == 0) {
                                char sname[128];
                                sscanf(line,"%*s %s\n",sname);
                                solidName[isolid] = sname;
                                isolid++;
                        }
                }
        }

        // checks...
        if (iv != 3) {
                cerr << "    > ERROR: problem in number of vertices" << endl;
                throw(-1);
        }

        if (ifacet != nf) {
                cerr << "    > ERROR: mismatch between facets counted/read: " << nf << "/" << ifacet << endl;
                throw(-1);
        }

        if (isolid != ns) {
                cerr << "    > ERROR: mismatch between solids counted/read: " << ns << "/" << isolid << endl;
                throw(-1);
        }
}

/// stuff for Binary STLs


int Stl::check_stl_file(const char * filename) {

  int stl_type=-1;
  unsigned char chtest[128];
  unsigned int i;
  FILE * fp;
  //long file_size;

  /* Open the file */
  fp = fopen(filename, "r");
  if(fp == NULL) {
    cerr << " > Couldn't open/find " << filename <<  "  for reading\n      verify path and filetype" << endl;
    throw(0);
  }
  // Find size of file
  fseek(fp, 0, SEEK_END);
  //file_size = ftell(fp);

  // Check for binary or ASCII file
  fseek(fp, HEADER_SIZE, SEEK_SET);
  fread(chtest, sizeof(chtest), 1, fp);
  stl_type = 0;
  for(i = 0; i < sizeof(chtest); i++)
    {
      if(chtest[i] > 127)
        {
          stl_type = 1;
          break;
        }
    }
  fclose(fp);

  return stl_type;

}


float Stl::stl_get_little_float(FILE *fp)
{
  union
  {
    int   int_value;
    float float_value;
  } value;

  value.int_value  =  fgetc(fp) & 0xFF;
  value.int_value |= (fgetc(fp) & 0xFF) << 0x08;
  value.int_value |= (fgetc(fp) & 0xFF) << 0x10;
  value.int_value |= (fgetc(fp) & 0xFF) << 0x18;
  return(value.float_value);
}

int Stl::stl_get_little_int(FILE *fp)
{

  int value;
  value  =  fgetc(fp) & 0xFF;
  value |= (fgetc(fp) & 0xFF) << 0x08;
  value |= (fgetc(fp) & 0xFF) << 0x10;
  value |= (fgetc(fp) & 0xFF) << 0x18;
  return(value);
}

int Stl::stl_get_number_of_facet(const char * filename)
{

  long           file_size;
  //int            header_num_facets;
  int            num_facets;
  FILE * fp;
  char this_header[81];

  /* Open the file */
  fp = fopen(filename, "rb");
  if(fp == NULL) {
    cerr << "Couldn't open " << filename <<  "  for reading" << endl;
    throw(0);
  }
  // Find size of file
  fseek(fp, 0, SEEK_END);
  file_size = ftell(fp);

  /* Test if the STL file has the right size  */
  if(((file_size - HEADER_SIZE) % SIZEOF_STL_FACET != 0)
     || (file_size < STL_MIN_FILE_SIZE))
    {
      cerr << "The file has the wrong size " << endl;
      throw(-1);
    }

  num_facets = (file_size - HEADER_SIZE) / SIZEOF_STL_FACET;

  /* Read the header */
  fread(this_header, LABEL_SIZE, 1, fp);
  this_header[80] = '\0';

  /* Read the int following the header.  This should contain # of facets */
  //header_num_facets = stl_get_little_int(fp);
  //if(num_facets != header_num_facets)
  //  {
  //    cout << "Warning: File size doesn't match number of facets in the header " << endl;
  //  }
  fclose(fp);

  return num_facets;
}

void Stl::initFromFileBinary(const char * filename) {

  int i,iv;
  //char extra[2];
  FILE *fp;


  nf = stl_get_number_of_facet(filename);
  nv = 3*nf;
  // allocate...
  normal = new float[nf][3];
  v_of_f = new int[nf][3];
  xv     = new float[nv][3];

  ns = 1; //binary files cannot list multiple solids
  sof    = new int[nf];
  solidName.resize(ns);

  string t_fname = filename;
  string delimeter = "/";
  size_t pos = 0;
  while ((pos = t_fname.find(delimeter)) != std::string::npos) {
    string token = t_fname.substr(0,pos);
    t_fname.erase(0, pos + delimeter.length());
  }

  solidName[0] = t_fname.c_str();

  /* Open the file */
  fp = fopen(filename, "rb");
  if(fp == NULL) {
    cerr << "Couldn't open " << filename <<  "  for reading" << endl;
    throw(20);
  }
  //fseek(fp, HEADER_SIZE, SEEK_SET);
  //fread(this_header, LABEL_SIZE, 1, fp);
  fseek(fp, HEADER_SIZE, SEEK_SET);

  for(i = 0; i < nf; i++)
    {
      normal[i][0] =  stl_get_little_float(fp);
      normal[i][1] =  stl_get_little_float(fp);
      normal[i][2] =  stl_get_little_float(fp);
      iv = 0;
      xv[3*(i)+iv][0] = stl_get_little_float(fp);
      xv[3*(i)+iv][1] = stl_get_little_float(fp);
      xv[3*(i)+iv][2] = stl_get_little_float(fp);
      v_of_f[i][iv] = 3*(i) + iv;
      iv++;
      xv[3*(i)+iv][0] = stl_get_little_float(fp);
      xv[3*(i)+iv][1] = stl_get_little_float(fp);
      xv[3*(i)+iv][2] = stl_get_little_float(fp);
      v_of_f[i][iv] = 3*(i) + iv;
      iv++;
      xv[3*(i)+iv][0] = stl_get_little_float(fp);
      xv[3*(i)+iv][1] = stl_get_little_float(fp);
      xv[3*(i)+iv][2] = stl_get_little_float(fp);
      v_of_f[i][iv] = 3*(i) + iv;
      //extra[0] = fgetc(fp);
      //extra[1] = fgetc(fp);
      fgetc(fp);
      fgetc(fp);
      //cout << " normal " << normal[i][0] << " " << normal[i][1] << " " << normal[i][2] << endl;
      //cout << "extra " << extra[0] << extra[1] << endl;

      sof[i] = 0; //only one solid for binary files
    }

  fclose(fp);

}

/// stuff for ASCII STLs

void Stl::initFromFile(const char * filename) {
  cout << " > reading ascii formatted file " << filename << endl;
  FILE * fp;
  if ( (fp = fopen(filename,"r")) == NULL ) {
    cerr << " > Error: stl file " << filename << " cannot be opened." << endl;
    throw(20);
  }

  asciiCountObjects(fp);

  // allocate...
  normal = new float[nf][3];
  v_of_f = new int[nf][3];
  xv     = new float[nv][3];
  sof    = new int[nf];
  solidName.resize(ns);

  asciiReadObjects(fp);

  fclose(fp);
}

void Stl::repair() {
  cout << " > Stl::repair()" << endl;

  const bool verbose = true;
  int * fa_flag = new int[nf];

  // now walk the faces in order to flip everyone into a specified orientation...
  for (int ifa = 0; ifa < nf ; ++ifa)
    fa_flag[ifa] = 0;

  std::vector<int> flipFaGroups;
  // read which face groups require normal flipping
  if (Param * param = Params::getParam("STL_FLIP_GROUPS")) {
    int iarg = 0;
    while (iarg < param->size()) {
      flipFaGroups.push_back(param->getInt(iarg++));
    }
  }


  int group = 0;
  while ( true ) {

    ++group;


    stack< pair<int,int> > faces;
    int ifa_seed = -1;
    for (int ifa = 0 ; ifa < nf ; ++ifa) {
      if ( fa_flag[ifa] ==0 ) {
        ifa_seed = ifa;
        fa_flag[ifa] = group;

        // flip seed normal for user defined face groups
        if (std::find(flipFaGroups.begin(), flipFaGroups.end(), group) != flipFaGroups.end()) {
                cout << "    > flipping normals for face group " << group << endl;
                swap(v_of_f[ifa][0], v_of_f[ifa][1]);
                swap(nb_of_f[ifa][1], nb_of_f[ifa][2]);
        }

        break;
      }
    }//ifa

    if ( ifa_seed == -1) {
      --group;
      break; // done
    }

    if ( verbose) {
      cout << "    > Starting with ifa seed : " << ifa_seed << " in zone "
           << solidName[sof[ifa_seed]] << endl;
    }

    // we need to push the nbrs of the seed face onto the stack.
    for (int i =0 ; i <3 ; ++i) {
      const int ifa_nbr = nb_of_f[ifa_seed][i];
      if ( ifa_nbr >= 0 ) {
        faces.push(pair<int,int>(ifa_nbr,ifa_seed));
      }
    }//i

    while ( faces.size() > 0 ) {

      pair<int,int> front = faces.top(); faces.pop();
      const int ifa       = front.first;
      const int ifa_nbr   = front.second;
      assert( ifa != -1);

      if ( fa_flag[ifa] != 0 ) {
        assert(fa_flag[ifa] == group );
        continue;
      }

      // need to decide if we need to flip the orientation
      // of this face based on the edge alignment with ifa_nbr
      int ii, jj=0;
      for (ii=0; ii < 3 ; ++ii) {
        int ino0 = v_of_f[ifa][ii];
        int ino1 = v_of_f[ifa][(ii+1)%3];

        const int this_nbr = nb_of_f[ifa][ii];
        if ( ifa_nbr == this_nbr ) {

          for (jj =0; jj < 3 ; ++jj) {
            if ( (v_of_f[ifa_nbr][jj] == ino0) && (v_of_f[ifa_nbr][(jj+1)%3] ==ino1) ) {
              // the edge should appear in the opposite orientation...
              jj = -jj-1;
              break;
            } else if ( (v_of_f[ifa_nbr][jj] == ino1) && (v_of_f[ifa_nbr][(jj+1)%3] == ino0) ) {
              break;
            }
          }
          assert(jj != 3);
          break;
        }
      }//ii


      if ( jj < 0) {
        // the face orientation needs to be flipped.
        jj = ii;
        assert( (jj >= 0) && (jj < 3)) ;

        //cout << " flip ifa : " << ifa << endl;

        swap(v_of_f[ifa][jj], v_of_f[ifa][(jj+1)%3]);
        swap(nb_of_f[ifa][(jj+1)%3], nb_of_f[ifa][(jj+2)%3]);

      }

      // now add the faces..
      for (int i=0; i < 3 ; ++i) {
        const int this_nbr = nb_of_f[ifa][i];
        if ( this_nbr >= 0 ) {
          if ( fa_flag[this_nbr] ==0 ) {
            faces.push(pair<int,int>(this_nbr,ifa));
          } else {
            assert( fa_flag[this_nbr] == group );
          }
        }
      }

      // and complete this face.
      fa_flag[ifa] = group;
    }//stack.size() >0

  }//while(true)

  if ( verbose )
    cout << "    > Number of face groups : " << group << endl;

  for (int igr =1; igr <= group ; ++igr) {

    // compute the volume for each group..
    double vol = 0.0;
    double gcl[3] = {0.0,0.0,0.0};
    int nfa_gr = 0;
    for (int ifa = 0; ifa < nf ; ++ifa) {
      if ( fa_flag[ifa] == igr ) {
        ++nfa_gr;
        double xv0[3], xv1[3], xv2[3];
        for (int i =0; i < 3 ; ++i) {
          xv0[i] = double(xv[v_of_f[ifa][0]][i]);
          xv1[i] = double(xv[v_of_f[ifa][1]][i]);
          xv2[i] = double(xv[v_of_f[ifa][2]][i]);
        }
        const double normal[3] = TRI_NORMAL_2(xv0, xv1, xv2);
        const double dx[3] = DIFF(xv[v_of_f[ifa][0]],xv[0]); // just use a common point for ALL tris

        for (int i =0; i < 3 ; ++i)
          gcl[i] += 0.5*normal[i];

        vol += DOT_PRODUCT(normal,dx)/6.0;
        fa_flag[ifa] = -1;
      }
    }//ifa

    cout << "      igroup[" << igr-1 << "] nfa, vol, gcl: "
      << nfa_gr << "    " << vol << "   " << COUT_VEC(gcl) << endl;
  }//igr

  // check everyone was found.
  for (int ifa =0; ifa < nf ; ++ifa)
    assert( fa_flag[ifa] == -1);

  delete[] fa_flag;

}//repair

void Stl::compress(float tol) {

  cout << " > Stl::compress(), tol=" << tol << endl;
  // the nb_of_f contains the facet nb for each edge, where the
  // edges are v0-v1, v1-v2, v2-v0. At the end of processing this
  // will either contain a facet number or -1, so it can be used
  // to determine where the stl's holes are, if any

  // push the stl facets into a float bounding box ADT...
  float (*bbmin)[3] = new float[nf][3];
  float (*bbmax)[3] = new float[nf][3];

  int iv;
  for (int ifacet = 0; ifacet < nf; ifacet++) {
    iv = v_of_f[ifacet][0];

    for (int i = 0; i < 3; i++) {
      bbmin[ifacet][i] = xv[iv][i];
      bbmax[ifacet][i] = xv[iv][i];
    }

    for (int j = 1; j < 3; j++) {
      int iv = v_of_f[ifacet][j];
      for (int i = 0; i < 3; i++) {
        bbmin[ifacet][i] = min(xv[iv][i],bbmin[ifacet][i]);
        bbmax[ifacet][i] = max(xv[iv][i],bbmax[ifacet][i]);
      }
    }
  }

  // push the adt into a bounding box...
  fadt = new Adt<float>(nf,bbmin,bbmax);

  // now stored in fadt...
  delete[] bbmin;
  delete[] bbmax;

  int nbb,bbList[ADT_LIST_MAX];

  nb_of_f = new int[nf][3];

  // use nb_of_f to store facet pairings,
  // where -1 indicates no pairing found...
  for (int ifacet = 0; ifacet < nf; ifacet++) {
    for (int e = 0; e < 3; e++) {
      nb_of_f[ifacet][e] = -1;
    }
  }

  // initially v_flag[iv] contains iv (all vertices assumed unique)...
  int * v_flag = new int[nv];
  for (int iv = 0; iv < nv; iv++)
    v_flag[iv] = iv;

  // now cycle through the edges and connect...
  float max_mergedTol = 0.0;    // max tolerance for inexact merges
  int mergedWithTolCount = 0;   // count inexact merges
  int mergeFlipped = 0;         // edges merged in flipped orientation
  int mergeExpected = 0;         // edges merged in expected orientation
  float min_unMergedTol = 1e20; // min tolerance tried for unmerged edges
  int unmatched_edge_count = 0; // count unmerged edges

  int nedge = 0;
  for (int ifacet = 0; ifacet < nf; ifacet++) {
    // unpaired edges only? - may have to rethink for certain tol > 0 cases...

    // define facet tolerance based on 5% of max edge length
    float maxLength = 0.0;
    float eLength;
    for (int iedge = 0; iedge < 3; iedge++) {
      // get this edge's vertices...
      const int iv0 = v_of_f[ifacet][iedge];
      const int iv1 = v_of_f[ifacet][(iedge+1)%3];
      float dx[3];
      for (int i = 0; i < 3; i++) {
        dx[i] = xv[iv1][i] - xv[iv0][i];
      }
      eLength = sqrt(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]);
      maxLength = max(maxLength,eLength);
    }

    // not a good idea -- changes tol for good...
    /*
    float tol_e = 0.05 * maxLength;
    // switch to tol facet tolerance if user did not provide any tolerance
    if (tol < 1.0e-15)
      tol = tol_e;
    */

    float tol2 = tol*tol;

    // loop through edges, use iv0 to build vertex bounding box based on tolerance
    for (int iedge = 0; iedge < 3; iedge++) if (nb_of_f[ifacet][iedge] == -1) {
      nedge++;
      // get this edge's vertices...
      int iv0 = v_of_f[ifacet][iedge];
      int iv1 = v_of_f[ifacet][(iedge+1)%3];

      // use iv0's bounding box to build the list of candidates...
      float this_bbmin[3],this_bbmax[3];

      for (int i = 0; i < 3; i++) {
        this_bbmin[i] = xv[iv0][i] - tol;
        this_bbmax[i] = xv[iv0][i] + tol;
      }
      // use the adt to get the candidate facets with overlap...
      fadt->buildListForBBox(nbb,bbList,this_bbmin,this_bbmax);

      // cycle through the list and find the matching edge (if any)...
      int found = 0; // check to make sure we find ourself
      int ifacet2_match = -1;
      int iedge2_match = -1,iv0_match = -1,iv1_match = -1; 
      bool edgeFlip = false;
      float min_d2 = 0.0;
      for (int ibb = 0; ibb < nbb; ibb++) {
        if (bbList[ibb] == ifacet) {
          found = 1;
        }
        else {
          // cycle through the edges of this facet and look for matches...
          int ifacet2 = bbList[ibb];
          for (int iedge2 = 0; iedge2 < 3; iedge2++) if (nb_of_f[ifacet2][iedge2] == -1) {

            // get iedge2's vertices...
            const int iv0_2 = v_of_f[ifacet2][iedge2];
            const int iv1_2 = v_of_f[ifacet2][(iedge2+1)%3];

            // ************************************************
            // look for a match in the expected orientation...
            // ************************************************

            // check that iedge2's second vertex within tolerance
            float this_d2_0 = 0.0;
            for (int i = 0; i < 3; i++) {
              float dx = xv[iv0][i] - xv[iv1_2][i];
              this_d2_0 += dx*dx;
            }

            if (this_d2_0 <= tol2) {
              // check that iedge2's first vertex is within tolerance
              float this_d2_1 = 0.0;
              for (int i = 0; i < 3; i++) {
                float dx = xv[iv1][i] - xv[iv0_2][i];
                this_d2_1 += dx*dx;
              }
              if (this_d2_1 <= tol2) {
                float this_d2 = max(this_d2_0,this_d2_1);

                if ((ifacet2_match == -1)||(this_d2 < min_d2)) {
                  ifacet2_match = ifacet2;
                  iedge2_match = iedge2;
                  min_d2 = this_d2;
                  iv0_match = iv1_2;
                  iv1_match = iv0_2;
                  edgeFlip = false;
                }
              }
            }

            // ************************************************
            // look for a match in the flipped orientation...
            // ************************************************

            // check that iedge2's first vertex within tolerance
            this_d2_0 = 0.0;
            for (int i = 0; i < 3; i++) {
              float dx = xv[iv0][i] - xv[iv0_2][i];
              this_d2_0 += dx*dx;
            }

            if (this_d2_0 <= tol2) {
              // check that iedge2's second vertex is within tolerance
              float this_d2_1 = 0.0;
              for (int i = 0; i < 3; i++) {
                float dx = xv[iv1][i] - xv[iv1_2][i];
                this_d2_1 += dx*dx;
              }
              if (this_d2_1 <= tol2) {
                float this_d2 = max(this_d2_0,this_d2_1);

                if ((ifacet2_match == -1)||(this_d2 < min_d2)) {
                  ifacet2_match = ifacet2;
                  iedge2_match = iedge2;
                  min_d2 = this_d2;
                  iv0_match = iv0_2;
                  iv1_match = iv1_2;
                  edgeFlip = true;
                }
              }
            }

          }
        }
      }
      max_mergedTol = max(max_mergedTol,min_d2);

      // check...
      if (found == 0) {
        cerr << "Error: did not find ourself in fadt." << endl;
        throw(-1);
      }

      if (ifacet2_match >= 0) {

        (edgeFlip) ? mergeFlipped++ : mergeExpected++;
        if (min_d2 > 0.0) mergedWithTolCount++;

        // set the connection in both ifacet and ifacet2_match...
        nb_of_f[ifacet][iedge] = ifacet2_match;
        nb_of_f[ifacet2_match][iedge2_match] = ifacet;

        // iv0 and iv0_match are the same vertex...
        // push down in both lists until we get to a unique vertex...
        while (v_flag[iv0] != iv0)
          iv0 = v_flag[iv0];
        while (v_flag[iv0_match] != iv0_match)
          iv0_match = v_flag[iv0_match];
        if (iv0 < iv0_match) {
          v_flag[iv0_match] = iv0;
        }
        else if (iv0 > iv0_match) {
          v_flag[iv0] = iv0_match;
        }

        // iv1 and iv1_match are the same vertex...
        // push down in both lists until we get to a unique vertex...
        while (v_flag[iv1] != iv1)
          iv1 = v_flag[iv1];
        while (v_flag[iv1_match] != iv1_match)
          iv1_match = v_flag[iv1_match];
        if (iv1 < iv1_match) {
          v_flag[iv1_match] = iv1;
        }
        else if (iv1 > iv1_match) {
          v_flag[iv1] = iv1_match;
        }

      }
      else {
        unmatched_edge_count++;
        min_unMergedTol = min(min_unMergedTol,tol);
      }

    }

  }


  // and compress...
  for (int iv = 0; iv < nv; iv++) {
    if (v_flag[iv] != iv) {
      // push this v_flag down to a unique index...
      int iv2 = v_flag[iv];
      while (v_flag[iv2] != iv2)
        iv2 = v_flag[iv2];
      v_flag[iv] = iv2;
    }
  }

  int nv_new = 0;
  for (int iv = 0; iv < nv; iv++) {
    if (v_flag[iv] == iv) {
      if (nv_new != iv) {
        // copy down the unique vertex location...
        xv[nv_new][0] = xv[iv][0];
        xv[nv_new][1] = xv[iv][1];
        xv[nv_new][2] = xv[iv][2];
      }
      v_flag[iv] = nv_new++;
    }
    else {
      v_flag[iv] = v_flag[v_flag[iv]]; // get the index...
    }
  }

  cout << "   vertices:" << endl;
  cout << "    > compressed from " << nv << " to " << nv_new << endl;
  cout << "   edges:" << endl;
  cout << "    > tolerance-based merges: " << mergedWithTolCount << " (max tolerance: " << sqrt(max_mergedTol) << ")" << endl;
  cout << "    > standard orientation merges: " << mergeExpected << endl;
  cout << "    > flipped orientation merges: " << mergeFlipped << endl;
  if (unmatched_edge_count > 0) {
    cout << "    > unmatched_edges: " << unmatched_edge_count << " (minimum tolerance: " << min_unMergedTol << ")" << endl;
  }
#ifdef DEBUG
  cout << "Euler: nv_new - nedge + nf - 2 (expect zero)  = " << nv_new - nedge + nf - 2 << endl;
#endif

  // compress the v_of_f...
  for (int ifacet = 0; ifacet < nf; ifacet++) {
    for (int i = 0; i < 3; i++) {
      v_of_f[ifacet][i] = v_flag[v_of_f[ifacet][i]];
    }
  }

  delete[] v_flag;

  // clear xv memory...
  float (*xv_tmp)[3] = new float[nv_new][3];
  for (int iv = 0; iv < nv_new; iv++) {
    for (int i = 0; i < 3; i++) {
      xv_tmp[iv][i] = xv[iv][i];
    }
  }

  delete[] xv;
  xv = xv_tmp;

  // and take the new nv...
  nv = nv_new;
}

void Stl::writeTecplot(const char * filename) {

#ifdef DEBUG
  cout << "Stl::writeTecplot: " << filename << "..." << endl;
#endif

  FILE * fp;
  if ( (fp=fopen(filename,"w"))==NULL ) {
    cerr << "Error: cannot open file " << filename << endl;
    throw(-1);
  }

  fprintf(fp,"TITLE = \"writeTecplot\"\n");
  fprintf(fp,"VARIABLES = \"X\"\n");
  fprintf(fp,"\"Y\"\n");
  fprintf(fp,"\"Z\"\n");

  char dummy[128];
  sprintf(dummy,"faces");
  writeTecplotZone(fp,dummy);

  fclose(fp);

}

void Stl::writeTecplotZone(FILE * fp,char * zonename) {

  fprintf(fp,"ZONE T=\"%s\"\n",zonename);
  fprintf(fp,"N=%d, E=%d, F=FEPOINT, ET=TRIANGLE\n",getNv(),getNf());

  int iv;
  for (iv = 0; iv < getNv(); iv++) {
    fprintf(fp,"%lf %lf %lf\n",xv[iv][0],xv[iv][1],xv[iv][2]);
  }

  // connectivity...
  int ifa;
  for (ifa = 0; ifa < getNf(); ifa++) {
    fprintf(fp,"%d %d %d \n",v_of_f[ifa][0]+1,v_of_f[ifa][1]+1,v_of_f[ifa][2]+1);
  }

}

void Stl::writeFlaggedFacetsTecplot(const char * filename,int flag) {

#ifdef DEBUG
  std::cout << "Stl::writeFlaggedFacetsTecplot: " << filename << " flag = " << flag << "..." << std::endl;
#endif

  FILE * fp;
  if ( (fp=fopen(filename,"w"))==NULL ) {
    cerr << "Error: cannot open file " << filename << endl;
    throw(-1);
  }

  fprintf(fp,"TITLE = \"writeFlaggedFacetsTecplot\"\n");
  fprintf(fp,"VARIABLES = \"X\"\n");
  fprintf(fp,"\"Y\"\n");
  fprintf(fp,"\"Z\"\n");

  char dummy[128];
  sprintf(dummy,"flagged_faces");
  writeFlaggedFacetsTecplotZone(fp,dummy,flag);

  fclose(fp);

}

void Stl::writeFlaggedFacetsTecplotZone(FILE * fp,char * zonename,int flag) {

  assert(int(facetFlag.size()) == nf);

  fprintf(fp,"ZONE T=\"%s\"\n",zonename);

  // count faces and vertices...
  int fcount = 0;
  int ifa;
  for (ifa = 0; ifa < getNf(); ifa++) {
    if (facetFlag[ifa] == flag)
      fcount++;
  }

  fprintf(fp,"N=%d, E=%d, F=FEPOINT, ET=TRIANGLE\n",getNv(),fcount);

  int iv;
  for (iv = 0; iv < getNv(); iv++) {
    fprintf(fp,"%lf %lf %lf\n",xv[iv][0],xv[iv][1],xv[iv][2]);
  }

  // connectivity...
  for (ifa = 0; ifa < getNf(); ifa++) {
    if (facetFlag[ifa] == flag)
      fprintf(fp,"%d %d %d \n",v_of_f[ifa][0]+1,v_of_f[ifa][1]+1,v_of_f[ifa][2]+1);
  }

}
