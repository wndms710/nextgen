#include "SimpleSurface.hpp"

void SimpleSurface::writeTecplot(const string& filename) const {

  COUT1("Surface::writeTecplot(" << filename << ")");

  IntFlag st_flag(nst);
  st_flag.setAll(1);  // currently write entire surface
  writeSelectedFacesByZoneToTecplot(filename,st_flag);

}

void SimpleSurface::writeFlaggedZonesToTecplot(const string& filename) {
  COUT2("SimpleSurface::writeFlaggedZonesToTecplot("<<filename<<")");

  if (zone_flag.count() == 0) {
    CWARN("no zones were flagged for export; skipping");
    return; // no tris were flagged
  }

  st_flag.setLength(nst);
  st_flag.setAll(0);

  FOR_IST {
    if (zone_flag[znost[ist]]) st_flag[ist] = 1;
  }
  writeSelectedFacesByZoneToTecplot(filename,st_flag);
}

void SimpleSurface::writeSelectedFacesByZoneToTecplot(const string& filename, const IntFlag& st_flag,const bool b_flag) const {

  COUT2("SimpleSurface::writeSelectedFacesByZoneToTecplot("<<filename<<")");

  FILE * fp = fopen(filename.c_str(),"w");
  assert(fp != NULL);

  fprintf(fp,"TITLE = \"%s\"\n",filename.c_str());
  fprintf(fp,"VARIABLES = \"X\"\n");
  fprintf(fp,"\"Y\"\n");
  fprintf(fp,"\"Z\"\n");
  if (b_flag) fprintf(fp,"\"st_flag_sum\"\n");

  cout << "nst: " << nst << " nsp: " << nsp << " zoneVec.size(): " << zoneVec.size() << endl;

  IntFlag sp_flag(nsp);
  sp_flag.setAll(-1);
  IntFlag _zone_flag( zoneVec.size() );
  _zone_flag.setAll(0);
  // flag nodes and zones being dumped

  double * sp_flag_value = NULL;
  if (b_flag) {
    sp_flag_value = new double[nsp];
    for (int isp=0; isp<nsp; ++isp) sp_flag_value[isp] = 0.0;
  }

  int nst_selected = 0;
  for (int ist = 0; ist < nst; ++ist) {
    if (st_flag[ist] != 0) {
      ++nst_selected;
      if (_zone_flag[znost[ist]] == 0) _zone_flag[znost[ist]] = 1;

      if (b_flag) {
        for (int i=0; i<3; ++i) {
          const int isp = spost[ist][i];
          sp_flag_value[isp] += double(st_flag[ist]);
        }
      }
    }
  }

  // loop zones and write points and tris
  for (int izone = 0; izone < int(zoneVec.size()); ++izone) {
    if (_zone_flag[izone] > 0) {
      for (int isp = 0; isp < nsp; ++isp) sp_flag[isp] = -1;

      int nst_selected = 0;
      for (int ist = 0; ist < nst; ++ist) {
        if (st_flag[ist] != 0 && znost[ist] == izone) {
          ++nst_selected;
          FOR_I3 sp_flag[spost[ist][i]] = 0;
        }
      }

      int nsp_selected = 0;
      for (int isp = 0; isp < nsp; ++isp) {
        if (sp_flag[isp] == 0) {
          sp_flag[isp] = nsp_selected++;
        }
      }

      // zone header
      fprintf(fp,"ZONE T=\"%s\"\n",zoneVec[izone].getName().c_str());
      fprintf(fp,"N=%d, E=%d, F=FEPOINT, ET=TRIANGLE\n",nsp_selected,nst_selected);

      if (b_flag) {
        for (int isp = 0; isp < nsp; ++isp) {
          if (sp_flag[isp] >= 0) {
            fprintf(fp,"%18.15le %18.15le %18.15le %18.15le\n",xsp[isp][0],xsp[isp][1],xsp[isp][2],sp_flag_value[isp]);
          }
        }
      }
      else {
        for (int isp = 0; isp < nsp; ++isp) {
          if (sp_flag[isp] >= 0) {
            fprintf(fp,"%18.15le %18.15le %18.15le\n",xsp[isp][0],xsp[isp][1],xsp[isp][2]);
          }
        }
      }


      for (int ist = 0; ist < nst; ++ist) {
        if (st_flag[ist] != 0 && znost[ist] == izone) {
          fprintf(fp,"%d %d %d\n",
          sp_flag[spost[ist][0]]+1,
          sp_flag[spost[ist][1]]+1,
          sp_flag[spost[ist][2]]+1);  //tecplot file is 1-indexed
        }
      }
    }
  }
  DELETE(sp_flag_value);
  fclose(fp);
}

void SimpleSurface::writeSelectedFacesByZoneToTecplot(const string& filename, const IntFlag& st_flag,const double x0[3]) const {

  // this one centers the output about the passed x0...

  cout << "SimpleSurface::writeSelectedFacesByZoneToTecplot(" << filename << ") x0: " << COUT_VEC(x0) << endl;

  FILE * fp = fopen(filename.c_str(),"w");
  assert(fp != NULL);

  fprintf(fp,"TITLE = \"%s\"\n",filename.c_str());
  fprintf(fp,"VARIABLES = \"X\"\n");
  fprintf(fp,"\"Y\"\n");
  fprintf(fp,"\"Z\"\n");

  cout << "nst: " << nst << " nsp: " << nsp << " zoneVec.size(): " << zoneVec.size() << endl;

  IntFlag sp_flag(nsp);
  sp_flag.setAll(-1);
  IntFlag _zone_flag( zoneVec.size() );
  _zone_flag.setAll(0);
  // flag nodes and zones being dumped

  int nst_selected = 0;
  for (int ist = 0; ist < nst; ++ist) {
    if (st_flag[ist] != 0) {
      ++nst_selected;
      if (_zone_flag[znost[ist]] == 0) _zone_flag[znost[ist]] = 1;
      FOR_I3 {
        const int isp = spost[ist][i];
        sp_flag[isp] = 0;
      }
    }
  }

  // loop zones and write points and tris
  for (int izone = 0; izone < int(zoneVec.size()); ++izone) {
    if (_zone_flag[izone] > 0) {
      for (int isp = 0; isp < nsp; ++isp) sp_flag[isp] = -1;

      int nst_selected = 0;
      for (int ist = 0; ist < nst; ++ist) {
        if (st_flag[ist] != 0 && znost[ist] == izone) {
          ++nst_selected;
          FOR_I3 sp_flag[spost[ist][i]] = 0;
        }
      }

      int nsp_selected = 0;
      for (int isp = 0; isp < nsp; ++isp) {
        if (sp_flag[isp] == 0) {
          sp_flag[isp] = nsp_selected++;
        }
      }

      // zone header
      fprintf(fp,"ZONE T=\"%s\"\n",zoneVec[izone].getName().c_str());
      fprintf(fp,"N=%d, E=%d, F=FEPOINT, ET=TRIANGLE\n",nsp_selected,nst_selected);

      for (int isp = 0; isp < nsp; ++isp) {
        if (sp_flag[isp] >= 0) {
          fprintf(fp,"%18.15le %18.15le %18.15le\n",xsp[isp][0]-x0[0],xsp[isp][1]-x0[1],xsp[isp][2]-x0[2]);
        }
      }

      for (int ist = 0; ist < nst; ++ist) {
        if (st_flag[ist] != 0 && znost[ist] == izone) {
          fprintf(fp,"%d %d %d\n",
          sp_flag[spost[ist][0]]+1,
          sp_flag[spost[ist][1]]+1,
          sp_flag[spost[ist][2]]+1);  //tecplot file is 1-indexed
        }
      }
    }
  }

  fclose(fp);

}

void SimpleSurface::writeSpDataToTecplot(const string& filename,const double *sp_buf) {

  FILE * fp;
  fp = fopen(filename.c_str(),"w");
  assert(fp != NULL);
  fprintf(fp,"VARIABLES = \"X\"\n");
  fprintf(fp,"\"Y\"\n");
  fprintf(fp,"\"Z\"\n");
  fprintf(fp,"\"VAR\"\n");
  fprintf(fp,"ZONE T=\"SURFACE\"\n");
  fprintf(fp,"N=%d, E=%d, F=FEPOINT, ET=TRIANGLE\n",nst*3,nst);
  for (int ist = 0; ist < nst; ++ist) {
    const double * const x0 = xsp[spost[ist][0]];
    const double * const x1 = xsp[spost[ist][1]];
    const double * const x2 = xsp[spost[ist][2]];
    const double d0 = sp_buf[spost[ist][0]];
    const double d1 = sp_buf[spost[ist][1]];
    const double d2 = sp_buf[spost[ist][2]];
    fprintf(fp,"%18.15le %18.15le %18.15le %18.15le\n",x0[0],x0[1],x0[2],d0);
    fprintf(fp,"%18.15le %18.15le %18.15le %18.15le\n",x1[0],x1[1],x1[2],d1);
    fprintf(fp,"%18.15le %18.15le %18.15le %18.15le\n",x2[0],x2[1],x2[2],d2);
  }
  for (int ist = 0; ist < nst; ++ist) {
    fprintf(fp,"%d %d %d\n",ist*3+1,ist*3+2,ist*3+3);
  }
  fclose(fp);
}

void SimpleSurface::writeSpDataToTecplot(const string& filename,const double (*sp_buf3)[3]) {

  FILE * fp;
  fp = fopen(filename.c_str(),"w");
  assert(fp != NULL);
  fprintf(fp,"VARIABLES = \"X\"\n");
  fprintf(fp,"\"Y\"\n");
  fprintf(fp,"\"Z\"\n");
  fprintf(fp,"\"VAR_X\"\n");
  fprintf(fp,"\"VAR_Y\"\n");
  fprintf(fp,"\"VAR_Z\"\n");
  fprintf(fp,"ZONE T=\"SURFACE\"\n");
  fprintf(fp,"N=%d, E=%d, F=FEPOINT, ET=TRIANGLE\n",nst*3,nst);
  for (int ist = 0; ist < nst; ++ist) {
    const double * const x0 = xsp[spost[ist][0]];
    const double * const x1 = xsp[spost[ist][1]];
    const double * const x2 = xsp[spost[ist][2]];
    double d0[3]; FOR_I3 d0[i] = sp_buf3[spost[ist][0]][i];
    double d1[3]; FOR_I3 d1[i] = sp_buf3[spost[ist][1]][i];
    double d2[3]; FOR_I3 d2[i] = sp_buf3[spost[ist][2]][i];
    fprintf(fp,"%18.15le %18.15le %18.15le %18.15le %18.15le %18.15le\n",x0[0],x0[1],x0[2],d0[0],d0[1],d0[2]);
    fprintf(fp,"%18.15le %18.15le %18.15le %18.15le %18.15le %18.15le\n",x1[0],x1[1],x1[2],d1[0],d1[1],d1[2]);
    fprintf(fp,"%18.15le %18.15le %18.15le %18.15le %18.15le %18.15le\n",x2[0],x2[1],x2[2],d2[0],d2[1],d2[2]);
  }
  for (int ist = 0; ist < nst; ++ist) {
    fprintf(fp,"%d %d %d\n",ist*3+1,ist*3+2,ist*3+3);
  }
  fclose(fp);
}

int SimpleSurface::addTecplot(const string& filename) {

  FILE * fp = NULL;
  const int file_err = MiscUtils::openFile(&fp,filename,"rb");
  if (file_err != 0) return file_err;
  fclose(fp);

  ifstream ifile;
  ifile.open(filename.c_str());

  assert(ifile.is_open());

  const int ss_nsp0 = nsp;
  const int ss_nst0 = nst;

  string line;
  vector<string> tokens;

  vector<NewNode> newNodesVec;
  vector<NewTri> newTrisVec;

  COUT1(" > adding zones:");
  while (!ifile.eof()) {

    getline(ifile,line);
    tokens.clear();
    MiscUtils::tokenizeString(tokens,line," =\t,");

    if ( tokens.empty() || (tokens[0] != "ZONE")) continue;

    assert(tokens[0] == "ZONE");
    assert(tokens[1] == "T");
    const int izone = zoneVec.size();
    const int isz = nsz;
    zoneVec.push_back(SurfaceZone(tokens[2]));
    ++nsz;

    COUT2("    > " << zoneVec.back().getName());

    const int sp_offset = ss_nsp0 + newNodesVec.size();

    // zone metadata on next line
    getline(ifile,line);
    tokens.clear();
    MiscUtils::tokenizeString(tokens,line," =\t,");

    assert(tokens[0] == "N");
    const int zone_nsp = atoi(tokens[1].c_str());

    assert(tokens[2] == "E");
    const int zone_nst = atoi(tokens[3].c_str());

    assert(tokens[4] == "F");
    assert(tokens[5] == "FEPOINT");

    assert(tokens[6] == "ET");
    assert(tokens[7] == "TRIANGLE");

    for (int isp = 0; isp < zone_nsp; ++isp) {
      getline(ifile,line);
      double this_xsp[3];
      sscanf(line.c_str(),"%lf %lf %lf",&this_xsp[0],&this_xsp[1],&this_xsp[2]);
      newNodesVec.push_back(NewNode(this_xsp));
      // cout << "just got line: \"" << line << "\" x: " << COUT_VEC(xsp[isp]) << endl;
    }

    for (int ist = 0; ist < zone_nst; ++ist) {
      getline(ifile,line);
      int this_spost[3];
      sscanf(line.c_str(),"%d %d %d",&this_spost[0],&this_spost[1],&this_spost[2]);
      FOR_I3 this_spost[i] += sp_offset-1; // tecplot is 1-indexed
      newTrisVec.push_back(NewTri(this_spost,izone,isz));
    }

    // reached end of zone, should be eof or next zone....
  }
  ifile.close();

  // allocate and populate new values
  nsp += newNodesVec.size();
  growNspData(nsp,ss_nsp0);

  int count = 0;
  for (int isp = 0, nsp_new = newNodesVec.size(); isp < nsp_new; ++isp) {
    FOR_I3 xsp[isp+ss_nsp0][i] = newNodesVec[isp].xsp[i];
    ++count;
  }
  assert(count+ss_nsp0 == nsp);

  nst += newTrisVec.size();
  growNstData(nst,ss_nst0);

  count = 0;
  for (int ist = 0, nst_new = newTrisVec.size(); ist < nst_new; ++ist) {
    FOR_I3 spost[ist+ss_nst0][i] = newTrisVec[ist].spost[i];
    znost[ist+ss_nst0] = newTrisVec[ist].znost;
    ++count;
  }
  assert(count+ss_nst0 == nst);

  if (ss_nsp0 > 0 && ss_nst0 > 0) {
    COUT1(" > nsp=" << nsp << " (+" << newNodesVec.size() << ")");
    COUT1(" > nst=" << nst << " (+" << newTrisVec.size() << ")");
  }
  else {
    COUT1(" > nsp=" << nsp);
    COUT1(" > nst=" << nst);
  }

  return 0;
}
