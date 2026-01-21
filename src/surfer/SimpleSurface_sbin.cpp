#include "SimpleSurface.hpp"
#include "ByteSwap.hpp"
#include "Utils.hpp"

void SimpleSurface::writeBinary(const string& filename) const {

  COUT1("SimpleSurface::writeBinary("<<filename<<")");

  MiscUtils::mkdir_for_file(filename);
  string tmp_filename = MiscUtils::makeTmpPrefix(filename);
  FILE * fp = fopen(tmp_filename.c_str(),"wb");
  assert(fp != NULL);

  const int version = 1;
  cout << " > sbin version: " << version << endl;
  fwrite(&version,sizeof(int),1,fp);
  const int count = zoneVec.size();
  fwrite(&count,sizeof(int),1,fp);
  cout << " > writing " << count << " zones:" << endl;
  for (int izone = 0, limit = zoneVec.size(); izone < limit; ++izone) {
    const int length = zoneVec[izone].getName().length();
    fwrite(&length,sizeof(int),1,fp);
    fwrite(zoneVec[izone].getName().c_str(),sizeof(char),length,fp);
    cout << "    > \"" << zoneVec[izone].getName() << "\"" << endl;
  }
  cout << " > nsp: " << nsp << endl;
  fwrite(&nsp,sizeof(int),1,fp);
  fwrite(xsp,sizeof(double),nsp*3,fp);

  cout << " > nst: " << nst << endl;
  fwrite(&nst,sizeof(int),1,fp);
  fwrite(spost,sizeof(int),nst*3,fp);
  fwrite(znost,sizeof(int),nst,fp);

  // include periodic data if present...
  if (pbi) {
    const int npt = PeriodicData::periodicTransformVec.size();
    cout << " > periodic transforms: " << npt << endl;
    fwrite(&npt,sizeof(int),1,fp);
    for (int ipt = 0; ipt < npt; ++ipt)
      PeriodicData::periodicTransformVec[ipt].writeBinary(fp);
    // write periodic bits for periodic zones only...
    int npz = 0;
    for (int izone = 0, limit = zoneVec.size(); izone < limit; ++izone)
      if (zoneVec[izone].getPeriodicBits())
        ++npz;
    cout << " > npz: " << npz << endl;
    fwrite(&npz,sizeof(int),1,fp);
    int * periodic_zone_array = new int[npz];
    uint2 * periodic_zone_bits = new uint2[npz];
    int ipz = 0;
    for (int izone = 0, limit = zoneVec.size(); izone < limit; ++izone) {
      if (zoneVec[izone].getPeriodicBits()) {
        periodic_zone_array[ipz] = izone;
        periodic_zone_bits[ipz]  = zoneVec[izone].getPeriodicBits();
        ++ipz;
      }
    }
    assert(ipz == npz);
    fwrite(periodic_zone_array,sizeof(int),npz,fp);
    fwrite(periodic_zone_bits,sizeof(uint2),npz,fp);
    delete[] periodic_zone_array;
    delete[] periodic_zone_bits;
    // only write the pbi's that contain transform data...
    int npbi = 0;
    for (int isp = 0; isp < nsp; ++isp)
      if (pbi[isp] != uint8(isp))
        ++npbi;
    cout << " > npbi: " << npbi << endl;
    fwrite(&npbi,sizeof(int),1,fp);
    int * isp_array = new int[npbi];
    int * isp_p_array = new int[npbi];
    uint2 * isp_bits_array = new uint2[npbi];
    int ipbi = 0;
    for (int isp = 0; isp < nsp; ++isp) {
      if (pbi[isp] != uint8(isp)) {
        isp_array[ipbi]      = isp;
        isp_p_array[ipbi]    = int(pbi[isp]&MASK_52BITS);
        isp_bits_array[ipbi] = uint2(pbi[isp]>>52);
        ++ipbi;
      }
    }
    assert(ipbi == npbi);
    fwrite(isp_array,sizeof(int),npbi,fp);
    fwrite(isp_p_array,sizeof(int),npbi,fp);
    fwrite(isp_bits_array,sizeof(uint2),npbi,fp);
    delete[] isp_array;
    delete[] isp_p_array;
    delete[] isp_bits_array;
  }
  fclose(fp);

  remove(filename.c_str());
  rename(tmp_filename.c_str(),filename.c_str());
}

void SimpleSurface::writeSelectedTrisToBinary(const string& filename,const IntFlag& st_flag) const {

  COUT1("Surface::writeSelectedTrisToBinary("<<filename<<")");

  MiscUtils::mkdir_for_file(filename);
  string tmp_filename = MiscUtils::makeTmpPrefix(filename);
  FILE * fp = fopen(tmp_filename.c_str(),"wb");
  assert(fp != NULL);

  const int version = 1;
  fwrite(&version,sizeof(int),1,fp);

  int * sp_flag = new int[nsp];
  for (int isp=0; isp < nsp; ++isp) sp_flag[isp] = -1;

  // flag nodes and zones being dumped
  vector<SurfaceZone> newZoneVec;
  IntFlag zone_flag(zoneVec.size());
  zone_flag.setAll(-1);
  int nst_selected = 0;
  int nz_selected = 0;
  for (int ist=0; ist < nst; ++ist) {
    if (st_flag[ist]) {
      ++nst_selected;
      if (zone_flag[znost[ist]] == -1) {
        zone_flag[znost[ist]] = nz_selected;
        newZoneVec.push_back(SurfaceZone(zoneVec[znost[ist]]));
        ++nz_selected;
      }
      FOR_I3 {
        const int isp = spost[ist][i];
        sp_flag[isp] = 0;
      }
    }
  }

  int zone_count = newZoneVec.size();

  fwrite(&zone_count,sizeof(int),1,fp);
  for (uint izone = 0; izone < newZoneVec.size(); ++izone) {
    const int length = newZoneVec[izone].getName().length();
    COUT2(" > tris from zone: " << newZoneVec[izone].getName());
    fwrite(&length,sizeof(int),1,fp);
    fwrite(newZoneVec[izone].getName().c_str(),sizeof(char),length,fp);
  }

  // count and then write nodes
  int nsp_selected = 0;
  for (int isp=0; isp < nsp; ++isp) {
    if (sp_flag[isp] == 0) {
      sp_flag[isp] = nsp_selected;
      ++nsp_selected;
    }
  }
  fwrite(&nsp_selected,sizeof(int),1,fp);

  for (int isp=0; isp < nsp; ++isp) {
    if (sp_flag[isp] >= 0) {
      FOR_I3 fwrite(&xsp[isp][i],sizeof(double),1,fp);
    }
  }

  fwrite(&nst_selected,sizeof(int),1,fp);
  for (int ist=0; ist < nst; ++ist) {
    if (st_flag[ist]) {
      FOR_I3 fwrite(&sp_flag[spost[ist][i]],sizeof(int),1,fp);
    }
  }
  for (int ist=0; ist < nst; ++ist) {
    if (st_flag[ist]) {
      const int znost_ist = zone_flag[znost[ist]];
      fwrite(&znost_ist,sizeof(int),1,fp);
    }
  }

  fclose(fp);
  DELETE(sp_flag);

  remove(filename.c_str());
  rename(tmp_filename.c_str(),filename.c_str());
}

int SimpleSurface::addBinary(const string& filename) {

  FILE * fp = NULL;
  const int file_err = MiscUtils::openFile(&fp,filename,"rb");
  if (file_err != 0) return file_err;

  bool byte_swap = false;
  int version;
  fread(&version,sizeof(int),1,fp);
  if ((version < 1)||(version > 10)) {
    version = ByteSwap::byteSwap(version);
    if ((version < 1)||(version > 10)) {
      cerr << "Error: file does not start as expected: " << filename << endl;
      return -1;
    }
    else {
      byte_swap = true;
    }
  }
  COUT2(" > sbin version: " << version);
  assert(version == 1);

  // zone reading
  int n_new_zones;
  fread(&n_new_zones,sizeof(int),1,fp);
  if (byte_swap) n_new_zones = ByteSwap::byteSwap(n_new_zones);

  IntFlag new_zone_flag(n_new_zones);

  cout << " > adding " << n_new_zones << " zones:" << endl;
  const int zone_name_len_max = 128;
  char cbuf[zone_name_len_max];
  for (int izone = 0; izone < n_new_zones; ++izone) {
    int length;
    fread(&length,sizeof(int),1,fp);
    if (byte_swap) length = ByteSwap::byteSwap(length);
    assert(length+1 <= zone_name_len_max);
    fread(cbuf,sizeof(char),length,fp); // no byteswap for char
    cbuf[length] = '\0';
    new_zone_flag[izone] = zoneVec.size();
    zoneVec.push_back(SurfaceZone(cbuf));
    cout << "    > \"" << zoneVec[new_zone_flag[izone]].getName() << "\"" << endl;
  }

  // point reading
  const int ss_nsp0 = nsp;
  fread(&nsp,sizeof(int),1,fp);
  if (byte_swap) nsp = ByteSwap::byteSwap(nsp);

  if (ss_nsp0 == 0) cout << " > nsp: " << nsp << endl;  // just the new points
  else {
    nsp += ss_nsp0;
    cout << " > nsp: " << nsp << " (+" << (nsp-ss_nsp0) << ")" << endl;
  }

  growNspData(nsp,ss_nsp0);

  fread(xsp[ss_nsp0],sizeof(double),(nsp-ss_nsp0)*3,fp);
  if (byte_swap) ByteSwap::byteSwap(xsp[ss_nsp0],(nsp-ss_nsp0));

  // tri information
  const int ss_nst0 = nst;

  fread(&nst,sizeof(int),1,fp);
  if (byte_swap) nst = ByteSwap::byteSwap(nst);

  if (ss_nst0 == 0) cout << " > nst: " << nst << endl;
  else {
    nst += ss_nst0;
    cout << " > nst: " << nst << " (+" << (nst-ss_nst0) << ")" << endl;
  }

  growNstData(nst,ss_nst0);

  fread(spost[ss_nst0],sizeof(int),(nst-ss_nst0)*3,fp);
  if (byte_swap) ByteSwap::byteSwap(spost[ss_nst0],(nst-ss_nst0));

  fread(&znost[ss_nst0],sizeof(int),(nst-ss_nst0),fp);
  if (byte_swap) ByteSwap::byteSwap(&znost[ss_nst0],(nst-ss_nst0));

  if (ss_nst0 > 0) {
    // if adding surface, properly offset node & zone indices
    for (int ist = ss_nst0; ist < nst; ++ist) {
      FOR_I3 spost[ist][i] += ss_nsp0;
      znost[ist] = new_zone_flag[znost[ist]];
    }
  }

  // ==============================================================
  // recall that the sbin file can now include periodicity...
  // ==============================================================

  int npt;
  int ierr = fread(&npt,sizeof(int),1,fp);
  if (ierr == 1) {

    // supporting single file read for now
    // I don't really know how to combine periodicTransforms from
    // multiple sbins (we are limited to npt = 3). TODO
    assert(PeriodicData::periodicTransformVec.empty());
    assert(pbi == NULL);

    if (byte_swap) npt = ByteSwap::byteSwap(npt);

    // this sbin has periodic information...
    cout << " > sbin has " << npt << " periodic transforms..." << endl;
    PeriodicData::periodicTransformVec.resize(npt);
    for (int ipt = 0; ipt < npt; ++ipt) {
      PeriodicData::periodicTransformVec[ipt].readBinary(fp,byte_swap);
      cout << "    > " << ipt << ": "; PeriodicData::periodicTransformVec[ipt].dump();
    }

    // read periodic bits for periodic zones only...
    int npz;
    fread(&npz,sizeof(int),1,fp);
    if (byte_swap) npz = ByteSwap::byteSwap(npz);
    //assert(npz >= 2*npt); // surface can hold edge periodicity, so there may be no periodic zones

    int * periodic_zone_array = new int[npz];
    uint2 * periodic_zone_bits = new uint2[npz];

    fread(periodic_zone_array,sizeof(int),npz,fp);
    if (byte_swap) ByteSwap::byteSwap(periodic_zone_array,npz);
    fread(periodic_zone_bits,sizeof(uint2),npz,fp);
    if (byte_swap) ByteSwap::byteSwap(periodic_zone_bits,npz);

    for (int ipz = 0; ipz < npz; ++ipz) {
      const int izone = periodic_zone_array[ipz]; assert((izone >= 0)&&(izone < int(zoneVec.size())));
      zoneVec[izone].setPeriodicBits(periodic_zone_bits[ipz]);

      if (zoneVec[izone].getPeriodicBits() == 1) zoneVec[izone].setBC("PER1_A");
      else if (zoneVec[izone].getPeriodicBits() == 2) zoneVec[izone].setBC("PER1_B");
      else if (zoneVec[izone].getPeriodicBits() == 4) zoneVec[izone].setBC("PER2_A");
      else if (zoneVec[izone].getPeriodicBits() == 8) zoneVec[izone].setBC("PER2_B");
      else if (zoneVec[izone].getPeriodicBits() == 16) zoneVec[izone].setBC("PER3_A");
      else if (zoneVec[izone].getPeriodicBits() == 32) zoneVec[izone].setBC("PER3_B");
    }

    delete[] periodic_zone_array;
    delete[] periodic_zone_bits;

    // only read the pbi's that contain transform data...
    int npbi;
    fread(&npbi,sizeof(int),1,fp); assert(npbi > 0);
    if (byte_swap) npbi = ByteSwap::byteSwap(npbi);

    int * isp_array = new int[npbi];
    int * isp_p_array = new int[npbi];
    uint2 * isp_bits_array = new uint2[npbi];

    fread(isp_array,sizeof(int),npbi,fp);
    if (byte_swap) ByteSwap::byteSwap(isp_array,npbi);
    fread(isp_p_array,sizeof(int),npbi,fp);
    if (byte_swap) ByteSwap::byteSwap(isp_p_array,npbi);
    fread(isp_bits_array,sizeof(uint2),npbi,fp);
    if (byte_swap) ByteSwap::byteSwap(isp_bits_array,npbi);

    pbi = new uint8[nsp];
    for (int isp = 0; isp < nsp; ++isp)
      pbi[isp] = uint8(isp);

    for (int ipbi = 0; ipbi < npbi; ++ipbi) {
      const int isp = isp_array[ipbi];
      const int isp_parent = isp_p_array[ipbi]; assert(isp_parent != isp);
      const uint2 bits = isp_bits_array[ipbi]; assert(bits != 0);
      assert(pbi[isp] == uint8(isp));
      pbi[isp] = BitUtils::packPbiHash(bits,isp_parent);
    }

    delete[] isp_array;
    delete[] isp_p_array;
    delete[] isp_bits_array;

  }

  fclose(fp);

  cout << " > done read" << endl;

  return 0;
}

void SimpleSurface::writeSubzonesToBinary(const string& filename,const vector<int>& subzonesVec) {

  flagTrisFromSubzoneVec(subzonesVec);  // flagged tris in st_flag
  writeSelectedTrisToBinary(filename,st_flag);

}

int SimpleSurface::addPartSurface(const string& filename) {

  // this routine strips the surface from a part file...
  // this should follow the part binary format decribed in code in PartData.cpp read/writeBinary...

  // TODO: byte_swap not implemented yet. Also, no periodicity

  FILE * fp = fopen(filename.c_str(),"rb");
  if (fp == NULL) return -1;
  
  int ibuf[5];
  fread(ibuf,sizeof(int),5,fp);
  cout << " > part file version: " << ibuf[1] << endl;
  if (ibuf[1] != 1) {
    cout << "Error: file does not start as expected \"" << filename << "\"" << endl;
    fclose(fp);
    return -1;
  }
  
  bool b_surf    = (ibuf[2] == 1);
  bool b_ff_surf = (ibuf[3] == 1);
  bool b_pts     = (ibuf[4] == 1);
  cout << " > b_surf: " << b_surf << ", b_ff_surf: " << b_ff_surf << ", b_pts: " << b_pts << endl;
  
  if (!b_surf) {
    cout << "Error: Part file does not contain a surface." << endl;
    fclose(fp);
    return -1;
  }
  
  // read the surface...
  fread(ibuf,sizeof(int),3,fp);
  cout << " > surface nzones: " << ibuf[0] << " nsp: " << ibuf[1] << " nst: " << ibuf[2] << endl;
  
  int n_new_zones = ibuf[0];
  IntFlag new_zone_flag(n_new_zones);

  const int zone_name_len_max = 128;
  char cbuf[zone_name_len_max];
  for (int izone = 0; izone < n_new_zones; ++izone) {
    int length;
    fread(&length,sizeof(int),1,fp);
    //if (byte_swap) length = ByteSwap::byteSwap(length);
    assert(length+1 <= zone_name_len_max);
    fread(cbuf,sizeof(char),length,fp); // no byteswap for char
    cbuf[length] = '\0';
    new_zone_flag[izone] = zoneVec.size();
    zoneVec.push_back(SurfaceZone(cbuf));
    cout << "    > \"" << zoneVec[new_zone_flag[izone]].getName() << "\"" << endl;
  }

  // point reading
  const int ss_nsp0 = nsp;
  nsp = ibuf[1];
    
  if (ss_nsp0 == 0) cout << " > nsp: " << nsp << endl;  // just the new points
  else {
    nsp += ss_nsp0;
    cout << " > nsp: " << nsp << " (+" << (nsp-ss_nsp0) << ")" << endl;
  }

  growNspData(nsp,ss_nsp0);
    
  fread(xsp[ss_nsp0],sizeof(double),(nsp-ss_nsp0)*3,fp);
  //if (byte_swap) ByteSwap::byteSwap(xsp[ss_nsp0],(nsp-ss_nsp0));
    
  // tri information
  const int ss_nst0 = nst;
  nst = ibuf[2];
  //if (byte_swap) nst = ByteSwap::byteSwap(nst);
    
  if (ss_nst0 == 0) cout << " > nst: " << nst << endl;
  else {
    nst += ss_nst0;
    cout << " > nst: " << nst << " (+" << (nst-ss_nst0) << ")" << endl;
  }

  growNstData(nst,ss_nst0);
    
  fread(spost[ss_nst0],sizeof(int),(nst-ss_nst0)*3,fp);
  //if (byte_swap) ByteSwap::byteSwap(spost[ss_nst0],(nst-ss_nst0));
    
  fread(&znost[ss_nst0],sizeof(int),(nst-ss_nst0),fp);
  //if (byte_swap) ByteSwap::byteSwap(&znost[ss_nst0],(nst-ss_nst0));
    
  if (ss_nst0 > 0) {
    // if adding surface, properly offset node & zone indices
    for (int ist = ss_nst0; ist < nst; ++ist) {
      FOR_I3 spost[ist][i] += ss_nsp0;
      znost[ist] = new_zone_flag[znost[ist]];
    }
  }
    
  fclose(fp);

  return 0;
}

    
