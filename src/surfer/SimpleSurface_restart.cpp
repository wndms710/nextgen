#include "SimpleSurface.hpp"
#include "Defs.hpp"
#include "ByteSwap.hpp"
#include "Utils.hpp"

int SimpleSurface::addRestartSurface(const string& filename) {

  FILE * fp = NULL;
  const int file_err = MiscUtils::openFile(&fp,filename,"rb");
  if (file_err != 0) return file_err;

  int byte_swap = 0;
  int itmp[2];
  fread(itmp,sizeof(int),2,fp);
  if (itmp[0] != UGP_IO_MAGIC_NUMBER) {
    ByteSwap::byteSwap(itmp,2);
    if (itmp[0] != UGP_IO_MAGIC_NUMBER) {
      cerr << "Error: file does not start as expected."<< endl;
      return 1; //error
    }
    cout << "File requires byte swapping." << endl;
    byte_swap = 1;
  }

  const int io_version = itmp[1];
  cout << " > io version: " << io_version << endl;
  assert(io_version >= 4);

  // look for sbin in headers, use if present
  // if not present, build surface from restart boundary faces

  //bool b_foundHash = false;
  int8 offset = sizeof(int)*2;
  Header header;
  int8 nno,nbf=-1,nfa,ncv,nno_p,nno_pb=-1;
  int8 nfa_zone[27],noofa_zone[27]; // periodicity
  int8 noobf_s;
  vector<pair<string,int> > zoneNameAndCount;
  int8 zone_count_check = 0;

  int8 xno_offset = -1;
  int8 noobf_offset = -1;

  bool b_foundSurface = false;
  const int ss_nsp0 = nsp;
  const int ss_nst0 = nst;
  int * periodic_zone_array = NULL;
  uint2 * periodic_zone_bits = NULL;
  int npz = 0;

  int done = 0;
  while (done != 1) {

    fseek(fp,offset,SEEK_SET);

    fread(&header,sizeof(Header),1,fp);
    if (byte_swap) ByteSwap::byteSwapHeader(&header,1);

    cout << "header \"" << header.name << "\" id: " << header.id << endl;

    //Since expr parsing added to surfer, no longer use lowercase for all map keys
    string headerNameKey = header.name;
    //transform(headerNameKey.begin(),headerNameKey.end(), headerNameKey.begin(), ::tolower);

    switch (header.id) {
      case UGP_IO_SURFACE:
      {
        cout << "Found surface in restart, nsp=" << header.idata[0] << ", nst=" << header.idata[1] << endl;

        // point reading
        nsp = header.idata[0];
        nsp += ss_nsp0;

        growNspData(nsp,ss_nsp0);

        fread(xsp[ss_nsp0],sizeof(double),(nsp-ss_nsp0)*3,fp);
        if (byte_swap) ByteSwap::byteSwap(xsp[ss_nsp0],(nsp-ss_nsp0));

        // tri information
        nst = header.idata[1];
        nst += ss_nst0;

        growNstData(nst,ss_nst0);

        fread(spost[ss_nst0],sizeof(int),(nst-ss_nst0)*3,fp);
        if (byte_swap) ByteSwap::byteSwap(spost[ss_nst0],(nst-ss_nst0));

        fread(&znost[ss_nst0],sizeof(int),(nst-ss_nst0),fp);
        if (byte_swap) ByteSwap::byteSwap(&znost[ss_nst0],(nst-ss_nst0));

        if (ss_nst0 > 0) {
          // if adding surface, properly offset node indices
          for (int ist = ss_nst0; ist < nst; ++ist) {
            FOR_I3 spost[ist][i] += ss_nsp0;
          }
        }

        b_foundSurface = true;
        //not done...
        //continue to end of the file to get zone info and
        //store R1 and R2 offsets for imaging
        break;
      }
      case UGP_IO_SURFACE_PERIODIC_INFO:
      {
        // only allow one periodic read for now
        assert(pbi == NULL);
        int npt = header.idata[0];
        npz = header.idata[1];
        int npbi = header.idata[2];

        cout << " Periodic surface has " << npt << " periodic transforms, " << npz << " periodic zones and " << npbi << " periodic nodes. " << endl;

        // populate transform...
        assert(PeriodicData::periodicTransformVec.empty()); // if you hit this, your existing surfaces have periodicity, and this must be smarter:
        // i.e. look for common transforms, or clear periodicity if no commonality can be found...
        for (int ipt = 0; ipt < npt; ++ipt) {
          PeriodicData::periodicTransformVec.push_back(PeriodicTransform(header.idata[3+ipt],header.rdata[ipt*3],header.rdata[ipt*3+1],header.rdata[ipt*3+2]));
        }

        // for now just store these (they will be thrown into zoneVec soon)...
        periodic_zone_array = new int[npz];
        periodic_zone_bits = new uint2[npz];
        fread(periodic_zone_array,sizeof(int),npz,fp);
        if (byte_swap) ByteSwap::byteSwap(periodic_zone_array,npz);
        fread(periodic_zone_bits,sizeof(uint2),npz,fp);
        if (byte_swap) ByteSwap::byteSwap(periodic_zone_bits,npz);
        //delete[] periodic_zone_array;
        //delete[] periodic_zone_bits;

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
        for (int isp = ss_nsp0; isp < nsp; ++isp)
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

        break;

      }
      case UGP_IO_NO_FA_BF_CV_COUNTS:
        nno    = ByteSwap::getInt8FromLswMswPair(header.idata+0);
        nbf    = ByteSwap::getInt8FromLswMswPair(header.idata+2);
        nfa    = ByteSwap::getInt8FromLswMswPair(header.idata+4);
        ncv    = ByteSwap::getInt8FromLswMswPair(header.idata+6);
        nno_p  = ByteSwap::getInt8FromLswMswPair(header.idata+8); // first nodes are periodic
        nno_pb = ByteSwap::getInt8FromLswMswPair(header.idata+10); // and boundary

        cout << " > global nno, nbf, nfa, ncv: " << nno << " " << nbf << " "
        << nfa << " " << ncv << " nno_p, nno_pb: " << nno_p << " " << nno_pb << endl;

        // periodicity...
        fread(nfa_zone,sizeof(int8),27,fp);
        fread(noofa_zone,sizeof(int8),27,fp);
        break;
      case UGP_IO_NOOBF_I_AND_V_INT8:
        assert(ByteSwap::getInt8FromLswMswPair(header.idata+0) == nbf);
        noobf_s = ByteSwap::getInt8FromLswMswPair(header.idata+2);

        cout << " > noobf_s: " << noobf_s << endl;
        cout << " > boundary data memory estimate: " <<
        double(sizeof(int)*(nbf+1)+sizeof(int8)*noobf_s+sizeof(double)*3*nno_pb)/1.0E+6 << " [MB]" << endl;

        noobf_offset = offset;

        break;

      case UGP_IO_BF_ZONE_HEADER:
        {
          const int index = header.idata[1];
          const int8 nbf_zone = ByteSwap::getInt8FromLswMswPair(header.idata+2);
          const int8 noobf_zone = ByteSwap::getInt8FromLswMswPair(header.idata+4);
          const int8 ibf_begin = ByteSwap::getInt8FromLswMswPair(header.idata+6);
          const int8 noobf_begin = ByteSwap::getInt8FromLswMswPair(header.idata+8);
          UNUSED(noobf_begin);
          UNUSED(noobf_zone);
          //cout << "index: " << index << " nbf_zone: " << nbf_zone << " noobf_zone: " << noobf_zone << " ibf_begin: " << ibf_begin << " noobf_begin: " << noobf_begin << endl;
          assert(index == int(zoneNameAndCount.size()));
          zoneNameAndCount.push_back(pair<string,int>(header.name,nbf_zone));
          assert(ibf_begin == zone_count_check);
          zone_count_check += nbf_zone;
        }
        break;

      case UGP_IO_X_NO:
        assert(ByteSwap::getInt8FromLswMswPair(header.idata+0) == nno);
        xno_offset = offset;
        break;

      case UGP_IO_EOF:
        done = 1;
    }

    offset += header.skip;

  }

  if (b_foundSurface) {

    // finish zone setup
    if (!pbi) {
      // if the surface isn't periodic just use the bf zone info...

      IntFlag new_zone_flag(zoneNameAndCount.size());

      cout << " > adding zones:" << endl;
      for (int izone = 0, nzn=zoneNameAndCount.size(); izone < nzn; ++izone) {
        // check if this zone already exists...
        new_zone_flag[izone] = zoneVec.size();
        zoneVec.push_back(SurfaceZone(zoneNameAndCount[izone].first));
        cout << "    > \"" << zoneVec[izone].getName() << "\"" << endl;
      }

      if (ss_nst0 > 0) {
        // if adding surface, properly offset zone indices
        for (int ist = ss_nst0; ist < nst; ++ist) {
          znost[ist] = new_zone_flag[znost[ist]];
        }
      }
    }
    else {
      // need to add in periodic zone info to bf zone info to keep SimpleSurface happy...
      // periodic zones are at the end, so this is

      assert(pbi);
      assert(periodic_zone_array);
      assert(periodic_zone_bits);

      IntFlag new_zone_flag(zoneNameAndCount.size()+npz);
      cout << " > adding zones:" << endl;
      for (int izone = 0, nzn=zoneNameAndCount.size()+npz; izone < nzn; ++izone) {
        // check if this zone already exists...
        new_zone_flag[izone] = zoneVec.size();
        if (izone < (nzn-npz)) {
          zoneVec.push_back(SurfaceZone(zoneNameAndCount[izone].first));
          cout << "    > \"" << zoneVec[izone].getName() << "\"" << endl;
        }
        else {
          // add nameless surface zone (named next)...
          zoneVec.push_back(SurfaceZone());
        }
      }

      if (ss_nst0 > 0) {
        // if adding surface, properly offset zone indices
        for (int ist = ss_nst0; ist < nst; ++ist) {
          znost[ist] = new_zone_flag[znost[ist]];
        }
      }

      for (int ipz = 0; ipz < npz; ++ipz) {
        const int izone = periodic_zone_array[ipz]; assert((izone >= int(zoneNameAndCount.size()))&&(izone < int(zoneVec.size())));
        zoneVec[izone].setPeriodicBits(periodic_zone_bits[ipz]);
        if (zoneVec[izone].getPeriodicBits() == 1) {
          zoneVec[izone].setName("PER1_A");
          zoneVec[izone].setBC("PER1_A");
        }
        else if (zoneVec[izone].getPeriodicBits() == 2) {
          zoneVec[izone].setName("PER1_B");
          zoneVec[izone].setBC("PER1_B");
        }
        else if (zoneVec[izone].getPeriodicBits() == 4) {
          zoneVec[izone].setName("PER2_A");
          zoneVec[izone].setBC("PER2_A");
        }
        else if (zoneVec[izone].getPeriodicBits() == 8) {
          zoneVec[izone].setName("PER2_B");
          zoneVec[izone].setBC("PER2_B");
        }
        else if (zoneVec[izone].getPeriodicBits() == 16) {
          zoneVec[izone].setName("PER3_A");
          zoneVec[izone].setBC("PER3_A");
        }
        else if (zoneVec[izone].getPeriodicBits() == 32) {
          zoneVec[izone].setName("PER3_B");
          zoneVec[izone].setBC("PER3_B");
        }
        cout << "    > \"" << zoneVec[izone].getName() << "\"" << endl;
      }
      delete[] periodic_zone_array;
      delete[] periodic_zone_bits;

    }

  }
  else {

    // no surface was found be we should have offsets for
    // boundary face and x_no data
    assert(nno_pb >= 0);
    assert(nbf >= 0);
    assert(zone_count_check == nbf);

    int * noobf_i = NULL;
    int8 * noobf_v = NULL;
    double (*x_no)[3] = NULL;

    assert(noobf_offset>=0);
    fseek(fp,noobf_offset,SEEK_SET);
    fread(&header,sizeof(Header),1,fp);
    if (byte_swap) ByteSwap::byteSwapHeader(&header,1);
    assert(noobf_i == NULL);
    noobf_i = new int[nbf+1];
    fread(noobf_i+1,sizeof(int),nbf,fp);
    noobf_i[0] = 0;
    for (int ibf = 0; ibf < nbf; ++ibf) noobf_i[ibf+1] += noobf_i[ibf];
    assert(noobf_i[nbf] == noobf_s); // check

    assert(noobf_v == NULL);
    noobf_v = new int8[noobf_s];
    fread(noobf_v,sizeof(int8),noobf_s,fp);
    for (int nob = 0; nob < noobf_s; ++nob) assert((noobf_v[nob] >= 0)&&(noobf_v[nob] < nno_pb)); // check

    assert(xno_offset>=0);
    fseek(fp,xno_offset,SEEK_SET);
    fread(&header,sizeof(Header),1,fp);
    if (byte_swap) ByteSwap::byteSwapHeader(&header,1);
    assert(x_no == NULL);
    x_no = new double[nno_pb][3];
    fread(x_no,sizeof(double),nno_pb*3,fp);

    assert(noobf_i != NULL);
    assert(noobf_v != NULL);
    assert(x_no != NULL);

    // now add to the surface...

    const int ss_nsp0 = nsp;
    nsp += nno_pb+nbf; // one new surface point at the center of each boundary face

    const int ss_nst0 = nst;
    nst += noobf_s; // one new tri per edge link in noobf_i/v

    // allocate new sizes for nsp/nst surface objects
    growNspData(nsp,ss_nsp0);
    growNstData(nst,ss_nst0);

    // set new nodes up to the first boundary node...

    for (int ino = 0; ino < nno_pb; ++ino) {
      FOR_I3 xsp[ss_nsp0+ino][i] = x_no[ino][i];
    }

    // loop on noobf_i/v...

    int nsp_bf = ss_nsp0+nno_pb;
    int ist = ss_nst0;
    int ibf_end = 0;
    assert(zoneVec.empty()); // assume no previous zones
    for (int ii = 0,ii_end=zoneNameAndCount.size(); ii < ii_end; ++ii) {
      const int izone = ii;
      zoneVec.push_back(SurfaceZone(zoneNameAndCount[ii].first));
      const int ibf_start = ibf_end;
      ibf_end += zoneNameAndCount[ii].second;
      for (int ibf = ibf_start; ibf < ibf_end; ++ibf) {
        const int isp_bf = nsp_bf++;
        FOR_I3 xsp[isp_bf][i] = 0.f;
        double wgt = 0.f;
        int ino1 = noobf_v[noobf_i[ibf+1]-1];
        //int ist_prev = ist + noobf_i[ibf+1] - noobf_i[ibf] - 1;
        for (int nob = noobf_i[ibf]; nob != noobf_i[ibf+1]; ++nob) {
          const int ino0 = ino1;
          ino1 = noobf_v[nob];
          const double this_wgt = DIST(x_no[ino0],x_no[ino1]);
          FOR_I3 xsp[isp_bf][i] += this_wgt*(x_no[ino0][i]+x_no[ino1][i]);
          wgt += this_wgt;
          spost[ist][0] = ss_nsp0+ino0;
          spost[ist][1] = ss_nsp0+ino1;
          spost[ist][2] = isp_bf;
          znost[ist] = izone;
          //ist_prev = ist;
          ++ist;
        }
        // normalize xsp[isp_bf]...
        assert(wgt > 0.f);
        FOR_I3 xsp[isp_bf][i] /= 2.f*wgt;
      }
    }
    assert(ibf_end == nbf);
    assert(ist == nst);
    assert(nsp_bf == nsp);

    delete[] x_no;
    delete[] noobf_v;
    delete[] noobf_i;

  }
  fclose(fp);

  cout << " > done read" << endl;

  // build subzone data structures
  szost.resize(nst);  // will persist previous szost values, so only add new nst values
  for (int ist=ss_nst0; ist < nst; ++ist) szost[ist] = znost[ist];

  szozn_i.resize(zoneVec.size()+1);
  szozn_i[0] = 0;
  for (int izn=1,end=zoneVec.size(); izn<=end; ++izn) szozn_i[izn] = szozn_i[izn-1] + 1;
  nsz = szozn_i[zoneVec.size()];

  if (b_foundSurface && nst) return 0;  // surface found
  else return -1;

}


