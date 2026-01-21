#include "LesImageMapper.hpp"

inline bool icase_comp(unsigned char _a, unsigned char _b) {
  return std::tolower(_a) == std::tolower(_b);
}

inline bool equals_ci(const string& s1, const string &s2) {
  if ( s1.length() == s2.length())
    return std::equal(s1.begin(),s1.end(),s2.begin(),icase_comp);
  else
    return false;
}

void LesImageMapper::init(const string &_filename,const string &_filename_2) {

  assert(rfp == NULL);
  assert(rfp_eof == true);
  const int file_err = MiscUtils::openFile(&rfp,_filename,"rb");
  if (file_err != 0) throw(file_err);
  // rfp = fopen(_filename.c_str(),"rb");
  // if (rfp == NULL) {
  //   cerr << "\n\nError: could not open les file: " << _filename << "\n\n" << endl;
  //   throw(-1);
  // }
  rfp_eof = false; // used to indicate that we have hit the eof in rfp...

  rfp_byte_swap = false;
  int itmp[2];
  fread(itmp,sizeof(int),2,rfp);
  if (itmp[0] != UGP_IO_MAGIC_NUMBER) {
    ByteSwap::byteSwap(itmp,2);
    if (itmp[0] != UGP_IO_MAGIC_NUMBER) {
      cerr << "Error: file does not start as expected. aborting."<< endl;
      throw(-1);
    }
    cout << "File requires byte swapping." << endl;
    rfp_byte_swap = true;
  }
  assert(!rfp_byte_swap); // if you hit this, we need to work byte swapping through the routine more carefully
  io_version = itmp[1];
  assert(io_version>=3); //io_versions 3 and 4 need lengthscales added manually...

  rfp_offset = sizeof(int)*2;

  //init value of vvbboxPair
  vvbboxPair = make_pair(0,int8(-1));

  //Grab metadata elements (step,time,hashid...) from restart file header for image metadata
  if (!setRfpHashId()){
    cout << "Warning: restart hash id not found" << endl;
  }

  if (_filename_2 == "") {
    if (!getRfpI0Value("step",step)){
      cout << "Warning: step not found in " << _filename << endl;
    }

    if (!getRfpD0Value("time",time)){
      cout << "Warning: time not found in " << _filename << endl;
    }
  }

  if (io_version==3)
    getNcvGlobalIoV3(); //ncv_global is in its own header, insure it is retrieved.

  if (_filename_2 != "") {
    assert(io_version >= 5);
    assert(sfp == NULL);
    assert(sfp_eof == true);
    const int file_err2 = MiscUtils::openFile(&sfp,_filename_2,"rb");
    if (file_err2 != 0) throw(file_err2);
    // sfp = fopen(_filename_2.c_str(),"rb");
    // if (sfp == NULL) {
    //   cerr << "\n\nError: could not open sles file: " << _filename_2 << "\n\n" << endl;
    //   throw(-1);
    // }
    sfp_eof = false; // used to indicate that we have hit the eof in rfp...

    sfp_byte_swap = false;
    int itmp[2];
    fread(itmp,sizeof(int),2,sfp);
    if (itmp[0] != UGP_IO_MAGIC_NUMBER+1) {
      ByteSwap::byteSwap(itmp,2);
      if (itmp[0] != UGP_IO_MAGIC_NUMBER+1) {
        cerr << "Error: file does not start as expected. aborting."<< endl;
        throw(-1);
      }
      cout << "File requires byte swapping." << endl;
      sfp_byte_swap = true;
    }
    assert(!sfp_byte_swap); // if you hit this, we need to work byte swapping through the routine more carefully
    //io_version = itmp[1];
    //assert(io_version>=3); //io_versions 3 and 4 need lengthscales added manually...

    sfp_offset = sizeof(int)*2;

    //init value of vvbboxPair
    //vvbboxPair = make_pair(0,int8(-1));

    //Grab metadata elements (step,time,hasid...) from solution file header for image metadata
    if (!setSfpHashId()){
      cout << "Warning: sles hash id not found" << endl;
    }
    if (!getSfpI0Value("step",step)){
      cout << "Warning: step not found in " << _filename_2 << endl;
    }

    if (!getSfpD0Value("time",time)){
      cout << "Warning: time not found in " << _filename_2 << endl;
    }
  }
}

double LesImageMapper::getTime(){
  return time;
}

string LesImageMapper::buildCvImageMapFilename(const CtiCanvas * canvas, const float xp[3], const float np[3]) {

  assert(canvas!=NULL);
  //Build unique filename for CvImageMap
    //Includes Geometry xp, np defining the plane
    //Canvas x0, e0, e1, e2 defining the view and
    //Image size defining the projection (as well as e0,e1,e2 being scaled by ni/width)
    //TBD: for completeness a unique identifier for the mesh should be included,
    //     similar to the RestartHashId but should be constant across start/stops and snapshots


  //Compute np_origin, which is the normal vector from origin to the plane,
  //don't just use xp and np since many xp values with a single np could represent the same plane...
  float mag_np = MAG(np); assert(mag_np>0.0);
  float unit_np[3];
  FOR_I3 unit_np[i] = np[i]/mag_np;
  float xp_dot_np = DOT_PRODUCT(xp,unit_np);
  float np_origin[3];
  FOR_I3 np_origin[i] = xp_dot_np*unit_np[i];

  float x0[3];  //simulation coordinates, center of the view
  float e0[3];  //NOTE, not unit vectors, have length ni/width
  float e1[3];
  float e2[3];
  canvas->getX0(x0);
  canvas->getE0(e0);
  canvas->getE1(e1);
  canvas->getE2(e2);
  float ni = canvas->getNi();
  float nj = canvas->getNj();

  //copy memory to unsigned char: np_origin, x0, e0, e1, e2, ni, nj
  //size in bytes:
  int numBytes = (3+3+3+3+3+1+1)*sizeof(float);  //17*4bytes=68bytes
  unsigned char * memPtr = new unsigned char[numBytes];
  memcpy(memPtr+0*sizeof(float), np_origin, 3*sizeof(float));
  memcpy(memPtr+3*sizeof(float), x0, 3*sizeof(float));
  memcpy(memPtr+6*sizeof(float), e0, 3*sizeof(float));
  memcpy(memPtr+9*sizeof(float), e1, 3*sizeof(float));
  memcpy(memPtr+12*sizeof(float), e2, 3*sizeof(float));
  memcpy(memPtr+15*sizeof(float), &ni, 1*sizeof(float));
  memcpy(memPtr+16*sizeof(float), &nj, 1*sizeof(float));

  //FNV Hash to build a filename
  const unsigned int Prime = 0x01000193; //   16777619
  const unsigned int Seed  = 0x811C9DC5; // 2166136261
  unsigned int hash = Seed;
  unsigned char * tempPtr = memPtr;

  while (numBytes--)
    hash = (*tempPtr++ ^ hash) * Prime;

  delete[] memPtr;

  std::stringstream ss;
  ss << "cims/cim_" << std::setw(10) << std::setfill('0') << hash << ".bin";
  MiscUtils::mkdir_for_file(ss.str());  //create folder
  return ss.str();
}

void LesImageMapper::getPointsInsideGeom(vector<SimplePoint>& pointVec,SimpleGeom * geom) {

  double geom_bbox[6];
  geom->getBoundingBox(geom_bbox);

  int8 offset = getRfpCvD2Offset("x_vv");
  if (offset == -1) {
    cerr << "Cannot find x_vv data" << endl;
    assert(0);
  }
  assert(got_ncv_global);
  my_fseek(rfp,offset);

  const int read_size = 10000000; //chunk size left in place in case bounding box data not present
  double (*buf)[3] = new double[min(int(ncv_global),read_size)][3];

  int nbb;
  int8 offset_bb = getRfpVVBBoxOffset("vv_bbox", nbb);
  if (offset_bb==-1) { //vv_bbox was not found...just read all x_vv data
    cout << "vv_bbox not found, chunking through all nodal data" << endl;
    nbb = 1;
  }

  cout << "reading x_vv..";
  int8 index0 = 0;
  for (int ibb = 0; ibb<nbb; ++ibb){
    //int8 this_offset_xvv = offset_xvv;
    //int8 this_offset_rvv = offset_rvv;
    int ncv_bbox = ncv_global; //will be reasssigned in checkVVBBoxIntersection
    if ((offset_bb < 0) || checkVVBBoxBBoxIntersection(ncv_bbox,ibb,geom_bbox)) { //check offset_bb in case there is no vv_bbox record
      int index0_bbox = 0;
      while (index0_bbox < ncv_bbox){
        cout << "."<<ibb<<".";
        cout.flush();
        const int this_read_size = min(read_size,ncv_bbox-index0_bbox);
	my_fseek(rfp,offset+(index0+index0_bbox)*sizeof(double)*3);
        fread(buf,sizeof(double),this_read_size*3,rfp);
	for (int icv = 0; icv < this_read_size; ++icv) {
	  if (geom->pointIsInside(buf[icv])) {
	    pointVec.push_back(SimplePoint(buf[icv]));
	  }
	} //for (int icv = 0; icv < this_read_size; ++icv)
        index0_bbox += this_read_size;
      } //while (index0_bbox < ncv_bbox)
      assert(index0_bbox == ncv_bbox);
    } // if (...
    index0 += ncv_bbox;
  } // for (int ibb = 0; ibb<nbb; ++ibb)
  cout << "done" << endl;

  delete[] buf;

}

void LesImageMapper::buildCvPlaneMap(CvPlaneMap& cvPlaneMap, const string &lengthscaleName) {

  const int read_size = 10000000; //chunk size left in place in case bounding box data not present

  //r_vv is required for imaging
  double * r_vv=NULL;
  int8 offset_rvv = getRfpCvD1Offset(lengthscaleName);
  if (offset_rvv == -1)
    CERR(lengthscaleName << " lengthscale not found, use Charles POST to add it to your restart file");
  r_vv = new double[min(int(ncv_global),read_size)];

  int8 offset = getRfpCvD2Offset("x_vv");
  if (offset == -1) {
    cerr << "Cannot find x_vv data" << endl;
    assert(0);
  }
  assert(got_ncv_global);
  my_fseek(rfp,offset);

  double (*buf)[3] = new double[min(int(ncv_global),read_size)][3];

  int nbb;
  int8 offset_bb = getRfpVVBBoxOffset("vv_bbox", nbb);
  if (offset_bb==-1) { //vv_bbox was not found...just read all x_vv data
    cout << "vv_bbox not found, chunking through all nodal data" << endl;
    nbb = 1;
  }

  cout << "reading x_vv for CvImageMap";
  int8 index0 = 0;
  for (int ibb = 0; ibb<nbb; ++ibb){
    //int8 this_offset_xvv = offset_xvv;
    //int8 this_offset_rvv = offset_rvv;
    int ncv_bbox = ncv_global; //will be reasssigned in checkVVBBoxIntersection
    if (offset_bb < 0 || checkVVBBoxIntersection(ncv_bbox,ibb,cvPlaneMap.xp,cvPlaneMap.np)){ //check offset_bb in case there is no vv_bbox record
      int index0_bbox = 0;
      while (index0_bbox < ncv_bbox){
        cout << "."<<ibb<<".";
        cout.flush();
        const int this_read_size = min(read_size,ncv_bbox-index0_bbox);

        my_fseek(rfp,offset_rvv+(index0+index0_bbox)*sizeof(double));
        fread(r_vv,sizeof(double),this_read_size,rfp);
        my_fseek(rfp,offset+(index0+index0_bbox)*sizeof(double)*3);
        fread(buf,sizeof(double),this_read_size*3,rfp);

        for (int icv = 0; icv < this_read_size; ++icv) {
          // calculate the normal distance to the plane.
          const double dist =
            (buf[icv][0]-cvPlaneMap.xp[0])*cvPlaneMap.np[0] +
            (buf[icv][1]-cvPlaneMap.xp[1])*cvPlaneMap.np[1] +
            (buf[icv][2]-cvPlaneMap.xp[2])*cvPlaneMap.np[2];
          // const double dist_c = dist*maskCanvas->getNi()/maskCanvas->getWidth();
          const double rs2 = r_vv[icv]*r_vv[icv] - dist*dist;
          if (rs2 > 0.0) {
            cvPlaneMap.push_back(index0 + index0_bbox + int8(icv),buf[icv],r_vv[icv]);
          }
        }
        index0_bbox += this_read_size;
      } //while (index0_bbox < ncv_bbox)
      assert(index0_bbox == ncv_bbox);
    } // if (checkVVBBoxIntersectionAndSetOffsets(offset_xvv,offset_rvv,ncv_bbox,xp,np))
    index0 += ncv_bbox;
  }
  if (r_vv!=NULL)
    delete[] r_vv;
  delete[] buf;
}



int LesImageMapper::buildCvImageMap(CvImageMap& cim, const CtiCanvas * maskCanvas, CvPlaneMap& cvPlaneMap, const string &lengthscaleName) {

  assert(maskCanvas!=NULL);
  using namespace ByteSwap;

  //Work in simulation coordinates
  //1. Test x_vv proximity to plane with r_vv
  //2. If point is near plane, determine it's influence on the canvas
  //3. Convert to canvas coordinates using canvas x0,e0,e1,e2
  //4. Map canvas pixels on to plane, compute distance from CV
  //5. Try to limit the number of pixels iterated over per CV
  //    circular region of influce on plane project to an ellipse on the 2d canvas


  //Canvas vectors
  float x0[3];  //simulation coordinates
  float e0[3];  //NOTE, not unit vectors, have length ni/width
  float e1[3];
  float e2[3];
  maskCanvas->getX0(x0);
  maskCanvas->getE0(e0);
  maskCanvas->getE1(e1);
  maskCanvas->getE2(e2);

  //plane center and normal in canvas coordinates
  //use this later for depth testing CV-to plane pixel location distance
  float xp_c[3] = { ( (cvPlaneMap.xp[0]-x0[0])*e0[0] + (cvPlaneMap.xp[1]-x0[1])*e0[1] + (cvPlaneMap.xp[2]-x0[2])*e0[2] ),
                    ( (cvPlaneMap.xp[0]-x0[0])*e1[0] + (cvPlaneMap.xp[1]-x0[1])*e1[1] + (cvPlaneMap.xp[2]-x0[2])*e1[2] ),
                    ( (cvPlaneMap.xp[0]-x0[0])*e2[0] + (cvPlaneMap.xp[1]-x0[1])*e2[1] + (cvPlaneMap.xp[2]-x0[2])*e2[2] )};
  float unit_np_c[3] = { ( (cvPlaneMap.np[0])*e0[0] + (cvPlaneMap.np[1])*e0[1] + (cvPlaneMap.np[2])*e0[2] ),
                       ( (cvPlaneMap.np[0])*e1[0] + (cvPlaneMap.np[1])*e1[1] + (cvPlaneMap.np[2])*e1[2] ),
                       ( (cvPlaneMap.np[0])*e2[0] + (cvPlaneMap.np[1])*e2[1] + (cvPlaneMap.np[2])*e2[2] )};
  float np_c_mag = MAG(unit_np_c);
  FOR_I3 unit_np_c[i] /= np_c_mag;

  //canvas ei, ej unit vectors projected onto data plane
  //use these later to try to limit canvas pixels iterated on per cv
  double ei_proj_unit_c[3] = {-unit_np_c[2],0.0,unit_np_c[0]};
  double ei_proj_unit_c_mag = MAG(ei_proj_unit_c);
  FOR_I3 ei_proj_unit_c[i] /= ei_proj_unit_c_mag;
  double ej_proj_unit_c[3] = {0.0,-unit_np_c[2],unit_np_c[1]};
  double ej_proj_unit_c_mag = MAG(ej_proj_unit_c);
  FOR_I3 ej_proj_unit_c[i] /= ej_proj_unit_c_mag;

//  FOR_I3{
//   cout << " unit_np_c" << i << " " << unit_np_c[i] << " ei_proj_unit_c" << i << " " << ei_proj_unit_c[i] << " ej_proj_unit_c" << i << " " << ej_proj_unit_c[i] <<endl;
//  }

  const int image_size = maskCanvas->getNi()*maskCanvas->getNj();

  int8 * cv_buf = new int8[image_size]; // holds the global icv associated with the closest cv, otherwise -1 or -2
  if (d2_buf){
    for (int image_ij = 0; image_ij < image_size; ++image_ij){
      d2_buf[image_ij] = 0.0f;
    }
  }else{
    d2_buf = new float[image_size]; // holds the dist2 of the closest cv, if there is one
  }

  int dap_count = 0;
  for (int image_ij = 0; image_ij < image_size; ++image_ij){
    if (maskCanvas->isMaskPixel(image_ij)){  //apply mask from canvas
      cv_buf[image_ij] = -1;
      ++dap_count;
    }
    else
      cv_buf[image_ij] = -2;
  }


  if (!cvPlaneMap.checkMapFlag()){
    COUT1("buildCvPlaneMap()");
    buildCvPlaneMap(cvPlaneMap,lengthscaleName);
    cout << "done" << endl;
    cvPlaneMap.setMapFlag();
  }
  else{
   COUT1("using existing cvPlaneMap()");
  }

  for (int icv = 0; icv < cvPlaneMap.size(); ++icv){

    //project x_vv onto the plane, this is the center of circle of influence cv has on the plane
    //const double xc_s[3] = {buf[icv][0] - dist*unit_np[0],
    //                        buf[icv][1] - dist*unit_np[1],
    //                        buf[icv][2] - dist*unit_np[2]};
    double xc_s[3];
    cvPlaneMap.getXcs(xc_s,icv);

    // figure out a range of pixels that this point could affect...
    // xc,yc,zc the center of the circle of influence on the plane in canvas coordinates
    const double xc = ( (xc_s[0]-x0[0])*e0[0] + (xc_s[1]-x0[1])*e0[1] + (xc_s[2]-x0[2])*e0[2] );
    const double yc = ( (xc_s[0]-x0[0])*e1[0] + (xc_s[1]-x0[1])*e1[1] + (xc_s[2]-x0[2])*e1[2] );
    // const double zc = ( (xc_s[0]-x0[0])*e2[0] + (xc_s[1]-x0[1])*e2[1] + (xc_s[2]-x0[2])*e2[2] );

    // and the radius of the circle of influence in canvas dimensions...
    const double rc = sqrt(cvPlaneMap.getRs2(icv))*maskCanvas->getNi()/maskCanvas->getWidth();

    const double i_max = xc+rc;
    const double i_min = xc-rc;
    const double j_max = yc+rc;
    const double j_min = yc-rc;

    // test if projected circle falls on canvas
    if ((i_max >= 0)&&(i_min < maskCanvas->getNi())&&(j_max >= 0)&&(j_min < maskCanvas->getNj())) {
      //for all rows of pixels between j_min and j_max, loop over some columns
      for (int jp = max(0,int(ceil(j_min))); jp <= min(maskCanvas->getNj()-1,int(floor(j_max))); ++jp) {
        const double dyc = (double(jp) - yc);
        const double i_min_local = xc-sqrt(rc*rc-dyc*dyc);
        const double i_max_local = xc+sqrt(rc*rc-dyc*dyc);

        for (int ip = max(0,int(ceil(i_min_local))); ip <= min(maskCanvas->getNi()-1,int(floor(i_max_local)));  ++ip) {
          const int index = jp*maskCanvas->getNi()+ip;
          //if cv_buf we are in the data mask, compute the cv-plane distance
          //if the distance is less than already set store it
          if (cv_buf[index] != -2){
            //build cv to plane pixel depth using components already calculated...
            //normal distance from cv to plane dist_c plus in place distance to current pixel mapped to 3d plane location

            //    long form depth check, distance from current pixel to cv in canvas units
            //plane depth coordinate for pixel ip, jp
            double kp = xp_c[2] - (unit_np_c[0]*(ip-xp_c[0]) +
                                   unit_np_c[1]*(jp-xp_c[1]) )
                                  /unit_np_c[2];
            //cv location in canvas coordinates
            double x_vv[3];
            cvPlaneMap.getXvv(x_vv,icv);

            const double xcv_c = ( (x_vv[0]-x0[0])*e0[0] + (x_vv[1]-x0[1])*e0[1] + (x_vv[2]-x0[2])*e0[2] );
            const double ycv_c = ( (x_vv[0]-x0[0])*e1[0] + (x_vv[1]-x0[1])*e1[1] + (x_vv[2]-x0[2])*e1[2] );
            const double zcv_c = ( (x_vv[0]-x0[0])*e2[0] + (x_vv[1]-x0[1])*e2[1] + (x_vv[2]-x0[2])*e2[2] );
            double this_d2 = (ip-xcv_c)*(ip-xcv_c)+(jp-ycv_c)*(jp-ycv_c)+(kp-zcv_c)*(kp-zcv_c);


            if ((cv_buf[index] == -1)||((cv_buf[index] >= 0)&&(this_d2 < d2_buf[index]))) {
              cv_buf[index] = cvPlaneMap.getIcvGlobal(icv);
              d2_buf[index] = (float) this_d2;
            }
          }
        }
      }
    }
  }
  cout << "building data-to-pixel map..." << endl;

  // =====================================================================
  // now organize the pixel info in terms of ordered cv data...
  // =====================================================================

  {

    int active_image_size = 0;
    for (int index = 0; index < image_size; ++index)
      if (cv_buf[index] >= 0)
	++active_image_size;

    vector< pair<int8,unsigned int> > cvIndexPair(active_image_size);
    active_image_size = 0;
    for (int index = 0; index < image_size; ++index) {
      if (cv_buf[index] >= 0) {
	cvIndexPair[active_image_size].first = cv_buf[index];
	cvIndexPair[active_image_size].second = static_cast<unsigned int>(index);
	++active_image_size;
      }
    }

    sort(cvIndexPair.begin(),cvIndexPair.end());

    // try to see if shift works
    //unsigned int test = ((unsigned int)((1<<CIM_REPEAT_BITS)-1))<<(32-CIM_REPEAT_BITS);
    //cout << "test: " << test << endl;

    for (int iter = 0; iter < 2; ++iter) {
      cim.ncv_unique = 0;
      int8 first_prev = -1;
      unsigned int second_prev = 0;
      unsigned int repeat_count = 0;
      int poucv_s = 0;
      for (int ii = 0,ii_max=cvIndexPair.size(); ii < ii_max; ++ii) {
	if (cvIndexPair[ii].first != first_prev) {
	  // this is a new value...
	  assert(cvIndexPair[ii].first > first_prev);
	  if (iter == 1) {
	    // first, complete the last entry before moving on...
	    if (repeat_count > 0) {
	      cim.poucv_v[poucv_s-1] |= (repeat_count<<(32-CIM_REPEAT_BITS));
	      assert((cim.poucv_v[poucv_s-1]>>(32-CIM_REPEAT_BITS)) == repeat_count); // check extraction
	    }
	    // set the next entry details...
	    cim.ucv_v[cim.ncv_unique] = cvIndexPair[ii].first;
	    cim.poucv_i[cim.ncv_unique] = poucv_s;
	    cim.poucv_v[poucv_s] = cvIndexPair[ii].second; // the first pixel index
	  }
	  first_prev = cvIndexPair[ii].first;
	  second_prev = cvIndexPair[ii].second;
	  ++cim.ncv_unique;
	  ++poucv_s;
	  repeat_count = 0;
	}
	else if (second_prev+1 == cvIndexPair[ii].second) {
	  ++second_prev;
	  if (repeat_count == (1<<CIM_REPEAT_BITS)-1) {
	    // we hit the maximum repeat count. can't do this, so start a new poucv_s...
	    if (iter == 1) {
	      cim.poucv_v[poucv_s-1] |= (repeat_count<<(32-CIM_REPEAT_BITS));
	      assert((cim.poucv_v[poucv_s-1]>>(32-CIM_REPEAT_BITS)) == repeat_count); // check extraction
	      cim.poucv_v[poucv_s] = second_prev;
	    }
	    //cout << "image looped" << endl;
	    ++poucv_s;
	    repeat_count = 0;
	  }
	  else {
	    ++repeat_count;
	  }
	}
	else {
	  // we hit a pixel that is not contiguous...
	  if (iter == 1) {
	    if (repeat_count > 0) {
	      cim.poucv_v[poucv_s-1] |= (repeat_count<<(32-CIM_REPEAT_BITS));
	      assert((cim.poucv_v[poucv_s-1]>>(32-CIM_REPEAT_BITS)) == repeat_count); // check extraction
	    }
	    cim.poucv_v[poucv_s] = cvIndexPair[ii].second; // the first pixel index
	  }
	  second_prev = cvIndexPair[ii].second;
	  ++poucv_s;
	  repeat_count = 0;
	}
      }
      if (iter == 0) {
	cout << "active_image_size: " << active_image_size << " poucv_s: " << poucv_s << " ncv_unique: " << cim.ncv_unique << endl;
	cim.ucv_v = new int8[cim.ncv_unique];
	cim.poucv_i = new int[cim.ncv_unique+1];
	cim.poucv_v = new unsigned int[poucv_s];
      }
      else {
	if (repeat_count > 0) {
	  // on the second time through, complete the last entry if it repeats...
	  cim.poucv_v[poucv_s-1] |= (repeat_count<<(32-CIM_REPEAT_BITS));
	  assert((cim.poucv_v[poucv_s-1]>>(32-CIM_REPEAT_BITS)) == repeat_count); // check extraction
	}
	cim.poucv_i[cim.ncv_unique] = poucv_s;
      }
    }

    // ====================================================
    // now check the compressed map...
    // TODO skip this check once things are working?
    // ====================================================

    int8 * cv_buf_check = new int8[image_size];
    for (int i = 0; i < image_size; ++i){
      cv_buf_check[i] = -1;
    }

    for (int icv = 0; icv < cim.ncv_unique; ++icv) {
      const int8 icv_global = cim.ucv_v[icv];
      for (int poc = cim.poucv_i[icv]; poc != cim.poucv_i[icv+1]; ++poc) {
	unsigned int index = cim.poucv_v[poc] & ((1<<(32-CIM_REPEAT_BITS))-1);
	unsigned int repeat_count = cim.poucv_v[poc]>>(32-CIM_REPEAT_BITS);
	assert( (index | (repeat_count<<(32-CIM_REPEAT_BITS))) == cim.poucv_v[poc]);
	for (unsigned int offset = 0; offset <= repeat_count; ++offset) {
	  cv_buf_check[index+offset] = icv_global;
	}
      }
    }

    for (int i = 0; i < image_size; ++i) {
      //cout << "i,cv_buf[i],cv_buf_check[i]: " << i << " " << cv_buf[i] << " " << cv_buf_check[i] << endl;
      //getchar();
      if (cv_buf[i] >= 0) {
	assert(cv_buf[i] == cv_buf_check[i]);
      }
    }

    delete[] cv_buf_check;

  }

  delete[] cv_buf;

  return 1;
}

string LesImageMapper::buildCvImageMapFilenameSurface(const CtiCanvas * canvas, const float bbox[6], const int &nst, const std::vector< PlaneData<float> > &vBlankPlaneData, const IntFlag &show_sz_flag) {

  assert(canvas!=NULL);
  //Build unique filename for CvImageMap
    //Include surface identifiers: nst and bbox
    //Canvas x0, e0, e1, e2 defining the view and
    //Image size defining the projection (as well as e0,e1,e2 being scaled by ni/width)
    //blanking and hidden zones also impact displayed surfer

  float x0[3];  //simulation coordinates, center of the view
  float e0[3];  //NOTE, not unit vectors, have length ni/width
  float e1[3];
  float e2[3];
  canvas->getX0(x0);
  canvas->getE0(e0);
  canvas->getE1(e1);
  canvas->getE2(e2);
  float ni = canvas->getNi();
  float nj = canvas->getNj();

  int bitFlagLength = show_sz_flag.getLength()/8;
  if (show_sz_flag.getLength()%8!=0)
    ++bitFlagLength;

  //copy memory to unsigned char: nst, bbox, x0, e0, e1, e2, ni, nj
  //size in bytes:
  int numBytes = 1*sizeof(int) + (6+3+3+3+3+1+1+vBlankPlaneData.size()*3)*sizeof(float) + bitFlagLength*sizeof(unsigned char);
  unsigned char * memPtr = new unsigned char[numBytes];
  unsigned int offset = 0;
  memcpy(memPtr+offset, &nst, 1*sizeof(int));
  offset+=1*sizeof(int);
  memcpy(memPtr+offset, bbox, 6*sizeof(float));
  offset+=6*sizeof(float);
  memcpy(memPtr+offset, x0, 3*sizeof(float));
  offset+=3*sizeof(float);
  memcpy(memPtr+offset, e0, 3*sizeof(float));
  offset+=3*sizeof(float);
  memcpy(memPtr+offset, e1, 3*sizeof(float));
  offset+=3*sizeof(float);
  memcpy(memPtr+offset, e2, 3*sizeof(float));
  offset+=3*sizeof(float);
  memcpy(memPtr+offset, &ni, 1*sizeof(float));
  offset+=1*sizeof(float);
  memcpy(memPtr+offset, &nj, 1*sizeof(float));
  offset+=1*sizeof(float);

  //blanking planes...
  for (int ip=0,ip_max=vBlankPlaneData.size(); ip<ip_max; ip++){
    //Compute np_origin, which is the normal vector from origin to the blanking plane,
    //don't just use center and normal since many xp values with a single np could represent the same plane...
    float mag_np = MAG(vBlankPlaneData[ip].normal); assert(mag_np>0.0);
    float unit_np[3];
    FOR_I3 unit_np[i] = vBlankPlaneData[ip].normal[i]/mag_np;
    float xp_dot_np = DOT_PRODUCT(vBlankPlaneData[ip].center,unit_np);
    float np_origin[3];
    FOR_I3 np_origin[i] = xp_dot_np*unit_np[i];
    memcpy(memPtr+offset,np_origin,3*sizeof(float));
    offset+=3*sizeof(float);
  }

  //hidden zones, build a flag with bits set
  //for hidden zones
  unsigned char * bitFlag = new unsigned char[bitFlagLength];
  int iz = 0;
  for (int ibf = 0; ibf<bitFlagLength; ++ibf){
    bitFlag[ibf] = '\0';
    for (int ib = 0; ib<8; ++ib){
      if (iz<show_sz_flag.getLength()){
        if (show_sz_flag[iz])
          bitFlag[ibf] |= (1 << ib);
      }
      else break;

      ++iz;
    }
  }
  memcpy(memPtr+offset,bitFlag,bitFlagLength*sizeof(unsigned char));
  delete[] bitFlag;
  offset+=bitFlagLength*sizeof(unsigned char);
  assert(int(offset)==numBytes);

  //FNV Hash to build a filename
  const unsigned int Prime = 0x01000193; //   16777619
  const unsigned int Seed  = 0x811C9DC5; // 2166136261
  unsigned int hash = Seed;
  unsigned char * tempPtr = memPtr;

  while (numBytes--)
    hash = (*tempPtr++ ^ hash) * Prime;

  delete[] memPtr;

  std::stringstream ss;
  ss << "cims/cimS_" << std::setw(10) << std::setfill('0') << hash << ".bin";
  MiscUtils::mkdir_for_file(ss.str());  //create folder
  return ss.str();
}

int LesImageMapper::buildCvImageMapSurface(CvImageMap& cim, const CtiCanvas * _canvas, const IntFlag &bfzone_flag, const string &lengthscaleName) {

  assert(_canvas!=NULL);
  using namespace ByteSwap;

  //Canvas vectors
  float x0[3];  //simulation coordinates
  float e0[3];  //NOTE, not unit vectors, have length ni/width
  float e1[3];
  float e2[3];
  _canvas->getX0(x0);
  _canvas->getE0(e0);
  _canvas->getE1(e1);
  _canvas->getE2(e2);

  const int image_size = _canvas->getNi()*_canvas->getNj();

  int8 * cv_buf = new int8[image_size]; // holds the global icv associated with the closest cv, otherwise -1 or -2
  if (d2_buf){
    for (int image_ij = 0; image_ij < image_size; ++image_ij){
      d2_buf[image_ij] = 0.0f;
    }
  }else{
    d2_buf = new float[image_size]; // holds the dist2 of the closest cv, if there is one
  }


  for (int image_ij = 0; image_ij < image_size; ++image_ij){
    if (_canvas->isSurfacePixel(image_ij)){  //apply mask from canvas
      cv_buf[image_ij] = -1;
    }
    else
      cv_buf[image_ij] = -2;
  }

  const int read_size = 1000000; //chunk size left in place in case bounding box data not present

  int8 * cvobf=NULL;
  int8 offset_cvobf = getRfpCvobfOffsetAndNbfGlobal();
  if (offset_cvobf==-1)
    CERR("Unable to build surface data map, can't find cvobf_global");

  int8 offset_rvv = getRfpCvD1Offset(lengthscaleName);
  if (offset_rvv == -1){
    CERR(lengthscaleName << " lengthscale not found, use Charles POST to add it to your restart file");
    assert(0);
  }
  int8 offset_xvv = getRfpCvD2Offset("x_vv");
  if (offset_xvv == -1) {
    cerr << "Cannot find x_vv data" << endl;
    assert(0);
  }
  assert(got_ncv_global);

  cout << "reading cvobf";  //TODO DAP - read only zones flagged in bfzone_flag
  if (getRfpNbfZnOffsetVec()){
    for (int izn = 0,nzn=nbfZnOffsetVec.size(); izn < nzn; ++izn){
      if (bfzone_flag[izn]==1){
        cout << "; zone " << izn;
        int8 index0 = 0;
        int8 nbf_zn_global = nbfZnOffsetVec[izn].first;
        while (index0 < nbf_zn_global){
          cvobf = new int8[min(int(nbf_zn_global),read_size)];

          cout <<".";
          cout.flush();
          const int this_read_size = min(read_size,int(nbf_zn_global-index0));

          my_fseek(rfp,offset_cvobf+(nbfZnOffsetVec[izn].second+index0)*sizeof(int8));
          fread(cvobf,sizeof(int8),this_read_size,rfp);

          set<int8> cv_list; // build unique, ordered list of cv's
          for (int ibf = 0; ibf<this_read_size; ++ibf){
            int8 icv = cvobf[ibf];
//            if (cv_list.find(icv) == cv_list.end()) {
            cv_list.insert(icv);  //set insertion only accepts unique values by default...
//            }
          }
          delete[] cvobf;

          //sparse read x_vv and r_vv based on cv_list
          const int max_buf_size = 1000000;
          double * r_vv = new double[max_buf_size];
          double (*x_vv)[3] = new double[max_buf_size][3];

          set<int8>::iterator it_end=cv_list.begin();
          while (it_end != cv_list.end()) {
            set<int8>::iterator it_begin = it_end++;
            while ((it_end != cv_list.end())&&( (*it_end)-(*it_begin) < max_buf_size) ){
              ++it_end;
            }

            set<int8>::iterator it_end_m1 = it_end;
            --it_end_m1;

            //int ncv_sparse = (*(&(*it_end)-1))-(*it_begin)+1;
            int ncv_sparse = (*it_end_m1)-(*it_begin)+1;

            // seek to it_begin...
            my_fseek(rfp,offset_rvv+(*it_begin)*sizeof(double));
            fread(r_vv,sizeof(double),ncv_sparse,rfp);
            my_fseek(rfp,offset_xvv+(*it_begin)*sizeof(double)*3);
            fread(x_vv,sizeof(double),ncv_sparse*3,rfp);

            // we've read chunks of x_vv and r_vv, next pick
            // only the values from our cv_list and build map
            for (set<int8>::iterator it = it_begin; it != it_end; ++it){
              int icv_sparse = (*it)-(*it_begin);

              // figure out a range of pixels that this point could affect...
              // x_vv_c == x_vv in canvas coordinates
              double x_vv_c[3];
              x_vv_c[0] = ( (x_vv[icv_sparse][0]-x0[0])*e0[0] + (x_vv[icv_sparse][1]-x0[1])*e0[1] + (x_vv[icv_sparse][2]-x0[2])*e0[2] );
              x_vv_c[1] = ( (x_vv[icv_sparse][0]-x0[0])*e1[0] + (x_vv[icv_sparse][1]-x0[1])*e1[1] + (x_vv[icv_sparse][2]-x0[2])*e1[2] );
              x_vv_c[2] = ( (x_vv[icv_sparse][0]-x0[0])*e2[0] + (x_vv[icv_sparse][1]-x0[1])*e2[1] + (x_vv[icv_sparse][2]-x0[2])*e2[2] );

              // and the radius of the sphere of influence in canvas dimensions...
              const double r_vv_c = r_vv[icv_sparse]*_canvas->getNi()/_canvas->getWidth();

              // max canvas extents of projected circle of influence
              const double i_max = x_vv_c[0]+r_vv_c;
              const double i_min = x_vv_c[0]-r_vv_c;
              const double j_max = x_vv_c[1]+r_vv_c;
              const double j_min = x_vv_c[1]-r_vv_c;

              // test if projected circle falls on canvas
              if ((i_max >= 0)&&(i_min < _canvas->getNi())&&(j_max >= 0)&&(j_min < _canvas->getNj())) {
                //for all rows of pixels between j_min and j_max, loop over some columns
                for (int jp = max(0,int(ceil(j_min))); jp <= min(_canvas->getNj()-1,int(floor(j_max))); ++jp) {
                  const double dyc = (double(jp) - x_vv_c[1]);
                  const double i_min_local = x_vv_c[0]-sqrt(r_vv_c*r_vv_c-dyc*dyc);
                  const double i_max_local = x_vv_c[0]+sqrt(r_vv_c*r_vv_c-dyc*dyc);

                  for (int ip = max(0,int(ceil(i_min_local))); ip <= min(_canvas->getNi()-1,int(floor(i_max_local)));  ++ip) {
                    const int index = jp*_canvas->getNi()+ip;
                    //if cv_buf we are in the data mask, compute the cv-surface pixel distance
                    //if the distance is less than already set store it
                    if (cv_buf[index] != -2){

                      //distance from current surface pixel to cv in canvas units
                      float kp = _canvas->getDepth(ip,jp);

                      //cv location in canvas coordinates
                      double this_d2 = (ip-x_vv_c[0])*(ip-x_vv_c[0])+(jp-x_vv_c[1])*(jp-x_vv_c[1])+(kp-x_vv_c[2])*(kp-x_vv_c[2]);

                      if ((cv_buf[index] == -1)||((cv_buf[index] >= 0)&&(this_d2 < d2_buf[index]))) {
                        cv_buf[index] = *it;
                        d2_buf[index] = (float) this_d2;
                      }
                    }
                  }
                }
              }
            }
          } //while (it_end != cv_list.end())
          index0 += this_read_size;
          delete[] x_vv;
          delete[] r_vv;
        } //while(index0 < nbf_zn_global)
      } //if (bfzone_flag[izn]==1)
    } //for (int izn = 0; izn < nbfZnOffsetVec.size(); ++izn){
  } //if (getRfpNbfZnOffsetVec()){
  cout << "done" << endl;
  cout << "building data-to-pixel map..." << endl;

  // =====================================================================
  // now organize the pixel info in terms of ordered cv data...
  // =====================================================================

  {

    int active_image_size = 0;
    for (int index = 0; index < image_size; ++index)
      if (cv_buf[index] >= 0)
	++active_image_size;

    vector< pair<int8,unsigned int> > cvIndexPair(active_image_size);
    active_image_size = 0;
    for (int index = 0; index < image_size; ++index) {
      if (cv_buf[index] >= 0) {
	cvIndexPair[active_image_size].first = cv_buf[index];
	cvIndexPair[active_image_size].second = static_cast<unsigned int>(index);
	++active_image_size;
      }
    }

    sort(cvIndexPair.begin(),cvIndexPair.end());

    // try to see if shift works
    //unsigned int test = ((unsigned int)((1<<CIM_REPEAT_BITS)-1))<<(32-CIM_REPEAT_BITS);
    //cout << "test: " << test << endl;

    for (int iter = 0; iter < 2; ++iter) {
      cim.ncv_unique = 0;
      int8 first_prev = -1;
      unsigned int second_prev = 0;
      unsigned int repeat_count = 0;
      int poucv_s = 0;
      for (int ii = 0,ii_max=cvIndexPair.size(); ii < ii_max; ++ii) {
	if (cvIndexPair[ii].first != first_prev) {
	  // this is a new value...
	  assert(cvIndexPair[ii].first > first_prev);
	  if (iter == 1) {
	    // first, complete the last entry before moving on...
	    if (repeat_count > 0) {
	      cim.poucv_v[poucv_s-1] |= (repeat_count<<(32-CIM_REPEAT_BITS));
	      assert((cim.poucv_v[poucv_s-1]>>(32-CIM_REPEAT_BITS)) == repeat_count); // check extraction
	    }
	    // set the next entry details...
	    cim.ucv_v[cim.ncv_unique] = cvIndexPair[ii].first;
	    cim.poucv_i[cim.ncv_unique] = poucv_s;
	    cim.poucv_v[poucv_s] = cvIndexPair[ii].second; // the first pixel index
	  }
	  first_prev = cvIndexPair[ii].first;
	  second_prev = cvIndexPair[ii].second;
	  ++cim.ncv_unique;
	  ++poucv_s;
	  repeat_count = 0;
	}
	else if (second_prev+1 == cvIndexPair[ii].second) {
	  ++second_prev;
	  if (repeat_count == (1<<CIM_REPEAT_BITS)-1) {
	    // we hit the maximum repeat count. can't do this, so start a new poucv_s...
	    if (iter == 1) {
	      cim.poucv_v[poucv_s-1] |= (repeat_count<<(32-CIM_REPEAT_BITS));
	      assert((cim.poucv_v[poucv_s-1]>>(32-CIM_REPEAT_BITS)) == repeat_count); // check extraction
	      cim.poucv_v[poucv_s] = second_prev;
	    }
	    //cout << "image looped" << endl;
	    ++poucv_s;
	    repeat_count = 0;
	  }
	  else {
	    ++repeat_count;
	  }
	}
	else {
	  // we hit a pixel that is not contiguous...
	  if (iter == 1) {
	    if (repeat_count > 0) {
	      cim.poucv_v[poucv_s-1] |= (repeat_count<<(32-CIM_REPEAT_BITS));
	      assert((cim.poucv_v[poucv_s-1]>>(32-CIM_REPEAT_BITS)) == repeat_count); // check extraction
	    }
	    cim.poucv_v[poucv_s] = cvIndexPair[ii].second; // the first pixel index
	  }
	  second_prev = cvIndexPair[ii].second;
	  ++poucv_s;
	  repeat_count = 0;
	}
      }
      if (iter == 0) {
	cout << "active_image_size: " << active_image_size << " poucv_s: " << poucv_s << " ncv_unique: " << cim.ncv_unique << endl;
	cim.ucv_v = new int8[cim.ncv_unique];
	cim.poucv_i = new int[cim.ncv_unique+1];
	cim.poucv_v = new unsigned int[poucv_s];
      }
      else {
	if (repeat_count > 0) {
	  // on the second time through, complete the last entry if it repeats...
	  cim.poucv_v[poucv_s-1] |= (repeat_count<<(32-CIM_REPEAT_BITS));
	  assert((cim.poucv_v[poucv_s-1]>>(32-CIM_REPEAT_BITS)) == repeat_count); // check extraction
	}
	cim.poucv_i[cim.ncv_unique] = poucv_s;
      }
    }

    // ====================================================
    // now check the compressed map...
    // TODO skip this check once things are working?
    // ====================================================

    int8 * cv_buf_check = new int8[image_size];
    for (int i = 0; i < image_size; ++i){
      cv_buf_check[i] = -1;
    }

    for (int icv = 0; icv < cim.ncv_unique; ++icv) {
      const int8 icv_global = cim.ucv_v[icv];
      for (int poc = cim.poucv_i[icv]; poc != cim.poucv_i[icv+1]; ++poc) {
	unsigned int index = cim.poucv_v[poc] & ((1<<(32-CIM_REPEAT_BITS))-1);
	unsigned int repeat_count = cim.poucv_v[poc]>>(32-CIM_REPEAT_BITS);
	assert( (index | (repeat_count<<(32-CIM_REPEAT_BITS))) == cim.poucv_v[poc]);
	for (unsigned int offset = 0; offset <= repeat_count; ++offset) {
	  cv_buf_check[index+offset] = icv_global;
	}
      }
    }

    for (int i = 0; i < image_size; ++i) {
      //cout << "i,cv_buf[i],cv_buf_check[i]: " << i << " " << cv_buf[i] << " " << cv_buf_check[i] << endl;
      //getchar();
      if (cv_buf[i] >= 0) {
	assert(cv_buf[i] == cv_buf_check[i]);
      }
    }

    delete[] cv_buf_check;

  }

  delete[] cv_buf;

  return 1;
}


//If cim was read and mesh imaging was requested an additional buffer needs to be computed.
//If the cim was just built, the buffer was already computed.
//Point to this buffer with mesh_buf, so that d2_buf can be used by a single LesImageMapper
//to compute additional cims. (For a given canvas view and surface field and plane may
//be compute using the same LesImageMapper instance)
void LesImageMapper::computeCvPixelD2(const CvImageMap& cim, const CtiCanvas * canvas, const float xp[3], const float np[3]){

  assert(mesh_buf==NULL);
  if (d2_buf){ //reassign d2_buf data to mesh_buf pointer
    mesh_buf = d2_buf;
    d2_buf = NULL;
  }
  else{  //compute d2_buf data but store in mesh_buf

    const int image_size = canvas->getNi()*canvas->getNj();
    mesh_buf = new float[image_size]; // holds the dist2 of the closest cv, if there is one

    float x0[3];  //simulation coordinates
    float e0[3];  //NOTE, not unit vectors, have length ni/width
    float e1[3];
    float e2[3];
    canvas->getX0(x0);
    canvas->getE0(e0);
    canvas->getE1(e1);
    canvas->getE2(e2);

    //Data Plane vectors, simulation coordinates
    float mag_np = MAG(np); assert(mag_np>0.0);
    float unit_np[3];
    FOR_I3 unit_np[i] = np[i]/mag_np;

    //plane center and normal in canvas coordinates
    //use this later for depth testing CV-to plane pixel location distance
    float xp_c[3] = { ( (xp[0]-x0[0])*e0[0] + (xp[1]-x0[1])*e0[1] + (xp[2]-x0[2])*e0[2] ),
                      ( (xp[0]-x0[0])*e1[0] + (xp[1]-x0[1])*e1[1] + (xp[2]-x0[2])*e1[2] ),
                      ( (xp[0]-x0[0])*e2[0] + (xp[1]-x0[1])*e2[1] + (xp[2]-x0[2])*e2[2] )};
    float unit_np_c[3] = { ( (unit_np[0])*e0[0] + (unit_np[1])*e0[1] + (unit_np[2])*e0[2] ),
                           ( (unit_np[0])*e1[0] + (unit_np[1])*e1[1] + (unit_np[2])*e1[2] ),
                           ( (unit_np[0])*e2[0] + (unit_np[1])*e2[1] + (unit_np[2])*e2[2] )};
    float np_c_mag = MAG(unit_np_c);
    FOR_I3 unit_np_c[i] /= np_c_mag;

    double (*x_vv)[3] = new double[cim.ncv_unique][3];
    sparseReadCvR2(x_vv,cim,"x_vv");

    for (int icv = 0; icv < cim.ncv_unique; ++icv) {
      // const int8 icv_global = cim.ucv_v[icv];
      //compute x_vv in canvas coordinates
      float x_vv_c[3] = { (float) ( (x_vv[icv][0]-x0[0])*e0[0] + (x_vv[icv][1]-x0[1])*e0[1] + (x_vv[icv][2]-x0[2])*e0[2] ),
                          (float) ( (x_vv[icv][0]-x0[0])*e1[0] + (x_vv[icv][1]-x0[1])*e1[1] + (x_vv[icv][2]-x0[2])*e1[2] ),
                          (float) ( (x_vv[icv][0]-x0[0])*e2[0] + (x_vv[icv][1]-x0[1])*e2[1] + (x_vv[icv][2]-x0[2])*e2[2] )};

      for (int poc = cim.poucv_i[icv]; poc != cim.poucv_i[icv+1]; ++poc) {
        unsigned int index = cim.poucv_v[poc] & ((1<<(32-CIM_REPEAT_BITS))-1);
        unsigned int repeat_count = cim.poucv_v[poc]>>(32-CIM_REPEAT_BITS);
        assert( (index | (repeat_count<<(32-CIM_REPEAT_BITS))) == cim.poucv_v[poc]);
        for (unsigned int offset = 0; offset <= repeat_count; ++offset) {
          const int image_ij = index+offset;

          const int jp = (int) image_ij/canvas->getNi();
          const int ip = image_ij%canvas->getNi();
          double kp = xp_c[2] - (unit_np_c[0]*(ip-xp_c[0]) +
                                 unit_np_c[1]*(jp-xp_c[1]) )
                                /unit_np_c[2];

          mesh_buf[image_ij] = (ip-x_vv_c[0])*(ip-x_vv_c[0])+(jp-x_vv_c[1])*(jp-x_vv_c[1])+(kp-x_vv_c[2])*(kp-x_vv_c[2]);
        }
      }
    }
    delete[] x_vv;
  }
}

void LesImageMapper::sparseReadCvR1(double * var,const CvImageMap& cim,const string &img_var) {
  int8 data_offset = getRfpCvD1Offset(img_var);
  sparseReadCvR1(var,cim,img_var,data_offset);
}


void LesImageMapper::sparseReadCvR1(double * var,const CvImageMap& cim,const string &img_var, const int8 &data_offset) {

  const int max_buf_size = 100000;
  cout << "sparseReadCvR1, max_buf_size=" << max_buf_size << "...";
  cout.flush();

  double * rbuf = new double[max_buf_size];

  int icv_end = 0;
  while (icv_end < cim.ncv_unique) {
    const int icv_begin = icv_end++;
    while ((icv_end < cim.ncv_unique)&&(cim.ucv_v[icv_end]-cim.ucv_v[icv_begin] < max_buf_size))
      ++icv_end;
    // seek to icv_begin...
    my_fseek(rfp,data_offset+cim.ucv_v[icv_begin]*sizeof(double));
    fread(rbuf,sizeof(double),cim.ucv_v[icv_end-1]-cim.ucv_v[icv_begin]+1,rfp);
    for (int icv = icv_begin; icv != icv_end; ++icv){
      var[icv] = rbuf[cim.ucv_v[icv]-cim.ucv_v[icv_begin]];
    }
  }

  delete[] rbuf;
  cout << "done" << endl;
  //MiscUtils::dumpRange(var,cim.ncv_unique,"SparseVar");

}

void LesImageMapper::sparseReadCvR2(double (*varR2)[3],const CvImageMap& cim,const string &img_var) {
  int8 data_offset = getRfpCvD2Offset(img_var);
  sparseReadCvR2(varR2,cim,img_var,data_offset);
}

void LesImageMapper::sparseReadCvR2(double (*varR2)[3],const CvImageMap& cim,const string &img_var,const int8 &data_offset) {

  const int max_buf_size = 100000;
  cout << "sparseReadCvR2, max_buf_size=" << 3*max_buf_size << "...";
  cout.flush();

  double * rbuf = new double[3*max_buf_size];

  int icv_end = 0;
  while (icv_end < cim.ncv_unique) {
    const int icv_begin = icv_end++;
    while ((icv_end < cim.ncv_unique)&&(cim.ucv_v[icv_end]-cim.ucv_v[icv_begin] < max_buf_size))
      ++icv_end;
    // seek to icv_begin...
    my_fseek(rfp,data_offset+cim.ucv_v[icv_begin]*sizeof(double)*3);
    fread(rbuf,sizeof(double),3*(cim.ucv_v[icv_end-1]-cim.ucv_v[icv_begin]+1),rfp);
    for (int icv = icv_begin; icv != icv_end; ++icv){
      FOR_I3 varR2[icv][i] = rbuf[3*(cim.ucv_v[icv]-cim.ucv_v[icv_begin])+i];
    }
  }

  delete[] rbuf;
  cout << "done" << endl;
  //MiscUtils::dumpRange(varR2,cim.ncv_unique,"SparseVar");

}

void LesImageMapper::sparseReadCvR2Mag(double *varR2Mag,const CvImageMap& cim,const string &img_var,const int8 &data_offset) {

  const int max_buf_size = 100000;
  cout << "sparseReadCvR2Mag, max_buf_size=" << 3*max_buf_size << "...";
  cout.flush();

  double * rbuf = new double[3*max_buf_size];

  int icv_end = 0;
  while (icv_end < cim.ncv_unique) {
    const int icv_begin = icv_end++;
    while ((icv_end < cim.ncv_unique)&&(cim.ucv_v[icv_end]-cim.ucv_v[icv_begin] < max_buf_size))
      ++icv_end;
    // seek to icv_begin...
    my_fseek(rfp,data_offset+cim.ucv_v[icv_begin]*sizeof(double)*3);
    fread(rbuf,sizeof(double),3*(cim.ucv_v[icv_end-1]-cim.ucv_v[icv_begin]+1),rfp);
    for (int icv = icv_begin; icv != icv_end; ++icv){
      double r2_x = rbuf[3*(cim.ucv_v[icv]-cim.ucv_v[icv_begin])+0];
      double r2_y = rbuf[3*(cim.ucv_v[icv]-cim.ucv_v[icv_begin])+1];
      double r2_z = rbuf[3*(cim.ucv_v[icv]-cim.ucv_v[icv_begin])+2];
      varR2Mag[icv] = sqrt(r2_x*r2_x+r2_y*r2_y+r2_z*r2_z);
    }
  }

  delete[] rbuf;
  cout << "done" << endl;
  //MiscUtils::dumpRange(var,cim.ncv_unique,"SparseVar");

}


void LesImageMapper::my_fseek(FILE * fp, int8 offset) {

  // not very efficient, but works on BG
  int8 cur_offset = offset;
  int8 max_offset = 1000000000ll;

  fseek(fp,0,SEEK_SET);
  while (cur_offset > max_offset) {
    fseek(fp,max_offset,SEEK_CUR);
    cur_offset -= max_offset;
  }
  fseek(fp,cur_offset,SEEK_CUR);
}

int8 LesImageMapper::atoi8(const string& var) {
  std::istringstream iss(var);
  int8 value;
  if ((iss >> value >> std::dec).fail()) {
    //cerr << "\nError: cannot evaluate string to int8: \"" << var << "\"" << endl;
    throw(-1);
  }
  return value;
}

void LesImageMapper::checkContiguous(FILE * fp,const int n) {
  cout << "# Reading and checking " << n << " contiguous ints";
  const int read_size = 10000000;
  int *buf = new int[min(n,read_size)];
  int8 index0 = 0;
  while (index0 < n) {
    cout << ".";
    cout.flush();
    const int this_read_size = min(read_size,int(n-index0));
    fread(buf,sizeof(int),this_read_size,fp);
    for (int i = 0; i < this_read_size; ++i) {
      if (buf[i] != index0+i) {
	cerr << "\nError: value at position " << index0+i << " does not match: " << buf[i];
	throw(-1);
      }
    }
    index0 += this_read_size;
  }
  delete[] buf;
  cout << "OK" << endl;
}

int8 LesImageMapper::advanceRfp(const int ugp_io_tag,const string& name) {

  using namespace ByteSwap;

  //cout << "advanceRfp: looking for header: " << name << endl;

  Header header;
  while (!rfp_eof) {

    my_fseek(rfp,rfp_offset);

    fread(&header,sizeof(Header),1,rfp);
    if (rfp_byte_swap)
      byteSwapHeader(&header,1);
    assert(header.skip > 0);

    //cout << "advanceRfp: got header: " << header.name << endl;

    //use lowercase for all map keys
    string headerNameCI = header.name;
    transform(headerNameCI.begin(),headerNameCI.end(), headerNameCI.begin(), ::tolower);

    if (header.id == UGP_IO_EOF) {
      rfp_eof = true;
    }
    else if (header.id == UGP_IO_NO_FA_CV_COUNTS){
      if (io_version==3){
        got_ncv_global = true;
        ncv_global = getInt8FromLswMswPair(header.idata+4);
        cout << "UGP_IO_NO_FA_CV_COUNTS " << header.name << " ncv_global " << ncv_global << endl;
      }
    }
    else if (header.id == UGP_IO_CV_D1) {
      cout << "CV_D1 " << header.name << endl;
      // set or check ncv_global...
      if (!got_ncv_global&&io_version>3) {
	got_ncv_global = true;
        ncv_global = getInt8FromLswMswPair(header.idata);
      }
      // map this location...
      cvD1OffsetMap[headerNameCI] = rfp_offset + sizeof(Header);
    }
    else if (header.id == UGP_IO_CV_D2) {
      cout << "CV_D2 " << header.name << endl;
      // set or check ncv_global...
      if (!got_ncv_global&&io_version>3) {
	got_ncv_global = true;
        ncv_global = getInt8FromLswMswPair(header.idata);
      }
      // map this cv-d2 location...
      cvD2OffsetMap[headerNameCI] = rfp_offset + sizeof(Header);
    }
    else if (header.id == UGP_IO_I0) {
      cout << "I0: " << header.name << "=" << header.idata[0] << endl;
      i0ValueMap[headerNameCI] = header.idata[0];
    }
    else if (header.id == UGP_IO_D0) {
      cout << "R0: " << header.name << "=" << header.rdata[0] << endl;
      d0ValueMap[headerNameCI] = header.rdata[0];
    }
    else if (header.id == UGP_IO_HASHIDS) {
      RestartHashUtilities::clearHashes();
      RestartHashUtilities::readHashFromRestartHeader_s(rfp,header);
      encodedHash = RestartHashUtilities::myHash.getAsciiHash();
      bSetEncodedHash=true;
      cout << header.name << " mles " << encodedHash << endl;
    }
    else if (header.id == UGP_IO_VV_BBOX) {
      cout << "VV_BBOX: " << header.name << " " << header.idata[0] << endl;
      assert(header.idata[0] >= 1); // should have a bounding box
      vvbboxPair = make_pair(header.idata[0],rfp_offset+sizeof(Header));
    }
    else if (header.id == UGP_IO_CVOBF_INT8) {
      nbf_global = getInt8FromLswMswPair(header.idata);
      got_nbf_global = true;
      cvobf_offset = rfp_offset + sizeof(Header);
      cout << "CVOBF: " << header.name << " nbf_global " << nbf_global << endl;
    }
    else if (header.id == UGP_IO_BF_ZONE_HEADER){
        int8 rfp_offset_first = rfp_offset;
      do{
        const int8 nbf_zn_global   = ByteSwap::getInt8FromLswMswPair(header.idata+2);
        int8 offset_zn_local=0;
        if (nbfZnOffsetVec.size()>0){
          offset_zn_local += nbfZnOffsetVec.back().first;
          offset_zn_local += nbfZnOffsetVec.back().second;
        }
        nbfZnOffsetVec.push_back(make_pair(nbf_zn_global,offset_zn_local));
        cout << "BF_ZONE: " << header.name << " nbf_zn_global " << nbf_zn_global << endl;
        //advance to next record and read if it's another bf zone
        rfp_offset += header.skip;
        my_fseek(rfp,rfp_offset);
        fread(&header,sizeof(Header),1,rfp);
        if (rfp_byte_swap)
          byteSwapHeader(&header,1);
        assert(header.skip > 0);
      } while (header.id == UGP_IO_BF_ZONE_HEADER);

      if (ugp_io_tag == UGP_IO_BF_ZONE_HEADER) //we were looking for UGP_IO_BF_ZONE_HEADER, return offset
        return rfp_offset_first + sizeof(Header);
      else  //we have another header and need to return control to the main loop with rfp_offset just after final UGP_IO_BF_ZONE_HEADER header
        rfp_offset -= header.skip; //the offset is correct for reading the header again, just negate the skip applied outside the loop.
    }

    rfp_offset += header.skip;

    // ===========================================
    // look for the requested header and return if
    // matched...
    // ===========================================

    if ((header.id == ugp_io_tag)&&(equals_ci(headerNameCI,name))) {
      return rfp_offset - header.skip + sizeof(Header);
    }

  }

  // must not have found it...
  return -1;

}

bool LesImageMapper::setRfpHashId(){
  if (bSetEncodedHash){
    return true;
  }

  // it was not stored in the global var
  // advance rfp if possible
  if (rfp_eof)
    return false;

  int8 offset = advanceRfp(UGP_IO_HASHIDS,"UGP_IO_HASHIDS");
  if (offset < 0){
    //header was not found
    return false;
  }
  else{
    //header was found and set to encodedHash
    return true;
  }
}

bool LesImageMapper::getRfpI0Value(const string& name,int &i0Value){

  // first check the map...
  string nameLC = name;
  transform(nameLC.begin(),nameLC.end(), nameLC.begin(), ::tolower);

  map<const string,int>::iterator iter = i0ValueMap.find(nameLC);
  if (iter != i0ValueMap.end()){
    i0Value = iter->second;
    return true;
  }

  // it was not in the map. If the map is complete, then offset cannot be completed...

  if (rfp_eof)
    return false;

  // otherwise, advance the rfp and look for the

  int8 offset = advanceRfp(UGP_IO_I0,nameLC);
  if (offset < 0){
    //header was not found
    return false;
  }
  else{
    //header was found and map set, grab value
    i0Value = i0ValueMap[nameLC];
    return true;
  }

}

bool LesImageMapper::getRfpD0Value(const string& name,double &d0Value){

  // first check the map...
  string nameLC = name;
  transform(nameLC.begin(),nameLC.end(), nameLC.begin(), ::tolower);

  map<const string,double>::iterator iter = d0ValueMap.find(nameLC);
  if (iter != d0ValueMap.end()){
    d0Value = iter->second;
    return true;
  }

  // it was not in the map. If the map is complete, then offset cannot be completed...

  if (rfp_eof)
    return false;

  // otherwise, advance the rfp and look for the

  int8 offset = advanceRfp(UGP_IO_D0,nameLC);
  if (offset < 0){
    //header was not found
    return false;
  }
  else{
    //header was found and map set, grab value
    d0Value = d0ValueMap[nameLC];
    return true;
  }


}

int8 LesImageMapper::getRfpCvD1Offset(const string& name) {

  assert(rfp != NULL);

  // first check the map...
  string nameLC = name;
  transform(nameLC.begin(),nameLC.end(), nameLC.begin(), ::tolower);

  map<const string,int8>::iterator iter = cvD1OffsetMap.find(nameLC);
  if (iter != cvD1OffsetMap.end())
    return iter->second;

  // it was not in the map. If the map is complete, then offset cannot be completed...

  if (rfp_eof)
    return -1;

  // otherwise, advance the rfp and look for the

  return advanceRfp(UGP_IO_CV_D1,nameLC);

}

int8 LesImageMapper::getRfpCvD2Offset(const string& name) {

  assert(rfp != NULL);

  // first check the map...
  string nameLC = name;
  transform(nameLC.begin(),nameLC.end(), nameLC.begin(), ::tolower);


  map<const string,int8>::iterator iter = cvD2OffsetMap.find(nameLC);
  if (iter != cvD2OffsetMap.end())
    return iter->second;

  // it was not in the map. If the map is complete, then offset cannot be completed...

  if (rfp_eof)
    return -1;

  // otherwise, advance the rfp and look for the

  return advanceRfp(UGP_IO_CV_D2,nameLC);

}

int8 LesImageMapper::getRfpVVBBoxOffset(const string& name, int &nbb) {

  assert(rfp != NULL);

  // first check the map...
  if (vvbboxPair.second >= 0){
    assert(vvbboxPair.first > 0);
    nbb = vvbboxPair.first;
    return vvbboxPair.second;
  }

  // it was not in the map. If the map is complete, then offset cannot be completed...

  if (rfp_eof)
    return -1;

  // otherwise, advance the rfp and look for the
  string nameLC = name;
  transform(nameLC.begin(),nameLC.end(), nameLC.begin(), ::tolower);

  int8 theOffset = advanceRfp(UGP_IO_VV_BBOX,nameLC);
  if (theOffset>=0){
    nbb = vvbboxPair.first;
  }
  return theOffset;
}

void LesImageMapper::getNcvGlobalIoV3() {

  if (!got_ncv_global){
    assert(rfp != NULL);
    if (rfp_eof)
      CERR("Unable to retrieve ncv_global for io_version 3 restart file");

    cout << "Getting ncv_global for io_version 3 restart file" << endl;

    advanceRfp(UGP_IO_NO_FA_CV_COUNTS,"no_fa_cv_counts");
  }
  if (!got_ncv_global)
    CERR("Unable to retrieve ncv_global for io_version 3 restart file");

}

int8 LesImageMapper::getRfpCvobfOffsetAndNbfGlobal() {

  if (!got_nbf_global){
    assert(rfp != NULL);
    if (rfp_eof)
      return -1;

    return advanceRfp(UGP_IO_CVOBF_INT8,"cvobf_global");
  }
  else{
    return  cvobf_offset;
  }

}

bool LesImageMapper::getRfpNbfZnOffsetVec() {
  if (nbfZnOffsetVec.size()==0){
    assert(rfp != NULL);
    if (rfp_eof)
      return false;

    return (advanceRfp(UGP_IO_BF_ZONE_HEADER,"bfZone")>-1);
  }
  else
    return true;
}

//assumes a vvbbox record as been found and set
bool LesImageMapper::checkVVBBoxIntersection(int &ncv_bbox,const int &ibb,const float xp[3], const float np[3]) {
  //read bounding box data once
  if (ncv_bbox_buf == NULL) {
    int nbb = vvbboxPair.first;
    int8 offset = vvbboxPair.second;
    if (offset < 0) //no vvbbox record, shouldn't get here
      assert(0);
    ncv_bbox_buf = new int[nbb];
    my_fseek(rfp,offset);
    fread(ncv_bbox_buf,sizeof(int),nbb,rfp);
    assert(vv_bbox_buf == NULL);
    vv_bbox_buf = new double[nbb][6];
    fread(vv_bbox_buf,sizeof(double),nbb*6,rfp);
  }
  ncv_bbox = ncv_bbox_buf[ibb];

  // Convert AABB to center-extents representation
  double bbox_c[3], bbox_e[3];
  FOR_I3 bbox_c[i] = 0.5*(vv_bbox_buf[ibb][i] + vv_bbox_buf[ibb][i+3]); //center of bbox
  FOR_I3 bbox_e[i] = vv_bbox_buf[ibb][i+3] - bbox_c[i];        //positive extents

  // Compute the projection interval radius of bbox onto L(t) = bbox_c + t * np
  double r = bbox_e[0]*fabs(np[0]) + bbox_e[1]*fabs(np[1]) + bbox_e[2]*fabs(np[2]);

  // Compute distance of box center from plane
  double xpToBboxC[3] = DIFF(bbox_c,xp);
  double s = DOT_PRODUCT(xpToBboxC,np);

  // Intersection occurs when distance s falls within [-r,+r] interval
  return fabs(s) <= r;

}

//assumes a vvbbox record as been found and set
bool LesImageMapper::checkVVBBoxBBoxIntersection(int &ncv_bbox,const int &ibb,const double bbox[6]) {
  // DP: why are we re-reading this record above?? - read once and store data -- same for above routine...
  if (ncv_bbox_buf == NULL) {
    int nbb = vvbboxPair.first;
    int8 offset = vvbboxPair.second;
    if (offset < 0) //no vvbbox record, shouldn't get here
      assert(0);
    ncv_bbox_buf = new int[nbb];
    my_fseek(rfp,offset);
    fread(ncv_bbox_buf,sizeof(int),nbb,rfp);
    assert(vv_bbox_buf == NULL);
    vv_bbox_buf = new double[nbb][6];
    fread(vv_bbox_buf,sizeof(double),nbb*6,rfp);
  }
  ncv_bbox = ncv_bbox_buf[ibb];
  // note that vv_bbox is xmin,ymin,zmin,xmax,ymax,zmax...
  if ((bbox[0] > vv_bbox_buf[ibb][3])||(bbox[1] < vv_bbox_buf[ibb][0])||
      (bbox[2] > vv_bbox_buf[ibb][4])||(bbox[3] < vv_bbox_buf[ibb][1])||
      (bbox[4] > vv_bbox_buf[ibb][5])||(bbox[5] < vv_bbox_buf[ibb][2]))
    return false;
  return true;
}

// Functions below added to handle mles/sles paradigm --

void LesImageMapper::sparseReadCvR1_2(double * var,const CvImageMap& cim,const string &img_var, const int8 &data_offset) {

  const int max_buf_size = 100000;
  cout << "sparseReadCvR1_2, max_buf_size=" << max_buf_size << "...";
  cout.flush();

  double * rbuf = new double[max_buf_size];

  int icv_end = 0;
  while (icv_end < cim.ncv_unique) {
    const int icv_begin = icv_end++;
    while ((icv_end < cim.ncv_unique)&&(cim.ucv_v[icv_end]-cim.ucv_v[icv_begin] < max_buf_size))
      ++icv_end;
    // seek to icv_begin...
    my_fseek(sfp,data_offset+cim.ucv_v[icv_begin]*sizeof(double));
    fread(rbuf,sizeof(double),cim.ucv_v[icv_end-1]-cim.ucv_v[icv_begin]+1,sfp);
    for (int icv = icv_begin; icv != icv_end; ++icv){
      var[icv] = rbuf[cim.ucv_v[icv]-cim.ucv_v[icv_begin]];
    }
  }

  delete[] rbuf;
  cout << "done" << endl;
  //MiscUtils::dumpRange(var,cim.ncv_unique,"SparseVar");

}

void LesImageMapper::sparseReadCvR2_2(double (*varR2)[3],const CvImageMap& cim,const string &img_var,const int8 &data_offset) {

  const int max_buf_size = 100000;
  cout << "sparseReadCvR2, max_buf_size=" << 3*max_buf_size << "...";
  cout.flush();

  double * rbuf = new double[3*max_buf_size];

  int icv_end = 0;
  while (icv_end < cim.ncv_unique) {
    const int icv_begin = icv_end++;
    while ((icv_end < cim.ncv_unique)&&(cim.ucv_v[icv_end]-cim.ucv_v[icv_begin] < max_buf_size))
      ++icv_end;
    // seek to icv_begin...
    my_fseek(sfp,data_offset+cim.ucv_v[icv_begin]*sizeof(double)*3);
    fread(rbuf,sizeof(double),3*(cim.ucv_v[icv_end-1]-cim.ucv_v[icv_begin]+1),sfp);
    for (int icv = icv_begin; icv != icv_end; ++icv){
      FOR_I3 varR2[icv][i] = rbuf[3*(cim.ucv_v[icv]-cim.ucv_v[icv_begin])+i];
    }
  }

  delete[] rbuf;
  cout << "done" << endl;
  //MiscUtils::dumpRange(varR2,cim.ncv_unique,"SparseVar");

}

void LesImageMapper::sparseReadCvR2Mag_2(double *varR2Mag,const CvImageMap& cim,const string &img_var,const int8 &data_offset) {

  const int max_buf_size = 100000;
  cout << "sparseReadCvR2Mag_2, max_buf_size=" << 3*max_buf_size << "...";
  cout.flush();

  double * rbuf = new double[3*max_buf_size];

  int icv_end = 0;
  while (icv_end < cim.ncv_unique) {
    const int icv_begin = icv_end++;
    while ((icv_end < cim.ncv_unique)&&(cim.ucv_v[icv_end]-cim.ucv_v[icv_begin] < max_buf_size))
      ++icv_end;
    // seek to icv_begin...
    my_fseek(sfp,data_offset+cim.ucv_v[icv_begin]*sizeof(double)*3);
    fread(rbuf,sizeof(double),3*(cim.ucv_v[icv_end-1]-cim.ucv_v[icv_begin]+1),sfp);
    for (int icv = icv_begin; icv != icv_end; ++icv){
      double r2_x = rbuf[3*(cim.ucv_v[icv]-cim.ucv_v[icv_begin])+0];
      double r2_y = rbuf[3*(cim.ucv_v[icv]-cim.ucv_v[icv_begin])+1];
      double r2_z = rbuf[3*(cim.ucv_v[icv]-cim.ucv_v[icv_begin])+2];
      varR2Mag[icv] = sqrt(r2_x*r2_x+r2_y*r2_y+r2_z*r2_z);
    }
  }

  delete[] rbuf;
  cout << "done" << endl;
  //MiscUtils::dumpRange(varR2Mag,cim.ncv_unique,"SparseVar");

}

int8 LesImageMapper::advanceSfp(const int ugp_io_tag,const string& name) {

  using namespace ByteSwap;

  //cout << "advanceSfp: looking for header: " << name << endl;

  Header header;
  while (!sfp_eof) {

    my_fseek(sfp,sfp_offset);

    fread(&header,sizeof(Header),1,sfp);
    if (sfp_byte_swap)
      byteSwapHeader(&header,1);
    assert(header.skip > 0);

    //cout << "advanceSfp: got header: " << header.name << endl;

    //use lowercase for all map keys
    string headerNameCI = header.name;
    transform(headerNameCI.begin(),headerNameCI.end(), headerNameCI.begin(), ::tolower);

    if (header.id == UGP_IO_EOF) {
      sfp_eof = true;
    }
    else if (header.id == UGP_IO_CV_D1) {
      cout << "CV_D1 " << header.name << endl;
      // set or check ncv_global...
      if (!got_ncv_global) {
	got_ncv_global = true;
        ncv_global = getInt8FromLswMswPair(header.idata);
      }
      // map this location...
      cvD1OffsetMap_2[headerNameCI] = sfp_offset + sizeof(Header);
    }
    else if (header.id == UGP_IO_CV_D2) {
      cout << "CV_D2 " << header.name << endl;
      // set or check ncv_global...
      if (!got_ncv_global) {
	got_ncv_global = true;
        ncv_global = getInt8FromLswMswPair(header.idata);
      }
      // map this cv-d2 location...
      cvD2OffsetMap_2[headerNameCI] = sfp_offset + sizeof(Header);
    }
    else if (header.id == UGP_IO_I0) {
      cout << "I0: " << header.name << "=" << header.idata[0] << endl;
      i0ValueMap_2[headerNameCI] = header.idata[0];
    }
    else if (header.id == UGP_IO_D0) {
      cout << "R0: " << header.name << "=" << header.rdata[0] << endl;
      d0ValueMap_2[headerNameCI] = header.rdata[0];
    }
    else if (header.id == UGP_IO_HASHIDS) {
      RestartHashUtilities::clearHashes();
      RestartHashUtilities::readHashFromRestartHeader_s(sfp,header);
      encodedHash_2 = RestartHashUtilities::myHash.getAsciiHash();
      string encodedParentHash = RestartHashUtilities::mlesHash.getAsciiHash();
      bSetEncodedHash_2 = true;
      cout << header.name << " sles " << encodedHash_2 << endl;
      if (bSetEncodedHash)
        assert(encodedParentHash == encodedHash); // sles matches with mles??? TODO
    }

    sfp_offset += header.skip;

    // ===========================================
    // look for the requested header and return if
    // matched...
    // ===========================================

    if ((header.id == ugp_io_tag)&&(equals_ci(headerNameCI,name))) {
      return sfp_offset - header.skip + sizeof(Header);
    }

  }

  // must not have found it...
  return -1;

}

bool LesImageMapper::getSfpI0Value(const string& name,int &i0Value){

  // first check the map...
  string nameLC = name;
  transform(nameLC.begin(),nameLC.end(), nameLC.begin(), ::tolower);

  map<const string,int>::iterator iter = i0ValueMap_2.find(nameLC);
  if (iter != i0ValueMap_2.end()){
    i0Value = iter->second;
    return true;
  }

  // it was not in the map. If the map is complete, then offset cannot be completed...

  if (sfp_eof)
    return false;

  // otherwise, advance the rfp and look for the

  int8 offset = advanceSfp(UGP_IO_I0,nameLC);
  if (offset < 0){
    //header was not found
    return false;
  }
  else{
    //header was found and map set, grab value
    i0Value = i0ValueMap_2[nameLC];
    return true;
  }

}

bool LesImageMapper::getSfpD0Value(const string& name,double &d0Value){

  // first check the map...
  string nameLC = name;
  transform(nameLC.begin(),nameLC.end(), nameLC.begin(), ::tolower);

  map<const string,double>::iterator iter = d0ValueMap_2.find(nameLC);
  if (iter != d0ValueMap_2.end()){
    d0Value = iter->second;
    return true;
  }

  // it was not in the map. If the map is complete, then offset cannot be completed...

  if (sfp_eof)
    return false;

  // otherwise, advance the sfp and look for the

  int8 offset = advanceSfp(UGP_IO_D0,nameLC);
  if (offset < 0){
    //header was not found
    return false;
  }
  else{
    //header was found and map set, grab value
    d0Value = d0ValueMap_2[nameLC];
    return true;
  }


}

int8 LesImageMapper::getSfpCvD1Offset(const string& name) {

  assert(sfp != NULL);

  // first check the map...
  string nameLC = name;
  transform(nameLC.begin(),nameLC.end(), nameLC.begin(), ::tolower);

  map<const string,int8>::iterator iter = cvD1OffsetMap_2.find(nameLC);
  if (iter != cvD1OffsetMap_2.end())
    return iter->second;

  // it was not in the map. If the map is complete, then offset cannot be completed...

  if (sfp_eof)
    return -1;

  // otherwise, advance the sfp and look for the

  return advanceSfp(UGP_IO_CV_D1,nameLC);

}

int8 LesImageMapper::getSfpCvD2Offset(const string& name) {

  assert(sfp != NULL);

  // first check the map...
  string nameLC = name;
  transform(nameLC.begin(),nameLC.end(), nameLC.begin(), ::tolower);

  map<const string,int8>::iterator iter = cvD2OffsetMap_2.find(nameLC);
  if (iter != cvD2OffsetMap_2.end())
    return iter->second;

  // it was not in the map. If the map is complete, then offset cannot be completed...

  if (sfp_eof)
    return -1;

  // otherwise, advance the sfp and look for the

  return advanceSfp(UGP_IO_CV_D2,nameLC);

}

bool LesImageMapper::setSfpHashId(){
  if (bSetEncodedHash_2){
    return true;
  }

  // it was not stored in the global var
  // advance sfp if possible
  if (sfp_eof)
    return false;

  int8 offset = advanceSfp(UGP_IO_HASHIDS,"UGP_IO_HASHIDS");
  if (offset < 0){
    //header was not found
    return false;
  }
  else{
    //header was found and set to encodedHash
    return true;
  }
}
