#ifndef _CV_IMAGE_MAP_
#define _CV_IMAGE_MAP_

// the I,J location of a pixel is converted to a 1D pixel index
// using this function. This means the compression (described next)
// is along lines in the x-direction.
//
//#define PINDEX(I,J) ((J)*_parser.IMAGE_NX+(I))
 // grab the step number for image naming
      // grab the step number for image naming

// the number of repeat bits help to compress the image data
// by repeating contiguous pixels up to 2^CIM_REPEAT_BITS times,
// however repeat bits take away from the image indexing, so
// limit the maximum size of the image allowed. Note that if
// a particular value is asociated with more than 2^CIM_REPEAT_BITS+1
// pixels in a row, this is fine. The data will have multiple
// poucv_v entries, each (except perhaps the last) with 2^CIM_REPEAT_BITS+1
// pixels.
//
// in choosing the size here, 10000x10000 is a BIG image, and requires
// 27 bits to store the pixel indexing (2^27 = 134217728 > 10000*10000),
// and leaves 5 bits for pixel repeating (2^5 = 32).
//
// so, we recommend CIM_REPEAT_BITS=5

#define CIM_VERSION 1
#define CIM_REPEAT_BITS 5

class CvImageMap {
public:
  int ncv_unique;
  int8 * ucv_v; // unique-cv
  int * poucv_i; // pixel-of-unique-cv_i
  unsigned int * poucv_v; // pixel-of-unique-cv_v, + count of pixels after the first using top CIM_REPEAT_BITS bits

  CvImageMap() {
    ncv_unique = 0;
    ucv_v = NULL;
    poucv_i = NULL;
    poucv_v = NULL;
  }
  ~CvImageMap() {
    DELETE(ucv_v);
    DELETE(poucv_i);
    DELETE(poucv_v);
  }

  unsigned int getIndex(const int &poc) const{
    return poucv_v[poc] & ((1<<(32-CIM_REPEAT_BITS))-1);
  }

  unsigned int getRepeatCount(const int &poc) const{
    return poucv_v[poc]>>(32-CIM_REPEAT_BITS);
  }


  int write(const char * filename) {

    // ----------------------------------------
    // write a binary imap file.
    // returns -1 if a problem, 0 otherwise
    // ----------------------------------------

    FILE * fp = fopen(filename,"wb");
    if (fp == NULL)
      return -1;

    int ibuf[3];
    ibuf[0] = CIM_VERSION;
    ibuf[1] = CIM_REPEAT_BITS;
    ibuf[2] = ncv_unique;
    fwrite(ibuf,sizeof(int),3,fp);


    // recall...
    /*
      cim.ucv_v = new int8[cim.ncv_unique];
      cim.poucv_i = new int[cim.ncv_unique+1];
      cim.poucv_v = new unsigned int[poucv_s];
    */

    fwrite(ucv_v,sizeof(int8),ncv_unique,fp);

    // convert the poucv_i offset array into a count per ucv -
    // this will allow

    const int poucv_s = poucv_i[ncv_unique];

    for (int icv = ncv_unique; icv > 0; --icv)
      poucv_i[icv] -= poucv_i[icv-1];

    fwrite(poucv_i+1,sizeof(int),ncv_unique,fp);

    // and return poucv_i...

    assert(poucv_i[0] == 0);
    for (int icv = 0; icv < ncv_unique; ++icv)
      poucv_i[icv+1] += poucv_i[icv];

    // check...

    assert(poucv_s == poucv_i[ncv_unique]);

    // and write poucv_v...

    fwrite(poucv_v,sizeof(unsigned int),poucv_s,fp);

    fclose(fp);

    cout << " > Wrote CvImageMap with filename: " << filename << endl;
    return 0;
  }

  int read(const char * filename) {

    // ----------------------------------------
    // reads the binary cim file.
    // returns -1 if a problem, 0 otherwise
    // ----------------------------------------
    FILE * fp = NULL;
    const int file_err = MiscUtils::openFile(&fp,string(filename),"rb");
    if (file_err != 0) return -1;
    // FILE * fp = fopen(filename,"rb");
    // if (fp == NULL)
    //   return -1;

    int ibuf[3];
    fread(ibuf,sizeof(int),3,fp);
    if ((ibuf[0] != CIM_VERSION)||(ibuf[1] != CIM_REPEAT_BITS)){
      cout << "Warning: problem with binary cim format" << endl;
      return -1;
    }
    ncv_unique = ibuf[2];

    assert(ucv_v == NULL);
    ucv_v = new int8[ncv_unique];
    fread(ucv_v,sizeof(int8),ncv_unique,fp);

    assert(poucv_i == NULL);
    poucv_i = new int[ncv_unique+1];
    fread(poucv_i+1,sizeof(int),ncv_unique,fp);

    // and make poucv_i CSR format...

    poucv_i[0] = 0;
    for (int icv = 0; icv < ncv_unique; ++icv)
      poucv_i[icv+1] += poucv_i[icv];

    // get size...

    const int poucv_s = poucv_i[ncv_unique];
    poucv_v = new unsigned int[poucv_s];
    fread(poucv_v,sizeof(unsigned int),poucv_s,fp);

    fclose(fp);
    cout << " > Found and read CvImageMap with filename: " << filename << endl;
    return 0;

  }



};


#endif
