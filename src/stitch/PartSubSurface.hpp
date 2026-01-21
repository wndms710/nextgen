#ifndef _PART_SUB_SURFACE_HPP_
#define _PART_SUB_SURFACE_HPP_

// similar to the SubSurface class, but this one stores 

#include "Adt.hpp"

class PartSubSurface {
public:

  int8 *ipart_ist_and_bits; // st_global-and-bits-of-st -- back reference to full surface...
  int8 *ipart_isp_and_bits; // sp_global-and-bits-of-sp -- back reference to full surface...
  
  //int8 nsp_global,nst_global; // TODO: why?

  double bbmin[3];
  double bbmax[3];

  int nsp,nst,nse;
  double (*xsp)[3]; // local node copy - note that periodicity is no longer important
  int (*spost)[3]; // local sp-of-st
  int (*seost)[3]; // local se-of-st
  int *znost;

  Adt<double> * adt;

  int *sp_flag;
  int *se_flag;

  PartSubSurface() {
    
    //nsp_global = nst_global = 0;
    nsp = nst = nse = 0;
    
    xsp = NULL;
    ipart_ist_and_bits = NULL;
    ipart_isp_and_bits = NULL;
    spost = NULL;
    seost = NULL;
    znost = NULL;
    
    adt = NULL;

    sp_flag = NULL; // this is used for fast node matching
    se_flag = NULL;

  }

  void clear() {

    //nsp_global = nst_global = 0;
    nsp = nst = nse = 0;

    DELETE(xsp);
    DELETE(ipart_ist_and_bits);
    DELETE(ipart_isp_and_bits);
    DELETE(spost);
    DELETE(seost);
    DELETE(znost);
    
    if (adt != NULL) {
      delete adt;
      adt = NULL;
    }

    DELETE(sp_flag);
    DELETE(se_flag);

  }

  ~PartSubSurface() {

    clear();
    
  }
  
  void setIpartIstAndBits(const int ipart,const int ist,const int bits,const int ist_ss) {
    assert(ipart_ist_and_bits);
    assert((ipart >= 0)&&(ipart < PartData::partVec.size()));
    assert((ist >= 0)&&(ist < PartData::partVec[ipart]->surface->nst));
    assert((bits >= 0)&&(bits < (1<<6)));
    assert((ist_ss >= 0)&&(ist_ss < nst));
    ipart_ist_and_bits[ist_ss] = (int8(bits)<<(52)) | (int8(ipart)<<(32)) | ist;
    
    // check...
    int ipart_check,ist_check,bits_check;
    getIpartIstAndBits(ipart_check,ist_check,bits_check,ist_ss);
    assert(ipart_check == ipart);
    assert(ist_check == ist);
    assert(bits_check == bits);
  }
  
  void getIpartIstAndBits(int& ipart,int& ist,int& bits,const int ist_ss) const {
    assert(ipart_ist_and_bits);
    assert((ist_ss >= 0)&&(ist_ss < nst));
    bits = (ipart_ist_and_bits[ist_ss]>>(52));
    assert((bits >= 0)&&(bits < (1<<6)));
    ipart = ((ipart_ist_and_bits[ist_ss]>>(32))&MASK_20BITS);
    assert((ipart >= 0)&&(ipart < PartData::partVec.size()));
    ist = (ipart_ist_and_bits[ist_ss]&MASK_32BITS);
    assert((ist >= 0)&&(ist < PartData::partVec[ipart]->surface->nst));
  }

  void setIpartIspAndBits(const int ipart,const int isp,const int bits,const int isp_ss) {
    assert(ipart_isp_and_bits);
    assert((ipart >= 0)&&(ipart < PartData::partVec.size()));
    assert((isp >= 0)&&(isp < PartData::partVec[ipart]->surface->nsp));
    assert((bits >= 0)&&(bits < (1<<6)));
    assert((isp_ss >= 0)&&(isp_ss < nsp));
    ipart_isp_and_bits[isp_ss] = (int8(bits)<<(52)) | (int8(ipart)<<(32)) | isp;
    
    // check...
    int ipart_check,isp_check,bits_check;
    getIpartIspAndBits(ipart_check,isp_check,bits_check,isp_ss);
    assert(ipart_check == ipart);
    assert(isp_check == isp);
    assert(bits_check == bits);
  }
  
  void getIpartIspAndBits(int& ipart,int& isp,int& bits,const int isp_ss) const {
    assert(ipart_isp_and_bits);
    assert((isp_ss >= 0)&&(isp_ss < nsp));
    bits = (ipart_isp_and_bits[isp_ss]>>(52));
    assert((bits >= 0)&&(bits < (1<<6)));
    ipart = ((ipart_isp_and_bits[isp_ss]>>(32))&MASK_20BITS);
    assert((ipart >= 0)&&(ipart < PartData::partVec.size()));
    isp = (ipart_isp_and_bits[isp_ss]&MASK_32BITS);
    assert((isp >= 0)&&(isp < PartData::partVec[ipart]->surface->nsp));
  }

  void writeTecplot() {

    // not parallel...
    
    char filename[128];
    sprintf(filename,"partSubSurface.%06d.dat",mpi_rank);
    
    FILE * fp = fopen(filename,"w");
    fprintf(fp,"TITLE = \"%s\"\n",filename);
    fprintf(fp,"VARIABLES = \"X\"\n");
    fprintf(fp,"\"Y\"\n");
    fprintf(fp,"\"Z\"\n");
    
    fprintf(fp,"ZONE T=\"part sub surface\"\n");
    fprintf(fp,"N=%d, E=%d, F=FEPOINT, ET=TRIANGLE\n",nsp,nst);
    
    // data...
    for (int isp = 0; isp < nsp; ++isp)
      fprintf(fp,"%lf %lf %lf\n",xsp[isp][0],xsp[isp][1],xsp[isp][2]);
    
    for (int ist = 0; ist < nst; ++ist)
      fprintf(fp,"%d %d %d\n",
	      spost[ist][0]+1,
	      spost[ist][1]+1,
	      spost[ist][2]+1);
    
    fclose(fp);
    
  }

};

#endif

