#include "Common.hpp"
#include "PngImage.hpp"
#include "MiscUtils.hpp"

void writeSampleData(unsigned char (*rgb)[3],PngDataChunk &zoNe,PngDataChunk &dpTh,PngDataChunk &daTa,const float range[2],const int ni, const int nj, const int i_c,const int j_c, const double &ang0, const int r2_max,const int zoneId,const int dpthId){
  unsigned char zoNe_uc[2], dpTh_uc[2], daTa_uc;
  zoNe_uc[0] = (uint2)zoneId & 255; 
  zoNe_uc[1] = (uint2)zoneId >> 8;
  dpTh_uc[0] = (uint2)dpthId & 255; //arbitrary value
  dpTh_uc[1] = (uint2)dpthId >> 8;

  int mod_ang = (int) floor(ang0/(M_PI*2.0));
  double ang00 = ang0;
  if (mod_ang>0)
    ang00 -= mod_ang*2.0*M_PI;

  //cout << "ang00 " << ang00 << " ang0 " << ang0 << " mod_ang " << mod_ang << " (int) floor(ang0/(M_PI*2.0)) " << (int) floor(ang0/(M_PI*2.0)) << endl;

  for (int ij = 0; ij<ni*nj; ++ij){
    int i = ij%ni;
    int j = ij/ni;
    if (((i-i_c)*(i-i_c)+(j-j_c)*(j-j_c))<r2_max){
      float ang = atan2((j-j_c),-(i-i_c)) + M_PI; //0:2pi
      if (ang<ang00){
        ang += 2.0*M_PI - ang00;
      }
      else{
        ang -= ang00;
      }
      float phi = (ang - range[0])/(range[1]-range[0]);
      rgb[ij][0] = (unsigned char) (phi*255 + 0.5);
      rgb[ij][1] = (unsigned char) (phi*255 + 0.5);
      rgb[ij][2] = (unsigned char) (phi*255 + 0.5);
      dpTh.set(ij,dpTh_uc);
      zoNe.set(ij,zoNe_uc);
      daTa_uc = (unsigned char) (phi*255 + 0.5);
      daTa.set(ij,&daTa_uc);
    }
  }
}

int writeSingleReferenceImage(int argc,char * argv[]){
  cout << "writeImage.exe writeSingleReferenceImage()..." << endl;

  string filename;
  if (argc>2){
    filename = argv[2];
  }
  else{
    cerr << "No image output filename specified." << endl;
    return 1;
  }

  const int ni = 400;
  const int nj = 200;

  //fill buffers with data
  unsigned char (*rgb)[3] = new unsigned char[ni*nj][3];
  PngDataChunk dpTh("dpTh", ni, nj, 2);
  unsigned char dpTh_uc[2] = {(unsigned char)255,(unsigned char)255};
  PngDataChunk zoNe("zoNe", ni, nj, 2);
  unsigned char zoNe_uc[2] = {(unsigned char)255,(unsigned char)255};
  PngDataChunk daTa("daTa", ni, nj, 1);
  unsigned char daTa_uc = (unsigned char)0;

  //fill with background color
  for (int ij = 0; ij<ni*nj; ++ij){
    rgb[ij][0] = 73;
    rgb[ij][1] = 175;
    rgb[ij][2] = 205;
    dpTh.set(ij,dpTh_uc);
    zoNe.set(ij,zoNe_uc);
    daTa.set(ij,&daTa_uc);
  }

  ImageMetadata imd;
  //let width = ni;
  imd.setLengthScale(1.0);

  //fill with grayscale
  int r2_max = 2025;
  float range[2] = {0.0,2.0*M_PI};
  writeSampleData(rgb,zoNe,dpTh,daTa,range,ni,nj,3*ni/4,1*nj/4,0.0,r2_max,65534,10); //plane data
  imd.setColorMapName("GRAYSCALE_RGB");
  imd.setVarId("angle0");
  imd.setRangeMin(range[0]);
  imd.setRangeMax(range[1]);

  writeSampleData(rgb,zoNe,dpTh,daTa,range,ni,nj,1*ni/4,1*nj/4,M_PI/2.0,r2_max,65533,20); //particle data
  imd.setColorMapParticlesName("GRAYSCALE_RGB");
  imd.setVarOnParticleId("angle1");
  imd.setRangeOnParticleMin(range[0]);
  imd.setRangeOnParticleMax(range[1]);

  writeSampleData(rgb,zoNe,dpTh,daTa,range,ni,nj,1*ni/4,3*nj/4,M_PI,r2_max,65532,20); //iso data
  imd.setColorMapIsoName("GRAYSCALE_RGB");
  imd.setVarOnIsoId("angle2");
  imd.setRangeOnIsoMin(range[0]);
  imd.setRangeOnIsoMax(range[1]);

  writeSampleData(rgb,zoNe,dpTh,daTa,range,ni,nj,3*ni/4,3*nj/4,3*M_PI/2.0,r2_max,65531,20); //surface data  32768-65531
  imd.setColorMapSurfaceName("GRAYSCALE_RGB");
  imd.setVarOnSurfaceId("angle3");
  imd.setRangeOnSurfaceMin(range[0]);
  imd.setRangeOnSurfaceMax(range[1]);



  PngImage png(ni,nj,rgb);


  imd.transformMat[0]  = 1.0;
  imd.transformMat[1]  = 0.0;
  imd.transformMat[2]  = 0.0;
  imd.transformMat[3]  = 0.0;
  imd.transformMat[4]  = 0.0;
  imd.transformMat[5]  = 1.0;
  imd.transformMat[6]  = 0.0;
  imd.transformMat[7]  = 0.0;
  imd.transformMat[8]  = 0.0;
  imd.transformMat[9]  = 0.0;
  imd.transformMat[10] = 0.0;
  imd.transformMat[11] = 0.0;
  imd.transformMat[12] = 0.0;
  imd.transformMat[13] = 0.0;
  imd.transformMat[14] = 0.0;
  imd.transformMat[15] = 1.0;

  imd.setRgbMode(0);
  imd.setCamDepth(0.0);

  png.setMetadata(&imd);

  dpTh.compressChunk();
  png.addPngDataChunk(&dpTh);

  zoNe.compressChunk();
  png.addPngDataChunk(&zoNe);

  daTa.compressChunk();
  png.addPngDataChunk(&daTa);

  cout << "writing image " << filename << endl;
  png.write(filename.c_str());

  //cleanup buffers
  delete[] rgb; rgb = NULL;
 
  return 0;

}

int writeTemporalReferenceImages(){

  cout << "writeImage.exe writeTemporalReferenceImage()..." << endl;

  const int ni = 400;
  const int nj = 200;

  //store the mean here
  int * buf_cnt = new int[ni*nj];
  double * buf_avg = new double[ni*nj];
  double * buf_rms = new double[ni*nj];
  double * buf_min = new double[ni*nj];
  double * buf_max = new double[ni*nj];

  //fill buffers with data
  unsigned char (*rgb)[3] = new unsigned char[ni*nj][3];
  PngDataChunk dpTh("dpTh", ni, nj, 2);
  unsigned char dpTh_uc[2] = {(unsigned char)255,(unsigned char)255};
  PngDataChunk zoNe("zoNe", ni, nj, 2);
  unsigned char zoNe_uc[2] = {(unsigned char)255,(unsigned char)255};
  PngDataChunk daTa("daTa", ni, nj, 1);
  unsigned char daTa_uc = (unsigned char)0;

  //fill with background color
  for (int ij = 0; ij<ni*nj; ++ij){
    rgb[ij][0] = 73;
    rgb[ij][1] = 175;
    rgb[ij][2] = 205;
    dpTh.set(ij,dpTh_uc);
    zoNe.set(ij,zoNe_uc);
    daTa.set(ij,&daTa_uc);

    buf_cnt[ij] = 0;
    buf_avg[ij] = 0.0;
    buf_rms[ij] = 0.0;
    buf_min[ij] = HUGE_VAL;
    buf_max[ij] = -HUGE_VAL;
  }
 
  //Setup image template
  float range[2] = {0.0,2.0*M_PI};
  ImageMetadata imd;
  //let width = ni;
  imd.setLengthScale(1.0);
  imd.setColorMapName("GRAYSCALE_RGB");
  imd.setVarId("angle0");
  imd.setRangeMin(range[0]);
  imd.setRangeMax(range[1]);
  imd.setColorMapParticlesName("GRAYSCALE_RGB");
  imd.setVarOnParticleId("angle1");
  imd.setRangeOnParticleMin(range[0]);
  imd.setRangeOnParticleMax(range[1]);
  imd.setColorMapIsoName("GRAYSCALE_RGB");
  imd.setVarOnIsoId("angle2");
  imd.setRangeOnIsoMin(range[0]);
  imd.setRangeOnIsoMax(range[1]);
  imd.setColorMapSurfaceName("GRAYSCALE_RGB");
  imd.setVarOnSurfaceId("angle3");
  imd.setRangeOnSurfaceMin(range[0]);
  imd.setRangeOnSurfaceMax(range[1]);

  imd.transformMat[0]  = 1.0;
  imd.transformMat[1]  = 0.0;
  imd.transformMat[2]  = 0.0;
  imd.transformMat[3]  = 0.0;
  imd.transformMat[4]  = 0.0;
  imd.transformMat[5]  = 1.0;
  imd.transformMat[6]  = 0.0;
  imd.transformMat[7]  = 0.0;
  imd.transformMat[8]  = 0.0;
  imd.transformMat[9]  = 0.0;
  imd.transformMat[10] = 0.0;
  imd.transformMat[11] = 0.0;
  imd.transformMat[12] = 0.0;
  imd.transformMat[13] = 0.0;
  imd.transformMat[14] = 0.0;
  imd.transformMat[15] = 1.0;
  imd.setRgbMode(0);
  imd.setCamDepth(0.0);

  //End setup image template

  string prefix = "images_ref/temporal_";
  for (int iImg = 11; iImg>=0; --iImg){

    stringstream ss; 
    ss << prefix << std::setfill('0') << std::setw(2) << iImg << ".png";

    double ang_t = iImg * M_PI/6.0;

    //fill with grayscale
    int r2_max = 2025;
    writeSampleData(rgb,zoNe,dpTh,daTa,range,ni,nj,3*ni/4,1*nj/4,0.0+ang_t,r2_max,65534,10); //plane data
    writeSampleData(rgb,zoNe,dpTh,daTa,range,ni,nj,1*ni/4,1*nj/4,M_PI/2.0+ang_t,r2_max,65533,20); //particle data
    writeSampleData(rgb,zoNe,dpTh,daTa,range,ni,nj,1*ni/4,3*nj/4,M_PI+ang_t,r2_max,65532,20); //iso data
    writeSampleData(rgb,zoNe,dpTh,daTa,range,ni,nj,3*ni/4,3*nj/4,3*M_PI/2.0+ang_t,r2_max,65531,20); //surface data  32768-65531

    PngImage png(ni,nj,rgb);
    png.setMetadata(&imd);

    dpTh.compressChunk();
    png.addPngDataChunk(&dpTh);

    zoNe.compressChunk();
    png.addPngDataChunk(&zoNe);

    daTa.compressChunk();
    png.addPngDataChunk(&daTa);

    cout << "writing temporal image " << ss.str() << endl;
    MiscUtils::mkdir_for_file(ss.str());
    png.write(ss.str().c_str());

    //update avg, rms
    for (int ij = 0; ij<ni*nj; ++ij){
      unsigned char zoNe_uc[2];
      zoNe.get(ij,zoNe_uc);
      uint2 zoNe_us = (zoNe_uc[1] << 8 | zoNe_uc[0]); //little endian  
      if (zoNe_us<65535){
        uint8_t daTa_i8;
        daTa.get(ij,&daTa_i8);
        double val = daTa_i8/255.0 * (range[1]-range[0]) + range[0];

        ++buf_cnt[ij];
        buf_avg[ij] += val;
        buf_rms[ij] += val*val;

        if (val < buf_min[ij]) buf_min[ij] = val;
        if (val > buf_max[ij]) buf_max[ij] = val;
      } 
    }

  }
 
  //finalize stats
  for (int ij = 0; ij<ni*nj; ++ij){
    if (buf_cnt[ij]>0){
      buf_avg[ij] /= buf_cnt[ij];
      buf_rms[ij] = sqrt(max(0.0,buf_rms[ij]/buf_cnt[ij]-buf_avg[ij]*buf_avg[ij]));
    }
  }
  //hijack buf_cnt to flag pixels where lighting would not be reconstructed
  for (int ij = 0; ij<ni*nj; ++ij){
    if (buf_cnt[ij]>0){
      int pixel = rgb[ij][0];
      buf_cnt[ij] = pixel;
    }
    else
      buf_cnt[ij] = -1;
  }


  //compute ranges
  double r0min[4] = {HUGE_VAL,HUGE_VAL,HUGE_VAL,HUGE_VAL};
  double r1min[4] = {HUGE_VAL,HUGE_VAL,HUGE_VAL,HUGE_VAL};
  double r2min[4] = {HUGE_VAL,HUGE_VAL,HUGE_VAL,HUGE_VAL};
  double r3min[4] = {HUGE_VAL,HUGE_VAL,HUGE_VAL,HUGE_VAL};
  double r0max[4] = {-HUGE_VAL,-HUGE_VAL,-HUGE_VAL,-HUGE_VAL};
  double r1max[4] = {-HUGE_VAL,-HUGE_VAL,-HUGE_VAL,-HUGE_VAL};
  double r2max[4] = {-HUGE_VAL,-HUGE_VAL,-HUGE_VAL,-HUGE_VAL};
  double r3max[4] = {-HUGE_VAL,-HUGE_VAL,-HUGE_VAL,-HUGE_VAL};
  for (int ij = 0; ij<ni*nj; ++ij){
    zoNe.get(ij,zoNe_uc);
    uint2 zoNe_us = (zoNe_uc[1] << 8 | zoNe_uc[0]); //little endian  
    if (zoNe_us==65534){
      if (buf_avg[ij]<r0min[0])
        r0min[0] = buf_avg[ij];
      if (buf_rms[ij]<r0min[1])
        r0min[1] = buf_rms[ij];
      if (buf_min[ij]<r0min[2])
        r0min[2] = buf_min[ij];
      if (buf_max[ij]<r0min[3])
        r0min[3] = buf_max[ij];

      if (buf_avg[ij]>r0max[0])
        r0max[0] = buf_avg[ij];
      if (buf_rms[ij]>r0max[1])
        r0max[1] = buf_rms[ij];
      if (buf_min[ij]>r0max[2])
        r0max[2] = buf_min[ij];
      if (buf_max[ij]>r0max[3])
        r0max[3] = buf_max[ij];
    }
    else if (zoNe_us==65533){
      if (buf_avg[ij]<r1min[0])
        r1min[0] = buf_avg[ij];
      if (buf_rms[ij]<r1min[1])
        r1min[1] = buf_rms[ij];
      if (buf_min[ij]<r1min[2])
        r1min[2] = buf_min[ij];
      if (buf_max[ij]<r1min[3])
        r1min[3] = buf_max[ij];

      if (buf_avg[ij]>r1max[0])
        r1max[0] = buf_avg[ij];
      if (buf_rms[ij]>r1max[1])
        r1max[1] = buf_rms[ij];
      if (buf_min[ij]>r1max[2])
        r1max[2] = buf_min[ij];
      if (buf_max[ij]>r1max[3])
        r1max[3] = buf_max[ij];
    }
    else if (zoNe_us==65532){
      if (buf_avg[ij]<r2min[0])
        r2min[0] = buf_avg[ij];
      if (buf_rms[ij]<r2min[1])
        r2min[1] = buf_rms[ij];
      if (buf_min[ij]<r2min[2])
        r2min[2] = buf_min[ij];
      if (buf_max[ij]<r2min[3])
        r2min[3] = buf_max[ij];

      if (buf_avg[ij]>r2max[0])
        r2max[0] = buf_avg[ij];
      if (buf_rms[ij]>r2max[1])
        r2max[1] = buf_rms[ij];
      if (buf_min[ij]>r2max[2])
        r2max[2] = buf_min[ij];
      if (buf_max[ij]>r2max[3])
        r2max[3] = buf_max[ij];
    }
    else if (zoNe_us==65531){
      if (buf_avg[ij]<r3min[0])
        r3min[0] = buf_avg[ij];
      if (buf_rms[ij]<r3min[1])
        r3min[1] = buf_rms[ij];
      if (buf_min[ij]<r3min[2])
        r3min[2] = buf_min[ij];
      if (buf_max[ij]<r3min[3])
        r3min[3] = buf_max[ij];

      if (buf_avg[ij]>r3max[0])
        r3max[0] = buf_avg[ij];
      if (buf_rms[ij]>r3max[1])
        r3max[1] = buf_rms[ij];
      if (buf_min[ij]>r3max[2])
        r3max[2] = buf_min[ij];
      if (buf_max[ij]>r3max[3])
        r3max[3] = buf_max[ij];
    }
  }

  dpTh.compressChunk();
  zoNe.compressChunk();

  //write avg
  for (int ij = 0; ij<ni*nj; ++ij){
    if (buf_cnt[ij]>=0){
      zoNe.get(ij,zoNe_uc);
      uint2 zoNe_us = (zoNe_uc[1] << 8 | zoNe_uc[0]); //little endian  
      double phi;
      if (zoNe_us==65534){
        phi = (buf_avg[ij] - r0min[0])/(r0max[0]-r0min[0]);
      }
      else if (zoNe_us==65533){
        phi = (buf_avg[ij] - r1min[0])/(r1max[0]-r1min[0]);
      }
      else if (zoNe_us==65532){
        phi = (buf_avg[ij] - r2min[0])/(r2max[0]-r2min[0]);
      }
      else if (zoNe_us==65531){
        phi = (buf_avg[ij] - r3min[0])/(r3max[0]-r3min[0]);
      }
      else{
        phi = 0;
        assert(0);
      }

      if (buf_cnt[ij]==0 && (zoNe_us!=65534)){
        FOR_I3 rgb[ij][i] = 0;  //light cannot be reconstructed
      }
      else{
        FOR_I3 rgb[ij][i] = (unsigned char) (phi*255 + 0.5);
      }
      unsigned char daTa_uc = (unsigned char) (phi*255 + 0.5);
      daTa.set(ij,&daTa_uc);
    }
  }
  PngImage png_avg(ni,nj,rgb);
  daTa.compressChunk();
  png_avg.addPngDataChunk(&daTa);
  png_avg.addPngDataChunk(&dpTh);
  png_avg.addPngDataChunk(&zoNe);
  imd.setVarId("angle0_avg");
  imd.setVarOnParticleId("angle1_avg");
  imd.setVarOnIsoId("angle2_avg");
  imd.setVarOnSurfaceId("angle3_avg");
  imd.setRangeMin(r0min[0]);
  imd.setRangeMax(r0max[0]);
  imd.setRangeOnParticleMin(r1min[0]);
  imd.setRangeOnParticleMax(r1max[0]);
  imd.setRangeOnIsoMin(r2min[0]);
  imd.setRangeOnIsoMax(r2max[0]);
  imd.setRangeOnSurfaceMin(r3min[0]);
  imd.setRangeOnSurfaceMax(r3max[0]);
 
  png_avg.setMetadata(&imd);

  cout << "writing temporal image avg images_ref/stats_avg.png" << endl;
  png_avg.write("images_ref/stats_avg.png");

  //write rms
  for (int ij = 0; ij<ni*nj; ++ij){
    if (buf_cnt[ij]>=0){
      zoNe.get(ij,zoNe_uc);
      uint2 zoNe_us = (zoNe_uc[1] << 8 | zoNe_uc[0]); //little endian  
      double phi;
      if (zoNe_us==65534){
        phi = (buf_rms[ij] - r0min[1])/(r0max[1]-r0min[1]);
      }
      else if (zoNe_us==65533){
        phi = (buf_rms[ij] - r1min[1])/(r1max[1]-r1min[1]);
      }
      else if (zoNe_us==65532){
        phi = (buf_rms[ij] - r2min[1])/(r2max[1]-r2min[1]);
      }
      else if (zoNe_us==65531){
        phi = (buf_rms[ij] - r3min[1])/(r3max[1]-r3min[1]);
      }
      else{
        phi = 0;
        assert(0);
      }

      if (buf_cnt[ij]==0 && (zoNe_us!=65534)){
        FOR_I3 rgb[ij][i] = 0;  //light cannot be reconstructed
      }
      else{
        FOR_I3 rgb[ij][i] = (unsigned char) (phi*255 + 0.5);
      }
      unsigned char daTa_uc = (unsigned char) (phi*255 + 0.5);
      daTa.set(ij,&daTa_uc);
    }
  }
  PngImage png_rms(ni,nj,rgb);
  daTa.compressChunk();
  png_rms.addPngDataChunk(&daTa);
  png_rms.addPngDataChunk(&dpTh);
  png_rms.addPngDataChunk(&zoNe);
  imd.setVarId("angle0_rms");
  imd.setVarOnParticleId("angle1_rms");
  imd.setVarOnIsoId("angle2_rms");
  imd.setVarOnSurfaceId("angle3_rms");
  imd.setRangeMin(r0min[1]);
  imd.setRangeMax(r0max[1]);
  imd.setRangeOnParticleMin(r1min[1]);
  imd.setRangeOnParticleMax(r1max[1]);
  imd.setRangeOnIsoMin(r2min[1]);
  imd.setRangeOnIsoMax(r2max[1]);
  imd.setRangeOnSurfaceMin(r3min[1]);
  imd.setRangeOnSurfaceMax(r3max[1]);
 
  png_rms.setMetadata(&imd);

  cout << "writing temporal image avg images_ref/stats_rms.png" << endl;
  png_rms.write("images_ref/stats_rms.png");

  //write min
  for (int ij = 0; ij<ni*nj; ++ij){
    if (buf_cnt[ij]>=0){
      zoNe.get(ij,zoNe_uc);
      uint2 zoNe_us = (zoNe_uc[1] << 8 | zoNe_uc[0]); //little endian  
      double phi;
      if (zoNe_us==65534){
        phi = (buf_min[ij] - r0min[2])/(r0max[2]-r0min[2]);
      }
      else if (zoNe_us==65533){
        phi = (buf_min[ij] - r1min[2])/(r1max[2]-r1min[2]);
      }
      else if (zoNe_us==65532){
        phi = (buf_min[ij] - r2min[2])/(r2max[2]-r2min[2]);
      }
      else if (zoNe_us==65531){
        phi = (buf_min[ij] - r3min[2])/(r3max[2]-r3min[2]);
      }
      else{
        phi = 0;
        assert(0);
      }

      if (buf_cnt[ij]==0 && (zoNe_us!=65534)){
        FOR_I3 rgb[ij][i] = 0;  //light cannot be reconstructed
      }
      else{
        FOR_I3 rgb[ij][i] = (unsigned char) (phi*255 + 0.5);
      }
      unsigned char daTa_uc = (unsigned char) (phi*255 + 0.5);
      daTa.set(ij,&daTa_uc);
    }
  }
  PngImage png_min(ni,nj,rgb);
  daTa.compressChunk();
  png_min.addPngDataChunk(&daTa);
  png_min.addPngDataChunk(&dpTh);
  png_min.addPngDataChunk(&zoNe);
  imd.setVarId("angle0_min");
  imd.setVarOnParticleId("angle1_min");
  imd.setVarOnIsoId("angle2_min");
  imd.setVarOnSurfaceId("angle3_min");
  imd.setRangeMin(r0min[2]);
  imd.setRangeMax(r0max[2]);
  imd.setRangeOnParticleMin(r1min[2]);
  imd.setRangeOnParticleMax(r1max[2]);
  imd.setRangeOnIsoMin(r2min[2]);
  imd.setRangeOnIsoMax(r2max[2]);
  imd.setRangeOnSurfaceMin(r3min[2]);
  imd.setRangeOnSurfaceMax(r3max[2]);
 
  png_min.setMetadata(&imd);

  cout << "writing temporal image avg images_ref/stats_min.png" << endl;
  png_min.write("images_ref/stats_min.png");

  //write max
  for (int ij = 0; ij<ni*nj; ++ij){
    if (buf_cnt[ij]>=0){
      zoNe.get(ij,zoNe_uc);
      uint2 zoNe_us = (zoNe_uc[1] << 8 | zoNe_uc[0]); //little endian  
      double phi;
      if (zoNe_us==65534){
        phi = (buf_max[ij] - r0min[3])/(r0max[3]-r0min[3]);
      }
      else if (zoNe_us==65533){
        phi = (buf_max[ij] - r1min[3])/(r1max[3]-r1min[3]);
      }
      else if (zoNe_us==65532){
        phi = (buf_max[ij] - r2min[3])/(r2max[3]-r2min[3]);
      }
      else if (zoNe_us==65531){
        phi = (buf_max[ij] - r3min[3])/(r3max[3]-r3min[3]);
      }
      else{
        phi = 0;
        assert(0);
      }

      if (buf_cnt[ij]==0 && (zoNe_us!=65534)){
        FOR_I3 rgb[ij][i] = 0;  //light cannot be reconstructed
      }
      else{
        FOR_I3 rgb[ij][i] = (unsigned char) (phi*255 + 0.5);
      }
      unsigned char daTa_uc = (unsigned char) (phi*255 + 0.5);
      daTa.set(ij,&daTa_uc);
    }
  }
  PngImage png_max(ni,nj,rgb);
  daTa.compressChunk();
  png_max.addPngDataChunk(&daTa);
  png_max.addPngDataChunk(&dpTh);
  png_max.addPngDataChunk(&zoNe);
  imd.setVarId("angle0_max");
  imd.setVarOnParticleId("angle1_max");
  imd.setVarOnIsoId("angle2_max");
  imd.setVarOnSurfaceId("angle3_max");
  imd.setRangeMin(r0min[3]);
  imd.setRangeMax(r0max[3]);
  imd.setRangeOnParticleMin(r1min[3]);
  imd.setRangeOnParticleMax(r1max[3]);
  imd.setRangeOnIsoMin(r2min[3]);
  imd.setRangeOnIsoMax(r2max[3]);
  imd.setRangeOnSurfaceMin(r3min[3]);
  imd.setRangeOnSurfaceMax(r3max[3]);
 
  png_max.setMetadata(&imd);
 
  cout << "writing temporal image avg images_ref/stats_max.png" << endl;
  png_max.write("images_ref/stats_max.png");

  //cleanup buffers
  delete[] rgb; rgb = NULL;
  delete[] buf_cnt; buf_cnt = NULL;
  delete[] buf_avg; buf_avg = NULL;
  delete[] buf_rms; buf_rms = NULL;
  delete[] buf_min; buf_min = NULL;
  delete[] buf_max; buf_max = NULL;

  return 0;
}

int main(int argc,char * argv[]) {
 
  cout << "writeImage.exe main() ..." << endl;

  string flag;
  if (argc>1){
    flag = argv[1];
  }
  else{
    cerr << "No image output flag specified" << endl;
    return 1;
  }
 
  if (flag=="--SINGLE"){
    return writeSingleReferenceImage(argc,argv);
  }
  else if (flag=="--TEMPORAL"){
    return writeTemporalReferenceImages();
  }
  else {
    cerr << "Unrecognized image output flag" << endl;
    return 1;
  }

  cout << "Regression Logic Error" << endl; //shouldn't get here
  return 1; 

}
