#ifndef PHASEAVERAGING_HPP
#define PHASEAVERAGING_HPP

#include <string>
#include <cstdio>
#include "PngData.hpp"

// Helper functions for ImageTools Phase Averaging Routines

void writePhaseBinnedImages(const int mybin, const int dataFlag,
                            const double *buf, const double *count_buf,  const double *depth_buf,
                            const string& name, const string& varId, const string& varNm, const std::string& outputPathAndPrefix,
                            const int imgNx, const int imgNy, ImageMetadata *metadata,
                            const double range[],          const bool b_range,
                            const double rangeSurface[],   const bool b_rangeSurface,
                            const double rangeParticles[], const bool b_rangeParticles,
                            const double rangeIso[],       const bool b_rangeIso)
{
  assert(imgNx > 0 && imgNy > 0);
  PngData *imgOut = new PngData(imgNx, imgNy);

  // Build metadata
  ImageMetadata imdOut;
  imdOut.setTime(double(mybin));
  imdOut.transformMat[0]  = metadata->transformMat[0];
  imdOut.transformMat[1]  = metadata->transformMat[1];
  imdOut.transformMat[2]  = metadata->transformMat[2];
  imdOut.transformMat[3]  = metadata->transformMat[3];
  imdOut.transformMat[4]  = metadata->transformMat[4];
  imdOut.transformMat[5]  = metadata->transformMat[5];
  imdOut.transformMat[6]  = metadata->transformMat[6];
  imgOut->b_data = imgOut->b_dataSurface = imgOut->b_dataParticles = imgOut->b_dataIso = false;
  if(dataFlag == PngData::VOLUME_DATA_PIXEL)
  {
    imgOut->b_data = true;
    imgOut->setColorMap(metadata->getColorMapName());
    imdOut.setVarId(varId + "_" + varNm);
  }
  else if(dataFlag == PngData::SURFACE_DATA_PIXEL)
  {
    imgOut->b_dataSurface = true;
    imgOut->setColorMapSurface(metadata->getSurfColorMapName());
    imdOut.setVarOnSurfaceId(varId + "_" + varNm);
  }
  else if(dataFlag==PngData::PARTICLE_DATA_PIXEL)
  {
    imgOut->b_dataParticles = true;
    imgOut->setColorMapParticles(metadata->getPartColorMapName());
    imdOut.setVarOnParticleId(varId + "_" + varNm);
  }
  else if(dataFlag == PngData::ISO_DATA_PIXEL)
  {
    imgOut->b_dataIso = true;
    imgOut->setColorMapIso(metadata->getIsoColorMapName());
    imdOut.setVarOnIsoId(varId + "_" + varNm);
  }

  // Build image data
  int nPx = imgOut->getNx()*imgOut->getNy(); imgOut->npx = 0;
  for(int iPx=0; iPx < nPx; ++iPx)
  {
    if(count_buf[iPx] > 0)
    {
      imgOut->pixel_flag[iPx] = dataFlag; imgOut->npx++;
      imgOut->pixel_data[iPx] = buf[iPx];
      imgOut->setDepth(iPx, uint2(depth_buf[iPx]));
    }
    else
    {
      imgOut->pixel_flag[iPx] = PngData::BACKGROUND_PIXEL; // turn off this pixel
    }
  }

  char filename[256];
  sprintf(filename,"%s.%s.%s-%03d.png",outputPathAndPrefix.c_str(),name.c_str(),varNm.c_str(),mybin);
  cout << " > writing phase " << varNm << "(datatype " << name <<") to \"" << filename << "\"..." << endl;
  imgOut->rescaleRangesToData();
  if(b_range)          imgOut->setRange(range);
  if(b_rangeSurface)   imgOut->setRangeSurface(rangeSurface);
  if(b_rangeParticles) imgOut->setRangeParticles(rangeParticles);
  if(b_rangeIso)       imgOut->setRangeIso(rangeIso);
  imgOut->initializeWrite(imdOut);
  imgOut->write(filename);

  delete imgOut;
}

// Set the cyclic period index for the bin buffers for phase averaging
void setBinBufIndex(const int bindex, const int nbin, const int nbuf, std::vector<int>& bin2buf, std::vector<int>& buf2bin)
{
  bin2buf.clear(); bin2buf.resize(nbin); buf2bin.clear(); buf2bin.resize(nbuf);
  for(int ibin=0; ibin < nbin; ++ibin)
  {
    bin2buf[ibin] = (ibin-bindex)+(nbuf/2);
    if(bin2buf[ibin] < 0)
      bin2buf[ibin] = (+nbin+ibin-bindex)+(nbuf/2);
    else if(bin2buf[ibin] >= nbuf)
      bin2buf[ibin] = (-nbin+ibin-bindex)+(nbuf/2);
    if(bin2buf[ibin] >= 0 && bin2buf[ibin] < nbuf)
      buf2bin[bin2buf[ibin]] = ibin;
    else
      bin2buf[ibin] = -1;
  }
}

// Different phase-averaging methods can be used to compute various statistics of an image pixel with calcPxBinAvg
enum PhaseMethod { BIN, WEIGHT, LINEAR, SPLINE, FOURIER };

double calcPxBinAvg(const int px_buf, const double px_time, const int nbuf,
                    const double *px_mean_buf, const double *px_count_buf, const double *px_time_buf,
                    const PhaseMethod pmethod)
{
  double val = px_mean_buf[px_buf];
  if(pmethod == LINEAR)
  {
    if(     px_time < px_time_buf[px_buf] && (px_buf-1 >= 0 && px_buf-1 < nbuf) && px_count_buf[px_buf-1] > 0)
      val = px_mean_buf[px_buf-1] +
              (px_time - px_time_buf[px_buf-1])*
              (px_mean_buf[px_buf  ] - px_mean_buf[px_buf-1])/(px_time_buf[px_buf  ] - px_time_buf[px_buf-1]);
    else if(px_time > px_time_buf[px_buf] && (px_buf+1 >= 0 && px_buf+1 < nbuf) && px_count_buf[px_buf+1] > 0)
      val = px_mean_buf[px_buf  ] +
              (px_time - px_time_buf[px_buf  ])*
              (px_mean_buf[px_buf+1] - px_mean_buf[px_buf  ])/(px_time_buf[px_buf+1] - px_time_buf[px_buf  ]);
  }

  return val;
}

double img_getTime(const string& filename)
{
  double time = 0.0;
  PngImage img; img.read(filename.c_str(), false);
  if(img.isMetadataSet())
  {
    time = (img.getMetadata())->getTime();
  }

  return time;
}



#endif
