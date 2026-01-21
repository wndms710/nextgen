#ifndef IMAGETOOLS_HPP
#define IMAGETOOLS_HPP

#ifdef WITH_X264
#include "x264.h"
#endif

#include "PngCache.hpp"
#include "PngData.hpp"
#include "MiscUtils.hpp"
#include "PhaseAveraging.hpp"

#include "CtiType.hpp"

void addTextCentered(unsigned char (*rgb)[3],const int nx,const int ny,const string& theText,const int i0,const int j0,const int ifs) {

  // this is a simple example of how to add black-on-white text to an image.
  
  assert((ifs >= 0)&&(ifs < cti_type_n));
  
  const int length = theText.length();
  const int ishift = (length*cti_type_width[ifs])/2;
  
  for (int ii = 0; ii < theText.length(); ++ii) {
    const int ichar = int(theText[ii])-32;
    assert((ichar >= 0)&&(ichar < cti_type_count));
    for (int di = 0; di < cti_type_width[ifs]; ++di) {
      for (int dj = 0; dj < cti_type_height[ifs]; ++dj) {
        const int i = i0+ii*cti_type_width[ifs]+di-ishift;
        const int j = j0+dj;
        if ((i >= 0)&&(i < nx)&&(j >= 0)&&(j < ny)) {
          const int ipx = j*nx + i;
          // we could blend here, but for now, just push in black on white text...
          const unsigned char uc = cti_type_data_ptr[ifs][ichar*cti_type_width[ifs]*cti_type_height[ifs]+dj*cti_type_width[ifs]+di];
          rgb[ipx][0] = uc;
          rgb[ipx][1] = uc;
          rgb[ipx][2] = uc;
        }
      }
    }
  }
  
}

class ImageText {
  // this class is used to manage image text that can be added
  // using the TEXT param (currently only used in tile)...
  // for example:
  // TEXT <i> <j> [SIZE <ifs>] The text to add
  // TEXT FRAC <i-frac> <j-frac> SIZE <ifs> The text to add
public:
  string str;
  int i,j,ifs;
  bool b_frac; // used to indicate "frac" mode: i = ifrac*nx, j = jfrac*ny
  double ifrac,jfrac;
  ImageText() {
    i = -1;
    j = -1;
    b_frac = false;
    ifrac = -1.0;
    jfrac = -1.0;
    ifs = -1;
  }
};

void processTextParams(vector<ImageText>& textVec) {
  FOR_PARAM_MATCHING("TEXT") {
    // expect lines like... 
    // TEXT <i> <j> [SIZE <ifs>] The text to add
    // TEXT FRAC <i-frac> <j-frac> SIZE <ifs> The text to add
    // to better support this in the future, consider moving this out of this 
    // routine, supporting double-quote string parsing in the params, and
    // moving to a surfer-like interpreted param model here in ping...
    textVec.push_back(ImageText());
    ImageText * it = &textVec.back();
    int iarg = 0;
    if (param->getString(iarg) == "FRAC") {
      ++iarg;
      it->b_frac = true;
      it->ifrac = param->getDouble(iarg++);
      it->jfrac = param->getDouble(iarg++);
    }
    else {
      it->i = param->getInt(iarg++);
      it->j = param->getInt(iarg++);
    }
    if (param->getString(iarg) == "SIZE") {
      ++iarg;
      it->ifs = param->getInt(iarg++);
    }
    else {
      it->ifs = 8; // good choice?
    }
    // the rest is the text. Since we cannot currently process "a text entry"
    // as a single string param in the params, just rebuild it here with spaces.
    // This will mess up multi-spaces or the use of "=" but works ok for now...
    bool first = true;
    while (iarg < param->size()) {
      if (first) {
        it->str = param->getString(iarg++);
        first = false;
      }
      else {
        it->str = it->str + " " + param->getString(iarg++);
      }
    }
  }
}

class ImageTools {

protected:

  PngCache ** pngCaches;
  PngCache * pngCache;

public:

  template <typename T> std::string NumberToString (T Number) {
    std::ostringstream ss;
    ss << Number;
    return ss.str();
  }

  ImageTools() {
    pngCaches = NULL;
    pngCache = NULL;
  }

  ~ImageTools(){
    if (pngCache!=NULL){
      delete pngCache;
      pngCache = NULL;
    }
    if (pngCaches != NULL) {
      DELETE(pngCaches);
      pngCaches = NULL;
    }
    ColorMap::clearColorMaps();
  }

  // -------------------------------
  // BEGIN Parse
  // -------------------------------
  /*
    void parseAll(const string &paramName){
    FOR_PARAM_MATCHING(paramName)
    parse(&(*param));
    }

    void parse(const string &paramName){
    if(Param * param = getParam(paramName))
    parse(param);
    }
  */

  //Call the matching function for a given param
  void parse(Param * param){
    COUT2("Working on param \"" << param->str() << "\"");
    const string keyword = param->getName();
    if ((keyword == "IMAGE") || (keyword == "IMAGES")) {
      itpImages(param);
    }
    else if ((keyword == "IMAGELIST") || (keyword == "IMAGE_LIST")) {
      itpImagelist(param);
    }
    else if ((keyword == "ADD_COLORMAP") || (keyword=="ADD_COLORMAP")){
      itpAddColorMap(param);
    }
    else if ((keyword == "SPATIALSTATS") || (keyword == "SPATIAL_STATS")) {
      itpSpatialStats(param);
    }
    else if (keyword == "SPATIAL_STATS_DISK") {
      itpSpatialStatsDisk(param);
    }
    else if (keyword == "RADIAL_PROFILE") {
      itpRadialProfile(param);
    }
    else if (keyword == "ANNULAR_INTEGRATION_PROBE") {
      itpAnnularIntegrationProbe(param);
    }
    else if (keyword == "TUMBLE") {
      itpTumble(param);
    }
    // NOTE: could condense the following data writers into just WRITE_DATA
    // and add optional parameter for N the number of image types...
    else if ((keyword == "WRITE_VELOCITY")||(keyword == "WRITE_VECTOR")) {
      itpWriteVelocity(param);
    }
    else if ((keyword == "WRITE_SCALAR")||(keyword == "WRITE_DATA")) {
      itpWriteScalar(param);
    }
    else if (keyword == "VECTOR") {
      itpVector(param);
    }
    else if (keyword == "CROP") {
      itpCrop(param);
    }
    else if ((keyword == "INTEGRATIONPROBE") || (keyword == "INTEGRATION_PROBE")) {
      itpIntegrationProbe(param);
    }
    else if ((keyword == "RADIALPROBE") || (keyword == "RADIAL_PROBE")) {
      itpRadialProbe(param);
    }
    else if ((keyword == "TEMPORALSTATS") || (keyword == "TEMPORAL_STATS") || (keyword == "TS")) {
      itpTemporalStats(param);
    }
    else if ((keyword == "TEMPORALSTATSCROSS") || (keyword == "TEMPORAL_STATS_CROSS") || (keyword == "TSC")) {
      itpTemporalStatsCross(param);
    }
    else if ((keyword == "TEMPORALSTATSMOVING") || (keyword == "TEMPORAL_STATS_MOVING") || (keyword == "TSM")) {
      itpTemporalStatsMoving(param);
    }
    else if ((keyword == "TEMPORALSTATSPHASE") || (keyword == "TEMPORAL_STATS_PHASE") || (keyword == "TSP")) {
      itpTemporalStatsPhase(param);
    }
    else if ((keyword == "TEMPORALSTATSPHASELOCK") || (keyword == "TEMPORAL_STATS_PHASE_LOCK") || (keyword == "TSPL")) {
      itpTemporalStatsPhaseLock(param);
    }
    else if (keyword == "DIFF") {
      itpDiff(param);
    }
    else if (keyword == "DIFF_STATS") {
      itpDiffStats(param);
    }
    else if (keyword == "RESCALE") {
      itpRescale(param);
    }
    else if (keyword == "RECOLOR") {
      itpRecolor(param);
    }
    else if (keyword == "MARK_MEDIAL_AXIS"){
      itpMarkMedialAxis(param);
    }
    else if (keyword == "TILE") {
      itpTile(param);
    }
    else if (keyword == "PROBE") {
      itpProbe(param);
    }
    else if (keyword == "NPROBE") {
      processNprobe(param);
    }
    else if (keyword == "PROJECTED_SURFACE_AREA"){
      itpProjectedSurfaceArea(param);
    }
    else if ((keyword == "COLLAPSEHEIGHT") || (keyword == "COLLAPSE_HEIGHT")) {
      itpCollapseHeight(param);
    }
    else if ((keyword == "COLLAPSEWIDTH") || (keyword == "COLLAPSE_WIDTH")) {
      itpCollapseWidth(param);
    }
    else if ((keyword == "DUMPDATA") || (keyword == "DUMP_DATA")) {
      itpDumpData(param);
    }
    else if ((keyword == "DUMPCOLUMN") || (keyword == "DUMP_COLUMN")) {
      itpDumpColumn(param);
    }
    else if ((keyword == "1DPDF") || (keyword == "1D_PDF")) {
      itp1DPdf(param);
    }
    else if ((keyword == "2DPDF") || (keyword == "2D_PDF")) {
      itp2DPdf(param);
    }
    else if ((keyword == "DIFFSET") || (keyword == "DIFF_SET")) {
      itpDiffSet(param);
    }
    else if ((keyword == "DISKAVERAGE") || (keyword == "DISCAVERAGE") || (keyword == "DISK_AVERAGE") || (keyword == "DISC_AVERAGE")) {
      itpDiscAverage(param);
    }
    else if ((keyword == "NOTDISKAVERAGE") || (keyword == "NOTDISCAVERAGE") || (keyword == "NOT_DISK_AVERAGE") || (keyword == "NOT_DISC_AVERAGE")) {
      itpNotDiscAverage(param);
    }
    else if (keyword == "ISOCONTOUR") {
      itpIsocontours(param, 0);
    }
    else if (keyword == "ISOCONTOURS") {
      itpIsocontours(param, 1);
    }
    else if (keyword == "SUPERPOSE") {
      itpSuperpose(param);
    }
    else if (keyword == "INFO"){ //print image metadata information
      itpInfo(param);
    }
    else if (keyword == "ENCODE"){
      itpEncode(param);
    }
    else if (keyword == "GLOBAL_MASK"){
      itpGlobalMask(param);
    }
    else if ((keyword == "CTILICENSE") || (keyword == "CTI_LICENSE")) {
      COUT1("license location set by input file");
    }
    else if (keyword == "VERBOSE") {
      COUT1("log verbosity set");
    }
    else if (keyword == "HELP") {
      itpHelp(param);
    }
    else {
      CERR("unrecognized ping param \"" << keyword << "\"");
    }
    COUT2("Done working on param \'" << param->str() << "\'");
  }
  // -------------------------------
  // END Parse
  // -------------------------------

  void itpImages(Param * param){
    if (pngCache != NULL)
      delete pngCache;

    //  remaining tokens in the parameter should be image filenames...
    pngCache = new PngCache();
    int iImg = 0;
    while (iImg < param->size()) {
      pngCache->imageList.push_back(param->getString(iImg));
      ++iImg;
    }
    if (iImg > 0){
      pngCache->init();
    }

    if (pngCache->imageList.empty()) {
      CERR("pngCache is empty. Param IMAGES requires one or image filenames.");
    }
  }//itpImages()

  void itpImagelist(Param * param) {
    if (pngCache != NULL) delete pngCache;

    if (param->size() != 1) {
      CERR("Param IMAGELIST requires the name of a file containing the list of images");
    }

    string imageListName = param->getString(0);
    pngCache = new PngCache(imageListName);

    if (pngCache->imageList.empty()) {
      CERR("pngCache is empty. Check format of list file \"" << imageListName << "\"");
    }
  }//itpImageList

  void itpAddColorMap(Param * param){
    COUT1("ADD_COLORMAP...");
    //Syntax COLORMAP FILE=<filename.dat>

    if (param->size() < 2 || param->getString(0)!="FILE"){
      CERR("problem parsing ADD_COLORMAP; specify ADD_COLORMAP FILE=<filename.dat>");
    }
    string filename = param->getString(1);

    ColorMap::addUserDefinedColorMap(filename);

  }

  void makeMask(int * mask,const int mask_i0,const int mask_j0,const int mask_ni,const int mask_nj,const vector<pair<int,int> >& mask_ij) {
    // for each pixel in the mask, count the intersections with the polygons...
    for (int i0 = 0; i0 < mask_ni; ++i0) {
      for (int j0 = 0; j0 < mask_nj; ++j0) {
        int count = 0;
        const int i = i0 + mask_i0;
        const int j = j0 + mask_j0;
        // now loop through polylines and count intersections for line from 0,0 to i,j...
        for (int ii = 0; ii < mask_ij.size(); ++ii) {
          const int is0 = mask_ij[ii].first;
          const int js0 = mask_ij[ii].second;
          int is1,js1;
          if (ii+1 < mask_ij.size()) {
            is1 = mask_ij[ii+1].first;
            js1 = mask_ij[ii+1].second;
          }
          else {
            // loop back to first...
            is1 = mask_ij[0].first;
            js1 = mask_ij[0].second;
          }
          if ((min(is0,is1) <= i)&&(min(js0,js1) <= j)) {
            // intersection is possible...
            const int8 As0 = int8(j)*int8(is0) - int8(i)*int8(js0);
            const int8 As1 = int8(j)*int8(is1) - int8(i)*int8(js1);
            if ((As0 == 0)&&(As1 == 0)) {
              cout << "zero area!" << endl;
              assert(0);
            }
            else if (((As0 <= 0)&&(As1 >= 0))||((As0 >= 0)&&(As1 <= 0))) {
              // Areas are opposite sign, so potential intersection...
              const int A0 = int8(js1-js0)*int8(is0) - int8(is1-is0)*int8(js0);
              const int A1 = int8(js1-js0)*int8(is0-i) - int8(is1-is0)*int8(js0-j);
              if ((A0 == 0)&&(A1 == 0)) {
                cout << "zero area!" << endl;
                assert(0);
              }
              else if (((A0 <= 0)&&(A1 >= 0))||((A0 >= 0)&&(A1 <= 0))) {
                // Areas are also opposite sign, so this is an intersection...
                if (As0 != 0) ++count;
                if (As1 != 0) ++count;
              }
            }
          }
        }
        // we add intersections in pairs, so they will always be multiples of 2...
        assert(count%2 == 0);
        // but is divisible by 4, then this is out...
        if (count%4 == 0) {
          mask[i0*mask_nj+j0] = 0;
        }
        else {
          mask[i0*mask_nj+j0] = 1;
        }
      }
    }
        
    // write the mask to check...
    /*
      {
      FILE * fp = fopen("mask.dat","w");
      for (int i0 = 0; i0 < mask_ni; ++i0) {
      for (int j0 = 0; j0 < mask_nj; ++j0) {
      i = i0 + mask_i0;
      j = j0 + mask_j0;
      fprintf(fp,"%d %d %d\n",i,j,mask[i0*mask_nj+j0]);
      }
      }
      fclose(fp);
      }
      cout << "take a look..." << endl;
      getchar();
    */

  }
        
  void itpSpatialStats(Param * param){

    COUT1("SPATIAL_STATS...");

    if ((!pngCache)||(pngCache->imageList.size() == 0)) {
      CERR("cannot compute SPATIAL_STATS: pngCache is empty.");
    }

    // a MASK <filename> can be specified where the mask contains
    // a series of image coordinates (i,j) forming a closed polygon
    
    int mask_i0,mask_j0;
    int mask_ni,mask_nj;
    int * mask = NULL;
    
    string outputName = "spatial_stats";

    // look for other param's...
    int iarg = 0;
    while (iarg < param->size()) {
      string token = param->getString(iarg++);
      if (token == "MASK") {
        if (mpi_rank == 0) cout << " > MASK" << endl;
        token = param->getString(iarg++);
        // read the mask file: expect i0 j0\n i1 j1\n ...
        FILE * fp = fopen(token.c_str(),"r");
        if (fp == NULL) {
          CERR("cannot find or open mask file: " << token);
        }
        vector<pair<int,int> > mask_ij; // sequence of image coordinates defining closed polygonal mask
        mask_i0 = TWO_BILLION;
        mask_j0 = TWO_BILLION;
        mask_ni = -TWO_BILLION;
        mask_nj = -TWO_BILLION;
        char line[128];
        int i,j;
        while (fgets(line,128,fp) != NULL) {
          const int n = sscanf(line,"%d %d\n",&i,&j);
          if (n == 2) {
            cout << " > i,j: " << i << " " << j << endl;
            mask_ij.push_back(pair<int,int>(i,j));
            mask_i0 = min(mask_i0,i);
            mask_j0 = min(mask_j0,j);
            mask_ni = max(mask_ni,i);
            mask_nj = max(mask_nj,j);
          }
        }
        fclose(fp);
        if (mask_ij.empty()) {
          CERR("mask file format problem: " << token);
        }
        mask_ni -= mask_i0-1;
        mask_nj -= mask_j0-1;
        cout << " > building pixel mask..." << endl;
        mask = new int[mask_ni*mask_nj];
        makeMask(mask,mask_i0,mask_j0,mask_ni,mask_nj,mask_ij);
      }
      else {
        // assume this is the name...
        outputName = token;
        if (mpi_rank == 0) cout << " > output prefix will be: " << outputName << endl;
      }
    }
    
    //COUT1("Spatial Image Stats: Name, Min, Max, Avg, RMS");
    size_t fwidth = 13;
    map<string,ofstream*> ofileMap;
    for (int iI = 0; iI<pngCache->size(); ++iI){
      PngData * img = pngCache->getImage(iI);

      vector<pixel_type> pixelTypeVec;
      img->getPixelTypes(pixelTypeVec);

      for (int ipt = 0; ipt<pixelTypeVec.size(); ipt++){
        if (ofileMap.find(pixelTypeVec[ipt].name)==ofileMap.end()){ //first time we found this pixel type...setup output
          stringstream myOutputName;
          myOutputName << outputName << "." << pixelTypeVec[ipt].name << ".dat";

          cout << " > writing spatial stats for " << pixelTypeVec[ipt].name << " data to file \"" << outputName << "." << pixelTypeVec[ipt].name << ".dat" << "\"..." << endl;

          MiscUtils::mkdir_for_file(myOutputName.str());
          ofileMap[pixelTypeVec[ipt].name] = new ofstream();
          ofstream& ofmap = *(ofileMap[pixelTypeVec[ipt].name]);
          ofmap.open(myOutputName.str().c_str());
          if (!ofmap.is_open()) {
            CERR("could not open file: " << myOutputName.str());
          }

          //  the parameter...
          ofmap << "# " << param->str() << "\n";
          ofmap << "# " << pixelTypeVec[ipt].name << " data pixels only" << "\n";

          // the variable name...
          ofmap << "# " <<
            setw(fwidth-2) << "1:image" <<
            setw(fwidth) << "2:area" <<
            setw(fwidth) << "3:centroid-x" <<
            setw(fwidth) << "4:centroid-y" <<
            setw(fwidth) << "5:centroid-z" <<
            setw(fwidth) << "6:avg" <<
            setw(fwidth) << "7:rms" <<
            setw(fwidth) << "8:min" <<
            setw(fwidth) << "9:max" ;
          ofmap << endl;

        }

        double stats[8]; //x,y,z,min,max,avg,rms,area
        int count;
        if (mask != NULL) count = img->computeStatsMask(stats,pixelTypeVec[ipt].flag,mask_i0,mask_j0,mask_ni,mask_nj,mask);
        else  count = img->computeStats(stats,pixelTypeVec[ipt].flag);
        
        COUT1("image " << iI << " of " << pngCache->size() << ": pixel count: " << count << " area: " << stats[7] << " centroid: " << COUT_VEC(stats) << " avg: " << stats[5] << " rms: " << stats[6]);
        
        img->lock = false;
        ofstream& ofmap = *(ofileMap[pixelTypeVec[ipt].name]);
        ofmap << "  " <<
          setw(fwidth-2) << iI <<
          setw(fwidth) << stats[7] <<
          setw(fwidth) << stats[0] <<
          setw(fwidth) << stats[1] <<
          setw(fwidth) << stats[2] <<
          setw(fwidth) << stats[5] << // note order switch to match CONDITIONAL_PROBE output
          setw(fwidth) << stats[6] <<
          setw(fwidth) << stats[3] <<
          setw(fwidth) << stats[4] ;
        ofmap << endl;
      }
    }//(iI < pngCache->size())
    
    if (mask != NULL) delete[] mask;
    
    //close files
    for (map<string,ofstream*>::iterator it=ofileMap.begin(); it!=ofileMap.end(); ++it)
      { it->second->close(); delete it->second; }
    
  }//itpSpatialStats()
  
  void itpSpatialStatsDisk(Param * param){
    
    COUT1("SPATIAL_STATS_DISK...");
    
    if ((!pngCache)||(pngCache->imageList.size() == 0)) {
      CERR("cannot compute SPATIAL_STATS_DISK: no IMAGES");
    }
    
    string outputPrefix = "spatial_stats_disk";
    
    bool b_xc = false;
    double xc[3];
    
    bool b_r = false;
    double r;

    bool b_check = false;
    
    // look for other param's...
    int iarg = 0;
    while (iarg < param->size()) {
      string token = param->getString(iarg++);
      if ((token == "X")||(token == "XP")) {
	FOR_I3 xc[i] = param->getDouble(iarg++);
	b_xc = true;
      }
      else if ((token == "R")||(token == "RP")) {
	r = param->getDouble(iarg++);
	b_r = true;
      }
      else if (token == "CHECK") {
	b_check = true;
      }
      else if (token == "NAME") {
	outputPrefix = param->getString(iarg++);
      }
      else {
        // assume this is the name...
        outputPrefix = token;
      }
    }
    
    cout << " > output prefix will be: " << outputPrefix << endl;
    
    if (!b_xc) {
      CERR("cannot compute SPATIAL_STATS_DISK: missing X <x> <y> <z>");
    }
    if (!b_r) {
      CERR("cannot compute SPATIAL_STATS_DISK: missing R <r>");
    }
    
    MiscUtils::mkdir_for_file(outputPrefix+"/test.png");

    size_t fwidth = 14;
    ofstream ofp;
    ofp.open((outputPrefix+".dat").c_str(),std::ofstream::out|std::ofstream::trunc);
    if (!ofp.is_open()) {
      CERR("could not open file: " << outputPrefix << ".dat");
    }
    
    //  the parameter...
    ofp << "# " << param->str() << "\n";
    // variable names...
    ofp << "# " <<
      setw(fwidth-2) << "1:image-index" <<
      setw(fwidth) << "2:avg" <<
      setw(fwidth) << "3:rms" << 
      setw(fwidth) << "4:min" <<
      setw(fwidth) << "5:max" << endl;
    
    for (int iI = 0; iI<pngCache->size(); ++iI){
      PngData * img = pngCache->getImage(iI);
      cout << " > image " << iI << " of " << pngCache->size() << ": "; 
      cout.flush();
      const int nx = img->getNx();
      const int ny = img->getNy();
      int count = 0;
      double sum = 0.0;
      double sum2 = 0.0;
      double vmin = HUGE_VAL;
      double vmax = -HUGE_VAL;
      double xP[3];
      for (int ipx = 0; ipx < nx*ny; ipx++) {
        if (img->pixel_flag[ipx] == PngData::VOLUME_DATA_PIXEL) {
          const int iP = ipx%nx;
          const int jP = int(ipx/nx);
	  img->convertPxtoXp(iP,jP,xP);
	  if (DIST2(xc,xP) <= r*r) {
	    ++count;
            const double var = img->pixel_data[ipx];
	    sum += var;
	    sum2 += var*var;
	    vmin = min(vmin,var);
	    vmax = max(vmax,var);
	    if (b_check) img->pixel_flag[ipx] = PngData::FLAG_RED;
	  }
	}
      }
      if (count > 0) {
	const double avg = sum/double(count);
	const double rms = sqrt(sum2/double(count)-avg*avg);
	cout << count << " pixels in disk, avg,rms,min,max: " << avg << " " << rms << " " << vmin << " " << vmax << endl;
	ofp << "  " <<
	  setw(fwidth-2) << iI <<
	  setw(fwidth) << avg << 
	  setw(fwidth) << rms << 
	  setw(fwidth) << vmin << 
	  setw(fwidth) << vmax << endl;
      }
      else {
	cout << "WARNING: no volume data pixels in disk" << endl;
      }
      if (b_check) {
	img->initializeWrite(); // necessary?
	char filename[128];
	sprintf(filename,"%s/check_disk.%06d.png",outputPrefix.c_str(),iI);
	img->write(filename);
      }
      img->lock = false;
    }

    ofp.close();
    
  }

  void itpRadialProfile(Param * param) {

    // the radial profile is useful for azimuthal averaging image data: for
    // example, to build the outlet temperature profile for axial combustor 
    // simulations...

    COUT1("RADIAL_PROFILE...");
    
    if ((!pngCache)||(pngCache->imageList.size() == 0)) {
      CERR("cannot compute RADIAL_PROFILE: no IMAGES");
    }
    
    string outputPrefix = "radial_profile";
    
    bool b_wgt = false; // if true, the second half of the IMAGES is assumed to contain a wgt (e.g. rho*ux) 
    double xc[3] = { 0.0, 0.0, 0.0 };
    int nbin = 50;

    // look for other param's...
    int iarg = 0;
    while (iarg < param->size()) {
      string token = param->getString(iarg++);
      if (token == "WGT") {
	b_wgt = true;
      }
      else if (token == "NBIN") {
	nbin = param->getInt(iarg++);
      }
      else if (token == "NAME") {
	outputPrefix = param->getString(iarg++);
      }
      else {
        // assume this is the name...
        outputPrefix = token;
      }
    }
    
    cout << " > output prefix will be: " << outputPrefix << endl;
    
    // decide what "r" is and build the bins by inspecting the first image...
    PngData * img = pngCache->getImage(0);
    const int nx = img->getNx();
    const int ny = img->getNy();
    double xP[3];
    double bbmin[3] = { HUGE_VAL, HUGE_VAL, HUGE_VAL };
    double bbmax[3] = { -HUGE_VAL, -HUGE_VAL, -HUGE_VAL };
    double rmin[3] = { HUGE_VAL, HUGE_VAL, HUGE_VAL };
    double rmax[3] = { -HUGE_VAL, -HUGE_VAL, -HUGE_VAL };
    for (int ipx = 0; ipx < nx*ny; ipx++) {
      if (img->pixel_flag[ipx] == PngData::VOLUME_DATA_PIXEL) {
	const int iP = ipx%nx;
	const int jP = int(ipx/nx);
	img->convertPxtoXp(iP,jP,xP);
	FOR_I3 bbmin[i] = min(bbmin[i],xP[i]);
	FOR_I3 bbmax[i] = max(bbmax[i],xP[i]);
	const double dx[3] = DIFF(xP,xc);
	FOR_I3 {
	  const double r = sqrt(dx[(i+1)%3]*dx[(i+1)%3]+dx[(i+2)%3]*dx[(i+2)%3]);
	  rmin[i] = min(rmin[i],r);
	  rmax[i] = max(rmax[i],r);
	}
      }
    }

    // figure out which plane the images are in. This is used to calculate 
    // the 2D radius for the radial bins...
    
    int id = -1;
    if (bbmax[0]-bbmin[0] < min(bbmax[1]-bbmin[1],bbmax[2]-bbmin[2])) {
      id = 0;
    }
    else if (bbmax[1]-bbmin[1] < bbmax[2]-bbmin[2]) {
      id = 1;
    }
    else {
      id = 2;
    }
    
    // expand the range slightly...

    double dr = rmax[id]-rmin[id];
    rmax[id] += 1.0E-4*dr;
    rmin[id] -= 1.0E-4*dr;
    rmin[id] = max(0.0,rmin[id]);

    cout << " > images in radial plane with id: " << id << " radial range: " << rmin[id] << " to " << rmax[id] << endl; 
    
    double * sum_wgt = new double[nbin];
    double * sum_wgt_var = new double[nbin];
    double * sum_wgt_var2 = new double[nbin];
    double * var_min = new double[nbin];
    double * var_max = new double[nbin];
    for (int ibin = 0; ibin < nbin; ++ibin) {
      sum_wgt[ibin] = sum_wgt_var[ibin] = sum_wgt_var2[ibin] = 0.0;
      var_min[ibin] = HUGE_VAL;
      var_max[ibin] = -HUGE_VAL;
    }
    
    if (b_wgt) {
      
      assert(pngCache->size()%2 == 0);
      
      double s_cum = 0.0;
      double sv_cum = 0.0;
      for (int iI = 0; iI<pngCache->size()/2; ++iI){
	PngData * img = pngCache->getImage(iI);
	PngData * wimg = pngCache->getImage(iI+pngCache->size()/2);
	cout << " > working on image pair " << iI << " of " << pngCache->size()/2 << ", time: " << img->getTime() << " var: " << pngCache->imageList[iI] << ", wgt: " << pngCache->imageList[iI+pngCache->size()/2] << ", "; 
	cout.flush();
	const int nx = img->getNx();
	const int ny = img->getNy();
	int count = 0;
	double s = 0.0;
	double sv = 0.0;
	for (int ipx = 0; ipx < nx*ny; ipx++) {
	  if (img->pixel_flag[ipx] == PngData::VOLUME_DATA_PIXEL) {
	    assert(wimg->pixel_flag[ipx] == PngData::VOLUME_DATA_PIXEL);
	    const int iP = ipx%nx;
	    const int jP = int(ipx/nx);
	    img->convertPxtoXp(iP,jP,xP);
	    const double dx[3] = DIFF(xP,xc);
	    const double r = sqrt(dx[(id+1)%3]*dx[(id+1)%3]+dx[(id+2)%3]*dx[(id+2)%3]);
	    if ((r < rmin[id])||(r > rmax[id])) {
	      cout << "\nWARNING: this images has r outside range of first image: r = " << r << ", skipping. (Need to implement dynamic binning for variable geometry data)" << endl;
	      continue;
	    }
	    const int ibin = int((r-rmin[id])/(rmax[id]-rmin[id])*double(nbin));
	    assert((ibin >= 0)&&(ibin < nbin));
	    const double wgt = wimg->pixel_data[ipx];
	    sum_wgt[ibin] += wgt;
	    //assert(wgt > 0.0);
	    const double var = img->pixel_data[ipx];
	    sum_wgt_var[ibin] += wgt*var;
	    sum_wgt_var2[ibin] += wgt*var*var;
	    var_min[ibin] = min(var_min[ibin],var);
	    var_max[ibin] = max(var_max[ibin],var);
	    ++count;
	    s += wgt;
	    sv += wgt*var;
	  }
	}
	s_cum += s;
	sv_cum += sv;
	cout << "VOLUME_DATA pixels: " << count << " avg: " << sv/s << " cumulative avg: " << sv_cum/s_cum << endl;
	img->lock = false;
	wimg->lock = false;
      }
      
    }
    else {

      double s_cum = 0.0;
      double sv_cum = 0.0;
      for (int iI = 0; iI<pngCache->size(); ++iI){
	PngData * img = pngCache->getImage(iI);
	cout << " > working on image " << iI << " of " << pngCache->size() << ", time: " << img->getTime() << " var: " << pngCache->imageList[iI] << ", "; 
	cout.flush();
	const int nx = img->getNx();
	const int ny = img->getNy();
	int count = 0;
	double s = 0.0;
	double sv = 0.0;
	for (int ipx = 0; ipx < nx*ny; ipx++) {
	  if (img->pixel_flag[ipx] == PngData::VOLUME_DATA_PIXEL) {
	    const int iP = ipx%nx;
	    const int jP = int(ipx/nx);
	    img->convertPxtoXp(iP,jP,xP);
	    const double dx[3] = DIFF(xP,xc);
	    const double r = sqrt(dx[(id+1)%3]*dx[(id+1)%3]+dx[(id+2)%3]*dx[(id+2)%3]);
	    if ((r < rmin[id])||(r > rmax[id])) {
	      cout << "\nWARNING: this images has r outside range of first image: r = " << r << ", skipping. (Need to implement dynamic binning for variable geometry data)" << endl;
	      continue;
	    }
	    const int ibin = int((r-rmin[id])/(rmax[id]-rmin[id])*double(nbin));
	    assert((ibin >= 0)&&(ibin < nbin));
	    sum_wgt[ibin] += 1.0;
	    const double var = img->pixel_data[ipx];
	    sum_wgt_var[ibin] += var;
	    sum_wgt_var2[ibin] += var*var;
	    var_min[ibin] = min(var_min[ibin],var);
	    var_max[ibin] = max(var_max[ibin],var);
	    ++count;
	    s += 1.0;
	    sv += var;
	  }
	}
	s_cum += s;
	sv_cum += sv;
	cout << "VOLUME_DATA pixels: " << count << " avg: " << sv/s << " cumulative avg: " << sv_cum/s_cum << endl;
	img->lock = false;
      }

    }

    // write profile...

    cout << " > writing profile to " << outputPrefix << ".dat..." << endl;
    
    MiscUtils::mkdir_for_file(outputPrefix+".dat");
    size_t fwidth = 14;
    ofstream ofp;
    ofp.open((outputPrefix+".dat").c_str(),std::ofstream::out|std::ofstream::trunc);
    if (!ofp.is_open()) {
      CERR("could not open file: " << outputPrefix << ".dat");
    }
    
    //  the parameter...
    ofp << "# " << param->str() << "\n";
    
    // variable names...
    ofp << "# " <<
      setw(fwidth-2) << "1:bin" <<
      setw(fwidth) << "2:wgt" <<
      setw(fwidth) << "3:r" <<
      setw(fwidth) << "4:frac" <<
      setw(fwidth) << "5:avg" << 
      setw(fwidth) << "6:rms" << 
      setw(fwidth) << "7:min" << 
      setw(fwidth) << "8:max" << endl;
    
    for (int ibin = 0; ibin < nbin; ++ibin) {
      if (sum_wgt[ibin] > 0.0) {
	const double frac = (double(ibin)+0.5)/double(nbin);
	const double r = rmin[id] + frac*(rmax[id]-rmin[id]);
	const double avg = sum_wgt_var[ibin]/sum_wgt[ibin];
	const double rms = sqrt(max(0.0,sum_wgt_var2[ibin]/sum_wgt[ibin]-avg*avg));
	ofp << "  " <<
	  setw(fwidth-2) << ibin <<
	  setw(fwidth) << sum_wgt[ibin] << 
	  setw(fwidth) << r << 
	  setw(fwidth) << frac << 
	  setw(fwidth) << avg << 
	  setw(fwidth) << rms << 
	  setw(fwidth) << var_min[ibin] << 
	  setw(fwidth) << var_max[ibin] << endl;
      }
    }
    
    ofp.close();
    
    delete[] sum_wgt;
    delete[] sum_wgt_var;
    delete[] sum_wgt_var2;
    delete[] var_min;
    delete[] var_max;
    
  }

  void itpAnnularIntegrationProbe(Param * param) {

    // the radial profile is useful for azimuthal averaging image data: for
    // example, to build the outlet temperature profile for axial combustor 
    // simulations...
    
    COUT1("ANNULAR_INTEGRATION_PROBE...");
    
    if ((!pngCache)||(pngCache->imageList.size() == 0)) {
      CERR("please use param IMAGES <image.png> [...] to set the images");
    }
    
    string outputPrefix = "annular_integration_probe";
    
    double xc[3] = { 0.0, 0.0, 0.0 };
    
    double rmin,rmax;
    bool b_rmin = false;
    bool b_rmax = false;
    bool b_check = false;
    int clocking = 0;
    
    int nbin = -1;
    
    // look for other param's...
    int iarg = 0;
    while (iarg < param->size()) {
      string token = param->getString(iarg++);
      if (token == "NBIN") {
	nbin = param->getInt(iarg++);
      }
      else if (token == "RMIN") {
	rmin = param->getDouble(iarg++);
	b_rmin = true;
      }
      else if (token == "RMAX") {
	rmax = param->getDouble(iarg++);
	b_rmax = true;
      }
      else if (token == "CLOCKING") {
	clocking = param->getInt(iarg++);
      }
      else if (token == "CHECK") {
	b_check = true;
      }
      else if (token == "NAME") {
	outputPrefix = param->getString(iarg++);
      }
      else {
        // assume this is the name...
        outputPrefix = token;
      }
    }

    assert(nbin > 0);
    assert((clocking == 0)||(clocking == 1));
    
    cout << " > output prefix will be: " << outputPrefix << endl;

    // step 1: try to get the image orientation from metadata...
    // for now, insist it is axis aligned...
    PngData * img = pngCache->getImage(0);
    ImageMetadata * metadata = img->getMetadata();
    assert(metadata);
    int i0 = -1;
    int i1 = -1;
    if (fabs(metadata->transformMat[0+0]) > max(fabs(metadata->transformMat[0+1]),fabs(metadata->transformMat[0+2]))) {
      i0 = 0;
    }
    else if (fabs(metadata->transformMat[0+1]) > fabs(metadata->transformMat[0+2])) {
      i0 = 1;
    }
    else { 
      i0 = 2;
    }
    if (fabs(metadata->transformMat[4+0]) > max(fabs(metadata->transformMat[4+1]),fabs(metadata->transformMat[4+2]))) {
      i1 = 0;
    }
    else if (fabs(metadata->transformMat[4+1]) > fabs(metadata->transformMat[4+2])) {
      i1 = 1;
    }
    else { 
      i1 = 2;
    }
    assert(i0 >= 0);
    assert(i1 >= 0);
    cout << " > images coordinate directions appear to be: " << i0 << " " << i1 << endl;
    
    // step 2: has the user requested CHECK? if so, write out an images with the probe areas colored..
    MiscUtils::mkdir_for_file(outputPrefix+".dat"); // also serves to mkdir for the check image if requested
    if (b_check) {
      const int nx = img->getNx();
      const int ny = img->getNy();
      double xp[3];
      for (int ipx = 0; ipx < nx*ny; ipx++) {
        if (img->pixel_flag[ipx] == PngData::VOLUME_DATA_PIXEL) {
          const int ip = ipx%nx;
          const int jp = int(ipx/nx);
	  img->convertPxtoXp(ip,jp,xp);
	  double r = sqrt((xp[i0]-xc[i0])*(xp[i0]-xc[i0]) + (xp[i1]-xc[i1])*(xp[i1]-xc[i1]));
	  if (b_rmin && (r < rmin))
	    continue;
	  if (b_rmax && (r > rmax))
	    continue;
	  double angle = atan2(xp[i1]-xc[i1],xp[i0]-xc[i0])*180.0/M_PI;
	  if (clocking) angle -= 180.0/double(nbin); // clocking adjustment
	  if (angle < 0.0) angle += 360.0;
	  assert(angle >= 0.0);
	  int ibin = angle/(360.0/double(nbin));
	  assert((ibin >= 0)&&(ibin < nbin));
	  if (ibin%2 == 0) {
	    img->pixel_flag[ipx] = PngData::FLAG_GREEN;
	  }
	  else {
	    img->pixel_flag[ipx] = PngData::FLAG_RED;
	  }
	}
	else {
	  // otherwise, make sure there is no other red or green in the image
	  // that would confuse the meaning of the check image... 
	  assert(img->pixel_flag[ipx] != PngData::FLAG_GREEN);
	  assert(img->pixel_flag[ipx] != PngData::FLAG_RED);
	}
      }
      img->initializeWrite(); // necessary?
      char filename[128];
      sprintf(filename,"%s.png",outputPrefix.c_str());
      img->write(filename);
      cout << "\n\n***********************************************************************\n" << 
	"check image written to " << filename << "\nRemove CHECK token to process\n" << 
	"***********************************************************************\n\n" << endl;
      throw(0);
    }
    
    MiscUtils::mkdir_for_file(outputPrefix+".dat");
    size_t fwidth = 16;
    ofstream ofp;
    ofp.open((outputPrefix+".dat").c_str(),std::ofstream::out|std::ofstream::trunc);
    if (!ofp.is_open()) {
      CERR("could not open file: " << outputPrefix << ".dat");
    }
    
    //  the parameter...
    ofp << "# " << param->str() << "\n";
    
    // column names...
    ofp << "# 1:step 2:time";
    for (int ibin = 0; ibin < nbin; ++ibin)
      ofp << " " << ibin+3 << ":probe" << ibin;
    ofp << endl;

    double *avg = new double[nbin];
    double *count = new double[nbin];
    for (int iI = 0; iI<pngCache->size(); ++iI){
      for (int ibin = 0; ibin < nbin; ++ibin) {
	avg[ibin] = 0.0;
	count[ibin] = 0.0;
      }
      img = pngCache->getImage(iI);
      int step = MiscUtils::getIndexFromFilename(pngCache->imageList[iI]);
      cout << " > working on image " << iI << " of " << pngCache->size() << ", index: " << step << ", time: " << img->getTime() << ", file: " << pngCache->imageList[iI] << endl;
      const int nx = img->getNx();
      const int ny = img->getNy();
      double xp[3];
      for (int ipx = 0; ipx < nx*ny; ipx++) {
	if (img->pixel_flag[ipx] == PngData::VOLUME_DATA_PIXEL) {
          const int ip = ipx%nx;
          const int jp = int(ipx/nx);
	  img->convertPxtoXp(ip,jp,xp);
	  double r = sqrt((xp[i0]-xc[i0])*(xp[i0]-xc[i0]) + (xp[i1]-xc[i1])*(xp[i1]-xc[i1]));
	  if (b_rmin && (r < rmin))
	    continue;
	  if (b_rmax && (r > rmax))
	    continue;
	  double angle = atan2(xp[i1]-xc[i1],xp[i0]-xc[i0])*180.0/M_PI;
	  if (clocking) angle -= 180.0/double(nbin); // clocking adjustment
	  if (angle < 0.0) angle += 360.0;
	  assert(angle >= 0.0);
	  int ibin = angle/(360.0/double(nbin));
	  assert((ibin >= 0)&&(ibin < nbin));
	  avg[ibin] += img->pixel_data[ipx];
	  count[ibin] += 1.0;
	}
      }
      img->lock = false;
      ofp << step << " " << setprecision(8) << img->getTime();
      for (int ibin = 0; ibin < nbin; ++ibin) {
	ofp << " " << avg[ibin]/count[ibin];
      }
      ofp << endl;
    }
    ofp.close();
    delete[] avg;
    delete[] count;

    cout << " > probe data written to file " << outputPrefix+".dat" << endl;
    
  }
  
  void itpTumble(Param * param){

    COUT1("TUMBLE...");

    if ((!pngCache)||(pngCache->imageList.size() == 0)) {
      CERR("cannot compute TUMBLE: pngCache is empty.");
    }

    // also, we need an even number of images...
    if (pngCache->size()%2 != 0) {
      CERR("cannot compute TUMBLE: expecting even image count.");
    }

    const int nImgPair = pngCache->size()/2;

    // a MASK <filename> can be specified where the mask contains
    // a series of image coordinates (i,j) forming a closed polygon
    
    int mask_i0,mask_j0;
    int mask_ni,mask_nj;
    int * mask = NULL;
    
    string outputPrefix = "tumble";
    
    bool b_omega = false;
    double omega;

    bool b_range = false;
    double range[2];

    // look for other param's...
    int iarg = 0;
    while (iarg < param->size()) {
      string token = param->getString(iarg++);
      if (token == "MASK") {
        if (mpi_rank == 0) cout << " > MASK" << endl;
        token = param->getString(iarg++);
        // read the mask file: expect i0 j0\n i1 j1\n ...
        FILE * fp = fopen(token.c_str(),"r");
        if (fp == NULL) {
          CERR("cannot find or open mask file: " << token);
        }
        vector<pair<int,int> > mask_ij; // sequence of image coordinates defining closed polygonal mask
        mask_i0 = TWO_BILLION;
        mask_j0 = TWO_BILLION;
        mask_ni = -TWO_BILLION;
        mask_nj = -TWO_BILLION;
        char line[128];
        int i,j;
        while (fgets(line,128,fp) != NULL) {
          const int n = sscanf(line,"%d %d\n",&i,&j);
          if (n == 2) {
            cout << " > i,j: " << i << " " << j << endl;
            mask_ij.push_back(pair<int,int>(i,j));
            mask_i0 = min(mask_i0,i);
            mask_j0 = min(mask_j0,j);
            mask_ni = max(mask_ni,i);
            mask_nj = max(mask_nj,j);
          }
        }
        fclose(fp);
        if (mask_ij.empty()) {
          CERR("mask file format problem: " << token);
        }
        mask_ni -= mask_i0-1;
        mask_nj -= mask_j0-1;
        cout << " > building pixel mask..." << endl;
        mask = new int[mask_ni*mask_nj];
        makeMask(mask,mask_i0,mask_j0,mask_ni,mask_nj,mask_ij);
      }
      else if (token == "RPT") {
        const double rpt = param->getDouble(iarg++);
        omega = 2.0*M_PI*rpt;
        if (mpi_rank == 0) cout << " > RPT " << rpt << " is OMEGA " << omega << endl;
        b_omega = true;
      }
      else if (token == "RPM") {
        const double rpm = param->getDouble(iarg++);
        omega = 2.0*M_PI*rpm/60.0;
        if (mpi_rank == 0) cout << " > RPM " << rpm << " is OMEGA " << omega << endl;
        b_omega = true;
      }
      else if (token == "OMEGA") {
        omega = param->getDouble(iarg++);
        if (mpi_rank == 0) cout << " > OMEGA " << omega << endl;
        b_omega = true;
      }
      else if (token == "RANGE") {
        range[0] = param->getDouble(iarg++);
        range[1] = param->getDouble(iarg++);
        if (mpi_rank == 0) cout << " > RANGE " << range[0] << " " << range[1] << endl;
        b_range = true;
      }
      else {
        // assume this is the name...
        outputPrefix = token;
        if (mpi_rank == 0) cout << " > output file will be: " << outputPrefix << ".dat" << endl;
      }
    }

    if (mask == NULL) {
      CERR("expecting MASK for tumble calc");
    }
    
    if (!b_omega) {
      CERR("omega was not specified. It is used to normalize the tumble. Add RPT <rpt>, or RPM <rpm>, or OMEGA <omega>");
    }
    
    MiscUtils::mkdir_for_file(outputPrefix+"/test.png");
    
    size_t fwidth = 14;
    ofstream ofp;
    ofp.open((outputPrefix+".dat").c_str(), std::ofstream::out | std::ofstream::trunc);
    if (!ofp.is_open()) {
      CERR("could not open file: " << outputPrefix << ".dat");
    }
    
    //  the parameter...
    ofp << "# " << param->str() << "\n";
    // also record omega...
    ofp << "# omega = " << omega << "\n";
    // variable names...
    ofp << "# " <<
      setw(fwidth-2) << "1:image-pair" <<
      setw(fwidth) << "2:centroid-x" <<
      setw(fwidth) << "3:centroid-y" <<
      setw(fwidth) << "4:centroid-z" <<
      setw(fwidth) << "5:tumble" <<
      setw(fwidth) << "6:moi" << endl;
    
    // loop on image pairs...
    for (int iP = 0; iP < nImgPair; ++iP) {
      
      PngData * uimg = pngCache->getImage(iP);
      PngData * vimg = pngCache->getImage(iP+nImgPair);

      const int nx = uimg->getNx();
      const int ny = uimg->getNy();
      if ((nx != vimg->getNx())||(ny != vimg->getNy())) {
        CERR("image pair differs in size. Cannot compute tumble.");
      }

      // pixels masked with mask[i0*mask_nj+j0] == 1 are included in stats...
      
      double iAvg = 0.0;
      double jAvg = 0.0;
      double xc[3] = { 0.0, 0.0, 0.0 };
      int count = 0;
      double xP[3];
      for (int ipx =0; ipx < nx*ny; ipx++) {
        if (uimg->pixel_flag[ipx] != vimg->pixel_flag[ipx]) {
          CERR("image pair differs in pixel_flag. Cannot compute tumble.");
        }
        if (uimg->pixel_flag[ipx] == PngData::VOLUME_DATA_PIXEL) {
          const int iP = ipx%nx;
          const int jP = (int)ipx/nx;
          if ((iP >= mask_i0)&&(iP < mask_i0+mask_ni)&&(jP >= mask_j0)&&(jP < mask_j0+mask_nj)&&(mask[(iP-mask_i0)*mask_nj+(jP-mask_j0)] == 1)) {
            iAvg += double(iP);
            jAvg += double(jP);
            uimg->convertPxtoXp(iP,jP,xP);
            FOR_I3 xc[i] += xP[i];
            ++count;
          }
          else {
            // VOLUME pixels outside of the mask are flagged GRAY here...
            uimg->pixel_flag[ipx] = PngData::FLAG_GRAY;
          }
        }
      }
      
      double tumble = 0.0;
      double tumble_range[2] = { HUGE_VAL, -HUGE_VAL };
      double moi = 0.0;
      if (count > 0) {
        
        //image geom centroid in pixel coordinates
        iAvg /= double(count);
        jAvg /= double(count);
        FOR_I3 xc[i] /= double(count);
        
        // recall ip,jp,depth to xp conversion...
        //xp[i] = ip*metadata->transformMat[0+i] +
        //        jp*metadata->transformMat[4+i] +
        //        dpTh_us*metadata->transformMat[8+i] +
        
        ImageMetadata * metadata = uimg->getMetadata();
        assert(metadata);
        
        double ipixel_length = 0.0;
        if (fabs(metadata->transformMat[0+0]) > max(fabs(metadata->transformMat[0+1]),fabs(metadata->transformMat[0+2]))) {
          if ((metadata->transformMat[0+0] < 0.0)&&
              (fabs(metadata->transformMat[0+1]) < -1.0E-6*metadata->transformMat[0+0])&&
              (fabs(metadata->transformMat[0+2]) < -1.0E-6*metadata->transformMat[0+0])) {
            // i == -x...
            ipixel_length = metadata->transformMat[0+0];
          }
        }

        // get the j-pixel length from the transform...
        double jpixel_length = 0.0;
        if (fabs(metadata->transformMat[4+2]) > max(fabs(metadata->transformMat[4+0]),fabs(metadata->transformMat[4+1]))) {
          if ((metadata->transformMat[4+2] < 0.0)&&
              (fabs(metadata->transformMat[4+0]) < -1.0E-6*metadata->transformMat[4+2])&&
              (fabs(metadata->transformMat[4+1]) < -1.0E-6*metadata->transformMat[4+2])) {
            // j == -z...
            jpixel_length = metadata->transformMat[4+2];
          }
        }
        
        // if pixel_length was not set above, you need to add a case to the above...
        if ((ipixel_length == 0.0)||(jpixel_length == 0.0)) {
          FOR_I3 {
            cout << "i = " << i << endl;
            cout << " > metadata->transformMat[0+i]: " << metadata->transformMat[0+i] << endl;
            cout << " > metadata->transformMat[4+i]: " << metadata->transformMat[4+i] << endl;
            cout << " > metadata->transformMat[8+i]: " << metadata->transformMat[8+i] << endl;
          }
          CERR("this image orientation not implemented for tumble");
        }
        
        for (int ipx =0; ipx < nx*ny; ipx++) {
          if (uimg->pixel_flag[ipx] == PngData::VOLUME_DATA_PIXEL) {
            const int iP = ipx%nx;
            const int jP = (int)ipx/nx;
            // no need to check the mask again. We set everyone outside the mask to red above...
            const double u = uimg->pixel_data[ipx];
            const double v = vimg->pixel_data[ipx];
            const double dx = (double(iP)-iAvg)*ipixel_length;
            const double dy = (double(jP)-jAvg)*jpixel_length;
            const double this_tumble = u*dy - v*dx;
            tumble_range[0] = min(tumble_range[0],this_tumble);
            tumble_range[1] = max(tumble_range[1],this_tumble);
            uimg->pixel_data[ipx] = this_tumble;
            tumble += this_tumble;
            moi += dx*dx + dy*dy;
          }
        }

        // we are done with vimg...
        vimg->lock = false;
        
        // and normalize and write the tumble from uimg...
        for (int ipx =0; ipx < nx*ny; ipx++) {
          if (uimg->pixel_flag[ipx] == PngData::VOLUME_DATA_PIXEL) {
            uimg->pixel_data[ipx] *= double(count)/(moi*omega);
          }
        }

        // for reporting below...
        FOR_I2 tumble_range[i] *= double(count)/(moi*omega);

        // no need to touch the metadata: set the new variable name like this...
        uimg->varId = "tumble";
        if (b_range) {
          uimg->setRange(range);
        }
        else {
          uimg->rescaleRangesToData();
        }
        uimg->initializeWrite();

        char filename[128];
        sprintf(filename,"%s/tumble.%06d.png",outputPrefix.c_str(),iP);
        
        uimg->write(filename);
        uimg->lock = false;
        
        cout << " > image pair: " << iP << " of " << nImgPair << ": pixel count: " << count << " centroid: " << COUT_VEC(xc) << 
          " tumble: " << tumble/(moi*omega) << " tumble img: " << filename << " range: " << tumble_range[0] << " " << tumble_range[1] << endl;
        
      }
      else {
        
        cout << " > image pair: " << iP << " of " << nImgPair << ": no VOLUME_DATA found for tumble calc" << endl;
        
      }
      
      ofp << "  " <<
        setw(fwidth-2) << iP <<
        setw(fwidth) << xc[0] <<
        setw(fwidth) << xc[1] <<
        setw(fwidth) << xc[2] <<
        setw(fwidth) << tumble/(moi*omega) << 
        setw(fwidth) << moi << endl;
      
    }//(iP < nImgPair)
    
    ofp.close();

    if (mask != NULL) delete[] mask;

    cout << " > done TUMBLE" << endl;
    
  }//itpTumble()
  
  void itpWriteVelocity(Param * param){

    COUT1("WRITE_VECTOR...");

    if ((!pngCache)||(pngCache->imageList.size() == 0)) {
      CERR("cannot process WRITE_VELOCITY: pngCache is empty.");
    }
    
    // a MASK <filename> can be specified where the mask contains
    // a series of image coordinates (i,j) forming a closed polygon
    
    int mask_i0,mask_j0;
    int mask_ni,mask_nj;
    int * mask = NULL;
    
    string outputPrefix = "velocity";
    int starting_index = 0;
    
    int stride = 1;

    bool b_uw_only = false;

    // look for other param's...
    int iarg = 0;
    while (iarg < param->size()) {
      string token = param->getString(iarg++);
      if (token == "MASK") {
        if (mpi_rank == 0) cout << " > MASK" << endl;
        token = param->getString(iarg++);
        // read the mask file: expect i0 j0\n i1 j1\n ...
        FILE * fp = fopen(token.c_str(),"r");
        if (fp == NULL) {
          CERR("cannot find or open mask file: " << token);
        }
        vector<pair<int,int> > mask_ij; // sequence of image coordinates defining closed polygonal mask
        mask_i0 = TWO_BILLION;
        mask_j0 = TWO_BILLION;
        mask_ni = -TWO_BILLION;
        mask_nj = -TWO_BILLION;
        char line[128];
        int i,j;
        while (fgets(line,128,fp) != NULL) {
          const int n = sscanf(line,"%d %d\n",&i,&j);
          if (n == 2) {
            cout << " > i,j: " << i << " " << j << endl;
            mask_ij.push_back(pair<int,int>(i,j));
            mask_i0 = min(mask_i0,i);
            mask_j0 = min(mask_j0,j);
            mask_ni = max(mask_ni,i);
            mask_nj = max(mask_nj,j);
          }
        }
        fclose(fp);
        if (mask_ij.empty()) {
          CERR("mask file format problem: " << token);
        }
        mask_ni -= mask_i0-1;
        mask_nj -= mask_j0-1;
        cout << " > building pixel mask..." << endl;
        mask = new int[mask_ni*mask_nj];
        makeMask(mask,mask_i0,mask_j0,mask_ni,mask_nj,mask_ij);
      }
      else if (token == "STARTING_INDEX") {
	starting_index = param->getInt(iarg++);
      }
      else if (token == "STRIDE") {
	stride = param->getInt(iarg++);
      }
      else if (token == "UW_ONLY") {
	b_uw_only = true;
      }
      else if (token == "NAME") {
	assert(outputPrefix == "velocity");
	outputPrefix = param->getString(iarg++);
      }
      else {
        // assume this is the name...
	cout << "assuming this is the name: " << token << ". Please explicitly specify the NAME param." << endl;
	assert(outputPrefix == "velocity");
	outputPrefix = token;
      }
    }

    if (mask == NULL) {
      COUT1("no MASK will be applied");
    }

    if (mpi_rank == 0) cout << " > output files will be: " << outputPrefix << "/uvw.xxxxxx.dat, and start with index: " << starting_index << endl;
    
    MiscUtils::mkdir_for_file(outputPrefix+"/test.png");
    
    int nImgSet;
    if (b_uw_only) {
      assert(pngCache->size()%2 == 0);
      nImgSet = pngCache->size()/2;
    }
    else {
      assert(pngCache->size()%3 == 0);
      nImgSet = pngCache->size()/3;
    }
    
    PngData * uimg = NULL;
    PngData * vimg = NULL;
    PngData * wimg = NULL;
    
    int nx,ny;
    for (int ii = 0; ii < nImgSet; ++ii) {

      char filename[128];
      sprintf(filename,"%s/uvw.%06d.dat",outputPrefix.c_str(),ii+starting_index);
      
      ofstream ofp;
      ofp.open(filename, std::ofstream::out | std::ofstream::trunc);
      if (!ofp.is_open()) {
	CERR("could not open file: " << filename);
      }
    
      //  the parameter...
      ofp << "# " << param->str() << "\n";
      ofp << "# image-set: " << ii << " time: " << pngCache->getImage(ii)->getTime() << endl;;

      // variable names...
      ofp << "# " <<
	setw(14) << "1:x" <<
	setw(14) << "2:y" <<
	setw(14) << "3:z" <<
	setw(14) << "4:velocity-x" <<
	setw(14) << "5:velocity-y" <<
	setw(14) << "6:velocity-z" << endl;
    
      // loop on image trios...
      if (b_uw_only) {
	uimg = pngCache->getImage(ii);
	wimg = pngCache->getImage(ii+nImgSet);
      }
      else {
	uimg = pngCache->getImage(ii);
	vimg = pngCache->getImage(ii+nImgSet);
	wimg = pngCache->getImage(ii+nImgSet*2);
      }

      assert(uimg);
      nx = uimg->getNx();
      ny = uimg->getNy();
      
      // check other image sizes...
      if (vimg && ((nx != vimg->getNx())||(ny != vimg->getNy()))) {
        CERR("v image differs in size. Cannot process WRITE_VELOCITY.");
      }
      if (wimg && ((nx != wimg->getNx())||(ny != wimg->getNy()))) {
        CERR("w image differs in size. Cannot process WRITE_VELOCITY.");
      }
      
      int count = 0;
      assert(stride >= 1);
      assert(mask == NULL);
      for (int iP0 = 0; iP0 < nx; iP0 += stride) {
	for (int jP0 = 0; jP0 < ny; jP0 += stride) {
	  double xPsum[3] = { 0.0, 0.0, 0.0 };
	  double uPsum[3] = { 0.0, 0.0, 0.0 };
	  int wgt = 0;
	  for (int diP = 0; diP < stride; ++diP) {
	    const int iP = iP0+diP;
	    if (iP >= nx)
	      continue;
	    for (int djP = 0; djP < stride; ++djP) {
	      const int jP = jP0+djP;
	      if (jP >= ny)
		continue;
	      const int ipx = jP*nx+iP;
	      if (uimg->pixel_flag[ipx] == PngData::VOLUME_DATA_PIXEL) {
		assert((vimg == NULL)||(vimg->pixel_flag[ipx] == PngData::VOLUME_DATA_PIXEL));
		assert((wimg == NULL)||(wimg->pixel_flag[ipx] == PngData::VOLUME_DATA_PIXEL));
		// if mask != NULL, pixels masked with mask[i0*mask_nj+j0] == 1 are included...
		if ((mask == NULL)||((mask != NULL)&&(iP >= mask_i0)&&(iP < mask_i0+mask_ni)&&(jP >= mask_j0)&&(jP < mask_j0+mask_nj)&&(mask[(iP-mask_i0)*mask_nj+(jP-mask_j0)] == 1))) {
		  ++wgt;
		  double xP[3];
		  uimg->convertPxtoXp(iP,jP,xP);
		  FOR_I3 xPsum[i] += xP[i];
		  uPsum[0] += uimg->pixel_data[ipx];
		  if (vimg) uPsum[1] += vimg->pixel_data[ipx];
		  if (wimg) uPsum[2] += wimg->pixel_data[ipx];
		}
	      }
	    }
	  }
	  if (wgt > 0) {
	    ++count;
	    const double inv_wgt = 1.0/double(wgt);
	    ofp << "  " <<
	      setw(14) << xPsum[0]*inv_wgt <<
	      setw(14) << xPsum[1]*inv_wgt <<
	      setw(14) << xPsum[2]*inv_wgt <<
	      setw(14) << uPsum[0]*inv_wgt << 
	      setw(14) << uPsum[1]*inv_wgt <<
	      setw(14) << uPsum[2]*inv_wgt << endl;
	  }
	}
      }

      ofp.close();
      cout << " > vector count: " << count << endl;
      
      // we are done with these images...
      if (uimg) uimg->lock = false;
      if (vimg) vimg->lock = false;
      if (wimg) wimg->lock = false;
      
    }

    if (mask != NULL) delete[] mask;

    cout << " > done WRITE_VELOCITY" << endl;
    
  }//itpWriteVelocity()

  void itpWriteScalar(Param * param){

    COUT1("WRITE_SCALAR...");

    if ((!pngCache)||(pngCache->imageList.size() == 0)) {
      CERR("cannot process WRITE_SCALAR: pngCache is empty.");
    }

    // a MASK <filename> can be specified where the mask contains
    // a series of image coordinates (i,j) forming a closed polygon
    
    int mask_i0,mask_j0;
    int mask_ni,mask_nj;
    int * mask = NULL;
    
    string outputPrefix = "scalar";
    int starting_index = 0;
    
    int stride = 1;

    // look for other param's...
    int iarg = 0;
    while (iarg < param->size()) {
      string token = param->getString(iarg++);
      if (token == "MASK") {
        if (mpi_rank == 0) cout << " > MASK" << endl;
        token = param->getString(iarg++);
        // read the mask file: expect i0 j0\n i1 j1\n ...
        FILE * fp = fopen(token.c_str(),"r");
        if (fp == NULL) {
          CERR("cannot find or open mask file: " << token);
        }
        vector<pair<int,int> > mask_ij; // sequence of image coordinates defining closed polygonal mask
        mask_i0 = TWO_BILLION;
        mask_j0 = TWO_BILLION;
        mask_ni = -TWO_BILLION;
        mask_nj = -TWO_BILLION;
        char line[128];
        int i,j;
        while (fgets(line,128,fp) != NULL) {
          const int n = sscanf(line,"%d %d\n",&i,&j);
          if (n == 2) {
            cout << " > i,j: " << i << " " << j << endl;
            mask_ij.push_back(pair<int,int>(i,j));
            mask_i0 = min(mask_i0,i);
            mask_j0 = min(mask_j0,j);
            mask_ni = max(mask_ni,i);
            mask_nj = max(mask_nj,j);
          }
        }
        fclose(fp);
        if (mask_ij.empty()) {
          CERR("mask file format problem: " << token);
        }
        mask_ni -= mask_i0-1;
        mask_nj -= mask_j0-1;
        cout << " > building pixel mask..." << endl;
        mask = new int[mask_ni*mask_nj];
        makeMask(mask,mask_i0,mask_j0,mask_ni,mask_nj,mask_ij);
      }
      else if (token == "STARTING_INDEX") {
	starting_index = param->getInt(iarg++);
      }
      else if (token == "STRIDE") {
	stride = param->getInt(iarg++);
      }
      else if (token == "NAME") {
	assert(outputPrefix == "scalar");
	outputPrefix = param->getString(iarg++);
      }
      else {
        // assume this is the name...
	cout << "assuming this is the name: " << token << ". Please explicitly specify the NAME param." << endl;
	assert(outputPrefix == "scalar");
	outputPrefix = token;
      }
    }

    if (mask == NULL) {
      COUT1("no MASK will be applied");
    }

    if (mpi_rank == 0) cout << " > output files will be: " << outputPrefix << "/scalar.xxxxxx.dat, and start with index: " << starting_index << endl;
    
    MiscUtils::mkdir_for_file(outputPrefix+"/test.png");
    
    int nImg = pngCache->size();
    
    for (int ii = 0; ii < nImg; ++ii) {

      char filename[128];
      sprintf(filename,"%s/scalar.%06d.dat",outputPrefix.c_str(),ii+starting_index);
      
      ofstream ofp;
      ofp.open(filename, std::ofstream::out | std::ofstream::trunc);
      if (!ofp.is_open()) {
	CERR("could not open file: " << filename);
      }
    
      //  the parameter...
      ofp << "# " << param->str() << "\n";
      ofp << "# image-set: " << ii << " time: " << pngCache->getImage(ii)->getTime() << endl;;

      // variable names...
      ofp << "# " <<
	setw(14) << "1:x" <<
	setw(14) << "2:y" <<
	setw(14) << "3:z" <<
	setw(14) << "4:scalar" << endl;
      
      PngData * img = pngCache->getImage(ii);
      const int nx = img->getNx();
      const int ny = img->getNy();
      
      int count = 0;
      assert(stride >= 1);
      assert(mask == NULL);
      for (int iP0 = 0; iP0 < nx; iP0 += stride) {
	for (int jP0 = 0; jP0 < ny; jP0 += stride) {
	  double xPsum[3] = { 0.0, 0.0, 0.0 };
	  double ssum = 0.0;
	  int wgt = 0;
	  for (int diP = 0; diP < stride; ++diP) {
	    const int iP = iP0+diP;
	    if (iP >= nx)
	      continue;
	    for (int djP = 0; djP < stride; ++djP) {
	      const int jP = jP0+djP;
	      if (jP >= ny)
		continue;
	      const int ipx = jP*nx+iP;
	      if (img->pixel_flag[ipx] == PngData::VOLUME_DATA_PIXEL) {
		// if mask != NULL, pixels masked with mask[i0*mask_nj+j0] == 1 are included...
		if ((mask == NULL)||((mask != NULL)&&(iP >= mask_i0)&&(iP < mask_i0+mask_ni)&&(jP >= mask_j0)&&(jP < mask_j0+mask_nj)&&(mask[(iP-mask_i0)*mask_nj+(jP-mask_j0)] == 1))) {
		  ++wgt;
		  double xP[3];
		  img->convertPxtoXp(iP,jP,xP);
		  FOR_I3 xPsum[i] += xP[i];
		  ssum += img->pixel_data[ipx];
		}
	      }
	    }
	  }
	  if (wgt > 0) {
	    ++count;
	    const double inv_wgt = 1.0/double(wgt);
	    ofp << "  " <<
	      setw(14) << xPsum[0]*inv_wgt <<
	      setw(14) << xPsum[1]*inv_wgt <<
	      setw(14) << xPsum[2]*inv_wgt <<
	      setw(14) << ssum*inv_wgt << endl;
	  }
	}
      }

      ofp.close();
      cout << " > scalar count: " << count << endl;
      
      // we are done with these images...
      assert(img);
      img->lock = false;
      
    }
    
    if (mask != NULL) delete[] mask;
    
    cout << " > done WRITE_SCALAR" << endl;
    
  }//itpWriteScalar()

  void itpVector(Param * param){

    COUT1("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX VECTOR...");

    if ((!pngCache)||(pngCache->imageList.size() == 0)) {
      CERR("cannot compute vectors: pngCache is empty.");
    }

    // also, we need an even number of images...
    if (pngCache->size()%2 != 0) {
      CERR("cannot compute vectors: expecting even image count.");
    }
    
    double spacing = 12.0; // default spacing of vectors in pixels
    const double min_spacing = 5.0; // set an appropriate min
    
    double scale;
    bool b_scale = false;
    
    const int nImgPair = pngCache->size()/2;

    string outputPrefix = "vector";
    
    // look for other param's...
    int iarg = 0;
    while (iarg < param->size()) {
      string token = param->getString(iarg++);
      if (token == "SPACING") {
	spacing = param->getDouble(iarg++);
	if (spacing < min_spacing) {
	  spacing = min_spacing;
	  if (mpi_rank == 0) cout << " > vector SPACING set to " << spacing << " pixels (minimum)" << endl;
	}
	else {
	  if (mpi_rank == 0) cout << " > vector SPACING set to " << spacing << " pixels" << endl;
	}
      }
      else if (token == "SCALE") {
	scale = param->getDouble(iarg++);
	b_scale = true;
	if (mpi_rank == 0) cout << " > vector SCALE set to " << scale << endl;
      }
      else if (token == "PREFIX") {
        // assume this is the name...
        outputPrefix = param->getString(iarg++);
        if (mpi_rank == 0) cout << " > PREFIX will be: " << outputPrefix << endl;
      }
    }
    
    MiscUtils::mkdir_for_file(outputPrefix);

    // we need an array of points that sample the image velocities using 
    // a local kernel. these are stored on a grid that is coarser than the 
    // pixel i,j.
    
    double (*uv)[3] = NULL; // contains u,v,wgt
    int uv_size = 0;
    int sx = -1;
    int sy = -1;
    
    // loop on image pairs...
    for (int iP = 0; iP < nImgPair; ++iP) {
      
      PngData * uimg = pngCache->getImage(iP);
      PngData * vimg = pngCache->getImage(iP+nImgPair);
      
      const int nx = uimg->getNx();
      const int ny = uimg->getNy();
      if ((nx != vimg->getNx())||(ny != vimg->getNy())) {
        CERR("image pair differs in size. Cannot compute vectors.");
      }

      // set sx,sy using the image size and spacing, and ensure the velocity vector
      // grid is large enough for this image pair...
      
      sx = int(double(nx)/spacing);
      sy = int(double(ny)/spacing);
      if (uv_size < sx*sy) {
	if (uv != NULL) delete[] uv;
	uv_size = sx*sy;
	uv = new double[uv_size][3];
      }
      
      // use a Gaussian kernel to 

      const int spacing_int = int(spacing);
      for (int j = 0; j < sy; ++j) {
	for (int i = 0; i < sx; ++i) {
	  const int ij = j*sx+i;
	  FOR_K3 uv[ij][k] = 0.0;
	  // now loop through the image pixels weighting the velocity...
	  const double ip0 = (double(i)+0.5)*spacing;
	  const double jp0 = (double(j)+0.5)*spacing;
	  for (int dip = -spacing_int; dip <= spacing_int; ++dip) {
	    for (int djp = -spacing_int; djp <= spacing_int; ++djp) {
	      double rp2 = double(dip*dip+djp*djp)/(spacing*spacing);
	      if (rp2 < 0.999) {
		const int ip = int(ip0 + double(dip));
		const int jp = int(jp0 + double(djp));
		if ((ip >= 0)&&(ip < nx)&&(jp >= 0)&&(jp < ny)) {
		  const int ipx = jp*nx + ip;
		  if ((uimg->pixel_flag[ipx] == PngData::VOLUME_DATA_PIXEL)&&(vimg->pixel_flag[ipx] == PngData::VOLUME_DATA_PIXEL)) {
		    const double wgt = exp(-1.0/(1.0-rp2));
		    uv[ij][0] += wgt*uimg->pixel_data[ipx];
		    uv[ij][1] += wgt*vimg->pixel_data[ipx];
		    uv[ij][2] += wgt; // store the wgt 
		  }
		}
	      }
	    }
	  }
	}
      }
      
      // and compute the velocity magnitude...
      double mag_max = 0.0;
      for (int j = 0; j < sy; ++j) {
	for (int i = 0; i < sx; ++i) {
	  const int ij = j*sx+i;
	  if (uv[ij][2] > 0.0) {
	    // normalize...
	    uv[ij][0] /= uv[ij][2];
	    uv[ij][1] /= uv[ij][2];
	    const double mag = uv[ij][0]*uv[ij][0] + uv[ij][1]*uv[ij][1];
	    mag_max = max(mag_max,mag);
	  }
	}
      }
      mag_max = sqrt(mag_max);
      cout << " > max in-plane velocity magnitude: " << mag_max << endl;

      // if the user has NOT specified a scaling, then use 0.5*mag_max as the scale...
      if (!b_scale) {
	scale = 0.5*mag_max;
	cout << " > setting SCALE for this image pair to " << scale << endl;
      }

      // we are done with vimg...
      vimg->lock = false;
      
      // we are going to use the uimg to output the vectors, because 
      // it already has the geometry...
      
      for (int ipx = 0; ipx < nx*ny; ipx++) {
	if (uimg->pixel_flag[ipx] == PngData::VOLUME_DATA_PIXEL) {
	  uimg->pixel_data[ipx] = 0.0;
	}
      }
      
      // start and end of vector...
      double x0[2];
      double x1[2];
      for (int j = 0; j < sy; ++j) {
	for (int i = 0; i < sx; ++i) {
	  const int ij = j*sx+i;
	  if (uv[ij][2] > 0.0) {
	    const double u = uv[ij][0];
	    const double v = uv[ij][1];
	    const double mag = sqrt(u*u + v*v);
	    if (mag > 0.0) {
	      const double theta = atan2(v,u);
	      const double cos_theta = cos(theta);
	      const double sin_theta = sin(theta);
	      // now loop through pixels and compute the min signed distance to the arrow...
	      const double ip0 = (double(i)+0.5)*spacing;
	      const double jp0 = (double(j)+0.5)*spacing;
	      // line segment from x0 to x1...
	      x0[0] = ip0 + 0.5*cos_theta*spacing*mag/scale;
	      x0[1] = jp0 - 0.5*sin_theta*spacing*mag/scale;
	      x1[0] = ip0 - 0.5*cos_theta*spacing*mag/scale;
	      x1[1] = jp0 + 0.5*sin_theta*spacing*mag/scale;
	      // draw points along the arrow shaft...
	      for (int ii = 0; ii <= 20; ++ii) {
		const double wgt = double(ii)/20.0;
		const int ip = int(wgt*x1[0] + (1.0-wgt)*x0[0]);
		const int jp = int(wgt*x1[1] + (1.0-wgt)*x0[1]);
		if ((ip >= 0)&&(ip < nx)&&(jp >= 0)&&(jp < ny)) {
		  const int ipx = jp*nx + ip;
		  if (uimg->pixel_flag[ipx] == PngData::VOLUME_DATA_PIXEL) {
		    uimg->pixel_data[ipx] = 1.0;
		  }
		}
	      }
	      // and draw the arrow head...
	      for (int dip = -1; dip <= 1; ++dip) {
		for (int djp = -1; djp <= 1; ++djp) {
		  const int ip = int(x1[0]+dip);
		  const int jp = int(x1[1]+djp);
		  if ((ip >= 0)&&(ip < nx)&&(jp >= 0)&&(jp < ny)) {
		    const int ipx = jp*nx + ip;
		    if (uimg->pixel_flag[ipx] == PngData::VOLUME_DATA_PIXEL) {
		      uimg->pixel_data[ipx] = 1.0;
		    }
		  }
		}
	      }
	    }
	  }
	}
      }

      // no need to touch the metadata: set the new variable name like this...
      uimg->varId = "vector";
      uimg->rescaleRangesToData();
      uimg->initializeWrite();
      
      char filename[128];
      sprintf(filename,"%s.%06d.png",outputPrefix.c_str(),iP);
      
      uimg->write(filename);
      uimg->lock = false;
      
      cout << " > image pair: " << iP << " of " << nImgPair << ": vector img: " << filename << endl;
      
    }//(iP < nImgPair)
    
    if (uv != NULL) delete[] uv;

    cout << " > done VECTOR" << endl;
    
  }//itpVector()

  void itpCrop(Param * param){

    COUT1("CROP...");

    if ((!pngCache)||(pngCache->imageList.size() == 0)) {
      CERR("cannot CROP: pngCache is empty.");
    }
    
    string outPath = "crop";

    assert(param->size() >= 4);
    const double i0 = param->getInt(0);
    const double j0 = param->getInt(1);
    const double i1 = param->getInt(2);
    const double j1 = param->getInt(3);
    assert(i1 > i0); // etc...
    
    cout << " > i0,j0: " << i0 << " " << j0 << endl;
    cout << " > i1,j1: " << i1 << " " << j1 << endl;

    if (param->size() > 4) {
      outPath = param->getString(4);
    }
    
    cout << " > cropped files will be written to path: " << outPath << endl;

    MiscUtils::mkdir_for_file(outPath+"/test.png");
    for (int iI = 0; iI<pngCache->size(); ++iI){
      PngData * img = pngCache->getImage(iI);
      img->crop(i0,j0,i1,j1);
      string filename = outPath + "/" + MiscUtils::getFilenameNoPath(pngCache->imageList[iI]);
      cout << " > cropping image " << iI << " of " << pngCache->size() << " to output image: " << filename << endl;
      img->initializeWrite();
      img->write(filename.c_str());
      img->lock = false;
    }
     
  } //itpCrop
  
  void itpSuperpose(Param * param) {
    //  Allocate and set unflagged values
    bool    b_thresh      = false;
    bool    b_box         = false;
    bool    b_window      = false;
    bool    b_write       = false;
    bool    b_prefix      = false;
    double  t_lo          = HUGE_VAL;
    double  t_hi          = HUGE_VAL;
    double  box[6]        = {HUGE_VAL, -HUGE_VAL, HUGE_VAL, -HUGE_VAL, HUGE_VAL, -HUGE_VAL};
    double  wind[4]       = {HUGE_VAL, -HUGE_VAL, HUGE_VAL, -HUGE_VAL};
    string  prefix        = "";
    string list[2]        ={"", ""};

    //  Get parameters from the input file
    int iarg = 0;
    cout << "about to parse\n";
    cout << "size = " << param->size() << "\n";
    while (iarg < param->size()) {
      const string token = param->getString(iarg++);
      cout << "Parsing token " << token << "\n";
      if (token == "THRESHOLD") {
        b_thresh = true;
        t_lo = param->getDouble(iarg++);
        t_hi = param->getDouble(iarg++);
      }//"THRESHOLD"
      else if (token == "BOX") {
        b_box = true;
        for (int i = 0; i < 6; i++)
          box[i] = param->getDouble(iarg++);
      }//"BOX"
      else if (token == "WINDOW") {
        b_window = true;
        for (int i = 0; i < 4; i++)
          wind[i] = param->getInt(iarg++);
      }//"WINDOW"
      else if (token == "WRITE_TO_FILE")
        b_write = true;
      else if (token == "PREFIX") {
        b_prefix = true;
        prefix = param->getString(iarg++);
      }//"PREFIX"
      else if (token == "LIST_1") {
        list[0] = param->getString(iarg++);
      }//"LIST_1"
      else if (token == "LIST_2") {
        list[1] = param->getString(iarg++);
      }//"LIST_2"
      else {
        cout << "Unrecognized input argument:" << token << "\n";
        throw(-1);
      }//Unrecognized
    }//(iarg < param->size())

    ////////////////////////////////////////////////////////////////////////////////////////
    //  Error check and interpret the IO bounds
    ////////////////////////////////////////////////////////////////////////////////////////
    if (b_thresh == false) {
      cerr << "No threshold values set for contour work\n";
      throw(-1);
    }//(!b_thresh)
    else if (t_lo >= t_hi) {
      cerr << "Threshold values (" << t_lo << ", " << t_hi << ") don't make sense\n";
      throw(-1);
    }//(t_lo >= t_hi)
    else
      cout << "Pixels with values between " << t_lo << " and " << t_hi << "will be flagged!\n";

    if ((list[0] == "") || (list[1] == "")) {
      cerr << "must provide two image lists for superposition\n";
      cerr << "LIST_1, the list of images that will serve as the superpose canvas and \n";
      cerr << "LIST_2, the list of images that will provide the contour information for superposition!\n";
      cerr << "The current input file gives LIST_1 = " << list[0] << " and LIST_2 = " << list[1] << "\n";
      throw(-1);
    }//(list is incorrect);
    else  buildTwoPngCaches(list);

    if ((pngCaches[0]->size() != pngCaches[1]->size()) || (pngCaches[0]->size() <= 0) || (pngCaches[1]->size() <= 0)) {
      cerr << "something is wrong with the pngCaches!\n";
      cerr << "cache 0 contains " << pngCaches[0]->size() << " images while cache 1 contains " << pngCaches[1]->size() << "\n";
      throw(-1);
    }

    PngData * mask = pngCaches[0]->getImage(0);

    //  Convert the spatial bounds provided by the BOX command into corresponding pixel bounds
    if (b_box) {
      double i_i[8];
      double i_j[8];
      int count = 0;
      for (int i_x = 0; i_x < 2; i_x ++ ) {
        for (int i_y = 0; i_y < 2; i_y ++) {
          for (int i_z = 0; i_z < 2; i_z ++) {
            double x_in[3] = {box[i_x], box[i_y +2], box[i_z + 4]};
            double x_out[3], val;
            int ij[2];
            mask->probePoint(x_in, x_out, ij, val);
            i_i[count] = ij[0];
            i_j[count] = ij[1];
            count++;
          }//(i_z)
        }//(i_y)
      }//(i_x)
      for (count = 0; count < 8; count++) {
        if (wind[0] > i_i[count]) wind[0] = i_i[count];
        if (wind[1] < i_i[count]) wind[1] = i_i[count];
        if (wind[2] > i_j[count]) wind[2] = i_j[count];
        if (wind[3] < i_j[count]) wind[3] = i_j[count];
      }//(count < 8)
    }//(b_box)
    else if (!b_window){
      wind[0] = 0;
      wind[1] = mask->getNx();
      wind[2] = 0;
      wind[3] = mask->getNy();
      cout << "All pixels are eligible for thresholding\n";
    }

    cout << "Contour pixel bounds are: " << wind[0] << " < i < " << wind[1] << ", " << wind[2] << " < j < " << wind[3] << "\n";

    if (!b_prefix) prefix = "superpose";
    cout << "Outputs will be named: " << prefix << ".LOREM_IPSUM.png\n";

    for (int iI = 0; iI < pngCaches[0]->size(); ++iI) {
      cout << "working on images "
           << pngCaches[0]->imageList[iI] << "and"
           << pngCaches[1]->imageList[iI] << "\n";

      string flaggedFile = prefix + "." + pngCache->imageList[iI] + ".flagged.txt";

      ofstream flaggedStream;
      if (b_write) {
        flaggedStream.open(flaggedFile.c_str());
        cout << "dumping flagged locations into " << flaggedFile.c_str() << "\n";
      }//(b_write)
      int color = PngData::FLAG_RED2;
      PngData * img0 = pngCaches[0]->getImage(iI);
      PngData * img1 = pngCaches[1]->getImage(iI);
      //    Make sure both images are the same size
      assert(img0->getNx() == img1->getNx());
      assert(img0->getNy() == img1->getNy());
      int ipx = -1;
      for (int jp = 0; jp < img0->getNy(); ++jp) {//loop through y pixels
        for (int ip = 0; ip < img0->getNx(); ++ip) {//loop through x pixels
          ipx++;
          if ((img1->pixel_flag[ipx] >= 0)         //pixel contains data
              && (ip > wind[0]) && (ip < wind[1])    //pixel falls within x bounds
              && (jp > wind[2]) && (jp < wind[3])    //pixel falls within y bounds
              && (img1->pixel_data[ipx] > t_lo)
              && (img1->pixel_data[ipx] < t_hi)) {         //value is flaggable
            img0->pixel_flag[ipx] = color;
            if (b_write) {
              double flaggedX[3];
              img1->convertPxtoXp(ip, jp, flaggedX);
              flaggedStream << flaggedX[0] << "\t" << flaggedX[1] << "\t" << flaggedX[2] << "\n";
            }//(b_write)
          }//(flaggable)
        }//(ip < getNx()))
      }//(jp < getNy())
      if (b_write) flaggedStream.close();
      string contoured = prefix + "." + pngCaches[0]->imageList[iI];
      img0->initializeWrite();
      img0->write(contoured.c_str());
      img0->lock = false;
      img1->lock = false;
    }//(iI < pngCache->size())
  }//itpSuperpose()

  void dumpExtrema(PngData * img) {
    double x_test[3];
    img->convertPxtoXp(0, 0, x_test);
    cout << "0 0 (" << x_test[0] << ", " << x_test[1] << ", " << x_test[2] << ")\n";
    img->convertPxtoXp(img->getNx()-1, 0, x_test);
    cout << img->getNx() << " 0 (" << x_test[0] << ", " << x_test[1] << ", " << x_test[2] << ")\n";
    img->convertPxtoXp(0, img->getNy()-1, x_test);
    cout << "0 " << img->getNy() <<  "(" << x_test[0] << ", " << x_test[1] << ", " << x_test[2] << ")\n";
    img->convertPxtoXp(img->getNx()-1, img->getNy()-1, x_test);
    cout << img->getNx() << " " << img->getNy() << " (" << x_test[0] << ", " << x_test[1] << ", " << x_test[2] << ")\n";
  }//dumpExtrema()

  void itpIsocontours(Param * param, int mode) {
    // mode = 0 ISOCONTOUR   -> single contour level on multiple images
    // mode = 1 ISOCONTOURS  -> multiple contours levels on single image
    if      (mode == 0) cout << "ISOCONTOUR...\n";
    else if (mode == 1) cout << "ISOCONTOURS...\n";

    if ((!pngCache) || (pngCache->imageList.size() == 0)) {
      cout << "cannot compute ISOCONT: pngCache is empty.\n";
      throw(-1);
    }
    else
      cout << "pngCache checks out!\n";

    // Allocate and set unflagged valeus
    bool    b_thresh      = false;
    bool    b_box         = false;
    bool    b_window      = false;
    bool    b_write       = false;
    bool    b_write_mask  = false;
    double  t_lo          = HUGE_VAL;
    double  t_hi          = HUGE_VAL;
    double  box[6]        = {HUGE_VAL, -HUGE_VAL, HUGE_VAL, -HUGE_VAL, HUGE_VAL, -HUGE_VAL};
    double  wind[4]       = {HUGE_VAL, -HUGE_VAL, HUGE_VAL, -HUGE_VAL};

    // Get parameters from the input file
    int iarg = 0;
    while (iarg < param->size()) {
      const string token = param->getString(iarg++);
      cout << "Parsing token " << token << "\n";
      if (token == "THRESHOLD") {
        b_thresh = true;
        t_lo = param->getDouble(iarg++);
        t_hi = param->getDouble(iarg++);
      }//"THRESHOLD"
      else if (token == "BOX") {
        b_box = true;
        for (int i = 0; i < 6; i++)
          box[i] = param->getDouble(iarg++);
      }//"BOX"
      else if (token == "WINDOW") {
        b_window = true;
        for (int i = 0; i < 4; i++)
          wind[i] = param->getInt(iarg++);
      }//"WINDOW"
      else if (token == "WRITE_TO_FILE")
        b_write = true;
      else if (token == "WRITE_MASK") {
        b_write_mask = true;
      }//"WRITE_MASK"
      else {
        cout << "Unrecognized input argument:" << token << "\n";
        throw(-1);
      }//Unrecognized
    }//(iarg < param->size())

    ////////////////////////////////////////////////////////////////////////////////////////
    //  Build the geometry + bahama blue mask
    ////////////////////////////////////////////////////////////////////////////////////////
    PngData * mask = new PngData();
    mask->read(pngCache->imageList[0].c_str());
    mask->finalizeRead();
    int npx = mask->getNx() * mask->getNy();
    //  Null out all of the data
    for (int ipx = 0; ipx < npx; ++ipx)
      if (mask->pixel_flag[ipx] >= 0)
        mask->pixel_flag[ipx] = PngData::FLAG_WHITE;

    if (mode == 0) {
      MiscUtils::mkdir_for_file("./isocontour");
      cout << "building dir isocontour" << endl;
    }
    else if (mode == 1) {
      MiscUtils::mkdir_for_file("./isocontours");
      cout << "building dir isocontours" << endl;
    }
    else {
      cout << "crashing" << endl;
      throw(-1);
    }

    if (b_write_mask) {
      mask->initializeWrite();
      if (mode == 0) mask->write("./isocontour/mask.png");
      else if (mode == 1) mask->write("./isocontours/mask.png");
    }//(b_write_mask)

    //  Dump extrema
    if (1 == 1) dumpExtrema(mask);
    ////////////////////////////////////////////////////////////////////////////////////////
    //  Error check and interpret the IO bounds
    ////////////////////////////////////////////////////////////////////////////////////////
    if (b_thresh == false) {
      cerr << "No threshold values set for contour work\n";
      throw(-1);
    }//(!b_thresh)
    else if (t_lo >= t_hi) {
      cerr << "Threshold values ( " << t_lo << ", " << t_hi << " ) don't make sense\n";
      throw(-1);
    }//(t_lo >= t_hi)
    else
      cout << "Pixels with values between " << t_lo << " and " << t_hi << " will be flagged!\n";

    //  Dump extrema
    if (1 == 2) dumpExtrema(mask);

    //  Convert the spatial bounds provided by the BOX command into corresponding pixel bounds
    if (b_box) {
      double i_i[8];
      double i_j[8];
      int count = 0;
      for (int i_x = 0; i_x < 2; i_x ++ ) {
        for (int i_y = 0; i_y < 2; i_y ++) {
          for (int i_z = 0; i_z < 2; i_z ++) {
            double x_in[3] = {box[i_x], box[i_y +2], box[i_z + 4]};
            double x_out[3], val;
            int ij[2];
            mask->probePoint(x_in, x_out, ij, val);
            i_i[count] = ij[0];
            i_j[count] = ij[1];
            count++;
          }//(i_z)
        }//(i_y)
      }//(i_x)
      for (count = 0; count < 8; count++) {
        if (wind[0] > i_i[count]) wind[0] = i_i[count];
        if (wind[1] < i_i[count]) wind[1] = i_i[count];
        if (wind[2] > i_j[count]) wind[2] = i_j[count];
        if (wind[3] < i_j[count]) wind[3] = i_j[count];
      }//(count < 8)
    }//(b_box)
    else if (!b_window) {
      wind[0] = 0;
      wind[1] = mask->getNx();
      wind[2] = 0;
      wind[3] = mask->getNy();
      cout << "All pixels are eligible for thresholding\n";
    }//(!b_window)

    cout << "Contour pixel bounds are: " << wind[0] << " < i < " << wind[1] << ", " << wind[2] << " < j < " << wind[3] << "\n";

    ////////////////////////////////////////////////////////////////////////////////////////
    //  Single contour on multiple images
    ////////////////////////////////////////////////////////////////////////////////////////
    if (mode == 0) {
      for (int iI = 0; iI < pngCache->size(); ++iI) {
        cout << "working on image " << pngCache->imageList[iI] << "\n";
        string tmp_in = pngCache->imageList[iI];
        string tmp_out = tmp_in.substr(tmp_in.find_last_of("\\/") + 1, tmp_in.length());

        string flaggedFile = "./isocontours/isoc." + tmp_out + ".flagged.txt";
        ofstream flaggedStream;
        if (b_write) {
          flaggedStream.open(flaggedFile.c_str());
          cout << "dumping flagged locations into " << flaggedFile.c_str() << "\n";
        }//(b_write)
        int color = PngData::FLAG_RED2;
        PngData * img = pngCache->getImage(iI);
        assert(npx == img->getNx() * img->getNy());
        int ipx = -1;
        for (int jp = 0; jp < img->getNy(); ++jp) {//loop through y pixels
          for (int ip = 0; ip < img->getNx(); ++ip) {//loop through x pixels
            ipx++;
            if ((img->pixel_flag[ipx] >= 0)          //pixel contains data
                && (ip > wind[0]) && (ip < wind[1])    //pixel falls within x bounds
                && (jp > wind[2]) && (jp < wind[3])    //pixel falls within y bounds
                && (img->pixel_data[ipx] > t_lo)
                && (img->pixel_data[ipx] < t_hi)) {            //value is flaggable
              img->pixel_flag[ipx] = color;
              if (b_write) {
                double flaggedX[3];
                img->convertPxtoXp(ip, jp, flaggedX);
                flaggedStream << flaggedX[0] << "\t" << flaggedX[1] << "\t" << flaggedX[2] << "\n";
              }//(b_write)
            }//(flaggable)
            else //Take value from mask
              img->pixel_flag[ipx] = mask->pixel_flag[ipx];
          }//(jp < getNy())
        }//(ip < getNx()))
        if (b_write) flaggedStream.close();
        string contoured = "./isocontour/isoc." + tmp_out;
        img->initializeWrite();
        img->write(contoured.c_str());
        img->lock = false;
      }//(iI < pngCache->size())
    }//(mode == 0)

    ////////////////////////////////////////////////////////////////////////////////////////
    //  Multiple contours on single image
    ////////////////////////////////////////////////////////////////////////////////////////
    else if (mode == 1) {
      if (pngCache->size() > 10)
        cout << "WARNING: Isocontours can only show 10 levels at most. Only the first 10 will be shown\n";

      string tmp_in = pngCache->imageList[0];
      string tmp_out = tmp_in.substr(tmp_in.find_last_of("\\/") + 1, tmp_in.length());
      for (int iI = 0; iI < min(pngCache->size(),10); ++iI) {
        cout << "working on image " << pngCache->imageList[iI] << "\n";
        string flaggedFile = "";
        ofstream flaggedStream;
        if (b_write) {
          flaggedStream.open(flaggedFile.c_str());
          cout << "dumping flagged locations into " << flaggedFile.c_str() << "\n";
        }//(b_write)
        int color = PngData::FLAG_RED2 - iI;
        PngData * img = pngCache->getImage(iI);
        assert(npx == img->getNx() * img->getNy());
        int ipx = -1;
        for (int jp = 0; jp < img->getNy(); ++jp) {//loop through y pixels
          for (int ip = 0; ip < img->getNx(); ++ip) {//loop through x pixels
            ipx++;
            if ((img->pixel_flag[ipx] >= 0)          //pixel contains data
                && (ip > wind[0]) && (ip < wind[1])    //pixel falls within x bounds
                && (jp > wind[2]) && (jp < wind[3])    //pixel falls within y bounds
                && (img->pixel_data[ipx] > t_lo)
                && (img->pixel_data[ipx] < t_hi)) {             //value is flaggable
              mask->pixel_flag[ipx] = color;
              if (b_write) {
                double flaggedX[3];
                img->convertPxtoXp(ip, jp, flaggedX);
                flaggedStream << flaggedX[0] << "\t" << flaggedX[1] << "\t" << flaggedX[2] << "\n";
              }//(b_write)
            }//(flagged)
          }//(ip < Nx)
        }//(jp < Ny)
        img->lock = false;
        if (b_write) flaggedStream.close();
      }//(iI < min(pngCache->size(),10))

      string contoured = "./isocontours/isoc." + tmp_out;
      mask->initializeWrite();
      mask->write(contoured.c_str());
      mask->lock = false;
    }//(mode == 1)

  }//itpIsocontours()

  //flag contiguous pixels to a given point,
  //compute stats on these pixels, per image.
  //Similar to SPATIAL_STATS but on a pixel subset of the image
  void itpIntegrationProbe(Param * param){

    COUT1("INTEGRATION_PROBE...");

    if ((!pngCache)||(pngCache->imageList.size() == 0)) {
      CERR("cannot compute INTEGRATION_PROBE: pngCache is empty.");
    }

    int seed_i, seed_j, seed_pixel;
    if (param->size() > 2 && param->getString(0)=="PIXEL"){
      seed_i = param->getInt(1);
      seed_j = param->getInt(2);
    }
    else{
      CERR("cannot compute INTEGRATION_PROBE: specify a seed pixel with INTEGRATION_PROBE PIXEL <i> <j>");
    }

    string outputName;
    if (param->size() > 3) {
      outputName = param->getString(3);
    }
    else {
      outputName = "integration_probe.dat";
    }

    cout << " > writing integration probe to file \"" << outputName << "\"..." << endl;

    ofstream ofile;
    size_t fwidth = 13;
    MiscUtils::mkdir_for_file(outputName);
    ofile.open(outputName.c_str());
    if (!ofile.is_open()) {
      CERR("could not open file: " << outputName);
    }

    // the parameter...
    ofile << "# " << param->str() << "\n";

    // the variable name...
    ofile << "# " <<
      setw(fwidth-2) << "1:image" <<
      setw(fwidth) << "2:area" <<
      setw(fwidth) << "3:centroid-x" <<
      setw(fwidth) << "4:centroid-y" <<
      setw(fwidth) << "5:centroid-z" <<
      setw(fwidth) << "6:avg" <<
      setw(fwidth) << "7:rms" <<
      setw(fwidth) << "8:min" <<
      setw(fwidth) << "9:max" ;
    ofile << endl;

    //COUT1("Spatial Image Stats: Name, Min, Max, Avg, RMS");
    for (int iI = 0; iI<pngCache->size(); ++iI){
      PngData * img = pngCache->getImage(iI);

      seed_pixel = img->getNx()*seed_j + seed_i;
      if (seed_pixel >= img->getNx()*img->getNy()) {
        CERR("Seed pixel value out of range");
      }

      img->flagContiguousPixels(seed_pixel);

      double stats[8]; //x,y,z,min,max,avg,rms,area
      img->computeStats(stats,PngData::FLAG_GREEN);
      //COUT1("Image " << iI << " " << stats[0] << ", " << stats[1] << ", " << stats[2] << ", " << stats[3]);
      img->lock = false;
      ofile << "  " <<
        setw(fwidth-2) << iI <<
        setw(fwidth) << stats[7] <<
        setw(fwidth) << stats[0] <<
        setw(fwidth) << stats[1] <<
        setw(fwidth) << stats[2] <<
        setw(fwidth) << stats[5] << // note order switch to match CONDITIONAL_PROBE output
        setw(fwidth) << stats[6] <<
        setw(fwidth) << stats[3] <<
        setw(fwidth) << stats[4] ;
      ofile << endl;

      string outPath = "./images_cropped";
      COUT1("Writing cropped files to " << outPath);

      MiscUtils::mkdir_for_file(outPath+"/test.png");

      //build output path and filename
      //use original filename in a different folder
      vector<string> inputNameTokens;
      tokenizeString(inputNameTokens,pngCache->imageList[iI],"/");
      string outputPathAndFilename = outPath + "/";
      outputPathAndFilename += inputNameTokens[inputNameTokens.size()-1];

      img->rescaleRangesToData();
      img->initializeWrite();
      img->write(outputPathAndFilename.c_str());
      img->lock = false;
    }//(iI < pngCache->size())
    ofile.close();
  }//itpIntegrationProbe()


  //Report radial profile of azimuthally averaged data
  void itpRadialProbe(Param * param){

    // see newer routine RADIAL_PROFILE for outlet profile stats.

    COUT1("RADIAL_PROBE...");

    if ((!pngCache)||(pngCache->imageList.size() == 0)) {
      CERR("cannot compute RADIAL_PROBE: pngCache is empty.");
    }
    pngCache->clearImageCache(); //If multiple RADIAL_PROBES are requested, must re-read the images
                                 //to clear flagging of marked pixels.

    int seed_i, seed_j, seed_pixel;
    double radius;
    if (param->size() > 3 && param->getString(0)=="PIXEL"){
      seed_i = param->getInt(1);
      seed_j = param->getInt(2);
      radius = param->getDouble(3);
    }
    else{
      CERR("cannot compute RADIAL_PROBE: specify a seed pixel with RADIAL_PROBE PIXEL <i> <j> <radius>");
    }
    if (radius<0.0){
      CERR("RADIAL_PROBE radius must be positive: " << radius);
    }
    stringstream ss;
    ss << "i" << seed_i << "j" << seed_j;

    for (int iI = 0; iI<pngCache->size(); ++iI){
      PngData * img = pngCache->getImage(iI);

      COUT1(" > working on image " << iI << " of " << pngCache->size() << " \"" << pngCache->imageList[iI] << "\"...");

      seed_pixel = img->getNx()*seed_j + seed_i;
      if (seed_pixel >= img->getNx()*img->getNy()){
        COUT1("Seed location pixel out of image dimension size, skipping.");
        continue;
      }

      vector<pixel_type> pixelTypeVec;
      img->getPixelTypes(pixelTypeVec);
      if (pixelTypeVec.size()<1){
        COUT1("No data pixels found in image, skipping.");
        continue;
      }

      int ipt0 = -1;
      for (int ipt=0; ipt<pixelTypeVec.size(); ++ipt){
        if (img->pixel_flag[seed_pixel] == pixelTypeVec[ipt].flag){
          COUT1(" > seed location pixel has data type " << pixelTypeVec[ipt].name);
          ipt0 = ipt;
        }
      }
      if (ipt0<0){
        COUT1("Seed location does not fall on a data pixel, skipping");
        continue;
      }


      vector<string> inputNameTokens;
      tokenizeString(inputNameTokens,pngCache->imageList[iI],"/");
      string temp = inputNameTokens[inputNameTokens.size()-1];
      inputNameTokens.clear();
      tokenizeString(inputNameTokens,temp,".");
      string outputDat = "./";
      string outputPrefix = "";
      for (int ii=0;ii<inputNameTokens.size()-1;++ii){ //execlude .png
        outputPrefix += inputNameTokens[ii];
        outputPrefix += ".";
      }
      outputPrefix += ss.str();
      outputDat = "./" + outputPrefix + ".dat";

      cout << " > writing radial probe to file \"" << outputDat << "\"..." << endl;

      ofstream ofile;
      size_t fwidth = 13;
      MiscUtils::mkdir_for_file(outputDat);
      ofile.open(outputDat.c_str());
      if (!ofile.is_open()) {
        CERR("could not open file: " << outputDat);
      }

      // the parameter...
      ofile << "# " << param->str() << "\n";
      ofile << "# data type - " << pixelTypeVec[ipt0].name << "; varId - " << pixelTypeVec[ipt0].varId << endl;

      // the variable name...
      ofile << "# " <<
        setw(fwidth-2) << "1:index" <<
        setw(fwidth) << "2:R (pxl.)" <<
        setw(fwidth) << "3:R (Sim.)" <<
        setw(fwidth) << "4:avg" <<
        setw(fwidth) << "5:rms" <<
        setw(fwidth) << "6:min" <<
        setw(fwidth) << "7:max" ;
      ofile << endl;

      double * no_data = new double[(img->getNx()+1)*(img->getNy()+1)];
      img->buildNodalData(no_data,img->pixel_flag[seed_pixel]);
      double rad = 0.5;
      int index = 0;

      //store probe region, marking in the loop impacts
      //the interpolation, update pixel_flag before
      //writing image
      int * flag_buf = new int[img->getNx()*img->getNy()];
      for (int ipx=0; ipx<img->getNx()*img->getNy(); ++ipx)
        flag_buf[ipx] = PngData::FLAG_BLANK;

      while (rad<=radius){
        int pixel_flag_val = PngData::FLAG_RED;
        if (rad<=5)
          pixel_flag_val = PngData::FLAG_GREEN; //mark center in green

        double stats[6]; //rad, rad, avg, rms, min, max
        double theta = 0.0;
        double min=HUGE_VAL;
        double max=-HUGE_VAL;
        double avg_buf=0.0;
        double rms_buf=0.0;
        int count = 0;
        while (theta < 2.0*M_PI){
          double ipf = seed_i + (double) rad*cos(theta);
          double jpf = seed_j + (double) rad*sin(theta);
          double val;
          if (img->interpFromNodalData(val,ipf,jpf,no_data,img->pixel_flag[seed_pixel])){
            int nn = ((int) floor(jpf+0.5))*img->getNx() + (int)floor(ipf+0.5);
            flag_buf[nn] = pixel_flag_val;
            if (val<min)
              min = val;
            if (val>max)
              max = val;
            avg_buf += val;
            rms_buf += val*val;
            ++count;
          }
          theta += 0.5/(rad);
        }
        stats[0] = rad;
        stats[1] = rad*img->getLengthScale();
        stats[2] = avg_buf/(double)count;
        stats[3] = sqrt(rms_buf/(double)count - stats[2]*stats[2]);
        stats[4] = min;
        stats[5] = max;

        ofile << "  " <<
          setw(fwidth-2) << index++ <<
          setw(fwidth) << stats[0] <<
          setw(fwidth) << stats[1] <<
          setw(fwidth) << stats[2] <<
          setw(fwidth) << stats[3] <<
          setw(fwidth) << stats[4] <<
          setw(fwidth) << stats[5] ;
        ofile << endl;
        rad += 0.5;
      }
      ofile.close();
      delete[] no_data;

      string outPath = "./images_flagged";
      COUT1("Writing flagged pixels to " << outPath);

      MiscUtils::mkdir_for_file(outPath+"/test.png");

      //build output path and filename
      //use original filename in a different folder
      string outputPathAndFilename = outPath + "/";
      outputPathAndFilename += outputPrefix;
      outputPathAndFilename += ".png";

      //mark probe center green, full region red
      for (int ipx=0; ipx<img->getNx()*img->getNy(); ++ipx){
        if (flag_buf[ipx] != PngData::FLAG_BLANK){
          img->pixel_flag[ipx] = flag_buf[ipx];
        }
      }
      delete[] flag_buf;

      img->rescaleRangesToData();
      img->initializeWrite();
      img->write(outputPathAndFilename.c_str());
      img->lock = false;
    }
  }//itpRadialProbe()

  //Initial implementation of temporal stats on an
  //image sequence.  Currently assumes the
  //that all images capture the same geometry.
  void itpTemporalStats(Param * param){

    COUT1("TEMPORAL_STATS...");

    if ((!pngCache)||(pngCache->imageList.size() == 0)) {
      CERR("cannot compute TEMPORAL_STATS: pngCache is empty.");
    }

    string outputPathAndPrefix;
    string outputPath = "images_stats/";

    if (param->size() > 0) {
      // assume this string is an output path...
      const string token = param->getString(0);
      outputPath = param->getString(0)+"/";
    }

    // try to figure out an output prefix from the input filename...
    outputPathAndPrefix = outputPath + getPrefix(pngCache->imageList[0]);
    COUT1(" > output file path and prefix will be \"" << outputPathAndPrefix << "\"");
    MiscUtils::mkdir_for_file(outputPath+"/test.png");

    PngData * imgOut = new PngData();

    imgOut->read(pngCache->imageList[0].c_str());
    imgOut->finalizeRead();

    vector<pixel_type> pixelTypeVec;
    imgOut->getPixelTypes(pixelTypeVec);
    if (pixelTypeVec.size()<1){
      CERR("no data pixels found in first image");
    }

    int nPx = imgOut->getNx()*imgOut->getNy();

    double * count_buf = new double[nPx];
    double * avg_buf = new double[nPx];
    double * rms_buf = new double[nPx];
    double * min_buf = new double[nPx];
    double * max_buf = new double[nPx];
    for (int iPx = 0; iPx<nPx; ++iPx){
      count_buf[iPx] = 0.0;
      avg_buf[iPx] = 0.0;
      rms_buf[iPx] = 0.0;
      min_buf[iPx] = HUGE_VAL;
      max_buf[iPx] = -HUGE_VAL;
    }

    const int nDataTypes = pixelTypeVec.size();
    map<int,double> all_img_count;
    map<int,double> all_img_avg;
    map<int,double> all_img_rms;
    map<int,double> all_img_min;
    map<int,double> all_img_max;
    for (int iDataTypes=0;iDataTypes<nDataTypes;++iDataTypes){
      all_img_count[pixelTypeVec[iDataTypes].flag] = 0.0;
      all_img_avg[pixelTypeVec[iDataTypes].flag] = 0.0;
      all_img_rms[pixelTypeVec[iDataTypes].flag] = 0.0;
      all_img_min[pixelTypeVec[iDataTypes].flag] = HUGE_VAL;
      all_img_max[pixelTypeVec[iDataTypes].flag] = -HUGE_VAL;
    }
    map<int,double> img_count;
    map<int,double> img_avg;
    map<int,double> img_rms;
    map<int,double> img_min;
    map<int,double> img_max;

    // for output formatting...

    for (int iI = 0; iI<pngCache->size(); ++iI){

      COUT1(" > working on image " << iI << " of " << pngCache->size() << " \"" << pngCache->imageList[iI] << "\"...");

      PngData * img = pngCache->getImage(iI);
      assert(img->getNx() == imgOut->getNx()); // if you hit one of these, images are different sizes -- need to handle in some future version
      assert(img->getNy() == imgOut->getNy());

      for (int iDataTypes=0;iDataTypes<nDataTypes;++iDataTypes){
        img_count[pixelTypeVec[iDataTypes].flag] = 0.0;
        img_avg[pixelTypeVec[iDataTypes].flag] = 0.0;
        img_rms[pixelTypeVec[iDataTypes].flag] = 0.0;
        img_min[pixelTypeVec[iDataTypes].flag] = HUGE_VAL;
        img_max[pixelTypeVec[iDataTypes].flag] = -HUGE_VAL;
      }

      for (int iPx = 0; iPx<nPx; ++iPx) {
        if ( imgOut->pixel_flag[iPx] != img->pixel_flag[iPx] ) {
          CERR("images have different blanking. Consider TEMPORAL_STATS_MOVING or TEMPORAL_STATS_PHASE");
        }
        if (imgOut->pixel_flag[iPx] >= PngData::VOLUME_DATA_PIXEL) {
          const double val = img->pixel_data[iPx];

          count_buf[iPx] += 1.0; // every pixel gets equal weight for now
          avg_buf[iPx]  += val;
          rms_buf[iPx]   += val*val; // for now use rms to store sum(value^2)

          img_count[img->pixel_flag[iPx]]      += 1.0;
          img_avg[img->pixel_flag[iPx]]       += val;
          img_rms[img->pixel_flag[iPx]]        += val*val;

          if (val < min_buf[iPx]) min_buf[iPx] = val;
          if (val > max_buf[iPx]) max_buf[iPx] = val;

          if (val < img_min[img->pixel_flag[iPx]]) img_min[img->pixel_flag[iPx]] = val;
          if (val > img_max[img->pixel_flag[iPx]]) img_max[img->pixel_flag[iPx]] = val;
        }//(imgOut->pixel_flag[iPx] >= 0)
      }//for (iPx < nPx)
      img->lock = false;

      for (int iDataTypes=0;iDataTypes<nDataTypes;++iDataTypes){
        int dataFlag = pixelTypeVec[iDataTypes].flag;
        all_img_count[dataFlag] += img_count[dataFlag];
        all_img_avg[dataFlag]  += img_avg[dataFlag];
        all_img_rms[dataFlag]   += img_rms[dataFlag];

        if (img_min[dataFlag] < all_img_min[dataFlag]) all_img_min[dataFlag] = img_min[dataFlag];
        if (img_max[dataFlag] > all_img_max[dataFlag]) all_img_max[dataFlag] = img_max[dataFlag];


        if (img_count[dataFlag] == 0.0) cout << " > image \"" << pngCache->imageList[iI] << "\" has no active " << pixelTypeVec[iDataTypes].name << " data pixels." << endl;
        else {
          img_avg[dataFlag] /= img_count[dataFlag];
          img_rms[dataFlag] = sqrt(max(0.0,img_rms[dataFlag]/img_count[dataFlag]-img_avg[dataFlag]*img_avg[dataFlag]));
          cout << " > image \"" << pngCache->imageList[iI] << "\" - " << pixelTypeVec[iDataTypes].name << " data - avg: " << img_avg[dataFlag] << " rms: " << img_rms[dataFlag] << " min: " << img_min[dataFlag] << " max: " << img_max[dataFlag] << endl;
        }
      }

    }//for (iI < iI<pngCache->size())

    for (int iDataTypes=0;iDataTypes<nDataTypes;++iDataTypes){
      int dataFlag = pixelTypeVec[iDataTypes].flag;

      if (all_img_count[dataFlag] == 0.0) cout << " > all images have no active " << pixelTypeVec[iDataTypes].name << " data pixels." << endl;
      else {
        all_img_avg[dataFlag] /= all_img_count[dataFlag];
        all_img_rms[dataFlag] = sqrt(max(0.0,all_img_rms[dataFlag]/all_img_count[dataFlag]-all_img_avg[dataFlag]*all_img_avg[dataFlag]));
        cout << " > all images - " << pixelTypeVec[iDataTypes].name << " data - avg: " << all_img_avg[dataFlag] << " rms: " << all_img_rms[dataFlag] << " min: " << all_img_min[dataFlag] << " max: " << all_img_max[dataFlag] << endl;
      }
    }

    for (int iPx = 0; iPx<nPx; ++iPx){
      if ( imgOut->pixel_flag[iPx] >= PngData::VOLUME_DATA_PIXEL ) {

        const double avg = avg_buf[iPx]/count_buf[iPx];
        const double rms = sqrt(max(0.0,rms_buf[iPx]/count_buf[iPx]-avg*avg));
        rms_buf[iPx] = rms;

        imgOut->pixel_data[iPx] = avg;
      }
    }

    string avg_str = outputPathAndPrefix+".avg.png";
    cout << " > writing avg to \"" << avg_str << "\"..." << endl;
    string baseVarId = imgOut->varId;
    string baseVarIdSurface = imgOut->varIdSurface;
    string baseVarIdParticles = imgOut->varIdParticles;
    string baseVarIdIso = imgOut->varIdIso;
    imgOut->varId = baseVarId+"_avg";
    imgOut->varIdSurface = baseVarIdSurface+"_avg";
    imgOut->varIdParticles = baseVarIdParticles+"_avg";
    imgOut->varIdIso = baseVarIdIso+"_avg";

    imgOut->rescaleRangesToData();
    imgOut->initializeWrite();
    imgOut->write(avg_str.c_str());

    //re-init imgOut; initializeWrite() changes the rgb and data chunks used
    //to construct the lighting estimate.  To get the same lighting estimate
    //for avg, rms, min and max images must grab the original buffers
    delete imgOut;
    imgOut = new PngData();
    imgOut->read(pngCache->imageList[0].c_str());
    imgOut->finalizeRead();

    for (int iPx = 0; iPx<nPx; ++iPx){
      if ( imgOut->pixel_flag[iPx] >= PngData::VOLUME_DATA_PIXEL ) {
        imgOut->pixel_data[iPx] = rms_buf[iPx];
      }
    }
    string rms_str = outputPathAndPrefix+".rms.png";
    cout << " > writing rms to \"" << rms_str << "\"..." << endl;

    imgOut->varId = baseVarId+"_rms";
    imgOut->varIdSurface = baseVarIdSurface+"_rms";
    imgOut->varIdParticles = baseVarIdParticles+"_rms";
    imgOut->varIdIso = baseVarIdIso+"_rms";

    imgOut->rescaleRangesToData();
    imgOut->initializeWrite();
    imgOut->write(rms_str.c_str());

    //re-init imgOut; initializeWrite() changes the rgb and data chunks used
    //to construct the lighting estimate.  To get the same lighting estimate
    //for avg, rms, min and max images must grab the original buffers
    delete imgOut;
    imgOut = new PngData();
    imgOut->read(pngCache->imageList[0].c_str());
    imgOut->finalizeRead();

    for (int iPx = 0; iPx<nPx; ++iPx){
      if ( imgOut->pixel_flag[iPx] >= PngData::VOLUME_DATA_PIXEL ) {
        imgOut->pixel_data[iPx] = min_buf[iPx];
      }
    }
    string min_str = outputPathAndPrefix+".min.png";
    cout << " > writing min to \"" << min_str << "\"..." << endl;

    imgOut->varId = baseVarId+"_min";
    imgOut->varIdSurface = baseVarIdSurface+"_min";
    imgOut->varIdParticles = baseVarIdParticles+"_min";
    imgOut->varIdIso = baseVarIdIso+"_min";

    imgOut->rescaleRangesToData();
    imgOut->initializeWrite();
    imgOut->write(min_str.c_str());

    //re-init imgOut; initializeWrite() changes the rgb and data chunks used
    //to construct the lighting estimate.  To get the same lighting estimate
    //for avg, rms, min and max images must grab the original buffers
    delete imgOut;
    imgOut = new PngData();
    imgOut->read(pngCache->imageList[0].c_str());
    imgOut->finalizeRead();

    for (int iPx = 0; iPx<nPx; ++iPx){
      if ( imgOut->pixel_flag[iPx] >= PngData::VOLUME_DATA_PIXEL ) {
        imgOut->pixel_data[iPx] = max_buf[iPx];
      }
    }
    string max_str = outputPathAndPrefix+".max.png";
    cout << " > writing max to \"" << max_str << "\"..." << endl;

    imgOut->varId = baseVarId+"_max";
    imgOut->varIdSurface = baseVarIdSurface+"_max";
    imgOut->varIdParticles = baseVarIdParticles+"_max";
    imgOut->varIdIso = baseVarIdIso+"_max";

    imgOut->rescaleRangesToData();
    imgOut->initializeWrite();
    imgOut->write(max_str.c_str());

    delete imgOut;
    delete[] count_buf;
    delete[] avg_buf;
    delete[] rms_buf;
    delete[] min_buf;
    delete[] max_buf;

  }//itpTemporalStats()

  void itpTemporalStatsCross(Param * param){

    COUT1("TEMPORAL_STATS_CROSS...");

    int imageCount = pngCache->imageList.size();

    if ((!pngCache)||(imageCount == 0)) {
      CERR("cannot compute TEMPORAL_STATS_CROSS: pngCache is empty.");
    }
    else if (imageCount%2!=0) {
      CERR("cannot compute TEMPORAL_STATS_CROSS, pngCache image count must be even.");
    }

    string outputPath = "images_stats/";

    if (param->size() > 0) {
      // assume this string is an output path...
      const string token = param->getString(0);
      outputPath = param->getString(0)+"/";
    }

    // try to figure out an output prefix from the input filename...
    string prefix0 = getPrefix(pngCache->imageList[0]);
    string prefix1 = getPrefix(pngCache->imageList[imageCount/2]);
    string outputPathAndPrefix0 = outputPath + prefix0;
    string outputPathAndPrefix1 = outputPath + prefix1;
    stringstream crossPathAndPrefixSS;
    crossPathAndPrefixSS << outputPath << prefix0 << "_" << prefix1;
    COUT1(" > the image path and prefixes to be used are \n" << 
          "     \"" << outputPathAndPrefix0 << "\"\n" << 
          "     \"" << outputPathAndPrefix1 << "\"\n" << 
          "     \"" << crossPathAndPrefixSS.str() << "\"");
    MiscUtils::mkdir_for_file(outputPath+"/test.png");

    PngData * imgOut0 = new PngData();
    PngData * imgOut1 = new PngData();

    imgOut0->read(pngCache->imageList[0].c_str());
    imgOut0->finalizeRead();
    imgOut1->read(pngCache->imageList[imageCount/2].c_str());
    imgOut1->finalizeRead();

    vector<pixel_type> pixelTypeVec;
    imgOut0->getPixelTypes(pixelTypeVec);
    if (pixelTypeVec.size()<1){
      CERR("no data pixels found in first image type");
    }

    vector<pixel_type> pixelTypeVec1;
    imgOut1->getPixelTypes(pixelTypeVec1);
    if (pixelTypeVec1.size()<1){
      CERR("no data pixels found in second image type");
    }

    int nPx = imgOut0->getNx()*imgOut0->getNy();
    //int nPx1 = imgOut1->getNx()*imgOut1->getNy();
    if (imgOut0->getNx()!=imgOut1->getNx()||
        imgOut0->getNy()!=imgOut1->getNy()){
      CERR("Image dimensions must match");
    }

    double * count_buf = new double[nPx];
    double * avg0_buf = new double[nPx];
    double * rms0_buf = new double[nPx];
    double * avg1_buf = new double[nPx];
    double * rms1_buf = new double[nPx];
    double * cross_buf = new double[nPx];

    for (int iPx = 0; iPx<nPx; ++iPx){
      count_buf[iPx] = 0.0;
      avg0_buf[iPx] = 0.0;
      rms0_buf[iPx] = 0.0;
      avg1_buf[iPx] = 0.0;
      rms1_buf[iPx] = 0.0;
      cross_buf[iPx] = 0.0;
    }

    const int nDataTypes = pixelTypeVec.size();
    const int nDataTypes1 = pixelTypeVec1.size();
    if (nDataTypes!=nDataTypes1) {
      CERR("Image datatype counts must match");
    }
    for (int iDataTypes=0;iDataTypes<nDataTypes;++iDataTypes){
      if (pixelTypeVec[iDataTypes].flag != pixelTypeVec1[iDataTypes].flag){
        CERR("Image datatypes must match");
      }
    }

    map<int,double> all_img_count;
    map<int,double> all_img_avg0;
    map<int,double> all_img_rms0;
    map<int,double> all_img_avg1;
    map<int,double> all_img_rms1;
    map<int,double> all_img_cross;

    for (int iDataTypes=0;iDataTypes<nDataTypes;++iDataTypes){
      all_img_count[pixelTypeVec[iDataTypes].flag] = 0.0;
      all_img_avg0[pixelTypeVec[iDataTypes].flag] = 0.0;
      all_img_rms0[pixelTypeVec[iDataTypes].flag] = 0.0;
      all_img_avg1[pixelTypeVec[iDataTypes].flag] = 0.0;
      all_img_rms1[pixelTypeVec[iDataTypes].flag] = 0.0;
      all_img_cross[pixelTypeVec[iDataTypes].flag] = 0.0;
    }
    map<int,double> img_count;
    map<int,double> img_avg0;
    map<int,double> img_rms0;
    map<int,double> img_avg1;
    map<int,double> img_rms1;
    map<int,double> img_cross;

    // for output formatting...

    for (int iI = 0; iI<imageCount/2; ++iI){
      int jI = imageCount/2 + iI;

      COUT1(" > working on image pair " << iI << " of " << imageCount/2 << " \"" << pngCache->imageList[iI] << "\" and \"" << pngCache->imageList[jI] << "\"");

      PngData * img0 = pngCache->getImage(iI);
      PngData * img1 = pngCache->getImage(jI);

      assert(img0->getNx() == imgOut0->getNx()); // if you hit one of these, images are different sizes -- need to handle in some future version
      assert(img0->getNy() == imgOut0->getNy());
      assert(img1->getNx() == imgOut1->getNx()); 
      assert(img1->getNy() == imgOut1->getNy());

      for (int iDataTypes=0;iDataTypes<nDataTypes;++iDataTypes){
        img_count[pixelTypeVec[iDataTypes].flag] = 0.0;
        img_avg0[pixelTypeVec[iDataTypes].flag] = 0.0;
        img_rms0[pixelTypeVec[iDataTypes].flag] = 0.0;
        img_avg1[pixelTypeVec[iDataTypes].flag] = 0.0;
        img_rms1[pixelTypeVec[iDataTypes].flag] = 0.0;
        img_cross[pixelTypeVec[iDataTypes].flag] = 0.0;
      }

      for (int iPx = 0; iPx<nPx; ++iPx) {
        if ( imgOut0->pixel_flag[iPx] != img0->pixel_flag[iPx] ) {
          CERR("image0 has different blanking");
        }
        if ( imgOut1->pixel_flag[iPx] != img1->pixel_flag[iPx] ) {
          CERR("image1 has different blanking");
        }

        if (imgOut0->pixel_flag[iPx] >= PngData::VOLUME_DATA_PIXEL) {
          assert(imgOut1->pixel_flag[iPx] >= PngData::VOLUME_DATA_PIXEL);
          const double val0 = img0->pixel_data[iPx];
          const double val1 = img1->pixel_data[iPx];

          count_buf[iPx] += 1.0; // every pixel gets equal weight for now
          avg0_buf[iPx]  += val0;
          rms0_buf[iPx]   += val0*val0; // for now use rms to store sum(value^2)
          avg1_buf[iPx]  += val1;
          rms1_buf[iPx]   += val1*val1; // for now use rms to store sum(value^2)
          cross_buf[iPx]   += val0*val1; // for now use cross to store sum(value0*value1)

          img_count[img0->pixel_flag[iPx]]       += 1.0;
          img_avg0[img0->pixel_flag[iPx]]       += val0;
          img_rms0[img0->pixel_flag[iPx]]        += val0*val0;
          img_avg1[img1->pixel_flag[iPx]]       += val1;
          img_rms1[img1->pixel_flag[iPx]]        += val1*val1;
          img_cross[img0->pixel_flag[iPx]]       += val0*val1;

        }//(imgOut->pixel_flag[iPx] >= 0)
      }//for (iPx < nPx)
      img0->lock = false;
      img1->lock = false;

      for (int iDataTypes=0;iDataTypes<nDataTypes;++iDataTypes){
        int dataFlag = pixelTypeVec[iDataTypes].flag;
        all_img_count[dataFlag] += img_count[dataFlag];
        all_img_avg0[dataFlag]  += img_avg0[dataFlag];
        all_img_avg1[dataFlag]  += img_avg1[dataFlag];
        all_img_rms0[dataFlag]   += img_rms0[dataFlag];
        all_img_rms1[dataFlag]   += img_rms1[dataFlag];
        all_img_cross[dataFlag]   += img_cross[dataFlag];

        if (img_count[dataFlag] == 0.0) cout << " > image \"" << pngCache->imageList[iI] << "\" has no active " << pixelTypeVec[iDataTypes].name << " data pixels." << endl;
        else {
          img_avg0[dataFlag] /= img_count[dataFlag];
          img_avg1[dataFlag] /= img_count[dataFlag];
          img_rms0[dataFlag] = sqrt(max(0.0,img_rms0[dataFlag]/img_count[dataFlag]-img_avg0[dataFlag]*img_avg0[dataFlag]));
          img_rms1[dataFlag] = sqrt(max(0.0,img_rms1[dataFlag]/img_count[dataFlag]-img_avg1[dataFlag]*img_avg1[dataFlag]));
          img_cross[dataFlag] = img_cross[dataFlag]/img_count[dataFlag]-img_avg0[dataFlag]*img_avg1[dataFlag];
          cout << " > image \"" << pngCache->imageList[iI] << "\" - " << pixelTypeVec[iDataTypes].name << " data - avg0: " << img_avg0[dataFlag] 
                                                                                                       << " avg1: " << img_avg1[dataFlag]
                                                                                                       << " rms0: " << img_rms0[dataFlag]
                                                                                                       << " rms1: " << img_rms1[dataFlag] 
                                                                                                       << " cross: " << img_cross[dataFlag] << endl;
        }
      }

    }//for (iI < iI<pngCache->size())

    for (int iDataTypes=0;iDataTypes<nDataTypes;++iDataTypes){
      int dataFlag = pixelTypeVec[iDataTypes].flag;

      if (all_img_count[dataFlag] == 0.0) cout << " > all images have no active " << pixelTypeVec[iDataTypes].name << " data pixels." << endl;
      else {
        all_img_avg0[dataFlag] /= all_img_count[dataFlag];
        all_img_avg1[dataFlag] /= all_img_count[dataFlag];
        all_img_rms0[dataFlag] = sqrt(max(0.0,all_img_rms0[dataFlag]/all_img_count[dataFlag]-all_img_avg0[dataFlag]*all_img_avg0[dataFlag]));
        all_img_rms1[dataFlag] = sqrt(max(0.0,all_img_rms1[dataFlag]/all_img_count[dataFlag]-all_img_avg1[dataFlag]*all_img_avg1[dataFlag]));
        all_img_cross[dataFlag] = all_img_cross[dataFlag]/all_img_count[dataFlag]-all_img_avg0[dataFlag]*all_img_avg1[dataFlag];
        cout << " > all images - " << pixelTypeVec[iDataTypes].name << " data - avg0: " << all_img_avg0[dataFlag] 
                                                                    << " avg1: " << all_img_avg1[dataFlag] 
                                                                    << " rms0: " << all_img_rms0[dataFlag] 
                                                                    << " rms1: " << all_img_rms1[dataFlag] 
                                                                    << " cross: " << all_img_cross[dataFlag] << endl;
      }
    }

    for (int iPx = 0; iPx<nPx; ++iPx){
      if ( imgOut0->pixel_flag[iPx] >= PngData::VOLUME_DATA_PIXEL ) {
        const double avg0 = avg0_buf[iPx]/count_buf[iPx];
        const double avg1 = avg1_buf[iPx]/count_buf[iPx];
        const double rms0 = sqrt(max(0.0,rms0_buf[iPx]/count_buf[iPx]-avg0*avg0));
        const double rms1 = sqrt(max(0.0,rms1_buf[iPx]/count_buf[iPx]-avg1*avg1));
        const double cross = cross_buf[iPx]/count_buf[iPx]-avg0*avg1;
        rms0_buf[iPx] = rms0;
        rms1_buf[iPx] = rms1;
        cross_buf[iPx] = cross;

        imgOut0->pixel_data[iPx] = avg0;
        imgOut1->pixel_data[iPx] = avg1;
      }
    }

    string avg0_str = outputPathAndPrefix0+".avg0.png";
    cout << " > writing avg0 to \"" << avg0_str << "\"..." << endl;
    string baseVarId0 = imgOut0->varId;
    string baseVarIdSurface0 = imgOut0->varIdSurface;
    string baseVarIdParticles0 = imgOut0->varIdParticles;
    string baseVarIdIso0 = imgOut0->varIdIso;
    imgOut0->varId = baseVarId0+"_avg";
    imgOut0->varIdSurface = baseVarIdSurface0+"_avg";
    imgOut0->varIdParticles = baseVarIdParticles0+"_avg";
    imgOut0->varIdIso = baseVarIdIso0+"_avg";

    imgOut0->rescaleRangesToData();
    imgOut0->initializeWrite();
    imgOut0->write(avg0_str.c_str());
    delete imgOut0;

    //re-init imgOut0; initializeWrite() changes the rgb and data chunks used
    //to construct the lighting estimate.  To get the same lighting estimate
    //for avg, rms, min and max images must grab the original buffers
    imgOut0 = new PngData();
    imgOut0->read(pngCache->imageList[0].c_str());
    imgOut0->finalizeRead();

    for (int iPx = 0; iPx<nPx; ++iPx){
      if ( imgOut0->pixel_flag[iPx] >= PngData::VOLUME_DATA_PIXEL ) {
        imgOut0->pixel_data[iPx] = rms0_buf[iPx];
      }
    }
    string rms0_str = outputPathAndPrefix0+".rms0.png";
    cout << " > writing rms0 to \"" << rms0_str << "\"..." << endl;

    imgOut0->varId = baseVarId0+"_rms";
    imgOut0->varIdSurface = baseVarIdSurface0+"_rms";
    imgOut0->varIdParticles = baseVarIdParticles0+"_rms";
    imgOut0->varIdIso = baseVarIdIso0+"_rms";

    imgOut0->rescaleRangesToData();
    imgOut0->initializeWrite();
    imgOut0->write(rms0_str.c_str());
    delete imgOut0;

    string avg1_str = outputPathAndPrefix1+".avg1.png";
    cout << " > writing avg1 to \"" << avg1_str << "\"..." << endl;
    string baseVarId1 = imgOut1->varId;
    string baseVarIdSurface1 = imgOut1->varIdSurface;
    string baseVarIdParticles1 = imgOut1->varIdParticles;
    string baseVarIdIso1 = imgOut1->varIdIso;
    imgOut1->varId = baseVarId1+"_avg";
    imgOut1->varIdSurface = baseVarIdSurface1+"_avg";
    imgOut1->varIdParticles = baseVarIdParticles1+"_avg";
    imgOut1->varIdIso = baseVarIdIso1+"_avg";

    imgOut1->rescaleRangesToData();
    imgOut1->initializeWrite();
    imgOut1->write(avg1_str.c_str());
    delete imgOut1;

    //re-init imgOut1; initializeWrite() changes the rgb and data chunks used
    //to construct the lighting estimate.  To get the same lighting estimate
    //for avg, rms, min and max images must grab the original buffers
    imgOut1 = new PngData();
    imgOut1->read(pngCache->imageList[0].c_str());
    imgOut1->finalizeRead();

    for (int iPx = 0; iPx<nPx; ++iPx){
      if ( imgOut1->pixel_flag[iPx] >= PngData::VOLUME_DATA_PIXEL ) {
        imgOut1->pixel_data[iPx] = rms1_buf[iPx];
      }
    }
    string rms1_str = outputPathAndPrefix1+".rms1.png";
    cout << " > writing rms1 to \"" << rms1_str << "\"..." << endl;

    imgOut1->varId = baseVarId1+"_rms";
    imgOut1->varIdSurface = baseVarIdSurface1+"_rms";
    imgOut1->varIdParticles = baseVarIdParticles1+"_rms";
    imgOut1->varIdIso = baseVarIdIso1+"_rms";

    imgOut1->rescaleRangesToData();
    imgOut1->initializeWrite();
    imgOut1->write(rms1_str.c_str());
    delete imgOut1;

    //re-init imgOut1; initializeWrite() changes the rgb and data chunks used
    //to construct the lighting estimate.  To get the same lighting estimate
    //for avg, rms, crossimages must grab the original buffers
    imgOut1 = new PngData();
    imgOut1->read(pngCache->imageList[0].c_str());
    imgOut1->finalizeRead();

    for (int iPx = 0; iPx<nPx; ++iPx){
      if ( imgOut1->pixel_flag[iPx] >= PngData::VOLUME_DATA_PIXEL ) {
        imgOut1->pixel_data[iPx] = cross_buf[iPx];
      }
    }
    string cross_str = crossPathAndPrefixSS.str()+".cross.png";
    cout << " > writing cross " << " to \"" << cross_str << "\"..." << endl;

    stringstream ss;
    ss << baseVarId0 << "_" << baseVarId1 << "_cross";
    imgOut1->varId = ss.str();
    ss.str("");
    ss << baseVarIdSurface0 << "_" << baseVarIdSurface1 << "_cross";
    imgOut1->varIdSurface = ss.str();
    ss.str("");
    ss << baseVarIdParticles0 << "_" << baseVarIdParticles1 << "_cross";
    imgOut1->varIdParticles = ss.str();
    ss.str("");
    ss << baseVarIdIso0 << "_" << baseVarIdIso1 << "_cross";
    imgOut1->varIdIso = ss.str();
    ss.str("");

    imgOut1->rescaleRangesToData();
    imgOut1->initializeWrite();
    imgOut1->write(cross_str.c_str());
    delete imgOut1;



    delete[] count_buf;
    delete[] avg0_buf;
    delete[] avg1_buf;
    delete[] rms0_buf;
    delete[] rms1_buf;
    delete[] cross_buf;

  }//itpTemporalStatsCross()



  void itpTemporalStatsMoving(Param * param){

    COUT1("TEMPORAL_STATS_MOVING...");

    if ((!pngCache)||(pngCache->imageList.size() == 0)) {
      CERR("cannot compute TEMPORAL_STATS_MOVING: pngCache is empty.");
    }

    string outputPath = "images_stats/";

    int iarg = 0;
    while (iarg < param->size()) {
      const string token = param->getString(iarg++);
      // assume this string is a output path...
      outputPath = token+"/";
      cout << " > Assuming output path \"" << outputPath << endl;
    }

    string outputPathAndPrefix = outputPath + getPrefix(pngCache->imageList[0]);
    COUT1(" > output file path and prefix will be \"" << outputPathAndPrefix << "\"");
    MiscUtils::mkdir_for_file(outputPath+"/test.png");

    PngData * imgOut = pngCache->getImage(0);
    int nx = imgOut->getNx();
    int ny = imgOut->getNy();
    int nPx = nx*ny;
    imgOut->lock = false;

    map<int,double*> depth_buf;
    map<int,double*> count_buf;
    map<int,double*> avg_buf;
    map<int,double*> rms_buf;

    map<int,double> all_img_count;
    map<int,double> all_img_avg;
    map<int,double> all_img_rms;
    vector<pixel_type> pixelTypeVecGlobal;

    map<int,double> img_count;
    map<int,double> img_avg;
    map<int,double> img_rms;

    for (int iI = 0; iI<pngCache->size(); ++iI){

      COUT1(" > working on image " << iI << " of " << pngCache->size() << " \"" << pngCache->imageList[iI] << "\"...");

      PngData * img = pngCache->getImage(iI);
      assert(img->getNx() == nx); // if you hit one of these, images are different sizes -- need to handle in some future version
      assert(img->getNy() == ny);

      vector<pixel_type> pixelTypeVec;
      img->getPixelTypes(pixelTypeVec);

      img_count.clear();
      img_avg.clear();
      img_rms.clear();
      for (int ipt=0;ipt<pixelTypeVec.size();++ipt){
        img_count[pixelTypeVec[ipt].flag] = 0.0;
        img_avg[pixelTypeVec[ipt].flag] = 0.0;
        img_rms[pixelTypeVec[ipt].flag] = 0.0;
        if (count_buf.find(pixelTypeVec[ipt].flag)==count_buf.end()){
          //first time we have encountered this data type, allocate buffer
          depth_buf[pixelTypeVec[ipt].flag] = new double[nPx];
          count_buf[pixelTypeVec[ipt].flag] = new double[nPx];
          avg_buf[pixelTypeVec[ipt].flag] = new double[nPx];
          rms_buf[pixelTypeVec[ipt].flag] = new double[nPx];
          for (int iPx = 0; iPx<nPx; ++iPx){
            depth_buf[pixelTypeVec[ipt].flag][iPx] = 0.0;
            count_buf[pixelTypeVec[ipt].flag][iPx] = 0.0;
            avg_buf[pixelTypeVec[ipt].flag][iPx] = 0.0;
            rms_buf[pixelTypeVec[ipt].flag][iPx] = 0.0;
          }
          pixelTypeVecGlobal.push_back(pixelTypeVec[ipt]);
        }
        else{
          //this data type has appeared in other images, check
          //that we are working with the same variable name
          int iptg;
          for (iptg=0; iptg<pixelTypeVecGlobal.size(); ++iptg){
            if (pixelTypeVecGlobal[iptg].flag == pixelTypeVec[ipt].flag){
              if (pixelTypeVecGlobal[iptg].varId !=  pixelTypeVec[ipt].varId){
                CERR("Images contain " << pixelTypeVecGlobal[iptg].name << " data with different variables: " << pixelTypeVecGlobal[iptg].varId << " and " << pixelTypeVec[ipt].varId);
              }
              else{
                break;
              }
            }
          }
          assert(iptg<pixelTypeVecGlobal.size()); //global vec should have a matching data type flag
        }
      }

      for (int iPx = 0; iPx<nPx; ++iPx) {
        if (img->pixel_flag[iPx] >= PngData::VOLUME_DATA_PIXEL) {
          const double val = img->pixel_data[iPx];

          depth_buf[img->pixel_flag[iPx]][iPx] += double(img->getDepth(iPx));
          count_buf[img->pixel_flag[iPx]][iPx] += 1.0; // every pixel gets equal weight for now
          avg_buf[img->pixel_flag[iPx]][iPx]  += val;
          rms_buf[img->pixel_flag[iPx]][iPx]   += val*val; // for now use rms to store sum(value^2)

          img_count[img->pixel_flag[iPx]]      += 1.0;
          img_avg[img->pixel_flag[iPx]]       += val;
          img_rms[img->pixel_flag[iPx]]        += val*val;

        }//(img->pixel_flag[iPx] >= PngData::VOLUME_DATA_PIXEL)
      }//for (iPx < nPx)
      img->lock = false;

      for (int ipt=0; ipt<pixelTypeVec.size(); ++ipt){
        const int dataFlag = pixelTypeVec[ipt].flag;
        if (all_img_count.find(dataFlag)!=all_img_count.end()){
          all_img_count[dataFlag] += img_count[dataFlag];
          all_img_avg[dataFlag]  += img_avg[dataFlag];
          all_img_rms[dataFlag]   += img_rms[dataFlag];
        }
        else{ //first time we found this data type in an image...
          all_img_count[dataFlag] = img_count[dataFlag];
          all_img_avg[dataFlag]  = img_avg[dataFlag];
          all_img_rms[dataFlag]   = img_rms[dataFlag];
        }

        assert(img_count[dataFlag]>0); //pixelTypeVec should only include entries when data pixels are present...
        img_avg[dataFlag] /= img_count[dataFlag];
        img_rms[dataFlag] = sqrt(max(0.0,img_rms[dataFlag]/img_count[dataFlag]-img_avg[dataFlag]*img_avg[dataFlag]));
        cout << " > image \"" << pngCache->imageList[iI] << "\" - " << pixelTypeVec[ipt].name << " data - avg: " << img_avg[dataFlag] << " rms: " << img_rms[dataFlag] << endl;
      }

    }//for (iI < iI<pngCache->size())

    ImageMetadata * imd0 = pngCache->getMetadata(0);

    for (int ipt=0; ipt<pixelTypeVecGlobal.size(); ++ipt){
      const int dataFlag = pixelTypeVecGlobal[ipt].flag;

      if (all_img_count[dataFlag] == 0.0) {
        cout << " > all images have no active " <<  pixelTypeVecGlobal[ipt].name << " data pixels." << endl;
        assert(0); //should not get here...
      }
      else {
        all_img_avg[dataFlag] /= all_img_count[dataFlag];
        all_img_rms[dataFlag] = sqrt(max(0.0,all_img_rms[dataFlag]/all_img_count[dataFlag]-all_img_avg[dataFlag]*all_img_avg[dataFlag]));
        cout << " > all images - " << pixelTypeVecGlobal[ipt].name << " data - avg: " << all_img_avg[dataFlag] << " rms: " << all_img_rms[dataFlag] << endl;
      }

      imgOut = new PngData(nx,ny);

      for (int iPx = 0; iPx<nPx; ++iPx){
        if (count_buf[dataFlag][iPx] > 0) {
          imgOut->pixel_flag[iPx] = dataFlag;
          const double avg = avg_buf[dataFlag][iPx]/count_buf[dataFlag][iPx];
          const double rms = sqrt(max(0.0,rms_buf[dataFlag][iPx]/count_buf[dataFlag][iPx]-avg*avg));
          rms_buf[dataFlag][iPx] = rms;
          imgOut->pixel_data[iPx] = avg;
          //set average depth
          const double depth_avg = depth_buf[dataFlag][iPx]/count_buf[dataFlag][iPx];
          imgOut->setDepth(iPx,uint2 (depth_avg));
        }
        //else no data at this pixel, imgOut already set to background
      }

      //build metadata
      ImageMetadata imdOut;
      imdOut.transformMat[0]  = imd0->transformMat[0];
      imdOut.transformMat[1]  = imd0->transformMat[1];
      imdOut.transformMat[2]  = imd0->transformMat[2];
      imdOut.transformMat[3]  = imd0->transformMat[3];
      imdOut.transformMat[4]  = imd0->transformMat[4];
      imdOut.transformMat[5]  = imd0->transformMat[5];
      imdOut.transformMat[6]  = imd0->transformMat[6];
      imdOut.transformMat[7]  = imd0->transformMat[7];
      imdOut.transformMat[8]  = imd0->transformMat[8];
      imdOut.transformMat[9]  = imd0->transformMat[9];
      imdOut.transformMat[10] = imd0->transformMat[10];
      imdOut.transformMat[11] = imd0->transformMat[11];
      imdOut.transformMat[12] = imd0->transformMat[12];
      imdOut.transformMat[13] = imd0->transformMat[13];
      imdOut.transformMat[14] = imd0->transformMat[14];
      imdOut.transformMat[15] = imd0->transformMat[15];

      imdOut.setLengthScale(imd0->getLengthScale());
      imdOut.setRgbMode(imd0->getRgbMode());
      imdOut.setCamDepth(imd0->getCamDepth());
      imdOut.setTime(0.0);
      imdOut.setWriteImageParam(param->str());

      //reset image data fields
      imgOut->b_data = false;
      imgOut->b_dataSurface = false;
      imgOut->b_dataParticles = false;
      imgOut->b_dataIso = false;

      //avg image...
      if (pixelTypeVecGlobal[ipt].flag==PngData::VOLUME_DATA_PIXEL){
        imgOut->b_data = true;
        imgOut->setColorMap(imd0->getColorMapName());
        imdOut.setVarId(pixelTypeVecGlobal[ipt].varId + "_avg");
      }
      else if (pixelTypeVecGlobal[ipt].flag==PngData::SURFACE_DATA_PIXEL){
        imgOut->b_dataSurface = true;
        imgOut->setColorMapSurface(imd0->getSurfColorMapName());
        imdOut.setVarOnSurfaceId(pixelTypeVecGlobal[ipt].varId + "_avg");
      }
      else if (pixelTypeVecGlobal[ipt].flag==PngData::PARTICLE_DATA_PIXEL){
        imgOut->b_dataParticles = true;
        imgOut->setColorMapParticles(imd0->getPartColorMapName());
        imdOut.setVarOnParticleId(pixelTypeVecGlobal[ipt].varId + "_avg");
      }
      else if (pixelTypeVecGlobal[ipt].flag==PngData::ISO_DATA_PIXEL){
        imgOut->b_dataIso = true;
        imgOut->setColorMapIso(imd0->getIsoColorMapName());
        imdOut.setVarOnIsoId(pixelTypeVecGlobal[ipt].varId + "_avg");
      }

      char filename[128];
      sprintf(filename,"%s.%s.avg.png",outputPathAndPrefix.c_str(),pixelTypeVecGlobal[ipt].name.c_str());
      cout << " > writing avg (datatype " << pixelTypeVecGlobal[ipt].name << ") to \"" << filename << "\"..." << endl;

      imgOut->rescaleRangesToData();
      imgOut->initializeWrite(imdOut);
      imgOut->write(filename);


      //rms image...
      if (pixelTypeVecGlobal[ipt].flag==PngData::VOLUME_DATA_PIXEL){
        imdOut.setVarId(pixelTypeVecGlobal[ipt].varId + "_rms");
      }
      else if (pixelTypeVecGlobal[ipt].flag==PngData::SURFACE_DATA_PIXEL){
        imdOut.setVarOnSurfaceId(pixelTypeVecGlobal[ipt].varId + "_rms");
      }
      else if (pixelTypeVecGlobal[ipt].flag==PngData::PARTICLE_DATA_PIXEL){
        imdOut.setVarOnParticleId(pixelTypeVecGlobal[ipt].varId + "_rms");
      }
      else if (pixelTypeVecGlobal[ipt].flag==PngData::ISO_DATA_PIXEL){
        imdOut.setVarOnIsoId(pixelTypeVecGlobal[ipt].varId + "_rms");
      }

      char filename_rms[128];
      sprintf(filename_rms,"%s.%s.rms.png",outputPathAndPrefix.c_str(),pixelTypeVecGlobal[ipt].name.c_str());

      cout << " > writing rms (datatype " << pixelTypeVecGlobal[ipt].name << ") to \"" << filename_rms << "\"..." << endl;

      for (int iPx = 0; iPx<nPx; ++iPx){
        if (count_buf[dataFlag][iPx] > 0) {
          imgOut->pixel_data[iPx] = rms_buf[dataFlag][iPx];
        }
      }
      imgOut->rescaleRangesToData();
      imgOut->initializeWrite(imdOut);
      imgOut->write(filename_rms);

      // count image...
      if (pixelTypeVecGlobal[ipt].flag==PngData::VOLUME_DATA_PIXEL){
        imdOut.setVarId(pixelTypeVecGlobal[ipt].varId + "_count");
      }
      else if (pixelTypeVecGlobal[ipt].flag==PngData::SURFACE_DATA_PIXEL){
        imdOut.setVarOnSurfaceId(pixelTypeVecGlobal[ipt].varId + "_count");
      }
      else if (pixelTypeVecGlobal[ipt].flag==PngData::PARTICLE_DATA_PIXEL){
        imdOut.setVarOnParticleId(pixelTypeVecGlobal[ipt].varId + "_count");
      }
      else if (pixelTypeVecGlobal[ipt].flag==PngData::ISO_DATA_PIXEL){
        imdOut.setVarOnIsoId(pixelTypeVecGlobal[ipt].varId + "_count");
      }

      char filename_count[128];
      sprintf(filename_count,"%s.%s.count.png",outputPathAndPrefix.c_str(),pixelTypeVecGlobal[ipt].name.c_str());

      cout << " > writing count (datatype " << pixelTypeVecGlobal[ipt].name << ") to \"" << filename_count << "\"..." << endl;

      for (int iPx = 0; iPx<nPx; ++iPx){
        if (count_buf[dataFlag][iPx] > 0) {
          imgOut->pixel_data[iPx] = count_buf[dataFlag][iPx];
        }
      }
      imgOut->rescaleRangesToData();
      imgOut->initializeWrite(imdOut);
      imgOut->write(filename_count);
      
      delete imgOut;
      imgOut = NULL;

      //clean up memory for each data type...
      delete[] depth_buf[dataFlag];
      delete[] count_buf[dataFlag];
      delete[] avg_buf[dataFlag];
      delete[] rms_buf[dataFlag];
    }

    delete imgOut;
  }//itpTemporalStatsMoving()

  ///////////////////////////////////////////////////////////////////////////

  void itpTemporalStatsPhase(Param * param){

    COUT1("TEMPORAL_STATS_PHASE...");

    if ((!pngCache)||(pngCache->imageList.size() == 0)) {
      CERR("cannot compute TEMPORAL_STATS_PHASE: pngCache is empty.");
    }

    string outputPath = "images_stats/";
    bool b_period = false;
    double period = -1;
    bool b_nbin = false;
    int nbin = 1;
    bool b_imgdt= false;
    double imgdt = -1.0;

    int iarg = 0;
    while (iarg < param->size()) {
      const string token = param->getString(iarg++);
      if (token == "PERIOD") {
        b_period = true;
        period = param->getDouble(iarg++);
	cout << " > got PERIOD " << period << endl;
      }
      else if ((token == "NBIN")||(token == "NBINS")||(token == "BINS")) {
        b_nbin = true;
        nbin = param->getInt(iarg++);
	cout << " > got NBIN " << nbin << endl;
      }
      else if (token == "IMGDT") {
        b_imgdt = true;
        imgdt = param->getDouble(iarg++);
	cout << " > got IMGDT " << imgdt << endl;
      }
      else {
	// assume this string is a output path...
	outputPath = token+"/";
	cout << " > Assuming output path \"" << outputPath << endl;
      }
    }

    string outputPathAndPrefix = outputPath + getPrefix(pngCache->imageList[0]);
    COUT1(" > output file path and prefix will be \"" << outputPathAndPrefix << "\"");
    MiscUtils::mkdir_for_file(outputPath+"/test.png");

    if (!b_period || !b_nbin) {
      CERR("Phase-averaging images requires PERIOD=<double> NBIN=<int>");
    }

    PngData * imgOut = pngCache->getImage(0);
    int nx = imgOut->getNx();
    int ny = imgOut->getNy();
    int nPx = nx*ny;
    imgOut->lock = false;

    map<int,double*> depth_buf;
    map<int,double*> count_buf;
    map<int,double*> avg_buf;
    map<int,double*> rms_buf;

    map<int,double> all_img_count;
    map<int,double> all_img_avg;
    map<int,double> all_img_rms;
    vector<pixel_type> pixelTypeVecGlobal;

    map<int,double> img_count;
    map<int,double> img_avg;
    map<int,double> img_rms;

    for (int iI = 0; iI<pngCache->size(); ++iI){

      COUT1(" > working on image \"" << pngCache->imageList[iI] << "\"");

      PngData * img = pngCache->getImage(iI);
      assert(img->getNx() == nx); // if you hit one of these, images are different sizes -- need to handle in some future version
      assert(img->getNy() == ny);

      vector<pixel_type> pixelTypeVec;
      img->getPixelTypes(pixelTypeVec);

      img_count.clear();
      img_avg.clear();
      img_rms.clear();
      for (int ipt=0;ipt<pixelTypeVec.size();++ipt){
        img_count[pixelTypeVec[ipt].flag] = 0.0;
        img_avg[pixelTypeVec[ipt].flag] = 0.0;
        img_rms[pixelTypeVec[ipt].flag] = 0.0;
        if (count_buf.find(pixelTypeVec[ipt].flag)==count_buf.end()){
          //first time we have encountered this data type, allocate buffer
          depth_buf[pixelTypeVec[ipt].flag] = new double[nPx*nbin];
          count_buf[pixelTypeVec[ipt].flag] = new double[nPx*nbin];
          avg_buf[pixelTypeVec[ipt].flag] = new double[nPx*nbin];
          rms_buf[pixelTypeVec[ipt].flag] = new double[nPx*nbin];
          for (int iPx = 0; iPx<nPx*nbin; ++iPx){
            depth_buf[pixelTypeVec[ipt].flag][iPx] = 0.0;
            count_buf[pixelTypeVec[ipt].flag][iPx] = 0.0;
            avg_buf[pixelTypeVec[ipt].flag][iPx] = 0.0;
            rms_buf[pixelTypeVec[ipt].flag][iPx] = 0.0;
          }
          pixelTypeVecGlobal.push_back(pixelTypeVec[ipt]);
        }
        else{
          //this data type has appeared in other images, check
          //that we are working with the same variable name
          int iptg;
          for (iptg=0; iptg<pixelTypeVecGlobal.size(); ++iptg){
            if (pixelTypeVecGlobal[iptg].flag == pixelTypeVec[ipt].flag){
              if (pixelTypeVecGlobal[iptg].varId !=  pixelTypeVec[ipt].varId){
                CERR("Images contain " << pixelTypeVecGlobal[iptg].name << " data with different variables: " << pixelTypeVecGlobal[iptg].varId << " and " << pixelTypeVec[ipt].varId);
              }
              else{
                break;
              }
            }
          }
          assert(iptg<pixelTypeVecGlobal.size()); //global vec should have a matching data type flag
        }
      }

      double time;
      if (b_imgdt)
        time = imgdt*iI;
      else
        time = img->getTime();

      while (time >= period)
        time -= period;
      const int ibin = (int)floor(time/period*double(nbin));
      assert((ibin >= 0)&&(ibin < nbin));
      cout << " > got time: " << time << " time/period*nbin: " << time/period*double(nbin) << endl;

      for (int iPx = 0; iPx<nPx; ++iPx) {
        if (img->pixel_flag[iPx] >= PngData::VOLUME_DATA_PIXEL) {
          const double val = img->pixel_data[iPx];

          depth_buf[img->pixel_flag[iPx]][iPx+nPx*ibin] += double(img->getDepth(iPx));
          count_buf[img->pixel_flag[iPx]][iPx+nPx*ibin] += 1.0; // every pixel gets equal weight for now
          avg_buf[img->pixel_flag[iPx]][iPx+nPx*ibin]  += val;
          rms_buf[img->pixel_flag[iPx]][iPx+nPx*ibin]   += val*val; // for now use rms to store sum(value^2)

          img_count[img->pixel_flag[iPx]]      += 1.0;
          img_avg[img->pixel_flag[iPx]]       += val;
          img_rms[img->pixel_flag[iPx]]        += val*val;

        }//(img->pixel_flag[iPx] >= PngData::VOLUME_DATA_PIXEL)
      }//for (iPx < nPx)
      img->lock = false;

      for (int ipt=0; ipt<pixelTypeVec.size(); ++ipt){
        const int dataFlag = pixelTypeVec[ipt].flag;
        if (all_img_count.find(dataFlag)!=all_img_count.end()){
          all_img_count[dataFlag] += img_count[dataFlag];
          all_img_avg[dataFlag]  += img_avg[dataFlag];
          all_img_rms[dataFlag]   += img_rms[dataFlag];
        }
        else{ //first time we found this data type in an image...
          all_img_count[dataFlag] = img_count[dataFlag];
          all_img_avg[dataFlag]  = img_avg[dataFlag];
          all_img_rms[dataFlag]   = img_rms[dataFlag];
        }

        assert(img_count[dataFlag]>0); //pixelTypeVec should only include entries when data pixels are present...
        img_avg[dataFlag] /= img_count[dataFlag];
        img_rms[dataFlag] = sqrt(max(0.0,img_rms[dataFlag]/img_count[dataFlag]-img_avg[dataFlag]*img_avg[dataFlag]));
        cout << " > image \"" << pngCache->imageList[iI] << "\" - " << pixelTypeVec[ipt].name << " data - avg: " << img_avg[dataFlag] << " rms: " << img_rms[dataFlag] << endl;
      }

    }//for (iI < iI<pngCache->size())

    ImageMetadata * imd0 = pngCache->getMetadata(0);

    for (int ipt=0; ipt<pixelTypeVecGlobal.size(); ++ipt){
      const int dataFlag = pixelTypeVecGlobal[ipt].flag;

      if (all_img_count[dataFlag] == 0.0) {
        cout << " > all images have no active " <<  pixelTypeVecGlobal[ipt].name << " data pixels." << endl;
        assert(0); //should not get here...
      }
      else {
        all_img_avg[dataFlag] /= all_img_count[dataFlag];
        all_img_rms[dataFlag] = sqrt(max(0.0,all_img_rms[dataFlag]/all_img_count[dataFlag]-all_img_avg[dataFlag]*all_img_avg[dataFlag]));
        cout << " > all images - " << pixelTypeVecGlobal[ipt].name << " data - avg: " << all_img_avg[dataFlag] << " rms: " << all_img_rms[dataFlag] << endl;
      }

      for (int ibin = 0; ibin < nbin; ++ibin) {

        imgOut = new PngData(nx,ny);

        for (int iPx = 0; iPx<nPx; ++iPx){
          if (count_buf[dataFlag][iPx+nPx*ibin] > 0) {
            imgOut->pixel_flag[iPx] = dataFlag;
            const double avg = avg_buf[dataFlag][iPx+nPx*ibin]/count_buf[dataFlag][iPx+nPx*ibin];
            const double rms = sqrt(max(0.0,rms_buf[dataFlag][iPx+nPx*ibin]/count_buf[dataFlag][iPx+nPx*ibin]-avg*avg));
            rms_buf[dataFlag][iPx+nPx*ibin] = rms;
            imgOut->pixel_data[iPx] = avg;
            //set average depth
            const double depth_avg = depth_buf[dataFlag][iPx+nPx*ibin]/count_buf[dataFlag][iPx+nPx*ibin];
            imgOut->setDepth(iPx,uint2 (depth_avg));
          }
          //else no data at this pixel, imgOut already set to background
        }

        //build metadata
        ImageMetadata imdOut;
        imdOut.transformMat[0]  = imd0->transformMat[0];
        imdOut.transformMat[1]  = imd0->transformMat[1];
        imdOut.transformMat[2]  = imd0->transformMat[2];
        imdOut.transformMat[3]  = imd0->transformMat[3];
        imdOut.transformMat[4]  = imd0->transformMat[4];
        imdOut.transformMat[5]  = imd0->transformMat[5];
        imdOut.transformMat[6]  = imd0->transformMat[6];
        imdOut.transformMat[7]  = imd0->transformMat[7];
        imdOut.transformMat[8]  = imd0->transformMat[8];
        imdOut.transformMat[9]  = imd0->transformMat[9];
        imdOut.transformMat[10] = imd0->transformMat[10];
        imdOut.transformMat[11] = imd0->transformMat[11];
        imdOut.transformMat[12] = imd0->transformMat[12];
        imdOut.transformMat[13] = imd0->transformMat[13];
        imdOut.transformMat[14] = imd0->transformMat[14];
        imdOut.transformMat[15] = imd0->transformMat[15];

        imdOut.setLengthScale(imd0->getLengthScale());
        imdOut.setRgbMode(imd0->getRgbMode());
        imdOut.setCamDepth(imd0->getCamDepth());
        imdOut.setTime(double(ibin)/double(nbin)*period);
        imdOut.setWriteImageParam(param->str());

        //reset image data fields
        imgOut->b_data = false;
        imgOut->b_dataSurface = false;
        imgOut->b_dataParticles = false;
        imgOut->b_dataIso = false;

        //avg image...
        if (pixelTypeVecGlobal[ipt].flag==PngData::VOLUME_DATA_PIXEL){
          imgOut->b_data = true;
          imgOut->setColorMap(imd0->getColorMapName());
          imdOut.setVarId(pixelTypeVecGlobal[ipt].varId + "_avg");
        }
        else if (pixelTypeVecGlobal[ipt].flag==PngData::SURFACE_DATA_PIXEL){
          imgOut->b_dataSurface = true;
          imgOut->setColorMapSurface(imd0->getSurfColorMapName());
          imdOut.setVarOnSurfaceId(pixelTypeVecGlobal[ipt].varId + "_avg");
        }
        else if (pixelTypeVecGlobal[ipt].flag==PngData::PARTICLE_DATA_PIXEL){
          imgOut->b_dataParticles = true;
          imgOut->setColorMapParticles(imd0->getPartColorMapName());
          imdOut.setVarOnParticleId(pixelTypeVecGlobal[ipt].varId + "_avg");
        }
        else if (pixelTypeVecGlobal[ipt].flag==PngData::ISO_DATA_PIXEL){
          imgOut->b_dataIso = true;
          imgOut->setColorMapIso(imd0->getIsoColorMapName());
          imdOut.setVarOnIsoId(pixelTypeVecGlobal[ipt].varId + "_avg");
        }

        char filename[128];
        sprintf(filename,"%s.%s.avg-%03d.png",outputPathAndPrefix.c_str(),pixelTypeVecGlobal[ipt].name.c_str(),ibin);
        cout << " > writing phase avg (datatype " << pixelTypeVecGlobal[ipt].name << ") to \"" << filename << "\"..." << endl;

        imgOut->rescaleRangesToData();
        imgOut->initializeWrite(imdOut);
        imgOut->write(filename);


        //rms image...
        if (pixelTypeVecGlobal[ipt].flag==PngData::VOLUME_DATA_PIXEL){
          imdOut.setVarId(pixelTypeVecGlobal[ipt].varId + "_rms");
        }
        else if (pixelTypeVecGlobal[ipt].flag==PngData::SURFACE_DATA_PIXEL){
          imdOut.setVarOnSurfaceId(pixelTypeVecGlobal[ipt].varId + "_rms");
        }
        else if (pixelTypeVecGlobal[ipt].flag==PngData::PARTICLE_DATA_PIXEL){
          imdOut.setVarOnParticleId(pixelTypeVecGlobal[ipt].varId + "_rms");
        }
        else if (pixelTypeVecGlobal[ipt].flag==PngData::ISO_DATA_PIXEL){
          imdOut.setVarOnIsoId(pixelTypeVecGlobal[ipt].varId + "_rms");
        }

        char filename_rms[128];
        sprintf(filename_rms,"%s.%s.rms-%03d.png",outputPathAndPrefix.c_str(),pixelTypeVecGlobal[ipt].name.c_str(),ibin);

        cout << " > writing phase rms (datatype " << pixelTypeVecGlobal[ipt].name << ") to \"" << filename_rms << "\"..." << endl;

        for (int iPx = 0; iPx<nPx; ++iPx){
          if (count_buf[dataFlag][iPx+nPx*ibin] > 0) {
            imgOut->pixel_data[iPx] = rms_buf[dataFlag][iPx+nPx*ibin];
          }
        }
        imgOut->rescaleRangesToData();
        imgOut->initializeWrite(imdOut);
        imgOut->write(filename_rms);

        delete imgOut;
        imgOut = NULL;
      }
      //clean up memory for each data type...
      delete[] depth_buf[dataFlag];
      delete[] count_buf[dataFlag];
      delete[] avg_buf[dataFlag];
      delete[] rms_buf[dataFlag];
    }

    delete imgOut;
  }//itpTemporalStatsPhase()

  void itpTemporalStatsPhaseLock(Param * param){

    COUT1("TEMPORAL_STATS_PHASE_LOCK...");

    if ((!pngCache)||(pngCache->imageList.size() == 0)) {
      CERR("cannot compute TEMPORAL_STATS_PHASE_LOCK: pngCache is empty.");
    }

    string outputPath    = "images_stats/";
    bool b_period        = false; double period = 1;
    bool b_nbin          = false; int    nbin   = 1;
    bool b_imgdt         = false; double imgdt  = 1.0;
    bool b_bindex        = false; int    bindex = 0;
    PhaseMethod pmethod  = BIN;   int    nbuf   = 1;
    bool b_lock          = false; string lockstr("nodes.lock"); double time0 = 0.0;

    bool b_range_avg            = false; double range_avg[2];
    bool b_range_rms            = false; double range_rms[2];
    bool b_rangeSurface_avg     = false; double rangeSurface_avg[2];
    bool b_rangeSurface_rms     = false; double rangeSurface_rms[2];
    bool b_rangeParticles_avg   = false; double rangeParticles_avg[2];
    bool b_rangeParticles_rms   = false; double rangeParticles_rms[2];
    bool b_rangeIso_avg         = false; double rangeIso_avg[2];
    bool b_rangeIso_rms         = false; double rangeIso_rms[2];

    int iarg = 0;
    while(iarg < param->size())
      {
        const string token = param->getString(iarg++);
        if(token == "PERIOD")
          {
            b_period = true;
            period = param->getDouble(iarg++);
            cout << " > got PERIOD " << period << endl;
          }
        else if((token == "NBIN")||(token == "NBINS")||(token == "BINS"))
          {
            b_nbin = true;
            nbin = param->getInt(iarg++);
            cout << " > got NBIN " << nbin << endl;
          }
        else if((token == "IMGDT"))
          {
            b_imgdt = true;
            imgdt = param->getDouble(iarg++);
            cout << " > got " << token << " " << imgdt << endl;
          }
        else if((token == "BINDEX"))
          {
            b_bindex = true;
            bindex = param->getInt(iarg++);
            cout << " > got " << token << " " << bindex << endl;
          }
        else if((token == "BIN") || (token == "WEIGHT") || (token == "SPLINE") || (token == "FOURIER"))
          {
            pmethod = BIN;
            nbuf = 1;
          }
        else if((token == "LINEAR") || (token == "LIN"))
          {
            pmethod = LINEAR;
            nbuf = 3; // Local linear interpolation between nearest neighbor bins
          }
        else if((token == "LOCK") || (token == "LOCK_FILE"))
          {
            b_lock = true;
            lockstr = param->getString(iarg++);
          }
        else if((token == "START_TIME"))
          {
            time0 = param->getDouble(iarg++);
            cout << " > got START_TIME " << time0 << endl;
          }
        else if((token == "RANGE_AVG")||(token == "AVG_RANGE")||(token == "RANGE_AVG")||(token == "AVG_RANGE"))
          {
            b_range_avg = true;
            range_avg[0] = param->getDouble(iarg++);
            range_avg[1] = param->getDouble(iarg++);
            cout << " > got RANGE_AVG " << range_avg[0] << " " << range_avg[1] << endl;
          }
        else if((token == "RANGE_RMS")||(token == "RMS_RANGE"))
          {
            b_range_rms = true;
            range_rms[0] = param->getDouble(iarg++);
            range_rms[1] = param->getDouble(iarg++);
            cout << " > got RANGE_RMS " << range_rms[0] << " " << range_rms[1] << endl;
          }
        else if((token == "RANGE_SURFACE_AVG")||(token == "AVG_RANGE_SURFACE")||(token == "RANGE_SURFACE_AVG")||(token == "AVG_RANGE_SURFACE"))
          {
            b_rangeSurface_avg = true;
            rangeSurface_avg[0] = param->getDouble(iarg++);
            rangeSurface_avg[1] = param->getDouble(iarg++);
            cout << " > got RANGE_SURFACE_AVG " << rangeSurface_avg[0] << " " << rangeSurface_avg[1] << endl;
          }
        else if((token == "RANGE_SURFACE_RMS")||(token == "RMS_RANGE_SURFACE"))
          {
            b_rangeSurface_rms = true;
            rangeSurface_rms[0] = param->getDouble(iarg++);
            rangeSurface_rms[1] = param->getDouble(iarg++);
            cout << " > got RANGE_SURFACE_RMS " << rangeSurface_rms[0] << " " << rangeSurface_rms[1] << endl;
          }
        else if((token == "RANGE_PARTICLES_AVG")||(token == "AVG_RANGE_PARTICLES")||(token == "RANGE_PARTICLES_AVG")||(token == "AVG_RANGE_PARTICLES"))
          {
            b_rangeParticles_avg = true;
            rangeParticles_avg[0] = param->getDouble(iarg++);
            rangeParticles_avg[1] = param->getDouble(iarg++);
            cout << " > got RANGE_PARTICLES_AVG " << rangeParticles_avg[0] << " " << rangeParticles_avg[1] << endl;
          }
        else if((token == "RANGE_PARTICLES_RMS")||(token == "RMS_RANGE_PARTICLES"))
          {
            b_rangeParticles_rms = true;
            rangeParticles_rms[0] = param->getDouble(iarg++);
            rangeParticles_rms[1] = param->getDouble(iarg++);
            cout << " > got RANGE_PARTICLES_RMS " << rangeParticles_rms[0] << " " << rangeParticles_rms[1] << endl;
          }
        else if((token == "RANGE_ISO_AVG")||(token == "AVG_RANGE_ISO")||(token == "RANGE_ISO_AVG")||(token == "AVG_RANGE_ISO"))
          {
            b_rangeIso_avg = true;
            rangeIso_avg[0] = param->getDouble(iarg++);
            rangeIso_avg[1] = param->getDouble(iarg++);
            cout << " > got RANGE_ISO_AVG " << rangeIso_avg[0] << " " << rangeIso_avg[1] << endl;
          }
        else if((token == "RANGE_ISO_RMS")||(token == "RMS_RANGE_ISO"))
          {
            b_rangeIso_rms = true;
            rangeIso_rms[0] = param->getDouble(iarg++);
            rangeIso_rms[1] = param->getDouble(iarg++);
            cout << " > got RANGE_ISO_RMS " << rangeIso_rms[0] << " " << rangeIso_rms[1] << endl;
          }
        else
          {
            // assume this string is a output path...
  	    outputPath = token+"/";
            cout << " > Assuming output path \"" << outputPath << endl;
          }
      }

    string outputPathAndPrefix = outputPath + getPrefix(pngCache->imageList[0]);
    COUT1(" > output file path and prefix will be \"" << outputPathAndPrefix << "\"");
    MiscUtils::mkdir_for_file(outputPath+"/test.png");

    if(!b_imgdt){
      COUT1(" > time increment between images will be read from image metadata");
    }
    else{
      COUT1(" > time increment between images will be \"" << imgdt << "\"");
    }

    if(!b_period || !b_nbin)
      { CERR("Phase-averaging images requires PERIOD=<double> NBIN=<int>"); }

    vector<int> set_bindex;
    if(b_bindex == false)
      {
        COUT1(" > no bin index (BINDEX) provided; assuming all bins are binned");
        set_bindex.resize(nbin); for(int i=0; i < nbin; ++i) set_bindex[i] = i;
      }
    else if(bindex < 0 || bindex >= nbin)
      { CERR("Invalid BINDEX=" << bindex << " specified; change so that 0 <= BINDEX < NBIN"); }
    else
      { set_bindex.resize(1); set_bindex[0] = bindex; }

    // If no BINDEX set, then process all BINDEX values; otherwise, just process the given BINDEX
    for(std::vector<int>::iterator mybindexit = set_bindex.begin(); mybindexit != set_bindex.end(); ++mybindexit)
      {// set_bindex loop begin
        bindex = *mybindexit;

        PngData * imgOut = pngCache->getImage(0);
        int nx = imgOut->getNx();
        int ny = imgOut->getNy();
        int nPx = nx*ny;
        imgOut->lock = false;

        map<int,double*> depth_buf;
        map<int,double*> count_buf;
        map<int,double*> avg_buf;
        map<int,double*> rms_buf;

        map<int,double> all_img_count;
        map<int,double> all_img_avg;
        map<int,double> all_img_rms;
        vector<pixel_type> pixelTypeVecGlobal;

        map<int,double> img_count;
        map<int,double> img_avg;
        map<int,double> img_rms;

        int nptime = 0; double *ptime = NULL;
        double myperiod = period; double pbeg = 0.0;
        if(b_lock)
          {
            cout << "Reading lock file " << lockstr << " ..." << endl;
            // Setup phase lock temporal structure
            // Lock file (lockstr) is binary with # of nodes
            // followed the by the list of nodal time values
            FILE *myfh = NULL;
            myfh=fopen(lockstr.c_str(),"rb");
            if(myfh == NULL) { CERR("Warning: Unable to open lock file " << lockstr); }
            fread(&nptime,sizeof(int),1,myfh);
            assert(nptime > 0);
            ptime = new double[nptime];
            fread(ptime,sizeof(double),nptime,myfh);
            fclose(myfh);
            if(nptime <= 1)
              { CERR("Invalid lock file with " << nptime << " nodes (node count < 2)."); }

            // If PERIOD is not positive, then compute the average period
            // from the lock file and use that value instead (could set the
            // the PERIOD arguments to be optional if lock file is given)
            if(period <= 0)
              {
                period = 0;
                for(int i=0; i < nptime-1; ++i) period += ptime[i+1]-ptime[i];
                period /= double(nptime-1);
                cout << "Reset PERIOD to " << period << "." << endl;
              }
          }

        // Setup cyclic buffer index into image bins
        std::vector< std::vector< std::pair<int,double> > > buf2img(nbuf); // Holds the images/times associated with each buffer
        std::vector< std::pair<int,double> >::iterator it_buf2img;
        std::vector<int> bin2buf(nbin); std::vector<int> buf2bin(nbuf);
        setBinBufIndex(bindex, nbin, nbuf, bin2buf, buf2bin); const int bin2buf_bindex = bin2buf[bindex];

        // Initial loop through the image sequence to compute initial estimates of bin-average/rms statistics
        for (int iI = 0; iI<pngCache->size(); ++iI){

          //COUT1(" > working on image \"" << pngCache->imageList[iI] << "\"");

          double time = time0;
          if (b_imgdt)
            {
              // If IMGDT set, increment time from time0 by imgdt for each consecutive image
              time = time0 + imgdt*iI;
            }
          else
            {
              // If no IMGDT set, read from image; if TIME not set in image, time = time0;
              time = time0 + img_getTime(pngCache->imageList[iI]);
            }

          myperiod = period;
          if(!b_lock)
            {
              // The constant period given on input is used as no lock file provided
              // Time and period are assumed non-negative
              while (time >= myperiod)
                time -= myperiod;
            }
          else
            {
              // Determine the local time offset from the local node-to-node period defined by the lock file
              int pbeg_found = 0; // Was a starting index for a given cycle/period found?
              if(time >= ptime[0]-period && time < ptime[0]) // The average period is assumed to extend before/after the beg/end nodes
                { pbeg = ptime[0]-period; myperiod=period; pbeg_found = 1; }
              else if(time >= ptime[nptime-1] && time < ptime[nptime-1]+period)
                { pbeg = ptime[nptime-1]; myperiod=period; pbeg_found = 1; }
              else if(nptime > 1)
                {
                  // Compute the local node-to-node period defined by the lock file and set the time offset
                  for(int p=0; p < nptime-1; ++p)
                    if(time >= ptime[p] && time < ptime[p+1])
                      { pbeg = ptime[p]; myperiod = ptime[p+1]-ptime[p]; pbeg_found = 1; }
                }
              if(pbeg_found) time -= pbeg; // Subtract the time offset for this local cycle period
              else
                { CERR("Could not find time " << time << " for image " << pngCache->imageList[iI] << " in lock file " << lockstr); }
            }

          const int ibin = (int)floor(time/myperiod*double(nbin));
          assert((ibin >= 0)&&(ibin < nbin));

          const int ibuf = bin2buf[ibin];
          if(ibuf < 0) continue;

          cout << " > working on image " << iI+1 << " of " << pngCache->size() << " \"" << pngCache->imageList[iI] << "\"..." << endl;
          cout << " > got time: " << time << " time/period*nbin: " << time/myperiod*double(nbin) << endl;
          buf2img[ibuf].push_back(std::pair<int,double>(iI,time/myperiod*double(nbin)));

          PngData * img = pngCache->getImage(iI);
          assert(img->getNx() == nx); // if you hit one of these, images are different sizes -- need to handle in some future version
          assert(img->getNy() == ny);

          vector<pixel_type> pixelTypeVec;
          img->getPixelTypes(pixelTypeVec);

          img_count.clear();
          img_avg.clear();
          img_rms.clear();
          for (int ipt=0;ipt<pixelTypeVec.size();++ipt){
            img_count[pixelTypeVec[ipt].flag] = 0.0;
            img_avg[pixelTypeVec[ipt].flag] = 0.0;
            img_rms[pixelTypeVec[ipt].flag] = 0.0;
            if (count_buf.find(pixelTypeVec[ipt].flag)==count_buf.end()){
              //first time we have encountered this data type, allocate buffer
              depth_buf[pixelTypeVec[ipt].flag] = new double[nPx*nbuf];
              count_buf[pixelTypeVec[ipt].flag] = new double[nPx*nbuf];
              avg_buf[pixelTypeVec[ipt].flag]  = new double[nPx*nbuf];
              rms_buf[pixelTypeVec[ipt].flag]   = new double[nPx*nbuf];
              for (int iPx = 0; iPx<nPx*nbuf; ++iPx){
                depth_buf[pixelTypeVec[ipt].flag][iPx] = 0.0;
                count_buf[pixelTypeVec[ipt].flag][iPx] = 0.0;
                avg_buf[pixelTypeVec[ipt].flag][iPx]  = 0.0;
                rms_buf[pixelTypeVec[ipt].flag][iPx]   = 0.0;
              }
              pixelTypeVecGlobal.push_back(pixelTypeVec[ipt]);
            }
            else{
              //this data type has appeared in other images, check
              //that we are working with the same variable name
              int iptg;
              for (iptg=0; iptg<pixelTypeVecGlobal.size(); ++iptg){
                if (pixelTypeVecGlobal[iptg].flag == pixelTypeVec[ipt].flag){
                  if (pixelTypeVecGlobal[iptg].varId !=  pixelTypeVec[ipt].varId){
                    CERR("Images contain " << pixelTypeVecGlobal[iptg].name << " data with different variables: "
                         << pixelTypeVecGlobal[iptg].varId << " and " << pixelTypeVec[ipt].varId);
                  }
                  else{
                    break;
                  }
                }
              }
              assert(iptg<pixelTypeVecGlobal.size()); //global vec should have a matching data type flag
            }
          }

          for (int iPx = 0; iPx<nPx; ++iPx) {
            if (img->pixel_flag[iPx] >= PngData::VOLUME_DATA_PIXEL) {
              const double val = img->pixel_data[iPx];

              depth_buf[img->pixel_flag[iPx]][iPx+nPx*ibuf] += double(img->getDepth(iPx));
              count_buf[img->pixel_flag[iPx]][iPx+nPx*ibuf] += 1.0; // every pixel gets equal weight for now
              avg_buf[ img->pixel_flag[iPx]][iPx+nPx*ibuf] += val;
              rms_buf[  img->pixel_flag[iPx]][iPx+nPx*ibuf] += val*val; // for now use rms to store sum(value^2)

              img_count[img->pixel_flag[iPx]]      += 1.0;
              img_avg[img->pixel_flag[iPx]]       += val;
              img_rms[img->pixel_flag[iPx]]        += val*val;
            }//(img->pixel_flag[iPx] >= PngData::VOLUME_DATA_PIXEL)
          }//for (iPx < nPx)
          img->lock = false;

          for (int ipt=0; ipt<pixelTypeVec.size(); ++ipt){
            const int dataFlag = pixelTypeVec[ipt].flag;
            if (all_img_count.find(dataFlag)!=all_img_count.end()){
              all_img_count[dataFlag] += img_count[dataFlag];
              all_img_avg[dataFlag]  += img_avg[dataFlag];
              all_img_rms[dataFlag]   += img_rms[dataFlag];
            }
            else{ //first time we found this data type in an image...
              all_img_count[dataFlag] = img_count[dataFlag];
              all_img_avg[dataFlag]  = img_avg[dataFlag];
              all_img_rms[dataFlag]   = img_rms[dataFlag];
            }

            assert(img_count[dataFlag]>0); //pixelTypeVec should only include entries when data pixels are present...
            img_avg[dataFlag] /= img_count[dataFlag];
            img_rms[dataFlag] = sqrt(max(0.0,img_rms[dataFlag]/img_count[dataFlag]-img_avg[dataFlag]*img_avg[dataFlag]));
            cout << " > image \"" << pngCache->imageList[iI] << "\" - " << pixelTypeVec[ipt].name << " data - avg: " << img_avg[dataFlag] << " rms: " << img_rms[dataFlag] << endl;
          }
        }//for (iI < iI<pngCache->size())

        for(int ipt=0; ipt<pixelTypeVecGlobal.size(); ++ipt)
          {
            const int dataFlag = pixelTypeVecGlobal[ipt].flag;
            if (all_img_count[dataFlag] == 0.0) {
              cout << " > all images have no active " <<  pixelTypeVecGlobal[ipt].name << " data pixels." << endl;
              assert(0); //should not get here...
            }
            else {
              all_img_avg[dataFlag] /= all_img_count[dataFlag];
              all_img_rms[dataFlag] = sqrt(max(0.0,all_img_rms[dataFlag]/all_img_count[dataFlag]-all_img_avg[dataFlag]*all_img_avg[dataFlag]));
              cout << " > all images - " << pixelTypeVecGlobal[ipt].name << " data - avg: " << all_img_avg[dataFlag] << " rms: " << all_img_rms[dataFlag] << endl;
            }
          }

        // Compute the bin-averaged values for each pixel
        for(int ipt=0; ipt < pixelTypeVecGlobal.size(); ++ipt)
          {
            const int dataFlag = pixelTypeVecGlobal[ipt].flag;
            for(int ibuf=0; ibuf < nbuf; ++ibuf)
              {
                for(int iPx=0; iPx < nPx; ++iPx)
                  {
                    if(count_buf[dataFlag][iPx+nPx*ibuf] > 0)
                      {
                        avg_buf[dataFlag][iPx+nPx*ibuf]  /= count_buf[dataFlag][iPx+nPx*ibuf];
                        depth_buf[dataFlag][iPx+nPx*ibuf] /= count_buf[dataFlag][iPx+nPx*ibuf];
                      }
                  }
              }
          }

        // Compute the second-order statistics using the selected phase-averaging method
        if(pmethod == BIN)
          {
            for(int ipt=0; ipt < pixelTypeVecGlobal.size(); ++ipt)
              {
                const int dataFlag = pixelTypeVecGlobal[ipt].flag;
                for(int iPx=0; iPx < nPx; ++iPx)
                  {
                    if(count_buf[dataFlag][iPx+nPx*(bin2buf_bindex)] > 0)
                      {
                        const double avg = avg_buf[dataFlag][iPx+nPx*bin2buf_bindex];
                        const double rms = sqrt(max(0.0,rms_buf[dataFlag][iPx+nPx*bin2buf_bindex]/count_buf[dataFlag][iPx+nPx*bin2buf_bindex]-avg*avg));
                        rms_buf[dataFlag][iPx+nPx*bin2buf_bindex] = rms;
                      }
                  }
              }
          }
        else
          {
            // Reset rms buffers for each data type
            for(int ipt=0; ipt < pixelTypeVecGlobal.size(); ++ipt) {
              const int dataFlag = pixelTypeVecGlobal[ipt].flag;
              for(int iPx=0; iPx < nPx; ++iPx) rms_buf[dataFlag][iPx+nPx*bin2buf_bindex] = 0.0;
            }

            // Read all buffer images to compute the total sum of the squared deviation
            double * px_count_buf = new double[nbuf]; double * px_avg_buf  = new double[nbuf]; double * px_time_buf  = new double[nbuf];
            for(it_buf2img=buf2img[bin2buf_bindex].begin(); it_buf2img != buf2img[bin2buf_bindex].end(); ++it_buf2img)
              {
                const int iI = (*it_buf2img).first; const double imgTimeIndex = (*it_buf2img).second;
                PngData * img = pngCache->getImage(iI);
                assert(img->getNx() == nx); // if you hit one of these, images are different sizes -- need to handle in some future version
                assert(img->getNy() == ny);

                // Set time index for this buffer image (shift time by 1/2 bin for bin average)
                //const double imgbuf = bin2buf_bindex + (imgTimeIndex-bindex) - 0.5;
                const int imgbin = (int)floor(imgTimeIndex); assert(bin2buf[imgbin] >= 0);
                const double imgbuf = bin2buf[imgbin] + (imgTimeIndex-imgbin) - 0.5;
                for(int ibuf=0; ibuf < nbuf; ++ibuf) px_time_buf[ibuf] = bin2buf[buf2bin[ibuf]];

                for(int iPx=0; iPx < nPx; ++iPx)
                  {
                    if(img->pixel_flag[iPx] >= PngData::VOLUME_DATA_PIXEL)
                      {
                        const int dataFlag = img->pixel_flag[iPx];
                        assert(count_buf[dataFlag][iPx+nPx*(bin2buf_bindex)] > 0);
                        const double val = img->pixel_data[iPx];

                        // Store pixel data from buf neighbors
                        for(int ibuf=0; ibuf < nbuf; ++ibuf) {
                          px_avg_buf[ibuf]  =  avg_buf[dataFlag][iPx+nPx*ibuf];
                          px_count_buf[ibuf] = count_buf[dataFlag][iPx+nPx*ibuf];
                        }

                        // Compute average pixel at image time from buf neighbors
                        const double iavg = calcPxBinAvg(bin2buf_bindex,imgbuf,nbuf,px_avg_buf,px_count_buf,px_time_buf,pmethod);

                        // Add squared deviation from the pixel average at this image time
                        rms_buf[dataFlag][iPx+nPx*bin2buf_bindex] += (val-iavg)*(val-iavg);
                      }
                  }
                img->lock = false;
              }
            delete[] px_count_buf; delete[] px_avg_buf; delete[] px_time_buf;

            // Compute rms for each data type
            for(int ipt=0; ipt < pixelTypeVecGlobal.size(); ++ipt)
              {
                const int dataFlag = pixelTypeVecGlobal[ipt].flag;
                for(int iPx=0; iPx < nPx; ++iPx)
                  {
                    if(count_buf[dataFlag][iPx+nPx*(bin2buf_bindex)] > 0)
                      {
                        // Compute square root of the expectation of the squared deviation from the pixel averages within the bindex buf bin
                        const double rms = sqrt(max(0.0,rms_buf[dataFlag][iPx+nPx*bin2buf_bindex]/count_buf[dataFlag][iPx+nPx*bin2buf_bindex]));
                        rms_buf[dataFlag][iPx+nPx*bin2buf_bindex] = rms;
                      }
                  }
              }
          }

        // Write the avg/rms data
        for(int ipt=0; ipt < pixelTypeVecGlobal.size(); ++ipt)
          {
            int dataFlag = pixelTypeVecGlobal[ipt].flag;

            // Write the bin-averaged average images for bindex for each pixel type
            writePhaseBinnedImages(bindex, dataFlag,
                                   avg_buf[dataFlag]  +nPx*bin2buf_bindex,
                                   count_buf[dataFlag] +nPx*bin2buf_bindex,
                                   depth_buf[dataFlag] +nPx*bin2buf_bindex,
                                   pixelTypeVecGlobal[ipt].name,pixelTypeVecGlobal[ipt].varId,"avg",
                                   outputPathAndPrefix,nx,ny,pngCache->getMetadata(0),
                                   range_avg,         b_range_avg,
                                   rangeSurface_avg,  b_rangeSurface_avg,
                                   rangeParticles_avg,b_rangeParticles_avg,
                                   rangeIso_avg,      b_rangeIso_avg);

            // Write the bin-averaged rms images for bindex for each pixel type
            writePhaseBinnedImages(bindex, dataFlag,
                                   rms_buf[dataFlag]   +nPx*bin2buf_bindex,
                                   count_buf[dataFlag] +nPx*bin2buf_bindex,
                                   depth_buf[dataFlag] +nPx*bin2buf_bindex,
                                   pixelTypeVecGlobal[ipt].name,pixelTypeVecGlobal[ipt].varId,"rms",
                                   outputPathAndPrefix,nx,ny,pngCache->getMetadata(0),
                                   range_rms,         b_range_rms,
                                   rangeSurface_rms,  b_rangeSurface_rms,
                                   rangeParticles_rms,b_rangeParticles_rms,
                                   rangeIso_rms,      b_rangeIso_rms);

            // Clean up memory for each data type...
            delete[] depth_buf[dataFlag];
            delete[] count_buf[dataFlag];
            delete[] avg_buf[dataFlag];
            delete[] rms_buf[dataFlag];
          }

        if(ptime != NULL) delete ptime;

      }// set_bindex loop end
  }//itpTemporalStatsPhaseLock()


  void itpDiff(Param * param){

    if ((!pngCache)||(pngCache->imageList.size() == 0)) {
      CERR("cannot compute DIFF: pngCache is empty.");
    }

    if (pngCache->imageList.size() == 1){
      CERR("cannot compute DIFF: pngCache has only one image.");
    }

    //Don't allow diffing of more images than we can cache in memory,
    //this is defined by pngCache->ncache.  It will allow the
    //global diff range to be computed prior to writing output and
    //without multiple reads per image.
    if (pngCache->imageList.size() > pngCache->ncache){
      CERR("cannot compute DIFF: total number of images must be <= " << pngCache->ncache );
    }

    string diffPath = "images_diff";
    bool b_signed = false;
    int ip = 0;
    if (param->size() > 0) {
      if (param->getString(ip)=="SIGNED"){
        b_signed = true;
        ++ip;
      }
      if (param->size() > ip){
        diffPath = param->getString(ip);
      }
    }
    COUT1("Writing files to " << diffPath << "\n");
    MiscUtils::mkdir_for_file(diffPath+"/test.png");

    PngData * imgBase = pngCache->getImage(0);
    PngData * imgDiff;

    double globalRangeVolume[2]    = {HUGE_VAL,-HUGE_VAL};
    double globalRangeSurface[2]   = {HUGE_VAL,-HUGE_VAL};
    double globalRangeParticles[2] = {HUGE_VAL,-HUGE_VAL};
    double globalRangeIso[2]       = {HUGE_VAL,-HUGE_VAL};

    vector<pixel_type> pixelTypeVec;
    imgBase->getPixelTypes(pixelTypeVec);

    vector<pair<bool,string> > msgVec;
    for (int iI = 1; iI<pngCache->size(); ++iI){
      imgDiff = pngCache->getImage(iI);
      string msg;
      if(imgBase->diff(imgDiff,b_signed,msg)){ //the diff is stored in imgDiff
        msgVec.push_back(make_pair(true,msg));
        for (int ipt=0; ipt<pixelTypeVec.size(); ++ipt){
          double imgRMin, imgRMax;
          switch (pixelTypeVec[ipt].flag){
          case PngData::VOLUME_DATA_PIXEL:
            imgRMin = imgDiff->getRangeMin();
            imgRMax = imgDiff->getRangeMax();
            if (imgRMin<globalRangeVolume[0])
              globalRangeVolume[0] = imgRMin;
            if (imgRMax>globalRangeVolume[1])
              globalRangeVolume[1] = imgRMax;
            break;
          case PngData::SURFACE_DATA_PIXEL:
            imgRMin = imgDiff->getRangeSurfaceMin();
            imgRMax = imgDiff->getRangeSurfaceMax();
            if (imgRMin<globalRangeSurface[0])
              globalRangeSurface[0] = imgRMin;
            if (imgRMax>globalRangeSurface[1])
              globalRangeSurface[1] = imgRMax;
            break;
          case PngData::PARTICLE_DATA_PIXEL:
            imgRMin = imgDiff->getRangeParticlesMin();
            imgRMax = imgDiff->getRangeParticlesMax();
            if (imgRMin<globalRangeParticles[0])
              globalRangeParticles[0] = imgRMin;
            if (imgRMax>globalRangeParticles[1])
              globalRangeParticles[1] = imgRMax;
            break;
          case PngData::ISO_DATA_PIXEL:
            imgRMin = imgDiff->getRangeIsoMin();
            imgRMax = imgDiff->getRangeIsoMax();
            if (imgRMin<globalRangeIso[0])
              globalRangeIso[0] = imgRMin;
            if (imgRMax>globalRangeIso[1])
              globalRangeIso[1] = imgRMax;
            break;
          default:
            assert(0); //shouldn't get here, diff will return false
          }
        }
      }
      else
        msgVec.push_back(make_pair(false,msg));
    }
    imgBase->lock = false; //release base image from cache

    //write diff images with the global range
    for (int iI = 1; iI<pngCache->size(); ++iI){
      COUT1("\n > Image " << pngCache->imageList[iI]);
      COUT1(" >   " << msgVec[iI-1].second);

      if (msgVec[iI-1].first){ //diff successfully produced output image

        imgDiff = pngCache->getImage(iI);

        //build output path and filename
        //use original filename in a different folder
        vector<string> inputNameTokens;
        tokenizeString(inputNameTokens,pngCache->imageList[iI],"/");
        string outputPathAndFilename = diffPath + "/diff.";
        outputPathAndFilename += inputNameTokens[inputNameTokens.size()-1];

        cout << " >   writing diff to \"" << outputPathAndFilename << "\"..." << endl;

        for (int ipt=0; ipt<pixelTypeVec.size(); ++ipt){
          switch (pixelTypeVec[ipt].flag){
          case PngData::VOLUME_DATA_PIXEL:
            imgDiff->setRange(globalRangeVolume);
            break;
          case PngData::SURFACE_DATA_PIXEL:
            imgDiff->setRangeSurface(globalRangeSurface);
            break;
          case PngData::PARTICLE_DATA_PIXEL:
            imgDiff->setRangeParticles(globalRangeParticles);
            break;
          case PngData::ISO_DATA_PIXEL:
            imgDiff->setRangeIso(globalRangeIso);
            break;
          default:
            assert(0); //shouldn't get here, diff will return false
          }
        }

        imgDiff->initializeWrite();
        imgDiff->write(outputPathAndFilename.c_str());

      }
      imgDiff->lock = false; //release all images from cache
    }
  }//itpDiff()

  void itpDiffStats(Param * param) {

    if ((!pngCache)||(pngCache->imageList.size() == 0)) {
      CERR("cannot compute DIFF_STATS: pngCache is empty.");
    }

    if (pngCache->imageList.size() == 1){
      CERR("cannot compute DIFF_STATS: pngCache has only one image.");
    }

    if (pngCache->imageList.size() > 2){
      CWARN("can only compute stats on difference between two images; using first two and ignoring others");
    }

    int iarg = 0;
    string metric = "AE";
    double tol = -1.0;
    bool force_image = false;
    while (iarg < param->size()) {
      const string token=param->getString(iarg++);
      if (token == "METRIC") {
        metric = param->getString(iarg++);
      }
      else if (token == "FORCE_IMAGE") {
        force_image = true;
      }
      else if (token == "DIFF_TOL") {
        tol = param->getDouble(iarg++);
      }
      else {
        CWARN("unrecognized DIFF_STATS parameter \"" << token << "\"; skipping");
      }
    }

    // use reference filename for prefix
    string outputPrefix = getPrefix(pngCache->imageList[0]);
    COUT2(" > output file prefix will be \"" << outputPrefix << "\"");

    PngData * imgBase = pngCache->getImage(0,false);
    PngData * imgOther = pngCache->getImage(1,false);
    imgBase->diffWithStats(imgOther,metric,outputPrefix,tol,force_image);  //the diff-image is stored in imgOther

    //release all images from cache
    imgBase->lock = false;
    imgOther->lock = false;
  }  //itpDiffStats()


  void itpRescale(Param * param){

    COUT1("RESCALE...");

    if ((!pngCache)||(pngCache->imageList.size() == 0)) {
      CERR("cannot RESCALE: pngCache is empty.");
    }

    string rescalePath = "images_rescaled";

    double rangeVolume[2] = {0.0f,0.0f};
    double rangeSurface[2] = {0.0f,0.0f};
    double rangeParticles[2] = {0.0f,0.0f};
    double rangeIso[2] = {0.0f,0.0f};
    bool b_rangeVolume = false;
    bool b_rangeSurface = false;
    bool b_rangeParticles = false;
    bool b_rangeIso = false;

    if (param->size()==0){
      CERR("Invalid RESCALE syntax, examples: \n" <<
           " RESCALE <min> <max> \n" <<
           " RESCALE GLOBAL [TIGHT,LOOSE] \n" <<
           " RESCALE VOLUME <min> <max> SURFACE <min> <max> PARTICLES <min> <max> ISO <min> <max>");
    }

    //1. look for simple syntax:
    //RESCALE <min> <max>

    bool b_simpleRange = false;
    if (param->size()>=2){

      //if not two doubles, look for alternate syntax
      bool b_foundNumber0, b_foundNumber1;
      b_foundNumber0 = from_string<double>(rangeVolume[0],param->getString(0),std::dec);
      b_foundNumber1 = from_string<double>(rangeVolume[1],param->getString(1),std::dec);

      if (b_foundNumber0&&b_foundNumber1){
        b_simpleRange = true;
        //later, any found data types will be rescaled if they have the same varId
        b_rangeVolume = true;
        b_rangeSurface = true;
        b_rangeParticles = true;
        b_rangeIso = true;
        rangeSurface[0]   = rangeVolume[0];
        rangeSurface[1]   = rangeVolume[1];
        rangeParticles[0] = rangeVolume[0];
        rangeParticles[1] = rangeVolume[1];
        rangeIso[0]       = rangeVolume[0];
        rangeIso[1]       = rangeVolume[1];


        if (param->size()>=4){
          assert(param->getString(2)=="PATH");
          rescalePath = param->getString(3);
        }
      }
    }

    //2. look for "advanced" syntax
    bool bTightRescale = false;
    bool bLooseRescale = false;

    if (!b_simpleRange){
      int iarg = 0;
      while (iarg < param->size()) {
        const string token = param->getString(iarg++);
        if (token == "GLOBAL") {
          string type = param->getString(iarg++);
          if (type=="TIGHT")
            bTightRescale = true;
          else if (type=="LOOSE")
            bLooseRescale = true;
          else
            CERR("Invalid GLOBAL rescale type -> " << type << "\n" <<
                 "  valid types are TIGHT and LOOSE");
        }
        else if (token == "VOLUME"){
          b_rangeVolume = true;
          rangeVolume[0] = param->getFloat(iarg++);
          rangeVolume[1] = param->getFloat(iarg++);
        }
        else if (token == "SURFACE"){
          b_rangeSurface = true;
          rangeSurface[0] = param->getFloat(iarg++);
          rangeSurface[1] = param->getFloat(iarg++);
        }
        else if (token == "PARTICLES"){
          b_rangeParticles = true;
          rangeParticles[0] = param->getFloat(iarg++);
          rangeParticles[1] = param->getFloat(iarg++);
        }
        else if (token == "ISO"){
          b_rangeIso = true;
          rangeIso[0] = param->getFloat(iarg++);
          rangeIso[1] = param->getFloat(iarg++);
        }
        else if (token == "PATH"){
          rescalePath = param->getString(iarg++);
        }
        else {
          CERR("Invalid RESCALE syntax, examples: \n" <<
               " RESCALE <min> <max> \n" <<
               " RESCALE GLOBAL [TIGHT,LOOSE] \n" <<
               " RESCALE VOLUME <min> <max> SURFACE <min> <max> PARTICLES <min> <max> ISO <min> <max>");
        }
      }
    }

    COUT1("Writing files to " << rescalePath << "\n");
    MiscUtils::mkdir_for_file(rescalePath+"/test.png");

    //3. Finished parsing, prep ranges if global option requested

    if (bTightRescale||bLooseRescale){
      assert(!b_simpleRange);
      double flip[2];
      if (bTightRescale){
        flip[0] = -1.0f; flip[1] = 1.0f;
      }
      else{
        flip[0] = 1.0f; flip[1] = -1.0f;
      }
      rangeVolume[0]    = HUGE_VAL;
      rangeVolume[1]    = HUGE_VAL;
      rangeSurface[0]   = HUGE_VAL;
      rangeSurface[1]   = HUGE_VAL;
      rangeParticles[0] = HUGE_VAL;
      rangeParticles[1] = HUGE_VAL;
      rangeIso[0]       = HUGE_VAL;
      rangeIso[1]       = HUGE_VAL;
      b_rangeVolume = false;
      b_rangeSurface = false;
      b_rangeParticles = false;
      b_rangeIso = false;

      for (int iI = 0; iI<pngCache->size(); ++iI) {
        ImageMetadata * imd = pngCache->getMetadata(iI);
        if (imd){
          double rmin, rmax;
          if (imd->hasText("VAR") && imd->getVarId() != "OFF" &&
              imd->hasText("RANGE_MAX") && imd->hasText("RANGE_MIN")){
            b_rangeVolume = true;
            rmin = flip[0]*imd->getRangeMin();
            rmax = flip[1]*imd->getRangeMax();
            if (rmin < rangeVolume[0])
              rangeVolume[0] = rmin;
            if (rmax < rangeVolume[1])
              rangeVolume[1] = rmax;
          }

          if (imd->hasText("VAR_ON_SURFACE") && imd->getSurfVarId() != "OFF" &&
              imd->hasText("RANGE_ON_SURFACE_MAX") && imd->hasText("RANGE_ON_SURFACE_MIN")){
            b_rangeSurface = true;
            rmin = flip[0]*imd->getSurfRangeMin();
            rmax = flip[1]*imd->getSurfRangeMax();
            if (rmin < rangeSurface[0])
              rangeSurface[0] = rmin;
            if (rmax < rangeSurface[1])
              rangeSurface[1] = rmax;
          }

          if (imd->hasText("VAR_ON_PARTICLE") && imd->getPartVarId() != "OFF" &&
              imd->hasText("RANGE_ON_PARTICLE_MAX") && imd->hasText("RANGE_ON_PARTICLE_MIN")){
            b_rangeParticles = true;
            rmin = flip[0]*imd->getPartRangeMin();
            rmax = flip[1]*imd->getPartRangeMax();
            if (rmin < rangeParticles[0])
              rangeParticles[0] = rmin;
            if (rmax < rangeParticles[1])
              rangeParticles[1] = rmax;
          }

          if (imd->hasText("VAR_ON_ISO") && imd->getIsoVarId() != "OFF" &&
              imd->hasText("RANGE_ON_ISO_MAX") && imd->hasText("RANGE_ON_ISO_MIN")){
            b_rangeIso = true;
            rmin = flip[0]*imd->getIsoRangeMin();
            rmax = flip[1]*imd->getIsoRangeMax();
            if (rmin < rangeIso[0])
              rangeIso[0] = rmin;
            if (rmax < rangeIso[1])
              rangeIso[1] = rmax;
          }
        }
      }
      rangeVolume[0]    *= flip[0];
      rangeVolume[1]    *= flip[1];
      rangeSurface[0]   *= flip[0];
      rangeSurface[1]   *= flip[1];
      rangeParticles[0] *= flip[0];
      rangeParticles[1] *= flip[1];
      rangeIso[0]       *= flip[0];
      rangeIso[1]       *= flip[1];

    }

    //3. Range specification complete...summarize
    if (b_simpleRange){
      COUT1("Attempting to rescale image data to the range " << rangeVolume[0] << ":" << rangeVolume[1]);
      COUT1("  If images contain multiple data types (i.e. volume and surface),\n" <<
            "  data will only be rescaled if the contour variable is identical\n" <<
            "  for all types (i.e. only pressure is visualized).  If your image\n" <<
            "  visualizes more than one variable, please specify a rescale range\n" <<
            "  explicitly for each data type. Use parameter HELP for syntax info.");
    }
    else{
      if (bTightRescale){
        COUT1("Using TIGHT rescale range computed from image metadata,\n" <<
              "  for the following data types:");
      }
      else if (bLooseRescale){
        COUT1("Using LOOSE rescale range computed from image metadata,\n"
              "  for the following data types:");
      }
      else{
        COUT1("Rescaling requested for the following data types:")
          }

      if (b_rangeVolume){
        COUT1("Volume data " << rangeVolume[0] << ":" << rangeVolume[1]);
      }
      if (b_rangeSurface){
        COUT1("Surface data " << rangeSurface[0] << ":" << rangeSurface[1]);
      }
      if (b_rangeParticles){
        COUT1("Particles data " << rangeParticles[0] << ":" << rangeParticles[1]);
      }
      if (b_rangeIso){
        COUT1("Iso data " << rangeIso[0] << ":" << rangeIso[1]);
      }
      COUT1("");

      if (bTightRescale){
        if (b_rangeVolume&&(rangeVolume[0]>=rangeVolume[1]))
          CERR("invalid TIGHT range for Volume data, min>=max");
        if (b_rangeSurface&&(rangeSurface[0]>=rangeSurface[1]))
          CERR("invalid TIGHT range for Surface data, min>=max");
        if (b_rangeParticles&&(rangeParticles[0]>=rangeParticles[1]))
          CERR("invalid TIGHT range for Particles data, min>=max");
        if (b_rangeIso&&(rangeIso[0]>=rangeIso[1]))
          CERR("invalid TIGHT range for Iso data, min>=max");

      }

    }

    //4. Rescale, if b_simpleRange test for data uniformity
    for (int iI = 0; iI<pngCache->size(); ++iI) {
      PngData * img = pngCache->getImage(iI);
      vector<pixel_type> pixelTypeVec;
      img->getPixelTypes(pixelTypeVec);

      if (pixelTypeVec.size()>0){
        bool b_varMismatch = false;
        COUT1(" > Image " << pngCache->imageList[iI] << " with ranges:");
        string varId = pixelTypeVec[0].varId;
        for (int ipt=0; ipt<pixelTypeVec.size(); ++ipt){
          COUT1(" >   " << pixelTypeVec[ipt].name << " " << pixelTypeVec[ipt].varId << " " << pixelTypeVec[ipt].rmin << ":" << pixelTypeVec[ipt].rmax);
          if (b_simpleRange && varId!=pixelTypeVec[ipt].varId){
            b_varMismatch = true;
          }
        }
        if (b_varMismatch){
          COUT1(" > Warning - skipping rescale of image with multiple contour variables - you must specify multiple ranges");
        }
        else{
          for (int ipt=0; ipt<pixelTypeVec.size(); ++ipt){
            switch (pixelTypeVec[ipt].flag){
            case PngData::VOLUME_DATA_PIXEL:
              if (b_rangeVolume){
                img->setRange(rangeVolume);
              }
              break;
            case PngData::SURFACE_DATA_PIXEL:
              if (b_rangeSurface){
                img->setRangeSurface(rangeSurface);
              }
              break;
            case PngData::PARTICLE_DATA_PIXEL:
              if (b_rangeParticles){
                img->setRangeParticles(rangeParticles);
              }
              break;
            case PngData::ISO_DATA_PIXEL:
              if (b_rangeIso){
                img->setRangeIso(rangeIso);
              }
              break;
            default:
              COUT1("Warning - unrecognized pixel data type " << pixelTypeVec[ipt].flag);
            }
          }

          //build output path and filename
          //use original filename in a different folder
          vector<string> inputNameTokens;
          tokenizeString(inputNameTokens,pngCache->imageList[iI],"/");
          string outputPathAndFilename = rescalePath + "/";
          outputPathAndFilename += inputNameTokens[inputNameTokens.size()-1];

          img->initializeWrite();
          img->write(outputPathAndFilename.c_str());
        }
      }
      img->lock = false;
    }
  }//itpRescale()

  void itpRecolor(Param * param){

    if ((!pngCache)||(pngCache->imageList.size() == 0)) {
      CERR("cannot RECOLOR: pngCache is empty.");
    }

    string colorMapName;
    string colorMapSurfaceName;
    string colorMapParticlesName;
    string colorMapIsoName;
    bool b_colorMapName          = false;
    bool b_colorMapSurfaceName   = false;
    bool b_colorMapParticlesName = false;
    bool b_colorMapIsoName       = false;

    bool b_legend                = false;
    bool b_legend_surface        = false;
    bool b_legend_particles      = false;
    bool b_legend_iso            = false;

    double legend_x = 0;
    double legend_y = 0;
    double legend_surface_x = 0;
    double legend_surface_y = 1;
    double legend_particles_x = 1;
    double legend_particles_y = 0;
    double legend_iso_x = 1;
    double legend_iso_y = 1;

    bool b_background = false;
    unsigned char background_r,background_g,background_b;

    bool b_checktime = false;
    
    string recolorPath = "images_recolored";
    ColorMap colorMap;

    stringstream validColorMapNames;
    ColorMap::getValidColorMapNames(validColorMapNames);
    //string validColorMapNames = "GRAYSCALE, HOT_METAL, BLUE_RED_RAINBOW, COOL_WARM, RAINBOW_DESATURATED, INVERTED_GRAYSCALE, INVERTED_HOT_METAL, INVERTED_BLUE_RED_RAINBOW, INVERTED_COOL_WARM, INVERTED_RAINBOW_DESATURATED";

    int iarg = 0;
    while (iarg < param->size()) {
      string token = param->getString(iarg++);
      if ((token == "COLORMAP")||(token == "COLORMAP_VOLUME")||(token == "COLORMAP_PLANE")) {
        colorMapName = param->getString(iarg++);
        if (!ColorMap::checkColorMapName(colorMapName)){
          CERR("cannot RECOLOR: " << colorMapName << " is not a valid colormap name, valid names: " <<  validColorMapNames.str());
        }
        b_colorMapName          = true;
      }
      else if (token == "COLORMAP_SURFACE"){
        colorMapSurfaceName = param->getString(iarg++);
        if (!ColorMap::checkColorMapName(colorMapSurfaceName)){
          CERR("cannot RECOLOR: " << colorMapSurfaceName << " is not a valid colormap name, valid names: " << validColorMapNames.str());
        }
        b_colorMapSurfaceName   = true;
      }
      else if (token == "COLORMAP_PARTICLES"){
        colorMapParticlesName = param->getString(iarg++);
        if (!ColorMap::checkColorMapName(colorMapParticlesName)){
          CERR("cannot RECOLOR: " << colorMapParticlesName << " is not a valid colormap name, valid names: " << validColorMapNames.str());
        }
        b_colorMapParticlesName = true;
      }
      else if (token == "COLORMAP_ISO"){
        colorMapIsoName = param->getString(iarg++);
        if (!ColorMap::checkColorMapName(colorMapIsoName)){
          CERR("cannot RECOLOR: " << colorMapIsoName << " is not a valid colormap name, valid names: " << validColorMapNames.str());
        }
        b_colorMapIsoName       = true;
      }
      else if (token == "PATH"){
        recolorPath = param->getString(iarg++);
      }
      else if (token == "CHECKTIME"){
        b_checktime = true;
        cout << " > CHECKTIME mode: will check file times and only do new work" << endl;
      }
      else if (token == "BACKGROUND"){
        b_background = true;
        token = param->getString(iarg++);
        if (token == "WHITE") {
          background_r = 255;
          background_g = 255;
          background_b = 255;
          cout << " > BACKGROUND will be WHITE" << endl;
        }
        else if (token == "BLACK") {
          background_r = 0;
          background_g = 0;
          background_b = 0;
          cout << " > BACKGROUND will be BLACK" << endl;
        }
        else {
          CERR("unrecognized BACKGROUND color: " << token);
        }
      }
      else if (token == "LEGEND"){
        b_legend = true;
        int i_test = iarg;
        if (i_test < param->size() && param->getString(i_test++) == "IJ"){
          legend_x = param->getDouble(i_test++);
          legend_y = param->getDouble(i_test++);
          iarg = i_test;
        }
      }
      else if (token == "LEGEND_SURFACE"){
        b_legend_surface = true;
        int i_test = iarg;
        if (i_test < param->size() && param->getString(i_test++) == "IJ"){
          legend_surface_x = param->getDouble(i_test++);
          legend_surface_y = param->getDouble(i_test++);
          iarg = i_test;
        }
      }
      else if (token == "LEGEND_PARTICLES"){
        b_legend_particles = true;
        int i_test = iarg;
        if (i_test < param->size() && param->getString(i_test++) == "IJ"){
          legend_particles_x = param->getDouble(i_test++);
          legend_particles_y = param->getDouble(i_test++);
          iarg = i_test;
        }
      }
      else if (token == "LEGEND_ISO"){
        b_legend_iso = true;
        int i_test = iarg;
        if (i_test < param->size() && param->getString(i_test++) == "IJ"){
          legend_iso_x = param->getDouble(i_test++);
          legend_iso_y = param->getDouble(i_test++);
          iarg = i_test;
        }
      }
      else {
        CERR("cannot RECOLOR: unrecognized token " << token << "\n" <<
             "  use syntax 'RECOLOR COLORMAP <colormap name> [COLORMAP_SURFACE <colormap name>]\n" <<
             "[COLORMAP_PARTICLES <colormap name>] [COLORMAP_ISO <colormap name>] \n" <<
             "[PATH <path>] [LEGEND]'\n" <<
             "  valid colormap names: " << validColorMapNames.str());
      }
    }

    COUT1("Recoloring images...");
    
    if (b_colorMapName) COUT1("Volume data to " << colorMapName);
    if (b_colorMapSurfaceName) COUT1("Surface data to " << colorMapSurfaceName);
    if (b_colorMapParticlesName) COUT1("Particle data to " << colorMapParticlesName);
    if (b_colorMapIsoName) COUT1("Iso data to " << colorMapIsoName);

    if(b_legend) COUT1("Adding volume data legend...");
    if(b_legend_surface) COUT1("Adding surface data legend...");
    if(b_legend_particles) COUT1("Adding particle data legend...");
    if(b_legend_iso) COUT1("Adding iso data legend...");

    COUT1("Writing files to " << recolorPath);

    MiscUtils::mkdir_for_file(recolorPath+"/test.png");

    for (int iI = 0; iI<pngCache->size(); ++iI) {

      //use original filename in a different folder
      string outputPathAndFilename = recolorPath + "/" + MiscUtils::getFilenameNoPath(pngCache->getImageName(iI));
      
      cout << " > recoloring image " << iI << " of " << pngCache->size() << " " << pngCache->getImageName(iI) << " to " << outputPathAndFilename << "...";
      
      // if in "CHECKTIME" mode, only do work if the dates are out of order: like a Makefile
      if (b_checktime) {
        if (MiscUtils::fileExists(outputPathAndFilename)) {
          struct stat t_stat_out;
          stat(outputPathAndFilename.c_str(),&t_stat_out);
          struct stat t_stat_in;
          stat(pngCache->getImageName(iI).c_str(),&t_stat_in);
          //cout << "infile time: " << t_stat_in.st_ctime << " outfile time: " << t_stat_out.st_ctime << endl;
          if (t_stat_out.st_ctime > t_stat_in.st_ctime) {
            cout << "skipping. Turn off CHECKTIME mode to recolor this image" << endl;
            continue;
          }
        }
      }
      
      PngData * img = pngCache->getImage(iI,false);
      
      if (b_background) img->setBackgroundColor(background_r,background_g,background_b);

      if (b_colorMapName) img->setColorMap(colorMapName);
      if (b_colorMapSurfaceName) img->setColorMapSurface(colorMapSurfaceName);
      if (b_colorMapParticlesName) img->setColorMapParticles(colorMapParticlesName);
      if (b_colorMapIsoName) img->setColorMapIso(colorMapIsoName);

      if(b_legend) img->setLegend(legend_x,legend_y);
      if(b_legend_surface) img->setLegendSurface(legend_surface_x,legend_surface_y);
      if(b_legend_particles) img->setLegendParticles(legend_particles_x,legend_particles_y);
      if(b_legend_iso) img->setLegendIso(legend_iso_x,legend_iso_y);

      img->initializeWriteRecolor(b_colorMapName,b_colorMapSurfaceName,b_colorMapParticlesName,b_colorMapIsoName);

      // write image to a tmp file...
      string tmpPathAndFilename = recolorPath + "/." + MiscUtils::getFilenameNoPath(pngCache->getImageName(iI));
      img->write(tmpPathAndFilename.c_str());
      img->lock = false;
      
      // then rename...
      rename(tmpPathAndFilename.c_str(),outputPathAndFilename.c_str());
      
      cout << "done" << endl;
    
    }

  }//itpRecolor()


  void itpMarkMedialAxis(Param * param){

    if ((!pngCache)||(pngCache->imageList.size() == 0)) {
      CERR("cannot MARK_MEDIAL_AXIS: pngCache is empty.");
    }

    string outputPath = "output_mark_medial_axis";

    int iarg = 0;
    string zonesA = "";
    string zonesB = "";
    vector<uint2> zonesAVec;
    vector<uint2> zonesBVec;
    int line_start_ij[2] = {0,-1};
    while (iarg < param->size()) {
      string token = param->getString(iarg++);
      if (token == "PATH"){
        outputPath = param->getString(iarg++);
      }
      else if (token == "ZONES_A"){
        zonesA = param->getString(iarg++);
        MiscUtils::splitCsv(zonesAVec,zonesA);
      }
      else if (token == "ZONES_B"){
        zonesB = param->getString(iarg++);
        MiscUtils::splitCsv(zonesBVec,zonesB);
      }
      else if (token == "LINE_START_IJ"){
        line_start_ij[0] = param->getInt(iarg++);
        line_start_ij[1] = param->getInt(iarg++);
      }
      else {
        CERR("cannot MARK_MEDIAL_AXIS: unrecognized token " << token << "\n" <<
             "  use syntax 'MARK_MEDIAL_AXIS'");
      }
    }

    COUT1("MARK_MEDIAL_AXIS...files will be written to " << outputPath);
    
    MiscUtils::mkdir_for_file(outputPath+"/test.png");

    for (int iI = 0; iI<pngCache->size(); ++iI) {

      //use original filename in a different folder
      string outputPathAndFilename = outputPath + "/" + MiscUtils::getFilenameNoPath(pngCache->getImageName(iI));
      
      cout << " > working on image " << iI+1 << " of " << pngCache->size() << ". Input image: " << pngCache->getImageName(iI) << "; Output image: " << outputPathAndFilename << endl;
      
      PngData * img = pngCache->getImage(iI);

      if (zonesAVec.size()==0||zonesBVec.size()==0){
        vector<uint2> zoneIds;
        img->getZoneIds(zoneIds);
        vector<string> zoneNames;
        img->getZoneNames(zoneNames);

        stringstream zoneSS;
        for (int iz=0; iz<zoneNames.size(); ++iz){
          zoneSS << "\n" << zoneIds[iz] << " " << zoneNames[iz];
        }
        CERR("cannot MARK_MEDIAL_AXIS, use syntax MARK_MEDIAL_AXIS ZONES_A zone0,zone2,...,zone5 ZONES_B zone1,...,zone6\n" <<
               "where valid zone names are: " << zoneSS.str());
      }

      cout << "    flagContiguousZoneIdPixels()...";
      img->flagContiguousZoneIdPixels(line_start_ij[0],line_start_ij[1],PngData::ISO_DATA_PIXEL);
      cout << "done" << endl;
      cout << "    computeDistanceFromFlaggedPixels()...";
      img->computeDistanceFromFlaggedPixels();
      cout << "done" << endl;
      cout << "    flagZoneGroup()...";
      img->flagZoneGroup(zonesAVec,PngData::FLAG_GREEN);
      img->flagZoneGroup(zonesBVec,PngData::FLAG_BLUE);
      cout << "done" << endl;
      cout << "    flagEnclosedPixels()...";
      img->flagEnclosedPixels(PngData::FLAG_RED,PngData::FLAG_GREEN);
      img->flagEnclosedPixels(PngData::FLAG_RED,PngData::FLAG_BLUE);
      cout << "done" << endl;
      cout << "    diffuseZoneGroups()...";
      img->diffuseZoneGroups();
      cout << "done" << endl;

      cout << "    ";
      vector<int> pixelVec;
      img->extractLineSegment(pixelVec,PngData::FLAG_YELLOW,line_start_ij[0],line_start_ij[1],true);

      // write pixel information to a file

      string logPathAndFilename = outputPath + "/farfield_points.dat";
      ofstream ofile;
      size_t fwidth = 13;
      MiscUtils::mkdir_for_file(logPathAndFilename);
      ofile.open(logPathAndFilename.c_str());
      if (!ofile.is_open()) {
        CERR("could not open file: " << logPathAndFilename);
      }
      // the parameter...
      ofile << "# " << param->str() << "\n";
      ofile << "#1:segment_id" << setw(fwidth) << "2:image-i" << setw(fwidth) << "3:image-j" << setw(fwidth) << "4:xp" << setw(fwidth) << "5:yp" << setw(fwidth) << "6:zp" << endl;

      bool b_foundFlag;
      int segment_id = 0;
      do {
        for (int ii = 0; ii < pixelVec.size(); ++ii){
          const int ipx = pixelVec[ii];
          double xp[3];
          const int ip = ipx%img->getNx();
          const int jp = ipx/img->getNx();
          img->convertPxtoXp(ip,jp,xp);
          ofile << 
            setw(fwidth) << segment_id <<
            setw(fwidth) << ip <<
            setw(fwidth) << jp <<
            setw(fwidth) << xp[0] <<
            setw(fwidth) << xp[1] <<
            setw(fwidth) << xp[2] <<
            endl;
        }
        pixelVec.clear();

        b_foundFlag=false;
        for (int ipx=0; ipx < img->getNx()*img->getNy(); ++ipx) {
          if (img->pixel_flag[ipx] == PngData::FLAG_YELLOW){
            b_foundFlag = true;
            int i_new = ipx%(img->getNx());
            int j_new = ipx/(img->getNx());
            //cout << "IJ " << i_new << " " << j_new << endl;

            cout << "    ";
            img->extractLineSegment(pixelVec,PngData::FLAG_YELLOW,i_new,j_new,true);
            ++segment_id;
            break;
          }
        }
      } while (b_foundFlag);
      ofile.close();

      //write out an image
      cout << "    writing reference image...";
      img->rescaleRangesToData();
      img->initializeWrite();

      // write image to a tmp file...
      string tmpPathAndFilename = outputPath + "/." + MiscUtils::getFilenameNoPath(pngCache->getImageName(iI));
      img->write(tmpPathAndFilename.c_str());
      img->lock = false;
      
      // then rename...
      rename(tmpPathAndFilename.c_str(),outputPathAndFilename.c_str());
      cout << "done" << endl;
      ///
       COUT1("END MARK_MEDIAL_AXIS");;

      //cout << "DEPTH CHECK" << endl;
      //int ip = 663;
      //for (int jp = 0; jp<img->getNy(); ++jp){
      //  cout << img->getDepth(jp*img->getNx()+ip) << endl;
      //}
    }

  }

  void itpTile(Param * param) {

    // ===========================================================================
    // TILE <ncol> <nrow> CROP <i0> <j0> <i1> <j1> PATH <output-path> CHECKTIME
    // ===========================================================================
    
    COUT1("TILE...");
    
    if ((!pngCache)||(pngCache->imageList.size() == 0)) {
      CERR("cannot TILE: pngCache is empty.");
    }

    if (param->size() < 2) {
      CERR("Expect TILE <ncol> <nrow> [...]");
    }
    
    const int ncol = param->getInt(0);
    const int nrow = param->getInt(1);
    
    COUT1(" > size: " << ncol << " across by " << nrow << " down");

    const int nPanel = nrow*ncol; 
    if (nPanel <= 1) {
      CERR("TILE is for multiple rows and/or columns");
    }
    
    // also, we need an even number of images...
    if (pngCache->size()%nPanel != 0) {
      CERR("cannot TILE: expecting image count to be a multiple of nrowxncol: " << nPanel);
    }
    
    const int nImg = pngCache->size()/nPanel;

    COUT1(" > number of tiled images: " << nImg);

    bool b_crop = false;
    int crop[4]; // i0,j0,i1,j1

    string outputPath = "tile";
    
    bool b_checktime = false;
    bool b_index = false; int  i_index = -1;
    
    int iarg = 2;
    while (iarg < param->size()) {
      string token = param->getString(iarg++);
      if (token == "CROP") {
        crop[0] = param->getInt(iarg++);
        crop[1] = param->getInt(iarg++);
        crop[2] = param->getInt(iarg++);
        crop[3] = param->getInt(iarg++);
        b_crop = true;
        COUT1(" > CROP: " << crop[0] << " " << crop[1] << " to " << crop[2] << " " << crop[3]);
      }
      else if (token == "PATH") {
        outputPath = param->getString(iarg++);
        COUT1(" > PATH: " << outputPath);
      }
      else if (token == "CHECKTIME") {
        b_checktime = true;
        cout << " > CHECKTIME mode: will check file times and only do new work" << endl;
      }
      else if (token == "INDEX") {
        b_index = true;
        i_index = param->getInt(iarg++);
        cout << " > INDEX mode: process only image index" << i_index << endl;
      }
      else {
        CERR("unrecognized TILE token: " << token);
      }
    }

    MiscUtils::mkdir_for_file(outputPath+"/test.png");
    
    int nx,ny;
    if (b_crop) {
      nx = crop[2]-crop[0]+1;
      ny = crop[3]-crop[1]+1;
    }
    else {
      // get the size from image 0 and assume/insist it is the same...
      PngData * img = pngCache->getImage(0,false); // false here is for minimal read
      nx = img->getNx();
      ny = img->getNy();
    }

    COUT1(" > output images will be: nx,ny: " << nx*ncol << " " << ny*nrow);

    unsigned char (*rgb)[3] = new unsigned char[nx*ncol*ny*nrow][3];

    // parse any TEXT lines...
    vector<ImageText> textVec;
    processTextParams(textVec);
    
    // we know the final image size so set the i,j from any ifrac,jfrac... 
    for (int ii = 0; ii < textVec.size(); ++ii) {
      if (textVec[ii].b_frac) {
        textVec[ii].i = int(textVec[ii].ifrac*double(nx*ncol));
        textVec[ii].j = int(textVec[ii].jfrac*double(ny*nrow));
      }
    }
    
    for (int iI = 0; iI < nImg; ++iI) {

      char outfilename[128];
      const int iIout = (b_index) ? i_index : iI;
      sprintf(outfilename,"%s/tile.%06d.png",outputPath.c_str(),iIout);
      
      cout << " > building tiled image " << iI << " of " << nImg << ": " << outfilename << "..." << endl;
      
      if (b_checktime) {
        bool b_doit = true;
        if (MiscUtils::fileExists(outfilename)) {
          b_doit = false;
          struct stat t_stat_out;
          stat(outfilename,&t_stat_out);
          int offset = 0;
          for (int irow = 0; irow < nrow; ++irow) {
            for (int icol = 0; icol < ncol; ++icol) {
              struct stat t_stat_in;
              stat(pngCache->getImageName(iI+offset*nImg).c_str(),&t_stat_in);
              // if we find any input image that is older than the output, we have to doit...
              if (t_stat_in.st_ctime >= t_stat_out.st_ctime) {
                cout << "   > time of input file: " << pngCache->getImageName(iI+offset*nImg) << " is newer: rebuilding image" << endl; 
                b_doit = true;
                break;
              }
              else {
                cout << "   > time of input file: " << pngCache->getImageName(iI+offset*nImg) << " is older" << endl; 
              }
              ++offset;
            }
            if (b_doit)
              break;
          }
        }
        if (!b_doit) {
          cout << "   > skipping. Turn off CHECKTIME mode to rebuild this tile" << endl;
          continue;
        }
      }
      
      int offset = 0;
      for (int irow = 0; irow < nrow; ++irow) {
        for (int icol = 0; icol < ncol; ++icol) {
          PngData * img = pngCache->getImage(iI+offset*nImg,false); // false here is for minimal read
          cout << "   > adding tile at row " << irow << " and col " << icol << ": " << pngCache->getImageName(iI+offset*nImg) << endl;
          const int this_nx = img->getNx();
          const int this_ny = img->getNy();
          if(this_nx <= 0 || this_ny <= 0) {
            // Insert blank tiles for non-existent images
            CWARN("      > blank tile inserted for unopened image file");
          }
          else {
            if (b_crop) {
              if (crop[0] >= crop[2] || crop[1] >= crop[3] ||
                  crop[0] <= 0       || crop[1] <= 0       ||
                  crop[2] >= this_nx || crop[3] >= this_ny  ) {
                CERR("Invalid CROP window for image");
              }
              for (int ip = crop[0]; ip <= crop[2]; ++ip) {
                assert((ip >= 0)&&(ip < this_nx));
                for (int jp = crop[1]; jp <= crop[3]; ++jp) {
                  assert((jp >= 0)&&(jp < this_ny));
                  const int ipx = jp*this_nx + ip;
                  // ipx is J*NX + I...
                  const int ipx_out = (irow*ny+(jp-crop[1]))*ncol*nx + (icol*nx+(ip-crop[0]));
                  img->getRGB(rgb[ipx_out],ipx);
                }
              }
            }
            else {
              if ((this_nx != nx)||(this_ny != ny)) {
                CERR("In current TILE implementation all images must be the same size. Try the CROP option");
              }
              for (int jp = 0; jp < ny; ++jp) {
                for (int ip = 0; ip < nx; ++ip) {
                  const int ipx = jp*nx + ip;
                  // ipx is J*NX + I...
                  const int ipx_out = (irow*ny+jp)*ncol*nx + (icol*nx+ip);
                  img->getRGB(rgb[ipx_out],ipx);
                }
              }
            }
          }
          img->lock = false;
          ++offset;
        }
      }

      // we are about to write rgb. Add text here...

      for (int ii = 0; ii < textVec.size(); ++ii) {
        cout << "   > adding centered text \"" << textVec[ii].str << "\" at " << textVec[ii].i << " " << textVec[ii].j << endl;
        addTextCentered(rgb,ncol*nx,nrow*ny,textVec[ii].str,textVec[ii].i,textVec[ii].j,textVec[ii].ifs);
      }
      
      // this constructor takes the rgb and will not delete in its destructor...
      PngImage img_out(ncol*nx,nrow*ny,rgb);

      // use a temporary filename to write...
      char tmpfilename[128];
      sprintf(tmpfilename,"%s/.tile.%06d.png",outputPath.c_str(),iIout);
      img_out.write(tmpfilename);
      
      // and rename to final output file...
      rename(tmpfilename,outfilename);

    }

    delete[] rgb;

  }//itpTile()

  void itpProbe(Param * param){

    if ((!pngCache)||(pngCache->imageList.size() == 0)) {
      CERR("cannot compute PROBE: pngCache is empty.");
    }

    double *xp = NULL;
    int *ij_probe = NULL;
    int npoints;
    bool b_xp = false;
    bool b_ij = false;

    int iarg = 0;
    string token = param->getString(iarg++);
    if (token == "GEOM") {
      token = param->getString(iarg++);
      if (token == "POINT") {
        npoints = 1;
        token = param->getString(iarg++);
        if (token == "XP") {
          b_xp = true;
          xp = new double[3];
          xp[0] = param->getDouble(iarg++);
          xp[1] = param->getDouble(iarg++);
          xp[2] = param->getDouble(iarg++);
        }
        else if (token == "IJ") {
          b_ij = true;
          ij_probe = new int[2];
          ij_probe[0] = param->getInt(iarg++);
          ij_probe[1] = param->getInt(iarg++);
        }
        else {
          CERR("Invalid POINT PROBE syntax. Sample usage:\n" <<
	       " PROBE GEOM POINT XP <x> <y> <z> [<outputfile>]\n" <<
	       " PROBE GEOM POINT IJ <pixel_i> <pixel_j> [<outputfile>]");
        }
      }
      else if (token == "LINE") {
        double x_ends[6];
        int ij_ends[4];
        bool ioerr = false;
        token = param->getString(iarg++);
        if (token == "XP") {
          b_xp = true;
          for (int i = 0; i < 6; i++)
            x_ends[i] = param->getDouble(iarg++);
        }
        else if (token == "IJ") {
          b_ij = true;
          for (int i = 0; i < 4; i++) {
            ij_ends[i] = param->getInt(iarg++);
          }
        }
        else ioerr = true;
        token = param->getString(iarg++);
        if (token == "NPOINTS") {
          npoints = param->getInt(iarg++);
          if (b_xp) {
            xp = new double[npoints*3];
            for (int i = 0; i < npoints; i++) {
              for (int j = 0; j < 3; j++) {
                xp[3*i+j] = x_ends[j] + (double)(i-0)/(double)(npoints - 1)*(x_ends[3+j] - x_ends[j]);
              }
            }
          }
          else if (b_ij) {
            ij_probe = new int[npoints*2];
            for (int i = 0; i < npoints; i++) {
              ij_probe[2*i+0] = ij_ends[0] + (int)(double(i - 0)/double(npoints - 1)*double(ij_ends[2] - ij_ends[0]));
	      ij_probe[2*i+1] = ij_ends[1] + (int)(double(i - 0)/double(npoints - 1)*double(ij_ends[3] - ij_ends[1]));
            }
          }
        }
        else ioerr = true;
        if (ioerr == true) {
          CERR("Invalid LINE PROBE syntax. Sample usage:\n" <<
	       " PROBE GEOM LINE XP <x0> <y0> <z0> <x1> <y1> <z1> NPOINTS <npoints> [<outputfile>]\n" <<
	       " PROBE GEOM LINE IJ <pixel_i0> <pixel_j0> <pixel_i1> <pixel_j1> NPOINTS <npoints> [<outputfile>]");
        }
      }
      else CERR("UNRECOGNIZED PROBE GEOMETRY. \nCurrent allowed values are POINT and LINE\n");
    }
    else {
      CERR("MUST PROVIDE PROBE GEOMETRY.");
    }

    string outputName = "probe_point.dat";
    if (iarg < param->size())
      outputName = param->getString(iarg++);

    cout << " > writing points to file \"" << outputName << "\"..." << endl;

    ofstream ofile;
    size_t fwidth = 13;
    MiscUtils::mkdir_for_file(outputName);
    ofile.open(outputName.c_str());
    if (!ofile.is_open()) {
      CERR("could not open file: " << outputName);
    }

    // the parameter...
    ofile << "# " << param->str() << "\n";

    double * var = new double[npoints*pngCache->size()];
    int * pixel_type = new int[npoints*pngCache->size()];
    double * time = NULL;
    if (npoints == 1)
      time = new double[pngCache->size()];

    double * xp_out = new double[pngCache->size()*npoints*3];
    int * ij_probe_out = new int[pngCache->size()*npoints*2];

    for (int iI = 0; iI<pngCache->size(); ++iI){
      cout << " > working on image: " << iI << " of " << pngCache->size() << "...";
      PngData * img = pngCache->getImage(iI);
      if (b_xp){
        for (int i = 0; i < npoints; i++) {
          pixel_type[iI*npoints+i] = img->probePoint(xp+i*3,xp_out+3*(iI*npoints+i),ij_probe_out+2*(iI*npoints+i),var[iI*npoints+i]);
        }
      }
      else if (b_ij) {
        for (int i = 0; i < npoints; i++) {
          pixel_type[iI*npoints+i] = img->convertPxtoXp(ij_probe[i*2],ij_probe[i*2+1],xp_out+3*(iI*npoints+i));
          img->probeIJ(ij_probe[i*2], ij_probe[i*2+1], var[iI*npoints+i]);
          ij_probe_out[2*(iI*npoints+i)  ] = ij_probe[i*2  ];
          ij_probe_out[2*(iI*npoints+i)+1] = ij_probe[i*2+1];
        }
      }
      else{
        assert(0);
      }

      // for the case of one point, report the value to stdout...
      if (npoints == 1){
        time[iI] = img->getTime();
        cout << " " << var[iI];
      }
      cout << endl;
      img->lock = false;
    }
    // for single points, we report the time trace of points through the images
    // for multiple points, we assume the user wants x,y,z,stats in columns.
    if (npoints == 1) {
      // when writing just one point, we want the time history of the point as a column
      // single-point header...
      ofile << "# Pixel Type Legend " << getPixelTypeLegendStr() << endl;
      ofile << "#" << setw(fwidth-5) << "1:image" << setw(fwidth) << "2:time" << setw(fwidth) << "3:pixel_type" << setw(fwidth) << "4:image-i" << setw(fwidth) << "5:image-j" << setw(fwidth) << "6:xp" << setw(fwidth) << "7:yp" << setw(fwidth) << "8:zp" << setw(fwidth) << "9:var" << setw(fwidth) << "10:avg" << setw(fwidth) << "11:rms" << endl;
      double sum = 0.0;
      double sum2 = 0.0;
      double sum_wgt = 0.0;
      for (int iI = 0; iI<pngCache->size(); ++iI){
        sum += var[iI];
        sum2 += var[iI]*var[iI];
        sum_wgt += 1.0;
        ofile <<
          setw(fwidth-4) << iI <<
          setw(fwidth) << time[iI] <<
          setw(fwidth) << pixel_type[iI] <<
          setw(fwidth) << ij_probe_out[iI*2  ] <<
          setw(fwidth) << ij_probe_out[iI*2+1] <<
          setw(fwidth) << xp_out[iI*3  ] <<
          setw(fwidth) << xp_out[iI*3+1] <<
          setw(fwidth) << xp_out[iI*3+2] <<
          setw(fwidth) << var[iI] <<
          setw(fwidth) << sum/sum_wgt <<
          setw(fwidth) << sqrt(max(0.0,sum2/sum_wgt - sum/sum_wgt*sum/sum_wgt)) <<
          endl;
      }
      delete[] time;
      time = NULL;
    }
    else {
      // multi-point header...
      if (b_xp){
        ofile << "# Probing by constant simulation coordinate, image-ij reflects pixel location in image-0" << endl;
      }
      else if (b_ij){
        ofile << "# Probing by constant image pixel ij, simulation coordinates reflect location in image-0" << endl;
      }
      else{
        assert(0);
      }
      ofile << "# Pixel Type Legend " << getPixelTypeLegendStr() << " 99:Mixed" << endl;
      ofile << "#" << setw(fwidth-5) << "1:npt" << setw(fwidth) << "2:pixel_type" << setw(fwidth) << "3:image-i" << setw(fwidth) << "4:image-j" << setw(fwidth) << "5:xp" << setw(fwidth) << "6:yp" << setw(fwidth) << "7:zp" << setw(fwidth) << "8:avg" << setw(fwidth) << "9:rms" << setw(fwidth) << "10:avg(first-half)" << setw(fwidth) << "11:rms(first-half)" << setw(fwidth) << "12:avg(middle-half)" << setw(fwidth) << "13:rms(middle-half)" << setw(fwidth) << "14:avg(last-half)" << setw(fwidth) << "15:rms(last-half)";
      //ofile << "# 1:npt 2:pixel_type 3:xp 4:yp 5:zp 6:avg 7:rms 8:avg(first-half) 9:rms(first-half) 10:avg(middle-half) 11:rms(middle-half) 12:avg(last-half) 13:rms(last-half)";
      for (int iI = 0; iI<pngCache->size(); ++iI) {
        ofile << setw(fwidth) << 16+iI << ":image-" << iI;
      }
      ofile << endl;
      for (int i = 0; i < npoints; i++) {
        double avg[4]  = { 0.0, 0.0, 0.0, 0.0 };
        double rms[4]   = { 0.0, 0.0, 0.0, 0.0 };
        double count[4] = { 0.0, 0.0, 0.0, 0.0 };
        int ptype0 = pixel_type[i];
        int ptype = ptype0;
        for (int iI = 0; iI<pngCache->size(); ++iI){
          if (ptype0!=pixel_type[iI*npoints+i]){
            ptype = 99;
            cout << "Warning, image pixel data type at point " << i << " varies across image set" << endl;
          }
          // [0] stores the stats for the full set of images...
          avg[0]  += var[iI*npoints+i];
          rms[0]   += var[iI*npoints+i]*var[iI*npoints+i];
          count[0] += 1.0;
          // we also store the avg/rms over ranges to test statistical steadiness of the images...
          if (iI < pngCache->size()/2) {
            // [1] is the first half...
            avg[1]  += var[iI*npoints+i];
            rms[1]   += var[iI*npoints+i]*var[iI*npoints+i];
            count[1] += 1.0;
          }
          else {
            // [3] is the last half...
            avg[3]  += var[iI*npoints+i];
            rms[3]   += var[iI*npoints+i]*var[iI*npoints+i];
            count[3] += 1.0;
          }
          if ((iI >= pngCache->size()/4)&&(iI < pngCache->size()*3/4))  {
            // [2] is the middle half (note overlaps others)...
            avg[2]  += var[iI*npoints+i];
            rms[2]   += var[iI*npoints+i]*var[iI*npoints+i];
            count[2] += 1.0;
          }
        }
        ofile <<
          setw(fwidth) << i <<
          setw(fwidth) << ptype << //check earlier that this is the same across images, 99 if not
          setw(fwidth) << ij_probe_out[2*i  ] <<
          setw(fwidth) << ij_probe_out[2*i+1] <<
          setw(fwidth) << xp_out[3*i  ] <<
          setw(fwidth) << xp_out[3*i+1] <<
          setw(fwidth) << xp_out[3*i+2];
        FOR_J4 {
          if (count[j] > 0.0) {
            avg[j] /= count[j];
            rms[j] = sqrt(max(0.0,rms[j]/count[j] - avg[j]*avg[j]));
          }
          ofile <<
            setw(fwidth) << avg[j] <<
            setw(fwidth) << rms[j];
        }
        for (int iI = 0; iI<pngCache->size(); ++iI) {
          ofile <<
            setw(fwidth) << var[iI*npoints+i];
        }
        ofile << endl;
      }
    }
    delete[] var;
    delete[] pixel_type;
    ofile.close();

    if (xp != NULL)           delete[] xp;
    if (ij_probe != NULL)     delete[] ij_probe;
    if (time !=  NULL)        delete[] time;
    if (xp_out!=NULL)         delete[] xp_out;
    if (ij_probe_out != NULL) delete[] ij_probe_out;

  }//itpProbe()

  void processNprobe(Param * param){

    if ((!pngCache)||(pngCache->imageList.size() == 0)) {
      CERR("cannot compute PROBE: pngCache is empty.");
    }

    vector<double> xVec,yVec,zVec,targetVec;
    vector<int> binVec;
    bool b_with_target = false;
    bool b_with_bin = false;
    int bin_max = 0;

    int iarg = 0;
    while (iarg < param->size()) {
      string token = param->getString(iarg++);
      if (token == "WITH_TARGET") {
	b_with_target = true;
      }
      else if (token == "WITH_BIN") {
	b_with_bin = true;
      }
      else {
	// assuming x,y,z and target if WITH_TARGET
	xVec.push_back(param->getDouble(iarg-1));
	yVec.push_back(param->getDouble(iarg++));
	zVec.push_back(param->getDouble(iarg++));
	if (b_with_bin) {
	  const int bin = param->getInt(iarg++);
	  binVec.push_back(bin);
	  bin_max = max(bin_max,bin);
	}
	if (b_with_target) targetVec.push_back(param->getDouble(iarg++));
      }
    }

    const int np = xVec.size();
    cout << " > got np=" << np << endl;

    if (b_with_target) {
      double avg = 0.0;
      double rms = 0.0;
      for (int ip = 0; ip < np; ++ip) {
	avg += targetVec[ip];
	rms += targetVec[ip]*targetVec[ip];
      }
      avg /= double(targetVec.size());
      rms = sqrt(max(0.0,rms/double(targetVec.size())-avg*avg));
      cout << " > target avg=" << avg << " rms=" << rms << endl;
    }

    double * var = new double[np*pngCache->size()];
    for (int iI = 0; iI<pngCache->size(); ++iI) {
      cout << " > working on image: " << iI << " of " << pngCache->size() << "..." << endl;
      PngData * img = pngCache->getImage(iI);
      for (int ip = 0; ip < np; ++ip) {
	double xp[3] = { xVec[ip], yVec[ip], zVec[ip] };
	int ij[2];
        double xout[3];
	img->probePoint(xp,xout,ij,var[iI*np+ip]);
      }
      img->lock = false;
    }

    double * bin_sum = NULL;
    int * bin_count = NULL;
    if (b_with_bin) {
      bin_sum = new double[bin_max+1];
      bin_count = new int[bin_max+1];
      for (int bin = 0; bin <= bin_max; ++bin) {
	bin_sum[bin] = 0.0;
	bin_count[bin] = 0;
      }
    }

    double dtarget = 0.0;
    double dtarget2 = 0.0;
    for (int ip = 0; ip < np; ++ip) {
      double avg = 0.0;
      double rms  = 0.0;
      double count = 0.0;
      for (int iI = 0; iI < pngCache->size(); ++iI){
	avg += var[iI*np+ip];
	rms  += var[iI*np+ip]*var[iI*np+ip];
	count += 1.0;
	// and bin...
	if (b_with_bin) {
	  //bin_sum[binVec[ip]] += targetVec[ip];
	  bin_sum[binVec[ip]] += var[iI*np+ip];
	  bin_count[binVec[ip]] += 1;
	}
      }
      // the avg value for all images at this probe...
      avg /= count;
      cout << "NPROBE point " << ip << " avg: " << avg;
      if (b_with_target) {
	dtarget += avg - targetVec[ip];
	dtarget2 += (avg - targetVec[ip])*(avg - targetVec[ip]);
	cout << " target: " << targetVec[ip] << " error: " << avg-targetVec[ip];
      }
      cout << endl;
    }

    cout << "L2 error: " << sqrt(dtarget2/double(np)) <<
      " avg error: " << dtarget/double(np) <<
      " rms error: "  << sqrt(max(0.0,dtarget2/double(np) - dtarget/double(np)*dtarget/double(np))) << endl;

    if (b_with_bin) {
      cout << "bin averages: " << endl;
      for (int bin = 0; bin <= bin_max; ++bin) {
	if (bin_count[bin] > 0) {
	  cout << bin << " " << bin_sum[bin]/double(bin_count[bin]) << endl;
	}
      }
    }

  }//processNprobe()


  void itpProjectedSurfaceArea(Param * param){
    COUT1("PROJECTED_SURFACE_AREA");

    if ((!pngCache)||(pngCache->imageList.size() == 0)) {
      CERR("cannot compute PROJECTED_SURFACE_AREA: pngCache is empty.");
    }

    bool b_markSurfaces = false;
    if (param->size()>0&&(param->getString(0)=="WRITE"))
      b_markSurfaces = true;

    for (int iI = 0; iI<pngCache->size(); ++iI){
      PngData * img = pngCache->getImage(iI);

      double surface_area = img->computeSurfaceArea(b_markSurfaces);

      cout << " > Image \"" << pngCache->imageList[iI] << "\"";
      if (b_markSurfaces){
        string filename = buildFilenameFromImage("proj_area_green.png",pngCache->imageList[iI]);
        img->rescaleRangesToData();
        img->initializeWrite();
        img->write(filename.c_str());
        cout << ", flagged surface written to \"" << filename << "\"";
      }
      cout << ", projected surface area = " << surface_area << endl;
      img->lock = false;
    }

  }//itpProjectedSurfaceArea

  void itpCollapseHeight(Param * param){

    COUT1("COLLAPSE_HEIGHT...");

    if ((!pngCache)||(pngCache->imageList.size() == 0)) {
      CERR("cannot COLLAPSE_HEIGHT: pngCache is empty.");
    }

    int pixelType = PngData::BACKGROUND_PIXEL;

    bool b_between = false;
    int j_f = 0;
    int j_l = -1;

    int iarg = 0;
    while (iarg < param->size()) {
      const string token = param->getString(iarg++);
      if (token == "PIXEL_TYPE") {
        //If specified, limit output to a single pixel type
        pixelType = param->getInt(iarg++);
        cout << " > Including only pixels of type " << pixelType << endl;
      }
      else if (token == "BETWEEN"){
        b_between = true;
        j_f = param->getInt(iarg++);
        j_l = param->getInt(iarg++);
        if (j_f > j_l || j_f < 0){
          CERR("Row bounds in pixels must be greater than or equal to zero with the smaller value first");
        }
        cout << "limiting collapse to rows between " << j_f << " and " << j_l << ", inclusive" << endl;
      }
    }

    ofstream ofile;
    size_t fwidth = 13;

    for (int iI = 0; iI<pngCache->size(); ++iI){

      cout << " > working on image \"" << pngCache->imageList[iI] << "\"" << endl;

      string outputDat = buildFilenameFromImage("collapse_height",pngCache->imageList[iI]) + ".dat";
      ofile.open(outputDat.c_str());
      assert(ofile.is_open());

      ofile << "# " << param->str() << endl;
      ofile << "# Pixel Type Legend " << getPixelTypeLegendStr() << " 99:Mixed" << endl;
      ofile << setw(fwidth-2) << "# 1:i-pixel" <<
        setw(fwidth-5) << "2:type" <<
        setw(fwidth) << "3:x-avg" <<
        setw(fwidth) << "4:y-avg" <<
        setw(fwidth) << "5:z-avg" <<
        setw(fwidth) << "6:average" <<
        setw(fwidth) << "7:min" <<
        setw(fwidth) << "8:max" << endl;

      PngData * img = pngCache->getImage(iI);

      if (b_between){
        if(j_l>=img->getNy()){
          cout << "Warning: max row bound " <<  j_l << " out of range, setting to " << img->getNy()-1 << endl;
          j_l = img->getNy()-1;
        }
      }
      else{
        j_f=0;
        j_l=img->getNy()-1;
      }

      for (int ip = 0; ip < img->getNx(); ++ip) {
        double avg = 0.0;
        int count = 0;
        double minval =  HUGE_VAL;
        double maxval = -HUGE_VAL;
        bool b_first = true;
        int pixelTypeOut = -1;//init w nonsense
        double xp_sum[3] = {0.0,0.0,0.0};
        for (int jp = j_f; jp <= j_l; ++jp) {
          const int ipx = jp*img->getNx() + ip;
          if ((pixelType!=PngData::BACKGROUND_PIXEL && img->pixel_flag[ipx] == pixelType) ||
              (pixelType==PngData::BACKGROUND_PIXEL && img->pixel_flag[ipx] >= PngData::VOLUME_DATA_PIXEL)) {
            avg += img->pixel_data[ipx];
            minval = min(minval,img->pixel_data[ipx]);
            maxval = max(maxval,img->pixel_data[ipx]);
            double xp[3];
            img->convertPxtoXp(ip,jp,xp);
            FOR_I3 xp_sum[i] += xp[i];
            ++count;
            if (b_first){
              b_first = false;
              pixelTypeOut=img->pixel_flag[ipx];
            }
            else if (pixelTypeOut!=img->pixel_flag[ipx]){
              pixelTypeOut=99;
            }
          }
        }
        if (count > 0) {
          ofile << setw(fwidth-2) << ip <<
            setw(fwidth-5) << pixelTypeOut <<
	    setw(fwidth) << xp_sum[0]/double(count) <<
	    setw(fwidth) << xp_sum[1]/double(count) <<
	    setw(fwidth) << xp_sum[2]/double(count) <<
            setw(fwidth) << avg/double(count) <<
            setw(fwidth) << minval <<
            setw(fwidth) << maxval << endl;
        }
      }

      ofile.close();

      img->lock = false;

    }

  }//itpCollapseHeight()


  void itpCollapseWidth(Param * param){

    COUT1("COLLAPSE_WIDTH...");

    if ((!pngCache)||(pngCache->imageList.size() == 0)) {
      CERR("cannot COLLAPSE_WIDTH: pngCache is empty.");
    }

    int pixelType = PngData::BACKGROUND_PIXEL;

    bool b_between = false;
    int i_f = 0;
    int i_l = -1;

    int iarg = 0;
    while (iarg < param->size()) {
      const string token = param->getString(iarg++);
      if (token == "PIXEL_TYPE") {
        //If specified, limit output to a single pixel type
        pixelType = param->getInt(iarg++);
        cout << " > Including only pixels of type " << pixelType << endl;
      }
      else if (token == "BETWEEN"){
        b_between = true;
        i_f = param->getInt(iarg++);
        i_l = param->getInt(iarg++);
        if (i_f > i_l || i_f < 0){
          CERR("Column bounds in pixels must be greater than or equal to zero with the smaller value first");
        }
        cout << "limiting collapse to columns between " << i_f << " and " << i_l << ", inclusive" << endl;
      }
    }

    ofstream ofile;
    size_t fwidth = 13;

    for (int iI = 0; iI<pngCache->size(); ++iI){

      cout << " > working on image \"" << pngCache->imageList[iI] << "\"" << endl;

      string outputDat = buildFilenameFromImage("collapse_width",pngCache->imageList[iI]) + ".dat";
      ofile.open(outputDat.c_str());
      assert(ofile.is_open());

      ofile << "# " << param->str() << endl;
      ofile << "# Pixel Type Legend " << getPixelTypeLegendStr() << " 99:Mixed" << endl;
      ofile << setw(fwidth-2) << "# 1:j-pixel" <<
        setw(fwidth-5) << "2:type" <<
        setw(fwidth) << "3:x-avg" <<
        setw(fwidth) << "4:y-avg" <<
        setw(fwidth) << "5:z-avg" <<
        setw(fwidth) << "6:average" <<
        setw(fwidth) << "7:min" <<
        setw(fwidth) << "8:max" << endl;

      PngData * img = pngCache->getImage(iI);

      if (b_between){
        if(i_l>=img->getNx()){
          cout << "Warning: max column bound " <<  i_l << " out of range, setting to " << img->getNx()-1 << endl;
          i_l = img->getNx()-1;
        }
      }
      else{
        i_f=0;
        i_l=img->getNx()-1;
      }

      for (int jp = 0; jp < img->getNy(); ++jp) {
	double avg = 0.0;
	int count = 0;
	double minval =  HUGE_VAL;
	double maxval = -HUGE_VAL;
        bool b_first = true;
        int pixelTypeOut = -1;// init w nonsense
        double xp_sum[3] = {0.0,0.0,0.0};
	for (int ip = i_f; ip <= i_l; ++ip) {
	  const int ipx = jp*img->getNx() + ip;
	  if ((pixelType!=PngData::BACKGROUND_PIXEL && img->pixel_flag[ipx] == pixelType) ||
              (pixelType==PngData::BACKGROUND_PIXEL && img->pixel_flag[ipx] >= PngData::VOLUME_DATA_PIXEL)) {
	    avg += img->pixel_data[ipx];
	    minval = min(minval,img->pixel_data[ipx]);
	    maxval = max(maxval,img->pixel_data[ipx]);
            double xp[3];
            img->convertPxtoXp(ip,jp,xp);
            FOR_I3 xp_sum[i] += xp[i];
	    ++count;
            if (b_first){
              b_first = false;
              pixelTypeOut=img->pixel_flag[ipx];
            }
            else if (pixelTypeOut!=img->pixel_flag[ipx]){
              pixelTypeOut=99;
            }
	  }
	}
	if (count > 0) {
	  ofile << setw(fwidth-2) << jp <<
            setw(fwidth-5) << pixelTypeOut <<
	    setw(fwidth) << xp_sum[0]/double(count) <<
	    setw(fwidth) << xp_sum[1]/double(count) <<
	    setw(fwidth) << xp_sum[2]/double(count) <<
	    setw(fwidth) << avg/double(count) <<
	    setw(fwidth) << minval <<
	    setw(fwidth) << maxval << endl;
	}
      }

      ofile.close();

      img->lock = false;

    }

  }//itpCollapseWidth()


  // JOB Dump pixels to data utility
  void itpDumpData(Param * param) {
    if ((!pngCache)||(pngCache->imageList.size() == 0)) {
      CERR("cannot DUMP_DATA: pngCache is empty.");
    }

    ofstream ofile;
    for (int iI = 0; iI<pngCache->size(); ++iI) {
      cout << " working on image \"" << pngCache->imageList[iI] << "\"" << endl;

      string outputDat = buildFilenameFromImage("raw_values",pngCache->imageList[iI]) + ".dat";
      ofile.open(outputDat.c_str());
      assert(ofile.is_open());

      size_t fwidth = 13;
      ofile << "# Pixel Type Legend " << getPixelTypeLegendStr() << endl;
      ofile << "#"
            << setw(fwidth) << "1:count"
            << setw(fwidth) << "2:pixel_type"
            << setw(fwidth) << "3:x-pixel"
            << setw(fwidth) << "4:y-pixel"
            << setw(fwidth) << "5:x_physical"
            << setw(fwidth) << "6:y_physical"
            << setw(fwidth) << "7:z_physical"
            << setw(fwidth) << "8:value" << endl;

      PngData * img = pngCache->getImage(iI);

      double xp[3];
      int count = 0;
      for (int ip = 0; ip < img->getNx(); ++ip) {
        for (int jp = 0; jp < img->getNy(); ++jp) {
          const int ipx = jp*img->getNx() +ip;
          if (img->pixel_flag[ipx] >= PngData::VOLUME_DATA_PIXEL){
            img->convertPxtoXp(ip,jp,xp);
            ofile << " " << setw(fwidth) << count++
                  << setw(fwidth) << img->pixel_flag[ipx]
                  << setw(fwidth) << ip
                  << setw(fwidth) << jp
                  << setw(fwidth) << xp[0]
                  << setw(fwidth) << xp[1]
                  << setw(fwidth) << xp[2]
                  << setw(fwidth) << img->pixel_data[ipx] << "\n";
          }
        }//for jp
      }//for ip
      ofile.close();
      img->lock = false;
    }//for iI
  }//itpDumpData()

  void itpDumpColumn(Param * param) {
    if ((!pngCache)||(pngCache->imageList.size() == 0)) {
      CERR("cannot DUMP_COLUMN: pngCache is empty.");
    }

    int theColumn;
    if (param->size()>0){
      theColumn = param->getInt(0);
    }
    else{
      CERR("cannot DUMP_COLUMN: please specify a column index");
    }

    ofstream ofile;
    for (int iI = 0; iI<pngCache->size(); ++iI) {

      PngData * img = pngCache->getImage(iI);
      if (theColumn>=0&&theColumn<img->getNx())
        cout << " working on image \"" << pngCache->imageList[iI] << "\"" << endl;
      else
        cout << " skipping image, column out of range, valid range: 0 " << img->getNx()-1 << endl;

      stringstream ss;
      ss << "column_" << theColumn;
      string outputDat = "./" + buildFilenameFromImage(ss.str(),pngCache->imageList[iI]) + ".dat";
      ofile.open(outputDat.c_str());
      assert(ofile.is_open());

      size_t fwidth = 13;
      ofile << "#"
	    << setw(fwidth) << "1:count"
	    << setw(fwidth) << "2:x-pixel"
	    << setw(fwidth) << "3:y-pixel"
	    << setw(fwidth) << "4:x_physical"
	    << setw(fwidth) << "5:y_physical"
	    << setw(fwidth) << "6:z_physical"
	    << setw(fwidth) << "7:value" << endl;

      double xp[3];
      int count = 0;
      for (int jp = 0; jp < img->getNy(); ++jp) {
        const int ipx = jp*img->getNx() +theColumn;
        if (img->pixel_flag[ipx] >= 0) {
          ++count;
          img->convertPxtoXp(theColumn,jp,xp);
          ofile
            << setw(fwidth) << count
            << setw(fwidth) << theColumn
            << setw(fwidth) << jp
            << setw(fwidth) << xp[0]
            << setw(fwidth) << xp[1]
            << setw(fwidth) << xp[2]
            << setw(fwidth) << img->pixel_data[ipx] << endl;
        }
      }//for jp
      ofile.close();
      img->lock = false;
    }//for iI
  }//itpDumpColumn()


  // 1D pdf utility
  void itp1DPdf(Param * param) {

    // get params
    int iarg = 0;
    double x_min = param->getDouble(iarg++);
    double x_max = param->getDouble(iarg++);
    int    nbins = param->getInt(iarg++);
    const double x_small = abs(x_max-x_min)*1e-15;

    string scale = "LINEAR";
    bool x_logscale = false;

    if (iarg < param->size()) {
      string token = param->getString(iarg++);
      if ((token == "LOG") || (token == "LOGSCALE")) {
        // assume MIN/MAX are same convention as pixel values (i.e. not logscale)
        x_min = log10(max(x_min,x_small));
        x_max = log10(max(x_max,x_small));
        x_logscale = true;
        scale = "LOG";
      }
    }
    double x_delta = (x_max - x_min)/double(nbins);
    cout << "\n building 1D PDF... " << endl;
    cout << " > xmin  = " << x_min << endl;
    cout << " > xmax  = " << x_max << endl;
    cout << " > nbins = " << nbins << endl;
    cout << " > scale = " << scale << endl << endl;

    // bin the data

    long int n_bg  = 0; // background
    long int n_or  = 0; // out-of-range
    long int n_ok  = 0; // valid
    long int n_tot = 0; // total

    long int * bincount = new long int [nbins];
    for (int i = 0; i < nbins; ++i) bincount[i] = 0;

    if ((!pngCache)||(pngCache->imageList.size() == 0))
      CERR("cannot compute 1D_PDF: pngCache is empty.");

    for (int iI = 0; iI < pngCache->size(); ++iI) {
      cout << " reading image: " << pngCache->imageList[iI];
      PngData * img = pngCache->getImage(iI);
      cout << " --> size = (" << img->getNx() << " x " << img->getNy() <<
        "), range = [" << img->getRangeMin() << ":" << img->getRangeMax() << "]" << endl;

      for (int ip = 0; ip < img->getNx(); ++ip) {
        for (int jp = 0; jp < img->getNy(); ++jp) {
          const int ipx = jp*img->getNx() +ip;
          n_tot++;
          if (img->pixel_flag[ipx] > PngData::SURFACE_PIXEL) {
            double val;
            if (x_logscale)
              val = log10(max(img->pixel_data[ipx],x_small));
            else
              val = img->pixel_data[ipx];

            if ((val >= x_min) && (val < x_max)) {//Pixel is on [xmin, xmax)
              int mybin = int(floor((val - x_min)/x_delta));
              assert(mybin > -1); assert(mybin < nbins);
              bincount[mybin]++;
              n_ok++; 
            }//
            //JOB: The last bin includes the right edge of the interval...
            else if (val == x_max) {//pixel equals xmax
              bincount[nbins-1]++;
              n_ok++;
            }//(right edge)
            else // pixel value is out of range
              n_or++;
          }//(contains data)
          else n_bg++;
        }//jp
      }//ip
      img->lock = false;
    }//iI

    // dump global histogram
    ofstream ofile;
    ofile.open("histogram.dat");
    assert(ofile.is_open());
    ofile << "# 1D PDF: xmin = " << x_min << ", xmax = " << x_max << ", nbins = " << nbins << ", scale = " << scale << endl;
    ofile << "# pixel breakdown..." << endl;
    ofile << "# > background   : "  << n_bg  << " (" << double(n_bg)/double(n_tot)*100.0 << "%)" << endl;
    ofile << "# > out-of-range : "  << n_or  << " (" << double(n_or)/double(n_tot)*100.0 << "%)" << endl;
    ofile << "# > valid        : "  << n_ok  << " (" << double(n_ok)/double(n_tot)*100.0 << "%)" << endl;
    ofile << "# > total        : "  << n_tot << endl;
    ofile << "# 1:bin, 2:x(center), 3:pdf, 4:cdf" << endl;
    double cdf = 0.0;
    for (int i=0; i<nbins; ++i) {
      cdf += double(bincount[i])/double(n_ok);
      double xc = x_min + (double(i)+0.5)/double(nbins)*(x_max-x_min);
      ofile << i << "  " << xc << " " << bincount[i]/double(n_ok) << " " << cdf << endl;
    }
    ofile.close();
    delete[] bincount;
  }

  // JOB clunky item to allow for two sets of images to enable 2D pdf generation
  void buildTwoPngCaches(string * imageListName){
    pngCaches    = new PngCache*[2];
    cout << "imageListName[0] = " << imageListName[0] << "\n";
    cout << "imageListName[1] = " << imageListName[1] << "\n";
    for (int i = 0; i < 2; i++) {
      pngCaches[i] = new PngCache(imageListName[i]);
      if (pngCaches[i]->imageList.empty())
        CERR("The pngCache for series " << i << " is empty. Check format of list file \"" << imageListName[i] << "\"");
    }//for (i < 2)
  }

  // JOB Diff two series of images
  void itpDiffSet(Param * param) {
    string list[2];
    string otptprefix;
    int iarg = 0;

    list[0] = param->getString(iarg++);
    list[1] = param->getString(iarg++);
    otptprefix = param->getString(iarg++);
    const double range_min = param->getFloat(iarg++);
    const double range_max = param->getFloat(iarg++);
    double data_range[2];
    data_range[0] = range_min;
    data_range[1] = range_max;

    buildTwoPngCaches(list);

    for (int iI = 0; iI<pngCaches[0]->size(); ++iI) {
      string img1name = pngCaches[0]->imageList[iI];
      string img2name = pngCaches[1]->imageList[iI];
      string otptname;
      if (iI < 10) otptname = otptprefix + ".00" +  NumberToString(iI) + ".diff.png";
      else if (iI < 100) otptname = otptprefix + ".0" +  NumberToString(iI) + ".diff.png";
      else otptname = otptprefix + "." +  NumberToString(iI) + ".diff.png";
      cout << img1name << "\n" << img2name << "\n" << otptname << "\n";

      PngData * img1 = pngCaches[0]->getImage(iI);
      PngData * img2 = pngCaches[1]->getImage(iI);
      img1->diff2(img2, data_range);

      img2->setRange(data_range);
      img2->initializeWrite();
      img2->write(otptname.c_str());
      img1->lock = false; //release all images from cache
      img2->lock = false;
    }//(iI < nImg)
  }//itpDiffSet()

  // JOB 2D Pdf generator - currently just writes to text, should write to png or at least binary?
  void itp2DPdf(Param * param) {
    //  Allocation + initialization
    bool x_logscale;
    bool y_logscale;
    string token;
    string list[2];
    long int **bincount;
    int iarg = 0;
    ofstream ofile;

    //  Set binning params
    list[0] = param->getString(iarg++);
    double x_min = param->getDouble(iarg++);
    double x_max = param->getDouble(iarg++);
    int  x_nbins = param->getInt(iarg++);
    double x_delta = (x_max - x_min)/(double)(x_nbins);
    token = param->getString(iarg++);
    if ((token == "LINEAR") || (token == "LINSCALE") || (token == "LINEARSCALE"))  x_logscale = false;
    else if ((token == "LOG") || (token == "LOGSCALE"))                           x_logscale = true;
    else {
      CERR("NEED TO SPECIFY WHETHER 2D PDF SCALE IS LINEAR OR LOG")
        throw(-2);
    }
    list[1] = param->getString(iarg++);
    double y_min = param->getDouble(iarg++);
    double y_max = param->getDouble(iarg++);
    int  y_nbins = param->getInt(iarg++);
    double y_delta = (y_max - y_min)/(double)(y_nbins);
    token = param->getString(iarg++);
    if ((token == "LINEAR") || (token == "LINSCALE") || (token == "LINEARSCALE")) y_logscale = false;
    else if ((token == "LOG") || (token == "LOGSCALE"))                           y_logscale = true;
    else {
      CERR("NEED TO SPECIFY WHETHER 2D PDF SCALE IS LINEAR OR LOG")
        throw(-2);
    }
    string outprefix = param->getString(iarg++);

    buildTwoPngCaches(list);

    cout << "Building 2D Histogram with x_min = " << x_min << ", x_max = " << x_max << "and nbins = " << x_nbins << "\n";
    cout << "Building 2D Histogram with y_min = " << y_min << ", y_max = " << y_max << "and nbins = " << y_nbins << "\n";

    bincount = new long int*[x_nbins];
    for (int i = 0; i < x_nbins; ++i) bincount[i] = new long int[y_nbins];

    //  Bin the data
    for (int iI = 0; iI<pngCaches[0]->size(); ++iI) {
      for (int i = 0; i < x_nbins; ++i)
        for (int j = 0; j < y_nbins; ++j)
          bincount[i][j] = 0;
      cout << " working on images \"" << pngCaches[0]->imageList[iI] << " and " << pngCaches[1]->imageList[iI] << endl;
      string file1 = pngCaches[0]->imageList[iI];
      string file2 = pngCaches[1]->imageList[iI];
      string outfile;
      if (iI < 10) outfile = outprefix + ".00" + NumberToString(iI) + ".2Dpdf";
      else if (iI < 100) outfile = outprefix + ".0" + NumberToString(iI) + ".2Dpdf";
      else outfile = outprefix + "." + NumberToString(iI) + ".2Dpdf";
      cout << outfile << "\n";
      ofile.open(outfile.c_str());
      assert(ofile.is_open());
      PngData * img1 = pngCaches[0]->getImage(iI);
      PngData * img2 = pngCaches[1]->getImage(iI);
      cout << "Field range is X:" << img1->getRangeMin() << " to " << img1->getRangeMax() << ", Y: " << img2->getRangeMin() <<" to " << img2->getRangeMax() << "\n";
      if (!x_logscale) {
        if ((img1 -> getRangeMin() < x_min) || (img1->getRangeMax() > x_max))
          cout << "WARNING: Some pixel values lie outside the 2D Pdf x-range!\n";
      }
      else {
        if ((log10(img1->getRangeMin()) < x_min) || (log10(img1->getRangeMax()) > x_max))
          cout << "WARNING: Some pixel values lie outside the 2D Pdf x-range!\n";
      }
      if (!y_logscale) {
        if ((img1 -> getRangeMin() < y_min) || (img1->getRangeMax() > y_max))
          cout << "WARNING: Some pixel values lie outside the 2D Pdf y-range!\n";
      }
      else {
        if ((log10(img1->getRangeMin()) < y_min) || (log10(img1->getRangeMax()) > y_max))
          cout << "WARNING: Some pixel values lie outside the 2D Pdf y-range!\n";
      }
      int n_rejects = 0;
      int n_bb = 0;
      int n_good = 0;
      int n_total = 0;
      int n_x_oob = 0;
      int n_y_oob = 0;
      for (int ip = 0; ip < img1->getNx(); ++ip) {
        for (int jp = 0; jp < img1->getNy(); ++jp) {
          n_total++;
          const int ipx = jp*img1->getNx() +ip;
          // Remove Bahama blue
          if ((img1->pixel_flag[ipx] >= 0) && (img2->pixel_flag[ipx] >= 0)) {
            int x_mybin = -1;
            int y_mybin = -1;
            if (x_logscale) {
              if ((log10(img1->pixel_data[ipx]) > x_min) && (log10(img1->pixel_data[ipx]) < x_max))
                x_mybin = (int)floor((log10(img1->pixel_data[ipx]) - x_min)/x_delta);
              else {
                n_rejects++;
                n_x_oob++;
              }
            }//x_logscale
            else {//x_linscale
              if ((img1->pixel_data[ipx] > x_min) && (img1->pixel_data[ipx] < x_max))
                x_mybin = (int)floor((img1->pixel_data[ipx] - x_min)/x_delta);
              else {
                n_rejects++;
                n_x_oob++;
              }
            }//x_linscale
            if (y_logscale) {
              if ((log10(img2->pixel_data[ipx]) > y_min) && (log10(img2->pixel_data[ipx]) < y_max))
                y_mybin = (int)floor((log10(img2->pixel_data[ipx]) - y_min)/y_delta);
              else {
                n_rejects++;
                n_y_oob++;
              }
            }//y_logscale
            else {
              if ((img2->pixel_data[ipx] > y_min) && (img2->pixel_data[ipx] < y_max))
                y_mybin = (int)floor((img2->pixel_data[ipx] - y_min)/y_delta);
              else {
                n_rejects++;
                n_y_oob++;
              }
            }//y_linscale
            if ((x_mybin > -1) && (y_mybin > -1)) {
              bincount[x_mybin][y_mybin]++;
              n_good++;
            }//(valid x and y indices)
          }//(pixel_flags are valid)
          else n_bb++;
        }//for jp
      }//for ip
      cout << "Done binning: n_valid = " << n_good << ", n_bahama_blue = " << n_bb << ", n_out_of_bounds = " << n_rejects << ", n_total = "  << n_total << "\n";
      if (n_rejects > 0) cout << "n_x_oob = " << n_x_oob << ", n_y_oob = " << n_y_oob << "\n";

      //    Dump the histogram
      //    header
      ofile << "#\t" << x_min << "\t" << x_max << "\t" << x_nbins << "\t" << y_min << "\t" << y_max << "\t" << y_nbins << "\n";
      //    body
      for (int x_bin = 0; x_bin < x_nbins; ++x_bin) {
        for (int y_bin = 0; y_bin < y_nbins; ++y_bin) {
          ofile << x_bin << "\t" << y_bin << "\t" << bincount[x_bin][y_bin] << "\n" << flush;
        }//for y_bin
      }//for x_bin

      ofile.close();
      img1->lock = false;
      img2->lock = false;
    }//for iI
    DELETE(bincount);
  }//itp2DPdf()

  void itpDiscAverage(Param * param) {//Statistics along the major axis of a (multiple) cylinder(s)

    COUT1("DISC_AVERAGE...");

    if ((!pngCache)||(pngCache->imageList.size() == 0)) {
      CERR("cannot compute DISC_AVERAGE: pngCache is empty.");
    }
    pngCache->clearImageCache(); //If multiple RADIAL_PROBES are requested, must re-read the images
                                 //to clear flagging of marked pixels.

    //  Parse the param file to get probe quantity, locations, and dump name
    int arg = 0;
    string outfilename = param->getString(arg++);
    cout << "Dumping data to file " << outfilename << "\n";
    int const nprobes = param->getInt(arg++);
    int seed_i[nprobes];
    int seed_j[nprobes];
    double radius[nprobes];
    for (int probe = 0; probe < nprobes; probe++) {
      seed_i[probe] = param->getInt(arg++);
      seed_j[probe] = param->getInt(arg++);
      radius[probe] = param->getDouble(arg++);
      if (radius[probe]<0.0) CERR("DISC_AVERAGE radius must be positive: " << radius);
      cout << "Computing a disc average probe centered at pixel: " << seed_i[probe] << ", " << seed_j[probe] << " with radius " << radius[probe] << "\n";
    }//for probe

    //  Create the output file
    ofstream ofile;
    size_t fwidth = 13;
    ofile.open(outfilename.c_str());
    if (!ofile.is_open()) CERR("could not open file: " << outfilename);

    //  Dump the output header
    for (int probe = 0; probe < nprobes; probe++) ofile << "Probe " << probe << " , pix = " << seed_i[probe] << ", " << seed_j[probe] << ", radius = " << radius[probe] << "\n";
    ofile << "# " << setw(fwidth-2) << "1:index";
    for (int probe = 0; probe < nprobes; probe++) {
      //    the variable name...
      ofile << setw(fwidth) << "2:avg" << setw(fwidth) << "3:rms" << setw(fwidth) << "4:min" << setw(fwidth) << "5:max" ;
    }//for (probe)
    ofile << "\n";;

    //  Loop over the images and compute the statistics
    double stats[nprobes*4];
    int count[nprobes];
    for (int iI = 0; iI<pngCache->size(); ++iI){
      PngData * img = pngCache->getImage(iI);
      //    Initialize the buffers for computing statistics
      cout << " working on image \"" << pngCache->imageList[iI] << "\n";
      for (int probe = 0; probe < nprobes; probe++) {
        count[probe]   = 0;
        stats[4*probe+0] = 0;//avg
        stats[4*probe+1] = 0;//rms
        stats[4*probe+2] = HUGE_VAL;//min
        stats[4*probe+3] = -HUGE_VAL;//max
      }//probe

      //    Loop over the pixels to determine whether or not they're included in each probe region
      cout << "Nx = " << img->getNx() << ", Ny = " << img->getNy() << "\n";
      for (int pix_x = 0; pix_x < img->getNx(); pix_x++) {
        for (int pix_y = 0; pix_y < img->getNy(); pix_y++) {
          int pix = img->getNx()*pix_y + pix_x;
	  //        Exclude bahama blue
          if (img->pixel_flag[pix] >= 0) {
            for (int probe = 0; probe < nprobes; probe++) {
              double delta_x = (double)pix_x - (double)seed_i[probe];
              double delta_y = (double)pix_y - (double)seed_j[probe];

              double my_r = sqrt(delta_x*delta_x + delta_y*delta_y);
              if (my_r < (double) radius[probe]) {
		//              Valid point; augment the statistics
                count[probe]++;
                stats[4*probe + 0] += img->pixel_data[pix];//avg
                stats[4*probe + 1] += img->pixel_data[pix]*img->pixel_data[pix];//rms
                if (img->pixel_data[pix] < stats[4*probe + 2]) stats[4*probe + 2] = img->pixel_data[pix];//min
                if (img->pixel_data[pix] > stats[4*probe + 3]) stats[4*probe + 3] = img->pixel_data[pix];//max
              }//(my_r < radius[probe])
            }//probe
          }//(pixel_flag >= 0)
        }//pix_y
      }//pix_x

      ofile << setw(fwidth-2) << iI << " ";
      for (int probe = 0; probe < nprobes; probe++) {
        cout << "Probe " << probe << " resulted in " << count[probe] << " flagged pixels\n";
	//      Final steps to convert sum->avg; variance->rms
        stats[4*probe+0] /= count[probe];//avg
        stats[4*probe+1] = sqrt(stats[4*probe+1]/count[probe] - stats[probe+0]*stats[4*probe+0]);
	//      Dump the stats to the ouptut file
        ofile << setw(fwidth) << stats[4*probe+0] << setw(fwidth) << stats[4*probe+1] <<\
	  setw(fwidth) << stats[4*probe+2] << setw(fwidth) << stats[4*probe+3] << "\t";
      }//probe
      ofile << "\n";

      img->lock = false;
    }//iI

    //  Close the output file
    ofile.close();
  }//itpDiscAverage()

  void itpNotDiscAverage(Param * param) {//Statistics along the major axis of a (multiple) cylinder(s)

    COUT1("DISC_AVERAGE...");

    if ((!pngCache)||(pngCache->imageList.size() == 0)) {
      CERR("cannot compute DISC_AVERAGE: pngCache is empty.");
    }
    pngCache->clearImageCache(); //If multiple RADIAL_PROBES are requested, must re-read the images
                                 //to clear flagging of marked pixels.

    //  Parse the param file to get probe quantity, locations, and dump name
    int arg = 0;
    string outfilename = param->getString(arg++);
    cout << "Dumping data to file " << outfilename << "\n";
    int const nprobes = param->getInt(arg++);
    int seed_i[nprobes];
    int seed_j[nprobes];
    double radius[nprobes];
    for (int probe = 0; probe < nprobes; probe++) {
      seed_i[probe] = param->getInt(arg++);
      seed_j[probe] = param->getInt(arg++);
      radius[probe] = param->getDouble(arg++);
      if (radius[probe]<0.0) CERR("DISC_AVERAGE radius must be positive: " << radius);
      cout << "Computing a disc average probe centered at pixel: " << seed_i[probe] << ", " << seed_j[probe] << " with radius " << radius[probe] << "\n";
    }//for probe

    //  Create the output file
    ofstream ofile;
    size_t fwidth = 13;
    ofile.open(outfilename.c_str());
    if (!ofile.is_open()) CERR("could not open file: " << outfilename);

    //  Dump the output header
    for (int probe = 0; probe < nprobes; probe++) ofile << "Probe " << probe << " , pix = " << seed_i[probe] << ", " << seed_j[probe] << ", radius = " << radius[probe] << "\n";
    ofile << "# " << setw(fwidth-2) << "1:index";
    for (int probe = 0; probe < nprobes; probe++) {
      //    the variable name...
      ofile << setw(fwidth) << "2:avg" << setw(fwidth) << "3:rms" << setw(fwidth) << "4:min" << setw(fwidth) << "5:max" ;
    }//for (probe)
    ofile << "\n";;

    //  Loop over the images and compute the statistics
    double stats[nprobes*4];
    int count[nprobes];
    for (int iI = 0; iI<pngCache->size(); ++iI){
      PngData * img = pngCache->getImage(iI);
      //    Initialize the buffers for computing statistics
      cout << " working on image \"" << pngCache->imageList[iI] << "\n";
      for (int probe = 0; probe < nprobes; probe++) {
        count[probe]   = 0;
        stats[4*probe+0] = 0;//avg
        stats[4*probe+1] = 0;//rms
        stats[4*probe+2] = HUGE_VAL;//min
        stats[4*probe+3] = -HUGE_VAL;//max
      }//probe

      //    Loop over the pixels to determine whether or not they're included in each probe region
      cout << "Nx = " << img->getNx() << ", Ny = " << img->getNy() << "\n";
      for (int pix_x = 0; pix_x < img->getNx(); pix_x++) {
        for (int pix_y = 0; pix_y < img->getNy(); pix_y++) {
          int pix = img->getNx()*pix_y + pix_x;
	  //        Exclude bahama blue
          if (img->pixel_flag[pix] >= 0) {
            for (int probe = 0; probe < nprobes; probe++) {
              double delta_x = (double)pix_x - (double)seed_i[probe];
              double delta_y = (double)pix_y - (double)seed_j[probe];

              double my_r = sqrt(delta_x*delta_x + delta_y*delta_y);
              if (my_r > (double) radius[probe]) {
		//              Valid point; augment the statistics
                count[probe]++;
                stats[4*probe + 0] += img->pixel_data[pix];//avg
                stats[4*probe + 1] += img->pixel_data[pix]*img->pixel_data[pix];//rms
                if (img->pixel_data[pix] < stats[4*probe + 2]) stats[4*probe + 2] = img->pixel_data[pix];//min
                if (img->pixel_data[pix] > stats[4*probe + 3]) stats[4*probe + 3] = img->pixel_data[pix];//max
              }//(my_r < radius[probe])
            }//probe
          }//(pixel_flag >= 0)
        }//pix_y
      }//pix_x

      ofile << setw(fwidth-2) << iI << " ";
      for (int probe = 0; probe < nprobes; probe++) {
        cout << "Probe " << probe << " resulted in " << count[probe] << " flagged pixels\n";
	//      Final steps to convert sum->avg; variance->rms
        stats[4*probe+0] /= count[probe];//avg
        stats[4*probe+1] = sqrt(stats[4*probe+1]/count[probe] - stats[probe+0]*stats[4*probe+0]);
	//      Dump the stats to the ouptut file
        ofile << setw(fwidth) << stats[4*probe+0] << setw(fwidth) << stats[4*probe+1] <<\
	  setw(fwidth) << stats[4*probe+2] << setw(fwidth) << stats[4*probe+3] << "\t";
      }//probe
      ofile << "\n";

      img->lock = false;
    }//iI

    //  Close the output file
    ofile.close();
  }//itpDiscAverage()
  //~JOB


  void itpInfo(Param * param){
    COUT1("Processing param \"INFO\"...");

    if ((!pngCache)||(pngCache->imageList.size() == 0)) {
      CERR("cannot report image INFO: pngCache is empty.");
    }

    for (int iI = 0; iI<pngCache->size(); ++iI){
      cout << endl;
      cout << " > working on image \"" << pngCache->imageList[iI] << "\"" << endl;
      PngData * img = pngCache->getImage(iI);
      img->dumpInfo(cout);
    }


  }

  void itpEncode(Param * param){
    COUT1("Processing param \"ENCODE\"...");

    if ((!pngCache)||(pngCache->imageList.size() == 0)) {
      CERR("cannot ENCODE images: pngCache is empty.");
    }
#ifndef WITH_X264

    CERR("Compile with -DWITH_X264 to use ping_s.exe --ENCODE \n" <<
         "to produce animations from png files.  Requires libx264,\n" <<
         "https://www.videolan.org/developers/x264.html");

#else
    //setup options...
    int iarg = 0;
    int fps = 20;
    string prefix = getPrefix(pngCache->imageList[0]);

    while (iarg < param->size()) {
      const string token = param->getString(iarg++);
      if (token == "FPS") {
        fps = param->getDouble(iarg++);
      }
      else if (token == "NAME") {
        prefix = param->getString(iarg++);
      }
    }

    stringstream filename_ss;
    filename_ss << prefix << "_" << fps << "fps.264";
    COUT1("Encoding images to H264 file: " << filename_ss.str());
    MiscUtils::mkdir_for_file(filename_ss.str());

    FILE * fp;
    fp = fopen (filename_ss.str().c_str(), "wb");

    PngData * img0 = pngCache->getImage(0,false);
    const int img0_nx = img0->getNx();
    const int img0_ny = img0->getNy();


    x264_param_t param264;
    x264_param_default_preset(&param264, "veryfast", "zerolatency");
    param264.i_threads = 1;
    param264.i_width = img0->getNx();
    param264.i_height = img0->getNy();
    param264.i_fps_num = fps;
    param264.i_fps_den = 1;
    // Intra refres:
    param264.i_keyint_max = fps;
    param264.b_intra_refresh = 1;
    //Rate control:
    param264.rc.i_rc_method = X264_RC_CRF;
    param264.rc.f_rf_constant = 25;
    param264.rc.f_rf_constant_max = 35;
    //For streaming:
    param264.b_repeat_headers = 1;
    param264.b_annexb = 1;
    x264_param_apply_profile(&param264, "baseline");

    //initialize encoder
    x264_t* encoder = x264_encoder_open(&param264);
    x264_picture_t pic_in, pic_out;

    x264_picture_alloc(&pic_in, X264_CSP_I420, img0->getNx(), img0->getNy());

    for (int iI = 0; iI<pngCache->size(); ++iI){
      cout << "working on image " << pngCache->imageList[iI] << endl;

      //read image but skip finalize read, we
      //will work directly with the rgb data.
      PngData * img = pngCache->getImage(iI,false);

      if (img->getNx()!=img0_nx || img->getNy()!=img0_ny){
        CERR("  This image is a difference size " << img->getNx() << "x" << img->getNy() << "than the previous image " << img0_nx << "x" << img0_ny);
      }

      img->rgbToYuv420p(pic_in.img.plane[0]);

      x264_nal_t* nals;
      int i_nals;
      int i_frame_size = x264_encoder_encode(encoder, &nals, &i_nals, &pic_in, &pic_out);
      if( i_frame_size < 0 ){
        CERR("Problem encoding movie");
      }
      else if( i_frame_size ) {
        if( !fwrite( nals->p_payload, i_frame_size, 1, fp ) )
          CERR("Problem writing movie");
      }

      img->lock = false;
    }
    fclose(fp);
    x264_encoder_close(encoder);
    x264_picture_clean(&pic_in);

#endif
  }

  //Blank pixels that do not share the same pixel type across
  //a global set of images. Images must be the same size.
  //Initial application is to prep image sets with moving geometry
  //for imageDmd/Pod
  void itpGlobalMask(Param * param){
 
    if ((!pngCache)||(pngCache->imageList.size() == 0)) {
      CERR("cannot perform GLOBAL_MASK: pngCache is empty.");
    }
    else{
      cout << " > GLOBAL_MASK" << endl;
      cout << " >   building mask..." << endl;
    }

    string maskPath = "images_mask";

    //read first image, setup mask
    PngData * img = pngCache->getImage(0);
    const int nx = img->getNx();
    const int ny = img->getNy();
    int * pixel_flag_mask = new int[nx*ny];
    for (int ipx =0; ipx < nx*ny ; ++ipx) {
      pixel_flag_mask[ipx] = img->pixel_flag[ipx];
    }
    img->lock = false;
    //set non-matching pixels to background
    for (int iI = 1; iI<pngCache->size(); ++iI) {
      img = pngCache->getImage(iI);
      if (nx!=img->getNx()||ny!=img->getNy()){
        CERR("images must be the same size. " << pngCache->imageList[0] << " is " << nx << "x" << ny << ", " << pngCache->imageList[iI] << " is " << img->getNx() << "x" << img->getNy());
      }
      for (int ipx =0; ipx < nx*ny ; ++ipx) {
        if(pixel_flag_mask[ipx] != img->pixel_flag[ipx]){
          pixel_flag_mask[ipx] = PngData::BACKGROUND_PIXEL;
        }
      }
      img->lock = false;
    }
    //the mask was built successfully, now process (blank) the images
    COUT1(" >   writing files to " << maskPath);
    MiscUtils::mkdir_for_file(maskPath+"/test.png");

    cout << " >   working on image 0...";
    cout.flush();
    for (int iI = 0; iI<pngCache->size(); ++iI) {
      if (iI%10==0&&iI>0){
        cout << iI << "...";
        cout.flush();
      }
      img = pngCache->getImage(iI);
      img->blankWithMask(pixel_flag_mask);

      //build output path and filename
      //use original filename in a different folder
      vector<string> inputNameTokens;
      tokenizeString(inputNameTokens,pngCache->imageList[iI],"/");
      string outputPathAndFilename = maskPath + "/";
      outputPathAndFilename += inputNameTokens[inputNameTokens.size()-1];

      img->initializeWrite();
      img->write(outputPathAndFilename.c_str());
      img->lock = false;
    }
    cout << "done" << endl;

  }

  string getPathAndPrefix(const string& filename) {

    // here we assume that the filename has the form
    // <path>/<prefix>.<index>.<ext>

    size_t pos_path = filename.find_last_of('/');
    if (pos_path == string::npos)
      pos_path = 0;

    size_t pos00 = pos_path;
    size_t pos0 = pos_path;
    size_t pos = pos_path;
    while (pos != string::npos) {
      pos00 = pos0;
      pos0 = pos;
      pos = filename.find_first_of('.',pos0+1);
    }

    if (pos00 == pos_path) {
      // filename is there appears to be only one "." in the image, so just
      // return the path as the prefix. It will put the avg, rms, etc
      // into the specified directory.
      if (pos00 == 0) {
        // filename must be "<prefix>.<ext>" and no path. Just return an empty
        // string and the output images will be simply "avg.png", etc. in the
        // current directory...
        return "";
      }
      else {
        // filename must be "<path>/<prefix>.<ext>". Just return the path up to
        // and INCLUDING the "/". This will put any output images into
        // "<path>/avg.png" etc...
        return filename.substr(0,pos00+1);
      }
    }
    else {
      // file is as expected with format "<path>/<prefix>.<index>.<ext>" or
      // perhaps "<prefix>.<index>.<ext>". return "<path>/<prefix>."...
      return filename.substr(0,pos00+1);
    }
  }

  string getPrefix(const string& filename) {

    // here we assume that the filename has the form
    // <path>/<prefix>.<index>.<ext>

    size_t pos_path = filename.find_last_of('/');
    if (pos_path == string::npos) pos_path = 0;
    else ++pos_path;  // start at char after last slash

    size_t prefix_end = filename.find_first_of('.',pos_path);
    return filename.substr(pos_path,(prefix_end-pos_path));
  }


  void tokenizeString(vector<string> &tokens, const string &str, const string &delimiters) {
    string::size_type lastPos = str.find_first_not_of(delimiters, 0); // skip delimiters at beginning
    string::size_type pos = str.find_first_of(delimiters, lastPos); // find first "non-delimiter"

    while (string::npos != pos || string::npos != lastPos) {
      tokens.push_back(str.substr(lastPos, pos - lastPos)); // add token to the vector<string>
      lastPos = str.find_first_not_of(delimiters, pos); // skip delimiters
      pos = str.find_first_of(delimiters, lastPos); // find next "non-delimiter"
    }
  }


  //An alternative to writing to the original image location...
  //strips the leading filepath and trailing .png from an image filename
  //the output can be used to write ascii files to the ping_s run
  //directory or another user specified location
  string buildFilenameFromImage(const string &label, const string &imageFilename){
    vector<string> inputNameTokens;
    tokenizeString(inputNameTokens,imageFilename,"/");
    string temp = inputNameTokens[inputNameTokens.size()-1];
    inputNameTokens.clear();
    tokenizeString(inputNameTokens,temp,".");
    string outputPrefix = "";
    for (int ii=0;ii<inputNameTokens.size()-1;++ii){ //execlude .png
      outputPrefix += inputNameTokens[ii];
      outputPrefix += ".";
    }
    outputPrefix += label;
    return outputPrefix;
  }

  string getPixelTypeLegendStr(){
    stringstream ss;
    ss << PngData::BACKGROUND_PIXEL << ":background ";
    ss << PngData::SURFACE_PIXEL << ":surface ";
    ss << PngData::VOLUME_DATA_PIXEL << ":volume_data ";
    ss << PngData::SURFACE_DATA_PIXEL << ":surface_data ";
    ss << PngData::PARTICLE_DATA_PIXEL << ":particle_data ";
    ss << PngData::ISO_DATA_PIXEL << ":iso_data";
    return ss.str();
  }


  /*
    void setPathAndPrefixAndSuffix(string& pathAndPrefix,string& suffix,const string& filename) {

    // here we assume that the filename has the form
    // <path>/<prefix>.<index>.<ext>

    size_t pos_path = filename.find_last_of('/');
    if (pos_path == string::npos)
    pos_path = 0;

    size_t pos00 = pos_path;
    size_t pos0 = pos_path;
    size_t pos = pos_path;
    while (pos != string::npos) {
    pos00 = pos0;
    pos0 = pos;
    pos = filename.find_first_of('.',pos0+1);
    }

    if (pos00 == pos_path) {
    // there appears to be only one "." in the image, so just
    // return the path as the prefix. It will put the avg, rms, etc
    // into the specified directory.
    if (pos00 == 0) {
    // filename must be "<prefix>.<ext>" and no path. Just return an empty
    // string and the output images will be simply "avg.png", etc. in the
    // current directory...
    pathAndPrefix = "";
    suffix = filename.substr(pos0
    return "";
    }
    else {
    // filename must be "<path>/<prefix>.<ext>". Just return the path up to
    // and INCLUDING the "/". This will put any output images into
    // "<path>/avg.png" etc...
    return filename.substr(0,pos00+1);
    }
    }
    else {
    // file is as expected with format "<path>/<prefix>.<index>.<ext>" or
    // perhaps "<prefix>.<index>.<ext>". return "<path>/<prefix>."...
    return filename.substr(0,pos00+1);
    }

    }
  */

  void itpHelp(Param * param) {

    if (mpi_rank == 0) {

      cout <<
        "\n==========================================================================================================\n" <<
        "ping HELP\n\n" <<
        " IMAGES <file1> [<file2>...]               # Specify image filenames directly in the parameter line\n" <<
        "                                           # This option can parse wildcards when used at the command line.\n\n" <<
        " IMAGELIST <filename>                      # Specify a file containing the list of image files\n\n" <<
	" ANNULAR_INTEGRATION_PROBE NBIN 18 RMIN 0.17 RMAX 0.245 [CLOCKING] [CHECK]\n" <<
	"                                           # annular integration probe capability. Adding the CHECK option\n" <<
	"                                           # writes out an image that illustrates where the probe integrations\n" <<
	"                                           # are located. CLOCKING toggles the clocking.\n\n" << 
        " ADD_COLORMAP FILE <colormap.dat>          # Allows for the creation of a user-defined colormap. Once added\n" <<
        "                                           # it may be used for recoloring as any available colormap would be.\n" <<
        "                                           # <colormap.dat> is an ASCII text file with the following format:\n" <<
        "                                           # row 1           - colormap name           e.g. MYCOLORMAP\n" <<
        "                                           # row 2           - color space 1 or 255    e.g. 255 \n" <<
        "                                           # row 3           - # of control points Ncp e.g. 5\n" <<
        "                                           # rows 4 to Ncp+4 - <phi> <r> <g> <b>       e.g. 0.2 73 10 25\n" <<
        "                                           # The rows containing control point colors should be arranged\n" <<
        "                                           # in order of increasing <phi> where 0 <= <phi> <= 1. RGB values\n" <<
        "                                           # should be formatted according to the colorspace indicated.  They\n" <<
        "                                           # should vary from a minimum of 0 to a maximum of 1 or 255.  The \n" <<
        "                                           # colormap name specified may be used to recolor any datatype in the\n" <<
        "                                           # RECOLOR command. \n\n" <<
        " SPATIAL_STATS [MASK <maskfile>] [<outputfile>] \n" <<
        "                                           # Computes spatial averages of individual images and writes\n" <<
        "                                           # ascii table to the specified output file. A separate output file is\n" <<
        "                                           # produced for each data type in the image (e.g. volume, surface). The\n" <<
        "                                           # default filenames are \"spatial_stats.<data_type>.dat\". If the MASK option\n" << 
        "                                           # is specified, <maskfile> is a list of pixel (i,j) coordinates that define\n" <<
        "                                           # a polygon (not necessarily convex) that contains the region of interest. \n\n" <<
        " RADIAL_PROFILE [WGT] [NBIN <n>] [NAME <outputprefix>]\n" <<
        "                                           # Computes azimuthally-averaged radial profile (e.g. for combustor outlet temperature).\n" <<
        "                                           # Automatically determines the image plane and computes 2d radius appropriately.\n" << 
        "                                           # If WGT is specified, it is assumed that 2*N images have been specified, where the\n" << 
	"                                           # first N are the variable to be averaged, and the second N are the weight (for example,\n" << 
	"                                           # outlet Temperature can be weighted by outlet-normal momentum). Volume data only.\n\n" <<
        " SPATIAL_STATS_DISK X <x> <y> <z> R <r> [<outputfile>] \n" <<
        "                                           # Computes spatial averages of individual images and writes\n" <<
        "                                           # ascii table to the specified output file. Volume data only.\n\n" <<
        " TUMBLE MASK <maskfile> [RPM <rpm>] [RPT <rpt>] [OMEGA <omega>] [RANGE <min> <max>] [<outputprefix>] \n" <<
        "                                           # Computes the normalized tumble using planar image pairs that contain images of the\n" << 
        "                                           # i-aligned and j-aligned components of velocity. It is assumed that IMAGES or IMAGELIST\n" << 
        "                                           # has all i-images first, followed by all j-images. Images must be matched.\n" <<
        "                                           # Only volume (i.e. plane-cut) data is supported. The default prefix is \"tumble\".\n" <<
        "                                           # MASK is required to use TUMBLE, where <maskfile> is a list of pixel (i,j) coordinates that define\n" <<
        "                                           # a polygon (not necessarily convex) that contains the region of interest. The\n" <<
        "                                           # optional parameter RANGE sets the range of the output tumble images.\n\n" << 
	" WRITE_VELOCITY [MASK <maskfile>] [STARTING_INDEX <index>] [STRIDE <int>] [UW_ONLY] [NAME <outputprefix>] \n" <<
        "                                           # Reads planar image trios that contain images of the x, y, and z\n" << 
        "                                           # components of velocity, and writes data files containing the physical\n" << 
	"                                           # coordinates of the pixels and all three velocity components. It is assumed that IMAGES\n" << 
        "                                           # or IMAGELIST has all x-images first, followed by all y-images, and then all z-images,\n" << 
        "                                           # unless the UW_ONLY flag is specified, then all x and all z only.\n" << 
	"                                           # Images must be matched. Only volume (i.e. plane-cut) data is supported.\n" <<
        "                                           # The default prefix is \"velocity\", and the default file names are \"uvw.xxxxxx.dat\".\n" <<
	"                                           # The files start with the starting index if provided, otherwise, start with 0. \n" <<
        "                                           # MASK is optional, where <maskfile> is a list of pixel (i,j) coordinates that define\n" <<
        "                                           # a polygon (not necessarily convex) that contains the region of interest. STRIDE defaults\n" <<
	"                                           # to 1 meaning every pixel. STRIDE > 1 will average the pixels in the STRIDExSTRIDE group.\n\n" << 
	" WRITE_SCALAR [MASK <maskfile>] [STARTING_INDEX <index>] [STRIDE <int>] [NAME <outputprefix>] \n" <<
	"                                           # writes image data to files with columns x,y,z,value. Similar to\n" <<
	"                                           # options in WRITE_VELOCITY above.\n\n" << 
        " TEMPORAL_STATS [<output_path>] \n" <<
        "                                           # Computes temporal stats of images and writes output to\n" <<
        "                                           # new image files prefix.avg.png, prefix.rms.png, prefix.min.png,\n" <<
        "                                           # and prefix.max.png. The stats image filenames are based on input file\n" <<
        "                                           # prefixes.  The default output path if not specified is images_stats.\n" <<
        "                                           # NOTE: all images must contain the same pixel layout, i.e. pixel (i,j)\n" <<
        "                                           # must be of the same type (surface vs volume data, etc) for all images.\n" <<
        "                                           # Otherwise, consider TEMPORAL_STATS_MOVING or TEMPORAL_STATS_PHASE.\n\n" <<
        " TEMPORAL_STATS_CROSS [<output_path>] \n" <<
        "                                           # Similar to TEMPORAL_STATS but computes a cross correlation between images.\n" <<
        "                                           # This commands expects two groups of images passed to the IMAGES parameter,\n" <<
        "                                           # i.e. IMAGES <image_set0.*.png> <image_set1.*.png>. The images in set0 and set1\n" <<
        "                                           # should be ordered identically in time and be of the same length. \n" <<
        "                                           # NOTE: all images must contain the same pixel layout, i.e. pixel (i,j)\n" <<
        "                                           # must be of the same type (surface vs volume data, etc) for all images.\n" <<
        "                                           # A avg and rms image is reported for each image set and a cross correlation\n" <<
        "                                           # image is also written.\n\n" <<
        " TEMPORAL_STATS_MOVING [<output_path>] \n" <<
        "                                           # Similar to TEMPORAL_STATS but allows the pixel layout to vary between\n" <<
        "                                           # images. Will be identical to TEMPORAL_STATS_PHASE with a single bin\n" <<
        "                                           # and a uniform image weighting (time interval). Avg and rms images are\n" <<
        "                                           # written for each data type found in the image set (volume, surface,\n"  <<
        "                                           # iso, particle).\n\n" <<
        " TEMPORAL_STATS_PHASE PERIOD <time> NBIN <count> [IMGDT <time interval>] [<output path>]\n" <<
        "                                           # Computes phase-averaged temporal stats of images and writes output to\n" <<
        "                                           # new image files prefix.avg-000.png, prefix.avg-001.png, etc.\n" <<
        "                                           # The prefix is computed from input filenames and the default output\n" <<
        "                                           # path is images_stats. The time interval between images will be computed\n" <<
        "                                           # automatically from the image TIME metadata. If not present the user should\n" <<
        "                                           # specify an IMGDT.\n\n" <<
        " TEMPORAL_STATS_PHASE_LOCK PERIOD <time> NBIN <count> [BINDEX <bin index>] [IMGDT <time interval>] [LOCK <lock file>] [PHASE_TYPE] [<output path>]\n" <<
        "                                           # Computes phase-averaged temporal stats of images and writes output to\n" <<
        "                                           # new image files prefix.avg-000.png, prefix.avg-001.png, etc.\n" <<
        "                                           # The prefix is computed from input filenames and the default output\n" <<
        "                                           # path is images_stats. The time interval between images will be computed\n" <<
        "                                           # automatically from the image TIME metadata. If not present the user should\n" <<
        "                                           # specify an IMGDT.  User can set a lock file to specify the nodes in time\n" <<
        "                                           # where the signal cycle changes.  This overrides the PERIOD.  If PERIOD is\n" <<
        "                                           # not positive (or not specified), PERIOD is computed from the lock file. The\n" <<
        "                                           # lock file must have the number of nodes followed by the nodal time values.\n" <<
        "                                           # PHASE_TYPE specifies the type of bin-averaging for second-order statistics.\n" <<
        "                                           # Current options are BIN (std) or LIN (for linear interpolation between bins).\n" <<
        "                                           # User can process a specific bin [0..NBIN] by setting a BINDEX value. If no\n" <<
        "                                           # BINDEX value given, all bins are computed.  Compute bins in parallel with BINDEX.\n" <<
        "                                           # User can set ranges with RANGE_[/SURFACE/PARTICLES/ISO]_[AVG/RMS] like RESCALE.\n\n" <<
        " DIFF [SIGNED] [<path>] \n" <<
        "                                           # Computes the diff between the first image in the list and \n" <<
        "                                           # subsequent images.  By default the absolute value of the diff \n" <<
        "                                           # is reported. Include the keyword SIGNED to preserve the sign\n" <<
        "                                           # of the diff relative to the first image.\n" <<
        "                                           # Red regions indicate missing geometry and green regions\n" <<
        "                                           # indicate additional geometry as compared to the baseline image,\n" <<
        "                                           # respectively. Blue regions indicate a data type mismatch. \n" <<
        "                                           # Output images are named \"path/diff.<image_prefix>.png\" where\n" <<
        "                                           # image_prefix is from the baseline image if not specified. The\n" <<
        "                                           # output image data range is set to the global min/max diff from\n" <<
        "                                           # all compared images.\n\n" <<
        " RESCALE <min> <max> [PATH <path>] \n" <<
        " RESCALE GLOBAL [TIGHT or LOOSE] [PATH <path>] \n" <<
        " RESCALE [VOLUME <min> <max>] [SURFACE <min> <max>] [PARTICLES <min> <max>] [ISO <min> <max>] [PATH <path>] \n" <<
        "                                           # Rescale the variable(s) range(s) represented by image pixel color.\n" <<
        "                                           # Rescale supports multiple data types with multiple variables in a .\n" <<
        "                                           # single image. These may be rescaled explicity to an exact range \n" <<
        "                                           # by using the data type keyword followed by <min> and <max>. If all \n" <<
        "                                           # data types show the same variable, omitting the data type\n" <<
        "                                           # keyword will rescale all data types. The GLOBAL keyword computes \n" <<
        "                                           # a range from image metadata for each data type.  TIGHT uses the \n" <<
        "                                           # global minimum range while LOOSE uses the global maximum range. A \n" <<
        "                                           # custom output folder may be specified with PATH <path>.  The \n" <<
        "                                           # default output directory is <run_directoy>/images_rescaled\n\n" <<
        " RECOLOR [COLORMAP <colormap name>] [COLORMAP_SURFACE <colormap name>] \n"
        "         [COLORMAP_PARTICLES <colormap name>] [COLORMAP_ISO <colormap name>] [PATH <path>] \n" <<
        "         [LEGEND [IJ <ifrac> <jfrac>]] [LEGEND_SURFACE [IJ <ifrac> <jfrac>]] \n" <<
        "         [LEGEND_PARTICLES [IJ <ifrac> <jfrac>]] [LEGEND_ISO [IJ <ifrac> <jfrac>]] \n" <<
        "         [BACKGROUND WHITE|BLACK] [CHECKTIME] \n" <<
        "                                           # Recolor an existing GRAYSCALE image. Other image colormaps cannot\n" <<
        "                                           # be recolored.  Valid colormap names are HOT_METAL, BLUE_RED_RAINBOW,\n" <<
        "                                           # COOL_WARM, and RAINBOW_DESATURATED. To invert a given colormap scheme\n" <<
        "                                           # preface any of the previous names with \"INVERTED_\". Volume (planar),\n" <<
        "                                           # surface, particles, and iso data is recolored using the separate\n" <<
        "                                           # COLORMAP keywords above. The default output directory is\n" <<
        "                                           # <run_directoy>/images_recolored.  A custom output folder may be\n" <<
        "                                           # specified with PATH <path>. A legend may optionally be\n" <<
        "                                           # included for any data type with the LEGEND keyword. To change the\n" <<
        "                                           # default legend location specify a <ifrac> and <jfrac> from 0 to 1.\n" <<
        "                                           # The location is measured from the upper left corner.\n\n" <<
        " CROP <i0> <j0> <i1> <j1> [PATH <path>]\n" <<
        "                                           # Crop images in the current image cache. <i0> <j0> and <i1> <j1> are in image\n" <<
        "                                           # pixel coordinates, i.e. (0,0) is at the top left of the image. The default path is \"crop\"\n" <<
        " ENCODE [FPS <frames-per-second>] [NAME <filename-prefix>]\n" <<
        "                                           # Encode all images in the current list to a movie using the x264 codec.\n" <<
        "                                           # By default the filename will be constructed from the prefix of the\n" <<
        "                                           # first image in the list and the framerate. The framerate in frames per \n" <<
        "                                           # second can be set with FPS, if omitted the default is 20fps. The default\n" <<
        "                                           # output file prefix can be overridden with NAME.\n" <<
        " PROJECTED_SURFACE_AREA [WRITE] \n" <<
        "                                           # Compute the boundary surface area projected in the image normal \n" <<
        "                                           # direction.  All boundary surfaces in the images are included, surfaces\n" <<
        "                                           # with data are excluded. To visualize the flagged region include the\n" <<
        "                                           # WRITE keyword to generate a sample image with the region in green.\n\n" <<
        " PROBE GEOM POINT XP <x> <y> <z> [<outputfile>] \n" <<
        "                                           # Find a point in the image closest to the requested point and\n" <<
        "                                           # return its value. The simulation coordinates and image pixel \n" <<
        "                                           # of the located point is also output for reference.\n\n" <<
        " PROBE GEOM POINT IJ <i> <j> [<outputfile>] \n" <<
        "                                           # Report the value at the requested image i,j point (standard image \n" <<
        "                                           # indexing with top-left = 0,0)\n\n" <<
        " PROBE GEOM LINE XP <x0> <y0> <z0> <x1> <y1> <z1> \n" <<
        "       NPOINTS <npoints> [<outputfile>]           \n" <<
        "                                           # Output the value at <npoints> along a line running from point\n" <<
        "                                           # <x0> <y0> <z0> to pixel <x1> <y1> <z1>. Output is directed to \n" <<
        "                                           # either <outputfile> or, if omitted, ''probe_output.dat''.\n\n" <<
        " PROBE GEOM LINE IJ <i0> <j0> <i1> <j1>          \n" <<
        "       NPOINTS <npoints> [<outputfile>]          \n" <<
        "                                           # Output the value at <npoints> along a line running from pixel\n" <<
        "                                           # <i0> <j0> to pixel <i1> <j1>. Output is directed to either\n" <<
        "                                           # <outputfile> or, if omitted, ''probe_output.dat''.\n\n" <<
        " INTEGRATION_PROBE PIXEL <i> <j> [<outputfile>] \n" <<
        "                                           # Compute spatial stats on an image sequence with output similar\n" <<
        "                                           # to that from the SPATIAL_STATS option.  Stats are computed \n" <<
        "                                           # only on the pixels contiguous to input pixel (and of the same \n" <<
        "                                           # data type).  Output images are written with included pixels in \n" <<
        "                                           # green.\n\n" <<
        " RADIAL_PROBE PIXEL <i> <j> <r> \n" <<
        "                                           # Compute a radial profile centered at i,j with radius r. \n" <<
        "                                           # The seed point i,j must fall on an image pixel with data and \n" <<
        "                                           # only the data type of the seed point will be considered. \n" <<
        "                                           # Azimuthal stats are computed at each radial step.\n" <<
        "                                           # Output images are marked with the center in green and \n" <<
        "                                           # profile area in red.\n\n" <<
        " COLLAPSE_HEIGHT [PIXEL_TYPE <p>] [BETWEEN <j_f> <j_l>]\n" <<
        "                                           # For each image in IMAGES, calculate the average, min, max in\n" <<
        "                                           # each column of pixels. The range of pixels included in each column\n" <<
        "                                           # may be limited with the BETWEEN keyword followed by the index of the\n" <<
        "                                           # first and last row to include.  By default pixels of all data types are\n" <<
        "                                           # included.  Use the PIXEL_TYPE keyword to limit output to a single\n" <<
        "                                           # type where p is one integer from "<<PngData::VOLUME_DATA_PIXEL<<":Volume, "<<PngData::SURFACE_DATA_PIXEL<<":Surface, "<<PngData::PARTICLE_DATA_PIXEL<<":Particle,\n" <<
        "                                           # "<<PngData::ISO_DATA_PIXEL<<":Iso.\n\n" <<
        " COLLAPSE_WIDTH [PIXEL_TYPE <p>] [BETWEEN <i_f> <i_l>] \n" <<
        "                                           # Identical behavior to COLLAPSE_HEIGHT except the stats are performed\n" <<
        "                                           # across rows or portions of rows of pixels.\n\n" <<
        " DUMP_DATA \n" <<
        "                                           # Writes the pixel values to an ASCII file as a brute force means\n" <<
        "                                           # of accessing the image data.\n\n" <<
        " DUMP_COLUMN <i> \n" <<
        "                                           # Writes a single column of pixels from each image in IMAGES to an \n" <<
        "                                           # ASCII file. Specify the column number in pixels.\n\n" <<
        " 1D_PDF <xmin> <xmax> <nbins> [LOGSCALE]\n" <<
        "                                           # Creates a 1-dimensional pdf from the pixel values.\n\n" <<
        " 2D_PDF <xlist.txt> <xmin> <xmax> <n_xbins> \n" <<
        "        <ylist.txt> <ymin> <ymax> <n_ybins> \n" <<
        "        <output_prefix> \n" <<
        "                                           # Creates a 2-dimensional pdf for each pair of images contained\n" <<
        "                                           # <xlist.txt> and <ylist.txt>.\n\n" <<
        " DIFF_SET <list1.txt> <list2.txt> <output> <min> <max> \n" <<
        "                                           # Creates a set of images showing the difference in values\n" <<
        "                                           # <list1.txt> and <list2.txt> with colorbar scaled to <min>\n" <<
        "                                           # and <max>.\n\n" <<
        " ISOCONTOUR THRESHOLD <lo> <hi>            # Color pixels falling within a given range on multiple images\n\n" <<
        " ISOCONTOURS                               # Document me!!\n\n" <<
        " SUPERPOSE                                 # Document me!!\n\n" <<
        " GLOBAL_MASK \n" <<
        "                                           # Blank pixel locations in each image that do not share a common\n" <<
        "                                           # pixel type across all images.  Writes blanked images to a new folder\n" <<
        "                                           # named images_mask.\n\n" <<
        " INFO \n" <<
        "                                           # For each image in IMAGES, display a summary of image metadata.\n" <<
        "                                           # The summary reports values of standard png chunks teXt and ztXt\n" <<
        "                                           # as well as the presence of custom Cascade data chunks\n\n" <<
        " TILE <ncol> <nrow> CROP <i0> <j0> <i1> <j1> [PATH <output-path>] [CHECKTIME] [INDEX] \n" <<
        "                                           # build tiled images from input images. Used to synchronize images\n" <<
        "                                           # for rendering into movies. Any images metadata is removed from the\n" << 
        "                                           # final tiled images.\n" << 
        "==========================================================================================================\n";
    }//(root)
  }//itpHelp()
};//class ImageTools()


#endif
