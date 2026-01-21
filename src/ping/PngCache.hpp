#ifndef __PNG_CACHE__
#define __PNG_CACHE__

#include "PngData.hpp"
#include <fstream>
#include <vector>
#include <deque>
#include <sstream>

struct PngDataComp {
  int n;
  explicit PngDataComp(int i) : n(i) { }
  inline bool operator()(const PngData* png) const { return png->id == n; }
};//struct PngDataComp()

class PngCache {

public:
  string infile ;
 
  //int nImg;
//images..
  int ncache;
  deque<PngData*> imageCache;
  vector<string> imageList;
  
  PngCache(){}

  PngCache(const string& _infile) : infile(_infile) {
//  everybody needs to keep a list of images. 
    ifstream fin(infile.c_str());
    if (fin.fail()) {
      cerr << "problem opening IMAGE_LIST " << infile << ".\nEither the file doesn't exist, you don't have permission to open it, or something else went wrong. Quitting here...\n";
      throw(-1);
    }
    while ( true ) {
      string tmp; 
	  fin >> tmp;
      if ( fin.eof() )
        break;
      imageList.push_back(tmp);
    }
    fin.close();

    init();
  }

  void init() { 
    // number of images that each rank can hold in memory
    ncache = 15; // must be >= 2 

    /*
    nImg = imageList.size(); // number of snapshots..

    if ( mpi_rank == 0 ) {  
      cout << " > nImg = " << nImg  << endl;
    }
    */
    
    imageCache.clear(); 

  }

  int size() const {
    return imageList.size();
  }
  
  string getImageName(const int i) const {
    return imageList[i];
  }

  //any time we are working with the image as data, the finalize read should be performed...
  PngData* getImage(const int i, const bool b_doFinalizeRead=true) { 

    deque<PngData*>::iterator it = find_if(imageCache.begin(), imageCache.end(), PngDataComp(i)); 
 
    if ( it != imageCache.end() ) { 
      COUT2(" > retrieving image " << imageList[i] << " from cache");
      (*it)->lock = true;
      return *it; 
    } 
    else {
      COUT2(" > retrieving image " << imageList[i] << " from disk");
      if ( imageCache.size() == ncache) { 
        // cache is full, remove the first unlocked element (first in-first out)
        for(it = imageCache.begin(); it != imageCache.end(); ++it) { 
          if (!((*it)->lock)) 
            break;
        }//it
        assert( it != imageCache.end()); // uh oh, the whole cache is locked...
        delete *it; 
        imageCache.erase(it); 
      }
     
      // when we get the image, we will lock it, but 
      // the user is responsible for releasing the lock
      PngData* png = new PngData(); 
      png->read(imageList[i].c_str());
      if (b_doFinalizeRead) 
        png->finalizeRead();
      png->id = i;
      png->lock = true; 
      imageCache.push_back(png);

      assert( imageCache.size() <= ncache); 

      return png;
    }
        
    // shouldn't get here..
    return NULL; 
  } 

  ImageMetadata* getMetadata(const int i) { 

    deque<PngData*>::iterator it = find_if(imageCache.begin(), imageCache.end(), PngDataComp(i)); 
 
    if ( it != imageCache.end() ) { 
      COUT2(" > retrieving metadata from image " << imageList[i] << " from cache");
      return (*it)->getMetadata();
    } else {
      COUT2(" > retrieving metadata from image " << imageList[i] << " from disk");
      //if the image is not in the cache, we will only read metadata from disk
      //therefore this image will not be added to the cache
     
      PngImage* png = new PngImage(); 
      png->read(imageList[i].c_str(),false); //read only image metadata
      return png->getMetadata();
    }
        
    // shouldn't get here..
    return NULL; 
  } 

  void clearImageCache(){
    for (int i=0; i < imageCache.size(); ++i) {
      if(imageCache[i]!=NULL)	{
        delete imageCache[i];
        imageCache[i] = NULL;
      }
    }
    imageCache.clear();
  }

  ~PngCache() {
      
    imageList.clear();
      
    for (int i=0; i < imageCache.size(); ++i) {
      if(imageCache[i]!=NULL)	{
        delete imageCache[i];
        imageCache[i] = NULL;
      }
    }
    
    imageCache.clear();
    
  }
};//class pngCache()

#endif // defined(__PNG_ANALYSIS__)
