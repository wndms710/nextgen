
#include "PngData.hpp"
#include <fstream>
#include <ctime> 
#include <vector> 



int main(const int argc, const char* argv[]) { 

  // get a sequence of images... 
  vector<string> imageList; 
  ifstream fin(argv[1]); 
  while(true) { 

    string tmp; fin >> tmp; 
    if ( fin.eof()) 
      break; 

    imageList.push_back(tmp);
  }
  fin.close();

  const int i_probe = atoi(argv[2]); 
  const int j_probe = atoi(argv[3]); 

  // check to make sure that we have located 
  // the right probe point... 
  int ipx_probe = -1; 
  { 
    PngData * png = new PngData(); 
    png->read(imageList[0].c_str(),true); 
    png->finalizeRead(); 
    png->buffer = new unsigned char[png->getNx()*png->getNy()*png->getStride()]; 

    vector<int> colorMe; 
    for (int i_width=-3; i_width < 3; ++i_width) { 
      for (int j_width = -3; j_width < 3; ++j_width) { 
        const int ii = i_probe + i_width; 
        const int jj = (png->getNy()-j_probe) + j_width; 
        const int kk = (png->getNy()-jj-1)*png->getNx() + ii;

        if ( (i_width ==0) && (j_width ==0) ) 
          ipx_probe = kk;

        colorMe.push_back(kk); 
      }
    }

    for (int ipx = 0; ipx < png->getNx()*png->getNy(); ++ipx) { 
      if ( std::find(colorMe.begin(), colorMe.end(), ipx) != colorMe.end()) { 
        png->buffer[png->getStride()*ipx+0] = 0; 
        png->buffer[png->getStride()*ipx+1] = 0; 
        png->buffer[png->getStride()*ipx+2] = 0;
      } else { 
        png->buffer[png->getStride()*ipx+0] = 255;
        png->buffer[png->getStride()*ipx+1] = 255;
        png->buffer[png->getStride()*ipx+2] = 255;
      }
    }//ipx

    png->setDepthFromPixelFlag();
    png->write("checkme.png"); 
    delete png; 
  }

  // now probe...
  const double dt_samp = atof(argv[4]); 
  for (int i=0; i < imageList.size(); ++i ) { 
    PngData * png = new PngData(); 
    png->read(imageList[i].c_str(),true); 
    png->finalizeRead(); 

    assert( png->pixel_flag[ipx_probe] >= 0 ); 
    const double time = dt_samp * double(i); 
    cout << " XXXX " << time << "    " << png->getPhi(ipx_probe) << endl; 

    delete png; 
  }


  return 0; 
} 
