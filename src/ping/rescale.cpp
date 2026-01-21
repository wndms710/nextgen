
#include <iostream> 
#include <fstream>
#include "PngImage.hpp" 

class PngRescale : public PngImage { 

public: 
  void rescale_gray(const double phi_lo, const double phi_hi) { 
    const int n = nx*ny;
    for (int i =0; i < n; ++i) {
      assert(stride == 3);
      if ( (buffer[i*stride+0] != 73) || (buffer[i*stride+1] != 175) || (buffer[i*stride+2] != 205)) {
        // rescale anything that isnt bahama blue...
        // (phi - phi_lo)/(phi_hi-phi_lo)*255 = gray
        double val     = double(buffer[i*stride+0])/255.0*(range[1]-range[0]) + range[0];
        double g_new   = min( 255.0, max( 0.0, (val - phi_lo)/(phi_hi-phi_lo) * 255.0));
        for (int j=0; j < stride; ++j)  
          buffer[i*stride+j] = (unsigned char)(g_new);
      }
    }

    // reset the range after the buffer is rescaled..
    range[0] = phi_lo;
    range[1] = phi_hi;
  }
};

void batch_rescale(const char* infile,const double phi_lo, const double phi_hi) { 

  ifstream fin(infile); 
  while (true) { 

    string img; fin >> img;  
    if ( fin.eof()) break;

    string fname;
    size_t pos   = img.find_last_of("/"); 
    if ( pos == string::npos) { 
      fname = img; 
    } else { 
      fname = img.substr(pos+1,img.length());
    }
    
    cout << "Working on: " << fname << endl;
    PngRescale png; 
    png.read(img.c_str(),true); 
    png.rescale_gray(phi_lo,phi_hi);
    string fout = "r_" + fname;
    png.write(fout.c_str()); 
  }

  fin.close(); 
}

int main(int argc, char* argv[] ) {

  if ( argc != 4) { 
    cout << "Usage: rescale.exe <image-list-text-file> <range-lo> <range-hi>" << endl;
    return 1;
  }

  batch_rescale(argv[1],atof(argv[2]),atof(argv[3])); 
  return 0; 
}




