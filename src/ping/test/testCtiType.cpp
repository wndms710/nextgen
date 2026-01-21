//#include "CTI.hpp"
//using namespace CTI;
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <assert.h>
#include <vector>

#include "PngImage.hpp"
#include "../CtiType.hpp"

/*
 * ***********************************************************
 *  Generate images using image_tools/charset/index.html. On
 *  this page, one sets a desired character height and the
 *  actual character height is measured.  The final reported
 *  height value includes a row of pixel white space at the
 *  top and bottom.  The reported width is accurate though
 *  several charters bleed over (#,w,W,R).  The image saved
 *  from the web browser is a PNG with colorspace srgb (one
 *  line of red is included at the end of the image to ensure
 *  it is not converted to pure grayscale).  For compatibility
 *  with PngImage read(), convert the colorspace to rgb i.e
 *  covert cti_type_20_09.png -colorspace cti_type_20_09-rgb.png.
 *  After this conversion the image can be read by this
 *  function.  Save the output to the header CtiType.hpp.
 *  Heights between 18px and 35px are included. This begins
 *  at the minimum Legend text size. Beyond 35px fonts are
 *  scaled by integer factors to get as close to the requested
 *  size without going over. i.e. the current default text
 *  size for a 2000x2000 image is 45px.  The 22px height font
 *  set will be scaled to 44px and used at that size in the
 *  legend box.
 * ***********************************************************
 */
void buildCtiTypeHpp() {
  int nfonts = 18;
  string imgstrs[18] = {"./cti_type_20_09-rgb.png",
                        "./cti_type_21_10-rgb.png",
                        "./cti_type_22_10-rgb.png",
                        "./cti_type_23_11-rgb.png",
                        "./cti_type_24_12-rgb.png",
                        "./cti_type_25_12-rgb.png",
                        "./cti_type_26_12-rgb.png",
                        "./cti_type_27_13-rgb.png",
                        "./cti_type_28_14-rgb.png",
                        "./cti_type_29_15-rgb.png",
                        "./cti_type_30_15-rgb.png",
                        "./cti_type_31_15-rgb.png",
                        "./cti_type_32_16-rgb.png",
                        "./cti_type_33_16-rgb.png",
                        "./cti_type_34_17-rgb.png",
                        "./cti_type_35_18-rgb.png",
                        "./cti_type_36_18-rgb.png",
                        "./cti_type_37_19-rgb.png"};
   
  int h0 = 18;
  int w_arr[18] = {9,10,10,11,12,12,12,13,14,15,15,15,16,16,17,18,18,19}; 
  
  cout << "#ifndef CTITYPE_HPP" << endl;
  cout << "#define CTITYPE_HPP" << endl;
  cout << "static int cti_type_n = " << nfonts << ";" << endl;
  cout << "static int cti_type_height[" << nfonts << "]" << " = {" << h0;
  for (int i = 1; i<nfonts; ++i) 
    cout << "," << h0+i; 
  cout <<  "};" << endl;
  cout << "static int cti_type_width[" << nfonts << "]" << " = {" << w_arr[0];
  for (int i = 1; i<nfonts; ++i)
    cout << "," << w_arr[i];
  cout <<  "};" << endl;
  cout << "static int cti_type_count = " << 126-32+1 << ";" << endl;
  
  for (int iI = 0; iI<nfonts; ++iI){
    int nci = w_arr[iI];
    int ncj = h0+iI;
    int i0 = nci;
    int j0 = 1;
  
    PngImage img;
    img.read(imgstrs[iI].c_str()); 
    
    cout << "static unsigned char cti_type_data_" << ncj << "[" << 126-32+1 << "]["<< nci*ncj <<"] = {" << "\n{";
    
    for (int ii = 32; ii <= 126; ++ii) {
  
      for (int j = 0; j < ncj; ++j) {
        for (int i = 0; i < nci; ++i) {
          unsigned char r,g,b;
          img.get(r,g,b,i+i0,j+j0);
          assert(r==g);
          assert(g==b);
                
          // evenly spaced so you can see it...
          //if (r <= 9)
          //  cout << "00" << int(r) << ",";
          //else if (r <= 99) 
          //  cout << "0" << int(r) << ",";
          //else
          //  cout << int(r) << ",";
                
          // or use this one when actually writing the CtiType.hpp,
          // because the leading zeros imply octal?...
          cout << int(r) << ",";      
  
        }
        cout << endl;
      }
      if (ii<126){
        cout << "},\n{";
      }
      else {
        cout << "}};" << endl;
      }
        
      //// check surroundings are 255...
      //bool problem = false;
      //for (int j = -1; j < ncj+1; ++j) {
      //  for (int i = -1; i < nci+1; ++i) {
      //    if ((i == -1)||(i == nci)||(j == -1)||(j == ncj)) {
      //      unsigned char r,g,b;
      //      img.get(r,g,b,i+i0,j+j0);
      //      if ((r!=255)||(g!=255)||(b!=255)) {
      //        cout << "surround problem: " << i << " " << j << endl;
      //        problem = true;
      //      }
      //    }
      //  }
      //}
      //if (problem)
      //  getchar();
        
      i0 += nci*2;
    }
  }
  cout << "static unsigned char * cti_type_data_ptr[" << nfonts << "]" << " = {cti_type_data_" << h0 << "[0]"; 
  for (int i = 1; i<nfonts; ++i)
    cout << ",cti_type_data_" << h0+i << "[0]";
  cout <<  "};" << endl;
  cout <<  "#endif" << endl;

}

void writeCtiTypeHpp(){
  cout << "writeCtiTypeHpp()..." << endl;
  int imghgt = 0;
  for (int iI=0;iI<cti_type_n;++iI)
    imghgt+=cti_type_height[iI];

  unsigned char (*buffer)[3] = new unsigned char[cti_type_width[cti_type_n-1]*95*imghgt][3];
  PngImage img(cti_type_width[cti_type_n-1]*95,imghgt,buffer);

  cout << "writing type..." << endl;
  int i = 0;
  int j = 0;
  for (int iI=0; iI<cti_type_n; ++iI){
    for (int iT = 0; iT < cti_type_count; ++iT) {
      for (int di = 0; di < cti_type_width[iI]; ++di) {
        for (int dj = 0; dj < cti_type_height[iI]; ++dj) {
          const unsigned char rgb = cti_type_data_ptr[iI][iT*cti_type_width[iI]*cti_type_height[iI]+dj*cti_type_width[iI]+di];
          img.set(i+di,j+dj,rgb,rgb,rgb);
        }
      }
      i += cti_type_width[iI];
    }
    i = 0;
    j += cti_type_height[iI];
  }

  img.write("test_type.png");

  delete[] buffer;
}
 
int main(int argc, char * argv[]) {
 
  if (argc < 2){
    cout << "Run testCtiType_s.exe with one of the following flags: " << endl;
    cout << " $ ./testCtiType_s.exe -b #build CtiType.hpp" << endl;
    cout << " $ ./testCtiType_s.exe -w #write test text from CtiType.hpp" << endl;
  }
  else if (strcmp(argv[1],"-b")==0){
    buildCtiTypeHpp();
  }
  else if (strcmp(argv[1],"-w")==0){
    writeCtiTypeHpp();
  }
  else{
    cout << "Invalid argument: " << argv[1] << endl;
    cout << "Run testCtiType_s.exe with one of the following flags: " << endl;
    cout << " $ ./testCtiType_s.exe -b #build CtiType.hpp" << endl;
    cout << " $ ./testCtiType_s.exe -w #write test text from CtiType.hpp" << endl;
  }

  return 0;
}


#ifdef testTypeJUNK
  void testType(Param * param) {
    
    cout << "testType()" << endl;
    
    // =============================================================
    // step 1: print the ASCII chars in some monospaced font...
    // =============================================================

    /*
      {
      cout << endl;
      for (int i = 32; i <= 126; ++i) {
      unsigned char v = i;
      cout << " " << v;
      }
      cout << endl;
      cout << endl;
      getchar();
      }
    */

    // =============================================================
    // step 2: Use this one to produce the CtiType.hpp header...
    // =============================================================

    
      {
      int nci = 18+2;
      int ncj = 35-2;
      int i0 = 17;
      int j0 = 1;
      img.read("test/cti_type_35_18-rgb.png"); 



      cout << img.getNx() << " " << img.getNy() << endl;

      cout << "static int cti_type_width = " << nci << ";" << endl;
      cout << "static int cti_type_height = " << ncj << ";" << endl;
      cout << "static int cti_type_count = " << 126-32+1 << ";" << endl;
      
      cout << "static char cti_type_data[" << 126-32+1 << "]["<< nci*ncj <<"] = {" <<
      "\n{";
      
      for (int ii = 32; ii <= 126; ++ii) {
	
        //cout << "i0,j0: " << i0 << " " << j0 << endl;

        for (int j = 0; j < ncj; ++j) {
          for (int i = 0; i < nci; ++i) {
            unsigned char r,g,b;
            img.get(r,g,b,i+i0,j+j0);
            assert(r==g);
            assert(g==b);
                  
            // evenly spaced so you can see it...
            if (r <= 9)
              cout << "00" << int(r) << ",";
            else if (r <= 99) 
              cout << "0" << int(r) << ",";
            else
              cout << int(r) << ",";
                  
            // or use this one when actually writing the CtiType.hpp,
            // because the leading zeros imply octal?...
            //cout << int(r) << ",";      

          }
          cout << endl;
        }
        cout << "},\n{";
          
        // check surroundings are 255...
        bool problem = false;
        for (int j = -1; j < ncj+1; ++j) {
          for (int i = -1; i < nci+1; ++i) {
            if ((i == -1)||(i == nci)||(j == -1)||(j == ncj)) {
              unsigned char r,g,b;
              img.get(r,g,b,i+i0,j+j0);
              if ((r!=255)||(g!=255)||(b!=255)) {
                cout << "surround problem: " << i << " " << j << endl;
                problem = true;
              }
            }
          }
        }
        if (problem)
          getchar();
          
        i0 += (nci-2)*2;
      }
      
      }
    }
    assert(0);
    
    
    
    // =============================================================
    // step 3: illustrate how the type works...
    // =============================================================
    
    {
    
      int imghgt = 0;
      for (int iI=0;iI<cti_type_n;++iI)
        imghgt+=cti_type_height[iI];

      PngData * img = new PngData(cti_type_width[cti_type_n-1]*95,imghgt);
    
      cout << "writing type..." << endl;
      int i = 0;
      int j = 0;
      for (int iI=0; iI<cti_type_n; ++iI){
        //static int cti_type_n = 16
        //static int cti_type_height[16] = {18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33};
        //static int cti_type_width[16] = {9,10,10,11,12,12,12,13,14,15,15,15,16,16,17,18};
        //static int cti_type_count = 95;


        for (int iT = 0; iT < cti_type_count; ++iT) {
          for (int di = 0; di < cti_type_width[iI]; ++di) {
            for (int dj = 0; dj < cti_type_height[iI]; ++dj) {
              const unsigned char rgb = cti_type_data_ptr[iI][iT*cti_type_width[iI]*cti_type_height[iI]+dj*cti_type_width[iI]+di];
              img->set(i+di,j+dj,rgb,rgb,rgb);
            }
          }
          i += cti_type_width[iI];
        }
        i = 0;
        j += cti_type_height[iI];
      }
    
      img->write("test_type.png");

      delete img;
    
    }
    
    cout << " > look at test_type.png" << endl;

  }

  void addBlackTextLeftJustified(PngData * img,const string& text,const int i0,const int j0) {
    
    // TODO: should check about overwriting image limits. For now, just assume we won't...
    // also, we are assuming there are no specuial chars in text (e.g. tab, newline, etc)...
    
    int i = i0 - cti_type_width*text.length();
    int j = j0 - cti_type_height/2;
    for (int ii = 0; ii < text.length(); ++ii) {
      const int ichar = int(text[ii])-32;
      assert((ichar >= 0)&&(ichar < cti_type_count));
      for (int di = 0; di < cti_type_width; ++di) {
	for (int dj = 0; dj < cti_type_height; ++dj) {
	  unsigned char r,g,b;
	  img->get(r,g,b,i+di,j+dj);
	  const unsigned char rgb = cti_type_data[ichar][dj*cti_type_width+di];
	  unsigned char r2 = (unsigned char)(double(rgb)/255.0*double(r));
	  unsigned char g2 = (unsigned char)(double(rgb)/255.0*double(g));
	  unsigned char b2 = (unsigned char)(double(rgb)/255.0*double(b));
	  img->set(i+di,j+dj,r2,g2,b2);
	}
      }
      i += cti_type_width;
    }
    
  }
#endif


