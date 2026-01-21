#ifndef PDFCONVOLUTION_HPP
#define PDFCONVOLUTION_HPP 

#include "Common.hpp"
#include "CTI.hpp" 
#include "MiscUtils.hpp" 

using namespace CTI;
using namespace MiscUtils;


namespace PdfConvolution { 
  
  enum PDFType { 
    BETA_PDF , 
    UNKNOWN_PDF
  }; 


  extern double * x_pdf ;    // x-values of the pdf, arg P(x)  
  extern int n_pts ; 
  extern string x_name ;     // name of the indep coordinate 
  extern int * idx ;         // sorted indicies of indep coordinate 
  extern double * pdf ;      // P(x) 
  extern double mean ; 
  extern double var ; 
  extern double invLx ; 
  extern double xmin ; 
  extern PDFType pdfType ; 

  void buildIndepCoord( int _n, const double * _x, const char * _name, const double _xmin, const double _xmax) ; 
  void updateBetaPDF(const string& name, const double _mean, const double _var) ; 
  double convolveWithPdf(const string& varName, const double * var , const int _n) ;  
  void updatePdf(const string& name, const double _mean, const double _var) ; 
  void clear() ; 


}; 






#endif 
