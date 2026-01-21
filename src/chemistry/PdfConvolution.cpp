#include "PdfConvolution.hpp" 

namespace PdfConvolution { 

  double * x_pdf = NULL;    // x-values of the pdf, arg P(x) 
  int n_pts      = 0   ; 
  string x_name ;           // name of the indep coordinate 
  int * idx      = NULL;    // sorted indicies of indep coordinate 
  double * pdf   = NULL;    // P(x) 
  double mean    = 0.0 ; 
  double var     = -1.0 ;   // invalid initial state ..  
  double invLx   = -1.0 ; 
  double xmin    =  0.0 ; 
  PDFType pdfType = BETA_PDF ; 

  void buildIndepCoord( const int _n, const double * _x, const char * _name, const double _xmin, const double _xmax) { 
    x_name = string(_name); 
    resize ( x_pdf, n_pts, _n) ; 
    resize ( idx  , n_pts, _n) ; 
    resize ( pdf  , n_pts, _n) ; 
    n_pts  = _n ; 


    // 
    // normalize the xcoordinate so that its [0,1]
    // 
    invLx         = 1.0/ (_xmax - _xmin) ; 
    xmin          = _xmin ; 
    bool isSorted = true ; 


    for (int i = 0 ; i < n_pts ; ++i) { 
      x_pdf[i] = (_x[i] - _xmin) * invLx ; 
      idx[i]   = i ; //initial ordering

      if ( i >= 1 )  
	isSorted &= (x_pdf[i] >= x_pdf[i-1]); 
    }


    // 
    // now sort the coordinate by value (store in idx) 
    // 
    if ( !isSorted) 
      quickSortByValue (idx, x_pdf, n_pts) ; 
 

  }//buildIndepCoords




  // 
  // this is roughly the same version that was coded 
  // in the original Flamelet.hpp 
  // 

  void updateBetaPDF(const string& name, const double _mean, const double _var) { 
    
    // rescale .. 
    double tmpMean = (_mean - xmin) * invLx ; 
    double tmpVar  = _var * invLx * invLx ; 


    // check if state has changed .. 
    const double TOL = 1.0e-10; 
    //if ( fabs(mean-tmpMean) > TOL || fabs(var-tmpVar) > TOL || 
    //	 x_name != name ) { 
    if ( true) { 
      
      mean = tmpMean ; 
      var  = tmpVar  ; 
      assert ( n_pts > 0 ) ; 

      for (int i = 0 ; i < n_pts ; ++i) 
	pdf[i] = 0.0 ; 
      
      
      if ( mean < TOL ) {
	
	pdf[idx[0]] = 1.0;
	return;
	
      }
      else if ( (1.0-mean) < TOL ) {
	
	pdf[idx[n_pts-1]] = 1.0;
	return;
	
      }      
      else if ( var <= TOL ) {
	
	int count = 0;
	for (int i=0 ; i<n_pts-1; ++i) { 
	  int ii = idx[i+1]; 
	  if ( mean < x_pdf[ii] ) 
	    break; 
	  else 
	    ++count ; 
	}


	if ( count == n_pts-1 ) 
	  CERR("mean is not in the x_pdf range") ; 
	
	assert ( x_pdf[idx[count+1]] > x_pdf[idx[count]] ); 
	assert ( mean >= x_pdf[idx[count]] ) ; 


	double wgt = ( x_pdf[idx[count+1]] - mean ) / 
	  ( x_pdf[idx[count+1]] - x_pdf[idx[count]] ) ; 

	pdf[idx[count]]    = wgt ; 
	pdf[idx[count+1]]  = 1.0- wgt; 
	assert ( wgt <= 1.0 && wgt >= 0.0); 
	return ; 
      }
      // Impossible cases (i.e., var > mean*(1-mean)): two delta pdf at x_pdf=0 and x_pdf=1 
      // note that we already checked for  var < mean*(1-mean)   
      // XXXX why do we treat impossible cases??? 
      else if ( var >= mean * (1.0 - mean)) {

	pdf[idx[0]]       = 1.0 - mean ; 
	pdf[idx[n_pts-1]] = mean ; 
	return; 
      }
      else {
	
	double a      = mean * (mean * (1.0 - mean) / var - 1.0);
	double b      = a / mean - a;

	assert( a != 0.0 );
	assert( b != 0.0 );
	double factor = lgamma(a+b) - lgamma(a) - lgamma(b);

	// Left BC: explicit integration
	double dx = 0.5 * (x_pdf[idx[1]] - x_pdf[idx[0]] ) ; 
	if ( dx < TOL ) 
	  pdf[idx[0]] = 0.0;
	else {
	  double tmp = a * log(dx) + factor;
	  pdf[idx[0]] = exp(tmp) / a;
	}
	// Right BC: explicit integration
	dx = 0.5 * (x_pdf[idx[n_pts-1]] - x_pdf[idx[n_pts-2]] ); 
	if ( dx < TOL ) 
	  pdf[idx[n_pts-1]] = 0.0;
	else {
	  double tmp               = b * log(dx) + factor;
	  pdf[idx[n_pts-1]] = exp(tmp) / b;
	}
	// other points
	for (int i=1; i<n_pts-1; ++i) {
	  dx     = 0.5 * (x_pdf[idx[i+1]] - x_pdf[idx[i-1]] ) ; 
	  if ( dx < TOL )
	    pdf[idx[i]] = 0.0;
	  else {
	    double tmp = (a - 1.0) * log(x_pdf[idx[i]]) + 
	      (b - 1.0) * log(1.0 - x_pdf[idx[i]]) + factor;
	    pdf[idx[i]] = exp(tmp) * dx;
	  }
	}

      }

      // normalize pdf
      double sum = 0.0;
      for  (int i=0; i<n_pts; ++i)
        sum += pdf[i];
      assert(sum>0.0);
      double one_over_sum = 1.0 / sum;
      for  (int i=0; i<n_pts; ++i)
        pdf[i] *= one_over_sum;
      
    } //  if ((mean != pdfMean) || (var != pdfVar))
  
  }//updateBetaPDF

  
  double convolveWithPdf(const string& varName, const double * var , const int _n) { 
    assert ( _n == n_pts) ; 
    double sum = 0.0 ; 

    if ( varName == "rho") { 
      for (int i = 0 ; i < n_pts ; ++i) 
	sum += pdf[i] / var[i] ; 
      return 1.0/sum ; 
    }
    else { 
      for (int i = 0 ; i < n_pts ; ++i) 
	sum += pdf[i] * var[i] ; 
      return sum ;  
    } 
    
  }//convolvewithpdf


  
  
  void updatePdf(const string& name, const double _mean, const double _var) { 
    
    switch ( pdfType) { 
    case BETA_PDF : 
      updateBetaPDF(name,_mean,_var); 
      break;  

    default: 
      assert(0); 
    }

  }//updatePdf 




  // 
  // in case you need to bring down the pdf 
  // interface and bring it back up .. 
  // 
  void clear() { 
    DELETE( x_pdf); x_pdf = NULL; 
    DELETE( pdf ) ; pdf = NULL; 
    DELETE( idx)  ; idx = NULL; 
  }




}; 







