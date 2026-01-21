#include <fstream>
#include <sstream>

// ==============================
// one dimensional profile reader
// desinged to read and store a 
// two colomn ASCII file
// ==============================
class oneDProfile {

private:
  int n;
  double (*data)[2];

public:

  oneDProfile() {    

    n = 0;
    data = NULL;

  }

  void read(const string& filename){

    if ( mpi_rank == 0 )
      cout << "reading 1D profile from file: " << filename << endl;
    // first read the file to get the size
    ifstream fp;
    fp.open(filename.c_str());
    if ( fp.fail()){
      cerr << "Error: could not find file " << filename << endl;
      throw(0);
    }
    string line;
    int linec = 0;
    while( !fp.eof() ){
      getline(fp,line);
      double data;
      stringstream ss(line);
      while(ss>> data ) {
	linec += 1;
      }
    }
    fp.close();
    n = linec/2;
    data = new double[n][2];
    ifstream fp2;
    fp2.open(filename.c_str());
    if ( fp2.fail()){
      cerr << "Error: could not find file " << filename << endl;
      throw(0);
    }
    int myn = 0;
    while( !fp2.eof() ){
      getline(fp2,line);
      stringstream ss(line);
      int i = 0;
      while(ss>> data[myn/2][i] ) {
	++i;
	++myn;
      }
    }
    fp2.close();
    assert(myn/2 == n);

    double min_y = 1000000000.0;
    double max_y = -1000000000.0;
    for ( int i=1; i<n; ++i ) {
      min_y = min(min_y,data[i][1]);
      max_y = max(max_y,data[i][1]);
    }
    if ( mpi_rank == 0 ) {
      cout << " > found    " << n << " points" << endl;
      cout << " > x range: (" << data[0][0] << "," << data[n-1][0] << ")" << endl;
      cout << " > y range: (" <<  min_y     << "," << max_y        << ")" << endl;
    }
  }


  double getData(const double x) {

    //if ( (x<data[0][0]) || (x>data[n-1][0]) ) {
    //  cerr << "Error in oneDProfile: data is out of bound " << endl;
    //  throw(0);
    //}

    if ( x <= data[0][0]   ) return(data[0][1]);
    if ( x >= data[n-1][0] ) return(data[n-1][1]);

    for ( int i=1; i<n; ++i ){
      if ( x <= data[i][0] ){
	assert(i>0);
	assert((data[i][0]-data[i-1][0])>0.0);
	double w = -(x-data[i][0])/(data[i][0]-data[i-1][0]);
	return( w*data[i-1][1]+(1.0-w)*data[i][1]);
      }
    }
  
    cerr << "Error: should not be here" << endl;
    throw(0);
  }

  double getMean(){

    double mean   = 0.0;
    double lenght = 0.0;
      for ( int i=1; i<n; ++i ) {
	mean   += 0.5 * (data[i][1] + data[i-1][1]) * (data[i][0] - data[i-1][0]);
	lenght += data[i][0] - data[i-1][0];
      }
      return(mean/lenght);
    
  }

  
  ~oneDProfile(){
    
    if ( data != NULL ) delete[] data;
    
  }

};
