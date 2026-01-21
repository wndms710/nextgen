
#include "MiscUtils.hpp"
#include "tomcrypt.hpp"
#include <fcntl.h>

namespace MiscUtils {

  void nullTerminate(char *name) {
    // f90 strings stored in binary files are space-terminated. This
    // adds a NULL termination after the last non-space character.
    int i = strlen(name)-1;

    while ((i > 0) && (name[i-1] == ' ')) {
      i--;
    }
    name[i] = '\0';
  }

  void nullTerminate(char *name, int len) {
    // f90 strings stored in binary files are space-terminated. This
    // adds a NULL termination after the last non-space character.
    int i = len;

    while ((i > 0) && (name[i-1] == ' ')) {
      i--;
    }
    name[i] = '\0';
  }

  string toLowerCase(const string word) {
    // copy string first and replace inline; assures memory is allocated
    // transform doesn't necessarily push properly into output string if empty
    string word_lc = word;
    std::transform(word_lc.begin(),word_lc.end(),word_lc.begin(),::tolower);
    return word_lc;
  }

  string toUpperCase(const string word) {
    // copy string first and replace inline; assures memory is allocated
    // transform doesn't necessarily push properly into output string if empty
    string word_uc = word;
    std::transform(word_uc.begin(),word_uc.end(),word_uc.begin(),::toupper);
    return word_uc;
  }

  void calcUniformDist(int * &xod, const int nx, const int ndist) {
    if (xod == NULL) xod = new int[ndist+1];
    xod[0] = 0;
    for (int id = 1; id <= ndist; id++) {
      xod[id] = (int)((double)id/(double)ndist*(double)nx + 0.5);
      assert(xod[id] >= xod[id-1]);
    }
    assert(xod[ndist] == nx);
  }

  void calcUniformDist(int8 * &xod, const int8 nx, const int ndist) {
    if  (xod == NULL) xod = new int8[ndist+1];
    xod[0] = 0;
    for (int id = 1; id <= ndist; ++id) {
      xod[id] = (int8)((double)id/(double)ndist*(double)nx + 0.5);
      assert(xod[id] >= xod[id-1]);
    }
    assert(xod[ndist] == nx);
    // make sure no local count exceeds TWO_BILLION...
    for (int id = 0; id < ndist; ++id) {
      assert(xod[id+1]-xod[id] < TWO_BILLION);
    }
  }

  void calcUniformDistAlign8(int8 * &xod, const int8 nx, const int ndist) {

    // this version builds a uniform dist and ensures alignment of all
    // offsets along a multiple of 8. It is used in stitch's nodal
    // periodicity calculations...

    if  (xod == NULL) xod = new int8[ndist+1];
    xod[0] = 0;
    for (int id = 1; id <= ndist; ++id) {
      xod[id] = (int8)((double)id/(double)ndist*(double)nx + 0.5);
      const int8 remainder = xod[id]%8;
      if (remainder) {
	xod[id] += 8 - remainder;
	assert(xod[id]%8 == 0);
      }
      xod[id] = min(xod[id],nx);
    }
    assert(xod[ndist] == nx);

    // make sure no local count exceeds TWO_BILLION...
    for (int id = 0; id < ndist; ++id) {
      assert(xod[id+1]-xod[id] < TWO_BILLION);
    }

  }

  void calcThresholdDist(int* &xod, const int nx, const int ndist, const int nmin) {
    if ( xod == NULL) xod = new int[ndist+1];
    int nunif = nx/ndist ;
    if ( nunif >= nmin ) calcUniformDist(xod,nx,ndist) ;
    else {
      int tmp = 0 ;
      int ix  = 0 ;
      while ( tmp < nx ) {
        assert ( ix < ndist ) ;
        xod[ix++] = tmp ;
        tmp      += nmin ;
      }
      for (int ii = ix ; ii <= ndist ; ++ii)
        xod[ii] = nx ;
    }
  }

  void calcThresholdDist(int8* &xod, const int8 nx, const int ndist, const int nmin) {
    if ( xod == NULL) xod = new int8[ndist+1];
    int nunif = nx/ndist ;
    if ( nunif >= nmin ) calcUniformDist(xod,nx,ndist) ;
    else {
      int8 tmp = 0 ;
      int ix  = 0 ;
      while ( tmp < nx ) {
        assert ( ix < ndist ) ;
        xod[ix++] = tmp ;
        tmp      += nmin ;
      }
      for (int ii = ix ; ii <= ndist ; ++ii)
        xod[ii] = nx ;
    }
  }

  void calcUniformDistNoNew(int8 * xod, const int8 nx, const int ndist) {
    assert(xod != NULL);
    xod[0] = 0;
    for (int id = 1; id <= ndist; id++) {
      xod[id] = (int8)((double)id/(double)ndist*(double)nx + 0.5);
    }
  }

  void buildUniformXora(int * &Xora, const int nX_global) {
    calcUniformDist(Xora,nX_global,mpi_size);
  }

  void buildUniformXora(int8 * &Xora, const int8 nX_global) {
    calcUniformDist(Xora,nX_global,mpi_size);
  }

  /*
    void buildUniformXora(int64_t * &Xora, const int64_t nX_global) {
    calcUniformDist(Xora,nX_global,mpi_size);
    }
  */

  void buildXora(int * &Xora,int nX_local) {

    // check...
    int8 nX_i8 = (int8)nX_local;
    int8 nX_global_i8;
    MPI_Allreduce(&nX_i8,&nX_global_i8,1,MPI_INT8,MPI_SUM,mpi_comm);
    if (nX_global_i8 > TWO_BILLION)
      CERR("Error: nX_global_i8 > TWO_BILLION: " << nX_global_i8);

    if (Xora == NULL) Xora = new int[mpi_size+1];
    MPI_Allgather(&nX_local,1,MPI_INT,Xora+1,1,MPI_INT,mpi_comm);
    Xora[0] = 0;
    for (int i = 0; i < mpi_size; i++)
      Xora[i+1] += Xora[i];
  }

  void buildXora(int8 * &Xora,int nX_local) {

    if (Xora == NULL) Xora = new int8[mpi_size+1];
    int8 nX_i8 = (int8)nX_local;
    MPI_Allgather(&nX_i8,1,MPI_INT8,Xora+1,1,MPI_INT8,mpi_comm);
    Xora[0] = 0;
    for (int i = 0; i < mpi_size; i++)
      Xora[i+1] += Xora[i];

  }

  void buildXora(int8 * &Xora,int8 nX_i8) {

    if (Xora == NULL) Xora = new int8[mpi_size+1];
    MPI_Allgather(&nX_i8,1,MPI_INT8,Xora+1,1,MPI_INT8,mpi_comm);
    Xora[0] = 0;
    for (int i = 0; i < mpi_size; i++)
      Xora[i+1] += Xora[i];

  }

  void copyXora(int * &xora,const int8 * xora_global) {

    if (xora == NULL) xora = new int[mpi_size+1];
    assert(xora_global[0] == 0);
    assert(xora_global[mpi_size] < TWO_BILLION);
    xora[0] = 0;
    FOR_RANK xora[rank+1] = (int)xora_global[rank+1];

  }

  int getRankInXora(const int ix,const int * xora) {

    // bracket...
    int left_rank = 0;
    int right_rank = mpi_size;
    while ((right_rank - left_rank) > 1) {
      const int middle_rank = (left_rank + right_rank)/2;   // equivalent to floor..
      if (ix >= xora[middle_rank])
	left_rank = middle_rank;
      else
	right_rank = middle_rank;
    }
    //assert(left_rank == rank);
    return(left_rank);

  }

  int getIntervalBisection(const double value,const double * const valueArray,const int n) {

    // double valueArray[n] contains monotonically increasing values. When value
    // is considered inside interval "i" when
    // value >= valueArray[i], and
    // value < valueArray[i+1]
    // ...

    assert(n > 1);
    assert(value >= valueArray[0]);
    assert(value <= valueArray[n-1]);

    int left = 0;
    int right = n-1;
    while ((right - left) > 1) {
      const int middle = (left + right)/2;   // equivalent to floor..
      if (value >= valueArray[middle])
	left = middle;
      else
	right = middle;
    }
    return(left);

  }

  int getIdxInXoso(const int ix,const int * xoso,const int nn) {

    // bracket...
    int left_rank = 0;
    int right_rank = nn;
    while ((right_rank - left_rank) > 1) {
      const int middle_rank = (left_rank + right_rank)/2;   // equivalent to floor..
      if (ix >= xoso[middle_rank])
	left_rank = middle_rank;
      else
	right_rank = middle_rank;
    }
    //assert(left_rank == rank);
    return(left_rank);

  }

  int getRankInXora(const int8 ix,const int8 * xora) {

    // x-of-rank defines the range on each rank of a stiped variable x
    // this should be written faster some day: Paul...

    /*
      static bool first = true;
      if (first) {
      if (mpi_rank == 0)
      cout << "XXXXXXX getRankInXora should be removed eventually" << endl;
      first = false;
      }

      assert(xora[0] == 0);

      assert((ix >= xora[0])&&(ix < xora[mpi_size]));
      int rank = 0;
      while (ix >= xora[rank+1])
      ++rank;

      //return(rank);
      */

    // bracket...
    int left_rank = 0;
    int right_rank = mpi_size;
    while ((right_rank - left_rank) > 1) {
      const int middle_rank = (left_rank + right_rank)/2;   // equivalent to floor..
      if (ix >= xora[middle_rank])
	left_rank = middle_rank;
      else
	right_rank = middle_rank;
    }
    //assert(left_rank == rank);
    return(left_rank);

  }

  void dumpRange(const int * s, const int n, const string& message) {

    int my_buf[2] = { 1000000000, 1000000000 };
    for (int i = 0; i < n; ++i) {
      my_buf[0] = min(my_buf[0],s[i]);
      my_buf[1] = min(my_buf[1],-s[i]);
    }

    int buf[2];
    MPI_Reduce(my_buf, buf, 2, MPI_INT, MPI_MIN, 0, mpi_comm);
    if (mpi_rank == 0)
      cout << " > dumpRange: "<< message << ", "<< buf[0]<< ":"<< -buf[1]<< endl;

  }

  void dumpRange(const int8 * s, const int n, const string& message) {

    int8 my_buf[2] = { (int8(1)<<60), (int8(1)<<60) };
    for (int i = 0; i < n; ++i) {
      my_buf[0] = min(my_buf[0],s[i]);
      my_buf[1] = min(my_buf[1],-s[i]);
    }

    int8 buf[2];
    MPI_Reduce(my_buf, buf, 2, MPI_INT8, MPI_MIN, 0, mpi_comm);
    if (mpi_rank == 0)
      cout << " > dumpRange: "<< message << ", "<< buf[0]<< ":"<< -buf[1]<< endl;

  }

  void dumpRange(const int (*v)[3], const int n, const string& message) {

    int my_buf[6];
    for (int j = 0; j < 6; j++)
      my_buf[j] = ~(1<<31); // something large

    for (int i = 0; i < n; i++) {
      for (int j = 0; j < 3; j++) {
        my_buf[2*j] = min(my_buf[2*j], v[i][j]);
        my_buf[2*j+1] = min(my_buf[2*j+1], -v[i][j]);
	// nan check...
	if ( v[i][j] != v[i][j] ) {
	  cerr << "Error: vector " << message << " failed NAN check." << endl;
	  throw(-1);
	}
      }
    }

    int buf[6];
    MPI_Reduce(my_buf, buf, 6, MPI_INT, MPI_MIN, 0, mpi_comm);
    if (mpi_rank == 0) {
      cout << " > dumpRange: "<< message;
      for (int j = 0; j < 3; j++) {
        cout << ", "<< j << ": "<< buf[2*j]<< ":"<< -buf[2*j+1];
      }
      cout << endl;
    }

    /*
    // HACK: also dump a check-sum-like quantity...
    {
    FOR_J3 {
    assert( mpi_size == 1 );
    double check = 0.0;
    for (int i = 0; i < n; i++) check += v[i][j]*v[i][j];
    cout << " > check: " << check << endl;
    }
    }
    */

  }

  void dumpRange(const double * s, const int n,const string& message,MPI_Comm& comm) {

    double my_buf[2];
    for (int j = 0; j < 2; j++)
      my_buf[j] = HUGE_VAL; // something large

    for (int i = 0; i < n; i++) {
      my_buf[0] = min(my_buf[0], s[i]);
      my_buf[1] = min(my_buf[1], -s[i]);
      // nan check...
      if ( s[i] != s[i] ) {
	cerr << "Error: scalar " << message << " failed NAN check." << endl;
	throw(-1);
      }
    }

    double buf[2];
    MPI_Reduce(my_buf, buf, 2, MPI_DOUBLE, MPI_MIN, 0, comm);
    if (mpi_rank == 0)
      cout << " > dumpRange: "<< message << ", "<< buf[0]<< ":"<< -buf[1]<< endl;

    /*
    // HACK: also dump a check-sum-like quantity...
    {
    assert( mpi_size == 1 );
    double check = 0.0;
    for (int i = 0; i < n; i++) check += s[i]*s[i];
    cout << " > check: " << check << endl;
    }
    */

  }

  void dumpRange(const double (*v)[3], const int n, const string& message,MPI_Comm& comm) {

    double my_buf[6];
    for (int j = 0; j < 6; j++)
      my_buf[j] = HUGE_VAL; // something large

    for (int i = 0; i < n; i++) {
      for (int j = 0; j < 3; j++) {
        my_buf[2*j] = min(my_buf[2*j], v[i][j]);
        my_buf[2*j+1] = min(my_buf[2*j+1], -v[i][j]);
	// nan check...
	if ( v[i][j] != v[i][j] ) {
	  cerr << "Error: vector " << message << " failed NAN check." << endl;
	  throw(-1);
	}
      }
    }

    double buf[6];
    MPI_Reduce(my_buf, buf, 6, MPI_DOUBLE, MPI_MIN, 0, comm);
    if (mpi_rank == 0) {
      cout << " > dumpRange: "<< message;
      for (int j = 0; j < 3; j++) {
        cout << ", "<< j << ": "<< buf[2*j]<< ":"<< -buf[2*j+1];
      }
      cout << endl;
    }

    /*
    // HACK: also dump a check-sum-like quantity...
    {
    FOR_J3 {
    assert( mpi_size == 1 );
    double check = 0.0;
    for (int i = 0; i < n; i++) check += v[i][j]*v[i][j];
    cout << " > check: " << check << endl;
    }
    }
    */
  }

  void dumpRange(const float * s, const int n, const string& message, MPI_Comm& comm) {

    float my_buf[2];
    for (int j = 0; j < 2; j++)
      my_buf[j] = HUGE_VALF; // something large

    for (int i = 0; i < n; i++) {
      my_buf[0] = min(my_buf[0], s[i]);
      my_buf[1] = min(my_buf[1], -s[i]);
      // nan check...
      if ( s[i] != s[i] ) {
        cerr << "Error: scalar " << message << " failed NAN check." << endl;
        throw(-1);
      }
    }

    float buf[2];
    MPI_Reduce(my_buf, buf, 2, MPI_FLOAT, MPI_MIN, 0, comm);
    if (mpi_rank == 0)
      cout << " > dumpRange: "<< message << ", "<< buf[0]<< ":"<< -buf[1]<< endl;

    /*
    // HACK: also dump a check-sum-like quantity...
    {
    assert( mpi_size == 1 );
    double check = 0.0;
    for (int i = 0; i < n; i++) check += s[i]*s[i];
    cout << " > check: " << check << endl;
    }
    */

  }

  void setRange(const double (*v)[3], const int n, double (&buf)[3][2]) {

    double my_buf[6], g_buf[6];
    for (int j = 0; j < 6; j++)
      my_buf[j] = HUGE_VAL; // something large

    for (int i = 0; i < n; i++) {
      for (int j = 0; j < 3; j++) {
        my_buf[2*j]   = min(my_buf[2*j], v[i][j]);
        my_buf[2*j+1] = min(my_buf[2*j+1], -v[i][j]);
      }
    }

    MPI_Allreduce(my_buf, g_buf, 6, MPI_DOUBLE, MPI_MIN, mpi_comm);

    for(int i=0; i<3; ++i)
      for(int j=0; j<2; ++j)
        buf[i][j] = (j==1? -1.0 : 1.0)*g_buf[i*2+j];
  }

  void dumpRange(const double (*t)[3][3], const int n, const string& message) {

    double my_buf[18];
    for (int j = 0; j < 18; j++)
      my_buf[j] = HUGE_VAL; // something large

    for (int i = 0; i < n; i++) {
      for (int j = 0; j < 3; j++) {
	for (int k = 0; k < 3; k++) {
	  my_buf[2*(3*j+k)] = min(my_buf[2*(3*j+k)], t[i][j][k]);
	  my_buf[2*(3*j+k)+1] = min(my_buf[2*(3*j+k)+1], -t[i][j][k]);
	  // nan check...
	  assert( t[i][j][k] == t[i][j][k] );
	}
      }
    }

    double buf[18];
    MPI_Reduce(my_buf, buf, 18, MPI_DOUBLE, MPI_MIN, 0, mpi_comm);
    if (mpi_rank == 0) {
      for (int j = 0; j < 3; j++) {
	cout << " > dumpRange: "<< message;
	for (int k = 0; k < 3; k++) {
	  cout << ", " << j << "," << k << ": "<< buf[2*(3*j+k)]<< ":"<< -buf[2*(3*j+k)+1];
	}
	cout << endl;
      }
    }

  }

  void dumpRange(const int value, const string& message) {

    int my_minmax[2] = { value, -value };
    int minmax[2];
    MPI_Reduce(my_minmax,minmax,2,MPI_INT,MPI_MIN,0,mpi_comm);

    int8 my_sum = value;
    int8 sum;
    MPI_Reduce(&my_sum,&sum,1,MPI_INT8,MPI_SUM,0,mpi_comm);

    if (mpi_rank == 0) cout << " > dumpRange: " << message << ", min,avg,max: " << minmax[0] << " " << sum/mpi_size << " " << -minmax[1] << endl;

  }

  void dumpHist0(const double * var,const int n,const string& message) {

    assert(mpi_rank == 0);

    // determine the range...

    assert(n > 0);
    double var_min = var[0];
    double var_max = var[0];
    for (int i = 0; i < n; ++i) {
      assert( var[i] == var[i] ); // nan check
      var_min = min(var_min,var[i]);
      var_max = max(var_max,var[i]);
    }

    // expand very slightly...

    double eps = 1.0E-6*(var_max-var_min);
    var_min -= eps;
    var_max += eps;

    // build historgram...

    const int nbin = 120;
    int count[nbin];
    for (int ib = 0; ib < nbin; ++ib)
      count[ib] = 0;

    for (int i = 0; i < n; ++i) {
      int ib = (int)((double)nbin*(var[i]-var_min)/(var_max-var_min));
      assert((ib >= 0)&&(ib < nbin));
      count[ib] += 1;
    }

    int count_max = 0;
    for (int ib = 0; ib < nbin; ++ib)
      count_max = max(count_max,count[ib]);

    // and print it...

    cout << endl;
    for (int ib = 0; ib < nbin; ++ib)
      cout << "-";
    cout << endl;
    cout << "Historgram: " << message << ", range: " << var_min << " " << var_max << ", samples: " << n << endl;
    int nrows = 40;
    for (int ir = nrows-1; ir >= 0; --ir) {
      for (int ib = 0; ib < nbin; ++ib) {
	if (count[ib]*nrows >= count_max*ir)
	  cout << "*";
	else
	  cout << " ";
      }
      cout << endl;
    }
    for (int ib = 0; ib < nbin; ++ib)
      cout << "-";
    cout << endl;

    /*
      cout << "press any key for tecplot stuff" << endl;
      getchar();
      for (int ib = 0; ib < nbin; ++ib) {
      double x = ( var_min*((double)ib+0.5) + var_max*((double)(nbin-ib)+0.5) )/(double)nbin;
      cout << "HISTOGRAM " << x << " " <<
      (double)count[ib]/(double)count_max*0.4 << " " <<
      exp(-x*x/2.0)/sqrt(2.0*M_PI) << endl;
      }
      throw(0);
    */

  }

  void dumpHist0(const double * var,const int n,const int nbin,const string& message, vector<double>& x, vector<double>& y, bool verbose) {

    assert(mpi_rank == 0);

    // determine the range...

    assert(n > 0);
    double var_min = var[0];
    double var_max = var[0];
    for (int i = 0; i < n; ++i) {
      assert( var[i] == var[i] ); // nan check
      var_min = min(var_min,var[i]);
      var_max = max(var_max,var[i]);
    }

    // expand very slightly...

    //double eps = 1.0E-6*(var_max-var_min);
    double eps = 1.0E-10;
    var_min -= eps;
    var_max += eps;

    // build historgram...

    //const int nbin = 120;
    int count[nbin];
    for (int ib = 0; ib < nbin; ++ib)
      count[ib] = 0;

    for (int i = 0; i < n; ++i) {
      int ib = (int)((double)nbin*(var[i]-var_min)/(var_max-var_min));
      assert((ib >= 0)&&(ib < nbin));
      count[ib] += 1;
    }

    assert(x.size()==0);
    assert(y.size()==0);
    int sum_count = 0;
    for (int ib = 0; ib < nbin; ++ib) {
      x.push_back(var_min + double(ib)*(var_max-var_min)/double(nbin));
      sum_count += count[ib];
    }
    double dx = x[1] - x[0];
    for (int ib = 0; ib < nbin; ++ib) {
      y.push_back(double(count[ib])/double(sum_count)/dx);
    }

    if (not verbose) return;

    int count_max = 0;
    for (int ib = 0; ib < nbin; ++ib)
      count_max = max(count_max,count[ib]);

    // and print it...

    cout << endl;
    for (int ib = 0; ib < nbin; ++ib)
      cout << "-";
    cout << endl;
    cout << "Historgram: " << message << ", range: " << var_min << " " << var_max << ", samples: " << n << ", nbins: " << nbin << endl;
    int nrows = 40;
    for (int ir = nrows-1; ir >= 0; --ir) {
      for (int ib = 0; ib < nbin; ++ib) {
	if (count[ib]*nrows >= count_max*ir)
	  cout << "*";
	else
	  cout << " ";
      }
      cout << endl;
    }
    for (int ib = 0; ib < nbin; ++ib)
      cout << "-";
    cout << endl;

    /*
      cout << "press any key for tecplot stuff" << endl;
      getchar();
      for (int ib = 0; ib < nbin; ++ib) {
      double x = ( var_min*((double)ib+0.5) + var_max*((double)(nbin-ib)+0.5) )/(double)nbin;
      cout << "HISTOGRAM " << x << " " <<
      (double)count[ib]/(double)count_max*0.4 << " " <<
      exp(-x*x/2.0)/sqrt(2.0*M_PI) << endl;
      }
      throw(0);
    */

  }

  void dumpSum(const int * s, const int n, const string& message) {

    int my_buf = 0;
    for (int i = 0; i < n; i++) {
      // nan check...
      if ( s[i] != s[i] ) {
	cerr << "Error: scalar " << message << " failed NAN check." << endl;
	throw(-1);
      }
      my_buf += s[i];
    }

    int buf;
    MPI_Reduce(&my_buf, &buf, 1, MPI_INT, MPI_SUM, 0, mpi_comm);
    if (mpi_rank == 0)
      cout << message << ": " << buf << endl;

  }

  void dumpSum(const double * s, const int n, const string& message) {

    double my_buf = 0.0;
    for (int i = 0; i < n; i++) {
      // nan check...
      if ( s[i] != s[i] ) {
	cerr << "Error: scalar " << message << " failed NAN check." << endl;
	throw(-1);
      }
      my_buf += s[i];
    }

    double buf;
    MPI_Reduce(&my_buf, &buf, 1, MPI_DOUBLE, MPI_SUM, 0, mpi_comm);
    if (mpi_rank == 0)
      cout << message << ": " << buf << endl;

  }

  void dumpSum(const double (*v)[3], const int n, const string& message) {

    double my_buf[3] = { 0.0, 0.0, 0.0 };
    for (int i = 0; i < n; i++) {
      // nan check...
      FOR_J3 if ( v[i][j] != v[i][j] ) {
	cerr << "Error: vector " << message << " failed NAN check." << endl;
	throw(-1);
      }
      FOR_J3 my_buf[j] += v[i][j];
    }

    double buf[3];
    MPI_Reduce(my_buf, buf, 3, MPI_DOUBLE, MPI_SUM, 0, mpi_comm);
    if (mpi_rank == 0)
      cout << message << ": " << buf[0] << " " << buf[1] << " " << buf[2] << endl;
  }

  void dumpBins(const int * s, const int n, const string& message) {

    int my_buf[2];
    for (int j = 0; j < 2; j++)
      my_buf[j] = 1000000000; // something large

    for (int i = 0; i < n; i++) {
      my_buf[0] = min(my_buf[0], s[i]);
      my_buf[1] = min(my_buf[1], -s[i]);
    }

    int buf[2];
    MPI_Allreduce(my_buf, buf, 2, MPI_INT, MPI_MIN, mpi_comm);

    int nbins = -buf[1] - buf[0]+ 1;
    //assert( (nbins > 0)&&(nbins < 100));

    int * my_bin_count = new int[nbins];
    for (int i = 0; i < nbins; i++)
      my_bin_count[i] = 0;

    for (int i = 0; i < n; i++)
      my_bin_count[s[i]-buf[0]] += 1;

    int * bin_count = new int[nbins];
    MPI_Reduce(my_bin_count, bin_count, nbins, MPI_INT, MPI_SUM, 0, mpi_comm);
    if (mpi_rank == 0) {
      cout << " > dumpScalarBins: " << message << endl;
      for (int i = 0; i < nbins; i++)
        cout << "   > value: "<< i+buf[0] << ", count: "<< bin_count[i] << endl;
    }

    delete[] my_bin_count;
    delete[] bin_count;

  }

  void dumpNodeAvgReport(double value,const string& message) {

    double node_data[3];
    MPI_Reduce(&value,node_data,1,MPI_DOUBLE,MPI_SUM,0,mpi_comm_shared);

    double my_min_max[2] = { value, -value };
    MPI_Reduce(my_min_max,node_data+1,2,MPI_DOUBLE,MPI_MIN,0,mpi_comm_shared);

    if (mpi_rank_shared == 0) {
      node_data[0] /= double(mpi_size_shared);
      double * root_data = NULL;
      if (mpi_rank == 0) root_data = new double[mpi_size_internode*3];
      MPI_Gather(node_data,3,MPI_DOUBLE,root_data,3,MPI_DOUBLE,0,mpi_comm_internode);
      if (mpi_rank == 0) {
	cout << " > dumpNodeAvgReport: " << message << " (node min avg max)" << endl;
	for (int i = 0; i < mpi_size_internode; ++i)
	  cout << " >> " << i << " " << root_data[3*i+1] << " " << root_data[3*i ] << " " << -root_data[3*i+2] << endl;
	delete[] root_data;
      }
    }

  }

  void dumpNodeSumReport(double value,const string& message) {

    double node_data;
    MPI_Reduce(&value,&node_data,1,MPI_DOUBLE,MPI_SUM,0,mpi_comm_shared);

    if (mpi_rank_shared == 0) {
      double * root_data = NULL;
      if (mpi_rank == 0) root_data = new double[mpi_size_internode];
      MPI_Gather(&node_data,1,MPI_DOUBLE,root_data,1,MPI_DOUBLE,0,mpi_comm_internode);
      if (mpi_rank == 0) {
	cout << " > dumpNodeSumReport: " << message << " (node sum)" << endl;
	for (int i = 0; i < mpi_size_internode; ++i)
	  cout << " >> " << i << " " << root_data[i] << endl;
	delete[] root_data;
      }
    }

  }

  void reorder_csr(int * x_i, int * x_v, int * order, const int n) {

    // backup the old data...
    int * x_i_old = new int[n+1];
    for (int i = 0; i <= n; i++)
      x_i_old[i] = x_i[i];

    int * x_v_old = new int[x_i[n]];
    for (int j = 0; j < x_i[n]; j++)
      x_v_old[j] = x_v[j];

    // put the counts into x_i...
    for (int i = 0; i < n; i++) {
      int i_new = order[i];
      x_i[i_new+1] = x_i_old[i+1] - x_i_old[i]; // put the count
    }

    // build new CSR...
    x_i[0] = 0;
    for (int i = 0; i < n; i++)
      x_i[i+1] += x_i[i];

    // copy in the values...
    for (int i = 0; i < n; i++) {
      int i_new = order[i];
      int j_new_f = x_i[i_new];
      int j_f = x_i_old[i];
      int j_l = x_i_old[i+1]-1;
      for (int j = j_f; j <= j_l; j++) {
        x_v[j-j_f+j_new_f] = x_v_old[j];
      }
    }

    delete[] x_i_old;
    delete[] x_v_old;
  }

  void reorder_csr(int * x_i, int * x_v1, int * x_v2, int * order, const int n) {

    // backup the old data...
    int * x_i_old = new int[n+1];
    for (int i = 0; i <= n; i++)
      x_i_old[i] = x_i[i];

    int * x_v1_old = new int[x_i[n]];
    int * x_v2_old = new int[x_i[n]];
    for (int j = 0; j < x_i[n]; j++) {
      x_v1_old[j] = x_v1[j];
      x_v2_old[j] = x_v2[j];
    }

    // put the counts into x_i...
    for (int i = 0; i < n; i++) {
      int i_new = order[i];
      x_i[i_new+1] = x_i_old[i+1] - x_i_old[i]; // put the count
    }

    // build new CSR...
    x_i[0] = 0;
    for (int i = 0; i < n; i++)
      x_i[i+1] += x_i[i];

    // copy in the values...
    for (int i = 0; i < n; i++) {
      int i_new = order[i];
      int j_new_f = x_i[i_new];
      int j_f = x_i_old[i];
      int j_l = x_i_old[i+1]-1;
      for (int j = j_f; j <= j_l; j++) {
        x_v1[j-j_f+j_new_f] = x_v1_old[j];
        x_v2[j-j_f+j_new_f] = x_v2_old[j];
      }
    }

    delete[] x_i_old;
    delete[] x_v1_old;
    delete[] x_v2_old;

  }

  void reorder(double * var,const int * order,const int n) {

    double * var_copy = new double[n];
    for (int i = 0; i < n; ++i)
      var_copy[i] = var[i];
    for (int i = 0; i < n; ++i) {
      int i_new = order[i];
      assert( (i_new >= 0)&&(i_new < n) );
      var[i_new] = var_copy[i];
    }
    delete[] var_copy;

  }

  void reorder(double (*var)[3],const int * order,const int n) {

    double * var_copy = new double[n];
    for (int j = 0; j < 3; ++j) {
      for (int i = 0; i < n; ++i)
	var_copy[i] = var[i][j];
      for (int i = 0; i < n; ++i) {
	int i_new = order[i];
	assert( (i_new >= 0)&&(i_new < n) );
	var[i_new][j] = var_copy[i];
      }
    }
    delete[] var_copy;

  }

  void reorder(double (*var)[3][3],const int * order,const int n) {

    double * var_copy = new double[n];
    for (int k = 0; k < 3; ++k) {
      for (int j = 0; j < 3; ++j) {
	for (int i = 0; i < n; ++i)
	  var_copy[i] = var[i][j][k];
	for (int i = 0; i < n; ++i) {
	  int i_new = order[i];
	  assert( (i_new >= 0)&&(i_new < n) );
	  var[i_new][j][k] = var_copy[i];
	}
      }
    }
    delete[] var_copy;

  }

  void reorder(int * v,const int * order,const int n) {

    int * v_copy = new int[n];
    for (int i = 0; i < n; ++i)
      v_copy[i] = v[i];
    for (int i = 0; i < n; ++i) {
      int i_new = order[i];
      assert( (i_new >= 0)&&(i_new < n) );
      v[i_new] = v_copy[i];
    }
    delete[] v_copy;

  }

  void reorder(int8 * v,const int * order,const int n) {

    int8 * v_copy = new int8[n];
    for (int i = 0; i < n; ++i)
      v_copy[i] = v[i];
    for (int i = 0; i < n; ++i) {
      int i_new = order[i];
      assert( (i_new >= 0)&&(i_new < n) );
      v[i_new] = v_copy[i];
    }
    delete[] v_copy;

  }

  void reorder(int (*v2)[2],const int * order,const int n) {

    int * v_copy = new int[n];
    for (int j = 0; j < 2; ++j) {
      for (int i = 0; i < n; ++i)
	v_copy[i] = v2[i][j];
      for (int i = 0; i < n; ++i) {
	int i_new = order[i];
	assert( (i_new >= 0)&&(i_new < n) );
	v2[i_new][j] = v_copy[i];
      }
    }
    delete[] v_copy;

  }

  int checkHash(const string &str)
  {
    size_t i = 0;
    if (str.length() == 0)     return 1;

    // ignore whitespace and tab
    while (str[i] == ' ' || str[i] == '\t')
      i++;

    if (str.at(i) == '#')      return 1;
    {
      int hash_position = (int)str.find('#');

      if (hash_position == -1)    return str.length();    //return length of string to be used for the parameter
      else                        return hash_position;   //only string until hash in line appears will be used
    }
  }

  bool parseFunction(string& function,string& args,const string& str) {

    // a function is written like:
    // MAG(U)

    //cout << "str: \"" << str << "\", length=" << str.length() << endl;

    string::size_type leftBracket = str.find_first_of('(',0);
    if (leftBracket == string::npos)
      return(false);

    int level = 1;
    string::size_type nextBracket = leftBracket;
    while (1) {
      nextBracket = str.find_first_of("()",nextBracket+1);
      if (nextBracket == string::npos) {
	cerr << "bracket mismatch in function: \"" << str << "\"" << endl;
	throw(0);
      }
      if (str[nextBracket] == '(') {
	level += 1;
      }
      else {
	// this closes the current level...
	assert(str[nextBracket] == ')');
	level -= 1;
      }
      if (level == 0) {
	if (nextBracket == str.length()-1) {
	  cout << "done" << endl;
	  break;
	}
	else {
	  cerr << "Error: cannot parse complex function yet: " << str << endl;
	  throw(0);
	}
      }
    }

    function = str.substr(0,leftBracket);
    args = str.substr(leftBracket+1,nextBracket-leftBracket-1);

    /*
      cout << "function=\"" << function << "\"" << endl;
      cout << "args=\"" << args << "\"" << endl;
      getchar();
    */

    return(true);

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

  void splitCsv(vector<string>& nameVec,const string& names) {
    // splits a set of comma-delimited names into a vector of name strings
    istringstream iss(names);
    string subtoken;
    while (std::getline(iss,subtoken,','))
      nameVec.push_back(subtoken);
  }

  void splitCsv(set<string>& nameSet,const string& names) {
    // splits a set of comma-delimited names into a set of name strings
    istringstream iss(names);
    string subtoken;
    while (std::getline(iss,subtoken,','))
      nameSet.insert(subtoken);
  }

  bool splitCsv(vector<uint2>& uint2Vec,const string& csv) {
    // splits a set of comma-delimited csv into a vector of uint2s
    istringstream iss(csv);
    string subtoken;
    while (std::getline(iss,subtoken,',')) {
      uint2 value;
      if (from_string<uint2>(value,subtoken,std::dec)) {
        uint2Vec.push_back(value);
      }
      else {
        return false;
      }
    }
    return true;
  }

  bool splitCsv(vector<int>& intVec,const string& csv) {
    // splits a set of comma-delimited csv into a vector of ints
    istringstream iss(csv);
    string subtoken;
    while (std::getline(iss,subtoken,',')) {
      int value;
      if (from_string<int>(value,subtoken,std::dec)) {
        intVec.push_back(value);
      }
      else {
        return false;
      }
    }
    return true;
  }

  bool splitCsv(set<int>& intSet,const string& names) {
    // splits a set of comma-delimited csv into a set of ints
    istringstream iss(names);
    string subtoken;
    while (std::getline(iss,subtoken,',')) {
      int value;
      if (from_string<int>(value,subtoken,std::dec)) {
        intSet.insert(value);
      }
      else {
        return false;
      }
    }
    return true;
  }

  bool splitCsv(vector<double>& doubleVec,const string& csv) {
    // splits a set of comma-delimited csv into a vector of doubles
    istringstream iss(csv);
    string subtoken;
    while (std::getline(iss,subtoken,',')) {
      double value;
      if (from_string<double>(value,subtoken,std::dec)) {
        doubleVec.push_back(value);
      }
      else {
        return false;
      }
    }
    return true;
  }

  bool splitCsv(set<double>& doubleSet,const string& names) {
    // splits a set of comma-delimited csv into a set of doubles
    istringstream iss(names);
    string subtoken;
    while (std::getline(iss,subtoken,',')) {
      double value;
      if (from_string<double>(value,subtoken,std::dec)) {
        doubleSet.insert(value);
      }
      else {
        return false;
      }
    }
    return true;
  }

  string getFileNameExtension(const string &str, const string &delimiters) {
    vector<string> tokens;
    tokenizeString(tokens, str, delimiters);
    return tokens[tokens.size()-1];
  }

  void buildTilePrefix(char * tilePrefix, const char * prefix, const int tile[3]){
    // it is the calling routines responsibility to make sure tilePrefix is long enough...
    std::string str(prefix);
    size_t found = str.find_last_of("//");
    if ((int)found > -1) {
      sprintf(tilePrefix, "%s/%d_%d_%d_%s", str.substr(0, found).c_str(), tile[0], tile[1], tile[2], str.substr(found+1).c_str());
    }
    else {
      sprintf(tilePrefix, "%d_%d_%d_%s", tile[0], tile[1], tile[2], prefix);
    }
  }

  void buildUnindexedFilename(char * filename,const char * prefix,const char * suffix) {
    // if is the calling routines responsibility to make sure filename is long enough...
    sprintf(filename,"%s.%s",prefix,suffix);
  }


  void buildIndexedFilename(char * filename,const char * prefix,const int index,const char * suffix) {
    // if is the calling routines responsibility to make sure filename is long enough...
    sprintf(filename,"%s.%08d.%s",prefix,index,suffix);
  }

  void buildIndexedDirAndFilename(char * filename,const char * dir_prefix,const int dir_index,const int file_index,const char * file_suffix) {
    // if is the calling routines responsibility to make sure filename is long enough...
    sprintf(filename,"%s.%08d/%08d.%s",dir_prefix,dir_index,file_index,file_suffix);
  }

  bool fileExists(const string& filename) {
    return( fileExists(filename.c_str()) );
  }

  bool fileExists(const char * filename) {
    struct stat stFileInfo;
    bool blnReturn;
    int intStat;
    // Attempt to get the file attributes
    intStat = stat(filename,&stFileInfo);
    if(intStat == 0) {
      // We were able to get the file attributes
      // so the file obviously exists.
      blnReturn = true;
    } else {
      // We were not able to get the file attributes.
      // This may mean that we don't have permission to
      // access the folder which contains this file. If you
      // need to do that level of checking, lookup the
      // return values of stat which will give you
      // more details on why stat failed.
      blnReturn = false;
    }
    return(blnReturn);
  }

  int getIndexFromFilename(const string& filename) {

    // here we assume the step index is the last sequence of
    // numbers in the passed filename separated by special
    // chars such as "_./". If no such set of chars can be found,
    // this routine returns -1...

    int index = -1;

    string::size_type pos2 = 0;
    while (1) {

      string::size_type pos1 = filename.find_first_of("_./",pos2);
      if (pos1 == string::npos)
	break;
      pos2 = filename.find_first_not_of("0123456789",pos1+1);
      if (pos2 == string::npos)
	break;

      // we may have a numerical range...
      if (pos2 > pos1+1) {
	index = 0;
	for (int pos = pos1+1,pos_max=pos2; pos < pos_max; ++pos) {
	  index = index*10 + int(filename.at(pos)-'0');
	}
      }

    }

    return index; // could be -1, or the last range of numbers in the name

  }

  int openFile(FILE ** fp,const string& filename,const string& mode) {
    const bool exists = fileExists(filename);
    if (!exists) {
      if (mpi_rank == 0) cout << "\nWARNING: could not find file: " << filename << " at specified location\n" << endl;
      return -1;
    }
    else {
      *fp = fopen(filename.c_str(),mode.c_str());
      if (fp == NULL) { // ME why is this not dereferenced?
        if (mpi_rank == 0) cout << "\nWARNING: problem opening file: " << filename << "\n" << endl;
        fclose(*fp);
        return -2;
      }
    }

    return 0;  // success
  }

  int read2DAsciiTable(double (* &xy_vals)[2], const string filename,const int skip) {
    // read file
    ifstream ifp(filename.c_str());
    if (ifp.fail()) {
      CWARN("was unable to open file " << filename <<"; check path or access permisions");
      return 0;
    }
    string line;
    if (mpi_rank == 0) cout << " > reading 2D table from file \"" << filename << "\"" << endl;
    vector<pair<double,double> > table_vals;
    uint line_count = 0;
    uint valid_count = 0;
    while(!ifp.eof()) {
      getline(ifp,line);
      ++line_count;
      // cout << "line: " << line << endl;
      if (line.empty()) continue;

      char * token = strtok((char *)(line.c_str()),"\t ,");
      if (token[0] == '#') {
        continue;  // skip comment lines
      }
      else {
        int ncols = 0;
        double col1,col2 = 0.0;
        while (token != NULL && ncols<2) {
          if (ncols == 0) col1 = atof(token);
          else if (ncols == 1) col2 = atof(token);
          ++ncols;

          token = strtok(NULL,"\t ,");
        }
        table_vals.push_back(pair<double,double>(col1,col2));
        if ((ncols > 2) && (mpi_rank == 0)) {
          cout << "Warning: row " << (line_count-1) << " had too many elements; only the first two were read" << endl;
        }

        if (ncols < 2) {
          if (mpi_rank == 0) cout << "Warning: row " << (line_count-1) << " had too few elements; ignoring this line" << endl;
          table_vals.pop_back();
        }
        else {
          ++valid_count;
        }
      }
    }

    if (table_vals.size()) {
      if (mpi_rank == 0) cout << " > found " << table_vals.size() << " valid rows in table" << endl;
      const int table_size = int(table_vals.size())/skip;
      xy_vals = new double[table_size][2];
      for (int i=0,limit=table_vals.size(); i<limit; i+=skip) {
        xy_vals[i][0] = table_vals[i].first;
        xy_vals[i][1] = table_vals[i].second;
      }
    }

    return table_vals.size();
  }

  int read3DAsciiTable(double (* &vals)[3], const string filename,const int skip) {
    // read file
    ifstream ifp(filename.c_str());
    if (ifp.fail()) {
      CWARN("was unable to open file " << filename <<"; check path or access permisions");
      return 0;
    }
    string line;
    if (mpi_rank == 0) cout << " > reading 3D table from file \"" << filename << "\"" << endl;
    vector<Triple<double,double,double> > table_vals;
    uint line_count = 0;
    uint valid_count = 0;
    while(!ifp.eof()) {
      ++line_count;
      getline(ifp,line);
      // cout << "line: " << line << endl;
      if (line.empty()) continue;

      char * token = strtok((char *)(line.c_str()),"\t ,");
      if (token[0] == '#') {
        continue;  // skip comment lines
      }
      else {
        int ncols = 0;
        double col1,col2,col3 = 0.0;
        while (token != NULL && ncols<3) {
          if (ncols == 0) col1 = atof(token);
          else if (ncols == 1) col2 = atof(token);
          else if (ncols == 2) col3 = atof(token);
          ++ncols;

          token = strtok(NULL,"\t ,");
        }
        table_vals.push_back(Triple<double,double,double> (col1,col2,col3));
        if ((ncols > 3) && (mpi_rank == 0)) {
          cout << "Warning: row " << (line_count-1) << " had too many elements; only the first three were read" << endl;
        }

        if (ncols < 3) {
          if (mpi_rank == 0) cout << "Warning: row " << (line_count-1) << " had too few elements; ignoring this line" << endl;
          table_vals.pop_back();
        }
        else {
          ++valid_count;
        }
      }
    }

    if (table_vals.size()) {
      if (mpi_rank == 0) cout << " > found " << table_vals.size() << " valid rows in table" << endl;
      const int table_size = int(table_vals.size())/skip;
      vals = new double[table_size][3];
      for (int i=0,limit=table_vals.size(); i<limit; i+=skip) {
        vals[i][0] = table_vals[i].first;
        vals[i][1] = table_vals[i].second;
        vals[i][2] = table_vals[i].third;
      }
    }

    return table_vals.size();
  }

  void resize(bool * &scalar,const int n_new) {
    if (scalar != NULL) {
      delete[] scalar;
      scalar = NULL;
    }
    resize(scalar,0,n_new);
  }

  void resize(bool * &scalar,const int n_old,const int n_new) {
    if (n_old == 0) {
      //assert(scalar == NULL);
      scalar = new bool[n_new];
    }
    else {
      assert(n_old > 0);
      assert(scalar != NULL);
      bool * tmp = new bool[n_new];
      for (int i = 0; i < n_old; ++i)
	tmp[i] = scalar[i];
      delete[] scalar;
      scalar = tmp;
    }
  }

  void resize(int * &scalar,const int n_new) {
    if (scalar != NULL) {
      delete[] scalar;
      scalar = NULL;
    }
    resize(scalar,0,n_new);
  }

  void resize(int * &scalar,const int n_old,const int n_new) {
    if (n_old == 0) {
      //assert(scalar == NULL);
      scalar = new int[n_new];
    }
    else {
      assert(n_old > 0);
      assert(scalar != NULL);
      int * tmp = new int[n_new];
      for (int i = 0; i < n_old; ++i)
	tmp[i] = scalar[i];
      delete[] scalar;
      scalar = tmp;
    }
  }

  void resize(int8 * &scalar,const int n_new) {
    if (scalar != NULL) {
      delete[] scalar;
      scalar = NULL;
    }
    resize(scalar,0,n_new);
  }

  void resize(int8 * &scalar,const int n_old,const int n_new) {
    if (n_old == 0) {
      //assert(scalar == NULL);
      scalar = new int8[n_new];
    }
    else {
      assert(n_old > 0);
      assert(scalar != NULL);
      int8 * tmp = new int8[n_new];
      for (int i = 0; i < n_old; ++i)
	tmp[i] = scalar[i];
      delete[] scalar;
      scalar = tmp;
    }
  }

  void resize(uint8 * &scalar,const int n_new) {
    if (scalar != NULL) {
      delete[] scalar;
      scalar = NULL;
    }
    resize(scalar,0,n_new);
  }

  void resize(uint8 * &scalar,const int n_old,const int n_new) {
    if (n_old == 0) {
      //assert(scalar == NULL);
      scalar = new uint8[n_new];
    }
    else {
      assert(n_old > 0);
      assert(scalar != NULL);
      uint8 * tmp = new uint8[n_new];
      for (int i = 0; i < n_old; ++i)
	tmp[i] = scalar[i];
      delete[] scalar;
      scalar = tmp;
    }
  }

  void resize(int (* &scalar)[2],const int n_new) {
    resize(scalar,0,n_new);
  }

  void resize(int (* &scalar)[2],const int n_old,const int n_new) {
    if (n_old == 0) {
      //assert(scalar == NULL);
      scalar = new int[n_new][2];
    }
    else {
      assert(n_old > 0);
      assert(scalar != NULL);
      int (*tmp)[2] = new int[n_new][2];
      for (int i = 0; i < n_old; ++i) {
	tmp[i][0] = scalar[i][0];
	tmp[i][1] = scalar[i][1];
      }
      delete[] scalar;
      scalar = tmp;
    }
  }

  void resize(int (* &scalar)[4],const int n_new) {
    if (scalar != NULL)
      delete[] scalar;
    scalar = new int[n_new][4];
  }

  void resize(int (* &scalar)[8],const int n_new) {
    if (scalar != NULL)
      delete[] scalar;
    scalar = new int[n_new][8];
  }

  void resize(int (* &scalar)[4],const int n_old,const int n_new) {
    if (n_old == 0) {
      //assert(scalar == NULL);
      scalar = new int[n_new][4];
    }
    else {
      assert(n_old > 0);
      assert(scalar != NULL);
      int (*tmp)[4] = new int[n_new][4];
      for (int i = 0; i < n_old; ++i) {
	for (int j = 0; j < 4; ++j) {
	  tmp[i][j] = scalar[i][j];
	}
      }
      delete[] scalar;
      scalar = tmp;
    }
  }

  void resize(int (* &scalar)[4][2],const int n_old,const int n_new) {
    if (n_old == 0) {
      //assert(scalar == NULL);
      scalar = new int[n_new][4][2];
    }
    else {
      assert(n_old > 0);
      assert(scalar != NULL);
      int (*tmp)[4][2] = new int[n_new][4][2];
      for (int i = 0; i < n_old; ++i) {
	for (int j = 0; j < 4; ++j) {
	  for (int k = 0; k < 2; ++k) {
	    tmp[i][j][k] = scalar[i][j][k];
	  }
	}
      }
      delete[] scalar;
      scalar = tmp;
    }
  }

  void resize(int (* &scalar)[6],const int n_old,const int n_new) {
    if (n_old == 0) {
      //assert(scalar == NULL);
      scalar = new int[n_new][6];
    }
    else {
      assert(n_old > 0);
      assert(scalar != NULL);
      int (*tmp)[6] = new int[n_new][6];
      for (int i = 0; i < n_old; ++i) {
	for (int j = 0; j < 6; ++j) {
	  tmp[i][j] = scalar[i][j];
	}
      }
      delete[] scalar;
      scalar = tmp;
    }
  }

  void resize(int (* &scalar)[6][4],const int n_old,const int n_new) {
    if (n_old == 0) {
      //assert(scalar == NULL);
      scalar = new int[n_new][6][4];
    }
    else {
      assert(n_old > 0);
      assert(scalar != NULL);
      int (*tmp)[6][4] = new int[n_new][6][4];
      for (int i = 0; i < n_old; ++i) {
	for (int j = 0; j < 6; ++j) {
	  for (int k = 0; k < 4; ++k) {
	    tmp[i][j][k] = scalar[i][j][k];
	  }
	}
      }
      delete[] scalar;
      scalar = tmp;
    }
  }

  void resize(int (* &scalar)[8],const int n_old,const int n_new) {
    if (n_old == 0) {
      //assert(scalar == NULL);
      scalar = new int[n_new][8];
    }
    else {
      assert(n_old > 0);
      assert(scalar != NULL);
      int (*tmp)[8] = new int[n_new][8];
      for (int i = 0; i < n_old; ++i) {
	for (int j = 0; j < 8; ++j) {
	  tmp[i][j] = scalar[i][j];
	}
      }
      delete[] scalar;
      scalar = tmp;
    }
  }

  void resize(int (* &scalar)[8][4],const int n_old,const int n_new) {
    if (n_old == 0) {
      //assert(scalar == NULL);
      scalar = new int[n_new][8][4];
    }
    else {
      assert(n_old > 0);
      assert(scalar != NULL);
      int (*tmp)[8][4] = new int[n_new][8][4];
      for (int i = 0; i < n_old; ++i) {
	for (int j = 0; j < 8; ++j) {
	  for (int k = 0; k < 4; ++k) {
	    tmp[i][j][k] = scalar[i][j][k];
	  }
	}
      }
      delete[] scalar;
      scalar = tmp;
    }
  }

  void resize(double * &scalar,const int n_new) {
    if (scalar != NULL) {
      delete[] scalar;
      scalar = NULL;
    }
    resize(scalar,0,n_new);
  }

  void resize(double * &scalar,const int n_old,const int n_new) {
    if (n_old == 0) {
      //assert(scalar == NULL);
      scalar = new double[n_new];
    }
    else {
      assert(n_old > 0);
      assert(scalar != NULL);
      double * tmp = new double[n_new];
      for (int i = 0; i < n_old; ++i)
	tmp[i] = scalar[i];
      delete[] scalar;
      scalar = tmp;
    }
  }

  void resize(double (*&d3)[3],const int n_new) {
    if (d3 != NULL) {
      delete[] d3;
      d3 = NULL;
    }
    resize(d3,0,n_new);
  }

  void resize(double (*&d3)[3],const int n_old,const int n_new) {
    if (n_old == 0) {
      //assert(d3 == NULL);
      d3 = new double[n_new][3];
    }
    else if (n_new > n_old) {
      double (*dtmp)[3] = new double[n_new][3];
      for (int i = 0; i < n_old; ++i) {
	dtmp[i][0] = d3[i][0];
	dtmp[i][1] = d3[i][1];
	dtmp[i][2] = d3[i][2];
      }
      delete[] d3;
      d3 = dtmp;
    }
  }

  void resize(double (*&d3x3)[3][3],const int n_new) {
    if (d3x3 != NULL) {
      delete[] d3x3;
      d3x3 = NULL;
    }
    resize(d3x3,0,n_new);
  }

  void resize(double (*&d3x3)[3][3],const int n_old,const int n_new) {
    if (n_old == 0) {
      //assert(d3x3 == NULL);
      d3x3 = new double[n_new][3][3];
    }
    else if (n_new > n_old) {
      double (*dtmp)[3][3] = new double[n_new][3][3];
      for (int i = 0; i < n_old; ++i)
	for (int j = 0; j < 3; ++j)
	  for (int k = 0; k < 3; ++k)
	    dtmp[i][j][k] = d3x3[i][j][k];
      delete[] d3x3;
      d3x3 = dtmp;
    }
  }

  double uniformRand(const double rmin,const double rmax) {

    return( rmin + (double)rand()/(double)RAND_MAX*(rmax-rmin) );

  }

  // box-muller normally distributed random numbers
  double randn() {

    double u, v,s ;
    u = 0.0; v= 0.0;
    int done = 0 ;
    while ( done == 0 ) {
      u = 2.0 * double( rand()) / double(RAND_MAX) - 1.0;
      v = 2.0 * double( rand()) / double(RAND_MAX) - 1.0;
      s = u*u + v*v ;
      if ( s > 0.0 && s < 1.0 )
        done = 1;
    }

    return u * sqrt( -2.0 * log(s)/ s) ;
  }

  int rangeMin(const int * values,const int n) {

    if (n <= 0)
      return(-(1<<30));

    int min_value = values[0];
    for (int i = 1; i < n; ++i)
      min_value = min( min_value, values[i] );

    return(min_value);

  }

  int rangeMax(const int * values,const int n) {

    if (n <= 0)
      return((1<<30));

    int max_value = values[0];
    for (int i = 1; i < n; ++i)
      max_value = max( max_value, values[i] );

    return(max_value);

  }

  void removeMean(double * var,const double * norm,const int n) {

    double my_buf[2] = { 0.0, 0.0 };
    for (int i = 0; i < n; ++i) {
      my_buf[0] += var[i]*norm[i];
      my_buf[1] += norm[i];
    }

    double buf[2];
    MPI_Allreduce(my_buf,buf,2,MPI_DOUBLE,MPI_SUM,mpi_comm);
    buf[0] /= buf[1];

    for (int i = 0; i < n; ++i)
      var[i] -= buf[0];


  }

  double calcMean(double * var,const double * norm,const int n) {

    double my_buf[2] = { 0.0, 0.0 };
    for (int i = 0; i < n; ++i) {
      my_buf[0] += var[i]*norm[i];
      my_buf[1] += norm[i];
    }

    double buf[2];
    MPI_Allreduce(my_buf,buf,2,MPI_DOUBLE,MPI_SUM,mpi_comm);

    return(buf[0]/buf[1]);

  }

  void newArray2d(double ** &data,const int ni,const int nj) {
    assert( ni >= 0 );
    assert( nj >= 0 );
    assert( data == NULL );
    if ((ni > 0)&&(nj > 0)) {
      data = new double*[ni];
      data[0] = new double[ni*nj]; // defines data[0][0]
      for (int i = 1; i < ni; ++i)
	data[i] = data[0] + nj*i;
    }
  }

  double ** newArray2d(const int ni,const int nj) {
    double ** data = NULL;
    newArray2d(data,ni,nj);
    return data;
  }

  void deleteArray2d(double ** &data) {
    if (data != NULL) {
      delete[] data[0];
      delete[] data;
      data = NULL;
    }
  }

  // so here is a pretty interesting use of const: check out the
  // second arg...

  void transposeArray2d(double ** pt,const double * const * p,const int * iora,const int * jora) {

    // data in p is distributed in iora[] and full jora[mpi_size]...

    const int ni_local = iora[mpi_rank+1]-iora[mpi_rank];
    const int ni = iora[mpi_size];

    const int nj_local = jora[mpi_rank+1]-jora[mpi_rank];
    const int nj = jora[mpi_size];

    // pack the send-side stuff...

    double * send_buffer = new double[ni_local*nj];
    int * send_count = new int[mpi_size];
    int ii = 0;
    for (int rank = 0; rank < mpi_size; ++rank) {
      send_count[rank] = ni_local*(jora[rank+1]-jora[rank]);
      for (int j = jora[rank]; j != jora[rank+1]; ++j)
	for (int i = 0; i < ni_local; ++i)
	  send_buffer[ii++] = p[i][j];
    }
    assert( ii == ni_local*nj );

    int * send_disp = new int[mpi_size];
    send_disp[0] = 0;
    for (int rank = 1; rank < mpi_size; ++rank)
      send_disp[rank] = send_disp[rank-1] + send_count[rank-1];

    // the recv-side stuff...

    double * recv_buffer = new double[nj_local*ni];
    int * recv_count = new int[mpi_size];
    for (int rank = 0; rank < mpi_size; ++rank)
      recv_count[rank] = nj_local*(iora[rank+1]-iora[rank]);

    int * recv_disp = new int[mpi_size];
    recv_disp[0] = 0;
    for (int rank = 1; rank < mpi_size; ++rank)
      recv_disp[rank] = recv_disp[rank-1] + recv_count[rank-1];

    // exchange...

    MPI_Alltoallv(send_buffer,send_count,send_disp,MPI_DOUBLE,
		  recv_buffer,recv_count,recv_disp,MPI_DOUBLE,mpi_comm);

    // clear send stuff...

    delete[] send_buffer;
    delete[] send_count;
    delete[] send_disp;
    delete[] recv_count;
    delete[] recv_disp;

    // now unpack into the transpose...

    int jj = 0;
    for (int rank = 0; rank < mpi_size; ++rank)
      for (int j = 0; j < nj_local; ++j)
	for (int i = iora[rank]; i != iora[rank+1]; ++i)
	  pt[j][i] = recv_buffer[jj++];
    assert( jj == nj_local*ni );

    // cleanup...

    delete[] recv_buffer;

  }

  void newComplexArray2d(complex<double> ** &data,const int ni,const int nj) {
    assert( ni >= 0 );
    assert( nj >= 0 );
    assert( data == NULL );
    if ((ni > 0)&&(nj > 0)) {
      data = new complex<double>*[ni];
      data[0] = new complex<double>[ni*nj]; // defines data[0][0]
      for (int i = 1; i < ni; ++i)
	data[i] = data[0] + nj*i;
    }
  }

  complex<double> ** newComplexArray2d(const int ni,const int nj) {
    complex<double> ** data = NULL;
    newComplexArray2d(data,ni,nj);
    return data;
  }

  void deleteComplexArray2d(complex<double> ** &data) {
    if (data != NULL) {
      delete[] data[0];
      delete[] data;
      data = NULL;
    }
  }

  void mkdir_for_file(const string& filename) {

    // should only call this from one rank-- not necessarily 0 though...
    //assert(mpi_rank == 0);

    string dir;
    int err = 0;
    size_t pos_ = 0;
    for (int i=0; i < std::count(filename.begin(),filename.end(),'/'); ++i) {
      size_t pos = filename.find('/',pos_);
      pos_ = pos+1;
      dir = filename.substr(0,pos);
      if (pos != 0) err = mkdir(dir.c_str(),0777);
    }
    if (err != 0 && errno != EEXIST)
      cout << "Warning: unable to create directory: " << dir << endl;

  }

  void mkdir_for_file_collective(const string& filename, const int write_rank) {

    // collective call: all ranks must call. default behavior is write_rank == 0.

    if (mpi_rank == write_rank) {

      mkdir_for_file(filename);

    }

    // this has to be synchronous, so no one tries to open the file before it's directory is created...

    MPI_Barrier(mpi_comm);

  }

  string makeTmpPrefix(const string& prefix) {
    // adds a "." to the passed filename for use a temporary file...
    size_t last_slash = prefix.find_last_of('/');
    if (last_slash != string::npos) {
      return prefix.substr(0,last_slash+1) + "." + prefix.substr(last_slash+1);
    }
    else {
      return "." + prefix;
    }
  }

  string getFilenameNoPath(const string& filename_plus_path) {
    size_t last_slash = filename_plus_path.find_last_of('/');
    if (last_slash != string::npos) {
      return filename_plus_path.substr(last_slash+1);
    }
    else {
      return filename_plus_path;
    }
  }

  string addPrefixToFilename(const string& prefix,const string& filename) {
    // adds "prefix" to the passed filename...
    // only challenge here is the additional directory names in the path that we do not want to mess up...
    const size_t last_slash = filename.find_last_of('/');
    if (last_slash != string::npos) {
      return filename.substr(0,last_slash+1) + prefix + filename.substr(last_slash+1);
    }
    else {
      return prefix + filename;
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

  // Erase all Occurrences of given substring from main string.
  void eraseAllSubStr(string& mainStr, const string& toErase) {
    size_t pos = string::npos;

    // Search for the substring in string in a loop untill nothing is found
    while ((pos  = mainStr.find(toErase) )!= string::npos) {
      // If found then erase it from string
      mainStr.erase(pos, toErase.length());
    }
  }

  // void replaceOldSubStringWithNewSubString(string& str,const string& oldSubString,const string& newSubString) {
  //   int iter = 0;
  //   size_t np;
  //   while ((np = str.find(oldSubString)) != string::npos) {
  //     assert(++iter < 10);
  //     string before = str.substr(0,np);
  //     string after  = str.substr(np+oldSubString.size(),string::npos);
  //     if (after.size() == 0) {
  //       if (before.size() == 0) {
  //         str = newSubString;
  //       }
  //       else {
  //         str = before+newSubString;
  //       }
  //       break;
  //     }
  //     else if (before.size() == 0) {
  //       str = newSubString+after;
  //     }
  //     else {
  //       str = before+newSubString+after;
  //     }
  //     np += newSubString.size() - oldSubString.size();
  //   }
  // }

  int gaussj(double a[][3], int n, double b[][3], int m) {

    //  Object : Linear equation solution by Gauss-Jordan elimination.
    //           a(1:n,1:n) is an input matrix stored in an array of physical
    //           dimensions np by np. b(1:n,1:m) is an input matrix containing
    //           the m right-hand side vectors, stored in an array of physical
    //           dimensions np by nm. On output, a(1:n,1:n) is replaced by its
    //           matrix inverse, and b(1:n,1:m) is replaced by the corresponding
    //           set of solution vectors.
    //           Parameter: NMAX is the largest anticipated value of n.
    //
    // Numerical Recipes (http://www.nr.com/)
    //                   (http://www.library.cornell.edu/nr/cbookfpdf.html)

    const int nmax = 3;
    assert(n <= nmax);
    assert(m <= 3);
    int indxc[nmax], indxr[nmax], ipiv[nmax];
    double big,dum,pivinv;

    int icol=0;
    int irow=0;
    for(int j=0; j<n; ++j) ipiv[j]=0;

    for(int i=0; i<n; ++i){
      big=0.0;
      for(int j=0; j<n; ++j){
	if(ipiv[j] != 1){
	  for(int k=0; k<n; ++k){
	    if(ipiv[k] == 0){
	      if(fabs(a[j][k]) >= big){
		big = fabs(a[j][k]);
		irow = j;
		icol = k;
	      }
	    }
	  }
	}
      }
      ipiv[icol] += 1;
      if(irow != icol){
	for(int l=0; l<n; ++l){
	  dum = a[irow][l];
	  a[irow][l] = a[icol][l];
	  a[icol][l] = dum;
	}
	for(int l=0; l<m; ++l){
	  dum = b[irow][l];
	  b[irow][l] = b[icol][l];
	  b[icol][l] = dum;
	}
      }
      indxr[i] = irow;
      indxc[i] = icol;
      if(a[icol][icol] == 0.0){
	cout << "Singular Matrix in gaussj" << endl;
	return(1);
      }
      pivinv = 1.0/a[icol][icol];
      a[icol][icol] = 1.0;
      for(int l=0; l<n; ++l) a[icol][l] *= pivinv;
      for(int l=0; l<m; ++l) b[icol][l] *= pivinv;

      for(int ll=0; ll<n; ++ll){
	if(ll != icol){
	  dum = a[ll][icol];
	  a[ll][icol] = 0.0;
	  for(int l=0; l<n; ++l) a[ll][l] -= a[icol][l]*dum;
	  for(int l=0; l<m; ++l) b[ll][l] -= b[icol][l]*dum;
	}
      }
    }

    for(int l=n-1; l>=0; --l){
      if(indxr[l] != indxc[l]){
	for(int k=0; k<n; ++k){
	  dum = a[k][indxr[l]];
	  a[k][indxr[l]] = a[k][indxc[l]];
	  a[k][indxc[l]] = dum;
	}
      }
    }

    return(0);
  }

  bool isNan(const double& value) {
    return(value != value);
  }

  bool compareFirstSecondIntIntPair(const std::pair<int,int>& a,const std::pair<int,int>& b) {
    return( (a.first < b.first) || ((a.first == b.first)&&(a.second < b.second)) );
  }

  bool compareSecondIntDoublePair(const std::pair<int,double>& a,const std::pair<int,double>& b) {
    return( a.second < b.second );
  }

  bool startsWith(const string& bigString,const string& smallString) {
    return( bigString.compare(0, smallString.length(), smallString) == 0 );
  }

  bool strcmp_wildcard(const string& str,const string& str_wildcard) {
    // returns true if str matches str_wildcard where str_wildcard
    // can contain one "*" character
    const size_t nstart = str_wildcard.find("*"); // ignore
    if (nstart == string::npos) {
      //cout << "no wildcard. Just compare" << endl;
      return str == str_wildcard;
    }
    //cout << "nstart: " << nstart << " str_wildcard.length(): " << str_wildcard.length() << endl;
    if (nstart > 0) {
      // there are some charaters at the front of str_wildcard that
      // we need to compare...
      if (nstart > str.length()) {
	// it is not possible to match...
	//cout << "start of str_wildcard is too long to match!" << endl;
	return false;
      }
      else if (str.compare(0,nstart,str_wildcard,0,nstart) != 0) {
	//cout << "start does not match!" << endl;
	return false;
      }
      //cout << "start matches" << endl;
    }
    // now the back...
    if (nstart+1 < str_wildcard.length()) {
      const size_t nend = str_wildcard.length()-(nstart+1);
      //cout << "some stuff at the back: nend: " << nend << endl;
      if (nstart + nend > str.length()) {
	// it is not possible to match...
	//cout << "end of str_wildcard is too long to match!" << endl;
	return false;
      }
      else if (str.compare(str.length()-nend,nend,str_wildcard,nstart+1,nend) != 0) {
	//cout << "end does not match!" << endl;
	return false;
      }
      //cout << "end matches" << endl;
    }
    return true;
  }

  bool regexMatchesString(const string& regex,const string& str) {

    string::size_type starPos = regex.find_first_of("*",0);

    // note: here we check for % sign instead of *, because * gets replaced
    // in the command line unless the user is careful to use a \*

    if (starPos == string::npos)
      starPos = regex.find_first_of("%",0);

    if (starPos == string::npos)
      return(regex == str);
    else if (starPos == 0)
      return true;
    else
      return(str.compare(0,starPos,regex,0,starPos) == 0);

  }

  // inline and template probably .. XXX
  void quickSortByValue( int * idx, int *vals, int n) {

    if ( n < 2 ) return ;

    int p      = vals[idx[n/2]];
    int * l  = idx ;
    int * r  = idx + n -1 ;

    while ( l <= r) {
      if ( vals[*l] < p ) {
	l++ ;
	continue ;
      }
      if ( vals[*r] > p ) {
	r-- ;
	continue ;
      }


      // else swap
      int t = *l ;
      *l++  = *r ;
      *r--  =  t ;
    }//while l<=r

    // recurse
    quickSortByValue( idx, vals, r-idx + 1) ;
    quickSortByValue( l  , vals, idx+n - l) ;
  } // quickSortByValue ..

  void quickSortByValue( int * idx, double *vals, int n) {

    if ( n < 2 ) return ;

    double p      = vals[idx[n/2]];
    int * l  = idx ;
    int * r  = idx + n -1 ;

    while ( l <= r) {
      if ( vals[*l] < p ) {
	l++ ;
	continue ;
      }
      if ( vals[*r] > p ) {
	r-- ;
	continue ;
      }


      // else swap
      int t = *l ;
      *l++  = *r ;
      *r--  =  t ;
    }//while l<=r

    // recurse
    quickSortByValue( idx, vals, r-idx + 1) ;
    quickSortByValue( l  , vals, idx+n - l) ;
  } // quickSortByValue ..

  // ==============================
  // eigenvalue stuff
  // ==============================
#define n 3

  static double hypot2(double x, double y) {
    return sqrt(x*x+y*y);
  }

  // Symmetric Householder reduction to tridiagonal form.
  static void tred2(double V[n][n], double d[n], double e[n]) {

    //  This is derived from the Algol procedures tred2 by
    //  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
    //  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
    //  Fortran subroutine in EISPACK.

    for (int j = 0; j < n; j++) {
      d[j] = V[n-1][j];
    }

    // Householder reduction to tridiagonal form.

    for (int i = n-1; i > 0; i--) {

      // Scale to avoid under/overflow.

      double scale = 0.0;
      double h = 0.0;
      for (int k = 0; k < i; k++) {
	scale = scale + fabs(d[k]);
      }
      if (scale == 0.0) {
	e[i] = d[i-1];
	for (int j = 0; j < i; j++) {
	  d[j] = V[i-1][j];
	  V[i][j] = 0.0;
	  V[j][i] = 0.0;
	}
      } else {

	// Generate Householder vector.

	for (int k = 0; k < i; k++) {
	  d[k] /= scale;
	  h += d[k] * d[k];
	}
	double f = d[i-1];
	double g = sqrt(h);
	if (f > 0) {
	  g = -g;
	}
	e[i] = scale * g;
	h = h - f * g;
	d[i-1] = f - g;
	for (int j = 0; j < i; j++) {
	  e[j] = 0.0;
	}

	// Apply similarity transformation to remaining columns.

	for (int j = 0; j < i; j++) {
	  f = d[j];
	  V[j][i] = f;
	  g = e[j] + V[j][j] * f;
	  for (int k = j+1; k <= i-1; k++) {
	    g += V[k][j] * d[k];
	    e[k] += V[k][j] * f;
	  }
	  e[j] = g;
	}
	f = 0.0;
	for (int j = 0; j < i; j++) {
	  e[j] /= h;
	  f += e[j] * d[j];
	}
	double hh = f / (h + h);
	for (int j = 0; j < i; j++) {
	  e[j] -= hh * d[j];
	}
	for (int j = 0; j < i; j++) {
	  f = d[j];
	  g = e[j];
	  for (int k = j; k <= i-1; k++) {
	    V[k][j] -= (f * e[k] + g * d[k]);
	  }
	  d[j] = V[i-1][j];
	  V[i][j] = 0.0;
	}
      }
      d[i] = h;
    }

    // Accumulate transformations.

    for (int i = 0; i < n-1; i++) {
      V[n-1][i] = V[i][i];
      V[i][i] = 1.0;
      double h = d[i+1];
      if (h != 0.0) {
	for (int k = 0; k <= i; k++) {
	  d[k] = V[k][i+1] / h;
	}
	for (int j = 0; j <= i; j++) {
	  double g = 0.0;
	  for (int k = 0; k <= i; k++) {
	    g += V[k][i+1] * V[k][j];
	  }
	  for (int k = 0; k <= i; k++) {
	    V[k][j] -= g * d[k];
	  }
	}
      }
      for (int k = 0; k <= i; k++) {
	V[k][i+1] = 0.0;
      }
    }
    for (int j = 0; j < n; j++) {
      d[j] = V[n-1][j];
      V[n-1][j] = 0.0;
    }
    V[n-1][n-1] = 1.0;
    e[0] = 0.0;
  }

  // Symmetric tridiagonal QL algorithm.
#define FAST_MAX(a, b) ((a)>(b)?(a):(b))

  static void tql2(double V[n][n], double d[n], double e[n]) {

    //  This is derived from the Algol procedures tql2, by
    //  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
    //  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
    //  Fortran subroutine in EISPACK.

    for (int i = 1; i < n; i++) {
      e[i-1] = e[i];
    }
    e[n-1] = 0.0;

    double f = 0.0;
    double tst1 = 0.0;
    double eps = pow(2.0,-52.0);
    for (int l = 0; l < n; l++) {

      // Find small subdiagonal element

      tst1 = FAST_MAX(tst1,fabs(d[l]) + fabs(e[l]));
      int m = l;
      while (m < n) {
	if (fabs(e[m]) <= eps*tst1) {
	  break;
	}
	m++;
      }

      // If m == l, d[l] is an eigenvalue,
      // otherwise, iterate.

      if (m > l) {
	int iter = 0;
	do {
	  iter = iter + 1;  // (Could check iteration count here.)

	  // Compute implicit shift

	  double g = d[l];
	  double p = (d[l+1] - g) / (2.0 * e[l]);
	  double r = hypot2(p,1.0);
	  if (p < 0) {
	    r = -r;
	  }
	  d[l] = e[l] / (p + r);
	  d[l+1] = e[l] * (p + r);
	  double dl1 = d[l+1];
	  double h = g - d[l];
	  for (int i = l+2; i < n; i++) {
	    d[i] -= h;
	  }
	  f = f + h;

	  // Implicit QL transformation.

	  p = d[m];
	  double c = 1.0;
	  double c2 = c;
	  double c3 = c;
	  double el1 = e[l+1];
	  double s = 0.0;
	  double s2 = 0.0;
	  for (int i = m-1; i >= l; i--) {
	    c3 = c2;
	    c2 = c;
	    s2 = s;
	    g = c * e[i];
	    h = c * p;
	    r = hypot2(p,e[i]);
	    e[i+1] = s * r;
	    s = e[i] / r;
	    c = p / r;
	    p = c * d[i] - s * g;
	    d[i+1] = h + s * (c * g + s * d[i]);

	    // Accumulate transformation.

	    for (int k = 0; k < n; k++) {
	      h = V[k][i+1];
	      V[k][i+1] = s * V[k][i] + c * h;
	      V[k][i] = c * V[k][i] - s * h;
	    }
	  }
	  p = -s * s2 * c3 * el1 * e[l] / dl1;
	  e[l] = s * p;
	  d[l] = c * p;

	  // Check for convergence.

	} while (fabs(e[l]) > eps*tst1);
      }
      d[l] = d[l] + f;
      e[l] = 0.0;
    }

    // Sort eigenvalues and corresponding vectors.

    for (int i = 0; i < n-1; i++) {
      int k = i;
      double p = d[i];
      for (int j = i+1; j < n; j++) {
	if (d[j] < p) {
	  k = j;
	  p = d[j];
	}
      }
      if (k != i) {
	d[k] = d[i];
	d[i] = p;
	for (int j = 0; j < n; j++) {
	  p = V[j][i];
	  V[j][i] = V[j][k];
	  V[j][k] = p;
	}
      }
    }
  }

  void eigenDecomposition(double A[n][n], double V[n][n], double d[n]) {
    double e[n];
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
	V[i][j] = A[i][j];
      }
    }
    tred2(V, d, e);
    tql2(V, d, e);
  }

#undef n

  void transpose(double A[3][3]) {
    {
      const double tmp = A[0][1];
      A[0][1] = A[1][0];
      A[1][0] = tmp;
    }
    {
      const double tmp = A[0][2];
      A[0][2] = A[2][0];
      A[2][0] = tmp;
    }
    {
      const double tmp = A[1][2];
      A[1][2] = A[2][1];
      A[2][1] = tmp;
    }
  }



  void setBit(int &i, int ibit) {
    i |= (1 << ibit);
  }

  void clearBit( int &i, int ibit) {
    i &= ~(1 << ibit) ;
  }

  int setNextBit(int &i) {
    for (int ibit = 0 ; ibit < 32; ++ibit) {
      int mask = 1 << ibit ;
      if  ( (i & mask) == 0 ) {
        i |= mask ;
        return ibit;
      }
    }

    return -1;
  }

  int clearLastSetBit( int& i ) {
    for (int ibit = 31 ; ibit >= 0 ; --ibit ) {
      int mask = 1 << ibit ;
      if ( i & mask) {
        i &= ~mask ;
        return ibit;
      }
    }//ibit

    return -1;
  }

  int nSetBits(const int &i) {
    int n = 0 ;
    for (int ibit = 0 ; ibit < 32 ; ++ibit) {
      if ( i & (1<<ibit)) ++n ;
    }
    return n;
  }

  void getBestE1E2FromE0(double e1[3],double e2[3],const double e0[3]) {

    // returns unit vectors e1,e2. Note that e0 need not be a unit vector.

    // recall
    //CROSS_PRODUCT(A,B) { (A)[1]*(B)[2]-(A)[2]*(B)[1] , (A)[2]*(B)[0]-(A)[0]*(B)[2] , (A)[0]*(B)[1]-(A)[1]*(B)[0] }

    if (fabs(e0[0]) < min(fabs(e0[1]),fabs(e0[2]))) {
      // x-component of e0 is the smallest, e1 is e0 x (1,0,0)...
      const double mag = sqrt(e0[1]*e0[1] + e0[2]*e0[2]);
      assert(mag > 0.0);
      e1[0] = 0.0;
      e1[1] = -e0[2]/mag;
      e1[2] = e0[1]/mag;
      // e2 is e0 x e1...
      e2[0] = e0[1]*e1[2]-e0[2]*e1[1];
      e2[1] = -e0[0]*e1[2];
      e2[2] = e0[0]*e1[1];
    }
    else if (fabs(e0[1]) < fabs(e0[2])) {
      // y component is the smallest: e1 = e0 x (0,1,0)...
      const double mag = sqrt(e0[0]*e0[0] + e0[2]*e0[2]);
      assert(mag > 0.0);
      e1[0] = e0[2]/mag;
      e1[1] = 0.0;
      e1[2] = -e0[0]/mag;
      e2[0] = e0[1]*e1[2];
      e2[1] = e0[2]*e1[0]-e0[0]*e1[2];
      e2[2] = -e0[1]*e1[0];
    }
    else {
      // z component is the smallest: e1 = e0 x (0,0,1)...
      const double mag = sqrt(e0[0]*e0[0] + e0[1]*e0[1]);
      assert(mag > 0.0);
      e1[0] = -e0[1]/mag;
      e1[1] = e0[0]/mag;
      e1[2] = 0.0;
      e2[0] = -e0[2]*e1[1];
      e2[1] = e0[2]*e1[0];
      e2[2] = e0[0]*e1[1]-e0[1]*e1[0];
    }

    // this should only be necessary if e0 was not a unit normal:
    const double mag = MAG(e2);
    assert(mag > 0.0);
    FOR_I3 e2[i] /= mag;

  }

  // ========================================================
  // sometimes these Point/Tri functions are used to cycle through a set of
  // tris and get the closest. It may be possible to do this
  // much faster by eliminating tris that are guaranteed to be further
  // away than the current d2. For now we use the calling process
  // to do this.
  // ========================================================

  double calcAngle(const double x0[3],const double x1[3],const double x2[3]) {
    // returns the angle in radians between 2 connected line segments x0->x1, and x1->x2. if
    // one or both line segments have zero length, this returns 0.0...
    const double v1[3] = DIFF(x1,x0);
    const double mag_v1 = MAG(v1);
    const double v2[3] = DIFF(x2,x1);
    const double mag_v2 = MAG(v2);
    if ((mag_v1 <= 0.0)||(mag_v2 <= 0.0))
      return 0.0;
    return acos(min(1.0,DOT_PRODUCT(v1,v2)/mag_v1/mag_v2));
  }

  double getCosAngleBetweenTris(const double v00[3],const double v01[3],const double v02[3],const double v10[3],const double v11[3],const double v12[3]) {
    // assumes user has checked before passing in v* that magnitude is non-zero
    double n_0[3] = TRI_NORMAL_2(v00,v01,v02);
    NORMALIZE(n_0);
    double n_1[3] = TRI_NORMAL_2(v10,v11,v12);
    NORMALIZE(n_1);
    return DOT_PRODUCT(n_0,n_1);
  }

  double getPointToLineDist(const double xp[3],const double x0[3],const double n0[3]) {
    return sqrt(getPointToLineDist2(xp,x0,n0));
  }

  double getPointToLineDist2(const double xp[3],const double x0[3],const double n0[3]) {
    const double dxp[3] = DIFF(xp,x0);
    const double dp = DOT_PRODUCT(n0,dxp);
    const double dx2 = DOT_PRODUCT(n0,n0);
    return( dxp[0]*dxp[0] + dxp[1]*dxp[1] + dxp[2]*dxp[2] - dp*dp/dx2 );
  }

  double getPointToPlaneDist(const double xp[3],const double x0[3],const double n0[3]) {
    return (xp[0]-x0[0])*n0[0] + (xp[1]-x0[1])*n0[1] + (xp[2]-x0[2])*n0[2];
  }

  double getPointToEdgeDist2(const double xp[3],const double v0[3],const double v1[3]) {

    // returns the square of the euclidean distance from point xp to edge v0->v1.

    const double dx[3] = DIFF(v1,v0);
    const double dxp[3] = DIFF(xp,v0);
    const double dp = DOT_PRODUCT(dx,dxp);
    if (dp <= 0.0) {
      // we are closest to the first point v0. Note this includes the case when dx is zero, and/or dxp is zero...
      return( dxp[0]*dxp[0] + dxp[1]*dxp[1] + dxp[2]*dxp[2] );
    }
    else {
      const double dx2 = DOT_PRODUCT(dx,dx);
      if (dp >= dx2) {
	// we are closest to the second point v1...
	return ( (xp[0]-v1[0])*(xp[0]-v1[0]) +
		 (xp[1]-v1[1])*(xp[1]-v1[1]) +
		 (xp[2]-v1[2])*(xp[2]-v1[2]) );
      }
      else {
        // can get slightly negative due to machine precision errors...
	return max(0.0, dxp[0]*dxp[0] + dxp[1]*dxp[1] + dxp[2]*dxp[2] - dp*dp/dx2);
      }
    }

  }

  void getClosestPointOnEdgeRobust(double xc[3],const double xp[3],const double v0[3],const double v1[3]) {

    // returns the square of the euclidean distance from point xp to edge v0->v1.

    const double dx[3] = DIFF(v1,v0);
    const double dxp[3] = DIFF(xp,v0);
    const double dp = DOT_PRODUCT(dx,dxp);
    const double dx2 = DOT_PRODUCT(dx,dx);
    const double eps = 1.0E-15*dx2;
    if (dp <= eps) {
      // we are closest to the first point v0. Note this includes the case when dx is zero, and/or dxp is zero...
      FOR_I3 xc[i] = v0[i];
    }
    else {
      if (dp >= dx2-eps) {
	// we are closest to the second point v1...
        FOR_I3 xc[i] = v1[i];
      }
      else {
        FOR_I3 xc[i] = dp/dx2*v1[i] + (1.0-dp/dx2)*v0[i];
      }
    }

  }

  void getPointToEdgeDir(double dir[3],const double xp[3],const double v0[3],const double v1[3]) {

    // euclidean dist weighted direction from point xp to edge v0->v1

    const double dx[3] = DIFF(v1,v0);
    const double dxp[3] = DIFF(xp,v0);
    const double dp = DOT_PRODUCT(dx,dxp);
    if (dp <= 0.0) {
      // we are closest to the first point v0. Note this includes the case when dx is zero, and/or dxp is zero...
      FOR_I3 dir[i] = -dxp[i];
    }
    else {
      const double dx2 = DOT_PRODUCT(dx,dx);
      if (dp >= dx2) {
	// we are closest to the second point v1...
        FOR_I3 dir[i] = v1[i]-xp[i];
      }
      else {
        FOR_I3 dir[i] = v0[i]+dp/dx2*dx[i]-xp[i]; // x_proj-xp[i]
      }
    }

  }

  /*
    double getPointToTriDist2(const double xp[3],const double v0[3],const double v1[3],const double v2[3]) {

    double xc[3];
    getClosestPointOnTri(xc,xp,v0,v1,v2);

    //return distance squared
    return ( (xp[0]-xc[0])*(xp[0]-xc[0]) +
    (xp[1]-xc[1])*(xp[1]-xc[1]) +
    (xp[2]-xc[2])*(xp[2]-xc[2]) );

    }
  */

  double getEdgeToEdgeDist2(const double v00[3],const double v01[3], const double v10[3],const double v11[3]) {
    // get minimum distance between edges (v00,v01) and (v10,v11)
    const double dx0[3] = DIFF(v01,v00);
    const double dx1[3] = DIFF(v11,v10);
    const double dx01[3] = DIFF(v10,v00);

    // DENOM = dx0 * dx0 * dy1 * dy1 + dx0 * dx0 * dz1 * dz1 - 2 * dx0 * dx1 * dy0 * dy1 - 2 * dx0 * dx1 * dz0 * dz1 + dx1 * dx1 * dy0 * dy0 + dx1 * dx1 * dz0 * dz0 + dy0 * dy0 * dz1 * dz1 - 2 * dy0 * dy1 * dz0 * dz1 + dy1 * dy1 * dz0 * dz0;cg = dx0 * dx0 * dy1 * dy1 + dx0 * dx0 * dz1 * dz1 - 2 * dx0 * dx1 * dy0 * dy1 - 2 * dx0 * dx1 * dz0 * dz1 + dx1 * dx1 * dy0 * dy0 + dx1 * dx1 * dz0 * dz0 + dy0 * dy0 * dz1 * dz1 - 2 * dy0 * dy1 * dz0 * dz1 + dy1 * dy1 * dz0 * dz0;
    // A = (dx0 * dy1 * dy1 + dx0 * dz1 * dz1 - dx1 * dy0 * dy1 - dx1 * dz0 * dz1) * dx01 + (-dx0 * dx1 * dy1 + dx1 * dx1 * dy0 + dy0 * dz1 * dz1 - dy1 * dz0 * dz1) * dy01 + (-dx0 * dx1 * dz1 + dx1 * dx1 * dz0 - dy0 * dy1 * dz1 + dy1 * dy1 * dz0) * dz01;
    // B = (dx0 * dy0 * dy1 + dx0 * dz0 * dz1 - dx1 * dy0 * dy0 - dx1 * dz0 * dz0) * dx01 + (-dx0 * dx0 * dy1 + dx0 * dx1 * dy0 + dy0 * dz0 * dz1 - dy1 * dz0 * dz0) * dy01 + (-dx0 * dx0 * dz1 + dx0 * dx1 * dz0 - dy0 * dy0 * dz1 + dy0 * dy1 * dz0) * dz01;

    const double denom = dx0[0] * dx0[0] * dx1[1] * dx1[1] + dx0[0] * dx0[0] * dx1[2] * dx1[2] - 2 * dx0[0] * dx1[0] * dx0[1] * dx1[1] - 2 * dx0[0] * dx1[0] * dx0[2] * dx1[2] + dx1[0] * dx1[0] * dx0[1] * dx0[1] + dx1[0] * dx1[0] * dx0[2] * dx0[2] + dx0[1] * dx0[1] * dx1[2] * dx1[2] - 2 * dx0[1] * dx1[1] * dx0[2] * dx1[2] + dx1[1] * dx1[1] * dx0[2] * dx0[2];

    if (fabs(denom) < 1.0e-14) {
      // indicates edges are essentially parallel, so return min node-node distance
      double d2_a = DIST2(v10,v00);
      double d2_b = DIST2(v11,v00);
      return min(d2_a,d2_b);
    }

    const double one_o_denom = 1.0/denom;

    double A = (dx0[0] * dx1[1] * dx1[1] + dx0[0] * dx1[2] * dx1[2] - dx1[0] * dx0[1] * dx1[1] - dx1[0] * dx0[2] * dx1[2]) * dx01[0] + (-dx0[0] * dx1[0] * dx1[1] + dx1[0] * dx1[0] * dx0[1] + dx0[1] * dx1[2] * dx1[2] - dx1[1] * dx0[2] * dx1[2]) * dx01[1] + (-dx0[0] * dx1[0] * dx1[2] + dx1[0] * dx1[0] * dx0[2] - dx0[1] * dx1[1] * dx1[2] + dx1[1] * dx1[1] * dx0[2]) * dx01[2];

    A *= one_o_denom;

    double B = (dx0[0] * dx0[1] * dx1[1] + dx0[0] * dx0[2] * dx1[2] - dx1[0] * dx0[1] * dx0[1] - dx1[0] * dx0[2] * dx0[2]) * dx01[0] + (-dx0[0] * dx0[0] * dx1[1] + dx0[0] * dx1[0] * dx0[1] + dx0[1] * dx0[2] * dx1[2] - dx1[1] * dx0[2] * dx0[2]) * dx01[1] + (-dx0[0] * dx0[0] * dx1[2] + dx0[0] * dx1[0] * dx0[2] - dx0[1] * dx0[1] * dx1[2] + dx0[1] * dx1[1] * dx0[2]) * dx01[2];

    B *= one_o_denom;

    // A/B are fractions along respective edge axes where edges are closest
    if ((A > 0.0 && A < 1.0) && (B > 0.0 && B < 1.0)) {

      double xp_a[3];
      FOR_I3 xp_a[i] = v00[i] + A*dx0[i];
      double xp_b[3];
      FOR_I3 xp_b[i] = v10[i] + B*dx1[i];

      return DIST2(xp_a,xp_b);
    }
    else if (A <= 0.0) {
      return getPointToEdgeDist2(v00,v10,v11);
    }
    else if (A >= 1.0) {
      return getPointToEdgeDist2(v01,v10,v11);
    }
    else if (B <= 0.0) {
      return getPointToEdgeDist2(v10,v00,v01);
    }
    else {
      assert(B >= 1.0);
      return getPointToEdgeDist2(v11,v00,v01);
    }
  }

  double getPointToTriDist2(const double xp[3],const double v0[3],const double v1[3],const double v2[3]) {

    double xc[3];
    getClosestPointOnTriRobust(xc,xp,v0,v1,v2);

    return ( (xp[0]-xc[0])*(xp[0]-xc[0]) +
	     (xp[1]-xc[1])*(xp[1]-xc[1]) +
	     (xp[2]-xc[2])*(xp[2]-xc[2]) );
  }

  void getTriUnitNormal(double (&un)[3],const double v0[3],const double v1[3],const double v2[3]) {
    const double dx1[3] = DIFF(v1,v0);
    const double dx2[3] = DIFF(v2,v0);
    un[0] = CROSS_PRODUCT_0(dx1,dx2);
    un[1] = CROSS_PRODUCT_1(dx1,dx2);
    un[2] = CROSS_PRODUCT_2(dx1,dx2);
    NORMALIZE(un);
  }

  void getClosestPointOnTriRobust(double xc[3],const double xp[3],const double v0[3],const double v1[3],const double v2[3],const bool debug) {

    // this one robust to bad tris...
    // should eventually modify the below routine in the same way...

    // compute Facet Triangle Edges: T(s,t) = v0 + s*edge0 + t*edge1
    const double edge0[3] = DIFF(v1,v0);
    const double edge1[3] = DIFF(v2,v0);

    // compute coefficients of distance squared function, Q(s,t)
    const double qa = DOT_PRODUCT(edge0,edge0);
    const double qc = DOT_PRODUCT(edge1,edge1);
    const double d12 = DIST2(v2,v1);
    const double l2 = max(qa,max(qc,d12));
    const double eps = 1.0E-15*l2;

    if (debug) cout << "qa,qc,d12,eps: " << qa << " " << qc << " " << d12 << " " << eps << endl;

    if ((qa <= eps)&&(qc <= eps)) {
      // this is a zero tri: all points must be the same, so return distance to the
      // first point...
      FOR_I3 xc[i] = v0[i];
      return;
    }
    else if ((qa <= eps)&&(d12 <= eps)) {
      FOR_I3 xc[i] = v1[i];
      return;
    }
    else if ((qc <= eps)&&(d12 <= eps)) {
      FOR_I3 xc[i] = v2[i];
      return;
    }

    const double qb = DOT_PRODUCT(edge0,edge1);

    // compute coordinates Qmin = Q(sbar/det,tbar/det)
    double det = qa*qc - qb*qb;

    if (debug) cout << "qb,det,eps2: " << qb << " " << det << " " << eps*l2 << endl;

    if (det <= eps*l2) {
      if (d12 > max(qa,qc)) {
	getClosestPointOnEdgeRobust(xc,xp,v1,v2);
      }
      else if (qa >= qc) {
	getClosestPointOnEdgeRobust(xc,xp,v0,v1);
      }
      else {
        getClosestPointOnEdgeRobust(xc,xp,v0,v2);
      }
      return;
    }

    double diff_v0xp[3];
    FOR_I3 diff_v0xp[i] = (v0[i]) - (xp[i]);
    const double qd = DOT_PRODUCT(edge0,diff_v0xp);
    const double qe = DOT_PRODUCT(edge1,diff_v0xp);
    //const double qf = DOT_PRODUCT(diff_v0xp,diff_v0xp);
    double sbar = qb*qe - qc*qd;
    double tbar = qb*qd - qa*qe;

    /***********************************
     * Identify region containing xp
     *        t
     *      \2|
     *       \|
     *        \
     *        |\
     *        | \
     *      3 | 0\  1
     *     ---|---\-----s
     *        |    \
     *      4 | 5   \  6
     ************************************/
    int region = -1;
    if (sbar + tbar <= det){
      if (sbar < 0.0){
        if (tbar < 0.0) region = 4;
        else region = 3;
      }
      else if (tbar < 0.0) region = 5;
      else region = 0;
    }
    else{
      if (sbar < 0.0) region = 2;
      else if (tbar < 0.0) region = 6;
      else region = 1;
    }
    assert(region>-1);

    // facet coordinates that give min distance to xp
    // will be assigned to sbar, tbar

    double invDet, numer, denom, tmp0, tmp1;
    switch (region) {
    case 0:
      invDet = 1.0/det;
      sbar *= invDet;
      tbar *= invDet;
      break;
    case 1:
      numer = (qc+qe-qb-qd);
      if (numer <= 0.0) sbar = 0.0;
      else {
	denom = qa - 2.0*qb + qc; // positive
	sbar = ( numer >= denom ? 1.0 : numer/denom );
      }
      tbar = 1.0 - sbar;
      break;
    case 2:
      tmp0 = qb + qd;
      tmp1 = qc + qe;
      if (tmp1 > tmp0){ // minimum on edge s+t = 1
	numer = tmp1 - tmp0;
	denom = qa - 2.0*qb + qc;
	sbar = ( numer >= denom ? 1.0 : numer/denom);
	tbar = 1.0-sbar;
      }
      else { // minimum on edge s = 0
	sbar = 0;
	tbar = ( tmp1 <= 0.0 ? 1.0 : ( qe >= 0.0 ? 0.0 : -qe/qc ) );
      }
      break;
    case 3:
      sbar = 0;
      tbar = ( qe >= 0.0 ? 0.0 : ( -qe >= qc ? 1.0 : -qe/qc ) );
      break;
    case 4:
      if (qe < 0.0) { // minimum on edge s = 0
	sbar = 0.0;
	tbar = ( qc <= -qe ? 1.0 : -qe/qc );
      }
      else { // minimum on edge t = 0
	tbar = 0.0;
	sbar = ( qa <= -qd ? 1.0 : ( qd >= 0.0 ? 0.0 : -qd/qa ) );
      }
      break;
    case 5:
      tbar = 0.0;
      sbar = ( qd >= 0.0 ? 0.0 : ( -qd >= qa ? 1.0 : -qd/qa ) );
      break;
    case 6:
      tmp0 = qa + qd;
      tmp1 = qb + qe;
      if (tmp0 > tmp1){ // minimum on edge s+t = 1
	numer = tmp0 - tmp1;
	denom = qa - 2.0*qb + qc;
	tbar = ( numer >= denom ? 1.0 : numer/denom );
	sbar = 1.0-tbar;
      }
      else{ // minimum on edge t = 0
	tbar = 0.0;
	sbar = ( tmp0 <= 0.0 ? 1.0 : ( qd >= 0.0 ? 0.0 : -qd/qa ) );
      }
      break;
    }

    if (debug) cout << "region,qd,qe,sbar,tbar,tmp0,tmp1,numer,denom: " << region << qd << " " << qe << " " << sbar << " " << tbar << " " << tmp0 << " " << tmp1 << " " << numer << " " << denom << endl;

    // return point coordinates in physical space
    FOR_I3 xc[i] = sbar*edge0[i] + tbar*edge1[i] + v0[i];

  }

  void getClosestPointOnTriInterior(double xc[3],const double xp[3],const double v0[3],const double v1[3],const double v2[3]) {

    const double eps = 1.0E-8;

    // NOTE: same as the Robust version above, but biased at the end
    // so it is never exactly along the tri edges when possible.

    //Compute Facet Triangle Edges: T(s,t) = v0 + s*edge0 + t*edge1
    const double edge0[3] = DIFF(v1,v0);
    const double edge1[3] = DIFF(v2,v0);

    //Compute coefficients of distance squared function, Q(s,t)
    const double qa = DOT_PRODUCT(edge0,edge0);
    const double qc = DOT_PRODUCT(edge1,edge1);

    if ((qa <= 0.0)&&(qc <= 0.0)) {
      // this is a zero tri: all points must be the same, so return the first point...
      FOR_I3 xc[i] = v0[i];
      return;
    }

    const double qb = DOT_PRODUCT(edge0,edge1);
    const double det = qa*qc - qb*qb;

    if (det <= 0.0) {
      const double edge12[3] = DIFF(v2,v1);
      const double edge12_d2 = DOT_PRODUCT(edge12,edge12);
      if (edge12_d2 > max(qa,qc)) {
	// closest point is along v1,v2...
	const double dxp[3] = DIFF(xp,v1);
	const double dp = DOT_PRODUCT(edge12,dxp);
	if (dp <= 0.0) {
	  // we are closest to the first point v1...
	  FOR_I3 xc[i] = (1.0-2.0*eps)*v1[i] + eps*(v0[i] + v2[i]);
	}
	else if (dp >= edge12_d2) {
	  // we are closest to the second point v2...
	  FOR_I3 xc[i] = (1.0-2.0*eps)*v2[i] + eps*(v0[i] + v1[i]);
	}
	else {
	  const double w1 = max(eps,1.0-dp/edge12_d2);
	  const double w2 = max(eps,dp/edge12_d2);
	  const double w0 = eps;
	  const double wsum = w0 + w1 + w2;
	  FOR_I3 xc[i] = (w0*v0[i] + w1*v1[i] + w2*v2[i])/wsum;
	}
      }
      else {
	const double dxp[3] = DIFF(xp,v0);
	if (qa >= qc) {
	  // closest point is along v0,v1...
	  const double dp = DOT_PRODUCT(edge0,dxp);
	  if (dp <= 0.0) {
	    // we are closest to the first point v0...
	    FOR_I3 xc[i] = (1.0-2.0*eps)*v0[i] + eps*(v1[i] + v2[i]);
	  }
	  else if (dp >= qa) {
	    // we are closest to the second point v1...
	    FOR_I3 xc[i] = (1.0-2.0*eps)*v1[i] + eps*(v0[i] + v2[i]);
	  }
	  else {
	    const double w0 = max(eps,1.0-dp/qa);
	    const double w1 = max(eps,dp/qa);
	    const double w2 = eps;
	    const double wsum = w0 + w1 + w2;
	    FOR_I3 xc[i] = (w0*v0[i] + w1*v1[i] + w2*v2[i])/wsum;
	  }
	}
	else {
	  // closest point is along v0,v2...
	  const double dp = DOT_PRODUCT(edge1,dxp);
	  if (dp <= 0.0) {
	    // we are closest to the first point v0...
	    FOR_I3 xc[i] = (1.0-2.0*eps)*v0[i] + eps*(v1[i] + v2[i]);
	  }
	  else if (dp >= qc) {
	    // we are closest to the second point v1...
	    FOR_I3 xc[i] = (1.0-2.0*eps)*v1[i] + eps*(v0[i] + v2[i]);
	  }
	  else {
	    const double w0 = max(eps,1.0-dp/qc);
	    const double w2 = max(eps,dp/qc);
	    const double w1 = eps;
	    const double wsum = w0 + w1 + w2;
	    FOR_I3 xc[i] = (w0*v0[i] + w1*v1[i] + w2*v2[i])/wsum;
	  }
	}
      }
      return;
    }

    const double diff_v0xp[3] = DIFF(v0,xp);
    const double qd = DOT_PRODUCT(edge0,diff_v0xp);
    const double qe = DOT_PRODUCT(edge1,diff_v0xp);
    //const double qf = DOT_PRODUCT(diff_v0xp,diff_v0xp);
    double sbar = qb*qe - qc*qd;
    double tbar = qb*qd - qa*qe;

    /***********************************
     * Identify region containing xp
     *        t
     *      \2|
     *       \|
     *        \
     *        |\
     *        | \
     *      3 | 0\  1
     *     ---|---\-----s
     *        |    \
     *      4 | 5   \  6
     ************************************/
    int region = -1;
    if (sbar + tbar <= det){
      if (sbar<0){
        if (tbar<0) region = 4;
        else region = 3;
      }
      else if (tbar<0) region = 5;
      else region = 0;
    }
    else{
      if (sbar<0) region = 2;
      else if (tbar<0) region = 6;
      else region = 1;
    }
    assert(region>-1);

    //Facet coordinates that give min distance to xp
    //will be assigned to sbar, tbar

    double invDet, numer, denom, tmp0, tmp1;
    switch (region) {
    case 0:
      invDet = 1.0/det;
      sbar *= invDet;
      tbar *= invDet;
      break;
    case 1:
      numer = (qc+qe-qb-qd);
      if (numer<=0) sbar = 0;
      else {
	denom = qa - 2*qb + qc; //positive
	sbar = ( numer >= denom ? 1 : numer/denom );
      }
      tbar = 1 - sbar;
      break;
    case 2:
      tmp0 = qb + qd;
      tmp1 = qc + qe;
      if (tmp1 > tmp0){ //minimum on edge s+t = 1
	numer = tmp1 - tmp0;
	denom = qa - 2*qb + qc;
	sbar = ( numer >= denom ? 1 : numer/denom);
	tbar = 1-sbar;
      }
      else { //minimum on edge s = 0
	sbar = 0;
	tbar = ( tmp1 <= 0 ? 1 : ( qe >= 0 ? 0 : -qe/qc ) );
      }
      break;
    case 3:
      sbar = 0;
      tbar = ( qe >= 0 ? 0 : ( -qe >= qc ? 1 : -qe/qc ) );
      break;
    case 4:
      if (qe < 0){ //minimum on edge s = 0
	sbar = 0;
	tbar = ( qc <= -qe ? 1 : -qe/qc );
      }
      else { //minimum on edge t = 0
	tbar = 0;
	sbar = ( qa <= -qd ? 1 : ( qd >= 0 ? 0 : -qd/qa ) );
      }
      break;
    case 5:
      tbar = 0;
      sbar = ( qd >= 0 ? 0 : ( -qd >= qa ? 1 : -qd/qa ) );
      break;
    case 6:
      tmp0 = qa + qd;
      tmp1 = qb + qe;
      if (tmp0 > tmp1){//minimum on edge s+t = 1
	numer = tmp0 - tmp1;
	denom = qa - 2*qb + qc;
	tbar = ( numer >= denom ? 1 : numer/denom );
	sbar = 1-tbar;
      }
      else{ //minimum on edge t = 0
	tbar = 0;
	sbar = ( tmp0 <= 0 ? 1 : ( qd >= 0 ? 0 : -qd/qa ) );
      }
      break;
    }

    //Return point coordinates in physical space
    // here we bias away from the corner and edges slightly to ensure we
    // truly get the closest tri...

    const double w0 = max(eps,1.0-sbar-tbar);
    const double w1 = max(eps,sbar);
    const double w2 = max(eps,tbar);
    const double wsum = w0 + w1 + w2;
    FOR_I3 xc[i] = (w0*v0[i] + w1*v1[i] + w2*v2[i])/wsum;

  }

  bool getClosestPointOnTri(double xc[3],const double xp[3],const double v0[3],const double v1[3],const double v2[3]) {

    // ==============================================================================
    // returns the closes point in xc, and the distance squared in returned double
    // the return value is true when the closest point is in region 0 -- i.e. the
    // closest point is formally inside the tri, not on the edge or corner...
    // ==============================================================================

    //Compute Facet Triangle Edges: T(s,t) = v0 + s*edge0 + t*edge1
    const double edge0[3] = DIFF(v1,v0);
    const double edge1[3] = DIFF(v2,v0);

    //Compute coefficients of distance squared function, Q(s,t)
    double diff_v0xp[3];
    double qa, qb, qc, qd, qe; //, qf;

    FOR_I3 diff_v0xp[i] = (v0[i]) - (xp[i]);
    qa = DOT_PRODUCT(edge0,edge0);
    qb = DOT_PRODUCT(edge0,edge1);
    qc = DOT_PRODUCT(edge1,edge1);
    qd = DOT_PRODUCT(edge0,diff_v0xp);
    qe = DOT_PRODUCT(edge1,diff_v0xp);
    //qf = DOT_PRODUCT(diff_v0xp,diff_v0xp);

    //Compute coordinates Qmin = Q(sbar/det,tbar/det)
    double det, sbar, tbar;
    det = qa*qc - qb*qb;
    if (det <= 0.0) {
      // edges of facet are colinear...
      // SEE above routine that returns d2 for a more robust version of this.
      // F. Ham Jan 2017
      cout << "Warning: edges of facet are potentially colinear: " << det << endl;
      cout << " > edge0: " << edge0[0] << " " << edge0[1] << " " << edge0[2] << endl;
      cout << " > edge1: " << edge1[0] << " " << edge1[1] << " " << edge1[2] << endl;
      const double n[3] = TRI_NORMAL_2(v0,v1,v2);
      double fa_area = 0.5*MAG(n);
      cout << " > facet area: " << fa_area << endl;
      if (fa_area > 1.0e-20) {
	cout << "    : based on the area we are going to keep this tri and not consider it collapsed" << endl;
      }
      else {
	cout << "Error! : based on the area we are going to consider it collapsed" << endl;
        throw(0);
      }
    }
    sbar = qb*qe - qc*qd;
    tbar = qb*qd - qa*qe;

    /***********************************
     * Identify region containing xp
     *        t
     *      \2|
     *       \|
     *        \
     *        |\
     *        | \
     *      3 | 0\  1
     *     ---|---\-----s
     *        |    \
     *      4 | 5   \  6
     ************************************/
    int region = -1;
    if (sbar + tbar <= det){
      if (sbar<0){
        if (tbar<0) region = 4;
        else region = 3;
      }
      else if (tbar<0) region = 5;
      else region = 0;
    }
    else{
      if (sbar<0) region = 2;
      else if (tbar<0) region = 6;
      else region = 1;
    }
    assert(region>-1);
    //Facet coordinates that give min distance to xp
    //will be assigned to sbar, tbar
    double invDet, numer, denom, tmp0, tmp1;
    switch (region) {
    case 0:
      invDet = 1/det;
      sbar *= invDet;
      tbar *= invDet;
      break;
    case 1:
      numer = (qc+qe-qb-qd);
      if (numer<=0) sbar = 0;
      else {
	denom = qa - 2*qb + qc; //positive
	sbar = ( numer >= denom ? 1 : numer/denom );
      }
      tbar = 1 - sbar;
      break;
    case 2:
      tmp0 = qb + qd;
      tmp1 = qc + qe;
      if (tmp1 > tmp0){ //minimum on edge s+t = 1
	numer = tmp1 - tmp0;
	denom = qa - 2*qb + qc;
	sbar = ( numer >= denom ? 1 : numer/denom);
	tbar = 1-sbar;
      }
      else { //minimum on edge s = 0
	sbar = 0;
	tbar = ( tmp1 <= 0 ? 1 : ( qe >= 0 ? 0 : -qe/qc ) );
      }
      break;
    case 3:
      sbar = 0;
      tbar = ( qe >= 0 ? 0 : ( -qe >= qc ? 1 : -qe/qc ) );
      break;
    case 4:
      if (qe < 0){ //minimum on edge s = 0
	sbar = 0;
	tbar = ( qc <= -qe ? 1 : -qe/qc );
      }
      else { //minimum on edge t = 0
	tbar = 0;
	sbar = ( qa <= -qd ? 1 : ( qd >= 0 ? 0 : -qd/qa ) );
      }
      break;
    case 5:
      tbar = 0;
      sbar = ( qd >= 0 ? 0 : ( -qd >= qa ? 1 : -qd/qa ) );
      break;
    case 6:
      tmp0 = qa + qd;
      tmp1 = qb + qe;
      if (tmp0 > tmp1){//minimum on edge s+t = 1
	numer = tmp0 - tmp1;
	denom = qa - 2*qb + qc;
	tbar = ( numer >= denom ? 1 : numer/denom );
	sbar = 1-tbar;
      }
      else{ //minimum on edge t = 0
	tbar = 0;
	sbar = ( tmp0 <= 0 ? 1 : ( qd >= 0 ? 0 : -qd/qa ) );
      }
      break;
    }

    //Return point coordinates in physical space
    FOR_I3 xc[i] = sbar*edge0[i] + tbar*edge1[i] + v0[i];

    return(region==0);

  }

  double spacingFunctionDoubleSided(const double xi,const double stretch,const double x0,const double x1) {
    assert((xi >= 0.0)&&(xi <= 1.0));
    if (xi <= 0.5) return MiscUtils::spacingFunction(2.0*xi,stretch,x0,x0+0.5*(x1-x0));
    else return MiscUtils::spacingFunction(2.0*(1.0-xi),stretch,x1,x0+0.5*(x1-x0));
  }

  double spacingFunction(const double xi,const double stretch,const double x0,const double x1) {
    assert((xi >= 0.0)&&(xi <= 1.0));
    return( x0 + (x1-x0)*xi*( 2.0*stretch + (1.0-stretch)*xi )/(stretch + 1.0) );
  }


  /*
   * given a point and normal, and potentially a hint for primary orthogonal vector, get the orthogonal plane vectors
   */
  void getOrthogonalVectors(double (&radDir1)[3], double (&radDir2)[3],const double xc[3],const double np[3], const double hint[3]) {
    double dir1[3] = {hint[0],hint[1],hint[2]};
    if (MAG(hint) == 0.0) {
      // user hasn't specified a hint direction; so compute
      int dir1Index = 0;
      if ( ( std::abs(np[1]) <= std::abs(np[0]) ) && ( std::abs(np[1]) <= std::abs(np[2]) ) ) {
        dir1Index = 1;
      }
      if ( ( std::abs(np[2]) <= std::abs(np[0]) ) && ( std::abs(np[2]) <= std::abs(np[1]) ) ) {
        dir1Index = 2;
      }
      if ( dir1Index == 0 ) { dir1[0] = 1.0; dir1[1] = 0.0; dir1[2] = 0.0; };
      if ( dir1Index == 1 ) { dir1[0] = 0.0; dir1[1] = 1.0; dir1[2] = 0.0; };
      if ( dir1Index == 2 ) { dir1[0] = 0.0; dir1[1] = 0.0; dir1[2] = 1.0; };
    }
    double _radDir1[3] = CROSS_PRODUCT( np , dir1 );
    double mag = DOT_PRODUCT( _radDir1 , _radDir1 );
    FOR_I3 {
      _radDir1[i] /= sqrt( mag );
      radDir1[i] = _radDir1[i];
    };
    double _radDir2[3] = CROSS_PRODUCT( _radDir1 , np );
    mag = DOT_PRODUCT( _radDir2 , _radDir2 );
    FOR_I3 {
      _radDir2[i] /= - sqrt( mag );
      radDir2[i] = _radDir2[i];
    };

    if (MAG(hint) != 0.0) {
      FOR_I3 {
	const double tmp = radDir1[i];
	radDir1[i] = radDir2[i];
        radDir2[i] = -tmp;
      }
    }
  }

  /*
   * given the centroid and normal, generate points on the circle with a specified radius
   */
  void createCirclePts(double (* xp)[3],const int index0,const double xc[3],const double np[3],const double radius, const int n,const bool stagger) {

    double radDir1[3];
    double radDir2[3];
    const double primary_dir[3] = {0.0,0.0,0.0};

    MiscUtils::getOrthogonalVectors(radDir1,radDir2,xc,np,primary_dir);
    // populate circle of points into xp at specified index (clockwise ordering)
    double offset = 0.0;
    if (stagger) offset = 0.5;
    for (int i = 0; i < n; ++i) {
      const double theta0 = (double(i)+offset)/double(n)*2.0*M_PI;
      FOR_J3 {
        xp[index0+i][j] = xc[j] + radius * radDir1[j] * cos(theta0) + radius * radDir2[j] * sin(theta0);
      }
    }
  }

  /*
   * given the centroid and normal and major/minor axis radii, generate points on the ellipse
   */
  void createEllipsePts(double (* xp)[3],const int index0,const double xc[3],const double np[3],const double rM, const double rm,const double primary_dir[3], const int n) {

    double radDir1[3];
    double radDir2[3];

    MiscUtils::getOrthogonalVectors(radDir1,radDir2,xc,np,primary_dir);
    // populate ellipse with points into xp at specified index (clockwise ordering)
    for (int i = 0; i < n; ++i) {
      const double theta0 = (double(i)+0.5)/double(n)*2.0*M_PI;
      const double radius = rM*rm/(sqrt(rm*rm*cos(theta0)*cos(theta0) + rM*rM*sin(theta0)*sin(theta0)));

      FOR_J3 {
        xp[index0+i][j] = xc[j] + radius * radDir1[j] * cos(theta0) + radius * radDir2[j] * sin(theta0);
      }
    }
  }

  /*
   * given the centroid and normal, convert 2D (xy profile) to 3D points in the plane defined by pt/normal
   */
  void create3DFrom2DPtAndNormal(double (* xp)[3],const int index0,const double xc[3],const double np[3],double (*xy_vals)[2],const int n) {
    int dir1Index = 0;
    if ( ( std::abs(np[1]) <= std::abs(np[0]) ) && ( std::abs(np[1]) <= std::abs(np[2]) ) ) {
      dir1Index = 1;
    }
    if ( ( std::abs(np[2]) <= std::abs(np[0]) ) && ( std::abs(np[2]) <= std::abs(np[1]) ) ) {
      dir1Index = 2;
    }
    double dir1[3];
    if ( dir1Index == 0 ) { dir1[0] = 1.0; dir1[1] = 0.0; dir1[2] = 0.0; };
    if ( dir1Index == 1 ) { dir1[0] = 0.0; dir1[1] = 1.0; dir1[2] = 0.0; };
    if ( dir1Index == 2 ) { dir1[0] = 0.0; dir1[1] = 0.0; dir1[2] = 1.0; };

    // create 3D vectors of in-plane x,y (e0,e1)
    double e0[3] = CROSS_PRODUCT( np , dir1 );
    double mag = DOT_PRODUCT( e0 , e0 );
    FOR_I3 {
      e0[i] /= sqrt( mag );
    };
    double e1[3] = CROSS_PRODUCT( e0 , np );
    mag = DOT_PRODUCT( e1 , e1 );
    FOR_I3 { e1[i] /= - sqrt( mag ); };

    // populate profile
    for (int i = 0; i < n; ++i) {
      FOR_J3 {
        xp[index0+i][j] = xc[j] - xy_vals[i][0]*e0[j] - xy_vals[i][1]*e1[j];
      }
    }
  }



  int getClosestPoint(const double xp_check[3],const double (* const xp)[3],const int np) {

    double d2_closest = HUGE_VAL;
    int ip_closest = -1;
    for (int ip = 0; ip < np; ++ip) {
      const double this_d2 = DIST2(xp[ip],xp_check);
      if (this_d2 < d2_closest) {
	d2_closest = this_d2;
	ip_closest = ip;
      }
    }

    // we need to get the minimum value AND the winning rank: use
    // MPI's MPI_MINLOC capability...

    DoubleInt myDi,di;

    myDi.this_double = d2_closest;
    myDi.this_int = mpi_rank;

    MPI_Allreduce(&myDi,&di,1,MPI_DOUBLE_INT,MPI_MINLOC,mpi_comm);

    if (di.this_int == mpi_rank) {
      assert(ip_closest >= 0);
      return ip_closest;
    }

    // we do not own the closest...
    return -1;

  }

  int lineSphereIntersections(float (&intersections)[2][3], const float line_xc[3], const float _line_n[3],const float sphere_xc[3], const float sphere_r) {
    float line_n[3] = {_line_n[0],_line_n[1],_line_n[2]};
    NORMALIZE(line_n);
    const float offset[3] = DIFF(line_xc,sphere_xc);
    const float term1 = -DOT_PRODUCT(line_n,offset);

    float value = term1*term1 - DOT_PRODUCT(offset,offset) + sphere_r*sphere_r;

    if (value < 0) return 0;
    else if (value == 0.0) {
      FOR_I3 intersections[0][i] = line_xc[i] + term1*line_n[i];
      return 1;
    }
    else {
      float sqrt_value = sqrt(value);
      FOR_I3 intersections[0][i] = line_xc[i] + (term1 + sqrt_value)*line_n[i];
      FOR_I3 intersections[1][i] = line_xc[i] + (term1 - sqrt_value)*line_n[i];
      return 2;
    }
  }

  void computePlaneCoeffsFromPoints(double (&N0)[3], double& d0, const double v0[3],const double v1[3],const double v2[3]) {
    // compute plane equation coeffs (N,d) where N*x + d = 0 for points x on the plane
    // combine info from all three edge pairs for more robust normal/offset

    double d2[3];
    d2[0] = DIST2(v0,v1);
    d2[1] = DIST2(v0,v2);
    d2[2] = DIST2(v1,v2);

    int d2_min_index = 0;
    for (int i=1; i<3; ++i) {
      if (d2[i] < d2[d2_min_index]) d2_min_index = i;
    }

    // take cross product of largest vectors to determine plane normal
    if (d2_min_index == 0) {
      const double d_v1[3] = DIFF(v0,v2);
      const double d_v2[3] = DIFF(v1,v2);
      double N0_temp[3] = CROSS_PRODUCT(d_v1,d_v2);
      FOR_I3 N0[i] = N0_temp[i];
      FOR_I3 N0_temp[i] *= -1.0;  // make negative to get displacement from origin
      d0 = DOT_PRODUCT(N0_temp,v2);
    }
    else if (d2_min_index == 1) {
      const double d_v1[3] = DIFF(v2,v1);
      const double d_v2[3] = DIFF(v0,v1);
      double N0_temp[3] = CROSS_PRODUCT(d_v1,d_v2);
      FOR_I3 N0[i] = N0_temp[i];
      FOR_I3 N0_temp[i] *= -1.0;  // make negative to get displacement from origin
      d0 = DOT_PRODUCT(N0_temp,v1);
    }
    else if (d2_min_index == 2) {
      const double d_v1[3] = DIFF(v2,v0);
      const double d_v2[3] = DIFF(v1,v0);
      double N0_temp[3] = CROSS_PRODUCT(d_v1,d_v2);
      FOR_I3 N0[i] = N0_temp[i];
      FOR_I3 N0_temp[i] *= -1.0;  // make negative to get displacement from origin
      d0 = DOT_PRODUCT(N0_temp,v0);
    }

  }

  int computeLinePlaneIntersection(double (&x_intersection)[3],const double xl[3], const double nl[3],const double xp[3],const double np[3]) {
    // check denominator value, make sure line isn't parallel to plane
    if (MAG(nl) <= 0.0 || MAG(np) <= 0.0) {
      // CWARN("invalid normal for plane or line");
      return -1;  // invalid
    }

    //double x1_line[3];
    //FOR_I3 x1_line[i] = xl[i]+nl[i];
    //const double denom_diff[3] = DIFF(x1_line,xl);
    //const double denom_val=DOT_PRODUCT(np,denom_diff);
    const double denom_val=DOT_PRODUCT(np,nl);

    if (fabs(denom_val) == 0.0) {
      // CWARN("line is parallel to plane-of-interest; no intersection");
      return 1;  // parallel
    }

    const double numer_diff[3] = DIFF(xp,xl);
    const double numer_val=DOT_PRODUCT(np,numer_diff);

    const double fac = numer_val/denom_val;

    //FOR_I3 x_intersection[i] = xl[i] + fac*denom_diff[i];
    FOR_I3 x_intersection[i] = xl[i] + fac*nl[i];
    return 0;  // intersection found

  }

  int howDoesLinePierceTri(const double v0[3],const double v1[3],const double v2[3],const double a[3],const double b[3], const double dist_tol) {
    // assume (for now) user has checked to see if either a or b lies on tri, so no imprint possible

    // categorize as either:
    // 0: does not intersect triangle
    // 1: does intersect triangle interior
    // -1,-2,-3: node at which line intersects
    // 2,3,4: value-1 is edge at which line intersects
    const double d2_tol = dist_tol*dist_tol;
    double d2[3];
    d2[0] = getEdgeToEdgeDist2(a,b,v0,v1);
    d2[1] = getEdgeToEdgeDist2(a,b,v1,v2);
    d2[2] = getEdgeToEdgeDist2(a,b,v2,v0);
    if (d2[0] <= d2_tol && d2[1] <= d2_tol) return -2;
    if (d2[0] <= d2_tol && d2[2] <= d2_tol) return -1;
    if (d2[1] <= d2_tol && d2[2] <= d2_tol) return -3;

    for (int i=0; i<3; ++i) {
      if (d2[i] <= d2_tol) return (i+2);
    }

    // above checks for edge/node intersection
    // we are sufficiently far enough away that volume should be adequate check (and should be non-zero)

    double tet_vol[3];
    tet_vol[0] = SIGNED_TET_VOLUME_6(a,v0,v1,b);
    tet_vol[1] = SIGNED_TET_VOLUME_6(a,v1,v2,b);
    tet_vol[2] = SIGNED_TET_VOLUME_6(a,v2,v0,b);

    int countPos = 0;
    int countNeg = 0;
    int countZero = 0;
    for (int i=0; i<3; ++i) {
      if (tet_vol[i] < 0.0) {
        ++countNeg;
      }
      else if (tet_vol[i] > 0.0) {
        ++countPos;
      }
      else {
        ++countZero;
      }
    }
    // assert(countZero == 0);

    //NOTE can get nearly coplanar and countZero == 3, but edge-edge checks above would have found edge intersection. Do we need a check on point-in-tri here?

    if (countPos==3 || countNeg==3) return 1;
    else return 0;
  }

  // replace all occurances of "find" by "replace" in "str" and return it
  std::string replaceAll(const std::string str ,const std::string find ,const std::string replace, const size_t start) {
    using namespace std;
    string result = str.substr(0,start);
    size_t find_len = find.size();
    size_t pos,from=start;
    while ( string::npos != ( pos=str.find(find,from) ) ) {
      result.append( str, from, pos-from );
      result.append( replace );
      from = pos + find_len;
    }
    result.append( str, from , string::npos );
    return result;
  }

  // replace last occurance of "find" by "replace" in "str"
  std::string replaceLast(const std::string str ,const std::string find ,const std::string replace) {
    using namespace std;

    size_t found = str.find_last_of(find);

    if (found < string::npos) {
      size_t find_len = find.size();
      string result = str.substr(0,found);
      result.append(replace);
      result.append(str,found+find_len,string::npos);
      return result;
    }
    else {
      return str;
    }
  }

  std::string cleanFilename(const std::string& str, const size_t start) {
    string result;
    result = replaceAll(str,"/","_d_",start);  // replace division symbols from filename
    return result;
  }

  bool getBarycentricCoordinates(double lambda[4],const double xp[3],const double x0[3],
                                 const double x1[3], const double x2[3],const double x3[3]) {
    const double v = SIGNED_TET_VOLUME_6(x0,x1,x2,x3);
    assert(v != 0.0);
    const double inv_v = 1.0/v;
    lambda[0] = SIGNED_TET_VOLUME_6(xp,x1,x2,x3)*inv_v;
    lambda[1] = SIGNED_TET_VOLUME_6(x0,xp,x2,x3)*inv_v;
    lambda[2] = SIGNED_TET_VOLUME_6(x0,x1,xp,x3)*inv_v;
    lambda[3] = SIGNED_TET_VOLUME_6(x0,x1,x2,xp)*inv_v;

    bool b_inside = true;
    FOR_I4 {
      if ((lambda[i] < 0.0)||(lambda[i] > 1.0)) {
        b_inside = false;
        break;
      }
    }
    return b_inside;
  }

  // compile-time endianness (see Defs.hpp)
  int getEndianness() {
    union {
      int i;
      char c[sizeof(int)];
    } x;
    x.i = 1;
    if (x.c[0] == 1)
      return CTI_LITTLE_ENDIAN; // little
    else
      return CTI_BIG_ENDIAN;
  }
  //int getEndianness() {
  //  unsigned char EndianTest[2] = {1,0};
  //  short endian = *(short *)EndianTest;
  //  if (endian != CTI_LITTLE_ENDIAN) endian = CTI_BIG_ENDIAN;
  //  return endian;
  //}
  //int getEndianness() {
  //  int num = 1;
  //  if (*(char*)&num == 1)
  //    return CTI_LITTLE_ENDIAN;
  //  else
  //    return CTI_BIG_ENDIAN;
  //}

  int solveCgSerial(double * phi,const int ncv,const double * const A,const double * const rhs,
                    const int * const cvocv_i,const int * const cvocv_v,
                    const double zero,const int maxiter,const bool verbose) {

    // just use the pased-in initial condition in phi...

    // we need the following work arrays...

    double * res      = new double[ncv];
    double * v        = new double[ncv];
    double * p        = new double[ncv];
    double * inv_diag = new double[ncv];

    // initialize...
    for (int icv = 0; icv < ncv; ++icv)  {
      assert(cvocv_v[cvocv_i[icv]] == icv); // confirm diagonal first in cvocv_i/v CSR
      assert(A[cvocv_i[icv]] != 0.0); // diag cannot be zero for diag-preconditioned cg.
      inv_diag[icv] = 1.0/A[cvocv_i[icv]];
    }

    for (int icv = 0; icv < ncv; ++icv)
      p[icv] = 0.0;
    double rho = 1.0;

    // calculate the residual in rhs format...
    for (int icv = 0; icv < ncv; ++icv) {
      res[icv] = rhs[icv] - A[cvocv_i[icv]]*phi[icv];
      for (int coc = cvocv_i[icv]+1; coc < cvocv_i[icv+1]; ++coc) {
	const int icv_nbr = cvocv_v[coc];
	res[icv] -= A[coc]*phi[icv_nbr];
      }
    }

    // diagonal precon/compute normalized residual...

    double  res_max_prev = 0.0;
    for (int icv = 0; icv < ncv; ++icv) {
      v[icv] = res[icv]*inv_diag[icv];
      res_max_prev = max( res_max_prev, fabs(v[icv]) );
    }
    cout << " > initial res_max: " << res_max_prev << endl;

    int iter = 0;
    int done = 0;
    while (done == 0) {

      ++iter;

      assert(rho != 0.0);
      const double rho_prev = rho;
      //if (fabs(rho_prev) < 1.0E-20)
      //rho_prev = -1.0E-20; // -1.0E-20? seems to help

      rho = 0.0;
      for (int icv = 0; icv < ncv; ++icv)
	rho += res[icv]*v[icv];

      assert(rho_prev != 0.0);
      const double beta = rho/rho_prev;
      for (int icv = 0; icv < ncv; ++icv)
	p[icv] = v[icv] + beta*p[icv];

      // v = [Ap]{p}...
      for (int icv = 0; icv < ncv; ++icv) {
	v[icv] = A[cvocv_i[icv]]*p[icv];
	for (int coc = cvocv_i[icv]+1; coc < cvocv_i[icv+1]; ++coc) {
	  const int icv_nbr = cvocv_v[coc];
	  v[icv] += A[coc]*p[icv_nbr];
	}
      }

      double gamma = 0.0;
      for (int icv = 0; icv < ncv; ++icv)
	gamma += p[icv]*v[icv];
      assert(gamma != 0.0);
      assert(gamma == gamma); // nan check

      const double alpha = rho/gamma;

      // check if we are done...
      if (iter%3 == 0) {

	for (int icv = 0; icv < ncv; ++icv) {
          if (p[icv] != p[icv]) {
            cout << "nan at icv: " << icv << endl;
          }
	  phi[icv] += alpha*p[icv];
        }

	// recompute the residual...
	for (int icv = 0; icv < ncv; ++icv) {
	  res[icv] = rhs[icv] - A[cvocv_i[icv]]*phi[icv];
	  for (int coc = cvocv_i[icv]+1; coc < cvocv_i[icv+1]; ++coc) {
	    const int icv_nbr = cvocv_v[coc];
	    res[icv] -= A[coc]*phi[icv_nbr];
	  }
	}

	for (int icv = 0; icv < ncv; ++icv)
	  v[icv] = res[icv]*inv_diag[icv];

	// compute the max (L-infinity) normalized residual...
	double  res_max = 0.0;
	for (int icv = 0; icv < ncv; ++icv)
	  res_max = max( res_max, fabs(v[icv]) );
	// only share the last half of the convergence behaviour...
	if (verbose || (iter > maxiter/2))
	  cout << " > solveCgLocal iter " << iter << " res_max " << res_max << endl;
	if (res_max <= zero) {
	  cout << "-> Successfully converged error to " << res_max << endl;
	  done = 1;
	}
        else if (res_max > res_max_prev) {
	  cout << "-> Warning: res_max not converging (zero may be too small): " << res_max << endl;
	  //done = 1;
        }
        else if (iter > maxiter) {
	  cout << "Warning: solveCgLocal did not converge after " << maxiter <<
	    " iters, res_max: " << res_max << endl;
	  done = 2;
	}
        res_max_prev = res_max;

      }
      else {

	// update full phi...
	for (int icv = 0; icv < ncv; ++icv)
	  phi[icv] += alpha*p[icv];

	for (int icv = 0; icv < ncv; ++icv) {
	  // on the other iterations, use this approximation to update
	  // the unreduced residual...
	  res[icv] -= alpha*v[icv];
	  // still need to compute v, diag precon for next iteration...
	  v[icv] = res[icv]*inv_diag[icv];
	}

      }

      //MiscUtils::dumpRange(phi,ncv,"PHI RANGE");
      //getchar();

    }

    delete[] res;
    delete[] v;
    delete[] p;
    delete[] inv_diag;

    // let the calling routine know if we were successful...
    return done;

  }

  int solveGaussSidelSerial(double * phi,const int ncv,const double * const A,const double * const rhs,
                            const int * const cvocv_i,const int * const cvocv_v,
                            const double zero,const int maxiter,const double relax,const bool verbose) {

    // assume we come in with a consistent initial condition...

    // we need the following work arrays...

    int iter = 0;
    int done = 0;
    while (done == 0) {

      ++iter;

      double res_max = 0.0;
      for (int icv = 0; icv < ncv; ++icv) {
        double res = rhs[icv] - A[cvocv_i[icv]]*phi[icv];
        for (int coc = cvocv_i[icv]+1; coc < cvocv_i[icv+1]; ++coc) {
          const int icv_nbr = cvocv_v[coc];
          res -= A[coc]*phi[icv_nbr];
        }
        const double dphi = res/A[cvocv_i[icv]];
        res_max = max(res_max,fabs(dphi));
        phi[icv] += relax*dphi;
      }

      if (verbose || (iter > maxiter/2))
        cout << " > solveGaussSidelSerial iter " << iter << " res_max " << res_max << endl;
      if (res_max <= zero) {
        cout << "-> Successfully converged error to " << res_max << endl;
        done = 1;
      }

    }

    return done;

  }

  // fast extraction of column of data from a text file, with grep and tail options...

  int xcol(double * &data,int &n,char const *fname,const int column) {

    int fd = open(fname, O_RDONLY);
    if (fd == -1) {
      cout << "Error: xcol file open problem for file: " << fname << endl;
      return -1;
    }

    // Advise the kernel of our access pattern: FDADVICE_SEQUENTIAL
    // note that this does not seem to be available on all systems tested...
    //posix_fadvise(fd, 0, 0, 1);

    //static const int BUFFER_SIZE = 16*8;
    static const int BUFFER_SIZE = 16*1024; // size used by "wc -l"...

    char buf[BUFFER_SIZE*2+1];
    size_t bytes_read = read(fd, buf, BUFFER_SIZE*2);
    if (bytes_read <= 0) {
      cout << "Error: xcol read returned bytes_read: " << bytes_read << " for file: " << fname << endl;
      close(fd);
      return -1;
    }

    // experimented with managing memory more directly, but file i/o is so much
    // slower than everything else, so just use a vector, because it is simplest...
    vector<double> doubleVec;

    // if the file is smaller than BUFFER_SIZE*2, then make sure it has '\n'
    // at the end of the last line and set pos_eof to the position after the '\n'...
    int pos_eof = -1;
    if (bytes_read < BUFFER_SIZE*2) {
      if (buf[bytes_read-1] != '\n') {
        buf[bytes_read++] = '\n';
      }
      pos_eof = bytes_read;
    }

    int pos = 0;
    int pos_nl = 0;
    while (pos != pos_eof) {

      // if we are in the back half of the buffer and we are not
      // done reading, consider advancing the buffer...
      if ((pos >= BUFFER_SIZE)&&(bytes_read == BUFFER_SIZE*2)) {
        memcpy(buf,buf+BUFFER_SIZE,BUFFER_SIZE);
        bytes_read = read(fd,buf+BUFFER_SIZE,BUFFER_SIZE);
        if (bytes_read < BUFFER_SIZE) {
          if (bytes_read < 0) {
            cerr << "Error: xcol bytes_read: " << bytes_read << endl;
            assert(0);
          }
          if (buf[BUFFER_SIZE+bytes_read-1] != '\n') {
            buf[BUFFER_SIZE+bytes_read] = '\n';
            ++bytes_read;
          }
          pos_eof = BUFFER_SIZE+bytes_read;
        }
        bytes_read += BUFFER_SIZE;
        pos -= BUFFER_SIZE;
      }

      // use memchr to find the next newline...
      const char * nl = (char*)memchr(buf+pos,'\n',bytes_read-pos);
      if (nl == NULL) {
        // if you hit this error, it means the buffer is not long enough to read the file
        // until a newline is found. This should produce an error?...
        cout << "Error: xcol file does not have a newline in " << bytes_read-pos << " chars. Increase BUFFER_SIZE in xcol." << endl;
        assert(0);
      }
      pos_nl = int(nl-buf);
      assert(pos <= pos_nl);

      // skip any line with first char a hash...
      if (buf[pos] == '#') {
        pos = pos_nl + 1;
        continue;
      }

      // clear leading whitespace once...
      while ((pos < pos_nl)&&((buf[pos]==' ')||(buf[pos]=='\t')))
        ++pos;

      // advance through the first columns to get to column1...
      int ic = 1;
      while (ic < column) {
        // advance through non-whitespace...
        while ((pos < pos_nl)&&(buf[pos]!=' ')&&(buf[pos]!='\t')&&(buf[pos]!=','))
          ++pos;
        bool b_comma = false;
        // then advance through whitespace which may contain up to 1 comma...
        while ((pos < pos_nl)&&((buf[pos]==' ')||(buf[pos]=='\t')||(buf[pos]==','))) {
          if (buf[pos]==',') {
            // if we hit a second comma, then break without advancing pos...
            if (b_comma)
              break;
            b_comma = true;
          }
          ++pos;
        }
        ++ic;
      }

      if (pos == pos_nl) {
        // this row is missing the requested column,
        // so set to zero...
        doubleVec.push_back(0.0);
      }
      else {
        //otherwise parse as double...
        doubleVec.push_back(strtod(buf+pos,NULL));
      }

      // and advance pos to next line...
      pos = pos_nl+1;

    }

    close(fd);

    n = doubleVec.size();
    if (n > 0) {
      assert(data == NULL);
      data = new double[n];
      memcpy(data,&doubleVec.front(),sizeof(double)*n);
    }

    return 0;

  }

  int xcol2(double * &data1,double * &data2,int &n,char const *fname,const int column1,const int column2) {

    // for now, insist on increasing, but this can be easily fixed...
    assert(column1 < column2);

    int fd = open(fname, O_RDONLY);
    if (fd == -1) {
      cout << "Error: xcol file open problem for file: " << fname << endl;
      return -1;
    }

    // Advise the kernel of our access pattern: FDADVICE_SEQUENTIAL
    // note that this does not seem to be available on all systems tested...
    //posix_fadvise(fd, 0, 0, 1);

    //static const int BUFFER_SIZE = 16*8;
    static const int BUFFER_SIZE = 16*1024; // size used by "wc -l"...

    char buf[BUFFER_SIZE*2+1];
    size_t bytes_read = read(fd, buf, BUFFER_SIZE*2);
    if (bytes_read <= 0) {
      cout << "Error: xcol read returned bytes_read: " << bytes_read << " for file: " << fname << endl;
      close(fd);
      return -1;
    }

    // experimented with managing memory more directly, but file i/o is so much
    // slower than everything else, so just use a vector, because it is simplest...
    vector<double> doubleVec;

    // if the file is smaller than BUFFER_SIZE*2, then make sure it has '\n'
    // at the end of the last line and set pos_eof to the position after the '\n'...
    int pos_eof = -1;
    if (bytes_read < BUFFER_SIZE*2) {
      if (buf[bytes_read-1] != '\n') {
        buf[bytes_read++] = '\n';
      }
      pos_eof = bytes_read;
    }

    int pos = 0;
    int pos_nl = 0;
    while (pos != pos_eof) {

      // if we are in the back half of the buffer and we are not
      // done reading, consider advancing the buffer...
      if ((pos >= BUFFER_SIZE)&&(bytes_read == BUFFER_SIZE*2)) {
        memcpy(buf,buf+BUFFER_SIZE,BUFFER_SIZE);
        bytes_read = read(fd,buf+BUFFER_SIZE,BUFFER_SIZE);
        if (bytes_read < BUFFER_SIZE) {
          if (bytes_read < 0) {
            cerr << "Error: xcol bytes_read: " << bytes_read << endl;
            assert(0);
          }
          if (buf[BUFFER_SIZE+bytes_read-1] != '\n') {
            buf[BUFFER_SIZE+bytes_read] = '\n';
            ++bytes_read;
          }
          pos_eof = BUFFER_SIZE+bytes_read;
        }
        bytes_read += BUFFER_SIZE;
        pos -= BUFFER_SIZE;
      }

      // use memchr to find the next newline...
      const char * nl = (char*)memchr(buf+pos,'\n',bytes_read-pos);
      if (nl == NULL) {
        // if you hit this error, it means the buffer is not long enough to read the file
        // until a newline is found. This should produce an error?...
        cout << "Error: xcol file does not have a newline in " << bytes_read-pos << " chars. Increase BUFFER_SIZE in xcol." << endl;
        assert(0);
      }
      pos_nl = int(nl-buf);
      assert(pos <= pos_nl);

      // skip any line with first char a hash...
      if (buf[pos] == '#') {
        pos = pos_nl + 1;
        continue;
      }

      // clear leading whitespace once...
      while ((pos < pos_nl)&&((buf[pos]==' ')||(buf[pos]=='\t')))
        ++pos;

      // advance through the first columns to get to column1...
      int ic = 1;
      while (ic < column1) {
        // advance through non-whitespace...
        while ((pos < pos_nl)&&(buf[pos]!=' ')&&(buf[pos]!='\t')&&(buf[pos]!=','))
          ++pos;
        bool b_comma = false;
        // then advance through whitespace which may contain up to 1 comma...
        while ((pos < pos_nl)&&((buf[pos]==' ')||(buf[pos]=='\t')||(buf[pos]==','))) {
          if (buf[pos]==',') {
            // if we hit a second comma, then break without advancing pos...
            if (b_comma)
              break;
            b_comma = true;
          }
          ++pos;
        }
        ++ic;
      }

      if (pos == pos_nl) {
        // this row is missing the requested column,
        // so set to zero...
        doubleVec.push_back(0.0);
      }
      else {
        //otherwise parse as double...
        doubleVec.push_back(strtod(buf+pos,NULL));
      }

      while (ic < column2) {
        // advance through non-whitespace...
        while ((pos < pos_nl)&&(buf[pos]!=' ')&&(buf[pos]!='\t')&&(buf[pos]!=','))
          ++pos;
        bool b_comma = false;
        // then advance through whitespace which may contain up to 1 comma...
        while ((pos < pos_nl)&&((buf[pos]==' ')||(buf[pos]=='\t')||(buf[pos]==','))) {
          if (buf[pos]==',') {
            // if we hit a second comma, then break without advancing pos...
            if (b_comma)
              break;
            b_comma = true;
          }
          ++pos;
        }
        ++ic;
      }

      if (pos == pos_nl) {
        // this row is missing the requested column,
        // so set to zero...
        doubleVec.push_back(0.0);
      }
      else {
        //otherwise parse as double...
        doubleVec.push_back(strtod(buf+pos,NULL));
      }

      // and advance pos to next line...
      pos = pos_nl+1;

    }

    close(fd);

    n = doubleVec.size()/2;
    if (n > 0) {
      assert(data1 == NULL);
      data1 = new double[n];
      assert(data2 == NULL);
      data2 = new double[n];
      for (int i = 0; i < n; ++i) {
        data1[i] = doubleVec[2*i];
        data2[i] = doubleVec[2*i+1];
      }
    }

    return 0;

  }

  int xcol4(double * &data1,double * &data2,double * &data3,double * &data4,
	    int &n,char const *fname,
	    const int column1,const int column2,const int column3,const int column4) {

    // for now, insist on increasing, but this can be easily fixed...
    assert(column1 < column2);
    assert(column2 < column3);
    assert(column3 < column4);

    int fd = open(fname, O_RDONLY);
    if (fd == -1) {
      cout << "Error: xcol file open problem for file: " << fname << endl;
      return -1;
    }

    // Advise the kernel of our access pattern: FDADVICE_SEQUENTIAL
    // note that this does not seem to be available on all systems tested...
    //posix_fadvise(fd, 0, 0, 1);

    //static const int BUFFER_SIZE = 16*8;
    static const int BUFFER_SIZE = 16*1024; // size used by "wc -l"...

    char buf[BUFFER_SIZE*2+1];
    size_t bytes_read = read(fd, buf, BUFFER_SIZE*2);
    if (bytes_read <= 0) {
      cout << "Error: xcol read returned bytes_read: " << bytes_read << " for file: " << fname << endl;
      close(fd);
      return -1;
    }

    // experimented with managing memory more directly, but file i/o is so much
    // slower than everything else, so just use a vector, because it is simplest...
    vector<double> doubleVec;

    // if the file is smaller than BUFFER_SIZE*2, then make sure it has '\n'
    // at the end of the last line and set pos_eof to the position after the '\n'...
    int pos_eof = -1;
    if (bytes_read < BUFFER_SIZE*2) {
      if (buf[bytes_read-1] != '\n') {
        buf[bytes_read++] = '\n';
      }
      pos_eof = bytes_read;
    }

    int pos = 0;
    int pos_nl = 0;
    while (pos != pos_eof) {

      // if we are in the back half of the buffer and we are not
      // done reading, consider advancing the buffer...
      if ((pos >= BUFFER_SIZE)&&(bytes_read == BUFFER_SIZE*2)) {
        memcpy(buf,buf+BUFFER_SIZE,BUFFER_SIZE);
        bytes_read = read(fd,buf+BUFFER_SIZE,BUFFER_SIZE);
        if (bytes_read < BUFFER_SIZE) {
          if (bytes_read < 0) {
            cerr << "Error: xcol bytes_read: " << bytes_read << endl;
            assert(0);
          }
          if (buf[BUFFER_SIZE+bytes_read-1] != '\n') {
            buf[BUFFER_SIZE+bytes_read] = '\n';
            ++bytes_read;
          }
          pos_eof = BUFFER_SIZE+bytes_read;
        }
        bytes_read += BUFFER_SIZE;
        pos -= BUFFER_SIZE;
      }

      // use memchr to find the next newline...
      const char * nl = (char*)memchr(buf+pos,'\n',bytes_read-pos);
      if (nl == NULL) {
        // if you hit this error, it means the buffer is not long enough to read the file
        // until a newline is found. This should produce an error?...
        cout << "Error: xcol file does not have a newline in " << bytes_read-pos <<
	  " chars. Increase BUFFER_SIZE in xcol." << endl;
        assert(0);
      }
      pos_nl = int(nl-buf);
      assert(pos <= pos_nl);

      // skip any line with first char a hash...
      if (buf[pos] == '#') {
        pos = pos_nl + 1;
        continue;
      }

      // clear leading whitespace once...
      while ((pos < pos_nl)&&((buf[pos]==' ')||(buf[pos]=='\t')))
        ++pos;

      // advance through the first columns to get to column1...
      int ic = 1;
      while (ic < column1) {
        // advance through non-whitespace...
        while ((pos < pos_nl)&&(buf[pos]!=' ')&&(buf[pos]!='\t')&&(buf[pos]!=','))
          ++pos;
        bool b_comma = false;
        // then advance through whitespace which may contain up to 1 comma...
        while ((pos < pos_nl)&&((buf[pos]==' ')||(buf[pos]=='\t')||(buf[pos]==','))) {
          if (buf[pos]==',') {
            // if we hit a second comma, then break without advancing pos...
            if (b_comma)
              break;
            b_comma = true;
          }
          ++pos;
        }
        ++ic;
      }

      if (pos == pos_nl) {
        // this row is missing the requested column,
        // so set to zero...
        doubleVec.push_back(0.0);
      }
      else {
        //otherwise parse as double...
        doubleVec.push_back(strtod(buf+pos,NULL));
      }

      while (ic < column2) {
        // advance through non-whitespace...
        while ((pos < pos_nl)&&(buf[pos]!=' ')&&(buf[pos]!='\t')&&(buf[pos]!=','))
          ++pos;
        bool b_comma = false;
        // then advance through whitespace which may contain up to 1 comma...
        while ((pos < pos_nl)&&((buf[pos]==' ')||(buf[pos]=='\t')||(buf[pos]==','))) {
          if (buf[pos]==',') {
            // if we hit a second comma, then break without advancing pos...
            if (b_comma)
              break;
            b_comma = true;
          }
          ++pos;
        }
        ++ic;
      }

      if (pos == pos_nl) {
        // this row is missing the requested column,
        // so set to zero...
        doubleVec.push_back(0.0);
      }
      else {
        //otherwise parse as double...
        doubleVec.push_back(strtod(buf+pos,NULL));
      }

      while (ic < column3) {
        // advance through non-whitespace...
        while ((pos < pos_nl)&&(buf[pos]!=' ')&&(buf[pos]!='\t')&&(buf[pos]!=','))
          ++pos;
        bool b_comma = false;
        // then advance through whitespace which may contain up to 1 comma...
        while ((pos < pos_nl)&&((buf[pos]==' ')||(buf[pos]=='\t')||(buf[pos]==','))) {
          if (buf[pos]==',') {
            // if we hit a second comma, then break without advancing pos...
            if (b_comma)
              break;
            b_comma = true;
          }
          ++pos;
        }
        ++ic;
      }

      if (pos == pos_nl) {
        // this row is missing the requested column,
        // so set to zero...
        doubleVec.push_back(0.0);
      }
      else {
        //otherwise parse as double...
        doubleVec.push_back(strtod(buf+pos,NULL));
      }

      while (ic < column4) {
        // advance through non-whitespace...
        while ((pos < pos_nl)&&(buf[pos]!=' ')&&(buf[pos]!='\t')&&(buf[pos]!=','))
          ++pos;
        bool b_comma = false;
        // then advance through whitespace which may contain up to 1 comma...
        while ((pos < pos_nl)&&((buf[pos]==' ')||(buf[pos]=='\t')||(buf[pos]==','))) {
          if (buf[pos]==',') {
            // if we hit a second comma, then break without advancing pos...
            if (b_comma)
              break;
            b_comma = true;
          }
          ++pos;
        }
        ++ic;
      }

      if (pos == pos_nl) {
        // this row is missing the requested column,
        // so set to zero...
        doubleVec.push_back(0.0);
      }
      else {
        //otherwise parse as double...
        doubleVec.push_back(strtod(buf+pos,NULL));
      }

      // and advance pos to next line...
      pos = pos_nl+1;

    }

    close(fd);

    n = doubleVec.size()/4;
    if (n > 0) {
      assert(data1 == NULL);
      data1 = new double[n];
      assert(data2 == NULL);
      data2 = new double[n];
      assert(data3 == NULL);
      data3 = new double[n];
      assert(data4 == NULL);
      data4 = new double[n];
      for (int i = 0; i < n; ++i) {
        data1[i] = doubleVec[4*i];
        data2[i] = doubleVec[4*i+1];
        data3[i] = doubleVec[4*i+2];
        data4[i] = doubleVec[4*i+3];
      }
    }

    return 0;

  }

  int grepxcol(double * &data,int &n,char const *fname,char const *expr,const int column) {

    int fd = open(fname, O_RDONLY);
    if (fd == -1) {
      cout << "Error: grepxcol file open problem for file: " << fname << endl;
      return -1;
    }

    // Advise the kernel of our access pattern: FDADVICE_SEQUENTIAL
    // note that this does not seem to be available on all systems tested...
    //posix_fadvise(fd, 0, 0, 1);

    //static const int BUFFER_SIZE = 16*8;
    static const int BUFFER_SIZE = 16*1024; // size used by "wc -l"...

    char buf[BUFFER_SIZE*2+1];
    size_t bytes_read = read(fd, buf, BUFFER_SIZE*2);
    if (bytes_read <= 0) {
      cout << "Error: grepxcol read returned bytes_read: " << bytes_read << " for file: " << fname << endl;
      close(fd);
      return -1;
    }

    const int nexpr = strlen(expr);

    // experimented with managing memory more directly, but file i/o is so much
    // slower than everything else, so just use a vector, because it is simplest...
    vector<double> doubleVec;

    // if the file is smaller than BUFFER_SIZE*2, then make sure it has '\n'
    // at the end of the last line and set pos_eof to the position after the '\n'...
    int pos_eof = -1;
    if (bytes_read < BUFFER_SIZE*2) {
      if (buf[bytes_read-1] != '\n') {
        buf[bytes_read++] = '\n';
      }
      pos_eof = bytes_read;
    }

    int pos = 0;
    int pos_nl = 0;
    while (pos != pos_eof) {

      // if we are in the back half of the buffer and we are not
      // done reading, consider advancing the buffer...
      if ((pos >= BUFFER_SIZE)&&(bytes_read == BUFFER_SIZE*2)) {
        memcpy(buf,buf+BUFFER_SIZE,BUFFER_SIZE);
        bytes_read = read(fd,buf+BUFFER_SIZE,BUFFER_SIZE);
        if (bytes_read < BUFFER_SIZE) {
          if (bytes_read < 0) {
            cerr << "Error: grepxcol bytes_read: " << bytes_read << endl;
            assert(0);
          }
          if (buf[BUFFER_SIZE+bytes_read-1] != '\n') {
            buf[BUFFER_SIZE+bytes_read] = '\n';
            ++bytes_read;
          }
          pos_eof = BUFFER_SIZE+bytes_read;
        }
        bytes_read += BUFFER_SIZE;
        pos -= BUFFER_SIZE;
      }

      // use memchr to find the next newline...
      const char * nl = (char*)memchr(buf+pos,'\n',bytes_read-pos);
      if (nl == NULL) {
        // if you hit this error, it means the buffer is not long enough to read the file
        // until a newline is found. This should produce an error?...
        cout << "Error: grepxcol file does not have a newline in " << bytes_read-pos << " chars. Increase BUFFER_SIZE in xcol." << endl;
        assert(0);
      }
      pos_nl = int(nl-buf);

      // no point in checking for anything in a line shorter than "expr"...
      if (pos_nl-pos < nexpr) {
        pos = pos_nl + 1;
        continue;
      }

      // general grep on expression...
      char * p = buf+pos;
      while ((p = (char*)memchr(p,expr[0],(buf+pos_nl)-p-nexpr+1))) {
        if (memcmp(p,expr,nexpr) == 0)
          break;
        ++p;
        assert((buf+pos_nl)-p-nexpr+1 >= 0);
      }
      if (!p) {
        pos = pos_nl + 1;
        continue;
      }

      // clear leading whitespace once...
      while ((pos < pos_nl)&&((buf[pos]==' ')||(buf[pos]=='\t')))
        ++pos;

      // advance through the first columns to get to column1...
      int ic = 1;
      while (ic < column) {
        // advance through non-whitespace...
        while ((pos < pos_nl)&&(buf[pos]!=' ')&&(buf[pos]!='\t')&&(buf[pos]!=','))
          ++pos;
        bool b_comma = false;
        // then advance through whitespace which may contain up to 1 comma...
        while ((pos < pos_nl)&&((buf[pos]==' ')||(buf[pos]=='\t')||(buf[pos]==','))) {
          if (buf[pos]==',') {
            // if we hit a second comma, then break without advancing pos...
            if (b_comma)
              break;
            b_comma = true;
          }
          ++pos;
        }
        ++ic;
      }

      if (pos == pos_nl) {
        // this row is missing the requested column,
        // so set to zero...
        doubleVec.push_back(0.0);
      }
      else {
        // otherwise parse double...
        doubleVec.push_back(strtod(buf+pos,NULL));
      }

      // and advance pos to next line...
      pos = pos_nl+1;

    }

    close(fd);

    n = doubleVec.size();
    if (n > 0) {
      assert(data == NULL);
      data = new double[n];
      memcpy(data,&doubleVec.front(),sizeof(double)*n);
    }

    return 0;

  }

  int grepxcol2(double * &data1,double * &data2,int &n,char const *fname,char const *expr,const int column1,const int column2) {

    // for now, insist on increasing, but this can be easily fixed...
    assert(column1 < column2);

    int fd = open(fname, O_RDONLY);
    if (fd == -1) {
      cout << "Error: grepxcol file open problem for file: " << fname << endl;
      return -1;
    }

    // Advise the kernel of our access pattern: FDADVICE_SEQUENTIAL
    // note that this does not seem to be available on all systems tested...
    //posix_fadvise(fd, 0, 0, 1);

    //static const int BUFFER_SIZE = 16*8;
    static const int BUFFER_SIZE = 16*1024; // size used by "wc -l"...

    char buf[BUFFER_SIZE*2+1];
    size_t bytes_read = read(fd, buf, BUFFER_SIZE*2);
    if (bytes_read <= 0) {
      cout << "Error: grepxcol read returned bytes_read: " << bytes_read << " for file: " << fname << endl;
      close(fd);
      return -1;
    }

    const int nexpr = strlen(expr);

    // experimented with managing memory more directly, but file i/o is so much
    // slower than everything else, so just use a vector, because it is simplest...
    vector<double> doubleVec;

    // if the file is smaller than BUFFER_SIZE*2, then make sure it has '\n'
    // at the end of the last line and set pos_eof to the position after the '\n'...
    int pos_eof = -1;
    if (bytes_read < BUFFER_SIZE*2) {
      if (buf[bytes_read-1] != '\n') {
        buf[bytes_read++] = '\n';
      }
      pos_eof = bytes_read;
    }

    int pos = 0;
    int pos_nl = 0;
    while (pos != pos_eof) {

      // if we are in the back half of the buffer and we are not
      // done reading, consider advancing the buffer...
      if ((pos >= BUFFER_SIZE)&&(bytes_read == BUFFER_SIZE*2)) {
        memcpy(buf,buf+BUFFER_SIZE,BUFFER_SIZE);
        bytes_read = read(fd,buf+BUFFER_SIZE,BUFFER_SIZE);
        if (bytes_read < BUFFER_SIZE) {
          if (bytes_read < 0) {
            cerr << "Error: grepxcol bytes_read: " << bytes_read << endl;
            assert(0);
          }
          if (buf[BUFFER_SIZE+bytes_read-1] != '\n') {
            buf[BUFFER_SIZE+bytes_read] = '\n';
            ++bytes_read;
          }
          pos_eof = BUFFER_SIZE+bytes_read;
        }
        bytes_read += BUFFER_SIZE;
        pos -= BUFFER_SIZE;
      }

      // use memchr to find the next newline...
      const char * nl = (char*)memchr(buf+pos,'\n',bytes_read-pos);
      if (nl == NULL) {
        // if you hit this error, it means the buffer is not long enough to read the file
        // until a newline is found. This should produce an error?...
        cout << "Error: grepxcol file does not have a newline in " << bytes_read-pos << " chars. Increase BUFFER_SIZE in xcol." << endl;
        close(fd);
        return -1;
      }
      pos_nl = int(nl-buf);

      // no point in checking for anything in a line shorter than "expr"...
      if (pos_nl-pos < nexpr) {
        pos = pos_nl + 1;
        continue;
      }

      // general grep on expression...
      char * p = buf+pos;
      while ((p = (char*)memchr(p,expr[0],(buf+pos_nl)-p-nexpr+1))) {
        if (memcmp(p,expr,nexpr) == 0)
          break;
        ++p;
        assert((buf+pos_nl)-p-nexpr+1 >= 0);
      }
      if (!p) {
        pos = pos_nl + 1;
        continue;
      }

      // clear leading whitespace once...
      while ((pos < pos_nl)&&((buf[pos]==' ')||(buf[pos]=='\t')))
        ++pos;

      // advance through the first columns to get to column1...
      int ic = 1;
      while (ic < column1) {
        // advance through non-whitespace...
        while ((pos < pos_nl)&&(buf[pos]!=' ')&&(buf[pos]!='\t')&&(buf[pos]!=','))
          ++pos;
        bool b_comma = false;
        // then advance through whitespace which may contain up to 1 comma...
        while ((pos < pos_nl)&&((buf[pos]==' ')||(buf[pos]=='\t')||(buf[pos]==','))) {
          if (buf[pos]==',') {
            // if we hit a second comma, then break without advancing pos...
            if (b_comma)
              break;
            b_comma = true;
          }
          ++pos;
        }
        ++ic;
      }

      if (pos == pos_nl) {
        // this row is missing the requested column,
        // so set to zero...
        doubleVec.push_back(0.0);
      }
      else {
        // otherwise parse double...
        doubleVec.push_back(strtod(buf+pos,NULL));
      }

      // note: We do not change the parsing pos with 2nd arg of strtod because its returned position
      // does not necessaily obey the correct definition of whitespace separation. For example,
      // strtod will process the character sequence "17-D" as a double 17, returning the
      // end as the "-" char.

      while (ic < column2) {
        // advance through non-whitespace...
        while ((pos < pos_nl)&&(buf[pos]!=' ')&&(buf[pos]!='\t')&&(buf[pos]!=','))
          ++pos;
        bool b_comma = false;
        // then advance through whitespace which may contain up to 1 comma...
        while ((pos < pos_nl)&&((buf[pos]==' ')||(buf[pos]=='\t')||(buf[pos]==','))) {
          if (buf[pos]==',') {
            // if we hit a second comma, then break without advancing pos...
            if (b_comma)
              break;
            b_comma = true;
          }
          ++pos;
        }
        ++ic;
      }

      if (pos == pos_nl) {
        // this row is missing the requested column,
        // so set to zero...
        doubleVec.push_back(0.0);
      }
      else {
        // otherwise parse double...
        doubleVec.push_back(strtod(buf+pos,NULL));
      }

      // and advance pos to next line...
      pos = pos_nl+1;

    }

    close(fd);

    n = doubleVec.size()/2;
    if (n > 0) {
      assert(data1 == NULL);
      data1 = new double[n];
      assert(data2 == NULL);
      data2 = new double[n];
      for (int i = 0; i < n; ++i) {
        data1[i] = doubleVec[2*i];
        data2[i] = doubleVec[2*i+1];
      }
    }

    return 0;

  }

  int tailxcol(double * &data,int &n,char const *fname,const int column) {

    // n is input/output
    if (n <= 0) {
      cout << "Error: tailxcol requires n>0: " << n << endl;
      return -1;
    }

    int fd = open(fname, O_RDONLY);
    if (fd == -1) {
      cout << "Error: tailxcol file open problem for file: " << fname << endl;
      return -1;
    }

    struct stat stat_buf;
    if (fstat(fd, &stat_buf) != 0) {
      cout << "Error: tailxcol fstat problem for file: " << fname << endl;
      close(fd);
      return -1;
    }

    //static const int BUFFER_SIZE = 16*8;
    static const int BUFFER_SIZE = 16*1024; // size used by "wc -l"...

    char buf[BUFFER_SIZE*2+1];
    size_t bytes_read = read(fd, buf, BUFFER_SIZE*2);
    if (bytes_read <= 0) {
      cout << "Error: tailxcol read returned bytes_read: " << bytes_read << " for file: " << fname << endl;
      close(fd);
      return -1;
    }

    // experimented with managing memory more directly, but file i/o is so much
    // slower than everything else, so just use a vector, because it is simplest...
    vector<double> doubleVec;

    // if the file is smaller than BUFFER_SIZE*2, then make sure it has '\n'
    // at the end of the last line and set pos_eof to the position after the '\n'...
    int pos_eof = -1;
    if (bytes_read < BUFFER_SIZE*2) {
      if (buf[bytes_read-1] != '\n') {
        buf[bytes_read++] = '\n';
      }
      pos_eof = bytes_read;
    }

    // for prediction, read in 256 lines and record infor at 128 and 256. If we don't get to
    // 256, then we have simply read everything, and we adjust at the bottom ot htis routine...
    const int n1 = max(128,n/4);
    const int n2 = max(256,n/2);
    bool b_init = true;

    size_t offset = 0,offset1;
    int pos = 0;
    int pos_nl = 0;
    while (pos != pos_eof) {

      // if we are in the back half of the buffer and we are not
      // done reading, consider advancing the buffer...
      if ((pos >= BUFFER_SIZE)&&(bytes_read == BUFFER_SIZE*2)) {
        memcpy(buf,buf+BUFFER_SIZE,BUFFER_SIZE);
        offset += BUFFER_SIZE;
        bytes_read = read(fd,buf+BUFFER_SIZE,BUFFER_SIZE);
        if (bytes_read < BUFFER_SIZE) {
          if (bytes_read < 0) {
            cerr << "Error: tailxcol bytes_read: " << bytes_read << endl;
            assert(0);
          }
          if (buf[BUFFER_SIZE+bytes_read-1] != '\n') {
            buf[BUFFER_SIZE+bytes_read] = '\n';
            ++bytes_read;
          }
          pos_eof = BUFFER_SIZE+bytes_read;
        }
        bytes_read += BUFFER_SIZE;
        pos -= BUFFER_SIZE;
      }

      // use memchr to find the next newline...
      const char * nl = (char*)memchr(buf+pos,'\n',bytes_read-pos);
      if (nl == NULL) {
        // if you hit this error, it means the buffer is not long enough to read the file
        // until a newline is found. This should produce an error?...
        cout << "Error: tailxcol file does not have a newline in " << bytes_read-pos << " chars. Increase BUFFER_SIZE in xcol." << endl;
        assert(0);
      }
      pos_nl = int(nl-buf);

      // skip any line with first char a hash...
      if (buf[pos] == '#') {
        pos = pos_nl + 1;
        continue;
      }

      // clear leading whitespace once...
      while ((pos < pos_nl)&&((buf[pos]==' ')||(buf[pos]=='\t')))
        ++pos;

      // advance through the first columns to get to column1...
      int ic = 1;
      while (ic < column) {
        // advance through non-whitespace...
        while ((pos < pos_nl)&&(buf[pos]!=' ')&&(buf[pos]!='\t')&&(buf[pos]!=','))
          ++pos;
        bool b_comma = false;
        // then advance through whitespace which may contain up to 1 comma...
        while ((pos < pos_nl)&&((buf[pos]==' ')||(buf[pos]=='\t')||(buf[pos]==','))) {
          if (buf[pos]==',') {
            // if we hit a second comma, then break without advancing pos...
            if (b_comma)
              break;
            b_comma = true;
          }
          ++pos;
        }
        ++ic;
      }

      if (pos == pos_nl) {
        // this row is missing the requested column,
        // so set to zero...
        doubleVec.push_back(0.0);
      }
      else {
        // otherwise parse double...
        doubleVec.push_back(strtod(buf+pos,NULL));
      }

      // and advance pos to next line...
      pos = pos_nl+1;

      if (b_init) {
        if (doubleVec.size() == n1) {
          offset1 = offset+pos_nl;
        }
        else if (doubleVec.size() == n2) {
          // when we hit n2, do the prediction...
          const size_t tail_size = (n+max(64,n/8))*((offset+pos_nl-offset1)/(n2-n1)+1);
          if (tail_size+offset+bytes_read < stat_buf.st_size) {
            // the tail_size is less than the remaining part of the file, so advance and
            // read just the tail_size...
            size_t target = stat_buf.st_size-tail_size;
            lseek(fd,target,SEEK_SET);
            doubleVec.clear();
            assert(pos_eof == -1);
            bytes_read = read(fd, buf, BUFFER_SIZE*2);
            if (bytes_read < BUFFER_SIZE*2) {
              if (buf[bytes_read-1] != '\n') {
                buf[bytes_read++] = '\n';
              }
              pos_eof = bytes_read;
            }
            const char * nl = (char*)memchr(buf+pos,'\n',bytes_read-pos);
            assert(nl != NULL);
            pos = (nl-buf)+1; // advances to the start of the next line...
          }
          b_init = false;
        }
      }

    }

    close(fd);

    if (doubleVec.size() < n) {
      if (!b_init) {
        cout << "possible prediction problem: n requested: " << n << ", n tailed: " << doubleVec.size() << endl;
      }
      n = doubleVec.size();
      if (data == NULL) data = new double[n];
      memcpy(data,&doubleVec.front(),sizeof(double)*n);
    }
    else {
      const int dn = doubleVec.size() - n;
      if (data == NULL) data = new double[n];
      memcpy(data,&doubleVec.front()+dn,sizeof(double)*n);
    }

    return 0;

  }

  int tailxcol2(double * &data1,double * &data2,int &n,char const *fname,const int column1,const int column2) {

    // n is input/output
    if (n <= 0) {
      cout << "Error: tailxcol requires n>0: " << n << endl;
      return -1;
    }

    // for now, insist on increasing, but this can be easily fixed...
    assert(column1 < column2);

    int fd = open(fname, O_RDONLY);
    if (fd == -1) {
      cout << "Error: tailxcol file open problem for file: " << fname << endl;
      return -1;
    }

    struct stat stat_buf;
    if (fstat(fd, &stat_buf) != 0) {
      cout << "Error: tailxcol fstat problem for file: " << fname << endl;
      close(fd);
      return -1;
    }

    //static const int BUFFER_SIZE = 16*8;
    static const int BUFFER_SIZE = 16*1024; // size used by "wc -l"...

    char buf[BUFFER_SIZE*2+1];
    size_t bytes_read = read(fd, buf, BUFFER_SIZE*2);
    if (bytes_read <= 0) {
      cout << "Error: tailxcol read returned bytes_read: " << bytes_read << " for file: " << fname << endl;
      close(fd);
      return -1;
    }

    // experimented with managing memory more directly, but file i/o is so much
    // slower than everything else, so just use a vector, because it is simplest...
    vector<double> doubleVec;

    // if the file is smaller than BUFFER_SIZE*2, then make sure it has '\n'
    // at the end of the last line and set pos_eof to the position after the '\n'...
    int pos_eof = -1;
    if (bytes_read < BUFFER_SIZE*2) {
      if (buf[bytes_read-1] != '\n') {
        buf[bytes_read++] = '\n';
      }
      pos_eof = bytes_read;
    }

    // for prediction, read in 256 lines and record infor at 128 and 256. If we don't get to
    // 256, then we have simply read everything, and we adjust at the bottom ot htis routine...
    const int n1 = max(128,n/4);
    const int n2 = max(256,n/2);
    bool b_init = true;

    size_t offset = 0,offset1;
    int pos = 0;
    int pos_nl = 0;
    while (pos != pos_eof) {

      // if we are in the back half of the buffer and we are not
      // done reading, consider advancing the buffer...
      if ((pos >= BUFFER_SIZE)&&(bytes_read == BUFFER_SIZE*2)) {
        memcpy(buf,buf+BUFFER_SIZE,BUFFER_SIZE);
        offset += BUFFER_SIZE;
        bytes_read = read(fd,buf+BUFFER_SIZE,BUFFER_SIZE);
        if (bytes_read < BUFFER_SIZE) {
          if (bytes_read < 0) {
            cerr << "Error: tailxcol bytes_read: " << bytes_read << endl;
            assert(0);
          }
          if (buf[BUFFER_SIZE+bytes_read-1] != '\n') {
            buf[BUFFER_SIZE+bytes_read] = '\n';
            ++bytes_read;
          }
          pos_eof = BUFFER_SIZE+bytes_read;
        }
        bytes_read += BUFFER_SIZE;
        pos -= BUFFER_SIZE;
      }

      // use memchr to find the next newline...
      const char * nl = (char*)memchr(buf+pos,'\n',bytes_read-pos);
      if (nl == NULL) {
        // if you hit this error, it means the buffer is not long enough to read the file
        // until a newline is found. This should produce an error?...
        cout << "Error: tailxcol file does not have a newline in " << bytes_read-pos << " chars. Increase BUFFER_SIZE in xcol." << endl;
        assert(0);
      }
      pos_nl = int(nl-buf);

      // skip any line with first char a hash...
      if (buf[pos] == '#') {
        pos = pos_nl + 1;
        continue;
      }

      // clear leading whitespace once...
      while ((pos < pos_nl)&&((buf[pos]==' ')||(buf[pos]=='\t')))
        ++pos;

      // advance through the first columns to get to column1...
      int ic = 1;
      while (ic < column1) {
        // advance through non-whitespace...
        while ((pos < pos_nl)&&(buf[pos]!=' ')&&(buf[pos]!='\t')&&(buf[pos]!=','))
          ++pos;
        bool b_comma = false;
        // then advance through whitespace which may contain up to 1 comma...
        while ((pos < pos_nl)&&((buf[pos]==' ')||(buf[pos]=='\t')||(buf[pos]==','))) {
          if (buf[pos]==',') {
            // if we hit a second comma, then break without advancing pos...
            if (b_comma)
              break;
            b_comma = true;
          }
          ++pos;
        }
        ++ic;
      }

      if (pos == pos_nl) {
        // this row is missing the requested column,
        // so set to zero...
        doubleVec.push_back(0.0);
      }
      else {
        // otherwise parse double...
        doubleVec.push_back(strtod(buf+pos,NULL));
      }

      while (ic < column2) {
        // advance through non-whitespace...
        while ((pos < pos_nl)&&(buf[pos]!=' ')&&(buf[pos]!='\t')&&(buf[pos]!=','))
          ++pos;
        bool b_comma = false;
        // then advance through whitespace which may contain up to 1 comma...
        while ((pos < pos_nl)&&((buf[pos]==' ')||(buf[pos]=='\t')||(buf[pos]==','))) {
          if (buf[pos]==',') {
            // if we hit a second comma, then break without advancing pos...
            if (b_comma)
              break;
            b_comma = true;
          }
          ++pos;
        }
        ++ic;
      }

      if (pos == pos_nl) {
        // this row is missing the requested column,
        // so set to zero...
        doubleVec.push_back(0.0);
      }
      else {
        // otherwise parse double...
        doubleVec.push_back(strtod(buf+pos,NULL));
      }

      // and advance pos to next line...
      pos = pos_nl+1;

      if (b_init) {
        if (doubleVec.size() == n1*2) {
          offset1 = offset+pos_nl;
        }
        else if (doubleVec.size() == n2*2) {
          // when we hit n2, do the prediction...
          const size_t tail_size = (n+max(64,n/8))*((offset+pos_nl-offset1)/(n2-n1)+1);
          if (tail_size+offset+bytes_read < stat_buf.st_size) {
            // the tail_size is less than the remaining part of the file, so advance and
            // read just the tail_size...
            size_t target = stat_buf.st_size-tail_size;
            lseek(fd,target,SEEK_SET);
            doubleVec.clear();
            assert(pos_eof == -1);
            bytes_read = read(fd, buf, BUFFER_SIZE*2);
            if (bytes_read < BUFFER_SIZE*2) {
              if (buf[bytes_read-1] != '\n') {
                buf[bytes_read++] = '\n';
              }
              pos_eof = bytes_read;
            }
            const char * nl = (char*)memchr(buf+pos,'\n',bytes_read-pos);
            assert(nl != NULL);
            pos = (nl-buf)+1; // advances to the start of the next line...
          }
          b_init = false;
        }
      }

    }

    close(fd);

    if (doubleVec.size() < n*2) {
      if (!b_init) {
        cout << "possible prediction problem: n requested: " << n << ", n tailed: " << doubleVec.size()/2 << endl;
      }
      n = doubleVec.size()/2;
      if (data1 == NULL) data1 = new double[n];
      if (data2 == NULL) data2 = new double[n];
      for (int i = 0; i < n; ++i) {
        data1[i] = doubleVec[2*i];
        data2[i] = doubleVec[2*i+1];
      }
    }
    else {
      const int dn2 = doubleVec.size() - n*2;
      if (data1 == NULL) data1 = new double[n];
      if (data2 == NULL) data2 = new double[n];
      for (int i = 0; i < n; ++i) {
        data1[i] = doubleVec[dn2+2*i];
        data2[i] = doubleVec[dn2+2*i+1];
      }
    }

    return 0;

  }

  int tailgrepxcol(double * &data,int &n,char const *fname,char const *expr,const int column) {

    // n is input/output
    if (n <= 0) {
      cout << "Error: tailgrepxcol requires n>0: " << n << endl;
      return -1;
    }

    int fd = open(fname, O_RDONLY);
    if (fd == -1) {
      cout << "Error: tailgrepxcol file open problem for file: " << fname << endl;
      return -1;
    }

    struct stat stat_buf;
    if (fstat(fd, &stat_buf) != 0) {
      cout << "Error: tailgrepxcol fstat problem for file: " << fname << endl;
      close(fd);
      return -1;
    }

    //static const int BUFFER_SIZE = 16*8;
    static const int BUFFER_SIZE = 16*1024; // size used by "wc -l"...

    char buf[BUFFER_SIZE*2+1];
    size_t bytes_read = read(fd, buf, BUFFER_SIZE*2);
    if (bytes_read <= 0) {
      cout << "Error: tailgrepxcol read returned bytes_read: " << bytes_read << " for file: " << fname << endl;
      close(fd);
      return -1;
    }

    const int nexpr = strlen(expr);

    // experimented with managing memory more directly, but file i/o is so much
    // slower than everything else, so just use a vector, because it is simplest...
    vector<double> doubleVec;

    // if the file is smaller than BUFFER_SIZE*2, then make sure it has '\n'
    // at the end of the last line and set pos_eof to the position after the '\n'...
    int pos_eof = -1;
    if (bytes_read < BUFFER_SIZE*2) {
      if (buf[bytes_read-1] != '\n') {
        buf[bytes_read++] = '\n';
      }
      pos_eof = bytes_read;
    }

    // for prediction, read in 256 lines and record infor at 128 and 256. If we don't get to
    // 256, then we have simply read everything, and we adjust at the bottom ot htis routine...
    const int n1 = max(128,n/4);
    const int n2 = max(256,n/2);
    bool b_init = true;

    size_t offset = 0,offset1;
    int pos = 0;
    int pos_nl = 0;
    while (pos != pos_eof) {

      // if we are in the back half of the buffer and we are not
      // done reading, consider advancing the buffer...
      if ((pos >= BUFFER_SIZE)&&(bytes_read == BUFFER_SIZE*2)) {
        memcpy(buf,buf+BUFFER_SIZE,BUFFER_SIZE);
        offset += BUFFER_SIZE;
        bytes_read = read(fd,buf+BUFFER_SIZE,BUFFER_SIZE);
        if (bytes_read < BUFFER_SIZE) {
          if (bytes_read < 0) {
            cerr << "Error: tailgrepxcol bytes_read: " << bytes_read << endl;
            assert(0);
          }
          if (buf[BUFFER_SIZE+bytes_read-1] != '\n') {
            buf[BUFFER_SIZE+bytes_read] = '\n';
            ++bytes_read;
          }
          pos_eof = BUFFER_SIZE+bytes_read;
        }
        bytes_read += BUFFER_SIZE;
        pos -= BUFFER_SIZE;
      }

      // use memchr to find the next newline...
      const char * nl = (char*)memchr(buf+pos,'\n',bytes_read-pos);
      if (nl == NULL) {
        // if you hit this error, it means the buffer is not long enough to read the file
        // until a newline is found. This should produce an error?...
        cout << "Error: tailgrepxcol file does not have a newline in " << bytes_read-pos << " chars. Increase BUFFER_SIZE in xcol." << endl;
        assert(0);
      }
      pos_nl = int(nl-buf);

      // no point in checking for anything in a line shorter than "expr"...
      if (pos_nl-pos < nexpr) {
        pos = pos_nl + 1;
        continue;
      }

      // general grep on expression...
      char * p = buf+pos;
      while ((p = (char*)memchr(p,expr[0],(buf+pos_nl)-p-nexpr+1))) {
        if (memcmp(p,expr,nexpr) == 0)
          break;
        ++p;
        assert((buf+pos_nl)-p-nexpr+1 >= 0);
      }
      if (!p) {
        pos = pos_nl + 1;
        continue;
      }

      // clear leading whitespace once...
      while ((pos < pos_nl)&&((buf[pos]==' ')||(buf[pos]=='\t')))
        ++pos;

      // advance through the first columns to get to column1...
      int ic = 1;
      while (ic < column) {
        // advance through non-whitespace...
        while ((pos < pos_nl)&&(buf[pos]!=' ')&&(buf[pos]!='\t')&&(buf[pos]!=','))
          ++pos;
        bool b_comma = false;
        // then advance through whitespace which may contain up to 1 comma...
        while ((pos < pos_nl)&&((buf[pos]==' ')||(buf[pos]=='\t')||(buf[pos]==','))) {
          if (buf[pos]==',') {
            // if we hit a second comma, then break without advancing pos...
            if (b_comma)
              break;
            b_comma = true;
          }
          ++pos;
        }
        ++ic;
      }

      if (pos == pos_nl) {
        // this row is missing the requested column,
        // so set to zero...
        doubleVec.push_back(0.0);
      }
      else {
        // otherwise parse double...
        doubleVec.push_back(strtod(buf+pos,NULL));
      }

      // and advance pos to next line...
      pos = pos_nl+1;

      if (b_init) {
        if (doubleVec.size() == n1) {
          offset1 = offset+pos_nl;
        }
        else if (doubleVec.size() == n2) {
          // when we hit n2, do the prediction...
          const size_t tail_size = (n+max(64,n/8))*((offset+pos_nl-offset1)/(n2-n1)+1);
          if (tail_size+offset+bytes_read < stat_buf.st_size) {
            // the tail_size is less than the remaining part of the file, so advance and
            // read just the tail_size...
            size_t target = stat_buf.st_size-tail_size;
            lseek(fd,target,SEEK_SET);
            doubleVec.clear();
            assert(pos_eof == -1);
            bytes_read = read(fd, buf, BUFFER_SIZE*2);
            if (bytes_read < BUFFER_SIZE*2) {
              if (buf[bytes_read-1] != '\n') {
                buf[bytes_read++] = '\n';
              }
              pos_eof = bytes_read;
            }
            const char * nl = (char*)memchr(buf+pos,'\n',bytes_read-pos);
            assert(nl != NULL);
            pos = (nl-buf)+1; // advances to the start of the next line...
          }
          b_init = false;
        }
      }

    }

    close(fd);

    if (doubleVec.size() < n) {
      if (!b_init) {
        cout << "possible prediction problem: n requested: " << n << ", n tailed: " << doubleVec.size() << endl;
      }
      n = doubleVec.size();
      if (data == NULL) data = new double[n];
      memcpy(data,&doubleVec.front(),sizeof(double)*n);
    }
    else {
      const int dn = doubleVec.size() - n;
      if (data == NULL) data = new double[n];
      memcpy(data,&doubleVec.front()+dn,sizeof(double)*n);
    }

    return 0;

  }

  int tailgrepxcol2(double * &data1,double * &data2,int &n,char const *fname,char const *expr,const int column1,const int column2) {

    // n is input/output
    if (n <= 0) {
      cout << "Error: tailgrepxcol requires n>0: " << n << endl;
      return -1;
    }

    // for now, insist on increasing, but this can be easily fixed...
    assert(column1 < column2);

    int fd = open(fname, O_RDONLY);
    if (fd == -1) {
      cout << "Error: tailgrepxcol file open problem for file: " << fname << endl;
      return -1;
    }

    struct stat stat_buf;
    if (fstat(fd, &stat_buf) != 0) {
      cout << "Error: tailgrepxcol fstat problem for file: " << fname << endl;
      close(fd);
      return -1;
    }

    //static const int BUFFER_SIZE = 16*8;
    static const int BUFFER_SIZE = 16*1024; // size used by "wc -l"...

    char buf[BUFFER_SIZE*2+1];
    size_t bytes_read = read(fd, buf, BUFFER_SIZE*2);
    if (bytes_read <= 0) {
      cout << "Error: tailgrepxcol read returned bytes_read: " << bytes_read << " for file: " << fname << endl;
      close(fd);
      return -1;
    }

    const int nexpr = strlen(expr);

    // experimented with managing memory more directly, but file i/o is so much
    // slower than everything else, so just use a vector, because it is simplest...
    vector<double> doubleVec;

    // if the file is smaller than BUFFER_SIZE*2, then make sure it has '\n'
    // at the end of the last line and set pos_eof to the position after the '\n'...
    int pos_eof = -1;
    if (bytes_read < BUFFER_SIZE*2) {
      if (buf[bytes_read-1] != '\n') {
        buf[bytes_read++] = '\n';
      }
      pos_eof = bytes_read;
    }

    // for prediction, read in 256 lines and record infor at 128 and 256. If we don't get to
    // 256, then we have simply read everything, and we adjust at the bottom ot htis routine...
    const int n1 = max(128,n/4);
    const int n2 = max(256,n/2);
    bool b_init = true;

    size_t offset = 0,offset1;
    int pos = 0;
    int pos_nl = 0;
    while (pos != pos_eof) {

      // if we are in the back half of the buffer and we are not
      // done reading, consider advancing the buffer...
      if ((pos >= BUFFER_SIZE)&&(bytes_read == BUFFER_SIZE*2)) {
        memcpy(buf,buf+BUFFER_SIZE,BUFFER_SIZE);
        offset += BUFFER_SIZE;
        bytes_read = read(fd,buf+BUFFER_SIZE,BUFFER_SIZE);
        if (bytes_read < BUFFER_SIZE) {
          if (bytes_read < 0) {
            cerr << "Error: tailgrepxcol bytes_read: " << bytes_read << endl;
            assert(0);
          }
          if (buf[BUFFER_SIZE+bytes_read-1] != '\n') {
            buf[BUFFER_SIZE+bytes_read] = '\n';
            ++bytes_read;
          }
          pos_eof = BUFFER_SIZE+bytes_read;
        }
        bytes_read += BUFFER_SIZE;
        pos -= BUFFER_SIZE;
      }

      // use memchr to find the next newline...
      const char * nl = (char*)memchr(buf+pos,'\n',bytes_read-pos);
      if (nl == NULL) {
        // if you hit this error, it means the buffer is not long enough to read the file
        // until a newline is found. This should produce an error?...
        cout << "Error: tailgrepxcol file does not have a newline in " << bytes_read-pos << " chars. Increase BUFFER_SIZE in xcol." << endl;
        assert(0);
      }
      pos_nl = int(nl-buf);

      // no point in checking for anything in a line shorter than "expr"...
      if (pos_nl-pos < nexpr) {
        pos = pos_nl + 1;
        continue;
      }

      // general grep on expression...
      char * p = buf+pos;
      while ((p = (char*)memchr(p,expr[0],(buf+pos_nl)-p-nexpr+1))) {
        if (memcmp(p,expr,nexpr) == 0)
          break;
        ++p;
        assert((buf+pos_nl)-p-nexpr+1 >= 0);
      }
      if (!p) {
        pos = pos_nl + 1;
        continue;
      }

      // clear leading whitespace once...
      while ((pos < pos_nl)&&((buf[pos]==' ')||(buf[pos]=='\t')))
        ++pos;

      // advance through the first columns to get to column1...
      int ic = 1;
      while (ic < column1) {
        // advance through non-whitespace...
        while ((pos < pos_nl)&&(buf[pos]!=' ')&&(buf[pos]!='\t')&&(buf[pos]!=','))
          ++pos;
        bool b_comma = false;
        // then advance through whitespace which may contain up to 1 comma...
        while ((pos < pos_nl)&&((buf[pos]==' ')||(buf[pos]=='\t')||(buf[pos]==','))) {
          if (buf[pos]==',') {
            // if we hit a second comma, then break without advancing pos...
            if (b_comma)
              break;
            b_comma = true;
          }
          ++pos;
        }
        ++ic;
      }

      if (pos == pos_nl) {
        // this row is missing the requested column,
        // so set to zero...
        doubleVec.push_back(0.0);
      }
      else {
        // otherwise parse double...
        doubleVec.push_back(strtod(buf+pos,NULL));
      }

      while (ic < column2) {
        // advance through non-whitespace...
        while ((pos < pos_nl)&&(buf[pos]!=' ')&&(buf[pos]!='\t')&&(buf[pos]!=','))
          ++pos;
        bool b_comma = false;
        // then advance through whitespace which may contain up to 1 comma...
        while ((pos < pos_nl)&&((buf[pos]==' ')||(buf[pos]=='\t')||(buf[pos]==','))) {
          if (buf[pos]==',') {
            // if we hit a second comma, then break without advancing pos...
            if (b_comma)
              break;
            b_comma = true;
          }
          ++pos;
        }
        ++ic;
      }

      if (pos == pos_nl) {
        // this row is missing the requested column,
        // so set to zero...
        doubleVec.push_back(0.0);
      }
      else {
        // otherwise parse double...
        doubleVec.push_back(strtod(buf+pos,NULL));
      }

      // and advance pos to next line...
      pos = pos_nl+1;

      if (b_init) {
        if (doubleVec.size() == n1*2) {
          offset1 = offset+pos_nl;
        }
        else if (doubleVec.size() == n2*2) {
          // when we hit n2, do the prediction...
          const size_t tail_size = (n+max(64,n/8))*((offset+pos_nl-offset1)/(n2-n1)+1);
          if (tail_size+offset+bytes_read < stat_buf.st_size) {
            // the tail_size is less than the remaining part of the file, so advance and
            // read just the tail_size...
            size_t target = stat_buf.st_size-tail_size;
            lseek(fd,target,SEEK_SET);
            doubleVec.clear();
            assert(pos_eof == -1);
            bytes_read = read(fd, buf, BUFFER_SIZE*2);
            if (bytes_read < BUFFER_SIZE*2) {
              if (buf[bytes_read-1] != '\n') {
                buf[bytes_read++] = '\n';
              }
              pos_eof = bytes_read;
            }
            const char * nl = (char*)memchr(buf+pos,'\n',bytes_read-pos);
            assert(nl != NULL);
            pos = (nl-buf)+1; // advances to the start of the next line...
          }
          b_init = false;
        }
      }

    }

    close(fd);

    if (doubleVec.size() < n*2) {
      if (!b_init) {
        cout << "possible prediction problem: n requested: " << n << ", n tailed: " << doubleVec.size()/2 << endl;
      }
      n = doubleVec.size()/2;
      if (data1 == NULL) data1 = new double[n];
      if (data2 == NULL) data2 = new double[n];
      for (int i = 0; i < n; ++i) {
        data1[i] = doubleVec[2*i];
        data2[i] = doubleVec[2*i+1];
      }
    }
    else {
      const int dn2 = doubleVec.size() - n*2;
      if (data1 == NULL) data1 = new double[n];
      if (data2 == NULL) data2 = new double[n];
      for (int i = 0; i < n; ++i) {
        data1[i] = doubleVec[dn2+2*i];
        data2[i] = doubleVec[dn2+2*i+1];
      }
    }

    return 0;

  }

  // ======================================================================
  // ROBUST (hopefully) double precision tri-tri intersection routines...
  // ======================================================================

  void getTriTriX(double x[3],double *x0[3],double *x1[3],const int kind,int * idata,double * ddata) {

    switch (kind) {
    case NODE_NODE_INT:
      {
	double xi0[3]; FOR_I3 xi0[i] = x0[idata[0]][i];
	double xi1[3]; FOR_I3 xi1[i] = x1[idata[1]][i];
	if (DIST(xi0,xi1) > 1.0E-12) {
          cout << "Warning: NODE_NODE_INT: DIST(xi0,xi1) > 1.0E-12: " << DIST(xi0,xi1) << endl;
        }
      }
      FOR_I3 x[i] = 0.5*(x0[idata[0]][i] + x1[idata[1]][i]);
      return;
    case EDGE_TRI_INT:
      // idata[0] contains the index of the edge in tri0...
      // store wgt wrt the head of the edge so sorted intersections
      // will go from tail to head...
      FOR_I3 x[i] = ddata[0]*x0[(idata[0]+1)%3][i] + (1.0-ddata[0])*x0[idata[0]][i];
      return;
    case TRI_EDGE_INT:
      // idata[0] contains the index of the edge in tri1...
      FOR_I3 x[i] = ddata[0]*x1[(idata[0]+1)%3][i] + (1.0-ddata[0])*x1[idata[0]][i];
      return;
    case EDGE_EDGE_INT:
      // idata[0] contains the index of the edge in tri0...
      {
	double xi0[3]; FOR_I3 xi0[i] = ddata[0]*x0[(idata[0]+1)%3][i] + (1.0-ddata[0])*x0[idata[0]][i];
	double xi1[3]; FOR_I3 xi1[i] = ddata[1]*x1[(idata[1]+1)%3][i] + (1.0-ddata[1])*x1[idata[1]][i];
	if (DIST(xi0,xi1) > 1.0E-12) {
          cout << "Warning: EDGE_EDGE_INT: DIST(xi0,xi1) > 1.0E-12: " << DIST(xi0,xi1) << endl;
        }
      }
      FOR_I3 x[i] = 0.5*( ddata[0]*x0[(idata[0]+1)%3][i] + (1.0-ddata[0])*x0[idata[0]][i] +
			  ddata[1]*x1[(idata[1]+1)%3][i] + (1.0-ddata[1])*x1[idata[1]][i] );
      return;
    case NODE_EDGE_INT:
      {
        double xi1[3]; FOR_I3 xi1[i] = ddata[0]*x1[(idata[1]+1)%3][i] + (1.0-ddata[0])*x1[idata[1]][i];
        if (DIST(x0[idata[0]],xi1) > 1.0E-12) {
          cout << "Warning: NODE_EDGE_INT: DIST(x0[idata[0]],xi1) > 1.0E-12: " << DIST(x0[idata[0]],xi1) << endl;
        }
      }
      // return the node...
      FOR_I3 x[i] = x0[idata[0]][i];
      return;
    case EDGE_NODE_INT:
      {
        double xi0[3]; FOR_I3 xi0[i] = ddata[0]*x0[(idata[0]+1)%3][i] + (1.0-ddata[0])*x0[idata[0]][i];
        if (DIST(xi0,x1[idata[1]]) > 1.0E-12) {
          cout << "Warning: EDGE_NODE_INT: DIST(xi0,x1[idata[1]]) > 1.0E-12: " << DIST(xi0,x1[idata[1]]) << endl;
          assert(0);
        }
      }
      // return the node...
      FOR_I3 x[i] = x1[idata[1]][i];
      return;
    default:
      assert(0); // not handled yet
    }

  }


  void writeTriTriBin(const int index,double *x0[3],double *x1[3]) {

    char filename[128];
    sprintf(filename,"tritri.%06d.bin",index);
    FILE * fp = fopen(filename,"wb");
    FOR_I3 fwrite(x0[i],sizeof(double),3,fp);
    FOR_I3 fwrite(x1[i],sizeof(double),3,fp);
    fclose(fp);

  }

  void writeTriTriBin(const string& filename,double *x0[3],double *x1[3]) {

    FILE * fp = fopen(filename.c_str(),"wb");
    FOR_I3 fwrite(x0[i],sizeof(double),3,fp);
    FOR_I3 fwrite(x1[i],sizeof(double),3,fp);
    fclose(fp);

  }

  void writeTriTriBinDebug(const int index,double *x0[3],double *x1[3]) {

    char filename[128];
    sprintf(filename,"debug/tritri.%06d.bin",index);
    FILE * fp = fopen(filename,"wb");
    FOR_I3 fwrite(x0[i],sizeof(double),3,fp);
    FOR_I3 fwrite(x1[i],sizeof(double),3,fp);
    fclose(fp);

  }

  void readTriTriBin(double *x0[3],double *x1[3],const string& filename) {

    FILE * fp = fopen(filename.c_str(),"rb");
    FOR_I3 fread(x0[i],sizeof(double),3,fp);
    FOR_I3 fread(x1[i],sizeof(double),3,fp);
    fclose(fp);

  }

  void writeTriTri(const int index,double *x0[3],double *x1[3]) {

    char filename[128];
    sprintf(filename,"tri0.%06d.dat",index);
    cout << " > writing ASCII tri file: " << filename << endl;

    FILE * fp = fopen(filename,"w");
    // write x0,x1,x2,x0 so tri closes on itself...
    fprintf(fp,"%18.15e %18.15e %18.15e\n",x0[0][0],x0[0][1],x0[0][2]);
    fprintf(fp,"%18.15e %18.15e %18.15e\n",x0[1][0],x0[1][1],x0[1][2]);
    fprintf(fp,"%18.15e %18.15e %18.15e\n",x0[2][0],x0[2][1],x0[2][2]);
    fprintf(fp,"%18.15e %18.15e %18.15e\n",x0[0][0],x0[0][1],x0[0][2]);
    fclose(fp);

    sprintf(filename,"tri1.%06d.dat",index);
    cout << " > writing ASCII tri file: " << filename << endl;

    fp = fopen(filename,"w");
    fprintf(fp,"%18.15e %18.15e %18.15e\n",x1[0][0],x1[0][1],x1[0][2]);
    fprintf(fp,"%18.15e %18.15e %18.15e\n",x1[1][0],x1[1][1],x1[1][2]);
    fprintf(fp,"%18.15e %18.15e %18.15e\n",x1[2][0],x1[2][1],x1[2][2]);
    fprintf(fp,"%18.15e %18.15e %18.15e\n",x1[0][0],x1[0][1],x1[0][2]);
    fclose(fp);

  }

  void writeEdge(FILE * fp,const double x0[3],const double x1[3]) {
    const int n = 50;
    for (int i = 0; i <= n; ++i) {
      const double wgt = pow(double(i)/double(n),2);
      fprintf(fp,"%18.16e %18.16e %18.16e\n",
              wgt*x1[0]+(1.0-wgt)*x0[0],
              wgt*x1[1]+(1.0-wgt)*x0[1],
              wgt*x1[2]+(1.0-wgt)*x0[2]);
    }
  }

  void writeEdge(const string& filename,const double x0[3],const double x1[3]) {
    FILE * fp = fopen(filename.c_str(),"w");
    assert(fp);
    writeEdge(fp,x0,x1);
    fclose(fp);
  }

  void writeEdge(const int ied,const double x0[3],const double x1[3]) {
    char filename[128];
    sprintf(filename,"edge.%06d.dat",ied);
    writeEdge(filename,x0,x1);
  }

  set<uint> caseSet; // used to store tritri cases and prevent re-output of CODE for a case already output

  int calcTriTriIntersection(int idata[6],double ddata[4],double *x0[3],double *x1[3]) {

    // start with a diagonal bbox check...

    // +x+y...
    if ((min(x0[0][0]+x0[0][1],min(x0[1][0]+x0[1][1],x0[2][0]+x0[2][1])) >
         max(x1[0][0]+x1[0][1],max(x1[1][0]+x1[1][1],x1[2][0]+x1[2][1]))) ||
        (min(x1[0][0]+x1[0][1],min(x1[1][0]+x1[1][1],x1[2][0]+x1[2][1])) >
         max(x0[0][0]+x0[0][1],max(x0[1][0]+x0[1][1],x0[2][0]+x0[2][1]))))
      return 0;

    // +x-y...
    if ((min(x0[0][0]-x0[0][1],min(x0[1][0]-x0[1][1],x0[2][0]-x0[2][1])) >
         max(x1[0][0]-x1[0][1],max(x1[1][0]-x1[1][1],x1[2][0]-x1[2][1]))) ||
        (min(x1[0][0]-x1[0][1],min(x1[1][0]-x1[1][1],x1[2][0]-x1[2][1])) >
         max(x0[0][0]-x0[0][1],max(x0[1][0]-x0[1][1],x0[2][0]-x0[2][1]))))
      return 0;

    // +x+z...
    if ((min(x0[0][0]+x0[0][2],min(x0[1][0]+x0[1][2],x0[2][0]+x0[2][2])) >
         max(x1[0][0]+x1[0][2],max(x1[1][0]+x1[1][2],x1[2][0]+x1[2][2]))) ||
        (min(x1[0][0]+x1[0][2],min(x1[1][0]+x1[1][2],x1[2][0]+x1[2][2])) >
         max(x0[0][0]+x0[0][2],max(x0[1][0]+x0[1][2],x0[2][0]+x0[2][2]))))
      return 0;

    // +x-z...
    if ((min(x0[0][0]-x0[0][2],min(x0[1][0]-x0[1][2],x0[2][0]-x0[2][2])) >
         max(x1[0][0]-x1[0][2],max(x1[1][0]-x1[1][2],x1[2][0]-x1[2][2]))) ||
        (min(x1[0][0]-x1[0][2],min(x1[1][0]-x1[1][2],x1[2][0]-x1[2][2])) >
         max(x0[0][0]-x0[0][2],max(x0[1][0]-x0[1][2],x0[2][0]-x0[2][2]))))
      return 0;

    // +y+z...
    if ((min(x0[0][1]+x0[0][2],min(x0[1][1]+x0[1][2],x0[2][1]+x0[2][2])) >
         max(x1[0][1]+x1[0][2],max(x1[1][1]+x1[1][2],x1[2][1]+x1[2][2]))) ||
        (min(x1[0][1]+x1[0][2],min(x1[1][1]+x1[1][2],x1[2][1]+x1[2][2])) >
         max(x0[0][1]+x0[0][2],max(x0[1][1]+x0[1][2],x0[2][1]+x0[2][2]))))
      return 0;

    // +y-z...
    if ((min(x0[0][1]-x0[0][2],min(x0[1][1]-x0[1][2],x0[2][1]-x0[2][2])) >
         max(x1[0][1]-x1[0][2],max(x1[1][1]-x1[1][2],x1[2][1]-x1[2][2]))) ||
        (min(x1[0][1]-x1[0][2],min(x1[1][1]-x1[1][2],x1[2][1]-x1[2][2])) >
         max(x0[0][1]-x0[0][2],max(x0[1][1]-x0[1][2],x0[2][1]-x0[2][2]))))
      return 0;

    // +x+y+z...
    if ((min(x0[0][0]+x0[0][1]+x0[0][2],min(x0[1][0]+x0[1][1]+x0[1][2],x0[2][0]+x0[2][1]+x0[2][2])) >
         max(x1[0][0]+x1[0][1]+x1[0][2],max(x1[1][0]+x1[1][1]+x1[1][2],x1[2][0]+x1[2][1]+x1[2][2]))) ||
        (min(x1[0][0]+x1[0][1]+x1[0][2],min(x1[1][0]+x1[1][1]+x1[1][2],x1[2][0]+x1[2][1]+x1[2][2])) >
         max(x0[0][0]+x0[0][1]+x0[0][2],max(x0[1][0]+x0[1][1]+x0[1][2],x0[2][0]+x0[2][1]+x0[2][2]))))
      return 0;

    // +x+y-z...
    if ((min(x0[0][0]+x0[0][1]-x0[0][2],min(x0[1][0]+x0[1][1]-x0[1][2],x0[2][0]+x0[2][1]-x0[2][2])) >
         max(x1[0][0]+x1[0][1]-x1[0][2],max(x1[1][0]+x1[1][1]-x1[1][2],x1[2][0]+x1[2][1]-x1[2][2]))) ||
        (min(x1[0][0]+x1[0][1]-x1[0][2],min(x1[1][0]+x1[1][1]-x1[1][2],x1[2][0]+x1[2][1]-x1[2][2])) >
         max(x0[0][0]+x0[0][1]-x0[0][2],max(x0[1][0]+x0[1][1]-x0[1][2],x0[2][0]+x0[2][1]-x0[2][2]))))
      return 0;

    // +x-y+z...
    if ((min(x0[0][0]-x0[0][1]+x0[0][2],min(x0[1][0]-x0[1][1]+x0[1][2],x0[2][0]-x0[2][1]+x0[2][2])) >
         max(x1[0][0]-x1[0][1]+x1[0][2],max(x1[1][0]-x1[1][1]+x1[1][2],x1[2][0]-x1[2][1]+x1[2][2]))) ||
        (min(x1[0][0]-x1[0][1]+x1[0][2],min(x1[1][0]-x1[1][1]+x1[1][2],x1[2][0]-x1[2][1]+x1[2][2])) >
         max(x0[0][0]-x0[0][1]+x0[0][2],max(x0[1][0]-x0[1][1]+x0[1][2],x0[2][0]-x0[2][1]+x0[2][2]))))
      return 0;

    // +x-y-z...
    if ((min(x0[0][0]-x0[0][1]-x0[0][2],min(x0[1][0]-x0[1][1]-x0[1][2],x0[2][0]-x0[2][1]-x0[2][2])) >
         max(x1[0][0]-x1[0][1]-x1[0][2],max(x1[1][0]-x1[1][1]-x1[1][2],x1[2][0]-x1[2][1]-x1[2][2]))) ||
        (min(x1[0][0]-x1[0][1]-x1[0][2],min(x1[1][0]-x1[1][1]-x1[1][2],x1[2][0]-x1[2][1]-x1[2][2])) >
         max(x0[0][0]-x0[0][1]-x0[0][2],max(x0[1][0]-x0[1][1]-x0[1][2],x0[2][0]-x0[2][1]-x0[2][2]))))
      return 0;

    uint value = 0;
    uint offset = 1;

    int p0[3] = { 0, 0, 0 };

    // TODO: change this from a vector to a scalar -- this will require a complete redo of the
    // codegen, but could be worth it!...

    double vp0[3][3];
    FOR_I3 {
      // the volume of the tet from a particular x0 to tri x1...
      vp0[i][0] = SIGNED_TET_VOLUME_6(x0[i],x1[0],x1[1],x1[2]);
      vp0[i][1] = SIGNED_TET_VOLUME_6(x0[i],x1[1],x1[2],x1[0]);
      vp0[i][2] = SIGNED_TET_VOLUME_6(x0[i],x1[2],x1[0],x1[1]);
      if ((vp0[i][0] > 0.0)&&(vp0[i][1] > 0.0)&&(vp0[i][2] > 0.0)) p0[i] = 1;
      else if ((vp0[i][0] < 0.0)&&(vp0[i][1] < 0.0)&&(vp0[i][2] < 0.0)) p0[i] = -1;
      value += offset*(p0[i]+1);
      offset *= 3;
    }

    // tri0 is entirely on one side of tri1...
    if ( ((p0[0] == 1)&&(p0[1] == 1)&&(p0[2] == 1)) || ((p0[0] == -1)&&(p0[1] == -1)&&(p0[2] == -1)) )
      return 0;

    if ((p0[0] == 0)&&(p0[1] == 0)&&(p0[2] == 0)) {
      //cout << "coplanar tris! skipping" << endl;
      return -1;
    }

    int p1[3] = { 0, 0, 0 };
    double vp1[3][3];
    FOR_I3 {
      // the volume of the tet from a particular x0 to tri x1...
      vp1[i][0] = SIGNED_TET_VOLUME_6(x1[i],x0[0],x0[1],x0[2]);
      vp1[i][1] = SIGNED_TET_VOLUME_6(x1[i],x0[1],x0[2],x0[0]);
      vp1[i][2] = SIGNED_TET_VOLUME_6(x1[i],x0[2],x0[0],x0[1]);
      if ((vp1[i][0] > 0.0)&&(vp1[i][1] > 0.0)&&(vp1[i][2] > 0.0)) p1[i] = 1;
      else if ((vp1[i][0] < 0.0)&&(vp1[i][1] < 0.0)&&(vp1[i][2] < 0.0)) p1[i] = -1;
      value += offset*(p1[i]+1);
      offset *= 3;
    }

    // tri1 is entirely on one side of tri0...
    if ( ((p1[0] == 1)&&(p1[1] == 1)&&(p1[2] == 1)) || ((p1[0] == -1)&&(p1[1] == -1)&&(p1[2] == -1)) )
      return 0;

    if ((p1[0] == 0)&&(p1[1] == 0)&&(p1[2] == 0)) {
      //cout << "coplanar tris! skipping" << endl;
      return -1;
    }

    // consider edges of tri0 against tri1...

    int e0[3][3] = { { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 } };
    double e0d[3][3] = { { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 } };
    FOR_I3 {
      if (p0[i] != p0[(i+1)%3]) {
        FOR_J3 {
          const double v0 = SIGNED_TET_VOLUME_6(x0[i],x0[(i+1)%3],x1[j],x1[(j+1)%3]);
          const double v1 = SIGNED_TET_VOLUME_6(x0[(i+1)%3],x0[i],x1[(j+1)%3],x1[j]);
          const double v2 = SIGNED_TET_VOLUME_6(x1[j],x1[(j+1)%3],x0[i],x0[(i+1)%3]);
          const double v3 = SIGNED_TET_VOLUME_6(x1[(j+1)%3],x1[j],x0[(i+1)%3],x0[i]);
          if ((v0 > 0.0)&&(v1 > 0.0)&&(v2 > 0.0)&&(v3 > 0.0)) {
            e0[i][j] = 1;
            e0d[i][j] = v0 + v1 + v2 + v3;
          }
          else if ((v0 < 0.0)&&(v1 < 0.0)&&(v2 < 0.0)&&(v3 < 0.0)) {
            e0[i][j] = -1;
            e0d[i][j] = v0 + v1 + v2 + v3;
          }
          value += offset*(e0[i][j]+1);
          offset *= 3;
        }
      }
    }

    // consider edges of tri1 against tri0...

    int e1[3][3] = { { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 } };
    double e1d[3][3] = { { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 } };
    FOR_I3 {
      if (p1[i] != p1[(i+1)%3]) {
        FOR_J3 {
          const double v0 = SIGNED_TET_VOLUME_6(x1[i],x1[(i+1)%3],x0[j],x0[(j+1)%3]);
          const double v1 = SIGNED_TET_VOLUME_6(x1[(i+1)%3],x1[i],x0[(j+1)%3],x0[j]);
          const double v2 = SIGNED_TET_VOLUME_6(x0[j],x0[(j+1)%3],x1[i],x1[(i+1)%3]);
          const double v3 = SIGNED_TET_VOLUME_6(x0[(j+1)%3],x0[j],x1[(i+1)%3],x1[i]);
          if ((v0 > 0.0)&&(v1 > 0.0)&&(v2 > 0.0)&&(v3 > 0.0)) {
            e1[i][j] = 1;
            e1d[i][j] = v0 + v1 + v2 + v3;
          }
          else if ((v0 < 0.0)&&(v1 < 0.0)&&(v2 < 0.0)&&(v3 < 0.0)) {
            e1[i][j] = -1;
            e1d[i][j] = v0 + v1 + v2 + v3;
          }
          value += offset*(e1[i][j]+1);
          offset *= 3;
        }
      }
    }

    switch (value) {

#include "tritri_rand.1.hpp"
#include "tritri_rand.2.hpp"
#include "tcc.1.hpp"
#include "tcc.2.hpp"
#include "tcc.3.hpp"
#include "tcc.4.hpp"
#include "node_node_int.hpp"
#include "node_edge_int.hpp"

    }

    cout << "========================= Working on =======================" << endl;
    FOR_I3 cout << "p0: " << p0[i] << " " << p0[(i+1)%3] << " e0" << i << " " << COUT_VEC(e0[i]) << endl;
    FOR_I3 cout << "p1: " << p1[i] << " " << p1[(i+1)%3] << " e1" << i << " " << COUT_VEC(e1[i]) << endl;

    vector<TriTriInt> triTriIntVec;

    FOR_I3 {
      if (p0[i] != p0[(i+1)%3]) {

        // this edge of tri 0 may have an intersection...
        if ( ((e0[i][0] >= 0)&&(e0[i][1] >= 0)&&(e0[i][2] >= 0)) || ((e0[i][0] <= 0)&&(e0[i][1] <= 0)&&(e0[i][2] <= 0)) ) {

          // check orientation of the edge through the tri: should relate to the sign of the edge volumes...
          if (p0[(i+1)%3]-p0[i] > 0) {
            if (!(((e0[i][0] <= 0)&&(e0[i][1] <= 0)&&(e0[i][2] <= 0)))) {
              cout << "Warning: e0 sign problem plus, tri 0 edge " << i << endl;
              continue;
            }
          }
          else {
            assert(p0[(i+1)%3]-p0[i] < 0);
            if (!(((e0[i][0] >= 0)&&(e0[i][1] >= 0)&&(e0[i][2] >= 0)))) {
              cout << "Warning: e0 sign problem minus, tri 0 edge " << i << endl;
              continue;
            }
          }

          // add an intersection...
          if ( ((e0[i][0] > 0)&&(e0[i][1] > 0)&&(e0[i][2] > 0)) || ((e0[i][0] < 0)&&(e0[i][1] < 0)&&(e0[i][2] < 0)) ) {

            // this edge0 is away from the tri1 edges...

            if (p0[i]*p0[(i+1)%3] == -1) {
              // this edge0 has its nodes on either side of tri1...
              stringstream ss;
              ss << "(vp0["<<i<<"][0]+vp0["<<i<<"][1]+vp0["<<i<<"][2])/(vp0["<<i<<"][0]+vp0["<<i<<"][1]+vp0["<<i<<"][2]-vp0["<<((i+1)%3)<<"][0]-vp0["<<((i+1)%3)<<"][1]-vp0["<<((i+1)%3)<<"][2])";
              triTriIntVec.push_back(TriTriInt(EDGE_TRI_INT,i,(vp0[i][0]+vp0[i][1]+vp0[i][2])/(vp0[i][0]+vp0[i][1]+vp0[i][2]-vp0[(i+1)%3][0]-vp0[(i+1)%3][1]-vp0[(i+1)%3][2]),ss.str()));
            }
            else {
              writeTriTriBin(mpi_rank,x0,x1);
              writeTriTri(mpi_rank,x0,x1);
              assert(0);
            }
          }
          else {

            // we are on an edge or a node...
            if ( (e0[i][0] == 0)&&((p0[(i+1)%3]-p0[i])*e0[i][1] < 0)&&((p0[(i+1)%3]-p0[i])*e0[i][2] < 0) ) {
              // edge-edge with tri1's edge0...

              // here we use the edge tet volumes to compute the weight...
              stringstream ss1;
              ss1 << "e0d["<<i<<"][2]/(e0d["<<i<<"][1]+e0d["<<i<<"][2])";
              const double wgt1 = e0d[i][2]/(e0d[i][1]+e0d[i][2]);

              if (p0[i]*p0[(i+1)%3] == -1) {
                // this edge0 has its nodes on either side of tri1...
                stringstream ss0;
                ss0 << "(vp0["<<i<<"][0]+vp0["<<i<<"][1]+vp0["<<i<<"][2])/(vp0["<<i<<"][0]+vp0["<<i<<"][1]+vp0["<<i<<"][2]-vp0["<<((i+1)%3)<<"][0]-vp0["<<((i+1)%3)<<"][1]-vp0["<<((i+1)%3)<<"][2])";
                triTriIntVec.push_back(TriTriInt(EDGE_EDGE_INT,i,0,
                                                 (vp0[i][0]+vp0[i][1]+vp0[i][2])/(vp0[i][0]+vp0[i][1]+vp0[i][2]-vp0[(i+1)%3][0]-vp0[(i+1)%3][1]-vp0[(i+1)%3][2]),ss0.str(),
                                                 wgt1,ss1.str()));
              }
              else if (p0[i] == 0) {
                assert(p0[(i+1)%3] != 0);
                // we are on node i of edge0...
                triTriIntVec.push_back(TriTriInt(NODE_EDGE_INT,i,0,
                                                 wgt1,ss1.str()));
              }
              else {
                assert(p0[(i+1)%3] == 0);
                //assert(p0[i] != 0); // guaranteed -- see previous "else if"
                triTriIntVec.push_back(TriTriInt(NODE_EDGE_INT,(i+1)%3,0,
                                                 wgt1,ss1.str()));
              }

            }
            else if ( (e0[i][1] == 0)&&((p0[(i+1)%3]-p0[i])*e0[i][2] < 0)&&((p0[(i+1)%3]-p0[i])*e0[i][0] < 0) ) {
              // edge-edge with tri1's edge1...

              // here we use the edge tet volumes to compute the weight...
              stringstream ss1;
              ss1 << "e0d["<<i<<"][0]/(e0d["<<i<<"][2]+e0d["<<i<<"][0])";
              const double wgt1 = e0d[i][0]/(e0d[i][2]+e0d[i][0]);

              if (p0[i]*p0[(i+1)%3] == -1) {
                // this edge0 has its nodes on either side of tri1...
                stringstream ss0;
                ss0 << "(vp0["<<i<<"][0]+vp0["<<i<<"][1]+vp0["<<i<<"][2])/(vp0["<<i<<"][0]+vp0["<<i<<"][1]+vp0["<<i<<"][2]-vp0["<<((i+1)%3)<<"][0]-vp0["<<((i+1)%3)<<"][1]-vp0["<<((i+1)%3)<<"][2])";
                triTriIntVec.push_back(TriTriInt(EDGE_EDGE_INT,i,1,
                                                 (vp0[i][0]+vp0[i][1]+vp0[i][2])/(vp0[i][0]+vp0[i][1]+vp0[i][2]-vp0[(i+1)%3][0]-vp0[(i+1)%3][1]-vp0[(i+1)%3][2]),ss0.str(),
                                                 wgt1,ss1.str()));
              }
              else if (p0[i] == 0) {
                assert(p0[(i+1)%3] != 0);
                // we are on node i of edge0...
                triTriIntVec.push_back(TriTriInt(NODE_EDGE_INT,i,1,
                                                 wgt1,ss1.str()));
              }
              else {
                assert(p0[(i+1)%3] == 0);
                //assert(p0[i] != 0); // guaranteed -- see previous "else if"
                triTriIntVec.push_back(TriTriInt(NODE_EDGE_INT,(i+1)%3,1,
                                                 wgt1,ss1.str()));
              }

            }
            else if ( (e0[i][2] == 0)&&((p0[(i+1)%3]-p0[i])*e0[i][0] < 0)&&((p0[(i+1)%3]-p0[i])*e0[i][1] < 0) ) {
              // edge-edge with tri1's edge2...

              // here we use the edge tet volumes to compute the weight...
              stringstream ss1;
              ss1 << "e0d["<<i<<"][1]/(e0d["<<i<<"][0]+e0d["<<i<<"][1])";
              const double wgt1 = e0d[i][1]/(e0d[i][0]+e0d[i][1]);

              if (p0[i]*p0[(i+1)%3] == -1) {
                // this edge0 has its nodes on either side of tri1...
                stringstream ss0;
                ss0 << "(vp0["<<i<<"][0]+vp0["<<i<<"][1]+vp0["<<i<<"][2])/(vp0["<<i<<"][0]+vp0["<<i<<"][1]+vp0["<<i<<"][2]-vp0["<<((i+1)%3)<<"][0]-vp0["<<((i+1)%3)<<"][1]-vp0["<<((i+1)%3)<<"][2])";
                triTriIntVec.push_back(TriTriInt(EDGE_EDGE_INT,i,2,
                                                 (vp0[i][0]+vp0[i][1]+vp0[i][2])/(vp0[i][0]+vp0[i][1]+vp0[i][2]-vp0[(i+1)%3][0]-vp0[(i+1)%3][1]-vp0[(i+1)%3][2]),ss0.str(),
                                                 wgt1,ss1.str()));
              }
              else if (p0[i] == 0) {
                assert(p0[(i+1)%3] != 0);
                // we are on node i of edge0...
                triTriIntVec.push_back(TriTriInt(NODE_EDGE_INT,i,2,
                                                 wgt1,ss1.str()));
              }
              else {
                assert(p0[(i+1)%3] == 0);
                //assert(p0[i] != 0); // guaranteed -- see previous "else if"
                triTriIntVec.push_back(TriTriInt(NODE_EDGE_INT,(i+1)%3,2,
                                                 wgt1,ss1.str()));
              }

            }
            else if ( (e0[i][0] == 0)&&(e0[i][1] == 0)&&((p0[(i+1)%3]-p0[i])*e0[i][2] < 0) ) {

              // edge "i" of tri0 is intersecting node on tri1:
              if (p0[i] == 0) {
                cout << "tri0 node " << i << " and tri1 node 1" << endl;
                triTriIntVec.push_back(TriTriInt(NODE_NODE_INT,i,1));
              }
              else if (p0[(i+1)%3] == 0) {
                cout << "tri0 node " << ((i+1)%3) << " and tri1 node 1" << endl;
                triTriIntVec.push_back(TriTriInt(NODE_NODE_INT,(i+1)%3,1));
              }
              else {
                stringstream ss0;
                ss0 << "(vp0["<<i<<"][0]+vp0["<<i<<"][1]+vp0["<<i<<"][2])/(vp0["<<i<<"][0]+vp0["<<i<<"][1]+vp0["<<i<<"][2]-vp0["<<((i+1)%3)<<"][0]-vp0["<<((i+1)%3)<<"][1]-vp0["<<((i+1)%3)<<"][2])";
                triTriIntVec.push_back(TriTriInt(EDGE_NODE_INT,i,1,
                                                 (vp0[i][0]+vp0[i][1]+vp0[i][2])/(vp0[i][0]+vp0[i][1]+vp0[i][2]-vp0[(i+1)%3][0]-vp0[(i+1)%3][1]-vp0[(i+1)%3][2]),ss0.str()));
              }

            }
            else if ( (e0[i][1] == 0)&&(e0[i][2] == 0)&&((p0[(i+1)%3]-p0[i])*e0[i][0] < 0) ) {

              // edge "i" of tri0 is intersecting node on tri1:
              if (p0[i] == 0) {
                cout << "tri0 node " << i << " and tri1 node 2" << endl;
                triTriIntVec.push_back(TriTriInt(NODE_NODE_INT,i,2));
              }
              else if (p0[(i+1)%3] == 0) {
                cout << "tri0 node " << ((i+1)%3) << " and tri1 node 2" << endl;
                triTriIntVec.push_back(TriTriInt(NODE_NODE_INT,(i+1)%3,2));
              }
              else {
                stringstream ss0;
                ss0 << "(vp0["<<i<<"][0]+vp0["<<i<<"][1]+vp0["<<i<<"][2])/(vp0["<<i<<"][0]+vp0["<<i<<"][1]+vp0["<<i<<"][2]-vp0["<<((i+1)%3)<<"][0]-vp0["<<((i+1)%3)<<"][1]-vp0["<<((i+1)%3)<<"][2])";
                triTriIntVec.push_back(TriTriInt(EDGE_NODE_INT,i,2,
                                                 (vp0[i][0]+vp0[i][1]+vp0[i][2])/(vp0[i][0]+vp0[i][1]+vp0[i][2]-vp0[(i+1)%3][0]-vp0[(i+1)%3][1]-vp0[(i+1)%3][2]),ss0.str()));
              }

            }
            else if ( (e0[i][2] == 0)&&(e0[i][0] == 0)&&((p0[(i+1)%3]-p0[i])*e0[i][1] < 0) ) {

              // edge "i" of tri0 is intersecting node on tri1:
              if (p0[i] == 0) {
                cout << "tri0 node " << i << " and tri1 node 0" << endl;
                triTriIntVec.push_back(TriTriInt(NODE_NODE_INT,i,0));
              }
              else if (p0[(i+1)%3] == 0) {
                cout << "tri0 node " << ((i+1)%3) << " and tri1 node 0" << endl;
                triTriIntVec.push_back(TriTriInt(NODE_NODE_INT,(i+1)%3,0));
              }
              else {
                stringstream ss0;
                ss0 << "(vp0["<<i<<"][0]+vp0["<<i<<"][1]+vp0["<<i<<"][2])/(vp0["<<i<<"][0]+vp0["<<i<<"][1]+vp0["<<i<<"][2]-vp0["<<((i+1)%3)<<"][0]-vp0["<<((i+1)%3)<<"][1]-vp0["<<((i+1)%3)<<"][2])";
                triTriIntVec.push_back(TriTriInt(EDGE_NODE_INT,i,0,
                                                 (vp0[i][0]+vp0[i][1]+vp0[i][2])/(vp0[i][0]+vp0[i][1]+vp0[i][2]-vp0[(i+1)%3][0]-vp0[(i+1)%3][1]-vp0[(i+1)%3][2]),ss0.str()));
              }

            }
            else {

              if ((e0[i][0] == 0)&&(e0[i][1] == 0)&&(e0[i][2] == 0)) {
                cout << "Warning: e0 all zero dispite edge endpoint diff. Skipping" << endl;
                continue;
              }

              cout << "failing at i: " << i << endl;
              cout << "e0[i][0]: " << e0[i][0] << " e0[i][1]: " << e0[i][1] << " e0[i][2]: " << e0[i][2] << endl;
              writeTriTriBin(mpi_rank,x0,x1);
              writeTriTri(mpi_rank,x0,x1);
              assert(0);
            }

          }

        }
      }
    }

    FOR_I3 {
      if (p1[i] != p1[(i+1)%3]) {

        if ( ((e1[i][0] >= 0)&&(e1[i][1] >= 0)&&(e1[i][2] >= 0)) || ((e1[i][0] <= 0)&&(e1[i][1] <= 0)&&(e1[i][2] <= 0)) ) {

          // check orientation of the edge through the tri: should relate to the sign of the edge volumes...
          if (p1[(i+1)%3]-p1[i] > 0) {
            if (!(((e1[i][0] <= 0)&&(e1[i][1] <= 0)&&(e1[i][2] <= 0)))) {
              cout << "Warning: e1 sign problem plus, tri 1 edge " << i << endl;
              continue;
            }
          }
          else {
            assert(p1[(i+1)%3]-p1[i] < 0);
            if (!(((e1[i][0] >= 0)&&(e1[i][1] >= 0)&&(e1[i][2] >= 0)))) {
              cout << "Warning: e1 sign problem minus, tri 1 edge " << i << endl;
              continue;
            }
          }

          // add an intersection...

          if ( ((e1[i][0] > 0)&&(e1[i][1] > 0)&&(e1[i][2] > 0)) || ((e1[i][0] < 0)&&(e1[i][1] < 0)&&(e1[i][2] < 0)) ) {

            // this edge1 is away from the tri0 edges...

            if (p1[i]*p1[(i+1)%3] == -1) {
              // this edge1 has its nodes on either side of tri0...
              stringstream ss;
              ss << "(vp1["<<i<<"][0]+vp1["<<i<<"][1]+vp1["<<i<<"][2])/(vp1["<<i<<"][0]+vp1["<<i<<"][1]+vp1["<<i<<"][2]-vp1["<<((i+1)%3)<<"][0]-vp1["<<((i+1)%3)<<"][1]-vp1["<<((i+1)%3)<<"][2])";
              triTriIntVec.push_back(TriTriInt(TRI_EDGE_INT,i,(vp1[i][0]+vp1[i][1]+vp1[i][2])/(vp1[i][0]+vp1[i][1]+vp1[i][2]-vp1[(i+1)%3][0]-vp1[(i+1)%3][1]-vp1[(i+1)%3][2]),ss.str()));
            }
            else {
              assert(0);
            }

          }
          else {

            // we are on an edge or a node...
            if ( (e1[i][0] == 0)&&((p1[(i+1)%3]-p1[i])*e1[i][1] < 0)&&((p1[(i+1)%3]-p1[i])*e1[i][2] < 0) ) {
              // edge-edge with tri0's edge0...

              // here we use the edge tet volumes to compute the weight...
              stringstream ss0;
              ss0 << "e1d["<<i<<"][2]/(e1d["<<i<<"][1]+e1d["<<i<<"][2])";
              const double wgt0 = e1d[i][2]/(e1d[i][1]+e1d[i][2]);

              if (p1[i]*p1[(i+1)%3] == -1) {
                // this edge1 has its nodes on either side of tri0...
                stringstream ss1;
                ss1 << "(vp1["<<i<<"][0]+vp1["<<i<<"][1]+vp1["<<i<<"][2])/(vp1["<<i<<"][0]+vp1["<<i<<"][1]+vp1["<<i<<"][2]-vp1["<<((i+1)%3)<<"][0]-vp1["<<((i+1)%3)<<"][1]-vp1["<<((i+1)%3)<<"][2])";
                triTriIntVec.push_back(TriTriInt(EDGE_EDGE_INT,0,i,
                                                 wgt0,ss0.str(),
                                                 (vp1[i][0]+vp1[i][1]+vp1[i][2])/(vp1[i][0]+vp1[i][1]+vp1[i][2]-vp1[(i+1)%3][0]-vp1[(i+1)%3][1]-vp1[(i+1)%3][2]),ss1.str()));
              }
              else if (p1[i] == 0) {
                assert(p1[(i+1)%3] != 0);
                // we are on node i of edge1...
                triTriIntVec.push_back(TriTriInt(EDGE_NODE_INT,0,i,
                                                 wgt0,ss0.str()));
              }
              else {
                assert(p1[(i+1)%3] == 0);
                //assert(p1[i] != 0); // guaranteed -- see previous "else if"
                triTriIntVec.push_back(TriTriInt(EDGE_NODE_INT,0,(i+1)%3,
                                                 wgt0,ss0.str()));
              }

            }
            else if ( (e1[i][1] == 0)&&((p1[(i+1)%3]-p1[i])*e1[i][2] < 0)&&((p1[(i+1)%3]-p1[i])*e1[i][0] < 0) ) {
              // edge-edge with tri0's edge1...

              // here we use the edge tet volumes to compute the weight...
              stringstream ss0;
              ss0 << "e1d["<<i<<"][0]/(e1d["<<i<<"][2]+e1d["<<i<<"][0])";
              const double wgt0 = e1d[i][0]/(e1d[i][2]+e1d[i][0]);

              if (p1[i]*p1[(i+1)%3] == -1) {
                // this edge1 has its nodes on either side of tri0...
                stringstream ss1;
                ss1 << "(vp1["<<i<<"][0]+vp1["<<i<<"][1]+vp1["<<i<<"][2])/(vp1["<<i<<"][0]+vp1["<<i<<"][1]+vp1["<<i<<"][2]-vp1["<<((i+1)%3)<<"][0]-vp1["<<((i+1)%3)<<"][1]-vp1["<<((i+1)%3)<<"][2])";
                triTriIntVec.push_back(TriTriInt(EDGE_EDGE_INT,1,i,
                                                 wgt0,ss0.str(),
                                                 (vp1[i][0]+vp1[i][1]+vp1[i][2])/(vp1[i][0]+vp1[i][1]+vp1[i][2]-vp1[(i+1)%3][0]-vp1[(i+1)%3][1]-vp1[(i+1)%3][2]),ss1.str()));
              }
              else if (p1[i] == 0) {
                assert(p1[(i+1)%3] != 0);
                // we are on node i of edge1...
                triTriIntVec.push_back(TriTriInt(EDGE_NODE_INT,1,i,
                                                 wgt0,ss0.str()));
              }
              else {
                assert(p1[(i+1)%3] == 0);
                //assert(p1[i] != 0); // guaranteed -- see previous "else if"
                triTriIntVec.push_back(TriTriInt(EDGE_NODE_INT,1,(i+1)%3,
                                                 wgt0,ss0.str()));
              }

            }
            else if ( (e1[i][2] == 0)&&((p1[(i+1)%3]-p1[i])*e1[i][0] < 0)&&((p1[(i+1)%3]-p1[i])*e1[i][1] < 0) ) {
              // edge-edge with tri0's edge2...

              // here we use the edge tet volumes to compute the weight...
              stringstream ss0;
              ss0 << "e1d["<<i<<"][1]/(e1d["<<i<<"][0]+e1d["<<i<<"][1])";
              const double wgt0 = e1d[i][1]/(e1d[i][0]+e1d[i][1]);

              if (p1[i]*p1[(i+1)%3] == -1) {
                // this edge1 has its nodes on either side of tri0...
                stringstream ss1;
                ss1 << "(vp1["<<i<<"][0]+vp1["<<i<<"][1]+vp1["<<i<<"][2])/(vp1["<<i<<"][0]+vp1["<<i<<"][1]+vp1["<<i<<"][2]-vp1["<<((i+1)%3)<<"][0]-vp1["<<((i+1)%3)<<"][1]-vp1["<<((i+1)%3)<<"][2])";
                triTriIntVec.push_back(TriTriInt(EDGE_EDGE_INT,2,i,
                                                 wgt0,ss0.str(),
                                                 (vp1[i][0]+vp1[i][1]+vp1[i][2])/(vp1[i][0]+vp1[i][1]+vp1[i][2]-vp1[(i+1)%3][0]-vp1[(i+1)%3][1]-vp1[(i+1)%3][2]),ss1.str()));
              }
              else if (p1[i] == 0) {
                assert(p1[(i+1)%3] != 0);
                // we are on node i of edge1...
                triTriIntVec.push_back(TriTriInt(EDGE_NODE_INT,2,i,
                                                 wgt0,ss0.str()));
              }
              else {
                assert(p1[(i+1)%3] == 0);
                //assert(p1[i] != 0); // guaranteed -- see previous "else if"
                triTriIntVec.push_back(TriTriInt(EDGE_NODE_INT,2,(i+1)%3,
                                                 wgt0,ss0.str()));
              }

            }
            else if ( (e1[i][0] == 0)&&(e1[i][1] == 0)&&((p1[(i+1)%3]-p1[i])*e1[i][2] < 0) ) {

              // edge "i" of tri1 is intersecting node on tri0:
              if (p1[i] == 0) {
                cout << "tri1 node " << i << " and tri0 node 1" << endl;
                triTriIntVec.push_back(TriTriInt(NODE_NODE_INT,1,i));
              }
              else if (p1[(i+1)%3] == 0) {
                cout << "tri1 node " << ((i+1)%3) << " and tri0 node 1" << endl;
                triTriIntVec.push_back(TriTriInt(NODE_NODE_INT,1,(i+1)%3));
              }
              else {
                stringstream ss1;
                ss1 << "(vp1["<<i<<"][0]+vp1["<<i<<"][1]+vp1["<<i<<"][2])/(vp1["<<i<<"][0]+vp1["<<i<<"][1]+vp1["<<i<<"][2]-vp1["<<((i+1)%3)<<"][0]-vp1["<<((i+1)%3)<<"][1]-vp1["<<((i+1)%3)<<"][2])";
                triTriIntVec.push_back(TriTriInt(NODE_EDGE_INT,1,i,
                                                 (vp1[i][0]+vp1[i][1]+vp1[i][2])/(vp1[i][0]+vp1[i][1]+vp1[i][2]-vp1[(i+1)%3][0]-vp1[(i+1)%3][1]-vp1[(i+1)%3][2]),ss1.str()));
              }

            }
            else if ( (e1[i][1] == 0)&&(e1[i][2] == 0)&&((p1[(i+1)%3]-p1[i])*e1[i][0] < 0) ) {

              // edge "i" of tri1 is intersecting node on tri0:
              if (p1[i] == 0) {
                cout << "tri1 node " << i << " and tri0 node 2" << endl;
                triTriIntVec.push_back(TriTriInt(NODE_NODE_INT,2,i));
              }
              else if (p1[(i+1)%3] == 0) {
                cout << "tri1 node " << ((i+1)%3) << " and tri0 node 2" << endl;
                triTriIntVec.push_back(TriTriInt(NODE_NODE_INT,2,(i+1)%3));
              }
              else {
                stringstream ss1;
                ss1 << "(vp1["<<i<<"][0]+vp1["<<i<<"][1]+vp1["<<i<<"][2])/(vp1["<<i<<"][0]+vp1["<<i<<"][1]+vp1["<<i<<"][2]-vp1["<<((i+1)%3)<<"][0]-vp1["<<((i+1)%3)<<"][1]-vp1["<<((i+1)%3)<<"][2])";
                triTriIntVec.push_back(TriTriInt(NODE_EDGE_INT,2,i,
                                                 (vp1[i][0]+vp1[i][1]+vp1[i][2])/(vp1[i][0]+vp1[i][1]+vp1[i][2]-vp1[(i+1)%3][0]-vp1[(i+1)%3][1]-vp1[(i+1)%3][2]),ss1.str()));
              }

            }
            else if ( (e1[i][2] == 0)&&(e1[i][0] == 0)&&((p1[(i+1)%3]-p1[i])*e1[i][1] < 0) ) {

              // edge "i" of tri1 is intersecting node on tri0:
              if (p1[i] == 0) {
                cout << "tri1 node " << i << " and tri0 node 0" << endl;
                triTriIntVec.push_back(TriTriInt(NODE_NODE_INT,0,i));
              }
              else if (p1[(i+1)%3] == 0) {
                cout << "tri1 node " << ((i+1)%3) << " and tri0 node 0" << endl;
                triTriIntVec.push_back(TriTriInt(NODE_NODE_INT,0,(i+1)%3));
              }
              else {
                stringstream ss1;
                ss1 << "(vp1["<<i<<"][0]+vp1["<<i<<"][1]+vp1["<<i<<"][2])/(vp1["<<i<<"][0]+vp1["<<i<<"][1]+vp1["<<i<<"][2]-vp1["<<((i+1)%3)<<"][0]-vp1["<<((i+1)%3)<<"][1]-vp1["<<((i+1)%3)<<"][2])";
                triTriIntVec.push_back(TriTriInt(NODE_EDGE_INT,0,i,
                                                 (vp1[i][0]+vp1[i][1]+vp1[i][2])/(vp1[i][0]+vp1[i][1]+vp1[i][2]-vp1[(i+1)%3][0]-vp1[(i+1)%3][1]-vp1[(i+1)%3][2]),ss1.str()));
              }

            }
            else {

              if ((e1[i][0] == 0)&&(e1[i][1] == 0)&&(e1[i][2] == 0)) {
                cout << "Warning: e1 all zero dispite edge endpoint diff. Skipping" << endl;
                continue;
              }

              cout << "failing at i: " << i << endl;
              cout << "e1[i][0]: " << e1[i][0] << " e1[i][1]: " << e1[i][1] << " e1[i][2]: " << e1[i][2] << endl;

              writeTriTriBin(mpi_rank,x0,x1);
              writeTriTri(mpi_rank,x0,x1);
              assert(0);

            }

          }

        }
      }
    }

    if (caseSet.find(value) == caseSet.end()) {

      cout << "WARNING: value: " << value << " was not found: triTriIntVec.size(): " << triTriIntVec.size() << endl;

      FOR_I3 cout << "p0: " << p0[i] << " " << p0[(i+1)%3] << " e0" << i << " " << COUT_VEC(e0[i]) << endl;
      FOR_I3 cout << "p1: " << p1[i] << " " << p1[(i+1)%3] << " e1" << i << " " << COUT_VEC(e1[i]) << endl;

      caseSet.insert(value);

      cout << "triTriIntVec.size(): " << triTriIntVec.size() << endl;
      for (int ii = 0; ii < triTriIntVec.size(); ++ii)
        triTriIntVec[ii].dump(" > ");

      // look for duplicates in the triTriIntVec...

      for (int ii = 0; ii < triTriIntVec.size(); ++ii) {
        if (triTriIntVec[ii].kind != MATCHED_INT) { // skip any already matched!
          for (int jj = ii+1; jj < triTriIntVec.size(); ++jj) {
            if (triTriIntVec[ii] == triTriIntVec[jj]) {
              cout << "GOT DUPLICATE triTriIntVec: " << endl;
              triTriIntVec[ii].dump(" > ii: ");
              triTriIntVec[jj].dump(" > jj: ");
              // confirm that the type is edge-edge and the doubles exactly match...
              // when we hit edge-node intersections, we will revisit...
              ++triTriIntVec[ii].match_count;
              if (triTriIntVec[ii].kind == EDGE_EDGE_INT) {
                assert(triTriIntVec[ii].ddata[0] == triTriIntVec[jj].ddata[0]);
                assert(triTriIntVec[ii].ddata[1] == triTriIntVec[jj].ddata[1]);
                triTriIntVec[jj].kind = MATCHED_INT;
                triTriIntVec[jj].idata[0] = ii;
              }
              else if (triTriIntVec[ii].kind == NODE_NODE_INT) {
                triTriIntVec[jj].kind = MATCHED_INT;
                assert(triTriIntVec[ii].idata[0] == triTriIntVec[jj].idata[0]);
                assert(triTriIntVec[ii].idata[1] == triTriIntVec[jj].idata[1]);
              }
              else if (triTriIntVec[ii].kind == EDGE_NODE_INT) {
                triTriIntVec[jj].kind = MATCHED_INT;
                assert(triTriIntVec[ii].idata[0] == triTriIntVec[jj].idata[0]);
                assert(triTriIntVec[ii].idata[1] == triTriIntVec[jj].idata[1]);
                cout << "triTriIntVec[ii].ddata[0]: " << triTriIntVec[ii].ddata[0] << " triTriIntVec[jj].ddata[0]: " << triTriIntVec[jj].ddata[0] <<
                  " diff: " << triTriIntVec[ii].ddata[0]-triTriIntVec[jj].ddata[0] << endl;
                //assert(triTriIntVec[ii].ddata[0] == triTriIntVec[jj].ddata[0]);
              }
              else if (triTriIntVec[ii].kind == NODE_EDGE_INT) {
                triTriIntVec[jj].kind = MATCHED_INT;
                assert(triTriIntVec[ii].idata[0] == triTriIntVec[jj].idata[0]);
                assert(triTriIntVec[ii].idata[1] == triTriIntVec[jj].idata[1]);
                cout << "triTriIntVec[ii].ddata[0]: " << triTriIntVec[ii].ddata[0] << " triTriIntVec[jj].ddata[0]: " << triTriIntVec[jj].ddata[0] <<
                  " diff: " << triTriIntVec[ii].ddata[0]-triTriIntVec[jj].ddata[0] << endl;
                //assert(triTriIntVec[ii].ddata[0] == triTriIntVec[jj].ddata[0]);
              }
              else {
                assert(0);
              }
            }
          }
        }
      }

      int ivec[10];
      int count = 0;
      for (int ii = 0; ii < triTriIntVec.size(); ++ii) {
        if (triTriIntVec[ii].kind != MATCHED_INT) {
          if ((triTriIntVec[ii].kind == NODE_NODE_INT)&&(triTriIntVec[ii].match_count == 0)) {
            triTriIntVec[ii].dump("Warning: skipping inconsistent NODE_NODE_INT");
          }
          else {
            ivec[count++] = ii;
          }
        }
      }

      if (count == 0) {

        cout << "CODE:" << endl;
        cout << "CODE: case " << value << ":" << endl;
        cout << "CODE:   return 0;" << endl;

      }
      else if (count == 1) {

        cout << "ONE INTERSECTION" << endl;

        cout << "CODE:" << endl;
        cout << "CODE: case " << value << ":" << endl;
        cout << "CODE:   idata[0] = " << triTriIntVec[ivec[0]].kindString() << ";" << endl;
        if (triTriIntVec[ivec[0]].kind == EDGE_EDGE_INT) {
          cout << "CODE:   idata[1] = " << triTriIntVec[ivec[0]].idata[0] << "; // edge0 index on tri0" << endl;
          cout << "CODE:   idata[2] = " << triTriIntVec[ivec[0]].idata[1] << "; // edge1 index on tri1" << endl;
          cout << "CODE:   idata[3] = 0; // unused" << endl;
          cout << "CODE:   idata[4] = 0; // unused" << endl;
          cout << "CODE:   ddata[0] = " << triTriIntVec[ivec[0]].str[0] << "; // edge0 wgt" << endl;
          cout << "CODE:   ddata[1] = " << triTriIntVec[ivec[0]].str[1] << "; // edge1 wgt" << endl;
          cout << "CODE:   ddata[2] = 0.0; // unused" << endl;
          cout << "CODE:   ddata[3] = 0.0; // unused" << endl;
        }
        else if (triTriIntVec[ivec[0]].kind == NODE_NODE_INT) {
          cout << "CODE:   idata[1] = " << triTriIntVec[ivec[0]].idata[0] << "; // node0 index on tri0" << endl;
          cout << "CODE:   idata[2] = " << triTriIntVec[ivec[0]].idata[1] << "; // node1 index on tri1" << endl;
          cout << "CODE:   idata[3] = 0; // unused" << endl;
          cout << "CODE:   idata[4] = 0; // unused" << endl;
          cout << "CODE:   ddata[0] = 0.0; // unused" << endl;
          cout << "CODE:   ddata[1] = 0.0; // unused" << endl;
          cout << "CODE:   ddata[2] = 0.0; // unused" << endl;
          cout << "CODE:   ddata[3] = 0.0; // unused" << endl;
        }
        else if (triTriIntVec[ivec[0]].kind == EDGE_NODE_INT) {
          cout << "CODE:   idata[1] = " << triTriIntVec[ivec[0]].idata[0] << "; // edge0 index on tri0" << endl;
          cout << "CODE:   idata[2] = " << triTriIntVec[ivec[0]].idata[1] << "; // node1 index on tri1" << endl;
          cout << "CODE:   idata[3] = 0; // unused" << endl;
          cout << "CODE:   idata[4] = 0; // unused" << endl;
          cout << "CODE:   ddata[0] = " << triTriIntVec[ivec[0]].str[0] << "; // edge0 wgt" << endl;
          cout << "CODE:   ddata[1] = 0.0; // unused" << endl;
          cout << "CODE:   ddata[2] = 0.0; // unused" << endl;
          cout << "CODE:   ddata[3] = 0.0; // unused" << endl;
        }
        else if (triTriIntVec[ivec[0]].kind == NODE_EDGE_INT) {
          cout << "CODE:   idata[1] = " << triTriIntVec[ivec[0]].idata[0] << "; // node0 index on tri0" << endl;
          cout << "CODE:   idata[2] = " << triTriIntVec[ivec[0]].idata[1] << "; // edge1 index on tri1" << endl;
          cout << "CODE:   idata[3] = 0; // unused" << endl;
          cout << "CODE:   idata[4] = 0; // unused" << endl;
          cout << "CODE:   ddata[0] = " << triTriIntVec[ivec[0]].str[0] << "; // edge1 wgt" << endl;
          cout << "CODE:   ddata[1] = 0.0; // unused" << endl;
          cout << "CODE:   ddata[2] = 0.0; // unused" << endl;
          cout << "CODE:   ddata[3] = 0.0; // unused" << endl;
        }
        else {
          cout << "unhandled 1-intersection problem!" << endl;
          writeTriTriBin(mpi_rank,x0,x1);
          writeTriTri(mpi_rank,x0,x1);
          assert(0);
        }
        cout << "CODE:   return 1;" << endl;

      }
      else if (count == 2) {

        cout << "TWO INTERSECTIONS" << endl;

        double xi0[3];
        triTriIntVec[ivec[0]].getX(xi0,x0,x1);
        double xi1[3];
        triTriIntVec[ivec[1]].getX(xi1,x0,x1);

        //fp = fopen("int.dat","w");
        //writeEdge(fp,xi0,xi1);
        //fclose(fp);

        const double dxi[3] = DIFF(xi1,xi0);
        const double dxi_mag = MAG(dxi);
        const double normal0[3] = TRI_NORMAL_2(x0[0],x0[1],x0[2]);
        const double normal1[3] = TRI_NORMAL_2(x1[0],x1[1],x1[2]);
        const double cp[3] = CROSS_PRODUCT(normal0,normal1);
        const double cp_mag = MAG(cp);
        const double dp = DOT_PRODUCT(dxi,cp);
        if (dp > 0.0) {
          if (fabs(1.0-dp/(dxi_mag*cp_mag)) > 1.0E-3) {
            cout << "bad alignment!" << endl;
            assert(0);
          }
        }
        else if (dp < 0.0) {
          if (fabs(1.0+dp/(dxi_mag*cp_mag)) > 1.0E-3) {
            cout << "bad flipped alignment!" << endl;
            assert(0);
          }
          // flip order...
          const int tmp = ivec[0];
          ivec[0] = ivec[1];
          ivec[1] = tmp;
        }
        cout << "CODE:" << endl;
        cout << "CODE: case " << value << ":" << endl;
        int ii = 0;
        int id = 0;
        FOR_I2 {
          cout << "CODE:   // intersection " << i << "..." << endl;
          cout << "CODE:   idata["<<ii<<"] = " << triTriIntVec[ivec[i]].kindString() << ";" << endl;
          ++ii;
          if (triTriIntVec[ivec[i]].kind == EDGE_TRI_INT) {
            cout << "CODE:   idata["<<ii<<"] = " << triTriIntVec[ivec[i]].idata[0] << "; // edge index on tri0" << endl;
            ++ii;
            cout << "CODE:   idata["<<ii<<"] = 0; // unused" << endl;
            ++ii;
            cout << "CODE:   ddata["<<id<<"] = " << triTriIntVec[ivec[i]].str[0] << "; // edge wgt" << endl;
            ++id;
            cout << "CODE:   ddata["<<id<<"] = 0.0; // unused" << endl;
            ++id;
          }
          else if (triTriIntVec[ivec[i]].kind == TRI_EDGE_INT) {
            cout << "CODE:   idata["<<ii<<"] = " << triTriIntVec[ivec[i]].idata[0] << "; // edge index on tri1" << endl;
            ++ii;
            cout << "CODE:   idata["<<ii<<"] = 0; // unused" << endl;
            ++ii;
            cout << "CODE:   ddata["<<id<<"] = " << triTriIntVec[ivec[i]].str[0] << "; // edge wgt" << endl;
            ++id;
            cout << "CODE:   ddata["<<id<<"] = 0.0; // unused" << endl;
            ++id;
          }
          else if (triTriIntVec[ivec[i]].kind == EDGE_EDGE_INT) {
            cout << "CODE:   idata["<<ii<<"] = " << triTriIntVec[ivec[i]].idata[0] << "; // edge0 index on tri0" << endl;
            ++ii;
            cout << "CODE:   idata["<<ii<<"] = " << triTriIntVec[ivec[i]].idata[1] << "; // edge1 index on tri1" << endl;
            ++ii;
            cout << "CODE:   ddata["<<id<<"] = " << triTriIntVec[ivec[i]].str[0] << "; // edge0 wgt" << endl;
            ++id;
            cout << "CODE:   ddata["<<id<<"] = " << triTriIntVec[ivec[i]].str[1] << "; // edge1 wgt" << endl;
            ++id;
          }
          else if (triTriIntVec[ivec[i]].kind == NODE_NODE_INT) {
            cout << "CODE:   idata["<<ii<<"] = " << triTriIntVec[ivec[i]].idata[0] << "; // node0 index on tri0" << endl;
            ++ii;
            cout << "CODE:   idata["<<ii<<"] = " << triTriIntVec[ivec[i]].idata[1] << "; // node1 index on tri1" << endl;
            ++ii;
            cout << "CODE:   ddata["<<id<<"] = 0.0; // unused" << endl;
            ++id;
            cout << "CODE:   ddata["<<id<<"] = 0.0; // unused" << endl;
            ++id;
          }
          else if (triTriIntVec[ivec[i]].kind == EDGE_NODE_INT) {
            cout << "CODE:   idata["<<ii<<"] = " << triTriIntVec[ivec[i]].idata[0] << "; // edge0 index on tri0" << endl;
            ++ii;
            cout << "CODE:   idata["<<ii<<"] = " << triTriIntVec[ivec[i]].idata[1] << "; // node1 index on tri1" << endl;
            ++ii;
            cout << "CODE:   ddata["<<id<<"] = " << triTriIntVec[ivec[i]].str[0] << "; // edge0 wgt" << endl;
            ++id;
            cout << "CODE:   ddata["<<id<<"] = 0.0; // unused" << endl;
            ++id;
          }
          else if (triTriIntVec[ivec[i]].kind == NODE_EDGE_INT) {
            cout << "CODE:   idata["<<ii<<"] = " << triTriIntVec[ivec[i]].idata[0] << "; // node0 index on tri0" << endl;
            ++ii;
            cout << "CODE:   idata["<<ii<<"] = " << triTriIntVec[ivec[i]].idata[1] << "; // edge1 index on tri1" << endl;
            ++ii;
            cout << "CODE:   ddata["<<id<<"] = " << triTriIntVec[ivec[i]].str[0] << "; // edge1 wgt" << endl;
            ++id;
            cout << "CODE:   ddata["<<id<<"] = 0.0; // unused" << endl;
            ++id;
          }
          else {
            cout << "unhandled 2-intersection problem!" << endl;
            writeTriTriBin(mpi_rank,x0,x1);
            writeTriTri(mpi_rank,x0,x1);
            assert(0);
          }
        }
        cout << "CODE:   return 2;" << endl;

      }
      else {
        cout << "unhandled: triTriIntVec.size(): " << triTriIntVec.size() << endl;
        FILE * fp = fopen("xint.dat","w");
        for (int ii = 0; ii < triTriIntVec.size(); ++ii) {
          if (triTriIntVec[ii].kind != MATCHED_INT) { // skip any already matched!
            double xint[3];
            triTriIntVec[ii].getX(xint,x0,x1);
            fprintf(fp,"%18.15e %18.15e %18.15e\n",xint[0],xint[1],xint[2]);
          }
        }
        fclose(fp);
        writeTriTriBin(mpi_rank,x0,x1);
        writeTriTri(mpi_rank,x0,x1);
        assert(0);
      }

      cout << "take a look!" << endl;
      writeTriTriBin(mpi_rank,x0,x1);
      writeTriTri(mpi_rank,x0,x1);
      //getchar();

    }

    // means we do not have an entry...
    return -2;

  }

  void checkTriTri() {

    // randomly generate tris and run the above TriTriIntersection routine...

    cout << "checkTriTri()" << endl;

    double x0[3][3];
    double x1[3][3];
    int idata[6];
    double ddata[4];

    int8 iter = 0;
    while (1) {

      ++iter;
      if (iter%10000 == 0)
        cout << " > done " << iter << " tri pairs. New intersections discovered: " << caseSet.size() << endl;

      FOR_I3 FOR_J3 x0[i][j] = double(rand())/double(RAND_MAX)-0.5;
      FOR_I3 FOR_J3 x1[i][j] = double(rand())/double(RAND_MAX)-0.5;

      double * x0p[3]; FOR_I3 x0p[i] = x0[i];
      double * x1p[3]; FOR_I3 x1p[i] = x1[i];

      const int ierr = calcTriTriIntersection(idata,ddata,x0p,x1p);
      if (ierr == -2) {
        // run the suite of permutations...
        FOR_I3 FOR_J3 {
          x0p[0] = x0[i];
          x0p[1] = x0[(i+1)%3];
          x0p[2] = x0[(i+2)%3];
          x1p[0] = x1[j];
          x1p[1] = x1[(j+1)%3];
          x1p[2] = x1[(j+2)%3];
          calcTriTriIntersection(idata,ddata,x0p,x1p);
          x1p[0] = x0[i];
          x1p[1] = x0[(i+1)%3];
          x1p[2] = x0[(i+2)%3];
          x0p[0] = x1[j];
          x0p[1] = x1[(j+1)%3];
          x0p[2] = x1[(j+2)%3];
          calcTriTriIntersection(idata,ddata,x0p,x1p);
          x0p[0] = x0[(i+2)%3];
          x0p[1] = x0[(i+1)%3];
          x0p[2] = x0[i];
          x1p[0] = x1[j];
          x1p[1] = x1[(j+1)%3];
          x1p[2] = x1[(j+2)%3];
          calcTriTriIntersection(idata,ddata,x0p,x1p);
          x0p[0] = x0[(i+2)%3];
          x0p[1] = x0[(i+1)%3];
          x0p[2] = x0[i];
          x1p[0] = x1[(j+2)%3];
          x1p[1] = x1[(j+1)%3];
          x1p[2] = x1[j];
          calcTriTriIntersection(idata,ddata,x0p,x1p);
          x1p[0] = x0[(i+2)%3];
          x1p[1] = x0[(i+1)%3];
          x1p[2] = x0[i];
          x0p[0] = x1[j];
          x0p[1] = x1[(j+1)%3];
          x0p[2] = x1[(j+2)%3];
          calcTriTriIntersection(idata,ddata,x0p,x1p);
          x1p[0] = x0[(i+2)%3];
          x1p[1] = x0[(i+1)%3];
          x1p[2] = x0[i];
          x0p[0] = x1[(j+2)%3];
          x0p[1] = x1[(j+1)%3];
          x0p[2] = x1[j];
          calcTriTriIntersection(idata,ddata,x0p,x1p);
        }

      }

    }

  }

}
