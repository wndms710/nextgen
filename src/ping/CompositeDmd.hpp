
#include "DmdNew.hpp"
#include <complex>
#include <vector>

using namespace std;

class SnapshotData { 
public: 

  // define the snapshot inner product

  virtual double dotp(const int i, const int j) = 0;

  // create a new snapshot (complex mode) .. 

  virtual void create_and_zero_mode(const string& name) = 0;

  // add to a named mode shape with snapshot "i" 

  virtual void axpy(const string& name, const double c_r, const double c_i, const int i) = 0;
  
  void axpy(const string& name, const complex<double> coeff, const int i) { 
    
    axpy(name,coeff.real(),coeff.imag(),i);

  }

  // mode normalization 

  virtual void normalize_mode(const string& name) = 0;

  // write out the mode in some unspecified format -- 
  // passing an additional parameter to allow finding some dummy information as necessary..

  virtual void write_mode(const string& mode_name, const int itmp) = 0;

  // cache blocking for snapshot data 

  virtual int getBlocksize() const = 0;
};


class CompositeDmd { 

public: 

  // we can build the snapshot correlation matrix against any unique 
  // integer hashes, but for now we will use the step indices... 

  int istart;
  int istep; 
  int iend;
  int ns;
  double * A_corr;
  double dt_samp;
  

  // dmd output --- corresponding to the forward Koopman operator 

  int nnz;
  double * lambda_r; 
  double * lambda_i;
  complex<double>* vr_imag;
  double * Wk;
  double * amp;
 
  // dmd output -- corresponding to the adjoint Koopman operator 

  double * lambda_r_adj; 
  double * lambda_i_adj;
  complex<double>* vrh_imag;

  SnapshotData* data;

  // frequencies of interest to find ... 

  vector<double> mode_freq;

  CompositeDmd(const int _istart, const int _istep, const int _iend, SnapshotData* _data, 
               const double _dt_samp) { 

    istart  = _istart;
    istep   = _istep ;
    iend    = _iend;
    data    = _data;
    dt_samp = _dt_samp;

    assert((iend - istart)%istep == 0);
    ns      = (iend - istart)/istep + 1; // inclusive at the end ... 
    A_corr  = NULL;

    assert( data != NULL);

    nnz      = -1;
    lambda_r = NULL;
    lambda_i = NULL;
    vr_imag  = NULL;
    Wk       = NULL;
    amp      = NULL;

    // the eigenvalues of the adjoint are redundant -- but we can keep 
    // them for double checking that everything is consistently done ... 

    lambda_r_adj = NULL;
    lambda_i_adj = NULL;
    vrh_imag     = NULL;

  }

  ~CompositeDmd() { 
   
    // its the responsibility of the calling process to delete its 
    // snapshot data manager .. 

    if ( lambda_r != NULL) delete[] lambda_r;
    if ( lambda_i != NULL) delete[] lambda_i;
    if ( vr_imag != NULL)  delete[] vr_imag;
    if ( Wk != NULL)       delete[] Wk;
    if ( amp != NULL)      delete[] amp;

    if ( lambda_r_adj != NULL) delete[] lambda_r_adj;
    if ( lambda_i_adj != NULL) delete[] lambda_i_adj;
    if ( vrh_imag     != NULL) delete[] vrh_imag;

    if ( A_corr != NULL)   delete[] A_corr;

  }

  void writeMatrix(const char* filename, const double* A) {
    int nd = 0;
    FILE * fp = fopen(filename, "wb");
    fwrite(&ns, sizeof(int), 1, fp); 
    fwrite(&nd, sizeof(int), 1, fp);
    fwrite(A, sizeof(double), ns*ns, fp); 
    fclose(fp); 
  }//writeMatrix 


  void build_correlation_matrix() { 

    assert( A_corr == NULL);
    A_corr = new double[ns*ns];

    if ( FILE* fp = fopen("A.dat", "rb") ) { 

      if ( mpi_rank == 0 )
        cout << " found a cached correlation matrix file " << endl;

      // found the matrix ... we'll read from the file
      int c_ns, c_nd;
      fread(&c_ns,sizeof(int),1,fp); assert( c_ns == ns);
      fread(&c_nd,sizeof(int),1,fp); assert( c_nd == 0);
      fread(A_corr, sizeof(double), ns*ns, fp);
      fclose(fp);

      const double dp = data->dotp(istart,istart);

    } else { 


      // note that the correlation matrix is indeed symmetric so we'll populate 
      // the lower triangle.  storing the matrix in col major...

      if ( mpi_rank == 0 ) { 
        cout << " building correlation matrix ... ";
        cout.flush();
      }

      /*
      for (int j = 0; j < ns; ++j) {
        
        cout << "." ; cout.flush();
        
        for (int i = 0; i <= j; ++i) { 
          
          const int ii = istart + istep*i;
          const int jj = istart + istep*j;
          
          A_corr[(j)*ns+i] = data->dotp(ii,jj);
          
          // copy down the value into the upper triangle .. 
          
          if ( i != j) { 
            A_corr[(i)*ns+j] = A_corr[(j)*ns+i];
          } 
        }
      }
      */

      
      double * my_A = new double[ns*ns];

      for (int i = 0; i < ns*ns; ++i) { 
        my_A[i] = 0.0;
        A_corr[i] = 0.0;
      }

      // do some manual cache blocking -- this helps if each dot product
      // needs to load snapshot data in .. 

      const int bs = max(data->getBlocksize(),2);

      vector<Block*> blocks; 
      
      for (int j=0; j < ns ; j += bs ) { 
        for (int i =j; i < ns ; i += bs ) { 
          
          Block* b  = new Block(); 
          b->istart = i; 
          b->iend   = min(i+bs, ns); // NOT inclusive of the end.
          b->jstart = j;
          b->jend   = min(j+bs, ns);
          b->rank   = -1;
          blocks.push_back(b);
        }
      }

      // in parallel .. round-robin the distribution of the blocks .. 
    
      for (int j=0; j < blocks.size(); ++j ) { 
        int r = j%mpi_size; 
        blocks[j]->rank = r;
      }
      
      for (int ibl = 0, lim = blocks.size(); ibl < lim; ++ibl) { 

        if ( blocks[ibl]->rank == mpi_rank) { 

          Block * b = blocks[ibl]; 
          if ( mpi_rank == 0 ) {  
            cout << "." ; 
            cout.flush(); 
          }
          
          for (int j = b->jstart; j < b->jend; ++j) { 
            for (int i = b->istart; i < b->iend ; ++i) {
              
              const int ii = istart + istep*i;
              const int jj = istart + istep*j;
              
              my_A[(j)*ns+i] = data->dotp(ii,jj);
              
            }
          }
        }
      }
      
      
      MPI_Allreduce(my_A, A_corr, ns*ns, MPI_DOUBLE, MPI_SUM, mpi_comm);
     
      // fill out the lower triangle..
     
      for (int j=0; j < ns ; ++j) { 
        for (int i=0; i < j ; ++i) { 
          A_corr[(j)*ns+i] = A_corr[(i)*ns+j];
        }
      }

      for (int ibl=0; ibl < blocks.size(); ++ibl) 
        delete blocks[ibl];
      
      delete[] my_A;

      if ( mpi_rank == 0 ) { 
        cout << " ok. " << endl;
        writeMatrix("A.dat", A_corr);
      }
    }

  } 
 
  void extract_mode(const int imo) { 

    //complex<double> tt(lambda_r[imo],lambda_i[imo]);
    //const double freq  =  arg(tt)/(2.0*M_PI*dt_samp);
    
    const double freq = lambda_i[imo];

    string mode_name; 
    { 
      char name[128]; sprintf(name, "mode.%g", freq);
      mode_name = string(name);
    }
    
    
    assert( Wk != NULL);
    assert( vr_imag != NULL);
    
    data->create_and_zero_mode(mode_name);
    
    const int nt = ns -1;
    
    for (int i = 0; i < nt; ++i) { 
      
      double a_r = 0.0;
      double a_i = 0.0;
      
      for (int j = 0; j < nnz; ++j) { 
        a_r += Wk[(j)*nt + i]* vr_imag[(imo)*nnz+j].real();
        a_i += Wk[(j)*nt + i]* vr_imag[(imo)*nnz+j].imag();
      }
      
      // these are a linear combination of the Y snapshots ... see notes
      
      const int ii = istart + istep*(i+1); assert( ii <= iend);
      data->axpy(mode_name,a_r,a_i,ii);
      
    }
    
    data->normalize_mode(mode_name);
    
    data->write_mode(mode_name,istart);
    
  } 

  void extract_adjoint_mode(const int imo) { 
    
    complex<double> tt(lambda_r_adj[imo],lambda_i_adj[imo]);
    const double freq  =  arg(tt)/(2.0*M_PI*dt_samp);
    
    string adj_mode_name; 
    { 
      char name[128]; sprintf(name, "adj_mode.%g", freq);
      adj_mode_name = string(name);
    }
    
    assert( Wk != NULL);
    assert( vrh_imag != NULL);
    
    data->create_and_zero_mode(adj_mode_name);
    
    const int nt = ns -1;
    
    for (int i = 0; i < nt; ++i) { 
      
      double a_r = 0.0;
      double a_i = 0.0;
      
      for (int j = 0; j < nnz; ++j) { 
        a_r += Wk[(j)*nt + i]* vrh_imag[(imo)*nnz+j].real();
        a_i += Wk[(j)*nt + i]* vrh_imag[(imo)*nnz+j].imag();
      }
      
      // these are a linear combination of the X snapshots ... see notes
      
      const int ii = istart + istep*(i); assert( ii < iend);
      data->axpy(adj_mode_name,a_r,a_i,ii);
      
    }
    
    data->normalize_mode(adj_mode_name);
    
    data->write_mode(adj_mode_name,istart);
  } 
 
  int findMode(const double& freq, const double* eig_r, const double* eig_i, const bool logscale = true) { 

    double min_dist = 1.0e+16;
    int ii          = -1;

    for (int i = 0; i < nnz; ++i) { 

      double f;
      if ( logscale ) { 

        f = eig_i[i];

      } else { 

        complex<double> tt(eig_r[i],eig_i[i]);
        f      = arg(tt)/(2.0*M_PI*dt_samp);
      }
      
      double this_dist    = abs(f - freq);
      if ( this_dist < min_dist) { 
        ii = i;
        min_dist = this_dist;
      }
    }

    if ( logscale) { 
      cout << " found freq: " << eig_i[ii] << endl;

    } else { 

      complex<double> tt(eig_r[ii],eig_i[ii]);
      double f      = arg(tt)/(2.0*M_PI*dt_samp);
      cout << " found not-lgscale freq: " << f << endl;

    }


    assert( (ii >= 0)&&(ii < nnz));
    return ii;
  }

  void run() { 

    build_correlation_matrix();

    // for now -- aside from the correlation matrix, everything else is done in serial

    if ( mpi_rank == 0 ) { 

      cout << " computing dmd decomposition .. " << endl;
      
      assert( dt_samp > 0.0);
      nnz = compute_dmd(ns,A_corr,lambda_r,lambda_i,Wk,vr_imag,amp,dt_samp);
      
      assert ( amp == NULL);
     
      compute_dmd_amp(amp, ns, A_corr, Wk, lambda_r, lambda_i, vr_imag, dt_samp, nnz);

      report_amplitude(amp);

      delete[] amp; amp = NULL;

      // compute the adjoint information --- this is just one more eigenvalue solve ..
      
      compute_dmd_adjoint(ns,A_corr,Wk,lambda_r_adj,lambda_i_adj,vrh_imag,dt_samp,nnz);
     

      // compute the amplitudes of the adjoint modes .. 

      double * adj_amp = NULL; 

      compute_dmd_adjoint_amp(adj_amp, ns, A_corr, Wk, lambda_r_adj, 
                              lambda_i_adj, vrh_imag, dt_samp, nnz);

      report_adj_amplitude(adj_amp);

      delete[] adj_amp;

      // extract the modes and their adjoints ... calling process needs to 
      // have set the frequency variables at some point .. 
      
      //mode_freq.push_back(215.0);
      mode_freq.push_back(5467.0);
      mode_freq.push_back(186.0);
      mode_freq.push_back(1488.0);
 
      for (int ii = 0, lim = mode_freq.size(); ii < lim; ++ii) { 
        
        cout << " working on mode at freq : " << mode_freq[ii] << endl;
        
        // extract the forward mode .. 
        
        const int imo_f = findMode(mode_freq[ii],lambda_r,lambda_i,true); 
        extract_mode(imo_f);
        
        // extract the adjoint mode ... 
        
        const int imo_a = findMode(mode_freq[ii],lambda_r_adj,lambda_i_adj,false);
        extract_adjoint_mode(imo_a);
        
      }
      
      cout << " done ... " << endl;
    
    }
  
    MPI_Barrier(mpi_comm);
  
  }

  void report_amplitude(double* amp) { 

    vector<pair<double,double> > zipped;

    for (int i = 0; i < nnz; ++i) {

      // these eigenvalues are already in lgoscale format... 

      zipped.push_back(pair<double,double>(amp[i],lambda_i[i]));
    }
    
    sort(zipped.begin(), zipped.end());

    cout << " --- reporting dmd amplitudes --- " << endl;

    for (int i = 0; i < nnz; ++i) { 

      cout << "ii, freq, amp: " << i << "   " << zipped[i].second << "    " << zipped[i].first << endl;
    }

  } 


  void report_adj_amplitude(double* adj_amp) { 

    vector<pair<double,double> > zipped;

    for (int i = 0; i < nnz; ++i) {

      complex<double> tt(lambda_r_adj[i],lambda_i_adj[i]);
      const double freq  =  arg(tt)/(2.0*M_PI*dt_samp);
      zipped.push_back(pair<double,double>(adj_amp[i],freq));
    }
    
    sort(zipped.begin(), zipped.end());

    cout << " --- reporting adj amplitudes --- " << endl;

    for (int i = 0; i < nnz; ++i) { 

      cout << "i, freq, adj_amp: " << i << "   " << zipped[i].second << "    " << zipped[i].first << endl;
    }

  } 
};

