
#include "PngAnalysis.hpp"
#include "CompositeDmd.hpp"

template<int n>
float dp_aux2(const float *__restrict__ t1, const float *__restrict__ t2) { 
  
  float sum = 0.0;
  for (int i = 0; i <n; ++i) 
    sum += t1[i]*t2[i];
  
  return sum;
}

class ImageCache { 

public:
  int ncache;
  deque<PngData*> image_cache;
  string prefix;

  ImageCache(const string& img_prefix) { 

    ncache = 15;
    image_cache.clear();
    prefix = img_prefix;

  }

  ~ImageCache() { 
    for (deque<PngData*>::iterator it = image_cache.begin(); it != image_cache.end(); ++it) 
      delete *it;
  }

  PngData* getImage(const int ii) { 

    deque<PngData*>::iterator it = find_if(image_cache.begin(), image_cache.end(), PngDataComp(ii)); 
    
    if ( it != image_cache.end() ) { 
      (*it)->lock = true;
      return *it; 
    } else {

      if ( image_cache.size() == ncache) { 
        // cache is full, remove the first unlocked element (first in-first out)
        for(it = image_cache.begin(); it != image_cache.end(); ++it) { 
          if (!((*it)->lock)) 
            break;
        }//it

        assert( it != image_cache.end()); // uh oh, the whole cache is locked...
        delete *it; 
        image_cache.erase(it); 
      }
     
      // when we get the image, we will lock it, but 
      // the user is responsible for releasing the lock
      
      int crop_start[2] = {0, 0};
      int crop_end[2]   = {-1, -1};
      
      PngData* png = new PngData(crop_start, crop_end);

      char filename[128]; sprintf(filename, "%s.%08d.png", prefix.c_str(), ii);
      png->read(filename);
      png->finalizeRead();
      png->id   = ii;
      png->lock = true; 
      image_cache.push_back(png);

      assert( image_cache.size() <= ncache); 
      
      return png;
    }
    
    // shouldnt get here..
    return NULL; 
  }
};


class ComplexModeM { 
public: 

  int n; 
  int m;
  vector<double*> phi_vec;

  ComplexModeM(const int _n, const int _m) { 

    n = _n;
    m = _m;

    // we need 2m entries to account for the real
    // and imaginary components -- otherwise store
    // complex data types... 

    for (int i = 0; i < 2*m; ++i) { 

      double * phi = new double[n];
      for (int i =0; i < n; ++i) 
        phi[i] = 0.0;

      phi_vec.push_back(phi);
    }

  }

  ~ComplexModeM() { 

    for (int i = 0; i < 2*m; ++i) 
      delete[] phi_vec[i];

  }
};


class ComplexMode2 { 
public: 

  int n;
  double * p_r;
  double * p_i;
  double * z_r;
  double * z_i;

  ComplexMode2(const int _n) { 

    n   = _n;
    p_r = new double[n];
    p_i = new double[n];
    z_r = new double[n];
    z_i = new double[n];

    for (int i = 0; i < n; ++i) { 
      p_r[i] = 0.0;
      p_i[i] = 0.0;
      z_r[i] = 0.0;
      z_i[i] = 0.0;
    }

  }

  ~ComplexMode2() { 

    delete[] p_r;
    delete[] p_i;
    delete[] z_r;
    delete[] z_i;

  }
};
    

class ImageSnapshotM : public SnapshotData { 
public: 

  int m;
  vector<ImageCache*> caches;
  vector<string> data_names;
  vector<double> d_ref;
  int n_size;

  map<string,ComplexModeM*> mode_map;

  ImageSnapshotM(vector<string>& _data_names, vector<string>& _img_prefix, vector<double>& _d_ref) { 

    assert( _data_names.size() == _d_ref.size());
    assert( _data_names.size() == _img_prefix.size());
    
    m = int(_data_names.size());

    for (int i = 0; i < m; ++i) { 

      caches.push_back(new ImageCache(_img_prefix[i]));
      data_names.push_back(_data_names[i]);
      d_ref.push_back(_d_ref[i]);

    }

    n_size = -1;
    mode_map.clear();
  }


  ~ImageSnapshotM() { 

    for (int i =0; i < m; ++i) 
      delete caches[i];

    for (map<string,ComplexModeM*>::iterator it = mode_map.begin(); it != mode_map.end(); ++it) 
      delete it->second;

  }

  int getBlocksize() const { 

    // assumes that all of the caches have the same size ... 

    return caches[0]->ncache/2;
  }

  // because the data buffers are float -- in order to ensure there is sufficient precision
  // we are chunking the dot product..

  double dp_aux(const float *__restrict__ p1, const float *__restrict__ p2, const int nn) { 

    double sum           = 0.0;
    const int chunk_size = 512;

    for (int i = 0; i < nn; i += chunk_size) { 

      if ( nn-i > chunk_size ) { 

        float tmp_sum = dp_aux2<chunk_size>(&(p1[i]),&(p2[i]));
        sum += double(tmp_sum);

      } else { 

        // cleanup the remainder manually ... 

        for (int j = i; j < nn; ++j) { 
          sum += double(p1[i])*double(p2[i]);
        }
      }
    }

    return sum;
  }

  double dotp(const int i, const int j) { 

    double sum = 0.0;
    int npx    = -1;

    for (int ivar = 0; ivar < m; ++ivar) { 

      PngData * png_i = caches[ivar]->getImage(i);
      PngData * png_j = caches[ivar]->getImage(j);

      assert( png_i && png_j);

      const int nn = png_i->getNx() * png_j->getNy();

      if ( n_size == -1) { 
        n_size = nn;
      } else { 
        assert( n_size == nn);
      }

      if ( npx == -1) { 
        npx = png_i->npx;
        assert( npx == png_j->npx);
      } else { 
        assert( npx == png_i->npx);
        assert( npx == png_j->npx);
      }

      assert( nn == (png_j->getNx() * png_j->getNy()));

      const double tmp_sum = dp_aux(png_i->data,png_j->data, nn)/(d_ref[ivar]*d_ref[ivar]);
      png_i->lock = false;
      png_j->lock = false;

      sum += tmp_sum;
    }
      
    assert( npx > 0);
    return sum/double(npx);
  }


  void create_and_zero_mode(const string& name) { 

    map<string,ComplexModeM*>::iterator it = mode_map.find(name);

    if ( it == mode_map.end()) { 

      assert( n_size > 0); // should have been set; currently done in the dot product ...

      ComplexModeM * c = new ComplexModeM(n_size,m);
      mode_map[name]   = c;

    } else { 

      assert(0); // this mode already exists -- something is wrong

    }
  }


  void axpy(const string& name, const double c_r, const double c_i, const int ii) { 


    ComplexModeM* c = mode_map[name]; assert(c);
    assert( n_size > 0);

    for (int ivar = 0; ivar < m; ++ivar) { 

      PngData * png = caches[ivar]->getImage(ii); assert(png);

      for (int i = 0; i < n_size; ++i) { 

        c->phi_vec[2*ivar+0][i] += c_r * double(png->data[i])/d_ref[ivar];
        c->phi_vec[2*ivar+1][i] += c_i * double(png->data[i])/d_ref[ivar];

      }

      png->lock = false;
    }

  }

  void normalize_mode(const string& name) { 

    ComplexModeM* c = mode_map[name]; assert( c);
    double norm     = 0.0;

    assert( n_size > 0);

    for (int ivar = 0; ivar < m; ++ivar) { 

      for (int i = 0; i < n_size; ++i) { 
        
        norm += c->phi_vec[2*ivar+0][i] * c->phi_vec[2*ivar+0][i];
        norm += c->phi_vec[2*ivar+1][i] * c->phi_vec[2*ivar+1][i];

      }
    }

    norm = 1.0/sqrt(norm);

    for (int ivar = 0; ivar < m; ++ivar) { 

      for (int i = 0; i < n_size; ++i) { 

        c->phi_vec[2*ivar+0][i] *= norm;
        c->phi_vec[2*ivar+1][i] *= norm;
      }
    }
  }

  template<typename T>
  void dumpRange(T range[2], const T * buf, const int nn, const char* msg) { 
    
    range[0] = 1.0e+16; 
    range[1] = -1.0e+16; 
    
    for (int i =0; i < nn ; ++i) { 
      range[0] = min(range[0], buf[i]); 
      range[1] = max(range[1], buf[i]); 
    }
    
    cout << "dumpRange, " << msg << " = " << range[0] << ": " << range[1] << endl;
  }
  

  void write_mode(const string& mode_name, const int ii) { 

    ComplexModeM* c = mode_map[mode_name]; assert(c);

    for (int ivar = 0; ivar < m; ++ivar) { 

      char filename[128];

      // XXX why arent these in a default constructor somewhere?
      
      int crop_start[2] = {0, 0};
      int crop_end[2]   = {-1, -1};
      
      PngData * png = new PngData(crop_start, crop_end);
      
      
      sprintf(filename, "%s.%08d.png", caches[ivar]->prefix.c_str(), ii);
      png->read(filename);
      png->finalizeRead();

      double range[2];
      range[0] = 1.0e+16;
      range[1] = -1.0e+16;

      for (int i = 0; i < n_size; ++i) { 

        range[0] = fmin(range[0], c->phi_vec[2*ivar+0][i]);
        range[0] = fmin(range[0], c->phi_vec[2*ivar+1][i]);

        range[1] = fmax(range[1], c->phi_vec[2*ivar+0][i]);
        range[1] = fmax(range[1], c->phi_vec[2*ivar+1][i]);

      }

      // the real part of this variable .. 

      sprintf(filename, "%s_%s.real.png", mode_name.c_str(), data_names[ivar].c_str());
      for (int i = 0; i < n_size; ++i) { 
        png->data[i] = float(c->phi_vec[2*ivar+0][i]);
      }

      png->initializeWrite(range[0],range[1]);
      png->write(filename);


      // the imaginary part of this variable ... 

      sprintf(filename, "%s_%s.imag.png", mode_name.c_str(), data_names[ivar].c_str());
      for (int i = 0; i < n_size; ++i) { 
        png->data[i] = float(c->phi_vec[2*ivar+1][i]);
      }

      png->initializeWrite(range[0],range[1]);
      png->write(filename);

      // magnitude of this variable .. 

      sprintf(filename, "%s_%s.mag.png", mode_name.c_str(), data_names[ivar].c_str());
      for (int i = 0; i < n_size; ++i) { 
        png->data[i] = float(sqrt( c->phi_vec[2*ivar+0][i]*c->phi_vec[2*ivar+0][i] + 
                                   c->phi_vec[2*ivar+1][i]*c->phi_vec[2*ivar+1][i]   ));
      }

      double mag_range[2] = {0.0, sqrt(range[0]*range[0] + range[1]*range[1])};
      png->initializeWrite(mag_range[0], mag_range[1]);
      png->write(filename);

      delete png;
    }
  }
};

class ImageSnapshot2 : public SnapshotData { 
public:

  // this stores a double stacked variable based on p and Z ... 

  ImageCache * p_cache;
  ImageCache * Z_cache;

  double p_ref, Z_ref;
  int n_size;

  map<string,ComplexMode2*> mode_map;

  ImageSnapshot2() { 

    p_cache = new ImageCache("./images/p");
    Z_cache = new ImageCache("./images/z");

    // XXX hack 
    p_ref = 8.0e5;
    Z_ref = 0.1;

    n_size = -1;
    mode_map.clear();
  }

  ~ImageSnapshot2() { 

    delete p_cache;
    delete Z_cache;
    
    for (map<string,ComplexMode2*>::iterator it = mode_map.begin(); it != mode_map.end(); ++it) 
      delete it->second;

  }


  int getBlocksize() const { 
    return p_cache->ncache/2;
  }

  double dp_aux(const float *__restrict__ p1, const float *__restrict__ p2, const int nn) { 

    double sum           = 0.0;
    const int chunk_size = 512;

    for (int i = 0; i < nn; i += chunk_size) { 

      if ( nn-i > chunk_size ) { 

        float tmp_sum = dp_aux2<chunk_size>(&(p1[i]),&(p2[i]));
        sum += double(tmp_sum);

      } else { 

        // cleanup the remainder manually ... 

        for (int j = i; j < nn; ++j) { 
          sum += double(p1[i])*double(p2[i]);
        }
      }
    }

    return sum;
  }


  double dotp(const int i, const int j) { 

    PngData * ppng_i = p_cache->getImage(i);
    PngData * zpng_i = Z_cache->getImage(i);

    PngData * ppng_j = p_cache->getImage(j);
    PngData * zpng_j = Z_cache->getImage(j);

    assert( ppng_i && ppng_j);
    assert( zpng_j && zpng_j);

    // recall that the data buffer is equal to zero for non-data locations

    const int nn     = ppng_i->getNx() * ppng_i->getNy();
    
    if ( n_size == -1) { 
      n_size = nn; 
    } else { 
      assert( n_size == nn);
    }

    assert( nn == (ppng_j->getNx() * ppng_j->getNy()));
    assert( nn == (zpng_i->getNx() * zpng_i->getNy()));
    assert( nn == (zpng_j->getNx() * zpng_j->getNy()));

    const double sum1 = dp_aux(ppng_i->data, ppng_j->data, nn)/(p_ref*p_ref);
    const double sum2 = dp_aux(zpng_i->data, zpng_j->data, nn)/(Z_ref*Z_ref);

    ppng_i->lock = false;
    ppng_j->lock = false;
    zpng_i->lock = false;
    zpng_j->lock = false;

    return (sum1 + sum2)/(double(ppng_i->npx));
  }


  // create a new snapshot (complex mode) .. 

  void create_and_zero_mode(const string& name) { 

    map<string,ComplexMode2*>::iterator it = mode_map.find(name);

    if ( it == mode_map.end()) { 

      // currently n_size is being set in the dot product -- but should be done earlier... 

      ComplexMode2 * c = new ComplexMode2(n_size);
      mode_map[name]  = c;

    } else { 

      assert(0); // this mode already exists -- something is wrong.. 

    }
  }

  // add to a named mode shape with snapshot "i" 

  void axpy(const string& name, const double c_r, const double c_i, const int ii) { 

    ComplexMode2 * c = mode_map[name]; assert(c);
    PngData * ppng   = p_cache->getImage(ii); assert(ppng);
    PngData * zpng   = Z_cache->getImage(ii); assert(zpng);

    assert( n_size == ppng->getNx() * ppng->getNy());
    assert( n_size == zpng->getNx() * zpng->getNy());

    for (int i =0 ; i < n_size; ++i) { 
      c->p_r[i] += c_r * double(ppng->data[i])/p_ref;
      c->p_i[i] += c_i * double(ppng->data[i])/p_ref;
      c->z_r[i] += c_r * double(zpng->data[i])/Z_ref;
      c->z_i[i] += c_i * double(zpng->data[i])/Z_ref;
    }

    ppng->lock = false;
    zpng->lock = false;
  } 
  
  // mode normalization 

  void normalize_mode(const string& name) { 

    ComplexMode2 * c = mode_map[name]; assert(c);
    double norm      = 0.0;

    assert( n_size > 0);
    for (int i =0; i < n_size; ++i) { 
      norm += c->p_r[i]*c->p_r[i] + c->p_i[i]*c->p_i[i];
      norm += c->z_r[i]*c->z_r[i] + c->z_i[i]*c->z_i[i];
    }

    norm = 1.0/sqrt(norm);

    for (int i = 0; i < n_size; ++i) { 
      c->p_r[i] *= norm;
      c->p_i[i] *= norm;
      c->z_r[i] *= norm;
      c->z_i[i] *= norm;
    }
  }

  template<typename T>
  void dumpRange(T range[2], const T * buf, const int nn, const char* msg) { 

    range[0] = 1.0e+16; 
    range[1] = -1.0e+16; 
    
    for (int i =0; i < nn ; ++i) { 
      range[0] = min(range[0], buf[i]); 
      range[1] = max(range[1], buf[i]); 
    }
    
    cout << "dumpRange, " << msg << " = " << range[0] << ": " << range[1] << endl;
  }
  
  // write out the mode in some unspecified format 

  void write_mode(const string& mode_name, const int ii) { 

    // need to grab a finite index to help start the process -- stored in ii

    ComplexMode2* c = mode_map[mode_name]; assert(c);

    int crop_start[2] = {0, 0};
    int crop_end[2]   = {-1, -1};
 
    PngData* ppng = new PngData(crop_start, crop_end);
    PngData* zpng = new PngData(crop_start, crop_end);

    char filename[128]; 
    
    sprintf(filename, "%s.%08d.png", p_cache->prefix.c_str(), ii);
    ppng->read(filename);
    ppng->finalizeRead();

    sprintf(filename, "%s.%08d.png", Z_cache->prefix.c_str(), ii);
    zpng->read(filename);
    zpng->finalizeRead();

    double zrange[2];
    zrange[0] = 1.0e+16;
    zrange[1] = -1.0e+16;

    double prange[2];
    prange[0] = 1.0e+16;
    prange[1] = -1.0e+16;

    cout << " > n_size : "  << n_size << endl;

    for (int i =0; i < n_size; ++i) {

      zrange[0] = fmin(zrange[0], c->z_r[i]);
      zrange[0] = fmin(zrange[0], c->z_i[i]);
      zrange[1] = fmax(zrange[1], c->z_r[i]);
      zrange[1] = fmax(zrange[1], c->z_i[i]);

      prange[0] = fmin(prange[0], c->p_r[i]);
      prange[0] = fmin(prange[0], c->p_i[i]);
      prange[1] = fmax(prange[1], c->p_r[i]);
      prange[1] = fmax(prange[1], c->p_i[i]);

    }

    // safety in the ragen .. 

    
    cout << " p range : " << prange[0] << "   " << prange[1] << endl;
    cout << " z range : " << zrange[0] << "   " << zrange[1] << endl;



    // p real... 

    sprintf(filename,"%s_p.real.png",mode_name.c_str());
    for (int i =0; i < n_size; ++i) { 
      ppng->data[i] = float(c->p_r[i]);
    }
  
    float tmp_buf[2]; dumpRange(tmp_buf,ppng->data,n_size,"check data buf");
    ppng->initializeWrite(prange[0],prange[1]);
    ppng->write(filename);

    // p imag ... 

    sprintf(filename,"%s_p.imag.png",mode_name.c_str());
    for (int i =0; i < n_size; ++i) 
      ppng->data[i] = float(c->p_i[i]);
    
    ppng->initializeWrite(prange[0],prange[1]);
    ppng->write(filename);

    // p mag .. 

    sprintf(filename,"%s_p.mag.png",mode_name.c_str());
    for (int i =0; i < n_size; ++i) 
      ppng->data[i] = float ( sqrt(c->p_r[i]*c->p_r[i] + c->p_i[i]*c->p_i[i])); 
    

    double pmag_range[2] = {0.0, sqrt(prange[0]*prange[0] + prange[1]*prange[1])};
    ppng->initializeWrite(pmag_range[0],pmag_range[1]);
    ppng->write(filename);

    // z real... 

    sprintf(filename,"%s_z.real.png",mode_name.c_str());
    for (int i =0; i < n_size; ++i) 
      zpng->data[i] = float(c->z_r[i]);

    zpng->initializeWrite(zrange[0],zrange[1]);
    zpng->write(filename);

    // z imag ... 

    sprintf(filename,"%s_z.imag.png",mode_name.c_str());
    for (int i =0; i < n_size; ++i) 
      zpng->data[i] = float(c->z_i[i]);
    
    zpng->initializeWrite(zrange[0],zrange[1]);
    zpng->write(filename);

    // z mag .. 

    sprintf(filename,"%s_z.mag.png",mode_name.c_str());
    for (int i =0; i < n_size; ++i) 
      zpng->data[i] = float ( sqrt(c->z_r[i]*c->z_r[i] + c->z_i[i]*c->z_i[i])); 
    
    double zmag_range[2] = {0.0, sqrt(zrange[0]*zrange[0] + zrange[1]*zrange[1])};
    zpng->initializeWrite(zmag_range[0],zmag_range[1]);
    zpng->write(filename);


    delete ppng;
    delete zpng;

  }
};


int main(int argc, char* argv[]) { 

  initMpiEnv(argc,argv);

  //ImageSnapshot2* data = new ImageSnapshot2();

  // stack p, Z

  /*
  vector<string> data_names; data_names.resize(2);
  vector<string> img_prefix; img_prefix.resize(2);
  vector<double> d_ref;      d_ref.resize(2);

  data_names[0] = "p"; img_prefix[0] = "./images/p"; d_ref[0] = 8.0e+5;
  data_names[1] = "z"; img_prefix[1] = "./images/z"; d_ref[1] = 0.1;
  */

  // stacked p, T, Z data ... 

  
  vector<string> data_names; data_names.resize(3);
  vector<string> img_prefix; img_prefix.resize(3);
  vector<double> d_ref;      d_ref.resize(3);

  data_names[0] = "p"; img_prefix[0] = "./images_50b/p"; d_ref[0] = 8.0e+5;
  data_names[1] = "z"; img_prefix[1] = "./images_50b/z"; d_ref[1] = 0.1;
  //data_names[2] = "T"; img_prefix[2] = "./images/t"; d_ref[2] = 2000.0;
  data_names[2] = "u"; img_prefix[2] = "./images_50b/u2"; d_ref[2] = 250.0;

  ImageSnapshotM* data = new ImageSnapshotM(data_names,img_prefix,d_ref);
  
  // XXX hack -- this should be parsed ... 

  //int istart = 3283000;
  //int istart = 3383000;
  //int istep  = 1000;
  //int iend   = 3482000;
  //double dt_samp = 1.4e-04;

  int istart = 6120000;
  //int istart = 3500000;
  int istep  = 50;
  int iend   = 6189950;
  //int iend   = 3750000;
  double dt_samp = 50.0*1.5e-07;

  CompositeDmd cdmd(istart,istep,iend,data,dt_samp);
  cdmd.run();

  delete data;
  return 0;
}

