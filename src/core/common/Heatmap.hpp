#ifndef _HEATMAP_HPP_
#define _HEATMAP_HPP_

enum HeatmapAction {
  HEATMAP_SUM = 1,
  HEATMAP_MAX = 2,
};

class Heatmap {
public:

  int write_rank; // rank to build png on (default to 0, updated by StaticSolver)

private:

  int action;
  
  int ni,nj; // image pixel size: width,height
  int nI,nJ; // block size (image divided into 8x8 blocks)

  bool b_x_range;
  double x_range[2];
  
  bool b_y_range;
  double y_range[2];

  class HeatmapBlockData {
  public:
    double wgt[64]; // 8x8 set of weights
    bool flag[64];
    HeatmapBlockData() {
      for (int ij = 0; ij < 64; ++ij) wgt[ij] = 0.0;
      for (int ij = 0; ij < 64; ++ij) flag[ij] = false; // i.e. NOT active
    }
  };

  // a local version of the image and depth
  HeatmapBlockData ** hbd;

  double * buf_double;
  uint2 * buf_uint2;
  int buf_size;
  
  //ImageMetadata * imd; ???

public:

  //CtiCanvas() : imd(NULL){};
  
  Heatmap(const int ni_,const int nj_) {
    
    COUT1("Heatmap()"); 
 
    // default action = sum

    action = HEATMAP_SUM;
    
    // set the size...
    
    ni = ni_;
    nj = nj_;
    
    nI = ((ni-1)>>3)+1; assert((nI >= 1)&&(nI < 2000));
    nJ = ((nj-1)>>3)+1; assert((nJ >= 1)&&(nJ < 2000));

    // clear range...
    
    b_x_range = false;
    b_y_range = false;

    // allocate hbd...
    
    hbd = new HeatmapBlockData*[nI*nJ];
    for (int ii = 0; ii < nI*nJ; ++ii) hbd[ii] = NULL;
    
    // and NULLIFY reduced buffers...
    
    write_rank = -1;
    buf_double = NULL;
    buf_uint2 = NULL;
    buf_size = 0;
    
  }
  
  ~Heatmap() {

    COUT1("~Heatmap()"); 
    
    if (hbd) {
      for (int ii = 0; ii < nI*nJ; ++ii) {
        if (hbd[ii] != NULL) delete hbd[ii];
      }
      delete[] hbd;
    }

    DELETE(buf_double);
    DELETE(buf_uint2);
    
  }

  void setAction(const int action_) {
    
    action = action_;
    
  }
  
  void addData(double * x,double * y,double * wgt,const int n) {
    
    // add x,y data to the heat map...
    
    if ((!b_x_range)&&(!b_y_range)) {
      double my_buf[4] = {1.0E+20,1.0E+20,1.0E+20,1.0E+20};
      for (int i = 0; i < n; ++i) {
	my_buf[0] = min(my_buf[0],x[i]);
	my_buf[1] = min(my_buf[1],-x[i]);
	my_buf[2] = min(my_buf[2],y[i]);
	my_buf[3] = min(my_buf[3],-y[i]);
      }
      double buf[4];
      MPI_Allreduce(my_buf,buf,4,MPI_DOUBLE,MPI_MIN,mpi_comm);
      x_range[0] = buf[0];
      x_range[1] = -buf[1] - 1.0E-6*(buf[1]+buf[0]); // add a little tolerance to the top
      b_x_range = true;
      y_range[0] = buf[2];
      y_range[1] = -buf[3] - 1.0E-6*(buf[3]+buf[2]); // add a little tolerance to the top
      b_y_range = true;
      if (mpi_rank == 0) cout << " > heatmap range set based on first addData: x: " << x_range[0] << " " << x_range[1] << ", y: " << y_range[0] << " " << y_range[1] << endl;
    }
    else if (!b_x_range) {
      assert(0); // do this later
    }
    else if (!b_y_range) {
      assert(0); // do this later
    }
    
    // now add the weights to the data...
    
    switch (action) {
    case HEATMAP_SUM: 
      for (int i = 0; i < n; ++i) {
	const double hm_i = (x[i] - x_range[0])/(x_range[1] - x_range[0])*double(ni);
	const double hm_j = (y[i] - y_range[0])/(y_range[1] - y_range[0])*double(nj);
	assert((hm_i >= 0)&&(hm_i < ni));
	assert((hm_j >= 0)&&(hm_j < nj));
	if ((hm_i >= 0)&&(hm_i < ni)&&(hm_j >= 0)&&(hm_j < nj)) {
	  addWgt((int)hm_i,(int)hm_j,wgt[i]);
	}
      }
      break;
    case HEATMAP_MAX:
      for (int i = 0; i < n; ++i) {
	const double hm_i = (x[i] - x_range[0])/(x_range[1] - x_range[0])*double(ni);
	const double hm_j = (y[i] - y_range[0])/(y_range[1] - y_range[0])*double(nj);
	assert((hm_i >= 0)&&(hm_i < ni));
	assert((hm_j >= 0)&&(hm_j < nj));
	if ((hm_i >= 0)&&(hm_i < ni)&&(hm_j >= 0)&&(hm_j < nj)) {
	  maxWgt((int)hm_i,(int)hm_j,wgt[i]);
	}
      }
      break;
    default:
      assert(0);
    }

  }

  void writeTecplot(const string& filename) {

    COUT1("Heatmap::writeTecplot: " << filename);
    
    setWriteRank(0);
    reduceData();
    
    if (mpi_rank == write_rank) {
      
      double *wgt = new double[ni*nj];
      for (int ij = 0; ij < ni*nj; ++ij)
	wgt[ij] = 0.0;

      switch (action) {
      case HEATMAP_SUM:
	for (int ibuf = 0; ibuf < buf_size; ++ibuf) {
	  const int i = buf_uint2[ibuf*2  ]; assert((i >= 0)&&(i < ni));
	  const int j = buf_uint2[ibuf*2+1]; assert((j >= 0)&&(j < nj));
	  wgt[j*ni+i] += buf_double[ibuf];
	}
	break;
      case HEATMAP_MAX:
	for (int ibuf = 0; ibuf < buf_size; ++ibuf) {
	  const int i = buf_uint2[ibuf*2  ]; assert((i >= 0)&&(i < ni));
	  const int j = buf_uint2[ibuf*2+1]; assert((j >= 0)&&(j < nj));
	  wgt[j*ni+i] = max(wgt[j*ni+i],buf_double[ibuf]);
	}
	break;
      default:
	assert(0);
      }
      
      FILE * fp = fopen(filename.c_str(),"w");
      fprintf(fp,"TITLE = \"heatmap\"\n");
      fprintf(fp,"VARIABLES = \"X\"\n\"Y\"\n\"Z\"\n");
      fprintf(fp,"ZONE T=\"heatmap\"\n");
      fprintf(fp,"I=%d, J=%d, K=1, ZONETYPE=Ordered DATAPACKING=POINT\n",ni,nj);
      for (int j = 0; j < nj; ++j) {
	const double y = y_range[0] + (double(j)+0.5)/double(nj)*(y_range[1]-y_range[0]);
	for (int i = 0; i < ni; ++i) {
	  const double x = x_range[0] + (double(i)+0.5)/double(ni)*(x_range[1]-x_range[0]);
	  fprintf(fp,"%18.16e %18.16e %18.16e\n",x,y,wgt[j*ni+i]);
	}
      }
      fclose(fp);
      
      delete[] wgt;

    }

  }
  
private:
  
  void addWgt(const int i,const int j,const double wgt) {
    assert((i >= 0)&&(i < ni));
    assert((j >= 0)&&(j < nj));
    const int IJ = (j>>3)*nI + (i>>3);
    assert((IJ >= 0)&&(IJ < nI*nJ));
    if (hbd[IJ] == NULL) hbd[IJ] = new HeatmapBlockData();
    const int ij = (j&7)*8+(i&7);
    assert((ij >= 0)&&(ij < 64));
    hbd[IJ]->flag[ij] = true;
    hbd[IJ]->wgt[ij] += wgt;
  }
  
  void maxWgt(const int i,const int j,const double wgt) {
    assert((i >= 0)&&(i < ni));
    assert((j >= 0)&&(j < nj));
    const int IJ = (j>>3)*nI + (i>>3);
    assert((IJ >= 0)&&(IJ < nI*nJ));
    if (hbd[IJ] == NULL) hbd[IJ] = new HeatmapBlockData();
    const int ij = (j&7)*8+(i&7);
    assert((ij >= 0)&&(ij < 64));
    if (!hbd[IJ]->flag[ij]) {
      hbd[IJ]->flag[ij] = true;
      hbd[IJ]->wgt[ij] = wgt;
    }
    else if (wgt > hbd[IJ]->wgt[ij]) {
      hbd[IJ]->wgt[ij] = wgt;
    }
  }
  
  void setWriteRank(const int write_rank_) {
    write_rank = write_rank_;
    assert((write_rank >= 0)&&(write_rank < mpi_size));
    assert (buf_double==NULL);
    assert (buf_uint2==NULL);
    assert (buf_size==0);
  }
  
  void reduceData() {
    
    assert(hbd);
    assert((write_rank >= 0)&&(write_rank < mpi_size));
    assert (buf_double==NULL);
    assert (buf_uint2==NULL);
    assert (buf_size==0);
      
    // count first...

    int my_buf_size = 0;
    for (int J = 0; J < nJ; ++J) {
      for (int I = 0; I < nI; ++I) {
	const int IJ = J*nI+I;
	if (hbd[IJ]) {
	  for (int j = 0; j < 8; ++j) {
	    for (int i = 0; i < 8; ++i) {
	      const int ij = j*8+i;
	      if (hbd[IJ]->flag[ij]) {
		++my_buf_size;
	      }
	    }
	  }
	}
      }
    }
      
    double * my_buf_double = new double[my_buf_size]; // wgt
    uint2 * my_buf_uint2 = new uint2[my_buf_size*2]; // i,j
      
    my_buf_size = 0;
    for (int J = 0; J < nJ; ++J) {
      for (int I = 0; I < nI; ++I) {
	const int IJ = J*nI+I;
	if (hbd[IJ]) {
	  for (int j = 0; j < 8; ++j) {
	    for (int i = 0; i < 8; ++i) {
	      const int ij = j*8+i;
	      if (hbd[IJ]->flag[ij]) {
		// float data...
		my_buf_double[my_buf_size]    = hbd[IJ]->wgt[ij];
		my_buf_uint2[my_buf_size*2  ] = (I<<3)+i;
		my_buf_uint2[my_buf_size*2+1] = (J<<3)+j;
		// increment...
		++my_buf_size;
	      }
	    }
	  }
	}
      }
    }

    // now gather at the write process (not necessarily rank==0, but set this for now)...
      
    if (mpi_size == 1) {
	
      // in serial mode, just take copies of the buffers...
      assert(write_rank == 0);
      buf_size = my_buf_size;
      buf_double = my_buf_double;
      buf_uint2 = my_buf_uint2;

    }
    else {
	
      int * count = NULL;
      if (mpi_rank == write_rank) count = new int[mpi_size];
      MPI_Gather(&my_buf_size,1,MPI_INT,count,1,MPI_INT,write_rank,mpi_comm);

      int * disp = NULL;
      if (mpi_rank == write_rank) {
	disp = new int[mpi_size];
	disp[0] = 0;
	for (int rank = 1; rank < mpi_size; ++rank)
	  disp[rank] = disp[rank-1] + count[rank-1];
	buf_size = disp[mpi_size-1] + count[mpi_size-1];
	buf_double = new double[buf_size];
	buf_uint2 = new uint2[buf_size*2];
      }
	
      MPI_Gatherv(my_buf_double,my_buf_size,MPI_DOUBLE,buf_double,count,disp,MPI_DOUBLE,write_rank,mpi_comm);
      delete[] my_buf_double;
	
      if (mpi_rank == write_rank) {
	FOR_RANK {
	  count[rank] *= 2;
	  disp[rank] *= 2;
	}
      }
	
      MPI_Gatherv(my_buf_uint2,my_buf_size*2,MPI_UINT2,buf_uint2,count,disp,MPI_UINT2,write_rank,mpi_comm);
      delete[] my_buf_uint2;
      
      if (mpi_rank == write_rank) {
	delete[] count;
	delete[] disp;
      }
	
    }

  }
  
};

#endif
