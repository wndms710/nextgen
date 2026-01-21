#include "ChemtableGPU.hpp"

__device__
inline int findIndexMap(const double& x, indexMap * idxMap) { 
  const double xi = fmin(idxMap->xmax-idxMap->xsmall,
      fmax(idxMap->xmin+idxMap->xsmall,x));
  const int j = floor((xi-idxMap->xmin)*idxMap->diloc);
  return floor(idxMap->iloc[j]+(idxMap->iloc[j+1]-idxMap->iloc[j])*
      (xi-idxMap->xloc[j])*idxMap->diloc);
} 

__device__ 
inline void cubicInterpWtsD(double* w, const double* inv_denom, const double * x, const double xi) { 
  const double d0=xi-x[0];
  const double d1=xi-x[1];
  const double d2=xi-x[2];
  const double d3=xi-x[3];
  w[0]=d1*d2*d3*inv_denom[0];
  w[1]=d2*d3*d0*inv_denom[1];
  w[2]=d3*d0*d1*inv_denom[2];
  w[3]=d0*d1*d2*inv_denom[3];
}

__global__
void lookupSpecialKernel(double * rout, const double * y1, const double * y2a, 
    const double * y2b, const double * data, const double* x1d, const double* x2d, const double* x2lowerd, 
    const double* x2upperd, const double (*invDenom1d)[4], const double (*invDenom2d)[4], indexMap* idxMap1d, 
    indexMap* idxMap2d, const int n1, const int n2, const int n) { 
  
  const int i = blockIdx.x * blockDim.x + threadIdx.x; 
  if ( i < n ) { 
    double w1[4], w2[4], w2_2[4]; 
    const double Y1 = min(x1d[n1-1],max(x1d[0],y1[i])); 
    
    const int i1 = min(n1-4, max(0, findIndexMap(Y1,idxMap1d)-1));
    cubicInterpWtsD(w1,invDenom1d[i1],&x1d[i1],Y1);
    
    double x2l = 0.0, x2u = 0.0;
    for (int k =0; k < 4; ++k) 
      x2l += w1[k]*x2lowerd[i1+k];
    
    for (int k=0; k < 4; ++k) 
      x2u += w1[k]*x2upperd[i1+k];
    
    const double den = x2u - x2l + 1.0e-16; 
    const double Y2   = min(1.0,max(0.0,(y2a[i]-x2l)/den));
    const double Y2_2 = min(1.0,max(0.0,(y2b[i]-x2l)/den));
    
    const int i2   = min(n2-4,max(0,findIndexMap(Y2, idxMap2d)-1)); 
    const int i2_2 = min(n2-4,max(0,findIndexMap(Y2_2,idxMap2d)-1)); 
    
    cubicInterpWtsD(w2,invDenom2d[i2],&x2d[i2],Y2);
    cubicInterpWtsD(w2_2,invDenom2d[i2_2],&x2d[i2_2],Y2_2); 
    
    double r_1[4] = {0.0,0.0,0.0,0.0}; 
    double r_2[4] = {0.0,0.0,0.0,0.0}; 
    
    for (int k =0; k < 4 ; ++k) { 
      const int ii       = (i1+k)*n2; 
      const int offset_1 = ii+i2;   
      const int offset_2 = ii+i2_2; 
      
      for (int l=0; l < 4 ; ++l)  
        r_1[l] += w1[k]* data[offset_1+l] ; 
      
      for (int l=0; l < 4 ; ++l) 
        r_2[l] += w1[k]* data[offset_2+l]; 
    }//k
    
    for (int k=0; k < 4 ; ++k) 
      r_1[k] *= w2[k]; 
    
    for (int k =0; k < 4 ; ++k) 
      r_2[k] *= w2_2[k]; 
    
    const double sum1 = (r_1[0] + r_1[1]) + (r_1[2] + r_1[3]); 
    const double sum2 = (r_2[0] + r_2[1]) + (r_2[2] + r_2[3]); 
    rout[i] = fabs(sum1-sum2) ; 
  }//i
}//lookupSpecialKernel.. 

__global__
void lookupReducedVectorKernel(double** routd, const double* y1, 
    double** datav, const double* x1d, const double (*invDenom1d)[4],  
    indexMap* idxMap1d, const int ndata, const int n1, const int n2, const int n) { 

  const int i = blockIdx.x * blockDim.x + threadIdx.x; 
  if ( i < n ) { 
    double w1[4]; 
    const double Y1 = min(x1d[n1-1],max(x1d[0],y1[i])); 
    
    const int i1 = min(n1-4, max(0, findIndexMap(Y1,idxMap1d)-1));
    cubicInterpWtsD(w1,invDenom1d[i1],&x1d[i1],Y1);

    for (int j =0; j < ndata ; ++j) { 
      const double * data = datav[j];
      double val          = 0.0;
      // the data vector is 2d but it doesn't depend on the 
      // second coord, so we'll use idx 0.  going forward, 
      // these data values should be precomputed in 1d...
      for (int k =0 ; k < 4 ; ++k) { 
        const int offset = (i1+k)*n2; 
        val += w1[k]*data[offset]; 
      }
      routd[j][i] = val;
    }
  }
}

__global__
void lookupCubicPtrVectorKernel(double** routd, const double* y1, const double *y2, 
    double** datav, const double* x1d, const double* x2d, const double* x2lowerd, 
    const double* x2upperd, const double (*invDenom1d)[4], const double (*invDenom2d)[4], 
    indexMap* idxMap1d, indexMap* idxMap2d, const int ndata, const int n1, const int n2, const int n) { 

  const int i =blockDim.x * blockIdx.x + threadIdx.x;  
  if ( i < n ) { 
    double w1[4], w2[4];
    const double Y1 = min(x1d[n1-1],max(x1d[0], y1[i])); 
    const int i1 = min(n1-4, max(0, findIndexMap(Y1,idxMap1d)-1));
    cubicInterpWtsD(w1,invDenom1d[i1],&x1d[i1],Y1);
    
    double x2l = 0.0, x2u = 0.0;
    for (int k =0; k < 4; ++k) 
      x2l += w1[k]*x2lowerd[i1+k];
    
    for (int k=0; k < 4; ++k) 
      x2u += w1[k]*x2upperd[i1+k];
    
    const double den = x2u - x2l + 1.0e-16; 
    const double Y2   = min(1.0,max(0.0,(y2[i]-x2l)/den));
    
    const int i2 = min(n2-4,max(0,findIndexMap(Y2, idxMap2d)-1));
    cubicInterpWtsD(w2,invDenom2d[i2],&x2d[i2],Y2); 
  
    for (int j =0; j < ndata; ++j) { 
      const double * data = datav[j];
      double r[4]   = {0.0,0.0,0.0,0.0};
      
      for (int k =0; k<4; ++k) { 
        double tmp[4]; 
        const int offset = (i1+k)*n2+i2; 
        for (int l=0; l<4; ++l) 
          tmp[l] = w1[k]*data[offset+l];
        
        for (int l=0; l<4; ++l) 
          r[l] += tmp[l];
      }//k
      
      for (int k=0; k<4; ++k) 
        r[k] *= w2[k];
      
      routd[j][i] = (r[0]+r[1]) +(r[2]+r[3]); 
    }
  }
} 

__global__
void computeFaceReductionKernel(double * Z_fa_d,double * C0_fa_d,double * C1_fa_d,
                                const double * Z_cv_d, const double * C_cv_d, const int (*cvofa_d)[2],const int n) { 
  const int i =blockDim.x * blockIdx.x + threadIdx.x;  
  if ( i < n ) {
    const int i0 = cvofa_d[i][0];
    const int i1 = cvofa_d[i][1];
    Z_fa_d[i]    = 0.5*(Z_cv_d[i0] + Z_cv_d[i1]); 
    C0_fa_d[i]   = C_cv_d[i0];
    C1_fa_d[i]   = C_cv_d[i1];
  }
} 

void computeFaceReduction(double* Z_fa_d, double* C0_fa_d, double * C1_fa_d, 
                          const double* Z_cv_d, const double* C_cv_d, const int (*cvofa_d)[2], const int n) { 
  const int block_size = 256;
  const int grid_size  = (n+block_size-1)/block_size;
  computeFaceReductionKernel<<<grid_size,block_size>>>(Z_fa_d,C0_fa_d,C1_fa_d,
                                                       Z_cv_d,C_cv_d,cvofa_d,n); 
}

void CartesianChemtable2dGpu::lookupSpecial(double* rout, const string& name, 
                                            const double * y1d, const double * y2ad, 
                                            const double * y2bd, const int n) { 
  const double * datad = deviceVars[name]; //device ptr
  double * tmp = NULL; 
  cudaMalloc((void**)&tmp,n*sizeof(double)); 
  cudaBuffers[name] = tmp;

  const int block_size  = 256;
  const int grid_size   = (n+block_size-1)/block_size; 
  lookupSpecialKernel<<<grid_size,block_size>>>(tmp, y1d, y2ad, 
      y2bd, datad, x1d, x2d, x2lowerd, x2upperd, invDenom1d, invDenom2d, idxMap1d, idxMap2d, n1, n2, n);
  
  // in order to support async operations on the host, we will not copy 
  // back the result of this operation until later...
}

void CartesianChemtable2dGpu::lookupCubicPtrVector(vector<double*>& routPtrVec, const vector<string>& nameVec,
                                                   const double *y1, const double *y2,const int n) {
  const int ndata = routPtrVec.size();
  
  // need to allocate the cuda buffers...
  double** ptrs = new double*[ndata]; 
  for (int j =0; j < ndata; ++j) { 
    double * tmp = NULL;
    cudaMalloc((void**)&tmp,n*sizeof(double));
    cudaBuffers[nameVec[j]] = tmp;
    ptrs[j]                 = tmp;
  }

  double ** routd = NULL; 
  cudaMalloc((void**)&routd, ndata*sizeof(double*));
  cudaMemcpy(routd,ptrs,ndata*sizeof(double*),cudaMemcpyHostToDevice); 
  
  for (int j =0; j < ndata; ++j) 
    ptrs[j] = deviceVars[nameVec[j]];
  
  double ** datad = NULL;
  cudaMalloc((void**)&datad, ndata*sizeof(double*));
  cudaMemcpy(datad,ptrs,ndata*sizeof(double*),cudaMemcpyHostToDevice);
  
  delete[] ptrs;

  int block_size = 256;
  int grid_size  = (n+block_size-1)/block_size;

  lookupCubicPtrVectorKernel<<<grid_size,block_size>>>(routd, y1, y2, datad, x1d, 
      x2d, x2lowerd, x2upperd, invDenom1d, invDenom2d, idxMap1d, idxMap2d, ndata, n1, n2, n);
 
  // we have this buffered so you don't have to copy the data back yet, 
  // but for simplicity, we'll copy it back now... 
  for (int j=0; j < ndata; ++j) { 
    double * d_data = cudaBuffers[nameVec[j]];
    cudaMemcpy(routPtrVec[j],d_data,n*sizeof(double),cudaMemcpyDeviceToHost);
    cudaFree(d_data);
    cudaBuffers.erase(nameVec[j]);
  }
  
  cudaFree(routd);
  cudaFree(datad);
}

void CartesianChemtable2dGpu::lookupDataVecReduced(vector<double*>& routPtrVec, const vector<string>& nameVec, 
                                                   const double * y1, const int n) { 
  
  const int ndata = routPtrVec.size(); 
  double ** ptrs = new double*[ndata];
  for (int j =0; j < ndata; ++j) { 
    double * tmp = NULL;
    cudaMalloc((void**)&tmp,n*sizeof(double));
    cudaBuffers[nameVec[j]] = tmp;
    ptrs[j]                 = tmp;
  }

  double ** routd = NULL; 
  cudaMalloc((void**)&routd, ndata*sizeof(double*));
  cudaMemcpy(routd,ptrs,ndata*sizeof(double*),cudaMemcpyHostToDevice); 
  
  for (int j =0; j < ndata; ++j)  {
      ptrs[j] = deviceVars[nameVec[j]];
  }

  double ** datad = NULL;
  cudaMalloc((void**)&datad, ndata*sizeof(double*));
  cudaMemcpy(datad,ptrs,ndata*sizeof(double*),cudaMemcpyHostToDevice);
  
  delete[] ptrs;

  int block_size = 256;
  int grid_size  = (n+block_size-1)/block_size;

  lookupReducedVectorKernel<<<grid_size,block_size>>>(routd, y1, datad, x1d, invDenom1d, 
                                                      idxMap1d,ndata,n1,n2,n); 

  // we have this buffered so you don't have to copy the data back yet, 
  // but for simplicity, we'll copy it back now... 
  for (int j=0; j < ndata; ++j) { 
    double * d_data = cudaBuffers[nameVec[j]];
    cudaMemcpy(routPtrVec[j],d_data,n*sizeof(double),cudaMemcpyDeviceToHost);
    cudaFree(d_data);
    cudaBuffers.erase(nameVec[j]);
  }
  
  cudaFree(routd);
  cudaFree(datad);
}


void initChemtableGpu(AbstractChemtable2D * &chemtable, const string& tablename) {
  
  string tabletype = getChemtableType(tablename);
  COUT1(" > initializing 2D table: " << tablename << " with the type: " << tabletype);
  
  if ( (tabletype == "VIDA_PREMIXED_FPV_CART2D") ||
       (tabletype == "CHARLES_PREMIXED_FPV_CART2D")||(tabletype == "PREMIXED"))
    chemtable = new CartesianChemtable2dGpu(tablename);
  else
    CERR("incompatible 2D table type " <<  tabletype);
} 

void deleteChemtableGpu(AbstractChemtable2D * &chemtable) {
  if ( chemtable != NULL ) delete chemtable;
}

