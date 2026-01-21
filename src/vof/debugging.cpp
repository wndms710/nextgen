
/* 
 * svdcomp - SVD decomposition routine. 
 * Takes an mxn matrix a and decomposes it into udv, where u,v are
 * left and right orthogonal transformation matrices, and d is a 
 * diagonal matrix of singular values.
 *
 * This routine is adapted from svdecomp.c in XLISP-STAT 2.1 which is 
 * adapted by Luke Tierney and David Betz.
 *
 * the now dead xlisp-stat package seems to have been distributed
 * under some sort of BSD license.
 *
 * Input to dsvd is as follows:
 *   a = mxn matrix to be decomposed, gets overwritten with u
 *   m = row dimension of a
 *   n = column dimension of a
 *   w = returns the vector of singular values of a
 *   v = returns the right orthogonal transformation matrix
*/

#pragma once

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define COPY_SIGN(a, b) copysign(a, b)

static inline double PYTHAG(double a, double b)
{
  double at = fabs(a), bt = fabs(b), ct, result;

  if (at > bt)       { ct = bt / at; result = at * sqrt(1.0 + ct * ct); }
  else if (bt > 0.0) { ct = at / bt; result = bt * sqrt(1.0 + ct * ct); }
  else result = 0.0;
  return result;
}


// decompose (m >= n)
//      n             n               n
//   |      |      |     |   n     |     |
// m |  a   |  = m |  u  | diag(w) | v^t | n
//   |      |      |     |         |     |
//
// where the data layout of a (in) and u (out) is strided by str for every row
static inline int dsvd(
    double *a,    // input matrix a[j*str + i] is j-th row and i-th column. will be overwritten by u
    int m,        // number of rows of a and u
    int n,        // number of cols of a and u
    int str,      // row stride of a and u
    double *w,    // output singular values w[n]
    double *v)    // output v matrix v[n*n]
{
  int flag, i, its, j, jj, k, l, nm;
  double c, f, h, s, x, y, z;
  double anorm = 0.0, g = 0.0, scale = 0.0;

  if (m < n) 
  {
    fprintf(stderr, "[svd] #rows must be >= #cols \n");
    return(0);
  }

  double *rv1 = new double[n];

  /* Householder reduction to bidiagonal form */
  for (i = 0; i < n; i++) 
  {
    /* left-hand reduction */
    l = i + 1;
    rv1[i] = scale * g;
    g = s = scale = 0.0;
    if (i < m) 
    {
      for (k = i; k < m; k++) 
        scale += fabs(a[k*str+i]);
      if (scale) 
      {
        for (k = i; k < m; k++) 
        {
          a[k*str+i] = a[k*str+i]/scale;
          s += a[k*str+i] * a[k*str+i];
        }
        f = a[i*str+i];
        g = -COPY_SIGN(sqrt(s), f);
        h = f * g - s;
        a[i*str+i] = f - g;
        if (i != n - 1) 
        {
          for (j = l; j < n; j++) 
          {
            for (s = 0.0, k = i; k < m; k++) 
              s += a[k*str+i] * a[k*str+j];
            f = s / h;
            for (k = i; k < m; k++) 
              a[k*str+j] += f * a[k*str+i];
          }
        }
        for (k = i; k < m; k++) 
          a[k*str+i] = a[k*str+i]*scale;
      }
    }
    w[i] = scale * g;

    /* right-hand reduction */
    g = s = scale = 0.0;
    if (i < m && i != n - 1) 
    {
      for (k = l; k < n; k++) 
        scale += fabs(a[i*str+k]);
      if (scale) 
      {
        for (k = l; k < n; k++) 
        {
          a[i*str+k] = a[i*str+k]/scale;
          s += a[i*str+k] * a[i*str+k];
        }
        f = a[i*str+l];
        g = -COPY_SIGN(sqrt(s), f);
        h = f * g - s;
        a[i*str+l] = f - g;
        for (k = l; k < n; k++) 
          rv1[k] = a[i*str+k] / h;
        if (i != m - 1) 
        {
          for (j = l; j < m; j++) 
          {
            for (s = 0.0, k = l; k < n; k++) 
              s += a[j*str+k] * a[i*str+k];
            for (k = l; k < n; k++) 
              a[j*str+k] += s * rv1[k];
          }
        }
        for (k = l; k < n; k++) 
          a[i*str+k] = a[i*str+k]*scale;
      }
    }
    anorm = MAX(anorm, (fabs(w[i]) + fabs(rv1[i])));
  }

  /* accumulate the right-hand transformation */
  for (i = n - 1; i >= 0; i--) 
  {
    if (i < n - 1) 
    {
      if (g) 
      {
        for (j = l; j < n; j++)
          v[j*n+i] = a[i*str+j] / a[i*str+l] / g;
        /* double division to avoid underflow */
        for (j = l; j < n; j++) 
        {
          for (s = 0.0, k = l; k < n; k++) 
            s += a[i*str+k] * v[k*n+j];
          for (k = l; k < n; k++) 
            v[k*n+j] += s * v[k*n+i];
        }
      }
      for (j = l; j < n; j++) 
        v[i*n+j] = v[j*n+i] = 0.0;
    }
    v[i*n+i] = 1.0;
    g = rv1[i];
    l = i;
  }

  /* accumulate the left-hand transformation */
  for (i = n - 1; i >= 0; i--) 
  {
    l = i + 1;
    g = w[i];
    if (i < n - 1) 
      for (j = l; j < n; j++) 
        a[i*str+j] = 0.0;
    if (g) 
    {
      g = 1.0 / g;
      if (i != n - 1) 
      {
        for (j = l; j < n; j++) 
        {
          for (s = 0.0, k = l; k < m; k++) 
            s += a[k*str+i] * a[k*str+j];
          f = (s / a[i*str+i]) * g;
          for (k = i; k < m; k++) 
            a[k*str+j] += f * a[k*str+i];
        }
      }
      for (j = i; j < m; j++) 
        a[j*str+i] = a[j*str+i]*g;
    }
    else 
    {
      for (j = i; j < m; j++) 
        a[j*str+i] = 0.0;
    }
    ++a[i*str+i];
  }

  /* diagonalize the bidiagonal form */
  for (k = n - 1; k >= 0; k--) 
  {                             /* loop over singular values */
    for (its = 0; its < 30; its++) 
    {                         /* loop over allowed iterations */
      flag = 1;
      for (l = k; l >= 0; l--) 
      {                     /* test for splitting */
        nm = MAX(0, l - 1);
        if (fabs(rv1[l]) + anorm == anorm) 
        {
          flag = 0;
          break;
        }
        if (l == 0 || fabs(w[nm]) + anorm == anorm)
          break;
      }
      if (flag) 
      {
        c = 0.0;
        s = 1.0;
        for (i = l; i <= k; i++) 
        {
          f = s * rv1[i];
          if (fabs(f) + anorm != anorm) 
          {
            g = w[i];
            h = PYTHAG(f, g);
            w[i] = h; 
            h = 1.0 / h;
            c = g * h;
            s = (- f * h);
            for (j = 0; j < m; j++) 
            {
              y = a[j*str+nm];
              z = a[j*str+i];
              a[j*str+nm] = y * c + z * s;
              a[j*str+i]  = z * c - y * s;
            }
          }
        }
      }
      z = w[k];
      if (l == k) 
      {                  /* convergence */
        if (z < 0.0) 
        {              /* make singular value nonnegative */
          w[k] = -z;
          for (j = 0; j < n; j++) 
            v[j*n+k] = (-v[j*n+k]);
        }
        break;
      }
      if (its >= 30) {
        fprintf(stderr, "[svd] no convergence after 30,000! iterations\n");
        delete[] rv1;
        return 0;
      }

      /* shift from bottom 2 x 2 minor */
      x = w[l];
      nm = k - 1;
      y = w[nm];
      g = rv1[nm];
      h = rv1[k];
      f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
      g = PYTHAG(f, 1.0);
      f = ((x - z) * (x + z) + h * ((y / (f + COPY_SIGN(g, f))) - h)) / x;

      /* next QR transformation */
      c = s = 1.0;
      for (j = l; j <= nm; j++) 
      {
        i = j + 1;
        g = rv1[i];
        y = w[i];
        h = s * g;
        g = c * g;
        z = PYTHAG(f, h);
        rv1[j] = z;
        c = f / z;
        s = h / z;
        f = x * c + g * s;
        g = g * c - x * s;
        h = y * s;
        y = y * c;
        for (jj = 0; jj < n; jj++) 
        {
          x = v[jj*n+j];
          z = v[jj*n+i];
          v[jj*n+j] = x * c + z * s;
          v[jj*n+i] = z * c - x * s;
        }
        z = PYTHAG(f, h);
        w[j] = z;
        if (z) 
        {
          z = 1.0 / z;
          c = f * z;
          s = h * z;
        }
        f = (c * g) + (s * y);
        x = (c * y) - (s * g);
        for (jj = 0; jj < m; jj++) 
        {
          y = a[jj*str+j];
          z = a[jj*str+i];
          a[jj*str+j] = y * c + z * s;
          a[jj*str+i] = z * c - y * s;
        }
      }
      rv1[l] = 0.0;
      rv1[k] = f;
      w[k] = x;
    }
  }
  delete[] rv1;
  return 1;
}
  void solveSvd(double x[6], const double A[6][6], const double b[6]) {

    double U[36], V[36], sigma[6];
    FOR_I6 FOR_J6 U[i*6+j] = A[i][j];
    dsvd(U,6,6,6,sigma,V);
    
    double sigma_inv[6];
    const double sigma_epsilon = 1.0E-16;
    FOR_I6 {
      assert( sigma[i] >= 0.0);
      if ( sigma[i] < sigma_epsilon) { 
        sigma_inv[i] = 0.0; // cant trust the gradient in this dir.. 
      } 
      else { 
        sigma_inv[i] = 1.0/sigma[i];
      }
    }

    FOR_I6 FOR_J6 V[i*6+j] *= sigma_inv[j];
    double tmp[6];
    FOR_I6 {
      tmp[i] = 0.0;
      FOR_J6 tmp[i] += U[j*6+i]*b[j];
    }
    FOR_I6 {
      x[i] = 0.0;
      FOR_J6 x[i] += V[i*6+j]*tmp[j];
    }

  }

  void solveGaussElimationWithPivoting(double x[6], double A[6][6], double b[6]) {

    // forward reduction step...
 
    int rowx = 0; // keep count of the row interchanges
    for (int k = 0; k < 5; ++k) {
      double amax = A[k][k];
      int m = k;
      for (int i = k+1; i < 6; ++i) { // Find the row with largest pivot 
        double xfac = fabs(A[i][k]);
        if (xfac > amax) { 
          amax = xfac; 
          m = i;
        }
      }
      if (m != k) { // Row interchanges
        rowx = rowx + 1;
        double temp1 = b[k];
        b[k] = b[m];
        b[m] = temp1;
        for (int j = k; j < 6; ++j) {
          double temp = A[k][j];
          A[k][j] = A[m][j];
          A[m][j] = temp;
        }
      }
      for (int i = k+1; i < 6; ++i) {
        double xfac = A[i][k]/A[k][k];
        for (int j = k+1; j < 6; ++j) {
          A[i][j] = A[i][j]-xfac*A[k][j];
        }
        b[i] = b[i]-xfac*b[k];
      }
    }

    // back substition step...

    for (int j = 0; j < 6; ++j) {
      int k = 5-j;
      x[k] = b[k];
      for (int i = k+1; i < 6; ++i) {
        x[k] = x[k]-A[k][i]*x[i];
      }
      x[k] = x[k]/A[k][k];
    }

  }

  double computeDistToPlicPoly(const double xp[3],const int icv) {

    double dist = d_max;
    list<double>::iterator it = plicPointList[icv].begin();
    double x0[3],x1[3]; 
    FOR_I3 {
      x1[i] = *it;
      ++it;
    }
    while (it != plicPointList[icv].end()) {
      FOR_I3 {
        x0[i] = x1[i];
        x1[i] = *it;
        ++it;
      }
      dist = min(dist,sqrt(getPointToTriDist2(xp,x0,x1,plic_xc[icv])));
    }
    it = plicPointList[icv].begin();
    FOR_I3 {
      x0[i] = x1[i];
      x1[i] = *it;
      ++it;
    }
    dist = min(dist,sqrt(getPointToTriDist2(xp,x0,x1,plic_xc[icv])));
    return dist;
  }

  void writePlicOFF(const int icv) {
    if (plicPointList[icv].size() == 0) return;
    
    const int ntri = plicPointList[icv].size()/3;
    char filename[128];
    sprintf(filename,"tris.%04d.%06d.off",mpi_rank,icv);
    FILE * fp = fopen(filename,"w");
    fprintf(fp,"OFF\n");
    fprintf(fp,"%i %i %i\n",ntri+1,ntri,0);
    fprintf(fp,"%18.15le %18.15le %18.15le\n",plic_xc[icv][0],plic_xc[icv][1],plic_xc[icv][2]);
    list<double>::iterator it = plicPointList[icv].begin();
    double x0[3]; 
    while (it != plicPointList[icv].end()) {
      FOR_I3 {
        x0[i] = *it;
        ++it;
      }
      fprintf(fp,"%18.15le %18.15le %18.15le\n",x0[0],x0[1],x0[2]);
    }
    it = plicPointList[icv].begin();
    int i0,i1 = 1;
    double x1[3]; 
    FOR_I3 {
      x1[i] = *it;
      ++it;
    }
    while (it != plicPointList[icv].end()) {
      FOR_I3 {
        x0[i] = x1[i];
        x1[i] = *it;
        ++it;
      }
      i0 = i1;
      ++i1;
      fprintf(fp,"%i %i %i %i \n",3,i0,i1,0); 
    }
    i0 = i1;
    i1 = 1;
    fprintf(fp,"%i %i %i %i \n",3,i0,i1,0); 
    fclose(fp);
  }

  double computePointListAreaAndCentroid(double xc[3],list<double> points) {
    double area = 0.0;
    FOR_I3 xc[i] = 0.0;
    if (points.size() > 0) {
      list<double>::iterator it = points.begin();
      double x0[3],x1[3]; 
      FOR_I3 {
        x1[i] = *it;
        ++it;
      }
      const double xp[3] = {x1[0],x1[1],x1[2]};
      while (it != points.end()) {
        FOR_I3 {
          x0[i] = x1[i];
          x1[i] = *it;
          ++it;
        }
        const double nA[3] = TRI_NORMAL_2(x0,x1,xp);
        const double mag = MAG(nA);
        area += mag;
        FOR_I3 xc[i] += mag*(x0[i]+x1[i]+xp[i]);
        
      }
      it = points.begin();
      FOR_I3 {
        x0[i] = x1[i];
        x1[i] = *it;
        ++it;
      }
      const double nA[3] = TRI_NORMAL_2(x0,x1,xp);
      const double mag = MAG(nA);
      area += mag;
      FOR_I3 xc[i] += mag*(x0[i]+x1[i]+xp[i]);

      FOR_I3 xc[i] /= 3.0*area;
      area *= 0.5;
    }
    return area;
  }

  void intersectCvWithPlic(list<double>& points, const int icv) {

    map<const int,pair<int,int> > faceEdgePairMap;
    vector<double> unordered_points(3*(fpnpocv_i[icv+1]-fpnpocv_i[icv]));
    //cout << "icv: " << icv << endl;
    for (int fnoc = fpnpocv_i[icv]; fnoc != fpnpocv_i[icv+1]; ++fnoc) {
      const int ifa_l = fpnpocv_v[fnoc][0]; 
      const int ifa_r = fpnpocv_v[fnoc][1]; 
      const int ino0 = fpnpocv_v[fnoc][2];
      const int ino1 = fpnpocv_v[fnoc][3];

      const double g0 = g[icv] - ((x_no[ino0][0]-x_cv[icv][0])*n[icv][0] +
                                  (x_no[ino0][1]-x_cv[icv][1])*n[icv][1] +
                                  (x_no[ino0][2]-x_cv[icv][2])*n[icv][2]);
      const double g1 = g[icv] - ((x_no[ino1][0]-x_cv[icv][0])*n[icv][0] +
                                  (x_no[ino1][1]-x_cv[icv][1])*n[icv][1] +
                                  (x_no[ino1][2]-x_cv[icv][2])*n[icv][2]);
      
      // TODO node/PLIC intersection, g = 0...
      if (g1*g0 < 0.0) {
        //cout << "g1g0: " << g1*g0 << " " << g1 << " " << g0 << " " << COUT_VEC(x_no[ino0]) << " " << COUT_VEC(x_no[ino1]) << endl;
        // add to map...
        map<const int,pair<int,int> >::iterator it = faceEdgePairMap.find(ifa_l);
        //cout << g0 << " " << g1 << " " << it->second.first << " " << it->second.second << " " << ifa_l << " " << ifa_r << " " << fnoc-fpnpocv_i[icv] << endl;
        if (it != faceEdgePairMap.end()) {
          assert(it->second.first >= 0);
          //if (it->second.second != -1) {
          //  for (int ii = 0; ii < 3*(fnoc-fpnpocv_i[icv]); ii += 3) {
          //    if (unordered_points[ii] < HUGE_VAL) 
          //      cout << unordered_points[ii+0] << " " << unordered_points[ii+1] << " " << unordered_points[ii+2] << endl;
          //  }
          //  cout.flush();
          //}
          assert(it->second.second == -1);
          it->second.second = fnoc-fpnpocv_i[icv]; // vector index...
        }
        else {
          faceEdgePairMap[ifa_l] = pair<int,int>(fnoc-fpnpocv_i[icv],-1);
        }
        it = faceEdgePairMap.find(ifa_r);
        if (it != faceEdgePairMap.end()) {
          assert(it->second.first >= 0);
          assert(it->second.second == -1);
          it->second.second = fnoc-fpnpocv_i[icv]; // vector index...
        }
        else {
          faceEdgePairMap[ifa_r] = pair<int,int>(fnoc-fpnpocv_i[icv],-1);
        }

        // plic intersects edge, so store intersection point...
        const double factor = g0/(g1-g0); 
        FOR_I3 {
          unordered_points[3*(fnoc-fpnpocv_i[icv])+i] = x_no[ino0][i]-factor*(x_no[ino1][i]-x_no[ino0][i]);
        }
      }
      else {
        FOR_I3 {
          unordered_points[3*(fnoc-fpnpocv_i[icv])+i] = HUGE_VAL;
        }
      }
    }

    // check that each intersected face has two edges...
    for (map<const int,pair<int,int> >::iterator it = faceEdgePairMap.begin(); it != faceEdgePairMap.end(); ++it) {
      //cout << it->second.first << " " << it->second.second << endl;
      assert(it->second.first != -1);
      assert(it->second.second != -1);
    }
    //cout << "size: " << faceEdgePairMap.size() << endl;
    
    // build the walk...
    if (faceEdgePairMap.size() > 0) {
      map<const int,pair<int,int> >::iterator it = faceEdgePairMap.begin();
      int fnoc0_local_first = it->second.first;
      int fnoc0_local = it->second.first;
      int fnoc1_local = it->second.second;
      const int nfa_intersected = faceEdgePairMap.size();
      int iter = 0;
      while (1) {
        //cout << "1: " << fnoc0_local << " " << fnoc1_local << " " << fnoc0_local_first << endl; cout.flush();
        iter++;
        FOR_I3 points.push_back(unordered_points[3*fnoc0_local+i]);
        const int fnoc1 = fnoc1_local+fpnpocv_i[icv];
        const int ifa_l = fpnpocv_v[fnoc1][0];
        const int ifa_r = fpnpocv_v[fnoc1][1];

        // next face...
        if (ifa_l != it->first) {
          //cout << ifa_l << " " << ifa_r << " " << it->first << endl;
          assert(ifa_r == it->first);
          it = faceEdgePairMap.find(ifa_l);
          assert(it != faceEdgePairMap.end());
        }
        else {
          assert(ifa_l == it->first);
          it = faceEdgePairMap.find(ifa_r);
          assert(it != faceEdgePairMap.end());
        }

        // next edge... 
        fnoc0_local = fnoc1_local;
        if ( it->second.first == fnoc1_local ) {
          fnoc1_local = it->second.second;
        }
        else { 
          assert( it->second.second == fnoc1_local );
          fnoc1_local = it->second.first;
        }
        //cout << "2: " << fnoc0_local << " " << fnoc1_local << " " << fnoc0_local_first << endl;

        if (fnoc0_local == fnoc0_local_first)
          break;
      }
      assert(iter == nfa_intersected);
    }
    // check...
    /*
    {
      cout << vof[icv] << " " << COUT_VEC(n[icv]) << endl;
      list<double>::iterator it = plicPointList[icv].begin();
      while (it != plicPointList[icv].end()) {
        FOR_I3 {
          cout << *it << ",";
          ++it;
        }
        cout << endl;
      }
    }
    */
  }
  void calcCurvatureParaboloidFittingMoments() {

    COUT1("calcCurvature()...");
    assert(mpi_size == 1); // don't have communicator for plic yet... add one if you get this to work

    buildPlicPolys();
    
    // calc kappa by fitting paraboloid...

    double A[6][6], p[6], b[6];
    double R[3][3], t[3];
    FOR_ICV {

      kappa[icv] = 0.0; 
      
      if (cv_flag[icv] >= 1) {

        // count number of well-defined nbr n's 
        int np = 0;
        for (int coc = cvocv_i[icv]; coc != cvocv_i[icv+1]; ++coc) {
          const int icv_nbr = cvocv_v[coc];
          if (cv_flag[icv_nbr] >= 1 && (DOT_PRODUCT(n[icv],n[icv_nbr]) > 0.0) ) {
            np++;
          }
        }
        if (np >= 6) {

          // transform coordinate system to be aligned with n and centered on plic centroid...

          buildTransform(R,t,icv);

          FOR_I6 FOR_J6 A[i][j] = 0.0;
          FOR_I6 b[i] = 0.0;

          for (int coc = cvocv_i[icv]; coc != cvocv_i[icv+1]; ++coc) {
            const int icv_nbr = cvocv_v[coc];
            if (cv_flag[icv_nbr] >= 1 && (DOT_PRODUCT(n[icv],n[icv_nbr]) > 0.0) ) {
              //cout << "icv " << icv << endl;
              //writePlicOFF(icv_nbr);

              // transform interface point and normal...
              double xp[3]; FOR_I3 xp[i] = plic_xc[icv_nbr][i];
              //double xp[3]; FOR_I3 xp[i] = g[icv_nbr]*n[icv_nbr][i]+x_cv[icv_nbr][i];
              //cout << "x_int: " << COUT_VEC(xp) << " -> ";
              applyTransform(xp,R,t);
              //cout << COUT_VEC(xp) << endl;;
              FOR_I3 xp[i] /= r_vv[icv];
              double np[3]; FOR_I3 np[i] = n[icv_nbr][i];
              //cout << "n_int: " << COUT_VEC(np) << " -> ";
              applyRotation(np,R);
              //cout << COUT_VEC(np) << endl;;
              double ds[3] = {DOT_PRODUCT(np,xp)/np[0],-np[1]/np[0],-np[2]/np[0]};
              //cout << "ds: " << COUT_VEC(ds) << endl;

              double as[6] = {0.0,0.0,0.0,0.0,0.0,0.0};
              assert(plicPointList[icv_nbr].size()%3 == 0);
              //cout << plicPointList[icv_nbr].size() << endl;
              if (plicPointList[icv_nbr].size() > 0) {
                list<double>::iterator it = plicPointList[icv_nbr].begin();
                double x0[3],x1[3]; 
                FOR_I3 {
                  x1[i] = *it;
                  ++it;
                }
                applyTransform(x1,R,t);
                FOR_I3 x1[i] /= r_vv[icv];
                //cout << COUT_VEC(x1) << endl;
                while (it != plicPointList[icv_nbr].end()) {
                  FOR_I3 {
                    x0[i] = x1[i];
                    x1[i] = *it;
                    ++it;
                  }
                  applyTransform(x1,R,t);
                  //cout << COUT_VEC(x1) << endl;
                  FOR_I3 x1[i] /= r_vv[icv];
                  const double area_2 = x0[1]*x1[2]-x1[1]*x0[2];
                  //cout << "area_2: " << area_2 << endl;
                  as[0] += area_2;                                                         // int_plic 1 dA
                  as[1] += area_2*(x0[1]+x1[1]);                                           // int_plic x dA
                  as[2] += area_2*(x0[2]+x1[2]);                                           // int_plic y dA
                  as[3] += area_2*(x0[1]*x0[1]+x0[1]*x1[1]+x1[1]*x1[1]);                   // int_plic x^2 dA
                  as[4] += area_2*(2.0*(x0[1]*x0[2]+x1[1]*x1[2])+x0[1]*x1[2]+x1[1]*x0[2]); // int_plic xy dA
                  as[5] += area_2*(x0[2]*x0[2]+x0[2]*x1[2]+x1[2]*x1[2]);                   // int_plic y^2 dA
                  //FOR_I6 cout << as[i] << " ";
                }
                it = plicPointList[icv].begin();
                FOR_I3 {
                  x0[i] = x1[i];
                  x1[i] = *it;
                  ++it;
                }
                applyTransform(x1,R,t);
                //cout << COUT_VEC(x1) << endl;
                FOR_I3 x1[i] /= r_vv[icv];
                const double area_2 = x0[1]*x1[2]-x1[1]*x0[2];
                //cout << "area_2: " << area_2 << endl;
                as[0] += area_2;                                                         // int_plic 1 dA
                as[1] += area_2*(x0[1]+x1[1]);                                           // int_plic x dA
                as[2] += area_2*(x0[2]+x1[2]);                                           // int_plic y dA
                as[3] += area_2*(x0[1]*x0[1]+x0[1]*x1[1]+x1[1]*x1[1]);                   // int_plic x^2 dA
                as[4] += area_2*(2.0*(x0[1]*x0[2]+x1[1]*x1[2])+x0[1]*x1[2]+x1[1]*x0[2]); // int_plic xy dA
                as[5] += area_2*(x0[2]*x0[2]+x0[2]*x1[2]+x1[2]*x1[2]);                   // int_plic y^2 dA
                //cout << endl;
                //getchar();
                as[0] /= 2.0;
                as[1] /= 6.0;
                as[2] /= 6.0;
                as[3] /= 12.0;
                as[4] /= 24.0;
                as[5] /= 12.0;
                const double bs = DOT_PRODUCT(ds,as);
                //cout << "bs: " << endl;
                //cout << ds[0] << " " << ds[1] << " " << ds[2] << " " << bs << endl;

                // add samples to matrix and rhs...
                const double wgt = 1.0;
                //const double wgt = vof[icv_nbr]*(1.0-vof[icv_nbr]);
                //const double wgt = plic_area[icv_nbr];
                FOR_I6 FOR_J6 A[i][j] += wgt*as[i]*as[j];
                FOR_J6 b[j] += wgt*bs*as[j];
                //cout << "A: " << endl;
                //FOR_I6 {
                //  FOR_J6 cout << A[i][j] << " ";
                //  cout << endl;
                //}
                //cout << "b: " << endl;
                //FOR_J6 cout << b[j]<< endl;
              }
            }
          }
          //getchar();

          // solve for polynomial coefficients: z = p0*x^2+p1*x*y+p2*y^2+p3*x+p4*y+p5...
          //solveGaussElimationWithPivoting(p,A,b);
          solveSvd(p,A,b);
        //  FOR_I6 cout << " " << p[i]; cout << endl;

          // calc kappa...
          const double den = 1.0+p[1]*p[1]+p[2]*p[2];
          const double this_kappa = -2.0*(p[3]*(1.0+p[2]*p[2])+p[5]*(1.0+p[1]*p[1])-p[4]*p[1]*p[2])/sqrt(den*den*den);
          if (this_kappa == this_kappa) 
            kappa[icv] = this_kappa/r_vv[icv];

         // cout << "np kappa " << np << " " <<  this_kappa << endl;
        }
      }
    }
    updateCvData(kappa);
    dumpRange(kappa,ncv,"kappa");
    
    // when cell nears empty or full, the accuracy of this method degrades, so average...

    double *kappa0 = new double[ncv_g];
    FOR_ICV_G kappa0[icv] = kappa[icv];
    FOR_ICV {
      if (cv_flag[icv] >= 1.0) {
        // count number of well-defined nbr n's 
        int np = 0;
        for (int coc = cvocv_i[icv]; coc != cvocv_i[icv+1]; ++coc) {
          const int icv_nbr = cvocv_v[coc];
          if (cv_flag[icv_nbr] >= 1 && (DOT_PRODUCT(n[icv],n[icv_nbr]) > 0.0) ) {
            np++;
          }
        }
        if ( vof[icv] < 0.1 || vof[icv] > 0.9 || np < 6 ) { 
          double wgt_sum = 0.0;
          double new_kappa = 0.0;
          for (int coc = cvocv_i[icv]+1; coc != cvocv_i[icv+1]; ++coc) {
            const int icv_nbr = cvocv_v[coc];
            if (cv_flag[icv_nbr] >= 1) {
              const double wgt = 1.0;
              //const double wgt = plic_area[icv_nbr];
              new_kappa += wgt*kappa0[icv_nbr];
              wgt_sum += wgt;
            }
          }
          if (wgt_sum > 0.0) 
            kappa[icv] = new_kappa/wgt_sum;
        }
      }
    }
    delete[] kappa0;
  
    // bound by cell length scale...
    FOR_ICV {
      if (kappa[icv] > 0.0) {
        kappa[icv] = min(kappa[icv],1.0/r_vv[icv]);
      }
      else if (kappa[icv] < 0.0) {
        kappa[icv] = max(kappa[icv],-1.0/r_vv[icv]);
      }
    }
    updateCvData(kappa);

  }

  void calcCurvatureParaboloidFittingCentroids() {

    COUT1("calcCurvature()...");

    buildPlicPolys();
    
    // calc kappa by fitting paraboloid...

    double A[6][6], p[6], b[6];
    double R[3][3], t[3];
    FOR_ICV {

      kappa[icv] = 0.0; 
      
      if (cv_flag[icv] >= 1) {

        // count number of well-defined nbr n's 
        int np = 0;
        for (int coc2 = cvocv2_i[icv]; coc2 != cvocv2_i[icv+1]; ++coc2) {
          const int icv_nbr = cvocv2_v[coc2];
          if (cv_flag[icv_nbr] >= 1 && (DOT_PRODUCT(n[icv],n[icv_nbr]) > 0.0) ) {
            np++;
          }
        }
        if (np >= 6) {

          // transform coordinate system to be aligned with n and centered on x...

          //buildTransform(R,t,g[icv],n[icv],x_cv[icv]);
          buildTransform(R,t,icv);

          FOR_I6 FOR_J6 A[i][j] = 0.0;
          FOR_I6 b[i] = 0.0;
          for (int coc2 = cvocv2_i[icv]; coc2 != cvocv2_i[icv+1]; ++coc2) {
            const int icv_nbr = cvocv2_v[coc2];
            if (cv_flag[icv_nbr] >= 1 && (DOT_PRODUCT(n[icv],n[icv_nbr]) > 0.0) ) {

              //double xp[3]; FOR_I3 xp[i] = g[icv_nbr]*n[icv_nbr][i]+x_cv[icv_nbr][i];
              double xp[3]; FOR_I3 xp[i] = plic_xc[icv_nbr][i];

              // transform point...
              applyTransform(xp,R,t);
              FOR_I3 xp[i] /= r_vv[icv];

              // weight point...
              // f*(1-f) ~ 5/8*A^3 + 3/16*A^2, so you could depress this weight w/ sqrt or cube root to make it grow less w/ PLIC A
              //const double wgt = 1.0;
              const double wgt = plic_area[icv_nbr];
              //const double wgt = sqrt(vof[icv_nbr]*(1.0-vof[icv_nbr]));
              ///const double wgt = vof[icv_nbr]*(1.0-vof[icv_nbr]);

              // add samples to matrix and rhs...
              const double x = xp[1];
              const double x2 = x*x;
              const double x3 = x2*x;
              const double x4 = x3*x;
              const double y = xp[2];
              const double y2 = y*y;
              const double y3 = y2*y;
              const double y4 = y3*y;
              const double xy = x*y;
              const double x2y = x*xy;
              const double x3y = x*x2y;
              const double x2y2 = x2y*y;
              const double xy2 = xy*y;
              const double xy3 = xy2*y;
              const double z = xp[0];
              A[0][0] += wgt*x4;   A[0][1] += wgt*x3y;  A[0][2] += wgt*x2y2; A[0][3] += wgt*x3;  A[0][4] += wgt*x2y; A[0][5] += wgt*x2;  b[0] += wgt*z*x2;
              A[1][0] += wgt*x3y;  A[1][1] += wgt*x2y2; A[1][2] += wgt*xy3;  A[1][3] += wgt*x2y; A[1][4] += wgt*xy2; A[1][5] += wgt*xy;  b[1] += wgt*z*xy;
              A[2][0] += wgt*x2y2; A[2][1] += wgt*xy3;  A[2][2] += wgt*y4;   A[2][3] += wgt*xy2; A[2][4] += wgt*y3;  A[2][5] += wgt*y2;  b[2] += wgt*z*y2;
              A[3][0] += wgt*x3;   A[3][1] += wgt*x2y;  A[3][2] += wgt*xy2;  A[3][3] += wgt*x2;  A[3][4] += wgt*xy;  A[3][5] += wgt*x;   b[3] += wgt*z*x;
              A[4][0] += wgt*x2y;  A[4][1] += wgt*xy2;  A[4][2] += wgt*y3;   A[4][3] += wgt*xy;  A[4][4] += wgt*y2;  A[4][5] += wgt*y;   b[4] += wgt*z*y;
              A[5][0] += wgt*x2;   A[5][1] += wgt*xy;   A[5][2] += wgt*y2;   A[5][3] += wgt*x;   A[5][4] += wgt*y;   A[5][5] += wgt*1.0; b[5] += wgt*z;
            }
          }

          // solve for polynomial coefficients: z = p0*x^2+p1*x*y+p2*y^2+p3*x+p4*y+p5...
          //solveGaussElimationWithPivoting(p,A,b);
          solveSvd(p,A,b);

          // calc kappa...
          const double den = 1.0+p[3]*p[3]+p[4]*p[4];
          const double this_kappa = -2.0*(p[0]*(1.0+p[4]*p[4])+p[2]*(1.0+p[3]*p[3])-p[1]*p[3]*p[4])/sqrt(den*den*den);
          if (this_kappa == this_kappa) 
            kappa[icv] = this_kappa/r_vv[icv];

        }
      }
    }
    updateCvData(kappa);

    // when cell nears empty or full, or the normal is angled close to 45, the accuracy of this method degrades, so average...

    double *kappa0 = new double[ncv_g];
    FOR_ICV_G kappa0[icv] = kappa[icv];
    FOR_ICV {
      if (cv_flag[icv] >= 1.0) {
        //double wgt_avg = 0.0;
        int np = 0;
        for (int coc = cvocv_i[icv]; coc != cvocv_i[icv+1]; ++coc) {
          const int icv_nbr = cvocv_v[coc];
          if (cv_flag[icv_nbr] >= 1 && (DOT_PRODUCT(n[icv],n[icv_nbr]) > 0.0) ) {
            //wgt_avg += vof[icv_nbr]*(1.0-vof[icv_nbr]);
            np++;
          }
        }
        //wgt_avg /= (double)np;
        //const double theta = acos(MAX3(abs(n[icv][0]),abs(n[icv][1]),abs(n[icv][2])));
        //cout << "HELP: vof,theta,wgt_avg,kappa: " << vof[icv] << ", " << theta << ", " << wgt_avg << ", " << kappa[icv] << endl; 
        if ( vof[icv] < 0.1 || vof[icv] > 0.9 || np < 6) { //|| theta > 0.55 || wgt_avg < 0.16) {
          double wgt_sum = 0.0;
          double new_kappa = 0.0;
          for (int coc = cvocv_i[icv]+1; coc != cvocv_i[icv+1]; ++coc) {
            const int icv_nbr = cvocv_v[coc];
            if (cv_flag[icv_nbr] >= 1) {
            //const double theta = acos(MAX3(abs(n[icv][0]),abs(n[icv][1]),abs(n[icv][2])));
            //if (cv_flag[icv_nbr] >= 1 && (DOT_PRODUCT(n[icv],n[icv_nbr]) > 0.0) ) {
              //if ( vof[icv] > 0.2 || vof[icv] < 0.8 || theta < 0.55 ) {
                const double wgt = plic_area[icv_nbr];
                new_kappa += wgt*kappa0[icv_nbr];
                wgt_sum += wgt;
              //}
            }
          }
          if (wgt_sum > 0.0) 
            kappa[icv] = new_kappa/wgt_sum;
        }
      }
    }
    delete[] kappa0;
    
    // bound by cell length scale...
    FOR_ICV {
      if (kappa[icv] > 0.0) {
        kappa[icv] = min(kappa[icv],1.0/r_vv[icv]);
      }
      else if (kappa[icv] < 0.0) {
        kappa[icv] = max(kappa[icv],-1.0/r_vv[icv]);
      }
    }
    updateCvData(kappa);

    dumpRange(kappa,ncv,"kappa");

  }

