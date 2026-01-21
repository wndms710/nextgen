#include "SplineStuff.hpp"

namespace SplineStuff {
  
  void CubicSpline::init(const double (*xp)[3],const int np,const bool b_periodic) {

    double (*rhs)[3] = new double[np][3];

    if (b_periodic) {

      this->b_periodic = true;

      double *dx = new double[np];
      for (int ip = 0; ip < np-1; ++ip)
        dx[ip] = DIST(xp[ip],xp[ip+1]);
      dx[np-1] = DIST(xp[0],xp[np-1]);

      FOR_I3 rhs[0][i] = (xp[1][i]-xp[0][i])/dx[0] - (xp[0][i]-xp[np-1][i])/dx[np-1];
      for (int ip = 1; ip < np-1; ++ip)
        FOR_I3 rhs[ip][i] = (xp[ip+1][i]-xp[ip][i])/dx[ip] - (xp[ip][i]-xp[ip-1][i])/dx[ip-1];
      FOR_I3 rhs[np-1][i] = (xp[0][i]-xp[np-1][i])/dx[np-1] - (xp[np-1][i]-xp[np-2][i])/dx[np-2];
      
      // This uses the Sherman-Morrison formula to split out the cyclic components using a decomp: A -> A'+ u(x)v
      // u_0 = gamma, u_n-1 = c_n-1, v_0 = 1, v_n-1 = a_0/gamma, other terms 0
      
      const double gamma = -(dx[np-1]+dx[0])/3.0; // free parameter, -A_00 recomended
      
      // note that A'_00 = A_00 - gamma and A'_n-1,n-1 = A_n-1,n-1 - a_0*c_n-1/gamma

      // solve A' z = u and A' y = rhs using Thomas algo...

      double *cp = new double[np-1];
      double *z = new double[np];
      assert(dp == NULL); 
      dp = new double[np][3]; // holds y for now
      
      // forward sweep to get cp and dp...
      cp[0] = dx[0]/6.0/((dx[np-1]+dx[0])/3.0-gamma);
      z[0] = gamma/((dx[np-1]+dx[0])/3.0-gamma);
      FOR_I3 dp[0][i] = rhs[0][i]/((dx[np-1]+dx[0])/3.0-gamma);
      for (int ip = 1; ip < np-1; ++ip) {
        cp[ip] = dx[ip]/6.0/((dx[ip]+dx[ip-1])/3.0-dx[ip-1]/6.0*cp[ip-1]);
        z[ip] = (0.0-dx[ip-1]/6.0*z[ip-1])/((dx[ip]+dx[ip-1])/3.0-dx[ip-1]/6.0*cp[ip-1]);
        FOR_I3 dp[ip][i] = (rhs[ip][i]-dx[ip-1]/6.0*dp[ip-1][i])/((dx[ip]+dx[ip-1])/3.0-dx[ip-1]/6.0*cp[ip-1]);
      }
      z[np-1] = (dx[np-1]/6.0-dx[np-2]/6.0*z[np-2])/((dx[np-1]+dx[np-2])/3.0-dx[np-1]/6.0*dx[np-1]/6.0/gamma-dx[np-2]/6.0*cp[np-2]);
      FOR_I3 dp[np-1][i] = (rhs[np-1][i]-dx[np-2]/6.0*dp[np-2][i])/((dx[np-1]+dx[np-2])/3.0-dx[np-1]/6.0*dx[np-1]/6.0/gamma-dx[np-2]/6.0*cp[np-2]);

      // back substitution (right into dp)...
      for (int ip = np-2; ip >= 0; --ip) {
        z[ip] -= cp[ip]*z[ip+1];
        FOR_I3 dp[ip][i] -= cp[ip]*dp[ip+1][i];
      }
      
      // finish up with dp = y - (v.y)/(1+v.z) z...
      double vy[3]; 
      FOR_I3 vy[i] = dp[0][i]+dx[np-1]/6.0*dp[np-1][i]/gamma; 
      const double vz = z[0]+dx[np-1]/6.0*z[np-1]/gamma;  
      assert(vz != -1.0);
      FOR_I3 vy[i] /= 1.0+vz;
      for (int ip = 0; ip <= np-1; ++ip)
        FOR_I3 dp[ip][i] -= vy[i]*z[ip];

      delete[] z;
      delete[] cp;
      delete[] dx;

    }
    else {

      // notes: this 

      //double alpha = 1.0;

      double *dx = new double[np-1];
      for (int ip = 0; ip < np-1; ++ip)
        dx[ip] = DIST(xp[ip],xp[ip+1]);
      
      //FOR_I3 rhs[0][i] = -(xp[1][i]-xp[0][i])/dx[0];
      // endpoint "runout"...
      FOR_I3 rhs[0][i] = 0.0;

      for (int ip = 1; ip < np-1; ++ip)
        FOR_I3 rhs[ip][i] = (xp[ip+1][i]-xp[ip][i])/dx[ip] - (xp[ip][i]-xp[ip-1][i])/dx[ip-1];
      
      //FOR_I3 rhs[np-1][i] = (xp[np-1][i]-xp[np-2][i])/dx[np-2];
      // endpoint "runout"...
      FOR_I3 rhs[np-1][i] = 0.0;
      
      //const double a = 1.0, b = 4.0, c = 1.0;   
      double *cp = new double[np-1];
      assert(dp == NULL); 
      dp = new double[np][3];
      
      // forward sweep to get cp and dp...
      //cp[0] = 0.5/(1.0+alpha); // 2.0*c/b
      // runout hack...
      cp[0] = 0.0;
      FOR_I3 dp[0][i] = 3.0*rhs[0][i]/dx[0];
      for (int ip = 1; ip < np-1; ++ip) {
        cp[ip] = dx[ip]/(2.0*(dx[ip]+dx[ip-1])-dx[ip-1]*cp[ip-1]);
        FOR_I3 dp[ip][i] = (6.0*rhs[ip][i]-dx[ip-1]*dp[ip-1][i])/(2.0*(dx[ip]+dx[ip-1])-dx[ip-1]*cp[ip-1]);
      }
      //FOR_I3 dp[np-1][i] = (6.0*rhs[np-1][i]-dx[np-2]*dp[np-2][i])/(2.0*dx[np-2]-dx[np-2]*cp[np-2]);
      // runout hack...
      FOR_I3 dp[np-1][i] = 0.0;
      
      delete[] dx;

      // back substitution (right into dp)...
      for (int ip = np-2; ip >= 0; --ip) 
        FOR_I3 dp[ip][i] -= cp[ip]*dp[ip+1][i];
      
      delete[] cp;
    }
      
    // copy the points. reuse the rhs memory here...
    assert(this->xp == NULL);
    this->xp = rhs;
    for (int ip = 0; ip < np; ++ip)
      FOR_I3 this->xp[ip][i] = xp[ip][i];
    this->np = np;

    // compute the length of each interval using quadrature...
    const double tmp1 = sqrt(5.0-2.0*sqrt(10.0/7.0))/6.0;
    const double tmp2 = sqrt(5.0+2.0*sqrt(10.0/7.0))/6.0;
    /*
    const double s_gauss5[5] = {
      0.5,
      0.5-tmp1,
      0.5+tmp1,
      0.5-tmp2,
      0.5+tmp2
    };
    */
    s_gauss5[0] = 0.5;
    s_gauss5[1] = 0.5-tmp1;
    s_gauss5[2] = 0.5+tmp1;
    s_gauss5[3] = 0.5-tmp2;
    s_gauss5[4] = 0.5+tmp2;    
    
    const double tmp3 = (322.0+13.0*sqrt(70.0))/1800.0;    
    const double tmp4 = (322.0-13.0*sqrt(70.0))/1800.0;    
    /*
    const double wgt_gauss5[5] = { 
      0.5*128.0/225.0,
      tmp3,
      tmp3,
      tmp4,
      tmp4
    };
    */

    wgt_gauss5[0] = 0.5*128.0/225.0;
    wgt_gauss5[1] = tmp3;
    wgt_gauss5[2] = tmp3;
    wgt_gauss5[3] = tmp4;
    wgt_gauss5[4] = tmp4;    
    
    assert(sp == NULL);
    if (b_periodic) {
      sp = new double[np+1];
      sp[0] = 0.0;
      double dxdt[3];
      for (int ip0 = 0; ip0 < np; ++ip0) {
        double l = 0.0;
        FOR_I5 {
          getDxdt(dxdt,s_gauss5[i],ip0);
          l += wgt_gauss5[i]*MAG(dxdt);
        }
        sp[ip0+1] = sp[ip0] + l;
      }
      // set the member length...
      length = sp[np];
      
      // finally, normalize "sp"...
      for (int ip = 1; ip < np; ++ip)
        sp[ip] /= length;
      sp[np] = 1.0;
    }
    else { 
      sp = new double[np];
      sp[0] = 0.0;
      double dxdt[3];
      for (int ip0 = 0; ip0 < np-1; ++ip0) {
        double l = 0.0;
        FOR_I5 {
          getDxdt(dxdt,s_gauss5[i],ip0);
          l += wgt_gauss5[i]*MAG(dxdt);
        }
        sp[ip0+1] = sp[ip0] + l;
      }
      // set the member length...
      length = sp[np-1];
      
      // finally, normalize "sp"...
      for (int ip = 1; ip < np-1; ++ip)
        sp[ip] /= length;
      sp[np-1] = 1.0;
    }
    
    //cout << " > length: " << length << endl;
    
  }

  void CubicSpline::getXp(double xpi[3],const int i) const {
    assert(xp);
    assert((i >= 0)&&(i < np));
    FOR_J3 xpi[j] = xp[i][j];
  }

  double CubicSpline::getSp(const int i) const {
    assert(sp);
    assert((i >= 0)&&(i < np));
    return sp[i];
  }
  
  void CubicSpline::getX(double x[3],const double s) const {
    if (s <= 0.0) {
      if (s < 0.0) cout << "Warning: CubicSpline::getX has s out of range: " << s << endl;
      FOR_I3 x[i] = xp[0][i];
    }
    else if ((s >= 1.0)&&(!b_periodic)) {
      if (s > 1.0) cout << "Warning: CubicSpline::getX has s out of range: " << s << endl;
      FOR_I3 x[i] = xp[np-1][i];
    }
    else {
      int left = 0;
      int right = np-1;
      if (b_periodic)
        ++right;
      while ((right - left) > 1) {
        const int middle = (left + right)/2;   // equivalent to floor..
        if (s >= sp[middle])
          left = middle;
        else
          right = middle;
      }
      // left is the final interval...
      const double starget = (s-sp[left])/(sp[left+1]-sp[left]);

      // TODO: it is a crying shame to discard this very good initial guess
      // and just do bisection. If this becomes a rate-limiting step for use,
      // move to a better solver...
      //double t = starget;
      const double tol = 1.0E-12;

      double t0 = 0.0;
      double t1 = 1.0;
      double dxdt[3],l0,l1;
      while (t1-t0 > tol) {
        const double t = 0.5*(t0+t1);
        // get the length on either side of this point...
        l0 = 0.0;
        l1 = 0.0;
        FOR_I5 {
          getDxdt(dxdt,s_gauss5[i]*t,left);
          l0 += wgt_gauss5[i]*MAG(dxdt);
          getDxdt(dxdt,t+s_gauss5[i]*(1.0-t),left);
          l1 += wgt_gauss5[i]*MAG(dxdt);
        }
        l0 *= t;
        l1 *= (1.0-t);
        if (l0/(l0+l1) < starget) 
          t0 = t;
        else
          t1 = t;
      }
      // return the x at this converged "t"...
      getX(x,0.5*(t0+t1),left);
      
    }
  }

  void CubicSpline::getXAndDx(double x[3],double dx[3],const double s) const {
    if (s <= 0.0) {
      if (s < 0.0) cout << "Warning: CubicSpline::getX has s out of range: " << s << endl;
      FOR_I3 x[i] = xp[0][i];
      getDxdt(dx,0.0,0);
    }
    else if ((s >= 1.0)&&(!b_periodic)) {
      if (s > 1.0) cout << "Warning: CubicSpline::getX has s out of range: " << s << endl;
      FOR_I3 x[i] = xp[np-1][i];
      getDxdt(dx,1.0,np-1);
    }
    else {
      int left = 0;
      int right = np-1;
      if (b_periodic)
        ++right;
      while ((right - left) > 1) {
        const int middle = (left + right)/2;   // equivalent to floor..
        if (s >= sp[middle])
          left = middle;
        else
          right = middle;
      }
      // left is the final interval...
      const double starget = (s-sp[left])/(sp[left+1]-sp[left]);

      // TODO: it is a crying shame to discard this very good initial guess
      // and just do bisection. If this becomes a rate-limiting step for use,
      // move to a better solver...
      //double t = starget;
      const double tol = 1.0E-12;

      double t0 = 0.0;
      double t1 = 1.0;
      double dxdt[3],l0,l1;
      while (t1-t0 > tol) {
        const double t = 0.5*(t0+t1);
        // get the length on either side of this point...
        l0 = 0.0;
        l1 = 0.0;
        FOR_I5 {
          getDxdt(dxdt,s_gauss5[i]*t,left);
          l0 += wgt_gauss5[i]*MAG(dxdt);
          getDxdt(dxdt,t+s_gauss5[i]*(1.0-t),left);
          l1 += wgt_gauss5[i]*MAG(dxdt);
        }
        l0 *= t;
        l1 *= (1.0-t);
        if (l0/(l0+l1) < starget) 
          t0 = t;
        else
          t1 = t;
      }
      // return the x and dx at this converged "t"...
      getX(x,0.5*(t0+t1),left);
      getDxdt(dx,0.5*(t0+t1),left);
      
    }
  }


  
  void CubicSpline::getX(double x[3],const double t,const int ip0) const {
    if (b_periodic) {
      assert(ip0 <= np-1);
      if (ip0 < np-1) {
        const double dx = DIST(xp[ip0],xp[ip0+1]);
        FOR_I3 x[i] = dx*dx/6.0*(dp[ip0][i]*((1.0-t)*(1.0-t)*(1.0-t)-(1.0-t)) + dp[ip0+1][i]*(t*t*t-t)) + xp[ip0][i]*(1.0-t) + xp[ip0+1][i]*t;
      }
      else {
        const double dx = DIST(xp[np-1],xp[0]);
        FOR_I3 x[i] = dx*dx/6.0*(dp[np-1][i]*((1.0-t)*(1.0-t)*(1.0-t)-(1.0-t)) + dp[0][i]*(t*t*t-t)) + xp[np-1][i]*(1.0-t) + xp[0][i]*t;
      }
    }
    else {
      assert(ip0 < np-1);
      const double dx = DIST(xp[ip0],xp[ip0+1]);
      FOR_I3 x[i] = dx*dx/6.0*(dp[ip0][i]*((1.0-t)*(1.0-t)*(1.0-t)-(1.0-t)) + dp[ip0+1][i]*(t*t*t-t)) + xp[ip0][i]*(1.0-t) + xp[ip0+1][i]*t;
    }
  }
  
  void CubicSpline::getDxdt(double dxdt[3],const double t,const int ip0) const {
    if (b_periodic) {
      assert(ip0 <= np-1);
      if (ip0 < np-1) {
        const double dx = DIST(xp[ip0],xp[ip0+1]);
        FOR_I3 dxdt[i] = dx*dx/6.0*(dp[ip0][i]*(6.0*t-3.0*t*t-2.0) + dp[ip0+1][i]*(3.0*t*t-1.0)) - xp[ip0][i] + xp[ip0+1][i];
      }
      else {
        const double dx = DIST(xp[np-1],xp[0]);
        FOR_I3 dxdt[i] = dx*dx/6.0*(dp[np-1][i]*(6.0*t-3.0*t*t-2.0) + dp[0][i]*(3.0*t*t-1.0)) - xp[np-1][i] + xp[0][i];
      }
    }
    else {
      assert(ip0 < np-1);
      const double dx = DIST(xp[ip0],xp[ip0+1]);
      FOR_I3 dxdt[i] = dx*dx/6.0*(dp[ip0][i]*(6.0*t-3.0*t*t-2.0) + dp[ip0+1][i]*(3.0*t*t-1.0)) - xp[ip0][i] + xp[ip0+1][i];
    }
  }
  
  void CubicSplineXy::init(const double *x, const double *y, const int n) {

    this->np = n;
    this->xp = new double [n];
    this->yp = new double [n];
    this->dx = new double [n-1];
    this->dp = new double [n];
    for (int i = 0; i < n; i++) { 
      this->xp[i] = x[i];
      this->yp[i] = y[i];
    }

    if (n == 2) {
      for (int i = 0; i < n; i++) {
        this->xp[i] = x[i];
        this->yp[i] = y[i];
        this->dp[i] = 0.0;
      }
      dx[0] = xp[1] - xp[0];
      return;
    }

    double * rhs = new double [n-2];
    double * cp = new double [n-3];

    double * a = new double [n-2];
    double * b = new double [n-2];
    double * c = new double [n-2];

    dp[0] = 0.0;
    dp[n-1] = 0.0;
    dx[0] = x[1] - x[0];
    dx[n-2] = x[n-1] - x[n-2];
    for (int i = 0; i < n-2; i++) {
      dx[i+1] = x[i+2] - x[i+1];
      rhs[i] = (y[i+2]-y[i+1])/dx[i+1] - (y[i+1]-y[i])/dx[i];
      b[i] = (dx[i]+dx[i+1])/3.0;
      a[i] = (dx[i])/6.0;
      c[i] = (dx[i+1])/6.0;
    }
    a[0] = 0.0;
    c[n-3] = 0.0;

    cp[0] = c[0]/b[0];
    dp[1] = rhs[0]/b[0];
    for (int i = 1; i < n-3; i++) {
      cp[i] = c[i]/(b[i]-a[i]*cp[i-1]);
      dp[i+1] = (rhs[i]-a[i]*dp[i])/(b[i]-a[i]*cp[i-1]);
    }
    dp[n-2] = (rhs[n-3]-a[n-3]*dp[n-3])/(b[n-3]-a[n-3]*cp[n-4]);

    for (int i = n-3; i > 0; i--) {
      dp[i] -= cp[i-1]*dp[i+1];
    }
      
    delete [] rhs;
    delete [] cp;
    delete [] a;
    delete [] b;
    delete [] c;
  }

  int CubicSplineXy::getIp(const double x) const {
    int left = 0;
    int right = np-1;
    while ((right - left)>1) {
      int middle = (left + right)/2;
      if (x > xp[middle]) {
        left = middle;
      } else {
        right = middle;
      }
    }
    return left;
  }

  double CubicSplineXy::getY(const double x) const {
    if (x <= xp[0]) {
      return yp[0];
    } else if (x >= xp[np-1]) {
      return yp[np-1];
    }

    int ip = getIp(x);
    double ret = dp[ip]/6.*( (xp[ip+1]-x)*(xp[ip+1]-x)*(xp[ip+1]-x)/dx[ip] - dx[ip]*(xp[ip+1]-x) ) +
                 dp[ip+1]/6.*( (x-xp[ip])*(x-xp[ip])*(x-xp[ip])/dx[ip] - dx[ip]*(x-xp[ip])) +
                 yp[ip]*(xp[ip+1]-x)/dx[ip] + yp[ip+1]*(x-xp[ip])/dx[ip];
    return ret;
  }

}
