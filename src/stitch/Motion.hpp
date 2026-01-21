#ifndef _MOTION_HPP_
#define _MOTION_HPP_

class AbstractMotion {
public:
  virtual ~AbstractMotion() {}
  virtual void move(double (*x)[3],const int n,const double time,const double dt) = 0;
  virtual void move0(double (*x)[3],const int n,const double time) { }
  virtual void setCentroid(const double (*x)) = 0; // we should have a default implementation of this
};

class ConstantMotion : public AbstractMotion {
private:
  double v[3];
public:
  ConstantMotion(const double v[3]) {
    FOR_I3 this->v[i] = v[i];
  }
  ConstantMotion(const double vx,const double vy,const double vz) {
    this->v[0] = vx;
    this->v[1] = vy;
    this->v[2] = vz;
  }
  void move(double (*x)[3],const int n,const double time,const double dt) {
    for (int i = 0; i < n; ++i) {
      FOR_J3 x[i][j] += dt*v[j];
    }
  }
  void setCentroid(const double (*x)) {
  }
};

class ConstantOmegaZMotion : public AbstractMotion {
private:
  double v[3];
  double omega_z;
  double xc[3];
public:
  ConstantOmegaZMotion(const double v[3],const double omega_z) {
    FOR_I3 this->v[i] = v[i];
    this->omega_z = omega_z;
    FOR_I3 xc[i] = 0.0;
  }
  ConstantOmegaZMotion(const double vx,const double vy,const double vz,const double omega_z) {
    this->v[0] = vx;
    this->v[1] = vy;
    this->v[2] = vz;
    this->omega_z = omega_z;
    FOR_I3 xc[i] = 0.0;
  }
  void move(double (*x)[3],const int n,const double time,const double dt) {
    const double dtheta = omega_z*dt;
    const double cos_dtheta = cos(dtheta);
    const double sin_dtheta = sin(dtheta);
    const double R[3][3] = { 
      { cos_dtheta, -sin_dtheta, 0.0 },
      { sin_dtheta, cos_dtheta,  0.0 },
      { 0.0,        0.0,         1.0 }
    };
    for (int ip = 0; ip < n; ++ip) {
      const double xr[3] = {
	x[ip][0] - xc[0],
	x[ip][1] - xc[1],
	x[ip][2] - xc[2]
      };
      FOR_J3 x[ip][j] = DOT_PRODUCT(R[j],xr) + xc[j] + dt*v[j];
    }
  }
  void setCentroid(const double (*x)) {
    FOR_I3 this->xc[i] = x[i];
  }
};

class LinearMotion : public AbstractMotion {
private:
  vector<double> timeVec;
  vector<double> xVec;
  vector<double> yVec;
  vector<double> zVec;
  void lookupPosition(double x[3],const double time) {
    assert(timeVec.size() >= 2);
    if (time <= timeVec[0]) {
      x[0] = xVec[0];
      x[1] = yVec[0];
      x[2] = zVec[0];
    }
    else if (time >= timeVec.back()) {
      x[0] = xVec.back();
      x[1] = yVec.back();
      x[2] = zVec.back();
    }
    else {
      // find the interval that contains time...
      int i1;
      for (i1 = 1; i1 < timeVec.size(); ++i1) {
	if (timeVec[i1] > time)
	  break;
      }
      assert(i1 < timeVec.size());
      const double w0 = (timeVec[i1]-time)/(timeVec[i1]-timeVec[i1-1]);
      assert((w0 >= 0.0)&&(w0 <= 1.0));
      x[0] = w0*xVec[i1-1] + (1.0-w0)*xVec[i1];
      x[1] = w0*yVec[i1-1] + (1.0-w0)*yVec[i1];
      x[2] = w0*zVec[i1-1] + (1.0-w0)*zVec[i1];
    }
  }
public:
  void addEntry(const double time,const double x[3]) {
    timeVec.push_back(time);
    xVec.push_back(x[0]);
    yVec.push_back(x[1]);
    zVec.push_back(x[2]);
  }
  void move(double (*x)[3],const int n,const double time,const double dt) {
    double x0[3];
    lookupPosition(x0,time-dt);
    double dx[3];
    lookupPosition(dx,time);
    FOR_I3 dx[i] -= x0[i];
    for (int ip = 0; ip < n; ++ip) 
      FOR_I3 x[ip][i] += dx[i];
  }

  void setCentroid(const double (*x)) {
    // nothing to do
  }
};
  
// BOUNCE: constant velocity with bouncing 
class BounceProfile : public AbstractMotion {
private:
  double v[3];
  double xc[3];
public:
  BounceProfile(const double v[3]) : AbstractMotion() {
    FOR_I3 this->v[i] = v[i];
    xc[0] = 0.0;
    xc[1] = 0.0;
    xc[2] = 0.0;
  }
  void move(double (*x)[3],const int n,const double time,const double dt) {
    const double xmax = 2.0;
    const double xmin = -1.0;
    const double ymax = 1.0;
    const double ymin = -1.0;
    const double r = 0.25;

    double dx[3];
    
    // Update location
    FOR_I3 dx[i] = v[i]*dt;
    for (int ip = 0; ip < n; ++ip) 
      FOR_I3 x[ip][i] += dx[i];

    // Update centroid (used for velocity update below)
    FOR_I3 xc[i] += dx[i];
    
    // Update velocity
    if ((v[0] > 0.0)&&(xc[0] >= xmax-r)) {
      v[0] = -v[0];
    }
    else if ((v[0] < 0.0)&&(xc[0] < xmin+r)) {
      v[0] = -v[0];
    }
    else if ((v[1] > 0.0)&&(xc[1] > ymax-r)) {
      v[1] = -v[1];
    }
    else if ((v[1] < 0.0)&&(xc[1] < ymin+r)) {
      v[1] = -v[1];
    }

  }
  void setCentroid(const double (*x)) {
    FOR_I3 this->xc[i] = x[i];
  }

};

// Engine profile - valve lifts and piston motion from files
class PistonMotion : public AbstractMotion {
private:
  double radius;
  double conrod;
  double omega; //360deg*rpm/60
  double start_angle;
  double cycle_angle;
  double nx[3];

public:
  PistonMotion(double stroke, double conrod, double rpm, double nx[3]) : AbstractMotion() {
    this->radius = 0.5*stroke;
    this->conrod = conrod;
    this->omega = 360.0*rpm/60.0; //degree
    FOR_I3 this->nx[i] = nx[i];
    // default angle for one full cycle is 720 degree ...
    this->start_angle = 360.0; // default start angle (deg) when time is zero...
    if (checkParam("START_ANGLE")) this->start_angle = getDoubleParam("START_ANGLE");
    this->cycle_angle = 720; 
  }
  
  double pistonPosition(const double rad) {
    double lrod = this->conrod*this->conrod-this->radius*this->radius*sin(rad)*sin(rad);
    assert(lrod >=0.0);
    return(fabs(this->radius*cos(rad)+sqrt(lrod)-(this->radius+this->conrod)));
  }

  void move0(double (*x)[3],const int n,const double time) {
    double crank_angle  = this->omega*time +this->start_angle; // default start is 360 deg
    double dx=pistonPosition(crank_angle*2.0*M_PI/360.0);
    //cout << "piston crank=" << " " << crank_angle << " " << dx << endl;
    for (int ip = 0; ip < n; ++ip) 
      FOR_I3 x[ip][i] += dx*this->nx[i];
  }
  
  void move(double (*x)[3],const int n,const double time,const double dt) {
    double crank_angle0  = this->omega*(time-dt) + this->start_angle; // default start is 360 deg
    double crank_angle  = this->omega*time +this->start_angle; // default start is 360 deg
    double x0=pistonPosition(crank_angle0*2.0*M_PI/360.0);
    double dx=pistonPosition(crank_angle*2.0*M_PI/360.0);
    //cout << "piston crank=" << " " << crank_angle << " " << dx << endl;
    dx -= x0;
    for (int ip = 0; ip < n; ++ip) 
      FOR_I3 x[ip][i] += dx*this->nx[i];
  }
  
  void setCentroid(const double (*x)) {
    // nothing to do
  }
};


// Amazon testing - motion for 2D cylinder (x = Ax sin(wt) and/or y = Ay sin(wt))
class CylinderMotion : public AbstractMotion {
private:
  double Ax; //amplitude in x motion
  double Ay; //amplitude in y motion
  double omega; 

public:
  CylinderMotion(double Ax, double Ay, double omega) : AbstractMotion() {
    this->Ax = Ax;   
    this->Ay = Ay;
    this->omega = omega; 
  }
  
  void move(double (*x)[3],const int n,const double time,const double dt) {
    double dx=this->Ax*(sin(this->omega*time)-sin(this->omega*(time-dt)));
    double dy=this->Ay*(sin(this->omega*time)-sin(this->omega*(time-dt)));
    for (int ip = 0; ip < n; ++ip) {
      x[ip][0] += dx;
      x[ip][1] += dy;
    }
  }
  
  void setCentroid(const double (*x)) {
    // nothing to do
  }
};

class ValveMotion : public AbstractMotion {
private:
  int n;
  double (*data)[2];
  double nx[3];
  bool set_cycle;
  double cycle_angle; // 720 deg by default
  double omega; // deg/sec: 360*rpm/60
  double start_angle; // 360 deg by default

public:
  ValveMotion(string filename, double rpm, double nx[3]) : AbstractMotion() {
    // Read from the file
    FOR_I3 this->nx[i] = nx[i];
    this->omega = 360.0*rpm/60.0;
    // default angle for one full cycle is 720 degree ...
    this->set_cycle = true;
    this->cycle_angle = 720; 
    this->start_angle = 360; //360 deg by default...
    if (checkParam("START_ANGLE")) this->start_angle = getDoubleParam("START_ANGLE");
    readProfile(filename);
  }
  
  void move0(double (*x)[3],const int n,const double time) {
    double crank_angle  = this->omega*time + this->start_angle; // default start is 360 deg
    double lift = getData(crank_angle); // deg
    double dx[3];
    FOR_I3 dx[i] = lift*this->nx[i];
    for (int ip = 0; ip < n; ++ip) 
      FOR_I3 x[ip][i] += dx[i];
  }

  void move(double (*x)[3],const int n,const double time,const double dt) {
    double crank_angle0 = this->omega*(time-dt) + this->start_angle;
    double crank_angle  = this->omega*time + this->start_angle; // default start is 360 deg
    double lift0 = getData(crank_angle0);
    double dlift = getData(crank_angle)-lift0;
    //cout <<"CRANK=" << " " << crank_angle << " " << dlift << endl;
    double dx[3];
    FOR_I3 dx[i] = dlift*this->nx[i];
    for (int ip = 0; ip < n; ++ip) 
      FOR_I3 x[ip][i] += dx[i];
  }

  void setCentroid(const double (*x)) {
    // nothing to do
  }

private:

  void readProfile(const string& filename){

    if ( mpi_rank == 0 )
      cout << "reading 1D valve-lift profile from " << filename << endl;
    // first read the file to get the size
    ifstream fp;
    fp.open(filename.c_str());
    if ( fp.fail()){
      CERR("Error: could not find engine profile: " << filename);
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

    double min_lift = 1.0E10;
    double max_lift = -1.0E10;
    for ( int i=1; i<n; ++i ) {
      min_lift = min(min_lift,data[i][1]);
      max_lift = max(max_lift,data[i][1]);
    }
    if ( mpi_rank == 0 ) {
      cout << " > found    " << n << " points" << " " << "in "<< filename << endl;
      cout << " > crank angle range: (" << data[0][0] << "," << data[n-1][0] << ")" << endl;
      cout << " > lift range: (" <<  min_lift     << "," << max_lift  << ")" << endl;
    }
  }

  double getData(double x) {

    // cycle periodic ?
    if (x < 0.0) { 
      if (this->set_cycle) x=x-floor(x/this->cycle_angle)*this->cycle_angle;
      else return(this->data[0][1]);
    }
    else if (x > this->cycle_angle) {
      if (this->set_cycle) x=x-floor(x/this->cycle_angle)*this->cycle_angle;
      else return(this->data[n-1][1]);
    }

    for ( int i=0; i<n; ++i ){
      if ( x <= this->data[i][0] ){
        if (i>0) {
          assert((this->data[i][0]-this->data[i-1][0])>0.0);
          double w = (this->data[i][0]-x)/(this->data[i][0]-this->data[i-1][0]);
          return( w*this->data[i-1][1]+(1.0-w)*this->data[i][1]);
        }
        else if (i==0) {
          assert(x>=this->data[n-1][0]-this->cycle_angle);
          assert(this->data[0][0]-(this->data[n-1][0]-this->cycle_angle)>0.0);
          double w = (this->data[0][0]-x)/(this->data[0][0]-(this->data[n-1][0]-this->cycle_angle));
          return( w*this->data[n-1][1]+(1.0-w)*this->data[0][1]);
        }
      }
      else if (x > this->data[n-1][0]) {
        assert(this->data[0][0]+this->cycle_angle-this->data[n-1][0]>0.0);
        double w = (this->data[0][0]+this->cycle_angle-x)/(this->data[0][0]+this->cycle_angle-this->data[n-1][0]);
        return( w*this->data[n-1][1]+(1.0-w)*this->data[0][1]);
      }
    }
  
    CERR("Error: should not be here");
  }



};
 
// FLUENT PROFILE MOTION
class FluentProfile : public AbstractMotion {
private:
  string part_name;
  string fluent_profile;
  string motion_dir;
  string pos_calc;
  double **motiondata;
  int ncol;
  string *col_name;
  
  int t_index;
  int p_index[3];
  int ntimes;
  int loop;
  double end_time;
  int first_call;

public:
  FluentProfile(string fluent_profile, string part_name, string motion_dir) : AbstractMotion() {
    
    assert(0); // use move0 dummy

    // Read in the profile and setup data
    this->part_name = part_name;
    this->motion_dir = motion_dir;
    this->fluent_profile = fluent_profile;
    t_index = -1;
    FOR_I3 p_index[i] = -1;
    pos_calc = "LINEAR";
    end_time = 0.0;
    loop = -1;
    first_call = -1;
  }

  // For Fluent profile, the displacement from the profile file is computed by comparing the positions at the current time, and the previous time
  // The previous time is computed as time-dt, rather than storing the last time the routine is called. This should work even for the first timestep
  // and negative times are ok in the interpolation routine below

  void setCentroid(const double (*x)) {
    // nothing to do for now
  }

  void move(double (*x)[3],const int n,const double time,const double dt) {
    // Read in the profile file if not done yet
    if (first_call == -1) {
      readFluentProfile(fluent_profile);
      first_call = 0;
    }
    
    double localtime, localtime_p;
    localtime = time;
    localtime_p = time - dt;
    if (loop == 1) {
      localtime = fmod(localtime,end_time);
      localtime_p = fmod(localtime_p,end_time);
    }

    // Do nothing if both the previous and new times are outside profile. This already accounts for loop here
    // due to the time correction above
    if (localtime_p >= motiondata[t_index][ntimes-1]) {
      assert(localtime > motiondata[t_index][ntimes-1]); // monotonic time
      return;
    }

    if (localtime <= motiondata[t_index][0]) {
      assert(localtime_p < motiondata[t_index][0]); // monotonic time
      return;
    }

    // At least one of the times should be inside the profile now
    // Find an index for previous and new times
    int index_p, index;
    index_p = index = -1;

    // Get locations at these times now
    double x_n[3]; FOR_I3 x_n[i] = 0.0;
    double x_p[3]; FOR_I3 x_p[i] = 0.0;

    // If old time is before the profile, assume the first profile location as start
    if (localtime_p <= motiondata[t_index][0]) {
      FOR_I3 {
	if (p_index[i] != -1)
	  x_p[i] = motiondata[p_index[i]][0];
      }
    } 
    else {
      // Update to the time index less than or equal to the localtime_p.
      for (int i=ntimes-1; i>=0; i--) {
	if (motiondata[t_index][i] <= localtime_p) {
	  index_p = i;
	  break;
	}
      }
    }

    // If current time position is after the profile ending, assume the first profile location as startCopy over last value into current if new time is at end of profile
    if (localtime >= motiondata[t_index][ntimes-1]) {
      FOR_I3 {
	if (p_index[i] != -1)
	  x_n[i] = motiondata[p_index[i]][ntimes-1];
      }
    }
    else {
      // Update to the time index less than or equal to the localtime_p.
      for (int i=ntimes-1; i>=0; i--) {
	if (motiondata[t_index][i] <= localtime) {
	  index = i;
	  break;
	}
      }
    }

    assert(index_p != ntimes-1); // double check that the previous time is not last
    //int loop_interp;

    // Set locations for previous/new time
    FOR_I3 {
      if (p_index[i] != -1) {
	if (index_p != -1) {
	  // Simple linear interpolation for now
	  if (pos_calc.compare("LINEAR")==0)
	    linearInterpolate(x_p[i], localtime_p, p_index[i], index_p);
	}

	if (index != -1) {
	  // Simple linear interpolation for now
	  if (pos_calc.compare("LINEAR")==0)
	    linearInterpolate(x_n[i], localtime, p_index[i], index);
	}
      }
    }

    // Increment x by the dx
    double dx[3];
    FOR_I3 dx[i] = x_n[i] - x_p[i];
    
    for (int ip = 0; ip < n; ++ip) 
      FOR_I3 x[ip][i] += dx[i];

    // COUT1("x at previous time for part " << part_name << " = " << COUT_VEC(x_p) << " at time " << localtime_p);
    // COUT1("x at next time for part " << part_name << " = " << COUT_VEC(x_n) << " at time " << localtime);
    COUT1("DX for part " << part_name << " = " << COUT_VEC(dx) << " at time " << time);

    return;
  }

private:
  
  void linearInterpolate(double &value, double localtime, int v_index, int index) {
    if (index == ntimes - 1) {
      value =  motiondata[t_index][ntimes-1]; // assume we are at the end 
      return;
    }
    double last_time = motiondata[t_index][index];
    double next_time = motiondata[t_index][index+1];
    value = (localtime-last_time)*(motiondata[v_index][index+1]-motiondata[v_index][index]);
    value = value/(next_time-last_time)+motiondata[v_index][index];
  }

  // Read Fluent Motion Profile in this routine
  void readFluentProfile(string filename) {
    
    // File variable
    ifstream fh;
    
    const int max_chars_line = 100;
    const int max_tokens_line = 40;
    const char DELIMITER[] = " \t";
    string line;
    
    fh.open(filename.c_str());    
    if (!fh.good()) {
      COUT1("\n\nError: could not open file: " << filename << "\n\n");
      throw(-1);
    }
    else
      {
        COUT1("Reading Fluent Profile file " << filename << " for object " << part_name);
      }
    
    // read top two lines
    int i;
    for (i=0; i<2; i++)
      {
        char buf[max_chars_line];
        fh.getline(buf, max_chars_line);
        int n = 0;
        const char* token[max_tokens_line] = {};
        token[0]= strtok(buf, DELIMITER);
        if (token[0])
          {
            for (n=1; n<max_tokens_line; n++)
              {
                token[n] = strtok(0,DELIMITER);
                if (!token[n]) break;
              }
          }
    
        if (i==0) {
          ncol = atoi(token[1]);
          ntimes = atoi(token[2]);
          loop = atoi(token[3]);
          col_name = new string[ncol];
          COUT1("N col, N rows, Loop : " << ncol << ", " << ntimes << ", " << loop);
        } 
        else if (i==1)
          {
            if (n!=ncol)
              {
                cerr << "Incorrect number of cell array names. Check";
                throw(-1);
              }
            for (int j=0; j<ncol; j++)
              {
                col_name[j] = token[j];
              }
          }        
      }

    // Allocate arrays
    motiondata = new double*[ncol];
    for (i=0; i<ncol; i++)
      motiondata[i] = new double[ntimes];

    i = 0;
    while(!fh.eof())
      {
        char buf[max_chars_line];
        fh.getline(buf, max_chars_line);
        int n = 0;
        const char* token[max_tokens_line] = {};
        token[0]= strtok(buf, DELIMITER);
        int tcount = 0;
        if (token[0])
          {
            tcount++;
            for (n=1; n<max_tokens_line; n++)
              {
                token[n] = strtok(0,DELIMITER);
                if (!token[n]) break;
                tcount++;
              }
          }
        if (tcount==ncol && i<ntimes)
          {
            for (int j=0; j<ncol; j++)
              {
                motiondata[j][i] = atof(token[j]);
              }
          }
        i++;
      } 

    if (i-1<ntimes)
      {
	cerr << "Incorrect number of time rows. Expected : " << ntimes << ", Actual : " << i-1 << ". Check" << endl;
	throw(-1);
      }
     
    // Find relevant strings in col name
    // Assume time is the first one
    COUT1("Time : " << col_name[0]);
    dumpRange(&motiondata[0][0], ntimes, "Time Min/Max");
    t_index = 0;
    end_time = motiondata[0][ntimes-1];

    if (motion_dir.compare("X")==0) {
      assert(ncol>1); p_index[0] = 1; }
    else if (motion_dir.compare("Y")==0) {
      assert(ncol>1); p_index[1] = 1; }
    else if (motion_dir.compare("Z")==0) {
      assert(ncol>1); p_index[2] = 1; }
    else if (motion_dir.compare("XY")==0) {
      assert(ncol>2); p_index[0] = 1; p_index[1] = 2; }
    else if (motion_dir.compare("XZ")==0) {
      assert(ncol>2); p_index[0] = 1; p_index[2] = 2; }
    else if (motion_dir.compare("YZ")==0) {
      assert(ncol>2); p_index[1] = 1; p_index[2] = 1; }
    else if (motion_dir.compare("XYZ")==0) {
      assert(ncol>3); p_index[0] = 1; p_index[1] = 2;  p_index[2] = 3; }
    else {
      cerr <<  "Motion type " << motion_dir << " not implemented yet" << endl;
      throw(-1); 
    }
    if (p_index[0] != -1) {
      COUT1("X Displacement  : " << col_name[p_index[0]]);
      dumpRange(&motiondata[p_index[0]][0], ntimes, "X Displacement Min/Max");
    }
    if (p_index[1] != -1) {
      COUT1("Y Displacement  : " << col_name[p_index[1]]);
      dumpRange(&motiondata[p_index[1]][0], ntimes, "Y Displacement Min/Max");
    }
    if (p_index[2] != -1) {
      COUT1("Z Displacement  : " << col_name[p_index[2]]);
      dumpRange(&motiondata[p_index[2]][0], ntimes, "Z Displacement Min/Max");
    }
     
    fh.close();
     
    if (t_index==-1) {
      cerr << "\n Time index not defined . Please check.\n\n" << endl;
      throw(-1);
    }
    
    bool found = false;
    FOR_I3 {
      if (p_index[i]!=1) 
        found = true;
    }
    
    if (!found) {
      cerr << "\n No position X or Y or Z found here . Please check.\n\n" << endl;
      throw(-1);
    }
  }
  
};


/*
class StarfishMotion : public AbstractMotion {
private:
  int n;
public:
  StarfishMotion(const int n) {
    this->n = n;
  }
  void move(double (*x)[3],const int n,const double time,const double dt) {
    for (int ip = 0; ip < n; ++ip) {
      const double theta = atan2(x[ip][1],x[ip][0]);
      const double r = 
      

      FOR_J3 x[i][j] += dt*v[j];
    }
  }
};
*/
    
    
#endif
