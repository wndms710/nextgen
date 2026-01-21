#include "LesPost.hpp"

class MyLesPost : public LesPost {
public:

  int ncv;
  double * vol_cv;
  double (*x_cv)[3];
  double * p;
  double * q;
  double * e_src;
  
  double * p_sum;
  double * q_sum;
  double * pq_sum;
  double * e_src_sum;
  double * pe_src_sum;
  int count;
  
  MyLesPost() {
    
    vol_cv = NULL;
    x_cv = NULL;

    p = NULL;
    q = NULL;
    e_src = NULL;

    p_sum = NULL;
    q_sum = NULL;
    pq_sum = NULL;
    e_src_sum = NULL;
    pe_src_sum = NULL;
    
  }

  ~MyLesPost() {
    
    DELETE(vol_cv);
    DELETE(x_cv);
    
    DELETE(p);
    DELETE(q);
    DELETE(e_src);

    DELETE(p_sum);
    DELETE(q_sum);
    DELETE(pq_sum);
    DELETE(e_src_sum);
    DELETE(pe_src_sum);
    
  }
  
  void initialHook() {

    cout << "initialHook()..." << endl;
    
    string mles = getStringParam("RESTART");
    
    // read the volume...
    cout << " > reading vol_cv..." << endl;
    SerialLesReader::readCvDN(vol_cv,ncv,"vol_cv",mles);
    assert(vol_cv != NULL);
    cout << " > got ncv: " << ncv << endl;
    
    // read x_cv...
    cout << " > reading x_cv..." << endl;
    x_cv = new double[ncv][3];
    SerialLesReader::readCvDN3(x_cv,ncv,"x_cv",mles);
    assert(x_cv != NULL);

    // allocate arrays used in temporalHook... 
    p = new double[ncv];
    q = new double[ncv];
    e_src = new double[ncv];

    p_sum = new double[ncv];
    q_sum = new double[ncv];
    pq_sum = new double[ncv];
    e_src_sum = new double[ncv];  
    pe_src_sum = new double[ncv];
    
    FOR_ICV {
      p_sum[icv] = 0.0;
      q_sum[icv] = 0.0;
      pq_sum[icv] = 0.0;
      e_src_sum[icv] = 0.0;
      pe_src_sum[icv] = 0.0;
    }
    count = 0;
    
  }
  
  void temporalHook(const string& snapshot,const int step) {
    
    cout << "temporalHook(): snapshot \"" << snapshot << "\" step " << step << "..." << endl;
    
    SerialLesReader::readCvDN(p,ncv,"p",snapshot);
    SerialLesReader::readCvDN(q,ncv,"heat_release()",snapshot);
    SerialLesReader::readCvDN(e_src,ncv,"e_src",snapshot);
    
    // output the pressure integrated of the first inch of the can...
    double p_vol_sum = 0.0;
    double q_vol_sum = 0.0;
    double e_src_vol_sum = 0.0;
    double vol_sum = 0.0;
    double total_vol_sum = 0.0;
    FOR_ICV {
      if ((x_cv[icv][2] >= 0)&&(x_cv[icv][2] < 0.0254)) {
	p_vol_sum += p[icv]*vol_cv[icv];
	vol_sum += vol_cv[icv];
      }
      // and q and e_src everywhere...
      q_vol_sum += q[icv]*vol_cv[icv];
      e_src_vol_sum += e_src[icv]*vol_cv[icv];
      total_vol_sum += vol_cv[icv];
    }
    assert(vol_sum > 0.0);
    cout << setprecision(10) << "PCAP,Q,ESRC " << step << " " << p_vol_sum/vol_sum << " " << q_vol_sum/total_vol_sum << " " << e_src_vol_sum/total_vol_sum << endl;

    FOR_ICV {
      p_sum[icv]       += p[icv];
      q_sum[icv]       += q[icv];
      pq_sum[icv]      += p[icv]*q[icv];
      e_src_sum[icv]   += e_src[icv];
      pe_src_sum[icv]  += p[icv]*e_src[icv];
    }
    ++count;
    
    /*
    const double xmax = 4.7;
    double vol_sum[2] = { 0.0, 0.0 };
    double vol_p_sum[2] = { 0.0, 0.0 };
    double vol_q_sum[2] = { 0.0, 0.0 };
    double vol_pq_sum[2] = { 0.0, 0.0 };
    FOR_ICV {
      if (x_cv[icv][0] < xmax) {
	const double theta = atan2(x_cv[icv][2],x_cv[icv][1]);
	if (theta*180.0/M_PI < 101.25) {
	  // can 1...
	  vol_sum[0] += vol_cv[icv];
	  vol_p_sum[0] += vol_cv[icv]*p[icv];
	  vol_q_sum[0] += vol_cv[icv]*q[icv];
	  vol_pq_sum[0] += vol_cv[icv]*p[icv]*q[icv];
	}
	else {
	  vol_sum[1] += vol_cv[icv];
	  vol_p_sum[1] += vol_cv[icv]*p[icv];
	  vol_q_sum[1] += vol_cv[icv]*q[icv];
	  vol_pq_sum[1] += vol_cv[icv]*p[icv]*q[icv];
	}
      }
    }
    
    cout << setprecision(10) << "CAN1 " << vol_sum[0] << " " << vol_p_sum[0]/vol_sum[0] << " " << vol_q_sum[0] << " " << vol_pq_sum[0] << endl;
    cout << setprecision(10) << "CAN2 " << vol_sum[1] << " " << vol_p_sum[1]/vol_sum[1] << " " << vol_q_sum[1] << " " << vol_pq_sum[1] << endl;
    */
    
  }

  void finalHook() {

    cout << "finalHook()..." << endl;
    
    double rg1 = 0.0;
    double rg2 = 0.0;

    FOR_ICV {
      const double p_avg = p_sum[icv]/double(count);
      const double q_avg = q_sum[icv]/double(count);
      const double pq = p_avg*q_avg*double(count) + pq_sum[icv] - p_avg*q_sum[icv] - q_avg*p_sum[icv];
      rg1 += vol_cv[icv]*pq;
      const double e_src_avg = e_src_sum[icv]/double(count);
      const double pe_src = p_avg*e_src_avg*double(count) + pe_src_sum[icv] - p_avg*e_src_sum[icv] - e_src_avg*p_sum[icv];
      rg2 += vol_cv[icv]*pe_src;
    }

    cout << setprecision(10) << "RG: " << rg1 << " " << rg2 << endl;
     
  }
  
};

class MyLesPost2 : public LesPost {
public:

  int ncv;
  double * vol_cv;
  double (*x_cv)[3];
  double * p;
  double (*u)[3];
  
  MyLesPost2() {
    
    vol_cv = NULL;
    x_cv = NULL;
    
    p = NULL;
    u = NULL;
    
  }
  
  ~MyLesPost2() {
    
    DELETE(vol_cv);
    DELETE(x_cv);
    
    DELETE(p);
    DELETE(u);
    
  }
  
  void initialHook() {

    cout << "initialHook()..." << endl;
    
    string mles = getStringParam("RESTART");
    
    // read the volume...
    cout << " > reading vol_cv..." << endl;
    SerialLesReader::readCvDN(vol_cv,ncv,"vol_cv",mles);
    assert(vol_cv != NULL);
    cout << " > got ncv: " << ncv << endl;
    
    // read x_cv...
    cout << " > reading x_cv..." << endl;
    x_cv = new double[ncv][3];
    SerialLesReader::readCvDN3(x_cv,ncv,"x_cv",mles);
    assert(x_cv != NULL);

    // allocate arrays used in temporalHook... 
    p = new double[ncv];
    u = new double[ncv][3];
    
  }
  
  void temporalHook(const string& snapshot,const int step) {
    
    cout << "temporalHook(): snapshot \"" << snapshot << "\" step " << step << "..." << endl;
    
    SerialLesReader::readCvDN(p,ncv,"p",snapshot);
    SerialLesReader::readCvDN3(u,ncv,"u",snapshot);
    
    /*
    // plug z locations: 32" cases...
    const double zp0 = 0.815;
    const double zp1 = 0.825;
    const double zp20 = 0.924;
    const double zp21 = 0.952;
    */

    // plug z locations: 37" cases...
    const double zp0 = 0.942;
    const double zp1 = 0.952;
    const double zp20 = 1.051;
    const double zp21 = 1.079;

    double p_plug_sum = 0.0;
    double uz_plug_sum = 0.0;
    double plug_sum = 0.0; 

    double p_plug2_sum = 0.0;
    double uz_plug2_sum = 0.0;
    double plug2_sum = 0.0;

    double p_vol_sum = 0.0;
    double vol_sum = 0.0;
    
    FOR_ICV {
      if ((x_cv[icv][2] >= zp0)&&(x_cv[icv][2] <= zp1)) {
	p_plug_sum += p[icv]*vol_cv[icv];
	uz_plug_sum += u[icv][2]*vol_cv[icv];
	plug_sum += vol_cv[icv];
      }
      if ((x_cv[icv][2] >= zp20)&&(x_cv[icv][2] <= zp21)) {
	p_plug2_sum += p[icv]*vol_cv[icv];
	uz_plug2_sum += u[icv][2]*vol_cv[icv];
	plug2_sum += vol_cv[icv];
      }
      
      // also cap pressure
      if ((x_cv[icv][2] >= 0)&&(x_cv[icv][2] < 0.0254)) {
	p_vol_sum += p[icv]*vol_cv[icv];
	vol_sum += vol_cv[icv];
      }
      
    }
    
    cout << setprecision(10) << "PUPUPCAP " << step << " " << p_plug_sum/plug_sum << " " << uz_plug_sum/plug_sum << " " << p_plug2_sum/plug2_sum << " " << uz_plug2_sum/plug2_sum << " " << p_vol_sum/vol_sum << endl;

  }

  void finalHook() {

    cout << "finalHook()..." << endl;

  }
  
};
  
int main(int argc, char * argv[]) {

  try {

    CTI_Init(argc,argv,"lespost.in");

    {
      
      //MyLesPost lespost;
      MyLesPost2 lespost;
      lespost.init(); // calls initialHook
      lespost.run(); // does SNAPSHOT loop if any with temporalHook, then finalHook
    
    }

    CTI_Finalize();

  }
  catch (int e) {
    if (e == 0)
      CTI_Finalize();
    else
      CTI_Abort();
  }
  catch(...) {
    CTI_Abort();
  }
  
  return(0);

}
