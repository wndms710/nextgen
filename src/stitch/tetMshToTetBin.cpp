#include "CTI.hpp"
using namespace CTI;

#include "FluentReader.hpp"
#include "Adt.hpp"

class FluentMshPlus : public FluentMsh {
public:

  int zone;
  int *no_flag;
  int nte;
  int (*noote)[4];
  
  FluentMshPlus() {
    no_flag = NULL; 
    nte = 0;
    noote = NULL;
  }
  
  FluentMshPlus(const string& filename) : FluentMsh(filename) {
    no_flag = NULL; 
    nte = 0;
    noote = NULL;
    buildTets(nte,noote);
    cout << " > got nte: " << nte << endl;
  }

  virtual ~FluentMshPlus() {
    DELETE(noote);
    DELETE(no_flag);
  }

  bool checkVolume() {
    double vol_sum = 0.0;
    double vol_min = HUGE_VAL;
    double vol_max = -HUGE_VAL;
    for (int ite = 0; ite < nte; ++ite) {
      // flip first 2 nodes?...
      //const int tmp = noote[ite][0];
      //noote[ite][0] = noote[ite][1];
      //noote[ite][1] = tmp;
      // and check volume...
      const int ino0 = noote[ite][0];
      const int ino1 = noote[ite][1];
      const int ino2 = noote[ite][2];
      const int ino3 = noote[ite][3];
      const double vol = SIGNED_TET_VOLUME_6(x_no[ino0],x_no[ino1],x_no[ino2],x_no[ino3])/6.0;
      vol_sum += vol;
      vol_min = min(vol_min,vol);
      vol_max = max(vol_max,vol);
    }
    // we require positive volumes...
    cout << " > total volume: " << vol_sum << ", individual tet volume (min,avg,max): " << 
      vol_min << " " << vol_sum/double(nte) << " " << vol_max << endl;
    return (vol_sum <= 0.0);
  }

  void flipVolume() {
    for (int ite = 0; ite < nte; ++ite) {
      // flip first 2 nodes...
      const int tmp = noote[ite][0];
      noote[ite][0] = noote[ite][1];
      noote[ite][1] = tmp;  
    }
  }
};

void processMsh() {
  
  // cap transformations, surface and tet writing...
  // Feb 2020
  
  vector<FluentMshPlus*> mshVec;
  FOR_PARAM_MATCHING("MSH") {

    string msh_filename = param->getString();
    
    cout << "\n===================================================================\n" << 
      " > working on msh file: " << msh_filename << "..." << endl;
    
    FluentMshPlus * msh = new FluentMshPlus(msh_filename);
    
    msh->zone = 0;
    
    // parse remaining args...
    int iarg = 1;
    while (iarg < param->size()) {
      string token = MiscUtils::toUpperCase(param->getString(iarg++));
      if (token == "ZONE") {
	msh->zone = param->getInt(iarg++);
	cout << " > ZONE set to " << msh->zone << endl;
      }
      else if (token == "SCALE") {
	const double factor = param->getDouble(iarg++);
	cout << " > applying SCALE " << factor << endl;
	for (int ino = 0; ino < msh->nno; ++ino) {
	  FOR_I3 msh->x_no[ino][i] *= factor;
	}	
      }
      else if (token == "TRANSLATE") {
	double dx[3];
	FOR_I3 dx[i] = param->getDouble(iarg++);
	cout << " > applying TRANSLATE " << dx[0] << " " << dx[1] << " " << dx[2] << endl;
	for (int ino = 0; ino < msh->nno; ++ino) {
	  FOR_I3 msh->x_no[ino][i] += dx[i];
	}	
      }
      else {
	CERR("unrecognized MSH token: " << token);
      }
    }

    // write the full surface mesh out. Delete zones in surfer to get the mesh to 
    // replace the old cap with... 
    msh->writeSurfaceSbin(msh_filename+".sbin");

    double bbmin[3] = { HUGE_VAL, HUGE_VAL, HUGE_VAL };
    double bbmax[3] = { -HUGE_VAL, -HUGE_VAL, -HUGE_VAL };
    for (int ino = 0; ino < msh->nno; ++ino) {
      FOR_I3 bbmin[i] = min(bbmin[i],msh->x_no[ino][i]);
      FOR_I3 bbmax[i] = max(bbmax[i],msh->x_no[ino][i]);
    }
    cout << " > bbox: " << COUT_VEC(bbmin) << " " << COUT_VEC(bbmax) << endl;

    bool is_neg = msh->checkVolume();
   
    if (is_neg) {
      COUT1(" ! negative volume was detected; flipping orientation of tets to get positive volume");
      // surface is flipped for some reason; reorder noote?
      msh->flipVolume();
      is_neg = msh->checkVolume();
      assert(is_neg == false);
    }

    // and add to the mshVec...
    mshVec.push_back(msh);
    
  }
  
  FOR_PARAM_MATCHING("TET_BIN") {

    string tet_bin_filename = param->getString();
    
    cout << "\n===================================================================\n" << 
      " > working on tet_bin file: " << tet_bin_filename << "..." << endl;
    
    // hack: for now, build an empty 
    FluentMshPlus * msh = new FluentMshPlus();
    msh->zone = 0;
    
    // read the tet bin format: it is trivial...
    bool b_zones = false;
    FILE *fp = fopen(tet_bin_filename.c_str(),"rb");
    assert(fp);
    fread(&msh->nno,sizeof(int),1,fp);
    if (msh->nno == -1) {
      // a leading -1 is used to indicate a tet_bin file with zone info. Note:
      // can't let this go on for too long. Need to make a proper file that 
      // is based on records and can journal where it came from...
      b_zones = true;
      fread(&msh->nno,sizeof(int),1,fp);
    }
    fread(&msh->nte,sizeof(int),1,fp);
    cout << "   > nno: " << msh->nno << " nte: " << msh->nte << " b_zones: " << b_zones << endl;
    assert(msh->x_no == NULL);
    msh->x_no = new double[msh->nno][3];
    fread(msh->x_no,sizeof(double),msh->nno*3,fp);
    assert(msh->noote == NULL);
    msh->noote = new int[msh->nte][4];
    fread(msh->noote,sizeof(int),msh->nte*4,fp);
    // skip the zone data at the end for now...
    assert(!b_zones);
    fclose(fp);

    msh->zone = 0;
    
    // parse remaining args...
    int iarg = 1;
    while (iarg < param->size()) {
      string token = MiscUtils::toUpperCase(param->getString(iarg++));
      if (token == "ZONE") {
	msh->zone = param->getInt(iarg++);
	cout << " > ZONE set to " << msh->zone << endl;
      }
      else if (token == "SCALE") {
	const double factor = param->getDouble(iarg++);
	cout << " > applying SCALE " << factor << endl;
	for (int ino = 0; ino < msh->nno; ++ino) {
	  FOR_I3 msh->x_no[ino][i] *= factor;
	}	
      }
      else if (token == "TRANSLATE") {
	double dx[3];
	FOR_I3 dx[i] = param->getDouble(iarg++);
	cout << " > applying TRANSLATE " << dx[0] << " " << dx[1] << " " << dx[2] << endl;
	for (int ino = 0; ino < msh->nno; ++ino) {
	  FOR_I3 msh->x_no[ino][i] += dx[i];
	}	
      }
      else {
	CERR("unrecognized MSH token: " << token);
      }
    }
    
    double bbmin[3] = { HUGE_VAL, HUGE_VAL, HUGE_VAL };
    double bbmax[3] = { -HUGE_VAL, -HUGE_VAL, -HUGE_VAL };
    for (int ino = 0; ino < msh->nno; ++ino) {
      FOR_I3 bbmin[i] = min(bbmin[i],msh->x_no[ino][i]);
      FOR_I3 bbmax[i] = max(bbmax[i],msh->x_no[ino][i]);
    }
    cout << " > bbox: " << COUT_VEC(bbmin) << " " << COUT_VEC(bbmax) << endl;
    
    msh->checkVolume();
    
    // and add to the mshVec...
    mshVec.push_back(msh);
    
  }
  
  cout << "\n===================================================================" << endl; 
  
  // now condense the nodes. First we need a length scale...
  
  int nno_max = 0;
  int nte = 0;
  double eps = HUGE_VAL;
  for (int imsh = 0, size = mshVec.size(); imsh < size; ++imsh) {
    FluentMshPlus* msh =  mshVec[imsh];
    nno_max += msh->nno;
    nte += msh->nte;
    for (int ite = 0; ite < msh->nte; ++ite) {
      for (int i0 = 0; i0 < 3; ++i0) {
	const int ino0 = msh->noote[ite][i0];
	for (int i1 = i0+1; i1 < 4; ++i1) {
	  const int ino1 = msh->noote[ite][i1];
	  const double d2 = DIST2(msh->x_no[ino0],msh->x_no[ino1]);
	  eps = min(eps,d2);
	}
      }
    }
  }
  eps = sqrt(eps);
  
  cout << " > min length scale: " << eps << endl;

  assert(mshVec.size() >= 1);
  double (*x_no)[3] = new double[nno_max][3];
  int nno = mshVec[0]->nno;
  for (int ino = 0; ino < nno; ++ino) {
    FOR_I3 x_no[ino][i] = mshVec[0]->x_no[ino][i];
  }
  for (int imsh = 1, size = mshVec.size(); imsh < size; ++imsh) {
    cout << " > condensing nodes from msh " << imsh << "..." << endl;
    FluentMshPlus* msh =  mshVec[imsh];
    assert(msh->no_flag == NULL);
    msh->no_flag = new int[msh->nno];
    // build an adt with the current x_no, and search for these xno's...
    Adt<double> * adt = new Adt<double>(nno,x_no,x_no);
    cout << " > adt built" << endl;
    vector<int> intVec;
    double d2_max = 0.0;
    int count = 0;
    for (int ino = 0; ino < msh->nno; ++ino) {
      //if (ino%100 == 0) cout << " > node " << ino << " out of " << msh->nno << endl;
      assert(intVec.empty());
      adt->buildListForSphere(intVec,msh->x_no[ino],0.1*eps);
      if (!intVec.empty()) {
	assert(intVec.size() == 1);
	const int ino_new = intVec[0];
	msh->no_flag[ino] = ino_new;
	const double d2 = DIST2(x_no[ino_new],msh->x_no[ino]);
	d2_max = max(d2_max,d2);
	++count;
	intVec.clear();
      }
      else {
	const int ino_new = nno++;
	msh->no_flag[ino] = ino_new;
	FOR_I3 x_no[ino_new][i] = msh->x_no[ino][i];
      }
    }
    delete adt;
    cout << " > matched " << count << " nodes with max dist: " << sqrt(d2_max) << endl; 
    // change the node references in noote...
    for (int ite = 0; ite < msh->nte; ++ite) {
      FOR_I4 {
	const int ino = msh->noote[ite][i];
	msh->noote[ite][i] = msh->no_flag[ino];
      }
    }
  }

  cout << " > final nno: " << nno << " compare to nno_max: " << nno_max << endl;

  // write the tets to a simple binary file...
  cout << " > writing cht.tet_bin with nno: " << nno << ", nte: " << nte << endl;
  FILE * fp = fopen("cht.tet_bin","wb");
  assert(fp != NULL);
  int minus_one = -1;
  fwrite(&minus_one,sizeof(int),1,fp); // use to indicate the version that includes the tet zone 
  fwrite(&nno,sizeof(int),1,fp);
  fwrite(&nte,sizeof(int),1,fp);
  fwrite(x_no,sizeof(double),nno*3,fp);
  
  // tet stuff...
  for (int imsh = 0, size = mshVec.size(); imsh < size; ++imsh) {
    FluentMshPlus* msh =  mshVec[imsh];
    fwrite(msh->noote,sizeof(int),msh->nte*4,fp);
  }
  
  // then the zone...
  for (int imsh = 0, size = mshVec.size(); imsh < size; ++imsh) {
    FluentMshPlus* msh =  mshVec[imsh];
    int * znote = new int[msh->nte];
    for (int ite = 0; ite < msh->nte; ++ite) {
      znote[ite] = msh->zone;
    }
    fwrite(znote,sizeof(int),msh->nte,fp);
    delete[] znote;
  }
  fclose(fp);
  
  cout << "done" << endl;
  
}

int main(int argc, char * argv[]) {

  try {
    
    CTI_Init(argc,argv,"tetMshToTetBin.in");

    {
      
      processMsh();

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

