void CuttableVoronoiData::doit2(const double Lxm_,const double Lxp_,const double Lym_,const double Lyp_,const double Lzm_,const double Lzp_,
                               const int * const ist_and_ipart_buf,const bool debug) {

  const double L = max(max(-Lxm_,Lxp_),max(max(-Lym_,Lyp_),max(-Lzm_,Lzp_)));
  assert(L == Lxp_);
  assert(L == Lyp_);
  assert(L == Lzp_);
  assert(L == -Lxm_);
  assert(L == -Lym_);
  assert(L == -Lzm_);
  
  //static int8 debug_count = 0;
  //const int debug_count_check = getIntParam("debug_count_check",-1); // 223,224,225,226
  
  if (debug) {
    writeTecplot(-1);
    cout << "CuttableVoronoiData::doit: nno_orig: " << nno << " ned_orig: " << ned << endl;
    cout << "original surface in tecplot -1" << endl;
  }
  
  // entering this routine with cube-cut surface patches...
  
  const int nno_orig = nno;
  const int ned_orig = ned;
  int nno_max = nno+8;
  int * no_flag = new int[nno_max]; // add space for up to 8 potential corner nodes (this allows reuse of no_flag)... 

  {

    FOR_INO no_flag[ino] = 0;
    FOR_IED {
      assert(faoed[ied][1] >= -1); 
      //if ((faoed[ied][0] >= -7)&&(faoed[ied][0] <= -2)) {
      if (faoed[ied][0] <= -2) {
	assert(faoed[ied][0] >= -7);
	// this is an edge on one of the cube faces. Set the corresponding bit in both its nodes:
	no_flag[nooed[ied][0]] |= (1<<(-faoed[ied][0]-2+6));
	no_flag[nooed[ied][1]] |= (1<<(-faoed[ied][0]-2));
        // to avoid problems downstream, use this as a chance to force the
        // cube limits to be exactly L...
        switch(faoed[ied][0]) {
        case X0_FACE:
          //cout << "err: " << fabs(x_no[nooed[ied][0]][0]-Lxm)/Lxm << endl;
          //cout << "err: " << fabs(x_no[nooed[ied][1]][0]-Lxm)/Lxm << endl;
          assert(fabs(x_no[nooed[ied][0]][0]+L) < 1.0E-12*L);
          //x_no[nooed[ied][0]][0] = Lxm;
          assert(fabs(x_no[nooed[ied][1]][0]+L) < 1.0E-12*L);
          //x_no[nooed[ied][1]][0] = Lxm;
          break;
        case X1_FACE:
          //cout << "err: " << fabs(x_no[nooed[ied][0]][0]-Lxp)/Lxp << endl;
          //cout << "err: " << fabs(x_no[nooed[ied][1]][0]-Lxp)/Lxp << endl;
          assert(fabs(x_no[nooed[ied][0]][0]-L) < 1.0E-12*L);
          //x_no[nooed[ied][0]][0] = Lxp;
          assert(fabs(x_no[nooed[ied][1]][0]-L) < 1.0E-12*L);
          //x_no[nooed[ied][1]][0] = Lxp;
          break;
        case Y0_FACE:
          //cout << "err: " << fabs(x_no[nooed[ied][0]][1]-Lym)/Lym << endl;
          //cout << "err: " << fabs(x_no[nooed[ied][1]][1]-Lym)/Lym << endl;
          assert(fabs(x_no[nooed[ied][0]][1]+L) < 1.0E-12*L);
          //x_no[nooed[ied][0]][1] = Lym;
          assert(fabs(x_no[nooed[ied][1]][1]+L) < 1.0E-12*L);
          //x_no[nooed[ied][1]][1] = Lym;
          break;
        case Y1_FACE:
          //cout << "err: " << fabs(x_no[nooed[ied][0]][1]-Lyp)/Lyp << endl;
          //cout << "err: " << fabs(x_no[nooed[ied][1]][1]-Lyp)/Lyp << endl;
          assert(fabs(x_no[nooed[ied][0]][1]-L) < 1.0E-12*L);
          //x_no[nooed[ied][0]][1] = Lyp;
          assert(fabs(x_no[nooed[ied][1]][1]-L) < 1.0E-12*L);
          //x_no[nooed[ied][1]][1] = Lyp;
          break;
        case Z0_FACE:
          //cout << "err: " << fabs(x_no[nooed[ied][0]][2]-Lzm)/Lzm << endl;
          //cout << "err: " << fabs(x_no[nooed[ied][1]][2]-Lzm)/Lzm << endl;
          assert(fabs(x_no[nooed[ied][0]][2]+L) < 1.0E-12*L);
          //x_no[nooed[ied][0]][2] = Lzm;
          assert(fabs(x_no[nooed[ied][1]][2]+L) < 1.0E-12*L);
          //x_no[nooed[ied][1]][2] = Lzm;
          break;
        case Z1_FACE:
          //cout << "err: " << fabs(x_no[nooed[ied][0]][2]-Lzp)/Lzp << endl;
          //cout << "err: " << fabs(x_no[nooed[ied][1]][2]-Lzp)/Lzp << endl;
          assert(fabs(x_no[nooed[ied][0]][2]-L) < 1.0E-12*L);
          //x_no[nooed[ied][0]][2] = Lzp;
          assert(fabs(x_no[nooed[ied][1]][2]-L) < 1.0E-12*L);
          //x_no[nooed[ied][1]][2] = Lzp;
          break;
        default:
          assert(0);
        }
      }
    }

    //getchar();
    
    // store surface-cube edge intersections in a vec, one for each of the 
    // 12 cube edges. We use a vec here so we can sort along the edge if there
    // are multiple intersections...
    // edge 0 faces x0,y0 nodes 0,4
    // edge 1 faces y1,x0 nodes 2,6
    // edge 2 faces z0,x0 nodes 0,2
    // edge 3 faces x0,z1 nodes 4,6
    // edge 4 faces y0,x1 nodes 1,5
    // edge 5 faces x1,y1 nodes 3,7
    // edge 6 faces x1,z0 nodes 1,3
    // edge 7 faces z1,x1 nodes 5,7
    // edge 8 faces y0,z0 nodes 0,1
    // edge 9 faces z1,y0 nodes 4,5
    // edge 10 faces z0,y1 nodes 2,3
    // edge 11 faces y1,z1 nodes 6,7

    vector<pair<double,int> > cutCubeNodeVec[12];
    bool got_ccn = false;
    
    for (int ino = 0; ino < nno_orig; ++ino) {
      if (no_flag[ino] == ((1<<(-X0_FACE-2))|(1<<(-Y0_FACE-2+6)))) {
	cutCubeNodeVec[0].push_back(pair<double,int>(x_no[ino][2],ino));
	got_ccn = true;
      }
      else if (no_flag[ino] == ((1<<(-X0_FACE-2+6))|(1<<(-Y0_FACE-2)))) {
	cutCubeNodeVec[0].push_back(pair<double,int>(x_no[ino][2],-ino-1));
	got_ccn = true;
      }
      else if (no_flag[ino] == ((1<<(-Y1_FACE-2))|(1<<(-X0_FACE-2+6)))) {
	cutCubeNodeVec[1].push_back(pair<double,int>(x_no[ino][2],ino));
	got_ccn = true;
      }
      else if (no_flag[ino] == ((1<<(-Y1_FACE-2+6))|(1<<(-X0_FACE-2)))) {
	cutCubeNodeVec[1].push_back(pair<double,int>(x_no[ino][2],-ino-1));
	got_ccn = true;
      }
      else if (no_flag[ino] == ((1<<(-Z0_FACE-2))|(1<<(-X0_FACE-2+6)))) {
	cutCubeNodeVec[2].push_back(pair<double,int>(x_no[ino][1],ino));
	got_ccn = true;
      }
      else if (no_flag[ino] == ((1<<(-Z0_FACE-2+6))|(1<<(-X0_FACE-2)))) {
	cutCubeNodeVec[2].push_back(pair<double,int>(x_no[ino][1],-ino-1));
	got_ccn = true;
      }
      else if (no_flag[ino] == ((1<<(-X0_FACE-2))|(1<<(-Z1_FACE-2+6)))) {
	cutCubeNodeVec[3].push_back(pair<double,int>(x_no[ino][1],ino));
	got_ccn = true;
      }
      else if (no_flag[ino] == ((1<<(-X0_FACE-2+6))|(1<<(-Z1_FACE-2)))) {
	cutCubeNodeVec[3].push_back(pair<double,int>(x_no[ino][1],-ino-1));
	got_ccn = true;
      }
      else if (no_flag[ino] == ((1<<(-Y0_FACE-2))|(1<<(-X1_FACE-2+6)))) {
	cutCubeNodeVec[4].push_back(pair<double,int>(x_no[ino][2],ino));
	got_ccn = true;
      }
      else if (no_flag[ino] == ((1<<(-Y0_FACE-2+6))|(1<<(-X1_FACE-2)))) {
	cutCubeNodeVec[4].push_back(pair<double,int>(x_no[ino][2],-ino-1));
	got_ccn = true;
      }
      else if (no_flag[ino] == ((1<<(-X1_FACE-2))|(1<<(-Y1_FACE-2+6)))) {
	cutCubeNodeVec[5].push_back(pair<double,int>(x_no[ino][2],ino));
	got_ccn = true;
      }
      else if (no_flag[ino] == ((1<<(-X1_FACE-2+6))|(1<<(-Y1_FACE-2)))) {
	cutCubeNodeVec[5].push_back(pair<double,int>(x_no[ino][2],-ino-1));
	got_ccn = true;
      }
      else if (no_flag[ino] == ((1<<(-X1_FACE-2))|(1<<(-Z0_FACE-2+6)))) {
	cutCubeNodeVec[6].push_back(pair<double,int>(x_no[ino][1],ino));
	got_ccn = true;
      }
      else if (no_flag[ino] == ((1<<(-X1_FACE-2+6))|(1<<(-Z0_FACE-2)))) {
	cutCubeNodeVec[6].push_back(pair<double,int>(x_no[ino][1],-ino-1));
	got_ccn = true;
      }
      else if (no_flag[ino] == ((1<<(-Z1_FACE-2))|(1<<(-X1_FACE-2+6)))) {
	cutCubeNodeVec[7].push_back(pair<double,int>(x_no[ino][1],ino));
	got_ccn = true;
      }
      else if (no_flag[ino] == ((1<<(-Z1_FACE-2+6))|(1<<(-X1_FACE-2)))) {
	cutCubeNodeVec[7].push_back(pair<double,int>(x_no[ino][1],-ino-1));
	got_ccn = true;
      }
      else if (no_flag[ino] == ((1<<(-Y0_FACE-2))|(1<<(-Z0_FACE-2+6)))) {
	cutCubeNodeVec[8].push_back(pair<double,int>(x_no[ino][0],ino));
	got_ccn = true;
      }
      else if (no_flag[ino] == ((1<<(-Y0_FACE-2+6))|(1<<(-Z0_FACE-2)))) {
	cutCubeNodeVec[8].push_back(pair<double,int>(x_no[ino][0],-ino-1));
	got_ccn = true;
      }
      else if (no_flag[ino] == ((1<<(-Z1_FACE-2))|(1<<(-Y0_FACE-2+6)))) {
	cutCubeNodeVec[9].push_back(pair<double,int>(x_no[ino][0],ino));
	got_ccn = true;
      }
      else if (no_flag[ino] == ((1<<(-Z1_FACE-2+6))|(1<<(-Y0_FACE-2)))) {
	cutCubeNodeVec[9].push_back(pair<double,int>(x_no[ino][0],-ino-1));
	got_ccn = true;
      }
      else if (no_flag[ino] == ((1<<(-Z0_FACE-2))|(1<<(-Y1_FACE-2+6)))) {
	cutCubeNodeVec[10].push_back(pair<double,int>(x_no[ino][0],ino));
	got_ccn = true;
      }
      else if (no_flag[ino] == ((1<<(-Z0_FACE-2+6))|(1<<(-Y1_FACE-2)))) {
	cutCubeNodeVec[10].push_back(pair<double,int>(x_no[ino][0],-ino-1));
	got_ccn = true;
      }
      else if (no_flag[ino] == ((1<<(-Y1_FACE-2))|(1<<(-Z1_FACE-2+6)))) {
	cutCubeNodeVec[11].push_back(pair<double,int>(x_no[ino][0],ino));
	got_ccn = true;
      }
      else if (no_flag[ino] == ((1<<(-Y1_FACE-2+6))|(1<<(-Z1_FACE-2)))) {
	cutCubeNodeVec[11].push_back(pair<double,int>(x_no[ino][0],-ino-1));
	got_ccn = true;
      }
    }

    /*
      FILE * fp = fopen("edge_intersection.dat","w");
      for (int i = 0; i < 12; ++i) {
      cout << "cutCubeNodeVec[" << i << "].size(): " << cutCubeNodeVec[i].size() << endl;
      for (int j = 0; j < cutCubeNodeVec[i].size(); ++j) {
      const double L = cutCubeNodeVec[i][j].first;
      int ino = cutCubeNodeVec[i][j].second;
      cout << " > L,ino: " << L << " " << ino << endl;
      if (ino < 0)
      ino = -ino-1;
      fprintf(fp,"%18.15e %18.15e %18.15e\n",x_no[ino][0],x_no[ino][1],x_no[ino][2]);
      }
      }
      fclose(fp);
      cout << "checkout edge_intersection.dat" << endl; 
      getchar();
    */

    if (!got_ccn) {

      // there were no intersections of the cut surface with any of the 12 edges...
      
      // if we didn't find any cube intersections, then it may still be 
      // necessary to include the cube...

      double face_area[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
      bool got_edge = false;
      
      FOR_IED {
        if (faoed[ied][0] < -1) {
          // even one edge will be enough for us to determine the surface orientation...
          got_edge = true;
          const int ic = -faoed[ied][0]-2;
          assert((ic >= 0)&&(ic < 6)); // -x,+x,-y,+y,-z,+z
          const int ino0 = nooed[ied][0]; assert(ino0 >= 0);
          const int ino1 = nooed[ied][1]; assert(ino1 >= 0);
          switch (ic) {
          case 0:
            // surface on -x...
            //assert( fabs(x_no[ino0][0] - Lxm) < 1.0E-10 );
            //assert( fabs(x_no[ino1][0] - Lxm) < 1.0E-10 );
            face_area[0] += x_no[ino0][1]*x_no[ino1][2] - x_no[ino0][2]*x_no[ino1][1];
            break;
          case 1:
            // surface on +x...
            //assert( fabs(x_no[ino0][0] - Lxp) < 1.0E-10 );
            //assert( fabs(x_no[ino1][0] - Lxp) < 1.0E-10 );
            face_area[1] -= x_no[ino0][1]*x_no[ino1][2] - x_no[ino0][2]*x_no[ino1][1];
            break;
          case 2:
            // surface on -y...
            //assert( fabs(x_no[ino0][1] - Lym) < 1.0E-10 );
            //assert( fabs(x_no[ino1][1] - Lym) < 1.0E-10 );
            face_area[2] += x_no[ino0][2]*x_no[ino1][0] - x_no[ino0][0]*x_no[ino1][2];
            break;
          case 3:
            // surface on +y...
            //assert( fabs(x_no[ino0][1] - Lyp) < 1.0E-10 );
            //assert( fabs(x_no[ino1][1] - Lyp) < 1.0E-10 );
            face_area[3] -= x_no[ino0][2]*x_no[ino1][0] - x_no[ino0][0]*x_no[ino1][2];
            break;
          case 4:
            // surface on -z...
            //assert( fabs(x_no[ino0][2] - Lzm) < 1.0E-10 );
            //assert( fabs(x_no[ino1][2] - Lzm) < 1.0E-10 );
            face_area[4] += x_no[ino0][0]*x_no[ino1][1] - x_no[ino0][1]*x_no[ino1][0];
            break;
          case 5:
            // surface on +z...
            //assert( fabs(x_no[ino0][2] - Lzp) < 1.0E-10 );
            //assert( fabs(x_no[ino1][2] - Lzp) < 1.0E-10 );
            face_area[5] -= x_no[ino0][0]*x_no[ino1][1] - x_no[ino0][1]*x_no[ino1][0];
            break;
          default:
            assert(0);
          }
        }
      }
      
      if (!got_edge) {

        double V_sum = 0.0;
        map<const int,double*> xMap;
        FOR_IED {
          map<const int,double*>::iterator iter = xMap.find(faoed[ied][0]);
          if (iter == xMap.end()) {
            xMap.insert(pair<int,double*>(faoed[ied][0],x_no[nooed[ied][0]]));
          }
          else if (x_no[nooed[ied][1]] != iter->second) {
            V_sum += CROSS_DOT(iter->second,x_no[nooed[ied][0]],x_no[nooed[ied][1]]);
          }
          iter = xMap.find(faoed[ied][1]);
          if (iter == xMap.end()) {
            xMap.insert(pair<int,double*>(faoed[ied][1],x_no[nooed[ied][1]]));
          }
          else if (x_no[nooed[ied][0]] != iter->second) {
            V_sum += CROSS_DOT(iter->second,x_no[nooed[ied][1]],x_no[nooed[ied][0]]);
          }
        }

        if (V_sum < 0.0) {
        
          // here we have grabbed an isolated void. Add the cube around it and continue...

          assert(nno == nno_orig);
          int ino;
          ino = new_node(); x_no[ino][0] = -L; x_no[ino][1] = -L; x_no[ino][2] = -L;
          ino = new_node(); x_no[ino][0] = L; x_no[ino][1] = -L; x_no[ino][2] = -L;
          ino = new_node(); x_no[ino][0] = -L; x_no[ino][1] = L; x_no[ino][2] = -L; 
          ino = new_node(); x_no[ino][0] = L; x_no[ino][1] = L; x_no[ino][2] = -L; 
          ino = new_node(); x_no[ino][0] = -L; x_no[ino][1] = -L; x_no[ino][2] = L;
          ino = new_node(); x_no[ino][0] = L; x_no[ino][1] = -L; x_no[ino][2] = L;
          ino = new_node(); x_no[ino][0] = -L; x_no[ino][1] = L; x_no[ino][2] = L; 
          ino = new_node(); x_no[ino][0] = L; x_no[ino][1] = L; x_no[ino][2] = L;
        
          assert(ned == ned_orig);
          int ied;
          // x-edges...
          ied = new_edge(); nooed[ied][0] = nno_orig+0; nooed[ied][1] = nno_orig+1; faoed[ied][0] = Y0_FACE; faoed[ied][1] = Z0_FACE;
          ied = new_edge(); nooed[ied][0] = nno_orig+2; nooed[ied][1] = nno_orig+3; faoed[ied][0] = Z0_FACE; faoed[ied][1] = Y1_FACE;
          ied = new_edge(); nooed[ied][0] = nno_orig+4; nooed[ied][1] = nno_orig+5; faoed[ied][0] = Z1_FACE; faoed[ied][1] = Y0_FACE;
          ied = new_edge(); nooed[ied][0] = nno_orig+6; nooed[ied][1] = nno_orig+7; faoed[ied][0] = Y1_FACE; faoed[ied][1] = Z1_FACE;
          // y-edges...
          ied = new_edge(); nooed[ied][0] = nno_orig+0; nooed[ied][1] = nno_orig+2; faoed[ied][0] = Z0_FACE; faoed[ied][1] = X0_FACE;
          ied = new_edge(); nooed[ied][0] = nno_orig+1; nooed[ied][1] = nno_orig+3; faoed[ied][0] = X1_FACE; faoed[ied][1] = Z0_FACE;
          ied = new_edge(); nooed[ied][0] = nno_orig+4; nooed[ied][1] = nno_orig+6; faoed[ied][0] = X0_FACE; faoed[ied][1] = Z1_FACE;
          ied = new_edge(); nooed[ied][0] = nno_orig+5; nooed[ied][1] = nno_orig+7; faoed[ied][0] = Z1_FACE; faoed[ied][1] = X1_FACE;
          // z-edges...
          ied = new_edge(); nooed[ied][0] = nno_orig+0; nooed[ied][1] = nno_orig+4; faoed[ied][0] = X0_FACE; faoed[ied][1] = Y0_FACE;
          ied = new_edge(); nooed[ied][0] = nno_orig+1; nooed[ied][1] = nno_orig+5; faoed[ied][0] = Y0_FACE; faoed[ied][1] = X1_FACE;
          ied = new_edge(); nooed[ied][0] = nno_orig+2; nooed[ied][1] = nno_orig+6; faoed[ied][0] = Y1_FACE; faoed[ied][1] = X0_FACE;
          ied = new_edge(); nooed[ied][0] = nno_orig+3; nooed[ied][1] = nno_orig+7; faoed[ied][0] = X1_FACE; faoed[ied][1] = Y1_FACE;
        
        }
        
        // no need to add the outer cube: a positive V_sum is simply a very coarse grid...
        
      }
      else if ( (face_area[0] >= 0.0) && (face_area[1] >= 0.0) && (face_area[2] >= 0.0) &&
                (face_area[3] >= 0.0) && (face_area[4] >= 0.0) && (face_area[5] >= 0.0) ) {

        // the surface is penetrating into the cube. We need to add the full outer cube...
        //cout << "add all cube nodes and edges" << endl;

        assert(nno == nno_orig);
        int ino;
        ino = new_node(); x_no[ino][0] = -L; x_no[ino][1] = -L; x_no[ino][2] = -L;
        ino = new_node(); x_no[ino][0] = L; x_no[ino][1] = -L; x_no[ino][2] = -L;
        ino = new_node(); x_no[ino][0] = -L; x_no[ino][1] = L; x_no[ino][2] = -L; 
        ino = new_node(); x_no[ino][0] = L; x_no[ino][1] = L; x_no[ino][2] = -L; 
        ino = new_node(); x_no[ino][0] = -L; x_no[ino][1] = -L; x_no[ino][2] = L;
        ino = new_node(); x_no[ino][0] = L; x_no[ino][1] = -L; x_no[ino][2] = L;
        ino = new_node(); x_no[ino][0] = -L; x_no[ino][1] = L; x_no[ino][2] = L; 
        ino = new_node(); x_no[ino][0] = L; x_no[ino][1] = L; x_no[ino][2] = L;
        
        assert(ned == ned_orig);
        int ied;
        // x-edges...
        ied = new_edge(); nooed[ied][0] = nno_orig+0; nooed[ied][1] = nno_orig+1; faoed[ied][0] = Y0_FACE; faoed[ied][1] = Z0_FACE;
        ied = new_edge(); nooed[ied][0] = nno_orig+2; nooed[ied][1] = nno_orig+3; faoed[ied][0] = Z0_FACE; faoed[ied][1] = Y1_FACE;
        ied = new_edge(); nooed[ied][0] = nno_orig+4; nooed[ied][1] = nno_orig+5; faoed[ied][0] = Z1_FACE; faoed[ied][1] = Y0_FACE;
        ied = new_edge(); nooed[ied][0] = nno_orig+6; nooed[ied][1] = nno_orig+7; faoed[ied][0] = Y1_FACE; faoed[ied][1] = Z1_FACE;
        // y-edges...
        ied = new_edge(); nooed[ied][0] = nno_orig+0; nooed[ied][1] = nno_orig+2; faoed[ied][0] = Z0_FACE; faoed[ied][1] = X0_FACE;
        ied = new_edge(); nooed[ied][0] = nno_orig+1; nooed[ied][1] = nno_orig+3; faoed[ied][0] = X1_FACE; faoed[ied][1] = Z0_FACE;
        ied = new_edge(); nooed[ied][0] = nno_orig+4; nooed[ied][1] = nno_orig+6; faoed[ied][0] = X0_FACE; faoed[ied][1] = Z1_FACE;
        ied = new_edge(); nooed[ied][0] = nno_orig+5; nooed[ied][1] = nno_orig+7; faoed[ied][0] = Z1_FACE; faoed[ied][1] = X1_FACE;
        // z-edges...
        ied = new_edge(); nooed[ied][0] = nno_orig+0; nooed[ied][1] = nno_orig+4; faoed[ied][0] = X0_FACE; faoed[ied][1] = Y0_FACE;
        ied = new_edge(); nooed[ied][0] = nno_orig+1; nooed[ied][1] = nno_orig+5; faoed[ied][0] = Y0_FACE; faoed[ied][1] = X1_FACE;
        ied = new_edge(); nooed[ied][0] = nno_orig+2; nooed[ied][1] = nno_orig+6; faoed[ied][0] = Y1_FACE; faoed[ied][1] = X0_FACE;
        ied = new_edge(); nooed[ied][0] = nno_orig+3; nooed[ied][1] = nno_orig+7; faoed[ied][0] = X1_FACE; faoed[ied][1] = Y1_FACE;
      
      }
      else {
        
        // no need to add the outer cube, but 
        // make sure all areas have the expected sign...

        if (!( (face_area[0] <= 0.0) && (face_area[1] <= 0.0) && (face_area[2] <= 0.0) &&
               (face_area[3] <= 0.0) && (face_area[4] <= 0.0) && (face_area[5] <= 0.0) ) ) {
          cout << "Error: CuttableVoronoiData: failed face_area check: " << 
            face_area[0] << " " << 
            face_area[1] << " " << 
            face_area[2] << " " << 
            face_area[3] << " " << 
            face_area[4] << " " << 
            face_area[5] << endl;
          throw(-1);
        }
        
      }

    }
    else {

      // set a corner_flag from 1 to 0 (in) or -1 (out) based on the orientation of
      // the edge intersections...
      
      int corner_flag[8];
      FOR_I8 corner_flag[i] = 1;
    
      // edge 0 faces x0,y0 nodes 0,4
#define IED 0
#define INO0 0
#define INO1 4
#include "corner_flag1.hpp"
#undef IED
#undef INO0
#undef INO1
      // edge 1 faces y1,x0 nodes 2,6
#define IED 1
#define INO0 2
#define INO1 6
#include "corner_flag1.hpp"
#undef IED
#undef INO0
#undef INO1
      // edge 2 faces z0,x0 nodes 0,2
#define IED 2
#define INO0 0
#define INO1 2
#include "corner_flag1.hpp"
#undef IED
#undef INO0
#undef INO1
      // edge 3 faces x0,z1 nodes 4,6
#define IED 3
#define INO0 4
#define INO1 6
#include "corner_flag1.hpp"
#undef IED
#undef INO0
#undef INO1
      // edge 4 faces y0,x1 nodes 1,5
#define IED 4
#define INO0 1
#define INO1 5
#include "corner_flag1.hpp"
#undef IED
#undef INO0
#undef INO1
      // edge 5 faces x1,y1 nodes 3,7
#define IED 5
#define INO0 3
#define INO1 7
#include "corner_flag1.hpp"
#undef IED
#undef INO0
#undef INO1
      // edge 6 faces x1,z0 nodes 1,3
#define IED 6
#define INO0 1
#define INO1 3
#include "corner_flag1.hpp"
#undef IED
#undef INO0
#undef INO1
      // edge 7 faces z1,x1 nodes 5,7
#define IED 7
#define INO0 5
#define INO1 7
#include "corner_flag1.hpp"
#undef IED
#undef INO0
#undef INO1
      // edge 8 faces y0,z0 nodes 0,1
#define IED 8
#define INO0 0
#define INO1 1
#include "corner_flag1.hpp"
#undef IED
#undef INO0
#undef INO1
      // edge 9 faces z1,y0 nodes 4,5
#define IED 9
#define INO0 4
#define INO1 5
#include "corner_flag1.hpp"
#undef IED
#undef INO0
#undef INO1
      // edge 10 faces z0,y1 nodes 2,3
#define IED 10
#define INO0 2
#define INO1 3
#include "corner_flag1.hpp"
#undef IED
#undef INO0
#undef INO1
      // edge 11 faces y1,z1 nodes 6,7
#define IED 11
#define INO0 6
#define INO1 7
#include "corner_flag1.hpp"
#undef IED
#undef INO0
#undef INO1

      // and add the corner nodes (we may need some of these)...
      resize_nno(nno_orig+8);
      assert(nno == nno_orig+8);
      x_no[nno_orig  ][0] = -L; x_no[nno_orig  ][1] = -L; x_no[nno_orig  ][2] = -L;
      x_no[nno_orig+1][0] = L; x_no[nno_orig+1][1] = -L; x_no[nno_orig+1][2] = -L;
      x_no[nno_orig+2][0] = -L; x_no[nno_orig+2][1] = L; x_no[nno_orig+2][2] = -L; 
      x_no[nno_orig+3][0] = L; x_no[nno_orig+3][1] = L; x_no[nno_orig+3][2] = -L; 
      x_no[nno_orig+4][0] = -L; x_no[nno_orig+4][1] = -L; x_no[nno_orig+4][2] = L;
      x_no[nno_orig+5][0] = L; x_no[nno_orig+5][1] = -L; x_no[nno_orig+5][2] = L;
      x_no[nno_orig+6][0] = -L; x_no[nno_orig+6][1] = L; x_no[nno_orig+6][2] = L; 
      x_no[nno_orig+7][0] = L; x_no[nno_orig+7][1] = L; x_no[nno_orig+7][2] = L;
      
      // and then loop edges to set all nodes...
      int done = 0;
      while (done == 0) {
	done = 1; // assume done
	// edge 0 faces x0,y0 nodes 0,4
#define IED 0
#define INO0 0
#define INO1 4
#include "corner_flag2.hpp"
#undef IED
#undef INO0
#undef INO1
	// edge 1 faces y1,x0 nodes 2,6
#define IED 1
#define INO0 2
#define INO1 6
#include "corner_flag2.hpp"
#undef IED
#undef INO0
#undef INO1
	// edge 2 faces z0,x0 nodes 0,2
#define IED 2
#define INO0 0
#define INO1 2
#include "corner_flag2.hpp"
#undef IED
#undef INO0
#undef INO1
	// edge 3 faces x0,z1 nodes 4,6
#define IED 3
#define INO0 4
#define INO1 6
#include "corner_flag2.hpp"
#undef IED
#undef INO0
#undef INO1
	// edge 4 faces y0,x1 nodes 1,5
#define IED 4
#define INO0 1
#define INO1 5
#include "corner_flag2.hpp"
#undef IED
#undef INO0
#undef INO1
	// edge 5 faces x1,y1 nodes 3,7
#define IED 5
#define INO0 3
#define INO1 7
#include "corner_flag2.hpp"
#undef IED
#undef INO0
#undef INO1
	// edge 6 faces x1,z0 nodes 1,3
#define IED 6
#define INO0 1
#define INO1 3
#include "corner_flag2.hpp"
#undef IED
#undef INO0
#undef INO1
	// edge 7 faces z1,x1 nodes 5,7
#define IED 7
#define INO0 5
#define INO1 7
#include "corner_flag2.hpp"
#undef IED
#undef INO0
#undef INO1
	// edge 8 faces y0,z0 nodes 0,1
#define IED 8
#define INO0 0
#define INO1 1
#include "corner_flag2.hpp"
#undef IED
#undef INO0
#undef INO1
	// edge 9 faces z1,y0 nodes 4,5
#define IED 9
#define INO0 4
#define INO1 5
#include "corner_flag2.hpp"
#undef IED
#undef INO0
#undef INO1
	// edge 10 faces z0,y1 nodes 2,3
#define IED 10
#define INO0 2
#define INO1 3
#include "corner_flag2.hpp"
#undef IED
#undef INO0
#undef INO1
	// edge 11 faces y1,z1 nodes 6,7
#define IED 11
#define INO0 6
#define INO1 7
#include "corner_flag2.hpp"
#undef IED
#undef INO0
#undef INO1
      }
    
      // at this point, the corners that will remain are flagged with 0. Use this to add 
      // all required edges along the sides of the cube...

      // edge 0 faces x0,y0 nodes 0,4
#define IED 0
#define INO0 0
#define INO1 4
#define IFA0 X0_FACE
#define IFA1 Y0_FACE
#include "corner_flag3.hpp"
#undef IED
#undef INO0
#undef INO1
#undef IFA0
#undef IFA1
      // edge 1 faces y1,x0 nodes 2,6
#define IED 1
#define INO0 2
#define INO1 6
#define IFA0 Y1_FACE
#define IFA1 X0_FACE
#include "corner_flag3.hpp"
#undef IED
#undef INO0
#undef INO1
#undef IFA0
#undef IFA1
      // edge 2 faces z0,x0 nodes 0,2
#define IED 2
#define INO0 0
#define INO1 2
#define IFA0 Z0_FACE
#define IFA1 X0_FACE
#include "corner_flag3.hpp"
#undef IED
#undef INO0
#undef INO1
#undef IFA0
#undef IFA1
      // edge 3 faces x0,z1 nodes 4,6
#define IED 3
#define INO0 4
#define INO1 6
#define IFA0 X0_FACE
#define IFA1 Z1_FACE
#include "corner_flag3.hpp"
#undef IED
#undef INO0
#undef INO1
#undef IFA0
#undef IFA1
      // edge 4 faces y0,x1 nodes 1,5
#define IED 4
#define INO0 1
#define INO1 5
#define IFA0 Y0_FACE
#define IFA1 X1_FACE
#include "corner_flag3.hpp"
#undef IED
#undef INO0
#undef INO1
#undef IFA0
#undef IFA1
      // edge 5 faces x1,y1 nodes 3,7
#define IED 5
#define INO0 3
#define INO1 7
#define IFA0 X1_FACE
#define IFA1 Y1_FACE
#include "corner_flag3.hpp"
#undef IED
#undef INO0
#undef INO1
#undef IFA0
#undef IFA1
      // edge 6 faces x1,z0 nodes 1,3
#define IED 6
#define INO0 1
#define INO1 3
#define IFA0 X1_FACE
#define IFA1 Z0_FACE
#include "corner_flag3.hpp"
#undef IED
#undef INO0
#undef INO1
#undef IFA0
#undef IFA1
      // edge 7 faces z1,x1 nodes 5,7
#define IED 7
#define INO0 5
#define INO1 7
#define IFA0 Z1_FACE
#define IFA1 X1_FACE
#include "corner_flag3.hpp"
#undef IED
#undef INO0
#undef INO1
#undef IFA0
#undef IFA1
      // edge 8 faces y0,z0 nodes 0,1
#define IED 8
#define INO0 0
#define INO1 1
#define IFA0 Y0_FACE
#define IFA1 Z0_FACE
#include "corner_flag3.hpp"
#undef IED
#undef INO0
#undef INO1
#undef IFA0
#undef IFA1
      // edge 9 faces z1,y0 nodes 4,5
#define IED 9
#define INO0 4
#define INO1 5
#define IFA0 Z1_FACE
#define IFA1 Y0_FACE
#include "corner_flag3.hpp"
#undef IED
#undef INO0
#undef INO1
#undef IFA0
#undef IFA1
      // edge 10 faces z0,y1 nodes 2,3
#define IED 10
#define INO0 2
#define INO1 3
#define IFA0 Z0_FACE
#define IFA1 Y1_FACE
#include "corner_flag3.hpp"
#undef IED
#undef INO0
#undef INO1
#undef IFA0
#undef IFA1
      // edge 11 faces y1,z1 nodes 6,7
#define IED 11
#define INO0 6
#define INO1 7
#define IFA0 Y1_FACE
#define IFA1 Z1_FACE
#include "corner_flag3.hpp"
#undef IED
#undef INO0
#undef INO1
#undef IFA0
#undef IFA1

      // we must have added something...
      assert(ned > ned_orig);

      /*
	FILE * fp = fopen("corners.dat","w");
	FOR_I8 {
	cout << "corner_flag i: " << i << " " << corner_flag[i] << endl;
	if (corner_flag[i] == 0) {
	const int ino = nno_orig+i;
	fprintf(fp,"%18.15le %18.15le %18.15le\n",x_no[ino][0],x_no[ino][1],x_no[ino][2]);
	}
	else {
	assert(corner_flag[i] == -1);
	}
	}
	fclose(fp);
	writeTecplot(-2);
	cout << "checkout corners.dat" << endl;
	getchar();
      */
    
    }

  }
    
  // at this point we have all the outer edges necessary. If there are multiple parts
  // involved, we may need to process intersections, so do this next...
    
  if (debug) {
    writeTecplot(-2);
    cout << "original surface + cube edges in tecplot -2" << endl;
  }
  
  // do a fast check to see if we need to do a costly intersections between different parts...
  set<int> intSet;
  for (int ied = 0; ied < ned_orig; ++ied) {
    if (faoed[ied][1] >= 0) {
      const int ipart = (ist_and_ipart_buf[faoed[ied][1]+1]>>6); // recall first entry is ist_global, second is ipart/bits
      intSet.insert(ipart);
    }
    if (faoed[ied][0] >= 0) {
      const int ipart = (ist_and_ipart_buf[faoed[ied][0]+1]>>6); // recall first entry is ist_global, second is ipart/bits
      intSet.insert(ipart);
    }
  }
    
  // ===================================================================================
  // solve for intersections between parts, introducing new edges and 
  // nodes so that parts are walkable...
  // ===================================================================================
  
  if (debug) cout << "got intSet.size() (number of unique parts): " << intSet.size() << endl;

  if (intSet.size() == 1) {

    // only one part. For now, assume it is manifold and take everything. In the 
    // future we could have identified certain problems and cleaned the 
    // geometry.
    // the only thing to do here is to remove any nodes from the 8 that 
    // were added for the cube that are not touched by edges...

    if (nno > nno_orig) {
      
      assert(nno == nno_orig+8);
      
      int no_flag[8];
      FOR_I8 no_flag[i] = -1;
      
      // loop on new edges and mark any nodes touched by these edges...
      
      for (int ied = ned_orig; ied < ned; ++ied) {
	FOR_I2 {
	  if (nooed[ied][i] >= nno_orig) {
	    assert(nooed[ied][i] < nno_orig+8);
	    no_flag[nooed[ied][i]-nno_orig] = 0;
	  }
	}
      }
      
      nno = nno_orig;
      FOR_I8 {
	if (no_flag[i] == 0) {
	  no_flag[i] = nno;
	  if (nno < nno_orig+i) {
	    FOR_J3 x_no[nno][j] = x_no[nno_orig+i][j];
	  }
	  ++nno;
	}
      }
      assert(nno <= nno_orig+8); 
      
      for (int ied = ned_orig; ied < ned; ++ied) {
	FOR_I2 {
	  if (nooed[ied][i] >= nno_orig) {
	    const int ino_old = nooed[ied][i];
	    nooed[ied][i] = no_flag[ino_old-nno_orig];
	  }
	}
      }
      
    }
    
  }
  else {

    assert(intSet.size() > 1);
    
    // multi-part routines -- expensive!...
     
    const int L_int = (1<<19); 
    //cout << "L_int: " << L_int << endl;

    // the basis for this intersection routine is the edge-polygon intersection routine, except on
    // cube edges, where we use edge-edge type intersections...
      
    // here we are going to intersect the parts successively:
    // 1+2 => (1+2)
    // (1+2)+3 => ((1+2)+3)
    // etc...

    vector<pair<int,int> > facePairVec;
    vector<pair<int,int> > faceEdgeVec0;
    vector<pair<int,int> > faceEdgeVec1;
    vector<int> edofa_i0;
    vector<int> edofa_i1;
    vector<int> intVec;
    vector<IntersectionData> intersectionVec;
    vector<IntersectionData> localIntersectionVec; // used for a particular ifa0/ifa1 pair
    vector<pair<int,int> > intersectionEdgeRange;
    vector<int> intersectionNodeNode;
    vector<pair<pair<int,double>,int> > cutEdgeDataVec;

    assert(nno <= nno_max);
    int (*x_no_int)[3] = new int[nno_max][3];
      
    for (set<int>::iterator ipart_iter = intSet.begin(); ipart_iter != intSet.end(); ++ipart_iter) {
	
      // skip the first...
      if (ipart_iter == intSet.begin())
	continue;
	
      // make sure arrays are big enough. When nodes are added with multiple intersecting parts,
      // we may need some fancier memory management here...
      assert(nno <= nno_max);
      FOR_INO no_flag[ino] = -1;
	
      // we are going to intersect part1 with everything lower...
      const int ipart1 = *ipart_iter;
      if (!(ipart1 >= 1)) {
        cout << "unexpected ipart1: " << ipart1 << endl;
        writeTecplot(-13);
        throw(-1);
      }
      assert(ipart1 >= 1);
	
      assert(intersectionVec.empty());
      assert(faceEdgeVec0.empty());
      assert(faceEdgeVec1.empty());
      for (int ied = 0; ied < ned; ++ied) {
	// faoed[ied][0,1] >= 0 indicates a tri. It is possible that 
	// faoed[ied][1] == -1, in which case it is an exterior edge of the geometry
	// and should be cut away by the intersection process.
	// faoed[ied][0] can be negative, (-2..-7) indicating an edge along one of the cube faces...
	if (faoed[ied][0] >= 0) {
	  const int ipart = (ist_and_ipart_buf[faoed[ied][0]+1]>>6); // recall first entry is ist_global, second is ipart/bits
	  if (ipart < ipart1) {
	    faceEdgeVec0.push_back(pair<int,int>(faoed[ied][0],ied));
	    // and put ipart into no_flag...
	    no_flag[nooed[ied][0]] = ipart1-1;
	    no_flag[nooed[ied][1]] = ipart1-1;
	  }
	  else if (ipart == ipart1) {
	    faceEdgeVec1.push_back(pair<int,int>(faoed[ied][0],ied));
	    // and put ipart into no_flag...
	    no_flag[nooed[ied][0]] = ipart1;
	    no_flag[nooed[ied][1]] = ipart1;
	  }
	}
	if (faoed[ied][1] >= 0) {
	  const int ipart = (ist_and_ipart_buf[faoed[ied][1]+1]>>6); // recall first entry is ist_global, second is ipart/bits
	  if (ipart < ipart1) {
	    faceEdgeVec0.push_back(pair<int,int>(faoed[ied][1],-ied-1));
	    // and put ipart into no_flag...
	    no_flag[nooed[ied][0]] = ipart1-1;
	    no_flag[nooed[ied][1]] = ipart1-1;
	  }
	  else if (ipart == ipart1) {
	    faceEdgeVec1.push_back(pair<int,int>(faoed[ied][1],-ied-1));
	    // and put ipart into no_flag...
	    no_flag[nooed[ied][0]] = ipart1;
	    no_flag[nooed[ied][1]] = ipart1;
	  }
	}
      }

      // compute the integer coords only where needed...
      FOR_INO {
	if (no_flag[ino] == ipart1) {
	  FOR_I3 {
	    x_no_int[ino][i] = (int)floor(x_no[ino][i]/L*double(L_int)+0.5);
	    // HACK: try even/odd?... best to avoid this so we can use this 
	    // routine (or something like it) to repair overlapping/intersecting
	    // tris on bad surfaces some day...
	    //if (x_no_int[ino][i]&1) x_no_int[ino][i] &= ~1;
	  }
	}
	else if (no_flag[ino] == ipart1-1) {
	  FOR_I3 {
	    x_no_int[ino][i] = (int)floor(x_no[ino][i]/L*double(L_int)+0.5);
	    // HACK: try even/odd?...
	    //if (!(x_no_int[ino][i]&1)) x_no_int[ino][i] |= 1;
	  }
        }
      }
      
      // HACK: stretch out the nodes that formally touch the cube edges
      // based on their faoed. This avoids certain problematic edge-edge 
      // intersections between boundary and internal edges that should not be allowed...
      
      // NOT in this version. The bounds have been selected to specifically avoid
      // other crap near the cut cube surfaces...

      FOR_IED {
	if (faoed[ied][0] == X0_FACE) {
          FOR_I2 {
            const int ino = nooed[ied][i];
            if ((no_flag[ino] == ipart1)||(no_flag[ino] == ipart1-1)) {
              assert(x_no_int[ino][0] == -L_int);
            }
          }
	}
	else if (faoed[ied][0] == X1_FACE) {
          FOR_I2 {
            const int ino = nooed[ied][i];
            if ((no_flag[ino] == ipart1)||(no_flag[ino] == ipart1-1)) {
              assert(x_no_int[ino][0] == L_int);
            }
          }
        }
	else if (faoed[ied][0] == Y0_FACE) {
          FOR_I2 {
            const int ino = nooed[ied][i];
            if ((no_flag[ino] == ipart1)||(no_flag[ino] == ipart1-1)) {
              assert(x_no_int[ino][1] == -L_int);
            }
          }
	}
	else if (faoed[ied][0] == Y1_FACE) {
          FOR_I2 {
            const int ino = nooed[ied][i];
            if ((no_flag[ino] == ipart1)||(no_flag[ino] == ipart1-1)) {
              assert(x_no_int[ino][1] == L_int);
            }
          }
        }
	else if (faoed[ied][0] == Z0_FACE) {
          FOR_I2 {
            const int ino = nooed[ied][i];
            if ((no_flag[ino] == ipart1)||(no_flag[ino] == ipart1-1)) {
              assert(x_no_int[ino][2] == -L_int);
            }
          }
	}
	else if (faoed[ied][0] == Z1_FACE) {
          FOR_I2 {
            const int ino = nooed[ied][i];
            if ((no_flag[ino] == ipart1)||(no_flag[ino] == ipart1-1)) {
              assert(x_no_int[ino][2] == L_int);
            }
          }
        }
      }
      
      // for fast lookup, we will need an adt...
      sort(faceEdgeVec0.begin(),faceEdgeVec0.end());
      assert(edofa_i0.empty());
      edofa_i0.push_back(0);
      for (int fe0 = 1; fe0 < faceEdgeVec0.size(); ++fe0) {
	if (faceEdgeVec0[fe0].first != faceEdgeVec0[fe0-1].first)
	  edofa_i0.push_back(fe0);
      }
      const int nfa0 = edofa_i0.size();
      edofa_i0.push_back(faceEdgeVec0.size());

      int * fa0_flag = new int[nfa0];
      for (int ifa0 = 0; ifa0 < nfa0; ++ifa0) 
	fa0_flag[ifa0] = -2;
	
      // put these in an adt...
      // include the first node in each edge to expand the bbox. The first
      // node depends on the edge orientation, i.e. faceEdgeVec0[fe0].second's sign...  
      // NOTE: is this risky? makes assumtions about connectivity of surface
      // face cutting. Should be OK...
      int (*bbmin_fa0)[3] = new int[nfa0][3];
      int (*bbmax_fa0)[3] = new int[nfa0][3];
      for (int ifa0 = 0; ifa0 < nfa0; ++ifa0) {
	FOR_I3 bbmin_fa0[ifa0][i] = L_int+1;
	FOR_I3 bbmax_fa0[ifa0][i] = -L_int-1;
	for (int fe0 = edofa_i0[ifa0]; fe0 != edofa_i0[ifa0+1]; ++fe0) {
	  if (faceEdgeVec0[fe0].second >= 0) {
	    // faceEdgeVec0[fe0].second stores the edge index normally...
	    const int ino = nooed[faceEdgeVec0[fe0].second][0]; // node 0
	    FOR_I3 {
	      bbmin_fa0[ifa0][i] = min(bbmin_fa0[ifa0][i],x_no_int[ino][i]);
	      bbmax_fa0[ifa0][i] = max(bbmax_fa0[ifa0][i],x_no_int[ino][i]);
	    }
	  }
	  else {
	    // faceEdgeVec0[fe0].second stores the edge index using -1 indexing...
	    const int ino = nooed[-faceEdgeVec0[fe0].second-1][1]; // node 1
	    FOR_I3 {
	      bbmin_fa0[ifa0][i] = min(bbmin_fa0[ifa0][i],x_no_int[ino][i]);
	      bbmax_fa0[ifa0][i] = max(bbmax_fa0[ifa0][i],x_no_int[ino][i]);
	    }
	  }
	}
      }
	
      Adt<int> * adt = new Adt<int>(nfa0,bbmin_fa0,bbmax_fa0);
      delete[] bbmin_fa0; bbmin_fa0 = NULL;
      delete[] bbmax_fa0; bbmax_fa0 = NULL;
      
      sort(faceEdgeVec1.begin(),faceEdgeVec1.end());
      assert(edofa_i1.empty());
      edofa_i1.push_back(0);
      for (int fe1 = 1; fe1 < faceEdgeVec1.size(); ++fe1) {
	if (faceEdgeVec1[fe1].first != faceEdgeVec1[fe1-1].first)
	  edofa_i1.push_back(fe1);
      }
      const int nfa1 = edofa_i1.size();
      edofa_i1.push_back(faceEdgeVec1.size());
	
      int * fa1_flag = new int[nfa1];
      for (int ifa1 = 0; ifa1 < nfa1; ++ifa1) 
	fa1_flag[ifa1] = -2;
      
      // figure out what faces intersect based on bbox check...
	
      assert(facePairVec.empty());
      for (int ifa1 = 0; ifa1 < nfa1; ++ifa1) {
	// we need the bbox of this face to query the adt...
	int bbmin_fa1[3] = { L_int+1, L_int+1, L_int+1 };
	int bbmax_fa1[3] = { -L_int-1, -L_int-1, -L_int-1 };
	for (int fe1 = edofa_i1[ifa1]; fe1 != edofa_i1[ifa1+1]; ++fe1) {
	  if (faceEdgeVec1[fe1].second >= 0) {
	    // faceEdgeVec1[fe1].second stores the edge index normally...
	    const int ino = nooed[faceEdgeVec1[fe1].second][0]; // node 0
	    FOR_I3 {
	      bbmin_fa1[i] = min(bbmin_fa1[i],x_no_int[ino][i]);
	      bbmax_fa1[i] = max(bbmax_fa1[i],x_no_int[ino][i]);
	    }
	  }
	  else {
	    // faceEdgeVec1[fe1].second stores the edge index using -1 indexing...
	    const int ino = nooed[-faceEdgeVec1[fe1].second-1][1]; // node 1
	    FOR_I3 {
	      bbmin_fa1[i] = min(bbmin_fa1[i],x_no_int[ino][i]);
	      bbmax_fa1[i] = max(bbmax_fa1[i],x_no_int[ino][i]);
	    }
	  }
	}
	// now grab the faces may intersect based on the bbox overlap...
	assert(intVec.empty());
	adt->buildListForBBox(intVec,bbmin_fa1,bbmax_fa1);
	if (!intVec.empty()) {
	  for (int ii = 0; ii < intVec.size(); ++ii) {
	    const int ifa0 = intVec[ii];
	    facePairVec.push_back(pair<int,int>(ifa0,ifa1));
	    // also set the flags so we can triangulate polygon faces with 
	    // new edges...
	    fa0_flag[ifa0] = -1; 
	    fa1_flag[ifa1] = -1; 
	  }
	  intVec.clear();
	}
      }
      
      // done with adt...
      
      delete adt; adt = NULL;
      
      // now ensure all faces that will be compared can be described as triangles. We
      // will need to break up faces that are more than 3 edges. For now these routines
      // assume these faces are simple convex loops, but this is not strictly the case
      // for certain multiple cut surface problems (e.g. tube cuts large tri, then the 
      // resulting intersected surface needs to cut something again)...
      
      for (int ifa0 = 0; ifa0 < nfa0; ++ifa0) {
	if (fa0_flag[ifa0] == -1) {
	  fa0_flag[ifa0] = triangulateFaceForDoit(ifa0,faceEdgeVec0,edofa_i0);
	}
      }
      
      for (int ifa1 = 0; ifa1 < nfa1; ++ifa1) {
	if (fa1_flag[ifa1] == -1) {
	  fa1_flag[ifa1] = triangulateFaceForDoit(ifa1,faceEdgeVec1,edofa_i1);
	}
      }

      if (debug) {
	writeTecplot(-13);
	cout << "original surface + cube edges + split polygonal faces in tecplot -13" << endl;
      }
      
      //cout << "about to loop on face pairs: " << facePairVec.size() << endl;

      // =====================================
      // now loop on face pairs...
      // =====================================
      
      const int ned_old = ned;
      for (int ifp = 0; ifp < facePairVec.size(); ++ifp) {

	const int ifa0 = facePairVec[ifp].first;
	vector<int> ied0Vec;
	packTriEdgeVecForDoit(ied0Vec,ifa0,faceEdgeVec0,edofa_i0,fa0_flag[ifa0]);
	assert(ied0Vec.size()%3 == 0);

	const int ifa1 = facePairVec[ifp].second;
	vector<int> ied1Vec;
	packTriEdgeVecForDoit(ied1Vec,ifa1,faceEdgeVec1,edofa_i1,fa1_flag[ifa1]);
	assert(ied1Vec.size()%3 == 0);
	
	for (int i0 = 0; i0 < ied0Vec.size(); i0 += 3) {
          
	  int edotr0[3];
          FOR_I3 edotr0[i] = ied0Vec[i0+i]; // leave as a signed edge

	  int nootr0[3];
	  int faoed0;
	  if (ied0Vec[i0] >= 0) {
	    nootr0[0] = nooed[ied0Vec[i0]][0];
	    nootr0[1] = nooed[ied0Vec[i0]][1];
	    faoed0 = faoed[ied0Vec[i0]][0];
	  }
	  else {
	    nootr0[0] = nooed[-ied0Vec[i0]-1][1];
	    nootr0[1] = nooed[-ied0Vec[i0]-1][0];
	    faoed0 = faoed[-ied0Vec[i0]-1][1];
	  }
	  if (ied0Vec[i0+1] >= 0) {
	    assert(nootr0[1] == nooed[ied0Vec[i0+1]][0]);
	    nootr0[2] = nooed[ied0Vec[i0+1]][1];
	    assert(faoed0 == faoed[ied0Vec[i0+1]][0]);
	  }
	  else {
	    assert(nootr0[1] == nooed[-ied0Vec[i0+1]-1][1]);
	    nootr0[2] = nooed[-ied0Vec[i0+1]-1][0];
	    assert(faoed0 == faoed[-ied0Vec[i0+1]-1][1]);
	  }
	  if (ied0Vec[i0+2] >= 0) {
	    assert(nootr0[2] == nooed[ied0Vec[i0+2]][0]);
	    assert(nootr0[0] == nooed[ied0Vec[i0+2]][1]);
	    assert(faoed0 == faoed[ied0Vec[i0+2]][0]);
	  }
	  else {
	    assert(nootr0[2] == nooed[-ied0Vec[i0+2]-1][1]);
	    assert(nootr0[0] == nooed[-ied0Vec[i0+2]-1][0]);
	    assert(faoed0 == faoed[-ied0Vec[i0+2]-1][1]);
	  }
	  
	  // modify faoed0 with a bit asociated with the tri. This prevents the 
	  // collapse of certain intersections...
	  
          assert(i0/3 <= 7);
	  faoed0 |= ((i0/3)<<28);
	  
	  // store the bbox for the tri...
	  
	  int bbmin0[3]; FOR_I3 bbmin0[i] = min(x_no_int[nootr0[0]][i],min(x_no_int[nootr0[1]][i],x_no_int[nootr0[2]][i]));
	  int bbmax0[3]; FOR_I3 bbmax0[i] = max(x_no_int[nootr0[0]][i],max(x_no_int[nootr0[1]][i],x_no_int[nootr0[2]][i]));
	  
	  for (int i1 = 0; i1 < ied1Vec.size(); i1 += 3) {
            
	    int edotr1[3];
            FOR_J3 edotr1[j] = ied1Vec[i1+j]; // leave as a signed edge
            
	    int nootr1[3];
	    int faoed1;
	    if (ied1Vec[i1] >= 0) {
	      nootr1[0] = nooed[ied1Vec[i1]][0];
	      nootr1[1] = nooed[ied1Vec[i1]][1];
	      faoed1 = faoed[ied1Vec[i1]][0];
	    }
	    else {
	      nootr1[0] = nooed[-ied1Vec[i1]-1][1];
	      nootr1[1] = nooed[-ied1Vec[i1]-1][0];
	      faoed1 = faoed[-ied1Vec[i1]-1][1];
	    }
	    if (ied1Vec[i1+1] >= 0) {
	      assert(nootr1[1] == nooed[ied1Vec[i1+1]][0]);
	      nootr1[2] = nooed[ied1Vec[i1+1]][1];
	      assert(faoed1 == faoed[ied1Vec[i1+1]][0]);
	    }
	    else {
	      assert(nootr1[1] == nooed[-ied1Vec[i1+1]-1][1]);
	      nootr1[2] = nooed[-ied1Vec[i1+1]-1][0];
	      assert(faoed1 == faoed[-ied1Vec[i1+1]-1][1]);
	    }
	    if (ied1Vec[i1+2] >= 0) {
	      assert(nootr1[2] == nooed[ied1Vec[i1+2]][0]);
	      assert(nootr1[0] == nooed[ied1Vec[i1+2]][1]);
	      assert(faoed1 == faoed[ied1Vec[i1+2]][0]);
	    }
	    else {
	      assert(nootr1[2] == nooed[-ied1Vec[i1+2]-1][1]);
	      assert(nootr1[0] == nooed[-ied1Vec[i1+2]-1][0]);
	      assert(faoed1 == faoed[-ied1Vec[i1+2]-1][1]);
	    }
	    
	    // modify faoed1 with a bit asociated with the tri. This prevents the 
	    // collapse of certain intersections...
	    
            assert(i1/3 <= 7);
	    faoed1 |= ((i1/3)<<28);
	    
	    // be sure face references are different...
	    
	    assert(faoed0 != faoed1);

	    // by convention, use "j" for face 1-related 3-indexing...
	    
	    int bbmin1[3]; FOR_J3 bbmin1[j] = min(x_no_int[nootr1[0]][j],min(x_no_int[nootr1[1]][j],x_no_int[nootr1[2]][j]));
	    int bbmax1[3]; FOR_J3 bbmax1[j] = max(x_no_int[nootr1[0]][j],max(x_no_int[nootr1[1]][j],x_no_int[nootr1[2]][j]));
	    
	    // build a local intersectionVec for every tri pair. This makes the 
	    // segment addition logic as simple as possible...
            
            //++debug_count;
            if (debug) cout << "\n\nface ifa0,i0: " << ifa0 << " " << i0 << " ifa1,i1: " << ifa1 << " " << i1 << endl; //" debug_count: " << debug_count << endl;
            
            if (debug) {
              //if (debug_count == debug_count_check) {
              
              /*
	      cout << "WRITING faces for debug_count_check: " << debug_count_check << 
                ", edges in ifa0: " << edofa_i0[ifa0+1]-edofa_i0[ifa0] << 
                " ifa1: " << edofa_i1[ifa1+1]-edofa_i1[ifa1] << endl;
              */

	      FILE * fp = fopen("ifa0.dat","w");
	      int ino_first = -1;
	      for (int fe0 = edofa_i0[ifa0]; fe0 != edofa_i0[ifa0+1]; ++fe0) {
                if (faceEdgeVec0[fe0].second >= 0) {
                  const int ino = nooed[faceEdgeVec0[fe0].second][0];
                  fprintf(fp,"%18.15le %18.15le %18.15le\n",x_no[ino][0],x_no[ino][1],x_no[ino][2]);
                  if (ino_first == -1)
                    ino_first = ino;
                }
                else {
                  const int ino = nooed[-faceEdgeVec0[fe0].second-1][1];
                  fprintf(fp,"%18.15le %18.15le %18.15le\n",x_no[ino][0],x_no[ino][1],x_no[ino][2]);
                  if (ino_first == -1)
                    ino_first = ino;
                }
	      }
	      // rewrite the first node to close the face loop...
	      assert(ino_first != -1);
	      fprintf(fp,"%18.15le %18.15le %18.15le\n",x_no[ino_first][0],x_no[ino_first][1],x_no[ino_first][2]);
	      fclose(fp);
	      
	      fp = fopen("tri0.dat","w");
	      FOR_I3 {
                const int ino = nootr0[i];
                fprintf(fp,"%18.15le %18.15le %18.15le\n",x_no[ino][0],x_no[ino][1],x_no[ino][2]);
	      }
	      int ino = nootr0[0];
	      fprintf(fp,"%18.15le %18.15le %18.15le\n",x_no[ino][0],x_no[ino][1],x_no[ino][2]);
	      fclose(fp);
	      
	      fp = fopen("ifa1.dat","w");
	      ino_first = -1;
	      for (int fe1 = edofa_i1[ifa1]; fe1 != edofa_i1[ifa1+1]; ++fe1) {
                if (faceEdgeVec1[fe1].second >= 0) {
                  const int ino = nooed[faceEdgeVec1[fe1].second][0];
                  fprintf(fp,"%18.15le %18.15le %18.15le\n",x_no[ino][0],x_no[ino][1],x_no[ino][2]);
                  if (ino_first == -1)
                    ino_first = ino;
                }
                else {
                  const int ino = nooed[-faceEdgeVec1[fe1].second-1][1];
                  fprintf(fp,"%18.15le %18.15le %18.15le\n",x_no[ino][0],x_no[ino][1],x_no[ino][2]);
                  if (ino_first == -1)
                    ino_first = ino;
                }
	      }
	      assert(ino_first != -1);
	      fprintf(fp,"%18.15le %18.15le %18.15le\n",x_no[ino_first][0],x_no[ino_first][1],x_no[ino_first][2]);
	      fclose(fp);
              
	      fp = fopen("tri1.dat","w");
	      FOR_J3 {
                const int ino = nootr1[j];
                fprintf(fp,"%18.15le %18.15le %18.15le\n",x_no[ino][0],x_no[ino][1],x_no[ino][2]);
	      }
	      ino = nootr1[0];
	      fprintf(fp,"%18.15le %18.15le %18.15le\n",x_no[ino][0],x_no[ino][1],x_no[ino][2]);
	      fclose(fp);
	      
	      cout << "files written: ifa0.dat and ifa1.dat, tri0.dat and tri1.dat..." << endl;
              
            }
	    
	    assert(localIntersectionVec.empty());

	    // ==================================================
	    // start by looping on ifa0/i0's 3 edges...
	    // ==================================================
            FOR_I3 { 
              
              //cout << "ifa0/i0 LOOP i: " << i << endl;
              
	      const int ied0 = max(edotr0[i],-edotr0[i]-1); assert(ied0 >= 0); // this edge index already accounts for sign
	      // -----------------------------------------------------------------------------------
	      // if this edge is a cut-cube edge (i.e. on the cut cube surface), then we can use
	      // 2D intersection routines...
	      // -----------------------------------------------------------------------------------
	      if (faoed[ied0][0] < 0) {
		// this is a boundary edge. See if it intersects with any other boundary edges...
		//assert(faceEdgeVec0[fe0].second < 0); // should have been a negative edge
		assert((faoed[ied0][0] <= -2)&&(faoed[ied0][0] >= -7));
		// look for an edge associated with ifa1 that sits on the same cube face...
		FOR_J3 {
		  if ((edotr1[j] < 0)&&(faoed[-edotr1[j]-1][0] == faoed[ied0][0])) {
		    const int ied1 = -edotr1[j]-1;
		    int8 A00,A01,A10,A11;
		    switch (faoed[ied0][0]) {
		    case X0_FACE:
		    case X1_FACE:
		      // X's should be the same!
		      assert(x_no_int[nooed[ied0][0]][0] == x_no_int[nooed[ied0][1]][0]);
		      assert(x_no_int[nooed[ied1][0]][0] == x_no_int[nooed[ied1][1]][0]);
		      assert(x_no_int[nooed[ied1][0]][0] == x_no_int[nooed[ied0][0]][0]);
		      // A00 means ied0 area on the 0-node side...
		      A00 = int8(x_no_int[nooed[ied1][0]][1]-x_no_int[nooed[ied0][0]][1])*int8(x_no_int[nooed[ied1][1]][2]-x_no_int[nooed[ied0][0]][2]) - 
			int8(x_no_int[nooed[ied1][0]][2]-x_no_int[nooed[ied0][0]][2])*int8(x_no_int[nooed[ied1][1]][1]-x_no_int[nooed[ied0][0]][1]);
		      A01 = int8(x_no_int[nooed[ied1][1]][1]-x_no_int[nooed[ied0][1]][1])*int8(x_no_int[nooed[ied1][0]][2]-x_no_int[nooed[ied0][1]][2]) - 
			int8(x_no_int[nooed[ied1][1]][2]-x_no_int[nooed[ied0][1]][2])*int8(x_no_int[nooed[ied1][0]][1]-x_no_int[nooed[ied0][1]][1]);
		      A10 = int8(x_no_int[nooed[ied0][0]][1]-x_no_int[nooed[ied1][0]][1])*int8(x_no_int[nooed[ied0][1]][2]-x_no_int[nooed[ied1][0]][2]) - 
			int8(x_no_int[nooed[ied0][0]][2]-x_no_int[nooed[ied1][0]][2])*int8(x_no_int[nooed[ied0][1]][1]-x_no_int[nooed[ied1][0]][1]);
		      A11 = int8(x_no_int[nooed[ied0][1]][1]-x_no_int[nooed[ied1][1]][1])*int8(x_no_int[nooed[ied0][0]][2]-x_no_int[nooed[ied1][1]][2]) - 
			int8(x_no_int[nooed[ied0][1]][2]-x_no_int[nooed[ied1][1]][2])*int8(x_no_int[nooed[ied0][0]][1]-x_no_int[nooed[ied1][1]][1]);
		      break;
		    case Y0_FACE:
		    case Y1_FACE:
		      // Y's should be the same!
		      assert(x_no_int[nooed[ied0][0]][1] == x_no_int[nooed[ied0][1]][1]);
		      assert(x_no_int[nooed[ied1][0]][1] == x_no_int[nooed[ied1][1]][1]);
		      assert(x_no_int[nooed[ied1][0]][1] == x_no_int[nooed[ied0][0]][1]);
		      A00 = int8(x_no_int[nooed[ied1][0]][2]-x_no_int[nooed[ied0][0]][2])*int8(x_no_int[nooed[ied1][1]][0]-x_no_int[nooed[ied0][0]][0]) - 
			int8(x_no_int[nooed[ied1][0]][0]-x_no_int[nooed[ied0][0]][0])*int8(x_no_int[nooed[ied1][1]][2]-x_no_int[nooed[ied0][0]][2]);
		      A01 = int8(x_no_int[nooed[ied1][1]][2]-x_no_int[nooed[ied0][1]][2])*int8(x_no_int[nooed[ied1][0]][0]-x_no_int[nooed[ied0][1]][0]) - 
			int8(x_no_int[nooed[ied1][1]][0]-x_no_int[nooed[ied0][1]][0])*int8(x_no_int[nooed[ied1][0]][2]-x_no_int[nooed[ied0][1]][2]);
		      A10 = int8(x_no_int[nooed[ied0][0]][2]-x_no_int[nooed[ied1][0]][2])*int8(x_no_int[nooed[ied0][1]][0]-x_no_int[nooed[ied1][0]][0]) - 
			int8(x_no_int[nooed[ied0][0]][0]-x_no_int[nooed[ied1][0]][0])*int8(x_no_int[nooed[ied0][1]][2]-x_no_int[nooed[ied1][0]][2]);
		      A11 = int8(x_no_int[nooed[ied0][1]][2]-x_no_int[nooed[ied1][1]][2])*int8(x_no_int[nooed[ied0][0]][0]-x_no_int[nooed[ied1][1]][0]) - 
			int8(x_no_int[nooed[ied0][1]][0]-x_no_int[nooed[ied1][1]][0])*int8(x_no_int[nooed[ied0][0]][2]-x_no_int[nooed[ied1][1]][2]);
		      break;
		    case Z0_FACE:
		    case Z1_FACE:
		      // Z's should be the same!
		      assert(x_no_int[nooed[ied0][0]][2] == x_no_int[nooed[ied0][1]][2]);
		      assert(x_no_int[nooed[ied1][0]][2] == x_no_int[nooed[ied1][1]][2]);
		      assert(x_no_int[nooed[ied1][0]][2] == x_no_int[nooed[ied0][0]][2]);
		      A00 = int8(x_no_int[nooed[ied1][0]][0]-x_no_int[nooed[ied0][0]][0])*int8(x_no_int[nooed[ied1][1]][1]-x_no_int[nooed[ied0][0]][1]) - 
			int8(x_no_int[nooed[ied1][0]][1]-x_no_int[nooed[ied0][0]][1])*int8(x_no_int[nooed[ied1][1]][0]-x_no_int[nooed[ied0][0]][0]);
		      A01 = int8(x_no_int[nooed[ied1][1]][0]-x_no_int[nooed[ied0][1]][0])*int8(x_no_int[nooed[ied1][0]][1]-x_no_int[nooed[ied0][1]][1]) - 
			int8(x_no_int[nooed[ied1][1]][1]-x_no_int[nooed[ied0][1]][1])*int8(x_no_int[nooed[ied1][0]][0]-x_no_int[nooed[ied0][1]][0]);
		      A10 = int8(x_no_int[nooed[ied0][0]][0]-x_no_int[nooed[ied1][0]][0])*int8(x_no_int[nooed[ied0][1]][1]-x_no_int[nooed[ied1][0]][1]) - 
			int8(x_no_int[nooed[ied0][0]][1]-x_no_int[nooed[ied1][0]][1])*int8(x_no_int[nooed[ied0][1]][0]-x_no_int[nooed[ied1][0]][0]);
		      A11 = int8(x_no_int[nooed[ied0][1]][0]-x_no_int[nooed[ied1][1]][0])*int8(x_no_int[nooed[ied0][0]][1]-x_no_int[nooed[ied1][1]][1]) - 
			int8(x_no_int[nooed[ied0][1]][1]-x_no_int[nooed[ied1][1]][1])*int8(x_no_int[nooed[ied0][0]][0]-x_no_int[nooed[ied1][1]][0]);
		      break;
		    default:
		      assert(0);
		    }
		    //cout << "A00: " << A00 << " A01: " << A01 << " A00+A01: " << A00+A01 << " A10: " << A10 << " A11: " << A11 << " A10+A11: " << A10+A11 << endl;
		    if ( ((A00 >= 0)&&(A01 >= 0)) || ((A00 <= 0)&&(A01 <= 0)) ) {
		      if ( ((A10 >= 0)&&(A11 >= 0)) || ((A10 <= 0)&&(A11 <= 0)) ) {
			// make sure one of either is non-zero...
			assert((A00 != 0)||(A01 != 0));
			assert((A10 != 0)||(A11 != 0));
			/*
			  double wgt = double(A00)/double(A00+A01);
			  double xi[3]; FOR_I3 xi[i] = wgt*x_no[nooed[ied0][1]][i] + (1.0-wgt)*x_no[nooed[ied0][0]][i];
			  double wgt_check = double(A10)/double(A10+A11);
			  double xi_check[3]; FOR_I3 xi_check[i] = wgt_check*x_no[nooed[ied1][1]][i] + (1.0-wgt_check)*x_no[nooed[ied1][0]][i];
			  cout << "check distance: should be small): " << DIST(xi,xi_check) << endl;
			  FILE * fp = fopen("xi.dat","w");
			  fprintf(fp,"%18.15le %18.15le %18.15le\n",xi[0],xi[1],xi[2]);
			  fclose(fp);
			*/
			if (A00 == 0) {
			  // we are at node 0 of ied0...
			  if (A10 == 0) {
			    // we are at node 0 of ied1...
			    localIntersectionVec.push_back(IntersectionData(NODE_NODE_INTERSECTION,nooed[ied0][0],nooed[ied1][0]));
			  }
			  else if (A11 == 0) {
			    // we are at node 1 of ied1...
			    localIntersectionVec.push_back(IntersectionData(NODE_NODE_INTERSECTION,nooed[ied0][0],nooed[ied1][1]));
			  }
			  else {
			    // we are somewhere along ied1...
			    localIntersectionVec.push_back(IntersectionData(NODE_EDGE_INTERSECTION,nooed[ied0][0],ied1,double(A10)/double(A10+A11)));
			  }
			}
			else if (A01 == 0) {
			  // we are at node 1 of ied0...
			  if (A10 == 0) {
			    // we are at node 0 of ied1...
			    localIntersectionVec.push_back(IntersectionData(NODE_NODE_INTERSECTION,nooed[ied0][1],nooed[ied1][0]));
			  }
			  else if (A11 == 0) {
			    // we are at node 1 of ied1...
			    localIntersectionVec.push_back(IntersectionData(NODE_NODE_INTERSECTION,nooed[ied0][1],nooed[ied1][1]));
			  }
			  else {
			    // we are somewhere along ied1...
			    localIntersectionVec.push_back(IntersectionData(NODE_EDGE_INTERSECTION,nooed[ied0][1],ied1,double(A10)/double(A10+A11)));
			  }
			}
			else {
			  // we are somewhere along ied0...
			  if (A10 == 0) {
			    // we are at node 0 of ied1...
			    localIntersectionVec.push_back(IntersectionData(NODE_EDGE_INTERSECTION,nooed[ied1][0],ied0,double(A00)/double(A00+A01)));
			  }
			  else if (A11 == 0) {
			    // we are at node 1 of ied1...
			    localIntersectionVec.push_back(IntersectionData(NODE_EDGE_INTERSECTION,nooed[ied1][1],ied0,double(A00)/double(A00+A01)));
			  }
			  else {
			    // we are somewhere along ied1...
			    localIntersectionVec.push_back(IntersectionData(EDGE_EDGE_INTERSECTION,ied0,ied1,double(A00)/double(A00+A01),double(A10)/double(A10+A11)));
			  }
			}
		      }
		    }
		  }
		}
	      }
	      // -----------------------------------------------------------------------------------
	      // otherwise, look for 3d intersection if the edge bbox overlaps/touches the bbox of tri1...
	      // -----------------------------------------------------------------------------------
	      else if ((min(x_no_int[nooed[ied0][0]][0],x_no_int[nooed[ied0][1]][0]) <= bbmax1[0])&&
		       (max(x_no_int[nooed[ied0][0]][0],x_no_int[nooed[ied0][1]][0]) >= bbmin1[0])&&
		       (min(x_no_int[nooed[ied0][0]][1],x_no_int[nooed[ied0][1]][1]) <= bbmax1[1])&&
		       (max(x_no_int[nooed[ied0][0]][1],x_no_int[nooed[ied0][1]][1]) >= bbmin1[1])&&
		       (min(x_no_int[nooed[ied0][0]][2],x_no_int[nooed[ied0][1]][2]) <= bbmax1[2])&&
		       (max(x_no_int[nooed[ied0][0]][2],x_no_int[nooed[ied0][1]][2]) >= bbmin1[2])) {
		// ied0 against ifa1's tri i1 in 3D...
		// here we compute the tet volumes so that, when the edge is on either side, the volumes have the same sign...
		const int8 tet_vol0 = SIGNED_TET_VOLUME_6_INT8(x_no_int[nooed[ied0][0]],x_no_int[nootr1[0]],x_no_int[nootr1[1]],x_no_int[nootr1[2]]);
		const int8 tet_vol1 = SIGNED_TET_VOLUME_6_INT8(x_no_int[nootr1[0]],x_no_int[nooed[ied0][1]],x_no_int[nootr1[1]],x_no_int[nootr1[2]]);
		if (tet_vol0+tet_vol1 == 0) {
		  // we should be able to skip this case, unless both are zero, meaning that a 2D intersection
		  // is possible. 
		  if (debug) {
		    if ((tet_vol0 == 0)&&(tet_vol1 == 0)) {
		      cout << "XXXXXXXXXXXX skipping possible 2D intersection: L1606: ied0: " << ied0 << " ifa1: " << ifa1 << " i1: " << i1 << endl;
		      getchar();
		      /*
			FILE * fp = fopen("face.dat","w");
			for (int fe1 = edofa_i1[ifa1]; fe1 != edofa_i1[ifa1+1]; ++fe1) {
			const int ied = max(faceEdgeVec1[fe1].second,-faceEdgeVec1[fe1].second-1);
			addEdgeToFp(fp,ied);
			}
			fclose(fp);
			fp = fopen("edge.dat","w");
			addEdgeToFp(fp,ied0);
			fclose(fp);
			cout << "see face_edge.lay..." << endl;
			getchar();
		      */
		    }
		  }
		}
		else if ((tet_vol0 >= 0)&&(tet_vol1 >= 0)) {
		  // this is a potenially a positive crossing of the edge...
		  const int8 tet_vol_no0 = SIGNED_TET_VOLUME_6_INT8(x_no_int[nooed[ied0][0]],x_no_int[nooed[ied0][1]],x_no_int[nootr1[1]],x_no_int[nootr1[2]]);
		  if ((tet_vol_no0 >= 0)&&(tet_vol_no0 <= tet_vol0+tet_vol1)) {
		    const int8 tet_vol_no1 = SIGNED_TET_VOLUME_6_INT8(x_no_int[nooed[ied0][0]],x_no_int[nooed[ied0][1]],x_no_int[nootr1[2]],x_no_int[nootr1[0]]);
		    if ((tet_vol_no1 >= 0)&&(tet_vol_no0+tet_vol_no1 <= tet_vol0+tet_vol1)) {
		      const int8 tet_vol_no2 = SIGNED_TET_VOLUME_6_INT8(x_no_int[nooed[ied0][0]],x_no_int[nooed[ied0][1]],x_no_int[nootr1[0]],x_no_int[nootr1[1]]);
                      if (tet_vol_no2 >= 0) {
			assert(tet_vol_no0+tet_vol_no1+tet_vol_no2 == tet_vol0+tet_vol1);
                        addEdgeTriIntersectionForDoit(localIntersectionVec,ied0,nootr1,edotr1,faoed1,tet_vol0,tet_vol1,tet_vol_no0,tet_vol_no1,tet_vol_no2);
                      }
                    }
		  }
		}
		else if ((tet_vol0 <= 0)&&(tet_vol1 <= 0)) {
		  // this is a potenially a positive crossing of the edge...
		  const int8 tet_vol_no0 = SIGNED_TET_VOLUME_6_INT8(x_no_int[nooed[ied0][0]],x_no_int[nooed[ied0][1]],x_no_int[nootr1[1]],x_no_int[nootr1[2]]);
		  if ((tet_vol_no0 <= 0)&&(tet_vol_no0 >= tet_vol0+tet_vol1)) {
		    const int8 tet_vol_no1 = SIGNED_TET_VOLUME_6_INT8(x_no_int[nooed[ied0][0]],x_no_int[nooed[ied0][1]],x_no_int[nootr1[2]],x_no_int[nootr1[0]]);
		    if ((tet_vol_no1 <= 0)&&(tet_vol_no0+tet_vol_no1 >= tet_vol0+tet_vol1)) {
		      const int8 tet_vol_no2 = SIGNED_TET_VOLUME_6_INT8(x_no_int[nooed[ied0][0]],x_no_int[nooed[ied0][1]],x_no_int[nootr1[0]],x_no_int[nootr1[1]]);
		      if (tet_vol_no2 <= 0) {
			assert(tet_vol_no0+tet_vol_no1+tet_vol_no2 == tet_vol0+tet_vol1);
                        addEdgeTriIntersectionForDoit(localIntersectionVec,ied0,nootr1,edotr1,faoed1,tet_vol0,tet_vol1,tet_vol_no0,tet_vol_no1,tet_vol_no2);
                      }
                    }
                  }
                }
	      }
	    }
            
	    // ==================================================
	    // second: edges of ifa1/i1 tri against ifa0/i0 tri...
	    // ==================================================
	    FOR_J3 {
	      const int ied1 = max(edotr1[j],-edotr1[j]-1); assert(ied1 >= 0); // this edge index already accounts for sign
	      // no need to reconsider the 2D case here. Only need to handle once (above)...
	      if ((faoed[ied1][0] >= 0)&&
		  (min(x_no_int[nooed[ied1][0]][0],x_no_int[nooed[ied1][1]][0]) <= bbmax0[0])&&
		  (max(x_no_int[nooed[ied1][0]][0],x_no_int[nooed[ied1][1]][0]) >= bbmin0[0])&&
		  (min(x_no_int[nooed[ied1][0]][1],x_no_int[nooed[ied1][1]][1]) <= bbmax0[1])&&
		  (max(x_no_int[nooed[ied1][0]][1],x_no_int[nooed[ied1][1]][1]) >= bbmin0[1])&&
		  (min(x_no_int[nooed[ied1][0]][2],x_no_int[nooed[ied1][1]][2]) <= bbmax0[2])&&
		  (max(x_no_int[nooed[ied1][0]][2],x_no_int[nooed[ied1][1]][2]) >= bbmin0[2])) {
		// ied1 against ifa0's tri i0 in 3D...
		// here we compute the tet volumes so that, when the edge is on either side, the volumes have the same sign...
		const int8 tet_vol0 = SIGNED_TET_VOLUME_6_INT8(x_no_int[nooed[ied1][0]],x_no_int[nootr0[0]],x_no_int[nootr0[1]],x_no_int[nootr0[2]]);
		const int8 tet_vol1 = SIGNED_TET_VOLUME_6_INT8(x_no_int[nootr0[0]],x_no_int[nooed[ied1][1]],x_no_int[nootr0[1]],x_no_int[nootr0[2]]);
		if (tet_vol0+tet_vol1 == 0) {
		  if (debug) {
		    if ((tet_vol0 == 0)&&(tet_vol1 == 0)) {
		      cout << "XXXXXXXXXXXX skipping possible 2D intersection: L1850: ied1: " << ied1 << " ifa0: " << ifa0 << " i0: " << i0 << endl;
		      getchar();
		    }
		  }
		}
		else if ((tet_vol0 >= 0)&&(tet_vol1 >= 0)) {
		  // this is a potenially a positive crossing of the edge...
		  const int8 tet_vol_no0 = SIGNED_TET_VOLUME_6_INT8(x_no_int[nooed[ied1][0]],x_no_int[nooed[ied1][1]],x_no_int[nootr0[1]],x_no_int[nootr0[2]]);
		  if ((tet_vol_no0 >= 0)&&(tet_vol_no0 <= tet_vol0+tet_vol1)) {
		    const int8 tet_vol_no1 = SIGNED_TET_VOLUME_6_INT8(x_no_int[nooed[ied1][0]],x_no_int[nooed[ied1][1]],x_no_int[nootr0[2]],x_no_int[nootr0[0]]);
		    if ((tet_vol_no1 >= 0)&&(tet_vol_no0+tet_vol_no1 <= tet_vol0+tet_vol1)) {
		      const int8 tet_vol_no2 = SIGNED_TET_VOLUME_6_INT8(x_no_int[nooed[ied1][0]],x_no_int[nooed[ied1][1]],x_no_int[nootr0[0]],x_no_int[nootr0[1]]);
		      if (tet_vol_no2 >= 0) {
			assert(tet_vol_no0+tet_vol_no1+tet_vol_no2 == tet_vol0+tet_vol1);
                        // use the same routine exploiting symmetry...
                        addEdgeTriIntersectionForDoit(localIntersectionVec,ied1,nootr0,edotr0,faoed0,tet_vol0,tet_vol1,tet_vol_no0,tet_vol_no1,tet_vol_no2);
                      }
                    }
		  }
		}
		else if ((tet_vol0 <= 0)&&(tet_vol1 <= 0)) {
		  // this is a potenially a positive crossing of the edge...
		  const int8 tet_vol_no0 = SIGNED_TET_VOLUME_6_INT8(x_no_int[nooed[ied1][0]],x_no_int[nooed[ied1][1]],x_no_int[nootr0[1]],x_no_int[nootr0[2]]);
		  if ((tet_vol_no0 <= 0)&&(tet_vol_no0 >= tet_vol0+tet_vol1)) {
		    const int8 tet_vol_no1 = SIGNED_TET_VOLUME_6_INT8(x_no_int[nooed[ied1][0]],x_no_int[nooed[ied1][1]],x_no_int[nootr0[2]],x_no_int[nootr0[0]]);
		    if ((tet_vol_no1 <= 0)&&(tet_vol_no0+tet_vol_no1 >= tet_vol0+tet_vol1)) {
		      const int8 tet_vol_no2 = SIGNED_TET_VOLUME_6_INT8(x_no_int[nooed[ied1][0]],x_no_int[nooed[ied1][1]],x_no_int[nootr0[0]],x_no_int[nootr0[1]]);
		      if (tet_vol_no2 <= 0) {
			assert(tet_vol_no0+tet_vol_no1+tet_vol_no2 == tet_vol0+tet_vol1);
                        addEdgeTriIntersectionForDoit(localIntersectionVec,ied1,nootr0,edotr0,faoed0,tet_vol0,tet_vol1,tet_vol_no0,tet_vol_no1,tet_vol_no2);
                      }
		    }
		  }
		}
	      }
	    }
            
	    if (debug) cout << "ifa0,ifa1: " << ifa0 << " " << ifa1 << " localIntersectionVec.size() : " << localIntersectionVec.size() << endl;
	    
	    if (!localIntersectionVec.empty()) {
              
	      if (!(localIntersectionVec.size() >= 2)) {
		cout << "unexpected localIntersectionVec.size(): " << localIntersectionVec.size() << endl;
		throw(-1);
	      }

	      sort(localIntersectionVec.begin(),localIntersectionVec.end());
	      vector<pair<double,int> > diVec;
	      if (debug) cout << "unique INTERSECTIONS between faces ifa0,i0: " << ifa0 << " " << i0 << " AND ifa1,i1: " << ifa1 << " " << i1 << endl;
	      for (int ii = 0; ii < localIntersectionVec.size(); ++ii) {
		if ((ii == 0)||(localIntersectionVec[ii] != localIntersectionVec[ii-1])) {
		  if (debug) {
		    cout << " > Adding unique intersection: ";
		    localIntersectionVec[ii].dump();
		  }
		  diVec.push_back(pair<double,int>(0.0,ii));
		}
		else {
		  if (debug) {
		    cout << " > Skipping duplicate intersection: ";
		    localIntersectionVec[ii].dump();
		  }
		}
	      }
	      
	      // must be 1 or 2 unique intersections...
	      if (diVec.size() == 1) {
		// must be a node touching down, or two tris touching at an edge ...
		assert((localIntersectionVec[diVec[0].second].kind == NODE_NODE_INTERSECTION) ||
		       (localIntersectionVec[diVec[0].second].kind == NODE_EDGE_INTERSECTION) ||
		       (localIntersectionVec[diVec[0].second].kind == NODE_FACE_INTERSECTION) ||
		       (localIntersectionVec[diVec[0].second].kind == EDGE_EDGE_INTERSECTION));
		// if its a node-node, then add it...
		// we have to do this because cases connected through a single node can fail the
		// walking algo below unless we include this intersectino to condense the common 
		// node and also ensure it is in (see intersectionNodeNode vec below)...
		if (localIntersectionVec[diVec[0].second].kind == NODE_NODE_INTERSECTION) {
		  intersectionVec.push_back(localIntersectionVec[diVec[0].second]);
		}
	      }
	      else {

		assert(diVec.size() == 2); // now that we are using tris, this needs to be exactly 2
		
		// TODO: fix this to use tri normals eventually...
		
		// add a new edge to the cvd. This edge will use a special treatment for its node
		// locations: they will refer to the new nodes associated with intersectionVec index.
		const int8 normal_tr0[3] = TRI_NORMAL_2_INT8(x_no_int[nootr0[0]],x_no_int[nootr0[1]],x_no_int[nootr0[2]]);
		const int8 normal_tr1[3] = TRI_NORMAL_2_INT8(x_no_int[nootr1[0]],x_no_int[nootr1[1]],x_no_int[nootr1[2]]);
		const double normal_tr0_d[3] = { (double)normal_tr0[0], (double)normal_tr0[1], (double)normal_tr0[2] };
		const double normal_tr1_d[3] = { (double)normal_tr1[0], (double)normal_tr1[1], (double)normal_tr1[2] };
		
		// cross the face normals: this gives a direction along the intersection...
		const double cp[3] = CROSS_PRODUCT(normal_tr0_d,normal_tr1_d);
		
		// note: cp can be zero, but in this case, we add the edge and "hope" that the 
		// small edge check code below solves any orientation issues.

		double xi[3];
		for (int ii = 0; ii < diVec.size(); ++ii) {
		  localIntersectionVec[diVec[ii].second].calcXi(xi,x_no,nooed);
		  assert(diVec[ii].first == 0.0);
		  diVec[ii].first = DOT_PRODUCT(xi,cp); ///cp_mag;
		}
		sort(diVec.begin(),diVec.end());

		// connect the two furthest together here...
		const int ied = new_edge();
		faoed[ied][0] = faoed0 & ~(7<<28);
		faoed[ied][1] = faoed1 & ~(7<<28);
		nooed[ied][0] = nno+intersectionVec.size();
		nooed[ied][1] = nno+intersectionVec.size()+1;

		//if (debug) cout << "added edge: " << ied << " with debug_count: " << debug_count << endl;
		
		if (debug) {
		  cout << "After cp sort, intersections addedin this order: " << endl;
		  localIntersectionVec[diVec.front().second].dump();
		  localIntersectionVec[diVec.back().second].dump();
                  const double cp_mag = MAG(cp);
                  cout << " > segment length in xi: " << (diVec[1].first-diVec[0].first)/cp_mag << endl;
                }

		// and add the two furthest to intersectionVec...
		// now that there are only 2, this does not leave anyone untouched...
		intersectionVec.push_back(localIntersectionVec[diVec.front().second]);
		intersectionVec.push_back(localIntersectionVec[diVec.back().second]);
		
		/*
		  if (debug_count == debug_count_check) {
		  FILE * fp = fopen("xi.dat","w");
		  for (int ii = 0; ii < diVec.size(); ++ii) {
		  localIntersectionVec[diVec[ii].second].calcXi(xi,x_no,nooed);
		  fprintf(fp,"%18.15le %18.15le %18.15le\n",xi[0],xi[1],xi[2]);
		  }
		  fclose(fp);
		  }
		*/
		
		/*
		  FILE * fp = fopen("xi.dat","w");
		  fprintf(fp,"%18.15le %18.15le %18.15le\n",xi0[0],xi0[1],xi0[2]);
		  fprintf(fp,"%18.15le %18.15le %18.15le\n",xi1[0],xi1[1],xi1[2]);
		  fclose(fp);
		
		  fp = fopen("face0.dat","w");
		  for (int fe0 = edofa_i0[ifa0]; fe0 != edofa_i0[ifa0+1]; ++fe0) {
		  const int ied = max(faceEdgeVec0[fe0].second,-faceEdgeVec0[fe0].second-1);
		  addEdgeToFp(fp,ied);
		  }
		  fclose(fp);
		
		  fp = fopen("face1.dat","w");
		  for (int fe1 = edofa_i1[ifa1]; fe1 != edofa_i1[ifa1+1]; ++fe1) {
		  const int ied = max(faceEdgeVec1[fe1].second,-faceEdgeVec1[fe1].second-1);
		  addEdgeToFp(fp,ied);
		  }
		  fclose(fp);
		*/
		  
	      }

	      localIntersectionVec.clear();
	    }

	    /*
	      if (debug_count == debug_count_check) {
	      cout << "how was that?" << endl;
	      getchar();
	      }
	    */
	    	    
	  }
	}
      }

      delete[] fa0_flag; fa0_flag = NULL;
      delete[] fa1_flag; fa1_flag = NULL;
      
      facePairVec.clear();
      faceEdgeVec0.clear();
      faceEdgeVec1.clear();
      edofa_i0.clear();
      edofa_i1.clear();
      
      // store the range of intersection edges...
      intersectionEdgeRange.push_back(pair<int,int>(ned_old,ned));
	
      if (!intersectionVec.empty()) {
	  
	// assume each intersection is associated with a new node...
	for (int ii = 0; ii < intersectionVec.size(); ++ii) {
	  intersectionVec[ii].ino = nno+ii;
	}
	sort(intersectionVec.begin(),intersectionVec.end());

	// reallocate no_flag if necessary...
	if (nno+intersectionVec.size() > nno_max) {
	  nno_max = nno+intersectionVec.size();
          assert(no_flag);
	  delete[] no_flag;
	  no_flag = new int[nno_max];
          assert(x_no_int);
          delete[] x_no_int;
          x_no_int = new int[nno_max][3];
	}
	  
	for (int ino = 0; ino < nno; ++ino)
	  no_flag[ino] = ino;
	for (int ino = nno; ino < nno+intersectionVec.size(); ++ino)
	  no_flag[ino] = -1;
    
	// add the nodes...
	assert(cutEdgeDataVec.empty());
	int ii_prev = -1;
	for (int ii = 0; ii < intersectionVec.size(); ++ii) {
	  if ((ii_prev == -1)||
	      (intersectionVec[ii].kind != intersectionVec[ii_prev].kind)||
	      (intersectionVec[ii].idata[0] != intersectionVec[ii_prev].idata[0])||
	      (intersectionVec[ii].idata[1] != intersectionVec[ii_prev].idata[1])) {
	    ii_prev = ii;
	    switch (intersectionVec[ii].kind) {
	    case NODE_NODE_INTERSECTION:
	      {
		// combine nodes in idata[0] and idata[1]...
		int ino0 = no_flag[intersectionVec[ii].idata[0]];
		while (ino0 != no_flag[ino0])
		  ino0 = no_flag[ino0];
		int ino1 = no_flag[intersectionVec[ii].idata[1]];
		while (ino1 != no_flag[ino1])
		  ino1 = no_flag[ino1];
		no_flag[intersectionVec[ii].ino] = no_flag[ino0] = no_flag[ino1] = min(ino0,ino1);
		// and add this to a special vec of node-node intersections...
		intersectionNodeNode.push_back(min(ino0,ino1));
	      }
	      break;
	    case NODE_EDGE_INTERSECTION:
	      {
		// use the node to define the x...
		no_flag[intersectionVec[ii].ino] = intersectionVec[ii].idata[0];
		const int ied = intersectionVec[ii].idata[1];
		const double wgt = intersectionVec[ii].ddata[1];
		cutEdgeDataVec.push_back(pair<pair<int,double>,int>(pair<int,double>(ied,wgt),intersectionVec[ii].idata[0]));
	      }
	      break;
	    case NODE_FACE_INTERSECTION:
	      {
		// use the node to define the x...
		no_flag[intersectionVec[ii].ino] = intersectionVec[ii].idata[0];
	      }
	      break;
	    case EDGE_EDGE_INTERSECTION:
	      {
		// this requires a new node. Use the simple average of the two edge intersections...
		const int ied0 = intersectionVec[ii].idata[0];
		const double wgt0 = intersectionVec[ii].ddata[0];
		const int ied1 = intersectionVec[ii].idata[1];
		const double wgt1 = intersectionVec[ii].ddata[1];
		const int ino = new_node();
		FOR_I3 x_no[ino][i] = 0.5*( wgt0*x_no[nooed[ied0][1]][i] + (1.0-wgt0)*x_no[nooed[ied0][0]][i] +
					    wgt1*x_no[nooed[ied1][1]][i] + (1.0-wgt1)*x_no[nooed[ied1][0]][i] );
		no_flag[intersectionVec[ii].ino] = ino;
		cutEdgeDataVec.push_back(pair<pair<int,double>,int>(pair<int,double>(ied0,wgt0),ino));
		cutEdgeDataVec.push_back(pair<pair<int,double>,int>(pair<int,double>(ied1,wgt1),ino));
	      }
	      break;
	    case EDGE_FACE_INTERSECTION:
	      { 
		const int ied = intersectionVec[ii].idata[0];
		const double wgt = intersectionVec[ii].ddata[0];
		const int ino = new_node();
		FOR_I3 x_no[ino][i] = wgt*x_no[nooed[ied][1]][i] + (1.0-wgt)*x_no[nooed[ied][0]][i];
		no_flag[intersectionVec[ii].ino] = ino;
		cutEdgeDataVec.push_back(pair<pair<int,double>,int>(pair<int,double>(ied,wgt),ino));
	      }
	      break;
	    default:
	      assert(0);
	    }
	  }
	  else {
	    assert(ii_prev != -1);
	    no_flag[intersectionVec[ii].ino] = no_flag[intersectionVec[ii_prev].ino];
	  }
	}
	intersectionVec.clear();

	// now update nooed for ALL edges to the no_flag value...
	for (int ied = 0; ied < ned; ++ied) {
	  FOR_I2 {
	    const int ino = nooed[ied][i];
	    assert((no_flag[ino] >= 0)&&(no_flag[ino] < nno));
	    nooed[ied][i] = no_flag[ino];
	  }
	}
	  
	// and cut the edges that have intersections...
	sort(cutEdgeDataVec.begin(),cutEdgeDataVec.end());
	int ied_current = -1;
	int ino_current = -1;
	for (int ii = 0; ii < cutEdgeDataVec.size(); ++ii) {
	  const int ied = cutEdgeDataVec[ii].first.first;
	  if (ied != ied_current) {
	    // this is the first of a new edge. Complete the old edge first...
	    if (ied_current >= 0) {
	      assert(ino_current >= 0);
	      nooed[ied_current][0] = ino_current;
	    }
	    ied_current = ied;
	    ino_current = nooed[ied][0];
	  }
	  //const double s = cutEdgeDataVec[ii].first.second;
	  const int ino_new = cutEdgeDataVec[ii].second;
	  assert(ino_new != ino_current);
	  const int ied_new = new_edge();
	  nooed[ied_new][0] = ino_current;
	  nooed[ied_new][1] = ino_new;
	  faoed[ied_new][0] = faoed[ied_current][0];
	  faoed[ied_new][1] = faoed[ied_current][1];
	  ino_current = ino_new;
	}
	cutEdgeDataVec.clear();
	// complete last edge (the original edge)...
	if (ied_current >= 0) {
	  assert(ino_current >= 0);
	  nooed[ied_current][0] = ino_current;
	}

      }

    }
    
    delete[] x_no_int; x_no_int = NULL;
    
    if (debug) {
      writeTecplot(-3);

      cout << "full geom with intersections in tecplot -3" << endl;

      FILE * fp = fopen("edges.dat","w");
      for (int ii = 0; ii < intersectionEdgeRange.size(); ++ii) {
	for (int ied = intersectionEdgeRange[ii].first; ied != intersectionEdgeRange[ii].second; ++ied) {
	  cout << " > intersection edge: " << ii << " length: " << DIST(x_no[nooed[ied][0]],x_no[nooed[ied][1]]) << endl; 
	  addEdgeToFp(fp,ied);
	}
      }
      fclose(fp);
      cout << "edges.dat contains intersection edges" << endl;

      // HACK: look for edge duplicates. I think this can happen...
    
      bool got_matching_edge = false;
      fp = fopen("matching_edges.dat","w");
      for (int ied1 = 0; ied1 < ned; ++ied1) {
	for (int ied2 = ied1+1; ied2 < ned; ++ied2) {
	  if ((nooed[ied2][0] == nooed[ied1][0])&&(nooed[ied2][1] == nooed[ied1][1])) {
	    cout << "EDGE MATCH: " << ied1 << " " << ied2 << endl;
	    addEdgeToFp(fp,ied1);
	    addEdgeToFp(fp,ied2);
	    got_matching_edge = true;
	  }
	  else if ((nooed[ied2][0] == nooed[ied1][1])&&(nooed[ied2][1] == nooed[ied1][0])) {
	    cout << "EDGE MATCH: " << ied1 << " " << ied2 << endl;
	    addEdgeToFp(fp,ied1);
	    addEdgeToFp(fp,ied2);
	    got_matching_edge = true;
	  }
	}
      }
      fclose(fp);
      
      if (got_matching_edge) {
	cout << "matching_edges.dat contains matching intersection edges" << endl;
	getchar();
      }

    }

    //cout << "TAKE A LOOK at -3!" << endl;
    //getchar();
    
    assert(nno <= nno_max);
    FOR_INO no_flag[ino] = 0;

    // first, use the link concept to ensure faoed are properly aligned on all edges...
    
    int loop_or_line_count = 0;
    for (int ii = 0; ii < intersectionEdgeRange.size(); ++ii) {
      // count...
      // here we use the no_flag to store the edge index in the first 15 and second 15 bits...
      for (int ied = intersectionEdgeRange[ii].first; ied != intersectionEdgeRange[ii].second; ++ied) {
	assert((ied >= 0)&&(ied < 32767));
	FOR_I2 {
	  const int ino = nooed[ied][i];
	  if (no_flag[ino] == 0) {
	    no_flag[ino] = (ied+1); // skip 0 and use 1,2,3... so zero can have meaning
	    assert((no_flag[ino]&32767) == (ied+1));
	  }
	  else {
	    assert((no_flag[ino]&32767) != 0);
	    if (!((no_flag[ino]>>15) == 0))
              throw(-1);
	    assert((no_flag[ino]>>15) == 0);
	    no_flag[ino] |= (ied+1)<<15;
	    assert((no_flag[ino]&1073709056) == ((ied+1)<<15));
	  }
	}
      }
      vector<int> singleNodeVec;
      FOR_INO {
	// look for only one edge reference in the first position of the no_flag...
	if ((no_flag[ino] != 0)&&((no_flag[ino]&1073709056) == 0)) {
	  singleNodeVec.push_back(ino);
	}
      }
      // now form line(s) and/or loop(s)...
      while (1) {
	++loop_or_line_count;
	// get a starting node...
	int ino_start = -1;
	int ino_prev;
	int flip_count = 0;
        double flip_length = 0.0;
        double total_length = 0.0;
	bool prev_ied_orientation = true; // true = 0,1 false = 1,0
	bool prev_ied_flipped = false;
	assert(intVec.empty()); // use intVec as an edgeVec
	for (int ii = 0; ii < singleNodeVec.size(); ++ii) {
	  const int ino = singleNodeVec[ii];
	  if (no_flag[ino] != 0) {
	    ino_start = ino;
	    // this should just be a single node...
	    assert((no_flag[ino]&32767) != 0);
	    assert((no_flag[ino]&1073709056) == 0);
	    const int ied = no_flag[ino]-1;
	    intVec.push_back(ied);
            const double this_length = DIST(x_no[nooed[ied][0]],x_no[nooed[ied][1]]);
            total_length += this_length;
            no_flag[ino] = 0;
	    if (nooed[ied][0] == ino_start) {
	      prev_ied_orientation = true;
	      ino_prev = nooed[ied][1];
	    }
	    else {
	      assert(nooed[ied][1] == ino_start);
	      prev_ied_orientation = false;
	      ino_prev = nooed[ied][0];
	    }
            assert(prev_ied_flipped == false);
	    // now remove ied from ino_prev. It could be the first or second...
	    if ((no_flag[ino_prev]>>15) == (ied+1)) {
	      no_flag[ino_prev] &= 32767;
	      assert(no_flag[ino_prev] != 0);
	    }
	    else {
	      assert((no_flag[ino_prev]&32767) == (ied+1));
	      no_flag[ino_prev] = (no_flag[ino_prev]>>15);
	      if (no_flag[ino_prev] == 0) {
		ino_start = -1;
		continue;
	      }
	    }
	    break;
	  }
	}
	// if we did not get one, then check if there is an internal 
	// intersection associated with loops...
	if (ino_start == -1) {
	  FOR_INO {
	    if (no_flag[ino] > 0) {
	      ino_start = ino;
	      // anything left should have both available...
	      assert((no_flag[ino]&32767) != 0); // first 15 bits
	      assert((no_flag[ino]&1073709056) != 0); // top 15 bits
	      // take the top edge to start...
	      const int ied = (no_flag[ino]>>15)-1;
	      intVec.push_back(ied);
              const double this_length = DIST(x_no[nooed[ied][0]],x_no[nooed[ied][1]]);
              total_length += this_length;
	      no_flag[ino] &= 32767;
	      if (nooed[ied][0] == ino_start) {
		prev_ied_orientation = true;
		ino_prev = nooed[ied][1];
	      }
	      else {
		assert(nooed[ied][1] == ino_start);
		prev_ied_orientation = false;
		ino_prev = nooed[ied][0];
	      }
              assert(prev_ied_flipped == false);
	      // now remove ied from ino_prev. It could be the first or second...
	      if ((no_flag[ino_prev]>>15) == (ied+1)) {
		no_flag[ino_prev] &= 32767;
		assert(no_flag[ino_prev] != 0);
	      }
	      else {
		assert((no_flag[ino_prev]&32767) == (ied+1));
		no_flag[ino_prev] = (no_flag[ino_prev]>>15);
		assert(no_flag[ino_prev] != 0);
	      }
	      break;
	    }
	  }
	}
	if (ino_start == -1)
	  break;
	// we have an ino_start as the beginning of a loop. Now 
	// cycle through the loop until we come back to the start, or
	// reach another terminal point...
	while (1) {
	  assert(no_flag[ino_prev] != 0);
	  assert((no_flag[ino_prev]&1073709056) == 0); // top 15 bits
	  const int ied = no_flag[ino_prev]-1;
	  no_flag[ino_prev] = 0;
	  intVec.push_back(ied);
          const double this_length = DIST(x_no[nooed[ied][0]],x_no[nooed[ied][1]]);
          total_length += this_length;
	  if (nooed[ied][0] == ino_prev) {
	    // this edge is aligned with the node direction...
	    if (prev_ied_orientation && prev_ied_flipped) {
	      const int tmp = faoed[ied][0];
	      faoed[ied][0] = faoed[ied][1];
	      faoed[ied][1] = tmp;
	      prev_ied_flipped = true;
	      ++flip_count;
              flip_length += this_length;
            }
	    else if ((!prev_ied_orientation) && (!prev_ied_flipped)) {
	      const int tmp = faoed[ied][0];
	      faoed[ied][0] = faoed[ied][1];
	      faoed[ied][1] = tmp;
	      prev_ied_flipped = true;
	      ++flip_count;
              flip_length += this_length;
	    }
	    else {
	      // the other 2 combinations do not require us to flip
	      prev_ied_flipped = false;
	    }
	    prev_ied_orientation = true;
	    ino_prev = nooed[ied][1];
	  }
	  else {
	    assert(nooed[ied][1] == ino_prev);
	    // this edge is mis-aligned with the node direction...
	    if ((!prev_ied_orientation) && prev_ied_flipped) {
	      const int tmp = faoed[ied][0];
	      faoed[ied][0] = faoed[ied][1];
	      faoed[ied][1] = tmp;
	      prev_ied_flipped = true;
	      ++flip_count;
              flip_length += this_length;
	    }
	    else if (prev_ied_orientation && (!prev_ied_flipped)) {
	      const int tmp = faoed[ied][0];
	      faoed[ied][0] = faoed[ied][1];
	      faoed[ied][1] = tmp;
	      prev_ied_flipped = true;
	      ++flip_count;
              flip_length += this_length;
	    }
	    else {
	      // the other 2 combinations do not require us to flip
	      prev_ied_flipped = false;
	    }
	    prev_ied_orientation = false;
	    ino_prev = nooed[ied][0];
	  }
	  // now remove ied from ino_prev. It could be the first or second...
	  if ((no_flag[ino_prev]>>15) == (ied+1)) {
	    no_flag[ino_prev] &= 32767;
	    assert(no_flag[ino_prev] != 0);
	  }
	  else {
	    assert((no_flag[ino_prev]&32767) == (ied+1));
	    no_flag[ino_prev] = (no_flag[ino_prev]>>15);
	    if (no_flag[ino_prev] == 0) {
	      break;
	    }
	  }
	}
	// if we flipped more than half the edges, then flip them all again...
	if (debug) cout << " > considering flip of loop or line: " << loop_or_line_count << 
                     " flip_count: " << flip_count << " intVec.size(): " << intVec.size() << 
                     " flip_length/total_length: " << flip_length/total_length << endl;
	//if ((flip_count*2 > intVec.size())||(loop_or_line_count == 1)) {
	if (flip_length*2.0 >= total_length) {
	  if (flip_count*2 == intVec.size()) cout << "WARNING: flip_count exactly half of intVec.size(): " << intVec.size() << endl;
	  for (int ii = 0; ii < intVec.size(); ++ii) {
	    const int ied = intVec[ii];
	    const int tmp = faoed[ied][0];
	    faoed[ied][0] = faoed[ied][1];
	    faoed[ied][1] = tmp;
	  }
	}
	intVec.clear();
      }
    }
    
    // ========================================================
    // with all edges in the correct order, process
    // in or out...
    // ========================================================

    // this is a check: above should have used each edge at each node and 
    // returned no_flag to 0...
    FOR_INO assert(no_flag[ino] == 0);

    // node-node intersections are in. These are treated specially because a node-node
    // intersection does not necessarily imply an edge intersection. i.e. surfaces can 
    // intersect at a point...
    for (int ii = 0; ii < intersectionNodeNode.size(); ++ii) {
      const int ino = intersectionNodeNode[ii];
      no_flag[ino] = -2;
    }
    
    // here we use the face orientation and the orientation of the intersection faces 
    // to determine some out points. Basically put the intersection edges into a set
    // in both directions around the face, then try and add edges that touch these
    // nodes along the intersections. If you cannot add the edge, then it must be
    // inconsistent with the face orientation set by the intersections, and thus out...
    
    // hack sort based on edge length, and introduce edges largest to smallest. 
    // smallest edges tend to have problems. Eventually when the combinatorial logic
    // is worked out, this can be removed.

    int * ed_flag = new int[ned];
    FOR_IED ed_flag[ied] = 0;
    
    // used to sort this: now we don't have to...
    /*
      vector<pair<double,int> > diVec;
      for (int ii = 0; ii < intersectionEdgeRange.size(); ++ii) {
      for (int ied = intersectionEdgeRange[ii].first; ied != intersectionEdgeRange[ii].second; ++ied) {
      diVec.push_back(pair<double,int>(DIST2(x_no[nooed[ied][0]],x_no[nooed[ied][1]]),ied));
      }
      }
      sort(diVec.begin(),diVec.end());
    */

    set<pair<int,int> > forwardFaNoSet;
    set<pair<int,int> > backwardFaNoSet;
    for (int ii = 0; ii < intersectionEdgeRange.size(); ++ii) {
      for (int ied = intersectionEdgeRange[ii].first; ied != intersectionEdgeRange[ii].second; ++ied) {
	ed_flag[ied] = 1; // all intersections are in...
	no_flag[nooed[ied][0]] = -2;
	no_flag[nooed[ied][1]] = -2;
	bool f0 = (forwardFaNoSet.find(pair<int,int>(faoed[ied][0],nooed[ied][0])) == forwardFaNoSet.end());
	bool f1 = (forwardFaNoSet.find(pair<int,int>(faoed[ied][1],nooed[ied][1])) == forwardFaNoSet.end());
	bool b0 = (backwardFaNoSet.find(pair<int,int>(faoed[ied][0],nooed[ied][1])) == backwardFaNoSet.end());
	bool b1 = (backwardFaNoSet.find(pair<int,int>(faoed[ied][1],nooed[ied][0])) == backwardFaNoSet.end());
	if (f0 && f1 && b0 && b1) {
	  // looks good...
	  forwardFaNoSet.insert(pair<int,int>(faoed[ied][0],nooed[ied][0]));
	  forwardFaNoSet.insert(pair<int,int>(faoed[ied][1],nooed[ied][1]));
	  backwardFaNoSet.insert(pair<int,int>(faoed[ied][0],nooed[ied][1]));
	  backwardFaNoSet.insert(pair<int,int>(faoed[ied][1],nooed[ied][0]));
	}
	else {
	  // one or more failures. We must have introduced the edge in the wrong 
	  // orientation above (during the dp check). This should no longer happen unless the loop
	  // logic above if flawed...
	  cout << "Error: intersection edge orientation FAILURE ied: " << ied << " f0: " << f0 << " f1: " << f1 << " b0: " << b0 << " b1: " << b1 << ", trying flip logic..." << endl;
	  throw(-1);
	}
      }
    }
    
    /*
    // old code when combinatorial logic is working...
    for (int ii = 0; ii < intersectionEdgeRange.size(); ++ii) {
    for (int ied = intersectionEdgeRange[ii].first; ied != intersectionEdgeRange[ii].second; ++ied) {
    ed_flag[ied] = 1; // all intersections are in...
    no_flag[nooed[ied][0]] = -2;
    no_flag[nooed[ied][1]] = -2;
    // faoed[0]...
    set<pair<int,int> >::iterator iter = forwardFaNoSet.find(pair<int,int>(faoed[ied][0],nooed[ied][0]));
    if (iter != forwardFaNoSet.end()) {
    }
    forwardFaNoSet.insert(pair<int,int>(faoed[ied][0],nooed[ied][0]));
    iter = backwardFaNoSet.find(pair<int,int>(faoed[ied][0],nooed[ied][1]));
    if (iter != backwardFaNoSet.end()) {
    cout << "Error: L2036" << endl;
    throw(-1);
    }
    backwardFaNoSet.insert(pair<int,int>(faoed[ied][0],nooed[ied][1]));
    // faoed[1]...
    iter = forwardFaNoSet.find(pair<int,int>(faoed[ied][1],nooed[ied][1]));
    if (iter != forwardFaNoSet.end()) {
    cout << "Error: L2043" << endl;
    throw(-1);
    }
    forwardFaNoSet.insert(pair<int,int>(faoed[ied][1],nooed[ied][1]));
    iter = backwardFaNoSet.find(pair<int,int>(faoed[ied][1],nooed[ied][0]));
    if (iter != backwardFaNoSet.end()) {
    cout << "Error: L2049" << endl;
    throw(-1);
    }
    backwardFaNoSet.insert(pair<int,int>(faoed[ied][1],nooed[ied][0]));
    }
    }
    */

    FOR_IED {
      if (ed_flag[ied] == 0) {
	// skip edges that were added to faces to make them tris. These
	// are eliminated below...
	if (faoed[ied][0] == faoed[ied][1])
	  continue;
	if (no_flag[nooed[ied][0]] == -2) {
	  if ( (forwardFaNoSet.find(pair<int,int>(faoed[ied][0],nooed[ied][0])) != forwardFaNoSet.end()) ||
	       (backwardFaNoSet.find(pair<int,int>(faoed[ied][1],nooed[ied][0])) != backwardFaNoSet.end()) ) {
	    ed_flag[ied] = -1;
	    if (no_flag[nooed[ied][1]] >= 0)
	      no_flag[nooed[ied][1]] = -1;
	  }
	}
	else if (no_flag[nooed[ied][0]] >= 0) {
	  // only count the cube edge nodes that touch a given node...
	  ++no_flag[nooed[ied][0]];
	}
	if (no_flag[nooed[ied][1]] == -2) {
	  if ( (forwardFaNoSet.find(pair<int,int>(faoed[ied][1],nooed[ied][1])) != forwardFaNoSet.end()) ||
	       (backwardFaNoSet.find(pair<int,int>(faoed[ied][0],nooed[ied][1])) != backwardFaNoSet.end()) ) {
	    ed_flag[ied] = -1;
	    if (no_flag[nooed[ied][0]] >= 0)
	      no_flag[nooed[ied][0]] = -1;
	  }
	}
	else if (no_flag[nooed[ied][1]] >= 0) {
	  // only count the cube edge nodes that touch a given node...
	  ++no_flag[nooed[ied][1]];
	}
      }
    }
    
    // convert counts where only 2 edges touch a node to out...
    
    FOR_INO {
      if (no_flag[ino] == 2) {
	no_flag[ino] = -1;
      }
      else if (no_flag[ino] > 0) {
	no_flag[ino] = 0;
      }
    }
    
    // HACK: skip counting logic...
    /*
      FOR_INO {
      if (no_flag[ino] > 0) {
      no_flag[ino] = 0;
      }
      }
    */
    
    if (debug) {
      FILE * fp = fopen("nodes0.dat","w");
      FOR_INO {
	fprintf(fp,"%18.15le %18.15le %18.15le %d\n",x_no[ino][0],x_no[ino][1],x_no[ino][2],no_flag[ino]);
      }
      fclose(fp);
      cout << "checkout nodes0.dat" << endl;
      //getchar();
    }

    // now no_flag has:
    // -2 : marks the intersection line
    // -1 : out
    //  0 : TBD
	
    int step = 0;
    int done = 0;
    while (done == 0) {

      ++step;
      if (debug) cout << "step : " << step << endl;
      
      done = 1;
      FOR_IED {
	if (no_flag[nooed[ied][0]] == -1) {
	  ed_flag[ied] = -1;
	  if (no_flag[nooed[ied][1]] == 0) {
	    done = 0;
	    no_flag[nooed[ied][1]] = -1;
	  }
	}
	else if (no_flag[nooed[ied][1]] == -1) {
	  ed_flag[ied] = -1;
	  if (no_flag[nooed[ied][0]] == 0) {
	    done = 0;
	    no_flag[nooed[ied][0]] = -1;
	  }
	}
      }
    }
    
    if (debug) {
      FILE * fp = fopen("nodes.dat","w");
      FOR_INO {
	fprintf(fp,"%18.15le %18.15le %18.15le %d\n",x_no[ino][0],x_no[ino][1],x_no[ino][2],no_flag[ino]);
      }
      fclose(fp);
    }

    // at this point, we keep nodes with no_flag 0 and -2, and discard -1,
    // and keep edges with 0,1, and discard -1
    
    FOR_IED {
      if (faoed[ied][0] == faoed[ied][1])
	continue;
      if (ed_flag[ied] >= 0) {
	// this edge is in! 
	// do some error checking...
	if (faoed[ied][1] == -1) {
	  cout << "ERROR: open edge not cut by other parts." << endl;
	  throw(-1);
	}
	no_flag[nooed[ied][0]] = -3;
	no_flag[nooed[ied][1]] = -3;
      }
    }
    const int nno_old = nno;
    nno = 0;
    for (int ino = 0; ino < nno_old; ++ino) {
      if (no_flag[ino] == -3) {
	const int ino_new = nno++;
	FOR_I3 x_no[ino_new][i] = x_no[ino][i];
	no_flag[ino] = ino_new;
      }
      else {
	no_flag[ino] = -1;
      }
    }
    const int ned_old = ned;
    ned = 0;
    for (int ied = 0; ied < ned_old; ++ied) {
      if ((faoed[ied][0] != faoed[ied][1])&&(ed_flag[ied] >= 0)) {
	const int ied_new = ned++;
	nooed[ied_new][0] = no_flag[nooed[ied][0]]; assert(nooed[ied_new][0] >= 0);
	nooed[ied_new][1] = no_flag[nooed[ied][1]]; assert(nooed[ied_new][1] >= 0);
	faoed[ied_new][0] = faoed[ied][0];
	faoed[ied_new][1] = faoed[ied][1];
      }
    }
    delete[] ed_flag;
    
  }
  
  if (debug) {
    writeTecplot(-4);
    cout << "final geom in tecplot -4" << endl;
    //getchar();
  }
  
  delete[] no_flag;

  //cout << "LOOK AT -4" << endl;
  //getchar();
    
}
