#include "bcprofile/ProfileReader.hpp"

void InflowTurbulenceIG::setStatistics() {
  Param * p = getParam(getName());
  //assert(p->getString() == "INFLOW_TURB"); // could be HOOK
  int iarg = 1;
  while (iarg < p->size()) {
    string token = p->getString(iarg++);
    //    cout<< token <<endl;

    //  first is the profile
    if (token == "PROFILE") {
      token = p->getString(iarg++);
      if (mpi_rank == 0)
        cout << " > PROFILE: " << token << endl;

      //    stats from a hook
      //    PROFILE HOOK
      if (token == "HOOK") {
        this->inflowTurbulenceStatsHook(um_it, Rd_it,Rod_it, Tm_it, pm_it);
      }//"HOOK"

      //    constant properties across the zone given tks
      //    PROFILE CONSTANT_UTKETP <ux> <uy> <uz> <tke> <T> <p>
      //    NOTE: turbulent kinetic enery (tke) is defined as:
      //          tke = 1/2 * (u'u'+v'v'+w'w')
      else if ( token == "CONSTANT_UTKETP" ) {
        double my_u_bc[3];
        my_u_bc[0] = p->getDouble(iarg++);
        my_u_bc[1] = p->getDouble(iarg++);
        my_u_bc[2] = p->getDouble(iarg++);
        double tke = p->getDouble(iarg++);
        double Tm  = p->getDouble(iarg++);
        double pm  = p->getDouble(iarg++);
        if (mpi_rank == 0)
          cout << "  > using constant mean velocity, constant turbulent kinetic energy, and constant T and P: (ux,uy,uz) = (" <<
            my_u_bc[0] << "," << my_u_bc[1] << "," << my_u_bc[2] << "), tke = " << tke <<", (T,p) = (" <<
            Tm << "," << pm << ")" << endl;
        double ms  = 2.0/3.0 * tke;
        for (int ibf =0; ibf < zone_ptr->nbf; ++ibf) {
          FOR_I3 um_it[ibf][i]  = my_u_bc[i];
          FOR_I3 Rd_it[ibf][i]  = ms;
          FOR_I3 Rod_it[ibf][i] = 0.0;
          Tm_it[ibf]            = Tm;
          pm_it[ibf]            = pm;
        }//(ibf < nbf)
      }//"CONTSTANT UTKETP"
      else if (token == "TRIP_CONSTANT_TKE_ADIABATIC") {
        double tke = p->getDouble(iarg++) ;
        if (mpi_rank == 0 ) {
          cout << " > using constant tke = " << tke << endl ;
        }//(root)

        double ms = 2.0/3.0 * tke ;
        for (int ibf =0; ibf < zone_ptr->nbf; ++ibf) {
          FOR_I3 {
            um_it[ibf][i]  = 0.0; // wall conditions ..
            Rd_it[ibf][i]  = ms ;
            Rod_it[ibf][i] = 0.0 ;
          }//FOR_I3
          Tm_it[ibf] = 0.0; // trigger an adiabatic treatment of this trip...
          pm_it[ibf] = 0.0 ; // "                                             "


        }//(ibf < nbf)
        //      TODO fazone->flag |= INFLOW_TURB_T_BIT | INFLOW_TURB_P_BIT ; /// use the local t, p....
      }//"TRIP_CONSTANT_TKE_ADIABATIC"
      //    constant properties across the zone given mean square quantities
      //    PROFILE CONSTANT_UMSTP <ux> <uy> <uz> <u'u'> <v'v'> <w'w'> <T> <p>
      else if ( token == "CONSTANT_UMSTP" ) {
        double my_u_bc[3];
        my_u_bc[0] = p->getDouble(iarg++);
        my_u_bc[1] = p->getDouble(iarg++);
        my_u_bc[2] = p->getDouble(iarg++);
        double my_uu_bc[3];
        my_uu_bc[0] = p->getDouble(iarg++);
        my_uu_bc[1] = p->getDouble(iarg++);
        my_uu_bc[2] = p->getDouble(iarg++);
        double Tm   = p->getDouble(iarg++);
        double pm   = p->getDouble(iarg++);
        if (mpi_rank == 0)
          cout << "  > using constant mean velocity and constant turbulent mean squares: (ux,uy,uz) = (" <<
            my_u_bc[0] << "," << my_u_bc[1] << "," << my_u_bc[2] << "), (uxux,uyuy,uzuz) = (" <<
            my_uu_bc[0] << "," << my_uu_bc[1] << "," << my_uu_bc[2] << "), (T,p) = (" <<
            Tm << "," << pm << ")" << endl;
        for (int ibf =0; ibf < zone_ptr->nbf; ++ibf) {
          FOR_I3 um_it[ibf][i]  = my_u_bc[i];
          FOR_I3 Rd_it[ibf][i]  = my_uu_bc[i];
          FOR_I3 Rod_it[ibf][i] = 0.0;
          Tm_it[ibf]            = Tm;
          pm_it[ibf]            = pm;
        }//(ibf < nbf)
      }//"CONSTANT_UMSTP"

      //    constant properties across the zone given Reynolds stresses
      //    PROFILE CONSTANT_UREYTP <ux> <uy> <uz> <u'u'> <v'v'> <w'w'> <v'w'> <u'w'> <u'v'> <T> <p>
      //    and
      //    Rd  = [<u'u'> <v'v'> <w'w'>]
      //    Rod = [<v'w'> <u'w'> <u'v'>]
      else if ( token == "CONSTANT_UREYTP" ) {
        double my_u_bc[3];
        my_u_bc[0] = p->getDouble(iarg++);
        my_u_bc[1] = p->getDouble(iarg++);
        my_u_bc[2] = p->getDouble(iarg++);
        double my_uu_bc[3];
        my_uu_bc[0] = p->getDouble(iarg++);
        my_uu_bc[1] = p->getDouble(iarg++);
        my_uu_bc[2] = p->getDouble(iarg++);
        double my_uv_bc[3];
        my_uv_bc[0] = p->getDouble(iarg++);
        my_uv_bc[1] = p->getDouble(iarg++);
        my_uv_bc[2] = p->getDouble(iarg++);
        double Tm   = p->getDouble(iarg++);
        double pm   = p->getDouble(iarg++);
        if (mpi_rank == 0)
          cout << "  > using constant mean velocity and constant Reynolds stresses: (ux,uy,uz) = (" <<
            my_u_bc[0] << "," << my_u_bc[1] << "," << my_u_bc[2] << "), (uxux,uyuy,uzuz) = (" <<
            my_uu_bc[0] << "," << my_uu_bc[1] << "," << my_uu_bc[2] << "), (uxuy,uxuz,uyuz) = (" <<
            my_uv_bc[0] << "," << my_uv_bc[1] << "," << my_uv_bc[2] << "), (T,p) = (" <<
            Tm << "," << pm << ")" << endl;
        for (int ibf =0; ibf < zone_ptr->nbf; ++ibf) {
          FOR_I3 um_it[ibf][i]  = my_u_bc[i];
          FOR_I3 Rd_it[ibf][i]  = my_uu_bc[i];
          FOR_I3 Rod_it[ibf][i] = my_uv_bc[i];
          Tm_it[ibf]            = Tm;
          pm_it[ibf]            = pm;
        }//(ibf < nbf)
      }//"CONSTANT_UREYTP"

      //    inflow profile from a fluent file based on k-e model
      //    PROFILE FLUENT_UTKETP
      //    NOTE: turbulent kinetic enery (tke) is defined as:
      //            tke = 1/2 * (u'u'+v'v'+w'w')
      else if ( token == "FLUENT_UTKETP" ) {
        //      filename
        token = p->getString(iarg++);
        bool use_mean = false;
        if ((iarg < p->size()) && (p->getString(iarg) == "FROM_MEAN")) {
          use_mean = true;
          ++iarg;
        }
        COUT2("  > using profiles of velocity, tke, p, and T from Fluent RANS file: " << token);
        //      read the profile
        FluentBPReader profile(token);
        //      old comment:
        //      set the points, they are contiguous in memory for now, but
        //      they might not always be ...
        //      new comment: Keeping this here. Don't think it is needed.
        double (*x_fa_zone)[3] = new double[zone_ptr->nbf][3];
        for (int ibf =0; ibf < zone_ptr->nbf; ++ibf) {
          FOR_I3 x_fa_zone[ibf][i] = zone_ptr->x_bf[ibf][i];
        }
        profile.setPoints(x_fa_zone,zone_ptr->x_global,zone_ptr->nbf);
        delete[] x_fa_zone;

        // if using mean values, only apply variable name switching to non-turbulence quantities for now (assume these are written without "mean-" prefix from Fluent...?

        profile.ensureVar("turb-kinetic-energy");
        for (int ibf =0; ibf < zone_ptr->nbf; ++ibf) {
          Rd_it[ibf][0] = 2.0/3.0 * profile.getData(ibf,"turb-kinetic-energy");
          Rd_it[ibf][1] = Rd_it[ibf][0];
          Rd_it[ibf][2] = Rd_it[ibf][0];
          //        set off diagonals to zero
          FOR_I3 Rod_it[ibf][i] = 0.0;
        }//(ibf < nbf)

        // now read state variables from "mean-" listed variables
        profile.setUseMean(use_mean);
        profile.ensureVar("x-velocity");
        profile.ensureVar("y-velocity");
        profile.ensureVar("z-velocity");
        profile.ensureVar("absolute-pressure");
        profile.ensureVar("temperature");

        //      get the profile
        for (int ibf =0; ibf < zone_ptr->nbf; ++ibf) {
          um_it[ibf][0] = profile.getData(ibf,"x-velocity");
          um_it[ibf][1] = profile.getData(ibf,"y-velocity");
          um_it[ibf][2] = profile.getData(ibf,"z-velocity");
          pm_it[ibf]    = profile.getData(ibf,"absolute-pressure");
          Tm_it[ibf]    = profile.getData(ibf,"temperature");
        }//(ibf < nbf)
      }//"FLUENT_UTKETP"
      else if ( token == "FLUENT_UREYTP" ) {
        //      filename
        token = p->getString(iarg++);
        bool use_mean = false;
        if ((iarg < p->size()) && (p->getString(iarg) == "FROM_MEAN")) {
          use_mean = true;
          ++iarg;
        }
        COUT2("  > using profiles of velocity, Reynolds stresses, p, and T from Fluent RANS file: " << token);
        //      read the profile
        FluentBPReader profile(token);
        //      set the points
        //      old comment:
        //      set the points, they are contiguous in memory for now, but
        //      they might not always be ...
        //      new comment: Keeping this here. Don't think it is needed.
        double (*x_fa_zone)[3] = new double[zone_ptr->nbf][3];
        for (int ibf =0; ibf < zone_ptr->nbf; ++ibf) {
          FOR_I3 x_fa_zone[ibf][i] = zone_ptr->x_bf[ibf][i];
        }//(ibf < nbf)
        profile.setPoints(x_fa_zone,zone_ptr->x_global,zone_ptr->nbf);
        delete[] x_fa_zone;

        profile.ensureVar("uu-reynolds-stress");
        profile.ensureVar("vv-reynolds-stress");
        profile.ensureVar("ww-reynolds-stress");
        profile.ensureVar("uv-reynolds-stress");
        profile.ensureVar("vw-reynolds-stress");
        profile.ensureVar("uw-reynolds-stress");
        //      get the profile
        for (int ibf =0; ibf < zone_ptr->nbf; ++ibf) {
          //        diagonal
          Rd_it[ibf][0]  = profile.getData(ibf,"uu-reynolds-stress");
          Rd_it[ibf][1]  = profile.getData(ibf,"vv-reynolds-stress");
          Rd_it[ibf][2]  = profile.getData(ibf,"ww-reynolds-stress");
          //        off-diagonal
          Rod_it[ibf][0] = profile.getData(ibf,"vw-reynolds-stress");
          Rod_it[ibf][1] = profile.getData(ibf,"uw-reynolds-stress");
          Rod_it[ibf][2] = profile.getData(ibf,"uv-reynolds-stress");
        }//(ibf < nbf)

        profile.setUseMean(use_mean);
        profile.ensureVar("x-velocity");
        profile.ensureVar("y-velocity");
        profile.ensureVar("z-velocity");
        profile.ensureVar("absolute-pressure");
        profile.ensureVar("temperature");
        //      get the profile
        for (int ibf =0; ibf < zone_ptr->nbf; ++ibf) {
          um_it[ibf][0]  = profile.getData(ibf,"x-velocity");
          um_it[ibf][1]  = profile.getData(ibf,"y-velocity");
          um_it[ibf][2]  = profile.getData(ibf,"z-velocity");
          pm_it[ibf]     = profile.getData(ibf,"absolute-pressure");
          Tm_it[ibf]     = profile.getData(ibf,"temperature");
        }//(ibf < nbf)
      }//"FLUENT_UREYTP"

      //    inflow profile from a fluent file based on k-e model
      //    PROFILE FLUENT_UTKE_CONSTANT_TP <filename> <T> <p>
      //    NOTE: turbulent kinetic enery (tke) is defined as:
      //            tke = 1/2 * (u'u'+v'v'+w'w')
      else if ( token == "FLUENT_UTKE_CONSTANT_TP" ) {
        //      filename
        token = p->getString(iarg++);
        double Tm   = p->getDouble(iarg++);
        double pm   = p->getDouble(iarg++);
        bool use_mean = false;
        if ((iarg < p->size()) && (p->getString(iarg) == "FROM_MEAN")) {
          use_mean = true;
          ++iarg;
        }
        COUT2("  > using velocity and tke profiles from Fluent RANS file: " << token  << " and setting (T,p) = (" <<
              Tm << "," << pm << ")");
        //      read the profile
        FluentBPReader profile(token);
        //      old comment:
        //      set the points, they are contiguous in memory for now, but
        //      they might not always be ...
        //      new comment: Keeping this here. Don't think it is needed.
        double (*x_fa_zone)[3] = new double[zone_ptr->nbf][3];
        for (int ibf =0; ibf < zone_ptr->nbf; ++ibf) {
          FOR_I3 x_fa_zone[ibf][i] = zone_ptr->x_bf[ibf][i];
        }//(ibf < nbf)
        profile.setPoints(x_fa_zone,zone_ptr->x_global,zone_ptr->nbf);
        delete[] x_fa_zone;
        //      check for mean velocity and turbulent kinetic energy

        profile.ensureVar("turb-kinetic-energy");
        //      get the profile
        for (int ibf =0; ibf < zone_ptr->nbf; ++ibf) {
          Rd_it[ibf][0] = 2.0/3.0 * profile.getData(ibf,"turb-kinetic-energy");
          Rd_it[ibf][1] = Rd_it[ibf][0];
          Rd_it[ibf][2] = Rd_it[ibf][0];
          //        set off diagonals to zero
          FOR_I3 Rod_it[ibf][i] = 0.0;
        }//(ibf < nbf)

        profile.setUseMean(use_mean);
        profile.ensureVar("x-velocity");
        profile.ensureVar("y-velocity");
        profile.ensureVar("z-velocity");
        for (int ibf =0; ibf < zone_ptr->nbf; ++ibf) {
          um_it[ibf][0] = profile.getData(ibf,"x-velocity");
          um_it[ibf][1] = profile.getData(ibf,"y-velocity");
          um_it[ibf][2] = profile.getData(ibf,"z-velocity");
          Tm_it[ibf]    = Tm;
          pm_it[ibf]    = pm;
        }//(ibf < nbf)
      }//"FLUENT_UTKE_CONSTANT_TP"
      else if ( token == "FLUENT_UREY_CONSTANT_TP" ) {
        //      filename
        token = p->getString(iarg++);
        double Tm   = p->getDouble(iarg++);
        double pm   = p->getDouble(iarg++);
        bool use_mean = false;
        if ((iarg < p->size()) && (p->getString(iarg) == "FROM_MEAN")) {
          use_mean = true;
          ++iarg;
        }
        COUT2("  > using velocity and Reynolds stresses profiles from Fluent RANS file: " << token
              << " and setting (T,p) = (" << Tm << "," << pm << ")");
        //      read the profile
        FluentBPReader profile(token);
        //      set the points
        //      old comment:
        //      set the points, they are contiguous in memory for now, but
        //      they might not always be ...
        //      new comment: Keeping this here. Don't think it is needed.
        double (*x_fa_zone)[3] = new double[zone_ptr->nbf][3];
        for (int ibf =0; ibf < zone_ptr->nbf; ++ibf) {
          FOR_I3 x_fa_zone[ibf][i] = zone_ptr->x_bf[ibf][i];
        }//(ibf < nbf)
        profile.setPoints(x_fa_zone,zone_ptr->x_global,zone_ptr->nbf);
        delete[] x_fa_zone;

        //      diagonal
        profile.ensureVar("uu-reynolds-stress");
        profile.ensureVar("vv-reynolds-stress");
        profile.ensureVar("ww-reynolds-stress");
        //      off-diagonal;
        profile.ensureVar("uv-reynolds-stress");
        profile.ensureVar("vw-reynolds-stress");
        profile.ensureVar("uw-reynolds-stress");
        //      get the profile
        for (int ibf =0; ibf < zone_ptr->nbf; ++ibf) {
          //        diagonal
          Rd_it[ibf][0]  = profile.getData(ibf,"uu-reynolds-stress");
          Rd_it[ibf][1]  = profile.getData(ibf,"vv-reynolds-stress");
          Rd_it[ibf][2]  = profile.getData(ibf,"ww-reynolds-stress");
          //        off-diagonal
          Rod_it[ibf][0] = profile.getData(ibf,"vw-reynolds-stress");
          Rod_it[ibf][1] = profile.getData(ibf,"uw-reynolds-stress");
          Rod_it[ibf][2] = profile.getData(ibf,"uv-reynolds-stress");
        }//(ibf < nbf)

        //      check for mean velocity
        profile.setUseMean(use_mean);
        profile.ensureVar("x-velocity");
        profile.ensureVar("y-velocity");
        profile.ensureVar("z-velocity");

        for (int ibf =0; ibf < zone_ptr->nbf; ++ibf) {
          um_it[ibf][0]  = profile.getData(ibf,"x-velocity");
          um_it[ibf][1]  = profile.getData(ibf,"y-velocity");
          um_it[ibf][2]  = profile.getData(ibf,"z-velocity");
          Tm_it[ibf]     = Tm;
          pm_it[ibf]     = pm;
        }//(ibf < nbf)
      }//"FLUENT_UREY_CONSTANT_TP"
      else {
        cerr << "Error: unrecognized PROFILE: " << token << endl;
        throw(0);
      } // end of profile
    }//"PROFILE"
    else if(token=="L_MAX"){ // the maximum lengthscale in the inflow zone
      max_length= p->getDouble(iarg++);
      COUT2("  > maximum lengt scale " <<token <<" set to: " << max_length);
    }//"L_MAX"
    else if(token=="L"){
      token=p->getString(iarg++);
      if(token=="HOOK"){
        COUT2("  > In plane length scale set by user in a HOOK ");
        in_plane_length=-1.0;
      }//"HOOK"
      else{
        in_plane_length= p->getDouble(iarg-1);
        COUT2("  > Constant in plane length scale set to: " << in_plane_length );
      }//!"HOOK"
    }//"L"
    else if ( token == "UN") {
      Un = p->getDouble(iarg++) ;
      COUT2("  > Constant convective velocity " <<token <<" set to: " << Un );
    }//"UN"
    else if(token=="TOL"){
      inflow_zero = p->getDouble(iarg++);
      COUT2("  > Inflow turbulence solver tolerance " <<token <<" set to: " << inflow_zero );
    }//"TOL"
    else if(token=="MAX_ITER"){
      inflow_maxiter= p->getInt(iarg++);
      COUT2("  > Inflow turbulence solver interation max " <<token <<" set to: " << inflow_maxiter );
    }//"MAX_ITER"
    else if(token=="RESET"){
      reset=true;
      COUT2("  > Resetting inflow turbulence solver statistitics..." );
    }//"RESET"
    else {
      cerr << "Error: unrecognized token: " << token << endl;
      throw(0);
    }//UNRECOGNIZED
  }//(i < p->size())

  //Now compute the scalings
  double my_mean[3]={0.0, 0.0, 0.0};
  for (int ibf =0; ibf < zone_ptr->nbf; ++ibf) {
    double normal[3];
    double um[3];
    FOR_I3{
      um[i]=um_it[ibf][i];
      normal[i]=zone_ptr->n_bf[ibf][i];
    }//FOR_I3
    double u_dot_A=DOT_PRODUCT(um,normal);
    double this_rho=pm_it[ibf]/(solver->R_gas*Tm_it[ibf]);
    double this_sos= sqrt(solver->gamma*pm_it[ibf]/this_rho);
    assert(u_dot_A<=0.0);
    my_mean[0] += zone_ptr->area_bf[ibf];
    my_mean[1] -= u_dot_A;
    my_mean[2] -= u_dot_A/this_sos;
  }//(ibf < nbf)
  double mean[3];
  MPI_Allreduce(my_mean,mean,3,MPI_DOUBLE,MPI_SUM,mpi_comm);
  Umean_it = mean[1]/mean[0];
  Ma_it    = mean[2]/mean[0];
}//setStatistics()

#define FOR_ICV_ for (int icv_ = 0; icv_ < ncv_; ++icv_)
#define FOR_IFA_ for (int ifa_ = 0; ifa_ < nfa_; ++ifa_)

//update the boundary state and provide out of plane correlations
//TODO give the user a separate option to specify the out of plane length scale
void InflowTurbulenceIG::updateBC(double (*FacePtrNew)[3],double (*FacePtrOld)[3]) {
  //Here we update the Boundary condition

  //convective time scale set based on u_mean and the in plane length scale or a constant convective velocity as specified by the user
  double gm=solver->gamma;
  double gm1=gm-1.0;
  for (int ibf =0; ibf < zone_ptr->nbf; ++ibf) {
    double length=length_scale[ibf];// Should add some factor here to get the ratio between out of plane and inplane length scales. Alternatively the user could specify...
    double Tper;
    if(Un==1.0E+200){ // this means that the user has not specified the convective velocity
      double um[3];
      double normal[3];
      FOR_I3{
        um[i]=um_it[ibf][i];
        normal[i]=zone_ptr->n_bf[ibf][i];
      }
      NORMALIZE(normal);
      double vel=DOT_PRODUCT(um,normal);
      Tper  = length/abs(vel);
    }//(Un == 1.0E+200)
    else{
      Tper=length/Un;
    }//(Un == 1.0E+200)
    double cf_old = exp(-0.5*M_PI*solver->dt/Tper) ;
    double cf_new = sqrt( 1.0 - exp(-M_PI*solver->dt/Tper)) ;

    double u_prime=0.0;
    FOR_I3 {
      FacePtrOld[ibf][i]=cf_new*FacePtrNew[ibf][i]+cf_old*FacePtrOld[ibf][i];
      u_prime+=FacePtrOld[ibf][i]*FacePtrOld[ibf][i];
    }//FOR_I3

    //  thermodynamic quantities
    u_prime = sqrt(u_prime/3.0);
    double T_prime = 0.0; //-gm1 * Ma_it*Ma_it * u_prime/Umean_it * Tm_it[ibf];
    //  limit Tprime to 90% of Tmean
    T_prime = max(T_prime,-0.9*Tm_it[ibf]);  // to have a positive temperature
    T_prime = min(T_prime,0.9*Tm_it[ibf]);   // to have a positive density

    T_bc[ibf]=T_prime+Tm_it[ibf]; // we have the temperature fluctuations here
    p_bc[ibf]=pm_it[ibf];
    if(pm_it[ibf]==0.0){ // This should only happen if we are using the turbulence trip optionand we treat it as adiabatic...ish
      T_bc[ibf]=solver->T[zone_ptr->cvobf[ibf]];
      p_bc[ibf]=solver->p[zone_ptr->cvobf[ibf]];
    }//(pm_it[ibf] == 0.0)
    FOR_I3{
      bf_state[ibf].u[i] = FacePtrOld[ibf][i]+um_it[ibf][i];
    }//FOR_I3
    const double p_it = p_bc[ibf]; // - 0.5*DOT_PRODUCT(FacePtrOld[ibf],FacePtrOld[ibf])*p_bc[ibf]/(solver->R_gas*T_bc[ibf]); //solver->p[zone_ptr->cvobf[ibf]];
    //bf_state[ibf].sp_vol = solver->R_gas*T_bc[ibf]/p_bc[ibf];
    bf_state[ibf].sp_vol = solver->R_gas*T_bc[ibf]/p_it;
    //bf_state[ibf].p      = p_bc[ibf];
    bf_state[ibf].p      = p_it;
    bf_state[ibf].h      = solver->R_gas*T_bc[ibf]*gm/gm1;
  }//(ibf < nbf)
}//updateBc()

void InflowTurbulenceIG::query(const string& param_str) const {
  double my_buf[5] = {0.0,0.0,0.0,0.0,0.0};
  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
    //const int icv = zone_ptr->cvobf[ibf];
    if (zone_ptr->area_bf[ibf] > 0.0) {
      my_buf[0] += zone_ptr->area_bf[ibf];
      // normal is pointing out, so take -ve dot product...
      const double un = -DOT_PRODUCT(bf_state[ibf].u,zone_ptr->n_bf[ibf])/zone_ptr->area_bf[ibf];
      my_buf[1] += un*zone_ptr->area_bf[ibf];
      my_buf[2] += un*un*zone_ptr->area_bf[ibf];
      my_buf[3] += solver->p[zone_ptr->cvobf[ibf]]*zone_ptr->area_bf[ibf];
      my_buf[4] += solver->p[zone_ptr->cvobf[ibf]]*solver->p[zone_ptr->cvobf[ibf]]*zone_ptr->area_bf[ibf];
    }
  }
  double buf[5];
  MPI_Reduce(my_buf,buf,5,MPI_DOUBLE,MPI_SUM,0,mpi_comm);
  if ( mpi_rank == 0 ) {
    const double un_avg = buf[1]/buf[0];
    const double p_avg = buf[3]/buf[0];
    cout << "[QUERY_BC] " << zone_ptr->getName() << " time: "  <<
      solver->time <<
      " un_avg, un_rms: " << un_avg << " " << sqrt(max(0.0,buf[2]/buf[0]-un_avg*un_avg)) <<
      " p_cv_avg, p_cv_rms: " << p_avg << " " << sqrt(max(0.0,buf[4]/buf[0]-p_avg*p_avg)) << endl;
  }
}

int InflowTurbulenceIG::solveVecBoundaryCvCg(double (*sol)[3],const double * const A,double (*rhs)[3], double * AddedDiag,double zero, int maxiter ) {
  bool verbose = false;
  if (solver->step%solver->check_interval==0){
    verbose=true;
  }

  //we need the following work arrays...
  double (*res)[3]      = new double[ncv_][3];
  double (*v)[3]        = new double[ncv_][3];
  double (*p)[3]      = new double[ncv_g_][3];
  double * inv_diag = new double[ncv_];

  double rho[3];
  FOR_I3 rho[i]=1.0;
  for (int icv_ = 0; icv_ < ncv_g_; ++icv_) {
    FOR_I3{
      p[icv_][i]=0.0;
    }//FOR_I3
  }//FOR_ICV_G_

  for (int icv_ = 0; icv_ < ncv_; ++icv_) {
    FOR_I3{
      res[icv_][i] = rhs[icv_][i]-AddedDiag[icv_]*sol[icv_][i];
      inv_diag[icv_]=AddedDiag[icv_];
    }//FOR_I3
  }//FOR_ICV_
  for(int ifa_=0;ifa_<nfa_i_;ifa_++){
    int icv0_=cvofa_[ifa_][0];
    int icv1_=cvofa_[ifa_][1];
    inv_diag[icv0_]+=A[ifa_];
    inv_diag[icv1_]+=A[ifa_];
    FOR_I3{
      double add=A[ifa_]*(sol[icv1_][i]-sol[icv0_][i]);
      res[icv0_][i]+=add;
      res[icv1_][i]-=add;
    }//FOR_I3
  }//(ifa < nfa_i)

  for(int ifa_=nfa_i_;ifa_<nfa_;ifa_++){  // special treatment needed for ghosts
    int icv0_=cvofa_[ifa_][0];
    int icv1_=cvofa_[ifa_][1];
    assert(icv1_>=ncv_&&icv0_<ncv_);
    inv_diag[icv0_]+=A[ifa_];
    FOR_I3{
      double add=A[ifa_]*(sol[icv1_][i]-sol[icv0_][i]);
      res[icv0_][i]+=add;
    }//FOR_I3
  }//(ifa < nfa_)

  for (int icv_ = 0; icv_ < ncv_; ++icv_) {
    inv_diag[icv_] = 1.0/(inv_diag[icv_] );
  }
  //diagonal precon/compute normalized residual...
  for (int icv_ = 0; icv_ < ncv_; ++icv_) {
    FOR_I3{
      v[icv_][i] = res[icv_][i]*inv_diag[icv_];
    }//FOR_I3
  }//FOR_ICV_
  int iter = 0;
  int done = 0;
  while (done == 0) {
    ++iter;

    double rho_prev[3];
    FOR_I3 rho_prev[i]=rho[i];

    FOR_I3{
      if (fabs(rho_prev[i]) < 1.0E-20)
        rho_prev[i] = -1.0E-20; // -1.0E-20? seems to help
    }//FOR_I3
    double my_rho[3] = {0.0,0.0,0.0};
    for (int icv_ = 0; icv_ < ncv_; ++icv_) {
      FOR_I3{
        my_rho[i] += res[icv_][i]*v[icv_][i];
      }//FOR_I3
    }//FOR_ICV_
    MPI_Allreduce(my_rho,rho,3,MPI_DOUBLE,MPI_SUM,mpi_comm);

    double beta[3];
    FOR_I3{
      beta[i]=rho[i]/rho_prev[i];
    }//FOR_I3

    for (int icv_ = 0; icv_ < ncv_; ++icv_){
      FOR_I3 p[icv_][i] = v[icv_][i] + beta[i]*p[icv_][i];
    }//FOR_ICV_

    updateCvData_(p);

    //  v = [Ap]{p}... Matrix vector product is done using faces.
    for (int icv_ = 0; icv_ < ncv_; ++icv_){
      FOR_I3 v[icv_][i]=AddedDiag[icv_]*p[icv_][i];
    }
    for(int ifa_=0;ifa_<nfa_i_;ifa_++){
      int icv0_=cvofa_[ifa_][0];
      int icv1_=cvofa_[ifa_][1];
      FOR_I3{
        double add=A[ifa_]*(p[icv0_][i]-p[icv1_][i]);
        v[icv0_][i]+=add;
        v[icv1_][i]-=add;
      }//FOR_I3
    }//(ifa_ < nfa_i_)
    for(int ifa_=nfa_i_;ifa_<nfa_;ifa_++){  // special treatment needed for ghosts
      int icv0_=cvofa_[ifa_][0];
      int icv1_=cvofa_[ifa_][1];
      assert(icv1_>=ncv_&&icv0_<ncv_);
      FOR_I3{
        double add=A[ifa_]*(p[icv0_][i]-p[icv1_][i]);
        v[icv0_][i]+=add;
      }//FOR_I3
    }//(ifa < nfa_)

    double my_gamma[3] = {0.0,0.0,0.0};
    for (int icv_ = 0; icv_ < ncv_; ++icv_){
      FOR_I3 {
        my_gamma[i] += p[icv_][i]*v[icv_][i];
      }//FOR_I3
    }//FOR_ICV_

    double gamma[3];
    MPI_Allreduce(my_gamma,gamma,3,MPI_DOUBLE,MPI_SUM,mpi_comm);
    FOR_I3{
      if (fabs(gamma[i]) < 1.0E-20)
        gamma[i] = 1.0E-20;
    }//FOR_I3
    const double alpha[3]={rho[0]/gamma[0],rho[1]/gamma[1],rho[2]/gamma[2]};

    //  check if we are done...
    if (iter%3 == 0) {
      for (int icv_ = 0; icv_ < ncv_; ++icv_) {
        FOR_I3 sol[icv_][i] += alpha[i]*p[icv_][i];
      }//FOR_ICV_
      updateCvData_(sol);
      //    recompute the residual...
      for (int icv_ = 0; icv_ < ncv_; ++icv_) {
        FOR_I3{
          res[icv_][i] = rhs[icv_][i]-AddedDiag[icv_]*sol[icv_][i];
        }//FOR_I3
      }//FOR_ICV_
      for(int ifa_=0;ifa_<nfa_i_;ifa_++){
        int icv0_=cvofa_[ifa_][0];
        int icv1_=cvofa_[ifa_][1];
        FOR_I3{
          double add=A[ifa_]*(sol[icv1_][i]-sol[icv0_][i]);
          res[icv0_][i]+=add;
          res[icv1_][i]-=add;
        }//FOR_I3
      }//FOR_IFA_
      for(int ifa_=nfa_i_;ifa_<nfa_;ifa_++){  // special treatment needed for ghosts
        int icv0_=cvofa_[ifa_][0];
        int icv1_=cvofa_[ifa_][1];
        assert(icv1_>=ncv_&&icv0_<ncv_);
        FOR_I3{
          double add=A[ifa_]*(sol[icv1_][i]-sol[icv0_][i]);
          res[icv0_][i]+=add;
        }//FOR_I3
      }//FOR_IFA_B

      for (int icv_ = 0; icv_ < ncv_; ++icv_) {
        FOR_I3{
          v[icv_][i] = res[icv_][i]*inv_diag[icv_];
        }//FOR_I3
      }//FOR_ICV_
      double  my_res_max = 0.0;
      for (int icv_ = 0; icv_ < ncv_; ++icv_) {
        FOR_I3{
          my_res_max = max( my_res_max, fabs(v[icv_][i]) );
        }//FOR_I3
      }//FOR_ICV_
      double res_max=0.0;
      MPI_Reduce(&my_res_max, &res_max, 1, MPI_DOUBLE, MPI_MAX, 0, mpi_comm);

      if (mpi_rank == 0) {
        //      only share the last half of the convergence behaviour...
        if (verbose || (iter > maxiter/2))
          cout << " > solveVecBoundaryCvCg for iter " << iter << " res_max " << res_max << endl;
        if (res_max <= zero) {
          if(verbose) cout << "-> Successfully converged error to " << res_max << endl;
          done = 1;
        }//(res_max < = zero)
        else if (iter > maxiter && verbose) {
          cout << "Warning: solveVecBoundaryCvCg did not converge after " << maxiter <<
            " iters, res_max: " << res_max << endl;
          done = 2;
        }//(iter > maxiter && verbose)
      }//(root)
      MPI_Bcast(&done, 1, MPI_INT, 0, mpi_comm);
    }//(iter%3 == 0)
    else {
      //    update full phi including ghosts...
      for (int icv_ = 0; icv_ < ncv_g_; ++icv_){
        FOR_I3 sol[icv_][i] += alpha[i]*p[icv_][i];
      }//FOR_ICV_G_

      for (int icv_ = 0; icv_ < ncv_; ++icv_) {
        //      on the other iterations, use this approximation to update
        //      the unreduced residual...
        FOR_I3 res[icv_][i] -= alpha[i]*v[icv_][i];
        //      still need to compute v, diag precon for next iteration...
        FOR_I3 v[icv_][i] = res[icv_][i]*inv_diag[icv_];
      }//FOR_ICV_
    }//(iter%3 != 0)
  }//(done == 0)
  delete[] res;
  delete[] v;
  delete[] p;
  delete[] inv_diag;

  //let the calling routine know if we were successful...
  return( done == 1 );
}//solveVecBoundaryCvCg()

///////////////////////////////////////////////////////////////////////////////////
//            REDUNDANT CODE
///////////////////////////////////////////////////////////////////////////////////
void InflowTurbulenceIG::updateInflowTurbulence(double (*FacePtr)[3]){
  //  updates the random field used in the inflow turbulence calculation and pushes it to the boundary
  double (*RandomField)[3]=new double[zone_ptr->nbf][3];
  updateRandomField(RandomField);
  //  removes the mean and sets variance to one for the boundary field
  removeMeanAndSetVarianceOne(RandomField);
  //  cholesky decomposition
  computeCholesky(RandomField);
  //  update the boundary conditions. The turbulent field is stored in u_bc// should rename
  updateBC(RandomField,FacePtr);
  delete [] RandomField;
}//updateInflowTurbulence()

void InflowTurbulenceIG::buildITBoundaryStructures() {
  COUT1("Building connectivity system..."); // need  to figure out how to do output when rank 0 isn't involved
  //  allocte a cv flag
  int ncv=solver->ncv;
  int ncv_g=solver->ncv_g;
  int* cv_flag=new int[solver->ncv_g];
  FOR_ICV_G cv_flag[icv]=0;
  //  step one. Flag all boundary cvs that are adjecent to a wall
  for(vector<IdealGasBc*>::iterator it = solver->bcs.begin(); it != solver->bcs.end(); ++it){
    string s=(*it)->getName();
    Param* p=getParam(s);
    //    cout<<p->getString()<<endl;
    if(p->getString()=="WALL_ADIABATIC"||p->getString()=="WALL_ISOTHERMAL"||p->getString()=="WM_ALG_ADIABATIC"||p->getString() == "WM_ALG_ISOTHERMAL"){
      //      TODO This needs to recognize all wall zones!!
      //      cout<<(*it)->getName()<<endl;
      BfZone* zone_ptr_it=(*it)->zone_ptr;
      for (int ibf =0; ibf < zone_ptr_it->nbf; ++ibf) {
        cv_flag[zone_ptr_it->cvobf[ibf]]=1;
      }//(ibf < nbf)
    }//"WALL_ZONES"
  }///(it in bcs)
  for (int ibf =0; ibf < zone_ptr->nbf; ++ibf) {
    cv_flag[zone_ptr->cvobf[ibf]]=cv_flag[zone_ptr->cvobf[ibf]]+2;
  }
  //  Now all cvs that only are attached to a wall are 1, all cvs that are attached to a wall and the inflow
  //  turbulence zone are 3 and all other inflow turbulence faces are 2. This does not take into account that other-non-neighbour walls may be important
  for (int icv = 0; icv < solver->ncv; ++icv){
    if(cv_flag[icv]>1){
      bool found=false;
      int bfa_f=solver->bfocv_i[icv];
      int bfa_l=solver->bfocv_i[icv+1];
      for(int ifa=bfa_f;ifa<bfa_l;ifa++){
        int bifa=solver->bfocv_v[ifa];
        if(solver->zone_bf[bifa]==zone_ptr->index){
          found=true;
        }//(solver->zone_bf[bifa] == zone_ptr->index)
      }//(ifa < bfa_l)
      assert(found);
    }//(cv_flag[icv] > 1)
    else{
      bool found=false;
      int bfa_f=solver->bfocv_i[icv];
      int bfa_l=solver->bfocv_i[icv+1]-1;
      for(int ifa=bfa_f;ifa<bfa_l;ifa++){
        int bifa=solver->bfocv_v[ifa];
        if(solver->zone_bf[bifa]==zone_ptr->index){
          found=true;
        }//(solver->zone_bf[bifa] == zone_ptr->index)
      }//(ifa < bfa_l)
      assert(!found);
    }//(cv_flag[icv] <= 1)
  }//FOR_ICV

  //  update the flag
  solver->updateCvData(cv_flag);
  //  Reaching this point means that the correct cvs have been tagged.

  //  we loop over the boundary faces and set the wall distance to -2.0 in the bfs that are in a cv next to the wall and a -1.0 everywhere else.
  //  The lengthscale will be set later
  for (int ibf =0; ibf < zone_ptr->nbf; ++ibf) {
    if(cv_flag[zone_ptr->cvobf[ibf]]==3)  length_scale[ibf]=-2.0;
    else                                  length_scale[ibf]=-1.0;
  }//(ibf < nbf)

  ncv_ = 0;
  FOR_ICV {
    if(cv_flag[icv]>1)  cv_flag[icv] = ncv_++;
    else                cv_flag[icv] = -1;
  }//FOR_ICV

  ncv_g_ = ncv_;
  for (int icv = solver->ncv; icv < solver->ncv_g; ++icv) {
    if(cv_flag[icv]>1){
      cv_flag[icv] = ncv_g_++;
    }//(cv_flag[icv] > 1)
    else {
      cv_flag[icv] = -1;
    }//(cv_flag[icv] <= 1)
  }//FOR_ICV

  activeCvs = new int[ncv_g_];
  int icv_ = 0;
  FOR_ICV_G {
    if(cv_flag[icv] >= 0){
      assert(cv_flag[icv] == icv_);
      activeCvs[icv_]=icv;
      icv_++;
    }//(cv_flag[icv] >= 0)
  }//FOR_ICV_G
  assert(icv_ == ncv_g_);

  //  We also store the face indices of the faces between the involved cvs
  //  first we need to know how many there are. (should be the number of neighbours-self divided by two)
  nfa_i_ = 0;
  for (int ifa = 0; ifa < solver->nfa_i; ++ifa) {
    if((cv_flag[solver->cvofa[ifa][0]]>=0)&&(cv_flag[solver->cvofa[ifa][1]]>=0)){
      nfa_i_++;
    }//(cv_flag >= 0)
  }//(ifa < nfa_i)
  nfa_ = nfa_i_;
  for (int ifa = solver->nfa_i; ifa < solver->nfa; ++ifa) {
    if((cv_flag[solver->cvofa[ifa][0]]>=0)&&(cv_flag[solver->cvofa[ifa][1]]>=0)){
      nfa_++;
    }//(cv_flag >= 0)
  }//(ifa < nfa)

  activeFaces = new int[nfa_];
  assert(cvofa_ == NULL);
  cvofa_ = new int[nfa_][2];
  int ifa_ = 0;
  for (int ifa = 0; ifa < solver->nfa; ++ifa) {
    if((cv_flag[solver->cvofa[ifa][0]]>=0)&&(cv_flag[solver->cvofa[ifa][1]]>=0)){
      cvofa_[ifa_][0] = cv_flag[solver->cvofa[ifa][0]];
      cvofa_[ifa_][1] = cv_flag[solver->cvofa[ifa][1]];
      activeFaces[ifa_] = ifa;
      ++ifa_;
    }//(cf_flag >= 0)
  }//(ifa < nfa)
  assert(ifa_ == nfa_);

  //  we also build a local cvobf
  assert(cvobf_ == NULL);
  cvobf_ = new int[zone_ptr->nbf];
  for (int ibf =0; ibf < zone_ptr->nbf; ++ibf) {
    cvobf_[ibf]=cv_flag[zone_ptr->cvobf[ibf]];
  }//(ibf < nbf)

  //  Clean up
  delete [] cv_flag;
}//buildITBoundaryStructures()

//This routine loops over all bfs and computes the distance based on all known near wll bfs. It may become exceedingly slow....
void InflowTurbulenceIG::setLengthscale() {
  COUT2(" > set lengthscale for inflow turbulence"); //TODO should change this to the active communicator once this is in place

  const int nbf = zone_ptr->nbf;

  if (in_plane_length==1.0E+200){ // A constant lengthscale is not specified and we calculate the distance to the wall with a maximum threshold L_MAX
    double * wall_dist = new double[nbf];  // compute for everybody to reuse machinery

    int ncv=solver->ncv;
    double * cv_flag=new double[ncv];
    FOR_ICV cv_flag[icv]=-2.0;

    FOR_IBF {
      cv_flag[zone_ptr->cvobf[ibf]] = -1;
    }//(ibf < nbf)

    //    Now we have flagged the cvs and we can loop over the boundary conditions to compute the distance to the cvs adjecent to the walls
    for(vector<IdealGasBc*>::const_iterator it = solver->bcs.begin(); it != solver->bcs.end(); ++it) {
      string s = (*it)->getName();
      Param* p = getParam(s);
      if ((p->getString()=="WALL_ADIABATIC") || (p->getString()=="WALL_ISOTHERMAL") || (p->getString()=="WM_ALG_ADIABATIC") || (p->getString()=="WM_ALG_ISOTHERMAL")) {
        // TODO This must recognize all wall bcs. How?
        BfZone* zone_ptr_it=(*it)->zone_ptr;
        for (int ibf=0; ibf < zone_ptr_it->nbf; ++ibf) {
          const int icv = zone_ptr_it->cvobf[ibf];

          if (cv_flag[icv] == -1.0) {
            // first time hit, identify as a corner cell
            cv_flag[icv] = HUGE_VAL;
          }

          if (cv_flag[icv] > 0.0) {
            // we are at a zone corner cv, so compute wall-distance and store minimum in cv_flag
            const double dist = DIST(zone_ptr_it->x_bf[ibf],solver->x_cv[icv]);
            if (dist < cv_flag[icv]) {
              cv_flag[icv] = dist;
            }
          }
        }
      }
    }  // bcs iterator
    // flagged cvs ( > 0) have their min wall distance stored
    solver->computeBfDistanceFromFlaggedCvs(wall_dist,zone_ptr,cv_flag);
    DELETE(cv_flag);

    //    Now we have all the boundary points and their distance to the wall.
    //    The next step is to loop over all boundary faces and set the length
    FOR_IBF {
      if (length_scale[ibf] < 0.0) {
        // is there a condition where some are set but not others?
        length_scale[ibf] = min(wall_dist[ibf],max_length);  // up to L_MAX
      }
      assert(length_scale[ibf] != 1.0E+200);
    } // FOR_IBF
    //    And we are done
    DELETE(wall_dist);
  }//(in_plane_length=1.0E+200)
  else if (in_plane_length > 0) { //there is a constant lengthscale
    FOR_IBF {
      length_scale[ibf] = in_plane_length;
    }
  }//(in_plane_length > 0)
  else this->inflowTurbulenceLengthscaleHook(length_scale);
}//setLengthScale()

void InflowTurbulenceIG::buildFilterOperators() {
  assert(A_fa_ == NULL);      A_fa_ = new double[nfa_];
  assert(A_diag_cv_ == NULL); A_diag_cv_ = new double[ncv_];

  //  get the length scale from the boundary
  double *lengthScale_     = new double[ncv_g_];
  for (int ibf =0; ibf < zone_ptr->nbf; ++ibf)
    lengthScale_[cvobf_[ibf]]=length_scale[ibf];
  updateCvData_(lengthScale_);

  double max_ls_over_dx_local=0.0;
  //  loop to compute l over dx and
  for (int ifa_ = 0; ifa_ < nfa_; ++ifa_) {
    const int icv0_ = cvofa_[ifa_][0]; assert((icv0_ >= 0)&&(icv0_ < ncv_g_));
    const int icv1_ = cvofa_[ifa_][1]; assert((icv1_ >= 0)&&(icv1_ < ncv_g_));
    double ls=0.5*(lengthScale_[icv0_]+lengthScale_[icv1_]);
    double x0[3];
    double x1[3];
    FOR_I3 {
      x0[i]=solver->x_cv[activeCvs[icv0_]][i];
      x1[i]=solver->x_cv[activeCvs[icv1_]][i];
    }//FOR_I3
    double dx=DIST(x1,x0);
    A_fa_[ifa_]=ls/dx;
    if(A_fa_[ifa_]>max_ls_over_dx_local)
      max_ls_over_dx_local=A_fa_[ifa_];
  }//(ifa_ < nfa_)

  double max_ls_over_dx;
  MPI_Allreduce(&max_ls_over_dx_local, &max_ls_over_dx, 1, MPI_DOUBLE, MPI_MAX, mpi_comm);

  FOR_ICV_ A_diag_cv_[icv_] = 0.0;

  for (int ifa_ = 0; ifa_ < nfa_i_; ++ifa_) {
    const int ifa = activeFaces[ifa_]; assert((ifa >= 0)&&(ifa < solver->nfa_i));
    const double area = MAG(solver->n_fa[ifa]);
    const int icv0_ = cvofa_[ifa_][0]; assert((icv0_ >= 0)&&(icv0_ < ncv_));
    const int icv1_ = cvofa_[ifa_][1]; assert((icv1_ >= 0)&&(icv1_ < ncv_));
    double ls=0.5*(lengthScale_[icv0_]+lengthScale_[icv1_]);
    double Ld=ls/A_fa_[ifa_];
    A_fa_[ifa_]=ls*(ls/Ld)*area;
    A_diag_cv_[icv0_] += area;
    A_diag_cv_[icv1_] += area;
  }//(ifa_ < nfa_i_)

  for (int ifa_ = nfa_i_; ifa_ < nfa_; ++ifa_) {
    const int ifa = activeFaces[ifa_]; assert((ifa >= solver->nfa_i)&&(ifa < solver->nfa));
    const double area = MAG(solver->n_fa[ifa]);
    const int icv0_ = cvofa_[ifa_][0]; assert((icv0_ >= 0)&&(icv0_ < ncv_));
    const int icv1_ = cvofa_[ifa_][1]; assert((icv1_ >= ncv_)&&(icv1_ < ncv_g_));
    double ls=0.5*(lengthScale_[icv0_]+lengthScale_[icv1_]);
    double Ld=ls/A_fa_[ifa_];
    A_fa_[ifa_]=ls*(ls/Ld)*area;
    A_diag_cv_[icv0_] +=area;
  }//(ifa_ < nfa_)

  delete[] lengthScale_;
}//buildFilterOperators()

void InflowTurbulenceIG::initiateStats() { // this must be smart enough to recognize a restart
  if (!zone_ptr->checkDataFlag("u_bc")) {  // if the velocity isn't read in it must be set to 0 and the statistics should be reset
    reset=true;
    for (int ibf =0; ibf < zone_ptr->nbf; ++ibf) {
      FOR_I3 u_bc[ibf][i]=0.0;
    }//(ibf < zone_ptr->nbf)
  }//("u_bc")
  if(solver->step==0 || reset){
    for (int ibf =0; ibf < zone_ptr->nbf; ++ibf) {
      FOR_I3{
        avg[ibf][i]=0.0;
        rms[ibf][i]=0.0;
      }//FOR_I3
    }//(ibf < nbf)
    IturbWgt    = 0.0;
  }//(solver->step >= 0 || reset)
  else{
    assert(IturbWgt!=0.0);
  }
}//initiateStats()

void InflowTurbulenceIG::updateRandomField(double (*FacePtr)[3]) {
  //  populate random noise
  double (*Rno_)[3]     = new double[ncv_g_][3];
  FOR_ICV_ FOR_I3
    Rno_[icv_][i] = max(-3.0,min(3.0,randn())); // clip at +/- 3? sigma

  //  initialize solution
  double (*sol_)[3]     = new double[ncv_g_][3];
  for(int icv_=0; icv_<ncv_g_;++icv_)
    FOR_I3
      sol_[icv_][i]=0.0;

  //  rescale the source term to remove the mean
  double my_sum[3]={0.0, 0.0, 0.0};
  for (int icv_ = 0; icv_ < ncv_; ++icv_) {
    FOR_I3{
      Rno_[icv_][i]=Rno_[icv_][i]*A_diag_cv_[icv_];
      my_sum[i]=my_sum[i]+Rno_[icv_][i];
    }//FOR_I3
  }//FOR_ICV_
  double sum[3];
  MPI_Allreduce(my_sum,sum,3,MPI_DOUBLE,MPI_SUM,mpi_comm);
  int tot_ncv_;
  MPI_Allreduce(&ncv_, &tot_ncv_, 1, MPI_INT, MPI_SUM, mpi_comm);
  double deduction[3]={ sum[0]/double(tot_ncv_), sum[1]/double(tot_ncv_), sum[2]/double(tot_ncv_)};

  for (int icv_ = 0; icv_ < ncv_; ++icv_)
    FOR_I3
      Rno_[icv_][i]=Rno_[icv_][i]-deduction[i];

  //  solve equation
  solveVecBoundaryCvCg(sol_,A_fa_,Rno_,A_diag_cv_,inflow_zero,inflow_maxiter); // TODO returns an error code; should be handled

  //  push data to the boundary
  getCvVectorData_(FacePtr,sol_);

  //  clean up
  delete [] sol_;
  delete [] Rno_;
}//updateRandomField()

void InflowTurbulenceIG::removeMeanAndSetVarianceOne(double (*FacePtr)[3]) {
  //  I am unsure if these normalizations are the right thing to do. Statistics are local so the variance and mean are pointwise properties

  //  need mean = 0  and var 1 in time for each point for the Cholesky decoposition to work. This will decorrelate the signal slightly....

  //  first we update the statistics...
  const double this_wgt = solver->dt;
  const double old_wgt  = IturbWgt;
  IturbWgt                  += this_wgt;

  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
    FOR_I3 {
      rms[ibf][i] = (old_wgt * (rms[ibf][i] * rms[ibf][i] + avg[ibf][i] * avg[ibf][i]) + this_wgt * (FacePtr[ibf][i] * FacePtr[ibf][i])) / IturbWgt;
      avg[ibf][i] = (old_wgt * avg[ibf][i] + this_wgt * FacePtr[ibf][i]) / IturbWgt;
    }//FOR_I3
  }//(ibf < nbf)

  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf)
    FOR_I3
      rms[ibf][i] = sqrt(abs(rms[ibf][i] - avg[ibf][i] * avg[ibf][i]));

  //  Now we have updated statistics and we can correct the signal
  for (int ibf = 0; ibf < zone_ptr->nbf; ++ibf) {
    FOR_I3{
      double a=rms[ibf][i];
      //      for the first steps the rms may be zero, so we need to check
      if(a==0.0) a=1.0;
      FacePtr[ibf][i]=(FacePtr[ibf][i]-avg[ibf][i])/a;
    }//FOR_I3
  }//(ibf < nbf)
}//removeMeanAndSetVariance()

void InflowTurbulenceIG::computeCholesky(double (*FacePtr)[3]) {
  //  TODO This routine should be merged into updateBC to gain some speed
  //  build the Cholesky decomposition of the Reynolds stresses
  //    REY = |<u'u'>   0      0   |
  //          |<u'v'> <v'v'>   0   |
  //          |<u'w'> <v'w'> <w'w'>|
  //  where
  //    Rd  = [<u'u'> <v'v'> <w'w'>]
  //    Rod = [<v'w'> <u'w'> <u'v'>]
  for (int ibf =0; ibf < zone_ptr->nbf; ++ibf) {
    assert(Rd_it[ibf][0]>=0.0);
    double a11 =  sqrt(Rd_it[ibf][0]);
    double a21 = 0.0;
    if ( a11 > 0.0 ) a21 = Rod_it[ibf][2] / a11;
    double a22 = sqrt( max(Rd_it[ibf][1]-a21*a21,0.0) );
    double a31 = 0.0;
    if ( a11> 0.0 ) a31 = Rod_it[ibf][1] / a11;
    double a32 = 0.0;
    if ( a22 > 0.0 ) a32 =  (Rod_it[ibf][0] - a21*a31) / a22;
    double a33 = sqrt(  max(Rd_it[ibf][2]-a31*a31-a32*a32,0.0) );

    //    velocity component
    double my_u_bc[3];
    my_u_bc[0]  = a11 * FacePtr[ibf][0];
    my_u_bc[1]  = a21 * FacePtr[ibf][0] + a22 * FacePtr[ibf][1];
    my_u_bc[2]  = a31 * FacePtr[ibf][0] + a32 * FacePtr[ibf][1] + a33 * FacePtr[ibf][2];

    FOR_I3 FacePtr[ibf][i]=my_u_bc[i];
  }//(ibf < nbf)
}//computeCholesky()

void InflowTurbulenceIG::updateCvData_(double (*var_)[3]) {
  double (*var)[3] = new double[solver->ncv_g][3];

  //  for debugging, fill the whole thisn with something that should bring down the house...
  for (int icv = 0; icv < solver->ncv_g; ++icv) FOR_I3 var[icv][i] = 1.0E+200;
  FOR_ICV_ {
    const int icv = activeCvs[icv_]; assert((icv >= 0)&&(icv < solver->ncv));
    FOR_I3 var[icv][i] = var_[icv_][i];
  }//FOR_ICV_
  solver->updateCvData(var);
  for (int icv_ = ncv_; icv_ < ncv_g_; ++icv_) {
    const int icv = activeCvs[icv_]; assert((icv >= solver->ncv)&&(icv < solver->ncv_g));
    FOR_I3 {
      var_[icv_][i] = var[icv][i];
      assert(var_[icv_][i] != 1.0E+200);
    }//FOR_I3
  }//FOR_ICV_

  delete[] var;
}//updateCvData()

//TODO This routine should be updated to use a local communicator to speed up the inflow turbulence
void InflowTurbulenceIG::updateCvData_(double *var_) {
  double (*var) = new double[solver->ncv_g];

  //  for debugging, fill the whole thisn with something that should bring down the house...
  for (int icv = 0; icv < solver->ncv_g; ++icv) var[icv] = 1.0E+200;

  FOR_ICV_ {
    const int icv = activeCvs[icv_]; assert((icv >= 0)&&(icv < solver->ncv));
    var[icv] = var_[icv_];
  }//FOR_ICV_
  solver->updateCvData(var);
  for (int icv_ = ncv_; icv_ < ncv_g_; ++icv_) {
    const int icv = activeCvs[icv_]; assert((icv >= solver->ncv)&&(icv < solver->ncv_g));
    var_[icv_] = var[icv];
    assert(var_[icv_] != 1.0E+200);
  }//FOR_ICV_G

  delete[] var;
}//updateCvData()

//two pointers as input. Moves data from the cv pointer to the right face.
void InflowTurbulenceIG::getCvData_(double *FacePtr, double *CvPtr_){
  for (int ibf =0; ibf < zone_ptr->nbf; ++ibf) {
    FacePtr[ibf]=CvPtr_[cvobf_[ibf]];
  }//(ibf < nbf)
}//getCvData_

//two pointers as input. Moves data from the cv pointer to the right face.
void InflowTurbulenceIG::getCvVectorData_(double (*FacePtr)[3], double (*CvPtr)[3]){
  for (int ibf =0; ibf < zone_ptr->nbf; ++ibf) {
    FOR_I3{
      FacePtr[ibf][i]=CvPtr[cvobf_[ibf]][i];
    }//FOR_I3
  }//(ibf < nbf)
}//getCvVectorData()
