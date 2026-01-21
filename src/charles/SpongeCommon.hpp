#ifndef SPONGECOMMON_HPP
#define SPONGECOMMON_HPP

inline void reportSpongeError(const string msg) {
  CERR( " SPONGE bc syntax error: \n" \
        << msg << "\n" \
        << " SPONGE TYPE  <L_X,L_Y,L_Z,L_RX,L_RZ> RELAX_TIME <t> STRENGTH <s> LENGTH <l> [EPS_P <eps_p>]" );
}


inline void parseSpongeParam(Param* param, SpongeType& sponge_type, double& sponge_strength,
                             double& eps_p_sponge, double& t_relax, double& sponge_length,
                             double& sponge_pressure) { 


  int iarg = 1;
  while ( iarg < param->size()) { 
    string token = param->getString(iarg++);
    if ( token == "RELAX_TIME") { 
      t_relax = param->getDouble(iarg++);
    } else if ( token == "STRENGTH") { 
      sponge_strength = param->getDouble(iarg++);
    } else if ( token == "TYPE" ) { 
      token = param->getString(iarg++);
      if ( token == "L_X") { 
        sponge_type = SPONGE_LX;
      } else if ( token == "L_Y") {
        sponge_type = SPONGE_LY;
      } else if ( token == "L_Z") {
        sponge_type = SPONGE_LZ;
      } else if ( token == "L_RX") { 
        sponge_type = SPONGE_LRX;
      } else if ( token == "L_RZ") { 
        sponge_type = SPONGE_LRZ; 
      } else { 
        CERR( " > unrecognized sponge type : " << token); 
      }
    } else if ( token == "LENGTH") {
      sponge_length = param->getDouble(iarg++);
    } else if ( token == "EPS_P") { 
      eps_p_sponge = param->getDouble(iarg++);
    } else if ( token == "SPONGE_PRESSURE") { 
      sponge_pressure = param->getDouble(iarg++);
    } else { 
      CERR( " > unrecognized sponge token : " << token);
    }
  }

  // error check the form of the boundary condition... 
  if ( eps_p_sponge < 0.0) { 
    reportSpongeError(" > eps_p : negative value specified"); 
  } else if ( sponge_length <= 0.0) { 
    reportSpongeError(" > LENGTH: nonpositive value specified or unspecified");
  } else if ( sponge_strength < 0.0) { 
    reportSpongeError(" > STRENGTH: negative value specified or unspecified");
  } else if ( t_relax <= 0.0) { 
    reportSpongeError(" > RELAX_TIME: nonpositive value specified or unspecified");
  } else if ( sponge_pressure < 0.0) { 
    reportSpongeError(" > SPONGE_PRESSURE: nonpositive value specified");
  } else if ( sponge_type == SPONGE_UNDEFINED) { 
    reportSpongeError( " > TYPE: unspecified");
  }
}

#endif
