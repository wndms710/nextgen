#define CATCH_CONFIG_RUNNER
#include "../../catch1/catch.hpp"
#include "../CTI.hpp"

int main( int argc, char* argv[] )
{
  // global setup...
  CTI::CTI_Init(argc, argv,"tests/assets/testExchanger.in");

  int result = Catch::Session().run( argc, argv );

  // global clean-up...
  CTI::CTI_Finalize();

  return ( result < 0xff ? result : 0xff );
}
