#define CATCH_CONFIG_RUNNER
#include "../../catch1/catch.hpp"
#include "../MpiStuff.hpp"
//#include "../CTI.hpp"

int main( int argc, char* argv[] )
{
  // global setup...
  MPI_Init(&argc, &argv);
  MpiStuff::initMpiStuff();

  // if using CTI.hpp
  //CTI::CTI_Init(argc, argv,"test.in");

  int result = Catch::Session().run( argc, argv );

  // global clean-up...
  MpiStuff::destroyMpiStuff();
  MPI_Finalize();

  // if using CTI.hpp
  //CTI::CTI_Finalize();

  return ( result < 0xff ? result : 0xff );
}
