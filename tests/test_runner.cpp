#define CATCH_CONFIG_RUNNER
#include "../src/core/catch1/catch.hpp"

// serial test runner, but should be able to call parallel executables from sysytem()

int main( int argc, char* argv[] )
{
  int result = Catch::Session().run( argc, argv );
  return ( result < 0xff ? result : 0xff );
}
