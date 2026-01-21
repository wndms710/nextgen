#include "../src/core/catch1/catch.hpp"
#include <stdlib.h>  // for system()
#include <sys/stat.h>  // for stat()
#include <string>

// fast snippet to check for file accessibility, i.e., existence
inline bool file_exists (const std::string& name) {
  struct stat buffer;   
  return (stat (name.c_str(), &buffer) == 0); 
}


TEST_CASE( "check for command processor","[runner]" ) {
  REQUIRE(system(NULL));  // otherwise system calls not available

  // success, let's clean up any potential files we will create (non-logs)
  if (file_exists("sector.periodic.sbin")) system("rm sector.periodic.sbin");
  if (file_exists("sector.periodic.mles")) system("rm sector.periodic.mles");
  if (file_exists("result.sles")) system("rm result.sles");
}

// use surfer to build periodic surface
// then stitch to create a simple mesh file
TEST_CASE( "sector: periodic cyl_x","[runner]" ) {

  // return value from system calls is the same as the executable's return value; so 0 on success...
  SECTION( "surfer: add periodicity" ) {
    REQUIRE(system("../bin/surfer.exe -i assets/sector_per_cyl_x/surfer.in 2>&1 > log.surfer.tmp") == 0);
    REQUIRE(file_exists("log.surfer.tmp"));
    REQUIRE(file_exists("sector.periodic.sbin"));
  }

  SECTION( "stitch: generate a periodic mesh" ) {
    REQUIRE(system("../bin/stitch.exe -i assets/sector_per_cyl_x/stitch.in 2>&1 > log.stitch.tmp") == 0);
    REQUIRE(file_exists("log.stitch.tmp"));
    REQUIRE(file_exists("sector.periodic.mles"));
  }

  SECTION( "charles: run a zero velocity, all walls domain" ) {
    REQUIRE(system("../bin/charles.exe -i assets/sector_per_cyl_x/charles.in 2>&1 > log.charles.tmp") == 0);
    REQUIRE(file_exists("log.charles.tmp"));
  }

  SECTION("cleanup") {
    // if we got here we can simply cleanup
    system("rm -f sector.periodic.sbin sector.periodic.mles result.sles log.*.tmp");
    system("rm -rf images_tmp");
  }
}

TEST_CASE( "sector with orphans: periodic cyl_x","[runner][orphan]" ) {

  // return value from system calls is the same as the executable's return value; so 0 on success...
  SECTION( "surfer: add periodicity" ) {
    REQUIRE(system("../bin/surfer.exe -i assets/sector_per_cyl_x/surfer.sheet.in 2>&1 > log.surfer.tmp") == 0);
    REQUIRE(file_exists("log.surfer.tmp"));
    REQUIRE(file_exists("sector.periodic.sbin"));
  }

  SECTION( "stitch: generate a periodic mesh" ) {
    REQUIRE(system("../bin/stitch.exe -i assets/sector_per_cyl_x/stitch.in 2>&1 > log.stitch.tmp") == 0);
    REQUIRE(file_exists("log.stitch.tmp"));
    REQUIRE(file_exists("sector.periodic.mles"));
  }

  SECTION( "charles: run a zero velocity, all walls domain" ) {
    REQUIRE(system("../bin/charles.exe -i assets/sector_per_cyl_x/charles.in 2>&1 > log.charles.tmp") == 0);
    REQUIRE(file_exists("log.charles.tmp"));
  }

  SECTION("cleanup") {
    // if we got here we can simply cleanup
    system("rm -f sector.periodic.sbin sector.periodic.mles result.sles log.*.tmp");
    system("rm -rf images_tmp");
  }
}
