#include "../StripedMesh.hpp"
#include "../../catch1/catch.hpp"

// build extended faces test

TEST_CASE( "striped mesh","[core]" ) {
  StripedMesh sm;

  SECTION("file does not exist") {
    CHECK_THROWS(sm.read_mles("restart.does.not.exist.mles"));
  }
  SECTION("file does exist") {
    CHECK_NOTHROW(sm.read_mles("tests/assets/orifice.mles"));

    CHECK(sm.nno_global == 46348);
    CHECK(sm.nbf_global == 693);
    CHECK(sm.nfa_global == 5654);
    CHECK(sm.nfa_p_global == 0);
    CHECK(sm.ncv_global == 1056);
    CHECK(sm.nno_p_global == 0);
    CHECK(sm.nno_pb_global == 42383);
    CHECK(sm.nef_global == 0);
  }
}
