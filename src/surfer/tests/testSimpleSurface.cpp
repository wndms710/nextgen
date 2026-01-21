#include "../SimpleSurface.hpp"
#include "../../core/catch1/catch.hpp"

TEST_CASE("SimpleSurface","[surfer]") {

  SimpleSurface ss;

  Param p("SURF");
  p.tokens.push_back("SIMPLE_BOX");
  p.tokens.push_back("0");
  p.tokens.push_back("1");
  p.tokens.push_back("0");
  p.tokens.push_back("1");
  p.tokens.push_back("0");
  p.tokens.push_back("1");
  
  ss.init(&p);

  double area,vol;
  ss.getAreaAndVolume(vol,area);
  CHECK(area == 6.0);
  CHECK(vol == 1.0);
}
