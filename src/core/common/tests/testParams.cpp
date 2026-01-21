#include "../../catch1/catch.hpp"
#include "../Params.hpp"

TEST_CASE("param initialization","[params]") {
  Param p;
  p.setName("hello");
  p.tokens.push_back("world");
  p.tokens.push_back("we");
  p.tokens.push_back("hardly");
  p.tokens.push_back("knew");
  p.tokens.push_back("ye");

  CHECK(p.getName() == "hello");
  CHECK(p.str() == "hello world we hardly knew ye");
  CHECK(p.size() == 5);
}
