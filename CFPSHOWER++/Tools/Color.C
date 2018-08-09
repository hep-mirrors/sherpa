#include "CFPSHOWER++/Tools/Color.H"

using namespace CFPSHOWER;
using namespace std;

std::ostream & CFPSHOWER::operator<<(std::ostream & s,const Color & color) {
  s<<"{"<<color[0]<<", "<<color[1]<<"}";
  return s;
}
