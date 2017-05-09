#include "AHADIC++/Tools/Z_Selector.H"
#include "AHADIC++/Tools/Hadronisation_Parameters.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Message.H"

using namespace AHADIC;
using namespace ATOOLS;


Z_Selector::Z_Selector() {}

void Z_Selector::Init() {}

Z_Selector::~Z_Selector() {}

double Z_Selector::operator()(const double & zmin,const double & zmax) {
  double z(-1.);
  do {
    z = zmin+ran->Get()*(zmax-zmin);
  } while (WeightFunction(z)<ran->Get());
  return z;
}

double Z_Selector::WeightFunction(const double & z) {
  return 2.*(sqr(z)+sqr(1.-z));
}

