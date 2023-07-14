#include "AHADIC++/Tools/Z_Selector.H"
#include "AHADIC++/Tools/Splitter_Base.H"
#include "AHADIC++/Tools/Hadronisation_Parameters.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Message.H"

using namespace AHADIC;
using namespace ATOOLS;


Z_Selector::Z_Selector() {}

void Z_Selector::Init(Splitter_Base * splitterbase) {
  p_splitterbase = splitterbase;
}

Z_Selector::~Z_Selector() {}

double Z_Selector::operator()(const double & zmin,const double & zmax,
			      const unsigned int & cnt) {
  if (p_splitterbase==NULL) return zmin+ran->Get()*(zmax-zmin);
  double z(-1.);
  do {
    z = zmin+ran->Get()*(zmax-zmin);
    auto wgt = p_splitterbase->WeightFunction(z,zmin,zmax,cnt);
    if((wgt<ran->Get())) {
      // add acceptance weight
      break;
    }
    else {
      // add rejection weight
      continue;
    }
  } while (true);
  return z;
}


