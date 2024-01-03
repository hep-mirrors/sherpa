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
  double z, z_range {zmax-zmin};
  do {
    z = zmin+ran->Get()*z_range;
    auto sel_wgt = p_splitterbase->WeightFunction(z,zmin,zmax,cnt);
    if((sel_wgt<ran->Get())) {
      p_splitterbase->z_rejected(sel_wgt, z,zmin,zmax,cnt);
      continue;
    }
    else {
      // add rejection weight
      p_splitterbase->z_accepted(sel_wgt, z,zmin,zmax,cnt);
      break;
    }
  } while (true);
  return z;
}


