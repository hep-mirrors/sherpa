#include "NEUTRINOS++/Tools/Propagator_Maps.H"
#include "ATOOLS/Phys/Flavour_Tags.H"
#include "ATOOLS/Org/Message.H"
#include <iomanip>

namespace NEUTRINOS {
  Propagator_Maps * ffprops(NULL);
}

using namespace NEUTRINOS;
using namespace ATOOLS;
using namespace std;

Propagator_Maps::Propagator_Maps() {}

Propagator_Maps::~Propagator_Maps() {}

Propagator_Base * 
Propagator_Maps::GetProp(kf_code & prop,prop_type::code & prop_type) {
  double mass = Flavour(prop).Mass();
  double width = Flavour(prop).Width();
  msg_Out() << mass << " " << width << "\n";
  prop_info info = prop_info(prop_type, mass, width);
  switch (prop_type) {
  case prop_type::none:
    return new Dummy_Prop(info);
  case prop_type::massless:
    return new Massless_Prop(info);
  case prop_type::massive:
    return new Massive_Prop(info);
  case prop_type::unstable:
    return new Unstable_Prop(info);

  case prop_type::unknown:
    msg_Error()<<"Error in "<<METHOD<<": found unknown propagator type.\n"
    <<"    Will return Dummy propagator for "
    <<Flavour(prop)<<"\n";
    return new Dummy_Prop(info);
  }
  return new Unstable_Prop(info);
}

void Propagator_Maps::Initialize() {
}