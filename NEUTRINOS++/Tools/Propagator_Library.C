#include "NEUTRINOS++/Tools/Propagator_Library.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Phys/Flavour.H"
#include <iomanip>

using namespace NEUTRINOS;
using namespace ATOOLS;
using namespace std;

std::ostream & NEUTRINOS::operator<<(std::ostream & s,const prop_type::code & type) {
  if (type==prop_type::none)            s<<setw(18)<<"none";
  if (type==prop_type::massless)        s<<setw(18)<<"massless";
  if (type==prop_type::massive)         s<<setw(18)<<"massive";
  if (type==prop_type::verymassive)     s<<setw(18)<<"verymassive";
  if (type==prop_type::unstable)        s<<setw(18)<<"unstable";
  if (type==prop_type::unknown)         s<<setw(18)<<"unknown";
  return s;
}

/////////////////////////////////////////////////////////////////////////////////////////////
//
/////////////////////////////////////////////////////////////////////////////////////////////
Dummy_Prop::Dummy_Prop(const prop_info & info) :
  Propagator_Base("None", info) {}
Complex Dummy_Prop::Calc(const double & p2) { 
  return Complex(0.0,-1.0); 
}

Massless_Prop::Massless_Prop(const prop_info & info) :
  Propagator_Base("Massless", info) {}
Complex Massless_Prop::Calc(const double & p2) { 
  return Complex(0.0,-1.0 / (p2)); 
}

Massive_Prop::Massive_Prop(const prop_info & info) :
  Propagator_Base("Massive", info),
  mass(info.m_mass) {}
Complex Massive_Prop::Calc(const double & p2) { 
  return Complex(0.0,-1.0 / (p2 - mass*mass)); 
}

VeryMassive_Prop::VeryMassive_Prop(const prop_info & info) :
  Propagator_Base("VeryMassive", info),
  mass(info.m_mass) {}
Complex VeryMassive_Prop::Calc(const double & p2) { 
  if ( mass == 0.0 ) {
    return Complex(0.0,-1.0 / p2);
  } else {
    return Complex(0.0,1.0 / (mass*mass)); 
  }
}

Unstable_Prop::Unstable_Prop(const prop_info & info) :
  Propagator_Base("Unstable", info),
  mass(info.m_mass), width(info.m_width) {}
Complex Unstable_Prop::Calc(const double & p2) { 
  double X = p2 - mass*mass;
  double Y = mass*width;
  return Complex(-Y/(X*X+Y*Y),-X/(X*X+Y*Y)); 
}