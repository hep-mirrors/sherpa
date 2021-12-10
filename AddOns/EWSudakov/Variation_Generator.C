#include "AddOns/EWSudakov/Variation_Generator.H"

#include <ostream>

using namespace ATOOLS;
using namespace EWSud;

using Base = Hard_Process_Variation_Generator_Base;
using Args = Hard_Process_Variation_Generator_Arguments;

Variation_Generator::Variation_Generator(const Args& args):
  Base {args}
{}

DECLARE_GETTER(Variation_Generator, "EWSud", Base, Args);

Base* ATOOLS::Getter<Base, Args, Variation_Generator>::
operator()(const Args& args) const
{
  return new Variation_Generator(args);
}

void ATOOLS::Getter<Base, Args, Variation_Generator>::
PrintInfo(std::ostream& str, const size_t width) const
{ 
  str << "Generator for EWSud(akov) variations"; 
}

