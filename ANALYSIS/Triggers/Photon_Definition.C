#include "Photon_Definition.H"
#include "MyStrStream.H"
#include "Particle_Qualifier.H"
#include <iomanip>

using namespace ANALYSIS;

DECLARE_GETTER(Photon_Definition_Getter,"Photon",
 	       Object_Definition_Base,Argument_Matrix);

void Photon_Definition_Getter::PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"{\n"
     <<std::setw(width+7)<<" "<<"blabla\n"
     <<std::setw(width+4)<<" "<<"}";
}

Object_Definition_Base * 
Photon_Definition_Getter::operator()(const Argument_Matrix &parameters) const
{
}
