#include "ATOOLS/Phys/Simple_Polarisation_Info.H"

namespace ATOOLS {

std::ostream& operator<< (std::ostream& s, const Simple_Polarisation_Info & polinfo)
{
  return s<<"("<<polinfo.GetPol()<<","<<polinfo.Info()<<","<<polinfo.Angle()<<")";
}

}
