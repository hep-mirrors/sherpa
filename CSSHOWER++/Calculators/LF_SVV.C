#include "CSSHOWER++/Showers/Single_Splitting_Function.H"

namespace CSSHOWER {
  
  class LF_SVV: public SF_Lorentz {
  public:
  public:

    inline LF_SVV(const SF_Key &key)
    {
    }

  };

}

using namespace CSSHOWER;

DECLARE_GETTER(LF_SVV_Getter,"Gab",SF_Lorentz,SF_Key);

SF_Lorentz *LF_SVV_Getter::operator()
  (const Parameter_Type &args) const
{
  return new LF_SVV(args);
}

void LF_SVV_Getter::PrintInfo
(std::ostream &str,const size_t width) const
{
  str<<"svv lorentz functions";
}
