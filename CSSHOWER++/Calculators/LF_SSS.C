#include "CSSHOWER++/Showers/Splitting_Function_Base.H"

namespace CSSHOWER {
  
  class LF_SSS: public SF_Lorentz {
  public:
  public:

    inline LF_SSS(const SF_Key &key): SF_Lorentz(key.p_ms)
    {
    }

  };

}

using namespace CSSHOWER;

DECLARE_GETTER(LF_SSS_Getter,"SSS",SF_Lorentz,SF_Key);

SF_Lorentz *LF_SSS_Getter::operator()
  (const Parameter_Type &args) const
{
  return new LF_SSS(args);
}

void LF_SSS_Getter::PrintInfo
(std::ostream &str,const size_t width) const
{
  str<<"sss lorentz functions";
}
