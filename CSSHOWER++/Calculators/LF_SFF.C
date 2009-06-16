#include "CSSHOWER++/Showers/Single_Splitting_Function.H"

namespace CSSHOWER {
  
  class LF_SFF: public SF_Lorentz {
  public:
  public:

    inline LF_SFF(const SF_Key &key)
    {
    }

  };

}

using namespace CSSHOWER;

DECLARE_GETTER(LF_SFF_Getter,"FFS",SF_Lorentz,SF_Key);

SF_Lorentz *LF_SFF_Getter::operator()
  (const Parameter_Type &args) const
{
  return new LF_SFF(args);
}

void LF_SFF_Getter::PrintInfo
(std::ostream &str,const size_t width) const
{
  str<<"sff lorentz functions";
}
