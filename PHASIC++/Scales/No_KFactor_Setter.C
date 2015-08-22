#include "PHASIC++/Scales/KFactor_Setter_Base.H"

namespace PHASIC {

  class No_KFactor_Setter: public KFactor_Setter_Base {
  public:

    No_KFactor_Setter(const KFactor_Setter_Arguments &args);

    double KFactor(const int mode=0);
    double KFactor(const ATOOLS::NLO_subevt& evt);

  };// end of class No_KFactor_Setter

}// end of namespace PHASIC

using namespace PHASIC;
using namespace ATOOLS;

DECLARE_GETTER(No_KFactor_Setter,"NO",
	       KFactor_Setter_Base,KFactor_Setter_Arguments);

KFactor_Setter_Base *ATOOLS::Getter
<KFactor_Setter_Base,KFactor_Setter_Arguments,No_KFactor_Setter>::
operator()(const KFactor_Setter_Arguments &args) const
{
  return new No_KFactor_Setter(args);
}

void ATOOLS::Getter<KFactor_Setter_Base,KFactor_Setter_Arguments,
		    No_KFactor_Setter>::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"no K-factor scheme\n";
}

No_KFactor_Setter::No_KFactor_Setter
(const KFactor_Setter_Arguments &args):
  KFactor_Setter_Base(args) {}

double No_KFactor_Setter::KFactor(const int mode) 
{
  return m_weight=1.0;
}

double No_KFactor_Setter::KFactor(const ATOOLS::NLO_subevt& evt) 
{
  return KFactor();
}

