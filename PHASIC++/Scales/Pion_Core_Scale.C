#include "PHASIC++/Scales/Scale_Setter_Base.H"

#include "MODEL/Main/Running_AlphaS.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"

namespace PHASIC {

  class Pion_Core_Scale : public Scale_Setter_Base {
  public:

    Pion_Core_Scale(const PHASIC::Scale_Setter_Arguments &args);

    double Calculate(const ATOOLS::Vec4D_Vector &p,
           const size_t &mode);

  };// end of class Scale_Setter_Base

}// end of namespace PHASIC

using namespace PHASIC;
using namespace ATOOLS;

Pion_Core_Scale::Pion_Core_Scale(const PHASIC::Scale_Setter_Arguments &args):
Scale_Setter_Base(args)
{

}

double Pion_Core_Scale::Calculate(const ATOOLS::Vec4D_Vector &p,
           const size_t &mode)
{
  if(p.size()!=4){
    msg_Error()<<"PionForm is for e+e- -> pi+pi- only"<<std::endl;
  }
  double Q2 = (p[0]+p[1]).Abs2();
  msg_Debugging()<<"Pion Scale Q2: "<<Q2<<std::endl;
  return Q2;
}

DECLARE_GETTER(Pion_Core_Scale,"Pion",
          Scale_Setter_Base,Scale_Setter_Arguments);

Scale_Setter_Base *ATOOLS::Getter
<Scale_Setter_Base,Scale_Setter_Arguments,Pion_Core_Scale>::
operator()(const Scale_Setter_Arguments &args) const
{
  return new Pion_Core_Scale(args);
}

void ATOOLS::Getter<Scale_Setter_Base,Scale_Setter_Arguments,
            Pion_Core_Scale>::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"Pion core scale"; 
}
