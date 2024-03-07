#include "ALARIC++/Tools/Splitting.H"

#include "ALARIC++/Tools/Amplitude.H"
#include "ATOOLS/Org/STL_Tools.H"
#include "ATOOLS/Org/Message.H"

using namespace ALARIC;
using namespace ATOOLS;

void Splitting::SetType()
{
  m_type=p_c->Beam()?1:0;
}

namespace ALARIC {

  std::ostream &operator<<(std::ostream &s,const Splitting &p)
  {
    s<<"["<<(p.p_c?p.p_c->Id():0)<<"<->"<<(p.p_s?p.p_s->Id():0)
     <<"](c="<<p.m_cm<<",kfac="<<p.m_kfac
     <<",soft="<<p.m_soft<<",mec="<<p.m_mec
     <<"){t="<<p.m_t<<",z="<<p.m_z<<",phi="<<p.m_phi
     <<",s2="<<p.m_s<<",z2="<<p.m_z2<<",phi2="<<p.m_phi2<<"|";
    if (p.m_ci.size()) s<<p.m_ci[0]<<p.m_cj[0];
    for (size_t i(1);i<p.m_ci.size();++i) s<<";"<<p.m_ci[i]<<p.m_cj[i];
    return s<<"|"<<p.m_h[0]<<","<<p.m_h[1]<<","<<p.m_h[2]<<"}"<<p.m_rcl;
  }

}
