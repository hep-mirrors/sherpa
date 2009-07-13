#include "CSSHOWER++/Tools/Parton.H"
#include "ATOOLS/Phys/Cluster_Leg.H"

using namespace CSSHOWER;
using namespace std;

namespace CSSHOWER {
  std::ostream& operator<<(std::ostream& str, const Parton &part) {
    str<<"  Parton "<<&part<<" ("<<part.m_stat<<")["<<ATOOLS::ID(part.m_id)
       <<"]: "<<part.m_flav<<" : "<<part.m_mom
       <<" ("<<part.GetFlow(1)<<","<<part.GetFlow(2)<<")"
       <<"["<<part.GetRFlow(1)<<","<<part.GetRFlow(2)<<"]"<<endl;
    if (part.m_pst==pst::IS)      str<<"     (Initial state parton)";
    else if (part.m_pst==pst::FS) str<<"     (Final state parton)  ";
    else                     str<<"                           ";
    str<<"  Colour partners ("
       <<part.p_left<<","<<part.p_right<<")"<<endl;
    str<<"  k_T start : "<<sqrt(part.m_kt_start);
    str<<"  k_T test : "<<sqrt(part.m_kt_test);
    str<<"  k_T veto : "<<sqrt(part.m_kt_veto)<<"("<<sqrt(part.m_kt_max)<<")";
    str<<"  x_B : "<<part.m_xBj<<std::endl;
    if (part.p_prev || part.p_next || 
	part.m_kt_prev!=0.0 || part.m_kt_next!=0.0) {
      if (part.p_prev) str<<"  P="<<part.p_prev;
      if (part.p_next) str<<"  N="<<part.p_next;
      str<<"  k_T prev : "<<sqrt(part.m_kt_prev);
      str<<"  k_T next : "<<sqrt(part.m_kt_next);
      str<<std::endl;
    }
    return str;
  }
}

Parton *Parton::FollowUp()
{
  if (p_next) return p_next->FollowUp();
  return this;
}

void Parton::UpdateDaughters()
{
  if (this==NULL || p_next==NULL) return;
  msg_Indent();
  msg_Debugging()<<METHOD<<"("<<this<<") {\n";
  p_next->SetMomentum(m_mom);
  msg_Debugging()<<*p_next;
  p_next->UpdateDaughters();
  msg_Debugging()<<"}\n";
}

