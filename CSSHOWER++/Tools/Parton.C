#include "CSSHOWER++/Tools/Parton.H"
#include "ATOOLS/Phys/Cluster_Leg.H"
#include "ATOOLS/Math/MathTools.H"

using namespace CSSHOWER;
using namespace ATOOLS;
using namespace std;

namespace CSSHOWER {
  std::ostream& operator<<(std::ostream& str, const Parton &part) {
    str<<"  Parton "<<&part<<" ("<<part.m_stat<<"|"
       <<part.m_kin<<")["<<ATOOLS::ID(part.m_id)
       <<"]: "<<part.m_flav<<" : "<<part.m_mom
       <<" ("<<part.GetFlow(1)<<","<<part.GetFlow(2)<<")"
       <<"["<<part.GetRFlow(1)<<","<<part.GetRFlow(2)<<"]"<<endl;
    if (part.m_pst==pst::IS)      str<<"     (Initial state parton)";
    else if (part.m_pst==pst::FS) str<<"     (Final state parton)  ";
    else                     str<<"                           ";
    str<<"  Colour partners ("
       <<part.p_left<<","<<part.p_right<<"), t_min = "
       <<sqrt(dabs(part.m_t_min))<<endl;
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

void Parton::DeleteAll()
{
  if (p_next) p_next->DeleteAll();
  delete this;
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

void Parton::UpdateNewDaughters()
{
  if (this==NULL || p_next==NULL) return;
  msg_Indent();
  msg_Debugging()<<METHOD<<"("<<this<<") {\n";
  p_next->SetMomentum(m_mom);
  for (int n(1);n<=2;++n) {
    p_next->SetFlow(n,GetFlow(n));
    p_next->SetMEFlow(n,GetMEFlow(n));
  }
  p_next->SetStart(m_kt_start);
  p_next->SetKtMax(m_kt_max);
  p_next->SetVeto(m_kt_veto);
  p_next->SetKtPrev(m_kt_prev);
  p_next->SetKtPrev(m_kin);
  msg_Debugging()<<*p_next;
  p_next->UpdateNewDaughters();
  msg_Debugging()<<"}\n";
}

void Parton::UpdateColours()
{
  if (this==NULL || p_next==NULL) return;
  msg_Indent();
  msg_Debugging()<<METHOD<<"("<<this<<") {\n";
  for (int n(1);n<=2;++n) {
    p_next->SetFlow(n,GetFlow(n));
    p_next->SetMEFlow(n,GetMEFlow(n));
  }
  msg_Debugging()<<*p_next;
  p_next->UpdateColours();
  msg_Debugging()<<"}\n";
}

double Parton::Weight(const double &scale) 
{
  double weight=1.0;
  for (size_t i(0);i<m_weights.size();++i)
    if (m_weights[i].first>scale) weight*=m_weights[i].second;
    else break;
  return weight;
}
    
