#include "CSSHOWER++/Tools/Parton.H"
#include "CSSHOWER++/Tools/Singlet.H"
#include "ATOOLS/Phys/Cluster_Leg.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Exception.H"

using namespace CSSHOWER;
using namespace ATOOLS;
using namespace std;

namespace CSSHOWER {
  std::ostream& operator<<(std::ostream& str, const Parton &part) {
    str<<"  Parton "<<&part<<" ("<<part.p_sing<<"), stat="
       <<part.m_stat<<", kin="<<part.m_kin<<", kscheme="<<part.m_kscheme
       <<", col="<<part.m_col<<" ["<<ATOOLS::ID(part.m_id)
       <<"]: "<<part.m_flav<<" : "<<part.m_mom
       <<" "<<sqrt(dabs(part.m_mom.Abs2()))<<" "<<sqrt(dabs(part.Mass2()))
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
    if (part.p_prev || part.p_next) {
      if (part.p_prev) str<<"  P="<<part.p_prev;
      if (part.p_next) str<<"  N="<<part.p_next;
      str<<std::endl;
    }
    if (part.m_fixspec!=Vec4D())
      str<<"  fix spec : "<<part.m_fixspec<<"\n";
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
  msg_IODebugging()<<METHOD<<"("<<this<<") {\n";
  p_next->SetMomentum(m_mom);
  msg_IODebugging()<<*p_next;
  p_next->UpdateDaughters();
  msg_IODebugging()<<"}\n";
}

void Parton::UpdateNewDaughters(Parton *ref)
{
  if (this==NULL || p_next==NULL) return;
  ref=ref->GetSing()->GetLeft();
  if (ref==NULL) THROW(fatal_error,"Internal error");
  msg_Indent();
  msg_IODebugging()<<METHOD<<"("<<this<<") {\n";
  p_next->SetMomentum(m_mom);
  for (int n(1);n<=2;++n) {
    p_next->SetFlow(n,GetFlow(n));
    p_next->SetMEFlow(n,GetMEFlow(n));
  }
  p_next->SetStart(m_kt_start);
  p_next->SetKtMax(m_kt_max);
  p_next->SetVeto(m_kt_veto);
  p_next->SetId(m_id);
  msg_IODebugging()<<*p_next;
  p_next->UpdateNewDaughters(ref);
  msg_IODebugging()<<"}\n";
}

void Parton::UpdateColours()
{
  if (this==NULL || p_next==NULL) return;
  msg_Indent();
  msg_IODebugging()<<METHOD<<"("<<this<<") {\n";
  for (int n(1);n<=2;++n) {
    p_next->SetFlow(n,GetFlow(n));
    p_next->SetMEFlow(n,GetMEFlow(n));
  }
  msg_IODebugging()<<*p_next;
  p_next->UpdateColours();
  msg_IODebugging()<<"}\n";
}

double Parton::Weight(const double &scale)
{
  double weight=1.0;
  for (size_t i(0);i<m_weights.size();++i)
    if (m_weights[i].first>scale) weight*=m_weights[i].second;
    else break;
  return weight;
}

void Parton::SetLeftOf(Parton * part) {
  //msg_Out()<<"   "<<METHOD
  //	   <<" for this ["<<m_flow.Code(1)<<", "<<m_flow.Code(2)<<"], "
  //	   <<"part ["<<part->GetFlow(1)<<", "<<part->GetFlow(2)<<"]";
  if (!m_connected) {
    //msg_Out()<<" --> NOT CONNECTED.\n";
    return;
  }
  //if (p_left) {
  //msg_Out()<<", left ["<<p_left->GetFlow(1)<<", "<<p_left->GetFlow(2)<<"]";
  //}
  //msg_Out()<<".\n";
  part->SetLeft(p_left);
  if (p_left) p_left->SetRight(part);
}

void Parton::SetRightOf(Parton * part) {
  //msg_Out()<<"   "<<METHOD
  //	   <<" for this ["<<m_flow.Code(1)<<", "<<m_flow.Code(2)<<"], "
  //	   <<"part ["<<part->GetFlow(1)<<", "<<part->GetFlow(2)<<"]";
  if (!m_connected) {
    //msg_Out()<<" --> NOT CONNECTED.\n";
    return;
  }
  //if (p_right) {
  //  msg_Out()<<", right ["<<p_right->GetFlow(1)<<", "
  //<<p_right->GetFlow(2)<<"]";
  //}
  //msg_Out()<<".\n";
  part->SetRight(p_right);
  if (p_right) p_right->SetLeft(part);
}

