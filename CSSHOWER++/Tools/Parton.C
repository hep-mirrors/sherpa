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
    if (part.p_invleft || part.p_invright) {
      str<<"   Inverse connections [";
      if (part.p_invleft) str<<ATOOLS::ID(part.p_invleft->Id());
      else str<<"(--)";
      str<<", ";
      if (part.p_invright) str<<ATOOLS::ID(part.p_invright->Id());
      else str<<"(--)";
      str<<"]"<<endl;
    }
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
      str<<"  fix spec : "<<part.m_fixspec<<", oldp : "<<part.OldMomentum()
	 <<" "<<(IsEqual(part.OldMomentum(),part.Momentum(),1.0e-6)?
		 "ok":"*** ERROR ***")<<"\n";
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

bool Parton::Splits()
{
  if (this==NULL) return false;
  if (this==p_sing->GetSplit()) return true;
  return p_next->Splits();
}

void Parton::UpdateDaughters()
{
  if (this==NULL || p_next==NULL) return;
  msg_Indent();
  msg_IODebugging()<<METHOD<<"("<<this<<") {\n";
  p_next->SetMomentum(m_mom);
  p_next->SetFlavour(m_flav);
  p_next->SetMass2(m_t);
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
  p_next->SetFlavour(m_flav);
  p_next->SetMass2(m_t);
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
  if (this==NULL) return;
  msg_IODebugging()<<METHOD<<"("<<this<<"): ("
		   <<GetMEFlow(1)<<","<<GetMEFlow(2)<<") -> ("
		   <<GetFlow(1)<<","<<GetFlow(2)<<") {\n";
  {
    msg_Indent();
    if (p_sing==NULL) THROW(fatal_error,"Cannot update flow");
    p_left=p_right=NULL;
    int f1(GetFlow(1)), f2(GetFlow(2));
    for (PLiter pit(p_sing->begin());pit!=p_sing->end();++pit) {
      if (f1 && f1==(*pit)->GetFlow(2)) (p_left=*pit)->SetRight(this);
      if (f2 && f2==(*pit)->GetFlow(1)) (p_right=*pit)->SetLeft(this);
    }
    msg_IODebugging()<<*this;
    if (this==p_sing->GetSplit() ||
	(p_prev && p_prev==p_sing->GetSplit())) {
      Parton *ds[2]={p_sing->GetLeft(),p_sing->GetRight()};
      for (int i(0);i<2;++i)
	for (int n(1);n<=2;++n)
	  if (ds[i]->GetFlow(n)==GetMEFlow(n)) {
	    ds[i]->SetFlow(n,GetFlow(n));
	    ds[i]->UpdateColours();
	  }
    }
    else if (p_next) {
      for (int n(1);n<=2;++n) p_next->SetFlow(n,GetFlow(n));
      p_next->UpdateColours();
    }
    for (int n(1);n<=2;++n) SetMEFlow(n,GetFlow(n));
  }
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

