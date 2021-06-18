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
    str<<"  Parton ["<<ATOOLS::ID(part.m_id)<<"], stat="
       <<part.m_stat<<", kin="<<part.m_kin<<", kscheme="<<part.m_kscheme
       <<", col="<<part.m_col<<" : "<<part.m_flav<<" : "<<part.m_mom
       <<" "<<sqrt(dabs(part.m_mom.Abs2()))<<" "<<sqrt(dabs(part.Mass2()))
       <<" ("<<part.GetFlow(1)<<","<<part.GetFlow(2)<<")"<<endl;
    if (part.m_pst==pst::IS)      str<<"     (Initial state parton)";
    else if (part.m_pst==pst::FS) str<<"     (Final state parton)  ";
    else                     str<<"                           ";
    str<<"  Colour partners ("
       <<(part.p_left?ATOOLS::ID(part.p_left->m_id):vector<int>())<<","
       <<(part.p_right?ATOOLS::ID(part.p_right->m_id):vector<int>())<<")"<<endl;
    if (part.m_kt_soft[0]<std::numeric_limits<double>::max() ||
	part.m_kt_soft[1]<std::numeric_limits<double>::max()) {
      str<<"  k_T left : "<<sqrt(part.KtSoft(0))<<", k_T right : "<<sqrt(part.KtSoft(1))<<endl;
    }
    str<<"  k_T start : "<<sqrt(part.m_kt_start);
    str<<"  k_T test : "<<sqrt(part.m_kt_test);
    str<<"  k_T veto : "<<sqrt(part.m_kt_veto);
    str<<"  x_B : "<<part.m_xBj<<std::endl;
    if (part.p_prev || part.p_next) {
      if (part.p_prev) str<<"  P="<<part.p_prev;
      if (part.p_next) str<<"  N="<<part.p_next;
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
  msg_IODebugging()<<METHOD<<"("<<this<<") {\n";
  if (this==p_sing->GetSplit()) {
    msg_Indent();
    Parton *left(p_sing->GetLeft()), *right(p_sing->GetRight());
    Vec4D pl(left->Momentum()), pr(right->Momentum());
    Poincare oldcms(pl+pr), newcms(m_mom);
    oldcms.Boost(pl);
    oldcms.Boost(pr);
    newcms.BoostBack(pl);
    newcms.BoostBack(pr);
    if (dabs(pl[0])>dabs(pr[0])) {
      left->SetMomentum(m_mom-pr);
      right->SetMomentum(pr);
    }
    else {
      left->SetMomentum(pl);
      right->SetMomentum(m_mom-pl);
    }
    left->UpdateDaughters();
    right->UpdateDaughters();
  }
  if (p_next) {
    msg_Indent();
    if (p_next==p_next->GetSing()->GetSplit() &&
	m_flav!=p_next->GetFlavour())
      THROW(fatal_error,"invalid flavor change");
    p_next->SetMomentum(m_mom);
    p_next->SetFlavour(m_flav);
    p_next->SetMass2(m_t);
    msg_IODebugging()<<*p_next;
    p_next->UpdateDaughters();
  }
  msg_IODebugging()<<"}\n";
}

void Parton::UpdateColours(int newr,int newa)
{
  msg_IODebugging()<<METHOD<<"("<<this<<"): ("
		   <<newr<<","<<newa<<") {\n";
  {
    msg_Indent();
    SetFlow(1,newr);
    SetFlow(2,newa);
    p_left=p_right=NULL;
    int f1(GetFlow(1)), f2(GetFlow(2));
    for (PLiter pit(p_sing->begin());pit!=p_sing->end();++pit) {
      if (f1 && f1==(*pit)->GetFlow(2)) (p_left=*pit)->SetRight(this);
      if (f2 && f2==(*pit)->GetFlow(1)) (p_right=*pit)->SetLeft(this);
    }
    msg_IODebugging()<<*this;
    if (p_next) p_next->UpdateColours(newr,newa);
    if (this==p_sing->GetSplit()) {
      int oldr(GetFlow(1)), olda(GetFlow(2));
      Parton *left(p_sing->GetLeft()), *right(p_sing->GetRight());
      if (left->GetFlow(1)==oldr) {
	if (left->GetFlow(2)==olda)
	  left->UpdateColours(newr,newa);
	else left->UpdateColours(newr,left->GetFlow(2));
      }
      else if (left->GetFlow(2)==olda) {
	left->UpdateColours(left->GetFlow(1),newa);
      }
      if (right->GetFlow(1)==oldr) {
	if (right->GetFlow(2)==olda)
	  right->UpdateColours(newr,newa);
	else right->UpdateColours(newr,right->GetFlow(2));
      }
      else if (right->GetFlow(2)==olda) {
	right->UpdateColours(right->GetFlow(1),newa);
      }
    }
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

void Parton::SetLeftOf(Parton * part)
{
  part->SetLeft(p_left);
  if (p_left) p_left->SetRight(part);
}

void Parton::SetRightOf(Parton * part)
{
  part->SetRight(p_right);
  if (p_right) p_right->SetLeft(part);
}

