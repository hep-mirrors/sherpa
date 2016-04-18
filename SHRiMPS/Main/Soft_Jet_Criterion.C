#include "SHRiMPS/Main/Soft_Jet_Criterion.H"
#include "ATOOLS/Org/Message.H"

using namespace SHRIMPS;
using namespace ATOOLS;
using namespace std;

double Soft_Jet_Criterion::PTij2(ATOOLS::Vec4D pi,ATOOLS::Vec4D pj) {
  double pti2  = pi.PPerp2();
  double ptj2  = pj.PPerp2();    
  return Min(pti2,ptj2)*(cosh(pi.Eta()-pj.Eta())-cos(pi.Phi()-pj.Phi()));
}
  
size_t Soft_Jet_Criterion::FindCombination(ATOOLS::Cluster_Amplitude *ampl) {
  ClusterLeg_Vector newlegs(ampl->Legs());
  if (newlegs.size()<m_reflegs.size()+1) return 0;
  msg_Out()<<METHOD<<".\n";
  size_t  winner(0),comp;
  Flavour lastflav(newlegs[newlegs.size()-1]->Flav().Bar());
  Vec4D   lastmom(newlegs[newlegs.size()-1]->Mom());
  for (size_t i=2;i<newlegs.size()-1;i++) {
    // only combinable flavours
    if (newlegs[i]->Flav()!=lastflav) continue;
    // colours to project on reference partons
    size_t cinew = newlegs[i]->Col().m_i; 
    size_t cjnew = newlegs[i]->Col().m_j; 
    bool   hit   = false;
    size_t maxj  = Min(i,m_reflegs.size()-1);
    for (size_t j=2;j<=maxj;j++) {
      if (m_reflegs[j]->Col().m_i==cinew && m_reflegs[j]->Col().m_j==cjnew) {
	double pt2 = PTij2(lastmom,newlegs[i]->Mom());
	if (pt2<m_pt2) {
	  m_pt2  = pt2;
	  winner = j;
	  comp   = i;
	  break;
	}
      }
    }
  }
  if (winner==0) {
    msg_Out()<<"Uncombinable step.\n";
  }
  else {
    msg_Out()<<"Winning pair: "<<(newlegs.size()-1)<<" + "<<comp<<", "
	     <<"pt^2 = "<<m_pt2<<" vs. "<<m_kt2veto[m_reflegs[winner]]<<" from "
	     <<(lastmom+newlegs[comp]->Mom())<<".\n";
  }
  return winner;
}

bool Soft_Jet_Criterion::Jets(ATOOLS::Cluster_Amplitude *ampl,int mode) {
  m_pt2  = 1.e12;
  ClusterLeg_Vector newlegs(ampl->Legs());
  for (size_t i=2;i<newlegs.size();i++) {
    msg_Out()<<"--> "<<i<<" "<<newlegs[i]->Flav()<<" "
	     <<"["<<newlegs[i]->Col().m_i<<" "<<newlegs[i]->Col().m_j<<"], "
	     <<newlegs[i]->Mom()<<".\n";
  }
  msg_Out()<<"\n";
  for (size_t i=2;i<m_reflegs.size();i++) {
    msg_Out()<<"--> "<<i<<" "<<m_reflegs[i]->Flav()<<" "
	     <<"["<<m_reflegs[i]->Col().m_i<<" "<<m_reflegs[i]->Col().m_j<<"], "
	     <<m_reflegs[i]->Mom()<<".\n";
  }
  size_t winner(FindCombination(ampl));
  if (m_kt2veto.find(m_reflegs[winner])==m_kt2veto.end()) {
    msg_Out()<<" reconstructed non-core parton.\n";
  }
  else if (winner>0 && m_pt2>m_kt2veto[m_reflegs[winner]]) {
    msg_Out()<<METHOD<<": "<<winner<<" yields pt^2 = "<<m_pt2
	     <<" vs. "<<m_kt2veto[m_reflegs[winner]]<<".\n";
    Output();
    return true;
  }
  return false;
}

void Soft_Jet_Criterion::Output() {
  msg_Out()<<"------------------------------------------------------\n"
	   <<METHOD<<" for reference ampl ["<<p_refampl<<"], "
	   <<m_kt2veto.size()<<" entries:\n";
  for (map<Cluster_Leg *,double>::iterator vetit=m_kt2veto.begin();
       vetit!=m_kt2veto.end();vetit++) {
    msg_Out()<<vetit->first->Id()<<": ktveto = "<<sqrt(vetit->second)<<"\n";
  }
  msg_Out()<<"------------------------------------------------------\n";
}

JF::JF() : PHASIC::Jet_Finder(), m_ycut(1.) {}
