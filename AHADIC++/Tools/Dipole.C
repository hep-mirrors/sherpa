#include "Dipole.H"

using namespace AHADIC;
using namespace ATOOLS;

Dipole::Dipole(SP(Proto_Particle) trip,SP(Proto_Particle) anti) 
{
  ++s_cnt;
  p_triplet    =trip;
  p_antitriplet=anti;
  m_mustdecay=(p_triplet->m_flav.IsGluon() || p_antitriplet->m_flav.IsGluon()); 
  m_mass2   =(p_triplet->m_mom+p_antitriplet->m_mom).Abs2();
  m_massbar2=ATOOLS::sqr(sqrt(m_mass2)-(p_triplet->m_mass+p_antitriplet->m_mass));
}

void Dipole::Update() 
{
  if (p_triplet!=NULL && p_antitriplet!=NULL) {
    m_mustdecay=(p_triplet->m_flav.IsGluon() || p_antitriplet->m_flav.IsGluon());
    m_mass2   =(p_triplet->m_mom+p_antitriplet->m_mom).Abs2();
    m_massbar2= ATOOLS::sqr(sqrt(m_mass2)-(p_triplet->m_mass+p_antitriplet->m_mass));
  }
}

void Dipole::Output() 
{
  //msg_Out()<<"--- Dipole["<<this<<"] ("<<p_triplet<<" "<<p_antitriplet<<" : ";
  msg_Out()<<"--- Dipole[";
  if (p_triplet!=NULL) msg_Out()<<p_triplet->m_flav;
  else msg_Out()<<" no flav ";
  msg_Out()<<", ";
  if (p_antitriplet!=NULL) msg_Out()<<p_antitriplet->m_flav;
  else msg_Out()<<" no flav ";
  msg_Out()<<"], "<<sqrt(m_mass2)<<", decay = "<<m_mustdecay<<") ---"<<std::endl<<"--- ";
  if (p_triplet!=NULL)	
    msg_Out()<<p_triplet->m_mom<<" "
	     <<sqrt(ATOOLS::Max(p_triplet->m_mom.Abs2(),0.));
  else msg_Out()<<" XXX ";
  msg_Out()<<" + "; 
  if (p_antitriplet!=NULL)	
    msg_Out()<<p_antitriplet->m_mom<<" "
	     <<sqrt(ATOOLS::Max(p_antitriplet->m_mom.Abs2(),0.));
  else msg_Out()<<" XXX ";
  msg_Out()<<" ---"<<std::endl; 
}
