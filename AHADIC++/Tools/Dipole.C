#include "AHADIC++/Tools/Dipole.H"

using namespace AHADIC;
using namespace ATOOLS;


namespace AHADIC {
  std::list<Dipole *>     Dipole::s_actives;
  std::list<DipoleList *> DipoleList::s_actives;
}


Dipole::Dipole(Proto_Particle * trip,Proto_Particle * anti) 
{
  ++s_cnt;
  p_triplet    =trip;
  p_antitriplet=anti;
  m_mustdecay=(p_triplet->m_flav.IsGluon() || p_antitriplet->m_flav.IsGluon()); 
  m_mass2   =(p_triplet->m_mom+p_antitriplet->m_mom).Abs2();
  m_massbar2=ATOOLS::sqr(sqrt(m_mass2)-(p_triplet->m_mass+p_antitriplet->m_mass));

  s_actives.push_back(this);
}

Dipole::~Dipole() 
{ 
  --s_cnt;
  s_actives.remove(this);
}

bool Dipole::CheckConsistency(std::ostream & s,std::string method) {
  bool passed(dabs(m_mass2-Momentum().Abs2())<1.e-8);
  if (!passed) {
    s<<"Error in "<<METHOD<<" called by "<<method<<":"<<std::endl
     <<"   Masses and momenta not consistent for dipole "
     <<"("<<p_triplet->m_flav<<", "<<p_antitriplet->m_flav<<", "
     <<"mass^2 = "<<m_mass2<<" vs. "<<Momentum()<<" ("<<Momentum().Abs2()<<")"<<std::endl;
  }
  if (p_triplet)     passed = passed && p_triplet->CheckConsistency(s,method);
  if (p_antitriplet) passed = passed && p_antitriplet->CheckConsistency(s,method);
  return passed;   
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
  msg_Out()<<"], (mass = "<<sqrt(m_mass2)<<", decay = "<<m_mustdecay<<") ---"<<std::endl<<"--- ";
  if (p_triplet!=NULL)	
    msg_Out()<<p_triplet->m_mom<<" "
	     <<sqrt(ATOOLS::Max(p_triplet->m_mom.Abs2(),0.))
	     <<" ("<<p_triplet->m_mom.Abs2()<<") ";
  else msg_Out()<<" XXX ";
  msg_Out()<<" + "; 
  if (p_antitriplet!=NULL)	
    msg_Out()<<p_antitriplet->m_mom<<" "
	     <<sqrt(ATOOLS::Max(p_antitriplet->m_mom.Abs2(),0.))
	     <<" ("<<p_antitriplet->m_mom.Abs2()<<") ";
  else msg_Out()<<" XXX ";
  msg_Out()<<" ---"<<std::endl; 
}
