#include "Particle_Selector.H"
#include "Primitive_Analysis.H"

using namespace ANALYSIS;
using namespace ATOOLS;

Particle_Selector::Particle_Selector(const std::string & inlistname,const std::string & outlistname, int mode)  : 
  m_inlistname(inlistname),m_outlistname(outlistname),
  m_mode(mode)
{
  m_splitt_flag = false;
  m_name        = std::string("ParticleSelector_")+outlistname;
  
  if (m_mode==0) {
    if (outlistname=="ChargedHadron") m_mode=1;
    else if (outlistname=="NeutralHadron") m_mode=2;
    else if (outlistname=="Hadron") m_mode=3;
    else if (outlistname=="ChargedParticle") m_mode=4;
    else if (outlistname=="ChargedPion") m_mode=11;
    else if (outlistname=="ChargedKaon") m_mode=12;
    else if (outlistname=="ProtonAntiproton") m_mode=13;
    else {
      msg.Error()<<"ERROR in Particle_Selector: unknown particle qualifier "<<outlistname<<std::endl;
      m_mode=4;
    }
  }

  switch (m_mode) {
  case 1:  p_qualifier = new Is_Charged_Hadron(); break;
  case 2:  p_qualifier = new Is_Neutral_Hadron(); break;
  case 3:  p_qualifier = new Is_Hadron(); break;
  case 4:  p_qualifier = new Is_Charged(); break;
  case 5:  p_qualifier = new Is_Neutral(); break;
  case 11: p_qualifier = new Is_Charged_Pion(); break;
  case 12: p_qualifier = new Is_Charged_Kaon(); break;
  case 13: p_qualifier = new Is_Proton_Antiproton(); break;
  default:
    msg.Error()<<"ERROR in Particle_Selector: unknown particle qualifier "<<m_mode<<std::endl;
    p_qualifier = new Is_Charged(); break;
  }
}

void Particle_Selector::CreateParticleList()
{
  Particle_List * pl_in = p_ana->GetParticleList(m_inlistname);
  if (pl_in==NULL) {
    msg.Out()<<"WARNING in Particle_Selector::Evaluate : particle list "<<m_inlistname<<" not found "<<std::endl;
    return;
  }
  
  Particle_List * pl = new Particle_List;
  copy_if(pl_in->begin(),pl_in->end(),
	  back_inserter(*pl),*p_qualifier);
  
  p_ana->AddParticleList(m_outlistname,pl);
}

void Particle_Selector::Evaluate(const ATOOLS::Blob_List & ,double weight, int ncout)
{
  CreateParticleList();
}

Primitive_Observable_Base * Particle_Selector::Copy() const
{
  return new Particle_Selector(m_inlistname,m_outlistname,m_mode);
}

Particle_Selector::~Particle_Selector()
{
  if (p_qualifier) delete p_qualifier;
}
