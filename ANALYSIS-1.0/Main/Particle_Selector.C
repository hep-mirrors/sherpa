#include "Particle_Selector.H"
#include "Primitive_Analysis.H"

using namespace ANALYSIS;
using namespace ATOOLS;

Particle_Selector::Particle_Selector(const std::string & inlistname1,
				     const std::string & inlistname2,
				     const std::string & outlistname, int mode)  : 
  m_inlistname1(inlistname1),m_inlistname2(inlistname2),m_outlistname(outlistname),
  m_mode(mode)
{
  m_splitt_flag = false;
  m_name        = std::string("ParticleSelector_")+outlistname;
  
  if (m_mode==0) {
    if (outlistname=="ChargedHadron") m_mode=1;
    else if (outlistname=="NeutralHadron") m_mode=2;
    else if (outlistname=="Hadron") m_mode=3;
    else if (outlistname=="ChargedParticle") m_mode=4;
    else if (outlistname=="NeutralParticle") m_mode=5;
    else if (outlistname=="ChargedPion") m_mode=11;
    else if (outlistname=="ChargedKaon") m_mode=12;
    else if (outlistname=="ProtonAntiproton") m_mode=13;
    else if (outlistname=="Parton") m_mode=21;
    else if (outlistname=="All") m_mode=22;
    else if (outlistname=="NeutralPion") m_mode=101;
    else if (outlistname=="NeutralKaon") m_mode=102;
    else if (outlistname=="ChargedKStar") m_mode=103;
    else if (outlistname=="NeutralKStar") m_mode=104;
    else if (outlistname=="Eta") m_mode=105;
    else if (outlistname=="Rho0") m_mode=106;
    else if (outlistname=="Omega") m_mode=107;
    else if (outlistname=="EtaPrime") m_mode=108;
    else if (outlistname=="Phi") m_mode=109;
    else if (outlistname=="Lambda") m_mode=110;
    else if (outlistname=="ChargedSigma") m_mode=111;
    else if (outlistname=="ChargedXi") m_mode=112;
    else if (outlistname=="NeutralXi") m_mode=113;
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
  case 21: p_qualifier = new Is_Parton(); break;
  case 22: p_qualifier = new Is_There(); break;

  case 101: p_qualifier = new Is_Neutral_Pion(); break;
  case 102: p_qualifier = new Is_Neutral_Kaon(); break; 
  case 103: p_qualifier = new Is_Charged_KStar(); break; 
  case 104: p_qualifier = new Is_Neutral_KStar(); break; 
  case 105: p_qualifier = new Is_Eta(); break;
  case 106: p_qualifier = new Is_Rho0(); break; 
  case 107: p_qualifier = new Is_Omega(); break;
  case 108: p_qualifier = new Is_EtaPrime(); break; 
  case 109: p_qualifier = new Is_Phi(); break; 
  case 110: p_qualifier = new Is_Lambda(); break; 
  case 111: p_qualifier = new Is_Charged_Sigma(); break;
  case 112: p_qualifier = new Is_Charged_Xi(); break; 
  case 113: p_qualifier = new Is_Neutral_Xi(); break; 
  default:
    msg.Error()<<"ERROR in Particle_Selector: unknown particle qualifier "<<m_mode<<std::endl;
    p_qualifier = new Is_Charged(); break;
  }
}

void Particle_Selector::CreateParticleList()
{
  Particle_List * pl_in = NULL;
  if (m_mode<100) pl_in = p_ana->GetParticleList(m_inlistname1);
  else pl_in = p_ana->GetParticleList(m_inlistname2);
  if (pl_in==NULL) {
    msg.Out()<<"WARNING in Particle_Selector::Evaluate : particle list ";
    if (m_mode<100) msg.Out()<<m_inlistname1;
    else msg.Out()<<m_inlistname2;
    msg.Out()<<" not found "<<std::endl;
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
  return new Particle_Selector(m_inlistname1,m_inlistname2,m_outlistname,m_mode);
}

Particle_Selector::~Particle_Selector()
{
  if (p_qualifier) delete p_qualifier;
}
