#include "Particle_Selector.H"
#include "Primitive_Analysis.H"
#include "Particle_Qualifier.H"
#include "MyStrStream.H"

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
  
  p_qualifier = 
    ATOOLS::Particle_Qualifier_Getter::GetObject(outlistname,"");
  if (!p_qualifier) {
    msg.Error()<<"ERROR in Particle_Selector: unknown particle qualifier "<<m_mode<<std::endl;
    p_qualifier = new Is_Charged();
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
