#include "MI_Base.H"
#include "Message.H"
#include "Particle.H"

using namespace AMISIC;

MI_Base::NameMIBaseMap MI_Base::m_bases=MI_Base::NameMIBaseMap();

long int MI_Base::m_particlecounter=0;

bool MI_Base::m_stophard=true;
bool MI_Base::m_stopsoft=true;

MI_Base::MI_Base(std::string _m_name,TypeID _m_type,unsigned int _n_parameter):
  m_name(_m_name),
  m_type(_m_type),
  n_parameter(_n_parameter),
  p_blob(NULL),
  p_xs(NULL)
{
  for (NameMIBaseMapIterator nbit=m_bases.begin();nbit!=m_bases.end();++nbit) {
    if (nbit->first==m_name) {
      ATOOLS::msg.Error()<<"MI_Base::MI_Base("<<m_name<<","<<m_type<<"): "
			 <<"MI_Base already exists!"<<std::endl
			 <<"   Run cannot continue."<<std::endl;
      exit(120);
    }
  }
  if (m_type==Unknown) {
    ATOOLS::msg.Error()<<"MI_Base::MI_Base("<<m_name<<","<<m_type<<"): "
		       <<"Base has no type!"<<std::endl
		       <<"   Run cannot continue."<<std::endl;
    exit(121);
  }
  m_start = new double[n_parameter];
  m_stop = new double[n_parameter];
  m_last = new double[n_parameter];
  p_blob = new ATOOLS::Blob();
  m_bases[m_name]=this;
}

MI_Base::~MI_Base()
{
  for (NameMIBaseMapIterator nbit=m_bases.begin();nbit!=m_bases.end();++nbit) {
    if (nbit->first==m_name) {
      m_bases.erase(nbit--);
      break;
    }
  }
  delete p_blob;
  delete [] m_start;
  delete [] m_stop;
  delete [] m_last;
}

void MI_Base::UpdateAll(const MI_Base *mibase)
{
  for (NameMIBaseMapIterator nbit=m_bases.begin();nbit!=m_bases.end();++nbit) {
    nbit->second->Update(mibase);
  }  
}

void MI_Base::Update(const MI_Base *mibase)
{
  ATOOLS::msg.Error()<<"MI_Base::Update("<<mibase<<"): "
		     <<"Virtual method called!"<<std::endl;
  return;
}

bool MI_Base::Initialize()
{
  ATOOLS::msg.Error()<<"MI_Base::Initialize(): "
		     <<"Virtual method called!"<<std::endl;
  return false;
}

void MI_Base::Reset()
{
  ATOOLS::msg.Error()<<"MI_Base::Reset(): "
		     <<"Virtual method called!"<<std::endl;
  return;
}

bool MI_Base::DiceOrderingParameter()
{
  ATOOLS::msg.Error()<<"MI_Base::DiceOrderingParameter(): "
		     <<"Virtual method called!"<<std::endl;
  return false;
}

bool MI_Base::DiceProcess()
{
  ATOOLS::msg.Error()<<"MI_Base::DiceProcess(): "
		     <<"Virtual method called!"<<std::endl;
  return false;
}

void MI_Base::ResetAll()
{
  m_particlecounter=0;
  for (NameMIBaseMapIterator nbit=m_bases.begin();nbit!=m_bases.end();++nbit) {
    nbit->second->Reset();
  }  
}

bool MI_Base::CreateBlob(ATOOLS::Blob *blob)
{
  if (blob==NULL) {
    ATOOLS::msg.Error()<<"MI_Base::CreateBlob(..): "
		       <<"Blob is not initialized!"<<std::endl
		       <<"   Cannot proceed in filling."<<std::endl;
    return false;
  }
  if (p_blob==NULL) {
    ATOOLS::msg.Error()<<"MI_Base::CreateBlob(..): "
		       <<"Did not select any blob yet!"<<std::endl
		       <<"   Cannot proceed in filling."<<std::endl;
    return false;
  }
  bool _m_dicedprocess=m_dicedprocess;
  m_dicedprocess=false;
  ATOOLS::Particle *particle;
  for (unsigned int i=0;i<(unsigned int)p_blob->NInP();++i) {
    particle = new ATOOLS::Particle(m_particlecounter++,
				    p_blob->InParticle(i)->Flav(),
				    p_blob->InParticle(i)->Momentum());
    particle->SetFlow(1,p_blob->InParticle(i)->GetFlow(1));
    particle->SetFlow(2,p_blob->InParticle(i)->GetFlow(2));
    particle->SetStatus(1);
    blob->AddToInParticles(particle);
  }
  for (unsigned int i=0;i<(unsigned int)p_blob->NOutP();++i) {
    particle = new ATOOLS::Particle(m_particlecounter++,
				    p_blob->OutParticle(i)->Flav(),
				    p_blob->OutParticle(i)->Momentum());
    particle->SetFlow(1,p_blob->OutParticle(i)->GetFlow(1));
    particle->SetFlow(2,p_blob->OutParticle(i)->GetFlow(2));
    particle->SetStatus(1);
    blob->AddToOutParticles(particle);
  }
  p_blob->DeleteOwnedParticles();
  return _m_dicedprocess;
}

std::string MI_Base::TypeToString(TypeID type)
{
  switch (type) {
  case HardEvent: return std::string("Hard Event");
  case SoftEvent: return std::string("Soft Event");
  case Unknown  : return std::string("Unknown");
  }
  return std::string("Unknown");
}

MI_Base::TypeID MI_Base::StringToType(std::string type)
{
  if (type==std::string("Hard Event")) return HardEvent;
  if (type==std::string("Soft Event")) return HardEvent;
  return Unknown;
}

MI_None::MI_None(TypeID _m_type):
  MI_Base(TypeToString(_m_type),_m_type) {}

MI_None::~MI_None() {}

void MI_None::Update(const MI_Base *mibase) 
{
  return;
}

bool MI_None::Initialize()
{
  return true;
}

void MI_None::Reset()
{
  return;
}

bool MI_None::DiceOrderingParameter()
{
  return true;
}

bool MI_None::DiceProcess()
{
  m_dicedprocess=false;
  return true;
}

