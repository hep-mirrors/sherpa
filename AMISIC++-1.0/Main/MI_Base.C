#include "MI_Base.H"
#include "Message.H"
#include "Particle.H"

#ifdef PROFILE__MI_Base
#include "prof.hh"
#else 
#define PROFILE_HERE
#endif
#ifdef USING__Sherpa
#include "Matrix_Element_Handler.H"
#endif

using namespace AMISIC;

MI_Base::NameMIBaseMap MI_Base::s_bases=MI_Base::NameMIBaseMap();
std::vector<ATOOLS::Remnant_Info*> 
MI_Base::s_remnanthandlers=std::vector<ATOOLS::Remnant_Info*>(2,(ATOOLS::Remnant_Info*)NULL);

bool MI_Base::m_stophard=true;
bool MI_Base::m_stopsoft=true;

MI_Base::MI_Base(std::string _m_name,TypeID _m_type,unsigned int _m_nparameter,
		 unsigned int infiles,unsigned int outfiles):
  File_IO_Base(infiles,outfiles),
  m_name(_m_name),
  m_type(_m_type),
  m_nparameter(_m_nparameter),
  p_blob(NULL),
#ifdef USING__Sherpa
  p_mehandler(NULL),
#endif
  p_xs(NULL)
{
  for (NameMIBaseMapIterator nbit=s_bases.begin();nbit!=s_bases.end();++nbit) {
    if (nbit->first==m_name) {
      throw(ATOOLS::Exception(ATOOLS::ex::fatal_error,"MI_Base already exists!",
			      "MI_Base","MI_Base"));
    }
  }
  if (m_type==Unknown) {
    throw(ATOOLS::Exception(ATOOLS::ex::fatal_error,"MI base has no type!",
			    "MI_Base","MI_Base"));
  }
  m_start = new double[m_nparameter];
  m_stop = new double[m_nparameter];
  m_last = new double[m_nparameter];
  p_blob = new ATOOLS::Blob();
  s_bases[m_name]=this;
}

MI_Base::~MI_Base()
{
  for (NameMIBaseMapIterator nbit=s_bases.begin();nbit!=s_bases.end();++nbit) {
    if (nbit->first==m_name) {
      s_bases.erase(nbit--);
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
  PROFILE_HERE;
  for (NameMIBaseMapIterator nbit=s_bases.begin();nbit!=s_bases.end();++nbit) {
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
  PROFILE_HERE;
  for (NameMIBaseMapIterator nbit=s_bases.begin();nbit!=s_bases.end();++nbit) {
    nbit->second->Reset();
  }  
}

bool MI_Base::CreateBlob(ATOOLS::Blob *blob)
{
  PROFILE_HERE;
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
    particle = new ATOOLS::Particle(-1,p_blob->InParticle(i)->Flav(),
				    p_blob->InParticle(i)->Momentum());
    particle->SetFlow(1,p_blob->InParticle(i)->GetFlow(1));
    particle->SetFlow(2,p_blob->InParticle(i)->GetFlow(2));
    particle->SetNumber(1);
    particle->SetStatus(1);
    particle->SetInfo('G');
    blob->AddToInParticles(particle);
  }
  for (unsigned int i=0;i<(unsigned int)p_blob->NOutP();++i) {
    particle = new ATOOLS::Particle(-1,p_blob->OutParticle(i)->Flav(),
				    p_blob->OutParticle(i)->Momentum());
    particle->SetFlow(1,p_blob->OutParticle(i)->GetFlow(1));
    particle->SetFlow(2,p_blob->OutParticle(i)->GetFlow(2));
    particle->SetNumber(1);
    particle->SetStatus(1);
    particle->SetInfo('H');
    blob->AddToOutParticles(particle);
  }
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

MI_None::~MI_None() 
{
}

void MI_None::Update(const MI_Base *mibase) 
{
  return;
}

bool MI_None::Initialize()
{
#ifdef USING__Sherpa
  p_mehandler = new SHERPA::Matrix_Element_Handler();
#endif
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

