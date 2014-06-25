#include "AMISIC++/Main/MI_Base.H"

#include "ATOOLS/Phys/Particle.H"
#include "ATOOLS/Org/Message.H"

using namespace AMISIC;

MI_Base::String_MI_Base_Map *MI_Base::s_bases;

bool MI_Base::s_stophard=true;
bool MI_Base::s_stopsoft=true;
bool MI_Base::s_cleaned=true;

MI_Base *MI_Base::s_hard=NULL;
MI_Base *MI_Base::s_soft=NULL;

MI_Base::MI_Base(std::string name,TypeID type,unsigned int nparameter,
		 unsigned int infiles,unsigned int outfiles):
  File_IO_Base(infiles,outfiles),
  m_name(name),
  m_type(type),
  m_start(NULL), m_stop(NULL), m_last(NULL),
  m_nparameter(nparameter),
  p_xs(NULL)
{
  static bool initialized=false;
  if (!initialized) {
    s_bases = new String_MI_Base_Map();
    initialized=true;
  }
  for (String_MI_Base_Map::iterator nbit=s_bases->begin();
       nbit!=s_bases->end();++nbit) {
    if (nbit->first==m_name) {
      THROW(fatal_error,"MI_Base already exists!");
    }
  }
  if (m_type==Unknown) {
    THROW(fatal_error,"MI base has no type!");
  }
  if (m_nparameter>0) {
    m_start = new double[m_nparameter];
    m_stop = new double[m_nparameter];
    m_last = new double[m_nparameter];
  }
  (*s_bases)[m_name]=this;
  switch (m_type) {
  case SoftEvent: 
    if (m_name!=TypeToString(type)+" None") s_soft=this; 
    break;
  case HardEvent: 
    if (m_name!=TypeToString(type)+" None") s_hard=this; 
    break;
  default: 
    break;
  }
}

MI_Base::~MI_Base()
{
  for (String_MI_Base_Map::iterator nbit=s_bases->begin();
       nbit!=s_bases->end();++nbit) {
    if (nbit->first==m_name) {
      s_bases->erase(nbit);
      break;
    }
  }
  if (m_nparameter>0) {
    delete [] m_start;
    delete [] m_stop;
    delete [] m_last;
  }
}

void MI_Base::UpdateAll(const MI_Base *mibase)
{
  for (String_MI_Base_Map::iterator nbit=s_bases->begin();
       nbit!=s_bases->end();++nbit) {
    nbit->second->Update(mibase);
  }  
}

void MI_Base::Update(const MI_Base *mibase)
{
  msg_Error()<<"MI_Base::Update("<<mibase<<"): "
		     <<"Virtual method called!"<<std::endl;
  return;
}

bool MI_Base::Initialize()
{
  msg_Error()<<"MI_Base::Initialize(): "
		     <<"Virtual method called!"<<std::endl;
  return false;
}

void MI_Base::Reset()
{
  msg_Error()<<"MI_Base::Reset(): "
		     <<"Virtual method called!"<<std::endl;
  return;
}

void MI_Base::CleanUp()
{
  for (String_MI_Base_Map::iterator nbit=s_bases->begin();
       nbit!=s_bases->end();++nbit) {
    nbit->second->m_inparticles.Clear();
    nbit->second->m_outparticles.Clear();
  }
  s_stophard=false;
  s_stopsoft=false;
  s_cleaned=true;
}

bool MI_Base::VetoProcess(ATOOLS::Blob *blob)
{
  msg_Error()<<"MI_Base::VetoProcess(): "
		     <<"Virtual method called!"<<std::endl;
  return false;
}

bool MI_Base::GenerateProcess()
{
  msg_Error()<<"MI_Base::GenerateProcess(): "
		     <<"Virtual method called!"<<std::endl;
  return false;
}

void MI_Base::ResetAll()
{
  for (String_MI_Base_Map::iterator nbit=s_bases->begin();
       nbit!=s_bases->end();++nbit) {
    nbit->second->Reset();
  }  
}

bool MI_Base::FillBlob(ATOOLS::Blob *blob)
{
  if (blob==NULL) {
    msg_Error()<<"MI_Base::FillBlob(..): "
		       <<"Blob is not initialized!"<<std::endl
		       <<"   Cannot proceed in filling."<<std::endl;
    return false;
  }
  if (!m_generatedprocess) return false;
  if (m_inparticles.empty()) {
    msg_Error()<<"MI_Base::FillBlob(..): "
		       <<"Did not create any process yet!"<<std::endl
		       <<"   Cannot proceed in filling."<<std::endl;
    return false;
  }
  bool generatedprocess=m_generatedprocess;
  m_generatedprocess=false;
  if (m_type==HardEvent) {
    blob->SetType(ATOOLS::btp::Hard_Collision);
    blob->SetStatus(ATOOLS::blob_status::needs_showers &
		    ATOOLS::blob_status::needs_beams &
		    ATOOLS::blob_status::needs_hadronization);
  }
  else {
    blob->SetType(ATOOLS::btp::Soft_Collision);
    blob->SetStatus(ATOOLS::blob_status::needs_beams &
		    ATOOLS::blob_status::needs_hadronization);
  }
  ATOOLS::Particle *particle;
  for (size_t i=0;i<m_inparticles.size();++i) {
    particle = new ATOOLS::Particle(-1,m_inparticles[i]->Flav(),
				    m_inparticles[i]->Momentum());
    particle->SetFlow(1,m_inparticles[i]->GetFlow(1));
    particle->SetFlow(2,m_inparticles[i]->GetFlow(2));
    particle->SetNumber(1);
    particle->SetStatus(ATOOLS::part_status::active);
    particle->SetInfo('G');
    blob->AddToInParticles(particle);
  }
  for (size_t i=0;i<m_outparticles.size();++i) {
    particle = new ATOOLS::Particle(-1,m_outparticles[i]->Flav(),
				    m_outparticles[i]->Momentum());
    particle->SetFlow(1,m_outparticles[i]->GetFlow(1));
    particle->SetFlow(2,m_outparticles[i]->GetFlow(2));
    particle->SetNumber(1);
    particle->SetStatus(ATOOLS::part_status::active);
    particle->SetInfo('H');
    blob->AddToOutParticles(particle);
  }
  return generatedprocess;
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

bool MI_Base::StopGeneration(TypeID type)
{ 
  switch (type) {
  case HardEvent: return s_stophard;
  case SoftEvent: return s_stopsoft;
  case Unknown: return s_stophard&&s_stopsoft;
  }
  return true;
}

void MI_Base::SetStopGeneration(TypeID type,const bool stop)
{ 
  switch (type) {
  case HardEvent: 
    s_stophard=stop; 
    break;
  case SoftEvent: 
    s_stopsoft=stop; 
    break;
  case Unknown: 
    s_stophard=s_stopsoft=stop;
  }
}

MI_None::MI_None(TypeID type):
  MI_Base(TypeToString(type)+" None",type) {}

MI_None::~MI_None() 
{
}

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
  switch (Type()) {
  case HardEvent: 
    s_stophard=true; 
    break;
  case SoftEvent: 
    s_stopsoft=true; 
    break;
  case Unknown: 
    THROW(fatal_error,"No type");
  }
  return;
}

bool MI_None::VetoProcess(ATOOLS::Blob *blob)
{
  return false;
}

bool MI_None::GenerateProcess()
{
  m_generatedprocess=false;
  return true;
}

