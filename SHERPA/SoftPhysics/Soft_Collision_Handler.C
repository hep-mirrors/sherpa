#include "SHERPA/SoftPhysics/Soft_Collision_Handler.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "SHRiMPS/Main/Shrimps.H"
#include "AMISIC++/Main/Amisic.H"

#ifdef PROFILE__all
#define PROFILE__Soft_Collision_Handler
#endif
#ifdef PROFILE__Soft_Collision_Handler
#include "prof.hh" 
#else
#define PROFILE_HERE
#endif

using namespace SHERPA;
using namespace ATOOLS;
using namespace std;

Soft_Collision_Handler::
Soft_Collision_Handler(AMISIC::Amisic * amisic,SHRIMPS::Shrimps * shrimps,
		       const bool bunch_rescatter) :
  m_bunch_rescatter(bunch_rescatter),
  m_mode(scmode::none),
  p_shrimps(NULL), p_amisic(NULL)
{
  Settings& s = Settings::GetMainSettings();
  m_dir     = s.GetPath();
  m_scmodel = (m_bunch_rescatter ?
	       s["BEAM_RESCATTERING"].SetDefault("None").UseNoneReplacements().Get<string>() :
	       s["SOFT_COLLISIONS"].SetDefault("None").UseNoneReplacements().Get<string>() );
  msg_Out()<<METHOD<<"("<<this<<"), model = "<<m_scmodel<<"\n";
  if (m_bunch_rescatter) {
    msg_Out()<<METHOD<<" for bunch rescattering: amisic = "<<amisic<<", shrimps = "<<shrimps<<".\n";
  }
  if (m_scmodel==string("Shrimps")) {
    m_mode    = scmode::shrimps;
    p_shrimps = shrimps;
    exh->AddTerminatorObject(this);
    return;
  }
  else if (m_scmodel==string("Amisic")) {
    m_mode    = scmode::amisic;
    p_amisic  = amisic;
    exh->AddTerminatorObject(this);
    return;
  }
  else if (m_scmodel==string("None")) return;
  THROW(critical_error,"Soft_Collision model not implemented.");
}
   
Soft_Collision_Handler::~Soft_Collision_Handler() 
{
  exh->RemoveTerminatorObject(this);
}

void Soft_Collision_Handler::CleanUp() {
  switch (m_mode) {
  case scmode::shrimps: 
    p_shrimps->CleanUp();
    break;
  case scmode::amisic:
    p_amisic->CleanUpMinBias();
    break;
  case scmode::none:
  default:
    break;
  }
} 

void Soft_Collision_Handler::PrepareTerminate() {}

ATOOLS::Return_Value::code
Soft_Collision_Handler::GenerateMinimumBiasEvent(ATOOLS::Blob_List* blobs)
{
  PROFILE_HERE;
  int outcome(-1);
  switch (m_mode) {
  case scmode::shrimps: 
    outcome = p_shrimps->InitMinBiasEvent(blobs);
    break;
  case scmode::amisic: 
    outcome = p_amisic->InitMinBiasEvent();
    break;
  case scmode::none:
    outcome = 0;
    break;
  default:
    break;
  }
  switch (outcome) {
  case 1:  return Return_Value::Success;
  case 0:  return Return_Value::Nothing;
  default: break;
  }
  msg_Tracking()<<"Error in "<<METHOD<<":\n"
		<<"   Did not manage to produce a Minimum Bias event with "<<m_scmodel<<".\n";
  return Return_Value::New_Event;
}

ATOOLS::Return_Value::code
Soft_Collision_Handler::GenerateBunchRescatter(ATOOLS::Blob_List * blobs) {
  int outcome(-1);
  switch (m_mode) {
  case scmode::shrimps: 
    msg_Error()<<METHOD<<" not yet available for SHRiMPS.  Will exit the run.\n";
    exit(1);
  case scmode::amisic:
    outcome = p_amisic->InitRescatterEvent();
    break;
  case scmode::none:
    outcome = 0;
    break;
  default:
    break;
  }
  switch (outcome) {
  case 1:  return Return_Value::Success;
  case 0:  return Return_Value::Nothing;
  default: break;
  }  
  return Return_Value::Nothing;
}

void Soft_Collision_Handler::SetPosition(const size_t & beam,const Vec4D & pos) {
  switch (m_mode) {
  case scmode::shrimps: 
    msg_Error()<<METHOD<<" not yet available for SHRiMPS.  Will exit the run.\n";
    exit(1);
  case scmode::amisic:
    p_amisic->SetPosition(beam,pos);
    break;
  case scmode::none:
  default:
    break;
  }
}

Cluster_Amplitude *Soft_Collision_Handler::ClusterConfiguration(Blob *const blob)
{
  return p_shrimps->ClusterConfiguration(blob);
}


