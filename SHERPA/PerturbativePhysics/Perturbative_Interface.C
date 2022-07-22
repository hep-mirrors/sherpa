#include "SHERPA/PerturbativePhysics/Perturbative_Interface.H"

#include "SHERPA/PerturbativePhysics/Shower_Handler.H"
#include "SHERPA/PerturbativePhysics/Matrix_Element_Handler.H"
#include "SHERPA/Single_Events/Decay_Handler_Base.H"
#include "SHERPA/PerturbativePhysics/Hard_Decay_Handler.H"
#include "SHERPA/PerturbativePhysics/MI_Handler.H"
#include "SHERPA/SoftPhysics/Soft_Collision_Handler.H"
#include "PDF/Main/Shower_Base.H"
#include "ATOOLS/Phys/Cluster_Amplitude.H"
#include "ATOOLS/Phys/Variations.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "PHASIC++/Process/Single_Process.H"
#include "PHASIC++/Process/ME_Generator_Base.H"
#include "PHASIC++/Process/MCatNLO_Process.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "ATOOLS/Math/Random.H"

#include <cassert>

using namespace SHERPA;
using namespace PHASIC;
using namespace ATOOLS;

std::ostream & SHERPA::operator<<(std::ostream & s,
                                  const SHERPA::mets_bbar_mode::code & mm)
{
  if (mm==mets_bbar_mode::none)         s<<"none";
  if (mm&mets_bbar_mode::enabled)       s<<"Enabled";
  if (mm&mets_bbar_mode::lowestmulti)   s<<"LowestMulti";
  if (mm&mets_bbar_mode::exclcluster)   s<<"ExclCluster";
  return s;
}

std::istream &SHERPA::operator>>(std::istream &s,
                                 SHERPA::mets_bbar_mode::code &mm)
{
  std::string tag;
  getline(s,tag);
  mm=mets_bbar_mode::none;
  if (tag.find("Enabled")!=std::string::npos)
    mm|=mets_bbar_mode::enabled;
  if (tag.find("LowestMulti")!=std::string::npos)
    mm|=mets_bbar_mode::lowestmulti;
  if (tag.find("ExclCluster")!=std::string::npos)
    mm|=mets_bbar_mode::exclcluster;
  return s;
}

Perturbative_Interface::Perturbative_Interface(Matrix_Element_Handler *const meh,
                                               Hard_Decay_Handler*const dec,
                                               Shower_Handler *const psh):
  p_me(meh), p_dec(dec), p_mi(NULL), p_hd(NULL), p_sc(NULL), p_shower(psh),
  p_ampl(NULL), m_globalkfac(0.), m_maxkfac(0.),
  m_bbarmode(mets_bbar_mode::none)
{
  Settings& s = Settings::GetMainSettings();
  m_bbarmode = s["METS_BBAR_MODE"].SetDefault(mets_bbar_mode::enabled|
                                              mets_bbar_mode::exclcluster)
                                  .Get<mets_bbar_mode::code>();
  m_globalkfac = s["GLOBAL_KFAC"].SetDefault(0.0).Get<double>();
  m_maxkfac = s["MENLOPS_MAX_KFAC"].SetDefault(10.0).Get<double>();
}

Perturbative_Interface::Perturbative_Interface
(MI_Handler *const mi,Shower_Handler *const psh):
  p_me(NULL), p_mi(mi), p_hd(NULL), p_sc(NULL), p_shower(psh),
  p_ampl(NULL) {}

Perturbative_Interface::Perturbative_Interface
(Decay_Handler_Base *const hdh,Shower_Handler *const psh):
  p_me(NULL), p_mi(NULL), p_hd(hdh), p_sc(NULL), p_shower(psh),
  p_ampl(NULL) {}

Perturbative_Interface::Perturbative_Interface
(Soft_Collision_Handler *const sch,Shower_Handler *const psh):
  p_me(NULL), p_mi(NULL), p_hd(NULL), p_sc(sch), p_shower(psh),
  p_ampl(NULL) {}

Perturbative_Interface::~Perturbative_Interface() 
{
  if (p_ampl) {
    Cluster_Amplitude *campl(p_ampl);
    while (campl->Prev()) campl=campl->Prev();
    campl->Delete();
  }
}

Return_Value::code
Perturbative_Interface::DefineInitialConditions(ATOOLS::Blob* blob,
                                                ATOOLS::Blob_List* bloblist)
{
  if (blob==NULL) {
    msg_Error()<<METHOD<<"(): Signal process not found."<<std::endl;
    return Return_Value::Error;
  }
  p_hard=blob;
  p_bloblist=bloblist;
  if (!p_shower->On() ||
      (p_me && p_me->Process()->Info().m_nlomode==nlo_mode::fixedorder)) {
    m_weightsmap.Clear();
    m_lkfweightsmap.Clear();
    return Return_Value::Success;
  }
  if (p_ampl) {
    Cluster_Amplitude *campl(p_ampl);
    while (campl->Prev()) campl=campl->Prev();
    campl->Delete();
    p_ampl=NULL;
  }
  p_shower->CleanUp();
  //msg_Out()<<METHOD<<"(MI = "<<p_mi<<", SC = "<<p_sc<<", "
  //	   <<bloblist->size()<<" blobs:\n"<<(*bloblist)<<"\n"
  //	   <<"Relevant blob = \n"<<(*blob)<<"\n";
  //msg_Out()<<"*** "<<METHOD<<"("<<bloblist->size()<<" blobs, id = "<<blob->Id()<<").\n";
  msg_Indent();
  if (p_mi) {
    p_ampl=p_mi->ClusterConfiguration(blob);
    if (p_ampl==NULL) return Return_Value::Retry_Event;
    if (p_ampl->Leg(0)->Mom()[3]*p_ampl->Leg(1)->Mom()[3]>0.0) {
      msg_Tracking()<<METHOD<<"(): Invalid beams. Retry event."<<std::endl;
      return Return_Value::Retry_Event;
    }
    p_mi->SetMassMode(1);
    int stat(p_mi->ShiftMasses(p_ampl));
    if (stat<0) {
      msg_Error()<<METHOD<<"(): MI Mass shift failed. Reject event.\n"
		 <<(*blob)<<"\n";
      exit(1);
      return Return_Value::Retry_Event;
    }
    if (stat==1) {
      stat=p_mi->Shower()->GetShower()->GetClusterDefinitions()->ReCluster(p_ampl);
      if (stat!=1) {
	msg_Error()<<METHOD<<"(): MI Reclustering failed. Reject event.\n"
		   <<(*blob)<<"\n";
	exit(1);
	return Return_Value::Retry_Event;
      }
    }
    if (!p_shower->GetShower()->PrepareShower(p_ampl))
      return Return_Value::New_Event;
    return Return_Value::Success;
  }
  if (p_hd) {
    p_ampl=p_hd->ClusterConfiguration(blob);
    if (!p_shower->GetShower()->PrepareShower(p_ampl))
      return Return_Value::New_Event;
    return Return_Value::Success;
  }
  if (p_sc) {
    //msg_Out()<<METHOD<<" in sc-mode.\n";
    p_sc->SetClusterDefinitions(p_shower->GetShower()->GetClusterDefinitions());
    p_ampl=p_sc->ClusterConfiguration(blob);
    if (p_ampl==NULL) {
      msg_Error()<<METHOD<<": Soft_Collision_Handler has no amplitude.\n";
      return Return_Value::New_Event;
    }
    //msg_Out()<<"  Amplitude has been constructed:\n"<<(*p_ampl)<<"\n";
    if (!p_shower->GetShower()->PrepareShower(p_ampl,true)) {
      msg_Error()<<METHOD<<": could not prepare shower.\n"; 
      return Return_Value::New_Event;
    }
    return Return_Value::Success;
  }
  assert(p_me != NULL);
  p_ampl=p_me->Process()->Get<Single_Process>()->Cluster
    (p_me->Process()->Integrator()->Momenta());
  if (p_ampl==NULL) {
    msg_Error()<<METHOD<<"(): Clustering failed. Reject event."<<std::endl;
    return Return_Value::New_Event;
  }
  if (p_ampl->Leg(0)->Mom()[3]*p_ampl->Leg(1)->Mom()[3]>0.0) {
    msg_Tracking()<<METHOD<<"(): Invalid beams. Reject event."<<std::endl;
    return Return_Value::New_Event;
  }
  m_weightsmap.Clear();
  m_lkfweightsmap.Clear();
  if (p_me->Process()->Info().m_ckkw&1) {
    if ((m_bbarmode&mets_bbar_mode::enabled) && p_me->HasNLO() &&
        p_me->Process()->Parent()->Info().m_fi.NLOType()==nlo_type::lo) {
      if (!LocalKFactor(p_ampl)) DEBUG_INFO("Process not found");
    }
  }
  p_me->Process()->Generator()->SetMassMode(1);
  size_t cmax(0);
  for (size_t i(0);i<p_ampl->Legs().size();++i)
    cmax=Max(cmax,(size_t)p_ampl->Leg(i)->Col().m_i);
  while (Flow::Counter()<cmax);
  p_me->Process()->Parent()->SetRBMap(p_ampl);
  if (p_dec) {
    p_dec->SetCluster(p_me->Shower()->GetShower()->GetClusterDefinitions());
    if (!p_dec->DefineInitialConditions(p_ampl, blob)) {
      msg_Tracking()<<METHOD<<"(): Decay clustering failed. Reject event."<<std::endl;
      return Return_Value::Retry_Event;
    }
  }
  while (p_ampl->Prev()) p_ampl=p_ampl->Prev();
  if (p_me->Process()->Info().m_ckkw&1) {
    auto wgtmap = (*p_hard)["WeightsMap"]->Get<Weights_Map>();
    wgtmap *= m_lkfweightsmap;
    p_hard->AddData("WeightsMap",new Blob_Data<Weights_Map>(wgtmap));
    if (p_me->EventGenerationMode()!=0) {
      const auto disc = ran->Get();
      const auto abswgt = std::abs(m_lkfweightsmap.Nominal());
      if (abswgt < disc) {
        return Return_Value::New_Event;
      }
      m_lkfweightsmap /= Min(1.0, abswgt);
      wgtmap /= Min(1.0, abswgt);
    }
    p_hard->AddData("WeightsMap",new Blob_Data<Weights_Map>(wgtmap));
  }
  if (!p_shower->GetShower()->PrepareShower(p_ampl)) 
    return Return_Value::New_Event;
  return Return_Value::Success;
}

bool Perturbative_Interface::LocalKFactor(ATOOLS::Cluster_Amplitude* ampl)
{
  if (m_globalkfac) {
    m_lkfweightsmap *= m_globalkfac;
    return true;
  }
  DEBUG_FUNC(ampl->Legs().size());
  Process_Vector procs(p_me->AllProcesses());
  Process_Base::SortFlavours(ampl);
  while (ampl->Next()!=NULL) {
    ampl=ampl->Next();
    Process_Base::SortFlavours(ampl);
    if (m_bbarmode&mets_bbar_mode::lowestmulti) {
      if (ampl->Next()) continue;
      Single_Process *proc(ampl->Proc<Single_Process>());
      if (ampl->Legs().size()-ampl->NIn()>
	  proc->Info().m_fi.NMinExternal()) break;
    }
    for (int i=procs.size()-1; i>=0; --i) {
      MCatNLO_Process* mcnloproc=dynamic_cast<MCatNLO_Process*>(procs[i]);
      if (mcnloproc) {
        Weights_Map K = mcnloproc->LocalKFactor(*ampl);
	if (K.Nominal()==0.0 || dabs(K.Nominal())>m_maxkfac) continue;
        m_lkfweightsmap *= K;
        return true;
      }
    }
  }
  // no process found along ampl
  return false;
}

bool Perturbative_Interface::FillBlobs()
{
  if (p_hard==NULL) return false;
  Blob *sblob = new Blob();
  sblob->SetType(btp::Shower);
  sblob->SetStatus(blob_status::needs_showers);
  sblob->SetId();
  sblob->SetPosition(p_hard->Position());
  if (p_shower->On()) {
    if (!p_hd) {
      for (int i(0);i<p_hard->NInP();++i)
	sblob->AddToOutParticles(p_hard->InParticle(i));
      for (size_t j(0);j<p_bloblist->size();++j) {
        Blob *cb((*p_bloblist)[j]);
        if (cb->Has(blob_status::needs_showers))
          for (int i(0);i<cb->NOutP();++i)
            if (cb->OutParticle(i)->DecayBlob()==NULL)
              sblob->AddToInParticles(cb->OutParticle(i));
      }
    }
    else {
      for (int i(0);i<p_hard->NOutP();++i) {
	if (!(p_hard->OutParticle(i)->GetFlow(1)==0 &&
	      p_hard->OutParticle(i)->GetFlow(2)==0))
	  sblob->AddToInParticles(p_hard->OutParticle(i));
      }
    }
  }
  p_bloblist->push_back(sblob);
  p_shower->FillBlobs(p_bloblist);
  //msg_Out()<<METHOD<<" has filled blob:\n"<<(*sblob)<<"\n";
  return true;
}

int Perturbative_Interface::PerformShowers()
{
  //msg_Out()<<"*** "<<METHOD<<".\n";
  PDF::Shower_Base *csh(p_shower->GetShower());
  int stat=csh->PerformShowers();
  if ((*p_hard)["WeightsMap"]!=NULL) {
    m_weightsmap["PS"]      = csh->WeightsMap()["PS"];
    m_weightsmap["PS_QCUT"] = csh->WeightsMap()["PS_QCUT"];
    auto wgtmap = (*p_hard)["WeightsMap"]->Get<Weights_Map>();
    wgtmap["PS"]      *= m_weightsmap["PS"];
    wgtmap["PS_QCUT"] *= m_weightsmap["PS_QCUT"];
    p_hard->AddData("WeightsMap",new Blob_Data<Weights_Map>(wgtmap));
  }
  return stat;
}

int Perturbative_Interface::PerformDecayShowers()
{ 
  return p_shower->GetShower()->PerformDecayShowers(); 
}

void Perturbative_Interface::CleanUp()
{
  if (p_me && p_me->Process())
    p_me->Process()->Generator()->SetMassMode(0);
  if (p_mi && p_mi->Process())
    p_mi->Process()->Generator()->SetMassMode(0);
  p_shower->CleanUp();
}

