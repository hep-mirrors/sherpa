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
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Math/Random.H"

#include <cassert>

using namespace SHERPA;
using namespace PHASIC;
using namespace ATOOLS;

Perturbative_Interface::Perturbative_Interface
(Matrix_Element_Handler *const meh,Hard_Decay_Handler*const dec,Shower_Handler *const psh):
  p_me(meh), p_dec(dec), p_mi(NULL), p_hd(NULL), p_sc(NULL), p_shower(psh),
  p_ampl(NULL), m_cmode(0), p_localkfactorvarweights(NULL)
{
  Data_Reader read(" ",";","!","=");
  read.AddComment("#");
  read.SetInputPath(p_me->Path());
  read.SetInputFile(p_me->File());
  m_cmode=ToType<int>(rpa->gen.Variable("METS_CLUSTER_MODE"));
  m_bbarmode=read.GetValue<int>("METS_BBAR_MODE",1);
  m_globalkfac=read.GetValue<double>("GLOBAL_KFAC",0.);
  m_maxkfac=read.GetValue<double>("MENLOPS_MAX_KFAC",10.0);
}

Perturbative_Interface::Perturbative_Interface
(MI_Handler *const mi,Shower_Handler *const psh):
  p_me(NULL), p_mi(mi), p_hd(NULL), p_sc(NULL), p_shower(psh),
  p_ampl(NULL), m_cmode(0), p_localkfactorvarweights(NULL) {}

Perturbative_Interface::Perturbative_Interface
(Decay_Handler_Base *const hdh,Shower_Handler *const psh):
  p_me(NULL), p_mi(NULL), p_hd(hdh), p_sc(NULL), p_shower(psh),
  p_ampl(NULL), m_cmode(0), p_localkfactorvarweights(NULL) {}

Perturbative_Interface::Perturbative_Interface
(Soft_Collision_Handler *const sch,Shower_Handler *const psh):
  p_me(NULL), p_mi(NULL), p_hd(NULL), p_sc(sch), p_shower(psh),
  p_ampl(NULL), m_cmode(0), p_localkfactorvarweights(NULL)  {}

Perturbative_Interface::~Perturbative_Interface() 
{
  if (p_ampl) {
    Cluster_Amplitude *campl(p_ampl);
    while (campl->Prev()) campl=campl->Prev();
    campl->Delete();
  }
  if (p_localkfactorvarweights) {
    delete p_localkfactorvarweights;
  }
}

Return_Value::code Perturbative_Interface::
DefineInitialConditions(ATOOLS::Blob *blob) 
{
  if (blob==NULL) {
    msg_Error()<<METHOD<<"(): Signal process not found."<<std::endl;
    return Return_Value::Error;
  }
  p_hard=blob;
  if (!p_shower->On() || (p_me && p_me->Process()->Info().m_nlomode==1)) {
    m_weight=1.0;
    return Return_Value::Success;
  }
  if (p_ampl) {
    Cluster_Amplitude *campl(p_ampl);
    while (campl->Prev()) campl=campl->Prev();
    campl->Delete();
    p_ampl=NULL;
  }
  p_shower->CleanUp();
  msg_Indent();
  if (p_mi) {
    p_ampl=p_mi->ClusterConfiguration();
    p_mi->Process()->Generator()->SetMassMode(1);
    int stat(p_mi->Process()->Generator()->ShiftMasses(p_ampl));
    if (stat<0) {
      msg_Tracking()<<METHOD<<"(): MI Mass shift failed. Reject event."<<std::endl;
      return Return_Value::Retry_Event;
    }
    if (stat==1) {
      stat=p_mi->Shower()->GetShower()->
	GetClusterDefinitions()->ReCluster(p_ampl);
      if (stat!=1) {
	msg_Tracking()<<METHOD<<"(): MI Reclustering failed. Reject event.\n";
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
    p_sc->SetClusterDefinitions(p_shower->GetShower()->GetClusterDefinitions());
    p_ampl=p_sc->ClusterConfiguration(blob);
    if (p_ampl==NULL) {
      msg_Out()<<METHOD<<": Soft_Collision_Handler has no amplitude.\n";
      return Return_Value::New_Event;
    }
    if (!p_shower->GetShower()->PrepareShower(p_ampl,true)) {
      msg_Out()<<METHOD<<": could not prepare shower.\n"; 
      return Return_Value::New_Event;
    }
    return Return_Value::Success;
  }
  assert(p_me != NULL);
  p_ampl=p_me->Process()->Get<Single_Process>()->Cluster
    (p_me->Process()->Integrator()->Momenta(),m_cmode);
  if (p_ampl==NULL) return Return_Value::New_Event;
  if (p_ampl->MS()==NULL)
    p_ampl=p_me->Process()->Get<Single_Process>()->Cluster
      (p_me->Process()->Integrator()->Momenta(),m_cmode|256);
  if (p_ampl==NULL) return Return_Value::New_Event;
  m_weight=1.0;
  if (p_localkfactorvarweights) {
    delete p_localkfactorvarweights;
    p_localkfactorvarweights = NULL;
  }
  if (p_me->Process()->Info().m_ckkw&1) {
    if ((m_bbarmode&1) && p_me->HasNLO() &&
        p_me->Process()->Parent()->Info().m_fi.NLOType()==nlo_type::lo) {
      Cluster_Amplitude *oampl=p_me->Process()->
	Get<Single_Process>()->Cluster
	(p_me->Process()->Integrator()->Momenta(),m_cmode);
      if (!LocalKFactor(oampl)) {
	DEBUG_INFO("didn't find process using original amplitude");
	if (m_bbarmode&4) {
	  Cluster_Amplitude *ampl=p_me->Process()->
	    Get<Single_Process>()->Cluster
	    (p_me->Process()->Integrator()->Momenta(),m_cmode|16|256|512);
	  while (ampl->Prev()) ampl=ampl->Prev();
	  if (!LocalKFactor(ampl))
	    DEBUG_INFO("didn't find process using exclusive clustering");
	  ampl->Delete();
	}
      }
      while (oampl->Prev()) oampl=oampl->Prev();
      oampl->Delete();
    }
  }
  p_me->Process()->Generator()->SetMassMode(1);
  int stat(p_me->Process()->Generator()->ShiftMasses(p_ampl));
  if (stat<0) {
    msg_Tracking()<<METHOD<<"(): ME Mass shift failed. Reject event."<<std::endl;
    return Return_Value::New_Event;
  }
  if (stat==1) {
    stat=p_me->Shower()->GetShower()->
      GetClusterDefinitions()->ReCluster(p_ampl);
    if (stat!=1) {
      msg_Tracking()<<METHOD<<"(): ME Reclustering failed. Reject event."<<std::endl;
      return Return_Value::New_Event;
    }
  }
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
    Cluster_Amplitude *ampl(p_ampl);
    while (ampl->Prev()) ampl=ampl->Prev();
    int stat(p_me->Process()->Generator()->ShiftMasses(ampl));
    if (stat<0) {
      msg_Tracking()<<METHOD<<"(): DH Mass shift failed. Reject event."<<std::endl;
      return Return_Value::Retry_Event;
    }
    if (stat==1) {
      stat=p_me->Shower()->GetShower()->
	GetClusterDefinitions()->ReCluster(ampl);
      if (stat!=1) {
	msg_Tracking()<<METHOD<<"(): DH Reclustering failed. Reject event."<<std::endl;
	return Return_Value::Retry_Event;
      }
    }
  }
  while (p_ampl->Prev()) p_ampl=p_ampl->Prev();
  if (p_me->Process()->Info().m_ckkw&1) {
    blob->AddData("Sud_Weight",new Blob_Data<double>(m_weight));
    if (p_me->EventGenerationMode()!=0) {
      if (m_weight>=ran->Get()) {
        if (m_weight < 1.0) {
          *p_localkfactorvarweights *= 1.0 / m_weight;
          m_weight = 1.0;
        }
      } else {
        return Return_Value::New_Event;
      }
    }
    Blob_Data_Base *winfo((*blob)["Weight"]);
    if (!winfo) THROW(fatal_error,"No weight information in signal blob");
    double meweight(winfo->Get<double>());
    blob->AddData("Weight",new Blob_Data<double>(meweight*m_weight));
    // also update reweighting weights
    Blob_Data_Base *vws((*blob)["Variation_Weights"]);
    if (vws) {
      if (p_localkfactorvarweights) {
        vws->Get<Variation_Weights>() *= *p_localkfactorvarweights;
      } else {
        vws->Get<Variation_Weights>() *= m_weight;
      }
    }
  }
  if (!p_shower->GetShower()->PrepareShower(p_ampl)) 
    return Return_Value::New_Event;
  return Return_Value::Success;
}

bool Perturbative_Interface::LocalKFactor(ATOOLS::Cluster_Amplitude* ampl)
{
  if (p_localkfactorvarweights) {
    delete p_localkfactorvarweights;
    p_localkfactorvarweights = NULL;
  }
  if (m_globalkfac) {
    m_weight*=m_globalkfac;
    return true;
  }
  DEBUG_FUNC(ampl->Legs().size());
  Process_Vector procs(p_me->AllProcesses());
  Process_Base::SortFlavours(ampl);
  if (p_hard) {
    Blob_Data_Base *vws((*p_hard)["Variation_Weights"]);
    if (vws) {
      if (p_localkfactorvarweights) {
        delete p_localkfactorvarweights;
      }
      Variations *variations = vws->Get<Variation_Weights>().GetVariations();
      p_localkfactorvarweights = new Variation_Weights(variations);
    }
  }
  while (ampl->Next()!=NULL) {
    ampl=ampl->Next();
    if (ampl->Next() && (m_bbarmode&2)) continue;
    Process_Base::SortFlavours(ampl);
    for (size_t i=0; i<procs.size(); ++i) {
      if (p_localkfactorvarweights) p_localkfactorvarweights->Reset();
      MCatNLO_Process* mcnloproc=dynamic_cast<MCatNLO_Process*>(procs[i]);
      if (mcnloproc) {
        if (mcnloproc->VariationWeights()) THROW(fatal_error, "Variation weights already set.");
        mcnloproc->SetVariationWeights(p_localkfactorvarweights);
	double K(mcnloproc->LocalKFactor(*ampl));
        mcnloproc->SetVariationWeights(NULL);
	if (K==0.0 || dabs(K)>m_maxkfac) continue;
	m_weight*=K;
	return true;
      }
    }
  }
  // no process found along ampl
  if (p_localkfactorvarweights) {
    delete p_localkfactorvarweights;
    p_localkfactorvarweights = NULL;
  }
  return false;
}

bool Perturbative_Interface::FillBlobs(ATOOLS::Blob_List *blobs)
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
      for (size_t j(0);j<blobs->size();++j) {
        Blob *cb((*blobs)[j]);
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
  blobs->push_back(sblob);
  p_shower->FillBlobs(blobs); 
  return true;
}

int Perturbative_Interface::PerformShowers()
{
  // see if the event has any weight
  Blob_Data_Base *winfo((*p_hard)["Weight"]);
  if (!winfo) THROW(fatal_error,"No weight information in signal blob");
  double meweight(winfo->Get<double>());
  if (meweight==0.0) return 0;

  PDF::Shower_Base *csh(p_shower->GetShower());

  // look for reweightings and set up the shower accordingly
  Blob_Data_Base *blob_data_base((*p_hard)["Variation_Weights"]);
  if (blob_data_base) {
    csh->SetVariationWeights(&blob_data_base->Get<Variation_Weights>());
  }

  int stat=csh->PerformShowers();
  double weight=csh->Weight();
  p_hard->AddData("Shower_Weight",new Blob_Data<double>(weight));
  p_hard->AddData("Weight",new Blob_Data<double>(meweight*weight));
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

