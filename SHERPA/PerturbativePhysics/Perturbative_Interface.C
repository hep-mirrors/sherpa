#include "SHERPA/PerturbativePhysics/Perturbative_Interface.H"

#include "SHERPA/PerturbativePhysics/Shower_Handler.H"
#include "SHERPA/PerturbativePhysics/Matrix_Element_Handler.H"
#include "SHERPA/Single_Events/Decay_Handler_Base.H"
#include "SHERPA/PerturbativePhysics/Hard_Decay_Handler.H"
#include "SHERPA/PerturbativePhysics/MI_Handler.H"
#include "SHERPA/SoftPhysics/Soft_Collision_Handler.H"
#include "PDF/Main/Shower_Base.H"
#include "ATOOLS/Phys/Cluster_Amplitude.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "PHASIC++/Process/Single_Process.H"
#include "PHASIC++/Process/ME_Generator_Base.H"
#include "PHASIC++/Process/MCatNLO_Process.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Math/Random.H"

using namespace SHERPA;
using namespace PHASIC;
using namespace ATOOLS;

Perturbative_Interface::Perturbative_Interface
(Matrix_Element_Handler *const meh,Hard_Decay_Handler*const dec,Shower_Handler *const psh):
  p_me(meh), p_dec(dec), p_mi(NULL), p_hd(NULL), p_sc(NULL), p_shower(psh),
  p_ampl(NULL), m_cmode(0)
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
  p_ampl(NULL), m_cmode(0) {}

Perturbative_Interface::Perturbative_Interface
(Decay_Handler_Base *const hdh,Shower_Handler *const psh):
  p_me(NULL), p_mi(NULL), p_hd(hdh), p_sc(NULL), p_shower(psh),
  p_ampl(NULL), m_cmode(0) {}

Perturbative_Interface::Perturbative_Interface
(Soft_Collision_Handler *const sch,Shower_Handler *const psh):
  p_me(NULL), p_mi(NULL), p_hd(NULL), p_sc(sch), p_shower(psh),
  p_ampl(NULL), m_cmode(0)  {}

Perturbative_Interface::~Perturbative_Interface() 
{
  if (p_ampl) p_ampl->Delete();
}

bool Perturbative_Interface::
SetColours(Cluster_Amplitude *ampl,Blob *const bl)
{
  if (!bl) return false;
  if (!ampl) {
    Particle * part;
    for (size_t i(0);i<bl->NInP();++i) {
      part = bl->InParticle(i); 
      for (size_t j(1);j<3;j++) {
	part->SetFlow(j,Max(0,int(part->GetFlow(j)-100)));
      }
    }
    for (size_t i(0);i<bl->NOutP();++i) {
      part = bl->OutParticle(i); 
      for (size_t j(1);j<3;j++) {
	part->SetFlow(j,Max(0,int(part->GetFlow(j)-100)));
      }
    }
    return true;
  }
  while (ampl->Prev()) ampl=ampl->Prev();
  for (size_t i(0);i<ampl->Legs().size();++i) {
    Cluster_Leg *cl(ampl->Leg(i));
    if (i<ampl->NIn()) {
      bl->InParticle(i)->SetFlow(1,Max(0,cl->Col().m_j-100));
      bl->InParticle(i)->SetFlow(2,Max(0,cl->Col().m_i-100));
    }
    else {
      bl->OutParticle(i-ampl->NIn())->SetFlow(1,Max(0,cl->Col().m_i-100));
      bl->OutParticle(i-ampl->NIn())->SetFlow(2,Max(0,cl->Col().m_j-100));
    }
  }
  return true;
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
    p_ampl->Delete();
    p_ampl=NULL;
  }
  p_shower->CleanUp();
  msg_Indent();
  if (p_mi) {
    p_ampl=p_mi->ClusterConfiguration();
    if (!SetColours(p_ampl,blob)) return Return_Value::New_Event;
    p_mi->Process()->Generator()->SetMassMode(1);
    int stat(p_mi->Process()->Generator()->ShiftMasses(p_ampl));
    if (stat<0) {
      msg_Tracking()<<METHOD<<"(): Mass shift failed. Reject event."<<std::endl;
      return Return_Value::New_Event;
    }
    if (stat==1) {
      stat=p_mi->Shower()->GetShower()->
	GetClusterDefinitions()->ReCluster(p_ampl);
      if (stat!=1) {
	msg_Tracking()<<METHOD<<"(): Reclustering failed. Reject event.\n";
	return Return_Value::New_Event;
      }
    }
    if (!p_shower->GetShower()->PrepareShower(p_ampl))
      return Return_Value::New_Event;
    return Return_Value::Success;
  }
  if (p_hd) {
    p_ampl=p_hd->ClusterConfiguration(blob);
    if (!SetColours(p_ampl,blob)) return Return_Value::New_Event;
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
    if (!SetColours(NULL,blob)) {
      msg_Out()<<METHOD<<": could not set colours.\n"; 
      return Return_Value::New_Event;
    }
    if (!p_shower->GetShower()->PrepareShower(p_ampl,true)) {
      msg_Out()<<METHOD<<": could not prepare shower.\n"; 
      return Return_Value::New_Event;
    }
    return Return_Value::Success;
  }
  p_ampl=p_me->Process()->Get<Single_Process>()->Cluster(m_cmode);
  if (p_ampl==NULL) return Return_Value::New_Event;
  p_me->Process()->Generator()->SetMassMode(1);
  int stat(p_me->Process()->Generator()->ShiftMasses(p_ampl));
  if (stat<0) {
    msg_Tracking()<<METHOD<<"(): Mass shift failed. Reject event."<<std::endl;
    return Return_Value::New_Event;
  }
  if (stat==1) {
    stat=p_me->Shower()->GetShower()->
      GetClusterDefinitions()->ReCluster(p_ampl);
    if (stat!=1) {
      msg_Tracking()<<METHOD<<"(): Reclustering failed. Reject event."<<std::endl;
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
    p_dec->DefineInitialConditions(p_ampl, blob);
  }
  while (p_ampl->Prev()) p_ampl=p_ampl->Prev();
  //if (!SetColours(p_ampl,blob)) return Return_Value::New_Event;
  //if (!p_dec->SetColours(p_ampl,blob)) return Return_Value::New_Event;
  m_weight=1.0;
  if (p_me->Process()->Info().m_ckkw&1) {
    if ((m_bbarmode&1) && p_me->HasNLO() &&
        p_me->Process()->Parent()->Info().m_fi.NLOType()==nlo_type::lo) {
        if (!LocalKFactor(p_ampl)) {
          DEBUG_INFO("didn't find PowProc in exclusively clustered amplitude");
	}
    }
    blob->AddData("Sud_Weight",new Blob_Data<double>(m_weight));
    if (p_me->EventGenerationMode()!=0) {
      if (m_weight>=ran->Get()) m_weight=Max(1.0,m_weight);
      else return Return_Value::New_Event;
    }
    Blob_Data_Base *winfo((*blob)["Weight"]);
    if (!winfo) THROW(fatal_error,"No weight information in signal blob");
    double meweight(winfo->Get<double>());
    blob->AddData("Weight",new Blob_Data<double>(meweight*m_weight));
  }
  if (!p_shower->GetShower()->PrepareShower(p_ampl)) 
    return Return_Value::New_Event;
  return Return_Value::Success;
}

bool Perturbative_Interface::LocalKFactor(ATOOLS::Cluster_Amplitude* ampl)
{
  if (m_globalkfac) {
    m_weight*=m_globalkfac;
    return true;
  }
  DEBUG_FUNC(ampl->Legs().size());
  Process_Vector procs(p_me->AllProcesses());
  Process_Base::SortFlavours(ampl);
  while (ampl->Next()!=NULL) {
    ampl=ampl->Next();
    if (ampl->Next() && (m_bbarmode&2)) continue;
    Process_Base::SortFlavours(ampl);
    for (size_t i=0; i<procs.size(); ++i) {
      MCatNLO_Process* mcnloproc=dynamic_cast<MCatNLO_Process*>(procs[i]);
      if (mcnloproc) {
	double K(mcnloproc->LocalKFactor(*ampl));
	if (K==0.0 || dabs(K)>m_maxkfac) continue;
	m_weight*=K;
	return true;
      }
    }
  }
  // no process found along ampl
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
      for (int i(0);i<p_hard->NOutP();++i)
        sblob->AddToInParticles(p_hard->OutParticle(i));
    }
  }
  blobs->push_back(sblob);
  p_shower->FillBlobs(blobs); 
  return true;
}

int Perturbative_Interface::PerformShowers()
{
  Blob_Data_Base *winfo((*p_hard)["Weight"]);
  if (!winfo) THROW(fatal_error,"No weight information in signal blob");
  double meweight(winfo->Get<double>());
  if (meweight==0.0) return 0;
  PDF::Shower_Base *csh(p_shower->GetShower());
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

