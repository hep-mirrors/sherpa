#include "Multiple_Interactions.H"

#include "Amisic.H"
#include "My_Limits.H"
#include "Remnant_Base.H"
#include "ISR_Handler.H"
#include "Jet_Finder.H"
#include "Run_Parameter.H"

#ifdef PROFILE__all
#define PROFILE__Multiple_Interactions
#endif
#ifdef PROFILE__Multiple_Interactions
#include "prof.hh"
#else 
#define PROFILE_HERE 
#endif

#ifdef ROOT_SUPPORT
#define ANALYSE__Multiple_Interactions
#include "My_Root.H"
#include "TH2D.h"
#endif

using namespace SHERPA;
using namespace AMISIC;
using namespace EXTRAXS;
using namespace ATOOLS;

Multiple_Interactions::Multiple_Interactions(MI_Handler *mihandler):
  p_mihandler(mihandler)
{
  m_name = std::string("Multiple_Interactions:")+p_mihandler->Name();
  m_type = eph::Perturbative;
  m_ecms = sqrt(p_mihandler->ISRHandler()->Pole());
  p_remnants[0]=GET_OBJECT(Remnant_Base,"Remnant_Base_0");
  p_remnants[1]=GET_OBJECT(Remnant_Base,"Remnant_Base_1");
  if (p_remnants[0]==NULL || p_remnants[1]==NULL) {
    THROW(fatal_error,"No beam remnant handler found.");
  }
}

Multiple_Interactions::~Multiple_Interactions() 
{
}

bool Multiple_Interactions::CheckBlobList(ATOOLS::Blob_List *const bloblist) 
{
  p_bloblist=bloblist;
  if (!p_bloblist->FourMomentumConservation()) {
    msg.Error()<<"Multiple_Interactions::CheckBlobList(..): "
	       <<"Retry event "<<rpa.gen.NumberOfDicedEvents()<<std::endl;
    p_bloblist->Clear();
  }
  for (size_t i=0;i<2;++i) {
    m_emax[i]=p_remnants[i]->BeamEnergy();
    p_remnants[i]->QuickClear();
  }
  for (Blob_List::const_iterator bit=bloblist->begin();
       bit!=bloblist->end();++bit) {
    if ((*bit)->Type()==btp::Hard_Collision ||
	(*bit)->Type()==btp::Signal_Process) 
      if ((*bit)->Status()!=0) return false;
  }
  Blob *signal=bloblist->FindFirst(btp::Signal_Process);
//     if (!m_diced && VetoHardProcess(signal)) {
//       m_ptmax=0.0;
//       return true;
//     }
  switch (p_mihandler->ScaleScheme()) {
  case 1: {
    Blob_List shower=bloblist->
      FindConnected(m_diced?bloblist->FindFirst(btp::Hard_Collision):
		    signal);
    int njets=2;
    if (!m_diced) {
      Blob_Data_Base *xsinfo=
	(*bloblist->FindFirst(btp::ME_PS_Interface_FS))["Core_Process"];
      if (xsinfo==NULL) THROW(critical_error,"No merging information.");
      XS_Base *xs=xsinfo->Get<XS_Base*>();
      if (xs!=NULL) if (xs->NStrong()<4) njets=1;
    }
    Particle_List jets=shower.ExtractLooseParticles(1);
    Jet_Finder finder(p_mihandler->YCut(),4);
    finder.ConstructJets(&jets,njets,true);
    double ptmax=0.0;
    for (Particle_List::const_iterator pit=jets.begin();
	 pit!=jets.end();++pit)
      ptmax=Max(ptmax,(*pit)->Momentum().PPerp());
    jets.Clear();
    m_ptmax=Min(m_ptmax,ptmax);
    break;
  }
  default: THROW(not_implemented,"Wrong mi scale scheme.");
  }
  Blob_List isr=bloblist->Find(btp::IS_Shower);
  for (Blob_List::const_iterator iit=isr.begin();iit!=isr.end();++iit) {
    m_emax[(*iit)->Beam()]-=(*iit)->InParticle(0)->Momentum()[0];
    p_mihandler->ISRHandler()->
      Extract((*iit)->InParticle(0)->Flav(),
	      (*iit)->InParticle(0)->Momentum()[0],(*iit)->Beam());
    if (!p_remnants[(*iit)->Beam()]->Extract((*iit)->InParticle(0))) {
      msg_Tracking()<<"Multiple_Interactions::CheckBlobList(..): "
		    <<"Cannot extract parton from hadron. \n"
		    <<*(*iit)->InParticle(0)<<std::endl;
      p_bloblist->DeleteConnected(*iit);
      if (bloblist->empty()) {
	Blob *blob = new Blob();
	blob->SetType(btp::Signal_Process);
	blob->SetStatus(-1);
	blob->SetId();
	blob->SetStatus(2);
	bloblist->push_back(blob);	  
      }
      return true;
    } 
  }
#ifdef ANALYSE__Multiple_Interactions
  if (!m_diced) {
    static TH2D *ept=NULL;
    if (ept==NULL) {
      ept = new TH2D("erem_pt","erem_pt",100,0.0,rpa.gen.Ecms(),
		     100,0.0,rpa.gen.Ecms());
      MYROOT::myroot->AddObject(ept,"erem_pt");
    }
    ept->Fill(2.0*m_ptmax,m_emax[0]+m_emax[1],
	      (*signal)["ME_Weight"]->Get<double>());
    static TH1D *mu0=NULL;
    if (mu0==NULL) {
      mu0 = new TH1D("mu_mi","mu_mi",200,0.0,log10(rpa.gen.Ecms()/2.0));
      MYROOT::myroot->AddObject(mu0,"mu_mi");
    }
    mu0->Fill(log10(m_ptmax),
	      (*signal)["ME_Weight"]->Get<double>());
  }
#endif
  return m_ptmax!=std::numeric_limits<double>::max();
}

bool Multiple_Interactions::Treat(ATOOLS::Blob_List *bloblist,double &weight)
{
  PROFILE_HERE;
  if (p_mihandler->Type()==MI_Handler::None ||
      MI_Base::StopGeneration()) return false;
  if (bloblist->empty()) {
    msg.Error()<<"Multiple_Interactions::Treat(): "
		       <<"Incoming blob list is empty!"<<std::endl;
    return false;
  }
  if (!CheckBlobList(bloblist)) return false;
  m_diced=true;
  p_mihandler->SetScaleMax(m_ptmax,0);
  p_mihandler->SetScaleMin(p_mihandler->ScaleMin(0),0);
  p_mihandler->SetScaleMax(m_emax[0],2);
  p_mihandler->SetScaleMax(m_emax[1],3);
  Blob *blob = new Blob();
  blob->AddData("MI_Scale",new Blob_Data<double>(m_ptmax));
  p_mihandler->Reset();
  bool success=false;
  if (!m_vetoed && m_ptmax>p_mihandler->ScaleMin(0)) {
    success=p_mihandler->GenerateHardProcess(blob);
  }
  else if (m_vetoed) {
    success=p_mihandler->GenerateSoftProcess(blob);
    // dummy settings for unweighted analysis
    blob->SetType(btp::Signal_Process);
    blob->SetTypeSpec("Soft UE");
    blob->SetStatus(5);
    blob->AddData("ME_Weight",new Blob_Data<double>(1.0));
    blob->AddData("ME_NumberOfTrials",new Blob_Data<int>(1));
  }
  if (success) {
    blob->SetId(bloblist->size());
    m_ptmax=blob->OutParticle(0)->Momentum().PPerp();
    for (size_t i=0;i<(size_t)blob->NInP();++i) {
      p_mihandler->ISRHandler()->Reset(i);
      if (!p_remnants[i]->Extract(blob->InParticle(i))) {
	msg_Tracking()<<"Multiple_Interactions::Treat(..): "
		      <<"Cannot extract parton from hadron. \n"
		      <<*blob->InParticle(0)<<std::endl;
	delete blob;
	return true;
      }
    }
    bloblist->push_back(blob);
    return true;
  }
  delete blob;
  p_mihandler->ISRHandler()->Reset(0);
  p_mihandler->ISRHandler()->Reset(1);
  return false;
}

bool Multiple_Interactions::VetoHardProcess(ATOOLS::Blob *const blob)
{
  if (p_mihandler->VetoHardProcess(blob)) {
    p_bloblist->DeleteConnected(blob);
    return m_vetoed=true;
  }
  return false;
}

void Multiple_Interactions::Finish(const std::string &resultpath) 
{
}

void Multiple_Interactions::CleanUp() 
{
  p_mihandler->CleanUp();
  m_ptmax=std::numeric_limits<double>::max();
  m_vetoed=false;
  m_diced=false;
}
