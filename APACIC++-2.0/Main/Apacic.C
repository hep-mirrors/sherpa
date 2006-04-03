#include "Apacic.H"

#include "Blob.H"
#include "Random.H"
#include "Run_Parameter.H"
#include "Veto_Info.H"

#ifdef PROFILE__all
#define PROFILE__Apacic
#endif
#ifdef PROFILE__Apacic
#include "prof.hh"
#else 
#define PROFILE_HERE
#define PROFILE_LOCAL(LOCALNAME)
#endif

using namespace APACIC;
using namespace PDF;
using namespace ATOOLS;

Apacic::Apacic(ISR_Handler *const isr,MODEL::Model_Base *const model,
	       ATOOLS::Jet_Finder *const jf,Data_Read *const dataread):
  p_inishower(NULL), p_finshower(NULL), 
  p_jetveto(NULL),
  p_initrees(NULL), p_fintree(NULL),
  m_showers(false)
{
  Splitting_Function::SetKFactorScheme
    (dataread->GetValue<int>("S_KFACTOR_SCHEME",1));        
  m_fsron=bool(dataread->GetValue<int>("FSR_SHOWER",1));
  m_isron=bool(dataread->GetValue<int>("ISR_SHOWER",1));
  if ((rpa.gen.Beam1().IsHadron() || rpa.gen.Beam2().IsHadron())
      && (m_fsron^m_isron)) 
    THROW(fatal_error,"Shower must be enabled for hadronic initial state.");
  jv::mode jvm((jv::mode)dataread->GetValue<int>("JET_VETO_SCHEME",11));
  if (m_fsron) {
    p_fintree   = new Tree();
    p_finshower = new Final_State_Shower(model,jf,dataread);
    p_jetveto = new Jet_Veto(jf,p_finshower->Kinematics());
    p_jetveto->SetLoseJetMode(dataread->GetValue<int>("LOSE_JET_MODE",0));
    p_jetveto->SetMode(jvm);
    p_jetveto->SetFSTree(p_fintree);
    p_finshower->SetJetVeto(p_jetveto);
    m_showers=true;
  }
  if (m_isron) {
    if (isr->On()==0) THROW(fatal_error,"ISR must be enabled.");
    p_initrees  = new Tree*[2];
    for (int i=0;i<2;i++) p_initrees[i] = new Tree();
    p_inishower = new Initial_State_Shower(isr,jf,p_finshower,model,dataread);
    if (m_fsron) p_jetveto->SetISTrees(p_initrees);
    m_showers=true;
  }
}
  
Apacic::~Apacic() 
{
  if (p_fintree) delete p_fintree;
  if (p_initrees) {
    for (short unsigned int i(0);i<2;++i) delete p_initrees[i];
    delete p_initrees;
  }  
  if (p_inishower) delete p_inishower;
  if (p_finshower) delete p_finshower;
  delete p_jetveto;
}

int Apacic::PerformShowers(const int &jetveto,const int &losejv,
			   const double &x1,const double &x2,
			   const double &ycut) 
{
  PROFILE_HERE;
  if (!m_showers) return 1;
  p_jetveto->SetYCut(ycut);
  p_jetveto->SetJetVeto(jetveto);
  p_jetveto->SetLoseJetVeto(losejv);
  if (m_isron) {
    p_initrees[0]->Store();
    p_initrees[1]->Store();
  }
  if (m_fsron) p_fintree->Store();
  static size_t m_maxtrials(100);
  static double rej(0.0), cnt(0.0);
  ++cnt;
  size_t trials(0);
  for (;trials<m_maxtrials;++trials) {
  if (msg.LevelIsDebugging()) {
    msg.Out()<<"Apacic::PerformShowers : Before showering."<<std::endl;
    OutputTrees();
  }
  static double accu(sqrt(rpa.gen.Accu()));
  Vec4D::SetAccu(accu);
  if (m_fsron) {
    m_cms=PrepareFSR();
    switch (p_finshower->PerformShower(p_fintree,jetveto)) {
    case -1:
      msg.Debugging()<<"Lose jet veto in FSR Shower.\n";
      m_last_ljv=true;
      Vec4D::ResetAccu();
      return -1;
    case 0:
      msg.Error()<<METHOD<<"(..): FSR shower failure. Abort."<<std::endl;
      Vec4D::ResetAccu();
      return 0;
    }
    p_finshower->SetAllColours(p_fintree->GetRoot());
  }
  if (m_isron) {
    p_inishower->InitShowerPT(p_initrees[0]->GetRoot()->maxpt2);
    switch (p_inishower->PerformShower(p_initrees,p_fintree,jetveto)) {
    case -1:
      msg.Debugging()<<"Lose jet veto in ISR Shower.\n";
      m_last_ljv=true;
      Vec4D::ResetAccu();
      return -1;
    case 0:
      msg.Error()<<METHOD<<"(..): ISR shower failure. Abort."<<std::endl;
      Vec4D::ResetAccu();
      return 0;
    }
  }
  switch (p_jetveto->TestKinematics(2)) {
  case 1:
    msg_Debugging()<<"passed\n";
    trials=2*m_maxtrials;
    break;
  case 0:
    msg_Debugging()<<"Jet veto\n";
    if (m_isron) {
      p_initrees[0]->Restore();
      p_initrees[1]->Restore();
    }
    if (m_fsron) p_fintree->Restore();
    break;
  case -1: 
    msg_Debugging()<<"Lose jet veto\n";
    if (m_isron) {
      p_initrees[0]->ClearStore();
      p_initrees[1]->ClearStore();
    }
    if (m_fsron) p_fintree->ClearStore();
    m_last_ljv=true;
    Vec4D::ResetAccu();
    return -1;
  }
  }
  if (m_isron) {
    p_initrees[0]->ClearStore();
    p_initrees[1]->ClearStore();
  }
  if (m_fsron) p_fintree->ClearStore();
  Vec4D::ResetAccu();
  if (trials==m_maxtrials) {
    ++rej;
    if (rej/cnt>0.1) 
      msg.Error()<<METHOD<<"(): rej. rate = "<<rej/cnt<<".\n";
    // m_maxtrials*=2;
    return -1;
  }
  BoostInLab();
  p_fintree->CheckMomentumConservation();
  if (m_isron) {
    p_initrees[0]->CheckMomentumConservation();
    p_initrees[1]->CheckMomentumConservation();
  }
  msg_Debugging()<<"kinematics check passed"<<std::endl;
  int number(0);
  Vec4D sum_fs(p_finshower->GetMomentum(p_fintree->GetRoot(),number));
  Vec4D::ResetAccu();
  if (number<0) {
    msg.Error()<<METHOD<<"(..): Four Momentum not conserved. Abort."
	       <<std::endl;
    return 0;
  }
#ifdef USING__Veto_Info
  if (m_last_ljv) 
    p_finshower->Sudakov()->Vetos(0).front()=
      p_finshower->Sudakov()->Vetos(0).front()|svc::lj_veto;
#endif
  m_last_ljv=false;
  if (msg.LevelIsDebugging()) {
    msg.Out()<<"Apacic::PerformShowers : After showering."<<std::endl;
    OutputTrees();
  }
  return 1;
}

void Apacic::PrepareTrees() 
{
  if (m_fsron) p_fintree->Reset();
  if (m_isron) for (int i=0;i<2;i++) p_initrees[i]->Reset();
}

void Apacic::BoostInLab() 
{
  if (!m_fsron) return;
  if (m_isron) {
    p_inishower->BoostFS();
  }
  else {
    Poincare lab(m_cms);
    lab.Invert();
    p_fintree->BoRo(lab);
  }
}

bool Apacic::ExtractPartons(const bool ini,const bool fin,
			    Blob_List *const bl,Particle_List *const pl) 
{
  if (fin) {
    if (p_fintree->CheckStructure(true)) 
      p_finshower->ExtractPartons(p_fintree->GetRoot(),0,bl,pl);
    else return false;
  }
  if (ini) {
    for (int i=0;i<2;i++) {
      if (p_initrees[i]->CheckStructure(true)) 
	p_inishower->ExtractPartons(p_initrees[i]->GetInitiator(),i,0,bl,pl);
      else return false;
    }
  }
  return true;
}

ATOOLS::Vec4D Apacic::PrepareFSR() 
{
  ATOOLS::Vec4D cms(p_fintree->GetRoot()->part->Momentum());
  ATOOLS::Poincare cmsb(cms);
  p_fintree->BoRo(cmsb);
  return cms;
}

void Apacic::SetMaxJetNumber(const size_t &maxjets) 
{ 
  if (p_jetveto!=NULL) p_jetveto->SetMaxJets(maxjets); 
}

void Apacic::SetJetvetoPt2(const double &q2i, const double &q2f)
{ 
  if (p_jetveto!=NULL) p_jetveto->SetJetPT2(q2f); 
}

void Apacic::OutputTrees() 
{
  if (m_fsron) p_finshower->OutputTree(p_fintree);
  if (m_isron) {
    p_inishower->OutputTree(p_initrees[0]);
    p_inishower->OutputTree(p_initrees[1]);
  }
}

