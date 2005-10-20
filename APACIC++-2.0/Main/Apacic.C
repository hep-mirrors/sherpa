#include "Apacic.H"

#include "Blob.H"
#include "Random.H"
#include "Run_Parameter.H"

using namespace APACIC;
using namespace PDF;
using namespace ATOOLS;


Apacic::Apacic(ISR_Handler *const isr,MODEL::Model_Base *const model,
	       ATOOLS::Jet_Finder *const jf,Data_Read *const dataread):
  p_inishower(NULL), p_finshower(NULL), 
  p_jetveto(NULL),
  p_initrees(NULL), p_fintree(NULL)
{
  Splitting_Function::SetKFactorScheme
    (dataread->GetValue<int>("S_KFACTOR_SCHEME",1));        
  m_fsron = bool(dataread->GetValue<int>("FSR_SHOWER",1));
  m_isron = bool(dataread->GetValue<int>("ISR_SHOWER",1));
  jv::mode mlm(dataread->GetValue<int>("MLM",0)==1?jv::mlm:jv::none);
  if (m_fsron) {
    p_fintree   = new Tree();
    p_finshower = new Final_State_Shower(model,jf,dataread);
    p_jetveto = new Jet_Veto(jf,p_finshower->Kinematics());
    p_jetveto->SetMode(mlm|jv::final);
    p_jetveto->SetFSTree(p_fintree);
    p_finshower->SetJetVeto(p_jetveto);
    m_showers=true;
  }
  if (m_isron) {
    if (isr->On()==0) {
      msg.Error()<<METHOD<<" Error: ISR must be enabled for ISR Shower."<<std::endl;
      abort();
    }
    p_initrees  = new Tree*[2];
    for (int i=0;i<2;i++) p_initrees[i] = new Tree();
    p_inishower = new Initial_State_Shower(isr,jf,p_finshower,model,dataread);
    m_showers=true;
    if (m_fsron) {
      p_jetveto->SetMode(p_jetveto->Mode()|jv::initial);
      p_jetveto->SetISTrees(p_initrees);
    }
  }
}
  
Apacic::~Apacic() 
{
  if (p_fintree)     { delete p_fintree; p_fintree = 0; }
  if (p_initrees)    {
    for (int i=0;i<2;i++) { delete p_initrees[i]; p_initrees[i] = 0; }
    delete p_initrees;p_initrees = 0;
  }  
  if (p_inishower)   { delete p_inishower; p_inishower = 0; }
  if (p_finshower)   { delete p_finshower; p_finshower = 0; }
}

int Apacic::PerformShowers(const int &jetveto,const int &losejv,
			   const double &x1,const double &x2,
			   const double &ycut) {
  if (!m_showers) return 1;
  if (msg.LevelIsDebugging()) {
    msg.Out()<<"############################################################"<<std::endl
	     <<"Apacic::PerformShowers : Before showering."<<std::endl;
    OutputTrees();
  }
  if (m_fsron) {
    Vec4D cms(PrepareFSR());
    if (!m_isron) p_jetveto->SetMode(jv::final);
    else p_jetveto->SetMode(jv::initial|jv::final);
    int fsrstatus(p_finshower->PerformShower(p_fintree,jetveto));
    switch (fsrstatus) {
    case -1:
      msg.Debugging()<<"Apacic::PerformShowers : "<<std::endl
		     <<"   Lose jet veto in FSR Shower."<<std::endl;
      return -1;
    case 0:
      msg.Error()<<"WARNING in Apacic::PerformShowers : "<<std::endl
		 <<"   FSR Shower failure. Abort."<<std::endl;
      //msg.SetLevel(0);
      return 0;
    }
    p_finshower->SetAllColours(p_fintree->GetRoot());
  }

  if (m_isron) {
    p_inishower->InitShowerPT(p_initrees[0]->GetRoot()->maxpt2); // ???
    int isrstatus(p_inishower->PerformShower(p_initrees,p_fintree,jetveto));
    switch (isrstatus) {
    case -1:
      msg.Debugging()<<"Apacic::PerformShowers : "<<std::endl
		     <<"   Lose jet veto in ISR Shower."<<std::endl;
      return -1;
    case 0:
      msg.Error()<<"WARNING in Apacic::PerformShowers : "<<std::endl
		 <<"   ISR Shower failure. Abort."<<std::endl;
      return 0;
    }
  }

  BoostInLab();

  p_fintree->CheckMomentumConservation();
  if (m_isron) {
    p_initrees[0]->CheckMomentumConservation();
    p_initrees[1]->CheckMomentumConservation();
  }

  switch (p_jetveto->TestKinematics(1)) {
  case 0:
    msg.Error()<<"WARNING in Apacic::PerformShowers : "<<std::endl
	       <<"   ISR Shower failure. Abort."<<std::endl;
    return 0;
  case -1: 
    msg_Debugging()<<"kinematics vetoed\n";
    return -1;
  }
  msg_Debugging()<<"kinematics check passed"<<std::endl;
  
  int number(0);
  Vec4D sum_fs=p_finshower->GetMomentum(p_fintree->GetRoot(),number);
  if (number<0) {
    msg.Error()<<"ERROR in Apacic::PerformShowers : "<<std::endl
	       <<"   Four Momentum not conserved, try new event."<<std::endl;
    //msg.SetLevel(0);
    return 0;
  }
  if (msg.LevelIsDebugging()) {
    msg.Out()<<"Apacic::PerformShowers : After showering."<<std::endl;
    OutputTrees();
    msg.Out()<<"############################################################"<<std::endl;
  }
  //msg.SetLevel(0);
  return 1;
}

bool Apacic::ExtractPartons(const bool ini,const bool fin,
			    Blob_List *const bl,Particle_List *const pl) {
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

//   Blob_List ble(bl->Find(btp::IS_Shower|btp::FS_Shower));
//   Particle_List jets(ble.ExtractLooseParticles(1));
//   PRINT_INFO(jets.size());
//   p_jetveto->JetFinder()->ConstructJets(&jets,
// 					p_jetveto->JetFinder()->Ycut(),true);
//   if (jets.size()>0) {
//     for (size_t i=0;i<jets.size();++i)
//       PRINT_INFO(i<<"-rate "<<jets[i]->Momentum().PPerp());
//   }
//   jets.Clear();

  return true;
}

void Apacic::BoostInLab() 
{
  if (!m_fsron) return;
  if (m_isron) {
    p_inishower->BoostFS();
  }
  else {
    Vec4D cms(p_fintree->GetRoot()->part->Momentum());
    cms=Vec4D(cms[0],-1.*Vec3D(cms));
    Poincare lab(cms);
    p_fintree->BoRo(lab);
  }
}

void Apacic::PrepareTrees() 
{
  if (m_fsron) { p_fintree->Reset(); } 
  if (m_isron) for (int i=0;i<2;i++) p_initrees[i]->Reset();
}

void Apacic::OutputTrees() 
{
  if (m_fsron) p_finshower->OutputTree(p_fintree);
  if (m_isron) {
    p_inishower->OutputTree(p_initrees[0]);
    p_inishower->OutputTree(p_initrees[1]);
  }
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

