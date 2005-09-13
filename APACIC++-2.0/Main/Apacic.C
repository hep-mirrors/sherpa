#include "Apacic.H"

#include "Blob.H"
#include "Random.H"
#include "Run_Parameter.H"

using namespace APACIC;
using namespace PDF;
using namespace ATOOLS;


Apacic::Apacic(ISR_Handler * isr,MODEL::Model_Base * model,
	       ATOOLS::Jet_Finder *jf,Data_Read * dataread):
  p_inishower(NULL), p_finshower(NULL), p_initrees(NULL), p_fintree(NULL)
{
  Splitting_Function::SetKFactorScheme
    (dataread->GetValue<int>("S_KFACTOR_SCHEME",1));        
  m_fsron = bool(dataread->GetValue<int>("FSR_SHOWER",1));
  m_isron = bool(dataread->GetValue<int>("ISR_SHOWER",1));
  if (m_fsron) {
    p_fintree   = new Tree();
    p_finshower = new Final_State_Shower(model,jf,dataread);
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

int Apacic::PerformShowers(int jetveto,int losejv,double x1,double x2, double ycut) {
  if (!m_showers) return 1;
  //msg.SetLevel(15);
  if (msg.LevelIsDebugging()) {
    msg.Out()<<"############################################################"<<std::endl
	     <<"Apacic::PerformShowers : Before showering."<<std::endl;
    OutputTrees();
  }
  if (m_fsron) {
    Vec4D cms(PrepareFSR());
    int fsrstatus(p_finshower->PerformShower(p_fintree,jetveto));
    switch (fsrstatus) {
    case -1:
      msg.Debugging()<<"Apacic::PerformShowers : "<<std::endl
		     <<"   Lose jet veto failed in FSR Shower, try new event."<<std::endl;
      return -1;
    case 0:
      msg.Error()<<"WARNING in Apacic::PerformShowers : "<<std::endl
		 <<"   FSR Shower did not work out, try new event."<<std::endl;
      //msg.SetLevel(0);
      return 0;
    }
    p_finshower->SetAllColours(p_fintree->GetRoot());
  }

  if (m_isron) {
    p_inishower->InitShowerPT(p_initrees[0]->GetRoot()->maxpt2); // ???
    int isrstatus(p_inishower->PerformShower(p_initrees,jetveto));
    switch (isrstatus) {
    case -1:
      msg.Debugging()<<"Apacic::PerformShowers : "<<std::endl
		     <<"   Lose jet veto failed in ISR Shower, try new event."<<std::endl;
    case 0:
      msg.Error()<<"WARNING in Apacic::PerformShowers : "<<std::endl
		 <<"   ISR Shower did not work out, try new event."<<std::endl;
      return isrstatus;
    }
  }

  BoostInCorrectSystem();
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

bool Apacic::ExtractPartons(bool ini,bool fin,Blob_List * bl,Particle_List * pl) {
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

void Apacic::BoostInCorrectSystem() {
  if (!m_fsron) return;
  if (!m_isron) {
    Vec4D cms(p_fintree->GetRoot()->part->Momentum());
    cms=Vec4D(cms[0],-1.*Vec3D(cms));
    Poincare lab(cms);
    p_fintree->BoRo(lab);
  }
  else {
    Vec4D mom1(p_initrees[0]->GetRoot()->part->Momentum());
    Vec4D mom2(p_initrees[1]->GetRoot()->part->Momentum());
    
    Vec4D vl(mom1[0]+mom2[0], -1.*Vec3D(mom1+mom2));
    Poincare lab(vl);
    lab.BoostBack(mom1);
    Poincare rot(Vec4D::ZVEC,mom1);
    p_fintree->BoRo(rot);
    p_fintree->BoRo(lab);
  }
}
