#include "Hard_Interface.H"
#include "Initial_State_Shower.H"
#include "Final_State_Shower.H"
#include "Tree.H"

#include "ISR_Handler.H"

#include "Run_Parameter.H"
#include "Random.H"

using namespace APACIC;
using namespace ISR;
using namespace AMATOOLS;
using namespace APHYTOOLS;
using namespace AORGTOOLS;


Hard_Interface::Hard_Interface(ISR_Handler * _isr,int _maxjetnumber,
			       bool _isron,bool _fsron,Data_Read * _dataread):
  m_isron(_isron), m_fsron(_fsron), m_showers(_isron||_fsron),
  p_fintree(NULL), p_finshower(NULL), p_initrees(NULL), p_inishower(NULL)
{
  msg.Debugging()<<"Passed isr : "<<_isr<<std::endl;
  if (m_fsron) {
    p_fintree   = new Tree();
    p_finshower = new Final_State_Shower(_dataread);
  }
  if (m_isron) {
    p_initrees  = new Tree*[2];
    for (int i=0;i<2;i++) p_initrees[i] = new Tree();
    p_inishower = new Initial_State_Shower(_isr,p_finshower,_dataread);
  }
}
  
Hard_Interface::~Hard_Interface() 
{
  if (p_fintree)     { delete p_fintree; p_fintree = 0; }
  if (p_initrees)    {
    for (int i=0;i<2;i++) { delete p_initrees[i]; p_initrees[i] = 0; }
    delete p_initrees;p_initrees = 0;
  }  
  if (p_inishower)   { delete p_inishower; p_inishower = 0; }
  if (p_finshower)   { delete p_finshower; p_finshower = 0; }

  msg.Tracking()<<"Out Hard_Interface::~Hard_Interface :"<<std::endl;
  msg.Tracking()<<"+++++++++++++++++++++++++++++++++++++++++++++++"<<std::endl;
}

void Hard_Interface::PrepareTrees() {
  msg.Tracking()<<"In Hard_Interface::PrepareTrees : Reset trees."<<std::endl;
  if (m_fsron) p_fintree->Reset();    
  if (m_isron) for (int i=0;i<2;i++) p_initrees[i]->Reset();
}

int Hard_Interface::PerformShowers(bool ini,bool fin,bool jetveto) {
  msg.Debugging()<<"In Hard_Interface::PerformShowers("<<ini<<","<<fin<<","<<jetveto<<")."<<endl;
  if (!m_showers) return 1;
  if (m_fsron) {
    int fsrstatus = p_finshower->PerformShower(p_fintree,jetveto);
    if (fsrstatus!=1) return fsrstatus;
    msg.Tracking()<<"Final State Shower successful !"<<std::endl
		  <<"Has to be boosted into lab frame !"<<std::endl;
    p_finshower->SetAllColours(p_fintree->GetRoot());
    if (rpa.gen.Tracking()) p_finshower->OutputTree(p_fintree);
  }
  if (m_isron) {
    p_inishower->InitShowerPT(p_initrees[0]->GetRoot()->maxpt2);
    if (!(p_inishower->PerformShower(p_initrees,jetveto))) return 0;
    msg.Tracking()<<"Initial State Shower successful !"<<std::endl
		  <<"Has to be boosted into lab frame !"<<std::endl;
  }
  return 1;
}

bool Hard_Interface::ExtractPartons(bool ini,bool fin,Blob_List * bl,Parton_List * pl) {
  if (fin) p_finshower->ExtractPartons(p_fintree->GetRoot(),0,bl,pl);
  if (ini) {
    for (int i=0;i<2;i++) {
      p_inishower->ExtractPartons(p_initrees[i]->GetInitiator(),i,0,bl,pl);
    }
  }
  return 1;
}

