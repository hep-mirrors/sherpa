#include "Hard_Interface.H"
#include "Run_Parameter.H"
#include "Random.H"

using namespace APACIC;
using namespace AMEGIC;
using namespace EXTRAXS;
using namespace ISR;
using namespace AMATOOLS;
using namespace APHYTOOLS;
using namespace AORGTOOLS;


Hard_Interface::Hard_Interface(ISR_Handler *& _isr,int _maxjetnumber,bool outside) :
  isr(_isr), maxjetnumber(_maxjetnumber)
{
  msg.Debugging()<<"Passed isr : "<<isr<<std::endl;
  fin_tree                           = new Tree();
  ini_trees                          = new Tree*[2];
  for (int i=0;i<2;i++) ini_trees[i] = new Tree();

  fin_shower   = new Final_State_Shower();
  ini_shower   = new Initial_State_Shower(isr,fin_shower);
  _isr=isr;
  tools        = new Interface_Tools(ini_trees,fin_tree);


  if (outside) two2two = 0;
  else         two2two = new XS_Group(2,2,"Core process");
}
  
Hard_Interface::~Hard_Interface() 
{
  if (fin_tree)     { delete fin_tree;fin_tree = 0; }
  if (ini_trees)    {
    for (int i=0;i<2;i++) {
      delete ini_trees[i];ini_trees[i]=0;
    }
    delete ini_trees;ini_trees = 0;
  }  
  if (ini_shower)   { delete ini_shower;ini_shower = 0; }
  if (fin_shower)   { delete fin_shower;fin_shower = 0; }
  if (tools)        { delete tools; tools = 0; }
  if (two2two)      { delete two2two; two2two = 0; }

  msg.Tracking()<<"Out Hard_Interface::~Hard_Interface :"<<std::endl;
  msg.Tracking()<<"+++++++++++++++++++++++++++++++++++++++++++++++"<<std::endl;
}




void Hard_Interface::PrepareTrees() {
  msg.Tracking()<<"In Hard_Interface::PrepareTrees : Reset trees."<<std::endl;
  fin_tree->Reset();    
  for (int i=0;i<2;i++) ini_trees[i]->Reset();
}




bool Hard_Interface::ReduceBlob(Blob * blob) 
{
  if (!blob) return 0;
  blob->BoostInCMS();

  nin  = blob->NInP();
  nout = blob->NOutP();

  if (nin!=2) 
    msg.Error()<<"ERROR in Hard_Interface::ReduceBlob(Blob * blob)."<<std::endl
	       <<"   No idea how to handle blobs with "<<nin<<" incoming legs."<<std::endl;
  if (nout>2) {}
  if (nout==2) return InitColours(blob);
}

bool Hard_Interface::InitColours(Blob * blob) {
  Flavour * fl = new Flavour[4];
  vec4d   * p  = new vec4d[4];
  for (int i=0;i<2;i++) {
    fl[i]   = blob->InParton(i)->flav(); 
    p[i]    = blob->InParton(i)->momentum(); 
    fl[i+2] = blob->OutParton(i)->flav();
    p[i+2]  = blob->OutParton(i)->momentum(); 
  } 
  XS_Base * xs; 
  if (!(xsselector->FindInGroup(two2two,xs,2,2,fl))) {
    xs = xsselector->GetXS(2,2,fl);
    two2two->Add(xs);
  }
  if (!(xs->SetColours(p))) return 0;
  for (int j=0;j<2;j++) {
    for (int i=0;i<2;i++) {
      blob->InParton(i)->set_flow(j+1,xs->Colours()[i][j]);
      blob->OutParton(i)->set_flow(j+1,xs->Colours()[i+2][j]);
    }
  }
  return InitialConditions(blob,xs);
}

bool Hard_Interface::InitialConditions(Blob * blob,XS_Base * xs) {
  double scale = xs->Scale();
  double th1,th2;
  th1   = tools->ColourAngle(blob->InParton(0),blob);
  th2   = tools->ColourAngle(blob->InParton(1),blob);
  tools->InitializeIncoming(blob,scale,th1,th2);
  th1   = tools->ColourAngle(blob->OutParton(0),blob);
  th2   = tools->ColourAngle(blob->OutParton(1),blob);
  tools->InitializeOutGoing(blob,scale,th1,th2);
  return 1;
}




int Hard_Interface::PerformShowers(bool ini,bool fin) {
  bool jetveto = 1;   // *AS* default is one !
  int stat=1;
  if (nout==maxjetnumber) jetveto = 0; 
  if (fin) {
    stat=fin_shower->PerformShower(fin_tree,jetveto);
    if (!(stat)) return 0;
    msg.Tracking()<<"Final State Shower successful !"<<std::endl;
    msg.Tracking()<<"Has to be boosted into lab frame !"<<std::endl;
    fin_shower->SetColours(fin_tree->GetRoot());
    if (rpa.gen.Tracking()) fin_shower->OutputTree(fin_tree);
  }
  if (ini) {
    ini_shower->InitShowerPT(ini_trees[0]->GetRoot()->maxpt2);
    if (!(ini_shower->PerformShower(ini_trees,jetveto))) return 0;
    msg.Tracking()<<"Initial State Shower successful !"<<std::endl;
    msg.Tracking()<<"Has to be boosted into lab frame !"<<std::endl;
  }
  return stat;
}


bool Hard_Interface::ExtractPartons(bool ini,bool fin,Blob_List *bl,Parton_List * pl) {
  if (fin) {
    if (bl) fin_shower->ExtractPartons(fin_tree->GetRoot(),0,bl,pl);
    else    fin_shower->ExtractPartons(fin_tree->GetRoot(),pl);
    msg.Tracking()<<std::endl<<std::endl;
  }
  if (ini) {
    for (int i=0;i<2;i++) {
      msg.Tracking()<<" Now Extract Partons for incoming leg "<<i<<std::endl;
      ini_shower->ExtractPartons(ini_trees[i]->GetInitiator(),0,bl,pl,i);
    }
  }
  return 1;
}



void Hard_Interface::FinalStats() {
  msg.Debugging()<<"Initialised "<<two2two->Size()
		 <<" new XS for colours :"<<std::endl;
  for (int i=0;i<two2two->Size();i++)
    msg.Debugging()<<"   "<<(*two2two)[i]->Name()<<std::endl;
}

