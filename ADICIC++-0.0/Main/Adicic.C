//bof
//Version: 2 ADICIC++-0.0/2004/09/10

//Implementation of Adicic.H.



#include "Message.H"
#include "Particle_List.H"
#include "Paraminit.H"
#include "Adicic.H"





using namespace std;
using namespace ATOOLS;
using namespace ADICIC;





//=============================================================================



Adicic::Adicic() : m_total(0), m_fail(0), p_handler(NULL), m_cascade() {
  cout<<om::green;/////////////////////////////////////////////////////////////
  cout<<"Dipole Flavour Initialization ......."<<endl;/////////////////////////
  cout<<"  "<<Dipole_Flavour_Init::Status()<<endl;/////////////////////////////
  cout<<"  "<<Dipole_Flavour_Init::DoIt()<<endl;
  cout<<"  "<<Dipole_Flavour_Init::DoIt()<<endl;///////////////////////////////
  cout<<"  "<<Dipole_Flavour_Init::Status()<<endl;/////////////////////////////
  cout<<"Dipole Parameter Initialization .......\n\n";/////////////////////////
  cout<<"  Dip Param Init?: "<<Dipole_Parameter_Init::Status()<<endl;//////////
  cout<<"\n  Do Init: "<<Dipole_Parameter_Init::DoIt()<<"\n";
  cout<<"  Do Init: "<<Dipole_Parameter_Init::DoIt()<<"\n";////////////////////
  cout<<"  Dip Param Init?: "<<Dipole_Parameter_Init::Status()<<"\n\n";////////
  cout<<"Sudakov Calculator Initialization ......."<<endl;/////////////////////
  cout<<"  Running?="<<Sudakov_Calculator::IsAlphaSRunning()<<endl;////////////
  cout<<"  MinScale="<<Sudakov_Calculator::MinOfK2t()<<endl;///////////////////
  cout<<"  MaxScale="<<Sudakov_Calculator::MaxOfK2t()<<endl;///////////////////
  cout<<"  ASFix="<<Sudakov_Calculator::AlphaSFix()<<endl;/////////////////////
  cout<<"  ASApp="<<Sudakov_Calculator::AlphaSApprox()<<endl;//////////////////
  cout<<"  ASCor="<<Sudakov_Calculator::AlphaSCorr(700.0)<<endl;///////////////
  cout<<"  ASCor="<<Sudakov_Calculator::AlphaSCorr(8100.0)<<endl;//////////////
  cout<<om::brown;/////////////////////////////////////////////////////////////
  cout<<"  SudakovInit(NULL)="<<Sudakov_Calculator::Init(NULL)<<endl;
  cout<<om::green;/////////////////////////////////////////////////////////////
  cout<<"  Running?="<<Sudakov_Calculator::IsAlphaSRunning()<<endl;////////////
  cout<<"  ASFix="<<Sudakov_Calculator::AlphaSFix()<<endl;/////////////////////
  cout<<"  ASApp="<<Sudakov_Calculator::AlphaSApprox()<<endl;//////////////////
  cout<<"  ASCor="<<Sudakov_Calculator::AlphaSCorr(700.0)<<endl;///////////////
  cout<<"  ASCor="<<Sudakov_Calculator::AlphaSCorr(8100.0)<<endl;//////////////
  cout<<"Dipole Initializations ....... finished."<<endl;//////////////////////
  Dipole_Handler::ShowCalcBox();///////////////////////////////////////////////
  cout<<om::reset;/////////////////////////////////////////////////////////////
  //assert(0);

  p_handler=new Cascade_Handler();
  assert(p_handler);

}





Adicic::Adicic(MODEL::Model_Base* pmod)
  : m_total(0), m_fail(0), p_handler(NULL), m_cascade() {
  cout<<om::green;/////////////////////////////////////////////////////////////
  cout<<"Dipole Flavour Initialization ......."<<endl;/////////////////////////
  cout<<"  "<<Dipole_Flavour_Init::Status()<<endl;/////////////////////////////
  cout<<"  "<<Dipole_Flavour_Init::DoIt()<<endl;
  cout<<"  "<<Dipole_Flavour_Init::DoIt()<<endl;///////////////////////////////
  cout<<"  "<<Dipole_Flavour_Init::Status()<<endl;/////////////////////////////
  cout<<"Dipole Parameter Initialization .......\n\n";/////////////////////////
  cout<<"  Dip Param Init?: "<<Dipole_Parameter_Init::Status()<<endl;//////////
  cout<<"\n  Do Init: "<<Dipole_Parameter_Init::DoIt()<<"\n";
  cout<<"  Do Init: "<<Dipole_Parameter_Init::DoIt()<<"\n";////////////////////
  cout<<"  Dip Param Init?: "<<Dipole_Parameter_Init::Status()<<"\n\n";////////
  cout<<"Sudakov Calculator Initialization ......."<<endl;/////////////////////
  cout<<"  Running?="<<Sudakov_Calculator::IsAlphaSRunning()<<endl;////////////
  cout<<"  MinScale="<<Sudakov_Calculator::MinOfK2t()<<endl;///////////////////
  cout<<"  MaxScale="<<Sudakov_Calculator::MaxOfK2t()<<endl;///////////////////
  cout<<"  ASFix="<<Sudakov_Calculator::AlphaSFix()<<endl;/////////////////////
  cout<<"  ASApp="<<Sudakov_Calculator::AlphaSApprox()<<endl;//////////////////
  cout<<"  ASCor="<<Sudakov_Calculator::AlphaSCorr(700.0)<<endl;///////////////
  cout<<"  ASCor="<<Sudakov_Calculator::AlphaSCorr(8100.0)<<endl;//////////////
  cout<<om::red;///////////////////////////////////////////////////////////////
  //cout<<"SudakovInit="<<Sudakov_Calculator::Init(NULL)<<endl;
  //cout<<"SudakovInit="<<Sudakov_Calculator::Init(pmod)<<endl;
  //cout<<"SudakovInit="<<Sudakov_Calculator::Init(NULL)<<endl;
  cout<<"  SudakovInit(MODEL)="<<Sudakov_Calculator::Init(pmod)<<endl;
  cout<<om::green;/////////////////////////////////////////////////////////////
  cout<<"  Running?="<<Sudakov_Calculator::IsAlphaSRunning()<<endl;////////////
  cout<<"  ASFix="<<Sudakov_Calculator::AlphaSFix()<<endl;/////////////////////
  cout<<"  ASApp="<<Sudakov_Calculator::AlphaSApprox()<<endl;//////////////////
  cout<<"  ASCor="<<Sudakov_Calculator::AlphaSCorr(700.0)<<endl;///////////////
  cout<<"  ASCor="<<Sudakov_Calculator::AlphaSCorr(8100.0)<<endl;//////////////
  cout<<"Dipole Initializations ....... finished."<<endl;//////////////////////
  Dipole_Handler::ShowCalcBox();///////////////////////////////////////////////
  cout<<om::reset;/////////////////////////////////////////////////////////////
  //assert(0);

  p_handler=new Cascade_Handler();
  assert(p_handler);

}





Adicic::~Adicic() {
  m_cascade|0;
  p_handler->PrintCounter();
  cout<<"Total number = "<<m_total<<"\tTotal number of failures = "<<m_fail;
  cout<<"   ("<<1.0*m_fail/m_total<<")."<<endl;
  if(p_handler) delete p_handler;
}



//=============================================================================



int Adicic::PerformShowers() {
  ++m_total;
  assert(m_cascade|(*p_handler));
  if(p_handler->EvolveCascade());
  else ++m_fail;
  m_cascade|0;
  return 1;
}





bool Adicic::ExtractPartons(ATOOLS::Blob_List* blobs, ATOOLS::Particle_List*) {

  if(m_total==1 || m_total%1000==0)
    cout<<"                "<<om::red<<om::blackbg<<m_total<<om::reset<<"\n";

  int num;
  Blob* blob;
  list<Particle_List> lists;

#ifdef ADICIC_OUTPUT
  cout<<"Blob list......: "<<blobs->size()<<endl;//////////////////////////////
  cout<<"Particle number: "<<ATOOLS::Particle::Counter()<<endl;////////////////
  cout<<"Current  number: "<<ATOOLS::Particle::CurrentNumber()<<"\n"<<endl;////
#endif

  for(Blob_List::iterator blit=blobs->begin(); blit!=blobs->end(); ++blit) {

    if((*blit)->Type()==btp::FS_Shower) {
      if(!m_cascade.ExtractPartons(lists)) {
	msg.Error()<<__PRETTY_FUNCTION__<<endl;
	return false;
      }
      assert(lists.size()==1);/////////////////////////////////////////////////
      Particle_List& plist=lists.front();
      blob=*blit;
#ifdef ADICIC_OUTPUT
      cout<<om::greenbg<<"Initial 2nd blob:"<<om::reset<<" \n"<<*blob<<endl;///
#endif
      for(int i=0; i<blob->NInP(); ++i) blob->InParticle(i)->SetStatus(2);
      num=Max(blob->InParticle(0)->Number(),blob->InParticle(1)->Number());
      num=num+1-plist.front()->Number();
      for(Particle_List::iterator piter=plist.begin(); piter!=plist.end(); ++piter) {
	blob->AddToOutParticles(*piter);
	(*piter)->SetInfo('F');
	(*piter)->SetNumber(-((*piter)->Number()+num/*-2*/));
      }
#ifdef ADICIC_OUTPUT
      cout<<om::greenbg<<"Final 2nd blob:"<<om::reset<<" \n"<<*blob<<endl;/////
      cout<<"Particle number: "<<ATOOLS::Particle::Counter()<<endl;////////////
      cout<<"Current  number: "<<ATOOLS::Particle::CurrentNumber()<<"\n"<<endl;
#endif
      return true;
    }

#ifdef ADICIC_OUTPUT
    cout<<om::brownbg<<"Initial 1st blob:"<<om::reset<<" \n"<<**blit<<endl;////
#endif
    blob=*blit;
    for(int i=0; i<blob->NInP(); ++i)
      blob->InParticle(i)->SetNumber(-(1+i));
    for(int i=0; i<blob->NOutP(); ++i)
      blob->OutParticle(i)->SetNumber(-(blob->NInP()+1+i));
#ifdef ADICIC_OUTPUT
    cout<<om::brownbg<<"Final 1st blob:"<<om::reset<<" \n"<<**blit<<endl;//////
#endif
  }

  return false;

}



//=============================================================================





//eof
