//bof
//Version: 2 ADICIC++-0.0/2004/08/10

//Implementation of Adicic.H.



#include "Message.H"
#include "Particle_List.H"
#include "Chain_Handler.H"
#include "Paraminit.H"
#include "Adicic.H"





using namespace std;
using namespace ATOOLS;
using namespace ADICIC;





//=============================================================================



Adicic::Adicic() : p_chain(new Chain()) {
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
}





Adicic::Adicic(MODEL::Model_Base* pmod) : p_chain(new Chain()) {
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
}





Adicic::~Adicic() {
  if(p_chain) { delete p_chain; p_chain=NULL;}
}



//=============================================================================



int Adicic::PerformShowers() {
  Chain_Handler chainhandler(*p_chain);
  while(chainhandler.EvolveChainByOneStep()) {
    Vec4D test;
    p_chain->CheckMomentumConservation(test);
  }
  return 1;
}





bool Adicic::ExtractPartons(ATOOLS::Blob_List* blobs, ATOOLS::Particle_List*) {

  static unsigned count=0;
  ++count;
  if(count==1 || count%1000==0) {
    cout<<"\n                "<<om::red<<om::blackbg<<count<<om::reset<<"\n\n";
    for(int i=0; i<100000000; ++i);
  }

  int num;
  Blob* blob;
  Particle_List plist;

#ifdef ADICIC_OUTPUT
  cout<<"Blob list......: "<<blobs->size()<<endl;//////////////////////////////
  cout<<"Particle number: "<<ATOOLS::Particle::Counter()<<endl;////////////////
  cout<<"Current  number: "<<ATOOLS::Particle::CurrentNumber()<<"\n"<<endl;////
#endif

  for(Blob_Iterator blit=blobs->begin(); blit!=blobs->end(); ++blit) {

    if((*blit)->Type()==btp::FS_Shower) {
      if(!p_chain->ExtractPartons(plist)) {
	msg.Error()<<"Error in Adicic::ExtractPartons."<<endl;
	return false;
      }
      blob=*blit;
#ifdef ADICIC_OUTPUT
      cout<<om::greenbg<<"Initial 2nd blob:"<<om::reset<<" \n"<<*blob<<endl;///
#endif
      for(int i=0; i<blob->NInP(); ++i) blob->InParticle(i)->SetStatus(2);
      num=Max(blob->InParticle(0)->Number(),blob->InParticle(1)->Number());
      num=num+1-plist.front()->Number();
      for(Particle_Iterator piter=plist.begin(); piter!=plist.end(); ++piter) {
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
