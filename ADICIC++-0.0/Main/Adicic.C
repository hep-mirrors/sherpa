//bof
//Version: 1 ADICIC++-0.0/2004/06/08

//Implementation of Adicic.H.



#include "Message.H"
#include "Particle_List.H"
#include "Chain_Handler.H"
#include "Adicic.H"





using namespace std;
using namespace ATOOLS;
using namespace ADICIC;





//=============================================================================



Adicic::Adicic() : p_chain(new Chain()) {
  cout<<om::green;/////////////////////////////////////////////////////////////
  cout<<"Running?="<<Sudakov_Calculator::IsAlphaSRunning()<<endl;//////////////
  cout<<"MinScale="<<Sudakov_Calculator::MinOfK2t()<<endl;/////////////////////
  cout<<"MaxScale="<<Sudakov_Calculator::MaxOfK2t()<<endl;/////////////////////
  cout<<"ASFix="<<Sudakov_Calculator::AlphaSFix()<<endl;///////////////////////
  cout<<"ASApp="<<Sudakov_Calculator::AlphaSApprox()<<endl;////////////////////
  cout<<"ASCor="<<Sudakov_Calculator::AlphaSCorr(700.0)<<endl;/////////////////
  cout<<"ASCor="<<Sudakov_Calculator::AlphaSCorr(8100.0)<<endl;////////////////
  cout<<om::brown;/////////////////////////////////////////////////////////////
  cout<<"SudakovInit="<<Sudakov_Calculator::Init(NULL)<<endl;
  cout<<om::green;/////////////////////////////////////////////////////////////
  cout<<"Running?="<<Sudakov_Calculator::IsAlphaSRunning()<<endl;//////////////
  cout<<"ASFix="<<Sudakov_Calculator::AlphaSFix()<<endl;///////////////////////
  cout<<"ASApp="<<Sudakov_Calculator::AlphaSApprox()<<endl;////////////////////
  cout<<"ASCor="<<Sudakov_Calculator::AlphaSCorr(700.0)<<endl;/////////////////
  cout<<"ASCor="<<Sudakov_Calculator::AlphaSCorr(8100.0)<<endl;////////////////
  cout<<om::reset;/////////////////////////////////////////////////////////////
}





Adicic::Adicic(MODEL::Model_Base* pmod) : p_chain(new Chain()) {
  cout<<om::green;/////////////////////////////////////////////////////////////
  cout<<"Running?="<<Sudakov_Calculator::IsAlphaSRunning()<<endl;//////////////
  cout<<"MinScale="<<Sudakov_Calculator::MinOfK2t()<<endl;/////////////////////
  cout<<"MaxScale="<<Sudakov_Calculator::MaxOfK2t()<<endl;/////////////////////
  cout<<"ASFix="<<Sudakov_Calculator::AlphaSFix()<<endl;///////////////////////
  cout<<"ASApp="<<Sudakov_Calculator::AlphaSApprox()<<endl;////////////////////
  cout<<"ASCor="<<Sudakov_Calculator::AlphaSCorr(700.0)<<endl;/////////////////
  cout<<"ASCor="<<Sudakov_Calculator::AlphaSCorr(8100.0)<<endl;////////////////
  cout<<om::red;///////////////////////////////////////////////////////////////
  //cout<<"SudakovInit="<<Sudakov_Calculator::Init(NULL)<<endl;
  //cout<<"SudakovInit="<<Sudakov_Calculator::Init(pmod)<<endl;
  //cout<<"SudakovInit="<<Sudakov_Calculator::Init(NULL)<<endl;
  cout<<"SudakovInit="<<Sudakov_Calculator::Init(pmod)<<endl;
  cout<<om::green;/////////////////////////////////////////////////////////////
  cout<<"Running?="<<Sudakov_Calculator::IsAlphaSRunning()<<endl;//////////////
  cout<<"ASFix="<<Sudakov_Calculator::AlphaSFix()<<endl;///////////////////////
  cout<<"ASApp="<<Sudakov_Calculator::AlphaSApprox()<<endl;////////////////////
  cout<<"ASCor="<<Sudakov_Calculator::AlphaSCorr(700.0)<<endl;/////////////////
  cout<<"ASCor="<<Sudakov_Calculator::AlphaSCorr(8100.0)<<endl;////////////////
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
  int num;
  Blob* blob;
  Particle_List plist;
  cout<<"Blob list:"<<" "<<blobs->size()<<endl;////////////////////////////////
  for(Blob_Iterator blit=blobs->begin(); blit!=blobs->end(); ++blit) {
    if((*blit)->Type()==btp::FS_Shower) {
      if(!p_chain->ExtractPartons(plist)) {
	msg.Error()<<"Error in Adicic::ExtractPartons."<<endl;
	return false;
      }
      blob=*blit;
      for(int i=0; i<blob->NInP(); ++i) blob->InParticle(i)->SetStatus(2);
      num=Max(blob->InParticle(0)->Number(),blob->InParticle(1)->Number());
      num=num+1-plist.front()->Number();
      for(Particle_Iterator piter=plist.begin(); piter!=plist.end(); ++piter) {
	blob->AddToOutParticles(*piter);
	(*piter)->SetInfo('F');
	(*piter)->SetNumber((*piter)->Number()+num);
      }
      cout<<om::greenbg<<"Final Blob:"<<om::reset<<" "<<endl<<blob<<endl;//////
      return true;
    }
    cout<<om::brownbg<<"First Blob:"<<om::reset<<" "<<endl<<*blit<<endl;///////
  }
  return false;
}



//=============================================================================





//eof
