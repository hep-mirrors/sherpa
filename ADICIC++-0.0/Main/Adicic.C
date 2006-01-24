//bof
//Version: 3 ADICIC++-0.0/2005/08/05

//Implementation of Adicic.H.



#include "Message.H"
#include "Particle_List.H"
#include "Paraminit.H"
#include "Adicic.H"





using namespace std;
using namespace ATOOLS;
using namespace ADICIC;





//=============================================================================



Adicic::Adicic()
  : m_parcountcorr(2),
    m_total(0), m_fail(0), m_noem(0), p_handler(NULL), m_cascade(),
    p_blist(NULL), l_paic() {
  cout<<om::brown;/////////////////////////////////////////////////////////////
  cout<<"Dipole Flavour Initialization .......\n  ";///////////////////////////
  Dipole_Flavour_Init dfi;
  cout<<"Dipole Parameter Initialization .......\n  ";/////////////////////////
  Dipole_Parameter_Init dpi;
  cout<<"Sudakov Calculator Initialization happened with NULL ......."<<endl;//
  cout<<"  Running?="<<Sudakov_Calculator::IsAlphaSRunning()<<endl;////////////
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





Adicic::Adicic(const string& path, MODEL::Model_Base* pmod)
  : m_parcountcorr(2),
    m_total(0), m_fail(0), m_noem(0), p_handler(NULL), m_cascade(),
    p_blist(NULL), l_paic() {
  cout<<om::green;/////////////////////////////////////////////////////////////
  cout<<"Dipole Flavour Initialization .......\n  ";///////////////////////////
  Dipole_Flavour_Init dfi;
  cout<<"Dipole Parameter Initialization .......\n  ";/////////////////////////
  Dipole_Parameter_Init dpi;
  cout<<"Sudakov Calculator Initialization so far ......."<<endl;//////////////
  cout<<"  Running?="<<Sudakov_Calculator::IsAlphaSRunning()<<endl;////////////
  cout<<"  ASApp="<<Sudakov_Calculator::AlphaSApprox()<<endl;//////////////////
  cout<<"  ASCor="<<Sudakov_Calculator::AlphaSCorr(700.0)<<endl;///////////////
  cout<<"  ASCor="<<Sudakov_Calculator::AlphaSCorr(8100.0)<<endl;//////////////
  cout<<om::red;///////////////////////////////////////////////////////////////
  //cout<<"SudakovInit="<<Sudakov_Calculator::AdjustEnvironment(pmod)<<endl;
  //cout<<"SudakovInit="<<Sudakov_Calculator::AdjustEnvironment(NULL)<<endl;
  cout<<"  SudakovInit(MODEL)="
      <<Sudakov_Calculator::AdjustEnvironment(path,pmod);
  cout<<endl<<om::green;///////////////////////////////////////////////////////
  cout<<"  Running?="<<Sudakov_Calculator::IsAlphaSRunning()<<endl;////////////
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
  cout<<"Total number="<<m_total
      <<".   Total number of no emissions="<<m_noem
      <<" ("<<1.0*m_noem/m_total<<")."
      <<"   Total number of failures="<<m_fail
      <<" ("<<1.0*m_fail/m_total<<")."<<endl;
  if(p_handler) delete p_handler;
}



//=============================================================================



int Adicic::PerformShowers() {
  ++m_total;
  //if(m_cascade.INumber()) m_cascade.Print();/////////////////////////////////
  assert(m_cascade|(*p_handler));
  if(p_handler->EvolveCascade()); else ++m_fail;
  m_cascade|0;
  if(m_cascade.DipoleNumber()==1) ++m_noem;
  //if(m_cascade.INumber()) m_cascade.Print();/////////////////////////////////
  return 1;
}





bool Adicic::ExtractPartons(ATOOLS::Blob_List* blobs) {

  assert(blobs);
  p_blist=blobs;

  if(m_total==1 || m_total%1000==0)
    cout<<"                "<<om::red<<om::blackbg<<m_total<<om::reset<<"\n";

  if(!m_cascade.ExtractPartons(l_paic)) {
    msg.Error()<<__PRETTY_FUNCTION__<<": Error: "
	       <<"Cannot extract partons from cascade!\n"<<endl;
    return false;
  }
  assert(l_paic.size()==1);////////////////////////////////////////////////////

#ifdef ADICIC_OUTPUT
  cout<<"Blob list......: "<<blobs->size()<<endl;//////////////////////////////
  cout<<"Particle number: "<<ATOOLS::Particle::Counter()<<endl;////////////////
  cout<<"Current  number: "<<ATOOLS::Particle::CurrentNumber()<<"\n"<<endl;////
#endif

  Blob* blob=NULL, * pis=NULL, * pisme=NULL, * pfs=NULL;
  for(Blob_List::iterator blit=p_blist->begin(); blit!=p_blist->end();
      ++blit) {
    switch((*blit)->Type()) {
    case btp::Signal_Process    : blob=*blit; break;
    case btp::IS_Shower         : pis=*blit; break;
    case btp::ME_PS_Interface_IS: pisme=*blit; break;
    case btp::FS_Shower         : pfs=*blit; break;
    default: assert(0);
    }
  }
  assert(blob);
  assert(pfs);

#ifdef ADICIC_OUTPUT
  cout<<om::brownbg<<"Initial 1st blob:"<<om::reset<<" \n"<<*blob<<endl;///////
#endif
  for(int i=0; i<blob->NInP(); ++i)
    blob->InParticle(i)->SetNumber(-(1+i));
  for(int i=0; i<blob->NOutP(); ++i)
    blob->OutParticle(i)->SetNumber(-(blob->NInP()+1+i));
#ifdef ADICIC_OUTPUT
  cout<<om::brownbg<<"Final 1st blob:"<<om::reset<<" \n"<<*blob<<endl;/////////
#endif

  bool result=false;
  if(m_cascade.INumber()==0) {
    assert(p_blist->size()==2);
    result=FExtract(pfs);
  }
  else if(m_cascade.INumber()==2) {
    assert(p_blist->size()==4);
    assert(pis);
    assert(pisme);
    result=IExtract(pis,pisme,pfs);
  }
  else {
    msg.Error()<<__PRETTY_FUNCTION__<<": Error: "
	       <<"Cannot handle this cascade!\n"<<endl;
    return false;
  }

#ifdef ADICIC_OUTPUT
  cout<<"Particle number: "<<ATOOLS::Particle::Counter()<<endl;////////////////
  cout<<"Current  number: "<<ATOOLS::Particle::CurrentNumber()<<"\n"<<endl;////
#endif

  return result;

}



//=============================================================================



bool Adicic::FExtract(Blob* pfs) {

  Particle_List& plist=l_paic.front();
#ifdef ADICIC_OUTPUT
  cout<<om::greenbg<<"Initial 2nd blob:"<<om::reset<<" \n"<<*pfs<<endl;////////
#endif
  for(int i=0; i<pfs->NInP(); ++i) pfs->InParticle(i)->SetStatus(2);
  int num=Max(pfs->InParticle(0)->Number(),pfs->InParticle(1)->Number());
  num=num+1-plist.front()->Number();
  for(Particle_List::iterator piter=plist.begin(); piter!=plist.end();
      ++piter) {
    pfs->AddToOutParticles(*piter);
    (*piter)->SetInfo('F');
    (*piter)->SetNumber(-((*piter)->Number()+num-m_parcountcorr));
  }
#ifdef ADICIC_OUTPUT
  cout<<om::greenbg<<"Final 2nd blob:"<<om::reset<<" \n"<<*pfs<<endl;//////////
#endif
  return true;

}





bool Adicic::IExtract(Blob* pis, Blob* pisme, Blob* pfs) {

  Particle_List& plist=l_paic.front();

#ifdef ADICIC_OUTPUT
  cout<<om::bluebg<<"Initial IS blob:"<<om::reset<<" \n"<<*pis<<endl;//////////
  cout<<om::redbg<<"Initial ISME blob:"<<om::reset<<" \n"<<*pisme<<endl;///////
  cout<<om::greenbg<<"Initial FS blob:"<<om::reset<<" \n"<<*pfs<<endl;/////////
#endif

  for(int i=0; i<pfs->NInP(); ++i) pfs->InParticle(i)->SetStatus(2);
  int num=Max(pfs->InParticle(0)->Number(),pfs->InParticle(1)->Number());
  num=num+1-plist.front()->Number();
  for(Particle_List::iterator piter=plist.begin(); piter!=plist.end();
      ++piter) {
    if((*piter)->Info()=='i') {
      pis->AddToInParticles(*piter);
      (*piter)->SetInfo('I');
      (*piter)->SetNumber(-((*piter)->Number()+num-2*m_parcountcorr));
    } else {
      pis->AddToOutParticles(*piter);
      (*piter)->SetInfo('F');
      (*piter)->SetNumber(-((*piter)->Number()+num-m_parcountcorr));
    }
  }
  if(pis->InParticle(0)->Momentum()[3]<pis->InParticle(1)->Momentum()[3])
    pis->SwapInParticles(0,1);

  Particle* par=new Particle(*pisme->OutParticle(0));
  Particle* per=new Particle(*pisme->OutParticle(1));
  pis->AddToOutParticles(par);
  pis->AddToOutParticles(per);
  pisme->AddToInParticles(par);
  pisme->AddToInParticles(per);

  Particle* pic=new Particle(*pfs->InParticle(0));
  Particle* poc=new Particle(*pfs->InParticle(1));
  pfs->AddToOutParticles(pic);
  pfs->AddToOutParticles(poc);
  pic->SetStatus(1);
  poc->SetStatus(1);

  const Vec4D PZ=par->Momentum()+per->Momentum();
  assert(PZ==pic->Momentum()+poc->Momentum());
  const Vec4D& PZprime=-1*m_cascade.Momentum();
  assert(dabs(PZ.Abs2()-PZprime.Abs2())<1.0e-7);
  Poincare fly(PZ);
  Poincare flyprime(PZprime);
  Vec4D p[4]={
    par->Momentum(), per->Momentum(), pic->Momentum(), poc->Momentum()};
  for(size_t i=0; i<4; ++i) {
    fly.Boost(p[i]); flyprime.BoostBack(p[i]);
  }
  par->SetMomentum(p[0]);
  per->SetMomentum(p[1]);
  pic->SetMomentum(p[2]);
  poc->SetMomentum(p[3]);

  pis->SetBeam(0);
  pisme->SetStatus(0);

#ifdef ADICIC_OUTPUT
  p[1]=pis->CheckMomentumConservation();
  p[2]=pisme->CheckMomentumConservation();
  p[3]=pfs->CheckMomentumConservation();
  cout<<om::bluebg<<"Final IS blob:"<<om::reset<<" \n"<<*pis<<p[1]<<endl;//////
  cout<<om::redbg<<"Final ISME blob:"<<om::reset<<" \n"<<*pisme<<p[2]<<endl;///
  cout<<om::greenbg<<"Final FS blob:"<<om::reset<<" \n"<<*pfs<<p[3]<<endl;/////
#endif

  return true;

}



//=============================================================================





//eof
