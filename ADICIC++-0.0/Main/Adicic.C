//bof
//Version: 4 ADICIC++-0.0/2005/08/02

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
  : m_parcountcorr(2), m_startnum(32767),
    m_total(0), m_fail(0), m_noem(0), p_handler(NULL), m_cascade(),
    p_blist(NULL), l_paic(),v_conn(), v_resi() {
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
  : m_parcountcorr(2), m_startnum(32767),
    m_total(0), m_fail(0), m_noem(0), p_handler(NULL), m_cascade(),
    p_blist(NULL), l_paic(), v_conn(), v_resi() {
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
  //if(m_cascade.ParticleNumber()>3) m_cascade.Print();////////////////////////
  return 1;
}





bool Adicic::ExtractPartons(ATOOLS::Blob_List* blobs) {

  assert(blobs);
  p_blist=blobs;

  if(m_total==1 || m_total%1000==0)
    cout<<"                "<<om::red<<om::blackbg<<m_total<<om::reset<<"\n";

  boolint res=m_cascade.ExtractPartons(l_paic);
  if(res.flag==false) {
    msg.Error()<<__PRETTY_FUNCTION__<<": Error: "
	       <<"Cannot extract partons from cascade!\n"<<endl;
    return false;
  }
  m_startnum=res.numr;
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
  cout<<om::brownbg<<"1st blob:"<<om::reset<<" \n"<<*blob<<endl;///////////////
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
  for(int i=0; i<pfs->NInP(); ++i)
    pfs->InParticle(i)->SetStatus(part_status::decayed);
  int num=Max(pfs->InParticle(0)->Number(),pfs->InParticle(1)->Number());
  num=num+1-m_startnum;    //-plist.front()->Number();
  for(Particle_List::iterator piter=plist.begin(); piter!=plist.end();
      ++piter) {
    assert((*piter)->Info()=='f');
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

  for(int i=0; i<pfs->NInP(); ++i)
    pfs->InParticle(i)->SetStatus(part_status::decayed);
  int num=Max(pfs->InParticle(0)->Number(),pfs->InParticle(1)->Number());
  //cout<<num<<":"<<plist.front()->Number()<<":"<<m_startnum<<":";
  num=num+3-m_startnum;    //cout<<num<<endl;

  for(Particle_List::iterator piter=plist.begin(); piter!=plist.end();
      ++piter) {
    switch((*piter)->Info()) {
    case 'f':
      //cout<<(*piter)->Number()<<endl;
      pis->AddToOutParticles(*piter);
      (*piter)->SetInfo('F');
      (*piter)->SetNumber(-((*piter)->Number()+num-m_parcountcorr));
      break;
    case 'i':
      //cout<<(*piter)->Number()<<endl;
      pis->AddToInParticles(*piter);
      (*piter)->SetInfo('I');
      (*piter)->SetNumber(-((*piter)->Number()+num-m_parcountcorr));
      break;
    case 'c':
      v_conn.push_back(*piter);
      //cout<<(**piter)<<endl;
      break;
    case 'r':
      v_resi.push_back(*piter);
      //cout<<(**piter)<<endl;
      break;
    default :
      cerr<<"\nMethod: "<<__PRETTY_FUNCTION__<<": "
	  <<"Warning: Particle carries non-expected information!\n"<<endl;
    }
  }
  assert(v_conn.empty() ^ v_resi.empty());
  if(pis->InParticle(0)->Momentum()[3]<pis->InParticle(1)->Momentum()[3])
    pis->SwapInParticles(0,1);

  Particle* par=new Particle(*pisme->OutParticle(0));
  Particle* per=new Particle(*pisme->OutParticle(1));
  pis->AddToOutParticles(par); pis->AddToOutParticles(per);
  pisme->AddToInParticles(par); pisme->AddToInParticles(per);
  const Vec4D  PZ=par->Momentum()+per->Momentum();
  const Vec4D& PZprime=-1*m_cascade.Momentum();
  assert(dabs(PZ.Abs2()-PZprime.Abs2())<1.0e-7);
  Poincare fly(PZ); Poincare flyprime(PZprime);
  Vec4D addup, eddup;
  addup=par->Momentum(); fly.Boost(addup); flyprime.BoostBack(addup);
  par->SetMomentum(addup);
  addup=per->Momentum(); fly.Boost(addup); flyprime.BoostBack(addup);
  per->SetMomentum(addup);

  addup=eddup=Vec4D(0.0,0.0,0.0,0.0);
  if(v_resi.empty()) {    //The correlated-particle treatment.
    assert(pfs->NInP()==int(v_conn.size()));
    const size_t num=v_conn.size();
    if(pfs->InParticle(0)->RefFlav()!=v_conn[0]->RefFlav()) {
      assert(pfs->InParticle(0)->RefFlav()==v_conn[num-1]->RefFlav());
      assert(pfs->InParticle(num-1)->RefFlav()==v_conn[0]->RefFlav());
      assert(pfs->InParticle(0)->Number()==v_conn[num-1]->Number());
      assert(pfs->InParticle(num-1)->Number()==v_conn[0]->Number());
      for(size_t i=0; i<num; ++i) {
	v_conn[num-1-i]->SetInfo(pfs->InParticle(i)->Info());
	v_conn[num-1-i]->SetStatus(part_status::active);
	pfs->AddToOutParticles(v_conn[num-1-i]);
	addup+=pfs->InParticle(i)->Momentum();
	eddup+=v_conn[num-1-i]->Momentum();
      }
    } else {
      assert(pfs->InParticle(num-1)->RefFlav()==v_conn[num-1]->RefFlav());
      assert(pfs->InParticle(0)->Number()==v_conn[0]->Number());
      assert(pfs->InParticle(num-1)->Number()==v_conn[num-1]->Number());
      for(size_t i=0; i<num; ++i) {
	v_conn[i]->SetInfo(pfs->InParticle(i)->Info());
	v_conn[i]->SetStatus(part_status::active);
	pfs->AddToOutParticles(v_conn[i]);
	addup+=pfs->InParticle(i)->Momentum();
	eddup+=v_conn[i]->Momentum();
      }
    }
    v_conn.clear();
    /*
    for(int i=0; i<pfs->NInP(); ++i) {    //For the purpose of comparison.
      Particle help(*pfs->InParticle(i));
      Vec4D temp=help.Momentum();
      fly.Boost(temp); flyprime.BoostBack(temp);
      help.SetMomentum(temp);
      help.SetStatus(part_status::active);
      cout<<help<<endl;
    }
    */
  }
  else {    //The residual-particle treatment.
    assert(v_resi.size()==1);
    assert(PZprime==v_resi[0]->Momentum());
    delete v_resi[0];
    v_resi.clear();
    for(int i=0; i<pfs->NInP(); ++i) {
      Particle* help=new Particle(*pfs->InParticle(i));
      pfs->AddToOutParticles(help);
      help->SetStatus(part_status::active);
      Vec4D temp=help->Momentum();
      addup+=temp;
      fly.Boost(temp); flyprime.BoostBack(temp);
      eddup+=temp;
      help->SetMomentum(temp);
      help->SetStatus(part_status::active);
      //cout<<*help<<endl;
    }
  }

  //addup-=PZ; eddup-=PZprime;
  //for(char i=0; i<4; ++i)
  //  assert(dabs(addup[i])<1.0e-10 && dabs(eddup[i])<1.0e-10);
  assert(PZ==addup);
  assert(PZprime==eddup);

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
