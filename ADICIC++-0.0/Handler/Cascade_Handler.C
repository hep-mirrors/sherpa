//bof
//Version: 2 ADICIC++-0.0/2004/09/10

//Implementation of Cascade_Handler.H.



#include <ioextra>
#include "Cascade_Handler.H"





using namespace std;
using namespace ATOOLS;
using namespace ADICIC;





//#include "....exa.cc"





//=============================================================================



//So far there is no static Cascade_Handler.
int Cascade_Handler::s_count=0;
const int& Cascade_Handler::InStore=Cascade_Handler::s_count;



//=============================================================================



Cascade_Handler::Cascade_Handler()
  : p_cas(NULL), m_chh(), v_count(), l_mit(), l_itt() {
  ++s_count;
  v_count.reserve(10);
  v_count.assign(10,0);
}





Cascade_Handler::Cascade_Handler(Cascade& cas)
  : p_cas(NULL), m_chh(), v_count(), l_mit(), l_itt() {

  ++s_count;

  v_count.reserve(10);
  v_count.assign(10,0);

  if(cas|*this) {
    if(cas.IsHandledBy(*this));
    else {
      cerr<<"\nBug: Wrong Cascade-Cascade_Handler connection emerged!\n";
      assert(cas.IsHandledBy(*this));
    }
  } else {
    cerr<<"\nMethod: "<<__PRETTY_FUNCTION__<<": "
	<<"Warning: Attaching Cascade failed!\n"<<endl;
  }

}





Cascade_Handler::~Cascade_Handler() {
  --s_count;
  v_count.clear();
  if(!p_cas) return;
  if(p_cas->IsHandledBy(*this)==false) {
    cerr<<"\nBug: Wrong Cascade-Cascade_Handler connection emerged!\n";
    assert(p_cas->IsHandledBy(*this));
  }
  FreeCascade();
  *p_cas|0;
}



//=============================================================================



void Cascade_Handler::PrintCounter() const {
  cout<<"\n"<<"======= "<<__PRETTY_FUNCTION__<<" =======\n";
  cout<<"Total number of cascade evolutions              = "
      <<v_count[Counter::total]<<"\n";
  cout<<"Total number of gluon  emissions in cascade     = "
      <<v_count[Counter::emission]<<"\t("
      <<1.0*v_count[Counter::emission]/v_count[Counter::total]<<").\n";
  cout<<"Total number of gluon  emissions -> gluons      = "
      <<v_count[Counter::gluon]<<"\t("
      <<1.0*v_count[Counter::gluon]/v_count[Counter::emission]<<").\n";
  cout<<"Total number of gluon splittings in cascade     = "
      <<v_count[Counter::splitting]<<"\t("
      <<1.0*v_count[Counter::splitting]/v_count[Counter::total]<<").\n";
  cout<<"Total number of gluon splittings -> ddbar pairs = "
      <<v_count[Counter::dquark]<<"\t("
      <<1.0*v_count[Counter::dquark]/v_count[Counter::splitting]<<").\n";
  cout<<"Total number of gluon splittings -> uubar pairs = "
      <<v_count[Counter::uquark]<<"\t("
      <<1.0*v_count[Counter::uquark]/v_count[Counter::splitting]<<").\n";
  cout<<"Total number of gluon splittings -> ssbar pairs = "
      <<v_count[Counter::squark]<<"\t("
      <<1.0*v_count[Counter::squark]/v_count[Counter::splitting]<<").\n";
  cout<<"Total number of gluon splittings -> ccbar pairs = "
      <<v_count[Counter::cquark]<<"\t("
      <<1.0*v_count[Counter::cquark]/v_count[Counter::splitting]<<").\n";
  cout<<"Total number of gluon splittings -> bbbar pairs = "
      <<v_count[Counter::bquark]<<"\t("
      <<1.0*v_count[Counter::bquark]/v_count[Counter::splitting]<<").\n";
  cout<<"Ratio of the # of splittings to the # of emissions = "
      <<1.0*v_count[Counter::splitting]/v_count[Counter::emission]<<" .\n";
  cout<<"==================================================================\n";
  assert(v_count[Counter::dquark]+v_count[Counter::uquark]+
	 v_count[Counter::squark]+v_count[Counter::cquark]+
	 v_count[Counter::bquark]==
	 v_count[Counter::splitting]);
}





const bool Cascade_Handler::EvolveCascadeByOneStep() {/////////////////////////
  return false;////////////////////////////////////////////////////////////////
}//////////////////////////////////////////////////////////////////////////////





const bool Cascade_Handler::EvolveCascade() {

  if(!p_cas) return false;
  if(p_cas->IsEmpty()) return false;

  m_typ=p_cas->CascadeType();

  if(p_cas->Status()==On && m_typ!=Cascade::incorrect);
  else return false;

  assert(!p_cas->IsEvolved());
  assert(!p_cas->HasBlockedChain());
  assert(p_cas->RootChainNumber());

  assert(!m_chh.IsDocked());
  assert(l_mit.empty() && l_itt.empty());

  //No testing of global parameters.

  m_old=p_cas->Momentum();

  {
    list<Cascade::Mirror>::iterator miter=p_cas->MirrorList().begin();
    for(list<Chain*>::iterator iter=p_cas->ChainPointerList().begin();
	iter!=p_cas->ChainPointerList().end(); ++iter) {
      if((*iter)->Status()==On) {
	l_mit.push_back(miter);
	l_itt.push_back(iter);
      }
      ++miter;
    }
  }

  do {
    if(EvolveCurrentChain());
    else {
      cerr<<"\nMethod: "<<__PRETTY_FUNCTION__<<": "
	  <<"Error: Evolution of current chain failed!\n"<<endl;
      return false;
    }
  }
  while(!l_itt.empty());

  ++v_count[Counter::total];

  //Check 4-momenta.
  Vec4D testcas;
  p_cas->CheckMomentumConservation(testcas);
  testcas+=(-1)*m_old;
  for(char i=0; i<4; ++i) assert(dabs(testcas[i])<1.0e-9);

  return true;

}



//=============================================================================



void Cascade_Handler::FreeCascade() {
  if(m_chh.IsDocked()) {
    assert(p_cas);
    for(list<Chain*>::const_iterator it=p_cas->ChainPointerList().begin();
	it!=p_cas->ChainPointerList().end(); ++it)
      **it|0;
  }
  l_mit.clear();
  l_itt.clear();
}





const bool Cascade_Handler::EvolveCurrentChain() {

  list<Cascade::Mirror>::iterator itmi=l_mit.front();
  list<Chain*>::iterator itca=l_itt.front();
  p_cuc=*itca;
  assert((*p_cuc)|m_chh);
  m_vec=p_cuc->Momentum();

#ifdef CASCADE_HANDLER_OUTPUT
  cout<<m_old<<"\t"<<m_old.Abs2()<<"\n";
  cout<<m_vec<<"\t"<<m_vec.Abs2()<<"\n";
#endif

  m_vec*=(-1.0);

  while(m_chh.EvolveChainByOneStep()) {

    assert(p_cuc->Status()==On);
    p_nec=NULL;
    m_chh.DecoupleNew(p_nec,m_kfc,f_bot);

    if(p_nec) {
      assert(p_nec->Status()==On);
      ++v_count[Counter::splitting];
      m_vec+=p_nec->Momentum();
      ++(*itmi).second;
      if(f_bot) {
	++itmi;
	itmi=p_cas->MirrorList().insert(itmi,Cascade::Mirror(p_cuc->Name,0));
	l_mit.push_back(itmi);
	--itmi;
	++itca;
	itca=p_cas->ChainPointerList().insert(itca,p_nec);
	l_itt.push_back(itca);
	--itca;
      } else {
	itmi=p_cas->MirrorList().insert(itmi,Cascade::Mirror(p_cuc->Name,0));
	l_mit.push_back(itmi);
	++itmi;
	itca=p_cas->ChainPointerList().insert(itca,p_nec);
	l_itt.push_back(itca);
	++itca;
      }
    } else {
      ++v_count[Counter::emission];
    }
 
    switch(m_kfc) {
    case kf::gluon: ++v_count[Counter::gluon]; break;
    case kf::d    : ++v_count[Counter::dquark]; break;
    case kf::u    : ++v_count[Counter::uquark]; break;
    case kf::s    : ++v_count[Counter::squark]; break;
    case kf::c    : ++v_count[Counter::cquark]; break;
    case kf::b    : ++v_count[Counter::bquark]; break;
    default       : assert(0);
    }

  }

  m_vec+=p_cuc->Momentum();
  p_cas->UpdateMomentum(1.0,m_vec);

  l_mit.remove(itmi);
  l_itt.remove(itca);
  p_cuc->SetStatus()=Off;
  (*p_cuc)|0;

#ifdef CASCADE_HANDLER_OUTPUT
  cout<<(*p_cas)<<"\n";
  Vec4D test;
  cout<<p_cas->CheckMomentumConservation(test)<<"\t"<<test<<"\n";
  cout<<p_cuc->CheckMomentumConservation(test)<<"\t"<<test<<"\n\n";
  cin>>enter; cout<<endl;
#endif

  return true;

}



//=============================================================================





//eof
