//bof
//Version: 4 ADICIC++-0.0/2006/07/02

//Implementation of Chain_Handler.H.



#include <stdlib.h>
#include <iomanip>
#include <ioextra>
#include "Poincare.H"
#include "Chain_Handler.H"
#include "Dipole_Parameter.H"





using namespace std;
using namespace ATOOLS;
using namespace ADICIC;





//#include "....exa.cc"





//=============================================================================



//ostream& ADICIC::operator<<(ostream& ost, const ADICIC::...&) {
//}



//=============================================================================



//So far there is no static Chain_Handler.
int Chain_Handler::s_count=0;
const int& Chain_Handler::InStore=Chain_Handler::s_count;



//=============================================================================



template
const bool Chain_Handler::FindTheDipole<Chain_Evolution_Strategy::Unknown>();
//template const bool
//Chain_Handler::FindTheDipole<Chain_Evolution_Strategy::Production>();
//template const bool
//Chain_Handler::FindTheDipole<Chain_Evolution_Strategy::Emission>();
template
const bool Chain_Handler::FindTheDipole<Chain_Evolution_Strategy::Mass>();



//=============================================================================



Chain_Handler::Chain_Handler()
  : f_below(false), m_code(), p_cix(NULL), m_box(),
    m_nii(-7), m_k2tcomp(0.0), p_cha(NULL),
    m_dh1(), m_dh2(),
    p_dhwait(&m_dh1), p_dhaciv(&m_dh2), p_dhtemp(NULL),
    i_fix(), i_run() {

  ++s_count;

  PresetCompScale();

#ifdef __GNUC__
#if __GNUC__ >2
  for(size_t i=0; ; ++i) {
    if(Chain_Evolution_Strategy::List[i]==Chain_Evolution_Strategy::stop)
      break;
    assert(Chain_Evolution_Strategy::List[i]>=0);
    if(v_finddip.size()<size_t(Chain_Evolution_Strategy::List[i])+1)
      v_finddip.resize(size_t(Chain_Evolution_Strategy::List[i])+1,NULL);
  }
  this->InitStrategies();
#endif
#endif

}





Chain_Handler::Chain_Handler(Chain& cha)
  : f_below(false), m_code(), p_cix(NULL), m_box(),
    m_nii(-7), m_k2tcomp(0.0), p_cha(NULL),
    m_dh1(), m_dh2(),
    p_dhwait(&m_dh1), p_dhaciv(&m_dh2), p_dhtemp(NULL),
    i_fix(), i_run() {

  ++s_count;

  PresetCompScale();

#ifdef __GNUC__
#if __GNUC__ >2
  for(size_t i=0; ; ++i) {
    if(Chain_Evolution_Strategy::List[i]==Chain_Evolution_Strategy::stop)
      break;
    assert(Chain_Evolution_Strategy::List[i]>=0);
    if(v_finddip.size()<size_t(Chain_Evolution_Strategy::List[i])+1)
      v_finddip.resize(size_t(Chain_Evolution_Strategy::List[i])+1,NULL);
  }
  this->InitStrategies();
#endif
#endif

  if(cha|*this) {
    if(cha.IsHandledBy(*this));
    else {
      cerr<<"\nBug: Wrong Chain-Chain_Handler connection emerged!\n";
      assert(cha.IsHandledBy(*this));
    }
  } else {
    cerr<<"\nMethod: "<<__PRETTY_FUNCTION__<<": "
	<<"Warning: Attaching Chain failed!\n"<<endl;
  }

}





Chain_Handler::~Chain_Handler() {

  --s_count;

  if(p_cix) delete p_cix;

  if(!p_cha) return;
  if(p_cha->IsHandledBy(*this)==false) {
    cerr<<"\nBug: Wrong Chain-Chain_Handler connection emerged!\n";
    assert(p_cha->IsHandledBy(*this));
  }

  FreeChain();
  *p_cha|0;

}



//=============================================================================



void Chain_Handler::DecoupleNew(Carrier& car) {
  assert(car.pCha==NULL);
  assert(car.Mup.GetItsVec().empty());
  car.ChaOrder=f_below; f_below=false;
  car.EmitFlav=m_code; m_code=ATOOLS::Flavour();
  car.pCha=p_cix; p_cix=NULL;
  if(m_box.Mup.GetItsVec().empty());
  else car.Mup.Swap(m_box.Mup);
}


void Chain_Handler::RemoveNewProducts() {    //Resets the news.
  f_below=false;
  m_code=ATOOLS::Flavour();
  if(p_cix) { delete p_cix; p_cix=NULL;}
  //m_box.Mup isn't really new, it was a condition, so here it's not changed.
  //It is handed over (to the cascade) in case other chains also have to be
  //corrected by applying the same transformations saved in the m_box.Mup.
}





const bool Chain_Handler::EvolveChainByOneStep() {

  if(!p_cha) return false;
  if(p_cha->IsEmpty()) return false;

  m_typ=p_cha->ChainType();
  this->PresetCompScale();

  if(p_cha->Status()==On &&
     m_typ!=Chain::incorrect &&
     p_cha->LastScale()>m_k2tcomp &&
     p_cha->ParticleNumber()<dpa.evo.ChainParticleLimit());//!!!!!!!!!!!!!!!!!!
  else return false;

  this->CleanUpDHs();
  this->RemoveNewProducts();    //No testing of global parameters.

  m_box.Reset();    //Note!
  //All the preceding results must have found their way into a carrying list or
  //been deleted before, otherwise this would be a potential memory leak.

  if(this->FindDipole()) return this->ModifyChain();

  return false;

}





const bool Chain_Handler::EvolveChain() {//////////////////////////////////////
  return false;////////////////////////////////////////////////////////////////
}//////////////////////////////////////////////////////////////////////////////





const bool Chain_Handler::CorrectChain(const Multipoincare& lmup) {
  if(!p_cha) return false;
  if(!p_cha->IsLine()) return false;
  if(lmup.Mode()!=rdt::iirecoil) return false;
  for(list<Particle*>::iterator it=p_cha->CorrParticlePointerList().begin();
      it!=p_cha->CorrParticlePointerList().end(); ++it) {
    Vec4D temp((*it)->Momentum());
    assert(lmup.Apply(temp));
    (*it)->SetMomentum(temp);
  }
  list<Dipole_Particle*> fins;
  assert(p_cha->BranchPointer()->Incoming()==false);
  fins.push_back(p_cha->BranchPointer());
  assert(p_cha->AntibranchPointer()->Incoming()==false);
  fins.push_back(p_cha->AntibranchPointer());
  for(list<Dipole::Glubranch*>::const_iterator
	it=p_cha->GlubranchPointerList().begin();
      it!=p_cha->GlubranchPointerList().end(); ++it) {
    assert((*it)->Incoming()==false);
    fins.push_back(*it);
  }
  for(list<Dipole_Particle*>::iterator it=fins.begin(); it!=fins.end(); ++it) {
    Vec4D temp=(*it)->Momentum();
    assert(lmup.Apply(temp));
    (*it)->SetMomentum(temp);    //Updates dipole mass!
    //But the chain momentum is still not up to date, I would guess!!!!!!!!!!!!
  }
  return false;//so long to see that it does nothing///////////////////////////
}



//=============================================================================



void Chain_Handler::InitStrategies() {
  //This must correspond to the settings made in Evolution_Strategy.hpp.
  assert(v_finddip.size()>Chain_Evolution_Strategy::Unknown);
  v_finddip[Chain_Evolution_Strategy::Unknown]=&ADICIC::Chain_Handler::
    FindTheDipole<Chain_Evolution_Strategy::Unknown>;
  assert(v_finddip.size()>Chain_Evolution_Strategy::Production);
  v_finddip[Chain_Evolution_Strategy::Production]=&ADICIC::Chain_Handler::
    FindTheDipole<Chain_Evolution_Strategy::Production>;
  assert(v_finddip.size()>Chain_Evolution_Strategy::Emission);
  v_finddip[Chain_Evolution_Strategy::Emission]=&ADICIC::Chain_Handler::
    FindTheDipole<Chain_Evolution_Strategy::Emission>;
  assert(v_finddip.size()>Chain_Evolution_Strategy::Mass);
  v_finddip[Chain_Evolution_Strategy::Mass]=&ADICIC::Chain_Handler::
    FindTheDipole<Chain_Evolution_Strategy::Mass>;
  //for(size_t i=0; i<v_finddip.size(); ++i)
  //  cout<<setw(7)<<i<<" : "<<v_finddip[i]<<endl;
  //abort();
}





void Chain_Handler::CleanUpDHs() {
  assert(!m_dh1.IsDocked() && !m_dh2.IsDocked());
  m_dh1.RemoveNewProducts();
  m_dh2.RemoveNewProducts();
  p_dhwait=&m_dh1;
  p_dhaciv=&m_dh2;
  p_dhtemp=NULL;
  i_fix=list<Dipole*>::iterator();
  i_run=list<Dipole*>::iterator();
}





void Chain_Handler::FreeChain() {
  if(m_dh1.IsDocked() || m_dh2.IsDocked()) {
    assert(p_cha);
    for(list<Dipole*>::const_iterator it=p_cha->DipolePointerList().begin();
	it!=p_cha->DipolePointerList().end(); ++it)
      **it|0;
  }
}



//=============================================================================



template<Chain_Evolution_Strategy::Type _Strategy>
const bool Chain_Handler::FindTheDipole() {
  cerr<<"\nMethod: "<<__PRETTY_FUNCTION__<<": "
      <<"Warning: Method has not been implemented yet!\n"<<endl;
  return false;
}





template<> const bool
Chain_Handler::FindTheDipole<Chain_Evolution_Strategy::Production>() {

  static const bool isfalse=ProductionStrategyInfo();    //Gives warning.

  i_fix=p_cha->DipolePointerList().end();
  i_run=p_cha->DipolePointerList().begin();
  bool spico=p_cha->ParticleNumber()<dpa.evo.ChainCorrelationLimit();

  for(; i_run!=p_cha->DipolePointerList().end(); ++i_run) {

    Dipole& dip=**i_run;
    dip.SetFactScale()=p_cha->FactScale();    //Updates it to the current muF.
    if(dip.SpinCorr() && !spico) dip.SetSpinCorr()=false;
    assert(dip|*p_dhaciv);

    if(p_dhaciv->InduceDipoleRadiation()==isfalse) {    //Removes the warning.
      dip|0; continue;
    }

    if(dip.EmitScale()>=p_cha->LastScale()) {
      bool result;
      do {
	dip.SetBootScale()=dip.EmitScale();
	dip.SetFactScale()=p_cha->FactScale();
	result=p_dhaciv->InduceDipoleRadiation();
	if(result==false) break;
      }
      while(dip.EmitScale()>=p_cha->LastScale());
      if(result==false) {
	dip.SetBootScale()=dip.ProdScale();
	dip|0; continue;
      }
    }

    if(dip.EmitScale()<m_k2tcomp) {
      dip.SetBootScale()=dip.ProdScale();
      //Correct muF is set at the beginning of the next emission round!
      dip|0; continue;
    }

    m_k2tcomp=dip.EmitScale();
    p_dhtemp=p_dhwait;
    p_dhwait=p_dhaciv;
    p_dhaciv=p_dhtemp;
    if(p_dhtemp->IsDocked()) {
      Dipole& loser=**i_fix;
      loser.SetBootScale()=loser.ProdScale();
      //Correct muF is set at the beginning of the next emission round!
      loser|0;
    }
    i_fix=i_run;

  }

  return bool(i_fix!=p_cha->DipolePointerList().end());

}





template<> const bool
Chain_Handler::FindTheDipole<Chain_Evolution_Strategy::Emission>() {

  static const bool isfalse=EmissionStrategyInfo();    //Gives warning.

  i_fix=p_cha->DipolePointerList().end();
  i_run=p_cha->DipolePointerList().begin();
  bool spico=p_cha->ParticleNumber()<dpa.evo.ChainCorrelationLimit();

  for(; i_run!=p_cha->DipolePointerList().end(); ++i_run) {

    Dipole& dip=**i_run;
    dip.SetFactScale()=p_cha->FactScale();    //Updates it to the current muF.
    if(dip.SpinCorr() && !spico) dip.SetSpinCorr()=false;
    assert(dip|*p_dhaciv);

    if(p_dhaciv->InduceDipoleRadiation()==isfalse) {    //Removes the warning.
      dip|0; continue;
    }

    assert(dip.EmitScale()<p_cha->LastScale());
    if(dip.EmitScale()<m_k2tcomp) { dip|0; continue;}

    m_k2tcomp=dip.EmitScale();
    p_dhtemp=p_dhwait; p_dhwait=p_dhaciv; p_dhaciv=p_dhtemp;
    if(p_dhtemp->IsDocked()) { Dipole& loser=**i_fix; loser|0;}
    i_fix=i_run;

  }

  if(i_fix==p_cha->DipolePointerList().end()) m_k2tcomp=0.0;
  i_run=p_cha->DipolePointerList().begin();
  for(; i_run!=p_cha->DipolePointerList().end(); ++i_run) {
    Dipole& dip=**i_run;
    if(i_run!=i_fix) dip.SetBootScale()=m_k2tcomp;
  }

  return bool(i_fix!=p_cha->DipolePointerList().end());

}



//=============================================================================



const bool Chain_Handler::ModifyChain() {

  p_win=*i_fix;    //cout<<*p_win<<endl;
  m_nii=-7;
  m_dm2=p_win->InvMass();
  m_emi=Vec4D();

  if(p_win->GetTopBranchPointer()->Incoming()) {
    m_nii=0;
    m_vec=-1*p_win->GetTopBranchPointer()->Momentum();
  } else {
    m_vec=p_win->GetTopBranchPointer()->Momentum();
  }
  if(p_win->GetBotBranchPointer()->Incoming()) {
    m_vec-=p_win->GetBotBranchPointer()->Momentum();
  } else {
    m_nii=-7;
    m_vec+=p_win->GetBotBranchPointer()->Momentum();
  }
  m_old=p_cha->Momentum();
  p_cha->UpdateMomentum(-1.0,m_vec);

#ifdef CHAIN_HANDLER_OUTPUT
  cout<<m_old<<"\t"<<m_old.Abs2()<<"\n";
  cout<<m_vec<<"\t"<<m_vec.Abs2()<<"\n";
  cout<<p_cha->Momentum()<<"\t"
      <<p_cha->Momentum().Abs2()<<"\t"
      <<p_cha->InvMass()<<endl;
#endif

  assert(p_dhwait->FinishDipoleRadiation());

  p_dhwait->DecoupleNew(m_box);
  if(m_box.pDip) {
    if(m_box.pQua) {
      assert(m_box.pAqu && m_box.pGlu);
      if(m_box.DipOrder) m_code=m_box.pAqu->Flav();
      else               m_code=m_box.pQua->Flav();
      (*p_win)|0;
      EmitQuark();
    } else {
      assert(!m_box.pAqu && m_box.pGlu);
      m_code=Flavour(kf::gluon);
      (*p_win)|0;
      EmitGluon();
    }
  } else {
    assert(m_box.pAqu && m_box.pQua && m_box.pGlu);
    m_code=m_box.pQua->Flav();
    assert(m_box.pAqu->Flav().Bar()==m_code);
    (*p_win)|0;
    if(m_typ==Chain::line) {
      if(m_box.DipOrder){
#ifdef CHAIN_HANDLER_OUTPUT
	cout<<p_win->GetTopBranchPointer()->Momentum()<<"\n";
	cout<<m_box.pAqu->Momentum()<<"\n";
	cout<<m_box.pQua->Momentum()<<"\n";
	cout<<p_win->GetTopBranchPointer()->Momentum()+
	  m_box.pAqu->Momentum()+m_box.pQua->Momentum()<<endl;
#endif
	BotSplitLine();
      }
      else {
#ifdef CHAIN_HANDLER_OUTPUT
	cout<<m_box.pAqu->Momentum()<<"\n";
	cout<<m_box.pQua->Momentum()<<"\n";
	cout<<p_win->GetBotBranchPointer()->Momentum()<<"\n";
	cout<<p_win->GetBotBranchPointer()->Momentum()+
	  m_box.pAqu->Momentum()+m_box.pQua->Momentum()<<endl;
#endif
	TopSplitLine();
      }
    } else {
      if(m_box.DipOrder) BotSplitRing();
      else               TopSplitRing();
    }
  }

  //Calculation of all factorization scale options:
  p_cha->FactScaleBox()[fascat::p2t]=p_win->EmitScale();
  p_cha->FactScaleBox()[fascat::m2t]=
    p_cha->FactScaleBox()[fascat::k2t]=m_emi.PPerp2();
  p_cha->FactScaleBox()[fascat::m2t]+=m_dm2;
  p_cha->FactScaleBox()[fascat::shat]=m_dm2*
    p_cha->FactScaleBox()[fascat::p2t]/p_cha->FactScaleBox()[fascat::k2t];

#ifdef CHAIN_HANDLER_OUTPUT
  cout<<p_cha->FactScale()<<" : "<<p_win->FactScale()<<" | "
      <<p_cha->FactScaleBox()[fascat::p2t]<<" : "
      <<p_cha->FactScaleBox()[fascat::k2t]<<"("
      <<sqrt(p_cha->FactScaleBox()[fascat::k2t])<<") : "
      <<p_cha->FactScaleBox()[fascat::m2t]<<"("
      <<sqrt(p_cha->FactScaleBox()[fascat::m2t])<<") : "
      <<p_cha->FactScaleBox()[fascat::shat]<<"("
      <<sqrt(p_cha->FactScaleBox()[fascat::shat])<<")"<<endl;
#endif

  return true;

}





void Chain_Handler::DistributeIIRecoil() {
  //Must be executed before new particles are introduced into the chain.

  //cout<<"Number of Poincares kept: "<<m_box.Mup.GetItsVec().size()
  //    <<"   using mode: "<<m_box.Mup.Mode()<<endl;

  assert(m_nii==0);
  assert(!p_cha->CorrParticlePointerList().empty());
  assert(m_box.Mup.GetItsVec().size());
  for(list<Particle*>::iterator it=p_cha->CorrParticlePointerList().begin();
      it!=p_cha->CorrParticlePointerList().end(); ++it) {
    Vec4D temp((*it)->Momentum());
    assert(m_box.Mup.Apply(temp));
    (*it)->SetMomentum(temp);
  }

  list<Dipole_Particle*> fins;
  if(p_cha->BranchPointer()->Incoming()==false) {
    ++m_nii; fins.push_back(p_cha->BranchPointer());}
  if(p_cha->AntibranchPointer()->Incoming()==false) {
    ++m_nii; fins.push_back(p_cha->AntibranchPointer());}
  for(list<Dipole::Glubranch*>::const_iterator
	it=p_cha->GlubranchPointerList().begin();
      it!=p_cha->GlubranchPointerList().end(); ++it)
    if((*it)->Incoming()==false) { ++m_nii; fins.push_back(*it);}
  if(m_nii==0) return;
  for(list<Dipole_Particle*>::iterator it=fins.begin(); it!=fins.end(); ++it) {
    Vec4D temp=(*it)->Momentum();
    m_vec-=temp;
    assert(m_box.Mup.Apply(temp));
    (*it)->SetMomentum(temp);    //Updates dipole mass!
    m_vec+=temp;
  }
}





void Chain_Handler::EmitGluon() {

  m_emi=m_vec=m_box.pGlu->Momentum();    //It is emitted, so it is an F gluon.

  if(m_box.DipOrder) {
    i_run=i_fix;
    ++i_run;
    i_run=p_cha->DipolePointerList().insert(i_run,m_box.pDip);
    //if(p_win->GetTopBranchPointer()->Incoming())
    //if(m_box.pDip->GetBotBranchPointer()->Incoming())
    if(m_nii>=0) {
      m_vec-=p_win->GetTopBranchPointer()->Momentum();
      m_vec-=m_box.pDip->GetBotBranchPointer()->Momentum();
      DistributeIIRecoil();    //-1*m_vec corresponds to new i+i-g!
    } else {
      m_vec+=p_win->GetTopBranchPointer()->Momentum();
      m_vec+=m_box.pDip->GetBotBranchPointer()->Momentum();
    }
  } else {
    i_run=p_cha->DipolePointerList().insert(i_fix,m_box.pDip);
    ++i_run;
    i_fix=i_run;
    //if(m_box.pDip->GetTopBranchPointer()->Incoming())
    //if(p_win->GetBotBranchPointer()->Incoming())
    if(m_nii>=0) {
      m_vec-=m_box.pDip->GetTopBranchPointer()->Momentum();
      m_vec-=p_win->GetBotBranchPointer()->Momentum();
      DistributeIIRecoil();    //-1*m_vec corresponds to new i+i-g!
    } else {
      m_vec+=m_box.pDip->GetTopBranchPointer()->Momentum();
      m_vec+=p_win->GetBotBranchPointer()->Momentum();
    }
  }

  p_cha->GlubranchPointerList().push_back(m_box.pGlu);
  p_cha->UpdateMomentum(1.0,m_vec);
  p_cha->SetLastScale()=m_k2tcomp;

  //Check 4-momenta.
  Vec4D testcha;
  assert(p_cha->CheckMomentumConservation(testcha));
  if(p_cha->INumber()==0) {
    testcha-=m_old;
    for(char i=0; i<4; ++i) assert(dabs(testcha[i])<1.0e-9);
  } else {
    assert(dabs(testcha.Abs2()-m_old.Abs2())<1.0e-7);
  }

}





void Chain_Handler::EmitQuark() {

  assert(m_nii>=0);
  m_vec=-1*m_box.pGlu->Momentum();    //This is the new initial(=I) gluon.
  //p_cha->GlubranchPointerList().push_back(m_box.pGlu);

  if(m_box.DipOrder) {
    i_run=i_fix;
    ++i_run;
    i_run=p_cha->DipolePointerList().insert(i_run,m_box.pDip);
    assert(p_win->GetTopBranchPointer()->Incoming());
    m_vec-=p_win->GetTopBranchPointer()->Momentum();
    assert(m_box.pDip->GetBotBranchPointer()->Incoming()==false);
    m_vec+=m_box.pDip->GetBotBranchPointer()->Momentum();
    DistributeIIRecoil();
    p_cha->AntibranchPointer()=m_box.pAqu;    //The chain's new Antibranch.
    delete m_box.pQua; m_box.pQua=NULL;    //Remove the old one.
    m_emi=m_box.pAqu->Momentum();
  } else {
    i_run=p_cha->DipolePointerList().insert(i_fix,m_box.pDip);
    ++i_run;
    i_fix=i_run;
    assert(m_box.pDip->GetTopBranchPointer()->Incoming()==false);
    m_vec+=m_box.pDip->GetTopBranchPointer()->Momentum();
    assert(p_win->GetBotBranchPointer()->Incoming());
    m_vec-=p_win->GetBotBranchPointer()->Momentum();
    DistributeIIRecoil();
    p_cha->BranchPointer()=m_box.pQua;    //The new Branch of the chain.
    delete m_box.pAqu; m_box.pAqu=NULL;    //Remove the old one.
    m_emi=m_box.pQua->Momentum();
  }

  p_cha->GlubranchPointerList().push_back(m_box.pGlu);
  p_cha->UpdateMomentum(1.0,m_vec);
  p_cha->SetLastScale()=m_k2tcomp;

  //Check 4-momenta.
  Vec4D testcha;
  assert(p_cha->CheckMomentumConservation(testcha));
  assert(dabs(testcha.Abs2()-m_old.Abs2())<1.0e-7);

}





void Chain_Handler::TopSplitLine() {

  //There are check asserts to find !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef CHAIN_HANDLER_OUTPUT
  cout<<__PRETTY_FUNCTION__;
  p_cha->Print();
#endif

  f_below=true;    //New Chain emerges below the old one.

  //m_box.pAqu is the new antiquark.
  //m_box.pQua is the new quark.
  //m_box.pGlu is the split gluon.

  p_cix=new Chain(); assert(p_cix);
  //p_cix->SetSource()=p_cha->Name;
  p_cix->m_memo=p_cha->m_name;
  p_cix->BranchPointer()=m_box.pQua;
  p_cix->AntibranchPointer()=p_cha->AntibranchPointer();
  p_cha->AntibranchPointer()=NULL;

  i_run=p_cha->DipolePointerList().end(); --i_run;

#ifdef CHAIN_HANDLER_OUTPUT
  cout<<"Type: "<<(*i_run)->IsType()<<"\n";
#endif

  if(i_run==i_fix) {    //Special case - only the winner dipole.
    assert(p_win->IsType()==Dipole::qqbar);//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    m_vec=m_box.pQua->Momentum();
    m_vec+=p_cix->AntibranchPointer()->Momentum();
    p_cix->UpdateMomentum(1.0,m_vec);
    p_cha->UpdateMomentum(1.0,m_box.pAqu->Momentum());
    p_cix->DipolePointerList().push_front(p_win);
    p_cha->DipolePointerList().pop_back();
    i_run=p_cha->DipolePointerList().end(); --i_run;
    assert(*i_fix==p_win);//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    p_cix->RootPointer()=p_win;
    if(p_cha->RootPointer()==p_win) p_cha->RootPointer()=*i_run;//?????????????
  } else {    //Last dipole.
    short root=0;
    Dipole_Particle* to;
    Dipole::Glubranch* topg;
    assert((*i_run)->IsType()==Dipole::gqbar);//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if(*i_run==p_cha->RootPointer()) ++root;
    to=(*i_run)->GetTopBranchPointer().operator->();
    assert(to->OrgType()==Nil);//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    topg=static_cast<Dipole::Glubranch*>(to);
    m_vec=topg->Momentum();
    p_cix->GlubranchPointerList().push_back(topg);
    p_cha->GlubranchPointerList().remove(topg);
    p_cix->DipolePointerList().push_back(*i_run);
    --i_run;
    while(i_run!=i_fix) {    //Inbetween dipoles.
      assert((*i_run)->IsType()==Dipole::gg);//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if(*i_run==p_cha->RootPointer()) ++root;
      to=(*i_run)->GetTopBranchPointer().operator->();
      assert(to->OrgType()==Nil);//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      topg=static_cast<Dipole::Glubranch*>(to);
      m_vec+=topg->Momentum();
      p_cix->GlubranchPointerList().push_front(topg);
      p_cha->GlubranchPointerList().remove(topg);
      p_cix->DipolePointerList().push_front(*i_run);
      --i_run;
    }    //Winner dipole.
    assert(p_win->IsType()==Dipole::qg);//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if(p_win==p_cha->RootPointer()) ++root;
    p_cix->UpdateMomentum(1.0,m_box.pQua->Momentum()+m_vec+
			  p_cix->AntibranchPointer()->Momentum());
    m_vec+=p_cix->AntibranchPointer()->Momentum();
    m_vec*=(-1);
    m_vec+=p_win->GetBotBranchPointer()->Momentum();
    m_vec+=m_box.pAqu->Momentum();
    p_cha->UpdateMomentum(1.0,m_vec);
    p_cix->DipolePointerList().push_front(p_win);
    i_run=p_cha->DipolePointerList().end();
    p_cha->DipolePointerList().erase(i_fix,i_run);
    i_run=p_cha->DipolePointerList().end(); --i_run;
    i_fix=p_cix->DipolePointerList().begin();
    assert(*i_fix==p_win);//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if(!root) p_cix->RootPointer()=p_win;
    else {
      assert(root==1);//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      p_cix->RootPointer()=p_cha->RootPointer();
      p_cha->RootPointer()=*i_run;//???????????????????????????????????????????
    }
  }

  (*i_run)->RenewBranch(*m_box.pAqu);
  p_cha->AntibranchPointer()=m_box.pAqu;
  p_cha->GlubranchPointerList().remove(m_box.pGlu);
  delete m_box.pGlu; m_box.pGlu=NULL;

  //Check 4-momenta.
  Vec4D testcha, testcix, testall;
#ifdef CHAIN_HANDLER_OUTPUT
  cout<<p_cha->CheckMomentumConservation(testcha)<<"\t"<<testcha<<"\n";
  cout<<p_cix->CheckMomentumConservation(testcix)<<"\t"<<testcix<<"\n";
#else
  p_cha->CheckMomentumConservation(testcha);
  p_cix->CheckMomentumConservation(testcix);
#endif
  testall=testcha+testcix+(-1)*m_old;
#ifdef CHAIN_HANDLER_OUTPUT
  cout<<"\t"<<testall<<"\n";
  cout<<"\t"<<p_cha->Momentum()+p_cix->Momentum()+(-1)*m_old<<"\n";
#endif
  testcha+=(-1)*p_cha->Momentum(); testcix+=(-1)*p_cix->Momentum();
#ifdef CHAIN_HANDLER_OUTPUT
  cout<<testcha<<"\n";
  cout<<testcix<<"\n";
#endif
  for(char i=0; i<4; ++i) assert(dabs(testcha[i])<1.0e-9 &&
				 dabs(testcix[i])<1.0e-9 &&
				 dabs(testall[i])<1.0e-9);

  p_cha->SetLastScale()=m_k2tcomp;
  p_cix->SetLastScale()=m_k2tcomp;
  p_cix->SetStatus()=On;

#ifdef CHAIN_HANDLER_OUTPUT
  p_cha->Print();
  p_cix->Print();
#endif

  //cout<<endl; cin>>enter; cout<<endl;////////////////////////////////////////

}





void Chain_Handler::BotSplitLine() {

  //There are check asserts to find!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef CHAIN_HANDLER_OUTPUT
  cout<<__PRETTY_FUNCTION__;
  p_cha->Print();
#endif

  f_below=false;    //New Chain is on top of the old one.

  //m_box.pAqu is the new antiquark.
  //m_box.pQua is the new quark.
  //m_box.pGlu is the split gluon.

  p_cix=new Chain(); assert(p_cix);
  //p_cix->SetSource()=p_cha->Name;
  p_cix->m_memo=p_cha->m_name;
  p_cix->BranchPointer()=p_cha->BranchPointer();
  p_cha->BranchPointer()=NULL;
  p_cix->AntibranchPointer()=m_box.pAqu;

  i_run=p_cha->DipolePointerList().begin();

#ifdef CHAIN_HANDLER_OUTPUT
  cout<<"Type: "<<(*i_run)->IsType()<<"\n";
#endif

  if(i_run==i_fix) {    //Special case - only the winner dipole.
    assert(p_win->IsType()==Dipole::qqbar);//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    m_vec=p_cix->BranchPointer()->Momentum();
    m_vec+=m_box.pAqu->Momentum();
    p_cix->UpdateMomentum(1.0,m_vec);
    p_cha->UpdateMomentum(1.0,m_box.pQua->Momentum());
    p_cix->DipolePointerList().push_front(p_win);
    p_cha->DipolePointerList().pop_front();
    i_run=p_cha->DipolePointerList().begin();
    assert(*i_fix==p_win);//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    p_cix->RootPointer()=p_win;
    if(p_cha->RootPointer()==p_win) p_cha->RootPointer()=*i_run;//?????????????
  } else {    //First dipole.
    short root=0;
    Dipole_Particle* bo;
    Dipole::Glubranch* botg;
    assert((*i_run)->IsType()==Dipole::qg);//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if(*i_run==p_cha->RootPointer()) ++root;
    bo=(*i_run)->GetBotBranchPointer().operator->();
    assert(bo->OrgType()==Nil);//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    botg=static_cast<Dipole::Glubranch*>(bo);
    m_vec=botg->Momentum();
    p_cix->GlubranchPointerList().push_front(botg);
    p_cha->GlubranchPointerList().remove(botg);
    p_cix->DipolePointerList().push_front(*i_run);
    ++i_run;
    while(i_run!=i_fix) {    //Inbetween dipoles.
      assert((*i_run)->IsType()==Dipole::gg);//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if(*i_run==p_cha->RootPointer()) ++root;
      bo=(*i_run)->GetBotBranchPointer().operator->();
      assert(bo->OrgType()==Nil);//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      botg=static_cast<Dipole::Glubranch*>(bo);
      m_vec+=botg->Momentum();
      p_cix->GlubranchPointerList().push_back(botg);
      p_cha->GlubranchPointerList().remove(botg);
      p_cix->DipolePointerList().push_back(*i_run);
      ++i_run;
    }    //Winner dipole.
    assert(p_win->IsType()==Dipole::gqbar);//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if(p_win==p_cha->RootPointer()) ++root;
    p_cix->UpdateMomentum(1.0,p_cix->BranchPointer()->Momentum()+
			  m_vec+m_box.pAqu->Momentum());
    m_vec+=p_cix->BranchPointer()->Momentum();
    m_vec*=(-1);
    m_vec+=p_win->GetTopBranchPointer()->Momentum();
    m_vec+=m_box.pQua->Momentum();
    p_cha->UpdateMomentum(1.0,m_vec);
    p_cix->DipolePointerList().push_back(p_win);
    i_run=p_cha->DipolePointerList().begin();
    i_run=p_cha->DipolePointerList().erase(i_run,++i_fix);
    i_fix=p_cix->DipolePointerList().end(); --i_fix;
    assert(i_run==p_cha->DipolePointerList().begin());//!!!!!!!!!!!!!!!!!!!!!!!
    assert(*i_fix==p_win);//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if(!root) p_cix->RootPointer()=p_win;
    else {
      assert(root==1);//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      p_cix->RootPointer()=p_cha->RootPointer();
      p_cha->RootPointer()=*i_run;//???????????????????????????????????????????
    }
  }

  (*i_run)->RenewBranch(*m_box.pQua);
  p_cha->BranchPointer()=m_box.pQua;
  p_cha->GlubranchPointerList().remove(m_box.pGlu);
  delete m_box.pGlu; m_box.pGlu=NULL;

  //Check 4-momenta.
  Vec4D testcix, testcha, testall;
#ifdef CHAIN_HANDLER_OUTPUT
  cout<<p_cix->CheckMomentumConservation(testcix)<<"\t"<<testcix<<"\n";
  cout<<p_cha->CheckMomentumConservation(testcha)<<"\t"<<testcha<<"\n";
#else
  p_cix->CheckMomentumConservation(testcix);
  p_cha->CheckMomentumConservation(testcha);
#endif
  testall=testcix+testcha+(-1)*m_old;
#ifdef CHAIN_HANDLER_OUTPUT
  cout<<"\t"<<testall<<"\n";
  cout<<"\t"<<p_cix->Momentum()+p_cha->Momentum()+(-1)*m_old<<"\n";
#endif
  testcix+=(-1)*p_cix->Momentum(); testcha+=(-1)*p_cha->Momentum();
#ifdef CHAIN_HANDLER_OUTPUT
  cout<<testcix<<"\n";
  cout<<testcha<<"\n";
#endif
  for(char i=0; i<4; ++i) assert(dabs(testcix[i])<1.0e-9 &&
				 dabs(testcha[i])<1.0e-9 &&
				 dabs(testall[i])<1.0e-9);

  p_cha->SetLastScale()=m_k2tcomp;
  p_cix->SetLastScale()=m_k2tcomp;
  p_cix->SetStatus()=On;

#ifdef CHAIN_HANDLER_OUTPUT
  p_cix->Print();
  p_cha->Print();
#endif

  //cout<<endl; cin>>enter; cout<<endl;////////////////////////////////////////

}





void Chain_Handler::TopSplitRing() {

  //There are check asserts to find!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef CHAIN_HANDLER_OUTPUT
  cout<<__PRETTY_FUNCTION__;
  p_cha->Print();
#endif

  f_below=false;    //Defined setting.
  p_cix=NULL;

  //m_box.pAqu is the new antiquark.
  //m_box.pQua is the new quark.
  //m_box.pGlu is the split gluon.

  assert(p_win->IsType()==Dipole::qg);//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  m_vec=m_box.pAqu->Momentum();
  m_vec+=m_box.pQua->Momentum();
  m_vec+=p_win->GetBotBranchPointer()->Momentum();
  p_cha->UpdateMomentum(1.0,m_vec);

  p_cha->BranchPointer()=m_box.pQua;
  p_cha->AntibranchPointer()=m_box.pAqu;

  if(p_cha->FirstGlubranchPointer()==m_box.pGlu) {
#ifdef CHAIN_HANDLER_OUTPUT
    cout<<"SPECIAL CASE\n";
#endif
    assert(p_win==p_cha->DipolePointerList().front());//!!!!!!!!!!!!!!!!!!!!!!!
    p_cha->FirstGlubranchPointer()=NULL;
    i_run=p_cha->DipolePointerList().end(); --i_run;
    (*i_run)->RenewBranch(*m_box.pAqu);
    p_cha->GlubranchPointerList().remove(m_box.pGlu);
    delete m_box.pGlu; m_box.pGlu=NULL;
  } else {
    i_run=i_fix; --i_run;
    (*i_run)->RenewBranch(*m_box.pAqu);
    p_cha->GlubranchPointerList().remove(m_box.pGlu);
    delete m_box.pGlu; m_box.pGlu=NULL;
    p_cha->DipolePointerList().splice(p_cha->DipolePointerList().end(),
				      p_cha->DipolePointerList(),
				      p_cha->DipolePointerList().begin(),
				      i_fix);
  }

  m_typ=p_cha->ChainType();
  assert(m_typ==Chain::line);//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  //Check 4-momenta.
  Vec4D testcha, testall;
  p_cha->CheckMomentumConservation(testcha);
  testall=testcha+(-1)*m_old;
  testcha+=(-1)*p_cha->Momentum();
#ifdef CHAIN_HANDLER_OUTPUT
  cout<<"\t"<<testall<<"\n";
  cout<<"\t"<<testcha<<"\n";
#endif
  for(char i=0; i<4; ++i) assert(dabs(testcha[i])<1.0e-9 &&
				 dabs(testall[i])<1.0e-9);

  p_cha->SetLastScale()=m_k2tcomp;

#ifdef CHAIN_HANDLER_OUTPUT
  p_cha->Print();
#endif

  //cout<<endl; cin>>enter; cout<<endl;////////////////////////////////////////

}





void Chain_Handler::BotSplitRing() {

  //There are check asserts to find!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef CHAIN_HANDLER_OUTPUT
  cout<<__PRETTY_FUNCTION__;
  p_cha->Print();
#endif

  f_below=false;    //Defined setting.
  p_cix=NULL;

  //m_box.pAqu is the new antiquark.
  //m_box.pQua is the new quark.
  //m_box.pGlu is the split gluon.

  assert(p_win->IsType()==Dipole::gqbar);//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  m_vec=p_win->GetTopBranchPointer()->Momentum();
  m_vec+=m_box.pAqu->Momentum();
  m_vec+=m_box.pQua->Momentum();
  p_cha->UpdateMomentum(1.0,m_vec);

  p_cha->BranchPointer()=m_box.pQua;
  p_cha->AntibranchPointer()=m_box.pAqu;

  if(p_cha->FirstGlubranchPointer()==m_box.pGlu) {
#ifdef CHAIN_HANDLER_OUTPUT
    cout<<"SPECIAL CASE\n";
#endif
    assert(p_win==p_cha->DipolePointerList().back());//!!!!!!!!!!!!!!!!!!!!!!!!
    p_cha->FirstGlubranchPointer()=NULL;
    i_run=p_cha->DipolePointerList().begin();
    (*i_run)->RenewBranch(*m_box.pQua);
    p_cha->GlubranchPointerList().remove(m_box.pGlu);
    delete m_box.pGlu; m_box.pGlu=NULL;
  } else {
    i_run=i_fix; ++i_run;
    (*i_run)->RenewBranch(*m_box.pQua);
    p_cha->GlubranchPointerList().remove(m_box.pGlu);
    delete m_box.pGlu; m_box.pGlu=NULL;
    p_cha->DipolePointerList().splice(p_cha->DipolePointerList().end(),
				      p_cha->DipolePointerList(),
				      p_cha->DipolePointerList().begin(),
				      i_run);
  }

  m_typ=p_cha->ChainType();
  assert(m_typ==Chain::line);//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  //Check 4-momenta.
  Vec4D testcha, testall;
  p_cha->CheckMomentumConservation(testcha);
  testall=testcha+(-1)*m_old;
  testcha+=(-1)*p_cha->Momentum();
#ifdef CHAIN_HANDLER_OUTPUT
  cout<<"\t"<<testall<<"\n";
  cout<<"\t"<<testcha<<"\n";
#endif

  for(char i=0; i<4; ++i) assert(dabs(testcha[i])<1.0e-9 &&
				 dabs(testall[i])<1.0e-9);

  p_cha->SetLastScale()=m_k2tcomp;

#ifdef CHAIN_HANDLER_OUTPUT
  p_cha->Print();
#endif

  //cout<<endl; cin>>enter; cout<<endl;////////////////////////////////////////

}



//=============================================================================





//eof
