//bof
//Version: 3 ADICIC++-0.0/2005/09/21

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
template
const bool Chain_Handler::FindTheDipole<Chain_Evolution_Strategy::Emission>();
template
const bool Chain_Handler::FindTheDipole<Chain_Evolution_Strategy::Mass>();



//=============================================================================



Chain_Handler::Chain_Handler()
  : f_below(false), m_code(), p_cix(NULL), p_rec(NULL),
    m_k2tcomp(0.0), p_cha(NULL),
    m_dh1(), m_dh2(),
    p_dhwait(&m_dh1), p_dhaciv(&m_dh2), p_dhtemp(NULL),
    i_fix(NULL), i_run(NULL) {

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
  : f_below(false), m_code(), p_cix(NULL), p_rec(NULL),
    m_k2tcomp(0.0), p_cha(NULL),
    m_dh1(), m_dh2(),
    p_dhwait(&m_dh1), p_dhaciv(&m_dh2), p_dhtemp(NULL),
    i_fix(NULL), i_run(NULL) {

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

  if(p_rec) delete p_rec;
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

  this->CleanUp();
  this->RemoveNewProducts();    //No testing of global parameters.

  if(this->FindDipole()) return this->ModifyChain();

  return false;

}





const bool Chain_Handler::EvolveChain() {//////////////////////////////////////
  return false;////////////////////////////////////////////////////////////////
}//////////////////////////////////////////////////////////////////////////////





const bool Chain_Handler::CorrectChain(const Recoil_Tool& reto) {
  if(!p_cha) return false;
  if(!p_cha->IsLine()) return false;
  if(reto.Mode()!=rdt::iirecoil) return false;
  Poincare& fly=reto.GetBoost();
  Poincare& flyprime=reto.GetBackBoost();
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
    fly.Boost(temp);
    flyprime.BoostBack(temp);
    (*it)->SetMomentum(temp);    //Updates dipole mass!
  }
  return false;////////////////////////////////////////////////////////////////
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





void Chain_Handler::CleanUp() {
  assert(!m_dh1.IsDocked() && !m_dh2.IsDocked());
  m_dh1.RemoveNewProducts();
  m_dh2.RemoveNewProducts();
  p_dhwait=&m_dh1;
  p_dhaciv=&m_dh2;
  p_dhtemp=NULL;
  i_fix=NULL;
  i_run=NULL;
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

  i_run=p_cha->DipolePointerList().begin();
  bool spico=p_cha->ParticleNumber()<dpa.evo.ChainCorrelationLimit();

  for(; i_run!=p_cha->DipolePointerList().end(); ++i_run) {

    Dipole& dip=**i_run;
    dip.SetSpinCorr()=spico;
    assert(dip|*p_dhaciv);

    if(p_dhaciv->InduceDipoleRadiation()==isfalse) {    //Removes the warning.
      dip|0; continue;
    }

    if(dip.EmitScale()>=p_cha->LastScale()) {
      bool result;
      do {
	dip.SetBootScale()=dip.EmitScale();
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
      dip|0; continue;
    }

    m_k2tcomp=dip.EmitScale();
    p_dhtemp=p_dhwait;
    p_dhwait=p_dhaciv;
    p_dhaciv=p_dhtemp;
    if(p_dhtemp->IsDocked()) {
      Dipole& loser=**i_fix;
      loser.SetBootScale()=loser.ProdScale();
      loser|0;
    }
    i_fix=i_run;

  }

  return bool(i_fix!=NULL);

}



//=============================================================================



const bool Chain_Handler::ModifyChain() {

  p_dw=NULL; p_gw=NULL; p_aw=NULL; p_bw=NULL;

  p_win=*i_fix;
  m_nii=-7;

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
  if(m_nii>=0) m_vii=-1*m_vec;

#ifdef CHAIN_HANDLER_OUTPUT
  cout<<m_old<<"\t"<<m_old.Abs2()<<"\n";
  cout<<m_vec<<"\t"<<m_vec.Abs2()<<"\n";
  cout<<p_cha->Momentum()<<"\t"
      <<p_cha->Momentum().Abs2()<<"\t"
      <<p_cha->InvMass()<<endl;
#endif

  assert(p_dhwait->FinishDipoleRadiation());

  p_dhwait->DecoupleNew(p_dw,p_gw,p_aw,p_bw,f_bot,m_rec);
  if(p_dw) {
    if(p_bw) {
      assert(p_aw && p_gw);
      if(f_bot) m_code=p_aw->Flav();
      else      m_code=p_bw->Flav();
      (*p_win)|0;
      EmitQuark();
    } else {
      assert(!p_aw && p_gw);
      m_code=Flavour(kf::gluon);
      (*p_win)|0;
      EmitGluon();
    }
  } else {
    assert(p_aw && p_bw && p_gw);
    m_code=p_bw->Flav();
    assert(p_aw->Flav().Bar()==m_code);
    (*p_win)|0;
    if(m_typ==Chain::line) {
      if(f_bot){
#ifdef CHAIN_HANDLER_OUTPUT
	cout<<p_win->GetTopBranchPointer()->Momentum()<<"\n";
	cout<<p_aw->Momentum()<<"\n";
	cout<<p_bw->Momentum()<<"\n";
	cout<<p_win->GetTopBranchPointer()->Momentum()+
	  p_aw->Momentum()+p_bw->Momentum()<<endl;
#endif
	BotSplitLine();
      }
      else {
#ifdef CHAIN_HANDLER_OUTPUT
	cout<<p_aw->Momentum()<<"\n";
	cout<<p_bw->Momentum()<<"\n";
	cout<<p_win->GetBotBranchPointer()->Momentum()<<"\n";
	cout<<p_win->GetBotBranchPointer()->Momentum()+
	  p_aw->Momentum()+p_bw->Momentum()<<endl;
#endif
	TopSplitLine();
      }
    } else {
      if(f_bot) BotSplitRing();
      else      TopSplitRing();
    }
  }

  return true;

}





void Chain_Handler::DistributeIIRecoil(const Vec4D& viiw) {
  //Must be executed before new particles are introduced to the chain.
  assert(m_nii==0 && p_rec==NULL);
  p_rec=new Recoil_Tool(rdt::iirecoil);
  assert(dabs(m_vii.Abs2()-viiw.Abs2())<1.0e-7);
  assert(m_vii.Abs2()>1.0e-12);
  Poincare* pfly=new Poincare(m_vii); assert(pfly);
  assert(p_rec->KeepThisBoost(*pfly));
  Poincare* pflyprime=new Poincare(viiw); assert(pflyprime);
  assert(p_rec->KeepThisBackBoost(*pflyprime));
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
    pfly->Boost(temp);
    pflyprime->BoostBack(temp);
    (*it)->SetMomentum(temp);    //Updates dipole mass!
    m_vec+=temp;
  }
}





void Chain_Handler::EmitGluon() {

  m_vec=p_gw->Momentum();    //This is emitted, so clearly it's an F gluon.

  if(f_bot) {
    i_run=i_fix;
    ++i_run;
    i_run=p_cha->DipolePointerList().insert(i_run,p_dw);
    //if(p_win->GetTopBranchPointer()->Incoming())
    //if(p_dw->GetBotBranchPointer()->Incoming())
    if(m_nii>=0) {
      m_vec-=p_win->GetTopBranchPointer()->Momentum();
      m_vec-=p_dw->GetBotBranchPointer()->Momentum();
      DistributeIIRecoil(-1*m_vec);    //-1*m_vec corresponds to an m_viinew!
    } else {
      m_vec+=p_win->GetTopBranchPointer()->Momentum();
      m_vec+=p_dw->GetBotBranchPointer()->Momentum();
    }
  } else {
    i_run=p_cha->DipolePointerList().insert(i_fix,p_dw);
    ++i_run;
    i_fix=i_run;
    //if(p_dw->GetTopBranchPointer()->Incoming())
    //if(p_win->GetBotBranchPointer()->Incoming())
    if(m_nii>=0) {
      m_vec-=p_dw->GetTopBranchPointer()->Momentum();
      m_vec-=p_win->GetBotBranchPointer()->Momentum();
      DistributeIIRecoil(-1*m_vec);    //-1*m_vec corresponds to an m_viinew!
    } else {
      m_vec+=p_dw->GetTopBranchPointer()->Momentum();
      m_vec+=p_win->GetBotBranchPointer()->Momentum();
    }
  }

  p_cha->GlubranchPointerList().push_back(p_gw);
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
  m_vec=-1*p_gw->Momentum();    //This is the new initial(=I) gluon.
  //p_cha->GlubranchPointerList().push_back(p_gw);

  if(f_bot) {
    i_run=i_fix;
    ++i_run;
    i_run=p_cha->DipolePointerList().insert(i_run,p_dw);
    assert(p_win->GetTopBranchPointer()->Incoming());
    m_vec-=p_win->GetTopBranchPointer()->Momentum();
    assert(p_dw->GetBotBranchPointer()->Incoming()==false);
    m_vec+=p_dw->GetBotBranchPointer()->Momentum();
    DistributeIIRecoil(-1*m_vec);
    p_cha->AntibranchPointer()=p_aw;    //The new Antibranch of the chain.
    delete p_bw; p_bw=NULL;    //Remove the old one.
  } else {
    i_run=p_cha->DipolePointerList().insert(i_fix,p_dw);
    ++i_run;
    i_fix=i_run;
    assert(p_dw->GetTopBranchPointer()->Incoming()==false);
    m_vec+=p_dw->GetTopBranchPointer()->Momentum();
    assert(p_win->GetBotBranchPointer()->Incoming());
    m_vec-=p_win->GetBotBranchPointer()->Momentum();
    DistributeIIRecoil(-1*m_vec);
    p_cha->BranchPointer()=p_bw;    //The new Branch of the chain.
    delete p_aw; p_aw=NULL;    //Remove the old one.
  }

  p_cha->GlubranchPointerList().push_back(p_gw);
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

  //p_aw is the new antiquark.
  //p_bw is the new quark.
  //p_gw is the split gluon.

  p_cix=new Chain(); assert(p_cix);
  p_cix->m_memo=p_cha->m_memo;
  p_cix->BranchPointer()=p_bw;
  p_cix->AntibranchPointer()=p_cha->AntibranchPointer();
  p_cha->AntibranchPointer()=NULL;

  i_run=p_cha->DipolePointerList().end(); --i_run;

#ifdef CHAIN_HANDLER_OUTPUT
  cout<<"Type: "<<(*i_run)->IsType()<<"\n";
#endif

  if(i_run==i_fix) {    //Special case - only the winner dipole.
    assert(p_win->IsType()==Dipole::qqbar);//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    m_vec=p_bw->Momentum();
    m_vec+=p_cix->AntibranchPointer()->Momentum();
    p_cix->UpdateMomentum(1.0,m_vec);
    p_cha->UpdateMomentum(1.0,p_aw->Momentum());
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
    p_cix->UpdateMomentum(1.0,p_bw->Momentum()+m_vec+
			  p_cix->AntibranchPointer()->Momentum());
    m_vec+=p_cix->AntibranchPointer()->Momentum();
    m_vec*=(-1);
    m_vec+=p_win->GetBotBranchPointer()->Momentum();
    m_vec+=p_aw->Momentum();
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

  (*i_run)->RenewBranch(*p_aw);
  p_cha->AntibranchPointer()=p_aw;
  p_cha->GlubranchPointerList().remove(p_gw);
  delete p_gw; p_gw=NULL;

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

  //p_aw is the new antiquark.
  //p_bw is the new quark.
  //p_gw is the split gluon.

  p_cix=new Chain(); assert(p_cix);
  p_cix->m_memo=p_cha->m_memo;
  p_cix->BranchPointer()=p_cha->BranchPointer();
  p_cha->BranchPointer()=NULL;
  p_cix->AntibranchPointer()=p_aw;

  i_run=p_cha->DipolePointerList().begin();

#ifdef CHAIN_HANDLER_OUTPUT
  cout<<"Type: "<<(*i_run)->IsType()<<"\n";
#endif

  if(i_run==i_fix) {    //Special case - only the winner dipole.
    assert(p_win->IsType()==Dipole::qqbar);//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    m_vec=p_cix->BranchPointer()->Momentum();
    m_vec+=p_aw->Momentum();
    p_cix->UpdateMomentum(1.0,m_vec);
    p_cha->UpdateMomentum(1.0,p_bw->Momentum());
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
			  m_vec+p_aw->Momentum());
    m_vec+=p_cix->BranchPointer()->Momentum();
    m_vec*=(-1);
    m_vec+=p_win->GetTopBranchPointer()->Momentum();
    m_vec+=p_bw->Momentum();
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

  (*i_run)->RenewBranch(*p_bw);
  p_cha->BranchPointer()=p_bw;
  p_cha->GlubranchPointerList().remove(p_gw);
  delete p_gw; p_gw=NULL;

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

  //p_aw is the new antiquark.
  //p_bw is the new quark.
  //p_gw is the split gluon.

  assert(p_win->IsType()==Dipole::qg);//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  m_vec=p_aw->Momentum();
  m_vec+=p_bw->Momentum();
  m_vec+=p_win->GetBotBranchPointer()->Momentum();
  p_cha->UpdateMomentum(1.0,m_vec);

  p_cha->BranchPointer()=p_bw;
  p_cha->AntibranchPointer()=p_aw;

  if(p_cha->FirstGlubranchPointer()==p_gw) {
#ifdef CHAIN_HANDLER_OUTPUT
    cout<<"SPECIAL CASE\n";
#endif
    assert(p_win==p_cha->DipolePointerList().front());//!!!!!!!!!!!!!!!!!!!!!!!
    p_cha->FirstGlubranchPointer()=NULL;
    i_run=p_cha->DipolePointerList().end(); --i_run;
    (*i_run)->RenewBranch(*p_aw);
    p_cha->GlubranchPointerList().remove(p_gw);
    delete p_gw; p_gw=NULL;
  } else {
    i_run=i_fix; --i_run;
    (*i_run)->RenewBranch(*p_aw);
    p_cha->GlubranchPointerList().remove(p_gw);
    delete p_gw; p_gw=NULL;
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

  //p_aw is the new antiquark.
  //p_bw is the new quark.
  //p_gw is the split gluon.

  assert(p_win->IsType()==Dipole::gqbar);//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  m_vec=p_win->GetTopBranchPointer()->Momentum();
  m_vec+=p_aw->Momentum();
  m_vec+=p_bw->Momentum();
  p_cha->UpdateMomentum(1.0,m_vec);

  p_cha->BranchPointer()=p_bw;
  p_cha->AntibranchPointer()=p_aw;

  if(p_cha->FirstGlubranchPointer()==p_gw) {
#ifdef CHAIN_HANDLER_OUTPUT
    cout<<"SPECIAL CASE\n";
#endif
    assert(p_win==p_cha->DipolePointerList().back());//!!!!!!!!!!!!!!!!!!!!!!!!
    p_cha->FirstGlubranchPointer()=NULL;
    i_run=p_cha->DipolePointerList().begin();
    (*i_run)->RenewBranch(*p_bw);
    p_cha->GlubranchPointerList().remove(p_gw);
    delete p_gw; p_gw=NULL;
  } else {
    i_run=i_fix; ++i_run;
    (*i_run)->RenewBranch(*p_bw);
    p_cha->GlubranchPointerList().remove(p_gw);
    delete p_gw; p_gw=NULL;
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
