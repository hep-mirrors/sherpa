//bof
//Version: 2 ADICIC++-0.0/2004/09/06

//Implementation of Dipole_Handler.H.



#include "Random.H"
#include "Poincare.H"
#include "Dipole_Handler.H"
#include "Dipole_Parameter.H"
#include "Sudakov_Calculator.H"
#include "Recoil_Calculator.H"





using namespace std;
using namespace ATOOLS;
using namespace ADICIC;





//#include "....exa.cc"





//=============================================================================



//ostream& ADICIC::operator<<(ostream& ost, const ADICIC::...&) {
//}



//=============================================================================



//So far there is no static Dipole_Handler.
int Dipole_Handler::s_count=0;
const int& Dipole_Handler::InStore=Dipole_Handler::s_count;

Dipole_Handler::Sudakovbox
Dipole_Handler::s_sumap=Dipole_Handler::Sudakovbox();
Dipole_Handler::Recoilbox
Dipole_Handler::s_remap=Dipole_Handler::Recoilbox();

const bool Dipole_Handler::sf_init=Dipole_Handler::AdjustCalcBox();



//=============================================================================



Dipole_Handler::Dipole_Handler()
  : m_key(Dipole::incorrect,Radiation::incorrect),
    //m_key(make_pair(Dipole::incorrect,Radiation::incorrect)),
    p_sudakov(NULL), p_recoil(NULL),
    p_dip(NULL),
    p_dix(NULL), p_ban(NULL), p_ati(NULL), p_glu(NULL),
    f_below(false), f_recoil(Nil), f_gate(0) {
  ++s_count;
}





Dipole_Handler::Dipole_Handler(Dipole& dip)
  : m_key(Dipole::incorrect,Radiation::incorrect),
    //m_key(make_pair(Dipole::incorrect,Radiation::incorrect)),
    p_sudakov(NULL), p_recoil(NULL),
    p_dip(NULL),
    p_dix(NULL), p_ban(NULL), p_ati(NULL), p_glu(NULL),
    f_below(false), f_recoil(Nil) {

  ++s_count;

  if(dip|*this) {
    if(dip.IsHandledBy(*this)); else {
      cerr<<"\nBug: Wrong Dipole-Dipole_Handler connection emerged!\n";
      assert(dip.IsHandledBy(*this));
    }
  } else {
    cerr<<"\nMethod: ADICIC::Dipole_Handler::Dipole_Handler(ADICIC::Dipole&): "
	<<"Warning: Attaching Dipole failed!\n"<<endl;
  }

}





Dipole_Handler::~Dipole_Handler() {

  --s_count;

  if(p_dix) {
    delete p_dix;
    if(p_glu) delete p_glu;
  } else {
    if(p_ban) delete p_ban;
    if(p_ati) delete p_ati;
  }

  assert(p_sudakov && p_dip || !p_sudakov && !p_recoil && !p_dip);

  if(!p_dip) return;
  if(p_dip->IsHandledBy(*this)==false) {
    cerr<<"\nBug: Wrong Dipole-Dipole_Handler connection emerged!\n";
    assert(p_dip->IsHandledBy(*this));
  }

  *p_dip|0;

}



//=============================================================================



void Dipole_Handler::ShowCalcBox() {    //Static.
  cout<<endl;
  cout<<"======================================="<<endl;
  cout<<"Calculator box for the Dipole_Handler's"<<endl;
  cout<<"---------------------------------------"<<endl;
  cout<<"Number of Sudakov_Calculators in store = "
      <<Sudakov_Calculator::InStore<<"."<<endl;
  cout<<"Mimic Ariadne? "<<Sudakov_Calculator::Ariadne<<"."<<endl;
  cout<<"Number of  Recoil_Calculators in store = "
      <<Recoil_Calculator::InStore<<"."<<endl;
  cout<<"---------------------------------------"<<endl;
  s_sumap[Dipole::qqbar]->Which();
  s_sumap[Dipole::qqbar]->ShowSpecification();
  s_sumap[Dipole::qg]->Which();
  s_sumap[Dipole::qg]->ShowSpecification();
  s_sumap[Dipole::gqbar]->Which();
  s_sumap[Dipole::gqbar]->ShowSpecification();
  s_sumap[Dipole::gg]->Which();
  s_sumap[Dipole::gg]->ShowSpecification();
  cout<<"---------------------------------------"<<endl;
  if(s_remap[Key(Dipole::qqbar,Radiation::gluon)])
    s_remap[Key(Dipole::qqbar,Radiation::gluon)]->Which();
  if(s_remap[Key(Dipole::qg,Radiation::gluon)])
    s_remap[Key(Dipole::qg,Radiation::gluon)]->Which();
  if(s_remap[Key(Dipole::gqbar,Radiation::gluon)])
    s_remap[Key(Dipole::gqbar,Radiation::gluon)]->Which();
  if(s_remap[Key(Dipole::gg,Radiation::gluon)])
    s_remap[Key(Dipole::gg,Radiation::gluon)]->Which();
  cout<<" -  -  -  -  -  -  -  -  -  -  -  -  - "<<endl;
  if(s_remap[Key(Dipole::qg,Radiation::qbot)])
    s_remap[Key(Dipole::qg,Radiation::qbot)]->Which();
  if(s_remap[Key(Dipole::gqbar,Radiation::qtop)])
    s_remap[Key(Dipole::gqbar,Radiation::qtop)]->Which();
  if(s_remap[Key(Dipole::gg,Radiation::qbot)])
    s_remap[Key(Dipole::gg,Radiation::qbot)]->Which();
  if(s_remap[Key(Dipole::gg,Radiation::qtop)])
    s_remap[Key(Dipole::gg,Radiation::qtop)]->Which();
  cout<<"======================================="<<endl;
}





const bool Dipole_Handler::AdjustCalcBox() {    //Static.

  //hier muss noch mehr gebohrt werden: group radiation type
  //sudakov strategy mit parametern handeln.

  static bool firsttime=true;
  static Calcbox box;    //Needed for a proper deletion at the very very end.

  if(firsttime) {

    firsttime=false;

    Dipole_Parameter::ForceFirstInit();

    const Radiation::Type raty=Radiation::g;///////////////////////////////////

    //Arrange the Sudakov's.
    box.p_sud[0]=new Sudakov_Group<Dipole::qqbar>(raty);
    box.p_sud[1]=new Sudakov_Group<Dipole::qg>(raty);
    box.p_sud[2]=new Sudakov_Group<Dipole::gqbar>(raty);
    box.p_sud[3]=new Sudakov_Group<Dipole::gg>(raty);
    for(short i=0; i<4; ++i) assert(box.p_sud[i]);

    //Fix the Sudakov map.
    s_sumap[Dipole::qqbar]=box.p_sud[0];
    s_sumap[Dipole::qg]   =box.p_sud[1];
    s_sumap[Dipole::gqbar]=box.p_sud[2];
    s_sumap[Dipole::gg]   =box.p_sud[3];

    //Establish the overall recoil strategy right now and here.
    if(raty>5) {
      box.p_rec[0]=new Recoil<Recoil_Strategy::Ret_qgqbar>;
      box.p_rec[1]=new Recoil<Recoil_Strategy::Ret_qgg>;
      box.p_rec[2]=new Recoil<Recoil_Strategy::Ret_ggqbar>;
      box.p_rec[3]=new Recoil<Recoil_Strategy::Ret_ggg>;
      for(short i=0; i<4; ++i) assert(box.p_rec[i]);
    }
    if(raty!=Radiation::g) {
      box.p_rec[4]=new Recoil<Recoil_Strategy::Ret_qqbarq>;
      box.p_rec[5]=new Recoil<Recoil_Strategy::Ret_qbarqqbar>;
      box.p_rec[6]=new Recoil<Recoil_Strategy::Ret_gqbarq>;
      box.p_rec[7]=new Recoil<Recoil_Strategy::Ret_qbarqg>;
      for(short i=4; i<8; ++i) assert(box.p_rec[i]);
    }

    //Fix the Recoil map.
    s_remap[Key(Dipole::qqbar,Radiation::gluon)] = box.p_rec[0];
    s_remap[Key(Dipole::qg,Radiation::gluon)]    = box.p_rec[1];
    s_remap[Key(Dipole::gqbar,Radiation::gluon)] = box.p_rec[2];
    s_remap[Key(Dipole::gg,Radiation::gluon)]    = box.p_rec[3];
    s_remap[Key(Dipole::qg,Radiation::qbot)]     = box.p_rec[4];
    s_remap[Key(Dipole::gqbar,Radiation::qtop)]  = box.p_rec[5];
    s_remap[Key(Dipole::gg,Radiation::qbot)]     = box.p_rec[6];
    s_remap[Key(Dipole::gg,Radiation::qtop)]     = box.p_rec[7];

    //Defined settings.
    s_remap[Key(Dipole::qqbar,Radiation::incorrect)] = NULL;
    s_remap[Key(Dipole::qg,Radiation::incorrect)]    = NULL;
    s_remap[Key(Dipole::gqbar,Radiation::incorrect)] = NULL;
    s_remap[Key(Dipole::gg,Radiation::incorrect)]    = NULL;

#ifdef DIPOLE_HANDLER_OUTPUT
    cout<<"ADICIC::Dipole_Handler: Calcbox is now initialized.\n";
#endif

    return true;

  }

  if(s_count) {
    cerr<<"\nStatic method: "
	<<"const bool ADICIC::Dipole_Handler::AdjustCalcBox(): "
	<<"Warning: Re-adjusting is not permitted "
	<<"since Dipole_Handler's are already present!\n"<<endl;
    return false;
  }

  const Radiation::Type newraty=Radiation::gduscb;/////////////////////////////

  for(short i=0; i<4; ++i) delete box.p_sud[i];
  box.p_sud[0]=new Sudakov_Group<Dipole::qqbar>(newraty);
  box.p_sud[1]=new Sudakov_Group<Dipole::qg>(newraty);
  box.p_sud[2]=new Sudakov_Group<Dipole::gqbar>(newraty);
  box.p_sud[3]=new Sudakov_Group<Dipole::gg>(newraty);
  for(short i=0; i<4; ++i) assert(box.p_sud[i]);
  //Re-fix the Sudakov map.
  s_sumap[Dipole::qqbar] = box.p_sud[0];
  s_sumap[Dipole::qg]    = box.p_sud[1];
  s_sumap[Dipole::gqbar] = box.p_sud[2];
  s_sumap[Dipole::gg]    = box.p_sud[3];

  for(short i=0; i<8; ++i)
    if(box.p_rec[i]) { delete box.p_rec[i]; box.p_rec[i]=NULL;}
  if(newraty>5) {
    box.p_rec[0]=
      ReadjustRecoilStrategy(Dipole_Parameter::RecoilStrategyQQbar());
    box.p_rec[1]=
      ReadjustRecoilStrategy(Dipole_Parameter::RecoilStrategyQG());
    box.p_rec[2]=
      ReadjustRecoilStrategy(Dipole_Parameter::RecoilStrategyGQbar());
    box.p_rec[3]=
      ReadjustRecoilStrategy(Dipole_Parameter::RecoilStrategyGG());
    for(short i=0; i<4; ++i) assert(box.p_rec[i]);
  }
  if(newraty!=Radiation::g) {
    box.p_rec[4]=
      ReadjustRecoilStrategy(Dipole_Parameter::RecoilStrategyGQbar());///////
    box.p_rec[5]=
      ReadjustRecoilStrategy(Dipole_Parameter::RecoilStrategyQG());//////////
    box.p_rec[6]=
      ReadjustRecoilStrategy(Dipole_Parameter::RecoilStrategyGQbar());///////
    box.p_rec[7]=
      ReadjustRecoilStrategy(Dipole_Parameter::RecoilStrategyQG());//////////
    for(short i=4; i<8; ++i) assert(box.p_rec[i]);
  }
  //Re-fix the Recoil map.
  s_remap[Key(Dipole::qqbar,Radiation::gluon)] = box.p_rec[0];
  s_remap[Key(Dipole::qg,Radiation::gluon)]    = box.p_rec[1];
  s_remap[Key(Dipole::gqbar,Radiation::gluon)] = box.p_rec[2];
  s_remap[Key(Dipole::gg,Radiation::gluon)]    = box.p_rec[3];
  s_remap[Key(Dipole::qg,Radiation::qbot)]     = box.p_rec[4];
  s_remap[Key(Dipole::gqbar,Radiation::qtop)]  = box.p_rec[5];
  s_remap[Key(Dipole::gg,Radiation::qbot)]     = box.p_rec[6];
  s_remap[Key(Dipole::gg,Radiation::qtop)]     = box.p_rec[7];

#ifdef DIPOLE_HANDLER_OUTPUT
  cout<<"ADICIC::Dipole_Handler: Calcbox is now re-initialized.\n";
#endif

  return true;

}





Recoil_Calculator* Dipole_Handler::ReadjustRecoilStrategy(const int s) {

  //Static.

  switch(s) {
  case  1: return new Recoil<Recoil_Strategy::FixDir1>;
  case  2: return new Recoil<Recoil_Strategy::Kleiss>;
  case  3: return new Recoil<Recoil_Strategy::FixDir3>;
  case  4: return new Recoil<Recoil_Strategy::MinimizePt>;
  case  5: return new Recoil<Recoil_Strategy::Lonnblad>;
  case  6: return new Recoil<Recoil_Strategy::OldAdicic>;
  case  7: return new Recoil<Recoil_Strategy::Test>;
  default: return new Recoil<Recoil_Strategy::Unknown>;
  }

}



//=============================================================================



void Dipole_Handler::RemoveNewProducts() {
  //Resets the news.
  f_below=false;
  f_recoil=Nil;
  if(p_dix) {
    delete p_dix; p_dix=NULL;
    if(p_glu) { delete p_glu; p_glu=NULL;}
  } else {
    if(p_ban) { delete p_ban; p_ban=NULL;}
    if(p_ati) { delete p_ati; p_ati=NULL;}
    p_glu=NULL;
  }
  assert(p_dix==NULL && p_glu==NULL && p_ban==NULL && p_ati==NULL);
}





void Dipole_Handler::ShowSudakov() const {
  cout<<endl;
  cout<<"=========================================="<<endl;
  cout<<"Sudakov_Calculator for this Dipole_Handler"<<endl;
  cout<<"------------------------------------------"<<endl;
  if(p_sudakov) { p_sudakov->Which(); p_sudakov->ShowSpecification();}
  else cout<<"Not initialized."<<endl;
  cout<<"Number of Sudakov_Calculators in store = "
      <<Sudakov_Calculator::InStore<<"."<<endl;
  cout<<"Mimic Ariadne? "<<Sudakov_Calculator::Ariadne<<"."<<endl;
  cout<<"=========================================="<<endl;
}





void Dipole_Handler::ShowRecoil() const {
  cout<<endl;
  cout<<"========================================="<<endl;
  cout<<"Recoil_Calculator for this Dipole_Handler"<<endl;
  cout<<"-----------------------------------------"<<endl;
  if(p_recoil) p_recoil->Which();
  else cout<<"Not initialized."<<endl;
  cout<<"Number of Recoil_Calculators in store = "
      <<Recoil_Calculator::InStore<<"."<<endl;
  cout<<"========================================="<<endl;
}





const bool Dipole_Handler::InduceDipoleRadiation() {

  f_gate=0;

  if(!p_dip) return false;

  if(p_dip->Status()==On && p_dip->PointerHandling()==0 &&
     p_dip->IsType()!=Dipole::incorrect);
  else { p_dip->SetEmitScale()=0.0; return false;}

  if(p_dix) {
    delete p_dix; p_dix=NULL;
    if(p_glu) { delete p_glu; p_glu=NULL;}
  } else {
    if(p_ban) { delete p_ban; p_ban=NULL;}
    if(p_ati) { delete p_ati; p_ati=NULL;}
    p_glu=NULL;
  }

  assert(p_dix==NULL && p_glu==NULL && p_ban==NULL && p_ati==NULL);

  //No testing of global parameters.
  //assert( p_dip->InvMass() > Sudakov_Calculator::MinOfK2t() );
  //assert( p_dip->ProdScale() > Sudakov_Calculator::MinOfK2t() );

  Sudakov_Strategy::Factorization factstrat;///////////////////////////////////

  if( p_sudakov->GenerateEfracsFor(*p_dip,factstrat) ) {
    p_sudakov->GetResult(m_sur);
    m_key.second=m_sur.Rad;
#ifdef DIPOLE_HANDLER_OUTPUT
    cout<<"(("<<m_key.first<<","<<m_key.second<<"))"<<endl;
#endif
    p_recoil=s_remap[m_key];
    p_dip->SetEmitScale()=m_sur.P2t;
    f_gate=p_dip->StateNumber;
  }
  else {
    p_dip->SetEmitScale()=Sudakov_Calculator::MinOfK2t();
    return false;
  }

#ifdef DIPOLE_HANDLER_OUTPUT
  cout<<"\ttransverse momentum and energy fractions:\n\t\t p2t=";
  cout<<m_sur.P2t<<endl;
  cout<<"\t\t x1="<<m_sur.X1<<endl;
  cout<<"\t\t x3="<<m_sur.X3<<endl;
#endif

  return true;

}





const bool Dipole_Handler::FinishDipoleRadiation() {

  //assert(p_dip);
  assert(f_gate);
  assert(f_gate==p_dip->StateNumber);

  if(m_sur.P2t!=p_dip->EmitScale() || this->Status()!=bool4::zero) {
    f_gate=0; return false;
  }

  assert((p_dip->TotP())[0] > 0.0);

  m_p1=p_dip->GetTopBranchPointer()->Momentum();
  m_p3=p_dip->GetBotBranchPointer()->Momentum();

  assert(m_p1.Abs2() > -1.0e-6);
  assert(m_p3.Abs2() > -1.0e-6);
  assert(m_p1[0] > 0.0);
  assert(m_p3[0] > 0.0);

  assert(GenerateMomenta());
  assert(GenerateSplitting());

  //Due to the dipole settings in GenerateSplitting,
  //actually the following should not be necessary.
  f_gate=0;

  return true;

}





const bool Dipole_Handler::ManageDipoleRadiation() {

  f_gate=0;

  if(!p_dip) return false;

  if(p_dip->Status()==On && p_dip->PointerHandling()==0 &&
     p_dip->IsType()!=Dipole::incorrect);
  else { p_dip->SetEmitScale()=0.0; return false;}

  if(p_dix) {
    delete p_dix; p_dix=NULL;
    if(p_glu) { delete p_glu; p_glu=NULL;}
  } else {
    if(p_ban) { delete p_ban; p_ban=NULL;}
    if(p_ati) { delete p_ati; p_ati=NULL;}
    p_glu=NULL;
  }

  assert(p_dix==NULL && p_glu==NULL && p_ban==NULL && p_ati==NULL);

  //No testing of global parameters.
  assert( p_dip->InvMass() > Sudakov_Calculator::MinOfK2t() );
  assert( p_dip->ProdScale() > Sudakov_Calculator::MinOfK2t() );

  Sudakov_Strategy::Factorization factstrat;///////////////////////////////////

  if( p_sudakov->GenerateEfracsFor(*p_dip,factstrat) ) {
    p_sudakov->GetResult(m_sur);
    m_key.second=m_sur.Rad;
#ifdef DIPOLE_HANDLER_OUTPUT
    cout<<"(("<<m_key.first<<","<<m_key.second<<"))"<<endl;
#endif
    p_recoil=s_remap[m_key];
    p_dip->SetEmitScale()=m_sur.P2t;
  }
  else {
    p_dip->SetEmitScale()=Sudakov_Calculator::MinOfK2t();
    return false;
  }

#ifdef DIPOLE_HANDLER_OUTPUT
  cout<<"\ttransverse momentum and energy fractions:\n\t\t p2t=";
  cout<<m_sur.P2t<<endl;
  cout<<"\t\t x1="<<m_sur.X1<<endl;
  cout<<"\t\t x3="<<m_sur.X3<<endl;
#endif

  assert((p_dip->TotP())[0] > 0.0);

  m_p1=p_dip->GetTopBranchPointer()->Momentum();
  m_p3=p_dip->GetBotBranchPointer()->Momentum();

  assert(m_p1.Abs2() > -1.0e-12);
  assert(m_p3.Abs2() > -1.0e-12);
  assert(m_p1[0] > 0.0);
  assert(m_p3[0] > 0.0);

  assert(GenerateMomenta());
  assert(GenerateSplitting());

  return true;

}



//=============================================================================



const bool Dipole_Handler::GenerateMomenta() {

  const Vec4D& Plab=p_dip->TotP();

#ifdef DIPOLE_HANDLER_OUTPUT
  cout<<"\tlab frame - before:\n\t\t P =";
  cout<<Plab<<"\t"<<p_dip->InvMass()<<"  "<<p_dip->Mass()<<endl;
  cout<<"\t\t q1="<<m_p1<<"\t "<<m_p1.Abs2()<<endl;
  cout<<"\t\t q3="<<m_p3<<"\t "<<m_p3.Abs2()<<endl;
#endif

  Recoil_Setup Iset;
  Iset.E2=sqrt(p_dip->InvMass());
  Iset.E1=0.5*Iset.E2*m_sur.X1;
  Iset.E3=0.5*Iset.E2*m_sur.X3;
  Iset.E2=Iset.E2-Iset.E1-Iset.E3;

  Vec4D& axis=m_p2;
  axis=m_p1;    //lab frame
  Poincare fly(Plab);
  fly.Boost(axis);    //This is always the initial cms frame axis.

  if(p_recoil->GenerateCmsMomenta(Iset,axis))
    p_recoil->GetResult(f_recoil,m_p1,m_p3);
  else return false;

#ifdef DIPOLE_HANDLER_OUTPUT
  cout<<"\tcms frame - after:\n";
  cout<<"\t\t p1="<<m_p1<<"\t "<<m_p1.Abs2()<<endl;
  cout<<"\t\t p3="<<m_p3<<"\t "<<m_p3.Abs2()<<endl;
#endif

  if(TEMP::CPTEST) CrossProductTest(axis);/////////////////////////////////////

  fly.BoostBack(m_p1);
  fly.BoostBack(m_p3);

  m_p2=Plab+(-1.0)*(m_p1+m_p3);

#ifdef DIPOLE_HANDLER_OUTPUT
  cout<<"\tlab frame - after:\n";
  cout<<"\t\t p1="<<m_p1<<"\t "<<m_p1.Abs2()<<endl;
  cout<<"\t\t p2="<<m_p2<<"\t "<<m_p2.Abs2()<<endl;
  cout<<"\t\t p3="<<m_p3<<"\t "<<m_p3.Abs2()<<endl;
#endif

  return true;

}





const bool Dipole_Handler::GenerateSplitting() {

  switch(m_sur.Rad) {

  case Radiation::gluon: {
    p_dip->GetTopBranchPointer()->SetMomentum(m_p1);
    p_dip->GetBotBranchPointer()->SetMomentum(m_p3);
    //That updates the dipole as well as the neighbouring ones.
    p_glu=new Dipole::Glubranch(m_p2);
    assert(p_glu); assert(p_glu->Flav()==Flavour(m_sur.Kfc));//////////////////
    p_dix=new Dipole(*p_dip);
    assert(p_dix);
    if(f_recoil==Positive) {
      f_below=false;
      p_dix->RenewBranch(false,*p_glu);    //dixbot
      p_dip->RenewBranch(true,*p_glu);    //diptop
    } else {
      f_below=true;
      p_dip->RenewBranch(false,*p_glu);    //dipbot
      p_dix->RenewBranch(true,*p_glu);    //dixtop
    }
    p_dix->SetSource()=p_dip->Name;
    p_dix->SetProdScale()=m_sur.P2t;
    p_dix->SetBootScale()=m_sur.P2t;
    p_dix->SetEmitScale()=m_sur.P2t;
    break;
  }

  case Radiation::qtop: {
    f_below=false;
    p_dip->GetBotBranchPointer()->SetMomentum(m_p3);
    //That updates the dipole and if existing the neighbouring one below.
    Dipole_Particle* topglu=p_dip->GetTopBranchPointer().operator->();
    assert(topglu->OrgType()==Nil);
    p_glu=static_cast<Dipole::Glubranch*>(topglu); assert(p_glu);
    p_ati=new Dipole::Antibranch(interface.antiq[m_sur.Kfc],m_p1);
    assert(p_ati); assert(p_ati->Flav()==Flavour(m_sur.Kfc,1));////////////////
    p_ban=new Dipole::Branch(interface.quark[m_sur.Kfc],m_p2);
    assert(p_ban); assert(p_ban->Flav()==Flavour(m_sur.Kfc));//////////////////
    p_dip->RenewBranch(*p_ban);    //dipbot
    break;
  }

  case Radiation::qbot: {
    f_below=true;
    p_dip->GetTopBranchPointer()->SetMomentum(m_p1);
    //That updates the dipole and if existing the neighbouring one above.
    Dipole_Particle* botglu=p_dip->GetBotBranchPointer().operator->();
    assert(botglu->OrgType()==Nil);
    p_glu=static_cast<Dipole::Glubranch*>(botglu); assert(p_glu);
    p_ati=new Dipole::Antibranch(interface.antiq[m_sur.Kfc],m_p2);
    assert(p_ati); assert(p_ati->Flav()==Flavour(m_sur.Kfc,1));////////////////
    p_ban=new Dipole::Branch(interface.quark[m_sur.Kfc],m_p3);
    assert(p_ban); assert(p_ban->Flav()==Flavour(m_sur.Kfc));//////////////////
    p_dip->RenewBranch(*p_ati);    //dipbot
    break;
  }

  default:
    assert(m_sur.Rad==Radiation::gluon || m_sur.Rad==Radiation::qtop ||
	   m_sur.Rad==Radiation::qbot);

  }

  p_dip->SetProdScale()=m_sur.P2t;
  p_dip->SetBootScale()=m_sur.P2t;
  p_dip->SetEmitScale()=m_sur.P2t;

  m_key.first=p_dip->IsType();
  m_key.second=Radiation::incorrect;
  p_sudakov=s_sumap[m_key.first];
  p_recoil=s_remap[m_key];

  return true;

}



//=============================================================================



void Dipole_Handler::CrossProductTest(const Vec4D& axis) const {
  Vec3D q1(m_p1);
  Vec3D q3(m_p3);
  Vec3D ax(axis);
  Vec3D B=cross(q1,q3);
  Vec3D A=cross(ax,q1);
  cout<<"  Cross product test."<<endl;
  cout<<"  ax="<<ax<<endl;
  cout<<"  q1="<<q1<<endl;
  cout<<"  q3="<<q3<<endl;
  cout<<"   A="<<A<<endl;
  cout<<"   B="<<B<<endl;
  cout<<"     ";
  for(char i=1; i<4; ++i) cout<<A[i]/B[i]<<" : ";
  cout<<"     : ";
  for(char i=1; i<4; ++i) cout<<B[i]/A[i]<<" : ";
  cout<<endl;
  cout<<"  +++++++++++++++++++"<<endl;
}



//=============================================================================



Dipole_Handler::Calcbox::Calcbox() {
  for(short i=0; i<4; ++i) { p_sud[i]=NULL; p_rec[i]=NULL;}
  for(short i=4; i<8; ++i) p_rec[i]=NULL;
}





Dipole_Handler::Calcbox::~Calcbox() {
  for(short i=0; i<4; ++i) {
    if(p_sud[i]) delete p_sud[i];
    if(p_rec[i]) delete p_rec[i];
  }
  for(short i=4; i<8; ++i) if(p_rec[i]) delete p_rec[i];
}



//=============================================================================





//eof
