//bof
//Version: 1 ADICIC++-0.0/2004/05/27

//Implementation of Chain_Handler.H.

///////////////////////////////////////////////////////////

#include "Random.H"
#include "Poincare.H"
#include "Dipole_Handler.H"
#include "Dipole_Handler.dat.cc"

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

Dipole_Handler::Calcbox Dipole_Handler::s_map=Dipole_Handler::Calcbox();
const bool Dipole_Handler::sf_init=Dipole_Handler::InitCalcBox();



//=============================================================================



Dipole_Handler::Dipole_Handler()
  : p_sudakov(NULL), p_recoil(NULL),
    p_dip(NULL),
    p_dix(NULL), p_ban(NULL), p_ati(NULL), p_glu(NULL) {
  ++s_count;
}





Dipole_Handler::Dipole_Handler(Dipole& dip)
  : p_sudakov(NULL), p_recoil(NULL),
    p_dip(NULL),
    p_dix(NULL), p_ban(NULL), p_ati(NULL), p_glu(NULL) {

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

  assert(p_sudakov && p_recoil && p_dip || !p_sudakov && !p_recoil && !p_dip);

  if(!p_dip) return;
  if(p_dip->IsHandledBy(*this)==false) {
    cerr<<"\nBug: Wrong Dipole-Dipole_Handler connection emerged!\n";
    assert(p_dip->IsHandledBy(*this));
  }

  *p_dip|0;

  if(p_dix) delete p_dix;
  if(p_ban) delete p_ban;
  if(p_ati) delete p_ati;
  if(p_glu) delete p_glu;

}





void Dipole_Handler::ShowCalcBox() {    //Static.
  cout<<endl;
  cout<<"======================================="<<endl;
  cout<<"Calculator box for the Dipole_Handler's"<<endl;
  cout<<"---------------------------------------"<<endl;
  cout<<"Number of Sudakov_Calculators in store = "
      <<Sudakov_Calculator::InStore<<"."<<endl;
  cout<<"Number of  Recoil_Calculators in store = "
      <<Recoil_Calculator::InStore<<"."<<endl;
  cout<<"---------------------------------------"<<endl;
  s_map[Dipole::qqbar]->p_sud->Which();
  s_map[Dipole::qqbar]->p_sud->ShowSpecification();
  s_map[Dipole::qg]->p_sud->Which();
  s_map[Dipole::qg]->p_sud->ShowSpecification();
  s_map[Dipole::gqbar]->p_sud->Which();
  s_map[Dipole::gqbar]->p_sud->ShowSpecification();
  s_map[Dipole::gg]->p_sud->Which();
  s_map[Dipole::gg]->p_sud->ShowSpecification();
  cout<<"---------------------------------------"<<endl;
  s_map[Dipole::qqbar]->p_rec->Which();
  s_map[Dipole::qg]->p_rec->Which();
  s_map[Dipole::gqbar]->p_rec->Which();
  s_map[Dipole::gg]->p_rec->Which();
  cout<<"======================================="<<endl;
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
  cout<<"=========================================="<<endl;
}





void Dipole_Handler::ShowRecoilStrategy() const {
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





const bool Dipole_Handler::InduceGluonEmission() {

  if(!p_dip) return false;

  if(p_dip->Status() && p_dip->PointerHandling()==0 &&
     p_dip->IsType()!=Dipole::incorrect);
  else return false;

  if(p_dix) { delete p_dix; p_dix=NULL;}
  if(p_ban) { delete p_ban; p_ban=NULL;}
  if(p_ati) { delete p_ati; p_ati=NULL;}
  if(p_glu) { delete p_glu; p_glu=NULL;}

  //No testing of global parameters.

  assert( p_dip->InvMass() > Sudakov_Calculator::MinOfK2t() );

  //assert(GenerateEfracs()); assert(TestEfracs());

  /*
  //OLD APPROACH.
  if( !GenerateEfracs() || !TestEfracs() ) {
    //cerr<<"\n"
    //<<"Method: const bool ADICIC::Dipole_Handler::InduceGluonEmission(): "
    //<<" Warning: Could not generate energy fractions!\n"<<endl;
    return false;
  }
  */

  //NEW APPROACH.
  if( p_sudakov->GenerateEfracsFor(*p_dip) ) {
    bool dummygsplit;
    p_sudakov->GetResult(dummygsplit,m_p2t,m_x1,m_x3);
  }
  else return false;

#ifdef DIPOLE_HANDLER_OUTPUT
  cout<<"\ttransverse momentum and energy fractions:\n\t\t p2t=";
  cout<<m_p2t<<endl;
  cout<<"\t\t x1="<<m_x1<<endl;
  cout<<"\t\t x3="<<m_x3<<endl;
#endif

  assert((p_dip->TotP())[0] > 0.0);

  m_p1=p_dip->GetTopBranchPointer()->Momentum();
  m_p3=p_dip->GetBotBranchPointer()->Momentum();
  assert(m_p1.Abs2() >= 0.0);
  assert(m_p3.Abs2() >= 0.0);
  assert(m_p1[0] > 0.0);
  assert(m_p3[0] > 0.0);

  assert(GenerateMomenta());
  assert(GenerateSplitting());

  return true;

}



//=============================================================================



const bool Dipole_Handler::InitCalcBox() {    //Static.

  static Calcpair qqpa, qgpa, gqpa, ggpa;

  if(Sudakov_Calculator::IsAlphaSRunning()==false) {

    //Arrange the Sudakov's.
    qqpa.p_sud=new Sudakov<Dipole::qqbar,Alpha_S_Fix>;
    qgpa.p_sud=new Sudakov<Dipole::qg,Alpha_S_Fix>;
    gqpa.p_sud=new Sudakov<Dipole::gqbar,Alpha_S_Fix>;
    ggpa.p_sud=new Sudakov<Dipole::gg,Alpha_S_Fix>;
    assert(qqpa.p_sud);
    assert(qgpa.p_sud);
    assert(gqpa.p_sud);
    assert(ggpa.p_sud);

    //Establish the overall recoil strategy right now and here.
    qqpa.p_rec=new Recoil<Kleiss_Strategy>;
    qgpa.p_rec=new Recoil<FixDir3_Strategy>;
    gqpa.p_rec=new Recoil<FixDir1_Strategy>;
    ggpa.p_rec=
      new Recoil<MinimizePt_Strategy>;
      //new Recoil<Lonnblad_Strategy>;
      //new Recoil<OldAdicic_Strategy>;
      //new Recoil<Test_Strategy>;
    assert(qqpa.p_rec);
    assert(qgpa.p_rec);
    assert(gqpa.p_rec);
    assert(ggpa.p_rec);

    //Fix the whole map - finishing arrangement of the calculator box.
    s_map[Dipole::qqbar] = &qqpa;
    s_map[Dipole::qg]    = &qgpa;
    s_map[Dipole::gqbar] = &gqpa;
    s_map[Dipole::gg]    = &ggpa;

    return true;

  }

  cerr<<"\nSorry :o( Option of using a running alpha_s ";
  cerr<<" has not been implemented yet.\n";
  assert(0);

}





const bool Dipole_Handler::GenerateMomenta() {

  const Vec4D& Plab=p_dip->TotP();
  //Preliminary approach - already done in EmitGluon():
  //m_p1=p_dip->GetTopBranchPointer()->Momentum();
  //m_p3=p_dip->GetBotBranchPointer()->Momentum();

#ifdef DIPOLE_HANDLER_OUTPUT
  cout<<"\tlab frame - before:\n\t\t P =";
  cout<<Plab<<"\t"<<p_dip->InvMass()<<"  "<<p_dip->Mass()<<endl;
  cout<<"\t\t q1="<<m_p1<<"\t "<<m_p1.Abs2()<<endl;
  cout<<"\t\t q3="<<m_p3<<"\t "<<m_p3.Abs2()<<endl;
#endif

  Recoil_Setup Iset;
  Iset.E2=sqrt(p_dip->InvMass());
  Iset.E1=0.5*Iset.E2*m_x1;
  Iset.E3=0.5*Iset.E2*m_x3;
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

  p_dip->GetTopBranchPointer()->SetMomentum(m_p1);
  p_dip->GetBotBranchPointer()->SetMomentum(m_p3);

  p_glu=new Dipole::Glubranch(m_p2); assert(p_glu);
  p_dix=new Dipole(*p_dip); assert(p_dix);

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
  p_dip->SetProdScale()=m_p2t;
  p_dip->SetEmitScale()=m_p2t;
  p_dix->SetProdScale()=m_p2t;
  p_dix->SetEmitScale()=m_p2t;

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



Dipole_Handler::Calcpair::~Calcpair() {
  if(p_sud) delete p_sud; if(p_rec) delete p_rec;
}



//=============================================================================



//=============================================================================





//eof
