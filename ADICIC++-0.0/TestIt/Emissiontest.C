//bof
//Version: 2 ADICIC++-0.0/2004/09/01

//Emissiontest.C - testing the first emission.



#include <typeinfo>
#include <cassert>
#include <iostream>
#include <ioextra>
#include <enumextra>
#include <mathextra>
#include "Message.H"
#include "Run_Parameter.H"
#include "Dipole.H"
#include "Dipole_Handler.H"
#include "Sudakov_Calculator.H"
#include "Recoil_Calculator.H"
#include "Chain.H"
#include "Dipole_Parameter.H"
#include "Paraminit.H"


#define EMISSIONTEST_OUTPUT EMISSIONTEST_OUTPUT
#undef  EMISSIONTEST_OUTPUT





using namespace std;
using namespace ATOOLS;
using namespace ADICIC;





class Lev1;

class Lev2 {
  friend class Lev1;
  class {
    friend class Lev2;
    int t;
  public:
    int& T() { t=333; return t;}
    int& Tt() { return t;}
  } lev3;
public:
  Lev2() { lev3.t=77;}
};

class Lev1 {
  Lev2 member;
public:
  void Print() { cout<<"In-class trick works: "<<member.lev3.T()<<endl;}
  void Vary() { member.lev3.T()=43; cout<<"Varying: "<<member.lev3.Tt()<<endl;}
  //void Test() { member.lev3.t=99;}    //does not work - as wished
};





int main() {

  //for(int q=0; q<=30000000; q++) factorial(56);
  //for(int q=0; q<=30000000; q++) Factorial<56>::RET;
  cout<<"RuntimeTest1: 53! = "
      <<factorial(53)<<" | "<<Factorial<53>::RET<<endl;

  //for(int q=0; q<=60000000; q++) power(5.67,22);
  //for(int q=0; q<=60000000; q++) power<22>(5.67);

  cout<<"RuntimeTest2: 5.67^22 = "
      <<power(5.67,22)<<" | "<<power<22>(5.67)<<endl;

  //abort();

  msg.SetModifiable(true);
  cout<<Dipole_Flavour_Init::Status()<<endl;
  cout<<Dipole_Flavour_Init::DoIt(true)<<endl;
  cout<<Dipole_Flavour_Init::DoIt()<<endl;
  cout<<Dipole_Flavour_Init::Status()<<endl;

  cout<<endl;

  Lev1 kkk;
  kkk.Print();
  kkk.Vary();

  Lev2 uuu;
  //uuu.lev3.t=11;    //does not work - as wished
  //uuu.lev3.T()=16;    //does not work - as wished

  cout<<endl;

  Chain cha;
  cha.Print();
  cout<<cha.MaxParticleNumber()<<endl;
  cout<<cha.MaxDipoleNumber()<<endl;

  cout<<endl;
  cout<<"=============================================================="<<endl;

  Dipole_Parameter::Show();
  Sudakov_Calculator::ShowParameters();
  Dipole_Handler::ShowCalcBox();

  cout<<endl; cin>>enter; cout<<endl;
  cout<<"=============================================================="<<endl;

  {

    Vec4D pl(45.0, 20.0,-5.0, 40.0);
    Vec4D pr(45.0,-20.0, 5.0,-40.0);
    Dipole D1;
    Dipole::Branch b0(info.quark.s,pl);
    Dipole::Antibranch a0(info.antiq.t,pl);
    Dipole D2(b0,a0,33);
    cout<<D1<<endl<<D2<<endl;
    D1.PrintTowers(); D2.PrintTowers();
    D1=D2;
    cout<<endl;

  }
  cout<<"=============================================================="<<endl;
  cout<<endl;

  cout<<"======================================="<<endl;
  cout<<" Testing the Dipole_Handler structure. "<<endl;
  cout<<"======================================="<<endl;

  {

    Vec4D pl(45.0, 20.0,-5.0, 40.0);
    Vec4D pr(45.0,-20.0, 5.0,-40.0);

    Dipole D1;
    Dipole::Branch b1(info.quark.u,pl);
    Dipole::Antibranch a1(info.antiq.u,pr);
    Dipole D2(b1,a1);
    cout<<D1<<endl<<D2<<endl;
    D1.PrintTowers(); D2.PrintTowers();
    cout<<"D1 handling="<<D1.IsHandled()<<endl;
    cout<<"D2 handling="<<D2.IsHandled()<<endl;
    D1|0; D2|0;
    cout<<"D1 handling="<<D1.IsHandled()<<endl;
    cout<<"D2 handling="<<D2.IsHandled()<<endl;

    {
      Dipole_Handler H1;
      H1.ShowSudakov();
      H1.ShowRecoil();
      Dipole_Handler H2(D2);
      H2.ShowSudakov();
      H2.ShowRecoil();

      cout<<"H1 docking="<<H1.IsDocked()<<endl;
      cout<<"H2 docking="<<H2.IsDocked()<<endl;
      cout<<"H2 docking D1? "<<H2.IsDockedAt(D1)<<endl;
      cout<<"H2 docking D2? "<<H2.IsDockedAt(D2)<<endl;
      cout<<"D2 handled by H1? "<<D2.IsHandledBy(H1)<<endl;
      cout<<"D2 handled by H2? "<<D2.IsHandledBy(H2)<<endl;
      cout<<endl;
      cout<<H1.DetachDipole(&D1)<<endl;
      cout<<H1.DetachDipole(&D2)<<endl;
      cout<<H2.DetachDipole(&D1)<<endl;
      cout<<H2.DetachDipole(&D2)<<endl;
      cout<<H1.AttachDipole(&D1)<<endl;
      cout<<H1.AttachDipole(&D2)<<endl;
      cout<<H2.AttachDipole(&D1)<<endl;
      cout<<H2.AttachDipole(&D2)<<endl;
      cout<<endl;

      Dipole_Handler H3;
      H3.ShowSudakov();
      H3.ShowRecoil();
      {
	Dipole E(b1,a1);
	E|H3;
	cout<<"H3 docking="<<H3.IsDocked()<<endl;
	H3.ShowSudakov();
	H3.ShowRecoil();
	cout<<E<<endl;
      }
      cout<<"H3 docking="<<H3.IsDocked()<<endl;
      cout<<endl;

      cout<<(D1|H3)<<endl;
      cout<<H3.DetachDipole(&D1)<<endl;
      cout<<H3.DetachDipole(&D2)<<endl;
      cout<<H3.AttachDipole(&D1)<<endl;
      cout<<H3.AttachDipole(&D2)<<endl;
      cout<<(D1|H1)<<endl;
      cout<<(D1|H2)<<endl;
      cout<<(D1|H3)<<endl;
      cout<<(D2|H1)<<endl;
      cout<<(D2|H2)<<endl;
      cout<<(D2|H3)<<endl;
      cout<<endl;
      D2|0;
      cout<<(D2|H3)<<endl;
      cout<<(D2|H1)<<endl;
      cout<<(D2|H2)<<endl;
      cout<<D1<<endl<<D2<<endl;
      cout<<"D1 handled by H1? "<<D1.IsHandledBy(H1)<<endl;
      cout<<"D1 handled by H2? "<<D1.IsHandledBy(H2)<<endl;
      cout<<"D1 handled by H3? "<<D1.IsHandledBy(H3)<<endl;
      cout<<"D2 handled by H1? "<<D2.IsHandledBy(H1)<<endl;
      cout<<"D2 handled by H2? "<<D2.IsHandledBy(H2)<<endl;
      cout<<"D2 handled by H3? "<<D2.IsHandledBy(H3)<<endl;

      Dipole_Handler::ShowCalcBox();
      cout<<"Number of Dipole_Handler's = "<<Dipole_Handler::InStore<<endl;
    }

    cout<<"D1 handling="<<D1.IsHandled()<<endl;
    cout<<"D2 handling="<<D2.IsHandled()<<endl;
    cout<<"Number of Dipole_Handler's = "<<Dipole_Handler::InStore<<endl;

  }

  cout<<om::greenbg;
  cout<<"====================================================================";
  cout<<om::reset<<endl;
  cout<<endl; cin>>enter; cout<<endl;
  cout<<om::greenbg;
  cout<<"====================================================================";
  cout<<om::reset<<endl;

  cout<<"============================="<<endl;
  cout<<" Testing the first emission. "<<endl;
  cout<<"============================="<<endl;

  {

    extern Run_Parameter ATOOLS::rpa;
    rpa.gen.SetEcms(90.0);

    Sudakov_Calculator::ShowParameters();
    Dipole_Handler::ShowCalcBox();
    cout<<"\nDip Param Init?: "<<Dipole_Parameter_Init::Status()<<endl;
    cout<<"\nDo Init: "<<Dipole_Parameter_Init::DoIt()<<endl;////////////////
    cout<<"Do Init: "<<Dipole_Parameter_Init::DoIt()<<endl;//////////////////
    cout<<"Dip Param Init?: "<<Dipole_Parameter_Init::Status()<<endl;

    cout<<"\nRunning?="<<Sudakov_Calculator::IsAlphaSRunning()<<endl;
    cout<<"MinScale="<<Sudakov_Calculator::MinOfK2t()<<endl;
    cout<<"MaxScale="<<Sudakov_Calculator::MaxOfK2t()<<endl;
    cout<<"NfFix="<<Sudakov_Calculator::NfFix()<<endl;
    cout<<"ASFix="<<Sudakov_Calculator::AlphaSFix()<<endl;
    cout<<"ASApp="<<Sudakov_Calculator::AlphaSApprox()<<endl;
    cout<<"ASCor="<<Sudakov_Calculator::AlphaSCorr(700.0)<<endl;
    cout<<"ASCor="<<Sudakov_Calculator::AlphaSCorr(8100.0)<<endl;
    cout<<"Nf(1.20)="<<Sudakov_Calculator::Nf(1.2)<<endl;
    cout<<"Nf(8100)="<<Sudakov_Calculator::Nf(8100.0)<<endl;

    cout<<"SudakovInit="<<Sudakov_Calculator::Init(NULL)<<endl;////////////////

    cout<<"Running?="<<Sudakov_Calculator::IsAlphaSRunning()<<endl;
    cout<<"ASFix="<<Sudakov_Calculator::AlphaSFix()<<endl;
    cout<<"ASApp="<<Sudakov_Calculator::AlphaSApprox()<<endl;
    cout<<"ASCor="<<Sudakov_Calculator::AlphaSCorr(700.0)<<endl;
    cout<<"ASCor="<<Sudakov_Calculator::AlphaSCorr(8100.0)<<endl;
    cout<<"Nf(1.20)="<<Sudakov_Calculator::Nf(1.2)<<endl;
    cout<<"Nf(8100)="<<Sudakov_Calculator::Nf(8100.0)<<endl;

    TEMP::CPTEST=true;

    Vec4D pl(45.0, 20.0,-5.0, 40.0);
    Vec4D pr(45.0,-20.0, 5.0,-40.0);

    //Dipole::Branch     b1(info.quark.u,pl);
    Dipole::Glubranch  b1(pl);
    Dipole::Antibranch a1(info.antiq.u,pr);
    //Dipole::Glubranch  a1(pr);

    Dipole* pDin=NULL;
    Dipole::Glubranch* pGlu=NULL;
    Dipole::Antibranch* pAti=NULL;
    Dipole::Branch* pBan=NULL;

    {

      unsigned trials=1;
      bool below;
      Trio recoil;
      Dipole Dip(b1,a1);
      Dipole_Handler H(Dip);

      H.ShowSudakov(); H.ShowRecoil();
      cout<<Dip<<endl; Dip.PrintTowers();

      //assert(H.ManageGluonEmission());// && H.ManageGluonEmission());
      while(H.ManageDipoleRadiation()==false) ++trials;
      cout<<"Trials for a successful radiation="<<trials<<endl;

      H.DecoupleNew(pDin,pGlu,pAti,pBan,below,recoil);

      if(pDin) {
	assert(!pAti && !pBan && pGlu);
	cout<<  "\t  p1="<<b1.Momentum()<<" \t "<<b1.Momentum().Abs2();
	cout<<"\n\t  p2="<<pGlu->Momentum()<<" \t "<<pGlu->Momentum().Abs2();
	cout<<"\n\t  p3="<<a1.Momentum()<<" \t "<<a1.Momentum().Abs2();
	cout<<"\n\t sum="<<b1.Momentum()+pGlu->Momentum()+a1.Momentum();
	cout<<"\n\t       ";
	for(char i=0; i<4; ++i) cout<<b1.Momentum()[i]/pl[i]<<" \t ";
	cout<<endl<<"\t       ";
	for(char i=0; i<4; ++i) cout<<pGlu->Momentum()[i]/pl[i]<<" \t ";
	cout<<endl<<"\t       ";
	for(char i=0; i<4; ++i) cout<<a1.Momentum()[i]/pl[i]<<" \t ";
	cout<<endl;
	const Dipole& Din=*pDin;
	cout<<"Recoil is "<<recoil
	    <<"; and, is order given as olddip newdip? "<<below<<endl;
	if(below) {
	  cout<<Dip<<endl<<Din<<endl;
	  Dip.PrintTowers(); Din.PrintTowers();
	} else {
	  cout<<Din<<endl<<Dip<<endl;
	  Din.PrintTowers(); Dip.PrintTowers();
	}
      } else {
	assert(pAti && pBan && pGlu);
	cout<<"Is top(0) or bottom(1) split? "<<below<<endl;
	if(below) {
	  cout<<  "\t  p1="<<b1.Momentum()<<" \t "<<b1.Momentum().Abs2();
	  cout<<"\n\t  p2="<<pAti->Momentum()<<" \t "<<pAti->Momentum().Abs2();
	  cout<<"\n\t  p3="<<pBan->Momentum()<<" \t "<<pBan->Momentum().Abs2();
	  cout<<"\n\t  pa="<<a1.Momentum()<<" \t "<<a1.Momentum().Abs2();
	  cout<<"\n\t  pg="<<pGlu->Momentum()<<" \t "<<pGlu->Momentum().Abs2();
	  cout<<"\n\t sum="<<b1.Momentum()+pAti->Momentum()+pBan->Momentum();
	  cout<<"\n\t       ";
	  for(char i=0; i<4; ++i) cout<<b1.Momentum()[i]/pl[i]<<" \t ";
	  cout<<endl<<"\t       ";
	  for(char i=0; i<4; ++i) cout<<pAti->Momentum()[i]/pl[i]<<" \t ";
	  cout<<endl<<"\t       ";
	  for(char i=0; i<4; ++i) cout<<pBan->Momentum()[i]/pl[i]<<" \t ";
	} else {
	  cout<<  "\t  p1="<<pAti->Momentum()<<" \t "<<pAti->Momentum().Abs2();
	  cout<<"\n\t  p2="<<pBan->Momentum()<<" \t "<<pBan->Momentum().Abs2();
	  cout<<"\n\t  p3="<<a1.Momentum()<<" \t "<<a1.Momentum().Abs2();
	  cout<<"\n\t  pb="<<b1.Momentum()<<" \t "<<b1.Momentum().Abs2();
	  cout<<"\n\t  pg="<<pGlu->Momentum()<<" \t "<<pGlu->Momentum().Abs2();
	  cout<<"\n\t sum="<<a1.Momentum()+pAti->Momentum()+pBan->Momentum();
	  cout<<"\n\t       ";
	  for(char i=0; i<4; ++i) cout<<pAti->Momentum()[i]/pl[i]<<" \t ";
	  cout<<endl<<"\t       ";
	  for(char i=0; i<4; ++i) cout<<pBan->Momentum()[i]/pl[i]<<" \t ";
	  cout<<endl<<"\t       ";
	  for(char i=0; i<4; ++i) cout<<a1.Momentum()[i]/pl[i]<<" \t ";
	}
	cout<<endl;
	cout<<"Recoil is "<<recoil<<endl;
	cout<<Dip<<endl;
	Dip.PrintTowers();
      }

      H.ShowSudakov(); H.ShowRecoil();

    }

    if(pDin) { delete pDin; if(pGlu) delete pGlu;}
    else { if(pBan) delete pBan; if(pAti) delete pAti;}

    TEMP::CPTEST=false;
  }

  cout<<"=============================================================="<<endl;
  cout<<"Test: 5^8="<<power<8>(5)<<endl;
  cout<<"Test: 8.123^7="<<power<7>(8.123)<<endl;
  cout<<"Test: 4.573^3="<<power<3>(4.573)<<endl;
  cout<<"Test: 7.000^2="<<power<2>(7.000);
  Dipole_Handler::ShowCalcBox();
  cout<<"Number of Dipole_Handler's = "<<Dipole_Handler::InStore<<endl;
  cout<<"=============================================================="<<endl;
  cout<<om::greenbg;
  cout<<"====================================================================";
  cout<<om::reset<<endl;
  cout<<endl; cin>>enter; cout<<endl;
  cout<<om::greenbg;
  cout<<"====================================================================";
  cout<<om::reset<<endl;

  {
    unsigned total=4000000; //total=0;
    unsigned count=0;
    unsigned gluons=0;
    unsigned quarks=0, cd=0, cu=0, cs=0, cc=0, cb=0;

    for(unsigned i=1; i<=total; ++i) {

      //cout<<endl; cin>>enter; cout<<endl;
      //cout<<"=>"<<i<<endl;
      bool control=!(i%50000);

      Vec4D pl(45.0, 20.0,-5.0, 40.0);
      Vec4D pr(45.0,-20.0, 5.0,-40.0);
      //Dipole::Branch     b1(info.quark.u,pl);
      Dipole::Glubranch  b1(pl);
      //Dipole::Antibranch a1(info.antiq.u,pr);
      Dipole::Glubranch  a1(pr);

      Dipole* pDin=NULL;
      Dipole::Branch*     pBan=NULL;
      Dipole::Antibranch* pAti=NULL;
      Dipole::Glubranch*  pGlu=NULL;

      {

	bool below;
	Trio recoil;
	Dipole Dip(b1,a1);
	Dipole_Handler H(Dip);

#ifdef EMISSIONTEST_OUTPUT
	cout<<Dip<<endl;
#endif

	if( H.InduceDipoleRadiation() && H.FinishDipoleRadiation() );
	else {
	  ++count;
	  H.DecoupleNew(pDin,pGlu,pAti,pBan,below,recoil);
	  assert(pBan==NULL); assert(pAti==NULL); assert(pGlu==NULL);
	  if(control) {//||1
	    cout<<i<<"\t\t"<<count<<"\t\t"<<Dip.ProdScale()<<"\t\t";
	    cout<<pGlu; cout<<",";
	    cout<<pBan; cout<<",";
	    cout<<pAti; cout<<endl;
	  }
	  continue;
	}
	H.DecoupleNew(pDin,pGlu,pAti,pBan,below,recoil);
	kf::code kfc;
	if(pDin) {
	  assert(!pAti && !pBan && pGlu);
	  kfc=kf::gluon;
	} else {
	  assert(pAti && pBan && pGlu);
	  kfc=pBan->Flav().Kfcode();
	  assert(pAti->Flav().Kfcode()==kfc);
	}
	switch(kfc) {
	case kf::gluon: ++gluons; break;
	case kf::d    : ++cd; break;
	case kf::u    : ++cu; break;
	case kf::s    : ++cs; break;
	case kf::c    : ++cc; break;
	case kf::b    : ++cb; break;
	default       : assert(0);
	}

	if(control) {//||1
	  cout<<i<<"\t\t"<<count<<"\t\t"<<Dip.ProdScale()<<"\t\t";
	  if(pGlu) cout<<pGlu->Flav(); else cout<<pGlu; cout<<",";
	  if(pBan) cout<<pBan->Flav(); else cout<<pBan; cout<<",";
	  if(pAti) cout<<pAti->Flav(); else cout<<pAti; cout<<"\t\t";
	  cout<<Flavour(kfc)<<","<<below<<","<<recoil<<endl;
	}

#ifdef EMISSIONTEST_OUTPUT
	if(pDin) {
	  const Dipole& Din=*pDin;
	  if(Dip.IsType()==Dipole::qg) cout<<Dip<<endl<<Din<<endl;
	  else cout<<Din<<endl<<Dip<<endl;
	}
#endif

      }

      if(pDin) { delete pDin; if(pGlu) delete pGlu;}
      else { if(pBan) delete pBan; if(pAti) delete pAti;}

    }

    quarks=cd+cu+cs+cc+cb;

    cout<<endl;
    cout<<"Total number of   gluons="<<gluons;
    cout<<"   ("<<1.0*gluons/total<<")."<<endl;
    cout<<"Total number of   quarks="<<quarks;
    cout<<"   ("<<1.0*quarks/total<<")."<<endl;
    cout<<"                       d="<<cd;
    cout<<"   ("<<1.0*cd/total<<")."<<endl;
    cout<<"                       u="<<cu;
    cout<<"   ("<<1.0*cu/total<<")."<<endl;
    cout<<"                       s="<<cs;
    cout<<"   ("<<1.0*cs/total<<")."<<endl;
    cout<<"                       c="<<cc;
    cout<<"   ("<<1.0*cc/total<<")."<<endl;
    cout<<"                       b="<<cb;
    cout<<"   ("<<1.0*cb/total<<")."<<endl;
    cout<<"Total number of failures="<<count;
    cout<<"   ("<<1.0*count/total<<")."<<endl;
    cout<<"Total number of   events="<<count+quarks+gluons;
    cout<<"   ("<<1.0*(count+quarks+gluons)/total<<")."<<endl;

  }

  cout<<om::greenbg;
  cout<<"====================================================================";
  cout<<om::reset<<endl;

  cout<<endl;

}





//eof
