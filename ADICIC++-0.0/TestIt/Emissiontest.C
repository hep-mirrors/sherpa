//bof
//Version: 1 ADICIC++-0.0/2004/06/08

//Emissiontest.C - testing the first emission.



#include <typeinfo>
#include <cassert>
#include <iostream>
#include <ioextra>
#include <enumextra>
#include "Message.H"
#include "Dipole.H"
#include "Dipole_Handler.H"
#include "Sudakov_Calculator.H"
#include "Recoil_Calculator.H"
#include "Chain.H"


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

  msg.SetModifiable(true);

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

  Dipole_Handler::ShowCalcBox();

  cout<<endl; cin>>enter; cout<<endl;
  cout<<"=============================================================="<<endl;

  {

    Vec4D pl(45.0, 20.0,-5.0, 40.0);
    Vec4D pr(45.0,-20.0, 5.0,-40.0);
    Dipole D1;
    Dipole::Branch b0(info.quark.s,pl);
    Dipole::Antibranch a0(info.antiquark.t,pl);
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
    Dipole::Antibranch a1(info.antiquark.u,pr);
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
      H1.ShowRecoilStrategy();
      Dipole_Handler H2(D2);
      H2.ShowSudakov();
      H2.ShowRecoilStrategy();

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
      H3.ShowRecoilStrategy();
      {
	Dipole E(b1,a1);
	E|H3;
	cout<<"H3 docking="<<H3.IsDocked()<<endl;
	H3.ShowSudakov();
	H3.ShowRecoilStrategy();
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

  cout<<"=============================================================="<<endl;
  cout<<endl;

  cout<<"============================="<<endl;
  cout<<" Testing the first emission. "<<endl;
  cout<<"============================="<<endl;

  {
    cout<<"Running?="<<Sudakov_Calculator::IsAlphaSRunning()<<endl;
    cout<<"MinScale="<<Sudakov_Calculator::MinOfK2t()<<endl;
    cout<<"MaxScale="<<Sudakov_Calculator::MaxOfK2t()<<endl;
    cout<<"ASFix="<<Sudakov_Calculator::AlphaSFix()<<endl;
    cout<<"ASApp="<<Sudakov_Calculator::AlphaSApprox()<<endl;
    cout<<"ASCor="<<Sudakov_Calculator::AlphaSCorr(700.0)<<endl;
    cout<<"ASCor="<<Sudakov_Calculator::AlphaSCorr(8100.0)<<endl;
    cout<<"SudakovInit="<<Sudakov_Calculator::Init(NULL)<<endl;
    cout<<"Running?="<<Sudakov_Calculator::IsAlphaSRunning()<<endl;
    cout<<"ASFix="<<Sudakov_Calculator::AlphaSFix()<<endl;
    cout<<"ASApp="<<Sudakov_Calculator::AlphaSApprox()<<endl;
    cout<<"ASCor="<<Sudakov_Calculator::AlphaSCorr(700.0)<<endl;
    cout<<"ASCor="<<Sudakov_Calculator::AlphaSCorr(8100.0)<<endl;

    TEMP::CPTEST=true;

    Vec4D pl(45.0, 20.0,-5.0, 40.0);
    Vec4D pr(45.0,-20.0, 5.0,-40.0);

    Dipole::Branch     b1(info.quark.u,pl);
    //Dipole::Glubranch  b1(pl);
    Dipole::Antibranch a1(info.antiquark.u,pr);
    //Dipole::Glubranch  a1(pr);

    Dipole* pDin=NULL;
    Dipole::Glubranch* pGlu=NULL;

    {

      bool below;
      Dipole Dip(b1,a1);
      Dipole_Handler H(Dip);

      H.ShowSudakov(); H.ShowRecoilStrategy();
      cout<<Dip<<endl; Dip.PrintTowers();

      assert(H.ManageGluonEmission());// && H.ManageGluonEmission());

      H.DecoupleNewDipole(pDin,below); assert(pDin);
      H.DecoupleGlubranch(pGlu); assert(pGlu);

      cout<<"\t  p1="<<b1.Momentum()<<" \t "<<b1.Momentum().Abs2()<<endl;
      cout<<"\t  p2="<<pGlu->Momentum()<<" \t "<<pGlu->Momentum().Abs2()<<endl;
      cout<<"\t  p3="<<a1.Momentum()<<" \t "<<a1.Momentum().Abs2()<<endl;
      cout<<"\t sum="<<b1.Momentum()+pGlu->Momentum()+a1.Momentum()<<endl;

      cout<<"\t       ";
      for(char i=0; i<4; ++i) cout<<b1.Momentum()[i]/pl[i]<<" \t ";
      cout<<endl<<"\t       ";
      for(char i=0; i<4; ++i) cout<<pGlu->Momentum()[i]/pl[i]<<" \t ";
      cout<<endl<<"\t       ";
      for(char i=0; i<4; ++i) cout<<a1.Momentum()[i]/pl[i]<<" \t ";
      cout<<endl;

      const Dipole& Din=*pDin;

      cout<<"Is order given as olddip newdip? "<<below<<endl;

      if(below) {
	cout<<Dip<<endl<<Din<<endl;
	Dip.PrintTowers(); Din.PrintTowers();
      } else {
	cout<<Din<<endl<<Dip<<endl;
	Din.PrintTowers(); Dip.PrintTowers();
      }
      H.ShowSudakov(); H.ShowRecoilStrategy();

    }

    delete pDin;
    delete pGlu;

    TEMP::CPTEST=false;
  }

  cout<<"=============================================================="<<endl;
  cout<<"Test: 8.123^7="<<power<7>(8.123)<<endl;
  cout<<"Test: 4.573^3="<<power<3>(4.573)<<endl;
  cout<<"Test: 7.000^2="<<power<2>(7.000);
  Dipole_Handler::ShowCalcBox();
  cout<<"Number of Dipole_Handler's = "<<Dipole_Handler::InStore<<endl;
  cout<<"=============================================================="<<endl;

  {
    unsigned total=4000000; //total=0;
    unsigned count=0;

    for(unsigned i=1; i<=total; ++i) {

      //cout<<endl; cin>>enter; cout<<endl;
      //cout<<"=>"<<i<<endl;
      bool control=!(i%50000);

      Vec4D pl(45.0, 20.0,-5.0, 40.0);
      Vec4D pr(45.0,-20.0, 5.0,-40.0);
      Dipole::Branch     b1(info.quark.u,pl);
      //Dipole::Glubranch  b1(pl);
      Dipole::Antibranch a1(info.antiquark.u,pr);
      //Dipole::Glubranch  a1(pr);

      Dipole* pDin=NULL;
      Dipole::Glubranch* pGlu=NULL;

      {

	bool dummy;
	Dipole Dip(b1,a1);
	Dipole_Handler H(Dip);

#ifdef EMISSIONTEST_OUTPUT
	cout<<Dip<<endl;
#endif

	if( H.InduceGluonEmission() && H.FinishGluonEmission() );
	else {
	  ++count;
	  if(control) cout<<i<<"\t\t"<<count<<"\t\t"<<Dip.ProdScale()<<endl;
	  continue;
	}

	if(control) cout<<i<<"\t\t"<<count<<"\t\t"<<Dip.ProdScale()<<endl;

	H.DecoupleNewDipole(pDin,dummy); assert(pDin);
	H.DecoupleGlubranch(pGlu); assert(pGlu);
	const Dipole& Din=*pDin;

#ifdef EMISSIONTEST_OUTPUT
	if(Dip.IsType()==Dipole::qg) cout<<Dip<<endl<<Din<<endl;
	else cout<<Din<<endl<<Dip<<endl;
#endif

      }

      delete pDin;
      delete pGlu;

    }

    cout<<endl<<"Total number of failures="<<count;
    cout<<"   ("<<1.0*count/total<<")."<<endl;
  }

  cout<<"=============================================================="<<endl;

  cout<<endl;

}





//eof
