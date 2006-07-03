//bof
//Version: 4 ADICIC++-0.0/2006/06/29

//Emissiontest.C - testing the first emission.



#include <typeinfo>
#include <cassert>
#include <iostream>
#include <ioextra>
#include <enumextra>
#include <mathextra>
#include <histoextra>
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
  Dipole_Flavour_Init dfi(true);

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
  Sudakov_Calculator::ShowEnvironment();
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

    Sudakov_Calculator::ShowEnvironment();
    Dipole_Handler::ShowCalcBox();
    Dipole_Parameter_Init dpi;

    cout<<"\nRunning?="<<Sudakov_Calculator::IsAlphaSRunning()<<endl;
    cout<<"NfFix="<<dpa.sud.NfFix()<<endl;
    cout<<"ASApp="<<Sudakov_Calculator::AlphaSApprox()<<endl;
    cout<<"ASCor="<<Sudakov_Calculator::AlphaSCorr(700.0)<<endl;
    cout<<"ASCor="<<Sudakov_Calculator::AlphaSCorr(8100.0)<<endl;
    cout<<"Nf(1.20)="<<Sudakov_Calculator::Nf(1.2)<<endl;
    cout<<"Nf(8100)="<<Sudakov_Calculator::Nf(8100.0)<<endl;

    cout<<"SudakovInit="<<Sudakov_Calculator::AdjustEnvironment()<<endl;

    cout<<"Running?="<<Sudakov_Calculator::IsAlphaSRunning()<<endl;
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

    Dipole_Handler::Carrier box;

    {

      unsigned trials=1;
      Dipole Dip(b1,a1);
      Dipole_Handler H(Dip);

      H.ShowSudakov(); H.ShowRecoil();
      cout<<Dip<<endl; Dip.PrintTowers();

      //At this stage factorization scale is correct.
      while(H.InduceDipoleRadiation()==false) ++trials;
      cout<<"Trials for a successful radiation="<<trials<<endl;
      assert(H.FinishDipoleRadiation());

      H.DecoupleNew(box);

      if(box.pDip) {
	assert(!box.pAqu && !box.pQua && box.pGlu);
	cout<<  "\t  p1="<<b1.Momentum()<<" \t "<<b1.Momentum().Abs2();
	cout<<"\n\t  p2="<<box.pGlu->Momentum()
	    <<" \t "<<box.pGlu->Momentum().Abs2();
	cout<<"\n\t  p3="<<a1.Momentum()<<" \t "<<a1.Momentum().Abs2();
	cout<<"\n\t sum="<<b1.Momentum()+box.pGlu->Momentum()+a1.Momentum();
	cout<<"\n\t       ";
	for(char i=0; i<4; ++i) cout<<b1.Momentum()[i]/pl[i]<<" \t ";
	cout<<endl<<"\t       ";
	for(char i=0; i<4; ++i) cout<<box.pGlu->Momentum()[i]/pl[i]<<" \t ";
	cout<<endl<<"\t       ";
	for(char i=0; i<4; ++i) cout<<a1.Momentum()[i]/pl[i]<<" \t ";
	cout<<endl;
	const Dipole& Din=*box.pDip;
	Dip.SetFactScale()=box.pDip->SetFactScale()=Dip.EmitScale();
	cout<<"Recoil is "<<box.RecoComp<<"; and, "
	    <<"is order given as olddip newdip? "<<box.DipOrder<<endl;
	if(box.DipOrder) {
	  cout<<Dip<<endl<<Din<<endl;
	  Dip.PrintTowers(); Din.PrintTowers();
	} else {
	  cout<<Din<<endl<<Dip<<endl;
	  Din.PrintTowers(); Dip.PrintTowers();
	}
      } else {
	assert(box.pAqu && box.pQua && box.pGlu);
	Dip.SetFactScale()=Dip.EmitScale();
	cout<<"Is top(0) or bottom(1) split? "<<box.DipOrder<<endl;
	if(box.DipOrder) {
	  cout<<  "\t  p1="<<b1.Momentum()<<" \t "<<b1.Momentum().Abs2();
	  cout<<"\n\t  p2="<<box.pAqu->Momentum()
	      <<" \t "<<box.pAqu->Momentum().Abs2();
	  cout<<"\n\t  p3="<<box.pQua->Momentum()
	      <<" \t "<<box.pQua->Momentum().Abs2();
	  cout<<"\n\t  pa="<<a1.Momentum()<<" \t "<<a1.Momentum().Abs2();
	  cout<<"\n\t  pg="<<box.pGlu->Momentum()
	      <<" \t "<<box.pGlu->Momentum().Abs2();
	  cout<<"\n\t sum="
	      <<b1.Momentum()+box.pAqu->Momentum()+box.pQua->Momentum();
	  cout<<"\n\t       ";
	  for(char i=0; i<4; ++i) cout<<b1.Momentum()[i]/pl[i]<<" \t ";
	  cout<<endl<<"\t       ";
	  for(char i=0; i<4; ++i) cout<<box.pAqu->Momentum()[i]/pl[i]<<" \t ";
	  cout<<endl<<"\t       ";
	  for(char i=0; i<4; ++i) cout<<box.pQua->Momentum()[i]/pl[i]<<" \t ";
	} else {
	  cout<<  "\t  p1="<<box.pAqu->Momentum()
	      <<" \t "<<box.pAqu->Momentum().Abs2();
	  cout<<"\n\t  p2="<<box.pQua->Momentum()
	      <<" \t "<<box.pQua->Momentum().Abs2();
	  cout<<"\n\t  p3="<<a1.Momentum()<<" \t "<<a1.Momentum().Abs2();
	  cout<<"\n\t  pb="<<b1.Momentum()<<" \t "<<b1.Momentum().Abs2();
	  cout<<"\n\t  pg="<<box.pGlu->Momentum()
	      <<" \t "<<box.pGlu->Momentum().Abs2();
	  cout<<"\n\t sum="
	      <<a1.Momentum()+box.pAqu->Momentum()+box.pQua->Momentum();
	  cout<<"\n\t       ";
	  for(char i=0; i<4; ++i) cout<<box.pAqu->Momentum()[i]/pl[i]<<" \t ";
	  cout<<endl<<"\t       ";
	  for(char i=0; i<4; ++i) cout<<box.pQua->Momentum()[i]/pl[i]<<" \t ";
	  cout<<endl<<"\t       ";
	  for(char i=0; i<4; ++i) cout<<a1.Momentum()[i]/pl[i]<<" \t ";
	}
	cout<<endl;
	cout<<"Recoil is "<<box.RecoComp<<endl;
	cout<<Dip<<endl;
	Dip.PrintTowers();
      }

      H.ShowSudakov(); H.ShowRecoil();

    }

    if(box.pDip) { delete box.pDip; if(box.pGlu) delete box.pGlu;}
    else { if(box.pQua) delete box.pQua; if(box.pAqu) delete box.pAqu;}

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

    //The Non STRICT_DIPOLE_HANDLER test.
    //It is helpful to use the DIPOLE_OUTPUT.

    Vec4D pl(45.0, 20.0,-5.0, 40.0);
    Vec4D pr(45.0,-20.0, 5.0,-40.0);
    //Dipole::Branch     b1(info.quark.u,pl);
    Dipole::Glubranch  b1(pl);
    //Dipole::Antibranch a1(info.antiq.u,pr);
    Dipole::Glubranch  a1(pr);

    Dipole Dip(b1,a1);
    {
      Dipole_Handler H(Dip);
      cout<<Dip<<endl<<endl;
      while(H.InduceDipoleRadiation()) {
	H.FinishDipoleRadiation();
	Dip.SetFactScale()=Dip.EmitScale();
	cout<<Dip<<endl;
      }
    }
    cout<<"ENDE: "<<Dip<<endl<<endl;
  }

  cout<<om::greenbg;
  cout<<"====================================================================";
  cout<<om::reset<<endl;
  cout<<endl; cin>>enter; cout<<endl;
  cout<<om::greenbg;
  cout<<"====================================================================";
  cout<<om::reset<<endl;

  {
    bool     xtest=false;
    unsigned total=4000000; //total=0;
    unsigned count=0;
    unsigned gluons=0;
    unsigned quarks=0, cd=0, cu=0, cs=0, cc=0, cb=0;

    Xhisto histo(40000);

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

      Dipole_Handler::Carrier box;

      {

	Dipole Dip(b1,a1);
	Dipole_Handler H(Dip);

#ifdef EMISSIONTEST_OUTPUT
	cout<<Dip<<endl<<endl;
#endif

	if( H.InduceDipoleRadiation(1,control,i==total) &&
	    H.FinishDipoleRadiation() );
	else {
	  ++count;
	  H.DecoupleNew(box);
	  assert(box.pQua==NULL);
	  assert(box.pAqu==NULL);
	  assert(box.pGlu==NULL);
	  if(control) {//||1
	    cout<<i<<"\t\t"<<count<<"\t\t"<<Dip.ProdScale()<<"\t\t";
	    cout<<box.pGlu; cout<<",";
	    cout<<box.pQua; cout<<",";
	    cout<<box.pAqu; cout<<endl;
	  }
	  continue;
	}
	H.DecoupleNew(box);
	kf::code kfc;
	if(box.pDip) {
	  assert(!box.pAqu && !box.pQua && box.pGlu);
	  Dip.SetFactScale()=box.pDip->SetFactScale()=Dip.EmitScale();
	  kfc=kf::gluon;
	  if(xtest) {
	    //x1x3 Test:
	    Multidouble m(2,0.0);
	    if(box.DipOrder) {
	      m[0]=Dip.GetTopBranchPointer()->Momentum()[0];
	      m[1]=box.pDip->GetBotBranchPointer()->Momentum()[0];
	    } else {
	      m[1]=Dip.GetBotBranchPointer()->Momentum()[0];
	      m[0]=box.pDip->GetTopBranchPointer()->Momentum()[0];
	    }
	    m[0]*=2/(pl[0]+pr[0]);
	    m[1]*=2/(pl[0]+pr[0]);
	    histo.Insert(m);
	  }
	} else {
	  assert(box.pAqu && box.pQua && box.pGlu);
	  Dip.SetFactScale()=Dip.EmitScale();
	  kfc=box.pQua->Flav().Kfcode();
	  assert(box.pAqu->Flav().Kfcode()==kfc);
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
	  if(box.pGlu) cout<<box.pGlu->Flav();
	  else cout<<box.pGlu; cout<<",";
	  if(box.pQua) cout<<box.pQua->Flav();
	  else cout<<box.pQua; cout<<",";
	  if(box.pAqu) cout<<box.pAqu->Flav();
	  else cout<<box.pAqu; cout<<"\t\t";
	  cout<<Flavour(kfc)<<","<<box.DipOrder<<","<<box.RecoComp<<endl;
	}

#ifdef EMISSIONTEST_OUTPUT
	if(box.pDip) {
	  const Dipole& Din=*box.pDip;
	  if(Dip.IsType()==Dipole::qg) cout<<Dip<<endl<<Din<<endl;
	  else cout<<Din<<endl<<Dip<<endl;
	}
#endif

      }

      if(box.pDip) { delete box.pDip; if(box.pGlu) delete box.pGlu;}
      else { if(box.pQua) delete box.pQua; if(box.pAqu) delete box.pAqu;}

    }

    //cout<<histo.Entries()<<endl;
    //histo.Reset();
    //Multidouble a(7,1.1);
    //Multidouble b(4,3.2);
    //histo.Insert(a);
    //histo.Insert(b);
    if(xtest) {
      string name("x1x3_emitest.dat");
      cout<<"Outputting "<<histo.Entries()
	  <<" Multidoubles to file "<<name<<".\n";
      histo.Output(name);
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
