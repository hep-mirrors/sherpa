//bof
//Version: 1 ADICIC++-0.0/2004/06/03

//Chaintest.C - testing the first chaining.



#include <typeinfo>
#include <cassert>
#include <iostream>
#include <ioextra>
#include <enumextra>
#include "Dipole.H"
#include "Dipole_Handler.H"
#include "Sudakov_Calculator.H"
#include "Chain.H"
#include "Chain_Handler.H"


#define CHAINTEST_OUTPUT CHAINTEST_OUTPUT
#undef  CHAINTEST_OUTPUT





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

  cout<<endl;

  Lev1 kkk;
  kkk.Print();
  kkk.Vary();

  Lev2 uuu;
  //uuu.lev3.t=11;    //does not work - as wished
  //uuu.lev3.T()=16;    //does not work - as wished

  cout<<endl;

  Dipole_Handler::ShowCalcBox();

  cout<<endl; cin>>enter; cout<<endl;
  cout<<"=============================================================="<<endl;
  cout<<endl;

  cout<<"======================================"<<endl;
  cout<<" Testing the chaining - Preparations. "<<endl;
  cout<<"======================================"<<endl;

  {

    Chain cha;
    cha.Print();
    cout<<cha.MaxParticleNumber()<<endl;
    cout<<cha.MaxDipoleNumber()<<endl;

    Vec4D pl(45.0, 20.0,-5.0, 40.0);
    Vec4D pr(45.0,-20.0, 5.0,-40.0);
    Dipole::Branch     b1(info.quark.u,pl);
    Dipole::Glubranch  g1(pl);
    Dipole::Antibranch a1(info.antiquark.u,pr);
    Dipole::Glubranch  g2(pr);

    {
      Chain ch1(b1,a1,Chain::Initiator::simple_epem);
      ch1.Print();
      cout<<ch1.MaxParticleNumber()<<endl;
      cout<<ch1.MaxDipoleNumber()<<endl;

      Chain ch2(g1,g1,Chain::Initiator::simple_epem);
      ch2.Print();
      cout<<ch2.MaxParticleNumber()<<endl;
      cout<<ch2.MaxDipoleNumber()<<endl;

      Chain ch3(g1,g2,Chain::Initiator::simple_epem);
      ch3.Print();
      cout<<ch3.MaxParticleNumber()<<endl;
      cout<<ch3.MaxDipoleNumber()<<endl;

      Chain ch4(ch3);//(ch1)
      ch4.Print();
      cout<<ch4.MaxParticleNumber()<<endl;
      cout<<ch4.MaxDipoleNumber()<<endl;

      ch2.Clear();
      ch3.Clear();
      ch2=ch4;
      ch3.Initialize(b1,a1);
      ch2.Print();
      ch3.Print();

      cout<<"Number of Dipole_Particle's = "<<Dipole_Particle::InStore<<endl;
      cout<<"Number of Dipole's          = "<<Dipole::InStore<<endl;
      cout<<"Number of Chain's           = "<<Chain::InStore<<endl;
    }

    cout<<"============================================================"<<endl;

    cout<<"Number of Dipole_Particle's = "<<Dipole_Particle::InStore<<endl;
    cout<<"Number of Dipole's          = "<<Dipole::InStore<<endl;
    cout<<"Number of Chain's           = "<<Chain::InStore<<endl;

    cout<<endl; cin>>enter; cout<<endl;

    Chain ch1(b1,a1,Chain::Initiator::simple_epem);
    cout<<ch1<<endl;
    cout<<"cha handling="<<cha.IsHandled()<<endl;
    cout<<"ch1 handling="<<ch1.IsHandled()<<endl;
    cha|0; ch1|0;
    cout<<"cha handling="<<cha.IsHandled()<<endl;
    cout<<"ch1 handling="<<ch1.IsHandled()<<endl;

    {
      Chain_Handler H1;
      Chain_Handler H2(ch1);

      cout<<"H1 docking="<<H1.IsDocked()<<endl;
      cout<<"H2 docking="<<H2.IsDocked()<<endl;
      cout<<"H2 docking cha? "<<H2.IsDockedAt(cha)<<endl;
      cout<<"H2 docking ch1? "<<H2.IsDockedAt(ch1)<<endl;
      cout<<"ch1 handled by H1? "<<ch1.IsHandledBy(H1)<<endl;
      cout<<"ch1 handled by H2? "<<ch1.IsHandledBy(H2)<<endl;
      cout<<endl;
      cout<<H1.DetachChain(&cha)<<endl;
      cout<<H1.DetachChain(&ch1)<<endl;
      cout<<H2.DetachChain(&cha)<<endl;
      cout<<H2.DetachChain(&ch1)<<endl;
      cout<<H1.AttachChain(&cha)<<endl;
      cout<<H1.AttachChain(&ch1)<<endl;
      cout<<H2.AttachChain(&cha)<<endl;
      cout<<H2.AttachChain(&ch1)<<endl;
      cout<<endl;

      Chain_Handler H3;
      {
	Chain E(g1,g2,Chain::Initiator::simple_epem);
	E|H3;
	cout<<"H3 docking="<<H3.IsDocked()<<endl;
	cout<<E.IsRing()<<endl;
	cout<<E<<endl;
      }
      cout<<"H3 docking="<<H3.IsDocked()<<endl;
      cout<<endl;

      cout<<(cha|H3)<<endl;
      cout<<H3.DetachChain(&cha)<<endl;
      cout<<H3.DetachChain(&ch1)<<endl;
      cout<<H3.AttachChain(&cha)<<endl;
      cout<<H3.AttachChain(&ch1)<<endl;
      cout<<(cha|H1)<<endl;
      cout<<(cha|H2)<<endl;
      cout<<(cha|H3)<<endl;
      cout<<(ch1|H1)<<endl;
      cout<<(ch1|H2)<<endl;
      cout<<(ch1|H3)<<endl;
      cout<<endl;
      ch1|0;
      cout<<(ch1|H3)<<endl;
      cout<<(ch1|H1)<<endl;
      cout<<(ch1|H2)<<endl;
      cout<<cha<<endl<<ch1<<endl;
      cout<<"cha handled by H1? "<<cha.IsHandledBy(H1)<<endl;
      cout<<"cha handled by H2? "<<cha.IsHandledBy(H2)<<endl;
      cout<<"cha handled by H3? "<<cha.IsHandledBy(H3)<<endl;
      cout<<"ch1 handled by H1? "<<ch1.IsHandledBy(H1)<<endl;
      cout<<"ch1 handled by H2? "<<ch1.IsHandledBy(H2)<<endl;
      cout<<"ch1 handled by H3? "<<ch1.IsHandledBy(H3)<<endl;

      cout<<cha.IsEmpty()<<endl;
      cout<<ch1.IsEmpty()<<endl;
      cout<<cha.IsRing()<<endl;
      cout<<ch1.IsRing()<<endl;
      cout<<cha.ChainType()<<endl;
      cout<<ch1.ChainType()<<endl;
      cout<<ch1.ChainRoot()<<endl;

      cout<<H1.CompScale()<<endl;
      cout<<H2.CompScale()<<endl;
      cout<<H3.CompScale()<<endl;
      H1.Reset();
      H2.Reset();
      H3.Reset();

      cout<<"Number of Chain_Handler's = "<<Chain_Handler::InStore<<endl;
    }

    cout<<"cha handling="<<cha.IsHandled()<<endl;
    cout<<"ch1 handling="<<ch1.IsHandled()<<endl;
    cout<<"Number of Chain_Handler's = "<<Chain_Handler::InStore<<endl;

  }

  cout<<"=============================================================="<<endl;
  cout<<endl; cin>>enter; cout<<endl;
  cout<<"=============================================================="<<endl;
  cout<<endl;

  cout<<"=================================="<<endl;
  cout<<" Testing the chaining stepwisely. "<<endl;
  cout<<"=================================="<<endl;

  {

    Vec4D pl(45.0, 20.0,-5.0, 40.0);
    Vec4D pr(45.0,-20.0, 5.0,-40.0);
    Dipole::Branch     b1(info.quark.u,pl);
    Dipole::Glubranch  g1(pl);
    Dipole::Antibranch a1(info.antiquark.u,pr);
    Dipole::Glubranch  g2(pr);

    Chain cha(b1,a1,Chain::Initiator::simple_epem);
    Chain_Handler H(cha);

    cha.Print();
    cout<<cha.MaxParticleNumber()<<endl;
    cout<<cha.MaxDipoleNumber()<<endl;
    cout<<endl;
    cout<<"\e[7m\e[31m                    \e[0m";
    cout<<"\e[1m\e[40m\e[33mSTART\e[0m";
    cout<<"\e[7m\e[31m                    \e[0m"<<endl;

    while(H.EvolveChainByOneStep()) {
      cout<<cha<<endl;
      Vec4D test;
      cha.CheckMomentumConservation(test);
      cout<<test<<endl<<endl;
    }

    cout<<"\e[7m\e[31m                    \e[0m";
    cout<<"\e[1m\e[40m\e[33mSTOP\e[0m";
    cout<<"\e[7m\e[31m                    \e[0m"<<endl;
    cha.Print();
    cout<<cha.MaxParticleNumber()<<endl;
    cout<<cha.MaxDipoleNumber()<<endl;

    /*
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
    */

  }

  cout<<"=============================================================="<<endl;
  cout<<"Test: 8.123^7="<<power<7>(8.123)<<endl;
  cout<<"Test: 4.573^3="<<power<3>(4.573)<<endl;
  cout<<"Test: 7.000^2="<<power<2>(7.000);
  Dipole_Handler::ShowCalcBox();
  cout<<"Number of Dipole_Handler's = "<<Dipole_Handler::InStore<<endl;
  cout<<"=============================================================="<<endl;

  {
    unsigned total=0; //total=0;
    unsigned count=0;

    for(unsigned i=1; i<=total; ++i) {

      //cout<<endl; cin>>enter; cout<<endl;
      //cout<<"=>"<<i<<endl;
      bool control=!(i%50000);

      Vec4D pl(45.0, 20.0,-5.0, 40.0);
      Vec4D pr(45.0,-20.0, 5.0,-40.0);
      //Dipole::Branch     b1(info.quark.u,pl);
      Dipole::Glubranch  b1(pl);
      //Dipole::Antibranch a1(info.antiquark.u,pr);
      Dipole::Glubranch  a1(pr);

      Dipole* pDin=NULL;
      Dipole::Glubranch* pGlu=NULL;

      {

	bool dummy;
	Dipole Dip(b1,a1);
	Dipole_Handler H(Dip);

#ifdef EMISSIONTEST_OUTPUT
	cout<<Dip<<endl;
#endif

	if(!H.InduceGluonEmission()) {
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
