//bof
//Version: 2 ADICIC++-0.0/2004/09/10

//Cascadetest.C - testing the first cascading.



#include <typeinfo>
#include <cassert>
#include <iostream>
#include <iomanip>
#include <ioextra>
#include <enumextra>
#include "Message.H"
#include "Dipole_Parameter.H"
#include "Paraminit.H"
#include "Dipole.H"
#include "Dipole_Handler.H"
#include "Sudakov_Calculator.H"
#include "Chain.H"
#include "Chain_Handler.H"
#include "Cascade.H"
#include "Cascade_Handler.H"


#define CASCADETEST_OUTPUT CASCADETEST_OUTPUT
#undef  CASCADETEST_OUTPUT





using namespace std;
using namespace ATOOLS;
using namespace ADICIC;





int main() {

  msg.SetModifiable(true);
  cout<<Dipole_Flavour_Init::Status()<<endl;
  cout<<Dipole_Flavour_Init::DoIt(true)<<endl;
  cout<<Dipole_Flavour_Init::DoIt()<<endl;
  cout<<Dipole_Flavour_Init::Status()<<endl;

  cout<<endl;
  cout<<cout.precision(6)<<endl;
  cout<<endl;

  Dipole_Parameter::Show();
  Sudakov_Calculator::ShowParameters();
  Sudakov_Calculator::AdjustParameters();
  Sudakov_Calculator::ShowParameters();
  Dipole_Handler::ShowCalcBox();
  Dipole_Handler::AdjustCalcBox();
  Dipole_Handler::AdjustCalcBox();
  Dipole_Handler::ShowCalcBox();
  Chain_Handler::ShowParameters();
  Chain_Handler::AdjustParameters();
  Chain_Handler::ShowParameters();

  cout<<endl;
  cout<<"=============================================================="<<endl;
  cout<<endl;

  cout<<"======================================="<<endl;
  cout<<" Testing the cascading - Preparations. "<<endl;
  cout<<"======================================="<<endl;

  cout<<endl; cin>>enter; cout<<endl;

  {//muss noch diese tests machen!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Chain cha;
    cha.Print();
    cout<<cha.MaxParticleNumber()<<endl;
    cout<<cha.MaxDipoleNumber()<<endl;

    Vec4D pl(45.0, 20.0,-5.0, 40.0);
    Vec4D pr(45.0,-20.0, 5.0,-40.0);
    Dipole::Branch     b1(info.quark.u,pl);
    Dipole::Glubranch  g1(pl);
    Dipole::Antibranch a1(info.antiq.u,pr);
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
    cout<<"cha handling(0)="<<cha.IsHandled()<<endl;
    cout<<"ch1 handling(0)="<<ch1.IsHandled()<<endl;
    cha|0; ch1|0;
    cout<<"cha handling(0)="<<cha.IsHandled()<<endl;
    cout<<"ch1 handling(0)="<<ch1.IsHandled()<<endl;

    {
      Chain_Handler H1;
      Chain_Handler H2(ch1);

      cout<<"H1 docking(0)="<<H1.IsDocked()<<endl;
      cout<<"H2 docking(1)="<<H2.IsDocked()<<endl;
      cout<<"H2 docking cha (0)? "<<H2.IsDockedAt(cha)<<endl;
      cout<<"H2 docking ch1 (1)? "<<H2.IsDockedAt(ch1)<<endl;
      cout<<"ch1 handled by H1 (0)? "<<ch1.IsHandledBy(H1)<<endl;
      cout<<"ch1 handled by H2 (1)? "<<ch1.IsHandledBy(H2)<<endl;
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
	cout<<"H3 docking(1)="<<H3.IsDocked()<<endl;
	cout<<E.IsRing()<<endl;
	cout<<E<<endl;
      }
      cout<<"H3 docking(0)="<<H3.IsDocked()<<endl;
      cout<<endl;

      cout<<"(1)"<<(cha|H3)<<endl;
      cout<<"(0)"<<H3.DetachChain(&cha)<<endl;
      cout<<"(0)"<<H3.DetachChain(&ch1)<<endl;
      cout<<"(0)"<<H3.AttachChain(&cha)<<endl;
      cout<<"(0)"<<H3.AttachChain(&ch1)<<endl;
      cout<<"(0)"<<(cha|H1)<<endl;
      cout<<"(0)"<<(cha|H2)<<endl;
      cout<<"(0)"<<(cha|H3)<<endl;
      cout<<"(0)"<<(ch1|H1)<<endl;
      cout<<"(0)"<<(ch1|H2)<<endl;
      cout<<"(0)"<<(ch1|H3)<<endl;
      cout<<endl;
      ch1|0;
      cout<<"(0)"<<(ch1|H3)<<endl;
      cout<<"(1)"<<(ch1|H1)<<endl;
      cout<<"(0)"<<(ch1|H2)<<endl;
      cout<<cha<<endl<<ch1<<endl;
      cout<<"cha handled by H1(0)? "<<cha.IsHandledBy(H1)<<endl;
      cout<<"cha handled by H2(0)? "<<cha.IsHandledBy(H2)<<endl;
      cout<<"cha handled by H3(1)? "<<cha.IsHandledBy(H3)<<endl;
      cout<<"ch1 handled by H1(1)? "<<ch1.IsHandledBy(H1)<<endl;
      cout<<"ch1 handled by H2(0)? "<<ch1.IsHandledBy(H2)<<endl;
      cout<<"ch1 handled by H3(0)? "<<ch1.IsHandledBy(H3)<<endl;

      cout<<"(1)"<<cha.IsEmpty()<<endl;
      cout<<"(0)"<<ch1.IsEmpty()<<endl;
      cout<<"(0)"<<cha.IsRing()<<endl;
      cout<<"(0)"<<ch1.IsRing()<<endl;
      cout<<"(inc.)"<<cha.ChainType()<<endl;
      cout<<"(line)"<<ch1.ChainType()<<endl;
      cout<<ch1.ChainRoot()<<endl;

      cout<<H1.CompScale()<<endl;
      cout<<H2.CompScale()<<endl;
      cout<<H3.CompScale()<<endl;
      H1.RemoveNewProducts();
      H2.RemoveNewProducts();
      H3.RemoveNewProducts();

      cout<<"Number of Chain_Handler's = "<<Chain_Handler::InStore<<endl;
    }

    cout<<"cha handling="<<cha.IsHandled()<<endl;
    cout<<"ch1 handling="<<ch1.IsHandled()<<endl;
    cout<<"Number of Chain_Handler's = "<<Chain_Handler::InStore<<endl;

  }

  cout<<"=============================================================="<<endl;
  cout<<endl;
  cout<<"=============================================================="<<endl;
  cout<<endl;

  cout<<"==================================="<<endl;
  cout<<" Testing the cascading stepwisely. "<<endl;
  cout<<"==================================="<<endl;

  cout<<endl; cin>>enter; cout<<endl;

  {

    cout<<"Running?="<<Sudakov_Calculator::IsAlphaSRunning()<<endl;
    cout<<"MinScale="<<Sudakov_Calculator::MinOfK2t()<<endl;
    cout<<"MaxScale="<<Sudakov_Calculator::MaxOfK2t()<<endl;
    cout<<"ASFix="<<Sudakov_Calculator::AlphaSFix()<<endl;
    cout<<"ASApp="<<Sudakov_Calculator::AlphaSApprox()<<endl;
    cout<<"ASCor(700)="<<Sudakov_Calculator::AlphaSCorr(700.0)<<endl;
    cout<<"ASCor(8100)="<<Sudakov_Calculator::AlphaSCorr(8100.0)<<endl;
    cout<<"SudakovInit="<<Sudakov_Calculator::Init(NULL)<<endl;
    cout<<"SudakovInit="<<Sudakov_Calculator::Init(NULL)<<endl;
    cout<<"Running?="<<Sudakov_Calculator::IsAlphaSRunning()<<endl;
    cout<<"ASFix="<<Sudakov_Calculator::AlphaSFix()<<endl;
    cout<<"ASApp="<<Sudakov_Calculator::AlphaSApprox()<<endl;
    cout<<"ASCor(700)="<<Sudakov_Calculator::AlphaSCorr(700.0)<<endl;
    cout<<"ASCor(8100)="<<Sudakov_Calculator::AlphaSCorr(8100.0)<<endl;

    Sudakov_Calculator::ShowParameters();
    Dipole_Handler::ShowCalcBox();
    Chain_Handler::ShowParameters();
    cout<<"\nDip Param Init?: "<<Dipole_Parameter_Init::Status()<<endl;
    //cout<<"\nDo Init: "<<Dipole_Parameter_Init::DoIt()<<endl;////////////////
    //cout<<"Do Init: "<<Dipole_Parameter_Init::DoIt()<<endl;//////////////////
    cout<<"Dip Param Init?: "<<Dipole_Parameter_Init::Status()<<endl;

    cout<<"\nRunning?="<<Sudakov_Calculator::IsAlphaSRunning()<<endl;
    cout<<"MinScale="<<Sudakov_Calculator::MinOfK2t()<<endl;
    cout<<"MaxScale="<<Sudakov_Calculator::MaxOfK2t()<<endl;
    cout<<"ASFix="<<Sudakov_Calculator::AlphaSFix()<<endl;
    cout<<"ASApp="<<Sudakov_Calculator::AlphaSApprox()<<endl;
    cout<<"ASCor(700)="<<Sudakov_Calculator::AlphaSCorr(700.0)<<endl;
    cout<<"ASCor(8100)="<<Sudakov_Calculator::AlphaSCorr(8100.0)<<endl;
    cout<<"SudakovInit="<<Sudakov_Calculator::Init(NULL)<<endl;
    cout<<"SudakovInit="<<Sudakov_Calculator::Init(NULL)<<endl;
    cout<<"Running?="<<Sudakov_Calculator::IsAlphaSRunning()<<endl;
    cout<<"ASFix="<<Sudakov_Calculator::AlphaSFix()<<endl;
    cout<<"ASApp="<<Sudakov_Calculator::AlphaSApprox()<<endl;
    cout<<"ASCor(700)="<<Sudakov_Calculator::AlphaSCorr(700.0)<<endl;
    cout<<"ASCor(8100)="<<Sudakov_Calculator::AlphaSCorr(8100.0)<<endl;

    cout<<"\n"<<Dipole_Parameter::ForceFirstInit()<<endl;
    //Dipole_Parameter::Reset();
    //Sudakov_Calculator::AdjustParameters();
    //Dipole_Handler::AdjustCalcBox();
    //Chain_Handler::AdjustParameters();

    cout<<endl; cin>>enter; cout<<endl;
    //abort();

    {
      unsigned total=3;//20000;
      unsigned fail=0;

      Cascade_Handler H;

      for(unsigned i=1; i<=total; ++i) {

	Vec4D pl(45.0, 20.0,-5.0, 40.0);
	Vec4D pr(45.0,-20.0, 5.0,-40.0);
	Dipole::Branch     b1(info.quark.u,pl);
	Dipole::Glubranch  g1(pl);
	Dipole::Antibranch a1(info.antiq.u,pr);
	Dipole::Glubranch  g2(pr);

	Cascade cas;

	cas.AddChain(b1,a1);
	cas.AddChain(b1,a1);
	//cas.AddChain(g1,g2);
	//cas.AddChain(g1,g2);

	Cascade cascop(cas);

	cas|H;

	cas.Print();
	cout<<cas.MaxChainNumber()<<endl;
	if(i==total) cascop.Print();
	cout<<endl;
	cout<<"\e[7m\e[31m                    \e[0m";
	cout<<"\e[1m\e[40m\e[33mSTART: "<<i<<"\e[0m";
	cout<<"\e[7m\e[31m                    \e[0m\n\n";

	if(H.EvolveCascade());
	else ++fail;

	Vec4D test;
	cout<<"Momentum conservation check = ";
	cout<<cas.CheckMomentumConservation(test)<<"\t"<<test<<"\n\n";

	cout<<"\e[7m\e[31m                    \e[0m";
	cout<<"\e[1m\e[40m\e[33mSTOP: "<<i<<"\e[0m";
	cout<<"\e[7m\e[31m                    \e[0m"<<endl;

	if(i==total) {
	  cas.Print();
	  //cascop.Print();
	  cout<<endl;
	  cout<<om::greenbg<<"Operator= test:"<<om::reset<<endl;
	  cascop=cas;
	  cascop.Print();
	  cout<<endl<<om::greenbg<<"Copy constructor test:"<<om::reset<<endl;
	  Cascade cascopy(cas);
	  cascopy.Print();
	  cout<<endl<<om::greenbg<<"Extracting test:"<<om::reset<<endl;
	  list<Particle_List> lists;
	  assert(cas.ExtractPartons(lists));
	  for(list<Particle_List>::iterator loc=lists.begin();
	      loc!=lists.end(); ++loc) {
	    cout<<(*loc);
	    for(Particle_List::iterator pit=(*loc).begin(); pit!=(*loc).end(); ++pit)
	      if(*pit) delete (*pit);
	  }
	}

	cout<<"\n\n"<<om::redbg<<"THIS IS NOW THE RESULT......."<<om::reset;
	cas.Print();
	cout<<cas.MaxChainNumber()<<"\n";
	cout<<endl;

	cas|0;    //Otherwise one gets a warning.

      }

      H.PrintCounter();
      cout<<"Total number of failures = "<<fail;
      cout<<"   ("<<1.0*fail/total<<")."<<endl;
    }

  }

  cout<<endl<<endl;
  cout<<cout.precision(6)<<endl<<endl<<om::greenbg;
  cout<<"====================================================================";
  cout<<om::reset<<endl;
  cout<<"Test: 8.123^7="<<power<7>(8.123)<<endl;
  cout<<"Test: 4.573^3="<<power<3>(4.573)<<endl;
  cout<<"Test: 7.000^2="<<power<2>(7.000);
  Dipole_Parameter::Show();
  Sudakov_Calculator::ShowParameters();
  Dipole_Handler::ShowCalcBox();
  cout<<"Number of Dipole_Handler's = "<<Dipole_Handler::InStore<<endl;
  Chain_Handler::ShowParameters();
  cout<<om::greenbg;
  cout<<"====================================================================";
  cout<<om::reset<<endl;

  {
    unsigned total=0; //total=0;
    unsigned count=0;

    for(unsigned i=1; i<=total; ++i) {

      bool control=!(i%50000);

      Vec4D pl(45.0, 20.0,-5.0, 40.0);
      Vec4D pr(45.0,-20.0, 5.0,-40.0);
      //Dipole::Branch     b1(info.quark.u,pl);
      Dipole::Glubranch  b1(pl);
      //Dipole::Antibranch a1(info.antiquark.u,pr);
      Dipole::Glubranch  a1(pr);

      Dipole* pDin=NULL;
      Dipole::Glubranch*  pGlu=NULL;
      Dipole::Branch*     pBan=NULL;
      Dipole::Antibranch* pAti=NULL;

      {

	bool bel;
	Trio rec;
	Dipole Dip(b1,a1);
	Dipole_Handler H(Dip);

	if(!H.InduceDipoleRadiation()) {
	  ++count;
	  if(control) cout<<i<<"\t\t"<<count<<"\t\t"<<Dip.ProdScale()<<endl;
	  continue;
	}

	if(control) cout<<i<<"\t\t"<<count<<"\t\t"<<Dip.ProdScale()<<endl;

	H.DecoupleNew(pDin,pGlu,pAti,pBan,bel,rec);

      }

      if(pDin) { delete pDin; if(pGlu) delete pGlu;}
      else { if(pBan) delete pBan; if(pAti) delete pAti;}

    }

    cout<<endl<<"Total number of failures="<<count;
    cout<<"   ("<<1.0*count/total<<")."<<endl;
  }

  cout<<"=============================================================="<<endl;

  cout<<endl;

}





//eof
