//bof
//Version: 4 ADICIC++-0.0/2006/05/23

//ICascadetest.C - testing the first cascading for the IS case.



#include <typeinfo>
#include <cassert>
#include <iostream>
#include <iomanip>
#include <ioextra>
#include <enumextra>
#include <mathextra>
#include "Message.H"
#include "Run_Parameter.H"
#include "Dipole_Parameter.H"
#include "Paraminit.H"
#include "Adicic.H"
#include "Dipole_Handler.H"
#include "Chain_Handler.H"
#include "Cascade_Handler.H"


#define ICASCADETEST_OUTPUT ICASCADETEST_OUTPUT
#undef  ICASCADETEST_OUTPUT





using namespace std;
using namespace ATOOLS;
using namespace ADICIC;





int main() {

  //cout<<endl;
  //cout<<cout.precision(6)<<endl;
  //cout<<endl;

  msg.SetModifiable(true);
  msg.SetLevel(15);
  Dipole_Flavour_Init dfi(true);

  dpa.Show();
  Sudakov_Calculator::ShowEnvironment();
  Dipole_Handler::ShowCalcBox();

  cout<<sqrt(-7.0)<<endl;

  //===========================================================================
  cout<<om::bluebg<<string(100,'+')<<om::reset<<endl;
  cout<<endl; cin>>enter; cout<<endl;
  cout<<om::bluebg<<string(100,'+')<<om::reset<<endl;
  //===========================================================================

  cout<<"====================================="<<endl;
  cout<<" Testing the I cascading stepwisely. "<<endl;
  cout<<"====================================="<<endl;

  extern Run_Parameter ATOOLS::rpa;
  rpa.gen.SetEcms(1960.0);
  Dipole_Parameter_Init dpi;

  cout<<"\nRunning?="<<Sudakov_Calculator::IsAlphaSRunning()<<endl;
  cout<<"NfFix="<<dpa.sud.NfFix()<<endl;
  cout<<"ASApp="<<Sudakov_Calculator::AlphaSApprox()<<endl;
  cout<<"ASCor(700)="<<Sudakov_Calculator::AlphaSCorr(700.0)<<endl;
  cout<<"ASCor(8100)="<<Sudakov_Calculator::AlphaSCorr(8100.0)<<endl;
  cout<<"Nf(1.20)="<<Sudakov_Calculator::Nf(1.2)<<endl;
  cout<<"Nf(8100)="<<Sudakov_Calculator::Nf(8100.0)<<endl;

  dpa.Show();

  Sudakov_Calculator::ShowEnvironment();
  Dipole_Handler::ShowCalcBox();

  //===========================================================================
  cout<<om::bluebg<<string(100,'+')<<om::reset<<endl;
  cout<<endl; cin>>enter; cout<<endl;
  cout<<om::bluebg<<string(100,'+')<<om::reset<<endl;
  //===========================================================================

  if(1) {

    unsigned total=33;//0000;//1000;//20000;
    unsigned noem=0;
    unsigned fail=0;

    Cascade_Handler H;

    for(unsigned i=1; i<=total; ++i) {

      //bool control=!(i%10000);

      //Vec4D pl(58.865, 0.0, 0.0, 58.865);    //
      //Vec4D pr(114.69, 0.0, 0.0,-114.69);
      //Vec4D pl(2.2285e+01, 0.0, 0.0,  2.2285e+01);    //
      //Vec4D pr(9.2635e+01, 0.0, 0.0, -9.2635e+01);
      Vec4D pl(4.5617e+01, 0.0, 0.0, -4.5617e+01);    //
      Vec4D pr(4.6770e+01, 0.0, 0.0, 4.6770e+01);    //pr[0]=pr[3]=78.91795;
      //Vec4D pl(7.1758e+00, 0.0, 0.0, 7.1758e+00);    //
      //Vec4D pr(3.4202e+02, 0.0, 0.0,-3.4202e+02);
      //Vec4D pl(1.3281e+01, 0.0, 0.0, 1.3281e+01);    //c
      //Vec4D pr(1.5894e+02, 0.0, 0.0,-1.5894e+02);    //anti-c

      Dipole::Antibranch a1(info.quark.d,pl);    //Incoming d quark!
      Dipole::Glubranch  g2(pr,true);
      Dipole::Branch     b1(info.antiq.d,pr);    //Incoming anti d quark!
      Dipole::Glubranch  g1(pl,true);

      Cascade cas;

      Dipole Dip(b1,a1);    //Convenient calculation.
      dpv.evo.SetChainParticleLimit(5);
      dpv.evo.SetChainCorrelationLimit(5);
      dpv.sud.SetMaxIIScale(Dip.InvMass());////////////////////////////////////
      //dpv.sud.SetMaxIIScale();///////////////////////////////////////////////
      double iscale=dpa.MaxIIInvScale(Dip.InvMass());

      cas.AddChain(b1,a1,NULL,iscale);
      //cas.AddChain(b1,a1);
      //cas.AddChain(g1,g2);
      //cas.AddChain(g1,g2);

      //Cascade cascop;    //Then the dipole and particle numbers are in order.
      Cascade cascop(cas);

      cas|H;
      //If you would like to shorten a bit, comment in:
      //if(H.EvolveCascade()); else ++fail; cas|0; continue;

      cout<<"HardScale="<<dpa.sud.MaxIIK2t()
	  <<"  "<<sqrt(dpa.sud.MaxIIK2t())<<endl;
      cout<<"-/SHatMax="<<dpa.MaxIISHat(Dip.InvMass())
	  <<"  "<<sqrt(dpa.MaxIISHat(Dip.InvMass()));
      cas.Print();
      cout<<cas.MaxChainNumber()<<endl;
      if(i==total) cascop.Print();
      cout<<endl;
      cout<<"\e[7m\e[31m                    \e[0m";
      cout<<"\e[1m\e[40m\e[33mSTART: "<<i<<"\e[0m";
      cout<<"\e[7m\e[31m                    \e[0m\n\n";

      if(H.EvolveCascade()); else ++fail;

      if(cas.DipoleNumber()==1) ++noem;
      Vec4D test;
      cout<<"Momentum conservation check = ";
      cout<<cas.CheckMomentumConservation(test)<<"\t"<<test<<"\n\n";

      cout<<"\e[7m\e[31m                    \e[0m";
      cout<<"\e[1m\e[40m\e[33mSTOP: "<<i<<"\e[0m";
      cout<<"\e[7m\e[31m                    \e[0m"<<endl;

      if(i==total) {
	cas.Print();
	cout<<"I Number = "<<cas.INumber()<<endl;
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
	assert(cas.ExtractPartons(lists).flag);
	for(list<Particle_List>::iterator loc=lists.begin();
	    loc!=lists.end(); ++loc) {
	  cout<<(*loc);
	  for(Particle_List::iterator pit=(*loc).begin();
	      pit!=(*loc).end(); ++pit)
	    if(*pit) delete (*pit);
	}
      }

      cout<<"\n\n"<<om::redbg<<"THIS IS NOW THE RESULT......."<<om::reset;
      cas.Print();
      cout<<cas.MaxChainNumber()<<"\n";
      cout<<"I Number = "<<cas.INumber()<<"\n";
      cout<<endl;

      cas|0;    //Otherwise one gets a warning.

    }

    H.PrintCounter();
    cout<<"Total number of no emissions = "<<noem;
    cout<<"   ("<<1.0*noem/total<<")."<<endl;
    cout<<"Total number of failures = "<<fail;
    cout<<"   ("<<1.0*fail/total<<")."<<endl;

  }

  //===========================================================================
  cout<<om::bluebg<<string(100,'+')<<om::reset<<endl;
  cout<<om::bluebg<<string(100,'+')<<om::reset<<endl<<endl;
  //===========================================================================

}





//eof
