//bof
//Version: 1 ADICIC++-0.0/2004/03/12

//Emissiontest.C - testing the first emission.



#include <typeinfo>
#include <cassert>
#include <iostream>
#include <ioextra>
#include <enumextra>
#include "Dipole.H"
#include "Dipole_Handler.H"


#define EMISSIONTEST_OUTPUT EMISSIONTEST_OUTPUT
#undef  EMISSIONTEST_OUTPUT





using namespace std;
using namespace ATOOLS;
using namespace ADICIC;





int main() {

  cout<<endl;
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
      Dipole_Handler H2(D2);

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
      {
	Dipole E(b1,a1);
	E|H3;
	cout<<"H3 docking="<<H3.IsDocked()<<endl;
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
    }

    cout<<"D1 handling="<<D1.IsHandled()<<endl;
    cout<<"D2 handling="<<D2.IsHandled()<<endl;

  }

  cout<<"=============================================================="<<endl;
  cout<<endl;

  cout<<"============================="<<endl;
  cout<<" Testing the first emission. "<<endl;
  cout<<"============================="<<endl;

  {
    Vec4D pl(45.0, 20.0,-5.0, 40.0);
    Vec4D pr(45.0,-20.0, 5.0,-40.0);

    Dipole::Branch b1(info.quark.u,pl);
    Dipole::Antibranch a1(info.antiquark.u,pr);

    Dipole* pDin=NULL;
    Dipole::Glubranch* pGlu=NULL;

    {

      Dipole Dip(b1,a1);
      Dipole_Handler H(Dip);

      cout<<Dip<<endl; Dip.PrintTowers();

      assert(H.InduceGluonEmission());

      H.DecoupleNewDipole(pDin); assert(pDin);
      H.DecoupleGlubranch(pGlu); assert(pGlu);

      const Dipole& Din=*pDin;

      if(Dip.IsType()==Dipole::qg) {
	cout<<Dip<<endl<<Din<<endl;
	Dip.PrintTowers(); Din.PrintTowers();
      } else {
	cout<<Din<<endl<<Dip<<endl;
	Din.PrintTowers(); Dip.PrintTowers();
      }

    }

    delete pDin;
    delete pGlu;

  }

  cout<<"=============================================================="<<endl;

  {
    unsigned total=4000000;
    unsigned count=0;

    for(unsigned i=1; i<=total; ++i) {

      //cout<<endl; cin>>enter; cout<<endl;
      //cout<<"=>"<<i<<endl;
      if(!(i%50000)) cout<<i<<"\t\t"<<count<<endl;

      Vec4D pl(45.0, 20.0,-5.0, 40.0);
      Vec4D pr(45.0,-20.0, 5.0,-40.0);
      Dipole::Branch b1(info.quark.u,pl);
      Dipole::Antibranch a1(info.antiquark.u,pr);
      Dipole* pDin=NULL;
      Dipole::Glubranch* pGlu=NULL;

      {

	Dipole Dip(b1,a1);
	Dipole_Handler H(Dip);

#ifdef EMISSIONTEST_OUTPUT
	cout<<Dip<<endl;
#endif

	if(!H.InduceGluonEmission()) { ++count; continue;}

	H.DecoupleNewDipole(pDin); assert(pDin);
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
