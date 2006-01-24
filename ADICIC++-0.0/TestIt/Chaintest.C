//bof
//Version: 3 ADICIC++-0.0/2005/09/21

//Chaintest.C - testing the first chaining.



#include <typeinfo>
#include <cassert>
#include <iostream>
#include <iomanip>
#include <ioextra>
#include <enumextra>
#include <mathextra>
#include "Message.H"
#include "Dipole_Parameter.H"
#include "Paraminit.H"
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

  msg.SetModifiable(true);
  Dipole_Flavour_Init dfi(true);

  cout<<endl;
  cout<<cout.precision(6)<<endl;
  cout<<endl;

  Lev1 kkk;
  kkk.Print();
  kkk.Vary();

  Lev2 uuu;
  //uuu.lev3.t=11;    //does not work - as wished
  //uuu.lev3.T()=16;    //does not work - as wished

  cout<<endl;

  Dipole_Handler::ListCalcBox();

  Dipole_Parameter::Show();
  Sudakov_Calculator::ShowEnvironment();
  //Sudakov_Calculator::AdjustEnvironment();
  //Sudakov_Calculator::ShowEnvironment();
  Dipole_Handler::ShowCalcBox();
  //Dipole_Handler::AdjustCalcBox();
  //Dipole_Handler::AdjustCalcBox();
  //Dipole_Handler::ShowCalcBox();

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
  cout<<endl; cin>>enter; cout<<endl;
  cout<<"=============================================================="<<endl;
  cout<<endl;

  cout<<"=================================="<<endl;
  cout<<" Testing the chaining stepwisely. "<<endl;
  cout<<"=================================="<<endl;

  {

    cout<<"Running?="<<Sudakov_Calculator::IsAlphaSRunning()<<endl;
    cout<<"ASApp="<<Sudakov_Calculator::AlphaSApprox()<<endl;
    cout<<"ASCor(700)="<<Sudakov_Calculator::AlphaSCorr(700.0)<<endl;
    cout<<"ASCor(8100)="<<Sudakov_Calculator::AlphaSCorr(8100.0)<<endl;
    cout<<"SudakovInit="<<Sudakov_Calculator::AdjustEnvironment()<<endl;
    cout<<"SudakovInit="<<Sudakov_Calculator::AdjustEnvironment()<<endl;
    cout<<"Running?="<<Sudakov_Calculator::IsAlphaSRunning()<<endl;
    cout<<"ASApp="<<Sudakov_Calculator::AlphaSApprox()<<endl;
    cout<<"ASCor(700)="<<Sudakov_Calculator::AlphaSCorr(700.0)<<endl;
    cout<<"ASCor(8100)="<<Sudakov_Calculator::AlphaSCorr(8100.0)<<endl;

    Sudakov_Calculator::ShowEnvironment();
    Dipole_Handler::ShowCalcBox();

    //Dipole_Parameter_Init dpi;

    cout<<"\nRunning?="<<Sudakov_Calculator::IsAlphaSRunning()<<endl;
    cout<<"ASApp="<<Sudakov_Calculator::AlphaSApprox()<<endl;
    cout<<"ASCor(700)="<<Sudakov_Calculator::AlphaSCorr(700.0)<<endl;
    cout<<"ASCor(8100)="<<Sudakov_Calculator::AlphaSCorr(8100.0)<<endl;
    cout<<"SudakovInit="<<Sudakov_Calculator::AdjustEnvironment()<<endl;
    cout<<"SudakovInit="<<Sudakov_Calculator::AdjustEnvironment()<<endl;
    cout<<"Running?="<<Sudakov_Calculator::IsAlphaSRunning()<<endl;
    cout<<"ASApp="<<Sudakov_Calculator::AlphaSApprox()<<endl;
    cout<<"ASCor(700)="<<Sudakov_Calculator::AlphaSCorr(700.0)<<endl;
    cout<<"ASCor(8100)="<<Sudakov_Calculator::AlphaSCorr(8100.0)<<endl;

    cout<<endl; cin>>enter; cout<<endl;
    //abort();

    {
      unsigned total=13;//20000;
      unsigned qcount=0, gcount=0, gluons=0, cd=0, cu=0, cs=0, cc=0, cb=0;

      for(unsigned i=1; i<=total; ++i) {

	Vec4D pl(45.0, 20.0,-5.0, 40.0);
	Vec4D pr(45.0,-20.0, 5.0,-40.0);
	Dipole::Branch     b1(info.quark.u,pl);
	Dipole::Glubranch  g1(pl);
	Dipole::Antibranch a1(info.antiq.u,pr);
	Dipole::Glubranch  g2(pr);

	Chain cha(b1,a1,Chain::Initiator::simple_epem);
	//Chain cha(g1,g2,Chain::Initiator::simple_epem);
	Chain chtest(cha);
	Chain_Handler H(cha);

	list<Chain*> chainlist;
	chainlist.push_front(&cha);
	list<Chain*>::iterator itcha=chainlist.begin();

	cha.Print();
	cout<<cha.MaxParticleNumber()<<endl;
	cout<<cha.MaxDipoleNumber()<<endl;
	if(i==total) chtest.Print();
	cout<<endl;
	cout<<"\e[7m\e[31m                    \e[0m";
	cout<<"\e[1m\e[40m\e[33mSTART: "<<i<<"\e[0m";
	cout<<"\e[7m\e[31m                    \e[0m"<<endl;

	while(H.EvolveChainByOneStep()) {
	  bool below;
	  Flavour fla;
	  Chain* pCha=NULL;
	  Recoil_Tool* pRct=NULL;
	  H.DecoupleNew(pRct,pCha,fla,below); assert(pRct==NULL);
	  if(pCha) {
	    ++qcount;
	    if(below) { ++itcha; itcha=chainlist.insert(itcha,pCha); --itcha;}
	    else      { itcha=chainlist.insert(itcha,pCha); ++itcha;}
	  } else {
	    ++gcount;
	  }
	  switch(fla.Kfcode()) {
	  case kf::gluon: ++gluons; break;
	  case kf::d    : ++cd; break;
	  case kf::u    : ++cu; break;
	  case kf::s    : ++cs; break;
	  case kf::c    : ++cc; break;
	  case kf::b    : ++cb; break;
	  default       : cout<<fla<<endl; assert(0);
	  }
	  cout<<cha<<endl;
	  Vec4D test;
	  cout<<cha.CheckMomentumConservation(test)<<"\t"<<test<<"\n\n";
	}

	cout<<"\e[7m\e[31m                    \e[0m";
	cout<<"\e[1m\e[40m\e[33mSTOP: "<<i<<"\e[0m";
	cout<<"\e[7m\e[31m                    \e[0m"<<endl;

	if(i==total) {
	  cha.Print();
	  //chtest.Print();
	  cout<<endl;
	  cout<<om::greenbg<<"Operator= test:"<<om::reset<<endl;
	  chtest=cha;
	  chtest.Print();
	  cout<<endl<<om::greenbg<<"Copy constructor test:"<<om::reset<<endl;
	  Chain chacopy(cha);
	  chacopy.Print();
	  cout<<endl<<om::greenbg<<"Extracting test:"<<om::reset<<endl;
	  Particle_List plist;
	  cha.ExtractPartons(plist);
	  cout<<plist<<endl;
	  for(Particle_List::iterator pit=plist.begin(); pit!=plist.end(); ++pit)
	    if(*pit) delete (*pit);
	}

	cout<<om::brownbg<<"THIS IS NOW THE RESULT......."<<om::reset;
	for(list<Chain*>::const_iterator cit=chainlist.begin();
	    cit!=chainlist.end(); ++cit) {
	  const Chain& chi=**cit;
	  if((*cit)==(&cha)) cout<<"This is the root:";
	  chi.Print();
	  cout<<chi.MaxParticleNumber()<<endl;
	  cout<<chi.MaxDipoleNumber()<<endl;
	}
	for(list<Chain*>::const_iterator cit=chainlist.begin();
	    cit!=chainlist.end(); ++cit) {
	  if((*cit)==(&cha)); else delete (*cit);
	}
	cout<<endl;

      }

      cout<<endl;
      cout<<"Total number of gluon  emissions (root chain)   = "<<gcount;
      cout<<"   ("<<1.0*gcount/total<<")   ["<<gluons<<"]."<<endl;
      cout<<"Total number of gluon splittings (extra chains) = "<<qcount;
      cout<<"   ("<<1.0*qcount/total<<")   ["<<cd+cu+cs+cc+cb<<"]."<<endl;
      cout<<"Total number of gluon splittings -> ddbar pairs = "<<cd;
      cout<<"   ("<<1.0*cd/qcount<<")."<<endl;
      cout<<"Total number of gluon splittings -> uubar pairs = "<<cu;
      cout<<"   ("<<1.0*cu/qcount<<")."<<endl;
      cout<<"Total number of gluon splittings -> ssbar pairs = "<<cs;
      cout<<"   ("<<1.0*cs/qcount<<")."<<endl;
      cout<<"Total number of gluon splittings -> ccbar pairs = "<<cc;
      cout<<"   ("<<1.0*cc/qcount<<")."<<endl;
      cout<<"Total number of gluon splittings -> bbbar pairs = "<<cb;
      cout<<"   ("<<1.0*cb/qcount<<")."<<endl;
      cout<<"Ratio of the # of splittings to the # of emissions = "
	  <<1.0*qcount/gcount<<" ."<<endl;
      //cout<<"Ratio of the # of produced q-qbar pairs to the # of gluons = "
      //  <<1.0*qcount/gcount<<" ."<<endl;
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
  Sudakov_Calculator::ShowEnvironment();
  Dipole_Handler::ShowCalcBox();
  cout<<"Number of Dipole_Handler's = "<<Dipole_Handler::InStore<<endl;
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

	bool  bel;
	xbool rec;
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
