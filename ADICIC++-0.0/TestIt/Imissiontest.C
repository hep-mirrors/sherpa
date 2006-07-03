//bof
//Version: 4 ADICIC++-0.0/2006/06/29

//Imissiontest.C - testing the first emission for the IS case.



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
#include "Adicic.H"


#define IMISSIONTEST_OUTPUT IMISSIONTEST_OUTPUT
#undef  IMISSIONTEST_OUTPUT





using namespace std;
using namespace ATOOLS;
using namespace ADICIC;





int main() {

  msg.SetModifiable(true);
  msg.SetLevel(15);
  Dipole_Flavour_Init dfi(true);

  dpa.Show(); Dipole_Handler::ShowCalcBox();

  cout<<endl; cin>>enter; cout<<endl;
  cout<<"=============================================================="<<endl;

  {

    Vec4D pl(45.0, 20.0,-5.0, 40.0);
    Vec4D pr(45.0,-20.0, 5.0,-40.0);
    Dipole D1;
    Dipole::Antibranch a0(info.quark.s,pl);
    Dipole::Branch     b0(info.antiq.t,pl);
    Dipole D2(b0,a0,33);    //II dipole!
    cout<<D1<<endl<<D2<<endl;
    D1.PrintTowers(); D2.PrintTowers();

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

  extern Run_Parameter ATOOLS::rpa;
  rpa.gen.SetEcms(1960.0);//90.0);
  Dipole_Parameter_Init dpi;

  cout<<"\nRunning?="<<Sudakov_Calculator::IsAlphaSRunning()<<endl;
  cout<<"NfFix="<<dpa.sud.NfFix()<<endl;
  cout<<"ASApp="<<Sudakov_Calculator::AlphaSApprox()<<endl;
  cout<<"ASCor="<<Sudakov_Calculator::AlphaSCorr(700.0)<<endl;
  cout<<"ASCor="<<Sudakov_Calculator::AlphaSCorr(8100.0)<<endl;
  cout<<"Nf(1.20)="<<Sudakov_Calculator::Nf(1.2)<<endl;
  cout<<"Nf(8100)="<<Sudakov_Calculator::Nf(8100.0)<<endl;

  dpa.Show();

  Sudakov_Calculator::ShowEnvironment();
  Dipole_Handler::ShowCalcBox();

  for(int i=0; i<10; ++i) {////////////////////////////////////////////////////

    TEMP::CPTEST=true;

    Vec4D pl(2.2285e+01, 0.0, 0.0,  2.2285e+01);
    Vec4D pr(9.2635e+01, 0.0, 0.0, -9.2635e+01);
    //Vec4D pl(4.5617e+01, 0.0, 0.0, -4.5617e+01);    //
    //Vec4D pr(4.6770e+01, 0.0, 0.0, 4.6770e+01);

    Dipole::Antibranch a1(info.quark.d,pl);
    //Dipole::Glubranch  a1(pl,true);
    //Dipole::Branch     b1(info.antiq.d,pr);
    Dipole::Glubranch  b1(pr,true);

    Dipole_Handler::Carrier box;

    {

      unsigned trials=1;
      Dipole Dip(b1,a1);
      dpv.sud.SetMaxIIScale(Dip.InvMass());////////////////////////////////////
      //Dip.SetEvolScales(dpa.MaxIIInvScale(Dip.InvMass()));
      Dip.SetEvolScales(dpa.HighestIIInvScale(Dip.InvMass()));
      Dipole Cop(Dip,true);
      Dipole_Handler H(Dip);
      H.ShowSudakov(); cout<<Dip<<endl; Dip.PrintTowers();

      while(H.InduceDipoleRadiation()==false) ++trials;
      cout<<"Trials for a successful radiation="<<trials<<endl;

      H.ShowRecoil();

      assert(H.FinishDipoleRadiation());
      //assert(H.FinishDipoleRadiation());    //as wanted, leads to abort.

      char stat=H.Status();
      H.DecoupleNew(box);

      if(stat=='g') {
	assert(!box.pAqu && !box.pQua && box.pGlu);
	Dip.SetFactScale()=box.pDip->SetFactScale()=Dip.EmitScale();
	cout<<  "\t        p1="
	    <<b1.Momentum()<<" \t "<<b1.Momentum().Abs2();
	cout<<"\n\t        p2="
	    <<box.pGlu->Momentum()<<" \t "<<box.pGlu->Momentum().Abs2();
	cout<<"\n\t        p3="
	    <<a1.Momentum()<<" \t "<<a1.Momentum().Abs2();
	Vec4D sum=b1.Momentum()-box.pGlu->Momentum()+a1.Momentum();
	cout<<"\n\tZ=p1-p2+p3="<<sum<<" \t "<<sum.Abs2();
	//cout<<"\n\t       ";
	//for(char i=0; i<4; ++i) cout<<b1.Momentum()[i]/pl[i]<<" \t ";
	//cout<<endl<<"\t       ";
	//for(char i=0; i<4; ++i) cout<<box.pGlu->Momentum()[i]/pl[i]<<" \t ";
	//cout<<endl<<"\t       ";
	//for(char i=0; i<4; ++i) cout<<a1.Momentum()[i]/pl[i]<<" \t ";
	cout<<endl;
	const Dipole& Din=*box.pDip;
	cout<<"Recoil is "<<box.RecoComp<<"; and, "
	    <<"is order given as olddip newdip? "<<box.DipOrder<<endl;
	cout<<Cop<<endl;
	if(box.DipOrder) {
	  cout<<Dip<<endl<<Din<<endl;
	  Dip.PrintTowers(); Din.PrintTowers();
	} else {
	  cout<<Din<<endl<<Dip<<endl;
	  Din.PrintTowers(); Dip.PrintTowers();
	}
      } else if(stat=='q') {
	assert(box.pDip && box.pAqu && box.pQua && box.pGlu);
	Dip.SetFactScale()=box.pDip->SetFactScale()=Dip.EmitScale();
	cout<<"Is q(0) or qbar(1) emission? "<<box.DipOrder<<endl;
	if(box.DipOrder==false) {
	  cout<<  "\t        p2="
	      <<box.pQua->Momentum()<<" \t "<<box.pQua->Momentum().Abs2();
	  cout<<"\n\t        p1="
	      <<box.pGlu->Momentum()<<" \t "<<box.pGlu->Momentum().Abs2();
	  cout<<"\n\t        p3="
	      <<a1.Momentum()<<" \t "<<a1.Momentum().Abs2();
	  Vec4D sum=box.pGlu->Momentum()-box.pQua->Momentum()+a1.Momentum();
	  cout<<"\n\tZ=p1-p2+p3="<<sum<<" \t "<<sum.Abs2();
	} else {
	  cout<<  "\t        p1="
	      <<b1.Momentum()<<" \t "<<b1.Momentum().Abs2();
	  cout<<"\n\t        p3="
	      <<box.pGlu->Momentum()<<" \t "<<box.pGlu->Momentum().Abs2();
	  cout<<"\n\t        p2="
	      <<box.pAqu->Momentum()<<" \t "<<box.pAqu->Momentum().Abs2();
	  Vec4D sum=box.pGlu->Momentum()-box.pAqu->Momentum()+b1.Momentum();
	  cout<<"\n\tZ=p1-p2+p3="<<sum<<" \t "<<sum.Abs2();
	}
	cout<<endl;
	const Dipole& Din=*box.pDip;
	cout<<"Recoil is "<<box.RecoComp<<endl;
	cout<<Cop<<endl;
	if(box.DipOrder) {
	  cout<<Dip<<endl<<Din<<endl;
	  Dip.PrintTowers(); Din.PrintTowers();
	} else {
	  cout<<Din<<endl<<Dip<<endl;
	  Din.PrintTowers(); Dip.PrintTowers();
	}
      } else abort();

      H.ShowSudakov(); H.ShowRecoil();

    }

    if(box.pDip) delete box.pDip;
    if(box.pGlu) delete box.pGlu;
    if(box.pQua && !box.DipOrder) delete box.pQua;
    if(box.pAqu &&  box.DipOrder) delete box.pAqu;

    TEMP::CPTEST=false;

  }

  cout<<om::greenbg;
  cout<<"====================================================================";
  cout<<om::reset<<endl;
  cout<<endl; cin>>enter; cout<<endl;
  cout<<om::greenbg;
  cout<<"====================================================================";
  cout<<om::reset<<endl;

  {
    unsigned total=2000;//000;
    unsigned count=0;
    unsigned gluons=0;
    unsigned quarks=0, cd=0, cu=0, cs=0, cc=0, cb=0;
    unsigned anti=0;

    Multiflavour   mufl(sf::stop,Flavour());
    Sudakov_Result sure; sure.Isr.resize(sr::stop,0.0);
    mufl[sf::plusini]=mufl[sf::plusfin]=Flavour(kf::u);
    mufl[sf::miusini]=mufl[sf::miusfin]=Flavour(kf::u,1);
    sure.Isr[sr::mdip]=90.0;
    sure.Isr[sr::fasc]=sqr(90.0);
    sure.Isr[sr::xpini]=0.04;
    sure.Isr[sr::xmini]=0.052712411;
    sure.Isr[sr::shat]=sqr(370.992);
    sure.Isr[sr::xpfin]=0.09;
    sure.Isr[sr::xmfin]=0.398083738;
    sure.Isr[sr::kt]=17.23;
    double pwgt=Sudakov_Calculator::PlusPDFCorr(mufl,sure);
    double mwgt=Sudakov_Calculator::MinusPDFCorr(mufl,sure);
    cout<<pwgt<<" : "<<mwgt<<endl;

    //abort();

    for(unsigned i=1; i<=total; ++i) {

      bool control=!(i%1000);//!(i%50000);

      //Vec4D pl(45.0, 20.0,-5.0, 40.0);    //
      //Vec4D pr(45.0,-20.0, 5.0,-40.0);
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

      Dipole::Antibranch a1(info.quark.u,pl);    //Incoming quark.
      //Dipole::Glubranch  a1(pr,true);
      Dipole::Branch     b1(info.antiq.u,pr);    //To get an II dipole.
      //Dipole::Glubranch  b1(pl,true);

      Dipole_Handler::Carrier box;

      {

	Dipole Dip(b1,a1);
	//Dipole_Particle* top=Dip.GetTopBranchPointer().operator->();
	//Dip.RenewBranch(*static_cast<Dipole::Branch*>(top));
	Dipole_Handler H(Dip);
	//H.ShowSudakov(); H.ShowRecoil();
	dpv.sud.SetMaxIIScale(Dip.InvMass());//////////////////////////////////
	//dpv.sud.SetMaxIIScale();/////////////////////////////////////////////
	Dip.SetEvolScales(dpa.MaxIIInvScale(Dip.InvMass()));
	//Dip.SetEvolScales(dpa.HighestIIInvScale(Dip.InvMass()));

#ifdef IMISSIONTEST_OUTPUT
	cout<<Dip<<endl;
	cout<<Dip.IsFF()<<endl;
	cout<<Dip.IsFI()<<endl;
	cout<<Dip.IsIF()<<endl;
	cout<<Dip.IsII()<<endl;
	cout<<dpa.sud.MaxIIK2t()<<"  "
	    <<sqrt(dpa.sud.MaxIIK2t())<<endl;
	cout<<dpa.MaxIISHat(Dip.InvMass())<<"  "
	    <<sqrt(dpa.MaxIISHat(Dip.InvMass()))<<endl;
#endif

	if(control) cout<<" "<<i<<":  ";
	if( H.InduceDipoleRadiation(1,control,i==total) &&
	    H.FinishDipoleRadiation() );
	else {
	  ++count;
	  H.DecoupleNew(box);
	  assert(box.pQua==NULL);
	  assert(box.pAqu==NULL);
	  assert(box.pGlu==NULL);
	  if(control) {//||1
	    cout<<" "<<i<<":  "<<count<<"\t\t"<<Dip.ProdScale()<<"\t\t";
	    cout<<box.pGlu; cout<<",";
	    cout<<box.pQua; cout<<",";
	    cout<<box.pAqu; cout<<endl;
	  }
	  continue;
	}
	char stat=H.Status();
	H.DecoupleNew(box);
	kf::code kfc;
	if(stat=='g') {
	  assert(!box.pAqu && !box.pQua && box.pGlu && box.pDip);
	  Dip.SetFactScale()=box.pDip->SetFactScale()=Dip.EmitScale();
	  kfc=kf::gluon;
	} else if(stat=='q') {
	  assert(box.pDip && box.pAqu && box.pQua && box.pGlu);
	  Dip.SetFactScale()=box.pDip->SetFactScale()=Dip.EmitScale();
	  kfc=box.pQua->Flav().Kfcode();
	  assert(box.pAqu->Flav().Kfcode()==kfc);
	  if(box.DipOrder) ++anti;
	} else abort();
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
	  cout<<" "<<i<<":  "<<count<<"\t\t"<<Dip.ProdScale()<<"\t\t";
	  if(box.pGlu) cout<<box.pGlu->Flav();
	  else cout<<box.pGlu; cout<<",";
	  if(box.pQua) cout<<box.pQua->Flav();
	  else cout<<box.pQua; cout<<",";
	  if(box.pAqu) cout<<box.pAqu->Flav();
	  else cout<<box.pAqu; cout<<"\t\t";
	  cout<<Flavour(kfc,box.DipOrder)<<","
	      <<box.DipOrder<<","
	      <<box.RecoComp<<endl;
	}

#ifdef IMISSIONTEST_OUTPUT
	if(box.pDip) {
	  const Dipole& Din=*box.pDip;
	  cout<<Dip<<endl<<Din<<endl;
	}
#endif

      }

      if(box.pDip) delete box.pDip;
      if(box.pGlu) delete box.pGlu;
      if(box.pQua && !box.DipOrder) delete box.pQua;
      if(box.pAqu &&  box.DipOrder) delete box.pAqu;

    }

    quarks=cd+cu+cs+cc+cb;

    cout<<endl;
    cout<<"Total number of   gluons="<<gluons;
    cout<<"   ("<<1.0*gluons/total<<")."<<endl;
    cout<<"Total number of   quarks="<<quarks;
    cout<<"   ("<<1.0*quarks/total<<")   anti="<<anti
	<<"   ("<<1.0*anti/quarks<<")."<<endl;
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
