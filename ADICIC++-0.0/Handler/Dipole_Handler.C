//bof
//Version: 4 ADICIC++-0.0/2006/06/14

//Implementation of Dipole_Handler.H.


#ifdef __GNUC__
#if __GNUC__ >2
#include <ios>
#endif
#endif
#include <iomanip>
#include "Random.H"
#include "Poincare.H"
#include "Dipole_Parameter.H"
#include "Sudakov_Group.H"
#include "IISudakov_Group.H"
#include "Recoil_Calculator.H"
#include "Dipole_Handler.H"
#include "Histogram.H"
#include <histoextra>





using namespace std;
using namespace ATOOLS;
using namespace ADICIC;





//#include "....exa.cc"





//=============================================================================



//ostream& ADICIC::operator<<(ostream& ost, const ADICIC::...&) {
//}



//=============================================================================



//So far nothing is adjusted. There is no static Dipole_Handler.
bool Dipole_Handler::sf_1stadjust=false;
int Dipole_Handler::s_count=0;
const int& Dipole_Handler::InStore=Dipole_Handler::s_count;

const rdt::code Dipole_Handler::s_mpcode=rdt::iirecoil;

Dipole_Handler::Sudakovbox
Dipole_Handler::s_sumap=Dipole_Handler::Sudakovbox();
Dipole_Handler::Recoilbox
Dipole_Handler::s_remap=Dipole_Handler::Recoilbox();



//=============================================================================



Dipole_Handler::Dipole_Handler()
  : m_key(Dipole::incorrect,Radiation::incorrect),
    //m_key(make_pair(Dipole::incorrect,Radiation::incorrect)),
    p_sudakov(NULL), p_recoil(NULL),
    p_dip(NULL),
    p_dix(NULL), p_ban(NULL), p_ati(NULL), p_glu(NULL),
    f_below(false), f_gate(0), m_sur(), m_rer(s_mpcode) {

  static int ini=sf_1stadjust || AdjustCalcBox();
  assert(ini==1);//////////////////////////////////////////////////////////////
  s_count+=ini;
}





Dipole_Handler::Dipole_Handler(Dipole& dip)
  : m_key(Dipole::incorrect,Radiation::incorrect),
    //m_key(make_pair(Dipole::incorrect,Radiation::incorrect)),
    p_sudakov(NULL), p_recoil(NULL),
    p_dip(NULL),
    p_dix(NULL), p_ban(NULL), p_ati(NULL), p_glu(NULL),
    f_below(false), f_gate(0), m_sur(), m_rer(s_mpcode) {

  //Since there is no static DH, the very first DH can only be obtained from
  //the standard constructor!

  ++s_count;

  if(dip|*this) {    //While processing, f_gate is set.
    if(dip.IsHandledBy(*this)); else {
      cerr<<"\nBug: Wrong Dipole-Dipole_Handler connection emerged!\n";
      assert(dip.IsHandledBy(*this));
    }
  } else {
    cerr<<"\nMethod: "<<__PRETTY_FUNCTION__<<": "
	<<"Warning: Attaching Dipole failed!\n"<<endl;
  }

}





Dipole_Handler::~Dipole_Handler() {

  --s_count;

  //In case decoupling has not yet carried out, see Note1!
  if(p_dix) {
    delete p_dix;
    if(p_glu) delete p_glu;
    if(p_ban) delete p_ban;
    if(p_ati) delete p_ati;
  } else {
    if(p_ban) delete p_ban;
    if(p_ati) delete p_ati;
  }

  //assert(p_sudakov && p_dip || !p_sudakov && !p_recoil && !p_dip);
  assert(p_dip || !p_sudakov && !p_recoil && !p_dip);

  if(!p_dip) return;
  if(p_dip->IsHandledBy(*this)==false) {
    cerr<<"\nBug: Wrong Dipole-Dipole_Handler connection emerged!\n";
    assert(p_dip->IsHandledBy(*this));
  }

  *p_dip|0;

}



//=============================================================================



void Dipole_Handler::ListCalcBox() {    //Static.
  cout<<endl;
  cout<<"==============================================="<<endl;
  cout<<"Sudakovs and Recoilers for the Dipole_Handler's"<<endl;
  cout<<"-----------------------------------------------"<<endl;
  for(Sudakovbox::const_iterator it=s_sumap.begin(); it!=s_sumap.end(); ++it) {
    cout<<" "<<setw(5)<<it->first<<"  :  ";
    if(it->second) it->second->Which();
    else cout<<"Not initialized."<<endl;
  }
  cout<<"-----------------------------------------------"<<endl;
  for(Recoilbox::const_iterator it=s_remap.begin(); it!=s_remap.end(); ++it) {
    cout<<" "<<setw(5)<<(it->first).first<<"  ,  "
	<<setw(5)<<(it->first).second<<"  :  ";
    if(it->second) {
      cout<<"(code="<<it->second->IsType()<<") ";
      it->second->Which();
    }
    else cout<<"Not initialized."<<endl;
  }
  cout<<"==============================================="<<endl;
}





void Dipole_Handler::ShowCalcBox() {    //Static.
  cout<<endl;
  cout<<"======================================="<<endl;
  cout<<"Calculator box for the Dipole_Handler's"<<endl;
  cout<<"---------------------------------------"<<endl;
  cout<<"Number of Sudakov_Calculators in store = "
      <<Sudakov_Calculator::InStore<<"."<<endl;
  cout<<"Mimic Ariadne? "<<Sudakov_Calculator::Ariadne<<"."<<endl;
  cout<<"Number of  Recoil_Calculators in store = "
      <<Recoil_Calculator::InStore<<"."<<endl;
  cout<<"---------------------------------------"<<endl;
  if(s_sumap[Dipole::qqbar]) {
    s_sumap[Dipole::qqbar]->Which();
    s_sumap[Dipole::qqbar]->ShowSpecification();}
  if(s_sumap[Dipole::qg]) {
    s_sumap[Dipole::qg]->Which();
    s_sumap[Dipole::qg]->ShowSpecification();}
  if(s_sumap[Dipole::gqbar]) {
    s_sumap[Dipole::gqbar]->Which();
    s_sumap[Dipole::gqbar]->ShowSpecification();}
  if(s_sumap[Dipole::gg]) {
    s_sumap[Dipole::gg]->Which();
    s_sumap[Dipole::gg]->ShowSpecification();}
  if(s_sumap[Dipole::iiqbarq]) {
    s_sumap[Dipole::iiqbarq]->Which();
    s_sumap[Dipole::iiqbarq]->ShowSpecification();}
  if(s_sumap[Dipole::iiqbarg]) {
    s_sumap[Dipole::iiqbarg]->Which();
    s_sumap[Dipole::iiqbarg]->ShowSpecification();}
  if(s_sumap[Dipole::iigq]) {
    s_sumap[Dipole::iigq]->Which();
    s_sumap[Dipole::iigq]->ShowSpecification();}
  if(s_sumap[Dipole::iigg]) {
    s_sumap[Dipole::iigg]->Which();
    s_sumap[Dipole::iigg]->ShowSpecification();}
  cout<<"---------------------------------------"<<endl;
  cout<<" -  -  FF g emit   -  -  -  -  -  -  - "<<endl;
  if(s_remap[Key(Dipole::qqbar,Radiation::gluon)])
    s_remap[Key(Dipole::qqbar,Radiation::gluon)]->Which();
  else cout<<"  --\n";
  if(s_remap[Key(Dipole::qg,Radiation::gluon)])
    s_remap[Key(Dipole::qg,Radiation::gluon)]->Which();
  else cout<<"  --\n";
  if(s_remap[Key(Dipole::gqbar,Radiation::gluon)])
    s_remap[Key(Dipole::gqbar,Radiation::gluon)]->Which();
  else cout<<"  --\n";
  if(s_remap[Key(Dipole::gg,Radiation::gluon)])
    s_remap[Key(Dipole::gg,Radiation::gluon)]->Which();
  else cout<<"  --\n";
  cout<<" -  -  II g emit   -  -  -  -  -  -  - "<<endl;
  if(s_remap[Key(Dipole::iiqbarq,Radiation::gluon)])
    s_remap[Key(Dipole::iiqbarq,Radiation::gluon)]->Which();
  else cout<<"  --\n";
  if(s_remap[Key(Dipole::iiqbarg,Radiation::gluon)])
    s_remap[Key(Dipole::iiqbarg,Radiation::gluon)]->Which();
  else cout<<"  --\n";
  if(s_remap[Key(Dipole::iigq,Radiation::gluon)])
    s_remap[Key(Dipole::iigq,Radiation::gluon)]->Which();
  else cout<<"  --\n";
  if(s_remap[Key(Dipole::iigg,Radiation::gluon)])
    s_remap[Key(Dipole::iigg,Radiation::gluon)]->Which();
  else cout<<"  --\n";
  cout<<" -  -  FF g split  -  -  -  -  -  -  - "<<endl;
  if(s_remap[Key(Dipole::qg,Radiation::qbarbot)])
    s_remap[Key(Dipole::qg,Radiation::qbarbot)]->Which();
  else cout<<"  --\n";
  if(s_remap[Key(Dipole::gqbar,Radiation::qtop)])
    s_remap[Key(Dipole::gqbar,Radiation::qtop)]->Which();
  else cout<<"  --\n";
  if(s_remap[Key(Dipole::gg,Radiation::qbarbot)])
    s_remap[Key(Dipole::gg,Radiation::qbarbot)]->Which();
  else cout<<"  --\n";
  if(s_remap[Key(Dipole::gg,Radiation::qtop)])
    s_remap[Key(Dipole::gg,Radiation::qtop)]->Which();
  else cout<<"  --\n";
  cout<<" -  -  II ig emit  -  -  -  -  -  -  - "<<endl;
  if(s_remap[Key(Dipole::iiqbarq,Radiation::qbarend)])
    s_remap[Key(Dipole::iiqbarq,Radiation::qbarend)]->Which();
  else cout<<"  --\n";
  if(s_remap[Key(Dipole::iiqbarq,Radiation::qfront)])
    s_remap[Key(Dipole::iiqbarq,Radiation::qfront)]->Which();
  else cout<<"  --\n";
  if(s_remap[Key(Dipole::iiqbarg,Radiation::qfront)])
    s_remap[Key(Dipole::iiqbarg,Radiation::qfront)]->Which();
  else cout<<"  --\n";
  if(s_remap[Key(Dipole::iigq,Radiation::qbarend)])
    s_remap[Key(Dipole::iigq,Radiation::qbarend)]->Which();
  else cout<<"  --\n";
  cout<<"======================================="<<endl;
}





const bool Dipole_Handler::AdjustCalcBox() {    //Static.

  static Calcbox box;    //Needed for a proper deletion at the very very end.

  //for(size_t i=0; i<box.v_psud.size(); ++i)cout<<i<<":"<<box.v_psud[i]<<"\n";
  //for(size_t i=0; i<box.v_prec.size(); ++i)cout<<i<<":"<<box.v_prec[i]<<"\n";

  if(s_count) {
    cerr<<"\nStatic method: "<<__PRETTY_FUNCTION__<<": "
	<<"Warning: Re-adjusting is not permitted "
	<<"since Dipole_Handler's are already present!\n"<<endl;
    return false;
  }

  const int             mode=dpa.kin.ShowerMode();
  const xbool           gspl=dpa.sud.GsplitRule();
  const Radiation::Type raty=dpa.sud.RadiationType();

  //Arrange the Sudakov_Group's.
  for(size_t i=0; i<box.v_psud.size(); ++i)
    if(box.v_psud[i]) { delete box.v_psud[i]; box.v_psud[i]=NULL;}
  box.v_psud.clear();
  if(mode & dsm::jff) {
    box.v_psud.push_back(new Sudakov_Group<Dipole::qqbar>(raty));
    box.v_psud.push_back(new Sudakov_Group<Dipole::qg>(raty));
    box.v_psud.push_back(new Sudakov_Group<Dipole::gqbar>(raty));
    box.v_psud.push_back(new Sudakov_Group<Dipole::gg>(raty));
    for(size_t i=0; i<4; ++i) assert(box.v_psud[i]);
  } else {
    for(size_t i=0; i<4; ++i) box.v_psud.push_back(NULL);
  }
  if(mode & dsm::jii) {
    box.v_psud.push_back(new IISudakov_Group<Dipole::iiqbarq>(raty));
    box.v_psud.push_back(new IISudakov_Group<Dipole::iiqbarg>(raty));
    box.v_psud.push_back(new IISudakov_Group<Dipole::iigq>(raty));
    box.v_psud.push_back(new IISudakov_Group<Dipole::iigg>(raty));
    for(size_t i=4; i<8; ++i) assert(box.v_psud[i]);
  } else {
    for(size_t i=4; i<8; ++i) box.v_psud.push_back(NULL);
  }
  //Fix the Sudakov map.
  s_sumap[Dipole::incorrect]=NULL;
  s_sumap[Dipole::qqbar]    =box.v_psud[0];
  s_sumap[Dipole::qg]       =box.v_psud[1];
  s_sumap[Dipole::gqbar]    =box.v_psud[2];
  s_sumap[Dipole::gg]       =box.v_psud[3];
  s_sumap[Dipole::iiqbarq]  =box.v_psud[4];
  s_sumap[Dipole::iiqbarg]  =box.v_psud[5];
  s_sumap[Dipole::iigq]     =box.v_psud[6];
  s_sumap[Dipole::iigg]     =box.v_psud[7];

  //Establish the overall recoil strategy right now and here.
  for(size_t i=0; i<box.v_prec.size(); ++i)
    if(box.v_prec[i]) { delete box.v_prec[i]; box.v_prec[i]=NULL;}
  box.v_prec.clear();
  //Which calculators are needed? (777...true, 999...false)
  std::vector<size_t> mpp(rl::stop,999);
  if(mode & dsm::jff) {
    mpp[rl::qag]=mpp[rl::qgg]=mpp[rl::gag]=mpp[rl::ggg]=777;
    mpp[rl::qga]=mpp[rl::gaq]=mpp[rl::gga]=mpp[rl::ggq]=777;
    if(raty<Radiation::g)
      mpp[rl::qag]=mpp[rl::qgg]=mpp[rl::gag]=mpp[rl::ggg]=999;
    else if(raty==Radiation::g)
      mpp[rl::qga]=mpp[rl::gaq]=mpp[rl::gga]=mpp[rl::ggq]=999;
    if(gspl==negative)
      mpp[rl::qga]=mpp[rl::gaq]=mpp[rl::gga]=mpp[rl::ggq]=999;
    else if(gspl==positive)
      mpp[rl::qag]=mpp[rl::qgg]=mpp[rl::gag]=mpp[rl::ggg]=999;
  }
  if(mode & dsm::jii) {
    mpp[rl::iiaqg]=mpp[rl::iiaqa]=mpp[rl::iiaqq]=777;
    mpp[rl::iiagg]=mpp[rl::iiagq]=777;
    mpp[rl::iigqg]=mpp[rl::iigqa]=777;
    mpp[rl::iiggg]=777;
    if(raty<Radiation::g)
      mpp[rl::iiaqg]=mpp[rl::iiagg]=mpp[rl::iigqg]=mpp[rl::iiggg]=999;
    else if(raty==Radiation::g)
      mpp[rl::iiaqa]=mpp[rl::iiaqq]=mpp[rl::iiagq]=mpp[rl::iigqa]=999;
    if(gspl==negative);
    else if(gspl==positive)
      mpp[rl::iiaqg]=mpp[rl::iiaqa]=mpp[rl::iiaqq]=
	mpp[rl::iiagg]=mpp[rl::iiagq]=
	mpp[rl::iigqg]=mpp[rl::iigqa]=
	mpp[rl::iiggg]=999;
  }
  std::vector<bool> gate;
  for(size_t i=0; ; ++i) {
    gate.push_back(false);
    for(size_t j=0; j<rl::stop; ++j) {
      if(mpp[j]==777) {
	if(dpa.kin.RecoilStrategy()[j]==Recoil_Strategy::stop) break;//(:o)//
	if(dpa.kin.RecoilStrategy()[j]==Recoil_Strategy::List[i]) {
	  gate.back()=true; mpp[j]=i;}
      }
    }
    ///cout<<"STRAT: "<<Recoil_Strategy::List[i]<<endl;
    if(Recoil_Strategy::List[i]==Recoil_Strategy::stop) break;
  }
  box.v_prec.resize(gate.size(),NULL);
  MakeRecos(gate,box.v_prec);
  ///cout<<gate.size()<<" == "<<box.v_prec.size()<<endl;
  ///for(size_t i=0; i<gate.size(); ++i) {
  ///  cout<<gate[i]<<" : "<<box.v_prec[i]<<" : ";
  ///  if(box.v_prec[i]) box.v_prec[i]->Which(); else cout<<"nil\n";}
  ///for(size_t j=0; j<mpp.size(); ++j) cout<<mpp[j]<<"\n";
  //Fix the Recoil map.
  if(mpp[rl::qag]>=777) s_remap[Key(Dipole::qqbar,Radiation::gluon)]=NULL;
  else s_remap[Key(Dipole::qqbar,Radiation::gluon)]=box.v_prec[mpp[rl::qag]];
  if(mpp[rl::qgg]>=777) s_remap[Key(Dipole::qg,Radiation::gluon)]   =NULL;
  else s_remap[Key(Dipole::qg,Radiation::gluon)]   =box.v_prec[mpp[rl::qgg]];
  if(mpp[rl::gag]>=777) s_remap[Key(Dipole::gqbar,Radiation::gluon)]=NULL;
  else s_remap[Key(Dipole::gqbar,Radiation::gluon)]=box.v_prec[mpp[rl::gag]];
  if(mpp[rl::ggg]>=777) s_remap[Key(Dipole::gg,Radiation::gluon)]   =NULL;
  else s_remap[Key(Dipole::gg,Radiation::gluon)]   =box.v_prec[mpp[rl::ggg]];
  //
  if(mpp[rl::qga]>=777) s_remap[Key(Dipole::qg,Radiation::qbarbot)] =NULL;
  else s_remap[Key(Dipole::qg,Radiation::qbarbot)] =box.v_prec[mpp[rl::qga]];
  if(mpp[rl::gaq]>=777) s_remap[Key(Dipole::gqbar,Radiation::qtop)] =NULL;
  else s_remap[Key(Dipole::gqbar,Radiation::qtop)] =box.v_prec[mpp[rl::gaq]];
  if(mpp[rl::gga]>=777) s_remap[Key(Dipole::gg,Radiation::qbarbot)] =NULL;
  else s_remap[Key(Dipole::gg,Radiation::qbarbot)] =box.v_prec[mpp[rl::gga]];
  if(mpp[rl::ggq]>=777) s_remap[Key(Dipole::gg,Radiation::qtop)]    =NULL;
  else s_remap[Key(Dipole::gg,Radiation::qtop)]    =box.v_prec[mpp[rl::ggq]];
  //
  if(mpp[rl::iiaqg]>=777) s_remap[Key(Dipole::iiqbarq,Radiation::gluon)]=NULL;
  else s_remap[Key(Dipole::iiqbarq,Radiation::gluon)]
	 =box.v_prec[mpp[rl::iiaqg]];
  if(mpp[rl::iiagg]>=777) s_remap[Key(Dipole::iiqbarg,Radiation::gluon)]=NULL;
  else s_remap[Key(Dipole::iiqbarg,Radiation::gluon)]
	 =box.v_prec[mpp[rl::iiagg]];
  if(mpp[rl::iigqg]>=777) s_remap[Key(Dipole::iigq,Radiation::gluon)]=NULL;
  else s_remap[Key(Dipole::iigq,Radiation::gluon)]=box.v_prec[mpp[rl::iigqg]];
  if(mpp[rl::iiggg]>=777) s_remap[Key(Dipole::iigg,Radiation::gluon)]=NULL;
  else s_remap[Key(Dipole::iigg,Radiation::gluon)]=box.v_prec[mpp[rl::iiggg]];
  //
  if(mpp[rl::iiaqa]>=777) s_remap[Key(Dipole::iiqbarq,Radiation::qbarend)]
			    =NULL;
  else s_remap[Key(Dipole::iiqbarq,Radiation::qbarend)]
	 =box.v_prec[mpp[rl::iiaqa]];
  if(mpp[rl::iiaqq]>=777) s_remap[Key(Dipole::iiqbarq,Radiation::qfront)]
			    =NULL;
  else s_remap[Key(Dipole::iiqbarq,Radiation::qfront)]
	 =box.v_prec[mpp[rl::iiaqq]];
  if(mpp[rl::iiagq]>=777) s_remap[Key(Dipole::iiqbarg,Radiation::qfront)]
			    =NULL;
  else s_remap[Key(Dipole::iiqbarg,Radiation::qfront)]
	 =box.v_prec[mpp[rl::iiagq]];
  if(mpp[rl::iigqa]>=777) s_remap[Key(Dipole::iigq,Radiation::qbarend)]
			    =NULL;
  else s_remap[Key(Dipole::iigq,Radiation::qbarend)]
	 =box.v_prec[mpp[rl::iigqa]];
  //Defined settings.
  s_remap[Key(Dipole::incorrect,Radiation::incorrect)] = NULL;
  s_remap[Key(Dipole::qqbar,Radiation::incorrect)]     = NULL;
  s_remap[Key(Dipole::qg,Radiation::incorrect)]        = NULL;
  s_remap[Key(Dipole::gqbar,Radiation::incorrect)]     = NULL;
  s_remap[Key(Dipole::gg,Radiation::incorrect)]        = NULL;
  s_remap[Key(Dipole::iiqbarq,Radiation::incorrect)]   = NULL;
  s_remap[Key(Dipole::iiqbarg,Radiation::incorrect)]   = NULL;
  s_remap[Key(Dipole::iigq,Radiation::incorrect)]      = NULL;
  s_remap[Key(Dipole::iigg,Radiation::incorrect)]      = NULL;
  //Check the recoil settings.
  for(Recoilbox::const_iterator it=s_remap.begin(); it!=s_remap.end(); ++it)
    if(it->second) it->second->TestKey(it->first);

#ifdef DIPOLE_PARAMETER_OUTPUT
  ListCalcBox();
  if(sf_1stadjust) cout<<"{ "<<__PRETTY_FUNCTION__<<" ... re-done }\n";
  else cout<<"{ "<<__PRETTY_FUNCTION__<<" ... done }\n";
  //abort();
#endif

  sf_1stadjust=true;

  return true;

}



//=============================================================================



void Dipole_Handler::DecoupleNew(Carrier& car) {
  assert(!(car.pDip || car.pGlu || car.pAqu || car.pQua));
  assert(car.Mup.GetItsVec().empty());
  car.pDip=p_dix; p_dix=NULL;
  car.pGlu=p_glu; p_glu=NULL;
  car.pAqu=p_ati; p_ati=NULL;
  car.pQua=p_ban; p_ban=NULL;
  car.DipOrder=f_below; f_below=false;
  car.RecoComp=m_rer.Poc; m_rer.Poc=both;
  if(m_rer.Mup.GetItsVec().empty());
  else car.Mup.Swap(m_rer.Mup);
}





void Dipole_Handler::RemoveNewProducts() {
  //Resets the news.
  f_below=false;
  m_rer.Reset(s_mpcode);
  //In case decoupling has not yet carried out, see Note1!
  if(p_dix) {
    delete p_dix; p_dix=NULL;
    if(p_glu) { delete p_glu; p_glu=NULL;}
    if(p_ban) { delete p_ban; p_ban=NULL;}
    if(p_ati) { delete p_ati; p_ati=NULL;}
  } else {
    if(p_ban) { delete p_ban; p_ban=NULL;}
    if(p_ati) { delete p_ati; p_ati=NULL;}
    p_glu=NULL;
  }
  assert(p_dix==NULL && p_glu==NULL && p_ban==NULL && p_ati==NULL);
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
  cout<<"Mimic Ariadne? "<<Sudakov_Calculator::Ariadne<<"."<<endl;
  cout<<"=========================================="<<endl;
}





void Dipole_Handler::ShowRecoil() const {
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





const bool Dipole_Handler::InduceDipoleRadiation(bool t, bool o, bool f) {

  static Histogram   histo(0,0,1000,200);//////////////////////////////////////
  static Xhisto      histu(40000);/////////////////////////////////////////////
  static Multidouble md(2,0.0);////////////////////////////////////////////////

  f_gate=0;

  if(!p_dip) return false;

  if(p_dip->Status()==On && p_dip->IsType()!=Dipole::incorrect);
  else { p_dip->SetEmitScale()=0.0; return false;}

  //Note1:
  //In case decoupling has not yet taken place, the following may lead the
  //loss of a Dipole_Particle of p_dip, however there is a mechanism to
  //compensate for this by creating a particle copy, which then physically
  //belongs to p_dip (this onsets the PointerHandling of p_dip).
  if(p_dix) {
    delete p_dix; p_dix=NULL;
    if(p_glu) { delete p_glu; p_glu=NULL;}
    if(p_ban) { delete p_ban; p_ban=NULL;}
    if(p_ati) { delete p_ati; p_ati=NULL;}
  } else {
    if(p_ban) { delete p_ban; p_ban=NULL;}
    if(p_ati) { delete p_ati; p_ati=NULL;}
    p_glu=NULL;
  }
  assert(p_dix==NULL && p_glu==NULL && p_ban==NULL && p_ati==NULL);

#ifdef STRICT_DIPOLE_HANDLER
  if(p_dip->PointerHandling()!=0) { p_dip->SetEmitScale()=0.0; return false;}
#endif

  //No testing of global parameters.
  //assert( p_dip->InvMass() > dpa.sud.MinK2t() );
  //assert( p_dip->ProdScale() > dpa.sud.MinK2t() );

  if(p_sudakov==NULL) { p_dip->SetEmitScale()=0.0; return false;}
  //assert(p_sudakov);

  //Kt Local-Analysis Sherpa-Analysis Comparison.
  //static size_t count=0;/////////////////////////////////////////////////////
  //static Histogram   histo(0,0,200,40);//////////////////////////////////////
  //++count;///////////////////////////////////////////////////////////////////
  //t=true;////////////////////////////////////////////////////////////////////
  //f=count==50000;////////////////////////////////////////////////////////////

  if( p_sudakov->GenerateVariablesFor(*p_dip,m_sur) ) {
    ///////////////////////////////////////////////////////////////////////////
    if(t) {
      if(o) m_sur.Print();
      md[0]=m_sur.X1;
      md[1]=m_sur.X3;
      histu.Insert(md);
      if(m_sur.Isr.size()==sr::stop) histo.Insert(m_sur.Isr[sr::kt]);
      else histo.Insert(0.0);
      if(f) {
	string name1("z_x1x3_test.dat");
	string name2("z_kt_test.dat");
	cout<<"Outputting "<<histu.Entries()
	    <<" Multidoubles to file "<<name1<<".\n";
	cout<<"Outputting the gluon kts to file "<<name2<<".\n";
	histu.Output(name1);
	histo.Output();
	histo.Finalize();
	histo.Output(name2);
      }
    }
    ///////////////////////////////////////////////////////////////////////////
    m_key.second=m_sur.Rad;
#ifdef DIPOLE_HANDLER_OUTPUT
    cout<<"(("<<m_key.first<<","<<m_key.second<<"))"<<endl;
#endif
    p_recoil=s_remap[m_key];
    p_dip->SetEmitScale()=m_sur.P2t;
    p_dip->SetFactScale()=dpa.evo.GetFactScaleFrom(m_sur);    //Suggestion.
    f_gate=p_dip->StateNumber;
  }
  else {
    ///////////////////////////////////////////////////////////////////////////
    if(t) {
      if(o) m_sur.Print();
      md[0]=m_sur.X1;
      md[1]=m_sur.X3;
      histu.Insert(md);
      histo.Insert(0.0);
      if(f) {
	string name1("z_x1x3_test.dat");
	string name2("z_kt_test.dat");
	cout<<"Outputting "<<histu.Entries()
	    <<" Multidoubles to file "<<name1<<".\n";
	cout<<"Outputting the gluon kts to file "<<name2<<".\n";
	histu.Output(name1);
	histo.Output();
	histo.Finalize();
	histo.Output(name2);
      }
    }
    ///////////////////////////////////////////////////////////////////////////
    if(p_dip->IsII()) p_dip->SetEmitScale()=dpa.sud.MinIIK2t();
    else p_dip->SetEmitScale()=dpa.sud.MinK2t();
    return false;
  }

#ifdef DIPOLE_HANDLER_OUTPUT
  cout<<"\ttransverse momentum and energy fractions:\n\t\t p2t=";
  cout<<m_sur.P2t<<endl;
  cout<<"\t\t x1="<<m_sur.X1<<endl;
  cout<<"\t\t x3="<<m_sur.X3<<endl;
#endif

  return true;

}





const bool Dipole_Handler::FinishDipoleRadiation() {

  assert(f_gate);    //Corresponds to: assert(p_dip) or no radiation.
  assert(f_gate==p_dip->StateNumber);

  if(m_sur.P2t!=p_dip->EmitScale() || this->Status()!='n') {
    f_gate=0; return false;
  }

  //assert((p_dip->TotP())[0] > 0.0);//////////////////////////////////////////

  m_rer.Reset(s_mpcode);
  m_rer.Vec.resize(3,Vec4D());
  m_rer.Vec[rr::p1]=p_dip->GetTopBranchPointer()->Momentum();
  m_rer.Vec[rr::p3]=p_dip->GetBotBranchPointer()->Momentum();

  if(dabs(m_rer.Vec[rr::p1].Abs2()) > 1.0e-9)
    cout<<" p1^2 ! -> "<<m_rer.Vec[rr::p1].Abs2()
	<<" \t"<<m_rer.Vec[rr::p1]<<endl;
  if(dabs(m_rer.Vec[rr::p3].Abs2()) > 1.0e-9)
    cout<<" p3^2 ! -> "<<m_rer.Vec[rr::p3].Abs2()
	<<" \t"<<m_rer.Vec[rr::p3]<<endl;

  assert(dabs(m_rer.Vec[rr::p1].Abs2()) < 1.0e-9);
  assert(dabs(m_rer.Vec[rr::p3].Abs2()) < 1.0e-9);
  assert(m_rer.Vec[rr::p1][0] > 0.0);
  assert(m_rer.Vec[rr::p3][0] > 0.0);

  assert(p_recoil);
  assert(p_recoil->GenerateMomenta(*p_dip,m_sur,m_rer));
  assert(GenerateSplitting());

  //Due to the dipole settings in GenerateSplitting,
  //actually the following should not be necessary.
  f_gate=0;

  return true;

}



//=============================================================================



const bool Dipole_Handler::GenerateSplitting() {

  //dpv.evo.SetFactScaleOffset(sqr(1.0));

  //So far, only the emission scale of p_dip has changed.
  //cout<<__PRETTY_FUNCTION__<<":\n"<<*p_dip<<endl;
  //m_sur.Print();

  switch(m_sur.Rad) {

  case Radiation::gluon: {
    p_dip->GetTopBranchPointer()->SetMomentum(m_rer.Vec[rr::p1]);
    p_dip->GetBotBranchPointer()->SetMomentum(m_rer.Vec[rr::p3]);
    //That updates the dipole as well as the neighbouring ones.
    p_glu=new Dipole::Glubranch(m_rer.Vec[rr::p2]); assert(p_glu);
    p_dix=new Dipole(*p_dip); assert(p_dix);    //Gives entry in copy memory.
    if(m_rer.Poc==front) {
      f_below=false;
      p_dix->RenewBranch(false,*p_glu);    //dixbot
      p_dip->RenewBranch(true,*p_glu);    //diptop
    } else {
      f_below=true;
      p_dip->RenewBranch(false,*p_glu);    //dipbot
      p_dix->RenewBranch(true,*p_glu);    //dixtop
    }
    p_dix->SetSource()=p_dip->Name;
    p_dix->SetSpinCorr()=p_dip->SpinCorr();
    p_dix->SetEvolScales(m_sur.P2t);
    p_dix->SetFactScale()=dpa.evo.GetFactScaleFrom(m_sur);
    break;
  }

  case Radiation::qtop: {
    f_below=false;
    p_dip->GetBotBranchPointer()->SetMomentum(m_rer.Vec[rr::p3]);
    //That updates the dipole and if existing the neighbouring one below.
    Dipole_Particle* topglu=p_dip->GetTopBranchPointer().operator->();
    assert(topglu->OrgType()==Nil);
    p_glu=static_cast<Dipole::Glubranch*>(topglu); assert(p_glu);
    p_ati=new Dipole::Antibranch(*m_sur.Sfc.Aqu,m_rer.Vec[rr::p1]);
    assert(p_ati);
    p_ban=new Dipole::Branch(*m_sur.Sfc.Qua,m_rer.Vec[rr::p2]);
    assert(p_ban);
    p_dip->RenewBranch(*p_ban);    //diptop
    break;
  }

  case Radiation::qbarbot: {
    f_below=true;
    p_dip->GetTopBranchPointer()->SetMomentum(m_rer.Vec[rr::p1]);
    //That updates the dipole and if existing the neighbouring one above.
    Dipole_Particle* botglu=p_dip->GetBotBranchPointer().operator->();
    assert(botglu->OrgType()==Nil);
    p_glu=static_cast<Dipole::Glubranch*>(botglu); assert(p_glu);
    p_ati=new Dipole::Antibranch(*m_sur.Sfc.Aqu,m_rer.Vec[rr::p2]);
    assert(p_ati);
    p_ban=new Dipole::Branch(*m_sur.Sfc.Qua,m_rer.Vec[rr::p3]);
    assert(p_ban);
    p_dip->RenewBranch(*p_ati);    //dipbot
    break;
  }

  case Radiation::qfront: {
    f_below=false;
    p_dip->GetBotBranchPointer()->SetMomentum(m_rer.Vec[rr::p3]);
    //That updates all dipoles carrying the initial quark.
    Dipole_Particle* top=p_dip->GetTopBranchPointer().operator->();
    assert(top->DipNum()==1);
    assert(top->OrgType()==Positive);
    p_ati=static_cast<Dipole::Antibranch*>(top); assert(p_ati);
    //Found initial antiquark!
    //Although being a Branch force its cast to an Antibranch!
    p_glu=new Dipole::Glubranch(m_rer.Vec[rr::p1],true); assert(p_glu);
    //Initial gluon!
    p_ban=new Dipole::Branch(*m_sur.Sfc.Qua,m_rer.Vec[rr::p2]);
    assert(p_ban);
    //Final quark!
    p_dip->RenewBranch(true,*p_glu);    //diptop
    p_dix=new Dipole(*p_ban,*p_glu); assert(p_dix);
    p_dix->SetSource()=p_dip->Name;
    p_dix->SetSpinCorr()=p_dip->SpinCorr();
    p_dix->SetEvolScales(m_sur.P2t);
    p_dix->SetFactScale()=dpa.evo.GetFactScaleFrom(m_sur);
    break;
  }

  case Radiation::qbarend: {
    f_below=true;
    p_dip->GetTopBranchPointer()->SetMomentum(m_rer.Vec[rr::p1]);
    //That updates all dipoles carrying the initial antiquark.
    Dipole_Particle* bot=p_dip->GetBotBranchPointer().operator->();
    assert(bot->DipNum()==1);
    assert(bot->OrgType()==Negative);
    p_ban=static_cast<Dipole::Branch*>(bot); assert(p_ban);
    //Found initial quark!
    //Although being an Antibranch force its cast to a Branch!
    p_glu=new Dipole::Glubranch(m_rer.Vec[rr::p3],true); assert(p_glu);
    //Initial gluon!
    p_ati=new Dipole::Antibranch(*m_sur.Sfc.Aqu,m_rer.Vec[rr::p2]);
    assert(p_ati);
    //Final antiquark!
    p_dip->RenewBranch(false,*p_glu);    //dipbot
    p_dix=new Dipole(*p_glu,*p_ati); assert(p_dix);
    p_dix->SetSource()=p_dip->Name;
    p_dix->SetSpinCorr()=p_dip->SpinCorr();
    p_dix->SetEvolScales(m_sur.P2t);
    p_dix->SetFactScale()=dpa.evo.GetFactScaleFrom(m_sur);
    break;
  }

  default:
    assert(m_sur.Rad==Radiation::gluon ||
	   m_sur.Rad==Radiation::qtop || m_sur.Rad==Radiation::qbarbot ||
	   m_sur.Rad==Radiation::qfront || m_sur.Rad==Radiation::qbarend);

  }

  p_dip->SetEvolScales(m_sur.P2t);
  p_dip->SetFactScale()=dpa.evo.GetFactScaleFrom(m_sur);    //Restatement.

  m_key.first=p_dip->IsType();
  m_key.second=Radiation::incorrect;
  p_sudakov=s_sumap[m_key.first];
  p_recoil=s_remap[m_key];

  //cout<<m_key.first<<endl;
  //cout<<m_key.second<<endl;
  //cout<<p_sudakov<<endl;
  //cout<<p_recoil<<endl;

  return true;

}



//=============================================================================



Dipole_Handler::Calcbox::Calcbox() : v_psud(), v_prec() {}





Dipole_Handler::Calcbox::~Calcbox() {
  for(vector<Sudakov_Calculator*>::iterator it=v_psud.begin();
      it!=v_psud.end(); ++it)
    if(*it) delete (*it);
  for(size_t i=0; i<v_prec.size(); ++i) cout<<" "<<i<<":"<<v_prec[i]<<"\n";////
  for(vector<Recoil_Calculator*>::iterator it=v_prec.begin();
      it!=v_prec.end(); ++it)
    if(*it) delete (*it);
  cout<<"~Calcbox done."<<endl;////////////////////////////////////////////////
}



//=============================================================================





//eof
