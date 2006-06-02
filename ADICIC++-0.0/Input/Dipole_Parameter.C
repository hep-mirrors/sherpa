//bof
//Version: 3 ADICIC++-0.0/2005/09/13

//Implementation of Dipole_Parameter.H.



#include <iostream>
#include "Evolution_Strategy.hpp"
#include "Dipole_Parameter.H"
#include "MathTools.H"




using namespace std;
using namespace ADICIC;





//#include "..."





//=============================================================================



bool            Dipole_Parameter::Sud::s_runalphas=false;
double          Dipole_Parameter::Sud::s_alphasfix=0.0;
unsigned short  Dipole_Parameter::Sud::s_nffix=0;
Radiation::Type Dipole_Parameter::Sud::s_radiatype=Radiation::g;
double          Dipole_Parameter::Sud::s_k2tmin=0.0;
double          Dipole_Parameter::Sud::s_k2tmax=0.0;
double          Dipole_Parameter::Sud::s_k2tiimin=0.0;
double          Dipole_Parameter::Sud::s_k2tiimax=0.0;
double          Dipole_Parameter::Sud::s_iieffexp=0.0;
//--------------------------------------------------------------
double          Dipole_Parameter::Sud::s_k2tiifac=0.0;
double          Dipole_Parameter::Sud::s_k2tiifixscale=0.0;
double          Dipole_Parameter::Sud::s_k2tiivarscale=0.0;



int Dipole_Parameter::Kin::s_dsmode=dsm::off;
//Explicit setup.
vector<Recoil_Strategy::Type> Dipole_Parameter::Kin::v_recostrat
=vector<Recoil_Strategy::Type>(rl::stop,Recoil_Strategy::stop);



//Explicit setup.
vector<Chain_Evolution_Strategy::Type> Dipole_Parameter::Evo::v_chevostrat
=vector<Chain_Evolution_Strategy::Type>(cel::stop,
					Chain_Evolution_Strategy::stop);
//------------------------------------------------------------------------
size_t Dipole_Parameter::Evo::s_chpartlim=0;
size_t Dipole_Parameter::Evo::s_chcorrlim=0;





//Initialize with meaningful values.
const bool Dipole_Parameter::sf_start=Dipole_Parameter::SetWithStatics();



//=============================================================================



Dipole_Parameter::Dipole_Parameter() {}

Dipole_Parameter::~Dipole_Parameter() {}



//=============================================================================



void Dipole_Parameter::Show() {    //Static.
  cout<<endl;
  cout<<"======================================================="<<endl;
  cout<<"         Current ADICIC parameter adjustments."<<endl;
  cout<<"-------------------------------------------------------"<<endl;
  cout<<"Enable running AlphaS   = "<<Sud::s_runalphas<<".\n";
  cout<<"Value for fixed AlphaS  = "<<Sud::s_alphasfix<<".\n";
  cout<<"Value for fixed Nf      = "<<Sud::s_nffix<<".\n";
  cout<<"Overall radiation type  = "<<Sud::s_radiatype<<".\n";
  cout<<"FF dipole shower cut-off scale = "<<Sud::s_k2tmin<<" GeV^2"
      <<"\t ["<<sqrt(Sud::s_k2tmin)<<"].\n";
  cout<<"FF dipole shower maximum scale = "<<Sud::s_k2tmax<<" GeV^2"
      <<"\t ["<<sqrt(Sud::s_k2tmax)<<"].\n";
  cout<<"II dipole shower cut-off scale = "<<Sud::s_k2tiimin<<" GeV^2"
      <<"\t ["<<sqrt(Sud::s_k2tiimin)<<"].\n";
  cout<<"II dipole shower maximum scale = "<<Sud::s_k2tiimax<<" GeV^2"
      <<"\t ["<<sqrt(Sud::s_k2tiimax)<<"]"
      <<"  ("<<Sud::s_k2tiifac
      <<", "<<Sud::s_k2tiifixscale
      <<", "<<Sud::s_k2tiivarscale<<").\n";
  cout<<"II dipole efficiency factor = "<<Sud::s_iieffexp<<".\n";
  cout<<"-------------------------------------------------------"<<endl;
  cout<<"Dipole shower mode = "<<Kin::s_dsmode<<".\n";
  cout<<"- - - - - - - - - - - - - - - - - - - - - - - - - - - -"<<endl;
  cout<<"Recoil strategy for q-qbar dipoles radiating g    = "
      <<Kin::v_recostrat[rl::qag]<<".\n";
  cout<<"Recoil strategy for q-g dipoles radiating g       = "
      <<Kin::v_recostrat[rl::qgg]<<".\n";
  cout<<"Recoil strategy for g-qbar dipoles radiating g    = "
      <<Kin::v_recostrat[rl::gag]<<".\n";
  cout<<"Recoil strategy for g-g dipoles radiating g       = "
      <<Kin::v_recostrat[rl::ggg]<<".\n";
  cout<<"Recoil strategy for ii qbar-q dipoles radiating g       = "
      <<Kin::v_recostrat[rl::iiaqg]<<".\n";
  cout<<"Recoil strategy for ii qbar-g dipoles radiating g       = "
      <<Kin::v_recostrat[rl::iiagg]<<".\n";
  cout<<"Recoil strategy for ii g-q dipoles radiating g          = "
      <<Kin::v_recostrat[rl::iigqg]<<".\n";
  cout<<"Recoil strategy for ii g-g dipoles radiating g          = "
      <<Kin::v_recostrat[rl::iiggg]<<".\n";
  cout<<"Recoil strategy for q-g dipoles radiating qbarbot = "
      <<Kin::v_recostrat[rl::qga]<<".\n";
  cout<<"Recoil strategy for g-qbar dipoles radiating qtop = "
      <<Kin::v_recostrat[rl::gaq]<<".\n";
  cout<<"Recoil strategy for g-g dipoles radiating qbarbot = "
      <<Kin::v_recostrat[rl::gga]<<".\n";
  cout<<"Recoil strategy for g-g dipoles radiating qtop    = "
      <<Kin::v_recostrat[rl::ggq]<<".\n";
  cout<<"Recoil strategy for ii qbar-q dipoles radiating qbarend = "
      <<Kin::v_recostrat[rl::iiaqa]<<".\n";
  cout<<"Recoil strategy for ii qbar-q dipoles radiating qfront  = "
      <<Kin::v_recostrat[rl::iiaqq]<<".\n";
  cout<<"Recoil strategy for ii qbar-g dipoles radiating qfront  = "
      <<Kin::v_recostrat[rl::iiagq]<<".\n";
  cout<<"Recoil strategy for ii g-q dipoles radiating qbarend    = "
      <<Kin::v_recostrat[rl::iigqa]<<".\n";
  cout<<"------------------------------------------------------"<<endl;
  cout<<"Chain evolution strategy = "
      <<Evo::v_chevostrat[cel::def]<<".\n";
  cout<<"Chain particle limit     = "
      <<Evo::s_chpartlim<<".\n";
  cout<<"Chain correlation limit  = "
      <<Evo::s_chcorrlim<<".\n";
  cout<<"======================================================"<<endl;
}





const bool Dipole_Parameter::Check(const double cmsen) {    //Static.

  if(0.001<Sud::s_alphasfix && Sud::s_alphasfix<1.0); else {
    cerr<<"\nParameter out of range: value of fixed AlphaS!\n";
    assert(0.001<Sud::s_alphasfix && Sud::s_alphasfix<1.0);
  }
  if(0.0<Sud::s_k2tmin && Sud::s_k2tmin<Sud::s_k2tmax); else {
    cerr<<"\nParameter out of range: FF dipole shower cut-off scale!\n";
    assert(0.0<Sud::s_k2tmin && Sud::s_k2tmin<Sud::s_k2tmax);
  }
  //Wait till Adicic.dat
  //if(Sud::s_k2tmax<=sqr(cmsen)); else {
  //  cerr<<"\nParameter out of range: FF dipole shower maximum scale!\n";
  //  assert(Sud::s_k2tmax<=sqr(cmsen));
  //}
  if(0.0<Sud::s_k2tiimin && Sud::s_k2tiimin<Sud::s_k2tiimax); else {
    cerr<<"\nParameter out of range: II dipole shower cut-off scale!\n";
    assert(0.0<Sud::s_k2tiimin && Sud::s_k2tiimin<Sud::s_k2tiimax);
  }
  //if(Sud::s_k2tiimax<=sqr(cmsen)/4); else {
  //  cerr<<"\nParameter out of range: II dipole shower maximum scale!\n";
  //  assert(Sud::s_k2tiimax<=sqr(cmsen)/4);
  //}
  if(Sud::s_iieffexp>=0.0); else {
    cerr<<"\nParameter out of range: II dipole efficiency exponent!\n";
    assert(Sud::s_iieffexp>=0.0);
  }


  if(Evo::v_chevostrat[cel::def]>=0); else {
    cerr<<"\nParameter out of range: chain evolution strategy!\n";
    assert(Evo::v_chevostrat[1]>=0);
  }
  if(Evo::v_chevostrat[cel::def]!=Chain_Evolution_Strategy::stop); else {
    cerr<<"\nParameter is unspecified: chain evolution strategy!\n";
    assert(Evo::v_chevostrat[1]!=Chain_Evolution_Strategy::stop);
  }

  //Surely more to check later.

  return true;

}



//=============================================================================



const bool Dipole_Parameter::SetWithStatics() {    //Static.
#ifdef DIPOLE_PARAMETER_OUTPUT
  //  cout<<"{ Executing ... "<<__PRETTY_FUNCTION__<<" }\n";
#endif
#include "Dipole_Parameter.prm.cc"
  return Check(2*sqrt(Sud::s_k2tmax));
}



//=============================================================================



const Dipole_Parameter ADICIC::dpa;    //"Static."
      Dipole_Parameter ADICIC::dpv;



//=============================================================================





//eof
