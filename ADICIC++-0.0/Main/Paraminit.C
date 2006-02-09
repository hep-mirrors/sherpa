//bof
//Version: 3 ADICIC++-0.0/2005/09/13

//Implementation of Paraminit.H.



#include "Run_Parameter.H"
#include "Dipole_Parameter.H"
#include "Sudakov_Calculator.H"
#include "Dipole_Handler.H"
#include "Paraminit.H"





using namespace std;
using namespace ATOOLS;
using namespace ADICIC;





//=============================================================================



const bool Dipole_Parameter_Init::Do() {    //Static.

#ifdef PARAMINIT_OUTPUT
  cout<<"{ Executing ... "<<__PRETTY_FUNCTION__<<" ...\n";
  Dipole_Parameter::Show();
#endif

  //###########################################################################
  //AlphaS treatment flag
  Dipole_Parameter::Sud::s_runalphas=true;//false;//true;
  //Coupling
  Dipole_Parameter::Sud::s_alphasfix=0.15;//0.12;
  //Number of quark flavours
  Dipole_Parameter::Sud::s_nffix=2;
  //Radiation type to build up the Sudakov groups
  Dipole_Parameter::Sud::s_radiatype=Radiation::gduscb;
  //GeV^2
  Dipole_Parameter::Sud::s_k2tmin=0.64;//1.0;
  Dipole_Parameter::Sud::s_k2tmax=8100.0;
  //GeV^2
  Dipole_Parameter::Sud::s_k2tiimin=1.0;//0.25;
  Dipole_Parameter::Sud::s_k2tiifac=12.005;//4.0;
  Dipole_Parameter::Sud::s_k2tiifixscale=8100.0;
  //II dipole efficiency factor
  Dipole_Parameter::Sud::s_iieffexp=1.0;
  //Dipole shower mode
  Dipole_Parameter::Kin::s_dsmode = dsm::iiff;
  //Recoil strategies (compare with Recoil_Strategy.hpp)
  //    qqbar, qg, gqbar, gg radiating g
  Dipole_Parameter::Kin::v_recostrat[rl::qag] = Recoil_Strategy::Kleiss;
  Dipole_Parameter::Kin::v_recostrat[rl::qgg] = Recoil_Strategy::FixDir3;
  Dipole_Parameter::Kin::v_recostrat[rl::gag] = Recoil_Strategy::FixDir1;
  Dipole_Parameter::Kin::v_recostrat[rl::ggg] = Recoil_Strategy::MinimizePt;
  //    iiqbarq, iiqbarg, iigq, iigg radiating g
  Dipole_Parameter::Kin::v_recostrat[rl::iiaqg] = Recoil_Strategy::Ktii;
  Dipole_Parameter::Kin::v_recostrat[rl::iiagg] = Recoil_Strategy::Ktii;
  Dipole_Parameter::Kin::v_recostrat[rl::iigqg] = Recoil_Strategy::Ktii;
  Dipole_Parameter::Kin::v_recostrat[rl::iiggg] = Recoil_Strategy::Ktii;
  //    qg radiating qbarbot, gqbar radiating qtop, gg radiating qbarbot, qtop
  Dipole_Parameter::Kin::v_recostrat[rl::qga] = Recoil_Strategy::FixDir1;
  Dipole_Parameter::Kin::v_recostrat[rl::gaq] = Recoil_Strategy::FixDir3;
  Dipole_Parameter::Kin::v_recostrat[rl::gga] = Recoil_Strategy::FixDir1;
  Dipole_Parameter::Kin::v_recostrat[rl::ggq] = Recoil_Strategy::FixDir3;
  //    iiqbarq radiating qbarend, qfront,
  //    iiqbarg radiating qfront, iigq radiating qbarend
  Dipole_Parameter::Kin::v_recostrat[rl::iiaqa] = Recoil_Strategy::Ktii;//stop;
  Dipole_Parameter::Kin::v_recostrat[rl::iiaqq] = Recoil_Strategy::Ktii;//stop;
  Dipole_Parameter::Kin::v_recostrat[rl::iiagq] = Recoil_Strategy::Ktii;//stop;
  Dipole_Parameter::Kin::v_recostrat[rl::iigqa] = Recoil_Strategy::Ktii;//Unknown;
  //Chain evolution strategy (compare with Evolution_Strategy.hpp)
  Dipole_Parameter::Evo::v_chevostrat[cel::def] =
    Chain_Evolution_Strategy::Production;
  //###########################################################################

  //Re-calculate the maximum k2t for II dipoles.
  dpv.sud.SetMaxIIScale(Dipole_Parameter::Sud::s_k2tiivarscale);

  //Dipole_Parameter::SetWithStatics();    //Test.

#ifdef PARAMINIT_OUTPUT
  Dipole_Parameter::Show();
#endif

  extern Run_Parameter ATOOLS::rpa;
  Dipole_Parameter::Check(rpa.gen.Ecms());
#ifdef PARAMINIT_OUTPUT
  cout<<"\n    Ecms="<<rpa.gen.Ecms()<<"\n\n";
#endif

  Sudakov_Calculator::AdjustEnvironment();
  Dipole_Handler::AdjustCalcBox();

#ifdef PARAMINIT_OUTPUT
  Sudakov_Calculator::ShowEnvironment();
  Dipole_Handler::ShowCalcBox();
  cout<<"}\n"<<endl;
#endif

  return true;

}





Dipole_Parameter_Init::Dipole_Parameter_Init() {
  static bool init=Do();
  assert(init);
}



//=============================================================================





//eof
