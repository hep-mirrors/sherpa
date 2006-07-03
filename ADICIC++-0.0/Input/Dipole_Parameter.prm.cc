//Globally defined parameter (static) variables of the ADICIC module.
//Preadjustments.





//AlphaS treatment flag
Sud::s_runalphas = true;

//Coupling
Sud::s_alphasfix = 0.12;

//Number of quark flavours
Sud::s_nffix = 3;

//Rule how to deal with gluon splittings:
//negative...all g-split procs are off,
//nil........no extra effect,
//positive...only and only g-split procs are allowed.
Sud::s_gsplit = nil;

//Radiation type to build up the Sudakov groups
Sud::s_radiatype = Radiation::g;

//FF dipoles, GeV^2
Sud::s_k2tmin = 1.0;
Sud::s_k2tmax = 8100.0;

//II dipoles, GeV^2 except for the factor
Sud::s_k2tiimin = 1.0;
Sud::s_k2tiifac = 1.0;
Sud::s_k2tiifixscale = 8100.0;
Sud::s_k2tiivarscale = -1.0;
Sud::s_k2tiimax = Sud::s_k2tiifac*Sud::s_k2tiifixscale;

//Efficiency factor to speed up the II dipole radiation calculation
Sud::s_iieffexp = 0.0;



//Dipole shower mode
Kin::s_dsmode = dsm::jff;
//Recoil strategies (compare with Recoil_Strategy.hpp)
//qqbar, qg, gqbar, gg radiating g
Kin::v_recostrat[rl::qag] = Recoil_Strategy::Kleiss;
Kin::v_recostrat[rl::qgg] = Recoil_Strategy::FixDir3;
Kin::v_recostrat[rl::gag] = Recoil_Strategy::FixDir1;
Kin::v_recostrat[rl::ggg] = Recoil_Strategy::Test;    //MinimizePt;
//iiqbarq, iiqbarg, iigq, iigg radiating g
Kin::v_recostrat[rl::iiaqg] = Recoil_Strategy::Unknown;
Kin::v_recostrat[rl::iiagg] = Recoil_Strategy::Unknown;
Kin::v_recostrat[rl::iigqg] = Recoil_Strategy::Unknown;
Kin::v_recostrat[rl::iiggg] = Recoil_Strategy::Unknown;
//qg radiating qbarbot, gqbar radiating qtop, gg radiating qbarbot, qtop
Kin::v_recostrat[rl::qga] = Recoil_Strategy::FixDir1;
Kin::v_recostrat[rl::gaq] = Recoil_Strategy::FixDir3;
Kin::v_recostrat[rl::gga] = Recoil_Strategy::FixDir1;
Kin::v_recostrat[rl::ggq] = Recoil_Strategy::FixDir3;
//iiqbarq radiating qbarend, qfront,
//iiqbarg radiating qfront, iigq radiating qbarend
Kin::v_recostrat[rl::iiaqa] = Recoil_Strategy::stop;
Kin::v_recostrat[rl::iiaqq] = Recoil_Strategy::stop;
Kin::v_recostrat[rl::iiagq] = Recoil_Strategy::stop;
Kin::v_recostrat[rl::iigqa] = Recoil_Strategy::stop;



//Factorization scale type (compare with Evolution_Strategy.hpp)
Evo::s_fascatype = fascat::p2t;
Evo::s_fmuf=Evo::s_fmur=1.0;
Evo::s_scaoffset = 0.0;
//Chain evolution strategy (compare with Evolution_Strategy.hpp)
Evo::v_chevostrat[cel::def] = Chain_Evolution_Strategy::Production;//Emission;
//Chain particle limit
Evo::s_chpartlim = 7777777;
//Chain correlation limit
Evo::s_chcorrlim = 7777777;
