//Globally defined parameter (static) variables of the ADICIC module.
//Preadjustments.





//AlphaS treatment flag
s_isalphasrun = true;

//coupling
s_alphasfix = 0.12;

//GeV^2
s_k2tmin = 1.0;
s_k2tmax = 8100.0;

//recoil strategies (compare with Recoil_Strategy.hpp)
s_restratqqbar = Recoil_Strategy::Label_qgqbar;
s_restratqg    = Recoil_Strategy::Label_qgg;
s_restratgqbar = Recoil_Strategy::Label_ggqbar;
s_restratgg    = Recoil_Strategy::Label_ggg;

//chain evolution strategy (compare with Evolution_Strategy.hpp)
s_chevolstrat = Chain_Evolution_Strategy::Label;




