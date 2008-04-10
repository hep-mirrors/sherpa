#include "Interaction_Model_EHC.H"
#include "MathTools.H"
#include "Message.H"
#include "Run_Parameter.H"
#include <stdio.h>


using namespace MODEL;
using namespace ATOOLS;
using namespace std;

Interaction_Model_EHC::Interaction_Model_EHC(MODEL::Model_Base * _model,
					   std::string _cplscheme,std::string _yukscheme) :
  Interaction_Model_Base("",_model,_cplscheme,_yukscheme)
{ 
  double Ecms2 = sqr(rpa.gen.Ecms());
  double hmass2 = sqr(Flavour(kf_h0).Mass());

  geff  = Kabbala(std::string("I_S"),ScalarConstant(std::string("HIGGS_GG_EFF")));
  ghgg  = Kabbala(std::string("ghgg"),ScalarConstant(std::string("Higgs_gg_fac"))*
		  ScalarFunction(std::string("alpha_S"),hmass2)/(2.*M_PI)/ScalarConstant(std::string("vev")));
  g1    = Kabbala(string("g_1"),
		  sqrt(4.*M_PI*ScalarFunction(std::string("alpha_QED"),Ecms2)));
  g2    = Kabbala(string("g_1/\\sin\\theta_W"), 
		  g1.Value()/sqrt(ScalarConstant(std::string("sin2_thetaW"))));
  g3  = Kabbala(string("g_3"),sqrt(4.*M_PI*ScalarFunction(std::string("alpha_S"),Ecms2)));
  sintW = Kabbala(std::string("\\sin\\theta_W"),
		  sqrt(ScalarConstant(std::string("sin2_thetaW"))));
  costW = Kabbala(std::string("\\cos\\theta_W"),
		  sqrt(1.-ScalarConstant(std::string("sin2_thetaW"))));
  PL    = Kabbala(string("P_L"),1.);
  PR    = Kabbala(string("P_R"),1.);
  M_I   = Kabbala(string("i"),Complex(0.,1.));
  root2 = Kabbala(string("\\sqrt{2}"),sqrt(2.));
  vev   = Kabbala(string("v_{EW}"),ScalarConstant(std::string("vev")));
}



void Interaction_Model_EHC::c_VVS(std::vector<Single_Vertex>& vertex,int& vanz)
{
  Flavour flh(kf_h0);
  Kabbala kcpl0,kcpl1;  
  Kabbala num_2 = Kabbala(string("2"),2.);  
 
 if (!flh.IsOn()) return;
  
  Flavour flav(kf_photon);
  // Photon h Photon
  if (flav.IsOn()) {
    vertex[vanz].in[0] = flav;
    vertex[vanz].in[1] = flh;
    vertex[vanz].in[2] = flav;
    
    kcpl0 = M_I*geff;
    kcpl1 = kcpl0;

    vertex[vanz].cpl[0]  = kcpl0;
    vertex[vanz].cpl[1]  = kcpl0;
    vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();

    
    vertex[vanz].Color.push_back(Color_Function(cf::None));;     

    
    vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("Triangle",LF_Key()));     
    vertex[vanz].Lorentz.back()->SetParticleArg(0,2);     

    vertex[vanz].on      = 1;
    vertex.push_back(Single_Vertex());vanz++;
  }

  Flavour flg(kf_gluon);
  // Gluon h Gluon
  if (flg.IsOn()) {
    vertex[vanz].in[0] = flg;
    vertex[vanz].in[1] = flh;
    vertex[vanz].in[2] = flg;
    
    kcpl0 = M_I*ghgg;
    kcpl1 = kcpl0;

    vertex[vanz].cpl[0]  = kcpl0;
    vertex[vanz].cpl[1]  = kcpl0;
    vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();

    
    vertex[vanz].Color.push_back(Color_Function(cf::G));;     
    vertex[vanz].Color.back().SetParticleArg(0,2);     
    vertex[vanz].Color.back().SetStringArg('0','2');     

    
    vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("Triangle",LF_Key()));     
    vertex[vanz].Lorentz.back()->SetParticleArg(0,2);     

    vertex[vanz].on      = 1;
    vertex.push_back(Single_Vertex());vanz++;
  }
  Flavour flsh(kf_shgluon);
  // gluon h shgluon
  if (flg.IsOn() && flsh.IsOn()) {
    vertex[vanz].in[2] = flg;
    vertex[vanz].in[0] = flsh;
    vertex[vanz].in[1] = flh;
    
    kcpl0 = M_I;
    kcpl1 = kcpl0;

    vertex[vanz].cpl[0]  = kcpl0;
    vertex[vanz].cpl[1]  = kcpl0;
    vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();

    
    vertex[vanz].Color.push_back(Color_Function(cf::G));;     
    vertex[vanz].Color.back().SetParticleArg(0,2);     
    vertex[vanz].Color.back().SetStringArg('0','2');     

    
    vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("C4GS",LF_Key()));     
    vertex[vanz].Lorentz.back()->SetParticleArg(0,2);     

    vertex[vanz].on      = 1;
    vertex[vanz].t       = -1;
    vertex.push_back(Single_Vertex());vanz++;
  }

}
 


void Interaction_Model_EHC::c_VVVV(std::vector<Single_Vertex>& vertex,int& vanz)
{
}
