#include "MODEL/Interaction_Models/Interaction_Model_QCD_Grav.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Run_Parameter.H"

using namespace MODEL;
using namespace ATOOLS;
using namespace std;

Interaction_Model_QCD_Grav::Interaction_Model_QCD_Grav(MODEL::Model_Base * _model,
					     std::string _cplscheme,std::string _yukscheme) :
  Interaction_Model_Base("",_model,_cplscheme,_yukscheme)
{
  g3    = Kabbala(string("g_3"),sqrt(4.*M_PI*ScalarFunction(std::string("alpha_S"),rpa.gen.CplScale())));
  PL    = Kabbala(string("P_L"),1.);
  PR    = Kabbala(string("P_R"),1.);
  M_I   = Kabbala(string("i"),Complex(0.,1.)); 
  kap   = Kabbala(string("kappa"),ScalarConstant(string("kappa")));
  om    = Kabbala(string("omega"),ScalarConstant(string("omega")));
  num2  = Kabbala(string("2"),2.);
  num15 = Kabbala(string("1.5"),1.5);
}

void Interaction_Model_QCD_Grav::c_VVT(std::vector<Single_Vertex>& vertex,int& vanz)
{
  Kabbala kcpl0,kcpl1; 
  
  Flavour flgraviton(kf_graviton);
  Flavour flgs(kf_gscalar);
  Flavour flav(kf_gluon);
  if(flgraviton.IsOn()){

    vertex[vanz].in[0] = flav;
    vertex[vanz].in[1] = flgraviton;
    vertex[vanz].in[2] = flav;
    
    kcpl0 = -M_I*kap;
    kcpl1 = Kabbala();
    
    vertex[vanz].cpl[0]  = kcpl0;
    vertex[vanz].cpl[1]  = kcpl1;
    vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();

    
    vertex[vanz].Color.push_back(Color_Function(cf::G));//GD;     
    vertex[vanz].Color.back().SetParticleArg(0,2);     
    vertex[vanz].Color.back().SetStringArg('0','2');     
    
    
    vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("VVT",LF_Key()));     
    vertex[vanz].Lorentz.back()->SetParticleArg(0,2,1);     

    vertex[vanz].on      = 1;
    vertex.push_back(Single_Vertex());vanz++;
  }

  // gluon gscalar gluon
  if (flav.IsOn()) {
    vertex[vanz].in[0] = flav;
    vertex[vanz].in[1] = flgs;
    vertex[vanz].in[2] = flav;
    
    kcpl0 = M_I*om*kap;
    kcpl1 = Kabbala();
  
    vertex[vanz].cpl[0]  = kcpl0;
    vertex[vanz].cpl[1]  = kcpl1;
    vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();

    
    vertex[vanz].Color.push_back(Color_Function(cf::G));     
    vertex[vanz].Color.back().SetParticleArg(0,2);     
    vertex[vanz].Color.back().SetStringArg('0','2');     
    
    
    vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("VVGS",LF_Key()));     
    vertex[vanz].Lorentz.back()->SetParticleArg(0,2,1);     

    vertex[vanz].on      = 1;
    vertex.push_back(Single_Vertex());vanz++;
  }
}

void Interaction_Model_QCD_Grav::c_FFVT(std::vector<Single_Vertex>& vertex,int& vanz)
{
  
  Flavour flgraviton(kf_graviton);
  Flavour flgs(kf_gscalar);
  Kabbala kcpl0,kcpl1;
  
  for (short int i=1;i<7;i++) {
    Flavour flav = Flavour((kf_code)(i));
    if (flav.Strong() && flav.IsOn()) {
      if(flgraviton.IsOn()) { 
	vertex[vanz].nleg    = 4;
	vertex[vanz].in[0] = flav;
	vertex[vanz].in[1] = Flavour(kf_gluon);
	vertex[vanz].in[2] = flav;
	vertex[vanz].in[3] = flgraviton;
	
	kcpl0 = g3*M_I*kap/num2;
	kcpl1 = kcpl0;
	
	vertex[vanz].cpl[0]  = kcpl0;
	vertex[vanz].cpl[1]  = kcpl1;
	vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
	
	
	vertex[vanz].Color.push_back(Color_Function(cf::T));     
	vertex[vanz].Color.back().SetParticleArg(1,2,0);     
	vertex[vanz].Color.back().SetStringArg('1','2','0');     
	
	
	vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("FFVT",LF_Key()));     
	vertex[vanz].Lorentz.back()->SetParticleArg(1,3);     
	
	vertex[vanz].on      = 1;
	vertex.push_back(Single_Vertex());vanz++;
      }
      if(flgs.IsOn()){
	vertex[vanz].nleg    = 4;
	vertex[vanz].in[0] = flav;
	vertex[vanz].in[1] = Flavour(kf_gluon);
	vertex[vanz].in[2] = flav;
	vertex[vanz].in[3] = flgs;
	
	kcpl0 = -g3*M_I*om*kap*num15;
	kcpl1 = kcpl0;
	
	vertex[vanz].cpl[0]  = kcpl0;
	vertex[vanz].cpl[1]  = kcpl1;
	vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
	
	
	vertex[vanz].Color.push_back(Color_Function(cf::T));     
	vertex[vanz].Color.back().SetParticleArg(1,2,0);     
	vertex[vanz].Color.back().SetStringArg('1','2','0');     

	
	vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("FFVGS",LF_Key()));     
	vertex[vanz].Lorentz.back()->SetParticleArg(1,3);     
	
	vertex[vanz].on      = 1;
	vertex.push_back(Single_Vertex());vanz++;
      }
    } 
  }
}

void Interaction_Model_QCD_Grav::c_VVVT(std::vector<Single_Vertex>& vertex,int& vanz)
{
  Kabbala kcpl0,kcpl1; 
  
  Flavour flgraviton(kf_graviton);
  if(!flgraviton.IsOn())return;

  vertex[vanz].nleg    = 4;
  for (short int i=0;i<3;i++)
    vertex[vanz].in[i] = Flavour(kf_gluon);
  vertex[vanz].in[3] = flgraviton;
  
  kcpl0 = g3*kap; 
  kcpl1 = kcpl0; 

  vertex[vanz].cpl[0]  = kcpl0;
  vertex[vanz].cpl[1]  = kcpl1;
  vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();

  
  vertex[vanz].Color.push_back(Color_Function(cf::F));     
  vertex[vanz].Color.back().SetParticleArg(0,2,1);     
  vertex[vanz].Color.back().SetStringArg('0','2','1');     

  
  vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("VVVT",LF_Key()));     
  vertex[vanz].Lorentz.back()->SetParticleArg(0,1,2,3);     

  vertex[vanz].on      = 1;
  vertex.push_back(Single_Vertex());vanz++;
}

















