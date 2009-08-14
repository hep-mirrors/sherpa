#include "MODEL/Interaction_Models/Interaction_Model_Higgs_SM.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include <stdio.h>


using namespace MODEL;
using namespace ATOOLS;
using namespace std;

Interaction_Model_Higgs_SM::Interaction_Model_Higgs_SM(MODEL::Model_Base * _model,
					   std::string _cplscheme,std::string _yukscheme) :
  Interaction_Model_Base("",_model,_cplscheme,_yukscheme)
{ 
  g1    = Kabbala(string("g_1"),
		  sqrt(4.*M_PI*ScalarFunction(std::string("alpha_QED"),rpa.gen.CplScale())));
  g2    = Kabbala(string("g_1/\\sin\\theta_W"), 
		  g1.Value()/sqrt(ScalarConstant(std::string("sin2_thetaW"))));
  sintW = Kabbala(std::string("\\sin\\theta_W"),
		  sqrt(ScalarConstant(std::string("sin2_thetaW"))));
  costW = Kabbala(std::string("\\cos\\theta_W"),
		  sqrt(1.-ScalarConstant(std::string("sin2_thetaW"))));
  PL    = Kabbala(string("P_L"),1.);
  PR    = Kabbala(string("P_R"),1.);
  M_I   = Kabbala(string("i"),Complex(0.,1.));
  root2 = Kabbala(string("\\sqrt{2}"),sqrt(2.));
  vev   = Kabbala(string("v_{Higgs_SM}"),ScalarConstant(std::string("vev")));
}

void Interaction_Model_Higgs_SM::c_FFS(std::vector<Single_Vertex>& vertex,int& vanz)
{
  Flavour flh(kf_h0);
  Kabbala kcpl0,kcpl1,M_h;
  if (!flh.IsOn()) return;

  for (short int i=1;i<17;i++) {
    if (i==7) i=11;
    Flavour flav = Flavour((kf_code)(i));
    if (flav.IsOn() && flav.IsFermion() && (flav.Yuk() > 0.)) {
      
      M_h = Kabbala(string("M_{")+flav.TexName()+string("}(m_h^2)"),
		    ScalarFunction(std::string("m")+std::string(flav.IDName()),sqr(flh.Mass(true))));

      kcpl0 = -M_I*M_h/vev;
      kcpl1 = kcpl0;
      
      if (!ATOOLS::IsZero(kcpl0.Value())) {
	vertex[vanz].in[0] = flav;
	vertex[vanz].in[1] = flh;
	vertex[vanz].in[2] = flav;

	vertex[vanz].cpl[0]  = kcpl0;
	vertex[vanz].cpl[1]  = kcpl1;
	vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();

	
	if (flav.Strong()) {
	  vertex[vanz].Color.push_back(Color_Function(cf::D));     
	  vertex[vanz].Color.back().SetParticleArg(0,2);     
	  vertex[vanz].Color.back().SetStringArg('0','2');     
	}
	else 
	  vertex[vanz].Color.push_back(Color_Function(cf::None));
	
	
	vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("FFS",LF_Key()));

	vertex[vanz].on      = 1;
	vertex.push_back(Single_Vertex());vanz++;
      }
    }
  }
}


void Interaction_Model_Higgs_SM::c_VVS(std::vector<Single_Vertex>& vertex,int& vanz)
{
  Flavour flh(kf_h0);
  if (!flh.IsOn()) return;

  Kabbala kcpl0,kcpl1;  
  Kabbala num_2 = Kabbala(string("2"),2.);    
  Flavour flWplus(kf_Wplus);
  // W h W
  if (flWplus.IsOn()) {
    vertex[vanz].in[0] = flWplus.Bar();
    vertex[vanz].in[1] = flh;
    vertex[vanz].in[2] = flWplus.Bar();
    
    kcpl0 = M_I*g2*flWplus.Yuk();
    kcpl1 = kcpl0;

    vertex[vanz].cpl[0]  = kcpl0;
    vertex[vanz].cpl[1]  = vertex[vanz].cpl[0];
    vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();

    
    vertex[vanz].Color.push_back(Color_Function(cf::None));     

    
    vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("Gab",LF_Key()));     
    vertex[vanz].Lorentz.back()->SetParticleArg(0,2);     

    vertex[vanz].on      = 1;
    vertex.push_back(Single_Vertex());vanz++;
  }

  Flavour flav = Flavour(kf_Z);
  // Z h Z
  if (flav.IsOn()) {
    vertex[vanz].in[0] = flav;
    vertex[vanz].in[1] = flh;
    vertex[vanz].in[2] = flav;
    
    kcpl0 = M_I*g2*flav.Yuk()/costW;
    kcpl1 = kcpl0;

    vertex[vanz].cpl[0]  = kcpl0;
    vertex[vanz].cpl[1]  = kcpl1;
    vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();

    
    vertex[vanz].Color.push_back(Color_Function(cf::None));     

    
    vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("Gab",LF_Key()));  
    vertex[vanz].Lorentz.back()->SetParticleArg(0,2);     

    vertex[vanz].on      = 1;
    vertex.push_back(Single_Vertex());vanz++;
  }
}


void Interaction_Model_Higgs_SM::c_SSS(std::vector<Single_Vertex>& vertex,int& vanz)
{
  Flavour flh = Flavour(kf_h0);
  Kabbala kcpl0,kcpl1,yuk;  
  Kabbala num_3 = Kabbala(string("3"),3.);  

  if (flh.IsOn()) {  
    vertex[vanz].in[0] = flh;
    vertex[vanz].in[1] = flh;
    vertex[vanz].in[2] = flh;

    yuk   = Kabbala(string("M_{")+flh.TexName()+string("}"),flh.Yuk());
    kcpl0 = -M_I*yuk*yuk*(num_3/vev);
    kcpl1 = kcpl0;
    
    vertex[vanz].cpl[0]  = kcpl0;
    vertex[vanz].cpl[1]  = kcpl1;
    vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();

    
    vertex[vanz].Color.push_back(Color_Function(cf::None));     
    
    
    vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("SSS",LF_Key()));     

    vertex[vanz].on      = 1;
    vertex.push_back(Single_Vertex());vanz++;
  }
}

void Interaction_Model_Higgs_SM::c_SSSS(std::vector<Single_Vertex>& vertex,int& vanz)
{
  Flavour flh = Flavour(kf_h0);
  Kabbala kcpl0,kcpl1,yuk;  
  Kabbala num_3 = Kabbala(string("3"),3.);  

  if (flh.IsOn()) {  
    vertex[vanz].in[0] = flh;
    vertex[vanz].in[1] = flh;
    vertex[vanz].in[2] = flh;
    vertex[vanz].in[3] = flh;

    vertex[vanz].nleg  = 4;  
    
    yuk   = Kabbala(string("M_{")+flh.TexName()+string("}"),flh.Yuk());
    kcpl0 = -M_I*yuk*yuk*(num_3/(vev*vev));
    kcpl1 = kcpl0;
    
    vertex[vanz].cpl[0]  = kcpl0;
    vertex[vanz].cpl[1]  = kcpl1;
    vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();

    
    vertex[vanz].Color.push_back(Color_Function(cf::None));     
    
    
    vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("SSSS",LF_Key()));     

    vertex[vanz].on      = 1;
    vertex.push_back(Single_Vertex());vanz++;
  }
}

void Interaction_Model_Higgs_SM::c_SSVV(std::vector<Single_Vertex>& vertex,int& vanz)
{
  Kabbala num_2 = Kabbala(string("2"),2.);  

  Flavour flavWplus(kf_Wplus);
  Flavour flavZ(kf_Z);
  Flavour flavh(kf_h0);
  Kabbala kcpl0,kcpl1;
  
  // h - Z - Z - h  
  if (flavZ.IsOn() && flavh.IsOn()) {
    vertex[vanz].in[0] = flavZ;
    vertex[vanz].in[1] = flavh;
    vertex[vanz].in[2] = flavh;
    vertex[vanz].in[3] = flavZ;
    
    vertex[vanz].nleg     = 4;
    
    kcpl0 = M_I*g2*g2/(costW*costW*num_2);
    kcpl1 = kcpl0;
    
    vertex[vanz].cpl[0]  = kcpl0;
    vertex[vanz].cpl[1]  = kcpl1;
    vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
    
    
    vertex[vanz].Color.push_back(Color_Function(cf::None));     
    
    
    vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("VVSS",LF_Key()));
    vertex[vanz].Lorentz.back()->SetParticleArg(0,3);     
    
    vertex[vanz].on      = 1;
    vertex.push_back(Single_Vertex());vanz++;
  }

  // h - W - W - h  
  if (flavWplus.IsOn() && flavh.IsOn()) {
    vertex[vanz].in[0] = flavWplus.Bar();
    vertex[vanz].in[1] = flavh;
    vertex[vanz].in[2] = flavh;
    vertex[vanz].in[3] = flavWplus.Bar();
    
    vertex[vanz].nleg     = 4;
    
    kcpl0 = M_I*g2*g2/num_2;
    kcpl1 = kcpl0;
    
    vertex[vanz].cpl[0]  = kcpl0;
    vertex[vanz].cpl[1]  = kcpl1;
    vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
    
    
    vertex[vanz].Color.push_back(Color_Function(cf::None));     
    
    
    vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("VVSS",LF_Key()));     
    vertex[vanz].Lorentz.back()->SetParticleArg(0,3);     
    
    vertex[vanz].on      = 1;
    vertex.push_back(Single_Vertex());vanz++;
  }
}
