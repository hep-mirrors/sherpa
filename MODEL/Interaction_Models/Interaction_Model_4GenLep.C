#include "MODEL/Interaction_Models/Interaction_Model_4GenLep.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include <stdio.h>


using namespace MODEL;
using namespace ATOOLS;
using namespace std;

DECLARE_GETTER(Interaction_Model_4GenLep_Getter,"SM+4thLF",
	       Interaction_Model_Base,Interaction_Model_Arguments);

Interaction_Model_Base *Interaction_Model_4GenLep_Getter::
operator()(const Interaction_Model_Arguments &args) const
{
  return new Interaction_Model_4GenLep
    (args.p_model,args.m_cplscheme,args.m_yukscheme);
}

void Interaction_Model_4GenLep_Getter::PrintInfo
(std::ostream &str,const size_t width) const
{ 
  str<<"The Standard Model + 4th generation leptons"; 
}

Interaction_Model_4GenLep::Interaction_Model_4GenLep(MODEL::Model_Base * _model,
					       std::string _cplscheme,std::string _yukscheme) :
  Interaction_Model_Base("FOURTH_GEN_LEPTONS",_model,_cplscheme,_yukscheme)
{ 
  p_mosm    = new Interaction_Model_SM(p_model,_cplscheme,_yukscheme); 

  g1    = Kabbala(string("g_1"),
		  sqrt(4.*M_PI*ScalarFunction(std::string("alpha_QED"),rpa->gen.CplScale())));
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
  vev   = Kabbala(string("v_{EW}"),ScalarConstant(std::string("vev")));
}

void Interaction_Model_4GenLep::c_FFV(std::vector<Single_Vertex>& vertex,int& vanz)
{
  p_mosm->c_FFV(vertex,vanz);

  Flavour flphoton(kf_photon);
  Flavour flZ(kf_Z);
  Flavour flWplus(kf_Wplus);
  for (short int i=17;i<19;i++) {
    Flavour flav1               = Flavour((kf_code)(i));
    Kabbala charge1             = Kabbala(string("Q_{")+flav1.TexName()+string("}"),flav1.Charge());
    Kabbala isoweak1            = Kabbala(string("T_{")+flav1.TexName()+string("}"),flav1.IsoWeak());

    Kabbala kcpl0,kcpl1;    
    if (flav1.IsOn()) {
      for (short int j=i;j<19;j++) {
	Flavour flav2           = Flavour((kf_code)(j));
	Kabbala charge2         = Kabbala(string("Q_{")+flav2.TexName()+string("}"),flav2.Charge());
	Kabbala isoweak2        = Kabbala(string("T_{")+ flav2.TexName()+string("}"),flav2.IsoWeak());	
	
	if (flav2.IsOn()) {
	  if (flav1==flav2) {
	    //photon
	    if (flphoton.IsOn()) {
	      kcpl0             = -g1*M_I*charge1;
	      kcpl1             = kcpl0;
	      if (!ATOOLS::IsZero(kcpl0.Value())) {
		vertex[vanz].in[0]   = flav1;
		vertex[vanz].in[1]   = Flavour(kf_photon);
		vertex[vanz].in[2]   = flav2;
		vertex[vanz].cpl[0]  = kcpl0;
		vertex[vanz].cpl[1]  = vertex[vanz].cpl[0];
		vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
	
		
		if (flav1.Strong()) {
		  vertex[vanz].Color.push_back(Color_Function(cf::D));     
		  vertex[vanz].Color.back().SetParticleArg(0,2);     
		  vertex[vanz].Color.back().SetStringArg('0','2');     
		}
		else 
		  vertex[vanz].Color.push_back(Color_Function(cf::None));

		
		vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("Gamma",LF_Key()));
		vertex[vanz].Lorentz.back()->SetParticleArg(1);     

		vertex[vanz].on      = 1;
		vertex.push_back(Single_Vertex());vanz++;
	      }
	    }
	    //Z
	    if (flZ.IsOn()) {
	      
	      kcpl0             = M_I/costW*charge1*sintW*sintW*g2;
	      kcpl1             = -M_I/costW*(isoweak1-charge1*sintW*sintW)*g2;
	     
	      vertex[vanz].in[0]     = flav1;
	      vertex[vanz].in[1]     = Flavour(kf_Z);
	      vertex[vanz].in[2]     = flav2;
	      vertex[vanz].cpl[0]    = kcpl0;
	      vertex[vanz].cpl[1]    = kcpl1;
	      vertex[vanz].Str       = (kcpl0*PR+kcpl1*PL).String();
	
	      
	      if (flav1.Strong()) {
		vertex[vanz].Color.push_back(Color_Function(cf::D));     
		vertex[vanz].Color.back().SetParticleArg(0,2);     
		vertex[vanz].Color.back().SetStringArg('0','2');     
	      }
	      else 
		  vertex[vanz].Color.push_back(Color_Function(cf::None));

	      
	      vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("Gamma",LF_Key()));
	      vertex[vanz].Lorentz.back()->SetParticleArg(1);     

	      vertex[vanz].on     = 1;
	      vertex.push_back(Single_Vertex());vanz++;
	    }
	  }
	  //W
	  if (flWplus.IsOn()) {
	    short int hit = 1;
	    Kabbala kcpl0,kcpl1;
	    kcpl0 = Kabbala(string("zero"),0.);
	    kcpl1 = Kabbala(string("1"),0.);

	    if (!((flav1.IsDowntype() && flav2.IsUptype()) ||
                  (flav2.IsDowntype() && flav1.IsUptype()))) hit = 0;
	    if ((flav1.IsLepton() && !flav2.IsLepton()) ||
		(flav1.IsQuark() && !flav2.IsQuark()) ) hit = 0;
	    if (hit==1) {
	      if (flav1.IsDowntype() && i>10 && j==i+1) 
		kcpl1 = -M_I/root2*g2;
	      if (!ATOOLS::IsZero(kcpl1.Value())) {
		vertex[vanz].in[1] = flWplus.Bar();
		if (flav1.IsDowntype()) {
		  vertex[vanz].in[0] = flav1;
		  vertex[vanz].in[2] = flav2;
		}
		else {
		  vertex[vanz].in[0] = flav2;
		  vertex[vanz].in[2] = flav1;
		}
		
		vertex[vanz].cpl[0]  = kcpl0;
		vertex[vanz].cpl[1]  = kcpl1;
		vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
	
		
		if (flav1.Strong()) {
		  vertex[vanz].Color.push_back(Color_Function(cf::D));     
		  vertex[vanz].Color.back().SetParticleArg(0,2);     
		  vertex[vanz].Color.back().SetStringArg('0','2');     
		}
		else 
		  vertex[vanz].Color.push_back(Color_Function(cf::None));

		
		vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("Gamma",LF_Key()));
		vertex[vanz].Lorentz.back()->SetParticleArg(1);     

		vertex[vanz].on      = 1;
		vertex.push_back(Single_Vertex());vanz++;
	      }
	    }
	  }
	}
      }
    }
  }
}

void Interaction_Model_4GenLep::c_VVV(std::vector<Single_Vertex>& vertex,int& vanz)
{
  p_mosm->c_VVV(vertex,vanz);
}
void Interaction_Model_4GenLep::c_VVVV(std::vector<Single_Vertex>& vertex,int& vanz)
{
  p_mosm->c_VVVV(vertex,vanz);
}

void Interaction_Model_4GenLep::c_FFS(std::vector<Single_Vertex>& vertex,int& vanz)  
{ 
  p_mosm->c_FFS(vertex,vanz);
  Flavour flh(kf_h0);
  Kabbala kcpl0,kcpl1,M_h;
  if (!flh.IsOn()) return;

  for (short int i=17;i<19;i++) {
    Flavour flav = Flavour((kf_code)(i));
    if (flav.IsOn() && flav.IsFermion() && (flav.Yuk() > 0.)) {
      
      M_h = Kabbala(string("M_{")+flav.TexName()+string("}"),flav.Mass(true));

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

void Interaction_Model_4GenLep::c_VVS(std::vector<Single_Vertex>& vertex,int& vanz)  
{ 
  p_mosm->c_VVS(vertex,vanz);
}

void Interaction_Model_4GenLep::c_SSS(std::vector<Single_Vertex>& vertex,int& vanz)  
{ 
  p_mosm->c_SSS(vertex,vanz);
}

void Interaction_Model_4GenLep::c_SSV(std::vector<Single_Vertex>& vertex,int& vanz) 
{ 
  p_mosm->c_SSV(vertex,vanz);
}

void Interaction_Model_4GenLep::c_SSVV(std::vector<Single_Vertex>& vertex,int& vanz) 
{ 
  p_mosm->c_SSVV(vertex,vanz); 
}

void Interaction_Model_4GenLep::c_SSSS(std::vector<Single_Vertex>& vertex,int& vanz) 
{ 
  p_mosm->c_SSSS(vertex,vanz);
}

Interaction_Model_4GenLep::~Interaction_Model_4GenLep()
{
  delete p_mosm;
}
