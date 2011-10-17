#include "MODEL/Interaction_Models/Interaction_Model_EW_Grav.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Run_Parameter.H"
//#include "ATOOLS/Math/Vector.H"
#include <stdio.h>

using namespace MODEL;
using namespace ATOOLS;
using namespace std;

Interaction_Model_EW_Grav::Interaction_Model_EW_Grav(MODEL::Model_Base * _model,
					   std::string _cplscheme,std::string _yukscheme) :
  Interaction_Model_Base("",_model,_cplscheme,_yukscheme)
{ 
  g1    = Kabbala(string("g_1"),
		  sqrt(4.*M_PI*ScalarFunction(std::string("alpha_QED"),rpa->gen.CplScale())));
  if(ScalarNumber(std::string("WidthScheme"))==0){
    g2       = Kabbala(string("g_1/\\sin\\theta_W"), 
		     g1.Value()/sqrt(ScalarConstant(std::string("sin2_thetaW"))));
    sintW    = Kabbala(std::string("\\sin\\theta_W"),
		     sqrt(ScalarConstant(std::string("sin2_thetaW"))));
    costW    = Kabbala(std::string("\\cos\\theta_W"),
		     sqrt(1.-ScalarConstant(std::string("sin2_thetaW"))));
    vev      = Kabbala(string("v_{EW}"),ScalarConstant(std::string("vev")));
  }else{
    g2       = Kabbala(string("g_1/\\sin\\theta_W"), 
		     g1.Value()/sqrt(ComplexConstant(std::string("csin2_thetaW"))));
    sintW    = Kabbala(std::string("\\sin\\theta_W"),
		     sqrt(ComplexConstant(std::string("csin2_thetaW"))));
    costW    = Kabbala(std::string("\\cos\\theta_W"),
		     sqrt(1.-ComplexConstant(std::string("csin2_thetaW"))));
    vev      = Kabbala(string("v_{EW}"),ComplexConstant(std::string("cvev")));
  }

  kap   = Kabbala(string("kappa"),ScalarConstant(string("kappa")));
  om    = Kabbala(string("omega"),ScalarConstant(string("omega")));
  
  PL    = Kabbala(string("P_L"),1.);
  PR    = Kabbala(string("P_R"),1.);
  M_I   = Kabbala(string("i"),Complex(0.,1.));
  root2 = Kabbala(string("\\sqrt{2}"),sqrt(2.));
  num2  = Kabbala(string("2"),2.);
  num4  = Kabbala(string("4"),4.);
  num15 = Kabbala(string("1.5"),1.5);

}

void Interaction_Model_EW_Grav::c_FFT(std::vector<Single_Vertex>& vertex,int& vanz)
{
  Flavour flgraviton(kf_graviton);
  Flavour flgs(kf_gscalar);
  for (short int i=1;i<17;i++) {
    if (i==7) i=11;
    Flavour flav1 = Flavour((kf_code)(i));
	
    if (flav1.IsOn()) {
      for (short int j=i;j<17;j++) {
	if (j==7) j=11;	
	Flavour flav2 = Flavour((kf_code)(j));
	
	if (flav2.IsOn()) {
	  if (flav1==flav2) {
	    
	    if (flgraviton.IsOn()) {
	      Kabbala mf;
	      if(ScalarNumber(std::string("WidthScheme"))==0){
	        mf = Kabbala(string("M_{")+flav1.TexName()+string("}"),flav1.Mass());
	      }else{
	        mf = Kabbala(string("M_{")+flav1.TexName()+string("}"),
	                           sqrt(sqr(flav1.Mass())-Complex(0.,1.)*flav1.Width()*flav1.Mass()));
	      }
	      Kabbala kcpl0,kcpl1;
	      kcpl0 = -M_I*kap/num4;
	      kcpl1 = num2*mf;
	      
	      vertex[vanz].in[0] = flav1;
	      vertex[vanz].in[1] = Flavour(kf_graviton);
	      vertex[vanz].in[2] = flav2;
	      vertex[vanz].cpl[0]  = kcpl0;
	      vertex[vanz].cpl[1]  = vertex[vanz].cpl[0];
	      vertex[vanz].Str     = (kcpl0*PR+kcpl0*PL).String();
	      vertex[vanz].cpl[2]  = kcpl1;
	      
	      

	      if (flav1.Strong()) {  
		vertex[vanz].Color.push_back(Color_Function(cf::D));     
		vertex[vanz].Color.back().SetParticleArg(0,2);     
		vertex[vanz].Color.back().SetStringArg('0','2');     
	      }
	      else 
		vertex[vanz].Color.push_back(Color_Function(cf::None));
	      
	      
	      vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("FFT",LF_Key()));     
	      vertex[vanz].Lorentz.back()->SetParticleArg(1);     

	      vertex[vanz].on      = 1;
	      vertex.push_back(Single_Vertex());vanz++;	      
	    }
	    //scalar graviton mode
	    if (flgs.IsOn()) {
	      Kabbala mf;
	      if(ScalarNumber(std::string("WidthScheme"))==0){
	        mf = Kabbala(string("M_{")+flav1.TexName()+string("}"),flav1.Mass());
	      }else{
	        mf = Kabbala(string("M_{")+flav1.TexName()+string("}"),
	                           sqrt(sqr(flav1.Mass())-Complex(0.,1.)*flav1.Width()*flav1.Mass()));
	      }
	      Kabbala kcpl0,kcpl1;
	      kcpl0 = M_I*om*kap*Kabbala(string("3/4"),.75);
	      kcpl1 = mf*Kabbala(string("8/3"),8./3.);
	      
	      vertex[vanz].in[0] = flav1;
	      vertex[vanz].in[1] = Flavour(kf_gscalar);
	      vertex[vanz].in[2] = flav2;
	      vertex[vanz].cpl[0] = kcpl0;
	      vertex[vanz].cpl[1] = kcpl0;
	      vertex[vanz].Str     = (kcpl0*PR+kcpl0*PL).String();
	      vertex[vanz].cpl[2] = kcpl1; 
	
	      
	      if (flav1.Strong()) {  
		vertex[vanz].Color.push_back(Color_Function(cf::D));     
		vertex[vanz].Color.back().SetParticleArg(0,2);     
		vertex[vanz].Color.back().SetStringArg('0','2');     
	      }
	      else 
		vertex[vanz].Color.push_back(Color_Function(cf::None));
	      
	      
	      vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("FFGS",LF_Key()));     
	      vertex[vanz].on     = 1;
	      vertex.push_back(Single_Vertex());vanz++;
	    }
	  }   
     	}
      }
    }
  }
}


void Interaction_Model_EW_Grav::c_FFVT(std::vector<Single_Vertex>& vertex,int& vanz)
{
  //return;
  Flavour flphoton(kf_photon);
  Flavour flZ(kf_Z);
  Flavour flWplus(kf_Wplus);
  Flavour flgraviton(kf_graviton);
  Flavour flgs(kf_gscalar);
  for (short int i=1;i<17;i++) {
    if (i==7) i=11;
    Flavour flav1 = Flavour((kf_code)(i));
    Kabbala charge1 = Kabbala(string("Q_{")+flav1.TexName()+string("}"),flav1.Charge());
    Kabbala isoweak1 = Kabbala(string("T_{")+flav1.TexName()+string("}"),flav1.IsoWeak());
    
    if (flav1.IsOn()) {
      for (short int j=i;j<17;j++) {
	if (j==7) j=11;	
	Flavour flav2 = Flavour((kf_code)(j));
	Kabbala charge2 = Kabbala(string("Q_{")+flav2.TexName()+string("}"),flav2.Charge());
	Kabbala isoweak2 = Kabbala(string("T_{")+ flav2.TexName()+string("}"),flav2.IsoWeak());	
	      
	if (flav2.IsOn()) {
	  if (flav1==flav2) {
	    //photon + graviton
	    if (flphoton.IsOn()&&flgraviton.IsOn()) {
	      
	      Kabbala kcpl0,kcpl1;
	      kcpl0 = g1*M_I*charge1*kap/num2;
	      kcpl1 = kcpl0;
	      
	      if (!ATOOLS::IsZero(charge1.Value())) {
		vertex[vanz].nleg     = 4;
		vertex[vanz].in[0] = flav1;
		vertex[vanz].in[1] = Flavour(kf_photon);
		vertex[vanz].in[2] = flav2;
		vertex[vanz].in[3] = flgraviton;
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
		
		
		vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("FFVT",LF_Key()));     
		vertex[vanz].Lorentz.back()->SetParticleArg(1,3);     

		vertex[vanz].on      = 1;
		vertex.push_back(Single_Vertex());vanz++;
	      }
	    }
	    //Z +graviton
	    if (flZ.IsOn()&&flgraviton.IsOn()) {
	      
	      Kabbala kcpl0,kcpl1;
	      kcpl0 = -M_I/costW*charge1*sintW*sintW*g2*kap/num2;
	      kcpl1 = M_I/costW*(isoweak1-charge1*sintW*sintW)*g2*kap/num2;
	      
	      vertex[vanz].nleg     = 4;
	      vertex[vanz].in[0] = flav1;
	      vertex[vanz].in[1] = Flavour(kf_Z);
	      vertex[vanz].in[2] = flav2;
	      vertex[vanz].in[3] = flgraviton;
	      vertex[vanz].cpl[0] = kcpl0;
	      vertex[vanz].cpl[1] = kcpl1;
	      vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
	
	      
	      if (flav1.Strong()) {
		vertex[vanz].Color.push_back(Color_Function(cf::D));     
		vertex[vanz].Color.back().SetParticleArg(0,2);     
		vertex[vanz].Color.back().SetStringArg('0','2');     
	      }
	      else 
		vertex[vanz].Color.push_back(Color_Function(cf::None));
		
	      
	      vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("FFVT",LF_Key()));     
	      vertex[vanz].Lorentz.back()->SetParticleArg(1,3);     

	      vertex[vanz].on     = 1;
	      vertex.push_back(Single_Vertex());vanz++;
	    }
	    //photon + gscalar
	    if (flphoton.IsOn()&&flgs.IsOn()) {
	      
	      Kabbala kcpl0,kcpl1;
	      kcpl0 = -g1*M_I*charge1*om*kap*num15;
	      kcpl1 = kcpl0;
	      
	      if (!ATOOLS::IsZero(charge1.Value())) {
		vertex[vanz].nleg     = 4;
		vertex[vanz].in[0] = flav1;
		vertex[vanz].in[1] = Flavour(kf_photon);
		vertex[vanz].in[2] = flav2;
		vertex[vanz].in[3] = flgs;
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
		
		
		vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("FFVGS",LF_Key()));     
		vertex[vanz].Lorentz.back()->SetParticleArg(1);     

		vertex[vanz].on      = 1;
		vertex.push_back(Single_Vertex());vanz++;
	      }
	    }
	    //Z +gscalar
	    if (flZ.IsOn()&&flgs.IsOn()) {
	      
	      Kabbala kcpl0,kcpl1;
	      kcpl0 = M_I/costW*charge1*sintW*sintW*g2*om*kap*num15;
	      kcpl1 = -M_I/costW*(isoweak1-charge1*sintW*sintW)*g2*om*kap*num15;
	      
	      vertex[vanz].nleg     = 4;
	      vertex[vanz].in[0] = flav1;
	      vertex[vanz].in[1] = Flavour(kf_Z);
	      vertex[vanz].in[2] = flav2;
	      vertex[vanz].in[3] = flgs;
	      vertex[vanz].cpl[0] = kcpl0;
	      vertex[vanz].cpl[1] = kcpl1;
	      vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
	
	      
	      if (flav1.Strong()) {
		vertex[vanz].Color.push_back(Color_Function(cf::D));     
		vertex[vanz].Color.back().SetParticleArg(0,2);     
		vertex[vanz].Color.back().SetStringArg('0','2');     
	      }
	      else 
		vertex[vanz].Color.push_back(Color_Function(cf::None));
	      
	      
	      vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("FFVGS",LF_Key()));     
	      vertex[vanz].Lorentz.back()->SetParticleArg(1);     

	      vertex[vanz].on     = 1;
	      vertex.push_back(Single_Vertex());vanz++;
	    }
	  }
	  //W + graviton
	  if (flWplus.IsOn()&&flgraviton.IsOn()) {
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
		kcpl1 = M_I/root2*g2*kap/num2;
	      if (i<7 && j<7) {
		if (flav1.IsDowntype())
		kcpl1 = M_I/root2*g2*K_CKM((i-1)/2,j/2-1)*kap/num2;
		else 	    
		kcpl1 = M_I/root2*g2*K_CKM(i/2-1,(j-1)/2)*kap/num2;		
	      }
	      if (!ATOOLS::IsZero(kcpl1.Value()/kap.Value())) {
		vertex[vanz].nleg     = 4;
		vertex[vanz].in[1] = flWplus.Bar();
		if (flav1.IsDowntype()) {
		  vertex[vanz].in[0] = flav1;
		  vertex[vanz].in[2] = flav2;
		}
		else {
		  vertex[vanz].in[0] = flav2;
		  vertex[vanz].in[2] = flav1;
		}
		vertex[vanz].in[3] = flgraviton;
		
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
	      
		
		vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("FFVT",LF_Key()));     
		vertex[vanz].Lorentz.back()->SetParticleArg(1,3);     

		vertex[vanz].on      = 1;
		vertex.push_back(Single_Vertex());vanz++;
	      }
	    }
	  }
	  //W + gscalar
	  if (flWplus.IsOn()&&flgs.IsOn()) {
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
		kcpl1 = -M_I/root2*g2*om*kap*num15;
	      if (i<7 && j<7) {
		if (flav1.IsDowntype())
		kcpl1 = -M_I/root2*g2*K_CKM((i-1)/2,j/2-1)*om*kap*num15;
		else 	    
		kcpl1 = -M_I/root2*g2*K_CKM(i/2-1,(j-1)/2)*om*kap*num15;		
	      }
	      if (!ATOOLS::IsZero(kcpl1.Value()/kap.Value())) {
		vertex[vanz].nleg     = 4;
		vertex[vanz].in[1] = flWplus.Bar();
		if (flav1.IsDowntype()) {
		  vertex[vanz].in[0] = flav1;
		  vertex[vanz].in[2] = flav2;
		}
		else {
		  vertex[vanz].in[0] = flav2;
		  vertex[vanz].in[2] = flav1;
		}
		vertex[vanz].in[3] = flgs;
		
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
	      
		
		vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("FFVGS",LF_Key()));     
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

void Interaction_Model_EW_Grav::c_VVT(std::vector<Single_Vertex>& vertex,int& vanz)
{
  Flavour flgraviton(kf_graviton);
  Flavour flgs(kf_gscalar);
  Kabbala kcpl0,kcpl1;  
  Flavour flav(kf_Wplus);

  if (flgraviton.IsOn()){  
    // W graviton W
    if (flav.IsOn()) {
      vertex[vanz].in[0] = flav;
      vertex[vanz].in[1] = flgraviton;
      vertex[vanz].in[2] = flav;
      
      kcpl0 = -M_I*kap;
      if(ScalarNumber(std::string("WidthScheme"))==0){
        kcpl1 = Kabbala(string("\\sqr(M_{")+flav.TexName()+string("})"),sqr(flav.Mass()));
      }else{
        kcpl1 = Kabbala(string("\\sqr(M_{")+flav.TexName()+string("})"),
                        sqr(flav.Mass())-Complex(0.,1.)*flav.Width()*flav.Mass());
      }
      
      vertex[vanz].cpl[0]  = kcpl0;
      vertex[vanz].cpl[1]  = kcpl1;
      vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();

      
      vertex[vanz].Color.push_back(Color_Function(cf::None));     
      
      
      vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("VVT",LF_Key()));     
      vertex[vanz].Lorentz.back()->SetParticleArg(0,2,1);     

      vertex[vanz].on      = 1;
      vertex.push_back(Single_Vertex());vanz++;
    }

    flav = Flavour(kf_Z);
    // Z graviton Z
    if (flav.IsOn()) {
      vertex[vanz].in[0] = flav;
      vertex[vanz].in[1] = flgraviton;
      vertex[vanz].in[2] = flav;
      
      kcpl0 = -M_I*kap;
      if(ScalarNumber(std::string("WidthScheme"))==0){
        kcpl1 = Kabbala(string("\\sqr(M_{")+flav.TexName()+string("})"),sqr(flav.Mass()));
      }else{
        kcpl1 = Kabbala(string("\\sqr(M_{")+flav.TexName()+string("})"),
                        sqr(flav.Mass())-Complex(0.,1.)*flav.Width()*flav.Mass());
      }
      vertex[vanz].cpl[0]  = kcpl0;
      vertex[vanz].cpl[1]  = kcpl1;
      vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();

      
      vertex[vanz].Color.push_back(Color_Function(cf::None));     
      
      
      vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("VVT",LF_Key()));     
      vertex[vanz].Lorentz.back()->SetParticleArg(0,2,1);     

      vertex[vanz].on      = 1;
      vertex.push_back(Single_Vertex());vanz++;
    }
    flav = Flavour(kf_photon);
    // photon graviton photon
    if (flav.IsOn()) {
      vertex[vanz].in[0] = flav;
      vertex[vanz].in[1] = flgraviton;
      vertex[vanz].in[2] = flav;
    
      kcpl0 = -M_I*kap;
      kcpl1 = Kabbala("",0);

      vertex[vanz].cpl[0]  = kcpl0;
      vertex[vanz].cpl[1]  = kcpl1;
      vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
      
      
      vertex[vanz].Color.push_back(Color_Function(cf::None));     

      
      vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("VVT",LF_Key()));     
      vertex[vanz].Lorentz.back()->SetParticleArg(0,2,1);     
      
      vertex[vanz].on      = 1;
      vertex.push_back(Single_Vertex());vanz++;
    }
  }
  if (!flgs.IsOn()) return;
  
  flav = Flavour(kf_Wplus);
  // W gscalar W
  if (flav.IsOn()) {
    vertex[vanz].in[0] = flav;
    vertex[vanz].in[1] = flgs;
    vertex[vanz].in[2] = flav;
    Kabbala ma;
    
    if(ScalarNumber(std::string("WidthScheme"))==0){
        ma = Kabbala(string("\\sqr(M_{")+flav.TexName()+string("})"),sqr(flav.Mass()));
      }else{
        ma = Kabbala(string("\\sqr(M_{")+flav.TexName()+string("})"),
                        sqr(flav.Mass())-Complex(0.,1.)*flav.Width()*flav.Mass());
      }
    kcpl0 = M_I*om*kap;
    kcpl1 = ma;
    
    vertex[vanz].cpl[0]  = kcpl0;
    vertex[vanz].cpl[1]  = kcpl1;
    vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();

    
    vertex[vanz].Color.push_back(Color_Function(cf::None));     

    
    vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("VVGS",LF_Key()));     
    vertex[vanz].Lorentz.back()->SetParticleArg(0,2,1);     

    vertex[vanz].on      = 1;
    vertex.push_back(Single_Vertex());vanz++;
  }

  flav = Flavour(kf_Z);
  // Z gscalar Z
  if (flav.IsOn()) {
    vertex[vanz].in[0] = flav;
    vertex[vanz].in[1] = flgs;
    vertex[vanz].in[2] = flav;
    Kabbala ma;
    
    if(ScalarNumber(std::string("WidthScheme"))==0){
        ma = Kabbala(string("\\sqr(M_{")+flav.TexName()+string("})"),sqr(flav.Mass()));
      }else{
        ma = Kabbala(string("\\sqr(M_{")+flav.TexName()+string("})"),
                        sqr(flav.Mass())-Complex(0.,1.)*flav.Width()*flav.Mass());
      }
    kcpl0 = M_I*om*kap;
    kcpl1 = ma;

    vertex[vanz].cpl[0]  = kcpl0;
    vertex[vanz].cpl[1]  = kcpl1;
    vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();

    
    vertex[vanz].Color.push_back(Color_Function(cf::None));     

    
    vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("VVGS",LF_Key()));     
    vertex[vanz].Lorentz.back()->SetParticleArg(0,2,1);     

    vertex[vanz].on      = 1;
    vertex.push_back(Single_Vertex());vanz++;
  }
  flav = Flavour(kf_photon);
  // photon gscalar photon
  if (flav.IsOn()) {
    vertex[vanz].in[0] = flav;
    vertex[vanz].in[1] = flgs;
    vertex[vanz].in[2] = flav;
    
    kcpl0 = M_I*om*kap;
    kcpl1 = Kabbala();

    vertex[vanz].cpl[0]  = kcpl0;
    vertex[vanz].cpl[1]  = kcpl1;
    vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();

    
    vertex[vanz].Color.push_back(Color_Function(cf::None));     

    
    vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("VVGS",LF_Key()));     
    vertex[vanz].Lorentz.back()->SetParticleArg(0,2,1);     

    vertex[vanz].on      = 1;
    vertex.push_back(Single_Vertex());vanz++;
  }
}

void Interaction_Model_EW_Grav::c_VVVT(std::vector<Single_Vertex>& vertex,int& vanz)
{
  Flavour flgraviton(kf_graviton);
  Flavour flav(kf_Wplus);
  Kabbala kcpl0,kcpl1,kcpl0_1,kcpl1_1,charge;
  charge = Kabbala(string("Q_{")+flav.TexName()+string("}"),flav.Charge());

  if (!flav.IsOn()) return;
  if (flgraviton.IsOn()){
    // photon WW graviton
    vertex[vanz].nleg     = 4;
    vertex[vanz].in[0] = flav;
    vertex[vanz].in[1] = Flavour(kf_photon);
    vertex[vanz].in[2] = flav;
    vertex[vanz].in[3] = flgraviton;
   
    kcpl0 = -M_I*g1*charge*kap;
    kcpl1 = kcpl0;

    vertex[vanz].cpl[0]  = kcpl0;
    vertex[vanz].cpl[1]  = vertex[vanz].cpl[0];
    vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
    
    
    vertex[vanz].Color.push_back(Color_Function(cf::None));     
    
    
    vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("VVVT",LF_Key()));     
    vertex[vanz].Lorentz.back()->SetParticleArg(0,1,2,3);     
    
    vertex[vanz].on      = 1;
    vertex.push_back(Single_Vertex());vanz++;
    
    // ZWW graviton
    vertex[vanz].nleg     = 4;
    vertex[vanz].in[0] = flav;
    vertex[vanz].in[1] = Flavour(kf_Z);
    vertex[vanz].in[2] = flav;
    vertex[vanz].in[3] = flgraviton;
    
    kcpl0 = -M_I*g2*charge*costW*kap;
    kcpl1 = kcpl0;
    
    vertex[vanz].cpl[0]  = kcpl0;
    vertex[vanz].cpl[1]  = kcpl1;
    vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
    
    
    vertex[vanz].Color.push_back(Color_Function(cf::None));     

    
    vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("VVVT",LF_Key()));     
    vertex[vanz].Lorentz.back()->SetParticleArg(0,1,2,3);     
    
    vertex[vanz].on      = 1;
    vertex.push_back(Single_Vertex());vanz++;
  }
}


void Interaction_Model_EW_Grav::c_SST(std::vector<Single_Vertex>& vertex,int& vanz)
{
  Kabbala kcpl0,kcpl1;

  Flavour flgraviton(kf_graviton);
  Flavour flgs(kf_gscalar);
  Flavour flh = Flavour(kf_h0);

 if (flh.IsOn()&&flgraviton.IsOn()) {  
    vertex[vanz].in[0] = flh;
    vertex[vanz].in[1] = flgraviton;
    vertex[vanz].in[2] = flh;

    kcpl0 = -M_I*kap;
    if(ScalarNumber(std::string("WidthScheme"))==0){
      kcpl1 = Kabbala(string("\\sqr(M_{")+flh.TexName()+string("})"),sqr(flh.Yuk()));
    }else{
      kcpl1 = Kabbala(string("\\sqr(M_{")+flh.TexName()+string("})"),
                      sqr(flh.Yuk())-Complex(0.,1.)*flh.Width()*flh.Yuk());
    }
    
    vertex[vanz].cpl[0]  = kcpl0;
    vertex[vanz].cpl[1]  = kcpl1;
    vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();

    
    vertex[vanz].Color.push_back(Color_Function(cf::None));     
    
    
    vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("SST",LF_Key()));     
    vertex[vanz].Lorentz.back()->SetParticleArg(0,2,1);     	

    vertex[vanz].on      = 1;
    vertex.push_back(Single_Vertex());vanz++;
  }

 if (flh.IsOn()&&flgs.IsOn()) {  
    vertex[vanz].in[0] = flh;
    vertex[vanz].in[1] = flgs;
    vertex[vanz].in[2] = flh;

    kcpl0 = -M_I*om*kap;
    if(ScalarNumber(std::string("WidthScheme"))==0){
      kcpl1 = num2*Kabbala(string("\\sqr(M_{")+flh.TexName()+string("})"),sqr(flh.Yuk()));
    }else{
      kcpl1 = num2*Kabbala(string("\\sqr(M_{")+flh.TexName()+string("})"),
                      sqr(flh.Yuk())-Complex(0.,1.)*flh.Width()*flh.Yuk());
    }
    
    vertex[vanz].cpl[0]  = kcpl0;
    vertex[vanz].cpl[1]  = kcpl1;
    vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();

    
    vertex[vanz].Color.push_back(Color_Function(cf::None));     
    
    
    vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("SSGS",LF_Key()));     
    vertex[vanz].Lorentz.back()->SetParticleArg(0,2);     	

    vertex[vanz].on      = 1;
    vertex.push_back(Single_Vertex());vanz++;
  }
}


void Interaction_Model_EW_Grav::c_SSST(std::vector<Single_Vertex>& vertex,int& vanz)
{
  Kabbala kcpl0,kcpl1,yuk;
  Kabbala num3  = Kabbala(string("3"),3.);

  Flavour flgraviton(kf_graviton);
  Flavour flgs(kf_gscalar);
  Flavour flh = Flavour(kf_h0);

 if (flh.IsOn()&&flgraviton.IsOn()) {  
    vertex[vanz].nleg     = 4;
    vertex[vanz].in[0] = flh;
    vertex[vanz].in[1] = flh;
    vertex[vanz].in[2] = flh;
    vertex[vanz].in[3] = flgraviton;


    if(ScalarNumber(std::string("WidthScheme"))==0){
      yuk = Kabbala(string("M_{")+flh.TexName()+string("}"),flh.Yuk());
    }else{
      yuk = Kabbala(string("M_{")+flh.TexName()+string("}"),
                      sqrt(sqr(flh.Yuk())-Complex(0.,1.)*flh.Width()*flh.Yuk()));
    }
    kcpl0 = -M_I*kap*yuk*yuk*(num3/vev);
    kcpl1 = kcpl0;
    //kcpl0=Kabbala("",320);
    //kcpl1=Kabbala("",321);
    
    vertex[vanz].cpl[0]  = kcpl0;
    vertex[vanz].cpl[1]  = kcpl1;
    vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();

    
    vertex[vanz].Color.push_back(Color_Function(cf::None));     
    
    
    vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("SSST",LF_Key()));     
    vertex[vanz].Lorentz.back()->SetParticleArg(3);     	

    vertex[vanz].on      = 1;
    vertex.push_back(Single_Vertex());vanz++;
  }

 if (flh.IsOn()&&flgs.IsOn()) {  
    vertex[vanz].nleg     = 4;
    vertex[vanz].in[0] = flh;
    vertex[vanz].in[1] = flh;
    vertex[vanz].in[2] = flh;
    vertex[vanz].in[3] = flgs;

    if(ScalarNumber(std::string("WidthScheme"))==0){
      yuk = Kabbala(string("M_{")+flh.TexName()+string("}"),flh.Yuk());
    }else{
      yuk = Kabbala(string("M_{")+flh.TexName()+string("}"),
                      sqrt(sqr(flh.Yuk())-Complex(0.,1.)*flh.Width()*flh.Yuk()));
    }
    kcpl0 = -M_I*num2*om*kap*yuk*yuk*(num3/vev);
    kcpl1 = kcpl0;
    //kcpl0=Kabbala("",300);
    //kcpl1=Kabbala("",0);
    
    vertex[vanz].cpl[0]  = kcpl0;
    vertex[vanz].cpl[1]  = kcpl1;
    vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();

    
    vertex[vanz].Color.push_back(Color_Function(cf::None));     
    
    
    vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("SSSS",LF_Key()));     

    vertex[vanz].on      = 1;
    vertex.push_back(Single_Vertex());vanz++;
  }
}

Kabbala Interaction_Model_EW_Grav::K_CKM(short int i,short int j)       
{   
  char hi[2],hj[2];
  sprintf(hi,"%i",i);
  sprintf(hj,"%i",j);
  return Kabbala(string("V_{")+string(hi)+string(hj)+string("}"),
		 ComplexMatrixElement(std::string("CKM"),i,j));
} 
  
Kabbala Interaction_Model_EW_Grav::conj_K_CKM(short int i,short int j)       
{   
  char hi[2],hj[2];
  sprintf(hi,"%i",i);
  sprintf(hj,"%i",j);
  return Kabbala(string("V^\\ti_{")+string(hi)+string(hj)+string("}"),
		 conj(ComplexMatrixElement(std::string("CKM"),i,j)));
} 
 
