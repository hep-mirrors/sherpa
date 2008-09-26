#include "Interaction_Model_EW.H"
#include "MathTools.H"
#include "Message.H"
#include "Run_Parameter.H"
#include <stdio.h>


using namespace MODEL;
using namespace ATOOLS;
using namespace std;

DECLARE_GETTER(Interaction_Model_EW_Getter,"pure_EW",
	       Interaction_Model_Base,Interaction_Model_Arguments);

Interaction_Model_Base *Interaction_Model_EW_Getter::
operator()(const Interaction_Model_Arguments &args) const
{
  return new Interaction_Model_EW
    (args.p_model,args.m_cplscheme,args.m_yukscheme);
}

void Interaction_Model_EW_Getter::PrintInfo
(std::ostream &str,const size_t width) const
{ 
  str<<"The SM electroweak sector only"; 
}

Interaction_Model_EW::Interaction_Model_EW(MODEL::Model_Base * _model,
					   std::string _cplscheme,std::string _yukscheme) :
  Interaction_Model_Base("pure_EW",_model,_cplscheme,_yukscheme)
{ 
  g1       = Kabbala(string("g_1"),
		     sqrt(4.*M_PI*ScalarFunction(std::string("alpha_QED"),rpa.gen.CplScale())));
  g2       = Kabbala(string("g_1/\\sin\\theta_W"), 
		     g1.Value()/sqrt(ScalarConstant(std::string("sin2_thetaW"))));
  sintW    = Kabbala(std::string("\\sin\\theta_W"),
		     sqrt(ScalarConstant(std::string("sin2_thetaW"))));
  costW    = Kabbala(std::string("\\cos\\theta_W"),
		     sqrt(1.-ScalarConstant(std::string("sin2_thetaW"))));

  PL       = Kabbala(string("P_L"),1.);
  PR       = Kabbala(string("P_R"),1.);
  M_I      = Kabbala(string("i"),Complex(0.,1.));
  root2    = Kabbala(string("\\sqrt{2}"),sqrt(2.));
  vev      = Kabbala(string("v_{EW}"),ScalarConstant(std::string("vev")));

  if (_model->GetScalarNumbers()->find("Extension")!=_model->GetScalarNumbers()->end()) {
    m_extension = (*_model->GetScalarNumbers())["Extension"];
  }
  if (m_extension==1) {
    tbW_f    = Kabbala(std::string("\\kappa_{tbW}"),
		       ScalarConstant(std::string("tbW_RelFactor")));
    tbW_cosa = Kabbala(std::string("\\cos\\theta_{tbW}"),
		       cos(ScalarConstant(std::string("tbW_Angle"))));
    tbW_sina = Kabbala(std::string("\\sin\\theta_{tbW}"),
		       sin(ScalarConstant(std::string("tbW_Angle"))));
  }
}

void Interaction_Model_EW::c_FFV(std::vector<Single_Vertex>& vertex,int & vanz)
{
  Flavour flphoton(kf_photon);
  Flavour flZ(kf_Z);
  Flavour flWplus(kf_Wplus);
  for (short int i=1;i<(m_extension==2?19:17);i++) {
    if (i==(m_extension==2?9:7)) i=11;
    Flavour flav1               = Flavour((kf_code)(i));
    Kabbala charge1             = Kabbala(string("Q_{")+flav1.TexName()+string("}"),flav1.Charge());
    Kabbala isoweak1            = Kabbala(string("T_{")+flav1.TexName()+string("}"),flav1.IsoWeak());

    Kabbala kcpl0,kcpl1;    
    if (flav1.IsOn()) {
      for (short int j=i;j<(m_extension==2?19:17);j++) {
	if (j==(m_extension==2?9:7)) j=11;
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
		  vertex[vanz].Color.push_back(Color_Function(cf::D));;     
		  vertex[vanz].Color.back().SetParticleArg(0,2);     
		  vertex[vanz].Color.back().SetStringArg('0','2');     
		}
		else 
		  vertex[vanz].Color.push_back(Color_Function(cf::None));;

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
		vertex[vanz].Color.push_back(Color_Function(cf::D));;     
		vertex[vanz].Color.back().SetParticleArg(0,2);     
		vertex[vanz].Color.back().SetStringArg('0','2');     
	      }
	      else 
		vertex[vanz].Color.push_back(Color_Function(cf::None));;

	      
	      vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("Gamma",LF_Key()));
	      vertex[vanz].Lorentz.back()->SetParticleArg(1);     

	      vertex[vanz].on     = 1;
	      vertex.push_back(Single_Vertex());vanz++;
	    }
	  }
	  //W
	  if (flWplus.IsOn()) {
	    bool hit = 1;
	    Kabbala kcpl0,kcpl1;
	    kcpl0 = Kabbala(string("zero"),0.);
	    kcpl1 = Kabbala(string("1"),0.);

	    if (!((flav1.IsDowntype() && flav2.IsUptype()) ||
                  (flav2.IsDowntype() && flav1.IsUptype()))) hit = false;
	    else if ((flav1.IsLepton() && !flav2.IsLepton()) ||
		     (flav1.IsQuark() && !flav2.IsQuark()) ) hit = false;
	    else if (flav1.IsDowntype() && i>10 && j>10) {
	      if (flav1.IsDowntype())
		kcpl1 = -M_I/root2*g2*K_L_CKM((j-10)/2-1,(i-11)/2);
	      else
		kcpl1 = -M_I/root2*g2*K_L_CKM((i-10)/2-1,(j-11)/2);
	    }
	    else if (i<10 && j<10) {
	      if (flav1.IsDowntype())
		kcpl1 = -M_I/root2*g2*K_CKM(j/2-1,(i-1)/2);
	      else
		kcpl1 = -M_I/root2*g2*K_CKM(i/2-1,(j-1)/2);
	      if (m_extension==1 && (tbW_f.Value()!=1. || tbW_cosa.Value()!=1.)) {
		if ((flav1==Flavour(kf_t) && flav2==Flavour(kf_b)) ||
		    (flav1==Flavour(kf_b) && flav2==Flavour(kf_t))) {
		  kcpl0 = kcpl1*tbW_sina*tbW_f;
		  kcpl1 = kcpl1*tbW_cosa*tbW_f;
                }
	      }
	    }
	    if (hit) {
	      if (!ATOOLS::IsZero(kcpl0.Value()) ||
		  !ATOOLS::IsZero(kcpl1.Value())) {
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
		  vertex[vanz].Color.push_back(Color_Function(cf::D));;     
		  vertex[vanz].Color.back().SetParticleArg(0,2);     
		  vertex[vanz].Color.back().SetStringArg('0','2');     
		}
		else 
		  vertex[vanz].Color.push_back(Color_Function(cf::None));;

		
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

void Interaction_Model_EW::c_VVV(std::vector<Single_Vertex>& vertex,int& vanz)
{
  Flavour flW(kf_Wplus);
  Flavour flZ(kf_Z);
  Flavour flP(kf_photon);
  Kabbala kcpl0,kcpl1,kcpl0_1,kcpl1_1,charge;

  charge = Kabbala(string("Q_{")+flW.TexName()+string("}"),flW.Charge());

  if (flW.IsOn() && flP.IsOn()) {

    // photon WW
    vertex[vanz].in[0] = flW;
    vertex[vanz].in[1] = Flavour(kf_photon);
    vertex[vanz].in[2] = flW;
    
    kcpl0 = M_I*g1*charge;
    kcpl1 = kcpl0;
    
    vertex[vanz].cpl[0]  = kcpl0;
    vertex[vanz].cpl[1]  = vertex[vanz].cpl[0];
    vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
    
    
    vertex[vanz].Color.push_back(Color_Function(cf::None));;     
    
    
    vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("Gauge3",LF_Key()));
    vertex[vanz].Lorentz.back()->SetParticleArg(0,1,2);     
    
    vertex[vanz].on      = 1;
    vertex.push_back(Single_Vertex());vanz++;
  }
  if (flZ.IsOn()) {
    
    // ZWW
    vertex[vanz].in[0] = flW;
    vertex[vanz].in[1] = Flavour(kf_Z);
    vertex[vanz].in[2] = flW;
    
    kcpl0 = M_I*g2*charge*costW;
    kcpl1 = kcpl0;
    
    vertex[vanz].cpl[0]  = kcpl0;
    vertex[vanz].cpl[1]  = kcpl1;
    vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
    
    
    vertex[vanz].Color.push_back(Color_Function(cf::None));;     
    
    
    vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("Gauge3",LF_Key()));
    vertex[vanz].Lorentz.back()->SetParticleArg(0,1,2);     
    
    vertex[vanz].on      = 1;
    vertex.push_back(Single_Vertex());vanz++;
  }
}

void Interaction_Model_EW::c_VVVV(std::vector<Single_Vertex>& vertex,int& vanz)
{
  Flavour flavWplus(kf_Wplus);
  Flavour flavZ(kf_Z);
  Flavour flavP(kf_photon);
  Kabbala kcpl0,kcpl1;
  
  // Ph - W - W - Ph  
  if (flavWplus.IsOn() && flavP.IsOn()) {
    vertex[vanz].in[0] = flavP;
    vertex[vanz].in[1] = flavWplus.Bar();
    vertex[vanz].in[2] = flavWplus;
    vertex[vanz].in[3] = flavP;
  
    vertex[vanz].nleg     = 4;

    kcpl0 = -M_I*g1*g1;
    kcpl1 = kcpl0;
    
    vertex[vanz].cpl[0]  = kcpl0;
    vertex[vanz].cpl[1]  = kcpl1;
    vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
  
    
    vertex[vanz].Color.push_back(Color_Function(cf::None));;     
    
    
    vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("Gauge4",LF_Key()));
    vertex[vanz].Lorentz.back()->SetParticleArg(0,3,1,2);     

    vertex[vanz].on      = 1;
    vertex.push_back(Single_Vertex());vanz++;
  }

  // Ph - W - W - Z  
  if (flavWplus.IsOn() && flavP.IsOn() && flavZ.IsOn()) {
    vertex[vanz].in[0] = flavP;
    vertex[vanz].in[1] = flavWplus.Bar();
    vertex[vanz].in[2] = flavWplus;
    vertex[vanz].in[3] = flavZ;

    vertex[vanz].nleg     = 4;  

    kcpl0 = -M_I*g1*g1*costW/sintW;
    kcpl1 = kcpl0;
  
    vertex[vanz].cpl[0]  = kcpl0;
    vertex[vanz].cpl[1]  = kcpl1;
    vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
  
    
    vertex[vanz].Color.push_back(Color_Function(cf::None));;     
    
    
    vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("Gauge4",LF_Key()));     
    vertex[vanz].Lorentz.back()->SetParticleArg(0,3,1,2);     

    vertex[vanz].on      = 1;
    vertex.push_back(Single_Vertex());vanz++;
  }
  // Z - W - W - Z  
  if (flavWplus.IsOn() && flavZ.IsOn()) {
    vertex[vanz].in[0] = flavZ;
    vertex[vanz].in[1] = flavWplus.Bar();
    vertex[vanz].in[2] = flavWplus;
    vertex[vanz].in[3] = flavZ;
  
    vertex[vanz].nleg     = 4;

    kcpl0 = -M_I*g1*g1*costW*costW/(sintW*sintW);
    kcpl1 = kcpl0;
  
    vertex[vanz].cpl[0]  = kcpl0;
    vertex[vanz].cpl[1]  = kcpl1;
    vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
  
    
    vertex[vanz].Color.push_back(Color_Function(cf::None));;     
    
    
    vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("Gauge4",LF_Key()));     
    vertex[vanz].Lorentz.back()->SetParticleArg(0,3,1,2);     

    vertex[vanz].on      = 1;
    vertex.push_back(Single_Vertex());vanz++;
  }
  
  // W - W - W - W  
  if (flavWplus.IsOn()) {
    vertex[vanz].in[0] = flavWplus;
    vertex[vanz].in[1] = flavWplus.Bar();
    vertex[vanz].in[2] = flavWplus;
    vertex[vanz].in[3] = flavWplus;
  
    vertex[vanz].nleg     = 4;

    kcpl0 = M_I*g1*g1/(sintW*sintW);
    kcpl1 = kcpl0;
  
    vertex[vanz].cpl[0]  = kcpl0;
    vertex[vanz].cpl[1]  = kcpl1;
    vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
  
    
    vertex[vanz].Color.push_back(Color_Function(cf::None));;     
    
    
    vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("Gauge4",LF_Key())); 
    vertex[vanz].Lorentz.back()->SetParticleArg(0,1,2,3);     

    vertex[vanz].on      = 1;
    vertex.push_back(Single_Vertex());vanz++;
  }
}

Kabbala Interaction_Model_EW::K_CKM(short int i,short int j)       
{   
  char hi[2],hj[2];
  sprintf(hi,"%i",i);
  sprintf(hj,"%i",j);
  return Kabbala(string("V_{")+string(hi)+string(hj)+string("}"),
		 ComplexMatrixElement(std::string("CKM"),i,j));
} 
  
Kabbala Interaction_Model_EW::K_L_CKM(short int i,short int j)       
{   
  char hi[2],hj[2];
  sprintf(hi,"%i",i);
  sprintf(hj,"%i",j);
  return Kabbala(string("V^L_{")+string(hi)+string(hj)+string("}"),
		 ComplexMatrixElement(std::string("L_CKM"),i,j));
} 
  
Kabbala Interaction_Model_EW::conj_K_CKM(short int i,short int j)       
{   
  char hi[2],hj[2];
  sprintf(hi,"%i",i);
  sprintf(hj,"%i",j);
  return Kabbala(string("V^\\ti_{")+string(hi)+string(hj)+string("}"),
		 conj(ComplexMatrixElement(std::string("CKM"),i,j)));
} 
 


