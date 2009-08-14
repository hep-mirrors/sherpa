#include "MODEL/Interaction_Models/Interaction_Model_AEW.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include <stdio.h>


using namespace MODEL;
using namespace ATOOLS;
using namespace std;

Interaction_Model_AEW::Interaction_Model_AEW(MODEL::Model_Base * _model,
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
  vev   = Kabbala(string("v_{EW}"),ScalarConstant(std::string("vev")));
  a4    = Kabbala(string("\\alpha_4"),ScalarConstant(std::string("Alpha_4")));
  a5    = Kabbala(string("\\alpha_5"),ScalarConstant(std::string("Alpha_5")));

  g1_p     = Kabbala(string("g_1^{gamma}"),ScalarConstant(std::string("g1_gamma")));
  kappa_p  = Kabbala(string("\\kappa^{gamma}"),ScalarConstant(std::string("kappa_gamma")));
  lambda_p = Kabbala(string("\\lambda^{gamma}"),ScalarConstant(std::string("lambda_gamma")));
  g4_p     = Kabbala(string("g_4^{gamma}"),ScalarConstant(std::string("g4_gamma")));
  g5_p     = Kabbala(string("g_5^{gamma}"),ScalarConstant(std::string("g5_gamma")));
  kappat_p = Kabbala(string("\\kappat^{gamma}"),ScalarConstant(std::string("kappat_gamma")));
  lambdat_p= Kabbala(string("\\lambdat^{gamma}"),ScalarConstant(std::string("lambdat_gamma")));
  g1_Z     = Kabbala(string("g_1^{Z}"),ScalarConstant(std::string("g1_Z")));
  kappa_Z  = Kabbala(string("\\kappa^{Z}"),ScalarConstant(std::string("kappa_Z")));
  lambda_Z = Kabbala(string("\\lambda^{Z}"),ScalarConstant(std::string("lambda_Z")));
  g4_Z     = Kabbala(string("g_4^{Z}"),ScalarConstant(std::string("g4_Z")));
  g5_Z     = Kabbala(string("g_5^{Z}"),ScalarConstant(std::string("g5_Z")));
  kappat_Z = Kabbala(string("\\kappat^{Z}"),ScalarConstant(std::string("kappat_Z")));
  lambdat_Z= Kabbala(string("\\lambdat^{Z}"),ScalarConstant(std::string("lambdat_Z")));

  f4_Z     = Kabbala(string("f_4^{Z}"),ScalarConstant(std::string("f4_Z")));
  f5_Z     = Kabbala(string("f_5^{Z}"),ScalarConstant(std::string("f5_Z")));
  f4_p     = Kabbala(string("f_4^{gamma}"),ScalarConstant(std::string("f4_gamma")));
  f5_p     = Kabbala(string("f_5^{gamma}"),ScalarConstant(std::string("f5_gamma")));
  
  h1_Z     = Kabbala(string("h_1^{Z}"),ScalarConstant(std::string("h1_Z")));
  h2_Z     = Kabbala(string("h_2^{Z}"),ScalarConstant(std::string("h2_Z")));
  h3_Z     = Kabbala(string("h_3^{Z}"),ScalarConstant(std::string("h3_Z")));
  h4_Z     = Kabbala(string("h_4^{Z}"),ScalarConstant(std::string("h4_Z")));
  h1_p     = Kabbala(string("h_1^{gamma}"),ScalarConstant(std::string("h1_gamma")));
  h2_p     = Kabbala(string("h_2^{gamma}"),ScalarConstant(std::string("h2_gamma")));
  h3_p     = Kabbala(string("h_3^{gamma}"),ScalarConstant(std::string("h3_gamma")));
  h4_p     = Kabbala(string("h_4^{gamma}"),ScalarConstant(std::string("h4_gamma")));
  
}

void Interaction_Model_AEW::c_FFV(std::vector<Single_Vertex>& vertex,int & vanz)
{
  Flavour flphoton(kf_photon);
  Flavour flZ(kf_Z);
  Flavour flWplus(kf_Wplus);
  for (short int i=1;i<17;i++) {
    if (i==7) i=11;
    Flavour flav1               = Flavour((kf_code)(i));
    Kabbala charge1             = Kabbala(string("Q_{")+flav1.TexName()+string("}"),flav1.Charge());
    Kabbala isoweak1            = Kabbala(string("T_{")+flav1.TexName()+string("}"),flav1.IsoWeak());

    Kabbala kcpl0,kcpl1;    
    if (flav1.IsOn()) {
      for (short int j=i;j<17;j++) {
	if (j==7) j=11;
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
	      if (i<7 && j<7) {
		if (flav1.IsDowntype())
		  kcpl1 = -M_I/root2*g2*K_CKM(j/2-1,(i-1)/2);
		else
		  kcpl1 = -M_I/root2*g2*K_CKM(i/2-1,(j-1)/2);
	      }
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

// Anomalous triple gauge couplings as defined in Nuclear Physics B282 (1987)253-307
// only terms ~f1,f2,f3,f4

void Interaction_Model_AEW::c_VVV(std::vector<Single_Vertex>& vertex,int& vanz)
{
  Flavour flWplus(kf_Wplus);
  Flavour flZ(kf_Z);
  Flavour flP(kf_photon);
  Kabbala kcpl,kcpl0,kcpl1,kcpl2,kcpl3,charge;
  Kabbala Wyuk = Kabbala(string("M_{")+flWplus.TexName()+string("}"),flWplus.Yuk());
  Kabbala Zyuk = Kabbala(string("M_{")+flZ.TexName()+string("}"),flZ.Yuk());

  charge = Kabbala(string("Q_{")+flWplus.TexName()+string("}"),flWplus.Charge());

  if (flWplus.IsOn() && flP.IsOn()) {

    // photon WW
    vertex[vanz].in[0] = flWplus;
    vertex[vanz].in[1] = flP;
    vertex[vanz].in[2] = flWplus;
    
    kcpl = M_I*g1*charge;
    kcpl0 = kcpl*g1_p;
    kcpl1 = kcpl*kappa_p;
    kcpl2 = kcpl*lambda_p/(Wyuk*Wyuk);
    kcpl3 = kcpl*M_I*g4_p;
        
    vertex[vanz].cpl[0]  = kcpl0;  //g1
    vertex[vanz].cpl[1]  = kcpl1;  //kappa
    vertex[vanz].cpl[2]  = kcpl2;  //lambda/m_W^2
    vertex[vanz].cpl[3]  = kcpl3;  //i*g4
    vertex[vanz].Str=kcpl0.String()+"|"+kcpl1.String()+"|"+kcpl2.String()+"|"+kcpl3.String();

    kcpl1 = kcpl*kappat_p;
    kcpl2 = kcpl*lambdat_p/(Wyuk*Wyuk);
    kcpl3 = kcpl*M_I*g5_p;
    vertex[vanz].cpl.push_back(kcpl3);  //i*g5
    vertex[vanz].cpl.push_back(kcpl1);  //kappat
    vertex[vanz].cpl.push_back(kcpl2);  //lambdat/m_W^2
    vertex[vanz].cpl.push_back(kcpl);  //overall coupling factor
    vertex[vanz].Str+=kcpl3.String()+"|"+kcpl1.String()+"|"+kcpl2.String();
    
    
    vertex[vanz].Color.push_back(Color_Function(cf::None));     
    
    
    vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("AGauge3",LF_Key()));
    vertex[vanz].Lorentz.back()->SetParticleArg(1,2,0);     
  
    vertex[vanz].on      = 1;
    vertex.push_back(Single_Vertex());vanz++;
  }
  if (flZ.IsOn()) {
 
    // ZWW
    vertex[vanz].in[0] = flWplus;
    vertex[vanz].in[1] = flZ;
    vertex[vanz].in[2] = flWplus;

    kcpl = M_I*g2*charge*costW;
    kcpl0 = kcpl*g1_Z;
    kcpl1 = kcpl*kappa_Z;
    kcpl2 = kcpl*lambda_Z/(Wyuk*Wyuk);
    kcpl3 = kcpl*M_I*g4_Z;
  
    vertex[vanz].cpl[0]  = kcpl0;  //g1
    vertex[vanz].cpl[1]  = kcpl1;  //kappa
    vertex[vanz].cpl[2]  = kcpl2;  //lambda/m_W^2
    vertex[vanz].cpl[3]  = kcpl3;  //i*g4

    kcpl1 = kcpl*kappat_Z;
    kcpl2 = kcpl*lambdat_Z/(Wyuk*Wyuk);
    kcpl3 = kcpl*M_I*g5_Z;
    vertex[vanz].cpl.push_back(kcpl3);  //i*g5
    vertex[vanz].cpl.push_back(kcpl1);  //kappat
    vertex[vanz].cpl.push_back(kcpl2);  //lambdat/m_W^2
    vertex[vanz].cpl.push_back(kcpl);  //overall coupling factor

    
    vertex[vanz].Color.push_back(Color_Function(cf::None));     
    
    
    vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("AGauge3",LF_Key()));
    vertex[vanz].Lorentz.back()->SetParticleArg(1,2,0);     

    vertex[vanz].on      = 1;
    vertex.push_back(Single_Vertex());vanz++;
  }
  

  if (flZ.IsOn()) {
    //ZZZ
    vertex[vanz].in[0] = flZ;
    vertex[vanz].in[1] = flZ;
    vertex[vanz].in[2] = flZ;

    kcpl = M_I*g1;
    kcpl0 = kcpl*f4_Z/(Zyuk*Zyuk);
    kcpl1 = kcpl*f5_Z/(Zyuk*Zyuk);
  
    vertex[vanz].cpl[0]  = kcpl0;  //f4
    vertex[vanz].cpl[1]  = kcpl1;  //f5

    vertex[vanz].Color.push_back(Color_Function(cf::None));     
       
    vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("AZZZ",LF_Key()));
    vertex[vanz].Lorentz.back()->SetParticleArg(1,2,0);     

    vertex[vanz].on      = 1;
    vertex.push_back(Single_Vertex());vanz++;
  
    // ZZ Photon
    vertex[vanz].in[0] = flZ;
    vertex[vanz].in[1] = flP;
    vertex[vanz].in[2] = flZ;

    kcpl = M_I*g1;
    kcpl0 = kcpl*f4_p/(Zyuk*Zyuk);
    kcpl1 = kcpl*f5_p/(Zyuk*Zyuk);
  
    vertex[vanz].cpl[0]  = kcpl0;  //f4
    vertex[vanz].cpl[1]  = kcpl1;  //f5

    kcpl0 = kcpl*h1_Z/(Zyuk*Zyuk);
    kcpl1 = kcpl*h2_Z/(Zyuk*Zyuk*Zyuk*Zyuk);
    kcpl2 = kcpl*h3_Z/(Zyuk*Zyuk);
    kcpl3 = kcpl*h4_Z/(Zyuk*Zyuk*Zyuk*Zyuk);

    vertex[vanz].cpl[2]  = kcpl0;       //h1
    vertex[vanz].cpl[3]  = kcpl1;       //h2
    vertex[vanz].cpl.push_back(kcpl2);  //h3
    vertex[vanz].cpl.push_back(kcpl3);  //h4
    vertex[vanz].cpl.push_back(kcpl0);  //h1
    vertex[vanz].cpl.push_back(kcpl1);  //h2

    vertex[vanz].Color.push_back(Color_Function(cf::None));     
       
    vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("AZZG",LF_Key()));
    vertex[vanz].Lorentz.back()->SetParticleArg(1,2,0);     

    vertex[vanz].on      = 1;
    vertex.push_back(Single_Vertex());vanz++;


    // Z Photon Photon
    vertex[vanz].in[0] = flP;
    vertex[vanz].in[1] = flP;
    vertex[vanz].in[2] = flZ;

    kcpl = M_I*g1;

    kcpl0 = kcpl*h1_p/(Zyuk*Zyuk);
    kcpl1 = kcpl*h2_p/(Zyuk*Zyuk*Zyuk*Zyuk);
    kcpl2 = kcpl*h3_p/(Zyuk*Zyuk);
    kcpl3 = kcpl*h4_p/(Zyuk*Zyuk*Zyuk*Zyuk);

    vertex[vanz].cpl[0]  = kcpl0;       //h1
    vertex[vanz].cpl[1]  = kcpl1;       //h2
    vertex[vanz].cpl[2]  = kcpl2;       //h3
    vertex[vanz].cpl[3]  = kcpl3;       //h4

    vertex[vanz].Color.push_back(Color_Function(cf::None));     
       
    vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("AZGG",LF_Key()));
    vertex[vanz].Lorentz.back()->SetParticleArg(1,2,0);     

    vertex[vanz].on      = 1;
    vertex.push_back(Single_Vertex());vanz++;
  }
}
//End anomalous triple gauge couplings

//Anomalous quadrupole Vertices defined in hep-ph/0001065

void Interaction_Model_AEW::c_VVVV(std::vector<Single_Vertex>& vertex,int& vanz)
{
  Flavour flavWplus(kf_Wplus);
  Flavour flavZ(kf_Z);
  Flavour flavP(kf_photon);
  Kabbala kcpl0,kcpl1;
  Kabbala num_2 = Kabbala(string("2"),2.);  
  Kabbala num_3 = Kabbala(string("3"),3.);  
  
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
  
    
    vertex[vanz].Color.push_back(Color_Function(cf::None));     
    
    
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
  
    
    vertex[vanz].Color.push_back(Color_Function(cf::None));     
    
    
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

    kcpl0 = -M_I*(g1*g1*costW*costW/(sintW*sintW)-a5/(costW*costW));
    kcpl1 = -M_I*(g1*g1*costW*costW/(sintW*sintW)+a4/(costW*costW));
  
    vertex[vanz].cpl[0]  = kcpl0;
    vertex[vanz].cpl[1]  = kcpl1;
    vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
  
    
    vertex[vanz].Color.push_back(Color_Function(cf::None));     
    
    
    vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("AGauge4",LF_Key()));     
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

    kcpl0 = M_I*(g1*g1/(sintW*sintW)+a4);
    kcpl1 = M_I*(g1*g1/(sintW*sintW)-a4-num_2*a5);
  
    vertex[vanz].cpl[0]  = kcpl0;
    vertex[vanz].cpl[1]  = kcpl1;
    vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
  
    
    vertex[vanz].Color.push_back(Color_Function(cf::None));     
    
    
    vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("AGauge4",LF_Key())); 
    vertex[vanz].Lorentz.back()->SetParticleArg(0,1,2,3);     

    vertex[vanz].on      = 1;
    vertex.push_back(Single_Vertex());vanz++;
  }

  // Z - Z - Z - Z  
  if (flavZ.IsOn()) {
    vertex[vanz].in[0] = flavZ;
    vertex[vanz].in[1] = flavZ;
    vertex[vanz].in[2] = flavZ;
    vertex[vanz].in[3] = flavZ;
  
    vertex[vanz].nleg     = 4;

    kcpl0 = M_I*num_3*(a4+a5)/(costW*costW*costW*costW);
    kcpl1 = Kabbala(string("0"),Complex(0.,0.));
  
    vertex[vanz].cpl[0]  = kcpl0;
    vertex[vanz].cpl[1]  = kcpl1;
    vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
  
    
    vertex[vanz].Color.push_back(Color_Function(cf::None));     
    
    
    vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("AGauge4",LF_Key())); 
    vertex[vanz].Lorentz.back()->SetParticleArg(0,1,2,3);     

    vertex[vanz].on      = 1;
    vertex.push_back(Single_Vertex());vanz++;
  }
}
//End anoumalous quadrpol vertices


Kabbala Interaction_Model_AEW::K_CKM(short int i,short int j)       
{   
  char hi[2],hj[2];
  sprintf(hi,"%i",i);
  sprintf(hj,"%i",j);
  return Kabbala(string("V_{")+string(hi)+string(hj)+string("}"),
		 ComplexMatrixElement(std::string("CKM"),i,j));
} 
  
Kabbala Interaction_Model_AEW::conj_K_CKM(short int i,short int j)       
{   
  char hi[2],hj[2];
  sprintf(hi,"%i",i);
  sprintf(hj,"%i",j);
  return Kabbala(string("V^\\ti_{")+string(hi)+string(hj)+string("}"),
		 conj(ComplexMatrixElement(std::string("CKM"),i,j)));
} 
 


