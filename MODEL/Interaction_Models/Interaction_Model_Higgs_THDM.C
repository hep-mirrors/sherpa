#include "MODEL/Interaction_Models/Interaction_Model_Higgs_THDM.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include <stdio.h>


using namespace MODEL;
using namespace ATOOLS;
using namespace std;

Interaction_Model_Higgs_THDM::Interaction_Model_Higgs_THDM(MODEL::Model_Base * _model,
						 std::string _cplscheme,std::string _yukscheme) :
  Interaction_Model_Base("",_model,_cplscheme,_yukscheme)
{ 
  g1     = Kabbala(string("g_1"),
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

  PL     = Kabbala(string("P_L"),1.);
  PR     = Kabbala(string("P_R"),1.);
  M_I    = Kabbala(string("i"),Complex(0.,1.));
  v1     = Kabbala(string("v_1"),
		   vev.Value() *
		   sqrt(1./(1.+sqr(ScalarConstant(std::string("tan(beta)"))))));
  v2     = Kabbala(string("v_2"),v1.Value()*ScalarConstant(std::string("tan(beta)")));
  root2  = Kabbala(string("\\sqrt{2}"),sqrt(2.));
  K_zero = Kabbala(string("zero"),0.);
  num_2  = Kabbala(string("2"),2.);    	
  num_3  = Kabbala(string("3"),3.);    	
  num_4  = Kabbala(string("4"),4.);    		
}


void Interaction_Model_Higgs_THDM::c_FFS(std::vector<Single_Vertex>& vertex,int& vanz)
{
  Flavour flHplus(kf_Hplus);
  Flavour flA0(kf_A0);
  
  // l(e)/q(d) -> h0/H0 + l(e)/q(d)
  
  for(short int i=25;i<36;i+=10) {
    
    Flavour flav = Flavour((kf_code)(i));
    Kabbala kcpl0,kcpl1;

    if(flav.IsOn()) {
      
      for(short int j=1;j<17;j++) {
	if (j==7) j=11;
	Flavour fl1 = Flavour((kf_code)(j)); 
       	if(fl1.IsOn() && fl1.IsFermion() && fl1.IsDowntype() && (fl1.Yuk() > 1.)) {
	  
	  vertex[vanz].in[0] = fl1; 
	  vertex[vanz].in[1] = flav;
	  vertex[vanz].in[2] = fl1; 
	  
	  kcpl0 = -M_I/v1*K_yuk(fl1,flav)*K_Z_R(0,(i-25)/10);
	  kcpl1 = kcpl0;

	  vertex[vanz].cpl[0]  = kcpl0;
	  vertex[vanz].cpl[1]  = kcpl1;
	  vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
	  vertex[vanz].on      = 1;

	  
 	  if (fl1.Strong()) {  
	    vertex[vanz].Color.push_back(Color_Function(cf::D));     
	    vertex[vanz].Color.back().SetParticleArg(0,2);     
	    vertex[vanz].Color.back().SetStringArg('0','2');     
	  }
	  else 
	    vertex[vanz].Color.push_back(Color_Function(cf::None));
	
	  
	  vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("FFS",LF_Key()));
	  	  
	  vertex.push_back(Single_Vertex());vanz++;
	  //checked FK & RK	  
	}
      }
      
      // q(u) -> h0/H0 + q(u)
      
      for(short int t=1;t<17;t++) {
	if (t==7) t=11;
	Flavour fl1 = Flavour((kf_code)(t)); 
	if(fl1.IsOn() && fl1.IsQuark() && fl1.IsUptype() && (fl1.Yuk() > 1.)) {
	  
	  vertex[vanz].in[0] = fl1; 
	  vertex[vanz].in[1] = flav;
	  vertex[vanz].in[2] = fl1; 
	  
	  kcpl0 = -M_I/v2*K_yuk(fl1,flav)*K_Z_R(1,(i-25)/10);
	  kcpl1 = kcpl0;
	  
	  vertex[vanz].cpl[0]  = kcpl0;
	  vertex[vanz].cpl[1]  = kcpl1;
	  vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
	  vertex[vanz].on      = 1;
	  
	  
	  vertex[vanz].Color.push_back(Color_Function(cf::D));     
	  vertex[vanz].Color.back().SetParticleArg(0,2);     
	  vertex[vanz].Color.back().SetStringArg('0','2');     
	  
	  
	  vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("FFS",LF_Key()));
	  vertex[vanz].Lorentz.back()->SetParticleArg(1);     
	  
	  vertex.push_back(Single_Vertex());vanz++;
	  //checked FK & RK
	}
      }
    }
  }
  
  // l(e)/q(d) -> A0 + l(e)/q(d)
  
  if(flA0.IsOn()) {

    Kabbala kcpl0,kcpl1;
 
    for(short int k=1;k<17;k++) {
      if (k==7) k=11;
      Flavour fl1 = Flavour((kf_code)(k)); 
      if(fl1.IsOn() && fl1.IsFermion() && fl1.IsDowntype() && (fl1.Yuk() > 1.)) {
		
	vertex[vanz].in[0] = fl1; 
	vertex[vanz].in[1] = flA0;
	vertex[vanz].in[2] = fl1; 


	kcpl0 = -K_yuk(fl1,flA0)/v1*K_Z_H(0,0);
	kcpl1 = -kcpl0;

	vertex[vanz].cpl[0]  = kcpl0;
	vertex[vanz].cpl[1]  = kcpl1;
	vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
	vertex[vanz].on      = 1;
	
	
	if (fl1.Strong()) {  
	  vertex[vanz].Color.push_back(Color_Function(cf::D));     
	  vertex[vanz].Color.back().SetParticleArg(0,2);     
	  vertex[vanz].Color.back().SetStringArg('0','2');     
	}
	else
	  vertex[vanz].Color.push_back(Color_Function(cf::None));
	
	
	vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("FFS",LF_Key()));     
	vertex[vanz].Lorentz.back()->SetParticleArg(1);     
	
	vertex.push_back(Single_Vertex());vanz++;
	//checked FK & RK	
      }
    }
    
    // q(u) -> A0 + q(u)
    
    for(short int z=1;z<7;z++) { 
      Flavour fl1 = Flavour((kf_code)(z)); 
      if(fl1.IsOn() && fl1.IsQuark() && fl1.IsUptype() && (fl1.Yuk() > 1.)) {
	
	vertex[vanz].in[0] = fl1; 
	vertex[vanz].in[1] = flA0;
	vertex[vanz].in[2] = fl1; 

	kcpl0 = -K_yuk(fl1,flA0)/v2*K_Z_H(1,0);
	kcpl1 = -kcpl0;

	vertex[vanz].cpl[0]  = kcpl0;
	vertex[vanz].cpl[1]  = kcpl1;
	vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
	vertex[vanz].on      = 1;
	
	
	vertex[vanz].Color.push_back(Color_Function(cf::D));     
	vertex[vanz].Color.back().SetParticleArg(0,2);     
	vertex[vanz].Color.back().SetStringArg('0','2');     
      
	
	vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("FFS",LF_Key()));    
	vertex[vanz].Lorentz.back()->SetParticleArg(1);     

	vertex.push_back(Single_Vertex());vanz++;
	//checked FK & RK
      }
    }
  } 
  
  if(flHplus.IsOn()) {
    
    // l(e) -> H- + nu(e)
    
    Kabbala kcpl0,kcpl1;

    for(short int i=11;i<17;i++) {
      Flavour fl1 = Flavour((kf_code)(i)); 
      if(fl1.IsOn() && fl1.IsLepton() && fl1.IsDowntype() && (fl1.Yuk() > 1.)) {
	Flavour fl2 = Flavour((kf_code)(i=i+1));
	if(fl2.IsOn() && fl2.IsLepton() && fl2.IsUptype() ) {
	  
	  vertex[vanz].in[0] = fl2;   
	  vertex[vanz].in[1] = flHplus;
	  vertex[vanz].in[2] = fl1;   
	  
	  kcpl1 = M_I/v1*root2*K_yuk(fl1,flHplus)*K_Z_H(0,0);
	  kcpl0 = K_zero;	 	  

	  vertex[vanz].cpl[0]  = kcpl0;
 	  vertex[vanz].cpl[1]  = kcpl1;
	  vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();	  
	  vertex[vanz].on      = 1;

	  
	  vertex[vanz].Color.push_back(Color_Function(cf::None));     
	  
	  
	  vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("FFS",LF_Key()));     
	  vertex[vanz].Lorentz.back()->SetParticleArg(1);     
	  
	  vertex.push_back(Single_Vertex());vanz++;
	  //checked FK & RK	  
	}
      }
    }  
    
    // q(d) -> H- + q(u) 
    
    for(short int j=1;j<7;j++) {
      Flavour fl1 = Flavour((kf_code)(j));
      Kabbala kcpl0,kcpl1;
      if(fl1.IsOn() && fl1.IsQuark() && fl1.IsDowntype()) {
	for(short int k=1;k<7;k++) {
	  Flavour fl2 = Flavour((kf_code)(k));
	  if(fl2.IsOn() && fl2.IsQuark() && fl2.IsUptype() && 
	     ((fl1.Yuk() > 1.) || (fl2.Yuk() > 1.)) ) {
	    
            int geni=(fl1.Kfcode()-1)/2; //downtype
	    int genj=(fl2.Kfcode()-2)/2; //uptype
	   
	    vertex[vanz].in[0] = fl2;
	    vertex[vanz].in[1] = flHplus;
	    vertex[vanz].in[2] = fl1;
	   
	    kcpl0 = M_I/v2*root2*K_yuk(fl2,flHplus)*K_Z_H(1,0)*K_CKM(geni,genj);
	    kcpl1 = M_I/v1*root2*K_yuk(fl1,flHplus)*K_Z_H(0,0)*K_CKM(geni,genj);
	    	   
	    vertex[vanz].cpl[0]  = kcpl0;
	    vertex[vanz].cpl[1]  = kcpl1;
	    vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();	  	    
	    vertex[vanz].on      = 1;

	    
	    vertex[vanz].Color.push_back(Color_Function(cf::D));     
	    vertex[vanz].Color.back().SetParticleArg(0,2);     
	    vertex[vanz].Color.back().SetStringArg('0','2');     
	    
	    
	    vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("FFS",LF_Key()));     
	    vertex[vanz].Lorentz.back()->SetParticleArg(1);     
	    
	    vertex.push_back(Single_Vertex());vanz++;
	    //checked FK & RK
	  }
	}
      }  
    }      
  }
}

void Interaction_Model_Higgs_THDM::c_VVS(std::vector<Single_Vertex>& vertex,int& vanz) 
{
    
  Flavour flWplus(kf_Wplus);
  Flavour flZ(kf_Z);
  
  Kabbala kcpl0,kcpl1;
  
  //W- -> h0,H0 + W-

  if(flWplus.IsOn()) {
    
    for(short int i=25; i<36;i+=10){
      Flavour flav = Flavour((kf_code)(i));
      if(flav.IsOn()) {
	
	vertex[vanz].in[0] = flWplus.Bar();
	vertex[vanz].in[1] = flav;
	vertex[vanz].in[2] = flWplus.Bar();

	kcpl0 = M_I/num_2*g2*g2*(v1*K_Z_R(0,(i-25)/10)+v2*K_Z_R(1,(i-25)/10));
	kcpl1 = kcpl0;

	vertex[vanz].cpl[0]  = kcpl0;
	vertex[vanz].cpl[1]  = kcpl1;
	vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();	  	    	
	vertex[vanz].on      = 1;
	
	
	vertex[vanz].Color.push_back(Color_Function(cf::None));     
	
	
	vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("Gab",LF_Key()));
	vertex[vanz].Lorentz.back()->SetParticleArg(0,2);     
	
	vertex.push_back(Single_Vertex());vanz++;
	//checked FK & RK
      }
    }
  }
  
  //  Z -> h0,H0 + Z  
  
  if(flZ.IsOn()){
    for(short int i=25; i<36;i+=10){
      Flavour flav = Flavour((kf_code)(i));
      if(flav.IsOn()) {
	
 	vertex[vanz].in[0] = flZ;
	vertex[vanz].in[1] = flav;
	vertex[vanz].in[2] = flZ;
	
	kcpl0 = M_I/(costW*costW*num_2)*g2*g2*(v1*K_Z_R(0,(i-25)/10)+v2*K_Z_R(1,(i-25)/10));
	kcpl1 = kcpl0;

	vertex[vanz].cpl[0]  = kcpl0;
	vertex[vanz].cpl[1]  = kcpl1;
	vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();	  	    	
	vertex[vanz].on      = 1;
	
	
	vertex[vanz].Color.push_back(Color_Function(cf::None));     
	
	
	vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("Gab",LF_Key()));
	vertex[vanz].Lorentz.back()->SetParticleArg(0,2);     
	
	vertex.push_back(Single_Vertex());vanz++;
	//checked FK & RK & SS
      }
    }  
  }
}


void Interaction_Model_Higgs_THDM::c_SSS(std::vector<Single_Vertex>& vertex,int& vanz) 
{
  Flavour flHplus(kf_Hplus);    
  Flavour flh0(kf_h0);
  Flavour flH0(kf_H0);
  Flavour flA0(kf_A0);
  Kabbala kcpl0,kcpl1;

  // h0,H0 -> H+ + H-
  
  if(flHplus.IsOn()) {  
    for(short int i=25; i<36;i+=10){
      Flavour flav = Flavour((kf_code)(i)); 
      if(flav.IsOn()) {
	
	vertex[vanz].in[0] = flHplus.Bar();
	vertex[vanz].in[1] = flav;
	vertex[vanz].in[2] = flHplus.Bar();
	
	kcpl0 = -M_I*g2*g2*(K_A_H(0,0)*K_B_R((i-25)/10)/(costW*costW*num_4)+
			    vev*K_A_P((i-25)/10,0)/num_2);
	kcpl1 = kcpl0;
	
	vertex[vanz].cpl[0]  = kcpl0;
	vertex[vanz].cpl[1]  = kcpl1;
	vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();	  	    		
	
	
	vertex[vanz].Color.push_back(Color_Function(cf::None)); 
	
	
	vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("SSS",LF_Key()));
	
	vertex[vanz].on      = 1;
	vertex.push_back(Single_Vertex());vanz++;
	//checked FK & RK
      }  
      
    }
  }
  
  // h0 -> h0/H0 + h0/H0 
  
  if(flh0.IsOn()) { 
    for(short int j=25; j<36;j+=10){
      Flavour flav2 = Flavour((kf_code)(j));
      if(flav2.IsOn()) {
	for(short int k=25; k<36;k+=10){
	  Flavour flav3 = Flavour((kf_code)(k));
	  if(flav3.IsOn()) {
	    if(flav2 != flH0 || flav3 != flh0) {
	      
	      vertex[vanz].in[0] = flh0;
	      vertex[vanz].in[1] = flav2;      
	      vertex[vanz].in[2] = flav3;
	      
	      kcpl0 = -M_I*g2*g2/(costW*costW*num_4)*(K_A_R(0,(j-25)/10)*K_B_R((k-25)/10) +
						      K_A_R((j-25)/10,(k-25)/10)*K_B_R(0) +
						      K_A_R((k-25)/10,0)*K_B_R((j-25)/10));
 
	      kcpl1 = kcpl0;

	      vertex[vanz].cpl[0]  = kcpl0;
	      vertex[vanz].cpl[1]  = kcpl1;
	      vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();	  	    		
	
	      
	      vertex[vanz].Color.push_back(Color_Function(cf::None)); 
    
	      
	      vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("SSS",LF_Key()));
	      
	      vertex[vanz].on      = 1;
	      vertex.push_back(Single_Vertex());vanz++;
	      //checked FK & RK
	      }
	    }
	  }
	} 
      }
    }
    //H0 -> H0 + H0
    
    if(flH0.IsOn()) { 
      vertex[vanz].in[0] = flH0;
      vertex[vanz].in[1] = flH0;      
      vertex[vanz].in[2] = flH0;
      
      kcpl0 = -M_I*g2*g2/(costW*costW*num_4)*(num_3*K_A_R(1,1)*K_B_R(1));
      kcpl1 = kcpl0;
      
      vertex[vanz].cpl[0]  = kcpl0;
      vertex[vanz].cpl[1]  = kcpl1;
      vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();	  	    		
      
      
      vertex[vanz].Color.push_back(Color_Function(cf::None)); 
      
      
      vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("SSS",LF_Key()));
      
      vertex[vanz].on      = 1;
      vertex.push_back(Single_Vertex());vanz++;
    }
      
  // A0 -> h0/H0 + A0
    
    if(flA0.IsOn()) {
      for(short int i=25;i<36;i+=10) {
	Flavour flav = Flavour((kf_code)(i));
	if(flav.IsOn()) {
	  
	vertex[vanz].in[0] = flA0;
	vertex[vanz].in[1] = flav;
	vertex[vanz].in[2] = flA0;
	
	kcpl0 = -M_I/(costW*costW*num_4)*g2*g2*K_A_H(0,0)*K_B_R((i-25)/10);
	kcpl1 = kcpl0;
	
	vertex[vanz].cpl[0]  = kcpl0;
	vertex[vanz].cpl[1]  = kcpl1;
	vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();	  	    		
	
	
	vertex[vanz].Color.push_back(Color_Function(cf::None)); 
	
	
	vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("SSS",LF_Key()));
	
	vertex[vanz].on      = 1;
	vertex.push_back(Single_Vertex());vanz++;
	//checked FK & RK	
      }
    }
  }
 
}


void Interaction_Model_Higgs_THDM::c_SSV(std::vector<Single_Vertex>& vertex,int& vanz)
{   
  Flavour flWplus(kf_Wplus);
  Flavour flZ(kf_Z);
  Flavour flHplus(kf_Hplus);
  Flavour flA0(kf_A0);
  Flavour flPhoton(kf_photon);
  Kabbala kcpl0,kcpl1;

  // A0 -> Z + h0/H0 
  
  if(flZ.IsOn()) {
    for(short int i=25; i<36;i+=10) {
      Flavour flav = Flavour((kf_code)(i)); 
      if(flav.IsOn() && flA0.IsOn()) {
	vertex[vanz].in[0] = flA0;
	vertex[vanz].in[1] = flZ;
	vertex[vanz].in[2] = flav;

	kcpl0 = g2/(costW*num_2)*K_A_M((i-25)/10,0);
	kcpl1 = kcpl0;

	vertex[vanz].cpl[0]  = kcpl0;
	vertex[vanz].cpl[1]  = kcpl1;
	vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();	  	    		
	
	
	vertex[vanz].Color.push_back(Color_Function(cf::None)); 
 
	
	vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("SSV",LF_Key()));
	vertex[vanz].Lorentz.back()->SetParticleArg(0,2,1);     	

	vertex[vanz ].on      = 1;
	vertex.push_back(Single_Vertex());vanz++;
	//checked FK & RK
      }  
    }
    
    // H- -> Z + H-
    
    if(flHplus.IsOn()) {
      vertex[vanz].in[0] = flHplus.Bar();
      vertex[vanz].in[1] = flZ;
      vertex[vanz].in[2] = flHplus.Bar();

      Complex cos2tw = (costW.Value()*costW.Value()-sintW.Value()*sintW.Value())/
	(2.*sintW.Value()*costW.Value());
      Kabbala cot2TW = Kabbala(string("cot2\\theta_W"),cos2tw);

      kcpl0 = M_I*g1*cot2TW;
      kcpl1 = kcpl0;

      vertex[vanz].cpl[0]  = kcpl0;
      vertex[vanz].cpl[1]  = kcpl1;
      vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();	  	    		
      
      
      vertex[vanz].Color.push_back(Color_Function(cf::None)); 
      
      
      vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("SSV",LF_Key()));
      vertex[vanz].Lorentz.back()->SetParticleArg(0,2,1);     	
      
      vertex[vanz].on      = 1;
      vertex.push_back(Single_Vertex());vanz++;
      //checked FK & RK
    }    
  }
  
  // H- -> Photon + H- 
  
  if(flPhoton.IsOn() && flHplus.IsOn()) {
    
    vertex[vanz].in[0] = flHplus.Bar();
    vertex[vanz].in[1] = flPhoton;
    vertex[vanz].in[2] = flHplus.Bar();
    
    kcpl0 = M_I*g1;
    kcpl1 = kcpl0;

    vertex[vanz].cpl[0]  = kcpl0;
    vertex[vanz].cpl[1]  = kcpl1;
    vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();	  	    		
    
    
    vertex[vanz].Color.push_back(Color_Function(cf::None)); 
    
    
    vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("SSV",LF_Key()));
    vertex[vanz].Lorentz.back()->SetParticleArg(0,2,1);     	
          
    vertex[vanz].on      = 1;
    vertex.push_back(Single_Vertex());vanz++;
    //checked FK & RK
  }    
  
  //H- -> W- + h0/H0  
  
  if(flWplus.IsOn()) {
    for(short int i=25; i<36;i+=10) {
      Flavour flav = Flavour((kf_code)(i)); 
      
      if(flav.IsOn() && flHplus.IsOn()) {
	
	vertex[vanz].in[0] = flHplus.Bar();
	vertex[vanz].in[1] = flWplus.Bar();  
	vertex[vanz].in[2] = flav;
	
	kcpl0 = -(M_I/num_2)*g2*K_A_M((i-25)/10,0);
	kcpl1 = kcpl0;

	vertex[vanz].cpl[0]  = kcpl0;
	vertex[vanz].cpl[1]  = kcpl1;
	vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();	  	    		

	
	vertex[vanz].Color.push_back(Color_Function(cf::None)); 
 
	
	vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("SSV",LF_Key()));
	vertex[vanz].Lorentz.back()->SetParticleArg(0,2,1);     	
	
	vertex[vanz].on      = 1;
	vertex.push_back(Single_Vertex());vanz++;
	//checked FK & RK
      }
    }
    
    //H- -> W- +  A0
    
    if(flHplus.IsOn() && flA0.IsOn()) {
      vertex[vanz].in[0] = flHplus.Bar();
      vertex[vanz].in[1] = flWplus.Bar();    
      vertex[vanz].in[2] = flA0;

      kcpl0 = -g2/num_2;
      kcpl1 = kcpl0;

      vertex[vanz].cpl[0]  = kcpl0;
      vertex[vanz].cpl[1]  = kcpl1;
      vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();	  	    		
      
      
      vertex[vanz].Color.push_back(Color_Function(cf::None)); 
 
      
      vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("SSV",LF_Key()));
      vertex[vanz].Lorentz.back()->SetParticleArg(0,2,1);     	
      
      vertex[vanz].on      = 1;
      vertex.push_back(Single_Vertex());vanz++;
      //checked FK & RK
    }
  } 
}

void Interaction_Model_Higgs_THDM::c_SSSS(std::vector<Single_Vertex>& vertex,int& vanz) 
{
  Flavour flHplus(kf_Hplus);    
  Flavour flA0(kf_A0);
  
  Kabbala kcpl0,kcpl1;

  // A0 - A0 - A0 - A0
  if (flA0.IsOn()) {
 	
    vertex[vanz].in[0] = flA0;
    vertex[vanz].in[1] = flA0;
    vertex[vanz].in[2] = flA0;
    vertex[vanz].in[3] = flA0;
    
    vertex[vanz].nleg     = 4;
    
    kcpl0 = -M_I*g2*g2/(num_4*costW*costW)*num_3*K_A_H(0,0)*K_A_H(0,0);
    kcpl1 = kcpl0;
    
    vertex[vanz].cpl[0]  = kcpl0;
    vertex[vanz].cpl[1]  = kcpl1;
    vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
    
    
    vertex[vanz].Color.push_back(Color_Function(cf::None));     
    
    
    vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("SSSS",LF_Key()));     
    
    vertex[vanz].on      = 1;
    vertex.push_back(Single_Vertex());vanz++;
  }
  // H- - H- - H+ - H-
  if (flHplus.IsOn()) {
 	
    vertex[vanz].in[0] = flHplus.Bar();
    vertex[vanz].in[1] = flHplus.Bar();
    vertex[vanz].in[2] = flHplus;
    vertex[vanz].in[3] = flHplus.Bar();
    
    vertex[vanz].nleg     = 4;
    
    kcpl0 = -M_I*g2*g2/(num_4*costW*costW)*num_2*K_A_H(0,0)*K_A_H(0,0);
    kcpl1 = kcpl0;
    
    vertex[vanz].cpl[0]  = kcpl0;
    vertex[vanz].cpl[1]  = kcpl1;
    vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
    
    
    vertex[vanz].Color.push_back(Color_Function(cf::None));     
    
    
    vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("SSSS",LF_Key()));     
    
    vertex[vanz].on      = 1;
    vertex.push_back(Single_Vertex());vanz++;
  }
  // A0 - H- - H+ - A0
  if (flA0.IsOn() && flHplus.IsOn()) {
    vertex[vanz].in[0] = flA0;
    vertex[vanz].in[1] = flHplus.Bar();
    vertex[vanz].in[2] = flHplus;
    vertex[vanz].in[3] = flA0;
    
    vertex[vanz].nleg     = 4;
    
    kcpl0 = -M_I*g2*g2/(num_4*costW*costW)*K_A_H(0,0)*K_A_H(0,0);
    kcpl1 = kcpl0;
    
    vertex[vanz].cpl[0]  = kcpl0;
    vertex[vanz].cpl[1]  = kcpl1;
    vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
    
    
    vertex[vanz].Color.push_back(Color_Function(cf::None));     
    
    
    vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("SSSS",LF_Key()));     
    
    vertex[vanz].on      = 1;
    vertex.push_back(Single_Vertex());vanz++;
  }

  for (int i=25;i<36;i+=10) {
    Flavour flav1 = Flavour((kf_code)(i));
    for (int j=i;j<36;j+=10) {
      Flavour flav2 = Flavour((kf_code)(j));
      if (flav1.IsOn() && flav2.IsOn()) {
	// h0/H0 - H- - H+ - h0/H0
	if (flHplus.IsOn()) {
	  
	  vertex[vanz].in[0] = flav1;
	  vertex[vanz].in[1] = flHplus.Bar();
	  vertex[vanz].in[2] = flHplus;
	  vertex[vanz].in[3] = flav2;
	  
	  vertex[vanz].nleg     = 4;
	  
	  kcpl0 = -M_I*g2*g2/num_4*(K_A_R((i-25)/10,(j-25)/10)*K_A_H(0,0)/costW*costW +
				    num_2*K_A_P((i-25)/10,0)*K_A_P((j-25)/10,0));
	  
	  kcpl1 = kcpl0;
	  
	  vertex[vanz].cpl[0]  = kcpl0;
	  vertex[vanz].cpl[1]  = kcpl1;
	  vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
	  
	  
	  vertex[vanz].Color.push_back(Color_Function(cf::None));     
	  
	  
	  vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("SSSS",LF_Key()));     
	  
	  vertex[vanz].on      = 1;
	  vertex.push_back(Single_Vertex());vanz++;
	}
	// A0 - h0/H0 - h0/H0 - A0
	if (flA0.IsOn()) {
	  
	  vertex[vanz].in[0] = flA0;
	  vertex[vanz].in[1] = flav1;
	  vertex[vanz].in[2] = flav2;
	  vertex[vanz].in[3] = flA0;
	  
	  vertex[vanz].nleg     = 4;
	  
	  kcpl0 = -M_I*g2*g2/(num_4*costW*costW)*K_A_R((i-25)/10,(j-25)/10)*K_A_H(0,0);
	  
	  kcpl1 = kcpl0;
	  
	  vertex[vanz].cpl[0]  = kcpl0;
	  vertex[vanz].cpl[1]  = kcpl1;
	  vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
	  
	  
	  vertex[vanz].Color.push_back(Color_Function(cf::None));     
	  
	  
	  vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("SSSS",LF_Key()));     
	  
	  vertex[vanz].on      = 1;
	  vertex.push_back(Single_Vertex());vanz++;
	}
	// h0/H0 - h0/H0 - h0/H0 - h0/H0 
	for (int k=25;k<36;k+=10) {
	  Flavour flav3 = Flavour((kf_code)(k));
	  for (int l=k;l<36;l+=10) {
	    Flavour flav4 = Flavour((kf_code)(l));
	    if (flav3.IsOn() && flav4.IsOn()) {
	      
	      vertex[vanz].in[0] = flav1;
	      vertex[vanz].in[1] = flav2;
	      vertex[vanz].in[2] = flav3;
	      vertex[vanz].in[3] = flav4;
	      
	      vertex[vanz].nleg     = 4;
	      
	      kcpl0 = -M_I*g2*g2/(num_4*costW*costW)*(K_A_R((i-25)/10,(j-25)/10)*K_A_R((k-25)/10,(l-25)/10)+
						      K_A_R((i-25)/10,(k-25)/10)*K_A_R((j-25)/10,(l-25)/10)+
						      K_A_R((i-25)/10,(l-25)/10)*K_A_R((j-25)/10,(k-25)/10));
	      
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
	}
      }
    }
  }
}

void Interaction_Model_Higgs_THDM::c_SSVV(std::vector<Single_Vertex>& vertex,int& vanz) 
{
    
  Flavour flWplus(kf_Wplus);
  Flavour flZ(kf_Z);
  Flavour flHplus(kf_Hplus);    
  Flavour flA0(kf_A0);
  Flavour flPhoton(kf_photon);
  
  Kabbala kcpl0,kcpl1;

  // Z - h0/H0/A0 - h0/H0/A0 - Z  
  if (flZ.IsOn()) {
    for(short int i=25; i<37;++i){
      if (i==26) i=35;
      Flavour flav = Flavour((kf_code)(i));
      if(flav.IsOn()) {
	
 	vertex[vanz].in[0] = flZ;
	vertex[vanz].in[1] = flav;
	vertex[vanz].in[2] = flav;
	vertex[vanz].in[3] = flZ;
  
	vertex[vanz].nleg     = 4;

	kcpl0 = M_I*g2*g2/(num_2*costW*costW);
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
  }
 // W - h0/H0/A0 - h0/H0/A0 - W  
  if (flZ.IsOn()) {
    for(short int i=25; i<37;++i){
      if (i==26) i=35;
      Flavour flav = Flavour((kf_code)(i));
      if(flav.IsOn()) {
	
 	vertex[vanz].in[0] = flWplus.Bar();
	vertex[vanz].in[1] = flav;
	vertex[vanz].in[2] = flav;
	vertex[vanz].in[3] = flWplus.Bar();
  
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
  }

  if (flWplus.IsOn() && flHplus.IsOn()) {
    //W -> H- - H- - W
    vertex[vanz].in[0] = flWplus.Bar();
    vertex[vanz].in[1] = flHplus;
    vertex[vanz].in[2] = flHplus.Bar();
    vertex[vanz].in[3] = flWplus.Bar();
    
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
    
    // W+ -> H+ + h0/H0 + Z/P 
    for (short int i=25;i<36;i+=10){
      Flavour flav = Flavour((kf_code)(i));
      if(flav.IsOn()) {
	if(flZ.IsOn()) {

	  vertex[vanz].in[0] = flWplus;
	  vertex[vanz].in[1] = flHplus;
	  vertex[vanz].in[2] = flav;
	  vertex[vanz].in[3] = flZ;
	  
	  vertex[vanz].nleg     = 4;
	  
	  kcpl0 = M_I*g1*g1/(num_2*costW)*K_A_M((i-25)/10,0);
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
	
	if(flPhoton.IsOn()) {

	  vertex[vanz].in[0] = flWplus;
	  vertex[vanz].in[1] = flHplus;
	  vertex[vanz].in[2] = flav;
	  vertex[vanz].in[3] = flPhoton;
  
	  vertex[vanz].nleg     = 4;
	  
	  kcpl0 = -M_I*g1*g1/(num_2*sintW)*K_A_M((i-25)/10,0);
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
    }
    // W+ -> H+ + A0 + Z/P   
    if (flA0.IsOn()) {
      if (flZ.IsOn()) {
	
 	vertex[vanz].in[0] = flWplus;
	vertex[vanz].in[1] = flHplus;
	vertex[vanz].in[2] = flA0;
	vertex[vanz].in[3] = flZ;
  
	vertex[vanz].nleg     = 4;
	
	kcpl0 = g1*g1/(num_2*costW);
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
      if (flPhoton.IsOn()) {
	
	vertex[vanz].in[0] = flWplus;
	vertex[vanz].in[1] = flHplus;
	vertex[vanz].in[2] = flA0;
	vertex[vanz].in[3] = flPhoton;
  
	vertex[vanz].nleg     = 4;
	
	kcpl0 = -g1*g1/(num_2*sintW);
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
  }
  //P/Z -> H- - H- - Z/P 
  if (flHplus.IsOn()) {
    if (flZ.IsOn()) {
      
      Kabbala cot2TW = Kabbala(string("cot2\\theta_W"),
			       (costW.Value()*costW.Value()-sintW.Value()*sintW.Value())/
			       (2.*sintW.Value()*costW.Value()));
      
      vertex[vanz].in[0] = flZ;
      vertex[vanz].in[1] = flHplus.Bar();
      vertex[vanz].in[2] = flHplus;
      vertex[vanz].in[3] = flZ;
      
      vertex[vanz].nleg     = 4;

      kcpl0 = M_I*g1*g1*num_2*cot2TW*cot2TW;
      kcpl1 = kcpl0;
      
      vertex[vanz].cpl[0]  = kcpl0;
      vertex[vanz].cpl[1]  = kcpl1;
      vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
      
      
      vertex[vanz].Color.push_back(Color_Function(cf::None));     
      
      
      vertex[vanz].Lorentz.push_back(LF_Getter::GetObject("VVSS",LF_Key()));
      vertex[vanz].Lorentz.back()->SetParticleArg(0,3);     
      
      vertex[vanz].on      = 1;
      vertex.push_back(Single_Vertex());vanz++;
      
      if(flPhoton.IsOn()) {

	vertex[vanz].in[0] = flZ;
	vertex[vanz].in[1] = flHplus.Bar();
	vertex[vanz].in[2] = flHplus;
	vertex[vanz].in[3] = flPhoton;
	
	vertex[vanz].nleg     = 4;
	
	kcpl0 = M_I*g1*g1*num_2*cot2TW;
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
    if(flPhoton.IsOn()) {
       vertex[vanz].in[0] = flPhoton;
       vertex[vanz].in[1] = flHplus.Bar();
       vertex[vanz].in[2] = flHplus;
       vertex[vanz].in[3] = flPhoton;
       
      vertex[vanz].nleg     = 4;
      
      kcpl0 = M_I*g1*g1*num_2;
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
}


Kabbala Interaction_Model_Higgs_THDM::K_CKM(short int i,short int j)       
{   
  char hi[2],hj[2];
  sprintf(hi,"%i",i);
  sprintf(hj,"%i",j);
  return Kabbala(string("V_{")+string(hi)+string(hj)+string("}"),
		 ComplexMatrixElement(std::string("CKM"),i,j));
} 
  
Kabbala Interaction_Model_Higgs_THDM::conj_K_CKM(short int i,short int j)       
{   
  char hi[2],hj[2];
  sprintf(hi,"%i",i);
  sprintf(hj,"%i",j);
  return Kabbala(string("V^\\ti_{")+string(hi)+string(hj)+string("}"),
		 conj(ComplexMatrixElement(std::string("CKM"),i,j)));
} 
 
Kabbala Interaction_Model_Higgs_THDM::K_yuk(Flavour fl,Flavour H) {
  std::string flname = std::string("m")+std::string(fl.IDName());
  double meff(p_model->ScalarFunction(flname,sqr(H.Mass())));
  if(ScalarNumber(std::string("WidthScheme"))==0){
    return Kabbala(string("M_{"+fl.TexName()+"}(m_"+H.IDName()+")"),meff);
  }else{
    return Kabbala(string("M_{"+fl.TexName()+"}(m_"+H.IDName()+")"),
                   sqrt(sqr(meff)-Complex(0.,1.)*fl.Width()*meff));
  }
}

Kabbala Interaction_Model_Higgs_THDM::K_yuk_sign(Flavour fl) {
  char hi[3];
  sprintf(hi,"%i",fl.MassSign());
  return Kabbala(string(hi),fl.MassSign());
}

Kabbala Interaction_Model_Higgs_THDM::K_Z_H(short int i,short int j) {   
  char hi[2],hj[2];
  sprintf(hi,"%i",i);
  sprintf(hj,"%i",j);
  return Kabbala(string("Z^{")+string(hi)+string(hj)+string("}_H"),
		 ComplexMatrixElement(std::string("Z_H"),i,j));
}     

Kabbala Interaction_Model_Higgs_THDM::K_Z_R(short int i,short int j) {   
  char hi[2],hj[2];
  sprintf(hi,"%i",i);
  sprintf(hj,"%i",j);
  return Kabbala(string("Z^{")+string(hi)+string(hj)+string("}_R"),
		 ComplexMatrixElement(std::string("Z_R"),i,j));
}  

Kabbala Interaction_Model_Higgs_THDM::K_A_M(short int i,short int j) {
  char hi[2],hj[2];
  sprintf(hi,"%i",i);
  sprintf(hj,"%i",j);
  return Kabbala(string("A^{")+string(hi)+string(hj)+string("}_M"),
		 ComplexMatrixElement(std::string("Z_R"),0,i) *
		 ComplexMatrixElement(std::string("Z_H"),0,j) -
		 ComplexMatrixElement(std::string("Z_R"),1,i) *
		 ComplexMatrixElement(std::string("Z_H"),1,j));
}  

Kabbala Interaction_Model_Higgs_THDM::K_A_H(short int i,short int j) {
  char hi[2],hj[2];
  sprintf(hi,"%i",i);
  sprintf(hj,"%i",j);
  return Kabbala(string("A^{")+string(hi)+string(hj)+string("}_H"),
		 ComplexMatrixElement(std::string("Z_H"),0,i) *
		 ComplexMatrixElement(std::string("Z_H"),0,j) -
		 ComplexMatrixElement(std::string("Z_H"),1,i) *
		 ComplexMatrixElement(std::string("Z_H"),1,j));
}  

Kabbala Interaction_Model_Higgs_THDM::K_A_P(short int i,short int j) {
  char hi[2],hj[2];
  sprintf(hi,"%i",i);
  sprintf(hj,"%i",j);
  return Kabbala(string("A^{")+string(hi)+string(hj)+string("}_P"),
		 ComplexMatrixElement(std::string("Z_R"),0,i) *
		 ComplexMatrixElement(std::string("Z_H"),1,j) +
		 ComplexMatrixElement(std::string("Z_R"),1,i) *
		 ComplexMatrixElement(std::string("Z_H"),0,j));
}  

Kabbala Interaction_Model_Higgs_THDM::K_A_R(short int i,short int j) {
  char hi[2],hj[2];
  sprintf(hi,"%i",i);
  sprintf(hj,"%i",j);
  return Kabbala(string("A^{")+string(hi)+string(hj)+string("}_R"),
		 ComplexMatrixElement(std::string("Z_R"),0,i) *
		 ComplexMatrixElement(std::string("Z_R"),0,j) -
		 ComplexMatrixElement(std::string("Z_R"),1,i) *
		 ComplexMatrixElement(std::string("Z_R"),1,j));
}  

Kabbala Interaction_Model_Higgs_THDM::K_B_R(short int i) {   
  char hi[2];
  sprintf(hi,"%i",i);
  return Kabbala(string("B^{")+string(hi)+string("}_R"),
		 v1.Value() * ComplexMatrixElement(std::string("Z_R"),0,i) -
		 v2.Value() * ComplexMatrixElement(std::string("Z_R"),1,i) );
}  





















