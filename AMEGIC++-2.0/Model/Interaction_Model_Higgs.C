#include "Interaction_Model_Higgs.H"
#include "MathTools.H"
#include "Message.H"
#include <stdio.h>


using namespace AMEGIC;
using namespace ATOOLS;
using namespace std;

Interaction_Model_Higgs::Interaction_Model_Higgs(MODEL::Model_Base * _model,
						 std::string _cplscheme,std::string _yukscheme) :
  Interaction_Model_Base(_model,_cplscheme,_yukscheme)
{ 
  double Ecms2 = sqr(rpa.gen.Ecms());

  g1     = Kabbala(string("g_1"),
		   sqrt(4.*M_PI*ScalarFunction(std::string("alpha_QED"),Ecms2)));
  g2     = Kabbala(string("g_1/\\sin\\theta_W"), 
		   g1.Value()/sqrt(p_model->ScalarConstant(std::string("sin2_thetaW"))));
  sintW  = Kabbala(std::string("\\sin\\theta_W"),
		   sqrt(p_model->ScalarConstant(std::string("sin2_thetaW"))));
  costW  = Kabbala(std::string("\\cos\\theta_W"),
		   sqrt(1.-p_model->ScalarConstant(std::string("sin2_thetaW"))));
  PL     = Kabbala(string("P_L"),1.);
  PR     = Kabbala(string("P_R"),1.);
  M_I    = Kabbala(string("i"),Complex(0.,1.));
  vev    = Kabbala(string("v_{EW}"),p_model->ScalarConstant(std::string("vev")));
  v1     = Kabbala(string("v_1"),
		   p_model->ScalarConstant(std::string("vev")) *
		   sqrt(1./(1.+sqr(p_model->ScalarConstant(std::string("tan(beta)"))))));
  v2     = Kabbala(string("v_2"),
		   p_model->ScalarConstant(std::string("vev")) *
		   p_model->ScalarConstant(std::string("tan(beta)")) *
		   sqrt(1./(1.+sqr(p_model->ScalarConstant(std::string("tan(beta)"))))));
  root2  = Kabbala(string("\\sqrt{2}"),sqrt(2.));
  K_zero = Kabbala(string("zero"),0.);
  num_2  = Kabbala(string("2"),2.);    	
  num_3  = Kabbala(string("3"),3.);    	
  num_4  = Kabbala(string("4"),4.);    		
}


void Interaction_Model_Higgs::c_FFS(Single_Vertex* vertex,int& vanz)
{
  Flavour flHmin(kf::Hmin);
  Flavour flA0(kf::A0);
  
  // l(e)/q(d) -> h0/H0 + l(e)/q(d)
  
  for(short int i=31;i<33;i++) {
    
    Flavour flav = Flavour(kf::code(i));
    Kabbala kcpl0,kcpl1;

    if(flav.IsOn()) {
      
      for(short int j=1;j<17;j++) {
	if (j==7) j=11;
	Flavour fl1 = Flavour(kf::code(j)); 
       	if(fl1.IsOn() && fl1.IsFermion() && fl1.IsDowntype() && (fl1.Yuk() > 1.)) {
	  
	  vertex[vanz].in[0] = fl1; 
	  vertex[vanz].in[1] = flav;
	  vertex[vanz].in[2] = fl1; 
	  
	  kcpl0 = -M_I/v1*K_yuk(fl1)*K_Z_R(0,i-31);
	  kcpl1 = kcpl0;

	  vertex[vanz].cpl[0]  = kcpl0.Value();
	  vertex[vanz].cpl[1]  = kcpl1.Value();
	  vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
	  vertex[vanz].cpl[2]  = 0.;vertex[vanz].cpl[3]  = 0.;
	  vertex[vanz].on      = 1;

	  vertex[vanz].ncf   = 1;
 	  if (fl1.Strong()) {  
	    vertex[vanz].Color = new Color_Function(cf::D);     
	    vertex[vanz].Color->SetParticleArg(0,2);     
	    vertex[vanz].Color->SetStringArg('0','2');     
	  }
	  else 
	    vertex[vanz].Color = new Color_Function(cf::None);
	
	  vertex[vanz].nlf     = 1;
	  vertex[vanz].Lorentz = new Lorentz_Function(lf::FFS);
	  	  
	  vanz++;
	  //checked FK & RK	  
	}
      }
    }
  
    // q(u) -> h0/H0 + q(u)
    
    for(short int t=1;t<17;t++) {
      if (t==7) t=11;
      Flavour fl1 = Flavour(kf::code(t)); 
      if(fl1.IsOn() && fl1.IsQuark() && fl1.IsUptype() && (fl1.Yuk() > 1.)) {
	
	vertex[vanz].in[0] = fl1; 
	vertex[vanz].in[1] = flav;
	vertex[vanz].in[2] = fl1; 

	kcpl0 = -M_I/v2*K_yuk(fl1)*K_Z_R(1,i-31);
	kcpl1 = kcpl0;

	vertex[vanz].cpl[0]  = kcpl0.Value();
	vertex[vanz].cpl[1]  = kcpl1.Value();
	vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
	vertex[vanz].cpl[2]  = 0.;vertex[vanz].cpl[3]  = 0.;
	vertex[vanz].on      = 1;
	
	vertex[vanz].ncf   = 1;
	vertex[vanz].Color = new Color_Function(cf::D);     
	vertex[vanz].Color->SetParticleArg(0,2);     
	vertex[vanz].Color->SetStringArg('0','2');     
	
	vertex[vanz].nlf     = 1;
	vertex[vanz].Lorentz = new Lorentz_Function(lf::FFS);
	vertex[vanz].Lorentz->SetParticleArg(1);     
	
	vanz++;
	//checked FK & RK
      }
    }
  }
  
  // l(e)/q(d) -> A0 + l(e)/q(d)
  
  if(flA0.IsOn()) {

    Kabbala kcpl0,kcpl1;
 
    for(short int k=1;k<17;k++) {
      if (k==7) k=11;
      Flavour fl1 = Flavour(kf::code(k)); 
      if(fl1.IsOn() && fl1.IsFermion() && fl1.IsDowntype() && (fl1.Yuk() > 1.)) {
		
	vertex[vanz].in[0] = fl1; 
	vertex[vanz].in[1] = flA0;
	vertex[vanz].in[2] = fl1; 


	kcpl0 = -K_yuk(fl1)/v1*K_Z_H(0,0);
	kcpl1 = -kcpl0;

	vertex[vanz].cpl[0]  = kcpl0.Value();
	vertex[vanz].cpl[1]  = kcpl1.Value();
	vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
	vertex[vanz].cpl[2]  = 0.;vertex[vanz].cpl[3]  = 0.;
	vertex[vanz].on      = 1;
	
	vertex[vanz].ncf   = 1;
	if (fl1.Strong()) {  
	  vertex[vanz].Color = new Color_Function(cf::D);     
	  vertex[vanz].Color->SetParticleArg(0,2);     
	  vertex[vanz].Color->SetStringArg('0','2');     
	}
	else
	  vertex[vanz].Color = new Color_Function(cf::None);
	
	vertex[vanz].nlf     = 1;
	vertex[vanz].Lorentz = new Lorentz_Function(lf::FFS);     
	vertex[vanz].Lorentz->SetParticleArg(1);     
	
	vanz++;
	//checked FK & RK	
      }
    }
    
    // q(u) -> A0 + q(u)
    
    for(short int z=1;z<7;z++) { 
      Flavour fl1 = Flavour(kf::code(z)); 
      if(fl1.IsOn() && fl1.IsQuark() && fl1.IsUptype() && (fl1.Yuk() > 1.)) {
	
	vertex[vanz].in[0] = fl1; 
	vertex[vanz].in[1] = flA0;
	vertex[vanz].in[2] = fl1; 

	kcpl0 = -K_yuk(fl1)/v2*K_Z_H(1,0);
	kcpl1 = -kcpl0;

	vertex[vanz].cpl[0]  = kcpl0.Value();
	vertex[vanz].cpl[1]  = kcpl1.Value();
	vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
	vertex[vanz].cpl[2]  = 0.;vertex[vanz].cpl[3]  = 0.;
	vertex[vanz].on      = 1;
	
	vertex[vanz].ncf   = 1;
	vertex[vanz].Color = new Color_Function(cf::D);     
	vertex[vanz].Color->SetParticleArg(0,2);     
	vertex[vanz].Color->SetStringArg('0','2');     
      
	vertex[vanz].nlf     = 1;
	vertex[vanz].Lorentz = new Lorentz_Function(lf::FFS);    
	vertex[vanz].Lorentz->SetParticleArg(1);     

	vanz++;
	//checked FK & RK
      }
    }
  } 
  
  if(flHmin.IsOn()) {
    
    // l(e) -> H- + nu(e)
    
    Kabbala kcpl0,kcpl1;

    for(short int i=11;i<17;i++) {
      Flavour fl1 = Flavour(kf::code(i)); 
      if(fl1.IsOn() && fl1.IsLepton() && fl1.IsDowntype() && (fl1.Yuk() > 1.)) {
	Flavour fl2 = Flavour(kf::code(i=i+1));
	if(fl2.IsOn() && fl2.IsLepton() && fl2.IsUptype() ) {
	  
	  vertex[vanz].in[0] = fl1;   
	  vertex[vanz].in[1] = flHmin;
	  vertex[vanz].in[2] = fl2;   

	  kcpl1 = M_I/v1*root2*K_yuk(fl1)*K_Z_H(0,0);
	  kcpl0 = K_zero;	 	  

	  vertex[vanz].cpl[0]  = kcpl0.Value();
 	  vertex[vanz].cpl[1]  = kcpl1.Value();
	  vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();	  
	  vertex[vanz].cpl[2]  = 0.;vertex[vanz].cpl[3]  = 0.;
	  vertex[vanz].on      = 1;

	  vertex[vanz].ncf   = 1;
	  vertex[vanz].Color = new Color_Function(cf::None);     
	  
	  vertex[vanz].nlf     = 1;
	  vertex[vanz].Lorentz = new Lorentz_Function(lf::FFS);     
	  vertex[vanz].Lorentz->SetParticleArg(1);     
	  
	  vanz++;
	  //checked FK & RK	  
	}
      }
    }  
    
    // q(d) -> H- + q(u) 
    
    for(short int j=1;j<7;j++) {
      Flavour fl1 = Flavour(kf::code(j));
      Kabbala kcpl0,kcpl1;
      if(fl1.IsOn() && fl1.IsQuark() && fl1.IsDowntype()) {
	for(short int k=1;k<7;k++) {
	  Flavour fl2 = Flavour(kf::code(k));
	  if(fl2.IsOn() && fl2.IsQuark() && fl2.IsUptype() && 
	     ((fl1.Yuk() > 1.) || (fl2.Yuk() > 1.)) ) {
	    
            int geni=(fl1.Kfcode()-1)/2; //downtype
	    int genj=(fl2.Kfcode()-2)/2; //uptype
	   
	    vertex[vanz].in[0] = fl1;
	    vertex[vanz].in[1] = flHmin;
	    vertex[vanz].in[2] = fl2;
	   
	    kcpl0 = M_I/v2*root2*K_yuk(fl2)*K_Z_H(1,0)*K_CKM(geni,genj);
	    kcpl1 = M_I/v1*root2*K_yuk(fl1)*K_Z_H(0,0)*K_CKM(geni,genj);
	    	   
	    vertex[vanz].cpl[0]  = kcpl0.Value();
	    vertex[vanz].cpl[1]  = kcpl1.Value();
	    vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();	  	    
	    vertex[vanz].cpl[2]  = 0.;vertex[vanz].cpl[3]  = 0.;
	    vertex[vanz].on      = 1;

	    vertex[vanz].ncf   = 1;
	    vertex[vanz].Color = new Color_Function(cf::D);     
	    vertex[vanz].Color->SetParticleArg(0,2);     
	    vertex[vanz].Color->SetStringArg('0','2');     
	    
	    vertex[vanz].nlf     = 1;
	    vertex[vanz].Lorentz = new Lorentz_Function(lf::FFS);     
	    vertex[vanz].Lorentz->SetParticleArg(1);     
	    
	    vanz++;
	    //checked FK & RK
	  }
	}
      }  
    }      
  }
}

void Interaction_Model_Higgs::c_VVS(Single_Vertex* vertex,int& vanz) 
{
    
  Flavour flW(kf::W);
  Flavour flZ(kf::Z);
  
  Kabbala kcpl0,kcpl1;
  
  //W- -> h0,H0 + W-

  if(flW.IsOn()) {
    
    for(short int i=31; i<33;i++){
      Flavour flav = Flavour(kf::code(i));
      if(flav.IsOn()) {
	
	vertex[vanz].in[0] = flW;
	vertex[vanz].in[1] = flav;
	vertex[vanz].in[2] = flW;

	kcpl0 = M_I/num_2*g2*g2*(v1*K_Z_R(0,i-31)+v2*K_Z_R(1,i-31));
	kcpl1 = kcpl0;

	vertex[vanz].cpl[0]  = kcpl0.Value();
	vertex[vanz].cpl[1]  = kcpl1.Value();
	vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();	  	    	
	vertex[vanz].cpl[2]  = 0.;vertex[vanz].cpl[3]  = 0.;
	vertex[vanz].on      = 1;
	
	vertex[vanz].ncf   = 1;
	vertex[vanz].Color = new Color_Function(cf::None);     
	
	vertex[vanz].nlf     = 1;
	vertex[vanz].Lorentz = new Lorentz_Function(lf::Gab);
	vertex[vanz].Lorentz->SetParticleArg(0,2);     
	
	vanz++;
	//checked FK & RK
      }
    }
  }
  
  //  Z -> h0,H0 + Z  
  
  if(flZ.IsOn()){
    for(short int i=31; i<33;i++){
      Flavour flav = Flavour(kf::code(i));
      if(flav.IsOn()) {
	
 	vertex[vanz].in[0] = flZ;
	vertex[vanz].in[1] = flav;
	vertex[vanz].in[2] = flZ;
	
	kcpl0 = M_I/(costW*costW*num_2)*g2*g2*(v1*K_Z_R(0,i-31)+v2*K_Z_R(1,i-31));
	kcpl1 = kcpl0;

	vertex[vanz].cpl[0]  = kcpl0.Value();
	vertex[vanz].cpl[1]  = kcpl1.Value();
	vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();	  	    	
	vertex[vanz].cpl[2]  = 0.;vertex[vanz].cpl[3]  = 0.;
	vertex[vanz].on      = 1;
	
	vertex[vanz].ncf   = 1;
	vertex[vanz].Color = new Color_Function(cf::None);     
	
	vertex[vanz].nlf     = 1;
	vertex[vanz].Lorentz = new Lorentz_Function(lf::Gab);
	vertex[vanz].Lorentz->SetParticleArg(0,2);     
	
	vanz++;
	//checked FK & RK & SS
      }
    }  
  }
}


void Interaction_Model_Higgs::c_SSS(Single_Vertex* vertex,int& vanz) 
{
  Flavour flHmin(kf::Hmin);    
  Flavour flh0(kf::h0);
  Flavour flH0(kf::H0);
  Flavour flA0(kf::A0);
  Kabbala kcpl0,kcpl1;

  // h0,H0 -> H+ + H-
  
  if(flHmin.IsOn()) {  
    for(short int i=31;i<33;i++) {
      Flavour flav = Flavour(kf::code(i)); 
      if(flav.IsOn()) {
	
	vertex[vanz].in[0] = flHmin;
	vertex[vanz].in[1] = flav;
	vertex[vanz].in[2] = flHmin;
	
	kcpl0 = -M_I*g2*g2*(K_A_H(0,0)*K_B_R(i-31)/(costW*costW*num_4)+
			    vev*K_A_P(i-31,0)/num_2);
	kcpl1 = kcpl0;
	
	vertex[vanz].cpl[0]  = kcpl0.Value();
	vertex[vanz].cpl[1]  = kcpl1.Value();
	vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();	  	    		
	vertex[vanz].cpl[2]  = 0.;vertex[vanz].cpl[3]  = 0.;
	
	vertex[vanz].ncf   = 1;
	vertex[vanz].Color = new Color_Function(cf::None); 
	
	vertex[vanz].nlf     = 1;
	vertex[vanz].Lorentz = new Lorentz_Function(lf::SSS);
	
	vertex[vanz].on      = 1;
	vanz++;
	//checked FK & RK
      }  
      
    }
  }
  
  // h0 -> h0/H0 + h0/H0 

    if(flh0.IsOn()) { 
      for(short int j=31;j<33;j++) {
	Flavour flav2 = Flavour(kf::code(j));
	if(flav2.IsOn()) {
	  for(short int k=31;k<33;k++) {
	    Flavour flav3 = Flavour(kf::code(k));
	    if(flav3.IsOn()) {
	      if(flav2 != flH0 || flav3 != flh0){

	      vertex[vanz].in[0] = flh0;
	      vertex[vanz].in[1] = flav2;      
	      vertex[vanz].in[2] = flav3;
	      
	      Kabbala C1 = K_Z_R(0,0)*K_Z_R(0,j-31)*K_Z_R(0,k-31);
	      Kabbala C2 = K_Z_R(1,0)*K_Z_R(1,j-31)*K_Z_R(1,k-31);
	      Kabbala C3 = K_Z_R(0,0)*K_Z_R(0,j-31)*K_Z_R(1,k-31);
	      Kabbala C4 = K_Z_R(0,0)*K_Z_R(1,j-31)*K_Z_R(0,k-31);
	      Kabbala C5 = K_Z_R(1,0)*K_Z_R(0,j-31)*K_Z_R(0,k-31);
	      Kabbala C6 = K_Z_R(1,0)*K_Z_R(1,j-31)*K_Z_R(0,k-31);
	      Kabbala C7 = K_Z_R(1,0)*K_Z_R(0,j-31)*K_Z_R(1,k-31);
	      Kabbala C8 = K_Z_R(0,0)*K_Z_R(1,j-31)*K_Z_R(1,k-31);
	      
	      kcpl0 = -M_I*g2*g2/(costW*costW*num_4)*
		(v1*(C1*num_3-C6-C7-C8)+v2*(C2*num_3-C3-C4-C5));
	      kcpl1 = kcpl0;

	      vertex[vanz].cpl[0]  = kcpl0.Value();
	      vertex[vanz].cpl[1]  = kcpl1.Value();
	      vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();	  	    		
	      vertex[vanz].cpl[2]  = 0.;vertex[vanz].cpl[3]  = 0.;
	
	      vertex[vanz].ncf   = 1;
	      vertex[vanz].Color = new Color_Function(cf::None); 
    
	      vertex[vanz].nlf     = 1;
	      vertex[vanz].Lorentz = new Lorentz_Function(lf::SSS);
	      
	      vertex[vanz].on      = 1;
	      vanz++;
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
      
      Kabbala C1 = K_Z_R(0,1)*K_Z_R(0,1)*K_Z_R(0,1);
      Kabbala C2 = K_Z_R(1,1)*K_Z_R(1,1)*K_Z_R(1,1);
      Kabbala C3 = K_Z_R(1,1)*K_Z_R(0,1)*K_Z_R(0,1);
      Kabbala C4 = K_Z_R(0,1)*K_Z_R(1,1)*K_Z_R(1,1);
      
      //noch einmal checken !!

      kcpl0 = -M_I*g2*g2/(costW*costW*num_4)*
	((C1-C4)*v1*num_3+(C2-C3)*v2*num_3);
      kcpl1 = kcpl0;

      vertex[vanz].cpl[0]  = kcpl0.Value();
      vertex[vanz].cpl[1]  = kcpl1.Value();
      vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();	  	    		
      vertex[vanz].cpl[2]  = 0.;vertex[vanz].cpl[3]  = 0.;
      
      vertex[vanz].ncf   = 1;
      vertex[vanz].Color = new Color_Function(cf::None); 
      
      vertex[vanz].nlf     = 1;
      vertex[vanz].Lorentz = new Lorentz_Function(lf::SSS);
      
      vertex[vanz].on      = 1;
      vanz++;
    }
      
  // A0 -> h0/H0 + A0
    
    if(flA0.IsOn()) {
      for(short int i=31;i<33;i++) {
	Flavour flav = Flavour(kf::code(i));
	if(flav.IsOn()) {
	  
	vertex[vanz].in[0] = flA0;
	vertex[vanz].in[1] = flav;
	vertex[vanz].in[2] = flA0;
	
	kcpl0 = -M_I/(costW*costW*num_4)*g2*g2*K_A_H(0,0)*K_B_R(i-31);
	kcpl1 = kcpl0;
	
	vertex[vanz].cpl[0]  = kcpl0.Value();
	vertex[vanz].cpl[1]  = kcpl1.Value();
	vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();	  	    		
	vertex[vanz].cpl[2]  = 0.;vertex[vanz].cpl[3]  = 0.;
	
	vertex[vanz].ncf   = 1;
	vertex[vanz].Color = new Color_Function(cf::None); 
	
	vertex[vanz].nlf     = 1;
	vertex[vanz].Lorentz = new Lorentz_Function(lf::SSS);
	
	vertex[vanz].on      = 1;
	vanz++;
	//checked FK & RK	
      }
    }
  }
 
}


void Interaction_Model_Higgs::c_SSV(Single_Vertex* vertex,int& vanz)
{   
  Flavour flW(kf::W);
  Flavour flZ(kf::Z);
  Flavour flHmin(kf::Hmin);
  Flavour flA0(kf::A0);
  Flavour flPhoton(kf::photon);
  Kabbala kcpl0,kcpl1;

  // A0 -> Z + h0/H0 
  
  if(flZ.IsOn()) {
    for(short int i=31; i<33;i++) {
      Flavour flav = Flavour(kf::code(i)); 
      if(flav.IsOn() && flA0.IsOn()) {
	vertex[vanz].in[0] = flA0;
	vertex[vanz].in[1] = flZ;
	vertex[vanz].in[2] = flav;

	kcpl0 = g2/(costW*num_2)*K_A_M(i-31,0);
	kcpl1 = kcpl0;

	vertex[vanz].cpl[0]  = kcpl0.Value();
	vertex[vanz].cpl[1]  = kcpl1.Value();
	vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();	  	    		
	vertex[vanz].cpl[2]  = 0.;vertex[vanz].cpl[3]  = 0.;
	
	vertex[vanz].ncf   = 1;
	vertex[vanz].Color = new Color_Function(cf::None); 
 
	vertex[vanz].nlf     = 1;
	vertex[vanz].Lorentz = new Lorentz_Function(lf::SSV);
	vertex[vanz].Lorentz->SetParticleArg(0,2,1);     	

	vertex[vanz ].on      = 1;
	vanz++;
	//checked FK & RK
      }  
    }
    
    // H- -> Z + H-
    
    if(flHmin.IsOn()) {
      vertex[vanz].in[0] = flHmin;
      vertex[vanz].in[1] = flZ;
      vertex[vanz].in[2] = flHmin;
     

      Complex cos2tw = (costW.Value()*costW.Value()-sintW.Value()*sintW.Value())/(2.*sintW.Value()*costW.Value());
      Kabbala cot2TW = Kabbala(string("cot2\\theta_W"),cos2tw);

      kcpl0 = M_I*g1*cot2TW;
      kcpl1 = kcpl0;

      vertex[vanz].cpl[0]  = kcpl0.Value();
      vertex[vanz].cpl[1]  = kcpl1.Value();
      vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();	  	    		
      vertex[vanz].cpl[2]  = 0.;vertex[vanz].cpl[3]  = 0.;
      
      vertex[vanz].ncf   = 1;
      vertex[vanz].Color = new Color_Function(cf::None); 
      
      vertex[vanz].nlf     = 1;
      vertex[vanz].Lorentz = new Lorentz_Function(lf::SSV);
      vertex[vanz].Lorentz->SetParticleArg(0,2,1);     	
      
      vertex[vanz].on      = 1;
      vanz++;
      //checked FK & RK
    }    
  }
  
  // H- -> Photon + H- 
  
  if(flPhoton.IsOn() && flHmin.IsOn()) {
    
    vertex[vanz].in[0] = flHmin;
    vertex[vanz].in[1] = flPhoton;
    vertex[vanz].in[2] = flHmin;
   
    kcpl0 = -M_I*g1;
    kcpl1 = kcpl0;

    vertex[vanz].cpl[0]  = kcpl0.Value();
    vertex[vanz].cpl[1]  = kcpl1.Value();
    vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();	  	    		
    vertex[vanz].cpl[2]  = 0.;vertex[vanz].cpl[3]  = 0.;
    
    vertex[vanz].ncf   = 1;
    vertex[vanz].Color = new Color_Function(cf::None); 
    
    vertex[vanz].nlf     = 1;
    vertex[vanz].Lorentz = new Lorentz_Function(lf::SSV);
    vertex[vanz].Lorentz->SetParticleArg(0,2,1);     	
          
    vertex[vanz].on      = 1;
    vanz++;
    //checked FK & RK
  }    
  
  //H- -> W- + h0/H0  
  
  if(flW.IsOn()) {
    for(short int i=31; i<33;i++) {
      Flavour flav = Flavour(kf::code(i)); 
      
      if(flav.IsOn() && flHmin.IsOn()) {
	
	vertex[vanz].in[0] = flHmin;
	vertex[vanz].in[1] = flW;  
	vertex[vanz].in[2] = flav;
	
	kcpl0 = -(M_I/num_2)*g2*K_A_M(i-31,0);
	kcpl1 = kcpl0;

	vertex[vanz].cpl[0]  = kcpl0.Value();
	vertex[vanz].cpl[1]  = kcpl1.Value();
	vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();	  	    		
	vertex[vanz].cpl[2]  = 0.;vertex[vanz].cpl[3]  = 0.;

	vertex[vanz].ncf   = 1;
	vertex[vanz].Color = new Color_Function(cf::None); 
 
	vertex[vanz].nlf     = 1;
	vertex[vanz].Lorentz = new Lorentz_Function(lf::SSV);
	vertex[vanz].Lorentz->SetParticleArg(0,2,1);     	
	
	vertex[vanz].on      = 1;
	vanz++;
	//checked FK & RK
      }
    }
    
    //H- -> W- +  A0
    
    if(flHmin.IsOn() && flA0.IsOn()) {
      vertex[vanz].in[0] = flHmin;
      vertex[vanz].in[1] = flW;    
      vertex[vanz].in[2] = flA0;

      //sollte so stimmen !

      kcpl0 = -g2/num_2;
      kcpl1 = kcpl0;

      vertex[vanz].cpl[0]  = kcpl0.Value();
      vertex[vanz].cpl[1]  = kcpl1.Value();
      vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();	  	    		
      vertex[vanz].cpl[2]  = 0.;vertex[vanz].cpl[3]  = 0.;
      
      vertex[vanz].ncf   = 1;
      vertex[vanz].Color = new Color_Function(cf::None); 
 
      vertex[vanz].nlf     = 1;
      vertex[vanz].Lorentz = new Lorentz_Function(lf::SSV);
      vertex[vanz].Lorentz->SetParticleArg(0,2,1);     	
      
      vertex[vanz].on      = 1;
      vanz++;
      //checked FK & RK
    }
  } 
}

void Interaction_Model_Higgs::c_SSVV(Single_Vertex* vertex,int& vanz) 
{
    
  Flavour flW(kf::W);
  Flavour flZ(kf::Z);
  Flavour flHmin(kf::Hmin);    
  Flavour flA0(kf::A0);
  Flavour flPhoton(kf::photon);
  
  Kabbala kcpl0,kcpl1;

  // Z - h0/H0/A0 - h0/H0/A0 - Z  
  if (flZ.IsOn()) {
    for(short int i=31; i<34;++i){
      Flavour flav = Flavour(kf::code(i));
      if(flav.IsOn()) {
	
 	vertex[vanz].in[0] = flZ;
	vertex[vanz].in[1] = flav;
	vertex[vanz].in[2] = flav;
	vertex[vanz].in[3] = flZ;
  
	vertex[vanz].nleg     = 4;

	kcpl0 = M_I*g2*g2/(num_2*costW*costW);
	kcpl1 = kcpl0;
    
	vertex[vanz].cpl[0]  = kcpl0.Value();
	vertex[vanz].cpl[1]  = kcpl1.Value();
	vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
	vertex[vanz].cpl[2]  = 0.;vertex[vanz].cpl[3]  = 0.;
  
	vertex[vanz].ncf   = 1;
	vertex[vanz].Color = new Color_Function(cf::None);     
	
	vertex[vanz].nlf     = 1;
	vertex[vanz].Lorentz = new Lorentz_Function(lf::VVSS);     
	vertex[vanz].Lorentz->SetParticleArg(0,3);     
	
	vertex[vanz].on      = 1;
	vanz++;
      }
    }
  }
 // W - h0/H0/A0 - h0/H0/A0 - W  
  if (flZ.IsOn()) {
    for(short int i=31; i<34;++i){
      Flavour flav = Flavour(kf::code(i));
      if(flav.IsOn()) {
	
 	vertex[vanz].in[0] = flW;
	vertex[vanz].in[1] = flav;
	vertex[vanz].in[2] = flav;
	vertex[vanz].in[3] = flW;
  
	vertex[vanz].nleg     = 4;

	kcpl0 = M_I*g2*g2/num_2;
	kcpl1 = kcpl0;
    
	vertex[vanz].cpl[0]  = kcpl0.Value();
	vertex[vanz].cpl[1]  = kcpl1.Value();
	vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
	vertex[vanz].cpl[2]  = 0.;vertex[vanz].cpl[3]  = 0.;
  
	vertex[vanz].ncf   = 1;
	vertex[vanz].Color = new Color_Function(cf::None);     
	
	vertex[vanz].nlf     = 1;
	vertex[vanz].Lorentz = new Lorentz_Function(lf::VVSS);     
	vertex[vanz].Lorentz->SetParticleArg(0,3);     
	
	vertex[vanz].on      = 1;
	vanz++;
      }
    }
  }

  if (flW.IsOn() && flHmin.IsOn()) {
    //W -> H- - H- - W
    vertex[vanz].in[0] = flW;
    vertex[vanz].in[1] = flHmin.Bar();
    vertex[vanz].in[2] = flHmin;
    vertex[vanz].in[3] = flW;
    
    vertex[vanz].nleg     = 4;

    kcpl0 = M_I*g2*g2/num_2;
    kcpl1 = kcpl0;
    
    vertex[vanz].cpl[0]  = kcpl0.Value();
    vertex[vanz].cpl[1]  = kcpl1.Value();
    vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
    vertex[vanz].cpl[2]  = 0.;vertex[vanz].cpl[3]  = 0.;
    
    vertex[vanz].ncf   = 1;
    vertex[vanz].Color = new Color_Function(cf::None);     
	
    vertex[vanz].nlf     = 1;
    vertex[vanz].Lorentz = new Lorentz_Function(lf::VVSS);
    vertex[vanz].Lorentz->SetParticleArg(0,3);     
    
    vertex[vanz].on      = 1;
    vanz++;
    
    // W+ -> H+ + h0/H0 + Z/P 
    for(short int i=31; i<33;++i){
      Flavour flav = Flavour(kf::code(i));
      if(flav.IsOn()) {
	if(flZ.IsOn()) {

 	vertex[vanz].in[0] = flW.Bar();
	vertex[vanz].in[1] = flHmin.Bar();
	vertex[vanz].in[2] = flav;
	vertex[vanz].in[3] = flZ;
  
	vertex[vanz].nleg     = 4;

	kcpl0 = M_I*g1*g1/(num_2*costW)*K_A_M(i-31,0);
	kcpl1 = kcpl0;
    
	vertex[vanz].cpl[0]  = kcpl0.Value();
	vertex[vanz].cpl[1]  = kcpl1.Value();
	vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
	vertex[vanz].cpl[2]  = 0.;vertex[vanz].cpl[3]  = 0.;
  
	vertex[vanz].ncf   = 1;
	vertex[vanz].Color = new Color_Function(cf::None);     
	
	vertex[vanz].nlf     = 1;
	vertex[vanz].Lorentz = new Lorentz_Function(lf::VVSS);
	vertex[vanz].Lorentz->SetParticleArg(0,3);     
	
	vertex[vanz].on      = 1;
	vanz++;
	}

	if(flPhoton.IsOn()) {

 	vertex[vanz].in[0] = flW.Bar();
	vertex[vanz].in[1] = flHmin.Bar();
	vertex[vanz].in[2] = flav;
	vertex[vanz].in[3] = flPhoton;
  
	vertex[vanz].nleg     = 4;

	kcpl0 = -M_I*g1*g1/(num_2*sintW)*K_A_M(i-31,0);
	kcpl1 = kcpl0;
    
	vertex[vanz].cpl[0]  = kcpl0.Value();
	vertex[vanz].cpl[1]  = kcpl1.Value();
	vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
	vertex[vanz].cpl[2]  = 0.;vertex[vanz].cpl[3]  = 0.;
  
	vertex[vanz].ncf   = 1;
	vertex[vanz].Color = new Color_Function(cf::None);     
	
	vertex[vanz].nlf     = 1;
	vertex[vanz].Lorentz = new Lorentz_Function(lf::VVSS);
	vertex[vanz].Lorentz->SetParticleArg(0,3);     
	
	vertex[vanz].on      = 1;
	vanz++;
	}
      }
    }
    // W+ -> H+ + A0 + Z/P   
    if (flA0.IsOn()) {
      if (flZ.IsOn()) {
	
 	vertex[vanz].in[0] = flW.Bar();
	vertex[vanz].in[1] = flHmin.Bar();
	vertex[vanz].in[2] = flA0;
	vertex[vanz].in[3] = flZ;
  
	vertex[vanz].nleg     = 4;
	
	kcpl0 = g1*g1/(num_2*costW);
	kcpl1 = kcpl0;
	
	vertex[vanz].cpl[0]  = kcpl0.Value();
	vertex[vanz].cpl[1]  = kcpl1.Value();
	vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
	vertex[vanz].cpl[2]  = 0.;vertex[vanz].cpl[3]  = 0.;
  
	vertex[vanz].ncf   = 1;
	vertex[vanz].Color = new Color_Function(cf::None);     
	
	vertex[vanz].nlf     = 1;
	vertex[vanz].Lorentz = new Lorentz_Function(lf::VVSS);
	vertex[vanz].Lorentz->SetParticleArg(0,3);     
	
	vertex[vanz].on      = 1;
	vanz++;
      }
      if (flPhoton.IsOn()) {
	
	vertex[vanz].in[0] = flW.Bar();
	vertex[vanz].in[1] = flHmin.Bar();
	vertex[vanz].in[2] = flA0;
	vertex[vanz].in[3] = flPhoton;
  
	vertex[vanz].nleg     = 4;
	
	kcpl0 = -g1*g1/(num_2*sintW);
	kcpl1 = kcpl0;
	
	vertex[vanz].cpl[0]  = kcpl0.Value();
	vertex[vanz].cpl[1]  = kcpl1.Value();
	vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
	vertex[vanz].cpl[2]  = 0.;vertex[vanz].cpl[3]  = 0.;
  
	vertex[vanz].ncf   = 1;
	vertex[vanz].Color = new Color_Function(cf::None);     
	
	vertex[vanz].nlf     = 1;
	vertex[vanz].Lorentz = new Lorentz_Function(lf::VVSS);
	vertex[vanz].Lorentz->SetParticleArg(0,3);     
	
	vertex[vanz].on      = 1;
	vanz++;
      }
    }
  }
  //P/Z -> H- - H- - Z/P 
  if (flHmin.IsOn()) {
    if (flZ.IsOn()) {
      
      Kabbala cot2TW = Kabbala(string("cot2\\theta_W"),
			       (costW.Value()*costW.Value()-sintW.Value()*sintW.Value())/(2.*sintW.Value()*sintW.Value()));
      
      vertex[vanz].in[0] = flZ;
      vertex[vanz].in[1] = flHmin;
      vertex[vanz].in[2] = flHmin.Bar();
      vertex[vanz].in[3] = flZ;
      
      vertex[vanz].nleg     = 4;
      
      kcpl0 = M_I*g1*g1*num_2*cot2TW*cot2TW;
      kcpl1 = kcpl0;
      
      vertex[vanz].cpl[0]  = kcpl0.Value();
      vertex[vanz].cpl[1]  = kcpl1.Value();
      vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
      vertex[vanz].cpl[2]  = 0.;vertex[vanz].cpl[3]  = 0.;
      
      vertex[vanz].ncf   = 1;
      vertex[vanz].Color = new Color_Function(cf::None);     
      
      vertex[vanz].nlf     = 1;
      vertex[vanz].Lorentz = new Lorentz_Function(lf::VVSS);
      vertex[vanz].Lorentz->SetParticleArg(0,3);     
      
      vertex[vanz].on      = 1;
      vanz++;
      
      if(flPhoton.IsOn()) {
	vertex[vanz].in[0] = flZ;
	vertex[vanz].in[1] = flHmin;
	vertex[vanz].in[2] = flHmin.Bar();
	vertex[vanz].in[3] = flPhoton;
	
	vertex[vanz].nleg     = 4;
	
	kcpl0 = M_I*g1*g1*num_2*cot2TW;
	kcpl1 = kcpl0;
	
	vertex[vanz].cpl[0]  = kcpl0.Value();
	vertex[vanz].cpl[1]  = kcpl1.Value();
	vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
	vertex[vanz].cpl[2]  = 0.;vertex[vanz].cpl[3]  = 0.;
  
	vertex[vanz].ncf   = 1;
	vertex[vanz].Color = new Color_Function(cf::None);     
	
	vertex[vanz].nlf     = 1;
	vertex[vanz].Lorentz = new Lorentz_Function(lf::VVSS);
	vertex[vanz].Lorentz->SetParticleArg(0,3);     
	
	vertex[vanz].on      = 1;
	vanz++;
	
      }
    }
    if(flPhoton.IsOn()) {
      vertex[vanz].in[0] = flPhoton;
      vertex[vanz].in[1] = flHmin;
      vertex[vanz].in[2] = flHmin.Bar();
      vertex[vanz].in[3] = flPhoton;
      
      vertex[vanz].nleg     = 4;
      
      kcpl0 = M_I*g1*g1*num_2;
      kcpl1 = kcpl0;
      
      vertex[vanz].cpl[0]  = kcpl0.Value();
      vertex[vanz].cpl[1]  = kcpl1.Value();
      vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
      vertex[vanz].cpl[2]  = 0.;vertex[vanz].cpl[3]  = 0.;
      
      vertex[vanz].ncf   = 1;
      vertex[vanz].Color = new Color_Function(cf::None);     
      
      vertex[vanz].nlf     = 1;
      vertex[vanz].Lorentz = new Lorentz_Function(lf::VVSS);
      vertex[vanz].Lorentz->SetParticleArg(0,3);     
      
      vertex[vanz].on      = 1;
      vanz++;
    }
  }    
}


Kabbala Interaction_Model_Higgs::K_CKM(short int i,short int j)       
{   
  char hi[2],hj[2];
  sprintf(hi,"%i",i);
  sprintf(hj,"%i",j);
  return Kabbala(string("V_{")+string(hi)+string(hj)+string("}"),
		 p_model->ComplexMatrixElement(std::string("CKM"),i,j));
} 
  
Kabbala Interaction_Model_Higgs::conj_K_CKM(short int i,short int j)       
{   
  char hi[2],hj[2];
  sprintf(hi,"%i",i);
  sprintf(hj,"%i",j);
  return Kabbala(string("V^\\ti_{")+string(hi)+string(hj)+string("}"),
		 conj(p_model->ComplexMatrixElement(std::string("CKM"),i,j)));
} 
 
inline Kabbala Interaction_Model_Higgs::K_yuk(Flavour fl) {
  return Kabbala(string("M_{"+fl.TexName()+"}"),fl.Yuk());
}

Kabbala Interaction_Model_Higgs::K_yuk_sign(Flavour fl) {
  char hi[3];
  sprintf(hi,"%i",fl.MassSign());
  return Kabbala(string(hi),fl.MassSign());
}

inline Kabbala Interaction_Model_Higgs::K_Z_H(short int i,short int j) {   
  char hi[2],hj[2];
  sprintf(hi,"%i",i);
  sprintf(hj,"%i",j);
  return Kabbala(string("Z^{")+string(hi)+string(hj)+string("}_H"),
		 p_model->ComplexMatrixElement(std::string("Z_H"),i,j));
}     

inline Kabbala Interaction_Model_Higgs::K_Z_R(short int i,short int j) {   
  char hi[2],hj[2];
  sprintf(hi,"%i",i);
  sprintf(hj,"%i",j);
  return Kabbala(string("Z^{")+string(hi)+string(hj)+string("}_R"),
		 p_model->ComplexMatrixElement(std::string("Z_R"),i,j));
}  

inline Kabbala Interaction_Model_Higgs::K_A_M(short int i,short int j) {
  char hi[2],hj[2];
  sprintf(hi,"%i",i);
  sprintf(hj,"%i",j);
  return Kabbala(string("A^{")+string(hi)+string(hj)+string("}_M"),
		 p_model->ComplexMatrixElement(std::string("Z_R"),0,i) *
		 p_model->ComplexMatrixElement(std::string("Z_H"),0,j) -
		 p_model->ComplexMatrixElement(std::string("Z_R"),1,i) *
		 p_model->ComplexMatrixElement(std::string("Z_H"),1,j));
}  

inline Kabbala Interaction_Model_Higgs::K_A_H(short int i,short int j) {
  char hi[2],hj[2];
  sprintf(hi,"%i",i);
  sprintf(hj,"%i",j);
  return Kabbala(string("A^{")+string(hi)+string(hj)+string("}_H"),
		 p_model->ComplexMatrixElement(std::string("Z_H"),0,i) *
		 p_model->ComplexMatrixElement(std::string("Z_H"),0,j) -
		 p_model->ComplexMatrixElement(std::string("Z_H"),1,i) *
		 p_model->ComplexMatrixElement(std::string("Z_H"),1,j));
}  

inline Kabbala Interaction_Model_Higgs::K_A_P(short int i,short int j) {
  char hi[2],hj[2];
  sprintf(hi,"%i",i);
  sprintf(hj,"%i",j);
  return Kabbala(string("A^{")+string(hi)+string(hj)+string("}_P"),
		 p_model->ComplexMatrixElement(std::string("Z_R"),0,i) *
		 p_model->ComplexMatrixElement(std::string("Z_H"),1,j) +
		 p_model->ComplexMatrixElement(std::string("Z_R"),1,i) *
		 p_model->ComplexMatrixElement(std::string("Z_H"),0,j));
}  

inline Kabbala Interaction_Model_Higgs::K_B_R(short int i) {   
  char hi[2];
  sprintf(hi,"%i",i);
  return Kabbala(string("B^{")+string(hi)+string("}_R"),
		 v1.Value() * p_model->ComplexMatrixElement(std::string("Z_R"),0,i) -
		 v2.Value() * p_model->ComplexMatrixElement(std::string("Z_R"),1,i) );
}  





















