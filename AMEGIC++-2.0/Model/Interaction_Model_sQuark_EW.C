#include "Interaction_Model_sQuark_EW.H"
#include "MathTools.H"
#include "Message.H"
#include <stdio.h>


using namespace AMEGIC;
using namespace AORGTOOLS;
using namespace APHYTOOLS;
using namespace AMATOOLS;
using namespace std;

Interaction_Model_sQuark_EW::Interaction_Model_sQuark_EW(MODEL::Model_Base * _model,
							 std::string _cplscheme,std::string _yukscheme) :
  Interaction_Model_Base(_model,_cplscheme,_yukscheme)
{ 
  double Ecms2 = sqr(rpa.gen.Ecms());

  g1       = Kabbala(string("g_1"),
		     sqrt(4.*M_PI*p_model->ScalarFunction(string("alpha_QED"),Ecms2)));
  g2       = Kabbala(string("g_1/\\sin\\theta_W"), 
		     g1.Value()/sqrt(p_model->ScalarConstant(string("sin2_thetaW"))));
  sintW    = Kabbala(string("\\sin\\theta_W"),
		     sqrt(p_model->ScalarConstant(string("sin2_thetaW"))));
  costW    = Kabbala(string("\\cos\\theta_W"),
		     sqrt(1.-p_model->ScalarConstant(string("sin2_thetaW"))));
  PL       = Kabbala(string("P_L"),1.);
  PR       = Kabbala(string("P_R"),1.);
  M_I      = Kabbala(string("i"),Complex(0.,1.));
  vev      = Kabbala(string("v_{EW}"),p_model->ScalarConstant(string("vev")));
  v1       = Kabbala(string("v_1"),
		     p_model->ScalarConstant(string("vev")) *
		     sqrt(1./(1.+sqr(p_model->ScalarConstant(string("tan(beta)"))))));
  v2       = Kabbala(string("v_2"),
		     p_model->ScalarConstant(string("vev")) *
		     p_model->ScalarConstant(string("tan(beta)")) *
		     sqrt(1./(1.+sqr(p_model->ScalarConstant(string("tan(beta)"))))));
  mu       = Kabbala(string("h"),p_model->ScalarConstant(string("mu")));
  conj_mu  = Kabbala(string("h"),p_model->ScalarConstant(string("mu")));
  K_zero   = Kabbala(string("zero"),0.);
  num_2    = Kabbala(string("2"),2.);    	
  num_3    = Kabbala(string("3"),3.);    	
  num_4    = Kabbala(string("4"),4.);    		
  root2    = Kabbala(string("\\sqrt{2}"),sqrt(2.));
  invroot2 = Kabbala(string("1/\\sqrt{2}"),sqrt(.5));
}


void Interaction_Model_sQuark_EW::c_FFS(Single_Vertex* vertex,int& vanz) {  
  Kabbala kcpl0,kcpl1; 
  
  //quark - squark - neutralino
  
  for (short int j=43;j<47;j++) {
    Flavour flneu = Flavour(kf::code(j));
    if (flneu.IsOn()) {
      //uptypes 
      for (short int k=2;k<7;k+=2) {
	Flavour flav1 = Flavour(kf::code(k));
	for (short int i=51;i<57;i++) {
	  Flavour flav2 = Flavour(kf::code(i));
	  if (flav1.IsOn() && flav2.IsOn()) {
	    vertex[vanz].in[0] = flav1;
	    vertex[vanz].in[1] = flav2;
	    vertex[vanz].in[2] = flneu;
	    
	    Kabbala K_u1 = Kabbala(string("u^1"),
				   (Flavour(kf::code(2)).Yuk())*sqrt(2.)/v2.Value());
	    Kabbala K_u2 = Kabbala(string("u^2"),
				   (Flavour(kf::code(4)).Yuk())*sqrt(2.)/v2.Value());
	    Kabbala K_u3 = Kabbala(string("u^3"),
				   (Flavour(kf::code(6)).Yuk())*sqrt(2.)/v2.Value());
	    
	    kcpl0 = M_I*(((g1*root2*num_2)/
			  (costW*num_3))*K_Z_U((k-2)/2+3,i-51)*K_Z_N(0,j-43)-
			 (K_u1*K_Z_U(0,i-51)+
			  K_u2*K_Z_U(1,i-51)+
			  K_u3*K_Z_U(2,i-51))*K_Z_N(3,j-43));
	    
	    kcpl1 = M_I*(-g2/(costW*root2)*K_Z_U((k-2)/2,i-51)*
			 (K_Z_N(0,j-43)*(sintW/num_3)+K_Z_N(1,j-43)*costW)-
		       (K_u1*K_Z_U(3,i-51)+
			K_u2*K_Z_U(4,i-51)+
			K_u3*K_Z_U(5,i-51))*K_Z_N(3,j-43));
	    
	    vertex[vanz].cpl[0] = kcpl0.Value();
	    vertex[vanz].cpl[1] = kcpl1.Value();
	    vertex[vanz].Str    = (kcpl0*PR+kcpl1*PL).String();
	    vertex[vanz].cpl[2]  = 0.;vertex[vanz].cpl[3]  = 0.;

	    vertex[vanz].ncf   = 1;
	    vertex[vanz].Color = new Color_Function; 
	    
	    vertex[vanz].Color->type       = cf::D;     
	    vertex[vanz].Color->SetParticleArg(0,1);     
	    vertex[vanz].Color->SetStringArg('0','1');     
	  
	    vertex[vanz].nlf     = 1;
	    vertex[vanz].Lorentz = new Lorentz_Function; 
	    
	    vertex[vanz].Lorentz->type       = lf::FFS;     
	    
	    vertex[vanz].on      = 1;
	    vanz++;
	  }
	}
      }
    }
  }

  for (short int j=43;j<47;j++) {
    Flavour flneu = Flavour(kf::code(j));
    if (flneu.IsOn()) {
      //downtypes 
      for (short int k=1;k<6;k+=2) {
	  Flavour flav1 = Flavour(kf::code(k));
	  for (short int i=61;i<67;i++) {
	    Flavour flav2 = Flavour(kf::code(i));
	    if (flav1.IsOn() && flav2.IsOn()) {
	      vertex[vanz].in[0] = flav1;
	      vertex[vanz].in[1] = flav2;
	      vertex[vanz].in[2] = flneu;

	      Kabbala K_dI = Kabbala(string("d^I"),
				     -flav1.Yuk()*sqrt(2.)/v1.Value());
	      Kabbala K_d1 = Kabbala(string("d^1"),
				     -Flavour(kf::code(1)).Yuk()*sqrt(2.)/v1.Value());
	      Kabbala K_d2 = Kabbala(string("d^2"),
				     -Flavour(kf::code(3)).Yuk()*sqrt(2.)/v1.Value());
	      Kabbala K_d3 = Kabbala(string("d^3"),
				     -Flavour(kf::code(5)).Yuk()*sqrt(2.)/v1.Value());

	      kcpl0 = M_I*((-(g1*root2)/
			    (costW*num_3))*K_Z_D((k-1)/2+3,i-61)*K_Z_N(0,j-43)+
			   (K_d1*K_Z_D(0,i-61)+
			    K_d2*K_Z_D(1,i-61)+
			    K_d3*K_Z_D(2,i-61))*K_Z_N(2,j-43));
	      
	      kcpl1 = M_I*(-g2/(costW*root2)*K_Z_D((k-1)/2,i-61)*
			   (K_Z_N(0,j-43)*(sintW/num_3)-K_Z_N(1,j-43)*costW)+
			   (K_d1*K_Z_D(3,i-61)+
			    K_d2*K_Z_D(4,i-61)+
			    K_d3*K_Z_D(5,i-61))*K_Z_N(2,j-43));
	      
	      vertex[vanz].cpl[0] = kcpl0.Value();
	      vertex[vanz].cpl[1] = kcpl1.Value();
	      vertex[vanz].Str    = (kcpl0*PR+kcpl1*PL).String();
	      vertex[vanz].cpl[2] = 0.;vertex[vanz].cpl[3]  = 0.;

	      vertex[vanz].ncf   = 1;
	      vertex[vanz].Color = new Color_Function; 
	      
	      vertex[vanz].Color->type       = cf::D;     
	      vertex[vanz].Color->SetParticleArg(0,1);     
	      vertex[vanz].Color->SetStringArg('0','1');     
	      
	      vertex[vanz].nlf     = 1;
	      vertex[vanz].Lorentz = new Lorentz_Function; 
	      
	      vertex[vanz].Lorentz->type       = lf::FFS;     
	    	
	      vertex[vanz].on     = 1;
	      vanz++;
	    }
	  }
      }
    }
  }

  //d-quark - Chargino - sup
  for (short int i=1;i<6;i+=2) {
    Flavour flav1 = Flavour(kf::code(i));
      
    Kabbala K_dI = Kabbala(string("d^I"),
			   -flav1.Yuk()/v1.Value()*sqrt(2.));
    
    for (short int j=41;j<43;j++) {
      Flavour flav2 = Flavour(kf::code(j));
      for (short int k=51;k<57;k++) {
	Flavour flav3 = Flavour(kf::code(k));
	if (flav1.IsOn() && flav2.IsOn() && flav3.IsOn()) {
	  vertex[vanz].in[0] = flav1;
	  vertex[vanz].in[1] = flav3;
	  vertex[vanz].in[2] = flav2;
	
	  Kabbala K_uI = Kabbala(string("u^I"),Flavour(kf::code(2*gen_sUp(flav3)+2)).Yuk()*
				 sqrt(2.)/v2.Value());

	  kcpl0 = -M_I*K_dI*K_Z_MI(1,j-41)*K_Z_U(gen_sUp(flav3),k-51)*
	    K_CKM((i-1)/2,gen_sUp(flav3));
	  
	  kcpl1 = M_I*(-g2*K_Z_PL(0,j-41)*K_Z_U(gen_sUp(flav3),k-51)+
		       K_uI*K_Z_PL(1,j-41)*K_Z_U(gen_sUp(flav3)+3,k-51))*
		       K_CKM((i-1)/2,gen_sUp(flav3));

	  vertex[vanz].cpl[0] = kcpl0.Value();
	  vertex[vanz].cpl[1] = kcpl1.Value();
	  vertex[vanz].Str    = (kcpl0*PR+kcpl1*PL).String();
	  vertex[vanz].cpl[2]  = 0.;vertex[vanz].cpl[3]  = 0.;
	  
	  vertex[vanz].ncf   = 1;
	  vertex[vanz].Color = new Color_Function; 
	
	  vertex[vanz].Color->type       = cf::D;     
	  vertex[vanz].Color->SetParticleArg(0,1);     
	  vertex[vanz].Color->SetStringArg('0','1');     
	  
	  vertex[vanz].nlf     = 1;
	  vertex[vanz].Lorentz = new Lorentz_Function; 
	      
	  vertex[vanz].Lorentz->type       = lf::FFS;     
	  
	  vertex[vanz].on      = 1;
	  vanz++;
	}
     }
    }
  }

  //u-quark - Chargino - sdown
  for (short int i=2;i<7;i+=2) {
    Flavour flav1 = Flavour(kf::code(i));
    
    Kabbala K_uJ = Kabbala(string("u^J"),flav1.Yuk()*sqrt(2.)/v2.Value());
    
    for (short int j=41;j<43;j++) {
      Flavour flav2 = Flavour(kf::code(j));
      for (short int k=61;k<67;k++) {
	Flavour flav3 = Flavour(kf::code(k));
	if (flav1.IsOn() && flav2.IsOn() && flav3.IsOn()) {
	  vertex[vanz].in[0] = flav1;
	  vertex[vanz].in[1] = flav3;
	  vertex[vanz].in[2] = flav2.Bar();

	  Kabbala K_dI = Kabbala(string("d^I"),
				 -Flavour(kf::code(2*gen_sDown(flav3)+1)).Yuk()*sqrt(2.)/
				 v1.Value());
	  
	  kcpl0 = M_I*K_uJ*K_Z_D(gen_sDown(flav3),k-61)*K_Z_PL(1,j-41)*
	    conj_K_CKM(gen_sDown(flav3),(i-2)/2);
	  
	  kcpl1 = -M_I*(g2*K_Z_D(gen_sDown(flav3),k-61)*K_Z_MI(0,j-41)+
			K_dI*K_Z_D(gen_sDown(flav3)+3,k-61)*K_Z_MI(1,j-41))*
			conj_K_CKM(gen_sDown(flav3),(i-2)/2);
	  
	  vertex[vanz].cpl[0] = kcpl0.Value();
	  vertex[vanz].cpl[1] = kcpl1.Value();
	  vertex[vanz].Str    = (kcpl0*PR+kcpl1*PL).String();
	  vertex[vanz].cpl[2] = 0.;vertex[vanz].cpl[3]  = 0.;
	  
	  vertex[vanz].ncf   = 1;
	  vertex[vanz].Color = new Color_Function; 
	  
	  vertex[vanz].Color->type       = cf::D;     
	  vertex[vanz].Color->SetParticleArg(0,1);     
	  vertex[vanz].Color->SetStringArg('0','1');     
	  
	  vertex[vanz].nlf     = 1;
	  vertex[vanz].Lorentz = new Lorentz_Function; 
	  
	  vertex[vanz].Lorentz->type       = lf::FFS;     
	  
	  vertex[vanz].on     = 1;
	  vanz++;
	}
      }
    } 
  }
}

void Interaction_Model_sQuark_EW::c_SSV(Single_Vertex* vertex,int& vanz)
{
  Kabbala kcpl0,kcpl1;
  
  //squark - Photon - squark
  
  Flavour flph = Flavour(kf::photon);
  if (flph.IsOn()) {
    //sUpypes
    for (short int i=51 ;i<57;i++) {
      Flavour flav = Flavour(kf::code(i));
      if (flav.IsOn()) {
	vertex[vanz].in[0] = flav;
	vertex[vanz].in[1] = flph;
	vertex[vanz].in[2] = flav;
	
	Kabbala charge = Kabbala(string("Q_{"+flav.TexName()+"}"),flav.Charge());
	
	//changed sign - -> +
	kcpl0 = M_I*charge*g1;
	kcpl1 = kcpl0;

	
	vertex[vanz].cpl[0]  = kcpl0.Value();
	vertex[vanz].cpl[1]  = kcpl1.Value();
	vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
	vertex[vanz].cpl[2]  = 0.;vertex[vanz].cpl[3]  = 0.;
	
	vertex[vanz].ncf   = 1;
	vertex[vanz].Color = new Color_Function; 
	
	vertex[vanz].Color->type       = cf::D;     
	vertex[vanz].Color->SetParticleArg(0,2);     
	vertex[vanz].Color->SetStringArg('0','2');     
	  
	vertex[vanz].nlf     = 1;
	vertex[vanz].Lorentz = new Lorentz_Function; 
	
	vertex[vanz].Lorentz->type       = lf::SSV;     
	vertex[vanz].Lorentz->SetParticleArg(0,2,1);     
	
	vertex[vanz].on      = 1;
	vanz++;
      }  
    }

    //sDowntypes
    for (short int i=61 ;i<67;i++) {
      Flavour flav = Flavour(kf::code(i));
      if (flav.IsOn()) {
	vertex[vanz].in[0] = flav;
	vertex[vanz].in[1] = flph;
	vertex[vanz].in[2] = flav;
	
	Kabbala charge = Kabbala(string("Q_{"+flav.TexName()+"}"),flav.Charge());

	kcpl0 = -M_I*charge*g1;
	kcpl1 = kcpl0;
	
	vertex[vanz].cpl[0]  = kcpl0.Value();
	vertex[vanz].cpl[1]  = kcpl1.Value();
	vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
	vertex[vanz].cpl[2]  = 0.;vertex[vanz].cpl[3]  = 0.;
	
	vertex[vanz].ncf   = 1;
	vertex[vanz].Color = new Color_Function; 
	
	vertex[vanz].Color->type       = cf::D;     
	vertex[vanz].Color->SetParticleArg(0,2);     
	vertex[vanz].Color->SetStringArg('0','2');     
	
	vertex[vanz].nlf     = 1;
	vertex[vanz].Lorentz = new Lorentz_Function; 
	
	vertex[vanz].Lorentz->type       = lf::SSV;     
	vertex[vanz].Lorentz->SetParticleArg(0,2,1);     
	
	vertex[vanz].on      = 1;
	vanz++;
      }
    }
  }   

  //squark - Z - squark

  Flavour flZ = Flavour(kf::Z);
  if (flZ.IsOn()) {
    //sUptypes
    for (short int i=51;i<57;i++) {
      Flavour flav1 = Flavour(kf::code(i));
      for (short int j=i;j<57;j++) {
	Flavour flav2 = Flavour(kf::code(j));
	if (flav1.IsOn() && flav2.IsOn()) {
	  
	  vertex[vanz].in[0] = flav1;
	  vertex[vanz].in[1] = flZ;
	  vertex[vanz].in[2] = flav2;
	  
	  Kabbala help = K_zero;
	  
	  if (i==j) {help = sintW*sintW*num_2/num_3;}  

	  kcpl0 = M_I*g2/costW*
	    (K_Z_U(gen_sUp(flav2),j-51)*K_Z_U(gen_sUp(flav2),i-51)/num_2-help);
	  kcpl1 = kcpl0;
	  
	  vertex[vanz].cpl[0]  = kcpl0.Value();
	  vertex[vanz].cpl[1]  = kcpl1.Value();
	  vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
	  vertex[vanz].cpl[2]  = 0.;vertex[vanz].cpl[3]  = 0.;
	
	  vertex[vanz].ncf   = 1;
	  vertex[vanz].Color = new Color_Function; 
	  
	  vertex[vanz].Color->type       = cf::D;     
	  vertex[vanz].Color->SetParticleArg(0,2);     
	  vertex[vanz].Color->SetStringArg('0','2');     
	  
	  vertex[vanz].nlf     = 1;
	  vertex[vanz].Lorentz = new Lorentz_Function; 
	  
	  vertex[vanz].Lorentz->type       = lf::SSV;     
	  vertex[vanz].Lorentz->SetParticleArg(0,2,1);     
	  
	  vertex[vanz].on      = 1;
	  vanz++;
	} 
      }
    }
    //sDowntypes
    for (short int i=61 ;i<67;i++) {
      Flavour flav1 = Flavour(kf::code(i));
      for (short int j=61 ;j<67;j++) {
	Flavour flav2 = Flavour(kf::code(j));
	if (flav1.IsOn() && flav2.IsOn()) {
	  
	  vertex[vanz].in[0] = flav1;
	  vertex[vanz].in[1] = flZ;
	  vertex[vanz].in[2] = flav2;
	  
	  Kabbala help = K_zero;
	  
	  if (i==j) { help = sintW*sintW/num_3; }  
	  
	  kcpl0 = M_I*g2/costW*
	    (K_Z_D(gen_sDown(flav2),j-61)*K_Z_D(gen_sDown(flav2),i-61)/num_2-help);
	  kcpl1 = kcpl0;

	  vertex[vanz].cpl[0]  = kcpl0.Value();
	  vertex[vanz].cpl[1]  = kcpl1.Value();
	  vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
	  vertex[vanz].cpl[2]  = 0.;vertex[vanz].cpl[3]  = 0.;
	  
	  vertex[vanz].ncf   = 1;
	  vertex[vanz].Color = new Color_Function; 
	
	  vertex[vanz].Color->type       = cf::D;     
	  vertex[vanz].Color->SetParticleArg(0,2);     
	  vertex[vanz].Color->SetStringArg('0','2');     
	  
	  vertex[vanz].nlf     = 1;
	  vertex[vanz].Lorentz = new Lorentz_Function; 
	
	  vertex[vanz].Lorentz->type       = lf::SSV;     
	  vertex[vanz].Lorentz->SetParticleArg(0,2,1);     
	  
	  vertex[vanz].on      = 1;
	  vanz++;
	}  
      }
    }
  }    
  
  //supquarks - W - sdownquarks
  Flavour flW = Flavour(kf::W);
  if (flW.IsOn()) {
    for (short int i=51;i<57;i++) {
      Flavour flav1 = Flavour(kf::code(i));
      for (short int j=61;j<67;j++) {
	Flavour flav2 = Flavour(kf::code(j));
	if (flav1.IsOn() && flav2.IsOn()) {
	  vertex[vanz].in[0] = flav2;
	  vertex[vanz].in[1] = flW;
	  vertex[vanz].in[2] = flav1;
	  
	  
	  kcpl0 = -M_I*g2*invroot2*(K_Z_D(gen_sDown(flav2),j-61)*K_Z_U(gen_sUp(flav1),i-51))*
	    K_CKM(gen_sDown(flav2),gen_sUp(flav1));
	  
	  kcpl1 = kcpl0;
	  
	  vertex[vanz].cpl[0]  = kcpl0.Value();
	  vertex[vanz].cpl[1]  = kcpl1.Value();
	  vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
	  vertex[vanz].cpl[2]  = 0.;vertex[vanz].cpl[3]  = 0.;
	  
	  vertex[vanz].ncf   = 1;
	  vertex[vanz].Color = new Color_Function; 
	  
	  vertex[vanz].Color->type       = cf::D;     
	  vertex[vanz].Color->SetParticleArg(0,2);     
	  vertex[vanz].Color->SetStringArg('0','2');     
	  
	  vertex[vanz].nlf     = 1;
	  vertex[vanz].Lorentz = new Lorentz_Function; 
	  
	  vertex[vanz].Lorentz->type       = lf::SSV;     
	  vertex[vanz].Lorentz->SetParticleArg(0,2,1);     
	  
	  vertex[vanz].on      = 1;
	  vanz++;
	} 
      }			
    }
  } 
}

void Interaction_Model_sQuark_EW::c_SSS(Single_Vertex* vertex,int& vanz)
{
  Kabbala kcpl0,kcpl1;
  Kabbala num_6 = Kabbala(string("6"),6.);
  
  //sQuarks - A0 - sQuarks
  
  Flavour flA0 = Flavour(kf::A0);
  if (flA0.IsOn()) {
    //uptypes
    for (short int i=51;i<57;i++) {
      Flavour flav1 = Flavour(kf::code(i));
      for (short int j=51;j<57;j++) {
	Flavour flav2 = Flavour(kf::code(j));
	if (flav1.IsOn() && flav2.IsOn() && i<=j) {
	  
	  Kabbala K_uI = Kabbala(string("u^I"),Flavour(kf::code(2*gen_sUp(flav1)+2)).Yuk()/
				 (v2).Value()*sqrt(2.));
	  
	  vertex[vanz].in[0] = flav1;
	  vertex[vanz].in[1] = flA0;
	  vertex[vanz].in[2] = flav2;
	  
	  kcpl0 = -(K_uI*K_Z_H(0,0)*(mu*K_Z_U(gen_sUp(flav1),i-51)*K_Z_U(gen_sUp(flav1)+3,j-51)-
				     conj_mu*K_Z_U(gen_sUp(flav1),j-51)*
				     K_Z_U(gen_sUp(flav1)+3,i-51))+
		    K_Z_H(1,0)*(K_u_S(gen_sUp(flav1),gen_sUp(flav2))*K_Z_U(gen_sUp(flav1),j-51)*
				K_Z_U(gen_sUp(flav2)+3,i-51)-K_u_S(gen_sUp(flav1),gen_sUp(flav2))*
				K_Z_U(gen_sUp(flav1),i-51)*K_Z_U(gen_sUp(flav2)+3,j-51))+
		    K_Z_H(0,0)*(K_w_S(gen_sUp(flav1),gen_sUp(flav2))*K_Z_U(gen_sUp(flav1),i-51)*
				K_Z_U(gen_sUp(flav2)+3,j-51)-K_w_S(gen_sUp(flav1),gen_sUp(flav2))*
				K_Z_U(gen_sUp(flav1),j-51)*K_Z_U(gen_sUp(flav2)+3,i-51)))
	    *invroot2;
	  
	  kcpl1 = kcpl0;
	  
	  vertex[vanz].cpl[0]  = kcpl0.Value();
	  vertex[vanz].cpl[1]  = kcpl1.Value();
	  vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
	  vertex[vanz].cpl[2]  = 0.;vertex[vanz].cpl[3]  = 0.;
	  
	  vertex[vanz].ncf   = 1;
	  vertex[vanz].Color = new Color_Function; 
	  
	  vertex[vanz].Color->type       = cf::D;     
	  vertex[vanz].Color->SetParticleArg(0,2);     
	  vertex[vanz].Color->SetStringArg('0','2');     
	  
	  vertex[vanz].nlf     = 1;
	  vertex[vanz].Lorentz = new Lorentz_Function; 
	  
	  vertex[vanz].Lorentz->type       = lf::SSS;     
	  
	  vertex[vanz].on      = 1;
	  vanz++;
	}
      }
    }
    //downtypes
    for (short int i=61;i<67;i++) {
      Flavour flav1 = Flavour(kf::code(i));
      for (short int j=61;j<67;j++) {
	Flavour flav2 = Flavour(kf::code(j));
	if (flav1.IsOn() && flav2.IsOn() && i<=j) {
	  
	  Kabbala K_dI = Kabbala(string("d^I"),
				 -Flavour(kf::code(2*gen_sDown(flav1)+1)).Yuk()/(v1).Value()*sqrt(2.));
	  
	  vertex[vanz].in[0] = flav1;
	  vertex[vanz].in[1] = flA0;
	  vertex[vanz].in[2] = flav2;
	  
	  kcpl0 = -(K_dI*K_Z_H(1,0)*(conj_mu*K_Z_D(gen_sDown(flav1),i-61)*
				     K_Z_D(gen_sDown(flav1)+3,j-61)-
				     mu*K_Z_D(gen_sDown(flav1),i-61)*
				     K_Z_D(gen_sDown(flav1)+3,j-61))+
		    K_Z_H(0,0)*(K_d_S(gen_sDown(flav1),gen_sDown(flav2))*
				K_Z_D(gen_sDown(flav1),j-61)*
				K_Z_D(gen_sDown(flav2)+3,i-61)-
				K_d_S(gen_sDown(flav1),gen_sDown(flav2))*
				K_Z_D(gen_sDown(flav1),i-61)*
				K_Z_D(gen_sDown(flav2)+3,j-61))+
		    K_Z_H(1,0)*(K_e_S(gen_sDown(flav1),gen_sDown(flav2))*
				K_Z_D(gen_sDown(flav1),j-61)*
				K_Z_D(gen_sDown(flav2)+3,i-61)-
				K_e_S(gen_sDown(flav1),gen_sDown(flav2))*
				K_Z_D(gen_sDown(flav1),i-61)*
				K_Z_D(gen_sDown(flav2)+3,j-61)))*invroot2;
	  
	  kcpl1 = kcpl0;
	  
	  vertex[vanz].cpl[0]  = kcpl0.Value();
	  vertex[vanz].cpl[1]  = kcpl1.Value();
	  vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
	  vertex[vanz].cpl[2]  = 0.;vertex[vanz].cpl[3]  = 0.;
	  
	  vertex[vanz].ncf   = 1;
	  vertex[vanz].Color = new Color_Function; 
      	  
	  vertex[vanz].Color->type       = cf::D;     
	  vertex[vanz].Color->SetParticleArg(0,2);     
	  vertex[vanz].Color->SetStringArg('0','2');     
	  
	  vertex[vanz].nlf     = 1;
	  vertex[vanz].Lorentz = new Lorentz_Function; 
	  
	  vertex[vanz].Lorentz->type       = lf::SSS;     
	  
	  vertex[vanz].on      = 1;
	  vanz++;
	}
      }
    }
  } 
  
  
  //sQuarks - h0/H0 - sQuarks
  
  for (short int k=31;k<33;k++) {
    Flavour flH = Flavour(kf::code(k)); 
    //uptypes  
    for (short int i=51;i<57;i++) {
      Flavour flav1 = Flavour(kf::code(i));
      for (short int j=51;j<57;j++) {
	Flavour flav2 =Flavour(kf::code(j));
	if(flav1.IsOn() && flav2.IsOn() && i<=j){
	  vertex[vanz].in[0] = flav1;
	  vertex[vanz].in[1] = flH;
	  vertex[vanz].in[2] = flav2;
	  
	  Kabbala help = K_zero;
	  
	  Kabbala K_uI = Kabbala(string("u^I"),Flavour(kf::code(2*gen_sUp(flav1)+2)).Yuk()/
				 (v2).Value()*sqrt(2.));
	  
	  Kabbala fac = Kabbala(string("\\frac{3-8sin^2\\theta_W}{4sin^2\\theta_W}"),
				(3.-8.*(sintW).Value()*(sintW).Value())/
				(4.*(sintW).Value()*(sintW).Value()));
	  
	  if (i==j) {help = Kabbala(string("1"),1.);}
	  
	  kcpl0 = M_I*(-g1*g1/(costW*costW*num_3)*K_B_R(k-31)*
		       (help+fac*K_Z_U(gen_sUp(flav1),i-51)*K_Z_U(gen_sUp(flav1),j-51))-
		       K_uI*K_uI*v2*K_Z_R(1,k-31)*
		       (K_Z_U(gen_sUp(flav1),i-51)*K_Z_U(gen_sUp(flav1),j-51)+
			K_Z_U(gen_sUp(flav1)+3,i-51)*K_Z_U(gen_sUp(flav1)+3,j-51))+
		       K_Z_R(1,k-31)*invroot2*K_u_S(gen_sUp(flav1),gen_sUp(flav2))*
		       (K_Z_U(gen_sUp(flav1),i-51)*K_Z_U(gen_sUp(flav2)+3,j-51)+
			K_Z_U(gen_sUp(flav1),j-51)*K_Z_U(gen_sUp(flav2)+3,i-51))+
		       K_Z_R(0,k-31)*invroot2*K_w_S(gen_sUp(flav1),gen_sUp(flav2))*
		       (K_Z_U(gen_sUp(flav1),i-51)*K_Z_U(gen_sUp(flav2)+3,j-51)+
			K_Z_U(gen_sUp(flav1),j-51)*K_Z_U(gen_sUp(flav2)+3,i-51))+
		       K_uI*K_Z_R(0,k-31)*invroot2*(conj_mu*K_Z_U(gen_sUp(flav1),j-51)*
						     K_Z_U(gen_sUp(flav1)+3,i-51)+
						     mu*K_Z_U(gen_sUp(flav1),i-51)*
						     K_Z_U(gen_sUp(flav2)+3,j-51)));
	  
	  kcpl1 = kcpl0;
	  
	  vertex[vanz].cpl[0]  = kcpl0.Value();
	  vertex[vanz].cpl[1]  = kcpl1.Value();
	  vertex[vanz].Str    = (kcpl0*PR+kcpl1*PL).String();
	  vertex[vanz].cpl[2]  = 0.;vertex[vanz].cpl[3]  = 0.;
	  
	  vertex[vanz].ncf   = 1;
	  vertex[vanz].Color = new Color_Function; 
      	  
	  vertex[vanz].Color->type       = cf::D;     
	  vertex[vanz].Color->SetParticleArg(0,2);     
	  vertex[vanz].Color->SetStringArg('0','2');     
	  
	  vertex[vanz].nlf     = 1;
	  vertex[vanz].Lorentz = new Lorentz_Function; 
	  
	  vertex[vanz].Lorentz->type       = lf::SSS;     
	  
	  vertex[vanz].on      = 1;
	  vanz++;
  	}
      }
    }
    
    //downtypes
    for (short int i=61;i<67;i++) {
      Flavour flav1 = Flavour(kf::code(i));
      for (short int j=61;j<67;j++) {
	Flavour flav2 =Flavour(kf::code(j));
	if(flav1.IsOn() && flav2.IsOn() && i<=j){
	  vertex[vanz].in[0] = flav1;
	  vertex[vanz].in[1] = flH;
	  vertex[vanz].in[2] = flav2;
	  
	  Kabbala help = K_zero;
	  
	  Kabbala K_dI = Kabbala(string("d^I"),
				 -Flavour(kf::code(2*gen_sDown(flav1)+1)).Yuk()/(v1).Value()*sqrt(2.));
	  
	  Kabbala fac = Kabbala(string("\\frac{3-4sin^2\\theta_W}{2sin^2\\theta_W}"),
				(3.-4.*(sintW).Value()*(sintW).Value())/
				(2.*(sintW).Value()*(sintW).Value()));
	  
	  if (i==j) {help = Kabbala(string("1"),1.);}
	  
	  kcpl0 = M_I*(-g1*g1/(costW*costW*num_6)*K_B_R(k-31)*
		       (help+fac*K_Z_D(gen_sDown(flav1),i-61)*K_Z_D(gen_sDown(flav1),j-61))-
		       K_dI*K_dI*v1*K_Z_R(0,k-31)*
		       (K_Z_D(gen_sDown(flav1),i-61)*K_Z_D(gen_sDown(flav1),j-61)+
			K_Z_D(gen_sDown(flav1)+3,i-61)*K_Z_D(gen_sDown(flav1)+3,j-61))+
		       K_Z_R(0,k-31)*invroot2*K_d_S(gen_sDown(flav1),gen_sDown(flav2))*
		       (K_Z_D(gen_sDown(flav1),j-61)*K_Z_D(gen_sDown(flav2)+3,i-61)+
			K_Z_D(gen_sDown(flav1),i-61)*K_Z_D(gen_sDown(flav2)+3,j-61))+
		       K_Z_R(1,k-31)*invroot2*K_e_S(gen_sDown(flav1),gen_sDown(flav2))*
		       (K_Z_D(gen_sDown(flav1),j-61)*K_Z_D(gen_sDown(flav2)+3,i-61)+
			K_Z_D(gen_sDown(flav1),i-61)*K_Z_D(gen_sDown(flav2)+3,j-61))-
		       K_dI*K_Z_R(1,k-31)*invroot2*(conj_mu*K_Z_D(gen_sDown(flav1),i-61)*
						     K_Z_D(gen_sDown(flav1)+3,j-61)+
						     mu*K_Z_D(gen_sDown(flav1),j-61)*
						     K_Z_D(gen_sDown(flav2)+3,i-61)));
	  
	  kcpl1 = kcpl0;
	  
	  vertex[vanz].cpl[0]  = kcpl0.Value();
	  vertex[vanz].cpl[1]  = kcpl1.Value();
	  vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
	  vertex[vanz].cpl[2]  = 0.;vertex[vanz].cpl[3]  = 0.;
	  
	  vertex[vanz].ncf   = 1;
	  vertex[vanz].Color = new Color_Function; 
      	  
	  vertex[vanz].Color->type       = cf::D;     
	  vertex[vanz].Color->SetParticleArg(0,2);     
	  vertex[vanz].Color->SetStringArg('0','2');     
	  
	  vertex[vanz].nlf     = 1;
	  vertex[vanz].Lorentz = new Lorentz_Function; 
	  
	  vertex[vanz].Lorentz->type       = lf::SSS;     
	  
	  vertex[vanz].on      = 1;
	  vanz++;
	  
	}
      }
    }
  }
  
  //sUp - Hmin - sDown
  
  Flavour flHmin = Flavour(kf::Hmin);
  if (flHmin.IsOn()) {
    for (short int i=51;i<57;i++) {
      Flavour flav1 = Flavour(kf::code(i));
      for (short int j=61;j<67;j++) {
	Flavour flav2 = Flavour(kf::code(j));
	vertex[vanz].in[0] = flav1.Bar();
	vertex[vanz].in[1] = flHmin.Bar();
	vertex[vanz].in[2] = flav2;
	
	Kabbala K_dI = Kabbala(string("d^I"),
			       -Flavour(kf::code(2*gen_sUp(flav1)+1)).Yuk()/(v1).Value()*sqrt(2.));
	
	Kabbala K_uJ = Kabbala(string("u^I"),Flavour(kf::code(2*gen_sDown(flav2)+2)).Yuk()/
			       (v2).Value()*sqrt(2.));
	
	Kabbala K_massW = Kabbala(string("M_W"),g1.Value()/2.*sqrt(v1.Value()*v1.Value()+
								   v2.Value()*v2.Value()));
	
	kcpl0 = M_I*((-(g2*g2)/num_2*(v1*K_Z_H(0,0)+v2*K_Z_H(1,0))+
		      v1*K_dI*K_dI*K_Z_H(0,0)+v2*K_uJ*K_uJ*K_Z_H(1,0))*invroot2*
		     conj_K_CKM(gen_sUp(flav1),gen_sDown(flav2))*
		     K_Z_D(gen_sUp(flav1),j-61)*K_Z_U(gen_sDown(flav2),i-51)-
		     sintW*K_massW*root2/g1*K_uJ*K_dI*
		     conj_K_CKM(gen_sUp(flav1),gen_sDown(flav2))*
		     K_Z_D(gen_sUp(flav1)+3,j-61)*K_Z_U(gen_sDown(flav2)+3,i-51)+
		     
		     (K_Z_H(0,0)*conj_mu*K_uJ*conj_K_CKM(gen_sUp(flav1),gen_sDown(flav2))+
		      (K_Z_H(0,0)*K_w_S(0,gen_sDown(flav2))-K_Z_H(1,0)*K_u_S(0,gen_sDown(flav2)))*
		      conj_K_CKM(gen_sUp(flav1),0)+
		      (K_Z_H(0,0)*K_w_S(1,gen_sDown(flav2))-K_Z_H(1,0)*K_u_S(1,gen_sDown(flav2)))*
		      conj_K_CKM(gen_sUp(flav1),1)+
		      (K_Z_H(0,0)*K_w_S(2,gen_sDown(flav2))-K_Z_H(1,0)*K_u_S(2,gen_sDown(flav2)))*
		      conj_K_CKM(gen_sUp(flav1),2))*
		     K_Z_U(gen_sDown(flav2)+3,i-51)*K_Z_D(gen_sUp(flav1),j-61)+
		     
		     ((K_Z_H(0,0)*K_d_S(0,gen_sUp(flav1))+K_Z_H(1,0)*K_e_S(0,gen_sUp(flav1)))*
		      conj_K_CKM(0,gen_sDown(flav2))+
		      (K_Z_H(0,0)*K_d_S(1,gen_sUp(flav1))+K_Z_H(1,0)*K_e_S(1,gen_sUp(flav1)))*
		      conj_K_CKM(1,gen_sDown(flav2))+
		      (K_Z_H(0,0)*K_d_S(2,gen_sUp(flav1))+K_Z_H(1,0)*K_e_S(2,gen_sUp(flav1)))*
		      conj_K_CKM(2,gen_sDown(flav2))*
		      K_Z_H(1,0)*mu*K_dI*conj_K_CKM(gen_sUp(flav1),gen_sDown(flav2)))*
		     K_Z_U(gen_sDown(flav2),i-51)*K_Z_D(gen_sUp(flav1)+3,j-61));
	
	kcpl1 = kcpl0;
	
	vertex[vanz].cpl[0]  = kcpl0.Value();
	vertex[vanz].cpl[1]  = kcpl1.Value();
	vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
	vertex[vanz].cpl[2]  = 0.;vertex[vanz].cpl[3]  = 0.;
	
	vertex[vanz].ncf   = 1;
	vertex[vanz].Color = new Color_Function; 
	
	vertex[vanz].Color->type       = cf::D;     
	vertex[vanz].Color->SetParticleArg(0,2);     
	vertex[vanz].Color->SetStringArg('0','2');     
	
	vertex[vanz].nlf     = 1;
	vertex[vanz].Lorentz = new Lorentz_Function; 
	
	vertex[vanz].Lorentz->type       = lf::SSS;     
	
	vertex[vanz].on      = 1;
	vanz++;
      }
    }
  }
}

void Interaction_Model_sQuark_EW::c_SSVV(Single_Vertex* vertex,int& vanz)
{
  Flavour flavW(kf::W);
  Flavour flavZ(kf::Z);
  Flavour flavPhoton(kf::photon);
  Kabbala kcpl0,kcpl1;
  
  // W - D - U - P/Z  
  for (short int i=51;i<57;i++) {
    Flavour flav1 = Flavour(kf::code(i));
    for (short int j=61;j<67;j++) {
      Flavour flav2 = Flavour(kf::code(j));
      if (flav1.IsOn() && flav2.IsOn()) {
	// W - D - U - P  
	if (flavW.IsOn()) {
	  if (flavPhoton.IsOn()) {
	    
	    vertex[vanz].in[0] = flavW;
	    vertex[vanz].in[1] = flav1.Bar();
	    vertex[vanz].in[2] = flav2;
	    vertex[vanz].in[3] = flavPhoton;
	    
	    vertex[vanz].nleg     = 4;
	    
	    //Attention, this is only a dummy coupling
	    kcpl0 = num_2;
	    kcpl1 = kcpl0;
	    
	    vertex[vanz].cpl[0]  = kcpl0.Value();
	    vertex[vanz].cpl[1]  = kcpl1.Value();
	    vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
	    vertex[vanz].cpl[2]  = 0.;vertex[vanz].cpl[3]  = 0.;
	    
	    vertex[vanz].ncf   = 1;
	    vertex[vanz].Color = new Color_Function; 
	    
	    vertex[vanz].Color->type       = cf::D;     
	    vertex[vanz].Color->SetParticleArg(1,2);     
	    vertex[vanz].Color->SetStringArg('1','2');     
	    
	    vertex[vanz].nlf     = 1;
	    vertex[vanz].Lorentz = new Lorentz_Function; 
	    
	    vertex[vanz].Lorentz->type = lf::VVSS;     
	    vertex[vanz].Lorentz->SetParticleArg(0,3);     
	    
	    //vertex[vanz].on      = 1;
	    vertex[vanz].on      = 0;
	    vanz++;
	  }
	  
	  if (flavZ.IsOn()) {
	    // W - D - U - Z  
	    
	    vertex[vanz].in[0] = flavW;
	    vertex[vanz].in[1] = flav1.Bar();
	    vertex[vanz].in[2] = flav2;
	    vertex[vanz].in[3] = flavZ;
	    
	    vertex[vanz].nleg     = 4;
	
	    //Attention, this is only a dummy coupling
	    kcpl0 = num_2;
	    kcpl1 = kcpl0;
	    
	    vertex[vanz].cpl[0]  = kcpl0.Value();
	    vertex[vanz].cpl[1]  = kcpl1.Value();
	    vertex[vanz].Str     = (kcpl0*PR+kcpl1*PL).String();
	    vertex[vanz].cpl[2]  = 0.;vertex[vanz].cpl[3]  = 0.;
	    
	    vertex[vanz].ncf   = 1;
	    vertex[vanz].Color = new Color_Function; 
	    
	    vertex[vanz].Color->type       = cf::D;     
	    vertex[vanz].Color->SetParticleArg(1,2);     
	    vertex[vanz].Color->SetStringArg('1','2');     
	    
	    vertex[vanz].nlf     = 1;
	    vertex[vanz].Lorentz = new Lorentz_Function; 
	    
	    vertex[vanz].Lorentz->type = lf::VVSS;     
	    vertex[vanz].Lorentz->SetParticleArg(0,3);     
	    
	    //vertex[vanz].on      = 1;
	    vertex[vanz].on      = 0;
	    vanz++;
	  }
	}
      }
    }
  }
}

Kabbala Interaction_Model_sQuark_EW::K_CKM(short int i,short int j)       
{   
  char hi[2],hj[2];
  sprintf(hi,"%i",i);
  sprintf(hj,"%i",j);
  return Kabbala(string("V_{")+string(hi)+string(hj)+string("}"),

		 p_model->ComplexMatrixElement(string("CKM"),i,j));
} 
  
Kabbala Interaction_Model_sQuark_EW::conj_K_CKM(short int i,short int j)       
{   
  char hi[2],hj[2];
  sprintf(hi,"%i",i);
  sprintf(hj,"%i",j);
  return Kabbala(string("V^\\ti_{")+string(hi)+string(hj)+string("}"),
		 conj(p_model->ComplexMatrixElement(string("CKM"),i,j)));
} 
 

Kabbala Interaction_Model_sQuark_EW::K_Z_D(short int i,short int j)       
{   
  char hi[2],hj[2];
  sprintf(hi,"%i",i);
  sprintf(hj,"%i",j);
  return Kabbala(string("Z^{")+string(hi)+string(hj)+string("}_D"),
		 p_model->ComplexMatrixElement(string("Z_d"),i,j));
}  

Kabbala Interaction_Model_sQuark_EW::K_Z_U(short int i,short int j)       
{   
  char hi[2],hj[2];
  sprintf(hi,"%i",i);
  sprintf(hj,"%i",j);
  return Kabbala(string("Z^{")+string(hi)+string(hj)+string("}_U"),
		 p_model->ComplexMatrixElement(string("Z_u"),i,j));
}  

Kabbala Interaction_Model_sQuark_EW::K_w_S(short int i,short int j)
{   
  char hi[2],hj[2];
  sprintf(hi,"%i",i);
  sprintf(hj,"%i",j);
  return Kabbala(string("w^{")+string(hi)+string(hj)+string("}_S"),
		 p_model->ComplexMatrixElement(string("ws"),i,j));

}  

Kabbala Interaction_Model_sQuark_EW::K_u_S(short int i,short int j)
{   
  char hi[2],hj[2];
  sprintf(hi,"%i",i);
  sprintf(hj,"%i",j);
  return Kabbala(string("u^{")+string(hi)+string(hj)+string("}_S"),
		 p_model->ComplexMatrixElement(string("us"),i,j));

}  

Kabbala Interaction_Model_sQuark_EW::K_e_S(short int i,short int j)
{   
  char hi[2],hj[2];
  sprintf(hi,"%i",i);
  sprintf(hj,"%i",j);
  return Kabbala(string("e^{")+string(hi)+string(hj)+string("}_S"),
		 p_model->ComplexMatrixElement(string("es"),i,j));
}  

Kabbala Interaction_Model_sQuark_EW::K_d_S(short int i,short int j)
{   
  char hi[2],hj[2];
  sprintf(hi,"%i",i);
  sprintf(hj,"%i",j);
  return Kabbala(string("d^{")+string(hi)+string(hj)+string("}_S"),
		 p_model->ComplexMatrixElement(string("ds"),i,j));
}

Kabbala Interaction_Model_sQuark_EW::K_yuk(Flavour fl) {
  return Kabbala(string("M_{"+fl.TexName()+"}"),fl.Yuk());
}

Kabbala Interaction_Model_sQuark_EW::K_yuk_sign(Flavour fl) {
  char hi[3];
  sprintf(hi,"%i",fl.MassSign());
  return Kabbala(string(hi),fl.MassSign());
}

Kabbala Interaction_Model_sQuark_EW::K_Z_H(short int i,short int j) {   
  char hi[2],hj[2];
  sprintf(hi,"%i",i);
  sprintf(hj,"%i",j);
  return Kabbala(string("Z^{")+string(hi)+string(hj)+string("}_H"),
		 p_model->ComplexMatrixElement(string("Z_H"),i,j));
}     

Kabbala Interaction_Model_sQuark_EW::K_Z_R(short int i,short int j) {   
  char hi[2],hj[2];
  sprintf(hi,"%i",i);
  sprintf(hj,"%i",j);
  return Kabbala(string("Z^{")+string(hi)+string(hj)+string("}_R"),
		 p_model->ComplexMatrixElement(string("Z_R"),i,j));
}  

Kabbala Interaction_Model_sQuark_EW::K_B_R(short int i) {   
  char hi[2];
  sprintf(hi,"%i",i);
  return Kabbala(string("B^{")+string(hi)+string("}_R"),
		 v1.Value() * p_model->ComplexMatrixElement(string("Z_R"),0,i) -
		 v2.Value() * p_model->ComplexMatrixElement(string("Z_R"),1,i) );
}  

Kabbala Interaction_Model_sQuark_EW::K_Z_PL(short int i,short int j)       
{   
  char hi[2];
  char hj[2];
  sprintf(hi,"%i",i);
  sprintf(hj,"%i",j);
  return Kabbala(string("Z^\\p_{")+string(hi)+string(hj)+string("}"),
		 p_model->ComplexMatrixElement(string("Z^+"),i,j));
}  

Kabbala Interaction_Model_sQuark_EW::K_Z_MI(short int i,short int j)       
{   
  char hi[2],hj[2];
  sprintf(hi,"%i",i);
  sprintf(hj,"%i",j);
  return Kabbala(string("Z^\\m_{")+string(hi)+string(hj)+string("}"),
		 p_model->ComplexMatrixElement(string("Z^-"),i,j));
}  

Kabbala Interaction_Model_sQuark_EW::K_Z_N(short int i,short int j)       
{   
  char hi[2],hj[2];
  sprintf(hi,"%i",i);
  sprintf(hj,"%i",j);
  return Kabbala(string("Z^{")+string(hi)+string(hj)+string("}_N"),
		 p_model->ComplexMatrixElement(string("Z^N"),i,j));
}  
//we use transposed convention !!! 

int Interaction_Model_sQuark_EW::gen_sUp(Flavour fl)
{
  int gen_sUp;

  if (fl.Kfcode() == 51 || fl.Kfcode() == 54)
    gen_sUp = 0;
  if (fl.Kfcode() == 52 || fl.Kfcode() == 55)
    gen_sUp = 1;
  if (fl.Kfcode() == 53 || fl.Kfcode() == 56)
    gen_sUp = 2;

  return gen_sUp;
}

int Interaction_Model_sQuark_EW::gen_sDown(Flavour fl)
{
  int gen_sDown;

  if (fl.Kfcode() == 61 || fl.Kfcode() == 64)
    gen_sDown = 0;
  if (fl.Kfcode() == 62 || fl.Kfcode() == 65)
    gen_sDown = 1;
  if (fl.Kfcode() == 63 || fl.Kfcode() == 66)
    gen_sDown = 2;

  return gen_sDown;
}
